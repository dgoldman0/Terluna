"""Independent continuous-density backward Monte Carlo spot checks.

Null-collision tracking; full scalar Rayleigh phase; Lambertian boundary.
A point Sun is used in checks safely away from disk contact. It has separate
spatial sampling from the deterministic solver and does not use its moment
or solar-transmittance lookup tables. Returns sampling standard errors.
"""
import math
import numpy as np
from numba import njit,prange
from solver import solar_columns,densities,EARTH,MOON,MOON_ZERO,SW,SF,O3

@njit(cache=True)
def rotate_direction(dx,dy,dz,c,phi):
    # Orthonormal tangent basis about direction d.
    if abs(dz)<.9:
        norm=math.sqrt(dx*dx+dy*dy);ux=-dy/norm;uy=dx/norm;uz=0.
    else:
        norm=math.sqrt(dy*dy+dz*dz);ux=0.;uy=-dz/norm;uz=dy/norm
    vx=dy*uz-dz*uy;vy=dz*ux-dx*uz;vz=dx*uy-dy*ux
    ss=math.sqrt(max(0.,1-c*c));cp=math.cos(phi);sp=math.sin(phi)
    return c*dx+ss*(cp*ux+sp*vx),c*dy+ss*(cp*uy+sp*vy),c*dz+ss*(cp*uz+sp*vz)

@njit(parallel=True,cache=True)
def monte_carlo(R,H,top,ozscale,beta,absorb,F,alb,sunangle,elevation,azimuth,ngroups,nper,seed):
    means=np.zeros(ngroups);squares=np.zeros(ngroups);trunc=np.zeros(ngroups)
    sx=math.cos(sunangle);sz=math.sin(sunangle)
    for group in prange(ngroups):
        np.random.seed(seed+group*103)
        for photon in range(nper):
            px=0.;py=0.;pz=R+1.7
            dx=math.cos(elevation)*math.cos(azimuth);dy=math.cos(elevation)*math.sin(azimuth);dz=math.sin(elevation)
            weight=1.;value=0.;events=0
            while weight>1e-12 and events<400:
                events+=1
                rr=px*px+py*py+pz*pz;r=math.sqrt(rr)
                rd=px*dx+py*dy+pz*dz;disc=rd*rd-rr+R*R
                ground=rd<0. and disc>0.
                if ground:end=max(0.,-rd-math.sqrt(disc))
                else:end=max(0.,-rd+math.sqrt(max(0.,rd*rd-rr+top*top)))
                rmin=R if ground else (math.sqrt(max(R*R,rr-rd*rd)) if rd<0 else r)
                major=beta*math.exp(-(rmin-R)/H)+absorb
                t=0.;hit=False
                while t<end:
                    t+=-math.log(max(1e-16,np.random.random()))/major
                    if t>=end:break
                    qx=px+t*dx;qy=py+t*dy;qz=pz+t*dz;qr=math.sqrt(qx*qx+qy*qy+qz*qz)
                    gas,oz=densities(qr-R,H,ozscale)
                    chi=beta*gas+absorb*oz
                    if np.random.random()<chi/major:
                        px=qx;py=qy;pz=qz;r=qr;hit=True
                        scatter=beta*gas/chi
                        ms=(px*sx+pz*sz)/r
                        cg,co=solar_columns(r,ms,R,top,H,ozscale,160)
                        if cg>=0:
                            nu=dx*sx+dz*sz
                            value+=weight*scatter*3./(16*math.pi)*(1+nu*nu)*F*math.exp(-beta*cg-absorb*co)
                        weight*=scatter
                        while True:
                            c=2*np.random.random()-1
                            if np.random.random()<.5*(1+c*c):break
                        dx,dy,dz=rotate_direction(dx,dy,dz,c,2*math.pi*np.random.random())
                        break
                if not hit:
                    if not ground:break
                    px+=dx*end;py+=dy*end;pz+=dz*end
                    r=math.sqrt(px*px+py*py+pz*pz)
                    nx=px/r;ny=py/r;nz=pz/r;ms=nx*sx+nz*sz
                    if ms>0:
                        cg,co=solar_columns(R+.01,ms,R,top,H,ozscale,160)
                        if cg>=0:value+=weight*alb/math.pi*ms*F*math.exp(-beta*cg-absorb*co)
                    weight*=alb
                    dx,dy,dz=rotate_direction(nx,ny,nz,math.sqrt(np.random.random()),2*math.pi*np.random.random())
                    px=nx*(R+.01);py=ny*(R+.01);pz=nz*(R+.01)
                if weight<.02:
                    survival=.2
                    if np.random.random()>survival:break
                    weight/=survival
            if events>=400:trunc[group]+=1
            means[group]+=value;squares[group]+=value*value
    return means.sum()/(ngroups*nper),math.sqrt(max(0.,squares.sum()/(ngroups*nper)-(means.sum()/(ngroups*nper))**2)/(ngroups*nper)),trunc.sum()

def check(atm,wavelength,sun_deg,elevation_deg,azimuth_deg,photons=60000,seed=82):
    beta=1.24062e-6*(wavelength/1000)**-4*atm.density_scale
    ab=float(np.interp(wavelength,SW,O3))*atm.ozone_du*2.687e20/(15000*atm.ozone_scale)
    F=float(np.interp(wavelength,SW,SF))
    mean,se,trunc=monte_carlo(atm.radius_m,atm.scale_height_m,atm.top,atm.ozone_scale,beta,ab,F,atm.ground_albedo,np.deg2rad(sun_deg),np.deg2rad(elevation_deg),np.deg2rad(azimuth_deg),40,photons//40,seed)
    return dict(world=atm.name,wavelength_nm=wavelength,sun_deg=sun_deg,elevation_deg=elevation_deg,azimuth_deg=azimuth_deg,photons=photons,mean=float(mean),standard_error=float(se),truncated_paths=int(trunc))
