"""Density-consistent spectral single-scattering reference and convergence records.

This first-order solver intentionally exposes scattering order and all spectral
inputs. It never substitutes its results silently into the shoreline renderer.
Run from any directory. Requires numpy, scipy and numba; no network access.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import time
import numpy as np
from numba import njit
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[2]
KB = 1.380649e-23
NREF = 101325 / (KB * 288.15)
DU = 2.687e20

@njit(cache=True)
def sample_air(z, rows):
    lo, hi = 0, rows.shape[0] - 1
    while hi - lo > 1:
        m = (lo + hi) // 2
        if rows[m, 0] <= z:
            lo = m
        else:
            hi = m
    t = min(1., max(0., (z - rows[lo, 0]) / (rows[hi, 0] - rows[lo, 0])))
    logp = rows[lo, 1] + t * (rows[hi, 1] - rows[lo, 1])
    temp = rows[lo, 2] + t * (rows[hi, 2] - rows[lo, 2])
    return math.exp(logp) / (KB * temp * NREF)

@njit(cache=True)
def ozone(z, settings):
    bottom, peak, top, amount = settings
    if z <= bottom or z >= top:
        return 0.
    frac = (z - bottom) / (peak - bottom) if z < peak else (top - z) / (top - peak)
    return 2 * frac / (top - bottom) * amount

@njit(cache=True)
def altitude(R, h, mu, s):
    r = R + h
    return ((2*R+h)*h + 2*r*mu*s + s*s) / (math.sqrt(r*r+2*r*mu*s+s*s)+R)

@njit(cache=True)
def ray_end(R, top, h, mu, length=1e99):
    r = R + h
    tangent2 = r*r*max(0., 1-mu*mu)
    end = -r*mu + math.sqrt(max(0., r*r*mu*mu+(top-h)*(2*R+top+h)))
    if mu < 0 and tangent2 < R*R:
        ground = max(0., -r*mu-math.sqrt(R*R-tangent2))
        if ground < min(end,length):
            return ground, True
    return min(end,length), False

@njit(cache=True)
def columns(R, top, h, mu, length, rows, oz, shells, nodes, weights):
    end, blocked = ray_end(R, top, h, mu, length)
    if blocked:
        return 0., 0., True
    points = np.empty(2*len(shells)+3)
    points[0], points[1], n = 0., end, 2
    r = R+h
    tangent = -r*mu
    if 0 < tangent < end:
        points[n] = tangent
        n += 1
    for z in shells:
        disc = r*r*mu*mu + (R+z-r)*(R+z+r)
        if disc >= 0:
            d = math.sqrt(disc)
            for s in (-r*mu-d, -r*mu+d):
                if 0 < s < end:
                    points[n] = s
                    n += 1
    pts = np.sort(points[:n])
    air, o3 = 0., 0.
    for i in range(1, n):
        mid, half = (pts[i]+pts[i-1])/2, (pts[i]-pts[i-1])/2
        if half < 1e-9:
            continue
        for k in range(len(nodes)):
            s = mid+half*nodes[k]
            z = max(0., min(top, altitude(R,h,mu,s)))
            air += half*weights[k]*sample_air(z,rows)
            o3 += half*weights[k]*ozone(z,oz)
    return air,o3,False

@njit(cache=True)
def sky_ray(R,top,h,elevation,azimuth,sun_elevation,rows,oz,shells,nodes,weights,beta,sigma,solar,steps):
    mu = math.sin(elevation)
    d = np.array([math.cos(elevation)*math.cos(azimuth),mu,math.cos(elevation)*math.sin(azimuth)])
    sun = np.array([math.cos(sun_elevation),math.sin(sun_elevation),0.])
    phase = 3*(1+np.dot(d,sun)**2)/(16*math.pi)
    end,blocked = ray_end(R,top,h,mu)
    out = np.zeros(len(beta))
    if blocked:
        return out
    tr = np.ones(len(beta))
    r = R+h
    for j in range(steps):
        a = end*(j/steps)**1.6
        b = end*((j+1)/steps)**1.6
        ds, s = b-a, (a+b)/2
        z = max(0.,min(top,altitude(R,h,mu,s)))
        air, o3 = sample_air(z,rows), ozone(z,oz)
        local_sun_mu = (r*sun[1]+s*np.dot(d,sun))/(R+z)
        sa,so,sb = columns(R,top,z,min(1.,max(-1.,local_sun_mu)),1e99,rows,oz,shells,nodes,weights)
        for k in range(len(beta)):
            scat = beta[k]*air
            ext = scat+sigma[k]*DU*o3
            opacity = -math.expm1(-ext*ds)
            if not sb and ext > 0:
                out[k] += tr[k]*opacity*(scat/ext)*solar[k]*phase*math.exp(-beta[k]*sa-sigma[k]*DU*so)
            tr[k] *= 1-opacity
    return out

@njit(cache=True)
def make_atlas(R,top,h,suns,el,az,rows,oz,shells,nodes,weights,beta,sigma,solar,steps):
    out = np.zeros((len(suns),len(el),len(az),len(beta)))
    for s in range(len(suns)):
        for i in range(len(el)):
            for j in range(len(az)):
                out[s,i,j] = sky_ray(R,top,h,el[i],az[j],suns[s],rows,oz,shells,nodes,weights,beta,sigma,solar,steps)
    return out

def unpack_profile(path: Path):
    raw=path.read_bytes()
    state=json.loads(raw)
    if state['schema']!='open-moon-atmospheric-profile/1':
        raise ValueError('Unsupported profile schema')
    rows=np.array([[v['z'],math.log(v['p']),v['T']] for v in state['rows']],dtype=np.float64)
    if not (np.isfinite(rows).all() and np.all(np.diff(rows[:,0])>0) and np.all(np.diff(rows[:,1])<0)):
        raise ValueError('Invalid shared profile grid')
    o=state['optics']['ozone']
    oz=np.array([o['bottom_m'],o['peak_m'],o['top_m'],o['column_DU']],dtype=np.float64)
    # Include every specified temperature/humidity break within the sounding.
    inputs=state['sounding']['inputs']
    xs=sorted(set(x[0] for k in ('temperature','humidity') for x in inputs[k] if x[0]<=2.2))
    heights=np.interp(xs,np.log(state['planet']['ps'])-rows[:,1],rows[:,0])
    shells=np.unique(np.concatenate((heights,[100,1000,5000,10000,25000,40000,80000,160000,320000,state['sounding']['end_m']],oz[:3])))
    shells=shells[(shells>0)&(shells<state['upper']['top_m'])]
    return state,rows,oz,shells,hashlib.sha256(raw).hexdigest()

def scipy_columns(state, rows, oz, shells, h, mu):
    """A separately expressed adaptive SciPy reference, avoiding production integrator."""
    R,top=state['planet']['radius'],state['upper']['top_m']
    r=R+h
    end,blocked=ray_end(R,top,h,mu)
    if blocked:
        return None
    breaks=[0.,float(end)]
    for z in np.unique(np.concatenate((shells, rows[:, 0]))):
        disc=r*r*mu*mu+(R+z-r)*(R+z+r)
        if disc>=0:
            for s in [-r*mu-math.sqrt(disc),-r*mu+math.sqrt(disc)]:
                if 0<s<end: breaks.append(s)
    if 0<-r*mu<end: breaks.append(-r*mu)
    ordered=sorted(set(breaks))
    breaks=[ordered[0]]
    for s in ordered[1:]:
        if s-breaks[-1]>1e-5: breaks.append(s)
    def f(s,which):
        z=((2*R+h)*h+2*r*mu*s+s*s)/(math.sqrt(r*r+s*s+2*r*mu*s)+R)
        if which==0:
            return math.exp(np.interp(z,rows[:,0],rows[:,1]))/(KB*np.interp(z,rows[:,0],rows[:,2])*NREF)
        b,p,t,amount=oz
        return amount*2/(t-b)*max(0.,min((z-b)/(p-b),(t-z)/(t-p)))
    return [sum(quad(lambda s:f(s,k),a,b,epsabs=1e-7 if k==0 else 1e-9,epsrel=2e-8,limit=200)[0] for a,b in zip(breaks[:-1],breaks[1:])) for k in (0,1)]

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--profiles',type=Path,default=ROOT/'data/a1/generated')
    ap.add_argument('--ids',nargs='*')
    ap.add_argument('--width',type=int,default=33)
    ap.add_argument('--height',type=int,default=25)
    ap.add_argument('--steps',type=int,default=96)
    ap.add_argument('--gauss',type=int,default=12)
    args=ap.parse_args()
    if min(args.width,args.height)<3 or args.steps<8 or args.gauss<4: ap.error('Resolution below supported minimum')
    manifest=json.loads((args.profiles/'profiles.json').read_text())
    spectral_path=ROOT/'data/a1/spectral-inputs.json'
    spec=json.loads(spectral_path.read_text())
    wavelengths=np.array(spec['wavelength_nm'],dtype=float)
    beta=1.24062e-6*(wavelengths/1000)**-4
    sigma,solar=np.array(spec['ozone_m2']),np.array(spec['solar_W_m2_nm'])
    cie=np.array(spec['cie1931_xyz']); rgb_matrix=np.array(spec['xyz_to_linear_srgb'])
    spectral_weights=np.full(len(wavelengths),10.);spectral_weights[[0,-1]]=5.
    rgb_weights=(cie.T*spectral_weights*683).T@rgb_matrix.T
    el=np.deg2rad(90*(np.arange(args.height)/(args.height-1))**2)
    az=np.linspace(0,math.pi,args.width)
    suns_deg=np.array([-6.,-2.,0.,2.,6.,15.,45.,90.]);suns=np.deg2rad(suns_deg)
    nodes,weights=np.polynomial.legendre.leggauss(args.gauss)
    high_n,high_w=np.polynomial.legendre.leggauss(24)
    start=time.monotonic()
    for record in manifest['profiles']:
        if args.ids and record['id'] not in args.ids: continue
        state,rows,oz,shells,sha=unpack_profile(args.profiles/record['file'])
        if sha!=record['profile_sha256']: raise ValueError('Profile hash mismatch')
        R,top=state['planet']['radius'],state['upper']['top_m']
        validation=[]
        for ref in record['independent_js_paths']:
            h,mu=ref['h_m'],ref['mu'];a,o,b=columns(R,top,h,mu,1e99,rows,oz,shells,high_n,high_w)
            independent=scipy_columns(state,rows,oz,shells,h,mu)
            error=None if b else [abs(a-independent[0])/max(1,independent[0]),abs(o-independent[1])/max(1,independent[1])]
            validation.append({'h_m':h,'mu':mu,'blocked':bool(b),'gauss24':[a,o] if not b else None,'scipy':independent,'js':[ref['air_m'],ref['ozone_DU']],'relative_to_scipy':error})
        spectra=make_atlas(R,top,2.,suns,el,az,rows,oz,shells,nodes,weights,beta,sigma,solar,args.steps)
        rgb=spectra@rgb_weights
        # Numerical radiance comparison keeps profile, spectrum and scattering order fixed.
        rays=[]
        for sun,e,a in [(45,90,0),(45,15,90),(6,1,0),(-2,5,0),(0,0,90),(90,30,180)]:
            x=sky_ray(R,top,2.,math.radians(e),math.radians(a),math.radians(sun),rows,oz,shells,nodes,weights,beta,sigma,solar,args.steps)
            y=sky_ray(R,top,2.,math.radians(e),math.radians(a),math.radians(sun),rows,oz,shells,high_n,high_w,beta,sigma,solar,args.steps*4)
            rays.append({'sun_deg':sun,'viewElevation_deg':e,'viewAzimuth_deg':a,'Y':float((x@rgb_weights)@[.2126,.7152,.0722]),'referenceY':float((y@rgb_weights)@[.2126,.7152,.0722]),'relativeSpectralL1':float(np.sum(np.abs(x-y))/max(1e-12,np.sum(y)))})
        direct=[]
        for s in suns:
            air,o3,b=columns(R,top,2.,math.sin(s),1e99,rows,oz,shells,high_n,high_w)
            v=np.zeros(len(beta)) if b else solar*np.exp(-beta*air-sigma*DU*o3)
            direct.append((v@rgb_weights).tolist())
        result={'schema':'open-moon-a1-clear-reference/1','id':record['id'],'profile_sha256':sha,
            'spectral_sha256':hashlib.sha256(spectral_path.read_bytes()).hexdigest(),
            'solver_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'model':{'scatteringOrders':1,'molecules':'dry-air-equivalent Rayleigh from total molecular number density','ozone':'selected geometric triangle, 233 K spectrum','surfaceBoundary':'black','aerosols':0,'refraction':False,'solarAngularExtent':'collimated beam for transport; display omits solar disk','observerHeight_m':2},
            'budgets':{'viewSteps':args.steps,'sunGaussOrder':args.gauss,'comparisonViewSteps':args.steps*4,'comparisonSunGaussOrder':24},
            'suns_deg':suns_deg.tolist(),'elevations_deg':np.rad2deg(el).tolist(),'azimuths_deg':np.rad2deg(az).tolist(),'radiance_rgb_units':'photopic-weighted linear sRGB, cd/m2-equivalent; negative gamut channels retained','direct_rgb_units':'photopic-weighted normal illuminance, linear sRGB lux-equivalent','radiance_rgb':np.round(rgb,8).tolist(),'direct_rgb':direct,
            'column_validation':validation,'ray_convergence':rays,'scope':'First-order optical reference. Higher scattering orders, spectral-bin convergence, climate validity and production-image parity remain separate gates.'}
        output=args.profiles/(record['id']+'.sky.json');output.write_text(json.dumps(result,separators=(',',':'))+'\n')
        np.savez_compressed(args.profiles/(record['id']+'.spectra.npz'),wavelength_nm=wavelengths,radiance=spectra,profile_sha256=sha)
        print(json.dumps({'id':record['id'],'elapsed_s':time.monotonic()-start,'maximum_column_relative_error':max((max(v['relative_to_scipy']) for v in validation if v['relative_to_scipy']),default=0),'maximum_ray_spectral_L1':max(v['relativeSpectralL1'] for v in rays)}),flush=True)

if __name__=='__main__': main()
