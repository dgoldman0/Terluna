"""Explicit, repeating finite-tile traffic in distinct meridional planes.

This is an inverse-dynamics construction, not a free-orbit propagation. Every
integer (plane, slot) identifies one tile at a common absolute ephemeris time.
The two alternating radii in a plane share a clock; a small smooth phase warp
registers their sunward shadows. Return arcs are part of the same closed path.
No ephemeris trajectory is shifted in time to manufacture another member.
"""
from dataclasses import dataclass

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.spatial import cKDTree

from shared import constants as K
from .active_formation import solar_frame
from .cycling import ray_geometry, visible_sun
from .fleet import beam_margins
from .optical import length, unit, arriving_ray_vectors


@dataclass(frozen=True)
class Layout:
    planes: int = 2304
    first_radius: float = 15e6
    radial_step: float = 4500.
    pitch: float = 9500.
    side: float = 10000.
    clear_side: float = 9890.
    protected_radius: float = 4*K.MOON_RADIUS
    normal_limit_deg: float = 17.


def phase_warp(theta, ratio):
    """Angle and first two derivatives; exact projected registration in service."""
    x = np.arctan2(np.sin(theta), np.cos(theta))
    z = np.clip((np.abs(x)-.65)/.30, 0, 1)
    w = 1-10*z**3+15*z**4-6*z**5
    wp = (-30*z*z+60*z**3-30*z**4)*np.sign(x)/.30
    wpp = (-60*z+180*z*z-120*z**3)/.30**2
    # The unwarped inner parity needs no inverse-sine branch at quadrature.
    changed = ratio < 1-1e-14
    rr = np.where(changed, ratio, .999)
    d = np.maximum(1-(rr*np.sin(x))**2, 1e-30)
    g = np.arcsin(rr*np.sin(x))-x
    gp = rr*np.cos(x)/np.sqrt(d)-1
    gpp = -rr*np.sin(x)/np.sqrt(d)+rr**3*np.sin(x)*np.cos(x)**2/d**1.5
    g, gp, gpp = [np.where(changed, v, 0.) for v in (g,gp,gpp)]
    return theta+w*g, 1+wp*g+w*gp, wpp*g+2*wp*gp+w*gpp


class Traffic:
    def __init__(self, env, days=30., layout=Layout()):
        self.env, self.layout = env, layout
        n = layout.planes
        if n < 4 or n % 2:
            raise ValueError('Use an even number of planes')
        self.alpha = np.arange(n)*np.pi/n
        # Cyclic neighbouring planes differ by at most two radius ranks.
        # A monotone assignment would leave a large parallax jump at its wrap.
        rank = np.r_[np.arange(0,n,2), np.arange(n-1,0,-2)]
        self.r0 = layout.first_radius+2*layout.radial_step*rank
        self.count = (2*np.ceil(2*np.pi*self.r0/layout.pitch/2)).astype(int)
        self.omega = np.sqrt(K.MOON_GM/(self.r0+layout.radial_step/2)**3)
        self.omega *= np.where(np.cos(self.alpha)>=0, 1., -1.)
        self.phase = 2*np.pi*np.mod(np.arange(n)*.6180339887498949,1)/self.count
        ts=np.arange(-1200,days*K.JULIAN_DAY+1801,600.)
        self.frames=CubicSpline(ts,np.array([solar_frame(env.at(t)) for t in ts]))
        # Remove the in-plane rotation already carried by the solar frame.
        # Otherwise the nominal lunar circular speed is counted twice in that
        # component, needlessly charging the controller for a radial mismatch.
        spins=[]
        for t in np.linspace(0,days*K.JULIAN_DAY,61):
            w=self.frames(t).T@self.frames(t,1)
            spins.append([w[2,1],w[0,2],w[1,0]])
        spin=np.mean(spins,axis=0)
        self.omega-= -np.sin(self.alpha)*spin[1]+np.cos(self.alpha)*spin[2]

    @property
    def inventory(self):
        return int(self.count.sum())

    def path(self, t, plane, slot=None, phase=None, parity=None):
        """Position, velocity and acceleration of each identified member.

        Explicit phase/parity is used only for spatial quadrature and the broad
        ray-search envelope; finite coverage always resolves an integer slot.
        """
        j=np.asarray(plane,dtype=int)
        if slot is not None:
            k=np.asarray(slot,dtype=int)%self.count[j]
            parity=k%2
            phase=2*np.pi*k/self.count[j]+self.phase[j]+self.omega[j]*t
        r=self.r0[j]+np.asarray(parity)*self.layout.radial_step
        g,gp,gpp=phase_warp(np.asarray(phase),self.r0[j]/r)
        ca,sa=np.cos(self.alpha[j]),np.sin(self.alpha[j])
        cg,sg=np.cos(g),np.sin(g)
        qr=r[...,None]*np.stack([cg,sg*ca,sg*sa],axis=-1)
        tangent=np.stack([-sg,cg*ca,cg*sa],axis=-1)
        vr=(r*gp*self.omega[j])[...,None]*tangent
        ar=-(gp*self.omega[j])[...,None]**2*qr+(r*gpp*self.omega[j]**2)[...,None]*tangent
        f,fd,fdd=(self.frames(t,o) for o in range(3))
        q=qr@f.T; v=vr@f.T+qr@fd.T
        a=ar@f.T+2*vr@fd.T+qr@fdd.T
        return np.concatenate([q,v],axis=-1),a,tangent@f.T

    def facets(self,t,plane,slot=None,phase=None,parity=None):
        state,acc,tangent=self.path(t,plane,slot,phase,parity)
        sample=self.env.at(t); q=state[...,:3]
        rays=ray_geometry(q,sample); e=-unit(rays['sun'])
        n=unit(q); n*=np.where(np.sum(n*e,axis=-1)>=0,1.,-1.)[...,None]
        original=n.copy(); changed=np.zeros(q.shape[:-1],bool)
        fallback=np.broadcast_to(self.frames(t)[:,1],q.shape)
        # Smooth out-of-plane detour around the Moon on the return hemisphere.
        # Without it, a nearest-cone projection flips the normal across the
        # lunar cone at the return pole. Earth handling below remains an ideal
        # pointing law, whose slew feasibility must be evaluated separately.
        out=e-2*np.sum(e*n,axis=-1)[...,None]*n
        moon=-unit(q)
        sep=np.arctan2(length(np.cross(out,moon)),np.sum(out*moon,axis=-1))
        old_margins=beam_margins(q,rays['sun'],n,sample,self.layout.protected_radius,self.layout.side)
        cone=sep-old_margins[...,0]+.003
        h=unit(np.cross(q,tangent))
        lift=1.1*np.tan(cone)*np.exp(-.5*(sep/cone)**2)
        z=np.clip((-np.sum(unit(q)*self.frames(t)[:,0],axis=-1)-.3)/.35,0,1)
        lift*=10*z**3-15*z**4+6*z**5
        candidate=unit(e-unit(out+lift[...,None]*h))
        n=np.where((lift>1e-10)[...,None],candidate,n)
        for _ in range(8):
            margins=beam_margins(q,rays['sun'],n,sample,self.layout.protected_radius,self.layout.side)
            if np.min(margins)>=.002-1e-12:break
            for h,target in enumerate((-q,sample['positions']['earth']-q)):
                bad=margins[...,h]<.002
                if not np.any(bad):continue
                changed|=bad
                outgoing=e-2*np.sum(e*n,axis=-1)[...,None]*n
                target=unit(target)
                away=outgoing*np.sum(outgoing*target,axis=-1)[...,None]-target
                # This ideal instantaneous pointing law has no slew certificate.
                away=np.where((length(away)<1e-9)[...,None],fallback,away)
                away=unit(away)
                delta=np.maximum(.002-margins[...,h],0)
                out=outgoing*np.cos(delta)[...,None]+away*np.sin(delta)[...,None]
                new=unit(e-out)
                n=np.where(bad[...,None],new,n)
        a=unit(tangent-n*np.sum(tangent*n,axis=-1)[...,None]); b=np.cross(n,a)
        margins=beam_margins(q,rays['sun'],n,sample,self.layout.protected_radius,self.layout.side)
        angle=np.arccos(np.clip(np.abs(np.sum(n*original,axis=-1)),0,1))
        bad=(angle>np.deg2rad(self.layout.normal_limit_deg))|(np.min(margins,axis=-1)<.0019)
        if np.any(bad):
            # In-plane alternatives bound the physical radial envelope even
            # when a grazing-incidence nearest-cone correction is ill-defined.
            flat_q=q.reshape(-1,3); flat_e=e.reshape(-1,3)
            ids=np.flatnonzero(bad); nn=original.reshape(-1,3)[ids]
            tt=unit(tangent.reshape(-1,3)[ids]-nn*np.sum(tangent.reshape(-1,3)[ids]*nn,axis=-1)[:,None])
            delta=np.deg2rad(np.arange(0,self.layout.normal_limit_deg+.0001,.1))
            delta=np.ravel(np.column_stack([delta,-delta]))
            opts=nn[:,None,:]*np.cos(delta)[None,:,None]+tt[:,None,:]*np.sin(delta)[None,:,None]
            margin=beam_margins(flat_q[ids,None,:],rays['sun'].reshape(-1,3)[ids,None,:],opts,
                sample,self.layout.protected_radius,self.layout.side)
            valid=np.min(margin,axis=-1)>=.002
            if not np.all(np.any(valid,axis=-1)):
                raise ValueError('No beam-safe attitude within the normal envelope')
            selection=valid.argmax(axis=-1)
            n.reshape(-1,3)[ids]=opts[np.arange(len(ids)),selection]
            a=unit(tangent-n*np.sum(tangent*n,axis=-1)[...,None]);b=np.cross(n,a)
            margins=beam_margins(q,rays['sun'],n,sample,self.layout.protected_radius,self.layout.side)
            angle=np.arccos(np.clip(np.abs(np.sum(n*original,axis=-1)),0,1))
        # Reflection geometry is invariant to n -> -n, but the momentum
        # vector must face the incident radiation, including fallback choices.
        sign=np.where(np.sum(n*e,axis=-1)>=0,1.,-1.)
        n*=sign[...,None];b*=sign[...,None]
        return dict(state=state,acc=acc,normal=n,a=a,b=b,rays=rays,
                    margins=margins,tilt=angle,changed=changed)

    def power_bounds(self,t,plane,phase,parity):
        """Force bracket valid for ANY mutual shadow fraction in [0,1].

        We do not assign uncomputed sunlight to mutually shadowed tiles.
        Gravity is centre point-mass; a finite-extent bound is reported apart.
        """
        d=self.facets(t,plane,phase=phase,parity=parity)
        q=d['state'][...,:3]; b=d['acc']-self.env.gravity(t,q)
        cos=np.maximum(np.sum(-unit(d['rays']['sun'])*d['normal'],axis=-1),0)
        lit=visible_sun(d['rays']['sun'],d['rays']['earth'],d['rays']['moon'])
        strength=(2*self.env.central_fraction*K.SOLAR_CONSTANT/(K.SPEED_OF_LIGHT*self.env.sigma)
                  *(K.AU/length(d['rays']['sun']))**2*cos**2*lit
                  *(self.layout.clear_side/self.layout.side)**2)
        sail=strength[...,None]*d['normal']
        lam=np.clip(np.sum(b*sail,axis=-1)/np.maximum(strength**2,1e-100),0,1)
        lo=length(b-lam[...,None]*sail)
        hi=np.maximum(length(b),length(b-sail))
        return lo,hi,d

    def separation_certificate(self):
        """Orientation-independent all-time separation for the defined paths.

        Different radial shells: reverse triangle inequality. Same shell:
        slots are two pitches apart; phase_warp derivative supplies a lower
        angular spacing. The derivative bound is independently tested/refined.
        """
        ratio=self.layout.first_radius/(self.layout.first_radius+self.layout.radial_step)
        # Mean-value bounds for asin(ratio*sin(theta)) on |theta|<=.95,
        # and max |quintic taper derivative|=1.875/.30. Outside that interval
        # the warp is identically zero. This covers every radius, not a grid.
        minimum=1-(1-ratio)*((1.875/.30)*np.tan(.95)+1/np.cos(.95)**2)
        distance=2*self.r0*np.sin(minimum*2*np.pi/self.count)
        size=np.sqrt(2)*self.layout.side
        same=float(distance.min()-size-100.)  # two independent 50 m errors
        radial_half=self.layout.side/np.sqrt(2)*np.sin(np.deg2rad(self.layout.normal_limit_deg))
        curvature=(self.layout.side**2/2)/(2*(self.layout.first_radius-radial_half))
        other=self.layout.radial_step-2*(radial_half+curvature)-100.
        return dict(same_shell_surface_clearance_m=same,
            different_shell_surface_clearance_m=other,
            required_clearance_m=100.,position_error_per_tile_m=50.,
            minimum_phase_derivative_bound=minimum,
            normal_envelope_deg=self.layout.normal_limit_deg,
            radial_half_extent_bound_m=radial_half+curvature,
            passes=bool(min(same,other)>=100),
            scope='All dates: different radial shells use the enforced normal envelope; same-shell pairs use enclosing spheres. Failure is inconclusive until an actual intersection is found')


def square_pair_separation(q,a,b,pairs,side=10000.):
    """Separating-axis test of actual, zero-thickness 3-D squares.

    Positive is a separating projection. Nonpositive on every axis indicates
    intersection; no clearance distance is inferred from a failed axis test.
    """
    i,j=pairs.T; delta=q[j]-q[i]
    axes_i=np.stack([a[i],b[i],np.cross(a[i],b[i])],axis=1)
    axes_j=np.stack([a[j],b[j],np.cross(a[j],b[j])],axis=1)
    axes=[axes_i[:,k] for k in range(3)]+[axes_j[:,k] for k in range(3)]
    axes += [np.cross(axes_i[:,k],axes_j[:,l]) for k in range(3) for l in range(3)]
    maximum=np.full(len(pairs),-np.inf)
    for v in axes:
        norm=length(v); u=v/np.maximum(norm,1e-100)[:,None]
        extent=side/2*(np.abs(np.sum(u*a[i],axis=-1))+np.abs(np.sum(u*b[i],axis=-1))+
                       np.abs(np.sum(u*a[j],axis=-1))+np.abs(np.sum(u*b[j],axis=-1)))
        sep=np.abs(np.sum(u*delta,axis=-1))-extent
        maximum=np.maximum(maximum,np.where(norm>1e-10,sep,-np.inf))
    return maximum


def pole_collisions(traffic,t,half_width=40000.):
    """Find actual snapshot intersections in a finite return-pole neighbourhood."""
    js=[];ks=[]
    for j,n in enumerate(traffic.count):
        centre=int(np.round((np.pi-traffic.phase[j]-traffic.omega[j]*t)*n/(2*np.pi)))
        h=int(np.ceil(half_width/(2*np.pi*traffic.r0[j]/n)))+1
        k=np.arange(centre-h,centre+h+1)%n
        js.extend([j]*len(k));ks.extend(k)
    js=np.asarray(js);ks=np.asarray(ks)
    d=traffic.facets(t,js,ks)
    pairs=cKDTree(d['state'][:,:3]).query_pairs(np.sqrt(2)*traffic.layout.side,output_type='ndarray')
    sep=square_pair_separation(d['state'][:,:3],d['a'],d['b'],pairs,traffic.layout.side)
    hit=pairs[sep<-1e-5]
    return dict(time_days=float(t/K.JULIAN_DAY),tiles_examined=len(js),
                actual_square_intersections=len(hit),
                example_ids=[[[int(js[i]),int(ks[i])],[int(js[j]),int(ks[j])]] for i,j in hit[:8]],
                scope='Exact static finite-square intersections in the stated return-pole subset; not an all-fleet pair census')


class RayIndex:
    """Conservative curved ribbon index, followed by exact finite-square tests."""
    def __init__(self,traffic,t,source):
        import shapely
        self.traffic,self.t,self.source=traffic,t,np.asarray(source)
        self.frame=traffic.frames(t); self.sample=traffic.env.at(t)
        self.source_local=self.source@self.frame
        j=np.repeat(np.arange(traffic.layout.planes),2)
        p=np.tile([0,1],traffic.layout.planes)
        # Fixed receiver-domain envelope, including solar parallax and motion.
        r=traffic.r0[j]+p*traffic.layout.radial_step
        extent=traffic.layout.protected_radius+K.SUN_RADIUS*r/length(source)+r*4e4/K.SPEED_OF_LIGHT+30000
        theta=np.arcsin(np.minimum(extent/r,.64))[:,None]*np.linspace(-1,1,33)
        plane=np.broadcast_to(j[:,None],theta.shape); parity=np.broadcast_to(p[:,None],theta.shape)
        state,_,_=traffic.path(t,plane,phase=theta,parity=parity)
        q=state[...,:3]; tau=np.maximum(q@self.frame[:,0],0)/K.SPEED_OF_LIGHT
        q=q-(state[...,3:]+self.sample['moon_v'])*tau[...,None]
        q=q@self.frame; s=self.source_local
        projected=s[1:]+(q[...,1:]-s[1:])*(s[0]/(s[0]-q[...,0]))[...,None]
        # Half diagonal + >= 250 m for curve interpolation, retardation and
        # finite-facet projection magnification in the supported radius range.
        self.buffer=np.sqrt(2)*traffic.layout.side/2+350.
        segments=np.stack([projected[:,:-1],projected[:,1:]],axis=2).reshape(-1,2,2)
        self.tree=shapely.STRtree(shapely.linestrings(segments))
        self.j,self.parity=j,p

    def intercept(self, receiver, return_ids=False, neighbour_slots=1):
        import shapely
        xy=np.asarray(receiver); n=len(xy); hit=np.zeros(n,bool); chosen=np.full((n,2),-1,int)
        points=shapely.points(xy)
        pairs=self.tree.query(points,predicate='dwithin',distance=self.buffer)
        if not pairs.size:return (hit,chosen) if return_ids else hit
        ri,segment=pairs
        keys=np.unique(ri*len(self.j)+segment//32)
        ri,li=keys//len(self.j),keys%len(self.j)
        j=self.j[li]; parity=self.parity[li]
        s=self.source_local; w=np.column_stack([np.zeros(n),xy])
        ray=unit(s-w); wr=w[ri]; dr=ray[ri]
        r=self.traffic.r0[j]+parity*self.traffic.layout.radial_step
        wd=np.sum(wr*dr,axis=-1)
        lam=-wd+np.sqrt(np.maximum(wd*wd+r*r-np.sum(wr*wr,axis=-1),0))
        qr=wr+lam[:,None]*dr
        # Undo the dominant Moon light-time motion to locate the right slot.
        qr+=self.sample['moon_v']@self.frame*(qr[:,0]/K.SPEED_OF_LIGHT)[:,None]
        radial=qr[:,1]*np.cos(self.traffic.alpha[j])+qr[:,2]*np.sin(self.traffic.alpha[j])
        theta=np.arcsin(np.clip(radial/self.traffic.r0[j],-.9,.9))
        near=np.round(((theta-self.traffic.phase[j]-self.traffic.omega[j]*self.t)*self.traffic.count[j]/(2*np.pi)-parity)/2).astype(int)*2+parity
        half=self.traffic.layout.clear_side/2
        for offset in range(-neighbour_slots,neighbour_slots+1):
            use=~hit[ri]
            if not np.any(use):break
            ii=ri[use]; jj=j[use]; kk=(near[use]+2*offset)%self.traffic.count[jj]
            d=self.traffic.facets(self.t,jj,kk)
            q=d['state'][:,:3]; tau=np.maximum(q@self.frame[:,0],0)/K.SPEED_OF_LIGHT
            q=q-(d['state'][:,3:]+self.sample['moon_v'])*tau[:,None]
            q=q@self.frame; normal=d['normal']@self.frame
            aa=d['a']@self.frame; bb=d['b']@self.frame
            denom=np.sum(ray[ii]*normal,axis=-1)
            distance=np.sum((q-w[ii])*normal,axis=-1)/np.where(np.abs(denom)>1e-12,denom,1e-100)
            point=w[ii]+distance[:,None]*ray[ii]-q
            yes=(distance>0)&(distance<length(s-w[ii]))&(np.abs(np.sum(point*aa,axis=-1))<=half)&(np.abs(np.sum(point*bb,axis=-1))<=half)
            hit[ii[yes]]=True; chosen[ii[yes]]=np.column_stack([jj[yes],kk[yes]])
        return (hit,chosen) if return_ids else hit


def sun_source(sample, frame, disk_xy):
    sun,_=arriving_ray_vectors(sample['positions']['sun'],sample['positions']['earth'],
        sample['sun_v'],sample['earth_v'],sample['sun_a'],sample['earth_a'])
    return sun+K.SUN_RADIUS*(disk_xy[0]*frame[:,1]+disk_xy[1]*frame[:,2])
