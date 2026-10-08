"""Release the prescribed fleet into actual-date gravitational trajectories.

Ideal cancellation of the incident reflected-band photon force is the only
translation control in this diagnostic. It is NOT a gap-repair controller.
Smooth phase/plane interpolation is validated against separately propagated
integer members before interpreting finite-tile geometry.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline, make_interp_spline

from shared import constants as K
from .active_global import Traffic
from .optical import length,unit


class NaturalTraffic(Traffic):
    def separation_certificate(self):
        raise ValueError('Prescribed-shell bounds do not apply to propagated natural traffic')

    def initial_states(self,plane,phase,parity):
        state,_,tangent=super().path(0,plane,phase=phase,parity=parity)
        # Preserve tile positions, but initialize actual lunar circular speeds.
        # The rejected Sun-following schedule's differentiated phase warp and
        # frame rotation are not a natural orbital velocity field.
        state[...,3:]=(np.sqrt(K.MOON_GM/length(state[...,:3]))*np.sign(self.omega[np.asarray(plane,dtype=int)]))[...,None]*unit(tangent)
        return state

    def propagate(self,days=7.,plane_stride=16,phase_nodes=128,max_step=1200.):
        n=self.layout.planes
        # Split both the radius-order cusp and the orbit-sense discontinuity.
        self.segments=[np.unique(np.r_[np.arange(0,n//2-1,plane_stride),n//2-1]),
                       np.array([n//2]),
                       np.unique(np.r_[np.arange(n//2+1,n-1,plane_stride),n-1])]
        self.grid_j=np.concatenate(self.segments).astype(int)
        self.phase_nodes=phase_nodes
        self.phases=np.arange(phase_nodes)*2*np.pi/phase_nodes
        j,par,ph=np.meshgrid(self.grid_j,[0,1],self.phases,indexing='ij')
        self.grid_shape=j.shape
        initial=self.initial_states(j,ph,par)
        def rhs(t,flat):
            y=flat.reshape(-1,6)
            return np.column_stack([y[:,3:],self.env.gravity(t,y[:,:3])]).ravel()
        self.times=np.arange(0,days*K.JULIAN_DAY+1,1800.)
        sol=solve_ivp(rhs,[0,self.times[-1]],initial.ravel(),t_eval=self.times,
            method='DOP853',rtol=2e-11,atol=np.tile([1e-4]*3+[1e-8]*3,initial.size//6),max_step=max_step)
        if not sol.success:raise RuntimeError(sol.message)
        self.states=sol.y.T.reshape(len(self.times),*self.grid_shape,6)
        self.nfev=sol.nfev;self.cached_time=None

    def cache(self,t):
        if self.cached_time==float(t):return
        if t<self.times[0] or t>self.times[-1]:raise ValueError('Time outside natural propagation')
        i=min(np.searchsorted(self.times,t,side='right')-1,len(self.times)-2)
        h=self.times[i+1]-self.times[i];z=(t-self.times[i])/h
        a=self.states[i];b=self.states[i+1]
        q=(2*z**3-3*z*z+1)*a[...,:3]+(z**3-2*z*z+z)*h*a[...,3:]+(-2*z**3+3*z*z)*b[...,:3]+(z**3-z*z)*h*b[...,3:]
        v=(6*z*z-6*z)/h*a[...,:3]+(3*z*z-4*z+1)*a[...,3:]+(-6*z*z+6*z)/h*b[...,:3]+(3*z*z-2*z)*b[...,3:]
        values=np.concatenate([q,v],axis=-1)
        # Remove the fast Kepler rotation from the VECTOR coordinates before
        # interpolating planes. Keep initial phase as the trajectory identity;
        # shifting the phase argument would move narrow initial features across
        # the interpolation axis. Every input is still an IVP state at date t.
        alpha=self.alpha[self.grid_j];f0=self.frames(0)
        m=np.cos(alpha)[:,None]*f0[:,1]+np.sin(alpha)[:,None]*f0[:,2]
        h=np.cross(f0[:,0],m)
        basis=np.stack([np.broadcast_to(f0[:,0],m.shape),m,h],axis=-1)
        local=np.einsum('jpkvi,jil->jpkvl',values.reshape(*values.shape[:-1],2,3),basis)
        angle=self.omega[self.grid_j]*t;c=np.cos(angle)[:,None,None,None];s=np.sin(angle)[:,None,None,None]
        x=local[...,0].copy();y=local[...,1].copy()
        local[...,0]=c*x+s*y;local[...,1]=-s*x+c*y
        values=local.reshape(values.shape)
        full=np.empty((self.layout.planes,2,self.phase_nodes,6));first=0
        for segment in self.segments:
            part=values[first:first+len(segment)];first+=len(segment)
            jj=np.arange(segment[0],segment[-1]+1)
            degree=min(5,len(segment)-1)
            full[jj]=make_interp_spline(segment,part,k=degree,axis=0)(jj) if len(segment)>1 else part
        ph=np.r_[self.phases,2*np.pi]
        periodic=np.concatenate([full,full[:,:,:1]],axis=2)
        self.coefficients=CubicSpline(ph,periodic,axis=2,bc_type='periodic').c
        self.cached_time=float(t)

    def path(self,t,plane,slot=None,phase=None,parity=None):
        self.cache(t);j=np.asarray(plane,dtype=int)
        if slot is not None:
            k=np.asarray(slot,dtype=int)%self.count[j];parity=k%2
            ph=2*np.pi*k/self.count[j]+self.phase[j]
        else:
            # API compatibility: explicit phase denotes the original clock at
            # date t. Its initial value identifies the propagated streamline.
            ph=np.asarray(phase)-self.omega[j]*t
        ph=ph%(2*np.pi);idx=np.floor(ph*self.phase_nodes/(2*np.pi)).astype(int)%self.phase_nodes
        z=ph-self.phases[idx];p=np.asarray(parity,dtype=int)
        c=self.coefficients[:,idx,j,p,:]
        state=((c[0]*z[...,None]+c[1])*z[...,None]+c[2])*z[...,None]+c[3]
        tangent=(3*c[0,...,:3]*z[...,None]+2*c[1,...,:3])*z[...,None]+c[2,...,:3]
        f0=self.frames(0);m=np.cos(self.alpha[j])[...,None]*f0[:,1]+np.sin(self.alpha[j])[...,None]*f0[:,2]
        h=np.cross(f0[:,0],m);angle=self.omega[j]*t
        u=np.cos(angle)[...,None]*f0[:,0]+np.sin(angle)[...,None]*m
        v=-np.sin(angle)[...,None]*f0[:,0]+np.cos(angle)[...,None]*m
        basis=np.stack([u,v,h],axis=-1)
        state=np.einsum('...vi,...li->...vl',state.reshape(*state.shape[:-1],2,3),basis).reshape(state.shape)
        tangent=unit(np.einsum('...i,...li->...l',tangent,basis))
        return state,self.env.gravity(t,state[...,:3]),tangent


class NaturalRayIndex:
    """Index propagated curve segments, then test actual integer square facets."""
    def __init__(self,fleet,t,source):
        import shapely
        self.fleet,self.t,self.source=fleet,t,source
        self.frame=fleet.frames(t);self.sample=fleet.env.at(t);self.s=source@self.frame
        phase=np.linspace(0,2*np.pi,257)
        j,p,ph=np.meshgrid(np.arange(fleet.layout.planes),[0,1],phase,indexing='ij')
        state,_,_=fleet.path(t,j,phase=ph+fleet.omega[j]*t,parity=p)
        q=state[...,:3];tau=np.maximum(q@self.frame[:,0],0)/K.SPEED_OF_LIGHT
        q=(q-(state[...,3:]+self.sample['moon_v'])*tau[...,None])@self.frame
        z=self.s[1:]+(q[...,1:]-self.s[1:])*(self.s[0]/(self.s[0]-q[...,0]))[...,None]
        seg=np.stack([z[:,:,:-1],z[:,:,1:]],axis=-2).reshape(-1,2,2)
        keep=((q[:,:,:-1,0]>0)&(q[:,:,1:,0]>0)).ravel()
        R=fleet.layout.protected_radius+15000.
        keep&=np.all(np.min(seg,axis=1)<=R,axis=1)&np.all(np.max(seg,axis=1)>=-R,axis=1)
        self.segment=seg[keep];self.labels=np.flatnonzero(keep)
        self.tree=shapely.STRtree(shapely.linestrings(self.segment))

    def intercept(self,xy):
        import shapely
        xy=np.asarray(xy);hit=np.zeros(len(xy),bool)
        pairs=self.tree.query(shapely.points(xy),predicate='dwithin',distance=self.fleet.layout.side/np.sqrt(2)+1500.)
        ri,si=pairs
        if not len(ri):return hit
        label=self.labels[si];j=label//512;par=(label//256)%2;step=label%256
        segment=self.segment[si];d=segment[:,1]-segment[:,0]
        fraction=np.clip(np.sum((xy[ri]-segment[:,0])*d,axis=1)/np.maximum(np.sum(d*d,axis=1),1e-100),0,1)
        ph=(step+fraction)*2*np.pi/256
        near=np.round(((ph-self.fleet.phase[j])*self.fleet.count[j]/(2*np.pi)-par)/2).astype(int)*2+par
        w=np.column_stack([np.zeros(len(xy)),xy]);ray=unit(self.s-w)
        for offset in (-1,0,1):
            take=~hit[ri]
            if not np.any(take):break
            ii=ri[take];jj=j[take];kk=(near[take]+2*offset)%self.fleet.count[jj]
            d=self.fleet.facets(self.t,jj,kk);q=d['state'][:,:3]
            tau=np.maximum(q@self.frame[:,0],0)/K.SPEED_OF_LIGHT
            q=(q-(d['state'][:,3:]+self.sample['moon_v'])*tau[:,None])@self.frame
            normal=d['normal']@self.frame;aa=d['a']@self.frame;bb=d['b']@self.frame
            den=np.sum(ray[ii]*normal,axis=1)
            distance=np.sum((q-w[ii])*normal,axis=1)/np.where(np.abs(den)>1e-12,den,1e-100)
            local=w[ii]+distance[:,None]*ray[ii]-q
            yes=(distance>0)&(distance<length(self.s-w[ii]))&(np.abs(np.sum(local*aa,axis=1))<=self.fleet.layout.clear_side/2)&(np.abs(np.sum(local*bb,axis=1))<=self.fleet.layout.clear_side/2)
            hit[ii[yes]]=True
        return hit
