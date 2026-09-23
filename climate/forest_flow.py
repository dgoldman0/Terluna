"""Spatial neutral mean-flow screen: Oseen advection and porous quadratic drag.

SI inputs, nondimensional velocity v=u/Uref. The prescribed eddy diffusivity is
nu_t/Uref = mixing_length_m. A lateral relaxation reservoir controls the imposed
ambient state. Fourier projection enforces incompressibility. Reflection imposes
free-slip impermeable ground/lid. This is neither a turbulence nor weather model.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
import numpy as np
from scipy.fft import fftn, ifftn, fftfreq
from scipy.sparse.linalg import LinearOperator, gmres


@dataclass(frozen=True)
class Grid:
    nx: int = 40
    ny: int = 32
    nz: int = 20
    lx: float = 2400.
    ly: float = 1920.
    top: float = 1200.
    mixing_length_m: float = 12.
    reservoir_fraction: float = .18
    reservoir_length_m: float = 150.
    tolerance: float = 3e-5
    max_outer: int = 45

    def validate(self):
        for name in ('nx', 'ny', 'nz', 'max_outer'):
            v = getattr(self, name)
            if type(v) is not int or v < 4: raise ValueError('Invalid '+name)
        for name in ('lx','ly','top','mixing_length_m','reservoir_length_m','tolerance'):
            v = getattr(self, name)
            if not np.isfinite(v) or v <= 0: raise ValueError('Invalid '+name)
        if not 0 < self.reservoir_fraction < .4: raise ValueError('Invalid reservoir')

    @property
    def shape(self): return self.nx, self.ny, self.nz
    @property
    def axes(self):
        return ((np.arange(self.nx)+.5)*self.lx/self.nx-self.lx/2,
                (np.arange(self.ny)+.5)*self.ly/self.ny-self.ly/2,
                (np.arange(self.nz)+.5)*self.top/self.nz)
    @property
    def cell_volume(self): return self.lx*self.ly*self.top/(self.nx*self.ny*self.nz)


class SpatialFlow:
    def __init__(self, grid=Grid(), direction_deg=0.):
        grid.validate()
        if not np.isfinite(direction_deg): raise ValueError('Finite wind direction required')
        self.grid = grid
        self.direction = np.array([np.cos(np.deg2rad(direction_deg)),
                                   np.sin(np.deg2rad(direction_deg)), 0.])
        self.shape = (3, grid.nx, grid.ny, 2*grid.nz)
        self.ambient = np.broadcast_to(self.direction[:,None,None,None],self.shape).copy()
        wave, lap = [], []
        for n, length in zip(self.shape[1:], (grid.lx, grid.ly, 2*grid.top)):
            k = 2*np.pi*fftfreq(n, d=length/n); lap.append(k*k)
            k = k.copy()
            if n % 2 == 0: k[n//2] = 0. # real collocated derivative at Nyquist
            wave.append(k)
        self.k = np.array(np.broadcast_arrays(wave[0][:,None,None],wave[1][None,:,None],wave[2][None,None,:]))
        self.k2 = np.sum(self.k*self.k,axis=0)
        self.nonzero = np.where(self.k2>0,self.k2,1.)
        lap2 = lap[0][:,None,None]+lap[1][None,:,None]+lap[2][None,None,:]
        self.linear = grid.mixing_length_m*lap2+1j*np.sum(self.direction[:,None,None,None]*self.k,axis=0)
        x,y,_ = grid.axes
        fx = np.clip((abs(x)/(grid.lx/2)-(1-2*grid.reservoir_fraction))/(2*grid.reservoir_fraction),0,1)
        fy = np.clip((abs(y)/(grid.ly/2)-(1-2*grid.reservoir_fraction))/(2*grid.reservoir_fraction),0,1)
        self.reservoir = np.broadcast_to(np.maximum(fx[:,None],fy[None,:])[:,:,None]**2/grid.reservoir_length_m,self.shape[1:])

    @staticmethod
    def fft(v): return fftn(v,axes=(-3,-2,-1),workers=1)
    @staticmethod
    def inv(v): return ifftn(v,axes=(-3,-2,-1),workers=1).real
    def project(self, vhat):
        return vhat-self.k*np.sum(self.k*vhat,axis=0)/self.nonzero
    def mirror(self, physical):
        if physical.ndim == 3: return np.concatenate((physical,physical[:,:,::-1]),axis=2)
        reflected = physical[:,:,:,::-1].copy(); reflected[2] *= -1
        return np.concatenate((physical,reflected),axis=3)

    def solve(self, drag, reference_speed=10., initial=None, progress=None):
        """drag.coefficient(v,Uref) returns Cd*area/volume in physical cells.

        Nonconvergence raises; returned residuals assess equations actually solved.
        Forces use exactly the same drag map as the momentum sink.
        """
        if not np.isfinite(reference_speed) or reference_speed <= 0:
            raise ValueError('Positive finite reference speed required')
        g = self.grid
        v = self.ambient.copy() if initial is None else self.mirror(np.asarray(initial,float))
        if v.shape != self.shape or not np.all(np.isfinite(v)): raise ValueError('Bad warm start')
        if not np.any(drag.coefficient(self.ambient[:,:,:,:g.nz], reference_speed)):
            physical=self.ambient[:,:,:,:g.nz].copy();forces,ledger=drag.forces(physical,reference_speed)
            return physical,forces,dict(status='EXACT_ZERO_DRAG_CONTROL',grid=asdict(g),reference_speed_m_s=reference_speed,
                direction_vector=self.direction.tolist(),outer_iterations=0,momentum_relative_residual=0.,
                mean_momentum_relative_residual=0.,max_divergence_per_m=0.,reservoir_velocity_departure_fraction=0.,
                max_speed_over_reference=1.,force_ledger=ledger,history=[],limitations=['prescribed_ambient_control'])
        rhs_hat = self.project(self.fft(self.reservoir*self.ambient))
        history=[]; converged=False; n = v.size
        for iteration in range(g.max_outer):
            physical = v[:,:,:,:g.nz]
            a = self.mirror(drag.coefficient(physical,reference_speed))
            beta = self.reservoir+.5*a*np.linalg.norm(v,axis=0)
            mean = float(beta.mean()); denom = self.linear+mean
            rhs = self.inv(rhs_hat/denom).ravel()
            def matvec(flat):
                u=flat.reshape(self.shape)
                return (u+self.inv(self.project(self.fft((beta-mean)*u))/denom)).ravel()
            op=LinearOperator((n,n),matvec=matvec,dtype=float)
            solution, info=gmres(op,rhs,x0=v.ravel(),rtol=min(1e-8,g.tolerance*.001),atol=1e-14,restart=25,maxiter=12)
            if info: raise RuntimeError(f'Oseen linear solve failed: {info}')
            new=solution.reshape(self.shape)
            # Relax only the nonlinear drag iteration. Fixed point solves full equation.
            v=.65*new+.35*v
            v=self.mirror(v[:,:,:,:g.nz])
            anew=self.mirror(drag.coefficient(v[:,:,:,:g.nz],reference_speed))
            friction=.5*anew*np.linalg.norm(v,axis=0)*v
            linear=self.inv(self.linear*self.fft(v))
            residual=linear+self.inv(self.project(self.fft(friction+self.reservoir*(v-self.ambient))))
            scale=max(np.linalg.norm(linear),np.linalg.norm(friction),1e-15)
            rel=float(np.linalg.norm(residual)/scale)
            change=float(np.max(abs(new-v)))
            history.append(dict(iteration=iteration+1,momentum_relative_residual=rel,velocity_change=change))
            if progress is not None: progress(history[-1])
            if rel < g.tolerance and change < 5*g.tolerance: converged=True;break
        if not converged: raise RuntimeError(f'Spatial flow unconverged: residual={rel:.4g}')
        vhat=self.fft(v)
        divergence=self.inv(1j*np.sum(self.k*vhat,axis=0))
        forces,ledger=drag.forces(v[:,:,:,:g.nz],reference_speed)
        demand=friction.mean(axis=(1,2,3))
        supplied=(self.reservoir*(self.ambient-v)).mean(axis=(1,2,3))
        mean_res=float(np.linalg.norm(demand-supplied)/max(np.linalg.norm(demand),1e-14))
        near_reservoir=self.reservoir[:,:,:g.nz]>.75/self.grid.reservoir_length_m
        reservoir_error=float(np.max(np.linalg.norm(v[:,:,:,:g.nz]-self.ambient[:,:,:,:g.nz],axis=0)[near_reservoir]))
        diagnostics=dict(status='CONDITIONAL_STEADY_OSEEN_SCREEN',grid=asdict(g),reference_speed_m_s=reference_speed,
            direction_vector=self.direction.tolist(),outer_iterations=len(history),momentum_relative_residual=rel,
            mean_momentum_relative_residual=mean_res,max_divergence_per_m=float(np.max(abs(divergence))),
            reservoir_velocity_departure_fraction=reservoir_error,max_speed_over_reference=float(np.max(np.linalg.norm(v,axis=0))),
            force_ledger=ledger,history=history,
            limitations=['prescribed_neutral_mixing','fixed_Oseen_advection','free_slip_ground_and_lid',
                         'periodic_box_with_lateral_reservoir','no_gusts_or_weather','upright_geometry_until_damage'])
        return v[:,:,:,:g.nz],forces,diagnostics
