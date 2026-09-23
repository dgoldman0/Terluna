"""Finite-patch, three-dimensional steady Oseen/porous-drag screening model.

Coordinates and velocities are nondimensionalized by Href and upstream Uref.
A MAC grid and a Neumann pressure projection enforce discrete incompressibility.
The streamwise transport speed is the imposed upstream speed (Oseen closure),
not the resolved velocity. Eddy viscosity is prescribed. This is a conditional
mean-flow experiment, not RANS/LES, a gust model, or a lunar weather forecast.
All walls are free slip; upstream/downstream normal velocity is prescribed.
Domain enlargement must test this finite wind-tunnel confinement.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
import math
import hashlib
import numpy as np
from scipy.fft import dctn, idctn


@dataclass(frozen=True)
class FlowConfig:
    nx: int = 72
    ny: int = 48
    nz: int = 24
    length_x_h: float = 12.0
    length_y_h: float = 8.0
    height_h: float = 3.0
    reference_height_m: float = 300.0
    eddy_viscosity_over_uh: float = .02
    tolerance: float = 2e-4
    max_steps: int = 6000
    cfl: float = .65
    advection_order: int = 2

    def validate(self):
        for n in (self.nx, self.ny, self.nz, self.max_steps):
            if not isinstance(n, int) or n < 4:
                raise ValueError('Integer grid sizes/step limit >=4 required')
        if self.advection_order not in (1,2):
            raise ValueError('Advection order must be 1 or 2')
        if not all(np.isfinite(v) for v in asdict(self).values()):
            raise ValueError('Finite flow parameters required')
        if min(self.length_x_h, self.length_y_h, self.height_h,
               self.reference_height_m, self.eddy_viscosity_over_uh,
               self.tolerance, self.cfl) <= 0 or self.cfl > .9:
            raise ValueError('Positive lengths, viscosity, tolerance; 0<CFL<=.9')


class PatchGrid:
    def __init__(self, config=FlowConfig()):
        config.validate()
        self.config = config
        self.shape = (config.nx, config.ny, config.nz)
        self.spacing = np.array([config.length_x_h/config.nx,
                                 config.length_y_h/config.ny,
                                 config.height_h/config.nz])
        # Forest centered horizontally. x is always the upstream wind direction.
        self.centers = tuple((np.arange(n)+.5)*d-(n*d/2 if a<2 else 0)
                             for a,(n,d) in enumerate(zip(self.shape,self.spacing)))
        self.cell_volume_h3 = float(np.prod(self.spacing))
        self.eigenvalues = sum((-4*np.sin(np.pi*np.arange(n)/(2*n))**2/d**2)
            .reshape(tuple(n if j==a else 1 for j in range(3)))
            for a,(n,d) in enumerate(zip(self.shape,self.spacing)))
        self.eigenvalues[0,0,0] = 1.

    def faces(self):
        nx,ny,nz=self.shape
        return [np.ones((nx+1,ny,nz)),np.zeros((nx,ny+1,nz)),np.zeros((nx,ny,nz+1))]

    @staticmethod
    def centered(velocity):
        u,v,w=velocity
        return np.stack(((u[1:]+u[:-1])/2,(v[:,1:]+v[:,:-1])/2,
                         (w[:,:,1:]+w[:,:,:-1])/2),axis=-1)

    def divergence(self, velocity):
        return sum(np.diff(v,axis=a)/self.spacing[a] for a,v in enumerate(velocity))

    def project(self, velocity):
        rhs=self.divergence(velocity)
        compatibility=float(rhs.mean())
        if abs(compatibility)>1e-10:
            raise ValueError('Incompatible imposed boundary flux')
        ph=dctn(rhs,type=2,norm='ortho')/self.eigenvalues
        ph[0,0,0]=0
        pressure=idctn(ph,type=2,norm='ortho')
        answer=[v.copy() for v in velocity]
        for a in range(3):
            section=[slice(None)]*3; section[a]=slice(1,-1)
            answer[a][tuple(section)]-=np.diff(pressure,axis=a)/self.spacing[a]
        return answer, pressure

    def _laplacian(self, a):
        p=np.pad(a,1,mode='edge')
        result=np.zeros_like(a)
        for dim,h in enumerate(self.spacing):
            lo=[slice(1,-1)]*3; hi=lo.copy()
            lo[dim]=slice(0,-2);hi[dim]=slice(2,None)
            result+=(p[tuple(lo)]-2*a+p[tuple(hi)])/h**2
        return result

    @staticmethod
    def _to_faces(a,axis):
        shape=list(a.shape);shape[axis]+=1
        r=np.zeros(shape)
        s=[slice(None)]*3;s[axis]=slice(1,-1)
        lo=[slice(None)]*3;hi=lo.copy();lo[axis]=slice(None,-1);hi[axis]=slice(1,None)
        r[tuple(s)]=(a[tuple(lo)]+a[tuple(hi)])/2
        return r

    @staticmethod
    def _boundary(velocity):
        velocity[0][0]=1.;velocity[0][-1]=1.
        velocity[1][:,0]=0.;velocity[1][:,-1]=0.
        velocity[2][:,:,0]=0.;velocity[2][:,:,-1]=0.

    def solve(self, drag_area_coefficient, initial=None):
        """Solve with dimensionless Cd*frontal-area-density*Href (no hidden LAI).

        Fixed crown shape/drag during one solve. Reconfiguration, when used,
        must update this field externally. A nonconverged solve is an error.
        """
        b=np.asarray(drag_area_coefficient,float)
        if b.shape!=self.shape or not np.all(np.isfinite(b)) or np.any(b<0):
            raise ValueError('Nonnegative finite drag on the complete grid required')
        c=self.config
        velocity=self.faces() if initial is None else [v.copy() for v in initial]
        if [v.shape for v in velocity] != [v.shape for v in self.faces()] or not all(np.all(np.isfinite(v)) for v in velocity):
            raise ValueError('Invalid initial velocity')
        self._boundary(velocity)
        velocity,_=self.project(velocity)
        nu=c.eddy_viscosity_over_uh
        def advection(v):
            if c.advection_order==1:
                return np.diff(v,axis=0,prepend=v[:1])/self.spacing[0]
            p=np.concatenate((v[:1],v[:1],v),axis=0)
            return (3*p[2:]-4*p[1:-1]+p[:-2])/(2*self.spacing[0])
        residual=math.inf
        pressure=np.zeros(self.shape)
        for step in range(1,c.max_steps+1):
            centered=self.centered(velocity)
            speed=np.linalg.norm(centered,axis=-1)
            drag=-.5*b[...,None]*speed[...,None]*centered
            dt=min(c.cfl/(2/self.spacing[0]+2*nu*np.sum(1/self.spacing**2)
                      +max(float(np.max(b*speed)),1e-9)),nu)
            predicted=[]
            for a,v in enumerate(velocity):
                # Positive-x Oseen transport. Domain refinement is required.
                adv=advection(v)
                force=self._to_faces(drag[...,a],a)
                predicted.append(v+dt*(-adv+nu*self._laplacian(v)+force))
            self._boundary(predicted)
            projected, phi=self.project(predicted)
            residual=max(float(np.max(abs(n-o)))/dt for n,o in zip(projected,velocity))
            velocity=projected;pressure=phi/dt
            if not np.isfinite(residual) or max(np.max(abs(v)) for v in velocity)>20:
                raise RuntimeError('Unbounded Oseen iterate')
            if residual<c.tolerance:
                break
        if residual>=c.tolerance:
            raise RuntimeError(f'Flow failed to converge: residual={residual:g}, steps={step}')
        centered=self.centered(velocity)
        speed=np.linalg.norm(centered,axis=-1)
        drag=.5*b[...,None]*speed[...,None]*centered
        integrated=drag.sum(axis=(0,1,2))*self.cell_volume_h3
        div=float(np.max(abs(self.divergence(velocity))))
        # Exact summed x-momentum balance of the discrete face equation. The
        # dual volumes omit the two prescribed x-boundary planes.
        u=velocity[0]
        transport=-advection(u)
        diffuse=nu*self._laplacian(u)
        force=-self._to_faces(drag[...,0],0)
        gradp=np.zeros_like(u);gradp[1:-1]=np.diff(pressure,axis=0)/self.spacing[0]
        ledger={name:float(a[1:-1].sum()*self.cell_volume_h3)
                for name,a in [('transport',transport),('diffusion',diffuse),
                               ('vegetation',force),('pressure',-gradp)]}
        ledger['residual']=sum(ledger.values())
        denom=max(abs(ledger['vegetation']),1e-12)
        diagnostic=dict(status='CONVERGED_CONDITIONAL_3D_OSEEN_MEAN_FLOW',steps=step,
            drag_coefficient_sha256=hashlib.sha256(b.astype('<f8').tobytes()).hexdigest(),
            max_projected_acceleration=residual,max_divergence=div,
            discrete_momentum=ledger,momentum_relative_residual=abs(ledger['residual'])/denom,
            integrated_drag_h2=integrated.tolist(),maximum_speed_ratio=float(speed.max()),
            max_vertical_speed_ratio=float(np.max(abs(centered[...,2]))),
            streamwise_upwind_numerical_viscosity_over_uh=float(self.spacing[0]/2) if c.advection_order==1 else 0.0,
            config=asdict(c),flow_model='Oseen transport; prescribed eddy viscosity; quadratic porous drag',
            density_closure='constant reference mass density; buoyancy omitted',
            boundary='fixed normal upstream/outlet flow; impermeable free-slip sides/top/ground',
            unresolved=['gusts','turbulence production','mean-flow self-advection','atmospheric stability','lunar circulation'])
        return FlowResult(self,velocity,centered,pressure,diagnostic)


@dataclass
class FlowResult:
    grid: PatchGrid
    velocity: list
    centered: np.ndarray
    pressure: np.ndarray
    diagnostics: dict
