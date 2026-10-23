"""Axisymmetric finite-volume diffusion around prescribed ice geometries.

Own conservative discretization of the quasistatic diffusion/attachment problem
in Libbrecht (2013, 2015), not a port of the author's cellular automaton. Faces
lie on mesh interfaces; the boundary half-cell resistance is retained. Growth
moves no boundaries. A prescribed open cavity is a geometry sensitivity, not an
inferred crystal profile. Lengths below are normalized by the outer radius.
"""
from __future__ import annotations
from dataclasses import dataclass
import math
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import spsolve
from scipy.optimize import brentq


@dataclass(frozen=True)
class Attachment:
    """Local flux p(s)=alpha(s)*s; all supersaturations are fractions."""
    kind: str
    amplitude: float = 1.
    critical: float = 0.
    mechanism: float = 10.

    def __post_init__(self):
        if self.kind not in ('constant', 'nucleation', 'harrington'):
            raise ValueError('Unknown attachment law')
        if not 0 < self.amplitude <= 1 or not math.isfinite(self.amplitude):
            raise ValueError('Attachment amplitude must lie in (0, 1]')
        if not math.isfinite(self.critical) or self.critical < 0:
            raise ValueError('Nonnegative finite nucleation threshold required')
        if self.kind == 'harrington' and (self.critical <= 0 or not 1 <= self.mechanism <= 15):
            raise ValueError('Harrington law requires positive critical and 1 <= m <= 15')

    def flux(self, s):
        s = np.asarray(s, dtype=float)
        if np.any(s < 0) or not np.all(np.isfinite(s)):
            raise ValueError('Nonnegative finite supersaturation required')
        if self.kind == 'constant' or self.critical == 0:
            return self.amplitude*s, np.full_like(s, self.amplitude)
        positive = s > 0
        p = np.zeros_like(s); derivative = np.zeros_like(s)
        if self.kind == 'nucleation':
            y = self.critical/s[positive]
            alpha = self.amplitude*np.exp(-y)
            p[positive] = alpha*s[positive]
            derivative[positive] = alpha*(1+y)
        else:
            y = np.exp(np.clip(self.mechanism*np.log(self.critical/s[positive]), -700, 700))
            alpha = np.ones_like(y); dp = np.ones_like(y)
            middle = (y >= 1e-5) & (y <= 30)
            large = y > 30
            alpha[middle] = np.tanh(y[middle])/y[middle]
            dp[middle] = (1+self.mechanism)*alpha[middle]-self.mechanism/np.cosh(y[middle])**2
            alpha[large] = 1/y[large]; dp[large] = (1+self.mechanism)/y[large]
            p[positive] = self.amplitude*alpha*s[positive]
            derivative[positive] = self.amplitude*dp
        return p, derivative


def surface_flux(center, distance, kinetic_length_ratio, law):
    """Eliminate the half-cell via s + distance*lambda*p(s) = center."""
    center, distance = np.broadcast_arrays(np.asarray(center, float), np.asarray(distance, float))
    if np.any(center < 0) or np.any(distance <= 0) or not np.all(np.isfinite(kinetic_length_ratio)) or np.any(np.asarray(kinetic_length_ratio) <= 0):
        raise ValueError('Invalid local diffusion inputs')
    lo = np.zeros_like(center); hi = center.copy(); s = .5*center
    k = distance*kinetic_length_ratio
    for _ in range(100):
        p, dp = law.flux(s)
        residual = s+k*p-center
        if np.all(np.abs(residual) <= 2e-14*np.maximum(center, 1e-300)):
            break
        done = np.abs(residual) <= 2e-14*np.maximum(center, 1e-300)
        lo = np.where(residual < 0, s, lo); hi = np.where(residual >= 0, s, hi)
        trial = s-residual/(1+k*dp)
        # A Newton proposal can alternate near opposite bracket ends on a
        # steep attachment transition. Require bracket contraction, and keep
        # individually converged roots fixed while other faces finish.
        safe = (trial > lo+.1*(hi-lo)) & (trial < hi-.1*(hi-lo))
        s = np.where(done, s, np.where(safe, trial, (lo+hi)/2))
    else:
        raise ArithmeticError('Local surface solve failed')
    p, dp = law.flux(s)
    return kinetic_length_ratio*p, kinetic_length_ratio*dp/(1+k*dp), s


def _axis(length, cells_per_radius, extent, breaks=()):
    # Include prescribed cavity boundaries exactly. The maximum inner-cell
    # width is 1/n; narrow rims receive at least four cells.
    cuts = sorted(set([0., float(length), *[float(x) for x in breaks if 0 < x < length]]))
    parts = [np.linspace(a, b, max(4, int(math.ceil((b-a)*cells_per_radius)))+1)[:-1]
             for a, b in zip(cuts[:-1], cuts[1:])]
    values = list(np.concatenate(parts))+[float(length)]
    step = 1/cells_per_radius
    while values[-1] < extent:
        values.append(min(extent, values[-1]+step)); step *= 1.18
    return np.asarray(values)


class Domain:
    """Upper half of a reflected cylinder, open-ended hollow cylinder or sphere.

    Extent is measured in maximum semiaxes. Outer supersaturation is prescribed
    at the finite cylindrical domain boundary; refine extent separately from
    mesh. No far-field humidity tuning or nominal-infinity substitution occurs.
    """
    def __init__(self, aspect=1., resolution=16, extent=12., *, cavity_radius=0., cavity_depth=0., shape='cylinder'):
        if not math.isfinite(aspect) or aspect <= 0 or int(resolution) != resolution or resolution < 4:
            raise ValueError('Positive aspect and integer resolution >= 4 required')
        if not math.isfinite(extent) or extent <= 2:
            raise ValueError('Finite outer extent > 2 required')
        if shape not in ('cylinder', 'sphere') or (shape == 'sphere' and aspect != 1):
            raise ValueError('Supported shapes are cylinders and unit spheres')
        if not 0 <= cavity_radius < 1 or not 0 <= cavity_depth <= aspect:
            raise ValueError('Cavity must lie inside cylinder envelope')
        if bool(cavity_radius) != bool(cavity_depth) or (shape == 'sphere' and cavity_radius):
            raise ValueError('Specify both cavity radius and depth for a cylinder')
        self.aspect = aspect; self.resolution = resolution; self.extent = extent
        self.cavity_radius = cavity_radius; self.cavity_depth = cavity_depth; self.shape = shape
        far = extent*max(1., aspect)
        self.rf = _axis(1., resolution, far, (cavity_radius,))
        self.zf = _axis(aspect, resolution, far, (aspect-cavity_depth,))
        self.r = (self.rf[1:]+self.rf[:-1])/2; self.z = (self.zf[1:]+self.zf[:-1])/2
        rr, zz = np.meshgrid(self.r, self.z, indexing='ij')
        if shape == 'sphere':
            solid = rr*rr+zz*zz < 1
        else:
            solid = (rr < 1) & (zz < aspect)
            solid &= ~((rr < cavity_radius) & (zz > aspect-cavity_depth))
        self.solid = solid
        index = np.full(solid.shape, -1, int); index[~solid] = np.arange(np.sum(~solid)); self.index = index
        n = int(np.sum(~solid)); diag = np.zeros(n); outer = np.zeros(n)
        rows=[]; cols=[]; vals=[]; ice=[]; edges=[]
        nr, nz = solid.shape
        for i, j in zip(*np.where(~solid)):
            p = index[i,j]
            for di, dj in ((-1,0),(1,0),(0,-1),(0,1)):
                k,l = i+di,j+dj
                if k < 0 or l < 0:
                    continue  # r=0 and z=0 reflection symmetry.
                radial = di != 0
                face_r = self.rf[i+(di>0)] if radial else self.r[i]
                face_z = self.zf[j+(dj>0)] if not radial else self.z[j]
                area = (2*math.pi*face_r*(self.zf[j+1]-self.zf[j]) if radial else
                        math.pi*(self.rf[i+1]**2-self.rf[i]**2))
                half = (abs(face_r-self.r[i]) if radial else abs(face_z-self.z[j]))
                if k >= nr or l >= nz:
                    conductance=area/half; diag[p]+=conductance; outer[p]+=conductance
                    edges.append((p, conductance, face_r, face_z))
                elif solid[k,l]:
                    ice.append((p, area, half, 0 if radial else 1, face_r, face_z))
                else:
                    distance=abs(self.r[k]-self.r[i]) if radial else abs(self.z[l]-self.z[j])
                    conductance=area/distance;diag[p]+=conductance
                    rows.append(p);cols.append(index[k,l]);vals.append(-conductance)
        rows.extend(range(n));cols.extend(range(n));vals.extend(diag)
        self.matrix = coo_matrix((vals,(rows,cols)),shape=(n,n)).tocsr()
        self.outer = outer; self.ice = np.asarray(ice); self.edges = np.asarray(edges)
        self.cells = n; self._capacity = None
        self.normal_projection = np.ones(len(self.ice))
        if shape == 'sphere':
            # Stair-step areas overestimate smooth-sphere area even as h -> 0.
            # Project the physical normal flux onto each coordinate face. This
            # is used only for the analytic smooth-sphere control; cylinder and
            # cavity interfaces already coincide with the mesh.
            radius = np.hypot(self.ice[:,4], self.ice[:,5])
            coordinate = np.where(self.ice[:,3] == 0, self.ice[:,4], self.ice[:,5])
            self.normal_projection = coordinate/radius
            # Integrate the known smooth-sphere area exactly. This geometric
            # quadrature normalization vanishes on refinement and is unrelated
            # to attachment, humidity or an observed growth rate.
            self.normal_projection *= 2*math.pi/np.dot(self.ice[:,1],self.normal_projection)

    def dirichlet(self):
        """Unit ice surface, zero finite outer boundary; full-crystal capacity."""
        p = self.ice[:,0].astype(int); g = self.ice[:,1]/self.ice[:,2]
        diagonal=np.bincount(p,weights=g,minlength=self.cells)
        x=spsolve(self.matrix+diags(diagonal),diagonal)
        flux=float(np.dot(g,1-x[p])); outer=float(np.dot(self.outer,x))
        self._capacity=flux
        return dict(half_conductance=flux, capacity_over_radius=2*flux/(4*math.pi),
                    relative_mass_residual=abs(flux-outer)/flux, cells=self.cells)

    def imposed_flux(self, prism=1., basal=0., face_flux=None):
        """Prescribed inward normal flux; returns finite-domain depletion field."""
        p=self.ice[:,0].astype(int)
        law=np.where(self.ice[:,3]==0,prism,basal) if face_flux is None else np.asarray(face_flux,float)
        if law.shape != p.shape or not np.all(np.isfinite(law)):
            raise ValueError('One finite normal flux per surface face required')
        rhs=np.bincount(p,weights=self.ice[:,1]*law,minlength=self.cells)
        x=spsolve(self.matrix,rhs)
        surface=x[p]+self.ice[:,2]*law
        q=float(np.sum(rhs)); far=float(np.dot(self.outer,x))
        return dict(surface=surface, half_flux=q, relative_mass_residual=abs(q-far)/max(abs(q),1e-30))

    def solve(self, supersaturation, kinetic_length_ratio, laws, *, tolerance=1e-10, heat_parameter=0.):
        """Solve local kinetics, optionally with isothermal-crystal latent heat.

        lambda=a*v_mol/(4*Dv). heat_parameter=beta*L*rho_eq*Dv/K, with
        beta=d ln(rho_eq)/dT. Crystal heating is spatially uniform; its external
        heat conductance uses this same geometry and finite outer boundary.
        This excludes substrate conduction, flow, jumps and surface diffusion.
        """
        if not math.isfinite(supersaturation) or supersaturation < 0 or not math.isfinite(heat_parameter) or heat_parameter < 0:
            raise ValueError('Nonnegative finite forcing and heat parameter required')
        if not math.isfinite(kinetic_length_ratio) or kinetic_length_ratio <= 0:
            raise ValueError('Positive finite kinetic length ratio required')
        if len(laws) != 2 or not math.isfinite(tolerance) or tolerance <= 0:
            raise ValueError('Two attachment laws and positive finite tolerance required')
        p=self.ice[:,0].astype(int); area=self.ice[:,1]; distance=self.ice[:,2]
        selections=[self.ice[:,3]==i for i in (0,1)]
        def fixed(effective):
            if effective == 0:
                return np.zeros(self.cells), np.zeros(len(p)), np.zeros(len(p)), 0., 0
            # Solve for depletion u = effective - s, so the far boundary is
            # zero. This avoids subtracting two large outer-boundary terms.
            x=np.zeros(self.cells)
            for iteration in range(70):
                q=np.zeros(len(p));dq=q.copy();surface=q.copy()
                for selected, law in zip(selections,laws):
                    if np.any(selected):
                        q[selected],dq[selected],surface[selected]=surface_flux(effective-x[p[selected]],distance[selected],kinetic_length_ratio*self.normal_projection[selected],law)
                residual=self.matrix@x-np.bincount(p,weights=area*q,minlength=self.cells)
                scale=max(float(np.dot(area,q)),1e-30)
                error=float(np.max(np.abs(residual)))/scale
                if error < tolerance:
                    return x,q,surface,error,iteration
                jacobian=self.matrix+diags(np.bincount(p,weights=area*dq,minlength=self.cells))
                change=spsolve(jacobian,-residual)
                fraction=1.
                if np.any(change < 0):
                    fraction=min(fraction,.99*float(np.min(-x[change<0]/change[change<0])))
                if np.any(change > 0):
                    fraction=min(fraction,.99*float(np.min((effective-x[change>0])/change[change>0])))
                x += fraction*change
                if np.any(x < 0) or np.any(x > effective):
                    raise ArithmeticError('Spatial depletion lies outside the supply range')
            raise ArithmeticError('Spatial nonlinear solve did not converge')
        effective=supersaturation
        if heat_parameter and supersaturation:
            if self._capacity is None: self.dirichlet()
            def balance(value):
                _,q,_,_,_=fixed(value)
                return value+heat_parameter*float(np.dot(area,q))/self._capacity-supersaturation
            effective=brentq(balance,0.,supersaturation,xtol=max(1e-15,supersaturation*1e-11),rtol=1e-11)
        x,q,surface,error,iterations=fixed(effective)
        total=float(np.dot(area,q));far=float(np.dot(self.outer,x))
        if abs(total-far)>max(total,1e-30)*1e-7:
            raise ArithmeticError('Spatial mass ledger failed')
        thermal_residual=(effective+heat_parameter*total/self._capacity-supersaturation if heat_parameter and supersaturation else 0.)
        return dict(half_flux=total, face_flux=q, surface_supersaturation=surface,
            cell_supersaturation=effective-x, effective_far_supersaturation=effective,
            thermal_supersaturation_shift=supersaturation-effective,
            thermal_balance_residual=thermal_residual, nonlinear_relative_residual=error,
            relative_mass_residual=abs(total-far)/max(total,1e-30), iterations=iterations,
            cells=self.cells, outer_extent=self.extent, resolution=self.resolution)
