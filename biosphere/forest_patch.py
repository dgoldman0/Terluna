"""Explicit forest geometry, matched-mass cases and conservative porous drag."""
from __future__ import annotations
from dataclasses import dataclass,replace
import math
import numpy as np
from scipy import sparse
from scipy.optimize import brentq
from biosphere.megaforest_mechanics import TreeConfig,RootConfig,StaticTree


@dataclass(frozen=True)
class Individual:
    key: str
    x: float
    y: float
    tree: TreeConfig
    roots: RootConfig
    alive: bool = True
    crown_remaining: float = 1.

    def effective(self):
        if not 0 < self.crown_remaining <= 1: raise ValueError('Invalid crown fraction')
        return replace(self.tree,crown_radius_m=self.tree.crown_radius_m*math.sqrt(self.crown_remaining),
                       crown_wood_to_stem_mass=self.tree.crown_wood_to_stem_mass*self.crown_remaining)


def biomass_proxy(tree):
    """Aboveground specified tissue/substrate mass; water and roots excluded."""
    f=tree.tip_diameter_fraction
    stem=tree.wood_density_kg_m3*math.pi*tree.base_diameter_m**2/4*tree.height_m*(1+f+f*f)/3
    return stem*(1+tree.crown_wood_to_stem_mass)+math.pi*tree.crown_radius_m**2*(
        tree.leaf_area_index*tree.leaf_mass_kg_m2_leaf+tree.epiphyte_mass_kg_m2_footprint)


def make_patch(shape='flat',nx=7,ny=7,spacing=90.,height=300.,matched_mass=True,elements=20,seed=2048):
    if shape not in ('flat','tapered','dome','ramp','irregular','gap'): raise ValueError('Unknown shape')
    if type(nx) is not int or type(ny) is not int or min(nx,ny)<2 or min(spacing,height)<=0: raise ValueError('Bad patch')
    rng=np.random.default_rng(seed); entries=[]
    for i in range(nx):
        for j in range(ny):
            x=(i-(nx-1)/2)*spacing;y=(j-(ny-1)/2)*spacing
            edge=min(i,j,nx-1-i,ny-1-j)
            h=1.
            if shape=='tapered': h=.55+.45*min(edge/2,1.)
            if shape=='dome':
                r2=(2*i/(nx-1)-1)**2+(2*j/(ny-1)-1)**2
                h=.55+.45*math.sqrt(max(0,1-r2/2))
            if shape=='ramp': h=.55+.45*i/(nx-1)
            if shape=='irregular': h=.6+.4*rng.random()
            if shape=='gap' and abs(x)<spacing*.6 and abs(y)<spacing*1.6: continue
            entries.append((f't{i:02d}-{j:02d}',x,y,h))
    if shape=='irregular':
        peak=max(e[3] for e in entries);entries=[(k,x,y,h/peak) for k,x,y,h in entries]
    def config(h):
        H=height*h
        return replace(TreeConfig(),height_m=H,base_diameter_m=H/35,crown_radius_m=.18*H,elements=elements)
    target=nx*ny*biomass_proxy(config(1))
    scale=brentq(lambda a:sum(biomass_proxy(config(e[3]*a)) for e in entries)-target,.2,4.) if matched_mass else 1.
    forest=[]
    for key,x,y,h in entries:
        tree=config(h*scale)
        roots=replace(RootConfig(),plate_radius_m=.1*tree.height_m,plate_depth_m=4*math.sqrt(tree.height_m/300))
        forest.append(Individual(key,x,y,tree,roots))
    return forest,dict(shape=shape,matched_mass=matched_mass,aboveground_mass_proxy_kg=sum(biomass_proxy(f.tree) for f in forest),
        target_mass_proxy_kg=target,matched_height_scale=scale,planting_area_m2=nx*ny*spacing**2,
        maximum_height_m=max(f.tree.height_m for f in forest),minimum_height_m=min(f.tree.height_m for f in forest),
        individual_count=len(forest),root_and_water_costs_excluded_from_matching=True)


def _horizontal_weights(x,y,cx,cy,radius):
    dx=x[1]-x[0];dy=y[1]-y[0]
    # Continuous triangular footprints for trunks; smooth porous crown footprints.
    if radius<.5*min(dx,dy):
        wx=np.maximum(0,1-abs(x-cx)/dx);wy=np.maximum(0,1-abs(y-cy)/dy);w=wx[:,None]*wy
    else:
        q=((x[:,None]-cx)/radius)**2+((y[None,:]-cy)/radius)**2
        w=np.where(q<2.25,np.exp(-2*q),0.)
    if w.sum()<=0: raise ValueError('Drag footprint outside domain')
    i,j=np.nonzero(w);return i,j,w[i,j]/w.sum()


class TreeDrag:
    """Same Cd-area quadrature for fluid momentum extraction and individual load."""
    def __init__(self,forest,grid,air):
        grid.validate();self.grid=grid;self.forest=forest
        if len({f.key for f in forest})!=len(forest): raise ValueError('Duplicate tree identifiers')
        x,y,z=grid.axes;dz=grid.top/grid.nz
        row=[];col=[];data=[];area=[];owner=[];kind=[];heights=[];exponent=[];starts=[];floors=[]
        gauss,gw=np.polynomial.legendre.leggauss(6)
        for n,f in enumerate(forest):
            f.tree.validate();f.roots.validate()
            if not f.alive: continue
            t=f.effective()
            if t.height_m>=grid.top or abs(f.x)+1.5*t.crown_radius_m>grid.lx*(.5-grid.reservoir_fraction) or abs(f.y)+1.5*t.crown_radius_m>grid.ly*(.5-grid.reservoir_fraction):
                raise ValueError('Tree intrudes into lid or lateral reservoir')
            for k in range(grid.nz):
                lo=k*dz;hi=min((k+1)*dz,t.height_m)
                if hi<=lo: continue
                zs=(lo+hi)/2+(hi-lo)/2*gauss
                s=(zs/t.height_m-t.crown_base_fraction)/(1-t.crown_base_fraction)
                shape=np.where((s>0)&(s<1),np.sqrt(np.maximum(0,4*s*(1-s))),0)
                widths=[t.base_diameter_m*(1-(1-t.tip_diameter_fraction)*zs/t.height_m)*t.stem_drag_coefficient,
                        2*t.crown_radius_m*shape*t.crown_frontal_fraction*t.crown_drag_coefficient]
                for typ,width in enumerate(widths):
                    a=float((width@gw)*(hi-lo)/2)
                    if a<=0: continue
                    h=float((width*zs)@gw/(width@gw))
                    s0=(h/t.height_m-t.crown_base_fraction)/(1-t.crown_base_fraction)
                    radius=t.crown_radius_m*math.sqrt(max(0,4*s0*(1-s0))) if typ else 0.
                    ii,jj,ww=_horizontal_weights(x,y,f.x,f.y,radius)
                    r=len(area);row.extend([r]*len(ii));col.extend(((ii*grid.ny+jj)*grid.nz+k).tolist());data.extend(ww)
                    area.append(a);owner.append(n);kind.append(typ);heights.append(h)
                    exponent.append(t.vogel_exponent if typ else 0.);starts.append(t.reconfiguration_start_m_s);floors.append(t.minimum_drag_fraction)
        self.W=sparse.csr_matrix((data,(row,col)),shape=(len(area),np.prod(grid.shape)))
        self.area=np.array(area);self.owner=np.array(owner,int);self.kind=np.array(kind,int);self.heights=np.array(heights)
        self.exponent=np.array(exponent);self.starts=np.array(starts);self.floors=np.array(floors)
        self.rho_z=air.sample(z)['density_kg_m3'];self.rho0=float(air.sample(0)['density_kg_m3'])
        self.rho=np.broadcast_to(self.rho_z,grid.shape).ravel()
        self.rigid=bool(np.all(self.exponent==0))
        self._fixed=np.asarray(self.W.T@self.area).reshape(grid.shape)/grid.cell_volume

    def response(self,v,speed):
        if self.rigid: return np.ones(len(self.area))
        rms=np.sqrt(np.maximum(0,self.W@(np.sum(v*v,axis=0).ravel())))*speed
        return np.maximum(self.floors,np.maximum(1,rms/self.starts)**self.exponent)

    def coefficient(self,v,speed):
        a=self._fixed if self.rigid else np.asarray(self.W.T@(self.area*self.response(v,speed))).reshape(self.grid.shape)/self.grid.cell_volume
        return a*self.rho.reshape(self.grid.shape)/self.rho0

    def forces(self,v,speed):
        drag=(np.linalg.norm(v,axis=0)[None]*v).reshape(3,-1)
        r=self.response(v,speed)
        forces=.5*speed**2*(self.W@(drag*self.rho).T)*(self.area*r)[:,None]
        fluid=(.5*self.rho0*speed**2*self.coefficient(v,speed)[None]*np.linalg.norm(v,axis=0)[None]*v).sum(axis=(1,2,3))*self.grid.cell_volume
        trees=forces.sum(axis=0)
        err=float(np.linalg.norm(fluid-trees)/max(np.linalg.norm(fluid),1e-15))
        return forces,dict(tree_force_n=trees.tolist(),fluid_sink_n=fluid.tolist(),relative_residual=err,
            summed_cd_area_m2=float(self.area.sum()),density_variation_fraction=float(np.ptp(self.rho_z)/self.rho0))
