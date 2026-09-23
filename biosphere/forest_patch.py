"""Individual-tree geometry, conserved drag deposition and forest experiments.

One list of trees supplies mechanical masses, projected drag, and contacts.
Voxelization preserves each tree's integrated projected area; numerical grid
convergence is still required to establish its spatial accuracy. Allometry and
patch shapes are hypotheses. Equal initial standing mass is a budget control,
not equal photosynthetic productivity or evidence of evolutionary selection.
"""
from __future__ import annotations
from dataclasses import dataclass, replace, asdict
import math
import hashlib
import numpy as np
from scipy.optimize import brentq
from biosphere.megaforest_mechanics import TreeConfig, RootConfig, StaticTree


@dataclass(frozen=True)
class PatchTree:
    identifier: str
    x_m: float
    y_m: float
    tree: TreeConfig
    roots: RootConfig
    row: int = 0
    column: int = 0
    remaining_crown: float = 1.
    alive: bool = True

    def effective_tree(self):
        """Crown-shedding hypothesis: radius contracts, mass per area preserved."""
        f=self.remaining_crown
        if not np.isfinite(f) or not 0<=f<=1:
            raise ValueError('Crown fraction in [0,1] required')
        t=self.tree
        if f==1: return t
        # f=0 retains a numerically nonzero crown radius with zero mass/drag.
        return replace(t,crown_radius_m=t.crown_radius_m*math.sqrt(max(f,1e-12)),
            crown_frontal_fraction=t.crown_frontal_fraction if f else 0.,
            crown_wood_to_stem_mass=t.crown_wood_to_stem_mass*f,
            leaf_area_index=t.leaf_area_index if f else 0.,
            epiphyte_mass_kg_m2_footprint=t.epiphyte_mass_kg_m2_footprint if f else 0.,
            retained_water_mm=t.retained_water_mm if f else 0.)


def make_patch(shape='flat', nrows=7, ncols=7, height_m=300., spacing_m=105.,
               direction_deg=0., matched_mass=True, air=None, soft_soil=True):
    from atmosphere.lower_air import LowerAir
    air=air or LowerAir()
    if shape not in ('flat','dome','ramp','irregular','gap'):
        raise ValueError('Unknown patch shape')
    if min(nrows,ncols)<2 or min(height_m,spacing_m)<=0 or not np.isfinite(direction_deg):
        raise ValueError('Invalid patch geometry')
    # Budget matching includes dry + wet standing mass; report foliar area apart.
    def build(scale,which):
        trees=[]
        for i in range(nrows):
            for j in range(ncols):
                a=(i-(nrows-1)/2)/max((nrows-1)/2,1)
                b=(j-(ncols-1)/2)/max((ncols-1)/2,1)
                if which=='gap' and abs(a)<.35 and abs(b)<.35: continue
                factor={'flat':1.,'gap':1.,'dome':.6+.4*max(0,1-max(abs(a),abs(b))**2),
                        'ramp':.6+.4*(a+1)/2,
                        'irregular':.8+.2*math.cos(2.1*i+1.3*j)}[which]
                h=height_m*scale*factor
                t=replace(TreeConfig(),height_m=h,base_diameter_m=h/35,
                    crown_radius_m=.18*h,elements=32)
                r=replace(RootConfig(),plate_radius_m=.1*h,plate_depth_m=4*math.sqrt(h/300))
                if soft_soil:r=replace(r,cohesion_pa=2000.,pore_pressure_ratio=.8)
                # Rotate forest, leaving the wind-tunnel inlet convention fixed.
                theta=math.radians(direction_deg)
                x=(i-(nrows-1)/2)*spacing_m;y=(j-(ncols-1)/2)*spacing_m
                trees.append(PatchTree(f'r{i:02d}c{j:02d}',x*math.cos(theta)-y*math.sin(theta),
                    x*math.sin(theta)+y*math.cos(theta),t,r,i,j))
        return trees
    def mass(trees):
        return sum(tree_mass(p.tree) for p in trees)
    target=mass(build(1,'flat'))
    scale=brentq(lambda s:mass(build(s,shape))-target,.3,3.) if matched_mass and shape!='flat' else 1.
    trees=build(scale,shape)
    return trees,dict(shape=shape,trees=len(trees),direction_deg=direction_deg,
        nrows=nrows,ncols=ncols,spacing_m=spacing_m,shape_height_scale=scale,
        reference_flat_mass_kg=target,initial_mass_kg=mass(trees),matched_mass=matched_mass,
        min_height_m=min(p.tree.height_m for p in trees),max_height_m=max(p.tree.height_m for p in trees),
        leaf_area_m2=sum(math.pi*p.tree.crown_radius_m**2*p.tree.leaf_area_index for p in trees),
        foliar_area_or_productivity_matched=False,evolution_simulated=False)


def tree_mass(t):
    # Exact tapered-frustum volume; prior beam masses use midpoint quadrature.
    stem=t.wood_density_kg_m3*math.pi*t.height_m*t.base_diameter_m**2*(1+t.tip_diameter_fraction+t.tip_diameter_fraction**2)/12
    crown=t.crown_wood_to_stem_mass*stem+math.pi*t.crown_radius_m**2*(
        t.leaf_area_index*t.leaf_mass_kg_m2_leaf+t.epiphyte_mass_kg_m2_footprint+t.retained_water_mm)
    return stem+crown


def crown_radius_at(tree,z):
    s=(z/tree.height_m-tree.crown_base_fraction)/(1-tree.crown_base_fraction)
    return tree.crown_radius_m*math.sqrt(max(0,4*s*(1-s))) if 0<s<1 else 0.


class PatchDrag:
    """Sparse per-tree voxel support. Cd*area shares exactly partition total drag."""
    def __init__(self, grid, trees):
        trees=tuple(trees)
        self.grid=grid;self.trees=trees;self.entries=[]
        self.coefficient=np.zeros(grid.shape)
        if len({p.identifier for p in trees})!=len(trees):
            raise ValueError('Unique tree identifiers required')
        H=grid.config.reference_height_m
        x,y,z=[a*H for a in grid.centers]
        dx,dy,dz=grid.spacing*H
        dvol=dx*dy*dz
        for p in trees:
            if not p.alive:continue
            if not np.isfinite(p.x_m+p.y_m): raise ValueError('Finite tree position required')
            t=p.effective_tree();t.validate();p.roots.validate()
            if t.vogel_exponent!=0:
                raise ValueError('Patch drag currently requires fixed Cd; streamlining is a future coupled closure')
            # Reject clipped geometry rather than renormalizing lost portions.
            if (p.x_m-t.crown_radius_m < x[0]-dx/2 or p.x_m+t.crown_radius_m>x[-1]+dx/2
                or p.y_m-t.crown_radius_m<y[0]-dy/2 or p.y_m+t.crown_radius_m>y[-1]+dy/2
                or t.height_m>z[-1]+dz/2):
                raise ValueError('Tree geometry extends outside the flow domain')
            idxs=[];vals=[];crowns=[];heights=[]
            # 8-point Gauss per layer for crown area. A narrow stem is deposited
            # with a normalized bilinear kernel, never discarded by cell misses.
            gx,gw=np.polynomial.legendre.leggauss(8)
            target=0.
            for k,zc in enumerate(z):
                lo=max(0,zc-dz/2);hi=min(t.height_m,zc+dz/2)
                if hi<=lo:continue
                zz=(lo+hi)/2+gx*(hi-lo)/2
                diam=t.base_diameter_m*(1-(1-t.tip_diameter_fraction)*zz/t.height_m)
                stemarea=float(np.dot(gw,diam)*(hi-lo)/2)
                # Split at crown base to retain area at coarse resolutions.
                cl=max(lo,t.crown_base_fraction*t.height_m)
                crownarea=0.
                if hi>cl:
                    zs=(cl+hi)/2+gx*(hi-cl)/2
                    width=np.array([2*crown_radius_at(t,h)*t.crown_frontal_fraction for h in zs])
                    crownarea=float(np.dot(gw,width)*(hi-cl)/2)
                rad=max(crown_radius_at(t,min(zc,.999*t.height_m)),.45*min(dx,dy))
                ix=np.flatnonzero(abs(x-p.x_m)<max(rad,dx)*1.4)
                iy=np.flatnonzero(abs(y-p.y_m)<max(rad,dy)*1.4)
                if not len(ix) or not len(iy): raise ValueError('Grid failed to cover tree support')
                xx,yy=np.meshgrid(x[ix]-p.x_m,y[iy]-p.y_m,indexing='ij')
                # Smooth compact footprint; explicit normalization preserves Cd A.
                weights=np.maximum(0,1-(xx/max(rad,dx))**2-(yy/max(rad,dy))**2)
                sw=np.maximum(0,1-abs(xx)/dx)*np.maximum(0,1-abs(yy)/dy)
                if weights.sum()==0: weights=sw.copy()
                weights/=weights.sum();sw/=sw.sum()
                area=t.crown_drag_coefficient*crownarea*weights+t.stem_drag_coefficient*stemarea*sw
                ca=t.crown_drag_coefficient*crownarea*weights
                ii,jj=np.meshgrid(ix,iy,indexing='ij')
                ids=np.ravel_multi_index((ii.ravel(),jj.ravel(),np.full(ii.size,k)),grid.shape)
                keep=area.ravel()>0
                idxs.extend(ids[keep]);vals.extend(area.ravel()[keep]*H/dvol)
                crowns.extend(ca.ravel()[keep]*H/dvol);heights.extend(np.full(sum(keep),(lo+hi)/2))
                target+=t.crown_drag_coefficient*crownarea+t.stem_drag_coefficient*stemarea
            ids=np.array(idxs,int);values=np.array(vals);crown=np.array(crowns)
            np.add.at(self.coefficient.ravel(),ids,values)
            actual=float(values.sum()*dvol/H)
            self.entries.append(dict(tree=p,indices=ids,coefficient=values,crown_coefficient=crown,
                z_m=np.array(heights),target_cd_area_m2=target,represented_cd_area_m2=actual))
        self.coefficient_sha256=hashlib.sha256(self.coefficient.astype('<f8').tobytes()).hexdigest()
        self.area_relative_error=max((abs(e['represented_cd_area_m2']-e['target_cd_area_m2'])/
                max(e['target_cd_area_m2'],1e-12) for e in self.entries),default=0.)

    def loads(self, result, upstream_speed_m_s, reference_density_kg_m3):
        if result.grid.config!=self.grid.config or result.diagnostics.get('drag_coefficient_sha256')!=self.coefficient_sha256:
            raise ValueError('Flow state does not match this geometry/grid; solve the new drag field')
        if not np.isfinite(upstream_speed_m_s) or upstream_speed_m_s<0 or not np.isfinite(reference_density_kg_m3) or reference_density_kg_m3<=0:
            raise ValueError('Nonnegative speed and positive density required')
        H=self.grid.config.reference_height_m
        q=.5*reference_density_kg_m3*upstream_speed_m_s**2
        v=result.centered.reshape((-1,3));mag=np.linalg.norm(v,axis=1)
        scale=q*H**2*self.grid.cell_volume_h3
        out={}
        total=np.zeros(3)
        for e in self.entries:
            ids=e['indices'];u=v[ids]
            f=scale*e['coefficient'][:,None]*mag[ids,None]*u
            fc=scale*e['crown_coefficient'][:,None]*mag[ids,None]*u
            total+=f.sum(axis=0)
            out[e['tree'].identifier]=dict(z_m=e['z_m'],force_n=f,crown_force_n=fc,
                total_force_n=f.sum(axis=0),base_wind_moment_nm=(f*e['z_m'][:,None]).sum(axis=0),
                drag_weighted_speed_ratio=float(np.average(mag[ids],weights=e['coefficient'])))
        bulk=np.asarray(result.diagnostics['integrated_drag_h2'])*reference_density_kg_m3*upstream_speed_m_s**2*H**2
        return out,dict(sum_tree_force_n=total.tolist(),bulk_force_n=bulk.tolist(),
            force_partition_relative_error=float(np.linalg.norm(total-bulk)/max(np.linalg.norm(bulk),1e-12)),
            area_partition_relative_error=self.area_relative_error)
