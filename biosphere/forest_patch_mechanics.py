"""Sparse two-plane beam stand, compliant foundations and optional interactions.

Reuses the established beam, self-weight and capacity functions. Crown contacts
are unilateral horizontal springs. Optional grafts are rotational spring/torque
proxies between foundations, with explicit finite breaking moments. These are
hypothetical mechanical connections, never assumed for ordinary neighboring trees.
Shared soil area is partitioned before foundation capacities are calculated.
Large deflection, collision dynamics, branch networks and plasticity remain open.
"""
from __future__ import annotations
from dataclasses import dataclass,asdict
import math
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.linalg import cholesky
from biosphere.megaforest_mechanics import StaticTree,beam_matrices
from biosphere.forest_patch import crown_radius_at


@dataclass(frozen=True)
class InteractionConfig:
    crown_contacts: bool = False
    root_grafts: bool = False
    foundation_rotation_at_capacity_rad: float = .025
    contact_stiffness_n_m: float = 100000.
    contact_break_force_n: float = 2500000.
    graft_stiffness_fraction: float = .2
    graft_break_moment_fraction: float = .3
    graft_range_m: float = 130.
    soil_cell_m: float = 5.
    max_contact_iterations: int = 60

    def validate(self):
        vals=[v for v in asdict(self).values() if not isinstance(v,bool)]
        if not all(np.isfinite(v) and v>0 for v in vals):
            raise ValueError('Positive finite interaction parameters required')
        if not isinstance(self.max_contact_iterations,int):raise ValueError('Integer iteration limit required')


def point_shape(z, height):
    """Cubic Hermite interpolation: displacement and slope DOFs."""
    if not np.isfinite(height) or height<0 or height>z[-1]+1e-8:
        raise ValueError('Attachment outside tree')
    i=min(np.searchsorted(z,height,side='right')-1,len(z)-2);i=max(0,i)
    l=z[i+1]-z[i];s=np.clip((height-z[i])/l,0,1)
    out=np.zeros(2*len(z))
    out[2*i:2*i+4]=[1-3*s*s+2*s**3,l*(s-2*s*s+s**3),3*s*s-2*s**3,l*(-s*s+s**3)]
    return out


def partition_soil(trees, cell_m):
    if not np.isfinite(cell_m) or cell_m<=0:raise ValueError('Positive finite soil cell size required')
    active=list(trees)  # Disturbed/failed root zones stay reserved during this event.
    if not active:return {},dict(unique_area_m2=0.,allocated_area_m2=0.,partition_error_m2=0.)
    xmin=min(p.x_m-p.roots.plate_radius_m for p in active);xmax=max(p.x_m+p.roots.plate_radius_m for p in active)
    ymin=min(p.y_m-p.roots.plate_radius_m for p in active);ymax=max(p.y_m+p.roots.plate_radius_m for p in active)
    x=np.arange(xmin+cell_m/2,xmax,cell_m);y=np.arange(ymin+cell_m/2,ymax,cell_m)
    xx,yy=np.meshgrid(x,y,indexing='ij');counts=np.zeros(xx.shape)
    masks=[]
    for p in active:
        m=(xx-p.x_m)**2+(yy-p.y_m)**2<=p.roots.plate_radius_m**2
        if not m.any():raise ValueError('Soil grid misses a root plate; refine soil_cell_m')
        masks.append(m);counts+=m
    allocation={};area=0.
    for p,m in zip(active,masks):
        allocated=float(np.sum(1/counts[m])*cell_m**2)
        original=float(m.sum()*cell_m**2)
        allocation[p.identifier]=dict(capacity_fraction=allocated/original,
            represented_plate_area_m2=original,allocated_area_m2=allocated)
        area+=allocated
    unique=float(np.sum(counts>0)*cell_m**2)
    return allocation,dict(unique_area_m2=unique,allocated_area_m2=area,
        partition_error_m2=area-unique,cell_m=cell_m,
        assumption='area-fraction capacity partition; soil stress field and anisotropic moment arms unresolved')


class PatchStructure:
    def __init__(self, air, trees, interactions=InteractionConfig(), broken_links=()):
        interactions.validate();self.config=interactions;all_trees=tuple(trees)
        if len({p.identifier for p in all_trees})!=len(all_trees):
            raise ValueError('Unique tree identifiers required')
        self.trees=[p for p in all_trees if p.alive]
        self.soil,self.soil_diagnostic=partition_soil(all_trees,interactions.soil_cell_m)
        self.models={};self.info={};blocks=[];offset=0
        self.unstable=[]
        for p in self.trees:
            t=p.effective_tree();m=StaticTree(air,t,p.roots);self.models[p.identifier]=m
            dmid=(m.diameter[:-1]+m.diameter[1:])/2
            k,kg=beam_matrices(m.z,t.young_modulus_pa*np.pi*dmid**4/64,
                              (m.axial[:-1]+m.axial[1:])/2)
            soilcap=m.anchorage['soil_capacity_nm']*self.soil[p.identifier]['capacity_fraction']
            rootk=soilcap/interactions.foundation_rotation_at_capacity_rad
            matrix=k-kg;matrix[1,1]+=rootk
            # Fix translation only. Scale rotations by height for conditioning.
            scale=np.ones(len(matrix)-1);scale[::2]=1/t.height_m
            reduced=matrix[1:,1:]*scale[:,None]*scale[None,:]
            try:cholesky(reduced,lower=True)
            except np.linalg.LinAlgError:self.unstable.append(p.identifier)
            n=len(scale)
            self.info[p.identifier]=dict(offset=offset,n=n,scale=scale,rootk=rootk,soilcap=soilcap,
                tissuecap=m.anchorage['tissue_capacity_nm'],tree=p)
            blocks.extend((sparse.csr_matrix(reduced),sparse.csr_matrix(reduced)));offset+=2*n
        self.size=offset
        self.base=sparse.block_diag(blocks,format='csr') if blocks else sparse.csr_matrix((0,0))
        self.links=[];broken=set(broken_links)
        if not self.unstable:
            for i,p in enumerate(self.trees):
                for q in self.trees[i+1:]:
                    dx=q.x_m-p.x_m;dy=q.y_m-p.y_m;distance=math.hypot(dx,dy)
                    if distance==0:raise ValueError('Trees require distinct positions')
                    normal=np.array([dx,dy])/distance
                    if interactions.crown_contacts and p.remaining_crown>0 and q.remaining_crown>0:
                        t=p.effective_tree();u=q.effective_tree()
                        lo=max(t.height_m*t.crown_base_fraction,u.height_m*u.crown_base_fraction)
                        hi=min(t.height_m,u.height_m)
                        if hi>lo:
                            h=(lo+hi)/2;r=crown_radius_at(t,h)+crown_radius_at(u,h)
                            # Nearby crowns can meet after displacement; distant
                            # pairs are outside the small-displacement horizon.
                            if distance<r+.2*min(t.height_m,u.height_m):
                                name=f'contact:{p.identifier}:{q.identifier}'
                                if name not in broken:
                                    a=self._point_vector(p.identifier,h,normal)-self._point_vector(q.identifier,h,normal)
                                    self.links.append(dict(id=name,kind='contact',a=a,gap=max(0,distance-r),
                                        stiffness=interactions.contact_stiffness_n_m,
                                        capacity=interactions.contact_break_force_n,
                                        first=p.identifier,second=q.identifier,height_m=h,normal=normal))
                    if interactions.root_grafts and distance<=interactions.graft_range_m:
                        for axis in (0,1):
                            name=f'graft{axis}:{p.identifier}:{q.identifier}'
                            if name in broken:continue
                            a=np.zeros(self.size)
                            for sign,tree in ((1,p),(-1,q)):
                                inf=self.info[tree.identifier]
                                a[inf['offset']+axis*inf['n']]=sign*inf['scale'][0]
                            ip,iq=self.info[p.identifier],self.info[q.identifier]
                            self.links.append(dict(id=name,kind='graft',a=a,gap=0.,
                                stiffness=interactions.graft_stiffness_fraction*min(ip['rootk'],iq['rootk']),
                                capacity=interactions.graft_break_moment_fraction*min(ip['tissuecap'],iq['tissuecap']),
                                first=p.identifier,second=q.identifier,height_m=0.,axis=axis,normal=normal))

    def _point_vector(self,identifier,h,direction):
        m=self.models[identifier];i=self.info[identifier];v=np.zeros(self.size)
        shape=point_shape(m.z,h)[1:]*i['scale']
        for axis in (0,1):
            o=i['offset']+axis*i['n'];v[o:o+i['n']]=shape*direction[axis]
        return v

    def solve(self,loads):
        if self.unstable:
            return dict(status='FOUNDATION_OR_SELF_WEIGHT_INSTABILITY',unstable=self.unstable,trees=[],links=[],domain_limited=True)
        if not self.trees:return dict(status='EMPTY_STAND',trees=[],links=[],domain_limited=False)
        rhs=np.zeros(self.size)
        for p in self.trees:
            for z,f in zip(loads[p.identifier]['z_m'],loads[p.identifier]['force_n']):
                rhs+=self._point_vector(p.identifier,z,f[:2])
        def stiffness(active):
            matrix=self.base.copy();f=rhs.copy()
            for n in active:
                e=self.links[n];a=sparse.csr_matrix(e['a'].reshape(1,-1))
                matrix=matrix+e['stiffness']*(a.T@a)
                f+=e['stiffness']*e['gap']*e['a']
            return matrix.tocsc(),f
        active={n for n,l in enumerate(self.links) if l['kind']=='graft'}
        last=None
        for iteration in range(self.config.max_contact_iterations):
            matrix,f=stiffness(active);displacements=spsolve(matrix,f)
            if not np.all(np.isfinite(displacements)):raise RuntimeError('Nonfinite structure solution')
            new={n for n,l in enumerate(self.links) if l['kind']=='graft' or l['a']@displacements>l['gap']+1e-9}
            if new==active:break
            if new==last:raise RuntimeError('Contact active set oscillated; reduce load/increase resolution')
            last=active;active=new
        else:raise RuntimeError('Contact solve failed to converge')
        residual=float(np.linalg.norm(matrix@displacements-f)/max(np.linalg.norm(f),1e-12))
        contactforces={p.identifier:[] for p in self.trees};grafttorques={p.identifier:np.zeros(2) for p in self.trees}
        linkrecords=[];net_contact=np.zeros(2);net_graft=np.zeros(2)
        for n,l in enumerate(self.links):
            force=l['stiffness']*(l['a']@displacements-l['gap']) if n in active else 0.
            if l['kind']=='contact':
                force=max(0,force)
                for sign,key in ((-1,l['first']),(1,l['second'])):
                    value=sign*force*l['normal'];contactforces[key].append((l['height_m'],value));net_contact+=value
            else:
                for sign,key in ((-1,l['first']),(1,l['second'])):
                    v=np.zeros(2);v[l['axis']]=sign*force;grafttorques[key]+=v;net_graft+=v
            linkrecords.append(dict(id=l['id'],kind=l['kind'],first=l['first'],second=l['second'],
                active=n in active,force_or_torque=force,utilization=abs(force)/l['capacity'],
                units='N' if l['kind']=='contact' else 'Nm'))
        records=[]
        for p in self.trees:
            key=p.identifier;m=self.models[key];inf=self.info[key];t=m.tree
            dofs=[]
            for a in (0,1):
                o=inf['offset']+a*inf['n'];dofs.append(np.r_[0.,displacements[o:o+inf['n']]*inf['scale']])
            dofs=np.stack(dofs,axis=-1);y=dofs[::2];slope=dofs[1::2]
            ym=(y[:-1]+y[1:])/2+m.dz[:,None]*(slope[:-1]-slope[1:])/8
            moment=np.zeros((len(m.z),2));load=loads[key]
            for z,f in zip(load['z_m'],load['force_n']):moment+=np.maximum(0,z-m.z)[:,None]*f[None,:2]
            for z,f in contactforces[key]:moment+=np.maximum(0,z-m.z)[:,None]*f
            for j,z in enumerate(m.z):
                mask=m.mid>z
                moment[j]+=(m.weights[mask,None]*(ym[mask]-y[j])).sum(axis=0)
            bend=np.linalg.norm(moment,axis=1)/m.section_modulus
            compression=m.axial/m.area
            stem=max(np.max((bend+compression)/t.compression_strength_pa),np.max(np.maximum(0,bend-compression)/t.bending_strength_pa))
            soilmoment=inf['rootk']*slope[0]
            root=max(np.linalg.norm(moment[0])/inf['tissuecap'],np.linalg.norm(soilmoment)/inf['soilcap'])
            crown_force=load['crown_force_n'].sum(axis=0)[:2]
            for z,f in contactforces[key]:crown_force+=f
            branch_wind=np.linalg.norm(crown_force)*m.span/(2*t.branch_count)
            branch=math.hypot(branch_wind,m.branch_gravity_moment)/(m.branch_section_modulus*t.branch_strength_pa)
            domain=max(float(np.max(np.linalg.norm(y,axis=1)))/(.1*t.height_m),
                float(np.max(np.linalg.norm(slope,axis=1)))/.2)
            records.append(dict(id=key,row=p.row,column=p.column,x_m=p.x_m,y_m=p.y_m,height_m=t.height_m,
                stem_utilization=float(stem),root_utilization=float(root),crown_utilization=float(branch),
                domain_utilization=domain,tip_displacement_m=float(np.linalg.norm(y[-1])),
                base_rotation_rad=float(np.linalg.norm(slope[0])),
                root_soil_moment_nm=float(np.linalg.norm(soilmoment)),base_moment_nm=float(np.linalg.norm(moment[0])),
                total_horizontal_drag_n=float(np.linalg.norm(load['total_force_n'][:2])),
                vertical_drag_n=float(load['total_force_n'][2]),
                drag_weighted_speed_ratio=load['drag_weighted_speed_ratio'],
                crown_fraction=p.remaining_crown,soil_capacity_fraction=self.soil[key]['capacity_fraction']))
        return dict(status='CONDITIONAL_LINEAR_INTERACTING_STAND',trees=records,links=linkrecords,
            contact_iterations=iteration+1,linear_equilibrium_relative_residual=residual,
            net_internal_contact_force_n=net_contact.tolist(),net_internal_graft_moment_nm=net_graft.tolist(),
            soil_partition=self.soil_diagnostic,domain_limited=any(r['domain_utilization']>1 for r in records),
            assumptions=['small-displacement two-plane beam with linearized self-weight',
                'quasi-static unilateral contacts; impact damping omitted',
                'finite-stiffness hypothetical root grafts',
                'vertical aerodynamic loads recorded but not applied to axial stability',
                'contact forces aggregated for branch onset; branch network unresolved'])
