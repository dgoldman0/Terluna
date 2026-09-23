"""Coupled small-displacement beams, compliant roots and optional neighbor links.

Reuses checkpoint beam/weight geometry. Contact is compression-only; clonal-root
links are optional finite-capacity rotational springs. Soil footprints share a
single area budget. Every strength/compliance parameter is hypothetical.
"""
from __future__ import annotations
from dataclasses import dataclass,asdict
import math
import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from scipy.sparse.linalg import spsolve
from biosphere.megaforest_mechanics import StaticTree,beam_matrices


@dataclass(frozen=True)
class Support:
    crown_mode: str = 'none'       # none, contact, bonded (proposed connection)
    crown_stiffness_n_m: float = 2e4
    crown_strength_n: float = 4e5
    contact_height_fraction: float = .70
    root_rotation_at_capacity_rad: float = .04
    root_link_nm_rad: float = 0.
    root_link_strength_nm: float = 1e8
    root_link_range_m: float = 120.
    soil_cell_m: float = 5.

    def validate(self):
        if self.crown_mode not in ('none','contact','bonded'): raise ValueError('Invalid support mode')
        for k,v in asdict(self).items():
            if k=='crown_mode':continue
            if not np.isfinite(v) or (v<0 if k=='root_link_nm_rad' else v<=0): raise ValueError('Invalid '+k)
        if not 0<self.contact_height_fraction<1: raise ValueError('Invalid contact level')


def share_soil(forest,cell_m=5.):
    """Allocate each horizontal soil cell once, equally among overlapping plates.

    Capacity is subsequently scaled by area fraction. Depth/load direction and
    soil constitutive interaction remain unresolved. Ownership is frozen during
    a damage sequence, so lost trees do not instantly donate soil to neighbors.
    """
    if not np.isfinite(cell_m) or cell_m<=0: raise ValueError('Invalid soil grid')
    if not forest: return np.array([]),dict(union_area_m2=0.,allocated_area_m2=0.)
    x0=min(f.x-f.roots.plate_radius_m for f in forest);x1=max(f.x+f.roots.plate_radius_m for f in forest)
    y0=min(f.y-f.roots.plate_radius_m for f in forest);y1=max(f.y+f.roots.plate_radius_m for f in forest)
    nx=max(1,int(np.ceil((x1-x0)/cell_m)));ny=max(1,int(np.ceil((y1-y0)/cell_m)))
    if nx*ny*len(forest)>100_000_000: raise ValueError('Soil grid too large')
    x=x0+(np.arange(nx)+.5)*cell_m;y=y0+(np.arange(ny)+.5)*cell_m
    masks=np.array([((x[:,None]-f.x)**2+(y[None,:]-f.y)**2)<f.roots.plate_radius_m**2 for f in forest])
    count=masks.sum(axis=0);weights=masks/np.maximum(count,1)
    raw=masks.sum(axis=(1,2));alloc=weights.sum(axis=(1,2))
    if np.any(raw==0): raise ValueError('Soil grid does not resolve a root plate')
    return alloc/raw,dict(union_area_m2=float(np.sum(count>0)*cell_m**2),allocated_area_m2=float(alloc.sum()*cell_m**2),
        cell_m=cell_m,overlapping_area_m2=float(np.sum(count>1)*cell_m**2),fractions=(alloc/raw).tolist())


def shape_row(z,height):
    """Cubic Hermite point-load / displacement row; base translation retained."""
    if not np.isfinite(height) or height < -1e-10 or height>z[-1]+1e-10: raise ValueError('Load height outside tree')
    i=min(len(z)-2,max(0,int(np.searchsorted(z,height,side='right')-1)))
    L=z[i+1]-z[i];s=np.clip((height-z[i])/L,0,1)
    row=np.zeros(2*len(z));row[2*i:2*i+4]=[1-3*s*s+2*s**3,L*(s-2*s*s+s**3),3*s*s-2*s**3,L*(-s*s+s**3)]
    return row


class StandMechanics:
    def __init__(self,forest,air,support=Support(),soil_fractions=None,broken_links=()):
        support.validate();self.forest=forest;self.support=support;self.broken=set(broken_links)
        self.soil_fractions,self.soil_ledger=share_soil(forest,support.soil_cell_m) if soil_fractions is None else (np.asarray(soil_fractions,float),{'ownership':'frozen_from_initial_patch'})
        if self.soil_fractions.shape!=(len(forest),) or np.any(self.soil_fractions<=0) or np.any(self.soil_fractions>1): raise ValueError('Invalid soil fractions')
        self.models={};self.offset={};self.matrices={};self.root_k={};self.cap={};blocks=[];offset=0;self.unstable=[]
        for i,f in enumerate(forest):
            if not f.alive:continue
            m=StaticTree(air,f.effective(),f.roots);self.models[i]=m
            t=m.tree;dm=(m.diameter[1:]+m.diameter[:-1])/2
            k,kg=beam_matrices(m.z,t.young_modulus_pa*math.pi*dm**4/64,(m.axial[1:]+m.axial[:-1])/2)
            cap=min(m.anchorage['tissue_capacity_nm'],m.anchorage['soil_capacity_nm']*self.soil_fractions[i])
            rk=cap/support.root_rotation_at_capacity_rad;kr=k[1:,1:].copy();kr[0,0]+=rk
            buck=float(eigh(kg[1:,1:],kr,eigvals_only=True,subset_by_index=[len(kr)-1,len(kr)-1])[0])
            if buck>=1:self.unstable.append(f.key)
            self.root_k[i]=rk;self.cap[i]=cap;self.matrices[i]=k-kg
            block=kr-kg[1:,1:];n=len(block);self.offset[i]=(offset,n);offset+=2*n
            blocks.extend([sparse.csr_matrix(block),sparse.csr_matrix(block)])
        self.K=sparse.block_diag(blocks,format='csc') if blocks else sparse.csc_matrix((0,0))
        self.links=[]
        active=list(self.models)
        for p,i in enumerate(active):
            fi=forest[i];ti=self.models[i].tree
            for j in active[p+1:]:
                fj=forest[j];tj=self.models[j].tree
                delta=np.array([fj.x-fi.x,fj.y-fi.y]);distance=np.linalg.norm(delta)
                if distance<=0: raise ValueError('Coincident tree centers')
                normal=delta/distance;h=support.contact_height_fraction*min(ti.height_m,tj.height_m)
                def radius(t):
                    s=(h/t.height_m-t.crown_base_fraction)/(1-t.crown_base_fraction)
                    return t.crown_radius_m*math.sqrt(max(0,4*s*(1-s))) if 0<s<1 else 0.
                radii=radius(ti)+radius(tj)
                if support.crown_mode!='none' and radii>0 and distance<radii+.15*min(ti.height_m,tj.height_m):
                    b=self._row(i,h,normal)-self._row(j,h,normal)
                    key=f'crown:{fi.key}:{fj.key}'
                    if key not in self.broken:
                        self.links.append(dict(key=key,kind='crown',i=i,j=j,b=b,k=support.crown_stiffness_n_m,
                            gap=max(0,float(distance-radii)) if support.crown_mode=='contact' else 0.,
                            unilateral=support.crown_mode=='contact',strength=support.crown_strength_n,h=h,normal=normal))
                if support.root_link_nm_rad>0 and distance<=support.root_link_range_m:
                    for axis in range(2):
                        b=np.zeros(self.K.shape[0]);oi,ni=self.offset[i];oj,nj=self.offset[j]
                        b[oi+axis*ni]=1;b[oj+axis*nj]=-1
                        key=f'root{axis}:{fi.key}:{fj.key}'
                        if key not in self.broken:
                            self.links.append(dict(key=key,kind='root',i=i,j=j,b=b,k=support.root_link_nm_rad,
                                gap=0.,unilateral=False,strength=support.root_link_strength_nm,h=0.,normal=np.eye(2)[axis]))
        self.B=sparse.csr_matrix(np.vstack([l['b'] for l in self.links])) if self.links else sparse.csr_matrix((0,self.K.shape[0]))
        self.lk=np.array([l['k'] for l in self.links]);self.gap=np.array([l['gap'] for l in self.links])
        self.unilateral=np.array([l['unilateral'] for l in self.links],bool)

    def _row(self,i,h,normal):
        offset,n=self.offset[i];b=np.zeros(self.K.shape[0]);v=shape_row(self.models[i].z,h)[1:]
        for axis in range(2):b[offset+axis*n:offset+(axis+1)*n]=normal[axis]*v
        return b

    def evaluate(self,drag,forces):
        if self.unstable: return dict(status='UNSTABLE_WITHOUT_PREENGAGED_SUPPORT',unstable_trees=self.unstable,trees=[],links=[])
        if forces.shape!=(len(drag.area),3) or not np.all(np.isfinite(forces)): raise ValueError('Invalid aerodynamic forces')
        if not self.models:return dict(status='EMPTY_STAND',trees=[],links=[],equilibrium_relative_residual=0.)
        load=np.zeros(self.K.shape[0]);points={i:[] for i in self.models};full_load={i:np.zeros((2,2*len(m.z))) for i,m in self.models.items()}
        for n,(i,h,force) in enumerate(zip(drag.owner,drag.heights,forces)):
            if i not in self.models: raise ValueError('Drag load applied to lost tree')
            row=shape_row(self.models[i].z,h);o,nfree=self.offset[i]
            for axis in range(2):
                full_load[i][axis]+=force[axis]*row
                load[o+axis*nfree:o+(axis+1)*nfree]+=force[axis]*row[1:]
            points[i].append((h,force[:2].copy(),int(drag.kind[n])))
        d=spsolve(self.K,load)
        def elong(q):
            r=self.B@q-self.gap
            return np.where(self.unilateral,np.maximum(r,0),r)
        def energy(q):
            e=elong(q);return .5*q@(self.K@q)-load@q+.5*np.sum(self.lk*e*e)
        iterations=0
        for iterations in range(40):
            e=elong(d);grad=self.K@d-load+self.B.T@(self.lk*e)
            res=float(np.linalg.norm(grad)/max(np.linalg.norm(load),1.))
            if res<2e-8:break
            active=(~self.unilateral)|(self.B@d>self.gap)
            ba=self.B[active];ka=self.lk[active]
            tangent=self.K+ba.T@sparse.diags(ka)@ba
            step=spsolve(tangent,-grad);alpha=1.;E=energy(d)
            while alpha>1e-7 and energy(d+alpha*step)>E+1e-4*alpha*(grad@step):alpha*=.5
            if alpha<=1e-7: raise RuntimeError('Contact energy line search failed')
            d+=alpha*step
        else:raise RuntimeError('Contact equilibrium did not converge')
        link_load=self.lk*elong(d);link_results=[];contact_branch={i:0. for i in self.models}
        for l,F in zip(self.links,link_load):
            link_results.append(dict(key=l['key'],kind=l['kind'],load=float(F),utilization=abs(float(F))/l['strength']))
            if l['kind']=='crown':
                for i,sign in ((l['i'],-1),(l['j'],1)):
                    force=sign*F*l['normal'];points[i].append((l['h'],force,2));contact_branch[i]+=abs(F)
                    full_load[i]+=force[:,None]*shape_row(self.models[i].z,l['h'])
            else:
                # External equal/opposite couples at the two root nodes.
                for i,sign in ((l['i'],-1),(l['j'],1)):
                    full_load[i][:,1]+=sign*F*l['normal']
        results=[];reactions=np.zeros(2)
        for i,m in self.models.items():
            o,n=self.offset[i];dofs=np.array([np.r_[0,d[o+a*n:o+(a+1)*n]] for a in range(2)])
            y=dofs[:,::2];slope=dofs[:,1::2]
            ym=(y[:,:-1]+y[:,1:])/2+m.dz*(slope[:,:-1]-slope[:,1:])/8
            moments=np.zeros_like(y)
            for h,force,_ in points[i]:moments+=force[:,None]*np.maximum(0,h-m.z)[None]
            for a in range(2):
                tail=np.r_[np.cumsum((m.weights*ym[a])[::-1])[::-1],0.]
                moments[a]+=tail-y[a]*m.axial
                reactions[a]+=(self.matrices[i]@dofs[a]-full_load[i][a])[0]
            bending=np.linalg.norm(moments,axis=0)/m.section_modulus;compression=m.axial/m.area;t=m.tree
            stem=float(np.max(np.maximum((bending+compression)/t.compression_strength_pa,np.maximum(0,bending-compression)/t.bending_strength_pa)))
            root=float(self.root_k[i]*np.linalg.norm(dofs[:,1])/self.cap[i])
            crownforce=np.sum([p[1] for p in points[i] if p[2]==1],axis=0) if any(p[2]==1 for p in points[i]) else np.zeros(2)
            branchwind=(np.linalg.norm(crownforce)+contact_branch[i])*m.span/(2*t.branch_count)
            crown=math.hypot(branchwind,m.branch_gravity_moment)/(m.branch_section_modulus*t.branch_strength_pa)
            domain=max(float(np.max(np.linalg.norm(y,axis=0)))/(.1*t.height_m),float(np.max(np.linalg.norm(slope,axis=0)))/.2)
            own=forces[drag.owner==i];vertical=float(abs(own[:,2].sum())/max(m.weights.sum(),1.))
            results.append(dict(key=self.forest[i].key,index=i,x_m=self.forest[i].x,y_m=self.forest[i].y,height_m=t.height_m,
                stem_utilization=stem,root_utilization=root,crown_utilization=crown,domain_utilization=domain,
                tip_displacement_m=float(np.linalg.norm(y[:,-1])),root_rotation_rad=float(np.linalg.norm(dofs[:,1])),
                base_stem_moment_nm=float(np.linalg.norm(moments[:,0])),anchorage_capacity_nm=self.cap[i],
                crown_contact_force_n=contact_branch[i],aerodynamic_force_x_n=float(own[:,0].sum()),
                aerodynamic_force_y_n=float(own[:,1].sum()),vertical_drag_to_weight=vertical,
                soil_share=float(self.soil_fractions[i]),crown_remaining=self.forest[i].crown_remaining,
                domain_valid=bool(domain<=1 and vertical<=.05)))
        total=forces[:,:2].sum(axis=0)
        force_res=float(np.linalg.norm(reactions+total)/max(np.linalg.norm(total),1.))
        return dict(status='COUPLED_SMALL_DISPLACEMENT_STAND',trees=results,links=link_results,
            equilibrium_relative_residual=res,horizontal_force_balance_relative_residual=force_res,
            contact_iterations=iterations,active_crown_contacts=sum(l['kind']=='crown' and abs(F)>1e-6 for l,F in zip(self.links,link_load)),
            assumptions=['hypothetical_root_and_link_properties','upright_aerodynamic_geometry',
                'vertical_aerodynamic_load_excluded_from_mechanics_5pct_guard','no_torsion_or_collisions','shared_soil_area_proxy'])
