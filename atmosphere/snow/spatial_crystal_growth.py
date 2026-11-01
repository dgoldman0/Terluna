"""Research prototype: conservative event growth on an axisymmetric voxel grid.

Published local attachment and optional 2015 terrace-distance rule. This is an
own implementation, not the original cellular automaton. Fixed physical cells,
finite outer boundary, no flow solver, stochastic nucleation or breakup.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.ndimage import label
from scipy.sparse import coo_matrix,diags
from scipy.sparse.linalg import spsolve
from atmosphere.snow.spatial_diffusion import Attachment,surface_flux,Domain
from atmosphere.snow.facet_growth import ICE_DENSITY


def exterior_air_masks(solid, reflect_top=False):
    """Face-connected exterior gas and sealed cells of the filled-cell grid.

    A sealed component has no ambient supply in the quasistatic problem.
    Excluding it from that solve does not turn its air or partial ice into ice.
    The moving solver separately bounds its finite trapped-vapour inventory.
    """
    gas = ~np.asarray(solid, bool)
    components, _ = label(gas)
    outer = components[-1, :]
    if not reflect_top:
        outer = np.r_[outer, components[:, -1]]
    connected = np.unique(outer[outer > 0])
    exterior = gas & np.isin(components, connected)
    return exterior, gas & ~exterior


class MaskDomain(Domain):
    def __init__(self,rf,zf,solid,reflect_top=False):
        self.rf=np.asarray(rf,float);self.zf=np.asarray(zf,float)
        self.r=(self.rf[:-1]+self.rf[1:])/2;self.z=(self.zf[:-1]+self.zf[1:])/2
        self.solid=np.asarray(solid,bool);self.shape='voxel';self._capacity=None
        self.aspect=float(max(self.z[np.any(solid,axis=0)])/max(self.r[np.any(solid,axis=1)]))
        self.resolution=1;self.extent=max(self.rf[-1],self.zf[-1]);self.mesh='voxel'
        self.exterior_air,self.closed_air=exterior_air_masks(self.solid,reflect_top)
        index=np.full(solid.shape,-1,int);index[self.exterior_air]=np.arange(np.sum(self.exterior_air));self.index=index;self.cells=int(np.sum(self.exterior_air))
        n=self.cells;diag=np.zeros(n);outer=np.zeros(n);rows=[];cols=[];vals=[];ice=[]
        # Interior face conductances, including the exact cylindrical areas.
        for axis in (0,1):
            if axis==0:
                left=index[:-1,:];right=index[1:,:]
                rr,zz=np.meshgrid(self.rf[1:-1],self.z,indexing='ij')
                area=2*math.pi*rr*np.diff(self.zf)[None,:]
                distance=np.diff(self.r)[:,None]*np.ones_like(area)
                half_left=np.broadcast_to((self.rf[1:-1]-self.r[:-1])[:,None],area.shape)
                half_right=np.broadcast_to((self.r[1:]-self.rf[1:-1])[:,None],area.shape)
            else:
                left=index[:,:-1];right=index[:,1:]
                rr,zz=np.meshgrid(self.r,self.zf[1:-1],indexing='ij')
                area=math.pi*np.diff(self.rf**2)[:,None]*np.ones_like(rr)
                distance=np.diff(self.z)[None,:]*np.ones_like(area)
                half_left=np.broadcast_to((self.zf[1:-1]-self.z[:-1])[None,:],area.shape)
                half_right=np.broadcast_to((self.z[1:]-self.zf[1:-1])[None,:],area.shape)
            both=(left>=0)&(right>=0);p=left[both];q=right[both];g=area[both]/distance[both]
            np.add.at(diag,p,g);np.add.at(diag,q,g)
            rows.extend([p,q]);cols.extend([q,p]);vals.extend([-g,-g])
            for ids,other,half in [(left,right,half_left),(right,left,half_right)]:
                select=(ids>=0)&(other<0)
                if np.any(select):ice.append(np.column_stack([ids[select],area[select],half[select],np.full(np.sum(select),axis),rr[select],zz[select]]))
        for axis in (0,1):
            if axis==1 and reflect_top:continue
            if axis==0:
                p=index[-1,:];area=2*math.pi*self.rf[-1]*np.diff(self.zf);half=self.rf[-1]-self.r[-1]
            else:
                p=index[:,-1];area=math.pi*np.diff(self.rf**2);half=self.zf[-1]-self.z[-1]
            if np.any(p<0):raise ValueError('Ice touched finite outer boundary')
            g=area/half;np.add.at(diag,p,g);np.add.at(outer,p,g)
        rows.append(np.arange(n));cols.append(np.arange(n));vals.append(diag)
        self.matrix=coo_matrix((np.concatenate(vals),(np.concatenate(rows),np.concatenate(cols))),shape=(n,n)).tocsr()
        self.outer=outer;self.ice=np.vstack(ice);self.normal_projection=np.ones(len(self.ice))

    def terrace_enhancement(self,n_surface_diffusion):
        """2015 rule: probability of reaching accessible pixel kink(s).

        n_SD = physical surface migration length / molecular step height,
        measured in grid pixels as in the source. This deliberately retains
        the source's resolution dependence; it is not a continuum limit.
        """
        if n_surface_diffusion<0:raise ValueError('Nonnegative surface diffusion scale required')
        p=self.ice[:,0].astype(int);orient=self.ice[:,3].astype(int)
        # Source neighbor rules assign rough uptake to every boundary pixel
        # with two or more ice neighbors, including a one-cell-wide gap.
        kink=np.bincount(p,minlength=self.cells)>=2
        enhance=np.where(kink[p],1.,0.)
        if n_surface_diffusion==0:return enhance
        ir,iz=np.where(self.exterior_air);cellr=self.r[ir];cellz=self.z[iz]
        for direction in (0,1):
            chosen=np.where(orient==direction)[0]
            coordinate=self.ice[:,4+direction]
            normal_sign=np.sign((cellr[p] if direction==0 else cellz[p])-coordinate)
            for position in np.unique(coordinate[chosen]):
                for sign in (-1,1):
                    group=chosen[(coordinate[chosen]==position)&(normal_sign[chosen]==sign)]
                    if not len(group):continue
                    tangential=self.ice[group,5-direction]
                    order=np.argsort(tangential);group=group[order];tangential=tangential[order]
                    # A gap means a disconnected terrace; no hopping through air.
                    cuts=np.r_[0,np.where(np.diff(tangential)>1.0000001)[0]+1,len(group)]
                    for lo,hi in zip(cuts[:-1],cuts[1:]):
                        ids=group[lo:hi];x=tangential[lo:hi];kinks=x[kink[p[ids]]]
                        if not len(kinks):continue
                        left=np.full(len(x),0.);right=left.copy()
                        for i,xx in enumerate(x):
                            ll=kinks[kinks<=xx];rr=kinks[kinks>=xx]
                            if len(ll):left[i]=math.exp(-(xx-ll[-1])/n_surface_diffusion)
                            if len(rr):right[i]=math.exp(-(rr[0]-xx)/n_surface_diffusion)
                        enhance[ids]=1-(1-left)*(1-right)
        return enhance

    def curvature_offsets(self,delta_over_cell):
        """Source-style outer-terrace curvature surrogate, kappa=2/width.

        Applies only to the outermost prism/basal terrace. This reproduces the
        stated plate/hollow-tip proxy, not a general curvature estimator.
        """
        offsets=np.zeros(len(self.ice))
        for axis in (0,1):
            chosen=np.where(self.ice[:,3]==axis)[0]
            if not len(chosen):continue
            normal=self.ice[chosen,4+axis];ids=chosen[normal==np.max(normal)]
            tangent=self.ice[ids,5-axis];order=np.argsort(tangent);ids=ids[order];tangent=tangent[order]
            cuts=np.r_[0,np.where(np.diff(tangent)>1.0000001)[0]+1,len(ids)]
            for lo,hi in zip(cuts[:-1],cuts[1:]):
                group=ids[lo:hi];offsets[group]=2*delta_over_cell/len(group)
        return offsets

    def coupled(self,si,lam,laws,heat_parameter=0.,guess=None,tolerance=1e-9):
        """Joint vapour/latent-heat Newton solve, one factorization per step."""
        if not math.isfinite(si) or si<0 or not math.isfinite(lam) or lam<=0 or not math.isfinite(heat_parameter) or heat_parameter<0:
            raise ValueError('Finite nonnegative forcing and heat, positive kinetic ratio required')
        if len(laws)!=2 or not math.isfinite(tolerance) or tolerance<=0:
            raise ValueError('Two attachment laws and positive tolerance required')
        p=self.ice[:,0].astype(int);area=self.ice[:,1];distance=self.ice[:,2]
        masks=[self.ice[:,3]==i for i in (0,1)]
        C=self.dirichlet()['half_conductance'] if heat_parameter else 1.
        h=0.;u=np.zeros(self.cells)
        if guess is not None:
            full,h0=guess;u=full[self.exterior_air].copy();h=float(h0)
            if si==0:u[:]=0.;h=0.
            elif np.any(u<0) or np.max(u)+h>=si:
                factor=.99*si/max(np.max(u)+h,si)
                u=np.maximum(u,0.)*factor;h=max(h,0.)*factor
        for it in range(50):
            q=np.zeros(len(p));dq=q.copy();ss=q.copy()
            for mask,law in zip(masks,laws):
                q[mask],dq[mask],ss[mask]=surface_flux(si-h-u[p[mask]],distance[mask],lam,law)
            source=np.bincount(p,weights=area*q,minlength=self.cells);Q=float(np.sum(source))
            residual=self.matrix@u-source;rh=h-heat_parameter*Q/C
            scale=max(Q,1e-30);error=max(np.max(np.abs(residual))/scale,abs(rh)/max(si,1e-30))
            if error<tolerance:
                full=np.zeros(self.solid.shape);full[self.exterior_air]=u
                mass_residual=abs(Q-np.dot(self.outer,u))/scale
                if mass_residual>2e-7:raise ArithmeticError('Water ledger failed')
                return dict(face_flux=q,surface=ss,half_flux=Q,heat_shift=h,guess=(full,h),iterations=it,
                    relative_mass_residual=float(mass_residual),relative_nonlinear_residual=float(error),capacity=C)
            g=np.bincount(p,weights=area*dq,minlength=self.cells)
            J=self.matrix+diags(g)
            two=spsolve(J,np.column_stack([-residual,g]))
            v,w=two[:,0],two[:,1]
            k=heat_parameter/C
            dh=(-rh-k*np.dot(g,v))/(1+k*np.sum(g)-k*np.dot(g,w))
            du=v-w*dh
            fraction=1.
            if dh<0:fraction=min(fraction,.99*(-h/dh))
            if np.any(du<0):fraction=min(fraction,.99*float(np.min(-u[du<0]/du[du<0])))
            together=du+dh;positive=together>0
            if np.any(positive):fraction=min(fraction,.99*float(np.min((si-h-u[positive])/together[positive])))
            u+=fraction*du;h+=fraction*dh
            # Tiny exact-equilibrium roundoff must not freeze the line search.
            u=np.maximum(u,0.);h=max(h,0.)
        raise ArithmeticError('Coupled diffusion/heat failed to converge')


class Combined:
    def __init__(self,law,enhancement,curvature=0.,kinetic_weight=1.):
        self.law=law;self.enhancement=np.asarray(enhancement)
        self.curvature=np.asarray(curvature);self.kinetic_weight=np.asarray(kinetic_weight)
    def flux(self,s):
        driving=np.maximum(0.,s-self.curvature)
        p,dp=self.law.flux(driving);a=self.enhancement
        return self.kinetic_weight*(a*driving+(1-a)*p),self.kinetic_weight*(a+(1-a)*dp)*(s>self.curvature)


def grid(core,extent):
    values=list(np.arange(core+1,dtype=float));step=1.
    while values[-1]<extent:
        values.append(min(extent,values[-1]+step));step*=1.18
    return np.array(values)


class GrowingCrystal:
    def __init__(self,dx,a,c,radial_core,vertical_core,extent=40.,cavity_radius=0.,cavity_depth=0.,adaptive=True):
        self.dx=dx;self.reflect_top=False;self.adaptive=adaptive;self.max_core=(radial_core,vertical_core)
        if adaptive:
            radial_core=min(radial_core,int(math.ceil(a/dx))+4)
            vertical_core=min(vertical_core,int(math.ceil(c/dx))+4)
        self.far=extent*max(a,c)/dx
        self.rf=grid(radial_core,self.far);self.zf=grid(vertical_core,self.far)
        rr,zz=np.meshgrid((self.rf[:-1]+self.rf[1:])/2,(self.zf[:-1]+self.zf[1:])/2,indexing='ij')
        solid=(rr*dx<a)&(zz*dx<c)
        if cavity_radius:solid&=~((rr*dx<cavity_radius)&(zz*dx>c-cavity_depth))
        self.solid=solid;self.fraction=solid.astype(float)
        self.volume=math.pi*np.diff(self.rf**2)[:,None]*np.diff(self.zf)[None,:]
        self.radial_core=radial_core;self.vertical_core=vertical_core
        self.time=0.;self.deposited=0.;self.initial_mass=self.mass();self.steps=0;self.last=None

    def mass(self):return float(2*ICE_DENSITY*self.dx**3*np.sum(self.volume*self.fraction))

    def expand(self):
        if not self.adaptive:return
        i,j=np.where(self.fraction>0)
        oldcores=(self.radial_core,self.vertical_core)
        cores=[oldcores[k] if max(ids)<oldcores[k]-3 else min(self.max_core[k],max(oldcores[k]+4,int(math.ceil(oldcores[k]*1.2)))) for k,ids in enumerate((i,j))]
        if tuple(cores)==oldcores:return
        rf=grid(cores[0],self.far);zf=grid(cores[1],self.far)
        fraction=np.zeros((len(rf)-1,len(zf)-1))
        fraction[:oldcores[0],:oldcores[1]]=self.fraction[:oldcores[0],:oldcores[1]]
        if np.any(self.fraction[oldcores[0]:,:]) or np.any(self.fraction[:,oldcores[1]:]):raise ArithmeticError('Regrid would lose partial ice')
        if self.last is not None:
            from scipy.interpolate import RegularGridInterpolator
            oldr=(self.rf[:-1]+self.rf[1:])/2;oldz=(self.zf[:-1]+self.zf[1:])/2
            interp=RegularGridInterpolator((oldr,oldz),self.last[0],bounds_error=False,fill_value=0.)
            rr,zz=np.meshgrid((rf[:-1]+rf[1:])/2,(zf[:-1]+zf[1:])/2,indexing='ij')
            self.last=(interp(np.column_stack([rr.ravel(),zz.ravel()])).reshape(rr.shape),self.last[1])
        self.rf=rf;self.zf=zf;self.fraction=fraction;self.solid=fraction>=1-1e-12
        self.volume=math.pi*np.diff(rf**2)[:,None]*np.diff(zf)[None,:]
        self.radial_core,self.vertical_core=cores

    def measured_seed(self,a,c,mass,tip_rim):
        """Axisymmetric hollow family constrained by initial mass and envelope.

        Inner radius squared = tip_radius^2 * (z/c)^power. The tip rim is a
        supplied unmeasured parameter; power follows from volume. Exact cell
        intersections conserve the supplied initial mass at every resolution.
        """
        if not 0<tip_rim<a:raise ValueError('Tip rim must lie within radius')
        f=mass/(2*math.pi*a*a*c*ICE_DENSITY);rt=a-tip_rim
        if not 0<f<1:raise ValueError('Mass must lie below the dense envelope')
        power=(rt/a)**2/(1-f)-1
        if power<=0:raise ValueError('Mass incompatible with this prescribed hollow profile')
        rlo=self.rf[:-1]*self.dx;rhi=self.rf[1:]*self.dx
        zlo=self.zf[:-1]*self.dx;zhi=self.zf[1:]*self.dx
        radial=np.maximum(0,np.minimum(rhi,a)**2-np.minimum(rlo,a)**2)
        axial=np.maximum(0,np.minimum(zhi,c)-np.minimum(zlo,c))
        actual=math.pi*radial[:,None]*axial[None,:]
        for i,(lo,hi) in enumerate(zip(rlo,rhi)):
            if lo>=rt:continue
            lower=c*(lo/rt)**(2/power)
            upper=c*(min(hi,rt)/rt)**(2/power)
            left=np.maximum(zlo,lower);right=np.minimum(np.minimum(zhi,c),upper)
            good=right>left
            integral=np.zeros_like(zlo)
            integral[good]=rt*rt/c**power*(right[good]**(power+1)-left[good]**(power+1))/(power+1)-lo*lo*(right[good]-left[good])
            full=np.maximum(0,np.minimum(zhi,c)-np.maximum(zlo,upper))*(min(hi,rt)**2-lo**2)
            # Beyond rt the cavity never fills the whole radial cell.
            actual[i]-=math.pi*(integral+full)
        fraction=actual/(self.volume*self.dx**3)
        if np.min(fraction)<-1e-10 or np.max(fraction)>1+1e-10:raise ArithmeticError('Seed cell overlap failed')
        self.fraction=np.clip(fraction,0.,1.);self.fraction[self.fraction>1-1e-12]=1.
        self.solid=self.fraction>=1.;self.initial_mass=self.mass();self.last=None
        if abs(self.initial_mass/mass-1)>1e-9:raise ArithmeticError('Initial measured mass not conserved')
        return dict(profile='inner_radius_squared_power_law',power=power,tip_rim_m=tip_rim,mass_fraction_of_cylinder=f)


    def state(self):
        i,j=np.where(self.solid)
        ii,jj=np.where(self.fraction>0)
        return dict(a_occupied_m=float(self.rf[max(ii)+1]*self.dx),c_occupied_m=float(self.zf[max(jj)+1]*self.dx),time_s=self.time,events=self.steps,a_m=float(self.rf[max(i)+1]*self.dx),c_m=float(self.zf[max(j)+1]*self.dx),
            mass_kg=self.mass(),mass_ledger_relative_residual=abs(self.mass()-self.initial_mass-self.deposited)/max(self.mass(),1e-30))

    def advance(self,air,T,si,laws,until, *, surface_diffusion_length=0.,molecular_step=0.3e-9,kinks=True,heating=True,ventilation=1.,tolerance=1e-9,capillary_length=0.):
        if not math.isfinite(until) or until<=self.time:
            raise ValueError('Finite target time after current time required')
        if any(not math.isfinite(x) or x<=0 for x in (T,ventilation,molecular_step)) or any(not math.isfinite(x) or x<0 for x in (surface_diffusion_length,capillary_length)):
            raise ValueError('Physical growth scales are out of range')
        self.expand()
        d=MaskDomain(self.rf,self.zf,self.solid,reflect_top=self.reflect_top)
        enhance=d.terrace_enhancement(surface_diffusion_length/molecular_step) if kinks else np.zeros(len(d.ice))
        curvature=d.curvature_offsets(capillary_length/self.dx)
        # Libbrecht 2015 eqs. (2)/(6): g=1/sqrt(N) prevents a 45-degree
        # boundary from receiving sqrt(2) times the isotropic kinetic uptake.
        # N>2 remains the source's approximate treatment of complex corners.
        neighbors=np.bincount(d.ice[:,0].astype(int),minlength=d.cells)
        weight=1/np.sqrt(neighbors[d.ice[:,0].astype(int)])
        combined=[Combined(law,enhance[d.ice[:,3]==k],curvature[d.ice[:,3]==k],weight[d.ice[:,3]==k]) for k,law in enumerate(laws)]
        rhoeq=air['ice_saturation_pressure_Pa']/(air['vapor_gas_constant_J_kg_K']*T)
        vth=air['vapor_molecular_speed_m_s']/4;D=air['diffusivity_m2_s']*ventilation;K=air['conductivity_W_m_K']*ventilation
        L=air['latent_sublimation_J_kg'];beta=(L/(air['vapor_gas_constant_J_kg_K']*T)-1)/T
        h=beta*L*rhoeq*D/K if heating else 0.
        ans=d.coupled(si,self.dx*vth/D,combined,h,guess=self.last,tolerance=tolerance)
        active=d.exterior_air
        # Growth-only stationary exterior: sealed air has no replenishment.
        # Freeze its partial ice and report an upper bound on ALL vapour there,
        # rather than supplying it from ambient humidity or filling it with ice.
        closed_volume=float(2*self.dx**3*np.sum(self.volume[d.closed_air]*(1-self.fraction[d.closed_air])))
        pocket=dict(closed_air_cells=int(np.sum(d.closed_air)),closed_gas_volume_m3=closed_volume,
            trapped_vapor_upper_bound_kg=rhoeq*(1+si)*closed_volume)
        if ans['half_flux']<=0:self.time=until;return dict(no_growth=True,**pocket)
        p=d.ice[:,0].astype(int)
        supply=np.bincount(p,weights=d.ice[:,1]*ans['face_flux'],minlength=d.cells)
        rate=rhoeq*D/(ICE_DENSITY*self.dx**2)*supply/self.volume[active]
        live=rate>0;remaining=1-self.fraction[active]
        # Compare inverse filling times, avoiding overflow for vanishingly
        # small attachment rates without imposing a physical coefficient floor.
        inverse_event=float(np.max(rate[live]/remaining[live]));dt=until-self.time
        if inverse_event*dt>1:dt=1/inverse_event
        updated=self.fraction[active]+dt*rate
        if np.max(updated)>1+1e-8:raise ArithmeticError('Growth event overfilled a cell')
        self.fraction[active]=np.minimum(updated,1.)
        self.deposited+=2*rhoeq*D*self.dx*ans['half_flux']*dt;self.time+=dt
        newly=(self.fraction>=1-1e-10)&~self.solid
        if np.any(newly):
            i,j=np.where(newly)
            if max(i)>=self.radial_core or max(j)>=self.vertical_core:raise ValueError('Crystal reached the uniform computational core')
            self.solid|=newly;self.fraction[newly]=1.;self.steps+=int(np.sum(newly))
        self.last=ans['guess']
        return dict(dt=dt,iterations=ans['iterations'],mass_residual=ans['relative_mass_residual'],
            nonlinear_residual=ans['relative_nonlinear_residual'],heat_shift=ans['heat_shift'],
            enhanced_faces=int(np.sum(enhance>0)),faces=len(enhance),cells=d.cells,**pocket)
