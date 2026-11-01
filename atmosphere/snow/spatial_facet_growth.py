"""Prescribed cavity transport with explicit facet step-source hypotheses.

Own stationary diagnostic. Uniform normal uptake is imposed on the outer side
and top rim; the cavity interior is inactive. Full-field supply checks can reject
that hypothesis. This is not a moving-interface or surface-diffusion model.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.optimize import brentq
from atmosphere.snow.spatial_diffusion import Attachment, Domain
from atmosphere.snow.facet_growth import ICE_DENSITY


def libbrecht_parameters(temperature_C, *, narrow=True, logarithm='log10'):
    """2023 M1 (narrow) or M2 broad limit, Fig. 1 and pp. 6-7.

    T in the printed fits is positive undercooling, as the figure axis shows;
    sigma0 is in percent. Base-ten log reproduces the plotted dip widths. The
    unlabelled printed 'log' motivates an explicit natural-log sensitivity.
    The plotted interval .5--50 C is a declared input range, not validation.
    """
    if not math.isfinite(temperature_C) or not -.5 >= temperature_C >= -50:
        raise ValueError('Published plotted interval is -50 to -0.5 C')
    if logarithm not in ('log10','ln'):
        raise ValueError('Declare log10 or ln')
    t=-temperature_C
    basal=(.02*t**1.75+.3)/100
    prism=(.015*t*t+.02*t**.6)/100
    if narrow:
        log=math.log10 if logarithm=='log10' else math.log
        basal*=1-.87*math.exp(-log(t/4.5)**2/.07)
        prism*=1-.95*math.exp(-log(t/14.4)**2/.06)
    return dict(prism_critical=prism,basal_critical=basal,amplitude=1.,
                temperature_C=temperature_C,narrow=narrow,logarithm=logarithm)


def libbrecht_laws(temperature_C, *, narrow=True, logarithm='log10'):
    p=libbrecht_parameters(temperature_C,narrow=narrow,logarithm=logarithm)
    return [Attachment('nucleation',1.,p['prism_critical']),Attachment('nucleation',1.,p['basal_critical'])]


def _flux(law,s):
    return float(law.flux(np.array([s]))[0][0])


def surface_balances(matrix, supersaturation, laws, tolerance=1e-12):
    """Bracketed two-source balance; permits different local attachment laws."""
    matrix=np.asarray(matrix,float)
    if matrix.shape!=(2,2) or not np.all(np.isfinite(matrix)) or np.any(matrix<0) or np.any(np.diag(matrix)<=0):
        raise ValueError('Finite nonnegative resistance matrix with positive diagonals required')
    if matrix[0,0]*matrix[1,1]<matrix[0,1]*matrix[1,0]*(1-1e-12):
        raise ValueError('Selected source positions do not give a monotone two-face balance')
    if not math.isfinite(supersaturation) or supersaturation<0 or len(laws)!=2:
        raise ValueError('Nonnegative supersaturation and two laws required')
    if not math.isfinite(tolerance) or tolerance<=0:
        raise ValueError('Positive finite root tolerance required')
    if supersaturation==0:return np.zeros(2)
    def diagonal(driving,resistance,law):
        if driving<=0:return 0.
        x=brentq(lambda f:f+resistance*_flux(law,driving*f)/driving-1.,0.,1.,xtol=tolerance,rtol=tolerance)
        return driving*x
    def first(sc):
        drive=supersaturation-matrix[0,1]*_flux(laws[1],sc)
        if drive < -supersaturation*1e-10:raise ArithmeticError('Source consumes more than available supply')
        return diagonal(max(0.,drive),matrix[0,0],laws[0])
    high=diagonal(supersaturation,matrix[1,1],laws[1])
    if matrix[0,1]*_flux(laws[1],high)>supersaturation:
        high=brentq(lambda f:matrix[0,1]*_flux(laws[1],high*f)/supersaturation-1.,0.,1.,xtol=tolerance,rtol=tolerance)*high
    def residual(f):
        sc=high*f;sa=first(sc)
        return (sc+matrix[1,0]*_flux(laws[0],sa)+matrix[1,1]*_flux(laws[1],sc))/supersaturation-1
    low_value,high_value=residual(0.),residual(1.)
    if low_value>1e-10 or high_value< -1e-10:
        raise ValueError('No nonnegative solution at selected step sources')
    x=0. if low_value>=0 else (1. if high_value<=0 else brentq(residual,0.,1.,xtol=tolerance,rtol=tolerance))
    sc=high*x;sa=first(sc);ss=np.array([sa,sc])
    p=np.array([_flux(law,s) for law,s in zip(laws,ss)])
    if np.max(np.abs(ss+matrix@p-supersaturation))/supersaturation>1e-8:
        raise ArithmeticError('Two-face source balance did not close')
    return ss


class FacetGeometry:
    def __init__(self,domain):
        if domain.shape!='cylinder':raise ValueError('Cylinder envelope required')
        self.domain=domain
        d=domain
        self.side=(d.ice[:,3]==0)&np.isclose(d.ice[:,4],1.,atol=1e-13,rtol=0.)
        self.top=(d.ice[:,3]==1)&np.isclose(d.ice[:,5],d.aspect,atol=1e-13,rtol=0.)
        self.masks=[self.side,self.top]
        self.areas=np.array([np.sum(d.ice[x,1]) for x in self.masks])
        self.capacity=d.dirichlet()['half_conductance']
        solutions=[d.imposed_flux(face_flux=x.astype(float)) for x in self.masks]
        self.maximum_mass_residual=max(x['relative_mass_residual'] for x in solutions)
        self.fields=np.array([x['surface'] for x in solutions]).T
        self.source_indices={}
        for name,mask,r,z in [('prism_middle',self.side,1.,0.),('basal_middle',self.top,(1+d.cavity_radius)/2,d.aspect),
                              ('prism_edge',self.side,1.,d.aspect),('basal_edge',self.top,1.,d.aspect)]:
            ids=np.where(mask)[0];self.source_indices[name]=int(ids[np.argmin((d.ice[ids,4]-r)**2+(d.ice[ids,5]-z)**2)])
        # The end-face "middle" is the centre of its radial width, not the
        # cavity floor. For solid disks the actual face centre is r=0.
        if not d.cavity_radius:
            ids=np.where(self.top)[0];self.source_indices['basal_middle']=int(ids[np.argmin(d.ice[ids,4])])

    def source_lengths(self,sources):
        if len(sources)!=2 or any(x not in ('edge','middle') for x in sources):raise ValueError('Two edge/middle step sources required')
        rows=[]
        # A common-edge field is sampled from the active face of each imposed
        # flux separately. Both approach the same physical corner on refinement.
        common=np.array([self.fields[self.source_indices['prism_edge'],0],self.fields[self.source_indices['basal_edge'],1]])
        for face,source in zip(('prism','basal'),sources):
            rows.append(common if source=='edge' else self.fields[self.source_indices[face+'_middle']])
        return np.array(rows)

    def solve(self,state,a,temperature,supersaturation,laws, *, sources=('edge','edge'), vapor_ventilation=1.,heat_ventilation=1., heating=True,tolerance=1e-12):
        if any(not math.isfinite(x) or x<=0 for x in (a,temperature,vapor_ventilation,heat_ventilation)):
            raise ValueError('Positive finite physical inputs required')
        rhoeq=state['ice_saturation_pressure_Pa']/(state['vapor_gas_constant_J_kg_K']*temperature)
        vth=state['vapor_molecular_speed_m_s']/4;latent=state['latent_sublimation_J_kg']
        beta=(latent/(state['vapor_gas_constant_J_kg_K']*temperature)-1)/temperature
        D=state['diffusivity_m2_s']*vapor_ventilation;K=state['conductivity_W_m_K']*heat_ventilation
        thermal=rhoeq*vth*a*latent*beta/K*self.areas/self.capacity if heating else np.zeros(2)
        matrix=a*vth/D*self.source_lengths(sources)+thermal[None,:]
        ss=surface_balances(matrix,supersaturation,laws,tolerance)
        p=np.array([_flux(law,s) for law,s in zip(laws,ss)])
        shift=float(np.dot(thermal,p));field=supersaturation-a*vth/D*(self.fields@p)-shift
        checks=[]
        for i,mask in enumerate(self.masks):
            minimum=float(np.min(field[mask]));required=p[i]/field[mask] if minimum>0 else None
            checks.append(dict(face=('prism','basal')[i],minimum_surface_supersaturation=minimum,
                maximum_required_efficiency=float(np.max(required)) if required is not None else None,
                passes_uniform_local_uptake_condition=bool((p[i]==0 and minimum>=0) or (required is not None and np.max(required)<=1+1e-9))))
        residual=ss+matrix@p-supersaturation
        return dict(prism_rate_m_s=float(rhoeq*vth*p[0]/ICE_DENSITY),basal_rate_m_s=float(rhoeq*vth*p[1]/ICE_DENSITY),
            mass_rate_kg_s=float(2*a*a*rhoeq*vth*np.dot(self.areas,p)),source_supersaturations=ss.tolist(),
            facet_rate_ratio=float(p[1]/p[0]) if p[0] else None,
            source_balance_relative_residual=float(np.max(abs(residual))/supersaturation) if supersaturation else 0.,
            temperature_increment_K=shift/beta,minimum_surface_supersaturation=float(np.min(field)),
            inactive_minimum_surface_supersaturation=float(np.min(field[~(self.side|self.top)])) if np.any(~(self.side|self.top)) else None,
            surface_negative_area_fraction=float(np.sum(self.domain.ice[field<0,1])/np.sum(self.domain.ice[:,1])),
            uniform_checks=checks,passes_uniform_local_uptake_conditions=all(x['passes_uniform_local_uptake_condition'] for x in checks),
            source_resistance_matrix=matrix.tolist(),half_active_areas_over_radius_squared=self.areas.tolist(),
            transport_relative_mass_residual=self.maximum_mass_residual,
            source_positions=list(sources),heating=heating,vapor_ventilation=vapor_ventilation,heat_ventilation=heat_ventilation)
