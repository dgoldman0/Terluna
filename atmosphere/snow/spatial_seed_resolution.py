"""Frozen initial-profile resolution test for the local warm growth treatment.

Surface fluxes from one event step are recovered on its saved initial domain
and checked against the independent integrated deposition ledger. Even if a
cell fills, the flux belongs to the initial shape. This resolves local uptake,
not a full history.
"""
from __future__ import annotations
import json
import math
from pathlib import Path
import time
import numpy as np
from atmosphere.snow import facet_attachment, facet_growth, flakes
from atmosphere.snow import takahashi_growth as lab, crystal_transport_replay as replay
from atmosphere.snow.spatial_crystal_growth import GrowingCrystal, MaskDomain, Combined
from atmosphere.snow.spatial_diffusion import surface_flux, Attachment
from atmosphere.snow.spatial_facet_growth import libbrecht_laws
from atmosphere.snow.spatial_crystal_assessment import file_id

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
RESULT=HERE/'results/spatial_seed_resolution.json'


def run():
    library=facet_attachment.AttachmentLibrary(facet_growth.BUILD)
    tc,age=-5.3,3.
    T=flakes.T0+tc
    obs=lab.observed(tc,age)
    a,c=obs['equatorial_dimension_m']/2,obs['polar_dimension_m']/2
    air=library.evaluate(.1,T,101000.,a,c,10.)
    si=flakes.liquid_saturation_pressure(T)/air['ice_saturation_pressure_Pa']-1
    flow=replay.air_state(tc)
    nu=flow['viscosity_Pa_s']/flow['density_kg_m3']
    fv=replay.jw_ventilation(lab.fall_speed(tc,age)*2*a/nu,.63,replay.branch_habit(tc,age))
    rhoeq=air['ice_saturation_pressure_Pa']/(air['vapor_gas_constant_J_kg_K']*T)
    D=air['diffusivity_m2_s']*fv
    nm,na,nc=replay.fitted_exponents(tc,age)
    target=np.array([na*a,nc*c,nm*obs['mass_kg']])/(60*age)
    rows=[]
    cases=[(dx,'M1') for dx in (1e-6,.5e-6,.25e-6,.125e-6,.0625e-6)]
    cases += [(.125e-6,mode) for mode in ('fast_basal','inactive_prism_fast_basal','fast_all')]
    for dx,mode in cases:
        start=time.monotonic()
        g=GrowingCrystal(dx,a,c,math.ceil(150e-6/dx),math.ceil(1000e-6/dx),extent=80.)
        seed=g.measured_seed(a,c,obs['mass_kg'],2e-6)
        initial=g.state()
        laws=libbrecht_laws(tc)
        if mode != 'M1':
            laws[1]=Attachment('constant',1.)
        if mode == 'inactive_prism_fast_basal':
            laws[0]=Attachment('constant',1e-12)
        if mode == 'fast_all':
            laws[0]=Attachment('constant',1.)
        d=MaskDomain(g.rf,g.zf,g.solid.copy())
        status=g.advance(air,T,si,laws,1e-5,kinks=False,
            capillary_length=1e-9,ventilation=fv)
        if not (np.array_equal(d.rf,g.rf) and np.array_equal(d.zf,g.zf)):
            raise ArithmeticError('Probe grid expanded unexpectedly')
        p=d.ice[:,0].astype(int);axis=d.ice[:,3].astype(int);area=d.ice[:,1]
        neighbors=np.bincount(p,minlength=d.cells)
        weight=1/np.sqrt(neighbors[p])
        curvature=d.curvature_offsets(1e-9/dx)
        lam=dx*(air['vapor_molecular_speed_m_s']/4)/D
        u,heat=g.last
        center=si-heat-u[d.exterior_air][p]
        q=np.zeros(len(p))
        for k,law in enumerate(laws):
            mask=axis==k
            boundary=Combined(law,0.,curvature[mask],weight[mask])
            q[mask]=surface_flux(center[mask],d.ice[mask,2],lam,boundary)[0]
        total=2*rhoeq*D*dx*float(np.dot(area,q))
        integrated=g.deposited/g.time
        residual=abs(total/integrated-1)
        if residual>1e-10:
            raise ArithmeticError('Recovered surface field differs from integrated uptake')
        rates=[]
        for k in (0,1):
            chosen=axis==k
            outer=chosen&(d.ice[:,4+k]==max(d.ice[chosen,4+k]))
            rates.append(rhoeq*D/(facet_growth.ICE_DENSITY*dx)*float(np.dot(area[outer],q[outer])/sum(area[outer])))
        rates.append(total)
        row=dict(dx_m=dx,kinetics=mode,seed=seed,filled_initial=initial,cells=d.cells,
            probe_elapsed_s=g.time,probe_filled_cell_events=g.steps,
            rates_a_c_mass=rates,ratios_to_fitted_rates=(np.array(rates)/target).tolist(),
            surface_to_integrated_mass_relative_residual=residual,step=status,
            wall_elapsed_s=time.monotonic()-start)
        rows.append(row)
        print(json.dumps(row),flush=True)
    return dict(schema='terluna.atmosphere.spatial-seed-resolution/1',
        evidence='local stationary flux and seed-discretization control; no trajectory or morphology validation',
        inputs=[file_id(p) for p in [Path(__file__).resolve(),HERE/'spatial_crystal_growth.py',
            HERE/'spatial_crystal_assessment.py',HERE/'spatial_diffusion.py',HERE/'spatial_facet_growth.py',
            HERE/'facet_attachment.py',HERE/'facet_growth.py',HERE/'takahashi_growth.py',
            HERE/'crystal_transport_replay.py',ROOT/'shared/constants.json']],
        native_build_identity=library.identity,temperature_K=T,pressure_Pa=101000.,
        ice_supersaturation=si,observed=obs,air=air,ventilation=fv,fitted_rates_a_c_mass=target.tolist(),cases=rows,
        limits=[
            'Prescribed 2 micrometre seed rim and exact fitted initial mass; partial boundary cells change the gas-interface envelope with resolution.',
            'Local M1 kinetics with the diagonal-area weight and supplied curvature proxy; no terrace completion or adjoining-facet transfer.',
            'The added controls set basal attachment to one, optionally prism attachment to one or an explicit near-inactive 1e-12. These are limiting hypotheses, not measured parameters or a coefficient floor.',
            'A fully absorbing local face does not bound adjoining-facet redistribution or a different three-dimensional shape.',
            'Reported side and top rates are area means on the outermost filled-cell terraces, not maximum-axis trajectory predictions.',
            'One frozen initial geometry cannot establish convergence of an evolving history or validate the unmeasured seed.'])


if __name__=='__main__':
    result=run()
    RESULT.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
