"""Conservative frozen-profile surface-transfer sensitivities and prism steps.

Unknown prism migration length, edge transmission, cavity and rim geometry are
varied explicitly. Agreement at a supplied geometry is not an autonomous fit.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from atmosphere.snow import facet_attachment, facet_growth, flakes
from atmosphere.snow import takahashi_growth as lab, crystal_transport_replay as replay
from atmosphere.snow.spatial_diffusion import Domain
from atmosphere.snow.spatial_facet_growth import libbrecht_laws
from atmosphere.snow.surface_transfer import transfer_rates, prism_spiral_rate

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
RESULT=HERE/'results/surface_transfer_assessment.json'


def identity(path):
    return dict(path=str(path.relative_to(ROOT)),sha256=facet_growth.sha(path))


def observed_state(library,age):
    tc=-5.3;T=flakes.T0+tc
    obs=lab.observed(tc,age)
    a,c=obs['equatorial_dimension_m']/2,obs['polar_dimension_m']/2
    air=library.evaluate(.1,T,101000.,a,c,10.)
    si=flakes.liquid_saturation_pressure(T)/air['ice_saturation_pressure_Pa']-1
    flow=replay.air_state(tc);nu=flow['viscosity_Pa_s']/flow['density_kg_m3']
    fv=replay.jw_ventilation(lab.fall_speed(tc,age)*2*a/nu,.63,replay.branch_habit(tc,age))
    nm,na,nc=replay.fitted_exponents(tc,age)
    fitted=np.array([na*a,nc*c,nm*obs['mass_kg']])/(60*age)
    return dict(age_min=age,temperature_C=tc,temperature_K=T,pressure_Pa=101000.,
        observed=obs,air=air,ice_supersaturation=si,ventilation=fv,
        radius_m=a,half_length_m=c,fitted_rates_a_c_mass=fitted.tolist())


def calculate(state,rim,xs,transmission,n=32,extent=80.,inner=True,capillary_length=0.,tolerance=1e-9):
    a,c=state['radius_m'],state['half_length_m']
    d=Domain(c/a,resolution=n,extent=extent,mesh='edge',cavity_radius=1-rim,cavity_depth=.5*c/a)
    row=dict(age_min=state['age_min'],rim_fraction=rim,migration_length_m=xs,
        transmission=transmission,resolution=n,extent=extent,inner_donor=inner,
        capillary_length_m=capillary_length,nonlinear_tolerance=tolerance)
    try:
        rates=transfer_rates(d,state['air'],a,state['temperature_K'],state['ice_supersaturation'],
            libbrecht_laws(state['temperature_C']),xs,transmission,inner_donor=inner,
            ventilation=state['ventilation'],capillary_length=capillary_length,tolerance=tolerance)
        fitted=state['fitted_rates_a_c_mass']
        row.update(complete=True,rates=rates,
            ratios_mean_prism_tip_mass=[rates['outer_prism_mean_normal_rate_m_s']/fitted[0],
                rates['tip_mean_normal_rate_m_s']/fitted[1],rates['mass_rate_kg_s']/fitted[2]],
            ratios_max_prism_tip=[rates['outer_prism_max_normal_rate_m_s']/fitted[0],
                rates['tip_max_normal_rate_m_s']/fitted[1]],
            fitted_tip_fraction_of_fitted_mass=facet_growth.ICE_DENSITY*rates['tip_area_m2']*fitted[1]/fitted[2])
    except (ValueError,ArithmeticError) as exc:
        row.update(complete=False,error=str(exc))
    return row


def spiral_diagnostics(states):
    rows=[]
    for state in states:
        for si in (.001,.005,.01,.02):
            for kappa in (1.5e-10,8.5e-10):
                for activation in (43e3,44e3,45e3):
                    r=prism_spiral_rate(state['temperature_K'],si,kappa,activation)
                    width_rate=state['fitted_rates_a_c_mass'][0]
                    rows.append(dict(age_min=state['age_min'],surface_supersaturation=si,
                        ledge_energy_J_m=kappa,activation_energy_J_mol=activation,**r,
                        normal_rate_to_fitted_width_rate=r['normal_rate_m_s']/width_rate,
                        step_spacing_needed_for_fitted_width_m=r['lateral_step_velocity_m_s']*.39e-9/width_rate))
    return dict(cases=rows,limits=[
        'Source: Miyamoto et al. 2022 prism spiral steps in nitrogen on a substrate, at low supersaturation.',
        'The two kappa values span the whole measured temperature range; they are not a temperature-specific confidence interval.',
        'The 43-45 kJ/mol sensitivity comes from the reported 44 +/- 1 kJ/mol activation fit.',
        'Surface supersaturation is supplied, not inferred from far-field humidity; no gas-phase diffusion correction is double counted.',
        'A large fitted width discrepancy does not exclude other step sources, dense vicinal terraces, nucleation or edge feeding.'])


def run():
    library=facet_attachment.AttachmentLibrary(facet_growth.BUILD)
    states=[observed_state(library,age) for age in (4.5,8.,27.)]
    cases=[]
    for state in states:
        for rim in (.03,.1,.3):
            for xs in (0.,14e-9,.1e-6,1e-6,5e-6,10e-6):
                for transmission in ((1.,) if xs==0 else (.1,1.)):
                    row=calculate(state,rim,xs,transmission,n=16)
                    cases.append(row)
            print(json.dumps(dict(age=state['age_min'],rim=rim,completed=len(cases),failures=sum(not r['complete'] for r in cases))),flush=True)
    refined=[]
    state=states[1]
    for xs in (0.,14e-9,1e-6,5e-6):
        for n in (16,32,64,128):
            refined.append(calculate(state,.1,xs,1.,n=n))
    boundaries=[calculate(state,.1,xs,1.,n=64,extent=extent)
        for xs in (0.,1e-6,5e-6) for extent in (40.,80.,160.)]
    donors=[calculate(state,.1,xs,1.,n=64,inner=inner)
        for xs in (1e-6,5e-6,10e-6) for inner in (False,True)]
    curvature=[calculate(state,rim,xs,1.,n=64,capillary_length=delta)
        for rim in (.03,.1) for xs in (0.,14e-9,1e-6,5e-6,10e-6)
        for delta in (.3e-9,1e-9)]
    roots=[calculate(state,.1,xs,1.,n=n,tolerance=tol)
        for xs in (0.,5e-6) for n in (16,32) for tol in (1e-9,1e-11)]
    allrows=cases+refined+boundaries+donors+curvature+roots
    return dict(schema='terluna.atmosphere.surface-transfer-assessment/1',
        evidence='conditional frozen-profile transfer and source-parameter diagnostics; no autonomous morphology validation',
        inputs=[identity(p) for p in [Path(__file__).resolve(),HERE/'surface_transfer.py',HERE/'surface_transport_sources.json',
            HERE/'spatial_crystal_growth.py',HERE/'spatial_diffusion.py',HERE/'spatial_facet_growth.py',
            HERE/'facet_attachment.py',HERE/'facet_growth.py',HERE/'takahashi_growth.py',HERE/'crystal_transport_replay.py',
            ROOT/'shared/constants.json']],
        states=states,cases=cases,refinement=refined,boundaries=boundaries,donor_sensitivity=donors,
        curvature_sensitivity=curvature,root_tolerance_sensitivity=roots,
        prism_spiral=spiral_diagnostics(states),
        completed_cases=sum(x['complete'] for x in allrows),failed_cases=sum(not x['complete'] for x in allrows),
        limits=[
            'Rims and half-depth cavities are supplied and unmeasured for this Earth history.',
            'An absorbing rough tip and exponential capture with supplied edge transmission define the hypothesis; the surface molecular field is not solved.',
            'A nanometre basal migration estimate cannot be transferred to prism faces as a measured parameter.',
            'The micrometre cases are deliberately broad sensitivities motivated by historical basal studies, not empirical prism bounds.',
            'Outer-prism mean motion differs from the maximum span; the maximum local rate is only an instantaneous envelope rate on this prescribed cylinder.',
            'Uniform redistribution over the annulus does not predict step propagation or rim shape.',
            'Main cases neglect capillarity. Separate 0.3 and 1 nm sensitivities use the source-style 2/rim_width curvature proxy and the same equilibrium shift for direct and donor-fed growth.',
            'The curvature sensitivity is growth-only and does not solve reversible surface diffusion or sublimation.',
            'Requested relative nonlinear tolerance is 1e-9, above measured graded-matrix roundoff; coarser 1e-11 controls distinguish solver tolerance from mesh error.',
            'Known whole-fit mass/axis inconsistencies and absent derivative errors prevent treating any snapshot agreement as validation.',
            'Nitrogen-pressure evidence is registered but no unreviewed pressure multiplier is implemented.'])


if __name__=='__main__':
    result=run()
    RESULT.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('completed_cases','failed_cases')},indent=2))
