"""Earth shape/rate falsification controls for explicit face transport.

Prescribed measurements and differentiated fits are distinguished from forward
predictions. No attachment parameter, rim width or initial shape is fitted.
"""
from __future__ import annotations
import argparse
import itertools
import json
import math
from pathlib import Path
import numpy as np
from atmosphere.snow import uniform_face_growth as model, facet_attachment as attachment
from atmosphere.snow import facet_growth as baseline, hollow_shape as shape
from atmosphere.snow import hollow_growth_assessment as hollow, facet_coupling
from atmosphere.snow import takahashi_growth as lab, crystal_transport_replay as replay, flakes

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
RESULT=HERE/'results/uniform_face_assessment.json'
SOURCES=HERE/'uniform_face_sources.json'


def solid_mass_budget(a,c,da,dc,mass_rate,density=baseline.ICE_DENSITY):
    """Hexagonal envelope budget; a is half maximum basal span, c half length.

    Active-area fractions are kinematic equivalents under assumed uniform normal
    advance. The measured maximum width need not describe a whole crystal side.
    A negative fraction rejects that geometric interpretation of the supplied
    fitted rates; it is not evidence of negative mass or sublimation.
    """
    if any(not math.isfinite(x) or x<=0 for x in (a,c,da,dc,mass_rate,density)):
        raise ValueError('Positive finite dimensions, rates and density required')
    prism=6*math.sqrt(3)*a*c*da*density
    basal=3*math.sqrt(3)*a*a*dc*density
    f=(mass_rate-prism)/basal
    return dict(prism_growth_mass_rate_kg_s=prism,basal_growth_mass_rate_kg_s=basal,
        full_envelope_mass_rate_to_fit=(prism+basal)/mass_rate,
        full_prism_face_fraction_of_mass_budget=prism/mass_rate,
        basal_active_fraction_with_full_prism=f,
        thin_tip_only_mass_budget_feasible=0<=f<=1,
        equivalent_rim_width_fraction=1-math.sqrt(1-f) if 0<=f<=1 else None,
        uniform_active_area_fraction=mass_rate/(prism+basal))


def compare_one(source,shapes,a,c,t,si,da,dc,dm,**kwargs):
    try:
        answer=model.solve(source,shapes,a,c,t,si,**kwargs)
        fine=model.solve(source,shapes,a,c,t,si,tolerance=2e-15,**kwargs)
    except (ValueError,ArithmeticError) as exc:
        return dict(complete=False,error=str(exc),calculation=None)
    values=('equatorial_rate_m_s','polar_rate_m_s','mass_rate_kg_s')
    return dict(complete=True,error=None,calculation=fine,
        tolerance_refinement_relative_change={key:abs(answer[key]/fine[key]-1) for key in values},
        equatorial_rate_ratio_to_fit=fine['equatorial_rate_m_s']/da,
        polar_rate_ratio_to_fit=fine['polar_rate_m_s']/dc,
        directional_rate_ratio_to_fit=fine['facet_rate_ratio']/(dc/da),
        mass_rate_ratio_to_fit=fine['mass_rate_kg_s']/dm)


def measured_shape_cases(library,shapes):
    cases=[];budgets=[]
    for tc,(end,habit,*_) in lab.FITS.items():
        t=flakes.T0+tc;air=replay.air_state(tc);nu=air['viscosity_Pa_s']/air['density_kg_m3']
        for age in (4.5,8.,.9*end):
            obs=lab.observed(tc,age);a=obs['equatorial_dimension_m']/2;c=obs['polar_dimension_m']/2
            nm,na,nc=replay.fitted_exponents(tc,age)
            da=na*a/(60*age);dc=nc*c/(60*age);dm=nm*obs['mass_kg']/(60*age)
            source=library.evaluate(.1,t,101000.,a,c,10.)
            si=flakes.liquid_saturation_pressure(t)/source['ice_saturation_pressure_Pa']-1
            common=dict(temperature_C=tc,age_min=age,habit=habit,pressure_Pa=101000.,
                ice_supersaturation=si,observed_fit=obs,
                fitted_equatorial_semiaxis_rate_m_s=da,fitted_polar_semiaxis_rate_m_s=dc,
                fitted_mass_rate_kg_s=dm,source_air_and_attachment=source)
            budgets.append(dict(temperature_C=tc,age_min=age,**solid_mass_budget(a,c,da,dc,dm)))
            re=lab.fall_speed(tc,age)*2*a/nu
            factors=[('diffusion',1.,None)]
            try:
                factors.append(('JW',replay.jw_ventilation(re,.63,replay.branch_habit(tc,age)),None))
            except ValueError as exc:
                factors.append(('JW',None,str(exc)))
            for (flow,factor,error),sources,mechanisms in itertools.product(factors,
                    itertools.product(('edge','middle'),repeat=2),((10.,10.),(1.,10.),(10.,1.))):
                result=compare_one(source,shapes,a,c,t,si,da,dc,dm,sources=sources,mechanisms=mechanisms,
                    vapor_ventilation=factor,heat_ventilation=factor) if error is None else dict(
                    complete=False,error=error,calculation=None)
                cases.append(dict(**common,flow=flow,ventilation=factor,sources=list(sources),
                    mechanisms=list(mechanisms),thin_basal_rim=None,**result))
            if tc==-5.3:
                for (flow,factor,error),rim in itertools.product(factors,(.01,.03,.1)):
                    result=compare_one(source,shapes,a,c,t,si,da,dc,dm,thin_basal_rim=rim,
                        vapor_ventilation=factor,heat_ventilation=factor) if error is None else dict(
                        complete=False,error=error,calculation=None)
                    cases.append(dict(**common,flow=flow,ventilation=factor,sources=['edge','edge'],
                        mechanisms=[10.,10.],thin_basal_rim=rim,**result))
    return dict(cases=cases,mass_budgets=budgets,
        limits=['Fitted derivatives are not directly measured instantaneous rates; no fitted-derivative error bars exist.',
                'The cylinder transport approximation is explicit; mass-budget geometry uses the laboratory hexagonal envelope separately.',
                'Uniform-face failure is diagnosed at edges and centres, not every point; passing these necessary conditions is insufficient validation.',
                'Changing m with unchanged critical supersaturation is a sensitivity, not a calibrated dislocation model.',
                'JW is a supplied common vapor/heat conductance factor using source Sc=0.63, not a resolved faceted flow field.',
                'Thin basal rims use NB96 C1 only at the common corner; their off-edge field, hollow volume and rim evolution are unavailable.'])


def cold_constraints(library,shapes):
    members=shape.archive_members();chamber=hollow.chamber_controls(members);cases=[]
    for experiment in chamber['cases']:
        _,data=hollow.parse_chamber(members[experiment['member']])
        source=library.evaluate(.48,experiment['temperature_K'],experiment['pressure_Pa'],50e-6,50e-6,10.)
        windows=[]
        for win in experiment['windows']:
            start,end=win['requested_interval_min']
            selected=data[(data[:,0]>=60*start)&(data[:,0]<60*end)]
            # If da=dc at all times, c-a is constant. Fit its slope directly;
            # sum of supplied axis errors gives a conservative per-point bound.
            difference=hollow.line_fit(selected[:,0],selected[:,2]-selected[:,1],selected[:,6]+selected[:,5])
            windows.append(dict(requested_interval_min=[start,end],axis_growth_ratio=win['fitted_directional_growth_ratio'],
                axis_growth_ratio_measurement_interval=win['directional_ratio_interval'],
                polar_minus_equatorial_slope_um_s=difference,
                equal_axis_rate_excluded_by_reported_error_bounds=difference['interval'][0]>0 or difference['interval'][1]<0,
                rim_fraction_median=win['rim_fraction_median']))
        calculations=[]
        for row in data:
            si=.2 if experiment['switch_elapsed_s'] is not None and row[0]>=experiment['switch_elapsed_s'] else .48
            a,c=row[1]*1e-6,row[2]*1e-6
            state=library.evaluate(si,experiment['temperature_K'],experiment['pressure_Pa'],a,c,10.)
            result=model.solve(state,shapes,a,c,experiment['temperature_K'],si)
            calculations.append(dict(elapsed_s=float(row[0]),predicted_directional_rate_ratio=result['facet_rate_ratio'],
                passes_sampled_uniform_growth_conditions=result['passes_sampled_uniform_growth_conditions'],
                rim_fraction=float(row[7]/row[1]),inside_declared_thin_rim_regime=bool(row[7]/row[1]<=.1)))
        cases.append(dict(member=experiment['member'],member_sha256=experiment['member_sha256'],
            characteristic_supersaturations=[source['characteristic_supersaturation_a'],source['characteristic_supersaturation_c']],
            windows=windows,states=calculations))
    return dict(cases=cases,source='Pinned 2021 getscrit explicitly sets sa=sc below -30 C.',
        invariant='With identical mechanisms and characteristics at a common corner, da=dc for every size, humidity and common transport multiplier, including the prescribed thin-basal-rim approximation.',
        limits=['This excludes the inherited equal-face/common-corner hypothesis for the resolved windows, not the original 2021 spheroidal solver with separate surface humidities.',
                'Substrate experiments can alter step sources and heat flow; they do not directly validate free-falling crystals.',
                'Reported error bounds exclude measurement perturbations only; forcing/systematic errors are not quantified.'])


def summarize(measured,cold):
    rows=measured['cases'];complete=[x for x in rows if x['complete']]
    solid=[x for x in complete if x['thin_basal_rim'] is None]
    return dict(requests=len(rows),completed=len(complete),rejected=len(rows)-len(complete),
        completed_solid_cases=len(solid),solid_cases_failing_necessary_uniform_growth=sum(
            not x['calculation']['passes_sampled_uniform_growth_conditions'] for x in solid),
        completed_rim_sensitivities=len(complete)-len(solid),
        maximum_balance_residual=max(max(map(abs,x['calculation']['relative_balance_residuals'])) for x in complete),
        maximum_rate_refinement_change=max(max(x['tolerance_refinement_relative_change'].values()) for x in complete),
        maximum_mean_free_path_over_semiaxis=max(x['calculation']['maximum_mean_free_path_over_semiaxis'] for x in complete),
        cold_measurements=sum(len(x['states']) for x in cold['cases']),
        cold_windows_excluding_equal_axis_speeds=sum(w['equal_axis_rate_excluded_by_reported_error_bounds'] for x in cold['cases'] for w in x['windows']))


def run(output=RESULT):
    sources=json.loads(SOURCES.read_text())
    for item in sources['retrieved']:
        path=ROOT/item['cache_path']
        if shape.sha(path)!=item['sha256'] or path.stat().st_size!=item['bytes']:
            raise ValueError('Original paper identity changed: '+str(path))
    replay.verify_sources();hollow.read_sources()
    paths=[Path(__file__),Path(model.__file__),Path(facet_coupling.__file__),Path(lab.__file__),
           Path(replay.__file__),Path(flakes.__file__),Path(baseline.__file__),Path(hollow.__file__),SOURCES]
    producer={str(p.relative_to(ROOT)):shape.sha(p) for p in paths};producer.update(shape.producer())
    library=attachment.AttachmentLibrary(baseline.BUILD);shapes=shape.ShapeLibrary()
    measured=measured_shape_cases(library,shapes);cold=cold_constraints(library,shapes)
    result=dict(schema='terluna.atmosphere.uniform-face-assessment/1',producer=producer,sources=sources,
        attachment_build=library.record if hasattr(library,'record') else attachment.verify_build(baseline.BUILD),
        geometry_build=shapes.identity,measured_shapes=measured,cold_constraints=cold,
        summary=summarize(measured,cold),
        evidence='Source-based continuum falsification and prescribed-rim sensitivities; no autonomous Earth calibration or adopted lunar snow distribution.')
    for name,digest in producer.items():
        if shape.sha(ROOT/name)!=digest:raise RuntimeError('Producer changed during assessment')
    Path(output).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,default=RESULT)
    args=parser.parse_args();answer=run(args.output);print(json.dumps(answer['summary'],indent=2))
