"""Matched growth histories and measured-shape controls for implicit kinetics.

The low-excess-vapor experiments use dense growth only. Archive comparisons are
numerical controls; the humidity diameters are laboratory evidence with unknown
seed shapes. Prescribed-shape rate comparisons avoid assuming those seed shapes.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import zipfile
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from atmosphere.snow import facet_attachment as attachment, facet_growth as baseline
from atmosphere.snow import facet_coupling as coupled, facet_coupling_audit as audit
from atmosphere.snow import growth_humidity_validation as humidity
from atmosphere.snow import crystal_transport_replay as replay, takahashi_growth as lab
from atmosphere.snow import flakes

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
RESULT=HERE/'results/facet_growth_controls.json'
SOURCES=HERE/'facet_growth_controls_sources.json'
DENSITY=baseline.ICE_DENSITY


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def integrate(library,temperature,pressure,si,a0,c0,mechanism,method,*,rtol=1e-8,ventilation=1.):
    """Fixed forcing, dense spheroid, independent mass integral, no fitted sizes."""
    if method not in ('semianalytic','implicit'):
        raise ValueError('Unknown attachment method')
    if method=='semianalytic' and ventilation!=1.:
        raise ValueError('Native approximation has no ventilation input')
    if any(not math.isfinite(x) or x<=0 for x in (temperature,pressure,si,a0,c0,rtol,ventilation)):
        raise ValueError('Positive finite integration inputs required')
    m0=DENSITY*4*math.pi/3*a0*a0*c0
    maximum_balance=0.

    def evaluate(a,c):
        nonlocal maximum_balance
        source=library.evaluate(si,temperature,pressure,a,c,mechanism)
        if method=='implicit':
            state=coupled.solve(source,a,c,temperature,si,mechanism,
                vapor_ventilation=ventilation,heat_ventilation=ventilation)
            rates=coupled.axis_rates(state,a,c,DENSITY)
            maximum_balance=max(maximum_balance,abs(state['relative_thermal_balance_residual']),
                                *map(abs,state['relative_surface_balance_residuals']))
            rates.update(mass_rate_kg_s=state['mass_rate_kg_s'],alpha_a=state['alpha_a'],alpha_c=state['alpha_c'])
        else:
            rates=baseline.mass_and_axis_rates(source,a,c,temperature,si,DENSITY)
            rates.update(alpha_a=source['alpha_a'],alpha_c=source['alpha_c'])
        return rates

    def rhs(seconds,y):
        a=a0*math.exp(y[0]);c=c0*math.exp(y[1]);rates=evaluate(a,c)
        return [rates['equatorial_rate_m_s']/a,rates['polar_rate_m_s']/c,rates['mass_rate_kg_s']/m0]

    solution=solve_ivp(rhs,(0.,600.),[0.,0.,0.],rtol=rtol,atol=rtol*.01,
                       t_eval=np.arange(0.,601.,60.),max_step=5.)
    rows=[]
    for seconds,y in zip(solution.t,solution.y.T):
        a=a0*math.exp(y[0]);c=c0*math.exp(y[1]);mass=DENSITY*4*math.pi/3*a*a*c
        rates=evaluate(a,c)
        rows.append(dict(age_s=float(seconds),equatorial_dimension_m=2*a,polar_dimension_m=2*c,
            diameter_mm=2*max(a,c)*1000,mass_kg=mass,
            integrated_deposition_kg=float(y[2])*m0,
            relative_mass_volume_residual=(mass-m0-float(y[2])*m0)/mass,
            alpha_a=rates['alpha_a'],alpha_c=rates['alpha_c']))
    complete=bool(solution.success and rows and rows[-1]['age_s']==600.)
    return dict(temperature_K=temperature,pressure_Pa=pressure,ice_supersaturation=si,
        initial_equatorial_radius_m=a0,initial_polar_radius_m=c0,initial_mass_kg=m0,
        mechanism_m=mechanism,method=method,ventilation=ventilation,rtol=rtol,
        deposition_density_kg_m3=DENSITY,complete=complete,solver_message=solution.message,
        function_evaluations=solution.nfev,trajectory=rows,endpoint=rows[-1] if complete else None,
        maximum_sampled_mass_volume_residual=max(abs(row['relative_mass_volume_residual']) for row in rows),
        maximum_rhs_balance_residual=maximum_balance)


def refined_pair(library,*args,**kwargs):
    coarse=integrate(library,*args,rtol=1e-6,**kwargs)
    fine=integrate(library,*args,rtol=1e-8,**kwargs)
    both=coarse['complete'] and fine['complete']
    metrics={key:abs(coarse['endpoint'][key]/fine['endpoint'][key]-1) if both else None
             for key in ('equatorial_dimension_m','polar_dimension_m','mass_kg')}
    return fine,dict(complete_pair=both,coarse_evaluations=coarse['function_evaluations'],
                     coarse_tolerance=1e-6,fine_tolerance=1e-8,relative_endpoint_changes=metrics)


def author_growth_controls(library):
    record=audit.verify_sources()
    item=next(x for x in record['retrieved'] if x['cache_path'].endswith('.zip'))
    with zipfile.ZipFile(ROOT/item['cache_path']) as archive:
        data,identity=audit.archive_array(archive,'Figure4.output')
        humidities,humidity_identity=audit.archive_array(archive,'Figure5b.output')
    if data.shape!=(100,13) or humidities.shape!=(100,7):
        raise ValueError('Unexpected archived growth table')
    cases=[]
    # Fixed logarithmically spaced low-supersaturation indices, selected before
    # comparing the solver. Only compact-growth conditions are included.
    for index in (0,20,40,60,75):
        row=data[index];si=float(row[0])/100
        source=library.evaluate(si,258.15,50000.,10e-6,10e-6,10.)
        excess=source['ice_saturation_pressure_Pa']/(source['vapor_gas_constant_J_kg_K']*258.15)*si*1000
        if excess>.05:
            raise ValueError('Compact control exceeds selected excess-vapor range')
        for mechanism,offset in ((10.,1),(1.,5),(3.,9)):
            for method in ('semianalytic','implicit'):
                fine,refinement=refined_pair(library,258.15,50000.,si,10e-6,10e-6,mechanism,method)
                reference_a=float(row[offset])*1e-6;reference_c=float(row[offset+1])*1e-6
                cases.append(dict(archive_index=index,archived_a_semiaxis_m=reference_a,
                    archived_c_semiaxis_m=reference_c,vapor_density_excess_g_m3=excess,
                    declared_characteristic_supersaturations=[float(x)/100 for x in humidities[index,5:7]],
                    calculation=fine,refinement=refinement,
                    equatorial_ratio_to_archived_accurate=fine['endpoint']['equatorial_dimension_m']/(2*reference_a) if fine['complete'] else None,
                    polar_ratio_to_archived_accurate=fine['endpoint']['polar_dimension_m']/(2*reference_c) if fine['complete'] else None))
    return dict(cases=cases,member_identity=identity,humidity_member_identity=humidity_identity,
        evidence='Numerical reference comparison; source forcing follows the later manuscript and archive labels. No Earth measurements in this subsection.',
        limits=['The exact archived executable and timestep are unavailable; complete numerical reproduction is unestablished.',
                'The separate static-curve archive has an unresolved mechanism mismatch; these growth tables are not silently substituted for it.',
                'Source 2019/2021 critical supersaturations retained; no parameter fitting. Dense spheroid uses 919 kg/m3 and no ventilation.'])


def humidity_controls(library):
    data=humidity.observations();cases=[]
    for datum in [row for row in data['cases'] if row['vapor_density_excess_g_m3']<=.05]:
        t=flakes.T0+datum['nominal_temperature_C']
        e=flakes.ice_saturation_pressure(t)+datum['vapor_density_excess_g_m3']*1e-3*flakes.R_VAPOUR*t
        reference=library.evaluate(.01,t,101000.,5e-6,5e-6,10.)
        si=e/reference['ice_saturation_pressure_Pa']-1
        for mechanism,phi,(method,ventilation) in itertools.product((1.,10.),(.1,.3,1.),
                (('semianalytic',1.),('implicit',1.),('implicit',2.))):
            fine,refinement=refined_pair(library,t,101000.,si,5e-6,phi*5e-6,mechanism,method,ventilation=ventilation)
            cases.append(dict(datum_id=datum['id'],observed_diameter_mm=datum['observed_diameter_mm'],
                vapor_density_excess_g_m3=datum['vapor_density_excess_g_m3'],nominal_temperature_C=datum['nominal_temperature_C'],
                initial_thickness_over_span=phi,imposed_vapor_pressure_Pa=e,calculation=fine,refinement=refinement,
                diameter_ratio_to_observed=fine['endpoint']['diameter_mm']/datum['observed_diameter_mm'] if fine['complete'] else None))
    return dict(cases=cases,source_observations=data,
        evidence='Earth diameter comparison with explicit unknown seed shape and flow scenarios; no selected initial-shape or mechanism fit.',
        limits=['Full 10 micrometre initial basal spans; varying thickness also changes initial mass.',
                'Constant dense ice restricted to observed excess vapor <= 0.05 g/m3; branching/hollowing absent.',
                'The factor-two conductance case is imposed; no observed fall-speed history supports it as an experimental reconstruction.',
                'Nominal temperature groups differ from some individual crystal temperatures by several tenths of a degree.',
                'No measured mass/thickness and no exact seed shape; matching a diameter alone cannot validate autonomous growth.',
                'Two selected cases are within the newer empirical morphology-onset range; the dense assumption is a restricted control, not guaranteed observed habit.'])


def measured_shape_rates(library):
    verification=replay.verify_sources();cases=[]
    for temperature,(end,habit,*_) in lab.FITS.items():
        t=flakes.T0+temperature;air=replay.air_state(temperature)
        nu=air['viscosity_Pa_s']/air['density_kg_m3']
        for age in (4.5,8.,.9*end):
            obs=lab.observed(temperature,age)
            nm,na,nc=replay.fitted_exponents(temperature,age)
            a=obs['equatorial_dimension_m']/2;c=obs['polar_dimension_m']/2
            source=library.evaluate(.1,t,101000.,a,c,10.)
            si=flakes.liquid_saturation_pressure(t)/source['ice_saturation_pressure_Pa']-1
            observed_q=nc*c/(na*a);observed_mass_rate=nm*obs['mass_kg']/(60*age)
            inferred_density=obs['mass_kg']/(4*math.pi/3*a*a*c)*nm/(2*na+nc)
            re=lab.fall_speed(temperature,age)*2*a/nu
            for method in ('diffusion','JW'):
                try:
                    factor=1. if method=='diffusion' else replay.jw_ventilation(re,.63,replay.branch_habit(temperature,age))
                    calculation=coupled.solve(source,a,c,t,si,10.,vapor_ventilation=factor,heat_ventilation=factor)
                    error=None
                except ValueError as exc:
                    calculation=None;factor=None;error=str(exc)
                cases.append(dict(temperature_C=temperature,age_min=age,habit=habit,method=method,
                    observed_fit=obs,observed_fitted_facet_rate_ratio=observed_q,
                    observed_fitted_mass_rate_kg_s=observed_mass_rate,
                    inferred_deposition_density_spheroid_kg_m3=inferred_density,
                    inferred_deposition_density_hexagonal_prism_kg_m3=inferred_density/replay.HEX_TO_SPHEROID,
                    inferred_spheroid_deposition_density_admissible=inferred_density<=DENSITY,
                    reynolds_basal=re,ventilation=factor,ventilation_schmidt=.63 if method=='JW' else None,
                    complete=calculation is not None,domain_error=error,calculation=calculation,
                    facet_rate_ratio_to_fit=(calculation['alpha_c']/calculation['alpha_a'])/observed_q if calculation else None,
                    mass_rate_ratio_to_fit=calculation['mass_rate_kg_s']/observed_mass_rate if calculation else None))
    return dict(cases=cases,source_verification=verification,
        evidence='Instantaneous rate diagnostic on supplied observed-fit shapes and falling speeds; derivatives of fits are not new directly measured rates.',
        limits=['Spheroid capacitance and critical-supersaturation approximation are retained; branches are not spatially resolved.',
                'JW source Sc = 0.63 is used and the same factor is assigned to heat and vapor; this is an explicit common-factor sensitivity.',
                'No extrapolation outside the chosen JW Reynolds interval and no unsupported habit substitution; rejected cases have no calculated endpoint.',
                'Planar thickness depends on mass/area, so shape and mass are not three independent measurements.',
                'Inferred deposition densities are diagnostics in both geometries, not adopted laws or silently clipped inputs.'])


def baseline_parity(empirical):
    """Require that moving the baseline into this integrator preserves all 24 cases."""
    path=baseline.RESULT;reference=json.loads(path.read_text())
    for name,digest in reference['producer'].items():
        if sha(ROOT/name)!=digest:
            raise ValueError('Frozen compact baseline producer differs: '+name)
    errors=[]
    for row in empirical['cases']:
        current=row['calculation']
        if current['method']!='semianalytic':continue
        old=next(item for item in reference['cases'] if item['datum_id']==row['datum_id']
            and item['mechanism_m']==current['mechanism_m']
            and item['initial_thickness_over_span']==row['initial_thickness_over_span'])
        errors.extend([abs(old['endpoint_diameter_mm']/current['endpoint']['diameter_mm']-1),
                       abs(old['endpoint_mass_kg']/current['endpoint']['mass_kg']-1)])
    if len(errors)!=48 or max(errors)>1e-12:
        raise ValueError('Semianalytic trajectory baseline was not preserved')
    return dict(cases=24,baseline_product_sha256=sha(path),maximum_relative_mass_or_diameter_difference=max(errors),
        evidence='Existing stagnant compact-crystal endpoints reproduced; empirical accuracy is separate.')


def morphology_screen():
    """Apply the 2024 measured onset range; no density law or automatic switch."""
    record=json.loads(SOURCES.read_text())
    for item in record['retrieved']:
        path=ROOT/item['cache_path']
        if sha(path)!=item['sha256'] or path.stat().st_size!=item['bytes']:
            raise ValueError('Morphology source identity differs')
    rows=[]
    for row in humidity.observations()['cases']:
        t=flakes.T0+row['nominal_temperature_C']
        liquid_excess=(flakes.liquid_saturation_pressure(t)-flakes.ice_saturation_pressure(t))/(flakes.R_VAPOUR*t)*1000
        excess=row['vapor_density_excess_g_m3']
        normalized=excess/liquid_excess
        category='below' if normalized<.2 else 'above' if normalized>.4 else 'within'
        rows.append(dict(datum_id=row['id'],nominal_temperature_C=row['nominal_temperature_C'],
            observed_vapor_density_excess_g_m3=excess,liquid_equilibrium_vapor_density_excess_g_m3=liquid_excess,
            relative_excess_vapor=normalized,onset_range_relative_excess=[.2,.4],
            central_onset_vapor_density_excess_g_m3=.27*liquid_excess,
            position_relative_to_published_onset_range=category))
    return dict(source_record=record,rows=rows,
        evidence='Calculated humidity position relative to a published empirical morphology-onset range; no measured transition time or deposition density is inferred.',
        limits=['The 2024 range describes substantial hollowing/branching across disparate experiments, not a universal sharp transition.',
                'Warm-regime experiments overlap these temperatures but seed shape, size, pressure and history still matter.',
                'The screen supplies no density-amplitude law, arm width or aggregate geometry; it is not installed in a lunar solver.'])


def calculate():
    library=attachment.AttachmentLibrary(baseline.BUILD);identity=library.identity
    sources=audit.verify_sources()
    numerical=author_growth_controls(library)
    empirical=humidity_controls(library)
    rates=measured_shape_rates(library)
    parity=baseline_parity(empirical)
    if attachment.verify_build(baseline.BUILD)!=identity:
        raise ValueError('Source helper changed during controls')
    integrations=numerical['cases']+empirical['cases']
    maximum_refinement=max(value for row in integrations for value in row['refinement']['relative_endpoint_changes'].values() if value is not None)
    paths=[Path(__file__),Path(coupled.__file__),Path(audit.__file__),Path(baseline.__file__),
           Path(attachment.__file__),attachment.DRIVER,Path(humidity.__file__),humidity.DATA,
           Path(flakes.__file__),Path(replay.__file__),Path(lab.__file__),audit.SOURCES,SOURCES]
    return dict(schema='terluna.atmosphere.facet-growth-controls/1',
        producer={str(path.relative_to(ROOT)):sha(path) for path in paths},source_verification=sources,
        attachment_build_record_sha256=sha(baseline.BUILD/'build.json'),
        evidence='Controlled numerical comparisons and restricted Earth validation; no autonomous branching or lunar population forecast.',
        author_numerical_controls=numerical,humidity_controls=empirical,measured_shape_rates=rates,
        morphology_screen=morphology_screen(),semianalytic_baseline_parity=parity,
        summary=dict(refined_trajectories=len(integrations),integration_executions=2*len(integrations),
            complete_trajectories=sum(row['calculation']['complete'] for row in integrations),
            maximum_relative_endpoint_refinement=maximum_refinement,
            maximum_sampled_mass_volume_residual=max(row['calculation']['maximum_sampled_mass_volume_residual'] for row in integrations),
            maximum_rhs_balance_residual=max(row['calculation']['maximum_rhs_balance_residual'] for row in integrations),
            complete_shape_rate_cases=sum(row['complete'] for row in rates['cases']),shape_rate_cases=len(rates['cases'])),
        numerics=dict(scipy_version=scipy.__version__,numpy_version=np.__version__,integrator='RK45 in logsemiaxes plus independently integrated deposition',
            relative_tolerances=[1e-6,1e-8],absolute_tolerance_factor=.01,max_step_s=5.,duration_s=600.))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=RESULT)
    args=p.parse_args();result=calculate()
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result['summary'],indent=2))

if __name__=='__main__':main()
