"""Conditional vapor-growth checks with prescribed laboratory dimensions/speeds.

This separates transport from autonomous habit evolution. It is a new calculation
following Westbrook & Heymsfield (2011), not an exact reproduction of their curves.
No temperature/pressure dependent habit law is calibrated here.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
from scipy.integrate import quad
from atmosphere.snow import flakes, takahashi_growth as lab
from atmosphere.rain import drops
from atmosphere.radiative_convective.thermodynamics import DryAir
from atmosphere.radiative_convective import thermodynamics
from shared.provenance import constants_used

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
RESULT=HERE/'results/crystal_transport_replay.json'
SOURCE_RECORD=HERE/'growth_physics_sources.json'
SOURCE_CACHE=ROOT/'research/runs/snow_scale/growth_physics_sources'
REQUIRED_SOURCES=('westbrook2008.pdf','westbrook_heymsfield2011.pdf','ji_wang1999.pdf','wang_ji1997.pdf')
PRESSURE=101000.  # Takahashi laboratory condition, not a new global constant.
ICE_DENSITY=919.  # Explicit McSnow comparison convention.
HEX_TO_SPHEROID=9*math.sqrt(3)/(4*math.pi)
JW_RANGES={'column':(.2,20.),'plate':(1.,120.),'broad_branch':(1.,120.)}
JW_SOURCE_SCHMIDT=.63  # Fixed value in Ji & Wang 1999, before using the published similarity variable.
TRANSPORT_CACHE=ROOT/'research/runs/snow_scale/mcsnow_ventilation'
TRANSPORT_TEMPERATURES=(-3.7,-5.3,-10.6,-14.4,-20.1)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_sources(source_dir=None,source_record=None,empirical_source=None):
    """Require exact formula and measurement sources; preserve failed retrievals."""
    source_dir=SOURCE_CACHE if source_dir is None else Path(source_dir)
    source_record=SOURCE_RECORD if source_record is None else Path(source_record)
    empirical_source=(lab.original.lunar.CACHE/lab.original.SOURCE_PDF
                      if empirical_source is None else Path(empirical_source))
    entries=json.loads(source_record.read_text())
    retrieved={entry['name']:entry for entry in entries if entry.get('status')=='retrieved'}
    verified=[]
    for name in REQUIRED_SOURCES:
        if name not in retrieved:raise ValueError('No retrieved formula-source identity: '+name)
        item=retrieved[name];path=source_dir/name
        if path.stat().st_size!=item['bytes'] or sha(path)!=item['sha256']:
            raise ValueError('Formula-source identity differs: '+name)
        verified.append(dict(name=name,url=item['url'],bytes=item['bytes'],sha256=item['sha256'],
            path=str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)))
    if sha(empirical_source)!=lab.SOURCE_SHA:raise ValueError('Takahashi measurement-source identity differs')
    return dict(formula_sources=verified,source_register_sha256=sha(source_record),
        empirical_source=dict(path=str(empirical_source.relative_to(ROOT)) if empirical_source.is_relative_to(ROOT) else str(empirical_source),
                              sha256=lab.SOURCE_SHA),
        evidence='Original source bytes verified against recorded identities; no unresolved download is substituted.')


def positive(*values):
    if any(not math.isfinite(v) or v<=0 for v in values):raise ValueError('Expected finite positive quantities')


def volume(full_a,full_c,geometry):
    positive(full_a,full_c)
    if geometry=='spheroid':return math.pi/6*full_a**2*full_c
    if geometry=='hexagonal_prism':return 3*math.sqrt(3)/8*full_a**2*full_c
    raise ValueError('Unidentified geometry')


def hex_capacitance(full_a,full_c):
    """Westbrook et al. 2008 Eq3; full basal span, full polar length, metres."""
    positive(full_a,full_c);phi=full_c/full_a
    if not .01<=phi<=50:raise ValueError('Outside sampled hexagonal-prism aspect range')
    return .58*(1+.95*phi**.75)*full_a/2


def branch_capacitance(full_a,arm_ratio,kind):
    """W08 Eqs7/8 at thickness/span=.01; explicit unmeasured arm-width scenario."""
    positive(full_a,arm_ratio)
    if not .02<=arm_ratio<=.35:raise ValueError('Outside selected published arm-width range')
    if kind=='stellar':coefficient,exponent=.38,4.7
    elif kind=='dendrite':coefficient,exponent=.30,5.8
    else:raise ValueError('Unknown branched crystal')
    return .596*(1-coefficient*math.exp(-exponent*arm_ratio))*full_a/2


def jw_ventilation(reynolds,schmidt,kind):
    """Ji & Wang 1999; Re uses basal width, including columns, never their length.

    No low-Re stitching or Reynolds extrapolation is performed. Reynolds validity
    does not certify another geometry. The source simulations held Sc=0.63; using
    another Sc through their X variable is a separate similarity assumption.
    """
    positive(reynolds,schmidt)
    if kind not in JW_RANGES:raise ValueError('Unknown ventilation habit')
    lower,upper=JW_RANGES[kind]
    if not lower<=reynolds<=upper:raise ValueError('Outside published Reynolds interval')
    x=schmidt**(1/3)*math.sqrt(reynolds)
    if kind=='column':
        y=x/4
        return 1-.00668*y+2.39402*y*y+.73409*y**3-.73911*y**4
    y=x/10
    if kind=='plate':return 1-.60420*y+2.79820*y*y+.31933*y**3-.06247*y**4
    return 1+.35463*y+3.55338*y*y


def air_state(temperature_C,pressure=PRESSURE):
    temperature=flakes.T0+temperature_C;positive(temperature,pressure)
    air=DryAir(pressure);e=flakes.liquid_saturation_pressure(temperature)
    if e>=pressure:raise ValueError('Vapor pressure exceeds total pressure')
    eta,dv,ka,_=drops.transport(temperature,pressure,e/pressure,air.fractions)
    rho=(pressure-e)*air.molar_mass/(flakes.GAS_CONSTANT*temperature)+e/(flakes.R_VAPOUR*temperature)
    return dict(temperature_K=temperature,pressure_Pa=pressure,vapor_pressure_Pa=e,
        ice_saturation_pressure_Pa=flakes.ice_saturation_pressure(temperature),
        density_kg_m3=rho,viscosity_Pa_s=eta,diffusivity_m2_s=dv,conductivity_W_m_K=ka)


def mass_rate(capacitance,ventilation,state):
    """Maxwell vapor/latent-heat resistances, with common enhancement f as WH2011."""
    positive(capacitance,ventilation)
    t=state['temperature_K'];es=state['ice_saturation_pressure_Pa']
    latent=drops.latent_heat(flakes.T0)+flakes.LATENT_FUSION
    resistance=(latent/(flakes.R_VAPOUR*t)-1)*latent/(state['conductivity_W_m_K']*t)
    resistance+=flakes.R_VAPOUR*t/(state['diffusivity_m2_s']*es)
    return 4*math.pi*capacitance*ventilation*(state['vapor_pressure_Pa']/es-1)/resistance


def branch_habit(temperature_C,age):
    if temperature_C in (-14.4,-16.5):return 'broad_branch'
    if temperature_C==-12.2:return 'plate' if age<10 else 'broad_branch'
    if temperature_C==-18.2:return 'plate'
    if temperature_C==-5.3:return 'column'
    return None


def capacitance_at(temperature_C,age,variant):
    obs=lab.observed(temperature_C,age);a=obs['equatorial_dimension_m'];c=obs['polar_dimension_m']
    if variant=='hexagonal_prism' or (temperature_C==-12.2 and age<10):return hex_capacitance(a,c)
    kind,ratio=variant.split(':')
    return branch_capacitance(a,float(ratio),kind)


def diagnostics(temperature_C,age,variant,ventilation,state,*,jw_schmidt=None):
    obs=lab.observed(temperature_C,age);a=obs['equatorial_dimension_m'];c=obs['polar_dimension_m']
    speed=lab.fall_speed(temperature_C,age);nu=state['viscosity_Pa_s']/state['density_kg_m3']
    sc=nu/state['diffusivity_m2_s'];reb=speed*a/nu;red=speed*max(a,c)/nu
    if jw_schmidt is not None and ventilation!='JW':raise ValueError('Schmidt sensitivity applies only to JW')
    vent_sc=sc if jw_schmidt is None else jw_schmidt
    if ventilation=='diffusion':factor=1.
    elif ventilation=='HPD':
        x=sc**(1/3)*math.sqrt(red);factor=1+.14*x*x if x<1 else .86+.28*x
    elif ventilation=='JW':factor=jw_ventilation(reb,vent_sc,branch_habit(temperature_C,age))
    else:raise ValueError('Unknown transport control')
    cap=capacitance_at(temperature_C,age,variant)
    return dict(age_min=age,capacitance_m=cap,ventilation=factor,reynolds_basal=reb,
        reynolds_maximum_dimension=red,schmidt=sc,ventilation_schmidt=vent_sc,
        jw_source_schmidt=JW_SOURCE_SCHMIDT if ventilation=='JW' else None,observed_fit_speed_m_s=speed,
        predicted_mass_rate_kg_s=mass_rate(cap,factor,state))


def geometry_audit():
    rows=[]
    for temperature,(end,habit,*_) in lab.FITS.items():
        for age in (3.,10.,float(end)):
            obs=lab.observed(temperature,age);a=obs['equatorial_dimension_m'];c=obs['polar_dimension_m'];m=obs['mass_kg']
            rows.append(dict(temperature_C=temperature,age_min=age,habit=habit,
                **{k:obs[k] for k in ('mass_kg','equatorial_dimension_m','polar_dimension_m')},
                densities_kg_m3={g:m/volume(a,c,g) for g in ('spheroid','hexagonal_prism')},
                admissible_at_native_density={g:m/volume(a,c,g)<=ICE_DENSITY for g in ('spheroid','hexagonal_prism')}))
    return dict(hexagonal_to_spheroid_volume_ratio=HEX_TO_SPHEROID,rows=rows,
        interpretation='Takahashi density uses a circumscribing hexagonal prism. Spheroid rejection is a mapping failure, not proof of impossible measured ice.',
        limits=['Independent fits can still disagree with bulk density under the original geometry.',
                'Thin planar thickness was inferred from mass/basal area; it is not an independent third observation.',
                'Do not replace spheroid volume while retaining its native axis reconstruction; this would alter the state.'])


def fitted_exponents(temperature_C,age):
    end,habit,*curves=lab.FITS[temperature_C];lab.observed(temperature_C,age)
    for curve in curves:
        if any(age==join for join,*_ in curve[:-1]):raise ValueError('Fit derivative undefined at a piecewise join')
    result=[]
    for curve in curves:
        selected=next((seg for seg in curve if age<seg[0]),curve[-1]);result.append(selected[2])
    return result


def derivative_audit():
    rows=[]
    for temperature in lab.FITS:
        obs=lab.observed(temperature,10.);nm,na,nc=fitted_exponents(temperature,10.)
        rho=obs['mass_kg']/volume(obs['equatorial_dimension_m'],obs['polar_dimension_m'],'spheroid')
        rows.append(dict(temperature_C=temperature,age_min=10.,mass_exponent=nm,equatorial_exponent=na,polar_exponent=nc,
            implied_constant_habit_growth_ratio=nc/na,implied_spheroid_deposition_density_kg_m3=rho*nm/(2*na+nc),
            mass_rate_kg_s=nm*obs['mass_kg']/600,
            equatorial_rate_m_s=na*obs['equatorial_dimension_m']/600,
            polar_rate_m_s=nc*obs['polar_dimension_m']/600))
    return dict(evidence='Derivatives of separate empirical power-law fits, not new directly measured growth laws.',rows=rows)


def replay(temperature,variant,ventilation,epsrel,*,jw_schmidt=None,state_override=None):
    state=air_state(temperature) if state_override is None else dict(state_override)
    end=float(lab.FITS[temperature][0]);initial=lab.observed(temperature,3.)['mass_kg']
    joins=sorted({3.,end,*[float(v[0]) for curve in lab.FITS[temperature][2:] for v in curve if 3<v[0]<end],
        *[float(v[0]) for v in lab.FALL_SPEED[temperature] if 3<v[0]<end],*([10.] if 3<10<end else [])})
    # Piecewise quadrature never differentiates independently fitted curves or crosses a join implicitly.
    error=None;mass=initial;quad_error=0.;records=[]
    kwargs={} if jw_schmidt is None else dict(jw_schmidt=jw_schmidt)
    try:
        records=[diagnostics(temperature,age,variant,ventilation,state,**kwargs) for age in sorted(set([3.,5.,10.,end])) if age<=end]
        for lo,hi in zip(joins,joins[1:]):
            amount,err=quad(lambda age:60*diagnostics(temperature,age,variant,ventilation,state,**kwargs)['predicted_mass_rate_kg_s'],
                            lo,hi,epsabs=1e-18,epsrel=epsrel,limit=100)
            mass+=amount;quad_error+=err
    except ValueError as exc:error=str(exc)
    return dict(temperature_C=temperature,capacitance_variant=variant,ventilation=ventilation,
        initial_age_min=3.,end_age_min=end,air=state,complete=error is None,domain_error=error,
        jw_schmidt_override=jw_schmidt,
        predicted_end_mass_kg=mass if error is None else None,
        observed_end_fit_mass_kg=lab.observed(temperature,end)['mass_kg'],
        mass_ratio_to_fit=mass/lab.observed(temperature,end)['mass_kg'] if error is None else None,
        quadrature_absolute_error_kg=quad_error if error is None else None,diagnostics=records,
        geometry_condition='Measured axes and fall speed prescribed; no independent shape or drag prediction.',
        branching_condition='Branched capacitance fixes thickness/span=0.01 and assumes an unmeasured arm-width ratio.' if ':' in variant else 'Solid enclosing hexagonal prism, not resolved dendritic geometry.')


def native_transport_inputs(cache=None):
    """Read verified native coefficients; this never executes or modifies McSnow."""
    from atmosphere.snow import mcsnow_ventilation_controls as control
    cache=TRANSPORT_CACHE if cache is None else Path(cache)
    build=control.verify_build(cache/'build');producer=control.producer()
    if build['producer']!=producer:raise ValueError('Native transport build producer differs')
    identities=[]
    for temperature in TRANSPORT_TEMPERATURES:
        name=f'm0_d0_igf2_sphere_T{abs(temperature):g}_dt0.05'.replace('.','p')
        path=cache/'runs'/name;record=json.loads((path/'run.json').read_text())
        if record['producer']!=producer or record['build_record_sha256']!=sha(cache/'build/build.json'):
            raise ValueError('Native transport execution identity differs')
        config=record['configuration']
        if not record['complete'] or config['mean_mode']!=0 or config['directional_mode']!=0 or record['temperature_C']!=temperature:
            raise ValueError('Native transport source case differs')
        for group in ('inputs','outputs'):
            for filename,digest in record[group].items():
                if sha(path/filename)!=digest:raise ValueError('Native transport source bytes differ: '+filename)
        if (path/'crystal.nml').read_text()!=control.native.namelist(config):raise ValueError('Native transport input differs')
        air,air_errors=control.native.read_csv(path/'air.csv',('temperature_K','pressure_Pa','vapor_pressure_Pa','density_kg_m3','dynamic_viscosity_Pa_s'))
        traces,trace_errors=control.native.read_csv(path/'growth_trace.csv',('diffusivity_raw_m2_s','conductivity_raw_W_m_K'))
        if air_errors or trace_errors or len(air)!=1 or not traces or air[0]!=record['air'] or traces!=record['growth_trace']:
            raise ValueError('Native transport diagnostic differs')
        state=air_state(temperature)
        for key in ('temperature_K','pressure_Pa','vapor_pressure_Pa'):
            if not math.isclose(air[0][key],state[key],rel_tol=1e-13,abs_tol=1e-10):raise ValueError('Native transport air forcing differs')
        trace=traces[0]
        for row in traces:
            for key in ('diffusivity_raw_m2_s','conductivity_raw_W_m_K'):
                if row[key]!=trace[key]:raise ValueError('Native transport coefficient changed in fixed air')
        identities.append(dict(temperature_C=temperature,
            run=str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path),run_sha256=sha(path/'run.json'),
            build_record_sha256=record['build_record_sha256'],producer=producer,inputs=record['inputs'],outputs=record['outputs'],
            native_diffusivity_m2_s=trace['diffusivity_raw_m2_s'],native_conductivity_W_m_K=trace['conductivity_raw_W_m_K'],
            native_density_kg_m3=air[0]['density_kg_m3'],native_viscosity_Pa_s=air[0]['dynamic_viscosity_Pa_s'],
            mixture_air=state,mixture_over_native_diffusivity=state['diffusivity_m2_s']/trace['diffusivity_raw_m2_s'],
            mixture_over_native_conductivity=state['conductivity_W_m_K']/trace['conductivity_raw_W_m_K']))
    return identities


def air_transport_sensitivity(cases):
    sources=native_transport_inputs();comparisons=[]
    for source in sources:
        temperature=source['temperature_C'];changed=dict(source['mixture_air'])
        changed['diffusivity_m2_s']=source['native_diffusivity_m2_s']
        changed['conductivity_W_m_K']=source['native_conductivity_W_m_K']
        for baseline in cases:
            if baseline['temperature_C']!=temperature:continue
            variant=baseline['capacitance_variant'];method=baseline['ventilation']
            coarse=replay(temperature,variant,method,1e-5,state_override=changed)
            fine=replay(temperature,variant,method,1e-9,state_override=changed)
            both=baseline['complete'] and fine['complete']
            initial=lab.observed(temperature,3.)['mass_kg']
            comparisons.append(dict(temperature_C=temperature,capacitance_variant=variant,ventilation=method,
                both_complete=both,mixture_mass_ratio_to_fit=baseline['mass_ratio_to_fit'],native_coefficients=fine,
                native_over_mixture_endpoint_mass=fine['predicted_end_mass_kg']/baseline['predicted_end_mass_kg'] if both else None,
                native_over_mixture_deposited_mass=(fine['predicted_end_mass_kg']-initial)/(baseline['predicted_end_mass_kg']-initial) if both else None,
                native_coefficients_quadrature_relative_change=abs(coarse['predicted_end_mass_kg']/fine['predicted_end_mass_kg']-1) if coarse['complete'] and fine['complete'] else None))
    return dict(evidence='Transport-coefficients-only control; coefficients from verified unchanged native McSnow source traces.',
        sources=sources,cases=comparisons,
        changed=['diffusivity_m2_s','conductivity_W_m_K'],
        preserved=['Murphy-Koop vapor and ice-saturation pressures','mixture density and viscosity','latent heat',
                   'empirical axes and fall speed','capacitance geometry','ventilation formula'],
        scope='Only raw diffusivity and conductivity are replaced. Schmidt number changes consistently with diffusivity; Reynolds number and density/viscosity remain unchanged. Native density/viscosity are recorded for comparison only. No source convention is selected as a repair.')


def calculate():
    sources=verify_sources()
    cases=[];refinements=[];schmidt_sensitivities=[]
    for temperature in lab.FITS:
        variants=['hexagonal_prism']
        if temperature in (-12.2,-14.4,-16.5):variants += [f'{kind}:{ratio}' for kind in ('stellar','dendrite') for ratio in (.02,.35)]
        methods=['diffusion','HPD']+(['JW'] if branch_habit(temperature,3.) else [])
        for variant in variants:
            for method in methods:
                coarse=replay(temperature,variant,method,1e-5);fine=replay(temperature,variant,method,1e-9)
                cases.append(fine)
                refinements.append(dict(temperature_C=temperature,capacitance_variant=variant,ventilation=method,
                    both_complete=coarse['complete'] and fine['complete'],
                    relative_endpoint_mass_change=abs(coarse['predicted_end_mass_kg']/fine['predicted_end_mass_kg']-1) if coarse['complete'] and fine['complete'] else None))
                if method=='JW':
                    fixed_coarse=replay(temperature,variant,method,1e-5,jw_schmidt=JW_SOURCE_SCHMIDT)
                    fixed_fine=replay(temperature,variant,method,1e-9,jw_schmidt=JW_SOURCE_SCHMIDT)
                    both=fine['complete'] and fixed_fine['complete']
                    schmidt_sensitivities.append(dict(temperature_C=temperature,capacitance_variant=variant,
                        source_schmidt=JW_SOURCE_SCHMIDT,calculated_mixture_schmidt=fine['air']['viscosity_Pa_s']/fine['air']['density_kg_m3']/fine['air']['diffusivity_m2_s'],
                        both_complete=both,calculated_schmidt_mass_ratio_to_fit=fine['mass_ratio_to_fit'],fixed_source_schmidt=fixed_fine,
                        fixed_over_calculated_endpoint_mass=fixed_fine['predicted_end_mass_kg']/fine['predicted_end_mass_kg'] if both else None,
                        fixed_schmidt_quadrature_relative_change=abs(fixed_coarse['predicted_end_mass_kg']/fixed_fine['predicted_end_mass_kg']-1) if fixed_coarse['complete'] and fixed_fine['complete'] else None,
                        scope='Only the Schmidt number inside the JW ventilation polynomial changes; capacitance, prescribed axes/speed and physical air resistances are identical.'))
    paths=[Path(__file__),Path(lab.__file__),Path(flakes.__file__),Path(drops.__file__),
           Path(thermodynamics.__file__),SOURCE_RECORD]
    return dict(schema='terluna.atmosphere.crystal-transport-replay/1',
        evidence='Conditional transport calculation against laboratory fits; no autonomous habit validation or lunar snowfall forecast.',
        producer={str(p.relative_to(ROOT)):sha(p) for p in paths},source_verification=sources,
        imported_constants=constants_used(paths),geometry=geometry_audit(),fit_derivatives=derivative_audit(),
        cases=cases,quadrature_refinements=refinements,jw_source_schmidt_sensitivity=schmidt_sensitivities,
        air_transport_sensitivity=air_transport_sensitivity(cases),
        protocol=['Fixed laboratory temperature and pressure; vapor pressure at Murphy-Koop liquid saturation.',
            'In-house mixture transport and latent heat; same f multiplies vapor and heat transport as published comparison.',
            'Observed fitted full axes and wind-tunnel fall speed imposed at every quadrature point.',
            'Each mass starts at the fitted three-minute mass; no assumed nucleation history.',
            'HPD uses maximum dimension; JW uses basal width even for columns, with explicit Reynolds limits.',
            'For -12.2 C, the branched variants switch from hexagonal-prism capacitance to the stated branch scenario at10min.',
            'Quadrature split at all mass/axis/fall-speed fit joins; tolerances1e-5 and1e-9 compared.'],
        limits=['Not exact reproduction of Westbrook-Heymsfield2011 numerical curves or their particular photographed branch geometry.',
            'Branched arm widths are unreported; scenarios are not measured confidence intervals.',
            'JW Reynolds validity is necessary but does not establish arbitrary-geometry validity; idealized source shapes remain an assumption.',
            'No JW Reynolds extrapolation or low-Re splice; invalid cases retain a reason and no endpoint.',
            'JW simulations used Schmidt number0.63; evaluating their published X polynomial at the calculated mixture Schmidt number assumes similarity outside that fixed-source condition.',
            'No facet attachment resistance, prognostic shape, finite reservoir, collision, or gravity transfer is validated here.'])


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=RESULT);args=p.parse_args()
    product=calculate();args.output.write_text(json.dumps(product,indent=2,allow_nan=False)+'\n')
    print(len(product['cases']),'conditional cases;',sum(c['complete'] for c in product['cases']),'complete')
if __name__=='__main__':main()
