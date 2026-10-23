"""Transfer the explicit plate-ventilation control to prescribed lunar conditions.

Native McSnow growth and drag are reused. The supply-limited descent is a steady,
monodisperse, collision-free column with an imposed cooling/moisture budget.
"""
from __future__ import annotations

import argparse
import difflib
import itertools
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import time

from shared.constants import MOON_SURFACE_GRAVITY
from atmosphere.snow import flakes
from atmosphere.snow import mcsnow_growth_trace as trace

model, native = trace.model, trace.native
ROOT = model.ROOT
CACHE = ROOT/'research/runs/snow_scale/mcsnow_lunar_growth'
DRIVER = ROOT/'atmosphere/snow/mcsnow_lunar_driver.f90'
PRODUCT = ROOT/'atmosphere/snow/results/mcsnow_lunar_growth.json'
MODES = ('native','plate_base')
STATUSES = {0:'duration_reached',1:'reached_freezing_level',2:'diameter_censored'}


def source_gravity():
    body = (model.BASELINE/'source/src/mo_atmo_types.f90').read_text()
    return float(re.search(r'grav\s*=\s*([\d.]+)_wp',body)[1])


def forcing():
    column,binding = flakes.moon_column()
    # Retain the exact previously executed lunar population forcing.
    path = ROOT/'atmosphere/snow/results/scale_lunar_comparison.json'
    product = json.loads(path.read_text())
    case = product['cases']['lunar_prescribed_p4']
    if binding['inputs']!=case['settings']['producer']['inputs']:
        raise ValueError('Lunar climate input differs from the identified prior experiment')
    for name,digest in binding['inputs'].items():
        if native.sha(ROOT/name)!=digest:
            raise ValueError('Lunar input hash mismatch')
    for row in case['settings']['profile']:
        if not math.isclose(column.state(column.height_of(row['t']))['p'],row['p'],rel_tol=1e-13):
            raise ValueError('Lunar pressure profile differs from the prior experiment')
    return column,binding,case,path


def patches(mode,world):
    if mode not in MODES or world not in ('earth','moon'):
        raise ValueError('Unknown source candidate/world')
    original = {name:(model.BASELINE/'source/src'/name).read_text()
                for name in ('mo_depo.f90','mo_atmo_types.f90')}
    dep = trace.source_patch(original['mo_depo.f90'],mode)
    dep = model.replace_once(dep,'  REAL(wp), PUBLIC :: growth_trace(22)=0._wp',
        '  REAL(wp), PUBLIC :: growth_trace(22)=0._wp\n'+
        '  REAL(wp), PUBLIC :: supply_rate=-1._wp, requested_rate=0._wp')
    dep = model.replace_once(dep,'  growth_trace(21)=dMdt',
        '  requested_rate=dMdt\n'+
        '  IF (supply_rate>=0._wp.AND.dMdt>0._wp) dMdt=MIN(dMdt,supply_rate)\n'+
        '  growth_trace(21)=dMdt')
    atmo = original['mo_atmo_types.f90']
    if world=='moon':
        atmo = model.replace_once(atmo,'    grav  = 9.81_wp,',f'    grav  = {MOON_SURFACE_GRAVITY:.17g}_wp,')
    return original,dict(zip(original,(dep,atmo)))


def build(mode,world):
    baseline = native.verify_build(model.BASELINE)
    path = CACHE/f'build_{mode}_{world}'
    if path.exists():
        return verify_build(path)
    original,changed = patches(mode,world)
    path.mkdir(parents=True)
    shutil.copytree(model.BASELINE/'source',path/'source')
    diff = ''
    for name,body in changed.items():
        (path/'source/src'/name).write_text(body)
        diff += ''.join(difflib.unified_diff(original[name].splitlines(True),body.splitlines(True),
                                            fromfile='native/'+name,tofile=mode+'_'+world+'/'+name))
    (path/'lunar.patch').write_text(diff)
    shutil.copyfile(DRIVER,path/DRIVER.name)
    for name in native.TABLE_NAMES:
        shutil.copyfile(model.BASELINE/name,path/name)
    compiler = baseline['compiler_path']
    sources = native.selected_sources(path/'source/src',path/DRIVER.name,compiler)
    objects,commands = [],[]
    with (path/'build.log').open('w') as log:
        for i,src in enumerate(sources):
            obj = path/f'{i:02d}_{src.stem}.o';objects.append(str(obj))
            cmd = [compiler,*native.FLAGS,'-I'+str(path/'source/src'),'-I'+str(path),'-J'+str(path),'-c',str(src),'-o',str(obj)]
            commands.append(cmd)
            subprocess.run(cmd,cwd=path,env=native.ENV,stdout=log,stderr=subprocess.STDOUT,check=True)
        cmd = [compiler,*native.FLAGS,*objects,'-o',str(path/'crystal')]
        commands.append(cmd)
        subprocess.run(cmd,cwd=path,env=native.ENV,stdout=log,stderr=subprocess.STDOUT,check=True)
    record = dict(schema='terluna.atmosphere.mcsnow-lunar-build/1',mode=mode,world=world,
        gravity_m_s2=MOON_SURFACE_GRAVITY if world=='moon' else source_gravity(),
        baseline_record_sha256=native.sha(model.BASELINE/'build.json'),source_identity=baseline['source_identity'],
        table_identity=baseline['table_identity'],
        source_files={str(p.relative_to(path/'source')):native.sha(p) for p in sorted((path/'source').rglob('*')) if p.is_file()},
        files_sha256={name:native.sha(path/name) for name in (DRIVER.name,'crystal','lunar.patch','build.log')},
        compiler_path=compiler,compiler_sha256=native.sha(compiler),compiler_version=baseline['compiler_version'],
        flags=native.FLAGS,commands=commands,
        evidence='Explicit mean-ventilation candidate, gravity substitution and optional externally supplied growth ceiling.')
    model.save(path/'build.json',record)
    return record


def verify_build(path):
    path = Path(path)
    record = json.loads((path/'build.json').read_text())
    native.verify_build(model.BASELINE)
    if native.sha(model.BASELINE/'build.json') != record['baseline_record_sha256']:
        raise ValueError('Baseline changed')
    for name,digest in record['files_sha256'].items():
        if native.sha(path/name) != digest:
            raise ValueError('Build artifact changed: '+name)
    for name,digest in record['source_files'].items():
        if native.sha(path/'source'/name) != digest:
            raise ValueError('Source changed: '+name)
    for name,item in record['table_identity'].items():
        native.pinned.verify((path/name).read_bytes(),item)
    if native.sha(DRIVER) != record['files_sha256'][DRIVER.name]:
        raise ValueError('Lunar caller changed')
    _,expected = patches(record['mode'],record['world'])
    if any((path/'source/src'/n).read_text()!=b for n,b in expected.items()):
        raise ValueError('Build does not implement its stated intervention')
    return record


def validate_output(rows,summary,air,steps,config,gravity):
    if not rows or len(summary)!=1 or len(air)!=1:
        raise ValueError('Missing trajectory, summary or air')
    final,initial = rows[-1],rows[0]
    status = summary[0]
    if status['status_code'] not in STATUSES:
        raise ValueError('Unknown termination status')
    if initial['elapsed_s']!=0 or final['elapsed_s']!=status['elapsed_s']:
        raise ValueError('Incorrect time coverage')
    if not math.isclose(air[0]['native_gravity_m_s2'],gravity,rel_tol=1e-14):
        raise ValueError('Gravity substitution differs')
    if any(b['elapsed_s']<=a['elapsed_s'] for a,b in zip(rows,rows[1:])):
        raise ValueError('Nonmonotone saved times')
    if status['status_code']==0 and final['elapsed_s']!=config['duration_s']:
        raise ValueError('Unfinished requested duration')
    if status['status_code']==1 and final['height_above_freezing_m']>1e-9:
        raise ValueError('Premature column exit')
    if status['status_code']==2 and final['diameter_m']<config.get('diameter_stop_m',1.):
        raise ValueError('Premature size censor')
    if config.get('max_mass_fraction',0)>0 and status['maximum_step_mass_fraction']>config['max_mass_fraction']*(1+1e-12):
        raise ValueError('Step control violated')
    for row in rows:
        if min(row[k] for k in (*model.QUANTITIES,'density_kg_m3'))<=0:
            raise ValueError('Invalid physical state')
        volume = math.pi/6*row['equatorial_dimension_m']**2*row['polar_dimension_m']
        if not math.isclose(row['ice_volume_m3'],volume,rel_tol=5e-13):
            raise ValueError('Axis/volume mapping differs')
        residual = row['mass_kg']-initial['mass_kg']-row['integrated_deposition_kg']
        if abs(residual)>1e-10*row['mass_kg'] or not math.isclose(
                residual,row['ice_mass_balance_residual_kg'],abs_tol=1e-14*row['mass_kg']):
            raise ValueError('Deposition account fails')
        if config.get('water_supply_kg_m2_s',-1)>=0:
            expected = supply_budget(config['water_supply_kg_m2_s'],config['number_flux_m2_s'],
                                     config['cloud_depth_m'],config['cloud_depth_m']-row['height_above_freezing_m'])
            if not math.isclose(row['allocated_water_kg'],expected,rel_tol=1e-9,abs_tol=1e-22):
                raise ValueError('Supply allocation differs')
            if row['integrated_deposition_kg']>expected+1e-9*max(expected,initial['mass_kg']):
                raise ValueError('Growth exceeds the integrated water allowance')
    for step in steps:
        if not math.isclose(step['state_elapsed_s']-step['evaluation_elapsed_s'],step['step_s'],rel_tol=1e-7,abs_tol=1e-10):
            raise ValueError('Incorrect source-evaluation time')
        if config.get('water_supply_kg_m2_s',-1)>=0:
            rate = config['water_supply_kg_m2_s']/config['cloud_depth_m']*step['fall_speed_m_s']/config['number_flux_m2_s']
            if not math.isclose(step['supply_rate_kg_s'],rate,rel_tol=1e-13):
                raise ValueError('Per-crystal supply rate differs from the number-flux account')
        expected = step['requested_rate_kg_s']
        if expected>0 and step['supply_rate_kg_s']>=0:
            expected = min(expected,step['supply_rate_kg_s'])
        # The caller diagnoses admitted rate by subtracting two rounded masses.
        # Near freezing, the clipped last step can add less than 1e-10 of the
        # existing mass; compare against its floating-point resolution as well
        # as the source rate. Output steps have exact post-state snapshots.
        post = next((r for r in rows if r['elapsed_s']==step['state_elapsed_s']),None)
        if post is not None:
            mass_scale = post['mass_kg']
        elif step['evaluation_elapsed_s']==0:
            mass_scale = initial['mass_kg']+abs(expected)*step['step_s']
        else:
            raise ValueError('Missing post-state for a source evaluation')
        roundoff = 2*math.ulp(mass_scale)/step['step_s']
        if not math.isclose(step['admitted_rate_kg_s'],expected,rel_tol=1e-12,abs_tol=roundoff):
            raise ValueError('Admitted source rate differs from the declared gate')


def analyze(path,config,returncode,gravity):
    rows,e1 = native.read_csv(path/'trajectory.csv',('elapsed_s',*model.QUANTITIES))
    summary,e2 = native.read_csv(path/'summary.csv',('elapsed_s','status_code','integrated_deposition_kg'))
    air,e3 = native.read_csv(path/'air.csv',('native_gravity_m_s2',))
    steps,e4 = native.read_csv(path/'step_trace.csv',('step_s','requested_rate_kg_s','admitted_rate_kg_s'))
    errors = e1+e2+e3+e4
    if returncode==0 and not errors:
        try:
            validate_output(rows,summary,air,steps,config,gravity)
        except ValueError as exc:
            errors.append(str(exc))
    return dict(snapshots=rows,summary=summary[0] if summary else None,air=air,steps=steps,
                validation_errors=errors,complete=returncode==0 and not errors)


def run(entry):
    mode,world,config = (entry[k] for k in ('mode','world','configuration'))
    path = CACHE/'runs'/entry['name']
    buildpath = CACHE/f'build_{mode}_{world}'
    identity = verify_build(buildpath)
    if path.exists():
        record = verify_run(path)
        if record['configuration']!=config or record['mode']!=mode or record['world']!=world:
            raise ValueError('Existing run has different inputs')
        return record
    path.mkdir(parents=True)
    (path/'crystal.nml').write_text(native.namelist(config))
    for name in native.TABLE_NAMES:
        shutil.copyfile(buildpath/name,path/name)
    start = time.monotonic()
    with (path/'run.log').open('w') as log:
        proc = subprocess.run([str(buildpath/'crystal'),'crystal.nml'],cwd=path,env=native.ENV,stdout=log,stderr=subprocess.STDOUT)
    record = dict(schema='terluna.atmosphere.mcsnow-lunar-run/1',mode=mode,world=world,configuration=config,
        build_record_sha256=native.sha(buildpath/'build.json'),returncode=proc.returncode,elapsed_wall_s=time.monotonic()-start,
        **analyze(path,config,proc.returncode,identity['gravity_m_s2']),
        file_sha256={p.name:native.sha(p) for p in sorted(path.iterdir()) if p.is_file()})
    model.save(path/'run.json',record)
    return record


def verify_run(path):
    path = Path(path)
    record = json.loads((path/'run.json').read_text())
    buildpath = CACHE/f"build_{record['mode']}_{record['world']}"
    identity = verify_build(buildpath)
    if native.sha(buildpath/'build.json')!=record['build_record_sha256']:
        raise ValueError('Run/build mismatch')
    for name,digest in record['file_sha256'].items():
        if native.sha(path/name)!=digest:
            raise ValueError('Raw file changed: '+name)
    if (path/'crystal.nml').read_text()!=native.namelist(record['configuration']):
        raise ValueError('Namelist changed')
    replay = analyze(path,record['configuration'],record['returncode'],identity['gravity_m_s2'])
    if any(record[k]!=v for k,v in replay.items()):
        raise ValueError('Saved analysis differs from raw output')
    return record


def supply_budget(water_flux,number_flux,depth,fallen):
    """Steady equal-number-flux column: allocated kg per crystal along a path."""
    if not all(math.isfinite(x) for x in (water_flux,number_flux,depth,fallen)):
        raise ValueError('Nonfinite supply input')
    if min(number_flux,depth)<=0 or water_flux<0 or not 0<=fallen<=depth:
        raise ValueError('Invalid supply input')
    return water_flux/number_flux*fallen/depth


def cases():
    column,binding,prior,path = forcing()
    depth = float(column.height_of(258.15)-column.z_freeze)
    ratio = MOON_SURFACE_GRAVITY/source_gravity()
    for Tc,world,pressure_kind,mode,dt in itertools.product((-14.4,-5.3),('earth','moon'),('lab','lunar'),MODES,(.03125,.0078125)):
        pressure = 101000. if pressure_kind=='lab' else float(column.state(column.height_of(Tc+273.15))['p'])
        config = model.configuration(Tc+273.15,dt=dt,duration=1800,pressure=pressure)
        yield dict(name=f'fixed_{Tc}_{world}_{pressure_kind}_{mode}_dt{dt}',group='fixed',mode=mode,world=world,
                   pressure_kind=pressure_kind,temperature_C=Tc,configuration=config)
    for mode in MODES:
        yield dict(name='parity_'+mode,group='parity',mode=mode,world='earth',
                   configuration=model.configuration(258.65,dt=.015625))
    # Exact source and allowance of the existing SCALE lunar comparison.
    settings = prior['settings']
    source_flux = settings['original_settings']['injected_number_flux_m2_s']
    water = settings['original_settings']['snowfall_mm_h']/3600.
    duration = settings['duration_s']
    descriptions = []
    for world,mode,supplied in itertools.product(('earth','moon'),MODES,(False,True)):
        descriptions.append((world,'shallow' if world=='earth' else 'deep',mode,supplied,4e-6,'main'))
    for world,depth_kind in (('earth','deep'),('moon','shallow')):
        descriptions.append((world,depth_kind,'plate_base',True,4e-6,'crossed'))
    for world,supplied in itertools.product(('earth','moon'),(False,True)):
        descriptions.append((world,'shallow' if world=='earth' else 'deep','plate_base',supplied,2e-4,'seed'))
    for world,depth_kind,mode,supplied,diameter,group in descriptions:
        factor = ratio if depth_kind=='shallow' else 1.
        for dt in (2.,.5):
            config = model.configuration(258.15,dt=dt,duration=1800.)
            config.update(column_mode=1,cloud_depth_m=depth*factor,freezing_pressure_Pa=column.p_freeze,
                pressure_exponent=float(column.exponent),lapse_K_m=float(column.lapse/factor),rh_water=1.,
                water_supply_kg_m2_s=water if supplied else -1.,number_flux_m2_s=source_flux,
                duration_s=duration,output_interval_s=600.,max_mass_fraction=.01 if dt==2 else .005,
                initial_equatorial_m=diameter,initial_polar_m=diameter,diameter_stop_m=1.)
            yield dict(name=f'fall_{group}_{world}_{depth_kind}_{mode}_{int(supplied)}_dt{dt}',group='fall',
                       control=group,world=world,depth_kind=depth_kind,mode=mode,supply_limited=supplied,configuration=config)
    for world in ('earth','moon'):
        config = model.configuration(273.15,dt=.125,duration=600.)
        config.update(column_mode=1,cloud_depth_m=1.,freezing_pressure_Pa=101000.,pressure_exponent=0.,
                      lapse_K_m=0.,growth_enabled=False,initial_equatorial_m=2e-4,initial_polar_m=2e-4)
        yield dict(name='drift_'+world,group='drift',mode='plate_base',world=world,configuration=config)
    for source_factor,dt in itertools.product((.1,.01),(2.,.5)):
        config = model.configuration(258.15,dt=dt,duration=1800.)
        config.update(column_mode=1,cloud_depth_m=depth,freezing_pressure_Pa=column.p_freeze,
            pressure_exponent=float(column.exponent),lapse_K_m=float(column.lapse),rh_water=1.,
            water_supply_kg_m2_s=water,number_flux_m2_s=source_flux*source_factor,
            duration_s=duration,output_interval_s=600.,max_mass_fraction=.01 if dt==2 else .005,
            diameter_stop_m=1.)
        yield dict(name=f'flux_{source_factor}_dt{dt}',group='flux',world='moon',mode='plate_base',
                   source_factor=source_factor,configuration=config)


def assemble():
    column,binding,prior,source_path = forcing()
    entries = list(cases())
    records,manifest,fixed,paths,parity,drift,refinements = {},[],[],[],[],[],[]
    quantities = (*model.QUANTITIES,'density_kg_m3','projected_area_m2')
    for entry in entries:
        path = CACHE/'runs'/entry['name']
        record = verify_run(path)
        if record['configuration']!=entry['configuration']:
            raise ValueError('Matrix configuration changed')
        records[entry['name']] = record
        metadata = {k:v for k,v in entry.items() if k!='configuration'}
        manifest.append(dict(**entry,path=str(path.relative_to(ROOT)),run_record_sha256=native.sha(path/'run.json'),
            complete=record['complete'],returncode=record['returncode'],validation_errors=record['validation_errors'],
            termination=STATUSES.get(record['summary']['status_code']) if record['summary'] else None))
        if not record['complete']:
            continue
        if entry['group']=='fixed':
            points = []
            for age in (180,600,1800):
                row = next((s for s in record['snapshots'] if s['elapsed_s']==age),None)
                points.append(dict(age_s=age,state=row))
            fixed.append(dict(**metadata,dt_s=entry['configuration']['dt_s'],points=points,
                              summary=record['summary'],last_state=record['snapshots'][-1]))
        elif entry['group'] in ('fall','flux'):
            final = record['snapshots'][-1]
            conf = entry['configuration']
            flux = conf['number_flux_m2_s']
            points = []
            for row in record['snapshots']:
                points.append(row)
            paths.append(dict(**metadata,dt_s=conf['dt_s'],summary=record['summary'],points=points,
                final_state=final,seed_dimension_m=conf['initial_equatorial_m'],
                number_flux_m2_s=flux,water_supply_kg_m2_s=conf['water_supply_kg_m2_s'],
                deposition_demand_mm_h=3600*flux*final['integrated_deposition_kg'],
                outgoing_ice_flux_mm_h=3600*flux*final['mass_kg'],
                latent_heat_removal_W_m2=flux*final['latent_heat_J'],
                interpretation='Fluxes describe the complete steady monodisperse passage only when the freezing level is reached. Censored cases give cumulative demand only.'))
        elif entry['group']=='parity':
            oldpath = trace.CACHE/'traces'/f"IGF2plus_{entry['mode']}_dt0.015625"
            old = trace.verify_run(oldpath)
            oldrows = {r['elapsed_s']:r for r in old['snapshots']}
            equal = all(all(row[k]==oldrows[row['elapsed_s']][k] for k in quantities) for row in record['snapshots'])
            if not equal:
                raise ValueError('Lunar caller changed the Earth growth control')
            parity.append(dict(**metadata,baseline_record_sha256=native.sha(oldpath/'run.json'),
                               compared_fields=quantities,exact_equality=equal))
        elif entry['group']=='drift':
            first,last = record['snapshots'][0],record['snapshots'][-1]
            expected = entry['configuration']['cloud_depth_m']/first['terminal_speed_m_s']
            if not math.isclose(last['elapsed_s'],expected,rel_tol=1e-12) or last['mass_kg']!=first['mass_kg']:
                raise ValueError('Constant-state settling control fails')
            drift.append(dict(**metadata,expected_exit_s=expected,actual_exit_s=last['elapsed_s'],
                              unchanged_mass=last['mass_kg']==first['mass_kg']))
    for entry in entries:
        dt = entry['configuration']['dt_s']
        if entry['group']=='fixed' and dt==.0078125:
            coarse_name = entry['name'].replace('dt0.0078125','dt0.03125')
        elif entry['group'] in ('fall','flux') and dt==.5:
            coarse_name = entry['name'].replace('dt0.5','dt2.0')
        else:
            continue
        a,b = records[coarse_name],records[entry['name']]
        if not (a['complete'] and b['complete']):
            continue
        av,bv = a['snapshots'][-1],b['snapshots'][-1]
        keys = (*model.QUANTITIES,'elapsed_s')
        refinements.append(dict(coarse=coarse_name,fine=entry['name'],
            same_termination=a['summary']['status_code']==b['summary']['status_code'],
            relative_change={k:bv[k]/av[k]-1 for k in keys}))
    profile = dict(gravity_m_s2=MOON_SURFACE_GRAVITY,earth_control_gravity_m_s2=source_gravity(),
        freezing_height_m=float(column.z_freeze),freezing_pressure_Pa=column.p_freeze,
        depth_m=float(column.height_of(258.15)-column.z_freeze),lapse_K_m=float(column.lapse),
        pressure_exponent=float(column.exponent),source_binding=binding,
        prior_population_product=dict(path=str(source_path.relative_to(ROOT)),sha256=native.sha(source_path)),
        source_flux_m2_s=prior['settings']['original_settings']['injected_number_flux_m2_s'],
        water_allowance_mm_h=prior['settings']['original_settings']['snowfall_mm_h'])
    return dict(schema='terluna.atmosphere.mcsnow-lunar-growth/1',
        producer={str(p.relative_to(ROOT)):native.sha(p) for p in
                  (Path(__file__),DRIVER,Path(trace.__file__),Path(model.__file__),Path(flakes.__file__),ROOT/'shared/constants.json')},
        evidence='Numerical conditional lunar transfer of the improved single-crystal control; no calibrated lunar snow distribution.',
        reading_rule='Separate fixed-age growth, unconstrained continuously saturated descent, and supply-limited steady monodisperse descent. The one-metre stop censors runaway calculations and is not a physical size limit. The imposed source flux sets the supply per crystal; seed histories are sensitivities. Compare complete paths at the freezing level, not at the ground.',
        profile=profile,builds={f'{mode}_{world}':verify_build(CACHE/f'build_{mode}_{world}')
                                for mode,world in itertools.product(MODES,('earth','moon'))},
        run_manifest=manifest,fixed_cases=fixed,falling_cases=paths,refinement=refinements,earth_parity=parity,drift_controls=drift,
        limitations=['Gravity is constant at the shared lunar surface value; the inherited air transport coefficients remain native McSnow.',
            'The fixed-age Earth controls retain native 9.81 m/s2 to preserve the identified benchmark.',
            'The hydrostatic comparison rescales depth with gravity; crossed depth/gravity controls prescribe nonhydrostatic air deliberately.',
            'Temperature and nominal water saturation are imposed. The supply gate represents an external deposition allowance, not a solved vapour reservoir.',
            'Number flux is constant along each monodisperse path: no aggregation, collisions, fragmentation, liquid collection or turbulent transport.',
            'Native latent heat release is recorded as required heat removal; no atmospheric energy feedback is integrated.',
            'The original Earth archive configuration and residual Earth mass/axis errors remain unresolved.',
            'No distribution quantiles, cloud occurrence probability, fracture survival or optical forecast follow from these trajectories.'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['run','assemble'])
    parser.add_argument('--group',action='append')
    args = parser.parse_args()
    if args.action=='run':
        for mode,world in itertools.product(MODES,('earth','moon')):
            build(mode,world)
        for entry in cases():
            if args.group and entry['group'] not in args.group:
                continue
            result = run(entry)
            print(entry['name'],result['complete'],result['summary'],flush=True)
    else:
        model.save(PRODUCT,assemble())


if __name__=='__main__':
    main()
