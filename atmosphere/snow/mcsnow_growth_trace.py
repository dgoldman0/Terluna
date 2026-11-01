"""Trace the identified 2024 crystal growth source and separate rate from shape.

The original experiment's executable remains unidentified. Source interventions
are diagnostic candidates, not accepted physics or fitted parameter choices.
"""
from __future__ import annotations

import argparse
import difflib
import hashlib
import itertools
import json
import math
from pathlib import Path
import shutil
import subprocess

import numpy as np
from atmosphere.snow import mcsnow_reproduction as model
from atmosphere.snow import mcsnow_reproduction_benchmarks as benchmark
from atmosphere.snow import mcsnow_ventilation_controls as controls

native = model.native
CACHE = model.CACHE / 'growth_diagnosis'
PRODUCT = model.ROOT / 'atmosphere/snow/results/benchmark_growth_trace.json'
FIELDS = controls.TRACE_FIELDS
MODES = {
    'native': 'Unchanged growth arithmetic, with source diagnostics.',
    'historical': 'The two documented historical mean-ventilation deletions.',
    'plate_base': 'Remove only the extra oblate mean-ventilation term; retain its equivalent-size base.',
    'sphere_radius': 'Use the physical envelope-volume equivalent sphere radius in directional ventilation.',
}
replace = model.replace_once


def source_patch(body, mode):
    if mode not in MODES:
        raise ValueError('Unknown diagnostic candidate')
    if hashlib.sha256(body.encode()).hexdigest() != model.DEPOSITION_SHA:
        raise ValueError('Unrecognized 2024 deposition source')
    if mode == 'historical':
        body = model.historical_deposition(body)
    if mode == 'plate_base':
        body = replace(body, '    fv      = fv + c2 / phi * X_equiv**(3./2.)',
                       '    ! Diagnostic: omit only the extra oblate mean-ventilation term.')
    body = replace(body, '  INTEGER, PARAMETER :: subcycles = 1',
                   '  INTEGER, PARAMETER :: subcycles = 1\n  REAL(wp), PUBLIC :: growth_trace(22)=0._wp')
    start = body.index('SUBROUTINE vapourDiffusion( sp,')
    end = body.index('END SUBROUTINE vapourDiffusion', start)
    part = body[start:end]
    part = replace(part, '  C = capacitance( sp, Re )', '''  C = capacitance( sp, Re )
  growth_trace=0._wp
  growth_trace(1)=sp%m_i
  growth_trace(4)=sp%phi
  growth_trace(5)=C
  growth_trace(6)=Re
  growth_trace(7)=Dv
  growth_trace(8)=Ka
  growth_trace(17)=vt
  growth_trace(18)=sqrtRe*cubrSc
  growth_trace(19)=sqrtRe*cubrPr''')
    anchor = '  Hi = (L/(Rv*atmo%T) - 1.0_wp)*(L*Dv*es)/(Ka*Rv*atmo%T**2)'
    part = replace(part, anchor, '''  growth_trace(9)=Dv
  growth_trace(10)=Ka
  growth_trace(11)=Dv/growth_trace(7)
  growth_trace(12)=Ka/growth_trace(8)
'''+anchor+'\n  growth_trace(22)=Hi')
    part = replace(part, '  ! sum heat for possible melting at surface temperature T_3',
                   '  growth_trace(21)=dMdt\n  ! sum heat for possible melting at surface temperature T_3')
    body = body[:start]+part+body[end:]
    start = body.index('SUBROUTINE habit_predict(')
    end = body.index('END SUBROUTINE habit_predict', start)
    part = body[start:end]
    part = replace(part, '  CALL spheroid_geometry(sp, r_a, r_c, V_tot)', '''  CALL spheroid_geometry(sp, r_a, r_c, V_tot)
  growth_trace(2)=2._wp*r_a
  growth_trace(3)=2._wp*r_c
  growth_trace(13)=IGF''')
    part = replace(part, '  dVidt    = dMdt / rho_depo', '''  growth_trace(15)=rho_depo
  growth_trace(16)=MERGE(1._wp,0._wp,branch_crit)
  dVidt    = dMdt / rho_depo''')
    anchor = '''  IGF_vent = IGF * ventilation_coef_IGF( ReSc, r_c, V_tot**z1_3, sp%phi ) /                       &
                   ventilation_coef_IGF( ReSc, r_a, V_tot**z1_3, sp%phi ) '''
    radius = '(z3_4pi*V_tot)**z1_3' if mode == 'sphere_radius' else 'V_tot**z1_3'
    part = replace(part, anchor, anchor.replace('V_tot**z1_3', radius)+
                   '\n  growth_trace(14)=IGF_vent\n  growth_trace(20)='+radius)
    return body[:start]+part+body[end:]


def trace_driver():
    body = model.reconstruction_driver(native.DRIVER.read_text())
    body = replace(body, '  use mo_depo, only: sp_depo', '  use mo_depo, only: sp_depo, growth_trace')
    body = replace(body, '  integer :: branch_unit', '  integer :: branch_unit, trace_unit\n  real(wp) :: old_branch')
    header = ['evaluation_elapsed_s', 'state_elapsed_s', *FIELDS]
    # Split the generated Fortran string at field boundaries to respect line limits.
    chunks = [','.join(header[i:i+3])+(',' if i+3<len(header) else '') for i in range(0,len(header),3)]
    statement = "  write(trace_unit,'(a)') "+'// &\n    '.join(repr(c) for c in chunks)
    body = replace(body, '  call write_state()\n  call branch_diagnostic(.true.)',
        "  open(newunit=trace_unit,file='growth_trace.csv',status='new',action='write')\n"+
        statement+'\n  old_branch=-1._wp\n  call write_state()\n  call branch_diagnostic(.true.)')
    body = replace(body, '    elapsed=k*dt_s', '''    elapsed=k*dt_s
    if (mod(k,output_steps)==0.or.k==1.or.k==nsteps.or.growth_trace(16)/=old_branch) then
      write(trace_unit,'(*(es25.17e3,:,","))') elapsed-dt_s,elapsed,growth_trace
    end if
    old_branch=growth_trace(16)''')
    return replace(body, '  close(branch_unit)', '  close(branch_unit)\n  close(trace_unit)')


def build(mode):
    path = CACHE / ('build_'+mode)
    baseline = native.verify_build(model.BASELINE)
    if path.exists():
        return verify_build(path)
    original = (model.BASELINE/'source/src/mo_depo.f90').read_text()
    changed = source_patch(original, mode)
    path.mkdir(parents=True)
    shutil.copytree(model.BASELINE/'source', path/'source')
    (path/'source/src/mo_depo.f90').write_text(changed)
    (path/native.DRIVER.name).write_text(trace_driver())
    for name in native.TABLE_NAMES:
        shutil.copyfile(model.BASELINE/name, path/name)
    diff = ''.join(difflib.unified_diff(original.splitlines(True), changed.splitlines(True),
                                      fromfile='native/mo_depo.f90', tofile=mode+'/mo_depo.f90'))
    diff += ''.join(difflib.unified_diff(model.reconstruction_driver(native.DRIVER.read_text()).splitlines(True),
        trace_driver().splitlines(True), fromfile='reconstruction/driver.f90', tofile='trace/driver.f90'))
    (path/'diagnostic.patch').write_text(diff)
    compiler = baseline['compiler_path']
    selected = native.selected_sources(path/'source/src', path/native.DRIVER.name, compiler)
    commands, objects = [], []
    with (path/'build.log').open('w') as log:
        for i, src in enumerate(selected):
            obj = path/f'{i:02d}_{src.stem}.o'; objects.append(str(obj))
            cmd = [compiler, *native.FLAGS, '-I'+str(path/'source/src'), '-I'+str(path),
                   '-J'+str(path), '-c', str(src), '-o', str(obj)]
            commands.append(cmd)
            subprocess.run(cmd, cwd=path, env=native.ENV, stdout=log, stderr=subprocess.STDOUT, check=True)
        cmd = [compiler, *native.FLAGS, *objects, '-o', str(path/'crystal')]
        commands.append(cmd)
        subprocess.run(cmd, cwd=path, env=native.ENV, stdout=log, stderr=subprocess.STDOUT, check=True)
    record = dict(schema='terluna.atmosphere.mcsnow-growth-trace-build/1', mode=mode,
        baseline_record_sha256=native.sha(model.BASELINE/'build.json'), source_identity=baseline['source_identity'],
        source_files={str(p.relative_to(path/'source')): native.sha(p)
                      for p in sorted((path/'source').rglob('*')) if p.is_file()},
        driver_sha256=native.sha(path/native.DRIVER.name), executable_sha256=native.sha(path/'crystal'),
        patch_sha256=native.sha(path/'diagnostic.patch'), build_log_sha256=native.sha(path/'build.log'),
        table_identity=baseline['table_identity'], compiler_path=compiler, compiler_sha256=native.sha(compiler),
        compiler_version=baseline['compiler_version'], flags=native.FLAGS, commands=commands,
        evidence=MODES[mode])
    model.save(path/'build.json', record)
    return record


def verify_build(path):
    path = Path(path)
    record = native.verify_build(path)
    native.verify_build(model.BASELINE)
    if native.sha(model.BASELINE/'build.json') != record['baseline_record_sha256']:
        raise ValueError('Baseline record changed')
    for name,key in [('diagnostic.patch','patch_sha256'), ('build.log','build_log_sha256')]:
        if native.sha(path/name) != record[key]:
            raise ValueError('Build artifact changed: '+name)
    expected = source_patch((model.BASELINE/'source/src/mo_depo.f90').read_text(), record['mode'])
    if (path/'source/src/mo_depo.f90').read_text() != expected or (path/native.DRIVER.name).read_text() != trace_driver():
        raise ValueError('Diagnostic source or caller differs')
    return record


def validate_trace(trace, snapshots, config):
    if not trace:
        raise ValueError('Missing source trace')
    dt = config['dt_s']
    times = [r['state_elapsed_s'] for r in trace]
    if times != sorted(set(times)) or times[0] != dt or times[-1] != config['duration_s']:
        raise ValueError('Incomplete or duplicate trace times')
    output_times = [s['elapsed_s'] for s in snapshots[1:]]
    if not set(output_times).issubset(times):
        raise ValueError('Missing saved source evaluation')
    states = {s['elapsed_s']: s for s in snapshots}
    for row in trace:
        if row['evaluation_elapsed_s'] != row['state_elapsed_s']-dt:
            raise ValueError('Incorrect source evaluation time')
        if not math.isclose(row['state_elapsed_s']/dt,round(row['state_elapsed_s']/dt),abs_tol=1e-9):
            raise ValueError('Source evaluation is off the integration grid')
        state = states.get(row['state_elapsed_s'])
        if row['branch_criterion_met'] not in (0.,1.):
            raise ValueError('Invalid branch state')
        if min(row[k] for k in FIELDS if k != 'branch_criterion_met') <= 0:
            raise ValueError('Nonpositive diagnostic')
        if state and not math.isclose(row['mass_kg']+dt*row['deposition_rate_kg_s'],state['mass_kg'],rel_tol=2e-13):
            raise ValueError('Trace/state deposition mismatch')
        if state:
            volume = math.pi/6*row['equatorial_dimension_m']**2*row['polar_dimension_m']
            expected_volume = volume+dt*row['deposition_rate_kg_s']/row['deposition_density_kg_m3']
            gamma = row['effective_growth_factor']
            expected_phi = row['aspect_ratio']*(expected_volume/volume)**((gamma-1)/(gamma+2))
            if not math.isclose(state['ice_volume_m3'],expected_volume,rel_tol=5e-13):
                raise ValueError('Trace/state volume mismatch')
            if not math.isclose(state['aspect_ratio'],expected_phi,rel_tol=5e-13):
                raise ValueError('Trace/state habit mismatch')
        for eff,raw,factor in [('diffusivity_effective_m2_s','diffusivity_raw_m2_s','ventilation_vapor'),
                               ('conductivity_effective_W_m_K','conductivity_raw_W_m_K','ventilation_heat')]:
            if not math.isclose(row[eff]/row[raw],row[factor],rel_tol=2e-14):
                raise ValueError('Trace ventilation mismatch')


def read_run(path, config, returncode):
    rows, air, errors = model.analyze_output(path, config, returncode)
    trace, trace_errors = native.read_csv(path/'growth_trace.csv', ('evaluation_elapsed_s','state_elapsed_s',*FIELDS))
    errors += trace_errors
    if not errors and returncode == 0:
        validate_trace(trace, rows, config)
    return dict(snapshots=rows, air=air, trace=trace, validation_errors=errors,
                complete=returncode == 0 and not errors)


def run(mode, name, config):
    path = CACHE/'traces'/name
    build_path = CACHE/('build_'+mode)
    verify_build(build_path)
    if path.exists():
        record = verify_run(path)
        if record['configuration'] != config or record['mode'] != mode:
            raise ValueError('Existing run uses different inputs')
        return record
    path.mkdir(parents=True)
    (path/'crystal.nml').write_text(native.namelist(config))
    for name in native.TABLE_NAMES:
        shutil.copyfile(build_path/name, path/name)
    with (path/'run.log').open('w') as log:
        proc = subprocess.run([str(build_path/'crystal'),'crystal.nml'],cwd=path,env=native.ENV,stdout=log,stderr=subprocess.STDOUT)
    record = dict(schema='terluna.atmosphere.mcsnow-growth-trace-run/1', mode=mode, configuration=config,
        build_record_sha256=native.sha(build_path/'build.json'), returncode=proc.returncode,
        **read_run(path,config,proc.returncode),
        file_sha256={p.name:native.sha(p) for p in sorted(path.iterdir()) if p.is_file()})
    model.save(path/'run.json',record)
    return record


def verify_run(path):
    path = Path(path)
    record = json.loads((path/'run.json').read_text())
    build_path = CACHE/('build_'+record['mode'])
    verify_build(build_path)
    if native.sha(build_path/'build.json') != record['build_record_sha256']:
        raise ValueError('Run/build identity changed')
    for name,digest in record['file_sha256'].items():
        if native.sha(path/name) != digest:
            raise ValueError('Raw run changed: '+name)
    if (path/'crystal.nml').read_text() != native.namelist(record['configuration']):
        raise ValueError('Namelist differs from configuration')
    replay = read_run(path, record['configuration'],record['returncode'])
    if any(record[k] != v for k,v in replay.items()):
        raise ValueError('Saved run differs from raw replay')
    return record


def matched_mass(rows, mass):
    """Interpolate only within a strictly growing, actually computed trajectory."""
    masses = np.array([r['mass_kg'] for r in rows])
    if not np.all(np.diff(masses)>0):
        raise ValueError('Mass must increase strictly')
    if mass < masses[0] or mass > masses[-1]:
        return None
    i = max(1, int(np.searchsorted(masses,mass)))
    lo,hi = rows[i-1],rows[i]
    fraction = (mass-lo['mass_kg'])/(hi['mass_kg']-lo['mass_kg'])
    return dict(state={k:lo[k]+fraction*(hi[k]-lo[k]) for k in lo},
                bracket_elapsed_s=[lo['elapsed_s'],hi['elapsed_s']])


def branch_mass(mass, volume, rho_ice, gamma):
    """Mass at the single solid-to-constant-porous-density transition, if assumed."""
    if not (0 < gamma < 1) or rho_ice <= 0 or min(mass,volume)<=0:
        raise ValueError('Requires a constant-temperature branching plate')
    return (mass-rho_ice*gamma*volume)/(1-gamma)


def cases():
    for variant,growth in benchmark.VARIANTS:
        for mode,dt in itertools.product(MODES, (1.,.0625,.03125,.015625)):
            config = model.configuration(258.65,dt=dt,growth=growth)
            config['output_interval_s'] = 1. if dt == 1. else .25
            yield dict(name=f'{variant}_{mode}_dt{dt}',mode=mode,variant=variant,group='seed',configuration=config)
        for age in (180,600):
            author = benchmark.target_row(model.author_data(variant=variant,duration=age),258.65)
            for mode in MODES:
                config = benchmark.archived_configuration(author)
                config.update(growth_mode=growth,duration_s=1.,dt_s=.015625,output_interval_s=1.)
                yield dict(name=f'fixed{age}_{variant}_{mode}',mode=mode,variant=variant,
                           group='fixed_state',author_age_s=age,configuration=config)
    # Continue the rate-control finding through the laboratory's measured interval.
    for dt in (.03125,.0078125):
        config = model.configuration(258.75,dt=dt,duration=1800.)
        config['output_interval_s'] = 1.
        yield dict(name=f'extended_plate_base_igf2_Tc-14.4_dt{dt}',mode='plate_base',variant='IGF2plus',
                   group='extension',configuration=config)
    for mode in ('native','plate_base'):
        config = model.configuration(258.65,dt=.015625)
        config['output_interval_s'] = .0625
        yield dict(name=f'interpolation_IGF2plus_{mode}',mode=mode,variant='IGF2plus',
                   group='interpolation',configuration=config)


def restart_cases():
    for variant,growth in benchmark.VARIANTS:
        author = benchmark.target_row(model.author_data(variant=variant,duration=180),258.65)
        for mode,dt in itertools.product(('native','historical'),(.0625,.015625)):
            config = benchmark.archived_configuration(author)
            config.update(growth_mode=growth,duration_s=420.,dt_s=dt)
            yield dict(name=f'restart180_{variant}_{mode}_dt{dt}',mode=mode,variant=variant,configuration=config)


def assemble():
    authors, sources = {}, {}
    for variant,age in itertools.product(dict(benchmark.VARIANTS),(180,600)):
        source = model.author_data(variant=variant,duration=age)
        authors[variant,age] = benchmark.target_row(source,258.65)
        sources[source['member']] = {k:v for k,v in source.items() if k != 'rows'}
    records,manifest = {},[]
    builds = {mode:verify_build(CACHE/('build_'+mode)) for mode in MODES}
    comparisons,matched,fixed,branches,selected_traces,extensions,parity = [],[],[],[],[],[],[]
    for entry in cases():
        path = CACHE/'traces'/entry['name']
        record = verify_run(path)
        config = entry['configuration']
        if record['configuration'] != config or record['mode'] != entry['mode']:
            raise ValueError('Matrix/configuration mismatch')
        if not record['complete']:
            raise ValueError('Incomplete diagnostic run: '+entry['name'])
        records[entry['name']] = record
        common = {k:entry[k] for k in ('name','mode','variant','group')}
        common['dt_s'] = config['dt_s']
        manifest.append(dict(**entry,path=str(path.relative_to(model.ROOT)),
            run_record_sha256=native.sha(path/'run.json'),complete=record['complete'],
            snapshots=len(record['snapshots']),source_evaluations=len(record['trace'])))
        if entry['group'] == 'seed':
            for age in (180,600):
                author = authors[entry['variant'],age]
                state = benchmark.endpoint(record,age)
                comparisons.append(dict(**common,age_s=age,state=state,comparison=model.compare(state,author)))
                match = matched_mass(record['snapshots'],author['values']['mass_kg'])
                matched.append(dict(**common,author_age_s=age,matched=match,
                    comparison=model.compare(match['state'],author) if match else None))
                selected_traces.append(dict(**common,age_s=age,
                    trace=next(r for r in record['trace'] if r['state_elapsed_s']==age)))
            first_branch = next(r for r in record['trace'] if r['branch_criterion_met'])
            final = record['snapshots'][-1]
            rho,gamma = record['air'][0]['native_ice_density_kg_m3'],record['air'][0]['inherent_growth_factor']
            author = authors[entry['variant'],600]
            am,av = (author['values'][k] for k in ('mass_kg','ice_volume_m3'))
            em,ev = (author['rounding_half_width'][k] for k in ('mass_kg','ice_volume_m3'))
            inferred = branch_mass(am,av,rho,gamma)
            branches.append(dict(**common,first_branch_evaluation=first_branch,
                native_final_inferred_branch_mass_kg=branch_mass(final['mass_kg'],final['ice_volume_m3'],rho,gamma),
                author_conditional_branch_mass_kg=inferred,
                author_printed_rounding_half_bound_kg=(em+rho*gamma*ev)/(1-gamma),
                native_to_author_conditional_mass_ratio=first_branch['mass_kg']/inferred,
                assumption='One irreversible solid-to-rho_i*Gamma deposition-density transition; same constant native Gamma as the archive.'))
            if entry['mode'] in ('native','historical'):
                oldname = (f"{entry['variant']}_{entry['mode']}_T258.65_sat1" if config['dt_s']==1 else
                           f"refine_{entry['variant']}_{entry['mode']}_T258.65_dt{config['dt_s']}")
                oldpath = model.CACHE/'runs'/oldname
                prior = model.verify_run(oldpath,model.CACHE/('build_'+entry['mode']))
                states = {r['elapsed_s']:r for r in record['snapshots']}
                equal = all(states[r['elapsed_s']]==r for r in prior['snapshots']) and prior['air']==record['air']
                if not equal:
                    raise ValueError('Trace changed the preceding trajectory: '+entry['name'])
                parity.append(dict(**common,baseline=str(oldpath.relative_to(model.ROOT)),
                    baseline_record_sha256=native.sha(oldpath/'run.json'),exact_common_snapshot_and_air_equality=equal))
        elif entry['group'] == 'fixed_state':
            fixed.append(dict(**common,author_age_s=entry['author_age_s'],trace=record['trace'][0]))
        elif entry['group'] == 'extension':
            points = []
            for age in (180,600,1200,1800):
                state = benchmark.endpoint(record,age)
                observed = benchmark.lab.observed(-14.4,age/60)
                points.append(dict(age_s=age,state=state,laboratory_fit=observed,
                    fit_ratios={k:state[k]/observed[k] for k in native.QUANTITIES}))
            extensions.append(dict(**common,temperature_C=-14.4,points=points))
    restarts = []
    for entry in restart_cases():
        path = CACHE/entry['name']
        record = model.verify_run(path,model.CACHE/('build_'+entry['mode']))
        if record['configuration'] != entry['configuration'] or not record['complete']:
            raise ValueError('Restart configuration/completion mismatch')
        restarts.append(dict(**entry,path=str(path.relative_to(model.ROOT)),
            run_record_sha256=native.sha(path/'run.json'),final_state=record['snapshots'][-1],
            comparison=model.compare(record['snapshots'][-1],authors[entry['variant'],600]),
            initial_age_s=180,final_age_s=600,complete=record['complete']))
    refinements,shape_paths,interpolation = [],[],[]
    for variant,mode in itertools.product(dict(benchmark.VARIANTS),MODES):
        coarse,fine = (records[f'{variant}_{mode}_dt{dt}'] for dt in (.03125,.015625))
        changes = []
        for age in (180,600):
            target = authors[variant,age]['values']['mass_kg']
            a,b = matched_mass(coarse['snapshots'],target),matched_mass(fine['snapshots'],target)
            changes.append(dict(author_age_s=age,
                endpoint_relative_change=benchmark.endpoint_change(benchmark.endpoint(coarse,age),benchmark.endpoint(fine,age)),
                matched_relative_change=benchmark.endpoint_change(a['state'],b['state']) if a and b else None,
                matched_time_change_s=b['state']['elapsed_s']-a['state']['elapsed_s'] if a and b else None))
        refinements.append(dict(variant=variant,mode=mode,coarse_dt_s=.03125,fine_dt_s=.015625,points=changes))
    for variant,dt in itertools.product(dict(benchmark.VARIANTS),(.0625,.015625)):
        mass = authors[variant,600]['values']['mass_kg']
        a,b = (matched_mass(records[f'{variant}_{mode}_dt{dt}']['snapshots'],mass) for mode in ('native','plate_base'))
        shape_paths.append(dict(variant=variant,dt_s=dt,mass_kg=mass,
            plate_base_relative_to_native={k:b['state'][k]/a['state'][k]-1 for k in
                ('equatorial_dimension_m','polar_dimension_m','ice_volume_m3','aspect_ratio')},
            native_age_s=a['state']['elapsed_s'],plate_base_age_s=b['state']['elapsed_s']))
    for mode in ('native','plate_base'):
        coarse,fine = records[f'IGF2plus_{mode}_dt0.015625'],records[f'interpolation_IGF2plus_{mode}']
        states = {r['elapsed_s']:r for r in fine['snapshots']}
        equal = all(states[r['elapsed_s']]==r for r in coarse['snapshots'])
        if not equal:
            raise ValueError('Output cadence changed the computed trajectory')
        points = []
        for age in (180,600):
            mass = authors['IGF2plus',age]['values']['mass_kg']
            a,b = (matched_mass(r['snapshots'],mass)['state'] for r in (coarse,fine))
            points.append(dict(author_age_s=age,relative_change=benchmark.endpoint_change(a,b),
                               time_change_s=b['elapsed_s']-a['elapsed_s']))
        interpolation.append(dict(mode=mode,coarse_output_interval_s=.25,fine_output_interval_s=.0625,
                                  common_snapshots_exactly_equal=equal,points=points))
    return dict(schema='terluna.atmosphere.mcsnow-growth-trace/1',
        producer_sha256={p.name:native.sha(p) for p in (Path(__file__),Path(model.__file__),Path(benchmark.__file__),Path(controls.__file__))},
        evidence='Actual native source terms, declared source interventions and conditional numerical-archive comparisons; no adopted repair.',
        source_archives=sources,builds=builds,modes=MODES,run_manifest=manifest,restarts=restarts,
        comparisons=comparisons,matched_mass=matched,branch_mass=branches,source_traces=selected_traces,
        fixed_state=fixed,refinement=refinements,interpolation=interpolation,shape_paths=shape_paths,extensions=extensions,parity=parity,
        limits=['Original figure driver/configuration remains unidentified.',
                'Earliest matching archive is 180 s; departure cannot be dated earlier from these outputs.',
                'Matched-mass states interpolate computed trajectories; they are not age-matched reproductions.',
                'Conditional branch mass assumes the native deposition-density law in the archive.',
                'Restarts use rounded numerical-archive states, not tracked laboratory crystals.',
                'Prescribed water-saturated air has no finite vapour reservoir, heat feedback or population processes.'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['run','assemble'])
    args = parser.parse_args()
    if args.action == 'run':
        for mode in MODES:
            build(mode)
        for entry in cases():
            result = run(entry['mode'],entry['name'],entry['configuration'])
            print(entry['name'],result['complete'],flush=True)
        for entry in restart_cases():
            result = model.run(model.CACHE/('build_'+entry['mode']),CACHE/entry['name'],entry['configuration'])
            print(entry['name'],result['complete'],flush=True)
    else:
        model.save(PRODUCT,assemble())


if __name__ == '__main__':
    main()
