"""Native collection rates for fixed laboratory and model crystal properties.

These are conditional pair calculations, not aggregation observations or a
population simulation. Growth remains disabled throughout every pair test.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path
import shutil
import subprocess

from atmosphere.snow import mcsnow_shape_diagnosis as shape
from atmosphere.snow import welss_aggregation_physics as independent

native = shape.native
model = shape.model
HERE = Path(__file__).resolve().parent
DRIVER = HERE / 'mcsnow_pair_driver.f90'
CACHE = shape.CACHE / 'pairs'
PRODUCT = HERE / 'results/mcsnow_pair_controls.json'


def build():
    path = CACHE / 'build'
    baseline = native.verify_build(model.BASELINE)
    if path.exists():
        return verify_build()
    path.mkdir(parents=True)
    shutil.copyfile(DRIVER, path/DRIVER.name)
    source = model.BASELINE/'source/src'
    compiler = baseline['compiler_path']
    selected = native.selected_sources(source, path/DRIVER.name, compiler)
    commands, objects = [], []
    with (path/'build.log').open('w') as log:
        for index, src in enumerate(selected):
            obj = path/f'{index:02d}_{src.stem}.o'
            objects.append(str(obj))
            command = [compiler, *native.FLAGS, '-I'+str(source), '-I'+str(path), '-J'+str(path),
                       '-c', str(src), '-o', str(obj)]
            commands.append(command)
            subprocess.run(command, cwd=path, env=native.ENV, stdout=log, stderr=subprocess.STDOUT, check=True)
        command = [compiler, *native.FLAGS, *objects, '-o', str(path/'pairs')]
        commands.append(command)
        subprocess.run(command, cwd=path, env=native.ENV, stdout=log, stderr=subprocess.STDOUT, check=True)
    record = dict(schema='terluna.atmosphere.mcsnow-pair-build/1',
                  baseline_record_sha256=native.sha(model.BASELINE/'build.json'),
                  source_identity=baseline['source_identity'],
                  compiled_sources={str(p):native.sha(p) for p in selected},
                  compiler_path=compiler, compiler_sha256=native.sha(compiler),
                  compiler_version=baseline['compiler_version'], flags=native.FLAGS,
                  commands=commands, physics_patches=[],
                  file_sha256={p.name:native.sha(p) for p in (path/DRIVER.name,path/'pairs',path/'build.log')})
    model.save(path/'build.json',record)
    return record


def verify_build():
    path = CACHE/'build'
    record = json.loads((path/'build.json').read_text())
    native.verify_build(model.BASELINE)
    if native.sha(model.BASELINE/'build.json') != record['baseline_record_sha256']:
        raise ValueError('Baseline identity differs')
    for filename, digest in record['compiled_sources'].items():
        if native.sha(filename) != digest:
            raise ValueError('Compiled source differs: '+filename)
    for filename, digest in record['file_sha256'].items():
        if native.sha(path/filename) != digest:
            raise ValueError('Build artifact differs: '+filename)
    if native.sha(DRIVER) != record['file_sha256'][DRIVER.name]:
        raise ValueError('Pair driver differs')
    return record


def particle(mass, a, c):
    if min(mass,a,c)<=0 or not all(math.isfinite(v) for v in (mass,a,c)):
        raise ValueError('Positive finite crystal properties required')
    return dict(mass_kg=mass, volume_m3=math.pi/6*a*a*c, phi=c/a)


def states():
    data = json.loads(shape.PRODUCT.read_text())
    rho = data['native_ice_density_kg_m3']
    rows = []
    for temperature in shape.lab.FITS:
        ages = shape.sample_ages(temperature)
        for mapping in shape.MAPPINGS:
            for age in ages:
                fit = shape.fit_state(temperature, age, mapping)
                if fit['density_kg_m3'] > rho:
                    continue
                rows.append(dict(temperature_C=temperature, age_min=age, family='fitted_'+mapping,
                                 state=particle(fit['mass_kg'],fit['model_a_m'],fit['model_c_m'])))
    for run in data['trajectories']:
        if run['dt_s'] != min(shape.TIMESTEPS):
            continue
        for row in run['comparisons']:
            value = row['age_matched']['values']
            # Undo only the reporting map: native collision geometry uses actual spheroid axes.
            scale = run['axis_scale']
            rows.append(dict(temperature_C=run['temperature_C'],age_min=row['age_min'],
                             family='model_'+run['initialization'],
                             state=particle(value['mass_kg'],value['equatorial_dimension_m']*scale,
                                            value['polar_dimension_m']*scale)))
    return rows


def cases():
    grouped = {}
    for row in states():
        grouped.setdefault((row['temperature_C'],row['family']),[]).append(row)
    rows = []
    for (temperature,family), entries in sorted(grouped.items()):
        for left,right in itertools.combinations_with_replacement(sorted(entries,key=lambda r:r['age_min']),2):
            for sticking in (0,1):
                rows.append(dict(group='shape',temperature_C=temperature,family=family,
                                 age_left_min=left['age_min'],age_right_min=right['age_min'],
                                 habit=1,sticking=sticking,left=left['state'],right=right['state']))
    # Independent implementation parity covers the selected older diagnostic-geometry closure.
    masses = (1e-12,1e-10,1e-8,1e-6)
    for temperature in (-5.3,-14.4,-20.1):
        for left,right in itertools.combinations_with_replacement(masses,2):
            rows.append(dict(group='independent_parity',temperature_C=temperature,
                             family='diagnostic_mass_geometry',habit=0,sticking=1,
                             left=dict(mass_kg=left,volume_m3=left/independent.w.RHOI,phi=1.),
                             right=dict(mass_kg=right,volume_m3=right/independent.w.RHOI,phi=1.)))
    for index,row in enumerate(rows):
        row['case_id']=index
    return rows


def input_text(rows):
    result=[]
    for row in rows:
        values=[row['case_id'],273.15+row['temperature_C'],101000.,row['habit'],row['sticking']]
        for side in ('left','right'):
            values.extend(row[side][k] for k in ('mass_kg','volume_m3','phi'))
        result.append(' '.join(format(v,'.17g') for v in values))
    return '\n'.join(result)+'\n'


def execute():
    build()
    rows = cases()
    path = CACHE/'execution'
    if path.exists():
        return verify_run(rows)
    path.mkdir()
    (path/'input.txt').write_text(input_text(rows))
    with (path/'pairs.csv').open('w') as out, (path/'run.log').open('w') as log:
        result=subprocess.run([str(CACHE/'build/pairs'),str(path/'input.txt')],cwd=path,env=native.ENV,
                              stdout=out,stderr=log)
    record=dict(returncode=result.returncode,build_record_sha256=native.sha(CACHE/'build/build.json'),
                shape_product_sha256=native.sha(shape.PRODUCT),
                file_sha256={name:native.sha(path/name) for name in ('input.txt','pairs.csv','run.log')})
    model.save(path/'run.json',record)
    return verify_run(rows)


def verify_run(inputs):
    verify_build()
    path=CACHE/'execution'
    record=json.loads((path/'run.json').read_text())
    if record['returncode'] != 0:
        raise ValueError('Failed pair execution retained; inspect run.log')
    if record['build_record_sha256'] != native.sha(CACHE/'build/build.json'):
        raise ValueError('Run/build mismatch')
    if record['shape_product_sha256'] != native.sha(shape.PRODUCT):
        raise ValueError('Growth input product differs')
    if (path/'input.txt').read_text() != input_text(inputs):
        raise ValueError('Pair input changed')
    for name,digest in record['file_sha256'].items():
        if native.sha(path/name) != digest:
            raise ValueError('Raw pair output differs: '+name)
    rows,errors=native.read_csv(path/'pairs.csv',('case_id','kernel_m3_s','collection_efficiency'))
    if errors or len(rows)!=len(inputs) or [r['case_id'] for r in rows]!=list(range(len(inputs))):
        raise ValueError('Incomplete or unordered pair output')
    for row,case in zip(rows,inputs):
        for side in ('left','right'):
            if not math.isclose(row['mass_'+side+'_kg'],case[side]['mass_kg'],rel_tol=1e-13):
                raise ValueError('Pair mass input differs')
        if min(row['collection_efficiency'],row['reverse_collection_efficiency'],row['kernel_m3_s'])<0:
            raise ValueError('Negative collection result')
        value=row['collection_efficiency']*row['encounter_area_m2']*row['relative_speed_m_s']
        if not math.isclose(value,row['kernel_m3_s'],rel_tol=2e-14,abs_tol=1e-30):
            raise ValueError('Encounter ledger does not reconstruct kernel')
    return record,rows


def assemble():
    inputs=cases()
    identity,outputs=verify_run(inputs)
    rows=[dict(**case,result=result) for case,result in zip(inputs,outputs)]
    parity=[]
    for row in rows:
        if row['group']!='independent_parity':
            continue
        result=row['result']
        state=dict(T=result['temperature_K'],rho=result['air_density_kg_m3'],eta=result['dynamic_viscosity_Pa_s'])
        comparison=independent.kernel(row['left']['mass_kg'],row['right']['mass_kg'],state)
        differences={}
        for python_key,native_key in (('K','kernel_m3_s'),('Ec','collection_efficiency'),
                                       ('encounter_area_m2','encounter_area_m2'),
                                       ('relative_speed_m_s','relative_speed_m_s')):
            value=float(comparison[python_key])
            if python_key=='Ec':
                value*=float(comparison['sticking'])
            expected=result[native_key]
            differences[native_key]=dict(independent=value,native=expected,
                relative_difference=abs(value/expected-1) if expected else None,
                absolute_difference=abs(value-expected))
        parity.append(dict(case_id=row['case_id'],differences=differences))
    comparisons=[]
    lookup={(r['temperature_C'],r['family'],r.get('age_left_min'),r.get('age_right_min'),r['sticking']):r for r in rows}
    for row in rows:
        if row['group']!='shape' or not row['family'].startswith('model_') or row['sticking']!=1:
            continue
        reference_families=('fitted_axes','fitted_volume') if row['family']=='model_sphere' else ('fitted_'+row['family'][6:],)
        for family in reference_families:
            ref=lookup.get((row['temperature_C'],family,row['age_left_min'],row['age_right_min'],1))
            if ref is None:
                continue
            a,b=row['result'],ref['result']
            comparisons.append(dict(temperature_C=row['temperature_C'],family=row['family'],reference=family,
                age_left_min=row['age_left_min'],age_right_min=row['age_right_min'],
                case_id=row['case_id'],reference_case_id=ref['case_id'],
                ratio={k:a[k]/b[k] if b[k] else None for k in
                       ('kernel_m3_s','collection_efficiency','encounter_area_m2','relative_speed_m_s')}))
    producer=[Path(__file__),DRIVER,Path(shape.__file__),Path(independent.__file__),Path(independent.w.__file__),
              Path(independent.boehm.__file__)]
    return dict(schema='terluna.atmosphere.mcsnow-pair-controls/1',
        evidence='Prescribed-property native collection kernels and independent implementation controls; no aggregation population or observed collection efficiency.',
        producer={str(p.relative_to(model.ROOT)):native.sha(p) for p in producer},
        shape_product_sha256=native.sha(shape.PRODUCT),execution_identity=identity,
        build_record_sha256=native.sha(CACHE/'build/build.json'),cases=rows,
        independent_parity=parity,model_to_prescribed_geometry=comparisons,
        limits=['Particles from different fitted ages are hypothetical coexisting pairs, not laboratory collision observations.',
                'Fitted geometry still uses native projected-area and drag closures; no measured projected-area history was supplied.',
                'Temperature-only sticking options zero (unity) and one are sensitivity settings, not inferred sticking efficiencies.',
                'Raw collection efficiency is not clipped at unity; it includes native hydrodynamic interception and sticking.',
                'Pure ice, Earth gravity, fixed native water saturation and 101000 Pa; no imposed turbulent speed spread.',
                'The outer ckt0 encounter area follows the native radius formula, including the prolate aspect correction.',
                'No evolving aggregate geometry, population source, water depletion, breakup or lunar forecast is established.'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('run','assemble'))
    args=parser.parse_args()
    if args.action=='run':
        _,rows=execute(); print(len(rows),'native pair evaluations')
    else:
        product=assemble(); model.save(PRODUCT,product)
        print(len(product['cases']),'pairs;',len(product['independent_parity']),'independent parity comparisons')


if __name__=='__main__':
    main()
