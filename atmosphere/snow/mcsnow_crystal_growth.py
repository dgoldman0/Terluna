"""Execute unmodified McSnow single-crystal growth in prescribed laboratory air.

The thin Fortran caller invokes native deposition, habit, shape and drag routines.
External archives/builds remain in scratch. These runs compare native configured
physics with empirical fits; they do not reproduce an archived paper executable.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tarfile
import time
import numpy as np
from atmosphere.snow import welss_aggregation_build as pinned
from atmosphere.snow import takahashi_growth as lab
from atmosphere.snow import scale_crystal_growth as scale

HERE=Path(__file__).resolve().parent
DRIVER=HERE/'mcsnow_crystal_driver.f90'
ROOT=Path(lab.__file__).resolve().parents[2]
MANIFEST=Path(lab.__file__).with_name('growth_model_sources.json')
TABLE_NAMES=('curve_inherent.txt','IGF_new_poly.txt')
FLAGS=['-O2','-g','-cpp','-fcheck=all','-ffpe-trap=invalid,zero,overflow','-fbacktrace','-finit-real=snan']
QUANTITIES=('mass_kg','equatorial_dimension_m','polar_dimension_m','diameter_m')
ENV={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def producer():
    return {p.name:sha(p) for p in (Path(__file__),DRIVER,Path(lab.__file__),
        Path(scale.__file__),Path(scale.lunar.flakes.__file__),Path(pinned.__file__),pinned.MANIFEST,MANIFEST)}

def table_sources(source):
    manifest=json.loads(pinned.MANIFEST.read_text());result={}
    for item in manifest['files']:
        if item['path'] in ['input/'+n for n in TABLE_NAMES]:
            p=Path(source)/item['path'];pinned.verify(p.read_bytes(),item);result[p.name]=dict(item)
    if len(result)!=2:raise ValueError('Missing pinned habit tables')
    return result

def selected_sources(source,driver,compiler):
    """Compile only transitive native module dependencies of the thin caller."""
    candidates=[p for p in Path(source).rglob('*') if p.suffix.lower()=='.f90' and p.name!='main.f90']+[driver]
    provides={};needs={}
    for p in candidates:
        body=subprocess.check_output([compiler,*FLAGS,'-E','-I'+str(source),str(p)],text=True,env=ENV)
        names=re.findall(r'^\s*module\s+(?!procedure\b|function\b|subroutine\b)(\w+)',body,re.I|re.M)
        for n in names:
            if n.lower() in provides:raise ValueError('Duplicate native module '+n)
            provides[n.lower()]=p
        needs[p]={n.lower() for n in re.findall(r'^\s*use\s+(\w+)',body,re.I|re.M)}
    result=[];visiting=set()
    def visit(p):
        if p in result:return
        if p in visiting:raise ValueError('Cyclic native dependencies')
        visiting.add(p)
        for name in sorted(needs[p]):
            if name in provides and provides[name]!=p:visit(provides[name])
        visiting.remove(p);result.append(p)
    visit(driver)
    return result

def build(output,version,public_source,archive=None,compiler='gfortran'):
    output=Path(output).resolve();public_source=Path(public_source).resolve()
    if output.exists():raise FileExistsError(output)
    own=producer();table_identity=table_sources(public_source)
    old=json.loads(pinned.MANIFEST.read_text())
    output.mkdir(parents=True)
    source=output/'source';source.mkdir()
    if version=='2024':
        for item in old['files']:
            if not item['path'].startswith('src/'):continue
            p=public_source/item['path'];pinned.verify(p.read_bytes(),item)
            q=source/item['path'];q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        identity=dict(version=version,revision=old['revision'],manifest_sha256=sha(pinned.MANIFEST))
    elif version=='2026':
        archive=Path(archive).resolve();meta=json.loads(MANIFEST.read_text())['mcsnow_stochastic_geometry']
        if sha(archive)!=meta['archive_sha256']:raise ValueError('Native archive identity differs')
        with tarfile.open(archive) as t:
            for member in t.getmembers():
                if not member.isfile() or not member.name.startswith('src/') or '/._' in member.name:continue
                if '..' in Path(member.name).parts:raise ValueError('Unsafe archive path')
                p=source/member.name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(t.extractfile(member).read())
        identity=dict(version=version,archive_sha256=sha(archive),archive_version=meta['archive_version'])
    else:raise ValueError('Unidentified native source version')
    for name in TABLE_NAMES:shutil.copyfile(public_source/'input'/name,output/name)
    shutil.copyfile(DRIVER,output/DRIVER.name)
    compiler_path=str(Path(shutil.which(compiler)).resolve())
    order=selected_sources(source/'src',output/DRIVER.name,compiler_path)
    commands=[];objects=[];start=time.monotonic()
    with (output/'build.log').open('w') as log:
        for i,p in enumerate(order):
            obj=output/f'{i:02d}_{p.stem}.o';objects.append(str(obj))
            cmd=[compiler_path,*FLAGS,'-I'+str(source/'src'),'-I'+str(output),'-J'+str(output),'-c',str(p),'-o',str(obj)]
            commands.append(cmd);subprocess.run(cmd,cwd=output,env=ENV,stdout=log,stderr=subprocess.STDOUT,check=True)
        cmd=[compiler_path,*FLAGS,*objects,'-o',str(output/'crystal')]
        commands.append(cmd);subprocess.run(cmd,cwd=output,env=ENV,stdout=log,stderr=subprocess.STDOUT,check=True)
    if own!=producer():raise RuntimeError('Harness producer changed while building')
    product=dict(schema='terluna.atmosphere.mcsnow-crystal-build/1',producer=own,source_identity=identity,
        source_files={str(p.relative_to(source)):sha(p) for p in sorted(source.rglob('*')) if p.is_file()},
        compiled_source_files=[str(p.relative_to(output)) for p in order],
        table_source_revision=old['revision'],table_identity=table_identity,physics_patches=[],
        compiler_path=compiler_path,compiler_sha256=sha(compiler_path),
        compiler_version=subprocess.check_output([compiler_path,'--version'],text=True).splitlines()[0],
        flags=FLAGS,commands=commands,executable_sha256=sha(output/'crystal'),
        driver_sha256=sha(output/DRIVER.name),build_log_sha256=sha(output/'build.log'),elapsed_s=time.monotonic()-start,
        precision='Native wp=selected_real_kind(12,307), mp=wp; double precision for mass, axes and shape.')
    (output/'build.json').write_text(json.dumps(product,indent=2,allow_nan=False)+'\n')
    return product

def initial(temperature_C,initialization,ice_density):
    if initialization=='fit':return lab.initial(temperature_C,ice_density)
    if initialization=='sphere':
        diameter=4e-6
        return dict(age_min=0.,mass_kg=math.pi/6*diameter**3*ice_density,
                    equatorial_dimension_m=diameter,polar_dimension_m=diameter,diameter_m=diameter,
                    density_kg_m3=ice_density,habit='Native dense sphere of 2 micrometre radius')
    raise ValueError('Unknown initialization')

def configuration(temperature_C,dt,growth_mode,initialization,*,ice_density):
    if not math.isfinite(dt) or dt<=0 or not math.isclose(60/dt,round(60/dt),rel_tol=0,abs_tol=1e-8):
        raise ValueError('Positive timestep must divide one minute')
    if growth_mode not in (1,2):raise ValueError('Unknown growth package')
    start=initial(temperature_C,initialization,ice_density)
    temperature=273.15+temperature_C
    return dict(temperature_K=temperature,pressure_Pa=scale.PRESSURE,
        vapor_pressure_Pa=float(scale.lunar.flakes.liquid_saturation_pressure(temperature)),
        dt_s=float(dt),duration_s=(lab.FITS[temperature_C][0]-start['age_min'])*60.,output_interval_s=60.,
        initial_mass_kg=start['mass_kg'],initial_equatorial_m=start['equatorial_dimension_m'],
        initial_polar_m=start['polar_dimension_m'],growth_mode=growth_mode,prolate_interpolation=True)

def namelist(config):
    def fmt(v):return '.true.' if v is True else '.false.' if v is False else format(v,'.17g')
    return '&crystal\n'+''.join('  '+k+' = '+fmt(v)+'\n' for k,v in config.items())+'/\n'

def verify_build(build_dir):
    build_dir=Path(build_dir);record=json.loads((build_dir/'build.json').read_text())
    if sha(build_dir/'crystal')!=record['executable_sha256']:raise ValueError('Executable changed')
    if sha(build_dir/DRIVER.name)!=record['driver_sha256']:raise ValueError('Native caller changed')
    for path,digest in record['source_files'].items():
        if sha(build_dir/'source'/path)!=digest:raise ValueError('Native source changed: '+path)
    for name,item in record['table_identity'].items():pinned.verify((build_dir/name).read_bytes(),item)
    return record

def read_csv(path,required):
    """Read a valid numeric prefix while retaining missing/truncated-output errors."""
    path=Path(path);errors=[];rows=[]
    if not path.exists():return rows,['Missing '+path.name]
    with path.open(newline='') as handle:
        reader=csv.reader(handle)
        names=next(reader,[])
        if len(set(names))!=len(names) or not set(required).issubset(names):
            return rows,['Invalid header in '+path.name]
        for line,values in enumerate(reader,start=2):
            if not values:continue
            if len(values)!=len(names):
                errors.append(f'Truncated or malformed {path.name} row {line}');break
            try:row={k:float(v) for k,v in zip(names,values)}
            except ValueError:
                errors.append(f'Nonnumeric {path.name} row {line}');break
            if not all(math.isfinite(v) for v in row.values()):
                errors.append(f'Nonfinite {path.name} row {line}');break
            rows.append(row)
    return rows,errors


def analyze_output(output,config,start,temperature_C,growth_mode,returncode):
    """Keep solver failure and output validity separate; never discard raw evidence."""
    output=Path(output)
    required=(*QUANTITIES,'elapsed_s','density_kg_m3','ice_volume_m3','aspect_ratio',
              'terminal_speed_m_s','ice_mass_balance_residual_kg')
    raw,errors=read_csv(output/'trajectory.csv',required)
    air_rows,air_errors=read_csv(output/'air.csv',('temperature_K','pressure_Pa','vapor_pressure_Pa'))
    errors.extend(air_errors)
    if len(air_rows)!=1:errors.append('Expected one prescribed-air state')
    air=air_rows[0] if air_rows else {}
    for key in ('temperature_K','pressure_Pa','vapor_pressure_Pa'):
        if key in air and not math.isclose(air[key],config[key],rel_tol=1e-13,abs_tol=1e-10):
            errors.append('Prescribed air differs: '+key)
    expected=np.arange(0,config['duration_s']+1,60.)
    rows=[]
    for index,row in enumerate(raw):
        if index>=len(expected) or row['elapsed_s']!=expected[index]:
            errors.append('Irregular trajectory timing');break
        if any(row[k]<=0 for k in (*QUANTITIES,'density_kg_m3','ice_volume_m3','aspect_ratio')):
            errors.append('Nonpositive crystal state');break
        if row['terminal_speed_m_s']<0:
            errors.append('Negative terminal speed');break
        row['age_min']=start['age_min']+row['elapsed_s']/60.
        if row['age_min']>=3:
            obs=lab.observed(temperature_C,row['age_min']);row['observed_fit']=obs
            row['ratio_to_fit']={k:row[k]/obs[k] for k in QUANTITIES}
            row['observed_fall_speed_m_s']=lab.fall_speed(temperature_C,row['age_min'])
            row['fall_speed_ratio_to_fit']=row['terminal_speed_m_s']/row['observed_fall_speed_m_s']
        rows.append(row)
    if not rows:errors.append('No valid crystal snapshots')
    else:
        for key in QUANTITIES:
            if not math.isclose(rows[0][key],start[key],rel_tol=1e-12):errors.append('Initial mapping differs: '+key)
    max_relative_ledger=max((abs(r['ice_mass_balance_residual_kg'])/r['mass_kg'] for r in rows),default=None)
    if max_relative_ledger is not None and max_relative_ledger>1e-10:errors.append('Ice growth accounting failed')
    for prev,row in zip(rows,rows[1:]):
        log_a=math.log(row['equatorial_dimension_m']/prev['equatorial_dimension_m'])
        row['interval_axis_log_growth_ratio']=math.log(row['polar_dimension_m']/prev['polar_dimension_m'])/log_a if log_a else None
    branches=[]
    if growth_mode==2:
        branches,branch_errors=read_csv(output/'branch_diagnostic.csv',('elapsed_s','branch_threshold_ratio','branch_criterion_met'))
        errors.extend(branch_errors)
    step_diagnostic=None
    if returncode==0:
        step_rows,step_errors=read_csv(output/'step_diagnostic.csv',('maximum_step_mass_fraction','steps'))
        errors.extend(step_errors)
        if len(step_rows)!=1:errors.append('Expected final step diagnostics')
        if step_rows:step_diagnostic=step_rows[0]
        if len(rows)!=len(expected):errors.append('Successful native exit with incomplete trajectory')
    complete=returncode==0 and len(rows)==len(expected) and not errors
    return dict(complete=complete,returncode=returncode,output_validation_errors=errors,
        snapshots=rows,air=air,maximum_relative_ice_deposition_residual=max_relative_ledger,
        branch_diagnostic=branches,step_diagnostic=step_diagnostic)


def run(build_dir,output,temperature_C,dt,growth_mode,initialization='fit'):
    build_dir=Path(build_dir).resolve();output=Path(output).resolve()
    if output.exists():raise FileExistsError(output)
    own=producer();build_record=verify_build(build_dir)
    if build_record['producer']!=own:raise ValueError('Build producer differs')
    constants=(build_dir/'source/src/mo_atmo_types.f90').read_text()
    density_match=re.search(r'\brhoi\s*=\s*([0-9.]+)_wp',constants,re.I)
    if not density_match:raise ValueError('Native ice density not found')
    ice_density=float(density_match[1]);start=initial(temperature_C,initialization,ice_density)
    config=configuration(temperature_C,dt,growth_mode,initialization,ice_density=ice_density)
    if sha(scale.lunar.CACHE/scale.SOURCE_PDF)!=lab.SOURCE_SHA:raise ValueError('Empirical source identity differs')
    output.mkdir(parents=True);(output/'crystal.nml').write_text(namelist(config))
    for name in TABLE_NAMES:shutil.copyfile(build_dir/name,output/name)
    inputs={p.name:sha(p) for p in output.iterdir() if p.is_file()}
    start_time=time.monotonic()
    with (output/'run.log').open('w') as log:
        completed=subprocess.run([str(build_dir/'crystal'),'crystal.nml'],cwd=output,env=ENV,stdout=log,stderr=subprocess.STDOUT)
    diagnostics=analyze_output(output,config,start,temperature_C,growth_mode,completed.returncode)
    if own!=producer() or verify_build(build_dir)!=build_record:raise RuntimeError('Producer/source identity changed during run')
    if inputs!={n:sha(output/n) for n in inputs}:raise RuntimeError('Input changed during execution')
    product=dict(schema='terluna.atmosphere.mcsnow-crystal-run/1',producer=own,build_record_sha256=sha(build_dir/'build.json'),
        **diagnostics,
        branch_diagnostic_scope='Recomputed revised-package 200 micrometre maximum-radius threshold; does not alter native growth. Original-package branching is not reconstructed.',
        build_identity=build_record,configuration=config,temperature_C=temperature_C,initialization=initialization,
        initial=start,growth_package='original IGF1' if growth_mode==1 else 'revised IGF2+ package',
        source_pdf_sha256=lab.SOURCE_SHA,inputs=inputs,
        outputs={p.name:sha(p) for p in output.iterdir() if p.is_file() and p.name not in inputs},
        elapsed_s=time.monotonic()-start_time,
        evidence='Unmodified pinned McSnow growth/drag/habit execution in fixed prescribed air; independent model comparison against fitted laboratory curves.',
        limits=['No nucleation, liquid collection, collisions, motion, vapour depletion or atmospheric heat feedback.',
                'Same imposed T,p,e as SCALE; native saturation, gravity and gas constants are retained and reported.',
                'Fit initialization compares processes after age 3 minutes; sphere initialization tests complete shape history from age zero.',
                'Native source differs from an identified paper executable; no exact published-figure reproduction is claimed.',
                'IGF2+ changes intrinsic table, small-crystal limiter and branching criterion together.'])
    (output/'run.json').write_text(json.dumps(product,indent=2,allow_nan=False)+'\n')
    return product

def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    b=sub.add_parser('build');b.add_argument('--output',type=Path,required=True);b.add_argument('--version',choices=('2024','2026'),required=True)
    b.add_argument('--public-source',type=Path,required=True);b.add_argument('--archive',type=Path)
    r=sub.add_parser('run');r.add_argument('--build',type=Path,required=True);r.add_argument('--output',type=Path,required=True)
    r.add_argument('--temperature',type=float,required=True);r.add_argument('--dt',type=float,default=1)
    r.add_argument('--growth-mode',type=int,choices=(1,2),required=True);r.add_argument('--initialization',choices=('fit','sphere'),default='fit')
    a=p.parse_args()
    if a.command=='build':result=build(a.output,a.version,a.public_source,a.archive);print(result['executable_sha256'])
    else:
        result=run(a.build,a.output,a.temperature,a.dt,a.growth_mode,a.initialization)
        print(json.dumps(dict(complete=result['complete'],returncode=result['returncode'],output_validation_errors=result['output_validation_errors'],last_snapshot=result['snapshots'][-1] if result['snapshots'] else None),indent=2))
if __name__=='__main__':main()
