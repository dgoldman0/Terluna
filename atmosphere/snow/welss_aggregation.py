"""Run and compare the declared-source McSnow author aggregation column.

This is an author-solver rerun, not an independently implemented collision model.
External executable, inputs, archive and raw output remain in scratch directories.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
import numpy as np
from scipy.io import netcdf_file
from atmosphere.snow import welss_aggregation_build as builder

HERE = Path(__file__).resolve().parent
MOMENTS = ('Nd','Md','Fn','Fm')
COLUMN = {name:col for start,suffix in [(2,''),(7,'_mm1'),(12,'_unr')] for col,name in enumerate((m+suffix for m in MOMENTS),start)}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def identity():
    return {p.name:sha(p) for p in [Path(__file__),HERE/'welss_aggregation_build.py',builder.MANIFEST]}


def configuration(aggregation, *, duration=36000., dtc=5., xi=1650):
    """Script1d_obl settings plus explicitly recorded inherited source defaults."""
    if duration <= 0 or dtc <= 0 or xi < 1 or duration % dtc != 0 or 600. % dtc != 0:
        raise ValueError('require positive integer multiplicity and time steps dividing output/duration')
    return dict(general_nml=dict(msg_level=6,out_type=1,tcase_id=111,ini_type=2,
        snowgeo_type=2,agggeo_type=1,rime_type=0,time_end=float(duration),dt_1dprof=600.,dt_2dprof=7200.,
        interpol_atmos=True,bnd_type=0,nucl_type=1,melt_type=1,shedding_type=1,diffusion_type=1,
        hydrobreakup_type=1,collbreakup_type=3,multiplication_type=0,mergeRatioLimit=0.,
        mergeNoTresh=1000000.,has_asol=False,perf_test=False),
        sp_nml=dict(iwc=1e-5,lwc0=0. if aggregation else .0002,nrp0=10000,xi0=int(xi),
        maxTotalNoSps=10000,dtc=float(dtc),idtgm=1,sigma0=.62,nclrate=16.5/xi,nclmass=4.8e-13,
        rm=1e-5,sp_vterm_id=3,agg_hist=-1,sp_habit=0,coll_kern_type=0,
        sp_kernel_id=2 if aggregation else 0,prolate_intp=True,sp_stick_effi=1,sp_frag=0,
        sp_kern_sig=0.,mu_gam=1.,nu_gam=1.),
        grid_nml=dict(dom_top=5000.,dom_bottom=0.,box_area=5.,nz=250),
        atmo_nml=dict(atmo_type=1,IGF_id=1,ssat=.05,h0=0.,h1=0.,h2=600.,nh1=2700.,nh2=3000.,
                      area_factor=1.1,atmofile='ecmwf_profile.txt'))


def namelist(config):
    def fmt(v):
        if isinstance(v,bool):return '.TRUE.' if v else '.FALSE.'
        if isinstance(v,str):return "'"+v.replace("'","''")+"'"
        return repr(v)
    return 'test_case_name = 1d_obl\n\n'+''.join('&'+group+'\n'+''.join('  '+k+' = '+fmt(v)+'\n' for k,v in values.items())+'/\n' for group,values in config.items())


def parse_profiles(path):
    times,frames,current=[],[],[]
    for line in Path(path).read_text().splitlines():
        if line.startswith('#"'):
            if current:frames.append(current);current=[]
            times.append(float(re.search(r'#"\s*([-+\d.]+)',line)[1])*60.)
        elif line.strip() and not line.startswith('#'):
            row=[float(v) for v in line.split()]
            if len(row)!=32:raise ValueError('profile row does not have32columns')
            current.append(row)
    if current:frames.append(current)
    a=np.asarray(frames,dtype=float);t=np.asarray(times)
    if a.ndim!=3 or a.shape[1:]!=(251,32) or len(a)!=len(t):raise ValueError('unexpected profile layout')
    if not np.isfinite(a).all() or np.any(a<0) or np.any(np.diff(t)<=0):raise ValueError('invalid profile values/times')
    if not np.all(a[:,:,0]==np.arange(251)*20.):raise ValueError('height bins differ from250-layer source')
    return t,a


def pathway_checks(run_dir, profiles):
    """Check saved process diagnostics and final particles, not unsaved stages."""
    run_dir=Path(run_dir)
    diag=np.loadtxt(run_dir/'hprofiles_diag.dat',comments='#')
    with (run_dir/'mass2fr.dat').open() as handle:
        next(handle)
        particles=np.loadtxt((line for line in handle if line.strip()),delimiter=',',ndmin=2)
    if diag.ndim!=2 or diag.shape[1]!=5 or particles.shape[1]!=15:
        raise ValueError('unexpected process/particle diagnostic layout')
    if not np.isfinite(diag).all() or not np.isfinite(particles).all():
        raise ValueError('nonfinite process or particle diagnostics')
    cols=[*range(17,21),*range(22,26),*range(27,31)]
    residual=particles[:,1]-particles[:,6]-particles[:,11]
    return dict(all_diagnostics_finite=True,
        maximum_saved_riming_rate_kg_s=float(np.max(np.abs(diag[:,3]))),
        maximum_saved_graupel_liquid_rimed_category_moment=float(np.max(np.abs(profiles[:,:,cols]))),
        final_superparticles=int(len(particles)),
        final_real_particle_number=float(np.sum(particles[:,13])),
        final_particle_mass_kg=float(np.sum(particles[:,13]*particles[:,1])),
        final_binned_mass_kg=float(np.sum(profiles[-1,:,3])*20.*5.),
        final_binned_mass_relative_rounding_difference=float(np.sum(profiles[-1,:,3])*100./np.sum(particles[:,13]*particles[:,1])-1.),
        final_maximum_rime_mass_kg=float(np.max(np.abs(particles[:,11]))),
        final_maximum_total_minus_pristine_minus_rime_mass_kg=float(np.max(np.abs(residual))),
        final_maximum_diameter_m=float(np.max(particles[:,4])),
        final_maximum_particle_mass_kg=float(np.max(particles[:,1])),
        final_maximum_monomer_count=float(np.max(particles[:,10])),
        scope='Saved600s process/category profiles and final particles; no assertion about every unsaved internal stage.')


def run(build_dir,output,*,aggregation,seed,duration=36000.,dtc=5.,xi=1650):
    build_dir,output=Path(build_dir).resolve(),Path(output).resolve()
    if output.exists():raise FileExistsError(output)
    if not -(2**31)<seed<2**31:raise ValueError('seed must fit Fortran default integer')
    own=identity();build_path=build_dir/'build.json';build=json.loads(build_path.read_text());source=Path(build['source_dir'])
    if build['producer']!={k:own[k] for k in ('welss_aggregation_build.py','welss_aggregation_sources.json')}:raise ValueError('builder identity differs')
    if sha(build_dir/'MCSNOW')!=build['executable_sha256']:raise ValueError('executable identity differs')
    manifest=json.loads(builder.MANIFEST.read_text())
    for item in manifest['files']:builder.verify((source/item['path']).read_bytes(),item)
    output.mkdir(parents=True)
    config=configuration(aggregation,duration=duration,dtc=dtc,xi=xi)
    (output/'input.nml').write_text(namelist(config))
    inputs={}
    for item in manifest['files']:
        if item['path'].startswith('input/'):
            p=source/item['path'];q=output/p.name;shutil.copyfile(p,q);inputs[q.name]=sha(q)
    inputs['input.nml']=sha(output/'input.nml')
    started=time.monotonic();env=os.environ.copy()
    env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    with (output/'run.log').open('w') as log:
        completed=subprocess.run([str(build_dir/'MCSNOW'),str(seed)],cwd=output,env=env,stdout=log,stderr=subprocess.STDOUT)
    elapsed=time.monotonic()-started
    if completed.returncode!=0:raise RuntimeError(f'McSnow exit{completed.returncode}; inspect{output}/run.log')
    t,a=parse_profiles(output/'hei2massdens.dat')
    if t[-1]!=duration:raise RuntimeError('author solver did not reach requested final output')
    if own!=identity() or sha(build_dir/'MCSNOW')!=build['executable_sha256']:raise RuntimeError('producer changed during run')
    if inputs!={name:sha(output/name) for name in inputs}:raise RuntimeError('run input changed')
    for item in manifest['files']:builder.verify((source/item['path']).read_bytes(),item)
    product=dict(schema='terluna.atmosphere.welss-author-run/1',producer=own,build=build,
        build_record_sha256=sha(build_path),configuration=config,seed=seed,elapsed_s=elapsed,
        aggregation=bool(aggregation),returncode=completed.returncode,inputs=inputs,
        outputs={p.name:sha(p) for p in sorted(output.iterdir()) if p.is_file() and p.name not in inputs},
        snapshot_times_s=t.tolist(),all_32_column_values_finite=True,
        pathway_checks=pathway_checks(output,a),
        evidence='Declared pinned author-model execution. This does not identify the2023 archive executable or independently validate aggregation physics.')
    (output/'run.json').write_text(json.dumps(product,indent=2,allow_nan=False)+'\n')
    return product


def read_archive(path):
    path=Path(path);manifest=json.loads(builder.MANIFEST.read_text());expected=manifest['archive']['members'][path.name]
    if sha(path)!=expected:raise ValueError('archive member checksum differs')
    with netcdf_file(path,mmap=False) as f:
        arrays={}
        for name in ['height',*COLUMN]:
            v=f.variables[name];a=v[:].astype(float)
            for attr in ('_FillValue','missing_value'):
                if hasattr(v,attr):a[a==getattr(v,attr)]=np.nan
            arrays[name]=a
        attrs={k:v.decode() if isinstance(v,bytes) else str(v) for k,v in f._attributes.items()}
    return arrays,dict(member=path.name,sha256=expected,attributes=attrs,**{k:v for k,v in manifest['archive'].items() if k!='members'})


def means(a, *, positive_only=False):
    values={}
    for name,col in COLUMN.items():
        x=a[:,:,col]
        if positive_only:
            count=(x>0).sum(axis=0);total=x.sum(axis=0)
            values[name]=np.divide(total,count,out=np.full(251,np.nan),where=count>0)
        else:values[name]=x.mean(axis=0)
    return values


def metrics(model,reference,z):
    valid=np.isfinite(model)&np.isfinite(reference)&(reference>0)&(z<3000.)
    result={}
    for label,mask in [('all_populated',valid),('below_nucleation',valid&(z<2700)),('nucleation_belt',valid&(z>=2700))]:
        if not np.any(mask):result[label]=None;continue
        d=model[mask]-reference[mask];r=reference[mask];rr=d/r;peak=np.argmax(np.abs(rr))
        result[label]=dict(count=int(mask.sum()),relative_L1=float(np.abs(d).sum()/r.sum()),relative_bias=float(d.sum()/r.sum()),
            maximum_absolute_relative_difference=float(np.max(np.abs(rr))),maximum_difference_height_m=float(z[mask][peak]))
    return result


def compare(run_dir,archive_path):
    run_dir=Path(run_dir);own=identity();run_path=run_dir/'run.json';record=json.loads(run_path.read_text())
    if record['producer']!=own:raise ValueError('run producer differs')
    for name,digest in {**record['inputs'],**record['outputs']}.items():
        if sha(run_dir/name)!=digest:raise ValueError('saved run file differs: '+name)
    t,a=parse_profiles(run_dir/'hei2massdens.dat');reference,archive=read_archive(archive_path);z=reference['height']
    if t[-1]!=36000 or not np.array_equal(t[t>=18000],np.arange(18000,36001,600)):raise ValueError('requires10h with31last5h snapshots')
    expected='ob_agg0.nc' if record['aggregation'] else 'ob_hp0.nc'
    if Path(archive_path).name!=expected:raise ValueError('archive case differs from run')
    modes={}
    for label,mask in [('last5h_31',t>=18000),('final4',t>=34200)]:
        raw=means(a[mask]);positive=means(a[mask],positive_only=True)
        modes[label]=dict(snapshot_times_s=t[mask].tolist(),snapshot_count=int(mask.sum()),
            raw_profiles={k:v.tolist() for k,v in raw.items()},
            positive_only_profiles={k:[None if not np.isfinite(x) else float(x) for x in v] for k,v in positive.items()},
            raw_comparison={k:metrics(v,reference[k],z) for k,v in raw.items()},
            positive_only_comparison={k:metrics(v,reference[k],z) for k,v in positive.items()})
    result=dict(schema='terluna.atmosphere.welss-author-comparison/1',producer=own,
        run_record_sha256=sha(run_path),run=record,reference=archive,height_m=z.tolist(),profiles=modes,
        evidence='Author-solver declared-configuration comparison to archived author population profiles; not independent solver validation or empirical validation.',
        reading_rule='Primary totals use31raw binned snapshots from18000through36000s. Four final snapshots and source-style zero-to-missing averages are explicit processing alternatives. No vertical smoothing.',
        caveats=['Pinned2024source may differ from2023archived executable.',
        'One random seed and prescribed discretization do not establish stochastic or numerical convergence.',
        'Totals are cell-average Nv and Nmv velocity moments, not boundary-crossing fluxes.',
        'Archive temporal window is supported but not uniquely identified by its31sample rounding lattice.',
        'Author output is rounded to three significant digits before arithmetic averaging.',
        'Fmono is omitted: pinned writer assigns a scalar to its full vector, so it is not used as a monomer-conservation diagnostic.'])
    if own!=identity() or sha(run_path)!=result['run_record_sha256'] or sha(archive_path)!=archive['sha256']:raise RuntimeError('comparison identity changed')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    r=sub.add_parser('run');r.add_argument('--build',type=Path,required=True);r.add_argument('--output',type=Path,required=True)
    r.add_argument('--aggregation',action='store_true');r.add_argument('--seed',type=int,default=20261021)
    r.add_argument('--duration',type=float,default=36000);r.add_argument('--dtc',type=float,default=5);r.add_argument('--xi',type=int,default=1650)
    c=sub.add_parser('compare');c.add_argument('--run',type=Path,required=True);c.add_argument('--archive',type=Path,required=True);c.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.command=='run':
        out=run(args.build,args.output,aggregation=args.aggregation,seed=args.seed,duration=args.duration,dtc=args.dtc,xi=args.xi)
        print(json.dumps({k:out[k] for k in ['aggregation','seed','elapsed_s','snapshot_times_s']},indent=2))
    else:
        out=compare(args.run,args.archive);args.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
        print(json.dumps({k:v['all_populated'] for k,v in out['profiles']['last5h_31']['raw_comparison'].items() if k in MOMENTS},indent=2))

if __name__=='__main__':main()
