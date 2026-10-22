"""Compare native SCALE snow against the complete prescribed lunar column.

The external-source adapter preserves SCALE's collision, drag and growth laws.
It supplies the original cloud profile, seed number flux and local water ceiling,
with explicit prescribed-air and evolving-reservoir experiments. Neither solver
is an empirical lunar snow forecast. Cached upstream files remain byte identified.
"""
from pathlib import Path
import argparse
import difflib
import hashlib
import json
import math
import os
import shutil
import signal
import subprocess
import time
import numpy as np
from atmosphere.snow import scale_column as base, scale_aggregation as boxes
from atmosphere.snow import scale_build, scale_collision_observer as observer
from atmosphere.snow import event_column as custom, flakes
from shared.provenance import constants_used

ROOT=base.REPO
CACHE=base.DEFAULT_CACHE
HOOK=Path(__file__).with_name('scale_lunar_hook.f90')
SOURCE='contrib/SDM/scale_atmos_phy_mp_sdm.F90'
SOURCE_SHA=scale_build.PATCH_IDENTITIES[SOURCE][1]
MODES={'prescribed':1,'thermal':2,'reservoir':3}
ARRIVAL_FIELDS=('time_s','horizontal_cell','multiplicity','mass_kg','dimension_m','density_kg_m3','aspect','monomers','phase')
LEDGER_FIELDS=('time_s','injected_kg_m2','deposited_kg_m2','vapor_source_kg_m2','precipitation_kg_m2','suspended_kg_m2','vapor_kg_m2','residual_kg_m2','max_absolute_residual_kg_m2','minimum_temperature_K','maximum_temperature_K','maximum_dimension_m','active_superparticles','peak_superparticles','removed_sensible_heat_J_m2')


def patch(text):
    if hashlib.sha256(text.encode()).hexdigest()!=SOURCE_SHA:
        raise ValueError('Unidentified SCALE source')
    text=text.replace('if( buffact < 1.0_RP .or. buffact > 1.0_RP ) then', 'if( abs(buffact-1.0_RP)>1e-12_RP ) then')
    state,procedures=HOOK.read_text().split('! PROCEDURES\n')
    state=state.split('! STATE\n')[1]
    text=text.replace('  subroutine ATMOS_PHY_MP_sdm_setup( MP_TYPE )', '  subroutine ATMOS_PHY_MP_sdm_setup( MP_TYPE )',1)
    # Module data before CONTAINS; helper procedures use the native module types.
    marker='contains\n'
    if text.count(marker)!=1: raise ValueError('Unexpected module structure')
    text=text.replace(marker,state+'\n'+marker+procedures+'\n',1)
    start=text.index('  subroutine sdm_calc(');end=text.index('  end subroutine sdm_calc',start)
    chunk=text[start:end]
    changes=[
      ('   real(RP) :: tmp_mink','   real(RP) :: tmp_mink\n   logical :: sl_active'),
      ('      ! Calculate super-droplets process.', '''      call sl_start(DENS,RHOT,QTRC,sd_num,sd_n,sd_liqice,sd_x,sd_y,sd_rk,sd_r,sd_asl,sdi,sl_active)
      if(.not.sl_active)return
      if(sl_mode>0)then
        sdm_dtevl=sl_h;sdm_dtcol=sl_h;sdm_dtadv=sl_h;sdm_dtmlt=sl_h;sdm_dtsbl=sl_h
      endif
      ! Calculate super-droplets process.'''),
      ('            call sdm_subldep(            &','            call sl_before_growth(sdi)\n            call sdm_subldep(            &'),
      ('            ! get density of solid-water after sublimation/deposition','            call sl_limit_growth(DENS,QTRC,sd_num,sd_n,sd_liqice,sd_ri,sd_rj,sd_rk,sdi,sdm_dtsbl)\n\n            ! get density of solid-water after sublimation/deposition'),
      ('            call sdm_subldep_updatefluid(RHOT,QTRC,DENS,crs_val1p,crs_val1c)','            if(sl_mode/=1) call sdm_subldep_updatefluid(RHOT,QTRC,DENS,crs_val1p,crs_val1c)'),
      ('      call sdm_sd2prec(dt,','      call sl_finish(DENS,RHOT,QTRC,sd_num,sd_n,sd_liqice,sd_r,sd_rk,sd_x,sd_y,sdi)\n      call sdm_sd2prec(dt,'),
    ]
    for old,new in changes:
        if chunk.count(old)!=1: raise ValueError('Ambiguous insertion: '+old)
        chunk=chunk.replace(old,new)
    return text[:start]+chunk+text[end:]


def collision_patch(text,collection='native'):
    result=observer.patch(text)
    if collection=='inhouse':
        old='             c_rate(tc) = 0.1_RP'
        if result.count(old)!=1: raise ValueError('Unexpected ice collection constant')
        new='             k=floor(sd_rk(icptc))+1; i=floor(sd_ri(icptc))+1; j=floor(sd_rj(icptc))+1\n             c_rate(tc) = exp(0.09_RP*(t_scale(k,i,j)-t0))'
        result=result.replace(old,new)
    elif collection!='native':raise ValueError('Unknown collection law')
    return result


def build(cache=CACHE,collection='native'):
    if collection not in ('native','inhouse'):raise ValueError('Unknown collection law')
    cache=Path(cache).resolve(); target=cache/('lunar' if collection=='native' else 'lunar_sticking')
    target.mkdir(exist_ok=False)
    work=target/'work'
    shutil.copytree(cache/'work',work,symlinks=True)
    patches=[]
    for rel,fn in [(SOURCE,patch),(observer.SOURCE_REL,lambda text:collision_patch(text,collection))]:
        p=work/rel;old=p.read_text();new=fn(old);p.write_text(new)
        diff=''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='a/'+rel,tofile='b/'+rel))
        dest=target/(p.name+'.patch');dest.write_text(diff)
        patches.append(dict(source=rel,original_sha256=hashlib.sha256(old.encode()).hexdigest(),patched_sha256=base.digest(p),patch_sha256=base.digest(dest)))
    deps=cache/'deps/usr'
    result=scale_build.build(work,'-I'+str(deps/'include'),str(deps/'lib/x86_64-linux-gnu/libnetcdff.so.7')+' /usr/lib/x86_64-linux-gnu/libnetcdf.so.19 -Wl,-rpath,'+str(deps/'lib/x86_64-linux-gnu')+' -lm -lz',workers=1)
    record=dict(collection=collection,patches=patches,build=result,hook_sha256=base.digest(HOOK),adapter_sha256=base.digest(Path(__file__)))
    boxes.save(target/'manifest.json',record)
    return record


def specification(layers=30):
    world,binding=flakes.moon_column()
    original=json.loads((ROOT/'atmosphere/snow/results/event_column.json').read_text())['groups']['moon_v4']['settings']
    depth=original['depth_m']; dz=depth/layers
    states=[world.state(float(world.z_freeze+(k+.5)*dz)) for k in range(layers)]
    return world,binding,original,states


def configure(name,mode='prescribed',packets=2,dt=60.,seed=1,layers=30,duration=None,collisions=True,cache=CACHE,capacity_per_layer=64,dx=1.,collection='native'):
    if mode not in MODES or packets<1 or int(packets)!=packets or dt not in (5.,10.,15.,30.,60.) or seed<1 or layers<1 or dx<=0:
        raise ValueError('Invalid lunar configuration')
    world,binding,original,states=specification(layers)
    depth=original['depth_m'];dz=depth/layers
    duration=original['total_time_s'] if duration is None else float(duration)
    if not math.isfinite(duration) or duration<=0: raise ValueError('Invalid duration')
    cache=Path(cache).resolve()
    if collection not in ('native','inhouse'):raise ValueError('Unknown collection law')
    build_directory='lunar' if collection=='native' else 'lunar_sticking'
    manifest=json.loads((cache/build_directory/'manifest.json').read_text())
    if manifest['hook_sha256']!=base.digest(HOOK): raise ValueError('Stale forcing executable')
    d=base.configure(cache,name,lunar=True,duration=60.,dt=dt,sd=capacity_per_layer,seed=seed)
    # Top guard cell accommodates the native restriction on its particle domain.
    for file in ('init.conf','run.conf'):
        text=(d/file).read_text()
        changes=dict(KMAX=layers+1,DX=dx,DY=dx,DZ=dz,BUFFER_DZ=dz)
        if file=='run.conf':
            changes.update(TIME_DURATION=math.ceil(duration/dt)*dt,TIME_DT=dt,TIME_DT_ATMOS_PHY_MP=dt,
              docondensation=False,domovement=True,domeltfreeze=(mode=='thermal'),dosublimation=True,doautoconversion=collisions,
              sdm_dtcmph=(dt,dt,dt,dt,dt),sdm_sdnmlvol=dx*dx*dz,sdm_zupper=depth,
              sdm_extbuf=0,sdm_dmpvar=0,HISTORY_DEFAULT_TINTERVAL=3600.,MONITOR_STEP_INTERVAL=1000000)
            text=text.replace('&PARAM_ATMOS_PHY_MP_SDM\n',"&PARAM_ATMOS_PHY_MP_SDM\n SD_IN_BASENAME = 'particles',\n")
        for k,v in changes.items():text=boxes.set_value(text,k,v)
        (d/file).write_text(text)
    # Native initializer needs a valid sounding even though the scenario hook
    # subsequently assigns the exact supplied cell T,p,e at its first step.
    rd=base.SOURCE_CONSTANTS['Rdry'];rv=base.SOURCE_CONSTANTS['Rvap'];cp=base.SOURCE_CONSTANTS['CPdry']
    def line(z,first=False):
        s=world.state(float(world.z_freeze+z));q=(rd/rv)*s['e']/(s['p']-(1-rd/rv)*s['e']);theta=s['t']*(1e5/s['p'])**(rd/cp)
        return f"{s['p']/100 if first else z:.15g} {theta:.15g} {q*1000:.15g}"+('' if first else ' 0 0')+'\n'
    (d/'env.txt').write_text(line(0,True)+''.join(line(float(z)) for z in np.linspace(dz/10,depth+4*dz,400)))
    n=8*layers*capacity_per_layer
    zeros=np.zeros(n)
    pop=dict(n=np.zeros(n,dtype=np.int64),rk=np.full(n,-999.999),x=zeros,y=zeros,z=np.full(n,-999.999),r=zeros,re=np.full(n,1e-4),rp=np.full(n,1e-4),rho=np.full(n,300.),phase=np.full(n,10,dtype=np.int16),monomers=np.ones(n,dtype=np.int32))
    boxes.write_restart(d/'particles.pe000000',pop)
    forcing=f'{MODES[mode]} {layers} {packets} {seed}\n'
    forcing+=' '.join(format(x,'.17g') for x in (depth,dx,original['injected_number_flux_m2_s'],original['snowfall_mm_h']/3600/depth,duration,duration/2))+'\n'
    forcing+=' '.join(format(x,'.17g') for x in (custom.GEOMETRY.mass_coefficient,custom.GEOMETRY.mass_exponent,*custom.ag.START_M))+'\n'
    forcing+=''.join(' '.join(format(s[k],'.17g') for k in ('t','p','e'))+'\n' for s in states)
    (d/'lunar_forcing.txt').write_text(forcing)
    settings=dict(name=name,mode=mode,packets_per_step_per_horizontal_cell=packets,dt_s=dt,seed=seed,layers=layers,
      duration_s=duration,gather_start_s=duration/2,depth_m=depth,dx_m=dx,horizontal_cells=8,capacity=n,collisions=collisions,
      original_settings=original,profile=[{k:float(v) for k,v in s.items()} for s in states],
      forcing='Uniform local 3 mm/h column ceiling, original seed flux and mass/dimension distribution, initially empty ice column. Prescribed mode shares the ceiling across native trial growth increments; reservoir modes add vapor up to that ceiling and liquid saturation.',
      unmatched='Native evolving porous spheroid geometry, Bohm drag and constant ice efficiency 0.1; spherical injected shapes with matched mass and maximum dimension. Ballistic settling; no chipping. Native dry-air thermodynamic constants. 8 columns of weighted particles, not the custom finite individual-particle volume.',
      build=manifest,build_directory=build_directory,collection=collection,producer=custom.identity(binding),constants=constants_used(['atmosphere/snow/scale_column.py']),
      adapter_sha256=base.digest(Path(__file__)))
    boxes.save(d/'settings.json',settings)
    boxes.save(d/'inputs.json',{p.name:base.digest(p) for p in d.iterdir() if p.is_file()})
    return d


def execute(d,timeout=3600):
    d=Path(d).resolve();cache=d.parents[1]
    settings=json.loads((d/'settings.json').read_text());manifest=settings['build']
    env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','GFORTRAN_UNBUFFERED_ALL':'y',
         'LD_LIBRARY_PATH':str(cache/'deps/usr/lib/x86_64-linux-gnu')+os.pathsep+os.environ.get('LD_LIBRARY_PATH','')}
    rows=[]
    for exe,conf in [('scale-les_init','init.conf'),('scale-les','run.conf')]:
        binary=cache/settings.get('build_directory','lunar')/'work'/scale_build.CASE_REL/exe
        if base.digest(binary)!=manifest['build']['executable_sha256'][exe]:raise ValueError('Binary identity mismatch')
        start=time.monotonic()
        with (d/(exe+'.log')).open('x') as log:
            p=subprocess.Popen([str(binary),conf],cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            try:rc=p.wait(timeout=timeout)
            except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait();rc='timeout'
        rows.append(dict(executable=exe,exit_code=rc,elapsed_s=time.monotonic()-start,sha256=base.digest(binary)))
        boxes.save(d/'execution.json',rows)
        print(d.name,exe,rc,round(rows[-1]['elapsed_s'],2),flush=True)
        if rc!=0:raise RuntimeError(f'Failed native run; retained {d}')
    return rows


def analyze(d):
    d=Path(d);s=json.loads((d/'settings.json').read_text());executions=json.loads((d/'execution.json').read_text())
    if len(executions)!=2 or any(r['exit_code']!=0 for r in executions):raise ValueError('Incomplete run')
    ledger=np.atleast_2d(np.loadtxt(d/'lunar_ledger.txt'))
    if not np.isclose(ledger[-1,0],s['duration_s'],rtol=0,atol=1e-6):raise ValueError('Incomplete physical duration')
    rows=np.fromfile(d/'lunar_arrivals.bin',dtype='>f8')
    if rows.size%len(ARRIVAL_FIELDS):raise ValueError('Incomplete arrival record')
    rows=rows.reshape(-1,len(ARRIVAL_FIELDS))
    kept=rows[rows[:,0]>=s['gather_start_s']]
    area=s['horizontal_cells']*s['dx_m']**2
    diam=custom.diameter_queries()
    def summary(a,footprint=area,period=s['duration_s']-s['gather_start_s']):
        w=a[:,2]*a[:,3]
        return dict(packets=len(a),number_m2=float(a[:,2].sum()/footprint),mass_kg_m2=float(w.sum()/footprint),
          precipitation_mm_h=float(w.sum()/footprint/period*3600),
          snowfall_mm_h=float(w[a[:,8]==10].sum()/footprint/period*3600),
          dimension_mass_quantiles_m=base.weighted_quantile(a[:,4],w),
          ice_dimension_mass_quantiles_m=base.weighted_quantile(a[a[:,8]==10,4],w[a[:,8]==10]),
          mass_quantiles_kg=base.weighted_quantile(a[:,3],w),
          inherited_geometry_dimension_mass_quantiles_m=base.weighted_quantile((a[:,3]/custom.GEOMETRY.mass_coefficient)**(1/custom.GEOMETRY.mass_exponent),w),
          density_mass_quantiles_kg_m3=base.weighted_quantile(a[:,5],w),aspect_mass_quantiles=base.weighted_quantile(a[:,6],w),
          aggregate_mass_fraction=float(w[a[:,7]>1].sum()/w.sum()) if w.sum() else None,
          ice_mass_fraction=float(w[a[:,8]==10].sum()/w.sum()) if w.sum() else None,
          cdf=[float(w[a[:,4]<=x].sum()/w.sum()) if w.sum() else None for x in diam],
          mass_fraction_above={str(x):float(w[a[:,4]>x].sum()/w.sum()) if w.sum() else None for x in (.01,.15,.6)})
    counters=np.atleast_2d(np.loadtxt(d/'collision_observer.txt')) if (d/'collision_observer.txt').exists() else np.zeros((0,8))
    final=dict(zip(LEDGER_FIELDS,map(float,ledger[-1])))
    initial_vapor=sum(v['e']/base.SOURCE_CONSTANTS['Rvap']/v['t'] for v in s['profile'])*s['depth_m']/s['layers']
    result=dict(schema='terluna.atmosphere.scale-lunar-case/1',analysis_sha256=base.digest(Path(__file__)),initial_vapor_kg_m2=initial_vapor,settings=s,execution=executions,
       ledger_fields=LEDGER_FIELDS,ledger=ledger.tolist(),final=final,
       max_relative_water_residual=final['max_absolute_residual_kg_m2']/(final['injected_kg_m2']+initial_vapor+(final['deposited_kg_m2'] if s['mode']=='prescribed' else final['vapor_source_kg_m2'])),
       arrival_fields=ARRIVAL_FIELDS,arrivals=summary(kept),query_diameter_m=diam.tolist(),
       horizontal_columns=[summary(kept[kept[:,1]==j],footprint=area/8) for j in range(8)],
       gathering_time_blocks=[summary(kept[(kept[:,0]>=lo)&(kept[:,0]<hi)],period=hi-lo) for lo,hi in zip(np.linspace(s['gather_start_s'],s['duration_s']+1e-6,5)[:-1],np.linspace(s['gather_start_s'],s['duration_s']+1e-6,5)[1:])],
       observer=dict(calls=len(counters),attempted=float(counters[:,2].sum()),omitted=float(counters[:,4].sum()),capped_pairs=int(counters[:,5].sum()),max_probability=float(counters[:,6].max()) if len(counters) else None),
       native_precipitation_vs_arrivals_absolute_error_kg_m2=None,
       files={p.name:base.digest(p) for p in d.iterdir() if p.is_file() and p.name!='analysis.json'},
       evidence='Numerical comparison under declared artificial forcing, not measured snow or a mechanically validated upper size. Prescribed runs change several native closures together; reservoir and thermal contrasts isolate atmospheric-feedback choices within SCALE.')
    boxes.save(d/'analysis.json',result)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=['build','run','analyze']);p.add_argument('--name');p.add_argument('--mode',choices=MODES,default='prescribed')
    p.add_argument('--packets',type=int,default=2);p.add_argument('--dt',type=float,default=60);p.add_argument('--seed',type=int,default=1)
    p.add_argument('--layers',type=int,default=30);p.add_argument('--duration',type=float);p.add_argument('--capacity',type=int,default=64)
    p.add_argument('--no-collisions',action='store_true');p.add_argument('--dx',type=float,default=1.)
    a=p.parse_args()
    if a.action=='build':print(json.dumps(build(),indent=2))
    elif a.action=='run':
        d=configure(a.name,a.mode,a.packets,a.dt,a.seed,a.layers,a.duration,not a.no_collisions,capacity_per_layer=a.capacity,dx=a.dx)
        execute(d);r=analyze(d);print(json.dumps({k:r[k] for k in ('final','arrivals','observer')},indent=2))
    else:
        r=analyze(CACHE/'runs'/a.name);print(json.dumps({k:r[k] for k in ('final','arrivals','observer')},indent=2))
