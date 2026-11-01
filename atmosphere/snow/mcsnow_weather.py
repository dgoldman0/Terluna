"""Native McSnow diagnostic crystals under resolved archived lunar vapour.

Crystals are dilute probes of imposed weather. They do not consume the archived
field or define an injected population. Source files and earlier products stay
unchanged; the matching saturation experiment is a humidity sensitivity.
"""
from __future__ import annotations

import argparse
import contextlib
import ctypes
import difflib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np

from atmosphere.snow import mcsnow_lunar_growth as lunar
from atmosphere.snow import moist_history as history
from shared.constants import CM1_DRY_AIR_GAS_CONSTANT, CM1_VAPOUR_GAS_CONSTANT

ROOT=lunar.ROOT
CACHE=ROOT/'research/runs/snow_scale/mcsnow_weather'
DRIVER=ROOT/'atmosphere/snow/mcsnow_weather_driver.f90'
PRODUCT=ROOT/'atmosphere/snow/results/mcsnow_weather.json'
native=lunar.native
OUTPUT_FIELDS=('diameter_m','equatorial_dimension_m','polar_dimension_m','density_kg_m3',
    'terminal_speed_m_s','deposited_mass_kg','latent_heat_J','peak_diameter_m','advanced_s',
    'substeps','saturation_water_Pa','saturation_ice_Pa','vapour_pressure_Pa','mass_residual_kg')
STATE_FIELDS=('mass_kg','ice_volume_m3','aspect_ratio','initial_temperature_K','initial_ice_supersaturation')
_INITIALIZED=set()


def save(path, value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def build(world='moon'):
    path=CACHE/('build_'+world)
    if path.exists():return verify_build(world)
    original,changed=lunar.patches('plate_base',world)
    baseline=native.verify_build(lunar.model.BASELINE)
    path.mkdir(parents=True)
    shutil.copytree(lunar.model.BASELINE/'source',path/'source')
    diff=''
    for name,body in changed.items():
        (path/'source/src'/name).write_text(body)
        diff+=''.join(difflib.unified_diff(original[name].splitlines(True),body.splitlines(True),
                                         fromfile='original/'+name,tofile=world+'/'+name))
    (path/'weather.patch').write_text(diff)
    shutil.copyfile(DRIVER,path/DRIVER.name)
    for name in native.TABLE_NAMES:shutil.copyfile(lunar.model.BASELINE/name,path/name)
    compiler=baseline['compiler_path'];flags=[*native.FLAGS,'-fPIC','-ffree-line-length-none']
    sources=native.selected_sources(path/'source/src',path/DRIVER.name,compiler)
    objects=[];commands=[]
    with (path/'build.log').open('w') as log:
        for i,source in enumerate(sources):
            obj=path/f'{i:02d}_{source.stem}.o';objects.append(str(obj))
            cmd=[compiler,*flags,'-I'+str(path/'source/src'),'-I'+str(path),'-J'+str(path),'-c',str(source),'-o',str(obj)]
            commands.append(cmd);subprocess.run(cmd,cwd=path,env=native.ENV,stdout=log,stderr=subprocess.STDOUT,check=True)
        cmd=[compiler,*flags,'-shared',*objects,'-o',str(path/'libweather.so')]
        commands.append(cmd);subprocess.run(cmd,cwd=path,env=native.ENV,stdout=log,stderr=subprocess.STDOUT,check=True)
    record=dict(schema='terluna.atmosphere.mcsnow-weather-build/1',world=world,
        baseline_record_sha256=native.sha(lunar.model.BASELINE/'build.json'),
        source_identity=baseline['source_identity'],table_identity=baseline['table_identity'],
        files={str(p.relative_to(path)):native.sha(p) for p in sorted(path.rglob('*'))
               if p.is_file() and p.suffix not in ('.o','.mod')},
        compiler_path=compiler,compiler_sha256=native.sha(compiler),compiler_version=baseline['compiler_version'],
        flags=flags,commands=commands)
    save(path/'build.json',record)
    return record


def verify_build(world):
    path=CACHE/('build_'+world);record=json.loads((path/'build.json').read_text())
    native.verify_build(lunar.model.BASELINE)
    if native.sha(lunar.model.BASELINE/'build.json')!=record['baseline_record_sha256']:raise ValueError('Baseline changed')
    for name,digest in record['files'].items():
        if native.sha(path/name)!=digest:raise ValueError('Build artifact changed: '+name)
    if native.sha(DRIVER)!=record['files'][DRIVER.name]:raise ValueError('Native interface changed')
    _,changed=lunar.patches('plate_base',world)
    if any((path/'source/src'/name).read_text()!=body for name,body in changed.items()):raise ValueError('Candidate differs')
    return record


class Library:
    """Serial batch calls; native module globals must not be used concurrently."""
    def __init__(self, world='moon'):
        self.identity=verify_build(world);self.world=world;path=CACHE/('build_'+world)
        self.lib=ctypes.CDLL(str(path/'libweather.so'))
        fp=np.ctypeslib.ndpointer(dtype=np.float64,flags='F_CONTIGUOUS')
        ip=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')
        self.lib.weather_evolve.argtypes=[ctypes.c_int,fp,fp,fp,ctypes.c_int,
            ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,fp,ip]
        self.lib.weather_evolve.restype=None
        self.lib.weather_constants.argtypes=[fp];self.lib.weather_constants.restype=None
        self.lib.weather_init.argtypes=[];self.lib.weather_init.restype=None
        old=Path.cwd()
        try:
            if str(path) not in _INITIALIZED:
                os.chdir(path);self.lib.weather_init();_INITIALIZED.add(str(path))
        finally:os.chdir(old)
        values=np.empty(4);self.lib.weather_constants(values)
        self.ice_density,self.Rd,self.Rv,self.gravity=values

    def seed(self, forcing, diameter):
        forcing=np.asarray(forcing)
        volume=math.pi/6*diameter**3
        state=np.empty((5,forcing.shape[1]),order='F')
        state[:]=np.array([volume*self.ice_density,volume,1.,0.,0.])[:,None]
        state[3]=forcing[0]
        diag,code=self.advance(state,forcing,0.)
        state[4]=forcing[2]/diag[11]-1
        return state

    def advance(self, state, forcing, duration, *, saturated=False,max_step=.5,max_fraction=.005,
                cutoff_diameter=2e-7,diameter_stop=1.):
        state=np.asarray(state)
        if state.ndim!=2 or state.shape[0]!=5 or not state.flags.f_contiguous:raise ValueError('Expected Fortran-contiguous five-field state')
        force=np.asfortranarray(forcing,dtype=float)
        if force.shape!=(3,state.shape[1]) or not np.all(np.isfinite(force)) or np.any(force<=0):raise ValueError('Invalid forcing')
        dt=np.array(np.broadcast_to(duration,(state.shape[1],)),dtype=float)
        if not np.all(np.isfinite(state)) or np.any(state[:3]<=0) or np.any(dt<0) or not np.all(np.isfinite(dt)):
            raise ValueError('Invalid state/time')
        if not (max_step>0 and 0<max_fraction<=1 and cutoff_diameter>0 and diameter_stop>0):
            raise ValueError('Invalid integration control')
        out=np.empty((len(OUTPUT_FIELDS),state.shape[1]),order='F');code=np.empty(state.shape[1],np.int32)
        if saturated==2:
            # Resolve an explicitly labelled instantaneous-condensation
            # sensitivity, without silently clipping the archived experiment.
            self.lib.weather_evolve(state.shape[1],state,force,np.zeros_like(dt),0,max_step,max_fraction,
                self.ice_density*math.pi/6*cutoff_diameter**3,diameter_stop,out,code)
            force=force.copy(order='F');force[2]=np.minimum(force[2],out[10])
        self.lib.weather_evolve(state.shape[1],state,force,dt,int(saturated),max_step,max_fraction,
            self.ice_density*math.pi/6*cutoff_diameter**3,diameter_stop,out,code)
        if np.any(code<0) or not np.all(np.isfinite(out)):raise ValueError('Native integration failed: '+str(code.tolist()))
        return out,code


def weather_air(raw):
    """Return native inputs from CM1 theta, total pressure and dry-air mixing ratio."""
    raw=np.asarray(raw)
    T=raw[:,0]*(raw[:,1]/history.t.CM1_P00)**(CM1_DRY_AIR_GAS_CONSTANT/history.t.CM1_CP)
    pv=raw[:,1]*raw[:,2]/(CM1_DRY_AIR_GAS_CONSTANT/CM1_VAPOUR_GAS_CONSTANT+raw[:,2])
    return np.asfortranarray(np.array([T,raw[:,1],pv]))


def prepare_weather():
    folder=CACHE/'weather';target=folder/'fields.npy'
    prior=json.loads(history.OUT.read_text());item=history.inventory(history.CASE)
    intended=next(r for r in prior['inventory'] if r['case']==history.CASE)
    if item['control_hashes']!=intended['controls']:raise ValueError('Archived controls changed')
    if target.exists():return verify_weather()
    folder.mkdir(parents=True)
    data,sources,diagnostics=history.load_window(item)
    if sources!=prior['source']['frames']:raise ValueError('Selected archived fields differ from prior identity')
    np.save(target,data,allow_pickle=False)
    record=dict(schema='terluna.atmosphere.mcsnow-weather-input/1',prior_history_sha256=native.sha(history.OUT),
        inventory=history.inventory_record(item),selected_levels=item['levels'].tolist(),
        z_m=item['z'][item['levels']].tolist(),dx_m=item['dx'],frames=sources,
        field_order=history.t.FIELDS,fields_sha256=native.sha(target),shape=list(data.shape),
        release_rule='Use the 354 existing distinct release cells; replace fixed drifts by native evolving settling.',
        releases=[r[:4] for r in prior['configuration']['seeds'] if r[4]==0.],
        diagnostics=diagnostics)
    save(folder/'input.json',record)
    return record


def verify_weather():
    folder=CACHE/'weather';record=json.loads((folder/'input.json').read_text())
    if native.sha(history.OUT)!=record['prior_history_sha256']:raise ValueError('Prior history changed')
    if native.sha(folder/'fields.npy')!=record['fields_sha256']:raise ValueError('Cached fields changed')
    current=history.inventory(history.CASE)
    if current['control_hashes']!=record['inventory']['controls']:raise ValueError('Archived controls changed')
    for entry in record['frames']:
        p=history.RUNS/history.CASE/f"cm1out_t{entry['number']:06d}_s.dat"
        st=p.stat()
        if st.st_size!=entry['bytes'] or st.st_mtime_ns!=entry['mtime_ns']:raise ValueError('Archive changed after caching')
    return record


STATUSES={0:'active',1:'small_ice_residual',2:'diameter_censored',3:'warm_exit',
          4:'cold_habit_range_exit',5:'lower_height_exit',6:'upper_height_exit',7:'archive_end'}


def variable_sample(data,positions,z,dx,times,stride):
    """Sample each particle at its own partial-step time without extrapolation."""
    result=np.empty((len(positions),len(history.t.FIELDS)))
    for seconds in np.unique(times):
        mask=times==seconds
        result[mask]=history.sample(data,positions[mask],z,dx,float(seconds),stride)
    return result


def thermal_fraction(before,after,low=243.15,high=273.15):
    return history.boundary_fraction(np.asarray(before),np.asarray(after),low,high)


def integrate(lib,data,record,*,diameter=4e-6,step_s=20.,stride=1,saturated=False,cutoff=2e-7,
              growth_step=.5,growth_fraction=.005):
    """Midpoint trajectory with adaptive native growth in midpoint air.

    Boundary and temperature crossings use linear segment fractions. Native
    sublimation/size exits can shorten a segment further. Prescribed weather
    already contains its original model's microphysical and transport response;
    these unweighted probes add no second population to its water budget.
    """
    start_wall=time.monotonic();z=np.asarray(record['z_m']);dx=record['dx_m']
    seeds=np.asarray(record['releases']);n=len(seeds);end=(len(data)-1)*600.
    if end%(600*stride) or 600%step_s or np.any(seeds[:,0]%step_s):raise ValueError('Time grids must align')
    position=seeds[:,1:].copy();initial_air=weather_air(variable_sample(data,position,z,dx,seeds[:,0],1))
    state=lib.seed(initial_air,diameter);initial_mass=state[0].copy()
    if saturated:
        initial_diag,_=lib.advance(state,initial_air,0.,saturated=saturated)
        state[4]=initial_diag[12]/initial_diag[11]-1
    status=np.zeros(n,np.int32);ages=np.zeros(n);peak=np.full(n,diameter)
    deposited=np.zeros(n);sublimated=np.zeros(n);heat=np.zeros(n);moist=np.zeros(n)
    max_balance=0.;total_substeps=0;rows=[];final_diag=np.zeros((14,n))
    row_fields=('id','window_s','age_s','x_m','y_m','height_m',*STATE_FIELDS[:3],
                *OUTPUT_FIELDS[:5],'temperature_K','pressure_Pa','vapour_pressure_Pa','ice_rh','status')
    def snapshot(ids,times,forcing,diag):
        for j,i in enumerate(ids):
            rows.append([float(i),float(times[j]),float(ages[i]),*position[i].tolist(),
                *state[:3,i].tolist(),*diag[:5,j].tolist(),*forcing[:2,j].tolist(),
                float(diag[12,j]),float(diag[12,j]/diag[11,j]),float(status[i])])
    for it in range(int(round(end/step_s))):
        seconds=it*step_s
        ids=np.flatnonzero((status==0)&(seeds[:,0]<=seconds))
        if not len(ids):continue
        old=position[ids].copy();s0=np.asfortranarray(state[:,ids]);oldraw=history.sample(data,old,z,dx,seconds,stride)
        f0=weather_air(oldraw);diag0,_=lib.advance(s0,f0,0.,saturated=saturated)
        born=seeds[ids,0]==seconds
        if born.any():snapshot(ids[born],np.full(born.sum(),seconds),f0[:,born],diag0[:,born])
        invalid=(f0[0]>=273.15)|(f0[0]<243.15)
        if invalid.any():
            bad=ids[invalid];status[bad]=np.where(f0[0,invalid]>=273.15,3,4);final_diag[:,bad]=diag0[:,invalid]
            snapshot(bad,np.full(len(bad),seconds),f0[:,invalid],diag0[:,invalid])
            keep=~invalid;ids=ids[keep];old=old[keep];s0=s0[:,keep];oldraw=oldraw[keep];f0=f0[:,keep];diag0=diag0[:,keep]
        if not len(ids):continue
        velocity0=oldraw[:,6:9].copy();velocity0[:,2]-=diag0[4]
        midpoint=old+.5*step_s*velocity0
        midheightfraction=history.boundary_fraction(old[:,2],midpoint[:,2],z[0],z[-1])*.5
        midpoint[:,2]=np.clip(midpoint[:,2],z[0],z[-1])
        midraw=history.sample(data,midpoint,z,dx,seconds+.5*step_s,stride);fm=weather_air(midraw)
        halfstate=np.asfortranarray(s0.copy());lib.advance(halfstate,f0,.5*step_s,saturated=saturated,cutoff_diameter=cutoff,
                                                        max_step=growth_step,max_fraction=growth_fraction)
        middiag,_=lib.advance(halfstate,fm,0.,saturated=saturated)
        velocity=midraw[:,6:9].copy();velocity[:,2]-=middiag[4]
        trial=old+step_s*velocity
        heightfraction=history.boundary_fraction(old[:,2],trial[:,2],z[0],z[-1])
        fraction=heightfraction.copy()
        midexit=midheightfraction<.5
        fraction[midexit]=midheightfraction[midexit]
        trial[midexit]=old[midexit]+step_s*velocity0[midexit]
        clipped=trial.copy();clipped[:,2]=np.clip(clipped[:,2],z[0],z[-1])
        fend=weather_air(history.sample(data,clipped,z,dx,seconds+step_s,stride))
        tempfraction=thermal_fraction(f0[0],fend[0])
        midtempfraction=thermal_fraction(f0[0],fm[0])*.5
        thermalcode=np.where(midtempfraction<.5,np.where(fm[0]>=273.15,3,4),np.where(fend[0]>=273.15,3,4))
        tempfraction=np.where(midtempfraction<.5,np.minimum(tempfraction,midtempfraction),tempfraction)
        fraction=np.minimum(fraction,tempfraction)
        effective=step_s*fraction
        middle=old+(trial-old)*fraction[:,None]/2
        middle[:,2]=np.clip(middle[:,2],z[0],z[-1])
        f=weather_air(variable_sample(data,middle,z,dx,seconds+effective/2,stride))
        # A linear crossing estimate can leave midpoint air outside the declared
        # native habit interval in a strongly curved resolved profile. Such a
        # path ends at the preceding accepted state, explicitly, not by clamping T.
        outside=(f[0]<243.15)|(f[0]>=273.15)
        if outside.any():
            bad=ids[outside];status[bad]=np.where(f[0,outside]>=273.15,3,4);final_diag[:,bad]=diag0[:,outside]
            snapshot(bad,np.full(len(bad),seconds),f0[:,outside],diag0[:,outside])
            keep=~outside
            ids=ids[keep];old=old[keep];s0=s0[:,keep];trial=trial[keep];fraction=fraction[keep]
            effective=effective[keep];f=f[:,keep];fend=fend[:,keep];fm=fm[:,keep];tempfraction=tempfraction[keep]
            thermalcode=thermalcode[keep]
        if not len(ids):continue
        updated=np.asfortranarray(s0.copy())
        growth,codes=lib.advance(updated,f,effective,saturated=saturated,cutoff_diameter=cutoff,
                                 max_step=growth_step,max_fraction=growth_fraction)
        ratio=np.divide(growth[8],effective,out=np.zeros_like(effective),where=effective>0)
        position[ids]=old+(trial-old)*(fraction*ratio)[:,None]
        position[ids,2]=np.clip(position[ids,2],z[0],z[-1])
        state[:,ids]=updated;ages[ids]+=growth[8];peak[ids]=np.maximum(peak[ids],growth[7])
        deposited[ids]+=np.maximum(growth[5],0);sublimated[ids]+=np.maximum(-growth[5],0);heat[ids]+=growth[6]
        moist[ids]+=growth[8]*(growth[12]>=growth[11]);total_substeps+=int(growth[9].sum())
        max_balance=max(max_balance,float(np.max(np.abs(growth[13])/np.maximum(updated[0],s0[0]))))
        terminal=fraction<1
        local_status=np.zeros(len(ids),np.int32)
        local_status[terminal]=np.where(trial[terminal,2]<z[0],5,6)
        thermal=terminal&(tempfraction<=fraction+1e-14)
        local_status[thermal]=thermalcode[thermal]
        local_status[codes>0]=codes[codes>0]
        if seconds+step_s==end:local_status[local_status==0]=7
        status[ids]=local_status
        actualtime=seconds+growth[8]
        finalforce=weather_air(variable_sample(data,position[ids],z,dx,actualtime,stride))
        final_diag[:,ids],_=lib.advance(np.asfortranarray(state[:,ids]),finalforce,0.,saturated=saturated)
        save_mask=(local_status>0)|((seconds+step_s)%600==0)
        if save_mask.any():snapshot(ids[save_mask],actualtime[save_mask],finalforce[:,save_mask],final_diag[:,ids[save_mask]])
        if it%int(3600/step_s)==0:print('weather hour',seconds/3600,'active',len(ids),'ended',int(np.count_nonzero(status)),flush=True)
    if np.any(status==0):raise ValueError('Unfinished probe')
    balance=state[0]-initial_mass-deposited+sublimated
    return dict(configuration=dict(diameter_m=diameter,step_s=step_s,field_stride=stride,saturated=saturated,
        cutoff_mass_equivalent_diameter_m=cutoff,growth_max_step_s=growth_step,growth_max_mass_fraction=growth_fraction,diameter_stop_m=1.),
        wall_s=time.monotonic()-start_wall,substeps=total_substeps,maximum_native_mass_balance_relative=max_balance,
        rows_columns=row_fields,rows=np.asarray(rows),states=state,final_diagnostics=final_diag,
        release=seeds,status=status,age_s=ages,peak_diameter_m=peak,deposited_mass_kg=deposited,
        sublimated_mass_kg=sublimated,latent_heat_J=heat,moist_time_s=moist,final_position_m=position,
        initial_mass_kg=initial_mass,mass_balance_residual_kg=balance)


def run(name,**options):
    path=CACHE/'runs'/name
    if path.exists():raise FileExistsError(path)
    record=verify_weather();lib=Library();data=np.load(CACHE/'weather/fields.npy',mmap_mode='r')
    path.mkdir(parents=True)
    source=Path(__file__).read_bytes();(path/'producer.py').write_bytes(source)
    save(path/'request.json',options)
    result=integrate(lib,data,record,**options)
    if Path(__file__).read_bytes()!=source:raise ValueError('Producer changed during execution')
    arrays={k:v for k,v in result.items() if isinstance(v,np.ndarray)}
    np.savez_compressed(path/'raw.npz',**arrays)
    meta={k:v for k,v in result.items() if not isinstance(v,np.ndarray)}
    meta.update(schema='terluna.atmosphere.mcsnow-weather-run/1',name=name,
        producer_sha256=native.sha(__file__),build_record_sha256=native.sha(CACHE/'build_moon/build.json'),
        weather_record_sha256=native.sha(CACHE/'weather/input.json'),raw_sha256=native.sha(path/'raw.npz'),
        terminations={text:int(np.sum(result['status']==code)) for code,text in STATUSES.items() if code},
        peak_dimension_stats_m=history.stats(result['peak_diameter_m']),age_stats_s=history.stats(result['age_s']))
    save(path/'run.json',meta)
    print(name,meta['terminations'],meta['peak_dimension_stats_m'],flush=True)
    return meta


def matrix_requests():
    cases=[]
    for size,diameter in (('small',4e-6),('large',2e-4)):
        for humidity,mode in (('actual',0),('saturated',1),('capped',2)):
            for step in (20.,10.):
                cases.append(dict(name=f'{humidity}_{size}_step{step:g}',
                                  options=dict(diameter=diameter,step_s=step,saturated=mode)))
        for stride in (2,3):
            cases.append(dict(name=f'actual_{size}_fields{stride*10}',
                              options=dict(diameter=diameter,step_s=10.,stride=stride)))
        cases.append(dict(name=f'actual_{size}_cutoff01',
                          options=dict(diameter=diameter,step_s=10.,cutoff=1e-7)))
        cases.append(dict(name=f'actual_{size}_growth_refined',
                          options=dict(diameter=diameter,step_s=10.,growth_step=.25,growth_fraction=.0025)))
    return cases


def run_matrix():
    cases=matrix_requests();save(CACHE/'matrix.json',cases)
    for case in cases:
        name=case['name'];print('START',name,flush=True)
        with (CACHE/(name+'.log')).open('w') as log,contextlib.redirect_stdout(log):
            result=run(name,**case['options'])
        print('DONE',name,result['terminations'],result['peak_dimension_stats_m'],flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['build','weather','run','matrix'])
    parser.add_argument('--name');parser.add_argument('--diameter',type=float,default=4e-6)
    parser.add_argument('--step',type=float,default=20.)
    parser.add_argument('--stride',type=int,default=1)
    parser.add_argument('--saturated',action='store_true')
    parser.add_argument('--cap-water',action='store_true')
    parser.add_argument('--cutoff',type=float,default=2e-7)
    args=parser.parse_args()
    if args.action=='build':
        for world in ('moon','earth'):build(world)
    elif args.action=='weather':prepare_weather()
    elif args.action=='matrix':run_matrix()
    else:run(args.name,diameter=args.diameter,step_s=args.step,stride=args.stride,
             saturated=2 if args.cap_water else int(args.saturated),cutoff=args.cutoff)


if __name__=='__main__':main()
