"""Finite-water SCALE-SDM column controls, using native vapor/heat feedback.

This is a feasibility adapter, not an empirical snow benchmark or storm forecast.
Air dynamics and momentum feedback are disabled. Water initially exists only in
vapor and hydrated aerosol; no source replenishes it during these controls.
The borrowed physics retains its own Earth closure constants and assumptions.
"""
from pathlib import Path
import argparse, json, math, os, subprocess, time, hashlib

import re
import signal
import numpy as np
from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY
from atmosphere.snow.scale_build import CASE_REL

REPO = Path(__file__).resolve().parents[2]
CONSTANTS = REPO / 'shared/constants.json'
DEFAULT_CACHE = REPO / 'research/runs/snow_scale'
# Pinned SCALE scale_const.F90: these are reference-solver closures, not new
# defaults for the terraformed atmosphere. Source identity is in the run product.
SOURCE_CONSTANTS = dict(Rdry=287.04, Rvap=461.46, CPdry=1004.64, DWATR=1000.)

def namelist(name, **items):
    def value(x):
        if isinstance(x, bool): return '.true.' if x else '.false.'
        if isinstance(x, str): return repr(x)
        if isinstance(x, (tuple,list)): return ', '.join(value(a) for a in x)
        if isinstance(x, float): return f'{x:.15e}'.replace('e','d')
        return str(x)
    return '&'+name+'\n'+''.join(f' {k} = {value(v)},\n' for k,v in items.items())+'/\n\n'

def configure(cache, name, lunar=False, duration=60., dt=.4, sd=64., move=False, deposition=True, collisions=True, seed=0):
    if not re.fullmatch(r'[a-z][a-z0-9_]*', name):
        raise ValueError('Case name must be a simple lowercase identifier')
    if not all(math.isfinite(x) and x > 0 for x in (duration, dt, sd)):
        raise ValueError('Duration, step and particle count must be positive and finite')
    if sd != int(sd) or seed < 0:
        raise ValueError('Particle count must be integral and seed nonnegative')
    if any(abs(x/dt-round(x/dt)) > 1e-7 for x in (60., duration)) or duration % 60.:
        raise ValueError('Use whole-minute durations and a step dividing 60 seconds')
    d=Path(cache).resolve()/'runs'/name
    d.mkdir(parents=True, exist_ok=False)
    g=MOON_SURFACE_GRAVITY if lunar else STANDARD_GRAVITY
    common=''
    for label,params in [
        ('PARAM_PRC',dict(PRC_NUM_X=1,PRC_NUM_Y=1)),
        ('PARAM_CONST',dict(CONST_GRAV=g)),
        ('PARAM_INDEX',dict(KMAX=16,IMAX=4,JMAX=2)),
        ('PARAM_GRID',dict(DX=125.,DY=125.,DZ=125.,BUFFER_DZ=1000.,BUFFFACT=1.)),
        ('PARAM_TRACER',dict(TRACER_TYPE='SDM')),
        ('PARAM_ATMOS',dict(ATMOS_DYN_TYPE='OFF',ATMOS_PHY_MP_TYPE='SDM')),
        ('PARAM_STATISTICS',dict(STATISTICS_checktotal=False,STATISTICS_use_globalcomm=False)),
    ]: common+=namelist(label,**params)
    start=dict(TIME_STARTDATE=(0,1,1,0,0,0),TIME_STARTMS=0.)
    init=common+namelist('PARAM_IO',IO_LOG_BASENAME='init_LOG')+namelist('PARAM_TIME',**start)
    init+=namelist('PARAM_ATMOS_VARS',ATMOS_RESTART_OUTPUT=True,ATMOS_RESTART_OUT_BASENAME='init')
    init+=namelist('PARAM_MKINIT',MKINIT_initname='SONDE_PERTURB',flg_sdm=True)
    init+=namelist('PARAM_MKINIT_SOUNDING',ENV_IN_SOUNDING_file='env.txt')
    init+=namelist('PARAM_MKINIT_SONDE_PERTURB',PERTURB_AMP=0.,BBL_THETA=0.)
    init+=namelist('PARAM_BUBBLE',BBL_CZ=1000.,BBL_CX=250.,BBL_CY=125.,BBL_RZ=1000.,BBL_RX=500.,BBL_RY=250.)
    init+=namelist('PARAM_SDMRANDOM',RANDOM_INIT_BASENAME='random_number_init',RANDOM_NUMBER_SEED=seed)
    init+=namelist('PARAM_RANDOM',RANDOM_FIX=True)
    (d/'init.conf').write_text(init)
    run=common+namelist('PARAM_IO',IO_LOG_BASENAME='run_LOG',IO_LOG_ALLNODE=False,IO_LOG_SUPPRESS=False)
    run+=namelist('PARAM_TIME',**start,TIME_DURATION=duration,TIME_DURATION_UNIT='SEC',TIME_DT=dt,
        TIME_DT_UNIT='SEC',TIME_DT_ATMOS_PHY_MP=dt,TIME_DT_ATMOS_PHY_MP_UNIT='SEC')
    run+=namelist('PARAM_ATMOS_VARS',ATMOS_RESTART_OUTPUT=False,
        ATMOS_RESTART_IN_BASENAME='init_00000000000.000',ATMOS_VARS_CHECKRANGE=False)
    run+=namelist('PARAM_ATMOS_PHY_MP_VARS',use_tendency_mp=False)
    run+=namelist('PARAM_ATMOS_PHY_MP',docondensation=True,doautoconversion=collisions,doprecipitation=True)
    run+=namelist('PARAM_ATMOS_PHY_MP_SDM',RANDOM_IN_BASENAME='random_number_init',
        RANDOM_OUT_BASENAME='',SD_OUT_BASENAME='',domovement=move,sdm_cold=True,domeltfreeze=True,
        dosublimation=deposition,sdm_rqc2qr=40e-6,sdm_dtcmph=(dt/4,dt/2,dt,dt,dt/4),sdm_rdnc=0.,
        sdm_sdnmlvol=125.**3,sdm_inisdnc=sd,sdm_zlower=0.,sdm_zupper=1875.,sdm_aslset=-1,
        sdm_fctr2multi=3.,sdm_aslmw=-999.,sdm_extbuf=15,sdm_aslfmrate=0.,sdm_aslfmdt=0.,sdm_aslfmsdnc=0.,
        sdm_colbrwn=0,sdm_colkrnl=2,sdm_mvexchg=0,sdm_nadjdt=0.,sdm_nadjvar=0,sdm_dmpvar=10,
        sdm_dmpitva=60.,sdm_dmpnskip=1,sdm_dmpitvb=60.,sdm_dmpitvl=60.,sdm_dmpsdsiz=0.)
    run+=namelist('PARAM_ATMOS_REFSTATE',ATMOS_REFSTATE_TYPE='INIT')
    run+=namelist('PARAM_ATMOS_BOUNDARY',ATMOS_BOUNDARY_TYPE='CONST',ATMOS_BOUNDARY_USE_VELZ=True,
        ATMOS_BOUNDARY_VALUE_VELZ=0.,ATMOS_BOUNDARY_TAUZ=10.)
    run+=namelist('PARAM_HISTORY',HISTORY_DEFAULT_BASENAME='history',HISTORY_DEFAULT_TINTERVAL=60.,
        HISTORY_DEFAULT_TUNIT='SEC',HISTORY_DEFAULT_TAVERAGE=False,HISTORY_DEFAULT_DATATYPE='REAL8',
        HISTORY_OUTPUT_STEP0=True)
    run+=namelist('PARAM_HIST',HIST_BND=False)
    for field in ('DENS','U','V','W','PT','T','PRES','RH','RHI','QV','QHYD_sd','QC_sd','QR_sd','QI_sd','QS_sd','QG_sd','PREC_sd','RAIN_sd','SNOW_sd'):
        run+=namelist('HISTITEM',item=field)
    run+=namelist('PARAM_MONITOR',MONITOR_STEP_INTERVAL=1000)
    run+=namelist('PARAM_USER',HEAT_TIME=0.,HEAT_RATE=0.)
    (d/'run.conf').write_text(run)
    # Approximate isothermal moist hydrostatic sounding; SCALE recomputes its
    # own consistent pressure/density. Record actual initialized T and RH.
    # Thermodynamic coefficients are the pinned SCALE closure constants.
    rd,rv,cp=(SOURCE_CONSTANTS[k] for k in ('Rdry','Rvap','CPdry'))
    temp=253.15
    # Deliberately prescribed initial supersaturation, not measured atmosphere.
    ew=6.112*math.exp(17.67*(-20.)/(-20.+243.5))*100
    def sounding(z):
        p=70000.*math.exp(-g*z/(rd*temp))
        q=1.03*(rd/rv)*ew/(p-(1-rd/rv)*ew)
        theta=temp*(100000./p)**(rd/cp)
        return theta,q
    theta,q=sounding(0.)
    lines=[f'700.0 {theta:.15g} {q*1000:.15g}\n']
    for z in range(25,3001,25):
        theta,q=sounding(z)
        lines.append(f'{z} {theta:.15g} {q*1000:.15g} 0.0 0.0\n')
    (d/'env.txt').write_text(''.join(lines))
    settings=dict(name=name,gravity_m_s2=g,duration_s=duration,dt_s=dt,sd_per_full_cell=sd,movement=move,
        deposition=deposition,collisions=collisions,seed=seed,grid=[16,4,2],spacing_m=[125.,125.,125.],
        supplied_profile='Approximate isothermal -20 C, 700 hPa base, 103% water RH; not a climate profile.')
    (d/'settings.json').write_text(json.dumps(settings,indent=2)+'\n')
    return d

def execute(d, timeout=500):
    d=Path(d).resolve()
    case=d.parents[1]/'work'/CASE_REL
    env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
    local_lib=d.parents[1]/'deps/usr/lib/x86_64-linux-gnu'
    if local_lib.is_dir():
        env['LD_LIBRARY_PATH']=str(local_lib)+os.pathsep+env.get('LD_LIBRARY_PATH','')
    runs=[]
    for exe,conf in [('scale-les_init','init.conf'),('scale-les','run.conf')]:
        t=time.monotonic()
        with (d/(exe+'.log')).open('x') as log:
            try:
                p=subprocess.Popen([str(case/exe),conf],cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
                rc=p.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                os.killpg(p.pid,signal.SIGKILL)
                p.wait()
                rc='timeout'
        runs.append(dict(executable=exe,sha256=hashlib.sha256((case/exe).read_bytes()).hexdigest(),exit_code=rc,elapsed_s=time.monotonic()-t))
        (d/'execution.json').write_text(json.dumps(runs,indent=2)+'\n')
        print(d.name,exe,rc,runs[-1]['elapsed_s'],flush=True)
        if rc!=0: break
    return runs

def _dataset(path):
    # Most repository checks do not need this optional output-reader dependency.
    from netCDF4 import Dataset
    return Dataset(path)


def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def weighted_quantile(x,w,quantiles=(.1,.5,.9,.99)):
    x,w=np.asarray(x),np.asarray(w)
    if x.shape != w.shape or not np.isfinite(x).all() or not np.isfinite(w).all() or (w<0).any():
        raise ValueError('Weighted quantiles require finite values and nonnegative weights')
    keep=w>0
    x,w=x[keep],w[keep]
    if not len(x): return [None for _ in quantiles]
    ix=np.argsort(x); x,w=x[ix],w[ix]
    c=np.cumsum(w)/w.sum()
    return [float(x[min(np.searchsorted(c,q),len(x)-1)]) for q in quantiles]

def particles(path, area):
    with _dataset(path) as nc:
        v={k:np.array(a[:]) for k,a in nc.variables.items()}
    valid=(v['sd_z']>=0)&(v['sd_n']>0)
    ice=valid&(v['sd_liqice']==10)
    liq=valid&(v['sd_liqice']==1)
    n=v['sd_n'][ice].astype(float)
    re,rp,rho=(v[k][ice] for k in ['sdi_re','sdi_rp','sdi_rho'])
    if any((~np.isfinite(x)).any() or (x<=0).any() for x in (re,rp,rho)):
        raise ValueError('Invalid live ice geometry in native SCALE output')
    m=4/3*np.pi*re**2*rp*rho
    size=2*np.maximum(re,rp)
    weight=n*m
    lm=float(np.sum(v['sd_n'][liq]*4/3*np.pi*v['sd_r'][liq]**3*SOURCE_CONSTANTS['DWATR']))
    im=float(weight.sum())
    return dict(file=path.name,sha256=digest(path),quantile_probabilities=[.1,.5,.9,.99],valid_superparticles=int(valid.sum()),ice_superparticles=int(ice.sum()),ice_number_per_m2=float(n.sum()/area),ice_kg_m2=im/area,liquid_kg_m2=lm/area,
                ice_max_dimension_number_quantiles_m=weighted_quantile(size,n),ice_max_dimension_mass_quantiles_m=weighted_quantile(size,weight),ice_max_dimension_max_m=float(size.max()) if len(size) else None,
                ice_aspect_ratio_mass_quantiles=weighted_quantile(rp/re,weight),ice_density_mass_quantiles_kg_m3=weighted_quantile(rho,weight),ice_fall_speed_mass_quantiles_m_s=weighted_quantile(v['sd_vz'][ice],weight),
                ice_aggregate_mass_fraction=float(weight[v['sdi_nmono'][ice]>1].sum()/im) if im else None,ice_rime_mass_fraction=float(np.sum(n*v['sdi_mrime'][ice])/im) if im else None,
                max_liquid_diameter_m=float(2*v['sd_r'][liq].max()) if liq.any() else None)


def water_history(rho, qv, qhyd, accumulated_precip_m, cell_volume, area):
    """Close water per horizontal area using simultaneous atmospheric records.

    Density varies as vapor leaves the gas, so integrating mixing ratios alone
    would introduce a spurious residual. PREC_sd is accumulation, not a rate.
    """
    rho,qv,qhyd,prec=(np.asarray(x) for x in (rho,qv,qhyd,accumulated_precip_m))
    if rho.ndim!=4 or qv.shape!=rho.shape or qhyd.shape!=rho.shape or prec.shape!=(rho.shape[0],*rho.shape[2:]):
        raise ValueError('Water history dimensions must be (time,z,y,x) and (time,y,x)')
    if any(not np.isfinite(x).all() for x in (rho,qv,qhyd,prec)):
        raise ValueError('Non-finite water history')
    if (rho<=0).any() or (qv<0).any() or (qhyd<0).any() or (prec<0).any() or cell_volume<=0 or area<=0:
        raise ValueError('Nonphysical water history or integration volume')
    vapor=(rho*qv).sum(axis=(1,2,3))*cell_volume/area
    hydro=(rho*qhyd).sum(axis=(1,2,3))*cell_volume/area
    precip=prec.mean(axis=(1,2))*SOURCE_CONSTANTS['DWATR']
    return vapor,hydro,precip,vapor+hydro+precip


def aggregation_snapshots(paths):
    """Detect represented aggregates at saved times; not an event counter."""
    rows=[]
    for path in paths:
        with _dataset(path) as nc:
            active=(np.array(nc['sd_z'][:])>=0)&(np.array(nc['sd_n'][:])>0)
            ice=active&(np.array(nc['sd_liqice'][:])==10)
            monomers=np.array(nc['sdi_nmono'][:])[ice]
        rows.append(dict(file=path.name,sha256=digest(path),
                         aggregate_superparticles=int(np.count_nonzero(monomers>1)),
                         maximum_monomers=int(monomers.max()) if len(monomers) else 0))
    return rows


def analyze(path):
    path=Path(path)
    settings=json.loads((path/'settings.json').read_text())
    dz,dx,dy=settings['spacing_m'];k,i,j=settings['grid']; area=i*j*dx*dy
    with _dataset(path/'history.pe000000.nc') as nc:
        a={k:np.array(v[:]) for k,v in nc.variables.items()}
    execution=json.loads((path/'execution.json').read_text())
    if len(execution)!=2 or any(x['exit_code']!=0 for x in execution):
        raise ValueError('Incomplete or failed SCALE execution')
    if not np.array_equal(a['time'],np.arange(0.,settings['duration_s']+1.,60.)):
        raise ValueError('Incomplete or unexpected history times')
    rho=a['DENS']; vol=dz*dx*dy
    integral=lambda v: np.sum(v,axis=(1,2,3))*vol/area
    vapor,hydro,precip,total=water_history(rho,a['QV'],a['QHYD_sd'],a['PREC_sd'],vol,area)
    ice=integral(rho*(a['QI_sd']+a['QS_sd']+a['QG_sd']))
    snow=np.mean(a['SNOW_sd'],axis=(1,2))*SOURCE_CONSTANTS['DWATR']
    dry=integral(rho*(1-a['QV']))
    temp=np.sum(rho*a['T'],axis=(1,2,3))/rho.sum(axis=(1,2,3))
    series=[dict(time_s=float(t),vapor_kg_m2=float(v),hydrometeor_kg_m2=float(h),ice_kg_m2=float(ic),precipitation_kg_m2=float(p),snowfall_kg_m2=float(s),total_water_kg_m2=float(tw),gas_mass_weighted_temperature_k=float(te),mean_water_rh_percent=float(rh),mean_ice_rh_percent=float(ri))
        for t,v,h,ic,p,s,tw,te,rh,ri in zip(a['time'],vapor,hydro,ice,precip,snow,total,temp,
            a['RH'].mean(axis=(1,2,3)),a['RHI'].mean(axis=(1,2,3)))]
    snapshots=sorted(path.glob('SD_output_NetCDF_*'))
    if len(snapshots)!=len(a['time']):
        raise ValueError('Missing particle snapshots for the declared output cadence')
    result=dict(settings=settings,execution=json.loads((path/'execution.json').read_text()),history_sha256=digest(path/'history.pe000000.nc'),
        water_ledger='sum(DENS * (QV + QHYD_sd) * cell_volume) / horizontal_area + 1000 * mean(PREC_sd). All fields from the same atmospheric history record. Baseline is first post-initialization record; initializer adds hydrated aerosol water before that baseline.',
        particle_timing='Native particle dumps precede sdm_calc; atmospheric hydrometeor histories follow it. Same nominal label differs by one main microphysics step. Particle masses are not mixed with history vapor in the water ledger.',
        max_relative_water_residual=float(np.max(np.abs(total-total[0]))/total[0]),max_absolute_water_residual_kg_m2=float(np.max(np.abs(total-total[0]))),max_relative_dry_air_residual=float(np.max(np.abs(dry-dry[0]))/dry[0]),
        temperature_change_k=float(temp[-1]-temp[0]),vapor_loss_kg_m2=float(vapor[0]-vapor[-1]),initial_temperature_range_k=[float(a['T'][0].min()),float(a['T'][0].max())],max_air_speed_m_s=float(np.sqrt(a['U']**2+a['V']**2+a['W']**2).max()),
        series=series,particle_snapshots=[particles(snapshots[0],area),particles(snapshots[-1],area)],
        aggregation_snapshot_scan=aggregation_snapshots(snapshots),
        input_sha256={name:digest(path/name) for name in ('init.conf','run.conf','env.txt','settings.json')})
    (path/'analysis.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return result


def summarize(cache, output):
    """Recompute the declared cases and bind them to preserved build evidence."""
    cache=Path(cache)
    provenance=json.loads((cache/'source_manifest.json').read_text())
    cases={name:analyze(cache/'runs'/name) for name in provenance['cases']}
    for result in cases.values():
        for run in result['execution']:
            if run['sha256'] not in provenance['allowed_executable_sha256'][run['executable']]:
                raise ValueError('Case executable is absent from the source/build manifest')
    base=cases['lunar_settling']
    def relative(case):
        row=cases[case]
        return dict(
            airborne_ice_percent=100*(row['series'][-1]['ice_kg_m2']/base['series'][-1]['ice_kg_m2']-1),
            accumulated_snow_percent=100*(row['series'][-1]['snowfall_kg_m2']/base['series'][-1]['snowfall_kg_m2']-1),
            mass_median_dimension_percent=100*(row['particle_snapshots'][-1]['ice_max_dimension_mass_quantiles_m'][1]/base['particle_snapshots'][-1]['ice_max_dimension_mass_quantiles_m'][1]-1),
            mass_p99_dimension_percent=100*(row['particle_snapshots'][-1]['ice_max_dimension_mass_quantiles_m'][3]/base['particle_snapshots'][-1]['ice_max_dimension_mass_quantiles_m'][3]-1))
    result=dict(schema='terluna.atmosphere.scale-column/1',
        producer=dict(path='atmosphere/snow/scale_column.py',sha256=digest(Path(__file__)),
                      build_adapter_sha256=digest(Path(__file__).with_name('scale_build.py'))),
        shared_constants_sha256=digest(CONSTANTS),source_manifest_sha256=digest(cache/'source_manifest.json'),
        provenance=provenance,cases=cases,
        comparisons={name:relative(name) for name in ('lunar_half_step','lunar_double_particles','lunar_seed1')},
        maximum_relative_water_residual=max(r['max_relative_water_residual'] for r in cases.values()),
        evidence='Executed finite-water mixed-phase columns with native vapor/heat feedback, fixed air, settling, riming and evolving spheroidal ice geometry. Numerical feasibility controls; no empirical validation or full-storm convergence.',
        reading_rule='Read only as the declared 2 km, 30 minute controls and 60 second startup checks. Particle snapshots precede atmospheric history by one microphysics step. Initial hydration water is included at the first history baseline. The saved snapshots do not establish behavior in an aggregation-heavy cloud; their narrow size distributions cannot qualify the historical 15–60 cm tail or its optics.')
    Path(output).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('name')
    parser.add_argument('--cache',type=Path,default=DEFAULT_CACHE)
    parser.add_argument('--lunar',action='store_true')
    parser.add_argument('--duration',type=float,default=1800.)
    parser.add_argument('--dt',type=float,default=.4)
    parser.add_argument('--sd',type=float,default=64.)
    parser.add_argument('--move',action='store_true')
    parser.add_argument('--no-deposition',action='store_true')
    parser.add_argument('--no-collisions',action='store_true')
    parser.add_argument('--seed',type=int,default=0)
    parser.add_argument('--timeout',type=float,default=500.)
    parser.add_argument('--analyze-only',action='store_true')
    parser.add_argument('--summarize-to',type=Path)
    args=parser.parse_args()
    if args.summarize_to:
        result=summarize(args.cache,args.summarize_to)
        print(f"Summarized {len(result['cases'])} declared cases in {args.summarize_to}")
        return
    if args.analyze_only:
        folder=args.cache/'runs'/args.name
    else:
        folder=configure(args.cache,args.name,args.lunar,args.duration,args.dt,args.sd,
                         args.move,not args.no_deposition,not args.no_collisions,args.seed)
        results=execute(folder,args.timeout)
        if len(results)!=2 or any(row['exit_code']!=0 for row in results):
            raise SystemExit('SCALE execution failed; preserved logs identify the cause')
    result=analyze(folder)
    print(json.dumps({key:result[key] for key in ('max_relative_water_residual',
        'temperature_change_k','vapor_loss_kg_m2','particle_snapshots')},indent=2))


if __name__=='__main__':
    main()
