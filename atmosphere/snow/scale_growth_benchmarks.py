"""Native SCALE empirical growth controls across representable crystal habits.

Use the same stationary eight-crystal driver as scale_crystal_growth, with full
published a and c tracked separately. Baseline physics is unchanged. Reject
three-minute fits that cannot be represented without exceeding bulk ice density.
"""
from pathlib import Path
import argparse
import json
import math
import re
import tempfile
import shutil
import numpy as np
from atmosphere.snow import scale_crystal_growth as first
from atmosphere.snow import takahashi_growth as lab
from shared.constants import STANDARD_GRAVITY

lunar=first.lunar
CACHE=lunar.CACHE
OUT=lunar.ROOT/'atmosphere/snow/results/scale_growth_benchmarks.json'
CONST_PATH='scalelib/src/common/scale_const.F90'
CONST_SHA='256841aff9d6358cb073a199d05424dc67f2c52a2cfe21666419b8a0ec26a8c2'
TEMPERATURES=(-3.7,-5.3,-10.6,-14.4,-20.1)
QUANTITIES=('mass_kg','equatorial_dimension_m','polar_dimension_m','diameter_m')


def ice_density(cache=CACHE):
    path=Path(cache)/'lunar/work'/CONST_PATH
    if lunar.base.digest(path)!=CONST_SHA:raise ValueError('Unidentified native constants')
    match=re.search(r'CONST_DICE\s*=\s*([0-9.]+)_RP',path.read_text())
    if not match:raise ValueError('Missing native bulk ice density')
    return float(match.group(1))


def case_name(temperature_C,dt,control=None):
    name='growth_sweep_'+str(abs(temperature_C)).replace('.','p')+'_dt'+format(dt,'g')
    return name if control is None else name+'_control'+str(control)


def configure(temperature_C,dt,cache=CACHE,control=None):
    cache=Path(cache)
    if lunar.base.digest(cache/first.SOURCE_PDF)!=lab.SOURCE_SHA:
        raise ValueError('Unidentified empirical source')
    if not math.isfinite(dt) or dt<=0 or 60/dt!=round(60/dt):
        raise ValueError('Use a positive timestep dividing one minute')
    start=lab.initial(temperature_C,ice_density(cache))
    duration=(lab.FITS[temperature_C][0]-3)*60.
    name=case_name(temperature_C,dt,control)
    d=lunar.configure(name,dt=5.,duration=duration,layers=4,capacity_per_layer=4,
                      collisions=False,cache=cache)
    for filename in ('init.conf','run.conf'):
        text=(d/filename).read_text();text=lunar.boxes.set_value(text,'CONST_GRAV',STANDARD_GRAVITY)
        if filename=='run.conf':
            for k,v in dict(TIME_DT=dt,TIME_DT_ATMOS_PHY_MP=dt,sdm_dtcmph=[dt]*5,
                           domovement=False,sdm_dmpvar=10,sdm_dmpitva=60.,
                           sdm_dmpitvb=60.,HISTORY_DEFAULT_TINTERVAL=60.).items():
                text=lunar.boxes.set_value(text,k,v)
        (d/filename).write_text(text)
    s=json.loads((d/'settings.json').read_text())
    depth=s['depth_m'];count=s['capacity'];zeros=np.zeros(count)
    pop=dict(n=np.zeros(count,dtype=np.int64),rk=np.full(count,-999.999),x=zeros.copy(),
             y=zeros.copy(),z=np.full(count,-999.999),r=zeros.copy(),
             re=np.full(count,start['equatorial_radius_m']),rp=np.full(count,start['polar_radius_m']),
             rho=np.full(count,start['density_kg_m3']),phase=np.full(count,10,dtype=np.int16),
             monomers=np.ones(count,dtype=np.int32))
    for idx,(i,j) in enumerate((i,j) for i in range(4) for j in range(2)):
        pop['n'][idx]=1;pop['rk'][idx]=2.5;pop['z'][idx]=depth/(2*s['layers'])
        pop['x'][idx]=i+.5;pop['y'][idx]=j+.5
    (d/'particles.pe000000').unlink();lunar.boxes.write_restart(d/'particles.pe000000',pop)
    temperature=273.15+temperature_C;pressure=first.PRESSURE
    vapor_pressure=float(lunar.flakes.liquid_saturation_pressure(temperature))
    lines=(d/'lunar_forcing.txt').read_text().splitlines()
    values=list(map(float,lines[1].split()));values[2]=0.;values[3]=1e-6;values[5]=0.
    lines[1]=' '.join(format(x,'.17g') for x in values)
    lines[3:]=[' '.join(format(x,'.17g') for x in (temperature,pressure,vapor_pressure))]*s['layers']
    (d/'lunar_forcing.txt').write_text('\n'.join(lines)+'\n')
    s.update(dt_s=dt,gravity_m_s2=STANDARD_GRAVITY,initial=start,temperature_C=temperature_C,
             temperature_K=temperature,pressure_Pa=pressure,vapor_pressure_Pa=vapor_pressure,
             source_pdf_sha256=lab.SOURCE_SHA,laboratory_start_age_min=3.,movement=False,
             initial_crystals=8,deposition_ceiling_kg_m3_s=values[3],gather_start_s=0.,
             profile=[dict(t=temperature,p=pressure,e=vapor_pressure)]*s['layers'],
             forcing='Stationary isolated crystals with native terminal-speed ventilation in fixed water-saturated air; no liquid, collisions, injection or transport.',
             unmatched='Empirical Table 2 fits at Earth gravity and 1010 hPa. Planar thickness is inferred, fits lack pointwise errors, piecewise fits may have discontinuous joins. No aggregate or lunar validation.',
             initial_water_kg_m2=start['mass_kg'],
             ledger_convention='Add known initial ice inventory to native empty-column source ledger.',
             producer_sha256=lunar.base.digest(Path(__file__)),lab_producer_sha256=lunar.base.digest(Path(lab.__file__)),
             native_constants_sha256=lunar.base.digest(cache/'lunar/work'/CONST_PATH))
    if control is not None:
        if control not in (0,1):raise ValueError('Unknown growth control')
        manifest=json.loads((cache/'growth_controls/manifest.json').read_text())
        s.update(build=manifest,build_directory='growth_controls',growth_control=control)
        (d/'growth_control.txt').write_text(str(control)+'\n')
    lunar.boxes.save(d/'settings.json',s)
    names=['init.conf','run.conf','env.txt','particles.pe000000','lunar_forcing.txt','settings.json']
    if control is not None:names.append('growth_control.txt')
    lunar.boxes.save(d/'inputs.json',{n:lunar.base.digest(d/n) for n in names})
    return d


def analyze(d):
    from netCDF4 import Dataset
    d=Path(d);s=json.loads((d/'settings.json').read_text());rows=[]
    hashes=json.loads((d/'inputs.json').read_text())
    if any(lunar.base.digest(d/n)!=h for n,h in hashes.items()):raise ValueError('Changed native inputs')
    execution=json.loads((d/'execution.json').read_text());expected=s['build']['build']['executable_sha256']
    if len(execution)!=2 or any(e['exit_code']!=0 or e['sha256']!=expected[e['executable']] for e in execution):
        raise ValueError('Incomplete or unidentified native execution')
    ledger=np.atleast_2d(np.loadtxt(d/'lunar_ledger.txt'));initial=s['initial_water_kg_m2']
    if not np.isfinite(ledger).all():raise ValueError('Nonfinite ledger')
    if abs(ledger[-1,0]-s['duration_s'])>1e-8:raise ValueError('Incomplete duration')
    for path in sorted(d.glob('SD_output_NetCDF_*.pe000000')):
        with Dataset(path) as nc:v={k:np.asarray(x[:]) for k,x in nc.variables.items()}
        live=(v['sd_n']>0)&(v['sd_z']>=0)
        if live.sum()!=8 or np.any(v['sd_n'][live]!=1) or np.any(v['sdi_nmono'][live]!=1):
            raise ValueError('Isolated population changed')
        a,c,rho=(v[k][live] for k in ('sdi_re','sdi_rp','sdi_rho'));mass=4*math.pi/3*a*a*c*rho
        for values in (a,c,rho,mass):
            if not np.isfinite(values).all() or np.any(values<=0):raise ValueError('Nonpositive or nonfinite ice state')
        speed=v['sd_vz'][live]
        if not np.isfinite(speed).all() or np.any(speed<0):raise ValueError('Invalid fall speed')
        for values in (a,c,rho,mass,speed):
            if np.ptp(values)>max(float(abs(values).max()),1e-30)*1e-12:raise ValueError('Identical crystals differ')
        elapsed=float(re.search(r'(\d{11}\.\d{3})',path.name).group(1));age=3+elapsed/60
        row=dict(age_min=age,elapsed_s=elapsed,mass_kg=float(mass.mean()),
                 equatorial_dimension_m=float(2*a.mean()),polar_dimension_m=float(2*c.mean()),
                 diameter_m=float((2*np.maximum(a,c)).mean()),density_kg_m3=float(rho.mean()),
                 terminal_speed_m_s=float(v['sd_vz'][live].mean()),sha256=lunar.base.digest(path))
        obs=lab.observed(s['temperature_C'],age);row['observed_fit']=obs
        row['observed_fall_speed_m_s']=lab.fall_speed(s['temperature_C'],age)
        row['fall_speed_ratio_to_fit']=row['terminal_speed_m_s']/row['observed_fall_speed_m_s'] if elapsed else None
        # Initial native fall speed predates our prescribed gas-state reset.
        row['ratio_to_fit']={k:row[k]/obs[k] for k in QUANTITIES}
        matching=np.flatnonzero(np.abs(ledger[:,0]-elapsed)<1e-8)
        if elapsed==0:ledger_mass=initial
        elif len(matching)==1:ledger_mass=ledger[matching[0],5]
        else:ledger_mass=None  # Native ledger is every 60 calls, snapshots every 60 seconds.
        row['ledger_mass_error_kg']=None if ledger_mass is None else abs(row['mass_kg']-ledger_mass)
        if ledger_mass is not None and row['ledger_mass_error_kg']>row['mass_kg']*1e-9:raise ValueError('Particle/ledger time mismatch')
        rows.append(row)
    for key in QUANTITIES:
        if not math.isclose(rows[0][key],s['initial'][key],rel_tol=1e-12):raise ValueError('Initial particle state differs from fitted input')
    if any(not np.isfinite(np.asarray([r[k] for r in rows])).all() for k in QUANTITIES):raise ValueError('Nonfinite growth trajectory')
    matched=sum(r['ledger_mass_error_kg'] is not None for r in rows)
    if matched<3 or rows[-1]['ledger_mass_error_kg'] is None:raise ValueError('Insufficient exact snapshot/ledger timing controls')
    ages=[r['age_min'] for r in rows]
    if ages!=list(range(3,lab.FITS[s['temperature_C']][0]+1)):raise ValueError('Incomplete snapshots')
    with Dataset(d/'history.pe000000.nc') as nc:
        fields={k:np.asarray(nc[k][:])[1:] for k in ('T','PRES','RH')}
    if any(not np.isfinite(v).all() for v in fields.values()):raise ValueError('Nonfinite atmospheric history')
    if np.max(abs(fields['T']-s['temperature_K']))>1e-9 or np.max(abs(fields['PRES']-s['pressure_Pa']))>1e-6:
        raise ValueError('Prescribed air state changed')
    if np.max(abs(fields['RH']-100))>.1:raise ValueError('Saturation convention discrepancy exceeds 0.1 percentage point')
    increment=np.diff(np.r_[0.,ledger[:,2]]);steps=np.diff(np.r_[0.,ledger[:,0]])
    budget=s['deposition_ceiling_kg_m3_s']*s['depth_m']/s['layers']*steps
    water=dict(ice_deposition_relative_residual=float(np.max(abs(ledger[:,5]-initial-ledger[:,2]))/ledger[-1,5]),
               corrected_maximum_absolute_residual_kg_m2=float(np.max(abs(ledger[:,7]-initial))),
               maximum_used_fraction_of_growth_ceiling=float(np.max(increment/budget)))
    if water['ice_deposition_relative_residual']>1e-11 or water['maximum_used_fraction_of_growth_ceiling']>1e-4:
        raise ValueError('Water accounting failed or deposition ceiling active')
    result=dict(settings=s,execution=execution,snapshots=rows,water=water,input_sha256=hashes,
                exactly_matched_snapshot_ledger_times=matched,
                air_ranges={k:[float(v.min()),float(v.max())] for k,v in fields.items()},
                ledger_sha256=lunar.base.digest(d/'lunar_ledger.txt'))
    lunar.boxes.save(d/'benchmark_analysis.json',result)
    return result


def verify_inputs(spec,d):
    with tempfile.TemporaryDirectory(prefix='terluna-growth-benchmark-') as tmp:
        cache=Path(tmp);(cache/'lunar').mkdir();(cache/'limitations_sources').mkdir()
        shutil.copyfile(CACHE/'lunar/manifest.json',cache/'lunar/manifest.json')
        (cache/'lunar/work').symlink_to(CACHE/'lunar/work',target_is_directory=True)
        (cache/first.SOURCE_PDF).symlink_to(CACHE/first.SOURCE_PDF)
        if spec.get('control') is not None:
            (cache/'growth_controls').mkdir()
            shutil.copyfile(CACHE/'growth_controls/manifest.json',cache/'growth_controls/manifest.json')
        copy=configure(cache=cache,**spec)
        names=['init.conf','run.conf','env.txt','particles.pe000000','lunar_forcing.txt']
        if spec.get('control') is not None:names.append('growth_control.txt')
        matches={n:lunar.base.digest(d/n)==lunar.base.digest(copy/n) for n in names}
        if not all(matches.values()):raise ValueError('Input reproduction failed')
    return matches


def native_parity(old,new):
    from netCDF4 import Dataset
    files=sorted(p.name for p in old.glob('SD_output_NetCDF_*.pe000000'))
    if files!=sorted(p.name for p in new.glob('SD_output_NetCDF_*.pe000000')):
        raise ValueError('Native snapshot file sets differ')
    counts={}
    for filename in files+['history.pe000000.nc']:
        with Dataset(old/filename) as a,Dataset(new/filename) as b:
            if set(a.variables)!=set(b.variables):raise ValueError('Native field sets differ')
            for field in a.variables:
                if not np.array_equal(np.asarray(a[field][:]),np.asarray(b[field][:]),equal_nan=True):
                    raise ValueError(f'Zero-change native parity failed: {filename} {field}')
            counts[filename]=len(a.variables)
    same=lunar.base.digest(old/'lunar_ledger.txt')==lunar.base.digest(new/'lunar_ledger.txt')
    if not same:raise ValueError('Zero-change ledger parity failed')
    return dict(ledger_bytes_identical=same,identical_native_array_counts=counts)


def run(controls=False):
    specs=([dict(temperature_C=t,dt=dt) for t in TEMPERATURES for dt in (5.,1.)]
           if not controls else [dict(temperature_C=-14.4,dt=dt,control=mode) for mode in (0,1) for dt in (5.,1.)])
    producers=CACHE/'growth_benchmarks_producers';producers.mkdir(exist_ok=True)
    for source in (Path(__file__),Path(lab.__file__)):
        frozen=producers/(lunar.base.digest(source)+'.py')
        if not frozen.exists():frozen.write_bytes(source.read_bytes())
        if lunar.base.digest(frozen)!=lunar.base.digest(source):raise ValueError('Changed frozen experiment producer')
    results={}
    for spec in specs:
        d=CACHE/'runs'/case_name(**spec)
        if not d.exists():d=configure(**spec)
        if not (d/'execution.json').exists():lunar.execute(d)
        result=analyze(d);result['physics_input_matches']=verify_inputs(spec,d)
        results[d.name]=result
        r=result['snapshots'][-1]
        print(d.name,{k:round(r['ratio_to_fit'][k],5) for k in QUANTITIES},flush=True)
    path=OUT if not controls else OUT.with_name('scale_growth_controls.json')
    rejection={}
    for t in lab.FITS:
        try:lab.initial(t,ice_density())
        except ValueError as e:rejection[str(t)]=str(e)
    refinements=[]
    for spec in specs:
        if spec['dt']!=1:continue
        fine=results[case_name(**spec)]['snapshots']
        coarse=results[case_name(**{**spec,'dt':5.})]['snapshots']
        refinements.append(dict(temperature_C=spec['temperature_C'],control=spec.get('control'),
            endpoint_relative_change={k:fine[-1][k]/coarse[-1][k]-1 for k in QUANTITIES},
            maximum_trajectory_relative_change={k:max(abs(a[k]/b[k]-1) for a,b in zip(fine,coarse)) for k in QUANTITIES}))
    parity={}
    if not controls:
        for dt in (5,1):
            old=CACHE/'runs'/f'growth_earth_observed_dt{dt}'
            new=CACHE/'runs'/case_name(-14.4,dt)
            parity[str(dt)]=native_parity(old,new)
    else:
        for dt in (5,1):
            old=CACHE/'runs'/case_name(-14.4,dt)
            new=CACHE/'runs'/case_name(-14.4,dt,0)
            parity[str(dt)]=native_parity(old,new)
    if not all(parity.values()):raise ValueError('Zero-change ledger parity failed')
    product=dict(schema='terluna.atmosphere.scale-growth-benchmarks/1',
        evidence='Empirical fit comparison and numerical/representation diagnostics; no calibrated growth repair or lunar forecast.',
        producer_sha256=lunar.base.digest(Path(__file__)),lab_producer_sha256=lunar.base.digest(Path(lab.__file__)),
        source_pdf_sha256=lab.SOURCE_SHA,initial_spheroid_rejections=rejection,
        baseline_ledger_parity=parity,refinements=refinements,cases=results)
    lunar.boxes.save(path,product)
    return product


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--controls',action='store_true')
    run(parser.parse_args().controls)
