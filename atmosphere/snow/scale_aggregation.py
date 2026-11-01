"""Exercise the pinned SCALE collision solver using native particle restarts.

The liquid control uses the published exponential/Golovin problem. Ice controls
borrow the X2M04 initial mass distribution and air state of Hashino & Tripoli
(2011), retaining SCALE's different collection and geometry closures. They are
numerical controls, not a reproduction of SHIPS or an observed snowstorm.
No collision, drag, or geometry code is replaced by this adapter.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import signal
import struct
import subprocess
import time

import numpy as np

from atmosphere.snow import scale_column as column
from atmosphere.snow.scale_build import CASE_REL
from shared.constants import STANDARD_GRAVITY, MOON_SURFACE_GRAVITY

REPO = Path(__file__).resolve().parents[2]
CACHE = REPO / 'research/runs/snow_scale'
GRID = (4, 4, 2)  # z, x, y; 32 stationary, noncommunicating particle cells
DX = 125.
CELL_VOLUME = DX**3
NCELL = math.prod(GRID)
# Published reference inputs, not new lunar atmosphere constants.
ICE_DENSITY = 300.
ICE_MEAN_MASS = 4 * math.pi / 3 * ICE_DENSITY * (200e-6)**3
ICE_WATER_CONTENT = 1e-4  # kg/m3; 10^-7 g/cm3 in Hashino & Tripoli (2011)
LIQUID_NUMBER = 2**23
LIQUID_WATER_CONTENT = 1e-3
GOLOVIN_B = 1500.  # native kernel: b * (volume_1 + volume_2)
EXPECTED_BINARY = {'scale-les': '88d9bb6f7040a9fc08afaf54cba81ae663390cd0a572ddcbd252879011eccbc2',
                   'scale-les_init': 'de580b9dbde137bcf34499ff91270850f49fd2bf0102714fb5324536dbd226c1'}


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def set_value(text, key, value):
    """Change exactly one existing scalar/list in our generated namelist."""
    encoded = column.namelist('x', **{key: value}).splitlines()[1]
    text, count = re.subn(r'^\s*'+re.escape(key)+r'\s*=.*$', encoded, text, flags=re.M|re.I)
    if count != 1:
        raise ValueError(f'Expected one namelist value for {key}, got {count}')
    return text


def initial_population(sd, kind, seed):
    """Equal-multiplicity midpoint quantiles; exact sampled mean mass in each cell.

    Quantile rescaling is recorded as a finite sampling approximation. Refinement
    changes that approximation; it never changes the specified bulk ice water.
    """
    if kind not in ('ice', 'golovin') or sd < 4 or sd != int(sd) or seed < 0:
        raise ValueError('Invalid initial population')
    mean = ICE_MEAN_MASS if kind == 'ice' else LIQUID_WATER_CONTENT / LIQUID_NUMBER
    concentration = ICE_WATER_CONTENT/mean if kind == 'ice' else LIQUID_NUMBER
    multiplicity = round(concentration * CELL_VOLUME / sd)
    if multiplicity < 1:
        raise ValueError('Insufficient physical particles for this sampling')
    concentration_actual = multiplicity*sd/CELL_VOLUME
    water = ICE_WATER_CONTENT if kind == 'ice' else LIQUID_WATER_CONTENT
    mass = -np.log1p(-(np.arange(sd)+.5)/sd)
    mass *= water/concentration_actual/mass.mean()
    rng = np.random.default_rng(seed)
    ix = np.array([(k,i,j) for k in range(GRID[0]) for i in range(GRID[1])
                   for j in range(GRID[2]) for _ in range(sd)])
    positions = (ix + rng.uniform(.05,.95,size=ix.shape))*DX
    masses = np.tile(mass, NCELL)
    radius = (3*masses/(4*math.pi*(ICE_DENSITY if kind=='ice' else 1000.)))**(1/3)
    n = len(masses)
    return dict(n=np.full(n,multiplicity,dtype=np.int64), rk=positions[:,0]/DX+2,
                x=positions[:,1],y=positions[:,2],z=positions[:,0],
                r=np.zeros(n) if kind=='ice' else radius,
                re=radius if kind=='ice' else np.zeros(n),
                rp=radius if kind=='ice' else np.zeros(n),
                rho=np.full(n,ICE_DENSITY) if kind=='ice' else np.zeros(n),
                phase=np.full(n,10 if kind=='ice' else 1,dtype=np.int16),
                monomers=np.full(n,1 if kind=='ice' else 0,dtype=np.int32),
                mass=masses)


def write_restart(path, pop):
    """The pinned gfortran build uses big-endian sequential records, 4-byte markers."""
    n=len(pop['n'])
    zeros=np.zeros(n)
    with Path(path).open('xb') as out:
        def record(data):
            out.write(struct.pack('>i',len(data))); out.write(data)
            out.write(struct.pack('>i',len(data)))
        record(struct.pack('>d5i',0.,8,8,n,1,1))
        fields=[(pop['n'],'>i8')]+[(pop[k],'>f8') for k in ('rk','x','y','z','r')]
        fields += [(zeros,'>f8')]*4  # u, v, fall speed, aerosol mass (one species)
        fields += [(pop['phase'],'>i2')]
        fields += [(pop[k],'>f8') for k in ('re','rp','rho')]
        fields += [(np.full(n,230.),'>f8'),(zeros,'>f8'),(pop['monomers'],'>i4')]
        for values,dtype in fields:
            record(np.asarray(values,dtype=dtype).tobytes())


def configure(cache, name, kind='ice', sd=512, dt=5., duration=3600., seed=0,
              lunar=False, collisions=True, observer=False):
    if not all(math.isfinite(x) and x > 0 for x in (dt,duration)):
        raise ValueError('Positive finite duration and step required')
    if kind not in ('ice','golovin') or duration % 300 or 300/dt != round(300/dt):
        raise ValueError('Use ice/golovin and positive steps dividing 300 s')
    if lunar and kind=='golovin':
        raise ValueError('The analytic kernel is gravity-independent')
    observer_record=None
    expected_binary=EXPECTED_BINARY
    if observer:
        from atmosphere.snow.scale_collision_observer import SOURCE_SHA256
        manifest=Path(cache)/'observer/observer.json'
        observer_record=json.loads(manifest.read_text())
        if observer_record['original_sha256']!=SOURCE_SHA256 or observer_record['build']['exit_code']!=0:
            raise ValueError('Unidentified observer build')
        expected_binary=observer_record['build']['executable_sha256']
    pop=initial_population(sd,kind,seed)
    d=column.configure(cache,name,lunar=lunar,duration=duration,dt=dt,sd=sd,seed=seed,
                       move=False,deposition=False,collisions=collisions)
    for filename in ('init.conf','run.conf'):
        text=(d/filename).read_text()
        text=set_value(text,'KMAX',GRID[0]+1)
        text=set_value(text,'BUFFER_DZ',0.)
        if filename=='run.conf':
            for k,v in dict(docondensation=False,domeltfreeze=not collisions,
                            sdm_colkrnl=0 if kind=='golovin' else 2,
                            sdm_dtcmph=[dt]*5,sdm_zupper=GRID[0]*DX,
                            sdm_extbuf=0,sdm_dmpitva=300.,sdm_dmpitvb=300.,
                            HISTORY_DEFAULT_TINTERVAL=300.).items():
                text=set_value(text,k,v)
            text=text.replace('&PARAM_ATMOS_PHY_MP_SDM\n',
                              "&PARAM_ATMOS_PHY_MP_SDM\n SD_IN_BASENAME = 'particles',\n")
        (d/filename).write_text(text)
    write_restart(d/'particles.pe000000',pop)
    settings=dict(name=name,kind=kind,sd_per_cell=sd,dt_s=dt,duration_s=duration,seed=seed,
        gravity_m_s2=MOON_SURFACE_GRAVITY if lunar else STANDARD_GRAVITY,
        collisions=collisions,grid_zxy=list(GRID),cell_volume_m3=CELL_VOLUME,
        air_temperature_K=268.15,air_pressure_Pa=100000.,vapor_mass_fraction=0.,
        initial_sampling='Identical exponential midpoint mass quantiles in every cell, rescaled to exact bulk water; equal integer multiplicity; seeded within-cell positions.',
        initial_mean_mass_kg=float(np.average(pop['mass'],weights=pop['n'])),
        initial_number_m3=float(pop['n'].sum()/(NCELL*CELL_VOLUME)),
        initial_water_kg_m3=float(np.dot(pop['n'].astype(float),pop['mass'])/(NCELL*CELL_VOLUME)),
        expected_binary_sha256=expected_binary,
        observer=observer,observer_build=observer_record,
        source_product_sha256=column.digest(Path(cache)/'source_manifest.json'))
    save(d/'settings.json',settings)
    save(d/'inputs.json',{p.name:column.digest(p) for p in d.iterdir() if p.is_file()})
    return d


def execute(d, timeout=600):
    d=Path(d).resolve()
    cache=d.parents[1]
    settings=json.loads((d/'settings.json').read_text())
    binaries=cache/('observer/work' if settings.get('observer') else 'work')/CASE_REL
    for name,expected in settings['expected_binary_sha256'].items():
        if column.digest(binaries/name)!=expected:
            raise ValueError(f'Unidentified binary: {name}')
    env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1',
         'GFORTRAN_UNBUFFERED_ALL':'y',
         'LD_LIBRARY_PATH':str(cache/'deps/usr/lib/x86_64-linux-gnu')+os.pathsep+os.environ.get('LD_LIBRARY_PATH','')}
    executions=[]
    for name,conf in [('scale-les_init','init.conf'),('scale-les','run.conf')]:
        start=time.monotonic()
        with (d/(name+'.log')).open('x') as log:
            proc=subprocess.Popen([str(binaries/name),conf],cwd=d,env=env,
                                  stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            try:
                rc=proc.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid,signal.SIGKILL); proc.wait(); rc='timeout'
        executions.append(dict(executable=name,exit_code=rc,sha256=column.digest(binaries/name),
                               elapsed_s=time.monotonic()-start))
        save(d/'execution.json',executions)
        print(d.name,name,rc,round(time.monotonic()-start,2),flush=True)
        if rc!=0:
            raise RuntimeError(f'Failed execution, logs retained: {d}')
        if name=='scale-les_init':
            # A box has uniform imposed air; no hydrostatic column is inferred.
            init=next(d.glob('init_*.nc'))
            before=column.digest(init)
            with column._dataset(init) as nc:
                old={k:np.array(nc[k][:]) for k in ('DENS','RHOT','QV')}
            from netCDF4 import Dataset
            with Dataset(init,'r+') as nc:
                nc['DENS'][:]=np.full(nc['DENS'].shape,100000./(column.SOURCE_CONSTANTS['Rdry']*268.15))
                nc['RHOT'][:]=np.full(nc['RHOT'].shape,100000./column.SOURCE_CONSTANTS['Rdry'])
                for k in ('QV','QC','QR','QI','QS','QG','MOMX','MOMY','MOMZ'):
                    nc[k][:]=np.zeros(nc[k].shape)
            save(d/'air_initialization.json',dict(before_sha256=before,after_sha256=column.digest(init),
                 imposed_temperature_K=268.15,imposed_pressure_Pa=100000.,
                 original_ranges={k:[float(v.min()),float(v.max())] for k,v in old.items()}))
    return executions


def sample(path,kind):
    with column._dataset(path) as nc:
        v={k:np.array(x[:]) for k,x in nc.variables.items()}
    keep=(v['sd_n']>0)&(v['sd_z']>=0)
    phase=10 if kind=='ice' else 1
    if not np.all(v['sd_liqice'][keep]==phase):
        raise ValueError('Unexpected phase in collision-only test')
    n=v['sd_n'][keep].astype(float)
    re,rp,rho=(v[k][keep] for k in ('sdi_re','sdi_rp','sdi_rho'))
    m=4*math.pi/3*(re**2*rp*rho if kind=='ice' else v['sd_r'][keep]**3*1000.)
    diameter=2*np.maximum(re,rp) if kind=='ice' else 2*v['sd_r'][keep]
    if not np.isfinite(m).all() or (m<=0).any():
        raise ValueError('Invalid particle mass')
    cell=(np.floor(v['sd_z'][keep]/DX).astype(int)*GRID[1]
          +np.floor(v['sd_x'][keep]/DX).astype(int))*GRID[2]+np.floor(v['sd_y'][keep]/DX).astype(int)
    if (cell<0).any() or (cell>=NCELL).any():
        raise ValueError('Particle outside box')
    moments=np.array([np.bincount(cell,weights=n*m**order,minlength=NCELL)/CELL_VOLUME for order in range(3)])
    agg=v['sdi_nmono'][keep]>1
    result=dict(file=path.name,sha256=column.digest(path),active_superparticles=int(keep.sum()),
        moments_per_m3=moments.mean(axis=1).tolist(),
        cell_moments_per_m3=moments.T.tolist(),
        mass_weighted_D_quantiles_m=column.weighted_quantile(diameter,n*m),
        number_weighted_D_quantiles_m=column.weighted_quantile(diameter,n),
        mean_D_m=float(np.dot(n,diameter)/n.sum()),max_D_m=float(diameter.max()),
        aggregate_mass_fraction=float(np.dot(n[agg],m[agg])/np.dot(n,m)) if kind=='ice' else None,
        monomer_units_m3=float(np.dot(n,v['sdi_nmono'][keep])/(CELL_VOLUME*NCELL)) if kind=='ice' else None,
        mass_weighted_density_kg_m3=column.weighted_quantile(rho,n*m) if kind=='ice' else None,
        mass_weighted_aspect=column.weighted_quantile(rp/re,n*m) if kind=='ice' else None)
    return result, (m,n,cell)


def golovin_moments(initial,time_s):
    n,m1,m2=initial
    rate=GOLOVIN_B*m1/1000.
    return [n*math.exp(-rate*time_s),m1,m2*math.exp(2*rate*time_s)]


def golovin_mass_cdf(mass,initial_number,water,time_s,points=24000):
    """Safranov-Golovin solution for a continuous exponential initial spectrum.

    Scaled I1 avoids overflow in the large-size tail. The integration is performed
    in dimensionless mass, with an independently checked normalization; no
    numerical renormalization conceals mass lost outside the integration domain.
    """
    from scipy.special import ive
    from scipy.integrate import cumulative_trapezoid
    mass=np.asarray(mass,dtype=float)
    if (mass<0).any() or not np.isfinite(mass).all() or time_s<0 or initial_number<=0 or water<=0:
        raise ValueError('Invalid analytic distribution inputs')
    tau=-math.expm1(-GOLOVIN_B*water/1000.*time_s)
    mean=water/initial_number
    if tau==0:
        y=mass/mean
        return 1-(1+y)*np.exp(-y),1.
    root=math.sqrt(tau)
    y=np.geomspace(1e-10,max(100.,50./(1-root)**2),points)
    density=(1-tau)/root*ive(1,2*y*root)*np.exp(-(1-root)**2*y)
    cdf=cumulative_trapezoid(density,y,initial=0.)
    if abs(cdf[-1]-1)>2e-6:
        raise ValueError(f'Analytic integration did not close mass: {cdf[-1]}')
    return np.interp(mass/mean,y,cdf,left=0.,right=cdf[-1]),float(cdf[-1])


def distribution_distance(m,n,expected_cdf):
    """Supremum CDF difference includes both sides of each empirical jump."""
    order=np.argsort(m)
    m,n=m[order],n[order]
    # Coincident masses must be one CDF jump, not order-dependent smaller steps.
    unique,ix=np.unique(m,return_index=True)
    weight=np.add.reduceat(m*n,ix)
    cumulative=np.cumsum(weight)/weight.sum()
    exact=expected_cdf(unique)
    return float(max(np.max(np.abs(cumulative-exact)),
                     np.max(np.abs(np.r_[0.,cumulative[:-1]]-exact))))


def analyze(d):
    d=Path(d)
    settings=json.loads((d/'settings.json').read_text())
    inputs=json.loads((d/'inputs.json').read_text())
    for name,expected_hash in inputs.items():
        if Path(name).name!=name or column.digest(d/name)!=expected_hash:
            raise ValueError(f'Changed or invalid run input: {name}')
    runs=json.loads((d/'execution.json').read_text())
    if len(runs)!=2 or any(r['exit_code']!=0 or r['sha256']!=settings['expected_binary_sha256'][r['executable']] for r in runs):
        raise ValueError('Incomplete or unidentified execution')
    paths=sorted(d.glob('SD_output_NetCDF_*.pe000000'))
    samples=[]
    for p in paths:
        row,(m,n,cell)=sample(p,settings['kind'])
        # Native file time is before microphysics; these snapshots have exact
        # elapsed collision time from the initial particle snapshot.
        matches=re.findall(r'[0-9]+\.[0-9]+',p.name)
        if len(matches)!=1:
            raise ValueError(f'Unknown particle snapshot name: {p.name}')
        row['time_s']=float(matches[0])
        if settings['kind']=='golovin':
            number=settings['initial_number_m3'];water=settings['initial_water_kg_m3']
            def cdf(mass):
                return golovin_mass_cdf(mass,number,water,row['time_s'])[0]
            row['analytic_mass_cdf_max_error']=distribution_distance(m,n,cdf)
        # Dimensionless fixed bins record distribution shape for cross-run
        # comparisons; underflow and overflow remain visible in the CDF ends.
        mean=settings['initial_mean_mass_kg']
        edges=mean*np.geomspace(1e-4,1e7,177)
        order=np.argsort(m);sm=m[order];sn=n[order]
        idx=np.searchsorted(sm,edges,side='right')
        cdf=np.r_[0.,np.cumsum(sm*sn)]/np.dot(sm,sn)
        row['mass_cdf']=cdf[idx].tolist()
        row['mass_cdf_edges_kg']=edges.tolist()
        samples.append(row)
    expected=np.arange(0,settings['duration_s']+1,300.)
    if len(samples)!=len(expected) or not np.array_equal([s['time_s'] for s in samples],expected):
        raise ValueError(f'Incomplete snapshot cadence: {[s["file"] for s in samples]}')
    initial=samples[0]['moments_per_m3']
    np.testing.assert_allclose(initial[:2],
        [settings['initial_number_m3'],settings['initial_water_kg_m3']],rtol=1e-12,atol=0.)
    m0=np.array(samples[0]['cell_moments_per_m3'])[:,1]
    mass_error=max(float(np.max(np.abs(np.array(s['cell_moments_per_m3'])[:,1]/m0-1))) for s in samples)
    with column._dataset(next(d.glob('history*.nc'))) as nc:
        history={k:np.array(nc[k][:]) for k in ('T','PRES','QV','PREC_sd','QHYD_sd','DENS')}
    air=dict(temperature_range_K=[float(history['T'].min()),float(history['T'].max())],
             pressure_range_Pa=[float(history['PRES'].min()),float(history['PRES'].max())],
             max_abs_vapor=float(np.abs(history['QV']).max()),max_precipitation_m=float(history['PREC_sd'].max()),
             max_bulk_water_relative_error=float(np.max(np.abs((history['QHYD_sd']*history['DENS'])[:, :GRID[0]]/initial[1]-1))))
    if not all(np.isfinite(v).all() for v in history.values()):
        raise ValueError('Non-finite atmospheric history')
    if settings['kind']=='golovin':
        for s in samples:
            exact=np.array(golovin_moments(initial,s['time_s']))
            cells=np.array(s['cell_moments_per_m3'])
            s['analytic_moments_per_m3']=exact.tolist()
            s['relative_moment_error']=(cells.mean(axis=0)/exact-1).tolist()
            s['cell_standard_error_relative']=(cells.std(axis=0,ddof=1)/math.sqrt(NCELL)/exact).tolist()
    result=dict(settings=settings,executions=runs,samples=samples,air=air,
        max_cell_mass_relative_error=mass_error,
        max_monomer_units_relative_error=max(abs(s['monomer_units_m3']/samples[0]['monomer_units_m3']-1) for s in samples) if settings['kind']=='ice' else None,
        input_hashes=json.loads((d/'inputs.json').read_text()),
        air_initialization=json.loads((d/'air_initialization.json').read_text()),
        producer_sha256=column.digest(Path(__file__)))
    result['mass_cdf_edges_kg']=samples[0]['mass_cdf_edges_kg']
    for row in samples:
        del row['mass_cdf_edges_kg']
    if settings.get('observer'):
        counter=np.loadtxt(d/'collision_observer.txt',ndmin=2)
        expected_calls=round(settings['duration_s']/settings['dt_s'])+1
        if counter.shape!=(expected_calls,8) or not np.array_equal(counter[:,0],np.arange(1,expected_calls+1)):
            raise ValueError('Incomplete observer trace')
        # The final native call occurs after the last saved particle snapshot.
        trace=counter[:-1]
        attempted,accepted,deficit=trace[:,2:5].sum(axis=0)
        actual_loss=(samples[0]['moments_per_m3'][0]-samples[-1]['moments_per_m3'][0])*NCELL*CELL_VOLUME
        if abs(accepted-actual_loss)>max(1.,1e-10*accepted):
            raise ValueError('Collision count does not match physical number loss')
        result['observer']=dict(sha256=column.digest(d/'collision_observer.txt'),
            calls_through_final_particle_snapshot=len(trace),excluded_final_call=counter[-1].tolist(),
            attempted_physical_collisions=float(attempted),accepted_physical_collisions=float(accepted),
            omitted_physical_collisions=float(deficit),omitted_fraction=float(deficit/attempted) if attempted else 0.,
            capped_pair_count=int(trace[:,5].sum()),max_probability=float(trace[:,6].max()),
            max_probability_over_capacity=float(trace[:,7].max()),
            number_loss_accounting_absolute_error=float(accepted-actual_loss))
    save(d/'aggregation_analysis.json',result)
    return result


def compare_distributions(cache,base,other,results):
    left=results[base]['samples'][-1];right=results[other]['samples'][-1]
    _,(lm,ln,_)=sample(Path(cache)/'runs'/base/left['file'],results[base]['settings']['kind'])
    _,(rm,rn,_)=sample(Path(cache)/'runs'/other/right['file'],results[other]['settings']['kind'])
    edges=np.unique(np.r_[lm,rm])
    def cdf(m,n):
        order=np.argsort(m);mass=m[order];weights=mass*n[order]
        return (np.r_[0.,np.cumsum(weights)]/weights.sum())[np.searchsorted(mass,edges,side='right')]
    lmom=np.array(left['moments_per_m3']);rmom=np.array(right['moments_per_m3'])
    return dict(base=base,other=other,
        relative_moment_change=(rmom/lmom-1).tolist(),
        relative_D_quantile_change=(np.array(right['mass_weighted_D_quantiles_m'])/left['mass_weighted_D_quantiles_m']-1).tolist(),
        mass_cdf_max_difference=float(np.max(np.abs(cdf(lm,ln)-cdf(rm,rn)))),
        aggregate_mass_fraction_difference=(right['aggregate_mass_fraction']-left['aggregate_mass_fraction']) if left['aggregate_mass_fraction'] is not None else None)


def summarize(cache,output):
    cache=Path(cache)
    declarations=[cache/'aggregation_suite_02.json',cache/'aggregation_followup_suite.json']
    names=[name for path in declarations for name,_ in json.loads(path.read_text())['cases']]
    if len(names)!=len(set(names)):
        raise ValueError('Repeated declared case')
    results={name:analyze(cache/'runs'/name) for name in names}
    pairs=[('ice_128','ice_512_02'),('ice_512_02','ice_2048'),
           ('ice_512_02','ice_half_step'),('ice_512_02','ice_512_seed1_02'),
           ('ice_four_hour_512','ice_four_hour'),('ice_four_hour','ice_four_hour_half'),
           ('ice_2048','lunar_ice_2048'),('golovin_8192_dt1','golovin_8192_half'),
           ('golovin_8192_dt1','golovin_8192_seed1')]
    references=REPO/'atmosphere/snow/scale_aggregation_sources.json'
    parity_path=cache/'aggregation_observer_parity.json'
    result=dict(schema='terluna.atmosphere.scale-aggregation/1',
        evidence='Native-solver numerical controls: analytic liquid coagulation and aggregation-only ice boxes. No observed snowstorm reproduction, coupled column convergence or giant-aggregate survival validation.',
        reading_rule='Read sizes only for the declared seeded finite-water boxes. The ice protocol borrows published initial mass and air conditions while retaining different SCALE efficiencies, drag and evolving geometry. Analytic convergence checks use continuous exponential reference distributions; initial quantile discretization and Monte Carlo sampling remain visible.',
        producer=dict(path='atmosphere/snow/scale_aggregation.py',sha256=column.digest(Path(__file__)),
            dependencies={name:column.digest(REPO/name) for name in
                ('atmosphere/snow/scale_column.py','atmosphere/snow/scale_build.py',
                 'atmosphere/snow/scale_collision_observer.py','shared/constants.json')}),
        references=dict(path=str(references.relative_to(REPO)),sha256=column.digest(references),
                        content=json.loads(references.read_text())),
        source_manifest_sha256=column.digest(cache/'source_manifest.json'),
        declarations={p.name:dict(sha256=column.digest(p),content=json.loads(p.read_text())) for p in declarations},
        cases=results,comparisons=[compare_distributions(cache,*pair,results) for pair in pairs],
        observer_parity=dict(sha256=column.digest(parity_path),content=json.loads(parity_path.read_text())),
        incomplete_attempts=['aggregation_smoke: vertical buffer exceeded small domain',
            'aggregation_smoke_02: restart count mismatched native empty-top-cell restriction',
            'aggregation_smoke_03: native run succeeded; initial analysis glob omitted extensionless particle outputs, corrected without rerunning',
            'ice_512 and ice_512_seed1: threaded NetCDF initialization crashed the Python driver; separate-process retries retained these directories'],
        unperformed=['Exact SHIPS size-distribution or Welss column reproduction with SCALE',
            'Adaptive collision stepping or spectral sampling repair',
            'Energy-budget closure and aggregation in a coupled lunar column',
            'Air motion, mechanical fracture and storm optical forecast'])
    save(output,result)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('name'); p.add_argument('--cache',type=Path,default=CACHE)
    p.add_argument('--kind',choices=['ice','golovin'],default='ice')
    p.add_argument('--sd',type=int,default=512); p.add_argument('--dt',type=float,default=5.)
    p.add_argument('--duration',type=float,default=3600.); p.add_argument('--seed',type=int,default=0)
    p.add_argument('--lunar',action='store_true'); p.add_argument('--no-collisions',action='store_true')
    p.add_argument('--observer',action='store_true')
    p.add_argument('--summarize-to',type=Path)
    p.add_argument('--analyze-only',action='store_true'); p.add_argument('--timeout',type=int,default=600)
    a=p.parse_args()
    if a.summarize_to:
        result=summarize(a.cache,a.summarize_to)
        print(f'Saved {len(result["cases"])} declared controls to {a.summarize_to}')
        return
    d=a.cache/'runs'/a.name
    if not a.analyze_only:
        d=configure(a.cache,a.name,a.kind,a.sd,a.dt,a.duration,a.seed,a.lunar,not a.no_collisions,a.observer)
        execute(d,a.timeout)
    result=analyze(d)
    print(json.dumps({k:result[k] for k in ('air','max_cell_mass_relative_error','max_monomer_units_relative_error')},indent=2))


if __name__=='__main__':
    main()
