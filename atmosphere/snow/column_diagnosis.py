"""Observe and ablate the frozen full-column calculation without changing it.

The diagnostic driver includes event_column.cpp and reuses its physical routines.
Input preparation and result conversion mirror the frozen interface so every
baseline scientific field can be compared at fixed seed. New files are explicitly
bound to the original source and do not invalidate the original products.
"""
from __future__ import annotations
import ctypes as ct
import hashlib
import math
import os
from pathlib import Path
import subprocess
import tempfile
import time
import numpy as np
from atmosphere.snow import event_column as baseline
from atmosphere.snow.event_column import ag,sf,GEOMETRY,parameters,diameter_queries

ROOT=baseline.ROOT
SOURCE=Path(__file__).with_suffix('.cpp')
SCHEMA='terluna.atmosphere.snow-column-diagnosis/1'
PTR=baseline.PTR
_LIBRARY=None
_BUILD=None
DIAGNOSTIC_COLUMNS=[
 'sample_time_s','empty_time_s','single_particle_time_s','vapour_limited_time_s',
 'unrestricted_demand_kg','admitted_demand_kg','empty_unused_supply_kg',
 'occupied_unused_supply_kg','number_time_s','number_squared_time_s',
 'effective_growth_number_time_s','largest_growth_share_time_s',
 'sticking_events','chipping_events','emitted_fragments']


def identity(binding):
    result=baseline.identity(binding)
    for path in (Path(__file__),SOURCE):
        result['files'][str(path.resolve().relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def library():
    global _LIBRARY,_BUILD
    if _LIBRARY is not None:return _LIBRARY
    compiler=subprocess.run(['g++','--version'],check=True,capture_output=True,text=True).stdout
    flags=['-std=c++17','-O3','-fPIC','-shared','-Wall','-Wextra','-Werror']
    sources={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in (SOURCE,baseline.SOURCE)}
    digest=hashlib.sha256((repr(sources)+compiler+repr(flags)).encode()).hexdigest()
    directory=Path(tempfile.gettempdir())/'terluna-snow-diagnosis-builds'/digest
    directory.mkdir(parents=True,exist_ok=True);target=directory/'column_diagnosis.so'
    if not target.exists():
        temporary=directory/f'build-{os.getpid()}.so'
        subprocess.run(['g++',*flags,str(SOURCE),'-o',str(temporary)],check=True,capture_output=True,text=True)
        temporary.replace(target)
    lib=ct.CDLL(str(target))
    lib.snow_column_diagnosis.argtypes=[PTR,PTR,PTR,PTR,ct.c_int,ct.c_uint64,PTR,PTR,ct.c_char_p]
    lib.snow_column_diagnosis.restype=ct.c_int
    _LIBRARY=lib
    _BUILD=dict(compiler=compiler.splitlines()[0],flags=flags,sources=sources,
                binary_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),build_key=digest)
    return lib


def column(world,*,top_c=-15.,number_per_m3=1e5,snowfall_mm_h=3.,layers=30,
           area_m2=1e-4,max_step_s=30.,fill=1.5,gather=1.5,breakup=True,seed=0,
           wall_limit_s=1800.,collisions=True,growth=True,settling=True,settling_scheme='upwind'):
    if settling_scheme not in ('upwind','ballistic'):raise ValueError('Unknown settling scheme')
    vals=(top_c,number_per_m3,snowfall_mm_h,area_m2,max_step_s,fill,gather,wall_limit_s)
    if not all(math.isfinite(x) for x in vals):raise ValueError('Finite parameters required')
    if top_c>=0 or number_per_m3<=0 or snowfall_mm_h<0 or min(area_m2,max_step_s,gather,wall_limit_s)<=0 or fill<0:
        raise ValueError('Invalid column forcing or duration')
    if not isinstance(layers,int) or layers<1:raise ValueError('Positive integer layer count required')
    if not isinstance(seed,int) or not 0<=seed<2**64:raise ValueError('Unsigned 64-bit seed required')
    geometry=GEOMETRY
    depth=world.height_of(sf.T0+top_c)-world.z_freeze;dz=depth/layers
    height=world.z_freeze+(np.arange(layers)+.5)*dz
    states=[world.state(float(z)) for z in height];p,air=parameters(world,states)
    seed_d=math.sqrt(ag.START_M[0]*ag.START_M[1])
    u=float(ag.speeds(np.array([seed_d]),np.array([geometry.mass(seed_d)]),world.gravity,
                      world.state(world.z_freeze+depth),geometry.area_ratio)[0][0])
    crossing=depth/u;gather_start=fill*crossing;end=(fill+gather)*crossing
    C=geometry.area_ratio*math.pi/4
    cfg=np.array([layers,area_m2,dz,max_step_s,gather_start,end,snowfall_mm_h/3600/depth,
                  number_per_m3*u,*ag.START_M,C,breakup,wall_limit_s,collisions,growth,settling,settling_scheme=='upwind'],float)
    diam=diameter_queries();query=np.r_[geometry.mass(diam),np.inf]
    out=np.zeros(24+4*query.size+4*layers+16);err=ct.create_string_buffer(1024)
    start=time.monotonic()
    trace=np.zeros((5,layers,len(DIAGNOSTIC_COLUMNS)))
    code=library().snow_column_diagnosis(p,air,cfg,query,query.size,seed,out,trace,err)
    if code:raise RuntimeError(err.value.decode())
    mass_names=('injected','deposited','freezing_exit','remaining','fractional_fragment_remaining','residual')
    ledger={k:float(v/area_m2) for k,v in zip(mass_names,out[:6])}
    ledger['relative_residual']=float(out[5]/(out[0]+out[1])) if out[0]+out[1] else 0.
    nres=out[6]+out[11]-out[7]-out[8]-out[9]-out[10]
    number=dict(zip(('injected','freezing_exit','remaining','fractional_fragment_remaining',
                     'coagulation_loss','chipping_production'),map(float,out[6:12]/area_m2)))
    number['residual']=float(nres/area_m2)
    number['relative_residual']=float(nres/(out[6]+out[11])) if out[6]+out[11] else 0.
    cumulative=out[24:24+4*query.size].reshape(-1,4)/area_m2
    profiles=out[24+4*query.size:24+4*query.size+4*layers].reshape(-1,4)/(area_m2*dz*gather*crossing)
    blocks=out[-16:].reshape(4,4)/area_m2
    total=cumulative[-1]
    summary=dict(arriving_number_m2=float(total[0]),snowfall_mm_h=float(total[1]*3600/(gather*crossing)),
                 mass_d50_m=float(out[21]) if total[0] else None,
                 mass_d90_m=float(out[22]) if total[0] else None,mass_d99_m=float(out[23]) if total[0] else None,
                 maximum_diameter_anywhere_m=float(out[18]),
                 arriving_mass_share_above={str(q):float(1-cumulative[np.searchsorted(diam,q),1]/total[1]) if total[1] else 0.
                                           for q in (.005,.01,.15,.6,1.,10.)})
    return dict(schema=SCHEMA,process_diagnostics=dict(columns=DIAGNOSTIC_COLUMNS,
        blocks=trace.tolist(),phase_names=['filling','gather_1','gather_2','gather_3','gather_4'],
        reading_rule='Time integrals sample both Heun stages at growth substeps. They are split-step diagnostics, not exact residence integrals. Demand and supply fields are extensive kg in the represented layer. Disabled growth has zero demand and supply; occupancy remains measured. Event counts are actual per-layer increments.'),settings=dict(top_c=top_c,number_per_m3=number_per_m3,
        snowfall_mm_h=snowfall_mm_h,layers=layers,area_m2=area_m2,layer_volume_m3=area_m2*dz,
        max_step_s=max_step_s,fill=fill,gather=gather,breakup=breakup,seed=seed,wall_limit_s=wall_limit_s,
        collisions=collisions,growth=growth,settling=settling,settling_scheme=settling_scheme,geometry=geometry.id,
        injection_reference_velocity_m_s=u,injected_number_flux_m2_s=number_per_m3*u,
        total_time_s=end,gathering_time_s=gather*crossing,depth_m=depth),
        summary=summary,mass_ledger_kg_m2=ledger,number_ledger_m2=number,
        query_diameter_m=[*map(float,diam),None],arriving_cumulative_moments=cumulative.tolist(),
        moment_columns=['number_m2','mass_kg_m2','second_mass_moment_kg2_m2','projected_area_m2_m2'],
        layer_height_m=height.tolist(),time_mean_layer_moments=profiles.tolist(),
        layer_moment_columns=['number_m3','mass_kg_m3','second_mass_moment_kg2_m3','projected_area_m2_m3'],
        arrival_time_blocks=dict(seconds_per_block=gather*crossing/4,moments=blocks.tolist()),
        diagnostics=dict(proposals=int(out[12]),chipping_events=int(out[13]),emitted_fragments=int(out[14]),
            peak_particles=int(out[15]),steps=int(out[16]),final_time_s=float(out[17]),maximum_thinning_probability=float(out[19]),
            maximum_fractional_fragment_mass_bound_kg_m2=float(layers*p[4]/area_m2),
            omitted_fractional_fragment_initial_demand_bound_kg_m2=float(out[20]/area_m2),
            wall_seconds=time.monotonic()-start,build=dict(_BUILD)),
        evidence='Continuous-mass finite-volume particle calculation with inherited geometry, drag, sticking, '
                 'chipping and liquid-saturated growth. Collision events are exact at fixed layer properties; '
                 'volume, vertical mixing, split-time, fractional-fragment allocation and sampling errors remain distinct. '
                 'No observed-storm calibration or giant-aggregate survival law.',
        reading_rule='The final null diameter means +infinity and retains all arrivals. Moment arrays are unnormalized '
                     'cumulative extensive quantities per square metre. Quantiles are per realization. Fractional '
                     'fragments retain water and number in local allocation remainders, participate after allocation '
                     'of a whole particle, and are bounded separately; no fragment water is discarded. '
                     'The demand bound excludes delayed-allocation feedback. Physical tail predictions require convergence '
                     'and suitable geometry, fracture and cloud exposure.')
