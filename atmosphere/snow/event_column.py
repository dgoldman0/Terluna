"""Continuous-mass event reference for closed aggregation and full snow columns.

The C++ kernel is compiled locally into a content-addressed temporary directory.
It implements exact thinning for collisions at fixed layer properties, continuous
mass growth, and either matched upwind or ballistic settling. Finite volume, operator splitting, layer
mixing, and fractional-fragment allocation are independent approximation axes.
The prescribed physical closures remain empirically unvalidated for giant snow.
"""
from __future__ import annotations
import ctypes as ct
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
import time
import numpy as np
from atmosphere.snow import aggregation as ag, flakes as sf
from atmosphere.snow.habits import GEOMETRIES
from shared.provenance import constants_used

ROOT = Path(__file__).resolve().parents[2]
SOURCE = Path(__file__).with_suffix('.cpp')
SCHEMA = 'terluna.atmosphere.snow-event-column/1'
GEOMETRY = GEOMETRIES['bf95_dmax_dendrites_area']
FILES = ('atmosphere/snow/event_column.py','atmosphere/snow/event_column.cpp',
         'atmosphere/snow/aggregation.py','atmosphere/snow/flakes.py',
         'atmosphere/snow/habits.py','atmosphere/rain/drops.py','atmosphere/rain/targets.py',
         'atmosphere/radiative_convective/thermodynamics.py')
PTR = np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS')
_LIBRARY = None
_BUILD = None


def identity(binding):
    sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
    inputs={p:sha(ROOT/p) for p in binding.get('inputs',{})}
    if inputs != binding.get('inputs',{}):raise RuntimeError('Column input changed')
    return dict(files={p:sha(ROOT/p) for p in FILES},inputs=inputs,
                constants=constants_used([p for p in FILES if p.endswith('.py')]))


def library():
    global _LIBRARY,_BUILD
    if _LIBRARY is not None:return _LIBRARY
    compiler=subprocess.run(['g++','--version'],check=True,capture_output=True,text=True).stdout
    flags=['-std=c++17','-O3','-fPIC','-shared','-Wall','-Wextra','-Werror']
    digest=hashlib.sha256(SOURCE.read_bytes()+compiler.encode()+repr(flags).encode()).hexdigest()
    directory=Path(tempfile.gettempdir())/'terluna-snow-event-builds'/digest
    directory.mkdir(parents=True,exist_ok=True)
    target=directory/'event_column.so'
    if not target.exists():
        temporary=directory/f'build-{os.getpid()}.so'
        subprocess.run(['g++',*flags,str(SOURCE),'-o',str(temporary)],check=True,capture_output=True,text=True)
        temporary.replace(target)
    lib=ct.CDLL(str(target))
    lib.snow_properties.argtypes=[PTR,PTR,PTR,ct.c_int,PTR,ct.c_char_p]
    lib.snow_fragment_loss.argtypes=[PTR,PTR,PTR,PTR,ct.c_int,PTR,ct.c_char_p]
    lib.snow_box.argtypes=[PTR,PTR,PTR,ct.c_int,ct.c_double,PTR,ct.c_int,PTR,ct.c_int,
                          ct.c_uint64,ct.c_int,ct.c_double,ct.c_double,PTR,PTR,PTR,ct.c_char_p]
    lib.snow_column.argtypes=[PTR,PTR,PTR,PTR,ct.c_int,ct.c_uint64,PTR,ct.c_char_p]
    for name in ('snow_properties','snow_fragment_loss','snow_box','snow_column'):getattr(lib,name).restype=ct.c_int
    _LIBRARY=lib
    _BUILD=dict(compiler=compiler.splitlines()[0],flags=flags,
                source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                binary_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),build_key=digest)
    return lib


def parameters(world,states,geometry=GEOMETRY):
    p=np.array([geometry.mass_coefficient,geometry.mass_exponent,sf.DELTA0,sf.C0,
                geometry.mass(ag.FRAGMENT_M),*ag.BREAKUP['size_m'],ag.BREAKUP['a_scale'],
                ag.BREAKUP['a_size'],ag.BREAKUP['c'],ag.BREAKUP['gamma'],ag.BREAKUP['most']],float)
    rows=[]
    for s in states:
        ls=ag.drops.latent_heat(sf.T0)+sf.LATENT_FUSION
        esi=sf.ice_saturation_pressure(s['t'])
        resistance=(ls/(sf.R_VAPOUR*s['t'])-1)*ls/(s['ka']*s['t'])+sf.R_VAPOUR*s['t']/(s['dv']*esi)
        rows.append([8*world.gravity*s['rho']/(math.pi*math.sqrt(geometry.area_ratio)*s['eta']**2),
                     s['eta']/s['rho'],max(0.,4*math.pi*geometry.capacitance_ratio*(s['e']/esi-1)/resistance),
                     (s['eta']/(s['rho']*s['dv']))**(1/3),math.exp(.09*(s['t']-sf.T0))])
    return p,np.ascontiguousarray(rows,float)


def _array(values,*,positive=False,ordered=False):
    result=np.ascontiguousarray(values,float)
    if result.ndim!=1 or not result.size or not np.isfinite(result).all():raise ValueError('Finite nonempty vector required')
    if positive and (result<=0).any():raise ValueError('Positive values required')
    if ordered and (np.diff(result)<0).any():raise ValueError('Ordered values required')
    return result


def properties(world,masses,temperature_c=-10.,geometry=GEOMETRY):
    masses=_array(masses,positive=True)
    p,a=parameters(world,[world.state(world.height_of(sf.T0+temperature_c))],geometry)
    out=np.zeros((masses.size,3));err=ct.create_string_buffer(1024)
    if library().snow_properties(p,a,masses,masses.size,out,err):raise RuntimeError(err.value.decode())
    return out


def box(world,masses,volume_m3,times_s,query_masses,*,temperature_c=-10.,seed=0,
        constant_kernel=None):
    masses=_array(masses,positive=True);times=_array(times_s,ordered=True)
    query=_array(query_masses,positive=True,ordered=True)
    if times[0]<0 or not math.isfinite(volume_m3) or volume_m3<=0:raise ValueError('Nonnegative times and positive volume required')
    if constant_kernel is not None and (not math.isfinite(constant_kernel) or constant_kernel<0):raise ValueError('Nonnegative finite kernel required')
    if not isinstance(seed,int) or not 0<=seed<2**64:raise ValueError('Unsigned 64-bit seed required')
    p,a=parameters(world,[world.state(world.height_of(sf.T0+temperature_c))])
    out=np.zeros((times.size,4));cdfs=np.zeros((times.size,2,query.size));info=np.zeros(3)
    err=ct.create_string_buffer(1024);start=time.monotonic()
    code=library().snow_box(p,a,masses,masses.size,volume_m3,times,times.size,query,query.size,
        seed,int(constant_kernel is not None),0 if constant_kernel is None else constant_kernel,
        GEOMETRY.area_ratio*math.pi/4,out,cdfs,info,err)
    if code:raise RuntimeError(err.value.decode())
    initial_mass=float(masses.sum())
    return dict(moments=out[:,:3]/volume_m3,maximum_mass=out[:,3],
                number_below=cdfs[:,0]/volume_m3,mass_cdf=cdfs[:,1]/initial_mass,
                relative_mass_residual=(out[:,1]-initial_mass)/initial_mass,
                proposals=int(info[0]),coagulations=int(info[1]),maximum_acceptance=float(info[2]),
                wall_seconds=time.monotonic()-start)


def diameter_queries():
    return np.unique(np.r_[np.geomspace(2e-5,100.,401),[.001,.003,.005,.01,.015,.03,.05,.15,.6,1.,10.]])


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
    code=library().snow_column(p,air,cfg,query,query.size,seed,out,err)
    if code:raise RuntimeError(err.value.decode())
    mass_names=('injected','deposited','freezing_exit','remaining','fractional_fragment_remaining','residual')
    ledger={k:float(v/area_m2) for k,v in zip(mass_names,out[:6])}
    ledger['relative_residual']=float(out[5]/(out[0]+out[1]))
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
    return dict(schema=SCHEMA,settings=dict(top_c=top_c,number_per_m3=number_per_m3,
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
