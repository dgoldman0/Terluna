"""Wave-resolving lunar shore profiles driven by incoming coastal spectra."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re

import numpy as np
from scipy.optimize import brentq
from scipy.signal import find_peaks

from climate.waves.coastal_checks import nesting_spectra
from climate.waves.model import ORIGIN, ROOT, sha256
from climate.waves.pilot import array_text, execute
from climate.waves.shore import SITES
from climate.waves.shore_checks import group_velocity
from climate.waves.swash_build import build
from geography.shore_profiles import numerical_profile
from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY

RUNS=ROOT/'research/runs/waves/runup'
PROFILES=ROOT/'research/runs/waves/shore_history/profiles.json'
SCENARIO=ROOT/'shared/scenarios/waves.json'
_scenario=json.loads(SCENARIO.read_text())
if _scenario['schema']!='terluna.scenario.waves/1':
    raise ValueError('Unsupported wave scenario')
WATER_DENSITY=_scenario['water_density_kg_m3']


def clock(seconds):
    milliseconds=round(seconds*1000)
    hour,remainder=divmod(milliseconds,3600000)
    minute,remainder=divmod(remainder,60000)
    second,ms=divmod(remainder,1000)
    return f'{hour:02d}{minute:02d}{second:02d}.{ms:03d}'


def spectral_integral(frequency,density,at):
    """Integral of a nonnegative piecewise-linear spectrum, zero outside it."""
    f,s=np.asarray(frequency),np.asarray(density)
    if f.ndim!=1 or s.shape!=f.shape or len(f)<2 or np.any(np.diff(f)<=0) or np.any(s<0):
        raise ValueError('Expected ordered frequencies and a nonnegative spectrum')
    cumulative=np.r_[0.,np.cumsum(np.diff(f)*(s[:-1]+s[1:])/2)]
    at=np.asarray(at)
    pos=np.clip(np.searchsorted(f,at)-1,0,len(f)-2)
    delta=np.clip(at,f[0],f[-1])-f[pos]
    return cumulative[pos]+delta*s[pos]+delta**2/2*(s[pos+1]-s[pos])/(f[pos+1]-f[pos])


def incoming_spectrum(path,station,*,date=None,layers=4,cycle_s=3600,gravity=MOON_SURFACE_GRAVITY,
                      frequency_ceiling=None,max_removed_variance=.01):
    """Collapse the incoming hemisphere while conserving shore-normal power."""
    path=Path(path).absolute()
    manifest=json.loads((path.parent/'product.json').read_text())
    if manifest.get('schema') not in ('terluna.climate.shore-case/1',
                                      'terluna.climate.coastal-history-case/1'):
        raise ValueError('Expected a completed coastal spectral product')
    if sha256(path)!=manifest['run']['output_sha256'][path.name]:
        raise ValueError('Coastal spectrum differs from its completed run record')
    source_input=path.parent/'INPUT'
    if sha256(source_input)!=manifest['run']['identity']['inputs']['INPUT']:
        raise ValueError('Coastal settings differ from their completed run record')
    setting=next(line for line in source_input.read_text().splitlines() if line.startswith('SET '))
    source_gravity=float(re.search(r'\bGRAV\s+([\d.eE+-]+)',setting)[1])
    source_density=float(re.search(r'\bRHO\s+([\d.eE+-]+)',setting)[1])
    if not (np.isclose(source_gravity,gravity,rtol=1e-10,atol=0) and
            np.isclose(source_density,WATER_DENSITY,rtol=1e-10,atol=0)):
        raise ValueError('Coastal and shoreline gravity or water density differ')
    source=nesting_spectra(path,{date}) if date is not None else nesting_spectra(path)
    if len(source['dates'])!=1:
        raise ValueError('Select exactly one coastal sea state')
    profiles=json.loads(PROFILES.read_text())
    if profiles.get('schema')!='terluna.geography.shore-profiles/1':
        raise ValueError('Expected a versioned shore-profile product')
    profile=profiles['profiles'][station]
    i=list(SITES).index(station)
    np.testing.assert_allclose(source['locations'][i],SITES[station],atol=1e-5,rtol=0)
    if not source['available'][0,i]:
        raise ValueError('The selected station has no spectrum')
    direction=np.deg2rad(source['direction'])
    unit=np.column_stack((np.cos(direction),np.sin(direction)))
    projection=np.maximum(unit@profile['normal_east_north'],0.)
    f=source['frequency']
    density=source['spectra'][0,i]@projection*(360/len(direction))
    total=float(np.trapezoid(density,f))
    if total<=0:
        raise ValueError('The selected shore has zero incoming variance')
    depth=profile['station_depth_m']
    cutoff=.9*layers/np.pi*np.sqrt(gravity/depth)
    if not 0<=max_removed_variance<1 or (frequency_ceiling is not None and
            (not np.isfinite(frequency_ceiling) or frequency_ceiling<=0)):
        raise ValueError('Expected a finite positive frequency ceiling and a variance bound below one')
    if frequency_ceiling is not None:
        cutoff=min(cutoff,frequency_ceiling)
    df=1/cycle_s
    frequency=np.arange(1,int(min(cutoff,f[-1])/df)+1)*df
    if len(frequency)==0:
        raise ValueError('The frequency ceiling leaves zero resolved components')
    edges=np.r_[frequency-df/2,frequency[-1]+df/2]
    variance=np.diff(spectral_integral(f,density,edges))
    target=variance/df
    retained=float(variance.sum())
    removed=1-retained/total
    if removed>max_removed_variance:
        raise ValueError(f'Vertical layers would discard {removed:.2%} of incoming variance')
    reference_power=float(np.trapezoid(density*group_velocity(f,depth,gravity),f))
    represented_power=float(np.sum(variance*group_velocity(frequency,depth,gravity)))
    text='SPEC1D\n'+array_text(np.column_stack((frequency,target)))
    record=dict(source_path=str(path.relative_to(ROOT)),source_sha256=sha256(path),
        source_record_sha256=sha256(path.parent/'product.json'),
        source_gravity_m_s2=source_gravity,source_water_density_kg_m3=source_density,
        source_date=None if source['dates'][0] is None else source['dates'][0].isoformat(),
        station=station,normal_east_north=profile['normal_east_north'],depth_m=depth,
        incoming_hs_m=4*np.sqrt(total),represented_hs_m=4*np.sqrt(retained),
        mean_period_s=float(total/np.trapezoid(f*density,f)),
        source_frequency_band_hz=[float(f[0]),float(f[-1])],
        source_edge_variance_fraction=dict(
            low=float(np.trapezoid(density[:3],f[:3])/total),
            high=float(np.trapezoid(density[-3:],f[-3:])/total)),
        removed_variance_fraction=removed,cutoff_frequency_hz=cutoff,
        explicit_frequency_ceiling_hz=frequency_ceiling,
        maximum_removed_variance_fraction=max_removed_variance,
        boundary_frequency_step_hz=df,cycle_s=cycle_s,components=len(frequency),
        source_variance_flux_m3_s=reference_power,represented_variance_flux_m3_s=represented_power,
        transport_relative_difference=abs(reference_power-represented_power)/reference_power,
        assumption='The 1D spectrum integrates E(f,theta)*max(n dot u,0) over direction. It preserves incoming normal power at the source depth. Alongshore variation and oblique-wave evolution require 2D wave-resolving calculations.')
    return text,record


def profile_files(x,depth,boundary,*,dx=.5,step=.025,layers=4,seconds=2400,
                  spinup=600,seed=12345678,gravity=MOON_SURFACE_GRAVITY,friction=0.,threshold=.005,
                  alpha=.6,beta=.3,regular=None,flat=False,water_density=WATER_DENSITY):
    x,depth=np.asarray(x),np.asarray(depth)
    if len(x)!=len(depth) or not np.allclose(np.diff(x),dx) or depth[0]<=0:
        raise ValueError('Expected a uniform wet-to-dry numerical profile')
    if not 0<=spinup<seconds or layers<1 or step<=0 or threshold<=0 or friction<0:
        raise ValueError('Expected a positive duration, time step, wet depth and layer count')
    n=len(x)-1
    boundary_command=(f"CON REG {regular[0]:.12g} {regular[1]:.12g}" if regular else "CON SPECFILE 'spectrum.dat'")
    text=f"""PROJECT 'Lunar runup' '001'
MODE DYN ONED
SET DEPMIN {threshold:.12g} SEED {seed} GRAV {gravity:.12g} RHOWAT {water_density:.12g}
OUTPUT OPTIONS TABLE 16
CGRID 0 0 0 {x[-1]:.12g} 0 {n} 0
VERT {layers}
INPGRID BOTTOM 0 0 0 {n} 0 {dx:.12g} 1
READINP BOTTOM 1 'depth.bot' 1 0 FREE
INIT ZERO
BOU SIDE W CCW BTYPE WEAK SMOO 120 SEC {boundary_command}
BOU SIDE E CCW BTYPE RADIATION
FRIC CONSTANT {friction:.12g}
BREAK {alpha:.12g} {beta:.12g}
NONHYDROSTATIC BOX PRECONDI ILU
DISCRET UPW MOM
TIMEI METH EXPL .1 .5
GROUP 'bed' 1 {n+1} 1 1
QUANT HS SETUP DUR {seconds-spinup:.12g} SEC
QUANT RUNUP DELRP {2*threshold:.12g}
TABLE 'bed' NOHEAD 'statistics.tbl' XP BOTL HS SETUP
TABLE 'bed' NOHEAD 'field.tbl' TSEC XP BOTL WATL DEP BRKP OUTPUT {clock(spinup)} 1 SEC
POINTS 'gauges' 25 0 50 0 75 0
TABLE 'gauges' NOHEAD 'gauges.tbl' TSEC XP WATL OUTPUT {clock(0)} .1 SEC
"""
    outputs=['statistics.tbl','field.tbl','gauges.tbl']
    if flat:
        text+=f'SPON E {min(x[-1]/4,40):.12g}\n'
    else:
        text+=f"TABLE 'NOGRID' NOHEAD 'runup.tbl' TSEC RUNUP OUTPUT {clock(spinup)} .1 SEC\n"
        outputs.append('runup.tbl')
    text+=f'COMPUTE {clock(0)} {step:.12g} SEC {clock(seconds)}\nSTOP\n'
    return dict(INPUT=text,**{'depth.bot':array_text(depth[None,:]),'spectrum.dat':boundary}),outputs


def controls(executable):
    """Check linear-wave dispersion and gravity similarity in the executable."""
    result={}
    x=np.arange(0,160.01,.5)
    for name,g in [('moon',MOON_SURFACE_GRAVITY),('four_gravity',4*MOON_SURFACE_GRAVITY),('earth',STANDARD_GRAVITY)]:
        scale=np.sqrt(MOON_SURFACE_GRAVITY/g)
        period=10*scale
        files,outputs=profile_files(x,np.full_like(x,5.),'',dx=.5,layers=3,
            step=.025*scale,seconds=300*scale,spinup=200*scale,gravity=g,regular=(.02,period),flat=True)
        # Scale the ramp and all output times with the physical clock.
        files['INPUT']=files['INPUT'].replace('SMOO 120 SEC',f'SMOO {20*scale:.12g} SEC')
        files['INPUT']=files['INPUT'].replace("OUTPUT 000000.000 .1 SEC",f"OUTPUT 000000.000 {.1*scale:.12g} SEC")
        directory=RUNS/f'control_{name}_ilu'
        run=execute(executable,directory,files,outputs,timeout_s=300,threads=1)
        data=np.loadtxt(directory/'gauges.tbl').reshape(-1,3,3)
        selected=data[:,0,0]>=200*scale
        t=data[selected,0,0]
        basis=np.column_stack((np.cos(2*np.pi*t/period),np.sin(2*np.pi*t/period),np.ones(len(t))))
        coefficient=np.linalg.lstsq(basis,data[selected,:,2],rcond=None)[0]
        amplitude=np.hypot(coefficient[0],coefficient[1])
        phase=np.unwrap(np.arctan2(coefficient[1],coefficient[0]))
        wave_number=brentq(lambda k:g*k*np.tanh(k*5)-(2*np.pi/period)**2,1e-8,10)
        # Choose the known dispersion branch between the gauges.
        expected=wave_number*25
        measured=np.mod(np.diff(phase),2*np.pi)
        measured+=2*np.pi*np.rint((expected-measured)/(2*np.pi))
        phase_error=np.abs(measured-expected)/expected
        result[name]=dict(gravity_m_s2=g,period_s=period,amplitudes_m=amplitude.tolist(),
            amplitude_relative_max=float(np.max(abs(amplitude/.01-1))),
            phase_relative_max=float(phase_error.max()),run=run,
            passed=bool(np.max(abs(amplitude/.01-1))<.05 and phase_error.max()<.02))
    moon=np.loadtxt(RUNS/'control_moon_ilu/gauges.tbl').reshape(-1,3,3)
    scaled=np.loadtxt(RUNS/'control_four_gravity_ilu/gauges.tbl').reshape(-1,3,3)
    np.testing.assert_allclose(moon[:,0,0],2*scaled[:,0,0],rtol=0,atol=1e-4)
    similarity=float(np.max(abs(moon[:,:,2]-scaled[:,:,2]))/.01)
    path=RUNS/'controls.json'
    path.write_text(json.dumps(dict(schema='terluna.climate.swash-controls/1',cases=result,
        gravity_similarity=dict(gravity_ratio=4,clock_ratio=.5,amplitude_normalized_max_error=similarity,
                                passed=similarity<.01),producer_sha256=sha256(Path(__file__))),indent=2)+'\n')
    print(json.dumps(dict(cases={k:{q:v[q] for q in ('amplitude_relative_max','phase_relative_max','passed')} for k,v in result.items()},
                          similarity=similarity)),flush=True)
    return result


def summarize(directory,settings,spectrum):
    r=np.loadtxt(directory/'runup.tbl')
    if not np.isfinite(r).all() or np.any(r[:,1]<-100):
        raise ValueError('Run-up output contains unavailable values')
    dt=float(np.median(np.diff(r[:,0])))
    peaks,_=find_peaks(r[:,1],distance=max(1,round(.35*spectrum['mean_period_s']/dt)),prominence=.01)
    profile=np.loadtxt(directory/'statistics.tbl')
    count=len(profile)
    field=np.loadtxt(directory/'field.tbl').reshape(-1,count,6)
    if not np.allclose(field[0,:,1:3],profile[:,:2],rtol=1e-6,atol=1e-5):
        raise ValueError('Profile output coordinates differ')
    flags=field[:,:,5]
    if np.any((flags!=0)&(flags!=1)):
        raise ValueError('Breaking flags must be sampled on native nodes')
    frequency=flags.mean(axis=0)
    active=frequency>=.01
    half=len(r)//2
    result=dict(runup_max_m=float(r[:,1].max()),runup_mean_m=float(r[:,1].mean()),
        runup_time_p98_m=float(np.percentile(r[:,1],98)),
        runup_peak_p98_m=float(np.percentile(r[peaks,1],98)) if len(peaks) else None,
        runup_peak_count=len(peaks),
        runup_p98_first_half_m=float(np.percentile(r[:half,1],98)),
        runup_p98_second_half_m=float(np.percentile(r[half:,1],98)),
        maximum_profile_hs_m=float(profile[:,2].max()),
        breaking_nodes=int(active.sum()),
        breaking_depth_range_m=[float(profile[active,1].min()),float(profile[active,1].max())] if active.any() else None,
        breaking_x_range_m=[float(profile[active,0].min()),float(profile[active,0].max())] if active.any() else None,
        duration_s=float(r[-1,0]-r[0,0]),
        runup_definition='Vertical waterline elevation above still water at the stated minimum wet depth. Time p98 is exceeded for 2% of sampled time. Peak p98 uses local maxima separated by at least 0.35 mean periods with 1 cm prominence.',
        breaking_definition='Native nodes flagged as breaking in at least 1% of the one-second samples during the reporting window.')
    return result


def run_case(executable,source,*,station='west_face',date=None,kind='rock',slope=50,dx=.5,
             step=.025,layers=4,seconds=2400,spinup=600,seed=12345678,friction=0.,name=None,
             threshold=.005,alpha=.6,buffer_m=100.,frequency_ceiling=None,max_removed_variance=.01):
    producer=Path(__file__).read_bytes()
    profiles=json.loads(PROFILES.read_text())
    x,depth,geometry=numerical_profile(profiles,station,dx=dx,kind=kind,slope_denominator=slope,buffer_m=buffer_m)
    text,spectrum=incoming_spectrum(source,station,date=date,layers=layers,
                                   frequency_ceiling=frequency_ceiling,max_removed_variance=max_removed_variance)
    files,outputs=profile_files(x,depth,text,dx=dx,step=step,layers=layers,seconds=seconds,
        spinup=spinup,seed=seed,friction=friction,threshold=threshold,alpha=alpha)
    name=name or f'{station}_{kind}{slope if kind=="slope" else ""}_dx{str(dx).replace(".","p")}_k{layers}'
    directory=RUNS/name
    run=execute(executable,directory,files,outputs,timeout_s=7200,threads=1)
    product_path=directory/'product.json'
    if product_path.exists():
        previous=json.loads(product_path.read_text())
        keys=['source_path','source_sha256','source_date','station']
        if 'source_record_sha256' in previous['boundary']:
            keys.append('source_record_sha256')
        if (previous['schema']!='terluna.climate.runup-case/1' or previous['run']!=run or
                previous['profile_product_sha256']!=sha256(PROFILES) or
                any(previous['boundary'][key]!=spectrum[key] for key in keys) or
                sha256(directory/'producer_source.py')!=previous['producer_sha256']):
            raise ValueError('Completed shore metadata differs from its inputs or source record')
        if 'water_scenario' in previous and previous['water_scenario']['sha256']!=sha256(SCENARIO):
            raise ValueError('Completed shore case used a different wave scenario')
        print(json.dumps(dict(name=name,cached=True)),flush=True)
        return previous
    settings=dict(geometry=geometry,step_s=step,layers=layers,seconds=seconds,spinup_s=spinup,
                  seed=seed,friction_coefficient=friction,wet_threshold_m=threshold,
                  gravity_m_s2=MOON_SURFACE_GRAVITY,water_density_kg_m3=WATER_DENSITY,
                  breaking_alpha=alpha,breaking_beta=.3)
    result=dict(schema='terluna.climate.runup-case/1',name=name,
        evidence='Non-hydrostatic wave-resolving simulation of an explicit lunar shore profile and one imposed coastal spectrum.',
        reading_rule='The directional sea is reduced to a power-equivalent normal-incidence spectrum. Rock profiles retain 118 m geographic information; smooth slopes are explicit beach scenarios. Wave phases, friction and breaking parameters remain model assumptions.',
        producer_sha256=hashlib.sha256(producer).hexdigest(),profile_product_sha256=sha256(PROFILES),
        water_scenario=dict(path=str(SCENARIO.relative_to(ROOT)),sha256=sha256(SCENARIO)),
        settings=settings,boundary=spectrum,run=run,statistics=summarize(directory,settings,spectrum))
    (directory/'producer_source.py').write_bytes(producer)
    product_path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(name=name,elapsed_wall_s=run['elapsed_wall_s'],statistics=result['statistics'])),flush=True)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('controls','case'))
    parser.add_argument('--source',type=Path)
    parser.add_argument('--station',default='west_face',choices=tuple(SITES)[:3])
    parser.add_argument('--hour',type=int)
    parser.add_argument('--kind',choices=('rock','slope'),default='rock')
    parser.add_argument('--slope',type=int,default=50)
    parser.add_argument('--dx',type=float,default=.5)
    parser.add_argument('--step',type=float,default=.025)
    parser.add_argument('--layers',type=int,default=4)
    parser.add_argument('--seconds',type=int,default=2400)
    parser.add_argument('--spinup',type=int,default=600)
    parser.add_argument('--seed',type=int,default=12345678)
    parser.add_argument('--friction',type=float,default=0.)
    parser.add_argument('--threshold',type=float,default=.005)
    parser.add_argument('--alpha',type=float,default=.6)
    parser.add_argument('--buffer',type=float,default=100.)
    parser.add_argument('--frequency-ceiling',type=float)
    parser.add_argument('--max-removed-variance',type=float,default=.01,
                        help='Explicit variance-loss bound for a labelled frequency-band sensitivity')
    parser.add_argument('--name')
    args=parser.parse_args()
    executable=build()
    if args.action=='controls':
        controls(executable)
    else:
        from datetime import timedelta
        if args.source is None:
            parser.error('--source is required for a shore case')
        run_case(executable,args.source,station=args.station,
            date=None if args.hour is None else ORIGIN+timedelta(hours=args.hour),kind=args.kind,
            slope=args.slope,dx=args.dx,step=args.step,layers=args.layers,seconds=args.seconds,
            spinup=args.spinup,seed=args.seed,friction=args.friction,name=args.name,
            threshold=args.threshold,alpha=args.alpha,buffer_m=args.buffer,
            frequency_ceiling=args.frequency_ceiling,max_removed_variance=args.max_removed_variance)


if __name__=='__main__':
    main()
