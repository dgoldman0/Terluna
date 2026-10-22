"""Replay one published single-crystal growth curve through native SCALE.

Takahashi et al. (1991), Table 2, -14.4 C: a and c are full lengths, mass
is in grams, and age is in minutes. Start at the measured three-minute fit,
not an invented nucleation state. Native terminal speed ventilates a stationary
particle in prescribed water-saturated air. This narrow empirical comparison
neither reproduces a snowfall nor validates lunar aggregation.
"""
from pathlib import Path
import argparse
import json
import math
import re
import numpy as np
from atmosphere.snow import scale_lunar as lunar
from shared.constants import STANDARD_GRAVITY, MOON_SURFACE_GRAVITY

CACHE=lunar.CACHE
OUT=lunar.ROOT/'atmosphere/snow/results/scale_crystal_growth.json'
TEMPERATURE=273.15-14.4
PRESSURE=101000.
START_MIN=3.
DURATION=27*60.
SOURCE_PDF='limitations_sources/takahashi1991.pdf'
SOURCE_SHA='2e5c6b492fcab7d1b3958be5a4c859fcd31174a1ac4865588ca48820bd4cd8eb'
CASES=[
    dict(name='growth_earth_observed_dt10',dt=10.,gravity=STANDARD_GRAVITY,shape='observed'),
    dict(name='growth_earth_observed_dt5',dt=5.,gravity=STANDARD_GRAVITY,shape='observed'),
    dict(name='growth_earth_observed_dt1',dt=1.,gravity=STANDARD_GRAVITY,shape='observed'),
    dict(name='growth_earth_sphere_dt5',dt=5.,gravity=STANDARD_GRAVITY,shape='sphere'),
    dict(name='growth_moon_observed_dt5',dt=5.,gravity=MOON_SURFACE_GRAVITY,shape='observed'),
]


def observed(age_min):
    """Original Table 2 coefficients converted to SI; fit valid through 30 min."""
    if not math.isfinite(age_min) or not 3<=age_min<=30:
        raise ValueError('This comparison uses the measured 3-30 minute interval')
    t=age_min
    mass=(6.4e-8*t**1.63 if t<4 else 2.3e-8*t**2.37)*1e-3
    diameter=(.14*t**.95 if t<4 else .13*t**1.02)*1e-3
    thickness=.0044*t**.52*1e-3
    return dict(age_min=t,mass_kg=mass,diameter_m=diameter,thickness_m=thickness)


def initial(shape):
    if shape not in ('observed','sphere'):raise ValueError('Unknown initial shape')
    obs=observed(START_MIN)
    a=obs['diameter_m']/2
    c=obs['thickness_m']/2 if shape=='observed' else a
    rho=obs['mass_kg']/(4*math.pi/3*a*a*c)
    return dict(**{**obs,'thickness_m':2*c},equatorial_radius_m=a,polar_radius_m=c,
                density_kg_m3=rho,observed_fit_at_start=obs)


def configure(name,dt,gravity,shape,cache=CACHE):
    cache=Path(cache)
    if lunar.base.digest(cache/SOURCE_PDF)!=SOURCE_SHA:
        raise ValueError('Unidentified laboratory source')
    if not math.isfinite(dt) or dt<=0 or 60/dt!=round(60/dt):
        raise ValueError('Use a positive timestep dividing one minute')
    start=initial(shape)
    d=lunar.configure(name,dt=5.,duration=DURATION,layers=4,capacity_per_layer=4,
                      collisions=False,cache=cache)
    for filename in ('init.conf','run.conf'):
        text=(d/filename).read_text();text=lunar.boxes.set_value(text,'CONST_GRAV',gravity)
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
    # Murphy-Koop saturation over liquid water, as in the existing column.
    vapor_pressure=float(lunar.flakes.liquid_saturation_pressure(TEMPERATURE))
    lines=(d/'lunar_forcing.txt').read_text().splitlines()
    values=list(map(float,lines[1].split()))
    values[2]=0. # no seed injection
    values[3]=1e-6 # kg m-3 s-1: generous ceiling, checked against actual demand
    values[5]=0.
    lines[1]=' '.join(format(x,'.17g') for x in values)
    lines[3:]=[' '.join(format(x,'.17g') for x in (TEMPERATURE,PRESSURE,vapor_pressure))]*s['layers']
    (d/'lunar_forcing.txt').write_text('\n'.join(lines)+'\n')
    s.update(dt_s=dt,gravity_m_s2=gravity,shape=shape,initial=start,temperature_K=TEMPERATURE,
             pressure_Pa=PRESSURE,vapor_pressure_Pa=vapor_pressure,source_pdf_sha256=SOURCE_SHA,
             laboratory_start_age_min=START_MIN,movement=False,initial_crystals=8,
             deposition_ceiling_kg_m3_s=values[3],gather_start_s=0.,
             profile=[dict(t=TEMPERATURE,p=PRESSURE,e=vapor_pressure)]*s['layers'],
             forcing='Stationary isolated crystals; terminal-speed ventilation, constant T,p,e, no injection, collisions, liquid, melting, or transport. Unlimited in practice deposition ceiling verified from ledger.',
             unmatched='Empirical fit at one Earth laboratory temperature, starting at three minutes. Initial mass and full axes map to a porous spheroid; observed thickness is inferred from mass and area. Spherical control preserves initial mass and maximum dimension. Moon run is an unvalidated gravity sensitivity at the same laboratory air state.',
             initial_water_kg_m2=start['mass_kg'],
             ledger_convention='Native hook assumes an empty initial population. Its raw residual therefore equals the known initial ice mass per area; independently add initial_water_kg_m2 to its source side.',
             growth_producer_sha256=lunar.base.digest(Path(__file__)))
    lunar.boxes.save(d/'settings.json',s)
    names=('init.conf','run.conf','env.txt','particles.pe000000','lunar_forcing.txt','settings.json')
    lunar.boxes.save(d/'inputs.json',{n:lunar.base.digest(d/n) for n in names})
    return d


def analyze(d):
    from netCDF4 import Dataset
    s=json.loads((d/'settings.json').read_text());rows=[]
    files=sorted(d.glob('SD_output_NetCDF_*.pe000000'))
    if not files:raise ValueError('No native particle snapshots')
    for path in files:
        with Dataset(path) as nc:
            v={k:np.asarray(x[:]) for k,x in nc.variables.items()}
        live=(v['sd_n']>0)&(v['sd_z']>=0)
        if live.sum()!=8 or np.any(v['sd_n'][live]!=1) or np.any(v['sdi_nmono'][live]!=1):
            raise ValueError('Population changed in an isolated growth test')
        a,c,rho=(v[k][live] for k in ('sdi_re','sdi_rp','sdi_rho'))
        mass=4*math.pi/3*a*a*c*rho
        # Native output is written before the microphysics update at its label.
        # The initial t=0 output is the initial population; later labels must be
        # checked against ledger mass to determine the actual completed growth.
        label=float(re.search(r'(\d{11}\.\d{3})',path.name).group(1))
        row=dict(file=path.name,sha256=lunar.base.digest(path),label_time_s=label,
                 mass_kg=float(mass.mean()),diameter_m=float((2*np.maximum(a,c)).mean()),
                 thickness_m=float((2*c).mean()),aspect=float((c/a).mean()),
                 density_kg_m3=float(rho.mean()),terminal_speed_m_s=float(v['sd_vz'][live].mean()),
                 maximum_relative_between_particle_difference=float(np.ptp(mass)/mass.mean()))
        rows.append(row)
    ledger=np.atleast_2d(np.loadtxt(d/'lunar_ledger.txt'))
    initial_mass=s['initial_water_kg_m2']
    matched=0
    for row in rows:
        elapsed=row['label_time_s']
        indices=np.flatnonzero(np.abs(ledger[:,0]-elapsed)<1e-8)
        error=None
        if elapsed==0:
            error=abs(row['mass_kg']-initial_mass)
        elif len(indices):
            error=abs(row['mass_kg']-ledger[indices[0],5]);matched+=1
        if error is not None and error>row['mass_kg']*1e-9:
            raise ValueError('Particle and ledger masses differ at matched times')
        row.update(elapsed_s=elapsed,age_min=START_MIN+elapsed/60.,ledger_mass_error_kg=error)
        obs=observed(row['age_min']);row['observed_fit']=obs
        row['relative_difference_from_fit']={k:row[k]/obs[k]-1 for k in ('mass_kg','diameter_m','thickness_m')}
    if matched<2 or rows[-1]['elapsed_s']!=DURATION:
        raise ValueError('Missing snapshot/ledger timing controls')
    with Dataset(d/'history.pe000000.nc') as nc:
        fields={k:np.asarray(nc[k][:]) for k in ('T','PRES','RH')}
    air={k:[float(v[1:].min()),float(v[1:].max())] for k,v in fields.items()}
    corrected=ledger[:,7]-initial_mass
    # One crystal per square metre in every horizontal cell. With identical
    # stationary crystals the per-area deposition equals each crystal's growth.
    increment=np.diff(np.r_[0.,ledger[:,2]])
    steps=np.diff(np.r_[0.,ledger[:,0]])
    budget=s['deposition_ceiling_kg_m3_s']*s['depth_m']/s['layers']*steps
    result=dict(settings=s,execution=json.loads((d/'execution.json').read_text()),
                snapshots=rows,air_ranges_after_first_step=air,
                water=dict(raw_residual_range_kg_m2=[float(ledger[:,7].min()),float(ledger[:,7].max())],
                           initial_water_kg_m2=initial_mass,
                           corrected_maximum_absolute_residual_kg_m2=float(np.abs(corrected).max()),
                           corrected_relative_residual=float(np.abs(corrected).max()/ledger[-1,5]),
                           ice_deposition_relative_residual=float(np.max(np.abs(ledger[:,5]-initial_mass-ledger[:,2]))/ledger[-1,5]),
                           maximum_used_fraction_of_growth_ceiling=float((increment/budget).max())),
                ledger_sha256=lunar.base.digest(d/'lunar_ledger.txt'),
                input_sha256=json.loads((d/'inputs.json').read_text()))
    lunar.boxes.save(d/'growth_analysis.json',result);return result


def assemble():
    results={s['name']:json.loads((CACHE/'runs'/s['name']/'growth_analysis.json').read_text()) for s in CASES}
    result=dict(schema='terluna.atmosphere.scale-crystal-growth/1',
                evidence='Native single-crystal comparison with one published wind-tunnel fit; empirical discrepancy, initialization and numerical refinement are distinct from snow-population validation.',
                source=dict(doi='10.2151/jmsj1965.69.1_15',table=2,page=24,temperature_C=-14.4,
                            original_pdf_sha256=SOURCE_SHA,age_min=[3,30],
                            units='Table mass g, full a/c lengths mm, age minutes; all outputs SI.',
                            limits='Published fits lack pointwise measurement errors. Planar thickness was inferred from mass and basal area. No aggregates or riming are represented by this selected curve. Corrigendum replaces Fig. 3, not Table 2.'),
                producer_sha256=lunar.base.digest(Path(__file__)),cases=results)
    lunar.boxes.save(OUT,result);return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['run','assemble'])
    args=parser.parse_args()
    if args.action=='run':
        for spec in CASES:
            d=CACHE/'runs'/spec['name']
            if not d.exists():d=configure(**spec)
            if not (d/'execution.json').exists():lunar.execute(d)
            analyze(d)
    assemble()
