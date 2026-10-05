"""The Open Moon's aerosol near the ground by region, and the air's conductivity and cloud nuclei it gives.

    python -m research.studies.open_moon_aerosol.run

writes results/open_moon_aerosol.json beside this file (schema terluna.research.open-moon-aerosol/1). It needs the GCM
climatology product (python -m climate.gcm.climatology A28_dim5_moon:20-29) and the conductivity column
(python -m atmosphere.electricity.conductivity).
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

from atmosphere.electricity import conductivity as cd
from research.studies.open_moon_aerosol import budget as bu
from research.studies.open_moon_aerosol import inputs as inp
from research.studies.open_moon_aerosol import particles as pa

HERE = Path(__file__).resolve().parent
RESULTS = HERE / 'results'
SCHEMA = 'terluna.research.open-moon-aerosol/1'
COLUMN = inp.ROOT / 'atmosphere' / 'electricity' / 'results' / 'conductivity_moon.json'
SUPERSATURATIONS = (0.001, 0.003, 0.01)
CURRENT_COLUMN = [pa.Mode(1.0e8, 1.0e-7, 1.01, cd.KAPPA, 'the column as it stands')]
EVIDENCE = ('A structured estimate. Each region\'s particles near the ground are a measured Earth analogue (the clean '
            'wet-season Amazon, the remote ocean by day, the Arctic summer) changed by lunar factors from the literature '
            'notes on the data drive (aerosol_marine.md, aerosol_biogenic.md, aerosol_dust_fog.md, aerosol_budget.md) '
            'and the project\'s own models: the regions\' areas from the GCM and geography, their weather near the ground '
            'from CM1\'s rings, the ionization from the conductivity column. Clean, central and loaded cases run from the '
            'fewest particles to the most. Nothing here is a measurement of the Open Moon; the design levers and the '
            'gaps are listed with the results.')
READING_RULE = ('Number concentrations per cm3; diameters in nm (dry medians of lognormal modes); conductivity in S/m '
                'at the ground (2 m), total of both signs; relaxation time eps0/sigma in minutes. Day and night are '
                'hour angles within and beyond 95 degrees of noon. The night values are means over the dark hours; '
                '"with fog" weights the clear and foggy conductivity by the share of hours in fog. Area shares are of '
                'the whole Moon. CCN at 0.1, 0.3 and 1 % supersaturation.')


def ground_air() -> dict:
    col = json.loads(COLUMN.read_text())['solar_minimum']
    return dict(q_m3_s=col['ion_pairs_per_cm3_s'][0] * 1e6, p_pa=col['pressure_hpa'][0] * 100.0)


def fog_conductivity(modes, q, p, t, case) -> float:
    """Conductivity in fog: the small ions go to the droplets and to the particles between them (at 99 % humidity)."""
    mu = cd.mobility_air(cd.MOBILITY_STANDARD[0], p, t)
    droplets = cd.droplet_attachment(t, mu, bu.FOG_DROPLETS[case] * bu.CM3, bu.FOG_RADIUS_M)
    sink = droplets + pa.ion_sink(modes, p, t, 0.99)
    return float(cd.conductivity_of(q, p, t, sink)[0])


def describe(modes) -> list:
    return [dict(name=m.name, number_cm3=m.number_m3 / bu.CM3, d_nm=m.d_m * 1e9, sigma=m.sigma, kappa=m.kappa)
            for m in modes if m.number_m3 > 0.0]


def state(region: str, case: str, phase: str, weather: dict, air: dict) -> dict:
    w = weather[phase]
    t = w['air_c'] + 273.15
    depth = max(weather['night']['mixed_layer_km'] * 1000.0, 50.0)
    modes = bu.modes(region, case, phase, depth)
    q, p = air['q_m3_s'], air['p_pa']
    clear = pa.conductivity(modes, q, p, t, min(w['relative_humidity'], 0.99))
    foggy = fog_conductivity(modes, q, p, t, case)
    mean = (1.0 - w['fog']) * clear + w['fog'] * foggy
    big = sum(m.number_m3 * 0.5 * float(math.erfc(np.log(100e-9 / m.d_m) / (np.sqrt(2.0) * np.log(m.sigma))))
              for m in modes) / bu.CM3
    return dict(modes=describe(modes), number_cm3=sum(m.number_m3 for m in modes) / bu.CM3, above_100nm_cm3=big,
                ccn_cm3={f'{s * 100:g}%': pa.ccn(modes, s, t) / bu.CM3 for s in SUPERSATURATIONS},
                ion_sink_s=pa.ion_sink(modes, p, t, min(w['relative_humidity'], 0.99)),
                conductivity_s_m=clear, in_fog_s_m=foggy, fog_share=w['fog'], with_fog_s_m=mean,
                relaxation_min=cd.EPSILON_0 / mean / 60.0, relative_humidity=w['relative_humidity'], air_c=w['air_c'])


def earth_checks() -> dict:
    """The same machinery on Earth's measured analogues, against measured ions and conductivity."""
    amazon = [pa.Mode(246e6, 70e-9, 1.6, 0.13), pa.Mode(145e6, 170e-9, 1.5, 0.21), pa.Mode(0.3e6, 1.6e-6, 1.6, 0.02)]
    out = dict(amazon_small_ions=dict(measured_cm3=[549.0, 856.0], source='Wimmer et al. 2018, inside the rainforest at '
                                      '97 % relative humidity (positive, negative; biogenic note Sect. 3)', model={}))
    for q in (2.5, 5.0, 10.0):
        n = cd.conductivity_of(q * 1e6, 100000.0, 298.0, pa.ion_sink(amazon, 100000.0, 298.0, 0.97))[1]
        out['amazon_small_ions']['model'][f'q {q:g} per cm3 per s'] = float(n) / bu.CM3
    marine = [pa.Mode(250e6, 45e-9, 1.5, 0.4), pa.Mode(200e6, 165e-9, 1.5, 0.5), pa.Mode(15e6, 0.4e-6, 2.0, 1.1)]
    out['marine_conductivity'] = dict(measured_s_m=[1.0e-14, 2.3e-14],
                                      source='Kamra et al. 1997 (1.1-2.3e-14 at 70-90 % RH) and Siingh et al. 2005 '
                                             '(1.0-1.7e-14 over the Arabian Sea), via humid_conductivity.md Sect. 4',
                                      model={})
    for q in (1.5, 2.0, 3.0):
        s = pa.conductivity(marine, q * 1e6, 101325.0, 300.0, 0.8)
        out['marine_conductivity']['model'][f'q {q:g} per cm3 per s'] = s
    return out


# design levers, each set in the central case alone (biogenic Sect. 7.6; marine Sect. 7; dust Sect. 8.4): the share of
# plants that emit isoprene (none, or every tree; isoprene is about two-fifths of the central forest's secondary organic
# aerosol, biogenic Sect. 7.2), the seas' dimethyl sulfide (1 or 10 nM against 3: marine day numbers 0.6 or 1.6 times),
# the plains' crust cover (complete, or bare grazed ground: dust 0 or 5 times), and the night decomposers' fungal bursts
# (none, or twice as often)
LEVERS = {
    'no isoprene emitters': dict(SOA_RATIO=dict(central=0.5 * 0.6)),
    'every tree an isoprene emitter': dict(SOA_RATIO=dict(central=0.5 * (0.6 + 0.4 * 2.5))),
    'seas at 1 nM dimethyl sulfide': dict(MARINE_SCALE=0.6),
    'seas at 10 nM dimethyl sulfide': dict(MARINE_SCALE=1.6),
    'plains fully crusted': dict(DUST_SCALE=0.0),
    'plains bare and grazed': dict(DUST_SCALE=5.0),
    'no night fungal bursts': dict(FUNGAL_SCALE=0.0),
    'night fungal bursts twice as often': dict(FUNGAL_SCALE=2.0),
}


def lever_effect(name: str, air: dict, areas: dict) -> dict:
    """Area-mean conductivity at the ground, day and night (with fog), with one lever set in the central case."""
    saved = dict(SOA_RATIO=dict(bu.SOA_RATIO), MARINE_DAY=dict(bu.MARINE_DAY), DUST_DUTY={k: dict(v) for k, v in
                 bu.DUST_DUTY.items()}, FUNGAL_SHARE={k: dict(v) for k, v in bu.FUNGAL_SHARE.items()})
    lever = LEVERS[name]
    try:
        if 'SOA_RATIO' in lever:
            bu.SOA_RATIO['central'] = lever['SOA_RATIO']['central']
        if 'MARINE_SCALE' in lever:
            n_a, d_a, n_c, d_c, spray = bu.MARINE_DAY['central']
            f = lever['MARINE_SCALE']
            bu.MARINE_DAY['central'] = (n_a * f, d_a, n_c * f, d_c, spray)
        if 'DUST_SCALE' in lever:
            for v in bu.DUST_DUTY.values():
                v['central'] *= lever['DUST_SCALE']
        if 'FUNGAL_SCALE' in lever:
            for v in bu.FUNGAL_SHARE.values():
                v['central'] *= lever['FUNGAL_SCALE']
        out = {}
        for phase in ('day', 'night'):
            out[phase] = sum(share * state(region, 'central', phase, inp.regime(region), air)['with_fog_s_m']
                             for region, share in areas.items())
        return out
    finally:
        bu.SOA_RATIO.clear(); bu.SOA_RATIO.update(saved['SOA_RATIO'])
        bu.MARINE_DAY.clear(); bu.MARINE_DAY.update(saved['MARINE_DAY'])
        bu.DUST_DUTY.clear(); bu.DUST_DUTY.update(saved['DUST_DUTY'])
        bu.FUNGAL_SHARE.clear(); bu.FUNGAL_SHARE.update(saved['FUNGAL_SHARE'])


def run() -> dict:
    air = ground_air()
    areas = inp.areas()
    regions = {}
    for region, share in areas.items():
        weather = inp.regime(region)
        regions[region] = dict(area_share=share, weather=weather,
                               cases={case: {phase: state(region, case, phase, weather, air) for phase in ('day', 'night')}
                                      for case in bu.CASES})
    summary = {}
    for case in bu.CASES:
        for phase in ('day', 'night'):
            vals = [(regions[r]['area_share'], regions[r]['cases'][case][phase]) for r in regions]
            summary[f'{case}_{phase}'] = dict(
                conductivity_s_m=sum(a * v['conductivity_s_m'] for a, v in vals),
                with_fog_s_m=sum(a * v['with_fog_s_m'] for a, v in vals),
                number_cm3=sum(a * v['number_cm3'] for a, v in vals),
                ccn_0_3_cm3=sum(a * v['ccn_cm3']['0.3%'] for a, v in vals))
    t_ground = 298.1
    current = dict(conductivity_s_m=pa.conductivity(CURRENT_COLUMN, air['q_m3_s'], air['p_pa'], t_ground, 0.84),
                   ccn_cm3={f'{s * 100:g}%': pa.ccn(CURRENT_COLUMN, s, t_ground) / bu.CM3 for s in SUPERSATURATIONS})
    land = areas['wet_land'] + areas['fog_desert'] + areas['polar_dry_land']
    so2 = {f'{dms:g}': bu.sulfur_dioxide_ppt(dms, areas['seas'], land, 5000.0) for dms in (0.4, 1.2, 3.9)}
    return dict(ground=dict(ion_pairs_per_cm3_s=air['q_m3_s'] / 1e6, pressure_hpa=air['p_pa'] / 100.0),
                regions=regions, area_mean=summary, current_column=current,
                levers={name: lever_effect(name, air, areas) for name in LEVERS},
                plains_sulfur_dioxide_ppt_by_sea_dms_umol_m2_day=so2, earth_checks=earth_checks())


def main(argv=None) -> int:
    result = run()
    digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]
    files = {f'{f.name}': digest(f) for f in sorted(HERE.glob('*.py'))}
    out = dict(schema=SCHEMA, evidence=EVIDENCE, reading_rule=READING_RULE,
               producer=dict(domain='research', files=files),
               inputs={str(p.relative_to(inp.ROOT)): digest(p) for p in (inp.CLIMATOLOGY, inp.DRAINAGE, inp.ATLAS, COLUMN)},
               **result)
    RESULTS.mkdir(parents=True, exist_ok=True)
    path = RESULTS / 'open_moon_aerosol.json'
    path.write_text(json.dumps(out, indent=1, default=float, allow_nan=False) + '\n')
    print(path)
    return 0


if __name__ == '__main__':
    sys.exit(main())
