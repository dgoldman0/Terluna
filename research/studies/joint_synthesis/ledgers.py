"""Ledgers across the five lines of work: the air by what it loses, gains and keeps; nitrogen; energy and heat.

    python -m research.studies.joint_synthesis.ledgers
    # -> results/ledgers.json

Each entry is read from a committed product or computed from one here, with its source. The atmosphere's mass is the
feasibility baseline's hydrostatic column (3.1e18 kg for 1.2 atm), so every loss rate carries its cycle time, the time
the loss and the resupply that matches it take to replace the whole atmosphere.
"""
from __future__ import annotations
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

from shared.constants import MOON_RADIUS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.joint-synthesis-ledgers/1'
P = lambda *parts: ROOT.joinpath(*parts)
OXYGEN = P('atmosphere', 'loss_response', 'results', 'oxygen_escape.json')
OCCURRENCE = P('atmosphere', 'electricity', 'results', 'occurrence_A28_dim5_moon.json')
NITROGEN = P('research', 'studies', 'atmospheric_electricity', 'results', 'nitrogen_oxides.json')
CANOPY = P('biosphere', 'canopy', 'results', 'canopy.json')
HEAT = P('research', 'studies', 'solar_shield_array', 'results', 'array_heat.json')
NIGHT = P('illumination', 'fleet_light', 'results', 'night_light.json')
OUT = HERE / 'results' / 'ledgers.json'
AIR_KG = 3.1e18                         # research/baselines/feasibility/reference.json, the 1.2 atm column
SECONDS_PER_YEAR = 365.25 * 86400.0
AIR_MOL = AIR_KG / 0.02897
H2_PPM = 0.53                           # the photochemistry's ground value (atmosphere/middle_atmosphere)
FLEET_RELEASE_PER_YEAR = (0.01, 0.1)    # share of the fleet's lift gas lost to the air each year, a bracket
NEW_NITROGEN_G_M2_YR = (1.0, 5.0)       # nitrogen a productive ecosystem must gain each year to replace its losses
EVIDENCE = ('Arithmetic on committed products of the five lines, with the brackets stated where a quantity has no '
            'product (the fleet\'s gas losses, an ecosystem\'s new nitrogen).')
READING_RULE = ('Rates per year unless labelled; cycle_years is the atmosphere\'s mass over the loss rate. Nitrogen in '
                'kg N per km2 a year. Power in W.')


def digest(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def cycle(kg_s):
    return AIR_KG / (kg_s * SECONDS_PER_YEAR)


def main(argv=None) -> int:
    area_m2 = 4 * np.pi * MOON_RADIUS ** 2
    # The air. The shield's cycle-mean losses for the warm titania case (requirements.md R1, final traced and
    # infrared-cooled screen), with and without the September magnets; the atoms and hydrogen from oxygen_escape.json.
    oxygen = json.loads(OXYGEN.read_text())
    h_escape = [x for v in oxygen['summary'].values() if isinstance(v, dict) and 'hydrogen_escape_kg_s' in v.get('titania_stack', {})
                for x in np.atleast_1d(v['titania_stack']['hydrogen_escape_kg_s'])]
    losses = dict(warm_case_no_magnets_kg_s=1.4, warm_case_september_magnets_kg_s=0.044, cooler_cases_kg_s=[0.61, 0.80])
    air = dict(mass_kg=AIR_KG,
               losses={k: dict(kg_s=v, cycle_years=([round(cycle(x), -8) for x in v] if isinstance(v, list) else round(cycle(v), -8)))
                       for k, v in losses.items()},
               hydrogen_escape_kg_s=[round(min(h_escape), 5), round(max(h_escape), 5)] if h_escape else None)
    h_air_kg = AIR_MOL * H2_PPM * 1e-6 * 0.002016
    escape_t_yr = [x * SECONDS_PER_YEAR / 1e3 for x in air['hydrogen_escape_kg_s']] if h_escape else [None]
    hydrogen = dict(air_inventory_t=round(h_air_kg / 1e3, -5), escape_t_per_year=[round(x, 0) for x in escape_t_yr],
                    fleet_inventory_t='31,000-54,000 in liners and about 17,000 in ferries (infrastructure review)',
                    release_bracket_per_year=FLEET_RELEASE_PER_YEAR,
                    release_t_per_year=[round(48000 * FLEET_RELEASE_PER_YEAR[0], -1), round(71000 * FLEET_RELEASE_PER_YEAR[1], -2)],
                    ppm_rise_per_mt_held=round(H2_PPM / (h_air_kg / 1e9), 5),
                    note=('Sinks: under the titania stack the design column has essentially no OH (atmosphere/middle_'
                          'atmosphere profiles), so released hydrogen leaves only through soils that take it up, as '
                          'Earth\'s do within about two years, or through escape, which is diffusion-limited and slow.'))
    # Nitrogen. Lightning carried over the Moon (atmosphere/electricity occurrence) against the new nitrogen a
    # productive ecosystem needs; the canopy model's equatorial stand for scale.
    occ, nit = json.loads(OCCURRENCE.read_text()), json.loads(NITROGEN.read_text())['estimate']
    ratio = occ['cases']['linear']['1']['ratio_to_box_moon_mean']
    lightning = {route: round(nit[route]['kg_n_per_km2_yr'] * ratio, 4) for route in ('per_flash', 'per_joule', 'per_metre')}
    need = [x * 1e3 for x in NEW_NITROGEN_G_M2_YR]                    # kg N per km2 a year
    canopy = json.loads(CANOPY.read_text())
    gpp = canopy['summaries']['moon']['0.0']['gpp_g_c_m2_per_24h']
    nitrogen = dict(lightning_moon_mean_kg_n_km2_yr=lightning, lightning_moon_total_kt_n_yr={k: round(v * area_m2 / 1e6 / 1e6, 2) for k, v in lightning.items()},
                    ecosystem_new_nitrogen_kg_n_km2_yr=need,
                    lightning_share_of_need=[round(min(lightning.values()) / need[1], 7), round(max(lightning.values()) / need[0], 5)],
                    canopy_gpp_g_c_m2_day_equator=gpp)
    # Energy and heat, from the Moon's absorbed sunlight to one city.
    heat = json.loads(HEAT.read_text())
    night = json.loads(NIGHT.read_text())
    absorbed = 293.3 * (1 - 0.28) * area_m2                            # design sunlight's global mean and the GCM's albedo
    energy = dict(moon_absorbed_sunlight_w=float(f'{absorbed:.3g}'), summit_metropolis_w=2.0e11,
                  summit_metropolis_local_w_m2=80.0, summit_metropolis_share_of_absorbed=float(f'{2e11 / absorbed:.2g}'),
                  magnets_refrigeration_w=[9.25e9, 1.48e10], film_plant_w=heat['film_plant']['power_TW'] if 'film_plant' in heat else None,
                  tiles_absorb_w='65-83 PW (array_heat.json)', fleet_night_half_intercepts_w=float(f"{night['night_half_luminous_flux_lm'] / 94.0:.3g}"),
                  heat_on_moon_per_tw_used_w_m2=heat['placement']['moon_W_m2_per_TW_used_on_moon'])
    product = dict(schema=SCHEMA, producer=dict(study='joint_synthesis', files={'research/studies/joint_synthesis/ledgers.py': digest(__file__)},
                                                inputs={str(p.relative_to(ROOT)): digest(p) for p in (OXYGEN, OCCURRENCE, NITROGEN, CANOPY, HEAT, NIGHT)}),
                   evidence=EVIDENCE, reading_rule=READING_RULE, air=air, hydrogen=hydrogen, nitrogen=nitrogen, energy=energy)
    OUT.write_text(json.dumps(product, indent=1, default=float) + '\n')
    print(json.dumps(dict(air=air, hydrogen={k: v for k, v in hydrogen.items() if k != 'note'}, nitrogen=nitrogen, energy=energy), indent=1, default=float)[:3500])
    return 0


if __name__ == '__main__':
    sys.exit(main())
