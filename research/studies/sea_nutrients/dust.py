"""The phosphorus that wind-blown dust brings the seas, from the Open Moon's own ground and winds.

    python -m research.studies.sea_nutrients.dust      # seconds

The phosphorus study carried an inherited range for dust, 0.7-28 mg of phosphorus per m2 a year, taken from Earth's
deposition and named there a sensitivity input. This module replaces it with a first lunar reckoning, from two things
the project has since worked out:

- the atlas round's sand study (research/studies/atlas_round/results/sand.json), which gives on the design climate
  how much sand a wind over bare loose ground carries in a year (the drift potential, m3 per metre of width, with
  the unresolved gusts as cases) and how much ground is bare and moving once plants and the kept crust hold the rest;
- the aerosol study's dust note (research drive, atmospheric-electricity/aerosol_dust_fog.md, section 8.2), which
  gives the dust a saltating surface raises as F = alpha Q, with Earth's alpha of 1e-5 to 1e-4 per metre for soils
  that hold fines and, as its central choice at lunar gravity, alpha scaled by gravity (0.165 of Earth's), so that
  the dust raised at a given friction velocity is about Earth's while the sand moved is 6.65 times Earth's.

The dust raised from a m2 of bare moving ground in a year is then alpha times the drift potential times the sand's
bulk density. It carries the soils' phosphorus and is taken to fall evenly over the Moon, since fine dust settles
six times slower here and travels far; near-source fallout on land would leave the seas less. Three things are not
in it and are named in the product: ground laid bare by people, given per percent of the arid and semi-arid land;
the same sum with Earth's alpha unscaled, the note's high scheme; and particles of living origin, which on Earth
carry about a tenth of the air's phosphorus and which the Open Moon's aerial life may make the larger share.
Writes results/dust.json.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from shared.constants import MOON_RADIUS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'dust.json'
SCHEMA = 'terluna.research.sea-nutrients-dust/1'
SAND = ROOT / 'research' / 'studies' / 'atlas_round' / 'results' / 'sand.json'
PHOSPHORUS = ROOT / 'resources' / 'results' / 'sea_phosphorus.json'
LAND_USE = ROOT / 'research' / 'studies' / 'atlas_round' / 'results' / 'land_use.json'
ALPHA_EARTH_PER_M = (1.0e-5, 1.0e-4)        # the dust note, 8.2: Earth's typical sandblasting efficiency (Kok et al. 2012;
#                                             González-Flórez et al. 2023 measured mostly 1e-5 to 1e-4 on a dry lake bed)
ALPHA_MOON_OVER_EARTH = 0.165               # the note's central choice: alpha follows gravity (Lu & Shao 1999; Shao 2004)
BULK_KG_M3 = 1600.0                         # the sand study's bulk density of loose sand
DRY_CLASSES = ('hyper-arid', 'arid', 'semi-arid')
GUSTS = dict(low='none', central='central', high='high')       # the sand study's cases for the wind inside a cell
MOVING = dict(low='young', central='central', high='opensand')  # its variants for the ground that is bare and moving


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def drift_potential(sand, case):
    """The mean drift potential over bare loose sand in the dry land (m3 per m of width a year), the classes weighed
    by their land."""
    groups = sand['winds']['over_sand'][case]['by_group']
    area = sand['land']['by_class_km2']
    w = np.array([area[c] for c in DRY_CLASSES], float)
    return float((np.array([groups[c]['drift_potential_m2_yr']['mean'] for c in DRY_CLASSES]) * w).sum() / w.sum())


def dust_raised(drift_m2_yr, alpha_per_m):
    """Dust raised from a m2 of bare moving ground in a year (kg): alpha times the sand carried past a metre of width."""
    return alpha_per_m * drift_m2_yr * BULK_KG_M3


def phosphorus_fall(dust_kg_yr, grade, moon_m2):
    """mg of phosphorus per m2 a year where dust of that phosphorus grade falls evenly over the Moon."""
    return 1e6 * dust_kg_yr * grade / moon_m2


def build():
    sand = json.loads(SAND.read_text())
    grade = json.loads(PHOSPHORUS.read_text())['land_loss']['land_soil_grade']
    moon = 4.0 * np.pi * MOON_RADIUS ** 2
    alpha = [ALPHA_MOON_OVER_EARTH * a for a in ALPHA_EARTH_PER_M]
    drift = {k: drift_potential(sand, case) for k, case in GUSTS.items()}
    moving = {k: sand['variants'][v]['moving_km2'] * 1e6 for k, v in MOVING.items()}

    def case(d, a, m, g):
        raised = dust_raised(drift[d], a)
        return dict(dust_kg_per_m2_of_moving_ground_yr=raised, dust_mt_yr=raised * moving[m] / 1e9,
                    phosphorus_mg_m2_yr=phosphorus_fall(raised * moving[m], g, moon))

    central = [case('central', alpha[0], 'central', grade[0]), case('central', alpha[1], 'central', grade[1])]
    outer = [case('low', alpha[0], 'low', grade[0]), case('high', alpha[1], 'high', grade[1])]
    # with no upkeep of the crust on sand, plants alone holding it: the sand study's open-sand variant
    untended = [case('central', alpha[0], 'high', grade[0]), case('central', alpha[1], 'high', grade[1])]
    dry_m2 = sum(sand['land']['by_class_km2'][c] for c in ('arid', 'semi-arid')) * 1e6
    bared = [phosphorus_fall(dust_raised(drift['central'], a) * 0.01 * dry_m2, g, moon) for a, g in zip(alpha, grade)]
    per_10000_km2 = [phosphorus_fall(dust_raised(drift['central'], a) * 1.0e10, g, moon) for a, g in zip(alpha, grade)]
    use = json.loads(LAND_USE.read_text())
    classes = use['inventory']['by_moisture_class']
    dry_all = sum(classes[c]['km2'] for c in ('hyper_arid', 'arid', 'semi_arid'))
    gardens = sum(use['variants']['A']['hybrid']['gardens_by_class'][c] for c in ('hyper_arid', 'arid', 'semi_arid'))
    inherited = json.loads(PHOSPHORUS.read_text())['seas']
    old = [1e6 * x / inherited['area_m2'] for x in inherited['dust_kg_yr_inherited_envelope']]        # kg to mg
    return dict(
        schema=SCHEMA,
        producer=dict(module='research.studies.sea_nutrients.dust', code_sha256=sha256(Path(__file__))),
        inputs={str(p.relative_to(ROOT)): sha256(p) for p in (SAND, PHOSPHORUS, LAND_USE)},
        evidence=('A first reckoning from the sand study\'s drift potentials and moving ground on the design climate '
                  'and the aerosol study\'s dust note: published relations fitted on Earth and carried to lunar '
                  'gravity, with no measurement of dust raised at 1.62 m/s2. The fallout is taken as even over the '
                  'Moon, and all of the dust\'s phosphorus as usable. Not validated.'),
        reading_rule=('natural: the phosphorus that dust from the ground the sand study finds bare and moving brings '
                      'each m2 of the Moon a year (mg), central as a range over the dust note\'s alpha, and outer '
                      'with the gusts, the alpha and the moving ground all at their low or all at their high ends. '
                      'per_percent_of_dry_land_bared: what each percent of the arid and semi-arid land laid bare and '
                      'loose would add in the central case. with_earths_alpha: the factor by which the note\'s high '
                      'scheme, alpha unscaled by gravity, multiplies every figure. natural.without_crust_upkeep: '
                      'the central case with the sand study\'s variant in which no crust is kept on sand and plants '
                      'alone hold it. per_10000_km2_bared: what each 10,000 km2 of dry ground left bare and loose '
                      'adds. dry_land: the hybrid gardens\' part of the dry land in the land-use screen. inherited: '
                      'the phosphorus study\'s earlier range, for comparison.'),
        setting=dict(alpha_per_m=alpha, alpha_earth_per_m=list(ALPHA_EARTH_PER_M), alpha_moon_over_earth=ALPHA_MOON_OVER_EARTH,
                     bulk_kg_m3=BULK_KG_M3, soil_phosphorus_grade=grade, moon_m2=moon,
                     drift_potential_m2_yr=drift, moving_ground_km2={k: v / 1e6 for k, v in moving.items()},
                     arid_and_semi_arid_km2=dry_m2 / 1e6),
        natural=dict(central=central, outer=outer,
                     without_crust_upkeep=untended,
                     phosphorus_mg_m2_yr=dict(central=[c['phosphorus_mg_m2_yr'] for c in central],
                                              outer=[c['phosphorus_mg_m2_yr'] for c in outer],
                                              without_crust_upkeep=[c['phosphorus_mg_m2_yr'] for c in untended])),
        per_percent_of_dry_land_bared_mg_m2_yr=bared,
        per_10000_km2_bared_mg_m2_yr=per_10000_km2,
        dry_land=dict(hyper_arid_arid_and_semi_arid_km2=dry_all, hybrid_gardens_km2=gardens,
                      hybrid_gardens_share=gardens / dry_all, wild_share_set=use['parameters']['register']['wild_shares'],
                      note='the land-use screen\'s variant A; the rest is wild land and land without use, with people\'s '
                           'fields and pasture a small part of it in the dry classes'),
        with_earths_alpha_factor=1.0 / ALPHA_MOON_OVER_EARTH,
        inherited_mg_m2_yr=old,
        not_counted=['particles of living origin: pollen, spores, fragments and the fall of aerial life',
                     'smoke', 'the share of dust\'s phosphorus that does not dissolve'])


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    out = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    n = out['natural']['phosphorus_mg_m2_yr']
    print(f'dust brings {n["central"][0]:.3f}-{n["central"][1]:.3f} mg P/m2/yr (outer {n["outer"][0]:.4f}-{n["outer"][1]:.2f}); '
          f'each percent of dry land bared adds {out["per_percent_of_dry_land_bared_mg_m2_yr"][0]:.2f}-'
          f'{out["per_percent_of_dry_land_bared_mg_m2_yr"][1]:.2f}; inherited {out["inherited_mg_m2_yr"]}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
