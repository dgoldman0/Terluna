#!/usr/bin/env python3
"""Carbon on the Open Moon: what plants need, how fast the land would take it, and what a settled Moon returns.

    python -m research.studies.atmospheric_co2.run            # writes research/studies/atmospheric_co2/results/atmospheric_co2.json

A literature synthesis with small calculations, not a carbon-cycle model (see README.md). It computes:

1. the design air's CO2: partial pressure, mole fraction, column mass and global inventory;
2. light-saturated C3 photosynthesis against CO2 partial pressure, from the biochemical model of
   Farquhar, von Caemmerer and Berry with the kinetics of Bernacchi et al. (2001) at a fixed ratio of
   internal to ambient CO2, and what the 1.2-atm total pressure costs at a fixed stomatal conductance;
3. drawdown: an Earth-like land biosphere's carbon against the atmosphere's inventory, and silicate
   weathering from the basalt law of Dessert et al. (2003), gross (to bicarbonate) and net (after
   carbonate precipitation returns half of it);
4. the heat that returning the net weathering sink would take by calcining carbonate.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
if __package__ in (None, ''):                  # executed as a file path rather than with -m
    sys.path.insert(0, str(ROOT))
from shared.constants import MOON_RADIUS, MOON_SURFACE_GRAVITY
from atmosphere.radiative_convective import thermodynamics as th
from atmosphere.radiative_convective.run import MOON_12

HERE = Path(__file__).resolve().parent
SCHEMA = 'terluna.research.atmospheric-co2/1'
DESIGN_CO2_PPM = 400.0                  # the value the atmosphere and climate runs use
EARTH_TODAY_PA = 41.0                   # about 405 ppm at one standard atmosphere, the reference
SECONDS_PER_YEAR = 3.15576e7

# Farquhar-von Caemmerer-Berry C3 leaf model. Bernacchi et al. (2001) give Gamma*, Kc (umol/mol) and
# Ko (mmol/mol) as exp(c - dH/RT), measured near 100 kPa, so they convert to partial pressures there.
# Capacities scale with the mean activation energies of Medlyn et al. (2002); respiration with
# Bernacchi's. Only ratios are reported, so the typical 25 C capacities chosen here matter little.
R_KJ = 0.008314
MEASUREMENT_PA = 100.0e3
KINETICS = {'gamma_star': (19.02, 37.83, 1e-6), 'kc': (38.05, 79.43, 1e-6), 'ko': (20.30, 36.38, 1e-3)}
VCMAX25, JMAX25, RD_FRACTION = 60.0, 100.0, 0.015          # umol m-2 s-1; respiration as a share of Vcmax
EA_VCMAX, EA_JMAX, EA_RD = 65.0, 50.0, 46.39                 # kJ/mol
CI_CA = 0.7                                                  # internal to ambient CO2, a typical C3 value

# Dessert et al. (2003), basaltic catchments: CO2 consumed (mol km-2 yr-1) = Rf * 323.44 * exp(0.0642 T)
# with runoff Rf in mm/yr and T in C. Silicate weathering takes two CO2 per Ca or Mg into bicarbonate;
# carbonate precipitation in the seas returns one, so the lasting sink is half the gross.
DESSERT_A, DESSERT_B = 323.44, 0.0642
NET_FRACTION = 0.5
GLASS_FACTOR = (1.0, 10.0)             # fine, glassy regolith weathers up to ~10x faster than catchment rock

# Earth's land carbon in plants and the top metre of soil (IPCC TAR WG1 Table 3.2), and accumulation
# rates (kg C m-2 yr-1): regrowing tropical forest above ground (Poorter et al. 2016), soil (Post & Kwon
# 2000), primary succession (Lichter 1998).
LAND_CARBON_KG_C_M2 = (15.0, 16.0)
ACCUMULATION = {'regrowing_tropical_forest': 0.31, 'soil_build_up': 0.033, 'primary_succession': 0.023}

# CaCO3 -> CaO + CO2 takes 178 kJ/mol at 25 C (standard enthalpies); real kilns need about twice that.
CALCINATION_KJ_MOL = 178.0

EVIDENCE = ('A literature synthesis with small calculations: plant responses and rates come from Earth '
            'experiments, catchments and inventories, applied to the Moon\'s design air (1.2 atm, O2 at '
            'Earth\'s partial pressure). The photosynthesis ratios come from a textbook leaf model at fixed '
            'CO2 ratio and light saturation, not from lunar plants. Weathering and biosphere rates are Earth '
            'ranges; lunar regolith lacks nitrogen and soils, and the ocean\'s uptake is not included. None of '
            'this is a carbon-cycle model or a measurement.')
READING_RULE = ('Use the partial pressures (Pa), not ppm, when comparing with plant data: O2 is at Earth\'s '
                'partial pressure, so the CO2 partial pressure is what plants respond to. Treat rates and '
                'times as orders of magnitude; gross weathering is the uptake into bicarbonate, net is what '
                'stays after carbonate forms in the seas.')


def design_air(co2_ppm=DESIGN_CO2_PPM):
    return th.earthlike_air(MOON_12['dry_pressure_pa'], co2_ppm)


def kinetics(t_c):
    """(Gamma*, Kc, Ko) in Pa at leaf temperature t_c (C)."""
    t = t_c + 273.15
    return tuple(math.exp(c - dh / (R_KJ * t)) * unit * MEASUREMENT_PA for c, dh, unit in KINETICS.values())


def assimilation(ci_pa, t_c, o2_pa):
    """Net light-saturated C3 assimilation (umol m-2 s-1): the lesser of the Rubisco- and RuBP-limited rates."""
    gamma, kc, ko = kinetics(t_c)
    scale = lambda ea: math.exp(ea / R_KJ * (1 / 298.15 - 1 / (t_c + 273.15)))
    vcmax, jmax = VCMAX25 * scale(EA_VCMAX), JMAX25 * scale(EA_JMAX)
    wc = vcmax * ci_pa / (ci_pa + kc * (1 + o2_pa / ko))
    wj = jmax * ci_pa / (4 * ci_pa + 8 * gamma)
    return min(wc, wj) * (1 - gamma / ci_pa) - RD_FRACTION * VCMAX25 * scale(EA_RD)


def relative_photosynthesis(ca_pa, t_c, o2_pa):
    """Assimilation at ambient CO2 ca_pa relative to Earth's today, both at the fixed ci/ca."""
    return assimilation(CI_CA * ca_pa, t_c, o2_pa) / assimilation(CI_CA * EARTH_TODAY_PA, t_c, o2_pa)


def _supply_limited(ca_pa, total_pa, conductance, t_c, o2_pa):
    """Assimilation where demand meets diffusive supply g (ca - ci)/P, at a fixed molar conductance."""
    lo, hi = 1e-3, ca_pa
    for _ in range(100):
        ci = 0.5 * (lo + hi)
        if assimilation(ci, t_c, o2_pa) > conductance * 1e6 * (ca_pa - ci) / total_pa:
            hi = ci
        else:
            lo = ci
    return assimilation(0.5 * (lo + hi), t_c, o2_pa)


def pressure_penalty(t_c, o2_pa, total_pa):
    """At a stomatal conductance giving ci/ca = 0.7 at Earth's 41 Pa and one atmosphere: assimilation at
    the design pressure relative to Earth's, and the CO2 partial pressure that restores Earth's rate."""
    earth_pa = 101325.0
    a_earth = assimilation(CI_CA * EARTH_TODAY_PA, t_c, o2_pa)
    g = a_earth * 1e-6 * earth_pa / (EARTH_TODAY_PA * (1 - CI_CA))
    ratio = _supply_limited(EARTH_TODAY_PA, total_pa, g, t_c, o2_pa) / a_earth
    lo, hi = EARTH_TODAY_PA, 2 * EARTH_TODAY_PA
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if _supply_limited(mid, total_pa, g, t_c, o2_pa) > a_earth:
            hi = mid
        else:
            lo = mid
    return ratio, 0.5 * (lo + hi)


def weathering_mol_km2_yr(runoff_mm_yr, t_c):
    """Gross CO2 consumption by basalt weathering (Dessert et al. 2003)."""
    return runoff_mm_yr * DESSERT_A * math.exp(DESSERT_B * t_c)


def results():
    air = design_air()
    p = air.surface_pressure_pa
    o2_pa = air.fractions['O2'] * p
    area = 4 * math.pi * MOON_RADIUS ** 2
    water = json.loads((ROOT / 'shared' / 'scenarios' / 'water.json').read_text())
    drainage = json.loads((ROOT / 'geography' / 'results' / 'drainage.json').read_text())
    land = (1 - water['surface_water_share']) * area
    kg_per_pa_m2 = th.MOLAR_MASS['CO2'] / air.molar_mass / MOON_SURFACE_GRAVITY
    per_pa_gt = kg_per_pa_m2 * area / 1e12
    co2_pa = DESIGN_CO2_PPM * 1e-6 * p
    inventory = co2_pa * kg_per_pa_m2 * area

    plants = {}
    for t_c in (25.0, 30.0, 35.0):
        gamma = kinetics(t_c)[0]
        ratio, restore = pressure_penalty(t_c, o2_pa, p)
        plants[f'{t_c:g}C'] = dict(
            compensation_gamma_star_pa=round(gamma, 2),
            relative_to_earth_today={f'{ca:g}': round(relative_photosynthesis(ca, t_c, o2_pa), 3)
                                     for ca in (18.0, 20.0, 28.0, 35.0, 41.0, 49.0, 60.0, 80.0)},
            design_pressure_rate_at_41pa=round(ratio, 3), partial_pressure_matching_earth_today_pa=round(restore, 1))

    runoff_mm = drainage['rivers']['discharge_to_sea_m3s'] * SECONDS_PER_YEAR / land * 1000.0
    weathering = []
    for runoff in (runoff_mm, 1000.0):
        for t_c in (20.0, 26.0, 30.0):
            base = weathering_mol_km2_yr(runoff, t_c)
            for factor in GLASS_FACTOR:
                gross = base * factor * th.MOLAR_MASS['CO2'] / 1e6 * land / 1e12           # Gt CO2 per year
                weathering.append(dict(runoff_mm_yr=round(runoff), temperature_c=t_c, fresh_glass_factor=factor,
                                       gross_mol_km2_yr=round(base * factor, -3), gross_gt_yr=round(gross, 2),
                                       net_gt_yr=round(gross * NET_FRACTION, 2),
                                       years_net_49_to_28_pa=round(21 * per_pa_gt / (gross * NET_FRACTION), -1)))
    land_carbon = [c * land * 44.0095 / 12.011 for c in LAND_CARBON_KG_C_M2]   # kg CO2 an Earth-like biosphere holds
    calcination_gw_per_gt = CALCINATION_KJ_MOL * 1e3 / th.MOLAR_MASS['CO2'] * 1e12 / SECONDS_PER_YEAR / 1e9

    return dict(
        schema=SCHEMA, evidence=EVIDENCE, reading_rule=READING_RULE,
        design_air=dict(pressure_pa=p, co2_ppm=DESIGN_CO2_PPM, co2_pa=round(co2_pa, 2), o2_pa=round(o2_pa),
                        ppm_per_pa=round(1e6 / p, 3), kg_co2_per_m2_per_pa=round(kg_per_pa_m2, 3),
                        gt_co2_per_pa=round(per_pa_gt, 1), inventory_kg_co2=float(f'{inventory:.3e}'),
                        land_share=round(1 - water['surface_water_share'], 3),
                        inventory_kg_co2_per_m2_land=round(inventory / land, 1)),
        photosynthesis=plants,
        biosphere=dict(earth_like_land_carbon_kg_co2=[float(f'{x:.3e}') for x in land_carbon],
                       share_of_inventory=[round(x / inventory, 2) for x in land_carbon],
                       years_to_accumulate={k: round(sum(LAND_CARBON_KG_C_M2) / 2 / v) for k, v in ACCUMULATION.items()}),
        weathering=dict(runoff_from=f"geography/results/drainage.json rivers ({drainage['rivers']['discharge_to_sea_m3s']} m3/s over "
                                    f"{round(1 - water['surface_water_share'], 2)} of the surface)",
                        cases=weathering),
        return_step=dict(calcination_kj_per_mol=CALCINATION_KJ_MOL, heat_gw_per_gt_co2_per_year=round(calcination_gw_per_gt)))


def run(out=HERE / 'results'):
    out.mkdir(parents=True, exist_ok=True)
    data = results()
    data['producer'] = dict(file='research/studies/atmospheric_co2/run.py',
                            sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16])
    (out / 'atmospheric_co2.json').write_text(json.dumps(data, indent=1) + '\n')
    return data


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--out', type=Path, default=HERE / 'results')
    args = parser.parse_args(argv)
    data = run(args.out)
    d = data['design_air']
    print(f"{d['co2_ppm']:g} ppm = {d['co2_pa']} Pa; {d['gt_co2_per_pa']} Gt CO2 per Pa; inventory {d['inventory_kg_co2']:.2e} kg")
    for t, p in data['photosynthesis'].items():
        print(t, p['relative_to_earth_today'], 'design-pressure rate', p['design_pressure_rate_at_41pa'],
              'matching Pa', p['partial_pressure_matching_earth_today_pa'])
    return 0


if __name__ == '__main__':
    sys.exit(main())
