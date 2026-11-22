"""Tests of the cement sizing: its own checks, its copies of this study's figures, and the figures the page states."""
from __future__ import annotations

import json
import math

from research.studies.rock_return import carriers as CA
from research.studies.rock_return import cement as CE
from research.studies.rock_return import run as RN
from shared.constants import MOON_SURFACE_GRAVITY


def test_the_modules_own_checks_pass():
    CE.self_check()


def test_the_figures_copied_from_the_study_are_the_studys():
    product = json.loads(RN.OUT.read_text())
    for case, gt in CE.ROCK_GT_YR.items():
        assert math.isclose(gt, product['job']['gt_yr'][case], rel_tol=1e-4)
    assert math.isclose(CE.RETURN_G_M2_YR, 1000.0 * product['job']['kg_m2_yr_central'], rel_tol=1e-4)
    assert math.isclose(CE.MEAN_LOWERING_M_MYR, product['job']['mean_lowering_m_per_myr'], rel_tol=1e-3)
    stone = product['staying_put']['the_rock_that_leaves_in_solution_mt_yr']
    assert all(math.isclose(CE.SEAS_STONE_MT_YR[k], stone[k], rel_tol=1e-5) for k in CE.SEAS_STONE_MT_YR)
    nitrogen = product['staying_put']['the_carriers_nitrogen_mt_yr']
    assert all(math.isclose(a, b, rel_tol=5e-3) for a, b in zip(CE.CARRIERS_N_MT_YR, nitrogen))
    assert math.isclose(CE.LOOSE_SOIL_KG_M3, product['setting']['loose_soil_kg_m3'])
    assert math.isclose(CE.MOON_G, MOON_SURFACE_GRAVITY, rel_tol=5e-4)
    lime = CE.carried_fractions(CE.ROCK_GT_YR['central'])
    assert math.isclose(lime['calcite'] + lime['magnesite'], product['staying_put']['lime_over_the_particles_by_mass']['central'], rel_tol=1e-3)


def test_the_urea_route_is_the_carriers_modules():
    assert math.isclose(CA.calcite_by_urea_kg(CE.carbon_for(100.087, 'urea_nitrogen')), 100.087, rel_tol=1e-9)


def test_what_the_write_up_says_of_cement():
    # rim rock: cohesion for each metre of thickness at lunar gravity, dry and saturated, 35 to 60 degrees
    dry = [CE.cohesion_needed_pa(1.0, s, rho) / 1e3 for s in (35.0, 45.0, 60.0) for rho in (1300.0, 1800.0)]
    wet = [CE.cohesion_needed_pa(1.0, s, rho, saturated=True) / 1e3 for s in (35.0, 45.0, 60.0) for rho in (1800.0, 2100.0)]
    assert 0.16 < min(dry) < 0.18 and 0.83 < max(dry) < 0.85 and 0.86 < min(wet) and max(wet) < 1.23
    assert 15.0 < CE.stable_thickness_m(10.0e3, 45.0, 1800.0) < 17.0                       # 10 kPa holds 16 m dry at 45 degrees
    assert 2.5 < CE.stable_thickness_m(10.0e3, 45.0, 1800.0, gravity=CE.EARTH_G) < 2.9     # against 2.7 m on Earth
    # the stock: what weathering frees in place and what the load carries
    here = CE.in_place_supply(CE.weathering_eq_m2_yr())
    assert 10.7 < here['calcite'] < 10.9 and 0.073 < here['calcite'] / CE.RETURN_G_M2_YR < 0.075
    assert 0.062 < CE.goethite_equivalent(CE.HIGHLAND_FEO) < 0.064
    # the need: a cement of 3% of the return
    assert 4.3 < CE.cement_need_g_m2_yr(0.03) < 4.5
    # water takes carbonate from wet ground: a 3% cement in the top 0.1 m under the humid land's water
    soil = CE.calcite_solubility_g_m3(CE.SOIL_PCO2_ATM)
    air = CE.calcite_solubility_g_m3(CE.AIR_PCO2_ATM)
    years = sorted(CE.cement_lifetime_yr(0.03, 0.1, rho, mm, sol) for rho in (1300.0, 1800.0) for mm in CE.HUMID_WATER_LEAVING_MM for sol in (air, soil))
    assert 20.0 < years[0] < 26.0 and 95.0 < years[-1] < 150.0
    # the binder for thinly spread rock
    low, high = (f * CE.RETURN_G_M2_YR for f in CE.ORGANIC_BINDER_SHARE)
    assert 1.4 < low < 1.5 and 7.2 < high < 7.4
