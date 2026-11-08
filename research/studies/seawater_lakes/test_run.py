import json
import math

import pytest

from research.studies.seawater_lakes import run as R


def test_seawater_density_matches_eos80_check_values():
    """UNESCO's check values for EOS-80 at one atmosphere (Fofonoff and Millard 1983)."""
    assert R.seawater_density(0.0, 5.0) == pytest.approx(999.96675, abs=1e-4)
    assert R.seawater_density(35.0, 5.0) == pytest.approx(1027.67547, abs=1e-4)
    assert R.seawater_density(35.0, 25.0) == pytest.approx(1023.34306, abs=1e-4)


def test_reference_salinity_converts_to_practical():
    assert R.practical_from_reference(35.16504) == pytest.approx(35.0)


def test_reduced_gravity_scales_with_gravity():
    rho_s, rho_f = R.seawater_density(35.0, 29.5), R.seawater_density(0.1, 29.5)
    ratio = R.reduced_gravity(R.G_MOON, rho_s, rho_f) / R.reduced_gravity(R.G_EARTH, rho_s, rho_f)
    assert ratio == pytest.approx(R.G_MOON / R.G_EARTH)
    assert R.G_MOON == pytest.approx(1.6242, abs=1e-4)


def test_arrest_discharge_is_densimetric_froude_one():
    gp, h = 0.042, 0.5
    q = R.arrest_discharge_per_width(gp, h)
    assert (q / h) / math.sqrt(gp * h) == pytest.approx(1.0)
    moon = R.arrest_discharge_per_width(gp, h)
    earth = R.arrest_discharge_per_width(gp * R.G_EARTH / R.G_MOON, h)
    assert moon / earth == pytest.approx(math.sqrt(R.G_MOON / R.G_EARTH))


def test_pool_drop_returns_the_jet_speed_and_scales_inversely_with_gravity():
    dz = R.pool_drop_for_speed(1.0, R.G_MOON)
    assert math.sqrt(2.0 * R.G_MOON * dz) == pytest.approx(1.0)
    assert dz / R.pool_drop_for_speed(1.0, R.G_EARTH) == pytest.approx(R.G_EARTH / R.G_MOON)


def test_pool_evaporation_per_lunar_day_to_a_year():
    per_lunar_day_mm = 62.0 + 5.6
    expected = per_lunar_day_mm / R.LUNAR_DAY_D * 365.25 * 1e-3
    assert R.pool_evaporation_m_yr(62.0, 5.6) == pytest.approx(expected)
    assert R.pool_evaporation_m_yr(62.0, 5.6, 2.0) == pytest.approx((124.0 + 5.6) / R.LUNAR_DAY_D * 0.36525)


def test_a_calcium_bicarbonate_feed_needs_no_bleed_or_dose():
    b = R.feed_budget(dict(alkalinity=1.68, ca=0.84, mg=0.0, na_k=0.0), 0.7, 20.0, 35.9)
    assert b['bleed_share_of_inflow'] == 0.0
    assert b['calcium_chloride_dose_t_yr_km2'] == 0.0
    assert b['lake_calcium_spent_yr'] is None
    assert b['calcite_from_feed_t_yr_km2'] == pytest.approx(0.84 * 0.7e6 * R.M_CACO3 * 1e-6)


def test_a_magnesium_feed_sets_the_bleed_and_the_calcium_draw():
    mg = 0.5
    b = R.feed_budget(dict(alkalinity=1.68 + 2 * mg, ca=0.84, mg=mg, na_k=0.0), 0.7, 20.0, 35.9)
    assert b['excess_alkalinity_meq_kg'] == pytest.approx(2 * mg)
    assert b['bleed_share_of_inflow'] == pytest.approx(mg / (R.TOLERANCE * R.SEAWATER['mg'] * 1e3))
    # the lake's calcium (10.28 mmol/kg over 20 m) falls by the excess over two each year
    assert b['lake_calcium_spent_yr'] == pytest.approx(R.SEAWATER['ca'] * 1e3 * 20.0 / (0.5 * 2 * mg * 0.7))
    assert b['salt_top_up_with_halite_returned_t_yr_km2'] < b['salt_top_up_t_yr_km2']


def test_liner_seepage_scales_with_gravity():
    moon = R.liner_seepage_m_yr(1e-10, R.G_MOON, 20.0, 1.0)
    earth = R.liner_seepage_m_yr(1e-10, R.G_EARTH, 20.0, 1.0)
    assert moon / earth == pytest.approx(R.G_MOON / R.G_EARTH)
    assert earth == pytest.approx(1e-10 * 21.0 * R.YEAR_S)


def test_wind_setup_is_inverse_in_gravity():
    m = R.wind_setup_m(10.0, 1e4, 20.0, R.G_MOON, 1022.0)
    e = R.wind_setup_m(10.0, 1e4, 20.0, R.G_EARTH, 1022.0)
    assert m / e == pytest.approx(R.G_EARTH / R.G_MOON)


def test_drift_check_parser():
    text = ('X\n  year 10: land: ts 300 evap -1.5 pr 2.2 alb 0.2 | sea: ts 303 evap -2.7 pr 1.6 alb 0.07\n'
            '  year 11: land: ts 300 evap -1.7 pr 2.0 alb 0.2 | sea: ts 303 evap -2.9 pr 1.8 alb 0.07\n')
    d = R.parse_drift_check(text)
    assert d['years'] == [10, 11]
    assert d['land_evaporation_mm_day'] == pytest.approx(1.6)
    assert d['sea_rain_mm_day'] == pytest.approx(1.7)


def test_the_recorded_turned_years_match_their_drift_check():
    """The carried means equal the drift check's on the research drive, where the drive is mounted."""
    if not R.DRIFT_CHECK.exists():
        pytest.skip('the research drive is not mounted')
    d = R.read_drift_check(R.DRIFT_CHECK)
    for k, v in R.TURNED_REFIT_YEARS.items():
        if k != 'source':
            assert d[k] == pytest.approx(v, abs=5e-4), k


def test_salt_per_cubic_kilometre_agrees_with_the_sea_chemistry():
    """The sea chemistry gives 36 Mt of salt for each km3 of seawater."""
    rho = R.seawater_density(R.practical_from_reference(R.SEAWATER_G_KG), R.LAKE_T_C)
    assert rho * R.SEAWATER_G_KG * 1e-3 == pytest.approx(36.0, abs=0.2)


def test_product_regenerates():
    """The committed product is what the runner writes from its inputs now."""
    built = json.loads(json.dumps(R.build()))
    assert json.loads(R.OUT.read_text()) == built
