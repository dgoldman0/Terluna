import hashlib
import json
import math
import os

import numpy as np
import pytest

from research.studies.river_deltas import physics as P
from research.studies.river_deltas import run as R

G_MOON, G_EARTH = R.MOON['g'], R.EARTH['g']


@pytest.fixture(scope='module')
def product():
    return json.loads(R.OUT.read_text())


# --- Physics ----------------------------------------------------------------------------------------------------------

def test_water_density_matches_tanaka_at_its_maximum_and_at_25_c():
    assert P.water_density(3.98) == pytest.approx(999.975, abs=1e-3)
    assert P.water_density(25.0) == pytest.approx(997.05, abs=0.02)


def test_settling_is_stokes_for_fine_grains_and_scales_with_gravity():
    d, nu = 2e-6, 1e-6
    stokes = 1.65 * G_EARTH * d * d / (18.0 * nu)
    assert P.settling_velocity(d, G_EARTH, 2650.0, nu) == pytest.approx(stokes, rel=1e-3)
    ratio = P.settling_velocity(d, G_MOON, 2650.0, nu) / P.settling_velocity(d, G_EARTH, 2650.0, nu)
    assert ratio == pytest.approx(G_MOON / G_EARTH, rel=1e-3)
    # Medium sand on Earth settles at about 2.5 cm/s (Ferguson & Church 2004's curve for natural grains).
    assert 0.02 < P.settling_velocity(2e-4, G_EARTH, 2650.0, 1e-6) < 0.03


def test_parker_gravel_reproduces_the_published_coefficients():
    q, d = 100.0, 0.03
    h, b, s = P.parker_gravel(q, d, G_EARTH)
    qhat = q / (math.sqrt(G_EARTH * d) * d * d)
    assert h == pytest.approx(0.382 * q ** 0.4 / G_EARTH ** 0.2)
    assert b == pytest.approx(4.63 * q ** 0.4 / G_EARTH ** 0.2 * qhat ** 0.0667)
    assert s == pytest.approx(0.101 * qhat ** -0.344)


def test_threshold_channel_holds_the_bank_stress_and_scales_as_stated():
    loose = dict(depth=np.array([2.0]), width=np.array([100.0]), velocity=np.array([1.0]), slope=np.array([1e-4]),
                 stress=np.array([4.0]))
    ch = P.threshold_channel(loose, np.array([0.1]), np.array([10.0]), np.array([1e9]))
    k = 12.0 / 4.0
    assert ch['stress'][0] == pytest.approx(12.0)
    assert ch['depth'][0] == pytest.approx(2.0 * k)
    assert ch['width'][0] == pytest.approx(100.0 * k ** -1.5)
    assert ch['width'][0] * ch['depth'][0] * ch['velocity'][0] == pytest.approx(200.0)    # continuity


def test_a_cohesive_bank_in_pascals_makes_lunar_channels_six_times_deeper_at_one_slope():
    """At one slope and one bank stress, depth is tau / (rho g S): six times Earth's at lunar gravity, while loose
    banks, whose stress falls with gravity, keep the depth."""
    tau, s = 12.0, 1e-4
    assert (tau / (1000.0 * G_MOON * s)) / (tau / (1000.0 * G_EARTH * s)) == pytest.approx(G_EARTH / G_MOON)
    nu, d, r = 1e-6, 2e-4, 1.7
    hs_moon = P.sand_bankfull_shields(P.particle_reynolds(d, G_MOON, r, nu)) * r * d
    hs_earth = P.sand_bankfull_shields(P.particle_reynolds(d, G_EARTH, r, nu)) * r * d
    assert 1.0 < hs_moon / hs_earth < (G_EARTH / G_MOON) ** 0.25 + 1e-9   # only through Rep^-1/2


def test_longshore_transport_holds_gravity_to_three_fifths():
    q_moon = P.longshore_max(1.0, 10.0, G_MOON, 2650.0)
    q_earth = P.longshore_max(1.0, 10.0, G_EARTH, 2650.0)
    assert q_moon / q_earth == pytest.approx((G_MOON / G_EARTH) ** 0.6)
    assert P.ANGLE_FACTOR == pytest.approx(0.4685, abs=2e-4)
    # A metre of 8 s swell on Earth moves of order a million cubic metres a year along a straight coast.
    assert 0.3 < 2.0 * q_earth * (8.0 / 10.0) ** 0.2 * R.YEAR_S / 1600.0 / 1e6 < 5.0


def test_wave_power_inverts_to_its_height():
    h = P.wave_height_from_power(1000.0 * G_MOON ** 2 * 1.5 ** 2 * 10.0 / (64.0 * math.pi), 10.0, G_MOON)
    assert h == pytest.approx(1.5)


def test_fetch_waves_at_lunar_gravity_come_near_the_wave_study():
    """Fully developed seas under the seas' median 2.6 m/s wind against SWAN's nearside mean of 1.39 m."""
    h, tp = P.fetch_waves(2.6, 1e6, G_MOON)
    assert h == pytest.approx(0.2433 * 2.6 ** 2 / G_MOON)
    assert 0.6 < h / 1.39 < 1.2
    assert 10.0 < tp < 15.0


def test_plunge_depth_and_rossby_radius():
    assert P.plunge_depth(2.0, 0.01) == pytest.approx((4.0 / (0.01 * 0.25)) ** (1.0 / 3.0))
    omega = R.MOON['omega']
    assert P.rossby_radius(1e-3, 4.0, 30.0, omega) == pytest.approx(math.sqrt(4e-3) / omega)


def test_fresh_river_water_needs_earths_sediment_load_to_plunge_into_seawater():
    """Mulder & Syvitski (1995): 36-43 kg/m3 of sediment for river water to sink into the sea."""
    assert 36.0 < R.seawater_plunge_threshold() < 43.0


def test_flocculation_share_rises_with_hardness_and_life():
    lo, hi = R.CCC_DIVALENT_MMOL
    assert R.flocculated_share(np.array([lo / 2.0]), np.array([0.0]), 0.0)[0] == 0.0
    assert R.flocculated_share(np.array([hi * 2.0]), np.array([0.0]), 0.0)[0] == 1.0
    assert R.flocculated_share(np.array([0.1]), np.array([0.0]), 0.6)[0] == pytest.approx(0.6)


def test_textures_run_from_sand_to_clay_and_regolith_has_little_clay():
    clay = {k: R.grain_fractions(k)[R.SIZES < 2e-6].sum() for k in R.SOILS}
    assert clay['sand'] < clay['sandy_loam'] < clay['loam'] < clay['clay_loam'] < clay['clay']
    assert clay['regolith'] < 0.02
    f = R.grain_fractions('regolith')
    assert 0.8 < f[(R.SIZES > 1e-5) & (R.SIZES < 1e-3)].sum() < 0.95      # Carrier: 90% between 10 um and 1 mm
    assert R.bank_critical_pa(R.grain_fractions('clay'), 1.0) > R.bank_critical_pa(R.grain_fractions('sand'), 1.0)


# --- The product ------------------------------------------------------------------------------------------------------

def test_product_binding(product):
    assert product['schema'] == R.SCHEMA
    for rel, digest in product['producer']['files'].items():
        assert hashlib.sha256((R.ROOT / rel).read_bytes()).hexdigest() == digest, rel
    for key in ('drainage', 'groundwater', 'cold_waters', 'atlas_bodies'):
        entry = product['inputs'][key]
        assert R.sha256(entry['path']) == entry['sha256'], key
    assert product['evidence'] and product['reading_rule']


def test_mouths_and_their_discharge_balance(product):
    c = product['checks']
    kinds = product['waters']['mouths_by_kind']
    assert sum(kinds.values()) == c['river_mouths']
    assert abs(c['common_level_discharge_m3s'] / c['drainage_discharge_to_sea_m3s'] - 1.0) < 0.02
    for row in c['largest_rivers'][:11]:
        assert row['screened_m3s'] == pytest.approx(row['drainage_m3s'], rel=0.01)


def test_shares_add_up(product):
    for state in ('young/regolith', 'mature/loam'):
        s = product['states'][state]['central']['all']
        assert sum(s['dominance_by_load'].values()) == pytest.approx(1.0, abs=1e-3)
        plume = s['plume_mean_flow_by_load']
        assert plume['floating'] + plume['level'] + plume['plunging'] == pytest.approx(1.0, abs=1e-3)
        assert plume['to_floor'] <= plume['plunging'] + 1e-9
        assert 0.0 <= s['keep_share_by_load'] <= s['keep_with_underflow_by_load'] <= 1.0


def test_young_moon_has_more_rock_and_weaker_banks_than_the_mature(product):
    young = product['states']['young/regolith']['central']['all']
    mature = product['states']['mature/loam']['central']['all']
    assert young['rock_mt_yr'] > mature['rock_mt_yr']
    assert young['bank_critical_pa'][1] < mature['bank_critical_pa'][1]
    assert young['pattern_by_mouth']['braided'] > mature['pattern_by_mouth']['braided']


def test_earth_control_reproduces_the_deltas_types(product):
    e = product['earth_deltas']
    for name in ('mississippi', 'sao_francisco', 'rhone_geneva'):
        assert e[name]['screened_type'] == e[name]['expected_type'], name
    for name, row in e.items():
        assert row['type_on_observed_channel'] == row['expected_type'], name
        assert row['screened_plume'] == row['expected_plume'], name
    assert e['mississippi']['bar_mode_on_observed_channel'] < 1.5


@pytest.mark.skipif(not os.environ.get('TERLUNA_RIVER_DELTAS_REGENERATE'), reason='set to rerun the screen (20 s)')
def test_regenerates(product):
    fresh = R.run()
    assert fresh['checks'] == product['checks']
    assert fresh['states']['mature/loam']['central']['all'] == product['states']['mature/loam']['central']['all']
