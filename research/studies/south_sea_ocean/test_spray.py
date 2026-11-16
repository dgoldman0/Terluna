import json
import os

import numpy as np
import pytest

from research.studies.south_sea_ocean import spray as SP
from shared.constants import JULIAN_DAY, MOON_RADIUS, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY


def test_the_source_reproduces_the_marine_notes_number_fluxes():
    # aerosol_marine.md, section 1.3: Gong (2003) times the warmth weight at 26 C, over r80 of 0.035-20 um
    number, mass, radius = SP.bins()
    table = {1: 146.0, 2: 1.56e3, 3: 6.2e3, 5: 3.54e4, 10: 3.76e5, 15: 1.5e6, 20: 4.0e6}
    for u, flux in table.items():
        assert 1.373 * u ** 3.41 * SP.sst_weight(26.0) * number.sum() == pytest.approx(flux, rel=0.01)
    assert mass.sum() > 0 and (np.diff(radius) > 0).all()
    assert mass[-3:].sum() / mass.sum() > 0.5                  # most of the salt rides on drops above 7 um of r80


def test_the_warmth_weight_is_held_above_the_fits_data():
    assert SP.sst_weight(26.0) == pytest.approx(1.453, abs=0.001)
    assert SP.sst_weight(31.4) == SP.sst_weight(30.0) == pytest.approx(2.13, abs=0.001)


def test_the_lunar_wind_is_carried_to_earths_at_equal_energy_into_the_waves():
    # aerosol_marine.md, section 1.5: friction velocities and their cube ratio at 5 and 10 m/s; section 7.1: U_eq
    assert SP.friction_velocity(5.0, STANDARD_GRAVITY, 1.18) == pytest.approx(0.168, abs=0.001)
    assert SP.friction_velocity(5.0, MOON_SURFACE_GRAVITY, 1.41) == pytest.approx(0.201, abs=0.001)
    assert SP.friction_velocity(10.0, MOON_SURFACE_GRAVITY, 1.41) == pytest.approx(0.484, abs=0.002)
    assert SP.equivalent_wind(5.0, MOON_SURFACE_GRAVITY, 1.41) == pytest.approx(6.13, abs=0.02)
    assert SP.equivalent_wind(10.0, MOON_SURFACE_GRAVITY, 1.41) == pytest.approx(12.71, abs=0.03)
    assert SP.equivalent_wind(7.0, STANDARD_GRAVITY, SP.EARTH_AIR) == pytest.approx(7.0, abs=0.02)


def test_the_laws_stand_in_their_order():
    u = np.array([2.0, 4.0])
    low, central, high = (SP.wind_factor(u, k) for k in ('low', 'central', 'high'))
    assert (low < central).all() and (central < high).all()
    assert low / central == pytest.approx(0.34) and high / central == pytest.approx(2.46)
    assert (SP.wind_factor(u, 'earth_law') < central).all()
    gentle = SP.wind_factor(u, 'gentle_wind_law')
    assert gentle[0] / central[0] > gentle[1] / central[1] > 1.0         # the gentler law gives more in light winds


def test_drops_settle_a_sixth_as_fast_at_lunar_gravity():
    assert SP.settling(5.0) == pytest.approx(5.9e-4, rel=0.02)
    assert SP.settling(10.0) == pytest.approx(2.4e-3, rel=0.03)
    earth = SP.settling(5.0, gravity=STANDARD_GRAVITY, pressure_pa=101325.0, temperature_k=293.0)
    assert SP.settling(5.0) / earth == pytest.approx(0.165, abs=0.002)


def test_the_grid_holds_the_whole_sphere_and_finds_its_cells():
    glat = np.array([60.0, 20.0, -20.0, -60.0])
    edges = SP.latitude_edges(glat)
    assert edges.tolist() == [90.0, 40.0, 0.0, -40.0, -90.0]
    assert SP.cell_area_m2(glat, 8).sum() * 8 == pytest.approx(4.0 * np.pi * MOON_RADIUS ** 2)
    x = SP.unit(np.array([89.0, 39.0, -1.0, -89.0]), np.array([0.0, 100.0, 359.0, 180.0]))
    i, j = SP.cell_of(x, edges, 0.0, 8)
    assert i.tolist() == [0, 1, 2, 3] and j.tolist() == [0, 2, 0, 4]


def test_a_parcel_moves_with_the_wind_and_crosses_the_pole():
    x = SP.unit(np.array([0.0]), np.array([10.0]))
    day = 86400.0
    east = SP.step(x, np.array([2.0]), np.array([0.0]), day)
    assert np.degrees(np.arctan2(east[0, 1], east[0, 0])) == pytest.approx(10.0 + np.degrees(2.0 * day / MOON_RADIUS),
                                                                           abs=0.01)
    assert np.linalg.norm(east) == pytest.approx(1.0)
    y = SP.unit(np.array([89.9]), np.array([0.0]))
    for _ in range(20):
        y = SP.step(y, np.array([0.0]), np.array([1.0]), 3600.0)             # northward at the start
    assert y[0, 2] < 1.0 and np.isfinite(y).all() and np.linalg.norm(y) == pytest.approx(1.0)


def still_air(rain_kg_m2_s, fate, water=True, writes=40):
    glat = np.array([60.0, 20.0, -20.0, -60.0])
    glon = np.arange(0.0, 360.0, 45.0)
    shape = (writes, 4, 8)
    emission = np.zeros((writes, 1, 2))
    emission[0, 0] = (3.0, 1.0)
    surface = np.full((4, 8), water)
    return SP.follow(emission, np.array([2]), np.array([5]), np.zeros(shape), np.zeros(shape),
                     np.full(shape, rain_kg_m2_s), surface, glat, glon, np.array([0.0, 0.002]), 1.4, fate,
                     max_years=400), emission


def test_every_gram_comes_down_and_in_still_air_where_it_rose():
    fate = dict(scavenging_ratio=300.0, deposition_m_s=0.0005)
    f, emission = still_air(2.0e-5, fate)
    total = f['wet'] + f['dry']
    assert total.sum() + f['still_airborne_kg'] == pytest.approx(emission.sum())
    assert total[2, 5] == pytest.approx(emission.sum()) and f['travel_m'] == pytest.approx(0.0, abs=1e-6)
    # each size leaves at its own rate: rain at W P / (rho H), the ground at (settling + deposition) / H
    k_wet = 300.0 * 2.0e-5 / (1.4 * SP.DEPTH_WATER_M)
    k_dry = (np.array([0.0, 0.002]) + 0.0005) / SP.DEPTH_WATER_M
    share = k_wet / (k_wet + k_dry)
    assert f['by_bin'][:, 0] / f['by_bin'].sum(axis=1) == pytest.approx(share, rel=2e-3)
    days = (emission[0, 0] / (k_wet + k_dry)).sum() / emission.sum() / JULIAN_DAY
    assert f['airborne'].sum() / emission.sum() / JULIAN_DAY == pytest.approx(days, rel=2e-3)


def test_dry_air_over_land_keeps_salt_aloft_longer_than_rain_over_water():
    fate = dict(scavenging_ratio=300.0, deposition_m_s=0.0005)
    rainy, emission = still_air(2.0e-5, fate, water=True)
    dry_land, _ = still_air(0.0, fate, water=False)
    assert dry_land['wet'].sum() == 0.0
    assert dry_land['airborne'].sum() > 3.0 * rainy['airborne'].sum()


def test_neighbours_join_across_the_seam_and_stop_at_the_poles():
    mask = np.zeros((3, 6), bool)
    mask[0, 0] = True
    near = SP.neighbours(mask)
    assert near[0, 5] and near[1, 5] and near[1, 1] and near[0, 1]
    assert not near[2].any() and near.sum() == 6


def test_the_product_is_current_and_its_salt_is_all_accounted_for():
    out = json.loads(SP.OUT.read_text())
    assert out['schema'] == SP.SCHEMA
    for path, digest in out['producer']['files'].items():
        assert SP.sha256(SP.ROOT / path) == digest, path
    assert out['chloride_share_of_salt'] == pytest.approx(0.5503, abs=0.0002)
    s = out['source']
    assert s['low']['salt_mt_yr'] < s['earth_law']['salt_mt_yr'] < s['central']['salt_mt_yr'] < s['high']['salt_mt_yr']
    assert s['central']['salt_g_m2_yr'] == pytest.approx(s['central']['salt_mt_yr'] * 1e12 / (out['sea']['area_km2'] * 1e6))
    assert s['central']['salt_g_m2_yr'] < 25.0 / 5.0           # well under the Earth ocean's mean of the lookup
    assert sum(out['sizes']['share_of_salt']) == pytest.approx(1.0)
    assert sum(b['share_of_salt'] for b in out['wind']['by_speed'].values()) == pytest.approx(1.0)
    days = []
    for f in out['fate']['cases'].values():
        assert f['deposited_over_emitted'] == pytest.approx(1.0, abs=1e-9)
        assert f['emitted_mt_yr'] == pytest.approx(s['central']['salt_mt_yr'])
        share = f['share_falling']
        assert (share['back_on_the_south_sea'] + sum(share['on_the_other_seas'].values())
                + share['on_land_and_lakes']) == pytest.approx(1.0)
        c = f['chloride_kg_ha_yr']
        assert c['on_the_south_sea'] > c['land_beside_the_sea'] > c['land_mean'] > c['desert_land_mean'] > 0
        days.append(f['mean_days_in_the_air'])
    assert days == sorted(days, reverse=True)                  # slow, central, fast
    lost = out['against_the_sea']
    assert lost['leaving_the_sea_for_good_mt_yr']['central'] < s['central']['salt_mt_yr']
    assert min(lost['years_to_carry_off_a_thousandth'].values()) > 1.0e4


@pytest.mark.skipif(os.environ.get('TERLUNA_SOUTH_SEA_REGENERATE') != '1' or not SP.RECORD.exists(),
                    reason='needs the 6.5-hour record on the research drive and netCDF4 (the ocean environment), '
                           'two minutes; set TERLUNA_SOUTH_SEA_REGENERATE=1')
def test_the_spray_regenerates():
    new = json.loads(json.dumps(SP.build()))
    old = json.loads(SP.OUT.read_text())
    for case in ('slow', 'central', 'fast'):
        new['fate']['cases'][case].pop('seconds'), old['fate']['cases'][case].pop('seconds')
    assert new == old
