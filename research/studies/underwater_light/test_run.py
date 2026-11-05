import hashlib
import json

import numpy as np
import pytest

from illumination.water_column import model as domain
from research.studies.underwater_light import optics as O
from research.studies.underwater_light import run as R
from shared.provenance import constants_changed

PRODUCT = R.OUT


def product():
    if not PRODUCT.exists():
        pytest.skip('results/underwater_light.json has not been built')
    return json.loads(PRODUCT.read_text())


# ==== Fresh water's optics =========================================================================================
def test_fresh_water_scattering_is_pure_waters():
    """Seawater at the table's 38.4 g/kg scatters 30% more than pure water; the 0.14 g/kg seas are pure water's."""
    assert O.scattering_factor(38.4) == pytest.approx(1.0)
    assert O.scattering_factor(0.0) == pytest.approx(1 / 1.30)
    w = np.array([400.0, 500.0, 600.0])
    _, b = O.pure_water(w)
    _, b_table = domain.seawater(w)
    assert np.allclose(b / b_table, (1 + 0.30 * 0.14 / 38.4) / 1.30)


def test_pure_water_absorption_keeps_pope_and_fry_and_replaces_the_ultraviolet():
    w = np.arange(420.0, 701.0, 5.0)
    assert np.allclose(O.pure_water(w)[0], domain.seawater(w)[0])
    assert O.pure_water(np.array([300.0]))[0][0] == pytest.approx(0.0147)          # IOCCG Table 1.1
    assert O.pure_water(np.array([350.0]))[0][0] == pytest.approx(0.0071)
    assert domain.seawater(np.array([300.0]))[0][0] > 9 * 0.0147                   # Smith and Baker's bound in the table
    assert np.all(O.pure_water(np.array([300.0, 350.0]), uv_absorption=False)[0] == 0)


def test_steeper_dissolved_organic_slope_only_below_the_knee():
    w = np.array([300.0, 340.0, 360.0, 440.0])
    base = O.dissolved_organic(w, 0.1)
    steep = O.dissolved_organic(w, 0.1, uv_slope=0.025)
    assert np.allclose(base[2:], steep[2:]) and np.all(steep[:2] > base[:2])
    assert base[3] == pytest.approx(0.1)


# ==== Surface and column ===========================================================================================
def test_fresnel_limits():
    assert O.fresnel_reflectance(0.0) == pytest.approx(((O.WATER_INDEX - 1) / (O.WATER_INDEX + 1)) ** 2, rel=1e-6)
    assert O.fresnel_reflectance(np.radians(89.99)) > 0.99
    # An isotropic sky over flat water: about 6.6% reflected.
    assert O.diffuse_transmission(8) == pytest.approx(0.934, abs=0.002)
    assert O.diffuse_transmission(8) == pytest.approx(O.diffuse_transmission(64), abs=1e-6)


def test_streams_share_an_isotropic_skys_irradiance():
    zenith, share = O.sky_streams(8)
    assert share.sum() == pytest.approx(1.0)
    mu = np.cos(np.radians(zenith))
    assert np.sum(share * mu) == pytest.approx(2 / 3, rel=1e-9)        # irradiance-weighted mean cosine, 2/3


def test_overhead_beam_follows_lees_attenuation():
    w = np.array([450.0, 550.0, 650.0])
    depth = np.array([0.0, 5.0, 20.0])
    col = O.Column(w, depth, chlorophyll_mg_m3=1.0, cdom_440=0.1)
    k = domain.diffuse_attenuation(col.iop['a'], col.iop['bb'], 0.0)
    expected = (1 - O.fresnel_reflectance(0.0)) * np.exp(-np.outer(k, depth))
    assert np.allclose(col.direct(0.0), expected)
    assert np.allclose(col.diffuse[:, 0], O.diffuse_transmission(8))
    # The diffuse light's slant streams fade faster than an overhead beam.
    assert np.all(col.diffuse[:, -1] / col.diffuse[:, 0] < col.direct(0.0)[:, -1] / col.direct(0.0)[:, 0])


def test_dna_action_spectrum_matches_temis():
    """Setlow's spectrum as TEMIS tabulates it (action_v20.dat), normalised at 300 nm."""
    w = np.array([290.0, 300.0, 310.0, 320.0, 340.0, 400.0])
    temis = np.array([7.93669, 1.0, 3.05838e-2, 9.35366e-4, 4.91196e-5, 3.05343e-5])
    assert np.allclose(O.dna_damage(w), temis, rtol=1e-4)


# ==== The day ======================================================================================================
def test_mean_sunlight_over_the_day():
    """The day's mean of sin(elevation) above the horizon is cos(latitude)/pi at equinox."""
    for lat in (0.0, 45.0, 75.0):
        _, e, f = R.day(lat)
        assert np.mean(np.maximum(np.sin(np.radians(e)), 0)) == pytest.approx(np.cos(np.radians(lat)) / np.pi, rel=1e-4)
        assert f.min() >= 0 and f.max() <= 1
    _, e, _ = R.day(0.0)
    assert np.sum(e > 0) * R.PERIOD_H / R.STEPS == pytest.approx(R.PERIOD_H / 2, abs=1.0)


def test_first_below_is_exact_for_an_exponential():
    z = np.linspace(0, 100, 201)
    assert R.first_below(z, np.exp(-0.05 * z), 0.01) == pytest.approx(np.log(100) / 0.05, rel=1e-9)
    assert np.isnan(R.first_below(z, np.exp(-0.001 * z), 0.01))
    assert R.first_below(z, np.full_like(z, 0.001), 0.01) == 0.0


def test_lit_area_counts_floor_shallower_than_the_limit():
    tables = dict(lat=np.array([10.0, -50.0]), area_km2=np.array([2.0, 3.0]),
                  sea=np.array([np.searchsorted(np.sort([1.0, 5.0, 50.0]), R.DEPTH_EDGES),
                                np.searchsorted(np.sort([2.0, 20.0]), R.DEPTH_EDGES)]))
    out = R.lit_area(tables, 'sea', lambda lat: np.full_like(lat, 10.0))
    assert out['0-15'] == pytest.approx(4.0) and out['45-60'] == pytest.approx(3.0) and out['all'] == pytest.approx(7.0)


# ==== The product ==================================================================================================
def test_product_header_and_provenance():
    p = product()
    assert p['schema'] == R.SCHEMA
    assert set(p) >= {'producer', 'inputs', 'evidence', 'reading_rule', 'units'}
    assert not constants_changed(p['producer']['constants'])
    for path in ('biosphere/living_water/waters.json', 'illumination/calendar/results/transfer.json',
                 'illumination/surface_light/results/surface_light.json'):
        assert p['inputs'][path] == hashlib.sha256((R.ROOT / path).read_bytes()).hexdigest()


def test_checks_against_the_sea_appearance_study_and_the_calendar():
    p = product()
    c = p['checks']
    for name, row in c['sea_appearance_45deg'].items():
        # Their own method reproduces their numbers; this study's surface, streams and fresh water stay within 3%.
        assert row['their_method_reproduced'] == pytest.approx(row['sea_appearance'], rel=0.005)
        assert row['this_study'] == pytest.approx(row['sea_appearance'], rel=0.03)
    assert c['twilight_diffuse_vs_calendar_max_relative'] < 0.02
    for name, eyes in p['eyes']['moon'].items():
        for e in eyes['full_earth'].values():
            assert e['surface_lux'] == pytest.approx(e['calendar_lux'], rel=0.01)


def test_shallow_sea_area_matches_the_aquatic_plan():
    p = product()
    assert p['floor_shallower_than_30m_km2']['sea']['all'] == pytest.approx(233000, rel=0.03)
    assert p['floor_total_km2']['sea']['all'] == pytest.approx(10.6e6, rel=0.01)


def test_zones_are_ordered_and_light_falls_with_depth():
    p = product()
    for world in ('moon', 'earth'):
        z = p['zones'][world]
        for lat in ('0', '45', '75'):
            depths = [z[w][lat]['euphotic_noon_1pct_m'] for w in ('river_mouth', 'productive_coast', 'crater_lake', 'open_sea')]
            assert depths == sorted(depths)
            row = z['open_sea'][lat]
            assert row['euphotic_noon_1pct_dusk_m'] < row['euphotic_noon_1pct_m'] < row['euphotic_noon_1pct_dawn_m']
            assert row['euphotic_noon_1pct_m'] < row['disphotic_end_m']
        for rows in p['par_profiles_umol_m2_s']['by_world'][world].values():
            for r in rows.values():
                assert np.all(np.diff(r['noon']) < 0) and np.all(np.diff(r['day_mean']) < 0)


def test_lit_floor_shrinks_as_the_threshold_rises():
    p = product()
    sea = p['floor']['sea']['open_sea']
    assert sea['seagrass']['all'] < sea['coral_reef']['all'] < sea['microphytobenthos']['all']
    assert sea['charophyte_9pct']['all'] < sea['charophyte_3pct']['all']
