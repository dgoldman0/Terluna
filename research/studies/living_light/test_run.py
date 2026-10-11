"""Bindings, regeneration, the energy and photon budgets, the limits and the skyglow model of the living-light study."""
import hashlib
import json
import math
import subprocess

import numpy as np
import pytest

from shared.provenance import constants_changed
from research.studies.living_light import run as m


@pytest.fixture(scope='module')
def product():
    return json.loads(m.OUT.read_text())


@pytest.fixture(scope='module')
def sky():
    return m.Sky(m.load_product(m.ROOT / 'illumination/calendar/results/transfer.json',
                                m.INPUTS['illumination/calendar/results/transfer.json']))


def _cross_branch_ready():
    return all(m.cross_branch_bytes(p)[0] is not None for p in m.CROSS_BRANCH)


# ---- Bindings and regeneration --------------------------------------------------------------------------------------
def test_stored_product_regenerates_exactly(product):
    if not _cross_branch_ready():
        pytest.skip('A bound cross-branch file is not readable from Git here')
    assert product['schema'] == m.SCHEMA
    assert product == json.loads(json.dumps(m.results()))


def test_producer_inputs_and_constants_are_current(product):
    p = product['producer']
    for path, h in {**p['files'], **p['imported_code'], **p['inputs']}.items():
        assert m.digest(m.ROOT / path) == h, path
    assert not constants_changed(p['constants'])
    assert set(p['inputs']) >= set(m.INPUTS) | {m.SOURCES}


@pytest.mark.parametrize('path', sorted(m.CROSS_BRANCH))
def test_cross_branch_inputs_bound_by_commit_blob_and_hash(path, product):
    raw, source = m.cross_branch_bytes(path)
    if raw is None:
        pytest.skip(f'{path} unavailable ({source})')
    spec = m.CROSS_BRANCH[path]
    assert hashlib.sha256(raw).hexdigest() == spec['sha256']
    assert source in ('in-tree', 'git object at the pinned commit')
    try:
        blob = subprocess.run(['git', '-C', str(m.ROOT), 'rev-parse', f"{spec['commit']}:{path}"], check=True,
                              capture_output=True, text=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        pytest.skip('git rev-parse unavailable')
    assert blob == spec['blob']
    recorded = product['producer']['cross_branch_inputs'][path]
    assert recorded['commit'] == spec['commit'] and recorded['sha256'] == spec['sha256']


def test_rejects_wrong_schema(tmp_path):
    p = tmp_path / 'wrong.json'
    p.write_text('{"schema": "wrong"}')
    with pytest.raises(ValueError):
        m.load_product(p, m.SCHEMA)


def test_natural_light_matches_the_daily_life_calendars(sky):
    """The calendar's evaluator called here reproduces the provisioning branch's 22 calendars to their print."""
    daily, _ = m.cross_branch_json('research/studies/daily_life/results/daily_life.json')
    if daily is None:
        pytest.skip('daily_life.json not readable')
    lv = daily['levels']
    for name, c in daily['light']['calendars'].items():
        r = m.night_light(sky, c['lat'], c['lon'], (lv['dusk'], lv['dark'], lv['scotopic']))
        assert r['darkest_lux'] == pytest.approx(c['darkest_lux'], rel=2e-3, abs=2e-4), name
        assert r['midnight_lux'] == pytest.approx(c['midnight_lux'], rel=2e-3, abs=2e-4), name
        hb = r['hours_below']
        assert hb[f"{lv['dusk']:g}"] == pytest.approx(c['below_dusk_h'], abs=0.1), name
        assert hb[f"{lv['dark']:g}"] == pytest.approx(c['below_dark_h'], abs=0.1), name
        assert hb[f"{lv['scotopic']:g}"] == pytest.approx(c['below_scotopic_h'], abs=0.1), name


def test_natural_light_matches_mains_places(sky):
    places, _ = m.cross_branch_json('research/studies/joint_synthesis/results/places.json')
    if places is None:
        pytest.skip('places.json not readable')
    for name, p in places['places'].items():
        r = m.night_light(sky, p['lat'], p['lon'], (0.1,))
        assert r['darkest_lux'] == pytest.approx(p['night']['natural_darkest_lux'], rel=5e-3, abs=2e-4), name
        assert r['midnight_lux'] == pytest.approx(p['night']['natural_midnight_lux'], rel=5e-3, abs=2e-4), name


# ---- Photometry, photons and energy ---------------------------------------------------------------------------------
def test_cie_tables_peak_where_they_should():
    assert m.LAMBDA_NM[np.argmax(m.V_PHOTOPIC)] == 555.0 and m.V_PHOTOPIC.max() == 1.0
    assert m.LAMBDA_NM[np.argmax(m.V_SCOTOPIC)] in (505.0, 510.0)
    narrow = m.per_photon(555.0, 0.5)
    assert narrow['lm_s_photopic'] / narrow['photon_j'] == pytest.approx(m.KM_PHOTOPIC, rel=1e-6)


def test_photons_for_one_lux_at_555_nm():
    """Independent of the tables' shape: at 555 nm one lux is 1/683 W/m2, about 4.09e15 photons a second per m2."""
    photon_j = m.PLANCK * m.SPEED_OF_LIGHT / 555e-9
    assert m.photons_for_lux(1.0, 555.0, 0.5) == pytest.approx(1.0 / m.KM_PHOTOPIC / photon_j, rel=1e-4)
    assert 4.0e15 < 1.0 / m.KM_PHOTOPIC / photon_j < 4.2e15


def test_photon_budget_closes():
    """Photons times their energy is the optical power; a closed canopy gives pi times its luminance."""
    for peak, fwhm in ((475.0, 40.0), (530.0, 80.0), (560.0, 70.0)):
        q = m.photons_for_lux(1.0, peak, fwhm)
        lum = m.lambert_luminance(q, peak, fwhm)['photopic']
        assert math.pi * lum * m.GEOMETRIES['closed canopy overhead']['view_factor'] == pytest.approx(1.0, rel=1e-9)
        p = m.per_photon(peak, fwhm)
        assert q * p['photon_j'] * m.KM_PHOTOPIC >= 1.0          # efficacy can only fall short of 683 lm/W


def test_pathway_costs_stand_above_the_thermodynamic_floor():
    for pathway in m.PATHWAYS:
        for i in (0, 1):
            per_photon = m.turnover_cost(pathway, i) / m.atp_values(i)['per_mol_C']
            assert per_photon > m.floor_mol_c_per_mol_photons(m.per_photon(530.0, 80.0)['photon_j'])


def test_hedge_view_factor_limits():
    assert m.hedge_view_factor(1.0, 1e9) == pytest.approx(1.0, abs=1e-6)
    assert m.hedge_view_factor(1.0, 1e-9) == pytest.approx(0.0, abs=1e-6)
    assert m.hedge_view_factor(1.0, 1.0) == pytest.approx(1 - 1 / math.sqrt(2))


# ---- Skyglow ----------------------------------------------------------------------------------------------------------
def test_skyglow_falls_with_walkers_slope_on_earth(product):
    """Earth's single-scattering curve falls about as distance to the -2.5 between 10 and 50 km (Walker's law)."""
    s = product['skyglow']['model']['earth_rayleigh_slope_10_50_km']
    assert -3.0 < s < -2.3


def test_skyglow_optical_depths_match_the_solved_column(product):
    model = product['skyglow']['model']
    assert model['moon_tau_550'] == pytest.approx(0.757, rel=0.03)
    assert model['earth_tau_550'] == pytest.approx(0.0974, rel=0.01)


def test_town_calibration_returns_one_against_itself():
    """A town read through Earth's curve in both worlds keeps its Earth area: the calibration closes."""
    d = m.SKY_DISTANCES_M
    curve = np.exp(-np.log(d) * 2.5)
    towns = m.Towns(dict(rayleigh=dict(earth=curve, moon=curve)), 1000.0)
    for people in m.TOWN_PEOPLE:
        assert towns.factor('rayleigh', people, 1.0, 1.0) == pytest.approx(1.0, rel=1e-6)
        assert towns.factor('rayleigh', people, 4.0, 1.0) < 1.0      # a brighter natural sky shrinks the zone
        assert towns.factor('rayleigh', people, 1.0, 0.25) < 1.0     # less light shrinks it too
        assert towns.factor('rayleigh', people, 1.0, 0.0) == 0.0


def test_moon_class_areas_partition_the_moon():
    c = m.moon_class_areas()
    assert sum(b['area_share'] for b in c.values()) == pytest.approx(1.0, abs=1e-9)
    for b in c.values():
        assert sum(b['classes'].values()) == pytest.approx(1.0, abs=1e-6)


def test_band_slices_lie_in_their_bands():
    for lo, hi in m.BANDS.values():
        qs = m.band_latitudes(lo, hi)
        assert len(qs) == m.BAND_SLICES and lo < qs[0] and qs == sorted(qs) and qs[-1] < hi


def test_class_weights_sum_to_one_at_each_latitude():
    totals = {}
    for band, lat, cls, e, lon, w in m.grid_points():
        totals[(round(lat, 6), cls)] = totals.get((round(lat, 6), cls), 0.0) + w
    assert all(abs(v - 1.0) < 1e-12 for v in totals.values())


def test_glow_cost_is_linear_in_light_and_time():
    a = m.glow_cost(1e15, 0.41, 'firefly', 100.0, 1)['g_c_m2']
    assert m.glow_cost(2e15, 0.41, 'firefly', 100.0, 1)['g_c_m2'] == pytest.approx(2 * a)
    assert m.glow_cost(1e15, 0.41, 'firefly', 200.0, 1)['g_c_m2'] == pytest.approx(2 * a)
    assert m.glow_cost(1e15, 0.205, 'firefly', 100.0, 1)['g_c_m2'] == pytest.approx(2 * a)


def test_waymark_contrast():
    assert m.marker_contrast(m.marker_needed_cd_m2(0.01, 1.0), 0.01) == pytest.approx(1.0)
    assert m.marker_contrast(1e-6, 1.0) < 0


# ---- The product's limits --------------------------------------------------------------------------------------------
def test_hours_stay_within_the_night(product):
    for g in product['natural']['grid']:
        hb = g['hours_below']
        values = [hb[f'{x:g}'] for x in m.LEVELS_LUX]
        assert all(0.0 <= v <= m.NIGHT_H + 1.0 for v in values)
        assert values == sorted(values)                      # more hours below a higher level


def test_areas_and_shares_hold_their_limits(product):
    for practice, designs in product['skyglow']['designs'].items():
        for design, d in designs.items():
            for home, (lo, hi) in d['per_person_m2_by_home'].items():
                assert 0.0 <= lo <= hi, (practice, design, home)
            for mix, rows in d['shares_of_dry_land'].items():
                for b, r in rows.items():
                    for key in ('all_on_the_land', 'half_aloft'):
                        lo, hi = r[key]
                        assert 0.0 <= lo <= hi <= 1.0
        for home in designs['shielded']['per_person_m2_by_home']:
            sh = designs['shielded']['per_person_m2_by_home'][home]
            tw = designs['twilight']['per_person_m2_by_home'][home]
            lv = designs['living']['per_person_m2_by_home'][home]
            assert lv[0] <= tw[0] + 1e-9 <= sh[0] + 2e-9 and lv[1] <= tw[1] + 1e-9 <= sh[1] + 2e-9


def test_emitter_photon_budgets_close(product):
    for r in product['emitters']['organisms']:
        p = m.per_photon(r['peak_nm'], r['fwhm_nm'])
        assert r['optical_w_m2'] == pytest.approx(r['photons_m2_s'] * p['photon_j'], rel=2e-3)
        assert r['path_lux']['closed canopy overhead'] == pytest.approx(math.pi * r['luminance_cd_m2']['photopic'],
                                                                        rel=2e-3)


def test_targets_cost_more_than_the_floor(product):
    for need, row in product['emitters']['targets'].items():
        for span in ('night', 'summit_hours'):
            floor = row['thermodynamic floor'][span]['g_c_m2'][0]
            for name, r in row.items():
                assert r[span]['g_c_m2'][0] >= floor * (1 - 2e-3), (need, name, span)


def test_every_cited_source_is_recorded():
    recorded = {s['id'] for s in json.loads((m.ROOT / m.SOURCES).read_text())['sources']}
    cited = {e['source'] for e in m.EMITTERS} | set(m.FIREFLY_FLASH['sources'])
    cited |= {x for f in m.FLASHES.values() for x in f['sources']}
    assert cited <= recorded, cited - recorded


def test_taxa_identifiers_are_registered(product):
    assert product['taxonomy']['all_registered']


def test_mesopic_luminance_limits():
    """CIE 191's MES2: photopic at 5 cd/m2 and above, scotopic at 0.005 and below, a light of S/P 1 unchanged."""
    assert m.mesopic_luminance(10.0, 25.0) == pytest.approx(10.0, rel=1e-6)
    assert m.mesopic_luminance(0.001, 0.004) == pytest.approx(0.004, rel=1e-6)
    for lp in (0.01, 0.1, 1.0):
        assert m.mesopic_luminance(lp, lp) == pytest.approx(lp, rel=1e-6)


def test_natural_floor_lies_within_earths_sky():
    lo, hi = m.natural_floor(0.764)
    assert 0.0 < lo < (1 - m.AIRGLOW_SHARE) * m.NATURAL_SKY_CD_M2 < hi == m.NATURAL_SKY_CD_M2


def test_never_below_means_the_darkest_moment_meets_the_level(product):
    """A reading counts as never below a need exactly when its darkest moment reaches the level (levels tabulated)."""
    for g in product['natural']['grid']:
        for need, spec in m.NEEDS.items():
            hours = g['hours_below'][f"{spec['lux']:g}"]
            if g['darkest_lux'] >= spec['lux'] * 1.001:
                assert hours == 0.0, (need, g['lat'], g['earth_class'])
            elif g['darkest_lux'] <= spec['lux'] * 0.999:
                assert hours > 0.0, (need, g['lat'], g['earth_class'])


def test_lit_shares_recompute(product):
    """The spread mix's lit share at 6 billion, all on the land, rebuilt from the per-person areas."""
    ways, _ = m.cross_branch_json('research/studies/ways_of_living/results/ways_of_living.json')
    if ways is None:
        pytest.skip('ways_of_living.json not readable')
    d = product['skyglow']['designs']['lunar_design']['shielded']
    shares, metro = m.mix_shares(ways, 'spread')
    dry = m.load_product(m.ROOT / 'biosphere/ecology/results/people.json',
                         m.INPUTS['biosphere/ecology/results/people.json'])['land']['dry_land_km2']
    for end in (0, 1):
        lit = metro * 1e9 * d['metropolis_m2_per_person'][end] / 1e6
        for h, w in shares.items():
            land = ways['settings'][h]['home_dry_land_km2']
            lit += min(land, (6.0 - metro) * w * 1e9 * d['per_person_m2_by_home'][h][end] / 1e6)
        assert lit / dry == pytest.approx(d['shares_of_dry_land']['spread']['6']['all_on_the_land'][end], rel=2e-3)


def test_delivered_share_scales_the_cost_and_spans_the_range(product):
    """Glow costs go as one over the share of the light that reaches the path, and the two shares bound the range."""
    lo_d, hi_d = m.DESIGN_GUESSES['delivered_share']
    for need, row in product['emitters']['targets'].items():
        for name, r in row.items():
            for span in ('night', 'summit_hours'):
                by = r[span].get('by_delivered_share')
                if by is None:
                    assert name == 'thermodynamic floor'
                    continue
                a, b = by[f'{hi_d:g}']['g_c_m2'], by[f'{lo_d:g}']['g_c_m2']
                assert b[0] == pytest.approx(a[0] * hi_d / lo_d, rel=2e-3) and b[1] == pytest.approx(a[1] * hi_d / lo_d,
                                                                                                    rel=2e-3)
                assert r[span]['g_c_m2'] == pytest.approx([a[0], b[1]], rel=2e-3)


def test_dark_belt_share_recomputes(product):
    ways, _ = m.cross_branch_json('research/studies/ways_of_living/results/ways_of_living.json')
    if ways is None:
        pytest.skip('ways_of_living.json not readable')
    life = product['night_life']
    land = {h: ways['settings'][h]['home_dry_land_km2'] for h in m.HOMES}
    belt = sum(land[h] * life['homes'][h]['share_with_any_hour_below_0_01_lux'] for h in m.HOMES)
    req = life['dark_night_requirement']
    assert req['dark_belt_km2'] == pytest.approx(belt, rel=2e-3)
    assert req['dark_belt_dry_land_share'] == pytest.approx(belt / sum(land.values()), rel=2e-3)
    assert req['dark_belt_darkest_lux'][1] < m.MELATONIN_LUX[0]


def test_natural_sky_ordering(product):
    """The far side's low latitudes have the darkest natural sky; twilight's high latitudes never let towns pollute."""
    by = product['skyglow']['natural_sky_by_band_and_class']
    assert by['equatorial']['never']['natural_zenith_over_earth'][0] == min(
        c['natural_zenith_over_earth'][0] for b in by.values() for c in b.values())
    assert all(not c['latitudes_where_towns_light_their_sky'] for c in by['high'].values())
    belt = product['skyglow']['ground_lux_at_the_polluted_boundary_in_the_dark_belt']
    edge = product['skyglow']['ground_lux_at_the_polluted_boundary']
    assert edge[0] <= belt[0] <= belt[1] <= edge[1]
