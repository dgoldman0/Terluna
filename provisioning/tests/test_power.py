"""The power screen: its binding and regeneration, the ITU-R models against their own Earth values, and the sums the
note states."""
import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np

from provisioning import power as P
from shared.provenance import constants_changed

ROOT = Path(__file__).resolve().parents[2]
PRODUCT = json.loads(P.OUT.read_text())
NOTE = (P.HERE / 'power.md').read_text()


def close(a, b, rel=2e-3):
    return math.isclose(a, b, rel_tol=rel, abs_tol=1e-12)


def test_product_binds_its_producer_inputs_and_constants():
    assert PRODUCT['schema'] == P.SCHEMA
    files, inputs = PRODUCT['producer']['files'], PRODUCT['producer']['inputs']
    assert set(files) == set(P.CODE) and set(inputs) == set(P.INPUTS)
    for path, h in {**files, **inputs}.items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == h, path
    assert not constants_changed(PRODUCT['producer']['constants'])
    for key in ('evidence', 'reading_rule', 'assumptions'):
        assert PRODUCT[key]


def test_a_fresh_run_matches_the_product():
    assert json.loads(P.dumps(P.results())) == PRODUCT


def test_the_absorption_model_reproduces_the_itu_earth_values():
    e = P.ITU_STANDARD['rho_g_m3'] * P.ITU_STANDARD['t_k'] / 216.7
    for f, (go, gw) in P.ITU_SPECIFIC.items():
        mo, mw = P.specific_attenuation(f, P.ITU_STANDARD['p_hpa'], e, P.ITU_STANDARD['t_k'])
        assert close(float(mo), go, 1e-9) and close(float(mw), gw, 1e-9), f
    edges = P.layer_edges(0.0, 100.0, imax=P.ITU_SLANT['layers'])
    ao, aw = P.slant_path(P.ITU_SLANT['f_ghz'], P.ITU_SLANT['elevation_deg'], P.earth_profile, P.EARTH_RADIUS / 1e3,
                          edges)
    assert close(ao + aw, P.ITU_SLANT['attenuation_db'], 1e-9)
    checks = PRODUCT['absorption']['earth_checks']
    assert checks['worst_relative_error'] < 1e-9 and abs(checks['slant']['relative_error']) < 1e-9
    # The 922 layers end where the validation file's do: the last layer's bottom at 99.457 km, 0.99966 km thick.
    assert close(edges[-2], 99.45702171642462, 1e-9) and close(edges[-1] - edges[-2], 0.99965968594378762, 1e-9)


def test_the_rain_model_reproduces_p838():
    for f, ref in P.P838_TABLE5.items():
        c = P.p838(f)
        for got, want in zip((c['kh'], c['ah'], c['kv'], c['av']), ref):
            assert close(got, want, 6e-4), f                    # the table prints four or five figures
    for f, ref in P.P838_VALIDATION.items():
        c = P.p838(f)
        for got, want in zip((c['kh'], c['ah'], c['kv'], c['av']), ref):
            assert close(got, want, 1e-9), f
    # Equations 4-5 at a 45 degree tilt: the mean k, and alpha weighted by k.
    c = P.p838(2.45)
    assert close(c['k'], (c['kh'] + c['kv']) / 2, 1e-12)
    assert close(c['alpha'], (c['kh'] * c['ah'] + c['kv'] * c['av']) / (c['kh'] + c['kv']), 1e-12)


def test_the_column_is_the_atmosphere_products():
    for ts, row in PRODUCT['absorption']['column_check'].items():
        for key, value in row['product'].items():
            assert close(row['here'][key], value, 1e-3), (ts, key)
    col = PRODUCT['absorption']['column']
    assert col['scenario'] == 'moon_1.2atm' and col['surface_temperature_k'] == 295.0
    comp = PRODUCT['absorption']['composition']
    assert close(comp['o2_strength_scale'], comp['o2_mole_fraction'] / P.EARTH_O2, 2e-3)
    assert close(comp['o2_mole_fraction'] * P.COLUMN['dry_pressure_pa'], 21227.0, 2e-3)    # Earth's partial pressure


def test_the_air_takes_more_lower_down_and_at_lower_elevations():
    order = ('90_deg', '60_deg', '45_deg', '30_deg', '20_deg', '10_deg', '5_deg')
    for f in ('2.45_ghz', '5.8_ghz'):
        rows = PRODUCT['absorption']['by_elevation'][f]
        for h, row in rows.items():
            db = [row[k]['db'] for k in order]
            assert db == sorted(db), (f, h)
        for k in order:
            db = [rows[h][k]['db'] for h in ('summit', 'towns', 'highland', 'lowland')]
            assert db == sorted(db) and len(set(db)) == 4, k
    a = PRODUCT['absorption']
    assert 6.0 < a['moon_over_earth_zenith']['2.45_ghz'] < 9.0
    # Oxygen at Earth's partial pressure takes less than at Earth's mole fraction, by its share of the strength.
    zen = a['by_elevation']['2.45_ghz']['lowland']['90_deg']
    wet = zen['db'] - zen['dry_db']
    earths = a['zenith_2_45_with_earths_oxygen_fraction_db'] - wet
    assert close(earths, zen['dry_db'] / a['composition']['o2_strength_scale'], 5e-3)


def test_the_link_sums():
    link = PRODUCT['absorption']['link']
    c = link['chain']
    assert close(c['tx'] * c['capture'] * c['rx'] * c['bus'], link['dc_efficiency_without_air'], 1e-3)
    assert close(1 / (c['capture'] * c['rx'] * c['bus']), link['heat_per_w_delivered_without_air'], 1e-3)
    best = link['with_the_best_tiles']
    for e in (0, 1):
        a = link['share_of_beams_lost'][e]
        assert close(link['heat_per_w_delivered_with_air'][e], link['heat_per_w_delivered_without_air'] / (1 - a))
        assert close(link['absorbed_in_the_air_tw'][e], link['intake_tw'][e] * a / (1 - a), 3e-3)
        assert close(link['dc_efficiency_with_air'][e], link['dc_efficiency_without_air'] * (1 - a))
        heat = link['array_heat_with_close_computing_tw']
        assert close(heat['with_air'][e], heat['synthesis'][e] + link['absorbed_in_the_air_tw'][e], 3e-3)
        assert close(link['rise_in_the_arrays_heat'][e], link['absorbed_in_the_air_tw'][e] / heat['synthesis'][e], 3e-3)
        assert best['share_of_beams_lost'][e] < a
    ledgers = json.loads((ROOT / P.LEDGERS).read_text())['array_power']
    assert np.allclose(link['intake_tw'], np.add(ledgers['sky_towns']['intake_tw'], ledgers['ground']['intake_tw']),
                       rtol=2e-3)


def test_the_ring_sums():
    r = PRODUCT['ring']
    for e in (0, 1):
        assert close(r['rf_installed_tw'][e], P.BEAMS_PER_RECEIVER * P.BEAM_RATING_SHARE * r['rf_in_use_tw'][e])
        assert close(r['electricity_at_the_tiles_tw'][e],
                     r['rf_in_use_tw'][e] / PRODUCT['absorption']['link']['chain']['tx'])
        assert close(r['transmitting_tiles'][e], r['rf_installed_tw'][e] * 1e3 / r['tile']['whole_tile_rf_gw'][e]
                     / P.TILE_UTILISATION)
        assert close(r['share_of_fleet'][e], r['transmitting_tiles'][e] / r['stack']['fleet_tiles'])
        assert close(r['collector_and_aperture_km2_installed'][e],
                     r['transmitting_tiles'][e] * r['tile']['clear_side_km'] ** 2)
        receivers = r['beams']['towns'] + r['beams']['ground_receivers'][e]
        assert close(r['beams']['total'][e], P.BEAMS_PER_RECEIVER * receivers)
    assert close(r['beams']['patch_gw'], 230.0 * (P.RECEIVER_SIDE_KM * 1e3) ** 2 / 1e9)
    assert close(r['tile']['sunlight_gw'], P.SOLAR_CONSTANT * r['tile']['clear_side_km'] ** 2 / 1e3)
    # The middle ring alone hides within asin(R/r) of its plane's poles; the stack keeps tiles over every horizon.
    assert close(r['middle_ring_alone_hidden_within_deg_of_its_poles'], math.degrees(math.asin(P.MOON_RADIUS / 2e7)))
    pole = r['by_latitude']['87.5']
    assert pole['tiles_above_horizon']['least'] > 5e6 and pole['share_of_time_with_no_usable_tile'] == 1.0
    for lat, v in r['by_latitude'].items():
        if float(lat) < 60.0:
            assert v['share_of_time_with_no_usable_tile'] == 0.0, lat
    # A tile overhead stays in view for the arc within acos(R/r) of the zenith, at the ring's speed over the surface.
    n = math.sqrt(P.MOON_GM / 2e7 ** 3) + P.OMEGA_SUN_FRAME
    assert close(r['overhead_pass_middle_ring']['horizon_h'], 2 * math.acos(P.MOON_RADIUS / 2e7) / n / 3600)
    # About the share of the fleet within acos(R/r) of a point's zenith stands above the equator's horizon.
    eq = r['by_latitude']['2.5']
    assert 0.40 < eq['tiles_above_horizon']['mean'] / r['stack']['fleet_tiles'] < 0.52
    assert eq['highest_tile_deg']['mean'] > 85.0


def test_the_ground_loads_and_receivers_add_up():
    g = PRODUCT['grid']
    comp = g['composition_billion']
    assert close(comp['metropolis'] + sum(comp['dry_highland'].values()) + comp['scattered'], 4.0, 1e-3)
    assert close(comp['scattered_0_30'] + comp['scattered_high_latitude'], comp['scattered'], 1e-3)
    ledgers = json.loads((ROOT / P.LEDGERS).read_text())['array_power']['ground']
    for e in (0, 1):
        assert close(sum(v['intake_tw'][e] for v in g['groups'].values()), ledgers['intake_tw'][e], 3e-3)
        km2 = sum(v['receiver_km2'][e] for v in g['groups'].values())
        assert close(km2, ledgers['receiver_km2_at_1978_peak'][e], 5e-3)
        assert close(g['receivers_over_settled_land'][e], g['receiver_km2'][e] / g['settled_km2'][e])
    x = g['exposure']
    for e in (0, 1):
        assert 0.0 < x['edge_buffer_km'][e] < 1.0
        assert close(x['closed_km2_per_patch'][e], (P.RECEIVER_SIDE_KM + 2 * x['edge_buffer_km'][e]) ** 2)
        assert close(x['closed_over_settled_land'][e], x['closed_land_km2'][e] / g['settled_km2'][e])
    # The buffer: where a Gaussian edge holding 95% within the whole tile's spot falls from 230 to 10 W/m2.
    spot = PRODUCT['ring']['apertures']['2.45_ghz']['whole_tile_spot_m'][0]
    w = spot / 2 / math.sqrt(-math.log(0.05) / 2)
    assert close(230.0 * 0.5 * math.erfc(math.sqrt(2) * x['edge_buffer_km'][0] * 1e3 / w), 10.0, 1e-2)
    for lat, line in g['poleward_lines'].items():
        assert close(line['km'], (float(lat) - (g['always_usable_to_deg'] - 2.5)) * P.KM_PER_DEG)
        assert close(line['line_loss'][0], 0.02 * line['km'] / 1000 + 0.02, 3e-3)
        assert close(line['line_loss'][1], 0.035 * line['km'] / 1000 + 0.02, 3e-3)


def test_the_reserves_add_up():
    f = PRODUCT['reserves']['ring_fault']
    hours = f['days'] * 24.0
    for e in (0, 1):
        assert close(f['critical_total_tw'][e], sum(v[e] for v in f['critical_tw'].values()), 3e-3)
        assert close(f['energy_twh'][e], f['critical_total_tw'][e] * hours)
        assert close(f['towns']['critical_tw'][e] + f['ground']['critical_tw'][e], f['critical_total_tw'][e], 3e-3)
        assert close(f['ground']['hydrogen_mt'][e], f['ground']['energy_twh'][e] / (P.FUEL_CELL * P.H2_LHV_KWH_KG))
        per_twh = 1.0 / (P.FUEL_CELL * P.H2_LHV_KWH_KG)                   # Mt of hydrogen a TWh delivered
        assert close(f['towns']['hydrogen_mt'][e], f['towns']['critical_tw'][e] * hours * per_twh)
        assert close(f['ground']['electricity_to_make_it_twh'][e], f['ground']['energy_twh'][e] / P.H2_ROUND_TRIP)
        assert close(f['ground']['pumped_hydro_km3_at_korolevs_head'][e], f['ground']['energy_twh'][e] / 1.35)
    t = PRODUCT['reserves']['single_tiles']
    lam, tau = t['tile_faults_per_year'], t['reassign_s'] / (P.HOURS_PER_YEAR * 3600)
    assert close(t['double_losses_per_receiver_year'], 3 * lam * 2 * lam * tau)
    assert close(PRODUCT['reserves']['storms']['rating_headroom_db'], 10 * math.log10(1.5))
    for who in ('towns', 'ground'):
        for row in PRODUCT['reserves']['reefs']['beam_columns'][who].values():
            for case in row['by_cover'].values():
                assert all(0.0 <= x < 1.0 for x in case['share_of_time_blocked_without_switching'])


def test_the_summit_sums():
    s = PRODUCT['summit']
    for e in (0, 1):
        d = s['delivered_tw']
        assert close(d['total'][e], d['residents'][e] + d['close_computing'][e] + d['food'][e], 3e-3)
        heat_per_w = PRODUCT['absorption']['link']['heat_per_w_delivered_without_air']
        assert close(s['intake_tw'][e], d['total'][e] * heat_per_w)
        assert close(s['receiver_km2'][e], s['intake_tw'][e] * 1e12 / 230.0 / 1e6)
        assert close(s['store_twh'][e], s['critical_tw'][e] * PRODUCT['reserves']['ring_fault']['days'] * 24.0)
        assert close(s['korolev_km3'][e], s['store_twh'][e] / 1.35)
    assert s['delivered_tw']['residents'][0] * 1e3 < s['brief_gw'] < s['delivered_tw']['residents'][1] * 1e3


def percent(x, digits=1):
    return f'{100 * x:.{digits}f}%'


def span(a, b, digits=1):
    return f'{100 * a:.{digits}f}–{100 * b:.{digits}f}%'


def test_the_note_states_the_products_figures():
    a, r, link = PRODUCT['absorption'], PRODUCT['ring'], PRODUCT['absorption']['link']
    zen = a['by_elevation']['2.45_ghz']['lowland']
    stated = [
        percent(zen['90_deg']['share_lost']), percent(zen['30_deg']['share_lost']),
        f"{zen['90_deg']['db']:.3f} dB",
        percent(a['by_elevation']['2.45_ghz']['towns']['90_deg']['share_lost']),
        span(link['with_the_best_tiles']['share_of_beams_lost'][0], link['share_of_beams_lost'][0]),
        f"{link['with_the_best_tiles']['heat_per_w_delivered'][0]:.2f}–"
        f"{link['heat_per_w_delivered_with_air'][0]:.2f} W",
        f"{a['moon_over_earth_zenith']['2.45_ghz']:.1f} times",
        f"{r['electricity_at_the_tiles_tw'][0]:.0f}–{r['electricity_at_the_tiles_tw'][1]:.0f} TW",
        span(*PRODUCT['grid']['receivers_over_settled_land'], digits=0),
        f"{PRODUCT['reserves']['ring_fault']['critical_total_tw'][0]:.1f}–"
        f"{PRODUCT['reserves']['ring_fault']['critical_total_tw'][1]:.1f} TW",
        f"{PRODUCT['summit']['delivered_tw']['total'][0]:.2f}–{PRODUCT['summit']['delivered_tw']['total'][1]:.2f} TW"]
    for text in stated:
        assert text in NOTE, text


def test_every_source_entry_is_complete_and_cited():
    sources = json.loads(P.SOURCES.read_text())['sources']
    for s in sources:
        assert set(s) == {'citation', 'doi_or_url', 'use', 'access'} and all(s.values()), s['citation']
    text = NOTE + Path(P.__file__).read_text()
    for s in sources:
        words = re.split(r'[ ,&(]+', s['citation'])
        if words[0] == 'ITU-R':
            key = re.search(r'P\.\d{3}-\d+', s['citation'] + ' ' + s['use']).group(0)
        elif words[0] == 'US':
            key = ' '.join(words[:2])
        elif words[0] == "Hydrogen's":
            key = 'heating value'
        else:
            key = words[0]
        assert key in text, s['citation']
