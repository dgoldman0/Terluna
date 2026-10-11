"""Conservation, limits and bindings of the food-systems product."""
import json
import math

import pytest

from shared.provenance import constants_changed
from research.studies.food_systems import mixes as mx
from research.studies.food_systems import parameters as par
from research.studies.food_systems import run as m
from research.studies.food_systems import systems as sy
from research.studies.food_systems import zones as zo


@pytest.fixture(scope='module')
def product():
    return json.loads(m.OUT.read_text())


def missing_inputs():
    missing = [p for p in zo.GRID_SCHEMAS if not (m.ROOT / p).exists()]
    for key in m.BOUND:
        try:
            m.bound_bytes(key)
        except zo.MissingInput as error:
            missing.append(str(error))
    return missing


@pytest.fixture(scope='module')
def fresh():
    missing = missing_inputs()
    if missing:
        pytest.skip(f'inputs absent: {missing}')
    return json.loads(m.dumps(m.results()))


def close(a, b, rel=1e-6, abs_=1e-12):
    return math.isclose(a, b, rel_tol=rel, abs_tol=abs_)


# ---- bindings ---------------------------------------------------------------------------------------------------

def test_product_regenerates(product, fresh):
    assert product['schema'] == m.SCHEMA
    assert product == fresh


def test_producer_and_input_hashes(product):
    producer = product['producer']
    for path, h in {**producer['files'], **producer['inputs']}.items():
        assert m.digest(m.ROOT / path) == h, path
    assert producer['sources'] == m.digest(m.HERE / 'sources.json')
    for path, h in producer['ignored_inputs'].items():
        if (m.ROOT / path).exists():
            assert m.digest(m.ROOT / path) == h, path
    assert set(producer['ignored_inputs']) == set(zo.GRID_SCHEMAS)
    assert not constants_changed(producer['constants'])


def test_bound_inputs_pinned(product):
    bound = product['producer']['bound_inputs']
    assert set(bound) == set(m.BOUND)
    for key, b in bound.items():
        pinned = m.BOUND[key]
        assert (b['path_at_integration'], b['branch'], b['commit'], b['sha256'], b['blob']) == \
            (pinned['path'], pinned['branch'], pinned['commit'], pinned['sha256'], pinned['blob'])
        try:
            raw = m.bound_bytes(key)
        except zo.MissingInput:
            pytest.skip(f'{key} is not readable here')
        if pinned['schema']:
            assert json.loads(raw)['schema'] == pinned['schema']


def test_rejects_wrong_schema(tmp_path, monkeypatch):
    bad = tmp_path / 'wrong.json'
    bad.write_text('{"schema": "wrong"}')
    monkeypatch.setattr(m, 'ROOT', tmp_path)
    with pytest.raises(ValueError):
        m.load_json('wrong.json', 'terluna.provisioning.population/1')
    with pytest.raises(zo.MissingInput):
        zo.load_npz(tmp_path, zo.ATLAS)


def test_every_parameter_has_a_source():
    sources = {s['key'] for s in json.loads((m.HERE / 'sources.json').read_text())['sources']}
    used = set()
    for key, entry in par.P.items():
        assert entry['source'], key
        for src in entry['source']:
            if src != 'design_guess':
                assert src in sources, (key, src)
                used.add(src)
    assert sources == used, sources ^ used


def test_ecology_basis_reproduces_the_bound_product(product):
    c = product['checks']['ecology_basis']
    for a, b in zip(c['plant_based_crop_m2_recomputed'], c['plant_based_crop_m2_bound']):
        assert close(a, b, rel=0.01)
    for a, b in zip(c['n2o_kg_n_recomputed'], c['n2o_kg_n_bound']):
        assert close(a, b, rel=0.01)


# ---- physics and conservation -----------------------------------------------------------------------------------

def test_heat_and_dimming_follow_power(product):
    h = product['basis']['heat']
    for name, sysd in product['systems'].items():
        for block in _blocks_with_heat(sysd):
            power, heat = block['power'], block['heat']
            for i in (0, 1):
                arr = heat['array']
                assert close(arr['heat_tw'][i], power[i] * h['heat_per_w']['array'], rel=1e-2)
                assert close(arr['w_m2'][i], arr['heat_tw'][i] * h['moon_w_m2_per_tw'], rel=1e-2)
                assert close(arr['added_dimming_percent'][i], arr['heat_tw'][i] * h['dimming_percent_per_tw'][i],
                             rel=1e-2)


def _blocks_with_heat(node):
    """(power, heat) pairs wherever a block carries both."""
    out = []
    if isinstance(node, dict):
        if 'heat' in node and isinstance(node['heat'], dict) and 'array' in node['heat']:
            p = node.get('power_tw')
            if isinstance(p, dict):
                p = p.get('total', p.get('nitrogen_fixation'))
            if isinstance(p, list):
                out.append(dict(power=p, heat=node['heat']))
        for v in node.values():
            out += _blocks_with_heat(v)
    return out


def test_lift_and_pipe_physics():
    # lifting a cubic metre a kilometre at the Moon's gravity at 100% efficiency
    assert close(sy.lift_kwh_m3(1.0, 1.0), 1000.0 * sy.MOON_SURFACE_GRAVITY * 1000.0 / 3.6e6)
    assert close(sy.lift_kwh_m3(2.0, 0.5), 4.0 * sy.lift_kwh_m3(1.0, 1.0))
    assert sy.pipe_kwh_m3(0.0, 'low') == 0.0
    assert close(sy.pipe_kwh_m3(200.0, 'low'), 2.0 * sy.pipe_kwh_m3(100.0, 'low'))
    assert sy.pipe_kwh_m3(100.0, 'high') > sy.pipe_kwh_m3(100.0, 'low')


def test_ends_are_ordered(product):
    """Every [low, high] cost pair in the systems and mixes runs from the low-cost end to the high-cost end for the
    quantities that are costs."""
    for name, s in product['systems'].items():
        if 'power_tw' in s and isinstance(s['power_tw'], list):
            assert s['power_tw'][0] <= s['power_tw'][1], name
    for mix, info in product['mixes'].items():
        for n, cases in info['cases'].items():
            for kept, c in cases.items():
                assert c['power_tw']['total'][0] <= c['power_tw']['total'][1] + 1e-12, (mix, n, kept)
                assert c['land_km2']['total'][0] <= c['land_km2']['total'][1] + 1e-6, (mix, n, kept)


def test_water_closes_in_glass(product):
    for name in ('fog_desert_glass_vented', 'highland_glass_vented', 'rain_belt_glass_vented',
                 'fog_desert_glass_closed'):
        g = product['systems'][name]
        for i in (0, 1):
            supplied = min(g['rain_on_glass_km3_year'][i], g['consumed_km3_year'][i]) + g['from_seas_km3_year'][i]
            assert close(supplied, g['consumed_km3_year'][i], rel=1e-2, abs_=1e-3), name
    closed = product['systems']['fog_desert_glass_closed']
    assert max(closed['consumed_km3_year']) < min(closed['transpiration_km3_year'])


def test_people_and_land_close_in_every_mix(product):
    dry = product['basis']['zones']['dry_land_km2']
    for mix, info in product['mixes'].items():
        for n, cases in info['cases'].items():
            for kept, c in cases.items():
                ppl = c['people_billion']
                for i in (0, 1):
                    assert close(ppl['aloft'][i] + ppl['surface'][i], ppl['lunar'][i])
                    zones = sum(c['land_km2'][z][i] for z in ('rain_belt', 'middle_wet', 'fog_desert',
                                                               'highlands', 'polar'))
                    assert close(zones, c['land_km2']['total'][i], rel=1e-2, abs_=1.0)
                    assert close(c['dry_land_share_used'][i], c['land_km2']['total'][i] / dry, rel=1e-2)


def test_phosphorus_balance_in_every_mix(product):
    for mix, info in product['mixes'].items():
        for n, cases in info['cases'].items():
            for kept, c in cases.items():
                p = c['phosphorus_mt_year']
                for i in (0, 1):
                    assert close(p['net_new_input'][i], p['lost'][i] + p['glass_discharge'][i] -
                                 p['drawn_from_seas'][i], rel=1e-2, abs_=1e-3)
                    assert p['lost'][i] <= p['eaten'][i]


def test_rain_fed_fill_respects_kept_land(product):
    bands = {tuple(b['abs_lat_deg']): b for b in product['basis']['zones']['bands']}
    wet = sum(bands[b]['dry_land_km2'] for b in ((0.0, 15.0), (15.0, 30.0), (30.0, 45.0)))
    for mix, info in product['mixes'].items():
        for n, cases in info['cases'].items():
            for kept, c in cases.items():
                k = int(kept.split('_')[1]) / 100.0
                if 'fields' in c:
                    for i in (0, 1):
                        assert c['fields']['field_km2'][i] <= wet * (1.0 - k) * (1 + 1e-6)


def _fog_basis(wind, fog_share=1.0, shape=1.2):
    zones = dict(zones=dict(fog_desert=dict(rain_mm_day=0.2, height_above_sea_km_mean=1.0,
                                            km_to_sea_or_great_lake_mean=50.0)))
    return dict(water=dict(pump_efficiency=0.7, air_extraction_l_per_kwh=[2.0, 4.0], earth_aqueducts_kwh_m3={}),
                fog=dict(rings=[dict(name='r', fog=[fog_share] * 36, wind=[wind] * 36, shape_whole_cycle=shape)],
                         lwc_g_m3=[0.1, 0.1], net_efficiency=[0.2, 0.2], stops_below_m_s=1.0, wind_rounding_m_s=0.1,
                         rayleigh_shape=2.0, earth_operational_l_m2_day=[3.0, 10.0],
                         night_land_evaporation_mm_h=[0.0, 0.0]),
                glass=dict(roof_capture=[0.9, 0.9]), dew=dict(earth_l_m2_night={}), zones=zones)


def test_weibull_partial_mean():
    """E[U; U > c] for Weibull winds matches a direct numerical integral, and the Rayleigh case its closed form."""
    for mean, shape, cut in ((0.8, 2.0, 1.0), (0.3, 1.2, 1.0), (2.5, 1.06, 1.0), (1.0, 1.5, 0.0)):
        lam = mean / math.gamma(1.0 + 1.0 / shape)
        steps = 200000
        du = (40.0 * lam - cut) / steps
        total = 0.0
        for i in range(steps):
            u = cut + (i + 0.5) * du
            total += u * (shape / lam) * (u / lam) ** (shape - 1.0) * math.exp(-(u / lam) ** shape) * du
        assert close(sy.weibull_above(mean, shape, cut), total, rel=1e-5, abs_=1e-9)
    sigma = 0.8 / math.sqrt(math.pi / 2.0)
    rayleigh = (1.0 * math.exp(-1.0 / (2.0 * sigma ** 2)) +
                sigma * math.sqrt(2.0 * math.pi) * 0.5 * math.erfc(1.0 / sigma / math.sqrt(2.0)))
    assert close(sy.weibull_above(0.8, 2.0, 1.0), rayleigh, rel=1e-9)


def test_fog_nets_need_wind():
    """Fog nets catch nothing without fog; with fog in still air the bracket's spread still finds a little wind
    above the cut-off; yields grow with the wind, and putting the fog in the windiest share of a bin never catches
    less than pairing it at random."""
    nothing = sy.fog_desert_water(_fog_basis(2.0, fog_share=0.0))['fog_nets']['yields_low_to_high']
    assert nothing['l_per_m2_year_all_wind'] == [0.0, 0.0]
    calm = sy.fog_desert_water(_fog_basis(0.0))['fog_nets']['yields_low_to_high']
    assert calm['l_per_m2_year_wind_above_cutoff'][0] == 0.0
    assert 0.0 < calm['l_per_m2_year_wind_above_cutoff'][1] < calm['l_per_m2_year_all_wind'][1]
    slow = sy.fog_desert_water(_fog_basis(0.5))['fog_nets']['yields_low_to_high']
    fast = sy.fog_desert_water(_fog_basis(2.0))['fog_nets']['yields_low_to_high']
    for key in ('l_per_m2_year_all_wind', 'l_per_m2_year_wind_above_cutoff'):
        assert fast[key][0] > slow[key][0] and fast[key][1] > slow[key][1]
    # fog all the time: every pairing gives the whole wind; at the low end the mean is rounded down by 0.05 m/s
    assert close(fast['l_per_m2_year_all_wind'][0], 0.1 * 1.95 * 0.2 * sy.SECONDS_PER_YEAR / 1e3, rel=1e-6)
    for f in (0.05, 0.3):
        for u in (0.3, 1.0, 2.5):
            rand = sy.bin_fog_flux(f, u, 1.3, 1.0, 'independent')
            windy = sy.bin_fog_flux(f, u, 1.3, 1.0, 'windiest')
            assert windy >= rand - 1e-15
            assert windy <= sy.weibull_above(u, 1.3, 1.0) + 1e-15


def test_vertical_farm_power_scales_with_light_use():
    """The lamp-lit diet's electricity is the diet's photons over the fixtures' efficacy and lighting's share, and
    vertical_farms computes it so."""
    kcal_day, g_mol, kcal_g, eff, share = 2500.0, 0.5, 3.7, 3.0, 0.9
    pair = lambda x: [x, x]
    basis = dict(diet=dict(intake_kcal_day=kcal_day, loss_factor=pair(1.0), variety_factor=pair(1.0),
                           crop_factor={'plant_rich': pair(1.0)}),
                 lamps=dict(staples=dict(potato=dict(edible_g_per_mol=pair(g_mol), kcal_per_g_dry=kcal_g, ppfd=800.0,
                                                     photoperiod_h=12.0)),
                            efficacy_umol_j=pair(eff), lighting_share=pair(share), layers=pair(5.0), floor_kg_m2=500.0),
                 heat=dict(heat_per_w=dict(array=1.0), moon_w_m2_per_tw=0.0264, dimming_percent_per_tw=pair(0.01)),
                 night=dict(hours=354.0, european_final_kwh=650.0), supply_land_km2_per_tw=dict(array=1.0))
    expected_kw = kcal_day / 86400.0 / (g_mol * kcal_g) / (eff * 1e-6) / share / 1e3
    for end in ('low', 'high'):
        vf = sy.vertical_farms(basis, 'plant_rich', end=end)['potato']
        assert close(vf['power_kw_per_resident'], expected_kw, rel=1e-9)
        assert close(vf['power_tw'], expected_kw, rel=1e-9)                 # TW a billion is kW a resident
    doubled = dict(basis, lamps=dict(basis['lamps'], efficacy_umol_j=pair(2.0 * eff)))
    assert close(sy.vertical_farms(doubled, 'plant_rich', end='low')['potato']['power_kw_per_resident'],
                 expected_kw / 2.0, rel=1e-9)


def test_mix_rules_cover_every_mix():
    assert set(mx.MIXES) == {'rain_fed_first', 'greenhouse_heavy', 'fog_desert_glass', 'power_heavy',
                             'fully_lamp_lit', 'power_made', 'sea_supported'}
    assert mx.POPULATIONS == (4.0, 6.0, 8.0)
    for n in mx.POPULATIONS:
        aloft = n * mx.ALOFT_SHARE
        assert 2.0 <= aloft <= 4.0


def test_feasibility_flags_follow_the_zones(product):
    """A case is feasible exactly when every zone it uses, grazing included, fits its land outside the kept share and
    its sea farms fit the seas shallower than 200 m; the feasible range ends where the first zone fills."""
    for mix, info in product['mixes'].items():
        for n, cases in info['cases'].items():
            for kept, c in cases.items():
                for i in (0, 1):
                    shares = [v[i] for v in c['zone_use_share'].values()]
                    assert c['feasible'][i] == all(x <= 1.0 + 1e-3 for x in shares), (mix, n, kept, i)
                t = c['feasible_to']
                if t is None:
                    assert not c['feasible'][0]
                    continue
                assert c['feasible'][0]
                if t == 1.0:
                    assert c['feasible'][1] and c['feasible_range'] is None
                    continue
                assert not c['feasible'][1] and 0.0 <= t < 1.0
                edge = max(v[1] for v in c['feasible_range']['zone_use_share'].values())
                assert 0.99 <= edge <= 1.0 + 1e-3, (mix, n, kept, edge)
                assert c['feasible_range']['feasible'] == [True, True]


def test_field_fill_spills_to_the_fog_desert():
    """When the wet bands' rain and watering cannot grow the need, the rest goes to watered fields in the fog desert,
    each growing what its sunlight allows, and the water follows the growth."""
    bands = {k: dict(dry_land_km2=1e6, highland_dry_land_km2=0.0) for k in ('0-15', '15-30', '30-45', '45-60')}
    zone = dict(dry_land_km2=1e6, height_above_sea_km_mean=1.0, km_to_sea_or_great_lake_mean=50.0)
    basis = dict(rain=dict(support_scaled={k: [0.2, 0.2] for k in bands},
                           watered_yield_cap={k: [0.5, 0.5] for k in bands}),
                 zones=dict(bands_by_key=bands, zones=dict(middle_wet=zone, rain_belt=zone, fog_desert=zone)),
                 water=dict(crop_need_mm_day=[5.0, 5.0], application_efficiency=[0.8, 0.8],
                            pump_efficiency=[0.8, 0.8], irrigation_pressure_kpa=[200.0, 200.0]))
    fill = mx.field_fill(basis, 3e6, 'low', 0.0, {})
    assert close(fill['rain_fed_eq_km2'], 0.6e6)
    assert close(fill['watered_eq_km2'], 0.9e6)
    assert close(fill['fog_desert_watered_km2'], 1.5e6 / 0.5)
    expected = (sy.yearly(5.0, 0.9e6) + sy.yearly(5.0, 1.5e6) * (1.0 - 0.2 / 0.5)) / 0.8
    assert close(fill['withdrawal_km3_year'], expected)

