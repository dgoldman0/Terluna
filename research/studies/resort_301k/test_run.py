"""Tests for the resort at 301 K across the domains: the product binds its code and inputs, its methane budget
reproduces the waters' methane study at the held pair and its brackets, its phosphorus room follows the lake kinds'
lines, its dimming is the heat levers' paired cases, its watering scale follows from its inputs, and the README's
figures are the product's."""
import hashlib
import json
import re

import pytest

from research.studies.resort_301k import run as R

PRODUCT = json.loads(R.OUT.read_text())
README = (R.HERE / 'README.md').read_text()


def pair(v, decimals=1, scale=1.0):
    return f'{scale * v[0]:.{decimals}f}–{scale * v[1]:.{decimals}f}'


def test_binds_its_code_and_inputs():
    assert PRODUCT['schema'] == R.SCHEMA
    producer = PRODUCT['producer']
    assert set(producer['files']) == set(R.CODE)
    assert set(producer['inputs']) == set(R.INPUTS)
    for path, digest in {**producer['files'], **producer['inputs']}.items():
        assert hashlib.sha256((R.ROOT / path).read_bytes()).hexdigest() == digest, path
    assert producer['constants']['MOON_SURFACE_GRAVITY'] == R.MOON_SURFACE_GRAVITY


def test_fresh_build_matches():
    fresh = {k: v for k, v in R.build().items() if k != 'producer'}
    stored = {k: v for k, v in PRODUCT.items() if k != 'producer'}
    assert json.loads(json.dumps(fresh)) == stored


def test_methane_reproduces_the_waters_study():
    data = R.load()
    air = data[R.METHANE]['air']
    theirs = air['lakes_and_seas']['traits']['oligotrophic']['central']['lakes_uncut']
    mine = PRODUCT['methane']['waters']['oligotrophic']['earth']
    assert mine['sources_tg_yr'] == pytest.approx(theirs['sources_tg_yr'], abs=1e-3)
    assert mine['steady_ppm_by_soils'] == pytest.approx(theirs['steady_ppm_by_soils'], abs=1e-2)
    assert mine['held']['removal_tg_yr'] == [0.0, 0.0]
    assert PRODUCT['methane']['share_0_8']['oligotrophic']['waters_tg_yr'] == \
        air['lakes_and_seas_share_tg_yr']['oligotrophic']['high']
    shih = PRODUCT['methane']['shih_bouffard']['oligotrophic']['earth']
    assert shih['waters_tg_yr'] == air['lakes_and_seas_shih_tg_yr']['oligotrophic']['earth']
    assert shih['held']['removal_tg_yr'][1] > 0
    for trophic, rows in PRODUCT['methane']['waters'].items():
        for community, row in rows.items():
            assert all(a >= b for a, b in zip(row['low']['removal_tg_yr'], row['held']['removal_tg_yr']))


def test_phosphorus_room_follows_the_kinds():
    p = PRODUCT['phosphorus']
    assert p['giants_share_of_lake_area'] == pytest.approx(0.70, abs=0.01)
    for s in R.SETTLING:
        v = p['by_settling'][s]
        oxic, anoxic = v['giants_oxic'], v['giants_anoxic']
        assert all(o > a for o, a in zip(oxic['room_kg_yr'], anoxic['room_kg_yr']))
        for land_g, kg in zip(p['land_bioavailable_g_m2_yr'], oxic['room_kg_yr']):
            want = sum(max(v['lines_g_m2_yr'][k] - land_g, 0.0) * p['kinds_area_km2'][k] * 1e3
                       for k in R.KINDS_BY_AREA)
            assert kg == pytest.approx(want)
    lines = [p['by_settling'][s]['giants_oxic']['area_weighted_line_g_m2_yr'] for s in R.SETTLING]
    assert lines == sorted(lines, reverse=True)


def test_dimming_is_the_levers_paired_cases():
    data = R.load()
    d = PRODUCT['dimming']
    assert d['resort_in_all_with_the_first_four_percent'] == \
        data[R.LEVERS]['ranking']['first_four_left_percent']['resort, 301 K']
    lo, hi = d['resort_in_all_with_the_first_four_percent']
    total = d['resort_in_all_percent']
    four = d['levers_first_four_points']
    assert total[0] - four[1] < lo < hi < total[1] - four[0]


def test_watering_scale():
    w = PRODUCT['watering']
    assert w['water_m3_yr_per_km2'] / 1e6 == pytest.approx(0.71, abs=0.01)
    for form, bands in w['bands'].items():
        suns = [b['sunlight_over_mean'] for b in bands.values()]
        assert suns == sorted(suns, reverse=True)
        for b in bands.values():
            assert b['per_1000_km2']['lift_m'] == b['median_depth_m']
    for row in w['rows']:
        for form in ('additive', 'log_linear'):
            x = row[form]
            assert x['water_km3_yr'] == pytest.approx(x['land_km2'] * w['water_m3_yr_per_km2'] / 1e9)
            assert x['power_gw'][0] < x['power_gw'][1] and x['albedo_points'][0] < x['albedo_points'][1]


def test_readme_figures():
    m = PRODUCT['methane']
    poor, meso = m['waters']['oligotrophic']['earth'], m['waters']['mesotrophic']['earth']
    assert f"{poor['waters_tg_yr']:.2f} Tg a year" in README
    assert f"{pair(poor['steady_ppm_by_soils'], 2)} ppm" in README
    low = poor['low']
    assert f"{pair(low['removal_tg_yr'])} Tg; {pair(low['air_passed_per_year_share'], 2, 100)}%" in README
    assert f"{pair(meso['held']['removal_tg_yr'])} Tg" in README
    share = m['share_0_8']['oligotrophic']
    assert f"{share['waters_tg_yr']:.2f} Tg a year" in README
    assert f"{pair(share['steady_ppm_by_soils'], 2)} ppm" in README
    assert f"{pair(share['low']['air_passed_per_year_share'], 2, 100)}%" in README
    shih = m['shih_bouffard']['oligotrophic']['earth']
    assert f"{shih['waters_tg_yr']:.2f} Tg" in README and f"{pair(shih['steady_ppm_by_soils'], 2)} ppm" in README
    assert f"up to {shih['held']['removal_tg_yr'][1]:.2f} Tg" in README
    years = poor['low']['contactors_turnover_yr'][::-1]
    assert f"{years[0]:.0f}–{years[1]:.0f} years" in README
    assert f"{pair(m['low_pair_n2o_decomposition_tw'], 2)} TW" in README
    p = PRODUCT['phosphorus']['by_settling']['aggregates']
    oxic, anoxic = p['giants_oxic']['per_person_g_yr'], p['giants_anoxic']['per_person_g_yr']
    for case in (oxic, anoxic):
        four, eight = case['4_billion'], case['8_billion']
        assert f"{four[1]:.0f}–{four[0]:.0f} g; {eight[1]:.0f}–{eight[0]:.0f} g" in README
    seas = PRODUCT['phosphorus']['seas']
    top = seas['at_1000_years_mg_m3']['rivers_and_dust_top']
    assert f"{top[0]:.0f}–{top[1]:.0f} mg/m³" in README
    assert f"{100 * seas['lakes_keep_share_of_the_land_loss'][0]:.0f}%" in README
    assert f"{seas['deepest_tenth_mg_l']['shih_bouffard_at_the_two_great_basins']['mesotrophic']:.2f}" in README
    assert f"{PRODUCT['methane']['seas_mesotrophic_tg_yr']:.1f} Tg a year" in README
    d = PRODUCT['dimming']
    assert f"{pair(d['resort_in_all_with_the_first_four_percent'])}%" in README
    assert f"{pair(d['first_four_share_of_the_resorts_dimming'], 0, 100)}%" in README
    w = PRODUCT['watering']['rows'][0]
    assert f"{w['additive']['water_km3_yr']:.0f}–{w['log_linear']['water_km3_yr']:.0f} km³" in README
    young = PRODUCT['watering']['deserts']['additive']['young_median_depth_m'] + \
        PRODUCT['watering']['deserts']['log_linear']['young_median_depth_m']
    assert f"{min(young):,.0f}–{max(young):,.0f} m" in README
