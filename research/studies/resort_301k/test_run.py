"""Tests for the resort at 301 K across the domains: the product binds its code and inputs, its methane budget
reproduces the waters' methane study at the held pair, its phosphorus room and watering scale follow from their inputs,
and the README's figures are the product's."""
import hashlib
import json
import re

import pytest

from research.studies.resort_301k import run as R

PRODUCT = json.loads(R.OUT.read_text())
README = (R.HERE / 'README.md').read_text()


def test_binds_its_code_and_inputs():
    assert PRODUCT['schema'] == R.SCHEMA
    producer = PRODUCT['producer']
    assert set(producer['files']) == set(R.CODE)
    assert set(producer['inputs']) == set(R.INPUTS)
    for path, digest in {**producer['files'], **producer['inputs']}.items():
        assert hashlib.sha256((R.ROOT / path).read_bytes()).hexdigest() == digest, path


def test_fresh_build_matches():
    fresh = {k: v for k, v in R.build().items() if k != 'producer'}
    stored = {k: v for k, v in PRODUCT.items() if k != 'producer'}
    assert json.loads(json.dumps(fresh)) == stored


def test_methane_reproduces_the_waters_study_at_the_held_pair():
    data = R.load()
    theirs = data[R.METHANE]['air']['lakes_and_seas']['traits']['oligotrophic']['central']['lakes_uncut']
    mine = PRODUCT['methane']['waters']['oligotrophic']['earth']
    assert mine['sources_tg_yr'] == pytest.approx(theirs['sources_tg_yr'], abs=1e-3)
    assert mine['steady_ppm_by_soils'] == pytest.approx(theirs['steady_ppm_by_soils'], abs=1e-2)
    assert mine['held']['removal_tg_yr'] == [0.0, 0.0]
    for trophic, rows in PRODUCT['methane']['waters'].items():
        for community, row in rows.items():
            low, held = row['low']['removal_tg_yr'], row['held']['removal_tg_yr']
            assert all(a >= b for a, b in zip(low, held)), (trophic, community)


def test_phosphorus_room_follows_the_lines():
    p = PRODUCT['phosphorus']
    lines = [p['room'][s]['oligotrophic_line_g_m2_yr'] for s in R.SETTLING]
    assert lines == sorted(lines, reverse=True)
    for s in R.SETTLING:
        room = p['room'][s]
        for land_g, kg in zip(p['land_bioavailable_g_m2_yr'], room['room_kg_yr']):
            assert kg == pytest.approx((room['oligotrophic_line_g_m2_yr'] - land_g) * p['lake_area_m2'] / 1e3)


def test_watering_scale():
    w = PRODUCT['watering']
    per_km2 = w['demand_mm_day'] * 1e-3 * 365.25 * 1e6 / 1e6        # million m3 a year for each km2
    assert per_km2 == pytest.approx(0.71, abs=0.01)
    for row in w['rows']:
        for form in ('additive', 'log_linear'):
            x = row[form]
            assert x['water_km3_yr'] == pytest.approx(x['land_km2'] * per_km2 * 1e6 / 1e9)
            assert x['power_gw'][0] < x['power_gw'][1] and x['albedo_points'][0] < x['albedo_points'][1]
    shares = [r['share_of_deserts'] for r in w['rows']]
    assert shares == sorted(shares)


def figures(pattern):
    return re.findall(pattern, README)


def test_readme_figures():
    m = PRODUCT['methane']['waters']
    poor, meso = m['oligotrophic']['earth'], m['mesotrophic']['earth']
    assert f"{poor['waters_tg_yr']:.2f} Tg a year" in README
    lo, hi = poor['steady_ppm_by_soils']
    assert f'{lo:.2f}–{hi:.2f} ppm' in README
    lo, hi = poor['low']['removal_tg_yr']
    assert f'{lo:.1f}–{hi:.1f} Tg' in README
    lo, hi = (100 * x for x in poor['low']['air_passed_per_year_share'])
    assert f'{lo:.2f}–{hi:.2f}%' in README
    lo, hi = meso['held']['removal_tg_yr']
    assert f'{lo:.1f}–{hi:.1f} Tg' in README
    d = PRODUCT['dimming']
    lo, hi = d['resort_in_all_with_the_first_four_percent']
    assert f'{lo:.1f}–{hi:.1f}%' in README
    w = PRODUCT['watering']['rows']
    one = w[0]
    assert f"{one['additive']['water_km3_yr']:.0f}–{one['log_linear']['water_km3_yr']:.0f} km³" in README
    p = PRODUCT['phosphorus']['room']['aggregates']['per_person_g_yr']
    assert f"{p['8_billion'][1]:.0f}–{p['8_billion'][0]:.0f} g" in README
