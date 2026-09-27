"""The equator weather page shows the climate domain's cloud-resolving products unchanged.

    python3 -m pytest visualization/equator-weather/tests
"""
import importlib.util
import json
from pathlib import Path
import re
import numpy as np
import pytest

HERE = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('equator_weather_build', HERE / 'build.py')
build = importlib.util.module_from_spec(spec)
spec.loader.exec_module(build)

if not (build.SUMMARY.exists() and build.FIELDS.exists()):
    pytest.skip('the ring products are missing; run climate/crm/ring_analysis.py', allow_module_level=True)


@pytest.fixture(scope='module')
def built():
    summary, fields, meta = build.load()
    page = build.build()
    data = json.loads(re.search(r'<script type="application/json" id="data">(.*?)</script>', page, re.S).group(1))
    return summary, fields, meta, page, data


def test_every_quoted_number_is_filled(built):
    page = built[3]
    assert '{{' not in page and '__DATA__' not in page


def test_local_time_composites_are_the_products(built):
    summary, _, _, _, data = built
    for surface in ('land', 'water'):
        for key, values in data['by_hour'][surface].items():
            assert values == summary['by_hour_angle'][surface][key]
    assert data['by_hour']['h'] == summary['by_hour_angle']['hour_angle_deg']


def test_sections_and_storm_round_trip(built):
    _, fields, _, _, data = built
    np.testing.assert_allclose(data['sections']['u'], fields['section_eastward_m_s'], atol=0.005)
    np.testing.assert_allclose(np.array(data['sections']['w']) / 100.0, fields['section_upward_m_s'], atol=5e-5)
    np.testing.assert_allclose(data['storm']['rain_mm_h'], fields['storm_rain_mm_h'], atol=0.05)
    assert data['storm']['rain_max_mm_h'] == pytest.approx(float(fields['storm_rain_max_mm_h']), abs=0.05)
    cond = np.maximum(fields['storm_condensate_g_kg'].astype(float), 1e-4)
    np.testing.assert_allclose(10 ** np.array(data['storm']['cond_log10']), cond, rtol=0.012)


def test_rain_map_code_keeps_each_cell_within_its_class():
    rain = np.array([0.0, 0.05, 0.1, 1.0, 12.0, 49.9, 80.0])
    code = build.code_rain(rain)
    assert code[0] == code[1] == 0 and code[-1] == 255
    edges = build.decode_rain(code)
    lo, hi = np.log10(build.RAIN_FLOOR), np.log10(build.RAIN_TOP)
    step = 10 ** ((hi - lo) / 254)
    wet = (rain >= build.RAIN_FLOOR) & (rain < build.RAIN_TOP)
    assert np.all(rain[wet] / edges[wet] <= step ** 1.5) and np.all(rain[wet] / edges[wet] >= step ** -0.5)


def test_text_values_follow_the_summary(built):
    summary = built[0]
    values = build.text_values(summary, 29.53059)
    assert values['width_max'] == f"{summary['storms']['width_km']['max']:,.0f}"
    assert values['rain_all'] == f"{summary['water_budget']['rain_mm_day']['all']:.1f}"
    assert build.crossing([-10.0, 0.0, 10.0], [-1.0, 1.0, 3.0], True) == pytest.approx(-5.0)


def test_page_states_the_products_and_their_hashes(built):
    _, _, meta, page, _ = built
    assert build.sha256(build.SUMMARY)[:16] in page and build.sha256(build.FIELDS)[:16] in page
    assert meta['schema'] in page
