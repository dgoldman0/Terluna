"""McSnow's vapour diffusivity in the lunar builds: the corrected sources, the old factor and the comparison product."""
import hashlib
import json

import pytest

from atmosphere.snow import mcsnow_aggregation as agg
from atmosphere.snow import mcsnow_diffusivity as dv
from atmosphere.snow import mcsnow_growth_aggregation as grow

RECORD = dv.HERE / 'mcsnow_diffusivity.md'


def _product():
    if not dv.PRODUCT.exists():
        pytest.skip('comparison product absent')
    return json.loads(dv.PRODUCT.read_text())


def test_corrected_lunar_sources_differ_only_in_the_diffusivity_reference():
    if not (grow.SOURCE / 'src/mo_atmo_types.f90').exists():
        pytest.skip('pinned McSnow source copy is not installed')
    differences = dv.source_differences()
    assert differences['earth_sources_unchanged'] and differences['lunar_changes_only_the_diffusivity']
    new = dv.fixed_sources('moon')
    for name in ('mo_atmo.f90', 'mo_cond.f90'):
        formula = [line for line in new[name].splitlines() if 'Dv = 0.211e-4_wp' in line]
        assert formula and all('p_dv /' in line and 'ps /' not in line for line in formula)
    old_types, new_types = dv._original_sources('moon')['mo_atmo_types.f90'], new['mo_atmo_types.f90']
    ps_line = f'    ps    = {agg.lunar_column()["ps"]:.10g}_wp,     &'
    assert ps_line in old_types and ps_line in new_types
    assert '    p_dv  = 101325._wp,' in new_types
    assert differences['changed_lines'] == _product()['build']['changed_lines']


def test_old_factor_is_the_freezing_level_pressure_over_one_atmosphere():
    product = _product()
    assert product['reference_pressure_pa'] == dv.REFERENCE_PA == 101325.
    assert product['old_factor'] == pytest.approx(agg.lunar_column()['ps'] / 101325., rel=1e-12)
    assert f"{product['old_factor']:.3f}" == '0.651'


def test_product_matches_its_producer_and_the_record_table():
    product = _product()
    for path, digest in product['producer']['files'].items():
        assert hashlib.sha256((dv.ROOT / path).read_bytes()).hexdigest().startswith(digest)
    assert set(product['cases']) == set(dv.CASES)
    text = ' '.join(RECORD.read_text().split())
    for name in ('col_q3_f0.01_base', 'col_sat_f0.01_base', 'col_q0.3_f1_base', 'col_deep_q0.3_f1_base'):
        old, new = (product['cases'][name][k]['arriving'] for k in ('old', 'corrected'))
        assert f"{old['water_mm_h']:.2f} → {new['water_mm_h']:.2f} mm/h" in text
        for key in ('d50_mm', 'd90_mm'):
            assert f'{old[key]:.0f} → {new[key]:.0f} mm' in text
    for name in ('col_q3_f0.01_base', 'col_sat_f0.01_base'):
        old, new = (product['cases'][name][k]['onset_depth_km']['10_cm'] for k in ('old', 'corrected'))
        assert f'{old} → {new} km' in text
    melted = product['cases']['col_q3_f0.01_base']['melting']['corrected']
    depths = sorted(melted[site]['depths_ice_gone_m']['99'] / 1e3 for site in ('land', 'sea'))
    assert f'{depths[0]:.2f}–{depths[1]:.2f} km' in text


def test_the_corrected_runs_take_the_corrected_lunar_build_and_keep_the_earth_build():
    """Within corrected() the lunar build and the runs are this module's, the Earth build stays the old cache's (the
    correction leaves its sources alone), and the melting record's raw outputs move beside the corrected runs."""
    old_cache = grow.CACHE
    old_earth = grow.build_path('earth')
    with dv.corrected():
        assert grow.CACHE == dv.CACHE
        assert grow.build_path('moon').parent == dv.CACHE
        assert grow.build_path('earth') == old_earth
    assert grow.CACHE == old_cache and grow.build_path('earth') == old_earth
    assert dv.MELTING_RAW.parent == dv.CACHE
