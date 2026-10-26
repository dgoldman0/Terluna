"""Depth-control identity and reuse of the frozen snow population calculation."""
import json

import pytest

from atmosphere.snow import flakes as sf, growth_depth as gd, population as pop
from atmosphere.snow.habits import GEOMETRIES
from shared.provenance import constants_changed
from shared.baselines import recorded_file


def test_depth_cases_isolate_cloud_top_without_changing_water_supply_or_geometry():
    rows = gd.cases()
    assert {(r['top_c'], r['snowfall_set_mm_h']) for r in rows} == {(-5., 1.), (-5., 3.), (-10., 1.), (-10., 3.)}
    assert len({r['id'] for r in rows}) == 4
    assert not ({r['id'] for r in rows} & {r['id'] for r in pop.cases()})
    assert all(r['geometry'] == gd.GEOMETRY and r['seed'] == 0 and r['crystals_per_litre'] == 100. for r in rows)
    moon, _ = sf.moon_column()
    depth = lambda top: moon.height_of(sf.T0 + top) - moon.z_freeze
    assert depth(-5.) / depth(-15.) == pytest.approx(1/3)
    assert depth(-10.) / depth(-15.) == pytest.approx(2/3)


def test_cache_invalidates_for_top_temperature_and_constants():
    spec = gd.cases()[0]
    args = dict(files={'runner.py': 'abc'}, inputs={'column.json': 'def'}, constants={'GAS_CONSTANT': 8.3})
    key = gd.cache_identity(spec, **args)
    assert gd.cache_identity(dict(spec, top_c=-10.), **args) != key
    assert gd.cache_identity(spec, **dict(args, constants={'GAS_CONSTANT': 8.4})) != key
    assert gd.cache_identity(spec, **dict(args, files={'runner.py': 'ghi'})) != key


def test_reused_population_run_conserves_mass_for_a_short_control():
    spec = dict(gd.cases()[0], id='earth_short_control', world='earth', top_c=-2.,
                bins=4, inject=4, fill=1., gather=.5)
    result = pop.run_case(spec, sf.earth_column(), GEOMETRIES[gd.GEOMETRY], {'inputs': {}})
    assert result['top_c'] == -2.
    assert abs(result['mass_ledger_kg_m2']['relative_residual']) < 1e-10
    assert result['freezing_level']['totals']['mass_flux_kg_m2_s'] * 3600 == pytest.approx(
        result['summary']['snowfall_mm_h'])
    assert result['layers'][-1]['temperature_k'] < result['layers'][0]['temperature_k']


def test_committed_depth_product_keeps_full_provenance_and_case_identity():
    if not gd.OUT.exists():
        pytest.skip('Growth-depth product not generated yet')
    product = json.loads(gd.OUT.read_text())
    assert product['schema'] == gd.SCHEMA
    assert not constants_changed(product['producer']['constants'])
    assert {r['id'] for r in product['runs']} == {r['id'] for r in gd.cases()}
    for path, expected in product['producer']['files'].items():
        assert pop.sha(recorded_file(path, expected)) == expected, path
    for path, expected in product['producer']['inputs'].items():
        assert pop.sha(recorded_file(path, expected)) == expected, path
    for row in product['runs']:
        assert row['cloud_depth_m'] > 0
        assert abs(row['mass_ledger_kg_m2']['relative_residual']) < 1e-10
