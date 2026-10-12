"""Bindings, regeneration and the coupling to the aerophyte product and the other branches' files."""
import hashlib
import json
import math
import subprocess

import pytest

from shared.provenance import constants_changed
from research.studies.sky_ecology import run as m


@pytest.fixture(scope='module')
def product():
    return json.loads(m.OUT.read_text())


def _cross_branch_ready():
    return all(m.cross_branch_bytes(p)[0] is not None for p in m.CROSS_BRANCH)


def test_stored_product_regenerates_exactly(product):
    if not _cross_branch_ready():
        pytest.skip('A bound cross-branch file is not readable from Git here')
    if product['producer']['optional_inputs'][m.ATLAS] and not (m.ROOT/m.ATLAS).exists():
        pytest.skip('The atlas grid (ignored input) is absent; the stored product was made with it')
    assert product['schema'] == m.SCHEMA
    assert product == json.loads(json.dumps(m.results()))


def test_producer_inputs_and_constants_are_current(product):
    p = product['producer']
    for path, h in {**p['files'], **p['imported_code'], **p['inputs']}.items():
        assert m.digest(m.ROOT/path) == h, path
    assert not constants_changed(p['constants'])
    assert set(p['inputs']) >= set(m.INPUTS) | {m.SOURCES}


def test_atlas_hash_when_present(product):
    h = product['producer']['optional_inputs'][m.ATLAS]
    if not (m.ROOT/m.ATLAS).exists():
        pytest.skip('The atlas grid (ignored input) is absent')
    assert h == m.digest(m.ROOT/m.ATLAS)


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


def test_cross_branch_json_schemas():
    for path, spec in m.CROSS_BRANCH.items():
        if 'schema' not in spec:
            continue
        data, source = m.cross_branch_json(path)
        if data is None:
            pytest.skip(f'{path} unavailable ({source})')
        assert data['schema'] == spec['schema']


def test_rejects_wrong_schema(tmp_path):
    p = tmp_path/'wrong.json'
    p.write_text('{"schema": "wrong"}')
    with pytest.raises(ValueError):
        m.load_product(p, m.SCHEMA)


def test_aerophyte_product_is_current_with_the_code_this_study_calls():
    aero = m.load_product(m.ROOT/'research/studies/aerophytes/results/aerophytes.json', 'terluna.research.aerophytes/1')
    for path, h in aero['producer']['files'].items():
        assert m.digest(m.ROOT/path) == h, path


def test_live_aerophyte_cases_reproduce_the_product_rows(product):
    aero = m.load_product(m.ROOT/'research/studies/aerophytes/results/aerophytes.json', 'terluna.research.aerophytes/1')
    rows = aero['size_limits']['10_km']['rows']
    live = {(r['cls'], r['fibre'], r['diameter_m'], r['gpp_c_kg_m2_year']): r for r in product['food_web']['production']
            if r['gas_route'] == 'photolysis'}
    pairs = [(('lesser', 'round', 100., 2.), 'round', 100.), (('round_giant', 'round_strong', 1000., 2.), 'round_strong', 1000.),
             (('sky_reef', 'module', 60., 2.), 'colony_module', 60.), (('round_giant', 'round_strong', 2000., 5.),
                                                                      'round_strong_gpp5', 2000.)]
    for key, name, d in pairs:
        row = next(r for r in rows[name] if r['diameter_m'] == d)
        assert math.isclose(live[key]['surplus_c'], row['surplus_c_kg_m2_year'], rel_tol=1e-9), key
        assert math.isclose(live[key]['structural_dry_kg_m2'], row['structural_dry_kg_m2'], rel_tol=1e-9), key


def test_taxonomy_identifiers_are_registered(product):
    assert product['taxonomy']['all_registered']
    taxa = json.loads((m.ROOT/'biosphere/ecology/results/taxa.json').read_text())
    ids = {e['id'] for k in ('groups', 'modules', 'communities') for e in taxa[k]}
    assert set(m.TAXA_USED) <= ids


def test_cover_basis(product):
    c = product['cover_basis']
    assert c['central'] == .01
    assert math.isclose(c['towns_shade_share_at_five_billion'], .018, rel_tol=1e-9)
    assert min(c['range']) <= .01 <= max(c['range'])


def test_layer_below_ten_km_and_wind_profile(product):
    assert .18 < product['layer_share_below_10_km'] < .2
    prof = product['wind_profile']
    assert prof['height_km'][0] == 0. and prof['u_m_s'][0] < 0 < prof['u_m_s'][prof['height_km'].index(10.)]


def test_headlines_hold_the_cross_checks(product):
    h = product['headlines']
    assert h['resources_reproduced']
    assert math.isclose(h['resources_published_kg'], 10_725_773.75, rel_tol=1e-6)
    assert 0 < h['underside_light_share']['raft'] < h['underside_light_share']['lone'] < 1
    assert 150 < h['hundred_metre_steady_tenants_t'] < 250        # the parent estimate's 150-200 t, checked
    assert .3 < h['ten_metre_living_share'] < .45                  # 'about two-fifths', checked
