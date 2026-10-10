"""The taxonomy register hash-binds its producer and the files it reads, lists the files it cites, and its names,
identifiers, ranks and parentage hold together."""
import copy
import hashlib
import json
import math
import re
from pathlib import Path

from biosphere.ecology import taxa
from shared.constants import MOON_RADIUS, SYNODIC_MONTH_DAYS
from shared.provenance import constants_changed

ROOT = Path(__file__).resolve().parents[2]
PRODUCT = json.loads(taxa.OUT.read_text(encoding='utf-8'))
ENTRIES = [*PRODUCT['modules'], *PRODUCT['groups'], *PRODUCT['communities']]
GROUPS = {g['id']: g for g in PRODUCT['groups']}
NODES = {b['name']: b for b in PRODUCT['backbone']}
SOURCES = json.loads(taxa.SOURCES.read_text(encoding='utf-8'))['sources']
DOC = (Path(taxa.HERE) / 'taxonomy.md').read_text(encoding='utf-8')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def ancestors(name):
    chain = []
    while name is not None:
        chain.append(name)
        name = NODES[name]['parent']
    return chain


def test_product_binds_producer_inputs_and_constants():
    assert PRODUCT['schema'] == taxa.SCHEMA
    for path, h in {**PRODUCT['producer']['files'], **PRODUCT['producer']['inputs']}.items():
        assert digest(ROOT / path) == h, f'{path} changed: rerun python -m biosphere.ecology.taxa'
    # only the files the producer reads are hash-bound
    assert sorted(PRODUCT['producer']['inputs']) == taxa.read_inputs() == [
        'biosphere/ecology/results/first_screen.json', 'biosphere/ecology/taxonomy_sources.json']
    assert set(PRODUCT['producer']['constants']) == {'MOON_RADIUS', 'SYNODIC_MONTH_DAYS'}
    assert not constants_changed(PRODUCT['producer']['constants'])


def test_cited_files_are_listed_by_path_and_exist():
    cited = PRODUCT['producer']['cited_files']
    assert isinstance(cited, list) and cited == sorted(set(cited)) == taxa.cited_files()
    assert not set(cited) & set(PRODUCT['producer']['inputs'])
    for path in cited:
        p = Path(path)
        assert not p.is_absolute() and '..' not in p.parts and (ROOT / p).is_file(), path
    entry_paths = {s['path'] for e in ENTRIES for s in e['sources']}
    assert entry_paths <= set(cited) | set(PRODUCT['producer']['inputs'])
    assert set(cited) <= entry_paths                               # nothing listed that no entry cites


def test_product_is_current_with_the_producer():
    assert json.loads(json.dumps(taxa.build(), ensure_ascii=False)) == PRODUCT


def test_register_passes_its_own_check():
    assert taxa.check(PRODUCT) == []


def test_schema():
    for key in ('schema', 'producer', 'evidence', 'reading_rule', 'system', 'backbone', 'modules', 'groups',
                'communities', 'aerophyte_names'):
        assert key in PRODUCT, key
    assert re.fullmatch(r'terluna\.biosphere\.ecology-taxa/\d+', PRODUCT['schema'])
    assert 'design' in PRODUCT['evidence'] and PRODUCT['reading_rule']
    fn = PRODUCT['system']['functional']
    for e in ENTRIES:
        kind = 'group' if e['kind'] in ('group', 'subgroup') else e['kind']
        for field in taxa.REQUIRED[kind]:
            assert field in e, (e['id'], field)
        assert e['status'] in PRODUCT['system']['statuses'] and e['scope'] in PRODUCT['system']['scopes']
        assert e['basis'] in PRODUCT['system']['bases']
        assert isinstance(e['working_name'], str) and e['working_name'] == e['working_name'].strip()
        assert isinstance(e['sources'], list) and isinstance(e['tests'], list)
    for g in PRODUCT['groups']:
        assert g['form'] in fn['forms'] or (g['form'] is None and g['form_candidates']), g['id']
        assert set(g['form_candidates']) <= set(fn['forms'])
        assert g['settings'] and set(g['settings']) <= set(fn['settings']), g['id']
        assert g['guilds'] and set(g['guilds']) <= set(fn['guilds']), g['id']
        assert g['night'] and set(g['night']) <= set(fn['night']), g['id']
        assert g['clock'] and set(g['clock']) <= set(fn['clocks']), g['id']
        assert g['realm'] in {fn['settings'][s]['realm'] for s in g['settings']}, g['id']
        assert g['degree'] in PRODUCT['system']['degrees']
        assert all(r['name'] and r['role'] for r in g['references']), g['id']
    # every status in the vocabulary is used, so the register shows the whole range of verdicts
    assert {e['status'] for e in ENTRIES} == set(PRODUCT['system']['statuses'])
    assert PRODUCT['system']['founding_register']['lineages'] == []


def test_names_and_identifiers_are_unique():
    ids = [e['id'] for e in ENTRIES] + [c['id'] for c in PRODUCT['aerophyte_names']['candidates']]
    assert len(ids) == len(set(ids))
    for e in ENTRIES:
        pattern = {'module': r'M\d\d', 'group': r'G\d\d', 'subgroup': r'G\d\d', 'community': r'C\d\d'}[e['kind']]
        assert re.fullmatch(pattern, e['id']), e['id']
    names = [e['working_name'].lower() for e in ENTRIES]
    assert len(names) == len(set(names))
    assert len(NODES) == len(PRODUCT['backbone'])
    per_candidate = [{f.lower() for part in ('group', 'giants', 'young') if c.get(part)
                      for f in (c[part].get('english'), c[part].get('plural'), c[part].get('scientific')) if f}
                     for c in PRODUCT['aerophyte_names']['candidates']]
    offered = [f for forms in per_candidate for f in forms]
    assert len(offered) == len(set(offered))                       # no name is offered by two candidates
    keys = [s['key'] for s in SOURCES]
    assert len(keys) == len(set(keys))


def test_every_entry_cites_an_existing_project_file():
    for e in ENTRIES:
        assert e['sources'], e['id']
        for s in e['sources']:
            path = Path(s['path'])
            assert not path.is_absolute() and '..' not in path.parts, s['path']
            assert (ROOT / path).is_file(), (e['id'], s['path'])
            assert (s['path'] in PRODUCT['producer']['cited_files']
                    or s['path'] in PRODUCT['producer']['inputs']), (e['id'], s['path'])


def test_backbone_ranks_descend_and_every_node_reaches_life():
    order = PRODUCT['system']['ranks']
    for b in PRODUCT['backbone']:
        chain = ancestors(b['name'])
        assert chain[-1] == 'Life' and len(chain) == len(set(chain)), b['name']
        if b['rank'] is not None:
            for a in chain[1:]:
                r = NODES[a]['rank']
                assert r is None or order.index(r) < order.index(b['rank']), (b['name'], a)


def test_groups_hang_from_the_backbone_with_consistent_degrees_and_parents():
    order = PRODUCT['system']['ranks']
    degrees = PRODUCT['system']['degrees']
    for g in GROUPS.values():
        assert g['anchor'] in NODES, g['id']
        rank = NODES[g['anchor']]['rank']
        floor = degrees[g['degree']]['anchor_at_or_above']
        assert rank is None or order.index(rank) <= order.index(floor), g['id']
        if g['kind'] == 'group':
            assert g['parent'] is None, g['id']
        else:
            parent = GROUPS[g['parent']]
            assert parent['kind'] == 'group' and parent['id'] < g['id'], g['id']
            assert parent['anchor'] in ancestors(g['anchor']), g['id']
            statuses = PRODUCT['system']['status_order']
            assert statuses.index(g['status']) <= statuses.index(parent['status']), g['id']
    # a line stays within one Earth species or a family of them; the aerophytes are a new lineage
    assert GROUPS['G22']['degree'] == 'line' and NODES[GROUPS['G22']['anchor']]['rank'] == 'species'
    assert GROUPS['G01']['degree'] == 'lineage'


def test_modules_and_communities_point_at_real_groups():
    modules = {m['id']: m for m in PRODUCT['modules']}
    used = {m for g in GROUPS.values() for m in g['modules']}
    assert used <= set(modules)
    for m in modules.values():
        assert m['id'] in used or m['status'] == 'set_aside', m['id']
        assert m['id'] not in used or m['status'] != 'set_aside', m['id']
    for c in PRODUCT['communities']:
        assert c['members'] and set(c['members']) <= set(GROUPS), c['id']


def test_carriers_carry_known_elements():
    for g in GROUPS.values():
        assert set(g['carries']) <= set(PRODUCT['system']['functional']['elements'])
        assert bool(g['carries']) <= ('carrier' in g['guilds']), g['id']
    carriers = {g['id'] for g in GROUPS.values() if 'carrier' in g['guilds']}
    assert {'G13', 'G14', 'G42'} <= carriers                       # migrants, colonies, run fish (H1)


def test_check_catches_broken_registers():
    broken = copy.deepcopy(PRODUCT)
    groups = {g['id']: g for g in broken['groups']}
    groups['G18']['anchor'] = 'Fungi'
    groups['G22']['degree'] = 'species'
    groups['G13']['working_name'] = 'aerophytes'
    groups['G03']['status'] = 'carried_forward'
    groups['G01']['sources'].append(dict(path='no/such/file.md', at=''))
    broken['backbone'].append(dict(name='Misplaced', rank='class', parent='Diptera', code=None,
                                   sources=['gbif_backbone']))
    broken['aerophyte_names']['scientific_name'] = 'N1'
    broken['producer']['inputs']['research/decisions.md'] = '0' * 16
    broken['producer']['cited_files'].append('no/such/record.md')
    problems = '\n'.join(taxa.check(broken))
    for expected in ('G18: anchor Fungi lies outside', 'G22: degree species needs an anchor',
                     'working name "aerophytes" used by G01 and G13', 'G03: status carried_forward above',
                     'cited file no/such/file.md does not exist', 'Misplaced (class) sits under Diptera',
                     'the author chooses', 'producer.inputs should hold exactly the files this module reads',
                     'producer.cited_files: no/such/record.md is not a file'):
        assert expected in problems, expected


def test_sources_register_is_complete():
    keys = {s['key'] for s in SOURCES}
    for s in SOURCES:
        for field in ('key', 'citation', 'doi_or_url', 'use', 'access'):
            assert s.get(field), (s.get('key'), field)
    system = PRODUCT['system']
    used = {k for item in (*system['layers'], *system['principles']) for k in item['sources']}
    used |= {k.split(':')[0] for rule in system['naming']['rules'] for k in rule['sources']}
    used |= {k for b in PRODUCT['backbone'] for k in b['sources']}
    used |= {k for e in ENTRIES for r in e.get('references', []) for k in r['lit']}
    used |= {k for c in PRODUCT['aerophyte_names']['candidates'] for k in c['sources']}
    used |= {r['source'] for c in PRODUCT['aerophyte_names']['candidates'] for r in c['roots']}
    assert used <= keys, used - keys


def test_cycle_follows_the_shared_constants():
    cycle = PRODUCT['system']['cycle']
    hours = SYNODIC_MONTH_DAYS * 24
    assert math.isclose(cycle['synodic_hours'], hours, abs_tol=0.01)
    assert math.isclose(cycle['half_cycle_hours'], hours / 2, abs_tol=0.01)
    for lat, v in cycle['dusk_ground_speed_m_s'].items():
        expected = 2 * math.pi * MOON_RADIUS * math.cos(math.radians(float(lat))) / (hours * 3600)
        assert math.isclose(v, expected, abs_tol=1e-3)
    # the ecology register's 15.4 km/h at the equator and 7.7 at 60 degrees; the aerophyte study's 4.30 m/s at 10 km
    assert round(cycle['dusk_ground_speed_km_h']['0'], 1) == 15.4
    assert round(cycle['dusk_ground_speed_km_h']['60'], 1) == 7.7
    assert round(cycle['sun_following_speed_m_s_at_10_km']['0'], 2) == 4.30


def test_size_classes_hold_their_examples():
    screen = json.loads(taxa.SCREEN.read_text(encoding='utf-8'))['settling']['particles']
    classes = PRODUCT['system']['functional']['size_classes']
    placed = set()
    for c in classes['aeroplankton']:
        lo, hi = c['range_um']
        for key, ex in c['examples'].items():
            assert lo <= ex['diameter_um'] < hi, key
            assert ex['residence_days_moon_5_km'] == screen[key]['residence_days']['moon_wet_land_day']
            assert ex['residence_days_earth_1_5_km'] == screen[key]['residence_days']['earth_day']
            placed.add(key)
    assert placed == set(screen)                                   # every screened particle has a class
    heights = [c['range_m'] for c in classes['trees']]
    assert all(a[1] == b[0] for a, b in zip(heights, heights[1:])) and heights[-1][1] is None


def test_aerophyte_scientific_name_waits_for_the_author():
    names = PRODUCT['aerophyte_names']
    assert names['scientific_name'] is None                        # the author chose the common names only
    assert names['common_name'] == GROUPS[names['group']]['working_name'] == 'aerophytes'
    assert names['giant_common_names'] == dict(group='G03', colony_form='sky reefs', round_form=None)
    assert GROUPS['G03']['parent'] == names['group'] and GROUPS['G03']['form'] == 'giant'
    preoccupied = {'holmia', 'pengia', 'peng', 'serenitas', 'nubecula', 'nubila', 'meteora', 'insularia'}
    for c in names['candidates']:
        assert c['group'].get('english') and c['giants'].get('english'), c['id']
        assert c['roots'] and all(r['form'] and r['meaning'] for r in c['roots']), c['id']
        for part in ('group', 'giants'):
            sci = c[part].get('scientific')
            assert sci is None or sci.lower() not in preoccupied, (c['id'], sci)
        assert 'IRMNG' in c['checks']['homonyms'] or 'irmng' in c['sources'], c['id']


def test_taxonomy_md_lists_every_entry_and_candidate():
    text = DOC.lower()
    for e in ENTRIES:
        assert e['id'] in DOC, e['id']
        assert e['working_name'].lower() in text, e['id']
    for c in PRODUCT['aerophyte_names']['candidates']:
        assert c['id'] in DOC and c['group']['plural'].lower() in text and c['giants']['plural'].lower() in text
