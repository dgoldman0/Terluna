"""How aerophytes join the aerial and global ecosystem: the food web, their tenants and phosphorus.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.sky_ecology.run    # ignored research/runs/sky_ecology/

One coupled screen in three components (food_web, tenants, phosphorus) on committed products:
- the coupled aerophyte model and its product (research/studies/aerophytes), whose code is called live
  for the representative bodies and whose product supplies the mass budgets by size;
- the aerial food-web screen, ecology's people screen (the lunar flyers, the land's P loss), its first
  screen (spore settling), the clear-sky surface-light product, the sky-ship air and the sky fleet;
- from other branches, read from Git at pinned commits and bound by SHA-256 until integration: the
  resources branch's phosphorus model and product, the provisioning branch's zonal-wind product and
  note, and its ways-of-living product (the sky towns' shade);
- the 28% atlas grid (an ignored input) for land and sea under the routes, when present.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import subprocess
from pathlib import Path

import numpy as np

from shared.provenance import constants_used
from research.studies.aerophytes import run as aero
from research.studies.aerophytes import structure as aero_structure
from research.studies.sky_ecology import food_web, phosphorus, tenants

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.sky-ecology/1'
OUT = HERE/'results'/'sky_ecology.json'
INPUTS = {
    'research/studies/aerophytes/results/aerophytes.json': 'terluna.research.aerophytes/1',
    'research/studies/aerial_ecology/results/aerial_ecology.json': 'terluna.research.aerial-ecology/1',
    'biosphere/ecology/results/people.json': 'terluna.biosphere.ecology-people/1',
    'biosphere/ecology/results/first_screen.json': 'terluna.biosphere.ecology-first-screen/1',
    'biosphere/ecology/results/taxa.json': 'terluna.biosphere.ecology-taxa/1',
    'illumination/surface_light/results/surface_light.json': 'terluna.illumination.surface-light/1',
    'research/studies/sky_ships/results/sky_ships.json': 'terluna.research.sky-ships/1',
    'research/studies/sky_fleet/results/sky_fleet.json': 'terluna.research.sky-fleet/1',
    'shared/scenarios/water.json': 'terluna.shared.water-scenario/1',
}
SOURCES = 'research/studies/sky_ecology/sources.json'
# Files the imported aerophyte structure component reads when this study calls it.
AEROPHYTE_READS = aero_structure.INPUT_FILES
COMPONENTS = ('food_web', 'tenants', 'phosphorus', 'stoichiometry')
AEROPHYTE_CODE = ('run', 'mechanics', 'envelope', 'biology', 'environment', 'photoperiod', 'navigation', 'trim_cycle',
                  'structure', 'growth', 'water', 'storms', 'sailing', 'gas_biology', 'resources')
# Other branches' files, bound by commit, blob and SHA-256 until the branches are integrated (then read in-tree).
RESOURCES_COMMIT = 'b47c452201a796a86aaf7dbedd9ae348f9ed2d8a'
WINDS_COMMIT = 'd9bc800bc2c9bbb1fa186b376710ed9e81c0919c'
WAYS_COMMIT = 'e364c08e7aa1b9cbe840f30ebc2a5af8132ac2a5'
CROSS_BRANCH = {
    'resources/results/phosphorus.json': dict(branch='domain/resources', commit=RESOURCES_COMMIT,
        blob='5f64d348f6c3a66911e724e3d813bd24b95dfc12',
        sha256='159462a220b6661e457d5adcde1590cf144090e90bb770f3a775d04d709fec8e',
        schema='terluna.resources.phosphorus/1'),
    'resources/phosphorus.py': dict(branch='domain/resources', commit=RESOURCES_COMMIT,
        blob='0b8568bc8af9c5ffc09e6baf9287bc86404cdf09',
        sha256='48f3ef951dc6d781f42de706df7652cbd53d945a67fbd8a6e857c253c9e81e43'),
    'resources/phosphorus.md': dict(branch='domain/resources', commit=RESOURCES_COMMIT,
        blob='030a5a80e4b6edab39c6a908d5bb519e7b31f721',
        sha256='afa8db4d5e8b698c6e4d8bb4b6d5e20f424f5edabd9144ccab7baf4c0a7af138'),
    'climate/results/gcm/zonal_winds_A28_dim5_moon.json': dict(branch='domain/provisioning', commit=WINDS_COMMIT,
        blob='c854c05285193ee493859fa087cb74d05903b5a9',
        sha256='a6a5d0cf96b42a1a18921f61566a638487cec3e602bd2ca752542a1c4a31affe',
        schema='terluna.climate.gcm-zonal-winds/1'),
    'climate/results/gcm/zonal_winds_A28_dim5_moon.npz': dict(branch='domain/provisioning', commit=WINDS_COMMIT,
        blob='20ebd86705a416c8bd18bf1b013dacef7644bc1c',
        sha256='c9b1a0ae8e8733b5d9fb1b60600aa2c43329c9b82bb86ccf4f1bb4a4a58a97e6'),
    'climate/gcm/zonal_winds.md': dict(branch='domain/provisioning', commit=WINDS_COMMIT,
        blob='fe893e2e172e7c22bd711d53444346e2b8db02ae',
        sha256='02bb50a2f2c31eae2172cd118e383f9f329a719576534d5f5412a399368a5e3f'),
    'research/studies/ways_of_living/results/ways_of_living.json': dict(branch='domain/provisioning', commit=WAYS_COMMIT,
        blob='8efea099bf293804c38fcab272a8c038218385ce',
        sha256='cbba4058819147f2d2bb2d02f1dae91320224b3e723e2820bd0f6ba95d64cd0a',
        schema='terluna.research.ways-of-living/1'),
}
# The atlas grid is an ignored input on the research drive; land and sea under the routes need it.
ATLAS = 'geography/products/atlas_28pct_4ppd.npz'
ATLAS_SCHEMA = 'terluna.geography.atlas-grid/1'
TAXA_USED = ('G01', 'G02', 'G03', 'G08', 'G11', 'G14', 'C01', 'C10', 'M18', 'M19', 'M20', 'M21', 'M27', 'M28')
TOWNS_ALOFT_BILLION = 5.          # provisioning's five billion aloft (ways of living, finding 5)


def digest(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_product(path, schema):
    data = json.loads(Path(path).read_text())
    if data.get('schema') != schema:
        raise ValueError(f'Unexpected schema for {path}: {data.get("schema")}')
    return data


def cross_branch_bytes(path):
    """A bound file from another branch: in-tree when its hash matches (after integration), else the Git object."""
    spec = CROSS_BRANCH[path]
    local = ROOT/path
    if local.exists():
        raw = local.read_bytes()
        if hashlib.sha256(raw).hexdigest() == spec['sha256']:
            return raw, 'in-tree'
    try:
        raw = subprocess.run(['git', '-C', str(ROOT), 'show', f"{spec['commit']}:{path}"], check=True,
                             capture_output=True).stdout
    except (OSError, subprocess.CalledProcessError):
        return None, 'unavailable'
    if hashlib.sha256(raw).hexdigest() != spec['sha256']:
        return None, 'hash mismatch'
    return raw, 'git object at the pinned commit'


def cross_branch_json(path):
    raw, source = cross_branch_bytes(path)
    if raw is None:
        return None, source
    data = json.loads(raw)
    schema = CROSS_BRANCH[path].get('schema')
    if schema and data.get('schema') != schema:
        raise ValueError(f'Unexpected schema for {path}')
    return data, source


def load_atlas():
    path = ROOT/ATLAS
    if not path.exists():
        return None, None
    data = np.load(path)
    if json.loads(str(data['metadata']))['schema'] != ATLAS_SCHEMA:
        raise ValueError('Unexpected schema for the atlas grid')
    return data, digest(path)


def wind_profile(winds, top_km=30.):
    """Area-weighted global-mean eastward wind by height from the zonal-wind product, where samples exist."""
    heights = winds['axes']['height_km']
    lat = np.radians(np.asarray(winds['axes']['lat_deg']))
    weight = np.cos(lat)
    out = dict(height_km=[], u_m_s=[])
    for i, h in enumerate(heights):
        if h > top_km:
            break
        row = np.array([np.nan if x is None else x for x in winds['ground_frame']['u_mean_m_s'][i]], float)
        ok = np.isfinite(row)
        if ok.any():
            out['height_km'].append(float(h))
            out['u_m_s'].append(float((row[ok]*weight[ok]).sum()/weight[ok].sum()))
    return out


def plain(x):
    """JSON-ready copy: numpy scalars and arrays to Python, tuples to lists."""
    if isinstance(x, dict):
        return {str(k): plain(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [plain(v) for v in x]
    if isinstance(x, np.ndarray):
        return plain(x.tolist())
    if isinstance(x, np.bool_):
        return bool(x)
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.floating):
        return float(x)
    if isinstance(x, float) and not math.isfinite(x):
        raise ValueError('non-finite value in the product')
    return x


def headlines(fw, tn, ph, towns_shade):
    """The numbers the README quotes, gathered in one place."""
    def lv(body, cover, key):
        return [r[key] for r in fw['levels'] if r['body'] == body and r['cover'] == cover]
    small = fw['reference_small_check']['with_ballonet_partition']
    eb = [r for r in fw['eaten_base'] if r['gas_route'] == 'photolysis']
    eb10 = [r for r in eb if r['diameter_m'] == 10.]
    eb_small = [r['edible_production_c_kg_m2_year'] for r in eb if r['diameter_m'] <= 15. and r['gpp_c_kg_m2_year'] == 2.]
    cap = {r['diameter_m']: r for r in tn['capacity']['round']}
    lone = next(z for z in tn['light']['cases'] if not z['raft'] and z['albedo'] == .1)['zones']
    raft = next(z for z in tn['light']['cases'] if z['raft'] and z['albedo'] == .1)['zones']
    tops = [z['level_top_par_umol_m2_s'] for z in tn['light']['cases'] if not z['raft']]
    chk = ph['resources_check']
    sea = ph['sea_to_land']
    reef = fw['chains']['sky_reef_60m']
    shares = {}
    for c in chk['cases']:
        shares.setdefault(c['state'], []).append(c['throughput_share_of_redfield'])
    return dict(
        cover=dict(central=.01, towns_shade_at_five_billion=towns_shade, range=fw['covers']),
        sky_biomass_mt_c_year_at_1pct=fw['sky']['mt_c_year_at_1pct'],
        reference_light_biomass_kg_c_m2_year={k: dict(growing=v['reference_light_biomass_c_kg_m2_year'],
                                                      mature=v['reference_light_mature_biomass_c_kg_m2_year'])
                                              for k, v in fw['by_class'].items()},
        growth_rate_per_year={k: v['growth_rate_per_year'] for k, v in fw['by_class'].items()},
        level2_production_mt_c_at_1pct=lv('sky_reef_60m', .01, 'level2_production_mt_c'),
        level2_carbon_limit_mt_c_at_1pct=lv('sky_reef_60m', .01, 'level2_production_carbon_limit_mt_c'),
        level2_biomass_mt_c_at_1pct=lv('sky_reef_60m', .01, 'level2_biomass_mt_c'),
        level3_biomass_mt_c_at_1pct=lv('sky_reef_60m', .01, 'level3_biomass_mt_c'),
        reef_p_limit_share=reef['p_corners']['share_range'],
        reef_p_subsidy_g_m2_year=reef['p_corners']['subsidy_range_g_m2_year'],
        ten_metre_living_share=small['living_share_of_reference'],
        ten_metre_gross_lift_kg_m2=small['gross_lift_kg_m2'],
        eaten_base_10m_edible_share=[r['edible_share_of_construction'] for r in eb10],
        eaten_base_10m_turnover_per_year=[r['turnover_per_year'] for r in eb10],
        eaten_base_8_to_15m_edible_production_kg_c_m2_year=[min(eb_small), max(eb_small)],
        filter_energy_return=fw['filters']['energy_return'],
        hundred_metre_steady_tenants_t=cap[100.]['steady_tenants_t'],
        hundred_metre_sudden_upper_bounds_t=[cap[100.]['sudden_trim_t'], cap[100.]['sudden_with_water_t']],
        level_top_par_umol_m2_s=[min(tops), max(tops)],
        underside_light_share=dict(lone=lone['under_cap']['share_of_level_top'], raft=raft['under_cap']['share_of_level_top']),
        resources_published_kg=chk['published']['additional_return_needed_p_kg_yr'],
        resources_reproduced=chk['reproduces'],
        aerophyte_throughput_share_of_redfield={s: [min(v), max(v)] for s, v in shares.items()},
        need_kg_year={f"{c['state']}_{c['retention_case']}_{c['edge']}": c['additional_return_needed_kg']
                      for c in chk['cases']},
        reef_recycled_shares=ph['recycled_shares'],
        sky_need_t_year_at_1pct={f"{r['state']}_{r['retention_case']}_{r['edge']}": r['need_t_year_at_1pct']
                                 for r in ph['sky_need']},
        forest_floor_need_albatross_scale_flyers={
            s: [next(n['flyers'][0] for n in ph['sky_need_flyers'] if n['state'] == s and n['edge'] == 'low'
                     and n['retention_case'] == 'earth_forest' and n['mass_kg'] < 5e3),
                next(n['flyers'][1] for n in ph['sky_need_flyers'] if n['state'] == s and n['edge'] == 'high'
                     and n['retention_case'] == 'earth_forest' and n['mass_kg'] < 5e3)] for s in ('mature', 'growing')},
        earth_seabirds_share_of_land_loss=sea['earth_seabirds_share_of_loss'],
        earth_seabirds_share_of_dissolved_loss=sea['earth_seabirds_share_of_dissolved'],
        land_share_of_routes=ph['land_share_used'],
        falls={f['body']: dict(p_t=f['p_t'], g_m2=f['g_m2_on_ground'], falls_per_year=f['falls_per_year'],
                               falls_per_year_at_growth_rate=f['falls_per_year_at_growth_rate'],
                               years_of_land_loss=f['years_of_land_loss']) for f in ph['falls']})


def results():
    """The whole product as a dict (deterministic; the stored product equals it)."""
    data = {p: load_product(ROOT/p, s) for p, s in INPUTS.items()}
    resources, res_source = cross_branch_json('resources/results/phosphorus.json')
    winds, winds_source = cross_branch_json('climate/results/gcm/zonal_winds_A28_dim5_moon.json')
    npz_raw, npz_source = cross_branch_bytes('climate/results/gcm/zonal_winds_A28_dim5_moon.npz')
    ways, ways_source = cross_branch_json('research/studies/ways_of_living/results/ways_of_living.json')
    sources = {}
    for path in ('resources/phosphorus.py', 'resources/phosphorus.md', 'climate/gcm/zonal_winds.md'):
        sources[path] = cross_branch_bytes(path)[1]
    if resources is None or winds is None or npz_raw is None or ways is None:
        raise RuntimeError('A bound cross-branch input is unreadable: '
                           f'resources {res_source}, winds {winds_source}, npz {npz_source}, ways {ways_source}')
    npz = np.load(io.BytesIO(npz_raw))
    atlas, atlas_hash = load_atlas()
    air_profile = data['research/studies/sky_ships/results/sky_ships.json']['air']
    air = {a['height_km']: a for a in air_profile}
    structural = aero_structure.evaluate(round_spare_pa=aero.SPARE_TRIM_PA, module_spare_pa=aero.SPARE_TRIM_PA_COLONY)
    people = data['biosphere/ecology/results/people.json']
    surface = data['illumination/surface_light/results/surface_light.json']
    aero_product = data['research/studies/aerophytes/results/aerophytes.json']
    water_share = data['shared/scenarios/water.json']['surface_water_share']
    profile = wind_profile(winds)
    fw = food_web.evaluate(air[10], structural, data['research/studies/aerial_ecology/results/aerial_ecology.json'],
                           people, profile, winds, air_profile, data['research/studies/sky_fleet/results/sky_fleet.json'])
    p_surface = surface['columns']['moon_1.2atm']['surface_pressure_pa']
    layer_share = (p_surface-air[10]['pressure_atm']*101325.)/p_surface
    first = data['biosphere/ecology/results/first_screen.json']
    gathering = {}
    g = winds['trajectories']['sequence']['long_run']['gathering']
    heights = winds['trajectories']['heights_km']
    for h in (2.5, 10., 20., 40.):
        k = heights.index(h)
        gathering[str(h)] = dict(days=g['days'], nearest_km_p50=[row[k] for row in g['nearest_km_p50']],
                                 within_cell_share=[row[k] for row in g['within_cell_share']])
    tn = tenants.evaluate(surface, air[10], structural, aero_product, fw['flyers'], fw['chains']['sky_reef_60m'],
                          layer_share, gathering, first['settling'])
    ph = phosphorus.evaluate(fw, resources, winds, npz, atlas, people, fw['flyers'], water_share)
    towns = ways['settings']['roaming_sky_towns']['per_billion']['projected_share_of_moon']*TOWNS_ALOFT_BILLION
    taxa = data['biosphere/ecology/results/taxa.json']
    registered = {e['id'] for key in ('groups', 'modules', 'communities') for e in taxa.get(key, [])}
    own = [HERE/'run.py']+[HERE/f'{c}.py' for c in COMPONENTS]
    aero_files = [ROOT/'research/studies/aerophytes'/f'{c}.py' for c in AEROPHYTE_CODE]
    inputs = list(INPUTS)+[SOURCES]+[p for p in AEROPHYTE_READS if p not in INPUTS]
    product = dict(
        schema=SCHEMA,
        producer=dict(
            study='sky_ecology',
            files={str(p.relative_to(ROOT)): digest(p) for p in own},
            imported_code={str(p.relative_to(ROOT)): digest(p) for p in aero_files},
            inputs={p: digest(ROOT/p) for p in inputs},
            cross_branch_inputs={p: dict(branch=s['branch'], commit=s['commit'], blob=s['blob'], sha256=s['sha256'],
                                         read_from=dict(**{'resources/results/phosphorus.json': res_source,
                                                           'climate/results/gcm/zonal_winds_A28_dim5_moon.json': winds_source,
                                                           'climate/results/gcm/zonal_winds_A28_dim5_moon.npz': npz_source,
                                                           'research/studies/ways_of_living/results/ways_of_living.json': ways_source},
                                                        **sources)[p],
                                         bind='Read in-tree once the branches are integrated; until then from Git at the commit.')
                                 for p, s in CROSS_BRANCH.items()},
            optional_inputs={ATLAS: atlas_hash},
            constants=constants_used(own+aero_files)),
        evidence=('A coupled screen: the coupled aerophyte model called live for representative bodies and its product\'s '
                  'mass budgets; Earth trophic efficiencies, allometries and precedents from the literature (sources.json); '
                  'the clear-sky surface-light product for the light on each surface; the resources branch\'s phosphorus '
                  'budget reproduced and redone; the design run\'s passive routes for where litter lands. No food web, '
                  'population, tenant community, phosphorus cycle or trajectory is simulated, and no organism exists.'),
        reading_rule=('Each component carries its own reading rule. Per-area values are per m² of projected aerophyte '
                      'area (a raft\'s hexagonal cell); sky-wide values take a cover, the share of the Moon\'s 4 pi R² '
                      'surface under aerophytes, 1% unless labelled. [low, high] pairs take every range in a calculation '
                      'at the end that bounds the quantity named; where the ends pull against each other the component '
                      'also gives every corner. Two states of the sky are kept apart: growing (every body spends its '
                      'surplus on new bodies) and mature (bodies hold their size, the dead are replaced). Values computed '
                      'here are screen values; arithmetic on cited values is derived; design guesses are named in each '
                      'component\'s design_guesses or RANGES.'),
        cover_basis=dict(central=.01, central_source='the aerial study\'s central case',
                         towns_shade_share_at_five_billion=towns, towns_source='provisioning ways of living, finding 5',
                         range=list(food_web.COVERS)),
        air_10_km=air[10], layer_share_below_10_km=layer_share, wind_profile=profile,
        food_web=fw, tenants=tn, phosphorus=ph,
        headlines=headlines(fw, tn, ph, towns),
        taxonomy=dict(identifiers_used=list(TAXA_USED), all_registered=all(t in registered for t in TAXA_USED)),
        no_claims=['No organism, tenant, flyer or community named here exists or has been tested.',
                   'Counts of flyers are the numbers a food supply feeds at their field metabolic rate, not populations.',
                   'Phosphorus contents, retention and roosting are design guesses; the budgets are accounting.'])
    return plain(product)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'research/runs/sky_ecology/sky_ecology.json')
    args = parser.parse_args()
    product = results()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(product, indent=1, allow_nan=False)+'\n')
    print(args.output)


if __name__ == '__main__':
    main()
