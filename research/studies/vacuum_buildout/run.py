"""A dated enclosed-build capacity path and finite-source sensitivity.

Consumes frozen provisioning responses; no cross-domain implementation import.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def population(year, q):
    u = max(0., min(1., (year-q['ramp_start_year']) /
                          (q['completion_year']-q['ramp_start_year'])))
    return q['seed_population']*math.exp(
        math.log(q['cap_population']/q['seed_population'])*(3*u*u-2*u*u*u))


def interpolate(a, b, fraction):
    if isinstance(a, dict):
        return {k:interpolate(v,b[k],fraction) for k,v in a.items()}
    if isinstance(a, bool) or isinstance(a, str) or a is None:
        return a if fraction < .5 else b
    if isinstance(a, (float,int)):
        return a+(b-a)*fraction
    return a


def response(product, name, growth, architecture='expanded', scenario='base'):
    grid=product['profiles'][name]['architectures'][architecture]['scenarios'][scenario]
    if growth < -1e-12 or growth > grid[-1]['growth_fraction']+1e-12:
        raise ValueError('Construction rate outside frozen engineering response')
    for lo,hi in zip(grid,grid[1:]):
        if growth <= hi['growth_fraction']:
            return interpolate(lo,hi,(growth-lo['growth_fraction'])/
                               (hi['growth_fraction']-lo['growth_fraction']))
    return grid[-1]


def inputs():
    q=json.loads((HERE/'buildout_inputs.json').read_text())
    for path,key in ((HERE/q['cap_source'],'cap_source_sha256'),
                     (HERE/'inputs/resource_inventory.json','inventory_sha256'),
                     (HERE/'inputs/engineering.json','engineering_sha256')):
        if digest(path)!=q[key]: raise ValueError('Unreviewed input: '+str(path))
    engineering=json.loads((HERE/'inputs/engineering.json').read_text())
    if engineering['schema']!='terluna.research.vacuum-buildout-engineering/1':
        raise ValueError('Unexpected engineering schema')
    return q,engineering,json.loads((HERE/'inputs/resource_inventory.json').read_text())


def trajectory(q, engineering, name, scenario='base'):
    rows=[]; lunar=source=0.
    for year in range(q['first_productive_year'],q['completion_year']+2):
        people=population(year,q); scale=people/q['founding_population']
        lead=q['commissioning_years']
        starts=population(year+lead,q)-population(year+lead-1,q)
        g=starts/people
        base=response(engineering,name,g)
        poor=max(0.,min(1.,(year-q['poor_feed_transition_start'])/
                       (q['poor_feed_transition_end']-q['poor_feed_transition_start'])))
        r=base if scenario=='base' else interpolate(
            base,response(engineering,name,g,scenario=scenario),poor)
        physical=r['physical']; ratios=dict(r['capacity_use'])
        jobs=(physical['technical_resident_jobs']+
              q['additional_technical_staff_per_founding_scale'])*scale
        ratios['technical_workforce']=jobs/(q['technical_workforce_fraction']*people)
        lunar+=physical['lunar_raw_kg_year']*scale
        source+=physical['source_feed_kg_year']*scale
        footprint=lunar/(q['regolith_density_kg_m3']*q['worked_depth_m'])/1e6
        rows.append(dict(year=year,population=people,scale=scale,
            construction_people=starts,construction_growth_fraction=g,
            poorer_feed_output_share=0. if scenario=='base' else poor,
            lunar_raw_kg_year=physical['lunar_raw_kg_year']*scale,
            source_feed_kg_year=physical['source_feed_kg_year']*scale,
            cumulative_lunar_raw_kg=lunar,cumulative_source_feed_kg=source,
            excavated_area_km2=footprint,
            prospective_area_km2=footprint/q['accessible_area_fraction'],
            source_process_W=physical['source_process_W']*scale,
            source_heat_W=physical['source_heat_W']*scale,
            lunar_process_W=physical['lunar_process_W']*scale,
            technical_resident_jobs=jobs,
            capacity_use=ratios,max_capacity_use=max(ratios.values()),
            feasible=max(ratios.values())<=1+1e-10))
    return rows


def product():
    q,e,inventory=inputs(); cases=[]
    groups=inventory['inventory']['groups']
    narrow=groups['dark_geometry_excluding_named_science_targets']['mass_at_assumed_density_kg']
    broad=groups['dark_with_diameter']['mass_at_assumed_density_kg']
    for name in e['profiles']:
        if name=='domestic_services': continue
        for scenario in ('base','lunar_stress','source_stress'):
            rows=trajectory(q,e,name,scenario)
            zero=response(e,name,0)['physical'];growth=response(e,name,.01)['physical']
            seed=q['seed_population']/q['founding_population']
            extra=lambda key:max(0.,growth[key]-zero[key])/.01*seed
            sg=e['base_source_inputs']
            seed_lunar=extra('lunar_raw_kg_year')
            seed_source=max(extra('source_water_kg_year')/(sg['water_grade']*sg['water_recovery']),
                extra('source_nitrogen_kg_year')/(sg['nitrogen_grade']*sg['nitrogen_recovery']),
                extra('source_virgin_metal_kg_year')/(sg['metal_grade']*sg['metal_recovery']))
            for row in rows:
                used=row['cumulative_source_feed_kg']
                row['narrow_screen_total_mass_multiple']=used/narrow
                row['broader_dark_catalogue_mass_fraction']=used/broad
                row['required_surveyed_body_mass_kg']=used/q['catalogue_extractable_fraction']
                row['source_feed_with_seed_proxy_kg']=used+seed_source
                row['lunar_feed_with_seed_proxy_kg']=row['cumulative_lunar_raw_kg']+seed_lunar
            cases.append(dict(strategy=name,scenario=scenario,annual=rows,
                seed_modular_feed_proxy=dict(lunar_kg=seed_lunar,source_kg=seed_source,
                    rule='Marginal initial stock from expanded base response; take maximum of separate source-product requirements. Initial survey operations and a separately simulated bootstrap transport campaign remain excluded.'),
                first_capacity_exceedance=next((r['year'] for r in rows if not r['feasible']),None),
                narrow_total_mass_crossing_year=next((r['year'] for r in rows if r['cumulative_source_feed_kg']>narrow),None),
                narrow_ten_percent_crossing_year=next((r['year'] for r in rows if r['cumulative_source_feed_kg']>narrow*q['catalogue_extractable_fraction']),None),
                peak_capacity_use=max(r['max_capacity_use'] for r in rows),
                checkpoints=[r for r in rows if r['year'] in q['checkpoints']]))
    return dict(schema='terluna.research.vacuum-buildout/1',
        evidence=q['evidence'],producer_sha256=digest(__file__),
        dependencies={p:digest(HERE/p) for p in ('buildout_inputs.json','inputs/engineering.json','inputs/resource_inventory.json')},
        reading_rule='Use the dated base path only where capacity_use <= 1. Stress paths retain unmet capacity explicitly; no qualified reserve or funding is implied. Feed totals start at first production and exclude preproduction seed extraction.',
        inventory_comparators=dict(narrow_catalogue_mass_kg=narrow,broader_dark_catalogue_mass_kg=broad,
            assumed_extractable_fraction=q['catalogue_extractable_fraction']),
        schedule=dict(first_productive_year=q['first_productive_year'],completion_year=q['completion_year'],
            cap_population=q['cap_population'],peak_construction_growth=max(r['construction_growth_fraction'] for r in cases[0]['annual']),
            checkpoints=[dict(year=y,population=population(y,q)) for y in q['checkpoints']]),
        cases=cases)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=HERE.parents[1]/'runs/vacuum_buildout/buildout.json')
    args=parser.parse_args(); result=product();args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result['schedule'],indent=2))


if __name__=='__main__': main()
