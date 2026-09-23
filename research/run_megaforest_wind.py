#!/usr/bin/env python3
"""Reproduce the conditional megaforest static-wind screen.

Run at repository root. Default output is research/runs/megaforest_wind.
An optional --a1-profile consumes an existing shared-profile export read-only.
Nothing in this runner writes under immersion/ or changes its source/assets.
"""
from __future__ import annotations
import argparse
import csv
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import platform
import sys
import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from atmosphere.lower_air import LowerAir, AirConfig, A1Profile
from climate.canopy_flow import CanopyFlow, CanopyConfig
from biosphere.megaforest_mechanics import StaticTree, TreeConfig, RootConfig

BASE_COMMIT = '9ad75cfee4fee5497cf6cce08947dca8cc3ef35f'
HEIGHTS = [100., 200., 300., 400., 500.]
DESIGNS = ['rigid', 'reconfiguring', 'wide_roots', 'thick_stem', 'porous_crown',
           'small_branches', 'combined', 'diameter_capped_8m']
SPEEDS = [0., 5., 10., 20., 30., 40., 60.]


def plain(v):
    if isinstance(v, np.ndarray): return v.tolist()
    if isinstance(v, np.generic): return v.item()
    raise TypeError(type(v).__name__)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, default=plain, allow_nan=False)+'\n')


def write_csv(path, rows):
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=keys, lineterminator='\n')
        writer.writeheader()
        for row in rows:
            writer.writerow({k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in row.items()})


def scenario(height, design, soft_soil=False, water_mm=30., elements=60):
    tree = replace(TreeConfig(), height_m=height, base_diameter_m=height/35,
                   crown_radius_m=.18*height, retained_water_mm=water_mm, elements=elements)
    root = replace(RootConfig(), plate_radius_m=.1*height, plate_depth_m=4*(height/300)**.5)
    if design in ('reconfiguring','combined'): tree=replace(tree,vogel_exponent=-.8)
    if design in ('wide_roots','combined'): root=replace(root,plate_radius_m=1.5*root.plate_radius_m)
    if design == 'thick_stem': tree=replace(tree,base_diameter_m=1.25*tree.base_diameter_m)
    if design == 'combined': tree=replace(tree,base_diameter_m=1.15*tree.base_diameter_m)
    if design == 'porous_crown': tree=replace(tree,crown_frontal_fraction=.15)
    if design == 'small_branches': tree=replace(tree,branch_diameter_fraction=.11)
    if design == 'diameter_capped_8m': tree=replace(tree,base_diameter_m=8.)
    if design not in DESIGNS: raise ValueError('Unknown trait scenario')
    if soft_soil: root=replace(root,cohesion_pa=2000.,pore_pressure_ratio=.8)
    return tree,root


def run(out, a1_profile=None):
    out = Path(out).resolve()
    if out == (ROOT/'immersion').resolve() or (ROOT/'immersion').resolve() in out.parents:
        raise ValueError('Output under immersion/ is outside this work scope')
    air = LowerAir() if a1_profile is None else A1Profile.from_file(a1_profile)
    out.mkdir(parents=True, exist_ok=True)
    envelopes, samples, profiles, configs, flow_checks = [], [], [], [], []
    for height in HEIGHTS:
        flow = CanopyFlow(air, replace(CanopyConfig(),height_m=height))
        emergent = CanopyFlow(air, replace(CanopyConfig(),height_m=.65*height))
        for label,f in [('canopy',flow),('emergent',emergent)]:
            flow_checks.append(dict(tree_height_m=height,exposure=label,**f.diagnostics))
            for z,u in zip(np.linspace(0, height,21),f.ratio(np.linspace(0,height,21))):
                profiles.append(dict(tree_height_m=height,exposure=label,
                    reference_height_m=f.config.height_m,z_m=float(z),wind_over_reference=float(u)))
        for design in DESIGNS:
            for soft in (False,True):
                for water in (0.,30.):
                    tree,root=scenario(height,design,soft,water)
                    for exposure,field in [('uniform',None),('canopy',flow),('emergent',emergent)]:
                        case_id=f'{int(height)}m-{design}-{exposure}-{"soft" if soft else "reference"}-{int(water)}mm'
                        ref_height=height if field is None else field.config.height_m
                        attrs=dict(case_id=case_id,height_m=height,design=design,exposure=exposure,
                            soil='soft' if soft else 'reference',water_mm=water,
                            reference_height_m=ref_height,base_diameter_m=tree.base_diameter_m,
                            root_plate_radius_m=root.plate_radius_m)
                        model=StaticTree(air,tree,root,None if field is None else field.ratio)
                        env=model.envelope()
                        ratio_top=1. if field is None else float(field.ratio(height))
                        env['tree_top_speed_at_first_threshold_m_s']=(None if env.get('first_threshold_m_s') is None
                            else ratio_top*env['first_threshold_m_s'])
                        envelopes.append({**attrs,**env})
                        configs.append(dict(**attrs,tree=asdict(tree),roots=asdict(root)))
                        for speed in SPEEDS:
                            samples.append({**attrs,**model.evaluate(speed)})
    # Same geometry refined independently. This is numerical evidence only.
    refinement=[]
    for height in (100.,300.,500.):
        for design in ('rigid','reconfiguring','diameter_capped_8m'):
            for elements in (30,60,120,240):
                tree,root=scenario(height,design,elements=elements)
                refinement.append(dict(height_m=height,design=design,elements=elements,
                    **StaticTree(air,tree,root).envelope()))
    max_refine=0.
    for row in refinement:
        if row['elements'] != 60: continue
        fine=next(r for r in refinement if r['height_m']==row['height_m'] and r['design']==row['design'] and r['elements']==240)
        for mode in ('root','stem','crown','domain'):
            a,b=row.get(mode+'_threshold_m_s'),fine.get(mode+'_threshold_m_s')
            if a is not None and b is not None and b>0:
                max_refine=max(max_refine,abs(a/b-1))
    references=[e for e in envelopes if e['height_m']==300 and e['exposure']=='uniform' and e['soil']=='reference' and e['water_mm']==30]
    sample_references=[s for s in samples if s['height_m']==300 and s['exposure']=='uniform'
        and s['soil']=='reference' and s['water_mm']==30 and s.get('reference_speed_m_s')==10]
    summary=dict(base_commit=BASE_COMMIT,stage='M0_M1_AND_CONDITIONAL_MEAN_CANOPY_SCREEN',
        envelope_cases=len(envelopes),steady_load_samples=len(samples),canopy_bvp_cases=len(flow_checks),
        refinement_cases=len(refinement),self_weight_buckling_cases=sum(e['status']=='SELF_WEIGHT_BUCKLING' for e in envelopes),
        zero_wind_damage_cases=sum(e.get('first_threshold_m_s')==0 for e in envelopes),
        right_censored_cases=sum(e.get('first_mode')=='right_censored' for e in envelopes),
        first_threshold_in_domain=sum(e.get('first_threshold_within_displacement_domain',False) for e in envelopes),
        first_threshold_outside_domain=sum(e.get('first_threshold_m_s') is not None and not e.get('first_threshold_within_displacement_domain',False) for e in envelopes),
        max_canopy_momentum_relative_residual=max(f['momentum_relative_residual'] for f in flow_checks),
        max_60_vs_240_element_threshold_relative_difference=max_refine,
        air=air.metadata,air_at_0_and_500m=air.sample(np.array([0.,500.])),
        reference_envelopes=references,reference_10m_s_samples=sample_references,
        scientific_status=dict(trait_calibration='hypothetical',root_soil_model='idealized_capacity_screen',
            wind_climate='unmodelled',gust_dynamics='unmodelled',evolution='unmodelled_trait_comparisons_only',
            biological_viability='unestablished',immersion='read_only_optional_input_no_changes'),
        runtime=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
        checks=dict(new_unit_tests='run separately; no status inferred here',
            full_repository_checks='not represented by this runner'))
    write_csv(out/'structural_envelopes.csv',envelopes)
    write_csv(out/'steady_load_samples.csv',samples)
    write_csv(out/'canopy_profiles.csv',profiles)
    write_csv(out/'reference_cases.csv',references)
    write_json(out/'configurations.json',dict(air=air.metadata,cases=configs))
    write_json(out/'grid_refinement.json',refinement)
    write_json(out/'canopy_checks.json',flow_checks)
    write_json(out/'summary.json',summary)
    files=['atmosphere/lower_air.py','climate/canopy_flow.py','biosphere/megaforest_mechanics.py',
           'research/run_megaforest_wind.py','tests/test_megaforest_wind.py']
    hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files}
    write_json(out/'source_hashes.json',hashes)
    print(json.dumps({k:summary[k] for k in ('envelope_cases','steady_load_samples','self_weight_buckling_cases',
        'first_threshold_in_domain','first_threshold_outside_domain','max_canopy_momentum_relative_residual',
        'max_60_vs_240_element_threshold_relative_difference')},indent=2))
    return summary


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=ROOT/'research/runs/megaforest_wind')
    parser.add_argument('--a1-profile',type=Path,help='Existing A1 JSON export; read only, no fallback')
    args=parser.parse_args()
    run(args.out,args.a1_profile)
