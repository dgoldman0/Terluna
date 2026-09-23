#!/usr/bin/env python3
"""Reproduce interacting-patch screens without writing under immersion/.

Default: five budget-matched shapes, directional cases, optional contacts/grafts,
refinement and domain checks, intact-tree maps, and a damage propagation test.
All wind speeds are imposed upstream steady values. Runtime and raw grids go in
research/runs/forest_patch; curated results are a separately identified snapshot.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,platform,sys
from dataclasses import asdict,replace
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from atmosphere.lower_air import LowerAir,A1Profile
from climate.forest_patch_flow import PatchGrid,FlowConfig
from biosphere.forest_patch import make_patch,PatchDrag,tree_mass
from biosphere.forest_patch_mechanics import PatchStructure,InteractionConfig
from biosphere.forest_damage import damage_sequence

SOURCE_FILES=['atmosphere/lower_air.py','biosphere/megaforest_mechanics.py','climate/forest_patch_flow.py','biosphere/forest_patch.py',
    'biosphere/forest_patch_mechanics.py','biosphere/forest_damage.py',
    'research/run_forest_patch.py','tests/test_forest_patch.py']


def plain(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x).__name__)


def write_json(path,value):
    path.write_text(json.dumps(value,indent=2,default=plain,allow_nan=False)+'\n')


def write_csv(path,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='') as stream:
        w=csv.DictWriter(stream,fieldnames=keys,lineterminator='\n');w.writeheader();w.writerows(rows)


def safe_output(path):
    path=Path(path).resolve();blocked=(ROOT/'immersion').resolve()
    if path==blocked or blocked in path.parents:
        raise ValueError('Writing under immersion/ is outside scope')
    return path


def summarize(result):
    rows=result['trees']
    if not rows:return dict(status=result['status'],domain_limited=result['domain_limited'])
    peak=max(rows,key=lambda r:max(r['stem_utilization'],r['root_utilization'],r['crown_utilization']))
    return dict(status=result['status'],domain_limited=result['domain_limited'],
        max_root_utilization=max(r['root_utilization'] for r in rows),
        max_stem_utilization=max(r['stem_utilization'] for r in rows),
        max_crown_utilization=max(r['crown_utilization'] for r in rows),
        max_domain_utilization=max(r['domain_utilization'] for r in rows),
        max_tip_displacement_m=max(r['tip_displacement_m'] for r in rows),
        min_drag_weighted_speed_ratio=min(r['drag_weighted_speed_ratio'] for r in rows),
        max_drag_weighted_speed_ratio=max(r['drag_weighted_speed_ratio'] for r in rows),
        peak_tree=peak['id'],active_contacts=sum(l['kind']=='contact' and l['active'] for l in result['links']),
        failed_links=sum(l['utilization']>=1 for l in result['links']),
        structural_residual=result['linear_equilibrium_relative_residual'])


def run(out,quick=False,a1_profile=None):
    out=safe_output(out);out.mkdir(parents=True,exist_ok=True)
    air=LowerAir() if a1_profile is None else A1Profile.from_file(a1_profile)
    rho=float(air.sample(0)['density_kg_m3']);cfg=FlowConfig()
    if quick:cfg=replace(cfg,nx=48,ny=32,nz=16)
    grid=PatchGrid(cfg)
    cases=[(s,d) for s in ('flat','dome','ramp','irregular') for d in (0.,90.,180.)]+[('gap',0.)]
    modes={'shelter':InteractionConfig(),'contacts':InteractionConfig(crown_contacts=True),
           'contacts_grafts':InteractionConfig(crown_contacts=True,root_grafts=True)}
    summaries=[];tree_rows=[];flows=[];configurations=[];cached={}
    for shape,direction in cases:
        key=f'{shape}-{int(direction)}';trees,meta=make_patch(shape,direction_deg=direction,air=air)
        drag=PatchDrag(grid,trees);flow=grid.solve(drag.coefficient)
        loads,ledger=drag.loads(flow,15.,rho)
        configurations.append(dict(case=key,geometry=meta,trees=[asdict(t) for t in trees]))
        flows.append(dict(case=key,**flow.diagnostics,load_ledger=ledger))
        for mode,interactions in modes.items():
            solved=PatchStructure(air,trees,interactions).solve(loads)
            summaries.append(dict(case=key,shape=shape,direction_deg=direction,interaction=mode,
                upstream_speed_m_s=15.,initial_mass_kg=meta['initial_mass_kg'],leaf_area_m2=meta['leaf_area_m2'],
                max_height_m=meta['max_height_m'],tree_count=len(trees),**summarize(solved)))
            tree_rows.extend(dict(case=key,interaction=mode,**r) for r in solved['trees'])
        if direction==0:
            np.savez_compressed(out/f'{key}_flow.npz',velocity_ratio=flow.centered.astype('float32'),
                pressure_ratio=flow.pressure.astype('float32'),coefficient=drag.coefficient.astype('float32'),
                x_m=grid.centers[0]*cfg.reference_height_m,y_m=grid.centers[1]*cfg.reference_height_m,
                z_m=grid.centers[2]*cfg.reference_height_m)
        cached[key]=(trees,flow,summaries[-3]);print('intact',key,flush=True)
    write_csv(out/'intact_summary.csv',summaries);write_csv(out/'tree_maps.csv',tree_rows)
    write_json(out/'flow_checks.json',flows);write_json(out/'configurations.json',configurations)
    # A one-wave control and a complete re-solved sequence use exactly the same
    # initial patch and imposed 22 m/s. These are not storm probabilities.
    damage={}
    for mode in ('shelter','contacts','contacts_grafts'):
        trees=cached['flat-0'][0]
        result=damage_sequence(air,grid,trees,[22.],modes[mode])
        damage[mode]=result;print('damage',mode,result['status'],result['standing_tree_count'],flush=True)
    control=damage_sequence(air,grid,cached['flat-0'][0],[22.],first_batch_only=True)
    write_json(out/'damage_sequences.json',damage);write_json(out/'first_batch_control.json',control)
    # Keep the same physical geometry when changing grid, enclosure or closure.
    checks=[]
    tests=[('coarse',replace(cfg,nx=48,ny=32,nz=16)),
           ('fine',replace(cfg,nx=96,ny=64,nz=32)),
           ('wide',replace(cfg,ny=cfg.ny*2,length_y_h=cfg.length_y_h*2)),
           ('tall',replace(cfg,nz=cfg.nz*2,height_h=cfg.height_h*2)),
           ('long',replace(cfg,nx=cfg.nx*2,length_x_h=cfg.length_x_h*2)),
           ('mixing_half',replace(cfg,eddy_viscosity_over_uh=.01)),
           ('mixing_double',replace(cfg,eddy_viscosity_over_uh=.04))]
    for shape in ('flat','dome'):
        trees=cached[f'{shape}-0'][0]
        base=cached[f'{shape}-0'][2]
        for label,c in tests:
            g=PatchGrid(c);drag=PatchDrag(g,trees);flow=g.solve(drag.coefficient)
            loads,ledger=drag.loads(flow,15.,rho);r=PatchStructure(air,trees).solve(loads);sr=summarize(r)
            checks.append(dict(shape=shape,variation=label,flow=flow.diagnostics,load_ledger=ledger,
                summary=sr,root_change_relative_to_main=sr['max_root_utilization']/base['max_root_utilization']-1))
            print('sensitivity',shape,label,flush=True)
    # Recompute the full failure cascade at finer spatial resolution.
    fine=PatchGrid(replace(cfg,nx=96,ny=64,nz=32))
    fine_damage=damage_sequence(air,fine,cached['flat-0'][0],[22.],modes['shelter'])
    write_json(out/'fine_damage.json',fine_damage)
    write_json(out/'sensitivity.json',checks)
    import scipy
    summary=dict(base_commit='b5b8cc50914633b6e22587b6d3743108b10ad6f4',
        stage='FINITE_PATCH_AERODYNAMIC_AND_ELASTIC_INTERACTION_SCREEN',
        intact_geometries=len(cases),intact_interaction_cases=len(summaries),trees_in_reference_patch=49,
        map_records=len(tree_rows),sensitivity_cases=len(checks),
        reference_grid=asdict(cfg),air=air.metadata,reference_density_kg_m3=rho,
        worst_mass_divergence=max(f['max_divergence'] for f in flows),
        worst_force_partition_error=max(f['load_ledger']['force_partition_relative_error'] for f in flows),
        max_intact_structural_residual=max(s.get('structural_residual',0) for s in summaries),
        first_batch_failed=49-control['standing_tree_count'],
        damage={k:dict(status=v['status'],failed=49-v['standing_tree_count'],waves=len(v['events']),
                       flow_solves=v['flow_solves']) for k,v in damage.items()},
        fine_damage=dict(status=fine_damage['status'],failed=49-fine_damage['standing_tree_count'],waves=len(fine_damage['events'])),
        max_sensitivity_abs_root_change=max(abs(c['root_change_relative_to_main']) for c in checks),
        scientific_status=dict(wind_climate='UNMODELLED',gusts='UNMODELLED',
            tree_traits='HYPOTHETICAL',root_links='OPTIONAL_HYPOTHETICAL_MECHANICS',
            evolutionary_outcome='UNMODELLED_DESIGN_COMPARISONS_ONLY',
            numerical_grid_limit='SENSITIVITY_TESTS_NOT_CONTINUUM_CERTIFICATION',
            whole_forest_viability='UNESTABLISHED',immersion='UNCHANGED'),
        runtime=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
        source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_FILES})
    write_json(out/'summary.json',summary)
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,default=ROOT/'research/runs/forest_patch')
    p.add_argument('--quick',action='store_true');p.add_argument('--a1-profile',type=Path)
    args=p.parse_args();print(json.dumps(run(args.out,args.quick,args.a1_profile),indent=2,default=plain))
