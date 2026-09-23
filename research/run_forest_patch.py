#!/usr/bin/env python3
"""Reproduce conditional interacting-forest screens; no writes under immersion/."""
from __future__ import annotations
import argparse,csv,json,sys,hashlib,platform
from pathlib import Path
from dataclasses import asdict,replace
import numpy as np
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from atmosphere.lower_air import LowerAir
from climate.forest_flow import Grid,SpatialFlow
from biosphere.forest_patch import make_patch,TreeDrag,biomass_proxy
from biosphere.forest_mechanics import StandMechanics,Support,share_soil
BASE='b5b8cc50914633b6e22587b6d3743108b10ad6f4'


def plain(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x).__name__)
def js(path,data):path.write_text(json.dumps(data,indent=2,default=plain,allow_nan=False)+'\n')
def csvout(path,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,keys,lineterminator='\n');w.writeheader()
        for r in rows:w.writerow({k:json.dumps(v,default=plain) if isinstance(v,(list,dict)) else v for k,v in r.items()})
def output_guard(path):
    path=Path(path).resolve();forbidden=(ROOT/'immersion').resolve()
    if path==forbidden or forbidden in path.parents:raise ValueError('immersion/ is outside scope')
    return path


def metrics(response):
    rows=response['trees']
    if not rows:return dict(status=response['status'],standing_trees=0)
    root=max(rows,key=lambda r:r['root_utilization']);stem=max(rows,key=lambda r:r['stem_utilization'])
    crown=max(rows,key=lambda r:r['crown_utilization']);domain=max(rows,key=lambda r:r['domain_utilization'])
    link=max((l['utilization'] for l in response['links']),default=0.)
    return dict(status=response['status'],standing_trees=len(rows),maximum_root_utilization=root['root_utilization'],
        maximum_stem_utilization=stem['stem_utilization'],maximum_crown_utilization=crown['crown_utilization'],
        maximum_link_utilization=link,maximum_domain_utilization=domain['domain_utilization'],
        maximum_vertical_drag_to_weight=max(r['vertical_drag_to_weight'] for r in rows),
        damage_utilization=max(root['root_utilization'],stem['stem_utilization'],crown['crown_utilization'],link),
        all_in_domain=all(r['domain_valid'] for r in rows),most_loaded_root=root['key'],
        active_crown_contacts=response['active_crown_contacts'],
        mean_root_utilization=float(np.mean([r['root_utilization'] for r in rows])),
        equilibrium_relative_residual=response['equilibrium_relative_residual'],
        horizontal_force_balance_relative_residual=response['horizontal_force_balance_relative_residual'])


def envelope(model,drag,F10,max_wind=60.):
    """First crossing on a scanned grid; later roots are intact-stand counterfactuals.

    Only allowed for speed-independent nondimensional drag. This prevents an
    invalid U^2 scaling when crown reconfiguration is enabled.
    """
    if not drag.rigid:raise ValueError('Envelope scaling requires fixed drag coefficients')
    cache={}
    def at(u):
        if u not in cache:cache[u]=metrics(model.evaluate(drag,F10*(u/10.)**2))
        return cache[u]
    if at(0).get('status')!='COUPLED_SMALL_DISPLACEMENT_STAND':return {'status':at(0)['status']}
    def crossing(key):
        old=0.
        if at(0)[key]>=1:return 0.
        for u in np.linspace(2,max_wind,30):
            if at(float(u))[key]>=1:return float(brentq(lambda v:at(v)[key]-1,old,float(u),xtol=1e-5))
            old=float(u)
        return None
    damage=crossing('damage_utilization');domain=crossing('maximum_domain_utilization')
    valid=damage is not None and (domain is None or damage<=domain) and at(damage)['maximum_vertical_drag_to_weight']<=.05
    return dict(status='CONDITIONAL_INTACT_STAND_ENVELOPE',first_nominal_damage_m_s=damage,
        displacement_domain_m_s=domain,damage_before_domain=valid,
        damage_state=None if damage is None else at(damage),speed_search_limit_m_s=max_wind)


def damage_sequence(forest,air,grid,support,winds,direction=0.,max_events=80):
    """Quasi-static prescribed removal / distal-crown-loss sequence.

    Links break first. Root/stem loss removes standing geometry; no falling-body
    impact or fallen-log aerodynamic resistance is included. Crown loss halves
    planform area once. Domain violations stop the sequence without inventing a
    survival outcome. Soil ownership stays at the initial allocation.
    """
    if any(f.tree.vogel_exponent!=0 for f in forest):raise ValueError('Damage ramp presently requires rigid drag')
    winds=np.asarray(winds,float)
    if winds.ndim!=1 or not len(winds) or np.any(~np.isfinite(winds)) or np.any(winds<0) or np.any(np.diff(winds)<=0):raise ValueError('Increasing nonnegative wind ramp required')
    state=list(forest);shares,soil=share_soil(state,support.soil_cell_m);broken=set();history=[];events=[];flows=[];stopped=None
    field=None;dirty=True;F10=None
    for speed in winds:
        round_id=0
        while True:
            if not any(f.alive for f in state):stopped='ALL_STANDING_TREES_REMOVED';break
            if dirty:
                drag=TreeDrag(state,grid,air);field,F10,diag=SpatialFlow(grid,direction).solve(drag,10.,initial=field)
                flows.append(diag);model=StandMechanics(state,air,support,shares,broken);dirty=False
            response=model.evaluate(drag,F10*(speed/10.)**2);m=metrics(response)
            history.append(dict(wind_m_s=float(speed),round=round_id,**m))
            if m['status']!='COUPLED_SMALL_DISPLACEMENT_STAND':stopped=m['status'];break
            if not m['all_in_domain']:stopped='DOMAIN_CENSORED';break
            failed_links=[l['key'] for l in response['links'] if l['utilization']>=1]
            if failed_links:
                broken.update(failed_links);events.append(dict(wind_m_s=float(speed),round=round_id,kind='link_loss',links=failed_links))
                model=StandMechanics(state,air,support,shares,broken)
            else:
                changed=[]
                for r in response['trees']:
                    i=r['index'];f=state[i]
                    if max(r['root_utilization'],r['stem_utilization'])>=1:
                        kind='root_loss' if r['root_utilization']>=r['stem_utilization'] else 'stem_loss'
                        state[i]=replace(f,alive=False);changed.append(dict(key=f.key,kind=kind))
                    elif r['crown_utilization']>=1:
                        if f.crown_remaining<1:stopped='CROWN_DAMAGE_RULE_EXHAUSTED';break
                        state[i]=replace(f,crown_remaining=.5);changed.append(dict(key=f.key,kind='prescribed_distal_crown_loss'))
                if changed:
                    events.append(dict(wind_m_s=float(speed),round=round_id,kind='geometry_loss',trees=changed));dirty=True
                elif stopped is None:break
            if stopped is not None:break
            round_id+=1
            if len(events)>=max_events:stopped='EVENT_CAP_UNRESOLVED';break
        if stopped is not None:break
    return dict(status=stopped or 'RAMP_COMPLETED',events=events,history=history,flow_solves=flows,initial_soil_budget=soil,
        final_state=[dict(key=f.key,alive=f.alive,crown_remaining=f.crown_remaining,x_m=f.x,y_m=f.y,height_m=f.tree.height_m) for f in state],
        final_standing_count=sum(f.alive for f in state),event_rule='simultaneous_exceedances_at_each_prescribed_wind_increment',
        collision_and_fallen_debris='unmodelled',evolution='unmodelled'),field


def extra_studies(out):
    air=LowerAir();grid=Grid();forest,meta=make_patch('dome')
    xedge=max(abs(t.x) for t in forest);yedge=max(abs(t.y) for t in forest)
    altered=[replace(t,roots=replace(t.roots,plate_radius_m=max(30.,t.roots.plate_radius_m),
        plate_depth_m=max(4.,t.roots.plate_depth_m))) if abs(t.x)==xedge or abs(t.y)==yedge else t for t in forest]
    d=TreeDrag(forest,grid,air);v,F,D=SpatialFlow(grid).solve(d)
    adaptation=[]
    for label,case in [('allometric_roots',forest),('retained_edge_root_size',altered)]:
        drag=TreeDrag(case,grid,air)
        # Root changes leave the aerodynamic geometry and Cd-area quadrature unchanged.
        if not np.allclose(d.area,drag.area,rtol=0,atol=0):raise RuntimeError('Adaptation changed drag geometry')
        model=StandMechanics(case,air);R=model.evaluate(drag,F*4)
        adaptation.append(dict(label=label,meta=meta,metrics_at_20m_s=metrics(R),response=R,
            envelope=envelope(model,drag,F),forest=[asdict(t) for t in case],
            summed_root_plate_volume_m3=sum(np.pi*t.roots.plate_radius_m**2*t.roots.plate_depth_m for t in case),
            root_tissue_construction_cost='unmodelled'))
    js(out/'edge_adaptation.json',dict(flow=D,cases=adaptation))
    forest,_=make_patch('flat');one=[forest[len(forest)//2]]
    d=TreeDrag(one,grid,air);v,F,D=SpatialFlow(grid).solve(d)
    js(out/'single_tree_control.json',dict(flow=D,response=StandMechanics(one,air).evaluate(d,F*4),
        description='Same compliant-root tree in a solitary porous-flow control; not the old uniform-wind threshold'))
    soft=[replace(t,roots=replace(t.roots,cohesion_pa=5000.,pore_pressure_ratio=.6)) for t in forest]
    variations=[]
    for step in (1.,.25):
        print('weaker-soil ramp',step,flush=True)
        result,field=damage_sequence(soft,air,grid,Support(),np.arange(0.,30.+step/2,step))
        variations.append(dict(wind_increment_m_s=step,**result))
    js(out/'weaker_soil_damage.json',dict(forest=[asdict(t) for t in soft],grid=asdict(grid),runs=variations,
        note='Changing soil capacity also changes foundation stiffness under the stated capacity/rotation prescription'))


def run(out,stage='all'):
    out=output_guard(out);out.mkdir(parents=True,exist_ok=True)
    air=LowerAir();g=Grid();summaries=[];trees=[];configs=[];flow_diags=[];envelopes=[]
    fields=out/'fields';fields.mkdir(exist_ok=True)
    supports={'independent':Support(),'contact':Support(crown_mode='contact'),
              'clonal_roots':Support(root_link_nm_rad=2e9)}
    if stage in ('all','cases'):
        for shape in ('flat','tapered','dome','ramp','irregular','gap'):
            for matched in (True,False):
                if shape=='flat' and not matched:continue
                forest,meta=make_patch(shape,matched_mass=matched)
                for direction in ((0.,180.) if shape=='ramp' else (0.,)):
                    case=f'{shape}-{"matched" if matched else "same_peak"}-{int(direction)}'
                    print('flow',case,flush=True)
                    drag=TreeDrag(forest,g,air);v,F,D=SpatialFlow(g,direction).solve(drag)
                    flow_diags.append(dict(case=case,**D));np.savez_compressed(fields/(case+'.npz'),velocity_ratio=v,forces_at_10m_s=F)
                    configs.append(dict(case=case,meta=meta,forest=[asdict(f) for f in forest],grid=asdict(g)))
                    for label,support in supports.items():
                        if label!='independent' and (not matched or shape not in ('flat','dome')):continue
                        model=StandMechanics(forest,air,support);R=model.evaluate(drag,F*4);m=metrics(R)
                        summaries.append(dict(case=case,support=label,wind_m_s=20.,**meta,**m))
                        trees.extend(dict(case=case,support=label,wind_m_s=20.,**r) for r in R['trees'])

                        if matched:
                            env=envelope(model,drag,F);envelopes.append(dict(case=case,support=label,**env))
        csvout(out/'shape_comparison.csv',summaries);csvout(out/'individual_trees.csv',trees)
        js(out/'envelopes.json',envelopes);js(out/'flow_checks.json',flow_diags);js(out/'configurations.json',configs)
    if stage in ('all','sensitivity'):
        # Grid and finite-domain sensitivities preserve geometry and physical viscosity.
        forest,_=make_patch();sensitivity=[]
        options=[('coarse',replace(g,nx=30,ny=24,nz=15)),('base',g),
            ('fine',replace(g,nx=60,ny=48,nz=30)),
            ('wide_box',replace(g,lx=3600.,ly=2880.,nx=60,ny=48)),
            ('high_lid',replace(g,top=1800.,nz=30)),
            ('low_mixing',replace(g,mixing_length_m=6.)),('high_mixing',replace(g,mixing_length_m=24.))]
        for label,gg in options:
            print('sensitivity',label,flush=True)
            drag=TreeDrag(forest,gg,air);v,F,D=SpatialFlow(gg).solve(drag)
            R=StandMechanics(forest,air).evaluate(drag,F*4)
            sensitivity.append(dict(label=label,grid=asdict(gg),metrics_at_20m_s=metrics(R),flow=D))
        # Longer narrow patch: fetch dependence, with lateral bypass still available.
        deep,meta=make_patch(nx=25,ny=5);gg=replace(g,lx=4200.,nx=70)
        drag=TreeDrag(deep,gg,air);v,F,D=SpatialFlow(gg).solve(drag)
        deep_result=StandMechanics(deep,air).evaluate(drag,F*4)
        js(out/'long_patch.json',dict(meta=meta,flow=D,response=deep_result))
        js(out/'sensitivity.json',sensitivity)
    if stage in ('all','damage'):
        # Disturbance at finite wind increments. Repeatability and increment sensitivity are separate checks.
        cascades=[]
        for shape,support in [('flat',Support()),('flat',Support(crown_mode='contact')),('dome',Support())]:
            f,_=make_patch(shape);print('cascade',shape,support.crown_mode,flush=True)
            result,field=damage_sequence(f,air,g,support,np.arange(0.,35.01,1.))
            cascades.append(dict(shape=shape,support=asdict(support),**result))
        js(out/'damage_sequences.json',cascades)
    if stage in ('all','extras'):extra_studies(out)
    files=['climate/forest_flow.py','biosphere/forest_patch.py','biosphere/forest_mechanics.py','research/run_forest_patch.py']
    js(out/'source_hashes.json',{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files})
    load=lambda name:json.loads((out/name).read_text()) if (out/name).exists() else []
    flow_diags=load('flow_checks.json');cascades=load('damage_sequences.json')
    summaries=list(csv.DictReader((out/'shape_comparison.csv').open())) if (out/'shape_comparison.csv').exists() else []
    trees=list(csv.DictReader((out/'individual_trees.csv').open())) if (out/'individual_trees.csv').exists() else []
    summary=dict(base_commit=BASE,flow_cases=len(flow_diags),structural_cases=len(summaries),individual_records=len(trees),
        max_flow_momentum_residual=max((d['momentum_relative_residual'] for d in flow_diags),default=None),
        max_force_exchange_residual=max((d['force_ledger']['relative_residual'] for d in flow_diags),default=None),
        max_structural_force_balance_residual=max((float(r.get('horizontal_force_balance_relative_residual',0)) for r in summaries),default=None),
        cascades=[dict(shape=c['shape'],support=c['support']['crown_mode'],status=c['status'],events=len(c['events']),standing=c['final_standing_count']) for c in cascades],
        scientific_status='UNCALIBRATED_NEUTRAL_STEADY_SPATIAL_REQUIREMENTS_SCREEN',
        evolution='not_implemented',gusts='not_implemented',immersion='unchanged',python=platform.python_version())
    js(out/'summary.json',summary);print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=ROOT/'research/runs/forest_patch')
    parser.add_argument('--stage',choices=['all','cases','sensitivity','damage','extras'],default='all')
    args=parser.parse_args();run(args.out,args.stage)
