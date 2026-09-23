"""Quasi-static damage propagation with geometry/airflow re-solves.

Each wave applies simultaneous overload decisions. Removed biomass enters an
explicit debris ledger, but debris collisions, blockage and regeneration are
not simulated. A large-displacement or unstable state stops interpretation.
"""
from __future__ import annotations
from dataclasses import replace
import numpy as np
from biosphere.forest_patch import PatchDrag,tree_mass
from biosphere.forest_patch_mechanics import PatchStructure,InteractionConfig


def damage_sequence(air,grid,initial_trees,speeds,interactions=InteractionConfig(),
                    max_waves=60,first_batch_only=False):
    speeds=np.asarray(speeds,float)
    if speeds.ndim!=1 or len(speeds)==0 or np.any(~np.isfinite(speeds)) or np.any(speeds<0) or np.any(np.diff(speeds)<0):
        raise ValueError('Finite nonnegative nondecreasing wind sequence required')
    trees=list(initial_trees);broken=set();events=[];stages=[];flow=None;drag=None;dirty=True
    original_mass=sum(tree_mass(p.effective_tree()) for p in trees if p.alive)
    rho=float(air.sample(0)['density_kg_m3']);flow_count=0
    status='COMPLETED';last_records=[]
    for speed in speeds:
        for wave in range(max_waves):
            if not any(p.alive for p in trees):break
            if dirty:
                drag=PatchDrag(grid,trees)
                flow=grid.solve(drag.coefficient,None if flow is None else flow.velocity)
                flow_count+=1;dirty=False
            loads,ledger=drag.loads(flow,float(speed),rho)
            model=PatchStructure(air,trees,interactions,broken)
            result=model.solve(loads);last_records=result['trees']
            stages.append(dict(speed_m_s=float(speed),wave=wave,standing=sum(p.alive for p in trees),
                force_ledger=ledger,flow=flow.diagnostics,structure_status=result['status'],
                domain_limited=result['domain_limited'],
                linear_residual=result.get('linear_equilibrium_relative_residual'),
                tree_records=last_records))
            if result['domain_limited']:
                status='MODEL_DOMAIN_LIMIT';break
            newbroken={l['id'] for l in result['links'] if l['utilization']>=1 and l['id'] not in broken}
            fatalities={r['id'] for r in result['trees'] if max(r['root_utilization'],r['stem_utilization'])>=1}
            shedding={r['id'] for r in result['trees'] if r['crown_utilization']>=1 and r['id'] not in fatalities}
            if not (newbroken or fatalities or shedding):break
            broken.update(newbroken)
            changes=[];newtrees=[]
            for p in trees:
                if p.identifier in fatalities:
                    p=replace(p,alive=False);changes.append(dict(id=p.identifier,kind='tree_failure'))
                elif p.identifier in shedding:
                    remaining=.35 if p.remaining_crown>.35 else 0.
                    p=replace(p,remaining_crown=remaining);changes.append(dict(id=p.identifier,kind='crown_shedding',remaining=remaining))
                newtrees.append(p)
            trees=newtrees
            events.append(dict(speed_m_s=float(speed),wave=wave,changes=changes,
                broken_links=sorted(newbroken),standing_after=sum(p.alive for p in trees)))
            dirty=bool(fatalities or shedding)
            if first_batch_only:status='FIRST_BATCH_CONTROL';break
        else:
            status='WAVE_LIMIT';break
        if status!='COMPLETED':break
    remaining=sum(tree_mass(p.effective_tree()) for p in trees if p.alive)
    debris=original_mass-remaining
    return dict(status=status,events=events,stages=stages,flow_solves=flow_count,
        original_tree_count=sum(p.alive for p in initial_trees),standing_tree_count=sum(p.alive for p in trees),
        crown_damaged_standing=sum(p.alive and p.remaining_crown<1 for p in trees),
        remaining_ids=[p.identifier for p in trees if p.alive],broken_links=sorted(broken),
        initial_mass_kg=original_mass,standing_mass_kg=remaining,debris_mass_kg=debris,
        mass_ledger_residual_kg=original_mass-remaining-debris,
        damage_order='simultaneous overloads within each quasi-static wave',
        limitations=['no gust or collision dynamics','no fallen-debris aerodynamic blockage',
            'no plasticity or post-buckling','soil capacity remains reserved for failed root zones',
            'undeformed flow geometry between damage events','fixed crown drag coefficient'])
