"""Conditional inventory/impulse trade at the bounded orbital pilot checkpoint."""
import json
from scipy.optimize import brentq
from shared.provenance import constants_used
from protection.dynamics.coverage_pilot import conditional_cost
from .cycling_search import HERE, ROOT, SCENARIO
from .fleet_run import digest


def main():
    path=HERE/'results/pilot_validation.json'
    validation=json.loads(path.read_text())
    n=validation['constraint_generation']['refined_solution']['inventory_tiles']
    rows=[]
    for inventory in [n,25e6,35e6]:
        def residual(dv):
            cost=conditional_cost(inventory,dv,SCENARIO['propulsion'],1.2,.1)
            return cost['propulsion_TW']-238.35 if cost['closed'] else 1e12
        rows.append(dict(tiles=inventory,break_even_dv_m_s_day=brentq(residual,0,100),
            three_m_s_day=conditional_cost(inventory,3.,SCENARIO['propulsion'],1.2,.1)))
    source_paths=['research/studies/solar_shield_array/pilot_budget.py',
        'protection/dynamics/coverage_pilot.py','research/studies/solar_shield_array/scenario.json']
    product=dict(schema='terluna.research.coverage-pilot-budget/1',
        producer=dict(source_hashes={p:digest(ROOT/p) for p in source_paths},
            constants=constants_used([p for p in source_paths if p.endswith('.py')]),
            input_sha256=digest(path)),
        evidence='Conditional scalar inventory/impulse accounting, with mass feedback; no achieved recurring operation',
        fixed_extra_mass_fraction=.2,peak_duty=.1,propulsion=SCENARIO['propulsion'],
        benchmark_TW=238.35,rows=rows,collection_capacity_W=None,delivered_power_W=None,
        scope='The first inventory is the restricted capacity estimate. The larger inventories are sensitivity scenarios with no coverage claim. Added mass has not been fed back into sail trajectories. Peak hardware capacity is a load rating, not electrical collection.')
    (HERE/'results/pilot_budget.json').write_text(json.dumps(product,indent=2)+'\n')
    print(json.dumps(rows),flush=True)


if __name__=='__main__':main()
