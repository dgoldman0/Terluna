"""Frozen faster-rollout orders for the combined reclamation comparison."""
from copy import deepcopy
import json
from pathlib import Path
from . import markets,bottleneck_orders
HERE=Path(__file__).resolve().parent

def product():
    cfg,q,d,b,_,_=markets.inputs();cases=[]
    for portfolio in ('services','propellant_logistics','diversified'):
        for growth in ('reference','wider'):
            cc=deepcopy(cfg);cc.update(commercial_capture_growth=.08,service_staff_fraction=.2)
            case=bottleneck_orders.scheduled_case(cc,q,d,b,portfolio,'breakthrough',growth,cfg['market_growth'][growth],4)
            name=case['name'];case['name']+=':combined_rollout'
            case['experiment']={'baseline_name':name,'variant':'combined_rollout','engineering_overrides':{}}
            cases.append(case)
    return {'schema':'terluna.humanities.enclosed-bottleneck-orders/1','producer':{'path':'research/studies/enclosed_returns/revised_orders.py','sha256':markets.digest(__file__)},
            'dependencies':{path:markets.digest(HERE/path) for path in ('markets.py','bottleneck_orders.py','enclosed_inputs.json',*cfg['inherited_inputs'])},
            'evidence':'Conditional orders under the previously studied faster supplier/service rollout; all new physical schedules require independent capacity checks before financial admission.',
            'rules':['Retain worldwide buyer budgets, prices, qualification, synergies, growth and residential schedule. Commercial offer growth rises from 4 to 8 percent annually; service staffing allowance from 10 to 20 percent of residents; service-package orders from one to four per year. These are planning interventions, not demonstrated productivity gains.',
                     'The physical producer separately combines these requests with recovery and refurbishment. No revenue from a failed physical schedule is admitted.'],
            'cases':cases}

if __name__=='__main__':
    p=HERE/'results/revised_orders.json';p.write_text(json.dumps(product(),sort_keys=True,separators=(',',':'),allow_nan=False)+'\n');print(p,markets.digest(p))
