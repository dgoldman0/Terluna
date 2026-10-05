"""Conditional electrical accounts on executed intervals and total resource use."""
import argparse
import json
import resource
import time

from shared.provenance import constants_used
from engineering.shield_power import routed_power
from .repair_screen import RUN,source_parents,write_product
from .fleet_run import digest
from .cycling_search import HERE,SCENARIO


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--earlier-summary-charge-s',type=float,default=0.)
    args=parser.parse_args()
    if args.earlier_summary_charge_s<0:raise ValueError('Negative resource charge')
    before=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    names=['repair_screen.json','repair_local_screen.json','repair_local_refined.json',
        'repair_trial_local_refined.json','repair_trial_local_refined_fine.json',
        'repair_trial_roll_60.json','repair_trial_roll_minus30.json',
        'repair_trade.json','repair_budget.json','repair_audit.json']
    p,sources=source_parents(names,['research/studies/solar_shield_array/repair_synthesis.py'])
    if not all(x['completed'] for x in p.values()):raise ValueError('Incomplete parent')
    out=dict(schema='terluna.research.return-repair-synthesis/1',
        producer=dict(source_hashes=sources,constants=constants_used(sources),
            inputs={n:digest(HERE/'results'/n) for n in names}),
        accepted_return=False,accepted_cycle=False,accepted_fleet=False,new_target_adopted=False,
        actual_delivered_power_W=None,full_cycle_energy_J=None,full_cycle_average_W=None,
        scope='Conversion/routing sensitivity of each actual executed prefix; temporal storage and power availability remain unspecified. Positive interval energy does not establish instantaneous self-sufficiency or full-cycle delivery.',
        trials=[])
    for trial in p['repair_budget.json']['trials']:
        full=trial['whole_executed_prefix'];duration=full['duration_s']
        ret=trial['evaluations'][0]['cost'];rt=trial['end_s']-trial['start_s']
        rb=trial['grazing_refinement'][-1]['optical']['redirected_band_optical_W']['mean_W']
        cases=[]
        for e in p['repair_trade.json']['conversion_cases']:
            row=routed_power(full['redirected_band_energy_J']/duration,e['capture'],e['conversion'],e['delivery'],
                full['propulsion_energy_J']/duration)
            row.update(case=e['label'],interval_net_delivered_potential_J=row['delivered_W']*duration,
                return_mean_gross_electric_W=rb*e['capture']*e['conversion'],
                return_mean_propulsion_W=ret['propulsion_energy_J']/rt,
                return_interval_energy_shortfall_J=max(ret['propulsion_energy_J']-rb*rt*e['capture']*e['conversion'],0.),
                propulsion_nonjet_loss_W=(1-SCENARIO['propulsion']['efficiency'])*full['propulsion_energy_J']/duration)
            cases.append(row)
        out['trials'].append(dict(source=trial['source'],duration_s=duration,electricity_sensitivities=cases))
    interrupted=(p['repair_budget.json'].get('accounting_continuation') or {}).get('prior_attempt_wall_s',0.)
    prior=10.+interrupted+args.earlier_summary_charge_s
    out['resources']=dict(initial_inspection_charge_s=10.,interrupted_accounting_charge_s=interrupted,
        earlier_summary_charge_s=args.earlier_summary_charge_s,
        completed_producer_wall_s=sum(x['elapsed_s'] for x in p.values()),
        numerical_wall_s_before_this_summary=prior+sum(x['elapsed_s'] for x in p.values()),
        ceiling_s=3600.,one_calculation_process=True,numerical_threads=1,address_space_cap_GiB=2,
        maximum_recorded_rss_MiB=max(x['max_rss_MiB'] for x in p.values()),
        raw_output_MiB=sum(f.stat().st_size for f in RUN.iterdir() if f.is_file())/2**20,
        raw_ceiling_MiB=512,checks_have_separate_allowance=True)
    if out['resources']['numerical_wall_s_before_this_summary']+time.monotonic()-before>=3600:
        raise ValueError('Numerical budget exceeded')
    out['completed']=True;write_product('repair_synthesis.json',out,sources,before)
    print(json.dumps(out['resources']),flush=True)


if __name__=='__main__':main()
