"""Quick gross build-out funding screen; no measured costs or demand forecast."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from .accounts import HERE, ROOT, producer, read_bound_inputs

INPUT = HERE/'buildout_inputs.json'


def results():
    spec = read_bound_inputs(INPUT,'terluna.humanities.cultural-buildout-inputs/1')
    sport = json.loads((ROOT/'humanities/economics/results/sports.json').read_text())
    entertainment = json.loads((ROOT/'humanities/economics/results/entertainment.json').read_text())
    c = sport['inputs']['touring_circuit']
    s = sport['touring_circuit']
    e = entertainment['cases']['expanded']
    extra = entertainment['inputs']['production']
    pay = spec['development_U_person_year']
    sports_base_pay = (s['capacity']['resident_athletes']*c['resident_athlete_pay_U_year']
        + c['resident_support_workers']*c['resident_support_pay_U_year']
        + c['visits_per_person_year']*(s['capacity']['visiting_athletes']*c['athlete_fee_U_visit']
                                      + c['visiting_support_workers']*c['support_fee_U_visit']))
    entertainment_base_pay = (sports_base_pay
        + extra['resident_workers']*entertainment['inputs']['resident_pay_U_year']
        + extra['visiting_workers']*extra['visits_year']*entertainment['inputs']['visitor_pay_U_visit'])
    pilot = spec['pilot']
    pilot_visits = pilot['people']*pilot['visits_each']
    pilot_travel = pilot_visits*c['return_fare_U']
    pilot_accommodation = pilot_visits*pilot['nights_each']*c['host_U_guest_night']
    circuit = spec['programmes']
    outputs = {}
    for name,row,base_pay,paid in (
        ('sports',s,sports_base_pay,s['support']['paid_viewers']),
        ('entertainment',e,entertainment_base_pay,e['audience']['paid_viewer_equivalents'])):
        fixed_purchases = sum(v['U_year'] for k,v in row['purchases'].items()
                              if k not in ('audience_acquisition','commercial_services'))
        cases = {}
        for case in spec['cases']:
            key = case['name']
            stages = {stage: v['person_years'][key]*pay+v['purchased_services_U'][key]
                      for stage,v in circuit.items()}
            stages['physical_pilot'] += pilot_travel+pilot_accommodation
            if name == 'sports':
                stages.pop('additional_entertainment_formats')
            development = sum(stages.values())
            duration = case['runway_months']/12
            uses = dict(
                prelaunch_development_U=development,
                installed_assets_U=row['installed_capital_U'],
                asset_contingency_U=row['installed_capital_U']*case['asset_contingency_fraction'],
                no_receipt_operating_runway_U=(base_pay+fixed_purchases)*duration,
                renewal_reserve_U=row['replacement_U_year']*duration,
                initial_paid_promotion_U=paid*case['paid_promotion_share']*spec['promotion_U_paid_customer'],
            )
            total = sum(uses.values())
            bond = sport['inputs']['finance']
            interest = bond['debt_U']*bond['real_rate']*case['construction_years']*case['average_debt_draw_fraction']
            fees = bond['debt_U']*bond['fee_fraction']
            debt_reserve = row['debt_service_U_year']*bond['reserve_years']
            cases[key] = dict(stages_U=stages, uses=uses,
                gross_uses_before_financing_U=total,
                development_person_years=sum(v['person_years'][key] for k,v in circuit.items()
                    if name=='entertainment' or k!='additional_entertainment_formats'),
                paid_promotion_share=case['paid_promotion_share'],
                blended_initial_promotion_U_customer=case['paid_promotion_share']*spec['promotion_U_paid_customer'],
                optional_inherited_debt=dict(principal_U=bond['debt_U'],
                    construction_interest_U=interest, issuance_fees_U=fees,
                    debt_reserve_U=debt_reserve,
                    other_funding_U=total+interest+fees+debt_reserve-bond['debt_U']))
        outputs[name] = dict(base_payroll_U_year=base_pay,
            fixed_purchases_U_year=fixed_purchases,
            paid_audience_target=paid,cases=cases)
    increments = {case['name']: {
        k:outputs['entertainment']['cases'][case['name']]['uses'][k]-outputs['sports']['cases'][case['name']]['uses'][k]
        for k in outputs['sports']['cases'][case['name']]['uses']}
        for case in spec['cases']}
    for row in increments.values():
        row['total_U'] = sum(row.values())
    exploratory = {case['name']:dict(
        person_years=spec['small_exploration']['person_years'][case['name']],
        purchased_U=spec['small_exploration']['purchased_U'][case['name']],
        total_U=spec['small_exploration']['person_years'][case['name']]*pay
                +spec['small_exploration']['purchased_U'][case['name']]) for case in spec['cases']}
    return dict(schema='terluna.humanities.cultural-buildout/1',
        producer=producer('humanities/economics/accounts.py','humanities/economics/buildout.py',
                          'humanities/economics/buildout_inputs.json'),
        evidence=spec['evidence'],reading_rule=spec['reading_rule'],
        units='Real U, person-years of effort, and months of funded operating exposure; no USD conversion.',
        inputs=spec,small_exploration=exploratory,programmes=outputs,
        entertainment_increment_after_sports=increments)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'research/runs/cultural_economics/buildout.json')
    args=parser.parse_args()
    data=results()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    print(args.output)


if __name__=='__main__':
    main()
