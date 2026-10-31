"""Preliminary enclosed-venue operating cash and project-debt scenarios.

All quantities are conditional. U is a real numeraire, not an adopted currency.
The cash account includes wages and cannot be interpreted as value added.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from shared.constants import JULIAN_YEAR_DAYS
from shared.provenance import ROOT, constants_used
from .finance import annuity_payment, debt_service_first_year

HERE = Path(__file__).resolve().parent
INPUT_PATH = HERE/'sports_finance_inputs.json'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def nonnegative(*values):
    if any(isinstance(v,bool) or not math.isfinite(v) or v<0 for v in values):
        raise ValueError('Finite nonnegative quantities required')


def fraction(value):
    nonnegative(value)
    if value>1:raise ValueError('Fraction must not exceed one')


def visitor_capacity(venue, stay_nights):
    nonnegative(*[venue[k] for k in ('guest_beds','occupancy_ceiling','team_guest_nights',
        'route_arrivals_year','training_completions_year','arena_days_year','arena_hours_day',
        'arena_simultaneous_participants','professional_participant_hours','guest_arena_hours')])
    fraction(venue['occupancy_ceiling'])
    if stay_nights<=0 or venue['guest_arena_hours']<=0:
        raise ValueError('Positive visit duration and arena use required')
    nights=venue['guest_beds']*JULIAN_YEAR_DAYS*venue['occupancy_ceiling']
    arena=venue['arena_days_year']*venue['arena_hours_day']*venue['arena_simultaneous_participants']
    if venue['team_guest_nights']>nights or venue['professional_participant_hours']>arena:
        raise ValueError('Professional use exceeds the planned facility allowance')
    caps=dict(beds=(nights-venue['team_guest_nights'])/stay_nights,
        travel=venue['route_arrivals_year'],training=venue['training_completions_year'],
        arena=(arena-venue['professional_participant_hours'])/venue['guest_arena_hours'])
    return dict(caps_arrivals_year=caps,limiting_arrivals_year=min(caps.values()),
        usable_guest_nights=nights-venue['team_guest_nights'],
        professional_guest_nights=venue['team_guest_nights'],
        professional_arena_participant_hours=venue['professional_participant_hours'])


def operating_case(spec, inherited, case):
    capacity=visitor_capacity(spec['venue'],inherited['stay_nights'])
    for key in ('visitor_demand_year','rights_receipts_U','sponsorship_receipts_U','capital_scale'):
        nonnegative(case[key])
    visitors=min(case['visitor_demand_year'],capacity['limiting_arrivals_year'])
    onsite=visitors*inherited['stay_nights']*inherited['destination_U_per_night']
    media=case['rights_receipts_U']+case['sponsorship_receipts_U']
    costs=spec['costs']
    for key in ('visitor_cash_cost_fraction','media_cash_cost_fraction','tax_allowance_fraction'):
        fraction(costs[key])
    visitor_capital=spec['venue']['guest_beds']*inherited['visitor_capital_U_per_bed']*case['capital_scale']
    production_capital=spec['production_capital_U']*case['capital_scale']
    capital=visitor_capital+production_capital
    replacement=visitor_capital/inherited['visitor_asset_life_years']+production_capital/spec['production_asset_life_years']
    variable=onsite*costs['visitor_cash_cost_fraction']+media*costs['media_cash_cost_fraction']
    fixed=capital*costs['fixed_operating_fraction_capital']+costs['core_programme_U_year']+costs['host_services_U_year']
    before_tax=onsite+media-variable-fixed-replacement
    tax=max(0,before_tax)*costs['tax_allowance_fraction']
    return dict(visitor_capacity=capacity,realized_visitors_year=visitors,
        visitor_nights_year=visitors*inherited['stay_nights'],onsite_receipts_U=onsite,
        rights_receipts_U=case['rights_receipts_U'],sponsorship_receipts_U=case['sponsorship_receipts_U'],
        total_receipts_U=onsite+media,variable_cash_cost_U=variable,fixed_cash_cost_U=fixed,
        installed_capital_U=capital,replacement_U=replacement,pre_tax_cash_U=before_tax,
        tax_allowance_U=tax,cash_available_for_debt_service_U=before_tax-tax)


def bond_capacity(cash,capital_and_working_cash,rate,years,coverage,reserve_years,fee_fraction,investor_demand):
    nonnegative(capital_and_working_cash,rate,reserve_years,investor_demand)
    fraction(fee_fraction)
    if not math.isfinite(cash) or not math.isfinite(coverage) or coverage<1:
        raise ValueError('Finite cash and coverage at least one required')
    a=annuity_payment(1,rate,years)
    retained=reserve_years*a+fee_fraction
    if retained>=1:raise ValueError('Reserve and fees consume all issue proceeds')
    limits=dict(operating_cash=max(0,cash)/(coverage*a),investor_demand=investor_demand,
                total_uses=capital_and_working_cash/(1-retained))
    debt=min(limits.values())
    service=debt*a;reserve=reserve_years*service;fees=fee_fraction*debt
    other=capital_and_working_cash+reserve+fees-debt
    return dict(principal_limits_U=limits,gross_principal_U=debt,annual_service_U=service,
        debt_service_reserve_U=reserve,issuance_cost_U=fees,net_proceeds_for_capital_and_working_cash_U=debt-reserve-fees,
        other_funding_required_U=other,rate_real=rate,years=years,coverage_target=coverage,
        first_year=debt_service_first_year(debt,rate,years),
        sources_minus_uses_U=debt+other-capital_and_working_cash-reserve-fees)


def fixed_bond_stress(cash,bond):
    """One independent stress year; opening debt reserve resets between cases."""
    service=bond['annual_service_U'];reserve=bond['debt_service_reserve_U']
    nonnegative(service,reserve)
    if not math.isfinite(cash):raise ValueError('Finite cash required')
    operating_gap=max(0,-cash)
    debt_gap=max(0,service-max(0,cash))
    draw=min(reserve,debt_gap)
    return dict(operating_coverage=cash/service if service else None,
        cash_after_debt_before_reserve_U=cash-service,
        operating_funding_gap_U=operating_gap,debt_gap_before_reserve_U=debt_gap,
        debt_reserve_draw_U=draw,closing_debt_reserve_U=reserve-draw,
        uncovered_debt_payment_U=debt_gap-draw,funding_gap_after_debt_reserve_before_working_cash_U=operating_gap+debt_gap-draw,
        cash_after_service_positive_U=max(0,cash-service))


def load_inputs(path=INPUT_PATH,root=ROOT):
    spec=json.loads(Path(path).read_text())
    if spec['schema']!='terluna.research.sports-finance-inputs/2':raise ValueError('Unexpected input schema')
    bound={}
    for name,identity in spec['products'].items():
        p=Path(root)/name
        if digest(p)!=identity['sha256']:raise ValueError(f'Input changed; review and rebind: {name}')
        data=json.loads(p.read_text())
        if data['schema']!=identity['schema']:raise ValueError(f'Unexpected schema: {name}')
        bound[name]=data
    old=next(iter(bound.values()))['scenarios']['reference'][0]
    inherited=dict(stay_nights=old['realized_visitor_nights']/old['realized_trips'],
        destination_U_per_night=old['tourism']['sales_U']/old['realized_visitor_nights'],
        visitor_capital_U_per_bed=old['tourism']['initial_capital_U']/old['visitor_bed_places'],
        visitor_asset_life_years=old['tourism']['initial_capital_U']/old['tourism']['cfc_and_replacement_U'])
    return spec,inherited


def results():
    s,old=load_inputs();f=s['bond'];ops={name:operating_case(s,old,c) for name,c in s['cases'].items()}
    base=ops['reference'];cash=base['cash_available_for_debt_service_U'];uses=base['installed_capital_U']+s['working_cash_U']
    def issue(rate,years,demand):
        return bond_capacity(cash,uses,rate,years,f['coverage'],f['reserve_years'],f['fee_fraction'],demand)
    bond=issue(f['real_rate'],f['years'],f['investor_demand_U'])
    cases={name:dict(operation=row,additional_upfront_capital_U=max(0,row['installed_capital_U']-base['installed_capital_U']),
        fixed_reference_bond=fixed_bond_stress(row['cash_available_for_debt_service_U'],bond)) for name,row in ops.items()}
    c=s['costs']
    if c['media_cash_cost_fraction']>=1 or c['tax_allowance_fraction']>=1:
        raise ValueError('Finite break-even examples require media cost and tax fractions below one')
    onsite_contribution=base['onsite_receipts_U']*(1-c['visitor_cash_cost_fraction'])
    fixed=base['fixed_cash_cost_U']+base['replacement_U']
    need_media=max(0,(fixed-onsite_contribution)/(1-c['media_cash_cost_fraction']))
    need_debt_media=max(0,(fixed+bond['annual_service_U']/(1-c['tax_allowance_fraction'])-onsite_contribution)/(1-c['media_cash_cost_fraction']))
    files=['research/studies/world_economy/sports_finance.py','research/studies/world_economy/finance.py',
           'research/studies/world_economy/sports_finance_inputs.json','research/studies/world_economy/sports_finance_sources.json']
    return dict(schema='terluna.research.sports-finance/2',evidence=s['evidence'],reading_rule=s['reading_rule'],
        assumptions=s,inherited_screen=old,
        producer=dict(files={p:digest(ROOT/p) for p in files},constants=constants_used(files)),
        reference_bond=bond,operating_scenarios=cases,
        break_even=dict(combined_media_and_sponsor_receipts_before_debt_U=need_media,
            combined_media_and_sponsor_receipts_to_pay_reference_debt_U=need_debt_media,
            rights_audience_equivalents_before_debt=max(0,need_media-base['sponsorship_receipts_U'])/s['rights_receipt_U_per_audience_equivalent'],
            rights_audience_equivalents_to_pay_debt=max(0,need_debt_media-base['sponsorship_receipts_U'])/s['rights_receipt_U_per_audience_equivalent']),
        financing_sensitivity={f'{rate:g}_{years}':issue(rate,years,f['investor_demand_U']) for rate,years in f['rate_term_cases']},
        investor_demand_sensitivity={str(demand):issue(f['real_rate'],f['years'],demand) for demand in f['demand_cases_U']})


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'research/runs/world_economy/sports_finance.json');args=p.parse_args()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(results(),indent=2,allow_nan=False)+'\n');print(args.output)


if __name__=='__main__':main()
