"""A small resident league selling to off-Moon audiences, without tourism."""
from __future__ import annotations

import argparse
from copy import deepcopy
import json
from pathlib import Path

from shared.constants import JULIAN_YEAR_DAYS

from .accounts import (HERE, ROOT, account, count, fraction, nonnegative,
                       payment, positive, producer, read_bound_inputs)

INPUT = HERE / 'sports_inputs.json'
SCHEMA = 'terluna.humanities.sports-economics/1'


def schedule(spec):
    c = count(spec['clubs'])
    if c < 2:
        raise ValueError('Competition needs at least two clubs')
    athletes = c * count(spec['roster'])
    matches = c * (c-1) // 2 * count(spec['meetings_per_pair']) + count(spec['playoff_matches'])
    training = athletes * positive(spec['training_hours_per_athlete']) / count(spec['training_group'])
    match_hours = matches * positive(spec['hours_per_match_slot'])
    other_hours = nonnegative(spec['maintenance_hours']) + nonnegative(spec['community_hours'])
    available = positive(spec['arena_days']) * positive(spec['arena_hours_per_day'])
    used = training + match_hours + other_hours
    if used > available:
        raise ValueError('Arena use exceeds annual scheduled hours')
    workers = athletes + count(spec['production_workers']) + count(spec['other_workers'])
    residents = count(spec['host_residents'])
    if workers > residents:
        raise ValueError('League workforce exceeds host population')
    return dict(clubs=c, athletes=athletes, workers=workers,
                host_residents=residents, workforce_share=workers/residents,
                matches_year=matches, training_arena_hours=training,
                competition_arena_hours=match_hours, occupied_arena_hours=used,
                available_arena_hours=available, arena_utilization=used/available)


def evaluate(spec, case):
    s = schedule(spec['season'])
    a = nonnegative(case['paying_audience_year'])
    y = positive(case['rights_yield_U'])
    retention = fraction(case['retention'])
    acquisition_cost = nonnegative(case['acquisition_U_per_new_customer'])
    new_customers = a * (1-retention)
    cost = spec['costs']
    wage_scale = positive(case.get('wage_scale', 1))
    athlete_pay = s['athletes'] * nonnegative(cost['athlete_compensation_U_year']) * wage_scale
    wages = athlete_pay + wage_scale * (
        spec['season']['production_workers'] * nonnegative(cost['production_compensation_U_year'])
        + spec['season']['other_workers'] * nonnegative(cost['other_compensation_U_year']))
    host = s['workers'] * nonnegative(cost['host_service_U_per_worker_year'])
    media = a*y
    sponsor = nonnegative(case['sponsor_U_year'])
    licensing = nonnegative(case['royalties_U_year'])
    membership = nonnegative(case.get('membership_U_year', 0))
    remote = media + sponsor + licensing + membership
    wage_floor = remote * fraction(case.get('payroll_revenue_floor_fraction', 0))
    if wage_floor > wages:
        athlete_pay *= wage_floor/wages
        wages = wage_floor
    variable_fraction = fraction(case['variable_cost_fraction'])
    receipts = {
        'media': dict(U_year=media, export_fraction=1),
        'sponsor': dict(U_year=sponsor, export_fraction=1),
        'game_and_merchandise_royalties': dict(U_year=licensing, export_fraction=1),
        'additional_membership_services': dict(U_year=membership, export_fraction=1),
        'local_admissions': dict(U_year=nonnegative(case['local_sales_U_year']), export_fraction=0),
    }
    purchases = {
        'host_services': dict(U_year=host, import_fraction=0),
        'facility_operations': dict(U_year=nonnegative(cost['facility_operations_U_year'])),
        'event_services': dict(U_year=s['matches_year']*nonnegative(cost['event_services_U_match'])),
        'customer_acquisition': dict(U_year=new_customers*acquisition_cost),
        'variable_commercial_services': dict(U_year=remote*variable_fraction),
    }
    # Import shares are scenario assumptions, recorded explicitly in the inputs.
    for name, row in purchases.items():
        row['import_fraction'] = fraction(cost['purchase_import_fractions'][name])
    cap = deepcopy(spec['capital'])
    cap['installed_U'] *= positive(case.get('capital_scale', 1))
    out = account(receipts, wages, purchases, cap, spec['finance'],
                  spec['tax_fraction'], spec['supplier_import_fraction'])
    club_grants = (media+sponsor)*fraction(spec['club_revenue_fraction'])
    club_costs = athlete_pay + host*s['athletes']/s['workers']
    out['clubs'] = dict(total_internal_grants_U_year=club_grants,
        grant_per_club_U_year=club_grants/s['clubs'],
        payroll_and_host_per_club_U_year=club_costs/s['clubs'],
        aggregate_surplus_before_tax_U_year=club_grants-club_costs,
        central_surplus_before_tax_U_year=out['net_operating_surplus_U_year']-(club_grants-club_costs))
    out['audience'] = dict(paying_customer_years=a, new_customers_to_replace_churn=new_customers,
        annual_acquisition_U=new_customers*acquisition_cost,
        initial_customer_acquisition_U=a*acquisition_cost,
        required_reach_people=a/(positive(spec['engagement_fraction'])*positive(spec['paid_conversion_fraction'])))
    for k in ('engagement_fraction', 'paid_conversion_fraction'):
        fraction(spec[k])
    out['schedule'] = s
    return out


def audience_threshold(spec, case, debt_coverage):
    """Hold sponsorship, royalties, price and capacity fixed; vary paid audience."""
    nonnegative(debt_coverage)
    if case.get('payroll_revenue_floor_fraction', 0):
        raise ValueError('Analytic threshold requires fixed payroll')
    zero = dict(case, paying_audience_year=0)
    base = evaluate(spec, zero)
    marginal = case['rights_yield_U']*(1-case['variable_cost_fraction']) - (
        1-case['retention'])*case['acquisition_U_per_new_customer']
    if marginal <= 0:
        return None
    tax = fraction(spec['tax_fraction'])
    if tax == 1:
        raise ValueError('Tax fraction must be below one for debt thresholds')
    target = debt_coverage*base['debt_service_U_year']/(1-tax)
    return max(0, (target-base['net_operating_surplus_U_year'])/marginal)


def fan_support(spec, row):
    """Reach is independent of host population; each payment has its own service.

    Membership buyers are a subset of paid viewers, receiving additional services.
    Gifts are transfers and fan investment is financing, both outside revenue.
    """
    fans = nonnegative(row['fans'])
    paid = fraction(row['paid_viewer_share'])
    members = fraction(row['membership_share'])
    if members > paid:
        raise ValueError('Membership cohort must be a subset of paid viewers')
    case = dict(spec['cases']['reference'], paying_audience_year=fans*paid,
        membership_U_year=fans*members*nonnegative(row['membership_owner_yield_U']),
        royalties_U_year=fans*nonnegative(row['royalty_U_per_fan']),
        sponsor_U_year=fans*nonnegative(row['sponsor_U_per_fan']),
        variable_cost_fraction=row['variable_cost_fraction'],
        wage_scale=row['wage_scale'], capital_scale=row['capital_scale'],
        payroll_revenue_floor_fraction=row['payroll_revenue_floor_fraction'])
    out = evaluate(spec, case)
    gifts = fans*fraction(row['patron_share'])*nonnegative(row['patron_gift_U'])
    investors = fans*fraction(row['founding_investor_share'])
    investment = investors*nonnegative(row['investment_per_fan_U'])
    out['support'] = dict(fans=fans, paid_viewers=fans*paid,
        membership_buyers=fans*members,
        earned_export_U_per_fan=out['exports_U_year']/fans if fans else 0,
        patron_transfers_U_year=gifts, founding_investors=investors,
        one_time_fan_investment_U=investment,
        cash_after_debt_and_patron_transfers_U_year=out['cash_after_debt_U_year']+gifts)
    return out


def touring_circuit(spec, circuit=None):
    """Four candidate competitions, overlapping fans and visiting-team blocks.

    Registered athletes can live on Earth between lunar appearances. This is
    an annual/peak capacity screen, not a qualified travel or match timetable.
    """
    c = circuit or spec['touring_circuit']
    n = len(c['competitions'])
    clubs = count(c['clubs_per_competition'])
    roster = count(c['roster'])
    resident_roster = count(c['resident_athletes_per_club'])
    if resident_roster > roster:
        raise ValueError('Resident roster exceeds club roster')
    athletes = n*clubs*roster
    resident_athletes = n*clubs*resident_roster
    visitors = athletes-resident_athletes
    crew = count(c['resident_support_workers'])
    visiting_crew = count(c['visiting_support_workers'])
    resident_workers = resident_athletes+crew
    rotations = count(c['visits_per_person_year'])
    stay = positive(c['nights_per_visit'])
    arrivals = (visitors+visiting_crew)*rotations
    nights = arrivals*stay
    beds = count(c['guest_beds'])
    bed_allowance = beds*JULIAN_YEAR_DAYS*fraction(c['guest_occupancy_ceiling'])
    venues = count(c['venues'])
    concurrent = count(c['simultaneous_competitions'])
    if concurrent > min(n, venues):
        raise ValueError('Concurrent competitions exceed venues')
    peak_guests = (visitors+visiting_crew)/n*concurrent
    games_each = clubs*(clubs-1)//2*count(c['meetings_per_pair'])+count(c['playoff_matches'])
    matches = n*games_each
    resident_training = resident_athletes*positive(c['resident_arena_training_hours'])/count(c['training_group'])
    visitor_training = visitors*rotations*positive(c['visitor_training_hours_per_visit'])/c['training_group']
    match_hours = matches*positive(c['hours_per_match_slot'])
    shared_hours = venues*(nonnegative(c['maintenance_hours_per_venue'])+nonnegative(c['community_hours_per_venue']))
    used = resident_training+visitor_training+match_hours+shared_hours
    available = venues*positive(c['arena_days'])*positive(c['arena_hours_per_day'])
    blocks = (n+concurrent-1)//concurrent*rotations
    if blocks*stay > c['arena_days']:
        raise ValueError('Visiting blocks exceed programme calendar')
    per_block = ((games_each+rotations-1)//rotations)*c['hours_per_match_slot'] + (
        visitors/n*c['visitor_training_hours_per_visit']/c['training_group'])
    if (nights > bed_allowance or peak_guests > beds or arrivals > c['arrival_capacity_year']
            or used > available or per_block > stay*c['arena_hours_per_day']):
        raise ValueError('Touring programme exceeds a capacity allowance')
    if resident_workers > count(c['host_residents']):
        raise ValueError('Resident workforce exceeds host population')
    populations = c['audience_groups']
    fan_by_sport = {name: 0 for name in c['competitions']}
    unique_fans = 0
    for group in populations:
        fans = nonnegative(group['fans'])
        names = group['competitions']
        if not names or len(names) != len(set(names)) or not set(names) <= set(fan_by_sport):
            raise ValueError('Invalid audience-group membership')
        unique_fans += fans
        for name in names:
            fan_by_sport[name] += fans
    support_row = dict(spec['fan_support_cases']['hundred_million_fans'], fans=unique_fans)
    support = fan_support(spec, support_row)
    receipts = support['receipts']
    remote = support['exports_U_year']
    resident_pay = resident_athletes*c['resident_athlete_pay_U_year']+crew*c['resident_support_pay_U_year']
    visiting_pay = rotations*(visitors*c['athlete_fee_U_visit']+visiting_crew*c['support_fee_U_visit'])
    base_pay = resident_pay+visiting_pay
    scale = max(1,remote*support_row['payroll_revenue_floor_fraction']/positive(base_pay))
    wages = base_pay*scale
    nonresident_pay = visiting_pay*scale
    purchases = {
        'resident_host_services': dict(U_year=resident_workers*c['host_U_resident_year'],import_fraction=0),
        'guest_host_services': dict(U_year=nights*c['host_U_guest_night'],import_fraction=0),
        'return_travel': dict(U_year=arrivals*c['return_fare_U'],import_fraction=1),
        'venue_operations': dict(U_year=venues*c['venue_operations_U_year'],import_fraction=0.4),
        'event_services': dict(U_year=matches*c['event_services_U_match'],import_fraction=0.5),
        'audience_acquisition': dict(U_year=support['audience']['annual_acquisition_U'],import_fraction=0.9),
        'commercial_services': dict(U_year=remote*support_row['variable_cost_fraction'],import_fraction=0.6),
    }
    capital = dict(spec['capital'], installed_U=venues*c['capital_U_venue']+beds*c['capital_U_guest_bed'])
    out = account(receipts,wages,purchases,capital,spec['finance'],spec['tax_fraction'],
                  spec['supplier_import_fraction'],nonresident_pay)
    out['capacity'] = dict(competitions=n, clubs=n*clubs, registered_athletes=athletes,
        resident_athletes=resident_athletes, visiting_athletes=visitors,
        resident_workers=resident_workers, host_residents=c['host_residents'],
        visitor_arrivals_year=arrivals, guest_nights_year=nights,
        guest_night_allowance_year=bed_allowance, peak_guests=peak_guests,
        matches_year=matches, occupied_arena_hours=used, available_arena_hours=available,
        block_days=blocks*stay, busiest_competition_block_hours=per_block,
        available_hours_per_competition_block=stay*c['arena_hours_per_day'])
    out['audience'] = dict(unique_fans=unique_fans, fans_by_competition=fan_by_sport,
        sum_of_sport_followings=sum(fan_by_sport.values()),
        reading_rule='Audience groups are mutually exclusive cohorts; their sport memberships overlap. One bundle is sold once per customer.')
    out['support'] = support['support']
    # Recalculate after-debt support cash for this larger physical programme.
    out['support']['cash_after_debt_and_patron_transfers_U_year'] = out['cash_after_debt_U_year']+out['support']['patron_transfers_U_year']
    return out


def results():
    spec = read_bound_inputs(INPUT, 'terluna.humanities.sports-economics-inputs/1')
    ref = spec['cases']['reference']
    cases = {name: evaluate(spec, row) for name, row in spec['cases'].items()}
    host_double = deepcopy(spec)
    host_double['costs']['host_service_U_per_worker_year'] *= 2
    automated = deepcopy(spec)
    automated['costs']['facility_operations_U_year'] += 1000000
    alternatives = dict(
        host_services_double=evaluate(host_double, ref),
        half_payroll_with_extra_million_services=evaluate(automated, dict(ref, wage_scale=0.5)),
        half_contractual_rights_retained=evaluate(spec, dict(ref, rights_yield_U=ref['rights_yield_U']/2)),
    )
    sensitivities = []
    for y in (10, 20, 40):
        for retention in (0.4, 0.75, 0.9):
            row = dict(ref, rights_yield_U=y, retention=retention)
            sensitivities.append(dict(yield_U=y, retention=retention,
                operating_audience=audience_threshold(spec, row, 0),
                debt_coverage_1_5_audience=audience_threshold(spec, row, 1.5)))
    legacy_path = 'research/studies/world_economy/results/sports_finance.json'
    legacy = json.loads((ROOT/legacy_path).read_text())
    bond = spec['finance']
    service = cases['reference']['debt_service_U_year']
    launch = dict(capital_U=spec['capital']['installed_U'], working_cash_U=bond['working_cash_U'],
        initial_customer_acquisition_U=cases['reference']['audience']['initial_customer_acquisition_U'],
        debt_reserve_U=service*bond['reserve_years'], debt_fees_U=bond['debt_U']*bond['fee_fraction'])
    launch['total_uses_U'] = sum(launch.values())
    launch['assumed_debt_U'] = bond['debt_U']
    launch['other_funding_required_U'] = launch['total_uses_U']-bond['debt_U']
    supported_debt = max(0,cases['reference']['cash_before_debt_U_year'])/(
        1.5*payment(1,bond['real_rate'],bond['years']))
    ad = spec['advertising_alternative']
    ad_yield = positive(ad['slots_per_hour'])*fraction(ad['fill_fraction'])*positive(ad['CPM_U'])/1000*fraction(ad['owner_share'])
    return dict(schema=SCHEMA, producer=producer('humanities/economics/accounts.py',
        'humanities/economics/sports.py', 'humanities/economics/sports_inputs.json',
        'humanities/economics/sports_sources.json'), evidence=spec['evidence'],
        units='Real U per year unless a field gives another unit; customer-years are annual average paid users.',
        reading_rule=spec['reading_rule'], inputs=spec, cases=cases,
        legacy_reference=dict(revenue_U_year=legacy['operating_scenarios']['reference']['operation']['total_receipts_U'],
            cash_before_debt_U_year=legacy['operating_scenarios']['reference']['operation']['cash_available_for_debt_service_U']),
        operating_sensitivities=alternatives,
        fan_support_cases={name: fan_support(spec, row) for name,row in spec['fan_support_cases'].items()},
        touring_circuit=touring_circuit(spec),
        touring_fare_sensitivity=[dict(return_fare_U=fare, account=touring_circuit(spec,dict(spec['touring_circuit'],return_fare_U=fare))) for fare in spec['touring_circuit']['fare_cases_U']],
        audience_thresholds=dict(operating=audience_threshold(spec, ref, 0),
            debt_service=audience_threshold(spec, ref, 1), coverage_1_5=audience_threshold(spec, ref, 1.5)),
        yield_retention_sensitivity=sensitivities, launch_funding=launch,
        cash_supported_debt_at_coverage_1_5_U=supported_debt,
        scale_requirements=[dict(target_receipts_U_year=t,
            audience_at_reference_yield=max(0,t-ref['sponsor_U_year']-ref['royalties_U_year'])/ref['rights_yield_U'])
            for t in (1e8, 1e9)],
        advertising_alternative=dict(owner_receipt_U_per_viewer_hour=ad_yield,
            hours_for_reference_rights_receipts=cases['reference']['receipts']['media']['U_year']/positive(ad_yield)),
        tax_and_finance='Single annual operating snapshot; tax is an allowance before interest deductions; debt schedule is year one; no subsidy or reserve release is revenue.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'research/runs/cultural_economics/sports.json')
    args = parser.parse_args()
    data = results()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
    print(args.output)


if __name__ == '__main__':
    main()
