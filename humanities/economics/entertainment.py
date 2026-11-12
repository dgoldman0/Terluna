"""Cultural-export portfolio with shared audiences, facilities and visiting artists."""
from __future__ import annotations

import argparse
from copy import deepcopy
import json
from pathlib import Path

from .accounts import (HERE, ROOT, account, count, fraction, nonnegative,
                       positive, producer, read_bound_inputs)

INPUT = HERE/'entertainment_inputs.json'
SPORTS_PRODUCT = 'humanities/economics/results/sports.json'


def audience(groups):
    """Disjoint cohorts carry possibly overlapping genre interests."""
    unique = 0
    by_genre = {}
    for row in groups:
        people = nonnegative(row['fans'])
        genres = row['genres']
        if not genres or len(set(genres)) != len(genres):
            raise ValueError('Genre memberships must be nonempty and unique')
        unique += people
        for genre in genres:
            by_genre[genre] = by_genre.get(genre,0)+people
    return dict(unique_fans=unique, fans_by_genre=by_genre,
                sum_genre_followings=sum(by_genre.values()))


def evaluate(spec, sport, case, fixed_payroll_U=0):
    src = sport['inputs']['touring_circuit']
    old = sport['touring_circuit']
    extra = spec['production'] if case['extra_programme'] else spec['empty_production']
    fans = audience(case['audience_groups'])
    n = fans['unique_fans']
    paid_share = fraction(case['paid_share'])
    member_share = fraction(spec['membership_share'])
    if member_share > paid_share:
        raise ValueError('Membership cohort exceeds paid cohort')
    audience_factor = fraction(case.get('earned_receipt_factor',1))
    paid = n*paid_share
    membership = n*member_share
    receipts = {
        'portfolio_viewing': dict(U_year=paid*nonnegative(case['owner_viewing_yield_U'])*audience_factor, export_fraction=1),
        'additional_membership_services': dict(U_year=membership*spec['membership_yield_U']*audience_factor, export_fraction=1),
        'game_merchandise_and_archive_royalties': dict(U_year=n*spec['royalty_U_fan']*audience_factor, export_fraction=1),
        'sponsor_services': dict(U_year=n*spec['sponsor_U_fan']*audience_factor, export_fraction=1),
        'separate_production_commissions': dict(U_year=extra['commission_receipts_U_year']*audience_factor,export_fraction=1),
        'local_admissions': deepcopy(old['receipts']['local_admissions']),
    }
    remote = sum(r['U_year'] for r in receipts.values() if r['export_fraction'])
    purchases = {name: deepcopy(row) for name,row in old['purchases'].items()
                 if name not in ('commercial_services','audience_acquisition')}
    purchases['commercial_services'] = dict(U_year=remote*spec['variable_cost_fraction'],import_fraction=0.6)
    purchases['audience_acquisition'] = dict(U_year=paid*(1-fraction(case['retention']))*spec['acquisition_U_customer'],import_fraction=0.9)
    extra_residents = nonnegative(extra['resident_workers'])
    visits = nonnegative(extra['visiting_workers'])*nonnegative(extra['visits_year'])
    nights = visits*nonnegative(extra['nights_visit'])
    purchases['additional_resident_host'] = dict(U_year=extra_residents*src['host_U_resident_year'],import_fraction=0)
    purchases['additional_guest_host'] = dict(U_year=nights*src['host_U_guest_night'],import_fraction=0)
    purchases['additional_return_travel'] = dict(U_year=visits*src['return_fare_U'],import_fraction=1)
    purchases['studio_operations'] = dict(U_year=nonnegative(extra['studio_operations_U_year']),import_fraction=0.4)
    capture_cost = 0
    studio_hours = 0
    arena_hours = 0
    work_hours = 0
    for row in extra['programmes']:
        units = count(row['units_year'])
        capture_cost += units*nonnegative(row['purchased_capture_U_unit'])
        studio_hours += units*nonnegative(row['studio_hours_unit'])
        arena_hours += units*nonnegative(row['arena_hours_unit'])
        work_hours += units*positive(row['crew'])*positive(row['hours_per_crew_member_unit'])
    purchases['additional_capture_services'] = dict(U_year=capture_cost,import_fraction=0.5)
    resident_base = old['capacity']['resident_athletes']*src['resident_athlete_pay_U_year'] + src['resident_support_workers']*src['resident_support_pay_U_year'] + extra_residents*spec['resident_pay_U_year']
    foreign_base = src['visits_per_person_year']*(old['capacity']['visiting_athletes']*src['athlete_fee_U_visit']+src['visiting_support_workers']*src['support_fee_U_visit']) + visits*spec['visitor_pay_U_visit']
    minimum_wages = resident_base+foreign_base
    wages = max(minimum_wages,remote*spec['payroll_revenue_fraction'],nonnegative(fixed_payroll_U))
    foreign_wages = wages*foreign_base/positive(minimum_wages)
    capital = dict(sport['inputs']['capital'], installed_U=old['installed_capital_U']+nonnegative(extra['studio_capital_U']))
    out = account(receipts,wages,purchases,capital,sport['inputs']['finance'],
        sport['inputs']['tax_fraction'],sport['inputs']['supplier_import_fraction'],foreign_wages)
    cap = old['capacity']
    total_hours = cap['occupied_arena_hours']+arena_hours
    total_nights = cap['guest_nights_year']+nights
    total_arrivals = cap['visitor_arrivals_year']+visits
    peak = cap['peak_guests']+extra['visiting_workers']
    residents = cap['resident_workers']+extra_residents
    studio_limit = nonnegative(extra['studio_available_hours_year'])
    work_limit = extra_residents*spec['work_hours_per_resident_year']
    if (total_hours > cap['available_arena_hours'] or total_nights > cap['guest_night_allowance_year']
            or total_arrivals > src['arrival_capacity_year'] or peak > src['guest_beds']
            or residents > src['host_residents'] or studio_hours > studio_limit or work_hours > work_limit):
        raise ValueError('Portfolio exceeds shared production capacity')
    out['audience'] = dict(**fans, paid_viewer_equivalents=paid,membership_buyers=membership)
    out['capacity'] = dict(resident_workers=residents, host_residents=src['host_residents'],
        arrivals_year=total_arrivals, arrival_capacity_year=src['arrival_capacity_year'],
        guest_nights_year=total_nights, guest_night_allowance_year=cap['guest_night_allowance_year'],
        peak_guests=peak, guest_beds=src['guest_beds'],
        arena_hours_year=total_hours, arena_hours_available=cap['available_arena_hours'],
        studio_hours_year=studio_hours, studio_hours_available=studio_limit,
        additional_crew_hours=work_hours, additional_crew_hours_available=work_limit)
    out['increment_against_sport'] = {k: out[k]-old[k] for k in (
        'exports_U_year','gross_value_added_U_year','cash_before_debt_U_year',
        'external_cash_after_financing_U_year','installed_capital_U')}
    out['patron_transfers_U_year'] = n*spec['patron_share']*spec['patron_gift_U']
    out['initial_paid_audience_acquisition_U'] = paid*spec['acquisition_U_customer']
    return out


def results():
    spec = read_bound_inputs(INPUT,'terluna.humanities.entertainment-inputs/1')
    sport = json.loads((ROOT/SPORTS_PRODUCT).read_text())
    cases = {name:evaluate(spec,sport,row) for name,row in spec['cases'].items()}
    outage = dict(spec['cases']['expanded'], earned_receipt_factor=0.5)
    cases['shared_revenue_shock'] = evaluate(spec,sport,outage,
        fixed_payroll_U=cases['expanded']['wages_U_year'])
    return dict(schema='terluna.humanities.entertainment-economics/1',
        producer=producer('humanities/economics/accounts.py','humanities/economics/entertainment.py',
            'humanities/economics/entertainment_inputs.json','humanities/economics/entertainment_sources.json'),
        evidence=spec['evidence'],reading_rule=spec['reading_rule'],
        units='Real U per year except capital, initial audience acquisition, people and physical capacities.',
        inputs=spec,cases=cases,
        sports_input_sha256=spec['products'][SPORTS_PRODUCT]['sha256'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'research/runs/cultural_economics/entertainment.json')
    args = parser.parse_args()
    data = results()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    print(args.output)


if __name__ == '__main__':
    main()
