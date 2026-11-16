"""Which lines of sale the system serves, at what price and how far, decided from its own costs and its residents.

The demand of a scenario says what customers would buy at a price. This module turns it into orders the system can
fill and is better off filling:

1. It sizes the industry for the residents alone, then for the residents and every line, and cuts every line alike
   in the years the residents could not staff them all. That is the first look.
2. It accounts each line at the first look: the industry is sized again without the line, and the difference is the
   line's cost and the jobs it needs.
3. It prices the lines whose price the system sets. Where buyers answer a price cut with more than a matching rise
   in quantity, the price that brings the seller most is the cost per unit times e/(e-1), e being the elasticity,
   under any ceiling the demand states; where they do not, the seller holds the ceiling, or for visits the price at
   which the published travel model is anchored.
4. It tries a line that costs more than it brings over the programme with a later start, at the first year of each
   later period, and drops it where no start pays. It gives the residents' free hands, year by year, to the lines
   kept, in the order of what a job in each brings.
5. It sizes the industry for those orders, lowers any year the residents still cannot staff, holds the metals
   recovered to the feed those orders mine, and accounts each served line once more against the final orders.

The residents' hands go to the lines of sale first and to making the working parts of machinery after, so the
industry model staffs that manufacture from what the final orders leave free. Within one account the staffing of
that manufacture is held as the reference orders set it, so that a difference is the line's alone.

Costs are compared at one real discount rate. A line's cost is the whole difference it makes to the industry, with
the equipment, the freight and the propellant it calls for and the financing's allowance on capital, so a line is
charged for what the system builds for it and credited with nothing it does not sell. Where the residents cannot
staff all the working parts the industry is able to make, a line is also charged for the parts its jobs would have
made, since those are then bought from Earth.
"""
from __future__ import annotations
import math

from provisioning.vacuum_industry import industry, model
from . import demand, programme as definition

RATE = 0.03                      # the real discount rate at which lines are compared
MARGIN = 0.97                    # the share of the workforce the allocation fills
PASSES = 6


def by_year(physical):
    return {row['year']: row for row in physical['annual']}


def yearly_cost(physical, allowance=0.0):
    """What the industry spends in each year, in dollars of 2041, with the founding deployment as one amount.

    The financing adds an allowance to capital, as a fraction of the machinery bought and built and of the founding
    deployment beyond its provision for resources; the same allowance is counted here.
    """
    costs = {}
    for y, row in by_year(physical).items():
        cash = row['cash']
        costs[y] = (cash['capital_total_USD2041'] + cash['operating_total_USD2041']
                    + allowance * (cash['capital_external_USD2041'] + cash['capital_local_USD2041']))
    seed = physical['seed']
    founding = seed['cost']['total_USD2041']
    return costs, founding + allowance * max(0.0, founding - seed['total_resource_provision_USD2041'])


def jobs(physical):
    return {y: row['technical_resident_jobs'] for y, row in by_year(physical).items()}


def feed(physical):
    return {y: row['physical']['source_feed_kg_year'] for y, row in by_year(physical).items()}


def period_of(periods, year):
    for index, (first, last) in enumerate(periods):
        if first <= year <= last:
            return index
    raise ValueError('Year outside the reporting periods')


class Decision:
    """The orders of one scenario case, with the industry that fills them."""

    def __init__(self, programme, ramp, inputs, scenario, coefficients, recovery, measures, periods, hotel, log=None, pricing='seller',
                 allowance=0.0):
        self.programme, self.ramp, self.inputs, self.scenario = programme, ramp, inputs, scenario
        self.coefficients, self.recovery, self.measures, self.periods = coefficients, recovery, measures, periods
        self.hotel = hotel               # what a visitor's bed costs in staff each year, beside the industry's account
        self.pricing = pricing
        self.allowance = allowance       # the financing's allowance on capital
        self.log = log or (lambda *a: None)
        self.solves = 0
        calendar = programme['calendar']
        self.years = range(calendar['first_productive_year'], calendar['last_year'] + 1)
        self.founding_years = calendar['first_productive_year']

    def orders(self, prices=None, scale=None, mined=None):
        return demand.book(self.programme, self.ramp, self.inputs, self.scenario, self.coefficients, prices, scale, mined)

    def industry(self, book, staffed=None):
        """Size the industry for a book. A stated staffing of local manufacture is held; none lets the model find it."""
        data = self.coefficients
        if staffed is not None:
            production = data['economy']['production']
            working = dict(production.get('local_working_parts', {}), staffed=staffed)
            data = dict(data, economy=dict(data['economy'], production=dict(production, local_working_parts=working)))
        self.solves += 1
        return industry.build(data, self.recovery, book, self.measures)

    @staticmethod
    def staffing(loops):
        """The staffing of local manufacture a sizing settled on, as a schedule that can be held."""
        found = loops.get('local_working_parts', {}).get('staffed_by_order_year')
        return {'by_year': [[int(y), v] for y, v in sorted(found.items(), key=lambda kv: int(kv[0]))]} if found else None

    def displaced(self, staffed):
        """What a job taken from making working parts costs in each year: the parts it would have made, bought instead.

        Where the residents cannot staff all the working parts the industry is able to make, a job given to a line of
        sale is a job that manufacture loses. A job-year makes parts worth its wage divided by the labour's share
        of their import price, so the loss is that value less the wage.
        """
        working = model.local_working_parts(self.coefficients)
        wage = self.coefficients['industrial']['economics']['loaded_staff_USD2026_year'] * self.coefficients['normalization']['USD2026_to_USD2041']
        out = {}
        for y in self.years:
            fraction = model.at(working['labour_fraction_of_import_price'], y)
            short = staffed is not None and model.at(working['share'], y) > 0 and model.at(staffed, y) < 1 - 1e-9
            out[y] = wage * model.labour_index(self.coefficients, y) * (1 / fraction - 1) if short and fraction > 0 else 0.0
        return out

    def none(self, lines=demand.LINES):
        return {line: {y: 0.0 for y in self.years} for line in lines}

    def without(self, scale, line):
        return {**(scale or {}), line: {y: 0.0 for y in self.years}}

    def discount(self, year):
        return (1 + RATE) ** -year

    def account(self, prices, scale, mined, book, physical, lines, staffed):
        """Each line's receipts, cost and jobs against the same orders without it."""
        rows = {row['year']: row for row in book['annual']}
        cost, founding = yearly_cost(physical, self.allowance)
        hands = jobs(physical)
        lost = self.displaced(staffed)
        out = {}
        for line in lines:
            if not any(rows[y]['lines'][line]['quantity'] for y in self.years):
                continue
            less_book = self.orders(prices, self.without(scale, line), mined)
            less, _ = self.industry(less_book, staffed)
            less_cost, less_founding = yearly_cost(less, self.allowance)
            less_hands = jobs(less)
            less_rows = {row['year']: row for row in less_book['annual']}
            years = {}
            for y in self.years:
                detail = rows[y]['lines'][line]
                beds = (rows[y]['visitor_beds'] - less_rows[y]['visitor_beds']) * self.hotel
                years[y] = dict(receipts=detail['receipts_USD2041'], quantity=detail['quantity'],
                                cost=cost[y] - less_cost[y] + detail['own_cost_USD2041'] + beds + max(0.0, hands[y] - less_hands[y]) * lost[y],
                                jobs=hands[y] - less_hands[y])
            out[line] = dict(years=years, founding=founding - less_founding, unit=rows[self.years[0]]['lines'][line]['unit'],
                             without=(less_cost, less_founding, less_hands))
        return out

    def summary(self, account):
        """Each line by period, and its worth over the programme at the comparison rate."""
        out = {}
        for line, record in account.items():
            periods = []
            for index, (first, last) in enumerate(self.periods):
                span = [record['years'][y] for y in self.years if first <= y <= last]
                if not span:
                    continue
                quantity = math.fsum(v['quantity'] for v in span)
                cost = math.fsum(v['cost'] for v in span) + (record['founding'] if index == 1 else 0.0)
                periods.append(dict(first_year=first, last_year=last, receipts_USD2041=math.fsum(v['receipts'] for v in span),
                                    cost_USD2041=cost, quantity=quantity, unit_cost_USD2041=cost / quantity if quantity else None,
                                    mean_jobs=math.fsum(v['jobs'] for v in span) / len(span)))
            worth = math.fsum((v['receipts'] - v['cost']) * self.discount(y) for y, v in record['years'].items())
            worth -= record['founding'] * self.discount(self.founding_years / 2)
            hands = math.fsum(v['jobs'] * self.discount(y) for y, v in record['years'].items())
            out[line] = dict(unit=record['unit'], periods=periods, worth_USD2041=worth, worth_per_job_USD2041=worth / hands if hands > 0 else None)
        return out

    def own_prices(self, summary):
        """The system's price for each line it prices, by year, from the line's cost per unit in each period."""
        prices = {}
        for line in demand.PRICED:
            if line not in summary:
                continue
            unit = {(p['first_year'], p['last_year']): p['unit_cost_USD2041'] for p in summary[line]['periods']}
            known = [v for v in unit.values() if v and v > 0]
            if not known:
                continue
            chosen = {}
            for y in self.years:
                cost = unit.get(tuple(self.periods[period_of(self.periods, y)]))
                price = demand.seller_price(self.inputs, self.scenario, line, cost if cost and cost > 0 else known[0],
                                            self.coefficients['normalization']['USD2026_to_USD2041'], self.pricing)
                if price is not None:
                    chosen[y] = price
            if chosen:
                prices[line] = chosen
        return prices

    def worth_at(self, summary, line, book):
        """A priced line's worth at the system's price, from its cost per unit at the first look."""
        unit = {(p['first_year'], p['last_year']): p['unit_cost_USD2041'] for p in summary[line]['periods']}
        known = [v for v in unit.values() if v and v > 0]
        total = 0.0
        for row in book['annual']:
            cost = unit.get(tuple(self.periods[period_of(self.periods, row['year'])]))
            cost = cost if cost and cost > 0 else (known[0] if known else 0.0)
            detail = row['lines'][line]
            total += (detail['price_USD2041'] - cost) * detail['quantity'] * self.discount(row['year'])
        return total

    def late_start(self, line, scale, mined, staffed, record):
        """Whether a line that loses when served from the opening pays if it starts with a later period, and from when.

        The industry is sized with the line left out until each later period's first year, against the same orders
        without the line. Returns the best worth and its first year, or nothing where no start pays.
        """
        less_cost, less_founding, less_hands = record['without']
        lost = self.displaced(staffed)
        best = (0.0, None)
        for start in [first for first, _ in self.periods[2:] if first <= self.programme['calendar']['completion_year']]:
            later = {**scale, line: {y: (scale[line][y] if y >= start else 0.0) for y in self.years}}
            book = self.orders(None, later, mined)
            physical, _ = self.industry(book, staffed)
            cost, founding = yearly_cost(physical, self.allowance)
            hands = jobs(physical)
            worth = -(founding - less_founding) * self.discount(self.founding_years / 2)
            for row in book['annual']:
                y, detail = row['year'], row['lines'][line]
                beds = row['visitor_beds'] * self.hotel if line == 'visits' else 0.0
                taken = max(0.0, hands[y] - less_hands[y]) * lost[y]
                worth += (detail['receipts_USD2041'] - detail['own_cost_USD2041'] - beds - taken - (cost[y] - less_cost[y])) * self.discount(y)
            if worth > best[0]:
                best = (worth, start)
        return best

    def allocate(self, summary, account, first_book, book, base_jobs, starts=None):
        """Free hands to the lines worth serving, in the order of what a job in each brings."""
        starts = starts or {}
        first = {row['year']: row for row in first_book['annual']}
        now = {row['year']: row for row in book['annual']}
        opening = {line: starts[line][1] if line in starts else self.years[0] for line in summary}
        worth = {line: (starts[line][0] if line in starts else self.worth_at(summary, line, book) if line in demand.PRICED
                        else summary[line]['worth_USD2041']) for line in summary}
        hands = {line: math.fsum(max(0.0, v['jobs']) * self.discount(y) for y, v in account[line]['years'].items() if y >= opening[line])
                 for line in summary}
        kept = [line for line in summary if worth[line] > 0]
        order = sorted(kept, key=lambda line: worth[line] / hands[line] if hands[line] > 0 else math.inf, reverse=True)
        scale = self.none()
        for y in self.years:
            free = MARGIN * model.WORKFORCE_SHARE * now[y]['population'] - base_jobs[y]
            for line in order:
                if y < opening[line]:
                    continue
                before = first[y]['lines'][line]['quantity']
                per_unit = max(0.0, account[line]['years'][y]['jobs']) / before if before else 0.0
                wanted = per_unit * now[y]['lines'][line]['quantity']
                share = 1.0 if wanted <= 0 else max(0.0, min(1.0, free / wanted))
                scale[line][y] = share
                free -= share * wanted
        return scale, order, [line for line in summary if line not in kept]

    def decide(self):
        residents = {y: definition.population(self.programme, self.ramp, y) for y in self.years}
        room = {y: MARGIN * model.WORKFORCE_SHARE * residents[y] for y in self.years}
        # The residents alone, and every line at its first price, with every working part bought: the hands each needs.
        base, _ = self.industry(self.orders(scale=self.none()), 0.0)
        base_jobs = jobs(base)
        mined = feed(base)
        every, _ = self.industry(self.orders(mined=mined), 0.0)
        wanted = jobs(every)
        alike = {y: max(0.0, min(1.0, (room[y] - base_jobs[y]) / max(1.0, wanted[y] - base_jobs[y]))) for y in self.years}
        first_scale = {line: dict(alike) for line in demand.LINES}
        # The first look: every line cut alike to what the residents can staff.
        first_book = self.orders(None, first_scale, mined)
        first, first_loops = self.industry(first_book)
        held = self.staffing(first_loops)
        first_account = self.account(None, first_scale, mined, first_book, first, demand.LINES, held)
        first_summary = self.summary(first_account)
        prices = self.own_prices(first_summary)
        priced_book = self.orders(prices, None, feed(first))
        # A line that loses when served from the opening may pay if it starts later.
        starts = {}
        for line in first_summary:
            worth = self.worth_at(first_summary, line, priced_book) if line in demand.PRICED else first_summary[line]['worth_USD2041']
            if worth <= 0:
                late = self.late_start(line, first_scale, mined, held, first_account[line])
                if late[1] is not None:
                    starts[line] = late
        mined = dict(feed(first))
        scale, order, dropped = self.allocate(first_summary, first_account, first_book, priced_book, base_jobs, starts)
        self.log('lines in order of worth per job:', order, 'dropped:', dropped, 'starting late:', {line: late[1] for line, late in starts.items()})
        for _ in range(PASSES):
            book = self.orders(prices, scale, mined)
            physical, loops = self.industry(book)
            use = {y: row['capacity_summary']['resident_workforce'] for y, row in by_year(physical).items()}
            over = [y for y in self.years if use[y] > 1.0]
            # The metals recovered cannot come from more feed than these orders mine.
            now = feed(physical)
            short = [row['year'] for row in book['annual'] if row['lines']['metals']['treated_feed_kg'] > now[row['year']] * (1 + 1e-9)]
            if not over and not short:
                break
            for y in short:
                mined[y] = now[y]
            self.log('years over the workforce:', len(over), 'worst', round(max(use.values()), 3), '| years short of feed for the metals:', len(short))
            hands = jobs(physical)
            for y in over:
                # The hands the lines hold beyond the residents' own work, cut to those that are free.
                held = max(1.0, hands[y] - base_jobs[y])
                cut = MARGIN * min(1.0, max(0.0, room[y] - base_jobs[y]) / held)
                for line in scale:
                    scale[line][y] *= cut
        served = [line for line in order if any(scale[line][y] > 0 for y in self.years)]
        final_account = self.account(prices, scale, mined, book, physical, served, self.staffing(loops))
        final = self.summary(final_account)
        return dict(book=book, physical=physical, loops=loops, prices=prices, scale=scale, order=order, dropped=dropped,
                    starts={line: late[1] for line, late in starts.items()}, late=starts,
                    first=first_summary, first_account=first_account, first_book=first_book, final=final, base_jobs=base_jobs,
                    solves=self.solves)
