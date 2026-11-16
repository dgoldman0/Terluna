"""How many residents the vacuum buildout carries when their number follows from its own production and financing.

usage: python3 -m research.studies.vacuum_economy.residents [--scenario NAME ...] [--jobs N] [--output PATH]   (in the repository root)

The programme's two ramps state how many residents are housed in each year, and the run tests whether the system pays
for them. The solved ramps of programme.json turn the question round. From the opening population in the first
productive year, the settlement grows period by period, in order, at the highest steady rate the system sustains. A
rate is sustained where, with no growth after its period,

  money   the opening funding is no more than on the constant-growth ramp, and the system needs no equity after it,
          under the study's financing rules with credit in U and with its surplus kept for growth;
  hands   the residents' own work, for themselves, for the newcomers' housing and supply and for the industry these
          need, fits within the workforce, and the lines of sale take the hands that are left;
  income  the receipts of years 201 to 220 pay the running costs and renewals of the population reached, so that the
          last population is carried by what the system earns and not by what it saved.

The housing of the first newcomers is ordered before opening and paid from the opening funding, so the first period
covers the years those orders serve, and the opening funding is held to what the constant-growth ramp needs: faster
early growth must be paid from the system's own receipts and credit. Growth after a period's rate is chosen only raises
costs that later periods must also meet, so each period's search starts from a schedule already sustained. The ramp `resident_led_beyond_plan` is this schedule; `resident_led` is the
same schedule held at the planned population from the year it reaches it, checked once more on its own.

The lines of sale, their prices and their first years are those mix.py chooses on the constant-growth ramp. Once a
schedule is found, the choice is made again on that schedule and the schedule is found once more with it. A trial
sizes the industry twice: for the residents alone with every working part bought, which counts the hands their own
work takes, and for the residents with the lines, the staffing of local manufacture held at the last choice's and
lowered in any year where it would take hands that are not free.

The water each schedule draws is set against the contained water of the wider near-Earth screen of the source survey,
at its 5th, 50th and 95th percentiles. It is reported and does not limit growth: the lines of sale draw on the same
water, and the study does not yet divide the pool between them and the residents.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import multiprocessing
from pathlib import Path
import pickle
import time

from provisioning.vacuum_industry import industry, model
from research.studies.enclosed_returns import finance, revised
from . import inputs, mix, programme as definition
from .inputs import HERE, ROOT, digest

SCHEMA = 'terluna.research.vacuum-economy-residents/1'
PRODUCT = HERE / 'results/residents.json'
RUNS = ROOT / 'research/runs/vacuum_economy_residents'
RAMP, LIMITED = 'resident_led_beyond_plan', 'resident_led'
REFERENCE = 'constant_growth'
SCENARIOS = ('low', 'central', 'high')
# The periods of steady growth, from the first year after opening to the completion year. The first holds the years
# whose newcomers' housing is ordered before opening.
PERIODS = ((9, 12), (13, 18), (19, 25), (26, 37), (38, 50), (51, 65), (66, 80), (81, 100), (101, 125), (126, 150), (151, 175),
           (176, 200))
# The share by which a trial's opening funding may exceed the cap's basis: the rounding of two sizings of the same
# opening differs by about 1e-4, and a point of growth in the first period adds about 8e-4.
OPENING_TOLERANCE = 1e-3
INCOME_TOLERANCE = 1e-6        # the share of the opening funding by which holdings may fall over years 201 to 220
ESCALATION = 4                  # passes, after the proportionate ones, that halve the lines of a year still over the workforce
HIGHEST = 0.30           # the highest yearly growth rate tried
RESOLUTION = 0.0005      # a period's search ends when its highest sustained and lowest unsustained rates are this close,
RELATIVE = 0.05          # or this share of each other,
TRIALS = 8               # or after this many trials
ROUNDS = 1               # schedules found in turn, each on the choice of sales made on the one before
TAIL = (201, 220)
EQUITY_TOLERANCE_USD2041 = 1e6


def dependencies():
    """Every file a solution reads: the coupled study's inputs and this module."""
    return inputs.dependencies((HERE / 'residents.py',))


def schedule(programme, rates, limit=None):
    """[[year, residents], ...] from first production to the last year, for one growth rate in each period.

    Periods without a rate have none; the population is held after the last period, and at the limit from the year it
    reaches it."""
    calendar = programme['calendar']
    first, last = calendar['first_productive_year'], calendar['last_year']
    if PERIODS[0][0] != first + 1 or PERIODS[-1][1] != calendar['completion_year'] or any(a[1] + 1 != b[0] for a, b in zip(PERIODS, PERIODS[1:])):
        raise ValueError('The periods of growth must run without gaps from the first year after opening to completion')
    if len(rates) > len(PERIODS) or any(r < 0 for r in rates):
        raise ValueError('One rate of growth, none negative, for each period at most')
    rate = {y: r for (a, b), r in zip(PERIODS, rates) for y in range(a, b + 1)}
    people, out = programme['population']['opening'], []
    for y in range(first, last + 1):
        if y > first:
            people *= 1 + rate.get(y, 0.0)
            if limit is not None:
                people = min(people, limit)
        out.append([y, people])
    return out


def money(path, financed, rules):
    """The least cash the system holds above the buffer the financing rules require, over all years, and the equity
    it is given after its opening funding."""
    start = path['annual'][0]['year']
    least = min(f['pre_topup_cash'] - r['operating_cost'] * max(rules['liquidity_buffer_years'], rules['initial_runway_years'] - (f['year'] - start))
                for f, r in zip(financed['annual'], path['annual']))
    return least, math.fsum(f['contribution'] for f in financed['annual'])


def holding(financed, year):
    """What the system holds at the end of a year: cash, reserves and custody, less its debt."""
    f = next(f for f in financed['annual'] if f['year'] == year)
    return f['cash'] + f['reserves'] + f['owner_custody'] - f['debt']


class Solver:
    """The resident-led schedule of one scenario of demand."""

    def __init__(self, scenario, log=None, state=None):
        s = inputs.sources()
        self.s, self.scenario = s, scenario
        self.log = log or (lambda *a: print(scenario, *a, flush=True))
        self.state = state
        self.case = dict(name='%s:%s' % (scenario, RAMP), scenario=scenario, ramp=RAMP, recovery='combined')
        self.coefficients = inputs.engineering(self.case)
        self.money = dict(s['money'], horizons=s['programme']['horizons'])
        last = s['programme']['calendar']['last_year']
        # While a schedule is sought the surplus is kept for growth: nothing is paid out before the last year.
        self.keeping = dict(self.money, distribution_policy='retain_to_completion', retain_through_year=last)
        self.hotel = s['money']['hotel_staff_per_bed'] * s['money']['hotel_staff_USD2041_year']
        self.planned = s['programme']['population']['planned']
        self.reference_opening = None   # the constant-growth ramp's opening funding, set by the first choice of sales
        self.cap = None                 # the opening funding no schedule may exceed, set before the first period's search
        self.trials = 0

    def decision(self, programme, ramp):
        return mix.Decision(programme, ramp, self.s['scenarios'], self.scenario, self.coefficients, self.s['recovery'], 'combined',
                            self.s['periods'], self.hotel, allowance=self.s['money']['capital_cost_allowance_fraction'])

    def programme(self, residents, ramp=RAMP):
        return definition.solved(self.s['programme'], ramp, residents)

    def choose(self, residents=None):
        """The system's choice of sales on the reference ramp, or on a solved schedule, with what a trial reads of it.

        On the reference ramp the choice also gives the opening funding that no solved schedule may exceed."""
        label = 'choice on ' + (REFERENCE if residents is None else 'schedule %s' % hashlib.sha256(json.dumps(residents).encode()).hexdigest()[:16])
        choice = self.state.kept_choice(label) if self.state is not None else None
        if choice is None:
            start = time.monotonic()
            programme = self.s['programme'] if residents is None else self.programme(residents)
            decided = self.decision(programme, REFERENCE if residents is None else RAMP).decide()
            if not revised.valid(decided['physical']):
                raise ValueError('The choice of sales is not physically resolved for ' + self.scenario)
            choice = {key: decided[key] for key in ('prices', 'first', 'first_account', 'first_book', 'late', 'order', 'dropped', 'starts')}
            choice['opening_funding_USD2041'] = finance.plan(decided['book'], decided['physical'], self.money)['seed']
            self.log(label, 'made in', round(time.monotonic() - start), 's:', choice['order'], 'dropped', choice['dropped'], 'late', choice['starts'])
            if self.state is not None:
                self.state.keep_choice(label, choice)
        if residents is None:
            # The constant-growth ramp's opening funding, the unit in which the margins of a trial are stated.
            self.reference_opening = choice['opening_funding_USD2041']
        return choice

    def fill(self, d, programme, ramp, choice):
        """The orders and the industry for one schedule: hands to the lines as the choice ranks them, then to parts.

        The industry is first sized for the residents alone with every working part bought, as mix.py sizes it, which
        gives the hands their own work takes; its ledger gives the hands that making the parts it is able to make would
        take. The lines take their hands as the choice's accounts give them per unit sold, and the parts are staffed
        from what is left; a year that still takes more hands than the workforce holds has its manufacture, and then its
        lines, lowered."""
        data = d.coefficients
        production = data['economy']['production']
        bought = dict(data, economy=dict(data['economy'], production=dict(production, local_working_parts=dict(
            production.get('local_working_parts', {}), staffed=0.0))))
        d.solves += 1
        base, _, ledger = industry.solve_once(bought, d.recovery, d.orders(scale=d.none()), d.measures, True)
        wage = data['industrial']['economics']['loaded_staff_USD2026_year']
        parts = {y: ledger.flows[y].get('local_working_labour_if_staffed', 0.0) / wage for y in d.years}
        base_jobs = mix.jobs(base)
        mined = dict(mix.feed(base))
        priced = d.orders(choice['prices'], None, mined)
        scale, _, _ = d.allocate(choice['first'], choice['first_account'], choice['first_book'], priced, base_jobs, choice['late'])
        residents = {y: definition.population(programme, ramp, y) for y in d.years}
        room = {y: mix.MARGIN * model.WORKFORCE_SHARE * residents[y] for y in d.years}
        first = {row['year']: row for row in choice['first_book']['annual']}
        now = {row['year']: row for row in priced['annual']}

        def line_jobs(y):
            total = 0.0
            for line in scale:
                before = first[y]['lines'][line]['quantity']
                if line in choice['first_account'] and before:
                    total += max(0.0, choice['first_account'][line]['years'][y]['jobs']) / before * now[y]['lines'][line]['quantity'] * scale[line][y]
            return total
        reach = {y: 1.0 if parts[y] <= 0 else max(0.0, min(1.0, (room[y] - base_jobs[y] - line_jobs(y)) / parts[y])) for y in d.years}
        # An order is worked in its own year and the next, so both years' hands must reach.
        staffing = {y: min(reach[y], reach.get(y + 1, reach[y])) for y in d.years}
        # Whether the residents' own work, for themselves and their newcomers, fits within the workforce in every year.
        own = min((room[y] - base_jobs[y]) / room[y] for y in d.years)
        for attempt in range(mix.PASSES + ESCALATION):
            book = d.orders(choice['prices'], scale, mined)
            physical, loops = d.industry(book, {'by_year': [[y, staffing[y]] for y in d.years]})
            rows = mix.by_year(physical)
            over = [y for y in d.years if rows[y]['capacity_summary']['resident_workforce'] > 1 + 1e-9]
            if not over:
                break
            hands = mix.jobs(physical)
            for y in over:
                made = rows[y]['working_parts_jobs']
                others = hands[y] - made
                if attempt >= mix.PASSES:
                    # Where the proportionate cuts have not fitted a year, its manufacture stops and its lines halve.
                    for k in (y - 1, y):
                        if k in staffing:
                            staffing[k] = 0.0
                    for line in scale:
                        scale[line][y] *= 0.5
                elif made > 0 and others < room[y]:
                    # Fewer parts made: an order is worked in its own year and the next, so both orders are lowered.
                    factor = (room[y] - others) / made
                    for k in (y - 1, y):
                        if k in staffing:
                            staffing[k] *= factor
                else:
                    for k in (y - 1, y):
                        if k in staffing:
                            staffing[k] = 0.0
                    held = max(1.0, others - base_jobs[y])
                    cut = mix.MARGIN * min(1.0, max(0.0, room[y] - base_jobs[y]) / held)
                    for line in scale:
                        scale[line][y] *= cut
        use = max(row['capacity_summary']['resident_workforce'] for row in physical['annual'])
        return book, physical, loops, use, own

    def trial(self, rates, choice, limit=None, ramp=RAMP):
        """Whether the schedule of these rates is sustained: the conditions that fail, and the margin of each.

        The margins are shares of the constant-growth ramp's opening funding (money: the least cash above the buffer,
        or the equity needed; opening: what the cap exceeds the opening funding by; income: the change in what the
        system holds over years 201 to 220) or of the workforce (hands: the room left by the residents' own work, or
        the excess of jobs over the workforce)."""
        start = time.monotonic()
        self.trials += 1
        residents = schedule(self.s['programme'], rates, limit)
        programme = self.programme(residents, ramp)
        d = self.decision(programme, ramp)
        book, physical, loops, use, own = self.fill(d, programme, ramp, choice)
        out = dict(rates=list(rates), residents_at={str(h): definition.population(programme, ramp, h) for h in (100, 200, 220)},
                   workforce_use=use, own_work_room=own)
        if own < 0 or not revised.valid(physical) or use > 1 + 1e-6:
            out.update(sustained=False, failed=['hands'], binding='hands', margin=min(own, 1 - use, -1e-9))
        else:
            path = finance.plan(book, physical, self.money)
            kept = finance.finance(path, self.keeping, 'equity')
            credit = finance.finance(path, self.keeping, 'U_credit', kept)
            least, equity = money(path, credit, self.keeping)
            tail = holding(credit, TAIL[1]) - holding(credit, TAIL[0] - 1)
            scale = self.reference_opening
            checks = dict(money=(equity <= EQUITY_TOLERANCE_USD2041, (least if equity <= EQUITY_TOLERANCE_USD2041 else -equity) / scale),
                          income=(tail >= -INCOME_TOLERANCE * scale, tail / scale), hands=(True, own))
            if self.cap is not None:
                checks['opening'] = (path['seed'] <= self.cap, (self.cap - path['seed']) / scale)
            failed = [k for k, (ok, _) in checks.items() if not ok]
            margins = {k: m for k, (_, m) in checks.items()}
            binding = min(failed or margins, key=margins.get)
            rows = {r['year']: r for r in path['annual']}
            out.update(sustained=not failed, failed=failed, binding=binding, margin=margins[binding], margins=margins,
                       equity_after_opening_USD2041=equity, least_cash_above_buffer_USD2041=least, holding_change_201_220_USD2041=tail,
                       opening_funding_USD2041=path['seed'], receipts_USD2041_year={str(h): rows[h]['gross'] for h in (100, 200, 220)},
                       water=inputs.water_drawn(physical, loops, self.coefficients, (100, 200, 220)))
        out['seconds'] = round(time.monotonic() - start, 1)
        self.log('trial %d' % self.trials, ' '.join('%.4f' % r for r in rates), 'sustained' if out['sustained'] else 'not sustained',
                 'binding', out['binding'], 'margin %.4g' % out['margin'], 'residents at 200 %.0f' % out['residents_at']['200'], out['seconds'], 's')
        return out

    def period(self, rates, index, choice, guess, step):
        """The highest sustained rate for one period, given the rates before it: a bracket, then false position.

        The schedule with no growth in this period is the one the previous period's search sustained, so the bracket
        starts there and widens from the guess by the given factor."""
        tried = []
        lo, at_lo, hi, at_hi = 0.0, None, None, None
        r = min(HIGHEST, guess) if guess > 0 else 0.02
        while hi is None and len(tried) < TRIALS:
            result = self.trial(rates + [r], choice)
            tried.append(result)
            if result['sustained']:
                lo, at_lo = r, result
                if r >= HIGHEST:
                    break
                r = min(HIGHEST, r * step)
            else:
                hi, at_hi = r, result
        # side counts the moves of one end in a row: negative for the sustained end, positive for the other.
        side = 0
        while hi is not None and hi - lo > max(RESOLUTION, RELATIVE * lo) and len(tried) < TRIALS:
            m_lo = at_lo['margin'] if at_lo else None
            m_hi = at_hi['margin']
            if m_lo is not None and m_lo > 0 > m_hi:
                # False position, halving the margin of an end that has stayed put while the other moved twice.
                w_lo = m_lo / 2 if side >= 2 else m_lo
                w_hi = m_hi / 2 if side <= -2 else m_hi
                r = lo + (hi - lo) * w_lo / (w_lo - w_hi)
            else:
                r = (lo + hi) / 2
            r = min(max(r, lo + 0.1 * (hi - lo)), hi - 0.1 * (hi - lo))
            result = self.trial(rates + [r], choice)
            tried.append(result)
            if result['sustained']:
                lo, at_lo, side = r, result, (side - 1 if side < 0 else -1)
            else:
                hi, at_hi, side = r, result, (side + 1 if side > 0 else 1)
        return lo, hi, at_lo, at_hi, tried

    def march(self, choice, warm=None, label=''):
        """The schedule found period by period on one choice of sales."""
        rates, record = [], []
        saved = self.state.get(label) if self.state is not None else None
        if saved:
            rates, record, self.cap = list(saved['rates']), list(saved['periods']), saved['opening_cap']
            self.log('resuming', label, 'after', len(rates), 'periods')
        if not rates:
            self.cap = None
            opening = self.trial([], choice)
            if not opening['sustained']:
                self.log('the opening settlement is not sustained without growth; no period grows')
                return [0.0] * len(PERIODS), [dict(first_year=a, last_year=b, rate=0.0, binding=opening['binding'], trials=0) for a, b in PERIODS], opening
            # The opening funding no schedule may exceed: that of the constant-growth ramp's rate in the first period,
            # found by the same trial as every other schedule, or that of no growth where that rate is not sustained.
            planned = self.trial([definition.growth_rate(self.s['programme'], REFERENCE, PERIODS[0][0])], choice)
            self.cap = planned.get('opening_funding_USD2041', opening['opening_funding_USD2041']) * (1 + OPENING_TOLERANCE)
        for index in range(len(rates), len(PERIODS)):
            guess = warm[index] if warm else (rates[-1] if rates else 0.033)
            lo, hi, at_lo, at_hi, tried = self.period(rates, index, choice, guess, 1.1 if warm else 2.0)
            a, b = PERIODS[index]
            record.append(dict(first_year=a, last_year=b, rate=lo, lowest_unsustained_rate=hi,
                               binding=(at_hi or at_lo)['binding'], failed=(at_hi or {}).get('failed', []), trials=len(tried)))
            rates.append(lo)
            if index == 0 and at_lo:
                # Later growth leaves the opening funding as the first period set it, but for rounding.
                self.cap = at_lo['opening_funding_USD2041'] * (1 + OPENING_TOLERANCE)
            self.log('period %d-%d: %.4f a year sustained, %s not; binding: %s' % (a, b, lo, 'none tried' if hi is None else '%.4f' % hi, record[-1]['binding']))
            if self.state is not None:
                self.state[label] = dict(rates=rates, periods=record, opening_cap=self.cap)
                self.state.save()
        return rates, record, None

    def solve(self):
        start = time.monotonic()
        choice = self.choose()
        rounds, warm, last = [], None, None
        for number in range(1, ROUNDS + 1):
            rates, record, refused = self.march(choice, warm, label='round %d' % number)
            rounds.append(dict(choice_made_on=REFERENCE if number == 1 else 'the schedule of round %d' % (number - 1), rates=rates,
                               periods=record, lines_in_order=choice['order'], lines_dropped=choice['dropped'], lines_starting_late=choice['starts']))
            if refused:
                last = refused
                break
            warm = rates
            if number < ROUNDS:
                choice = self.choose(schedule(self.s['programme'], rates))
        beyond = schedule(self.s['programme'], rates)
        limited = schedule(self.s['programme'], rates, self.planned)
        reached = next((y for y, v in beyond if v >= self.planned * (1 - 1e-12)), None)
        final = last or self.trial(rates, choice)
        check = self.trial(rates, choice, self.planned, LIMITED) if reached is not None else None
        return dict(scenario=self.scenario, schedules={RAMP: beyond, LIMITED: limited}, rounds=rounds, periods=rounds[-1]['periods'],
                    constant_growth_opening_funding_USD2041=self.reference_opening, opening_funding_cap_USD2041=self.cap,
                    planned_population_reached_year=reached, opening_settlement_sustained=refused is None,
                    final_trial=final, limited_trial=check, trials=self.trials, minutes=round((time.monotonic() - start) / 60, 1))


class State(dict):
    """The rates found so far for one scenario, and the choices of sales made, kept on disk so that an interrupted
    solution resumes. Rates are kept against every file the solution reads; choices against the study's inputs."""

    def __init__(self, scenario, key, inputs_key):
        super().__init__()
        self.path = RUNS / ('%s.json' % scenario)
        self.key, self.inputs_key = key, inputs_key
        if self.path.exists():
            saved = json.loads(self.path.read_text())
            if saved.get('key') == key:
                self.update(saved['rounds'])

    def save(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text(json.dumps(dict(key=self.key, rounds=dict(self)), indent=1))

    def choice_path(self, label):
        return self.path.with_name('%s.%s.pickle' % (self.path.stem, hashlib.sha256(label.encode()).hexdigest()[:12]))

    def kept_choice(self, label):
        """A choice of sales made before on the same inputs, or None."""
        path = self.choice_path(label)
        if not path.exists():
            return None
        saved = pickle.loads(path.read_bytes())
        return saved['choice'] if saved.get('key') == self.inputs_key and saved.get('label') == label else None

    def keep_choice(self, label, choice):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.choice_path(label).write_bytes(pickle.dumps(dict(key=self.inputs_key, label=label, choice=choice)))


def solve(scenario):
    def key(files):
        return hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest()
    return Solver(scenario, state=State(scenario, key(dependencies()), key(inputs.dependencies()))).solve()


def product(scenarios, jobs=1):
    if jobs > 1:
        with multiprocessing.Pool(min(jobs, len(scenarios))) as workers:
            solved = workers.map(solve, scenarios, chunksize=1)
    else:
        solved = [solve(name) for name in scenarios]
    return dict(schema=SCHEMA, producer=dict(path=str(Path(__file__).resolve().relative_to(ROOT)), sha256=digest(__file__)),
                dependencies=dependencies(),
                evidence='Conditional schedules. Each is the fastest growth, period by period, that the coupled study\'s own '
                         'accounts sustain under its stated financing rules, workforce limit and choice of sales; the '
                         'accounts themselves rest on the scenarios of demand and the engineering coefficients, whose '
                         'evidence is stated with them. The water drawn is compared with a conditional survey range, not '
                         'a qualified reserve.',
                units='Years of the programme, residents, yearly growth rates as fractions, kilograms, constant USD2041.',
                reading_rule='For each scenario of demand, the schedule of `%s` grows in each period at the rate given, the '
                             'highest found sustained to within the stated resolution; the lowest rate found unsustained and '
                             'the condition that failed there are given beside it. `%s` is the same schedule held at the '
                             'planned population from the year it reaches it. A trial gives, for one schedule, the conditions '
                             'that fail and the margin of each, as a share of the constant-growth ramp\'s opening funding: '
                             'money, the least cash above the buffer the financing rules require, or the equity needed after '
                             'the opening; opening, what the cap exceeds the trial\'s opening funding by, the cap being the '
                             'opening funding of the constant-growth rate in the first period while that period is sought, '
                             'and the first period\'s chosen opening funding after, each with %g of it added for the rounding '
                             'of the sizing; income, the change '
                             'in what the system holds over years 201 to 220. Hands is the least share of the workforce the '
                             'residents\' own work leaves, in any year; a schedule fails it where that work exceeds the '
                             'workforce, or where the lines and the making of parts could not be fitted to it. A trial also '
                             'gives its receipts and the contained water drawn from the sources to each horizon against the '
                             'wider near-Earth screen at its 5th, 50th and 95th percentiles. While a schedule is sought '
                             'nothing is paid out to the funders; the coupled run finances the solved schedules under the '
                             'study\'s ordinary rules.' % (RAMP, LIMITED, OPENING_TOLERANCE),
                periods=[list(p) for p in PERIODS], highest_rate_tried=HIGHEST, resolution=dict(absolute=RESOLUTION, relative=RELATIVE),
                scenarios={entry['scenario']: entry for entry in solved})


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--scenario', action='append', choices=SCENARIOS, help='a scenario of demand; repeat for several (default: all)')
    parser.add_argument('--output', type=Path, default=PRODUCT)
    parser.add_argument('--jobs', type=int, default=1, help='scenarios solved at once (default: one)')
    args = parser.parse_args()
    result = product(tuple(args.scenario or SCENARIOS), args.jobs)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=1, sort_keys=True, allow_nan=False) + '\n')
    print(args.output)


if __name__ == '__main__':
    main()
