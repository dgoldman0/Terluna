"""U, the money of the zone beyond Earth, beside the account with Earth.

The period model's account with Earth counts what Earth's customers pay and what the settlement pays Earth, in
dollars. This ledger adds U as one unit for the accounts (Terluna's later work splits it into several zones):

- trade moves toward U, a share settled in U rising over the decade after the opening; Earth's customers and sellers
  keep a week of that trade as working balances, without interest, and top them up as U's prices rise;
- paying in U changes Earth's purchases from the zone by the trade multiple, none by default (the evidence on shared
  monies measures partners who use one money; moving from the dollar to U puts an exchange between the zone and
  Earth), with a fall as the case of a thin market or a widely swinging U;
- Earth's investors hold U assets when the real return they expect, the yield and U's last change in value, exceeds
  what they earn on Earth, up to a ceiling in proportion to the zone's output, closing part of the gap each year by
  buying U with dollars or selling it back; their yield is paid in U, and the share of it they spend in the zone buys
  goods the economy charges their work, power and Earth inputs;
- the funders receive what the account leaves above the reserve, part of it in U, which they spend in the zone or
  convert over the years;
- U's value in dollars, the rate, is a result. Until the account with Earth is repaid it stays at par and the funders
  keep the reserve at its target. After that U floats: a year whose dollars leave the reserve below its minimum
  lowers U within the year, so the residents' purchases from Earth absorb the gap, down to the policy's floor, and
  below the floor the reserve itself falls short; a year that leaves the reserve above its target raises U back
  toward par and pays the funders the rest;
- the reserve is a basket: a share in gold and platinum-group metals, taken from the zone's own refining where it can
  be and otherwise bought, earning nothing, and the rest in Earth's currencies at Earth's safe real rate;
- five tests of a safe and healthy policy are checked every year from the opening.

Amounts in dollars are constant dollars of 2041 of purchasing power. Amounts in U are U of constant purchasing power
in the zone, each worth the year's rate in dollars; at par the two are equal. The zone's output is its residents' pay
at the rate and what it keeps of its sales to Earth.
"""
from __future__ import annotations
import math

from .designs import inputs, v

OPENING = 8
PAYMENT_YEARS = 5                   # the years of payments to Earth over which a month of payments is counted
EXPECTATION_YEARS = 5               # the years over which investors read U's change in value
LEVERS = ('settlement_share_from', 'settlement_share_to', 'settlement_years', 'real_yield', 'yield_follows_growth', 'inflation',
          'reserve_share_of_holdings', 'reserve_months_of_payments', 'funders_paid_in_U', 'reserve_metal_share')


def policy(name, overrides=None):
    """The levers of a named policy with any overrides, or None for the dollar account alone."""
    if name in (None, 'dollar'):
        return None
    p = {k: x for k, x in inputs()['u']['policies'][name].items() if k in LEVERS}
    p.update({k: x for k, x in (overrides or {}).items() if x is not None})
    return p


def holders(overrides=None):
    h = {k: v(x) for k, x in inputs()['u']['holders'].items()}
    h.update({k: x for k, x in (overrides or {}).items() if x is not None})
    return h


def reserve_terms():
    return {k: v(x) for k, x in inputs()['u']['reserve'].items()}


def tests():
    return {k: x for k, x in inputs()['u']['tests'].items() if k != 'basis'}


def settlement_share(p, y):
    """The share of the zone's trade with Earth settled in U in a programme year."""
    if p is None or y < OPENING:
        return 0.0
    f = min(1.0, (y - OPENING) / p['settlement_years'])
    return p['settlement_share_from'] + (p['settlement_share_to'] - p['settlement_share_from']) * f


def trade_multiple(p, h, y):
    """Earth's purchases from the zone as a multiple of what they would be in dollars alone."""
    return 1.0 + h['trade_effect_of_full_settlement'] * settlement_share(p, y) if p is not None else 1.0


class Ledger:
    """The U accounts of one run, stepped year by year with the economy. `floor` is the lowest rate U may fall to;
    with `credit_line` the funders lend the reserve what it still lacks of its minimum at the floor, as in the dollar
    account they cover a year's deficit, and their flows carry the loan."""

    def __init__(self, p, h, floor=1.0, credit_line=False):
        self.p, self.h, self.t, self.r = p, h, tests(), reserve_terms()
        self.floor, self.credit_line = floor, credit_line
        self.rate = 1.0               # dollars a U this year; the economy prices the year's local work at it
        self.rates = []               # the closing rate of each year, for the change in value investors expect
        self.repaid = False
        self.balances = 0.0           # Earth's working balances, U
        self.investment = 0.0         # Earth's investors' U assets, U
        self.funders_U = 0.0          # the funders' U, held apart from the investors', U
        self.funders_U_due = 0.0      # U paid to the funders last year, reaching them this year
        self.currency = 0.0           # the reserve's Earth currencies, dollars, the zone's own dollars among them
        self.metals = 0.0             # the reserve's gold and platinum-group metals, dollars
        self.own = 0.0                # dollars the zone kept when the funders took U: they back the funders' U, and
                                      # what the funders spend in the zone leaves them with the zone as its savings
        self.spent_last = 0.0         # dollars of U holders spent in the zone last year, charged in this year's economy
        self.offered = 0.0            # dollars of metals the reserve asked of this year's refining
        self.payments = []            # payments to Earth in each year since the opening, for the months of payments
        self.output = []
        self.funders_cumulative = 0.0

    def investment_target(self, output, expected):
        h = self.h
        excess = expected - h['earth_safe_real_rate'] - h['investment_threshold_spread']
        if excess <= 0 or output <= 0:
            return 0.0
        return h['investment_ceiling_share_of_output'] * (1 - math.exp(-excess / h['investment_spread_scale'])) * output

    def growth(self, n):
        """The zone's output growth a year over the last n years, or None before there are n years since the opening."""
        if len(self.output) <= n or self.output[-1 - n] <= 0:
            return None
        return (self.output[-1] / self.output[-1 - n]) ** (1 / n) - 1

    def metals_wanted(self):
        """The dollars of metals the reserve lacks of its share, which it asks first of this year's refining."""
        if self.p is None:
            return 0.0
        self.offered = max(0.0, self.p.get('reserve_metal_share', 0.0) * (self.currency + self.metals) - self.metals)
        return self.offered

    def step(self, y, receipts, paid_other, earth_pay_par, pay_par, kept_sales, metals_retained=0.0):
        """One year. `receipts` and `paid_other` are the account with Earth in dollars at this year's opening rate,
        without the residents' purchases from Earth, which are `earth_pay_par` at par and follow the rate; `pay_par` is
        the residents' pay and `kept_sales` what the zone keeps of its sales to Earth; `metals_retained` the dollars of
        metals the zone's refining kept for the reserve instead of selling."""
        p, h, t, rt = self.p, self.h, self.t, self.r
        e0 = self.rate
        share = settlement_share(p, y)
        inflation = p['inflation'] if y >= OPENING else 0.0
        output0 = pay_par * e0 + kept_sales
        self.output.append(output0)
        n = t['yield_over_growth_years']
        growth = self.growth(n) if y - OPENING >= n else None
        # The yield offered; a policy whose yield follows growth holds it within the zone's growth of the last 20 years.
        real_yield = p['real_yield']
        if p.get('yield_follows_growth') and growth is not None:
            real_yield = max(0.0, min(real_yield, growth))
        # Working balances: a week of the trade settled in U, held in U and topped up as U's prices rise.
        trade = receipts + paid_other + earth_pay_par * e0
        wanted = share * trade * h['working_balance_days'] / 365.25 / e0
        from_balances = (wanted - self.balances + inflation * self.balances) * e0
        self.balances = wanted
        # Investment holdings, bought with dollars: the yield is paid in U, and the share spent in the zone returns as
        # purchases of its goods, retiring that U; the investors expect the yield and U's average yearly change in value
        # over the last five years, and close part of the gap to what they want each year.
        chi = h['share_of_income_spent_in_zone']
        income = real_yield * self.investment
        spent = chi * income
        held = self.investment + income - spent
        past = self.rates[-EXPECTATION_YEARS:]
        change = (e0 / past[0]) ** (1 / len(past)) - 1 if past else 0.0
        target = self.investment_target(output0, real_yield + change) / e0
        moved = h['adjustment_per_year'] * (target - held)
        self.investment = held + moved
        from_investors = moved * e0
        # The funders' U: they spend a share of what reaches them in the zone and convert the rest over the years, each
        # conversion paid from the dollars the zone kept for them.
        spent_by_funders = chi * self.funders_U_due
        held_by_funders = self.funders_U + self.funders_U_due - spent_by_funders
        converted_U = h['adjustment_per_year'] * held_by_funders
        self.funders_U = held_by_funders - converted_U
        converted = converted_U * e0
        holdings_U = self.balances + self.investment + self.funders_U
        self.spent_last = (spent + spent_by_funders) * e0
        # Dollars: the account with Earth at the opening rate, Earth's purchases of U net of sales, the funders'
        # conversions and the currency's interest; metals kept from the zone's refining join the reserve directly.
        interest = h['earth_safe_real_rate'] * self.currency
        own_interest = h['earth_safe_real_rate'] * max(0.0, self.own)
        flow = receipts - paid_other - earth_pay_par * e0 + from_balances + from_investors + interest - converted
        self.own += own_interest - converted
        self.metals += metals_retained
        pool = self.currency + self.metals + flow
        # Months of payments are counted on the payments of the last five years since the opening, so that a year's
        # lump of equipment does not by itself double the reserve the zone must hold, and with the residents' purchases
        # from Earth at par, so that the reserve the zone must hold does not fall because U has fallen.
        at_par = paid_other + earth_pay_par
        recent = self.payments[-(PAYMENT_YEARS - 1):]
        usual = (sum(recent) + at_par) / (len(recent) + 1)
        if y >= OPENING:
            self.payments.append(at_par)
        if y < OPENING:
            target_level, minimum = 0.0, 0.0
        else:
            target_level = max(p['reserve_share_of_holdings'] * holdings_U * e0, p['reserve_months_of_payments'] / 12 * usual)
            minimum = max(t['reserve_share_of_holdings_min'] * holdings_U * e0, t['reserve_months_of_payments_min'] / 12 * usual)
        rate = e0
        funders, in_U = 0.0, 0.0
        level = target_level + max(0.0, self.own)
        if not self.repaid:
            # The build: U at par, the funders keeping the reserve at its target and taking what is above it.
            rate = 1.0
            funders = pool - level
        elif pool >= level:
            # Above the target: U back toward par first, then the funders.
            excess = pool - level
            lift = min(excess, (1.0 - e0) * earth_pay_par)
            rate = e0 + (lift / earth_pay_par if earth_pay_par > 0 else 0.0)
            funders = excess - lift
            pool -= lift
        elif pool < minimum:
            # Below the minimum: U falls within the year, the residents buying less from Earth, down to the floor.
            gap = minimum - pool
            drop = min(gap, max(0.0, e0 - self.floor) * earth_pay_par)
            rate = e0 - (drop / earth_pay_par if earth_pay_par > 0 else 0.0)
            pool += drop
            if self.credit_line and pool < minimum:
                funders = pool - minimum
        # Between the minimum and the target the reserve absorbs the year and rebuilds before anything is paid out.
        if funders > 0:
            in_U = p['funders_paid_in_U'] * funders
            self.own += in_U
            pool -= funders - in_U
        elif funders < 0:
            pool -= funders
        if pool < 0:                    # nothing left even at the floor: the funders cover it
            funders += pool
            pool = 0.0
        if self.own < 0:                # conversions beyond the zone's own dollars: the reserve's other dollars pay them
            self.own = 0.0
        # The basket: what the share lacked and this year's refining did not supply is bought against the currency; the
        # share this year's growth adds is asked of next year's refining; metals beyond the share are sold.
        target_metals = p.get('reserve_metal_share', 0.0) * pool
        bought = min(max(0.0, self.offered - metals_retained), max(0.0, target_metals - self.metals))
        self.offered = 0.0
        self.metals = min(self.metals + bought, target_metals)
        self.currency = pool - self.metals
        if self.currency < 0:
            self.metals += self.currency
            self.currency = 0.0
        reserve = self.currency + self.metals
        self.funders_U_due = in_U / rate
        self.funders_cumulative += funders
        if y >= OPENING and self.funders_cumulative >= 0 and y > OPENING:
            self.repaid = True
        self.rates.append(rate)
        self.rate = rate
        earth_pay = earth_pay_par * rate
        paid = paid_other + earth_pay
        output = pay_par * rate + kept_sales
        holdings = holdings_U * rate
        checks = {}
        if y >= OPENING:
            # The reserve and the run are read on the terms the year's policy was set on: Earth's holdings at the opening
            # rate and the usual payments at par.
            liquid = self.currency + (1 - rt['metal_haircut']) * self.metals
            checks = dict(
                inflation=inflation <= min(t['inflation_year_max'], t['inflation_decade_mean_max']),
                reserve=reserve >= minimum * (1 - 1e-9),
                run=liquid >= t['run_share_of_holdings'] * holdings_U * e0 * (1 - 1e-9),
                yield_within_growth=self.investment <= 0 or growth is None or real_yield <= max(0.0, growth) + 1e-12,
                holdings_within_output=holdings <= t['holdings_over_output_max'] * output)
        return dict(u_settlement_share=share, u_inflation=inflation, u_real_yield=real_yield, u_rate=rate, u_rate_opening=e0,
                    u_at_floor=rate <= self.floor + 1e-12 and self.repaid and self.floor < 1.0, u_repaid=self.repaid,
                    u_balances=self.balances, u_investment=self.investment, u_funders_U=self.funders_U, u_holdings=holdings_U,
                    u_holdings_USD=holdings, u_investment_target=target, u_reserve=reserve, u_reserve_currency=self.currency,
                    u_reserve_metals=self.metals, u_metals_retained=metals_retained, u_metals_bought=bought,
                    u_reserve_target=target_level, u_reserve_minimum=minimum, u_reserve_own=self.own, u_interest=interest,
                    u_dollars_to_funders=funders - in_U, u_funders_converted=converted, u_from_balances=from_balances,
                    u_from_investors=from_investors, u_yield_income=income * e0, u_spent_in_zone=self.spent_last,
                    u_paid_to_funders=in_U, u_output=output, u_growth_20=growth, u_earth_pay=earth_pay, u_paid=paid,
                    funders_flow=funders, funders_cumulative=self.funders_cumulative, u_tests=checks)
