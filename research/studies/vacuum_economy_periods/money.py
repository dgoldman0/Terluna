"""U, the money of the zone beyond Earth, beside the account with Earth.

The period model's account with Earth counts what Earth's customers pay and what the settlement pays Earth, in
dollars. This ledger adds U as one unit for the accounts (Terluna's later work splits it into several zones):

- trade moves toward U, a share settled in U rising over the decade after the opening; Earth's customers and sellers
  keep a week of that trade as working balances, without interest, and top them up as U's prices rise;
- the shared unit raises Earth's purchases from the zone in proportion to the share settled (the trade multiple,
  applied to the lines of sale in the economy's own year so that their power and hands follow);
- Earth's investors hold U assets, bonds or staking, when the real yield offered exceeds what they earn on Earth, up to
  a ceiling in proportion to the zone's output, closing part of the gap each year by buying U with dollars or converting
  it back; their yield is paid in U, and the share of it they spend in the zone costs the zone its own work, not dollars;
- the funders receive what the account leaves above the reserve, part of it in U, which joins Earth's holdings;
- a dollar reserve is kept at the policy's share of Earth's holdings and months of payments to Earth, the funders
  topping it up when the year's dollars fall short;
- five tests of a safe and healthy policy are checked every year from the opening.

Amounts are in constant dollars of 2041 of purchasing power. The zone's output is its residents' pay, the workforce at
the wage, the only domestic product the model counts.
"""
from __future__ import annotations
import math

from .designs import inputs, v

OPENING = 8
LEVERS = ('settlement_share_from', 'settlement_share_to', 'settlement_years', 'real_yield', 'yield_follows_growth', 'inflation',
          'reserve_share_of_holdings', 'reserve_months_of_payments', 'funders_paid_in_U')


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
    """The U accounts of one run, stepped year by year with the economy."""

    def __init__(self, p, h):
        self.p, self.h, self.t = p, h, tests()
        self.balances = 0.0           # Earth's working balances of U
        self.investment = 0.0         # Earth's investors' U assets, the funders' U among them
        self.reserve = 0.0            # the zone's dollar reserve: the floor and the zone's own dollars above it
        self.own = 0.0                # dollars the zone kept when the funders took U: they back the funders' U, and
                                      # what the funders spend in the zone leaves them with the zone as its savings
        self.funders_U = 0.0          # the funders' U, held apart from the investors' holdings
        self.funders_U_due = 0.0      # U paid to the funders last year, reaching them this year
        self.output = []
        self.funders_cumulative = 0.0

    def investment_target(self, output, real_yield):
        h = self.h
        excess = real_yield - h['earth_safe_real_rate'] - h['investment_threshold_spread']
        if excess <= 0 or output <= 0:
            return 0.0
        return h['investment_ceiling_share_of_output'] * (1 - math.exp(-excess / h['investment_spread_scale'])) * output

    def growth(self, n):
        """The zone's output growth a year over the last n years, or None before there are n years since the opening."""
        if len(self.output) <= n or self.output[-1 - n] <= 0:
            return None
        return (self.output[-1] / self.output[-1 - n]) ** (1 / n) - 1

    def step(self, y, receipts, paid, output):
        """One year: receipts and payments in dollars from the account with Earth, and the zone's output (the
        residents' pay and what the zone keeps of its sales to Earth)."""
        p, h, t = self.p, self.h, self.t
        share = settlement_share(p, y)
        inflation = p['inflation'] if y >= OPENING else 0.0
        self.output.append(output)
        n = t['yield_over_growth_years']
        growth = self.growth(n) if y - OPENING >= n else None
        # The yield offered; a policy whose yield follows growth holds it within the zone's growth of the last 20 years.
        real_yield = p['real_yield']
        if p.get('yield_follows_growth') and growth is not None:
            real_yield = max(0.0, min(real_yield, growth))
        # Working balances: a week of the trade settled in U, topped up as U's prices rise.
        wanted = share * (receipts + paid) * h['working_balance_days'] / 365.25
        from_balances = wanted - self.balances + inflation * self.balances
        self.balances = wanted
        # Investment holdings, bought with dollars: the yield is paid in U, and the share spent in the zone returns as
        # purchases of its goods, retiring that U; the investors close part of the gap to what they want each year.
        chi = h['share_of_income_spent_in_zone']
        income = real_yield * self.investment
        spent = chi * income
        held = self.investment + income - spent
        target = self.investment_target(output, real_yield)
        from_investors = h['adjustment_per_year'] * (target - held)
        self.investment = held + from_investors
        # The funders' U: they spend a share of what reaches them in the zone and convert the rest over the years,
        # each conversion paid from the dollars the zone kept for them.
        spent_by_funders = chi * self.funders_U_due
        held_by_funders = self.funders_U + self.funders_U_due - spent_by_funders
        converted = h['adjustment_per_year'] * held_by_funders
        self.funders_U = held_by_funders - converted
        holdings = self.balances + self.investment + self.funders_U
        # Dollars: the account with Earth, Earth's purchases of U net of conversions, and the reserve's interest; the
        # interest on the zone's own dollars stays with them.
        interest = h['earth_safe_real_rate'] * self.reserve
        own_interest = h['earth_safe_real_rate'] * max(0.0, self.own)
        before = self.reserve + receipts - paid + from_balances + from_investors + interest - converted
        self.own += own_interest - converted
        floor = 0.0 if y < OPENING else max(p['reserve_share_of_holdings'] * holdings, p['reserve_months_of_payments'] / 12 * paid)
        # What is left above the floor and the zone's own dollars goes to the funders, part of it in U, whose dollars
        # the zone keeps; a shortfall comes first from the zone's own dollars and then from the funders.
        if before >= floor + self.own:
            excess = before - floor - self.own
            funders, in_U = excess, p['funders_paid_in_U'] * excess
            self.own += in_U
        else:
            short = floor + self.own - before
            used = min(short, max(0.0, self.own))
            self.own -= used
            funders, in_U = -(short - used), 0.0
        if self.own < 0:                # conversions beyond the zone's own dollars: the funders cover them
            funders += self.own
            self.own = 0.0
        self.reserve = floor + self.own
        self.funders_U_due = in_U
        self.funders_cumulative += funders
        checks = {}
        if y >= OPENING:
            checks = dict(
                inflation=inflation <= min(t['inflation_year_max'], t['inflation_decade_mean_max']),
                reserve=self.reserve >= max(t['reserve_share_of_holdings_min'] * holdings, t['reserve_months_of_payments_min'] / 12 * paid) * (1 - 1e-12),
                run=self.reserve >= t['run_share_of_holdings'] * holdings * (1 - 1e-12),
                yield_within_growth=self.investment <= 0 or growth is None or real_yield <= max(0.0, growth) + 1e-12,
                holdings_within_output=holdings <= t['holdings_over_output_max'] * output)
        return dict(u_settlement_share=share, u_inflation=inflation, u_real_yield=real_yield, u_balances=self.balances, u_investment=self.investment,
                    u_holdings=holdings, u_investment_target=target, u_reserve=self.reserve, u_reserve_floor=floor, u_reserve_own=self.own,
                    u_interest=interest, u_dollars_to_funders=funders - in_U, u_funders_U=self.funders_U, u_funders_converted=converted,
                    u_from_balances=from_balances, u_from_investors=from_investors, u_yield_income=income,
                    u_spent_in_zone=spent + spent_by_funders, u_paid_to_funders=in_U, u_output=output, u_growth_20=growth,
                    funders_flow=funders, funders_cumulative=self.funders_cumulative, u_tests=checks)
