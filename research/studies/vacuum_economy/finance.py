"""The financing of the coupled study: the dated plan of a case's money, its financing under each arrangement, and
the measures reported at each horizon.

These are the staged studies' functions as they stood at e49fca09, taken whole from research/studies/enclosed_returns (finance.plan,
finance.finance, finance.completion_cash_requirement, finance.metrics, revised.measures, revised.sponsor_cash,
bottlenecks.carry_bridge) and research/studies/complex_returns (credit.debt_due, credit.repay, credit.price_screen,
consolidated.pv, consolidated.durable_recovery, consolidated.nonnegative_irr), each written one statement to a line
with its syntax unchanged, so that the study no longer reads those studies. `cfg` is the money section of the
programme definition.
"""
from __future__ import annotations
from copy import deepcopy
import math


# From research/studies/complex_returns/consolidated.py at e49fca09, pv.
def pv(values, rate):
    return math.fsum((v * (1 + rate) ** (-t) for t, v in enumerate(values)))

# From research/studies/complex_returns/consolidated.py at e49fca09, durable_recovery.
def durable_recovery(flows, rate=0.0):
    cumulative = []
    balance = 0.0
    for year, value in enumerate(flows):
        balance += value / (1 + rate) ** year
        cumulative.append(balance)
    running = math.inf
    result = None
    for year in reversed(range(len(cumulative))):
        running = min(running, cumulative[year])
        if running >= 0:
            result = year
    return result

# From research/studies/complex_returns/consolidated.py at e49fca09, nonnegative_irr.
def nonnegative_irr(flows):
    signs = [1 if x > 0 else -1 for x in flows if x]
    changes = sum((a != b for a, b in zip(signs, signs[1:])))
    if not signs or signs[0] != -1 or changes != 1 or (pv(flows, 0.0) < 0):
        return None
    lo, hi = (0.0, 1.0)
    while pv(flows, hi) > 0:
        hi *= 2
        if hi > 1024:
            return None
    # All actual roots fall below100%; avoid high powers at long horizons.
    for _ in range(70):
        mid = (lo + hi) / 2
        if pv(flows, mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2

# From research/studies/complex_returns/credit.py at e49fca09, debt_due.
def debt_due(loans, year):
    interest = principal = U_principal = external_interest = U_interest = supplier_principal = supplier_interest = 0.0
    for loan in loans:
        if loan['remaining'] <= 0 or year <= loan['year']:
            continue
        charge = loan['remaining'] * loan['rate']
        repay = min(loan['remaining'], loan['original'] / loan['term']) if year > loan['year'] + loan['grace'] else 0.0
        interest += charge
        principal += repay
        if loan['kind'] == 'U':
            U_principal += repay
            U_interest += charge
            supplier_principal += repay * loan.get('supplier_share', 0.0)
            supplier_interest += charge * loan.get('supplier_share', 0.0)
        else:
            external_interest += charge
    return dict(interest=interest, principal=principal, U_principal=U_principal, external_interest=external_interest, U_interest=U_interest, supplier_principal=supplier_principal, supplier_interest=supplier_interest)

# From research/studies/complex_returns/credit.py at e49fca09, repay.
def repay(loans, year):
    due = debt_due(loans, year)
    for loan in loans:
        if year > loan['year'] + loan['grace']:
            loan['remaining'] -= min(loan['remaining'], loan['original'] / loan['term'])
    return due

# From research/studies/complex_returns/credit.py at e49fca09, price_screen.
def price_screen(basket, local_spending, local_payout, gross, U_outstanding, U_escrow, contractor_escrow, initial_money_need, age, monetary, outside_U=0.0):
    planned = basket + local_spending + max(0.0, local_payout)
    local_need = monetary['domestic_holding_days'] / 365 * planned
    settlement = monetary['settlement_initial_share'] + (monetary['settlement_final_share'] - monetary['settlement_initial_share']) * min(1.0, age / monetary['settlement_transition_years'])
    customer_need = settlement * gross * monetary['customer_holding_days'] / 365
    money_growth = max(0.0, local_need + customer_need - initial_money_need)
    circulating_credit = max(0.0, U_outstanding - U_escrow - contractor_escrow - outside_U)
    excess = max(0.0, circulating_credit - money_growth)
    spending = planned + excess * 365 / monetary['domestic_holding_days']
    price = (spending / basket) ** (1 / (1 + monetary['supply_elasticity']))
    return dict(price=price, planned_spending=planned, transaction_money_growth=money_growth, circulating_credit=circulating_credit, excess_credit_balances=excess, local_transaction_balances=local_need, customer_balances=customer_need, other_money_creation=max(0.0, money_growth - circulating_credit), total_money_creation=max(money_growth, circulating_credit))

# From research/studies/enclosed_returns/bottlenecks.py at e49fca09, carry_bridge.
def carry_bridge(path, sim, cfg, horizon, rate):
    """Reconcile PV of real cash to payouts plus terminal net financial claims."""
    start = path['annual'][0]['year']
    assets = 0.0
    debt = 0.0
    carry = interest_vs_discount = 0.0
    fs = {f['year']: f for f in sim['annual']}
    for y in range(horizon + 1):
        if y:
            carry += rate * assets / (1 + rate) ** y
            if y in fs:
                interest_vs_discount += (fs[y]['interest'] - rate * debt) / (1 + rate) ** y
        if y == 0:
            assets += sim['initial_reserve']
        if y == start - 1:
            assets += path['initial_working']
        if y in fs:
            f = fs[y]
            assets = f['cash'] + f['reserves'] + f['owner_custody']
            debt = f['debt']
    metric = metrics(path, sim, cfg, horizon)['discounted'][str(rate)]
    wedge = metric['unlevered_cash_NPV'] - metric['sponsor_cash_and_terminal_NPV']
    residual = wedge - carry - interest_vs_discount
    if not math.isclose(residual, 0.0, abs_tol=0.25):
        raise ValueError(f'Financial carry bridge failed by {residual}')
    return dict(unlevered_cash_NPV=metric['unlevered_cash_NPV'], sponsor_distributed_cash_NPV=metric['sponsor_distributed_cash_NPV'], terminal_net_financial_claim_PV=metric['terminal_financial_claim_PV'], zero_real_yield_asset_carry_cost_PV=carry, interest_minus_discount_cost_on_debt_PV=interest_vs_discount, reconciliation_residual_USD2041=residual)

# From research/studies/enclosed_returns/finance.py at e49fca09, plan.
def plan(case, physical, cfg):
    if not physical['convergence']['converged'] or physical['first_capacity_shortfall'] is not None:
        raise ValueError('Physical schedule has unresolved requirements: ' + case['name'])
    requests = {r['year']: r for r in case['annual']}
    start = min(requests)
    by_year = {r['year']: r for r in physical['annual']}
    residual_stock = {y: max(0.0, r['service_installed_capital_USD2041'] - by_year[y]['cash']['identified_service_stock_value_USD2041']) for y, r in requests.items()}
    quote = {y: 0.0 for y in requests}
    preproduction_service = 0.0
    for y, r in requests.items():
        if y == start:
            continue
        addition = max(0.0, residual_stock[y] - residual_stock[y - 1])
        for workyear in (y - 2, y - 1):
            if workyear in quote:
                quote[workyear] += addition / 2
            else:
                preproduction_service += addition / 2
    first = requests[start]
    sc = physical['seed']['cost']
    seed_service = residual_stock[start] + preproduction_service
    initial_hotel_payroll = first['visitor_beds'] * cfg['hotel_staff_per_bed'] * cfg['hotel_staff_USD2041_year']
    seed_carry = (cfg['preopening_commercial_fixed_USD2041_year'] + initial_hotel_payroll) * start / 2
    seed_direct = sc['total_USD2041'] + seed_service + seed_carry
    seed_allowance = cfg.get('capital_cost_allowance_fraction', 0.0) * max(0.0, seed_direct - physical['seed']['total_resource_provision_USD2041'])
    seed_hardware = seed_direct + seed_allowance
    rows = []
    for pr in physical['annual']:
        r = requests[pr['year']]
        cash = pr['cash']
        year = pr['year']
        hotel_payroll = r['visitor_beds'] * cfg['hotel_staff_per_bed'] * cfg['hotel_staff_USD2041_year']
        # Value known assets against the package before dating its residual.
        # Actual replacement peaks therefore cannot lose overlap credits merely
        # because the unitemized quote uses a smooth annual allowance.
        service_cap = quote[year]
        service_renew = max(0.0, r['service_replacement_allowance_USD2041'] - cash['identified_service_mean_renewal_allowance_USD2041'])
        external_cap = cash['capital_external_USD2041'] + (1 - cfg['service_capital_local_fraction']) * service_cap
        local_cap = cash['capital_local_USD2041'] + cfg['service_capital_local_fraction'] * service_cap
        allowance = cfg.get('capital_cost_allowance_fraction', 0.0) * (external_cap + local_cap)
        external_cap += allowance
        external_op = sum((cash[k] for k in cfg['operating_external_keys'])) + r['commercial_operating_USD2041'] - r['commercial_payroll_USD2041'] + service_renew
        local_op = cash['operating_local_payroll_USD2041'] + r['commercial_payroll_USD2041'] + hotel_payroll
        total = cash['operating_total_USD2041'] + r['commercial_operating_USD2041'] + service_renew + hotel_payroll
        capital_provision = cash['capital_resource_provision_USD2041']
        if min(external_op, local_op, external_cap, local_cap) < -1:
            raise ValueError('Negative cash category')
        rows.append(dict(year=year, population=r['population'], service_kits=r['service_kit_count'], capital=external_cap + local_cap + capital_provision, external_capital=external_cap, local_capital=local_cap, capital_cost_allowance=allowance, service_residual_capital=service_cap, service_residual_renewal=service_renew, gross=r['gross_USD2041'], operating_cost=total, operating_available=r['gross_USD2041'] - total, external_operating=external_op, resource_provision=cash['resource_provision_USD2041'], capital_resource_provision=capital_provision, resident_payroll=cash['resident_payroll_USD2041'] + r['commercial_payroll_USD2041'] + hotel_payroll, resident_value=r['population'] * cfg['resident_service_value_USD2041_person_year'], domestic_transaction_basket=r['population'] * cfg['domestic_transaction_basket_USD2041_person_year'], revenue_by_sector=r['revenue_by_sector_USD2041'], maximum_capacity_use=pr['maximum_capacity_use']))
    working = rows[0]['operating_cost'] * cfg['initial_runway_years']
    for i, r in enumerate(rows):
        # The lender examines the remaining programme, including scheduled later
        # investment. This is not an attribution of marginal profit to residents.
        future = rows[i + 1:i + 1 + cfg['underwriting']['project_appraisal_years']]
        r['programme_NPV_screen'] = -r['capital'] + math.fsum(((v['operating_available'] - v['capital']) / (1 + cfg['underwriting']['project_discount_rate']) ** (j + 1) for j, v in enumerate(future)))
    return dict(name=case['name'], portfolio=case['portfolio'], market=case['market'], growth=case['market_growth_name'], seed_direct_cost=seed_direct, seed_capital_cost_allowance=seed_allowance, seed_hardware_and_carry=seed_hardware, seed_resource_provision=physical['seed']['total_resource_provision_USD2041'], seed_resident_payroll=physical['seed']['resident_payroll_USD2041'] + (cfg['preopening_commercial_payroll_USD2041_year'] + initial_hotel_payroll) * start / 2, seed_cost_scope=physical['seed'].get('evidence', ''), unpriced_deployment_requirements=physical['seed'].get('remaining_unpriced_requirements', []), seed_population=first['population'], initial_working=working, seed=seed_hardware + working, annual=rows)

# From research/studies/enclosed_returns/finance.py at e49fca09, completion_cash_requirement.
def completion_cash_requirement(path, index, loans, reserve, cfg, receipt_shock=None):
    """Cash needed now to meet the remaining stated schedule without new loans.

    This conditional forecast permits early distributions only after the worst
    future cumulative cash deficit and its required operating buffer are funded.
    Future sales are scenario inputs, not contracted certainty.
    """
    mon = cfg['monetary']
    shadow = deepcopy(loans)
    start = path['annual'][0]['year']
    cumulative = 0.0
    required = 0.0
    future_reserve = reserve
    for r in path['annual'][index + 1:]:
        y = r['year']
        due = repay(shadow, y)
        service = due['interest'] + due['principal']
        factor = receipt_shock['factor'] if receipt_shock and receipt_shock['start'] <= y <= receipt_shock['end'] else 1.0
        available = r['gross'] * factor - r['operating_cost']
        Udebt = sum((l['remaining'] for l in shadow if l['kind'] == 'U'))
        target = mon['base_reserve_years_core_imports'] * r['external_operating'] + (mon['U_credit_conversion_buffer_fraction'] + mon['U_credit_loss_capital_fraction']) * Udebt
        addition = min(max(0.0, target - future_reserve), mon['reserve_contribution_fraction'] * max(0.0, available - service))
        future_reserve += addition
        cumulative += available - service - addition - r['capital']
        buffer = r['operating_cost'] * max(cfg['liquidity_buffer_years'], cfg['initial_runway_years'] - (y - start))
        required = max(required, buffer - cumulative)
    return required

# From research/studies/enclosed_returns/finance.py at e49fca09, finance.
def finance(path, cfg, mode, control=None, receipt_shock=None):
    if cfg.get('cash_real_return', 0.0) != 0.0:
        raise ValueError('Nonzero cash investment return is not implemented')
    mon = cfg['monetary']
    uw = cfg['underwriting']
    start = path['annual'][0]['year']
    reserve = path['seed'] * mon['initial_reserve_fraction_of_seed']
    initial_reserve = reserve
    cash = path['initial_working']
    loans = []
    out = []
    initial_need = None
    previous_relative = 1.0
    custody = 0.0
    controls = {} if control is None else {r['year']: r for r in control['annual']}
    for i, r in enumerate(path['annual']):
        y = r['year']
        opening = cash
        due = repay(loans, y)
        service = due['interest'] + due['principal']
        factor = receipt_shock['factor'] if receipt_shock and receipt_shock['start'] <= y <= receipt_shock['end'] else 1.0
        gross = r['gross'] * factor
        available = gross - r['operating_cost']
        forecast_shock = receipt_shock if receipt_shock and y >= receipt_shock.get('known_from', start) else None
        debt = sum((l['remaining'] for l in loans))
        Udebt = sum((l['remaining'] for l in loans if l['kind'] == 'U'))
        outside = sum((l['remaining'] * l.get('supplier_share', 0.0) for l in loans if l['kind'] == 'U'))
        target = mon['base_reserve_years_core_imports'] * r['external_operating'] + (mon['U_credit_conversion_buffer_fraction'] + mon['U_credit_loss_capital_fraction']) * Udebt
        contribution = min(max(0.0, target - reserve), mon['reserve_contribution_fraction'] * max(0.0, available - service))
        reserve += contribution
        cash += available - service - contribution
        runway = max(cfg['liquidity_buffer_years'], cfg['initial_runway_years'] - (y - start))
        buffer = r['operating_cost'] * runway
        gap = max(0.0, r['capital'] - max(0.0, cash - buffer))
        annual_committed = sum((l['remaining'] * l['rate'] + l['original'] / l['term'] for l in loans if l['remaining'] > 0))
        headroom = max(0.0, max(0.0, available) * uw['pledgeable_operating_cash_fraction'] / uw['min_DSCR'] - annual_committed)
        dscr = available * uw['pledgeable_operating_cash_fraction'] / service if service else None
        lock = dscr is not None and dscr < uw['min_DSCR']
        draw = min(gap, uw['max_loan_fraction_of_budget'] * r['capital'], max(0.0, uw['max_debt_fraction_of_seed'] * path['seed'] - debt), headroom / (cfg['external_rate'] + 1 / cfg['loan_term_years'])) if mode != 'equity' and r['programme_NPV_screen'] > 0 and (not lock) else 0.0
        reserve_free = max(0.0, reserve - mon['base_reserve_years_core_imports'] * r['external_operating'])
        ucap = max(0.0, reserve_free / (mon['U_credit_conversion_buffer_fraction'] + mon['U_credit_loss_capital_fraction']) - Udebt)
        supplier = min(draw, ucap, uw['max_loan_fraction_of_budget'] * cfg['supplier_U_fraction'] * r['external_capital']) if mode == 'U_credit' else 0.0
        local = min(draw - supplier, max(0.0, ucap - supplier), uw['max_loan_fraction_of_budget'] * r['local_capital']) if mode == 'U_credit' else 0.0
        external = draw - supplier - local
        basket = max(1.0, r['domestic_transaction_basket'] * cfg.get('monetary_basket_multiplier', 1.0))
        local_flow = r['local_capital']
        if initial_need is None:
            initial_need = mon['domestic_holding_days'] / 365 * (basket + local_flow) + mon['settlement_initial_share'] * gross * mon['customer_holding_days'] / 365

        def pressure(local_draw, payout=0.0):
            return price_screen(basket, local_flow, payout, gross, Udebt + supplier + local_draw, 0.0, 0.0, initial_need, y - start, mon, outside_U=outside + supplier)
        pr = pressure(local)
        reference = controls.get(y, {}).get('price', pr['price'])
        ceiling = min(1 + mon['cumulative_relative_price_limit'], max(1.0, previous_relative) * (1 + mon['annual_relative_price_limit']))
        if pr['price'] / reference > ceiling:
            lo, hi = (0.0, local)
            for _ in range(45):
                mid = (lo + hi) / 2
                if pressure(mid)['price'] / reference <= ceiling:
                    lo = mid
                else:
                    hi = mid
            local = lo
            pr = pressure(local)
        draw = external + supplier + local
        for value, kind, rate, share in ((external, 'external', cfg['external_rate'], 0.0), (supplier, 'U', cfg['external_rate'], 1.0), (local, 'U', cfg['local_U_rate'], 0.0)):
            if value > 1e-06:
                loans.append(dict(year=y, original=value, remaining=value, rate=rate, term=cfg['loan_term_years'], grace=uw['principal_grace_years'], kind=kind, supplier_share=share))
        before_topup = cash + draw - r['capital']
        topup = max(0.0, buffer - before_topup)
        free_opening = max(0.0, opening - buffer)
        opgap = max(0.0, -free_opening - available)
        debtgap = max(0.0, service - max(0.0, free_opening + available))
        capgap = max(0.0, r['capital'] - max(0.0, free_opening + available - service - contribution + draw))
        purposes = dict(operating=min(topup, opgap), debt=min(max(0.0, topup - opgap), debtgap), capital=min(max(0.0, topup - opgap - debtgap), capgap))
        purposes['liquidity'] = topup - sum(purposes.values())
        cash = before_topup + topup
        nxt = path['annual'][min(i + 1, len(path['annual']) - 1)]
        nd = debt_due(loans, y + 1)
        ns = nd['interest'] + nd['principal']
        next_U = sum((l['remaining'] for l in loans if l['kind'] == 'U')) - nd['U_principal']
        next_factor = forecast_shock['factor'] if forecast_shock and forecast_shock['start'] <= y + 1 <= forecast_shock['end'] else 1.0
        next_available = nxt['gross'] * next_factor - nxt['operating_cost']
        nt = mon['base_reserve_years_core_imports'] * nxt['external_operating'] + (mon['U_credit_conversion_buffer_fraction'] + mon['U_credit_loss_capital_fraction']) * next_U
        nc = min(max(0.0, nt - reserve), mon['reserve_contribution_fraction'] * max(0.0, next_available - ns))
        nb = nxt['operating_cost'] * max(cfg['liquidity_buffer_years'], cfg['initial_runway_years'] - (y + 1 - start))
        retention = max(buffer, nb + nxt['capital'] + ns + nc - next_available)
        if cfg.get('distribution_policy') == 'completion_funded':
            retention = max(retention, completion_cash_requirement(path, i, loans, reserve, cfg, forecast_shock))
        hold = cfg.get('distribution_policy') == 'retain_to_completion' and y <= cfg['retain_through_year']
        dividend = 0.0 if lock or hold else max(0.0, cash - retention)
        cash -= dividend
        intended = dividend * mon['incremental_distribution_spending_fraction']
        releasable = custody + intended
        released = releasable
        pr = pressure(local, released)
        if control is None:
            reference = pr['price']
        elif pr['price'] / reference > ceiling:
            lo, hi = (0.0, released)
            for _ in range(45):
                mid = (lo + hi) / 2
                if pressure(local, mid)['price'] / reference <= ceiling:
                    lo = mid
                else:
                    hi = mid
            released = lo
            pr = pressure(local, released)
        custody = releasable - released
        relative = pr['price'] / reference
        previous_relative = relative
        residual = opening + available + draw + topup - service - contribution - r['capital'] - dividend - cash
        if abs(residual) > max(1.0, 1e-12 * max(opening, r['capital'])):
            raise ValueError('Cash identity failed')
        out.append(dict(year=y, contribution=topup, topup_by_purpose=purposes, pre_topup_cash=before_topup, cash=cash, reserves=reserve, reserve_contribution=contribution, retained_requirement=retention, dividend=dividend, owner_payout=dividend - intended + released, owner_custody=custody, local_payout=released, capital=r['capital'], operating_available=available, gross=gross, draw=draw, external_draw=external, supplier_U_draw=supplier, local_U_draw=local, debt=sum((l['remaining'] for l in loans)), U_debt=sum((l['remaining'] for l in loans if l['kind'] == 'U')), interest=due['interest'], principal=due['principal'], local_U_interest=due['U_interest'] - due['supplier_interest'], supplier_interest=due['supplier_interest'], external_interest=due['external_interest'], covenant_DSCR=dscr, covenant_lock=lock, price=pr['price'], relative_price=relative, price_breach=relative > ceiling + 1e-10, excess_credit_balances=pr['excess_credit_balances'], cash_residual=residual))
    first = next((r['year'] for r in out if r['contribution'] > 1.0), None)
    return dict(mode=mode, initial_reserve=initial_reserve, first_unfunded_year=first, self_financing_prescribed_schedule=first is None, annual=out)

# From research/studies/enclosed_returns/finance.py at e49fca09, metrics.
def metrics(path, sim, cfg, horizon):
    n = horizon + 1
    start = path['annual'][0]['year']
    sponsor = [0.0] * n
    unlevered = [0.0] * n
    exports = [0.0] * n
    resident = [0.0] * n
    costs = [0.0] * n
    for y in range(start):
        sponsor[y] = -path['seed_hardware_and_carry'] / start
        unlevered[y] = sponsor[y]
        costs[y] = (path['seed_hardware_and_carry'] - path['seed_resource_provision'] - path['seed_resident_payroll']) / start
        # Half occupancy is the same seed-work assumption used in provisioning.
        costs[y] += path['seed_population'] / 2 * cfg['resident_work_fraction'] * cfg['opportunity_USD2041_worker_year']
    sponsor[start - 1] -= path['initial_working']
    sponsor[0] -= sim['initial_reserve']
    records = []
    for r, f in zip(path['annual'], sim['annual']):
        y = r['year']
        if y > horizon:
            break
        records.append(f)
        sponsor[y] = f['owner_payout'] - f['contribution']
        unlevered[y] = f['operating_available'] - r['capital']
        exports[y] = f['gross']
        resident[y] = r['resident_value']
        opportunity = r['population'] * cfg['resident_work_fraction'] * cfg['opportunity_USD2041_worker_year']
        costs[y] = r['operating_cost'] + r['capital'] - r['resource_provision'] - r['resident_payroll'] + opportunity
    last = records[-1]
    end = next((r for r in path['annual'] if r['year'] == horizon))
    terminal = last['cash'] + last['reserves'] + last['owner_custody'] - last['debt']
    discounted = {}
    for rate in cfg['rates']:
        S = pv(sponsor, rate)
        E = pv(exports, rate)
        R = pv(resident, rate)
        C = pv(costs, rate)
        annuity = math.fsum(((1 + rate) ** (-y) for y in range(start, horizon + 1)))
        discounted[str(rate)] = dict(sponsor_distributed_cash_NPV=S, sponsor_cash_and_terminal_NPV=S + terminal / (1 + rate) ** horizon, terminal_financial_claim_PV=terminal / (1 + rate) ** horizon, unlevered_cash_NPV=pv(unlevered, rate), full_price_socioeconomic_NPV=E + R - C, half_price_socioeconomic_NPV=0.5 * (E + R) - C, export_only_full_price_NPV=E - C, full_price_benefit_cost_ratio=(E + R) / C, export_value_break_even_fraction=(C - R) / E if E else None, additional_annual_net_social_benefit=max(0.0, C - E - R) / annuity, additional_initial_cost_at_unlevered_break_even=max(0.0, pv(unlevered, rate)) / sum(((1 + rate) ** (-y) / start for y in range(start))), discounted_exports=E, discounted_resident_value=R, discounted_resource_cost=C, durable_distributed_cash_recovery_year=durable_recovery(sponsor, rate), valuation_sensitivity=[dict(export_fraction=x, resident_fraction=y, net_value=x * E + y * R - C) for x in (0.0, 0.5, 1.0) for y in (0.0, 0.5, 1.0)], local_interest_relief_and_bank_margin_loss=math.fsum((f['local_U_interest'] * (cfg['external_rate'] - cfg['local_U_rate']) / cfg['local_U_rate'] / (1 + rate) ** f['year'] for f in records)))
    funding = sum((f['contribution'] for f in records))
    first = next((f['year'] for f in records if f['contribution'] > 1), None)
    return dict(pooled_funder_IRR_if_unique=nonnegative_irr(sponsor), name=path['name'], portfolio=path['portfolio'], market=path['market'], growth=path['growth'], mode=sim['mode'], horizon=horizon, seed_cost_scope=path.get('seed_cost_scope', ''), unpriced_deployment_requirements=path.get('unpriced_deployment_requirements', []), population=end['population'], service_kits=end['service_kits'], seed=path['seed'], seed_direct_cost=path.get('seed_direct_cost', path['seed_hardware_and_carry']), seed_capital_cost_allowance=path.get('seed_capital_cost_allowance', 0.0), initial_reserve=sim['initial_reserve'], no_additional_equity_to_checkpoint=first is None, credit_supported_self_financing_to_checkpoint=first is None, first_unfunded_year=first, additional_equity_required=funding, equity_by_purpose={k: sum((f['topup_by_purpose'][k] for f in records)) for k in ('operating', 'debt', 'capital', 'liquidity')}, capital_orders=sum((r['capital'] for r in path['annual'] if r['year'] <= horizon)), owner_payouts=sum((f['owner_payout'] for f in records)), peak_debt=max((f['debt'] for f in records)), cumulative_U_draws=sum((f['local_U_draw'] + f['supplier_U_draw'] for f in records)), cumulative_external_draws=sum((f['external_draw'] for f in records)), minimum_covenant_DSCR=min((f['covenant_DSCR'] for f in records if f['covenant_DSCR'] is not None), default=None), price_breach_years=[f['year'] for f in records if f['price_breach']], maximum_relative_price=max((f['relative_price'] for f in records)), annual_gross=last['gross'], annual_operating_available=last['operating_available'], annual_capital=end['capital'], annual_post_capital_cash=last['operating_available'] - end['capital'], annual_operating_break_even_service_gross_at_half_margin=max(0.0, -last['operating_available']) / 0.5, annual_build_and_debt_break_even_service_gross_at_half_margin=max(0.0, end['capital'] + last['interest'] + last['principal'] + last['reserve_contribution'] - last['operating_available']) / 0.5, revenue_by_sector={k: v * last['gross'] / end['gross'] if end['gross'] else 0.0 for k, v in end['revenue_by_sector'].items()}, terminal_cash=last['cash'], terminal_reserves=last['reserves'], terminal_debt=last['debt'], terminal_custody=last['owner_custody'], discounted=discounted)

# From research/studies/enclosed_returns/revised.py at e49fca09, sponsor_cash.
def sponsor_cash(path, sim, horizon):
    start = path['annual'][0]['year']
    cash = [0.0] * (horizon + 1)
    for y in range(min(start, horizon + 1)):
        cash[y] = -path['seed_hardware_and_carry'] / start
    cash[0] -= sim['initial_reserve']
    if start - 1 <= horizon:
        cash[start - 1] -= path['initial_working']
    for r in sim['annual']:
        if r['year'] <= horizon:
            cash[r['year']] = r['owner_payout'] - r['contribution']
    return cash

# From research/studies/enclosed_returns/revised.py at e49fca09, measures.
def measures(path, sim, cfg, horizon):
    m = metrics(path, sim, cfg, horizon)
    rows = [r for r in sim['annual'] if r['year'] <= horizon]
    initial = path['seed'] + sim['initial_reserve']
    contributions = initial + m['additional_equity_required']
    cf = sponsor_cash(path, sim, horizon)
    m.update(initial_full_funding=initial, all_funder_contributions=contributions, cumulative_gross_receipts=sum((r['gross'] for r in rows)), undiscounted_payout_multiple=m['owner_payouts'] / contributions, first_payout_year=next((r['year'] for r in rows if r['owner_payout'] > 1), None), real_IRR_roundoff_screened=nonnegative_irr([0.0 if abs(v) < 1 else v for v in cf]), sustaining_years_with_negative_operating_cash=[r['year'] for r in rows if r['year'] >= 200 and r['operating_available'] < 0], evidence='Conditional funded-completion account; additional equity, if required, must be supplied. Cash distributions, retained financial claims and physical assets are separate valuations.')
    return m
