"""Bounded U funding-policy overlays, with a dated receivables facility ledger.

The pinned credit producer remains unchanged. Interest repricing holds its real
investment schedule fixed; wider supplier acceptance reruns that same producer.
Monthly payment stress pauses discretionary uses and conserves invoice receipts.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import json
import math
from pathlib import Path

from . import credit as c

HERE = c.HERE


def inputs():
    cfg = json.loads((HERE / 'policy_inputs.json').read_text())
    path = HERE / 'results/credit.json'
    if c.s.run.digest(path) != cfg['credit_result_sha256']:
        raise ValueError('Pinned credit result changed')
    saved = json.loads(path.read_text())
    if saved['producer_sha256'] != c.s.run.digest(HERE / 'credit.py'):
        raise ValueError('Credit producer changed')
    for name, sha in saved['dependencies'].items():
        if c.s.run.digest(HERE / name) != sha:
            raise ValueError('Credit dependency changed: ' + name)
    return cfg, saved, json.loads((HERE / 'credit_inputs.json').read_text())


def rate_relief(sim, horizon, reduction):
    rows = [r for r in sim['annual'] if r['year'] <= horizon and 'failure' not in r]
    # Existing contracts are not retroactively rewritten: this is an alternative
    # offer at each actual origination date, preserving its principal schedule.
    savings = [(r['year'], max(0., r['U_interest'] - r['supplier_interest']) * reduction / .03) for r in rows]
    nominal = math.fsum(value for _, value in savings)
    pv = math.fsum(value / 1.03 ** year for year, value in savings)
    return dict(rate_reduction=reduction, realized_borrower_interest_saving=nominal,
                borrower_interest_saving_NPV_3_percent=pv,
                bank_gross_income_reduction_NPV_3_percent=pv,
                consolidated_interest_transfer_NPV=0.,
                operating_less_capital_NPV_change=0.,
                affected_portion_coverage_capacity_increase=.08/(.08-reduction)-1)


def snapshot(row, profile, mon, seed):
    basket = (profile['domestic']['USD2041']['resident_consumption'] +
              profile['domestic']['USD2041']['domestic_capital_formation']) * row['scale']
    return dict(revenue=row['gross'], operating_cost=row['gross']-row['available'],
                core_imports=row['core_imports'], cash=max(0., row['cash']-row['FX_restricted']),
                debt_service=row['interest']+row['principal'],
                external_debt_service=row['external_interest']+row['principal']-row['U_principal'],
                reserves=row['reserves'], existing_U_loans=row['U_debt'],
                existing_convertible_U=max(0., row['U_debt']-row['outside_U']),
                basket=basket, price=row['price'], relative_price=row['relative_price'],
                holding_years=mon['domestic_holding_days']/365,
                elasticity=mon['supply_elasticity'], cumulative_limit=mon['cumulative_relative_price_limit'],
                seed=seed, baseline_annual_payout=row['cash_distribution_to_owners'],
                baseline_annual_new_orders=row['industrial_order']+row['service_order'])


def balance_allowance(snap, cfg):
    ratio = max(1., min(1+cfg['maximum_incremental_price_premium'],
                         (1+snap['cumulative_limit'])/snap['relative_price']))
    return snap['holding_years']*snap['basket']*snap['price']**(1+snap['elasticity'])*(ratio**(1+snap['elasticity'])-1)


def bridge(snap, cfg, policy, months, *, trace=False):
    """Independent monthly episode; all values in constant USD2041 equivalents.

    A bank starts with FX assets R, old U loans D, old liabilities D and equity R.
    A bridge draw adds a loan and claims; FX conversion retires the corresponding
    claims and spends R. Repayment restores FX or retires U claims pro rata.
    Term redemption uses matched external saving assets: no unlocked payout enters
    local circulation. Costs are accrued and settled at episode end, separately.
    """
    cash = snap['cash']; reserve = snap['reserves']; accessible = reserve*policy['FX_access']
    floor = cfg['future_core_import_floor_years']*snap['core_imports']
    cap = cfg['facility_fraction_of_seed']*snap['seed']
    allowance = balance_allowance(snap, cfg)
    ext_cost = snap['core_imports'] + cfg['unclassified_operating_external_fraction']*max(0., snap['operating_cost']-snap['core_imports'])
    obligation = (snap['operating_cost']+snap['debt_service'])/12
    external_debt_obligation = snap['external_debt_service']/12
    revenue = snap['revenue']/12
    # Four principal accounts: converted into outside money, outside U term
    # holdings, spendable local U, and locked local U saving claims.
    parts = dict(FX=0., outside_term=0., local=0., local_term=0.)
    queue = {}; receivable = collected = invoiced = paid = 0.
    fee = interest = term_interest = draws = repayments = 0.
    peak = peak_terms = peak_local = peak_conversion = 0.
    peak_local_terms = peak_outside_terms = term_redemptions = 0.
    peak_receivable = peak_at_receivable = 0.
    min_FX = accessible; max_identity = 0.; active_months = 0; first_draw = last_repay = None
    rows = []; fail = None; peak_parts = dict(parts)
    external = policy.get('external_funded', False)
    rate = cfg['external_bridge_rate'] if external else cfg['bridge_rate']
    total_months = 2*months + cfg['collection_tail_months']
    for month in range(1, total_months+1):
        opening_loan = sum(parts.values())
        interest += opening_loan*rate/12
        term_interest += (parts['outside_term']+parts['local_term'])*cfg['term_claim_rate']/12
        delayed = cfg['delayed_fraction']*revenue if month <= months else 0.
        due = queue.pop(month, 0.)
        receipts = revenue-delayed+due
        if delayed: queue[month+months] = queue.get(month+months, 0.)+delayed
        receivable += delayed-due; cash += receipts
        invoiced += revenue; collected += receipts
        peak_receivable = max(peak_receivable, receivable)
        need = max(0., obligation-cash)
        # Cash is assumed convertible at the inherited reference valuation.
        # It pays the external bill first; escrow is excluded from initial cash.
        missing_external_debt = max(0., external_debt_obligation-cash)
        missing_external_suppliers = max(0., ext_cost/12-max(0.,cash-external_debt_obligation))
        missing_external = min(need, missing_external_debt+missing_external_suppliers)
        missing_local = need-missing_external
        add = dict(FX=missing_external_debt+missing_external_suppliers*(1-policy['outside_U_share']),
                   outside_term=missing_external_suppliers*policy['outside_U_share'],
                   local=missing_local, local_term=0.)
        if external: add = dict(FX=need, outside_term=0., local=0., local_term=0.)
        if policy['term_absorption']:
            locked = max(0., parts['local']+add['local']-allowance)
            add['local'] -= locked; add['local_term'] += locked
        candidate = {k: parts[k]+add[k] for k in parts}
        loan = sum(candidate.values())
        reasons = []
        if need > 1:
            if not policy['facility']: reasons.append('no_facility')
            if snap['revenue'] <= snap['operating_cost']+snap['debt_service']:
                reasons.append('nonpositive_operating_coverage')
            if loan > (1-cfg['collateral_haircut'])*receivable+1: reasons.append('collateral')
            if loan > cap+1: reasons.append('facility_limit')
            if not external and policy['facility']:
                liquidity_floor = floor+cfg['conversion_buffer_fraction']*(snap['existing_convertible_U']+candidate['local'])
                if accessible-add['FX'] < liquidity_floor-1: reasons.append('external_liquidity')
                if reserve < cfg['capital_fraction_of_loans']*(snap['existing_U_loans']+loan)-1: reasons.append('capital')
                if candidate['local'] > allowance+1: reasons.append('price_balance_allowance')
        if reasons:
            fail = dict(month=month, payment_gap_before_partial_facility_draw=need, reasons=reasons,
                        required_total_loan=loan, receivables=receivable)
            break
        if need:
            parts = candidate; draws += need; fee += cfg['origination_fee']*need
            if not external: accessible -= add['FX']
            cash += need
            if first_draw is None: first_draw = month
        cash -= obligation; paid += obligation
        loan = sum(parts.values())
        if loan > peak:
            peak = loan; peak_parts = dict(parts); peak_at_receivable = receivable
        peak_terms = max(peak_terms, parts['local_term']+parts['outside_term'])
        peak_local_terms = max(peak_local_terms, parts['local_term'])
        peak_outside_terms = max(peak_outside_terms, parts['outside_term'])
        peak_local = max(peak_local, parts['local'])
        peak_conversion = max(peak_conversion, parts['FX'])
        if loan > 1: active_months += 1
        # Receipts above current operating bills sweep into principal repayment.
        repayment = min(cash, loan)
        if repayment > 0 and loan > 0:
            factor = repayment/loan
            if not external: accessible += parts['FX']*factor
            term_redemptions += (parts['outside_term']+parts['local_term'])*factor
            parts = {k: max(0., value*(1-factor)) for k, value in parts.items()}
            cash -= repayment; repayments += repayment; last_repay = month
        min_FX = min(min_FX, accessible)
        remaining = sum(parts.values())
        liabilities = remaining-parts['FX'] if not external else 0.
        fx_assets = reserve-parts['FX'] if not external else reserve
        identity = fx_assets+snap['existing_U_loans']+(remaining if not external else 0.)-snap['existing_U_loans']-liabilities-reserve
        max_identity = max(max_identity, abs(identity))
        cash_identity = snap['cash']+collected+draws-paid-repayments-cash
        max_identity = max(max_identity, abs(cash_identity), abs(invoiced-collected-receivable))
        if trace:
            rows.append(dict(month=month, receipts=receipts, receivables=receivable, cash=cash,
                             loan=remaining, accessible_FX=accessible, **parts))
    costs = interest+fee
    terminal_loan = sum(parts.values())
    if not fail and cash+1 < costs:
        fail = dict(month=total_months, unpaid_obligation=costs-cash, reasons=['financing_cost'])
    if not fail and terminal_loan > 1:
        fail = dict(month=total_months, unpaid_obligation=terminal_loan, reasons=['facility_not_repaid'])
    closing_cash = cash-costs if not fail else cash
    payout_exposure = snap['baseline_annual_payout']*total_months/12
    payout_delay_pv = (sum(snap['baseline_annual_payout']/12/1.03**(m/12)
                           for m in range(1,total_months+1))
                       - payout_exposure/1.03**(total_months/12))
    stress = []
    if peak > 0:
        for loss in cfg['loss_fractions']:
            recovery = min(peak, (1-loss)*peak_at_receivable)
            credit_loss = peak-recovery
            stress.append(dict(receivable_loss_fraction=loss, loan_loss=credit_loss,
                               U_fund_loan_loss=0. if external else credit_loss,
                               external_lender_loan_loss=credit_loss if external else 0.,
                               fund_equity_after_loss=reserve-(0. if external else credit_loss),
                               recapitalization_to_15_percent=max(0.,cfg['capital_fraction_of_loans']*snap['existing_U_loans']-(reserve-(0. if external else credit_loss))),
                               issued_U_claims_before_resolution=0. if external else peak-peak_parts['FX']))
    return dict(status='payment_shortfall' if fail else 'paid_and_repaid', failure=fail,
                opening_usable_cash=snap['cash'], opening_reserve_assets=reserve,
                annual_operating_cost=snap['operating_cost'], local_balance_allowance=allowance,
                peak_bridge=peak, gross_draws=draws, principal_repaid=repayments,
                peak_FX_use=peak_conversion, peak_term_claims=peak_terms, peak_local_balances=peak_local,
                peak_local_term_claims=peak_local_terms, peak_outside_term_claims=peak_outside_terms,
                term_principal_redeemed_into_external_savings=term_redemptions,
                peak_receivables=peak_receivable, minimum_accessible_FX=min_FX,
                bridge_active_months=active_months, first_draw_month=first_draw, final_repayment_month=last_repay,
                borrower_interest=interest, origination_fees=fee, borrower_financing_cost=costs,
                term_holder_interest=term_interest, intermediary_margin_before_costs=costs-term_interest,
                terminal_bridge=terminal_loan, terminal_receivables=receivable, closing_cash=closing_cash,
                max_accounting_residual=max_identity, collateral_loss_stress=stress,
                cashflow_freeze_months=total_months,
                baseline_payouts_deferred_over_episode=payout_exposure if not fail else None,
                baseline_orders_paused_over_episode=snap['baseline_annual_new_orders']*total_months/12 if not fail else None,
                payout_delay_PV_if_released_at_end_3_percent=payout_delay_pv if not fail else None,
                # Financial cost is conditional on this one episode, not expected loss.
                paused_baseline_payout_rate=snap['baseline_annual_payout'],
                paused_baseline_new_order_rate=snap['baseline_annual_new_orders'], rows=rows)


def product():
    cfg, saved, base = inputs(); q,d,b,e,t = c.s.inputs()
    funding = []; liquidity = []; sensitivity = []
    old = {(r['strategy'],r['market'],r['horizon']): r for r in saved['horizons'] if r['mode']=='international_U_mix'}
    for name in q['strategies']:
        for market in q['markets']:
            control = c.simulate(q,d,b,e,t,name,market,max(cfg['horizons']),base,'retained_earnings')
            sim = c.simulate(q,d,b,e,t,name,market,max(cfg['horizons']),base,'international_U_mix',control)
            wider = deepcopy(base)
            wider['financing']['international_U_mix']['outside_supplier_U_fraction'] = cfg['wider_capital_supplier_U_fraction']
            alt = c.simulate(q,d,b,e,t,name,market,max(cfg['horizons']),wider,'international_U_mix',control)
            rows = {r['year']:r for r in sim['annual'] if 'failure' not in r}
            for horizon in cfg['horizons']:
                ref = c.summary(sim,horizon,base,control)
                if ref != old[name,market,horizon]: raise ValueError('Credit baseline failed reproduction')
                new = c.summary(alt,horizon,wider,control)
                valid = ref['status']==new['status']=='continuing_conditional_path'
                funding.append(dict(strategy=name,market=market,horizon=horizon,baseline=ref,wider_supplier_U=new,
                    rate_relief=[rate_relief(sim,horizon,v) for v in cfg['local_rate_reductions']],
                    wider_rate_relief=[rate_relief(alt,horizon,v) for v in cfg['local_rate_reductions']],
                    wider_incremental_operating_less_capital_NPV=(new['incremental_operating_less_capital_NPV_at_3_percent']-ref['incremental_operating_less_capital_NPV_at_3_percent']) if valid and new['matched_comparison_complete'] and ref['matched_comparison_complete'] else None))
                if horizon not in rows:
                    liquidity.append(dict(strategy=name,market=market,horizon=horizon,status='unavailable_after_baseline_failure'))
                    continue
                r = rows[horizon]
                row = dict(r,gross=r['operation']['gross_USD2041'],available=r['operation']['available_USD2041'],core_imports=r['operation']['core_imports_USD2041'])
                snap = snapshot(row,c.s.profile(q,d,name),base['monetary'],sim['external_seed'])
                for months in cfg['liquidity']['outage_months']:
                    for policy, spec in cfg['policies'].items():
                        result = bridge(snap,cfg['liquidity'],spec,months)
                        liquidity.append(dict(strategy=name,market=market,horizon=horizon,months=months,policy=policy,**result))
                if market=='breakthrough' and horizon in (100,800):
                    for label, edits in [('all_unclassified_local',{'unclassified_operating_external_fraction':0.}),('all_unclassified_external',{'unclassified_operating_external_fraction':1.}),('half_invoices_delayed',{'delayed_fraction':.5}),('half_advance',{'collateral_haircut':.5}),('full_year_import_floor',{'future_core_import_floor_years':1.})]:
                        local = dict(cfg['liquidity'],**edits)
                        sensitivity.append(dict(strategy=name,market=market,horizon=horizon,variant=label,months=6,policy='U_term_50',**bridge(snap,local,cfg['policies']['U_term_50'],6)))
    dependencies = ['policy_inputs.json','credit.py','results/credit.json',*saved['dependencies']]
    return dict(schema='terluna.research.u-policy/1',evidence=cfg['evidence'],boundary=saved['boundary'],
                producer_sha256=c.s.run.digest(HERE/'policy.py'),dependencies={p:c.s.run.digest(HERE/p) for p in dependencies},
                reading_rule='Fixed investment rate relief is a borrower/bank transfer. Supplier acceptance can alter investment timing. Monthly liquidity episodes freeze activity at a dated checkpoint, pause discretionary uses, and require voluntary term holdings and actual outside settlement assets. No CPI or default probability is inferred.',
                funding=funding,liquidity=liquidity,sensitivities=sensitivity)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=c.s.run.ROOT/'research/runs/complex_returns/policy.json')
    args=parser.parse_args();args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(product(),indent=2,allow_nan=False)+'\n')
    print(args.output)
