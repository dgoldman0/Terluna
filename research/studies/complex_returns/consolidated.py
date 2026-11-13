"""Dated financial and socioeconomic return accounts for the latest U package.

Consumes unchanged investment producers. The new saving-release assumption is an
explicit incidence overlay; it does not reinvest interest savings or create output.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import json
import math
from pathlib import Path

from . import credit as c
from . import policy as p

HERE = c.HERE


def inputs():
    cfg = json.loads((HERE/'consolidated_inputs.json').read_text())
    path = HERE/'results/policy.json'
    if c.s.run.digest(path) != cfg['policy_result_sha256']:
        raise ValueError('Pinned policy result changed')
    policy = json.loads(path.read_text())
    if c.s.run.digest(HERE/'policy.py') != policy['producer_sha256']:
        raise ValueError('Policy producer changed')
    for name, sha in policy['dependencies'].items():
        if c.s.run.digest(HERE/name) != sha:
            raise ValueError('Policy dependency changed: '+name)
    pcfg, saved, base = p.inputs()
    return cfg, policy, pcfg, saved, base


def pv(values, rate):
    return math.fsum(v*(1+rate)**(-t) for t,v in enumerate(values))


def durable_recovery(flows, rate=0.):
    cumulative = []; balance = 0.
    for year, value in enumerate(flows):
        balance += value/(1+rate)**year; cumulative.append(balance)
    running = math.inf; result = None
    for year in reversed(range(len(cumulative))):
        running = min(running,cumulative[year])
        if running >= 0: result = year
    return result


def nonnegative_irr(flows):
    signs = [1 if x>0 else -1 for x in flows if x]
    changes = sum(a!=b for a,b in zip(signs,signs[1:]))
    if not signs or signs[0]!=-1 or changes!=1 or pv(flows,0.)<0:
        return None
    lo, hi = 0., 1.
    while pv(flows,hi)>0:
        hi *= 2
        if hi > 1024: return None
    # All actual roots fall below100%; avoid high powers at long horizons.
    for _ in range(70):
        mid = (lo+hi)/2
        if pv(flows,mid)>0: lo=mid
        else: hi=mid
    return (lo+hi)/2


def cashflows(sim, horizon, bootstrap, reduction=0.):
    n = horizon+1
    initial = [0.]*n; dividends = [0.]*n; relief = [0.]*n
    bank = [0.]*n; external = [0.]*n; suppliers = [0.]*n
    for year in range(min(bootstrap,n)):
        initial[year] = sim['initial_hardware_and_carry']/bootstrap
    initial[bootstrap-1] += sim['initial_working_cash']
    for row in sim['annual']:
        year = row['year']
        if year>horizon or 'failure' in row: continue
        local_interest = max(0.,row['U_interest']-row['supplier_interest'])
        relief[year] = local_interest*reduction/.03
        dividends[year] = row['cash_distribution_to_owners']
        bank[year] = local_interest-relief[year]
        external[year] = row['external_interest']+row['principal']-row['U_principal']-row['external_draw']
        suppliers[year] = row['supplier_interest']+row['supplier_principal']-row['supplier_credit_origination']
    equity = [a+b-z for a,b,z in zip(dividends,relief,initial)]
    sponsor = equity[:]; sponsor[0] -= sim['initial_reserve']
    return dict(initial=initial,dividends=dividends,relief=relief,bank_margin=bank,
                external_lender=external,supplier_claims=suppliers,equity=equity,sponsor=sponsor)


def embedded_resident_payroll(q,b,profile,name,row):
    scale=row['scale']; kits=row['service_kits']
    staff = b['cultural_platform']['staff']+kits*b['service_kit']['staff']
    staff += q['strategies'][name]['premium_share']*scale*q['premium']['staff_per_founding_scale']
    staff += bool(profile['repricing']['propellant_kg_year'])*scale*q['logistics']['support_staff_per_founding_scale']
    return staff*b['service_staff_USD2041_year']


def capital_reserve(profile,b,committed,starts):
    first = min(starts,max(0.,1-committed))
    commission = profile['commissioning']['additional_resource_provision_USD2041']
    g=profile['remaining_balanced_growth']
    growth=(g['at_limit']['external']['USD2041']['resource_provision']-profile['external']['USD2041']['resource_provision'])/g['fraction_of_residents_year']
    return b['marginal_capital_multiplier']*(first*commission+(starts-first)*(commission+growth))


def economic_components(sim,q,b,profile,name,horizon,bootstrap,old_inputs):
    n=horizon+1; out={k:[0.]*n for k in ('exports','resident_services','procurement_cost','resident_opportunity_cost','removed_resource_provision','removed_resident_payroll')}
    seed_reserve = b['seed_scale']*profile['commissioning']['additional_resource_provision_USD2041']
    bootstrap_reserve = b['seed_scale']*profile['external']['USD2041']['resource_provision']*b['bootstrap_carrying_fraction']*bootstrap
    initial = sim['initial_hardware_and_carry']-seed_reserve-bootstrap_reserve
    for year in range(bootstrap):
        out['procurement_cost'][year]=initial/bootstrap
        out['removed_resource_provision'][year]=(seed_reserve+bootstrap_reserve)/bootstrap
    committed=b['seed_scale']
    for row in sim['annual']:
        year=row['year']
        if year>horizon or 'failure' in row: continue
        op=row['operation']; wages=embedded_resident_payroll(q,b,profile,name,row)
        reserve=op['resource_provision_USD2041']+capital_reserve(profile,b,committed,row['starts'])
        out['exports'][year]=op['gross_USD2041']
        out['resident_services'][year]=profile['domestic']['USD2041']['resident_consumption']*row['scale']
        out['procurement_cost'][year]=op['total_cost_USD2041']+row['industrial_order']+row['service_order']-reserve-wages
        out['resident_opportunity_cost'][year]=profile['engineering']['population']*row['scale']*old_inputs['resident_work_fraction']*old_inputs['alternative_USD2041_worker_year']
        out['removed_resource_provision'][year]=reserve
        out['removed_resident_payroll'][year]=wages
        committed+=row['starts']
    return out


def financial_metrics(sim,horizon,bootstrap,rates,reduction=0.):
    f=cashflows(sim,horizon,bootstrap,reduction)
    rows=[r for r in sim['annual'] if r['year']<=horizon];last=rows[-1]
    failed='failure' in last
    terminal=None if failed else dict(
        net_project_liquid_assets=last['cash']+last['U_escrow']-last['supplier_prepay']-last['debt'],
        owner_custody=last['custody'],monetary_reserve_assets=last['reserves'],
        debt=last['debt'],supplier_prepayments_excluded=last['supplier_prepay'])
    discounted={}
    for rate in rates:
        cost=pv(f['initial'],rate)+sim['initial_reserve']
        paid=pv(f['dividends'],rate);relief=pv(f['relief'],rate)
        sponsor=pv(f['sponsor'],rate)
        terminal_value=None if terminal is None else sum(terminal[k] for k in ('net_project_liquid_assets','owner_custody','monetary_reserve_assets'))/(1+rate)**horizon
        factor=sum((1+rate)**(-year) for year in range(bootstrap,horizon+1))
        discounted[str(rate)]=dict(project_equity_cash_NPV=pv(f['equity'],rate),
            sponsor_contractual_cash_NPV=sponsor-relief,sponsor_cash_NPV_with_saving_release=sponsor,
            initial_sponsor_funding_PV=cost,payout_PV=paid,conditional_saving_release_PV=relief,
            payout_to_funding_multiple=(paid+relief)/cost,
            U_bank_gross_margin_PV=pv(f['bank_margin'],rate),
            U_bank_margin_reduction_PV=relief,external_lender_scheduled_cash_NPV=pv(f['external_lender'],rate),
            supplier_term_claim_scheduled_cash_NPV=pv(f['supplier_claims'],rate),
            terminal_book_claim_PV=terminal_value,
            sponsor_NPV_including_terminal_book_claims=None if terminal_value is None else sponsor+terminal_value,
            additional_annual_sponsor_cash_to_close_gap=None if failed else max(0.,-sponsor)/factor,
            discounted_recovery_year=durable_recovery(f['sponsor'],rate))
    return dict(status='failed_path_cash_to_failure' if failed else 'continuing_conditional_path',
                failure=last.get('failure'),original_seed=sim['external_seed'],initial_monetary_reserve=sim['initial_reserve'],
                first_productive_year=bootstrap,terminal_book_claims=terminal,
                cumulative_owner_payout=sum(f['dividends']),cumulative_saving_release=sum(f['relief']),
                sponsor_nondiscounted_cash_return=sum(f['sponsor']),
                sponsor_nonnegative_IRR=nonnegative_irr(f['sponsor']),
                undiscounted_recovery_year=durable_recovery(f['sponsor']),discounted=discounted)


def economic_metrics(components,rate,years,cfg):
    v={key:pv(value,rate) for key,value in components.items()}
    factor=sum(1/(1+rate)**year for year in years)
    valuations=[]
    for exports in cfg['economic_export_value_fractions']:
        for residents in cfg['resident_service_value_fractions']:
            benefit=exports*v['exports']+residents*v['resident_services']
            cost=v['procurement_cost']+v['resident_opportunity_cost']
            valuations.append(dict(export_value_fraction=exports,resident_value_fraction=residents,
                net_value=benefit-cost,benefit_cost_ratio=benefit/cost,
                additional_constant_net_benefit_to_close_gap=max(0.,cost-benefit)/factor))
    labour=[]
    for multiplier in cfg['resident_opportunity_cost_multipliers']:
        cost=v['procurement_cost']+multiplier*v['resident_opportunity_cost']
        labour.append(dict(opportunity_cost_multiplier=multiplier,
            full_price_proxy_net_value=v['exports']+v['resident_services']-cost))
    return dict(discounted_components=v,constant_annual_benefit_factor=factor,
                export_value_break_even_with_full_resident_value=(v['procurement_cost']+v['resident_opportunity_cost']-v['resident_services'])/v['exports'],
                valuation_grid=valuations,labour_sensitivity=labour)


def product():
    cfg,policy,pcfg,saved,base=inputs();q,d,b,e,t=c.s.inputs(); old,_,_=c.s.run.inputs()
    saved_credit={(r['strategy'],r['market'],r['mode'],r['horizon']):r for r in saved['horizons']}
    saved_policy={(r['strategy'],r['market'],r['horizon']):r for r in policy['funding']}
    historical=json.loads((HERE/'results/strategies.json').read_text())
    history={(r['strategy'],r['market'],r['horizon_elapsed_years']):r for r in historical['horizons']}
    financial=[];economic=[];reconciliation=[];liquidity=[];dated=[]
    for name in q['strategies']:
        profile=c.s.profile(q,d,name)
        for market in q['markets']:
            control=c.simulate(q,d,b,e,t,name,market,max(cfg['horizons']),base,'retained_earnings')
            simulations={'retained_earnings':control}
            for mode in ('mixed_credit','international_U_mix'):
                simulations[mode]=c.simulate(q,d,b,e,t,name,market,max(cfg['horizons']),base,mode,control)
            wider=deepcopy(base);wider['financing']['international_U_mix']['outside_supplier_U_fraction']=cfg['latest_capital_supplier_U_share']
            simulations['wider_U']=c.simulate(q,d,b,e,t,name,market,max(cfg['horizons']),wider,'international_U_mix',control)
            latest=simulations['wider_U'];simulations['latest_U_package']=latest
            for mode,sim in simulations.items():
                reduction=cfg['local_U_rate_reduction'] if mode=='latest_U_package' else 0.
                for horizon in cfg['horizons']:
                    source=c.summary(sim,horizon,wider if mode in ('wider_U','latest_U_package') else base,control if mode!='retained_earnings' else None)
                    expected=(saved_policy[name,market,horizon]['wider_supplier_U'] if mode in ('wider_U','latest_U_package') else saved_credit[name,market,mode,horizon])
                    if source!=expected: raise ValueError('Source financial reproduction failed')
                    metric=financial_metrics(sim,horizon,b['bootstrap_years'],cfg['real_discount_rates'],reduction)
                    if not math.isclose(metric['discounted']['0.03']['sponsor_contractual_cash_NPV'],source['project_and_reserve_funders_cash_NPV']['0.03'],rel_tol=1e-12,abs_tol=.01):
                        raise ValueError('Cash payout reconciliation failed')
                    end=next(r for r in reversed(sim['annual']) if r['year']<=horizon)
                    revenue=end['operation']['revenue_by_sector_USD2041']; services=sum(value for key,value in revenue.items() if key not in ('materials','source_products','propellant'))
                    financial.append(dict(strategy=name,market=market,mode=mode,horizon=horizon,**metric,
                        final_population=profile['engineering']['population']*end['scale'],
                        first_founding_scale_year=source['first_founding_scale_year'],
                        final_scale=end['scale'],final_service_kits=end['service_kits'],
                        annual_receipts=end['operation']['gross_USD2041'],annual_available_before_finance=end['operation']['available_USD2041'],
                        nonmaterial_and_service_receipt_fraction=services/end['operation']['gross_USD2041'] if end['operation']['gross_USD2041'] else 0.,
                        revenue_by_sector=revenue,cumulative_capital_orders=source['cumulative_capital_orders'],
                        peak_debt=source['peak_debt'],minimum_DSCR=source['minimum_DSCR'],
                        max_price_premium_percent=source['max_relative_price_change_percent'],price_breaches=source['price_limit_breach_years'],
                        incremental_operating_less_capital_NPV=source.get('incremental_operating_less_capital_NPV_at_3_percent',0.)))
            latest_rows={r['year']:r for r in latest['annual'] if 'failure' not in r}
            for horizon in cfg['horizons']:
                by_mode={r['mode']:r for r in financial if (r['strategy'],r['market'],r['horizon'])==(name,market,horizon)}
                orig=history[name,market,horizon]
                at3=lambda mode:by_mode[mode]['discounted']['0.03']['sponsor_cash_NPV_with_saving_release']
                historical_distributions=orig['NPV_USD2041']['0.03']-orig['terminal_liquid_USD2041']/1.03**horizon
                reconstruction=dict(strategy=name,market=market,horizon=horizon,
                    original_distribution_NPV=historical_distributions,
                    initial_monetary_reserve=-control['initial_reserve'],
                    matched_controller_change=at3('retained_earnings')-(historical_distributions-control['initial_reserve']),
                    commercial_U25_change=at3('international_U_mix')-at3('retained_earnings'),
                    wider_supplier_change=at3('wider_U')-at3('international_U_mix'),
                    conditional_saving_release_change=at3('latest_U_package')-at3('wider_U'),
                    latest_sponsor_cash_NPV=at3('latest_U_package'))
                reconciliation.append(reconstruction)
                if horizon not in latest_rows:
                    economic.append(dict(strategy=name,market=market,horizon=horizon,status='unavailable_after_baseline_failure'))
                    liquidity.append(dict(strategy=name,market=market,horizon=horizon,status='unavailable_after_baseline_failure'))
                    continue
                components=economic_components(latest,q,b,profile,name,horizon,b['bootstrap_years'],old)
                economic.append(dict(strategy=name,market=market,horizon=horizon,status='conditional_price_proxy',
                    rates={str(rate):economic_metrics(components,rate,range(b['bootstrap_years'],horizon+1),cfg) for rate in cfg['real_discount_rates']}))
                row=latest_rows[horizon];op=row['operation']
                snap=p.snapshot(dict(row,gross=op['gross_USD2041'],available=op['available_USD2041'],core_imports=op['core_imports_USD2041']),profile,base['monetary'],latest['external_seed'])
                relief=max(0.,row['U_interest']-row['supplier_interest'])*cfg['local_U_rate_reduction']/.03
                snap['debt_service']-=relief;snap['baseline_annual_payout']+=relief
                for months in pcfg['liquidity']['outage_months']:
                    for option in ('cash_only','U_term_50','U_term_50_half_FX_access','external_bridge'):
                        liquidity.append(dict(strategy=name,market=market,horizon=horizon,months=months,policy=option,
                            **p.bridge(snap,pcfg['liquidity'],pcfg['policies'][option],months)))
            # Annual products make the owner-saving rule and cash timing inspectable.
            f=cashflows(latest,max(cfg['horizons']),b['bootstrap_years'],cfg['local_U_rate_reduction'])
            dated.append(dict(strategy=name,market=market,initial_sponsor_funding=[v+(latest['initial_reserve'] if year==0 else 0.) for year,v in enumerate(f['initial'][:b['bootstrap_years']])],
                annual=[dict(year=r['year'],scale=r['scale'],service_kits=r['service_kits'],
                             contractual_payout=f['dividends'][r['year']],conditional_saving_release=f['relief'][r['year']],
                             U_bank_gross_margin=f['bank_margin'][r['year']],
                             gross=r['operation']['gross_USD2041'],operating_cost=r['operation']['total_cost_USD2041'],
                             capital_orders=r['industrial_order']+r['service_order']) for r in latest['annual'] if 'failure' not in r]))
    paths=['consolidated_inputs.json','policy.py','results/policy.json','returns_inputs.json',*policy['dependencies']]
    cap=c.s.staged.capital_account(e,t)
    return dict(schema='terluna.research.consolidated-returns/1',evidence=cfg['evidence'],boundary=policy['boundary'],
        producer_sha256=c.s.run.digest(HERE/'consolidated.py'),dependencies={key:c.s.run.digest(HERE/key) for key in paths},
        architecture=cap,reading_rule='Sponsor return is not global welfare. Interest relief is a transfer with an explicit external-saving release rule. Terminal book claims and dated liquidity episodes are separate. Economic purchase-price proxies omit unmeasured surplus and externalities; failed paths receive no future extrapolation.',
        financial=financial,economic=economic,reconciliation=reconciliation,liquidity=liquidity,latest_annual=dated,
        sustained_cap_operation=historical['cap_operation_checks'])


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=c.s.run.ROOT/'research/runs/complex_returns/consolidated.json')
    args=parser.parse_args();args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(product(),indent=2,allow_nan=False)+'\n');print(args.output)
