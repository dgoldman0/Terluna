"""Conditional annual GDP and GNI distributions for the economic zone served by U.

This is an accounting consumer. It neither adds customer demand nor changes the
physical schedule, currency price, ownership contracts or financing histories.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
from .extract import HERE, canonical, digest, read, write


def triangular_ppf(u, low, mode, high):
    u = np.asarray(u)
    if high == low:
        return np.full_like(u, low, dtype=float)
    split = (mode - low) / (high - low)
    return np.where(u < split, low + np.sqrt(u * (high - low) * (mode - low)),
                    high - np.sqrt((1 - u) * (high - low) * (high - mode)))


def draws(spec, n, seed):
    """Independent stratified marginal draws, reused to compare the same cases."""
    rng = np.random.default_rng(seed)
    return {key: triangular_ppf((rng.permutation(n) + rng.random(n)) / n, *bounds)
            for key, bounds in spec['triangular_priors'].items()}


def accounts(row, base, p, spec):
    f = p['sales_realization']
    price = p['current_import_price']
    freight = p['capital_freight_share'] * row['earth_launch']
    base_freight = p['capital_freight_share'] * base['earth_launch']
    routine = p['routine_repair_share'] * row['repair_parts']
    base_routine = p['routine_repair_share'] * base['repair_parts']
    # Fixed physical support remains in place even when fewer exports sell.
    # Only the stated variable part of purchased commercial services falls.
    commercial = row['commercial_imports'] * (1 - spec['commercial_variable_cost_share'] * (1 - f))
    current_imports = price * (row['physical_current_imports'] - freight + commercial) + routine
    exports = f * row['exports']
    local_capital = p['local_capital_value'] * row['local_capital_output']
    # Existing asset quotes omit delivery and the unallocated service package.
    # At the tapered endpoint, annual capitalized freight is the additional
    # depreciation proxy. It is varied with the same depreciation multiplier.
    depreciation = p['depreciation_multiplier'] * (row['depreciation'] + price*freight + row['service_renewal'])
    base_d = p['depreciation_multiplier'] * (base['depreciation'] + price*base_freight)
    housing_d = p['depreciation_multiplier'] * base['housing_depreciation']
    net_ratio = base['net_assets'] / base['depreciation']
    q = p['resident_support_share']
    resident_depreciation = housing_d + q * (base_d - housing_d)
    resident_net_assets = base['housing_net_assets'] + q * (base['net_assets'] - base['housing_net_assets'] + price*base_freight * net_ratio)
    resident_current_inputs = q * (price * (base['physical_current_imports'] - base_freight) + base_routine + base['core_payroll'])
    resident_services = resident_current_inputs + resident_depreciation + p['cost_valuation_return'] * resident_net_assets
    # Housing, food and shared resident services are represented as household
    # benefits in kind. Equal output and compensation entries do not add cash
    # receipts or profit. Work-only support stays intermediate consumption.
    cash_wages = row['resident_payroll']
    compensation = cash_wages + resident_services
    gdp = exports + resident_services + local_capital - current_imports
    imported_capital = (row['imported_growth_capital'] + row['imported_renewal'] - routine
                        + price*freight + row['service_renewal'] + .9 * row['service_capital'])
    investment = local_capital + imported_capital
    imports = current_imports + imported_capital
    expenditure_gdp = resident_services + investment + exports - imports
    net_profit = gdp - compensation - depreciation - row['local_interest'] - row['external_interest']
    operating_cash = (row['annual_operating_cash'] + (f - 1)*row['exports']
                      - (price*(row['physical_current_imports'] + commercial)
                         - row['physical_current_imports'] - row['commercial_imports']))
    income_gdp = compensation + depreciation + net_profit + row['local_interest'] + row['external_interest']
    return dict(gdp=gdp, expenditure_gdp=expenditure_gdp, income_gdp=income_gdp,
        exports=exports, intermediate_imports=current_imports, resident_services=resident_services,
        cash_compensation=np.zeros_like(gdp) + cash_wages, compensation=compensation,
        own_account_capital=local_capital, imported_capital=imported_capital, gross_capital_formation=investment,
        total_imports=imports, depreciation=depreciation, net_profit=net_profit,
        external_interest=np.zeros_like(gdp) + row['external_interest'],
        local_interest=np.zeros_like(gdp) + row['local_interest'], capitalized_freight=price*freight,
        routine_repair_inputs=routine,
        unallocated_resource_provision=np.zeros_like(gdp) + row['resource_provision'],
        year200_operating_cash=operating_cash)


def income(a, local_share, outside_income=0., foreign_direct_investment_share=1.):
    """Accrued FDI income includes retained earnings and signed losses.

The primary study treats outside equity as direct investment. A portfolio-only
ownership structure needs its own dividends and income attribution account.
"""
    if foreign_direct_investment_share != 1.:
        raise ValueError('Portfolio-equity earnings require a separate income account')
    foreign_equity_income = (1 - local_share) * a['net_profit']
    gni = a['gdp'] - a['external_interest'] - foreign_equity_income + outside_income
    return dict(gni=gni, gnp=gni, foreign_equity_income=foreign_equity_income,
                net_external_earned_income=gni - a['gdp'])


def weighted_quantile(x, w, probabilities):
    x, w = np.asarray(x), np.asarray(w)
    keep = w > 0
    x, w = x[keep], w[keep]
    if not len(x) or not np.all(np.isfinite(x)) or not np.all(np.isfinite(w)):
        raise ValueError('Empty or nonfinite distribution')
    order = np.argsort(x, kind='stable')
    x, w = x[order], w[order]
    cumulative = np.cumsum(w) / np.sum(w)
    # Inverse empirical CDF; no interpolation across unsupported scenario gaps.
    return x[np.minimum(np.searchsorted(cumulative, probabilities, side='left'), len(x) - 1)]


def statistics(values, weights, curve=False):
    values, weights = np.asarray(values), np.asarray(weights)
    weights = weights / weights.sum()
    ps = [.01, .05, .1, .25, .5, .75, .9, .95, .99]
    qs = weighted_quantile(values, weights, ps)
    mean = float(np.sum(weights * values))
    out = dict(mean=mean, standard_deviation=float(np.sqrt(np.sum(weights * (values - mean)**2))),
               quantiles={str(p): float(v) for p, v in zip(ps, qs)},
               probability_below={str(t): float(weights[values < t].sum()) for t in (0, 1e12, 5e12, 10e12, 20e12, 50e12)})
    if curve:
        ps = np.linspace(0, 1, 501)
        out['inverse_cdf'] = dict(probability=ps.tolist(), value=weighted_quantile(values, weights, ps).tolist())
    return out


def correlation(x, y, w):
    w = w / w.sum()
    a, b = x - np.sum(w*x), y - np.sum(w*y)
    d = np.sqrt(np.sum(w*a*a) * np.sum(w*b*b))
    return float(np.sum(w*a*b) / d) if d > 0 else None


def scenario_weight(row, prior, spec):
    return prior['market'][spec['market_order'].index(row['market'])] * prior['growth'][spec['growth_order'].index(row['growth'])]


def aggregate(selected, account_rows, weights, spec, actual_ownership=None, curve=True):
    weights = np.asarray(weights, dtype=float)
    weights /= weights.sum()
    n = np.asarray(account_rows[0]['gdp']).size
    w = np.repeat(weights / n, n)
    gdp = np.concatenate([a['gdp'] for a in account_rows])
    gni = {}
    ownership = spec['ownership_shares'] if actual_ownership is None else actual_ownership
    for label, share in ownership.items():
        values = np.concatenate([income(a, row['local_shares'][share] if isinstance(share, str) else share)['gni']
                                 for row, a in zip(selected, account_rows)])
        gni[label] = statistics(values, w, curve=curve)
        gni[label]['correlation_with_gdp'] = correlation(gdp, values, w)
    means = np.array([np.mean(a['gdp']) for a in account_rows])
    between = float(np.sum(weights * (means - np.sum(weights * means))**2))
    within = float(np.sum(weights * [np.var(a['gdp']) for a in account_rows]))
    return dict(cases=[dict(name=r['name'], probability=float(p)) for r, p in zip(selected, weights)],
                gdp=statistics(gdp, w, curve=curve), gni=gni,
                probability_year200_operating_deficit=float(w[np.concatenate([a['year200_operating_cash'] for a in account_rows]) < 0].sum()),
                gdp_variance=dict(between_market_cases=between, within_case=within,
                                 fraction_between_cases=between/(between+within) if between+within else 0.))


def verify_inputs(spec, data):
    if digest(HERE/spec['extracted_input']['path']) != spec['extracted_input']['sha256']:
        raise ValueError('Changed extracted input; reproduce and admit its new identity')
    if data['schema'] != spec['extracted_input']['schema']:
        raise ValueError('Wrong extracted input schema')
    if data['producer']['sha256'] != digest(HERE/'extract.py'):
        raise ValueError('Extracted input uses a different extractor')
    for item in spec['products'].values():
        if digest(HERE/item['path']) != item['sha256'] or data['dependencies'][item['path']] != item['sha256']:
            raise ValueError('Changed source: ' + item['path'])
    for path, value in spec['implementation_dependencies'].items():
        if digest(HERE/path) != value or data['implementation_dependencies'][path] != value:
            raise ValueError('Changed implementation: ' + path)
    names = {r['name'] for r in data['cases']}
    expected = {f'{p}:{m}:{g}' for p in spec['portfolios'] for m in spec['market_order'] for g in spec['growth_order']}
    if names != expected or len(data['cases']) != len(expected):
        raise ValueError('Expected all nine portfolio/market/growth combinations')
    for r in data['cases']:
        if not np.isclose(sum(r['revenue_by_sector'].values()), r['exports'], rtol=1e-12, atol=1):
            raise ValueError('Sector receipts do not reconcile')
        if not np.isclose(r['population'], data['cases'][0]['population']):
            raise ValueError('Different completion populations')
    if any(any(p < 0 or p > 1 for p in v[k]) or abs(sum(v[k])-1)>1e-12 for v in spec['scenario_priors'].values() for k in ('market','growth')):
        raise ValueError('Invalid prior probabilities')
    for key, bounds in spec['triangular_priors'].items():
        if len(bounds) != 3 or not np.all(np.isfinite(bounds)) or not bounds[0] <= bounds[1] <= bounds[2]:
            raise ValueError('Invalid triangular prior: ' + key)


def product(spec, data, n=None, seed=None):
    n = n or spec['samples_per_case']
    seed = spec['seed'] if seed is None else seed
    p = draws(spec, n, seed)
    modes = {k: np.array([v[1]]) for k,v in spec['triangular_priors'].items()}
    rows = data['cases']
    base = next(r for r in rows if r['name'] == spec['base_resident_case'])
    out, byname, centres = {}, {}, []
    controls = dict(max_production_expenditure_gap=0., max_production_income_gap=0., max_gni_identity_gap=0.)
    for r in rows:
        a = accounts(r, base, p, spec)
        byname[r['name']] = a
        c = accounts(r, base, modes, spec)
        centre = {k:float(v[0]) for k,v in c.items()}
        centre.update(name=r['name'], asset_method=r['asset_method'], resource_treatment_sensitivity={
            k:dict(base_gni=float(income(c,s)['gni'][0]),
                   all_provision_current_local_levy_gni=float(income(c,s)['gni'][0]+(1-s)*r['resource_provision']),
                   all_provision_external_current_service_gni=float(income(c,s)['gni'][0]-s*r['resource_provision']))
            for k,s in spec['ownership_shares'].items()}, gni={k:float(income(c,s)['gni'][0]) for k,s in spec['ownership_shares'].items()})
        centres.append(centre)
        controls['max_production_expenditure_gap'] = max(controls['max_production_expenditure_gap'], float(np.max(np.abs(a['gdp'] - a['expenditure_gdp']))))
        controls['max_production_income_gap'] = max(controls['max_production_income_gap'], float(np.max(np.abs(a['gdp'] - a['income_gdp']))))
        for s in spec['ownership_shares'].values():
            inc = income(a, s)
            independently = a['compensation'] + a['depreciation'] + a['local_interest'] + s*a['net_profit']
            controls['max_gni_identity_gap'] = max(controls['max_gni_identity_gap'], float(np.max(np.abs(inc['gni']-independently))))
    if max(controls.values()) > 1.:
        raise ValueError('Annual account does not close within one USD2041')
    for prior_name, prior in spec['scenario_priors'].items():
        out[prior_name] = {}
        for portfolio in spec['portfolios']:
            rs = [r for r in rows if r['portfolio'] == portfolio]
            ps = np.array([scenario_weight(r, prior, spec) for r in rs])
            screens = dict(funded_completion=rs,
                previously_self_financing_paths=[r for r in rs if r['self_financing'] and r['annual_operating_cash'] > 0 and r['annual_post_capital_cash'] > 0])
            result = {}
            for screen, selected in screens.items():
                if not selected:
                    result[screen] = None
                    continue
                weights = [scenario_weight(r, prior, spec) for r in selected]
                result[screen] = aggregate(selected, [byname[r['name']] for r in selected], weights, spec)
                result[screen]['retained_scenario_weight'] = float(sum(weights)/ps.sum())
            out[prior_name][portfolio] = result
    # Keep the genuinely computed ownership histories separate from hypothetical
    # endpoint shares. No extrapolation of vesting or buyouts into 63 other cases.
    actual = {}
    for portfolio in spec['portfolios']:
        selected = [r for r in rows if r['portfolio']==portfolio and r['local_shares']]
        weights = [scenario_weight(r, spec['scenario_priors']['balanced'], spec) for r in selected]
        actual[portfolio] = aggregate(selected, [byname[r['name']] for r in selected], weights, spec,
            actual_ownership={k:k for k in ('external_equity','community_purchase','negotiated_vesting')})
    return dict(schema='terluna.humanities.production-income/1', year=200,
        units='Annual constant USD2041 purchasing-power equivalents; no U exchange rate or nominal-currency forecast.',
        boundary='Resident producers and institutions throughout the planned off-Earth economic zone served by U, before Open Moon, the habitat array and outer-system development.',
        reading_rules=[
            'Treat each portfolio and ownership arrangement as a separate decision; do not pool their probabilities.',
            'All distributions assume the planned year-200 population and physical build. The previously self-financing filter describes source histories, not financing under every new draw.',
            'Scenario weights and triangular bounds are analyst assumptions; these intervals do not include every structural uncertainty or a probability of completion.',
            'GNP equals GNI in this reconciled model. Both include depreciation; neither is cash income, investor payout, wealth or currency supply.',
            'Outside equity is direct investment: accrued earnings include retained profits and signed losses. Portfolio-equity treatment requires a separate distribution policy.',
            'The actual-ownership comparison uses only the 18 matched source cases. Fixed local shares elsewhere are counterfactual endpoints.',
            'Figures, means and quantiles refer to the supplied production boundary and probability prior; they are not independently calibrated economic forecasts.'],
        evidence=spec['evidence'], producer=dict(path='research/studies/production_income/model.py', sha256=digest(__file__)),
        dependencies={p:digest(HERE/p) for p in ('income_inputs.json','extract.py','inputs/year200.json.gz','income_sources.json')},
        sampling=dict(method='Independent randomized stratified marginal sampling; common draws across scenarios and ownership cases.',
                      samples_per_case=n, seed=seed, scenarios=len(rows), draws_evaluated=n*len(rows)),
        calibration=data['calibration'], controls=controls, priors=spec['scenario_priors'], triangular_priors=spec['triangular_priors'],
        distributions=out, actual_ownership_conditional_breakthrough=actual, mode_accounts=centres)


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--samples',type=int)
    parser.add_argument('--seed',type=int)
    parser.add_argument('--output',default=str(HERE/'results/year200.json.gz'))
    args=parser.parse_args()
    spec=read(HERE/'income_inputs.json')
    data=read(HERE/'inputs/year200.json.gz')
    verify_inputs(spec,data)
    result=product(spec,data,args.samples,args.seed)
    write(args.output,result)
    print('Completed',result['sampling'],'; accounting controls:',result['controls'])
