"""Extract a year-200 account from admitted products; never resize industry."""
from __future__ import annotations
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n').encode()


def read(path):
    raw = Path(path).read_bytes()
    return json.loads(gzip.decompress(raw) if str(path).endswith('.gz') else raw)


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = canonical(value)
    path.write_bytes(gzip.compress(raw, mtime=0) if str(path).endswith('.gz') else raw)


def admitted(spec, name):
    item = spec['products'][name]
    path = HERE / item['path']
    if digest(path) != item['sha256']:
        raise ValueError('Source hash mismatch: ' + name)
    data = read(path)
    if data['schema'] != item['schema']:
        raise ValueError('Source schema mismatch: ' + name)
    return data


def year(rows, y=200):
    matches = [r for r in rows if r['year'] == y]
    if len(matches) != 1:
        raise ValueError('Expected one annual row at year ' + str(y))
    return matches[0]


def extract(spec):
    # Finance belongs to this research lane. Physical quantities arrive only in
    # pinned provisioning products, not through a provisioning implementation.
    from ..enclosed_returns import finance
    for path, expected in spec['implementation_dependencies'].items():
        if digest(HERE / path) != expected:
            raise ValueError('Changed finance dependency: ' + path)
    orders = admitted(spec, 'orders')
    assets = admitted(spec, 'assets')
    wealth = admitted(spec, 'wealth')
    returns = admitted(spec, 'returns')
    physical = admitted(spec, 'physical')
    flows = admitted(spec, 'flows')
    if assets['source_physical_sha256'] != spec['products']['physical']['sha256']:
        raise ValueError('Asset and physical products do not match')
    cfg = read(HERE / spec['finance_config'])
    pmap = {c['baseline_case']: c for c in physical['cases'] if c['variant'] == 'combined'}
    amap = {c['baseline_case']: c for c in assets['cases']}
    fmap = {c['baseline_case']: c for c in flows['cases']}
    rmap = {c['name']: c for c in returns['primary_cases'] if c['variant'] == 'combined'}
    wmap = {}
    for c in wealth['cases']:
        wmap.setdefault(c['baseline_case'], {})[c['ownership']] = year(c['annual'])['local_economic_ownership_fraction']
    if set(pmap) != {c['name'] for c in orders['cases']}:
        raise ValueError('Incomplete physical/order coverage')
    ratios, net_ratios, repair_ratios, asset_rows = [], [], [], {}
    for name, c in amap.items():
        pc = pmap[name]
        if hashlib.sha256(canonical(pc['physical_case'])).hexdigest() != c['source_case_sha256']:
            raise ValueError('Asset case mismatch: ' + name)
        if hashlib.sha256(canonical(pc)).hexdigest() != c['source_loop_case_sha256']:
            raise ValueError('Asset loop mismatch: ' + name)
        row = year(c['annual'])
        cash = year(pc['physical_case']['annual'])['cash']
        denominator = cash['operating_external_replacement_USD2041'] + cash['local_replacement_account_value_USD2041']
        ratios.append(row['totals']['consumption_of_fixed_capital_proxy_USD2041'] / denominator)
        net_ratios.append(row['totals']['net_remaining_service_value_USD2041'] / row['totals']['consumption_of_fixed_capital_proxy_USD2041'])
        repair_ratios.append(row['totals']['repair_parts_purchase_USD2041'] / cash['operating_external_replacement_USD2041'])
        asset_rows[name] = row
    # A bounded extrapolation for 63 orderbooks without a full asset export.
    # Its error against the 18 exported cases is reported, not hidden.
    ratio = float(np.median(ratios))
    loo = [float(np.median(np.delete(ratios, i))) / v - 1 for i, v in enumerate(ratios)]
    housing = [r['by_functional_sector']['settlement'] for r in asset_rows.values()]
    for key in ('consumption_of_fixed_capital_proxy_USD2041', 'net_remaining_service_value_USD2041'):
        if max(r[key] for r in housing) - min(r[key] for r in housing) > 1:
            raise ValueError('Housing stock cannot be shared across the fixed-population cases')
    result = []
    for c in orders['cases']:
        name = c['name']
        pc = pmap[name]
        if hashlib.sha256(canonical(c)).hexdigest() != pc['source_order_case_sha256']:
            raise ValueError('Orders and physical case mismatch: ' + name)
        pr = year(pc['physical_case']['annual'])
        request = year(c['annual'])
        path = finance.plan(c, pc['physical_case'], cfg)
        f = year(path['annual'])
        reference = next(r for r in rmap[name]['checkpoints'] if r['horizon'] == 200 and r['mode'] == 'U_credit')
        if f['gross'] != reference['annual_gross'] or f['operating_available'] != reference['annual_operating_available']:
            raise ValueError('Changed annual cash account: ' + name)
        if name in fmap:
            s = year(fmap[name]['annual'])
        else:
            eq = finance.finance(path, cfg, 'equity')
            sim = finance.finance(path, cfg, 'U_credit', eq)
            s = year(sim['annual'])
            if abs(s['debt'] - reference['terminal_debt']) > 1:
                raise ValueError('Changed financing account: ' + name)
        cash = pr['cash']
        renewal = cash['operating_external_replacement_USD2041'] + cash['local_replacement_account_value_USD2041']
        a = asset_rows.get(name)
        d = a['totals']['consumption_of_fixed_capital_proxy_USD2041'] if a else renewal * ratio
        n = a['totals']['net_remaining_service_value_USD2041'] if a else d * float(np.median(net_ratios))
        repairs = a['totals']['repair_parts_purchase_USD2041'] if a else cash['operating_external_replacement_USD2041'] * float(np.median(repair_ratios))
        result.append(dict(name=name, portfolio=c['portfolio'], market=c['market'], growth=c['market_growth_name'],
            population=pr['population'], exports=request['gross_USD2041'], revenue_by_sector=request['revenue_by_sector_USD2041'],
            physical_current_imports=cash['operating_external_logistics_materials_USD2041'] + cash['operating_external_services_USD2041'],
            earth_launch=min(cash['Earth_launch_USD2041'], cash['operating_external_logistics_materials_USD2041']),
            imported_renewal=cash['operating_external_replacement_USD2041'], repair_parts=repairs,
            imported_growth_capital=cash['capital_external_USD2041'],
            local_capital_output=cash['local_capital_account_value_USD2041'] + cash['local_replacement_account_value_USD2041'],
            core_payroll=cash['operating_local_payroll_USD2041'], resident_payroll=f['resident_payroll'],
            commercial_imports=request['commercial_operating_USD2041'] - request['commercial_payroll_USD2041'],
            service_renewal=f['service_residual_renewal'], service_capital=f['service_residual_capital'],
            resource_provision=f['resource_provision'], unallocated_capital_allowance=f['capital_cost_allowance'],
            depreciation=d, net_assets=n, housing_depreciation=housing[0]['consumption_of_fixed_capital_proxy_USD2041'],
            housing_net_assets=housing[0]['net_remaining_service_value_USD2041'], asset_method='exported' if a else 'renewal_proxy',
            external_interest=s['external_interest'] + s['supplier_interest'], local_interest=s['local_U_interest'],
            annual_operating_cash=f['operating_available'], annual_owner_payout=s['owner_payout'],
            self_financing=reference['no_additional_equity_to_checkpoint'],
            annual_post_capital_cash=reference['annual_post_capital_cash'],
            additional_equity=reference['additional_equity_required'], local_shares=wmap.get(name, {}),
            source_case_sha256=hashlib.sha256(canonical(pc['physical_case'])).hexdigest()))
        print('Year-200 account:', name, flush=True)
    return dict(schema='terluna.humanities.production-income-inputs/1',
        producer=dict(path='research/studies/production_income/extract.py', sha256=digest(__file__)),
        dependencies={item['path']: item['sha256'] for item in spec['products'].values()},
        implementation_dependencies=spec['implementation_dependencies'],
        calibration=dict(exported_cases=len(ratios), renewal_to_depreciation_median=ratio,
            observed_ratio_range=[min(ratios), max(ratios)],
            leave_one_out_relative_error_range=[min(loo), max(loo)],
            median_net_stock_to_depreciation=float(np.median(net_ratios)),
            median_repair_share=float(np.median(repair_ratios)),
            evidence='Internal approximation to engineering asset proxies; not empirical depreciation validation.'),
        cases=result)


if __name__ == '__main__':
    write(HERE / 'inputs/year200.json.gz', extract(read(HERE / 'income_inputs.json')))
