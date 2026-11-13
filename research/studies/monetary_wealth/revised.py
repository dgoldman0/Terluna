"""Revalue the revised physical/financial paths with unchanged ownership rules.

The consumer reads frozen products. It never imports provisioning code or adds
mineral potential, an equity claim and its underlying assets as separate wealth.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import math
from pathlib import Path

from . import model

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.humanities.revised-monetary-wealth/1'
STOCK_FIELDS = ('cash', 'reserves', 'owner_custody', 'debt', 'local_U_debt',
                'supplier_U_debt', 'external_debt')
FLOW_FIELDS = ('gross', 'resident_payroll', 'owner_payout', 'dividend',
               'outside_equity_contribution', 'local_U_draw', 'supplier_U_draw',
               'external_draw', 'local_U_principal', 'supplier_U_principal',
               'external_principal', 'local_U_interest', 'supplier_interest',
               'external_interest')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_product(spec, base=HERE):
    path = base / spec['path']
    if not spec.get('sha256'):
        raise ValueError('Unadmitted product: ' + spec['path'])
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != spec['sha256']:
        raise ValueError('Changed product: ' + spec['path'])
    result = json.loads(gzip.decompress(raw) if path.suffix == '.gz' else raw)
    if result.get('schema') != spec['schema']:
        raise ValueError('Unexpected product schema: ' + spec['path'])
    return result


def _close(a, b):
    return abs(a-b) <= max(1., 1e-11*max(abs(a), abs(b)))


def validate_finance(case, settings):
    rows = case['annual']
    if [r['year'] for r in rows] != list(range(settings['years'])):
        raise ValueError('Missing or repeated finance date: ' + case['name'])
    old = {kind: 0. for kind in ('local_U', 'supplier_U', 'external')}
    old_custody = 0.
    for row in rows:
        for key in STOCK_FIELDS + FLOW_FIELDS:
            value = row[key]
            if not isinstance(value, (int, float)) or not math.isfinite(value) or value < -1.:
                raise ValueError('Invalid finance amount: ' + case['name'] + '/' + key)
        if not _close(row['debt'], sum(row[k+'_debt'] for k in old)):
            raise ValueError('Debt sectors do not equal programme debt')
        for kind in old:
            closing = old[kind] + row[kind+'_draw'] - row[kind+'_principal']
            if not _close(row[kind+'_debt'], closing):
                raise ValueError('Debt stock/flow mismatch: ' + kind)
            old[kind] = row[kind+'_debt']
        # Custody belongs to the owners already entitled to a declared payout.
        if not _close(old_custody + row['dividend'] - row['owner_payout'], row['owner_custody']):
            raise ValueError('Payout and custody stocks do not reconcile')
        old_custody = row['owner_custody']


def matched_cases(finance, assets, cfg, settings):
    physical = finance.get('source_physical_sha256')
    if not physical or assets.get('source_physical_sha256') != physical:
        raise ValueError('Finance and asset products have different physical sources')
    indexed = {}
    for case in assets['cases']:
        key = (case['baseline_case'], case['variant'])
        if key in indexed:
            raise ValueError('Duplicate physical asset case')
        indexed[key] = case
    pairs = []
    seen = set()
    for case in finance['cases']:
        if case['market'] != settings['market'] or case['growth'] not in settings['growth_cases']:
            raise ValueError('Financial case outside the declared market/growth scope')
        key = (case.get('baseline_case', case['name']), case['variant'])
        if key[1] != cfg['variant'] or key in seen:
            raise ValueError('Unexpected or duplicate financial variant')
        seen.add(key)
        asset = indexed.get(key)
        if asset is None:
            raise ValueError('No matching asset stock: ' + repr(key))
        for identity in ('source_case_sha256', 'source_loop_case_sha256'):
            if not case.get(identity) or case[identity] != asset.get(identity):
                raise ValueError('Financial and asset case identities differ: '+identity)
        if [r['year'] for r in asset['annual']] != list(range(settings['years'])):
            raise ValueError('Missing or repeated asset date')
        if case['name'] != key[0]:
            raise ValueError('Finance names must retain the baseline case identity')
        for row in asset['annual']:
            for field in ('net_remaining_service_value_USD2041', 'construction_WIP_value_USD2041'):
                if not math.isfinite(row['totals'][field]) or row['totals'][field] < 0:
                    raise ValueError('Invalid included asset value')
        validate_finance(case, settings)
        pairs.append((case, asset))
    if len(pairs) != cfg['expected_financial_cases']:
        raise ValueError('Incomplete financial case set')
    if len({c['portfolio'] for c, _ in pairs}) != cfg['expected_portfolios']:
        raise ValueError('Incomplete portfolio set')
    for portfolio in {c['portfolio'] for c, _ in pairs}:
        if {c['growth'] for c, _ in pairs if c['portfolio'] == portfolio} != set(settings['growth_cases']):
            raise ValueError('Incomplete growth pair: ' + portfolio)
    return pairs


def checkpoint_summary(row):
    """Dated cash and terminal claims stay separate, including custody rights."""
    book = row['book_balance']
    local = book['consolidated_off_earth']
    sectors = book['sectors']
    cumulative = row['cumulative']
    outside = row['outside_investor_cash_returns']
    share = row['local_economic_ownership_fraction']
    own = book['own_funds_equity']
    outside_custody = row['claims']
    outside_custody = math.fsum(c['value'] for c in outside_custody
                               if c['instrument'] == 'custody'
                               and c['holder'] == 'outside_owners_customers')
    terminal_outside = (1-share)*own + outside_custody
    marked = row['cash_claim_balance']
    result = {
        'gross_receipts_cumulative': cumulative['gross_receipts'],
        'outside_capital_contributions_cumulative': cumulative['outside_contributions'],
        'programme_payouts_cumulative': cumulative['outside_payouts']+cumulative['local_payouts'],
        'outside_payouts_cumulative': cumulative['outside_payouts'],
        'outside_buyout_receipts_cumulative': cumulative['outside_buyout_receipts'],
        'local_payouts_cumulative': cumulative['local_payouts'],
        'household_net_saving_cumulative': cumulative['household_saving'],
        'community_purchase_spending_cumulative': cumulative['purchase_spending'],
        'negotiated_transfer_book_value_cumulative': cumulative['gift_book_value'],
        'local_economic_ownership_fraction': share,
        'included_real_assets': local['real_assets'],
        'net_external_position': local['net_international_investment_position'],
        'consolidated_local_net_wealth': local['net_worth'],
        'household_net_wealth': sectors['households']['net_worth'],
        'community_fund_net_wealth': sectors['community_fund']['net_worth'],
        'community_uninvested_financial_assets': row['community_uninvested_financial_assets'],
        'household_retained_financial_assets': row['household_retained_financial_assets'],
        'domestic_issuer_transaction_money': book['money']['resident_transaction_money']+book['money']['external_transaction_money'],
        'resident_transaction_money': book['money']['resident_transaction_money'],
        'outside_transaction_money': book['money']['external_transaction_money'],
        'bank_settlement_reserves': book['money']['resident_bank_reserves'],
        'monetary_base': book['money']['monetary_base'],
        'programme_cash': row['programme_cash'],
        'programme_reserves': row['programme_reserves'],
        'programme_debt': row['programme_debt'],
        'programme_local_U_debt': row['programme_local_U_debt'],
        'programme_supplier_U_debt': row['programme_supplier_U_debt'],
        'programme_external_debt': row['programme_external_debt'],
        'enterprise_own_funds_equity': own,
        'local_enterprise_equity_at_own_funds': share*own,
        'outside_enterprise_equity_at_own_funds': (1-share)*own,
        'outside_custody_claim': outside_custody,
        'outside_cash_NPV_at_3percent': outside['dated_cash_NPV_at_3percent'],
        'outside_cash_IRR_if_unique': outside['nonnegative_IRR_if_unique'],
        'outside_cash_plus_terminal_claim_NPV_at_3percent': outside['dated_cash_NPV_at_3percent']+terminal_outside/1.03**row['year'],
        'conditional_cash_claim_equity_quote': row['equity_quote']['quote'],
        'cash_claim_marked_consolidated_local_net_wealth': None if marked is None else marked['consolidated_off_earth']['net_worth'],
        'cash_claim_marked_community_net_wealth': None if marked is None else marked['sectors']['community_fund']['net_worth'],
        'cash_claim_marked_enterprise_residual': None if marked is None else marked['sectors']['enterprises']['net_worth'],
        'remaining_modelled_capital_calls': row['equity_quote']['remaining_capital_calls'],
    }
    return result


def baseline_checkpoints(product):
    return {(case['name'], case['ownership']): {
        year: checkpoint_summary(row) for year, row in case['checkpoints'].items()
    } for case in product['cases']}


def compare_path(path, baseline):
    old = baseline[(path['name'], path['ownership'])]
    comparisons = []
    for year, row in path['checkpoints'].items():
        current = checkpoint_summary(row)
        previous = old[year]
        delta = {k: None if current[k] is None or previous[k] is None else current[k]-previous[k]
                 for k in current}
        comparisons.append({'name': path['name'], 'portfolio': path['portfolio'],
                            'growth': path['growth'], 'ownership': path['ownership'],
                            'year': int(year), 'baseline': previous, 'revised': current,
                            'change': delta})
    return comparisons


def compact_path(path, finance_case):
    """Keep annual book accounts; full alternate marks and claims at checkpoints."""
    result = dict(path)
    rows = []
    checkpoints = set(map(int, path['checkpoints']))
    for source, flow in zip(path['annual'], finance_case['annual']):
        row = dict(source)
        if not _close(row['flows']['local_payouts']+row['flows']['outside_payouts'], flow['owner_payout']):
            raise ValueError('Ownership changed programme payouts')
        for exported, original in [('programme_cash', 'cash'), ('programme_reserves', 'reserves'), ('programme_debt', 'debt')]:
            if row[exported] != flow[original]:
                raise ValueError('Ownership changed programme balance')
        if row['flows']['gross_receipts'] != flow['gross'] or row['flows']['outside_contributions'] != flow['outside_equity_contribution']:
            raise ValueError('Ownership changed operating receipts or capital calls')
        row['cash_claim_equity_mark'] = row['equity_quote']['quote']
        if row['year'] not in checkpoints:
            for field in ('cash_claim_balance', 'physical_by_function', 'cumulative'):
                row.pop(field, None)
        rows.append(row)
    result['annual'] = rows
    result['checkpoints'] = {str(y): rows[y] for y in sorted(checkpoints)}
    result['variant'] = finance_case['variant']
    result['baseline_case'] = finance_case.get('baseline_case', finance_case['name'])
    result['physical_supply_qualified'] = False
    result['reading_rule'] = ('Annual book accounts retain every sector. Full cash-claim valuation accounts, cumulative totals and claim lists are retained at checkpoints. Counterparty closure does not qualify mining access or commissioning. Terminal equity includes enterprise cash and restricted reserves; those assets must not be added again to the equity claim.')
    return result


def load(config_path=HERE/'revised_inputs.json'):
    cfg = json.loads(Path(config_path).read_text())
    for path, expected in cfg['implementation_dependencies'].items():
        if sha(HERE/path) != expected:
            raise ValueError('Changed accounting dependency: '+path)
    settings = read_product(cfg['settings_product'])
    finance = read_product(cfg['finance_product'])
    if sha(ROOT/finance['producer']['path']) != finance['producer']['sha256']:
        raise ValueError('Changed revised finance producer')
    finance_base = HERE/cfg['finance_product']['dependency_base']
    for name, expected in finance['dependencies'].items():
        if sha(finance_base/name) != expected:
            raise ValueError('Changed revised finance dependency: '+name)
    if finance['source_physical_sha256'] not in finance['dependencies'].values():
        raise ValueError('Financial physical source is not pinned among its dependencies')
    policy = json.loads((HERE/cfg['finance_policy_path']).read_text())
    if settings['owner_custody_fraction'] != policy['monetary']['incremental_distribution_spending_fraction']:
        raise ValueError('Wealth custody allocation differs from the source finance policy')
    assets = read_product(cfg['asset_product'])
    # Reduce the 54-path historical product before retaining new annual results.
    old = read_product(cfg['baseline_product'])
    baseline = baseline_checkpoints(old)
    baseline_metadata = {k: old[k] for k in ('schema', 'producer', 'unvalued_assets')}
    del old
    return cfg, settings, finance, assets, baseline, baseline_metadata


def build(cfg, settings, finance, assets, baseline, baseline_metadata):
    pairs = matched_cases(finance, assets, cfg, settings)
    results = []
    comparisons = []
    for case, asset in pairs:
        for ownership in settings['ownership_scenarios']:
            path = model.path(case, asset, settings, ownership)
            comparisons.extend(compare_path(path, baseline))
            results.append(compact_path(path, case))
    expected = cfg['expected_financial_cases']*len(settings['ownership_scenarios'])
    if len(results) != expected:
        raise ValueError('Incomplete ownership comparison')
    return {
        'schema': SCHEMA,
        'producer': {'path': 'research/studies/monetary_wealth/revised.py', 'sha256': sha(__file__)},
        'source_physical_sha256': finance['source_physical_sha256'],
        'source_finance_schema': finance['schema'],
        'source_asset_schema': assets['schema'],
        'settings': settings,
        'boundary': 'Whole planned off-Earth enclosed zone, with Earth as an explicit external counterparty. No Open Moon, array or outer-system expansion.',
        'units': 'constant USD2041 reference purchasing-power values under the U accounting convention; no nominal U face or exchange rate is inferred',
        'evidence': 'Conditional balance sheets and ownership scenarios for revised physical and financial products. The accounting policies match the retained baseline. Reclamation and restored service life remain engineering scenarios; money holding, community ownership and valuation remain institutional assumptions.',
        'cases': results,
        'baseline_comparisons': comparisons,
        'historical_product': baseline_metadata,
        'unvalued_assets': list(dict.fromkeys(assets.get('unvalued_assets', [])+baseline_metadata['unvalued_assets'])),
        'rules': cfg['rules'],
        'controls': {'financial_cases': len(pairs), 'ownership_paths': len(results),
                     'years_per_path': settings['years'],
                     'identical_physical_source_for_cash_and_assets': True,
                     'programme_receipts_payouts_and_balances_preserved': True,
                     'baseline_financial_and_wealth_products_modified': False,
                     'terminal_equity_added_to_distributed_cash_returns': False},
    }


def product(config_path=HERE/'revised_inputs.json'):
    loaded = load(config_path)
    result = build(*loaded)
    cfg = loaded[0]
    dependencies = {'revised_inputs.json': sha(config_path)}
    dependencies.update({spec['path']: spec['sha256'] for key, spec in cfg.items()
                         if key.endswith('_product')})
    dependencies.update(cfg['implementation_dependencies'])
    result['dependencies'] = dependencies
    return result


def main():
    result = product()
    path = HERE/'results/revised.json.gz'
    raw = (json.dumps(result, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode()
    path.write_bytes(gzip.compress(raw, mtime=0))
    print(path, len(result['cases']), sha(path))


if __name__ == '__main__':
    main()
