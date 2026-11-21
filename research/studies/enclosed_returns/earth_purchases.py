"""Decompose the unchanged commercial orders into purchased mass and invoices.

All 81 primary orderbooks are alternatives. No additional coproduct sales,
refinery ramp, price forecast or financing revenue is introduced here.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    path = Path(path)
    return json.loads(gzip.decompress(path.read_bytes()) if path.suffix == '.gz' else path.read_bytes())


def inputs():
    manifest = read(HERE / 'earth_purchases_inputs.json')
    for path, expected in manifest['dependencies'].items():
        if digest(HERE / path) != expected:
            raise ValueError('Changed purchase-account dependency: ' + path)
    return (manifest, read(HERE / 'results/orders.json'),
            read(HERE / 'enclosed_inputs.json'),
            read(HERE / '../complex_returns/strategy_inputs.json'),
            read(HERE / '../complex_returns/inputs/strategy_profiles.json'),
            read(HERE / 'results/revised.json.gz'))


def close(a, b):
    if not math.isclose(a, b, rel_tol=2e-11, abs_tol=0.01):
        raise ValueError(f'Purchase-account mismatch: {a} != {b}')


def annual(row, strategy, cfg, q, inflation):
    """Recover accepted premium mass from the two-price materials invoice.

    The orders product preserves production mass and revenues but omits the
    accepted premium mass. Its invoice uniquely determines that omitted mass.
    A premium attempt that fails acceptance recovers only the fallback share.
    """
    delivered_yield = cfg['delivered_quality_yield'] * (1 - cfg['return_loss_probability'])
    goods = row['goods'] * delivered_yield
    source = row['source_goods'] * delivered_yield
    share = strategy['premium_share']
    fallback = q['premium']['fallback_recovery_fraction']
    ordinary_price = strategy['lunar_price_USD2026_kg'] * inflation
    premium_price = (q['premium']['equivalent_product_offer_USD2026_kg']
                     * q['premium']['willingness_to_pay_multiple']
                     * q['competition']['premium_competitor_price_factor'] * inflation)
    revenues = dict(row['revenue_by_sector_USD2041'])
    fallback_mass = goods * (1 - share + share * fallback)
    premium = 0.0
    if share and goods:
        premium = ((revenues['materials'] - fallback_mass * ordinary_price)
                   / (premium_price - fallback * ordinary_price))
        if premium < -0.01 or premium > goods * share + 0.01:
            raise ValueError('Accepted premium mass outside attempted mass')
        premium = min(goods * share, max(0.0, premium))
    standard = fallback_mass - fallback * premium
    sold = standard + premium + source
    close(standard * ordinary_price + premium * premium_price, revenues['materials'])
    close(source * strategy['source_price_USD2026_kg'] * inflation, revenues['source_products'])
    close(sum(revenues.values()), row['gross_USD2041'])
    mass = dict(lunar_standard_products=standard,
                lunar_qualified_premium_products=premium,
                asteroid_refined_products=source,
                Earth_purchased_goods=sold,
                customer_propellant_delivered_in_space=row['propellant_sales'],
                gross_export_goods_before_losses=row['goods'] + row['source_goods'],
                goods_not_sold_after_quality_return_and_qualification=row['goods'] + row['source_goods'] - sold)
    if any(not math.isfinite(v) or v < -0.01 for v in mass.values()):
        raise ValueError('Invalid purchased mass')
    return dict(year=row['year'], population=row['population'], mass_kg=mass,
                revenue_by_sector_USD2041=revenues, total_invoices_USD2041=row['gross_USD2041'])


def checkpoint(rows, horizon):
    selected = [r for r in rows if r['year'] <= horizon]
    if not selected or selected[-1]['year'] != horizon:
        raise ValueError('Missing purchase checkpoint')
    return dict(horizon=horizon, first_sales_year=selected[0]['year'],
                annual_rows_counted=len(selected),
                cumulative_mass_kg={k: math.fsum(r['mass_kg'][k] for r in selected)
                                    for k in selected[0]['mass_kg']},
                cumulative_revenue_by_sector_USD2041={k: math.fsum(r['revenue_by_sector_USD2041'][k] for r in selected)
                                                     for k in selected[0]['revenue_by_sector_USD2041']},
                cumulative_total_invoices_USD2041=math.fsum(r['total_invoices_USD2041'] for r in selected),
                annual_at_checkpoint=selected[-1])


def product():
    manifest, orders, cfg, q, profiles, finance = inputs()
    if orders['schema'] != 'terluna.humanities.enclosed-orders/1':
        raise ValueError('Unexpected orders schema')
    revised = {c['name']: c for c in finance['primary_cases'] if c['variant'] == 'combined'}
    cases = []
    for c in orders['cases']:
        years = [r['year'] for r in c['annual']]
        if years != list(range(cfg['first_productive_year'], cfg['last_year'] + 1)):
            raise ValueError('Missing or duplicate annual order')
        rows = [annual(r, q['strategies'][c['portfolio']], cfg, q, profiles['inflation_factor'])
                for r in c['annual']]
        points = [checkpoint(rows, h) for h in manifest['horizons']]
        for point in points:
            matches = [r for r in revised[c['name']]['checkpoints'] if r['horizon'] == point['horizon']]
            if len(matches) != 3:
                raise ValueError('Expected all three financing comparisons')
            for r in matches:
                close(point['cumulative_total_invoices_USD2041'], r['cumulative_gross_receipts'])
                close(point['annual_at_checkpoint']['total_invoices_USD2041'], r['annual_gross'])
        cases.append(dict(name=c['name'], portfolio=c['portfolio'], market=c['market'],
                          market_growth_name=c['market_growth_name'], market_growth=c['market_growth'],
                          label=q['strategies'][c['portfolio']]['label'], checkpoints=points))
    if len(cases) != 81 or len({c['name'] for c in cases}) != 81 or len(revised) != 81:
        raise ValueError('Expected 81 unique primary alternatives')
    return dict(schema='terluna.humanities.earth-purchases/1',
                producer=dict(path='research/studies/enclosed_returns/earth_purchases.py', sha256=digest(__file__)),
                dependencies={'earth_purchases_inputs.json': digest(HERE / 'earth_purchases_inputs.json'),
                              **manifest['dependencies']},
                evidence=manifest['evidence'], units=manifest['units'], rules=manifest['rules'],
                reconciliation=dict(financial_checkpoints_compared=81 * len(manifest['horizons']) * 3,
                                    financing_modes_add_no_purchase_revenue=True), cases=cases)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=HERE / 'results/earth_purchases.json.gz')
    args = parser.parse_args()
    data = product()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(gzip.compress((json.dumps(data, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n').encode(), mtime=0))
    print(f"{len(data['cases'])} purchase alternatives; {data['reconciliation']['financial_checkpoints_compared']} financial checkpoints reconciled")
    print(args.output)


if __name__ == '__main__':
    main()
