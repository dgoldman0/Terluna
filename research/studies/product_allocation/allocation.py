"""Product and processing-location screens on admitted resource/economic data.

These calculations identify price, freight, recovery and power thresholds. They
do not invent customer contracts, elemental domestic requirements or a globally
optimal production mix when those inputs are still unspecified.
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


def load():
    cfg = read(HERE / 'allocation_inputs.json')
    for path, expected in cfg['dependencies'].items():
        if digest(HERE / path) != expected:
            raise ValueError('Changed product-allocation dependency: ' + path)
    return cfg, read(HERE / 'inputs/assays.json'), read(HERE / 'inputs/flows.json'), read(HERE / 'inputs/recovery.json'), read(HERE / 'market_benchmarks.json'), read(HERE / 'inputs/basis.json')


def maximum_refining_cost(grade, earth_recovery, space_recovery, price,
                          freight, survival, earth_processing):
    """USD per kg of common starting concentrate, before either routing choice.

    Earth route: ship 1 kg; refine the surviving cargo on Earth.
    Space route: refine 1 kg; ship only its recovered target product.
    Both routes sell at the same target-product price. Freight is paid per kg
    dispatched; Earth treatment is paid per kg of concentrate that arrives.
    Common concentration costs cancel. Other coproducts remain outside this
    single-product threshold and must be reconciled in a joint plant account.
    """
    if not 0 < grade <= 1 or not 0 < earth_recovery <= 1 or not 0 < space_recovery <= 1 or not 0 < survival <= 1:
        raise ValueError('Invalid fraction in refining comparison')
    return (survival * grade * (space_recovery - earth_recovery) * price
            + freight * (1 - grade * space_recovery)
            + survival * earth_processing)


def fabrication_mass_penalty(active_fraction, freight):
    """Extra dispatch cost per kg of active material incorporated in a product."""
    if not 0 < active_fraction <= 1:
        raise ValueError('Invalid finished-product mass fraction')
    return (1 / active_fraction - 1) * freight


def qualified_product_surplus(price, freight, survival):
    """Maximum budget for all remaining costs per kg of dispatch-ready product."""
    return survival * price - freight


def product():
    cfg, assays, flows, recovery, market, basis = load()
    index = cfg['USD2025_to_USD2026_planning_index']
    survival = cfg['return_survival']
    prices = {p['id']: p['USD2025_per_product_kg'] * index for p in market['products']}
    rows = []
    for p in market['products']:
        price = prices[p['id']]
        rows.append(dict(id=p['id'], label=p['label'], product_basis=p['price_basis'],
                         USD2026_per_product_kg=price,
                         maximum_freight_before_any_other_cost_USD2026_kg=survival * price,
                         transport_screens=[dict(freight_USD2026_kg=f,
                             remaining_cost_budget_USD2026_per_kg_shipped=qualified_product_surplus(price, f, survival),
                             at_or_above_inherited_ground_floor=f >= cfg['inherited_ground_floor_USD2026_kg'])
                             for f in cfg['freight_USD2026_kg']]))

    # Treat measured composition and engineering recovery separately. Osmium
    # remains in the physical PGM total, with no invented price attached.
    normalized = []
    for profile in ('ci', 'kreep15386'):
        for level in ('low', 'middle', 'high'):
            process = math.prod(recovery['recovery_scenarios'][level].values())
            factor = recovery['grade_transfer_factors'][profile][('low', 'middle', 'high').index(level)]
            quantities = {e: cfg['reference_feed_kg'] * v['ppm'] * 1e-6 * factor * process
                          for e, v in assays['profiles'][profile]['elements'].items()}
            normalized.append(dict(profile=profile, scenario=level, feed_kg=cfg['reference_feed_kg'],
                                   element_kg=quantities,
                                   groups_kg={name: (sum(quantities[e] for e in elements)
                                                     if all(e in quantities for e in elements) else None)
                                              for name, elements in assays['group_definitions'].items()}))

    middle = next(r for r in normalized if r['profile'] == 'ci' and r['scenario'] == 'middle')
    circuit = []
    for price_factor in cfg['benchmark_price_multipliers']:
        for freight in cfg['freight_USD2026_kg']:
            contributions = {e: max(0., qualified_product_surplus(prices[e] * price_factor, freight, survival))
                             * middle['element_kg'][e] / cfg['reference_feed_kg']
                             for e in cfg['joint_circuit_priced_elements']}
            budget = sum(contributions.values())
            circuit.append(dict(price_multiplier=price_factor, freight_USD2026_kg=freight,
                                individual_contributions_USD2026_per_feed_kg=contributions,
                                maximum_joint_incremental_cost_USD2026_per_feed_kg=budget,
                                maximum_joint_incremental_cost_USD2041_per_feed_kg=budget * cfg['USD2026_to_USD2041']))

    raw = []
    for e in cfg['joint_circuit_priced_elements']:
        ppm = assays['profiles']['ci']['elements'][e]['ppm']
        r = math.prod(recovery['recovery_scenarios']['middle'].values())
        raw.append(dict(element=e, CI_grade_ppm=ppm,
                        raw_feed_freight_ceiling_USD2026_kg=ppm * 1e-6 * r * prices[e] * survival,
                        minimum_raw_grade_at_inherited_floor_fraction=cfg['inherited_ground_floor_USD2026_kg'] / (r * prices[e] * survival)))

    refining = [dict(target='Pt', concentrate_target_mass_fraction=g,
                     earth_recovery=cfg['refining']['earth_recovery'], space_recovery=r,
                     freight_USD2026_kg=f, earth_treatment_USD2026_per_concentrate_kg=e,
                     maximum_space_treatment_USD2026_per_concentrate_kg=maximum_refining_cost(
                         g, cfg['refining']['earth_recovery'], r, prices['Pt'], f, survival, e))
                for g in cfg['refining']['target_mass_fractions']
                for r in cfg['refining']['space_recoveries']
                for f in cfg['refining']['freight_USD2026_kg']
                for e in cfg['refining']['Earth_treatment_USD2026_per_concentrate_kg']]
    fabrication = [dict(product=p['product'], active_fraction=p['active_fraction'],
                        freight_USD2026_kg=f,
                        added_freight_USD2026_per_active_kg=fabrication_mass_penalty(p['active_fraction'], f))
                   for p in cfg['fabrication_examples'] for f in cfg['refining']['freight_USD2026_kg']]
    energy = [dict(architecture=p['architecture'], imported_mass_fraction=p['imported_mass_fraction'],
                   energy_kWh_per_kg=kwh,
                   extra_energy_cost_USD2026_per_product_kg=kwh * (p['levelized_USD2026_kWh'] - cfg['Earth_electricity_comparison_USD2026_kWh']))
              for p in basis['power'] for kwh in cfg['processing_energy_kWh_per_kg']]

    # Annual checkpoint quantities, not cumulative output and not committed sales.
    cases = {c['name']: c for c in flows['cases']}
    production = {p['id']: p['estimated_world_primary_production_2025_t'] * 1000
                  for p in market['products'] if p['estimated_world_primary_production_2025_t'] is not None
                  and p['id'] in middle['element_kg']}
    scale = []
    for case_name in cfg['physical_cases']:
        for row in cases[case_name]['checkpoints']:
            if row['year'] not in cfg['checkpoints']:
                continue
            for share in cfg['treated_feed_shares']:
                feed = row['source_feed_kg_year'] * share
                elements = {e: m * feed / cfg['reference_feed_kg'] for e, m in middle['element_kg'].items()}
                scale.append(dict(case=case_name, year=row['year'], treated_share=share,
                                  inherited_physical_feasible=row['physical_feasible'],
                                  original_source_feed_kg_year=row['source_feed_kg_year'], treated_feed_kg_year=feed,
                                  refined_element_kg_year=elements,
                                  ratio_to_2025_world_primary_output={e: elements[e] * survival / q for e, q in production.items()},
                                  receipts_booked_USD2041=0,
                                  reading_rule='Conditional recovery opportunity before domestic requirements, future stocks, priced refining capacity and product orders. Production ratios are scale diagnostics, not demand caps or market forecasts.'))

    purchase = read(HERE / '../enclosed_returns/results/earth_purchases.json.gz')
    invoices = []
    for c in purchase['cases']:
        if c['name'] in cfg['purchase_cases']:
            for point in c['checkpoints']:
                r = point['annual_at_checkpoint']
                invoices.append(dict(case=c['name'], year=point['horizon'],
                                     Earth_goods_kg_year=r['mass_kg']['Earth_purchased_goods'],
                                     Earth_goods_invoices_USD2041_year=r['revenue_by_sector_USD2041']['materials'] + r['revenue_by_sector_USD2041']['source_products']))
    return dict(schema='terluna.research.product-allocation/1',
                producer=dict(path='research/studies/product_allocation/allocation.py', sha256=digest(__file__)),
                dependencies={'allocation_inputs.json': digest(HERE / 'allocation_inputs.json'), **cfg['dependencies']},
                evidence=cfg['evidence'], boundary=cfg['boundary'],
                rules=cfg['rules'], market_products=rows, normalized_recovery=normalized,
                joint_circuit_budgets=circuit, raw_feed_screens=raw, refining_location=refining,
                fabrication_location=fabrication, power_cost_differences=energy,
                annual_opportunity=scale, existing_customer_requirements=invoices)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=HERE / 'results/allocation.json')
    args = parser.parse_args()
    data = product()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')
    print(f"{len(data['market_products'])} products; {len(data['refining_location'])} refining comparisons; {len(data['annual_opportunity'])} checkpoint opportunities")


if __name__ == '__main__':
    main()
