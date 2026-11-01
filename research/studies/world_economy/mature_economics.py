"""Component-level delivered costs and bounded recurring demand.

All monetary inputs are one consistent real unit U. One production attempt
targets one specified product unit; its accepted fraction determines saleable
output. Component prices include their own operation and capital service. The
extra annual site cost (for example a centrifuge service) must be disjoint.
These functions infer no cost advantage from gravity, market share from cost,
or new world GDP from replacement spending. They do not model project finance.
"""
from __future__ import annotations

import math
from numbers import Real


def _positive(name, value, *, zero=True):
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f'{name} must be numeric')
    result = float(value)
    if not math.isfinite(result) or result < 0 or (not zero and result == 0):
        raise ValueError(f'{name} must be finite and ' + ('nonnegative' if zero else 'positive'))
    return result


def _fraction(name, value, *, zero=True):
    value = _positive(name, value, zero=zero)
    if value > 1:
        raise ValueError(f'{name} must not exceed one')
    return value


def component_cost(base_attempt_cost, shares, factors, accepted_yield,
                   annual_output, annual_site_cost=0, delivery_cost_per_good=0):
    """Levelized cost per accepted unit at the supplied output quantity.

    Each share is a disjoint fraction of the reference attempt cost. A factor
    modifies only that component. For example pressure_gas=1.2 on a 1% share
    adds 0.2%, not 20%, to reference cost before yield. Shell and cycle factors
    belong to distinct shares. A low-g saving is an independent assumed factor.
    Costs of failed attempts are paid; no scrap credit is silently inferred.
    annual_output is qualified units/year, not gross attempts/year. No ramp,
    demand, equipment-capacity or financing feasibility is inferred here.
    """
    base = _positive('base_attempt_cost', base_attempt_cost, zero=False)
    output = _positive('annual_output', annual_output, zero=False)
    yield_fraction = _fraction('accepted_yield', accepted_yield, zero=False)
    annual_site = _positive('annual_site_cost', annual_site_cost)
    delivery = _positive('delivery_cost_per_good', delivery_cost_per_good)
    if not shares or set(factors)-set(shares):
        raise ValueError('Components required; factor names must match components')
    weights = {k:_fraction(f'{k} share', v) for k,v in shares.items()}
    if not math.isclose(math.fsum(weights.values()), 1, rel_tol=0, abs_tol=1e-12):
        raise ValueError('Disjoint component shares must sum to one')
    scales = {k:_positive(f'{k} factor', factors.get(k, 1)) for k in weights}
    attempt = {k:base*weights[k]*scales[k] for k in weights}
    qualified = {k:v/yield_fraction for k,v in attempt.items()}
    site_unit = annual_site/output
    factory = math.fsum(qualified.values())+site_unit
    return dict(base_attempt_cost_U=base, shares=weights, factors=scales,
        accepted_yield=yield_fraction, annual_qualified_output=output,
        required_annual_attempts=output/yield_fraction,
        rejected_annual_attempts=output*(1/yield_fraction-1),
        cost_per_attempt_U=math.fsum(attempt.values()),
        qualified_component_costs_U=qualified,
        annual_site_cost_U=annual_site,
        annual_site_cost_per_qualified_unit_U=site_unit,
        factory_cost_per_qualified_unit_U=factory,
        delivery_cost_per_qualified_unit_U=delivery,
        delivered_cost_per_qualified_unit_U=factory+delivery,
        annual_levelized_production_cost_U=factory*output,
        annual_delivery_cost_U=delivery*output,
        interpretation='Scenario cost at a fixed qualified output. Component capital service is included; no additional depreciation, debt or plant capital may be added without reconciling that boundary.')


def compare_locations(locations, reduction_component):
    """Compare Moon against both Earth and array and any additional recipes.

    Each input is a component_cost result for the same specified product and
    destination. Different accepted yields are allowed. One can include both
    array_freefall and array_rotating alongside an array case. All competitors
    are considered; no market share or realized selling price is inferred.

    Reduction thresholds change only the named lunar component, at fixed yield,
    output, annual site cost and delivery. Fractions use that component's original
    reference cost, not total factory cost. A threshold above one cannot be met
    by removing that component alone; negative thresholds need no reduction.
    """
    if not {'moon','earth','array'} <= set(locations):
        raise ValueError('Moon, Earth and array comparisons are required')
    moon = locations['moon']
    if reduction_component not in moon['shares']:
        raise ValueError('Unknown lunar reduction component')
    candidates = {k:_positive(f'{k} delivered cost',v['delivered_cost_per_qualified_unit_U'])
                  for k,v in locations.items() if k != 'moon'}
    cheapest = min(candidates, key=candidates.get)
    alternative = candidates[cheapest]
    lunar_factory = moon['factory_cost_per_qualified_unit_U']
    lunar_delivered = moon['delivered_cost_per_qualified_unit_U']
    reference_component = moon['base_attempt_cost_U']*moon['shares'][reduction_component]/moon['accepted_yield']
    additional_reduction = (lunar_delivered-alternative)/reference_component if reference_component else None
    total_reduction = (1-moon['factors'][reduction_component]+additional_reduction
                       if additional_reduction is not None else None)
    return dict(cheapest_alternative=cheapest,
        cheapest_alternative_delivered_cost_U=alternative,
        lunar_delivered_headroom_U=alternative-lunar_delivered,
        moon_strictly_cheaper_than_all_alternatives=lunar_delivered < alternative,
        maximum_lunar_delivery_cost_for_parity_U=alternative-lunar_factory,
        candidate_component=reduction_component,
        additional_reduction_fraction_of_reference_component_for_parity=additional_reduction,
        total_reduction_fraction_of_reference_component_for_parity=total_reduction,
        component_reduction_alone_can_reach_parity=(total_reduction <= 1 if total_reduction is not None else None),
        interpretation='A delivered-cost threshold, not a sale, share, product qualification or monopoly rent. Ownership does not change the operating producer territory.')


def renewal_demand(stock_units, service_life_years, eligible_fraction,
                   captured_share, qualified_capacity_year,
                   customer_budget_U_year, delivered_price_U_unit):
    """Bound annual sales from a maintained mature stock, with no stock growth.

    Lifetime, eligible product fraction, awarded/captured share, capacity and
    budget are explicit scenario assumptions. Replacement volume is already part
    of the mature economy; it does not create an additional GDP increment. The
    caller still must establish process qualification and delivered competition.
    Units may be kg, panels or other consistent physical/service units.
    """
    stock = _positive('stock_units', stock_units)
    life = _positive('service_life_years', service_life_years, zero=False)
    eligible = _fraction('eligible_fraction', eligible_fraction)
    capture = _fraction('captured_share', captured_share)
    capacity = _positive('qualified_capacity_year', qualified_capacity_year)
    budget = _positive('customer_budget_U_year', customer_budget_U_year)
    price = _positive('delivered_price_U_unit', delivered_price_U_unit, zero=False)
    renewal = stock/life
    eligible_orders = renewal*eligible
    awarded = eligible_orders*capture
    sold = min(awarded,capacity,budget/price)
    return dict(stock_units=stock,service_life_years=life,
        total_annual_renewal_units=renewal,
        eligible_annual_orders_units=eligible_orders,
        assumed_awarded_annual_orders_units=awarded,
        qualified_annual_capacity_units=capacity,
        budget_limited_annual_units=budget/price,
        realized_annual_units=sold, annual_customer_spending_U=sold*price,
        unused_customer_budget_U=budget-sold*price,
        remaining_eligible_orders_units=eligible_orders-sold,
        captured_share_assumption=capture,
        interpretation='Maintained-stock replacement with explicit award assumptions; spending overlaps existing replacement/CFC and is not added to mature GWP.')
