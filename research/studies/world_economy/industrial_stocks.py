"""A finite reserve compared with later purchases of the same service quantity.

All quantities use kg and one consistent real monetary unit U. This is a
deterministic sensitivity calculation, not a price forecast or a financing
model. Goods are acquired and prepositioned at t=0. In each subsequent year,
physical attrition affects both compatible and incompatible stock, compatibility
decay affects surviving compatible stock, and end-year demand draws compatible
stock. Annual storage is charged on opening physical stock. Unserved orders are
purchased from the alternative supplier at year end, preserving equal service.

Compatibility loss leaves recoverable physical material in storage. No recovery
process, salvage receipts, terminal liquidation or new GDP from resale is assumed.
The real discount rate represents the chosen comparison rate; actual borrowing,
interest payments, taxes, risk and buyer financing require separate accounts.
"""
from __future__ import annotations

import math
from numbers import Real


def _nonnegative(name, value):
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f'{name} must be a finite nonnegative number')
    value = float(value)
    if not math.isfinite(value) or value < 0:
        raise ValueError(f'{name} must be a finite nonnegative number')
    return value


def _fraction(name, value):
    value = _nonnegative(name, value)
    if value > 1:
        raise ValueError(f'{name} must lie between zero and one')
    return value


def _years(value, *, allow_zero=False):
    if isinstance(value, bool) or not isinstance(value, int) or value < (0 if allow_zero else 1):
        raise ValueError('horizon_years must be an integer in its allowed range')
    return value


def _demands(annual_demand, horizon_years):
    years = _years(horizon_years)
    if isinstance(annual_demand, Real):
        return [_nonnegative('annual_demand', annual_demand)] * years
    values = list(annual_demand)
    if len(values) != years:
        raise ValueError('Annual demand sequence must match horizon_years')
    return [_nonnegative('annual_demand', x) for x in values]


def required_initial_stock(annual_demand, horizon_years, physical_attrition=0,
                           compatibility_decay=0, initial_compatible_fraction=1):
    """Minimum starting kg that can serve every annual order without new stock.

    Backward recurrence implements q*C[t-1] - demand[t] = C[t], where
    q=(1-physical_attrition)*(1-compatibility_decay). This assumes homogeneous
    stock with one shared compatibility risk, not active sorting or refurbishment.
    Impossible or nonrepresentable requirements raise ValueError.
    """
    demand = _demands(annual_demand, horizon_years)
    attrition = _fraction('physical_attrition', physical_attrition)
    decay = _fraction('compatibility_decay', compatibility_decay)
    fraction = _fraction('initial_compatible_fraction', initial_compatible_fraction)
    if not any(demand):
        return 0.0
    survival = (1-attrition)*(1-decay)
    if survival == 0 or fraction == 0:
        raise ValueError('Positive demand cannot be served by a reserve with zero usable survival')
    required = 0.0
    for amount in reversed(demand):
        required = (required + amount) / survival
        if not math.isfinite(required):
            raise ValueError('Required initial stock exceeds finite numeric range')
    required /= fraction
    if not math.isfinite(required):
        raise ValueError('Required initial stock exceeds finite numeric range')
    return required


def evaluate_stock(initial_stock, annual_demand, horizon_years,
                   acquisition_cost_per_kg, future_buy_cost_per_kg,
                   annual_price_decline=0, storage_cost_per_kg_year=0,
                   storage_fraction_initial_acquisition_value=0,
                   real_discount_rate=0, preposition_cost_per_kg=0,
                   physical_attrition=0, compatibility_decay=0,
                   initial_compatible_fraction=1):
    """Return a closed physical ledger and equal-service present-value comparison.

    Future buy cost is the t=0 delivered price for the matching grade/product at
    the use site; year t price is p0*(1-decline)**t. Acquisition excludes initial
    prepositioning. Storage may be priced per opening physical kg-year OR as a
    fixed annual fraction of the initial acquisition value, not both. The latter
    remains charged for the whole horizon; it represents a contracted fixed cost.
    Terminal reference value is compatible kg times the alternative delivered
    replacement price. It is unrealized, excludes selling costs and is never a PV
    credit. Incompatible stock has unpriced recovery potential, not zero mass.
    """
    demand = _demands(annual_demand, horizon_years)
    stock = _nonnegative('initial_stock', initial_stock)
    acquisition = _nonnegative('acquisition_cost_per_kg', acquisition_cost_per_kg)
    alternative = _nonnegative('future_buy_cost_per_kg', future_buy_cost_per_kg)
    storage_unit = _nonnegative('storage_cost_per_kg_year', storage_cost_per_kg_year)
    storage_fraction = _nonnegative('storage_fraction_initial_acquisition_value', storage_fraction_initial_acquisition_value)
    if storage_unit and storage_fraction:
        raise ValueError('Select one storage charging basis')
    discount = _nonnegative('real_discount_rate', real_discount_rate)
    preposition = _nonnegative('preposition_cost_per_kg', preposition_cost_per_kg)
    attrition = _fraction('physical_attrition', physical_attrition)
    decay = _fraction('compatibility_decay', compatibility_decay)
    compatible_fraction = _fraction('initial_compatible_fraction', initial_compatible_fraction)
    decline = _fraction('annual_price_decline', annual_price_decline)
    settings = dict(initial_stock_kg=stock, annual_demand_kg=demand,
        horizon_years=horizon_years, acquisition_cost_U_per_kg=acquisition,
        future_buy_cost_t0_U_per_kg=alternative, annual_real_price_decline=decline,
        storage_cost_U_per_opening_physical_kg_year=storage_unit,
        annual_storage_fraction_initial_acquisition_value=storage_fraction,
        real_discount_rate=discount, preposition_cost_U_per_kg=preposition,
        annual_physical_attrition=attrition, annual_compatibility_decay=decay,
        initial_compatible_fraction=compatible_fraction)
    compatible = stock * compatible_fraction
    incompatible = stock - compatible
    rows = []
    for year, amount in enumerate(demand, 1):
        opening_compatible, opening_incompatible = compatible, incompatible
        opening = math.fsum((compatible, incompatible))
        lost_compatible = compatible * attrition
        lost_incompatible = incompatible * attrition
        compatible -= lost_compatible
        incompatible -= lost_incompatible
        newly_incompatible = compatible * decay
        compatible -= newly_incompatible
        incompatible += newly_incompatible
        usable_before_withdrawal = compatible
        withdrawal = min(amount, usable_before_withdrawal)
        compatible -= withdrawal
        shortfall = amount - withdrawal
        loss = math.fsum((lost_compatible, lost_incompatible))
        closing = math.fsum((compatible, incompatible))
        price = alternative * (1-decline)**year
        discount_factor = (1+discount)**(-year)
        storage_cost = (opening*storage_unit + stock*acquisition*storage_fraction)
        rows.append(dict(year=year, opening_physical_kg=opening,
            opening_compatible_kg=opening_compatible,
            opening_incompatible_kg=opening_incompatible,
            physical_loss_kg=loss, newly_incompatible_kg=newly_incompatible,
            compatible_before_withdrawal_kg=usable_before_withdrawal,
            demand_kg=amount, withdrawal_kg=withdrawal,
            just_in_time_shortfall_purchase_kg=shortfall,
            closing_physical_kg=closing, closing_compatible_kg=compatible,
            closing_incompatible_kg=incompatible,
            physical_closure_residual_kg=opening-loss-withdrawal-closing,
            demand_closure_residual_kg=amount-withdrawal-shortfall,
            alternative_delivered_price_U_per_kg=price,
            storage_cost_U=storage_cost,
            shortfall_purchase_U=shortfall*price,
            all_just_in_time_purchase_U=amount*price,
            pv_factor=discount_factor,
            pv_storage_U=storage_cost*discount_factor,
            pv_shortfall_purchase_U=shortfall*price*discount_factor,
            pv_all_just_in_time_purchase_U=amount*price*discount_factor))
    pv_storage = math.fsum(row['pv_storage_U'] for row in rows)
    pv_shortfall = math.fsum(row['pv_shortfall_purchase_U'] for row in rows)
    pv_jit = math.fsum(row['pv_all_just_in_time_purchase_U'] for row in rows)
    initial_purchase = stock*acquisition
    initial_delivery = stock*preposition
    pv_stock_strategy = math.fsum((initial_purchase, initial_delivery, pv_storage, pv_shortfall))
    # The cost coefficient includes fixed storage based on acquisition value.
    coefficient = stock*(1+storage_fraction*math.fsum(row['pv_factor'] for row in rows))
    pv_physical_storage = math.fsum(row['opening_physical_kg']*storage_unit*row['pv_factor'] for row in rows)
    breakeven = ((pv_jit-initial_delivery-pv_physical_storage-pv_shortfall)/coefficient
                 if coefficient else None)
    total_loss = math.fsum(row['physical_loss_kg'] for row in rows)
    total_withdrawal = math.fsum(row['withdrawal_kg'] for row in rows)
    terminal = math.fsum((compatible, incompatible))
    totals = dict(initial_purchase_U=initial_purchase, initial_preposition_U=initial_delivery,
        stock_withdrawals_kg=total_withdrawal,
        just_in_time_shortfall_purchases_kg=math.fsum(row['just_in_time_shortfall_purchase_kg'] for row in rows),
        physical_loss_kg=total_loss, terminal_physical_kg=terminal,
        terminal_compatible_kg=compatible, terminal_incompatible_kg=incompatible,
        physical_closure_residual_kg=stock-total_loss-total_withdrawal-terminal,
        pv_storage_U=pv_storage, pv_shortfall_purchases_U=pv_shortfall,
        pv_all_just_in_time_U=pv_jit, pv_stock_plus_shortfall_U=pv_stock_strategy,
        pv_saving_vs_just_in_time_U=pv_jit-pv_stock_strategy,
        break_even_acquisition_cost_U_per_kg=breakeven,
        terminal_compatible_market_reference_value_U=compatible*rows[-1]['alternative_delivered_price_U_per_kg'],
        terminal_liquidation_credit_in_pv_U=0.0,
        original_manufacture_gva_recounted_on_withdrawal_U=0.0)
    return dict(schema='terluna.research.industrial-stocks/1', settings=settings,
        annual=rows, totals=totals,
        conventions=dict(timing='t0 acquisition; annual opening-stock storage; attrition then compatibility decay then end-year orders',
            equal_service='JIT purchases cover every order not served by compatible stock; alternative supply and route capacity are assumed available',
            compatibility='Incompatible material remains physical stock; recovery and conversion are unpriced',
            terminal_value='Delivered replacement-price reference only; no identified buyer, salvage model or terminal PV credit',
            inventory_production='No GDP estimate; original manufacture is not counted again on withdrawal or resale',
            finance='PV comparison only; acquisition is paid at t0 and no free capital, loan or independent programme cash is inferred'))


def evaluate_reserve_horizons(annual_demand, reserve_horizons=(5, 20, 50),
                             analysis_years=50, **cost_and_loss_settings):
    """Compare candidate reserves over ONE shared analysis horizon and demand.

    A reserve is sized to supply its first N years, then later orders use JIT.
    Stranded incompatible material remains in storage for the analysis horizon.
    Returns the cheapest of the specified cases and a no-reserve comparator;
    this is a discrete scenario comparison, not a general optimum proof.
    """
    demands = _demands(annual_demand, analysis_years)
    horizons = list(reserve_horizons)
    if not horizons:
        raise ValueError('At least one reserve horizon is required')
    for value in horizons:
        _years(value, allow_zero=True)
        if value > analysis_years:
            raise ValueError('Reserve horizon cannot exceed analysis horizon')
    horizons = sorted(set([0] + horizons))
    cases = []
    for horizon in horizons:
        required = required_initial_stock(demands[:horizon], horizon,
            physical_attrition=cost_and_loss_settings.get('physical_attrition', 0),
            compatibility_decay=cost_and_loss_settings.get('compatibility_decay', 0),
            initial_compatible_fraction=cost_and_loss_settings.get('initial_compatible_fraction', 1)) if horizon else 0.0
        result = evaluate_stock(required, demands, analysis_years, **cost_and_loss_settings)
        cases.append(dict(reserve_years=horizon, result=result))
    cheapest = min(cases, key=lambda case: case['result']['totals']['pv_stock_plus_shortfall_U'])
    return dict(analysis_years=analysis_years, cases=cases,
        cheapest_tested_reserve_years=cheapest['reserve_years'],
        interpretation='Lowest PV among these discrete reserve choices for identical annual service; not an unrestricted optimum')
