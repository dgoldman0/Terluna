"""Six-route visitor spending, service output and external trade attribution.

This is a transparent service-only scenario screen, not a complete Tourism
Satellite Account or a demand forecast. Visitor status is assumed: interworld
commuters, migrants and people employed at the destination must be excluded by
the caller. Local amenity consumers remain a separate group. Their enjoyment
of free public or ecological amenities is not valued by this spending model.

All money uses one common real unit U; budgets and capacities are annual. Stay
nights and destination capacity use the same caller-chosen night unit. Routes
represent residents' round-trip itineraries, so access cost includes return
travel. Route capacities are allocated itinerary quotas, not a solved shared
vehicle/fleet constraint. Destination night capacities are available to these
visitors after other uses; annual utilization does not establish peak capacity.

The model assumes service receipts equal output at basic prices: product taxes,
retail merchandise, booking commissions and packages need separate treatment
before entering these inputs. Direct GVA plus domestic intermediate inputs
plus explicitly imported inputs equals service output. Carrier production
follows operator residence, independently of origin and destination. First-tier
import purchases are credited to their named supplier world; indirect supplier
GVA, further imported content, foreign profits, depreciation, induced effects
and fiscal/project cash collection are not estimated.

By default all output is attributed within existing accounts. No receipts or
GVA are added to baseline leisure, transport or world product. An optional
incrementality specification identifies additional fractions and displaced
GVA, but still returns only a conditional adjustment for the caller to assess.

Sources read 2026-10-31: UN Tourism Satellite Account RMF2008, paras2.6-2.8,
4.81-4.98 and Annex6 linkage/leakage discussion; ABS TSA2024-25 methodology,
visitor, direct/indirect output and GVA definitions. This was targeted reading,
not full-document review. Exact URLs and scope are returned in SOURCE_NOTES.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from numbers import Real


WORLDS = ("earth", "moon", "array")
COMPONENTS = ("visitor_services", "local_amenities", "transport")
SOURCE_NOTES = (
    {
        "url": "https://unstats.un.org/unsd/publication/Seriesf/SeriesF_80rev1e.pdf",
        "accessed": "2026-10-31",
        "scope": "Targeted read: paras2.6-2.8,4.81-4.98; Annex6 linkage and leakage discussion. Direct tourism GVA, local residents versus visitors, and supply-chain boundaries.",
    },
    {
        "url": "https://www.abs.gov.au/methodologies/tourism-satellite-account-methodology/2024-25",
        "accessed": "2026-10-31",
        "scope": "Targeted read: visitor/usual-environment, direct and indirect tourism output, and tourism GVA definitions. Tourism is already embedded in core national accounts.",
    },
)


def _keys(value, expected, name):
    if not isinstance(value, Mapping) or set(value) != set(expected):
        raise ValueError(f"{name} must contain exactly {sorted(expected)}")


def _number(value, name, *, positive=False, signed=False, share=False):
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite number")
    value = float(value)
    if not math.isfinite(value) or (not signed and value < 0) or (positive and value <= 0):
        raise ValueError(f"{name} is outside its allowed finite range")
    if share and value > 1:
        raise ValueError(f"{name} must not exceed one")
    return value


def _coefficients(value, producer, component):
    _keys(value, {"direct_gva_share", "imported_input_shares"}, f"{producer}.{component}")
    gva = _number(value["direct_gva_share"], "direct_gva_share", share=True)
    imports = value["imported_input_shares"]
    if not isinstance(imports, Mapping):
        raise ValueError("imported_input_shares must map suppliers to shares of output")
    clean = {}
    for supplier, fraction in imports.items():
        if supplier not in (*WORLDS, "outside") or supplier == producer:
            raise ValueError("Imported inputs require another world or 'outside' as supplier")
        clean[supplier] = _number(fraction, "imported input share", share=True)
    total_imports = math.fsum(clean.values())
    if gva + total_imports > 1 + 1e-12:
        raise ValueError("Direct GVA plus imported input shares exceed output")
    return {"direct_gva_share": gva, "imported_input_shares": clean,
            "domestic_input_share": max(0.0, 1 - gva - total_imports)}


def _finite_output(value):
    if isinstance(value, Mapping):
        for child in value.values():
            _finite_output(child)
    elif isinstance(value, (list, tuple)):
        for child in value:
            _finite_output(child)
    elif isinstance(value, Real) and not isinstance(value, bool) and not math.isfinite(value):
        raise ValueError("Tourism calculation exceeded finite numerical range")


def evaluate(worlds, routes, *, baseline_leisure_gva=None, incremental=None,
             required_net_exports=None):
    """Attribute six routes and local amenity spending without changing GDP.

    ``worlds[w]`` has outbound_budget_U, local_amenity_spend_U,
    inbound_visitor_night_capacity, and production. Production supplies one
    coefficient block for each of visitor_services, local_amenities, transport;
    each block specifies direct_gva_share and imported_input_shares by supplier.

    Each of six route dictionaries supplies origin, destination, budget_share,
    stay_nights, destination_spend_per_night_U, access_cost_per_trip_U,
    operator_shares and trip_capacity. Each origin's budget shares sum to one;
    operator shares also sum to one. Budget implies requested visits; route
    caps bind first, then destination nights are rationed proportionally across
    incoming routes. Unspent money is reported and is NOT automatically spent
    locally or rerouted. Fractional visits represent expected annual volumes.

    ``baseline_leisure_gva`` optionally supplies all three worlds' existing
    leisure-sector GVA. It is compared with visitor-service plus local-amenity
    GVA; transport is separate. This comparison does not establish an exact
    industrial mapping and never adds tourism to the baseline.

    ``incremental`` is off by default. If provided it must supply all three
    worlds and components, each with additional_fraction and displaced_gva_U.
    Potential incremental GVA is attributed direct GVA times that fraction,
    minus displaced GVA. Budgets, available factors and price responses remain
    assumptions; no multiplier or automatic GWP adjustment follows.

    ``required_net_exports`` optionally supplies the earlier annual accounts'
    signed net-export requirement for each world. The comparison uses net
    external tourism trade cash, including outbound visitor purchases, carrier
    residence and specified input imports. It is not a project profit, public
    tax receipt or additional credit to already-counted baseline production.
    """
    _keys(worlds, WORLDS, "worlds")
    model = {}
    for world, value in worlds.items():
        _keys(value, {"outbound_budget_U", "local_amenity_spend_U",
                      "inbound_visitor_night_capacity", "production"}, world)
        _keys(value["production"], COMPONENTS, world + ".production")
        model[world] = {
            field: _number(value[field], world + "." + field)
            for field in ("outbound_budget_U", "local_amenity_spend_U", "inbound_visitor_night_capacity")
        }
        model[world]["production"] = {
            component: _coefficients(value["production"][component], world, component)
            for component in COMPONENTS
        }
    if not isinstance(routes, (list, tuple)) or len(routes) != 6:
        raise ValueError("Exactly six directed interworld routes are required")
    expected_routes = {(a, b) for a in WORLDS for b in WORLDS if a != b}
    seen = set()
    route_rows = []
    route_fields = {"origin", "destination", "budget_share", "stay_nights",
                    "destination_spend_per_night_U", "access_cost_per_trip_U",
                    "operator_shares", "trip_capacity"}
    for route in routes:
        _keys(route, route_fields, "route")
        origin, destination = route["origin"], route["destination"]
        pair = (origin, destination)
        if pair not in expected_routes or pair in seen:
            raise ValueError("Routes must contain every distinct directed world pair once")
        seen.add(pair)
        share = _number(route["budget_share"], "budget_share", share=True)
        nights = _number(route["stay_nights"], "stay_nights", positive=True)
        nightly = _number(route["destination_spend_per_night_U"], "nightly spending")
        access = _number(route["access_cost_per_trip_U"], "access spending")
        capacity = _number(route["trip_capacity"], "trip_capacity")
        operator_shares = route["operator_shares"]
        if not isinstance(operator_shares, Mapping) or not operator_shares or set(operator_shares) - set(WORLDS):
            raise ValueError("operator_shares must name resident operators in the three worlds")
        operators = {w: _number(f, "operator share", share=True) for w, f in operator_shares.items()}
        if not math.isclose(math.fsum(operators.values()), 1, rel_tol=0, abs_tol=1e-10):
            raise ValueError("Operator shares must sum to one")
        trip_spend = nights * nightly + access
        if trip_spend <= 0 or not math.isfinite(trip_spend):
            raise ValueError("A budget-based route needs a positive finite trip price; free trips require a separate demand model")
        budget = model[origin]["outbound_budget_U"] * share
        requested = budget / trip_spend
        route_limited = min(requested, capacity)
        route_rows.append({
            "origin": origin, "destination": destination, "allocated_budget_U": budget,
            "budget_share": share, "stay_nights": nights, "destination_spend_per_night_U": nightly,
            "access_cost_per_trip_U": access, "total_spend_per_trip_U": trip_spend,
            "operator_shares": operators, "trip_capacity": capacity,
            "requested_visits": requested, "visits_after_route_capacity": route_limited,
            "route_capacity_rejected_visits": requested - route_limited,
        })
    for world in WORLDS:
        total = math.fsum(r["budget_share"] for r in route_rows if r["origin"] == world)
        if not math.isclose(total, 1, rel_tol=0, abs_tol=1e-10):
            raise ValueError("Each origin's route budget shares must sum to one")

    capacity = {}
    for world in WORLDS:
        requested_nights = math.fsum(r["visits_after_route_capacity"] * r["stay_nights"]
                                    for r in route_rows if r["destination"] == world)
        available = model[world]["inbound_visitor_night_capacity"]
        scale = min(1.0, available / requested_nights) if requested_nights else 1.0
        used = requested_nights * scale
        capacity[world] = {"visitor_nights_after_route_caps": requested_nights,
                           "available_visitor_nights": available, "destination_capacity_scale": scale,
                           "realized_visitor_nights": used,
                           "annual_capacity_utilization": used / available if available else None}

    production = {world: {component: {"gross_receipts_U": 0.0, "direct_gva_U": 0.0,
                                     "domestic_intermediate_inputs_U": 0.0,
                                     "imported_intermediate_inputs_U": 0.0}
                         for component in COMPONENTS} for world in WORLDS}
    flows = []

    def flow(payer, recipient, amount, category, kind, route_id=None):
        if amount:
            flows.append({"payer": payer, "recipient": recipient, "amount_U": amount,
                          "category": category, "kind": kind, "route": route_id})

    def serve(producer, component, receipts, route_id=None):
        coef = model[producer]["production"][component]
        row = production[producer][component]
        row["gross_receipts_U"] += receipts
        row["direct_gva_U"] += receipts * coef["direct_gva_share"]
        row["domestic_intermediate_inputs_U"] += receipts * coef["domestic_input_share"]
        imports = receipts * math.fsum(coef["imported_input_shares"].values())
        row["imported_intermediate_inputs_U"] += imports
        category = "local_amenities" if component == "local_amenities" else "tourism"
        for supplier, share in coef["imported_input_shares"].items():
            flow(producer, supplier, receipts * share, category, "imported_intermediate_input", route_id)

    origins = {w: {"outbound_budget_U": model[w]["outbound_budget_U"],
                   "realized_visits": 0.0, "resident_trip_spend_U": 0.0,
                   "outbound_destination_spend_U": 0.0,
                   "resident_operator_transport_spend_U": 0.0,
                   "foreign_operator_transport_spend_U": 0.0,
                   "local_amenity_spend_U": model[w]["local_amenity_spend_U"]}
               for w in WORLDS}
    for route in route_rows:
        origin, destination = route["origin"], route["destination"]
        route_id = origin + "_to_" + destination
        visits = route["visits_after_route_capacity"] * capacity[destination]["destination_capacity_scale"]
        destination_receipts = visits * route["stay_nights"] * route["destination_spend_per_night_U"]
        access_receipts = visits * route["access_cost_per_trip_U"]
        spend = destination_receipts + access_receipts
        route.update({"realized_visits": visits, "realized_visitor_nights": visits * route["stay_nights"],
                      "destination_capacity_scale": capacity[destination]["destination_capacity_scale"],
                      "destination_gross_receipts_U": destination_receipts,
                      "transport_gross_receipts_U": access_receipts,
                      "resident_trip_spend_U": spend,
                      "unspent_route_budget_U": max(0.0, route["allocated_budget_U"] - spend)})
        origin_row = origins[origin]
        origin_row["realized_visits"] += visits
        origin_row["resident_trip_spend_U"] += spend
        origin_row["outbound_destination_spend_U"] += destination_receipts
        flow(origin, destination, destination_receipts, "tourism", "visitor_service_purchase", route_id)
        serve(destination, "visitor_services", destination_receipts, route_id)
        for operator, share in route["operator_shares"].items():
            receipts = access_receipts * share
            origin_row["resident_operator_transport_spend_U" if operator == origin
                       else "foreign_operator_transport_spend_U"] += receipts
            flow(origin, operator, receipts, "tourism", "transport_purchase", route_id)
            serve(operator, "transport", receipts, route_id)
    for world, origin in origins.items():
        origin["unspent_outbound_budget_U"] = max(0.0, origin["outbound_budget_U"] - origin["resident_trip_spend_U"])
        origin["outbound_crossborder_visitor_spend_U"] = (origin["outbound_destination_spend_U"]
                                                         + origin["foreign_operator_transport_spend_U"])
        origin["outbound_budget_retained_domestically_U"] = (origin["unspent_outbound_budget_U"]
                                                              + origin["resident_operator_transport_spend_U"])
        flow(world, world, origin["local_amenity_spend_U"], "local_amenities", "resident_amenity_purchase")
        serve(world, "local_amenities", origin["local_amenity_spend_U"])

    trade = {}
    for world in (*WORLDS, "outside"):
        row = {}
        for category in ("tourism", "local_amenities"):
            credits = math.fsum(f["amount_U"] for f in flows if f["category"] == category
                                and f["recipient"] == world and f["payer"] != world)
            debits = math.fsum(f["amount_U"] for f in flows if f["category"] == category
                               and f["payer"] == world and f["recipient"] != world)
            row[category + "_external_receipts_U"] = credits
            row[category + "_external_payments_U"] = debits
            row[category + "_net_external_cash_U"] = credits - debits
        row["combined_net_external_cash_U"] = row["tourism_net_external_cash_U"] + row["local_amenities_net_external_cash_U"]
        trade[world] = row

    summaries = {}
    if baseline_leisure_gva is not None:
        _keys(baseline_leisure_gva, WORLDS, "baseline_leisure_gva")
    if required_net_exports is not None:
        _keys(required_net_exports, WORLDS, "required_net_exports")
    if incremental is not None:
        _keys(incremental, WORLDS, "incremental")
    for world in WORLDS:
        p = production[world]
        visitor_gva = p["visitor_services"]["direct_gva_U"]
        local_gva = p["local_amenities"]["direct_gva_U"]
        transport_gva = p["transport"]["direct_gva_U"]
        summary = {"production": p, "origin_budget": origins[world], "inbound_capacity": capacity[world],
                   "visitor_facing_direct_gva_U": visitor_gva,
                   "local_amenity_direct_gva_U": local_gva,
                   "tourism_transport_direct_gva_U": transport_gva,
                   "direct_tourism_gva_U": visitor_gva + transport_gva,
                   "leisure_facing_direct_gva_U": visitor_gva + local_gva,
                   "trade": trade[world], "baseline_gwp_adjustment_U": 0.0}
        if baseline_leisure_gva is not None:
            baseline = _number(baseline_leisure_gva[world], "baseline leisure GVA")
            measured = visitor_gva + local_gva
            summary["baseline_leisure_comparison"] = {
                "baseline_gva_U": baseline, "modeled_leisure_facing_gva_U": measured,
                "baseline_remainder_U": baseline - measured,
                "modeled_fraction_of_baseline": measured / baseline if baseline else None,
                "requires_reconciliation": measured > baseline,
                "treatment": "Attribution comparison only; transport is separate and no output is added to the baseline.",
            }
        if incremental is not None:
            _keys(incremental[world], COMPONENTS, world + ".incremental")
            component_changes = {}
            for component in COMPONENTS:
                assumption = incremental[world][component]
                _keys(assumption, {"additional_fraction", "displaced_gva_U"}, "incremental component")
                fraction = _number(assumption["additional_fraction"], "additional fraction", share=True)
                displaced = _number(assumption["displaced_gva_U"], "displaced GVA")
                component_changes[component] = {
                    "additional_fraction": fraction, "displaced_gva_U": displaced,
                    "potential_incremental_gva_U": p[component]["direct_gva_U"] * fraction - displaced,
                }
            summary["optional_incrementality"] = {
                "components": component_changes,
                "potential_incremental_gva_U": math.fsum(x["potential_incremental_gva_U"] for x in component_changes.values()),
                "status": "Conditional additional direct GVA less explicit displacement; no automatic GDP change or indirect multiplier. Net-product and financing effects need additional capital/cost accounts.",
            }
        if required_net_exports is not None:
            required = _number(required_net_exports[world], "required net exports", signed=True)
            tourism_net = trade[world]["tourism_net_external_cash_U"]
            combined_net = trade[world]["combined_net_external_cash_U"]
            summary["financing_comparison"] = {
                "required_net_exports_U": required,
                "tourism_net_external_cash_U": tourism_net,
                "tourism_fraction_of_positive_requirement": tourism_net / required if required > 0 else None,
                "remaining_net_exports_after_tourism_U": required - tourism_net,
                "remaining_net_exports_after_tourism_and_local_input_trade_U": required - combined_net,
                "status": "Contribution to an economy-wide trade requirement, not extra GDP, fiscal revenue or cash collectible by a project. Existing tourism trade must be reconciled with the rest of the baseline trade account.",
            }
        summaries[world] = summary
    partition_residuals = {
        world: math.fsum(p["gross_receipts_U"] - p["direct_gva_U"] - p["domestic_intermediate_inputs_U"]
                          - p["imported_intermediate_inputs_U"] for p in production[world].values())
        for world in WORLDS
    }
    result = {
        "schema": "terluna.world-economy.tourism-satellite-screen/1",
        "evidence": "Budget, price, service-coefficient and capacity scenarios. No demand estimation, market equilibrium, income ownership or empirical capacity validation.",
        "routes": route_rows, "worlds": summaries, "trade_flows": flows,
        "outside_trade": trade["outside"],
        "closure": {
            "service_output_partition_residuals_U": partition_residuals,
            "tourism_external_cash_including_outside_U": math.fsum(x["tourism_net_external_cash_U"] for x in trade.values()),
            "local_external_cash_including_outside_U": math.fsum(x["local_amenities_net_external_cash_U"] for x in trade.values()),
            "three_world_combined_net_external_cash_U": math.fsum(trade[w]["combined_net_external_cash_U"] for w in WORLDS),
        },
        "unspent_budget_rule": "Unused outbound funds are retained; they do not automatically create local amenity output. The domestic-retention field also includes purchases from resident transport operators before those operators' imported inputs.",
        "baseline_rule": "Satellite attribution only. Gross visitor receipts, direct tourism GVA and local amenity GVA are not additional to existing leisure/transport/world output. Optional incrementality is a separate scenario result.",
        "scope_limits": "Service-only basic-price screen; no product taxes, imported retail resale, commissions, foreign profit transfers, depreciation, project/fiscal cash flows, indirect GVA, induced multipliers, seasonality, crowding feedback or complete supply chains. First-tier imported-input trade is explicit; further supply-chain imports are omitted.",
        "sources": [dict(source) for source in SOURCE_NOTES],
    }
    _finite_output(result)
    return result
