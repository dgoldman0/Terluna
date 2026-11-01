"""Budget, capacity, residence and non-additive tourism accounting controls."""

import copy
import json

import pytest

from research.studies.world_economy import tourism


def _worlds():
    return {
        world: {
            "outbound_budget_U": 1000,
            "local_amenity_spend_U": 100,
            "inbound_visitor_night_capacity": 1000,
            "production": {component: {"direct_gva_share": 0.5, "imported_input_shares": {}}
                           for component in tourism.COMPONENTS},
        }
        for world in tourism.WORLDS
    }


def _routes():
    return [
        {"origin": origin, "destination": destination, "budget_share": 0.5,
         "stay_nights": 2, "destination_spend_per_night_U": 10,
         "access_cost_per_trip_U": 30, "operator_shares": {origin: 1}, "trip_capacity": 100}
        for origin in tourism.WORLDS for destination in tourism.WORLDS if origin != destination
    ]


def test_bidirectional_balanced_tourism_has_gross_receipts_but_no_net_trade_gain():
    result = tourism.evaluate(_worlds(), _routes())
    assert len(result["routes"]) == 6
    for world in tourism.WORLDS:
        row = result["worlds"][world]
        assert row["origin_budget"]["realized_visits"] == 20
        assert row["inbound_capacity"]["realized_visitor_nights"] == 40
        assert row["production"]["visitor_services"]["gross_receipts_U"] == 400
        assert row["direct_tourism_gva_U"] == 500
        assert row["local_amenity_direct_gva_U"] == 50
        assert row["trade"]["tourism_net_external_cash_U"] == 0
        assert row["baseline_gwp_adjustment_U"] == 0
        assert "optional_incrementality" not in row
    assert result["closure"]["tourism_external_cash_including_outside_U"] == 0


def test_carrier_residence_changes_trade_and_transport_gva_without_relocating_visitors():
    routes = _routes()
    for route in routes:
        route["operator_shares"] = {"array": 1}
    result = tourism.evaluate(_worlds(), routes)
    assert result["worlds"]["earth"]["trade"]["tourism_net_external_cash_U"] == -600
    assert result["worlds"]["moon"]["trade"]["tourism_net_external_cash_U"] == -600
    array = result["worlds"]["array"]
    assert array["trade"]["tourism_net_external_cash_U"] == 1200
    assert array["tourism_transport_direct_gva_U"] == 900
    assert array["visitor_facing_direct_gva_U"] == 200
    assert array["inbound_capacity"]["realized_visitor_nights"] == 40
    assert array["origin_budget"]["resident_operator_transport_spend_U"] == 600
    assert result["worlds"]["earth"]["origin_budget"]["foreign_operator_transport_spend_U"] == 600


def test_outbound_budget_includes_resident_transport_that_is_not_an_import():
    result = tourism.evaluate(_worlds(), _routes())
    origin = result["worlds"]["moon"]["origin_budget"]
    assert origin["resident_trip_spend_U"] == 1000
    assert origin["outbound_crossborder_visitor_spend_U"] == 400
    assert origin["resident_operator_transport_spend_U"] == 600
    assert origin["outbound_budget_retained_domestically_U"] == 600


def test_destination_capacity_rations_all_inbound_routes_and_keeps_unspent_budget():
    worlds = _worlds()
    worlds["moon"]["inbound_visitor_night_capacity"] = 10
    result = tourism.evaluate(worlds, _routes())
    moon_routes = [r for r in result["routes"] if r["destination"] == "moon"]
    assert [r["realized_visits"] for r in moon_routes] == [2.5, 2.5]
    assert result["worlds"]["moon"]["inbound_capacity"]["annual_capacity_utilization"] == 1
    for origin in ("earth", "array"):
        budget = result["worlds"][origin]["origin_budget"]
        assert budget["resident_trip_spend_U"] == 625
        assert budget["unspent_outbound_budget_U"] == 375
        assert budget["local_amenity_spend_U"] == 100
        assert budget["outbound_budget_retained_domestically_U"] == 750


def test_route_capacity_is_applied_before_night_capacity_without_rerouting():
    routes = _routes()
    for r in routes:
        if (r["origin"], r["destination"]) == ("earth", "moon"):
            r["trip_capacity"] = 0
    result = tourism.evaluate(_worlds(), routes)
    earth_routes = [r for r in result["routes"] if r["origin"] == "earth"]
    assert sum(r["realized_visits"] for r in earth_routes) == 10
    assert result["worlds"]["earth"]["origin_budget"]["unspent_outbound_budget_U"] == 500
    assert result["worlds"]["moon"]["inbound_capacity"]["realized_visitor_nights"] == 20


def test_first_tier_import_leakage_has_explicit_outside_counterpart_and_no_second_gva_deduction():
    worlds = _worlds()
    for world in tourism.WORLDS:
        worlds[world]["production"]["visitor_services"]["imported_input_shares"] = {"outside": 0.1}
        worlds[world]["production"]["local_amenities"]["imported_input_shares"] = {"outside": 0.2}
    result = tourism.evaluate(worlds, _routes())
    for world in tourism.WORLDS:
        row = result["worlds"][world]
        assert row["visitor_facing_direct_gva_U"] == 200
        assert row["production"]["visitor_services"]["imported_intermediate_inputs_U"] == 40
        assert row["production"]["visitor_services"]["domestic_intermediate_inputs_U"] == 160
        assert row["trade"]["tourism_net_external_cash_U"] == -40
        assert row["trade"]["local_amenities_net_external_cash_U"] == -20
    assert result["outside_trade"]["combined_net_external_cash_U"] == 180
    assert result["closure"]["three_world_combined_net_external_cash_U"] == -180
    assert result["closure"]["tourism_external_cash_including_outside_U"] == 0
    assert result["closure"]["local_external_cash_including_outside_U"] == 0


def test_named_input_supplier_receives_trade_credit_without_invented_indirect_gva():
    worlds = _worlds()
    worlds["moon"]["production"]["visitor_services"]["imported_input_shares"] = {"earth": 0.1}
    result = tourism.evaluate(worlds, _routes())
    assert result["worlds"]["earth"]["trade"]["tourism_net_external_cash_U"] == 40
    assert result["worlds"]["moon"]["trade"]["tourism_net_external_cash_U"] == -40
    assert result["worlds"]["earth"]["direct_tourism_gva_U"] == 500
    assert result["closure"]["three_world_combined_net_external_cash_U"] == 0


def test_leisure_baseline_comparison_excludes_transport_and_never_adds_receipts():
    result = tourism.evaluate(_worlds(), _routes(), baseline_leisure_gva=dict.fromkeys(tourism.WORLDS, 200))
    for row in result["worlds"].values():
        compare = row["baseline_leisure_comparison"]
        assert compare["modeled_leisure_facing_gva_U"] == 250
        assert compare["baseline_remainder_U"] == -50
        assert compare["requires_reconciliation"]
        assert row["baseline_gwp_adjustment_U"] == 0


def test_explicit_incrementality_accounts_for_displaced_gva_and_can_be_negative():
    assumptions = {
        world: {component: {"additional_fraction": 0.5, "displaced_gva_U": 20}
                for component in tourism.COMPONENTS}
        for world in tourism.WORLDS
    }
    assumptions["moon"]["visitor_services"]["displaced_gva_U"] = 1000
    result = tourism.evaluate(_worlds(), _routes(), incremental=assumptions)
    assert result["worlds"]["earth"]["optional_incrementality"]["potential_incremental_gva_U"] == 215
    assert result["worlds"]["moon"]["optional_incrementality"]["potential_incremental_gva_U"] == -765
    assert result["worlds"]["earth"]["baseline_gwp_adjustment_U"] == 0


def test_financing_comparison_uses_net_trade_not_gross_receipts_or_direct_gva():
    routes = _routes()
    for route in routes:
        route["operator_shares"] = {"moon": 1}
    result = tourism.evaluate(_worlds(), routes,
                              required_net_exports={"earth": -1000, "moon": 2000, "array": -1000})
    compare = result["worlds"]["moon"]["financing_comparison"]
    assert compare["tourism_net_external_cash_U"] == 1200
    assert compare["tourism_fraction_of_positive_requirement"] == 0.6
    assert compare["remaining_net_exports_after_tourism_U"] == 800
    assert result["worlds"]["earth"]["financing_comparison"]["tourism_fraction_of_positive_requirement"] is None


def test_zero_budgets_and_zero_night_capacity_have_finite_nonnegative_results():
    worlds = _worlds()
    for value in worlds.values():
        value["outbound_budget_U"] = 0
        value["inbound_visitor_night_capacity"] = 0
    result = tourism.evaluate(worlds, _routes())
    assert json.dumps(result, allow_nan=False)
    assert all(r["realized_visits"] == 0 for r in result["routes"])
    assert all(row["direct_tourism_gva_U"] == 0 for row in result["worlds"].values())
    assert all(row["local_amenity_direct_gva_U"] == 50 for row in result["worlds"].values())


def test_inputs_are_unchanged_and_output_partitions_close():
    worlds, routes = _worlds(), _routes()
    before = copy.deepcopy((worlds, routes))
    result = tourism.evaluate(worlds, routes)
    assert (worlds, routes) == before
    assert all(abs(v) < 1e-10 for v in result["closure"]["service_output_partition_residuals_U"].values())


@pytest.mark.parametrize("failure", ["missing_route", "duplicate_route", "budget_shares", "operator_shares", "gva_imports", "domestic_import", "negative_capacity", "nan_budget"])
def test_invalid_boundaries_shares_and_capacities_fail_explicitly(failure):
    worlds, routes = _worlds(), _routes()
    if failure == "missing_route": routes.pop()
    if failure == "duplicate_route": routes[-1] = copy.deepcopy(routes[0])
    if failure == "budget_shares": routes[0]["budget_share"] = 0.6
    if failure == "operator_shares": routes[0]["operator_shares"] = {"earth": 0.9}
    if failure == "gva_imports": worlds["moon"]["production"]["visitor_services"]["imported_input_shares"] = {"earth": 0.6}
    if failure == "domestic_import": worlds["moon"]["production"]["visitor_services"]["imported_input_shares"] = {"moon": 0.1}
    if failure == "negative_capacity": routes[0]["trip_capacity"] = -1
    if failure == "nan_budget": worlds["earth"]["outbound_budget_U"] = float("nan")
    with pytest.raises(ValueError):
        tourism.evaluate(worlds, routes)
