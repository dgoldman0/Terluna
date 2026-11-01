"""Conserve photons and carbon; carry the rain channel and the sky rules; distinguish failed communities and double counting."""
import json
import math

import pytest

from research.studies.aerial_productivity import run as A


def inputs():
    return dict(zip(("sky", "aero", "old", "surface", "previous", "aerial", "gases", "ledgers", "power", "p"), A.load()))


def test_depth_does_not_multiply_the_photon_supply():
    one = A.photon_column(0.2)
    four_levels = A.photon_column(sum([0.05] * 4))
    assert one == four_levels
    assert one["intercepted_fraction"] == pytest.approx(1 - math.exp(-0.2))
    assert one["mean_organism_light_fraction"] * 0.2 == pytest.approx(one["intercepted_fraction"])
    assert A.photon_column(100)["intercepted_fraction"] <= 1
    assert A.photon_column(100, 0.5)["intercepted_fraction"] <= 0.5


def test_darker_lower_layers_are_distinguished_from_the_mean():
    column = A.photon_column(0.1, 0.71)
    assert column["bottom_plane_mean_light_fraction"] < column["mean_organism_light_fraction"] < 1
    assert A.photon_column(0)["intercepted_fraction"] == 0
    with pytest.raises(ValueError):
        A.photon_column(-1)


def test_stratified_classes_share_one_photon_stream():
    shares = dict(round_giant=0.3, lesser=0.2, sky_reef=0.5)
    for c in (0, 0.01, 0.2, 1, 10):
        column = A.photon_column(c, 0.71)
        light = A.class_light(column, shares, list(shares), "stratified")
        assert sum(c * shares[k] * v["mean"] for k, v in light.items()) == pytest.approx(column["intercepted_fraction"])
        assert all(v["bottom"] <= v["mean"] + 1e-12 for v in light.values())
        assert light["round_giant"]["mean"] >= light["sky_reef"]["mean"]


def test_dark_bodies_cannot_keep_producing_or_replace_their_dead():
    sky = inputs()["sky"]
    row = next(r for r in sky["food_web"]["production"] if r["all_gates"])
    out = A.body_budget(row, 0, 0.003, 0.5, 0.8)
    assert not out["feasible_at_mean_light"]
    assert not out["photolytic_h2_area_time_feasible"]
    assert out["accounted_organic_production_c"] is None
    assert out["remaining_assimilate_c"] < 0
    assert out["assimilate_balance_residual_c"] == pytest.approx(0)


def test_prescribed_renewal_is_not_the_unallocated_surplus():
    sky = inputs()["sky"]
    row = next(r for r in sky["food_web"]["production"] if r["all_gates"])
    a = A.body_budget(row, 1, 0.001, 0, 0.8)
    b = A.body_budget(row, 1, 0.001, 0.5, 0.8)
    assert a["accounted_organic_production_c"] == b["accounted_organic_production_c"]
    assert b["extra_tissue_c"] == pytest.approx(0.4 * b["remaining_assimilate_c"])
    assert a["extra_tissue_c"] == 0
    assert a["remaining_assimilate_c"] == pytest.approx(row["surplus_c"] - 0.001 * row["construction_c_kg_m2"])


def test_water_limited_surface_cannot_produce_without_light():
    budgets = dict(rainfed=10, wet_margins=2, coasts=1, oligotrophic_water=1)
    e = dict(rainfed_light_potential=0.8, wet_margins=0.65, coasts=0.6, oligotrophic_water=0.1)
    zero = A.surface_loss(budgets, 40, A.photon_column(0), e, rain_fraction=0.5)
    assert all(v == 0 for v in zero.values())
    dark = A.surface_loss(budgets, 40, A.photon_column(1000), e, rain_fraction=0.5)
    assert dark["homogenized_gt_c_year"] == pytest.approx(14)
    assert dark["persistent_opaque_gt_c_year"] == pytest.approx(14)


def test_rain_channel_reads_the_gcm_dimming_branch():
    d = inputs()
    rain = A.rain_channel(d["gases"], d["old"])
    run = d["gases"]["runs"]["gases"]
    assert rain["dimming_fraction"] == pytest.approx(run["added_dimming_percent"] / 100)
    # About 0.12 Gt C a year of rain-fed growth at the 1.69% dimming: some 0.07 Gt C per percent.
    assert 0.06 < rain["gt_c_per_percent_interception"] < 0.08
    assert rain["growth_fraction_per_unit_interception"] * rain["rain_fed_reference_gt_c_year"] / 100 \
        == pytest.approx(rain["gt_c_per_percent_interception"])


def test_rain_channel_adds_its_cost_once_to_each_surface_response():
    budgets = dict(rainfed=10, wet_margins=2, coasts=1, oligotrophic_water=1)
    e = dict(rainfed_light_potential=0.8, wet_margins=0.65, coasts=0.6, oligotrophic_water=0.1)
    column = A.photon_column(0.1)
    f = column["intercepted_fraction"]
    dry = A.surface_loss(budgets, 40, column, e)
    wet = A.surface_loss(budgets, 40, column, e, rain_fraction=0.5)
    # Light binds nowhere on this water-limited land, so the homogenized response adds the whole rain share ...
    assert wet["homogenized_gt_c_year"] - dry["homogenized_gt_c_year"] == pytest.approx(10 * 0.5 * f)
    # ... and the persistent response charges the rain only on the ground outside the footprints.
    assert wet["persistent_opaque_gt_c_year"] - dry["persistent_opaque_gt_c_year"] == pytest.approx((1 - f) * 10 * 0.5 * f)
    # Where the light already holds the land below its reduced water, the rain adds nothing.
    shaded = A.surface_loss(budgets, 10.5, A.photon_column(0.5), e, rain_fraction=0.01)
    unrained = A.surface_loss(budgets, 10.5, A.photon_column(0.5), e)
    assert shaded["homogenized_gt_c_year"] == pytest.approx(unrained["homogenized_gt_c_year"])


def test_sky_rules_scale_their_one_percent_basis():
    d = inputs()
    towns = d["ledgers"]["heights"]["sky_towns"]
    at_one = A.sky_rules(0.01, d["ledgers"], d["power"], 0.5)
    assert at_one["reef_passes_per_town_day"] == pytest.approx(towns["reef_passes_per_town_day"]["0.01"])
    at_two = A.sky_rules(0.02, d["ledgers"], d["power"], 0.5)
    assert at_two["reef_passes_per_town_day"] == pytest.approx(towns["reef_passes_per_town_day"]["0.02"], rel=0.01)
    beams = d["power"]["reserves"]["reefs"]["beam_columns"]["towns"]["operating"]["by_cover"]
    assert at_two["beam_switches_per_receiver_day"]["towns"] == pytest.approx(beams["0.02"]["switches_per_receiver_day"], rel=0.01)
    with pytest.raises(ValueError):
        A.sky_rules(0.1, d["ledgers"], d["power"], 0.4)


def test_microbe_column_mass_is_fixed_by_optics_not_habitat_depth():
    m = inputs()["p"]["microbes"]
    a = A.microbial_screen(100, 0.45, m)
    b = A.microbial_screen(100, 0.45, {**m, "column_depth_m": m["column_depth_m"] * 2})
    for x, y in zip(a, b):
        assert x["stock_gt_c"] == y["stock_gt_c"]
        assert x["cell_synthesis_gt_c_year"] == y["cell_synthesis_gt_c_year"]
        assert x["dry_concentration_mg_m3"] == 2 * y["dry_concentration_mg_m3"]
        assert x["growth_minus_removal_per_day"]["30.0"] < 0
        assert x["growth_minus_removal_per_day"]["100.0"] > 0


def test_recorded_product_regenerates_and_counts_the_sky_once():
    product = A.results()
    assert json.loads(A.OUT.read_text()) == product
    assert set(product["producer"]["inputs"]) == set(A.INPUTS)
    assert "research/decisions.md" not in product["producer"]["inputs"]
    for case in product["scenarios"]:
        sky = case["sky_accounted_gt_c_year"]
        for kind, value in case["additional_world_gt_c_year"].items():
            if sky is None:
                assert value is None
            else:
                assert value == pytest.approx(sky - case["surface_loss"][kind])
                assert value <= case["additional_world_light_only_gt_c_year"][kind] + 1e-12
    assert any(not r["whole_mixture_feasible_at_mean_light"] for r in product["scenarios"])
    assert any(r["whole_mixture_feasible_at_mean_light"] and not r["all_classes_feasible_at_bottom_plane_mean"]
               for r in product["scenarios"])


def test_phosphorus_losses_precede_one_shared_fallout_credit():
    d = inputs()
    product = A.results()
    case = next(r for r in product["scenarios"] if r["arrangement"] == "stratified"
                and r["column"]["habitat_fraction"] == 1 and r["column"]["area_loading"] == 0.1)
    old_net = [r["need_t_year_at_1pct"] * 10 for r in d["sky"]["phosphorus"]["sky_need"]
               if r["state"] == "mature" and r["retention_case"] == "earth_forest"]
    assert all(a > b for a, b in zip(case["mature_p_resupply_before_fallout_t_year"], old_net))
    assert all(c < total / 1000 for c, total in zip(case["ideal_interception_of_shared_p_fallout_t_year"],
                                                  d["aerial"]["baseline"]["p_fallout_analogue_kg_year"]))


def test_failing_bottom_layers_are_not_credited_as_a_whole_community():
    cases = A.results()["scenarios"]
    rejected = [c for c in cases if c["whole_mixture_feasible_at_mean_light"]
                and not c["all_classes_feasible_at_bottom_plane_mean"]]
    assert rejected
    for case in rejected:
        assert case["mean_budget_only_diagnostic_gt_c_year"] > 0
        assert case["sky_accounted_gt_c_year"] is None
        assert all(v is None for v in case["additional_world_gt_c_year"].values())
