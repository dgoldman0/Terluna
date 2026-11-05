"""Check the eel screen's mass boundaries, units and reproducibility."""
import json

import pytest

from shared.constants import JULIAN_YEAR_DAYS
from shared.provenance import constants_changed
from . import eel_budget as model


@pytest.fixture
def params():
    return json.loads((model.HERE / "eel_parameters.json").read_text())


def test_reference_cohort_hand_account(params):
    c = model.cohort(params["cohort"])
    assert c["harvest_individuals"] == 4000
    assert c["glass_eels_required"] == pytest.approx(4000 / .9 / .8)
    assert c["eggs_before_fertilization_required"] == pytest.approx(4000 / .9 / .8 / .01 / .6 / .7)
    assert c["nursery_dry_feed_kg"] == pytest.approx(52)
    assert c["growout_dry_feed_kg"] == pytest.approx(1433 + 1 / 3)
    assert c["cooked_fish_kg"] == 520
    assert c["cooked_portions"] == 5200
    assert c["post_glass_dry_feed_per_kg_cooked"] == pytest.approx(1485.3333333333333 / 520)


def test_processing_and_nutrient_mass_close(params):
    c = model.cohort(params["cohort"])
    assert c["cooked_fish_kg"] + c["cooking_loss_kg"] + c["processing_coproducts_and_trim_kg"] == c["harvest_live_kg"]
    n = c["nutrients_kg"]
    assert n["nonretained_n"] == pytest.approx(85.3056)
    assert n["nonretained_p"] == pytest.approx(10.857777777777779)
    for element in ("n", "p"):
        assert n[f"feed_{element}"] + n[f"seed_{element}"] == pytest.approx(n[f"harvest_{element}"] + n[f"nonretained_{element}"])


def test_cohort_scales_but_intensity_does_not(params):
    a = model.cohort(params["cohort"])
    b = model.cohort(dict(params["cohort"], harvest_live_kg=2000))
    for k in ("glass_eels_required", "eggs_before_fertilization_required", "post_glass_dry_feed_kg", "cooked_fish_kg"):
        assert b[k] == pytest.approx(2 * a[k])
    assert b["post_glass_dry_feed_per_kg_cooked"] == a["post_glass_dry_feed_per_kg_cooked"]


def test_larval_bottleneck_changes_eggs_not_post_glass_feed(params):
    a = model.cohort(params["cohort"])
    b = model.cohort(dict(params["cohort"], larva_to_glass_survival=.001))
    assert b["eggs_before_fertilization_required"] == pytest.approx(10 * a["eggs_before_fertilization_required"])
    for k in ("glass_eels_required", "post_glass_dry_feed_kg", "cooked_fish_kg"):
        assert b[k] == a[k]


@pytest.mark.parametrize("key,value", [
    ("growout_survival", 0), ("nursery_survival", 1.01),
    ("harvest_live_kg", float("nan")), ("cut_yield_of_live", -1),
    ("cook_yield_of_cut", 1.1), ("glass_eel_individual_kg", .02),
    ("growout_survival", .001),  # Incoming biomass exceeds outgoing biomass.
    ("feed_crude_protein_fraction_dry", .001),  # Impossible N retention.
    ("feed_phosphorus_fraction_dry", .00001),  # Impossible P retention.
])
def test_reject_impossible_cohort_inputs(params, key, value):
    with pytest.raises(ValueError):
        model.cohort(dict(params["cohort"], **{key: value}))


def test_oxygen_units_and_nonoverlapping_demands(params):
    o = model.oxygen(params["oxygen"])
    assert o["fish_o2_kg_per_day"] == pytest.approx(4.8)
    assert o["nitrification_o2_kg_per_day"] == pytest.approx(1.371)
    assert o["other_biological_o2_kg_per_day"] == pytest.approx(2)
    assert o["total_o2_kg_per_day"] == pytest.approx(8.171)
    assert o["capacity_o2_kg_per_hour"] == pytest.approx(.5106875)
    assert o["nitrification_alkalinity_kg_caco3_per_day"] == pytest.approx(2.142)
    assert o["equivalent_flow_without_reoxygenation_m3_per_hour"] == pytest.approx(255.34375)
    assert o["water_o2_inventory_above_floor_kg"] == pytest.approx(.2)
    assert o["ideal_inventory_minutes_at_capacity_demand"] == pytest.approx(.2 / .5106875 * 60)


def test_oxygen_boundaries_and_margin(params):
    a = model.oxygen(params["oxygen"])
    b = model.oxygen(dict(params["oxygen"], feed_dry_kg_per_day=0))
    assert b["total_o2_kg_per_day"] == b["fish_o2_kg_per_day"]
    assert b["nitrification_alkalinity_kg_caco3_per_day"] == 0
    c = model.oxygen(dict(params["oxygen"], inlet_do_mg_per_l=9))
    assert c["equivalent_flow_without_reoxygenation_m3_per_hour"] == pytest.approx(a["equivalent_flow_without_reoxygenation_m3_per_hour"] / 2)
    assert c["ideal_inventory_minutes_at_capacity_demand"] == pytest.approx(a["ideal_inventory_minutes_at_capacity_demand"] * 2)


@pytest.mark.parametrize("key,value", [
    ("inlet_do_mg_per_l", 5), ("capacity_factor", .8),
    ("feed_dry_kg_per_day", -1), ("nitrified_fraction_of_tan", 1.1),
])
def test_reject_impossible_oxygen_inputs(params, key, value):
    with pytest.raises(ValueError):
        model.oxygen(dict(params["oxygen"], **{key: value}))


def test_salt_inventory_and_recovery_have_distinct_boundaries(params):
    a = model.salt_inventory(35, params["hatchery"])
    b = model.salt_inventory(17.5, params["hatchery"])
    c = model.salt_inventory(35, dict(params["hatchery"], discharged_salt_recovery_fraction=1))
    assert a["salt_inventory_kg"] == 3500
    assert a["makeup_salt_after_recovery_kg_per_day"] == 35
    assert b["salt_inventory_kg"] == a["salt_inventory_kg"] / 2
    assert c["salt_inventory_kg"] == a["salt_inventory_kg"]
    assert c["makeup_salt_after_recovery_kg_per_day"] == 0
    assert model.salt_inventory(0, params["hatchery"])["salt_inventory_kg"] == 0


def test_hatchery_annualization_is_not_batch_output(params):
    c = model.cohort(params["cohort"])
    h = model.hatchery(params["hatchery"], c)
    assert h["expected_glass_eels_per_occupied_batch_slot_year"] == pytest.approx(1000 * .8 * JULIAN_YEAR_DAYS / 300)
    assert h["batch_slots_rounded_up"] == 6
    p = dict(params["hatchery"], successful_batch_fraction=.4)
    assert model.hatchery(p, c)["batch_slots_continuous_equivalent"] == pytest.approx(2 * h["batch_slots_continuous_equivalent"])


def test_product_reproduces_and_pins_only_used_inputs():
    product = model.results()
    assert model.OUT.read_text() == model.dumps(product)
    assert not constants_changed(product["producer"]["constants"])
    assert product["producer"]["constants"] == {"JULIAN_YEAR_DAYS": JULIAN_YEAR_DAYS}
    for name, sha in product["producer"]["files"].items():
        assert model.digest(model.ROOT / name) == sha
    manifest = product["producer"]["environment_manifest"]
    assert model.digest(model.ROOT / manifest["path"]) == manifest["sha256"]
    sources = json.loads((model.HERE / "eel_sources.json").read_text())["sources"]
    assert len({s["id"] for s in sources}) == len(sources)
