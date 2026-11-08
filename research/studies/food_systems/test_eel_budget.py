"""Check the eel screen's mass boundaries, units, per-person cases and reproducibility."""
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


def test_older_factors_are_read_from_the_older_product():
    f = model.older_factors()
    basis = json.loads(model.OLDER.read_text())["basis"]
    assert f["n2o_n_share_of_feed_n"] == basis["aquaponics"]["n2o_share_of_n"]
    assert f["closed_system_kwh_per_kg_live"] == basis["aquaponics"]["kwh_per_kg_fish"]
    assert f["waters_spare_phosphorus_kt_year"] == basis["seas_phosphorus"]["all_waters_ceiling_kt_year"]
    assert f["n2o_tg_n_per_ppb"] == basis["gases"]["n2o_tg_n_per_ppb"]
    assert min(f["n2o_n_share_of_feed_n"]) < max(f["n2o_n_share_of_feed_n"])
    assert min(f["waters_spare_phosphorus_kt_year"]) < max(f["waters_spare_phosphorus_kt_year"])


def test_per_person_reference_case_by_hand(params):
    f = model.older_factors()
    c = model.cohort(params["cohort"])
    case = model.per_person(8, 1.0, c, f, {"pond": 100.0}, [1.0, 2.0])
    live = 8e9 / (.65 * .8)                                  # kg live a year for 1 kg cooked each
    feed = live * 1485.3333333333333 / 1000
    feed_n = feed * .48 / 6.25
    assert case["live_harvest_mt_year"] == pytest.approx(live / 1e9)
    assert case["post_glass_dry_feed_mt_year"] == pytest.approx(feed / 1e9)
    assert case["feed_phosphorus_kt_year"] == pytest.approx(feed * .01 / 1e6)
    assert case["harvest_phosphorus_kt_year"] == pytest.approx(live * .004 / 1e6)
    assert case["nonretained_phosphorus_kt_year"] == pytest.approx(167.0427, rel=1e-5)   # the review's 167 kt
    spare = f["waters_spare_phosphorus_kt_year"]
    assert case["nonretained_phosphorus_over_waters_spare"] == pytest.approx(
        [case["nonretained_phosphorus_kt_year"] / max(spare), case["nonretained_phosphorus_kt_year"] / min(spare)])
    assert case["n2o_kt_n_year"] == pytest.approx([feed_n * s / 1e6 for s in f["n2o_n_share_of_feed_n"]])
    assert case["n2o_ppb_per_year"] == pytest.approx([x / (f["n2o_tg_n_per_ppb"] * 1e3) for x in case["n2o_kt_n_year"]])
    kwh = live * f["closed_system_kwh_per_kg_live"]
    assert case["closed_system_electricity_twh_year"] == pytest.approx(kwh / 1e9)
    assert case["closed_system_mean_power_gw"] == pytest.approx(kwh / (24 * JULIAN_YEAR_DAYS) / 1e6)
    assert case["open_pond_methane_mt_year"]["pond"] == pytest.approx(live / 1e3 * 100 / 1e9)
    assert case["share_of_heritage_meat_allowance"] == [.5, 1.0]
    low, high = case["feed_crop_land_km2"]
    assert low == pytest.approx(feed * f["feed_kcal_per_g"] * 1e3 * min(f["variety_factor"])
                                / (max(f["crop_kcal_m2_day"]) * JULIAN_YEAR_DAYS) / 1e6)
    assert low < high


def test_per_person_cases_scale_with_people_and_portion(params):
    f = model.older_factors()
    c = model.cohort(params["cohort"])
    a, b, d = (model.per_person(people, kg, c, f, {"pond": 1.0}, [1, 2]) for people, kg in ((4, 1), (8, 1), (8, 2)))
    for k in ("live_harvest_mt_year", "feed_phosphorus_kt_year", "nonretained_phosphorus_kt_year",
              "closed_system_water_km3_year", "feed_crop_land_km2", "n2o_kt_n_year"):
        assert b[k] == pytest.approx([2 * x for x in a[k]] if isinstance(a[k], list) else 2 * a[k])
        assert d[k] == pytest.approx([2 * x for x in b[k]] if isinstance(b[k], list) else 2 * b[k])


@pytest.mark.parametrize("people,kg", [(0, 1), (8, 0), (float("nan"), 1), (8, -1)])
def test_reject_impossible_per_person_inputs(params, people, kg):
    with pytest.raises(ValueError):
        model.per_person(people, kg, model.cohort(params["cohort"]), model.older_factors(), {}, [1, 2])


def test_pond_methane_pairs_flux_with_its_yield():
    hours = 24 * JULIAN_YEAR_DAYS
    assert model.pond_methane_kg_per_t({"ch4_mg_m2_h": .48, "yield_t_ha": 13.12}) == pytest.approx(.48 * hours / 100 / 13.12)
    assert model.pond_methane_kg_per_t({"ch4_mg_m2_h": 16.17, "yield_t_ha": 2.25}) == pytest.approx(630.0, rel=1e-3)
    with pytest.raises(ValueError):
        model.pond_methane_kg_per_t({"ch4_mg_m2_h": 1.0, "yield_t_ha": 0})


def test_spawning_run_returns_little_phosphorus_with_its_young(params):
    r = model.spawning_run(params["cohort"], [.1, .01])
    assert r["phosphorus_carried_down_kg_per_t_spawners"] == pytest.approx(4.0)
    culture, *wild = r["cases"]
    assert culture["glass_eel_to_spawner_survival"] == pytest.approx(.72)
    assert culture["phosphorus_share_carried_back"] == pytest.approx(.0002 / (.72 * .25))
    assert [w["phosphorus_share_carried_back"] for w in wild] == pytest.approx([.008, .08])
    with pytest.raises(ValueError):
        model.spawning_run(params["cohort"], [0])


def test_product_covers_every_per_person_case():
    product = model.results()
    cases = product["per_person_cases"]
    assert len(cases) == 3 * 2 * 3
    assert {c["processing"] for c in cases} == {"efficient_example", "reference_example", "higher_resource_example"}


def test_product_reproduces_and_pins_only_used_inputs():
    product = model.results()
    assert model.OUT.read_text() == model.dumps(product)
    assert not constants_changed(product["producer"]["constants"])
    assert product["producer"]["constants"] == {"JULIAN_YEAR_DAYS": JULIAN_YEAR_DAYS}
    assert list(product["producer"]["inputs"]) == ["research/studies/food_systems/results/food_systems.json"]
    for name, sha in {**product["producer"]["files"], **product["producer"]["inputs"]}.items():
        assert model.digest(model.ROOT / name) == sha
    sources = json.loads((model.HERE / "eel_sources.json").read_text())["sources"]
    assert len({s["id"] for s in sources}) == len(sources)
