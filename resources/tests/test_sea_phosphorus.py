"""Closure, limits and bindings of the seas' phosphorus screen and its waters product.

Where a value could only agree with itself, the test recomputes it from the inputs by its own closed form instead."""
import json
import math
import os

import numpy as np
import pytest

from resources import sea_phosphorus as sp
from resources import sea_phosphorus_waters as sw
from shared.constants import MOON_GM, MOON_RADIUS, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY
from shared.provenance import constants_changed

PRODUCT = json.loads(sp.OUT.read_text())
WATERS = json.loads(sw.OUT.read_text())


def bound_available():
    try:
        for key in sp.BOUND:
            sp.bound_bytes(key)
    except sp.MissingInput:
        return False
    return True


needs_bound = pytest.mark.skipif(not bound_available(), reason="bound inputs are not in this tree")


# ---- Bindings -------------------------------------------------------------------------------------------------------
def test_producer_inputs_and_constants_are_pinned():
    assert PRODUCT["schema"] == sp.SCHEMA
    for path, sha in {**PRODUCT["producer"]["files"], **PRODUCT["producer"]["inputs"]}.items():
        assert sp.digest(sp.ROOT / path) == sha, path
    assert set(PRODUCT["producer"]["inputs"]) == set(sp.INPUTS)
    assert not constants_changed(PRODUCT["producer"]["constants"])
    bound = PRODUCT["producer"]["bound_inputs"]
    assert set(bound) == set(sp.BOUND)
    for key, b in bound.items():
        pinned = sp.BOUND[key]
        assert (b["path_at_integration"], b["branch"], b["commit"], b["sha256"]) == \
            (pinned["path"], pinned["branch"], pinned["commit"], pinned["sha256"])


@needs_bound
def test_product_regenerates_from_bound_inputs():
    assert PRODUCT == json.loads(json.dumps(sp.run()))


@needs_bound
def test_bound_values_match_their_sources():
    people = sp.load_bound("ecology_people")
    sky = sp.load_bound("sky_ecology")
    ways = sp.load_bound("ways_of_living")
    loss = PRODUCT["land_loss"]
    assert loss["total_kg_yr"] == [x * 1e9 for x in people["capacity"]["land_p_loss_tg"]]
    assert loss["dissolved_kg_yr"] == sky["phosphorus"]["sea_to_land"]["dissolved_loss_kg_year"]
    # The dissolved loss is 0.003-0.008 g/m2 over the dry land, as the sky-ecology study derives it.
    dry = people["land"]["dry_land_km2"] * 1e6
    assert np.allclose([x / dry * 1e3 for x in loss["dissolved_kg_yr"]], [0.003, 0.008], rtol=1e-3)
    assert PRODUCT["people"]["eaten_kg_per_person"] == people["needs"]["practices"]["lunar_design"]["p_eaten_kg"]
    # The design target's per-billion loss is the ecology screen's 1-12% of the land's loss (table 4).
    share = PRODUCT["people"]["ecology_design_target"]["share_of_land_loss_per_billion"]
    assert round(share[0], 2) == 0.01 and round(share[1], 2) == 0.12
    # Earth's seabirds return 2.4-6.7% of the loss, as the sky-ecology study has it.
    sea_to_land = sky["phosphorus"]["sea_to_land"]
    flyers = PRODUCT["natural_cycle"]["flyers"]
    assert np.allclose(flyers["earth_seabirds_share_of_loss"], sorted(sea_to_land["earth_seabirds_share_of_loss"]))
    # The fish P content reproduces the sky-ecology study's food P against the land's loss.
    target = next(t for t in sea_to_land["targets"] if t["target"] == "earth_seabird_scale")
    assert target["food_kg_wet_year"][0] * sp.FISH_P_WET / loss["total_kg_yr"][1] == pytest.approx(
        target["food_p_share_of_seas_input"][0], rel=1e-9)
    assert PRODUCT["seas"]["discharge_m3_yr"] == ways["moon"]["drainage_discharge_to_sea_km3_year"] * 1e9
    assert PRODUCT["routes"]["sorption"]["pump_efficiency"] == ways["assumptions"]["inherited"]["pump_efficiency"]
    assert "4–8 billion people on the Moon itself" in sp.load_bound("provisioning_decisions")


def test_register_and_box_model_values_match():
    register = sp.load("research/studies/lunar_cycle_ecology/results/lunar_cycle_ecology.json")
    stock = register["return_limb"]["rock_stock_top_metre"]
    p2o5 = sp.LUNAR_SOIL_P2O5_WT_PCT
    ratio = lambda a, b: stock[a]["P_kg_m2"] / stock[b]["P_kg_m2"]       # noqa: E731
    assert ratio("mare_soil", "highland_soil") == pytest.approx(p2o5["mare"] / p2o5["highland"], rel=1e-3)
    assert ratio("kreep_rich", "mare_soil") == pytest.approx(p2o5["kreep_rich"] / p2o5["mare"], rel=1e-3)
    assert sp.P_IN_P2O5 == pytest.approx(2 * 30.974 / 141.94, rel=1e-4)
    box = sp.load("resources/results/phosphorus.json")
    edge = next(e for e in box["finite_reservoir_cases"][0]["edges"]
                if e["source"] == "sea" and e["destination"] == "sediment")
    assert PRODUCT["seas"]["residence_yr"]["box_model"] == 1 / edge["rate_per_year"] == 1000.0
    first = sp.load("resources/results/first_screen.json")
    assert PRODUCT["seas"]["volume_m3"] == first["inventories"]["water"]["seas_kg"] / 1000
    # Literature arithmetic: Xiao et al.'s 9,496 t from 1,250 km2; Earth's ocean mean near 2.25 uM.
    assert PRODUCT["routes"]["seaweed"]["uptake_kg_p_m2_yr"] == pytest.approx(7.6e-3, rel=1e-3)
    assert PRODUCT["seas"]["earth_ocean"]["mean_um"] == pytest.approx(2.25, abs=0.01)
    assert PRODUCT["statements"]["open_water"]["km3_per_tg_at_earth_ocean"] == pytest.approx(14339, rel=1e-3)
    assert sp.STRUVITE_PER_P == pytest.approx(245.41 / 30.974, rel=1e-4)


# ---- The waters product ---------------------------------------------------------------------------------------------
def test_waters_product_reproduces_drainage_and_stays_in_bounds():
    assert WATERS["schema"] == sw.SCHEMA
    drainage = json.loads(sw.DRAINAGE.read_text())
    for key, (mine, theirs) in WATERS["reproduced_drainage"].items():
        assert mine == theirs, key
    assert WATERS["reproduced_drainage"]["discharge_to_sea_m3s"][0] == drainage["rivers"]["discharge_to_sea_m3s"]
    for path, sha in {**WATERS["producer"]["files"], **WATERS["producer"]["inputs"]}.items():
        assert sw.digest(sw.ROOT / path) == sha, path
    for path, sha in WATERS["producer"]["ignored_inputs"].items():
        if (sw.ROOT / path).exists():
            assert sw.digest(sw.ROOT / path) == sha, path
    assert not constants_changed(WATERS["producer"]["constants"])
    water = WATERS["water_m3s"]
    assert water["land_runoff"] + water["lake_rain_less_evaporation"] == pytest.approx(water["to_common_level_water"])
    passage = list(WATERS["lake_passage_share_of_discharge"].values())
    land = list(WATERS["lake_passage_share_of_land_runoff"].values())
    for series in (passage, land):
        assert all(0 <= x <= 1 for x in series) and series == sorted(series, reverse=True)
    # Rain on the lakes counts as having crossed one, so the discharge shares sit above the land runoff's.
    assert all(d >= l for d, l in zip(passage, land))
    assert WATERS["land_runoff_reaching_common_level_m3s"] <= water["land_runoff"]
    lakes = WATERS["lakes"]
    # Joining across the seam merges split lakes; retention bodies are the spilling lakes, every one with an outflow.
    assert lakes["drainage_product"] - lakes["lakes_seam_joined"] == lakes["seam_joins"] > 0
    assert lakes["retention_bodies"] == lakes["retention_bodies_spilling"]
    assert lakes["closed_lakes"] + lakes["retention_bodies"] >= lakes["lakes_seam_joined"]
    assert lakes["closed_lake_km2_cleared_from_spilling_lakes"] <= lakes["closed_lake_km2"]
    drained_km2 = drainage["lakes"]["area_share"] * 4 * math.pi * (MOON_RADIUS / 1e3) ** 2
    assert lakes["spilling_lake_km2"] + lakes["closed_lake_km2"] == pytest.approx(drained_km2, rel=2e-3)
    # Retention cannot keep the runoff that crosses no lake: survival is at least the share that bypasses them.
    bypass = 1 - WATERS["lake_passage_share_of_land_runoff"]["tau_ge_0_yr"]
    assert bypass * 0.99 < WATERS["dissolved_p_survival"] < 1
    survival = WATERS["particulate_p_survival"]
    assert all(bypass * 0.99 < x < 1 for x in survival.values())
    for curve in sw.BRUNE_BETA:      # slower lunar settling lets more through on every curve
        assert survival[f"{curve}_moon"] > survival[f"{curve}_earth"]
    assert survival["upper_earth"] < survival["median_earth"] < survival["lower_earth"]
    assert WATERS["gravity_ratio"] == MOON_SURFACE_GRAVITY / STANDARD_GRAVITY
    # The per-body sums, computed by atlas label, add up to the tracer totals computed at the outlets.
    shares = sum(row["dissolved_p_share"] for row in WATERS["by_receiving_water"].values())
    assert shares == pytest.approx(WATERS["dissolved_p_survival"], rel=1e-9)
    discharge = sum(row["discharge_m3s"] for row in WATERS["by_receiving_water"].values())
    assert discharge == pytest.approx(water["to_common_level_water"], rel=1e-9)
    atlas = json.loads(sw.ATLAS.read_text())
    hyp = WATERS["hypsometry"]["all_common_level_water"]
    assert hyp["area_km2"] == pytest.approx(atlas["water_share"] * 4 * math.pi * (MOON_RADIUS / 1e3) ** 2, rel=2e-3)
    assert hyp["volume_km3"] * 1e9 * 1000 == pytest.approx(atlas["water_mass_kg"], rel=2e-3)
    assert sum(hyp["by_depth_m_km2"].values()) == pytest.approx(hyp["area_km2"])
    # The largest river ends in Orientale and the second in a 229 km2 crater water, as the note says.
    mouths = WATERS["largest_rivers_receiving_water"]
    assert mouths[0]["waters"][0]["name"] == "Mare Orientale"
    assert mouths[1]["waters"][0]["area_km2"] == pytest.approx(229.4, abs=0.5)


def test_lake_retention_relations():
    assert sw.lake_survival_dissolved(0.0) == 1.0
    assert sw.lake_survival_dissolved(1.0) == pytest.approx(0.5)
    assert sw.lake_survival_dissolved(100.0) == pytest.approx(1 / 11)
    # Brune's median curve in Siyam's form: C/I of 0.01 and 1 year trap 45% and 99%.
    beta = sw.BRUNE_BETA["median"]
    assert 1 - sw.lake_survival_particulate(0.01, beta) == pytest.approx(math.exp(-0.79), rel=1e-9)
    assert 1 - sw.lake_survival_particulate(1.0, beta) == pytest.approx(0.992, abs=1e-3)
    assert sw.lake_survival_particulate(0.0, beta) == 1.0
    assert sw.lake_survival_particulate(1.0, beta, 0.5) > sw.lake_survival_particulate(1.0, beta)
    for bad in ((-1.0, beta), (1.0, 0.0), (1.0, beta, 0.0), (1.0, beta, 1.5)):
        with pytest.raises(ValueError):
            sw.lake_survival_particulate(*bad)


@pytest.mark.skipif(not os.environ.get("TERLUNA_WATERS_REGENERATE")
                    or not all((sw.ROOT / p).exists() for p in sw.IGNORED),
                    reason="set TERLUNA_WATERS_REGENERATE=1 with the grid inputs restored (1-3 minutes, 4 GB)")
def test_waters_product_regenerates():
    assert WATERS == json.loads(json.dumps(sw.run()))


# ---- Phosphorus closes ----------------------------------------------------------------------------------------------
def test_lakes_and_seas_share_the_land_loss():
    loss, part = PRODUCT["land_loss"], PRODUCT["partition"]
    s_d = WATERS["dissolved_p_survival"]
    s_p = sorted(WATERS["particulate_p_survival"].values())
    for i in range(2):
        d = loss["dissolved_kg_yr"][i]
        assert part["to_seas_dissolved_kg_yr"][i] == pytest.approx(s_d * d)
        assert part["kept_by_lakes_dissolved_kg_yr"][i] == pytest.approx((1 - s_d) * d)
    particles = [t - d for t, d in zip(loss["total_kg_yr"], loss["dissolved_kg_yr"])]
    assert part["to_seas_particulate_kg_yr"] == pytest.approx([s_p[0] * particles[0], s_p[-1] * particles[1]])
    # Lakes and seas together hold the whole loss at each land-loss edge once their survival is matched.
    for i, s in ((0, s_p[-1]), (1, s_p[0])):
        seas = s_d * loss["dissolved_kg_yr"][i] + s * particles[i]
        assert seas + part["kept_by_lakes_kg_yr"][i] == pytest.approx(loss["total_kg_yr"][i])
    # The listed receiving waters, summed by name, account for what reaches the common level; waters taking under
    # 100 m3/s are left out of the list and hold under 0.1% of it.
    bodies = part["by_receiving_water"].values()
    for i in range(2):
        assert sum(b["dissolved_kg_yr"][i] for b in bodies) == pytest.approx(part["to_seas_dissolved_kg_yr"][i],
                                                                              rel=1e-3)
    assert part["kept_by_lakes_share_of_loss"][0] > 0.95 and part["to_seas_share_of_loss"][1] < 0.05
    implied, gates = loss["implied_erosion_kg_yr"], loss["conservation_erosion_kg_yr"]
    assert implied[0] < gates[1] and gates[0] < implied[1]


def test_sea_box_conserves_and_matches_closed_form():
    for b, tau in ((5e6, 1000.0), (5e7, 2e4), (0.0, 1e5)):
        box = sp.sea_box(b, tau, 1e16, initial_kg=1e9)
        for row in box["states"]:
            assert row["inventory_kg"] == pytest.approx(row["closed_form_kg"], rel=1e-9, abs=1e-3)
            assert row["inventory_kg"] + row["buried_kg"] == pytest.approx(1e9 + b * row["year"], rel=1e-9)
            assert row["unit_state"] == pytest.approx(1.0)
            assert row["inventory_kg"] >= 0 and row["buried_kg"] >= -1e-6
        inventories = [r["inventory_kg"] for r in box["states"]]
        if b * tau > 1e9:
            assert inventories == sorted(inventories) and inventories[-1] <= box["steady_kg"]
    m = sp.sea_matrix(1000.0)
    assert np.max(np.abs(m.sum(axis=0))) < 1e-15
    with pytest.raises(ValueError):
        sp.sea_matrix(0.0)
    with pytest.raises(ValueError):
        sp.sea_box(-1.0, 1000.0, 1e16)


def test_sea_tables_follow_the_closed_form():
    seas = PRODUCT["seas"]
    v = seas["volume_m3"]
    tau_lo, tau_hi = seas["residence_yr"]["box_model"], seas["residence_yr"]["earth_ocean_low"]

    def expect(b, tau, t):
        return b * tau * (1 - math.exp(-t / tau)) / v * 1e6
    for case, (lo, hi) in seas["input_cases_kg_yr"].items():
        table = seas["table_mg_m3"][case]
        assert table["1000_yr"] == pytest.approx([min(expect(lo, tau_lo, 1e3), expect(lo, tau_hi, 1e3)),
                                                  max(expect(hi, tau_lo, 1e3), expect(hi, tau_hi, 1e3))])
        assert table["settled_earth_ocean_low"] == pytest.approx([lo * tau_hi / v * 1e6, hi * tau_hi / v * 1e6])
    ocean = seas["earth_ocean"]["mean_mg_m3"]
    assert seas["table_mg_m3"]["dissolved_after_lakes"]["1000_yr"][1] < 0.02 * ocean
    # Every body's input is its share of the land's dissolved loss in the waters product.
    near = seas["bodies"]["near_side_sea"]
    share = WATERS["by_receiving_water"]["Oceanus Procellarum"]["dissolved_p_share"]
    assert near["input_cases_kg_yr"]["dissolved_after_lakes"] == pytest.approx(
        [share * x for x in PRODUCT["land_loss"]["dissolved_kg_yr"]])
    # Release from particles never lowers a body's input, and one sea taking everything holds the most.
    for key, body in seas["bodies"].items():
        cases = body["input_cases_kg_yr"]
        if "with_particle_release" in cases:
            assert all(r >= d for r, d in zip(cases["with_particle_release"], cases["dissolved_after_lakes"])), key
    takes_all = seas["bodies"]["near_side_sea_takes_all"]["table_mg_m3"]["dissolved_after_lakes"]["1000_yr"]
    assert takes_all[1] > near["table_mg_m3"]["dissolved_after_lakes"]["1000_yr"][1]
    atlas = sp.load("geography/results/atlas.json")
    named = sum(b["volume_km3"] for b in atlas["bodies"]) * 1e9
    assert seas["named_bodies_volume_share"] == pytest.approx(named / v)
    assert seas["named_bodies_factor_on_mixed"] == pytest.approx(v / named)
    lever = seas["delivered_lever"]["box_model"]["states"]
    assert lever[0]["mg_m3"] == pytest.approx(ocean)
    assert [r["mg_m3"] for r in lever] == sorted((r["mg_m3"] for r in lever), reverse=True)


def test_food_systems_match_their_closed_form():
    eaten = PRODUCT["people"]["eaten_kg_per_person"]
    loss_share, ret_x = sp.FOOD_LOSS_SHARE, sp.WASTEWATER_P_RECOVERY
    ret_w = sp.DESIGN_GUESSES["food_waste_p_return"]["range"]
    for name, system in PRODUCT["people"]["systems"].items():
        eff = system["efficiency"]
        for edge, e, w, eta, rx, rw in (("best", eaten[0], loss_share[0], eff[1], ret_x[1], ret_w[1]),
                                        ("worst", eaten[1], loss_share[1], eff[0], ret_x[0], ret_w[0])):
            new = e / ((1 - w) * eta) - rx * e - rw * e * w / (1 - w)
            row = system["per_person_kg"][edge]
            assert row["new"] == pytest.approx(new, rel=1e-12)
            assert row["new"] == pytest.approx(row["unrecovered_on_farms"] + row["lost_from_towns"], rel=1e-12)
        assert system["new_per_billion_kg"] == pytest.approx([x["new"] * 1e9 for x in system["per_person_kg"].values()])
    # Aquaponics returns 85-90% of the 55% of feed P that reaches the sludge.
    assert PRODUCT["people"]["systems"]["aquaponics"]["efficiency"] == pytest.approx([1 - 0.55 * 0.15, 1 - 0.55 * 0.10])
    systems = PRODUCT["people"]["systems"]
    assert systems["hydroponics"]["new_per_billion_kg"][1] < systems["field_crops"]["new_per_billion_kg"][1]
    # With nothing lost on farms, the design target's 0.044 kg needs towns to return 92.5%.
    w = loss_share[0] / (1 - loss_share[0])
    need = PRODUCT["people"]["ecology_design_target"]["town_return_needed_at_low_end"]
    assert (1 - need) * eaten[0] * (1 + w) == pytest.approx(0.044)
    assert sp.food_system(1.0, 0.0, 1.0, 1.0, 1.0)["new"] == 0.0
    for bad in ((1.0, 1.0, 0.5, 0.5, 0.5), (1.0, 0.2, 0.0, 0.5, 0.5), (1.0, 0.2, 0.5, 1.2, 0.5)):
        with pytest.raises(ValueError):
            sp.food_system(*bad)


def test_ledgers_close_at_every_corner():
    loss, part = PRODUCT["land_loss"], PRODUCT["partition"]
    carriers = PRODUCT["natural_cycle"]["decent"]["wild_dissolved_loss_kg_yr"]
    for case, people in (("design_target", PRODUCT["people"]["ecology_design_target"]["for_working_range_kg"]),
                         ("field_crops_all_unrecovered_to_waters",
                          PRODUCT["people"]["systems"]["field_crops"]["new_for_working_range_kg"])):
        shares = []
        for i in range(2):
            for j in range(2):
                into = loss["total_kg_yr"][i] + people[j]
                to_lakes = into * (1 - part["to_seas_share_of_loss"][i])
                assert to_lakes - people[j] - carriers[i] > 0
                shares.append(people[j] / to_lakes)
        assert PRODUCT["ledger"][case]["dredged_share_of_lake_input"] == pytest.approx([min(shares), max(shares)])
        for row in PRODUCT["ledger"][case]["corners"].values():
            assert abs(row["closure_kg_yr"]) <= 1e-9 * row["into_waters"]


# ---- Limits hold ----------------------------------------------------------------------------------------------------
def test_biological_routes_stay_under_their_waters_input():
    rt, part, seas = PRODUCT["routes"], PRODUCT["partition"], PRODUCT["seas"]
    ceiling = sp.DESIGN_GUESSES["extraction_ceiling_share"]["range"][0]
    inputs = dict(seas=seas["input_cases_kg_yr"]["dissolved_after_lakes"], lakes=part["kept_by_lakes_dissolved_kg_yr"])
    for water, b in inputs.items():
        for route in ("seaweed", "fish"):
            assert rt[route]["ceiling_kg_yr"][water]["base"] == pytest.approx([ceiling * x for x in b])
    assert rt["bioavailable_input_kg_yr"]["lakes"]["with_particle_release"][1] == pytest.approx(
        sp.RIVER_REACTIVE_SHARE[1] * part["kept_by_lakes_kg_yr"][1])
    need = rt["need_for_working_range_kg"]
    assert rt["waters_life_ceiling_kg_yr"]["base"][1] < need[0]
    assert rt["dredging"]["supply_kg_yr"]["lake_delta_particles"][0] > need[0]
    assert rt["dredging"]["supply_kg_yr"]["lake_delta_particles"] == part["kept_by_lakes_particulate_kg_yr"]
    for key, row in PRODUCT["natural_cycle"]["sharing"]["lakes"].items():
        assert row["left_to_food_webs_and_burial_kg_yr"] > 0, key
    # Carriers that return the whole dissolved loss eat more P than the lakes receive in solution.
    assert PRODUCT["natural_cycle"]["flyers"]["p_eaten_over_lakes_dissolved_input"][0] > 1


def test_sorption_energy_follows_concentration():
    s = PRODUCT["routes"]["sorption"]
    earth = s["earth_ocean_mean"]
    lunar = s["by_case"]["dissolved_after_lakes"]["at_1000_years"][0]
    ratio = earth["mg_m3"] / lunar["mg_m3"]
    assert lunar["water_m3_per_kg_p"][0] == pytest.approx(earth["water_m3_per_kg_p"][0] * ratio)
    assert lunar["pumping_mj_per_kg_p"][1] == pytest.approx(earth["pumping_mj_per_kg_p"][1] * ratio)
    assert lunar["minimum_work_mj_per_kg_p"] < lunar["pumping_mj_per_kg_p"][0] / 100
    assert lunar["minimum_work_mj_per_kg_p"] > earth["minimum_work_mj_per_kg_p"]
    with pytest.raises(ValueError):
        sp.minimum_work_mj_per_kg(0.0, 290.0, 1.0)


def test_dredging_masses_grades_and_fuel():
    d = PRODUCT["routes"]["dredging"]
    grade = PRODUCT["land_loss"]["soil_p_grade"]
    assert d["sediment_kg_per_kg_p"] == pytest.approx([1 / grade["kreep_rich"], 1 / grade["highland"]])
    lo, hi = sp.DESIGN_GUESSES["sediment_dry_bulk_t_m3"]["range"]
    assert d["volume_m3_per_kg_p"] == pytest.approx([d["sediment_kg_per_kg_p"][0] / (hi * 1e3),
                                                     d["sediment_kg_per_kg_p"][1] / (lo * 1e3)])
    # Hopper dredging, which lifts the sediment, sets the low end of the fuel.
    assert d["fuel_kg_per_m3"][0] == pytest.approx(0.5 * 0.820)
    g = MOON_GM / MOON_RADIUS**2
    assert d["lift_mj_per_kg_p_per_km"][0] == pytest.approx(d["sediment_kg_per_kg_p"][0] * g / 1e3)


def test_pair_and_fraction_guards():
    with pytest.raises(ValueError):
        sp.pair(2.0, 1.0)
    with pytest.raises(ValueError):
        sp.fraction(1.5)
