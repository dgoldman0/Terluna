"""The seas' phosphorus, people's phosphorus and the routes that carry it from the waters to farms and towns.

    python -m resources.sea_phosphorus           # results/sea_phosphorus.json
    python -m resources.sea_phosphorus_waters    # first, when the routing changes (ignored grid inputs)

Closed-form arithmetic and linear sea boxes on committed products and on products bound from the ecology and
provisioning domains (BOUND, read in-tree since the joint integration of 12 October 2026 and checked by SHA-256).
Each sea box is the
resources box model's sea and sediment boxes with its first-order burial, fed by the ecology product's land loss as
the routed lakes pass it on. Literature values are named where they are defined and recorded in
sea_phosphorus_sources.json; design guesses are collected in DESIGN_GUESSES.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np

from shared.constants import GAS_CONSTANT, JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_RADIUS, MOON_SURFACE_GRAVITY
from shared.provenance import constants_used
from resources.phosphorus import P_ATOMIC, evolve

ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent
OUT = HERE / "results/sea_phosphorus.json"
SOURCES = HERE / "sea_phosphorus_sources.json"
SCHEMA = "terluna.resources.sea-phosphorus/1"
YEAR_S = JULIAN_DAY * JULIAN_YEAR_DAYS
MOON_AREA_M2 = 4 * math.pi * MOON_RADIUS**2
CODE = ("resources/sea_phosphorus.py", "resources/phosphorus.py")

# Committed products read here: path -> schema.
INPUTS = {
    "resources/results/phosphorus.json": "terluna.resources.phosphorus/1",
    "resources/results/first_screen.json": "terluna.resources.first-screen/1",
    "resources/results/sea_phosphorus_waters.json": "terluna.resources.sea-phosphorus-waters/1",
    "geography/results/atlas.json": "terluna.geography.atlas/1",
    "geography/results/drainage.json": "terluna.geography.drainage/1",
    "research/studies/lunar_cycle_ecology/results/lunar_cycle_ecology.json": "terluna.research.lunar-cycle-ecology/1",
    "research/studies/conservation/results/gates.json": None,
    "resources/sea_phosphorus_sources.json": None,
}
# Inputs from the former domain branches, read in-tree since the joint integration, their hashes checked.
BOUND = {
    "ecology_people": dict(path="biosphere/ecology/results/people.json", schema="terluna.biosphere.ecology-people/1",
                           branch="main", commit=None,
                           sha256="758a8dcf68c6eb58082e5749e22c35e2da7ffd530f7483622e6c1b3011779ee8"),
    "sky_ecology": dict(path="research/studies/sky_ecology/results/sky_ecology.json",
                        schema="terluna.research.sky-ecology/1", branch="main", commit=None,
                        sha256="2d0f1f08d95cfb0f3c20f11b5097535c4d97b9cbb439fe76a5315378f819342a"),
    "sky_phosphorus_note": dict(path="research/studies/sky_ecology/phosphorus.md", schema=None,
                                branch="main", commit=None,
                                sha256="65ec921e7a5dea769a120406c54a22e5bfc47d3077fb06439408b0b9d05f84a2"),
    "living_water": dict(path="biosphere/living_water/waters.json", schema="terluna.biosphere.living-water/1",
                         branch="main", commit=None,
                         sha256="cda0594b8b59acfb469a39b007f977d616db60307d325e80af759242d2baecd1"),
    "ways_of_living": dict(path="research/studies/ways_of_living/results/ways_of_living.json",
                           schema="terluna.research.ways-of-living/1", branch="main", commit=None,
                           sha256="b55285cd89077ee8a2366e638b942e2ac0faa890d25fdf0c6f15ff8a3dffa65a"),
    "provisioning_decisions": dict(path="research/decisions.md", schema=None, branch="main",
                                   commit=None,
                                   sha256="b3255386680217a7ff78dc2b4731e42133589ba03390f6abf67c42bea6e2ea6a"),
}
# Bound at integration when provisioning's food-systems study is committed; until then the systems below stand alone.
PENDING = {
    "provisioning_food_systems": dict(
        branch="domain/provisioning", status="being built",
        replaces=["people.systems.*.efficiency", "people.food_loss_share", "people.waste_return",
                  "people.excreta_return", "the shares of each system in the diet"],
        rule="Bind the study's per-system P applied, recovered and lost per person; the ledger and routes "
             "then follow.")}

# Atomic weights (IUPAC standard values) for stoichiometric conversion, the convention of phosphorus.py.
O_ATOMIC, MG_ATOMIC, N_ATOMIC, H_ATOMIC = 15.999, 24.305, 14.007, 1.008
P_IN_P2O5 = 2 * P_ATOMIC / (2 * P_ATOMIC + 5 * O_ATOMIC)
STRUVITE_PER_P = (MG_ATOMIC + N_ATOMIC + P_ATOMIC + 16 * H_ATOMIC + 10 * O_ATOMIC) / P_ATOMIC  # MgNH4PO4.6H2O
MG_PER_P = MG_ATOMIC / P_ATOMIC

# Literature (sea_phosphorus_sources.json).
EARTH_OCEAN_P_MOL = 3e15              # Paytan & McLaughlin 2007: the ocean holds about 3e15 mol of P
EARTH_OCEAN_KM3 = 1.3324e9            # Charette & Smith 2010
EARTH_OCEAN_P_RESIDENCE_YR = (2e4, 1e5)   # Paytan & McLaughlin 2007: 20-100 kyr
RIVER_REACTIVE_SHARE = (0.1, 0.3)     # Compton et al. 2000, via Paytan & McLaughlin 2007: 10-30% of river P reactive
LUNAR_SOIL_P2O5_WT_PCT = dict(highland=0.11, mare=0.13, kreep_rich=0.51)  # Lunar Sourcebook Table 7.15 (register)
SEAWEED_CHINA = dict(area_km2=1250.0, harvest_t_dry=2004600.0, p_removed_t=9496.0,
                     p_dry_share=(0.0034, 0.0052))   # Xiao et al. 2017, Chinese seaweed farms in 2014
DREDGE_FUEL_KG_M3 = dict(trailing_suction_hopper_10_km=0.78, trailing_suction_hopper_27_km=1.15)  # CEDA 2022
DREDGE_FUEL_L_M3 = dict(trailing_suction_hopper=0.5, backhoe_barge=1.0)   # CEDA 2022, Table C3-1 campaign
DIESEL_MJ_KG = 43.0                   # IPCC 2006, net calorific value of gas/diesel oil, 43.0 TJ/Gg
DIESEL_KG_L = (0.820, 0.845)          # EN 590 density at 15 C
FISHING_FUEL_L_T = dict(small_pelagic=42.0, world_mean=489.0, demersal_fish=539.0,
                        crustaceans=1739.0)   # Parker et al. 2018, Table 1, 2011
EARTH_OCEAN_AIR_P_UMOL_M2_DAY = (0.08, 0.9)   # Paytan & McLaughlin 2007: Pacific low to Mediterranean
WORLD_LANDINGS_T = 81e6               # Parker et al. 2018, marine landings in 2011
FIELD_P_EFFICIENCY = (0.5, 0.9)       # Syers et al. 2008: balance-method recovery often 50-90%
RECIRCULATION_P_SAVING = (0.029, 0.176)  # Sanchez del Castillo et al. 2021, Table 3: P saved by reusing drainage
AQUAPONIC_FEED_P = dict(sludge=0.55, soluble=0.17, mineralised=(0.85, 0.90))  # Goddek et al. 2019 (Neto & Ostrensky)
WASTEWATER_P_RECOVERY = (0.6, 0.9)    # Egle et al. 2016: sludge-ash routes, share of the P the works receive
FOOD_LOSS_SHARE = (0.25, 0.33)        # FAO 2011, as the ecology people note takes it: a quarter to a third
INSECT_P_DRY = (0.0062, 0.013)        # Woods et al. 2004 (0.62-0.97%); Gratton et al. 2008 (1.3%), register H6
MIDGE_LANDING_M = (100.0, 200.0)      # Dreyer et al. 2015 (register H5): most midges land within 100-200 m
FISH_P_WET = 0.006                    # Otero et al. 2018's seabird diet, as the sky-ecology study takes it

# Design guesses, labelled once here and carried as ranges into the product.
DESIGN_GUESSES = dict(
    hydroponic_p_efficiency=dict(range=[0.8, 0.95], basis="closed recirculation with bleed and residue P returned; "
                                 "reusing the fifth of the solution that drained applied 3-18% less P than an open "
                                 "system"),
    food_waste_p_return=dict(range=[0.5, 0.9], basis="composted or digested food losses back on the fields"),
    sediment_dry_bulk_t_m3=dict(range=[0.5, 1.5], basis="fine lake muds to sandy deltas, in place"),
    sorption_capture=dict(range=[0.5, 0.9],
                          basis="share of the dissolved P a contactor takes from the water it treats"),
    sorption_pressure_drop_kpa=dict(range=[1.0, 10.0], basis="a low-resistance sorbent contactor"),
    sorption_concentrate_mol_l=dict(range=[1.0, 1.0], basis="the phosphate concentrate the minimum work runs to"),
    extraction_ceiling_share=dict(range=[0.3, 0.3], basis="largest share of a water body's bioavailable input taken "
                                  "by fish, seaweed or sorption; the rest stays with its food webs and carriers"),
)
YEARS = (0.0, 100.0, 500.0, 1000.0, 5000.0, 10000.0, 100000.0)
SEA_BODIES = dict(near_side_sea="Oceanus Procellarum", orientale="Mare Orientale", south_pole_aitken_sea="Mare Ingenii")
UNIT_C_TO_K = 273.15
EVIDENCE = ("A screen on committed and bound products: the ecology product's land loss split into its dissolved and "
            "particulate parts, passed through the rain-fed lakes by the waters product's rerun routing and Earth "
            "lake-retention relations, into linear sea boxes built on the resources box model; Earth figures from "
            "the literature for food systems, recovery, dredging, seaweed, fishing and the ocean. No sea chemistry, "
            "sediment diagenesis, lake ecology, food-web or carrier population is simulated, and no food system, "
            "route or population is chosen.")
READING_RULE = ("Phosphorus in kg of elemental P and flows per Julian year unless a key names Tg (1e9 kg), kt (1e6 kg) "
                "or g/m2. Pairs are [low, high] with every range taken at the end that bounds the quantity named; "
                "ledger corners pair each land-loss end with each population end. Concentrations are mg P per m3 of "
                "water (= ug/L) with uM beside them. Per-billion values are per 1e9 people on the Moon. 'screen' "
                "marks numbers this code computes, 'derived' arithmetic on cited values. Shares of the land's loss "
                "use the ecology product's 1.47-4.16 Tg a year. The waters product's lakes are connected bodies of "
                "the drainage product's fill-and-spill lakes; its outlets include every cell below the common water "
                "level, so 'all common-level water' bounds the named seas' input from above and their direct inflow "
                "from below. Bioavailable inputs: 'base' is dissolved P; 'with_particle_release' adds Earth's 10-30% "
                "reactive share of all river P, which salinity or anoxic beds release from particles.")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class MissingInput(RuntimeError):
    """A bound input is not in this tree."""


def bound_bytes(key):
    """The bound input's bytes: the in-tree file, which since the joint integration of 12 October 2026 holds every
    former branch's inputs; its hash is checked against the pin."""
    b = BOUND[key]
    local = ROOT / b["path"]
    if not local.exists():
        raise MissingInput(f"{b['path']} is not in this tree")
    raw = local.read_bytes()
    if hashlib.sha256(raw).hexdigest() != b["sha256"]:
        raise ValueError(f"{b['path']} does not match its pinned hash")
    return raw


def load_bound(key):
    raw = bound_bytes(key)
    if BOUND[key]["schema"] is None:
        return raw.decode()
    data = json.loads(raw)
    if data.get("schema") != BOUND[key]["schema"]:
        raise ValueError(f"Unexpected schema in bound input {key}")
    return data


def load(rel):
    data = json.loads((ROOT / rel).read_text())
    if INPUTS[rel] is not None and data.get("schema") != INPUTS[rel]:
        raise ValueError(f"Unexpected schema in {rel}")
    return data


def pair(lo, hi):
    if lo > hi:
        raise ValueError(f"Range runs backwards: {lo} > {hi}")
    return [float(lo), float(hi)]


def fraction(x):
    if not 0 <= x <= 1:
        raise ValueError("Fraction must be between zero and one")
    return x


# ---- 1. The seas' phosphorus --------------------------------------------------------------------------------------
def land_loss(people, sky, register, gates, waters):
    """The land's loss to the waters, its dissolved and particulate parts, and what it implies."""
    total = [x * 1e9 for x in people["capacity"]["land_p_loss_tg"]]
    dissolved = list(sky["phosphorus"]["sea_to_land"]["dissolved_loss_kg_year"])
    particulate = pair(total[0] - dissolved[0], total[1] - dissolved[1])
    dry_m2 = people["land"]["dry_land_km2"] * 1e6
    grade = {k: v / 100 * P_IN_P2O5 for k, v in LUNAR_SOIL_P2O5_WT_PCT.items()}
    land_grade = pair(grade["highland"], grade["mare"])          # the dry land's soils; KREEP-rich ground is local
    erosion = pair(particulate[0] / land_grade[1], particulate[1] / land_grade[0])
    stock = register["return_limb"]["rock_stock_top_metre"]
    top_metre = pair(stock["highland_soil"]["P_kg_m2"] * dry_m2, stock["mare_soil"]["P_kg_m2"] * dry_m2)
    runoff_basis = register["return_limb"]["runoff_m_per_yr_over_land"] * dry_m2
    land_runoff = waters["water_m3s"]["land_runoff"] * YEAR_S
    gate_rock = [v["rock_gt_per_year"] * 1e12 for v in gates["sediment_return"].values()]
    return dict(
        total_kg_yr=total, dissolved_kg_yr=dissolved, particulate_kg_yr=particulate,
        dissolved_share=pair(dissolved[1] / total[1], dissolved[0] / total[0]),
        per_m2_dry_land_g=pair(total[0] / dry_m2 * 1e3, total[1] / dry_m2 * 1e3),
        soil_p_grade=grade, land_soil_grade=land_grade,
        implied_erosion_kg_yr=erosion, conservation_erosion_kg_yr=pair(min(gate_rock), max(gate_rock)),
        top_metre_p_kg=top_metre,
        cycle_time_against_top_metre_yr=pair(top_metre[0] / total[1], top_metre[1] / total[0]),
        water_basis_m3_yr=dict(ecology_loss_basis=runoff_basis, land_runoff=land_runoff,
                               land_runoff_over_basis=land_runoff / runoff_basis,
                               loss_on_land_runoff_kg_yr=pair(total[0] * land_runoff / runoff_basis,
                                                              total[1] * land_runoff / runoff_basis)),
        marking="derived")


def partition(loss, waters, drainage, first):
    """Where the land's loss goes: kept by the rain-fed lakes, or passed to the common-level waters."""
    s_d = waters["dissolved_p_survival"]
    s_p = sorted(waters["particulate_p_survival"].values())
    d, p = loss["dissolved_kg_yr"], loss["particulate_kg_yr"]
    seas_d = pair(s_d * d[0], s_d * d[1])
    seas_p = pair(s_p[0] * p[0], s_p[-1] * p[1])
    seas = pair(seas_d[0] + seas_p[0], seas_d[1] + seas_p[1])
    lakes_d = pair(d[0] - seas_d[0], d[1] - seas_d[1])
    lakes_p = pair(p[0] - s_p[-1] * p[0], p[1] - s_p[0] * p[1])
    lakes = pair(lakes_d[0] + lakes_p[0], lakes_d[1] + lakes_p[1])
    to_level = sum(row["dissolved_p_share"] for row in waters["by_receiving_water"].values())
    named = {}
    for name, row in waters["by_receiving_water"].items():
        if name in ("small waters", "below the atlas grid") or row.get("discharge_m3s", 0) > 100:
            named[name] = dict(discharge_share=row["discharge_m3s"] / waters["water_m3s"]["to_common_level_water"],
                               dissolved_share_of_arrivals=row["dissolved_p_share"] / to_level,
                               dissolved_kg_yr=pair(row["dissolved_p_share"] * d[0], row["dissolved_p_share"] * d[1]),
                               particulate_kg_yr=pair(row["particulate_p_share_median_moon"] * p[0],
                                                      row["particulate_p_share_median_moon"] * p[1]))
    lake_area = MOON_AREA_M2 * drainage["lakes"]["area_share"]
    return dict(
        lakes=waters["lakes"],
        outflow_weighted_tau_yr=waters["outflow_weighted_tau_yr"],
        lake_passage_share_of_discharge=waters["lake_passage_share_of_discharge"],
        lake_passage_share_of_land_runoff=waters["lake_passage_share_of_land_runoff"],
        dissolved_survival=s_d, particulate_survival=pair(s_p[0], s_p[-1]),
        to_seas_dissolved_kg_yr=seas_d, to_seas_particulate_kg_yr=seas_p, to_seas_kg_yr=seas,
        to_seas_share_of_loss=pair(seas[0] / loss["total_kg_yr"][0], seas[1] / loss["total_kg_yr"][1]),
        kept_by_lakes_dissolved_kg_yr=lakes_d, kept_by_lakes_particulate_kg_yr=lakes_p, kept_by_lakes_kg_yr=lakes,
        kept_by_lakes_share_of_loss=pair(lakes[0] / loss["total_kg_yr"][0], lakes[1] / loss["total_kg_yr"][1]),
        kept_by_lakes_g_m2_lake=pair(lakes[0] / lake_area * 1e3, lakes[1] / lake_area * 1e3),
        lake_area_m2=lake_area, lake_volume_m3=first["inventories"]["water"]["lakes_kg"] / 1000,
        by_receiving_water=named,
        largest_rivers_receiving_water=waters["largest_rivers_receiving_water"],
        if_lakes_passed_all_dissolved_kg_yr=list(d),
        marking="screen")


def sea_matrix(tau_yr):
    """The box model's sea and sediment boxes with first-order burial, and a unit source as a third state."""
    if tau_yr <= 0:
        raise ValueError("Positive residence required")
    m = np.zeros((3, 3))
    m[0, 0], m[1, 0] = -1 / tau_yr, 1 / tau_yr       # sea water -> sediment
    return m


def sea_box(input_kg_yr, tau_yr, volume_m3, years=YEARS, initial_kg=0.0):
    """Inventory, concentration and burial of bioavailable P in one water, by matrix exponential and closed form."""
    if input_kg_yr < 0 or initial_kg < 0 or volume_m3 <= 0:
        raise ValueError("Invalid sea box input")
    m = sea_matrix(tau_yr)
    m[0, 2] = input_kg_yr                              # constant input carried by the unit state
    rows = []
    for t in years:
        sea, sediment, unit = evolve(np.array([initial_kg, 0.0, 1.0]), m, t)
        closed = input_kg_yr * tau_yr * (1 - math.exp(-t / tau_yr)) + initial_kg * math.exp(-t / tau_yr)
        mg_m3 = sea / volume_m3 * 1e6
        rows.append(dict(year=t, inventory_kg=float(sea), buried_kg=float(sediment), unit_state=float(unit),
                         closed_form_kg=closed, mg_m3=float(mg_m3), um=float(mg_m3 / P_ATOMIC)))
    steady = input_kg_yr * tau_yr
    return dict(input_kg_yr=input_kg_yr, tau_yr=tau_yr, states=rows, steady_kg=steady,
                steady_mg_m3=steady / volume_m3 * 1e6, steady_um=steady / volume_m3 * 1e6 / P_ATOMIC)


def table(runs, taus_shown=("box_model", "earth_ocean_low")):
    """mg/m3 by year over the residence times shown, and settled at each, as [low, high]."""
    rows = {}
    for year in (100.0, 1000.0, 10000.0):
        values = [row["mg_m3"] for t in taus_shown for b in runs[t] for row in b["states"] if row["year"] == year]
        rows[f"{year:g}_yr"] = pair(min(values), max(values))
    for t in runs:
        rows[f"settled_{t}"] = pair(runs[t][0]["steady_mg_m3"], runs[t][1]["steady_mg_m3"])
    return rows


def seas(part, loss, people, first, box_model, ways, atlas, waters):
    """The common-level waters' inventory and concentration through time, together and by body."""
    volume = first["inventories"]["water"]["seas_kg"] / 1000
    rate = next(e["rate_per_year"] for e in box_model["finite_reservoir_cases"][0]["edges"]
                if e["source"] == "sea" and e["destination"] == "sediment")
    taus = dict(box_model=1 / rate, earth_ocean_low=EARTH_OCEAN_P_RESIDENCE_YR[0],
                earth_ocean_high=EARTH_OCEAN_P_RESIDENCE_YR[1])
    sea_area = MOON_AREA_M2 * people["land"]["seas"]
    dust = [x / 1e3 * sea_area for x in people["aloft"]["dust_p_g_m2_yr"]]
    total_in = part["to_seas_kg_yr"]
    d_in = part["to_seas_dissolved_kg_yr"]
    ocean_air = [x * 1e-6 * P_ATOMIC * JULIAN_YEAR_DAYS for x in EARTH_OCEAN_AIR_P_UMOL_M2_DAY]   # g/m2 a year
    cases = dict(
        dissolved_after_lakes=d_in,
        dissolved_after_lakes_and_dust=pair(d_in[0] + dust[0], d_in[1] + dust[0]),
        dissolved_after_lakes_and_dust_top=pair(d_in[0] + dust[1], d_in[1] + dust[1]),
        with_particle_release=pair(max(d_in[0], RIVER_REACTIVE_SHARE[0] * total_in[0]),
                                   max(d_in[1], RIVER_REACTIVE_SHARE[1] * total_in[1])),
        lakes_pass_dissolved=list(loss["dissolved_kg_yr"]))
    runs = {c: {t: [sea_box(b, tau, volume) for b in rng] for t, tau in taus.items()} for c, rng in cases.items()}
    # Each body on the inflow the routing gives it directly; the pits below the atlas grid hold no atlas volume.
    vol = {b["water_names"][0]: b["volume_km3"] * 1e9 for b in atlas["bodies"] if b["water_names"]}
    named_volume = sum(b["volume_km3"] * 1e9 for b in atlas["bodies"])
    bodies = {}
    rows = waters["by_receiving_water"]
    d, p = loss["dissolved_kg_yr"], loss["particulate_kg_yr"]
    for key, name in list(SEA_BODIES.items()) + [("crater_waters", "small waters")]:
        v = vol[name] if name in vol else volume - named_volume
        share_d = rows.get(name, {}).get("dissolved_p_share", 0.0)
        share_p = rows.get(name, {}).get("particulate_p_share_median_moon", 0.0)
        body_d = pair(share_d * d[0], share_d * d[1])
        body_total = [body_d[0] + share_p * p[0], body_d[1] + share_p * p[1]]
        body_cases = dict(dissolved_after_lakes=body_d,
                          with_particle_release=pair(max(body_d[0], RIVER_REACTIVE_SHARE[0] * body_total[0]),
                                                     max(body_d[1], RIVER_REACTIVE_SHARE[1] * body_total[1])))
        discharge = rows.get(name, {}).get("discharge_m3s", 0.0) * YEAR_S
        body_runs = {c: {t: [sea_box(b, tau, v) for b in rng] for t, tau in taus.items()}
                     for c, rng in body_cases.items()}
        bodies[key] = dict(atlas_name=name, volume_m3=v, input_cases_kg_yr=body_cases,
                           water_residence_yr=(v / discharge) if discharge > 0 else None,
                           table_mg_m3={c: table(r) for c, r in body_runs.items()})
    # One sea taking everything that reaches the common level in solution, the Near-side Sea for scale.
    near = vol[SEA_BODIES["near_side_sea"]]
    takes_all = {t: [sea_box(b, tau, near) for b in d_in] for t, tau in taus.items()}
    bodies["near_side_sea_takes_all"] = dict(atlas_name=SEA_BODIES["near_side_sea"], volume_m3=near,
                                             input_cases_kg_yr=dict(dissolved_after_lakes=d_in),
                                             table_mg_m3=dict(dissolved_after_lakes=table(takes_all)))
    ocean_mol_m3 = EARTH_OCEAN_P_MOL / (EARTH_OCEAN_KM3 * 1e9)
    ocean_mg_m3 = ocean_mol_m3 * P_ATOMIC * 1e3
    lever = {t: sea_box(0.0, tau, volume, initial_kg=ocean_mg_m3 / 1e6 * volume) for t, tau in taus.items()}
    discharge = ways["moon"]["drainage_discharge_to_sea_km3_year"] * 1e9
    tables = {c: table(r) for c, r in runs.items()}
    return dict(
        volume_m3=volume, area_m2=sea_area, residence_yr=taus, input_cases_kg_yr=cases,
        named_bodies_volume_share=named_volume / volume, named_bodies_factor_on_mixed=volume / named_volume,
        dust_kg_yr_inherited_envelope=dust,
        earth_ocean_air_deposition=dict(g_m2_yr=ocean_air, kg_yr_over_these_waters=[x / 1e3 * sea_area
                                                                                    for x in ocean_air]),
        table_mg_m3=tables,
        at_1000_years_mg_m3={c: tables[c]["1000_yr"] for c in cases},
        bodies=bodies, runs=runs, delivered_lever=lever,
        earth_ocean=dict(mean_mol_m3=ocean_mol_m3, mean_mg_m3=ocean_mg_m3, mean_um=ocean_mol_m3 * 1e3),
        delivered_at_earth_ocean_mean_kg=ocean_mg_m3 / 1e6 * volume,
        water_residence_yr=volume / discharge, discharge_m3_yr=discharge,
        marking="screen")


def gathering(part, loss, waters, drainage, living):
    """Where the phosphorus gathers: lake deltas, the seas' river mouths and shelves, and the deeps."""
    rivers = drainage["rivers"]
    top = sum(rv["discharge_m3s"] for rv in rivers["largest"]) / rivers["discharge_to_sea_m3s"]
    hyp = waters["hypsometry"]

    def shelf(name):
        rows = hyp[name]["by_depth_m_km2"]
        return dict(area_km2=hyp[name]["area_km2"], under_50_m=rows["0-50"] / hyp[name]["area_km2"],
                    under_200_m=(rows["0-50"] + rows["50-200"]) / hyp[name]["area_km2"],
                    deeper_than_1000_m=(rows["1000-3000"] + rows["3000-inf"]) / hyp[name]["area_km2"])
    land_grade = loss["land_soil_grade"]
    lakes = part["kept_by_lakes_kg_yr"]
    particles = part["kept_by_lakes_particulate_kg_yr"]
    bulk = DESIGN_GUESSES["sediment_dry_bulk_t_m3"]["range"]
    deposit = pair(particles[0] / land_grade[1] / (bulk[1] * 1e3), particles[1] / land_grade[0] / (bulk[0] * 1e3))
    return dict(
        lakes=dict(kept_kg_yr=lakes, particles_kg_yr=particles, g_m2_lake=part["kept_by_lakes_g_m2_lake"],
                   sediment_m3_yr=deposit,
                   fill_time_yr=pair(part["lake_volume_m3"] / deposit[1], part["lake_volume_m3"] / deposit[0])),
        sea_river_mouths=dict(particulate_kg_yr=part["to_seas_particulate_kg_yr"], largest_twelve_share=top),
        shelves=dict(all_common_level_water=shelf("all_common_level_water"),
                     near_side_sea=shelf("Oceanus Procellarum"), south_pole_aitken_sea=shelf("Mare Ingenii")),
        design_waters_chlorophyll_mg_m3={k: living["waters"][k]["chlorophyll_mg_m3"]["range"]
                                         for k in ("open_sea", "productive_coast", "river_mouth")},
        marking="screen")


# ---- 2. People's phosphorus ---------------------------------------------------------------------------------------
def food_system(eaten, food_loss, efficiency, excreta_return, waste_return):
    """One person's P a year through a food system at steady state; the losses are what must come in new.

    For field crops, P the crops do not recover partly builds soil stores while the soils fill (Syers et al. 2008)
    and partly reaches the waters; both must be supplied."""
    for x in (food_loss, efficiency, excreta_return, waste_return):
        fraction(x)
    if eaten < 0 or efficiency == 0 or food_loss == 1:
        raise ValueError("Invalid food system")
    harvest = eaten / (1 - food_loss)
    applied = harvest / efficiency
    waste = harvest - eaten
    returned = excreta_return * eaten + waste_return * waste
    unrecovered_on_farms = applied - harvest
    lost_from_towns = (1 - excreta_return) * eaten + (1 - waste_return) * waste
    return dict(eaten=eaten, harvested=harvest, applied=applied, food_waste=waste, returned=returned,
                unrecovered_on_farms=unrecovered_on_farms, lost_from_towns=lost_from_towns, new=applied - returned)


def people_p(people, ways, decisions):
    """P per person in the lunar design's diet, by food system, and the new P per billion and for 4-8 billion."""
    design = people["needs"]["practices"]["lunar_design"]
    eaten = design["p_eaten_kg"]
    if "4–8 billion people on the Moon itself" not in decisions:
        raise ValueError("The working range is not in the bound decisions")
    population = (4e9, 8e9)
    feed = AQUAPONIC_FEED_P
    efficiency = dict(
        field_crops=list(FIELD_P_EFFICIENCY),
        hydroponics=DESIGN_GUESSES["hydroponic_p_efficiency"]["range"],
        aquaponics=pair(1 - feed["sludge"] * (1 - feed["mineralised"][0]),
                        1 - feed["sludge"] * (1 - feed["mineralised"][1])))
    waste_return = DESIGN_GUESSES["food_waste_p_return"]["range"]
    systems = {}
    natural = [x * 1e9 for x in people["capacity"]["land_p_loss_tg"]]
    for name, eff in efficiency.items():
        best = food_system(eaten[0], FOOD_LOSS_SHARE[0], eff[1], WASTEWATER_P_RECOVERY[1], waste_return[1])
        worst = food_system(eaten[1], FOOD_LOSS_SHARE[1], eff[0], WASTEWATER_P_RECOVERY[0], waste_return[0])
        per_billion = pair(best["new"] * 1e9, worst["new"] * 1e9)
        systems[name] = dict(efficiency=eff, per_person_kg=dict(best=best, worst=worst),
                             new_per_billion_kg=per_billion,
                             new_share_of_land_loss_per_billion=pair(per_billion[0] / natural[1],
                                                                     per_billion[1] / natural[0]),
                             new_for_working_range_kg=pair(per_billion[0] * population[0] / 1e9,
                                                           per_billion[1] * population[1] / 1e9))
    target = pair(design["p_lost_kg"][0] * 1e9, design["p_lost_kg"][1] * 1e9)
    # The design target's low end, 10% of what is eaten, met by towns alone with nothing lost on farms.
    w = FOOD_LOSS_SHARE[0] / (1 - FOOD_LOSS_SHARE[0])
    town_return_for_target = 1 - design["p_lost_kg"][0] / eaten[0] / (1 + w)
    return dict(
        diet=design["diet"], eaten_kg_per_person=eaten, food_loss_share=list(FOOD_LOSS_SHARE),
        excreta_return=list(WASTEWATER_P_RECOVERY), waste_return=waste_return, systems=systems,
        aquaponic_feed_p=dict(fish_and_other=1 - feed["sludge"] - feed["soluble"], soluble=feed["soluble"],
                              sludge=feed["sludge"], mineralised=list(feed["mineralised"]),
                              efficiency_sludge_discarded=1 - feed["sludge"]),
        ecology_design_target=dict(lost_kg_per_person=design["p_lost_kg"], per_billion_kg=target,
                                   share_of_land_loss_per_billion=pair(target[0] / natural[1], target[1] / natural[0]),
                                   for_working_range_kg=pair(target[0] * 4, target[1] * 8),
                                   town_return_needed_at_low_end=town_return_for_target),
        earth_practice_lost_kg_per_person=people["needs"]["practices"]["earth_practice"]["p_lost_kg"],
        excreta_per_billion_kg=pair(eaten[0] * 1e9, eaten[1] * 1e9),
        working_range_people=list(population), cropland_km2_per_billion=ways["per_billion_common"][
            "food_cropland_equivalent_km2"],
        recirculation_p_saving=list(RECIRCULATION_P_SAVING),
        pending_binding=PENDING["provisioning_food_systems"],
        marking="screen")


# ---- 3. Routes from the waters to people --------------------------------------------------------------------------
def minimum_work_mj_per_kg(mg_m3, kelvin, concentrate_mol_l):
    """Ideal-dilute work to raise dissolved P from mg_m3 to a concentrate, per kg of P (RT ln c_out/c)."""
    c = mg_m3 / P_ATOMIC / 1e3                         # mol/m3
    if c <= 0:
        raise ValueError("Positive concentration required")
    return GAS_CONSTANT * kelvin * math.log(concentrate_mol_l * 1e3 / c) / (P_ATOMIC / 1e3) / 1e6


def bioavailable(part, sea):
    """Dissolved input of each water, and the same with Earth's reactive share released from particles."""
    seas = sea["input_cases_kg_yr"]
    return dict(
        seas=dict(base=seas["dissolved_after_lakes"], with_particle_release=seas["with_particle_release"]),
        lakes=dict(base=part["kept_by_lakes_dissolved_kg_yr"],
                   with_particle_release=pair(max(part["kept_by_lakes_dissolved_kg_yr"][0],
                                                  RIVER_REACTIVE_SHARE[0] * part["kept_by_lakes_kg_yr"][0]),
                                              RIVER_REACTIVE_SHARE[1] * part["kept_by_lakes_kg_yr"][1])))


def routes(part, sea, loss, people_out, people, ways):
    need = people_out["ecology_design_target"]["for_working_range_kg"]
    ceiling = DESIGN_GUESSES["extraction_ceiling_share"]["range"][0]
    grade = loss["soil_p_grade"]
    bulk = DESIGN_GUESSES["sediment_dry_bulk_t_m3"]["range"]
    # A delta's sediment carries its catchment's soil, highland to KREEP-rich.
    sediment = pair(1 / grade["kreep_rich"], 1 / grade["highland"])           # kg sediment per kg P
    volume = pair(sediment[0] / (bulk[1] * 1e3), sediment[1] / (bulk[0] * 1e3))  # m3 per kg P
    fuel_m3 = pair(DREDGE_FUEL_L_M3["trailing_suction_hopper"] * DIESEL_KG_L[0],
                   DREDGE_FUEL_KG_M3["trailing_suction_hopper_27_km"])
    dredge_mj = pair(volume[0] * fuel_m3[0] * DIESEL_MJ_KG, volume[1] * fuel_m3[1] * DIESEL_MJ_KG)
    lift_km = pair(sediment[0] * MOON_SURFACE_GRAVITY * 1e3 / 1e6, sediment[1] * MOON_SURFACE_GRAVITY * 1e3 / 1e6)
    b = bioavailable(part, sea)
    ceil = {w: {c: pair(ceiling * v[0], ceiling * v[1]) for c, v in cases.items()} for w, cases in b.items()}
    deltas = part["kept_by_lakes_particulate_kg_yr"]
    dredging = dict(
        supply_kg_yr=dict(lake_delta_particles=deltas, lake_floors_from_solution=part["kept_by_lakes_dissolved_kg_yr"],
                          sea_mouth_particles=part["to_seas_particulate_kg_yr"]),
        sediment_kg_per_kg_p=sediment, volume_m3_per_kg_p=volume, fuel_kg_per_m3=fuel_m3,
        fuel_mj_per_kg_p=dredge_mj, lift_mj_per_kg_p_per_km=lift_km,
        for_working_range=dict(p_kg_yr=need, sediment_kg_yr=pair(need[0] * sediment[0], need[1] * sediment[1]),
                               volume_m3_yr=pair(need[0] * volume[0], need[1] * volume[1]),
                               fuel_power_w=pair(need[0] * dredge_mj[0] * 1e6 / YEAR_S,
                                                 need[1] * dredge_mj[1] * 1e6 / YEAR_S),
                               need_over_lake_delta_particles=pair(need[0] / deltas[1], need[1] / deltas[0])),
        marking="derived")
    p_per_m2 = SEAWEED_CHINA["p_removed_t"] * 1e3 / (SEAWEED_CHINA["area_km2"] * 1e6)
    dry_per_m2 = SEAWEED_CHINA["harvest_t_dry"] * 1e3 / (SEAWEED_CHINA["area_km2"] * 1e6)
    dry_per_p = pair(1 / SEAWEED_CHINA["p_dry_share"][1], 1 / SEAWEED_CHINA["p_dry_share"][0])
    seaweed = dict(
        uptake_kg_p_m2_yr=p_per_m2, harvest_kg_dry_m2_yr=dry_per_m2, dry_kg_per_kg_p=dry_per_p,
        ceiling_kg_yr=ceil,
        minimum_area_m2_at_ceiling={w: {c: [x / p_per_m2 for x in v] for c, v in cases.items()}
                                    for w, cases in ceil.items()},
        minimum_area_m2_for_working_range=pair(need[0] / p_per_m2, need[1] / p_per_m2),
        marking="derived")
    fish_per_p = 1 / FISH_P_WET
    fuel_l = pair(fish_per_p / 1e3 * FISHING_FUEL_L_T["small_pelagic"],
                  fish_per_p / 1e3 * FISHING_FUEL_L_T["demersal_fish"])
    fish = dict(
        fish_kg_per_kg_p=fish_per_p, fuel_l_per_kg_p=fuel_l, fishing_fuel_l_per_t=FISHING_FUEL_L_T,
        fuel_mj_per_kg_p=pair(fuel_l[0] * DIESEL_KG_L[0] * DIESEL_MJ_KG, fuel_l[1] * DIESEL_KG_L[1] * DIESEL_MJ_KG),
        ceiling_kg_yr=ceil,
        catch_kg_yr_at_ceiling={w: {c: [x * fish_per_p for x in v] for c, v in cases.items()}
                                for w, cases in ceil.items()},
        catch_kg_yr_for_working_range=pair(need[0] * fish_per_p, need[1] * fish_per_p),
        world_landings_kg_yr=WORLD_LANDINGS_T * 1e3,
        marking="derived")
    kelvin = people["sky"]["bands"]["low_sky"]["temperature_c"][0] + UNIT_C_TO_K
    capture = DESIGN_GUESSES["sorption_capture"]["range"]
    dp = DESIGN_GUESSES["sorption_pressure_drop_kpa"]["range"]
    pump = ways["assumptions"]["inherited"]["pump_efficiency"]
    out_mol_l = DESIGN_GUESSES["sorption_concentrate_mol_l"]["range"][0]

    def sorb(mg_m3):
        water = pair(1 / (mg_m3 / 1e6 * capture[1]), 1 / (mg_m3 / 1e6 * capture[0]))   # m3 per kg P
        return dict(mg_m3=mg_m3, water_m3_per_kg_p=water,
                    pumping_mj_per_kg_p=pair(water[0] * dp[0] * 1e3 / pump / 1e6, water[1] * dp[1] * 1e3 / pump / 1e6),
                    minimum_work_mj_per_kg_p=minimum_work_mj_per_kg(mg_m3, kelvin, out_mol_l),
                    km3_per_tg=1e9 / (mg_m3 / 1e6) / 1e9)
    at1000 = sea["at_1000_years_mg_m3"]
    sorption = dict(
        by_case={c: dict(at_1000_years=[sorb(x) for x in at1000[c]],
                         ceiling_kg_yr=pair(ceiling * sea["input_cases_kg_yr"][c][0],
                                            ceiling * sea["input_cases_kg_yr"][c][1]))
                 for c in ("dissolved_after_lakes", "lakes_pass_dissolved")},
        earth_ocean_mean=sorb(sea["earth_ocean"]["mean_mg_m3"]), kelvin=kelvin, pump_efficiency=pump,
        marking="derived")
    excreta = people_out["excreta_per_billion_kg"]
    wastewater = dict(
        recoverable_per_billion_kg=pair(excreta[0] * WASTEWATER_P_RECOVERY[0], excreta[1] * WASTEWATER_P_RECOVERY[1]),
        struvite_kg_per_kg_p=STRUVITE_PER_P, magnesium_kg_per_kg_p=MG_PER_P,
        magnesium_per_billion_kg=pair(excreta[0] * WASTEWATER_P_RECOVERY[0] * MG_PER_P,
                                      excreta[1] * WASTEWATER_P_RECOVERY[1] * MG_PER_P),
        marking="derived")
    return dict(need_for_working_range_kg=need, ceiling_share=ceiling, bioavailable_input_kg_yr=b,
                waters_life_ceiling_kg_yr={c: pair(ceil["seas"][c][0] + ceil["lakes"][c][0],
                                                   ceil["seas"][c][1] + ceil["lakes"][c][1])
                                           for c in ("base", "with_particle_release")},
                dredging=dredging, seaweed=seaweed, fish=fish, sorption=sorption, wastewater=wastewater)


def distribution(rt, ways, people_out):
    """What each route moves per kg of P to farms and towns, and how far the land lies from the waters."""
    d = ways["moon"]["land_distances_lunar_cycle_ecology"]
    crop = [x * 1e6 / 1e9 for x in people_out["cropland_km2_per_billion"]]          # m2 of cropland a person
    field = people_out["systems"]["field_crops"]["new_per_billion_kg"]
    new_m2 = pair(field[0] / 1e9 / crop[1], field[1] / 1e9 / crop[0])                # kg P per m2 of cropland a year
    sediment = rt["dredging"]["sediment_kg_per_kg_p"]
    bulk = DESIGN_GUESSES["sediment_dry_bulk_t_m3"]["range"]
    spread = dict(cropland_m2_per_person=crop, new_p_g_m2=[x * 1e3 for x in new_m2],
                  sediment_kg_m2=pair(new_m2[0] * sediment[0], new_m2[1] * sediment[1]),
                  sediment_mm=pair(new_m2[0] * sediment[0] / (bulk[1] * 1e3) * 1e3,
                                   new_m2[1] * sediment[1] / (bulk[0] * 1e3) * 1e3))
    carriers = dict(sediment=rt["dredging"]["sediment_kg_per_kg_p"], seaweed_dry=rt["seaweed"]["dry_kg_per_kg_p"],
                    fish=[rt["fish"]["fish_kg_per_kg_p"]] * 2, struvite=[STRUVITE_PER_P] * 2,
                    p2o5=[1 / P_IN_P2O5] * 2)
    lift = {k: [v[0] * MOON_SURFACE_GRAVITY * 1e3 / 1e6, v[1] * MOON_SURFACE_GRAVITY * 1e3 / 1e6]
            for k, v in carriers.items()}
    reach = dict(within_50_km=d["any_water"]["land_share_within"]["50_km"],
                 within_100_km=d["any_water"]["land_share_within"]["100_km"],
                 within_500_km=d["any_water"]["land_share_within"]["500_km"],
                 median_km_any_water=d["any_water"]["p50"], median_km_sea=d["sea"]["p50"])
    hauls = {k: [v[0] * reach["median_km_any_water"] / 1e3, v[1] * reach["median_km_any_water"] / 1e3]
             for k, v in carriers.items()}
    lakes = ways["moon"]["rain_fed_lakes"]
    return dict(carrier_kg_per_kg_p=carriers, lift_mj_per_kg_p_per_km=lift, land_reach=reach,
                tonne_km_per_kg_p_at_median_reach=hauls, sediment_on_field_crops=spread,
                lake_shore_km=lakes["big_lake_shore_km"] + lakes["small_lake_shore_km"],
                lakes=lakes["lakes"], marking="derived")


# ---- 4. The natural cycle -----------------------------------------------------------------------------------------
def natural(loss, part, rt, people, sky, register, ways):
    s = sky["phosphorus"]["sea_to_land"]
    total = loss["total_kg_yr"]
    dissolved = loss["dissolved_kg_yr"]
    targets = {t["target"]: t for t in s["targets"] if t["flyer"].startswith("wandering albatross")}
    lakes_in = part["kept_by_lakes_dissolved_kg_yr"]
    eaten = pair(targets["dissolved_loss"]["food_kg_wet_year"][0] * FISH_P_WET,
                 targets["dissolved_loss"]["food_kg_wet_year"][1] * FISH_P_WET)
    flyers = dict(
        earth_seabirds_kg_yr=s["earth_seabirds_kg_year"],
        earth_seabirds_share_of_loss=pair(s["earth_seabirds_kg_year"] / total[1],
                                          s["earth_seabirds_kg_year"] / total[0]),
        albatross_scale_for_dissolved_loss=targets["dissolved_loss"]["flyers"],
        p_eaten_for_dissolved_loss_kg_yr=eaten,
        p_eaten_over_lakes_dissolved_input=pair(eaten[0] / lakes_in[0], eaten[1] / lakes_in[1]),
        p_eaten_for_earth_seabird_scale_kg_yr=pair(targets["earth_seabird_scale"]["food_kg_wet_year"][0] * FISH_P_WET,
                                                   targets["earth_seabird_scale"]["food_kg_wet_year"][1] * FISH_P_WET))
    conveyor = register["return_limb"]["conveyor"]
    emergence = conveyor["myvatn_emergence_g_m2_yr"]
    lake_area = part["lake_area_m2"]
    insects_kg = pair(emergence[0] / 1e3 * lake_area * INSECT_P_DRY[0],
                      emergence[1] / 1e3 * lake_area * INSECT_P_DRY[1])
    lakes_shore = ways["moon"]["rain_fed_lakes"]
    shore_m = (lakes_shore["big_lake_shore_km"] + lakes_shore["small_lake_shore_km"]) * 1e3
    strip = pair(shore_m * MIDGE_LANDING_M[0], shore_m * MIDGE_LANDING_M[1])
    top_rate = [emergence[1] / 1e3 * lake_area * x for x in INSECT_P_DRY]
    insects = dict(lake_emergence_g_dry_m2=emergence, p_kg_yr=insects_kg, shore_strip_m2=strip,
                   g_p_m2_on_strip=pair(insects_kg[0] / strip[1] * 1e3, insects_kg[1] / strip[0] * 1e3),
                   myvatn_top_export_over_lakes_dissolved_input=pair(top_rate[0] / lakes_in[1],
                                                                     top_rate[1] / lakes_in[0]))
    run_fish = dict(fish_kg_per_kg_p=1 / FISH_P_WET,
                    fish_kg_yr_for_tenth_of_dissolved=pair(0.1 * dissolved[0] / FISH_P_WET,
                                                           0.1 * dissolved[1] / FISH_P_WET))
    falls = sky["phosphorus"]["falls"]
    falls_kg = pair(sum(f["p_t_year"][0] for f in falls) * 1e3, sum(f["p_t_year"][1] for f in falls) * 1e3)
    litter = [x * 1e3 for x in sky["phosphorus"]["deposition_flux_t_year"]["earth_forest"]]
    kept = people["capacity"]["kept_shares"]
    wild_dissolved = pair(min(kept.values()) * dissolved[0], max(kept.values()) * dissolved[1])
    flyers_needed = pair(min(kept.values()) * targets["dissolved_loss"]["flyers"][0],
                         max(kept.values()) * targets["dissolved_loss"]["flyers"][1])
    ceiling = rt["ceiling_share"]
    b = rt["bioavailable_input_kg_yr"]
    budget = dict(lakes={}, seas={c: dict(bioavailable_kg_yr=v, harvest_ceiling_kg_yr=pair(ceiling * v[0],
                                                                                         ceiling * v[1]))
                                  for c, v in b["seas"].items()})
    for i, edge in enumerate(("low", "high")):
        for case, rng in b["lakes"].items():
            for name, share in kept.items():
                carriers = share * dissolved[i]
                left = rng[i] - carriers - ceiling * rng[i]
                budget["lakes"][f"{edge}_{case}_{name}"] = dict(
                    bioavailable_kg_yr=rng[i], carriers_kg_yr=carriers, harvest_ceiling_kg_yr=ceiling * rng[i],
                    left_to_food_webs_and_burial_kg_yr=left, left_share=left / rng[i])
    return dict(
        flyers=flyers, lake_insects=insects, run_fish=run_fish,
        aerophyte_falls_kg_yr_at_1pct=falls_kg, litter_kg_yr_forest_floor_pools=litter,
        land_share_of_routes=sky["phosphorus"]["land_share_used"]["value"],
        return_ceiling_kg_yr=list(dissolved),
        return_ceiling_share_of_loss=pair(dissolved[1] / total[1], dissolved[0] / total[0]),
        return_ceiling_with_particle_release_kg_yr=pair(
            b["lakes"]["with_particle_release"][0] + b["seas"]["with_particle_release"][0],
            b["lakes"]["with_particle_release"][1] + b["seas"]["with_particle_release"][1]),
        decent=dict(kept_shares=kept, wild_dissolved_loss_kg_yr=wild_dissolved,
                    albatross_scale_flyers=flyers_needed,
                    lake_emergence_share_of_myvatn_top=pair(wild_dissolved[0] / top_rate[1],
                                                            wild_dissolved[1] / top_rate[0]),
                    extraction_from_deposits_kg_yr=rt["need_for_working_range_kg"],
                    harvest_ceiling_share=ceiling),
        sharing=budget, marking="derived")


# ---- Checks of the two statements and the ledger --------------------------------------------------------------------
def statements(part, sea, rt):
    ocean = sea["earth_ocean"]["mean_mg_m3"]
    per_tg_km3 = 1e9 / (ocean / 1e6) / 1e9
    flow_km3 = sea["discharge_m3_yr"] / 1e9
    at1000 = sea["at_1000_years_mg_m3"]
    lunar = {c: [1e9 / (x / 1e6) / 1e9 for x in reversed(v)] for c, v in at1000.items()}
    d_in = sea["input_cases_kg_yr"]["dissolved_after_lakes"]
    return dict(
        easy_route=dict(kept_by_lakes_share=part["kept_by_lakes_share_of_loss"],
                        to_seas_share=part["to_seas_share_of_loss"],
                        sea_mouth_deposits_kg_yr=part["to_seas_particulate_kg_yr"],
                        farm_ceiling_kg_yr=rt["seaweed"]["ceiling_kg_yr"]),
        open_water=dict(earth_ocean_mg_m3=ocean, km3_per_tg_at_earth_ocean=per_tg_km3, river_flow_km3_yr=flow_km3,
                        river_flows_per_tg=per_tg_km3 / flow_km3,
                        km3_per_tg_lunar_seas_1000_yr=lunar,
                        river_flows_per_tg_lunar_seas_1000_yr={c: [x / flow_km3 for x in v] for c, v in lunar.items()},
                        share_of_seas_volume_per_tg_1000_yr={c: [x * 1e9 / sea["volume_m3"] for x in v]
                                                             for c, v in lunar.items()},
                        teragram_over_dissolved_input=pair(1e9 / d_in[1], 1e9 / d_in[0]),
                        seas_volume_km3=sea["volume_m3"] / 1e9),
        marking="derived")


def ledger(land, lost_by_people, to_seas_share, carriers):
    """Steady-state P a year: the land, the waters, the people's loop and the carriers."""
    into_waters = land + lost_by_people
    to_seas = to_seas_share * into_waters
    to_lakes = into_waters - to_seas
    dredged = lost_by_people
    lake_burial = to_lakes - dredged - carriers
    if lake_burial < 0:
        raise ValueError("The lakes cannot supply this ledger")
    flows = dict(land_loss=land, people_losses=lost_by_people, into_waters=into_waters, to_lakes=to_lakes,
                 to_seas=to_seas, dredged_to_farms=dredged, carriers_to_wild_land=carriers,
                 lake_burial=lake_burial, sea_burial=to_seas)
    flows["closure_kg_yr"] = into_waters - (dredged + carriers + lake_burial + to_seas)
    flows["dredged_share_of_lake_input"] = dredged / to_lakes
    return flows


def ledgers(loss, part, folks, nat):
    """Every corner of the land loss against the population, for the design target and for field crops."""
    out = {}
    cases = dict(design_target=folks["ecology_design_target"]["for_working_range_kg"],
                 field_crops_all_unrecovered_to_waters=folks["systems"]["field_crops"]["new_for_working_range_kg"])
    for case, people_ends in cases.items():
        corners = {}
        for i, land_edge in enumerate(("low_loss", "high_loss")):
            for j, people_edge in enumerate(("4_billion", "8_billion")):
                corners[f"{land_edge}_{people_edge}"] = ledger(
                    loss["total_kg_yr"][i], people_ends[j], part["to_seas_share_of_loss"][i],
                    nat["decent"]["wild_dissolved_loss_kg_yr"][i])
        shares = [c["dredged_share_of_lake_input"] for c in corners.values()]
        out[case] = dict(corners=corners, dredged_share_of_lake_input=pair(min(shares), max(shares)))
    out["basis"] = ("Steady state: people's losses join the land's in the waters and are dredged back from the lake "
                    "deposits; carriers return the wild land's dissolved loss, 30% of the land kept wild at the low "
                    "land loss and half at the high; every other P is buried. The field-crop case sends all P the "
                    "crops do not recover to the waters, an upper bound while soils build stores.")
    return out


def run():
    data = {rel: load(rel) for rel in INPUTS}
    bound = {k: load_bound(k) for k in BOUND}
    people, sky = bound["ecology_people"], bound["sky_ecology"]
    ways, living = bound["ways_of_living"], bound["living_water"]
    waters = data["resources/results/sea_phosphorus_waters.json"]
    register = data["research/studies/lunar_cycle_ecology/results/lunar_cycle_ecology.json"]
    loss = land_loss(people, sky, register, data["research/studies/conservation/results/gates.json"], waters)
    part = partition(loss, waters, data["geography/results/drainage.json"], data["resources/results/first_screen.json"])
    sea = seas(part, loss, people, data["resources/results/first_screen.json"],
               data["resources/results/phosphorus.json"], ways, data["geography/results/atlas.json"], waters)
    gather = gathering(part, loss, waters, data["geography/results/drainage.json"], living)
    folks = people_p(people, ways, bound["provisioning_decisions"])
    rt = routes(part, sea, loss, folks, people, ways)
    dist = distribution(rt, ways, folks)
    nat = natural(loss, part, rt, people, sky, register, ways)
    producer = [ROOT / p for p in CODE]
    return dict(
        schema=SCHEMA,
        producer=dict(files={p: digest(ROOT / p) for p in CODE},
                      inputs={rel: digest(ROOT / rel) for rel in INPUTS},
                      bound_inputs={k: dict(path_at_integration=b["path"], branch=b["branch"], commit=b["commit"],
                                            sha256=b["sha256"], read_from="in-tree")
                                    for k, b in BOUND.items()},
                      pending_bindings=PENDING,
                      constants=constants_used(producer)),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        design_guesses=DESIGN_GUESSES,
        land_loss=loss, partition=part, seas=sea, gathering=gather, people=folks, routes=rt, distribution=dist,
        natural_cycle=nat, statements=statements(part, sea, rt), ledger=ledgers(loss, part, folks, nat))


def main():
    product = run()
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(product, indent=2, allow_nan=False) + "\n")
    print(OUT.relative_to(ROOT))


if __name__ == "__main__":
    main()
