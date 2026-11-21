"""The salt for an ocean-like sea: what it takes and what the planned imports carry.

Run ``python -m resources.sea_salt``. On 15 November 2026 the author asked to test the viability and the cost and
benefit of making the South Pole-Aitken Sea an Earth-like ocean. Lunar rock gives the water almost no chloride, so
the sea's salt is an import question. This screen sets the sea's chloride and sodium beside the chlorine that rides
in the feed the project already plans to process, the water and the air's nitrogen, by the kind of body the feed
comes from, and beside lunar rock and the near-Earth asteroids of the source-access study.

It is arithmetic on published grades and project inventories. The grades are group means of meteorites and one
comet's gas; no body the project would mine has been assayed for chlorine, and the processing that would keep the
chlorine is not designed here.
"""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

from shared.constants import GRAVITATIONAL_CONSTANT, JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_GM
from shared.provenance import constants_used

ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent
OUT = HERE / "results/sea_salt.json"
SOURCES = HERE / "sea_salt_sources.json"
SCREEN = HERE / "results/first_screen.json"
ATLAS = ROOT / "geography/results/atlas.json"
POTENTIAL = HERE / "source_access/results/potential.json.gz"
SCHEMA = "terluna.resources.sea-salt/1"
SEAS = {433: "near_side", 1417: "south_pole_aitken", 874: "smythii_marginis", 213: "humboldtianum",
        1296: "orientale", 714: "moscoviense"}
SEA = "south_pole_aitken"
YEAR_S = JULIAN_DAY * JULIAN_YEAR_DAYS
MOON_MASS = MOON_GM / GRAVITATIONAL_CONSTANT
ASTEROID_HORIZON = 200                  # the vacuum period's checkpoint the feed is read at, years
ASTEROID_POOLS = ("narrow", "wider", "all_no_don_quixote")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def span(low, high):
    return [min(low, high), max(low, high)]


def need(volume_km3, values, share=1.0):
    """Salt, chloride and sodium (kg) in a volume of water at a share of seawater's strength. The density is
    seawater's at every share, which overstates a half-strength sea's mass by about 1%."""
    sea = values["seawater_g_per_kg"]
    water_kg = volume_km3 * 1e9 * values["seawater_density_kg_m3"]["value"]
    return dict(volume_km3=volume_km3, share_of_seawater=share,
                salt_kg=water_kg * share * sea["total"] * 1e-3,
                cl_kg=water_kg * share * sea["cl"] * 1e-3, na_kg=water_kg * share * sea["na"] * 1e-3)


def riding_chlorine(values, screen):
    """Chlorine (kg) carried by the feed for the water and for the air's nitrogen, by the kind of body it comes
    from. Each entry gives low, central and high where the grade has them, else the cases named."""
    ci, comet = values["ci_rock"], values["comet_ice"]
    water = screen["inventories"]["water"]
    standing = water["seas_kg"] + water["lakes_kg"]
    with_crust = water["total_kg"]
    n2 = screen["inventories"]["n2_kg"]
    water_share = ci["h_ppm"]["value"] * 1e-6 * 18.015 / 2.016      # all of the rock's hydrogen counted as water
    cl = ci["cl_ppm"]
    ci_per_water = {k: cl[k] * 1e-6 / water_share for k in ("low", "central", "high")}
    ci_per_n = dict(low=cl["low"] / ci["n_ppm"]["high"], central=cl["central"] / ci["n_ppm"]["high"],
                    high=cl["high"] / ci["n_ppm"]["low"])
    disputed = ci["cl_ppm_disputed_low"]["value"]
    per_water = comet["cl_per_water_by_mass"]
    return dict(
        water_from_comet_like_ice=dict(
            basis="gas-phase HCl of comet 67P per kg of water; a lower bound on the nucleus",
            cl_per_kg_water=per_water,
            standing_water_kg=standing,
            cl_kg={k: per_water[k] * standing for k in ("low", "central", "high")},
            cl_kg_with_the_crusts_water=span(per_water["central"] * with_crust[0], per_water["central"] * with_crust[1]),
            cl_kg_at_the_solar_ceiling=comet["cl_per_water_solar_ceiling"]["value"] * standing),
        water_from_ci_like_rock=dict(
            basis="group-mean chlorine of CI chondrites per kg of the water their hydrogen makes",
            water_share_of_rock=water_share, cl_per_kg_water=ci_per_water,
            rock_kg=standing / water_share,
            cl_kg={k: ci_per_water[k] * standing for k in ("low", "central", "high")},
            cl_kg_at_the_disputed_low_grade=disputed * 1e-6 / water_share * standing),
        nitrogen_from_comet_like_ice=dict(
            basis="chlorine over nitrogen in comet 67P: gas over gas, and gas chlorine over gas-and-dust nitrogen",
            n2_kg=n2,
            cl_kg_gas_over_gas=comet["cl_per_n_gas"]["value"] * n2,
            cl_kg_floor=span(comet["cl_per_n_gas_over_gas_and_dust_nitrogen"]["low"] * n2,
                             comet["cl_per_n_gas_over_gas_and_dust_nitrogen"]["high"] * n2)),
        nitrogen_from_ci_like_rock=dict(
            basis="group-mean chlorine over nitrogen of CI chondrites",
            cl_per_kg_n=ci_per_n, rock_kg=span(n2 / (ci["n_ppm"]["high"] * 1e-6), n2 / (ci["n_ppm"]["low"] * 1e-6)),
            cl_kg={k: ci_per_n[k] * n2 for k in ("low", "central", "high")},
            cl_kg_at_the_disputed_low_grade=disputed / ci["n_ppm"]["high"] * n2))


def rock_for_the_chloride(cl_kg, na_kg, values):
    """CI-like rock that holds a sea's chloride, the sodium that comes with it, and the same from lunar ground."""
    ci, moon, project = values["ci_rock"], values["moon"], values["project"]
    rock = {k: cl_kg / (ci["cl_ppm"][k] * 1e-6) for k in ("low", "central", "high")}
    rock_disputed = cl_kg / (ci["cl_ppm_disputed_low"]["value"] * 1e-6)
    planned = project["small_bodies_processed_for_the_air_kg"]
    sodium = rock["central"] * ci["na_ppm"]["value"] * 1e-6
    soil_na = moon["na_weight_fraction_highland_soil"]
    lunar_cl = moon["cl_ppm_soil"]
    return dict(
        ci_like_rock_kg=dict(at_800_ppm=rock["high"], at_717_ppm=rock["central"], at_600_ppm=rock["low"],
                             at_the_disputed_115_ppm=rock_disputed),
        share_of_the_feed_planned_for_the_air=dict(
            planned_kg=[planned["low"], planned["high"]],
            central_grade=span(rock["central"] / planned["high"], rock["central"] / planned["low"]),
            disputed_low_grade=span(rock_disputed / planned["high"], rock_disputed / planned["low"])),
        share_of_the_kuiper_belts_small_bodies=rock["central"] / project["kuiper_belt_small_bodies_kg"]["value"],
        sodium_in_that_rock_kg=sodium, sodium_over_the_seas_need=sodium / na_kg,
        na_over_cl_in_ci_rock=ci["na_ppm"]["value"] / ci["cl_ppm"]["central"],
        na_over_cl_in_seawater=values["seawater_g_per_kg"]["na"] / values["seawater_g_per_kg"]["cl"],
        lunar_soil_for_the_sodium_kg=span(na_kg / soil_na["high"], na_kg / soil_na["low"]),
        lunar_ground_for_the_chloride_kg=span(cl_kg / (lunar_cl["high"] * 1e-6), cl_kg / (lunar_cl["low"] * 1e-6)),
        lunar_ground_for_the_chloride_share_of_the_moon=span(cl_kg / (lunar_cl["high"] * 1e-6) / MOON_MASS,
                                                             cl_kg / (lunar_cl["low"] * 1e-6) / MOON_MASS))


def near_earth_asteroids(potential, values, cl_need_kg):
    """Chlorine in the near-Earth asteroid feed and pools of the source-access study, if all of it were carbonaceous
    rock at the CM to CI group means and every atom were kept: an upper bound, since most near-Earth bodies are
    ordinary-chondrite-like at 87-126 ppm."""
    sea = values["seawater_g_per_kg"]
    cl_per_km3 = values["seawater_density_kg_m3"]["value"] * 1e9 * sea["cl"] * 1e-3
    grades = dict(cm=values["other_chondrites_cl_ppm"]["cm"], ci=values["ci_rock"]["cl_ppm"]["central"],
                  ordinary=values["other_chondrites_cl_ppm"]["l"])
    feed = next(b["ranges"]["source_feed_kg"] for b in potential["baseline_demand_ranges"]
                if b["horizon"] == ASTEROID_HORIZON)
    pools = {}
    for row in potential["comparisons"]:
        pool = row.get("pool") or row.get("population_pool")
        if row["available_share"] == 1.0 and pool in ASTEROID_POOLS and pool not in pools:
            pools[pool] = row["feed_mass"]["available_stock_kg_p05_p50_p95"]

    def chlorine(mass_kg):
        cl = {k: mass_kg * g * 1e-6 for k, g in grades.items()}
        return dict(cl_kg=cl, seawater_km3={k: v / cl_per_km3 for k, v in cl.items()},
                    share_of_the_seas_need={k: v / cl_need_kg for k, v in cl.items()})

    return dict(grades_ppm=grades, cl_kg_per_km3_of_seawater=cl_per_km3,
                feed_through_year_200=dict(feed_kg=feed, low=chlorine(feed[0]), high=chlorine(feed[1])),
                pools={k: dict(body_mass_kg_p05_p50_p95=v, at_the_median=chlorine(v[1])) for k, v in pools.items()})


def delivery(salt_kg, screen):
    """The salt beside the planned delivery of air and water (the first screen's stream over its build time)."""
    d = screen["delivery"]
    return dict(delivery_mass_kg=d["mass_kg"], delivery_flow_kg_s=d["flow_kg_s"], build_years=d["build_years"],
                salt_share_of_the_delivery=span(salt_kg / d["mass_kg"][1], salt_kg / d["mass_kg"][0]),
                years_of_the_stream=span(salt_kg / d["flow_kg_s"][1] / YEAR_S, salt_kg / d["flow_kg_s"][0] / YEAR_S),
                salt_share_of_the_standing_water=salt_kg / (screen["inventories"]["water"]["seas_kg"]
                                                            + screen["inventories"]["water"]["lakes_kg"]))


def run():
    sources = json.loads(SOURCES.read_text())
    values = sources["values"]
    screen = json.loads(SCREEN.read_text())
    atlas = json.loads(ATLAS.read_text())
    with gzip.open(POTENTIAL, "rt") as f:
        potential = json.load(f)
    volumes = {SEAS[b["label"]]: float(b["volume_km3"]) for b in atlas["bodies"] if b["label"] in SEAS}
    needs = {name: need(v, values) for name, v in volumes.items()}
    own = needs[SEA]
    riding = riding_chlorine(values, screen)

    def shares(entry):
        out = {}
        for key, value in entry.items():
            if not key.startswith("cl_kg"):
                continue
            label = "share_of_need" + key[len("cl_kg"):]
            if isinstance(value, dict):
                out[label] = {k: v / own["cl_kg"] for k, v in value.items()}
            elif isinstance(value, list):
                out[label] = [v / own["cl_kg"] for v in value]
            else:
                out[label] = value / own["cl_kg"]
        return out

    for entry in riding.values():
        entry.update(shares(entry))
    producer = Path(__file__)
    inputs = (SOURCES, SCREEN, ATLAS, POTENTIAL)
    return dict(
        schema=SCHEMA,
        producer=dict(files={str(producer.relative_to(ROOT)): digest(producer)},
                      inputs={str(p.relative_to(ROOT)): digest(p) for p in inputs},
                      constants=constants_used([producer])),
        evidence="Arithmetic on published grades and project inventories; nothing is modelled. Chlorine grades are "
                 "group means of CI chondrites (read in the compilations; one 2017 study, rejected by the compilers, "
                 "gives a grade six times lower and is carried as a disputed low case) and the gas of one comet, "
                 "67P, which bounds a nucleus's chlorine from below. No body the project would mine has been assayed "
                 "for chlorine. Shares assume every atom of chlorine in the processed feed is kept; the processing "
                 "that would keep it is not designed. Lunar soil chlorine is inferred from a general statement, not "
                 "measured values. The sea's need is its fill at seawater's strength, without the salt the crust "
                 "takes with the water it draws from the floor.",
        reading_rule="Masses in kg; shares are of the South Pole-Aitken Sea's chloride at seawater's strength unless "
                     "named otherwise. low, central and high follow the grade's range; a two-entry list is (low, "
                     "high). In riding, each route gives the chlorine carried by all of that import if it came "
                     "wholly from that kind of body; the routes are alternatives for one import and add across "
                     "imports. Near-Earth asteroid figures are upper bounds at carbonaceous grades.",
        sea=SEA, needs=needs,
        other_strengths={name: need(volumes[SEA], values, share) for name, share in (("half", 0.5), ("fifth", 0.2))},
        riding=riding,
        rock=rock_for_the_chloride(own["cl_kg"], own["na_kg"], values),
        near_earth_asteroids=near_earth_asteroids(potential, values, own["cl_kg"]),
        delivery=delivery(own["salt_kg"], screen))


def main():
    product = run()
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(product, indent=2, allow_nan=False) + "\n")
    print(OUT.relative_to(ROOT))


if __name__ == "__main__":
    main()
