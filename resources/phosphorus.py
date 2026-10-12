"""Conservative phosphorus return requirements and finite-reservoir sensitivities.

Run ``python -m resources.phosphorus``. This is an accounting/box-model screen;
its rates, stocks and biomass architecture are scenarios, not an evolved ecology.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_GM, MOON_RADIUS
from shared.provenance import constants_used

ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent
OUT = HERE / "results/phosphorus.json"
SOURCES = HERE / "phosphorus_sources.json"
POPULATION = ROOT / "shared/scenarios/population.json"
INDUSTRY = HERE / "results/industry.json"
SCHEMA = "terluna.resources.phosphorus/1"
YEAR_S = JULIAN_DAY * JULIAN_YEAR_DAYS
AREA_M2 = 4 * np.pi * MOON_RADIUS**2
# Atomic weights for stoichiometric conversion; same convention as industry.py.
C_ATOMIC = 12.011
P_ATOMIC = 30.974
BOXES = ("sky", "land", "sea", "sediment", "human_return", "unrecoverable")


def fraction(value):
    if not 0 <= value <= 1:
        raise ValueError("Fraction must be between zero and one")
    return value


def stock_and_flux(cp_molar=106.0, coverage=0.01, npp_g_c_m2_yr=1000.0,
                   carbon_dry_fraction=0.45, standing_dry_kg_m2=(1.0, 10.0)):
    """Coverage is projected local footprint / 4*pi*R^2, never pi*R^2."""
    fraction(coverage)
    if min(cp_molar, npp_g_c_m2_yr, carbon_dry_fraction, *standing_dry_kg_m2) <= 0:
        raise ValueError("Positive stoichiometry, productivity and stock required")
    cp_mass = cp_molar * C_ATOMIC / P_ATOMIC
    area = AREA_M2 * coverage
    carbon = area * npp_g_c_m2_yr / 1000
    stock = [area * s * carbon_dry_fraction / cp_mass for s in standing_dry_kg_m2]
    return dict(cp_molar=cp_molar, cp_mass=cp_mass, coverage=coverage,
                occupied_surface_footprint_m2=area, carbon_kg_yr=carbon,
                dry_kg_yr=carbon / carbon_dry_fraction, throughput_p_kg_yr=carbon / cp_mass,
                standing_p_kg=stock, stock_over_throughput_yr=[s / (carbon / cp_mass) for s in stock])


def budget(throughput, harvest, natural_recycling, food_return, upward):
    """r applies only to unharvested P; q includes capture, treatment and delivery."""
    for f in (harvest, natural_recycling, food_return):
        fraction(f)
    if min(throughput, upward) < 0:
        raise ValueError("Negative flow")
    harvest_p = harvest * throughput
    natural_down = (1 - harvest) * (1 - natural_recycling) * throughput
    food_unreturned = (1 - food_return) * harvest_p
    loss = natural_down + food_unreturned
    return dict(harvest_p_kg_yr=harvest_p, natural_downward_p_kg_yr=natural_down,
                direct_food_return_p_kg_yr=food_return * harvest_p,
                food_unreturned_to_sky_p_kg_yr=food_unreturned,
                sky_export_after_immediate_returns_p_kg_yr=loss,
                upward_from_surface_p_kg_yr=upward,
                additional_return_needed_p_kg_yr=max(0.0, loss - upward),
                sky_net_p_kg_yr=upward - loss)


def required_natural_recycling(throughput, harvest, food_return, upward):
    """Raw algebraic threshold; below zero means no recycling needed, above one impossible."""
    for f in (harvest, food_return):
        fraction(f)
    if throughput <= 0 or upward < 0 or harvest == 1:
        raise ValueError("Positive throughput and an unharvested share required")
    r = 1 - (upward / throughput - harvest * (1 - food_return)) / (1 - harvest)
    return dict(required_fraction=r, feasible_at_or_below_unity=r <= 1,
                minimum_food_return_at_perfect_natural_recycling=(
                    max(0.0, 1 - upward / (harvest * throughput)) if harvest else 0.0))


def transition_model(stock, throughput, upward_initial, natural_recycling=0.999,
                     food_return=0.99, sediment_return_years=None,
                     irreversible_fraction=0.0, source_stock_multiple=10.0):
    """dx/dt=M*x, kg P and years, exactly conservative including the loss box.

    Sky throughput scales with its remaining stock. Instant recycling stays in sky.
    Surface-source stock is finite. Human return has a 0.1-year residence time.
    Rates and reservoir sizes are illustrative, explicitly not lunar field values.
    """
    for f in (natural_recycling, food_return, irreversible_fraction):
        fraction(f)
    if min(stock, throughput, source_stock_multiple) <= 0 or upward_initial < 0:
        raise ValueError("Invalid stock or rate")
    if sediment_return_years is not None and sediment_return_years <= 0:
        raise ValueError("Positive sediment return timescale required")
    initial = np.array([stock, stock * source_stock_multiple * 0.72,
                        stock * source_stock_multiple * 0.28, 0.0, 0.0, 0.0])
    matrix = np.zeros((len(BOXES), len(BOXES)))
    edges = []

    def transfer(source, dest, rate):
        i, j = BOXES.index(source), BOXES.index(dest)
        matrix[j, i] += rate
        matrix[i, i] -= rate
        edges.append(dict(source=source, destination=dest, rate_per_year=rate))

    turnover = throughput / stock
    transfer("sky", "human_return", 0.05 * turnover)
    down = 0.95 * (1 - natural_recycling) * turnover
    transfer("sky", "land", down * 0.72)
    transfer("sky", "sea", down * 0.28)
    transfer("human_return", "sky", food_return / 0.1)
    transfer("human_return", "land", (1 - food_return) * 0.72 / 0.1)
    transfer("human_return", "sea", (1 - food_return) * 0.28 / 0.1)
    # A finite donor loses exactly the phosphorus that the recipient receives.
    for box in ("land", "sea"):
        rate = upward_initial / (stock * source_stock_multiple)
        transfer(box, "sky", rate * (1 - irreversible_fraction))
        transfer(box, "unrecoverable", rate * irreversible_fraction)
    transfer("land", "sea", 1 / 100.0)
    transfer("sea", "sediment", 1 / 1000.0)
    if sediment_return_years is not None:
        transfer("sediment", "sky", (1 - irreversible_fraction) / sediment_return_years)
        transfer("sediment", "unrecoverable", irreversible_fraction / sediment_return_years)
    return initial, matrix, edges


def evolve(initial, matrix, years):
    if years < 0:
        raise ValueError("Negative elapsed time")
    return expm(matrix * years) @ initial


def lifting(kg_p_yr, altitude_m, p_mass_fraction, efficiency=0.7):
    """Gravity-only work on transported mass; no credit for return recovery."""
    if min(kg_p_yr, altitude_m) < 0 or not 0 < p_mass_fraction <= 1 or not 0 < efficiency <= 1:
        raise ValueError("Invalid transport input")
    work = MOON_GM * (1 / MOON_RADIUS - 1 / (MOON_RADIUS + altitude_m))
    carrier = kg_p_yr / p_mass_fraction
    return dict(altitude_m=altitude_m, p_mass_fraction=p_mass_fraction,
                transported_kg_yr=carrier, ideal_j_kg=work,
                efficiency=efficiency, power_w=carrier * work / (YEAR_S * efficiency))


def population_account():
    """Use historical coefficients with new common populations; no capacity claim."""
    scenarios = json.loads(POPULATION.read_text())
    old = json.loads(INDUSTRY.read_text())
    coefficients = old["people"]
    stock_t = coefficients["stock_t_per_person"]
    lifetime = coefficients["lifetime_years"]
    sp = coefficients["sp413"]
    out = {}
    for name, case in scenarios["cases"].items():
        assert abs(sum(case[k] for k in ("earth", "lunar_surface", "lunar_aerial", "orbital", "elsewhere")) - case["total"]) < 1e-12
        lunar = (case["lunar_surface"] + case["lunar_aerial"]) * 1e9
        orbital = case["orbital"] * 1e9
        stocks = [lunar * x * 1e3 for x in stock_t]
        out[name] = dict(population_billion=case,
                        lunar_generic_material_proxy_kg=stocks,
                        lunar_generic_turnover_kg_yr=[stocks[0] / max(lifetime), stocks[1] / min(lifetime)],
                        orbital_shield_proxy_kg=orbital * sp["shield_t"] / sp["people"] * 1e3,
                        orbital_structure_proxy_kg=orbital * sp["structure_t"] / sp["people"] * 1e3)
    return dict(status="comparison only; no selected allocation or carrying capacity",
                coefficients_from="resources/results/industry.json",
                exclusions="Earth and elsewhere infrastructure; aerial gas, envelope and flight structure. Generic material proxies overlap settlement district inventories and must not be added wholesale to them. SP-413 shielding applies only to the orbital habitat comparator, not atmospheric residents.",
                cases=out)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run():
    central = stock_and_flux()
    throughput = central["throughput_p_kg_yr"]
    upward = [AREA_M2 * d / 1000 * 0.1 for d in (0.0007, 0.028)]
    requirements = []
    grid = []
    for q in (0.0, 0.9, 0.99, 0.999):
        for u in upward:
            requirements.append(dict(food_return=q, upward_p_kg_yr=u,
                                     **required_natural_recycling(throughput, 0.05, q, u)))
            for r in (0.99, 0.999, 0.9999):
                grid.append(dict(food_return=q, natural_recycling=r,
                                 **budget(throughput, 0.05, r, q, u)))
    finite = []
    for label, sediment_years, loss in (("no_sediment_return", None, 0),
                                       ("sediment_return_no_irreversible_loss", 10000, 0),
                                       ("sediment_return_0.1percent_loss_per_processed_pass", 10000, 0.001)):
        for u in upward:
            for stock in central["standing_p_kg"]:
                initial, matrix, edges = transition_model(stock, throughput, u,
                                                         sediment_return_years=sediment_years,
                                                         irreversible_fraction=loss)
                states = []
                for years in (0, 100, 1000, 10000, 1e6, 1e9):
                    state = evolve(initial, matrix, years)
                    states.append(dict(year=years, p_kg=dict(zip(BOXES, state.tolist())),
                                       sky_relative_to_initial=state[0] / stock,
                                       total_relative_error=(state.sum() - initial.sum()) / initial.sum()))
                finite.append(dict(case=label, initial_p_kg=dict(zip(BOXES, initial.tolist())),
                                   upward_initial_p_kg_yr=u, edges=edges, states=states))
    # Complete support of a constant target throughput would require this extra
    # transfer from a reserve. This is the reserve duration, not a biomass lifetime.
    reserves = []
    for row in grid:
        gap = row["additional_return_needed_p_kg_yr"]
        reserves.append(dict(food_return=row["food_return"], natural_recycling=row["natural_recycling"],
                             upward_p_kg_yr=row["upward_from_surface_p_kg_yr"],
                             extra_return_p_kg_yr=gap, reserve_p_kg=1e12,
                             reserve_duration_years=None if gap == 0 else 1e12 / gap,
                             reserve_required_for_1e9_years_kg=gap * 1e9))
    producer = Path(__file__)
    inputs = (SOURCES, POPULATION, INDUSTRY)
    return dict(schema=SCHEMA,
                producer=dict(files={str(producer.relative_to(ROOT)): digest(producer)},
                              inputs={str(p.relative_to(ROOT)): digest(p) for p in inputs},
                              constants=constants_used([producer])),
                evidence="Closed-form requirements and six-box linear conservation sensitivities. No actual closed ecology, ore reserve, capture efficiency, sediment return technology or sustained productivity is established.",
                reading_rule="All phosphorus in kg elemental P; annual flow kg P per Julian year; stock differs from throughput. Null reserve duration means no extra reserve drain under a fixed assumed upward flow, not infinite planetary sustainability. Deposited P is internal redistribution, not fresh off-world import. Box-model unrecoverable means outside the defined recovery system, not nuclear destruction.",
                assumptions=dict(coverage_fraction_of_lunar_surface=0.01, npp_g_c_m2_yr=1000,
                                 biomass_carbon_dry_fraction=0.45, harvest_fraction=0.05,
                                 deposition_scenario_g_p_m2_yr=[0.0007, 0.028], capture_and_availability_fraction=0.1,
                                 natural_recycling_grid=[0.99, 0.999, 0.9999], food_return_grid=[0, 0.9, 0.99, 0.999],
                                 finite_model="Land+sea donor stock=10 times initial sky stock; donor split72:28, not inferred from measured P geography. Runoff timescale100yr, sea burial1000yr, engineered sediment return10000yr, human handling0.1yr. Natural recycling.999, food return.99. No imports; sky throughput scales with stock. Loss fraction applies to processed upward/engineered return, not all P turnover."),
                central=central, stoichiometry=[stock_and_flux(cp) for cp in (78, 106, 137, 195, 500)],
                upward_p_kg_yr=upward, recycling_requirements=requirements, budget_grid=grid,
                constant_target_reserve_sensitivity=reserves, finite_reservoir_cases=finite,
                transport=[lifting(throughput * 0.05, z, concentration) for z in (10000, 20000)
                           for concentration in (0.1, 0.001, 1e-5)],
                population=population_account())


def main():
    product = run()
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(product, indent=2, allow_nan=False) + "\n")
    print(OUT.relative_to(ROOT))


if __name__ == "__main__":
    main()
