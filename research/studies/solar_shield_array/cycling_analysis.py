"""Verify raw cycling runs; export reproducible trajectories and inventory scales."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re

import numpy as np
from scipy.interpolate import CubicHermiteSpline

from shared import constants as K
from shared.provenance import constants_changed, constants_used
from protection.dynamics.ephemeris import Ephemeris
from protection.dynamics.cycling import (CyclingEnvironment, ray_geometry,
    service_coast_angles, service_fraction, outgoing_clearance, sail_basis)
from .cycling_search import ROOT, HERE, RUN, SCENARIO, CENTRAL, ReturnProblem
from .cycling_feedback import source_identity


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def compact_series_json(product):
    """Keep metadata readable and numeric series on individual lines."""
    pretty = json.dumps(product, indent=2, allow_nan=False)
    compact = re.sub(r"\[\n[\s0-9eE+.,\-]+\n\s*\]",
                     lambda m: json.dumps(json.loads(m[0]), separators=(",", ":")), pretty)
    if json.loads(compact) != product:
        raise ValueError("Formatting changed product values")
    return compact+"\n"


def read_raw(name):
    path = RUN/f"{name}.npz"
    metadata = json.loads(path.with_suffix(".json").read_text())
    if metadata["producer"]["source_hashes"] != source_identity():
        raise ValueError(f"Stale dynamics source: {name}")
    if constants_changed(metadata["producer"]["constants"]):
        raise ValueError(f"Stale constants: {name}")
    if digest(path) != metadata["raw_sha256"]:
        raise ValueError(f"Raw trace hash mismatch: {name}")
    with np.load(path) as data:
        trace = {k: data[k] for k in data.files}
    return trace, metadata


def spline(trace):
    return CubicHermiteSpline(trace["t"], trace["state"][:, :3], trace["state"][:, 3:], extrapolate=False)


def episodes(t, active):
    """Sample-resolved episodes; each boundary has at most one sample uncertainty."""
    changes = np.diff(np.r_[False, active, False].astype(int))
    starts, stops = np.flatnonzero(changes == 1), np.flatnonzero(changes == -1)-1
    return [[float(t[a]), float(t[b])] for a, b in zip(starts, stops)]


def time_average(values, t):
    """Trapezoidal time mean, including Boolean indicator plateaus.

    NumPy adds adjacent ordinates before division by two. Boolean addition
    preserves dtype (True + True is True), so cast before that addition.
    """
    return float(np.trapezoid(np.asarray(values, dtype=float), t)/(t[-1]-t[0]))


def analyze(trace, env, sample_step=300.):
    t = np.linspace(trace["t"][0], trace["t"][-1], int(np.ceil(np.ptp(trace["t"])/sample_step))+1)
    p = spline(trace)
    state = np.column_stack([p(t), p(t, 1)])
    sample = env.at(t)
    rays = ray_geometry(state[:, :3], sample)
    angles = service_coast_angles(state, rays, sample, SCENARIO["protected_radii"]*K.MOON_RADIUS)
    diagnostic = env.acceleration(t, state[:, :3], angles, diagnostics=True)
    full = dict(t=t, state=state, angles=angles, **diagnostic)
    useful = service_fraction(full, SCENARIO["protected_radii"]*K.MOON_RADIUS)
    duration = t[-1]-t[0]
    avg = lambda y: time_average(y, t)
    radius = np.linalg.norm(state[:, :3], axis=-1)
    windows = episodes(t, useful > .001)
    gaps = [b[0]-a[1] for a, b in zip(windows[:-1], windows[1:])]
    active = diagnostic["projected"] > 1e-7
    margin = radius/K.SPEED_OF_LIGHT*np.linalg.norm(sample["moon_v"], axis=-1)+50.
    clearance = outgoing_clearance(state[:, :3], rays, angles,
                                    SCENARIO["protected_radii"]*K.MOON_RADIUS+margin)
    e, b1, b2 = sail_basis(rays["sun"])
    alpha, clock = angles.T
    normal = e*np.cos(alpha)[:, None]+np.sin(alpha)[:, None]*(
        b1*np.cos(clock)[:, None]+b2*np.sin(clock)[:, None])
    rotation = np.arccos(np.clip(np.abs(np.sum(normal[:-1]*normal[1:], axis=1)), 0, 1))
    rate = np.rad2deg(rotation)/np.diff(t)
    feathering = active[:-1] | active[1:]
    # The vectorized geometry must reproduce the controller used in propagation.
    original_rays = ray_geometry(trace["state"][:, :3], env.at(trace["t"]))
    original_angles = service_coast_angles(trace["state"], original_rays, env.at(trace["t"]),
                                          SCENARIO["protected_radii"]*K.MOON_RADIUS)
    if np.max(np.abs(original_angles-trace["angles"])) > 1e-9:
        raise ValueError("Unmodelled fallback in the service-and-coast control")
    result = dict(
        completed=bool(trace["completed"]), duration_days=duration/K.JULIAN_DAY,
        moon_range_km=[float(radius.min()/1000), float(radius.max()/1000)],
        useful_projected_area_fraction=avg(useful),
        raw_service_time_fraction=avg(useful > .001), service_visits=len(windows),
        service_window_hours=[(b-a)/3600 for a, b in windows],
        maximum_between_visit_gap_days=max(gaps, default=0)/K.JULIAN_DAY,
        incident_eclipse_fraction=avg(1-diagnostic["visibility"]),
        sail_impulse_m_s=float(np.trapezoid(np.linalg.norm(diagnostic["sail"], axis=-1), t)),
        electric_impulse_m_s=0.,
        minimum_reflected_beam_clearance_deg=float(np.rad2deg(clearance[active].min())),
        maximum_commanded_plane_rate_deg_s=float(rate.max()),
        maximum_rate_while_deployed_deg_s=float(rate[feathering].max()),
        maximum_sunward_photon_acceleration_m_s2=float(np.max(np.sum(diagnostic["sail"]*(-e), axis=-1))),
        sampling_seconds=float(t[1]-t[0]),
        central_diverted_fraction_max=float(np.max(diagnostic["fraction"][diagnostic["central"]]))
        if np.any(diagnostic["central"]) else None,
    )
    # One-hour display samples plus the first ten days at analysis cadence.
    every = max(1, int(round(7200/sample_step)))
    result["display"] = dict(day=np.round(t[::every]/K.JULIAN_DAY, 6).tolist(),
        range_km=np.round(radius[::every]/1000, 2).tolist(), useful=np.round(useful[::every], 6).tolist())
    day_index = np.floor(t/K.JULIAN_DAY).astype(int)
    minimum = np.full(day_index.max()+1, np.inf)
    maximum = np.full(day_index.max()+1, -np.inf)
    np.minimum.at(minimum, day_index, radius/1000)
    np.maximum.at(maximum, day_index, radius/1000)
    result["daily_range_km"] = dict(day=np.arange(len(minimum)).tolist(),
        minimum=np.round(minimum, 3).tolist(), maximum=np.round(maximum, 3).tolist())
    keep = (t <= 10*K.JULIAN_DAY) & (np.arange(len(t)) % max(1, int(round(600/sample_step))) == 0)
    solar_coordinates = np.column_stack([rays["sunward"], np.sum(rays["transverse"]*b1, axis=1),
                                         np.sum(rays["transverse"]*b2, axis=1)])
    result["first_ten_days"] = dict(day=np.round(t[keep]/K.JULIAN_DAY, 6).tolist(),
        solar_coordinates_km=np.round(solar_coordinates[keep]/1000, 2).tolist(),
        useful=np.round(useful[keep], 6).tolist(), cone_deg=np.round(np.rad2deg(alpha[keep]), 4).tolist())
    return result


def replay_arc(env, archive):
    data = np.load(archive)
    problem = ReturnProblem(env, float(data["radius"]), float(data["duration"]),
                            float(data["start"]), int(data["intervals"]))
    problem.anchor = data["anchor"]
    vector = data["vector"]
    propagated = problem.propagate(vector, step=600., rtol=2e-11)
    gravity_only = problem.propagate(vector, step=600., rtol=2e-11, enabled=False)
    state, angles = problem.unpack(vector)
    defect = problem.residual(vector, False)
    dt = propagated["t"][-1]-propagated["t"][0]
    result = dict(
        radius_seed_km=problem.radius/1000, duration_days=dt/K.JULIAN_DAY,
        optimizer_status="Evaluation budget reached; accepted by independent propagation tolerances, not optimizer success",
        closure_position_m=float(np.linalg.norm(propagated["local_closure"][:3])),
        closure_velocity_m_s=float(np.linalg.norm(propagated["local_closure"][3:])),
        gravity_only_counterfactual_closure_position_m=float(np.linalg.norm(gravity_only["local_closure"][:3])),
        gravity_only_counterfactual_closure_velocity_m_s=float(np.linalg.norm(gravity_only["local_closure"][3:])),
        max_scaled_dynamics_defect=float(np.max(np.abs(defect[:6*problem.intervals]))),
        max_scaled_soft_anchor_residual=float(np.max(np.abs(defect[-3:]))),
        useful_projected_area_fraction=float(np.trapezoid(service_fraction(propagated, 4*K.MOON_RADIUS),
                                                         propagated["t"])/dt),
        electric_impulse_m_s=0.,
        minimum_reflected_beam_clearance_deg=float(np.rad2deg(np.min(outgoing_clearance(
            propagated["state"][:, :3], propagated, propagated["angles"], 4*K.MOON_RADIUS)))),
        initial_icrf_state_m_m_s=(state[0]*problem.scale).tolist(),
        control_times_s=problem.times.tolist(), control_angles_rad=angles.tolist(),
        boundary_rule="Same state in instantaneous Earth–Moon frames at two dates; no periodic ephemeris is assumed",
        raw_archive_sha256=digest(archive),
    )
    if result["closure_position_m"] > 50 or result["closure_velocity_m_s"] > .001:
        raise ValueError("Return arc failed the independent 50 m / 1 mm/s closure limits")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arc", type=Path, default=RUN/"refined.npz")
    args = parser.parse_args()
    scenario = json.loads((HERE/"cycling_scenario.json").read_text())
    eph = Ephemeris(epoch=SCENARIO["epoch_tdb"])
    maximum_days = max(c["days"] for c in scenario["cases"])
    samples = eph.sample(np.arange(-3600, maximum_days*K.JULIAN_DAY+7200, 1800.))
    eph.close()
    cases, raw_inputs, traces = {}, {}, {}
    for config in scenario["cases"]:
        name = config["name"]
        trace, meta = read_raw(name)
        for key, value in config.items():
            if key != "name" and meta["config"][key] != value:
                raise ValueError(f"Raw config mismatch: {name}/{key}")
        traces[name] = trace
        env = CyclingEnvironment(samples, SCENARIO["baseline_areal_mass_kg_m2"], CENTRAL, config["spectrum"])
        if config["mode"] == "service_coast":
            cases[name] = {"config": config, **analyze(trace, env)}
        else:
            cases[name] = {"config": config, **meta["summary"],
                           "reading_rule": "Rejected comparison; no long-term feasibility or inventory inferred"}
        raw_inputs[name] = dict(trace=meta["raw_sha256"], metadata=digest(RUN/f"{name}.json"))
        print(f"Analyzed {name}", flush=True)
    env = CyclingEnvironment(samples, .05, CENTRAL, "filter_only")
    reference, fine, perturbed = (traces[k] for k in ("inner_filter", "inner_filter_fine", "inner_filter_perturbed"))
    t = fine["t"]
    position_difference = np.linalg.norm(spline(reference)(t)-fine["state"][:, :3], axis=-1)
    perturbation_distance = np.linalg.norm(spline(reference)(perturbed["t"])-perturbed["state"][:, :3], axis=-1)
    finer_sampling = analyze(fine, env, 120.)
    coarse_fraction = cases["inner_filter"]["useful_projected_area_fraction"]
    fine_fraction = cases["inner_filter_fine"]["useful_projected_area_fraction"]
    convergence = dict(
        maximum_position_difference_m=float(position_difference.max()),
        final_position_difference_m=float(position_difference[-1]),
        useful_fraction_relative_difference=abs(fine_fraction/coarse_fraction-1),
        service_quadrature_300_to_120_s_relative_difference=abs(
            finer_sampling["useful_projected_area_fraction"]/fine_fraction-1),
        maximum_perturbation_separation_m=float(perturbation_distance.max()),
        final_perturbation_separation_m=float(perturbation_distance[-1]),
        perturbed_useful_fraction_relative_change=cases["inner_filter_perturbed"]["useful_projected_area_fraction"]/coarse_fraction-1,
        perturbation="1 m in ICRF x and 1 mm/s in ICRF y at initial epoch",
        interpretation="Bounded orbital response is separate from seam/formation stability; no attitude actuator dynamics were integrated",
    )
    if convergence["useful_fraction_relative_difference"] > .005 or convergence["service_quadrature_300_to_120_s_relative_difference"] > .005:
        raise ValueError("Useful-area statistic has not converged to 0.5%")
    holding_path = HERE/"results/holding.json"
    held = json.loads(holding_path.read_text())["aperture_quadrature"]["variable_distance_finer"]
    yearly_fuel = held["propellant_kg_s"]*K.JULIAN_DAY*K.JULIAN_YEAR_DAYS
    inventory = {}
    for name in ("inner_filter", "middle_filter", "middle_full"):
        duty = cases[name]["useful_projected_area_fraction"]
        inventory[name] = dict(reference_aperture_area_multiplier=1/duty,
            base_optical_inventory_kg=held["base_optical_mass_kg"]/duty,
            inventory_in_reference_propellant_years=held["base_optical_mass_kg"]/duty/yearly_fuel)
    arc = replay_arc(CyclingEnvironment(samples, .05, CENTRAL), args.arc)
    sources = source_identity()
    sources[str(Path(__file__).relative_to(ROOT))] = digest(__file__)
    sources[str((HERE/"cycling_scenario.json").relative_to(ROOT))] = digest(HERE/"cycling_scenario.json")
    # Keep a compact display for only the primary case, with daily ranges for others.
    for name, case in cases.items():
        if name != "inner_filter":
            case.pop("first_ten_days", None)
            if "display" in case:
                case["display"] = {k: v[::24] for k, v in case["display"].items()}
    product = dict(
        schema="terluna.research.solar-shield-cycling/1",
        producer=dict(runner=str(Path(__file__).relative_to(ROOT)), source_hashes=sources,
                      constants=constants_used(sources), raw_inputs=raw_inputs,
                      holding_product_sha256=digest(holding_path)),
        evidence="Numerically propagated ideal point-sail trajectories in DE440 point-mass dynamics; not empirical or materials validation",
        epoch_tdb=SCENARIO["epoch_tdb"], areal_mass_kg_m2=.05,
        reading_rule="Annual boundedness and shadow-service statistics concern individual tiles. No fleet population, spatial handover, collision avoidance or continuous UV coverage is demonstrated. The inventory is a reference-aperture equivalent under perfect area placement, not a solved fleet or cost payback. Full finite-tile retarded intersections, actual spectral/angular optics, attitude hardware, power collection and habitat loads remain open.",
        physical_model=dict(gravity="Moon plus Sun, Earth and DE440 planetary systems; subtract ephemeris Moon acceleration",
            actuator="Ideal spectral specular membrane; a=2 f S cos(alpha)^2 n/(c sigma); transmitted remainder unchanged",
            eclipses="Union of finite solar disk occultations by Earth and Moon",
            rays="Second-order retarded sources/blockers; target motion during tile-to-Moon flight; infinitesimal patch footprint",
            omitted=["Lunar/Earth multipoles", "relativity", "limb/spectral variation", "non-specular/absorptive optical force",
                     "attitude dynamics and hardware mass", "finite-tile shape and seams", "deployment and replacements"]),
        cases=cases, convergence=convergence, independently_replayed_return_arc=arc,
        reference_area_inventory=inventory,
        reference_holding=dict(aperture_area_m2=held["aperture_area_m2"], base_optical_mass_kg=held["base_optical_mass_kg"],
                               mean_power_TW=held["mean_power_TW"], annual_expended_propellant_kg=yearly_fuel),
        fleet_coverage_demonstrated=False,
        units=dict(position="m unless named km", velocity="m/s", area="m2", mass="kg", angle="radians unless named deg", time="seconds or TDB days as named"),
    )
    path = HERE/"results/cycling.json"
    path.write_text(compact_series_json(product))
    print(json.dumps({"product": str(path.relative_to(ROOT)), "convergence": convergence,
                      "inventory": inventory, "return_closure_m": arc["closure_position_m"]}), flush=True)


if __name__ == "__main__":
    main()
