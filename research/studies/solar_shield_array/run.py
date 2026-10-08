"""Full-ephemeris holding, trajectory search and mass/power closure.

Run from the repository root. Raw traces and the external kernel stay in runs/;
the compact result product carries source and input hashes.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize, brentq

from shared import constants as K
from shared.provenance import constants_used
from protection.dynamics.ephemeris import Ephemeris, DEFAULT_KERNEL
from protection.dynamics.model import geometry, evaluate, moon_acceleration_residual
from protection.dynamics.optical import length, optical_projection, axial_optical, covering_radius, solar_visibility

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCENARIO_PATH = HERE / "scenario.json"
SCENARIO = json.loads(SCENARIO_PATH.read_text())
CORE_DIVERTED = 1-SCENARIO["climate_flux_before_dimming_W_m2"]*SCENARIO["climate_sunlight_scale"]/K.SOLAR_CONSTANT


def constant_coeff(distance_km):
    c = np.zeros((3, 1+2*SCENARIO["harmonics"]))
    c[0, 0] = distance_km
    return c


def average(value, t):
    return float(np.trapezoid(value, t)/(t[-1]-t[0]))


def mixed_residual(case, sigma, core_fraction, mode="ideal"):
    if mode == "none":
        return case["required"], case["required"]
    if mode == "ideal":
        core = optical_projection(case["required"], case["sun"], sigma, CORE_DIVERTED, case["visibility"])[1]
        annulus = optical_projection(case["required"], case["sun"], sigma, 1, case["visibility"])[1]
    elif mode in ("absorb", "reflect"):
        coefficient = 1 if mode == "absorb" else 2
        core = case["required"]-axial_optical(case["sun"], sigma, coefficient*CORE_DIVERTED, case["visibility"])
        annulus = case["required"]-axial_optical(case["sun"], sigma, coefficient, case["visibility"])
    else:
        raise ValueError(mode)
    return core, annulus


def summarize(case, samples, sigma, mode="ideal", protected_radii=4):
    t = samples["t"]
    radius = float(case["covering_radius"].max()+SCENARIO["formation_margin_m"])
    # Conservative fixed central window, including the changing Sun's angular size.
    core_radius = float(covering_radius(case["coords"][:, 0], length(samples["positions"]["sun"]),
                                      K.MOON_RADIUS, length(case["coords"][:, 1:])).max()+case["light_time_allowance"].max())
    core_fraction = min((core_radius/radius)**2, 1)
    core, annulus = mixed_residual(case, sigma, core_fraction, mode)
    magnitudes = core_fraction*length(core)+(1-core_fraction)*length(annulus)
    mean = average(magnitudes, t)
    peak = float(magnitudes.max())
    req = average(length(case["required"]), t)
    return {"areal_mass_kg_m2": sigma, "optical_mode": mode,
            "protected_radii": protected_radii, "aperture_radius_km": radius/1000,
            "aperture_area_m2": float(np.pi*radius**2), "core_area_fraction": core_fraction,
            "mean_required_mm_s2": req*1000, "mean_residual_mm_s2": mean*1000,
            "peak_residual_mm_s2": peak*1000,
            "residual_delta_v_km_s": mean*(t[-1]-t[0])/1000,
            "solar_reduction_percent": (1-mean/req)*100,
            "sunward_thrust_fraction": average((case["sunward_required"] > 0).astype(float), t),
            "unlimited_sail_sunward_floor_mm_s2": average(np.maximum(case["sunward_required"], 0), t)*1000,
            "sunlit_fraction": average(case["visibility"], t),
            "distance_min_km": float(case["coords"][:, 0].min()/1000),
            "distance_max_km": float(case["coords"][:, 0].max()/1000),
            "transverse_max_km": float(length(case["coords"][:, 1:]).max()/1000),
            "earth_clearance_min_km": float(case["earth_clearance"].min()/1000)}


def closure(case, samples, base_sigma=0.05, mode="ideal", extra_sigma=0.0):
    """Re-solve optical authority as power hardware and propellant add mass."""
    prop = SCENARIO["propulsion"]
    ve, eta, kappa = prop["exhaust_velocity_m_s"], prop["efficiency"], prop["specific_power_W_kg"]
    cant = np.cos(np.deg2rad(prop["cant_deg"]))
    buffer = prop["propellant_buffer_days"]*K.JULIAN_DAY
    sigma = base_sigma+extra_sigma
    for count in range(100):
        summary = summarize(case, samples, sigma, mode)
        mean = summary["mean_residual_mm_s2"]/1000
        peak = summary["peak_residual_mm_s2"]/1000*prop["peak_margin_factor"]
        mdot_area = sigma*mean/(ve*cant)
        power_area = mdot_area*ve*ve/(2*eta)
        peak_area = sigma*peak*ve/(2*eta*cant)
        new_sigma = base_sigma+extra_sigma+peak_area/kappa+mdot_area*buffer
        if new_sigma > 1e5:
            return {"closed": False, "reason": "power and buffer mass diverge"}
        if abs(new_sigma-sigma) < 1e-11:
            break
        sigma = new_sigma
    else:
        return {"closed": False, "reason": "mass iteration did not converge"}
    area = summary["aperture_area_m2"]
    return {"closed": True, "iterations": count+1, "base_areal_mass_kg_m2": base_sigma,
            "extra_payload_kg_m2": extra_sigma, "total_areal_mass_kg_m2": sigma,
            "total_mass_kg": area*sigma, "mean_power_TW": area*power_area/1e12,
            "design_peak_power_TW": area*peak_area/1e12,
            "propellant_kg_s": area*mdot_area, "base_optical_mass_kg": area*base_sigma,
            "power_hardware_mass_kg": area*peak_area/kappa,
            "buffer_mass_kg": area*mdot_area*buffer,
            "gyr_propellant_kg": area*mdot_area*K.JULIAN_DAY*K.JULIAN_YEAR_DAYS*1e9,
            "summary": summary, "storage_included": False}


def search(samples, geo, off_axis_km=0, sigma=0.05, initial_coefficients=None):
    """Bounded three-harmonic search; no claim of global optimality."""
    n = geo["basis"].shape[1]
    initial = constant_coeff(78000) if initial_coefficients is None else np.asarray(initial_coefficients)
    active = 1 if off_axis_km == 0 else 3
    x0 = initial[:active].ravel()/1000
    lo, hi = np.asarray(SCENARIO["distance_bounds_km"])/1000
    bounds = [(lo, hi)]+[(-80, 80)]*(n-1)
    if active == 3:
        bounds += [(-off_axis_km/1000, off_axis_km/1000)]*(2*n)
    def unpack(x):
        c = np.zeros_like(initial); c[:active] = np.asarray(x).reshape(active, n)*1000
        return c
    def constraints(x):
        coords = geo["basis"]@unpack(x).T
        values = [coords[:, 0].min()/1000-lo, hi-coords[:, 0].max()/1000]
        if active == 3:
            values.append((off_axis_km-length(coords[:, 1:]).max())/1000)
        return np.array(values)
    reference_case = evaluate(samples, geo, initial, sigma=sigma)
    baseline = closure(reference_case, samples, base_sigma=sigma)["mean_power_TW"]
    def objective(x):
        coords = geo["basis"]@unpack(x).T
        if coords[:, 0].min() <= 0:
            return 1e5+abs(coords[:, 0].min())
        case = evaluate(samples, geo, unpack(x), sigma=sigma)
        value = closure(case, samples, base_sigma=sigma)
        return value["mean_power_TW"]/baseline if value["closed"] else 1e6
    # Independent phase profiles distinguish local convergence from a single start.
    starts = [x0.copy(), x0.copy()]
    starts[1][0] = 95; starts[1][1] = 10
    if active == 3:
        # Off-axis zero is a cusp in the aperture cost. Start inside the allowed
        # region as well as from the feasible, optimized on-axis trajectory.
        starts[1][n+1] = 0.1
    candidates = []
    for start in starts:
        result = minimize(objective, start, method="SLSQP", bounds=bounds,
                          constraints={"type": "ineq", "fun": constraints},
                          options={"maxiter": 140, "ftol": 2e-8, "eps": 1e-4})
        if constraints(result.x).min() < -1e-5:
            raise RuntimeError("Trajectory search violates its sampled distance/offset bounds")
        candidates.append(result)
    best = min(candidates, key=lambda x: x.fun)
    c = unpack(best.x)
    case = evaluate(samples, geo, c, sigma=sigma)
    return {"coefficients_km": c.tolist(), "off_axis_bound_km": off_axis_km,
            "success": bool(best.success), "message": str(best.message),
            "iterations": int(best.nit), "function_evaluations": int(best.nfev),
            "starts": [{"success": bool(r.success), "objective_ratio": float(r.fun),
                        "constraint_margin": float(constraints(r.x).min())} for r in candidates],
            "closure": closure(case, samples, base_sigma=sigma)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kernel", type=Path, default=DEFAULT_KERNEL)
    parser.add_argument("--out", type=Path, default=ROOT/"research/runs/solar_shield_array")
    parser.add_argument("--quick", action="store_true", help="pilot, without optimisation or nodal-span run")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    ephem = Ephemeris(args.kernel, SCENARIO["epoch_tdb"])
    start = time.monotonic()
    days = SCENARIO["primary_days"]
    t = np.linspace(0, days*K.JULIAN_DAY, int(days*24)+1)
    samples = ephem.sample(t); geo = geometry(samples, SCENARIO["harmonics"])
    baseline = evaluate(samples, geo, constant_coeff(SCENARIO["baseline_distance_km"]))
    result = {"schema": "terluna.protection.photogravitational-holding/1",
              "evidence": "DE440s inverse dynamics and optimistic optical-control bounds; finite-Sun Earth eclipses; point-mass Sun, Earth, Moon and planetary systems. Centre-force scans and area-quadrature final cases. A conservative first-order ray-flight allowance protects the moving Moon. Trajectory search is bounded, finite-dimensional and local.",
              "producer": {"runner": "research/studies/solar_shield_array/run.py",
                           "source_hashes": {}}, "ephemeris": ephem.manifest,
              "scenario": SCENARIO, "epoch_scale": "TDB", "primary_samples": len(t),
              "primary_step_s": float(t[1]-t[0]), "core_diverted_fraction": CORE_DIVERTED,
              "moon_point_mass_residual_max_m_s2": float(length(moon_acceleration_residual(samples)).max()),
              "baseline": {}, "fixed_distance_scan": [], "mass_sensitivity": []}
    for mode in ("none", "ideal", "absorb", "reflect"):
        result["baseline"][mode] = closure(baseline, samples, mode=mode)
    for sigma in (0.005, 0.01, 0.0259, 0.05, 0.1, 0.2, 1, 10):
        result["mass_sensitivity"].append({"sigma_kg_m2": sigma,
            "force_bound": summarize(baseline, samples, sigma),
            "closure": closure(baseline, samples, base_sigma=sigma)})
    for distance in np.arange(30000, 180001, 2000):
        case = evaluate(samples, geo, constant_coeff(float(distance)))
        result["fixed_distance_scan"].append({"distance_km": float(distance),
             "closure": closure(case, samples)})
    print("fixed-distance and mass scans complete", flush=True)
    if not args.quick:
        coarse_t = np.linspace(0, days*K.JULIAN_DAY, int(days*8)+1)
        coarse = ephem.sample(coarse_t); cg = geometry(coarse)
        result["trajectory_search"] = {}
        for label, offset in (("variable_distance", 0), ("variable_distance_offset", 2000)):
            print("search", label, flush=True)
            seed = None if offset == 0 else result["trajectory_search"]["variable_distance"]["coefficients_km"]
            result["trajectory_search"][label] = search(coarse, cg, offset, initial_coefficients=seed)
            c = np.array(result["trajectory_search"][label]["coefficients_km"])
            fine_case = evaluate(samples, geo, c)
            result["trajectory_search"][label]["fine_closure"] = closure(fine_case, samples)
            result["trajectory_search"][label]["fine_summary"] = summarize(fine_case, samples, 0.05)
        print("nodal-span screen", flush=True)
        years = SCENARIO["nodal_span_years"]
        nt = np.linspace(0, years*K.JULIAN_YEAR_DAYS*K.JULIAN_DAY, int(years*K.JULIAN_YEAR_DAYS*4)+1)
        ns = ephem.sample(nt); ng = geometry(ns)
        result["nodal_span"] = {"years": years, "samples": len(nt), "step_s": float(nt[1]-nt[0]), "cases": {}}
        trajectories = {"baseline": constant_coeff(78000),
                        **{k: np.array(v["coefficients_km"]) for k,v in result["trajectory_search"].items()}}
        for label, c in trajectories.items():
            case = evaluate(ns, ng, c)
            result["nodal_span"]["cases"][label] = closure(case, ns)
        print("spatial tile quadrature", flush=True)
        from .aperture import evaluate_aperture
        result["aperture_quadrature"] = {}
        for label in ("baseline", "variable_distance", "variable_distance_offset"):
            c = trajectories[label]
            result["aperture_quadrature"][label] = evaluate_aperture(coarse, cg, c)
        result["aperture_quadrature"]["baseline_without_solar"] = evaluate_aperture(coarse, cg, trajectories["baseline"], mode="none")
        result["aperture_quadrature"]["variable_distance_finer"] = evaluate_aperture(coarse, cg, trajectories["variable_distance"], rings=8, angles=24)
        result["aperture_quadrature"]["variable_distance_10g"] = evaluate_aperture(coarse, cg, trajectories["variable_distance"], base_sigma=0.01)
    np.savez_compressed(args.out/"annual_trace.npz", t=t, phase=geo["phase"],
                        required=baseline["required"], sunward=baseline["sunward_required"],
                        solar_visibility=baseline["visibility"])
    for relative in ("protection/dynamics/ephemeris.py", "protection/dynamics/model.py",
                     "protection/dynamics/optical.py", "research/studies/solar_shield_array/run.py", "research/studies/solar_shield_array/aperture.py", "research/studies/solar_shield_array/scenario.json"):
        result["producer"]["source_hashes"][relative] = hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()
    result["producer"]["constants"] = constants_used(result["producer"]["source_hashes"])
    result["runtime_s"] = time.monotonic()-start
    (args.out/"holding.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    ephem.close()
    print("wrote", args.out/"holding.json", "seconds", result["runtime_s"], flush=True)


if __name__ == "__main__":
    main()
