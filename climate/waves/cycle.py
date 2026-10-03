"""Measure the second solar cycle after a full cycle of wave spin-up.

Hourly SWAN samples define piecewise-linear height histories. Threshold
crossings and the two cycle boundaries are integrated on those histories.
"""
from __future__ import annotations

import argparse
from datetime import timedelta
import json
from pathlib import Path

import numpy as np

from climate.waves.model import ORIGIN, ROOT, sha256
from climate.waves.pilot_checks import spectrum_diagnostics
from shared.constants import SYNODIC_MONTH_DAYS

RUNS = ROOT / "research/runs/waves/cycle"
HERE = Path(__file__).resolve().parent
CYCLE_HOURS = SYNODIC_MONTH_DAYS * 24


def window(times, values, start, end):
    """Clip a sampled history to exact endpoints, interpolating only its edges."""
    times, values = np.asarray(times, dtype=float), np.asarray(values, dtype=float)
    if (times.ndim != 1 or len(times) < 2 or values.shape[0] != len(times)
            or not np.isfinite(times).all() or not np.isfinite(values).all()
            or np.any(np.diff(times) <= 0)):
        raise ValueError("A finite history with strictly increasing times is required")
    if not times[0] <= start < end <= times[-1]:
        raise ValueError("The complete reporting interval must lie inside the history")
    selected = np.r_[start, times[(times > start) & (times < end)], end]
    left = np.clip(np.searchsorted(times, selected, side="right") - 1, 0, len(times)-2)
    fraction = (selected-times[left]) / (times[left+1]-times[left])
    fraction = fraction.reshape((-1,) + (1,)*(values.ndim-1))
    return selected, values[left]*(1-fraction) + values[left+1]*fraction


def exceedance(times, values, threshold):
    """Duration and contiguous episodes at/above a threshold on a linear history."""
    times, values = np.asarray(times, float), np.asarray(values, float)
    if (times.ndim != 1 or values.shape != times.shape or len(times) < 2
            or not np.isfinite(times).all() or not np.isfinite(values).all()
            or not np.isfinite(threshold) or np.any(np.diff(times) <= 0)):
        raise ValueError("Threshold statistics require finite, ordered scalar samples")
    intervals = []
    for a, b, ha, hb in zip(times[:-1], times[1:], values[:-1], values[1:]):
        if ha < threshold and hb < threshold:
            continue
        left, right = a, b
        if ha < threshold:
            left = a+(b-a)*(threshold-ha)/(hb-ha)
        if hb < threshold:
            right = a+(b-a)*(threshold-ha)/(hb-ha)
        if right <= left:
            continue
        if intervals and abs(left-intervals[-1][1]) <= 1e-9:
            intervals[-1][1] = right
        else:
            intervals.append([left, right])
    duration = sum(b-a for a, b in intervals)
    return dict(hours=float(duration), fraction=float(duration/(times[-1]-times[0])),
                longest_episode_hours=float(max((b-a for a,b in intervals), default=0)),
                episode_count=len(intervals), intervals_hours=intervals,
                begins_above=bool(intervals and abs(intervals[0][0]-times[0]) <= 1e-9),
                ends_above=bool(intervals and abs(intervals[-1][1]-times[-1]) <= 1e-9))


def read_case(path):
    result = json.loads(path.read_text())
    if result["schema"] != "terluna.climate.weather-wave-case/1":
        raise ValueError("Unsupported wave-history product")
    for name, expected in {**result["run"]["identity"]["inputs"], **result["run"]["output_sha256"]}.items():
        if sha256(path.parent/name) != expected:
            raise ValueError(f"Wave-history bytes changed: {name}")
    return result


def summarize(path, *, cycle_hours=CYCLE_HOURS):
    # Preserve the checkout path when the run directory links to external storage.
    path = path.absolute()
    case = read_case(path)
    if case["forcing"]["first_snapshot_index"] != 0:
        raise ValueError("A cycle study begins at the first atmospheric snapshot")
    atmosphere = case["forcing"]["atmosphere_metadata"]
    config = atmosphere["source_configuration"]
    if config["model"].get("sun_clock") != "synodic":
        raise ValueError("Cycle analysis requires the corrected synodic solar clock")
    planet = config["planet"]
    solar_days = 1/(1/planet["rotationperiod"] - 1/planet["year"])
    clock_difference_s = 2*(SYNODIC_MONTH_DAYS-solar_days)*86400
    if abs(clock_difference_s) > 1:
        raise ValueError("The atmospheric solar period differs from the shared lunar cycle")
    cube = np.asarray(case["values"])
    if cube.ndim != 3 or cube.shape[2] != 9:
        raise ValueError("Expected the nine SWAN basin columns")
    times = cube[:,0,0]/3600
    np.testing.assert_allclose(cube[...,0]-times[:,None]*3600, 0, atol=.1, rtol=0)
    if not np.array_equal(times, np.arange(case["hours"]+1)):
        raise ValueError("Cycle statistics require complete hourly output")
    if np.any(cube[...,3] < 0):
        raise ValueError("Invalid significant wave height")
    # Interpolate the direction on its unwrapped angular path at the endpoints.
    unwrapped = cube.copy()
    unwrapped[...,7] = np.degrees(np.unwrap(np.radians(cube[...,7]), axis=0))
    t, values = window(times, unwrapped, cycle_hours, 2*cycle_hours)
    values[...,7] %= 360
    hs, weights = values[...,3], np.cos(np.radians(values[0,:,2]))
    weights /= weights.sum()
    dt = np.diff(t)
    mean_hs = np.sum(dt[:,None]*(hs[:-1]+hs[1:])*.5, axis=0)/cycle_hours
    # Exact integral of the square of a linearly interpolated height.
    mean_hs2 = np.sum(dt[:,None]*(hs[:-1]**2+hs[:-1]*hs[1:]+hs[1:]**2)/3, axis=0)/cycle_hours
    peak_time, peak_node = np.unravel_index(np.argmax(hs), hs.shape)
    offshore = case["offshore_index"]
    thresholds = []
    for level in (.5, 1., 2., 3.):
        nodes = [exceedance(t, hs[:,n], level) for n in range(hs.shape[1])]
        thresholds.append(dict(hs_m=level,
            node_hours=[n["hours"] for n in nodes],
            node_fraction=[n["fraction"] for n in nodes],
            node_longest_episode_hours=[n["longest_episode_hours"] for n in nodes],
            node_episode_count=[n["episode_count"] for n in nodes],
            node_begins_above=[n["begins_above"] for n in nodes],
            node_ends_above=[n["ends_above"] for n in nodes],
            area_time_fraction=float(np.dot(weights,[n["fraction"] for n in nodes])),
            sampled_area_fraction=((hs >= level)*weights).sum(axis=1).tolist(),
            offshore=nodes[offshore]))
    # Energy-weighted propagation direction: a low resultant marks variable
    # directions. Endpoint trapezoidal quadrature is stated in the product.
    angle = np.radians(values[...,7])
    vector = hs**2 * np.exp(1j*angle)
    integrated = np.sum(dt[:,None]*(vector[:-1]+vector[1:])*.5, axis=0)
    energy_weight = np.sum(dt[:,None]*(hs[:-1]**2+hs[1:]**2)*.5, axis=0)
    direction = np.degrees(np.angle(integrated)) % 360
    strength = np.divide(np.abs(integrated), energy_weight,
                         out=np.zeros_like(energy_weight), where=energy_weight > 0)
    stamp = (ORIGIN+timedelta(hours=case["hours"])).strftime("%Y%m%d.%H%M%S")
    result = dict(schema="terluna.climate.wave-cycle/1",
        evidence="A conditional Smythii–Marginis wave record driven by one GCM continuation; the first solar cycle initializes the waves and the second supplies the statistics.",
        reading_rule="Wave statistics use the second solar cycle; the forcing section retains the full input history, including spin-up. Threshold durations use piecewise-linear hourly heights and exact cycle endpoints. Episodes touching a reporting boundary are censored. Spatial fractions use cosine-latitude weights at resolved wet nodes. The stated climate, coarse geography and wave physics define this conditional cycle; long-term occurrence and coastal surf remain further calculations.",
        producer=dict(domain="climate", files={str(p.relative_to(ROOT)):sha256(p) for p in
            (Path(__file__), HERE/"weather.py", HERE/"pilot.py", ROOT/"shared/constants.json")}),
        inputs=dict(case_path=str(path.relative_to(ROOT)), case_sha256=sha256(path),
                    atmosphere_sha256=case["forcing"]["atmosphere_sha256"],atlas_sha256=case["atlas_sha256"]),
        case=dict(name=case["name"], step_s=case["step_s"], frequency_intervals=case["frequency_intervals"],
                  stride=case["stride"], hours=case["hours"], run=case["run"]),
        clock=dict(cycle_hours=cycle_hours, cycle_days=cycle_hours/24,
                   basis="Shared synodic month; elapsed time is measured from the first atmospheric snapshot",
                   phase_origin="First saved snapshot; the reporting boundary is independent of local sunrise",
                   atmospheric_solar_period_days=solar_days,
                   two_cycle_clock_difference_s=clock_difference_s,
                   spinup_interval_hours=[0,cycle_hours], report_interval_hours=[cycle_hours,2*cycle_hours],
                   reported_hours=cycle_hours, interpolated_boundary_hours=[float(t[0]),float(t[-1])],
                   first_snapshot_restart_elapsed_s=atmosphere["first_snapshot_restart_step"]*config["model"]["timestep_min"]*60),
        grid=dict(longitude_deg=case["longitude_deg"],latitude_deg=case["latitude_deg"],wet=case["wet"],
                  node_longitude_deg=values[0,:,1].tolist(),node_latitude_deg=values[0,:,2].tolist(),
                  depth_m=values[0,:,6].tolist(),area_weights=weights.tolist(),offshore_index=offshore),
        units=dict(height="m",period="s",time="Earth hours since first atmospheric snapshot",direction="degrees counterclockwise from east, pointing toward propagation"),
        forcing=case["forcing"],
        statistics=dict(peak_record=values[peak_time,peak_node].tolist(),max_hs_m=float(hs.max()),
                        peak_day_in_reporting_cycle=float((t[peak_time]-cycle_hours)/24),
                        node_max_hs_m=hs.max(axis=0).tolist(),node_mean_hs_m=mean_hs.tolist(),
                        node_rms_hs_m=np.sqrt(mean_hs2).tolist(),
                        area_time_mean_hs_m=float(np.dot(weights,mean_hs)),
                        node_energy_weighted_direction_deg=direction.tolist(),
                        node_direction_resultant=strength.tolist(),
                        direction_quadrature="Trapezoidal integration of Hs squared times the unit propagation vector at hourly samples and interpolated endpoints",
                        thresholds=thresholds),
        series=dict(time_hours=t.tolist(),max_hs_m=hs.max(axis=1).tolist(),
                    area_mean_hs_m=(hs*weights).sum(axis=1).tolist(),
                    offshore_columns=case["columns"],offshore=values[:,offshore].tolist()),
        spectral_edge_check=spectrum_diagnostics(path.parent/"offshore.spc",stamp),
        limitations=["One reporting cycle establishes occurrence within that cycle; inter-cycle and seasonal variability remain open.",
                     "The 1-degree basin grid retains unresolved islands, shores and underwater slopes.",
                     "The earlier 150-to-75-second comparison passed for metre-scale waves and failed the full-field criterion in small waves; cycle-specific refinement is recorded separately.",
                     "Final offshore spectral coverage has a separate check; evolving spectra across the basin need further sampling.",
                     "Lunar wind-input and dissipation calibration, currents, changing water levels and beach run-up remain open."])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, default=RUNS/"basin_dt150/product.json")
    parser.add_argument("--output", type=Path, default=HERE/"results/cycle.json")
    args = parser.parse_args()
    result = summarize(args.case)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,separators=(",",":"),allow_nan=False)+"\n")
    print(json.dumps(dict(output=str(args.output), max_hs_m=result["statistics"]["max_hs_m"],
                         reporting_hours=result["clock"]["reported_hours"]),indent=2))


if __name__ == "__main__":
    main()
