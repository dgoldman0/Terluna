"""Drive a lunar sea with recovered atmospheric momentum-transfer histories."""
from __future__ import annotations

import argparse
from datetime import timedelta
import json
from pathlib import Path

import numpy as np

from climate.waves.model import ORIGIN, ROOT, sha256
from climate.waves.pilot import BasinCase, array_text, basin_files, execute, load_basin, read_basin
from shared.constants import DATA

HERE = Path(__file__).resolve().parent
RUNS = ROOT / "research/runs/waves/weather"
DRAG = DATA["model_closures"]["swan_41_51_wu"]


def drag_coefficient(speed):
    speed = np.asarray(speed)
    if not np.isfinite(speed).all() or np.any(speed < 0):
        raise ValueError("Drag requires finite, nonnegative speed")
    return np.where(speed <= DRAG["transition_speed_m_s"], DRAG["low_speed_drag_coefficient"],
                    DRAG["high_speed_intercept"] + DRAG["high_speed_slope_s_m"] * speed)


def stress_equivalent_wind(eastward, northward, density):
    """Invert rho Cd(U) U², preserving the signed stress direction, including calm."""
    eastward, northward = np.broadcast_arrays(eastward, northward)
    if density <= 0 or not np.isfinite(density) or not np.isfinite(eastward + northward).all():
        raise ValueError("Stress and density must be finite, with positive density")
    stress = np.hypot(eastward, northward)
    # The low-speed coefficient supplies an upper speed bound because Cd grows.
    low, high = np.zeros_like(stress), np.sqrt(stress / (density * DRAG["low_speed_drag_coefficient"]))
    for _ in range(55):
        middle = (low + high) * .5
        below = density * drag_coefficient(middle) * middle**2 < stress
        low, high = np.where(below, middle, low), np.where(below, high, middle)
    speed = (low + high) * .5
    scale = np.divide(speed, stress, out=np.zeros_like(speed), where=stress > 0)
    return eastward * scale, northward * scale


def ocean_weights(longitude, latitude, water, target_lon, target_lat):
    """Bilinear weights renormalized over ocean corners; nearest ocean fills gaps.

    The regional target must lie inside the native coordinate range. Water is
    a static mask of cells that remain open water throughout the saved record.
    """
    longitude, latitude = np.asarray(longitude), np.asarray(latitude)
    if np.any(np.diff(longitude) <= 0) or np.any(np.diff(latitude) <= 0):
        raise ValueError("Ocean interpolation requires ascending coordinates")
    x, y = np.broadcast_arrays(target_lon, target_lat)
    if np.any((x < longitude[0]) | (x > longitude[-1]) | (y < latitude[0]) | (y > latitude[-1])):
        raise ValueError("Target lies outside the atmospheric grid")
    if water.shape != (len(latitude), len(longitude)) or not water.any():
        raise ValueError("An open-water source mask is required")
    ix = np.clip(np.searchsorted(longitude, x) - 1, 0, len(longitude) - 2)
    iy = np.clip(np.searchsorted(latitude, y) - 1, 0, len(latitude) - 2)
    a = (x - longitude[ix]) / (longitude[ix+1] - longitude[ix])
    b = (y - latitude[iy]) / (latitude[iy+1] - latitude[iy])
    indices = np.stack((iy*len(longitude)+ix, iy*len(longitude)+ix+1,
                        (iy+1)*len(longitude)+ix, (iy+1)*len(longitude)+ix+1), axis=-1)
    weights = np.stack(((1-a)*(1-b), a*(1-b), (1-a)*b, a*b), axis=-1)
    weights *= water.ravel()[indices]
    available = weights.sum(axis=-1)
    fallback = available == 0
    weights /= np.where(fallback, 1, available)[..., None]
    rows, cols = np.where(water)
    for place in np.ndindex(fallback.shape):
        if not fallback[place]:
            continue
        # Great-circle dot product chooses the closest open-water cell.
        yr, xr = np.radians(y[place]), np.radians(x[place])
        yr0, xr0 = np.radians(latitude[rows]), np.radians(longitude[cols])
        score = np.sin(yr)*np.sin(yr0) + np.cos(yr)*np.cos(yr0)*np.cos(xr0-xr)
        nearest = np.argmax(score)
        indices[place] = rows[nearest]*len(longitude)+cols[nearest]
        weights[place] = [1, 0, 0, 0]
    return indices, weights, fallback


def apply_weights(field, indices, weights):
    field = np.asarray(field)
    return np.sum(field.reshape(*field.shape[:-2], -1)[..., indices] * weights, axis=-1)


def forcing(atmosphere, basin, density, hours=168):
    """Select a week containing the largest basin-mean stress and retain its dates."""
    with np.load(atmosphere, allow_pickle=False) as d:
        meta = json.loads(d["metadata"].item())
        if meta["schema"] != "terluna.climate.gcm-wave-snapshots/1":
            raise ValueError("Unsupported atmospheric product")
        times = d["time_s"]
        lat, lon = d["latitude_deg"][::-1], d["longitude_deg"]
        mask = np.all((d["land_fraction"][:, ::-1] < .5) & (d["sea_ice_fraction"][:, ::-1] < .01), axis=0)
        xx, yy = np.meshgrid(basin["lon"], basin["lat"])
        indices, weights, fallback = ocean_weights(lon, lat, mask, xx, yy)
        tx = apply_weights(d["stress_eastward_Pa"][:, ::-1], indices, weights)
        ty = apply_weights(d["stress_northward_Pa"][:, ::-1], indices, weights)
        rho = apply_weights(d["stress_closure_density_kg_m3"][:, ::-1], indices, weights)
    if not np.allclose(np.diff(times), 10800) or hours % 3 or hours < 48 or hours*3600 > times[-1]:
        raise ValueError("A complete three-hourly record and multi-day window are required")
    wet = basin["wet"]
    area_weight = np.cos(np.radians(yy[wet]))
    mean_stress = np.average(np.hypot(tx, ty)[:, wet], weights=area_weight, axis=1)
    peak_index = int(np.argmax(mean_stress))
    first = int(np.clip(peak_index - 24//3, 0, len(times) - 1 - hours//3))
    last = first + hours//3
    # Interpolate stress vectors every 15 minutes before converting to SWAN
    # input. This keeps momentum interpolation close between native snapshots.
    target_time = np.arange(0, hours*3600 + 1, 900)
    interval = target_time / 10800 + first
    left = np.minimum(interval.astype(int), last - 1)
    alpha = (interval - left)[:, None, None]
    targets = [(1-alpha)*f[left] + alpha*f[left+1] for f in (tx, ty)]
    u, v = stress_equivalent_wind(*targets, density)
    # Midpoint error quantifies SWAN's linear interpolation of these inputs.
    um, vm = .5*(u[:-1]+u[1:]), .5*(v[:-1]+v[1:])
    magnitude = np.hypot(um, vm)
    computed = [density * drag_coefficient(magnitude) * magnitude * f for f in (um, vm)]
    wanted = [.5*(f[:-1]+f[1:]) for f in targets]
    error = np.hypot(computed[0]-wanted[0], computed[1]-wanted[1])[:, wet]
    reference = np.hypot(*wanted)[:, wet]
    xnative, ynative = np.meshgrid(lon, lat)
    source_distance = np.arccos(np.clip(
        np.sin(np.radians(yy[..., None]))*np.sin(np.radians(ynative.ravel()[indices])) +
        np.cos(np.radians(yy[..., None]))*np.cos(np.radians(ynative.ravel()[indices])) *
        np.cos(np.radians(xx[..., None]-xnative.ravel()[indices])), -1, 1)) * meta["radius_m"]
    used = (weights > 0) & wet[..., None]
    record = dict(schema="terluna.climate.swan-weather-forcing/1", atmosphere_sha256=sha256(atmosphere),
                  atmosphere_metadata=meta, air_density_kg_m3=density,
                  coupling="rho_swan * Cd_Wu(Ueq) * |Ueq| * (ueq,veq) = interpolated GCM surface stress vector",
                  interpretation="Stress-equivalent SWAN input; atmospheric 10 m wind diagnostics remain a separate surface-layer calculation.",
                  time_interpolation="Linear surface stress at 15-minute knots, then SWAN linear vector interpolation between knots",
                  spatial_interpolation="Bilinear ocean-corner weights renormalized over cells remaining ice-free; nearest open-water cell fills an empty stencil",
                  selection=f"{hours}-hour window containing the maximum area-weighted basin-mean stress; starts 24 hours before that peak where the archive permits",
                  first_snapshot_index=first, last_snapshot_index=last, hours=hours,
                  start_day_from_first_snapshot=float(times[first]/86400),
                  end_day_from_first_snapshot=float(times[last]/86400),
                  peak_mean_stress_day=float(times[peak_index]/86400),
                  global_snapshot_count=len(times), regional_open_water_source_cells=int(len(np.unique(indices[used]))),
                  wet_nodes=int(wet.sum()), nearest_ocean_fallback_wet_nodes=int(fallback[wet].sum()),
                  farthest_contributing_ocean_cell_km=float(source_distance[used].max()/1000),
                  stress_interpolation_midpoint_absolute_max_Pa=float(error.max()),
                  stress_interpolation_midpoint_relative_p99=float(np.percentile((error/np.maximum(reference, 1e-4))[reference >= 1e-4], 99)),
                  equivalent_speed_percentiles_m_s=dict(zip(("p10", "p50", "p90", "p99", "max"),
                        np.percentile(np.hypot(u,v)[:,wet], [10,50,90,99,100]).tolist())),
                  stress_closure_density_kg_m3=dict(mean=float(rho[first:last+1,wet].mean()),
                                                    min=float(rho[first:last+1,wet].min()), max=float(rho[first:last+1,wet].max())),
                  monthly_time_days=(times/86400).tolist(), monthly_basin_mean_stress_Pa=mean_stress.tolist())
    return dict(u=u, v=v, time_s=target_time, record=record)


def files(case, basin, inputs):
    result = basin_files(case, basin, 1.)
    lon, lat = basin["lon"], basin["lat"]
    start = ORIGIN.strftime("%Y%m%d.%H%M%S")
    end = (ORIGIN + timedelta(hours=case.hours)).strftime("%Y%m%d.%H%M%S")
    old = result["INPUT"]
    begin, finish = old.index("INPGRID WIND"), old.index("GEN3 KOMEN")
    wind_grid = f"""INPGRID WIND {lon[0]:.12g} {lat[0]:.12g} 0 {len(lon)-1} {len(lat)-1} {lon[1]-lon[0]:.12g} {lat[1]-lat[0]:.12g} &
 NONSTAT {start} 15 MIN {end}
READINP WIND 1 'wind.dat' 3 0 0 0 FREE
"""
    result["INPUT"] = old[:begin] + wind_grid + old[finish:]
    # Extra output checks the physical components seen by SWAN against the file.
    result["INPUT"] = result["INPUT"].replace("COMPUTE NONSTAT", f"TABLE 'COMPGRID' NOHEAD 'wind.tbl' TSEC XP YP WIND UFRI CDRAG &\n OUTPUT {start} 3 HR\nCOMPUTE NONSTAT")
    result["wind.dat"] = "".join(array_text(field) for pair in zip(inputs["u"], inputs["v"]) for field in pair)
    return result


def run_case(executable, atmosphere, name, stride=4, step_s=150, hours=168):
    build = json.loads((executable.parent.parent / "coupled_air.json").read_text())
    if sha256(executable) != build["build"]["executable_sha256"]:
        raise ValueError("SWAN executable differs from its density/build record")
    basin = load_basin(stride=stride)
    inputs = forcing(atmosphere, basin, build["air_density_kg_m3"], hours)
    case = BasinCase(name, stride=stride, step_s=step_s, frequency_intervals=48, hours=hours)
    record = execute(executable, RUNS / name, files(case, basin, inputs),
                     ("basin.tbl", "offshore.spc", "wind.tbl"), timeout_s=3600)
    cube = read_basin(RUNS / name / "basin.tbl", case, basin)
    wind = np.loadtxt(RUNS / name / "wind.tbl").reshape(hours//3 + 1, *basin["wet"].shape, 7)
    if not np.allclose(wind[:, basin["wet"], 3:5],
                       np.stack((inputs["u"][::12, basin["wet"]], inputs["v"][::12, basin["wet"]]), axis=-1),
                       rtol=8e-5, atol=1e-5):
        raise ValueError("SWAN wind components differ from the supplied forcing")
    # The initial UFRI field precedes the first wind-source calculation.
    speed = np.hypot(wind[1:, basin["wet"], 3], wind[1:, basin["wet"], 4])
    expected_ustar = np.sqrt(drag_coefficient(speed)) * speed
    ustar_error = np.abs(wind[1:, basin["wet"], 5] - expected_ustar)
    if not np.allclose(wind[1:, basin["wet"], 5], expected_ustar, rtol=1e-4, atol=1e-6):
        raise ValueError("SWAN friction velocity differs from the stress coupling")
    result = dict(schema="terluna.climate.weather-wave-case/1", name=name, stride=stride,
                  step_s=step_s, frequency_intervals=48, hours=hours, run=record,
                  forcing=inputs["record"], atlas_sha256=basin["source_sha256"],
                  longitude_deg=basin["lon"].tolist(), latitude_deg=basin["lat"].tolist(),
                  offshore_index=basin["offshore_index"], wet=basin["wet"].tolist(),
                  columns=["time_s","longitude_deg","latitude_deg","hs_m","peak_period_s","mean_period_s","depth_m","direction_deg","breaking_fraction"],
                  stress_coupling_ustar_absolute_max_m_s=float(ustar_error.max()),
                  values=cube.tolist())
    path = RUNS / name / "product.json"
    path.write_text(json.dumps(result, separators=(",", ":")) + "\n")
    return path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--executable", type=Path, required=True)
    p.add_argument("--atmosphere", type=Path, default=RUNS / "atmosphere.npz")
    p.add_argument("--name", default="basin_dt150")
    p.add_argument("--step", type=int, default=150)
    p.add_argument("--stride", type=int, default=4)
    p.add_argument("--hours", type=int, default=168)
    a = p.parse_args()
    print(run_case(a.executable, a.atmosphere, a.name, a.stride, a.step, a.hours))


if __name__ == "__main__":
    main()
