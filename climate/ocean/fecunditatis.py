"""Mare Fecunditatis's monthly tide through its straits, in MITgcm and SCHISM.

geography's linear check found Fecunditatis joined to the Near-side Sea by
narrow, shallow passages that set a 9.5-day exchange, close to the 9.61-day
tidal line, so that bed friction decides its tide (geography/README.md, The
monthly tide). Here both models run the whole Near-side Sea under the Moon's
tide from the ephemeris, with the passages at 16 pixels per degree and
quadratic bed friction:

- MITgcm on a latitude-longitude grid of 1/16 degree over Fecunditatis, its
  western channels and its passage to Crisium (33-66 E, 15 S-13 N), widening
  to 1/2 degree elsewhere;
- SCHISM on triangles over the same 16 px/deg terrain, 1/16 degree in the
  same box, widening to 1/2 degree beyond.

The tide is geography's degree-2 and degree-3 forcing with the seas'
equilibrium self-attraction (forcing.py), from 2026-01-01, ramped over five
days; the first 30 days spin the seas up and the next seven anomalistic
months (192.9 days) are analysed. Each chunk's snapshots (every three hours)
are reduced as it finishes and then deleted: basin-mean levels and volumes,
each anomalistic month's range at every water point, the largest speed and
the mean cube of speed (for the bed's dissipation), and the zones' largest
speeds through time. The Fecunditatis basin is the water joined to its deep
core below 160 m, with shallower water given to the nearest core, so its
boundary runs through the sills and its volume's change is the exchange
through the straits.

`python -m climate.ocean.fecunditatis grids` builds and describes the grids;
`... run mitgcm|schism [--drag 0.0025]` runs or resumes a case (`--smoke`: one
day in its own folder, timed, with its output through the reducer);
`... analyse` writes climate/results/ocean/fecunditatis.json.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from climate.ocean import forcing, lunar, mesh as meshes, mitgcm, schism, terrain

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = ROOT / "climate" / "results" / "ocean"
FINE = dict(lon=(33.0, 66.0), lat=(-15.0, 13.0))
DOMAIN = dict(lon=(-88.0, 76.0), lat=(-35.0, 74.0))
ZONES = {"western channels": dict(lon=(35.0, 44.5), lat=(-11.5, -0.5)),
         "Crisium passage": dict(lon=(53.0, 62.0), lat=(6.0, 11.5))}
OUTER_DEG = 0.5               # cells away from the fine box: half a degree, as SCHISM's coarsest triangles
CORE_DEPTH_M = 160.0          # deeper than both of Fecunditatis's sills (to the main sea about 110 m, to Crisium about 150 m)
SPINUP_DAYS = 30.0
MONTH_DAYS = 27.554550        # anomalistic month
MONTHS = 7
RAMP_DAYS = 5.0
FORCING_STEP_S = 10800.0
OUTPUT_S = 10800.0
CHUNK_DAYS = 10.0
DT_S = 300.0


def end_s():
    return (SPINUP_DAYS + MONTHS * MONTH_DAYS) * 86400.0


# ------------------------------------------------- the sea and its basins

def sea():
    """The 16 px/deg terrain and the Near-side Sea's connected water."""
    t = terrain.load()
    mask = terrain.connected_water(t["height"], t["level"], t["lat"], t["lon"], seeds=("Mare Imbrium", "Mare Fecunditatis"))
    return t, mask


def basins(t, mask, core_depth=CORE_DEPTH_M):
    """Basin labels on the terrain's pixels: deep cores below `core_depth`, every other water pixel given to the
    nearest core. Returns the labels and Fecunditatis's label."""
    from scipy import ndimage
    depth = t["level"] - t["height"]
    cores, _ = terrain.components(mask & (depth > core_depth))
    _, (ri, ci) = ndimage.distance_transform_edt(cores == 0, return_indices=True)
    labels = np.where(mask, cores[ri, ci], 0)
    la, lo = terrain.SEEDS["Mare Fecunditatis"]
    j = int(np.argmin(np.abs(t["lat"] - la)))
    i = int(np.argmin(np.abs(t["lon"] - lo)))
    return labels, int(labels[j, i])


def lookup(t, field, lat, lon):
    """Values of a pixel field at points (nearest pixel)."""
    j = np.clip(np.rint((np.asarray(lat) - t["lat"][0]) / (t["lat"][1] - t["lat"][0])).astype(int), 0, len(t["lat"]) - 1)
    i = np.clip(np.rint((np.asarray(lon) - t["lon"][0]) / (t["lon"][1] - t["lon"][0])).astype(int), 0, len(t["lon"]) - 1)
    return field[j, i]


# ------------------------------------------------- grids

def mitgcm_grid():
    t, mask = sea()
    xg = mitgcm.stretched_edges(*DOMAIN["lon"], FINE["lon"], 1 / 16, OUTER_DEG)
    yg = mitgcm.stretched_edges(*DOMAIN["lat"], FINE["lat"], 1 / 16, OUTER_DEG)
    if (len(xg) - 1) % 2:
        xg = np.r_[xg, xg[-1] + OUTER_DEG]
    if (len(yg) - 1) % 2:
        yg = np.r_[yg, yg[-1] + OUTER_DEG]
    depth = mitgcm.aggregate(t, mask, xg, yg)
    grid = mitgcm.Grid(xg=xg, yg=yg, depth=depth, name="Near-side Sea, 1/16 degree round Fecunditatis")
    seeds = []
    for name in ("Mare Imbrium", "Mare Fecunditatis"):
        la, lo = terrain.SEEDS[name]
        seeds.append((int(np.searchsorted(yg, la) - 1), int(np.searchsorted(xg, lo) - 1)))
    grid.depth = mitgcm.keep_connected(grid.depth, seeds)
    return grid


def schism_mesh(fine_level=3, outer_level=0):
    t, mask = sea()
    ny, nx = mask.shape
    by, bx = ny // meshes.BASE, nx // meshes.BASE
    lat_b = t["lat"][: by * meshes.BASE].reshape(by, meshes.BASE).mean(axis=1)
    lon_b = t["lon"][: bx * meshes.BASE].reshape(bx, meshes.BASE).mean(axis=1)
    inside = ((lat_b > FINE["lat"][0]) & (lat_b < FINE["lat"][1]))[:, None] & \
             ((lon_b > FINE["lon"][0]) & (lon_b < FINE["lon"][1]))[None, :]
    level = np.where(inside, fine_level, outer_level)
    d = t["level"] - t["height"]
    return meshes.build(d[: by * meshes.BASE, : bx * meshes.BASE], mask[: by * meshes.BASE, : bx * meshes.BASE],
                        level, t["lat"][: by * meshes.BASE], t["lon"][: bx * meshes.BASE])


def nodal_areas(m: meshes.Mesh):
    """Each node's share of the elements round it (a third of each), m^2."""
    areas = m.element_areas(lunar.RADIUS)
    out = np.zeros(m.nodes)
    for k in range(3):
        np.add.at(out, m.triangles[:, k], areas / 3.0)
    return out


# ------------------------------------------------- forcing

def forcing_series():
    seconds = np.arange(0.0, end_s() + FORCING_STEP_S, FORCING_STEP_S)
    return seconds, forcing.ephemeris_series(seconds) * forcing.ramp(seconds, RAMP_DAYS)[:, None]


# ------------------------------------------------- reductions

class Reducer:
    """Accumulates the analysis quantities from snapshots at water points (cells or nodes), saved after each chunk."""

    def __init__(self, path, lat, lon, area, basin, fec_label):
        self.path = Path(path)
        self.lat, self.lon, self.area = lat, lon, area
        self.fec = basin == fec_label
        self.zones = {name: (lon >= z["lon"][0]) & (lon <= z["lon"][1]) & (lat >= z["lat"][0]) & (lat <= z["lat"][1])
                      for name, z in ZONES.items()}
        n = len(lat)
        if self.path.exists():
            with np.load(self.path, allow_pickle=False) as s:
                self.state = {k: s[k] for k in s.files}
        else:
            self.state = dict(t=np.zeros(0), fec_mean=np.zeros(0), fec_volume=np.zeros(0), sea_mean=np.zeros(0),
                              month_max=np.full((MONTHS, n), -np.inf, np.float32),
                              month_min=np.full((MONTHS, n), np.inf, np.float32),
                              speed_max=np.zeros(n, np.float32), speed_cube_sum=np.zeros(n), count=np.zeros(1),
                              **{f"zone_speed_{k}": np.zeros(0) for k in range(len(ZONES))})

    def add(self, t_s, eta, speed):
        s = self.state
        if s["t"].size and t_s <= s["t"][-1] + 1:
            return                                      # already counted (a resumed chunk)
        s["t"] = np.r_[s["t"], t_s]
        s["fec_mean"] = np.r_[s["fec_mean"], np.sum(eta[self.fec] * self.area[self.fec]) / np.sum(self.area[self.fec])]
        s["fec_volume"] = np.r_[s["fec_volume"], np.sum(eta[self.fec] * self.area[self.fec])]
        s["sea_mean"] = np.r_[s["sea_mean"], np.sum(eta * self.area) / np.sum(self.area)]
        month = int((t_s / 86400.0 - SPINUP_DAYS) // MONTH_DAYS)
        if 0 <= month < MONTHS:
            s["month_max"][month] = np.maximum(s["month_max"][month], eta)
            s["month_min"][month] = np.minimum(s["month_min"][month], eta)
        s["speed_max"] = np.maximum(s["speed_max"], speed)
        s["speed_cube_sum"] = s["speed_cube_sum"] + speed.astype(float) ** 3
        s["count"] = s["count"] + 1
        for k, zone in enumerate(self.zones.values()):
            s[f"zone_speed_{k}"] = np.r_[s[f"zone_speed_{k}"], speed[zone].max() if zone.any() else 0.0]

    def save(self):
        tmp = self.path.with_suffix(".tmp.npz")
        np.savez_compressed(tmp, **self.state)
        tmp.replace(self.path)


def mitgcm_case(drag=lunar.DRAG, smoke=False):
    """The MITgcm case; `smoke` is one day in its own folder with snapshots from the start, to time the run and
    pass real output through the reducer."""
    grid = mitgcm_grid()
    name = "fecunditatis_smoke" if smoke else f"fecunditatis_cd{drag:g}"
    spec = mitgcm.tide_forcing_grid(grid, step=1.0)
    exp = mitgcm.Experiment(name=name, grid=grid, deltat=DT_S, nonlinear=True, nonlinear_surface=True, drag=drag, rotation=True,
                            tide=dict(spec, period_s=FORCING_STEP_S, file="tide.bin"), output_s=OUTPUT_S,
                            output_fields=("ETAN", "UVEL", "VVEL"), output_from_s=0.0 if smoke else SPINUP_DAYS * 86400.0,
                            implicit=0.6, viscosity_smag=2.2, tiles=(2, 2))
    return exp, spec


def run_mitgcm(drag=lunar.DRAG, until_s=None, smoke=False):
    exp, spec = mitgcm_case(drag, smoke)
    t, mask = sea()
    labels, fec = basins(t, mask)
    g = exp.grid
    if smoke or not (exp.folder / "tide.bin").exists():
        exp.write(fresh=smoke)
        seconds, coefs = forcing_series()
        forcing.write_mitgcm(exp.folder / "tide.bin.tmp", spec, seconds, coefs)
        (exp.folder / "tide.bin.tmp").replace(exp.folder / "tide.bin")
    wet = g.wet
    yc, xc = np.meshgrid(g.yc, g.xc, indexing="ij")
    reducer = Reducer(exp.folder / "reduced.npz", yc[wet], xc[wet], g.area[wet], lookup(t, labels, yc[wet], xc[wet]), fec)

    def after_chunk(e, first_it, last_it):
        for it in e.snapshots(first_it + 1, last_it):
            s = e.state(it)
            u = s["UVEL"]
            v = s["VVEL"]
            uc = 0.5 * (u + np.roll(u, -1, axis=1))     # centre values from the faces (land faces hold zero)
            vc = 0.5 * (v + np.roll(v, -1, axis=0))
            reducer.add(it * e.deltat, s["ETAN"][wet], np.hypot(uc, vc)[wet])
        reducer.save()
        for it in e.snapshots(first_it + 1, last_it):
            for suffix in ("data", "meta"):
                (e.folder / f"state.{it:010d}.{suffix}").unlink(missing_ok=True)

    exp.binary()                                        # built before the clock starts
    timed(exp.run, lambda: max([0] + exp.pickups()) * exp.deltat, 86400.0 if smoke else until_s or end_s(),
          CHUNK_DAYS * 86400.0, exp.deltat, after_chunk=after_chunk)
    return exp


def timed(run, progress_s, until_s, chunk_s, deltat, **kw):
    """Run to `until_s`, then print the model days and steps advanced, the wall time and the steps a minute
    (`progress_s()` gives the model time reached)."""
    import time
    begin, start = progress_s(), time.time()
    run(until_s, chunk_s, **kw)
    minutes, steps = (time.time() - start) / 60.0, (progress_s() - begin) / deltat
    print(f"day {begin / 86400:.2f} to {progress_s() / 86400:.2f}: {steps:.0f} steps in {minutes:.1f} min "
          f"({steps / max(minutes, 1e-9):.0f} steps a minute)")


def schism_case(drag=lunar.DRAG, smoke=False):
    m = schism_mesh()
    return schism.Experiment(name="fecunditatis_smoke" if smoke else f"fecunditatis_cd{drag:g}", mesh=m, dt=DT_S,
                             drag=drag, advection=True, rotation=True, output_s=OUTPUT_S,
                             file_s=(86400.0 if smoke else CHUNK_DAYS * 86400.0), compute_ranks=3)


def run_schism(drag=lunar.DRAG, until_s=None, smoke=False):
    exp = schism_case(drag, smoke)
    t, mask = sea()
    labels, fec = basins(t, mask)
    m = exp.mesh
    if smoke or not (exp.folder / "lunar_tide.th").exists():
        seconds, coefs = forcing_series()
        exp.write(seconds, coefs, fresh=smoke)
    reducer = Reducer(exp.folder / "reduced.npz", m.lat, m.lon, nodal_areas(m), lookup(t, labels, m.lat, m.lon), fec)
    first_s = 0.0 if smoke else SPINUP_DAYS * 86400.0

    def after_chunk(e, out, start_s, stop_s):
        import netCDF4
        for path in sorted(out.glob("out2d_*.nc"), key=lambda p: int(p.stem.split("_")[1])):
            with netCDF4.Dataset(path) as nc:
                times = np.array(nc["time"][:], float)
                eta = np.array(nc["elevation"][:], float)
                u = np.array(nc["depthAverageVelX"][:], float)
                v = np.array(nc["depthAverageVelY"][:], float)
            for k, ts in enumerate(times):
                if ts >= first_s - 1:
                    reducer.add(ts, eta[k], np.hypot(u[k], v[k]))
            path.unlink()
        reducer.save()

    timed(exp.run, lambda: exp.chunks_done()[-1]["end_s"] if exp.chunks_done() else 0.0,
          86400.0 if smoke else until_s or end_s(), CHUNK_DAYS * 86400.0, exp.dt, after_chunk=after_chunk)
    return exp


def describe_grids():
    t, mask = sea()
    labels, fec = basins(t, mask)
    px_area = terrain.cell_area(t["lat"], t["lon"])
    fec_area = float(px_area[labels == fec].sum())
    g = mitgcm_grid()
    m = schism_mesh()
    return dict(fecunditatis_basin_km2=fec_area / 1e6,
                mitgcm=dict(nx=g.nx, ny=g.ny, wet_cells=int(g.wet.sum()), water_km2=float(g.area[g.wet].sum() / 1e6)),
                schism=dict(nodes=m.nodes, elements=m.elements,
                            water_km2=float(m.element_areas(lunar.RADIUS).sum() / 1e6)),
                terrain_water_km2=float(px_area[mask & (t["level"] - t["height"] > 0)].sum() / 1e6))


# ------------------------------------------------- analysis

BANDS = {"monthly (25-35 d)": (25.0, 35.0), "fortnightly (13-16 d)": (13.0, 16.0), "nine-day (8.5-10.5 d)": (8.5, 10.5)}


def equilibrium_basin_mean(seconds, fec_label, t=None, mask=None):
    """Geography's equilibrium tide (each sea's volume held, self-attraction included) averaged over the
    Fecunditatis basin's atlas nodes, at model times."""
    from geography import tides
    if t is None:
        t, mask = sea()
    labels, _ = basins(t, mask)
    with np.load(ROOT / "geography/products/tides_28pct_4ppd.npz", allow_pickle=False) as product:
        lat, lon = product["lat_deg"][product["rows"]], product["lon_deg"][product["cols"]]
        basis, body = product["basis"].astype(float), product["labels"]
    lon = np.where(lon > 180, lon - 360, lon)
    inside = (lat > t["lat"][0]) & (lat < t["lat"][-1]) & (lon > t["lon"][0]) & (lon < t["lon"][-1])
    pick = np.zeros(len(lat), bool)
    pick[inside] = (body[inside] == 433) & (lookup(t, labels, lat[inside], lon[inside]) == fec_label)
    w = np.cos(np.radians(lat[pick]))
    mean_basis = basis[:, pick] @ w / w.sum()
    return forcing.ephemeris_series(seconds) @ mean_basis, int(pick.sum())


def band_gain(t_s, model, reference, band):
    """Complex gain of `model` against `reference` within a band of periods (days): sum M E* / sum |E|^2."""
    dt = t_s[1] - t_s[0]
    m = np.fft.rfft(model - model.mean())
    e = np.fft.rfft(reference - reference.mean())
    f = np.fft.rfftfreq(len(t_s), dt / 86400.0)
    sel = (f >= 1 / band[1]) & (f <= 1 / band[0])
    g = np.sum(m[sel] * np.conj(e[sel])) / np.sum(np.abs(e[sel]) ** 2)
    return dict(gain=float(abs(g)), lag_deg=float(-np.degrees(np.angle(g))),
                reference_rms_m=float(np.sqrt(2 * np.sum(np.abs(e[sel]) ** 2)) / len(t_s)))


def analyse_case(folder, area_of, fec_area_m2, eq_mean, eq_t):
    with np.load(Path(folder) / "reduced.npz", allow_pickle=False) as s:
        r = {k: s[k] for k in s.files}
    t_s = r["t"]
    if t_s.size < 10:
        return None
    eq = np.interp(t_s, eq_t, eq_mean)
    months = []
    for k in range(MONTHS):
        a = (SPINUP_DAYS + k * MONTH_DAYS) * 86400.0
        sel = (t_s >= a) & (t_s < a + MONTH_DAYS * 86400.0)
        if sel.sum() < 100:
            continue
        months.append(dict(month=k + 1, model_range_m=float(np.ptp(r["fec_mean"][sel])), equilibrium_range_m=float(np.ptp(eq[sel]))))
    q = np.gradient(r["fec_volume"], t_s)
    q_eq = np.gradient(eq * fec_area_m2, t_s)
    rho, n = lunar.WATER_DENSITY, float(r["count"][0])
    speed_cube_mean = r["speed_cube_sum"] / max(n, 1)
    zones = {}
    for k, name in enumerate(ZONES):
        z = r[f"zone_speed_{k}"]
        daily = z[: (len(z) // 8) * 8].reshape(-1, 8).max(axis=1)
        zones[name] = dict(largest_speed_m_s=float(z.max()), median_daily_largest_m_s=float(np.median(daily)),
                           p95_daily_largest_m_s=float(np.percentile(daily, 95)))
    return dict(snapshots=int(n), span_days=float((t_s[-1] - t_s[0]) / 86400.0),
                months=months,
                median_month_range=dict(model_m=float(np.median([m["model_range_m"] for m in months])),
                                        equilibrium_m=float(np.median([m["equilibrium_range_m"] for m in months]))),
                bands={name: band_gain(t_s, r["fec_mean"], eq, band) for name, band in BANDS.items()},
                exchange=dict(largest_m3_s=float(np.abs(q).max()), equilibrium_largest_m3_s=float(np.abs(q_eq).max()),
                              rms_m3_s=float(np.sqrt(np.mean(q ** 2))), equilibrium_rms_m3_s=float(np.sqrt(np.mean(q_eq ** 2)))),
                sea_mean_level_spread_m=float(np.ptp(r["sea_mean"])),
                zones=zones, speed_cube_mean=speed_cube_mean, rho=rho)


def analyse():
    t, mask = sea()
    labels, fec = basins(t, mask)
    px_area = terrain.cell_area(t["lat"], t["lon"])
    fec_area = float(px_area[labels == fec].sum())
    eq_t = np.arange(0.0, end_s() + 3600.0, 3600.0)
    eq_mean, eq_nodes = equilibrium_basin_mean(eq_t, fec, t, mask)
    record = dict(schema="terluna.climate.ocean-fecunditatis/1", fecunditatis_basin_km2=fec_area / 1e6,
                  equilibrium_nodes=eq_nodes, cases={})
    for model in ("mitgcm", "schism"):
        for folder in sorted((HERE / "runs" / model).glob("fecunditatis_cd*")):
            drag = float(folder.name.split("_cd")[1])
            if model == "mitgcm":
                grid = mitgcm_case(drag)[0].grid
                wet = grid.wet
                yc, xc = np.meshgrid(grid.yc, grid.xc, indexing="ij")
                lat, lon, area = yc[wet], xc[wet], grid.area[wet]
            else:
                m = schism_case(drag).mesh
                lat, lon, area = m.lat, m.lon, nodal_areas(m)
            out = analyse_case(folder, area, fec_area, eq_mean, eq_t)
            if out is None:
                continue
            basin = lookup(t, labels, lat, lon)
            power = out.pop("speed_cube_mean") * out.pop("rho") * drag * area
            out["bed_dissipation_w"] = dict(near_side_sea=float(power.sum()), fecunditatis_basin=float(power[basin == fec].sum()),
                                            **{name: float(power[(lon >= z["lon"][0]) & (lon <= z["lon"][1]) &
                                                                 (lat >= z["lat"][0]) & (lat <= z["lat"][1])].sum())
                                               for name, z in ZONES.items()})
            out["drag"] = drag
            record["cases"][f"{model}_cd{drag:g}"] = out
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "fecunditatis.json").write_text(json.dumps(record, indent=1) + "\n")
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("action", choices=["grids", "run", "analyse"])
    parser.add_argument("model", nargs="?", choices=["mitgcm", "schism"])
    parser.add_argument("--drag", type=float, default=lunar.DRAG)
    parser.add_argument("--days", type=float, help="stop at this model day (a partial run)")
    parser.add_argument("--smoke", action="store_true",
                        help="one day in runs/<model>/fecunditatis_smoke with output from the start: speed and reducer")
    args = parser.parse_args()
    if args.action == "grids":
        print(json.dumps(describe_grids(), indent=1))
    elif args.action == "run":
        until = args.days * 86400.0 if args.days else None
        (run_mitgcm if args.model == "mitgcm" else run_schism)(args.drag, until, smoke=args.smoke)
    else:
        print(json.dumps(analyse(), indent=1)[:4000])


if __name__ == "__main__":
    main()
