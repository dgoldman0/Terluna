"""Figures of Mare Fecunditatis's tide from the ocean models' reduced runs.

Reads climate/ocean/runs/<model>/fecunditatis_cd<Cd>/reduced.npz (written by
climate.ocean.fecunditatis) with the run's own grid or mesh, and geography's
equilibrium tide, and draws into out/ (kept out of Git):

- level.png: the Fecunditatis basin's mean level through the analysed months, each model against the equilibrium;
- range.png: each model's median monthly range round Fecunditatis, and the equilibrium's, on one colour scale;
- speed.png: each model's largest depth-averaged speed over the analysed months round Fecunditatis.

The maps show the models' own water points (MITgcm's cells, SCHISM's triangles); nothing is smoothed.
`python visualization/ocean-tides/render.py [--drag 0.0025]`
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from climate.ocean import fecunditatis as study

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
VIEW = dict(lon=(30.0, 68.0), lat=(-18.0, 16.0))


def load(model, drag):
    folder = study.HERE / "runs" / model / f"fecunditatis_cd{drag:g}"
    path = folder / "reduced.npz"
    if not path.exists():
        return None
    with np.load(path, allow_pickle=False) as s:
        r = {k: s[k] for k in s.files}
    if model == "mitgcm":
        with np.load(folder / "grid.npz") as g:
            grid = study.mitgcm.Grid(xg=g["xg"], yg=g["yg"], depth=g["depth"])
        r["grid"] = grid
    else:
        with np.load(folder / "mesh.npz") as m:
            r["mesh"] = {k: m[k] for k in m.files}
    return r


def field_on_grid(grid, values):
    out = np.full(grid.depth.shape, np.nan)
    out[grid.wet] = values
    return out


def draw(ax, run, values, **kw):
    if "grid" in run:
        g = run["grid"]
        return ax.pcolormesh(g.xg, g.yg, field_on_grid(g, values), shading="flat", **kw)
    m = run["mesh"]
    import matplotlib.tri as mtri
    tri = mtri.Triangulation(m["lon"], m["lat"], m["triangles"])
    return ax.tripcolor(tri, values, shading="gouraud", **kw)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--drag", type=float, default=study.lunar.DRAG)
    args = parser.parse_args()
    OUT.mkdir(exist_ok=True)
    runs = {m: r for m in ("mitgcm", "schism") if (r := load(m, args.drag)) is not None}
    if not runs:
        raise SystemExit("no reduced runs")
    t, mask = study.sea()
    _, fec = study.basins(t, mask)
    eq_t = np.arange(0.0, study.end_s() + 3600.0, 3600.0)
    eq_mean, _ = study.equilibrium_basin_mean(eq_t, fec, t, mask)

    fig, ax = plt.subplots(figsize=(11, 4.2))
    days = eq_t / 86400.0
    keep = days >= study.SPINUP_DAYS
    ax.plot(days[keep], eq_mean[keep], color="0.45", lw=1.4, label="equilibrium (geography)")
    for (model, r), colour in zip(runs.items(), ("#1f6fb4", "#d1495b")):
        ax.plot(r["t"] / 86400.0, r["fec_mean"], color=colour, lw=1.1, label=f"{model}, Cd {args.drag:g}")
    ax.set_xlabel("days from 2026-01-01")
    ax.set_ylabel("Fecunditatis basin mean level (m)")
    ax.legend(frameon=False, ncol=3)
    ax.grid(alpha=0.3)
    fig.savefig(OUT / "level.png", dpi=110, bbox_inches="tight")
    plt.close(fig)

    with np.load(study.ROOT / "geography/products/tides_28pct_4ppd.npz", allow_pickle=False) as p:
        lat_n, lon_n = p["lat_deg"][p["rows"]], p["lon_deg"][p["cols"]]
        eq_range = p["monthly_range_median_m"].astype(float)
    lon_n = np.where(lon_n > 180, lon_n - 360, lon_n)
    panels = len(runs) + 1
    vmax = 0.0
    ranges = {}
    for model, r in runs.items():
        months = r["month_max"] - r["month_min"]
        filled = np.isfinite(months).all(axis=1)
        ranges[model] = np.median(months[filled], axis=0) if filled.any() else np.full(months.shape[1], np.nan)
        vmax = max(vmax, np.nanpercentile(ranges[model], 99.5))
    fig, axes = plt.subplots(1, panels, figsize=(5.2 * panels, 5), constrained_layout=True)
    for ax, (model, r) in zip(axes, runs.items()):
        im = draw(ax, r, ranges[model], vmin=0, vmax=vmax, cmap="viridis")
        ax.set_title(f"{model}: median monthly range")
    sel = (lon_n > VIEW["lon"][0] - 1) & (lon_n < VIEW["lon"][1] + 1) & (lat_n > VIEW["lat"][0] - 1) & (lat_n < VIEW["lat"][1] + 1)
    axes[-1].scatter(lon_n[sel], lat_n[sel], c=eq_range[sel], s=4, marker="s", vmin=0, vmax=vmax, cmap="viridis")
    axes[-1].set_title("equilibrium: median monthly range")
    for ax in axes:
        ax.set_xlim(*VIEW["lon"]); ax.set_ylim(*VIEW["lat"]); ax.set_aspect(1); ax.set_facecolor("#9a8f80")
    fig.colorbar(im, ax=axes, shrink=0.8, label="metres")
    fig.savefig(OUT / "range.png", dpi=110)
    plt.close(fig)

    fig, axes = plt.subplots(1, len(runs), figsize=(5.6 * len(runs), 5), constrained_layout=True, squeeze=False)
    vmax = max(np.percentile(r["speed_max"], 99.9) for r in runs.values())
    for ax, (model, r) in zip(axes[0], runs.items()):
        im = draw(ax, r, r["speed_max"], vmin=0, vmax=vmax, cmap="magma")
        ax.set_title(f"{model}: largest depth-averaged speed")
        ax.set_xlim(*VIEW["lon"]); ax.set_ylim(*VIEW["lat"]); ax.set_aspect(1); ax.set_facecolor("#9a8f80")
        for name, z in study.ZONES.items():
            ax.add_patch(plt.Rectangle((z["lon"][0], z["lat"][0]), np.ptp(z["lon"]), np.ptp(z["lat"]), fill=False, ec="w", lw=0.8, ls="--"))
    fig.colorbar(im, ax=axes[0], shrink=0.8, label="m/s")
    fig.savefig(OUT / "speed.png", dpi=110)
    plt.close(fig)
    (OUT / "sources.json").write_text(json.dumps({m: str(study.HERE / "runs" / m / f"fecunditatis_cd{args.drag:g}" / "reduced.npz")
                                                  for m in runs}, indent=1) + "\n")
    print(OUT)


if __name__ == "__main__":
    main()
