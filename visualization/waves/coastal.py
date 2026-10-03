"""Display the coastal refinement and recovered geographic wind coverage."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(fig, output, inputs, display):
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=160, facecolor="white")
    plt.close(fig)
    record = dict(schema="terluna.visualization.coastal-waves/1",
                  evidence="Displays computed coastal wave fields and recovered atmospheric records with their original sampling and model assumptions.",
                  reading_rule="Read wave heights as responses to the prescribed wind episode. CM1 ring snapshots and GCM model-level means retain separate atmospheric realisations and time conventions.",
                  producer=dict(path="visualization/waves/coastal.py", sha256=digest(Path(__file__))),
                  inputs=[dict(path=str(p.relative_to(ROOT)), sha256=digest(p)) for p in inputs],
                  output=dict(path=output.name, sha256=digest(output)),
                  matplotlib_version=matplotlib.__version__, display=display)
    output.with_suffix(".json").write_text(json.dumps(record, indent=2)+"\n")


def coastal(output):
    path = ROOT/"climate/waves/results/coastal.json"
    product = json.loads(path.read_text())
    if product["schema"] != "terluna.climate.coastal-resolution/1":
        raise ValueError("Unsupported coastal product")
    cases = {c["name"]: c for c in product["cases"]}
    fig, axes = plt.subplots(2, 2, figsize=(11, 10))
    fig.subplots_adjust(top=.85, bottom=.16, left=.09, right=.87, hspace=.34, wspace=.28)
    fig.suptitle("Resolving eastern Smythii's coast", x=.09, y=.975, ha="left", fontsize=20, weight="bold")
    fig.text(.09, .929, "12 hours of 4.001 m/s eastward wind · lunar gravity · flooded LOLA terrain", fontsize=11)
    fig.text(.09, .9, "The basin supplies full directional spectra every 15 minutes. Each map shows computed nodes.", fontsize=10)
    for ax, stride, km in zip(axes.flat, (4, 2, 1), (7.6, 3.8, 1.9)):
        c = cases[f"coast_s{stride}_d36"]
        values, wet = np.asarray(c["final_values"]), np.asarray(c["wet"])
        hs = np.ma.array(values[..., 3], mask=~wet)
        mesh = ax.pcolormesh(c["longitude_deg"], c["latitude_deg"], hs, cmap="viridis", vmin=0, vmax=1.5, shading="nearest")
        np.testing.assert_array_equal(np.ma.compressed(mesh.get_array()), values[..., 3][wet])
        ax.set_title(f"{km:g} km wave grid ({c['wet_nodes']:,} wet nodes)")
    c, f = cases["coast_s2_d36"], cases["coast_s1_d36"]
    a, b = np.asarray(c["final_values"]), np.asarray(f["final_values"])[::2, ::2]
    common = np.asarray(c["wet"]) & np.asarray(f["wet"])[::2, ::2]
    difference = np.ma.array(b[..., 3]-a[..., 3], mask=~common)
    limit = max(.05, float(np.ceil(np.max(np.abs(difference))/.05)*.05))
    delta = axes[1, 1].pcolormesh(c["longitude_deg"], c["latitude_deg"], difference,
                                 cmap="RdBu_r", vmin=-limit, vmax=limit, shading="nearest")
    np.testing.assert_array_equal(np.ma.compressed(delta.get_array()), (b[..., 3]-a[..., 3])[common])
    axes[1, 1].set_title("Height change: 1.9 km minus 3.8 km")
    for ax in axes.flat:
        ax.set_facecolor("#e9e4d9")
        ax.set(xlabel="Longitude (°E)", ylabel="Latitude (°N)")
        ax.set_aspect("equal")
    fig.colorbar(mesh, cax=fig.add_axes([.91, .55, .017, .27]), label="Significant wave height (m)")
    fig.colorbar(delta, cax=fig.add_axes([.91, .18, .017, .27]), label="Height change (m)")
    check = product["checks"]["coast_s2_d36_to_coast_s1_d36"]
    direction = product["checks"]["coast_s2_d36_to_coast_s2_d72"]
    fig.text(.09, .105, f"Interior 3.8 → 1.9 km differences: maximum Hs {check['hs_relative_max']:.1%}; mean period {check['tm01_relative_max']:.1%}.", fontsize=10)
    fig.text(.09, .075, f"36 → 72 child directions at 3.8 km: maximum Hs {direction['hs_relative_max']:.1%}. Parent keeps 36 directions.", fontsize=10)
    fig.text(.09, .043, "Conditional wind episode. Beach slopes, breaking calibration and surf occurrence require further work.", fontsize=10)
    save(fig, output, [path], "Native node Hs at hour 12 with dry nodes masked. Differences use common wet nodes including boundaries; quoted errors use the source product's interior comparison.")


def coverage(output):
    ring_path = ROOT/"climate/waves/results/ring_coverage.json"
    gcm_path = ROOT/"climate/waves/results/gcm_coverage.json"
    rings, gcm = json.loads(ring_path.read_text()), json.loads(gcm_path.read_text())
    if rings["schema"] != "terluna.climate.cm1-wind-coverage/1" or gcm["schema"] != "terluna.climate.gcm-regional-wind-context/1":
        raise ValueError("Unsupported wind coverage")
    inputs = [ring_path, gcm_path]
    products = []
    for record in [*rings["cases"], gcm]:
        path = ROOT/record["product"]["path"]
        if digest(path) != record["product"]["sha256"]:
            raise ValueError("Wind product differs from its record")
        with np.load(path, allow_pickle=False) as f:
            products.append({k: f[k].copy() for k in f.files})
        inputs.append(path)
    fig = plt.figure(figsize=(12, 10))
    grid = fig.add_gridspec(2, 2, left=.09, right=.94, top=.84, bottom=.17, hspace=.48, wspace=.3)
    ax, regional, series = fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1]), fig.add_subplot(grid[1, :])
    fig.suptitle("Recovered winds: tracks and regional means", x=.09, y=.975, ha="left", fontsize=20, weight="bold")
    fig.text(.09, .928, "CM1: 237 three-hourly surface snapshots per ring · GCM: 1,200 regional mean fields", fontsize=11)
    fig.text(.09, .899, "Independent atmospheric runs; each product retains its own coordinates, height and time convention.", fontsize=10)
    colors = ("#2979a0", "#bd6033", "#7b5c9c")
    for product, record, color in zip(products, rings["cases"], colors):
        lon, lat = product["longitude_deg"], product["latitude_deg"]
        order = np.argsort(lon)
        ax.plot(lon[order], lat[order], color=color, lw=1.5, label=f"{record['tilt_deg']:g}° tilt, {record['node_deg']:g}° node")
    ax.add_patch(Rectangle((76, -9), 26, 31, fill=False, ec="#202b38", lw=1.5))
    ax.set(xlim=(0, 360), ylim=(-80, 80), xlabel="Longitude (°E)", ylabel="Latitude (°N)", title="Three recovered CM1 ring tracks")
    ax.legend(fontsize=8, loc="lower center", ncol=1, frameon=False)
    d = products[-1]
    xx, yy = np.meshgrid(d["longitude_deg"], d["latitude_deg"])
    u, v = d["eastward_m_s"][-1], d["northward_m_s"][-1]
    regional.scatter(xx, yy, c=np.where(d["water"][-1], "#bddce8", "#d4c8b4").ravel(), s=55)
    arrows = regional.quiver(xx, yy, u, v, scale=32, width=.005, color="#203449")
    np.testing.assert_array_equal(arrows.U, u.ravel())
    np.testing.assert_array_equal(arrows.V, v.ravel())
    regional.quiverkey(arrows, .83, 1.09, 4, "4 m/s", labelpos="E", coordinates="axes")
    regional.axhline(0, color=colors[0], lw=1.5)
    regional.add_patch(Rectangle((76, -9), 26, 31, fill=False, ec="#202b38", lw=1.5))
    regional.set(xlim=(70, 110), ylim=(-15, 30), xlabel="Longitude (°E)", ylabel="Latitude (°N)", title="GCM: final saved 71.5-hour mean")
    equator = products[0]
    distance = np.abs(equator["longitude_deg"]-90)
    column = np.argmin(np.where(equator["water"], distance, np.inf))
    for key, label, color in (("eastward_10m_m_s", "Eastward", colors[0]), ("northward_10m_m_s", "Northward", colors[1])):
        line, = series.plot(equator["time_s"]/86400, equator[key][:, column], label=label, color=color, lw=1.3)
        np.testing.assert_array_equal(line.get_ydata(), equator[key][:, column])
    series.axhline(0, lw=.6, color="#333333")
    series.set(xlabel="CM1 elapsed Earth days", ylabel="10 m wind component (m/s)", title=f"Equatorial water column at {equator['longitude_deg'][column]:.2f}°E")
    series.legend(frameon=False, ncol=2)
    series.grid(alpha=.2)
    fig.text(.09, .102, "Blue GCM nodes are water; tan nodes are land. Vectors are at the lowest model level, sigma 0.9833.", fontsize=10)
    fig.text(.09, .073, "Only the equatorial CM1 track crosses this sea. The tilted rings extend coverage to higher latitudes.", fontsize=10)
    fig.text(.09, .043, "Basin forcing needs hourly or three-hourly 10 m vectors over the full water area.", fontsize=10)
    save(fig, output, inputs, "Ring coordinates and one native CM1 column. Native GCM last-window vectors at the lowest model level; coloured nodes retain the GCM surface mask. Rectangles mark the Smythii–Marginis selection.")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--coverage-only", action="store_true")
    args = p.parse_args()
    coverage(HERE/"results/wind_coverage.png")
    if not args.coverage_only:
        coastal(HERE/"results/coastal.png")
