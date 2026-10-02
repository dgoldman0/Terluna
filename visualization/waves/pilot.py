"""Render maps and profiles from the conditional basin/slope data product."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def plot(product_path, output):
    product = json.loads(product_path.read_text())
    if product.get("schema") != "terluna.climate.wave-basin-pilot/1":
        raise ValueError("Unexpected basin product schema")
    expected = ["time_s", "longitude_deg", "latitude_deg", "hs_m", "peak_period_s", "mean_period_s", "depth_m", "direction_deg", "breaking_fraction"]
    if product["units"]["basin_columns"] != expected:
        raise ValueError("Unexpected basin columns")
    cases = {c["name"]: c for c in product["basin_cases"]}
    fields = product["basin_fields"]
    max_h = max(float(np.asarray(f["hourly_values"])[-1, :, 3].max()) for f in fields.values())
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.subplots_adjust(top=.86, bottom=.14, left=.09, right=.89, hspace=.38, wspace=.3)
    fig.suptitle("A wave pilot for Smythii–Marginis", x=.09, y=.975, ha="left", fontsize=20, weight="bold")
    fig.text(.09, .932, "Nominal 4.001 m/s winds for 12 hours  ·  lunar gravity  ·  air density from CM1", fontsize=11)
    spatial_error = product["checks"]["smythii_east_half_grid"]["hs_final_relative_max"]
    fig.text(.09, .903, f"28% water scenario. Local grid differences reach {spatial_error:.0%}; surf locations need finer coastal grids.", fontsize=10)
    colors = {"smythii_east": "#22628b", "smythii_turn": "#ba5a28"}
    for ax, (name, label) in zip(axes[0], (("smythii_east", "Wind stays eastward"),
                                         ("smythii_turn", "Wind turns north during hours 5–6"))):
        f, case = fields[name], cases[name]
        points, cube = np.asarray(f["points"]), np.asarray(f["hourly_values"])
        if cube.shape != (case["hours"] + 1, len(points), 9) or not np.isfinite(cube).all():
            raise ValueError("Incomplete plotted field")
        delta = case["stride"] / 4
        lon = np.arange(points[:, 0].min(), points[:, 0].max() + delta / 2, delta)
        lat = np.arange(points[:, 1].min(), points[:, 1].max() + delta / 2, delta)
        values = np.full((len(lat), len(lon)), np.nan)
        j = np.rint((points[:, 0] - lon[0]) / delta).astype(int)
        i = np.rint((points[:, 1] - lat[0]) / delta).astype(int)
        values[i, j] = cube[-1, :, 3]
        np.testing.assert_array_equal(values[i, j], cube[-1, :, 3])
        ax.set_facecolor("#eeeae2")
        mesh = ax.pcolormesh(lon, lat, values, shading="nearest", cmap="viridis", vmin=0, vmax=max_h)
        x, y = case["offshore_lon_lat"]
        ax.plot(x, y, "o", color="white", mec="#212c39", ms=6)
        ax.set(title=label, xlabel="Longitude (°E)", ylabel="Latitude (°N)")
        ax.set_aspect("equal")
        k = int(np.argmin(np.sum((points - [x, y])**2, axis=1)))
        axes[1, 0].plot(cube[:, k, 0] / 3600, cube[:, k, 3], marker="o", ms=3,
                        color=colors[name], label=label)
    cax = fig.add_axes([.915, .552, .017, .287])
    fig.colorbar(mesh, cax=cax, label="Significant wave height at hour 12 (m)")
    axes[1, 0].set(title="Offshore node marked on the maps", xlabel="Hours after initially calm water", ylabel="Significant wave height (m)")
    axes[1, 0].legend(frameon=False, fontsize=9)
    axes[1, 0].grid(alpha=.2)
    axes[1, 0].axvspan(0, 2, color="#c95d24", alpha=.07)
    for row in product["slopes"]:
        if row["cells"] != 400 or row["breaker_index"] != .73:
            continue
        data = np.asarray(row["rows"])
        line, = axes[1, 1].plot(data[:, 4], data[:, 1], label=f"Slope 1:{round(1 / row['slope'])}")
        np.testing.assert_array_equal(line.get_ydata(), data[:, 1])
    axes[1, 1].set(title="The offshore spectrum on assumed slopes", xlabel="Water depth (m)", ylabel="Significant wave height (m)", xlim=(15, 0))
    axes[1, 1].legend(frameon=False)
    axes[1, 1].grid(alpha=.2)
    fig.text(.09, .08, "Conditional responses. The vector turn briefly lowers wind speed to 2.83 m/s. Initial growth is timestep-sensitive.", fontsize=9)
    fig.text(.09, .054, "Slopes isolate shoaling and depth breaking. Breaking coefficients need calibration at lunar gravity.", fontsize=9)
    fig.text(.09, .029, "Surf frequency requires wind histories; crests and shoreline run-up require individual-wave simulations.", fontsize=10)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=160, facecolor="white")
    plt.close(fig)
    record = dict(schema="terluna.visualization.wave-pilot/1", evidence=product["evidence"], reading_rule=product["reading_rule"],
                  input=dict(path=str(product_path), sha256=digest(product_path), schema=product["schema"]),
                  producer=dict(path="visualization/waves/pilot.py", sha256=digest(Path(__file__))),
                  output=dict(path=output.name, sha256=digest(output)), matplotlib_version=matplotlib.__version__,
                  display="Map cells display native coarse SWAN node values. Straight lines connect hourly node values. Slopes show the gamma=0.73, 400-cell cases.")
    output.with_suffix(".json").write_text(json.dumps(record, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ROOT / "climate/waves/results/pilot.json")
    parser.add_argument("--output", type=Path, default=HERE / "results/pilot.png")
    args = parser.parse_args()
    plot(args.input, args.output)
    print(args.output)
