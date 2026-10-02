"""Display the climate wind-wave data product without recomputing its physics."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SCHEMA = "terluna.climate.wind-waves/1"
HOURS = (1, 3, 6, 12, 24, 48)
COLORS = {"earth": "#286AA6", "moon": "#BB6423"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_product(path: Path) -> tuple[dict, list[float], dict, list[dict]]:
    """Require the complete plotted subset and its stated scenario conditions."""
    product = json.loads(path.read_text())
    if product.get("schema") != SCHEMA:
        raise ValueError(f"Expected {SCHEMA}: {path}")
    required_units = dict(hs_m="m", mean_period_s="s", wind_m_s="m/s",
                          fetch_km="km", duration_h="hours")
    if any(product.get("units", {}).get(key) != unit for key, unit in required_units.items()):
        raise ValueError("Wind-wave product uses unexpected units")
    cases = {case["name"]: case for case in product["cases"]}
    if len(cases) != 6:
        raise ValueError("Expected three Earth and three Moon forcing cases")
    for case in cases.values():
        if case["depth_m"] != 1000 or case["length_m"] != 100_000 or case["seed_height_m"] != 0:
            raise ValueError("Figure assumptions require 1000 m depth, 100 km fetch and initially calm water")
    winds = sorted({case["wind_m_s"] for case in cases.values()})
    if len(winds) != 3:
        raise ValueError("Expected three wind magnitudes")
    selected = [row for row in product["rows"] if row["fetch_km"] == 100]
    series = {}
    for wind in winds:
        for world in ("earth", "moon"):
            matching = [case for case in cases.values()
                        if case["name"].startswith(world + "_") and case["wind_m_s"] == wind]
            if len(matching) != 1:
                raise ValueError("Missing or duplicate Earth/Moon wind case")
            case = matching[0]
            for build in ("stock", "patched"):
                rows = sorted((row for row in selected if row["case"] == case["name"]
                               and row["build"] == build), key=lambda row: row["duration_h"])
                if tuple(row["duration_h"] for row in rows) != HOURS:
                    raise ValueError("Plotted series is incomplete or has duplicate times")
                for row in rows:
                    if row["gravity_m_s2"] != case["gravity_m_s2"] or row["wind_m_s"] != wind:
                        raise ValueError("Plotted forcing differs from the case metadata")
                    if any(not math.isfinite(row[key]) or row[key] < 0
                           for key in ("hs_m", "mean_period_s")):
                        raise ValueError("Invalid wave height or period")
                series[wind, world, build] = rows
    if len(selected) != 72:
        raise ValueError("Unexpected extra rows at the plotted fetch")
    return product, winds, series, selected


def plot(input_path: Path, output_path: Path) -> dict:
    product, winds, series, selected = read_product(input_path)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.labelcolor": "#253244", "text.color": "#253244",
                         "xtick.color": "#425166", "ytick.color": "#425166"})
    figure, axes = plt.subplots(2, 3, figsize=(13.6, 8.7), sharex=True, sharey="row")
    figure.subplots_adjust(left=0.075, right=0.98, bottom=0.17, top=0.77,
                           hspace=0.18, wspace=0.13)
    figure.suptitle("Wave growth under Earth and lunar gravity", x=0.075, y=0.97,
                   ha="left", fontsize=20, weight="bold")
    figure.text(0.075, 0.918,
                "100 km fetch  ·  1,000 m water depth  ·  uniform steady wind  ·  initially calm",
                fontsize=12)
    legend = [Line2D([], [], color=COLORS[world], lw=2.6, label=label)
              for world, label in (("earth", "Earth gravity"), ("moon", "Lunar gravity"))]
    legend.extend([Line2D([], [], color="#626B76", lw=2.4, label="Gravity-scaled AGROW"),
                   Line2D([], [], color="#626B76", lw=1.7, ls="--", alpha=0.65,
                          label="Stock SWAN")])
    figure.legend(handles=legend, loc="upper left", bbox_to_anchor=(0.068, 0.89),
                  ncol=4, frameon=False, columnspacing=2.1, handlelength=2.6)
    for column, wind in enumerate(winds):
        axes[0, column].set_title(f"Wind  {wind:.3f} m/s", fontsize=13, pad=11)
        for row_index, quantity in enumerate(("hs_m", "mean_period_s")):
            axis = axes[row_index, column]
            for world in ("earth", "moon"):
                for build in ("stock", "patched"):
                    rows = series[wind, world, build]
                    line, = axis.plot([item["duration_h"] for item in rows],
                                      [item[quantity] for item in rows],
                                      color=COLORS[world],
                                      ls="--" if build == "stock" else "-",
                                      lw=1.7 if build == "stock" else 2.6,
                                      alpha=0.58 if build == "stock" else 1,
                                      marker=None if build == "stock" else "o",
                                      ms=4.2, zorder=2 if build == "stock" else 3)
                    if list(line.get_ydata()) != [item[quantity] for item in rows]:
                        raise RuntimeError("Rendered values differ from the source product")
            axis.set_xscale("log")
            axis.set_xlim(0.88, 55)
            axis.set_xticks(HOURS, [str(value) for value in HOURS])
            axis.minorticks_off()
            axis.grid(axis="y", color="#DAE1E9", lw=0.7)
            axis.spines["bottom"].set_color("#AAB5C2")
            axis.spines["left"].set_color("#AAB5C2")
            if row_index:
                axis.set_xlabel("Wind duration (hours)", labelpad=9)
    axes[0, 0].set_ylabel("Significant wave height Hₛ (m)", labelpad=10)
    axes[1, 0].set_ylabel("Mean period Tₘ₀₁ (s)", labelpad=10)
    for row_index, quantity in enumerate(("hs_m", "mean_period_s")):
        maximum = max(item[quantity] for item in selected)
        axes[row_index, 0].set_ylim(0, 1.08 * maximum)
    figure.text(0.075, 0.087,
                "Conditional numerical experiments. No incoming swell or currents. "
                "Wind statistics select forcing magnitudes only.", fontsize=10)
    figure.text(0.075, 0.052,
                "Solid lines use a gravity-scaled wind-wave onset term. "
                "Lunar wind-input physics remains empirically unvalidated.", fontsize=10)
    checks = product["checks"]
    numerical_note = ("Early growth is timestep-sensitive. Resolution checks cover only "
                      f"the 48 h lunar endpoint at {winds[1]:.3f} m/s.")
    if not (checks["endpoint_resolution_pass"] and checks["spectral_edges_pass"]
            and checks["gravity_similarity"]["patched_pass"]):
        numerical_note = "One or more numerical checks require review; consult the source product."
    figure.text(0.075, 0.018, numerical_note, color="#9C441E", fontsize=10)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=170, facecolor="white")
    plt.close(figure)
    record = dict(schema="terluna.visualization.wind-waves/1",
                  input=dict(path=str(input_path.resolve()), sha256=sha256(input_path),
                             schema=product["schema"], producer=product["producer"]),
                  source=dict(path="visualization/waves/plot.py", sha256=sha256(Path(__file__))),
                  output=dict(path=output_path.name, sha256=sha256(output_path)),
                  matplotlib_version=matplotlib.__version__, evidence=product["evidence"],
                  reading_rule=product["reading_rule"], checks=checks,
                  display=dict(fetch_km=100, hours=list(HOURS), horizontal_scale="logarithmic",
                               connection_between_samples="straight line in plotted coordinates",
                               numerical_note=numerical_note,
                               fields=["hs_m", "mean_period_s"]), plotted_rows=selected)
    output_path.with_suffix(".json").write_text(json.dumps(record, indent=2) + "\n")
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ROOT / "climate/waves/results/waves.json")
    parser.add_argument("--output", type=Path, default=HERE / "results/waves.png")
    args = parser.parse_args()
    plot(args.input, args.output)
    print(args.output)


if __name__ == "__main__":
    main()
