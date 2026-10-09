"""Render four conditional floater screens from the coupled JSON product only.

Run from the repository root:
    python visualization/floater_viability/plot.py

Physics stays in research/studies/floater_viability. Every plotted point is
exported with its product/source hashes. Generated images/data are ignored.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, ScalarFormatter

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
DEFAULT_PRODUCT = ROOT / "research/studies/floater_viability/results/floater_viability.json"
COLORS = ("#126A83", "#C16921", "#70479B", "#27816A")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def close(a: float, b: float) -> bool:
    return abs(a - b) < 1e-8


def select_data(product: dict) -> dict:
    if product.get("schema") != "terluna.research.floater-viability/1":
        raise ValueError("Expected coupled floater-viability/1 data product")
    structural = product.get("structural_size_scans")
    navigation = product.get("navigation")
    if not structural or not navigation:
        raise ValueError("Product needs structural_size_scans and navigation; regenerate the study product")
    strength_rows = [r for r in structural
                     if r["height_km"] == 10 and r["traits"]["architecture"] == "sphere"
                     and close(r["traits"]["living_dry_kg_m2"], 1.0)]
    if len(strength_rows) < 3:
        raise ValueError("Need at least three stress-allowable size scans for the reference sphere")
    chosen_materials = ("kurata_2022_cellophane", "kurata_2022_polypropylene",
                        "mfc_2025_incomplete_conditions", "wu_2014_chitin_incomplete_conditions")
    envelope = [r for r in product["envelope_comparators"] if r["comparator"] in chosen_materials]
    if not envelope or any("hydrogen_chemical_w_m2_projected" not in r for r in envelope):
        raise ValueError("Envelope product needs hydrogen_chemical_w_m2_projected on every comparator")
    water = [r for r in product["environment"]["carbon_water_requirements"]
             if r["height_km"] == 10 and r["co2_ppm"] == 400
             and close(r["assimilation_to_npp"], 1.0)
             and close(r["npp_g_c_m2_year"], 1000.0)]
    motion = [r for r in navigation["relative_speed_cases"]
              if r["shape"] in ("sphere_sensitivity", "streamlined_hypothesis")]
    if len(water) != 6 or len(motion) < 8:
        raise ValueError("Expected six water points and biological-input motion sensitivities")
    return {"structural": strength_rows, "envelope": envelope, "water": water,
            "motion": motion, "navigation_assumptions": navigation["source_assumptions"],
            "air_10km": navigation["air_10km"]}


def style_axis(ax, title, subtitle):
    ax.set_title(title, loc="left", fontsize=13.5, fontweight="bold", pad=43)
    ax.text(0, 1.035, subtitle, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=9, color="#4B5965", linespacing=1.35)
    ax.grid(axis="both", which="major", color="#DCE3E8", linewidth=0.7)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["bottom", "left"]].set_color("#A8B2BC")
    ax.tick_params(labelsize=9.5, color="#A8B2BC")


def render(product_path: Path, output_dir: Path) -> dict:
    product = json.loads(product_path.read_text())
    data = select_data(product)
    product_hash = digest(product_path)
    script_hash = digest(Path(__file__))
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.labelcolor": "#263845", "text.color": "#172E3D",
                         "figure.facecolor": "white", "axes.facecolor": "white",
                         "svg.hashsalt": "terluna-floater-viability-v1"})
    fig, axes = plt.subplots(2, 2, figsize=(14, 11.4))
    fig.subplots_adjust(left=.085, right=.97, bottom=.145, top=.82, wspace=.26, hspace=.69)
    fig.text(.085, .955, "Biological floaters: four coupled requirements", fontsize=22, weight="bold")
    fig.text(.085, .921, "Conditional material and physiology sensitivities • No integrated organism or flight system is validated",
             fontsize=11, color="#586975")

    ax = axes[0, 0]
    rows = sorted(data["structural"], key=lambda r: r["traits"]["allowable_mpa"])
    base = rows[0]["traits"]
    threshold = base["maximum_lift_utilization"]
    for i, scan in enumerate(rows):
        points = sorted(scan["rows"], key=lambda r: r["diameter_m"])
        ax.plot([r["diameter_m"] for r in points], [r["gas_fraction"] for r in points],
                "o-", ms=3.5, lw=2, color=COLORS[i % len(COLORS)],
                label=f"Allowable stress {scan['traits']['allowable_mpa']:g} MPa")
    ax.axhspan(.08, threshold, color="#E5F0EB", zorder=0)
    ax.axhline(threshold, color="#456E5D", lw=1.1, ls="--")
    ax.text(3.0, threshold * .83, f"Mass limit with {100*(1-threshold):g}% lift reserve", fontsize=8.5, color="#456E5D")
    ax.set(xscale="log", yscale="log", xlim=(2, 1000), ylim=(.08, 15),
           xlabel="Envelope diameter (m)", ylabel="Required fraction of full supported lift")
    ax.set_xticks([2, 10, 30, 100, 300, 1000])
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.set_yticks([.1, .3, .7, 1, 3, 10])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.legend(loc="upper right", fontsize=8.7, framealpha=.94, edgecolor="none")
    style_axis(ax, "A  Small and large spheres face different mass limits",
               f"10 km air; sphere; 1 kg dry living biomass/m² projected\n"
               f"Material density {base['material_density_kg_m3']:g} kg/m³; {base['relative_gust_m_s']:g} m/s relative gust; traits assumed")

    ax = axes[0, 1]
    materials = (("kurata_2022_cellophane", "Cellophane: 0.173 Barrer", "-"),
                 ("kurata_2022_polypropylene", "Polypropylene: 6.69 Barrer", "-"),
                 ("mfc_2025_incomplete_conditions", "MFC: 10.7 Barrer*", "--"),
                 ("wu_2014_chitin_incomplete_conditions", "Chitin: 0.024 Barrer*", "--"))
    for i, (name, label, ls) in enumerate(materials):
        points = sorted((r for r in data["envelope"] if r["comparator"] == name), key=lambda r: r["thickness_um"])
        ax.plot([r["thickness_um"] for r in points], [r["hydrogen_chemical_w_m2_projected"] for r in points],
                marker="o", ms=3.5, lw=2, ls=ls, color=COLORS[i], label=label)
    ax.set(xscale="log", yscale="log", xlim=(1, 1000),
           xlabel="Intact barrier thickness (µm)", ylabel="H₂ replacement chemical power (W/m² projected)")
    ax.set_xticks([1, 10, 100, 1000])
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.legend(loc="upper right", fontsize=8.5, framealpha=.94, edgecolor="none")
    style_axis(ax, "B  Thicker barriers trade lower leakage for mass",
               "100 kPa total; 98% H₂; full sphere, gas area/projected = 4\n"
               "LHV of lost H₂; conversion losses excluded. *Conditions incomplete")

    ax = axes[1, 0]
    leaf_temperatures = sorted(set(r["leaf_c"] for r in data["water"]))
    for i, leaf in enumerate(leaf_temperatures):
        points = sorted((r for r in data["water"] if close(r["leaf_c"], leaf)), key=lambda r: r["relative_humidity"])
        delta = leaf - points[0]["air_c"]
        ax.plot([100 * r["relative_humidity"] for r in points], [r["water_kg_m2_year"] for r in points],
                "o-", lw=2.3, ms=5, color=COLORS[i], label=f"Leaf {delta:+g} K above air")
        end = points[-1]
        ax.annotate(f"{end['water_kg_m2_year']:.1f}", (98, end["water_kg_m2_year"]),
                    xytext=(-8, 9), textcoords="offset points", ha="right", fontsize=9, color=COLORS[i])
    air = data["water"][0]["air_c"]
    ax.set(xlim=(58, 100), ylim=(0, None), xlabel="Ambient relative humidity (%)",
           ylabel="Water demand (kg/m² projected/year)")
    ax.set_xticks([60, 70, 80, 90, 98])
    ax.legend(loc="upper right", fontsize=9, frameon=False)
    style_axis(ax, "C  Warm leaves still lose water in humid air",
               f"10 km; air {air:g}°C; CO₂ 400 ppm; imposed leaf net carbon 1 kg/m²/year\n"
               "Diffusion requirement only; capture and nighttime losses not included")

    ax = axes[1, 1]
    for i, (shape, label) in enumerate((("sphere_sensitivity", "Sphere sensitivity"),
                                       ("streamlined_hypothesis", "Streamlined hypothesis"))):
        points = sorted((r for r in data["motion"] if r["shape"] == shape), key=lambda r: r["airspeed_m_s"])
        points = [r for r in points if r["airspeed_m_s"] > 0]
        ax.plot([r["airspeed_m_s"] for r in points], [r["mean_input_w_m2"] for r in points],
                "o-", lw=2.3, ms=4.5, color=COLORS[i],
                label=f"{label}, Cᴅ = {points[0]['drag_coefficient']:g}")
    gpp = data["navigation_assumptions"]["reference_gpp_chemical_equivalent_w_m2"]
    ax.axhline(gpp, ls="--", color="#697B88", lw=1.3)
    ax.text(.28, gpp * 1.19, "Reference gross carbon energy; not spare power", fontsize=8.5, color="#586975")
    eta = data["navigation_assumptions"]["biological_overall_efficiency"]
    rho = data["air_10km"]["density_kg_m3"]
    ax.set(xscale="log", yscale="log", xlim=(.25, 5),
           xlabel="Continuous speed relative to local air (m/s)",
           ylabel="Motion chemical input (W/m² projected)")
    ax.set_xticks([.25, .5, 1, 2, 5])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.legend(loc="upper left", fontsize=8.6, framealpha=.94, edgecolor="none")
    style_axis(ax, "D  Relative airspeed controls the motion cost",
               f"10 km; air density {rho:g} kg/m³; assumed overall efficiency {100*eta:g}%\n"
               "Frontal area/projected = 1; continuous motion; no proven wind corridor")

    fig.text(.085, .087, "Reading rule: panels use separate stated sensitivity conditions, not one demonstrated design. Lines join computed points; they are not fitted empirical curves.",
             fontsize=9, color="#435966")
    fig.text(.085, .062, "Passing a mass or energy screen leaves water supply, gas purity, seams, repair, reproduction, weather avoidance and lifecycle unresolved.",
             fontsize=9, color="#435966")
    fig.text(.085, .036, f"Terluna • floater-viability/1 • product SHA-256 {product_hash[:16]} • plot SHA-256 {script_hash[:16]}",
             fontsize=8.2, color="#72828D")

    stem = output_dir / "floater_viability"
    fig.savefig(stem.with_suffix(".png"), dpi=180, facecolor="white", metadata={"Software": "Terluna floater_viability plot.py"})
    fig.savefig(stem.with_suffix(".svg"), facecolor="white", metadata={"Date": None, "Creator": "Terluna floater_viability plot.py"})
    plt.close(fig)
    plotted_path = output_dir / "plotted_data.json"
    plotted_path.write_text(json.dumps({"schema": "terluna.visualization.floater-viability-points/1",
                                      "product_sha256": product_hash, "points": data}, indent=2) + "\n")
    # JSON round trip verifies exported values are exactly the selected model points.
    if json.loads(plotted_path.read_text())["points"] != data:
        raise AssertionError("Plotted-point round trip changed values")
    manifest = {"schema": "terluna.visualization.floater-viability-manifest/1",
                "product": str(product_path.relative_to(ROOT)), "product_sha256": product_hash,
                "plot_source_sha256": script_hash, "producer": product["producer"],
                "outputs": {p.name: digest(p) for p in (stem.with_suffix(".png"), stem.with_suffix(".svg"), plotted_path)},
                "evidence": product["evidence"],
                "validation": "Exact plotted-point JSON round trip; scientific rendering must also be visually inspected."}
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--product", type=Path, default=DEFAULT_PRODUCT)
    parser.add_argument("--output-dir", type=Path, default=HERE / "outputs")
    args = parser.parse_args()
    result = render(args.product.resolve(), args.output_dir.resolve())
    print(json.dumps({"product_sha256": result["product_sha256"], "outputs": list(result["outputs"])}, indent=2))
