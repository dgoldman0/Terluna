"""Render the floater study's screens from its coupled JSON product only.

Run from the repository root:
    python visualization/floater_viability/plot.py

Physics stays in research/studies/floater_viability. Two figures: the four
first-round requirements, and the giants' limits, ages, water and day length.
Every plotted point is exported with the product and plot hashes. Generated
images and data are ignored by Git.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, ScalarFormatter

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
DEFAULT_PRODUCT = ROOT / "research/studies/floater_viability/results/floater_viability.json"
SCHEMA = "terluna.research.floater-viability/2"
# Categorical slots 1-4 of the validated reference palette, in fixed order (light surface).
COLORS = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100")
SURFACE = "#fcfcfb"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
TEXT_MUTED = "#6f6e69"
GRID = "#e4e3df"
AXIS = "#b9b8b2"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def close(a: float, b: float) -> bool:
    return abs(a - b) < 1e-8


def select_requirements(product: dict) -> dict:
    """The first-round panels: structure, barrier, water demand and motion cost."""
    structural = product.get("structural_size_scans")
    navigation = product.get("navigation")
    if not structural or not navigation:
        raise ValueError("Product needs structural_size_scans and navigation; regenerate the study product")
    strength_rows = [r for r in structural
                     if r["height_km"] == 10 and r["traits"]["architecture"] == "sphere"
                     and close(r["traits"]["living_dry_kg_m2"], 1.0)]
    if len(strength_rows) < 3:
        raise ValueError("Need at least three stress-allowable size scans for the reference sphere")
    chosen = ("kurata_2022_cellophane", "kurata_2022_polypropylene",
              "mfc_2025_incomplete_conditions", "wu_2014_chitin_incomplete_conditions")
    envelope = [r for r in product["envelope_comparators"] if r["comparator"] in chosen]
    if not envelope or any("hydrogen_chemical_w_m2_projected" not in r for r in envelope):
        raise ValueError("Envelope product needs hydrogen_chemical_w_m2_projected on every comparator")
    water = [r for r in product["environment"]["carbon_water_requirements"]
             if r["height_km"] == 10 and r["co2_ppm"] == 400
             and close(r["assimilation_to_npp"], 1.0) and close(r["npp_g_c_m2_year"], 1000.0)]
    motion = [r for r in navigation["relative_speed_cases"]
              if r["shape"] in ("sphere_sensitivity", "streamlined_hypothesis")]
    if len(water) != 6 or len(motion) < 8:
        raise ValueError("Expected six water points and biological-input motion sensitivities")
    return {"structural": strength_rows, "envelope": envelope, "water": water, "motion": motion,
            "navigation_assumptions": navigation["source_assumptions"], "air_10km": navigation["air_10km"]}


def select_giants(product: dict) -> dict:
    """The giants' panels: carbon surplus by size, ages, rain against need, day length."""
    for key in ("size_limits", "ages", "water", "sailing"):
        if key not in product:
            raise ValueError(f"Product lacks {key}; regenerate the study product")
    rows = product["size_limits"]["10_km"]["rows"]
    size = {name: [dict(diameter_m=r["diameter_m"], surplus_c_kg_m2_year=r["surplus_c_kg_m2_year"],
                        all_gates=r["all_gates"]) for r in rows[name]]
            for name in ("round_collagen", "round", "round_strong", "colony_module")}
    ages = product["ages"]
    targets = ages["target_diameters_m"]

    def round_case(fibre, gpp, route):
        row = next(r for r in ages["round"] if r["fibre"] == fibre and close(r["gpp_c_kg_m2_year"], gpp)
                   and r["hydrogen_route"] == route and close(r["growth_share"], 1.0))
        years = [row["years_to_radius"][str(d/2)] for d in targets]
        return [dict(diameter_m=d, years=y) for d, y in zip(targets, years) if y is not None]

    def colony_case(gpp, route):
        row = next(r for r in ages["colony"] if close(r["gpp_c_kg_m2_year"], gpp)
                   and r["hydrogen_route"] == route and close(r["growth_share"], 1.0))
        return [dict(diameter_m=d, years=y) for d, y in zip(targets, row["years_to_area"].values())]

    age_series = {
        "Colony, photolytic gas, reference light": colony_case(2.0, "photolysis"),
        "Colony, fermented gas, reference light": colony_case(2.0, "fermentation_2_1"),
        "Round, 300 MPa fibre, canopy-like light": round_case("round_strong", 5.0, "photolysis"),
        "Round, 300 MPa fibre, reference light": round_case("round_strong", 2.0, "photolysis"),
    }
    water = product["water"]
    zonal = [dict(latitude_deg=l, rain_mm_day=r) for l, r in
             zip(water["zonal_rain"]["latitude_deg"], water["zonal_rain"]["rain_mm_day"])]
    needs = {n["name"]: n["water_kg_m2_day"] for n in water["needs"]}
    rows = product["sailing"]["day_length"]
    days = [dict(latitude_deg=r["latitude_deg"], height_km=r["height_km"], period_days=r["relative_period_days"])
            for r in rows]
    # A fixed place sees T = relative period x (v + u)/v; every row must give the same month.
    fixed = [r["relative_period_days"]*(r["sun_pace_m_s"] + r["eastward_m_s"])/r["sun_pace_m_s"] for r in rows]
    if not zonal or len(days) < 8 or any(not v for v in size.values()):
        raise ValueError("Giants' panels need size rows, zonal rain and day lengths")
    module = {r["module_diameter_m"] for r in ages["colony"]}
    if len(module) != 1:
        raise ValueError("Colony ages must share one founding module")
    if max(fixed) - min(fixed) > 1e-6:
        raise ValueError("Day-length rows disagree on the fixed-place period")
    return {"size": size, "ages": age_series, "zonal_rain": zonal,
            "needs": {k: needs[k] for k in ("reference", "concentrating_cool_leaves", "concentrating_humid_layer")},
            "day_length": days, "fixed_place_days": fixed[0], "module_diameter_m": module.pop(),
            "start_diameter_m": ages["start_diameter_m"],
            "freezing_height_km": product["environment"]["mean_freezing_height_km"]}


def style_axis(ax, title, subtitle):
    ax.set_title(title, loc="left", fontsize=13, fontweight="bold", pad=40, color=TEXT_PRIMARY)
    ax.text(0, 1.035, subtitle, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=8.8, color=TEXT_SECONDARY, linespacing=1.35)
    ax.grid(axis="both", which="major", color=GRID, linewidth=0.7)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["bottom", "left"]].set_color(AXIS)
    ax.tick_params(labelsize=9.2, color=AXIS, labelcolor=TEXT_SECONDARY)
    ax.xaxis.label.set_color(TEXT_SECONDARY)
    ax.yaxis.label.set_color(TEXT_SECONDARY)


def end_label(ax, x, y, text, color_index, label_y=None):
    """Direct label at a line's end, in text ink beside the series marker; label_y dodges neighbours."""
    ax.annotate(text, (x, y if label_y is None else label_y), xytext=(7, 0), textcoords="offset points",
                va="center", fontsize=8.3, color=TEXT_PRIMARY, annotation_clip=False)
    ax.plot([x], [y], "o", ms=5, color=COLORS[color_index], mec=SURFACE, mew=1.5, zorder=4)


def dodge(values, gap):
    """Label positions at least gap apart, keeping their order and moving upward from the lowest."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    placed = list(values)
    previous = None
    for i in order:
        if previous is not None and placed[i] < previous + gap:
            placed[i] = previous + gap
        previous = placed[i]
    return placed


def setup():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.labelcolor": TEXT_SECONDARY, "text.color": TEXT_PRIMARY,
                         "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                         "svg.hashsalt": "terluna-floater-viability-v2"})


def draw_requirements(data, product_hash, script_hash):
    fig, axes = plt.subplots(2, 2, figsize=(14, 11.4))
    fig.subplots_adjust(left=.085, right=.97, bottom=.145, top=.82, wspace=.26, hspace=.69)
    fig.text(.085, .955, "Biological floaters: four coupled requirements", fontsize=22, weight="bold")
    fig.text(.085, .921, "Material and physiology sensitivities from the first round of the study",
             fontsize=11, color=TEXT_SECONDARY)

    ax = axes[0, 0]
    rows = sorted(data["structural"], key=lambda r: r["traits"]["allowable_mpa"])
    base = rows[0]["traits"]
    threshold = base["maximum_lift_utilization"]
    for i, scan in enumerate(rows[:4]):
        points = sorted(scan["rows"], key=lambda r: r["diameter_m"])
        ax.plot([r["diameter_m"] for r in points], [r["gas_fraction"] for r in points],
                "o-", ms=4, lw=2, color=COLORS[i], label=f"Allowable stress {scan['traits']['allowable_mpa']:g} MPa")
    ax.axhline(threshold, color=TEXT_MUTED, lw=1.1, ls="--")
    ax.text(2.2, threshold*.93, f"Mass limit, {100*(1-threshold):g}% lift reserve", fontsize=8.5, va="top",
            color=TEXT_SECONDARY)
    ax.set(xscale="log", yscale="log", xlim=(2, 1000), ylim=(.08, 15),
           xlabel="Envelope diameter (m)", ylabel="Share of full lift needed")
    ax.set_xticks([2, 10, 30, 100, 300, 1000])
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.set_yticks([.1, .3, .7, 1, 3, 10])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.legend(loc="upper right", fontsize=8.6, framealpha=.94, edgecolor="none")
    style_axis(ax, "A  Small and large spheres meet different mass limits",
               f"10 km air; sphere; 1 kg living dry tissue per m²\n"
               f"Material {base['material_density_kg_m3']:g} kg/m³; {base['relative_gust_m_s']:g} m/s relative gust")

    ax = axes[0, 1]
    materials = (("kurata_2022_cellophane", "Cellophane, 0.173 Barrer", "-"),
                 ("kurata_2022_polypropylene", "Polypropylene, 6.69 Barrer", "-"),
                 ("mfc_2025_incomplete_conditions", "MFC, 10.7 Barrer*", "--"),
                 ("wu_2014_chitin_incomplete_conditions", "Chitin, 0.024 Barrer*", "--"))
    for i, (name, label, ls) in enumerate(materials):
        points = sorted((r for r in data["envelope"] if r["comparator"] == name), key=lambda r: r["thickness_um"])
        ax.plot([r["thickness_um"] for r in points], [r["hydrogen_chemical_w_m2_projected"] for r in points],
                marker="o", ms=4, lw=2, ls=ls, color=COLORS[i], label=label)
    ax.set(xscale="log", yscale="log", xlim=(1, 1000), xlabel="Intact barrier thickness (µm)",
           ylabel="Hydrogen replacement power (W/m² projected)")
    ax.set_xticks([1, 10, 100, 1000])
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.legend(loc="upper right", fontsize=8.5, framealpha=.94, edgecolor="none")
    style_axis(ax, "B  Thicker barriers trade gas loss for mass",
               "100 kPa total; 98% H₂; full sphere, gas area four times the projection\n"
               "Chemical energy of the lost hydrogen. *Measurement conditions incomplete")

    ax = axes[1, 0]
    for i, leaf in enumerate(sorted(set(r["leaf_c"] for r in data["water"]))):
        points = sorted((r for r in data["water"] if close(r["leaf_c"], leaf)), key=lambda r: r["relative_humidity"])
        delta = leaf - points[0]["air_c"]
        ax.plot([100*r["relative_humidity"] for r in points], [r["water_kg_m2_year"] for r in points],
                "o-", lw=2.2, ms=5, color=COLORS[i], label=f"Leaf {delta:+g} K against the air")
        end = points[-1]
        ax.annotate(f"{end['water_kg_m2_year']:.0f}", (98, end["water_kg_m2_year"]), xytext=(-8, 9),
                    textcoords="offset points", ha="right", fontsize=9, color=TEXT_PRIMARY)
    air = data["water"][0]["air_c"]
    ax.set(xlim=(58, 100), ylim=(0, None), xlabel="Air relative humidity (%)",
           ylabel="Water need (kg/m² projected a year)")
    ax.set_xticks([60, 70, 80, 90, 98])
    ax.legend(loc="upper right", fontsize=9, frameon=False)
    style_axis(ax, "C  Warm leaves still lose water in humid air",
               f"10 km; air {air:g} °C; CO₂ 400 ppm; net leaf carbon 1 kg/m² a year\nStomatal diffusion need")

    ax = axes[1, 1]
    for i, (shape, label) in enumerate((("sphere_sensitivity", "Sphere"), ("streamlined_hypothesis", "Streamlined"))):
        points = [r for r in sorted((r for r in data["motion"] if r["shape"] == shape), key=lambda r: r["airspeed_m_s"])
                  if r["airspeed_m_s"] > 0]
        ax.plot([r["airspeed_m_s"] for r in points], [r["mean_input_w_m2"] for r in points],
                "o-", lw=2.2, ms=4.5, color=COLORS[i], label=f"{label}, C_d {points[0]['drag_coefficient']:g}")
    gpp = data["navigation_assumptions"]["reference_gpp_chemical_equivalent_w_m2"]
    ax.axhline(gpp, ls="--", color=TEXT_MUTED, lw=1.2)
    ax.text(.28, gpp*1.19, "The reference's whole gross fixation", fontsize=8.5, color=TEXT_SECONDARY)
    eta = data["navigation_assumptions"]["biological_overall_efficiency"]
    rho = data["air_10km"]["density_kg_m3"]
    ax.set(xscale="log", yscale="log", xlim=(.25, 5), xlabel="Steady speed through the air (m/s)",
           ylabel="Chemical input for motion (W/m² projected)")
    ax.set_xticks([.25, .5, 1, 2, 5])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.legend(loc="upper left", fontsize=8.6, framealpha=.94, edgecolor="none")
    style_axis(ax, "D  Airspeed sets the cost of motion",
               f"10 km; air {rho:g} kg/m³; overall efficiency {100*eta:g}%\nFrontal area equal to the projection")

    fig.text(.085, .087, "Panels use separate stated conditions. Lines join computed points.", fontsize=9, color=TEXT_SECONDARY)
    fig.text(.085, .036, f"Terluna · floater-viability/2 · product SHA-256 {product_hash[:16]} · plot SHA-256 {script_hash[:16]}",
             fontsize=8.2, color=TEXT_MUTED)
    return fig


def draw_giants(data, product_hash, script_hash):
    fig, axes = plt.subplots(2, 2, figsize=(14, 11.6))
    fig.subplots_adjust(left=.085, right=.965, bottom=.12, top=.82, wspace=.3, hspace=.66)
    fig.text(.085, .955, "Giants of the Open Moon's sky", fontsize=22, weight="bold")
    fig.text(.085, .921, "How large a floater can be, how old a giant is, where rain waters it and how long its day is",
             fontsize=11, color=TEXT_SECONDARY)
    # One colour per body across panels A and B; line style separates growth cases.
    body_colour = {"round_collagen": 0, "round": 1, "round_strong": 2, "colony_module": 3}

    ax = axes[0, 0]
    names = (("round_collagen", "Round, collagen tendons"), ("round", "Round, cellulose fibre"),
             ("round_strong", "Round, 300 MPa fibre"), ("colony_module", "Colony module"))
    for key, label in names:
        i = body_colour[key]
        pts = data["size"][key]
        ax.plot([p["diameter_m"] for p in pts], [p["surplus_c_kg_m2_year"] for p in pts], "-", lw=2,
                color=COLORS[i], label=label)
        ok = [p for p in pts if p["all_gates"]]
        bad = [p for p in pts if not p["all_gates"]]
        ax.plot([p["diameter_m"] for p in ok], [p["surplus_c_kg_m2_year"] for p in ok], "o", ms=5.5,
                color=COLORS[i], mec=SURFACE, mew=1.2, zorder=3)
        ax.plot([p["diameter_m"] for p in bad], [p["surplus_c_kg_m2_year"] for p in bad], "o", ms=5.5,
                mfc=SURFACE, mec=COLORS[i], mew=1.4, zorder=3)
    ax.axhline(0, color=TEXT_MUTED, lw=1.1)
    ax.set(xscale="log", xlim=(18, 5600), ylim=(-3, 1.3), xlabel="Diameter (m; module diameter for the colony)",
           ylabel="Carbon left for growth (kg C/m² a year)")
    ax.set_xticks([20, 50, 100, 200, 500, 1000, 2000, 5000])
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.legend(loc="lower left", fontsize=8.5, framealpha=.94, edgecolor="none")
    style_axis(ax, "A  Tendon renewal ends the round giants' growth",
               "10 km; gross fixation 2 kg C/m² a year; gusts, spare pressure, gas cells and wet skins paid\n"
               "Filled points pass every gate; open points fail one (the lift reserve when small, mostly carbon when large)")

    ax = axes[0, 1]
    ends = {}
    for (label, key, ls), pts in zip((("Colony, photolytic gas, reference light", "colony_module", "-"),
                                      ("Colony, fermented gas, reference light", "colony_module", "--"),
                                      ("Round, 300 MPa fibre, canopy-like light", "round_strong", "-"),
                                      ("Round, 300 MPa fibre, reference light", "round_strong", "--")),
                                     data["ages"].values()):
        if not pts:
            continue
        i = body_colour[key]
        xs = [p["diameter_m"] for p in pts]
        ys = [p["years"] for p in pts]
        ax.plot(xs, ys, ls=ls, lw=2, color=COLORS[i], label=label)
        ax.plot(xs, ys, "o", ms=5, color=COLORS[i], mec=SURFACE, mew=1.2)
        ends.setdefault(xs[-1], []).append(ys[-1])
    last = max(ends)
    for x, ys in ends.items():
        logs = [math.log10(y) for y in ys]
        for y, placed in zip(ys, dodge(logs, .075)):
            inside = x < last  # a series that stops short is labelled above and left of its end
            ax.annotate(f"{y:.0f} years", (x, 10**placed), xytext=(-6, 9) if inside else (7, 0),
                        textcoords="offset points", ha="right" if inside else "left",
                        va="center", fontsize=8.3, color=TEXT_PRIMARY, annotation_clip=False)
    ax.set(xscale="log", yscale="log", xlim=(150, 3600), ylim=(8, 600), xlabel="Diameter (m)",
           ylabel="Years of growth")
    ax.set_xticks([200, 500, 1000, 2000])
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.set_yticks([10, 30, 100, 300])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.legend(loc="lower right", fontsize=8.2, framealpha=.94, edgecolor="none")
    style_axis(ax, "B  Colonies grow fast; round giants grow old",
               f"Round bodies from a {data['start_diameter_m']:g}-m juvenile, colonies from one {data['module_diameter_m']:g}-m module, all surplus to growth\n"
               "Gross fixation 2 or 5 (canopy-like) kg C/m² a year; no point where a body stops short")

    ax = axes[1, 0]
    pts = sorted(data["zonal_rain"], key=lambda p: p["latitude_deg"])
    ax.plot([p["latitude_deg"] for p in pts], [p["rain_mm_day"] for p in pts], "-", lw=2, color=COLORS[0],
            label="Zonal-mean rain, design climate")
    labels = (("reference", "Reference need"), ("concentrating_cool_leaves", "Concentrating, cool leaves"),
              ("concentrating_humid_layer", "Concentrating, humid layer"))
    for i, (key, label) in enumerate(labels, start=1):
        need = data["needs"][key]
        ax.axhline(need, color=COLORS[i], lw=1.6, ls="--", label=f"{label}, {need:.2g} kg/m² a day")
        ax.annotate(f"{need:.2g}", (88, need), xytext=(0, 4), textcoords="offset points", ha="right",
                    va="bottom", fontsize=8.3, color=TEXT_PRIMARY)
    ax.set(xlim=(-90, 90), yscale="log", ylim=(.004, 10), xlabel="Latitude (°)",
           ylabel="Water, mm (kg/m²) a day")
    ax.set_xticks([-90, -60, -30, 0, 30, 60, 90])
    ax.set_yticks([.01, .03, .1, .3, 1, 3, 10])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.legend(loc="lower center", fontsize=8.2, framealpha=.94, edgecolor="none")
    style_axis(ax, "C  Rain waters the tropics; thrift widens the belt",
               "Rain at full catch and the surface rate; needs at 10 km in air of 64% relative humidity\n"
               "Rain meets a need where its line is above it")

    ax = axes[1, 1]
    lats = (0., 30., 45., 60.)
    ends = []
    for i, lat in enumerate(lats):
        pts = sorted((p for p in data["day_length"] if close(p["latitude_deg"], lat)), key=lambda p: p["height_km"])
        xs = [p["height_km"] for p in pts]
        ys = [p["period_days"] for p in pts]
        ax.plot(xs, ys, "-", lw=2, color=COLORS[i], label=f"{lat:g}° latitude")
        ends.append((xs[-1], ys[-1]))
    for i, (lat, (x, y), y_label) in enumerate(zip(lats, ends, dodge([e[1] for e in ends], 1.6))):
        end_label(ax, x, y, f"{lat:g}°", i, y_label)
    fixed = data["fixed_place_days"]
    ax.axhline(fixed, color=TEXT_MUTED, lw=1.1, ls="--")
    ax.text(5.3, fixed + .9, f"A fixed place: one synodic month, {fixed:.1f} days", fontsize=8.5, color=TEXT_SECONDARY)
    freezing = data["freezing_height_km"]
    ax.axvline(freezing, color=TEXT_MUTED, lw=1, ls=":")
    ax.text(freezing + .6, 24, f"Mean freezing\nheight, {freezing:g} km", fontsize=8.5, color=TEXT_SECONDARY, va="top")
    ax.set(xlim=(4, 44), ylim=(0, 33), xlabel="Height riding the mean eastward wind (km)",
           ylabel="Days from noon to noon")
    ax.legend(loc="lower left", fontsize=8.5, framealpha=.94, edgecolor="none")
    style_axis(ax, "D  The eastward wind shortens day and night",
               "Design GCM layer-mean eastward wind, applied at every latitude\n"
               "The longest night is half the period")

    fig.text(.085, .062, "Lines join computed points from the coupled product; each panel has its own axes.", fontsize=9,
             color=TEXT_SECONDARY)
    fig.text(.085, .036, f"Terluna · floater-viability/2 · product SHA-256 {product_hash[:16]} · plot SHA-256 {script_hash[:16]}",
             fontsize=8.2, color=TEXT_MUTED)
    return fig


def render(product_path: Path, output_dir: Path) -> dict:
    product = json.loads(product_path.read_text())
    if product.get("schema") != SCHEMA:
        raise ValueError(f"Expected {SCHEMA}")
    data = {"requirements": select_requirements(product), "giants": select_giants(product)}
    product_hash = digest(product_path)
    script_hash = digest(Path(__file__))
    output_dir.mkdir(parents=True, exist_ok=True)
    setup()
    outputs = []
    for stem, draw, key in (("floater_viability", draw_requirements, "requirements"),
                            ("floater_giants", draw_giants, "giants")):
        fig = draw(data[key], product_hash, script_hash)
        base = output_dir / stem
        fig.savefig(base.with_suffix(".png"), dpi=180, facecolor=SURFACE,
                    metadata={"Software": "Terluna floater_viability plot.py"})
        fig.savefig(base.with_suffix(".svg"), facecolor=SURFACE,
                    metadata={"Date": None, "Creator": "Terluna floater_viability plot.py"})
        plt.close(fig)
        outputs += [base.with_suffix(".png"), base.with_suffix(".svg")]
    plotted_path = output_dir / "plotted_data.json"
    plotted_path.write_text(json.dumps({"schema": "terluna.visualization.floater-viability-points/2",
                                        "product_sha256": product_hash, "points": data}, indent=2) + "\n")
    if json.loads(plotted_path.read_text())["points"] != data:
        raise AssertionError("Plotted-point round trip changed values")
    manifest = {"schema": "terluna.visualization.floater-viability-manifest/2",
                "product": str(product_path.relative_to(ROOT)), "product_sha256": product_hash,
                "plot_source_sha256": script_hash, "producer": product["producer"],
                "outputs": {p.name: digest(p) for p in outputs + [plotted_path]},
                "evidence": product["evidence"],
                "validation": "Exact plotted-point JSON round trip; categorical slots 1-4 validated on the light surface; figures inspected by eye."}
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--product", type=Path, default=DEFAULT_PRODUCT)
    parser.add_argument("--output-dir", type=Path, default=HERE / "outputs")
    args = parser.parse_args()
    result = render(args.product.resolve(), args.output_dir.resolve())
    print(json.dumps({"product_sha256": result["product_sha256"], "outputs": list(result["outputs"])}, indent=2))
