"""Render the enclosed engineering screens without introducing new physics."""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
PRODUCT = ROOT / "research/studies/enclosed_feasibility/results/feasibility.json"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "research/runs/enclosed_feasibility")
    args = parser.parse_args()
    data = json.loads(PRODUCT.read_text())
    case = next(c for c in data["requirements"] if c["name"] == data["reference_case"])
    rows = case["checkpoints"]
    years = [r["year"] for r in rows]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.titleweight": "bold"})
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9.2))
    for key, label, colour in [("ground_firm_W", "Lunar firm electricity", "#267c9c"),
                               ("source_process_W", "Source processing", "#c26936")]:
        axes[0, 0].plot(years, [r[key] / 1e9 for r in rows], "o-", color=colour, label=label)
    axes[0, 0].set(yscale="log", ylabel="Electric power (GW)", xlabel="Programme year", title="Requirements continue beyond the population cap")
    axes[0, 0].legend(loc="upper left", frameon=False)
    axes[0, 1].plot(years, [r["source_feed_kg_year"] / 1e9 for r in rows], "o-", color="#267c9c", label="Operating checkpoints")
    seed = data["prestart"]["seed_average_feed_kg_year"] / 1e9
    axes[0, 1].plot([0, 7], [seed, seed], color="#c26936", lw=4, label="Construction average, years 0–7")
    axes[0, 1].annotate(f"{seed:.1f} Mt/year\nSequence unresolved", (4, seed), xytext=(33, 27),
                       textcoords="offset points", arrowprops={"arrowstyle": "-", "color": "#c26936"}, color="#a24f26")
    axes[0, 1].set(yscale="log", ylabel="Asteroid-source feed (million tonnes/year)", xlabel="Programme year", title="Opening demand does not describe deployment")
    axes[0, 1].legend(loc="lower right", frameon=False)
    beam = data["beam"]["ideal_focus_cases"][0]
    sites = beam["site_comparisons"]
    axes[1, 0].plot([r["delivered_MW"] for r in sites], [r["capture_at_flux_sized_radius"] * 100 for r in sites], "o-", color="#267c9c")
    axes[1, 0].axhline(80, linestyle="--", color="#a24f26", label="Assumed capture: 80%")
    axes[1, 0].set(xscale="log", ylim=(0, 100), ylabel="Ideal captured beam power (%)", xlabel="Delivered site load used to size receiver (MW)", title="Receiver area must satisfy diffraction")
    axes[1, 0].text(0.04, 0.95, "1.2 km aperture · 2.45 GHz · 10,000 km\nReceiver sized at 200 W/m² delivered electricity", transform=axes[1, 0].transAxes, va="top", fontsize=9)
    axes[1, 0].legend(loc="lower right", frameon=False)
    repair = sorted(data["refurbishment"]["cases"], key=lambda r: r["restored_life_fraction"])
    axes[1, 1].plot([r["restored_life_fraction"] * 100 for r in repair], [r["purchase_ratio"] * 100 for r in repair], "o-", color="#267c9c")
    axes[1, 1].axhline(100, color="#a24f26", linestyle="--", label="All-new equipment")
    axes[1, 1].set(xlabel="Service life restored (% of new lifetime)", ylabel="Material purchased per service-year (%)", title="Repair savings require durable rebuilt equipment", ylim=(55, 125))
    axes[1, 1].legend(loc="lower left", frameon=False)
    axes[1, 1].text(0.97, 0.96, "48% accepted for repair\n20% replacement parts\nWorkshops and downtime excluded", transform=axes[1, 1].transAxes, va="top", ha="right", fontsize=9)
    for ax in axes.flat:
        ax.grid(alpha=0.18)
    fig.suptitle("Engineering feasibility of the planned enclosed off-Earth programme", fontsize=16, x=0.06, ha="left")
    fig.text(0.06, 0.035, "Conditional requirements and analytical screens; no integrated system qualification. Lines join checkpoints, not annual simulations.\nScale panels use the strong-market, wider-growth propellant/logistics case. Year 0 = 2041; year 200 = 2241.", fontsize=9, color="#454545")
    fig.tight_layout(rect=(0.025, 0.09, 0.99, 0.95), h_pad=2.3)
    args.out.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf"):
        fig.savefig(args.out / f"engineering-screen.{suffix}", dpi=160)
    provenance = {"schema": "terluna.visualization.enclosed-feasibility/1",
                  "source_product": str(PRODUCT.relative_to(ROOT)),
                  "source_sha256": hashlib.sha256(PRODUCT.read_bytes()).hexdigest(),
                  "renderer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  "reference_case": data["reference_case"]}
    (args.out / "render-provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(args.out / "engineering-screen.png")


if __name__ == "__main__":
    main()
