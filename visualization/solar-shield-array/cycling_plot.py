"""Render the numerical cycling product; no orbit physics lives here."""
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT/"research/studies/solar_shield_array/results/cycling.json"
OUTPUT = Path(__file__).with_name("results")/"cycling.png"


def main():
    product = json.loads(SOURCE.read_text())
    cases = product["cases"]
    primary = cases["inner_filter"]
    first = primary["first_ten_days"]
    color = "#007f79"
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.titleweight": "bold"})
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.8))
    fig.subplots_adjust(top=.84, bottom=.13, wspace=.48, hspace=.42)
    fig.suptitle("Orbital cycling at 50 g/m²", x=.08, y=.975, ha="left", fontsize=21, fontweight="bold")
    fig.text(.08, .927, "Gravity returns the tile; the filter faces the Sun during useful passes and feathers between them.", fontsize=11)
    fig.text(.08, .898, "DE440 point-mass dynamics · ideal spectral mirror · no electric thrust · fleet coverage remains unsolved", fontsize=10, color="#555555")

    ax = axes[0, 0]
    q = np.array(first["solar_coordinates_km"])/1000
    useful = np.array(first["useful"])
    segments = np.stack([q[:-1, :2], q[1:, :2]], axis=1)
    ax.add_collection(LineCollection(segments, colors=np.where(useful[:-1] > .001, color, "#b6bec5"), linewidths=1.5))
    radius = product["producer"]["constants"]["MOON_RADIUS"]/1e6
    ax.add_patch(plt.Circle((0, 0), radius, color="#657582", zorder=3))
    ax.add_patch(plt.Circle((0, 0), 4*radius, fill=False, ls="--", color="#999999", lw=1))
    ax.set(xlim=(-18, 18), ylim=(-18, 18), aspect="equal", xlabel="Sunward position (1,000 km)",
           ylabel="Transverse position (1,000 km)", title="Returning motion · first 10 days")
    ax.annotate("Sunlight", xy=(9, 16), xytext=(16, 16), arrowprops={"arrowstyle": "->", "color": "#bf8b00"},
                color="#8f6800", ha="center", va="center")
    ax.text(.02, .03, "Teal: useful shadow service\nDashed circle: 4 lunar radii\nProjection omits the third coordinate", transform=ax.transAxes, fontsize=8,
            bbox={"facecolor": "white", "alpha": .85, "edgecolor": "none", "pad": 2})

    ax = axes[0, 1]
    long = cases["inner_filter_three_year"]["daily_range_km"]
    x = np.array(long["day"])/365.25
    low, high = np.array(long["minimum"])/1000, np.array(long["maximum"])/1000
    ax.fill_between(x, low, high, color=color, alpha=.2)
    ax.plot(x, low, color=color, lw=.8)
    ax.plot(x, high, color=color, lw=.8)
    ax.set(xlabel="Years after 2026-10-04 TDB", ylabel="Moon-centre distance (1,000 km)",
           title="Longer propagation · daily range")
    ax.grid(axis="y", alpha=.2)
    ax.text(.03, .05, "Shading shows daily minimum–maximum distance.",
            transform=ax.transAxes, fontsize=8, bbox={"facecolor": "white", "alpha": .85, "edgecolor": "none", "pad": 2})

    ax = axes[1, 0]
    ax.plot(first["day"], useful, color=color, lw=1.5)
    ax.fill_between(first["day"], 0, useful, color=color, alpha=.15)
    ax.set(xlim=(0, 10), ylim=(-.03, 1.08), xlabel="Days after epoch", ylabel="Useful projected area / physical area",
           title=f"Annual useful-area fraction: {100*primary['useful_projected_area_fraction']:.2f}%")
    ax.text(.02, .86, "One tile's service passes.\nOther tiles must cover each absence.", transform=ax.transAxes, fontsize=8,
            bbox={"facecolor": "white", "alpha": .85, "edgecolor": "none", "pad": 2})
    ax.grid(axis="y", alpha=.2)

    ax = axes[1, 1]
    labels = ["15,000 km seed\nFilter spectrum", "30,000 km seed\nFilter spectrum", "30,000 km seed\nFull outside core"]
    keys = ["inner_filter", "middle_filter", "middle_full"]
    inventory = product["reference_area_inventory"]
    values = [inventory[k]["base_optical_inventory_kg"]/1e13 for k in keys]
    bars = ax.barh(np.arange(3), values, color=[color, "#6f9eae", "#9ab5be"], height=.58)
    ax.set_yticks(np.arange(3), labels)
    ax.tick_params(axis="y", labelsize=9)
    ax.invert_yaxis()
    ax.set(xlabel="Base optical mass (10¹³ kg)", xlim=(0, max(values)*1.22), title="Equivalent reference-area inventory")
    for bar, key in zip(bars, keys):
        multiplier = inventory[key]["reference_aperture_area_multiplier"]
        ax.text(bar.get_width()+.15, bar.get_y()+bar.get_height()/2, f"{multiplier:.2f}× area", va="center", fontsize=9)
    ax.grid(axis="x", alpha=.2)
    ax.set_axisbelow(True)

    fig.text(.08, .065, "Inventory divides the held aperture's base area by the measured useful-area fraction and assumes perfect placement.", fontsize=9)
    fig.text(.08, .040, "It excludes spatial overlap, seams, collision avoidance, attitude/power hardware, habitats and replacements; it is not a complete fleet design.", fontsize=9)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=170, facecolor="white")
    plt.close(fig)
    metadata = {"source": str(SOURCE.relative_to(ROOT)), "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                "renderer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "image_sha256": hashlib.sha256(OUTPUT.read_bytes()).hexdigest()}
    OUTPUT.with_suffix(".json").write_text(json.dumps(metadata, indent=2)+"\n")
    print(OUTPUT)


if __name__ == "__main__":
    main()
