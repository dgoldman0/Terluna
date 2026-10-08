"""The sea-appearance study's lighting calendar and its map of where the Earth stands over the seas.

    python visualization/sea-appearance/lighting.py

Reads research/studies/sea_appearance/results/lighting_calendar.json and earth_over_seas.npz and the 28% atlas grid.
Writes ignored figures with provenance sidecars to results/.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.patheffects
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "visualization" / "atlas"))
import render as atlas  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"
CALENDAR = ROOT / "research/studies/sea_appearance/results/lighting_calendar.json"
MAP = ROOT / "research/studies/sea_appearance/results/earth_over_seas.npz"
ATLAS_GRID = ROOT / "geography/products/atlas_28pct_4ppd.npz"
ATLAS_JSON = ROOT / "geography/results/atlas.json"
INK, INK2, SURF = atlas.INK, atlas.INK2, atlas.SURF
SUN, EARTH, FLOOR = "#eb6834", "#2a78d6", "#8a8984"          # categorical slots 2 and 1, and a neutral
EARTH_RAMP = plt.get_cmap("cividis")      # dark where the Earth never rises, bright where it always stands
GREY, MESOPIC = 0.18, (0.005, 5.0)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_sidecar(image, inputs, checks):
    record = dict(schema="terluna.visualization.figure/1", image=image.name, image_sha256=digest(image),
                  renderer_sha256=digest(Path(__file__)), inputs={str(Path(p).relative_to(ROOT)): digest(p) for p in inputs},
                  checks=checks)
    image.with_suffix(".json").write_text(json.dumps(record, indent=1) + "\n")


def lux_for_luminance(cd_m2):
    return cd_m2 * math.pi / GREY


def calendar_figure(record):
    sites = record["sites"]
    t = np.array(record["time_hours"])
    days = (t - t[0]) / 24
    fig = plt.figure(figsize=(15.5, 13.2), facecolor=SURF)
    outer = fig.add_gridspec(3, 2, left=0.06, right=0.985, top=0.865, bottom=0.05, hspace=0.42, wspace=0.12)
    photopic, scotopic = lux_for_luminance(MESOPIC[1]), lux_for_luminance(MESOPIC[0])
    checks = {}
    for k, site in enumerate(sites):
        cell = outer[k // 2, k % 2].subgridspec(2, 1, height_ratios=[2.3, 1.0], hspace=0.08)
        ax, ax2 = fig.add_subplot(cell[0]), fig.add_subplot(cell[1])
        s = site["series"]
        sun, earth = np.array(s["ground_lux_from_sun"]), np.array(s["ground_lux_from_earth"])
        floor = np.full_like(sun, 1e-3)
        ax.axhspan(photopic, 1e6, color="#f6f3e8", lw=0, zorder=0)
        ax.axhspan(scotopic, photopic, color="#e8ecf2", lw=0, zorder=0)
        ax.axhspan(1e-5, scotopic, color="#d9d9e4", lw=0, zorder=0)
        ax.plot(days, np.maximum(sun, 1e-5), color=SUN, lw=1.6, label="Sky lit by the Sun")
        ax.plot(days, np.maximum(earth, 1e-5), color=EARTH, lw=1.6, label="The Earth")
        ax.plot(days, floor, color=FLOOR, lw=1.0, ls=(0, (4, 3)), label="Stars and airglow (placeholder)")
        total = sun + earth + 1e-3
        ax.plot(days, total, color=INK, lw=0.9, alpha=0.85, label="All light")
        ax.set_yscale("log"); ax.set_ylim(1e-4, 2e5); ax.set_xlim(0, 29.6)
        ax.set_yticks([1e-4, 1e-2, 1, 1e2, 1e4])
        ax.set_yticklabels(["0.0001", "0.01", "1", "100", "10,000"])
        ax.set_ylabel("Light on the\nground (lux)", fontsize=8, color=INK2)
        ax.set_title(site["name"], fontsize=9.5, color=INK, loc="left")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False); ax2.spines[side].set_visible(False)
        ax.tick_params(labelsize=7, colors=INK2); plt.setp(ax.get_xticklabels(), visible=False)
        if k == 0:
            ax.text(29.4, photopic * 3, "day vision", fontsize=6.8, color=INK2, ha="right", va="bottom")
            ax.text(29.4, math.sqrt(photopic * scotopic), "twilight vision", fontsize=6.8, color=INK2, ha="right", va="center")
            ax.text(29.4, scotopic / 3, "night vision", fontsize=6.8, color=INK2, ha="right", va="top")
        darkest = int(np.argmin(total))
        ax.plot(days[darkest], total[darkest], "o", ms=4, mfc=INK, mec="white", mew=0.8, zorder=5)
        right = days[darkest] < 20
        ax.annotate(f"darkest {total[darkest]:.2g} lux", (days[darkest], total[darkest]),
                    xytext=(8 if right else -8, -11), textcoords="offset points", ha="left" if right else "right",
                    fontsize=7, color=INK, path_effects=[matplotlib.patheffects.withStroke(linewidth=2.2, foreground="white")])
        ax2.axhline(0, color="#b9b8b2", lw=0.7)
        ax2.plot(days, s["sun_elevation_deg"], color=SUN, lw=1.3)
        ax2.plot(days, s["earth_elevation_deg"], color=EARTH, lw=1.3)
        ax2.set_ylim(-95, 95); ax2.set_yticks([-90, 0, 90]); ax2.set_xlim(0, 29.6)
        ax2.set_ylabel("Elevation (°)", fontsize=8, color=INK2)
        ax2.tick_params(labelsize=7, colors=INK2)
        if k >= 4:
            ax2.set_xlabel("Days from 7 February 2038", fontsize=8, color=INK2)
        checks[site["short"]] = dict(darkest_lux=float(total[darkest]), hours_by_mode=site["hours_by_mode"])
    handles = [plt.Line2D([], [], color=INK, lw=0.9), plt.Line2D([], [], color=SUN, lw=1.6),
               plt.Line2D([], [], color=EARTH, lw=1.6), plt.Line2D([], [], color=FLOOR, lw=1.0, ls=(0, (4, 3)))]
    fig.legend(handles, ["All light", "Sky lit by the Sun (clear sky)", "The Earth",
                         "Stars and airglow (placeholder, 0.001 lux)"],
               loc="upper left", bbox_to_anchor=(0.06, 0.915), ncol=4, frameon=False, fontsize=8)
    fig.text(0.06, 0.975, "Light through a month at six coasts", fontsize=13, color=INK, ha="left", va="top")
    note = ("Clear-sky horizontal illuminance hour by hour, 7 February to 8 March 2038, the tide month of the wave "
            "studies' shore month. The Sun's light comes from the solved spherical sky of the design atmosphere, which "
            "scatters light far past the horizon; the Earth's from a model spectrum of the whole Earth at its observed "
            "visual brightness (geometric albedo 0.242), seen through the same air. Bands: modes of vision for an 18% grey surface (CIE 191:2010). Lower panels: the Sun's and "
            "the Earth's elevation.")
    fig.text(0.06, 0.95, textwrap.fill(note, 215), fontsize=7.6, color=INK2, ha="left", va="top")
    path = OUT / "lighting_calendar.png"
    fig.savefig(path, dpi=160, facecolor=SURF)
    plt.close(fig)
    write_sidecar(path, [CALENDAR], checks)
    return path


def map_figure(record):
    with np.load(MAP, allow_pickle=False) as z:
        lat, lon, rows, cols = z["lat_deg"], z["lon_deg"], z["rows"], z["cols"],
        share = z["share_earth_up"].astype(float)
    with np.load(ATLAS_GRID, allow_pickle=False) as z:
        height, alat = z["height_m"].astype(float), z["lat_deg"]
        labels = z["water_label"]
    meta = json.loads(ATLAS_JSON.read_text())
    shade = atlas.shaded(height, alat, meta["grid"]["pixels_per_degree"])
    grid = np.full((len(lat), len(lon)), np.nan)
    grid[rows, cols] = share
    # Upsample the 1-degree shares to the atlas grid inside its water, so coasts keep their quarter-degree shape.
    ppd = meta["grid"]["pixels_per_degree"]
    fine = np.kron(np.nan_to_num(grid, nan=-1), np.ones((ppd, ppd)))
    water = labels > 0
    img = (0.74 + 0.18 * shade)[..., None] * np.array([1.0, 0.985, 0.95])     # light, slightly warm land
    known = water & (fine >= 0)
    img[known] = EARTH_RAMP(fine[known])[:, :3]
    img[water & ~known] = atlas.hexrgb("#d9e0e8")
    # Columns run from 0 E; show the map from 180 W to 180 E with the near side in the middle.
    half = img.shape[1] // 2
    img = np.concatenate([img[:, half:], img[:, :half]], axis=1)
    fig, ax = plt.subplots(figsize=(15.5, 8.7), facecolor=SURF)
    fig.subplots_adjust(left=0.05, right=0.98, top=0.83, bottom=0.14)
    ax.imshow(img, extent=(-180, 180, -90, 90), interpolation="bilinear")
    mask = np.concatenate([water[:, half:], water[:, :half]], axis=1).astype(float)
    ny, nx = mask.shape
    ax.contour(np.linspace(-180, 180, nx), np.linspace(90, -90, ny), mask, levels=[0.5], colors="#2f3237",
               linewidths=0.35)
    ax.set_xticks(range(-180, 181, 30)); ax.set_yticks(range(-90, 91, 30))
    ax.set_xticklabels([f"{abs(x)}°{'W' if x < 0 else 'E' if 0 < x < 180 else ''}" for x in range(-180, 181, 30)])
    ax.set_yticklabels([f"{abs(y)}°{'S' if y < 0 else 'N' if y > 0 else ''}" for y in range(-90, 91, 30)])
    ax.tick_params(labelsize=7.5, colors=INK2)
    for x in (-90, 90):
        ax.axvline(x, color="white", lw=0.8, ls=(0, (3, 3)))
    stroke = [matplotlib.patheffects.withStroke(linewidth=2.4, foreground="white")]
    ax.text(0, 84, "near side: the Earth always up", ha="center", fontsize=8.5, color=INK, path_effects=stroke)
    ax.text(150, 84, "far side: never", ha="center", fontsize=8.5, color=INK, path_effects=stroke)
    ax.text(-150, 84, "far side: never", ha="center", fontsize=8.5, color=INK, path_effects=stroke)
    names = {"Mare Smythii": (87, 2), "Mare Humboldtianum": (82, 57), "Mare Orientale": (-95, -20),
             "Mare Moscoviense": (148, 27), "Mare Ingenii": (164, -34), "Oceanus Procellarum": (-55, 15)}
    for name, (x, y) in names.items():
        ax.annotate(name.replace("Mare ", "").replace("Oceanus ", ""), (x, y), xytext=(6, 6), textcoords="offset points",
                    fontsize=7.5, color=INK, path_effects=[matplotlib.patheffects.withStroke(linewidth=2.2, foreground="white")])
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=mcolors.Normalize(0, 1), cmap=EARTH_RAMP), ax=ax,
                       orientation="horizontal", fraction=0.04, pad=0.08, aspect=50, ticks=[0, .25, .5, .75, 1])
    bar.ax.set_xticklabels(["never", "25%", "50%", "75%", "always"])
    bar.set_label("Share of 2026–2045 with the Earth's centre above the horizon", fontsize=8.5, color=INK2)
    bar.ax.tick_params(labelsize=7.5, colors=INK2); bar.outline.set_visible(False)
    fig.text(0.05, 0.965, "Where the Earth stands over the seas", fontsize=13, color=INK, ha="left", va="top")
    seas = {s["names"][0] if s["names"] else str(s["body"]): s for s in record["earth_over_seas"]["seas"]}
    smythii = seas["Mare Smythii"]["share_of_area_earth_rises_and_sets"]
    note = ("Every water cell at 1-degree spacing, from the Earth's topocentric elevation every six hours through "
            "2026–2045. The Earth stays fixed over the near side and wobbles about 8° with the libration, so the seas "
            f"along the limb see it rise and set through the month: {100 * smythii:.0f}% of Smythii–Marginis, as at its "
            "headland of the wave studies. Dashed lines: the mean limb at 90° W and 90° E.")
    fig.text(0.05, 0.93, textwrap.fill(note, 210), fontsize=7.8, color=INK2, ha="left", va="top")
    path = OUT / "earth_over_seas.png"
    fig.savefig(path, dpi=160, facecolor=SURF)
    plt.close(fig)
    write_sidecar(path, [CALENDAR, MAP, ATLAS_GRID], dict(water_cells=int(len(share)),
                  always=int(np.sum(share >= 1)), never=int(np.sum(share <= 0))))
    return path


def main():
    import matplotlib.patheffects  # noqa: F401
    OUT.mkdir(parents=True, exist_ok=True)
    record = json.loads(CALENDAR.read_text())
    print(calendar_figure(record))
    print(map_figure(record))


if __name__ == "__main__":
    main()
