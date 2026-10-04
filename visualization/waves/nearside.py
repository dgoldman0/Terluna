"""The nearside sea's waves and the seas' monthly tide, as seen from Earth.

    python -m visualization.waves.nearside           # every figure whose product exists
    python -m visualization.waves.nearside tides     # only the tide maps

Reads climate/waves/results/nearside.json, the geography tide product
(results/tides.json and its grid in products/) and the 28% atlas grid. Land is
the atlas's neutral shaded relief; the mapped quantity takes a one-hue orange
ramp, since the atlas's blues mean depth. The wave fields are computed on
1-degree nodes; the maps interpolate between nodes and clip to the atlas's
quarter-degree shoreline. Each image has a provenance sidecar with product and
renderer hashes and round-trip checks of the drawn values.
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
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "visualization" / "atlas"))
import render as atlas  # noqa: E402  (the atlas sheets' relief, projection and labels)

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"
NEARSIDE = ROOT / "climate/waves/results/nearside.json"
SMYTHII = ROOT / "climate/waves/results/cycle.json"
TIDES = ROOT / "geography/results/tides.json"
TIDE_GRID = ROOT / "geography/products/tides_28pct_4ppd.npz"
ATLAS_GRID = ROOT / "geography/products/atlas_28pct_4ppd.npz"
ATLAS_JSON = ROOT / "geography/results/atlas.json"
ORANGE = mcolors.LinearSegmentedColormap.from_list(
    "orange", ["#fdf1ea", "#f9d6c2", "#f4b48f", "#ef8d5e", "#eb6834", "#c94f1d", "#9b3a13", "#6b270b"])
SERIES = ("#2a78d6", "#eb6834", "#1baf7a")      # categorical slots 1-3, in order
OTHER_WATER = "#dde3ea"
INK, INK2, SURF = atlas.INK, atlas.INK2, atlas.SURF
HALO = [pe.withStroke(linewidth=2.4, foreground="white")]
N = 1100


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class Basemap:
    """Neutral relief and water labels on the atlas grid, sampled through the orthographic view."""

    def __init__(self):
        meta = json.loads(ATLAS_JSON.read_text())
        with np.load(ATLAS_GRID, allow_pickle=False) as z:
            self.height, self.lat, self.lon = z["height_m"].astype(float), z["lat_deg"], z["lon_deg"]
            self.labels = z["water_label"]
        self.ppd = meta["grid"]["pixels_per_degree"]
        self.level = meta["sea_level_m"]
        self.shade = atlas.shaded(self.height, self.lat, self.ppd)

    def view(self, lon0, n=N):
        x, y = np.meshgrid(np.linspace(-1, 1, n), np.linspace(1, -1, n))
        r2 = x**2 + y**2
        inside = r2 <= 1
        z = np.sqrt(np.clip(1 - r2, 0, 1))
        lat = np.degrees(np.arcsin(np.clip(y, -1, 1)))
        lon = lon0 + np.degrees(np.arctan2(x, z))
        i, j = atlas.sample_index(lat, lon, self.ppd, self.height.shape)
        img = np.ones((n, n, 3)) * atlas.hexrgb(SURF)
        gray = (0.74 + 0.22 * self.shade[i, j])[..., None] * np.ones(3)
        img[inside] = gray[inside]
        label = np.where(inside, self.labels[i, j], 0)
        return dict(img=img, inside=inside, lat=lat, lon=lon, label=label, rim=inside & (r2 > 0.994))


def fill_grid(values, valid):
    """Extend a field into invalid cells from the nearest valid cell, for interpolation near coasts."""
    _, (r, c) = ndimage.distance_transform_edt(~valid, return_indices=True)
    return values[r, c]


def sample(field, lat0, lon0, step, lat, lon, lon_period=True):
    """Bilinear sample of a regular grid (rows from lat0 north, cols from lon0 east) at points."""
    rows = (lat - lat0) / step
    lon_rel = np.mod(lon - lon0, 360.0) if lon_period else lon - lon0
    cols = lon_rel / step
    return ndimage.map_coordinates(field, [rows.ravel(), cols.ravel()], order=1, mode="nearest").reshape(lat.shape)


def project(lon0, lat, lon):
    return atlas.project(lon0, lat, lon)


def draw_field(ax, base_view, mask, values, cmap, norm):
    img = base_view["img"].copy()
    other = base_view["inside"] & (base_view["label"] > 0) & ~mask
    img[other] = atlas.hexrgb(OTHER_WATER)
    colors = cmap(norm(values[mask]))[:, :3]
    img[mask] = colors
    img[base_view["rim"]] = atlas.hexrgb("#8f8e8a")
    ax.imshow(img, extent=(-1, 1, -1, 1), interpolation="bilinear")
    ax.set_xlim(-1.02, 1.02); ax.set_ylim(-1.02, 1.02); ax.set_aspect("equal"); ax.axis("off")


def sea_labels(ax, lon0, keep=None, color=INK):
    for name, la, lo in atlas.SEA_LABELS:
        if keep is not None and name.replace("\n", " ") not in keep:
            continue
        x, y, vis = project(lon0, la, lo)
        if vis:
            ax.text(x, y, name, ha="center", va="center", fontsize=6.6, color=color, style="italic",
                    path_effects=HALO, zorder=5, linespacing=0.95)


def colorbar(fig, ax, cmap, norm, label, ticks=None):
    mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    bar = fig.colorbar(mappable, ax=ax, orientation="horizontal", fraction=0.045, pad=0.02, aspect=32, ticks=ticks)
    bar.set_label(label, fontsize=8, color=INK2)
    bar.ax.tick_params(labelsize=7, colors=INK2)
    bar.outline.set_visible(False)
    return bar


def write_sidecar(image, inputs, checks):
    record = dict(schema="terluna.visualization.figure/1", image=image.name, image_sha256=digest(image),
                  renderer_sha256=digest(Path(__file__)),
                  inputs={str(Path(p).relative_to(ROOT) if Path(p).is_relative_to(ROOT) else p): digest(p) for p in inputs},
                  checks=checks)
    image.with_suffix(".json").write_text(json.dumps(record, indent=1) + "\n")


# ------------------------------------------------------------------ waves

def wave_grid(product, key):
    """A node field on the sea's 1-degree crop, rows by latitude (ascending), with validity."""
    maps = product["maps"]
    rows, cols = np.array(maps["rows"]), np.array(maps["cols"])
    active = np.array(maps["active"], bool)
    shape = (rows.max() + 2, cols.max() + 2)
    grid, valid = np.zeros(shape), np.zeros(shape, bool)
    grid[rows[active], cols[active]] = np.array(maps[key])[active]
    valid[rows[active], cols[active]] = True
    lat0 = float(np.array(maps["lat"])[0] - rows[0] * 1.0)
    lon0 = float(np.array(maps["lon"])[0] - cols[0] * 1.0)
    return fill_grid(grid, valid), lat0, lon0


def wave_maps(base, product, smythii):
    view = base.view(0.0)
    mask = view["inside"] & (view["label"] == 433)
    fig = plt.figure(figsize=(13.5, 10.6), facecolor=SURF)
    grid = fig.add_gridspec(2, 2, height_ratios=[3.1, 1.25], hspace=0.16, wspace=0.04,
                            left=0.03, right=0.97, top=0.88, bottom=0.09)
    checks = {}
    stats = product["statistics"]
    top = math.ceil(2 * float(np.percentile(np.array(product["maps"]["mean_hs_m"]), 99))) / 2
    panels = [("mean_hs_m", "Mean significant wave height through the cycle", "Mean Hs (m)", 0.0, top),
              ("hours_above_1m", "Share of the cycle with Hs of at least 1 m", "Share of the cycle (%)", 0.0, 100.0)]
    for k, (key, title, label, vmin, vmax) in enumerate(panels):
        ax = fig.add_subplot(grid[0, k])
        field, lat0, lon0 = wave_grid(product, key)
        if key == "hours_above_1m":
            field = 100 * field / product["clock"]["reported_hours"]
        values = sample(field, lat0, lon0 - 0.0, 1.0, view["lat"], np.where(view["lon"] > 180, view["lon"] - 360, view["lon"]),
                        lon_period=False)
        norm = mcolors.Normalize(vmin, vmax)
        draw_field(ax, view, mask, values, ORANGE, norm)
        sea_labels(ax, 0.0, keep={"Oceanus Procellarum", "Mare Imbrium", "Mare Serenitatis", "Mare Crisium",
                                  "Mare Nectaris", "Mare Fecunditatis", "Mare Nubium", "Mare Humorum", "Mare Frigoris"})
        if key == "mean_hs_m":
            maps = product["maps"]
            lon, lat = np.array(maps["lon"]), np.array(maps["lat"])
            direction, strength = np.radians(maps["direction_deg"]), np.array(maps["direction_resultant"])
            pick = (np.mod(np.round(lon - 0.125), 2) == 0) & (np.mod(np.round(lat - 0.875), 2) == 0) & (strength >= 0.3)
            x0, y0, vis = project(0.0, lat[pick], lon[pick])
            # Keep arrows at least 0.055 apart on the page, so foreshortening near the limb does not crowd them.
            kept = []
            for k in np.flatnonzero(vis):
                if all((x0[k] - x0[j])**2 + (y0[k] - y0[j])**2 > 0.055**2 for j in kept):
                    kept.append(k)
            vis = np.zeros_like(vis)
            vis[kept] = True
            d = 0.6
            x1, y1, _ = project(0.0, lat[pick] + d * np.sin(direction[pick]),
                                lon[pick] + d * np.cos(direction[pick]) / np.cos(np.radians(lat[pick])))
            u, v = x1 - x0, y1 - y0
            scale = np.hypot(u, v)
            u, v = u / scale * 0.045 * strength[pick], v / scale * 0.045 * strength[pick]
            ax.quiver(x0[vis], y0[vis], u[vis], v[vis], angles="xy", scale_units="xy", scale=1, color=INK,
                      width=0.0018, headwidth=3.6, headlength=3.8, alpha=0.6, zorder=4)
            peak = stats["peak_lon_lat"]
            px, py, pvis = project(0.0, peak[1], peak[0])
            if pvis:
                ax.plot(px, py, marker="*", ms=11, mec="white", mfc=INK, zorder=6)
                ax.annotate(f"largest Hs {stats['max_hs_m']:.1f} m", (px, py), xytext=(10, 8), textcoords="offset points",
                            fontsize=7.5, color=INK, path_effects=HALO, zorder=6)
        ax.set_title(title, fontsize=10, color=INK, loc="left")
        colorbar(fig, ax, ORANGE, norm, label)
        # Round trip: the drawn field at each node centre equals the node value.
        node_lon, node_lat = np.array(product["maps"]["lon"]), np.array(product["maps"]["lat"])
        drawn = sample(field, lat0, lon0, 1.0, node_lat, node_lon, lon_period=False)
        target = np.array(product["maps"][key], float) * (100 / product["clock"]["reported_hours"] if key == "hours_above_1m" else 1)
        active = np.array(product["maps"]["active"], bool)
        checks[key] = dict(node_round_trip_max_abs=float(np.abs(drawn - target)[active].max()),
                           data_range=[float(target[active].min()), float(target[active].max())], colour_range=[vmin, vmax])
    ax = fig.add_subplot(grid[1, :])
    series = product["series"]
    t = (np.array(series["time_hours"]) - product["clock"]["report_interval_hours"][0]) / 24
    st = smythii["series"]
    ts = (np.array(st["time_hours"]) - smythii["clock"]["report_interval_hours"][0]) / 24
    lines = [(t, series["max_hs_m"], "Nearside sea: largest Hs anywhere", SERIES[0]),
             (t, series["area_mean_hs_m"], "Nearside sea: area mean", SERIES[1]),
             (ts, st["area_mean_hs_m"], "Smythii–Marginis: area mean", SERIES[2])]
    for x, y, name, color in lines:
        ax.plot(x, y, color=color, lw=1.8, label=name)
        ax.text(x[-1] + 0.25, y[-1], name.split(":")[0] if "Smythii" in name else name.split(": ")[1],
                fontsize=7.5, color=INK2, va="center")
    ax.set_xlim(0, 31.5)
    ax.set_ylim(0, max(max(series["max_hs_m"]), 1) * 1.08)
    ax.set_xlabel("Earth days into the second lunar cycle", fontsize=8.5, color=INK2)
    ax.set_ylabel("Significant wave height (m)", fontsize=8.5, color=INK2)
    ax.tick_params(labelsize=7.5, colors=INK2)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(axis="y", color="#e6e5e1", lw=0.6)
    ax.legend(loc="upper left", fontsize=7.5, frameon=False, ncol=3)
    fig.suptitle("Waves of the nearside sea through a lunar cycle, as seen from Earth", fontsize=13, color=INK,
                 x=0.03, ha="left", y=0.975)
    fig.text(0.03, 0.935,
             f"SWAN at lunar gravity on 1-degree nodes, driven by the GCM's surface stress; the first lunar cycle spins the "
             f"sea up and the second is shown. Largest Hs {stats['max_hs_m']:.2f} m; area-and-time mean "
             f"{stats['area_time_mean_hs_m']:.2f} m (Smythii–Marginis {product['smythii_marginis']['area_time_mean_hs_m']:.2f} m). "
             "Arrows: mean direction of travel, length by its steadiness. Maps interpolate between nodes inside the atlas shoreline.",
             fontsize=7.5, color=INK2, ha="left", va="top", wrap=True)
    path = OUT / "nearside_waves.png"
    fig.savefig(path, dpi=170, facecolor=SURF)
    plt.close(fig)
    write_sidecar(path, [NEARSIDE, SMYTHII, ATLAS_GRID], checks)
    return path


def place_labels(ax, anchors, texts, width=0.013, height=0.055):
    """Annotate anchors (data coordinates) at the first nearby offset that keeps labels apart."""
    boxes = []
    offsets = [(0.05, 0.03), (0.05, -0.08), (-0.05, 0.03), (-0.05, -0.08), (0.05, 0.12), (-0.05, 0.12),
               (0.05, -0.17), (-0.05, -0.17), (0.12, 0.2), (-0.12, -0.25)]
    for (x, y), text in zip(anchors, texts):
        w = width * max(len(line) for line in text.split("\n"))
        for dx, dy in offsets:
            left = x + dx if dx > 0 else x + dx - w
            box = (left, y + dy, left + w, y + dy + height)
            if all(box[2] < b[0] or box[0] > b[2] or box[3] < b[1] or box[1] > b[3] for b in boxes):
                break
        boxes.append(box)
        ax.annotate(text, (x, y), xytext=(box[0], box[1]), textcoords="data", fontsize=7.2, color=INK,
                    path_effects=HALO, zorder=7, va="bottom", ha="left",
                    arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6, shrinkA=0, shrinkB=2))


def coast_map(base, product):
    view = base.view(0.0)
    img = view["img"].copy()
    water = view["inside"] & (view["label"] > 0)
    img[water] = atlas.hexrgb(OTHER_WATER)
    maps, coast = product["maps"], product["coast"]
    nodes = np.array(maps["coastal_nodes"])
    power = np.array(maps["coastal_mean_power_W_m"])
    lon, lat = np.radians(np.array(maps["lon"])[nodes]), np.radians(np.array(maps["lat"])[nodes])
    node_unit = np.stack((np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)), axis=-1)
    # The sea's own shoreline band, three pixels wide, coloured by the nearest 1-degree coastal node.
    sea = view["inside"] & (view["label"] == 433)
    band = sea & ndimage.binary_dilation(view["inside"] & ~sea, iterations=3)
    la, lo = np.radians(view["lat"][band]), np.radians(view["lon"][band])
    pixel_unit = np.stack((np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)), axis=-1)
    from scipy.spatial import cKDTree
    distance, nearest = cKDTree(node_unit).query(pixel_unit)
    vmax = float(np.percentile(power, 99))
    norm = mcolors.Normalize(0, vmax)
    near = distance * atlas.MOON_RADIUS < 45e3          # within 1.5 node spacings of a computed coast
    coloured = np.zeros_like(band)
    coloured[band] = near
    img[coloured] = ORANGE(norm(power[nearest[near]]))[:, :3]
    img[view["rim"]] = atlas.hexrgb("#8f8e8a")
    fig, ax = plt.subplots(figsize=(9.6, 10.4), facecolor=SURF)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.86, bottom=0.08)
    ax.imshow(img, extent=(-1, 1, -1, 1), interpolation="bilinear")
    anchors, texts = [], []
    for entry in coast["strongest"][:7]:
        ex, ey, ev = project(0.0, entry["lon_lat"][1], entry["lon_lat"][0])
        if ev:
            anchors.append((float(ex), float(ey)))
            name = entry["near"]
            for plain, display in atlas.DISPLAY.items():
                name = name.replace(plain, display)
            texts.append(f"{name}\n{entry['mean_power_W_m']:,.0f} W/m")
    place_labels(ax, anchors, texts)
    ax.set_xlim(-1.02, 1.02); ax.set_ylim(-1.02, 1.02); ax.set_aspect("equal"); ax.axis("off")
    colorbar(fig, ax, ORANGE, norm, "Mean wave power travelling toward the coast (W per metre of coast)")
    fig.suptitle("Where the nearside sea's waves reach the coast", fontsize=13, color=INK, x=0.03, ha="left", y=0.975)
    note = ("Time-mean shoreward energy transport through the second lunar cycle, computed at 1-degree sea nodes "
            "beside land; the atlas shoreline takes the value of its nearest coastal node within 45 km, and inlets "
            "farther from a node stay uncoloured. Median coast "
            f"{coast['mean_power_W_m']['median']:.0f} W/m, 90th percentile {coast['mean_power_W_m']['p90']:.0f} W/m, "
            f"strongest {coast['mean_power_W_m']['max']:,.0f} W/m; colours stop at the 99th percentile.")
    fig.text(0.03, 0.935, textwrap.fill(note, 130), fontsize=7.8, color=INK2, ha="left", va="top")
    path = OUT / "nearside_coasts.png"
    fig.savefig(path, dpi=170, facecolor=SURF)
    plt.close(fig)
    write_sidecar(path, [NEARSIDE, ATLAS_GRID], dict(
        coastal_nodes=int(len(power)), shoreline_pixels=int(band.sum()), coloured_shoreline_pixels=int(near.sum()),
        coloured_to_node_distance_max_km=float(distance[near].max() * atlas.MOON_RADIUS / 1000),
        colour_max_W_m=vmax, data_max_W_m=float(power.max())))
    return path


# ------------------------------------------------------------------ tides

def tide_maps(base):
    product = json.loads(TIDES.read_text())
    with np.load(TIDE_GRID, allow_pickle=False) as z:
        rows, cols = z["rows"], z["cols"]
        monthly = z["monthly_range_median_m"].astype(float)
        lat, lon = z["lat_deg"], z["lon_deg"]
    grid, valid = np.zeros((len(lat), len(lon))), np.zeros((len(lat), len(lon)), bool)
    grid[rows, cols], valid[rows, cols] = monthly, True
    field = fill_grid(grid, valid)
    fig, axes = plt.subplots(1, 3, figsize=(18, 8.0), facecolor=SURF)
    fig.subplots_adjust(left=0.015, right=0.985, top=0.82, bottom=0.1, wspace=0.03)
    norm = mcolors.Normalize(0, 6.0)
    checks = {}
    views = (("Near side, facing Earth (0°)", 0.0), ("Eastern limb (90° E): Smythii–Marginis, Humboldtianum", 90.0),
             ("Far side (180°)", 180.0))
    for ax, (title, lon0) in zip(axes, views):
        view = base.view(lon0)
        mask = view["inside"] & (view["label"] > 0) & np.isin(view["label"], list(map(int, product["seas"])))
        values = sample(field[::-1], float(lat[-1]), float(lon[0]), 0.25, view["lat"], view["lon"])
        draw_field(ax, view, mask, values, ORANGE, norm)
        levels = np.arange(1, 7)
        shown = np.where(mask, values, np.nan)
        x = np.linspace(-1, 1, N)
        cs = ax.contour(x, x[::-1], shown, levels=levels, colors=INK, linewidths=0.7, alpha=0.85)
        ax.clabel(cs, fmt=lambda v: f"{v:g} m", fontsize=7, colors=INK)
        for station in product["stations"]:
            headland = station["name"].startswith("Eastern")
            if headland and lon0 != 90.0:
                continue
            sx, sy, sv = project(lon0, station["lon_lat"][1], station["lon_lat"][0])
            if sv and abs(math.cos(math.radians(station["lon_lat"][0] - lon0))) * math.cos(math.radians(station["lon_lat"][1])) > 0.2:
                name = "Smythii headland (wave studies)" if headland else station["name"].replace("Mare ", "").replace("Oceanus ", "")
                ax.plot(sx, sy, "o", ms=3.5, mfc=INK, mec="white", mew=0.8, zorder=6)
                ax.annotate(f"{name} {station['monthly_range_median_m']:.{2 if station['monthly_range_median_m'] < 1 else 1}f} m",
                            (sx, sy), xytext=(5, 4), textcoords="offset points", fontsize=7, color=INK,
                            path_effects=HALO, zorder=6)
        ax.set_title(title, fontsize=10, color=INK, loc="left")
        checks[title] = dict(drawn_max_m=float(np.nanmax(shown)))
    # Round trip: the drawn field at each station equals the station's typical monthly range.
    st = product["stations"]
    drawn = sample(field[::-1], float(lat[-1]), float(lon[0]), 0.25, np.array([s["lon_lat"][1] for s in st]),
                   np.array([s["lon_lat"][0] for s in st]))
    checks["station_round_trip_max_abs_m"] = float(np.max(np.abs(drawn - np.array([s["monthly_range_median_m"] for s in st]))))
    colorbar(fig, axes, ORANGE, norm, "Typical monthly tide range (m): highest minus lowest water in an anomalistic month, median over 19 years")
    fig.suptitle("The monthly tide of the Open Moon's seas", fontsize=13, color=INK, x=0.02, ha="left", y=0.975)
    seas = product["seas"]
    note = ("Equilibrium tide from the Earth and the Sun, 2026–2045, with each sea's volume held and the water's "
            "self-attraction. Typical monthly range over each sea's area: nearside "
            f"{seas['433']['monthly_range_median_m']['area_median']:.1f} m, South Pole–Aitken "
            f"{seas['1417']['monthly_range_median_m']['area_median']:.1f} m, Smythii–Marginis "
            f"{seas['874']['monthly_range_median_m']['area_median']:.2f} m. Inside Mare Fecunditatis the fortnightly "
            "lines run at about twice this equilibrium (geography/README.md). Contours every metre.")
    note = note.replace(" m,", "\u00a0m,").replace(" m.", "\u00a0m.")   # keep values with their units
    fig.text(0.02, 0.925, textwrap.fill(note, 175), fontsize=8, color=INK2, ha="left", va="top")
    path = OUT / "tides.png"
    fig.savefig(path, dpi=170, facecolor=SURF)
    plt.close(fig)
    write_sidecar(path, [TIDES, TIDE_GRID, ATLAS_GRID], checks)
    return path


def main():
    which = sys.argv[1:] or ["tides", "waves"]
    OUT.mkdir(parents=True, exist_ok=True)
    base = Basemap()
    if "tides" in which and TIDE_GRID.exists():
        print(tide_maps(base))
    if "waves" in which and NEARSIDE.exists():
        product = json.loads(NEARSIDE.read_text())
        print(wave_maps(base, product, json.loads(SMYTHII.read_text())))
        print(coast_map(base, product))


if __name__ == "__main__":
    main()
