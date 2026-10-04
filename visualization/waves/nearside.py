"""The nearside sea's waves and coasts, a month at four of its shores, and the seas' monthly tide.

    python -m visualization.waves.nearside           # every figure whose product exists
    python -m visualization.waves.nearside tides     # only the tide maps

Reads climate/waves/results/nearside.json and shore_month.json, the geography
tide product (results/tides.json and its grid in products/) and the 28% atlas
grid. The wave and coast maps use a Lambert azimuthal equal-area projection
centred on the nearside sea, so its western coasts keep their width; the tide
maps show the whole Moon in three orthographic views. Land is the atlas's
neutral shaded relief, and coastlines are the half-level contour of the
bilinearly interpolated water mask. Heights take a one-hue orange ramp and
periods a violet one, since the atlas's blues mean depth. The wave fields are
computed on 1-degree nodes; the maps interpolate between nodes and clip to the
atlas's quarter-degree shoreline. Each image has a provenance sidecar with
product and renderer hashes and round-trip checks of the drawn values.
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
import render as atlas  # noqa: E402  (the atlas sheets' relief, labels and palette)
from shared.constants import MOON_SURFACE_GRAVITY  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"
NEARSIDE = ROOT / "climate/waves/results/nearside.json"
SHORE_MONTH = ROOT / "climate/waves/results/shore_month.json"
SMYTHII = ROOT / "climate/waves/results/cycle.json"
TIDES = ROOT / "geography/results/tides.json"
TIDE_GRID = ROOT / "geography/products/tides_28pct_4ppd.npz"
ATLAS_GRID = ROOT / "geography/products/atlas_28pct_4ppd.npz"
ATLAS_JSON = ROOT / "geography/results/atlas.json"
ORANGE = mcolors.LinearSegmentedColormap.from_list(
    "orange", ["#fdf1ea", "#f9d6c2", "#f4b48f", "#ef8d5e", "#eb6834", "#c94f1d", "#9b3a13", "#6b270b"])
VIOLET = mcolors.LinearSegmentedColormap.from_list(
    "violet", ["#f4f0fa", "#ddd3ef", "#c2b1e1", "#a48fd1", "#866dbf", "#6a4fa8", "#4f348b", "#351c66"])
SERIES = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100")      # categorical slots 1-4, in order
OTHER_WATER = "#dde3ea"
SEA_FLAT = "#c4cfdb"                                       # the coast map's sea: water without a measure on it
COAST = "#56606b"
INK, INK2, SURF = atlas.INK, atlas.INK2, atlas.SURF
HALO = [pe.withStroke(linewidth=2.4, foreground="white")]
SEA = 433
CENTRE = (-6.0, 19.0)                                      # lon, lat at the middle of the nearside sea
NEARSIDE_MARIA = {"Oceanus Procellarum", "Mare Imbrium", "Mare Serenitatis", "Mare Crisium", "Mare Nectaris",
                  "Mare Fecunditatis", "Mare Nubium", "Mare Humorum", "Mare Frigoris", "Mare Cognitum"}
# The coast table of climate/waves/nearside.md, by the product's nearest landmark.
COAST_NAMES = {"Russell": "Western Oceanus Procellarum, by Russell",
               "Harding": "North-western Oceanus Procellarum, by Harding",
               "Cardanus": "Western Oceanus Procellarum, near 12° N",
               "Delisle": "South-western Mare Imbrium, by Delisle",
               "Gay Lussac": "Southern Mare Imbrium, at the eastern end of Montes Carpatus",
               "Atwood": "Eastern Mare Fecunditatis, by Atwood",
               "Jansen": "The Tranquillitatis plateau's north-western shore, by Jansen",
               "Rupes Toscanelli": "The Aristarchus plateau's north-eastern shore, by Rupes Toscanelli"}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def hexrgb(color):
    return atlas.hexrgb(color)


# ------------------------------------------------------------------ projections

def azimuthal_inverse(x, y, c, lon0, s0, c0):
    """Latitude and longitude (degrees, longitude in [-180, 180)) at angular distance c from the centre."""
    rho = np.hypot(x, y)
    sinc, cosc = np.sin(c), np.cos(c)
    safe = np.where(rho > 0, rho, 1.0)
    lat = np.arcsin(np.clip(cosc * s0 + np.where(rho > 0, y * sinc * c0 / safe, 0.0), -1, 1))
    lon = lon0 + np.degrees(np.arctan2(x * sinc, rho * c0 * cosc - y * s0 * sinc))
    return np.degrees(lat), (lon + 180.0) % 360.0 - 180.0


class Ortho:
    """The Moon seen from far away above (lon0, lat0); x east, y north, radius 1."""

    def __init__(self, lon0, lat0=0.0):
        self.lon0, self.s0, self.c0 = lon0, math.sin(math.radians(lat0)), math.cos(math.radians(lat0))

    def forward(self, lat, lon):
        p, dl = np.radians(lat), np.radians(np.asarray(lon, float) - self.lon0)
        cosc = self.s0 * np.sin(p) + self.c0 * np.cos(p) * np.cos(dl)
        return np.cos(p) * np.sin(dl), self.c0 * np.sin(p) - self.s0 * np.cos(p) * np.cos(dl), cosc > 0.05

    def inverse(self, x, y):
        rho = np.hypot(x, y)
        lat, lon = azimuthal_inverse(x, y, np.arcsin(np.clip(rho, 0, 1)), self.lon0, self.s0, self.c0)
        return lat, lon, rho <= 1


class Laea:
    """Lambert azimuthal equal-area projection about (lon0, lat0); x east, y north, unit sphere."""

    def __init__(self, lon0, lat0):
        self.lon0, self.s0, self.c0 = lon0, math.sin(math.radians(lat0)), math.cos(math.radians(lat0))

    def forward(self, lat, lon):
        p, dl = np.radians(lat), np.radians(np.asarray(lon, float) - self.lon0)
        cosc = self.s0 * np.sin(p) + self.c0 * np.cos(p) * np.cos(dl)
        k = np.sqrt(2 / np.maximum(1 + cosc, 1e-12))
        return k * np.cos(p) * np.sin(dl), k * (self.c0 * np.sin(p) - self.s0 * np.cos(p) * np.cos(dl)), cosc > -0.9

    def inverse(self, x, y):
        rho = np.hypot(x, y)
        lat, lon = azimuthal_inverse(x, y, 2 * np.arcsin(np.clip(rho / 2, 0, 1)), self.lon0, self.s0, self.c0)
        return lat, lon, rho <= 2


def unit_vectors(lat, lon):
    p, l = np.radians(lat), np.radians(lon)
    return np.stack((np.cos(p) * np.cos(l), np.cos(p) * np.sin(l), np.sin(p)), axis=-1)


# ------------------------------------------------------------------ the base map

class Basemap:
    """Neutral relief, water labels and smooth water masks on the atlas grid, sampled through a projection."""

    def __init__(self):
        meta = json.loads(ATLAS_JSON.read_text())
        with np.load(ATLAS_GRID, allow_pickle=False) as z:
            self.height, self.lat, self.lon = z["height_m"].astype(float), z["lat_deg"], z["lon_deg"]
            self.labels = z["water_label"]
        self.ppd = meta["grid"]["pixels_per_degree"]
        self.level = meta["sea_level_m"]
        self.shade = atlas.shaded(self.height, self.lat, self.ppd)
        self._shade = self._ready(self.shade)

    def _ready(self, grid):
        """Columns from 180 W, with two columns of wrap on each side, for bilinear sampling."""
        rolled = np.roll(grid, int(180 * self.ppd), axis=1).astype(np.float32)
        return np.concatenate([rolled[:, -2:], rolled, rolled[:, :2]], axis=1)

    def bilinear(self, ready, lat, lon):
        rows = (90 - lat) * self.ppd - 0.5
        cols = (lon + 180) * self.ppd - 0.5 + 2
        return ndimage.map_coordinates(ready, [rows.ravel(), cols.ravel()], order=1, mode="nearest").reshape(lat.shape)

    def plane(self, proj, extent, nx):
        x0, x1, y0, y1 = extent
        step = (x1 - x0) / nx
        ny = int(round((y1 - y0) / step))
        x = x0 + step * (np.arange(nx) + 0.5)
        y = y1 - step * (np.arange(ny) + 0.5)
        X, Y = np.meshgrid(x, y)
        lat, lon, inside = proj.inverse(X, Y)
        i, j = atlas.sample_index(lat, lon, self.ppd, self.labels.shape)
        return dict(proj=proj, x=x, y=y, X=X, Y=Y, lat=lat, lon=lon, inside=inside, extent=(x0, x1, y0, y1),
                    relief=self.bilinear(self._shade, lat, lon), label=np.where(inside, self.labels[i, j], 0))

    def soft(self, view, bodies):
        """The water mask of the given bodies, bilinearly interpolated: its 0.5 contour is the coastline."""
        return self.bilinear(self._ready(np.isin(self.labels, list(bodies))), view["lat"], view["lon"])


def compose(view, soft, values=None, cmap=None, norm=None, bodies=(SEA,), land=(0.74, 0.22), sea_color=None,
            other_water=OTHER_WATER):
    """The relief with other water flat, and the given bodies coloured by `values` (or a flat colour)."""
    img = np.repeat((land[0] + land[1] * view["relief"])[..., None], 3, axis=2)
    other = view["inside"] & (view["label"] > 0) & ~np.isin(view["label"], list(bodies))
    img[other] = hexrgb(other_water)
    sea = view["inside"] & (soft >= 0.5)
    img[sea] = hexrgb(sea_color) if values is None else cmap(norm(values[sea]))[:, :3]
    img[~view["inside"]] = hexrgb(SURF)
    return img, sea


def show(ax, view, img):
    ax.imshow(img, extent=view["extent"], interpolation="bilinear", zorder=1)
    ax.set_xlim(*view["extent"][:2]); ax.set_ylim(*view["extent"][2:]); ax.set_aspect("equal"); ax.axis("off")


def coastline(ax, view, soft, color=COAST, lw=0.5, zorder=3):
    ax.contour(view["X"], view["Y"], np.where(view["inside"], soft, 0.0), levels=[0.5], colors=color,
               linewidths=lw, zorder=zorder)


def sea_extent(base, proj, body=SEA, margin=0.025):
    rows, cols = np.nonzero(base.labels == body)
    x, y, _ = proj.forward(base.lat[rows], base.lon[cols])
    dx, dy = x.max() - x.min(), y.max() - y.min()
    return (x.min() - margin * dx, x.max() + margin * dx, y.min() - margin * dy, y.max() + margin * dy)


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


def sea_labels(ax, view, keep=None, size=6.6):
    x0, x1, y0, y1 = view["extent"]
    for name, la, lo in atlas.SEA_LABELS:
        if keep is not None and name.replace("\n", " ") not in keep:
            continue
        x, y, vis = view["proj"].forward(la, lo)
        if vis and x0 < x < x1 and y0 < y < y1:
            ax.text(x, y, name, ha="center", va="center", fontsize=size, color=INK, style="italic",
                    path_effects=HALO, zorder=6, linespacing=0.95)


def colorbar(fig, ax, cmap, norm, label, ticks=None, **kw):
    mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    options = dict(orientation="horizontal", fraction=0.05, pad=0.02, aspect=34)
    options.update(kw)
    bar = fig.colorbar(mappable, ax=ax, ticks=ticks, **options)
    bar.set_label(label, fontsize=8, color=INK2)
    bar.ax.tick_params(labelsize=7, colors=INK2)
    bar.outline.set_visible(False)
    return bar


def wavelength(period):
    return MOON_SURFACE_GRAVITY * np.asarray(period, float)**2 / (2 * math.pi)


def period_of(length):
    return np.sqrt(2 * math.pi * np.maximum(np.asarray(length, float), 0) / MOON_SURFACE_GRAVITY)


def write_sidecar(image, inputs, checks):
    record = dict(schema="terluna.visualization.figure/1", image=image.name, image_sha256=digest(image),
                  renderer_sha256=digest(Path(__file__)),
                  inputs={str(Path(p).relative_to(ROOT) if Path(p).is_relative_to(ROOT) else p): digest(p) for p in inputs},
                  checks=checks)
    image.with_suffix(".json").write_text(json.dumps(record, indent=1) + "\n")


def note(fig, x, y, text, width, size=7.8):
    fig.text(x, y, textwrap.fill(text, width), fontsize=size, color=INK2, ha="left", va="top")


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


def direction_arrows(ax, view, product, spacing=0.105, length=0.07):
    """Mean direction of travel at steady nodes, thinned to an even spacing on the page; length by steadiness."""
    maps = product["maps"]
    lon, lat = np.array(maps["lon"]), np.array(maps["lat"])
    direction, strength = np.radians(maps["direction_deg"]), np.array(maps["direction_resultant"])
    active = np.array(maps["active"], bool)
    proj = view["proj"]
    x0, y0, vis = proj.forward(lat, lon)
    candidates = np.flatnonzero(active & vis & (strength >= 0.3))
    candidates = candidates[np.argsort(-strength[candidates])]
    kept = []
    for k in candidates:
        if all((x0[k] - x0[j])**2 + (y0[k] - y0[j])**2 > spacing**2 for j in kept):
            kept.append(k)
    kept = np.array(kept)
    d = 0.5
    x1, y1, _ = proj.forward(lat[kept] + d * np.sin(direction[kept]),
                             lon[kept] + d * np.cos(direction[kept]) / np.cos(np.radians(lat[kept])))
    u, v = x1 - x0[kept], y1 - y0[kept]
    scale = np.hypot(u, v)
    u, v = u / scale * length * strength[kept], v / scale * length * strength[kept]
    ax.quiver(x0[kept] - u / 2, y0[kept] - v / 2, u, v, angles="xy", scale_units="xy", scale=1, color="#2b2a28",
              alpha=0.78, width=0.0026, headwidth=3.8, headlength=4.0, headaxislength=3.6, zorder=5)
    return len(kept)


def wave_maps(base, product, smythii):
    proj = Laea(*CENTRE)
    extent = sea_extent(base, proj)
    view = base.plane(proj, extent, nx=900)
    soft = base.soft(view, [SEA])
    width, margin, gap = 15.6, 0.3, 0.16
    map_w = (width - 2 * margin - 2 * gap) / 3
    map_h = map_w * (extent[3] - extent[2]) / (extent[1] - extent[0])
    series_h, below = 2.0, 0.62
    height = 1.4 + map_h + 1.15 + series_h + below + 0.35
    fig = plt.figure(figsize=(width, height), facecolor=SURF)
    box = lambda left, bottom, w, h: fig.add_axes([left / width, bottom / height, w / width, h / height])
    map_bottom = height - 1.4 - map_h
    checks = {}
    stats = product["statistics"]
    maps = product["maps"]
    active = np.array(maps["active"], bool)
    lon_w = np.where(view["lon"] > 180, view["lon"] - 360, view["lon"])
    height_norm = mcolors.Normalize(0.0, 4.0)
    period_norm = mcolors.Normalize(8.0, 17.0)
    panels = [("mean_hs_m", "Mean significant wave height", ORANGE, height_norm),
              ("p95_hs_m", "Height exceeded 5% of the time", ORANGE, height_norm),
              ("energy_mean_tm01_s", "Mean period of the wave energy", VIOLET, period_norm)]
    axes = []
    for k, (key, title, cmap, norm) in enumerate(panels):
        ax = box(margin + k * (map_w + gap), map_bottom, map_w, map_h)
        axes.append(ax)
        field, lat0, lon0 = wave_grid(product, key)
        values = sample(field, lat0, lon0, 1.0, view["lat"], lon_w, lon_period=False)
        img, _ = compose(view, soft, values, cmap, norm)
        show(ax, view, img)
        coastline(ax, view, soft, lw=0.45)
        sea_labels(ax, view, keep=NEARSIDE_MARIA, size=6.2)
        if key == "mean_hs_m":
            checks["direction_arrows"] = direction_arrows(ax, view, product)
        if key == "p95_hs_m":
            peak = stats["peak_lon_lat"]
            px, py, _ = proj.forward(peak[1], peak[0])
            ax.plot(px, py, marker="*", ms=10, mec="white", mfc=INK, zorder=7)
            ax.annotate(f"largest Hs {stats['max_hs_m']:.1f} m", (px, py), xytext=(8, 6), textcoords="offset points",
                        fontsize=7.2, color=INK, path_effects=HALO, zorder=7)
        ax.set_title(title, fontsize=10, color=INK, loc="left")
        node_lon, node_lat = np.array(maps["lon"]), np.array(maps["lat"])
        drawn = sample(field, lat0, lon0, 1.0, node_lat, node_lon, lon_period=False)
        target = np.array(maps[key], float)
        checks[key] = dict(node_round_trip_max_abs=float(np.abs(drawn - target)[active].max()),
                           data_range=[float(target[active].min()), float(target[active].max())],
                           colour_range=[norm.vmin, norm.vmax])
    bar_bottom = map_bottom - 0.62
    bars = [(box(margin + 0.4, bar_bottom, 2 * map_w + gap - 0.8, 0.13), ORANGE, height_norm,
             "Significant wave height (m)", np.arange(0, 4.5, 0.5), "max"),
            (box(margin + 2 * (map_w + gap) + 0.4, bar_bottom, map_w - 0.8, 0.13), VIOLET, period_norm,
             "Mean period Tm01, weighted by wave energy (s)", np.arange(8, 18, 1), "both")]
    for cax, cmap, norm, label, ticks, extend in bars:
        bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal",
                           ticks=ticks, extend=extend)
        bar.set_label(label, fontsize=8, color=INK2)
        bar.ax.tick_params(labelsize=7, colors=INK2)
        bar.outline.set_visible(False)
    top = bar.ax.secondary_xaxis("top", functions=(wavelength, period_of))
    top.set_ticks([20, 30, 40, 50, 60, 70])
    top.tick_params(labelsize=7, colors=INK2, length=2.5, pad=1)
    bar.ax.text(-0.01, 1.9, "wavelength (m)", transform=bar.ax.transAxes, fontsize=7, color=INK2, ha="right", va="center")
    ax = box(0.95, below, width - 0.95 - 1.75, series_h)
    series = product["series"]
    t = (np.array(series["time_hours"]) - product["clock"]["report_interval_hours"][0]) / 24
    st = smythii["series"]
    ts = (np.array(st["time_hours"]) - smythii["clock"]["report_interval_hours"][0]) / 24
    lines = [(t, series["max_hs_m"], "Nearside sea: largest Hs anywhere", "largest anywhere", SERIES[0]),
             (t, series["area_mean_hs_m"], "Nearside sea: mean over its area", "nearside mean", SERIES[1]),
             (ts, st["area_mean_hs_m"], "Smythii–Marginis: mean over its area", "Smythii–Marginis mean", SERIES[2])]
    for x, y, name, short, color in lines:
        ax.plot(x, y, color=color, lw=1.6, label=name)
        ax.plot([x[-1] + 0.2, x[-1] + 0.55], [y[-1], y[-1]], color=color, lw=1.6)
        ax.text(x[-1] + 0.7, y[-1], short, fontsize=7.2, color=INK2, va="center")
    ax.set_xlim(0, 32.5)
    ax.set_ylim(0, max(max(series["max_hs_m"]), 1) * 1.1)
    ax.set_xticks(np.arange(0, 30, 5))
    ax.set_xlabel("Earth days into the second lunar cycle", fontsize=8.5, color=INK2)
    ax.set_ylabel("Significant wave\nheight (m)", fontsize=8.5, color=INK2)
    ax.tick_params(labelsize=7.5, colors=INK2)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(axis="y", color="#e6e5e1", lw=0.6)
    ax.legend(loc="lower left", fontsize=7.4, frameon=False, ncol=3, bbox_to_anchor=(0, 1.0))
    fig.text(margin / width, 1 - 0.32 / height, "Waves of the nearside sea through a lunar cycle", fontsize=13,
             color=INK, ha="left", va="top")
    note(fig, margin / width, 1 - 0.62 / height,
         f"SWAN at lunar gravity on 1-degree nodes, driven by the GCM's surface stress; the first lunar cycle spins the "
         f"sea up and the second is shown. Area-and-time mean Hs {stats['area_time_mean_hs_m']:.2f} m "
         f"(Smythii–Marginis {product['smythii_marginis']['area_time_mean_hs_m']:.2f} m); largest "
         f"{stats['max_hs_m']:.2f} m. Arrows: the energy-weighted mean direction of travel, drawn where it holds "
         "(length by its steadiness, none below 0.3). Equal-area maps centred on the sea; fields interpolate between "
         "nodes inside the atlas shoreline.", 205, size=7.6)
    path = OUT / "nearside_waves.png"
    fig.savefig(path, dpi=170, facecolor=SURF)
    plt.close(fig)
    write_sidecar(path, [NEARSIDE, SMYTHII, ATLAS_GRID], checks)
    return path


def coast_map(base, product):
    proj = Laea(*CENTRE)
    extent = sea_extent(base, proj)
    map_w, side_w, margin = 11.2, 3.7, 0.25
    map_h = map_w * (extent[3] - extent[2]) / (extent[1] - extent[0])
    width, height = margin + map_w + 0.2 + side_w + margin, 1.15 + map_h + 0.95
    view = base.plane(proj, extent, nx=1800)
    soft = base.soft(view, [SEA])
    img, _ = compose(view, soft, land=(0.80, 0.17), sea_color=SEA_FLAT, other_water="#d9e0e8")
    maps, coast = product["maps"], product["coast"]
    nodes = np.array(maps["coastal_nodes"])
    power = np.array(maps["coastal_mean_power_W_m"])
    lat, lon = np.array(maps["lat"])[nodes], np.array(maps["lon"])[nodes]
    fig = plt.figure(figsize=(width, height), facecolor=SURF)
    box = lambda left, bottom, w, h: fig.add_axes([left / width, bottom / height, w / width, h / height])
    ax = box(margin, 0.95, map_w, map_h)
    show(ax, view, img)
    coastline(ax, view, soft, lw=0.45)
    # Each coastal node as a dot coloured and sized by its mean power, the strongest drawn last.
    norm = mcolors.LogNorm(5.0, 800.0)
    x, y, _ = proj.forward(lat, lon)
    order = np.argsort(power)
    size = 5 + 55 * np.sqrt(np.clip(power, 0, None) / power.max())
    ax.scatter(x[order], y[order], s=size[order], c=np.clip(power[order], norm.vmin, None), cmap=ORANGE, norm=norm,
               edgecolors="#3c3a36", linewidths=0.25, zorder=4)
    strongest = coast["strongest"][:8]
    listing = []
    sea_mask = base._ready(base.labels == SEA)
    for k, entry in enumerate(strongest, 1):
        sx, sy, _ = proj.forward(entry["lon_lat"][1], entry["lon_lat"][0])
        # The number goes on the land side: the direction in which the sea mask, 40 km out, is lowest.
        angles = np.radians(np.arange(0, 360, 30))
        probe_x, probe_y = sx + 0.023 * np.cos(angles), sy + 0.023 * np.sin(angles)
        probe_lat, probe_lon, _ = proj.inverse(probe_x, probe_y)
        land = angles[np.argmin(base.bilinear(sea_mask, probe_lat, probe_lon))]
        ax.annotate(str(k), (sx, sy), xytext=(15 * np.cos(land), 15 * np.sin(land)), textcoords="offset points",
                    fontsize=7.5, color=INK,
                    weight="bold", ha="center", va="center", zorder=7,
                    bbox=dict(boxstyle="circle,pad=0.25", fc="white", ec=INK, lw=0.7),
                    arrowprops=dict(arrowstyle="-", color=INK, lw=0.6, shrinkA=0, shrinkB=3))
        listing.append((k, COAST_NAMES.get(entry["near"], entry["near"]), entry["mean_power_W_m"]))
    sea_labels(ax, view, keep=NEARSIDE_MARIA, size=7.0)
    side = box(margin + map_w + 0.2, 0.95, side_w, map_h)
    side.set_xlim(0, 1); side.set_ylim(0, 1); side.axis("off")
    side.text(0, 0.98, "The strongest coasts", fontsize=10, color=INK, va="top")
    side.text(0, 0.945, "Mean power toward the coast; each at least\n300 km from a stronger one",
              fontsize=7.6, color=INK2, va="top", linespacing=1.3)
    for row, (k, name, value) in enumerate(listing):
        top = 0.85 - row * 0.094
        side.text(0.035, top - 0.013, str(k), fontsize=7.5, color=INK, weight="bold", ha="center", va="center",
                  bbox=dict(boxstyle="circle,pad=0.25", fc="white", ec=INK, lw=0.7))
        side.text(0.1, top, f"{value:,.0f} W/m", fontsize=8.4, color=INK, va="top", weight="bold")
        side.text(0.1, top - 0.03, textwrap.fill(name, 46), fontsize=7.6, color=INK2, va="top", linespacing=1.2)
    cax = box(margin + 1.2, 0.42, map_w - 2.4, 0.15)
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=ORANGE), cax=cax, orientation="horizontal", extend="both",
                       ticks=[5, 10, 20, 50, 100, 200, 500])
    bar.ax.set_xticklabels(["5", "10", "20", "50", "100", "200", "500"])
    bar.set_label("Mean wave power travelling toward the coast (W per metre of coast); larger dots carry more",
                  fontsize=8, color=INK2)
    bar.ax.tick_params(labelsize=7, colors=INK2)
    bar.outline.set_visible(False)
    fig.text(margin / width, 1 - 0.3 / height, "Where the nearside sea's waves reach the coast", fontsize=13,
             color=INK, ha="left", va="top")
    stats = coast["mean_power_W_m"]
    note(fig, margin / width, 1 - 0.6 / height,
         "Time-mean shoreward wave energy through the second lunar cycle at the 1-degree sea nodes beside land, each "
         "drawn as a dot: the exposure of about 30 km of coast. Median coast "
         f"{stats['median']:.0f} W/m, 90th percentile {stats['p90']:.0f} W/m, strongest {stats['max']:,.0f} W/m. "
         "Headlands and bays inside a node's stretch need nested grids (the shore study at Smythii). Equal-area map "
         "centred on the sea.", 200, size=7.8)
    path = OUT / "nearside_coasts.png"
    fig.savefig(path, dpi=170, facecolor=SURF)
    plt.close(fig)
    write_sidecar(path, [NEARSIDE, ATLAS_GRID], dict(
        coastal_nodes=int(len(power)), colour_range_W_m=[norm.vmin, norm.vmax],
        data_range_W_m=[float(power.min()), float(power.max())],
        listed=[dict(rank=k, name=n, mean_power_W_m=v) for k, n, v in listing]))
    return path


# ------------------------------------------------------------------ a month at the shores

def label_ends(ax, x, ys, texts, colors, gap):
    """Direct labels at the ends of lines, nudged apart vertically; ink text beside a stub in the series colour."""
    order = np.argsort(ys)
    placed = []
    for k in order:
        y = ys[k] if not placed else max(ys[k], placed[-1] + gap)
        placed.append(y)
        ax.plot([x, x + 0.35], [y, y], color=colors[k], lw=1.8, clip_on=False)
        ax.text(x + 0.5, y, texts[k], fontsize=7.2, color=INK2, va="center", clip_on=False)


def shore_month(base, product, month):
    shores = month["shores"]
    start = product["clock"]["report_interval_hours"][0]
    t = (np.array(month["time_hours"]) - start) / 24
    sun = np.array(month["subsolar_longitude_deg"])
    fig = plt.figure(figsize=(15.0, 10.0), facecolor=SURF)
    right = fig.add_gridspec(4, 1, height_ratios=[1, 1, 1, 0.34], hspace=0.12, left=0.40, right=0.885,
                             top=0.86, bottom=0.07)
    axes = [fig.add_subplot(right[k]) for k in range(4)]
    for ax in axes[1:]:
        ax.sharex(axes[0])
    names = [f"{k}  {s['short']}" for k, s in enumerate(shores, 1)]
    quantities = [("hs_m", "Significant wave\nheight (m)"), ("tm01_s", "Mean period (s)"),
                  ("tide_m", "Tide (m)")]
    for ax, (key, label) in zip(axes, quantities):
        ends = []
        for k, shore in enumerate(shores):
            y = np.array(shore[key], float)
            if key == "tm01_s":                       # the period of a few decimetres of sea says little
                y = np.where(np.array(shore["hs_m"]) >= 0.5, y, np.nan)
            ax.plot(t, y, color=SERIES[k], lw=1.35, label=names[k])
            ends.append(y[-1])
        if key != "tm01_s":
            label_ends(ax, t[-1] + 0.25, np.array(ends), names, SERIES,
                       gap=0.065 * (ax.get_ylim()[1] - ax.get_ylim()[0]))
        ax.set_ylabel(label, fontsize=8.5, color=INK2)
        ax.tick_params(labelsize=7.5, colors=INK2)
        ax.grid(axis="y", color="#e6e5e1", lw=0.6)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        plt.setp(ax.get_xticklabels(), visible=False)
    axes[0].set_ylim(bottom=0)
    axes[1].set_ylim(5, 21)
    axes[1].text(0.005, 0.97, "where Hs is at least 0.5 m", transform=axes[1].transAxes, fontsize=7, color=INK2,
                 va="top")
    second = axes[1].secondary_yaxis(1.0, functions=(wavelength, period_of))
    second.set_ticks([10, 25, 50, 75, 100])
    second.tick_params(labelsize=7, colors=INK2)
    second.set_ylabel("deep-water\nwavelength (m)", fontsize=7.5, color=INK2)
    axes[2].axhline(0, color="#b9b8b2", lw=0.7, zorder=0)
    axes[0].legend(loc="upper left", fontsize=7.4, frameon=False, ncol=4, bbox_to_anchor=(0, 1.2))
    # Lunar day and night at each shore: daylight where the Sun stands within 90 degrees of its longitude.
    night = axes[3]
    for k, shore in enumerate(shores):
        dark = np.abs(((shore["lon_lat"][0] - sun) + 180) % 360 - 180) > 90
        row = len(shores) - 1 - k
        night.fill_between(t, row + 0.12, row + 0.88, where=dark, color="#4a4a52", lw=0, step="mid")
        night.fill_between(t, row + 0.12, row + 0.88, where=~dark, color="#efeee9", lw=0, step="mid")
        night.text(-0.4, row + 0.5, names[k], fontsize=7, color=INK2, ha="right", va="center")
    night.set_ylim(0, len(shores)); night.set_yticks([])
    for side in ("top", "right", "left"):
        night.spines[side].set_visible(False)
    night.set_xlim(0, 31.0)
    night.set_xticks(np.arange(0, 30, 5))
    night.tick_params(labelsize=7.5, colors=INK2)
    night.set_xlabel("Earth days into the second lunar cycle", fontsize=8.5, color=INK2)
    night.text(31.0, len(shores) / 2, "lunar night\nshaded", fontsize=7, color=INK2, va="center")
    # The locator: the shores on the equal-area map.
    proj = Laea(*CENTRE)
    view = base.plane(proj, sea_extent(base, proj), nx=700)
    soft = base.soft(view, [SEA])
    img, _ = compose(view, soft, land=(0.80, 0.17), sea_color=SEA_FLAT, other_water="#d9e0e8")
    loc = fig.add_axes([0.015, 0.50, 0.30, 0.36])
    show(loc, view, img)
    coastline(loc, view, soft, lw=0.35)
    for k, shore in enumerate(shores):
        x, y, _ = proj.forward(shore["lon_lat"][1], shore["lon_lat"][0])
        loc.plot(x, y, "o", ms=11, mfc=SERIES[k], mec="white", mew=1.0, zorder=6)
        loc.text(x, y, str(k + 1), fontsize=7, color="white", ha="center", va="center", weight="bold", zorder=7)
    side = fig.add_axes([0.02, 0.07, 0.31, 0.40]); side.axis("off")
    for k, shore in enumerate(shores):
        y = 1.0 - k * 0.25
        side.plot(0.025, y - 0.02, "o", ms=11, mfc=SERIES[k], mec="white", mew=1.0, transform=side.transAxes)
        side.text(0.025, y - 0.02, str(k + 1), fontsize=7, color="white", ha="center", va="center", weight="bold",
                  transform=side.transAxes)
        side.text(0.08, y, shore["name"], fontsize=8, color=INK, va="top", transform=side.transAxes)
        typical = shore["tide_range_matched_months_m"]["median"]
        detail = (f"{shore['depth_m']:.0f} m deep; {shore['mean_power_W_m']:,.0f} W/m toward the coast on average "
                  f"({shore['reason']}). Tide range this month {shore['tide_range_this_month_m']:.1f} m, typical "
                  f"{typical:.1f} m.")
        side.text(0.08, y - 0.045, textwrap.fill(detail, 62), fontsize=7.1, color=INK2, va="top",
                  transform=side.transAxes, linespacing=1.2)
    tide = month["tide_month"]
    fig.suptitle("A month at four shores of the nearside sea", fontsize=13, color=INK, x=0.015, ha="left", y=0.975)
    first, last = tide["start_tt"][:10], tide["end_tt"][:10]
    note(fig, 0.015, 0.945,
         "Hourly sea states at four coastal 1-degree nodes through the nearside run's second lunar cycle, with the "
         "equilibrium tide at the same points. The GCM's month has no date, so the tide is a real one, "
         f"{first} to {last}: of the {tide['candidates']} months of 2026–2045 matched to the GCM's Sun at mid-month "
         f"(within {tide['sun_mismatch_max_deg']:.1f}° through the month), the one closest to each shore's typical "
         f"range (within {100 * tide['largest_departure_from_typical_range']:.0f}% at every shore). Air pressure moves "
         "the level by under 0.2 m and wind by about 1 cm more.", 215, size=7.6)
    path = OUT / "nearside_shores.png"
    fig.savefig(path, dpi=170, facecolor=SURF)
    plt.close(fig)
    write_sidecar(path, [NEARSIDE, SHORE_MONTH, ATLAS_GRID], dict(
        shores=[dict(name=s["name"], hs_max_m=max(s["hs_m"]), tide_range_m=float(np.ptp(s["tide_m"])),
                     night_share=float(np.mean(np.abs(((s["lon_lat"][0] - sun) + 180) % 360 - 180) > 90)))
                for s in shores]))
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
    bodies = [int(b) for b in product["seas"]]
    # Contours where a sea spans several metres of range; the smaller seas carry their values at stations.
    contoured = [int(b) for b, s in product["seas"].items() if s["area_km2"] > 1e6]
    fig, axes = plt.subplots(1, 3, figsize=(18, 8.0), facecolor=SURF)
    fig.subplots_adjust(left=0.015, right=0.985, top=0.82, bottom=0.1, wspace=0.03)
    norm = mcolors.Normalize(0, 6.0)
    checks = {}
    views = (("Near side, facing Earth (0°)", 0.0), ("Eastern limb (90° E): Smythii–Marginis, Humboldtianum", 90.0),
             ("Far side (180°)", 180.0))
    for ax, (title, lon0) in zip(axes, views):
        proj = Ortho(lon0)
        view = base.plane(proj, (-1.02, 1.02, -1.02, 1.02), nx=1100)
        soft = base.soft(view, bodies)
        values = sample(field[::-1], float(lat[-1]), float(lon[0]), 0.25, view["lat"], view["lon"])
        img, sea = compose(view, soft, values, ORANGE, norm, bodies=bodies)
        show(ax, view, img)
        ax.add_patch(plt.Circle((0, 0), 1.0, fill=False, ec="#8f8e8a", lw=0.8, zorder=2))
        coastline(ax, view, soft, lw=0.45)
        large = base.soft(view, contoured) >= 0.5
        shown = np.where(large & view["inside"], values, np.nan)
        cs = ax.contour(view["X"], view["Y"], shown, levels=np.arange(1, 7), colors=INK, linewidths=0.6, alpha=0.85,
                        zorder=4)
        ax.clabel(cs, fmt=lambda v: f"{v:g} m", fontsize=7, colors=INK)
        for station in product["stations"]:
            headland = station["name"].startswith("Eastern")
            if headland and lon0 != 90.0:
                continue
            sx, sy, sv = proj.forward(station["lon_lat"][1], station["lon_lat"][0])
            if sv and abs(math.cos(math.radians(station["lon_lat"][0] - lon0))) * math.cos(math.radians(station["lon_lat"][1])) > 0.2:
                name = "Smythii headland (wave studies)" if headland else station["name"].replace("Mare ", "").replace("Oceanus ", "")
                ax.plot(sx, sy, "o", ms=3.5, mfc=INK, mec="white", mew=0.8, zorder=6)
                ax.annotate(f"{name} {station['monthly_range_median_m']:.{2 if station['monthly_range_median_m'] < 1 else 1}f} m",
                            (sx, sy), xytext=(5, 4), textcoords="offset points", fontsize=7, color=INK,
                            path_effects=HALO, zorder=6)
        ax.set_title(title, fontsize=10, color=INK, loc="left")
        checks[title] = dict(drawn_max_m=float(np.nanmax(np.where(sea, values, np.nan))))
    # Round trip: the drawn field at each station equals the station's typical monthly range.
    st = product["stations"]
    drawn = sample(field[::-1], float(lat[-1]), float(lon[0]), 0.25, np.array([s["lon_lat"][1] for s in st]),
                   np.array([s["lon_lat"][0] for s in st]))
    checks["station_round_trip_max_abs_m"] = float(np.max(np.abs(drawn - np.array([s["monthly_range_median_m"] for s in st]))))
    colorbar(fig, axes, ORANGE, norm, "Typical monthly tide range (m): highest minus lowest water in an anomalistic month, median over 19 years")
    fig.suptitle("The monthly tide of the Open Moon's seas", fontsize=13, color=INK, x=0.02, ha="left", y=0.975)
    seas = product["seas"]
    text = ("Equilibrium tide from the Earth and the Sun, 2026–2045, with each sea's volume held and the water's "
            "self-attraction. Typical monthly range over each sea's area: nearside "
            f"{seas['433']['monthly_range_median_m']['area_median']:.1f} m, South Pole–Aitken "
            f"{seas['1417']['monthly_range_median_m']['area_median']:.1f} m, Smythii–Marginis "
            f"{seas['874']['monthly_range_median_m']['area_median']:.2f} m. Mare Fecunditatis's strait sets its "
            "fortnightly tide: twice this equilibrium with bed friction for slow flows, 0.8–1.05 times it with the "
            "friction of the metre-per-second flows the strait would carry (geography/README.md). Contours every "
            "metre in the two large seas.")
    text = text.replace(" m,", " m,").replace(" m.", " m.")   # keep values with their units
    note(fig, 0.02, 0.925, text, 175, size=8)
    path = OUT / "tides.png"
    fig.savefig(path, dpi=170, facecolor=SURF)
    plt.close(fig)
    write_sidecar(path, [TIDES, TIDE_GRID, ATLAS_GRID], checks)
    return path


def main():
    which = sys.argv[1:] or ["tides", "waves", "shores"]
    OUT.mkdir(parents=True, exist_ok=True)
    base = Basemap()
    if "tides" in which and TIDE_GRID.exists():
        print(tide_maps(base))
    if NEARSIDE.exists():
        product = json.loads(NEARSIDE.read_text())
        if "waves" in which:
            print(wave_maps(base, product, json.loads(SMYTHII.read_text())))
            print(coast_map(base, product))
        if "shores" in which and SHORE_MONTH.exists():
            print(shore_month(base, product, json.loads(SHORE_MONTH.read_text())))


if __name__ == "__main__":
    main()
