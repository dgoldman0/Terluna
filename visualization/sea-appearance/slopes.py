"""The sea's slopes: the short-wave spectrum at lunar and Earth gravity, and a month of slopes at six coasts.

    python visualization/sea-appearance/slopes.py

Reads research/studies/sea_appearance/results/sea_slopes.json and sea_slopes.npz, and evaluates
illumination/water_surface/short_waves.py for the spectra. Writes ignored figures with provenance sidecars to
results/.
"""
from __future__ import annotations

import json
import sys
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from lighting import INK, INK2, OUT, SURF, write_sidecar  # noqa: E402
from illumination.water_surface import short_waves as S  # noqa: E402
from shared.constants import (EXOPLASIM_CHARNOCK, EXOPLASIM_SEA_ROUGHNESS_FLOOR, EXOPLASIM_VON_KARMAN,  # noqa: E402
                              MOON_SURFACE_GRAVITY, STANDARD_GRAVITY)

SLOPES = ROOT / "research/studies/sea_appearance/results/sea_slopes.json"
SERIES = ROOT / "research/studies/sea_appearance/results/sea_slopes.npz"
EARTH, MOON = "#2a78d6", "#eb6834"          # categorical slots 1 and 2
GLASSY = "#ecebe6"
HALO = [matplotlib.patheffects.withStroke(linewidth=2.4, foreground="white")]


def friction_velocity(u10, gravity):
    """Invert the GCM's neutral surface layer (Charnock's roughness) for the friction velocity under a 10 m wind."""
    low, high = 1e-4, 2.0
    for _ in range(80):
        middle = (low + high) / 2
        wind = S.neutral_wind(middle, gravity, EXOPLASIM_CHARNOCK, EXOPLASIM_SEA_ROUGHNESS_FLOOR, EXOPLASIM_VON_KARMAN)
        low, high = (middle, high) if wind < u10 else (low, middle)
    return (low + high) / 2


def spectra_figure(product):
    fig, ax = plt.subplots(figsize=(11.5, 6.6), facecolor=SURF)
    fig.subplots_adjust(left=0.08, right=0.97, top=0.76, bottom=0.12)
    wavelength = np.geomspace(400, 0.002, 3000)
    k = 2 * np.pi / wavelength
    cutoff = product["coasts"][0]["cutoff_wavelength_m"]
    ax.axvspan(400, cutoff, color="#f2f1ec", lw=0, zorder=0)
    ax.text(np.sqrt(400 * cutoff), 2.2e-2, "the wave runs resolve these waves", ha="center", fontsize=8, color=INK2)
    ax.text(cutoff * 0.85, 2.2e-2, "the unified spectrum adds these", ha="left", fontsize=8, color=INK2)
    checks = {}
    for gravity, colour, world in ((MOON_SURFACE_GRAVITY, MOON, "Moon"), (STANDARD_GRAVITY, EARTH, "Earth")):
        k_m, _ = S.capillary_scales(gravity)
        for u10, style in ((3.0, (0, (5, 2.5))), (5.0, "-")):
            u_star = friction_velocity(u10, gravity)
            b = sum(S.curvature(k, u10, u_star, S.FULLY_DEVELOPED, gravity))
            ax.plot(wavelength, b, color=colour, lw=1.7 if style == "-" else 1.4, ls=style)
            checks[f"{world} {u10:g} m/s"] = dict(u_star_m_s=round(u_star, 4), peak_curvature=float(b.max()))
        ax.axvline(2 * np.pi / k_m, color=colour, lw=0.8, ls=(0, (1, 2)))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(400, 0.002); ax.set_ylim(1e-4, 0.03)
    ax.set_xticks([100, 10, 1, 0.1, 0.01, 0.001])
    ax.set_xticklabels(["100 m", "10 m", "1 m", "10 cm", "1 cm", "1 mm"])
    ax.set_xlabel("Wavelength", fontsize=9, color=INK2)
    ax.set_ylabel("Curvature spectrum B(k) = k³S(k)", fontsize=9, color=INK2)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(labelsize=8, colors=INK2)
    moon_m = 2 * np.pi / S.capillary_scales(MOON_SURFACE_GRAVITY)[0]
    earth_m = 2 * np.pi / S.capillary_scales(STANDARD_GRAVITY)[0]
    ax.text(moon_m * 1.08, 1.6e-4, f"lunar capillary\nscale {100 * moon_m:.1f} cm", ha="right", fontsize=7.5, color=MOON,
            path_effects=HALO)
    ax.text(earth_m * 0.93, 1.6e-4, f"Earth's\n{100 * earth_m:.1f} cm", ha="left", fontsize=7.5, color=EARTH,
            path_effects=HALO)
    handles = [plt.Line2D([], [], color=MOON, lw=1.7), plt.Line2D([], [], color=EARTH, lw=1.7),
               plt.Line2D([], [], color=INK2, lw=1.7), plt.Line2D([], [], color=INK2, lw=1.4, ls=(0, (5, 2.5)))]
    fig.legend(handles, ["Lunar gravity", "Earth's gravity", "10 m wind 5 m/s", "10 m wind 3 m/s"],
               loc="upper left", bbox_to_anchor=(0.08, 0.83), ncol=4, frameon=False, fontsize=8.5)
    fig.text(0.08, 0.965, "Short waves on the Moon and on Earth under the same wind", fontsize=13, color=INK,
             ha="left", va="top")
    note = ("The unified spectrum of Elfouhaily et al. (1997) for a fully developed sea, with gravity and surface "
            "tension explicit. At one sixth of Earth's gravity the gravity-capillary waves are 2.5 times as long, the "
            f"slowest ripples travel at {product['minimum_phase_speed_m_s']:.3f} m/s against "
            f"{product['earth_minimum_phase_speed_m_s']:.3f}, and the same wind drives them harder, so the short "
            "waves stand steeper. Each world's friction velocity comes from the same neutral surface layer with "
            "Charnock's roughness at its own gravity. 3 and 5 m/s span the coasts' median and 90th-percentile winds.")
    fig.text(0.08, 0.925, textwrap.fill(note, 175), fontsize=8, color=INK2, ha="left", va="top")
    path = OUT / "slope_spectra.png"
    fig.savefig(path, dpi=160, facecolor=SURF)
    plt.close(fig)
    write_sidecar(path, [SLOPES], checks)
    return path


def runs(mask):
    """Start and end indices of the True runs in a boolean series."""
    edges = np.flatnonzero(np.diff(np.r_[0, mask.astype(int), 0]))
    return edges.reshape(-1, 2)


def month_figure(product):
    with np.load(SERIES) as z:
        series = {k: z[k] for k in z.files}
    hours = series["hours"]
    days = (hours - hours[0]) / 24
    fig = plt.figure(figsize=(15.5, 11.0), facecolor=SURF)
    grid = fig.add_gridspec(3, 2, left=0.06, right=0.985, top=0.835, bottom=0.06, hspace=0.38, wspace=0.12)
    checks = {}
    for n, coast in enumerate(product["coasts"]):
        ax = fig.add_subplot(grid[n // 2, n % 2])
        key = coast["short"].replace(" ", "_").lower()
        short = series[f"{key}/short_covariance"].reshape(-1, 2, 2)
        short_mss = np.trace(short, axis1=1, axis2=2)
        glassy = series[f"{key}/short_level"] == 0
        for a, b in runs(glassy):
            ax.axvspan(days[a] - 1 / 48, days[b - 1] + 1 / 48, color=GLASSY, lw=0, zorder=0)
        ax.plot(days, series[f"{key}/earth_cox_munk_mss"], color=EARTH, lw=1.1, alpha=0.9, zorder=2)
        ax.plot(days, short_mss, color=MOON, lw=1.1, zorder=3)
        resolved_hours = series[f"{key}/resolved_hours"]
        if len(resolved_hours):
            at = np.searchsorted(hours, resolved_hours)
            resolved = series[f"{key}/resolved_covariance"].reshape(-1, 2, 2)
            total = np.trace(resolved + short[at], axis1=1, axis2=2)
            x = days[at]
            if len(at) > 100:
                ax.plot(x, total, color=INK, lw=1.2, zorder=4)
            else:
                ax.plot(x, total, "o", ms=5, mfc=INK, mec="white", mew=0.8, zorder=4)
        else:
            ax.text(0.98, 0.93, "No wave run: the resolved waves are unknown here", transform=ax.transAxes,
                    ha="right", va="top", fontsize=7.5, color=INK2)
        ax.set_xlim(0, 29.6); ax.set_ylim(0, 0.065)
        ax.set_title(coast["name"], fontsize=9.5, color=INK, loc="left")
        ax.set_ylabel("Mean square slope", fontsize=8, color=INK2)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.tick_params(labelsize=7, colors=INK2)
        if n >= 4:
            ax.set_xlabel("Days from 7 February 2038", fontsize=8, color=INK2)
        share = glassy.mean()
        ax.text(0.01, 0.93, f"no short waves {100 * share:.0f}% of hours", transform=ax.transAxes, fontsize=7.5,
                color=INK2, va="top")
        checks[coast["short"]] = dict(glassy_share=round(float(share), 4), total_median=coast["total_mss"] and
                                      coast["total_mss"]["p50"])
    handles = [plt.Line2D([], [], color=INK, marker="o", ms=5, mfc=INK, mec="white", lw=1.2),
               plt.Line2D([], [], color=MOON, lw=1.1), plt.Line2D([], [], color=EARTH, lw=1.1),
               plt.Rectangle((0, 0), 1, 1, color=GLASSY)]
    fig.legend(handles, ["All waves (resolved spectrum plus short waves)", "Short waves alone, shorter than 1 m",
                         "Earth's clean sea under the same wind (Cox and Munk 1954)",
                         "Hours too calm for short waves"],
               loc="upper left", bbox_to_anchor=(0.06, 0.885), ncol=4, frameon=False, fontsize=8)
    fig.text(0.06, 0.975, "How steep the sea's surface is through a month at six coasts", fontsize=13, color=INK,
             ha="left", va="top")
    note = ("Mean square slope of the sea surface, the quantity that sets how widely it scatters reflected light. "
            "Short waves come from the unified spectrum at lunar gravity under the GCM's surface stress, hour by hour; "
            "they vanish when the friction velocity falls below "
            f"{product['short_wave_threshold_u_star_m_s']:.3f} m/s. The longer waves come from the wave runs' spectra: "
            "every 48 hours at the four nearside shores (dots), hourly at the Smythii headland (line). The GCM's "
            "stress is a 170 km, three-hourly mean, without gusts, sea breezes or surface films.")
    fig.text(0.06, 0.945, textwrap.fill(note, 215), fontsize=7.8, color=INK2, ha="left", va="top")
    path = OUT / "sea_slopes.png"
    fig.savefig(path, dpi=160, facecolor=SURF)
    plt.close(fig)
    write_sidecar(path, [SLOPES, SERIES], checks)
    return path


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    product = json.loads(SLOPES.read_text())
    print(spectra_figure(product))
    print(month_figure(product))


if __name__ == "__main__":
    main()
