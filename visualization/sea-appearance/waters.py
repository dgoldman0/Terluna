"""The colour of the Open Moon's waters: reflectance spectra and swatches under the lunar and the Earth daylight.

    python visualization/sea-appearance/waters.py

Reads research/studies/sea_appearance/results/water_colours.json. Swatches show the light leaving each water at
one exposure, on a display with Earth-daylight (D65) white and no adaptation to the lunar daylight, so the
lunar sky's warmer light shows. Writes an ignored figure with a provenance sidecar to results/.
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
sys.path.insert(0, str(Path(__file__).resolve().parent))
from lighting import INK, INK2, OUT, SURF, write_sidecar  # noqa: E402

COLOURS = ROOT / "research/studies/sea_appearance/results/water_colours.json"
SERIES = {"open_sea": "#2a78d6", "productive_coast": "#1baf7a", "river_mouth": "#eb6834"}   # slots 1, 3, 2
SOIL_DASH = {"12001": "-", "10084": (0, (5, 2.5)), "62231": (0, (1.5, 1.8))}
LABELS = {("open_sea", "12001"): "Open sea", ("productive_coast", "12001"): "Productive coast, mare fines",
          ("productive_coast", "10084"): "Productive coast, high-titanium mare fines",
          ("productive_coast", "62231"): "Productive coast, highland fines", ("river_mouth", "12001"): "River mouth"}
XYZ_TO_SRGB = np.array([[3.2406, -1.5372, -0.4986], [-0.9689, 1.8758, 0.0415], [0.0557, -0.2040, 1.0570]])


def srgb(xyz, white):
    linear = np.clip(XYZ_TO_SRGB @ (np.asarray(xyz) / white), 0, 1)
    return np.where(linear <= 0.0031308, 12.92 * linear, 1.055 * linear ** (1 / 2.4) - 0.055)


def figure(product):
    cases = [c for c in product["cases"] if c["phase"] == "mean"]
    fig = plt.figure(figsize=(15.5, 9.2), facecolor=SURF)
    ax = fig.add_axes([0.06, 0.09, 0.40, 0.73])
    wavelength = np.array(product["display_wavelength_nm"])
    for c in cases:
        ax.plot(wavelength, np.array(c["rrs_display_per_sr"]) * 1000, color=SERIES[c["water"]],
                ls=SOIL_DASH[c["soil"]], lw=1.6)
    ax.set_xlim(400, 700); ax.set_ylim(0, None)
    ax.set_xlabel("Wavelength (nm)", fontsize=9, color=INK2)
    ax.set_ylabel("Remote-sensing reflectance (10⁻³ per steradian)", fontsize=9, color=INK2)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(labelsize=8, colors=INK2)
    handles = [plt.Line2D([], [], color=SERIES[c["water"]], ls=SOIL_DASH[c["soil"]], lw=1.6) for c in cases]
    ax.legend(handles, [LABELS[(c["water"], c["soil"])] for c in cases], loc="upper right", frameon=False, fontsize=8)
    # Swatches: one exposure for all, so brightness compares across waters and worlds.
    white = 1.15 * max(c[world]["xyz_cd_m2"][1] for c in cases for world in ("moon", "earth"))
    height, aspect = 0.098, 9.2 / 15.5
    columns = dict(moon=0.705, earth=0.855)
    for world, x in columns.items():
        title = "Moon's daylight" if world == "moon" else "Earth's daylight"
        fig.text(x, 0.80, f"{title}\n{product['daylight_illuminance_lux'][world]:,.0f} lux, Sun 45° up",
                 fontsize=8.5, color=INK, va="bottom")
    checks = {}
    for k, c in enumerate(cases):
        y = 0.68 - k * 0.14
        fig.text(0.53, y + height / 2, textwrap.fill(LABELS[(c["water"], c["soil"])], 24), fontsize=8.2,
                 color=INK, va="center")
        for world, x in columns.items():
            box = fig.add_axes([x, y, height * aspect, height])
            colour = srgb(c[world]["xyz_cd_m2"], white)
            box.set_facecolor(colour); box.set_xticks([]); box.set_yticks([])
            for spine in box.spines.values():
                spine.set_visible(False)
            fig.text(x + height * aspect + 0.008, y + height * 0.62, f"{c[world]['luminance_cd_m2']:.0f} cd/m²",
                     fontsize=7.6, color=INK2, va="center")
            fig.text(x + height * aspect + 0.008, y + height * 0.32,
                     f"1% light at {c[world]['depth_one_percent_photons_m']:.0f} m", fontsize=7.6, color=INK2,
                     va="center")
            checks[f"{world} {c['water']} {c['soil']}"] = dict(srgb=np.round(colour, 3).tolist(),
                                                              xy=c[world]["chromaticity_xy"])
    fig.text(0.06, 0.965, "The colour of the Open Moon's waters", fontsize=13, color=INK, ha="left", va="top")
    note = ("The light leaving the water itself, seen looking down past the reflected sky, for the biosphere's design "
            "guesses of what the seas carry: a nutrient-poor open sea (0.1 mg/m³ chlorophyll), a productive coast "
            "(2 mg/m³, organic matter, 3 g/m³ of regolith fines) and a river mouth (5 mg/m³, much organic matter, "
            "30 g/m³ of fines). Earth's seawater optics; fines from the finest fraction of Apollo soils. Swatches share "
            "one exposure on a display with Earth-daylight (D65) white, without adapting to the lunar daylight, which "
            "reaches the water warmer than Earth's: the Moon's tall air takes more blue out of the sunbeam. Depths are "
            "where photosynthetic light falls to 1%. Saturated blues beyond the display's range are clipped.")
    fig.text(0.06, 0.93, textwrap.fill(note, 205), fontsize=7.8, color=INK2, ha="left", va="top")
    path = OUT / "water_colours.png"
    fig.savefig(path, dpi=160, facecolor=SURF)
    plt.close(fig)
    write_sidecar(path, [COLOURS], dict(display_white_cd_m2=round(white, 2), swatches=checks))
    return path


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    print(figure(json.loads(COLOURS.read_text())))


if __name__ == "__main__":
    main()
