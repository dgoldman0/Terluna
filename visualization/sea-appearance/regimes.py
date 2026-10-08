"""How the seas look by regime: panoramas of sky and sea at four coasts through the month.

    python visualization/sea-appearance/regimes.py

Reads research/studies/sea_appearance/results/regimes.json and its panoramas on the research drive
(research/runs/sea_appearance/regimes.npz). Each strip is one moment: the whole horizon, 30 degrees above and
below, from an eye 2 m above the water, with the sea in every direction. Colours are as calculated, on a display
with Earth-daylight (D65) white, without adapting to the light of the moment; each strip has its own exposure,
printed as the luminance shown white, and brighter highlights roll off toward white. Writes ignored figures with
provenance sidecars to results/.
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
from waters import XYZ_TO_SRGB  # noqa: E402
sys.path.insert(0, str(ROOT))
from matplotlib.patches import Ellipse  # noqa: E402
from research.studies.sea_appearance import lighting as study_lighting  # noqa: E402
from shared.constants import AU, EARTH_MOON_DISTANCE, EARTH_RADIUS, SUN_RADIUS  # noqa: E402
SHOULDER = 0.8

REGIMES = ROOT / "research/studies/sea_appearance/results/regimes.json"
CALENDAR = ROOT / "research/studies/sea_appearance/results/lighting_calendar.json"
HALO = [matplotlib.patheffects.withStroke(linewidth=2.2, foreground="black")]
MOMENT_LABELS = {"noon": "The Sun at its highest", "low Sun": "The Sun 4° up",
                 "after sunset": "The Sun 4° below the horizon", "evening": "The Sun 15° below",
                 "deep evening": "The Sun 35° below", "darkest": "The darkest hour"}


def exposure(xyz, elevation):
    """Luminance shown white: 1.25 times the 95th percentile of the sky's, as a photographer exposes a bright sky."""
    return 1.25 * float(np.percentile(xyz[elevation > 0][..., 1], 95))


def encode(linear):
    linear = np.clip(linear, 0, 1)
    return np.where(linear <= 0.0031308, 12.92 * linear, 1.055 * linear ** (1 / 2.4) - 0.055)


def strip(xyz, white):
    """sRGB for a strip of XYZ radiance at one exposure. Luminance above 80% of white rolls off smoothly toward
    white, keeping its chromaticity, so the glitter keeps its shape instead of clipping flat."""
    y = xyz[..., 1] / white
    rolled = np.where(y <= SHOULDER, y, SHOULDER + (1 - SHOULDER) * (1 - np.exp(-(y - SHOULDER) / (1 - SHOULDER))))
    scale = np.where(y > 0, rolled / np.maximum(y, 1e-30), 0.0)
    return encode((xyz / white * scale[..., None]) @ XYZ_TO_SRGB.T)


class Disks:
    """The colours of the Sun's and the Earth's disks seen through the air, from the solved sky's direct beam."""

    def __init__(self):
        self.sky = study_lighting.Sky()
        self.earth = study_lighting.Earthlight(self.sky)
        self.xyz = np.load(study_lighting.SPHERICAL / "moon_1.2atm_standard.npz")["xyz"]

    def colour(self, body, elevation, phase=0.0):
        beam = np.array([np.interp(elevation, self.sky.beam_angles, self.sky.direct_channels[:, c])
                         for c in range(self.sky.direct_channels.shape[1])])
        if body == "earth":
            beam = beam * self.earth.weights(np.array([phase]), np.array([EARTH_MOON_DISTANCE]), np.ones(1))[0]
        rgb = np.clip((beam @ self.xyz) @ XYZ_TO_SRGB.T, 0, None)
        return encode(rgb / max(rgb.max(), 1e-30))


def figure(product, panoramas, coast, month_start, disks):
    scenes = [(k, s) for k, s in enumerate(product["scenes"]) if s["coast"] == coast]
    elevation, azimuth = panoramas["elevation_deg"], panoramas["azimuth_deg"]
    fig = plt.figure(figsize=(16, 1.6 + 2.62 * len(scenes)), facecolor=SURF)
    top = 1 - 1.2 / fig.get_figheight()
    height = 2.25 / fig.get_figheight()
    checks = {}
    for row, (k, s) in enumerate(scenes):
        xyz = panoramas["xyz"][k].astype(float)
        white = exposure(xyz, elevation)
        ax = fig.add_axes([0.05, top - (row + 1) * (height + 0.37 / fig.get_figheight()), 0.93, height])
        ax.imshow(strip(xyz, white), extent=(0, 360, elevation[-1] - 0.25, elevation[0] + 0.25), aspect="auto",
                  interpolation="bilinear")
        ax.axhline(0, color="white", lw=0.4, alpha=0.5)
        radii = dict(sun=np.degrees(SUN_RADIUS / AU), earth=np.degrees(np.arcsin(EARTH_RADIUS / EARTH_MOON_DISTANCE)))
        for body in ("sun", "earth"):
            e, a = s[body]["elevation_deg"], s[body]["azimuth_deg"]
            if 0 < e < 30:
                colour = disks.colour(body, e, s["earth"]["phase_angle_deg"])
                r = max(radii[body], 0.45)          # at least a pixel or two across
                ax.add_patch(Ellipse((a, e), 2 * r, 2 * r, facecolor=colour, edgecolor="none", zorder=3))
                ax.text(a + 2.5, e + 1.5, body.title(), color="white", fontsize=7.5, path_effects=HALO)
        ax.set_xticks(range(0, 361, 45))
        ax.set_xticklabels(["N", "NE", "E", "SE", "S", "SW", "W", "NW", "N"])
        ax.set_yticks([-30, -15, 0, 15, 30])
        ax.tick_params(labelsize=7, colors=INK2)
        for side in ax.spines.values():
            side.set_visible(False)
        label = (f"{MOMENT_LABELS[s['moment']]}  ·  day {(s['hour'] - month_start) / 24:.1f}"
                 f"  ·  Sun {s['sun']['elevation_deg']:+.1f}°, Earth {s['earth']['elevation_deg']:+.1f}°"
                 f" ({100 * s['earth']['lit_fraction']:.0f}% lit)  ·  {s['light_on_ground_lux']:,.3g} lux on the ground"
                 f"  ·  mean square slope {s['slopes']['mean_square']:.3f}  ·  white = {white:,.3g} cd/m²")
        ax.set_title(label, fontsize=8.5, color=INK, loc="left", pad=3)
        checks[s["moment"]] = dict(white_cd_m2=round(white, 6), hour=s["hour"])
    fig.text(0.05, 1 - 0.35 / fig.get_figheight(), f"How the sea looks through the month: {scenes[0][1]['coast']}",
             fontsize=13, color=INK, va="top")
    note = ("The whole horizon from an eye 2 m above the water, 30° above and below, with the sea in every direction "
            "(the coastline is left out). Clear sky from the solved atmosphere lit by the Sun and the Earth; the sea "
            "reflects it through that hour's wave slopes, with the glitter of both disks and the water's own light "
            "for a productive coast. Each strip has its own exposure, printed as the luminance shown white; brighter "
            "highlights roll off from 80% of white. The Sun's and the Earth's disks are drawn at their colour through "
            "the air, far brighter than white. Colours as calculated on a display with Earth-daylight (D65) white, "
            "unadapted; below "
            "a few cd/m² the eye sees colour ever more weakly, so the darkest strips show what a long photographic "
            "exposure records. Days count from 7 February 2038.")
    fig.text(0.05, 1 - 0.68 / fig.get_figheight(), textwrap.fill(note, 230), fontsize=7.8, color=INK2, va="top")
    path = OUT / f"regimes_{scenes[0][1]['coast'].replace(' ', '_').lower()}.png"
    fig.savefig(path, dpi=150, facecolor=SURF)
    plt.close(fig)
    write_sidecar(path, [REGIMES], checks)
    return path


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    product = json.loads(REGIMES.read_text())
    month_start = json.loads(CALENDAR.read_text())["time_hours"][0]
    with np.load(ROOT / product["panoramas"]["path"]) as z:
        panoramas = {k: z[k] for k in z.files}
    disks = Disks()
    for coast in dict.fromkeys(s["coast"] for s in product["scenes"]):
        print(figure(product, panoramas, coast, month_start, disks))


if __name__ == "__main__":
    main()
