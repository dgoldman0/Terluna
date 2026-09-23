#!/usr/bin/env python3
"""Bake the sky bodies: the site ephemeris, the bright stars and Earth's disk calibration.

Reads the illumination domain's ephemeris model and star catalogue and the site in
shared/scenarios/sites.json, and writes assets/sky/site-sky.json and assets/sky/stars.json.
The Earth imagery is checked against assets/earth/manifest.json; a gain is written to
assets/earth/calibration.json so the drawn disk's full-phase brightness matches the
geometric albedo used for earthlight.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
IMMERSION = HERE.parent
REPO = IMMERSION.parent
sys.path.insert(0, str(REPO))

import numpy as np  # noqa: E402
from PIL import Image  # noqa: E402

from illumination.ephemeris import Site, product  # noqa: E402
from shared import constants as K  # noqa: E402

SITE_ID = "development"
CLOUD_ALBEDO = 0.75  # linear albedo of opaque cloud in the disk shader
HAZE = [0.02, 0.035, 0.08]  # Earth's own Rayleigh sky, added to the disk albedo


def earth_calibration() -> dict:
    folder = IMMERSION / "assets/earth"
    manifest = json.loads((folder / "manifest.json").read_text())
    for name, spec in manifest["files"].items():
        digest = hashlib.sha256((folder / name).read_bytes()).hexdigest()
        if digest != spec["sha256"]:
            raise SystemExit(f"Earth imagery hash mismatch: {name}")
    rgb = np.asarray(Image.open(folder / "land_shallow_topo_2048.jpg").convert("RGB"), float) / 255
    cloud = np.asarray(Image.open(folder / "cloud_combined_2048.jpg").convert("L"), float)[..., None] / 255
    surface = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    albedo = surface * (1 - cloud) + CLOUD_ALBEDO * cloud
    weight = np.cos(np.linspace(np.pi / 2, -np.pi / 2, albedo.shape[0]))[:, None, None]
    mean = (albedo * weight).sum(axis=(0, 1)) / weight.sum() / albedo.shape[1]
    luminance = float(mean @ [0.2126, 0.7152, 0.0722]) + float(np.dot(HAZE, [0.2126, 0.7152, 0.0722]))
    # A Lambert sphere's geometric albedo is 2/3 of its (uniform) albedo.
    gain = 1.5 * K.EARTH_GEOMETRIC_ALBEDO / luminance
    return {"schema": "terluna.immersion.earth-calibration/1", "gain": gain, "cloud_albedo": CLOUD_ALBEDO,
            "haze": HAZE, "mean_albedo_before_gain": mean.tolist(),
            "target_geometric_albedo": K.EARTH_GEOMETRIC_ALBEDO,
            "note": "Global-mean calibration so the disk's full-phase brightness matches earthlight; informed, not measured."}


def main() -> int:
    spec = json.loads((REPO / "shared/scenarios/sites.json").read_text())["sites"][SITE_ID]
    sky = IMMERSION / "assets/sky"
    (sky / "site-sky.json").write_text(json.dumps(product(SITE_ID, Site.from_scenario(spec))))
    stars = json.loads((REPO / "illumination/stars/bright_stars.json").read_text())
    (sky / "stars.json").write_text(json.dumps(stars, separators=(",", ":")))
    calibration = earth_calibration()
    (IMMERSION / "assets/earth/calibration.json").write_text(json.dumps(calibration, indent=2) + "\n")
    print(f"site-sky ({SITE_ID}), {len(stars['stars'])} stars, Earth gain {calibration['gain']:.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
