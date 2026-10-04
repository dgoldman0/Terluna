"""Terrain round the sea-appearance study's rendering coasts, from pinned 118 m LOLA rows.

    python -m geography.coast_terrain --download     # restore the pinned row slices (about 350 MB)
    python -m geography.coast_terrain                # export each coast's local grid

Each coast takes a band of rows, 2.5 degrees either side of it, from a LOLA 256 pixel-per-degree radius tile of
the PDS Gridded Data Record, fetched by HTTP byte range as geography/shore_terrain.py fetches the Smythii
headland's rows (which this module reuses for that coast). Radii become heights above the degree-200 GRAIL
geoid (geography/topography.py) and then above the atlas's common sea level. The export resamples them onto a
local grid of metres east and north of the coast's point, by the azimuthal equidistant projection, for
renderers and other local uses; the tide is the caller's to add.

The PDS labels place pixel centres at latitude (LINE_PROJECTION_OFFSET - line) / 256 and east longitude
CENTER_LONGITUDE + (sample - SAMPLE_PROJECTION_OFFSET) / 256, lines and samples counted from zero.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
from scipy.ndimage import map_coordinates

from geography import shore_terrain
from geography.topography import Grid, geoid
from shared.constants import MOON_RADIUS

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = Path(__file__).with_name("coast_inputs.json")
OUTPUT = ROOT / "research/runs/sea_appearance/terrain"
BASE = "https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/lrolol_1xxx/data/lola_gdr/cylindrical/img/"
ATLAS = ROOT / "geography/products/atlas_28pct_4ppd.npz"
HALF_BAND_DEG = 2.5


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def manifest():
    return json.loads(MANIFEST.read_text())


def _fetch(url, path, size, headers=None, expect_range=None):
    partial = path.with_suffix(path.suffix + ".partial")
    if partial.exists():
        raise FileExistsError(f"Inspect the retained partial download: {partial}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with urlopen(Request(url, headers=headers or {}), timeout=90) as response:
        if expect_range and (response.status != 206 or response.headers.get("Content-Range") != expect_range):
            raise ValueError("The terrain server must return the exact requested byte range")
        with partial.open("xb") as stream:
            remaining = size
            while remaining:
                chunk = response.read(min(1024 ** 2, remaining))
                if not chunk:
                    raise ValueError("Incomplete terrain download")
                stream.write(chunk)
                remaining -= len(chunk)
            if response.read(1):
                raise ValueError("Terrain response exceeds the requested byte count")
    return partial


def restore(download=False):
    """Accept each recorded byte slice only with its exact label and hash."""
    for entry in manifest()["slices"]:
        for kind in ("label", "image"):
            path = ROOT / entry["local_" + kind]
            expected = entry["label_sha256" if kind == "label" else "selected_sha256"]
            size = entry["label_bytes" if kind == "label" else "selected_bytes"]
            if path.exists():
                if path.stat().st_size != size or digest(path) != expected:
                    raise ValueError(f"Pinned coastal input differs: {path}")
                continue
            if not download:
                raise FileNotFoundError(f"Restore the coastal input with --download: {path}")
            if kind == "label":
                partial = _fetch(entry["label_url"], path, size)
            else:
                first, last = entry["selected_byte_range"]
                partial = _fetch(entry["source_url"], path, size,
                                 {"Range": f"bytes={first}-{last}", "If-Range": entry["source_etag"]},
                                 f"bytes {first}-{last}/{entry['source_size_bytes']}")
            if digest(partial) != expected:
                raise ValueError(f"Terrain checksum differs; retained for inspection: {partial}")
            partial.rename(path)


def pin(name, tile, latitude_deg, directory):
    """Retrieve a new slice once and record its identity (used to create the manifest's entries)."""
    label_url, source_url = BASE + f"ldem_256_{tile}.lbl", BASE + f"ldem_256_{tile}.img"
    with urlopen(label_url, timeout=90) as response:
        label = response.read()
    text = label.decode()

    def value(key):
        line = next(x for x in text.splitlines() if x.strip().startswith(key + " ") or x.strip().startswith(key + "="))
        return line.split("=", 1)[1].split("<")[0].strip().strip('"')

    lines, samples = int(value("LINES")), int(value("LINE_SAMPLES"))
    row_bytes, ppd = int(value("RECORD_BYTES")), int(value("MAP_RESOLUTION"))
    lpo, spo = float(value("LINE_PROJECTION_OFFSET")), float(value("SAMPLE_PROJECTION_OFFSET"))
    if value("SAMPLE_TYPE") != "LSB_INTEGER" or int(value("SAMPLE_BITS")) != 16 or row_bytes != 2 * samples:
        raise ValueError("Unexpected LOLA tile layout")
    first_row = int(math.floor(lpo - (latitude_deg + HALF_BAND_DEG) * ppd))
    last_row = int(math.ceil(lpo - (latitude_deg - HALF_BAND_DEG) * ppd))
    first_row, last_row = max(first_row, 0), min(last_row, lines - 1)
    head = urlopen(Request(source_url, method="HEAD"), timeout=90).headers
    size_total = int(head["Content-Length"])
    if size_total != lines * row_bytes:
        raise ValueError("The tile's size differs from its label")
    first, last = first_row * row_bytes, (last_row + 1) * row_bytes - 1
    local_image = Path(directory) / f"ldem_256_{name}_rows_{first_row}_{last_row}.img"
    local_label = Path(directory) / f"ldem_256_{tile}.lbl"
    (ROOT / local_label).parent.mkdir(parents=True, exist_ok=True)
    (ROOT / local_label).write_bytes(label)
    partial = _fetch(source_url, ROOT / local_image, last - first + 1,
                     {"Range": f"bytes={first}-{last}", "If-Range": head["ETag"]},
                     f"bytes {first}-{last}/{size_total}")
    sha = digest(partial)
    partial.rename(ROOT / local_image)
    return dict(name=name, source_url=source_url, label_url=label_url, label_sha256=hashlib.sha256(label).hexdigest(),
                label_bytes=len(label), source_size_bytes=size_total, source_etag=head["ETag"],
                source_last_modified=head["Last-Modified"], selected_byte_range=[first, last], selected_sha256=sha,
                selected_bytes=last - first + 1, first_row=first_row, row_count=last_row - first_row + 1,
                row_bytes=row_bytes, lines=lines, samples=samples, pixels_per_degree=ppd,
                line_projection_offset=lpo, sample_projection_offset=spo,
                center_longitude_deg=float(value("CENTER_LONGITUDE")),
                latitude_bounds_deg=[float(value("MINIMUM_LATITUDE")), float(value("MAXIMUM_LATITUDE"))],
                longitude_bounds_deg=[float(value("WESTERNMOST_LONGITUDE")), float(value("EASTERNMOST_LONGITUDE"))],
                sample_type="<i2", scale_m=float(value("SCALING_FACTOR")), offset_m=float(value("OFFSET")),
                product_version=value("PRODUCT_VERSION_ID"), product_creation_date=value("PRODUCT_CREATION_TIME"),
                local_image=str(local_image), local_label=str(local_label))


def smythii_entry():
    """The shore studies' Smythii rows, described in this manifest's terms."""
    s = shore_terrain.restore()
    return dict(name="smythii_headland", first_row=s["first_row"], row_count=s["row_count"],
                samples=s["samples"], pixels_per_degree=s["pixels_per_degree"], line_projection_offset=23039.5,
                sample_projection_offset=46079.5, center_longitude_deg=180.0, scale_m=s["scale_m"],
                local_image=s["local_image"], selected_sha256=s["selected_sha256"],
                manifest=str(shore_terrain.MANIFEST.relative_to(ROOT)))


def slice_for(name):
    if name == "smythii_headland":
        return smythii_entry()
    return next(e for e in manifest()["slices"] if e["name"] == name)


def radius_heights(entry, latitude, longitude):
    """LOLA heights above 1737.4 km at (latitude, longitude) points, bilinear between pixel centres."""
    raw = np.memmap(ROOT / entry["local_image"], dtype="<i2", mode="r",
                    shape=(entry["row_count"], entry["samples"]))
    ppd = entry["pixels_per_degree"]
    row = entry["line_projection_offset"] - np.asarray(latitude) * ppd - entry["first_row"]
    col = entry["sample_projection_offset"] + (np.asarray(longitude) % 360 - entry["center_longitude_deg"]) * ppd
    if row.min() < 0 or row.max() > entry["row_count"] - 1 or col.min() < 0 or col.max() > entry["samples"] - 1:
        raise ValueError(f"The local grid leaves the restored rows of {entry['name']}")
    r0, r1 = int(np.floor(row.min())), int(np.ceil(row.max())) + 1
    c0, c1 = int(np.floor(col.min())), int(np.ceil(col.max())) + 1
    block = np.asarray(raw[r0:r1, c0:c1], dtype=float) * entry["scale_m"]
    return map_coordinates(block, [row - r0, col - c0], order=1, mode="nearest")


def destination(lon0, lat0, x, y):
    """Latitude and east longitude (degrees) of points x east and y north (m) of a point, azimuthal equidistant."""
    distance = np.hypot(x, y) / MOON_RADIUS
    azimuth = np.arctan2(x, y)
    phi1, lam1 = math.radians(lat0), math.radians(lon0)
    phi2 = np.arcsin(np.sin(phi1) * np.cos(distance) + np.cos(phi1) * np.sin(distance) * np.cos(azimuth))
    lam2 = lam1 + np.arctan2(np.sin(azimuth) * np.sin(distance) * np.cos(phi1),
                             np.cos(distance) - np.sin(phi1) * np.sin(phi2))
    return np.degrees(phi2), np.degrees(lam2) % 360


def local_grid(name, lon0, lat0, half_width_m=72e3, spacing_m=100.0):
    """Heights above the atlas sea level on a local grid round (lon0, lat0), with its metadata."""
    entry = slice_for(name)
    axis = np.arange(-half_width_m, half_width_m + spacing_m / 2, spacing_m)
    x, y = np.meshgrid(axis, axis, indexing="xy")
    lat, lon = destination(lon0, lat0, x, y)
    radius_height = radius_heights(entry, lat, lon)
    # The degree-200 geoid is smooth over ~50 km: evaluate it on a 0.05-degree grid and interpolate.
    step = 0.05
    glat = np.arange(lat.max() + step, lat.min() - step, -step)
    glon = np.arange(lon.min() - step, lon.max() + step, step)
    datum = geoid(Grid(glat, glon, int(round(1 / step))), max_degree=200)
    rows = (glat[0] - lat) / step
    cols = (lon - glon[0]) / step
    geoid_height = map_coordinates(datum, [rows, cols], order=1, mode="nearest")
    with np.load(ATLAS, allow_pickle=False) as atlas:
        meta = json.loads(atlas["metadata"].item())
    level = float(meta["sea_level_m"])
    height = radius_height - geoid_height - level
    metadata = dict(
        schema="terluna.geography.local-coast/1",
        evidence=("Archived LOLA altimetry raster (256 pixels per degree, about 118 m spacing, with interpolation "
                  "where coverage is incomplete) above the degree-200 GRAIL geoid, relative to the atlas's common "
                  "sea level at its selected water share. Static flooded rock: beaches, sediment, erosion and "
                  "vegetation are not modelled."),
        reading_rule=("height_m is metres above the atlas sea level on a grid of metres east (columns) and north "
                      "(rows, increasing) of the centre, by the azimuthal equidistant projection; negative heights "
                      "are under water at that level. Add the tide's level to the sea, not to the ground."),
        name=name, centre_lon_lat_deg=[lon0, lat0], spacing_m=spacing_m, half_width_m=half_width_m,
        sea_level_m=level, source_slice_sha256=entry["selected_sha256"],
        manifest_sha256=digest(MANIFEST) if name != "smythii_headland" else digest(shore_terrain.MANIFEST),
        atlas_sha256=digest(ATLAS),
        producer=dict(domain="geography", files={str(p.relative_to(ROOT)): digest(p) for p in
                      (Path(__file__), ROOT / "geography/topography.py", ROOT / "geography/shore_terrain.py")}))
    return dict(x_m=axis, y_m=axis, height_m=height.astype(np.float32), metadata=json.dumps(metadata))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--coasts", nargs="*", help="name lon lat triples; default: the study's four coasts")
    args = parser.parse_args()
    restore(args.download)
    calendar = json.loads((ROOT / "research/studies/sea_appearance/results/lighting_calendar.json").read_text())
    names = {"W Procellarum": "w_procellarum", "Smythii headland": "smythii_headland", "S Nubium": "s_nubium",
             "Ingenii coast": "ingenii_coast"}
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for site in calendar["sites"]:
        if site["short"] not in names:
            continue
        name = names[site["short"]]
        # The shore studies' Smythii rows start at the equator, 2.3 degrees south of the headland.
        grid = local_grid(name, *site["lon_lat"], half_width_m=69e3 if name == "smythii_headland" else 72e3)
        path = OUTPUT / f"{name}.npz"
        np.savez_compressed(path, **grid)
        h = grid["height_m"]
        print(name, path, "land share", round(float((h > 0).mean()), 3), "highest", round(float(h.max()), 1))


if __name__ == "__main__":
    main()
