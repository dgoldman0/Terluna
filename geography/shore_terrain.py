"""Restore a pinned LOLA row slice and export native 118 m coastal terrain."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
from scipy.ndimage import label

from geography import fetch_inputs
from geography.topography import Grid, geoid

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = Path(__file__).with_name("shore_inputs.json")
OUTPUT = ROOT/"research/runs/waves/shore"


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def restore(download=False):
    """Accept the recorded byte slice only, with its exact label and hash."""
    source = json.loads(MANIFEST.read_text())
    for kind in ("label", "image"):
        path = ROOT/source["local_"+kind]
        expected = source["label_sha256" if kind == "label" else "selected_sha256"]
        size = source["label_bytes" if kind == "label" else "selected_bytes"]
        if path.exists():
            if path.stat().st_size != size or digest(path) != expected:
                raise ValueError(f"Pinned coastal input differs: {path}")
            continue
        if not download:
            raise FileNotFoundError(f"Restore the coastal input with --download: {path}")
        path.parent.mkdir(parents=True, exist_ok=True)
        partial = path.with_suffix(path.suffix+".partial")
        if partial.exists():
            raise FileExistsError(f"Inspect the retained partial download: {partial}")
        headers = {}
        url = source["label_url"]
        if kind == "image":
            first, last = source["selected_byte_range"]
            headers = {"Range":f"bytes={first}-{last}", "If-Range":source["source_etag"]}
            url = source["source_url"]
        with urlopen(Request(url, headers=headers), timeout=90) as response:
            if kind == "image" and (response.status != 206 or response.headers.get("Content-Range") !=
                    f"bytes {first}-{last}/{source['source_size_bytes']}"):
                raise ValueError("The terrain server must return the exact requested byte range")
            with partial.open("xb") as stream:
                remaining = size
                while remaining:
                    chunk = response.read(min(1024**2, remaining))
                    if not chunk:
                        raise ValueError("Incomplete terrain download")
                    stream.write(chunk)
                    remaining -= len(chunk)
                if response.read(1):
                    raise ValueError("Terrain response exceeds the pinned byte count")
        if digest(partial) != expected:
            raise ValueError(f"Terrain checksum differs; retained for inspection: {partial}")
        partial.rename(path)
    return source


def native_indices(bounds, source):
    """Native pixel centres, increasing east and north, and their file indices."""
    west, south, width, height = bounds
    ppd = source["pixels_per_degree"]
    if not np.isfinite(bounds).all() or min(width, height) <= 0:
        raise ValueError("Finite positive coastal bounds are required")
    cells = np.array([width, height])*ppd
    if not np.allclose(cells, np.rint(cells), rtol=0, atol=1e-8):
        raise ValueError("Coastal extents must contain complete native cells")
    lon = west+np.arange(round(cells[0])+1)/ppd
    lat = south+np.arange(round(cells[1])+1)/ppd
    columns_float = lon*ppd-.5
    rows_float = (90-lat)*ppd-.5
    if not (np.allclose(columns_float,np.rint(columns_float),rtol=0,atol=1e-8)
            and np.allclose(rows_float,np.rint(rows_float),rtol=0,atol=1e-8)):
        raise ValueError("Coastal origins must coincide with native pixel centres")
    columns = np.rint(columns_float).astype(int)
    rows = np.rint(rows_float).astype(int)-source["first_row"]
    if columns.min() < 0 or columns.max() >= source["samples"] or rows.min() < 0 or rows.max() >= source["row_count"]:
        raise ValueError("Coastal sector lies outside the restored rows")
    return lon, lat, rows, columns


def flooded_component(physical_depth):
    """Select the four-connected water component with the most west-edge nodes."""
    components, _ = label(physical_depth > 0)
    west = components[:,0]
    counts = np.bincount(west[west > 0])
    if not len(counts):
        raise ValueError("The sea must meet the western boundary")
    wet = components == int(counts.argmax())
    return wet, np.where(wet, physical_depth, np.minimum(physical_depth, -1.))


def sector(bounds=(92.751953125, 1.751953125, 1., 1.)):
    source = restore()
    lon, lat, rows, columns = native_indices(bounds, source)
    raw = np.memmap(ROOT/source["local_image"], dtype=source["sample_type"], mode="r",
                    shape=(source["row_count"],source["samples"]))
    radius_height = np.asarray(raw[np.ix_(rows, columns)],dtype=float)*source["scale_m"]
    gravity = fetch_inputs.path("jggrx_0420a_sha.tab")
    expected = next(f for f in fetch_inputs.manifest()["files"] if f["name"] == gravity.name)
    if gravity.stat().st_size != expected["bytes"] or digest(gravity) != expected["sha256"]:
        raise ValueError("The GRAIL source differs from its pinned record")
    datum = geoid(Grid(lat,lon,source["pixels_per_degree"]),max_degree=200)
    ground = radius_height-datum
    atlas_path = ROOT/"geography/products/atlas_28pct_4ppd.npz"
    with np.load(atlas_path,allow_pickle=False) as atlas:
        atlas_meta = json.loads(atlas["metadata"].item())
    if atlas_meta["schema"] != "terluna.geography.atlas-grid/1":
        raise ValueError("Unsupported atlas schema")
    level = float(atlas_meta["sea_level_m"])
    physical_depth = level-ground
    wet, depth = flooded_component(physical_depth)
    metadata = dict(schema="terluna.geography.coastal-grid/1",
        evidence="Static flooded rock from the archived 256-pixel/degree LOLA raster and degree-200 GRAIL geoid. Native pixel spacing is about 118 m; the source interpolates gaps in altimetry coverage. Beach sediment and erosion remain separate models.",
        reading_rule="Longitude increases eastward and latitude northward at native pixel centres. The four-connected water component with the most western-edge nodes is selected; isolated water is held dry. All wave-grid refinements sample this same bottom grid. The inherited geoid calculation retains its mean-Earth/principal-axis frame approximation.",
        units=dict(longitude_deg="degrees east",latitude_deg="degrees north",depth_m="m",height_m="m above geoid"),
        sea_level_m=level,pixels_per_degree=source["pixels_per_degree"],bounds_deg=list(bounds),
        wet_nodes=int(wet.sum()),excluded_wet_nodes=int(((physical_depth>0)&~wet).sum()),
        source_sha256=dict(selected_lola_bytes=source["selected_sha256"],label=source["label_sha256"],
                           manifest=digest(MANIFEST),grail=digest(gravity),atlas=digest(atlas_path)),
        source=source,credit=source["credit"],
        producer=dict(domain="geography",files={str(p.relative_to(ROOT)):digest(p) for p in
            (Path(__file__),ROOT/"geography/topography.py",ROOT/"shared/constants.json")}))
    return dict(longitude_deg=lon,latitude_deg=lat,depth_m=depth,height_m=ground,
                radius_height_m=radius_height,geoid_height_m=datum,wet=wet,metadata=json.dumps(metadata))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download",action="store_true")
    parser.add_argument("--bounds",type=float,nargs=4,default=(92.751953125,1.751953125,1.,1.))
    parser.add_argument("--output",type=Path,default=OUTPUT/"terrain.npz")
    args = parser.parse_args()
    restore(args.download)
    data = sector(args.bounds)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(args.output,**data)
    print(json.dumps(dict(output=str(args.output),sha256=digest(args.output),
                         shape=data["depth_m"].shape,wet_nodes=int(data["wet"].sum())),indent=2))


if __name__ == "__main__":
    main()
