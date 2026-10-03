"""Export a flooded LOLA terrain sector for coastal wave calculations."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import label

from geography import fetch_inputs
from geography.topography import Grid, geoid

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT / "research/runs/waves/coastal/terrain.npz"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sector(west=91.03125, south=-1.03125, width=4., height=4., ppd=16):
    """Retain native pixel centres and the water component meeting the west edge."""
    if ppd != 16 or min(width, height) <= 0 or width > 10 or height > 10:
        raise ValueError("Use a positive sector at most ten degrees across at 16 pixels/degree")
    manifest = fetch_inputs.manifest()
    required = ("ldem_16.img", "jggrx_0420a_sha.tab")
    source_hashes = {}
    for name in required:
        path = fetch_inputs.path(name)
        entry = next(f for f in manifest["files"] if f["name"] == name)
        digest = sha256(path)
        if digest != entry["sha256"] or path.stat().st_size != entry["bytes"]:
            raise ValueError(f"Terrain source bytes differ: {name}")
        source_hashes[name] = digest
    lon = west + np.arange(round(width * ppd) + 1) / ppd
    lat = south + np.arange(round(height * ppd) + 1) / ppd
    cols = (lon * ppd - .5).round().astype(int)
    rows = ((90 - lat) * ppd - .5).round().astype(int)
    if not np.allclose((cols + .5) / ppd, lon) or not np.allclose(90 - (rows + .5) / ppd, lat):
        raise ValueError("Sector origin must coincide with native LOLA centres")
    raw = np.memmap(fetch_inputs.path("ldem_16.img"), dtype="<i2", mode="r", shape=(180 * ppd, 360 * ppd))
    radius_height = np.asarray(raw[np.ix_(rows, cols)], dtype=float) * .5
    ground = radius_height - geoid(Grid(lat, lon, ppd), max_degree=200)
    atlas_path = ROOT / "geography/products/atlas_28pct_4ppd.npz"
    with np.load(atlas_path, allow_pickle=False) as atlas:
        meta = json.loads(atlas["metadata"].item())
        if meta["schema"] != "terluna.geography.atlas-grid/1":
            raise ValueError("Unsupported atlas schema")
        level = float(meta["sea_level_m"])
    physical_depth = level - ground
    components, _ = label(physical_depth > 0, structure=np.ones((3, 3)))
    west_labels = components[:, 0]
    counts = np.bincount(west_labels[west_labels > 0])
    if not len(counts):
        raise ValueError("The selected sea must meet the western boundary")
    selected = int(counts.argmax())
    wet = components == selected
    depth = np.where(wet, physical_depth, np.minimum(physical_depth, -1.))
    metadata = dict(schema="terluna.geography.coastal-grid/1",
                    evidence="Flooded lunar rock from LOLA at 16 pixels/degree and the degree-200 GRAIL geoid with static Earth tide and rotation. Beach formation and sediment evolution remain open.",
                    reading_rule="Latitudes increase northward and longitudes eastward. Nodes retain native LOLA centres. Depth is the atlas sea level minus ground; the water component with the most western-edge nodes is selected. Other components are held dry. Geoid harmonics are evaluated directly on this sector.",
                    units=dict(longitude_deg="degrees east", latitude_deg="degrees north", depth_m="m", height_m="m above geoid"),
                    sea_level_m=level, pixels_per_degree=ppd, bounds_deg=[west, south, width, height],
                    wet_nodes=int(wet.sum()), excluded_wet_nodes=int(((physical_depth > 0) & ~wet).sum()),
                    source_sha256=dict(**source_hashes, atlas=sha256(atlas_path)),
                    credit={f["name"]: f["credit"] for f in manifest["files"] if f["name"] in required},
                    producer=dict(domain="geography", files={str(p.relative_to(ROOT)): sha256(p) for p in
                                  (Path(__file__), ROOT / "geography/topography.py", ROOT / "shared/constants.json")}))
    return dict(longitude_deg=lon, latitude_deg=lat, depth_m=depth, height_m=ground, wet=wet, metadata=json.dumps(metadata))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT)
    args = parser.parse_args()
    data = sector()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, **data)
    print(data["metadata"])


if __name__ == "__main__":
    main()
