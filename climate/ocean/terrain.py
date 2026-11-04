"""The Near-side Sea's floor at 16 pixels per degree, for the ocean models.

Heights are LOLA's 16 px/deg radii above the degree-200 GRAIL geoid with the
static Earth tide (geography/topography.py, as the atlas and the sea-oxygen
study use them), in a box round the Near-side Sea: 36 S to 76 N and 90 W to
78 E, longitudes running continuously across 0. The water stands at the
atlas's common level, -1,654.2 m at 28% cover. The product is written once to
products/nearside_terrain_16ppd.npz (schema terluna.climate.ocean-terrain/1,
about 20 MB, kept out of Git); building it takes a minute or two and about a
gigabyte of memory.

Water masks join cells through shared faces (four neighbours), the way water
flows between grid cells, and keep the parts connected to named seed points.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PRODUCT = HERE / "products" / "nearside_terrain_16ppd.npz"
ATLAS = ROOT / "geography" / "products" / "atlas_28pct_4ppd.npz"
SCHEMA = "terluna.climate.ocean-terrain/1"
PPD = 16
BOX = dict(lat=(-36.0, 76.0), lon=(-90.0, 78.0))   # degrees; longitudes east, continuous across 0

# Points inside the seas, for picking connected water (lat, lon east in -180..180).
SEEDS = {"Mare Imbrium": (36.0, -16.0), "Mare Fecunditatis": (-6.0, 54.0), "Mare Nectaris": (-15.0, 34.5),
         "Mare Crisium": (17.0, 59.0), "Mare Nubium": (-21.0, -17.0), "Mare Tranquillitatis": (8.0, 31.0)}


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sea_level() -> float:
    with np.load(ATLAS, allow_pickle=False) as atlas:
        return float(json.loads(atlas["metadata"].item())["sea_level_m"])


def build(path=PRODUCT):
    """Cut the box from the global 16 px/deg heights above the geoid and write the product."""
    from geography import topography
    height, _, grid = topography.height_above_geoid(PPD, 200)
    lat = grid.lat_deg                       # north first
    lon = np.where(grid.lon_deg > 180.0, grid.lon_deg - 360.0, grid.lon_deg)
    order = np.argsort(lon, kind="stable")
    lon, height = lon[order], height[:, order]
    rows = np.flatnonzero((lat > BOX["lat"][0]) & (lat < BOX["lat"][1]))[::-1]   # south first
    cols = np.flatnonzero((lon > BOX["lon"][0]) & (lon < BOX["lon"][1]))
    box = height[np.ix_(rows, cols)].astype(np.float32)
    inputs = {name: sha256(topography.fetch_inputs.path(name)) for name in ("ldem_16.img", "ldem_4.img", "jggrx_0420a_sha.tab")}
    metadata = dict(schema=SCHEMA, pixels_per_degree=PPD, box=BOX, geoid_degree=200, sea_level_m=sea_level(),
                    producer={str(p.relative_to(ROOT)): sha256(p) for p in (Path(__file__), ROOT / "geography/topography.py")},
                    inputs=inputs, atlas_sha256=sha256(ATLAS),
                    evidence="LOLA 16 px/deg radii above the GRAIL GL0420A geoid to degree 200 with the static Earth tide "
                             "and rotation (geography/topography.py); the geoid is evaluated at 4 px/deg and interpolated.",
                    reading_rule="height_m[j, i] at lat_deg[j] (south first), lon_deg[i] (east, -180..180), metres above "
                                 "the geoid; water depth is sea_level_m minus height where that is positive.")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, height_m=box, lat_deg=lat[rows], lon_deg=lon[cols], metadata=np.array(json.dumps(metadata)))
    return Path(path)


def load(path=PRODUCT):
    """The terrain box: lat (south first), lon, height (m above the geoid), the sea level and the metadata."""
    if not Path(path).exists():
        build(path)
    with np.load(path, allow_pickle=False) as source:
        metadata = json.loads(source["metadata"].item())
        if metadata["schema"] != SCHEMA:
            raise ValueError(f"{path}: unsupported schema {metadata['schema']}")
        return dict(lat=source["lat_deg"].astype(float), lon=source["lon_deg"].astype(float),
                    height=source["height_m"].astype(float), level=float(metadata["sea_level_m"]), metadata=metadata)


def components(mask):
    """Labels of the regions of `mask` joined through shared faces (0 outside), and their count."""
    from scipy import ndimage
    return ndimage.label(mask, structure=[[0, 1, 0], [1, 1, 1], [0, 1, 0]])


def connected_water(height, level, lat, lon, seeds=("Mare Imbrium",)):
    """Cells below `level` joined through faces to any of the named seed points."""
    labels, _ = components(height < level)
    keep = set()
    for name in seeds:
        s_lat, s_lon = SEEDS[name]
        j = int(np.argmin(np.abs(lat - s_lat)))
        i = int(np.argmin(np.abs(lon - s_lon)))
        if labels[j, i] == 0:
            raise ValueError(f"seed {name} lies on land at this resolution")
        keep.add(int(labels[j, i]))
    return np.isin(labels, list(keep))


def depth(terrain, mask=None):
    """Water depth (m, positive) where `mask` (default: everything below the level), zero elsewhere."""
    d = terrain["level"] - terrain["height"]
    wet = d > 0 if mask is None else mask & (d > 0)
    return np.where(wet, d, 0.0)


def cell_area(lat, lon):
    """Areas (m^2) of cells centred on a regular lat/lon grid of this resolution, on the Moon's sphere."""
    from shared.constants import MOON_RADIUS
    dlat = np.radians(abs(lat[1] - lat[0]))
    dlon = np.radians(abs(lon[1] - lon[0]))
    phi = np.radians(lat)
    band = MOON_RADIUS**2 * dlon * (np.sin(phi + dlat / 2) - np.sin(phi - dlat / 2))
    return np.repeat(band[:, None], len(lon), axis=1)
