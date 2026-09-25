"""Static water storage on lunar topography: level curves and the depression tree.

Two views of how much water the Moon can hold before it floods:

- Level curve ("bathtub"): for a common equipotential water level, the area below
  it and the water needed to fill it. It describes water standing at one level
  everywhere; isolated depressions fill only if supplied.
- Depression tree: every closed depression, with the level at which it spills
  into a neighbour, its capacity and area at that level, found by merging grid
  cells in order of height (a union-find merge tree). The major merges show which
  basins join into one sea as water rises, and what each stage holds.

Heights are above the geoid (topography.height_above_geoid). Volumes are exact
sums over spherical grid cells. This is hydrostatics only: rainfall, runoff,
evaporation, groundwater, erosion and crustal loading are outside it.
"""
from __future__ import annotations
from dataclasses import dataclass
import math
import numpy as np

from shared.constants import MOON_RADIUS

AREA = 4 * math.pi * MOON_RADIUS**2
WATER_DENSITY = 1000.0

# Approximate centres of large named basins (IAU/USGS nomenclature), used only to
# label computed depressions; the computed geometry does not depend on them.
NAMED_BASINS = {
    'South Pole-Aitken': (-53.0, 191.0), 'Imbrium': (33.0, 342.0), 'Serenitatis': (27.0, 18.0),
    'Crisium': (17.0, 59.0), 'Tranquillitatis': (8.5, 31.0), 'Nectaris': (-15.0, 35.0),
    'Fecunditatis': (-8.0, 51.0), 'Nubium': (-21.0, 343.0), 'Humorum': (-24.0, 321.0),
    'Procellarum': (18.0, 303.0), 'Orientale': (-20.0, 265.0), 'Moscoviense': (27.0, 147.0),
    'Smythii': (1.0, 87.0), 'Frigoris': (56.0, 1.0), 'Hertzsprung': (2.0, 231.0),
    'Apollo': (-36.0, 209.0), 'Korolev': (-4.0, 203.0), 'Mendel-Rydberg': (-50.0, 266.0),
    'Schrodinger': (-75.0, 132.0), 'Australe': (-47.0, 94.0), 'Marginis': (13.0, 86.0),
    'Humboldtianum': (57.0, 82.0), 'Coulomb-Sarton': (52.0, 237.0), 'Freundlich-Sharonov': (19.0, 175.0),
}


def level_curve(height, area, levels):
    """Flooded area fraction and water volume (m^3) at each common water level."""
    h = height.ravel()
    a = area.ravel()
    order = np.argsort(h)
    hs, as_ = h[order], a[order]
    cum_area = np.cumsum(as_)
    cum_ah = np.cumsum(as_ * hs)
    k = np.searchsorted(hs, levels, side='right')
    area_below = np.where(k > 0, cum_area[np.maximum(k - 1, 0)], 0.0)
    ah_below = np.where(k > 0, cum_ah[np.maximum(k - 1, 0)], 0.0)
    volume = np.asarray(levels) * area_below - ah_below
    return area_below / a.sum(), volume


def level_for_area(height, area, fraction):
    """Common water level whose flooded area is the given fraction of the Moon (the height of
    the cell at which the cumulative area, in order of height, reaches that fraction)."""
    h = height.ravel()
    order = np.argsort(h)
    cum = np.cumsum(area.ravel()[order])
    k = int(np.searchsorted(cum, fraction * cum[-1]))
    return float(h[order][min(k, h.size - 1)])


def level_for_volume(height, area, volume):
    """Common water level holding a given volume (m^3), by bisection on the level curve."""
    lo, hi = float(height.min()), float(height.max())
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if level_curve(height, area, [mid])[1][0] < volume:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


@dataclass
class Depression:
    """A closed depression at the level where it spills into a neighbour."""
    ident: int
    spill_level_m: float
    floor_m: float
    area_m2: float
    volume_m3: float
    floor_lat: float
    floor_lon: float
    children: list
    names: tuple = ()        # named basin centres inside it
    event: int = -1          # merge event at which it spilled

    @property
    def gel_m(self):
        """Global equivalent layer: volume spread over the whole Moon."""
        return self.volume_m3 / AREA

    @property
    def area_fraction(self):
        return self.area_m2 / AREA


def _neighbours(i, j, nlat, nlon):
    for di in (-1, 0, 1):
        ii = i + di
        if ii < 0 or ii >= nlat:
            continue
        for dj in (-1, 0, 1):
            if di == 0 and dj == 0:
                continue
            yield ii * nlon + (j + dj) % nlon
    if i == 0 or i == nlat - 1:            # across the pole
        yield i * nlon + (j + nlon // 2) % nlon


def depression_tree(height, area, lat, lon):
    """All closed depressions and their spill levels (union-find merge tree).

    Cells are added from the lowest up. A cell touching no processed cell starts a
    depression; a cell joining two or more depressions is their spill point, and
    each joining depression is recorded with its capacity at that level. The merged
    depression continues and is recorded in turn when it spills. Returns the list
    of recorded depressions and, per cell, the id of the depression it first
    belonged to (for maps).
    """
    nlat, nlon = height.shape
    h = height.ravel()
    a = area.ravel()
    n = h.size
    order = np.argsort(h, kind='stable')
    parent = np.full(n, -1, dtype=np.int64)
    comp_area = np.zeros(n)
    comp_ah = np.zeros(n)
    comp_floor = np.zeros(n, dtype=np.int64)
    comp_node = {}                       # root -> list of recorded child depressions
    records = []
    events = []
    centres = {}
    for name, (la, lo) in NAMED_BASINS.items():
        ii = int(np.argmin(np.abs(lat - la)))
        jj = int(np.argmin(np.abs(((lon - lo) + 180) % 360 - 180)))
        centres[name] = ii * nlon + jj

    def find(x):
        root = x
        while parent[root] != root:
            root = parent[root]
        while parent[x] != root:
            parent[x], x = root, parent[x]
        return root

    for idx in order:
        i, j = divmod(int(idx), nlon)
        roots = set()
        for nb in _neighbours(i, j, nlat, nlon):
            if parent[nb] >= 0:
                roots.add(find(nb))
        level = h[idx]
        if not roots:
            parent[idx] = idx
            comp_area[idx], comp_ah[idx], comp_floor[idx] = a[idx], a[idx] * level, idx
            comp_node[idx] = []
            continue
        if len(roots) > 1:
            event = len(events)
            events.append([])
            for r in roots:
                vol = comp_area[r] * level - comp_ah[r]
                f = comp_floor[r]
                inside = tuple(nm for nm, cell in centres.items() if parent[cell] >= 0 and find(cell) == r)
                records.append(Depression(len(records), float(level), float(h[f]), float(comp_area[r]),
                                          float(vol), float(lat[f // nlon]), float(lon[f % nlon]),
                                          comp_node.pop(r), inside, event))
                events[event].append(records[-1])
        roots = sorted(roots, key=lambda r: comp_area[r], reverse=True)
        keep = roots[0]
        kids = [records[-len(roots) + k] for k in range(len(roots))] if len(roots) > 1 else comp_node[keep]
        for r in roots[1:]:
            parent[r] = keep
            comp_area[keep] += comp_area[r]
            comp_ah[keep] += comp_ah[r]
            if h[comp_floor[r]] < h[comp_floor[keep]]:
                comp_floor[keep] = comp_floor[r]
        parent[idx] = keep
        comp_area[keep] += a[idx]
        comp_ah[keep] += a[idx] * level
        comp_node[keep] = kids
    return records


def major_merges(records, min_area_fraction=2e-3):
    """Merge events where at least two depressions above min_area_fraction meet.

    Returns (level, [depressions]) sorted by level: the points where substantial
    bodies of water standing in separate basins would join into one.
    """
    events = {}
    for d in records:
        events.setdefault(d.event, []).append(d)
    out = []
    for parts in events.values():
        big = [d for d in parts if d.area_fraction >= min_area_fraction]
        if len(big) >= 2:
            out.append((big[0].spill_level_m, sorted(big, key=lambda d: -d.area_m2)))
    return sorted(out, key=lambda e: e[0])


def seas(height, area, level, min_area_fraction=1e-3):
    """Separate water bodies below a common level (8-connected, wrapping in longitude)."""
    from scipy import ndimage
    wet = height < level
    lab, count = ndimage.label(wet, structure=np.ones((3, 3), bool))
    # Join labels across the 0/360 longitude seam and across each pole.
    pairs = [(lab[:, 0], lab[:, -1])]
    n = wet.shape[1]
    pairs += [(lab[0], np.roll(lab[0], n // 2)), (lab[-1], np.roll(lab[-1], n // 2))]
    uf = np.arange(count + 1)

    def root(x):
        while uf[x] != x:
            uf[x] = uf[uf[x]]
            x = uf[x]
        return x
    for a_, b_ in pairs:
        for x, y in zip(a_, b_):
            if x and y:
                rx, ry = root(x), root(y)
                if rx != ry:
                    uf[ry] = rx
    roots = np.array([root(x) for x in range(count + 1)])
    merged = roots[lab]
    sizes = np.bincount(merged.ravel(), weights=area.ravel(), minlength=count + 1)
    sizes[0] = 0.0
    big = np.sort(sizes[sizes > min_area_fraction * AREA])[::-1]
    return big / AREA, merged


def label(dep: Depression, radius_deg=25.0):
    """Nearest named basin within radius_deg of the depression's floor, else coordinates."""
    best, best_d = None, radius_deg
    la1, lo1 = math.radians(dep.floor_lat), math.radians(dep.floor_lon)
    for name, (la, lo) in NAMED_BASINS.items():
        la2, lo2 = math.radians(la), math.radians(lo)
        d = math.degrees(math.acos(max(-1.0, min(1.0, math.sin(la1) * math.sin(la2)
                                                 + math.cos(la1) * math.cos(la2) * math.cos(lo1 - lo2)))))
        if d < best_d:
            best, best_d = name, d
    return best or f'{dep.floor_lat:+.0f}, {dep.floor_lon:.0f}E'
