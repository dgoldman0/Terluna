"""The crust's sag under the seas and lakes, and what it does to the sea level and the straits.

    python -m geography.crustal_loading     # results/crustal_loading.json; grids on the research drive

The water's weight bends the crust down under the seas and the deep rain-fed lakes.
The response comes from love_numbers.py: the immediate elastic response of the
published interior, and the relaxed response of an elastic lithosphere over an
interior that has flowed, bracketed over the elastic thicknesses published for the
mare basins. The sea surface is an equipotential of the deformed Moon and its water,
found by the sea-level equation (Farrell & Clark 1976) with the shoreline free to
move (Mitrovica & Milne 2003): the water thickness is the sea surface less the
deformed floor wherever that is positive, one common level holds either today's
volume or 28% of the surface, and the load, the deformation and the shoreline are
iterated together. The lakes are a fixed load at the drainage product's depths.

Spherical harmonics are evaluated per order on the 4 px/deg grid (the load averaged
from 16 px/deg); the deformation is interpolated back to 16 px/deg, where levels,
areas, sills and lakes are measured. Degree 0 is a uniform shift, absorbed in the level.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import time
import numpy as np
from scipy import ndimage

from geography import basins as bs, fetch_inputs, love_numbers as ln, topography as tp
from shared.constants import (GAS_CONSTANT, GRAVITATIONAL_CONSTANT, MOON_GM, MOON_RADIUS, MOON_SURFACE_GRAVITY,
                              STANDARD_ATMOSPHERE)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESULT = HERE / "results/crustal_loading.json"
BULK = Path("/media/projectspace/terluna-research/crustal-loading")
ATLAS_JSON = HERE / "results/atlas.json"
DRAINAGE_JSON = HERE / "results/drainage.json"
DRAINAGE_GRID = HERE / "products/drainage_28pct_16ppd.npz"
UPTAKE_JSON = HERE / "results/groundwater_uptake.json"
SCHEMA = "terluna.geography.crustal-loading/1"

WATER_DENSITY = bs.WATER_DENSITY
MEAN_DENSITY = 3 * MOON_GM / (4 * math.pi * GRAVITATIONAL_CONSTANT * MOON_RADIUS**3)
SHARE = 0.28
TE_CASES_KM = (20.0, 40.0, 75.0)       # the published range under the mare basins (sources file)
THICK_TE_KM = 150.0                     # a colder, thicker lithosphere today: a sensitivity, not a source value
YOUNG_PA, POISSON = 1.0e11, 0.25
VISCOSITIES_PA_S = (1e20, 1e21, 1e22)
NMAX = 359
NMAX_CHECK = (179, 719)
TRANSFORM_PPD = 4
FINE_PPD = 16
TOLERANCE_M = 0.01
AIR_PRESSURE_ATM = 1.2                  # at the sea surface, design value
AIR_TEMPERATURE_K = 290.0               # mean temperature of the air column, for its scale height
AIR_MOLAR_MASS = 0.02896                # Earthlike air, kg/mol

# The strait study's zones (climate/ocean/fecunditatis.py) and the seeds of the basins they join.
ZONES = {"western channels": dict(lon=(35.0, 44.5), lat=(-11.5, -0.5)),
         "Crisium passage": dict(lon=(53.0, 62.0), lat=(6.0, 11.5))}
STRAITS = {
    "Fecunditatis to the Near-side Sea (western channels)":
        dict(box=dict(lon=(25.0, 58.0), lat=(-25.0, 5.0)), a=(-7.8, 53.7), b=(-15.2, 34.6)),
    "Fecunditatis to Crisium (Crisium passage)":
        dict(box=dict(lon=(45.0, 70.0), lat=(-20.0, 25.0)), a=(-7.8, 53.7), b=(16.2, 59.1)),
    "Tenebrarum sub-basin to the South Pole-Aitken Sea":
        dict(box=dict(lon=(60.0, 260.0), lat=(-90.0, -25.0)), a=None, b=None),
}
NEARSIDE_MARIA = {"Imbrium": (34.7, 345.1), "Procellarum": (18.0, 303.0), "Serenitatis": (27.0, 18.0),
                  "Crisium": (16.2, 59.1), "Nectaris": (-15.2, 34.6), "Fecunditatis": (-7.8, 53.7),
                  "Nubium": (-21.0, 343.0), "Humorum": (-24.0, 321.0), "Cognitum": (-10.5, 337.7),
                  "Frigoris": (56.0, 1.0)}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ------------------------------------------------------------- harmonics

class Filter:
    """Isotropic spherical-harmonic filters on a cell-centred latitude-longitude grid.

    Analysis weights each cell by its solid angle; synthesis evaluates the filtered
    coefficients per order, so no table of Legendre functions is held."""

    def __init__(self, nlat, nlon, nmax):
        if nmax > nlon // 2 - 1:
            raise ValueError("Degree above the grid's longitude resolution")
        self.nlat, self.nlon, self.nmax = nlat, nlon, nmax
        dlat, dlon = math.pi / nlat, 2 * math.pi / nlon
        lat = math.pi / 2 - (np.arange(nlat) + 0.5) * dlat
        self.t, self.u = np.sin(lat), np.cos(lat)
        self.weight = (np.sin(lat + dlat / 2) - np.sin(lat - dlat / 2)) * dlon / (4 * math.pi)
        m = np.arange(nlon // 2 + 1)
        self.shift = np.exp(-0.5j * m * dlon)          # the cells' half-cell offset in longitude

    def apply(self, field, factors):
        """Each row of `factors` (degree 0..nmax) applied to `field`; returns one field per row."""
        factors = np.atleast_2d(factors)
        F = np.fft.rfft(field, axis=1) * self.shift
        out = np.zeros((len(factors), self.nlat, self.nlon // 2 + 1), dtype=complex)
        for m in range(self.nmax + 1):
            p = tp.legendre_column(m, self.nmax, self.t, self.u)
            coeff = p @ (self.weight * F[:, m])
            for k, f in enumerate(factors):
                out[k, :, m] = (f[m:] * coeff) @ p
        out *= np.conj(self.shift)
        out[:, :, 0] *= 2.0
        return np.fft.irfft(out * (self.nlon / 2), n=self.nlon, axis=2)


def block_mean(field, factor):
    """Area-weighted means over factor x factor blocks of a cell-centred grid."""
    nlat, nlon = field.shape
    lat = 90.0 - (np.arange(nlat) + 0.5) * 180.0 / nlat
    w = np.cos(np.radians(lat)).reshape(nlat // factor, factor, 1)
    blocks = field.reshape(nlat // factor, factor, nlon // factor, factor).sum(axis=3)
    return (blocks * w).sum(axis=1) / (w.sum(axis=1) * factor)


def upsample(field, factor, chunk=360):
    """Bilinear interpolation of a coarse cell-centred grid onto one `factor` times finer."""
    nlat, nlon = field.shape
    pos = (np.arange(nlat * factor) + 0.5) / factor - 0.5
    i0 = np.clip(np.floor(pos).astype(int), 0, nlat - 1)
    i1 = np.clip(i0 + 1, 0, nlat - 1)
    fi = np.clip(pos - np.floor(pos), 0.0, 1.0)[:, None]
    pos = (np.arange(nlon * factor) + 0.5) / factor - 0.5
    j0 = np.floor(pos).astype(int)
    fj = (pos - j0)[None, :]
    out = np.empty((nlat * factor, nlon * factor))
    for start in range(0, nlat * factor, chunk):
        s = slice(start, start + chunk)
        rows = field[i0[s]] * (1 - fi[s]) + field[i1[s]] * fi[s]
        out[s] = rows[:, j0 % nlon] * (1 - fj) + rows[:, (j0 + 1) % nlon] * fj
    return out


# ------------------------------------------------------------- the sea-level equation

class Response:
    """Ground (U) and equipotential (N) heights from a surface load, for one set of load Love numbers."""

    def __init__(self, h, k, filt):
        n = np.arange(filt.nmax + 1)
        scale = 3.0 / ((2 * n + 1) * MEAN_DENSITY)
        u, e = h[:filt.nmax + 1] * scale, (1.0 + k[:filt.nmax + 1]) * scale
        u[0] = e[0] = 0.0
        self.factors = np.vstack([u, e])
        self.filt = filt

    def __call__(self, load_coarse):
        U, N = self.filt.apply(load_coarse, self.factors)
        return U, N


def level_for_volume(eff, row_area, volume, guess):
    level = guess
    for _ in range(30):
        depth = np.maximum(level - eff, 0.0)
        v = float(row_area @ depth.sum(axis=1))
        a = float(row_area @ np.count_nonzero(depth > 0, axis=1))
        step = (volume - v) / max(a, 1.0)
        level += step
        if abs(step) < 1e-4:
            break
    return level


def level_for_area(eff, row_area, fraction, guess):
    total = row_area.sum() * eff.shape[1]
    lo, hi = guess - 2000.0, guess + 2000.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if row_area @ np.count_nonzero(eff < mid, axis=1) < fraction * total:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def solve(height, row_area, response, fixed_load, level0, volume=None, fraction=None, iterations=60, dev=None):
    """Sea-level equation on the fine grid: returns the level, the effective heights' deformation
    (U - N on the coarse grid), U and N, and the iteration history."""
    coarse = (height.shape[0] // TRANSFORM_FACTOR, height.shape[1] // TRANSFORM_FACTOR)
    dev = np.zeros(coarse) if dev is None else dev
    level, history = level0, []
    for it in range(iterations):
        eff = height + upsample(dev, TRANSFORM_FACTOR)
        if fraction is not None:
            level = level_for_area(eff, row_area, fraction, level)
        else:
            level = level_for_volume(eff, row_area, volume, level)
        load = WATER_DENSITY * np.maximum(level - eff, 0.0)
        del eff
        if fixed_load is not None:
            load += fixed_load
        U, N = response(block_mean(load, TRANSFORM_FACTOR))
        del load
        change = float(np.max(np.abs(U - N - dev)))
        dev = U - N
        history.append(round(change, 4))
        if change < TOLERANCE_M:
            break
    else:
        raise RuntimeError(f"Sea-level equation did not converge: {history[-5:]}")
    return dict(level=level, dev=dev, U=U, N=N, history=history)


TRANSFORM_FACTOR = FINE_PPD // TRANSFORM_PPD


# ------------------------------------------------------------- measures on the fine grid

def grid_index(lat, lon, ppd=FINE_PPD):
    return int((90.0 - lat) * ppd), int((lon % 360.0) * ppd)


def box_slices(box, ppd=FINE_PPD):
    rows = slice(int((90.0 - box["lat"][1]) * ppd), int((90.0 - box["lat"][0]) * ppd))
    cols = slice(int((box["lon"][0] % 360.0) * ppd), int((box["lon"][1] % 360.0) * ppd) or None)
    return rows, cols


def deepest_near(height, lat, lon, radius_deg=3.0, ppd=FINE_PPD):
    i, j = grid_index(lat, lon, ppd)
    r = int(radius_deg * ppd)
    window = height[max(i - r, 0):i + r, max(j - r, 0):j + r]
    a, b = np.unravel_index(np.argmin(window), window.shape)
    return max(i - r, 0) + a, max(j - r, 0) + b


def pass_level(height, seeds, lo, hi, structure=None):
    """Lowest water level at which all seed cells lie in one body of water (bisection on
    connected components; four-connected unless a structure is given)."""
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        lab, _ = ndimage.label(height < mid, structure=structure)
        ids = {lab[s] for s in seeds}
        if len(ids) == 1 and 0 not in ids:
            hi = mid
        else:
            lo = mid
    return hi


def strait_sills(eff, level, seeds):
    """Sill depth below the sea surface for each strait: the level at which its two basins join."""
    out = {}
    for name, spec in STRAITS.items():
        rows, cols = box_slices(spec["box"])
        sub = eff[rows, cols]
        sa, sb = seeds[name]
        sa = (sa[0] - rows.start, sa[1] - cols.start)
        sb = (sb[0] - rows.start, sb[1] - cols.start)
        sill = pass_level(sub, [sa, sb], float(min(sub[sa], sub[sb])), level + 500.0)
        out[name] = dict(sill_level_m=round(sill, 2), sill_depth_m=round(level - sill, 2))
    return out


def nearside_join(eff, level, seeds):
    """Level at which every near-side mare seed lies in Imbrium's water, on a grid rolled so the
    near side is contiguous; the margin is the sea level above it."""
    shift = eff.shape[1] // 2
    rolled = np.roll(eff, shift, axis=1)
    rows = slice(int((90 - 75) * FINE_PPD), int((90 + 50) * FINE_PPD))
    cols = slice(int(80 * FINE_PPD), int(280 * FINE_PPD))        # 100 W to 100 E
    sub = rolled[rows, cols]
    pts = {k: (i - rows.start, (j + shift) % eff.shape[1] - cols.start) for k, (i, j) in seeds.items()}
    join = pass_level(sub, list(pts.values()), level - 1500.0, level + 500.0)
    lab, _ = ndimage.label(sub < join - 0.5)
    last = sorted(k for k, p in pts.items() if lab[p] != lab[pts["Imbrium"]])
    return dict(join_level_m=round(join, 2), margin_m=round(level - join, 2), last_to_join=last)


def shore_measures(height, eff, level0, level, row_area, cell_km):
    """Sea area, the relative sea-level change along today's shore and the shoreline's mean shift."""
    total = row_area.sum() * height.shape[1]
    wet0 = height < level0
    wet = eff < level
    area0 = float(row_area @ np.count_nonzero(wet0, axis=1))
    area1 = float(row_area @ np.count_nonzero(wet, axis=1))
    shore = wet0 & ~ndimage.binary_erosion(wet0, structure=np.ones((3, 3), bool), border_value=1)
    rsl = (level - eff) - (level0 - height)
    values = rsl[shore]
    lat = 90.0 - (np.nonzero(shore)[0] + 0.5) / FINE_PPD
    length_km = float(np.sum(cell_km * np.sqrt(0.5 * (1 + np.cos(np.radians(lat))**2))))
    return dict(sea_share=round(area1 / total, 5), sea_share_change=round((area1 - area0) / total, 5),
                land_share=round(1 - area1 / total, 5),
                shore_rsl_change_m=dict(median=round(float(np.median(values)), 2),
                                        p05=round(float(np.percentile(values, 5)), 2),
                                        p95=round(float(np.percentile(values, 95)), 2)),
                shoreline_mean_shift_km=round((area1 - area0) / 1e6 / max(length_km, 1.0), 3),
                shoreline_length_km=round(length_km, 0))


# ------------------------------------------------------------- loads

def lake_load(shape):
    with np.load(DRAINAGE_GRID, allow_pickle=False) as d:
        meta = json.loads(d["metadata"].item())
        cells, levels, spills = d["lake_cell"], d["lake_level_m"].astype(float), d["lake_spills"]
    if meta["grid"]["rows"] != shape[0] or meta["grid"]["cols"] != shape[1]:
        raise ValueError("Drainage grid does not match the heights")
    return cells, levels, spills, meta


def groundwater_load(height, level0, lake_cells, years="10000", case="central"):
    """Water held in the crust's pores (kg/m^2) by column type and latitude band at one time."""
    uptake = json.loads(UPTAKE_JSON.read_text())
    bands = uptake["cases"][case]["bands"]
    nlat, nlon = height.shape
    lat = np.abs(90.0 - (np.arange(nlat) + 0.5) / FINE_PPD)
    load = np.zeros_like(height)
    sea = height < level0
    lake = np.zeros(height.size, bool)
    lake[lake_cells] = True
    lake = lake.reshape(height.shape)
    for name, band in bands.items():
        lo, hi = (float(x) for x in name.split("-"))
        rows = (lat >= lo) & (lat < hi if hi < 90 else lat <= hi)
        held = {col: band["columns"][col]["by_year"][years]["held_m"] for col in band["columns"]}
        # A band without lakes in the uptake study takes its land column for any lake cells.
        block = np.where(sea[rows], held["sea_floor"], np.where(lake[rows], held.get("lake_floor", held["land"]), held["land"]))
        load[rows] = WATER_DENSITY * block
    return load


def air_load(height, level0):
    """Weight of an isothermal column of air (kg/m^2) over the ground or the sea surface."""
    g = MOON_SURFACE_GRAVITY
    scale = GAS_CONSTANT * AIR_TEMPERATURE_K / (AIR_MOLAR_MASS * g)
    surface = np.maximum(height, level0)
    return AIR_PRESSURE_ATM * STANDARD_ATMOSPHERE / g * np.exp(-(surface - level0) / scale), scale


# ------------------------------------------------------------- lakes

def lake_changes(height, dev_fine, U_fine, cells, levels, spills, row_area, count=12):
    """For the largest rain-fed lakes: the floor's sag against the outlet's, the spill level's change
    relative to the floor and the volume gained with the shoreline held."""
    mask = np.zeros(height.size, bool)
    mask[cells] = True
    level_grid = np.full(height.size, np.nan)
    level_grid[cells] = levels
    spill_grid = np.zeros(height.size, bool)
    spill_grid[cells] = spills
    mask, level_grid, spill_grid = (x.reshape(height.shape) for x in (mask, level_grid, spill_grid))
    lab, nlab = ndimage.label(mask, structure=np.ones((3, 3), bool))
    area = np.broadcast_to(row_area[:, None], height.shape)
    sizes = ndimage.sum(area, lab, index=np.arange(1, nlab + 1))
    order = np.argsort(sizes)[::-1][:count]
    objs = ndimage.find_objects(lab)
    out, total_gain = [], 0.0
    for idx in order:
        sl = objs[idx]
        sl = (slice(max(sl[0].start - 2, 0), sl[0].stop + 2), slice(max(sl[1].start - 2, 0), min(sl[1].stop + 2, height.shape[1])))
        inside = lab[sl] == idx + 1
        ring = ndimage.binary_dilation(inside, structure=np.ones((3, 3), bool)) & ~inside & ~mask[sl]
        h, d, a = height[sl], dev_fine[sl], area[sl]
        level = float(np.nanmedian(level_grid[sl][inside]))
        depth = np.where(inside, level - h, 0.0)
        deepest = np.unravel_index(np.argmax(depth), depth.shape)
        outlet = np.unravel_index(np.argmin(np.where(ring, h, np.inf)), h.shape)
        spill = bool(spill_grid[sl][inside].any())
        gain = float(np.sum(a[inside] * (d[outlet] - d[inside]))) if spill else None
        if gain is not None:
            total_gain += gain
        i, j = sl[0].start + deepest[0], sl[1].start + deepest[1]
        out.append(dict(centre=[round(90 - (i + 0.5) / FINE_PPD, 2), round(((j + 0.5) / FINE_PPD + 180) % 360 - 180, 2)],
                        area_km2=round(float(sizes[idx]) / 1e6, 0), max_depth_m=round(float(depth[deepest]), 0),
                        spills=spill, floor_fall_m=round(float(-U_fine[sl][deepest]), 1),
                        floor_fall_below_equipotential_m=round(float(-d[deepest]), 1),
                        outlet_fall_below_equipotential_m=round(float(-d[outlet]), 1),
                        spill_level_rise_over_floor_m=round(float(d[outlet] - d[deepest]), 1),
                        volume_gain_km3=None if gain is None else round(gain / 1e9, 0)))
    return out, total_gain / 1e9


# ------------------------------------------------------------- the run

def love_sets(nmax):
    """Load Love numbers for each case, and the validation of the interior's tidal response."""
    sets, notes = {}, {}
    comp = ln.interior_model(compressible=True)
    inc = ln.interior_model(compressible=False)
    sets["rigid"] = (np.zeros(nmax + 1), np.zeros(nmax + 1))
    h, k, rows = ln.elastic_spectrum(ln.Interior(comp), nmax)
    sets["elastic"] = (h, k)
    tides = {}
    for name, layers in (("compressible, Poisson's ratio 0.25", comp), ("incompressible, as the source", inc)):
        r = ln.elastic_love_numbers(ln.Interior(layers), 2)
        tides[name] = dict(k2=round(r["k"], 5), h2=round(r["h"], 5))
    # The same interior with the low-viscosity zone relaxed (its Maxwell time is about nine days).
    soft = [ln.Layer(l.top_m, l.density, 0.0 if i == 1 else l.rigidity, None if i == 1 else l.lame)
            for i, l in enumerate(comp)]
    r = ln.elastic_love_numbers(ln.Interior(soft), 2)
    tides["compressible, low-viscosity zone relaxed"] = dict(k2=round(r["k"], 5), h2=round(r["h"], 5))
    hs, ks, _ = ln.elastic_spectrum(ln.Interior(soft), nmax, samples=12)
    sets["elastic, low-viscosity zone relaxed"] = (hs, ks)
    mass, inertia = ln.mass_and_inertia(comp)
    notes["interior"] = dict(mass_kg=mass, mass_over_gm_g=round(mass * GRAVITATIONAL_CONSTANT / MOON_GM, 5),
                             moment_of_inertia=round(inertia, 5), tidal=tides,
                             elastic_load_love=[dict(n=n, h_load=round(rows[n]["h_load"], 5), k_load=round(rows[n]["k_load"], 6))
                                                for n in sorted(rows) if n in (1, 2, 3, 5, 10, 20, 50, 100, 200, 500) or n == max(rows)])
    for te in TE_CASES_KM + (THICK_TE_KM,):
        hr, kr = ln.relaxed_love_numbers(comp, nmax, te * 1e3, YOUNG_PA, POISSON)
        # The shell's own compression, absent from thin-shell flexure, is the elastic response;
        # adding it counts the deep interior's elastic part twice at the lowest degrees (a few percent).
        sets[f"relaxed, Te {te:g} km"] = (hr + h, kr + k)
        notes[f"relaxed, Te {te:g} km"] = [dict(n=n, h_load=round(float(hr[n] + h[n]), 4), k_load=round(float(kr[n] + k[n]), 4))
                                          for n in (1, 2, 3, 5, 10, 20, 30, 50, 100, 200) if n <= nmax]
    return sets, notes


def relaxation_times():
    mantle = ln.interior_model()[-2].density
    rows = []
    for te in TE_CASES_KM:
        for eta in VISCOSITIES_PA_S:
            rows.append(dict(te_km=te, viscosity_pa_s=eta, years={str(n): round(ln.relaxation_time_years(
                n, eta, te * 1e3, YOUNG_PA, POISSON, mantle, MOON_SURFACE_GRAVITY), -1) for n in (3, 5, 10, 20, 50)}))
    return rows


def run(nmax=NMAX, cases=None, extra=True, targets=("volume", "share")):
    started = time.time()
    height, _, grid = tp.height_above_geoid(FINE_PPD, 200)
    row_area = grid.cell_area_m2[:, 0].copy()
    atlas = json.loads(ATLAS_JSON.read_text())
    level0 = float(atlas["sea_level_m"])
    wet0 = height < level0
    volume0 = float(row_area @ np.maximum(level0 - height, 0.0).sum(axis=1))
    del wet0
    cells, levels, spills, _ = lake_load(height.shape)
    lakes = np.zeros(height.size)
    lakes[cells] = WATER_DENSITY * np.maximum(levels - height.ravel()[cells], 0.0)
    lakes = lakes.reshape(height.shape)
    filt = Filter(height.shape[0] // TRANSFORM_FACTOR, height.shape[1] // TRANSFORM_FACTOR, nmax)
    sets, notes = love_sets(max(nmax, 1))
    if cases:
        sets = {k: v for k, v in sets.items() if k in cases}
    seeds = {name: (deepest_near(height, *spec["a"]), deepest_near(height, *spec["b"]))
             for name, spec in STRAITS.items() if spec["a"]}
    spa_name = "Tenebrarum sub-basin to the South Pole-Aitken Sea"
    rows, cols = box_slices(dict(lat=(-86.0, -70.0), lon=(96.0, 161.0)))
    sub = height[rows, cols]
    i, j = np.unravel_index(np.argmin(sub), sub.shape)
    rows2, cols2 = box_slices(dict(lat=(-70.0, -30.0), lon=(150.0, 230.0)))
    sub2 = height[rows2, cols2]
    i2, j2 = np.unravel_index(np.argmin(sub2), sub2.shape)
    seeds[spa_name] = ((rows.start + i, cols.start + j), (rows2.start + i2, cols2.start + j2))
    mare_seeds = {k: deepest_near(height, *v) for k, v in NEARSIDE_MARIA.items()}
    # The Near-side Sea's deepest point: the lowest cell of the water joined to Imbrium at today's level.
    lab, _ = ndimage.label(height < level0)
    body = lab == lab[mare_seeds["Imbrium"]]
    nearside_deep = np.unravel_index(np.argmin(np.where(body, height, np.inf)), height.shape)
    del lab, body
    cell_km = 2 * math.pi * MOON_RADIUS / 1e3 / (360 * FINE_PPD)
    base = dict(sills=strait_sills(height, level0, seeds), nearside=nearside_join(height, level0, mare_seeds),
                shore=shore_measures(height, height, level0, level0, row_area, cell_km),
                volume_km3=round(volume0 / 1e9, 0), global_layer_m=round(volume0 / (4 * math.pi * MOON_RADIUS**2), 2))
    results, saved = {}, {}
    BULK.mkdir(parents=True, exist_ok=True)
    for name, (h, k) in sets.items():
        response = Response(h, k, filt)
        out = {}
        for target in targets:
            t0 = time.time()
            sol = solve(height, row_area, response, lakes, level0,
                        volume=volume0 if target == "volume" else None, fraction=SHARE if target == "share" else None)
            dev = upsample(sol["dev"], TRANSFORM_FACTOR)
            eff = height + dev
            level = sol["level"]
            depth0 = np.maximum(level0 - height, 0.0)
            vol = float(row_area @ np.maximum(level - eff, 0.0).sum(axis=1))
            sea = depth0 > 0
            U = upsample(sol["U"], TRANSFORM_FACTOR)
            sag = -U[sea]
            ratio = np.sum((sag * depth0[sea]) * np.broadcast_to(row_area[:, None], height.shape)[sea]) / \
                np.sum((depth0[sea]**2) * np.broadcast_to(row_area[:, None], height.shape)[sea])
            deep = {"Near-side Sea": nearside_deep, "South Pole-Aitken Sea": seeds[spa_name][1]}
            land = ~sea & (height < level0 + 3000.0)
            entry = dict(level_m=round(level, 2), level_change_m=round(level - level0, 2),
                         iterations=len(sol["history"]), seconds=round(time.time() - t0, 1),
                         volume_km3=round(vol / 1e9, 0), extra_volume_km3=round((vol - volume0) / 1e9, 0),
                         extra_global_layer_m=round((vol - volume0) / (4 * math.pi * MOON_RADIUS**2), 2),
                         extra_mass_kg=float((vol - volume0) * WATER_DENSITY),
                         **shore_measures(height, eff, level0, level, row_area, cell_km),
                         sag_m=dict(max=round(float(np.max(sag)), 1),
                                    depth_weighted_fraction=round(float(ratio), 4),
                                    at_deepest={k: dict(depth_m=round(float(depth0[v]), 0), sag_m=round(float(-U[v]), 1),
                                                        relative_m=round(float(-dev[v]), 1)) for k, v in deep.items()}),
                         land_uplift_max_m=round(float(np.max(U[land])), 2),
                         sills=strait_sills(eff, level, seeds), nearside=nearside_join(eff, level, mare_seeds))
            del U, sag, land, depth0, sea
            if target == "volume" and name.startswith(("elastic", "relaxed")):
                entry["lakes"], entry["lakes_volume_gain_km3"] = lake_changes(height, dev, upsample(sol["U"], TRANSFORM_FACTOR),
                                                                              cells, levels, spills, row_area)
                entry["lakes_volume_gain_km3"] = round(entry["lakes_volume_gain_km3"], 0)
            out[target] = entry
            saved[f"{name}|{target}"] = dict(U=sol["U"].astype(np.float32), N=sol["N"].astype(np.float32))
            del eff, dev
            print(f"{name:40s} {target:6s} level {level:9.2f}  sea {entry['sea_share']:.4f}  extra {entry['extra_volume_km3']:9.0f} km3  "
                  f"it {entry['iterations']:2d}  {entry['seconds']:.0f} s", flush=True)
        results[name] = out
    extra_loads = {}
    central = f"relaxed, Te {TE_CASES_KM[1]:g} km"
    if extra and central in sets:
        response = Response(*sets[central], filt)
        ref = results[central]["volume"]
        for label, load in (("groundwater, central, 10,000 years", groundwater_load(height, level0, cells)),
                            ("groundwater, central, 100,000 years", groundwater_load(height, level0, cells, "100000")),
                            ("air, 1.2 atm at the sea surface", air_load(height, level0)[0])):
            load = load - float(row_area @ load.sum(axis=1)) / (row_area.sum() * load.shape[1])
            sol = solve(height, row_area, response, lakes + load, level0, volume=volume0)
            eff = height + upsample(sol["dev"], TRANSFORM_FACTOR)
            shore = shore_measures(height, eff, level0, sol["level"], row_area, cell_km)
            sills = strait_sills(eff, sol["level"], seeds)
            extra_loads[label] = dict(
                load_range_kg_m2=[round(float(load.min()), 0), round(float(load.max()), 0)],
                sea_share_change=round(shore["sea_share"] - ref["sea_share"], 5),
                shore_rsl_change_m={k: round(shore["shore_rsl_change_m"][k] - ref["shore_rsl_change_m"][k], 2) for k in ("median", "p05", "p95")},
                sill_depth_change_m={k: round(sills[k]["sill_depth_m"] - ref["sills"][k]["sill_depth_m"], 2) for k in sills},
                nearside_margin_change_m=round(nearside_join(eff, sol["level"], mare_seeds)["margin_m"] - ref["nearside"]["margin_m"], 2))
            del eff
            print(f"extra load {label}: {extra_loads[label]}", flush=True)
        extra_loads["air_scale_height_km"] = round(air_load(height[:2], level0)[1] / 1e3, 1)
    np.savez_compressed(BULK / f"crustal_loading_nmax{nmax}.npz", lat_deg=90.0 - (np.arange(filt.nlat) + 0.5) / TRANSFORM_PPD,
                        lon_deg=(np.arange(filt.nlon) + 0.5) / TRANSFORM_PPD,
                        **{f"{k.replace(' ', '_').replace(',', '').replace('|', '__')}__{f}": v[f] for k, v in saved.items() for f in v})
    return dict(base=base, cases=results, love=notes, extra_loads=extra_loads, nmax=nmax,
                seconds=round(time.time() - started, 0), seeds={k: [list(map(int, a)), list(map(int, b))] for k, (a, b) in seeds.items()})


def producer():
    files = ["geography/love_numbers.py", "geography/crustal_loading.py", "geography/topography.py", "geography/basins.py"]
    return dict(domain="geography", files={f: sha256(ROOT / f) for f in files})


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--nmax", type=int, default=NMAX)
    parser.add_argument("--check", action="store_true", help="degree convergence at 179 and 719 for two cases")
    parser.add_argument("--out", type=Path, default=RESULT)
    args = parser.parse_args(argv)
    if fetch_inputs.restore()["status"] != "PASS":
        raise SystemExit("BLOCKED: topography inputs missing; run python -m geography.fetch_inputs --download")
    main_run = run(args.nmax)
    convergence = {}
    if args.check:
        keep = ("elastic", f"relaxed, Te {TE_CASES_KM[0]:g} km")
        for n in NMAX_CHECK:
            r = run(n, cases=keep, extra=False, targets=("volume",))
            convergence[str(n)] = {c: dict(level_change_m=r["cases"][c]["volume"]["level_change_m"],
                                           sea_share=r["cases"][c]["volume"]["sea_share"],
                                           sills={s: v["sill_depth_m"] for s, v in r["cases"][c]["volume"]["sills"].items()},
                                           sag_max_m=r["cases"][c]["volume"]["sag_m"]["max"]) for c in keep}
        convergence[str(args.nmax)] = {c: dict(level_change_m=main_run["cases"][c]["volume"]["level_change_m"],
                                               sea_share=main_run["cases"][c]["volume"]["sea_share"],
                                               sills={s: v["sill_depth_m"] for s, v in main_run["cases"][c]["volume"]["sills"].items()},
                                               sag_max_m=main_run["cases"][c]["volume"]["sag_m"]["max"]) for c in keep}
    write(main_run, convergence, args.out, args.nmax)


def write(main_run, convergence, out, nmax):
    inputs = [ATLAS_JSON, DRAINAGE_JSON, UPTAKE_JSON, ln.SOURCES, ROOT / "shared/constants.json"]
    ignored = [DRAINAGE_GRID] + [fetch_inputs.path(n) for n in ("ldem_16.img", "ldem_4.img", "jggrx_0420a_sha.tab")]
    meta = dict(
        schema=SCHEMA, producer=dict(**producer(), inputs={str(p.relative_to(ROOT)): sha256(p) for p in inputs},
                                     ignored_inputs={str(Path(p).relative_to(ROOT)): sha256(p) for p in ignored}),
        evidence=("A calculation on committed products and published values: load Love numbers of a published lunar interior "
                  "(Briaud et al. 2023's tidal fit: density and rigidity by layer, fluid core), elastic by integration of the static "
                  "self-gravitating equations, relaxed as a thin elastic shell over an interior in hydrostatic equilibrium with the "
                  "published effective elastic thicknesses of the mare basins; the sea-level equation with migrating shorelines solved "
                  "on LOLA at 16 px/deg with spherical harmonics to degree %d. Checked against closed forms (homogeneous sphere, "
                  "Boussinesq limit, fluid limit) and against the Moon's measured k2 and h2. Not validated against any observed "
                  "loading on the Moon, which has none of this size." % nmax),
        reading_rule=("Read 'elastic' as the response within years of the water's arrival and 'relaxed, Te' as the state the "
                      "mantle's flow approaches over the relaxation times given (thousands to hundreds of thousands of years). "
                      "'volume' holds today's water (the atlas's 28% inventory) and lets the shore move; 'share' holds 28% cover "
                      "and gives the water that takes. Sill depths are below the sea surface, four-connected cells at 16 px/deg. "
                      "Sag is the ground's radial fall (degrees 2 and up, and degree 1 in the solid Moon's frame); 'relative' is the "
                      "fall below the sea's equipotential. The products of geography are not regenerated: whether to move the "
                      "shorelines is a decision this informs."),
        settings=dict(nmax=nmax, transform_ppd=TRANSFORM_PPD, fine_ppd=FINE_PPD, water_density_kg_m3=WATER_DENSITY,
                      te_cases_km=list(TE_CASES_KM), thick_te_sensitivity_km=THICK_TE_KM, young_pa=YOUNG_PA, poisson=POISSON,
                      mean_density_kg_m3=round(MEAN_DENSITY, 2), tolerance_m=TOLERANCE_M,
                      air=dict(pressure_atm=AIR_PRESSURE_ATM, temperature_k=AIR_TEMPERATURE_K, molar_mass_kg_mol=AIR_MOLAR_MASS)),
        relaxation_times=relaxation_times(), convergence=convergence, **main_run)
    out.write_text(json.dumps(meta, indent=1, default=lambda x: x.item() if hasattr(x, "item") else str(x)) + "\n")
    print(f"wrote {out} in {main_run['seconds']:.0f} s")


if __name__ == "__main__":
    main()
