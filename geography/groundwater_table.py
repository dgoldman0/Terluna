"""The water table on land: where groundwater stands, where it comes out, and what the seas and lakes trade with it.

    python -m geography.groundwater_table           # results/groundwater_table.json, about two minutes, 2 GB

Groundwater's second stage. geography.groundwater_uptake follows how fast the cold, porous crust takes the delivered
water down columns; this module follows the water sideways. On the atlas averaged to 2 pixels a degree (15 km), about
the resolution of Earth's global groundwater models, it solves steady Dupuit-Forchheimer flow in the upper crust: heads
on land, the seas and the rain-fed lakes held at their levels, recharge from the design climate, and the crust's
transmissivity from stage 1's permeability, falling with depth below the ground so that a deep water table carries
little. Where groundwater would rise above the ground it comes out, as springs, seeps and wetlands do (a drain at the
surface, as global groundwater models treat seepage; de Graaf et al. 2015, Fan et al. 2013). Pseudo-transient Newton
iteration solves it: the 4 pixel/degree grid's steeper relief from cell to cell kept it from converging there, and
finer regional solutions belong with the atlas round's districts.

It solves three states for each of stage 1's cases: the mature Moon, once the deep crust has filled; and the young Moon
at 1,000 and 10,000 years, with stage 1's deep uptake under the land taken from the recharge. The results are the depth
to water, how much of the land is waterlogged or seeping, how much groundwater reaches the surface, the seas and the
lakes, and which lakes lose water to the ground and how much.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu

from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_RADIUS, MOON_SURFACE_GRAVITY

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = HERE / 'results' / 'groundwater_table.json'
GRIDS = ROOT / 'research' / 'runs' / 'resources' / 'groundwater' / 'groundwater_table_2ppd.npz'
SCHEMA = 'terluna.geography.groundwater-table/1'
INPUTS = ('geography/results/atlas.json', 'geography/results/groundwater_uptake.json')
IGNORED = ('climate/gcm/products/climatology_A28_dim5_moon.npz', 'geography/products/atlas_28pct_4ppd.npz',
           'geography/products/drainage_28pct_16ppd.npz')
CODE = ('geography/groundwater_table.py', 'geography/groundwater_uptake.py')
YEAR_S = JULIAN_DAY * JULIAN_YEAR_DAYS
PPD = 2                                  # the solved grid: the atlas's 4 pixels a degree averaged two by two
RHO_WATER = 1000.0
WATER_K = 290.0                          # the groundwater's temperature for its viscosity, near the floors'
DRAIN_FACTOR = 100.0                     # a seeping cell's drain conductance against its largest face conductance
DEPTH_CLASSES_M = (0.0, 2.0, 10.0, 50.0, 200.0, 1000.0, math.inf)
STATES = ('mature', 'young_10000_yr', 'young_1000_yr')   # the mature state first; the young ones start from it
SPANS = dict(young_1000_yr='100-1000', young_10000_yr='1000-10000')
# Newton's method to steady state.
SEEP_SMOOTHING_M = 0.5                   # the drain's softplus width at the ground
NEWTON_TOLERANCE = 1e-4                  # the cells' root-mean-square imbalance, against their recharge or the mean
NEWTON_WORST = 0.05                      # and the largest any one cell may keep
MAX_NEWTON = 200
MAX_HEAD_STEP_M = 500.0
SPECIFIC_YIELD = 0.1                     # the pseudo-time storage that damps the early steps (it sets the path only)
PTC_FIRST_YR = 1.0
PTC_GROWTH = 4.0
PTC_LAST_YR = 1e12
PTC_REJECT = 2.0
VERBOSE = False
EVIDENCE = ('Steady Dupuit-Forchheimer flow on the atlas averaged to 2 pixels a degree, with the seas and the drainage '
            'product\'s lakes as fixed heads, recharge from the design climate run on stage 1\'s terms, transmissivity '
            'from stage 1\'s permeability integrated below the water table at lunar gravity and 290 K, and seepage '
            'where the water table reaches the ground. The young states take stage 1\'s deep uptake under the land, '
            'band by band, from the recharge. Solved by pseudo-transient Newton iteration to a root-mean-square '
            'imbalance of 1e-4 of the recharge. Relief within a 15 km cell, evaporation from shallow water tables, '
            'flow through the deep crust and the lakes\' and rivers\' own beds are outside it.')
READING_RULE = ('Depth to water in m below the cell\'s mean ground, 0 where groundwater seeps out. Land shares are of '
                'the land outside the seas and the 2 pixel/degree lakes. Flows in km3 a year: recharge, seepage at '
                'the surface, groundwater into the seas and into the lakes (negative where the seas or lakes lose '
                'water to the ground). Lakes are 8-connected 2 pixel/degree lake cells; their exchange is m a year '
                'over their area, positive where groundwater flows in.')


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def log(msg, t0=[time.time()]):
    print(f'{time.time() - t0[0]:7.1f} s  {msg}', file=sys.stderr, flush=True)


def viscosity(t_k):
    return 2.414e-5 * 10.0 ** (247.8 / (t_k - 140.0))


def transmissivity(depth, k_1km, cap, exponent):
    """m2/s below a water table depth m down: the permeability holds its cap above the depth where the falling law
    meets it, and falls as depth^-exponent below; integrated at lunar gravity and the water's viscosity."""
    k_cap = cap * RHO_WATER * MOON_SURFACE_GRAVITY / viscosity(WATER_K)
    z_c = 1000.0 * (k_1km / cap) ** (1.0 / exponent)
    d = np.maximum(depth, 0.0)
    above = k_cap * (z_c - d) + k_cap * z_c / (exponent - 1.0)
    below = k_cap * z_c * (np.maximum(d, z_c) / z_c) ** (1.0 - exponent) / (exponent - 1.0)
    return np.where(d <= z_c, above, below)


# --- Grid and inputs ------------------------------------------------------------------------------------------------

def bilinear(field, lat_src, lon_src, lat, lon):
    """A T21 field on the atlas's grid (rows from the north, longitudes east)."""
    order = np.argsort(lat_src)
    ls, fs = lat_src[order], field[order]
    lon_ext = np.concatenate([lon_src - 360.0, lon_src, lon_src + 360.0])
    fs = np.concatenate([fs, fs, fs], axis=1)
    i = np.clip(np.searchsorted(ls, np.clip(lat, ls[0], ls[-1])) - 1, 0, ls.size - 2)
    ty = ((np.clip(lat, ls[0], ls[-1]) - ls[i]) / (ls[i + 1] - ls[i]))[:, None]
    j = np.clip(np.searchsorted(lon_ext, lon) - 1, 0, lon_ext.size - 2)
    tx = ((lon - lon_ext[j]) / (lon_ext[j + 1] - lon_ext[j]))[None, :]
    top = fs[i][:, j] * (1 - tx) + fs[i][:, j + 1] * tx
    bottom = fs[i + 1][:, j] * (1 - tx) + fs[i + 1][:, j + 1] * tx
    return top * (1 - ty) + bottom * ty


def setting():
    atlas = json.loads((ROOT / INPUTS[0]).read_text())
    z = np.load(ROOT / IGNORED[1])
    lat, lon = z['lat_deg'].astype(float), z['lon_deg'].astype(float)
    if not lat[0] > lat[-1]:
        raise ValueError('Atlas rows must run from the north')
    rows, cols = lat.size, lon.size
    ground = z['height_m'].astype(float)
    sea = z['water_label'] > 0
    level = float(atlas['sea_level_m'])
    # The drainage product's 16 pixel/degree lakes on this grid: a cell is lake where most of its sub-cells are.
    dr = np.load(ROOT / IGNORED[2])
    cells = dr['lake_cell'].astype(np.int64)
    r16, c16 = np.divmod(cells, 360 * 16)
    r4, c4 = r16 // 4, c16 // 4
    flat = r4 * cols + c4
    count = np.bincount(flat, minlength=rows * cols).reshape(rows, cols)
    level_sum = np.bincount(flat, weights=dr['lake_level_m'].astype(float), minlength=rows * cols).reshape(rows, cols)
    lake = (count >= 8) & ~sea
    lake_level = np.where(lake, level_sum / np.maximum(count, 1), np.nan)
    clim = dict(np.load(ROOT / IGNORED[0]))
    clim.pop('metadata')
    rain = np.maximum(bilinear(clim['pr_mm_day'], clim['lat'], clim['lon'], lat, lon), 0.0) * 365.25 / 1e3
    evap = np.maximum(bilinear(clim['evap_mm_day'], clim['lat'], clim['lon'], lat, lon), 0.0) * 365.25 / 1e3
    return dict(lat=lat, lon=lon, rows=rows, cols=cols, ground=ground, sea=sea, sea_level=level, lake=lake,
                lake_level=lake_level, rain=rain, evap=evap)


def recharge(s, share, arid_share):
    """m/s on each cell: a share of the rain, at most the rain less evaporation, at least the arid share."""
    r = np.minimum(share * s['rain'], np.maximum(s['rain'] - s['evap'], 0.0))
    return np.maximum(r, arid_share * s['rain']) / YEAR_S


# --- The solver -----------------------------------------------------------------------------------------------------

def geometry(s):
    rows, cols = s['rows'], s['cols']
    step = abs(s['lat'][0] - s['lat'][1])
    dphi = math.radians(step)
    dlam = math.radians(360.0 / cols)
    phi = np.radians(s['lat'])
    dy = MOON_RADIUS * dphi
    dx = MOON_RADIUS * np.cos(phi) * dlam                     # cell width at its centre
    phi_face = np.radians(s['lat'][:-1] - 0.5 * step)         # the north-south faces between rows
    dx_face = MOON_RADIUS * np.cos(phi_face) * dlam
    band = MOON_RADIUS ** 2 * dlam * (np.sin(phi + 0.5 * dphi) - np.sin(phi - 0.5 * dphi))
    area = band[:, None] * np.ones((1, cols))
    return dx, dy, dx_face, area


def conductivity_at(depth, k_1km, cap, exponent):
    """m/s at a depth below the ground: the derivative of the transmissivity with the water table's height."""
    k_cap = cap * RHO_WATER * MOON_SURFACE_GRAVITY / viscosity(WATER_K)
    z_c = 1000.0 * (k_1km / cap) ** (1.0 / exponent)
    d = np.maximum(depth, 0.0)
    return np.where(depth <= 0.0, 0.0, np.where(d <= z_c, k_cap, k_cap * (np.maximum(d, z_c) / z_c) ** -exponent))


def solve(s, rech, k_1km, cap, exponent, h0=None, seeping0=None):
    """Steady heads on land with fixed seas and lakes and seepage at the ground, by Newton's method on the cells'
    balances: recharge, the flows to their neighbours through the transmissivity below the water table, and a drain
    that lets water out where the water table reaches the ground (smoothed over SEEP_SMOOTHING_M so the balance is
    differentiable). The transmissivity's derivative with the water table's height is the conductivity there. Each
    step is backtracked until the balances' norm falls."""
    rows, cols = s['rows'], s['cols']
    n = rows * cols
    dx, dy, dx_face, area = geometry(s)
    fixed = (s['sea'] | s['lake']).ravel()
    fixed_h = np.where(s['sea'], s['sea_level'], s['lake_level']).ravel()
    ground = s['ground'].ravel()
    free = np.flatnonzero(~fixed)
    pos = np.full(n, -1, np.int64)
    pos[free] = np.arange(free.size)
    idx = np.arange(n).reshape(rows, cols)
    a_all = np.concatenate([idx.ravel(), idx[:-1].ravel()])
    b_all = np.concatenate([np.roll(idx, -1, axis=1).ravel(), idx[1:].ravel()])
    g_all = np.concatenate([np.repeat(dy / dx, cols), np.repeat(dx_face / dy, cols)])
    m = free.size
    q = rech.ravel() * area.ravel()                            # m3/s of recharge on each cell
    q[fixed] = 0.0
    t0 = transmissivity(np.zeros(n), k_1km, cap, exponent)
    c0 = t0[a_all] * g_all
    cdrain = np.zeros(n)
    np.maximum.at(cdrain, a_all, c0)
    np.maximum.at(cdrain, b_all, c0)
    cdrain *= DRAIN_FACTOR
    delta = SEEP_SMOOTHING_M
    scale = np.maximum(q, np.mean(q[~fixed]))                  # what a cell's imbalance is judged against

    def balance(h):
        d = ground - h
        t = transmissivity(d, k_1km, cap, exponent)
        ta, tb = t[a_all], t[b_all]
        c = 2.0 * ta * tb / (ta + tb) * g_all
        f = c * (h[a_all] - h[b_all])                          # m3/s from a to b
        x = (h - ground) / delta
        drain = cdrain * delta * np.logaddexp(0.0, x)          # softplus: the seepage
        r = q + np.bincount(b_all, weights=f, minlength=n) - np.bincount(a_all, weights=f, minlength=n) - drain
        r[fixed] = 0.0
        return r, d, t, c, f, drain, x

    h = np.where(fixed, fixed_h, ground) if h0 is None else np.where(fixed, fixed_h, h0.ravel())
    r, d, t, c, f, drain, x = balance(h)
    norm = float(np.sqrt(np.mean((r[free] / scale[free]) ** 2)))
    storage = SPECIFIC_YIELD * area.ravel() / YEAR_S            # m2 per year of pseudo-time
    dt, steps, worst = PTC_FIRST_YR, 0, math.inf
    for steps in range(MAX_NEWTON):
        worst = float(np.max(np.abs(r[free]) / scale[free]))
        if VERBOSE:
            log(f'  newton {steps}: dt {dt:.3g} yr, rms imbalance {norm:.3g}, worst {worst:.3g}, seeping '
                f'{int(((h > ground) & ~fixed).sum())}')
        if norm < NEWTON_TOLERANCE and worst < NEWTON_WORST:
            break
        kc = conductivity_at(d, k_1km, cap, exponent)          # dT/dh
        ta, tb = t[a_all], t[b_all]
        dc_da = 2.0 * tb * tb / (ta + tb) ** 2 * g_all * kc[a_all]
        dc_db = 2.0 * ta * ta / (ta + tb) ** 2 * g_all * kc[b_all]
        dh = h[a_all] - h[b_all]
        df_da = c + dh * dc_da                                  # d(flow a to b)/dh_a
        df_db = -c + dh * dc_db
        ddrain = cdrain / (1.0 + np.exp(-np.clip(x, -50.0, 50.0)))
        # Pseudo-transient Newton: storage over dt less the Jacobian of r (r_a -= f, r_b += f, r -= drain).
        rows_j = np.concatenate([a_all, a_all, b_all, b_all, np.arange(n)])
        cols_j = np.concatenate([a_all, b_all, a_all, b_all, np.arange(n)])
        vals = np.concatenate([df_da, df_db, -df_da, -df_db, ddrain + storage / dt])
        keep = ~fixed[rows_j] & ~fixed[cols_j]
        M = sp.csc_matrix((vals[keep], (pos[rows_j[keep]], pos[cols_j[keep]])), shape=(m, m))
        step = np.zeros(n)
        step[free] = splu(M, permc_spec='COLAMD').solve(r[free])
        step = np.clip(step, -MAX_HEAD_STEP_M, MAX_HEAD_STEP_M)
        trial = h + step
        r2, d2, t2, c2, f2, drain2, x2 = balance(trial)
        norm2 = float(np.sqrt(np.mean((r2[free] / scale[free]) ** 2)))
        if norm2 > PTC_REJECT * norm:
            dt *= 0.25                                          # the step made things worse: shorten it and retry
            continue
        if norm2 <= norm:
            dt = min(dt * PTC_GROWTH, PTC_LAST_YR)              # damping no longer needed as much: lengthen the step
        h, r, d, t, c, f, drain, x, norm = trial, r2, d2, t2, c2, f2, drain2, x2, norm2
    into = np.bincount(b_all, weights=f, minlength=n) - np.bincount(a_all, weights=f, minlength=n)
    seep = np.where(~fixed, drain, 0.0)
    return dict(h=h.reshape(rows, cols), seep=seep.reshape(rows, cols),
                seeping=((h > ground - delta) & ~fixed).reshape(rows, cols),
                into_fixed=np.where(fixed, into, 0.0).reshape(rows, cols), iterations=steps + 1, worst=worst,
                rms=norm, recharge=q.reshape(rows, cols))


def coarsen(s, f=2):
    """The setting on a grid f times coarser: means of the ground and climate, sea and lakes by majority."""
    rows, cols = s['rows'] // f, s['cols'] // f
    blk = lambda a: a.reshape(rows, f, cols, f)
    out = dict(rows=rows, cols=cols, sea_level=s['sea_level'], lat=s['lat'].reshape(rows, f).mean(1),
               lon=s['lon'].reshape(cols, f).mean(1), ground=blk(s['ground']).mean((1, 3)),
               rain=blk(s['rain']).mean((1, 3)), evap=blk(s['evap']).mean((1, 3)))
    out['sea'] = blk(s['sea']).sum((1, 3)) * 2 >= f * f
    lake_n = blk(s['lake']).sum((1, 3))
    out['lake'] = (lake_n * 2 >= f * f) & ~out['sea']
    level = blk(np.where(s['lake'], s['lake_level'], 0.0)).sum((1, 3)) / np.maximum(lake_n, 1)
    out['lake_level'] = np.where(out['lake'], level, np.nan)
    return out


# --- The states -----------------------------------------------------------------------------------------------------

def band_of(lat):
    from geography.groundwater_uptake import BANDS_DEG
    a = np.abs(lat)
    out = np.zeros(lat.shape, np.int64)
    for i, (lo, hi) in enumerate(BANDS_DEG):
        out[(a >= lo) & ((a < hi) if hi < 90.0 else (a <= hi))] = i
    return out


def lakes_of(s):
    from scipy import ndimage
    lab, n = ndimage.label(s['lake'], structure=np.ones((3, 3), bool))
    return lab, n


def summarise(s, r, state_rech, lake_lab, n_lakes):
    from geography.groundwater_uptake import BANDS_DEG
    dx, dy, dx_face, area = geometry(s)
    land = ~(s['sea'] | s['lake'])
    depth = np.where(land, np.maximum(s['ground'] - r['h'], 0.0), np.nan)
    depth = np.where(r['seeping'], 0.0, depth)
    a_land = area[land].sum()
    lat2 = s['lat'][:, None] * np.ones((1, s['cols']))
    bands = band_of(lat2)

    def shares(mask):
        a = area[mask]
        d = depth[mask]
        tot = a.sum()
        wet = r['seeping'][mask]
        out = {'seeping': float(a[wet].sum() / tot)}
        for lo, hi in zip(DEPTH_CLASSES_M[:-1], DEPTH_CLASSES_M[1:]):
            name = f'{lo:g}_to_{hi:g}_m' if math.isfinite(hi) else f'over_{lo:g}_m'
            out[name] = float(a[(d >= lo) & (d < hi) & ~wet].sum() / tot)
        return out
    km3 = lambda x: float(x) * YEAR_S / 1e9
    into = r['into_fixed']
    exchange = np.bincount(lake_lab.ravel(), weights=into.ravel(), minlength=n_lakes + 1)
    lake_area = np.bincount(lake_lab.ravel(), weights=(area * s['lake']).ravel(), minlength=n_lakes + 1)
    rate = exchange[1:] / np.maximum(lake_area[1:], 1.0) * YEAR_S        # m/yr over the lake, + gaining
    la = lake_area[1:] / 1e6
    by_band = {}
    for i, (lo, hi) in enumerate(BANDS_DEG):
        m = land & (bands == i)
        if not m.any():
            continue
        by_band[f'{lo:g}-{hi:g}'] = dict(
            land_km2=float(area[m].sum() / 1e6), recharge_mm_yr=float((state_rech * area)[m].sum() / area[m].sum()
                                                                     * YEAR_S * 1e3),
            seepage_km3_yr=km3((r['seep'] * 1.0)[m].sum()), depth_shares=shares(m),
            median_depth_m=float(np.median(depth[m])))
    losing = rate < -0.001
    gaining = rate > 0.001
    return dict(
        iterations=r['iterations'], rms_imbalance=r['rms'], worst_imbalance=r['worst'],
        land_km2=float(a_land / 1e6),
        recharge_km3_yr=km3((state_rech * area)[land].sum()),
        seepage_km3_yr=km3(r['seep'][land].sum()),
        into_seas_km3_yr=km3(into[s['sea']].sum()),
        into_lakes_km3_yr=km3(into[s['lake']].sum()),
        depth_shares=shares(land),
        median_depth_m=float(np.nanmedian(depth[land])),
        by_band=by_band,
        lakes=dict(count=int(n_lakes), losing=int(losing.sum()), gaining=int(gaining.sum()),
                   losing_area_km2=float(la[losing].sum()), gaining_area_km2=float(la[gaining].sum()),
                   loss_m_yr_percentiles={f'p{p}': float(np.percentile(-rate[losing], p)) if losing.any() else 0.0
                                          for p in (50, 90, 99)},
                   gain_m_yr_percentiles={f'p{p}': float(np.percentile(rate[gaining], p)) if gaining.any() else 0.0
                                          for p in (50, 90, 99)},
                   losing_faster_than={f'{x:g}_m_yr': dict(count=int((rate < -x).sum()),
                                                            area_km2=float(la[rate < -x].sum()))
                                       for x in (0.1, 1.0, 10.0)}),
        seas_losing_to_land_km2=float((area * (s['sea'] & (into < 0))).sum() / 1e6),
        seas_losing_to_land_km3_yr=km3(-into[s['sea'] & (into < 0)].sum()))


def results():
    from geography import groundwater_uptake as gu
    stage1 = json.loads((ROOT / INPUTS[1]).read_text())
    s = coarsen(setting())
    log('setting')
    lake_lab, n_lakes = lakes_of(s)
    lat2 = s['lat'][:, None] * np.ones((1, s['cols']))
    bands = band_of(lat2)
    band_names = [f'{lo:g}-{hi:g}' for lo, hi in gu.BANDS_DEG]
    out, grids = {}, {}
    for case in gu.CASES:
        p = stage1['cases'][case]['params']
        full = recharge(s, p['recharge_share'], gu.ARID_RECHARGE_SHARE)
        rows = {}
        h0 = None
        for state in STATES:
            if state == 'mature':
                rech = full
            else:
                frac = np.zeros(len(band_names))
                for i, name in enumerate(band_names):
                    col = stage1['cases'][case]['bands'][name]['columns']['land']
                    year = SPANS[state].split('-')[1]
                    frac[i] = min(col['by_year'][year]['uptake_m_yr_since_last'] / max(col['recharge_mm_yr'] / 1e3,
                                                                                      1e-12), 1.0)
                rech = full * (1.0 - frac[bands])
            r = solve(s, rech, p['perm_1km'], gu.PERM_CAP_M2, gu.PERM_EXPONENT, h0)
            log(f'{case} {state}: {r["iterations"]} iterations')
            rows[state] = summarise(s, r, rech, lake_lab, n_lakes)
            if state == 'mature':
                h0 = r['h']
            if case == 'central':
                land = ~(s['sea'] | s['lake'])
                grids[f'depth_m_{state}'] = np.where(land, np.maximum(s['ground'] - r['h'], 0.0), np.nan).astype(
                    np.float32)
                grids[f'seepage_mm_yr_{state}'] = (r['seep'] / geometry(s)[3] * YEAR_S * 1e3).astype(np.float32)
        out[case] = dict(perm_1km_m2=p['perm_1km'], recharge_share=p['recharge_share'], states=rows)
    grids.update(lake_label=lake_lab.astype(np.int32), lat_deg=s['lat'].astype(np.float32),
                 lon_deg=s['lon'].astype(np.float32))
    return dict(
        schema=SCHEMA, evidence=EVIDENCE, reading_rule=READING_RULE,
        setting=dict(grid_pixels_per_degree=PPD, lakes_on_grid=int(n_lakes),
                     lake_km2_on_grid=float((geometry(s)[3] * s['lake']).sum() / 1e6),
                     sea_km2_on_grid=float((geometry(s)[3] * s['sea']).sum() / 1e6), water_k=WATER_K,
                     drain_factor=DRAIN_FACTOR, transmissivity_m2_s_at_depth={
                         case: {f'{d:g} m': float(transmissivity(np.array(d), out[case]['perm_1km_m2'],
                                                                 gu.PERM_CAP_M2, gu.PERM_EXPONENT))
                                for d in (0.0, 100.0, 500.0, 2000.0)} for case in gu.CASES}),
        cases=out), grids


def main(argv=None) -> int:
    from geography.groundwater_uptake import rnd
    product, grids = results()
    product = rnd(product)
    product['producer'] = dict(
        domain='geography', files={p: sha256(ROOT / p) for p in CODE},
        inputs={p: sha256(ROOT / p) for p in INPUTS}, ignored_inputs={p: sha256(ROOT / p) for p in IGNORED})
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1, allow_nan=False) + '\n')
    GRIDS.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(GRIDS, metadata=json.dumps(dict(schema=SCHEMA + '-grids', reading_rule=READING_RULE)), **grids)
    for case, v in product['cases'].items():
        for state, row in v['states'].items():
            print(case, state, 'seeping', row['depth_shares']['seeping'], 'median depth', row['median_depth_m'])
    return 0


if __name__ == '__main__':
    sys.exit(main())
