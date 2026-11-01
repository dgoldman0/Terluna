"""Junction loads in a snowflake aggregate falling at terminal speed: a first contact-network calculation toward a
breaking law for lunar snow (the README's step 6, a breakup law in place of the load analogy).

    OPENBLAS_NUM_THREADS=1 python -m atmosphere.snow.breakup_network      # writes results/breakup_network.json

Monomers are thin circular plates of one diameter d0 and one mass, each with a random orientation. Hierarchical
ballistic cluster-cluster aggregation in 3D joins pairs of equal clusters, each turned at random, by moving one
straight toward the other with a random offset until they touch. The plates' circumscribing spheres (diameter d0)
decide the contact, so each merge makes exactly one contact and the contacts form a spanning tree.

A flake falls at a fixed orientation (its thinnest principal axis vertical; a random orientation is the sensitivity)
at terminal speed, so its whole drag equals its weight. The drag is shared over the plates by their projected areas
seen from below, on a grid of 0.1 d0 pixels: each pixel gives its area to the lowest plate over it (an opaque shadow,
the most uneven sharing), to the plates over it in shares (1-s)^rank (a partial shadow, s = 0.5), or to every plate
over it (no shadow: each plate's own projected area). The force and torque that weight and drag leave unbalanced are
carried by the flake's linear and angular acceleration (d'Alembert); every tree contact then carries the force and
the moment, about the contact point, of everything on one side of it. Turbulence is a vertical air velocity varying
linearly across the flake, by du over its span, which changes each plate's drag by 2 du_local / V.

The junction strength is fixed by Earth flakes of 4 and 6 cm at their limit; the lunar span with the same peak force,
or moment, follows from the fitted size dependence. Fall speed, mass laws and the turbulent increment come from the
package (turbulence_survival.load, Heymsfield and Westbrook 2010). This is a calculation on an idealised geometry and
drag rule; nothing in it is measured.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import resource
import sys
import time
from pathlib import Path

import numpy as np

from atmosphere.snow import flakes
from atmosphere.snow import turbulence_survival as ts

HERE = Path(__file__).resolve().parent
PRODUCT = HERE / 'results/breakup_network.json'
SCHEMA = 'terluna.atmosphere.snow-breakup-network/1'
SEED = 20261030
D0_M = 1.0e-3                    # m: the monomer (plate) diameter, the geometry's length unit
PIXEL = 0.1                      # z-buffer pixel, in plate diameters
OWN_INERTIA = 1.0 / 12.0         # a plate's own moment of inertia, m0 d0^2, averaged over axes (thin disk: 1/16-1/8)
SHADOWS = {'shadow': 1.0, 'partial': 0.5, 'own': 0.0}
ORIENTATIONS = ('broadside', 'random')
X_GRID = np.array([0.0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25, 0.3])   # du/V
FIT_MIN_N = (64, 256)
EARTH_SPANS_M = (0.04, 0.06)
DISSIPATION = (1e-4, 1e-3, 1e-2)
BOOT = 400


# ---------------------------------------------------------------- aggregation

def rotation(rng):
    q = rng.normal(size=4)
    w, x, y, z = q / np.linalg.norm(q)
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                     [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                     [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])


def close_pairs(pa, pb):
    """Pairs (i of A, j of B) whose centres lie within one diameter of each other in the horizontal plane."""
    na, nb = len(pa), len(pb)
    if na * nb <= 4096:
        i, j = np.divmod(np.arange(na * nb), nb)
    else:                                                    # a 2D grid of cells one diameter wide
        ca = np.floor(pa[:, :2]).astype(np.int64)
        cb = np.floor(pb[:, :2]).astype(np.int64)
        low = np.minimum(ca.min(0), cb.min(0)) - 2
        ca, cb = ca - low, cb - low
        width = int(max(ca[:, 1].max(), cb[:, 1].max())) + 3
        ka = ca[:, 0] * width + ca[:, 1]
        order = np.argsort(ka, kind='stable')
        ks = ka[order]
        ii, jj = [], []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                kb = (cb[:, 0] + dx) * width + cb[:, 1] + dy
                lo = np.searchsorted(ks, kb, 'left')
                cnt = np.searchsorted(ks, kb, 'right') - lo
                tot = int(cnt.sum())
                if tot:
                    first = np.repeat(lo - (np.cumsum(cnt) - cnt), cnt)
                    ii.append(order[first + np.arange(tot)])
                    jj.append(np.repeat(np.arange(nb), cnt))
        if not ii:
            return np.empty(0, np.int64), np.empty(0, np.int64), np.empty(0)
        i, j = np.concatenate(ii), np.concatenate(jj)
    d2 = (pa[i, 0] - pb[j, 0]) ** 2 + (pa[i, 1] - pb[j, 1]) ** 2
    keep = d2 < 1.0
    return i[keep], j[keep], d2[keep]


def merge(a, b, rng):
    """Ballistic merge: both clusters turned at random, B moved along +z from far below with a random offset until
    the first pair of plates' circumscribing spheres touch. Positions are kept about the centre of mass."""
    ra, rb = rotation(rng), rotation(rng)
    pa, va = a['p'] @ ra.T, a['n'] @ ra.T
    pb0, vb = b['p'] @ rb.T, b['n'] @ rb.T
    reach = math.sqrt((pa * pa).sum(1).max()) + math.sqrt((pb0 * pb0).sum(1).max()) + 1.0
    for _ in range(100000):
        rr, ph = reach * math.sqrt(rng.random()), 2.0 * math.pi * rng.random()
        pb = pb0 + np.array([rr * math.cos(ph), rr * math.sin(ph), 0.0])
        i, j, d2 = close_pairs(pa, pb)
        if i.size:
            t = pa[i, 2] - pb[j, 2] - np.sqrt(1.0 - d2)
            k = int(np.argmin(t))
            pb[:, 2] += t[k]
            p = np.vstack([pa, pb])
            na = len(pa)
            return dict(p=p - p.mean(0), n=np.vstack([va, vb]),
                        e=np.vstack([a['e'], b['e'] + na, [[i[k], na + j[k]]]]))
    raise RuntimeError('no contact found')


def monomers(count, rng):
    v = rng.normal(size=(count, 3))
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    return [dict(p=np.zeros((1, 3)), n=v[k:k + 1].copy(), e=np.zeros((0, 2), np.int64)) for k in range(count)]


# ---------------------------------------------------------------- loads

def tree(e, n):
    """Parent of every node (root 0) and the breadth-first levels of the contact tree."""
    parent = np.full(n, -1, np.int64)
    levels = [np.array([0])]
    if n == 1:
        return parent, levels
    src, dst = np.concatenate([e[:, 0], e[:, 1]]), np.concatenate([e[:, 1], e[:, 0]])
    o = np.argsort(src, kind='stable')
    src, dst = src[o], dst[o]
    start = np.searchsorted(src, np.arange(n + 1))
    seen = np.zeros(n, bool)
    seen[0] = True
    front = levels[0]
    while True:
        cnt = start[front + 1] - start[front]
        tot = int(cnt.sum())
        if not tot:
            break
        nb = dst[np.repeat(start[front] - (np.cumsum(cnt) - cnt), cnt) + np.arange(tot)]
        par = np.repeat(front, cnt)
        new = ~seen[nb]
        nb, par = nb[new], par[new]
        if not nb.size:
            break
        seen[nb] = True
        parent[nb] = par
        levels.append(nb)
        front = nb
    if not seen.all():
        raise RuntimeError('contact graph is not connected')
    return parent, levels


def subtree_sum(v, parent, levels):
    s = v.copy()
    for lv in reversed(levels[1:]):
        np.add.at(s, parent[lv], s[lv])
    return s


def span(p):
    """Maximum dimension: the largest centre-to-centre distance plus one plate diameter."""
    sq = (p * p).sum(1)
    best = 0.0
    for s in range(0, len(p), 512):
        best = max(best, float((sq[s:s + 512, None] + sq[None, :] - 2.0 * (p[s:s + 512] @ p.T)).max()))
    return math.sqrt(max(best, 0.0)) + 1.0


def orient(p, v, how, rng):
    if how == 'broadside':                     # thinnest principal axis vertical, the broad side facing the flow
        rot = np.linalg.eigh(p.T @ p)[1][:, ::-1].copy()
        if np.linalg.det(rot) < 0:
            rot[:, 0] = -rot[:, 0]
    else:
        rot = rotation(rng)
    return p @ rot, v @ rot


def projected_areas(p, v):
    """Area each plate presents to the flow from below under each shadow rule, and the union area."""
    n, r, h = len(p), 0.5, PIXEL
    m = int(math.ceil(r / h))
    off = np.arange(-m, m + 1)
    ox, oy = (a.ravel() for a in np.meshgrid(off, off, indexing='ij'))
    ix = np.floor(p[:, 0] / h).astype(np.int64)[:, None] + ox
    iy = np.floor(p[:, 1] / h).astype(np.int64)[:, None] + oy
    ux, uy = (ix + 0.5) * h - p[:, :1], (iy + 0.5) * h - p[:, 1:2]
    nz = v[:, 2:3]
    nz = np.where(np.abs(nz) < 1e-6, 1e-6, nz)
    uz = -(v[:, :1] * ux + v[:, 1:2] * uy) / nz
    inside = ux * ux + uy * uy + uz * uz <= r * r
    owner = np.nonzero(inside)[0]
    z = (p[:, 2:3] + uz)[inside]
    ix, iy = ix[inside], iy[inside]
    key = (ix - ix.min()) * (int(iy.max() - iy.min()) + 1) + (iy - iy.min())
    o = np.lexsort((z, key))
    key, owner = key[o], owner[o]
    first = np.ones(key.size, bool)
    first[1:] = key[1:] != key[:-1]
    idx = np.arange(key.size)
    rank = idx - np.maximum.accumulate(np.where(first, idx, 0))
    areas = {}
    for name, s in SHADOWS.items():
        if s == 0.0:
            areas[name] = math.pi * r * r * np.abs(v[:, 2])
        else:
            areas[name] = np.bincount(owner, weights=(1.0 - s) ** rank, minlength=n) * h * h
    return areas, float(first.sum()) * h * h


def equilibrate(p, f, inertia_inv):
    """Add the inertial forces and couples of the flake's linear and angular acceleration (unit plate masses)."""
    n = len(p)
    alpha = inertia_inv @ np.cross(p, f).sum(0)
    return f - f.sum(0) / n - np.cross(alpha, p), np.broadcast_to(-OWN_INERTIA * alpha, (n, 3))


def junction_loads(p, parent, levels, kids, drag, xk):
    """Force and moment at every contact (units m0 g and m0 g d0) from settling and from a unit du/V shear."""
    n = len(p)
    inv = np.linalg.inv((p * p).sum() * np.eye(3) - p.T @ p + n * OWN_INERTIA * np.eye(3))
    fs = np.zeros((n, 3))
    fs[:, 2] = drag - 1.0
    ft = np.zeros((n, 3))
    ft[:, 2] = 2.0 * drag * xk
    fs, cs = equilibrate(p, fs, inv)
    ft, ct = equilibrate(p, ft, inv)
    resid = max(np.abs(fs.sum(0)).max(), np.abs((np.cross(p, fs) + cs).sum(0)).max(),
                np.abs(ft.sum(0)).max(), np.abs((np.cross(p, ft) + ct).sum(0)).max())
    s = subtree_sum(np.hstack([fs, np.cross(p, fs) + cs, ft, np.cross(p, ft) + ct, np.ones((n, 1))]),
                    parent, levels)[kids]
    c = 0.5 * (p[kids] + p[parent[kids]])
    f_s, f_t = s[:, 0:3], s[:, 6:9]
    return f_s, s[:, 3:6] - np.cross(c, f_s), f_t, s[:, 9:12] - np.cross(c, f_t), s[:, 12], resid


def analyse(c, rng):
    p0 = c['p'] - c['p'].mean(0)
    n = len(p0)
    parent, levels = tree(c['e'], n)
    kids = np.flatnonzero(parent >= 0)
    d = span(p0)
    rec = dict(n=n, span=d, rg=float(math.sqrt((p0 * p0).sum(1).mean())), depth=len(levels), models={}, diag={})
    for how in ORIENTATIONS:
        p, v = orient(p0, c['n'], how, rng)
        areas, covered = projected_areas(p, v)
        ph = 2.0 * math.pi * rng.random()
        xk = (p[:, 0] * math.cos(ph) + p[:, 1] * math.sin(ph)) / d     # local du/V per unit du/V across the span
        rec['diag'][how] = dict(area_ratio=covered / (math.pi / 4.0 * d * d),
                                unshadowed_fraction=float(areas['shadow'].sum() / areas['own'].sum()),
                                thickness_ratio=float((p[:, 2].max() - p[:, 2].min() + 1.0) / d))
        for name in SHADOWS:
            drag = n * areas[name] / areas[name].sum()
            fs, ms, ft, mt, cnt, resid = junction_loads(p, parent, levels, kids, drag, xk)
            share = np.maximum(cnt, n - cnt) / n
            out = dict(F=[], M=[], shareF=[], shareM=[], resid=float(resid))
            for x in X_GRID:
                fa = np.maximum(np.linalg.norm(fs + x * ft, axis=1), np.linalg.norm(fs - x * ft, axis=1))
                ma = np.maximum(np.linalg.norm(ms + x * mt, axis=1), np.linalg.norm(ms - x * mt, axis=1))
                kf, km = int(np.argmax(fa)), int(np.argmax(ma))
                out['F'].append(float(fa[kf]))
                out['M'].append(float(ma[km]))
                out['shareF'].append(float(share[kf]))
                out['shareM'].append(float(share[km]))
            rec['models'][(how, name)] = out
    return rec


# ---------------------------------------------------------------- checks

def selftest(rng):
    """Geometry, equilibrium and tree-sum checks on one 128-plate aggregate."""
    cl = monomers(128, rng)
    while len(cl) > 1:
        cl = [merge(cl[2 * k], cl[2 * k + 1], rng) for k in range(len(cl) // 2)]
    c = cl[0]
    p, e, n = c['p'], c['e'], len(c['p'])
    dist = np.linalg.norm(p[e[:, 0]] - p[e[:, 1]], axis=1)
    pair = np.sqrt(((p[:, None, :] - p[None, :, :]) ** 2).sum(-1)) + 2.0 * np.eye(n)
    parent, levels = tree(e, n)
    kids = np.flatnonzero(parent >= 0)
    p2, v2 = orient(p, c['n'], 'broadside', rng)
    areas, _ = projected_areas(p2, v2)
    drag = n * areas['shadow'] / areas['shadow'].sum()
    xk = p2[:, 0] / span(p2)
    fs, ms, ft, mt, cnt, resid = junction_loads(p2, parent, levels, kids, drag, xk)
    # brute force: walk every node up to the root to find each subtree, and sum the loads directly
    inv = np.linalg.inv((p2 * p2).sum() * np.eye(3) - p2.T @ p2 + n * OWN_INERTIA * np.eye(3))
    f = np.zeros((n, 3))
    f[:, 2] = drag - 1.0
    feq, ceq = equilibrate(p2, f, inv)
    worst = 0.0
    for idx, k in enumerate(kids[:40]):
        member = np.zeros(n, bool)
        for q in range(n):
            a = q
            while a >= 0 and a != k:
                a = parent[a]
            member[q] = a == k
        cp = 0.5 * (p2[k] + p2[parent[k]])
        fb = feq[member].sum(0)
        mb = (np.cross(p2[member] - cp, feq[member]) + ceq[member]).sum(0)
        worst = max(worst, np.abs(fb - fs[idx]).max(), np.abs(mb - ms[idx]).max(), abs(member.sum() - cnt[idx]))
    fu, mu, _, _, _, _ = junction_loads(p2, parent, levels, kids, np.ones(n), xk)
    mu = np.hstack([fu, mu])                                 # uniform drag per unit weight leaves no load at all
    return dict(edges=int(len(e)), contact_distance_error=float(np.abs(dist - 1.0).max()),
                closest_pair=float(pair.min()), equilibrium_residual=float(resid),
                tree_sum_vs_brute_force=float(worst), uniform_drag_max_moment=float(np.abs(mu).max()),
                passed=bool(len(e) == n - 1 and np.abs(dist - 1.0).max() < 1e-9 and pair.min() > 1.0 - 1e-9
                            and resid < 1e-8 and worst < 1e-8 and np.abs(mu).max() < 1e-9))


# ---------------------------------------------------------------- statistics

def ols(x, y):
    a = np.vstack([np.ones_like(x), x]).T
    return np.linalg.lstsq(a, y, rcond=None)[0]


def size_table(recs, key, xi):
    """Per size: n, log span, and log of the normalised loads F/(m g) and M/(m g D) at du/V = X_GRID[xi]."""
    rows = {}
    for r in recs:
        m = r['models'][key]
        rows.setdefault(r['n'], []).append((math.log(r['span']), math.log(m['F'][xi] / r['n']),
                                            math.log(m['M'][xi] / (r['n'] * r['span']))))
    return {n: np.array(v) for n, v in sorted(rows.items())}


def fit(tab, min_n, rng, boot=BOOT):
    """Slopes against log span of log n, log F/(mg), log M/(mgD), fitted on per-size means, with bootstrap SDs."""
    ns = [n for n in tab if n >= min_n]

    def slopes(sel):
        x = np.array([sel[k][:, 0].mean() for k in range(len(ns))])
        out = [ols(x, np.log(np.array(ns, float)))]
        for col in (1, 2):
            out.append(ols(x, np.array([sel[k][:, col].mean() for k in range(len(ns))])))
        return out

    base = slopes([tab[n] for n in ns])
    draws = np.array([[s[1] for s in slopes([tab[n][rng.integers(0, len(tab[n]), len(tab[n]))] for n in ns])]
                      for _ in range(boot)])
    sd = draws.std(0) if boot else np.zeros(3)
    return dict(sizes=ns, mass_exponent=float(base[0][1]), mass_exponent_sd=float(sd[0]),
                force_norm_exponent=float(base[1][1]), force_norm_exponent_sd=float(sd[1]),
                force_norm_intercept=float(base[1][0]),
                moment_norm_exponent=float(base[2][1]), moment_norm_exponent_sd=float(sd[2]),
                moment_norm_intercept=float(base[2][0]),
                force_exponent=float(base[0][1] + base[1][1]),
                moment_exponent=float(base[0][1] + base[2][1] + 1.0))


# ---------------------------------------------------------------- main

def build_all(levels, realisations, per_size, pair_finals, rng, rng_a, log):
    recs, finals = [], []
    for real in range(realisations):
        t0 = time.process_time()
        cl = monomers(2 ** levels, rng)
        for k in range(1, levels + 1):
            cl = [merge(cl[2 * q], cl[2 * q + 1], rng) for q in range(len(cl) // 2)]
            if k >= 3:
                take = len(cl) if k == levels else min(len(cl), -(-per_size // realisations))
                for c in cl[:take]:
                    recs.append(analyse(c, rng_a))
        finals.append(cl[0])
        log(f'realisation {real + 1}/{realisations}: {time.process_time() - t0:.1f} s cpu')
    if pair_finals:
        for q in range(len(finals) // 2):
            t0 = time.process_time()
            recs.append(analyse(merge(finals[2 * q], finals[2 * q + 1], rng), rng_a))
            log(f'pair {q + 1}: {2 ** (levels + 1)} plates, {time.process_time() - t0:.1f} s cpu')
    return recs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--levels', type=int, default=12)
    ap.add_argument('--realisations', type=int, default=20)
    ap.add_argument('--per-size', type=int, default=40)
    ap.add_argument('--no-pairs', action='store_true', help='skip the 2^(levels+1) aggregates')
    ap.add_argument('--quick', action='store_true', help='print only; do not write the product')
    args = ap.parse_args()
    cpu0, wall0 = time.process_time(), time.time()
    log = lambda s: print(s, flush=True)
    rng, rng_a = np.random.default_rng(SEED), np.random.default_rng(SEED + 1)
    check = selftest(np.random.default_rng(SEED + 2))
    log(f'selftest {check}')
    recs = build_all(args.levels, args.realisations, args.per_size, not args.no_pairs, rng, rng_a, log)
    rng_b = np.random.default_rng(SEED + 3)

    # size table and fits
    keys = [(o, s) for o in ORIENTATIONS for s in SHADOWS]
    sizes = sorted({r['n'] for r in recs})
    per_size = []
    for n in sizes:
        rs = [r for r in recs if r['n'] == n]
        row = dict(n=n, count=len(rs), span_d0_median=float(np.median([r['span'] for r in rs])),
                   rg_d0_median=float(np.median([r['rg'] for r in rs])),
                   tree_depth_median=float(np.median([r['depth'] for r in rs])))
        for o in ORIENTATIONS:
            row[f'area_ratio_{o}'] = float(np.median([r['diag'][o]['area_ratio'] for r in rs]))
            row[f'unshadowed_fraction_{o}'] = float(np.median([r['diag'][o]['unshadowed_fraction'] for r in rs]))
            row[f'thickness_ratio_{o}'] = float(np.median([r['diag'][o]['thickness_ratio'] for r in rs]))
        for o, s in keys:
            m = [r['models'][(o, s)] for r in rs]
            fn = np.array([q['F'][0] / r['n'] for q, r in zip(m, rs)])
            mn = np.array([q['M'][0] / (r['n'] * r['span']) for q, r in zip(m, rs)])
            row[f'{o}/{s}'] = dict(force_over_mg_median=float(np.median(fn)),
                                   force_over_mg_quartiles=[float(q) for q in np.percentile(fn, [25, 75])],
                                   moment_over_mgD_median=float(np.median(mn)),
                                   moment_over_mgD_quartiles=[float(q) for q in np.percentile(mn, [25, 75])],
                                   largest_remnant_force_median=float(np.median([q['shareF'][0] for q in m])),
                                   largest_remnant_moment_median=float(np.median([q['shareM'][0] for q in m])))
        per_size.append(row)
        sh, own = row['broadside/shadow'], row['broadside/own']
        log(f"n {n:5d} x{len(rs):3d} span {row['span_d0_median']:6.1f} d0  area ratio "
            f"{row['area_ratio_broadside']:.2f}  F/mg shadow {sh['force_over_mg_median']:.3f} "
            f"own {own['force_over_mg_median']:.4f}  M/mgD shadow "
            f"{sh['moment_over_mgD_median']:.4f} own {own['moment_over_mgD_median']:.5f}"
            f"  unshadowed {row['unshadowed_fraction_broadside']:.2f}/{row['unshadowed_fraction_random']:.2f}"
            f"  thickness {row['thickness_ratio_broadside']:.2f}  remnant F/M shadow "
            f"{sh['largest_remnant_force_median']:.2f}/{sh['largest_remnant_moment_median']:.2f}")

    fits = {}
    for o, s in keys:
        for mn in FIT_MIN_N:
            fits[f'{o}/{s}/from{mn}'] = fit(size_table(recs, (o, s), 0), mn, rng_b)
    turb_fits = {}
    for s in SHADOWS:
        turb_fits[s] = [fit(size_table(recs, ('broadside', s), xi), FIT_MIN_N[0], rng_b, boot=0) for xi in
                        range(len(X_GRID))]

    # largest remnant at the first failure, sizes of 256 plates and up
    remnant = {}
    for o, s in keys:
        big = [r['models'][(o, s)] for r in recs if r['n'] >= 256]
        out = {}
        for crit, k in (('force', 'shareF'), ('moment', 'shareM')):
            for xi, tag in ((0, 'settling'), (int(np.argmin(np.abs(X_GRID - 0.1))), 'du_over_V_0.1')):
                sh = np.array([q[k][xi] for q in big])
                out[f'{crit}/{tag}'] = dict(median=float(np.median(sh)),
                                            quartiles=[float(q) for q in np.percentile(sh, [25, 75])],
                                            p10_p90=[float(q) for q in np.percentile(sh, [10, 90])],
                                            fraction_above_0p9=float(np.mean(sh > 0.9)), count=int(sh.size))
        remnant[f'{o}/{s}'] = out

    # Earth calibration and lunar spans
    moon, binding = flakes.moon_column()
    earth = flakes.earth_column()
    g_ratio = earth.gravity / moon.gravity
    laws = {'simulated': None, **{k: v for k, v in ts.MASS_LAWS.items()}}
    spans = {}
    for o, s in keys:
        for mn in FIT_MIN_N:
            fr = fits[f'{o}/{s}/from{mn}']
            rows = {}
            for law, ab in laws.items():
                b = fr['mass_exponent'] if ab is None else ab[1]
                rf = g_ratio ** (1.0 / (b + fr['force_norm_exponent']))
                rm = g_ratio ** (1.0 / (b + 1.0 + fr['moment_norm_exponent']))
                rows[law] = dict(mass_exponent=b, load_analogy_ratio=g_ratio ** (1.0 / b), force_ratio=rf,
                                 moment_ratio=rm, force_spans_cm=[100 * d * rf for d in EARTH_SPANS_M],
                                 moment_spans_cm=[100 * d * rm for d in EARTH_SPANS_M],
                                 load_analogy_spans_cm=[100 * d * g_ratio ** (1.0 / b) for d in EARTH_SPANS_M])
            spans[f'{o}/{s}/from{mn}'] = rows

    def phi(s, crit, d_m, x):
        """Normalised peak load at span d_m (m) and du/V = x, from the per-x power-law fits (broadside)."""
        col = 'force' if crit == 'force' else 'moment'
        c = np.array([f[f'{col}_norm_intercept'] for f in turb_fits[s]])
        e = np.array([f[f'{col}_norm_exponent'] for f in turb_fits[s]])
        ld = math.log(d_m / D0_M)
        return math.exp(np.interp(x, X_GRID, c) + np.interp(x, X_GRID, e) * ld)

    def load(world, d_m, eps, law, s, crit, form):
        q = ts.load(world, d_m, eps, law)
        x = q['velocity_increment_m_s'] / q['fall_speed_m_s']
        lever = d_m if crit == 'moment' else 1.0
        if form == 'network':
            return phi(s, crit, d_m, x) * q['weight_n'] * lever, x
        return phi(s, crit, d_m, 0.0) * q['weight_n'] * lever * (1.0 + 2.0 * x), x

    turbulence = []
    for s in SHADOWS:
        for law in ts.MASS_LAWS:
            for eps in DISSIPATION:
                for d_e in EARTH_SPANS_M:
                    for crit in ('force', 'moment'):
                        for form in ('network', 'whole_flake'):
                            target, x_e = load(earth, d_e, eps, law, s, crit, form)
                            lo, hi = d_e * 0.5, d_e * 40.0
                            for _ in range(80):
                                mid = math.sqrt(lo * hi)
                                if load(moon, mid, eps, law, s, crit, form)[0] < target:
                                    lo = mid
                                else:
                                    hi = mid
                            d_m = math.sqrt(lo * hi)
                            turbulence.append(dict(drag=s, mass_law=law, dissipation=eps, earth_span_cm=100 * d_e,
                                                   criterion=crit, form=form, moon_span_cm=100 * d_m,
                                                   ratio=d_m / d_e, du_over_V_earth=x_e,
                                                   du_over_V_moon=load(moon, d_m, eps, law, s, crit, form)[1],
                                                   earth_load=target))
    # settling-only lunar spans with the per-x fits at x = 0 (same machinery, for the table beside the turbulence)
    strengths = {}
    mitchell = 'mitchell1996_mixed_aggregates'
    for s in SHADOWS:
        rows = {}
        for d_e in EARTH_SPANS_M:
            q = ts.load(earth, d_e, DISSIPATION[0], mitchell)
            rows[f'{100 * d_e:.0f}cm'] = dict(weight_n=q['weight_n'],
                                              junction_force_n=phi(s, 'force', d_e, 0.0) * q['weight_n'],
                                              junction_moment_n_m=phi(s, 'moment', d_e, 0.0) * q['weight_n'] * d_e)
        strengths[s] = rows

    # derived quantities quoted in the note
    derived = dict(dense_limit={}, lever_arm_mm={}, reynolds={},
                   simulated_span_d0=[min(r['span_d0_median'] for r in per_size if r['n'] >= FIT_MIN_N[0]),
                                      max(r['span_d0_median'] for r in per_size)])
    derived['largest_simulated_span_cm'] = derived['simulated_span_d0'][1] * D0_M * 100.0
    for law, ab in ts.MASS_LAWS.items():          # loads in proportion to m g D: the densely shadowed limit
        r = g_ratio ** (1.0 / (ab[1] + 1.0))
        derived['dense_limit'][law] = dict(moment_ratio=r, moment_spans_cm=[100 * d * r for d in EARTH_SPANS_M])
    for s in SHADOWS:                             # peak moment over peak force at the calibration flakes
        derived['lever_arm_mm'][s] = {k: 1e3 * v['junction_moment_n_m'] / v['junction_force_n']
                                      for k, v in strengths[s].items()}
    moon_spans = [d / 100.0 for k in ('broadside/shadow/from64', 'broadside/own/from64')
                  for d in spans[k][mitchell]['moment_spans_cm']]
    for name, world, ds in (('earth', earth, EARTH_SPANS_M), ('moon', moon, moon_spans)):
        st = world.state(world.z_freeze)
        nu = st['eta'] / st['rho']
        for d in ds:
            q = ts.load(world, d, DISSIPATION[0], mitchell)
            derived['reynolds'][f'{name}_{100 * d:.1f}cm'] = dict(fall_speed_m_s=q['fall_speed_m_s'],
                                                                 flake=q['fall_speed_m_s'] * d / nu,
                                                                 plate=q['fall_speed_m_s'] * D0_M / nu)

    cpu = time.process_time() - cpu0
    usage = resource.getrusage(resource.RUSAGE_SELF)
    files = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(ts.__file__), Path(flakes.__file__))}
    product = dict(
        schema=SCHEMA,
        evidence=('Calculation on idealised aggregates (equal randomly oriented plates joined by single contacts from '
                  'hierarchical ballistic cluster-cluster aggregation) with geometric drag-sharing rules; no measured '
                  'junction strength, flow field or lunar breakup. The Earth calibration assumes 4-6 cm flakes sit at '
                  'their breaking point.'),
        reading_rule=('Loads are in units of the plate weight m0 g (force) and m0 g d0 (moment); normalised loads are '
                      'F/(m g) and M/(m g D) with m the aggregate mass and D its maximum dimension. "shadow" gives '
                      'each pixel of projected area to the lowest plate (opaque, the upper bound on uneven drag), '
                      '"own" to every plate over it (no shadowing, the lower bound), "partial" in shares 0.5^rank. '
                      'Ratios are lunar over Earth span at equal peak junction load; spans come from power-law fits '
                      'over sizes from N plates up, N named in each key (from64, from256).'),
        assumptions=['The aggregates are equal plates joined by single contacts, so the contacts form a tree; natural '
                     'flakes interlock and share a moment over several contacts.',
                     'The aggregates grow by joining equal clusters; natural flakes grow largely by sweeping much '
                     'smaller crystals and flakes, a history that changes the contact tree and may move the largest '
                     'remnant share.',
                     'The aggregates are sparser in projection than natural flakes (area ratio 0.08-0.14 against '
                     '0.35), which weakens their shadowing; the dense limit, moments in proportion to m g D, brackets '
                     'that side.',
                     'Drag is shared by geometric rules at flake Reynolds numbers of 3000-6700; wakes, flow through '
                     'the porous flake and edge effects are not resolved. The orientation is fixed and the loads are '
                     'instantaneous.',
                     'Every junction is equally strong at every flake size, and Earth flakes of 4-6 cm sit at their '
                     'breaking point; only the first failure is followed.'],
        inputs=dict(seed=SEED, levels=args.levels, realisations=args.realisations, per_size=args.per_size,
                    pairs=not args.no_pairs, plate_diameter_m=D0_M, pixel_d0=PIXEL, own_inertia_m0d0sq=OWN_INERTIA,
                    shadows=SHADOWS, orientations=list(ORIENTATIONS), du_over_V_grid=X_GRID.tolist(),
                    fit_min_n=list(FIT_MIN_N), earth_spans_m=list(EARTH_SPANS_M), dissipation=list(DISSIPATION),
                    gravity=dict(earth=earth.gravity, moon=moon.gravity, ratio=g_ratio),
                    air_density_freezing_level=dict(earth=earth.state(earth.z_freeze)['rho'],
                                                    moon=moon.state(moon.z_freeze)['rho']),
                    mass_laws={k: list(v) for k, v in ts.MASS_LAWS.items()},
                    observed_area_ratio_large_aggregates=flakes.AREA_RATIO,
                    climate_binding=binding.get('inputs', {}), imported_files_sha256=files,
                    numpy=np.__version__, python=sys.version.split()[0]),
        selftest=check, per_size=per_size, fits=fits, spans=spans, remnant=remnant, turbulence=turbulence,
        turbulence_fits={s: [dict(du_over_V=float(x), force_norm_exponent=f['force_norm_exponent'],
                                  force_norm_intercept=f['force_norm_intercept'],
                                  moment_norm_exponent=f['moment_norm_exponent'],
                                  moment_norm_intercept=f['moment_norm_intercept'])
                             for x, f in zip(X_GRID, turb_fits[s])] for s in SHADOWS},
        implied_strengths_mitchell=strengths, derived=derived,
        max_equilibrium_residual=max(q['resid'] for r in recs for q in r['models'].values()),
        compute=dict(cpu_seconds=cpu, wall_seconds=time.time() - wall0, max_rss_mb=usage.ru_maxrss / 1024.0),
        producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    for k in ('broadside/shadow/from64', 'broadside/partial/from64', 'broadside/own/from64',
              'random/shadow/from64', 'random/own/from64', 'broadside/shadow/from256', 'broadside/own/from256'):
        f = fits[k]
        log(f"{k}: b {f['mass_exponent']:.3f}±{f['mass_exponent_sd']:.3f}  F/mg exp {f['force_norm_exponent']:+.3f}"
            f"±{f['force_norm_exponent_sd']:.3f}  M/mgD exp {f['moment_norm_exponent']:+.3f}"
            f"±{f['moment_norm_exponent_sd']:.3f}  F~D^{f['force_exponent']:.2f} M~D^{f['moment_exponent']:.2f}")
        sp = spans[k]
        for law in ('simulated', mitchell):
            log(f"   {law}: ratio force {sp[law]['force_ratio']:.2f} moment {sp[law]['moment_ratio']:.2f} "
                f"(analogy {sp[law]['load_analogy_ratio']:.2f})")
    log(f'cpu {cpu:.1f} s, max rss {usage.ru_maxrss / 1024:.0f} MB')
    if not args.quick:
        PRODUCT.write_text(json.dumps(product, separators=(',', ':'), allow_nan=False) + '\n')
        log(f'wrote {PRODUCT}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
