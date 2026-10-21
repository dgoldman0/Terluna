#!/usr/bin/env python3
"""How large snowflakes grow in the sky winter: a column of falling snow that grows from the vapour, sticks together
and breaks apart, the Moon against Earth.

    python -m atmosphere.snow.aggregation        # a few minutes: results/aggregation.json

Crystals of 0.1-0.4 mm enter at the cloud top at a set number concentration and fall at their own speeds through a
winter column (flakes.Column) to the freezing level, where the snow arriving is counted once the column has filled.
The particles are super-particles (Shima and co-authors 2009), each standing for many real ones, at their own heights.
They grow from the vapour by flakes' diffusion law as far as the cloud supplies it: a snowfall of a set rate R
(kg/m2/s) needs the cloud to condense R spread through its depth, shared among the particles at each height by their
diffusional demand. They meet the others at their height, sweeping the area (A1^(1/2) + A2^(1/2))^2 of their projected
areas (Mitchell 1988) at their difference in speed. A meeting sticks with Lin, Farley and Orville's (1983)
efficiency, exp(0.09 (T - 273.15)). One that does not stick breaks fragments off the smaller particle by Phillips, Yano
and Khain's (2017) law for snow, as Sotiropoulou and co-authors set it out (Atmospheric Chemistry and Physics
Discussions, acp-2020-786, appendix B): F = a A (1 - exp(-(C K0 / (a A))^g)), K0 = m1 m2 / (m1 + m2) du^2 with du the
two speeds' difference, a = pi D^2 for the smaller particle's size D, held within the law's 0.5-5 mm,
A = 1.58e7 (1 + 1.33e-4 / D^1.5) m^-2, C = 7.08e6 x 3.5e-3, g = 0.5 and at most 100 fragments a collision. The
fragments are 0.1 mm crystals taken from the smaller particle, at most half its mass. Gravity reaches the result
through the fall speeds: how often particles meet, how hard, and how long they take to fall through the cloud.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np

from atmosphere.rain import drops
from atmosphere.snow import flakes as sf

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCHEMA = 'terluna.atmosphere.snow-aggregation/1'
OUT = HERE / 'results' / 'aggregation.json'
CODE = ('atmosphere/snow/aggregation.py', 'atmosphere/snow/flakes.py', 'atmosphere/rain/drops.py',
        'atmosphere/rain/targets.py')
START_M = (1.0e-4, 4.0e-4)                 # the crystals at the cloud top, log-uniform between these sizes
FRAGMENT_M = 1.0e-4                        # the fragments' size
BREAKUP = dict(a_scale=1.58e7, a_size=1.33e-4, c=7.08e6 * 3.5e-3, gamma=0.5, size_m=(5e-4, 5e-3), most=100.0)


def speeds(d, m, gravity, s, area_ratio=sf.AREA_RATIO):
    """Fall speeds (m/s) and Reynolds numbers of particles of sizes d and masses m (arrays), Heymsfield and
    Westbrook (2010) as flakes.fall_speed."""
    x = 8.0 * m * gravity * s['rho'] / (math.pi * math.sqrt(area_ratio) * s['eta'] ** 2)
    re = sf.DELTA0 ** 2 / 4.0 * (np.sqrt(1.0 + 4.0 * np.sqrt(x) / (sf.DELTA0 ** 2 * math.sqrt(sf.C0))) - 1.0) ** 2
    return s['eta'] * re / (s['rho'] * d), re


def deposition(d, re, s, capacitance_ratio=sf.CAPACITANCE):
    """kg/s each particle gains from the vapour, flakes.deposition_rate over arrays."""
    t = s['t']
    ls = drops.latent_heat(sf.T0) + sf.LATENT_FUSION
    esi = sf.ice_saturation_pressure(t)
    resistance = (ls / (sf.R_VAPOUR * t) - 1.0) * ls / (s['ka'] * t) + sf.R_VAPOUR * t / (s['dv'] * esi)
    x = (s['eta'] / (s['rho'] * s['dv'])) ** (1.0 / 3.0) * np.sqrt(re)
    vent = np.where(x < 1.0, 1.0 + 0.14 * x * x, 0.86 + 0.28 * x)
    return 4.0 * math.pi * capacitance_ratio * d * vent * (s['e'] / esi - 1.0) / resistance


def fragments(smaller_m, energy, mass_to_size=sf.flake_size):
    """Fragments a non-sticking collision breaks off the smaller particle (Phillips, Yano and Khain 2017, snow)."""
    d = np.clip(mass_to_size(smaller_m), *BREAKUP['size_m'])
    a_area = math.pi * d * d
    scale = BREAKUP['a_scale'] * (1.0 + BREAKUP['a_size'] / d ** 1.5) * a_area
    return np.minimum(scale * (1.0 - np.exp(-(BREAKUP['c'] * energy / scale) ** BREAKUP['gamma'])), BREAKUP['most'])


def collide(m, xi, u, d, idx, volume, temperature, dt, rng, breakup, m_frag,
            area_ratio=sf.AREA_RATIO, mass_to_size=sf.flake_size):
    """One step of collisions among the particles idx sharing one height (volume, m3): masses and multiplicities are
    updated in place; returns the fragments made (real particles)."""
    n = idx.size
    if n < 2:
        return 0.0
    order = rng.permutation(idx)
    j, k = order[: n // 2], order[n // 2: 2 * (n // 2)]
    du = np.abs(u[j] - u[k])
    kernel = area_ratio * math.pi / 4.0 * (d[j] + d[k]) ** 2 * du
    p = np.maximum(xi[j], xi[k]) * kernel * dt / volume * (n * (n - 1) / 2) / (n // 2)
    gamma = np.floor(p) + (rng.random(p.size) < p - np.floor(p))
    hit = gamma > 0
    if not hit.any():
        return 0.0
    j, k, gamma, du = j[hit], k[hit], gamma[hit], du[hit]
    swap = xi[j] < xi[k]
    a, b = np.where(swap, k, j), np.where(swap, j, k)                 # a has the larger multiplicity
    gt = np.minimum(gamma, np.floor(xi[a] / xi[b]))
    ok = gt >= 1
    a, b, gt, du = a[ok], b[ok], gt[ok], du[ok]
    stick = rng.random(a.size) < math.exp(0.09 * (temperature - sf.T0))
    # sticking (Shima and co-authors 2009)
    sa, sb, sg = a[stick], b[stick], gt[stick]
    merged = m[sb] + sg * m[sa]
    rest = xi[sa] - sg * xi[sb]
    split = rest <= 0
    xi[sa[~split]] = rest[~split]
    m[sb[~split]] = merged[~split]
    half = np.floor(xi[sb[split]] / 2.0)
    xi[sa[split]], xi[sb[split]] = half, xi[sb[split]] - half
    m[sa[split]] = merged[split]
    m[sb[split]] = merged[split]
    if not breakup:
        return 0.0
    ba, bb, bg, bu = a[~stick], b[~stick], gt[~stick], du[~stick]
    if ba.size == 0:
        return 0.0
    count = fragments(np.minimum(m[ba], m[bb]), m[ba] * m[bb] / (m[ba] + m[bb]) * bu ** 2,
                      mass_to_size=mass_to_size)
    b_small = m[bb] <= m[ba]                                           # each b particle met bg of a's
    lost_b = np.minimum(bg * count * m_frag, 0.5 * m[bb])
    lost_a = np.minimum(bg * xi[bb] * count * m_frag / xi[ba], 0.5 * m[ba])   # bg * xi[b] of a's met one each
    m[bb[b_small]] -= lost_b[b_small]
    m[ba[~b_small]] -= lost_a[~b_small]
    return float(np.sum(lost_b[b_small] * xi[bb[b_small]]) + np.sum(lost_a[~b_small] * xi[ba[~b_small]])) / m_frag


def column(world, top_c, number_per_m3, snowfall_mm_h, bins=30, inject=8, breakup=True, seed=0, fill=1.5, gather=1.5,
           geometry=None, observer=None, max_step_s=30.0, arrival_observer=None):
    """Snow falling from the cloud top (top_c, °C) to the freezing level: the masses and numbers of the particles
    arriving there over the gathering time once the column has filled, per square metre, and the time they cover."""
    rng = np.random.default_rng(seed)
    # Optional geometry and observation are used by the subsequent population study.
    # Omitting them preserves the original run's equations and random sequence.
    mass_of = sf.flake_mass if geometry is None else geometry.mass
    size_of = sf.flake_size if geometry is None else geometry.size
    area_ratio = sf.AREA_RATIO if geometry is None else geometry.area_ratio
    capacitance = sf.CAPACITANCE if geometry is None else geometry.capacitance_ratio
    g = world.gravity
    z_top, z_f = world.height_of(sf.T0 + top_c), world.z_freeze
    depth = z_top - z_f
    edges = np.linspace(z_f, z_top, bins + 1)
    dz = depth / bins
    states = [world.state(0.5 * (lo + hi)) for lo, hi in zip(edges[:-1], edges[1:])]
    production = snowfall_mm_h / 3600.0 / depth                         # kg/m3/s the cloud condenses
    m_frag = mass_of(FRAGMENT_M)
    s_top = states[-1]
    d_in = math.sqrt(START_M[0] * START_M[1])
    u_in = speeds(np.array([d_in]), np.array([mass_of(d_in)]), g, s_top, area_ratio)[0][0]
    dt = min(max_step_s, 0.5 * dz / max(1.5, 3.0 * u_in))
    crossing = depth / u_in
    m = np.empty(0)
    xi = np.empty(0)
    z = np.empty(0)
    arrived_m, arrived_xi = [], []
    t, made_total = 0.0, 0.0
    injected_mass = deposited_mass = exited_mass = discarded_fragment_mass = 0.0
    while t < (fill + gather) * crossing:
        d_new = np.exp(rng.uniform(math.log(START_M[0]), math.log(START_M[1]), inject))
        m = np.append(m, mass_of(d_new))
        xi = np.append(xi, np.full(inject, number_per_m3 * u_in * dt / inject))
        injected_mass += float(np.sum(m[-inject:] * xi[-inject:]))
        z = np.append(z, z_top - rng.uniform(0.0, u_in * dt, inject))
        layer = np.clip(((z - z_f) / dz).astype(int), 0, bins - 1)
        d = size_of(m)
        u = np.empty_like(m)
        broken = []
        for b in range(bins):
            here = layer == b
            if not here.any():
                continue
            s = states[b]
            u[here], re = speeds(d[here], m[here], g, s, area_ratio)
            demand = np.maximum(deposition(d[here], re, s, capacitance), 0.0) * dt
            total = float(np.sum(demand * xi[here]))
            supply = production * dz * dt
            addition = demand * (min(1.0, supply / total) if total > 0 else 0.0)
            deposited_mass += float(np.sum(addition * xi[here]))
            m[here] += addition
            made = collide(m, xi, u, d, np.flatnonzero(here), dz, s['t'], dt, rng, breakup, m_frag,
                           area_ratio=area_ratio, mass_to_size=size_of)
            if made >= 1.0:
                broken.append((made, edges[b] + 0.5 * dz))
            else:
                discarded_fragment_mass += made * m_frag
        for count, height in broken:
            made_total += count
            m = np.append(m, m_frag)
            xi = np.append(xi, count)
            z = np.append(z, height)
            u = np.append(u, 0.0)
        z = z - u * dt
        out = z < z_f
        exited_mass += float(np.sum(m[out] * xi[out]))
        if arrival_observer is not None and out.any():
            arrival_observer(m[out], xi[out], t, dt)
        if t >= fill * crossing and out.any():
            arrived_m.append(m[out])
            arrived_xi.append(xi[out])
        keep = (~out) & (xi > 0) & (m > 0)
        m, xi, z = m[keep], xi[keep], z[keep]
        if observer is not None and t >= fill * crossing:
            observer(m, xi, z, t, dt)
        t += dt
    return dict(masses=np.concatenate(arrived_m) if arrived_m else np.empty(0),
                numbers=np.concatenate(arrived_xi) if arrived_xi else np.empty(0),
                seconds=gather * crossing, crossing_hours=crossing / 3600.0, fragments=made_total / t,
                mass_ledger_kg_m2=dict(injected=injected_mass, deposited=deposited_mass, exited=exited_mass,
                                     remaining=float(np.sum(m * xi)), discarded_fragment=discarded_fragment_mass,
                                     residual=injected_mass + deposited_mass - exited_mass
                                              - float(np.sum(m * xi)) - discarded_fragment_mass))


def summary(result, mass_to_size=sf.flake_size):
    """The snow arriving at the freezing level: its sizes by mass (the median and the sizes below which 90% and 99% of
    it lies), the largest flake arriving at least once an hour on each square metre, the snowfall it makes (mm of
    water an hour), the flakes arriving each hour on a square metre and the crystals' crossing time."""
    m, n = result['masses'], result['numbers']
    if m.size == 0:
        return dict(median_by_mass_mm=None, d90_by_mass_mm=None, d99_by_mass_mm=None, largest_hourly_mm=None,
                    snowfall_mm_h=0.0, flakes_per_m2_hour=0.0, crossing_hours=result['crossing_hours'])
    d = mass_to_size(m)
    order = np.argsort(d)
    d, m, n = d[order], m[order], n[order]
    cum = np.cumsum(m * n) / np.sum(m * n)
    at = lambda q: float(1e3 * d[min(np.searchsorted(cum, q), d.size - 1)])
    hours = result['seconds'] / 3600.0
    hourly = d[np.cumsum(n[::-1])[::-1] / hours >= 1.0]
    return dict(median_by_mass_mm=at(0.5), d90_by_mass_mm=at(0.9), d99_by_mass_mm=at(0.99),
                largest_hourly_mm=float(1e3 * hourly.max()) if hourly.size else None,
                snowfall_mm_h=float(np.sum(m * n) / hours), flakes_per_m2_hour=float(np.sum(n) / hours),
                crossing_hours=result['crossing_hours'], fragments_per_m2_s=result['fragments'])


def analyse(top_c=-15.0, snowfalls=(0.3, 1.0, 3.0), numbers=(1e4, 1e5)):
    moon, binding = sf.moon_column()
    earth = sf.earth_column()
    runs = {}
    for world in (moon, earth):
        rows = []
        for rate in snowfalls:
            for n in numbers:
                for breakup in (True, False):
                    rows.append(dict(snowfall_set_mm_h=rate, crystals_per_litre=n / 1e3, breakup=breakup,
                                     **summary(column(world, top_c, n, rate, breakup=breakup))))
        runs[world.name] = rows
    files = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()[:16] for p in CODE}
    return dict(
        schema=SCHEMA, cloud_top_c=top_c, runs=runs, breakup_law=dict(BREAKUP, fragment_m=FRAGMENT_M),
        start_size_m=list(START_M),
        producer=dict(domain='atmosphere', files=files, inputs=binding['inputs']),
        evidence=('A one-dimensional super-particle model of snow in a cloud, with Earth-fitted laws for growth, fall, '
                  'sticking and fragmentation carried to lunar gravity and air through the fall speeds; one seed a '
                  'case; a still cloud whose condensation feeds a set snowfall. Sizes for comparing the two worlds, not '
                  'a forecast of single flakes.'),
        reading_rule=('Sizes are maximum dimensions of the snow arriving at the freezing level once the column has '
                      'filled: the median, 90th and 99th percentiles by mass and the largest flake arriving at least '
                      'once an hour on a square metre. snowfall_set_mm_h is the condensation the cloud is given, '
                      'snowfall_mm_h what arrives; crystals_per_litre is the number entering at the cloud top; breakup '
                      'false lets meetings that do not stick bounce whole.'))


def main() -> int:
    out = analyse()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1, default=float) + '\n')
    for name, rows in out['runs'].items():
        for r in rows:
            print(f"{name} {r['snowfall_set_mm_h']:g} mm/h {r['crystals_per_litre']:g}/L breakup {r['breakup']}: "
                  f"median {r['median_by_mass_mm']}, d90 {r['d90_by_mass_mm']}, d99 {r['d99_by_mass_mm']}, "
                  f"hourly largest {r['largest_hourly_mm']}, arriving {r['snowfall_mm_h']:.2f} mm/h")
    print(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
