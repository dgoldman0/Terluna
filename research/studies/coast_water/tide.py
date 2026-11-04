"""The monthly tide in the shores' groundwater: how far inland, and with what lag, the seas' monthly tide moves the
water table, and how much the tide on a sloping shore raises the mean water table inland.

Linear theory (Ferris 1951): a sea level A cos(wt) at a fixed shore drives the water table A exp(-kx) cos(wt - kx)
inland, with k = sqrt(w Sy / 2T), so the swing falls to 1/e in 1/k and lags kx/w. On a sloping shore the waterline
moves A cot(beta) up and down the face, and the face takes the sea in faster on the rise, over the wetted slope, than
it lets it out on the fall, through a seepage face. Nielsen (1990) gives the steady part to first order in
eps = kA cot(beta): the mean water table inland stands eps A / 2 above the mean sea; the finite amplitude adds
A^2 / 4D over an aquifer D thick.

The numerical check solves unconfined Boussinesq flow across the shore, Sy dh/dt = d/dx(T dh/dx) - U, on a uniform
grid by backward Euler: nodes under the sea are held at its level, so the shoreline moves on the slope; nodes whose
water table would rise above the ground seep, held at the ground; and the plants draw U = E wherever the water table
lies within the reach R below the ground (a ramp of RAMP_M across the edge), nothing elsewhere. The inland end is
closed, as the coast-water study's belts are. T is constant (the coast-water study's) unless an aquifer depth is
given, when it is T (1 + h/D). Without the tide and on a steady rise the solution is belts.width's belt.

Lunar gravity enters through T, stage 2's at lunar gravity: the same ground carries the tide inland sqrt(6) times
less far than on Earth at the same period. The specific yield is a drained pore volume and holds at any gravity for a
water table deeper than the soil's capillary length, 3-6 m at lunar gravity in geography's Gardner soils; within that
of the ground the fringe lowers the apparent yield and the tide reaches further.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.linalg import solve_banded

from shared.constants import ANOMALISTIC_MONTH_DAYS, JULIAN_DAY, SYNODIC_MONTH_DAYS

MONTH_S = ANOMALISTIC_MONTH_DAYS * JULIAN_DAY        # the leading line M', the eccentricity's month
LUNAR_DAY_S = SYNODIC_MONTH_DAYS * JULIAN_DAY
# Specific yield by ground: the crust's fractured and vesicular basalt and its developed soil, 5-10% (Hunt 1996's
# Oahu lavas, mostly 5-10%; Woessner and Poeter's fractured basalt 2-10%, vesicular 5-15%); sandy ground, the
# averages of fine to coarse sand, 21-27% (Johnson 1967).
SPECIFIC_YIELD = dict(slow=(0.05, 0.10), central=(0.05, 0.10), sandy=(0.21, 0.27))
RAMP_M = 0.2                    # the plants' draw fades out over this depth at the reach's edge
ROOTS_M = 3.0                   # the shallower roots of the coast-water study's 3-5 m


def damping(period_s, specific_yield, transmissivity):
    """Ferris's wavenumber k (1/m) for a tide of the given period: sqrt(w Sy / 2T)."""
    return math.sqrt(math.pi * specific_yield / (period_s * transmissivity))


def ferris(x, amplitude, period_s, specific_yield, transmissivity):
    """The swing's amplitude (m) at x (m) inland of a fixed shore and its lag behind the sea (s)."""
    k = damping(period_s, specific_yield, transmissivity)
    x = np.asarray(x, float)
    return amplitude * np.exp(-k * x), k * x * period_s / (2.0 * math.pi)


def nielsen_overheight(amplitude, k, slope):
    """The mean water table's rise inland (m) to first order in eps = kA/slope: eps A / 2 (Nielsen 1990)."""
    return 0.5 * k * np.asarray(amplitude, float) ** 2 / np.asarray(slope, float)


def finite_amplitude_overheight(amplitude, depth):
    """The finite-amplitude rise A^2 / 4D (m) over an aquifer depth D (Nielsen 1990; Philip 1973)."""
    return np.asarray(amplitude, float) ** 2 / (4.0 * depth)


def sinusoid(amplitude, period_s=MONTH_S, phase=0.0):
    return lambda t: amplitude * math.cos(2.0 * math.pi * t / period_s + phase)


def lines(amplitudes, periods_s, phases):
    """A sea level made of tide lines: sum of A cos(2 pi t / P + phi)."""
    a, p, f = (np.asarray(v, float) for v in (amplitudes, periods_s, phases))
    return lambda t: float(np.sum(a * np.cos(2.0 * math.pi * t / p + f)))


def lunar_day_draw(mean, day_share, period_s=LUNAR_DAY_S, phase=0.0):
    """The plants' draw through the lunar day: day_share of the mean's water drawn while the Sun is up (half the
    lunar day), the rest at night."""
    day, night = 2.0 * day_share * mean, 2.0 * (1.0 - day_share) * mean
    return lambda t: day if math.cos(2.0 * math.pi * t / period_s + phase) > 0.0 else night


def cross_shore(x, ground, specific_yield, transmissivity, sea, draw=0.0, reach=math.inf, dt=None, window=MONTH_S,
                start=None, tolerance=1e-3, min_windows=4, max_windows=400, depth=math.inf, probe=None,
                record=None):
    """Boussinesq flow across a shore to a repeating state. x (m) is uniform and increasing inland, ground (m) the
    ground's height above the mean sea; sea(t) the sea's level (m) and draw a rate (m/s) or a function of time. Each
    window (one period of the forcing) is averaged; the run stops once the window means at the probe (default the
    inland end) change so little that their geometric tail is under the tolerance (m).

    Returns the last window's mean, lowest and highest water table, the plants' mean draw (m/s per metre of ground),
    the water that entered from the sea and the seepage face and that left by them (m3 per metre of shore per window),
    the number of windows, and the last change at the probe; with record (node indices), the water table there at
    every step of the last window and the steps' times."""
    x = np.asarray(x, float)
    ground = np.asarray(ground, float)
    n = x.size
    dx = float(x[1] - x[0])
    dt = dt or window / 240.0
    steps = max(int(round(window / dt)), 1)
    dt = window / steps
    rate = draw if callable(draw) else (lambda t, e=float(draw): e)
    h = np.minimum(np.zeros(n) if start is None else np.array(start, float), np.maximum(ground, sea(0.0)))
    probe = n - 1 if probe is None else probe
    t, history, changes = 0.0, [], []
    for w in range(max_windows):
        total = np.zeros(n)
        low, high = np.full(n, np.inf), np.full(n, -np.inf)
        drawn = np.zeros(n)
        inflow = outflow = 0.0
        series, times = [], []
        for _ in range(steps):
            t += dt
            eta = sea(t)
            wet = ground <= eta
            h = np.where(wet, h, np.minimum(h, ground))        # the falling sea's water above the ground runs off
            flooded = wet & (h < ground)                     # land the rising sea now covers fills to the ground
            inflow += float(specific_yield * dx * (ground[flooded] - h[flooded]).sum())
            u = rate(t) * np.clip((reach - (ground - h)) / RAMP_M + 0.5, 0.0, 1.0)
            u[wet] = 0.0
            tf = transmissivity * (1.0 + 0.5 * (h[1:] + h[:-1]) / depth) if math.isfinite(depth) else \
                np.full(n - 1, transmissivity)
            a = tf / (dx * dx)
            fixed, value = wet.copy(), np.where(wet, eta, 0.0)
            for _ in range(8):
                ab = np.zeros((3, n))
                ab[1] = specific_yield / dt
                ab[1, :-1] += a
                ab[1, 1:] += a
                ab[0, 1:] = -a
                ab[2, :-1] = -a
                rhs = specific_yield / dt * h - u
                ab[1, fixed] = 1.0
                rhs[fixed] = value[fixed]
                up = np.flatnonzero(fixed)
                ab[0, up[up + 1 < n] + 1] = 0.0               # a fixed row's upper neighbour
                ab[2, up[up - 1 >= 0] - 1] = 0.0              # and its lower
                new = solve_banded((1, 1), ab, rhs)
                over = ~fixed & (new > ground)
                if not over.any():
                    break
                fixed |= over
                value[over] = ground[over]
            q = -tf * (new[1:] - new[:-1]) / dx               # landward flux at each face (m2/s)
            into = np.where(fixed[:-1] & ~fixed[1:], q, 0.0) - np.where(~fixed[:-1] & fixed[1:], q, 0.0)
            inflow += float(np.clip(into, 0.0, None).sum() * dt)
            outflow += float(-np.clip(into, None, 0.0).sum() * dt)
            u = np.where(fixed, 0.0, u)
            drawn += u
            h = new
            if record is not None:
                series.append(h[record].copy())
                times.append(t)
            total += h
            np.minimum(low, h, out=low)
            np.maximum(high, h, out=high)
        mean = total / steps
        history.append(mean[probe])
        if len(history) >= 2:
            changes.append(abs(history[-1] - history[-2]))
        if len(changes) >= 2 and w + 1 >= min_windows:
            r = changes[-1] / max(changes[-2], 1e-30)
            tail = changes[-1] * r / (1.0 - r) if r < 0.98 else math.inf
            if changes[-1] < 1e-6 or tail < tolerance:
                break
    return dict(mean=mean, low=low, high=high, draw=drawn / steps, inflow=inflow, outflow=outflow,
                windows=w + 1, last_change=changes[-1] if changes else None, dx=dx, dt=dt,
                series=np.array(series) if record is not None else None, times=np.array(times))


def lag(series, times, period_s, level):
    """The phase lag (s) of a recorded water table's fundamental behind the sea's, over one period."""
    w = 2.0 * math.pi / period_s
    z = np.exp(-1j * w * times)
    sea = np.array([level(t) for t in times])
    phase = np.angle((series - series.mean()) @ z) - np.angle((sea - sea.mean()) @ z)
    return float(((-phase) % (2.0 * math.pi)) / w)


def edge(x, depth, reach, start=0.0):
    """Where the depth to water (m) below the ground first exceeds the reach inland of start (m): the belt's edge,
    linearly between nodes; x[-1] if it never does."""
    sel = np.flatnonzero((x >= start) & (depth > reach))
    if sel.size == 0:
        return float(x[-1])
    i = int(sel[0])
    if i == 0 or x[i - 1] < start:
        return float(max(x[i], start))
    d0, d1 = depth[i - 1] - reach, depth[i] - reach
    return float(x[i - 1] + (x[i] - x[i - 1]) * (-d0) / (d1 - d0))


def planar_grid(slope, amplitude, k, extent, nodes_per_damping=20, nodes_per_excursion=12, floor=0.25):
    """A uniform grid across a planar shore rising at slope: from beneath the lowest water to extent (m) inland,
    fine enough for the damping length 1/k and the waterline's excursion amplitude/slope."""
    excursion = amplitude / slope if amplitude > 0 else 0.0
    dx = 1.0 / (nodes_per_damping * k)
    if excursion > 0:
        dx = min(dx, max(excursion / nodes_per_excursion, floor))
    x0 = -1.5 * excursion - 4.0 * dx
    x = np.arange(x0, extent + dx, dx)
    return x, slope * x


def overheight_run(eps, specific_yield=0.07, transmissivity=5.2e-4, amplitude=1.0, period_s=MONTH_S, steps=None,
                   tolerance=2e-4):
    """The no-draw sloping shore at a given eps: the mean water table's rise inland over the amplitude (eta/A), the
    swing's half-range there and at 1/k, against Nielsen's first order, eps/2."""
    k = damping(period_s, specific_yield, transmissivity)
    slope = k * amplitude / eps
    x, ground = planar_grid(slope, amplitude, k, 4.0 / k + 1.5 * amplitude / slope)
    guess = min(0.5 * eps, 0.8) * amplitude
    steps = steps or (480 if eps > 3 else 240)
    out = cross_shore(x, ground, specific_yield, transmissivity, sinusoid(amplitude, period_s),
                      dt=period_s / steps, window=period_s, start=np.minimum(np.full(x.size, guess), ground),
                      tolerance=tolerance * amplitude)
    inland = out['mean'][-1]
    at = np.searchsorted(x, 1.0 / k)
    swing = 0.5 * (out['high'] - out['low'])
    return dict(eps=eps, slope=slope, overheight_over_amplitude=inland / amplitude, nielsen=0.5 * eps,
                swing_at_one_over_k=float(swing[at] / amplitude), ferris_at_one_over_k=math.exp(-1.0),
                windows=out['windows'], last_change_m=out['last_change'], nodes=int(x.size), dx_m=out['dx'],
                exchange_m3_per_m=out['inflow'], exchange_linear_m3_per_m=2.0 * amplitude
                * math.sqrt(transmissivity * specific_yield * period_s / (2.0 * math.pi)))


def overheight(eps, table):
    """eta/A at eps from the numerical table (eps, eta/A), log-linear between rows; Nielsen's eps/2 below its
    first row and the last row's value beyond its end."""
    e = np.asarray(table[0], float)
    f = np.asarray(table[1], float)
    eps = np.asarray(eps, float)
    out = np.interp(np.log(np.maximum(eps, 1e-12)), np.log(e), f)
    return np.where(eps < e[0], 0.5 * eps * f[0] / (0.5 * e[0]), out)


def belt_run(slope, transmissivity, specific_yield, draw, reach, amplitude, period_s=MONTH_S, sea=None,
             draw_t=None, window=None, extent=None, tolerance=2e-3, max_windows=800, steps=240, start_rise=None,
             profile=None, depth=math.inf):
    """A planar shore's belt with and without the tide. Returns both runs' belt edges (m from the mean waterline) by
    the mean depth and by the deepest water table of the window, the plants' draw (m3 a year per km of shore), the
    swing's half-range at 25-400 m, the rise that belts.width's rule needs to give the same edge (s W + c W^2 - R),
    and the analytic belt without the tide."""
    from research.studies.coast_water import belts as B
    k = damping(period_s, specific_yield, transmissivity)
    c = draw / (2.0 * transmissivity)
    if slope is None:                                       # a real profile: its mean rise over the first km
        px_, pz = profile
        near = (px_ >= 0) & (px_ <= 1000.0)
        slope = max(float(np.interp(min(1000.0, px_[-1]), px_, pz)) / min(1000.0, px_[-1]), 1e-3) \
            if near.any() else 0.01
    w0 = float(B.width(slope, c, reach))
    rise = start_rise if start_rise is not None else 0.0
    w_guess = float(B.width(slope, c, reach + max(rise, 4.0 * amplitude)))
    extent = extent or (1.4 * w_guess + 6.0 / k)
    if profile is None:
        x, ground = planar_grid(slope, max(amplitude, 0.05), k, extent)
    else:
        x, ground = profile
    window = window or period_s
    # the belt's slowest mode decays over about 0.4 W^2 Sy / T; three of them before the geometric tail may stop it
    settle = max(4, int(math.ceil(1.2 * w_guess ** 2 * specific_yield / transmissivity / window)))
    rise = rise if start_rise is not None else min(0.5 * k * amplitude / max(slope, 1e-6), 1.0) * amplitude
    year_km = 365.25 * JULIAN_DAY * 1000.0
    out = {}
    for name, a in (('without', 0.0), ('with', amplitude)):
        level = (sea if (sea is not None and a > 0) else sinusoid(a, period_s))
        base = rise if a > 0 else 0.0
        first = base - c * (2.0 * w0 * np.clip(x, 0, w0) - np.clip(x, 0, w0) ** 2)
        r = cross_shore(x, ground, specific_yield, transmissivity, level, draw=draw_t or draw, reach=reach,
                        dt=window / steps, window=window, start=np.minimum(first, np.maximum(ground, 0.0)),
                        tolerance=tolerance, min_windows=min(settle, max_windows), max_windows=max_windows,
                        depth=depth)
        depth_mean = ground - r['mean']
        depth_low = ground - r['low']
        swing = 0.5 * (r['high'] - r['low'])
        edge_mean = edge(x, depth_mean, reach)
        at = lambda d: float(np.interp(d, x, swing))
        within_roots = (ground - r['high']) <= ROOTS_M
        out[name] = dict(edge_mean_m=edge_mean, edge_always_m=edge(x, depth_low, reach),
                         belt_length_m=float(((depth_mean <= reach) & (x >= 0)).sum() * r['dx']),
                         reaches_domain_end=bool(edge_mean >= 0.95 * x[-1]),
                         draw_share_reached_at_high_water=float(r['draw'][within_roots].sum()
                                                                / max(r['draw'].sum(), 1e-30)),
                         draw_m3_yr_per_km=float(r['draw'].sum() * r['dx'] * year_km),
                         swing_half_range_m={f'{d:g}': at(d) for d in (0, 25, 50, 100, 200, 400)},
                         effective_rise_m=float(slope * edge_mean + c * edge_mean ** 2 - reach),
                         exchange_m3_per_m_window=r['inflow'], outflow_m3_per_m_window=r['outflow'],
                         windows=r['windows'], last_change_m=r['last_change'])
    out['analytic_without_m'] = w0
    out['k_per_m'] = k
    out['eps'] = k * amplitude / slope
    out['nodes'] = int(x.size)
    return out
