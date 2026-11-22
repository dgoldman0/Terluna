"""A sea's oxygen as one column of levels over its own floor, with the three parts of a canopy colony added.

The column carries the oxygen spent since the water last left the surface (mg/L, the same as g/m3), level by level.
The surface level is held at none, since the air keeps it saturated. Between the levels the water mixes at a
diffusivity given for each level; each level takes the demand of the matter that decays in it; and the colony's parts
enter as terms of their own:

- the feeders take a share of the matter sinking past their tier, so that share of the demand below it never arises
  there, and the colony respires what it caught in the lit layer, where the air pays;
- the trunk sends surface water down and lets it out over a range of depths; the water it displaces rises, each level
  passing its water to the one above;
- the capsules release oxygen in the levels they are given.

A level stands for all the sea's water between two depths: its volume and the open area at its top and bottom follow
the sea's floor, so matter that sinks onto a slope is spent in the level where the slope lies. The sinking flux is
Martin's curve from a reference depth with the depth below it stretched (the south sea study's sinking_stretch): at a
stretch of 1 it is Earth's curve.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

MARTIN_B = 0.86                 # the south sea study's exponent (Martin et al. 1987, recalled there)
MARTIN_FROM_M = 100.0
YEAR_S = 31557600.0
G_O2_PER_MOL = 31.998
O2_PER_C = 1.0                  # mol of oxygen for each mol of carbon respired, as the south sea study counts its export
G_C_PER_MOL = 12.011


@dataclass(frozen=True)
class Column:
    """Levels from the surface down: their edges (m), each level's volume and the open area at each edge, per m2 of the
    sea's surface (so area[0] is 1), the diffusivity at each level's centre (m2/s) and saturation (mg/L)."""
    edge: np.ndarray
    area: np.ndarray
    volume: np.ndarray
    k: np.ndarray
    saturation: float

    @property
    def n(self):
        return len(self.volume)

    @property
    def mid(self):
        return 0.5 * (self.edge[:-1] + self.edge[1:])

    @property
    def dz(self):
        return np.diff(self.edge)

    def with_k(self, k):
        return Column(self.edge, self.area, self.volume, np.asarray(k, float), self.saturation)


def edges_from_centres(mid):
    edge = [0.0]
    for m in mid:
        edge.append(2.0 * m - edge[-1])
    return np.array(edge)


def from_level_volumes(mid, volume, k, saturation):
    """A column from each level's volume: a level's mean area is its volume over its thickness, the open area at an
    edge the mean of the two levels' it parts, the surface the top level's and the last edge closed by the floor."""
    edge = edges_from_centres(mid)
    volume = np.asarray(volume, float)
    n = int(np.flatnonzero(volume > 0).max()) + 1
    edge, volume = edge[:n + 1], volume[:n]
    mean = volume / np.diff(edge)
    area = np.concatenate([[mean[0]], 0.5 * (mean[:-1] + mean[1:]), [0.0]])
    return Column(edge, area / area[0], volume / area[0], np.asarray(k, float)[:n], float(saturation))


def from_areas(mid, depths, areas, k, saturation):
    """A column on the levels `mid` from a sea's area of water deeper than each of `depths`: the open area at each
    edge is read between the table's depths, and a level's volume is its thickness times the mean of its two edges'
    areas. Levels below the sea's deepest water are dropped."""
    edge = edges_from_centres(mid)
    area = np.interp(edge, depths, areas, right=0.0)
    n = int(np.flatnonzero(area[:-1] > 0).max()) + 1
    edge, area = edge[:n + 1], area[:n + 1].copy()
    area[-1] = 0.0
    volume = np.diff(edge) * 0.5 * (area[:-1] + area[1:])
    return Column(edge, area / area[0], volume / area[0], np.asarray(k, float)[:n], float(saturation))


def flux_share(z, stretch=1.0, b=MARTIN_B, z_ref=MARTIN_FROM_M):
    """The share of the export at z_ref still sinking at depth z."""
    below = np.maximum(np.asarray(z, float), z_ref)
    return ((z_ref + stretch * (below - z_ref)) / z_ref) ** -b


def demand(col, export_g_o2, stretch=1.0, catch=0.0, catch_below_m=250.0, z_ref=MARTIN_FROM_M):
    """Each level's demand below z_ref (g O2 per m3 a year) for an export of export_g_o2 (g O2 per m2 of open water at
    z_ref a year): what the sinking flux brings in through the level's top and does not carry out through its bottom,
    spent in its water or on the floor inside it. With `catch` the feeders take that share of the flux that reaches
    catch_below_m, and the flux below is that much smaller. Returns the demand and what the feeders took (g O2 per
    m2 of the sea's surface a year)."""
    top, bottom = col.edge[:-1], col.edge[1:]
    passed = lambda z: np.where(np.asarray(z, float) > catch_below_m, 1.0 - catch, 1.0)    # noqa: E731
    # a level that holds the tier's lower edge is cut at it, so that the flux is taken exactly there
    cut = np.unique(np.concatenate([col.edge, [z_ref, catch_below_m]]))
    cut = cut[(cut >= z_ref) & (cut <= col.edge[-1])]
    a = np.interp(cut, col.edge, col.area)
    carried = export_g_o2 * flux_share(cut, stretch, z_ref=z_ref) * a
    spent = np.zeros(col.n)
    for i in range(len(cut) - 1):
        level = int(np.searchsorted(col.edge, cut[i], side='right')) - 1
        factor = float(passed(0.5 * (cut[i] + cut[i + 1])))
        last = i == len(cut) - 2                                    # the floor closes the column: nothing leaves
        spent[level] += factor * (carried[i] - (0.0 if last else carried[i + 1]))
    a_catch = float(np.interp(catch_below_m, col.edge, col.area)) if catch_below_m < col.edge[-1] else 0.0
    caught = catch * export_g_o2 * float(flux_share(catch_below_m, stretch, z_ref=z_ref)) * a_catch
    return spent / col.volume, caught


def operator(col, trunk_m_yr=0.0, foot_m=(300.0, 500.0)):
    """The matrix L of d(spent)/dt = L spent + demand (per year): mixing between the levels, and the trunk's water,
    let out evenly by thickness over the levels whose centres lie in foot_m, with each level below the surface passing
    its water to the one above. trunk_m_yr is the water sent down each year over each m2 of the sea's surface."""
    n, mid = col.n, col.mid
    L = np.zeros((n, n))
    a, b = col.k[:-1], col.k[1:]
    ki = np.where((a > 0) & (b > 0), 2.0 * a * b / np.maximum(a + b, 1e-300), 0.0) * YEAR_S   # the weaker level limits
    for i in range(1, n):
        g = ki[i - 1] * col.area[i] / (mid[i] - mid[i - 1])          # m3 a year per m2 of surface across edge i
        L[i - 1, i - 1] -= g / col.volume[i - 1]
        L[i - 1, i] += g / col.volume[i - 1]
        L[i, i] -= g / col.volume[i]
        L[i, i - 1] += g / col.volume[i]
    if trunk_m_yr:
        foot = (mid >= foot_m[0]) & (mid < foot_m[1])
        q = np.where(foot, col.dz, 0.0)
        q = trunk_m_yr * q / q.sum()                                 # let out in each level
        rising = np.concatenate([np.cumsum(q[::-1])[::-1][1:], [0.0]])    # through each level's bottom edge
        for i in range(n):
            L[i, i] -= (rising[i] + q[i]) / col.volume[i]            # leaves upward with the level's own water
            if i + 1 < n and rising[i]:
                L[i, i + 1] += rising[i] / col.volume[i]
    L[0, :] = 0.0                                                    # the air keeps the surface
    return L


def integrate(col, net_demand, years, marks=(), L=None, dt=0.05, start=None, yearly=0):
    """The oxygen spent after `years` from `start` (none spent by default), in implicit steps of dt years, never more
    than the water holds. Returns the end state, the states at each of `marks` (years) and the first `yearly` years'
    means."""
    L = operator(col) if L is None else L
    n = col.n
    step = np.linalg.inv(np.eye(n) - dt * L)
    d = np.asarray(net_demand, float).copy()
    d[0] = 0.0
    u = np.zeros(n) if start is None else np.asarray(start, float).copy()
    at, means, acc = {}, [], np.zeros(n)
    per_year = int(round(1.0 / dt))
    for s in range(int(round(years / dt))):
        u = np.clip(step @ (u + dt * d), 0.0, col.saturation)
        if len(means) < yearly:
            acc += u
            if (s + 1) % per_year == 0:
                means.append(acc / per_year)
                acc = np.zeros(n)
        t = (s + 1) * dt
        for m in marks:
            if abs(t - m) < dt / 2:
                at[m] = u.copy()
    return u, at, means


def steady(col, net_demand, L):
    """The oxygen spent once nothing changes any more, where the water never runs out: the balance of L. Where the
    balance asks for more than the water holds the state is not this one, and the caller reads that as failing."""
    A = L.copy()
    d = np.asarray(net_demand, float).copy()
    A[0, :] = 0.0
    A[0, 0] = 1.0
    d[0] = 0.0
    return np.linalg.solve(A, -d)


def calibrated_share(col, run_demand, run_used, upto_m=300.0, fit_m=(90.0, 280.0)):
    """The share of the mean diffusivity above upto_m with which the column follows a tracer run's oxygen spent, year
    by year, on the levels whose centres lie in fit_m. A sea-wide mean of a mixing that comes in bursts and patches
    overstates what the mean water feels, so the share is under 1. Returns the share and the fit's root-mean-square
    error in the logarithm."""
    fit = (col.mid >= fit_m[0]) & (col.mid <= fit_m[1])
    target = np.array([np.asarray(y, float)[:col.n][fit] for y in run_used])
    best = None
    for c in np.geomspace(0.02, 1.0, 80):
        k = col.k.copy()
        k[col.mid < upto_m] *= c
        _, _, means = integrate(col.with_k(k), run_demand, len(run_used), yearly=len(run_used))
        got = np.array([m[fit] for m in means])
        err = float(np.sqrt(np.mean((np.log(got) - np.log(target)) ** 2)))
        if best is None or err < best[1]:
            best = (float(c), err)
    return best


def layer_mean(col, values, span):
    m = (col.mid >= span[0]) & (col.mid <= span[1])
    return float((values[m] * col.volume[m]).sum() / col.volume[m].sum()) if m.any() else None
