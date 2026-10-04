"""The monthly tide of the Open Moon's seas.

Earth's pull raises a tide on the Moon about 13 m high at the sub-Earth point.
Its permanent part lies in the atlas's equipotential (the GRAIL geoid with the
static Earth tide). What moves the water is the change through each month:
the Earth's distance varies with the orbit's eccentricity and the evection,
and the sub-Earth point wanders with the optical librations, so the tidal
bulge swells, shrinks and rocks back and forth. The Sun adds a tide of 7 cm.

This product computes that changing tide for every sea of the 28% atlas:

- Earth and Sun over the Moon from lunar_ephemeris.py, checked against
  JPL Horizons, every three hours from 2026 to 2045 (the 18.6-year nodal
  cycle).
- The equilibrium height of the degree-2 and degree-3 tidal potential, with
  the GRAIL and LOLA Love numbers for the solid Moon's own tide at degree 2.
- Each sea's mean level held, since its water cannot leave it, and the
  water's self-attraction (the gravity of the moved water itself), solved
  together for all seas through spherical harmonics to degree 60, on the
  atlas's quarter-degree nodes.
- Constituents: a harmonic analysis of the whole forcing, each line named by
  its combination of the lunar arguments.

tide_dynamics.py checks the equilibrium assumption against a dynamic
shallow-water solution for the leading constituents and gives each sea's
seiche periods. Run `python -m geography.tides` from the repository root;
it writes results/tides.json and the grid product in products/.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from geography import lunar_ephemeris as ephemeris
from shared.constants import (EARTH_GM, GRAVITATIONAL_CONSTANT, MOON_GM, MOON_LOVE_H2, MOON_LOVE_K2,
                              MOON_RADIUS, MOON_SURFACE_GRAVITY, SUN_GM)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
ATLAS_GRID = HERE / "products/atlas_28pct_4ppd.npz"
ATLAS_JSON = HERE / "results/atlas.json"
CHECK = HERE / "tides_ephemeris_check.csv"
WATER = ROOT / "shared/scenarios/waves.json"
PRODUCT = HERE / "products/tides_28pct_4ppd.npz"
RESULT = HERE / "results/tides.json"

STRIDE = 1                  # the atlas's native quarter-degree nodes
DEGREE = 60                 # spherical-harmonic degree of the self-attraction
JD_START = 2461041.5        # 2026-01-01 00:00 TT
YEARS = 19.0                # one 18.6-year nodal cycle and a little more
STEP_DAYS = 0.125
ANOMALISTIC_DAYS = 27.554550

MONOMIALS = ("1", "x", "y", "z", "xx", "yy", "zz", "xy", "xz", "yz",
             "xxx", "yyy", "zzz", "xxy", "xxz", "xyy", "yyz", "xzz", "yzz", "xyz")
DEGREE3 = [i for i, m in enumerate(MONOMIALS) if m in ("x", "y", "z") or len(m) == 3]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def moon_mean_density():
    return 3 * MOON_GM / (4 * math.pi * GRAVITATIONAL_CONSTANT * MOON_RADIUS**3)


# ---------------------------------------------------------------- forcing

def monomial_values(points):
    """Each monomial of MONOMIALS evaluated at unit vectors `points` (..., 3)."""
    x, y, z = points[..., 0], points[..., 1], points[..., 2]
    values = dict(x=x, y=y, z=z)
    out = []
    for name in MONOMIALS:
        out.append(np.ones_like(x) if name == "1" else np.prod([values[c] for c in name], axis=0))
    return np.stack(out, axis=-1)


def monomial_coefficients(direction, coefficient2, coefficient3):
    """Coefficients of the monomials for one body's degree-2 and degree-3 equilibrium height.

    height = c2 [3/2 (e.x)^2 - 1/2] + c3 [5/2 (e.x)^3 - 3/2 (e.x)], with e the
    body's direction and x the unit position, both in the Moon's frame.
    """
    ex, ey, ez = direction[..., 0], direction[..., 1], direction[..., 2]
    c2, c3 = coefficient2, coefficient3
    cube = dict(xxx=ex**3, yyy=ey**3, zzz=ez**3, xxy=3*ex*ex*ey, xxz=3*ex*ex*ez, xyy=3*ex*ey*ey,
                yyz=3*ey*ey*ez, xzz=3*ex*ez*ez, yzz=3*ey*ez*ez, xyz=6*ex*ey*ez)
    terms = {"1": -0.5*c2, "x": -1.5*c3*ex, "y": -1.5*c3*ey, "z": -1.5*c3*ez,
             "xx": 1.5*c2*ex*ex, "yy": 1.5*c2*ey*ey, "zz": 1.5*c2*ez*ez,
             "xy": 3*c2*ex*ey, "xz": 3*c2*ex*ez, "yz": 3*c2*ey*ez}
    terms.update({k: 2.5*c3*v for k, v in cube.items()})
    return np.stack([terms[name] for name in MONOMIALS], axis=-1)


def forcing_coefficients(jd_tt):
    """Monomial coefficients (time, monomial) of the equilibrium tide height, in metres."""
    geometry = ephemeris.earth_and_sun(jd_tt)
    love2 = 1 + MOON_LOVE_K2 - MOON_LOVE_H2
    g, r = MOON_SURFACE_GRAVITY, MOON_RADIUS
    total = 0.0
    for gm, direction, distance in ((EARTH_GM, geometry["earth"], geometry["earth_distance_m"]),
                                    (SUN_GM, geometry["sun"], geometry["sun_distance_m"])):
        c2 = love2 * gm * r**2 / (g * distance**3)
        c3 = gm * r**3 / (g * distance**4)  # degree 3 with the solid Moon held rigid
        total = total + monomial_coefficients(direction, c2, c3)
    return total


# ------------------------------------------------- spherical harmonics

def legendre(nmax, t):
    """Fully normalised associated Legendre functions P[n, m] at t = sin(latitude); 4-pi normalisation."""
    t = np.asarray(t, dtype=float)
    u = np.sqrt(1 - t**2)
    p = np.zeros((nmax + 1, nmax + 1) + t.shape)
    p[0, 0] = 1.0
    if nmax >= 1:
        p[1, 1] = math.sqrt(3) * u
    for m in range(2, nmax + 1):
        p[m, m] = math.sqrt((2*m + 1) / (2*m)) * u * p[m-1, m-1]
    for m in range(0, nmax):
        p[m+1, m] = math.sqrt(2*m + 3) * t * p[m, m]
    for m in range(0, nmax + 1):
        for n in range(m + 2, nmax + 1):
            a = math.sqrt((2*n - 1) * (2*n + 1) / ((n - m) * (n + m)))
            b = math.sqrt((2*n + 1) * (n + m - 1) * (n - m - 1) / ((n - m) * (n + m) * (2*n - 3)))
            p[n, m] = a * t * p[n-1, m] - b * p[n-2, m]
    return p


class Harmonics:
    """Analysis and synthesis on a regular cell-centred latitude-longitude grid."""

    def __init__(self, lat_deg, lon_deg, nmax):
        lat, lon = np.radians(lat_deg), np.radians(lon_deg)
        dlat, dlon = abs(lat[1] - lat[0]), abs(lon[1] - lon[0])
        self.nmax = nmax
        self.p = legendre(nmax, np.sin(lat))                       # (n, m, lat)
        m = np.arange(nmax + 1)
        self.cos, self.sin = np.cos(np.outer(lon, m)), np.sin(np.outer(lon, m))   # (lon, m)
        self.area = (np.sin(lat + dlat/2) - np.sin(lat - dlat/2))[:, None] * dlon   # steradians, (lat, 1)

    def analyse(self, field):
        weighted = field * self.area / (4 * math.pi)
        a = np.einsum("nml,lm->nm", self.p, weighted @ self.cos)
        b = np.einsum("nml,lm->nm", self.p, weighted @ self.sin)
        return a, b

    def synthesise(self, a, b):
        ca = np.einsum("nml,nm->lm", self.p, a)
        cb = np.einsum("nml,nm->lm", self.p, b)
        return ca @ self.cos.T + cb @ self.sin.T


# ------------------------------------------------- the equilibrium response

def load_seas(path=ATLAS_GRID, stride=STRIDE):
    """The atlas at every `stride`-th node: sea labels, depths and the grid."""
    with np.load(path, allow_pickle=False) as source:
        metadata = json.loads(source["metadata"].item())
        if metadata["schema"] != "terluna.geography.atlas-grid/1":
            raise ValueError("Unsupported atlas grid")
        labels = source["water_label"][::stride, ::stride]
        height = source["height_m"][::stride, ::stride]
        lat, lon = source["lat_deg"][::stride], source["lon_deg"][::stride]
    if not (np.allclose(np.diff(lon), lon[1] - lon[0]) and np.allclose(np.diff(lat), lat[1] - lat[0])):
        raise ValueError("Uneven atlas grid")
    return dict(labels=labels, depth=metadata["sea_level_m"] - height, lat=lat, lon=lon,
                level_m=metadata["sea_level_m"], source_sha256=sha256(path))


def sea_bodies(atlas_json=ATLAS_JSON):
    atlas = json.loads(Path(atlas_json).read_text())
    return {b["label"]: b for b in atlas["bodies"]}


class Response:
    """Each sea's level response to a forcing height: its mean held, its own water's attraction included."""

    def __init__(self, grid, bodies, nmax=DEGREE, water_density=None):
        self.grid = grid
        labels = grid["labels"]
        self.bodies = sorted(bodies)
        self.masks = [labels == b for b in self.bodies]
        self.wet = np.any(self.masks, axis=0)
        self.harmonics = Harmonics(grid["lat"], grid["lon"], nmax)
        area = np.broadcast_to(self.harmonics.area, labels.shape)
        self.area = area
        water = water_density or json.loads(WATER.read_text())["water_density_kg_m3"]
        n = np.arange(nmax + 1)
        factor = 3 * water / (moon_mean_density() * (2 * n + 1))
        factor[:2] = 0.0  # degree 0 vanishes with each sea's mass held; degree 1 is a shift of frame
        self.sal_factor = factor
        self.water_density = water

    def hold_means(self, field):
        out = np.where(self.wet, field, 0.0)
        for mask in self.masks:
            out[mask] -= np.sum(field[mask] * self.area[mask]) / np.sum(self.area[mask])
        return out

    def self_attraction(self, height):
        a, b = self.harmonics.analyse(np.where(self.wet, height, 0.0))
        return self.harmonics.synthesise(a * self.sal_factor[:, None], b * self.sal_factor[:, None])

    def solve(self, forcing, attraction=True, tolerance=1e-10, iterations=60):
        height = self.hold_means(forcing)
        if not attraction:
            return height, 0
        for k in range(1, iterations + 1):
            updated = self.hold_means(forcing + self.self_attraction(height))
            change = np.max(np.abs(updated - height)) / max(np.max(np.abs(updated)), 1e-30)
            height = updated
            if change < tolerance:
                return height, k
        raise RuntimeError("Self-attraction iteration did not converge")


def basis_responses(response):
    """The response to each monomial forcing, on every sea node: (monomial, lat, lon)."""
    lat, lon = np.radians(response.grid["lat"]), np.radians(response.grid["lon"])
    la, lo = np.meshgrid(lat, lon, indexing="ij")
    points = np.stack((np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)), axis=-1)
    values = monomial_values(points)
    with_sal, without, iterations = [], [], []
    for k in range(len(MONOMIALS)):
        h, n = response.solve(values[..., k])
        with_sal.append(h)
        without.append(response.solve(values[..., k], attraction=False)[0])
        iterations.append(n)
    return np.array(with_sal), np.array(without), values, iterations


# ------------------------------------------------- constituents

def argument_rates():
    """Rates of D, M, M', F and the node, in cycles per day, from the mean elements."""
    t = np.array([0.0, 1e-4])
    e = ephemeris.mean_elements(t)
    days = 1e-4 * 36525
    return {k: (e[k][1] - e[k][0]) / (2 * math.pi * days) for k in ("D", "M", "Mp", "F", "node")}


def identify(frequency, rates, tolerance):
    """The simplest combination of D, M, M', F and the node within `tolerance` of a frequency (cycles/day).

    Simplicity counts each multiple, the node's twice; the nearest frequency breaks ties.
    """
    names = ("D", "M", "Mp", "F", "node")
    best = None
    for d in range(-4, 5):
        for m in range(-2, 3):
            for mp in range(-3, 4):
                for f in range(-3, 4):
                    for node in (-1, 0, 1):
                        combo = (d, m, mp, f, node)
                        value = sum(c * rates[n] for c, n in zip(combo, names))
                        error = abs(value - frequency)
                        if value <= 0 or error > tolerance:
                            continue
                        key = (abs(d) + abs(m) + abs(mp) + abs(f) + 2 * abs(node), error)
                        if best is None or key < best[0]:
                            best = (key, combo, value)
    if best is None:
        return None
    (_, error), combo, value = best
    label = " ".join(f"{c:+d}{n}".replace("+1", "+").replace("-1", "-")
                     for c, n in zip(combo, ("D", "M", "M'", "F", "N")) if c).lstrip("+")
    return dict(argument=label, multiples=dict(zip(names, combo)), argument_frequency_cpd=value,
                frequency_error_cpd=error)


def constituents(time_days, coefficients, count=12):
    """The strongest spectral lines of the forcing, fitted by least squares, with their arguments."""
    rates = argument_rates()
    centred = coefficients - coefficients.mean(axis=0)
    spectrum = np.fft.rfft(centred * np.hanning(len(time_days))[:, None], axis=0)
    freq = np.fft.rfftfreq(len(time_days), d=time_days[1] - time_days[0])
    power = np.sum(np.abs(spectrum)**2, axis=1)
    lines = []
    band = (freq > 1 / 400) & (freq < 1 / 5)   # from the annual to five-day periods
    candidate = np.flatnonzero(band & (power >= np.roll(power, 1)) & (power >= np.roll(power, -1)))
    candidate = candidate[np.argsort(power[candidate])[::-1]]
    resolution = 1 / (time_days[-1] - time_days[0])
    for index in candidate[:4 * count]:
        named = identify(freq[index], rates, 1.5 * resolution)
        if named is None:
            continue
        exact = named["argument_frequency_cpd"]
        if any(abs(exact - l["frequency_cpd"]) < 1e-9 for l in lines):
            continue
        lines.append(dict(frequency_cpd=exact, period_days=1 / exact, argument=named["argument"],
                          multiples=named["multiples"]))
        if len(lines) == count:
            break
    # Least-squares fit of the mean and every line to every monomial coefficient.
    omega = 2 * math.pi * np.array([l["frequency_cpd"] for l in lines])
    design = np.column_stack([np.ones_like(time_days)] + [np.cos(o * time_days) for o in omega]
                             + [np.sin(o * time_days) for o in omega])
    solution, *_ = np.linalg.lstsq(design, coefficients, rcond=None)
    k = len(lines)
    complex_amplitude = solution[1:k+1] - 1j * solution[k+1:2*k+1]   # a cos(wt) + b sin(wt) = Re[(a - ib) e^{iwt}]
    residual = coefficients - design @ solution
    return lines, solution[0], complex_amplitude, float(np.sqrt(np.mean(residual**2)) / np.sqrt(np.mean(centred**2)))


# ------------------------------------------------- statistics

def node_ranges(coefficients, basis, wet, time_days, chunk=4000):
    """Extreme and median monthly ranges of the water level at every sea node."""
    nodes = basis[:, wet]                                   # (monomial, node)
    window = int(round(ANOMALISTIC_DAYS / (time_days[1] - time_days[0])))
    whole = (len(time_days) // window) * window
    high = np.full(nodes.shape[1], -np.inf)
    low = np.full(nodes.shape[1], np.inf)
    monthly = []
    for start in range(0, whole, window):
        level = coefficients[start:start + window] @ nodes   # (time, node)
        top, bottom = level.max(axis=0), level.min(axis=0)
        monthly.append(top - bottom)
        high, low = np.maximum(high, top), np.minimum(low, bottom)
    monthly = np.array(monthly)
    return dict(extreme=high - low, monthly_median=np.median(monthly, axis=0),
                monthly_p95=np.percentile(monthly, 95, axis=0), months=len(monthly))


def check_ephemeris(path=CHECK):
    import csv
    with open(path, encoding="ascii") as stream:
        records = list(csv.DictReader(line for line in stream if not line.startswith("#")))
    rows = {k: np.array([float(r[k]) for r in records]) for k in records[0] if k != "sample"}
    geometry = ephemeris.earth_and_sun(rows["jd_tt"])
    wrap = lambda d: (d + 180) % 360 - 180
    earth_lon = wrap(geometry["earth_lon_lat_deg"][0] - rows["sub_earth_lon_deg"])
    earth_lat = geometry["earth_lon_lat_deg"][1] - rows["sub_earth_lat_deg"]
    sun_lon = wrap(geometry["sun_lon_lat_deg"][0] - rows["sub_sun_lon_deg"])
    sun_lat = geometry["sun_lon_lat_deg"][1] - rows["sub_sun_lat_deg"]
    distance = geometry["earth_distance_m"] / 1000 - rows["earth_distance_km"]
    return dict(samples=len(records), span_jd_tt=[float(rows["jd_tt"].min()), float(rows["jd_tt"].max())],
                sub_earth_longitude_error_max_deg=float(np.abs(earth_lon).max()),
                sub_earth_latitude_error_max_deg=float(np.abs(earth_lat).max()),
                sub_solar_longitude_error_max_deg=float(np.abs(sun_lon).max()),
                sub_solar_latitude_error_max_deg=float(np.abs(sun_lat).max()),
                earth_distance_error_max_km=float(np.abs(distance).max()),
                earth_distance_error_rms_km=float(np.sqrt(np.mean(distance**2))),
                reading_rule="Differences from JPL Horizons (DE441): the omitted physical librations stay under 0.05 degrees.")


# ------------------------------------------------- the product

STATION_FEATURES = {433: ("Oceanus Procellarum", "Mare Imbrium", "Mare Frigoris", "Mare Serenitatis",
                          "Mare Tranquillitatis", "Mare Fecunditatis", "Mare Nubium", "Mare Humorum"),
                    874: ("Mare Smythii", "Mare Marginis"), 1417: ("Mare Ingenii",)}
HEADLAND = (93.4, 2.34)          # the eastern Smythii headland of the wave shore studies (lon, lat)
DYNAMIC_SEAS = (433, 1417, 874)
DYNAMIC_LINES = 8
LINES = 12


def unit_points(lat_rad, lon_rad):
    return np.stack((np.cos(lat_rad) * np.cos(lon_rad), np.cos(lat_rad) * np.sin(lon_rad), np.sin(lat_rad)), axis=-1)


def nearest(grid, mask, lon, lat):
    rows, cols = np.nonzero(mask)
    la, lo = np.radians(grid["lat"][rows]), np.radians(grid["lon"][cols])
    target = unit_points(math.radians(lat), math.radians(lon))
    k = int(np.argmax(unit_points(la, lo) @ target))
    distance = math.acos(min(1.0, float(unit_points(la[k], lo[k]) @ target))) * MOON_RADIUS / 1000
    return int(rows[k]), int(cols[k]), distance


def weighted_quantile(values, weights, q):
    order = np.argsort(values)
    cumulative = np.cumsum(weights[order])
    return float(values[order][np.searchsorted(cumulative, q * cumulative[-1])])


def station_record(name, body, row, col, offset_km, grid, coefficients, basis, amplitude, lines, time_days):
    series = coefficients @ basis[:, row, col]
    degree2 = coefficients.copy()
    degree2[:, DEGREE3] = 0.0
    series2 = degree2 @ basis[:, row, col]
    window = int(round(ANOMALISTIC_DAYS / STEP_DAYS))
    whole = (len(series) // window) * window
    months = series[:whole].reshape(-1, window)
    response = amplitude @ basis[:, row, col]
    return dict(name=name, body=body, lon_lat=[float(grid["lon"][col]), float(grid["lat"][row])],
                offset_km=offset_km, depth_m=float(grid["depth"][row, col]),
                monthly_range_median_m=float(np.median(np.ptp(months, axis=1))),
                monthly_range_p95_m=float(np.percentile(np.ptp(months, axis=1), 95)),
                extreme_range_m=float(np.ptp(series)), standard_deviation_m=float(series.std()),
                degree3_share_of_extreme_range_m=float(np.ptp(series) - np.ptp(series2)),
                lines=[dict(argument=l["argument"], period_days=l["period_days"], amplitude_m=float(abs(r)),
                            phase_deg=float(np.degrees(np.angle(r))))
                       for l, r in zip(lines, response)])


def dynamic_check(grid_native, amplitude, lines):
    from geography import tide_dynamics as dynamics
    from geography import nomenclature
    named = nomenclature.by_name(nomenclature.features())
    out = {}
    for body in DYNAMIC_SEAS:
        sea = dynamics.SeaGrid(grid_native["labels"], grid_native["depth"], grid_native["lat"], grid_native["lon"], body)
        values = monomial_values(unit_points(sea.lat, sea.lon))
        whole = dynamics.SeaGrid(grid_native["labels"], grid_native["depth"], grid_native["lat"], grid_native["lon"],
                                 body, main_body=False)
        record = dict(main_body_cells=sea.cells, body_cells=sea.body_cells,
                      main_body_area_share=float(sea.area.sum() / whole.area.sum()), lines=[])
        for line, amp in list(zip(lines, amplitude))[:DYNAMIC_LINES]:
            forcing = values @ amp
            result = dynamics.compare(sea, forcing, line["frequency_cpd"])
            result.update(argument=line["argument"], period_days=line["period_days"])
            record["lines"].append(result)
        if body == 433:
            centre = named["Mare Fecunditatis"]
            angle = np.arccos(np.clip(unit_points(sea.lat, sea.lon) @ unit_points(math.radians(centre["lat"]), math.radians(centre["lon"])), -1, 1))
            region = angle * MOON_RADIUS / 1000 < centre["diameter_km"] / 2
            fecund = []
            for line, amp in list(zip(lines, amplitude))[:DYNAMIC_LINES]:
                forcing = values @ amp
                dyn = dynamics.forced_response(sea, forcing, line["frequency_cpd"])
                eq = sea.equilibrium(forcing)
                fecund.append(dict(argument=line["argument"], period_days=line["period_days"],
                                   rms_amplitude_ratio=float(np.sqrt(np.sum(sea.area[region] * np.abs(dyn[region])**2) /
                                                                     np.sum(sea.area[region] * np.abs(eq[region])**2)))))
            # Inside Fecunditatis the shorter lines approach the sub-basin's exchange mode; bed friction sets the response.
            for entry, (line, amp) in zip(fecund, list(zip(lines, amplitude))[:DYNAMIC_LINES]):
                if line["period_days"] > 20:
                    continue
                forcing = values @ amp
                eq = sea.equilibrium(forcing)
                ratios = {}
                for factor in (10, 100):
                    dyn = dynamics.forced_response(sea, forcing, line["frequency_cpd"], friction=factor * dynamics.FRICTION_M_S)
                    ratios[f"friction_x{factor}"] = float(np.sqrt(np.sum(sea.area[region] * np.abs(dyn[region])**2) /
                                                                  np.sum(sea.area[region] * np.abs(eq[region])**2)))
                entry["rms_amplitude_ratio_with_more_friction"] = ratios
            record["mare_fecunditatis"] = dict(cells=int(region.sum()), radius_km=centre["diameter_km"] / 2, lines=fecund)
            # Sensitivity of the strongest fortnightly line to bed friction and rotation.
            fortnightly = max((x for x in zip(lines, amplitude) if x[0]["period_days"] < 20),
                              key=lambda x: np.max(np.abs(values @ x[1])))
            forcing = values @ fortnightly[1]
            record["sensitivity"] = dict(
                argument=fortnightly[0]["argument"],
                friction_x10=dynamics.compare(sea, forcing, fortnightly[0]["frequency_cpd"], friction=10 * dynamics.FRICTION_M_S),
                without_rotation=dynamics.compare(sea, forcing, fortnightly[0]["frequency_cpd"], rotation=False))
        out[str(body)] = record
    return out


def seiches(bodies=DYNAMIC_SEAS, strides=(4, 2, 1), count=5):
    from geography import tide_dynamics as dynamics
    from geography import nomenclature
    named = [f for f in nomenclature.features() if f["type"].split(",")[0] in ("Mare", "Oceanus", "Sinus", "Lacus", "Palus")]
    out = {}
    for body in bodies:
        record = {}
        for stride in strides:
            grid = load_seas(stride=stride)
            sea = dynamics.SeaGrid(grid["labels"], grid["depth"], grid["lat"], grid["lon"], body)
            periods, vectors = dynamics.seiche_periods(sea, count)
            entry = dict(grid_deg=0.25 * stride, cells=sea.cells, periods_hours=[float(p) for p in periods])
            if stride == 1:
                mode = vectors[:, 0] / np.max(np.abs(vectors[:, 0]))
                strong = np.abs(mode) > 0.5
                centre = unit_points(sea.lat[strong], sea.lon[strong]).mean(axis=0)
                centre /= np.linalg.norm(centre)
                feature = max(named, key=lambda f: float(unit_points(math.radians(f["lat"]), math.radians(f["lon"])) @ centre))
                lon = np.degrees(sea.lon[strong]); lon = np.where(lon > 180, lon - 360, lon)
                entry["gravest_mode"] = dict(nearest_feature=feature["name"], cell_share=float(strong.mean()),
                                             lon_range_deg=[float(lon.min()), float(lon.max())],
                                             lat_range_deg=[float(np.degrees(sea.lat[strong]).min()), float(np.degrees(sea.lat[strong]).max())])
            record[f"{0.25 * stride:g}"] = entry
        out[str(body)] = record
    return out


def build():
    ephemeris_check = check_ephemeris()
    grid = load_seas()
    bodies = sea_bodies()
    response = Response(grid, bodies.keys())
    basis, plain, _, iterations = basis_responses(response)
    time_days = np.arange(0.0, YEARS * 365.25, STEP_DAYS)
    coefficients = forcing_coefficients(JD_START + time_days)
    lines, mean, amplitude, residual = constituents(time_days, coefficients, count=LINES)
    wet = response.wet
    ranges = node_ranges(coefficients, basis, wet, time_days)
    ranges_plain = node_ranges(coefficients, plain, wet, time_days)
    labels, area = grid["labels"][wet], response.area[wet]
    seas = {}
    for body in sorted(bodies, key=lambda b: -bodies[b]["area_km2"]):
        m = labels == body
        w = area[m]
        seas[str(body)] = dict(
            names=bodies[body]["water_names"][:4], area_km2=bodies[body]["area_km2"], nodes=int(m.sum()),
            monthly_range_median_m=dict(area_median=weighted_quantile(ranges["monthly_median"][m], w, .5),
                                        area_p95=weighted_quantile(ranges["monthly_median"][m], w, .95),
                                        maximum=float(ranges["monthly_median"][m].max())),
            extreme_range_m=dict(area_median=weighted_quantile(ranges["extreme"][m], w, .5),
                                 maximum=float(ranges["extreme"][m].max())),
            self_attraction_gain=float(np.sum(w * ranges["monthly_median"][m]) / np.sum(w * ranges_plain["monthly_median"][m])))
    from geography import nomenclature
    named = nomenclature.by_name(nomenclature.features())
    stations = []
    for body, names in STATION_FEATURES.items():
        for name in names:
            row, col, offset = nearest(grid, grid["labels"] == body, named[name]["lon"], named[name]["lat"])
            stations.append(station_record(name, body, row, col, offset, grid, coefficients, basis, amplitude, lines, time_days))
    row, col, offset = nearest(grid, grid["labels"] == 874, *HEADLAND)
    stations.append(station_record("Eastern Smythii headland (wave shore studies)", 874, row, col, offset,
                                   grid, coefficients, basis, amplitude, lines, time_days))
    dynamics = dynamic_check(grid, amplitude, lines)
    modes = seiches()
    values = monomial_values(unit_points(*np.meshgrid(np.radians(grid["lat"]), np.radians(grid["lon"]), indexing="ij")))
    line_records = []
    for line, amp in zip(lines, amplitude):
        on_seas = (amp @ basis[:, wet])
        line_records.append(dict(line, forcing_max_m=float(np.max(np.abs(values @ amp))),
                                 sea_response_max_m=float(np.max(np.abs(on_seas)))))
    PRODUCT.parent.mkdir(parents=True, exist_ok=True)
    rows, cols = np.nonzero(wet)
    np.savez_compressed(PRODUCT, lat_deg=grid["lat"], lon_deg=grid["lon"], rows=rows.astype(np.int32),
                        cols=cols.astype(np.int32), labels=labels.astype(np.int32),
                        basis=basis[:, wet].astype(np.float32), monomials=np.array(MONOMIALS),
                        line_period_days=np.array([l["period_days"] for l in lines]),
                        line_amplitude=(amplitude @ basis[:, wet]).astype(np.complex64),
                        monthly_range_median_m=ranges["monthly_median"].astype(np.float32),
                        extreme_range_m=ranges["extreme"].astype(np.float32),
                        metadata=np.array(json.dumps(dict(schema="terluna.geography.tides-grid/1",
                            reading_rule="Water-level response on every sea node of the atlas (rows, cols into lat_deg, lon_deg). The level at time t is coefficients(t) @ basis, with coefficients from tides.forcing_coefficients; line_amplitude gives each constituent's complex amplitude (m), Re[A exp(i w t)] with t in days from 2026-01-01 00:00 TT."))))
    producer = [Path(__file__), HERE / "lunar_ephemeris.py", HERE / "tide_dynamics.py", ROOT / "shared/constants.json"]
    inputs = [ATLAS_GRID, ATLAS_JSON, CHECK, WATER]
    return dict(
        schema="terluna.geography.tides/1",
        producer=dict(domain="geography", files={str(p.relative_to(ROOT)): sha256(p) for p in producer}),
        inputs={str(Path(p).relative_to(ROOT)): sha256(p) for p in inputs},
        product=dict(path=str(PRODUCT.relative_to(ROOT)), sha256=sha256(PRODUCT)),
        evidence=("Equilibrium tide of the 28% atlas's seas from Earth and the Sun, 2026-2045, every three hours: "
                  "degree-2 and degree-3 potential, lunar Love numbers k2 and h2 at degree 2, each sea's volume held, "
                  "and the water's self-attraction to spherical-harmonic degree 60 with the solid Moon rigid under the load. "
                  "A linear shallow-water solution on each main sea's quarter-degree cells checks the equilibrium for the eight "
                  "leading lines."),
        reading_rule=("Heights are metres of water level relative to the 2026-2045 mean, which the atlas's static "
                      "equipotential already holds. The monthly range is the highest minus the lowest level within an "
                      "anomalistic month (27.55 days), its median taken over the 251 months; the extreme range spans the "
                      "whole record. Sea statistics weight nodes by area."),
        settings=dict(grid_deg=0.25 * STRIDE, harmonic_degree=DEGREE, love_numbers=dict(k2=MOON_LOVE_K2, h2=MOON_LOVE_H2),
                      degree3="Rigid solid Moon (Love factor 1)", water_density_kg_m3=response.water_density,
                      moon_mean_density_kg_m3=moon_mean_density(),
                      self_attraction_factor_degree2=float(response.sal_factor[2]),
                      self_attraction_iterations=iterations[1:], span_days=float(time_days[-1]), step_days=STEP_DAYS,
                      start_jd_tt=JD_START),
        ephemeris_check=ephemeris_check,
        constituents=dict(lines=line_records, residual_rms_fraction=residual,
                          reading_rule="Lines fitted together by least squares to the forcing's monomial coefficients; "
                                       "each is named by the simplest combination of the mean elongation D, solar anomaly M, "
                                       "lunar anomaly M', argument of latitude F and node N within the record's resolution."),
        seas=seas, stations=stations, dynamics=dynamics, seiches=modes)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", type=Path, default=RESULT)
    args = parser.parse_args()
    record = build()
    args.output.write_text(json.dumps(record, indent=1) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
