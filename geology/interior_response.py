"""How today's Moon yields under the Open Moon's water through time, and where its orientation settles.

    python -m geology.interior_response          # results/interior_response.json, results/moon_orientation.json

The present thermal state (thermal_state.py) gives the mantle's viscosity with depth for a central case and two
bracketing ones, a near side and a far side, and single-factor sensitivities. viscoelastic.py turns each profile
into load and tidal Love numbers from a year to a billion years, degrees 2-50 for the central and bracketing
profiles and degree 2 for the rest, and checks them against closed forms, the crustal-loading study's elastic and
thin-shell numbers and the Moon's monthly tide. orientation.py then follows the Moon's axes as the water's load
and the rotational bulge relax: the seas and lakes from the crustal-loading study's rigid geoid change applied at
once, the pores' water as the central uptake delivers it. The deep seas' sag is sized from the same response for
degrees 3-20 against the study's relaxed limit. Bulky arrays go to the research drive.

Inputs read by path: the crustal-loading study's grids (research drive) and results (domain/hydrosphere,
c2c8624a), geography's atlas and drainage products and the groundwater uptake. Their sha256 is recorded.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
import time
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
from scipy.special import lpmv

from geology import orientation as ori
from geology import thermal_state as ts
from geology import viscoelastic as ve
from shared.constants import MOON_RADIUS
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESULTS = HERE / "results"
DRIVE = Path("/media/projectspace/terluna-research/geology/interior-response")
CRUSTAL_NPZ = Path("/media/projectspace/terluna-research/crustal-loading/crustal_loading_nmax359.npz")
# The crustal-loading study reaches this tree at the next integration; until then it is read from domain/hydrosphere.
CRUSTAL_JSON = next((p for p in (ROOT / "geography/results/crustal_loading.json",
                                 Path("/home/kir/Documents/Projects/terluna-hydrosphere/geography/results/crustal_loading.json"))
                     if p.is_file()), ROOT / "geography/results/crustal_loading.json")
ATLAS = ROOT / "geography/products/atlas_28pct_4ppd.npz"
DRAINAGE = ROOT / "geography/products/drainage_28pct_16ppd.npz"
UPTAKE = ROOT / "geography/results/groundwater_uptake.json"
CODE = ("geology/interior_response.py", "geology/thermal_state.py", "geology/viscoelastic.py", "geology/orientation.py",
        "geology/interior_sources.json")
SCHEMA_RESPONSE = "terluna.geology.viscoelastic-response/1"
SCHEMA_ORIENTATION = "terluna.geology.moon-orientation/1"
TIMES = np.r_[0.0, np.logspace(0.0, 9.0, 91)]          # years
REPORT_YEARS = (1.0, 10.0, 100.0, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9)
DEGREES = tuple(range(2, 51))
SAG_DEGREES = tuple(range(3, 21))
LID_TE_KM = (25, 50, 75, 100, 150, 200, 300, 400, 600, 800)   # elastic lids over a fluid interior
FLEXURE_YOUNG = 100e9                                  # the flexure literature's moduli (crustal-loading sources)
FLEXURE_POISSON = 0.25
SPA_LABEL = 1417
MONTH_RATE = 2 * math.pi * 365.25 / 27.321661          # 1/yr
PLACES = {"Jura Island's centre (the CM1 box)": (45.0, -36.1), "Mare Imbrium": (34.7, -14.9),
          "Mare Serenitatis": (27.3, 18.4), "Mare Tranquillitatis": (8.3, 30.8), "Apollo 11": (0.674, 23.473),
          "Copernicus": (9.6, -20.1), "Tycho": (-43.3, -11.2), "Korolev": (-4.2, -157.4),
          "Hertzsprung": (1.4, -128.7), "Chang'e 4": (-45.44, 177.6), "South Pole-Aitken sea's centre": (-53.0, -169.0)}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rel(path):
    try:
        return str(Path(path).relative_to(ROOT))
    except ValueError:
        return str(path)


def moon_of(rheology):
    return ve.Moon([ve.Shell(*shell) for shell in ts.shells(rheology)])


def lid_moon(te_m, young=FLEXURE_YOUNG, poisson=FLEXURE_POISSON):
    """The interior's densities under an elastic lid `te_m` thick with the flexure literature's moduli, over an
    interior viscous enough only to flow: fluid on any time asked about. The lid is compressible; the fluid below is
    incompressible."""
    base = MOON_RADIUS - te_m
    rigidity = young / (2 * (1 + poisson))
    bulk = young / (3 * (1 - 2 * poisson))
    out, bottom = [], 0.0
    for top, density, mu, _ in ts.shells(ts.CENTRAL):
        if mu == 0:
            out.append(ve.Shell(top, density, 0.0))
        elif top <= base + 1.0:
            out.append(ve.Shell(top, density, mu, 1e10))
        else:
            if bottom < base - 1.0:
                out.append(ve.Shell(base, density, mu, 1e10))
            out.append(ve.Shell(top, density, rigidity, math.inf, bulk))
        bottom = top
    return ve.Moon(out)


# ------------------------------------------------------------- the response

def profile_record(rheology):
    depths = np.array([50, 60, 80, 100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1237]) * 1e3
    eta = ts.viscosity_pa_s(depths, rheology)
    t = ts.temperature_k(depths, rheology.profile, rheology.offset_k)
    lid = {f"{y:g}": round(ts.lid_depth_m(y, rheology) / 1e3) for y in (1e3, 1e4, 1e5, 1e6, 1e7)}
    # the depth above which a degree-2 load's viscous time (19/2) eta / (rho g a) exceeds the time
    moon = moon_of(rheology)
    z = np.linspace(47.4e3, 1237e3, 1200)
    viscous = 9.5 * ts.viscosity_pa_s(z, rheology) / moon.stress / ve.YEAR_S
    gravity_lid = {}
    for y in (1e3, 1e4, 1e5, 1e6, 1e7):
        weak = np.flatnonzero(viscous < y)
        gravity_lid[f"{y:g}"] = round(float(z[weak[0]]) / 1e3) if len(weak) else 1237
    return dict(label=rheology.label(), assumptions=asdict(rheology),
                depth_km=(depths / 1e3).round().tolist(), temperature_k=t.round().tolist(),
                log10_viscosity_pa_s=np.round(np.log10(eta), 2).tolist(),
                maxwell_lid_km=lid, degree2_viscous_lid_km=gravity_lid)


def degree_response(moon, n, kinds=("h_load", "k_load", "h", "k")):
    r = ve.response(moon, n, TIMES, kinds)
    return r


def equivalent_te(k_values, lid_curve):
    """Thickness (km) of the elastic lid (E = 100 GPa) over a fluid interior whose relaxed k'_2 equals each value."""
    te = np.array([te for te, _ in lid_curve], float)
    k = np.array([kv for _, kv in lid_curve])
    order = np.argsort(k)
    return np.interp(np.asarray(k_values), k[order], te[order], left=np.nan, right=np.nan)


# ------------------------------------------------------------- loads on the grids

def pore_tensors(lat, lon, height, level, uptake):
    """The pores' water's inertia (3 x 3) by uptake year, from the central bands on sea, lake and land columns."""
    drainage = np.load(DRAINAGE)
    cells = drainage["lake_cell"]
    lake = np.zeros(height.shape, bool)
    lake[(cells // 5760) // 4, (cells % 5760) // 4] = True
    sea = height < level
    lake &= ~sea
    abslat = np.abs(lat)[:, None] * np.ones((1, height.shape[1]))
    out = {}
    bands = uptake["cases"]["central"]["bands"]
    years = list(next(iter(bands.values()))["columns"]["land"]["by_year"].keys())
    for year in years:
        load = np.zeros_like(height)
        for name, band in bands.items():
            lo, hi = (float(x) for x in name.split("-"))
            rows = (abslat >= lo) & ((abslat < hi) if hi < 90 else (abslat <= hi))
            held = {c: band["columns"][c]["by_year"][year]["held_m"] for c in band["columns"]}
            value = np.where(sea, held["sea_floor"], np.where(lake, held.get("lake_floor", held["land"]), held["land"]))
            load[rows] = 1000.0 * value[rows]
        dI, mass = ori.inertia_of_load(load, lat, lon)
        out[float(year)] = (dI, mass)
    return out


def piecewise(tensors):
    """A load history linear between tabulated times, zero at t = 0 and constant after the last."""
    years = np.array([0.0] + sorted(tensors))
    stack = np.array([np.zeros((3, 3))] + [tensors[y] for y in sorted(tensors)])
    def load(t):
        if t >= years[-1]:
            return stack[-1]
        i = np.searchsorted(years, t, side="right") - 1
        w = (t - years[i]) / (years[i + 1] - years[i])
        return (1 - w) * stack[i] + w * stack[i + 1]
    return load


def real_harmonics(field, lat, lon, nmax):
    """Fully normalised real spherical-harmonic coefficients up to nmax, by area-weighted quadrature on the grid,
    and a synthesiser for chosen degrees."""
    weight, _ = ori.grid_geometry(lat, lon)
    colat = np.radians(90.0 - lat)
    x = np.cos(colat)
    phi = np.radians(lon)
    P = {}
    for m in range(nmax + 1):
        for n in range(m, nmax + 1):
            norm = math.sqrt((2 - (m == 0)) * (2 * n + 1) * math.factorial(n - m) / math.factorial(n + m))
            P[n, m] = norm * (-1) ** m * lpmv(m, n, x)
    cos_m = {m: np.cos(m * phi) for m in range(nmax + 1)}
    sin_m = {m: np.sin(m * phi) for m in range(nmax + 1)}
    fw = field * weight
    coef = {}
    for (n, m), p in P.items():
        coef[n, m] = (np.einsum("i,ij,j->", p, fw, cos_m[m]) / (4 * np.pi),
                      np.einsum("i,ij,j->", p, fw, sin_m[m]) / (4 * np.pi))
    def synth(scale):
        out = np.zeros_like(field)
        for (n, m), (c, s) in coef.items():
            if scale.get(n, 0.0) != 0.0:
                out += scale[n] * P[n, m][:, None] * (c * cos_m[m][None] + s * sin_m[m][None])
        return out
    return coef, synth


# ------------------------------------------------------------- main

def main():
    began = time.time()
    DRIVE.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(exist_ok=True)
    sources = json.loads((HERE / "interior_sources.json").read_text())
    ref = sources["crustal_loading_reference"]
    log = lambda *a: print(f"[{time.time() - began:6.0f} s]", *a, flush=True)

    # Validation of the solver against closed forms and the crustal-loading study
    checks = {}
    stack = ts.layers()
    compressible = ve.Moon([ve.Shell(t, d, mu, math.inf, None if mu == 0 else 5 * mu / 3) for t, d, mu in stack])
    incompressible = ve.Moon([ve.Shell(t, d, mu, math.inf) for t, d, mu in stack])
    rows = []
    for i, n in enumerate(ref["elastic_load"]["n"]):
        v = compressible.love(n, np.array([np.inf]))
        rows.append(dict(n=n, h_load=float(v["h_load"][0]), k_load=float(v["k_load"][0]),
                         study_h_load=ref["elastic_load"]["h_load"][i], study_k_load=ref["elastic_load"]["k_load"][i]))
    v2c, v2i = compressible.love(2, np.array([np.inf])), incompressible.love(2, np.array([np.inf]))
    checks["elastic_against_crustal_loading"] = dict(
        load=rows, tidal_compressible=dict(k2=float(v2c["k"][0]), h2=float(v2c["h"][0]), study=ref["tidal"]["compressible"]),
        tidal_incompressible=dict(k2=float(v2i["k"][0]), h2=float(v2i["h"][0]), study=ref["tidal"]["incompressible"]))
    log("elastic checks", rows[0], checks["elastic_against_crustal_loading"]["tidal_incompressible"])

    # Monthly tide of the central model: the low-viscosity zone is the only shell that relaxes within a month
    central_moon = moon_of(ts.CENTRAL)
    monthly = {}
    for lvz in ts.LVZ_VISCOSITY_PA_S:
        m = moon_of(replace(ts.CENTRAL, lvz_pa_s=lvz))
        v = m.love(2, np.array([1j * MONTH_RATE, np.inf]))
        monthly[f"{lvz:g}"] = dict(k2=abs(complex(v["k"][0])), h2=abs(complex(v["h"][0])),
                                   k2_q=abs(complex(v["k"][0])) / abs(complex(v["k"][0]).imag),
                                   elastic_k2=float(v["k"][1].real), elastic_h2=float(v["h"][1].real))
    checks["monthly_tide"] = dict(model=monthly, measured=dict(k2=0.02405, h2=0.0371, k2_park2025=0.024223),
                                  elastic_estimates_garcia2019=dict(k2=0.02294, h2=0.0353))
    log("monthly", monthly)

    # Elastic lids over a fluid interior: the thin-shell comparison and the equivalent-thickness curve
    lid_rows = {}
    for te in sorted(set(LID_TE_KM) | {20, 40}):
        m = lid_moon(te * 1e3)
        per = {}
        for n in (2, 3, 5, 10, 20):
            v = m.love(n, np.array([1e-13]))
            per[n] = dict(h_load=float(v["h_load"][0]), k_load=float(v["k_load"][0]), k=float(v["k"][0]))
        lid_rows[te] = per
    thin = ref["relaxed_thin_shell"]
    checks["lid_against_thin_shell"] = {
        f"te_{te}_km": [dict(n=n, k_load=lid_rows[te][n]["k_load"], h_load=lid_rows[te][n]["h_load"],
                             study_k_load=thin[f"te_{te}_km"]["k_load"][i], study_h_load=thin[f"te_{te}_km"]["h_load"][i])
                        for i, n in enumerate(thin["n"])] for te in (20, 40, 75, 150)}
    checks["secular_tidal_k2_against_keane_matsuyama"] = dict(
        lid_25_km=lid_rows[25][2]["k"], lid_50_km=lid_rows[50][2]["k"], keane_matsuyama=dict(lid_25_km=1.08, lid_50_km=0.86, none=1.44))
    lid_curve = [(te, lid_rows[te][2]["k_load"]) for te in sorted(lid_rows)]
    log("lids", {te: round(k, 4) for te, k in lid_curve})

    # Homogeneous Maxwell sphere: elastic, fluid and the single relaxation time through the fitted spectrum
    mu, rho, eta = 65.6e9, 3344.0, 1e21
    sphere = ve.Moon([ve.Shell(MOON_RADIUS, rho, mu, eta)])
    closed = ve.homogeneous_maxwell(2, mu / sphere.stress, eta / mu / ve.YEAR_S)
    tau = closed["tau_years"]
    r = ve.response(sphere, 2, np.array([0.0, 0.3 * tau, tau, 3 * tau]), ("k_load",))
    exact = closed["k_load"][1] + (closed["k_load"][0] - closed["k_load"][1]) * np.exp(-np.array([0.0, 0.3, 1.0, 3.0]))
    checks["homogeneous_maxwell"] = dict(tau_years=tau, fitted=r["k_load"]["series"].tolist(), exact=exact.tolist(),
                                         elastic=r["k_load"]["elastic"], fluid=closed["k_load"][1])

    # Profiles
    sides = {}
    for off in ts.NEAR_FAR_OFFSET_K:
        near, far = ts.near_far(ts.CENTRAL, off)
        sides[f"near_{off:g}"], sides[f"far_{off:g}"] = near, far
    main_cases = {"central": ts.CENTRAL, "soft": ts.SOFT, "stiff": ts.STIFF}
    sensitivity = {"profile_cool": replace(ts.CENTRAL, profile="cool"), "profile_warm": replace(ts.CENTRAL, profile="warm"),
                   "wet": replace(ts.CENTRAL, wet=True), "grain_1_mm": replace(ts.CENTRAL, grain_m=0.001),
                   "grain_10_mm": replace(ts.CENTRAL, grain_m=0.01), "stress_0.3_mpa": replace(ts.CENTRAL, stress_mpa=0.3),
                   "stress_3_mpa": replace(ts.CENTRAL, stress_mpa=3.0), "lvz_4.5e16": replace(ts.CENTRAL, lvz_pa_s=4.5e16),
                   "soft_near_75": replace(ts.SOFT, offset_k=75.0), "stiff_far_75": replace(ts.STIFF, offset_k=-75.0)}
    profiles = {name: profile_record(rh) for name, rh in {**main_cases, **sides, **sensitivity}.items()}

    # Love numbers through time
    arrays = {"times_years": TIMES, "degrees": np.array(DEGREES), "rates_per_year": ve.RATES_PER_YEAR}
    spectra = {}
    summary = {}
    for name, rh in main_cases.items():
        moon = central_moon if name == "central" else moon_of(rh)
        series = {kind: np.zeros((len(DEGREES), len(TIMES))) for kind in ("h_load", "k_load", "h", "k")}
        misfit = 0.0
        for i, n in enumerate(DEGREES):
            r = degree_response(moon, n)
            for kind in series:
                series[kind][i] = r[kind]["series"]
                misfit = max(misfit, r[kind]["misfit"])
            if n == 2:
                spectra[name] = {kind: dict(elastic=r[kind]["elastic"], amplitudes=r[kind]["amplitudes"].tolist(),
                                            slowest=r[kind]["slowest"]) for kind in ("k_load", "k", "h_load")}
        for kind in series:
            arrays[f"{name}__{kind}"] = series[kind]
        idx = [int(np.argmin(np.abs(TIMES - y))) for y in REPORT_YEARS]
        summary[name] = dict(
            spectrum_misfit=misfit,
            by_degree={str(n): {kind: dict(elastic=float(series[kind][i][0]),
                                           at_years={f"{y:g}": round(float(series[kind][i][j]), 5) for y, j in zip(REPORT_YEARS, idx)})
                                for kind in ("h_load", "k_load")}
                       for i, n in enumerate(DEGREES) if n in (2, 3, 5, 10, 20, 50)},
            tidal_k2={f"{y:g}": round(float(series["k"][0][j]), 4) for y, j in zip(REPORT_YEARS, idx)},
            equivalent_lid_te_km={f"{y:g}": (None if not np.isfinite(e) else round(float(e)))
                                  for y, e in zip(REPORT_YEARS, equivalent_te(series["k_load"][0][idx], lid_curve))})
        log("love", name, "k'_2:", summary[name]["by_degree"]["2"]["k_load"]["at_years"], "misfit", misfit)
    # Stehfest cross-check of the central degree 2
    steh = []
    for y in (10.0, 1e3, 1e5, 1e7):
        pts = ve.stehfest_points(y)
        vals = central_moon.love(2, pts)["k_load"].real
        fitted = float(ve.step_response(y, spectra["central"]["k_load"]["elastic"], np.array(spectra["central"]["k_load"]["amplitudes"]))[0])
        steh.append(dict(years=y, stehfest=ve.stehfest_step(vals), spectrum=fitted))
    checks["stehfest_central_k2"] = steh
    for name, rh in {**sides, **sensitivity}.items():
        r = degree_response(moon_of(rh), 2, ("k_load", "k"))
        spectra[name] = {kind: dict(elastic=r[kind]["elastic"], amplitudes=r[kind]["amplitudes"].tolist(),
                                    slowest=r[kind]["slowest"]) for kind in ("k_load", "k")}
        idx = [int(np.argmin(np.abs(TIMES - y))) for y in REPORT_YEARS]
        summary[name] = dict(k_load_2={f"{y:g}": round(float(r["k_load"]["series"][j]), 5) for y, j in zip(REPORT_YEARS, idx)},
                             tidal_k2={f"{y:g}": round(float(r["k"]["series"][j]), 4) for y, j in zip(REPORT_YEARS, idx)})
    log("sides and sensitivities done")

    # ---------------------------------------------------------------- orientation
    grids = np.load(CRUSTAL_NPZ)
    lat, lon = grids["lat_deg"], grids["lon_deg"]
    tensors = {case: ori.inertia_of_geoid(grids[f"{case}__share__N"].astype(float), lat, lon)
               for case in ("rigid", "elastic", "relaxed_Te_20_km", "relaxed_Te_40_km", "relaxed_Te_75_km", "relaxed_Te_150_km")}
    elastic_state = ori.describe(tensors["elastic"])
    frame = elastic_state["frame"]
    atlas = np.load(ATLAS)
    height = atlas["height_m"].astype(float)
    level = json.loads(atlas["metadata"].item())["sea_level_m"]
    assert np.allclose(atlas["lat_deg"], lat) and np.allclose(atlas["lon_deg"], lon)
    # The South Pole-Aitken sea's water in the rigid case: the sea surface rises with the case's level and N
    rigid_case = json.loads(CRUSTAL_JSON.read_text())["cases"]["rigid"]["share"] if CRUSTAL_JSON.exists() else None
    rigid_level = rigid_case["level_m"] if rigid_case else level + 28.91
    water = 1000.0 * np.clip(rigid_level + grids["rigid__share__N"].astype(float) - height, 0, None)
    water_mass = ori.inertia_of_load(water, lat, lon)[1]
    spa_tensor, spa_mass = ori.inertia_of_load(np.where(atlas["water_label"] == SPA_LABEL, water, 0.0), lat, lon)
    uptake = json.loads(UPTAKE.read_text())
    pores = pore_tensors(lat, lon, height, level, uptake)
    pore_load = piecewise({y: dI for y, (dI, _) in pores.items()})
    seas = tensors["rigid"]
    rest = seas - spa_tensor
    # The water added to hold 28% as the floors sag: the study's relaxed load less the rigid load, by Te
    extra = {te: tensors[f"relaxed_Te_{te}_km"] / (1 + ref["relaxed_thin_shell"][f"te_{te}_km"]["k_load"][0]) - seas
             for te in (20, 40, 75, 150)}

    def love_pair(name, kind):
        sp = spectra[name][kind]
        return sp["elastic"], np.array(sp["amplitudes"])

    def run(loads, tide_name, bulge=True):
        return ori.evolve(TIMES, loads, love_pair(tide_name, "k"), ve.RATES_PER_YEAR, bulge=bulge)

    def compact(rows):
        out = []
        for row in rows:
            if row["years"] == 0 or any(abs(math.log10(row["years"]) - math.log10(y)) < 1e-6 for y in REPORT_YEARS):
                out.append(dict(years=row["years"], pole_deg=round(row["pole_deg"], 3), pole_toward_lon_deg=round(row["pole_toward_lon_deg"], 1),
                                sub_earth_lat_deg=round(row["sub_earth_lat_deg"], 2), sub_earth_lon_deg=round(row["sub_earth_lon_deg"], 2),
                                face_move_deg=round(row["face_move_deg"], 2), b_minus_a_change=round(row["b_minus_a_change"], 4),
                                c_minus_a_change=round(row["c_minus_a_change"], 4), tilt_to_ecliptic_deg=round(row["tilt_to_ecliptic_deg"], 3)))
        return out

    # ---------------------------------------------------------------- the deep seas' sag, degrees 3-20
    seas_load = 1000.0 * np.clip(level - height, 0, None)
    weight, _ = ori.grid_geometry(lat, lon)
    deep = (level - height) > 2000.0
    land = height >= level
    coef, synth = real_harmonics(seas_load, lat, lon, max(SAG_DEGREES))
    mean_density = central_moon.mean_density

    def sag(h_by_degree):
        field = synth({n: h_by_degree[n] * 3.0 / ((2 * n + 1) * mean_density) for n in SAG_DEGREES})
        return float(np.sum(field * weight * deep) / np.sum(weight * deep) - np.sum(field * weight * land) / np.sum(weight * land))

    lid_h = {}
    for te in (20, 40, 75, 150):
        m = lid_moon(te * 1e3)
        lid_h[te] = {n: float(m.love(n, np.array([1e-13]))["h_load"][0]) for n in SAG_DEGREES}
    sag_lids = {te: sag(lid_h[te]) for te in lid_h}
    study_sag = ref["deep_sea_sag"]
    scale = study_sag["relaxed_te_40_m"] / -sag_lids[40]
    sag_series = {}
    for name in main_cases:
        h = arrays[f"{name}__h_load"]
        rows_ = {}
        for y in REPORT_YEARS:
            j = int(np.argmin(np.abs(TIMES - y)))
            rows_[f"{y:g}"] = round(-scale * sag({n: h[DEGREES.index(n)][j] for n in SAG_DEGREES}), 1)
        rows_["elastic"] = round(-scale * sag({n: h[DEGREES.index(n)][0] for n in SAG_DEGREES}), 1)
        sag_series[name] = rows_
    sizing = dict(degrees=list(SAG_DEGREES), scale_to_study=scale,
                  lids_degrees_3_20_m={f"te_{te}_km": round(-sag_lids[te], 1) for te in sag_lids},
                  lids_scaled_m={f"te_{te}_km": round(-scale * sag_lids[te], 1) for te in sag_lids},
                  study_m=study_sag, by_profile_m=sag_series)
    log("sag", sizing["lids_scaled_m"], sag_series["central"])
    log("water", water_mass, "SPA", spa_mass)

    const = lambda tensor: (lambda t: tensor)
    time_series = {}
    for name in main_cases:
        kl = love_pair(name, "k_load")
        time_series[name] = dict(
            seas_and_lakes=compact(run([(const(seas), kl)], name)),
            seas_lakes_and_pores=compact(run([(const(seas), kl), (pore_load, kl)], name)),
            seas_and_lakes_without_bulge=compact(run([(const(seas), kl)], name, bulge=False)))
        log("orientation", name, time_series[name]["seas_and_lakes"][-1])
    # The sea-level feedback, estimated: the extra water that holds 28% grows with the deep floors' sag, taken as the
    # study's relaxed extra water (Te 40 km) times the sag's progress toward the study's 490 m
    for name in main_cases:
        years = np.array([0.0] + list(REPORT_YEARS))
        progress = np.array([sag_series[name]["elastic"]] + [sag_series[name][f"{y:g}"] for y in REPORT_YEARS]) / study_sag["relaxed_te_40_m"]
        def extra_load(t, years=years, progress=progress):
            return extra[40] * np.interp(math.log10(1.0 + t), np.log10(1.0 + years), progress)
        kl = love_pair(name, "k_load")
        time_series[name]["seas_and_lakes_with_sea_level_feedback"] = compact(run([(const(seas), kl), (extra_load, kl)], name))
    for off in ts.NEAR_FAR_OFFSET_K:
        near, far = love_pair(f"near_{off:g}", "k_load"), love_pair(f"far_{off:g}", "k_load")
        time_series[f"spa_far_{off:g}_rest_near"] = dict(
            seas_and_lakes=compact(run([(const(spa_tensor), far), (const(rest), near)], "central")))
    for name in sensitivity:
        time_series[f"sensitivity_{name}"] = dict(seas_and_lakes=compact(run([(const(seas), love_pair(name, "k_load"))], name)))

    # The study's relaxed states against this lane's reading of them: (1 + k'_2) on the rigid load
    feedback = {}
    for te in (20, 40, 75, 150):
        k2 = ref["relaxed_thin_shell"][f"te_{te}_km"]["k_load"][0]
        mine = ori.describe((1 + k2) * seas)
        study = ori.describe(tensors[f"relaxed_Te_{te}_km"])
        feedback[f"te_{te}_km"] = dict(k_load_2=k2, pole_rigid_load_scaled=round(mine["pole_deg"], 3), pole_study=round(study["pole_deg"], 3))
    report_states = {case: {k: (round(v, 4) if isinstance(v, float) else None) for k, v in ori.describe(tensors[case]).items() if k != "frame"}
                     for case in tensors}
    places = {}
    for name, (la, lo) in PLACES.items():
        nl, nlo = ori.to_turned(la, lo, frame)
        places[name] = dict(today=[la, lo], turned=[round(float(nl), 2), round(float(nlo), 2)])

    # ---------------------------------------------------------------- write
    inputs = {rel(p): sha(p) for p in (CRUSTAL_NPZ, ATLAS, DRAINAGE, UPTAKE) if Path(p).exists()}
    if CRUSTAL_JSON.exists():
        inputs["domain/hydrosphere@c2c8624a:geography/results/crustal_loading.json"] = sha(CRUSTAL_JSON)
    producer = dict(domain="geology", files={f: sha(ROOT / f) for f in CODE}, inputs=inputs,
                    constants=constants_used([ROOT / f for f in CODE if f.endswith(".py")]))
    np.savez_compressed(DRIVE / "interior_response_love.npz", **arrays)
    drive_file = DRIVE / "interior_response_love.npz"
    response = dict(
        schema=SCHEMA_RESPONSE, producer=producer,
        evidence=("A calculation on published values: present-day temperature profiles (Kuskov & Kronrod 2009, Karato 2013, "
                  "Khan et al. 2014, digitised), Hirth & Kohlstedt's (2003) olivine flow laws with stated grain size, water and "
                  "stress, the crustal-loading study's interior (Briaud et al. 2023) and low-viscosity zone (Harada et al. 2014). "
                  "Love numbers of a self-gravitating, incompressible Maxwell Moon by the correspondence principle at real s and "
                  "a fitted relaxation spectrum, checked against closed forms, the crustal-loading study's elastic numbers, "
                  "elastic lids over a fluid interior and the monthly tide. Not validated against any observed loading of the "
                  "present Moon; the flow laws are measured on Earth's olivine."),
        reading_rule=("Read 'central' as the response of the present Moon's mid-range thermal state and 'soft' and 'stiff' as its "
                      "bounds; series are the response to a load applied at t = 0 (years). k'_n(t) multiplies a load's own "
                      "potential; the Moon-wide inertia of a load is (1 + k'_2(t)) times its rigid inertia. The tidal k_2(t) is "
                      "the bulge's re-formation per unit forcing. Spectra hold elastic + sum a_j (1 - exp(-r_j t)) on "
                      "rates_per_year. The full degree 2-50 tables are on the research drive."),
        settings=dict(times_years=[0.0, 1.0, 1e9], times_per_decade=10, degrees=[2, 50], fluid_fraction=ve.FLUID_FRACTION,
                      rates_per_year=dict(first=float(ve.RATES_PER_YEAR[0]), last=float(ve.RATES_PER_YEAR[-1]), count=len(ve.RATES_PER_YEAR)),
                      fit_s_per_year=dict(first=float(ve.FIT_S_PER_YEAR[0]), last=float(ve.FIT_S_PER_YEAR[-1]), count=len(ve.FIT_S_PER_YEAR)),
                      incompressible=True, crust="elastic", flexure_young_pa=FLEXURE_YOUNG),
        drive=dict(path=str(drive_file), sha256=sha(drive_file)),
        profiles=profiles, checks=checks, response=summary,
        spectra={name: {kind: dict(elastic=sp["elastic"], slowest=sp["slowest"],
                                   amplitudes=[round(a, 9) for a in sp["amplitudes"]]) for kind, sp in v.items()}
                 for name, v in spectra.items()},
        lid_curve_k_load_2={f"{te}": round(k, 5) for te, k in lid_curve},
        sag=sizing, seconds=round(time.time() - began))
    (RESULTS / "interior_response.json").write_text(json.dumps(response, indent=1) + "\n")
    orientation_product = dict(
        schema=SCHEMA_ORIENTATION, producer=producer,
        evidence=("A calculation: today's moments from lunar laser ranging (beta, gamma; Murphy 2013) and C/MR^2 (Williams et al. "
                  "2014); the water's degree-2 inertia from the crustal-loading study's geoid changes (seas held at 28% with the "
                  "rain-fed lakes) and the groundwater uptake's central bands; the response through time from "
                  "results/interior_response.json. The elastic frame reproduces the orientation report's elastic state."),
        reading_rule=("'elastic_frame' is the orientation set for the second set of full climate simulations: rows of "
                      "rotation_today_to_turned are the turned Moon's axes (toward Earth, the third, north) in today's body frame "
                      "(x toward 0 N 0 E, z north); a point's turned coordinates are R @ n(today). orientation.to_turned and "
                      "to_today convert latitude and longitude. Time series hold the Moon after the seas and lakes arrive at t = 0, "
                      "years after that; 'pole_deg' is the pole's move from today's pole."),
        elastic_frame=dict(**ori.frame_record(frame), pole_move_deg=round(elastic_state["pole_deg"], 3),
                           pole_toward_lon_deg=round(elastic_state["pole_toward_lon_deg"], 2),
                           face_move_deg=round(elastic_state["face_move_deg"], 3), places=places),
        moments=dict(beta=ori.BETA, gamma=ori.GAMMA, c_mr2=ori.C_MR2, radius_m=ori.R_C, A=ori.A_MOMENT, B=ori.B_MOMENT, C=ori.C_MOMENT,
                     fluid_love=ori.fluid_love(ori.C_MR2), q=ori.SPIN_Q),
        load_tensors_kg_m2=dict(**{f"{case}_share": t.tolist() for case, t in tensors.items()},
                                south_pole_aitken_rigid=spa_tensor.tolist(),
                                pores_by_year={f"{y:g}": dI.tolist() for y, (dI, _) in pores.items()}),
        load_masses_kg=dict(rigid_water=water_mass, south_pole_aitken=spa_mass,
                            pores_by_year={f"{y:g}": m for y, (_, m) in pores.items()}),
        report_states=report_states, relaxed_states_rigid_load_scaled=feedback, time_series=time_series,
        seconds=round(time.time() - began))
    (RESULTS / "moon_orientation.json").write_text(json.dumps(orientation_product, indent=1) + "\n")
    log("written")


if __name__ == "__main__":
    sys.exit(main())
