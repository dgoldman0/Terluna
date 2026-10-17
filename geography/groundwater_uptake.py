"""How fast the crust takes up the delivered water: heat and water through the cold, porous crust, column by column.

    python -m geography.groundwater_uptake          # results/groundwater_uptake.json, a few seconds

geography.groundwater counts the pore space the crust holds; this module asks how fast it fills. The airless Moon's
crust is cold near the surface and warm below. Its shattered megaregolith, 2-3 km thick under the highlands and 1 km
under the maria, conducts heat so poorly (about 0.2 W/m/K; Warren & Rasmussen 1987) that the endogenic heat flow of
13-17 mW/m2 warms it by 60-85 K a kilometre, from about 255 K below the reach of the month's heat wave at the equator,
and colder toward the poles as the fourth root of the cosine of latitude (Apollo 15 and 17; Diviner). Water that seeps
into rock below 273 K freezes and seals it, as permafrost does, so the crust takes water only where it has warmed.

Each column follows heat and water from the moment the planned climate and the water are in place to ten million
years:
- heat by implicit conduction, with the megaregolith's conductivity rising as air and then water fill its pores and
  its heat capacity rising with the water, the surface held at the water's or ground's temperature and the endogenic
  heat flow entering from below;
- water by Green-Ampt infiltration through permeability falling with depth (Manning & Ingebritsen 1999) at lunar
  gravity and the water's temperature, ponded under the seas and lakes and fed by recharge on land, stopping where
  the rock is still below freezing, with the floors' fine sediment building up as a seal.
The columns stand on the seas' floors, the rain-fed lakes' floors and the land, in six latitude bands, with the seas
over the maria and the lakes and land over the highlands. Moon-wide they give the water the crust holds through time,
how far it would lower the seas, and the top-up that holds them; a warm crust that never freezes the water is run
beside them, and an airless megaregolith of 0.5 W/m/K as a sensitivity.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import solve_banded

from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_RADIUS, MOON_SURFACE_GRAVITY
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = HERE / 'results' / 'groundwater_uptake.json'
SCHEMA = 'terluna.geography.groundwater-uptake/1'
INPUTS = ('geography/results/atlas.json', 'geography/results/drainage.json', 'geography/results/groundwater.json')
IGNORED = ('climate/gcm/products/climatology_A28_dim5_moon.npz', 'geography/products/atlas_28pct_4ppd.npz',
           'geography/products/drainage_28pct_16ppd.npz', 'geography/inputs/ldem_16.img',
           'geography/inputs/ldem_16.lbl', 'geography/inputs/ldem_4.img', 'geography/inputs/ldem_4.lbl',
           'geography/inputs/jggrx_0420a_sha.tab', 'geography/inputs/jggrx_0420a_sha.lbl')
CODE = ('geography/groundwater_uptake.py', 'geography/topography.py')
YEAR_S = JULIAN_DAY * JULIAN_YEAR_DAYS
FREEZING_K = 273.15
RHO_WATER = 1000.0
HEAT_WATER = 4.18e6                      # J/m3/K
HEAT_ROCK = 2.2e6                        # J/m3/K of solid anorthosite, about 2,800 kg/m3 at 0.78 kJ/kg/K
BANDS_DEG = ((0.0, 15.0), (15.0, 30.0), (30.0, 45.0), (45.0, 60.0), (60.0, 75.0), (75.0, 90.0))

# The airless crust (sources in groundwater_sources.json).
AIRLESS_EQUATOR_K = 255.0                # Apollo 15 and 17: 252 and 256 K at 26 and 20 N, below the month's wave
AIRLESS_SPREAD_K = (-5.0, 0.0, 5.0)
HEAT_FLOW_W_M2 = (0.013, 0.015, 0.017)   # regional means, about 25% below Apollo 15's 21 and Apollo 17's 16 mW/m2
K_MEGAREGOLITH_AIRLESS = 0.2             # W/m/K (Warren & Rasmussen 1987)
K_MEGAREGOLITH_GENTLE = 0.5              # the sensitivity: a better-conducting airless megaregolith
K_BEDROCK = 1.4                          # about seven times the megaregolith's (Warren & Rasmussen 1987)
MEGAREGOLITH_M = dict(highland=2500.0, mare=1000.0)   # 2-3 km under the highlands, 1 km under the maria
# The crust once air and water reach it (design guesses: fractured rock with gas, then water, in its cracks).
K_MEGAREGOLITH_AIR = (0.3, 0.5, 0.8)
K_MEGAREGOLITH_WET = (1.0, 1.4, 2.0)
# Permeability: Manning & Ingebritsen's (1999) mean crust, log k = -14 - 3.2 log z (z in km), shifted an order either
# way for the impact-shattered crust, and capped near the surface at volcanic rock's high end (Gleeson et al. 2011).
PERM_1KM_M2 = (1e-15, 1e-14, 1e-13)
PERM_EXPONENT = 3.2
PERM_CAP_M2 = 1e-12
PORE_BOTTOM_M = 20000.0                  # pores close near 20 km (GRAIL's 4% there; Clifford 1993)
HEAT_BOTTOM_M = 30000.0
# The waters' floors and the land's recharge.
SEA_BOTTOM_OFFSET_K = (-2.0, 0.0, 2.0)   # the seas' deep water near their coldest surface, the polar band's
LAKE_BOTTOM_OFFSET_K = (-4.0, -2.0, 0.0)  # a lake's deep water a little below its band's ground temperature
RECHARGE_SHARE = (0.1, 0.2, 0.3)         # of rain, at most rain less evaporation (Earth's mean 0.17; Moeck et al. 2020)
# The floors' fine sediment seals them as it builds up: the seas take their own deposits and the rock that reaches them,
# about 0.1-0.4 mm a year over their floors, and the lakes 0.15-2.4 mm a year of rock and their own accretion (the sea
# chemistry and terrain studies' totals spread over the floors); its permeability runs from lacustrine clay's to
# GLHYMPS's regional clay, 10^-15.8 m2 (Huscroft et al. 2018). Design guesses, slow to fast.
FLOOR_SEDIMENT_MM_YR = dict(sea=(0.36, 0.2, 0.08), lake=(2.4, 0.75, 0.15))
SEAL_PERM_M2 = (1e-18, 1e-17, 1.7e-16)
ARID_RECHARGE_SHARE = 0.01               # dry land still recharges 0.1-5% of its rain (Scanlon et al. 2006)
CASES = ('slow', 'central', 'fast')
STEPS_PER_DECADE = 40
FIRST_YEAR = 0.1
END_YEARS = 1e7
REPORT_YEARS = (10.0, 100.0, 1000.0, 1e4, 1e5, 1e6, 1e7)
EVIDENCE = ('One-dimensional columns of heat and water on committed products and published values: the airless '
            'crust\'s temperatures from Apollo and Diviner with Warren & Rasmussen\'s megaregolith; conduction solved '
            'implicitly with the crust\'s properties changing as air and water fill it; Green-Ampt infiltration '
            'through Manning & Ingebritsen\'s permeability at lunar gravity, the water stopping where the rock is '
            'below freezing and the floors\' sediment sealing them as it builds up. The air and water arrive at once '
            'with the planned climate; the build\'s own warming, delivery heat, lateral flow, the heat the inflowing '
            'water carries, hydrothermal circulation and trapped air are outside it. Earth\'s values are reference '
            'points; design guesses are named.')
READING_RULE = ('Years from the moment the planned climate and the water are in place. Depths in m below the floor or '
                'ground; held water in m of water a square metre of floor or ground; Moon-wide amounts as a global '
                'layer in m, as kg and as a share of the seas\' 1.07e19 kg. The seas\' fall is the held water spread '
                'over the seas\' area without top-up, and the share the seas would keep follows the atlas\'s levels. A '
                'frozen layer is the depth range below freezing; a column is thawed through once no part of it is.')


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rnd(x, sig=6):
    if isinstance(x, dict):
        return {k: rnd(v, sig) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [rnd(v, sig) for v in x]
    if isinstance(x, (np.floating, np.integer)):
        x = x.item()
    if isinstance(x, float) and not math.isfinite(x):
        return None
    if isinstance(x, bool) or not isinstance(x, float) or x == 0.0:
        return x
    return round(x, sig - 1 - int(math.floor(math.log10(abs(x)))))


# --- The column -------------------------------------------------------------------------------------------------------

def grid():
    """Cell edges and centres to HEAT_BOTTOM_M: half a metre at the top, growing 6% a cell to 250 m."""
    dz, z, step = [], 0.0, 0.5
    while z < HEAT_BOTTOM_M:
        dz.append(step)
        z += step
        step = min(step * 1.06, 250.0)
    dz = np.array(dz)
    edges = np.concatenate([[0.0], np.cumsum(dz)])
    return edges, 0.5 * (edges[:-1] + edges[1:]), dz


def viscosity(t_k):
    """Water's viscosity (Pa s), Vogel's fit, held within its range."""
    t = np.clip(t_k, 273.15, 573.0)
    return 2.414e-5 * 10.0 ** (247.8 / (t - 140.0))


def airless_temperature(z, top_k, heat_flow, k_mega, mega_m):
    """The airless crust's steady temperatures: steep through the megaregolith, gentle in the bedrock below."""
    base = top_k + heat_flow * mega_m / k_mega
    return np.where(z <= mega_m, top_k + heat_flow * z / k_mega, base + heat_flow * (z - mega_m) / K_BEDROCK)


def permeability(z, k_1km):
    return np.minimum(k_1km * (np.maximum(z, 1.0) / 1000.0) ** -PERM_EXPONENT, PERM_CAP_M2)


def heat_step(t, dt, k, c, dz, top_k, heat_flow):
    """One implicit step of conduction: Dirichlet temperature at the top, heat flow entering at the bottom."""
    n = t.size
    h = 0.5 * (dz[:-1] + dz[1:])
    k_face = h / (0.5 * dz[:-1] / k[:-1] + 0.5 * dz[1:] / k[1:])          # series conductance between centres
    g_face = k_face / h                                                    # W/m2/K
    g_top = k[0] / (0.5 * dz[0])
    cap = c * dz / dt                                                      # W/m2/K
    diag = cap.copy()
    diag[:-1] += g_face
    diag[1:] += g_face
    diag[0] += g_top
    rhs = cap * t
    rhs[0] += g_top * top_k
    rhs[-1] += heat_flow
    ab = np.zeros((3, n))
    ab[0, 1:] = -g_face
    ab[1] = diag
    ab[2, :-1] = -g_face
    return solve_banded((1, 1), ab, rhs)


def advance(front, dt, z_open, edges, phi, cond, head, supply, seal=0.0):
    """Green-Ampt: move the wetting front over dt (s), from front toward z_open, with the flux (head + depth) over the
    resistance above the front and of the floor's seal (s), or supply when that is smaller."""
    if front >= z_open:
        return front
    j0 = int(np.searchsorted(edges, front, side='right') - 1)
    stops = np.concatenate([[front], edges[j0 + 1:]])
    stops = stops[stops < z_open]
    stops = np.append(stops, z_open)
    seg = np.diff(stops)
    if not seg.size:
        return front
    mids = 0.5 * (stops[:-1] + stops[1:])
    cells = np.searchsorted(edges, mids, side='right') - 1
    resist = np.concatenate([[0.0], np.cumsum(np.diff(edges) / cond)])
    r_mid = seal + resist[cells] + (mids - edges[cells]) / cond[cells]
    flux = (head + mids) / r_mid
    if supply is not None:
        flux = np.minimum(flux, supply)
    tau = phi[cells] * seg / flux
    cum = np.cumsum(tau)
    i = int(np.searchsorted(cum, dt))
    if i >= tau.size:
        return float(stops[-1])
    before = cum[i - 1] if i else 0.0
    return float(stops[i] + seg[i] * (dt - before) / tau[i])


def column(p, top_k, head, supply, setting, warm=False, k_airless=K_MEGAREGOLITH_AIRLESS, sediment_m_yr=0.0):
    """Heat and water in one column through time. supply is the land's recharge (m/s) or None under water, where
    the floor's sediment builds up at sediment_m_yr and seals it."""
    edges, z, dz = grid()
    phi0, d = p['phi0'], p['porosity_e_folding_m']
    phi = phi0 * np.exp(-z / d) * (z < PORE_BOTTOM_M)
    mega = MEGAREGOLITH_M[setting]
    t = airless_temperature(z, p['airless_top_k'], p['heat_flow'], k_airless, mega)
    frozen0 = z[t < FREEZING_K]
    band0 = (float(frozen0.min() - 0.5 * dz[0]) if frozen0.size else 0.0,
             float(frozen0.max()) if frozen0.size else 0.0)
    years = np.concatenate([[0.0], 10 ** np.arange(math.log10(FIRST_YEAR), math.log10(END_YEARS) + 1e-9,
                                                    1.0 / STEPS_PER_DECADE)])
    front, thawed_through = 0.0, None
    rows = {}
    k_perm = permeability(z, p['perm_1km'])
    seal_cond = p['seal_perm'] * RHO_WATER * MOON_SURFACE_GRAVITY / float(viscosity(np.array(top_k)))
    for y0, y1 in zip(years[:-1], years[1:]):
        wet = z < front
        k = np.where(z <= mega, np.where(wet, p['k_wet'], p['k_air']), K_BEDROCK)
        c = (1.0 - phi) * HEAT_ROCK + np.where(wet, phi * HEAT_WATER, 0.0)
        t = heat_step(t, (y1 - y0) * YEAR_S, k, c, dz, top_k, p['heat_flow'])
        cold = z[(t < FREEZING_K) & (z < PORE_BOTTOM_M)]
        if warm or not cold.size:
            z_open = PORE_BOTTOM_M
            if thawed_through is None and not cold.size:
                thawed_through = y1
        else:
            z_open = float(edges[np.searchsorted(z, cold.min())])
        cond = k_perm * RHO_WATER * MOON_SURFACE_GRAVITY / viscosity(t)
        seal = sediment_m_yr * y1 / seal_cond
        front = max(front, advance(front, (y1 - y0) * YEAR_S, z_open, edges, np.maximum(phi, 1e-9), cond, head,
                                   supply, seal))
        held = phi0 * d * (1.0 - math.exp(-front / d))
        for r in REPORT_YEARS:
            if abs(y1 - r) <= 1e-9 * r:
                rows[f'{r:g}'] = dict(front_m=front, held_m=held, seal_m=sediment_m_yr * y1,
                                      frozen_top_m=float(cold.min()) if cold.size else None,
                                      frozen_bottom_m=float(cold.max()) if cold.size else None)
    prev_y, prev = 0.0, 0.0
    for r in REPORT_YEARS:
        row = rows[f'{r:g}']
        row['uptake_m_yr_since_last'] = (row['held_m'] - prev) / (r - prev_y)
        prev_y, prev = r, row['held_m']
    return dict(airless_frozen_layer_m=list(band0), thawed_through_yr=thawed_through, by_year=rows)


# --- The Moon --------------------------------------------------------------------------------------------------------

def bands_setting():
    """Areas, depths and temperatures of the seas, lakes and land by latitude band, from the committed grids."""
    from geography.topography import Grid, height_above_geoid
    from geography.drainage import PPD
    atlas = json.loads((ROOT / INPUTS[0]).read_text())
    clim = dict(np.load(ROOT / IGNORED[0]))
    clim.pop('metadata')
    lat_c = np.abs(clim['lat'])[:, None] * np.ones_like(clim['ts_k'])
    w = np.cos(np.radians(clim['lat']))[:, None] * np.ones_like(clim['ts_k'])
    land_c = clim['lsm'] >= 0.5
    z4 = np.load(ROOT / IGNORED[1])
    g4 = Grid(z4['lat_deg'], z4['lon_deg'], 4)
    a4 = g4.cell_area_m2
    lat4 = np.abs(g4.lat_deg)[:, None] * np.ones_like(a4)
    sea4 = z4['water_label'] > 0
    depth4 = np.maximum(atlas['sea_level_m'] - z4['height_m'], 0.0)
    dr = np.load(ROOT / IGNORED[2])
    cells = dr['lake_cell'].astype(np.int64)
    h16, _, _ = height_above_geoid(PPD, atlas['grid']['geoid_degree'])
    h16 = h16.ravel()
    rows16 = cells // (360 * PPD)
    lat16 = np.abs(90.0 - (rows16 + 0.5) / PPD)
    a16 = (MOON_RADIUS * math.radians(1.0 / PPD)) ** 2 * np.cos(np.radians(90.0 - (rows16 + 0.5) / PPD))
    lake_depth16 = np.maximum(dr['lake_level_m'] - h16[cells], 0.0)

    def within(lat, lo, hi):
        return (lat >= lo) & ((lat < hi) if hi < 90.0 else (lat <= hi))
    out = []
    for lo, hi in BANDS_DEG:
        mc = within(lat_c, lo, hi)
        ml = mc & land_c
        ms = mc & ~land_c
        mean = lambda a, m: float((a * w)[m].sum() / w[m].sum())
        m4 = within(lat4, lo, hi)
        s4 = m4 & sea4
        l16 = within(lat16, lo, hi)
        band_m2 = float(a4[m4].sum())
        sea_m2 = float(a4[s4].sum())
        lake_m2 = float(a16[l16].sum())
        rain = 365.25 * mean(clim['pr_mm_day'], ml) / 1e3
        evap = 365.25 * mean(clim['evap_mm_day'], ml) / 1e3
        out.append(dict(
            band=f'{lo:g}-{hi:g}', area_m2=band_m2, sea_m2=sea_m2, lake_m2=lake_m2,
            land_m2=max(band_m2 - sea_m2 - lake_m2, 0.0),
            sea_depth_m=float((depth4 * a4)[s4].sum() / max(sea_m2, 1.0)),
            lake_depth_m=float((lake_depth16 * a16)[l16].sum() / max(lake_m2, 1.0)) if lake_m2 else 0.0,
            ground_k=mean(clim['ts_k'], ml), sea_surface_k=mean(clim['ts_k'], ms),
            rain_m_yr=rain, evaporation_m_yr=evap, mid_lat=0.5 * (lo + hi)))
    return out, atlas


def params(case, phi0, d, k_airless=K_MEGAREGOLITH_AIRLESS):
    i = CASES.index(case)
    return dict(airless_spread_k=AIRLESS_SPREAD_K[i], heat_flow=HEAT_FLOW_W_M2[i], k_air=K_MEGAREGOLITH_AIR[i],
                k_wet=K_MEGAREGOLITH_WET[i], perm_1km=PERM_1KM_M2[i], sea_bottom_offset_k=SEA_BOTTOM_OFFSET_K[i],
                lake_bottom_offset_k=LAKE_BOTTOM_OFFSET_K[i], recharge_share=RECHARGE_SHARE[i], phi0=phi0,
                porosity_e_folding_m=d, k_airless=k_airless, seal_perm=SEAL_PERM_M2[i],
                sea_sediment_mm_yr=FLOOR_SEDIMENT_MM_YR['sea'][i], lake_sediment_mm_yr=FLOOR_SEDIMENT_MM_YR['lake'][i])


def run_case(case, bands, phi0, d, warm=False, k_airless=K_MEGAREGOLITH_AIRLESS):
    p0 = params(case, phi0, d, k_airless)
    coldest_sea = min(b['sea_surface_k'] for b in bands if b['sea_m2'] > 0)
    out = {}
    for b in bands:
        top = (AIRLESS_EQUATOR_K + p0['airless_spread_k']) * math.cos(math.radians(b['mid_lat'])) ** 0.25
        p = dict(p0, airless_top_k=top)
        recharge = min(p0['recharge_share'] * b['rain_m_yr'], max(b['rain_m_yr'] - b['evaporation_m_yr'], 0.0))
        recharge = max(recharge, ARID_RECHARGE_SHARE * b['rain_m_yr'])
        cols = {}
        if b['sea_m2'] > 0:
            cols['sea_floor'] = column(p, coldest_sea + p0['sea_bottom_offset_k'], b['sea_depth_m'], None, 'mare',
                                       warm, k_airless, p0['sea_sediment_mm_yr'] / 1e3)
        if b['lake_m2'] > 0:
            cols['lake_floor'] = column(p, b['ground_k'] + p0['lake_bottom_offset_k'], b['lake_depth_m'], None,
                                        'highland', warm, k_airless, p0['lake_sediment_mm_yr'] / 1e3)
        cols['land'] = column(p, b['ground_k'], 0.0, recharge / YEAR_S, 'highland', warm, k_airless)
        cols['land']['recharge_mm_yr'] = recharge * 1e3
        out[b['band']] = dict(airless_top_k=top, columns=cols)
    return dict(params=p0, bands=out)


def moon(result, bands, atlas):
    """The crust's held water Moon-wide, the seas' fall without top-up and the top-up that holds them."""
    area = 4.0 * math.pi * MOON_RADIUS ** 2
    sea_m2 = sum(b['sea_m2'] for b in bands)
    seas_kg = atlas['water_mass_kg']
    table = sorted(atlas['share_comparison'], key=lambda r: r['water_mass_kg'])
    masses = np.array([r['water_mass_kg'] for r in table])
    shares = np.array([r['share'] for r in table])
    levels = np.array([r['level_m'] for r in table])
    held = {}
    for y in REPORT_YEARS:
        key = f'{y:g}'
        m3 = 0.0
        for b in bands:
            cols = result['bands'][b['band']]['columns']
            for kind, a in (('sea_floor', b['sea_m2']), ('lake_floor', b['lake_m2']), ('land', b['land_m2'])):
                if kind in cols:
                    m3 += a * cols[kind]['by_year'][key]['held_m']
        kg = m3 * RHO_WATER
        left = seas_kg - kg
        within = left >= masses[0]
        held[key] = dict(global_layer_m=m3 / area, kg=kg, share_of_seas=kg / seas_kg,
                         seas_fall_m=m3 / sea_m2,
                         seas_share_without_top_up=float(np.interp(left, masses, shares)) if within else None,
                         seas_level_without_top_up_m=float(np.interp(left, masses, levels)) if within else None)
    spans = {}
    prev_y, prev = 0.0, 0.0
    for y in REPORT_YEARS:
        kg = held[f'{y:g}']['kg']
        spans[f'{prev_y:g}-{y:g}'] = (kg - prev) / RHO_WATER / 1e9 / (y - prev_y)
        prev_y, prev = y, kg
    return dict(held=held, top_up_km3_yr=spans)


def results():
    phi = json.loads((ROOT / INPUTS[2]).read_text())['porosity_profile']
    phi0, d = phi['phi0'], phi['e_folding_m']
    bands, atlas = bands_setting()
    cases = {}
    for case in CASES:
        r = run_case(case, bands, phi0, d)
        cases[case] = dict(r, moon=moon(r, bands, atlas))
    warm = run_case('central', bands, phi0, d, warm=True)
    gentle = run_case('central', bands, phi0, d, k_airless=K_MEGAREGOLITH_GENTLE)
    return dict(
        schema=SCHEMA, evidence=EVIDENCE, reading_rule=READING_RULE,
        setting=dict(bands=bands, porosity_phi0=phi0, porosity_e_folding_m=d, pore_bottom_m=PORE_BOTTOM_M,
                     seas_kg=atlas['water_mass_kg'], sea_level_m=atlas['sea_level_m']),
        cases=cases,
        warm_crust_central=dict(warm, moon=moon(warm, bands, atlas)),
        gentle_airless_gradient_central=dict(gentle, moon=moon(gentle, bands, atlas)))


def main(argv=None) -> int:
    product = rnd(results())
    product['producer'] = dict(
        domain='geography', files={p: sha256(ROOT / p) for p in CODE},
        inputs={p: sha256(ROOT / p) for p in INPUTS}, ignored_inputs={p: sha256(ROOT / p) for p in IGNORED},
        constants=constants_used([ROOT / p for p in CODE]))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1, allow_nan=False) + '\n')
    for case in CASES:
        h = product['cases'][case]['moon']['held']
        print(case, {y: round(v['global_layer_m'], 1) for y, v in h.items()})
    return 0


if __name__ == '__main__':
    sys.exit(main())
