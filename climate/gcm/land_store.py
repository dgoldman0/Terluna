"""The hollows' store in PlaSim's land scheme: stage (b) of the storm-water study (3 November 2026; the record is the
section "The hollows' store in PlaSim's land" in README.md beside this file, the design note
session-tools/agents/storm_water/stage_b_design.md on the research drive).

PlaSim holds each land cell's water in one bucket of wsmax (0.5 m) that sheds water only when full (landmod.f90's
wandr), so on the design years the drylands return all their rain to the air. The storm-water study
(research/studies/storm_water on domain/hydrosphere) worked out offline what the drylands' storms do instead, and its
product results/store_tables.json (schema terluna.research.storm-water-store-tables/1) carries what this scheme reads.
With nlstore 1, in every land cell the product gives a store, each step

1. the slopes shed runoff by infiltration excess: the step's convective rain falls on a wetted share dlswet of the
   cell, and the share of it that runs off comes from the product's table for the chosen surface against the rain's
   depth on that patch and the bucket's relative wetness, the study's Green and Ampt model at lunar gravity over one
   30-minute step; the large-scale rain falls on the whole cell through the same table;
2. the runoff splits by the cell's routing shares (phase 2a): into the hollows' store, to the rivers (PlaSim's own
   drunoff, which roffstep carries), below the roots into the slow store (the channel beds' and slopes' recharge),
   and the rest back into the bucket (the slopes' and beds' soak);
3. the store also takes the rain on its water area, and it spills to the rivers above its capacity;
4. it evaporates over its water area at the pools' own rate (the plan's variant b1, pool_evaporation): Penman's
   evaporation of open water with no heat stored, from the step's sunlight and longwave at the ground, the lowest
   level's air and PlaSim's 10 m wind, through PlaSim's own surface-transfer form over a smooth water surface. The
   evaporation joins the land's (devap), the lowest layer's vapour and the latent heat the cell's one skin loses, so
   the pools share the land's skin temperature;
5. it seeps through its floors at the chosen floor design's rate into a slow store that leaves the model, as the
   drylands' crust takes all their recharge (geography's groundwater rule); the slow store's total is kept.

The store's water area and seepage come from tables at common fills of the cell's hollows (each hollow holding the
same share of its capacity), interpolated linearly in the square root of the fill, which is exact for the paraboloid
bowls and gives nil area at nil volume. Without nlstore every run is as before: every edit only inserts lines, and the
inserted statements run only with nlstore 1. Water is conserved exactly each step: the rain equals the bucket's
change, its own evaporation and runoff, the water sent to the rivers, the store's change, its evaporation and the
slow store's gain (land_step, the tests).

The terms PlaSim keeps from Earth are left as they are and listed (EARTH_TERMS), with the wetted fraction, whose Earth
value of 0.3 the scheme carries as a parameter without settling it. LAND_STORE_PATCHES put the scheme into landmod.f90
alone; its year-file fields go out through PlaSim's own 2-d diagnostic slots (CODES), so no file another run's key
hashes changes. The runner takes it with --land-store (exoplasim_run.py).
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np

from shared.constants import EXOPLASIM_FREE_CONVECTION, STANDARD_ATMOSPHERE, STANDARD_GRAVITY

LAND_STORE_MARKER = "Terluna: the hollows' store"
PRODUCT_SCHEMA = 'terluna.research.storm-water-store-tables/1'
TABLE_VERSION = 1                 # landstore.dat's layout, checked by lstoreini
TABLE_NAME = 'landstore.dat'      # the run's copy of the tables, in PlaSim's working folder
FLOORS = ('clay_floored', 'sandy_banked', 'tight_seal')            # the floor designs, nlsfloor 1-3 (phase 2b's)
SURFACES = ('sandy_loam', 'loam', 'sealed_loam', 'biocrust')        # the slopes' surfaces, nlssurf 1-4 (phase 1's)
LIBRARIES = ('walnut_gulch', 'cm1_6km')                             # the storm libraries the routing shares follow
HOLLOWS = ('young', 'mature')                                       # the hollows' age (phase 2a's states)
SHARES = ('hollows', 'rivers', 'recharge')                          # of the slopes' runoff; the rest soaks back
NLON, NLAT = 64, 32                                                 # T21: the product's cells are j * NLON + i
NUGP = NLON * NLAT
NLSF, NLSD, NLSW = 12, 19, 11     # fills, rain depths and wetnesses of the product's grids (store_tables.py)
NLSOUT = 7
# The pools' surface and air (b1). Water's albedo and emissivity are phase 2b's (a property and the plan's reading);
# the roughness is that of calm water, about COARE 3.0's cap on the roughness for heat and vapour (1.15e-4 m) and
# Charnock's at lunar gravity for these winds (5e-5 to 1e-4 m), a design choice; the wind is PlaSim's 10 m wind
# (WIND10_PATCHES), at least CM1's least wind of 0.1 m/s as the gust patch has it; Louis's constants are PlaSim's
# vdiff_b, vdiff_c and vdiff_d (fluxmod.f90), Earth calibrations; the free convection is Miller, Beljaars and
# Palmer's tank law scaled to the planet as exoplasim_run.free_convection_setting scales it (pool_free_convection).
POOL = dict(albedo=0.07, emissivity=0.97, roughness_m=1.0e-4, wind_height_m=10.0, least_wind_m_s=0.1,
            louis_b=5.0, louis_c=5.0, louis_d=5.0)
STEFAN_BOLTZMANN_PLASIM = 5.67e-8   # radmod.f90's SBK, the radiation's own constant
KARMAN = 0.4                        # fluxmod.f90's vonkarman
PENMAN_ITERATIONS = 48              # halvings of the bracket: the water's temperature to 3e-13 K
PENMAN_RANGE_K = 40.0               # the bracket about the air's temperature (a guard that the drylands never reach)
WETTED_FRACTION_EARTH = 0.3         # the Met Office scheme's share of a cell under convective rain (JULES, Best et
                                    # al. 2011, recalled via the plan): an Earth-tuned term, carried and not settled
# The year-file fields: PlaSim's 2-d diagnostic slots 1-7 (dgp2d, which outdiag writes as codes 1-7 and nothing else
# fills), each the mean over the output interval; the runner sets ndiaggp2d to NLSOUT and names the codes.
CODES = {1: ('mrros', 'surface_runoff_by_infiltration_excess_from_the_slopes', 'm s-1'),
         2: ('lstore', 'water_in_the_hollows_store_over_the_cell', 'm'),
         3: ('lstorea', 'water_area_share_of_the_hollows_store', '1'),
         4: ('lstoree', 'evaporation_from_the_hollows_store_downward', 'm s-1'),
         5: ('lstores', 'seepage_from_the_hollows_store_into_the_slow_store', 'm s-1'),
         6: ('lstoreg', 'recharge_below_the_roots_from_the_slopes_and_channel_beds', 'm s-1'),
         7: ('lstorer', 'water_the_store_scheme_sends_to_the_rivers', 'm s-1')}
RUNOFF_NAME = 'mrro'                # PlaSim's own code 160: the local runoff with the rivers' convergence, as mkradv
                                    # leaves drunoff (the river store's tendency over land, the rivers' inflow at sea)
# PlaSim's own Earth-tuned terms the build leaves as they are (landmod.f90), with what the study's evidence says.
EARTH_TERMS = {
    'wsmax': dict(value_m=0.5, where='landmod.f90 WSMAX_EARTH, unchanged by the design run',
                  evidence='the drylands keep the bucket 0.4-8% full on the design years; the study reads it as the '
                           'plant-available water between the wilting point and field capacity of a soil of '
                           'unspecified depth, and its capacity sets how much rain the slopes soak before evaporating'),
    'drhsfull': dict(value=0.4, where='landmod.f90 drhsfull: the wetness factor reaches 1 at 0.4 wsmax',
                     evidence='an Earth calibration of the bucket\'s evaporation; the drylands\' factor is 0.05-0.08, '
                              'which with the 28 mm/day daytime potential sets their evaporation'),
    'wetted_fraction': dict(value=WETTED_FRACTION_EARTH, where='dlswet, this scheme',
                            evidence='Earth\'s schemes spread a step\'s convective rain over about 0.3 of a cell; '
                                     'the lunar storms\' cover at T21 is unmeasured until a dryland CM1 box runs, and '
                                     'the share decides how intense PlaSim\'s cell-mean rain is on the patch'),
    'green_ampt': dict(where='the product\'s runoff tables (infiltration.py)',
                       evidence='Rawls et al. (1983) class means and Rawls et al. (1990) seal at lunar gravity; the '
                                'surfaces are design choices'),
    'louis': dict(values=[5.0, 5.0, 5.0], where='POOL louis_b, louis_c, louis_d',
                  evidence='PlaSim\'s own Earth-calibrated stability factors, here over the pools'),
}


def pool_constants(c: dict, gravity_m_s2: float, pressure_pa: float, pool: dict = POOL) -> np.ndarray:
    """The constants lspool reads, in its order: PlaSim's gravity, gas constant, Rd/Rv, heat capacity of dry air,
    latent heat, the saturation formula's ra1, ra2, ra4 and tmelt, Stefan and Boltzmann's constant, and the pool's
    albedo, emissivity, roughness, wind height, least wind, free-convection constant and Louis's three constants."""
    return np.array([c['ga'], c['gascon'], c['rdbrv'], c['acpd'], c['alv'], c['ra1'], c['ra2'], c['ra4'],
                     c['tmelt'], STEFAN_BOLTZMANN_PLASIM, pool['albedo'], pool['emissivity'], pool['roughness_m'],
                     pool['wind_height_m'], pool['least_wind_m_s'], pool_free_convection(gravity_m_s2, pressure_pa),
                     pool['louis_b'], pool['louis_c'], pool['louis_d']], float)


def pool_free_convection(gravity_m_s2: float, pressure_pa: float) -> float:
    """Miller, Beljaars and Palmer's (1992) still-air constant over water, 0.0016 m s-1 K-1/3 at Earth's gravity and
    1013 hPa, carried to the planet as (g / p_s)^(1/3) (exoplasim_run.free_convection_setting's scaling, which phase
    2b's pools take): 0.00083 at lunar gravity and 1.2 atm."""
    return EXOPLASIM_FREE_CONVECTION * (gravity_m_s2 / STANDARD_GRAVITY * STANDARD_ATMOSPHERE / pressure_pa) ** (1.0 / 3.0)


PLASIM_CONSTANTS = dict(alv=2.5008e6, ra1=610.78, ra2=17.2693882, ra4=35.86, tmelt=273.16)   # pumamod's and p_exo's


def saturation_humidity(t, p, pc):
    """PlaSim's saturation humidity (landstep's dqs): rdbrv ra1 exp(ra2 (T - tmelt) / (T - ra4)) / p, over
    1 - (1/rdbrv - 1) of itself."""
    q = pc[2] * pc[5] * np.exp(pc[6] * (t - pc[8]) / (t - pc[7])) / p
    return q / (1.0 - (1.0 / pc[2] - 1.0) * q)


def pool_transfer(tw, ta, qa, p, s, pc):
    """The pool's transfer velocity (m/s) at water temperature tw: PlaSim's surface form (surflx) over smooth water,
    the wind s at the wind's height zw times the neutral coefficient (kappa / ln(zw / z0 + 1))^2 times Louis's factor
    of the bulk Richardson number between the air and the water's virtual temperatures, and at least the free
    convection off a pool warmer than its air."""
    epsv = 1.0 / pc[2] - 1.0
    cn = (KARMAN / np.log(pc[13] / pc[12] + 1.0)) ** 2
    qw = saturation_humidity(tw, p, pc)
    thva = ta * (1.0 + epsv * qa)
    thvw = tw * (1.0 + epsv * qw)
    ri = pc[0] * pc[13] * (thva - thvw) / (s * s * thva)
    unstable = 1.0 - 3.0 * pc[16] * ri / (1.0 + 3.0 * pc[16] * pc[17] * cn * np.sqrt(np.maximum(-ri, 0.0) * (pc[13] / pc[12] + 1.0)))
    stable = 1.0 / (1.0 + 3.0 * pc[16] * np.maximum(ri, 0.0) * np.sqrt(1.0 + pc[18] * np.maximum(ri, 0.0)))
    fh = np.where(ri <= 0.0, unstable, stable)
    return np.maximum(s * cn * fh, pc[15] * np.maximum(thvw - thva, 0.0) ** (1.0 / 3.0))


def pool_residual(tw, v, sw, lw, ta, qa, p, rho, pc):
    """The water's net radiation less its latent and sensible heat (W/m2) at tw with the transfer v."""
    rn = (1.0 - pc[10]) * sw + pc[11] * (lw - pc[9] * tw ** 4)
    return rn - rho * pc[4] * v * (saturation_humidity(tw, p, pc) - qa) - rho * pc[3] * v * (tw - ta)


def pool_evaporation(sw, lw, ta, qa, p, wind, pc, iterations=PENMAN_ITERATIONS):
    """Penman's evaporation of open water in its exact form, with no heat stored (lspool step for step): the water
    temperature at which the water's net radiation balances its latent and sensible heat with the transfer at that
    temperature, found by halving a bracket of PENMAN_RANGE_K about the air's temperature (the free convection's cube
    root gives the transfer a cusp where the water meets the air's virtual temperature, which defeats Newton's method).
    sw and lw are the sunlight and longwave down at the ground (W/m2), ta the air's temperature and qa its humidity,
    p the pressure (Pa), wind the wind at the pool's wind height (m/s), pc pool_constants. Returns the evaporation
    (m/s of water, upward positive), the water's temperature (K) and the transfer velocity (m/s).

    Phase 2b's pools evaporate by an energy balance with heat stored in the pool and its floor and COARE 3.0's transfer
    at lunar gravity on the wind with the boundary layer's gust; its no-storage limit (phase2b.penman) is this balance
    with its own transfer, and it agrees with the shallow pools' evaporation through the lunar day within a few per
    cent and with the 0.3 m pool's within about a tenth (stage_b_report.md)."""
    sw, lw, ta, qa, p, wind = np.broadcast_arrays(*(np.asarray(v, float) for v in (sw, lw, ta, qa, p, wind)))
    epsv = 1.0 / pc[2] - 1.0
    s = np.maximum(wind, pc[14])
    rho = p / (pc[1] * ta * (1.0 + epsv * qa))
    lo, hi = ta - PENMAN_RANGE_K, ta + PENMAN_RANGE_K
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        up = pool_residual(mid, pool_transfer(mid, ta, qa, p, s, pc), sw, lw, ta, qa, p, rho, pc) > 0.0
        lo, hi = np.where(up, mid, lo), np.where(up, hi, mid)
    tw = 0.5 * (lo + hi)
    v = pool_transfer(tw, ta, qa, p, s, pc)
    return rho * v * (saturation_humidity(tw, p, pc) - qa) / 1000.0, tw, v


def radiation_from_plasim(sw_net, albedo, lw_net, ts, sbk=STEFAN_BOLTZMANN_PLASIM):
    """Sunlight and longwave down at the ground from PlaSim's net fluxes there (dswfl and dlwfl at NLEP, the year
    files' rss and rls), its albedo and the ground's temperature, the ground emitting as a black body as PlaSim's
    radiation has it: (sw_net / (1 - albedo), lw_net + sbk ts^4)."""
    return (np.asarray(sw_net, float) / np.maximum(1.0 - np.asarray(albedo, float), 0.05),
            np.asarray(lw_net, float) + sbk * np.asarray(ts, float) ** 4)


# ---------------------------------------------------------------------------------------------------------------------
# The store's tables and one step.

def curve(fill, area, seep, cap, v):
    """The store's water area (share of the cell) and seepage (m/s over the cell) at volume v (m), from its tables at
    the common fills (lscurve): linear in the square root of the fill between the table's fills, nil at nil volume and
    held at the full store's values above its capacity. fill (NLSF,), area and seep (cell, NLSF), cap and v (cell,)."""
    cap, v = np.asarray(cap, float), np.asarray(v, float)
    on = (cap > 0.0) & (v > 0.0)
    f = np.where(on, np.minimum(v / np.where(cap > 0.0, cap, 1.0), 1.0), 0.0)
    zs = np.sqrt(f)
    root = np.sqrt(np.asarray(fill, float))
    k = np.clip(np.searchsorted(root, zs, side='left'), 1, root.size - 1)    # first fill whose root reaches zs
    w = (zs - root[k - 1]) / np.maximum(root[k] - root[k - 1], 1e-30)
    rows = np.arange(cap.size)
    fw = area[rows, k - 1] + w * (area[rows, k] - area[rows, k - 1])
    q = seep[rows, k - 1] + w * (seep[rows, k] - seep[rows, k - 1])
    return np.where(on, fw, 0.0), np.where(on, q, 0.0)


def runoff_share(depths, wets, table, d_mm, w):
    """The share of a step's rain that runs off (lsshare): bilinear in the rain's depth on its patch (mm) and the
    bucket's relative wetness, both held within the table's grids. table (NLSD, NLSW) for one surface."""
    d = np.clip(np.asarray(d_mm, float), depths[0], depths[-1])
    w = np.clip(np.asarray(w, float), wets[0], wets[-1])
    i = np.clip(np.searchsorted(depths, d, side='left') - 1, 0, depths.size - 2)    # depths[i] < d <= depths[i+1]
    j = np.clip(np.searchsorted(wets, w, side='left') - 1, 0, wets.size - 2)
    a = (d - depths[i]) / (depths[i + 1] - depths[i])
    b = (w - wets[j]) / (wets[j + 1] - wets[j])
    return ((1.0 - a) * (1.0 - b) * table[i, j] + a * (1.0 - b) * table[i + 1, j] + (1.0 - a) * b * table[i, j + 1]
            + a * b * table[i + 1, j + 1])


def store_step(dt, prc, prl, es, wetted, w_bucket, wmax, fw, v, cap, fill, area, seep, shares, depths, wets, table):
    """One step of the store in each cell (lsstep), after its evaporation es (m/s, upward positive, at most v / dt
    where positive) has joined the land's: the rain on its water area, the slopes' runoff by infiltration excess on the
    rest, its routing, the seepage at the step's start, the spill above capacity. prc and prl are the convective and
    large-scale rain (m/s), wetted the convective rain's wetted share, w_bucket and wmax the bucket's water and
    capacity (m), fw and v the store's area share and water at the step's start, shares (cell, 3) the hollows', rivers'
    and recharge's shares of the runoff. Returns the change to PlaSim's dwater (m/s) and the step's quantities."""
    rel = np.where(wmax > 0.0, np.minimum(np.maximum(w_bucket / np.where(wmax > 0.0, wmax, 1.0), 0.0), 1.0), 0.0)
    ps = fw * (prc + prl)
    dc = prc * dt * 1000.0 / max(wetted, 1e-6)
    dl = prl * dt * 1000.0
    roff = (runoff_share(depths, wets, table, dc, rel) * prc + runoff_share(depths, wets, table, dl, rel) * prl) * (1.0 - fw)
    inflow = ps + shares[:, 0] * roff
    _, q0 = curve(fill, area, seep, cap, v)
    vv = v - es * dt + inflow * dt
    seepage = np.minimum(q0, np.maximum(vv, 0.0) / dt)
    vv = vv - seepage * dt
    spill = np.maximum(vv - cap, 0.0) / dt
    vn = np.maximum(vv - spill * dt, 0.0)
    fwn, _ = curve(fill, area, seep, cap, vn)
    rivers = shares[:, 1] * roff
    recharge = shares[:, 2] * roff
    dwater = es - ps - (shares[:, 0] + shares[:, 1] + shares[:, 2]) * roff
    return dict(dwater=dwater, v=vn, fw=fwn, runoff=roff, seepage=seepage, recharge=recharge, spill=spill,
                rivers=rivers, inflow=inflow, rain_on_store=ps)


def store_evaporation(fw, v, pool_rate, dt, snow=None):
    """The store's evaporation this step (lstoreevap): its water area times the pools' rate (m/s, upward positive), at
    most its water where it evaporates; none under snow. Dew onto the store (a negative rate) is not limited."""
    es = np.asarray(fw, float) * np.asarray(pool_rate, float)
    es = np.where(es > 0.0, np.minimum(es, np.asarray(v, float) / dt), es)
    if snow is not None:
        es = np.where(np.asarray(snow, float) > 0.0, 0.0, es)
    return es


def land_step(dt, state, forcing, tables, wetted):
    """One land step of a store cell with PlaSim's bucket as soilstep and wandr have it (no snow): the store's
    evaporation, tands' dwater with the land's own evaporation devap_bucket (m/s, PlaSim's sign), the store's step and
    the bucket's. Returns the new state, the fluxes and the water budget's terms (m over the cell), which close
    exactly: rain = bucket change + bucket's evaporation + runoff to the rivers + store change + store's evaporation
    + slow store's gain."""
    w, v, fw = state['w'], state['v'], state['fw']
    es = store_evaporation(fw, v, forcing['pool_rate'], dt)
    devap = forcing['devap_bucket'] - es
    dwater = devap + forcing['prl'] + forcing['prc']
    s = store_step(dt, forcing['prc'], forcing['prl'], es, wetted, w, tables['wmax'], fw, v, tables['cap'],
                   tables['fill'], tables['area'], tables['seep'], tables['shares'], tables['depths'], tables['wets'],
                   tables['table'])
    dwater = dwater + s['dwater']
    wn = w + dt * dwater
    bucket_runoff = np.maximum(0.0, wn - tables['wmax']) / dt
    clamp = np.maximum(-wn, 0.0)                                   # water PlaSim's clamp at nil would create
    wn = np.maximum(np.minimum(tables['wmax'], wn), 0.0)
    drunoff = bucket_runoff + s['rivers'] + s['spill']
    budget = dict(rain=(forcing['prc'] + forcing['prl']) * dt, bucket_change=wn - w,
                  bucket_evaporation=-forcing['devap_bucket'] * dt, to_rivers=drunoff * dt, store_change=s['v'] - v,
                  store_evaporation=es * dt, slow_store=(s['seepage'] + s['recharge']) * dt, clamp=clamp)
    return dict(w=wn, v=s['v'], fw=s['fw']), dict(s, es=es, devap=devap, drunoff=drunoff), budget


# ---------------------------------------------------------------------------------------------------------------------
# The product and the run's input file.

def read_product(path) -> dict:
    """The storm-water study's store tables (store_tables.py), refused unless they carry the schema, the grids this
    patch reads and every floor, surface, library and hollow age."""
    path = Path(path)
    product = json.loads(path.read_text())
    if product.get('schema') != PRODUCT_SCHEMA:
        raise ValueError(f'{path} is not a {PRODUCT_SCHEMA} product')
    g = product['grids']
    if (len(g['fill']) != NLSF or len(g['rain_depth_mm']) != NLSD or len(g['wetness']) != NLSW
            or g['fill'][0] != 0.0 or g['fill'][-1] != 1.0 or g['rain_depth_mm'][0] != 0.0):
        raise ValueError(f'{path}: its grids are not the {NLSF} fills, {NLSD} depths and {NLSW} wetnesses this patch reads')
    if (tuple(product['floors']) != FLOORS or tuple(product['surfaces']) != SURFACES
            or tuple(product['libraries']) != LIBRARIES or tuple(product['hollows']) != HOLLOWS):
        raise ValueError(f'{path}: its floors, surfaces, libraries or hollow ages are not the patch\'s')
    if product['grid'] != dict(nlon=NLON, nlat=NLAT):
        raise ValueError(f'{path}: its cells are not on PlaSim\'s T21 grid')
    return product


def table_rows(product: dict, library: str, hollows: str) -> tuple:
    """The input file's tables: the grids, each surface's runoff shares (NLSD x NLSW), and per T21 point in PlaSim's
    order (north to south, from 0 E) the capacity (m), the water area at each fill, the seepage at each fill under each
    floor (m/s) and the three routing shares under each surface; nothing where the product has no store."""
    if library not in LIBRARIES or hollows not in HOLLOWS:
        raise ValueError(f'the library is one of {LIBRARIES} and the hollows one of {HOLLOWS}')
    g = product['grids']
    runoff = np.array([product['runoff_tables'][s] for s in SURFACES], float)        # (surface, depth, wetness)
    rows = np.zeros((NUGP, 1 + NLSF * (1 + len(FLOORS)) + len(SHARES) * len(SURFACES)))
    for c in product['cells']:
        k = int(c['j']) * NLON + int(c['i'])
        st = c['store'][hollows]
        if st['capacity_m'] <= 0.0:
            continue
        shares = [c['shares'][library][s][hollows][q] for s in SURFACES for q in SHARES]
        rows[k] = [st['capacity_m'], *st['area_share'], *[x for f in FLOORS for x in st['seepage_m_s'][f]], *shares]
    return np.array(g['fill'], float), np.array(g['rain_depth_mm'], float), np.array(g['wetness'], float), runoff, rows


def write_table(product: dict, library: str, hollows: str, path) -> Path:
    """landstore.dat: a header (layout, points, fills, depths, wetnesses, floors, surfaces), the grids, each surface's
    runoff table with the depth running fastest, as Fortran stores ptab(depth, wetness), then one line per T21 point,
    as lstoreini reads them."""
    fill, depths, wets, runoff, rows = table_rows(product, library, hollows)
    fmt = lambda xs: ' '.join(repr(float(x)) for x in xs)
    lines = [' '.join(str(n) for n in (TABLE_VERSION, NUGP, NLSF, NLSD, NLSW, len(FLOORS), len(SURFACES))),
             fmt(fill), fmt(depths), fmt(wets)]
    lines += [fmt(t.T.ravel()) for t in runoff]                  # Fortran's column order: depth fastest in ptab(d, w)
    lines += [fmt(r) for r in rows]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text('\n'.join(lines) + '\n')
    return path


def setting(product_path, floor: str, surface: str, library: str, hollows: str, wetted_fraction: float,
            gravity_m_s2: float, pressure_pa: float, root: Path | None = None) -> dict:
    """The store as the configuration records it: the product and its hash, the floor, the surface, the storm library
    and hollow age whose routing and storage it reads, the wetted fraction, the pools' constants and the output."""
    if floor not in FLOORS or surface not in SURFACES:
        raise ValueError(f'the floor is one of {FLOORS} and the surface one of {SURFACES}')
    if library not in LIBRARIES or hollows not in HOLLOWS:
        raise ValueError(f'the library is one of {LIBRARIES} and the hollows one of {HOLLOWS}')
    if not 0.0 < float(wetted_fraction) <= 1.0:
        raise ValueError('the wetted fraction lies above 0 and at most 1')
    path = Path(product_path)
    product = read_product(path)
    resolved = path.resolve()
    shown = str(resolved.relative_to(root)) if root is not None and resolved.is_relative_to(root) else str(resolved)
    return dict(product=shown, sha256=hashlib.sha256(path.read_bytes()).hexdigest()[:16], schema=PRODUCT_SCHEMA,
                producer=product['producer'].get('files', {}).get('research/studies/storm_water/store_tables.py'),
                basis=product.get('basis'), floor=floor, surface=surface, library=library, hollows=hollows,
                wetted_fraction=float(wetted_fraction),
                wetted_fraction_status=('an Earth-tuned term carried as a parameter, not settled: Earth\'s schemes '
                                        'take 0.3 (JULES, Best et al. 2011, recalled via the plan)'),
                pool=dict(POOL, free_convection_m_s_k3=pool_free_convection(gravity_m_s2, pressure_pa),
                          form="Penman's evaporation of open water with no heat stored, PlaSim's surface transfer over "
                               'smooth water on the 10 m wind (variant b1)'),
                output=dict({str(k): v[0] for k, v in CODES.items()}, **{'160': RUNOFF_NAME}),
                earth_terms_kept=dict(wsmax_m=EARTH_TERMS['wsmax']['value_m'],
                                      drhsfull=EARTH_TERMS['drhsfull']['value']),
                code='climate/gcm/land_store.py (pool_evaporation, store_step; lspool and lsstep in Fortran)')


def namelist_values(s: dict) -> dict:
    """landmod_namelist's entries for a store setting (setting's record)."""
    p = s['pool']
    return dict(NLSTORE='1', NLSFLOOR=str(FLOORS.index(s['floor']) + 1), NLSSURF=str(SURFACES.index(s['surface']) + 1),
                DLSWET=repr(float(s['wetted_fraction'])), DLSALBW=repr(float(p['albedo'])),
                DLSEMW=repr(float(p['emissivity'])), DLSZ0W=repr(float(p['roughness_m'])),
                DLSZW=repr(float(p['wind_height_m'])), DLSUMIN=repr(float(p['least_wind_m_s'])),
                DLSCFREE=repr(float(p['free_convection_m_s_k3'])), DLSLB=repr(float(p['louis_b'])),
                DLSLC=repr(float(p['louis_c'])), DLSLD=repr(float(p['louis_d'])),
                DLSSBK=repr(float(STEFAN_BOLTZMANN_PLASIM)))


def land_mask_cells(landmap: Path) -> np.ndarray:
    """The land points of a run's land mask (the runner's moon_landmask.sra, PlaSim's order)."""
    text = Path(landmap).read_text().split('\n')
    values = np.array(' '.join(text[1:]).split(), float)
    if values.size != NUGP:
        raise ValueError(f'{landmap} holds {values.size} points, not T21\'s {NUGP}')
    return values > 0.5


# ---------------------------------------------------------------------------------------------------------------------
# Fortran: landmod.f90's module variables, namelist, broadcasts, restart, steps and routines.

DECLARATIONS = f"""!     {LAND_STORE_MARKER} (climate/gcm/land_store.py): nlstore 1 adds it in the cells its table gives
      integer, parameter :: NLSF = {NLSF}, NLSD = {NLSD}, NLSW = {NLSW}   ! the table's fills, rain depths and wetnesses
      integer, parameter :: NLSFLOORS = {len(FLOORS)}, NLSSURFS = {len(SURFACES)}, NLSOUT = {NLSOUT}, NLSPC = 19
      integer, parameter :: LSVERSION = {TABLE_VERSION}   ! landstore.dat's layout
      integer :: nlstore  = 0      ! 1: the slopes' runoff, the hollows' store, its evaporation and seepage
      integer :: nlsfloor = 1      ! the floor design (1 clay seal, 2 sandy banks above 1 m, 3 a seal ten times tighter)
      integer :: nlssurf  = 1      ! the slopes' surface (1 sandy loam, 2 loam, 3 sealed loam, 4 biological crust)
      real :: dlswet   = {WETTED_FRACTION_EARTH!r}       ! the convective rain's wetted share of a cell (Earth's schemes' 0.3; open)
      real :: dlsalbw  = {POOL['albedo']!r}      ! the pools' albedo
      real :: dlsemw   = {POOL['emissivity']!r}      ! the pools' emissivity
      real :: dlsz0w   = {POOL['roughness_m']!r}    ! the pools' roughness (m)
      real :: dlszw    = {POOL['wind_height_m']!r}      ! the height of the wind over the pools (m): PlaSim's 10 m wind
      real :: dlsumin  = {POOL['least_wind_m_s']!r}       ! the least wind over the pools (m/s)
      real :: dlscfree = {EXOPLASIM_FREE_CONVECTION!r}    ! the free convection off a warm pool (m s-1 K-1/3; the runner sets the planet's)
      real :: dlslb = {POOL['louis_b']!r}, dlslc = {POOL['louis_c']!r}, dlsld = {POOL['louis_d']!r}   ! Louis's constants, PlaSim's vdiff_b, vdiff_c, vdiff_d
      real :: dlssbk = {STEFAN_BOLTZMANN_PLASIM!r}    ! Stefan and Boltzmann's constant, radmod's
      real :: dlsv(NHOR) = 0.      ! the store's water (m over the cell)
      real :: dlsslow(NHOR) = 0.   ! the slow store's total, which leaves the model (m over the cell)
      real :: dlsfw(NHOR) = 0.     ! the store's water area (share of the cell)
      real :: dlses(NHOR) = 0.     ! the store's evaporation in this step (m/s, upward)
      real :: dlscap(NHOR) = 0.    ! its capacity (m); none where the cell has no store
      real :: dlsarea(NHOR,NLSF) = 0., dlsseep(NHOR,NLSF) = 0., dlsshare(NHOR,3) = 0.
      real :: dlsfill(NLSF) = 0., dlsdgrid(NLSD) = 0., dlswgrid(NLSW) = 0., dlsrtab(NLSD,NLSW) = 0.
      real :: dlsflx(NHOR,NLSOUT) = 0.   ! this step's runoff, water, area, evaporation, seepage, recharge, rivers' share
      real :: dlsacc(NHOR,NLSOUT) = 0.   ! their sums since the last output
      integer :: nlsacc = 0
"""
NAMELIST_NAMES = ('nlstore', 'nlsfloor', 'nlssurf', 'dlswet', 'dlsalbw', 'dlsemw', 'dlsz0w', 'dlszw', 'dlsumin',
                  'dlscfree', 'dlslb', 'dlslc', 'dlsld', 'dlssbk')
INTEGER_NAMES = ('nlstore', 'nlsfloor', 'nlssurf')
_DWATCINI = "      real    :: dwatcini = 0           ! water content of soil (m)\n"
_NAMELIST = "     &                ,newsurf,rinifor,nwatcini,dwatcini,dgroundalb\n"
_BROADCAST = "      call mpbcrn(dgroundalb,2)\n"
_RESTART_GET = "       call mpgetgp('dalbcl2' ,dalbcl2 ,NHOR,    14)\n"
_RESTART_PUT = "      call mpputgp('dalbcl2' ,dalbcl2 ,NHOR,14)\n"
_LANDINI_END = "      return\n      end subroutine landini\n"
_SOILSTEP = "      call soilstep\n"
_WANDR = "      if(nlandw==1) then\n       call wandr\n      else\n"
_AFTER_WANDR = "      end subroutine wandr\n"

ROUTINES = """
!     ====================
!     SUBROUTINE LSTOREINI
!     ====================

      subroutine lstoreini
      use landmod
!
!     Terluna: the hollows' store's tables from landstore.dat (land_store.write_table, from the storm-water study's
!     product): the grids of fills, rain depths and wetnesses, the chosen surface's runoff shares (nlssurf), and per
!     point the store's capacity, its water area at each fill, its seepage under the chosen floor (nlsfloor) and the
!     chosen surface's routing shares. The store's area follows from the water the restart gave it (none on a cold
!     start). The store needs the prognostic land.
!
      real zcap(NUGP),zarea(NUGP,NLSF),zseep(NUGP,NLSF),zshare(NUGP,3)
      real zrow(1+NLSF*(1+NLSFLOORS)+3*NLSSURFS)
      real ztab(NLSD,NLSW)
      real zq
      integer ihead(7)
      integer jhor,jsurf,jf
!
      if (nlandt /= 1 .or. nlandw /= 1) then
       write(nud,*) '*** the hollows store needs the prognostic land (nlandt and nlandw 1) ***'
       stop 1
      endif
      if (nlsfloor < 1 .or. nlsfloor > NLSFLOORS .or. nlssurf < 1 .or. nlssurf > NLSSURFS) then
       write(nud,*) '*** the hollows store: nlsfloor or nlssurf outside its table ***'
       stop 1
      endif
      if (ndiaggp2d < NLSOUT) then
       write(nud,*) '*** the hollows store writes its means in 2-d diagnostic slots: ndiaggp2d must be at least ',NLSOUT
       stop 1
      endif
      if (mypid == NROOT) then
       open(73,file='landstore.dat',form='formatted',status='old')
       read(73,*) ihead
       if (ihead(1) /= LSVERSION .or. ihead(2) /= NUGP .or. ihead(3) /= NLSF .or. ihead(4) /= NLSD                     &
     &     .or. ihead(5) /= NLSW .or. ihead(6) /= NLSFLOORS .or. ihead(7) /= NLSSURFS) then
        write(nud,*) '*** landstore.dat does not match the hollows store patch: ',ihead
        stop 1
       endif
       read(73,*) dlsfill
       read(73,*) dlsdgrid
       read(73,*) dlswgrid
       do jsurf=1,NLSSURFS
        read(73,*) ztab
        if (jsurf == nlssurf) dlsrtab(:,:)=ztab(:,:)
       enddo
       jf=1+NLSF*(1+NLSFLOORS)+3*(nlssurf-1)
       do jhor=1,NUGP
        read(73,*) zrow
        zcap(jhor)=zrow(1)
        zarea(jhor,:)=zrow(2:1+NLSF)
        zseep(jhor,:)=zrow(2+NLSF*nlsfloor:1+NLSF*(nlsfloor+1))
        zshare(jhor,:)=zrow(jf+1:jf+3)
       enddo
       close(73)
       write(nud,*) '* Read the hollows store: floor ',nlsfloor,' surface ',nlssurf,' wetted share ',dlswet
      endif
      call mpbcrn(dlsfill,NLSF)
      call mpbcrn(dlsdgrid,NLSD)
      call mpbcrn(dlswgrid,NLSW)
      call mpbcrn(dlsrtab,NLSD*NLSW)
      call mpscgp(zcap,dlscap,1)
      call mpscgp(zarea,dlsarea,NLSF)
      call mpscgp(zseep,dlsseep,NLSF)
      call mpscgp(zshare,dlsshare,3)
      do jhor=1,NHOR
       if (dls(jhor) > 0. .and. dlscap(jhor) > 0.) then
        call lscurve(NLSF,dlsfill,dlsarea(jhor,:),dlsseep(jhor,:),dlscap(jhor),dlsv(jhor),dlsfw(jhor),zq)
       else
        dlsv(jhor)=0.
        dlsfw(jhor)=0.
       endif
      enddo
      dlsacc(:,:)=0.
      nlsacc=0
      return
      end subroutine lstoreini

!     =====================
!     SUBROUTINE LSTOREEVAP
!     =====================

      subroutine lstoreevap
      use landmod
!
!     Terluna: the store's evaporation in this step (variant b1 of the storm-water plan; land_store.store_evaporation):
!     its water area times Penman's evaporation of open water with no heat stored (lspool), from the sunlight and
!     longwave down at the ground, the lowest level's air carried to the ground and the 10 m wind, at most the store's
!     water; none under snow. It joins the land's evaporation, the lowest layer's vapour (and the diagnostic of the
!     surface's) and the latent heat the cell's one skin loses; lstorewat takes it from the store, not the bucket.
!
      real zc(NLSPC)
      real zsw,zlw,zta,zpool,ztw,zv,zes
      integer jhor
!
      zc(:)=(/ga,gascon,rdbrv,acpd,alv,ra1,ra2,ra4,tmelt,dlssbk,dlsalbw,dlsemw,dlsz0w,dlszw,dlsumin,dlscfree,            &
     &        dlslb,dlslc,dlsld/)
      dlses(:)=0.
      do jhor=1,NHOR
       if (dls(jhor) > 0. .and. dlscap(jhor) > 0. .and. dlsfw(jhor) > 0. .and. dsnowz(jhor) <= 0.) then
        zsw=dswfl(jhor,NLEP)/MAX(1.-dalb(jhor),0.05)
        zlw=dlwfl(jhor,NLEP)+dlssbk*dt(jhor,NLEP)**4
        zta=dt(jhor,NLEV)*sigma(NLEV)**(-akap)
        call lspool(zsw,zlw,zta,dq(jhor,NLEV),dp(jhor),SQRT(du10(jhor)**2+dv10(jhor)**2),zc,zpool,ztw,zv)
        zes=dlsfw(jhor)*zpool
        if (zes > 0.) zes=MIN(zes,dlsv(jhor)/deltsec)
        dlses(jhor)=zes
        devap(jhor)=devap(jhor)-zes
        dlhfl(jhor)=dlhfl(jhor)-zes*alv*1000.
        dqdt(jhor,NLEV)=dqdt(jhor,NLEV)+zes*1000.*ga/(dp(jhor)*dsigma(NLEV))
        if (ndiaggp == 1) dgp3d(jhor,NLEV,13)=dgp3d(jhor,NLEV,13)+zes
       endif
      enddo
      return
      end subroutine lstoreevap

!     ====================
!     SUBROUTINE LSTOREWAT
!     ====================

      subroutine lstorewat
      use landmod
!
!     Terluna: the store's water in this step, between tands and wandr (lsstep): the bucket keeps what the store does
!     not take, and gives back the store's evaporation, which devap carries; the slow store gains the seepage and the
!     recharge below the roots.
!
      real zdwat,zvn,zfwn,zroff,zseep,zrech,zspill,zriver,zinfl
      integer jhor
!
      do jhor=1,NHOR
       dlsflx(jhor,:)=0.
       if (dls(jhor) > 0. .and. dlscap(jhor) > 0.) then
        call lsstep(deltsec,dprc(jhor),dprl(jhor),dlses(jhor),dlswet,dwatc(jhor),dwmax(jhor),dlsfw(jhor),dlsv(jhor),     &
     &              dlscap(jhor),NLSF,dlsfill,dlsarea(jhor,:),dlsseep(jhor,:),dlsshare(jhor,:),NLSD,NLSW,dlsdgrid,        &
     &              dlswgrid,dlsrtab,zdwat,zvn,zfwn,zroff,zseep,zrech,zspill,zriver,zinfl)
        dwater(jhor)=dwater(jhor)+zdwat
        dlsv(jhor)=zvn
        dlsfw(jhor)=zfwn
        dlsslow(jhor)=dlsslow(jhor)+(zseep+zrech)*deltsec
        dlsflx(jhor,1)=zroff
        dlsflx(jhor,2)=zvn
        dlsflx(jhor,3)=zfwn
        dlsflx(jhor,4)=-dlses(jhor)
        dlsflx(jhor,5)=zseep
        dlsflx(jhor,6)=zrech
        dlsflx(jhor,7)=zriver+zspill
       endif
      enddo
      return
      end subroutine lstorewat

!     =====================
!     SUBROUTINE LSTOREROFF
!     =====================

      subroutine lstoreroff
      use landmod
!
!     Terluna: after wandr, the store scheme's water for the rivers (the runoff's rivers' share and the spill) joins
!     the bucket's runoff, which roffstep carries. The step's quantities join their sums, and PlaSim's 2-d diagnostic
!     slots 1 to NLSOUT hold their mean since the last output, which outdiag writes as codes 1 to NLSOUT
!     (land_store.CODES); the sums start again after a step that ends an output interval.
!
      integer jout
!
      where (dls(:) > 0. .and. dlscap(:) > 0.) drunoff(:)=drunoff(:)+dlsflx(:,7)
      dlsacc(:,:)=dlsacc(:,:)+dlsflx(:,:)
      nlsacc=nlsacc+1
      do jout=1,NLSOUT
       dgp2d(:,jout)=dlsacc(:,jout)/real(nlsacc)
      enddo
      if (nafter > 0) then
       if (mod(nhcstp,nafter) == 0) then
        dlsacc(:,:)=0.
        nlsacc=0
       endif
      endif
      return
      end subroutine lstoreroff

!     =================
!     SUBROUTINE LSPOOL
!     =================

      subroutine lspool(psw,plw,pta,pqa,pp,pu,pc,pevap,ptw,ptrans)
!
!     Terluna: Penman's evaporation of open water with no heat stored (land_store.pool_evaporation step for step): the
!     water temperature ptw (K) at which the water's net radiation balances its latent and sensible heat with the
!     transfer there, by halving a bracket of 40 K about the air's temperature 48 times, and there the evaporation
!     pevap (m/s of water, upward) and the transfer ptrans (m/s). The transfer is PlaSim's surface form over smooth
!     water: the wind pu (at least pc(15)) times (kappa/ln(z/z0+1))**2 and Louis's factor of the bulk Richardson number
!     between the air's and the water's virtual temperatures, at least the free convection off a warm pool. psw and
!     plw are the sunlight and longwave down (W/m2), pta and pqa the air's temperature (K) and humidity, pp the pressure
!     (Pa), pc the constants in land_store.pool_constants' order.
!
      real psw,plw,pta,pqa,pp,pu,pevap,ptw,ptrans
      real pc(19)
      real zepsv,zs,zrho,zcn,zlo,zhi,zmid
      integer jit
!
      zepsv=1./pc(3)-1.
      zs=MAX(pu,pc(15))
      zrho=pp/(pc(2)*pta*(1.+zepsv*pqa))
      zcn=(0.4/LOG(pc(14)/pc(13)+1.))**2
      zlo=pta-40.
      zhi=pta+40.
      do jit=1,48
       zmid=0.5*(zlo+zhi)
       if (lsres(zmid,lstrans(zmid)) > 0.) then
        zlo=zmid
       else
        zhi=zmid
       endif
      enddo
      ptw=0.5*(zlo+zhi)
      ptrans=lstrans(ptw)
      pevap=zrho*ptrans*(lsqs(ptw)-pqa)/1000.
      return
!
      contains
!
      real function lsqs(pt)
!     PlaSim's saturation humidity (landstep's dqs)
      real pt,zq
      zq=pc(3)*pc(6)*EXP(pc(7)*(pt-pc(9))/(pt-pc(8)))/pp
      lsqs=zq/(1.-(1./pc(3)-1.)*zq)
      end function lsqs
!
      real function lstrans(pt)
!     the transfer velocity at water temperature pt
      real pt,zqw,zthva,zthvw,zri,zfh
      zqw=lsqs(pt)
      zthva=pta*(1.+zepsv*pqa)
      zthvw=pt*(1.+zepsv*zqw)
      zri=pc(1)*pc(14)*(zthva-zthvw)/(zs*zs*zthva)
      if (zri <= 0.) then
       zfh=1.-3.*pc(17)*zri/(1.+3.*pc(17)*pc(18)*zcn*SQRT(-zri*(pc(14)/pc(13)+1.)))
      else
       zfh=1./(1.+3.*pc(17)*zri*SQRT(1.+pc(19)*zri))
      endif
      lstrans=MAX(zs*zcn*zfh,pc(16)*MAX(zthvw-zthva,0.)**(1./3.))
      end function lstrans
!
      real function lsres(pt,pv)
!     the water's net radiation less its latent and sensible heat (W/m2) at pt with the transfer pv
      real pt,pv,zrn
      zrn=(1.-pc(11))*psw+pc(12)*(plw-pc(10)*pt**4)
      lsres=zrn-zrho*pc(5)*pv*(lsqs(pt)-pqa)-zrho*pc(4)*pv*(pt-pta)
      end function lsres
!
      end subroutine lspool

!     ==================
!     SUBROUTINE LSCURVE
!     ==================

      subroutine lscurve(kf,pfill,parea,pseep,pcap,pv,pfw,pq)
!
!     Terluna: the store's water area pfw (share of the cell) and seepage pq (m/s) at its water pv (m) from its tables
!     at the common fills pfill (land_store.curve): linear in the square root of the fill, nil at nil water.
!
      integer kf
      real pfill(kf),parea(kf),pseep(kf),pcap,pv,pfw,pq
      real zf,zs,zw
      integer k
!
      pfw=0.
      pq=0.
      if (pcap <= 0. .or. pv <= 0.) return
      zf=MIN(pv/pcap,1.)
      zs=SQRT(zf)
      do k=2,kf
       if (SQRT(pfill(k)) >= zs) exit
      enddo
      k=MIN(k,kf)
      zw=(zs-SQRT(pfill(k-1)))/MAX(SQRT(pfill(k))-SQRT(pfill(k-1)),1.e-30)
      pfw=parea(k-1)+zw*(parea(k)-parea(k-1))
      pq=pseep(k-1)+zw*(pseep(k)-pseep(k-1))
      return
      end subroutine lscurve

!     ==================
!     SUBROUTINE LSSHARE
!     ==================

      real function lsshare(kd,kw,pdg,pwg,ptab,pd,pw)
!
!     Terluna: the share of a step's rain that runs off (land_store.runoff_share): bilinear in the rain's depth on its
!     patch pd (mm) and the bucket's relative wetness pw, both held within the table's grids.
!
      integer kd,kw
      real pdg(kd),pwg(kw),ptab(kd,kw),pd,pw
      real zd,zw,za,zb
      integer i,j
!
      zd=MIN(MAX(pd,pdg(1)),pdg(kd))
      zw=MIN(MAX(pw,pwg(1)),pwg(kw))
      do i=1,kd-1
       if (zd <= pdg(i+1)) exit
      enddo
      i=MIN(i,kd-1)
      do j=1,kw-1
       if (zw <= pwg(j+1)) exit
      enddo
      j=MIN(j,kw-1)
      za=(zd-pdg(i))/(pdg(i+1)-pdg(i))
      zb=(zw-pwg(j))/(pwg(j+1)-pwg(j))
      lsshare=(1.-za)*(1.-zb)*ptab(i,j)+za*(1.-zb)*ptab(i+1,j)+(1.-za)*zb*ptab(i,j+1)+za*zb*ptab(i+1,j+1)
      return
      end function lsshare

!     =================
!     SUBROUTINE LSSTEP
!     =================

      subroutine lsstep(pdt,pprc,pprl,pes,pwet,pw,pwmax,pfw,pv,pcap,kf,pfill,parea,pseep,pshare,kd,kw,pdg,pwg,ptab,       &
     &                  pdwat,pvn,pfwn,proff,pseepo,precho,pspill,priver,pinfl)
!
!     Terluna: one step of the store in a cell (land_store.store_step), after its evaporation pes (m/s, upward, at most
!     its water where positive) has joined the land's: the rain on its water area, the slopes' runoff by infiltration
!     excess on the rest (the convective rain on its wetted share pwet, the large-scale rain on the whole cell), the
!     runoff's routing by the shares pshare (hollows, rivers, recharge below the roots; the rest soaks back), the
!     seepage at the step's start, the spill above the capacity. pdwat is the change to the bucket's dwater (m/s).
!
      integer kf,kd,kw
      real pdt,pprc,pprl,pes,pwet,pw,pwmax,pfw,pv,pcap,pfill(kf),parea(kf),pseep(kf),pshare(3)
      real pdg(kd),pwg(kw),ptab(kd,kw)
      real pdwat,pvn,pfwn,proff,pseepo,precho,pspill,priver,pinfl
      real lsshare
      real zrel,zps,zdc,zdl,zq0,zv,zfw0,zq
!
      zrel=0.
      if (pwmax > 0.) zrel=MIN(MAX(pw/pwmax,0.),1.)
      zps=pfw*(pprc+pprl)
      zdc=pprc*pdt*1000./MAX(pwet,1.e-6)
      zdl=pprl*pdt*1000.
      proff=(lsshare(kd,kw,pdg,pwg,ptab,zdc,zrel)*pprc+lsshare(kd,kw,pdg,pwg,ptab,zdl,zrel)*pprl)*(1.-pfw)
      pinfl=zps+pshare(1)*proff
      call lscurve(kf,pfill,parea,pseep,pcap,pv,zfw0,zq0)
      zv=pv-pes*pdt+pinfl*pdt
      pseepo=MIN(zq0,MAX(zv,0.)/pdt)
      zv=zv-pseepo*pdt
      pspill=MAX(zv-pcap,0.)/pdt
      pvn=MAX(zv-pspill*pdt,0.)
      call lscurve(kf,pfill,parea,pseep,pcap,pvn,pfwn,zq)
      priver=pshare(2)*proff
      precho=pshare(3)*proff
      pdwat=pes-zps-(pshare(1)+pshare(2)+pshare(3))*proff
      return
      end subroutine lsstep
"""

LANDMOD_EDITS = [
    (_DWATCINI, _DWATCINI + DECLARATIONS),
    (_NAMELIST, _NAMELIST[:-1] + "  &\n"
                "     &                ," + ','.join(NAMELIST_NAMES[:7]) + "  &\n"
                "     &                ," + ','.join(NAMELIST_NAMES[7:]) + "\n"),
    (_BROADCAST, _BROADCAST + ''.join(f"      call mpbc{'i' if n in INTEGER_NAMES else 'r'}({n})"
                                      + (f"   ! {LAND_STORE_MARKER}\n" if i == 0 else "\n")
                                      for i, n in enumerate(NAMELIST_NAMES))),
    (_RESTART_GET, _RESTART_GET
     + f"       if (nlstore == 1) then   ! {LAND_STORE_MARKER}: its water, none in a restart without it\n"
       "        ilsread = 0\n"
       "        if (mypid == NROOT) call restart_has_array('dlstorev',ilsread)\n"
       "        call mpbci(ilsread)\n"
       "        if (ilsread == 1) then\n"
       "         call mpgetgp('dlstorev',dlsv   ,NHOR,     1)\n"
       "         call mpgetgp('dlstoreg',dlsslow,NHOR,     1)\n"
       "        endif\n"
       "       endif\n"),
    (_LANDINI_END, f"      if (nlstore == 1) call lstoreini   ! {LAND_STORE_MARKER}: its tables\n" + _LANDINI_END),
    (_RESTART_PUT, _RESTART_PUT
     + f"      if (nlstore == 1) then   ! {LAND_STORE_MARKER}\n"
       "       call mpputgp('dlstorev',dlsv   ,NHOR, 1)\n"
       "       call mpputgp('dlstoreg',dlsslow,NHOR, 1)\n"
       "      endif\n"),
    (_SOILSTEP, f"      if (nlstore == 1) call lstoreevap   ! {LAND_STORE_MARKER}: its evaporation joins the land's\n"
                + _SOILSTEP),
    (_WANDR, "      if(nlandw==1) then\n"
             f"       if (nlstore == 1) call lstorewat   ! {LAND_STORE_MARKER}: the slopes' runoff and the store\n"
             "       call wandr\n"
             "       if (nlstore == 1) call lstoreroff   ! its water for the rivers\n"
             "      else\n"),
    (_AFTER_WANDR, _AFTER_WANDR + ROUTINES),
]
# Only landmod.f90 changes: the year-file fields go out through PlaSim's own 2-d diagnostic slots, so outmod.f90, which
# every run's key hashes, stays as it is. The store's restart reads take restart_has_array (exoplasim_run's switch
# patch to restartmod.f90), which the runner applies with this one.
LAND_STORE_PATCHES = [('landmod.f90', LAND_STORE_MARKER, LANDMOD_EDITS)]
