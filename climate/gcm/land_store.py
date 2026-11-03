"""The hollows' store in PlaSim's land scheme: stage (b) of the storm-water study (3 November 2026; the record is the
section "The hollows' store in PlaSim's land" in README.md beside this file, the design note
session-tools/agents/storm_water/stage_b_design.md and the second pass's report stage_b2_report.md on the research
drive).

PlaSim holds each land cell's water in one bucket of wsmax (0.5 m) that sheds water only when full (landmod.f90's
wandr), so on the design years the drylands return all their rain to the air. The storm-water study
(research/studies/storm_water on domain/hydrosphere) worked out offline what the drylands' storms do instead, and its
product results/store_tables.json (schema terluna.research.storm-water-store-tables/2) carries what this scheme reads.
With nlstore 1, in every land cell the product gives a store, each step

1. the slopes shed runoff as the chosen storm library's storms would: each millimetre of the step's convective rain
   runs off in the share the product's table gives for the chosen library and surface against a running mean of the
   cell's convective rain, which stands for the 3-day mean the study's tables are built on, and the bucket's relative
   wetness (the study's phase 2a: the library's storms at random times within a 3-day mean, the soil redistributing
   and drying between them, Green and Ampt at lunar gravity); the large-scale rain runs off in the share its own table
   gives against its own running mean, the mean's large-scale rain falling as one steady rain of 2 mm/h as the study
   reads it; both running means are exponential, with a time constant of RUNNING_MEAN_DAYS;
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

The store's water area and seepage come from tables at shares of its capacity, with the cell's water shared among its
hollows by their own yearly inflow (each hollow at most full) for the chosen library, surface and hollow age, read
linearly in the square root of the share, which is exact for the paraboloid bowls and gives nil area at nil volume.
Where a run's restart holds no store, the store starts at the product's balanced volume for its floor, at which its
mean inflow meets its mean losses, and the running means at the record's mean rain. Without nlstore every run is as
before: every edit only inserts lines, and the inserted statements run only with nlstore 1. Water is conserved exactly
each step: the rain equals the bucket's change, its own evaporation and runoff, the water sent to the rivers, the
store's change, its evaporation and the slow store's gain (land_step, the tests).

The terms PlaSim keeps from Earth are left as they are and listed (EARTH_TERMS) with the Earth calibrations the scheme
carries, and the scheme's stated limits are recorded (LIMITS). LAND_STORE_PATCHES put the scheme into landmod.f90
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
PRODUCT_SCHEMA = 'terluna.research.storm-water-store-tables/2'
TABLE_VERSION = 2                 # landstore.dat's layout, checked by lstoreini
TABLE_NAME = 'landstore.dat'      # the run's copy of the tables, in PlaSim's working folder
FLOORS = ('clay_floored', 'sandy_banked', 'tight_seal')            # the floor designs, nlsfloor 1-3 (phase 2b's)
SURFACES = ('sandy_loam', 'loam', 'sealed_loam', 'biocrust')        # the slopes' surfaces, nlssurf 1-4 (phase 1's)
LIBRARIES = ('walnut_gulch', 'cm1_6km')                             # the storm libraries the runoff and routing follow
HOLLOWS = ('young', 'mature')                                       # the hollows' age (phase 2a's states)
SHARES = ('hollows', 'rivers', 'recharge')                          # of the slopes' runoff; the rest soaks back
NLON, NLAT = 64, 32                                                 # T21: the product's cells are j * NLON + i
NUGP = NLON * NLAT
NLSF, NLSR, NLSW = 12, 29, 13     # fills, rain means and wetnesses of the product's grids (store_tables.py)
NLSOUT = 7
NLSROW = 1 + NLSF * (1 + len(FLOORS)) + len(SHARES) + len(FLOORS)  # a surface's entries in a point's row of the table
MM_DAY_PER_M_S = 86400.0 * 1000.0
# The running means of the rain: exponential, one state per cell, standing for the 3-day means (a tenth of the lunar
# day) on which the product's runoff tables are built; the time constant is a choice of this scheme, tested against
# the 3-day means on PlaSim's own 5.5-hourly rain (stage_b2_report.md).
RUNNING_MEAN_DAYS = 1.5
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
# The year-file fields: PlaSim's 2-d diagnostic slots 1-7 (dgp2d, which outdiag writes as codes 1-7 and nothing else
# fills), each the mean over the output interval; the runner sets ndiaggp2d to NLSOUT and names the codes.
CODES = {1: ('mrros', 'surface_runoff_from_the_slopes_as_the_storm_librarys_storms', 'm s-1'),
         2: ('lstore', 'water_in_the_hollows_store_over_the_cell', 'm'),
         3: ('lstorea', 'water_area_share_of_the_hollows_store', '1'),
         4: ('lstoree', 'evaporation_from_the_hollows_store_downward', 'm s-1'),
         5: ('lstores', 'seepage_from_the_hollows_store_into_the_slow_store', 'm s-1'),
         6: ('lstoreg', 'recharge_below_the_roots_from_the_slopes_and_channel_beds', 'm s-1'),
         7: ('lstorer', 'water_the_store_scheme_sends_to_the_rivers', 'm s-1')}
RUNOFF_NAME = 'mrro'                # PlaSim's own code 160: the local runoff with the rivers' convergence, as mkradv
                                    # leaves drunoff (the river store's tendency over land, the rivers' inflow at sea)
# PlaSim's own Earth-tuned terms the build leaves as they are (landmod.f90), and the Earth calibrations the scheme
# carries, with what the study's evidence says.
EARTH_TERMS = {
    'wsmax': dict(value_m=0.5, where='landmod.f90 WSMAX_EARTH, unchanged by the design run',
                  evidence='the drylands keep the bucket 0.4-8% full on the design years; the study reads it as the '
                           'plant-available water between the wilting point and field capacity of a soil of '
                           'unspecified depth, and its capacity sets how much rain the slopes soak before evaporating; '
                           'kept for the first store run and decided with the soils, which the project chooses with '
                           'the biome work after a maintainability study'),
    'drhsfull': dict(value=0.4, where='landmod.f90 drhsfull: the wetness factor reaches 1 at 0.4 wsmax',
                     evidence='an Earth calibration of the bucket\'s evaporation; the drylands\' factor is 0.05-0.08, '
                              'which with the 28 mm/day daytime potential sets their evaporation; kept for the first '
                              'store run and decided with wsmax and the soils'),
    'storms': dict(where='the product\'s convective runoff tables and routing shares (storms.py, phase 2a)',
                   evidence='Walnut Gulch\'s measured point storms on Earth (a family fitted to Osborn 1983) and CM1\'s '
                            'model storms over wet equatorial ground at 6 km; neither is the drylands\' own storms'),
    'large_scale_rate': dict(value_mm_h=2.0, where='the product\'s large-scale runoff tables (run.LS_RATE_MM_H)',
                             evidence='an estimate of the large-scale rain\'s steady rate; it carries 0.2% of the '
                                      'drylands\' runoff or less'),
    'green_ampt': dict(where='the product\'s runoff tables (infiltration.py, wetting.py)',
                       evidence='Rawls et al. (1983) class means and Rawls et al. (1990) seal, and Brooks and Corey '
                                'curves for the redistribution between storms, at lunar gravity; the surfaces are '
                                'design choices'),
    'drying_between_storms': dict(value=1.26, where='the product\'s convective runoff tables (extract.PT_ALPHA)',
                                  evidence='the soil dries between a mean\'s storms at Priestley and Taylor\'s demand '
                                           'with Earth\'s coefficient of 1.26, weighted by the design years\' '
                                           'convective rain'),
    'louis': dict(values=[5.0, 5.0, 5.0], where='POOL louis_b, louis_c, louis_d',
                  evidence='PlaSim\'s own Earth-calibrated stability factors, here over the pools'),
    'free_convection': dict(where='POOL, pool_free_convection',
                            evidence='Miller, Beljaars and Palmer\'s (1992) still-air constant over Earth\'s water, '
                                     'carried to the planet as (g / p_s)^(1/3)'),
}
# The scheme's stated limits, which the configuration records.
LIMITS = {
    'sandy_banks': ('one store per cell leaves the sandy banks\' pools short by half or more of the evaporation of '
                    'phase 2b\'s pools under any sharing of the water among the hollows, because the store\'s water '
                    'rises above the banks\' 1 m stage, which drains it; no second store is built'),
    'gust': ('the pools evaporate on the 10 m wind without the boundary layer\'s gust, which PlaSim forms only inside '
             'surflx (fluxmod.f90, which every run\'s key hashes); their daytime evaporation runs about a tenth low'),
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
    shares of its capacity (lscurve): linear in the square root of the share between the table's, nil at nil volume
    and held at the full store's values above its capacity. fill (NLSF,), area and seep (cell, NLSF), cap and v
    (cell,)."""
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


def runoff_share(rgrid, wgrid, table, r_mm_day, w):
    """The slopes' runoff per unit of rain (lsshare): bilinear in the rain's running mean (mm/day) and the bucket's
    relative wetness, both held within the table's grids. table (NLSR, NLSW) for one surface (and library)."""
    r = np.clip(np.asarray(r_mm_day, float), rgrid[0], rgrid[-1])
    w = np.clip(np.asarray(w, float), wgrid[0], wgrid[-1])
    i = np.clip(np.searchsorted(rgrid, r, side='left') - 1, 0, rgrid.size - 2)    # rgrid[i] < r <= rgrid[i+1]
    j = np.clip(np.searchsorted(wgrid, w, side='left') - 1, 0, wgrid.size - 2)
    a = (r - rgrid[i]) / (rgrid[i + 1] - rgrid[i])
    b = (w - wgrid[j]) / (wgrid[j + 1] - wgrid[j])
    return ((1.0 - a) * (1.0 - b) * table[i, j] + a * (1.0 - b) * table[i + 1, j] + (1.0 - a) * b * table[i, j + 1]
            + a * b * table[i + 1, j + 1])


def running_mean(r, p, dt, tau_s):
    """An exponential running mean r of the rain p (m/s) after a step dt (s): the mean relaxes towards the step's rain
    with the time constant tau_s, exactly for rain steady through the step."""
    return r + (1.0 - math.exp(-dt / tau_s)) * (np.asarray(p, float) - r)


def store_step(dt, prc, prl, es, tau_s, rc, rl, w_bucket, wmax, fw, v, cap, fill, area, seep, shares, rgrid, wgrid,
               ctab, ltab):
    """One step of the store in each cell (lsstep), after its evaporation es (m/s, upward positive, at most v / dt
    where positive) has joined the land's: the running means of the rain take this step's, the rain on the store's
    water area, the slopes' runoff on the rest as the storm library's storms (the convective rain through ctab at its
    running mean, the large-scale rain through ltab at its own), its routing, the seepage at the step's start, the
    spill above capacity. prc and prl are the convective and large-scale rain (m/s), rc and rl their running means at
    the step's start, w_bucket and wmax the bucket's water and capacity (m), fw and v the store's area share and water
    at the step's start, shares (cell, 3) the hollows', rivers' and recharge's shares of the runoff. Returns the change
    to PlaSim's dwater (m/s) and the step's quantities."""
    rel = np.where(wmax > 0.0, np.minimum(np.maximum(w_bucket / np.where(wmax > 0.0, wmax, 1.0), 0.0), 1.0), 0.0)
    rcn = running_mean(rc, prc, dt, tau_s)
    rln = running_mean(rl, prl, dt, tau_s)
    ps = fw * (prc + prl)
    roff = (runoff_share(rgrid, wgrid, ctab, rcn * MM_DAY_PER_M_S, rel) * prc
            + runoff_share(rgrid, wgrid, ltab, rln * MM_DAY_PER_M_S, rel) * prl) * (1.0 - fw)
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
    return dict(dwater=dwater, v=vn, fw=fwn, rc=rcn, rl=rln, runoff=roff, seepage=seepage, recharge=recharge,
                spill=spill, rivers=rivers, inflow=inflow, rain_on_store=ps)


def store_evaporation(fw, v, pool_rate, dt, snow=None):
    """The store's evaporation this step (lstoreevap): its water area times the pools' rate (m/s, upward positive), at
    most its water where it evaporates; none under snow. Dew onto the store (a negative rate) is not limited."""
    es = np.asarray(fw, float) * np.asarray(pool_rate, float)
    es = np.where(es > 0.0, np.minimum(es, np.asarray(v, float) / dt), es)
    if snow is not None:
        es = np.where(np.asarray(snow, float) > 0.0, 0.0, es)
    return es


def start_state(tables, w, store=None, means=None):
    """The store's state as lstoreini leaves it: where the restart holds no store (store None), the balanced volume
    vbal, and where it holds no running means (means None), the record's mean rain rc0 and rl0; nothing in cells
    without a store. store is the restart's water, means its (rc, rl)."""
    on = tables['cap'] > 0.0
    v = np.where(on, tables['vbal'] if store is None else store, 0.0)
    rc, rl = (tables['rc0'], tables['rl0']) if means is None else means
    fw, _ = curve(tables['fill'], tables['area'], tables['seep'], tables['cap'], v)
    return dict(w=np.asarray(w, float), v=v, fw=fw, rc=np.where(on, rc, 0.0), rl=np.where(on, rl, 0.0))


def land_step(dt, state, forcing, tables, tau_s=RUNNING_MEAN_DAYS * 86400.0):
    """One land step of a store cell with PlaSim's bucket as soilstep and wandr have it (no snow): the store's
    evaporation, tands' dwater with the land's own evaporation devap_bucket (m/s, PlaSim's sign), the store's step and
    the bucket's. Returns the new state, the fluxes and the water budget's terms (m over the cell), which close
    exactly: rain = bucket change + bucket's evaporation + runoff to the rivers + store change + store's evaporation
    + slow store's gain."""
    w, v, fw = state['w'], state['v'], state['fw']
    es = store_evaporation(fw, v, forcing['pool_rate'], dt)
    devap = forcing['devap_bucket'] - es
    dwater = devap + forcing['prl'] + forcing['prc']
    s = store_step(dt, forcing['prc'], forcing['prl'], es, tau_s, state['rc'], state['rl'], w, tables['wmax'], fw, v,
                   tables['cap'], tables['fill'], tables['area'], tables['seep'], tables['shares'], tables['rgrid'],
                   tables['wgrid'], tables['ctab'], tables['ltab'])
    dwater = dwater + s['dwater']
    wn = w + dt * dwater
    bucket_runoff = np.maximum(0.0, wn - tables['wmax']) / dt
    clamp = np.maximum(-wn, 0.0)                                   # water PlaSim's clamp at nil would create
    wn = np.maximum(np.minimum(tables['wmax'], wn), 0.0)
    drunoff = bucket_runoff + s['rivers'] + s['spill']
    budget = dict(rain=(forcing['prc'] + forcing['prl']) * dt, bucket_change=wn - w,
                  bucket_evaporation=-forcing['devap_bucket'] * dt, to_rivers=drunoff * dt, store_change=s['v'] - v,
                  store_evaporation=es * dt, slow_store=(s['seepage'] + s['recharge']) * dt, clamp=clamp)
    return dict(w=wn, v=s['v'], fw=s['fw'], rc=s['rc'], rl=s['rl']), dict(s, es=es, devap=devap, drunoff=drunoff), budget


# ---------------------------------------------------------------------------------------------------------------------
# The product and the run's input file.

def read_product(path) -> dict:
    """The storm-water study's store tables (store_tables.py), refused unless they carry this patch's schema (version
    1, whose runoff tables read each step as a storm of its own and whose curves fill the hollows by their capacity,
    is refused by name), the grids this patch reads and every floor, surface, library and hollow age."""
    path = Path(path)
    product = json.loads(path.read_text())
    schema = product.get('schema')
    if schema == PRODUCT_SCHEMA.rsplit('/', 1)[0] + '/1':
        raise ValueError(f'{path} is a version 1 store-tables product (step-depth runoff tables, curves by capacity), '
                         f'which this patch does not read; regenerate it as {PRODUCT_SCHEMA}')
    if schema != PRODUCT_SCHEMA:
        raise ValueError(f'{path} is not a {PRODUCT_SCHEMA} product')
    g = product['grids']
    if (len(g['fill']) != NLSF or len(g['rain_mm_day']) != NLSR or len(g['wetness']) != NLSW
            or g['fill'][0] != 0.0 or g['fill'][-1] != 1.0 or g['rain_mm_day'][0] != 0.0
            or g['wetness'][0] != 0.0 or g['wetness'][-1] != 1.0):
        raise ValueError(f'{path}: its grids are not the {NLSF} fills, {NLSR} rain means and {NLSW} wetnesses this '
                         'patch reads')
    if (tuple(product['floors']) != FLOORS or tuple(product['surfaces']) != SURFACES
            or tuple(product['libraries']) != LIBRARIES or tuple(product['hollows']) != HOLLOWS):
        raise ValueError(f'{path}: its floors, surfaces, libraries or hollow ages are not the patch\'s')
    if product['grid'] != dict(nlon=NLON, nlat=NLAT):
        raise ValueError(f'{path}: its cells are not on PlaSim\'s T21 grid')
    return product


def table_rows(product: dict, library: str, hollows: str) -> tuple:
    """The input file's tables for one storm library and hollow age: the grids, each surface's convective and
    large-scale runoff tables (NLSR x NLSW), and per T21 point in PlaSim's order (north to south, from 0 E) the start
    of the running means of the convective and large-scale rain (m/s) and, for each surface, the capacity (m), the
    water area at each fill, the seepage at each fill under each floor (m/s), the three routing shares and the balanced
    volume under each floor (m); nothing where the product has no store."""
    if library not in LIBRARIES or hollows not in HOLLOWS:
        raise ValueError(f'the library is one of {LIBRARIES} and the hollows one of {HOLLOWS}')
    g = product['grids']
    ctab = np.array([product['convective_runoff'][library][s] for s in SURFACES], float)    # (surface, rain, wetness)
    ltab = np.array([product['large_scale_runoff'][s] for s in SURFACES], float)
    if ctab.shape[1:] != (NLSR, NLSW) or ltab.shape[1:] != (NLSR, NLSW):
        raise ValueError(f'the runoff tables are not {NLSR} rain means by {NLSW} wetnesses')
    rows = np.zeros((NUGP, 2 + len(SURFACES) * NLSROW))
    for c in product['cells']:
        k = int(c['j']) * NLON + int(c['i'])
        blocks = []
        for s in SURFACES:
            st = c['store'][library][s][hollows]
            if len(st['area_share']) != NLSF or any(len(st['seepage_m_s'][f]) != NLSF for f in FLOORS):
                raise ValueError(f'cell {c["j"]}, {c["i"]}: its curve is not {NLSF} fills')
            if st['capacity_m'] <= 0.0:
                blocks.append([0.0] * NLSROW)
                continue
            blocks.append([st['capacity_m'], *st['area_share'], *[x for f in FLOORS for x in st['seepage_m_s'][f]],
                           *[c['shares'][library][s][hollows][q] for q in SHARES],
                           *[st['balanced_m'][f] for f in FLOORS]])
        if not any(b[0] > 0.0 for b in blocks):
            continue
        rows[k] = [c['rain_m_s']['convective'], c['rain_m_s']['large_scale'], *[x for b in blocks for x in b]]
    return (np.array(g['fill'], float), np.array(g['rain_mm_day'], float), np.array(g['wetness'], float), ctab, ltab,
            rows)


def unpack_row(row, surface, floor):
    """A point's row as lstoreini reads it for one surface and floor: (rc0, rl0, cap, area, seep, shares, vbal)."""
    jf = 2 + NLSROW * SURFACES.index(surface)
    kf = FLOORS.index(floor) + 1
    jo = jf + 1 + NLSF * (1 + len(FLOORS))
    return (row[0], row[1], row[jf], row[jf + 1:jf + 1 + NLSF], row[jf + 1 + NLSF * kf:jf + 1 + NLSF * (kf + 1)],
            row[jo:jo + 3], row[jo + 3 + kf - 1])


def write_table(product: dict, library: str, hollows: str, path) -> Path:
    """landstore.dat: a header (layout, points, fills, rain means, wetnesses, floors, surfaces), the grids, each
    surface's convective runoff table for the library, then each surface's large-scale table, the rain mean running
    fastest as Fortran stores ptab(rain, wetness), then one line per T21 point, as lstoreini reads them."""
    fill, rgrid, wgrid, ctab, ltab, rows = table_rows(product, library, hollows)
    fmt = lambda xs: ' '.join(repr(float(x)) for x in xs)
    lines = [' '.join(str(n) for n in (TABLE_VERSION, NUGP, NLSF, NLSR, NLSW, len(FLOORS), len(SURFACES))),
             fmt(fill), fmt(rgrid), fmt(wgrid)]
    lines += [fmt(t.T.ravel()) for t in ctab]                    # Fortran's column order: rain fastest in ptab(r, w)
    lines += [fmt(t.T.ravel()) for t in ltab]
    lines += [fmt(r) for r in rows]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text('\n'.join(lines) + '\n')
    return path


def setting(product_path, floor: str, surface: str, library: str, hollows: str, gravity_m_s2: float,
            pressure_pa: float, root: Path | None = None) -> dict:
    """The store as the configuration records it: the product and its hash, the floor, the surface, the storm library
    and hollow age whose runoff, routing and storage it reads, the running means' form, the start, the pools'
    constants, the output, PlaSim's kept Earth terms and the scheme's stated limits."""
    if floor not in FLOORS or surface not in SURFACES:
        raise ValueError(f'the floor is one of {FLOORS} and the surface one of {SURFACES}')
    if library not in LIBRARIES or hollows not in HOLLOWS:
        raise ValueError(f'the library is one of {LIBRARIES} and the hollows one of {HOLLOWS}')
    path = Path(product_path)
    product = read_product(path)
    resolved = path.resolve()
    shown = str(resolved.relative_to(root)) if root is not None and resolved.is_relative_to(root) else str(resolved)
    return dict(product=shown, sha256=hashlib.sha256(path.read_bytes()).hexdigest()[:16], schema=PRODUCT_SCHEMA,
                producer=product['producer'].get('files', {}).get('research/studies/storm_water/store_tables.py'),
                basis=product.get('basis'), floor=floor, surface=surface, library=library, hollows=hollows,
                runoff=dict(reading='each millimetre of rain runs off in the share the product gives for the storm '
                                    'library and surface at a running mean of the cell\'s rain of its kind and the '
                                    'bucket\'s relative wetness (phase 2a\'s storms; the large-scale rain as a steady '
                                    '2 mm/h)',
                            running_mean='exponential', running_mean_days=RUNNING_MEAN_DAYS,
                            tables_built_on_mean_days=(product.get('basis') or {}).get('mean_days')),
                start=('where the restart holds no store, the product\'s balanced volume for the floor; where it holds '
                       'no running means, the record\'s mean rain'),
                pool=dict(POOL, free_convection_m_s_k3=pool_free_convection(gravity_m_s2, pressure_pa),
                          form="Penman's evaporation of open water with no heat stored, PlaSim's surface transfer over "
                               'smooth water on the 10 m wind (variant b1)'),
                output=dict({str(k): v[0] for k, v in CODES.items()}, **{'160': RUNOFF_NAME}),
                earth_terms_kept=dict(wsmax_m=EARTH_TERMS['wsmax']['value_m'],
                                      drhsfull=EARTH_TERMS['drhsfull']['value']),
                limits=dict(LIMITS),
                code='climate/gcm/land_store.py (pool_evaporation, store_step; lspool and lsstep in Fortran)')


def namelist_values(s: dict) -> dict:
    """landmod_namelist's entries for a store setting (setting's record)."""
    p = s['pool']
    return dict(NLSTORE='1', NLSFLOOR=str(FLOORS.index(s['floor']) + 1), NLSSURF=str(SURFACES.index(s['surface']) + 1),
                DLSTAU=repr(float(s['runoff']['running_mean_days']) * 86400.0), DLSALBW=repr(float(p['albedo'])),
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
      integer, parameter :: NLSF = {NLSF}, NLSR = {NLSR}, NLSW = {NLSW}   ! the table's fills, rain means and wetnesses
      integer, parameter :: NLSFLOORS = {len(FLOORS)}, NLSSURFS = {len(SURFACES)}, NLSOUT = {NLSOUT}, NLSPC = 19
      integer, parameter :: NLSROW = 1+NLSF*(1+NLSFLOORS)+3+NLSFLOORS   ! a surface's entries in a point's row
      integer, parameter :: LSVERSION = {TABLE_VERSION}   ! landstore.dat's layout
      integer :: nlstore  = 0      ! 1: the slopes' runoff, the hollows' store, its evaporation and seepage
      integer :: nlsfloor = 1      ! the floor design (1 clay seal, 2 sandy banks above 1 m, 3 a seal ten times tighter)
      integer :: nlssurf  = 1      ! the slopes' surface (1 sandy loam, 2 loam, 3 sealed loam, 4 biological crust)
      integer :: nlsread  = 0      ! 1: the restart held the store's water
      integer :: nlsrread = 0      ! 1: the restart held the running means of the rain
      real :: dlstau   = {RUNNING_MEAN_DAYS * 86400.0!r}   ! the rain's running means' time constant (s)
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
      real :: dlsrc(NHOR) = 0.     ! the running mean of the convective rain (m/s)
      real :: dlsrl(NHOR) = 0.     ! the running mean of the large-scale rain (m/s)
      real :: dlsvbal(NHOR) = 0., dlsrc0(NHOR) = 0., dlsrl0(NHOR) = 0.   ! the starts: balanced water, record's rain
      real :: dlsarea(NHOR,NLSF) = 0., dlsseep(NHOR,NLSF) = 0., dlsshare(NHOR,3) = 0.
      real :: dlsfill(NLSF) = 0., dlsrgrid(NLSR) = 0., dlswgrid(NLSW) = 0.
      real :: dlsctab(NLSR,NLSW) = 0., dlsltab(NLSR,NLSW) = 0.   ! runoff per unit of convective, large-scale rain
      real :: dlsflx(NHOR,NLSOUT) = 0.   ! this step's runoff, water, area, evaporation, seepage, recharge, rivers' share
      real :: dlsacc(NHOR,NLSOUT) = 0.   ! their sums since the last output
      integer :: nlsacc = 0
"""
NAMELIST_NAMES = ('nlstore', 'nlsfloor', 'nlssurf', 'dlstau', 'dlsalbw', 'dlsemw', 'dlsz0w', 'dlszw', 'dlsumin',
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
!     product for one storm library and hollow age): the grids of fills, rain means and wetnesses, the chosen
!     surface's convective and large-scale runoff tables (nlssurf), and per point the start of the rain's running
!     means and, for the chosen surface, the store's capacity, its water area at each fill, its seepage under the
!     chosen floor (nlsfloor), the routing shares and the balanced water under that floor. Where the restart held no
!     store, the store starts at its balanced water, and where it held no running means, they start at the record's
!     rain.
!     The store needs the prognostic land.
!
      real zcap(NUGP),zarea(NUGP,NLSF),zseep(NUGP,NLSF),zshare(NUGP,3),zvbal(NUGP),zrc0(NUGP),zrl0(NUGP)
      real zrow(2+NLSSURFS*NLSROW)
      real ztab(NLSR,NLSW)
      real zq
      integer ihead(7)
      integer jhor,jsurf,jf,jo
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
      if (dlstau <= 0.) then
       write(nud,*) '*** the hollows store: the running means need a time constant dlstau above nil ***'
       stop 1
      endif
      if (mypid == NROOT) then
       open(73,file='landstore.dat',form='formatted',status='old')
       read(73,*) ihead
       if (ihead(1) /= LSVERSION .or. ihead(2) /= NUGP .or. ihead(3) /= NLSF .or. ihead(4) /= NLSR                     &
     &     .or. ihead(5) /= NLSW .or. ihead(6) /= NLSFLOORS .or. ihead(7) /= NLSSURFS) then
        write(nud,*) '*** landstore.dat does not match the hollows store patch: ',ihead
        stop 1
       endif
       read(73,*) dlsfill
       read(73,*) dlsrgrid
       read(73,*) dlswgrid
       do jsurf=1,NLSSURFS
        read(73,*) ztab
        if (jsurf == nlssurf) dlsctab(:,:)=ztab(:,:)
       enddo
       do jsurf=1,NLSSURFS
        read(73,*) ztab
        if (jsurf == nlssurf) dlsltab(:,:)=ztab(:,:)
       enddo
       jf=2+NLSROW*(nlssurf-1)
       jo=jf+1+NLSF*(1+NLSFLOORS)
       do jhor=1,NUGP
        read(73,*) zrow
        zrc0(jhor)=zrow(1)
        zrl0(jhor)=zrow(2)
        zcap(jhor)=zrow(jf+1)
        zarea(jhor,:)=zrow(jf+2:jf+1+NLSF)
        zseep(jhor,:)=zrow(jf+2+NLSF*nlsfloor:jf+1+NLSF*(nlsfloor+1))
        zshare(jhor,:)=zrow(jo+1:jo+3)
        zvbal(jhor)=zrow(jo+3+nlsfloor)
       enddo
       close(73)
       write(nud,*) '* Read the hollows store: floor ',nlsfloor,' surface ',nlssurf,' running means (s) ',dlstau
      endif
      call mpbcrn(dlsfill,NLSF)
      call mpbcrn(dlsrgrid,NLSR)
      call mpbcrn(dlswgrid,NLSW)
      call mpbcrn(dlsctab,NLSR*NLSW)
      call mpbcrn(dlsltab,NLSR*NLSW)
      call mpscgp(zcap,dlscap,1)
      call mpscgp(zarea,dlsarea,NLSF)
      call mpscgp(zseep,dlsseep,NLSF)
      call mpscgp(zshare,dlsshare,3)
      call mpscgp(zvbal,dlsvbal,1)
      call mpscgp(zrc0,dlsrc0,1)
      call mpscgp(zrl0,dlsrl0,1)
      do jhor=1,NHOR
       if (dls(jhor) > 0. .and. dlscap(jhor) > 0.) then
        if (nlsread == 0) dlsv(jhor)=dlsvbal(jhor)
        if (nlsrread == 0) then
         dlsrc(jhor)=dlsrc0(jhor)
         dlsrl(jhor)=dlsrl0(jhor)
        endif
        call lscurve(NLSF,dlsfill,dlsarea(jhor,:),dlsseep(jhor,:),dlscap(jhor),dlsv(jhor),dlsfw(jhor),zq)
       else
        dlsv(jhor)=0.
        dlsfw(jhor)=0.
        dlsrc(jhor)=0.
        dlsrl(jhor)=0.
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
!     Terluna: the store's water in this step, between tands and wandr (lsstep): the rain's running means take the
!     step's rain, the bucket keeps what the store does not take, and gives back the store's evaporation, which devap
!     carries; the slow store gains the seepage and the recharge below the roots.
!
      real zdwat,zvn,zfwn,zrcn,zrln,zroff,zseep,zrech,zspill,zriver,zinfl
      integer jhor
!
      do jhor=1,NHOR
       dlsflx(jhor,:)=0.
       if (dls(jhor) > 0. .and. dlscap(jhor) > 0.) then
        call lsstep(deltsec,dprc(jhor),dprl(jhor),dlses(jhor),dlstau,dlsrc(jhor),dlsrl(jhor),dwatc(jhor),dwmax(jhor),   &
     &              dlsfw(jhor),dlsv(jhor),dlscap(jhor),NLSF,dlsfill,dlsarea(jhor,:),dlsseep(jhor,:),dlsshare(jhor,:),    &
     &              NLSR,NLSW,dlsrgrid,dlswgrid,dlsctab,dlsltab,zdwat,zvn,zfwn,zrcn,zrln,zroff,zseep,zrech,zspill,zriver, &
     &              zinfl)
        dwater(jhor)=dwater(jhor)+zdwat
        dlsv(jhor)=zvn
        dlsfw(jhor)=zfwn
        dlsrc(jhor)=zrcn
        dlsrl(jhor)=zrln
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
!     at shares pfill of its capacity (land_store.curve): linear in the square root of the share, nil at nil water.
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

      real function lsshare(kr,kw,prg,pwg,ptab,pr,pw)
!
!     Terluna: the slopes' runoff per unit of rain (land_store.runoff_share): bilinear in the rain's running mean pr
!     (mm/day) and the bucket's relative wetness pw, both held within the table's grids.
!
      integer kr,kw
      real prg(kr),pwg(kw),ptab(kr,kw),pr,pw
      real zr,zw,za,zb
      integer i,j
!
      zr=MIN(MAX(pr,prg(1)),prg(kr))
      zw=MIN(MAX(pw,pwg(1)),pwg(kw))
      do i=1,kr-1
       if (zr <= prg(i+1)) exit
      enddo
      i=MIN(i,kr-1)
      do j=1,kw-1
       if (zw <= pwg(j+1)) exit
      enddo
      j=MIN(j,kw-1)
      za=(zr-prg(i))/(prg(i+1)-prg(i))
      zb=(zw-pwg(j))/(pwg(j+1)-pwg(j))
      lsshare=(1.-za)*(1.-zb)*ptab(i,j)+za*(1.-zb)*ptab(i+1,j)+(1.-za)*zb*ptab(i,j+1)+za*zb*ptab(i+1,j+1)
      return
      end function lsshare

!     =================
!     SUBROUTINE LSSTEP
!     =================

      subroutine lsstep(pdt,pprc,pprl,pes,ptau,prc,prl,pw,pwmax,pfw,pv,pcap,kf,pfill,parea,pseep,pshare,kr,kw,prg,pwg,  &
     &                  pctab,pltab,pdwat,pvn,pfwn,prcn,prln,proff,pseepo,precho,pspill,priver,pinfl)
!
!     Terluna: one step of the store in a cell (land_store.store_step), after its evaporation pes (m/s, upward, at most
!     its water where positive) has joined the land's: the rain's running means prc and prl relax towards the step's
!     rain with the time constant ptau, the rain falls on the store's water area, and the slopes shed runoff on the
!     rest as the storm library's storms (each unit of convective rain through pctab at its running mean, of
!     large-scale rain through pltab at its own, against the bucket's relative wetness); the runoff's routing by the
!     shares pshare (hollows, rivers, recharge below the roots; the rest soaks back), the seepage at the step's start,
!     the spill above the capacity. pdwat is the change to the bucket's dwater (m/s).
!
      integer kf,kr,kw
      real pdt,pprc,pprl,pes,ptau,prc,prl,pw,pwmax,pfw,pv,pcap,pfill(kf),parea(kf),pseep(kf),pshare(3)
      real prg(kr),pwg(kw),pctab(kr,kw),pltab(kr,kw)
      real pdwat,pvn,pfwn,prcn,prln,proff,pseepo,precho,pspill,priver,pinfl
      real lsshare
      real zrel,za,zps,zq0,zv,zfw0,zq
!
      zrel=0.
      if (pwmax > 0.) zrel=MIN(MAX(pw/pwmax,0.),1.)
      za=1.-EXP(-pdt/ptau)
      prcn=prc+za*(pprc-prc)
      prln=prl+za*(pprl-prl)
      zps=pfw*(pprc+pprl)
      proff=(lsshare(kr,kw,prg,pwg,pctab,prcn*8.64e7,zrel)*pprc+lsshare(kr,kw,prg,pwg,pltab,prln*8.64e7,zrel)*pprl)    &
     &      *(1.-pfw)
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
     + f"       if (nlstore == 1) then   ! {LAND_STORE_MARKER}: its water and the rain's running means, if held\n"
       "        ilsread = 0\n"
       "        if (mypid == NROOT) call restart_has_array('dlstorev',ilsread)\n"
       "        call mpbci(ilsread)\n"
       "        if (ilsread == 1) then\n"
       "         call mpgetgp('dlstorev',dlsv   ,NHOR,     1)\n"
       "         call mpgetgp('dlstoreg',dlsslow,NHOR,     1)\n"
       "         nlsread = 1\n"
       "        endif\n"
       "        ilsread = 0\n"
       "        if (mypid == NROOT) call restart_has_array('dlstorec',ilsread)\n"
       "        call mpbci(ilsread)\n"
       "        if (ilsread == 1) then\n"
       "         call mpgetgp('dlstorec',dlsrc  ,NHOR,     1)\n"
       "         call mpgetgp('dlstorel',dlsrl  ,NHOR,     1)\n"
       "         nlsrread = 1\n"
       "        endif\n"
       "       endif\n"),
    (_LANDINI_END, f"      if (nlstore == 1) call lstoreini   ! {LAND_STORE_MARKER}: its tables and start\n" + _LANDINI_END),
    (_RESTART_PUT, _RESTART_PUT
     + f"      if (nlstore == 1) then   ! {LAND_STORE_MARKER}\n"
       "       call mpputgp('dlstorev',dlsv   ,NHOR, 1)\n"
       "       call mpputgp('dlstoreg',dlsslow,NHOR, 1)\n"
       "       call mpputgp('dlstorec',dlsrc  ,NHOR, 1)\n"
       "       call mpputgp('dlstorel',dlsrl  ,NHOR, 1)\n"
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
