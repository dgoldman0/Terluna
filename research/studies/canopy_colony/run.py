"""The canopy colony sized against the seas' oxygen: what the sea does unaided, what each part of the colony must do
for the oxygen to hold, what the colony gets back, and the same for the fresh seas.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.canopy_colony.run      # about a minute

It reads the south sea study's stored estimate (the sea's floor by level, its saturation and its demand), the paused
oxygen run's reading (inputs/south_sea_lunar_mixing_1deg.json, a copy of the south sea study's product on
domain/hydrosphere), this study's seas (seas.py) and currents (currents.py), and the sea oxygen study's air, and
writes results/canopy_colony.json. Nothing here runs a circulation model: the sea is the column of column.py, mixed
at the paused run's sea-wide diffusivity with the share above 300 m that makes the column follow that run's three
years.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from biosphere.ecology import lake as L
from research.studies.canopy_colony import column as CL
from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'canopy_colony.json'
SCHEMA = 'terluna.research.canopy-colony/1'
MIXING = HERE / 'inputs' / 'south_sea_lunar_mixing_1deg.json'
MIXING_SOURCE = 'research/studies/south_sea_ocean/results/lunar_mixing_1deg.json'   # on domain/hydrosphere, 18916584
ESTIMATE = ROOT / 'research/studies/south_sea_ocean/results/circulation_1deg.json'
SEA_OXYGEN = ROOT / 'research/studies/sea_oxygen/results/sea_oxygen.json'
SEAS = HERE / 'results' / 'seas.json'
CURRENTS = HERE / 'results' / 'currents.json'

# --- the sinking matter (the south sea study's correction of 22 November, repeated here until the branches join)
EARTH_FLUX_WATER_C = 7.5                # the water Martin's curve was measured in below 100 m, about 5-10 C (recalled)
SINKING_WATER_C = 28.0                  # the seas between 100 and 500 m
DECAY_Q10 = 2.0                         # the decay's rise for each 10 K (recalled, a usual value for microbes)
EARTH_SINKING_M_DAY = (20.0, 100.0)     # speeds that carry most of Earth's sinking flux (recalled)

# --- the cases
EXPORTS_G_C = dict(poorest=12.0, middling=64.2)          # beside the reference; g C per m2 a year below 100 m
HYPOXIC_MG_L = 2.0                      # the oxygen under which most fish fail (the sea oxygen study's threshold)
EELS_M = (270.0, 600.0)                 # the layer cool enough for the eels' eggs (the south sea study)
EELS_MG_L = 4.0                         # what the south sea study's spawning ground was sized to hold
MARKS = (100.0, 300.0, 1000.0)
UPPER_MIXING = (0.5, 2.0)               # the fitted share halved and doubled

# --- the colony's parts
CATCH = 0.5                             # the share of the flux at the tier's lower edge the feeders take: a target
CATCHES = (0.25, 0.5, 0.75)
CATCH_BELOW_M = 250.0                   # the feeders' tier lies at 150-250 m, in water the air renews
TRUNKS = {'short': (250.0, 500.0), 'long': (250.0, 900.0)}      # the depths over which the trunk lets its water out
CENTRAL_TRUNK = 'long'
TRUNK_MAX_M_YR = 60.0
CAPSULE_DELIVERY = 0.7                  # the share of capsules that release where they are needed (the brief's placeholder)

# --- physics and chemistry of the parts
NUSSELT = 48.0 / 11.0                   # laminar flow in a tube whose wall temperature changes evenly along it
HEAT_DIFFUSIVITY_M2_S = 1.49e-7         # water at 28 C (conductivity 0.61 W/m/K over 4.1 MJ/m3/K, recalled)
TRUNK_LAGS_K = (0.3, 0.1, 0.03)         # how much warmer than the water around it the trunk's water arrives
CHANNEL_BORE_M = 0.05
ENERGY_KJ_PER_MOL_C = 2870.0 / 6.0      # glucose respired (recalled)
C_TO_N, C_TO_P = 106.0 / 16.0, 106.0    # Redfield (recalled)
ELECTRONS_PER_O2_DELIVERED = 4.0        # two peroxides, two electrons each (O2 + 2 e- + 2 H+ -> H2O2)
ELECTRONS_PER_C_FIXED = 4.0
CALCITE_MT_YR = 150.0                   # the held seas' calcium, buried as calcite (the sea chemistry study)
CALCITE_G_PER_MOL = 100.087
REEF_MOL_M2_YR = (10.0, 100.0)          # carbonate laid down per m2 of reef a year (the brief's range, recalled)
SEAWATER_CA_MOL_M3 = 10.28 * 1.022      # 10.28 mmol/kg (Millero et al. 2008, the sea chemistry study's seawater)
FRESH_SEA_CA_MOL_M3 = 0.84              # the held fresh sea (the sea chemistry study)
CAO2_G_PER_MOL = 72.08
FEEDER_EFFICIENCIES = (0.1, 0.5)        # the share of the particles in the water passing a feeder that it takes (a range)
TISSUE_RESPIRATION_A_YEAR = (2.0, 7.0)  # a suspension feeder's own carbon respired a year at 28 C: cold-water corals'
#                                         0.1-0.5% a day at 8 C carried up by the decay's Q10 (recalled)
CANOPY_G_C_M2_YR = (300.0, 1000.0)      # carbon a lit canopy of symbiont-bearing members hands on, per m2 of canopy a
#                                         year (reefs' production, recalled; the long night is not counted)
P_G_PER_MOL = 30.974


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sinking_stretch(water_c=SINKING_WATER_C, gravity=MOON_SURFACE_GRAVITY, decay_q10=None):
    """How many times shorter a sinking particle's fall is than on Earth before the same share of it has decayed:
    Stokes's speed goes as gravity over the water's viscosity, and with decay_q10 the matter also decays faster in the
    warmer water. The south sea study's function of 22 November (circulation.sinking_stretch on domain/hydrosphere)."""
    stretch = STANDARD_GRAVITY / gravity * float(L.viscosity(water_c) / L.viscosity(EARTH_FLUX_WATER_C))
    if decay_q10:
        stretch *= decay_q10 ** ((water_c - EARTH_FLUX_WATER_C) / 10.0)
    return stretch


def sinking_forms():
    return {'earth': 1.0, 'moon': sinking_stretch(), 'moon_warm': sinking_stretch(decay_q10=DECAY_Q10)}


# ---------------------------------------------------------------- the south sea's column

def south_sea_column():
    """The south sea on the paused run's levels and floor, mixed at that run's sea-wide diffusivity of model year 53
    with the share above 300 m that follows its three years. Returns the column, the fit and the run's figures."""
    est = json.loads(ESTIMATE.read_text())['oxygen_estimate']
    run = json.loads(MIXING.read_text())
    mid = np.array(run['level_mid_m'])
    k = np.array([0.0 if v is None else v for v in run['years'][-1]['diffusivity_by_depth_m2_s']])
    k[0] = k[1]                                         # the run keeps none for its top level, which the air holds anyway
    raw = CL.from_level_volumes(mid, est['basin_volume_km3'], k, est['saturation_mg_l'])
    ocean = est['states']['oligotrophic']['ocean']
    study = (np.array(ocean['demand_g_m3_yr']['water']) + np.array(ocean['demand_g_m3_yr']['floor']))[:raw.n]
    used = [y['used_by_depth_mg_l'] for y in run['years']]
    share, err = CL.calibrated_share(raw, study, used)
    kc = raw.k.copy()
    kc[raw.mid < 300.0] *= share
    col = raw.with_k(kc)
    _, _, means = CL.integrate(col, study, len(used), yearly=len(used))
    fit = dict(share_of_the_mean_diffusivity_above_300_m=share, log_error=err,
               run_used_mg_l=[[None if v is None else float(v) for v in y[:col.n]] for y in used],
               column_used_mg_l=[[float(v) for v in m] for m in means])
    temperature = np.array(run['years'][-1]['warmth']['by_depth_c'][:col.n], float)
    return col, fit, dict(export_g_o2=float(ocean['export_g_m2_yr']), study_demand=study, raw_k=raw.k,
                          temperature_c=temperature, new_water=np.array(
                              [np.nan if v is None else v for v in run['years'][-1]['new_water_by_depth'][:col.n]]))


def oxygen_summary(col, spent):
    o2 = col.saturation - spent
    below = col.mid > CL.MARTIN_FROM_M
    lowest = int(np.flatnonzero(below)[np.argmin(o2[below])])
    short = np.flatnonzero(below & (o2 < HYPOXIC_MG_L))
    return dict(by_level_mg_l=[float(x) for x in o2], eels_layer_mg_l=CL.layer_mean(col, o2, EELS_M),
                lowest_mg_l=float(o2[lowest]), lowest_at_m=float(col.mid[lowest]),
                under_2_mg_l_between_m=[float(col.edge[short[0]]), float(col.edge[short[-1] + 1])] if short.size else None)


def placement(col, stretch):
    """Where the demand below 100 m falls: the share of the export still sinking in open water at each depth, and
    the shares spent in each band of depth over the sea's own floor."""
    carried = lambda z: float(CL.flux_share(float(z), stretch) * np.interp(z, col.edge, col.area))   # noqa: E731
    total = carried(CL.MARTIN_FROM_M)
    bands = {}
    for lo, hi in ((100, 250), (250, 400), (400, 850), (850, None)):
        bands[f'{lo}-{hi} m' if hi else f'below {lo} m'] = (carried(lo) - (carried(hi) if hi else 0.0)) / total
    d, _ = CL.demand(col, 1.0, stretch)
    floor = sum(float(CL.flux_share(0.5 * (col.edge[i] + col.edge[i + 1]), stretch)) * (col.area[i] - col.area[i + 1])
                for i in range(col.n) if col.edge[i] >= CL.MARTIN_FROM_M)
    bands['of it on the floor'] = floor / total
    return dict(still_sinking_at_m={str(z): float(CL.flux_share(float(z), stretch)) for z in (150, 200, 250, 300, 400, 500, 850, 2100)},
                share_of_demand=bands)


def unaided(col, export_g_o2, stretch):
    d, _ = CL.demand(col, export_g_o2, stretch)
    _, at, _ = CL.integrate(col, d, max(MARKS), marks=MARKS)
    return {f'{int(m)} years': oxygen_summary(col, at[m]) for m in MARKS}


# ---------------------------------------------------------------- the colony's parts on a column

def colony_terms(col, export_g_o2, stretch, catch, foot_m, capsules='need', trunk_m_yr=0.0):
    """The net demand with the feeders' catch and the capsules' oxygen, and the operator with the trunk. The capsules
    release below the trunk's foot in proportion to the demand left there: 'need' gives all of it, a number that many
    mol of O2 per m2 a year, None gives none. Returns the net demand, the operator and the amounts."""
    d, caught = CL.demand(col, export_g_o2, stretch, catch, CATCH_BELOW_M)
    below = col.mid >= foot_m[1]
    left = float((d * col.volume)[below].sum())                       # g O2 per m2 a year below the trunk's reach
    given = 0.0 if capsules is None else left if capsules == 'need' else min(float(capsules) * CL.G_O2_PER_MOL, left)
    net = d.copy()
    if left > 0:
        net[below] *= 1.0 - given / left
    return net, CL.operator(col, trunk_m_yr, foot_m), dict(caught_g_o2=caught, below_the_foot_g_o2=left, capsules_g_o2=given)


def steady_oxygen(col, export_g_o2, stretch, catch, foot_m, trunk_m_yr, capsules='need'):
    net, op, _ = colony_terms(col, export_g_o2, stretch, catch, foot_m, capsules, trunk_m_yr)
    return col.saturation - CL.steady(col, net, op)


def holds(col, o2, eels_mg_l=None):
    below = col.mid > CL.MARTIN_FROM_M
    ok = bool((o2[below] >= HYPOXIC_MG_L - 1e-9).all())
    return ok and (eels_mg_l is None or CL.layer_mean(col, o2, EELS_M) >= eels_mg_l - 1e-9)


def trunk_needed(col, export_g_o2, stretch, catch, foot_m, eels_mg_l=None):
    """The least water the trunk must send down (m a year over the sea) for the steady sea to keep 2 mg/L at every
    level below 100 m, and with eels_mg_l that much in the eels' layer too, the capsules meeting the demand below
    the trunk's foot. None where TRUNK_MAX_M_YR is not enough."""
    test = lambda w: holds(col, steady_oxygen(col, export_g_o2, stretch, catch, foot_m, w), eels_mg_l)   # noqa: E731
    if test(0.0):
        return 0.0
    if not test(TRUNK_MAX_M_YR):
        return None
    lo, hi = 0.0, TRUNK_MAX_M_YR
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if test(mid) else (mid, hi)
    return hi


def system(col, export_g_o2, stretch, catch=CATCH, trunk=CENTRAL_TRUNK):
    """One sized system: what each part must carry for the steady sea to hold, and the sea it then keeps."""
    foot = TRUNKS[trunk]
    _, _, amounts = colony_terms(col, export_g_o2, stretch, catch, foot)
    w2 = trunk_needed(col, export_g_o2, stretch, catch, foot)
    w4 = trunk_needed(col, export_g_o2, stretch, catch, foot, EELS_MG_L)
    out = dict(catch=catch, trunk=trunk, foot_m=list(foot),
               caught_mol_o2_m2_yr=amounts['caught_g_o2'] / CL.G_O2_PER_MOL,
               capsules_needed_mol_o2_m2_yr=amounts['below_the_foot_g_o2'] / CL.G_O2_PER_MOL,
               trunk_m_yr_for_2_mg_l_everywhere=w2, trunk_m_yr_with_the_eels_layer_at_4_mg_l=w4)
    if w4 is not None:
        out['steady_with_that_trunk'] = oxygen_summary(col, col.saturation - steady_oxygen(col, export_g_o2, stretch, catch, foot, w4))
    return out


def parts_removed(col, export_g_o2, stretch, sized):
    """The sized system run for a thousand years from a saturated sea, whole and with each part taken away."""
    foot, w = TRUNKS[sized['trunk']], sized['trunk_m_yr_with_the_eels_layer_at_4_mg_l']
    cases = (('whole', sized['catch'], w, 'need'), ('without the feeders', 0.0, w, sized['capsules_needed_mol_o2_m2_yr']),
             ('without the trunk', sized['catch'], 0.0, 'need'), ('without the capsules', sized['catch'], w, None),
             ('feeders alone', sized['catch'], 0.0, None))
    out = {}
    for name, catch, flow, capsules in cases:
        net, op, _ = colony_terms(col, export_g_o2, stretch, catch, foot, capsules, flow)
        _, at, _ = CL.integrate(col, net, max(MARKS), marks=MARKS, L=op)
        out[name] = {f'{int(m)} years': {k: v for k, v in oxygen_summary(col, at[m]).items() if k != 'by_level_mg_l'}
                     for m in MARKS[1:]}
    # without the capsules the water below the trunk's foot lives on its stock: how long until a level there is short
    net, op, _ = colony_terms(col, export_g_o2, stretch, sized['catch'], foot, None, w)
    below = col.mid >= foot[1]
    spans = tuple(float(x) for x in range(1000, 40001, 1000))
    _, at, _ = CL.integrate(col, net, spans[-1], marks=spans, L=op, dt=0.25)
    short = [t for t in spans if (col.saturation - at[t][below]).min() < HYPOXIC_MG_L]
    out['without the capsules']['first_millennium_with_under_2_mg_l_below_the_foot'] = short[0] if short else None
    return out


# ---------------------------------------------------------------- what the colony spends and gets

def trunk_engineering(col, flow_m_yr, foot_m, temperature_c, currents):
    """The trunk as thin-walled channels that give their water's warmth to the water rising round them. A channel
    whose water arrives `lag` warmer than its surroundings can carry Nu pi kappa lag over the temperature gradient,
    whatever its bore; the lag times the flow is the heat the trunk leaks downward."""
    top = float(temperature_c[col.mid < 50.0].mean())
    foot = (col.mid >= foot_m[0]) & (col.mid < foot_m[1])
    foot_c = float((temperature_c[foot] * col.volume[foot]).sum() / col.volume[foot].sum())
    length = float(0.5 * (foot_m[0] + min(foot_m[1], col.edge[-1])))
    gradient = (top - foot_c) / length
    rho_c = 4.1e6
    out = dict(surface_c=top, foot_c=foot_c, mean_length_m=length, litres_per_m2_a_day=flow_m_yr * 1000.0 / 365.25, lags={})
    for lag in TRUNK_LAGS_K:
        q = NUSSELT * np.pi * HEAT_DIFFUSIVITY_M2_S * lag / gradient                    # m3/s in one channel
        leak = rho_c * flow_m_yr / CL.YEAR_S * lag                                      # W per m2 of sea
        out['lags'][f'{lag:g} K'] = dict(
            channel_m3_a_day=q * 86400.0, channels_per_km2=flow_m_yr * 1e6 / (q * CL.YEAR_S),
            speed_in_a_5_cm_bore_m_s=q / (np.pi * (CHANNEL_BORE_M / 2) ** 2), heat_leak_w_m2=leak,
            foot_warms_k_a_century_if_nothing_carries_it_off=flow_m_yr * lag * 100.0 / float(col.dz[foot].sum()))
    out['uncooled_pipe_warms_the_foot_k_a_century'] = flow_m_yr * (top - foot_c) * 100.0 / float(col.dz[foot].sum())
    # the work of pushing the water down: it arrives lighter or heavier than the foot's water by its salt
    water = currents['surface_water_cooled_to_the_foot'][f'{foot_m[0]:g}-{foot_m[1]:g} m']
    flow_m_s = flow_m_yr / CL.YEAR_S
    heads = {k: max(-water['excess_density_kg_m3'][k], 0.0) * MOON_SURFACE_GRAVITY * length for k in ('median', 'p5')}
    out['pumping'] = dict(share_of_the_sea_where_it_sinks_by_itself=water['share_of_area_denser'],
                          head_pa=heads, work_j_m2_yr={k: v * flow_m_s * CL.YEAR_S for k, v in heads.items()},
                          head_per_metre_of_trunk_pa={k: v / length for k, v in heads.items()})
    return out


def feeders_net(catch, stretch, currents):
    """The feeding surface that takes `catch` of the sinking matter: a particle that sinks at w spends H / w in a tier
    H thick while the water carries it u H / w past the colony, so a share 1 - exp(-e F u / w) is taken by feeders of
    efficiency e whose surface facing the flow is F m2 for each m2 of sea."""
    past = currents['runs'][currents['central_run']]['water_past_a_drifting_canopy_m_s']['feeders']
    need = -np.log(1.0 - catch)
    out = dict(sinking_m_day=[v / stretch for v in EARTH_SINKING_M_DAY], water_past_the_feeders_m_s=dict(median=past['median'], p10=past['p10']))
    out['surface_times_efficiency_m2_per_m2'] = {
        f'{label} flow, {speed:g} m/day': need * (speed / 86400.0) / u
        for label, u in (('median', past['median']), ('slow', past['p10'])) for speed in out['sinking_m_day']}
    values = list(out['surface_times_efficiency_m2_per_m2'].values())
    out['facing_surface_share_of_the_sea'] = [min(values) / max(FEEDER_EFFICIENCIES), max(values) / min(FEEDER_EFFICIENCIES)]
    return out


def calcium(capsules_mol_o2, sea, all_rivers_m3s, stock_mol_m3):
    """The calcium the capsules take against what the sea's rivers bring, the held seas' calcite shared out by river
    flow, and against the sea's dissolved stock."""
    made = capsules_mol_o2 / CAPSULE_DELIVERY                                   # mol O2 carried by all capsules made
    ca = 2.0 * made
    per_m3 = CALCITE_MT_YR * 1e12 / CALCITE_G_PER_MOL / (all_rivers_m3s * CL.YEAR_S)
    rivers = per_m3 * sea['rivers_m_yr']                                        # mol Ca per m2 a year
    stock = stock_mol_m3 * sea['mean_depth_m']
    beyond = max(ca - rivers, 0.0)
    return dict(capsules_made_mol_o2_m2_yr=made, calcium_mol_m2_yr=ca, peroxide_g_m2_yr=ca * CAO2_G_PER_MOL,
                rivers_calcium_mol_m2_yr=rivers, rivers_calcium_mol_m3=per_m3, share_of_the_rivers_calcium=ca / rivers,
                oxygen_the_rivers_calcium_carries_mol_m2_yr=0.5 * rivers * CAPSULE_DELIVERY,
                dissolved_stock_mol_m2=stock, stock_drawn_percent_a_millennium=100.0 * beyond * 1000.0 / stock,
                canopy_share_of_the_sea_at_reef_rates=[ca / REEF_MOL_M2_YR[1], ca / REEF_MOL_M2_YR[0]],
                light_carbon_forgone_g_m2_yr=made * ELECTRONS_PER_O2_DELIVERED / ELECTRONS_PER_C_FIXED * CL.G_C_PER_MOL)


def loop(col, export_g_o2, stretch, sized, basis, currents, sea, all_rivers_m3s):
    """What passes between the colony and the sea in the sized system, per m2 of sea a year."""
    foot = TRUNKS[sized['trunk']]
    export_c = export_g_o2 / CL.G_O2_PER_MOL / CL.O2_PER_C                      # mol C
    caught_c = sized['caught_mol_o2_m2_yr'] / CL.O2_PER_C
    flux = lambda z: float(CL.flux_share(z, stretch))                           # noqa: E731
    a = lambda z: float(np.interp(z, col.edge, col.area))                       # noqa: E731
    d_none, _ = CL.demand(col, 1.0, stretch)
    d_catch, caught_share = CL.demand(col, 1.0, stretch, sized['catch'], CATCH_BELOW_M)
    reach = (col.mid >= foot[0]) & (col.mid < foot[1])
    below = col.mid >= foot[1]
    at_ref = a(CL.MARTIN_FROM_M)                                                # the open water at 100 m, where the export is counted
    where = lambda d: dict(above_the_trunk=float((d * col.volume)[col.mid < foot[0]].sum()) / at_ref,  # noqa: E731
                           in_the_trunks_reach=float((d * col.volume)[reach].sum()) / at_ref,
                           below_it=float((d * col.volume)[below].sum()) / at_ref)
    flow = sized['trunk_m_yr_with_the_eels_layer_at_4_mg_l']
    engineering = trunk_engineering(col, flow, foot, basis['temperature_c'], currents)
    pump = engineering['pumping']['work_j_m2_yr']
    capsules = calcium(sized['capsules_needed_mol_o2_m2_yr'], sea, all_rivers_m3s, SEAWATER_CA_MOL_M3 if sea['salt'] else FRESH_SEA_CA_MOL_M3)
    lost = lambda d: where(d)['below_it']                                       # noqa: E731
    export_p_mg = 1e3 * export_c / C_TO_P * P_G_PER_MOL                         # mg P per m2 a year leaving the lit layer
    unreturned = lambda d: where(d)['in_the_trunks_reach'] + where(d)['below_it']   # noqa: E731
    phosphorus = dict(
        export_mg_m2_yr=export_p_mg,
        lost_below_the_trunks_reach_mg_m2_yr=export_p_mg * lost(d_catch),
        without_the_colony_lost_below_250_m_mg_m2_yr=export_p_mg * unreturned(d_none),
        rivers_m_yr=sea['rivers_m_yr'],
        river_phosphorus_that_replaces_the_loss_mg_m3=export_p_mg * lost(d_catch) / sea['rivers_m_yr'],
        export_g_c_m2_yr_held_by_each_mg_m3_of_river_phosphorus=export_c * CL.G_C_PER_MOL
        / (export_p_mg * lost(d_catch) / sea['rivers_m_yr']))
    steady = sized['steady_with_that_trunk']['by_level_mg_l']
    spent_at_foot = col.saturation - float(np.mean(np.array(steady)[reach]))
    return dict(
        export=dict(g_c_m2_yr=export_c * CL.G_C_PER_MOL, mol_o2_m2_yr=export_g_o2 / CL.G_O2_PER_MOL,
                    nitrogen_mmol_m2_yr=1e3 * export_c / C_TO_N, phosphorus_mmol_m2_yr=1e3 * export_c / C_TO_P),
        catch=dict(share_of_the_export=caught_share / a(CL.MARTIN_FROM_M), g_c_m2_yr=caught_c * CL.G_C_PER_MOL,
                   energy_kj_m2_yr=caught_c * ENERGY_KJ_PER_MOL_C, nitrogen_mmol_m2_yr=1e3 * caught_c / C_TO_N,
                   phosphorus_mmol_m2_yr=1e3 * caught_c / C_TO_P,
                   oxygen_to_respire_it_mol_m2_yr=sized['caught_mol_o2_m2_yr'],
                   carried_up_as_food_l_m2_yr_at_10_g_c_per_l=caught_c * CL.G_C_PER_MOL / 10.0,
                   living_tissue_it_keeps_g_c_m2=[caught_c * CL.G_C_PER_MOL / x for x in TISSUE_RESPIRATION_A_YEAR[::-1]],
                   net=feeders_net(sized['catch'], stretch, currents)),
        trunk=dict(m_yr=flow, oxygen_carried_mol_m2_yr=flow * col.saturation / CL.G_O2_PER_MOL,
                   oxygen_given_up_at_the_foot_mol_m2_yr=flow * spent_at_foot / CL.G_O2_PER_MOL,
                   pumping_share_of_the_catchs_energy={k: v / 1e3 / (caught_c * ENERGY_KJ_PER_MOL_C) for k, v in pump.items()},
                   **engineering),
        capsules=dict(needed_mol_o2_m2_yr=sized['capsules_needed_mol_o2_m2_yr'], **capsules,
                      light_carbon_forgone_share_of_the_catch=capsules['light_carbon_forgone_g_m2_yr'] / (caught_c * CL.G_C_PER_MOL)),
        canopy=dict(share_of_the_sea=capsules['canopy_share_of_the_sea_at_reef_rates'],
                    carbon_handed_on_g_m2_yr=[capsules['canopy_share_of_the_sea_at_reef_rates'][0] * CANOPY_G_C_M2_YR[0],
                                              capsules['canopy_share_of_the_sea_at_reef_rates'][1] * CANOPY_G_C_M2_YR[1]]),
        nutrients=dict(regenerated_share_of_the_export=dict(without_the_colony=where(d_none), with_the_feeders=where(d_catch),
                                                            taken_by_the_feeders=caught_share / a(CL.MARTIN_FROM_M)),
                       still_sinking_at_the_tier=flux(CATCH_BELOW_M) * a(CATCH_BELOW_M), phosphorus=phosphorus))


# ---------------------------------------------------------------- the fresh seas

def fresh_saturation():
    setting = json.loads(SEA_OXYGEN.read_text())['setting']
    return float(L.oxygen_saturation(setting['seas_median_c'], setting['design_pressure_pa'] / 101325.0,
                                     setting['salinity_g_kg'], setting['o2_fraction'])), setting['seas_median_c']


def one_sea(name, sea, salt_col, export_g_o2, forms, saturation, all_rivers_m3s):
    """A fresh sea as still as the salt sea: its own floor on the salt sea's levels and mixing, fresh water's
    saturation. Unaided after 300 and 1,000 years, and the system that holds it."""
    table = sea['below']
    col = CL.from_areas(salt_col.mid, [b['depth_m'] for b in table] + [sea['deepest_m']],
                        [b['area_km2'] for b in table] + [0.0], salt_col.k, saturation)
    out = dict(mean_depth_m=sea['mean_depth_m'], levels=col.n, share_of_floor_above_250_m=1.0 - float(np.interp(250.0, col.edge, col.area)),
               share_of_floor_above_400_m=1.0 - float(np.interp(400.0, col.edge, col.area)), cases={})
    for form in ('moon', 'moon_warm'):
        runs = unaided(col, export_g_o2, forms[form])
        sized = system(col, export_g_o2, forms[form])
        entry = dict(unaided={k: {kk: vv for kk, vv in v.items() if kk != 'by_level_mg_l'} for k, v in runs.items() if k != '100 years'},
                     system={k: v for k, v in sized.items() if k != 'steady_with_that_trunk'})
        entry['calcium'] = calcium(sized['capsules_needed_mol_o2_m2_yr'], sea, all_rivers_m3s, FRESH_SEA_CA_MOL_M3)
        out['cases'][form] = entry
    return out


# ---------------------------------------------------------------- the product

def build(log=print):
    started = time.time()
    forms = sinking_forms()
    col, fit, basis = south_sea_column()
    reference = basis['export_g_o2']
    exports = {'nutrient-poor': reference, **{k: v * CL.O2_PER_C / CL.G_C_PER_MOL * CL.G_O2_PER_MOL for k, v in EXPORTS_G_C.items()}}
    seas = json.loads(SEAS.read_text())['seas']
    currents = json.loads(CURRENTS.read_text())
    currents['central_run'] = next(k for k in currents['runs'] if 'own mixing' in k)
    all_rivers = sum(s['rivers_m3s'] for s in seas.values())
    log(f'the column follows the paused run with {fit["share_of_the_mean_diffusivity_above_300_m"]:.3f} of the mean '
        f'diffusivity above 300 m (log error {fit["log_error"]:.2f})')

    south = dict(
        column=dict(level_mid_m=[float(x) for x in col.mid], level_edge_m=[float(x) for x in col.edge],
                    open_area_at_each_edge=[float(x) for x in col.area], mean_depth_m=float(col.volume.sum()),
                    diffusivity_m2_s=[float(x) for x in col.k], run_diffusivity_m2_s=[float(x) for x in basis['raw_k']],
                    saturation_mg_l=col.saturation, temperature_c=[float(x) for x in basis['temperature_c']], fit=fit),
        placement={form: placement(col, s) for form, s in forms.items()},
        unaided={form: {name: unaided(col, e, s) for name, e in exports.items()} for form, s in forms.items()})
    sensitivity = {}
    for factor in UPPER_MIXING:
        k = col.k.copy()
        k[col.mid < 300.0] *= factor
        runs = unaided(col.with_k(k), reference, forms['moon'])
        sensitivity[f'upper mixing times {factor:g}'] = {m: {kk: vv for kk, vv in v.items() if kk != 'by_level_mg_l'} for m, v in runs.items()}
    south['unaided_sensitivity_moon_nutrient_poor'] = sensitivity
    log(f'unaided cases done, {time.time() - started:.0f} s')

    systems = {}
    for form, s in forms.items():
        systems[form] = {}
        for name, e in exports.items():
            systems[form][name] = {f'catch {c:g}, {t} trunk': system(col, e, s, c, t) for c in CATCHES for t in TRUNKS}
            systems[form][name]['no catch, long trunk'] = system(col, e, s, 0.0, 'long')
    south['systems'] = systems
    central = system(col, reference, forms['moon'])
    south['central'] = dict(sinking='moon', export='nutrient-poor', sized=central,
                            thousand_years=parts_removed(col, reference, forms['moon'], central),
                            loop=loop(col, reference, forms['moon'], central, basis, currents, seas['south_pole_aitken'], all_rivers))
    south['central_at_other_exports'] = {
        name: {k: v for k, v in system(col, e, forms['moon']).items() if k != 'steady_with_that_trunk'} for name, e in exports.items()}
    log(f'systems sized, {time.time() - started:.0f} s')

    saturation, fresh_c = fresh_saturation()
    fresh = dict(saturation_mg_l=saturation, at_c=fresh_c,
                 seas={name: one_sea(name, sea, col, reference, forms, saturation, all_rivers)
                       for name, sea in seas.items() if not sea['salt']})
    return dict(
        schema=SCHEMA,
        producer=dict(module='research.studies.canopy_colony.run',
                      code_sha256=hashlib.sha256(Path(__file__).read_bytes() + (HERE / 'column.py').read_bytes()).hexdigest()),
        inputs={str(p.relative_to(ROOT)): sha256(p) for p in (MIXING, ESTIMATE, SEA_OXYGEN, SEAS, CURRENTS)},
        input_copies={str(MIXING.relative_to(ROOT)): dict(copy_of=MIXING_SOURCE, branch='domain/hydrosphere', commit='18916584')},
        evidence=('A sizing on a one-dimensional column of each sea, mixed at the sea-wide diffusivity of three model '
                  'years of one circulation run of the south sea at lunar gravity, with the share above 300 m fitted to '
                  'that run. The sinking matter follows Earth\'s curve shortened by Stokes\'s law; the feeders\' catch, '
                  'the capsules\' chemistry and their delivery are targets and placeholders that no organism or '
                  'experiment has shown. The fresh seas borrow the salt sea\'s mixing. Not validated.'),
        reading_rule=('south_sea.unaided: oxygen (mg/L) by level, in the eels\' layer and at its lowest below 100 m, '
                      'after 100, 300 and 1,000 years from a saturated sea, for each form of the sinking (Earth\'s '
                      'curve, the Moon\'s slower sinking, and that with the warm water\'s faster decay) and each '
                      'export. south_sea.systems: for each catch and trunk, the oxygen the capsules must release '
                      'below the trunk\'s foot and the water the trunk must send down for the steady sea to keep '
                      '2 mg/L at every level, and also 4 mg/L in the eels\' layer; None where no trunk of up to '
                      f'{TRUNK_MAX_M_YR:g} m a year does. south_sea.central: the central system, the sea after a '
                      'thousand years with each part removed, and what passes between the colony and the sea (loop). '
                      'fresh_seas: each fresh sea on its own floor with the salt sea\'s mixing, the upper bound of '
                      'its need. Amounts are per m2 of the sea\'s surface a year unless said.'),
        setting=dict(sinking=dict(stretch=forms, water_c=SINKING_WATER_C, earth_curve_water_c=EARTH_FLUX_WATER_C,
                                  decay_q10=DECAY_Q10, martin_b=CL.MARTIN_B, from_m=CL.MARTIN_FROM_M),
                     exports_g_o2_m2_yr=exports, exports_g_c_m2_yr={k: v / CL.G_O2_PER_MOL / CL.O2_PER_C * CL.G_C_PER_MOL for k, v in exports.items()},
                     hypoxic_mg_l=HYPOXIC_MG_L, eels_layer_m=list(EELS_M), eels_mg_l=EELS_MG_L,
                     catch=CATCH, catch_below_m=CATCH_BELOW_M, trunks={k: list(v) for k, v in TRUNKS.items()},
                     capsule_delivery=CAPSULE_DELIVERY, calcite_mt_yr=CALCITE_MT_YR, reef_mol_m2_yr=list(REEF_MOL_M2_YR),
                     all_rivers_into_the_seas_m3s=all_rivers),
        south_sea=south, fresh_seas=fresh, seconds=round(time.time() - started))


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    out = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print(f'wrote {OUT} in {out["seconds"]} s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
