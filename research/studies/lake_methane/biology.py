"""The methane-eating biology: the methanotrophs' kinetics and the anaerobic pathways' limits, from Earth's lakes and
cultures (sources.json), as design variables within those ranges, and what the Moon's rivers and lakes supply to them.

The Open Moon's organisms are engineered: Earth's methanotrophs give the reference points, and a designed community
may sit anywhere in the ranges measured on Earth, within physical and chemical limits. The aerobic sets are lanes for
the tracer (tracer.Lane's fields). The anaerobic pathways are each held to the supply of their electron acceptor:
iron(III) and manganese(IV) oxides settled from the rivers onto the floors, nitrate carried down from the oxic water,
and the lakes' sulphate; and to the rates measured where they work on Earth.
"""
from __future__ import annotations

import json
from pathlib import Path

from research.studies.lake_methane import tracer as T

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CHEMISTRY = ROOT / 'research' / 'studies' / 'sea_chemistry' / 'results' / 'sea_chemistry.json'
LAKES = ROOT / 'research' / 'studies' / 'lake_oxygen' / 'results' / 'lakes.json'
SEAS = ROOT / 'research' / 'studies' / 'sea_oxygen' / 'results' / 'seas.json'
DAY, YEAR = T.DAY, T.YEAR
CELL_C_PG = 0.42                          # pg C a methanotroph cell (Zimmermann et al. 2021: 0.26-0.57)
UM_C_PER_CELL_ML = CELL_C_PG * 1e-12 * 1e6 / 12.011 * 1e6 / 1e3   # uM of carbon for one cell per mL

# Aerobic methanotrophs. Earth's lake communities: growth from the doubling times measured at oxyclines (1.8-3.2 days
# at Rotsee and Zug, Oswald et al. 2015, 2016) carried to 20 C at the uptake's 0.49 eV, the methane half-saturation of
# lake communities (0.7-15 uM, Mayr et al. 2020a; 1.7 uM fitted at Rotsee's overturn, Zimmermann et al. 2021; Earth's
# lakes' mean 6.9 uM, Martinez-Cruz et al. 2015), oxygen's (sub-micromolar to 20 uM: Blees et al. 2014; Ren et al.
# 1997; Lidstrom & Somers 1984), the carbon conversion efficiency (0.13-0.16 in the field, 0.19-0.7 in cultures,
# Zimmermann et al. 2021) and a seed stock of 1e4 cells per mL (Rotsee's mixed layer before its overturn, 1-2e4,
# Mayr et al. 2020a; Zimmermann et al. 2021). The loss rate (grazing, lysis, settling) is a design guess.
EARTH = dict(mu20_per_h=0.5 / 24.0, yield_c=0.3, k_ch4_um=5.0, k_o2_um=5.0, loss_per_day=0.1,
             seed_c_um=1e4 * UM_C_PER_CELL_ML, activation_ev=0.49)
# The ends of Earth's ranges a design can reach for or must allow for.
FAST = dict(mu20_per_h=2.0 / 24.0, k_ch4_um=1.7, k_o2_um=1.0, seed_c_um=1e5 * UM_C_PER_CELL_ML)
POOR = dict(mu20_per_h=0.16 / 24.0, k_ch4_um=15.0, k_o2_um=20.0, seed_c_um=1e3 * UM_C_PER_CELL_ML, loss_per_day=0.3)

# (a) the oxycline: each kinetic parameter at the ends of Earth's ranges, the others central
KINETICS = {
    'central': {},
    'k_ch4_1.7': dict(k_ch4_um=1.7), 'k_ch4_15': dict(k_ch4_um=15.0),
    'k_o2_1': dict(k_o2_um=1.0), 'k_o2_20': dict(k_o2_um=20.0),
    'growth_0.16': dict(mu20_per_h=0.16 / 24.0), 'growth_2': dict(mu20_per_h=2.0 / 24.0),
    'yield_0.13': dict(yield_c=0.13), 'yield_0.6': dict(yield_c=0.6),
    'loss_0.03': dict(loss_per_day=0.03), 'loss_0.3': dict(loss_per_day=0.3),
    'activation_0.43': dict(activation_ev=0.43), 'activation_0.55': dict(activation_ev=0.55),
    'designed': FAST, 'poor': POOR, 'none': dict(aerobic=False),
}
# (b) the overturn: the standing stock the water keeps (cells per mL) against the response time (1 / growth)
STOCKS_CELLS_ML = (1e3, 1e4, 1e5, 1e6)
GROWTH_PER_DAY = (0.16, 0.5, 2.0)
SHARES = dict(low=0.3, central=0.5, high=0.8)   # lake_anoxia's methane share of the anoxic mineralization


def lane_sets(column, share='central'):
    """The tracer's lanes for one oxygen column: (a) the kinetic sets, (b) the stock and growth grid, each with the
    central methane share, and the central biology at the low and high shares."""
    lanes = []
    for name, kw in KINETICS.items():
        lanes.append(T.Lane(column, share=SHARES[share], label=f'kinetics/{name}', **dict(EARTH, **kw)))
    for cells in STOCKS_CELLS_ML:
        for mu in GROWTH_PER_DAY:
            kw = dict(EARTH, seed_c_um=cells * UM_C_PER_CELL_ML, mu20_per_h=mu / 24.0)
            lanes.append(T.Lane(column, share=SHARES[share], label=f'overturn/{cells:g}/{mu:g}', **kw))
    for s in ('low', 'high'):
        lanes.append(T.Lane(column, share=SHARES[s], label=f'share/{s}', **EARTH))
    return lanes


# --- the electron acceptors' supply --------------------------------------------------------------------------------

M_FE, M_MN = 55.845, 54.938
FE_PER_CH4, MN_PER_CH4, NO3_PER_CH4, SO4_PER_CH4 = 8.0, 4.0, 1.6, 1.0   # Fe(OH)3, MnO2, nitrate to N2 by nitrite,
#                                                                          sulphate
# Iron: the developed soils' iron, all ferric once weathered (the carbon cycle's 2 FeO + 1/2 O2 -> Fe2O3), from the
# highland soil's 5.1 wt% FeO to the high-titanium mare soil's 15.3 (Lunar Sourcebook, resources' literature), with
# the reactive share of Earth's river particulates, 0.43 of the total (Poulton & Raiswell 2002: 1.94 of 4.49 wt%).
FEO_WT = dict(highland=5.1, mare=15.3)
REACTIVE_SHARE = 0.43
# Manganese: lunar soils hold 0.07-0.21 wt% MnO (highland to mare); taken reactive at the same share.
MNO_WT = dict(highland=0.07, mare=0.21)
# Nitrate: Lake Superior's deep water, 26 uM today against about 5 a century ago (Sterner et al. 2007), the lakes'
# design range; nitrite-dependent oxidation at Lake Constance's measured actual rates, 31-437 umol m-2 d-1, and
# potential 660-4,890 (Deutzmann et al. 2014).
NITRATE_UM = (5.0, 26.0)
NDAMO_ACTUAL_UMOL_M2_D = (31.0, 437.0)
NDAMO_POTENTIAL_UMOL_M2_D = (660.0, 4890.0)
# Sulphate: the lakes' 1.4 uM (sea chemistry's export, as lake_anoxia takes it), renewed by the throughflow; the held
# fresh seas none (sea chemistry's managed sinks take it out of the rivers: 0.0 mmol/kg in its held run).
SULPHATE_UM = 1.4
SEAS_SULPHATE_UM = 0.0
THROUGHFLOW_M_YR = 1.8                     # lake_anoxia: the lakes pass 1.8 m of water a year
# Lake Matano's measured iron arithmetic (Crowe et al. 2011): iron burial 8.5 mmol m-2 d-1, five times the methane
# flux of 1.6, so that iron could carry at most 66% of the oxidation; anaerobic oxidation near the pycnocline at
# 0.5 uM a day, each iron atom recycled eighty times.
MATANO = dict(fe_burial_mmol_m2_d=8.5, ch4_flux_mmol_m2_d=1.6, fe_share_max=0.66, aom_um_d=0.5, recycling=80)


def acceptor_supply(st_chem=None, lake_area_km2=None, waters='lakes', seas_area_km2=None):
    """The acceptors each m2 of lake (or sea) floor receives a year and the methane each could oxidise (mol CH4 m-2
    a-1), at the land's latest and earliest erosion (sea chemistry) and the highland and mare soils. The lakes keep
    what of their rivers' particles does not pass them; the seas take the particles of the runoff that skips the
    lakes and those the lakes pass on. The held seas carry no sulphate (sea chemistry's managed sinks)."""
    chem = st_chem or json.loads(CHEMISTRY.read_text())
    s = chem['setting']
    erosion = sorted(s['erosion_gt_yr'])                                   # Gt a year, latest and earliest
    into = s['lake_share']
    passed = sorted(s['particle_pass'])
    if waters == 'lakes':
        area = (lake_area_km2 or json.loads(LAKES.read_text())['lakes']['lake_area_km2']) * 1e6
        reach = [into * (1.0 - p) for p in passed[::-1]]                   # the lakes keep what does not pass
    else:
        area = (seas_area_km2 or json.loads(SEAS.read_text())['seas']['area_km2']) * 1e6
        reach = [(1.0 - into) + into * p for p in passed]
    sediment = [e * 1e12 * r / area for e, r in zip(erosion, reach)]          # kg m-2 a-1
    out = dict(waters=waters, sediment_kg_m2_yr=sediment, erosion_gt_yr=erosion, iron={}, manganese={})
    for soil, feo in FEO_WT.items():
        fe = feo / 100.0 * M_FE / (M_FE + 15.999) * REACTIVE_SHARE                     # reactive Fe, kg/kg
        mn = MNO_WT[soil] / 100.0 * M_MN / (M_MN + 15.999) * REACTIVE_SHARE
        out['iron'][soil] = dict(reactive_wt_percent=100.0 * fe,
                                 fe_mol_m2_yr=[x * fe * 1000.0 / M_FE for x in sediment],
                                 ch4_mol_m2_yr=[x * fe * 1000.0 / M_FE / FE_PER_CH4 for x in sediment])
        out['manganese'][soil] = dict(reactive_wt_percent=100.0 * mn,
                                      ch4_mol_m2_yr=[x * mn * 1000.0 / M_MN / MN_PER_CH4 for x in sediment])
    so4 = SULPHATE_UM if waters == 'lakes' else SEAS_SULPHATE_UM
    out['sulphate'] = dict(um=so4, ch4_mol_m2_yr=so4 * 1e-3 * THROUGHFLOW_M_YR / SO4_PER_CH4 if waters == 'lakes'
                           else 0.0)
    out['nitrate'] = dict(um=list(NITRATE_UM),
                          ndamo_actual_ch4_mol_m2_yr=[x * 1e-6 * 365.25 for x in NDAMO_ACTUAL_UMOL_M2_D],
                          ndamo_potential_ch4_mol_m2_yr=[x * 1e-6 * 365.25 for x in NDAMO_POTENTIAL_UMOL_M2_D])
    out['matano'] = MATANO
    return out


def aom_lanes(column, supply, share='central'):
    """Lanes adding each anaerobic pathway to Earth's aerobic community, at its acceptor's supply: the floors' iron and
    manganese (at the latest and earliest erosion, highland soil, and the earliest with mare soil), nitrate-coupled
    oxidation in anoxic water at Lake Constance's measured actual rates, sulphate at its renewal, and all together."""
    def floor(x):
        return dict(aom_floor_mol_m2_s=x / YEAR)

    def water(x):
        return dict(aom_water_per_s=1.0 / DAY, aom_water_cap_mol_m2_s=x / YEAR)

    fe_h, fe_m = supply['iron']['highland']['ch4_mol_m2_yr'], supply['iron']['mare']['ch4_mol_m2_yr']
    mn_h, mn_m = supply['manganese']['highland']['ch4_mol_m2_yr'], supply['manganese']['mare']['ch4_mol_m2_yr']
    nd = supply['nitrate']['ndamo_actual_ch4_mol_m2_yr']
    so4 = supply['sulphate']['ch4_mol_m2_yr']
    if supply.get('waters') == 'seas':
        sets = {'all/central': dict(floor(fe_h[1] + mn_h[1]), **water((nd[0] * nd[1]) ** 0.5))}
        return [T.Lane(column, share=SHARES[share], label=f'aom/{k}', **dict(EARTH, **kw)) for k, kw in sets.items()]
    sets = {
        'iron/latest_highland': floor(fe_h[0] + mn_h[0]), 'iron/earliest_highland': floor(fe_h[1] + mn_h[1]),
        'iron/earliest_mare': floor(fe_m[1] + mn_m[1]),
        'nitrate/low': water(nd[0]), 'nitrate/high': water(nd[1]),
        'sulphate': water(so4),
        'all/central': dict(floor(fe_h[1] + mn_h[1]), **water((nd[0] * nd[1]) ** 0.5 + so4)),
        'all/high': dict(floor(fe_m[1] + mn_m[1]), **water(nd[1] + so4)),
    }
    return [T.Lane(column, share=SHARES[share], label=f'aom/{k}', **dict(EARTH, **kw)) for k, kw in sets.items()]
