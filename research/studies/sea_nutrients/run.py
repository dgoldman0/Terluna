"""The seas' phosphorus budget on the design surface: what enters each sea, what open water that supply sustains,
the oxygen its decay takes below the water the air renews, and the supply each sea can take.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.sea_nutrients.run      # about two minutes

It reads this study's loads (loads.py: the land's phosphorus through the design surface's lakes), the phosphorus
study's product (resources/results/sea_phosphorus.json: the land's loss, the dust's range and people's losses), and
the canopy colony study's seas and column (the same column of each sea, with its mixing and the Moon's sinking), and
writes results/sea_nutrients.json.

The budget's rule. A sea without an outlet that is layered as the south sea's circulation run has it keeps what it
receives. Phosphorus entering the lit layer is taken up, sinks as matter and is set free again where the matter
decays: above about 250 m the mixing brings it back within a year or so, and below it comes back over centuries or
longer. So the lit layer strips itself until the phosphorus its sinking matter loses below 250 m each year equals
what enters, and from then on

    export = supply / (the share of the export that passes 250 m), and
    the oxygen spent below 250 m = 106 mol for each mol of phosphorus supplied,

whatever the export is. The second line is the budget's centre: in the water the air does not renew, the oxygen
demand is set by the phosphorus supply and by nothing else, until that water has gathered enough phosphorus to
send it back up, which takes tens of thousands of years at these supplies.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from research.studies.canopy_colony import column as CL
from research.studies.canopy_colony import run as CC

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'sea_nutrients.json'
SCHEMA = 'terluna.research.sea-nutrients/1'
LOADS = HERE / 'results' / 'loads.json'
PHOSPHORUS = ROOT / 'resources' / 'results' / 'sea_phosphorus.json'
SEAS = CC.SEAS

P_G_PER_MOL = 30.974
C_TO_P = 106.0                              # Redfield, mol of carbon for each mol of phosphorus (recalled)
P_MG_PER_G_C = 1e3 * P_G_PER_MOL / (C_TO_P * CL.G_C_PER_MOL)        # 24.3 mg of phosphorus in each g of carbon sunk
O2_MMOL_PER_MG_P = C_TO_P * CL.O2_PER_C / P_G_PER_MOL               # 3.42 mmol of oxygen for each mg set free again
VENTILATED_TO_M = 250.0                     # the water above returns what decays in it within a year or so
HORIZONS = (300.0, 1000.0, 3000.0)
REACTIVE_PARTICLES = (0.1, 0.3)             # the share of river particles' phosphorus that salt water frees (Compton et al.
#                                             2000 in Paytan & McLaughlin 2007, the phosphorus study's range); none in fresh
PEOPLE = (4.0e9, 8.0e9)
SHALLOWS_M = (20, 50)
DELIVERED_MG_M3 = (1.0, 10.0, 70.0)         # phosphorus in the water the seas are filled with; 70 is Earth's ocean mean


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sea_column(name, sea, salt_col, fresh_saturation):
    if sea['salt']:
        return salt_col
    table = sea['below']
    return CL.from_areas(salt_col.mid, [b['depth_m'] for b in table] + [sea['deepest_m']],
                         [b['area_km2'] for b in table] + [0.0], salt_col.k, fresh_saturation)


def passing_shares(col, stretch):
    """Of the phosphorus in the matter leaving 100 m: the share that decays below 250 m with no colony, and with the
    canopy colony's central system the share its feeders take, the share that decays in its trunk's reach, which the
    trunk's rising water returns, and the share that decays below the reach."""
    carried = lambda z: float(CL.flux_share(float(z), stretch) * np.interp(z, col.edge, col.area))    # noqa: E731
    at_ref = carried(CL.MARTIN_FROM_M)
    foot = CC.TRUNKS[CC.CENTRAL_TRUNK][1]
    below_250 = carried(VENTILATED_TO_M) / at_ref
    below_foot = (carried(foot) if foot < col.edge[-1] else 0.0) / at_ref
    return dict(below_250_m=below_250, taken_by_the_feeders=CC.CATCH * below_250,
                returned_by_the_trunk=(1.0 - CC.CATCH) * (below_250 - below_foot),
                below_the_trunks_reach=(1.0 - CC.CATCH) * below_foot)


def lowest_oxygen(col, export_g_o2, stretch, years):
    d, _ = CL.demand(col, export_g_o2, stretch)
    spent, _, _ = CL.integrate(col, d, years, dt=0.25)
    return float((col.saturation - spent)[col.mid > CL.MARTIN_FROM_M].min())


def safe_export(col, stretch, years, top_g_o2=400.0):
    """The largest export (g of oxygen demand per m2 a year) with which the unaided sea keeps 2 mg/L at every level
    below 100 m through `years` from a saturated start."""
    lo, hi = 0.0, top_g_o2
    if lowest_oxygen(col, hi, stretch, years) >= CC.HYPOXIC_MG_L:
        return hi
    for _ in range(22):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if lowest_oxygen(col, mid, stretch, years) >= CC.HYPOXIC_MG_L else (lo, mid)
    return lo


def supply_cases(load, land, dust, salt):
    """What enters a sea's water each year, mg of phosphorus per m2 of sea, low and high."""
    area = load['area_m2']
    per_m2 = lambda share, kg: [1e6 * share * kg[0] / area, 1e6 * share * kg[1] / area]     # noqa: E731
    dissolved = per_m2(load['share_of_the_lands_dissolved_loss'] + load['dissolved_share_from_other_seas'], land['dissolved_kg_yr'])
    particles = [1e6 * load['share_of_the_lands_particles']['low'] * land['particulate_kg_yr'][0] / area,
                 1e6 * load['share_of_the_lands_particles']['high'] * land['particulate_kg_yr'][1] / area]
    freed = [particles[0] * REACTIVE_PARTICLES[0], particles[1] * REACTIVE_PARTICLES[1]] if salt else [0.0, 0.0]
    rivers = [dissolved[0] + freed[0], dissolved[1] + freed[1]]
    return dict(rivers_dissolved=dissolved, rivers_particles=particles, freed_from_particles=freed,
                dust=[1e3 * dust[0], 1e3 * dust[1]], rivers_alone=rivers,
                rivers_and_dust=[rivers[0] + 1e3 * dust[0], rivers[1] + 1e3 * dust[1]])


def what_a_supply_does(supply_mg, shares, safe):
    """For a supply (mg P per m2 a year): the export it sustains with no colony and with the canopy colony, the
    oxygen its phosphorus takes below 250 m each year, and the longest horizon through which the unaided sea keeps
    2 mg/L everywhere at that export."""
    export = supply_mg / (shares['below_250_m'] * P_MG_PER_G_C)
    held = [h for h in HORIZONS if export <= safe[f'{int(h)} years']['export_g_c_m2_yr'] + 1e-12]
    return dict(supply_mg_m2_yr=supply_mg, export_g_c_m2_yr=export,
                export_with_the_colony_g_c_m2_yr=supply_mg / (shares['below_the_trunks_reach'] * P_MG_PER_G_C)
                if shares['below_the_trunks_reach'] > 0 else None,
                oxygen_spent_below_250_m_mol_m2_yr=supply_mg * O2_MMOL_PER_MG_P / 1e3,
                unaided_sea_holds_through_years=max(held) if held else 0.0)


def one_sea(name, sea, load, col, stretch, land, dust):
    shares = passing_shares(col, stretch)
    safe = {}
    for years in HORIZONS:
        e = safe_export(col, stretch, years)
        e_c = e / CL.G_O2_PER_MOL / CL.O2_PER_C * CL.G_C_PER_MOL
        safe[f'{int(years)} years'] = dict(export_g_c_m2_yr=e_c, supply_mg_m2_yr=e_c * shares['below_250_m'] * P_MG_PER_G_C)
    supply = supply_cases(load, land, dust, sea['salt'])
    does = {case: {end: what_a_supply_does(supply[case][i], shares, safe) for i, end in enumerate(('low', 'high'))}
            for case in ('rivers_alone', 'rivers_and_dust')}
    area = sea['area_km2'] * 1e6
    deeper = {b['depth_m']: b['area_km2'] for b in sea['below']}
    shallows = {f'under {z} m': sea['area_km2'] - deeper.get(z, 0.0) for z in SHALLOWS_M}
    river_kg = [load['share_of_the_lands_dissolved_loss'] * land['dissolved_kg_yr'][i] for i in (0, 1)]
    particle_kg = [load['share_of_the_lands_particles'][c] * land['particulate_kg_yr'][i] for i, c in ((0, 'low'), (1, 'high'))]
    mouths = dict(dissolved_kt_yr=[x / 1e6 for x in river_kg], particles_kt_yr=[x / 1e6 for x in particle_kg],
                  floor_km2=shallows,
                  dissolved_g_per_m2_of_floor_under_20_m=[x * 1e3 / (shallows['under 20 m'] * 1e6) for x in river_kg],
                  all_river_phosphorus_g_per_m2_of_floor_under_20_m=[(a + b) * 1e3 / (shallows['under 20 m'] * 1e6)
                                                                     for a, b in zip(river_kg, particle_kg)])
    room = {}
    for horizon, s in safe.items():
        left = [s['supply_mg_m2_yr'] - supply['rivers_and_dust'][1], s['supply_mg_m2_yr'] - supply['rivers_and_dust'][0]]
        room[horizon] = dict(mg_m2_yr=left, kt_yr=[x * area / 1e12 for x in left])
    above = col.mid < VENTILATED_TO_M
    upper = float(col.volume[col.edge[1:] <= VENTILATED_TO_M].sum()
                  + sum((VENTILATED_TO_M - col.edge[i]) * 0.5 * (col.area[i] + np.interp(VENTILATED_TO_M, col.edge, col.area))
                        for i in range(col.n) if col.edge[i] < VENTILATED_TO_M < col.edge[i + 1]))
    still = (col.mid >= VENTILATED_TO_M) & (col.mid < 900.0)
    usable = float((col.volume[still] * (col.saturation - CC.HYPOXIC_MG_L)).sum()) / CL.G_O2_PER_MOL
    delivered = {f'{c:g} mg/m3': dict(lit_and_renewed_stock_mg_m2=c * upper,
                                      oxygen_it_takes_below_mol_m2=c * upper * O2_MMOL_PER_MG_P / 1e3,
                                      share_of_the_oxygen_above_2_mg_l_at_250_900_m=c * upper * O2_MMOL_PER_MG_P / 1e3 / usable)
                 for c in DELIVERED_MG_M3}
    return dict(salt=sea['salt'], area_km2=sea['area_km2'], mean_depth_m=sea['mean_depth_m'], saturation_mg_l=col.saturation,
                shares_of_the_exports_phosphorus=shares, supply_mg_m2_yr=supply, what_the_supply_does=does,
                supply_the_unaided_sea_can_take=safe, room_beside_rivers_and_dust=room, river_mouths=mouths,
                water_above_250_m_per_m2=upper, oxygen_above_2_mg_l_at_250_900_m_mol_m2=usable,
                delivered_water=delivered, levels_above_250_m=int(above.sum()))


def build(log=print):
    started = time.time()
    loads = json.loads(LOADS.read_text())
    phosphorus = json.loads(PHOSPHORUS.read_text())
    land = phosphorus['land_loss']
    seas = json.loads(SEAS.read_text())['seas']
    salt_col, fit, basis = CC.south_sea_column()
    fresh_saturation, _ = CC.fresh_saturation()
    stretch = CC.sinking_forms()['moon']
    reference_c = basis['export_g_o2'] / CL.G_O2_PER_MOL / CL.O2_PER_C * CL.G_C_PER_MOL
    dust = phosphorus['seas']['earth_ocean_air_deposition']['g_m2_yr']
    envelope_area = phosphorus['seas']['area_m2']
    inherited = [x / envelope_area * 1e3 for x in phosphorus['seas']['dust_kg_yr_inherited_envelope']]   # g/m2 a year
    out = {}
    for name, sea in seas.items():
        load = dict(loads['seas'][name], area_m2=sea['area_km2'] * 1e6)
        out[name] = one_sea(name, sea, load, sea_column(name, sea, salt_col, fresh_saturation), stretch, land, inherited)
        s = out[name]
        log(f'{name}: rivers {s["supply_mg_m2_yr"]["rivers_alone"][0]:.2f}-{s["supply_mg_m2_yr"]["rivers_alone"][1]:.2f} mg/m2/yr, '
            f'safe for 1,000 years {s["supply_the_unaided_sea_can_take"]["1000 years"]["supply_mg_m2_yr"]:.1f}, {time.time() - started:.0f} s')
    area = sum(s['area_km2'] for s in seas.values()) * 1e6
    total = lambda key: [sum(out[n]['supply_mg_m2_yr'][key][i] * seas[n]['area_km2'] * 1e6 for n in out) / 1e12 for i in (0, 1)]   # noqa: E731
    all_seas = dict(
        area_km2=area / 1e6,
        kt_yr=dict(rivers_dissolved=total('rivers_dissolved'), rivers_particles=total('rivers_particles'),
                   dust_inherited_envelope=total('dust'), dust_at_earths_oceans_rate=[dust[0] * area / 1e9, dust[1] * area / 1e9]),
        share_of_the_lands_loss_reaching_the_seas=dict(
            dissolved=sum(v['share_of_the_lands_dissolved_loss'] for v in loads['seas'].values()),
            particles={c: sum(v['share_of_the_lands_particles'][c] for v in loads['seas'].values()) for c in ('low', 'central', 'high')}),
        people=dict(mg_m2_yr_for_each_g_a_person_a_year_reaching_the_seas={f'{p / 1e9:g} billion': p * 1e3 / area for p in PEOPLE},
                    design_loss_g_a_person_a_year=[1e3 * x for x in phosphorus['people']['ecology_design_target']['lost_kg_per_person']],
                    through_the_lakes_share=sum(v['share_of_the_lands_dissolved_loss'] for v in loads['seas'].values())))
    return dict(
        schema=SCHEMA,
        producer=dict(module='research.studies.sea_nutrients.run',
                      code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),
        inputs={str(p.relative_to(ROOT)): sha256(p) for p in (LOADS, PHOSPHORUS, SEAS, CC.MIXING, CC.ESTIMATE, CC.SEA_OXYGEN)},
        evidence=('Arithmetic on the phosphorus study\'s loss from the land, this study\'s routing of it through the '
                  'design surface\'s lakes, and the canopy colony study\'s one-dimensional column of each sea. The '
                  'rule that a layered sea without an outlet spends 106 mol of oxygen below its renewed water for '
                  'each mol of phosphorus supplied assumes Redfield\'s matter and no return from below within the '
                  'horizon. The fresh seas borrow the salt sea\'s mixing. Not validated.'),
        reading_rule=('seas: for each sea, the phosphorus entering a year per m2 of sea (low and high: the rivers\' '
                      'dissolved part after the lakes, with what other seas pass on; the rivers\' particles, which '
                      'settle at the mouths, and the part of them salt water frees; and dust on the phosphorus '
                      'study\'s inherited envelope); the export each supply sustains with no colony and with the '
                      'canopy colony\'s central system; the oxygen it takes below 250 m; the supply the unaided sea '
                      'can take and keep 2 mg/L everywhere for 300, 1,000 and 3,000 years, and the room that leaves '
                      'beside rivers and dust; the river phosphorus against the floor shallower than 20 m; and what '
                      'the phosphorus in the delivered water does.'),
        setting=dict(p_mg_per_g_c=P_MG_PER_G_C, o2_mmol_per_mg_p=O2_MMOL_PER_MG_P, ventilated_to_m=VENTILATED_TO_M,
                     sinking_stretch=stretch, reference_export_g_c_m2_yr=reference_c,
                     supply_that_would_sustain_the_reference_export_mg_m2_yr={
                         n: reference_c * out[n]['shares_of_the_exports_phosphorus']['below_250_m'] * P_MG_PER_G_C for n in out},
                     land_loss_kg_yr=dict(dissolved=land['dissolved_kg_yr'], particulate=land['particulate_kg_yr']),
                     dust_g_m2_yr=dict(inherited_envelope=inherited, earths_oceans=dust),
                     reactive_particles=list(REACTIVE_PARTICLES), catch=CC.CATCH, trunk_foot_m=list(CC.TRUNKS[CC.CENTRAL_TRUNK])),
        all_seas=all_seas, seas=out, seconds=round(time.time() - started))


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    out = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print(f'wrote {OUT} in {out["seconds"]} s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
