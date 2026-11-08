"""Where the sea's dissolved load settles along the shores: the minerals the seas' water lays down as the plants and
the ground's evaporation take it, in the order its evaporation path gives them, and how fast they build a hardpan or a
sodic layer under the belts.

The plants draw the shore's groundwater from the capillary fringe and keep almost none of its salts, so each cubic
metre they transpire leaves its load in the root zone, where the water ends its path. The load is the sea's: the held
fresh sea of the sea chemistry study (calcium bicarbonate with silica, 0.14 g/kg), the left-alone sea as it passed
0.14 g/kg (with magnesium, sodium and potassium), and the held sea carrying the delivered water's own chloride
(resources' question 17). The evaporation path concentrates the water step by step at the root zone's CO2 and lays
down calcite at saturation (the carbon cycle's Plummer and Busenberg constants with Davies activities), amorphous
silica at its solubility (Gunnarsson and Arnorsson 2000, as the sea chemistry takes it), magnesium clays once the water
passes pH 9 (the sea chemistry's design threshold, one and a half silicons a magnesium), and keeps sodium, potassium and
chloride in solution: on the chemical divides of evaporating waters (Hardie and Eugster 1970) a water whose alkalinity
exceeds twice its calcium turns to soda after calcite, and the soda minerals form past the activities' range,
where only the masses are followed.
"""
from __future__ import annotations

import math

import numpy as np

from research.studies.carbon_cycle import run as cc
from research.studies.sea_chemistry import run as sc

SOIL_CO2_ATM = (4.86e-4, 0.01, 0.03)     # the air's 400 ppm at 1.216 bar, and a root zone's 1-3% (design range)
ROOT_ZONE_K = 298.15                     # the constants' 25 C, near the design climate's ground at 301 K
CEMENT_POROSITY = 0.3                    # the soil's pore space the cement fills (geography's shore screen)
ROOT_ZONE_M = 1.0                        # the depth over which the root zone gathers what the roots leave
BULK_DENSITY = 1500.0                    # kg/m3 of soil
CEC_MOL_KG = (0.05, 0.25)                # exchange capacity, 5-25 cmol(+)/kg, a loamy developed soil (recalled)
SODIC_ESP = 0.15                         # the exchangeable sodium share of a sodic soil (US Salinity Laboratory 1954)
SALINE_MOL_M3 = 40.0                     # soil water near 4 dS/m, the saline soil's threshold, as sodium chloride
CHLORIDE_G_KG = 0.28                     # the delivered ice's chloride (resources' question 17; sea chemistry)
DENSITY = dict(calcite=sc.MINERAL_DENSITY['calcite'], opal=sc.MINERAL_DENSITY['opal'],
               magnesium_clays=sc.MINERAL_DENSITY['magnesium_clays'], soda=sc.MINERAL_DENSITY['soda'],
               halite=2165.0)
MOLAR = dict(calcite=sc.MOLAR_MASS['calcite'], opal=sc.MOLAR_MASS['opal'],
             magnesium_clays=sc.MOLAR_MASS['magnesium_clays'], soda=0.5 * 105.99, halite=58.44)


def seas(chemistry):
    """The seas' waters (mol/kg, alkalinity eq/kg) from the sea chemistry product: the held fresh sea (managed sinks,
    central weathering, the billion-year row), the left-alone sea when it first passed 0.14 g/kg, and the held sea
    with the delivered chloride balanced by sodium."""
    def run(sinks):
        return next(r for r in chemistry['runs'] if r['corner'] == 'central' and r['delivery'] == 'fresh'
                    and r['sinks'] == sinks and r['sediment'] == 'kept_open')['rows']
    held = run('managed')[-1]
    natural = next(r for r in run('natural') if r['salinity_g_kg'] >= 0.14)

    def water(row, **extra):
        out = dict(ca=row['ca'], mg=row['mg'], na=row['na_k'], cl=row['cl'], so4=row['so4'], si=row['si'],
                   ta=row['alkalinity_meq_kg'], salinity_g_kg=row['salinity_g_kg'], year=row['year'])
        out.update(extra)
        return {k: (v / 1000.0 if k in ('ca', 'mg', 'na', 'cl', 'so4', 'si', 'ta') else v) for k, v in out.items()}
    cl = CHLORIDE_G_KG / sc.M_CL                                       # mol/kg
    held_w = water(held)
    return dict(held=held_w, left_alone=water(natural),
                held_with_chloride=dict(held_w, na=held_w['na'] + cl, cl=held_w['cl'] + cl,
                                        salinity_g_kg=held_w['salinity_g_kg'] + CHLORIDE_G_KG
                                        + cl * cc.M_NA))


def ionic(b):
    return 0.5 * (4.0 * b['ca'] + 4.0 * b['mg'] + b['na'] + b['cl'] + 4.0 * b['so4'] + max(b['ta'], 0.0)) / b['mass']


def evaporation_path(sea, pco2_atm, factors=None):
    """The water concentrated by the factors (each kilogram of sea water down to 1/factor kg), calcite and silica
    held at saturation and magnesium clays formed above pH 9 as it goes. Returns where each mineral starts (the
    concentration factor) and the share of the calcium, silica and magnesium each mineral holds at the factor where
    the Davies activities end (ionic strength 0.5); sodium, potassium and chloride stay dissolved to there."""
    factors = np.geomspace(1.0, 1e5, 401) if factors is None else factors
    c = cc.fresh_constants(ROOT_ZONE_K)
    b = dict(mass=1.0, ca=sea['ca'], mg=sea['mg'], na=sea['na'], cl=sea['cl'], so4=sea['so4'], si=sea['si'],
             ta=sea['ta'], dic=0.0)
    laid = dict(calcite=0.0, opal=0.0, magnesium_clays=0.0)
    onset, rows, end = {}, [], None
    for f in factors:
        b['mass'] = 1.0 / f
        if ionic(b) > sc.VALID_IONIC_STRENGTH:
            end = float(f)
            break
        box = dict(mass=b['mass'], ta=b['ta'], ca=b['ca'], mg=b['mg'], dic=0.0,
                   mono=b['na'] + b['cl'] + 4.0 * b['so4'])
        x = cc.precipitate(box, pco2_atm, c, 1.0) if b['ca'] > 0 else 0.0
        if x > 0:
            onset.setdefault('calcite', float(f))
            b['ca'], b['ta'] = box['ca'], box['ta']
            laid['calcite'] += x
        sp = cc.speciate(dict(box, ca=b['ca'], ta=b['ta']), pco2_atm, c)
        if b['mg'] > 0 and sp['ph'] > sc.PH_CLAY:
            onset.setdefault('magnesium_clays', float(f))
            m = min(b['mg'], b['si'] / sc.SI_PER_MG_CLAY, 0.5 * b['ta'])
            b['mg'] -= m
            b['ta'] -= 2.0 * m
            b['si'] -= sc.SI_PER_MG_CLAY * m
            laid['magnesium_clays'] += m
        limit = 10 ** sc.SILICA_LOG_K_25C * (1.0 + 10 ** (sp['ph'] - sc.SILICA_PKA)) * b['mass']
        if b['si'] > limit:
            onset.setdefault('opal', float(f))
            laid['opal'] += b['si'] - limit
            b['si'] = limit
        rows.append(dict(factor=float(f), ph=sp['ph'], omega=sp['omega'], ionic_strength=ionic(b),
                         ca_mmol_kg=b['ca'] / b['mass'] * 1e3, ta_meq_kg=b['ta'] / b['mass'] * 1e3,
                         si_mmol_kg=b['si'] / b['mass'] * 1e3))
    held = dict(calcite=laid['calcite'] / sea['ca'] if sea['ca'] else None,
                opal=laid['opal'] / sea['si'] if sea['si'] else None,
                magnesium_clays=laid['magnesium_clays'] / sea['mg'] if sea['mg'] else None)
    return dict(pco2_atm=pco2_atm, onset_factor=onset, share_laid_by_end=held, activities_end_factor=end,
                order=sorted(onset, key=onset.get), path=rows[::20])


def deposits(sea, draw_m_yr):
    """What the root zone gathers under a belt drawing draw_m_yr (m of water a year): every mineral of the sea's load
    in the end, the calcium as calcite, the silica as opal (less what the clays take), the magnesium as clays, the
    sodium and potassium with the remaining alkalinity as soda, the chloride as halite. Returns the solid's
    thickness and the cemented soil's (mm per century), and for the soluble salts the years to a sodic or a saline
    root zone."""
    w = draw_m_yr * 1000.0                                   # kg of water per m2 a year
    mol = dict(calcite=w * sea['ca'], magnesium_clays=w * sea['mg'],
               opal=w * max(sea['si'] - sc.SI_PER_MG_CLAY * sea['mg'], 0.0),
               soda=w * max(sea['na'] - sea['cl'], 0.0), halite=w * sea['cl'])
    solid = {k: v * MOLAR[k] / 1000.0 / DENSITY[k] * 1e5 for k, v in mol.items()}      # mm per century
    hardpan = solid['calcite'] + solid['opal'] + solid['magnesium_clays']
    sodium = w * sea['na']                                   # mol/m2 a year
    needed = [SODIC_ESP * cec * BULK_DENSITY * ROOT_ZONE_M for cec in CEC_MOL_KG]
    chloride = w * sea['cl']
    pore = CEMENT_POROSITY * ROOT_ZONE_M                     # m3 of soil water a m2 of root zone
    return dict(draw_m_yr=draw_m_yr, load_g_m2_yr=w * sea['salinity_g_kg'] / 1000.0 * 1000.0,
                mol_m2_yr=mol, solid_mm_per_century=solid,
                hardpan_solid_mm_per_century=hardpan, hardpan_cemented_mm_per_century=hardpan / CEMENT_POROSITY,
                sodium_mol_m2_per_century=sodium * 100.0,
                years_to_sodic_root_zone=[n / sodium if sodium > 0 else None for n in needed],
                years_to_saline_root_zone=(SALINE_MOL_M3 * pore / chloride) if chloride > 0 else None,
                halite_kg_m2_per_century=chloride * MOLAR['halite'] / 1000.0 * 100.0)


def leaching_requirement(sea, draw_m_yr, limit_mol_m3=SALINE_MOL_M3):
    """The water (m a year) that must pass below the roots to hold the root zone's chloride and soda at the saline
    threshold: the draw times the sea's soluble salt over the limit."""
    soluble = max(sea['na'], sea['cl']) * 1000.0             # mol/m3
    return draw_m_yr * soluble / limit_mol_m3
