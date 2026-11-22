"""Cement for the returned rock: what it has to do, what can supply it, and how long it takes.

    python -m research.studies.rock_return.cement      # checks itself and prints the tables, in under a second

The rock return study sets 146 g of rock down on each square metre of dry land a year and finds that the carriers'
own nitrogen cannot cement it. This module sizes the cement: the strength thinly spread rock and thick rim rock need
at lunar gravity, the rates at which Earth's organisms and soils lay binder and cement, the stock that weathering
and the seas' stone give, and the water that takes carbonate away. It holds the arithmetic as small pure functions.

It was worked out by a delegated agent on 22 November 2026 from the study's products and from outside sources read
mostly as search summaries (sources.json records how far each was read), and was read through and run before it was
adopted here. The project values it needs are repeated as constants, each with its source in a comment;
test_cement.py checks those copied from this study's own products against them. Binding the rest to their products
is still to do.

Every constant carries one of five labels:
  [project]   a committed project product or note, named, read on 22 November 2026
  [measured]  an Earth measurement, with its key in sources.json, where the depth of reading is recorded
  [computed]  arithmetic done here
  [recalled]  a value from memory, not looked up
  [guess]     a design guess
"""
from __future__ import annotations

import math

# ---- the job, from the study --------------------------------------------------------------------------------------
MOON_G = 1.624                    # m/s2 [project] the brief's value; shared/constants.json gives GM/R^2 = 1.62422
EARTH_G = 9.80665                 # m/s2 [project] shared/constants.json
ROCK_GT_YR = dict(slow=0.6467, central=3.2334, fast=6.4669)   # [project] rock_return.json job.gt_yr
RETURN_G_M2_YR = 146.107          # [project] rock_return.json job.kg_m2_yr_central
DRY_LAND_M2 = ROCK_GT_YR['central'] * 1e15 / RETURN_G_M2_YR   # [computed] 2.213e13 m2, the design surface's dry land
ROCK_KG_M3 = 2600.0               # [project] terrain study and carriers.py
LOOSE_SOIL_KG_M3 = 1300.0         # [project] rock_return.json setting.loose_soil_kg_m3
MEAN_LOWERING_M_MYR = 56.19       # [project] rock_return.json job.mean_lowering_m_per_myr
SEAS_STONE_MT_YR = dict(calcite=149.665, magnesite=94.8914, opal=281.581)   # [project] rock_return.json staying_put
CARRIERS_N_MT_YR = (0.194, 0.670)         # [project] albatross-scale carriers' nitrogen, rock_return.json staying_put
SMALL_CARRIERS_FOOD_MT_YR = (379.0, 381.0)  # [project] 60 g carriers' dry food, rock_return README
FOOD_NITROGEN = (0.01, 0.03)              # [project, a design guess there] rock_return.json setting.food_nitrogen
PLANT_GROWTH_G_M2_YR = (1440.0, 1900.0)   # [project] dry matter, rock_return README from the biomass comparison

# ---- slopes -------------------------------------------------------------------------------------------------------
FRICTION_DEG = 30.0               # [guess] the study's repose slope, taken as the cemented rock's friction angle too
DRY_BULK_KG_M3 = (1300.0, 1800.0)   # [guess] loose soil to a settled deposit
SATURATED_KG_M3 = (1800.0, 2100.0)  # [guess] the same deposits with their pores full of water
WATER_KG_M3 = 1000.0

# ---- the weathering that frees calcium, magnesium, silicon and iron in place -------------------------------------
DESSERT_A, DESSERT_B = 323.44, 0.0642   # [project] carbon_cycle/run.py: mol CO2 km-2 yr-1 per mm of runoff, exp(B T/C)
LAND_AIR_C = 294.3 - 273.15             # [project] carbon_cycle/run.py LAND_AIR_K
RUNOFF_MM_YR = 265.1                    # [project] carbon_cycle.json setting.runoff_mm_yr
CATION_SHARES = dict(ca=0.6501, mg=0.3284, mono=0.0214)   # [project] carbon_cycle.json setting.cation_shares, by eq
SILICA_MOL_PER_EQ = 0.87                # [project] sea_chemistry README, "0.87 mol of silica ... for each equivalent"
HIGHLAND_OXIDES = dict(SiO2=0.450, MgO=0.057, CaO=0.157, Na2O=0.0046, K2O=0.0017)   # [project] lunar_cycle_ecology
MARE_OXIDES = dict(SiO2=0.4305, MgO=0.0935, CaO=0.114, Na2O=0.00365, K2O=0.0016)    # run.py ROCK_WT_PCT (Sourcebook 7.15)
HIGHLAND_FEO = 0.051                    # [project] carbon_cycle/run.py HIGHLAND_FEO_WT (Apollo 16 soil)
MARE_FEO = 0.15                         # [project] the soils note, section 1.2, "about 15% FeO" (resources findings)
YOUNG_SOIL_FACTOR = (1.0, 10.0)         # [project, a design guess there] carbon_cycle FRESH_FACTOR, first centuries
M = dict(CaCO3=100.087, MgCO3=84.314, SiO2=60.084, CaO=56.077, MgO=40.304, Na2O=61.979, K2O=94.196, FeO=71.844,
         FeOOH=88.852, Fe=55.845, C=12.011, N=14.007)

# ---- water, by the design years' aridity class -------------------------------------------------------------------
# [project] the soils note on the research drive (session-tools/agents/soils/soils_maintenance.md, section 1.1):
# share of the land, rain in mm a year. Water leaving below the roots: the storm-water study's recharge, 1.3-17% of
# the rain on the drylands; in the humid land the rain less the demand (695 mm) down to the land's whole runoff
# gathered on the humid land (265.1 / 0.422 = 628 mm).
CLASSES = (('hyper-arid', 0.042, 19.0), ('arid', 0.191, 70.0), ('semi-arid', 0.251, 243.0),
           ('dry sub-humid', 0.094, 433.0), ('humid', 0.422, 1717.0))
DRYLAND_RECHARGE_SHARE = (0.013, 0.17)   # [project] storm-water study through the soils note
HUMID_WATER_LEAVING_MM = (628.0, 695.0)  # [computed] from the project figures above
AIR_PCO2_ATM = 48.6 / 101325.0           # [project] carbon_cycle README: 400 ppm, 48.6 Pa
SOIL_PCO2_ATM = 0.026                    # [project] springs README table: soil CO2 2.6% within 15 degrees

# ---- what Earth's organisms and soils are measured to do ---------------------------------------------------------
# Rates in g of cement per m2 a year unless stated; keys are in sources.json with the depth of reading.
TUFA_KG_M3 = (1000.0, 2000.0)            # [recalled] porous tufa and stromatolite carbonate, used to convert mm to g
MM_YR = dict(
    plitvice_barrier=(13.5, 13.5),       # [project] terrain study README
    stream_tufa=(1.0, 9.0),              # [measured] pedley_2009, river_piedra, pentecost: 1, 2, 3-9, 5-9 mm a year
    stromatolite=(0.1, 0.5),             # [measured] chivas_1990 (<= 0.4), collins_jahnert_2014 (< 0.1 to 0.5)
)
STONEWORT_G_M2 = (440.0, 690.0)          # [project] terrain study README (Pelechaty et al. 2013), a year's standing lime
PEDOGENIC_G_M2_YR = (1.0, 5.0)           # [measured] marion_1985: most field studies; schlesinger gives 1.0-3.5
PEDOGENIC_EXTREMES = (0.16, 20.0)        # [measured] bockheim_douglass_2006 low; calciaquoll at a wetland high
EARTHWORM_MOL_HA_YR = (18.0, 3139.0)     # [measured] lambkin_2011: calcite granules, 10-20 Lumbricus terrestris a m2
IROKO_KG_C_TREE_YR = (5.76, 12.5)        # [measured] 5.76: a thesis figure seen in a search summary; 12.5 = about 1 t of
                                         # mineral carbon round an 80-year-old stump (cailleau_2004 abstract) / 80 years
IROKO_GROUND_M2 = (200.0, 700.0)         # [guess] the ground one tree works, a crown 16-30 m across
BIOCRUST_G_C_M2_YR = (2.75, 16.0)        # [measured] elbert_2012 median 16 for dryland crusts; field values 2.75, 11.7
FERRICRETE_MM_KYR = (0.4, 17.0)          # [measured] goudie_1985 through fenske_2025: Uganda to the Liberian Shield
CALCRETE_MM_KYR = (0.5, 350.0)           # [measured] goudie_1985 through fenske_2025: pedogenic calcretes
EROSION_M_MYR = dict(canga=(0.08, 0.54), silcrete=(0.2, 0.2), arid_bedrock=(1.0, 5.1))   # [measured] shuster_2012,
                                         # monteiro_2018; central Australia; Namib and Australian outcrops
CANGA_TRIAL_MONTHS = 16.0                # [measured] gagen_2020: ten monthly feeds and six months of drying
CANGA_TRIAL_GLUCOSE_G_M3_MONTH = (108.0, 254.0)   # [computed from measured] 20-47 L of 30 mM glucose a month per m3
CRUSHED_CANGA_KG_M3 = (2000.0, 3000.0)   # [guess] bulk density of the trial's crushed canga
ROOT_COHESION_KPA = dict(clearcut=(0.0, 10.0), industrial_forest=(6.8, 23.2), natural_forest=(25.6, 94.3))
                                         # [measured] schmidt_2001 abstract; a 2022 study finds under 30% of these
WEAK_CEMENT_UCS_KPA = dict(weak=(0.0, 100.0), moderate=(100.0, 400.0))   # [measured] collins_sitar_2009 classes
MICP_ANCHORS = dict(at_contacts=(0.057, 1140.0), over_grains=(0.143, 1000.0))   # [measured, tentative] cheng_2013:
                                         # calcite per sand by mass and UCS in kPa, at low and at full saturation
MICP_THRESHOLD_KG_M3 = 60.0              # [measured, secondhand] whiffin_2007 through van Paassen: least calcite a m3
BARE_DRYLAND_LOSS_MM_YR = (0.13, 0.63)   # [project] storm-water README: bare soils, an upper bound
CRUSTED_LOSS_MM_YR = 0.013               # [project] storm-water README: the study's crust
NATIVE_COVER_LOSS_MM_YR = (0.013, 0.053)  # [measured, through the project] Montgomery 2007 median and mean
CRUST_CUT = (0.68, 0.85)                 # [measured, through the project] Eldridge et al. 2020: sediment cut by crusts
COVER_EXPONENT = (0.0435, 0.053)         # [measured, through the project] Elwell & Stocking 1976, per % of cover
# wind, from the project's sand study (research/studies/atlas_round/sand.md), which read Shao & Lu (2000) in Kok 2012
SHAO_LU_A, GRAIN_COHESION_N_M = 0.111, 3.0e-4
SURFACE_AIR_KG_M3, FELDSPAR_KG_M3 = 1.35, 2730.0
SAND_ROUGHNESS_M, VON_KARMAN = 1.0e-3, 0.4
GUSTS_M_S = (5.0, 9.0)                   # [project] the soils note from the desert analysis: storm gusts at 10 m
CRUST_WIND_EXPONENT = 0.0157             # [measured, through the project] Rodriguez-Caballero et al. 2022, per % cover
ORGANIC_BINDER_SHARE = (0.01, 0.05)      # [guess] organic matter of a bound mineral soil, by mass
SECRETION_NITROGEN = (0.03, 0.16)        # [guess] a mucus glycoprotein to a plain protein


# ---- 1. what the cement has to do ---------------------------------------------------------------------------------
def cohesion_needed_pa(thickness_m, slope_deg, bulk_kg_m3, gravity=MOON_G, friction_deg=FRICTION_DEG,
                       saturated=False, safety=1.0):
    """Cohesion on a plane parallel to the surface at a depth, for an infinite slope to stand with a safety factor.
    Dry: c = z g rho cos(b) (F sin(b) - cos(b) tan(phi)). Saturated with seepage down the slope and the water table at
    the surface: the weight is the saturated one and the friction acts on the buoyant weight."""
    b, phi = math.radians(slope_deg), math.radians(friction_deg)
    drive = safety * bulk_kg_m3 * math.sin(b)
    resist = ((bulk_kg_m3 - WATER_KG_M3) if saturated else bulk_kg_m3) * math.cos(b) * math.tan(phi)
    return max(thickness_m * gravity * math.cos(b) * (drive - resist), 0.0)


def stable_thickness_m(cohesion_pa, slope_deg, bulk_kg_m3, **kw):
    """The thickness a cohesion holds: the inverse of cohesion_needed_pa, which is linear in the thickness."""
    per_metre = cohesion_needed_pa(1.0, slope_deg, bulk_kg_m3, **kw)
    return math.inf if per_metre == 0 else cohesion_pa / per_metre


def culmann_height_m(cohesion_pa, slope_deg, bulk_kg_m3, gravity=MOON_G, friction_deg=FRICTION_DEG):
    """Culmann's critical height of a planar slope face: 4 c sin(b) cos(phi) / (g rho (1 - cos(b - phi)))."""
    b, phi = math.radians(slope_deg), math.radians(friction_deg)
    return 4.0 * cohesion_pa * math.sin(b) * math.cos(phi) / (gravity * bulk_kg_m3 * (1.0 - math.cos(b - phi)))


def ucs_from_cohesion(cohesion, friction_deg=FRICTION_DEG):
    """Unconfined compressive strength of a Mohr-Coulomb material: 2 c tan(45 + phi/2)."""
    return 2.0 * cohesion * math.tan(math.radians(45.0 + friction_deg / 2.0))


def cohesion_from_ucs(ucs, friction_deg=FRICTION_DEG):
    return ucs / (2.0 * math.tan(math.radians(45.0 + friction_deg / 2.0)))


def cement_for_cohesion(cohesion_kpa, settled_kg_m3=1600.0):
    """The mass fraction of calcite cement that gives a cohesion, between two readings of the measured anchors
    [guess on measured anchors]. Best: all of it at the grain contacts, a line through the origin and the anchor
    treated at low saturation. Worst: spread over the grains, nothing below the threshold of 60 kg a m3 and a line
    from there to the anchor treated at full saturation."""
    f1, ucs1 = MICP_ANCHORS['at_contacts']
    f2, ucs2 = MICP_ANCHORS['over_grains']
    best = cohesion_kpa * f1 / cohesion_from_ucs(ucs1)
    f0 = MICP_THRESHOLD_KG_M3 / settled_kg_m3
    worst = f0 + cohesion_kpa * (f2 - f0) / cohesion_from_ucs(ucs2)
    return best, worst


def threshold_friction_velocity(grain_m, grain_kg_m3=FELDSPAR_KG_M3, air_kg_m3=SURFACE_AIR_KG_M3, gravity=MOON_G):
    """Shao and Lu's friction velocity at which wind starts loose dry grains, as the project's sand study carries it."""
    return SHAO_LU_A * math.sqrt((grain_kg_m3 - air_kg_m3) / air_kg_m3 * gravity * grain_m
                                 + GRAIN_COHESION_N_M / (air_kg_m3 * grain_m))


def easiest_grain_m(grain_kg_m3=FELDSPAR_KG_M3, air_kg_m3=SURFACE_AIR_KG_M3, gravity=MOON_G):
    return math.sqrt(GRAIN_COHESION_N_M / ((grain_kg_m3 - air_kg_m3) * gravity))


def wind_at_10_m(friction_velocity, roughness_m=SAND_ROUGHNESS_M):
    return friction_velocity / VON_KARMAN * math.log(10.0 / roughness_m)


def crust_cover_for(multiplier):
    """The crust cover, in percent, that raises the wind threshold by a factor, on the mean fit e^(0.0157 cover)."""
    return math.log(multiplier) / CRUST_WIND_EXPONENT


# ---- 3. supply ----------------------------------------------------------------------------------------------------
def cement_need_g_m2_yr(fraction, concentration=1.0, returned=RETURN_G_M2_YR):
    """Cement a year for a mass fraction of the rock set down, where the rock of `concentration` m2 lands on one."""
    return fraction * returned * concentration


def carried_fractions(rock_gt_yr):
    """The seas' stone as mass fractions of the particles, if all of it is laid within reach and carried up."""
    return {k: v / (rock_gt_yr * 1000.0) for k, v in SEAS_STONE_MT_YR.items()}


def weathering_eq_m2_yr(runoff_mm=RUNOFF_MM_YR, temp_c=LAND_AIR_C):
    """Dessert's basalt law as the carbon cycle carries it: equivalents of cation charge freed per m2 a year."""
    return DESSERT_A * math.exp(DESSERT_B * temp_c) * runoff_mm / 1e6


def rock_g_per_eq(oxides):
    eq = 2.0 * oxides['CaO'] / M['CaO'] + 2.0 * oxides['MgO'] / M['MgO'] \
        + 2.0 * oxides['Na2O'] / M['Na2O'] + 2.0 * oxides['K2O'] / M['K2O']
    return 1.0 / eq


def in_place_supply(eq_m2_yr, feo=HIGHLAND_FEO, oxides=HIGHLAND_OXIDES):
    """What weathering frees on the spot, as the cement each element could make, g per m2 a year."""
    rock = eq_m2_yr * rock_g_per_eq(oxides)
    return dict(calcite=eq_m2_yr * CATION_SHARES['ca'] / 2.0 * M['CaCO3'],
                magnesite=eq_m2_yr * CATION_SHARES['mg'] / 2.0 * M['MgCO3'],
                opal=eq_m2_yr * SILICA_MOL_PER_EQ * M['SiO2'],
                goethite=rock * feo * M['FeOOH'] / M['FeO'], rock_dissolved=rock)


def goethite_equivalent(feo):
    """The iron of a rock as goethite, by mass of the rock: the iron cement the particles carry in themselves."""
    return feo * M['FeOOH'] / M['FeO']


def carbonate_constants(temp_k):
    """Plummer and Busenberg's (1982) fits for Henry's constant, carbonic acid and calcite [recalled; the 25 C values
    they give are checked in self_check]."""
    t, lg = temp_k, math.log10(temp_k)
    kh = 108.3865 + 0.01985076 * t - 6919.53 / t - 40.45154 * lg + 669365.0 / t ** 2
    k1 = -356.3094 - 0.06091964 * t + 21834.37 / t + 126.8339 * lg - 1684915.0 / t ** 2
    k2 = -107.8871 - 0.03252849 * t + 5151.79 / t + 38.92561 * lg - 563713.9 / t ** 2
    ks = -171.9065 - 0.077993 * t + 2839.319 / t + 71.595 * lg
    return tuple(10.0 ** x for x in (kh, k1, k2, ks))


def calcite_solubility_g_m3(pco2_atm, temp_k=LAND_AIR_C + 273.15):
    """Calcite that pure water holds at a CO2 pressure, open to the gas, with Davies activity coefficients."""
    kh, k1, k2, ks = carbonate_constants(temp_k)
    product = k1 * ks * kh * pco2_atm / k2            # a(Ca) a(HCO3)^2
    ca = (product / 4.0) ** (1.0 / 3.0)
    for _ in range(50):
        ionic = 3.0 * ca
        root = math.sqrt(ionic)
        g1 = 10.0 ** (-0.507 * (root / (1.0 + root) - 0.3 * ionic))
        ca = (product / (4.0 * g1 ** 4 * g1 ** 2)) ** (1.0 / 3.0)   # gamma2 = gamma1^4
    return ca * M['CaCO3'] * 1000.0


def carbonate_moved_g_m2_yr(water_mm, solubility_g_m3):
    """Calcite that a depth of water can dissolve and carry: off the ground if it leaves, to the contacts if it dries."""
    return water_mm * 1e-3 * solubility_g_m3


def cement_lifetime_yr(fraction, depth_m, bulk_kg_m3, water_mm, solubility_g_m3):
    """Years for water leaving at saturation to take a carbonate cement out of a layer."""
    return fraction * depth_m * bulk_kg_m3 * 1000.0 / carbonate_moved_g_m2_yr(water_mm, solubility_g_m3)


def carbon_for(cement_g, route):
    """Carbon an organism spends for a mass of cement, g: light-driven calcification fixes one carbon a calcite; the
    oxalate route burns one oxalate, two carbons, a calcite; iron takes one organic carbon for four irons reduced
    (cement_g is then goethite); the urea route takes two nitrogens a calcite (returned as g of nitrogen)."""
    if route == 'light':
        return cement_g / M['CaCO3'] * M['C']
    if route == 'oxalate':
        return cement_g / M['CaCO3'] * 2.0 * M['C']
    if route == 'iron':
        return cement_g / M['FeOOH'] / 4.0 * M['C']
    if route == 'urea_nitrogen':
        return cement_g / M['CaCO3'] * 2.0 * M['N']
    raise ValueError(route)


def secretion_share_limit(nitrogen_mt_yr, rock_gt_yr, nitrogen_content):
    """The largest binder, as a share of the load, that carriers could secrete if all their nitrogen went into it."""
    return nitrogen_mt_yr / (rock_gt_yr * 1000.0 * nitrogen_content)


# ---- 4. time ------------------------------------------------------------------------------------------------------
def deposit_mm_yr(concentration=1.0, bulk_kg_m3=LOOSE_SOIL_KG_M3, returned=RETURN_G_M2_YR):
    return returned * concentration / bulk_kg_m3


def fraction_kept_up(rate_g_m2_yr, concentration=1.0, carried=0.0, returned=RETURN_G_M2_YR):
    """The cement fraction a deposit ends with when a process lays cement at a rate per m2 while it builds."""
    return carried + rate_g_m2_yr / (returned * concentration)


def concentration_served(rate_g_m2_yr, fraction, carried=0.0, returned=RETURN_G_M2_YR):
    """How many m2 of land's rock one m2 of deposit can take while a process still cements it to a fraction."""
    need = max(fraction - carried, 0.0)
    return math.inf if need == 0 else rate_g_m2_yr / (returned * need)


def years_per_century_of_return(rate_g_m2_yr, fraction, returned=RETURN_G_M2_YR):
    """Years a process takes to cement a century of returned rock to a fraction; under 100 it keeps pace."""
    return 100.0 * returned * fraction / rate_g_m2_yr


def mm_to_g_m2_yr(mm_yr, bulk_kg_m3):
    return mm_yr * bulk_kg_m3


def process_rates():
    """Measured rates as g of carbonate per m2 a year, low and high, with the condition each needs."""
    tufa = lambda key: (mm_to_g_m2_yr(MM_YR[key][0], TUFA_KG_M3[0]), mm_to_g_m2_yr(MM_YR[key][1], TUFA_KG_M3[1]))
    iroko = (IROKO_KG_C_TREE_YR[0] * 1000.0 * M['CaCO3'] / M['C'] / IROKO_GROUND_M2[1],
             IROKO_KG_C_TREE_YR[1] * 1000.0 * M['CaCO3'] / M['C'] / IROKO_GROUND_M2[0])
    worms = tuple(x * M['CaCO3'] / 1e4 for x in EARTHWORM_MOL_HA_YR)
    crust = tuple(x * M['CaCO3'] / M['C'] for x in BIOCRUST_G_C_M2_YR)
    return [
        ('Soil carbonate, drylands (calcrete)', PEDOGENIC_G_M2_YR, 'dry land; calcium from dust and rain'),
        ('Earthworm calcite granules', worms, 'moist soil with calcium; granules, loose'),
        ('Dryland crust at one calcite a carbon (potential)', crust, 'wet hours in light; calcium at hand; not measured'),
        ('Oxalate tree (iroko), under the tree', iroko, 'wet tropics; calcium from rock, dust or rain'),
        ('Stromatolites', tufa('stromatolite'), 'standing water above saturation'),
        ('Stonewort meadow', STONEWORT_G_M2, 'hard-water lake, 0.5-6 m deep'),
        ('Stream tufa', tufa('stream_tufa'), 'flowing water above 2-5.6 times saturation'),
        ('Plitvice barriers', tufa('plitvice_barrier'), 'the same, at a waterfall'),
    ]


# ---- checks -------------------------------------------------------------------------------------------------------
def self_check():
    kh, k1, k2, ks = carbonate_constants(298.15)
    for got, want in ((kh, -1.468), (k1, -6.352), (k2, -10.329), (ks, -8.480)):
        assert abs(math.log10(got) - want) < 0.003, (math.log10(got), want)
    assert abs(threshold_friction_velocity(easiest_grain_m()) - 0.145) < 0.001      # the sand study's 0.145 m/s
    assert abs(easiest_grain_m() - 260e-6) < 2e-6                                    # and its 260 micrometres
    assert abs(weathering_eq_m2_yr() * 27.31e12 / 9.1e12 - 1.0) < 0.01               # the carbon cycle's 9.1e12 eq
    assert abs(cohesion_needed_pa(1.0, 45.0, 1300.0) - 446.1) < 0.5
    assert abs(cohesion_needed_pa(1.0, 45.0, 1300.0, gravity=EARTH_G)
               / cohesion_needed_pa(1.0, 45.0, 1300.0) - EARTH_G / MOON_G) < 1e-9
    assert cohesion_needed_pa(5.0, 30.0, 1500.0) < 1e-9                              # at repose nothing is needed
    assert cohesion_needed_pa(5.0, 20.0, 1500.0) == 0.0
    assert abs(stable_thickness_m(cohesion_needed_pa(7.0, 60.0, 1800.0), 60.0, 1800.0) - 7.0) < 1e-9
    assert abs(cohesion_needed_pa(1.0, 45.0, 2100.0, saturated=True) - 1190.0) < 1.0
    assert abs(calcite_solubility_g_m3(AIR_PCO2_ATM) - 59.0) < 1.0
    assert abs(calcite_solubility_g_m3(SOIL_PCO2_ATM) - 243.0) < 3.0
    here = in_place_supply(weathering_eq_m2_yr())
    assert abs(here['calcite'] - 10.8) < 0.1 and abs(here['magnesite'] - 4.6) < 0.1 and abs(here['opal'] - 17.4) < 0.1
    assert abs(here['calcite'] / M['CaCO3'] / (here['magnesite'] / M['MgCO3'])
               - CATION_SHARES['ca'] / CATION_SHARES['mg']) < 1e-9
    stone = carried_fractions(ROCK_GT_YR['central'])
    assert abs(stone['calcite'] - 0.0463) < 0.0002 and abs(stone['opal'] - 0.0871) < 0.0002
    assert abs(stone['calcite'] + stone['magnesite'] - 0.0756) < 0.0001              # rock_return.json's lime, central
    assert abs(years_per_century_of_return(0.03 * RETURN_G_M2_YR, 0.03) - 100.0) < 1e-9
    assert abs(fraction_kept_up(20.0, concentration_served(20.0, 0.03)) - 0.03) < 1e-12
    assert concentration_served(20.0, 0.03, carried=0.03) == math.inf
    assert abs(carbon_for(100.087, 'urea_nitrogen') - 28.014) < 1e-9
    assert abs(carbon_for(100.087, 'oxalate') - 2.0 * carbon_for(100.087, 'light')) < 1e-12


# ---- the tables ---------------------------------------------------------------------------------------------------
def num(x):
    """A figure to about three significant digits, without exponents."""
    a = abs(x)
    if a == 0:
        return '0'
    if a >= 100:
        return f'{x:,.0f}'
    if a >= 10:
        return f'{x:.1f}'
    if a >= 1:
        return f'{x:.2f}'.rstrip('0').rstrip('.') if a < 1.995 else f'{x:.1f}'
    if a >= 0.1:
        return f'{x:.2f}'
    return f'{x:.2g}'


def rng(pair, fmt=None):
    lo, hi = pair
    show = num if fmt is None else fmt.format
    return show(lo) if abs(hi - lo) < 1e-12 * max(abs(hi), 1.0) else f'{show(lo)}-{show(hi)}'


def main():
    self_check()
    out = print
    need = {f: cement_need_g_m2_yr(f) for f in (0.01, 0.03, 0.10)}

    out('1a. ORDINARY GROUND: what a bound surface has to hold')
    layer_mm = deposit_mm_yr()
    out(f'  a year of returned rock: {RETURN_G_M2_YR:.0f} g/m2, {layer_mm:.3f} mm as loose soil '
        f'({100 * layer_mm:.1f} mm a century; {100 * RETURN_G_M2_YR / ROCK_KG_M3:.2f} mm as rock)')
    lost = [layer_mm / x * 12.0 for x in BARE_DRYLAND_LOSS_MM_YR[::-1]]
    out(f'  left bare on dry land it is stripped at {rng(BARE_DRYLAND_LOSS_MM_YR)} mm/yr (upper bound): a year\'s layer '
        f'lasts {rng(lost, "{:.1f}")} months')
    d = easiest_grain_m()
    u = threshold_friction_velocity(d)
    w = wind_at_10_m(u)
    out(f'  loose grains start at u* = {u:.3f} m/s ({1e6 * d:.0f} um), a 10 m wind of {w:.1f} m/s; '
        f'Earth quartz {threshold_friction_velocity(easiest_grain_m(2650.0, 1.225, EARTH_G), 2650.0, 1.225, EARTH_G):.3f} m/s')
    mult = [g / w for g in GUSTS_M_S]
    out(f'  storm gusts of {rng(GUSTS_M_S)} m/s ask the threshold raised {rng(mult, "{:.1f}")} times: crust cover of '
        f'{rng([crust_cover_for(m) for m in mult], "{:.0f}")}%')
    cuts = (BARE_DRYLAND_LOSS_MM_YR[0] / NATIVE_COVER_LOSS_MM_YR[1], BARE_DRYLAND_LOSS_MM_YR[1] / NATIVE_COVER_LOSS_MM_YR[0])
    out(f'  to wear as Earth\'s ground under native cover ({rng(NATIVE_COVER_LOSS_MM_YR)} mm/yr) the bare loss has to '
        f'fall {rng(cuts, "{:.1f}")} times; Earth\'s crusts cut {rng([1 / (1 - c) for c in CRUST_CUT], "{:.1f}")} times, '
        f'60% plant cover {rng([math.exp(b * 60) for b in COVER_EXPONENT], "{:.0f}")} times')
    binder = [s * RETURN_G_M2_YR for s in ORGANIC_BINDER_SHARE]
    out(f'  organic binder at {rng([100 * s for s in ORGANIC_BINDER_SHARE])}% of the layer: {rng(binder, "{:.1f}")} g/m2/yr, '
        f'{rng([100 * binder[0] / PLANT_GROWTH_G_M2_YR[1], 100 * binder[1] / PLANT_GROWTH_G_M2_YR[0]], "{:.2f}")}% of the '
        f'plants\' growth')

    out('\n1b. RIM ROCK: cohesion for a layer to stand, kPa (friction 30 degrees)')
    out('  thickness  slope  dry, Moon      saturated with seepage  with safety 1.5, saturated  dry, Earth')
    for z in (1.0, 10.0, 100.0):
        for b in (35.0, 45.0, 60.0):
            dry = [cohesion_needed_pa(z, b, r) / 1e3 for r in DRY_BULK_KG_M3]
            wet = [cohesion_needed_pa(z, b, r, saturated=True) / 1e3 for r in SATURATED_KG_M3]
            safe = [cohesion_needed_pa(z, b, r, saturated=True, safety=1.5) / 1e3 for r in SATURATED_KG_M3]
            earth = [cohesion_needed_pa(z, b, r, gravity=EARTH_G) / 1e3 for r in DRY_BULK_KG_M3]
            out(f'  {z:6.0f} m   {b:.0f}    {rng(dry, "{:.2f}"):13s}  {rng(wet, "{:.2f}"):22s}  '
                f'{rng(safe, "{:.2f}"):26s}  {rng(earth, "{:.1f}")}')
    out('  the thickness a cohesion holds at 45 degrees (dry, 1,800 kg/m3; saturated, 2,100 kg/m3), Moon, and dry on Earth:')
    for c in (1.0, 10.0, 30.0, 100.0):
        out(f'    {c:5.0f} kPa: {stable_thickness_m(c * 1e3, 45.0, 1800.0):6.1f} m dry, '
            f'{stable_thickness_m(c * 1e3, 45.0, 2100.0, saturated=True):6.1f} m saturated; Earth '
            f'{stable_thickness_m(c * 1e3, 45.0, 1800.0, gravity=EARTH_G):5.1f} m; a vertical face '
            f'{culmann_height_m(c * 1e3, 90.0, 1800.0):5.1f} m (Earth {culmann_height_m(c * 1e3, 90.0, 1800.0, EARTH_G):.1f} m); '
            f'UCS {ucs_from_cohesion(c):.0f} kPa')
    out('  calcite cement for a cohesion, % of mass, at the contacts to spread over the grains [guess on anchors]:')
    for c in (1.0, 12.0, 30.0, 120.0, 290.0):
        best, worst = cement_for_cohesion(c)
        out(f'    {c:5.0f} kPa: {100 * best:.2f}-{100 * worst:.1f}%')

    out('\n3. SUPPLY, g per m2 of dry land a year (central erosion)')
    out(f'  need for a cement fraction of 1, 3, 10% of {RETURN_G_M2_YR:.0f} g: '
        + ', '.join(f'{v:.1f}' for v in need.values()))
    for case in ('slow', 'central', 'fast'):
        fr = carried_fractions(ROCK_GT_YR[case])
        out(f'  carried with the mud, {case} erosion: calcite {100 * fr["calcite"]:.1f}%, magnesite '
            f'{100 * fr["magnesite"]:.1f}%, opal {100 * fr["opal"]:.1f}% of the load'
            + (f' = {fr["calcite"] * RETURN_G_M2_YR:.1f}, {fr["magnesite"] * RETURN_G_M2_YR:.1f}, '
               f'{fr["opal"] * RETURN_G_M2_YR:.1f} g/m2/yr' if case == 'central' else ''))
    eq = weathering_eq_m2_yr()
    here = in_place_supply(eq)
    mare = in_place_supply(eq, MARE_FEO, MARE_OXIDES)
    out(f'  weathering in place at the catchment rate, {eq:.3f} eq/m2/yr, dissolving {here["rock_dissolved"]:.0f} g of '
        f'highland rock: calcite {here["calcite"]:.1f}, magnesite {here["magnesite"]:.1f}, opal {here["opal"]:.1f}, '
        f'goethite {here["goethite"]:.1f} (mare {mare["goethite"]:.1f}); young soils up to {YOUNG_SOIL_FACTOR[1]:.0f} times')
    out(f'  as shares of the returned rock: calcite {100 * here["calcite"] / RETURN_G_M2_YR:.1f}%, with magnesite '
        f'{100 * (here["calcite"] + here["magnesite"]) / RETURN_G_M2_YR:.1f}%, opal {100 * here["opal"] / RETURN_G_M2_YR:.1f}%')
    out(f'  iron the particles carry in themselves, as goethite: highland {100 * goethite_equivalent(HIGHLAND_FEO):.1f}% '
        f'({goethite_equivalent(HIGHLAND_FEO) * RETURN_G_M2_YR:.1f} g/m2/yr), mare {100 * goethite_equivalent(MARE_FEO):.1f}%')
    sol = (calcite_solubility_g_m3(AIR_PCO2_ATM), calcite_solubility_g_m3(SOIL_PCO2_ATM))
    out(f'  calcite a m3 of water holds: {sol[0]:.0f} g at the air\'s CO2, {sol[1]:.0f} g at the soil\'s 2.6%')
    supply = here['calcite'] + carried_fractions(ROCK_GT_YR['central'])['calcite'] * RETURN_G_M2_YR
    out(f'  one stock of calcium: weathering frees {here["calcite"] * DRY_LAND_M2 / 1e12:.0f} Mt/yr of calcite-equivalent on the '
        f'dry land, and {SEAS_STONE_MT_YR["calcite"]:.0f} Mt/yr of it reaches the seas; what is caught in place never '
        f'reaches the lakes. The most that can pass through the ground in a year, freed in place and carried back, is '
        f'{supply:.1f} g/m2/yr')
    out(f'  in-place calcium keeps 3% cement up to a concentration of {concentration_served(here["calcite"], 0.03):.1f}')
    out(f'  water: leaching takes that whole flux at {rng([supply / s * 1e3 for s in sol[::-1]], "{:.0f}")} mm/yr leaving '
        f'the soil')
    out('  class           land   rain   water leaving  calcite leached  moved to contacts by drying  land a m2 serves at 3%  3% lasts in 0.1 m')
    for name, share, rain in CLASSES:
        leaving = HUMID_WATER_LEAVING_MM if name == 'humid' else tuple(rain * s for s in DRYLAND_RECHARGE_SHARE)
        leached = (carbonate_moved_g_m2_yr(leaving[0], sol[0]), carbonate_moved_g_m2_yr(leaving[1], sol[1]))
        drying = (carbonate_moved_g_m2_yr(rain - leaving[1], sol[0]), carbonate_moved_g_m2_yr(rain - leaving[0], sol[1]))
        life = (cement_lifetime_yr(0.03, 0.1, LOOSE_SOIL_KG_M3, leaving[1], sol[1]),
                cement_lifetime_yr(0.03, 0.1, LOOSE_SOIL_KG_M3, leaving[0], sol[0]))
        serves = [x / need[0.03] for x in drying]
        out(f'  {name:14s} {100 * share:4.1f}%  {rain:5.0f}  {rng(leaving, "{:.1f}"):13s}  {rng(leached):15s}  '
            f'{rng(drying):27s}  {rng(serves):22s}  {rng(life)} yr')
    crust_kg = [0.010 * d for d in TUFA_KG_M3]
    crust_life = (crust_kg[0] * 1e3 / carbonate_moved_g_m2_yr(HUMID_WATER_LEAVING_MM[1], sol[1]),
                  crust_kg[1] * 1e3 / carbonate_moved_g_m2_yr(RUNOFF_MM_YR, sol[0]))
    out(f'  a 10 mm carbonate crust ({rng(crust_kg, "{:.0f}")} kg/m2) dissolves in {rng(crust_life)} years, from the humid '
        f'land\'s water at the soil\'s CO2 to the mean runoff at the air\'s')
    out('  cost to the organisms for 1, 3, 10% cement at the mean return, g/m2/yr:')
    for route, label in (('light', 'carbon fixed, one a calcite'), ('oxalate', 'carbon as oxalate'),
                         ('urea_nitrogen', 'nitrogen by the urea route, for comparison')):
        out(f'    {label}: ' + ', '.join(f'{carbon_for(v, route):.2f}' for v in need.values()))
    iron_g = goethite_equivalent(HIGHLAND_FEO) * RETURN_G_M2_YR
    out(f'    organic carbon to reduce all the particles\' iron once ({iron_g:.1f} g goethite): {carbon_for(iron_g, "iron"):.2f}')
    trial = [10.0 * g * 72.066 / 180.156 / k for g, k in zip(CANGA_TRIAL_GLUCOSE_G_M3_MONTH, CRUSHED_CANGA_KG_M3[::-1])]
    out(f'    the canga trial\'s dose: {rng(trial, "{:.2f}")} g C per kg cemented in {CANGA_TRIAL_MONTHS:.0f} months, '
        f'{rng([t * RETURN_G_M2_YR / 1e3 for t in trial], "{:.2f}")} g C/m2/yr for the mean return')
    out(f'  light: a dryland crust fixing {rng(BIOCRUST_G_C_M2_YR)} g C/m2/yr could lay '
        f'{rng([x * M["CaCO3"] / M["C"] for x in BIOCRUST_G_C_M2_YR], "{:.0f}")} g of calcite at one a carbon; calcium '
        f'supplies {supply:.1f}')
    for label, n_mt in (('albatross-scale carriers', CARRIERS_N_MT_YR),
                        ('60 g carriers', tuple(f * n for f, n in zip(SMALL_CARRIERS_FOOD_MT_YR, FOOD_NITROGEN)))):
        lim = (secretion_share_limit(n_mt[0], ROCK_GT_YR['central'], SECRETION_NITROGEN[1]),
               secretion_share_limit(n_mt[1], ROCK_GT_YR['central'], SECRETION_NITROGEN[0]))
        out(f'  a binder secreted by the {label}, with all their nitrogen ({rng(n_mt)} Mt/yr) in it: at most '
            f'{rng([100 * x for x in lim])}% of the load')
    for f in (0.01, 0.03, 0.10):
        mt = f * ROCK_GT_YR['central'] * 1000.0
        out(f'  all the land at {100 * f:.0f}% calcite: {mt:.0f} Mt/yr of cement holding {mt * M["C"] / M["CaCO3"] / 1000:.3f} Gt C/yr')

    out('\n2 and 4. RATES AND TIME: measured rates against the need')
    out(f'  need at 3%: {need[0.03]:.1f} g/m2/yr at the mean return, times the concentration on a rim')
    out('  process                                        g/m2/yr        years per century of return  land a m2 can serve at 3%')
    for name, rate, cond in process_rates():
        yrs = (years_per_century_of_return(rate[1], 0.03), years_per_century_of_return(rate[0], 0.03))
        served = (concentration_served(rate[0], 0.03), concentration_served(rate[1], 0.03))
        out(f'  {name:46s} {rng(rate):13s}  {rng(yrs):27s}  {rng(served)}   [{cond}]')
    out(f'  weathering in place, calcium: {here["calcite"]:.1f} g/m2/yr keeps {100 * fraction_kept_up(here["calcite"]):.1f}% '
        f'at the mean return, {100 * fraction_kept_up(here["calcite"], 10):.2f}% at ten times, '
        f'{100 * fraction_kept_up(here["calcite"], 100):.3f}% at a hundred; the lime in the load keeps '
        f'{100 * carried_fractions(ROCK_GT_YR["central"])["calcite"]:.1f}% at any concentration')
    for key, (lo, hi) in (('ferricrete', FERRICRETE_MM_KYR), ('pedogenic calcrete', CALCRETE_MM_KYR)):
        rock_mm_kyr = MEAN_LOWERING_M_MYR
        out(f'  {key} thickens {lo:g}-{hi:g} mm a thousand years against a return of {rock_mm_kyr:.0f} mm as rock: '
            f'{rock_mm_kyr / hi:.2g}-{rock_mm_kyr / lo:.3g} times too slow (below 1, fast enough)')
    out('  rim deposits by concentration (loose soil):')
    out('    land served   builds     1 m in     10 m in     100 m in    cement at 3%   in-place calcium gives')
    for n in (1.0, 3.0, 10.0, 30.0, 100.0):
        r = deposit_mm_yr(n)
        out(f'    {n:5.0f} m2   {r:6.2f} mm/yr  {1e3 / r:8,.0f} yr  {1e4 / r:8,.0f} yr  {1e5 / r:9,.0f} yr  '
            f'{cement_need_g_m2_yr(0.03, n):7.1f} g/m2/yr   {100 * fraction_kept_up(here["calcite"], n):.2f}%')
    out('  what lasts: years for 1 m of crust to wear away, and the return such ground needs against the mean lowering')
    for key, (lo, hi) in EROSION_M_MYR.items():
        out(f'    {key:13s} {rng((lo, hi))} m/Myr: {rng((1.0 / hi, 1.0 / lo))} million years a metre; '
            f'{rng((100 * lo / MEAN_LOWERING_M_MYR, 100 * hi / MEAN_LOWERING_M_MYR))}% of the mean return')
    out(f'  root cohesion of {ROOT_COHESION_KPA["industrial_forest"][0]}-{ROOT_COHESION_KPA["natural_forest"][1]} kPa holds '
        f'{stable_thickness_m(6.8e3, 45.0, 2100.0, saturated=True):.0f}-'
        f'{stable_thickness_m(94.3e3, 45.0, 2100.0, saturated=True):.0f} m at 45 degrees saturated, within the roots\' depth')


if __name__ == '__main__':
    main()
