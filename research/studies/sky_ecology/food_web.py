"""The sky's food web on aerophytes: production by class, trophic levels, flyers and the eaten base.

Pure functions plus evaluate(). Carbon in kg C, energy in J, areas in m² of projected aerophyte
area unless a name says otherwise. The aerophyte ledgers come from the coupled aerophyte model
(research/studies/aerophytes), called here for the representative bodies of each class; run.py
checks them against the committed aerophyte product. Literature values and design guesses enter
as (low, high) ranges named in RANGES, each with its source in sources.json.

Two states of the sky are kept apart. In a growing sky every body spends its surplus on new bodies, so
the cover grows; in a mature sky the bodies hold their size and number, deaths at the adult hazard are
replaced, and the surplus beyond that replacement pays for tenants' gas and food, domatia and sailing.
"""
from __future__ import annotations

import itertools
import math
from dataclasses import replace

import numpy as np

from shared.constants import JULIAN_YEAR_DAYS, MOON_RADIUS, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY
from research.studies.aerophytes import growth as aero_growth
from research.studies.aerophytes import run as aero
from research.studies.sky_ecology import stoichiometry as st

MOON_AREA_M2 = 4*math.pi*MOON_RADIUS**2
C_PER_DRY = st.C_PER_DRY             # carbon share of dry tissue, the aerophyte and aerial studies' convention
KEPT_PER_CONSTRUCTION = C_PER_DRY/aero_growth.CONSTRUCTION_C_PER_KG_DRY   # 0.809 of construction C stays in tissue
FOOD_J_PER_KG_DRY = 18e6             # the aerial study's feeding energy (aerial_ecology assumptions)
FOOD_J_PER_KG_C = FOOD_J_PER_KG_DRY/C_PER_DRY
SEA_DIET_J_PER_G_WET = 6.5e3         # seabird diets' energy content (Otero et al. 2018)
CARNIVORE_ASSIMILATION = .8          # seabirds' assimilation of animal food (Otero et al. 2018)
SEA_DIET_P_G_PER_G_WET = .006        # seabird diets' P content (Otero et al. 2018)
SUGAR_ASSIMILATION = .9              # design guess: the host's sugars and food bodies are absorbed almost whole
HAZARD_PER_YEAR = (.001, .005)       # adult death hazard (aerophyte growth note)
BIRD_FMR = (10.5, .681)              # FMR = 10.5 M^0.681 kJ/day, M in g (Nagy 2005, 95 bird species)
BIRD_FMR_RESIDUAL = (.3, 1.7)        # a species lies at 30-170% of the prediction (Nagy 2005)
SEABIRD_FMR = (9.2, .774)            # the seabird FMR in Otero et al.'s (2018) excretion model
BIRD_CARBON_PER_WET = .15            # design guess: 30% dry matter at 50% carbon
UPKEEP_EXPONENT = .72                # upkeep scales with mass to about 0.72 (Nagy 2005's mammals 0.734, birds 0.681)
KCAL_J = 4184.
INVERTEBRATE_PB = (.65, -.37)        # P/B = 0.65 M^-0.37 a year, M adult energy in kcal (Banse & Mosher 1980)
FLAKE_KG_M2 = (.05, .2)              # design guess: shed film and tissue flakes
PELLET_M = (1e-3, 3e-3)              # design guess: frass pellets, density 1,000 kg/m³
DRAG = dict(flake=1.2, pellet=.5)    # flat plate broadside, sphere

# (low, high) ranges: literature values, then design guesses (labelled in sources.json and the notes).
RANGES = dict(
    tree_grazed_share=(.12, .19),           # leaf herbivory 12-19% of foliar production in tropical forests (Metcalfe
                                            # et al. 2014); forests and shrublands lose 7.9% of above-ground NPP
                                            # (Cebrian 2004) and terrestrial plants a median 18% (Cyr & Pace 1993)
    eaten_share=(.30, .79),                 # aquatic macrophytes 30%, aquatic algae 79% (Cyr & Pace 1993 medians)
    herbivore_assimilation=(.3, .5),        # design guess anchored by insect herbivores respiring 21% of the leaf
                                            # tissue they eat (Metcalfe et al. 2014)
    poikilotherm_production_efficiency=(.2, .4),   # ~20% in poikilotherm populations (Golley 1968, in van der Meer
                                                   # 2021); 0.4 the high end
    carnivore_transfer=(.0432, .1594),      # mean 10.13% +- one s.d. of 5.81 over 140 estimates (Pauly & Christensen 1995)
    top_harvest_share=(.1, .5),             # design guess: share of the level below that great flyers take
    invertebrate_dry_mg=(1., 10.),          # design guess: the reef's grazers and predators are insect-sized
    glide_ratio=(15., 25.),                 # design guess within the sky fleet's 10 (paraglider) to 50 (sailplane)
    muscle_efficiency=(.2, .25),            # 0.25 the aerophyte navigation study's biological mover
)
COVERS = (.001, .003, .01, .02)             # 0.1%, 0.3%, the aerial study's 1%, and 2% beside the towns' 1.8%


def _check(**values):
    for k, v in values.items():
        if not (isinstance(v, (int, float)) and math.isfinite(v) and v >= 0):
            raise ValueError(f'{k} must be finite and nonnegative')


def _fraction(name, value):
    if not 0 <= value <= 1:
        raise ValueError(f'{name} must lie in [0, 1]')


def _photo_cost(gpp):
    return aero_growth.photo_carbon_per_hydrogen(gpp, aero.Traits().incident_cycle_mean_w_m2,
                                                 aero.Traits().assumed_full_solar_to_h2_lhv_efficiency)


# ---- producers -----------------------------------------------------------------------------------

def class_bodies(air, structural):
    """Representative coupled bodies at 10 km for the three classes, reference and canopy-like light.

    lesser (G02): round cellulose-class bodies 50-150 m; round giants (G03, round form): cellulose-class
    500 m, 300 MPa 1 km and, in canopy-like light only, 300 MPa 2 km; sky reefs (G03, colony form): one-
    layer rafts of 60-150 m modules. GPP 2 is the reference light and 5 the canopy-like light.
    """
    fam = aero.design_families(10, structural)
    specs = [('lesser', 'round', d) for d in (50., 100., 150.)]
    specs += [('round_giant', 'round', 500.), ('round_giant', 'round_strong', 1000.)]
    specs += [('sky_reef', 'module', d) for d in (60., 100., 150.)]
    out = []
    for gpp in (2., 5.):
        for cls, fibre, d in specs + ([('round_giant', 'round_strong', 2000.)] if gpp == 5. else []):
            if fibre == 'module':
                traits = replace(aero.colony_module(d/2, fam), gpp_c_kg_m2_year=gpp)
                case, t, cells = aero.settle_case(d, air, traits, cells=1)
            else:
                case, t, cells = aero.settle_case(d, air, replace(fam[fibre], gpp_c_kg_m2_year=gpp))
            out.append(dict(cls=cls, fibre=fibre, diameter_m=d, gpp_c_kg_m2_year=gpp, case=case, traits=t,
                            cells=cells))
    return out


def production(case, traits, carbon_per_kg_h2):
    """Biomass production per projected m² a year from the coupled carbon ledger, partitioned.

    The ledger: GPP after photolytic gas = maintenance + night gas substrate + food to tenants +
    renewal construction + surplus. Renewal construction keeps 0.45/0.556 of its carbon in tissue; the
    surplus builds new bodies, whose gas share (the growth component's construction costs at the
    stated carbon per kg of hydrogen) and reserve fix no biomass. The food to tenants is sugar the host
    gives away, counted apart from its biomass. biomass_c is the growing sky's (renewal and the new
    bodies); renewal_biomass_c is the renewed tissue and structure alone, to which a mature sky adds the
    replacement of its dead (mature_biomass).
    """
    c, m = case['carbon'], case['mass']
    area = math.pi*case['radius_m']**2*traits.cell_area_factor
    living = traits.living_dry_kg_m2
    community = living*traits.community_dry_host_fraction
    build = aero_growth.construction_per_area(case['hydrogen']['inventory_kg']/area, m['structural_dry_kg_m2'],
                                              living+community, case['storage']['minimum_stored_c_kg_m2'],
                                              carbon_per_kg_h2)
    biomass_share = (build['structure_c_kg_m2']+build['tissue_c_kg_m2'])/build['total_c_kg_m2']
    surplus = max(0., c['remaining_assimilate_c_kg_m2_year'])
    living_dry = living*traits.living_turnover_per_year
    structure_dry = c['turnover_dry_kg_m2_year']-living_dry
    renewal_c = c['turnover_construction_c_kg_m2_year']
    out = dict(
        gpp_c=c['gpp_before_h2_c_kg_m2_year'],
        photolytic_gas_not_fixed_c=c['gpp_before_h2_c_kg_m2_year']-c['gpp_after_h2_c_kg_m2_year'],
        maintenance_c=c['host_maintenance_c_kg_m2_year'],
        night_gas_substrate_c=c['night_h2_gross_substrate_c_kg_m2_year'],
        renewal_respiration_c=renewal_c*(1-KEPT_PER_CONSTRUCTION),
        green_tissue_c=living_dry*C_PER_DRY, green_tissue_dry_kg=living_dry,
        structure_c=structure_dry*C_PER_DRY, structure_dry_kg=structure_dry,
        food_to_tenants_c=c['consumer_food_c_kg_m2_year'],
        surplus_c=surplus,
        growth_biomass_c=surplus*biomass_share*KEPT_PER_CONSTRUCTION,
        growth_respiration_c=surplus*biomass_share*(1-KEPT_PER_CONSTRUCTION),
        growth_gas_and_reserve_c=surplus*(1-biomass_share),
        biomass_share_of_growth=biomass_share, gas_share_of_growth=build['gas_share'],
        construction_c_kg_m2=build['total_c_kg_m2'],
        standing_dry_kg_m2=m['structural_dry_kg_m2']+living+community,
        standing_host_c_kg_m2=(m['structural_dry_kg_m2']+living)*C_PER_DRY,
        living_dry_kg_m2=living, structural_dry_kg_m2=m['structural_dry_kg_m2'],
        tendon_dry_kg_m2=m['tendon_dry_kg_m2'])
    out['renewal_biomass_c'] = out['green_tissue_c']+out['structure_c']
    out['biomass_c'] = out['renewal_biomass_c']+out['growth_biomass_c']
    out['growth_rate_per_year'] = out['growth_biomass_c']/out['standing_host_c_kg_m2']
    out['ledger_residual_c'] = (c['gpp_after_h2_c_kg_m2_year']-out['maintenance_c']-out['night_gas_substrate_c']
                                - renewal_c-out['food_to_tenants_c']-c['remaining_assimilate_c_kg_m2_year'])
    return out


def mature_biomass(prod, hazard):
    """A mature sky's biomass production: renewal plus the replacement of the bodies that die at the hazard."""
    _fraction('hazard', hazard)
    return prod['renewal_biomass_c']+prod['standing_host_c_kg_m2']*hazard


def production_table(bodies, photolytic=True):
    """Production by class at both lights; growth gas costed by photolysis or 2.1 mol/mol fermentation."""
    rows = []
    for b in bodies:
        gpp = b['gpp_c_kg_m2_year']
        cost = _photo_cost(gpp) if photolytic else aero_growth.fermentation_carbon_per_hydrogen(2.1)
        gates = b['case']['gates']
        prod = production(b['case'], b['traits'], cost)
        prod['mature_biomass_c'] = [mature_biomass(prod, h) for h in HAZARD_PER_YEAR]
        rows.append(dict(cls=b['cls'], fibre=b['fibre'], diameter_m=b['diameter_m'], gpp_c_kg_m2_year=gpp,
                         gas_route='photolysis' if photolytic else 'fermentation_2_1', carbon_per_kg_h2=cost,
                         all_gates=bool(gates['mass_and_reserve'] and gates['solar_h2_allocation']
                                        and gates['positive_construction_surplus']), **prod))
    return rows


def sky_totals(per_m2_c, cover):
    """Mt C a year over a share `cover` of the Moon's surface."""
    _check(per_m2=per_m2_c)
    _fraction('cover', cover)
    return per_m2_c*cover*MOON_AREA_M2/1e9


# ---- trophic structure ---------------------------------------------------------------------------

def consumer_level(tissue_c, sugar_c, tissue_assimilation, sugar_assimilation, production_efficiency):
    """Consumption, assimilation, production, respiration and egestion of one consumer level, kg C.

    Tissue and sugar are assimilated at their own efficiencies. Closes exactly: consumption = egestion +
    assimilation; assimilation = production + respiration.
    """
    _check(tissue=tissue_c, sugar=sugar_c)
    for name, f in (('tissue assimilation', tissue_assimilation), ('sugar assimilation', sugar_assimilation),
                    ('production efficiency', production_efficiency)):
        _fraction(name, f)
    a = tissue_c*tissue_assimilation+sugar_c*sugar_assimilation
    p = a*production_efficiency
    food = tissue_c+sugar_c
    return dict(consumption_c=food, assimilation_c=a, production_c=p, respiration_c=a-p, egestion_c=food-a)


def invertebrate_pb(dry_mg):
    """Annual production/biomass of an invertebrate population from adult mass (Banse & Mosher 1980)."""
    _check(mass=dry_mg)
    if dry_mg <= 0:
        raise ValueError('positive mass required')
    kcal = dry_mg*1e-6*FOOD_J_PER_KG_DRY/KCAL_J
    a, b = INVERTEBRATE_PB
    return a*kcal**b


def trophic_chain(food_on_offer, edges='low'):
    """One projected m² of a class a year, with every range at the end that bounds production.

    food_on_offer: the green tissue (or edible bodies) the consumers can graze (grazable_c) and its P
    content (food_p_g_kg), the share they take (grazed_share), the sugar the host gives its tenants
    (food_to_tenants_c, without P) and what nobody grazes (shed_c). Level 2, the tenants that eat the host
    (grazers on its tissue, guards on its sugar), produces the lesser of the carbon limit and what the food's
    phosphorus can build (stoichiometry.consumer_p_balance at stoichiometry.edge_parameters, so that the low
    edge bounds production from below). Level 3 takes the carnivore transfer of that production; great flyers
    take a share of level 3. Litter is the ungrazed tissue and the shed structure with level 2's egestion.
    Unallocated assimilate after the P limit closes the accounting explicitly; its
    eventual respiration, excretion or storage is not solved or credited elsewhere.
    """
    i = 0 if edges == 'low' else 1
    r = {k: v[i] for k, v in RANGES.items()}
    grazed = food_on_offer['grazed_share']*food_on_offer['grazable_c']
    l2 = consumer_level(grazed, food_on_offer['food_to_tenants_c'], r['herbivore_assimilation'], SUGAR_ASSIMILATION,
                        r['poikilotherm_production_efficiency'])
    balance = st.consumer_p_balance(grazed, l2['production_c'], *st.edge_parameters(edges, food_on_offer['food_p_g_kg']))
    p2 = balance['p_limited_production_c']
    carbon_potential = l2['production_c']
    # P-limited tissue production must agree with the returned consumer ledger.
    # The fate of assimilate not built into tissue is unresolved: do not silently
    # assign it to respiration, litter, food for predators or retained biomass.
    l2 = dict(l2, production_c=p2, unallocated_assimilate_c=carbon_potential-p2)
    p3 = p2*r['carnivore_transfer']
    pb = invertebrate_pb(RANGES['invertebrate_dry_mg'][i])   # the low end takes the faster small bodies
    return dict(edges=edges, grazed_c=grazed, level2=l2, p_balance=balance, level2_production_c=p2,
                level3_production_c=p3, great_flyer_food_c=p3*r['top_harvest_share'], invertebrate_pb_per_year=pb,
                level2_biomass_c=p2/pb, level3_biomass_c=p3/pb,
                carbon_limited=dict(level2_production_c=carbon_potential,
                                    level3_production_c=carbon_potential*r['carnivore_transfer'],
                                    level2_biomass_c=carbon_potential/pb),
                litter_c=food_on_offer['grazable_c']-grazed+food_on_offer['shed_c']+l2['egestion_c'])


def p_limit_corners(food_on_offer):
    """The level-2 phosphorus limit's share of the carbon limit, and the P subsidy that lifts it to the carbon
    limit, over every corner of the ranges they depend on: grazed share, tissue assimilation, production
    efficiency, the food's P, P assimilation and the consumers' own P."""
    shares, subsidies = [], []
    for gs, a, pe, fp, pa, cp in itertools.product(food_on_offer['grazed_share_range'], RANGES['herbivore_assimilation'],
                                                   RANGES['poikilotherm_production_efficiency'],
                                                   food_on_offer['food_p_g_kg'], st.CONSUMER_P_ASSIMILATION,
                                                   st.TENANT_P_G_KG):
        grazed = gs*food_on_offer['grazable_c']
        level = consumer_level(grazed, food_on_offer['food_to_tenants_c'], a, SUGAR_ASSIMILATION, pe)
        b = st.consumer_p_balance(grazed, level['production_c'], fp, pa, cp)
        shares.append(b['p_limited_share_of_carbon_limit'])
        subsidies.append(b['subsidy_for_carbon_limit_g_m2_year'])
    shares = np.array(shares)
    return dict(share_range=[float(shares.min()), float(shares.max())], corners=int(shares.size),
                corners_below_half=int((shares < .5).sum()), median_share=float(np.median(shares)),
                subsidy_range_g_m2_year=[float(min(subsidies)), float(max(subsidies))])


def tree_offer(prod, edges='low'):
    """A long-lived body's offer: grazers take a tree canopy's share of its renewed green tissue."""
    i = 0 if edges == 'low' else 1
    return dict(grazable_c=prod['green_tissue_c'], grazed_share=RANGES['tree_grazed_share'][i],
                grazed_share_range=RANGES['tree_grazed_share'], food_p_g_kg=st.GREEN_P_G_KG,
                food_to_tenants_c=prod['food_to_tenants_c'], shed_c=prod['structure_c'])


def eaten_offer(row, edges='low'):
    """The eaten base's offer: grazers take a macrophyte's to an alga's share of the new bodies' living tissue;
    the regenerated-cellulose envelopes of the eaten bodies go to the decomposers with the litter."""
    i = 0 if edges == 'low' else 1
    return dict(grazable_c=row['edible_production_c_kg_m2_year'], grazed_share=RANGES['eaten_share'][i],
                grazed_share_range=RANGES['eaten_share'], food_p_g_kg=st.GREEN_P_G_KG,
                food_to_tenants_c=row['food_to_tenants_c_kg_m2_year'], shed_c=row['envelope_production_c_kg_m2_year'])


# ---- flyers --------------------------------------------------------------------------------------

def bird_fmr_j_day(mass_kg, equation=BIRD_FMR):
    """Field metabolic rate by allometry, J per day (M in g inside the equation)."""
    _check(mass=mass_kg)
    a, b = equation
    return a*(mass_kg*1e3)**b*1e3


def similarity_fmr_ratio(length_factor, basal_share, basal_exponent=UPKEEP_EXPONENT):
    """FMR of a flyer k times Earth's in every length under lunar gravity, over the Earth flyer's FMR.

    Flight work per km is weight over glide ratio, so at the same speed and hours it scales as k³ times
    the gravity ratio; upkeep scales with mass to `basal_exponent`. basal_share is the part of the Earth
    flyer's FMR that is not flight, a design input. An allometric exponent between the two factors'
    exponents gives a ratio between them, so agreement with the allometry is a consistency check.
    """
    _check(k=length_factor)
    _fraction('basal share', basal_share)
    g = MOON_SURFACE_GRAVITY/STANDARD_GRAVITY
    return basal_share*length_factor**(3*basal_exponent)+(1-basal_share)*length_factor**3*g


def flight_j_per_km(mass_kg, glide_ratio, muscle_efficiency, gravity=MOON_SURFACE_GRAVITY):
    """Metabolic energy per km of flight at a fixed glide (lift-to-drag) ratio."""
    _check(mass=mass_kg, ratio=glide_ratio)
    if glide_ratio <= 0 or not 0 < muscle_efficiency <= 1:
        raise ValueError('positive glide ratio and efficiency in (0, 1] required')
    return mass_kg*gravity*1e3/glide_ratio/muscle_efficiency


def body_spacing_m(projected_area_m2, cover):
    """Mean spacing of bodies of one size on a square lattice at a given sky cover."""
    _check(area=projected_area_m2)
    if not 0 < cover <= 1:
        raise ValueError('cover in (0, 1]')
    return math.sqrt(projected_area_m2/cover)


def flyer_table(lunar_flyers):
    """FMR, food and flight costs of the people product's lunar flyers (Earth's largest scaled six times)."""
    rows = []
    g = MOON_SURFACE_GRAVITY/STANDARD_GRAVITY
    for name, f in lunar_flyers.items():
        for i, m_t in enumerate(f['moon_mass_t']):
            m = m_t*1e3
            fmr = bird_fmr_j_day(m)
            earth = bird_fmr_j_day(f['earth_mass_kg'][i])
            k = f['moon_span_m'][i]/f['earth_span_m'][i]
            sea_food_g = fmr/(SEA_DIET_J_PER_G_WET*CARNIVORE_ASSIMILATION)
            hi_cost = flight_j_per_km(m, RANGES['glide_ratio'][0], RANGES['muscle_efficiency'][0])
            lo_cost = flight_j_per_km(m, RANGES['glide_ratio'][1], RANGES['muscle_efficiency'][1])
            rows.append(dict(
                flyer=name, edge=['low', 'high'][i], earth_mass_kg=f['earth_mass_kg'][i], mass_kg=m,
                span_m=f['moon_span_m'][i], length_factor=k,
                fmr_j_day=fmr, fmr_range_j_day=[fmr*BIRD_FMR_RESIDUAL[0], fmr*BIRD_FMR_RESIDUAL[1]],
                seabird_equation_fmr_j_day=bird_fmr_j_day(m, SEABIRD_FMR),
                earth_fmr_j_day=earth, allometric_ratio=fmr/earth,
                similarity_ratio=[similarity_fmr_ratio(k, s) for s in (.3, .7)],
                similarity_bounds=[min(k**3*g, k**(3*UPKEEP_EXPONENT)), max(k**3*g, k**(3*UPKEEP_EXPONENT))],
                sea_food_kg_wet_day=sea_food_g/1e3,
                p_excreted_kg_year=sea_food_g*SEA_DIET_P_G_PER_G_WET*JULIAN_YEAR_DAYS/1e3,
                flight_j_per_km=[lo_cost, hi_cost],
                flight_j_per_km_same_animal_on_earth=[lo_cost*STANDARD_GRAVITY/MOON_SURFACE_GRAVITY,
                                                      hi_cost*STANDARD_GRAVITY/MOON_SURFACE_GRAVITY],
                km_per_day_of_fmr=[fmr/hi_cost, fmr/lo_cost],
                biomass_c_kg=m*BIRD_CARBON_PER_WET))
    return rows


def flyers_fed(food_c_kg_year, assimilation, fmr_j_day):
    """How many flyers a supply of food carbon feeds at their FMR (their production is a percent or two)."""
    _check(food=food_c_kg_year, fmr=fmr_j_day)
    _fraction('assimilation', assimilation)
    return food_c_kg_year*FOOD_J_PER_KG_C*assimilation/(fmr_j_day*JULIAN_YEAR_DAYS)


def commute_share(mass_kg, spacing_m, fmr_j_day):
    """One flight between neighbouring bodies as a share of the flyer's daily FMR, [low, high]."""
    lo = flight_j_per_km(mass_kg, RANGES['glide_ratio'][1], RANGES['muscle_efficiency'][1])
    hi = flight_j_per_km(mass_kg, RANGES['glide_ratio'][0], RANGES['muscle_efficiency'][0])
    return [x*spacing_m/1e3/fmr_j_day for x in (lo, hi)]


# ---- grazing against filter feeding -----------------------------------------------------------

def filter_return(concentration_mg_m3, break_even_mg_m3):
    """Energy a powered filter feeder recovers per unit spent: linear in concentration, one at break-even."""
    _check(c=concentration_mg_m3, b=break_even_mg_m3)
    if break_even_mg_m3 <= 0:
        raise ValueError('positive break-even required')
    return concentration_mg_m3/break_even_mg_m3


def open_sky_food(aerial):
    """The aerial study's resident aeroplankton spread through its warm air, mg per m³ (wet mass)."""
    base = aerial['baseline']
    volume = base['warm_air_volume_km3']*1e9
    return [w*1e9/volume for w in base['microbial_wet_t_at_1pg_cell'].values()]


# ---- the eaten base: small aerophytes ----------------------------------------------------------------

def small_body(diameter_m, air, structural, free_water_kg_m2=.5, gpp=2., gate=.7):
    """The heaviest living tissue a small single-cell aerophyte carries with a 30% reserve.

    Storm gust and the colony modules' 100 Pa spare pressure; no wet fire-barrier skin (a burning
    small body is lost whole); half a kilogram of stored water per m². Bisection on living tissue.
    """
    fam = aero.design_families(10, structural)
    t0 = replace(fam['base'], name='small_aerophyte', skin_water_kg_m2=0., free_water_kg_m2=free_water_kg_m2,
                 bottom_overpressure_pa=aero.Traits().bottom_overpressure_pa+aero.SPARE_TRIM_PA_COLONY,
                 gpp_c_kg_m2_year=gpp)
    lo, hi = 0., 5.
    for _ in range(80):
        mid = (lo+hi)/2
        c = aero.coupled_case(diameter_m/2, air, replace(t0, living_dry_kg_m2=mid))
        lo, hi = (mid, hi) if c['mass']['required_lifting_mix_volume_fraction'] <= gate else (lo, mid)
    t = replace(t0, living_dry_kg_m2=lo)
    return aero.coupled_case(diameter_m/2, air, t), t


def reference_small_check(air):
    """The 10-m reference body with no stored water: its gross lift and the tissue it carries at 70%."""
    out = {}
    for label, kw in (('with_ballonet_partition', {}),
                      ('no_partition', dict(partition_count=0., gas_air_partitions=0.))):
        t0 = replace(aero.Traits(), free_water_kg_m2=0., **kw)
        lo, hi = 0., 3.
        for _ in range(80):
            mid = (lo+hi)/2
            c = aero.coupled_case(5., air, replace(t0, living_dry_kg_m2=mid))
            lo, hi = (mid, hi) if c['mass']['required_lifting_mix_volume_fraction'] <= .7 else (lo, mid)
        c = aero.coupled_case(5., air, replace(t0, living_dry_kg_m2=lo))
        out[label] = dict(living_share_of_reference=lo/aero.Traits().living_dry_kg_m2,
                          gross_lift_kg_m2=c['mass']['full_gross_supported_mass_kg']/(math.pi*25.),
                          two_thirds_d_delta_rho_kg_m2=2/3*10*air['lift_hydrogen_kg_m3'],
                          structural_dry_kg_m2=c['mass']['structural_dry_kg_m2'],
                          gas_fraction=c['mass']['required_lifting_mix_volume_fraction'])
    return out


def eaten_base(air, structural, diameters=(8., 10., 15., 20., 30., 40.)):
    """Design limits and turnover of small aerophytes that grazers eat whole.

    Turnover is the surplus over the whole body's construction carbon (gas, envelope, tissue, reserve): the
    share of the population that can be eaten each year with numbers held steady. The edible carbon is the
    living tissue kept in biomass. The envelope is regenerated-cellulose film, as the reference organism's,
    and goes to the decomposers; the gas is lost when a body is eaten; the starch store is under 1% of the
    construction carbon and is left out.
    """
    rows = []
    for d in diameters:
        for gpp in (2., 5.):
            c, t = small_body(d, air, structural, gpp=gpp)
            area = math.pi*(d/2)**2
            h2 = (c['hydrogen'].get('inventory_kg') or 0.)/area
            for route, cost in (('photolysis', _photo_cost(gpp)),
                                ('fermentation_2_1', aero_growth.fermentation_carbon_per_hydrogen(2.1))):
                build = aero_growth.construction_per_area(h2, c['mass']['structural_dry_kg_m2'],
                                                          t.living_dry_kg_m2*(1+t.community_dry_host_fraction),
                                                          c['storage']['minimum_stored_c_kg_m2'], cost)
                total = build['total_c_kg_m2']
                surplus = max(0., c['carbon']['remaining_assimilate_c_kg_m2_year'])
                edible = build['tissue_c_kg_m2']*KEPT_PER_CONSTRUCTION
                envelope = build['structure_c_kg_m2']*KEPT_PER_CONSTRUCTION
                turnover = surplus/total
                rows.append(dict(diameter_m=d, gpp_c_kg_m2_year=gpp, gas_route=route,
                                 living_dry_kg_m2=t.living_dry_kg_m2,
                                 living_share_of_reference=t.living_dry_kg_m2/aero.Traits().living_dry_kg_m2,
                                 structural_dry_kg_m2=c['mass']['structural_dry_kg_m2'],
                                 gas_fraction=c['mass']['required_lifting_mix_volume_fraction'],
                                 hydrogen_kg_m2=h2, surplus_c_kg_m2_year=surplus, construction_c_kg_m2=total,
                                 construction_shares=dict(
                                     gas=build['gas_c_kg_m2']/total, edible_tissue=edible/total, envelope=envelope/total,
                                     construction_respiration=(build['structure_c_kg_m2']+build['tissue_c_kg_m2'])
                                     * (1-KEPT_PER_CONSTRUCTION)/total,
                                     starch_store=build['reserve_c_kg_m2']/total),
                                 edible_c_kg_m2=edible, edible_share_of_construction=edible/total,
                                 envelope_c_kg_m2=envelope, turnover_per_year=turnover,
                                 edible_production_c_kg_m2_year=turnover*edible,
                                 envelope_production_c_kg_m2_year=turnover*envelope,
                                 food_to_tenants_c_kg_m2_year=c['carbon']['consumer_food_c_kg_m2_year']))
    return rows


def grazing_tolerance(prod, grazed_share):
    """Carbon a host spends regrowing grazed green tissue, and its share of the host's surplus."""
    _fraction('grazed share', grazed_share)
    regrow = grazed_share*prod['green_tissue_c']/KEPT_PER_CONSTRUCTION
    return dict(grazed_share=grazed_share, regrowth_c=regrow,
                share_of_surplus=regrow/prod['surplus_c'] if prod['surplus_c'] > 0 else None)


# ---- detritus: sky snow --------------------------------------------------------------------------

def fall_speed(air_density, areal_density_kg_m2=None, sphere_diameter_m=None, particle_density=1000.,
               gravity=MOON_SURFACE_GRAVITY):
    """Terminal fall speed of a flat flake (from its areal density) or a compact pellet (from its diameter)."""
    if air_density <= 0:
        raise ValueError('positive air density required')
    if areal_density_kg_m2 is not None:
        _check(sigma=areal_density_kg_m2)
        return math.sqrt(2*areal_density_kg_m2*gravity/(air_density*DRAG['flake']))
    _check(d=sphere_diameter_m)
    return math.sqrt(4*particle_density*gravity*sphere_diameter_m/(3*DRAG['pellet']*air_density))


def density_at(height_km, air_profile):
    """Air density at a height, log-linear between the sky-ship product's levels."""
    hs = np.array([a['height_km'] for a in air_profile], float)
    rho = np.array([a['density_kg_m3'] for a in air_profile], float)
    return float(np.exp(np.interp(height_km, hs, np.log(rho))))


def fall(start_km, air_profile, heights_km, u_m_s, areal_density_kg_m2=None, sphere_diameter_m=None, dz_km=.05):
    """Hours to fall from start_km to sea level, and the drift in a wind profile u(z), km.

    The fall speed follows the air's density at each height; the wind is linear between its levels and held
    at the end values beyond them.
    """
    if start_km <= 0:
        raise ValueError('positive start height required')
    n = max(1, int(round(start_km/dz_km)))
    z = (np.arange(n)+.5)*start_km/n
    v = np.array([fall_speed(density_at(h, air_profile), areal_density_kg_m2, sphere_diameter_m) for h in z])
    u = np.interp(z, np.asarray(heights_km, float), np.asarray(u_m_s, float))
    dt = start_km/n*1e3/v
    return dict(hours=float(dt.sum()/3600), drift_km=float((u*dt).sum()/1e3))


def drift_by_latitude(winds, air_profile, start_km, **particle):
    """The largest eastward and northward drift of a falling particle in each latitude's mean winds."""
    heights = np.array(winds['axes']['height_km'], float)
    out = {}
    for key, comp in (('eastward', 'u_mean_m_s'), ('northward', 'v_mean_m_s')):
        field = winds['ground_frame'][comp]
        drifts = []
        for j in range(len(winds['axes']['lat_deg'])):
            col = np.array([np.nan if row[j] is None else row[j] for row in field], float)
            ok = np.isfinite(col) & (heights <= start_km+2.5)
            if ok.sum() >= 2:
                drifts.append(fall(start_km, air_profile, heights[ok], col[ok], **particle)['drift_km'])
        out[f'largest_{key}_km'] = float(max(abs(x) for x in drifts))
    return out


def evaluate(air10, structural, aerial, people, wind_profile, winds, air_profile, sky_fleet):
    """The food web's numbers for the product."""
    bodies = class_bodies(air10, structural)
    photo = production_table(bodies, photolytic=True)
    ferm = production_table(bodies, photolytic=False)
    keys = ('biomass_c', 'renewal_biomass_c', 'food_to_tenants_c', 'green_tissue_c', 'structure_c', 'growth_biomass_c',
            'growth_rate_per_year')
    by_class = {}
    for cls in ('lesser', 'round_giant', 'sky_reef'):
        rows = [r for r in photo+ferm if r['cls'] == cls and r['all_gates']]
        entry = {f'{k}_kg_m2_year' if k != 'growth_rate_per_year' else k: [min(r[k] for r in rows), max(r[k] for r in rows)]
                 for k in keys}
        entry['mature_biomass_c_kg_m2_year'] = [min(r['mature_biomass_c'][0] for r in rows),
                                                max(r['mature_biomass_c'][1] for r in rows)]
        for light, gpp in (('reference', 2.), ('canopy', 5.)):
            sel = [r for r in rows if r['gpp_c_kg_m2_year'] == gpp]
            entry[f'{light}_light_biomass_c_kg_m2_year'] = [min(r['biomass_c'] for r in sel), max(r['biomass_c'] for r in sel)]
            entry[f'{light}_light_mature_biomass_c_kg_m2_year'] = [min(r['mature_biomass_c'][0] for r in sel),
                                                                   max(r['mature_biomass_c'][1] for r in sel)]
        entry['bodies'] = sorted({(r['fibre'], r['diameter_m']) for r in rows})
        by_class[cls] = entry
    totals = [dict(cover=cover, cls=cls,
                   biomass_mt_c_year=[sky_totals(x, cover) for x in v['biomass_c_kg_m2_year']],
                   reference_light_biomass_mt_c_year=[sky_totals(x, cover) for x in v['reference_light_biomass_c_kg_m2_year']],
                   reference_light_mature_biomass_mt_c_year=[sky_totals(x, cover) for x in
                                                             v['reference_light_mature_biomass_c_kg_m2_year']],
                   food_to_tenants_mt_c_year=[sky_totals(x, cover) for x in v['food_to_tenants_c_kg_m2_year']])
              for cover in COVERS for cls, v in by_class.items()]
    reference = [r for r in photo+ferm if r['all_gates'] and r['gpp_c_kg_m2_year'] == 2.]
    per_m2 = dict(growing=[min(r['biomass_c'] for r in reference), max(r['biomass_c'] for r in reference)],
                  mature=[min(r['mature_biomass_c'][0] for r in reference), max(r['mature_biomass_c'][1] for r in reference)],
                  food_to_tenants=[min(r['food_to_tenants_c'] for r in reference),
                                   max(r['food_to_tenants_c'] for r in reference)])
    sky = dict(basis='every class\'s bodies at the reference light (GPP 2), both gas routes, all gates passed',
               kg_c_m2_year=per_m2, mt_c_year_at_1pct={k: [sky_totals(x, .01) for x in v] for k, v in per_m2.items()})
    base = eaten_base(air10, structural)
    chains = {}
    for r in photo:
        key = {('lesser', 100.): 'lesser_100m', ('round_giant', 1000.): 'round_giant_1km',
               ('sky_reef', 60.): 'sky_reef_60m'}.get((r['cls'], r['diameter_m']))
        if key and r['gpp_c_kg_m2_year'] == 2.:
            chains[key] = dict(production=r, low=trophic_chain(tree_offer(r, 'low'), 'low'),
                               high=trophic_chain(tree_offer(r, 'high'), 'high'),
                               p_corners=p_limit_corners(tree_offer(r)),
                               tolerance=[grazing_tolerance(r, g) for g in RANGES['tree_grazed_share']])
    eb = next(r for r in base if r['diameter_m'] == 10. and r['gpp_c_kg_m2_year'] == 2. and r['gas_route'] == 'photolysis')
    chains['eaten_base_10m'] = dict(production=eb, low=trophic_chain(eaten_offer(eb, 'low'), 'low'),
                                    high=trophic_chain(eaten_offer(eb, 'high'), 'high'),
                                    p_corners=p_limit_corners(eaten_offer(eb)))
    levels = []
    for cover in COVERS:
        scale = cover*MOON_AREA_M2/1e9
        for key, ch in chains.items():
            prod = ch['production']
            if key == 'eaten_base_10m':
                growing = prod['edible_production_c_kg_m2_year']+prod['envelope_production_c_kg_m2_year']
                mature = None
            else:
                growing, mature = prod['biomass_c'], prod['mature_biomass_c']
            for edge in ('low', 'high'):
                t = ch[edge]
                levels.append(dict(cover=cover, body=key, edges=edge, growing_biomass_mt_c=growing*scale,
                                   mature_biomass_mt_c=None if mature is None else mature[0 if edge == 'low' else 1]*scale,
                                   level2_consumption_mt_c=t['level2']['consumption_c']*scale,
                                   level2_production_mt_c=t['level2_production_c']*scale,
                                   level2_production_carbon_limit_mt_c=t['carbon_limited']['level2_production_c']*scale,
                                   level2_biomass_mt_c=t['level2_biomass_c']*scale,
                                   level3_production_mt_c=t['level3_production_c']*scale,
                                   level3_biomass_mt_c=t['level3_biomass_c']*scale,
                                   great_flyer_food_mt_c=t['great_flyer_food_c']*scale,
                                   litter_mt_c=t['litter_c']*scale))
    flyers = flyer_table(people['lift']['flyers'])
    reef = chains['sky_reef_60m']
    spread = [1/BIRD_FMR_RESIDUAL[1], 1/BIRD_FMR_RESIDUAL[0]]
    fed = []
    for cover in COVERS:
        area = cover*MOON_AREA_M2
        for f in flyers:
            graze = [RANGES['tree_grazed_share'][i]*reef['production']['green_tissue_c']*area for i in (0, 1)]
            herb = [flyers_fed(graze[i], RANGES['herbivore_assimilation'][i], f['fmr_j_day']) for i in (0, 1)]
            pred = [flyers_fed(reef[e]['level2_production_c']*area*RANGES['top_harvest_share'][i],
                               CARNIVORE_ASSIMILATION, f['fmr_j_day']) for i, e in ((0, 'low'), (1, 'high'))]
            fed.append(dict(cover=cover, flyer=f['flyer'], edge=f['edge'], mass_kg=f['mass_kg'],
                            as_grazers=herb, as_grazers_biomass_mt_c=[n*f['biomass_c_kg']/1e9 for n in herb],
                            as_grazers_over_fmr_residual=[herb[0]*spread[0], herb[1]*spread[1]],
                            as_predators=pred, as_predators_biomass_mt_c=[n*f['biomass_c_kg']/1e9 for n in pred],
                            as_predators_over_fmr_residual=[pred[0]*spread[0], pred[1]*spread[1]]))
    commute = []
    for cover in COVERS:
        for f in flyers:
            if f['edge'] != 'low':
                continue
            for label, area in (('100m_body', math.pi*50**2), ('1km_giant', math.pi*500**2),
                                ('2km_reef', math.pi*1000**2)):
                s = body_spacing_m(area, cover)
                commute.append(dict(cover=cover, flyer=f['flyer'], bodies=label, spacing_m=s,
                                    one_flight_share_of_fmr_day=commute_share(f['mass_kg'], s, f['fmr_j_day'])))
    open_sky = open_sky_food(aerial)
    feeding10 = {str(r['drag_coefficient']): r['powered_break_even_mg_m3'] for r in aerial['feeding']
                 if r['speed_m_s'] == 10.}
    filters = dict(open_sky_aeroplankton_mg_m3=open_sky, break_even_mg_m3_at_10_m_s=feeding10,
                   energy_return={cd: [filter_return(c, b) for c in open_sky] for cd, b in feeding10.items()},
                   aerophyte_green_tissue_kg_dry_m2=aero.Traits().living_dry_kg_m2)
    hz, u = wind_profile['height_km'], wind_profile['u_m_s']
    snow = dict(flake_areal_density_kg_m2=list(FLAKE_KG_M2), pellet_diameter_m=list(PELLET_M),
                model_cell_km=winds['trajectories']['cell_km'], wind_profile=wind_profile)
    for start in (10., 20.):
        rho = density_at(start, air_profile)
        snow[f'from_{start:g}_km'] = dict(
            air_density_kg_m3=rho,
            flake_fall_m_s_at_start=[fall_speed(rho, areal_density_kg_m2=s) for s in FLAKE_KG_M2],
            pellet_fall_m_s_at_start=[fall_speed(rho, sphere_diameter_m=d) for d in PELLET_M],
            flakes=[fall(start, air_profile, hz, u, areal_density_kg_m2=s) for s in FLAKE_KG_M2],
            pellets=[fall(start, air_profile, hz, u, sphere_diameter_m=d) for d in PELLET_M],
            slowest_flakes_by_latitude=drift_by_latitude(winds, air_profile, start, areal_density_kg_m2=FLAKE_KG_M2[0]))
    gliders = {g['name']: g['earth_sea_level']['best_glide'] for g in sky_fleet['people_on_wings']['earth_gliders']}
    return dict(
        evidence='Screen on the coupled aerophyte ledgers (called live and checked against the aerophyte product), '
                 'Earth trophic efficiencies and allometries from the literature, and stated design guesses. No food '
                 'web, population or species is simulated.',
        reading_rule='kg C per m² of projected aerophyte area a year unless labelled; Mt C a year over a cover, the '
                     'share of the Moon\'s 4 pi R² surface under aerophytes. Biomass production in a growing sky counts '
                     'the renewed tissue and structure and the new bodies the surplus builds; in a mature sky the renewal '
                     'and the replacement of the bodies that die at the adult hazard; the sugar given to tenants is '
                     'counted apart. [low, high] pairs take every range at the end that bounds the quantity; p_corners '
                     'give the level-2 phosphorus limit over every corner of its ranges. FMR in J per day; counts of '
                     'flyers are the numbers a food supply feeds at their FMR, and over_fmr_residual widens them by a '
                     'species\' 30-170% spread about the allometry.',
        ranges={k: list(v) for k, v in RANGES.items()}, covers=list(COVERS), sugar_assimilation=SUGAR_ASSIMILATION,
        hazard_per_year=list(HAZARD_PER_YEAR), sky_fleet_best_glide=gliders,
        production=photo+ferm, by_class=by_class, sky=sky, totals=totals, chains=chains, levels=levels,
        eaten_base=base, reference_small_check=reference_small_check(air10),
        flyers=flyers, flyers_fed=fed, commute=commute, filters=filters, sky_snow=snow,
        earth=dict(wild_birds_mt_c=2., terrestrial_arthropods_mt_c=200., marine_arthropods_mt_c=1000.),
        unresolved=['Grazing shares, efficiencies and P/B are Earth values applied to designed organisms.',
                    'Great flyers\' diets, densities and behaviour are not modelled; the counts are energy budgets.',
                    'The food the host gives its tenants is the reference organism\'s 0.1 kg C per m² a year of sugar; '
                    'its chemistry and its tenants are open.',
                    'How a real sky divides between growing and mature bodies is open; the two states bound it.'])
