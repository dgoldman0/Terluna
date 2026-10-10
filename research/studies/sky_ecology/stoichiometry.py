"""Phosphorus contents shared by the food web and the phosphorus components, as (low, high) ranges.

Design guesses unless a source is named; the first two are the aerophyte growth note's, which warns
against one Redfield ratio for the whole organism.
"""
from __future__ import annotations

C_PER_DRY = .45                            # carbon share of dry plant tissue (aerophyte and aerial studies)
GREEN_P_G_KG = (1., 3.)                    # active green tissue
STRUCTURE_P_G_KG = (.05, .3)               # inert envelope, tendons and ties
RESERVE_P_G_KG = (0., .1)                  # starch night store: almost none
TENANT_P_G_KG = (6.2, 9.7)                 # insects 0.62-0.97% P of dry mass (Woods et al. 2004)
CONSUMER_C_PER_DRY = .5                    # design guess: animal tissue half carbon
CONSUMER_P_ASSIMILATION = (.5, .8)         # design guess: share of the P in food a consumer can build into itself
REDFIELD_C_TO_P_MASS = 106*12.011/30.974   # 41.1, the resources model's and the aerial study's Redfield mass ratio


def edge_parameters(edge, food_p_g_kg=GREEN_P_G_KG):
    """The food's P, the P assimilation and the consumers' own P content that bound production from below
    ('low': poor food, poor assimilation, P-rich consumers) or above ('high')."""
    if edge == 'low':
        return food_p_g_kg[0], CONSUMER_P_ASSIMILATION[0], TENANT_P_G_KG[1]
    if edge == 'high':
        return food_p_g_kg[1], CONSUMER_P_ASSIMILATION[1], TENANT_P_G_KG[0]
    raise ValueError("edge must be 'low' or 'high'")


def consumer_p_balance(grazed_c, production_c, food_p_g_kg, p_assimilation, consumer_p_g_kg):
    """Whether the food's phosphorus covers the production its carbon allows, per m² a year.

    The consumers eat green tissue (grazed_c, at food_p_g_kg) and the host's sugars (no P). Production at
    the carbon limit needs consumer_p_g_kg in bodies CONSUMER_C_PER_DRY carbon; what the food holds, built in
    at p_assimilation, sets a phosphorus limit. The shortfall is the P that droppings, caught dust and pollen
    must bring (the food that carries it is assimilated at the same efficiency) for the carbon-limited
    production.
    """
    if not (0 < p_assimilation <= 1 and consumer_p_g_kg > 0):
        raise ValueError('P assimilation in (0, 1] and a positive consumer P content required')
    food_p = grazed_c/C_PER_DRY*food_p_g_kg
    built_p = food_p*p_assimilation
    need_p = production_c/CONSUMER_C_PER_DRY*consumer_p_g_kg
    limited = min(production_c, built_p/consumer_p_g_kg*CONSUMER_C_PER_DRY)
    return dict(food_p_g_kg=food_p_g_kg, p_assimilation=p_assimilation, consumer_p_g_kg=consumer_p_g_kg,
                food_p_g_m2_year=food_p, buildable_p_g_m2_year=built_p,
                carbon_limited_need_p_g_m2_year=need_p, p_limited_production_c=limited,
                production_p_g_m2_year=limited/CONSUMER_C_PER_DRY*consumer_p_g_kg,
                p_limited_share_of_carbon_limit=limited/production_c if production_c else None,
                subsidy_for_carbon_limit_g_m2_year=max(0., need_p-built_p)/p_assimilation)
