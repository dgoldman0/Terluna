"""Phosphorus in the sky: by tissue, the sky's own loop, the resources check, the sea's return and the falls.

Pure functions plus evaluate(). Phosphorus in g per m² of projected aerophyte area (stocks) or g per m² a
year (flows) unless labelled; sky-wide flows in t a year. P contents, resorption, retention, roosting and
hazards are design guesses given as ranges (the aerophyte growth note's, extended here). The resources
branch's phosphorus product and the provisioning branch's zonal-wind product are read from Git at pinned
commits by run.py and passed in as data.

Two states of the sky (see food_web): a mature sky holds its stock, so what must come in equals what leaves;
a growing sky also builds new bodies, whose phosphorus is new stock that must come in. The need is what
leaves plus the change of the stock aloft, less the low supply of rain and dust.
"""
from __future__ import annotations

import math

import numpy as np

from shared.constants import MOON_RADIUS
from research.studies.sky_ecology.food_web import HAZARD_PER_YEAR, SEA_DIET_P_G_PER_G_WET
from research.studies.sky_ecology.stoichiometry import (C_PER_DRY, GREEN_P_G_KG, REDFIELD_C_TO_P_MASS, RESERVE_P_G_KG,
                                                        STRUCTURE_P_G_KG, TENANT_P_G_KG)

SEA_DIET_P = SEA_DIET_P_G_PER_G_WET*1e3   # g P per kg of wet food

MOON_AREA_M2 = 4*math.pi*MOON_RADIUS**2
STATES = ('mature', 'growing')

# Design guesses (with stoichiometry's contents; the growth note's are marked there), as (low, high).
RESORPTION = (.5, .8)                      # central 0.65, the mean of terrestrial leaves (Vergutz et al. 2012)
POOL_P_MG_L = (.05, 1.)                    # design guess: dissolved P in the pools
RECRUITS_PER_ADULT_YEAR = (.002, .05)      # juveniles released per adult a year (aerophyte growth note)
JUVENILE_P_KG = (2.4, 8.8)                 # a 50-m juvenile's P from its host (aerophyte growth note, derived)
DUST_P_G_M2_YEAR = (.0003, .8)             # particle capture (aerophyte growth note)
RAIN_P_G_M2_YEAR = (.004, .1)              # rain (aerophyte growth note)
HARVEST_SHARE = (.05, .2)                  # of renewed living tissue (aerophyte resources note)
FOOD_RETURN = .99                          # share of the harvest's P people return to the sky (the resources model's q)
ROOST_SHARE = (.3, .6)                     # design guess: share of a sea-feeding flyer's excretion at its roost
CAUGHT_SHARE = (.5, .9)                    # design guess: share of a roost's droppings the host keeps
DRAPE_FACTOR = (1., 4.)                    # design guess: a fallen body covers one to four times its footprint
FLYER_TAKE_SHARE = (.1, .5)                # the food web's top_harvest_share: great flyers' share of the grazers
EARTH_SEABIRD_P_KG_YEAR = 99e6             # seabird colonies' P (Otero et al. 2018)
EARTH_SEABIRD_FOOD_KG_YEAR = 70e9          # seabirds eat 70 Mt a year, 95% CI 55.9-83.7 Mt (Brooke 2004)
WHALE_FALL = dict(mass_t=40., carbon_g=2e6, area_under_m2=50., years_of_background=2000.,
                  hectare_years=(100., 200.))   # Smith & Baco 2003
RETENTION_CASES = dict(organism_only=None, earth_forest=.98, designed_pools=.999)  # see sky_loop; 0.98 is the
                                           # recycled share of plant P demand on land (Cleveland et al. 2013)
DISSOLVED_LOSS_G_M2_YEAR = (.003, .008)    # Earth rivers' dissolved P at the Moon's runoff (lunar-cycle ecology, H4)
CLASS_SHARES = dict(sky_reef=.5, round_giant=.3, lesser=.2)   # design guess: shares of the aerophyte cover
CLASS_BODIES = dict(sky_reef='sky_reef_60m', round_giant='round_giant_1km', lesser='lesser_100m')
RESOURCES_HARVEST, RESOURCES_FOOD_RETURN = .05, .99   # the resources model's h and q


def _check(**values):
    for k, v in values.items():
        if not (isinstance(v, (int, float)) and math.isfinite(v) and v >= 0):
            raise ValueError(f'{k} must be finite and nonnegative')


def _fraction(name, value):
    if not 0 <= value <= 1:
        raise ValueError(f'{name} must lie in [0, 1]')


# ---- stocks by tissue ----------------------------------------------------------------------------

def tissue_stock(prod, community_dry_kg_m2=.1, pool_water_l_m2=(1., 5.), reserve_dry_kg_m2=.05):
    """P held per projected m², by tissue, [low, high], beside a single Redfield ratio on the same biomass."""
    green = [prod['living_dry_kg_m2']*x for x in GREEN_P_G_KG]
    envelope = [(prod['structural_dry_kg_m2']-prod['tendon_dry_kg_m2'])*x for x in STRUCTURE_P_G_KG]
    tendon = [prod['tendon_dry_kg_m2']*x for x in STRUCTURE_P_G_KG]
    reserve = [reserve_dry_kg_m2*x for x in RESERVE_P_G_KG]
    tenants = [community_dry_kg_m2*x for x in TENANT_P_G_KG]
    pools = [pool_water_l_m2[0]*POOL_P_MG_L[0]/1e3, pool_water_l_m2[1]*POOL_P_MG_L[1]/1e3]
    total = [sum(v[i] for v in (green, envelope, tendon, reserve, tenants, pools)) for i in (0, 1)]
    dry = prod['living_dry_kg_m2']+prod['structural_dry_kg_m2']+community_dry_kg_m2+reserve_dry_kg_m2
    redfield = dry*C_PER_DRY/REDFIELD_C_TO_P_MASS*1e3
    return dict(green_g_m2=green, envelope_g_m2=envelope, tendon_g_m2=tendon, reserve_g_m2=reserve,
                tenants_g_m2=tenants, pools_g_m2=pools, total_g_m2=total, dry_kg_m2=dry,
                redfield_g_m2=redfield, share_of_redfield=[t/redfield for t in total])


def growth_p(prod, growth_c, edge='low'):
    """P in a year's new bodies, g per m², at the body's own mix of living tissue and structure."""
    _check(growth=growth_c)
    i = 0 if edge == 'low' else 1
    living, structural = prod['living_dry_kg_m2'], prod['structural_dry_kg_m2']
    mix = (living*GREEN_P_G_KG[i]+structural*STRUCTURE_P_G_KG[i])/(living+structural)
    return growth_c/C_PER_DRY*mix


def throughput(prod, stock, growth_c):
    """P through a year's biomass production, g per m², [low, high], in each state, beside Redfield's for the
    same biomass carbon.

    Renewal: the renewed green tissue and structure. A mature sky adds the replacement of the bodies that die
    at the hazard (their stock's P); a growing sky the new bodies (growth_c, kg C a year, [low, high]: the
    fermented-gas body at the low edge, the photolytic one at the high).
    """
    renewal = [prod['green_tissue_dry_kg']*g+prod['structure_dry_kg']*s for g, s in zip(GREEN_P_G_KG, STRUCTURE_P_G_KG)]
    growth = [growth_p(prod, growth_c[i], e) for i, e in enumerate(('low', 'high'))]
    replacement = [stock['total_g_m2'][i]*HAZARD_PER_YEAR[i] for i in (0, 1)]
    out = dict(renewal_g_m2_year=renewal, growth_g_m2_year=growth, replacement_g_m2_year=replacement)
    for state, extra, extra_c in (('growing', growth, list(growth_c)),
                                  ('mature', replacement, [prod['standing_host_c_kg_m2']*h for h in HAZARD_PER_YEAR])):
        total = [renewal[i]+extra[i] for i in (0, 1)]
        carbon = [prod['renewal_biomass_c']+extra_c[i] for i in (0, 1)]
        out[state] = dict(total_g_m2_year=total, biomass_c_kg_m2_year=carbon,
                          redfield_g_m2_year=[c/REDFIELD_C_TO_P_MASS*1e3 for c in carbon],
                          share_of_redfield=[total[i]/(carbon[i]/REDFIELD_C_TO_P_MASS*1e3) for i in (0, 1)])
    return out


# ---- the sky's own loop ------------------------------------------------------------------------

def sky_loop(prod, stock, grazed_share, retention=None, edge='low', grazer_p=0., state='mature', growth_c=0.):
    """One m² of a class a year: where its P goes and what must come in, for one retention case and state.

    The renewed green tissue and structure carry the year's renewal P. People harvest a share of the green
    tissue and return FOOD_RETURN of its P (the resources model's food return). Grazers eat a share before
    resorption and build grazer_p into themselves; great flyers take FLYER_TAKE_SHARE of the grazers and drop
    the share they do not leave at their roost elsewhere. What remains aloft (the rest) is recycled or falls
    as litter: with retention None nothing catches what is shed, so only the resorbed share of the shed green
    tissue is kept and the frass, the grazers' remains and the structure fall (the organism alone); with a
    number, that share of the rest is kept (resorption, pools, canopy soils) and the remainder falls. Whole
    bodies fall at the hazard with their stock. A mature sky replaces them, so its stock holds; a growing sky
    also builds growth_c of new bodies, its stock changing by their P less the falls. The need is the P that
    leaves plus the change of stock; the net need takes off rain and dust at their low ends.
    """
    if state not in STATES:
        raise ValueError(f'state must be one of {STATES}')
    i = 0 if edge == 'low' else 1
    _fraction('grazed share', grazed_share)
    resorb = RESORPTION[1-i]                               # the low edge keeps more
    green_p = prod['green_tissue_dry_kg']*GREEN_P_G_KG[i]
    struct_p = prod['structure_dry_kg']*STRUCTURE_P_G_KG[i]
    renewal = green_p+struct_p
    grazed = grazed_share*green_p
    if not 0 <= grazer_p <= grazed+1e-15:
        raise ValueError('grazers cannot build more P than they eat')
    harvest = HARVEST_SHARE[i]*green_p
    if grazed+harvest > green_p:
        raise ValueError('grazing and harvest exceed the green tissue')
    harvest_not_returned = harvest*(1-FOOD_RETURN)
    leaving = grazer_p*FLYER_TAKE_SHARE[i]*(1-ROOST_SHARE[1-i])
    rest = renewal-harvest-leaving
    if retention is None:
        recycled = resorb*(green_p-grazed-harvest)
    else:
        _fraction('retention', retention)
        recycled = retention*rest
    litter = rest-recycled
    falls = stock['total_g_m2'][i]*HAZARD_PER_YEAR[i]
    if state == 'mature':
        new_p, stock_change = falls, 0.
        biomass_c = prod['renewal_biomass_c']+prod['standing_host_c_kg_m2']*HAZARD_PER_YEAR[i]
    else:
        new_p = growth_p(prod, growth_c, edge)
        stock_change = new_p-falls
        biomass_c = prod['renewal_biomass_c']+growth_c
    through = renewal+new_p
    losses = litter+falls+harvest_not_returned+leaving
    need = losses+stock_change
    low_supply = RAIN_P_G_M2_YEAR[0]+DUST_P_G_M2_YEAR[0]
    return dict(edge=edge, state=state, retention=retention, renewal_g_m2_year=renewal,
                new_bodies_g_m2_year=new_p, throughput_g_m2_year=through, biomass_c_kg_m2_year=biomass_c,
                grazed_g_m2_year=grazed, grazers_built_g_m2_year=grazer_p, frass_g_m2_year=grazed-grazer_p,
                harvest_g_m2_year=harvest, harvest_not_returned_g_m2_year=harvest_not_returned,
                tenants_leaving_g_m2_year=leaving, recycled_g_m2_year=recycled, litter_g_m2_year=litter,
                falls_g_m2_year=falls, losses_g_m2_year=losses, stock_change_g_m2_year=stock_change,
                need_g_m2_year=need, net_need_g_m2_year=max(0., need-low_supply),
                rain_g_m2_year=list(RAIN_P_G_M2_YEAR), dust_g_m2_year=list(DUST_P_G_M2_YEAR),
                rain_and_dust_at_their_high_cover_it=RAIN_P_G_M2_YEAR[1]+DUST_P_G_M2_YEAR[1] >= need,
                retained_share_of_rest=recycled/rest if rest else None,
                recycled_share=1-(need-harvest_not_returned)/through if through else None)


def nursery_supply(adult_area_m2):
    """P a year per m² of adults that recruits bring from their nursery hosts (an explicit design guess)."""
    _check(area=adult_area_m2)
    lo = RECRUITS_PER_ADULT_YEAR[0]*JUVENILE_P_KG[0]*1e3/adult_area_m2
    hi = RECRUITS_PER_ADULT_YEAR[1]*JUVENILE_P_KG[1]*1e3/adult_area_m2
    return [lo, hi]


# ---- the resources check ---------------------------------------------------------------------------

def resources_budget(throughput_kg, harvest, natural_recycling, food_return, upward_kg):
    """The resources model's sky balance: export = [(1-h)(1-r) + h(1-q)] F, deficit = export - U (kg a year)."""
    for name, f in (('harvest', harvest), ('recycling', natural_recycling), ('food return', food_return)):
        _fraction(name, f)
    _check(throughput=throughput_kg, upward=upward_kg)
    natural = (1-harvest)*(1-natural_recycling)*throughput_kg
    food = harvest*(1-food_return)*throughput_kg
    export = natural+food
    return dict(throughput_kg=throughput_kg, natural_downward_kg=natural, food_unreturned_kg=food,
                export_kg=export, upward_kg=upward_kg, additional_return_needed_kg=max(0., export-upward_kg))


def equivalent_recycling(recycled_share, harvest=RESOURCES_HARVEST):
    """The resources model's r for a loop: its r acts on the unharvested (1 - h) of the throughput."""
    _fraction('recycled share', recycled_share)
    return max(0., 1-(1-recycled_share)/(1-harvest))


def resources_check(resources, loops):
    """The resources branch's published low case, reproduced, then redone with aerophyte stoichiometry.

    Their sky: 1% of the Moon at 1,000 g C/m² a year of biomass at Redfield C:P, 5% harvested, 99% of the
    harvest's P returned, 99.9% of the rest recycled, against the low upward supply. Here the same equation
    runs with the P an aerophyte's biomass holds per unit of carbon and the recycling each loop reaches
    (loops: (state, case, edge) -> a sky_loop result for the sky reef).
    """
    central = resources['central']
    upward = resources['upward_p_kg_yr']
    row = next(r for r in resources['budget_grid'] if r['food_return'] == .99 and r['natural_recycling'] == .999
               and abs(r['upward_from_surface_p_kg_yr']-min(upward)) < 1)
    redone = resources_budget(central['throughput_p_kg_yr'], RESOURCES_HARVEST, .999, RESOURCES_FOOD_RETURN, min(upward))
    area = central['occupied_surface_footprint_m2']
    carbon_basis = central['carbon_kg_yr']/area                    # 1 kg C/m² a year
    out = dict(published=row, reproduced=redone,
               reproduces=abs(redone['additional_return_needed_kg']-row['additional_return_needed_p_kg_yr'])
               < 1e-3*row['additional_return_needed_p_kg_yr'],
               redfield_throughput_g_m2_year=central['throughput_p_kg_yr']*1e3/area,
               basis='g P per kg C of biomass production (no sugars), as their 1,000 g C/m² is all dry biomass',
               cases=[])
    for (state, name, edge), loop in loops.items():
        per_c = loop['throughput_g_m2_year']/loop['biomass_c_kg_m2_year']   # g P per kg C of biomass
        f = per_c*carbon_basis*area/1e3                                     # kg P a year at their 1,000 g C/m²
        r = equivalent_recycling(loop['recycled_share'])
        b = resources_budget(f, RESOURCES_HARVEST, r, RESOURCES_FOOD_RETURN, min(upward))
        out['cases'].append(dict(state=state, edge=edge, retention_case=name, recycled_share=loop['recycled_share'],
                                 natural_recycling=r, **b, p_per_kg_c_g=per_c,
                                 throughput_share_of_redfield=f/central['throughput_p_kg_yr'],
                                 against_published=b['additional_return_needed_kg']/row['additional_return_needed_p_kg_yr']))
    return out


# ---- the sea's return through roosting flyers ----------------------------------------------------

def flyers_for_return(target_kg_year, p_excreted_kg_year, roost_share, caught_share=1.):
    """How many sea-feeding flyers deliver a target of P a year to their roosts (and the hosts keep)."""
    _check(target=target_kg_year, p=p_excreted_kg_year)
    _fraction('roost share', roost_share)
    _fraction('caught share', caught_share)
    per = p_excreted_kg_year*roost_share*caught_share
    if per <= 0:
        raise ValueError('a flyer must deliver some P')
    return target_kg_year/per


def sea_to_land(land_loss_tg, dry_land_km2, dissolved_g_m2_year, flyers, land_share_of_routes):
    """Earth's seabird return against the Moon's land loss, and what lunar flyers would have to carry and eat.

    A flyer returns P to the land through the droppings it leaves at roosts over land, directly or through
    its host's litter: its excretion times its roost share times the routes' share of time over land.
    """
    loss = [x*1e9 for x in land_loss_tg]
    dissolved = [x*dry_land_km2*1e6/1e3 for x in dissolved_g_m2_year]
    out = dict(land_loss_kg_year=loss, dissolved_loss_kg_year=dissolved, land_share_of_routes=land_share_of_routes,
               earth_seabirds_kg_year=EARTH_SEABIRD_P_KG_YEAR,
               earth_seabirds_share_of_loss=[EARTH_SEABIRD_P_KG_YEAR/loss[1], EARTH_SEABIRD_P_KG_YEAR/loss[0]],
               earth_seabirds_share_of_dissolved=[EARTH_SEABIRD_P_KG_YEAR/dissolved[1],
                                                  EARTH_SEABIRD_P_KG_YEAR/dissolved[0]],
               targets=[])
    targets = dict(earth_seabird_scale=[EARTH_SEABIRD_P_KG_YEAR]*2, dissolved_loss=dissolved,
                   tenth_of_land_loss=[.1*x for x in loss])
    for name, target in targets.items():
        for f in flyers:
            if f['edge'] != 'low':
                continue
            per = f['p_excreted_kg_year']
            n = [flyers_for_return(target[0], per, ROOST_SHARE[1]*land_share_of_routes),
                 flyers_for_return(target[1], per, ROOST_SHARE[0]*land_share_of_routes)]
            food = [x*f['sea_food_kg_wet_day']*365.25 for x in n]
            out['targets'].append(dict(target=name, target_kg_year=target, flyer=f['flyer'], mass_kg=f['mass_kg'],
                                       flyers=n, biomass_kg=[x*f['mass_kg'] for x in n], food_kg_wet_year=food,
                                       food_over_earth_seabirds=[x/EARTH_SEABIRD_FOOD_KG_YEAR for x in food],
                                       food_p_share_of_seas_input=[food[0]*SEA_DIET_P/1e3/loss[1],
                                                                   food[1]*SEA_DIET_P/1e3/loss[0]]))
    return out


# ---- where litter and falls land ---------------------------------------------------------------------

def band_areas(edges_deg):
    """Exact shares of a sphere's area between latitude edges."""
    e = np.radians(np.asarray(edges_deg, float))
    return (np.sin(e[1:])-np.sin(e[:-1]))/2


def band_occupancy(winds, height_km):
    """Share of passive aerophyte time by 10-degree latitude band at a height, years 5-10 (zonal-wind product),
    and its concentration over the band's exact share of the Moon's area."""
    t = winds['trajectories']
    k = t['heights_km'].index(height_km)
    occ = np.asarray(t['sequence']['long_run']['occupancy'][k]).sum(axis=1)
    edges = t['occupancy_lat_edges_deg']
    area = band_areas(edges)
    return dict(height_km=height_km, lat_edges_deg=edges, time_share=occ.tolist(), area_share=area.tolist(),
                concentration=(occ/area).tolist())


def subsolar_longitude(days, sun):
    """The Sun's longitude on the wind product's fitted track at a time in days from the run's start.

    The wind code centres output n at (n + 1/2) output intervals, so at t = 0 the Sun stands half a track step
    east of the track's start; it then moves west at 360 degrees a synodic month.
    """
    return (sun['track_start_lon_deg']+sun['track_step_deg']/2-360.*np.asarray(days, float)/sun['synodic_month_days']) % 360.


def route_positions(npz, winds, height_km, first_record=20):
    """Geographic positions of the long-run routes at a height, records from `first_record` on.

    Positions are hour angle (east of the subsolar point) and latitude; record k lies k x position_every_days
    from the start, and the longitude is the subsolar longitude plus the hour angle. Each parcel carries its
    release cell's area as a weight.
    """
    t = winds['trajectories']
    lats, hours = t['release_lat_deg'], t['release_hour_deg']
    k = t['heights_km'].index(height_km)
    per = len(lats)*len(hours)
    pos = np.asarray(npz['long_run_positions'], float)/100.
    sel = pos[first_record:, k*per:(k+1)*per, :]
    days = (np.arange(first_record, pos.shape[0])*t['position_every_days'])[:, None]
    lon = (subsolar_longitude(days, winds['sun'])+sel[:, :, 0]) % 360.
    lat = sel[:, :, 1]
    w = np.repeat(np.cos(np.radians(lats)), len(hours))[None, :]*np.ones_like(lat)
    return lon, lat, w


def over_water(lon, lat, atlas):
    """Whether each position lies over water in the 28% atlas (nearest 0.25-degree cell)."""
    water = np.asarray(atlas['water_label']) > 0
    lat_grid, lon_grid = np.asarray(atlas['lat_deg']), np.asarray(atlas['lon_deg'])
    i = np.clip(np.round((lat_grid[0]-lat)/(lat_grid[0]-lat_grid[1])).astype(int), 0, len(lat_grid)-1)
    j = np.round((lon-lon_grid[0])/(lon_grid[1]-lon_grid[0])).astype(int) % len(lon_grid)
    return water[i, j]


def land_share_along(lon, lat, weight, atlas):
    """Weighted share of route positions over land, in all and record by record."""
    land = ~over_water(lon, lat, atlas)
    by_record = (weight*land).sum(axis=1)/weight.sum(axis=1)
    return dict(all=float((weight*land).sum()/weight.sum()), by_record=by_record.tolist())


def band_sea_share(atlas, edges):
    """Sea share of each latitude band in the atlas, area-weighted."""
    water = np.asarray(atlas['water_label']) > 0
    lat = np.asarray(atlas['lat_deg'])
    cw = np.cos(np.radians(lat))
    out = []
    for lo, hi in zip(edges, edges[1:]):
        m = (lat >= lo) & (lat < hi)
        out.append(float((water[m]*cw[m, None]).sum()/(cw[m].sum()*water.shape[1])))
    return out


def band_land_share(occupancy, sea_share):
    """Land share under the routes with each band's time spread evenly over its longitudes."""
    t = np.asarray(occupancy['time_share'], float)
    return float((t*(1-np.asarray(sea_share, float))).sum()/t.sum())


def deposition(occupancy, flux_t_year, sea_share=None):
    """Mean deposition per m² of ground by latitude band when a sky-wide flux falls beneath passive routes."""
    edges = occupancy['lat_edges_deg']
    area = band_areas(edges)
    rows = []
    for b, (lo, hi) in enumerate(zip(edges, edges[1:])):
        g = flux_t_year*1e6*occupancy['time_share'][b]/(MOON_AREA_M2*area[b])
        rows.append(dict(lat_deg=[lo, hi], time_share=occupancy['time_share'][b], concentration=occupancy['concentration'][b],
                         g_m2_year=g, sea_share=None if sea_share is None else sea_share[b]))
    return rows


# ---- falls ---------------------------------------------------------------------------------------

def body_fall(stock_g_m2, diameter_m, area_factor=1.):
    """P in one fallen body and the P it lays on the ground it drapes, [low, high]."""
    area = math.pi*(diameter_m/2)**2*area_factor
    p = [s*area/1e6 for s in stock_g_m2]                             # t
    per_m2 = [stock_g_m2[0]/DRAPE_FACTOR[1], stock_g_m2[1]/DRAPE_FACTOR[0]]
    return dict(diameter_m=diameter_m, area_m2=area, p_t=p, g_m2_on_ground=per_m2)


def falls_table(stocks, land_loss_g_m2_year, cover=.01):
    """Falls of a 2-km sky reef, a 1-km round giant and a 150-m body: P, fertilisation, frequency.

    stocks: name -> (stock dict from tissue_stock, diameter m, class share of the cover (a design guess), the
    class's growth rate a year [low, high]). Falls come at the adult hazard; a cover held fixed while its
    bodies spend their whole surplus on new bodies would need falls at the growth rate instead.
    """
    rows = []
    for name, (stock, d, share, growth_rate) in stocks.items():
        f = body_fall(stock['total_g_m2'], d)
        bodies = cover*MOON_AREA_M2*share/f['area_m2']
        per_year = [bodies*h for h in HAZARD_PER_YEAR]
        f.update(body=name, dry_t=stock['dry_kg_m2']*f['area_m2']/1e3, class_share_of_cover=share,
                 bodies_at_cover=bodies, falls_per_year=per_year,
                 falls_per_year_at_growth_rate=[bodies*g for g in growth_rate],
                 p_t_year=[per_year[0]*f['p_t'][0], per_year[1]*f['p_t'][1]],
                 ground_km2_year=[per_year[0]*f['area_m2']*DRAPE_FACTOR[0]/1e6,
                                  per_year[1]*f['area_m2']*DRAPE_FACTOR[1]/1e6],
                 years_of_land_loss=[f['g_m2_on_ground'][0]/land_loss_g_m2_year[1],
                                     f['g_m2_on_ground'][1]/land_loss_g_m2_year[0]])
        rows.append(f)
    return rows


def evaluate(food_web, resources, winds, winds_npz, atlas, people, flyers, water_share):
    """The phosphorus numbers for the product."""
    prods = {(r['cls'], r['diameter_m'], r['gpp_c_kg_m2_year'], r['gas_route']): r for r in food_web['production']}
    keys = dict(sky_reef_60m=('sky_reef', 60.), round_giant_1km=('round_giant', 1000.), lesser_100m=('lesser', 100.),
                lesser_150m=('lesser', 150.))
    bodies = {k: prods[(c, d, 2., 'photolysis')] for k, (c, d) in keys.items()}
    fermented = {k: prods[(c, d, 2., 'fermentation_2_1')] for k, (c, d) in keys.items()}
    growth_c = {k: [fermented[k]['growth_biomass_c'], bodies[k]['growth_biomass_c']] for k in keys}
    growth_rate = {k: [fermented[k]['growth_rate_per_year'], bodies[k]['growth_rate_per_year']] for k in keys}
    stocks = {k: tissue_stock(p) for k, p in bodies.items()}
    flows = {k: throughput(p, stocks[k], growth_c[k]) for k, p in bodies.items()}
    grazed = food_web['ranges']['tree_grazed_share']
    chains = food_web['chains']
    consumers = {e: chains['sky_reef_60m'][e]['p_balance'] for e in ('low', 'high')}
    eaten = {e: chains['eaten_base_10m'][e]['p_balance'] for e in ('low', 'high')}
    loops = []
    for key in ('sky_reef_60m', 'round_giant_1km', 'lesser_100m'):
        for state in STATES:
            for case, r in RETENTION_CASES.items():
                for i, edge in enumerate(('low', 'high')):
                    grazer_p = chains[key][edge]['p_balance']['production_p_g_m2_year']
                    loop = sky_loop(bodies[key], stocks[key], grazed[i], retention=r, edge=edge, grazer_p=grazer_p,
                                    state=state, growth_c=growth_c[key][i])
                    loop.update(body=key, retention_case=case,
                                sky_need_t_year_if_all_cover=loop['net_need_g_m2_year']*.01*MOON_AREA_M2/1e6)
                    loops.append(loop)

    def loop_of(body, state, case, edge):
        return next(l for l in loops if l['body'] == body and l['state'] == state and l['retention_case'] == case
                    and l['edge'] == edge)
    sky_need = [dict(state=s, retention_case=c, edge=e,
                     need_t_year_at_1pct=sum(CLASS_SHARES[cls]*loop_of(b, s, c, e)['net_need_g_m2_year']
                                             for cls, b in CLASS_BODIES.items())*.01*MOON_AREA_M2/1e6)
                for s in STATES for c in RETENTION_CASES for e in ('low', 'high')]
    recycled = {s: {c: [loop_of('sky_reef_60m', s, c, e)['recycled_share'] for e in ('low', 'high')]
                    for c in RETENTION_CASES} for s in STATES}
    nursery = {k: nursery_supply(a) for k, a in (('lesser_100m', math.pi*50**2), ('round_giant_1km', math.pi*500**2),
                                                 ('sky_reef_2km', math.pi*1000**2))}
    check = resources_check(resources, {(s, c, e): loop_of('sky_reef_60m', s, c, e) for s in STATES
                                        for c in RETENTION_CASES for e in ('low', 'high')})
    # Where the sky's litter lands, from the routes passive bodies follow at 10 and 20 km.
    occ = {h: band_occupancy(winds, h) for h in (10., 20.)}
    sea_band = band_sea_share(atlas, occ[10.]['lat_edges_deg']) if atlas is not None else None
    land = {}
    if atlas is not None and winds_npz is not None:
        for h in (10., 20.):
            lon, lat, w = route_positions(winds_npz, winds, h)
            snap = land_share_along(lon, lat, w, atlas)
            land[str(h)] = dict(snapshots=snap['all'], snapshot_range=[min(snap['by_record']), max(snap['by_record'])],
                                snapshots_by_record=snap['by_record'], bands=band_land_share(occ[h], sea_band))
    if '10.0' in land:
        land_route = land['10.0']['bands']
        land_basis = ('routes at 10 km over the 28% atlas: each band\'s time spread over its longitudes (the snapshots of '
                      'the routes move with the Sun\'s place and are kept as a check)')
    else:
        land_route, land_basis = 1-water_share, 'the water scenario\'s land share; the atlas grid was absent'
    litter_t = {c: [sum(CLASS_SHARES[cls]*loop_of(b, 'mature', c, e)['litter_g_m2_year'] for cls, b in CLASS_BODIES.items())
                    * .01*MOON_AREA_M2/1e6 for e in ('low', 'high')] for c in RETENTION_CASES}
    dep = {str(h): dict(per_kt=deposition(occ[h], 1e3, sea_band),
                        **{c: [deposition(occ[h], t, sea_band) for t in v] for c, v in litter_t.items()})
           for h in (10., 20.)}
    g = winds['trajectories']['sequence']['long_run']['gathering']
    gathering = dict(days=g['days'], heights_km=winds['trajectories']['heights_km'],
                     nearest_km_p50=g['nearest_km_p50'], within_cell_share=g['within_cell_share'])
    cap = people['capacity']
    sea = sea_to_land(cap['land_p_loss_tg'], cap['dry_land_km2'], DISSOLVED_LOSS_G_M2_YEAR, flyers, land_route)
    def roosting(rows, need_of, basis):
        out = []
        for row in rows:
            for f in flyers:
                if f['edge'] != 'low':
                    continue
                need = need_of(row)
                n = ([flyers_for_return(need, f['p_excreted_kg_year'], ROOST_SHARE[1], CAUGHT_SHARE[1]),
                      flyers_for_return(need, f['p_excreted_kg_year'], ROOST_SHARE[0], CAUGHT_SHARE[0])]
                     if need > 0 else [0., 0.])
                food = [x*f['sea_food_kg_wet_day']*365.25 for x in n]
                out.append(dict(basis=basis, state=row['state'], edge=row['edge'], retention_case=row['retention_case'],
                                flyer=f['flyer'], mass_kg=f['mass_kg'], need_kg_year=need, flyers=n, food_kg_wet_year=food,
                                food_over_earth_seabirds=[x/EARTH_SEABIRD_FOOD_KG_YEAR for x in food]))
        return out
    need_rows = roosting(check['cases'], lambda r: r['additional_return_needed_kg'],
                         'the resources model\'s sky at 1,000 g C/m² a year, redone')
    sky_need_flyers = roosting(sky_need, lambda r: r['need_t_year_at_1pct']*1e3,
                               'the aerophytes\' own production, classes at their design shares of 1% cover')
    per_body = []
    areas = dict(sky_reef_60m=math.pi*1000.**2, round_giant_1km=math.pi*500.**2, lesser_100m=math.pi*50.**2)
    albatross = next(f for f in flyers if f['edge'] == 'low')
    for l in loops:
        per = albatross['p_excreted_kg_year']*1e3
        need = l['net_need_g_m2_year']*areas[l['body']]
        per_body.append(dict(body=l['body'] if l['body'] != 'sky_reef_60m' else 'sky_reef_2km', state=l['state'],
                             retention_case=l['retention_case'], edge=l['edge'], need_kg_year=need/1e3,
                             flyer=albatross['flyer'],
                             roosting_flyers=[need/(per*ROOST_SHARE[1]*CAUGHT_SHARE[1]),
                                              need/(per*ROOST_SHARE[0]*CAUGHT_SHARE[0])]))
    land_loss_g = [x*1e12/(cap['dry_land_km2']*1e6) for x in cap['land_p_loss_tg']]
    falls = falls_table({'sky_reef_2km': (stocks['sky_reef_60m'], 2000., CLASS_SHARES['sky_reef'], growth_rate['sky_reef_60m']),
                         'round_giant_1km': (stocks['round_giant_1km'], 1000., CLASS_SHARES['round_giant'],
                                             growth_rate['round_giant_1km']),
                         'lesser_150m': (stocks['lesser_150m'], 150., CLASS_SHARES['lesser'], growth_rate['lesser_150m'])},
                        land_loss_g)
    return dict(
        evidence='Accounting on the aerophyte ledgers with design-guess P contents, resorption, retention, roosting and '
                 'hazards; the resources branch\'s published budget reproduced and redone; deposition from the design '
                 'run\'s passive routes and the 28% atlas. No P cycle is simulated.',
        reading_rule='Stocks in g P per m² of projected aerophyte area; flows in g P per m² a year; sky-wide flows in t '
                     'or kg a year at 1% cover unless labelled; [low, high] take every range at the end that bounds the '
                     'need. States: mature (stock held, falls replaced) and growing (the surplus builds new bodies). '
                     'recycled_share is the share of the year\'s P throughput met by recycling aloft, the resources '
                     'model\'s r once taken onto its unharvested share. Deposition is g P per m² of ground a year, the '
                     'flux spread by the routes\' time in each 10-degree band and smeared over its longitudes.',
        design_guesses=dict(green_p_g_kg=list(GREEN_P_G_KG), structure_p_g_kg=list(STRUCTURE_P_G_KG),
                            reserve_p_g_kg=list(RESERVE_P_G_KG), resorption=list(RESORPTION),
                            pool_p_mg_l=list(POOL_P_MG_L), hazard_per_year=list(HAZARD_PER_YEAR),
                            recruits_per_adult_year=list(RECRUITS_PER_ADULT_YEAR), juvenile_p_kg=list(JUVENILE_P_KG),
                            dust_p_g_m2_year=list(DUST_P_G_M2_YEAR), rain_p_g_m2_year=list(RAIN_P_G_M2_YEAR),
                            harvest_share=list(HARVEST_SHARE), roost_share=list(ROOST_SHARE),
                            caught_share=list(CAUGHT_SHARE), drape_factor=list(DRAPE_FACTOR),
                            flyer_take_share=list(FLYER_TAKE_SHARE), food_return=FOOD_RETURN,
                            class_shares_of_cover=CLASS_SHARES),
        stocks=stocks, throughput=flows, consumers=dict(sky_reef_60m=consumers, eaten_base_10m=eaten),
        loops=loops, sky_need=sky_need, roosting_flyers_per_body=per_body, nursery_g_m2_year=nursery,
        recycled_shares=recycled, resources_check=check, flyer_need=need_rows, sky_need_flyers=sky_need_flyers,
        sea_to_land=sea,
        occupancy={str(h): v for h, v in occ.items()}, band_sea_share=sea_band, land_share_of_routes=land,
        land_share_used=dict(value=land_route, basis=land_basis),
        deposition=dep, deposition_flux_t_year=litter_t, gathering=gathering, falls=falls,
        land_loss_g_m2_year=land_loss_g, whale_fall=WHALE_FALL,
        unresolved=['P contents, resorption and retention of designed tissues are guesses until tissues exist.',
                    'Roosting shares, droppings caught and flyer numbers set the sea\'s return; no flyer population is modelled.',
                    'Routes are passive at fixed heights in 3-day means; sailing, storms and steering move bodies off them.',
                    'The atlas and the climate model are taken to share their zero meridian.'])
