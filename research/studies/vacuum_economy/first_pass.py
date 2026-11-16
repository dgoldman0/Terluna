"""The first pass: the staged studies' order books, drawn up for the coupled study as its fourth case.

The functions are the staged studies' own as they stood at e49fca09, taken whole from research/studies/enclosed_returns/markets.py (population,
annual_market, case) and research/studies/complex_returns/strategies.py (profile) and written one statement to a line;
the residents of each year are an argument in place of the staged module's own schedule, which stays the default.
They read first_pass.json, the values the staged order books read, frozen from the staged inputs.
"""
from __future__ import annotations
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def inputs():
    """The rules, strategies, profiles and service scenarios of the staged order books."""
    frozen = json.loads((HERE / 'first_pass.json').read_text())
    return tuple(frozen[part]['values'] for part in ('rules', 'strategies', 'profiles', 'services'))


# From research/studies/complex_returns/strategies.py at e49fca09, profile.
def profile(q, d, name):
    return d['profiles'][q['strategies'][name]['profile']]

# From research/studies/enclosed_returns/markets.py at e49fca09, population.
def population(year, cfg):
    u = max(0.0, min(1.0, (year - cfg['ramp_start_year']) / (cfg['completion_year'] - cfg['ramp_start_year'])))
    return cfg['seed_population'] * math.exp(math.log(cfg['cap_population'] / cfg['seed_population']) * (3 * u * u - 2 * u * u * u))

# From research/studies/enclosed_returns/markets.py at e49fca09, annual_market.
def annual_market(cfg, q, d, b, name, market, rate, year, kits, paid_previous, licenses_previous, people):
    p = profile(q, d, name)
    strategy = q['strategies'][name]
    m = q['markets'][market]
    I = d['inflation_factor']
    age = year - cfg['first_productive_year']
    # Demand is held at the completion-year level during the sustaining check.
    economic_age = min(age, cfg['completion_year'] - cfg['first_productive_year'])
    growth = (1 + rate) ** economic_age
    maturity = min(1.0, b['service_demand_initial_fraction'] * (1 + b['service_demand_annual_growth']) ** min(age, 100))
    capture = cfg['seed_commercial_fraction'] * (1 + cfg['commercial_capture_growth']) ** economic_age
    capture *= min(1.0, (age + 1) / b['production_ramp_years'])
    scen = b['scenarios'][market]
    premium = q['premium']
    syn = q['synergies']
    logistics = q['logistics']
    mr = d['market_reference']
    audience = scen['audience'] * maturity
    paid = audience * scen['paid_fraction']
    licenses = scen['licenses'] * maturity
    culture = paid * scen['paid_USD2041_year'] * growth + audience * scen['advertising_and_rights_USD2041_fan_year'] * growth
    license_sales = licenses * scen['license_USD2041_year'] * growth
    raw = p['input_arguments']['physical']
    goods = raw.get('goods', 0.0) * capture
    source = raw.get('source_goods', 0.0) * capture
    refused = []
    # A necessary variable-cost test, before any machinery or outbound freight.
    if strategy['lunar_price_USD2026_kg'] < cfg['Earth_ground_service_USD2026_per_payload_kg'] and (not strategy['premium_share']):
        if goods:
            refused.append('lunar_quote_below_Earth_ground_service_cost')
        goods = 0.0
    if strategy['source_price_USD2026_kg'] < cfg['Earth_ground_service_USD2026_per_payload_kg']:
        if source:
            refused.append('source_quote_below_Earth_ground_service_cost')
        source = 0.0
    yield_delivered = cfg['delivered_quality_yield'] * (1 - cfg['return_loss_probability'])
    qualified = age >= max(0, premium['qualification_years_after_first_production'] - (syn['premium_qualification_years_saved'] if kits else 0))
    accept = min(1.0, premium['accepted_fraction'] + (syn['premium_accepted_fraction_gain'] if kits else 0)) if qualified else 0.0
    premium_input = goods * yield_delivered * strategy['premium_share']
    premium_sold = min(premium_input * accept, m['premium_demand_kg_year'] * growth)
    standard = goods * yield_delivered - premium_input + (premium_input - premium_sold) * premium['fallback_recovery_fraction']
    premium_price = premium['equivalent_product_offer_USD2026_kg'] * premium['willingness_to_pay_multiple'] * q['competition']['premium_competitor_price_factor'] * I
    goods_sales = standard * strategy['lunar_price_USD2026_kg'] * I + premium_sold * premium_price
    source_sales = source * yield_delivered * strategy['source_price_USD2026_kg'] * I
    industrial_budget = m['industrial_budget_USD2041_year'] * growth
    realized = min(1.0, industrial_budget / max(1.0, goods_sales + source_sales))
    goods *= realized
    source *= realized
    goods_sales *= realized
    source_sales *= realized
    premium_sold *= realized
    premium_scale = goods / max(1.0, raw.get('goods', 0.0)) * strategy['premium_share']
    research_capacity = max(0.0, kits * b['service_kit']['research_contracts_year'] - premium_scale * syn['qualification_research_slots_per_founding_scale'])
    research = min(scen['research_contracts'] * growth * maturity, research_capacity) * scen['research_USD2041_contract']
    professional = min(scen['professional_clients'] * growth * maturity, kits * b['service_kit']['professional_clients_year']) * scen['professional_USD2041_client_year']
    mission_research = min(research * q['competition']['mission_share_of_research_invoices'], m['mission_budget_USD2041_year'] * growth)
    nonmission_research = research * (1 - q['competition']['mission_share_of_research_invoices'])
    fuel_offer = raw.get('propellant_sales', 0.0) * capture
    fuel = min(fuel_offer, m['funded_compatible_missions_year'] * growth * logistics['delivered_fuel_kg_per_mission'])
    competitor = mr['Earth_propellant_competitor_USD2026_kg'] * I * (q['competition']['Earth_fuel_price_floor_share_of_delivered_alternative'] + (1 - q['competition']['Earth_fuel_price_floor_share_of_delivered_alternative']) * q['competition']['Earth_launch_cost_factor'])
    fuel_price = min(mr['propellant_USD2026_kg'] * I, mr['Earth_propellant_price_fraction_ceiling'] * competitor)
    fuel_sales = fuel * fuel_price
    support = fuel / logistics['delivered_fuel_kg_per_mission'] * logistics['support_invoice_USD2041_per_mission']
    mission_factor = min(1.0, max(0.0, m['mission_budget_USD2041_year'] * growth - mission_research) / max(1.0, fuel_sales + support))
    fuel *= mission_factor
    fuel_sales *= mission_factor
    support *= mission_factor
    compute_offer = raw.get('compute', 0.0) * capture
    internal_need = syn['internal_compute_W_for_culture'] + kits * syn['internal_compute_W_per_service_kit']
    internal = min(compute_offer, internal_need)
    external_compute = min(max(0.0, compute_offer - internal) * m['compute_capacity_offtake_fraction'], m['compute_contracted_capacity_W'] * growth)
    compute_sales = external_compute * mr['compute_USD2026_per_installed_W_year'] * I * m['compute_price_factor']
    institution_factor = min(1.0, m['institutional_services_budget_USD2041_year'] * growth / max(1.0, nonmission_research + professional + compute_sales))
    research = mission_research + nonmission_research * institution_factor
    professional *= institution_factor
    compute_sales *= institution_factor
    external_compute *= institution_factor
    # Destination nights use a separately priced bed stock; access fares are absent.
    visit_offer = cfg['reference_visitors_year'] * min(capture, people(year) / cfg['founding_population'])
    nights = visit_offer * cfg['visitor_stay_nights']
    destination = nights * cfg['visitor_USD2026_night'] * I
    if premium_sold:
        license_sales *= 1 + syn['additional_separate_branded_license_fraction']
    household_factor = min(1.0, m['household_budget_USD2041_year'] * growth / max(1.0, culture + license_sales + destination))
    culture *= household_factor
    license_sales *= household_factor
    paid *= household_factor
    licenses *= household_factor
    destination *= household_factor
    nights *= household_factor
    visit_offer *= household_factor
    beds = math.ceil(nights / (cfg['year_days'] * cfg['visitor_occupancy']))
    revenues = dict(materials=goods_sales, source_products=source_sales, propellant=fuel_sales, mission_support=support, computing=compute_sales, destination=destination, sports=0.45 * culture, entertainment=0.55 * culture, licenses=license_sales, research=research, professional=professional)
    v = b['variable_cost_fractions']
    variables = culture * v['culture'] + license_sales * v['licenses'] + research * v['research'] + professional * v['professional']
    digital_saving = (culture * v['culture'] + license_sales * v['licenses']) * syn['avoided_external_digital_cost_fraction'] * internal / internal_need
    acquisition = max(0.0, paid - paid_previous * (1 - b['paid_churn'])) * b['paid_acquisition_USD2041'] + max(0.0, licenses - licenses_previous * (1 - b['license_renewal_churn'])) * b['license_acquisition_USD2041']
    jobs = b['cultural_platform']['staff'] + kits * b['service_kit']['staff']
    premium_jobs = premium_scale * premium['staff_per_founding_scale']
    fuel_scale = fuel / max(1.0, raw.get('propellant_sales', 0.0))
    logistic_jobs = fuel_scale * logistics['support_staff_per_founding_scale']
    service_payroll = (jobs + premium_jobs + logistic_jobs) * b['service_staff_USD2041_year']
    fixed = b['cultural_platform']['fixed_purchases_USD2041_year'] + kits * b['service_kit']['fixed_purchases_USD2041_year']
    quoted_stock = b['cultural_platform']['capital_USD2041'] + kits * b['service_kit']['capital_USD2041']
    allowance_renew = b['cultural_platform']['capital_USD2041'] / b['cultural_platform']['life_years'] + kits * b['service_kit']['capital_USD2041'] / b['service_kit']['life_years']
    supplemental = premium_scale * premium['annual_development_USD2041_per_founding_scale'] + goods * strategy['premium_share'] * premium['incremental_cost_USD2041_attempted_kg'] + support * logistics['support_variable_cost_fraction'] + fuel_scale * logistics['support_fixed_USD2041_per_founding_scale_year']
    extra_stock = premium_scale * premium['installed_capital_USD2041_per_founding_scale'] + fuel_scale * logistics['support_capital_USD2041_per_founding_scale']
    extra_renew = premium_scale * premium['installed_capital_USD2041_per_founding_scale'] / premium['life_years'] + fuel_scale * logistics['support_capital_USD2041_per_founding_scale'] / logistics['life_years']
    earth_commercial = compute_sales * mr['compute_Earth_service_fraction'] + destination * cfg['visitor_Earth_service_fraction']
    return dict(goods=goods, source_goods=source, propellant_sales=fuel, compute=internal + external_compute, service_power_W=b['cultural_platform']['power_W'] + kits * b['service_kit']['power_W'], service_jobs=jobs, premium_power_W=goods * strategy['premium_share'] * premium['additional_mean_W_per_attempted_kg_year'], premium_jobs=premium_jobs, logistics_jobs=logistic_jobs, visitor_beds=beds, visitor_nights_year=nights, visitors_year=visit_offer, service_kit_count=kits, service_installed_capital_USD2041=quoted_stock + extra_stock, service_replacement_allowance_USD2041=allowance_renew + extra_renew, commercial_operating_USD2041=variables - digital_saving + acquisition + service_payroll + fixed + supplemental + earth_commercial, commercial_payroll_USD2041=service_payroll, gross_USD2041=sum(revenues.values()), revenue_by_sector_USD2041=revenues, audience=audience, paid_accounts=paid, licenses=licenses, external_compute_W=external_compute, internal_compute_W=internal, demand_growth_factor=growth, mission_research_USD2041=mission_research, industrial_budget_USD2041=industrial_budget, household_budget_USD2041=m['household_budget_USD2041_year'] * growth, institutional_budget_USD2041=m['institutional_services_budget_USD2041_year'] * growth, mission_budget_USD2041=m['mission_budget_USD2041_year'] * growth, industrial_realization=realized, mission_realization=mission_factor, institutional_realization=institution_factor, household_realization=household_factor, refused_export_offers=refused, premium_qualified=qualified)

# From research/studies/enclosed_returns/markets.py at e49fca09, case.
def case(cfg, q, d, b, name, market, growth_name, rate, people=None):
    """The orders and residents of one staged order book. `people(year)` gives the residents of a year; the staged
    schedule where it is not given."""
    people = people or (lambda year: population(year, cfg))
    kits = 1
    pending = []
    paid = licenses = 0.0
    rows = []
    last_stock = 0.0
    for year in range(cfg['first_productive_year'], cfg['last_year'] + 1):
        age = year - cfg['first_productive_year']
        kits += sum((n for due, n in pending if due == year))
        pending = [(y, n) for y, n in pending if y > year]
        maturation = min(1.0, b['service_demand_initial_fraction'] * (1 + b['service_demand_annual_growth']) ** min(age, 100))
        factor = (1 + rate) ** min(age, cfg['completion_year'] - cfg['first_productive_year'])
        scenario = b['scenarios'][market]
        wanted = math.ceil(max(scenario['research_contracts'] / b['service_kit']['research_contracts_year'], scenario['professional_clients'] / b['service_kit']['professional_clients_year']) * maturation * factor)
        staff_limit = max(1, math.floor((cfg['service_staff_fraction'] * people(year) - b['cultural_platform']['staff']) / b['service_kit']['staff']))
        order = int(kits + sum((n for _, n in pending)) < min(wanted, staff_limit) and year <= cfg['completion_year'] - cfg['commissioning_years'])
        if order:
            pending.append((year + cfg['commissioning_years'], 1))
        a = annual_market(cfg, q, d, b, name, market, rate, year, kits, paid, licenses, people)
        paid = a['paid_accounts']
        licenses = a['licenses']
        stock = a['service_installed_capital_USD2041']
        a.update(year=year, population=people(year), commissioned_people=0.0 if year == cfg['first_productive_year'] else people(year) - people(year - 1), new_residential_people=people(year + 2) - people(year + 1), service_new_capital_USD2041=max(0.0, stock - last_stock), service_kit_orders=order)
        last_stock = stock
        rows.append(a)
    return dict(name=f'{name}:{market}:{growth_name}', portfolio=name, market=market, market_growth=rate, market_growth_name=growth_name, annual=rows)
