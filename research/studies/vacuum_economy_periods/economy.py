"""The vacuum economy year by year, from 2041 (year 0) to the horizon, 2341 (year 300), with the engineering of each year.

Residents follow a ramp of the programme or a schedule; past the planned 6.47 million they go on growing, the ramps at
their final rates, until the sites are full: the caves' cities and the eventual set of crater sites in towns. Technology
is held at its 2241 values after 2241.

Each year the settlement houses its newcomers in towns or cities built from their designs with that year's engineering;
renews what wears out; supplies firm power from the option that costs least; makes in its own industry what it can and
buys the rest from Earth; fetches water from the near-Earth sources, nearest first, with the tug that costs least;
sells what pays; and counts its hands against its workforce.

Stocks of capital (settlement structures and equipment, power plant by option, industrial plant by line, the tug fleet)
wear out geometrically, as national accounts treat capital: each year a share of each stock equal to one over its life
retires and is replaced at that year's engineering. Cohorts retiring all at once would echo their first build for two
centuries.

Money is the account with Earth: what the settlement pays Earth (equipment, supplies, power plant and other hardware
bought there and landed, fares of residents arriving, what partners and services on Earth keep, and the share of
residents' pay spent on Earth) against what Earth's customers pay it. The funders cover the shortfall and are repaid
from the surplus.
"""
from __future__ import annotations
from functools import lru_cache
import math

from research.studies.vacuum_economy import programme as definition
from . import designs, engineering, money, sales, technology as tech
from .designs import inputs, v

OPENING = 8
DECADE = 10                         # the window of the growth test: once repaid, no ten years together in deficit
FOUNDING_RESIDENTS = 250000.0
STREAMS = designs.STREAMS
RENEWED = ('bulk_kg', 'metal_kg', 'pressure_kg', 'equipment_kg')
LINES = ('handling', 'bulk', 'metal', 'manufactured', 'water', 'propellant')
FILL_YEARS = 10                     # the study's own: a town's lake and stocks are filled over its first ten years


@lru_cache(maxsize=4096)
def prices(tcase, date, wage):
    return engineering.Prices(tcase, date, wage)


def horizon():
    """The last year modelled: 300, 2341."""
    return int(v(inputs()['growth']['horizon_year']))


def defaults():
    s, pe, m, tr = inputs()['settlements'], inputs()['people'], inputs()['money'], inputs()['transport']
    return dict(pay_counted=v(m['pay_counted']), share_spent_on_Earth=v(m['share_of_pay_spent_on_Earth']),
                births=v(pe['births_per_resident_year']), fares_paid_by=v(pe['fares_paid_by']),
                allocation=v(s['allocation']), caves=v(s['caves_for_cities']), town_residents=v(s['town_residents']),
                reach=v(tr['reach_percentile']), share_who_work=v(pe['share_who_work']),
                water_set_aside=v(tr['water_set_aside_share']), past_plan=v(inputs()['growth']['past_plan']),
                u_policy=v(inputs()['u']['policy']), u_floor=v(inputs()['u']['rate']['floor']), u_real_yield=None,
                u_reserve_share=None, u_funders_paid_in_U=None, u_metal_share=None, u_trade_effect=None, u_spending_share=None,
                u_investor_premium=None, u_credit_line=False)


class Run:
    """One case: a scenario of demand, a technology case and a schedule of residents (a ramp of the programme, or a
    list of residents by year for growth as a result), with options that override the defaults of period_inputs.json.
    The run ends at `last_year` (the horizon by default) or, with `stop_at_exhaustion`, in the first year the water open
    to mining runs out."""

    def __init__(self, scenario, technology, ramp='constant_growth', schedule=None, last_year=None, stop_at_exhaustion=False, **options):
        self.scenario, self.tcase, self.ramp = scenario, technology, ramp
        self.opt = defaults()
        unknown = set(options) - set(self.opt)
        if unknown:
            raise ValueError('Unknown options: %s' % sorted(unknown))
        self.opt.update(options)
        self.programme = definition.load()
        self.schedule = schedule
        self.last_year = horizon() if last_year is None else int(last_year)
        self.stop_at_exhaustion = stop_at_exhaustion
        self.world = sales.world(scenario)

    def capacity(self):
        """Residents the sites can house: the caves' cities, the eventual set of crater sites in towns of the case's size,
        and the founding system's towns."""
        s = inputs()['settlements']
        return (self.opt['caves'] * s['founding_system']['residents_per_city'] + v(s['eventual_crater_sites']) * self.opt['town_residents']
                + FOUNDING_RESIDENTS * 0.2)

    def initial_town_capacity(self):
        """Residents the initial quota's towns house, with the founding system's: for reporting who lives beyond it."""
        return v(inputs()['settlements']['crater_sites']) * self.opt['town_residents'] + FOUNDING_RESIDENTS * 0.2

    def residents(self, y):
        """Residents in a year: the schedule's, or the ramp's; past the plan's completion year a ramp carries on at its
        final rate (or holds, with past_plan 'hold'). Never more than the sites house."""
        if y < OPENING:
            return 0.0
        if self.schedule is not None:
            n = self.schedule[min(y, len(self.schedule) - 1)]
        else:
            done = self.programme['calendar']['completion_year']
            n = definition.population(self.programme, self.ramp, min(y, done))
            if y > done and self.opt['past_plan'] == 'continue':
                rate = n / definition.population(self.programme, self.ramp, done - 1) - 1
                n *= (1 + rate) ** (y - done)
        return min(n, self.capacity())

    def wage(self, date):
        return v(inputs()['money']['wage_USD2041_year']) * self.world.income(date) / self.world.income(2041)

    def simulate(self):
        I, opt = inputs(), self.opt
        ind = I['industry']
        lines = ind['lines']
        lives = I['settlements']['renewal_life_years']
        found = I['settlements']['founding_system']
        sites = I['settlements']['town_sites_within_km_of_a_polar_ridge']
        town_km = sites['km'][next(i for i, s in enumerate(sites['share']) if s >= 0.5)]
        city_km = v(I['settlements']['city_km_from_a_polar_ridge'])
        city_cap = opt['caves'] * found['residents_per_city']
        town_cap = v(I['settlements']['eventual_crater_sites']) * opt['town_residents']
        initial_towns = self.initial_town_capacity()
        capacity = self.capacity()
        set_aside = opt['water_set_aside']
        # U beside the account with Earth, when the case names a policy: the year's local work is priced at U's rate, any
        # effect of paying in U scales the lines of sale below, goods that holders of U buy in the zone are charged to the
        # lines' costs, and the ledger keeps Earth's holdings, the rate, the reserve and the funders' flows.
        upol = money.policy(opt['u_policy'], dict(real_yield=opt['u_real_yield'], reserve_share_of_holdings=opt['u_reserve_share'],
                                                  funders_paid_in_U=opt['u_funders_paid_in_U'], reserve_metal_share=opt['u_metal_share']))
        uh = money.holders(dict(trade_effect_of_full_settlement=opt['u_trade_effect'], share_of_income_spent_in_zone=opt['u_spending_share'],
                                investment_threshold_spread=opt['u_investor_premium']))
        ledger = money.Ledger(upol, uh, floor=opt['u_floor'], credit_line=opt['u_credit_line']) if upol is not None else None
        lo, hi = v(I['founding']['build_years'])
        plant_life0 = v(ind['plant_design_life_years'])
        ramp_years = v(ind['local_manufacture_ramp_years'])
        grade = v(I['transport']['asteroid_water_grade'])
        land = I['transport']['depot_propellant_per_kg_landed_from_Earth']
        payload = v(I['transport']['launch_payload_kg'])
        pads_rate = v(I['transport']['launches_per_pad_year'])
        metal_per_kg = lines['manufactured']['metal_per_kg']

        housed = dict(city=0.0, town=0.0, founding=0.0)
        installed = dict.fromkeys(RENEWED, 0.0)      # the settlements' renewable stock, kg
        plant = dict.fromkeys(LINES, 0.0)            # industrial plant by line, kg
        power = {}                                   # firm capacity by option, kW
        power_life = {}
        fleet = 0.0                                  # the tug fleet's capacity, kg a year delivered
        drawn_contained = 0.0
        leg = engineering.leg_for(1.0, opt['reach'], set_aside)
        cumulative = 0.0
        feed_last, depot_cost_last, cargo_history = 0.0, None, []
        made_history, stocks_due, rows = [], {}, []
        lines_scale = 1.0
        for y in range(0, self.last_year + 1):
            date = 2041 + y
            t = tech.Year(self.tcase, date)
            wage = self.wage(date)
            # The residents' pay is in U; their work costs Earth's buyers its value at the year's opening rate.
            rate = ledger.rate if ledger is not None else 1.0
            P = prices(self.tcase, date, round(wage * rate, 2))
            N, N1 = self.residents(y), self.residents(y + 1)
            # Newcomers housed this year move in next year; the opening's residents are housed in the founding years.
            if y < OPENING:
                newcomers = self.residents(OPENING) / (hi - lo + 1) if lo <= y <= hi else 0.0
            else:
                newcomers = max(0.0, N1 - N)
            local_cap = 0.0 if y < OPENING else min(1.0, (y - OPENING) / ramp_years) * P.local_share

            # 1. Where the newcomers live, and what their housing takes to build: the founding system first, in its
            # design's proportion of cities and towns, then by the case's rule within the caves and crater sites.
            built = dict.fromkeys(STREAMS, 0.0)
            kinds = dict(city=0.0, town=0.0)
            n_found = min(newcomers, max(0.0, FOUNDING_RESIDENTS - housed['founding']))
            plan = []
            if n_found:
                c = n_found * found['cities'] * found['residents_per_city'] / FOUNDING_RESIDENTS
                plan += [('city', c, True), ('town', n_found - c, True)]
                housed['founding'] += n_found
            housing_short = 0.0
            if newcomers - n_found:
                alloc = self.allocate(newcomers - n_found, housed, city_cap, town_cap, P, date)
                housing_short = alloc.pop('beyond', 0.0)
                plan += [(k, n, k == 'city') for k, n in alloc.items()]
            for kind, n, cave in plan:
                if n:
                    b = designs.build(kind, self.tcase, date, cave)
                    for k in STREAMS:
                        built[k] += n * b[k]
                    kinds[kind] += n
            housed['city'] += kinds['city']
            housed['town'] += kinds['town']
            renewed = dict.fromkeys(STREAMS, 0.0)
            for k in RENEWED:
                life = lives[k] * (t.equipment_life_multiple if k == 'equipment_kg' else 1.0)
                renewed[k] = installed[k] / life
                installed[k] += built[k]

            # 2. What every resident needs this year.
            share_city = housed['city'] / max(1.0, housed['city'] + housed['town'])
            need_city, need_town = designs.needs('city', self.tcase, date), designs.needs('town', self.tcase, date)
            mix = lambda key: N * (share_city * need_city[key] + (1 - share_city) * need_town[key])
            # The first air and water of new housing, the towns' lakes among them, arrive over ten years.
            for due in range(y, y + FILL_YEARS):
                stocks_due[due] = stocks_due.get(due, 0.0) + built['volatiles_kg'] / FILL_YEARS
            water_need = mix('water_kg_year') + mix('air_kg_year') + stocks_due.pop(y, 0.0)
            supplies_kg = mix('earth_supplies_kg_year')
            goods_kg = mix('local_goods_kg_year') * v(I['settlements']['local_goods_mineral_share'])

            # 3. Sales, priced from the year's costs: water at the depot at last year's freight, a visitor's day at a
            # resident's power and supplies.
            depot_cost = depot_cost_last if depot_cost_last is not None else P.line_kg['water'] * 2
            stay_day = (need_city['firm_W'] / 1000.0 * P.supply_here['USD_kW_year'] + need_city['earth_supplies_kg_year'] * (P.supplies + P.landed)) / 365.25
            keep = ledger.metals_wanted() if ledger is not None else 0.0
            lines_year = sales.year(self.scenario, self.tcase, y, P, depot_cost, stay_day, feed_last, keep) if y >= OPENING else {}
            served = {k: x for k, x in lines_year.items() if x['served']}
            trade = money.trade_multiple(upol, uh, y)
            for x in served.values():
                for q in ('receipts', 'to_Earth', 'jobs', 'cost', 'kg', 'visitors', 'watts', 'power_W', 'hardware_kg_year', 'launched_kg'):
                    if q in x:
                        x[q] *= lines_scale * trade
            retained = served['metals']['retained_value'] if 'metals' in served else 0.0
            # What holders of U bought in the zone last year takes the lines' average Earth inputs, hands and power.
            zone_goods = dict(to_Earth=0.0, jobs=0.0, power_W=0.0)
            sold = sum(x['receipts'] for x in served.values())
            if ledger is not None and ledger.spent_last > 0 and sold > 0:
                zone_goods = {q: ledger.spent_last / sold * sum(x[q] for x in served.values()) for q in zone_goods}
            lift = P.lift()
            materials_kg = served['materials']['kg'] if 'materials' in served else 0.0
            propellant_sold = served['propellant']['kg'] if 'propellant' in served else 0.0
            computing = served.get('computing')
            computing_local_kg = computing['hardware_kg_year'] if computing and computing['hardware_local'] else 0.0
            computing_imported_kg = computing['hardware_kg_year'] if computing and not computing['hardware_local'] else 0.0
            sales_W = sum(x['power_W'] for x in served.values())
            if ledger is not None:
                sales_W += zone_goods['power_W']
            equipment_need = built['pressure_kg'] + built['equipment_kg'] + renewed['pressure_kg'] + renewed['equipment_kg']

            # 4. Power and the industry's work. The work asked of the industry from outside it (building and renewal,
            # everyday goods, the lines' hardware, the tugs and the power plant) sets the firm power and the plant; the
            # plant the industry builds for itself is made on the plant's margin. Power and the industry's energy are
            # solved together; the power hardware's share of the work is small, and the rounds converge.
            power_local = 0.0
            for _ in range(12):
                outside = local_cap * equipment_need + computing_local_kg + power_local
                production = dict(handling=built['handling_kg'], bulk=built['bulk_kg'] + renewed['bulk_kg'] + goods_kg,
                                  metal=built['metal_kg'] + renewed['metal_kg'] + metal_per_kg * outside + materials_kg,
                                  manufactured=outside, water=0.0, propellant=propellant_sold + materials_kg * lift['propellant_kg'])
                cargo = sum(cargo_history[-5:]) / max(1, len(cargo_history[-5:]))
                out = self.freight(water_need, production, P, leg, land, y, fleet, cargo)
                production['water'] = out['drawn_kg']
                production['propellant'] += out['propellant_kg']
                extra = out['fleet_local_kg'] + lift['launcher_kg_per_kg_year'] * materials_kg * 0.05
                production['manufactured'] += extra
                production['metal'] += metal_per_kg * extra
                # The sources and the depot carry their own power in their plant; the Moon's firm supply serves the rest.
                industry_kWh = sum(production[k] * self.kwh(lines[k], t) for k in ('handling', 'bulk', 'metal', 'manufactured')) + materials_kg * lift['kWh']
                firm_kW = (mix('firm_W') + sales_W) / 1000.0 + industry_kWh / engineering.HOURS * v(ind['industry_power_firm_factor'])
                pw = self.power(power, power_life, firm_kW, P, town_km, city_km, share_city, t)
                if abs(pw['local_kg'] - power_local) <= 1e-6 * max(1.0, pw['local_kg']):
                    break
                power_local = pw['local_kg']
            power, power_life = pw['stock'], pw['lives']
            fleet = out['fleet_after']
            outside_work = production['manufactured']
            pl = self.plant(plant, production, lines, local_cap, t, plant_life0, made_history, outside_work)
            plant = pl['stock']
            wanted = outside_work + pl['local_kg']
            goods_imported = max(0.0, wanted - pl['capacity_kg_year'])
            made = wanted - goods_imported
            production['manufactured'] = made
            production['metal'] += metal_per_kg * (made - outside_work)
            industry_kWh = sum(production[k] * self.kwh(lines[k], t) for k in ('handling', 'bulk', 'metal', 'manufactured')) + materials_kg * lift['kWh']
            made_history.append(outside_work)
            drawn_contained += out['contained_kg']
            feed_last = out['contained_kg'] / grade
            leg = engineering.leg_for(max(1.0, drawn_contained), opt['reach'], set_aside)
            # Once the water open to mining is exhausted the depot keeps last year's price for its sales; the case's later
            # years are no longer physical, and the summary names the year.
            if math.isfinite(out['water_kg_cost']):
                depot_cost_last = out['water_kg_cost']

            # 5. What is bought from Earth and what is paid for it.
            equipment_imported = (1 - local_cap) * equipment_need + goods_imported
            earth = dict(equipment=equipment_imported * P.imported_kg(), supplies=supplies_kg * (P.supplies + P.landed),
                         power=pw['to_Earth'], plant=pl['imported_kg'] * P.imported_kg(),
                         fleet=out['fleet_imported_kg'] * P.imported_kg() + out['drives_to_Earth'],
                         lines=sum(x['to_Earth'] for x in served.values()))
            if ledger is not None:
                earth['zone_goods'] = zone_goods['to_Earth']
            # People: arrivals from Earth, whose fares the settlement pays unless migrants pay their own.
            births = opt['births'] * N
            arrivals = self.residents(OPENING) if y == OPENING - 1 else (max(0.0, N1 - N - births) if y >= OPENING else 0.0)
            earth['fares'] = arrivals * P.fare() if opt['fares_paid_by'] == 'settlement' else 0.0

            # 6. Hands: the work the settlement does against the residents who work.
            jobs = dict(services=N * v(I['people']['service_jobs_per_resident']) / t.services_output_multiple,
                        operators=mix('operators'),
                        industry=sum(production[k] * lines[k]['jobs_per_Mt_year'] / 1e9 for k in LINES) / t.output_per_job,
                        lines=sum(x['jobs'] for x in served.values()))
            if ledger is not None:
                jobs['zone_goods'] = zone_goods['jobs']
            zone_hands = jobs.get('zone_goods', 0.0)
            workforce = N * opt['share_who_work'] if y >= OPENING else v(I['people']['founding_crew'])
            pay = (workforce if y >= OPENING else 0.0) * wage
            # The residents' purchases from Earth: their share of the pay at par, worth the rate in dollars.
            earth_pay_par = pay * (1.0 if opt['pay_counted'] == 'all' else opt['share_spent_on_Earth'])
            earth['pay'] = earth_pay_par * rate
            # Hands go to the settlements' own work and to the goods U's holders buy first and to the lines of sale next,
            # which were sized on last year's balance; the industry's work beyond the hands left is done from Earth by
            # teleoperation and paid there at Earth's wages.
            left = workforce - jobs['services'] - jobs['operators'] - jobs['lines'] - zone_hands
            tele = I['people']['teleoperation']
            teleoperated = max(0.0, jobs['industry'] - max(0.0, left)) if y >= OPENING else 0.0
            earth['teleoperation'] = teleoperated * wage * tele['Earth_pay_over_lunar_wage'] * tele['time_factor']
            if y >= OPENING:
                full = jobs['lines'] / lines_scale if lines_scale > 0 else jobs['lines']
                room = workforce - jobs['services'] - jobs['operators'] - zone_hands
                lines_scale = 1.0 if full <= 0 else min(1.0, max(0.0, room / full))

            receipts = sum(x['receipts'] for x in served.values())
            paid = sum(earth.values())
            u = None
            if ledger is not None:
                # The ledger sets the year's rate; the residents' purchases from Earth follow it.
                u = ledger.step(y, receipts, paid - earth['pay'], earth_pay_par, pay, receipts - earth['lines'], retained)
                earth['pay'] = u['u_earth_pay']
                paid = sum(earth.values())
            net = receipts - paid
            cumulative += net
            imported_kg = equipment_imported + supplies_kg + pw['imported_kg'] + pl['imported_kg'] + out['fleet_imported_kg'] + computing_imported_kg
            visitors = served['visits']['visitors'] if 'visits' in served else 0.0
            launches = (imported_kg + (arrivals + 2 * visitors) * v(I['people']['launched_kg_per_passenger'])) / payload
            cargo_history.append(imported_kg)
            row = dict(
                year=y, date=date, residents=N, newcomers=newcomers, arrivals=arrivals, births=births,
                housed_city=housed['city'], housed_town=housed['town'],
                housed_beyond_initial_quota=max(0.0, housed['town'] - initial_towns), sites_full=N >= capacity * (1 - 1e-9),
                receipts=receipts, paid_to_Earth=paid, to_Earth=earth, net=net, cumulative=cumulative,
                lines={k: dict(receipts=x['receipts'], to_Earth=x['to_Earth'], served=x['served'], jobs=x['jobs'],
                               **{q: x[q] for q in ('kg', 'visitors', 'price', 'watts') if q in x}) for k, x in lines_year.items()},
                production_kg=production, industry_GWh=industry_kWh / 1e6, firm_GW=firm_kW / 1e6,
                power_chosen=pw['chosen'], power_added_kW={o: sum(kw for k, kw, _ in pw['added'] if k == o) for o, _, _ in pw['added']},
                power_by_option_kW=dict(power),
                jobs=jobs, workforce=workforce, teleoperated_jobs=teleoperated,
                hands_short=max(0.0, jobs['services'] + jobs['operators'] - workforce) if y >= OPENING else 0.0,
                lines_served_share=min(1.0, lines_scale), housing_short=housing_short,
                local_share=local_cap, imported_t=imported_kg / 1000.0,
                made_for_kg=dict(settlements=local_cap * equipment_need, power=pw['local_kg'], plant=pl['local_kg'], fleet=out['fleet_local_kg'],
                                 computing=computing_local_kg, bought_finished=goods_imported),
                launches_year=launches, pads=launches / pads_rate,
                tug=out['kind'], leg_m_s=out['leg'], water_delivered_t=out['delivered_kg'] / 1000.0, contained_drawn_t=drawn_contained / 1000.0,
                exhausted=out['exhausted'], wage=wage, kWh_USD=P.kWh, fare=P.fare(), landed_USD_kg=P.landed,
                visitors=visitors, visitor_share_of_world=visitors / self.world.population(date), trade_multiple=trade)
            if u is not None:
                row.update(u)
            rows.append(row)
            if out['exhausted'] and self.stop_at_exhaustion:
                break
        return rows

    @staticmethod
    def kwh(line, t):
        k = line['kWh_kg']
        return getattr(t, k) if isinstance(k, str) else k

    def allocate(self, n, housed, city_cap, town_cap, P, date):
        """Newcomers after the founding system: by cost between towns and cities within the caves and the eventual set of
        crater sites, the initial quota's first, or by the rule the case states. The residents never outnumber what the
        sites house, so none is left beyond them."""
        rule = self.opt['allocation']
        city_room = max(0.0, city_cap - housed['city'])
        town_room = max(0.0, town_cap - housed['town'] + FOUNDING_RESIDENTS * 0.2)
        if rule == 'founding_mix':
            c = min(city_room, n * 0.8)
        elif rule == 'towns' or city_room <= 0:
            c = 0.0
        elif self.resident_cost('city', P, date) <= self.resident_cost('town', P, date):
            c = min(city_room, n)
        else:
            c = max(0.0, n - town_room) if town_room < n else 0.0
            c = min(c, city_room)
        town = n - c
        return dict(city=c, town=town, beyond=max(0.0, town - town_room))

    def resident_cost(self, kind, P, date):
        """What housing one resident of a kind costs, for choosing between them: local work at the lines' costs, the
        equipment's imported share, the first air and water, and the capital of the firm power the resident needs."""
        b = designs.build(kind, self.tcase, date, kind == 'city')
        n = designs.needs(kind, self.tcase, date)
        line = P.line_kg
        equip = b['pressure_kg'] + b['equipment_kg']
        return (b['handling_kg'] * line['handling'] + b['bulk_kg'] * line['bulk'] + b['metal_kg'] * line['metal']
                + equip * (P.local_share * line['manufactured'] + (1 - P.local_share) * P.imported_kg())
                + b['volatiles_kg'] * line['water'] * 3 + n['firm_W'] / 1000.0 * P.supply_here['USD_kW_year'] / engineering.crf(20))

    def power(self, stock, lives, firm_kW, P, town_km, city_km, share_city, t):
        """Firm capacity added this year: what the year's load needs beyond what still stands after a year's wear, from
        the cheapest option at each kind of site; the polar grid no more than the poles supply."""
        alive = {o: kw * (1 - 1 / lives[o]) for o, kw in stock.items()}
        gap = max(0.0, firm_kW - sum(alive.values()))
        added, to_earth, imported, local = [], 0.0, 0.0, 0.0
        polar = alive.get('polar_grid', 0.0)
        cap = tech.parameter('polar_grid', 'supply_cap_GW', self.tcase) * 1e6 if t.arrived('polar_grid') else 0.0
        new_lives = dict(lives)
        for site_km, share in ((city_km, share_city), (town_km, 1 - share_city)):
            need = gap * share
            if need <= 0:
                continue
            s = P.firm_supply(site_km)
            for o in sorted(s['options'], key=lambda k: s['options'][k]['USD_kW_year']):
                if need <= 0:
                    break
                take = min(need, max(0.0, cap - polar)) if o == 'polar_grid' else need
                if take <= 0:
                    continue
                if o == 'polar_grid':
                    polar += take
                x = s['options'][o]
                added.append((o, take, x['life']))
                new_lives[o] = x['life']
                imp = x['imported_kg'] * take
                imported += imp
                if o in ('fission', 'fusion'):
                    to_earth += x['capital'] * take
                else:
                    to_earth += imp * P.imported_kg()
                    local += max(0.0, x['capital'] * take - imp * P.imported_kg()) / max(1e-9, P.line_kg['manufactured'])
                need -= take
        after = dict(alive)
        for o, kw, _ in added:
            after[o] = after.get(o, 0.0) + kw
        chosen = {}
        for o, kw, _ in added:
            chosen[o] = chosen.get(o, 0.0) + kw
        return dict(added=added, to_Earth=to_earth, imported_kg=imported, local_kg=local, stock=after, lives=new_lives,
                    chosen=max(chosen, key=chosen.get) if chosen else None)

    def plant(self, stock, production, lines, local_cap, t, life0, history, outside):
        """Plant each line needs beyond what stands after a year's wear, sized a margin above the work asked of it. The
        bulk lines (handling, bulk products, metal, water and propellant) are sized to the year's work, built in the
        settlement's own factories as far as the doubling allows and bought for the rest, since their products cannot
        be bought. The factories that make equipment are sized to the larger of the year's outside work and the last
        five years' average: they start from a seed bought from Earth when local manufacture begins, then grow from their
        own output no faster than their doubling time, renewing what wears out; their capacity bounds what is made, and
        the rest is bought as finished equipment."""
        life = life0 * t.equipment_life_multiple
        grow = 2 ** (1 / t.industry_doubling_years) - 1
        margin = v(inputs()['industry']['plant_margin'])
        after, local_kg, imported_kg, capacity = {}, 0.0, 0.0, 0.0
        for k in LINES:
            before = stock[k]
            alive = before * (1 - 1 / life)
            room = max(0.0, before * (1 + grow) - alive)        # renewal of what wore out, and growth at the doubling
            ratio = lines[k]['plant_kg_per_kg_year']
            if k == 'manufactured':
                recent = history[-4:] + [outside]
                gap = max(0.0, max(sum(recent) / len(recent), outside) * ratio * margin - alive)
                if local_cap <= 0:
                    add = 0.0
                elif before <= 0:
                    add = gap                                    # the seed, bought
                    imported_kg += add
                else:
                    add = min(gap, room)
                    local_kg += add
                after[k] = alive + add
                capacity = after[k] / ratio
                continue
            gap = max(0.0, production[k] * ratio * margin - alive)
            own = min(gap, room) * local_cap if local_cap > 0 else 0.0
            local_kg += own
            imported_kg += gap - own
            after[k] = alive + gap
        return dict(stock=after, local_kg=local_kg, imported_kg=imported_kg, capacity_kg_year=capacity)

    def freight(self, water_need, production, P, leg, land, y, fleet, earth_cargo_last):
        """Water delivered to the depot for the year: the settlements' makeup and first stocks landed on the Moon, the
        propellant sold or burned lifting lunar material, and from the year the depot fuels it the descent of the cargo
        from Earth (at the last five years' average); fetched by the tug that costs least on the leg of the water being mined, nearest first. The
        fleet, in kilograms a year it can deliver, wears out geometrically and is added to where it falls short."""
        tr = inputs()['transport']
        ve = v(tr['chemical_isp_s']) * engineering.G0
        r = math.exp(v(tr['moon_to_depot_m_s']) / ve)
        descent = (1 + v(tr['moon_vehicle_kg_per_kg'])) * (r - 1)
        delivered = water_need * (1 + descent) + production['propellant']
        if y >= v(land['from_year']):
            delivered += earth_cargo_last * v(land)
        t = P.t
        alive = fleet * (1 - 1 / t.tug_life_years)
        none = dict(kind=None, leg=leg, delivered_kg=delivered, drawn_kg=0.0, contained_kg=0.0, propellant_kg=0.0,
                    fleet_local_kg=0.0, fleet_imported_kg=0.0, fleet_after=alive, drives_to_Earth=0.0,
                    water_kg_cost=float('inf'), exhausted=True)
        if leg is None:
            return none
        best = None
        for kind in ('chemical', 'electric', 'nuclear_electric'):
            x = engineering.tug(kind, leg, P)
            if x is None:
                continue
            fleet_kg_cost = P.local_share * P.line_kg['manufactured'] + (1 - P.local_share) * P.imported_kg()
            drive = 0.0
            if kind == 'electric':
                drive = x['power_W_per_kg_year'] * engineering.to2041(v(tr['electric_drive_USD2026_W'])) * t.equipment_price_multiple * engineering.crf(v(tr['electric_drive_life_years']))
            per_kg = (x['dry_kg_per_kg_year'] * fleet_kg_cost * engineering.crf(t.tug_life_years) + drive
                      + x['drawn_kg'] * P.line_kg['water'] + x['propellant_kg'] * P.line_kg['propellant'])
            if best is None or per_kg < best[1]:
                best = (x, per_kg, drive)
        if best is None:
            return none
        x, per_kg, drive = best
        drawn = delivered * x['drawn_kg']
        add = max(0.0, delivered - alive)
        add_kg = add * x['dry_kg_per_kg_year']
        local = add_kg * P.local_share if y >= OPENING else 0.0
        return dict(kind=x['kind'], leg=leg, delivered_kg=delivered, drawn_kg=drawn, contained_kg=drawn / t.water_recovered_from_rock,
                    propellant_kg=delivered * x['propellant_kg'], fleet_local_kg=local, fleet_imported_kg=add_kg - local,
                    fleet_after=alive + add, drives_to_Earth=delivered * drive, water_kg_cost=per_kg, exhausted=False)


def summary(rows, rates=(0.0, 0.01, 0.03, 0.07)):
    """The measures of one case: the funding needed (the deepest the cumulative account goes, over the whole run and by
    years 100, 200 and the last) and how much of it is needed by the opening; the net present value at each rate to years 100 and 200 and to the last year; the first year
    after the opening in which the account is repaid, the years in deficit after it and the ten-year windows after it
    whose years together are in deficit, and the first year from which it stays repaid; the year the plan's 6.47 million is reached and the year the sites are full; the residents at chosen
    years; the years short of hands; and the first year the water open to mining runs out. Under a policy of U the
    funders' measures are taken from their own flows, which carry Earth's purchases of U and the reserve, and the years
    failing each of the five tests are counted, with U's lowest rate and the years it stands at its floor."""
    nets = [r.get('funders_flow', r['net']) for r in rows]
    cum, total = [], 0.0
    for n in nets:
        total += n
        cum.append(total)
    last = rows[-1]['year']
    opening = -min(0.0, min(cum[:OPENING + 1]))
    # The payback year: the first year after the opening from which the cumulative account stays repaid to the last.
    payback = None
    for r, c in zip(reversed(rows), reversed(cum)):
        if c < 0:
            break
        if r['year'] >= OPENING:
            payback = r['year']
    repaid = next((r['year'] for r, c in zip(rows, cum) if r['year'] > OPENING and c >= 0), None)
    failing = {}
    for r in rows:
        for name, ok in r.get('u_tests', {}).items():
            failing.setdefault(name, 0)
            failing[name] += 0 if ok else 1
    horizons = sorted({h for h in (100, 200, last) if h <= last})
    npv = {str(rate): {str(h): sum(n / (1 + rate) ** y for y, n in enumerate(nets[:h + 1])) for h in horizons} for rate in rates}
    planned = definition.load()['population']['planned']
    return dict(funding_needed_USD2041=-min(0.0, min(cum)), opening_funding_USD2041=opening,
                funding_after_opening_USD2041=max(0.0, -min(cum) - opening),
                funding_needed_by_year_USD2041={str(h): -min(0.0, min(cum[:h + 1])) for h in horizons},
                NPV_USD2041=npv, last_year=last, payback_year=payback, repaid_year=repaid,
                deficit_years_after_repayment=0 if repaid is None else sum(1 for r, n in zip(rows, nets) if r['year'] > repaid and n < 0),
                u_test_failure_years=sum(1 for r in rows if not all(r.get('u_tests', {}).values())), u_test_failures=failing,
                deficit_decades_after_repayment=0 if repaid is None else sum(
                    1 for i in range(repaid + 1, len(rows) - DECADE + 1) if sum(nets[i:i + DECADE]) < 0),
                plan_reached_year=next((r['year'] for r in rows if r['residents'] >= planned * (1 - 1e-9)), None),
                sites_full_year=next((r['year'] for r in rows if r['sites_full']), None),
                residents={str(y): rows[y]['residents'] for y in (50, 100, 150, 200, 250, 300) if y <= last},
                hands_short_years=sum(1 for r in rows if r['hands_short'] > 0),
                exhausted_year=next((r['year'] for r in rows if r['exhausted']), None),
                **({} if 'u_rate' not in rows[0] else dict(u_rate_lowest=min(r['u_rate'] for r in rows),
                                                           u_years_at_floor=sum(1 for r in rows if r['u_at_floor']))))
