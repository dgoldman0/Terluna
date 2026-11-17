"""The engineering of each year, chosen by cost: the prices of what is bought from Earth and landed, the cost of work done
on the Moon by each line of the industry, the supply of firm power at a site, the cost of carrying people, the lift of
lunar material to the depot, and the tugs that fetch water from the near-Earth sources.

Costs here are the study's prices for choosing: what an option pays to Earth together with the hands it takes, valued at
the wage, and its plant and power at their own cost, annuitized at the programme's 3%. The account with Earth itself is
kept in economy.py, which counts only the money that crosses between Earth and the settlement.
"""
from __future__ import annotations
from functools import lru_cache
import gzip
import json
import math
from pathlib import Path

from . import technology as tech
from .designs import inputs, v

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
REACH = ROOT / 'resources/source_access/results/reach.json.gz'
G0 = 9.80665
HOURS = 8766.0
RATE = 0.03


def crf(years, rate=RATE):
    """The capital recovery factor: the yearly payment that repays one dollar over a life at the rate."""
    return rate / (1 - (1 + rate) ** -years)


def to2041(usd2026):
    return usd2026 * v(inputs()['money']['USD2026_to_USD2041'])


class Prices:
    """The prices of one year in one technology case, in dollars of 2041: launch and landing, what equipment and
    electronics cost on Earth, the wage, and from them the cost of work on the Moon by line and of firm power."""

    def __init__(self, case, year, wage):
        self.case, self.year, self.t = case, year, tech.Year(case, year)
        m = inputs()['money']
        t = self.t
        self.launch = to2041(t.launch_USD2026_kg)
        self.landed = self.launch * t.landed_multiple
        self.equipment = to2041(v(m['equipment_USD2026_kg'])) * t.equipment_price_multiple
        self.supplies = to2041(v(m['supplies_USD2026_kg']))
        self.wage = wage
        self.output_per_job = t.output_per_job
        self.local_share = 1 - self.import_share()
        self._solve()

    def import_share(self):
        """The share of new equipment's mass bought from Earth once local production runs: the roadmap's path, held at
        the electronics' share until chips are made on the Moon."""
        floor = tech.parameter('chips_on_moon', 'electronics_share_of_equipment_mass', self.case if self.case != 'static' else 'central')
        share = self.t.equipment_import_share
        return share if self.t.arrived('chips_on_moon') else max(share, floor)

    def imported_kg(self):
        """A kilogram of equipment bought on Earth and landed on the Moon."""
        return self.equipment + self.landed

    def _solve(self):
        """Cost of a kilogram from each line, and of a firm kilowatt-year, solved together: lines use power and plant,
        and power uses manufactured hardware."""
        lines = inputs()['industry']['lines']
        life = v(inputs()['industry']['plant_design_life_years']) * self.t.equipment_life_multiple
        kwh_price = 0.05                       # first guess, dollars a kWh
        cost = {k: 1.0 for k in lines}
        for _ in range(30):
            plant_kg_cost = self.local_share * cost['manufactured'] + (1 - self.local_share) * self.imported_kg()
            new = {}
            for name, line in lines.items():
                kwh = line['kWh_kg']
                kwh = getattr(self.t, kwh) if isinstance(kwh, str) else kwh
                jobs = line['jobs_per_Mt_year'] / 1e9 / self.output_per_job
                c = kwh * kwh_price + line['plant_kg_per_kg_year'] * plant_kg_cost * crf(life) + jobs * self.wage
                if name == 'manufactured':
                    c += line['metal_per_kg'] * cost['metal']
                new[name] = c
            cost = new
            supply = self.firm_supply(site_km=None, kg_cost=cost['manufactured'])
            new_price = supply['USD_kW_year'] / HOURS
            if abs(new_price - kwh_price) < 1e-6 * max(1.0, kwh_price):
                kwh_price = new_price
                break
            kwh_price = new_price
        self.kWh = kwh_price
        self.line_kg = cost
        self.supply_here = supply

    def firm_supply(self, site_km, kg_cost=None):
        """The cost of a firm kilowatt for a year at a site, by each option available in the year, and the cheapest.
        `site_km` is the site's distance from a polar ridge (None for a site the polar grid does not reach)."""
        t, p = self.t, inputs()['power']
        made = self.line_kg['manufactured'] if kg_cost is None else kg_cost
        bought = self.imported_kg()
        share = t.array_local_share
        life = v(p['array_design_life_years']) * t.equipment_life_multiple
        rt = t.storage_round_trip

        def arrays(kW):
            kg = kW * 1000.0 / t.array_W_kg
            return kg * (share * made + (1 - share) * bought)

        def storage(kWh):
            kg = kWh * 1000.0 / t.storage_Wh_kg
            s = v(p['storage_local_share'])
            return kg * (s * made + (1 - s) * bought)

        out = {}
        night = v(p['night_hours'])
        out['solar_with_storage'] = dict(capital=arrays(1 + 1 / rt) + storage(night), life=life,
                                         imported_kg=(1 + 1 / rt) * 1000 / t.array_W_kg * (1 - share) + night * 1000 / t.storage_Wh_kg * (1 - v(p['storage_local_share'])))
        # A station's watts per kilogram are grid watts with the beam's whole chain; the constellation carries five times
        # the load because each station serves a site only part of the time.
        station_kg = v(p['beam_constellation_factor']) / t.beam_station_W_kg * 1000.0
        out['beamed'] = dict(capital=station_kg * (share * made + (1 - share) * bought) + v(p['beam_receiver_kg_W']) * 1000 * made,
                             life=v(p['beam_station_life_years']) * t.equipment_life_multiple, imported_kg=station_kg * (1 - share))
        fission_kg = t.fission_kg_kWe
        price = t.fission_USD2024_kWe * v(inputs()['money']['USD2024_to_USD2026']) * v(inputs()['money']['USD2026_to_USD2041'])
        out['fission'] = dict(capital=price + fission_kg * self.landed, life=t.fission_life_years, imported_kg=fission_kg)
        if t.arrived('fusion'):
            kg = tech.parameter('fusion', 'kg_kWe', self.case)
            price = tech.parameter('fusion', 'USD2024_kWe', self.case) * v(inputs()['money']['USD2024_to_USD2026']) * v(inputs()['money']['USD2026_to_USD2041'])
            out['fusion'] = dict(capital=price + kg * self.landed, life=tech.parameter('fusion', 'life_years', self.case), imported_kg=kg)
        if t.arrived('polar_grid') and site_km is not None:
            # Lines designed for a fixed share lost: their conductor grows with the square of the length.
            dark = tech.parameter('polar_grid', 'ridge_darkness_hours', self.case)
            sent = 1 / (1 - tech.parameter('polar_grid', 'line_loss', self.case))
            line_kg = tech.parameter('polar_grid', 'line_kg_kW_km2', self.case) * site_km ** 2
            ridge_kW = sent * (1 + dark / (708 - dark) / rt)
            out['polar_grid'] = dict(capital=arrays(ridge_kW) + storage(sent * dark) + line_kg * made, life=life,
                                     imported_kg=ridge_kW * 1000 / t.array_W_kg * (1 - share))
        if t.arrived('equatorial_belt'):
            belt_km = 2730.0                    # a quarter of the equator, the farthest a site is from the lit half
            line_kg = tech.parameter('polar_grid', 'line_kg_kW_km2', self.case) * belt_km ** 2
            out['equatorial_belt'] = dict(capital=arrays(2.2) + line_kg * made, life=life, imported_kg=2.2 * 1000 / t.array_W_kg * (1 - share))
        upkeep = v(p['upkeep_share_year'])
        for o in out.values():
            o['USD_kW_year'] = o['capital'] * (crf(o['life']) + upkeep)
        best = min(out, key=lambda k: out[k]['USD_kW_year'])
        return dict(options=out, chosen=best, USD_kW_year=out[best]['USD_kW_year'])

    def fare(self):
        """A passenger carried one way between Earth and the Moon: the mass launched from Earth landed, and the cabin
        moved between low Earth orbit and the Moon at the landing's share of the cost."""
        pp = inputs()['people']
        return v(pp['launched_kg_per_passenger']) * self.landed + v(pp['cabin_kg_per_passenger']) * self.launch * (self.t.landed_multiple - 1)

    def lift(self):
        """Lifting a kilogram of lunar material to the depot: by rocket on depot propellant, or by mass driver once it
        has arrived. Returns the propellant per kilogram (rocket) and the energy (mass driver)."""
        tr = inputs()['transport']
        if self.t.arrived('mass_driver'):
            first = tech.arrival('mass_driver', self.case)
            a, b = tech.parameter('mass_driver', 'kWh_kg', self.case), tech.parameter('mass_driver', 'kWh_kg_2241', self.case)
            kwh = a + (b - a) * min(1.0, (self.year - first) / max(1, 2241 - first))
            launcher = tech.parameter('mass_driver', 'launcher_kg_per_kg_year', self.case)
            return dict(how='mass_driver', propellant_kg=0.0, kWh=kwh, launcher_kg_per_kg_year=launcher)
        ve = v(tr['chemical_isp_s']) * G0
        r = math.exp(v(tr['moon_to_depot_m_s']) / ve)
        return dict(how='rocket', propellant_kg=(1 + v(tr['moon_vehicle_kg_per_kg'])) * (r - 1), kWh=0.0, launcher_kg_per_kg_year=0.0)


@lru_cache(maxsize=None)
def reach(percentile):
    """The source survey's contained water within each depot leg, nearest first: the leg's upper edges in m/s and the
    water within each, in kilograms."""
    r = json.load(gzip.open(REACH))
    return r['depot_leg_upper_edges_m_s'], r['contained_water_within_kg'][percentile]


def leg_for(drawn_kg, percentile, set_aside=0.0):
    """The leg the tugs fly once the cumulative draw of contained water has reached `drawn_kg`, mining the nearest
    water first, with a share `set_aside` of each body's water left unmined; None once the water open to mining is
    exhausted. Within each band of the survey the leg lengthens in step with the water drawn from the band, from its
    lower edge to its upper, since the nearest water is mined first."""
    edges, within = reach(percentile)
    lower, below = max(0.0, 2 * edges[0] - edges[1]), 0.0
    for edge, w in zip(edges, within):
        usable = w * (1.0 - set_aside)
        if usable >= drawn_kg and usable > 0:
            if usable <= below:
                return edge
            return lower + (edge - lower) * max(0.0, drawn_kg - below) / (usable - below)
        lower, below = edge, usable
    return None


def usable_water_kg(percentile, set_aside=0.0):
    """All the survey's contained water open to mining at a percentile, after the share set aside."""
    return reach(percentile)[1][-1] * (1.0 - set_aside)


def tug(kind, leg_m_s, prices, limit=True):
    """One kind of tug on a leg: per kilogram of water delivered net to the depot, the water drawn at the source, the
    tug's dry mass per kilogram a year, its power, and whether it closes. Kinds: chemical (hydrogen and oxygen from the
    water), electric (solar arrays at the year's specific power, water propellant), nuclear_electric (once arrived).
    The practical limit on a tug's mass applies to building new tugs (`limit`); tugs already built fly on beyond it,
    delivering less, for as long as the leg closes at all."""
    tr = inputs()['transport']
    t = prices.t
    s, k = v(tr['tug_structure_kg_per_kg_payload']), v(tr['tug_tank_kg_per_kg_propellant'])
    trip = v(tr['tug_round_trip_years'])
    if kind == 'chemical':
        ve = v(tr['chemical_isp_s']) * G0
        r = math.exp(leg_m_s / ve)
        if k * (r - 1) >= 1:
            return None
        dry = (s + k * (r - 1)) / (1 - k * (r - 1))
        f_in, f_out = (dry + 1) * (r - 1), dry * (r - 1)
        power_W = 0.0
    elif kind in ('electric', 'nuclear_electric'):
        if kind == 'nuclear_electric' and not t.arrived('nuclear_electric'):
            return None
        isp = v(tr['electric_isp_s']) if kind == 'electric' else tech.parameter('nuclear_electric', 'isp_s', prices.case)
        ve = isp * G0
        r = math.exp(leg_m_s * v(tr['electric_leg_multiple']) / ve)
        powered_s = v(tr['electric_powered_years']) * 365.25 * 86400
        eff = v(tr['electric_efficiency'])
        kg_per_W = 1 / t.array_W_kg if kind == 'electric' else t.fission_kg_kWe / 1000.0
        dry = s
        for _ in range(60):
            f_in, f_out = (dry + 1) * (r - 1), dry * (r - 1)
            power_W = 0.5 * (f_in + f_out) / powered_s * ve ** 2 / eff
            new = s + k * f_in + power_W * kg_per_W
            if abs(new - dry) < 1e-9:
                break
            dry = new
            if dry > 1e3:
                return None
    else:
        raise ValueError(kind)
    if f_out >= 1:
        return None
    net = 1 - f_out
    # Near the end of a tug's reach its mass per kilogram delivered grows without bound; past the practical limit the
    # leg is out of that tug's reach, since the factories could not build such a fleet as fast as the water is wanted.
    if limit and dry * trip / net > v(tr['tug_max_dry_kg_per_kg_year']):
        return None
    return dict(kind=kind, drawn_kg=(1 + f_in) / net, dry_kg_per_kg_year=dry * trip / net, power_W_per_kg_year=power_W * trip / net,
                propellant_kg=(f_in + f_out) / net)
