"""The distribution of the Open Moon's flyers, from canopies and micro gliders to the largest ships: what flies at
each size under a sixth of Earth's gravity, how many of each the summit metropolis and its port need, how they share
the air by height, and the long-haul ship for the port's crown.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.sky_fleet.run   # results/sky_fleet.json

It reads the sky-ship study's product (the rigid-hull model fitted to four airships' weight statements, the flight
band's winds, the port's long-haul traffic), the summit tower study's products (the port's traffic, its canopies and
its site; the sized form's terminal and disks), the design run's air over the whole Moon (climate/gcm/global_winds.py)
and the cloud-resolving equatorial ring's daytime mixed layer and storms (climate/crm), and sizes each class with the
engineering flight models (engineering/flight). Closed-form sizing and a few small fixed-point loops: about a second on one core.
"""
from __future__ import annotations
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

from engineering.flight import buoyant as bu
from engineering.flight import winged as wi
from research.studies.sky_ships import run as ships
from shared.constants import MOON_RADIUS, SYNODIC_MONTH_DAYS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.sky-fleet/1'
SKY_SHIPS = ROOT / 'research' / 'studies' / 'sky_ships' / 'results' / 'sky_ships.json'
PORT = ships.PORT
FORM = ROOT / 'research' / 'studies' / 'summit_tower' / 'results' / 'summit_tower_form.json'
GLOBAL_WINDS = ships.GLOBAL_WINDS
RING = ships.RING

G_MOON, G_EARTH = ships.G_MOON, ships.G_EARTH
HALF_CIRCUMFERENCE_KM = np.pi * MOON_RADIUS / 1e3
MEAN_TRIP_KM = 0.5 * HALF_CIRCUMFERENCE_KM          # the mean great-circle distance between random points
EVIDENCE = ('First-order sizing of every class of flyer with the engineering flight models (a parabolic drag polar '
            'for gliders and muscle-powered wings, momentum theory and a wing for electric flyers that take off '
            'vertically, the rigid-hull model fitted to four airships\' weight statements for buoyant ships), in the '
            'design run\'s air, the flight band\'s winds and the cloud-resolving ring\'s daytime mixed layer. The '
            'reference flyers, people\'s power and the cities\' travel are published values (sources.json). The '
            'metropolis\'s travel is a planning case built on Earth\'s transit cities, not a model of the metropolis; '
            'no flyer has been designed.')
READING_RULE = ('Masses in kg or t as labelled, speeds in m/s, energy per passenger-km in Wh. Heights are above the '
                'model\'s sea level unless labelled "above ground" or "above the summit". The busy hour is the busiest '
                'hour of a working day. Ranges are the planning case\'s low and high values.')


def digest(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def load() -> dict:
    products = dict(ships=json.loads(SKY_SHIPS.read_text()), port=json.loads(PORT.read_text()),
                    form=json.loads(FORM.read_text()),
                    winds=json.loads(GLOBAL_WINDS.read_text()), ring=json.loads(RING.read_text()))
    assert products['ships']['schema'] == 'terluna.research.sky-ships/1'
    assert products['port']['schema'] == 'terluna.research.summit-tower/1'
    return products


def air(winds: dict, z_km: float) -> bu.Air:
    return ships.moon_air(winds, z_km)


# ---------------------------------------------------------------------------------------------------------------
# Where the flyers fly

CITY_DROP_KM = 4.0                      # the metropolis runs from the summit to about 4 km lower
CITY_AIR_KM = 1.0                       # the city's flyers: about a kilometre above the summit's ground
FERRY_FULL_KM = 3.0                     # sky ferries' cells full at their working heights' top, above the summit
REGIONAL_KM = (10.0, 15.0)              # the port's regional docks, above the summit
FAST_CRUISE_KM = 55.0                   # winged liners cruise above the band, in air like Earth's at 7.3 km
PLATFORM_KM = 70.0


# ---------------------------------------------------------------------------------------------------------------
# People on wings: canopies, micro gliders and muscle-powered wings

PERSON_KG = 80.0                          # a person with clothes, as in the summit tower study's canopies
LEAST_SINK_PER_GLIDE_SPEED = 2.0 / 3.0 ** 0.75    # least sink = 0.877 x (speed of best glide) / (best glide)


def fit_polar(span_m, area_m2, mass_kg, glide, glide_speed_m_s=None, least_sink_m_s=None, max_lift=1.4, density=1.225,
              gravity=G_EARTH):
    """The parasite drag area and span efficiency that give a glider its published best glide at its speed of best
    glide: L/D = 0.5 sqrt(pi e b^2 / f) fixes e / f, and V^4 = (m g)^2 / ((rho/2)^2 pi e b^2 f) fixes e f. Where only
    the least sink is published, the parabolic polar's least sink = 0.877 V / (L/D) gives the speed."""
    if glide_speed_m_s is None:
        glide_speed_m_s = least_sink_m_s * glide / LEAST_SINK_PER_GLIDE_SPEED
    ratio = 4.0 * glide ** 2 / (np.pi * span_m ** 2)
    product = (mass_kg * gravity) ** 2 / ((0.5 * density) ** 2 * np.pi * span_m ** 2 * glide_speed_m_s ** 4)
    return wi.Polar(span_m, area_m2, float(np.sqrt(product / ratio)), float(np.sqrt(product * ratio)), max_lift)


def glider_row(name: str, polar: wi.Polar, mass_kg: float, gravity: float, density: float, efficiency=None,
               power_w=None) -> dict:
    """A glider's (or wing's) least sink, best glide and their speeds; with a power and an efficiency, the level
    speed and the climb it holds on that power."""
    low = wi.min_power(polar, mass_kg, gravity, density, stall_margin=1.2)
    flat = wi.best_glide(polar, mass_kg, gravity, density, stall_margin=1.2)
    row = dict(name=name, span_m=round(polar.span_m, 1), area_m2=round(polar.area_m2, 1), mass_kg=round(mass_kg, 0),
               stall_m_s=round(float(polar.stall_speed_m_s(mass_kg * gravity, density)), 1),
               least_sink_m_s=round(float(low['sink_m_s']), 2), least_sink_speed_m_s=round(float(low['speed_m_s']), 1),
               best_glide=round(float(flat['lift_to_drag']), 1), best_glide_speed_m_s=round(float(flat['speed_m_s']), 1),
               least_power_w=round(float(low['power_w']), 0))
    if efficiency is not None:
        row['least_shaft_power_w'] = round(float(low['power_w']) / efficiency, 0)
        if power_w is not None:
            row.update(level_flight_on_power(polar, mass_kg, gravity, density, efficiency, power_w))
    return row


def level_flight_on_power(polar: wi.Polar, mass_kg, gravity, density, efficiency, power_w) -> dict:
    """At a person's sustained power: the fastest level speed it holds, and the climb rate at the speed of least
    power."""
    out = {}
    for label, p in power_w.items():
        v = np.linspace(float(wi.min_power(polar, mass_kg, gravity, density)['speed_m_s']), 60.0, 6000)
        need = wi.flight_point(polar, mass_kg, gravity, density, v, efficiency)['power_w']
        ok = need <= p
        low = wi.min_power(polar, mass_kg, gravity, density, efficiency)
        out[label] = dict(power_w=p, top_level_speed_m_s=round(float(v[ok].max()), 1) if ok.any() else None,
                          climb_m_s=round(max(0.0, (p - float(low['power_w'])) * efficiency / (mass_kg * gravity)), 2))
    return out


def soaring(ring: dict) -> dict:
    """The cloud-resolving ring's daytime mixed layer over land, by hour angle: how deep the thermals rise and for how
    much of the lunar day, before and after the afternoon storms."""
    b = ring['by_hour_angle']
    ha = np.asarray(b['hour_angle_deg'], dtype=float)
    land = b['land']
    mixed, rain = np.asarray(land['mixed_layer_km'], dtype=float), np.asarray(land['rain_mm_h'], dtype=float)
    base = np.asarray(land['cloud_base_lcl_km'], dtype=float)
    hours_per_degree = SYNODIC_MONTH_DAYS * 24.0 / 360.0
    step = float(np.median(np.diff(ha))) * hours_per_degree
    deep = mixed >= 2.0
    dry = rain < 0.01
    k = int(np.argmax(mixed))
    return dict(mixed_layer_max_km=round(float(mixed[k]), 1), at_hour_angle_deg=float(ha[k]),
                cloud_base_at_max_km=round(float(base[k]), 1),
                hours_mixed_layer_over_2_km=round(float(deep.sum() * step), 0),
                hours_dry_and_over_2_km=round(float((deep & dry).sum() * step), 0),
                night_mixed_layer_km=round(float(np.median(mixed[np.abs(ha) > 120.0])), 1),
                storm_hours_rain_over_0_1_mm_h=round(float((rain > 0.1).sum() * step), 0),
                land_wet_share=round(ring['rain_rate_mm_h']['land']['wet_share'], 3),
                rain_spells_per_lunar_day=round(ring['rain_spells']['land']['spells_per_lunar_day'], 2),
                rain_spell_hours=ring['rain_spells']['land']['duration_h'],
                hour_angle_deg=ha.tolist(), mixed_layer_km=mixed.tolist(), rain_mm_h=rain.round(3).tolist())


# ---------------------------------------------------------------------------------------------------------------
# Electric flyers that take off vertically: personal flyers, sky boats, air taxis, cargo drones

def calibrate_electric(evtols: dict, tech: wi.Electric) -> tuple:
    """The airframe share that gives each of Earth's electric flyers its published take-off mass for its payload,
    range, cruise speed and rotor disc; their mean becomes the technology's share."""
    shares = {}
    for name, e in evtols.items():
        lo, hi = 0.05, 0.9
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            t = wi.Electric(**{**tech.__dict__, 'airframe_share': mid,
                                'disc_loading_kg_m2': e['mtow_kg'] / e['disc_m2']})
            r = wi.electric_vtol(e['payload_kg'], G_EARTH, 1.225, e['cruise_m_s'], e['range_m'], e['hover_s'],
                                 e['climb_m'], t)
            if r is None or r['mass_kg'] > e['mtow_kg']:
                hi = mid
            else:
                lo = mid
        shares[name] = round(0.5 * (lo + hi), 3)
    return float(np.mean(list(shares.values()))), shares


def electric_row(name: str, spec: dict, tech: wi.Electric, gravity: float, density: float, city_trip_km: float) -> dict:
    r = wi.electric_vtol(spec['payload_kg'], gravity, density, spec['cruise_m_s'], spec['range_km'] * 1e3,
                         spec['hover_s'], spec['climb_m'], tech)
    if r is None:
        return dict(name=name, closes=False)
    trip = wi.vtol_trip(r['mass_kg'], gravity, density, r['disc_area_m2'], tech.lift_to_drag, spec['cruise_m_s'],
                        city_trip_km * 1e3, spec['hover_s'], spec['climb_m'], tech.figure_of_merit, tech.efficiency)
    per_km = float(trip['energy_j']) / city_trip_km / 3600.0
    seats = spec.get('seats')
    return dict(name=name, closes=True, payload_kg=spec['payload_kg'], seats=seats, mass_kg=round(r['mass_kg'], 0),
                battery_kg=round(float(r['battery_kg']), 0), drive_kg=round(float(r['drive_kg']), 0),
                disc_area_m2=round(float(r['disc_area_m2']), 1),
                rotors_across_m=round(float(np.sqrt(4.0 * r['disc_area_m2'] / (np.pi * spec['rotors']))), 2),
                hover_kw=round(float(r['hover_w']) / 1e3, 1), cruise_kw=round(float(r['cruise_w']) / 1e3, 1),
                cruise_m_s=spec['cruise_m_s'], range_km=spec['range_km'],
                city_trip_km=city_trip_km, city_trip_kwh=round(float(trip['energy_j']) / 3.6e6, 2),
                city_trip_wh_per_km=round(per_km, 1),
                wh_per_seat_km=round(per_km / seats, 1) if seats else None)


# Small buoyant ships on Earth keep this share of what they lift for payload (sources.json): the Zeppelin NT 1,900 of
# 8,050 kg, the Skyship 600 a disposable load of 2,590 kg on a maximum of 6,630 kg. Non-rigid ships of fineness
# 3.5-4.4 had drag coefficients of 0.035-0.049 on volume to the two-thirds in NACA's full-scale tests (Report 397).
SMALL_SHIP_PAYLOAD_SHARE = (1900.0 / 8050.0, 2590.0 / 6630.0)
SMALL_SHIP_DRAG = 0.04


def buoyant_boat(payload_kg: float, fly: bu.Air, speeds_m_s=(15.0, 20.0, 30.0), fineness=4.0) -> dict:
    """A sky boat held up by helium instead of rotors: the hull that carries the payload at the payload shares of
    Earth's small ships, and the energy per kilometre its drag takes at a few speeds."""
    lift = float(bu.net_lift_kg_m3(fly, bu.HELIUM, ships.PURITY))
    rows = []
    for share in SMALL_SHIP_PAYLOAD_SHARE:
        hull = bu.Hull.of_volume(payload_kg / share / lift, fineness)
        drag = {f'{v:.0f}': 0.5 * fly.density_kg_m3 * v ** 2 * SMALL_SHIP_DRAG * float(hull.reference_area_m2)
                for v in speeds_m_s}
        rows.append(dict(payload_share=round(share, 2), volume_m3=round(float(hull.volume_m3), -1),
                         length_m=round(float(hull.length_m), 1), diameter_m=round(float(hull.diameter_m), 1),
                         wh_per_km={k: round(d / 0.7 / 3.6, 0) for k, d in drag.items()}))
    return dict(payload_kg=payload_kg, lift_kg_m3=round(lift, 3), drag_coefficient=SMALL_SHIP_DRAG, hulls=rows)


# ---------------------------------------------------------------------------------------------------------------
# Buoyant ships: ferries, regional ships, long-haul liners, freighters

TOP_OVER_CRUISE = ships.SPEED_M_S / ships.CRUISE_M_S      # top speed, for strength and engines, over cruise (1.2)
PER_JOULE = 1.0 / (ships.FUEL_CELL_EFFICIENCY * ships.HYDROGEN_HEAT_J_KG * ships.TANK_SHARE)   # kg of fuel and tanks


def buoyant_design(factors: dict, cruise_m_s: float, gust_m_s: float) -> bu.Design:
    """The sky-ship study's modern rigid hull, on its fitted factors, built for a top speed 1.2 times the cruise."""
    return bu.Design(**{**ships.DESIGN_1930, **factors, **ships.MODERN, 'speed_m_s': TOP_OVER_CRUISE * cruise_m_s,
                        'gust_m_s': gust_m_s})


def ship(fly: bu.Air, full: bu.Air, factors: dict, length_m: float, cruise_m_s: float, route_km: float,
         passenger_kg: float, gust_m_s: float, hotel_w: float = ships.HOTEL_W_PER_PASSENGER, fineness=6.0) -> dict:
    """A rigid hydrogen ship flying in `fly`, its cells full in `full` (the top of its working heights), carrying
    people over a route with fuel cells: the payload after the route's fuel and the hotel load, and what it takes.
    At 30 m/s over the sky-ship study's 48-hour route it is that study's long-haul ship."""
    lift = float(bu.net_lift_kg_m3(full, bu.HYDROGEN, ships.PURITY))
    design = buoyant_design(factors, cruise_m_s, gust_m_s)
    hull = bu.Hull(float(length_m), fineness)
    r = bu.rigid(hull, fly, lift, design)
    cruise_w = float(bu.power_w(hull, fly, cruise_m_s, design))
    seconds = route_km * 1e3 / cruise_m_s
    fuel = cruise_w * seconds * PER_JOULE
    hotel = hotel_w * seconds * PER_JOULE
    n = max(0.0, (float(r['useful_kg']) - fuel) / (passenger_kg + hotel))
    fins = design.fin_area_share * hull.reference_area_m2 / 4.0
    h2 = float(hull.volume_m3) * ships.PURITY * float(bu.gas_density(full, bu.HYDROGEN, 1.0))
    energy = (cruise_w + hotel_w * n) / (n * cruise_m_s) / 3.6 if n > 0 else None
    return dict(length_m=float(length_m), diameter_m=round(float(hull.diameter_m), 1),
                fin_span_m=round(float(np.sqrt(design.fin_aspect_ratio * fins)), 1),
                volume_m3=float(f'{float(hull.volume_m3):.3g}'), gross_lift_t=round(float(r['gross_lift_kg']) / 1e3, 1),
                useful_share=round(float(r['useful_share']), 3), route_km=round(route_km, 0), cruise_m_s=cruise_m_s,
                fuel_t=round((fuel + hotel * n) / 1e3, 1), passengers=round(n, -1) if n >= 100 else round(n, 0),
                payload_t=round((float(r['useful_kg']) - fuel - hotel * n) / 1e3, 1),
                cruise_mw=round(cruise_w / 1e6, 2), wh_per_passenger_km=None if energy is None else round(energy, 1),
                hydrogen_t=round(h2 / 1e3, 1))


def mooring(fly: bu.Air, row: dict, factors: dict, design_wind: float, service_wind: float) -> dict:
    """Wind loads on a berthed ship: nose-on (weathervaned) in the design wind, side-on in the service wind."""
    hull = bu.Hull(row['length_m'], 6.0)
    design = buoyant_design(factors, row['cruise_m_s'], 0.0)
    q = lambda v: 0.5 * fly.density_kg_m3 * v ** 2
    side = 0.8 * hull.length_m * hull.diameter_m
    return dict(nose_on_design_mn=round(float(q(design_wind) * design.drag_coefficient * hull.reference_area_m2) / 1e6, 2),
                side_on_service_mn=round(float(q(service_wind) * ships.BEAM_DRAG * side) / 1e6, 2))


def network(passengers_per_hour, per_ship: float, trip_km: float, cruise_m_s: float, berth_h: float) -> dict:
    """The service a class of ship gives at a port: calls an hour (each unloads and loads a full ship), the
    destinations it can serve at a given frequency, the berths busy at once and the fleet cycling through it."""
    calls = [p / (2.0 * per_ship) for p in passengers_per_hour]
    trip_h = trip_km * 1e3 / cruise_m_s / 3600.0
    return dict(calls_per_hour=[round(c, 2) for c in calls],
                destinations_hourly=[round(c, 1) for c in calls],
                destinations_every_4_hours=[round(4.0 * c, 0) for c in calls],
                destinations_daily=[round(24.0 * c, 0) for c in calls],
                berths=[round(c * berth_h, 1) for c in calls],
                mean_trip_hours=round(trip_h, 1),
                fleet=[round(c * 2.0 * (trip_h + berth_h), -1) for c in calls])


CLEARANCE_M = 40.0                   # between the fins of ships stacked one above another
SWING_DEG = 40.0                     # how far berthed ships may swing together from the prevailing wind
ARM_BEYOND = 0.3                     # the arm runs this share of a ship's width past its last berth


def crown_berths(berths: int, row: dict, crown_width_m: float, lowest_m: float, top_m: float) -> dict:
    """Berths at the crown for one class: ships lie nose-in on the lee side of arms running north and south from the
    frame, noses into the steady west wind. Side by side along an arm they are spaced so that, swinging together
    with the wind, their fins meet only past SWING_DEG. Levels are stacked a hull, its fins and a clearance apart,
    from the lowest height a ship clears the terminal's roof up to the frame's top."""
    across = row['diameter_m'] + 2.0 * row['fin_span_m']
    pitch = across / np.cos(np.radians(SWING_DEG))
    level = across + CLEARANCE_M
    first = lowest_m + 0.5 * across + CLEARANCE_M
    fit = int(max(0.0, top_m - first) // level) + 1
    levels = int(max(1, min(fit, int(np.ceil(berths / 2.0)))))
    arms = 2 * levels if berths > 1 else 1
    per_arm = int(np.ceil(berths / arms))
    arm = pitch * (per_arm - 0.5) + ARM_BEYOND * across     # to the last berth, and a little past it
    return dict(berths=berths, pitch_m=round(pitch, 0), level_spacing_m=round(level, 0), levels=levels, arms=arms,
                ships_per_arm=per_arm, arm_length_m=round(arm, 0),
                levels_above_summit_m=[round(first + k * level, 0) for k in range(levels)],
                reach_from_axis_m=round(crown_width_m / 2.0 + arm, 0), footprint_east_m=round(row['length_m'], 0))


# ---------------------------------------------------------------------------------------------------------------
# Winged liners, freighters and high platforms

def winged_liner(fly: bu.Air, dock: bu.Air, seats: int, kg_per_seat: float, wing_loading_pa: float,
                 lift_coefficient: float, lift_to_drag: float, disc_loading_pa: float, efficiency=0.8,
                 figure_of_merit=0.7) -> dict:
    """A winged liner that takes off and lands vertically on rotors that tilt forward to cruise: its cruise speed at
    a wing loading and lift coefficient, its energy per passenger-kilometre with the airships' hotel load, its hover
    power at the dock and the wing it needs."""
    mass = seats * kg_per_seat
    weight = mass * fly.gravity_m_s2
    speed = float(wi.speed_for_lift(wing_loading_pa / fly.gravity_m_s2, fly.gravity_m_s2, fly.density_kg_m3,
                                    lift_coefficient))
    area = weight / wing_loading_pa
    hover = float(wi.hover_power_w(mass, dock.gravity_m_s2, dock.density_kg_m3, weight / disc_loading_pa,
                                   figure_of_merit))
    per_pkm = (kg_per_seat * fly.gravity_m_s2 / (lift_to_drag * efficiency) + ships.HOTEL_W_PER_PASSENGER / speed) / 3.6
    return dict(seats=seats, mass_t=round(mass / 1e3, 0), wing_area_m2=round(area, 0),
                span_at_aspect_ratio_10_m=round(float(np.sqrt(10.0 * area)), 0), cruise_m_s=round(speed, 0),
                wh_per_passenger_km=round(per_pkm, 1), hover_mw=round(hover / 1e6, 1),
                mean_trip_hours=round(MEAN_TRIP_KM * 1e3 / speed / 3600.0, 1),
                longest_trip_hours=round(HALF_CIRCUMFERENCE_KM * 1e3 / speed / 3600.0, 1))


def freighter(fly: bu.Air, full: bu.Air, factors: dict, length_m: float, cruise_m_s: float, route_km: float,
              gust_m_s: float) -> dict:
    """A rigid hydrogen freighter: its cargo after the route's fuel, and the energy per tonne-kilometre."""
    r = ship(fly, full, factors, length_m, cruise_m_s, route_km, 1000.0, gust_m_s, hotel_w=0.0)
    cargo = r['payload_t']
    return dict(length_m=length_m, cruise_m_s=cruise_m_s, gross_lift_t=r['gross_lift_t'], cargo_t=cargo,
                wh_per_tonne_km=round(r['cruise_mw'] * 1e6 / (cargo * cruise_m_s) / 3.6, 1) if cargo > 0 else None,
                hydrogen_t=r['hydrogen_t'])


def platform(high: bu.Air, winds: dict, z_km: float, length_m: float, fabric: bu.Fabric, payload_share=0.3) -> dict:
    """A pressure hull that holds station high above the flight band: its lift, envelope and the power to hold
    against the wind at its height (the 99th percentile and the highest of the design run's 3-day means)."""
    lift = float(bu.net_lift_kg_m3(high, bu.HYDROGEN, ships.PURITY))
    p99, top = ships.at(winds, 'speed_p99_m_s', z_km), ships.at(winds, 'speed_max_m_s', z_km)
    hull = bu.Hull(length_m, 4.0)
    design = bu.Design(speed_m_s=top * ships.GUST_FACTOR, gust_m_s=top * ships.GUST_FACTOR)
    env = bu.pressure_hull(hull, high, lift, fabric, design)
    hold = lambda v: float(bu.power_w(hull, high, v, design))
    return dict(height_km=z_km, length_m=length_m, volume_m3=float(f'{float(hull.volume_m3):.3g}'),
                gross_lift_t=round(float(env['gross_lift_kg']) / 1e3, 1),
                envelope_t=round(float(env['envelope_kg']) / 1e3, 1), fabric_use=round(float(env['fabric_use']), 2),
                payload_t_at_share=round(payload_share * float(env['gross_lift_kg']) / 1e3, 1),
                wind_p99_m_s=round(p99, 1), wind_highest_m_s=round(top, 1),
                hold_p99_kw=round(hold(p99) / 1e3, 0), hold_highest_kw=round(hold(top) / 1e3, 0))


# ---------------------------------------------------------------------------------------------------------------
# The metropolis: trips, fleets, and the sky at the busy hour

CLIMB_M_S = 5.0                          # electric flyers' climb and descent to and from their lanes


def flights(modes: dict, electric: dict, tech: wi.Electric, wing_speed_m_s: float, ferry: dict, gravity: float,
            density: float) -> dict:
    """Each flying mode's minutes aloft and energy per flight: electric flyers hover, climb to their lane, cruise
    the mode's trip and descend; personal wings fly the trip at their cruising speed on their own power; ferries
    carry a passenger at their energy per passenger-kilometre at 70% load."""
    out = {}
    for name, m in modes.items():
        if 'flyer' not in m:
            continue
        if m['flyer'] in electric:
            spec, row = ELECTRIC_CLASSES[m['flyer']], electric[m['flyer']]
            trip = wi.vtol_trip(row['mass_kg'], gravity, density, row['disc_area_m2'], tech.lift_to_drag,
                                spec['cruise_m_s'], m['km'] * 1e3, spec['hover_s'], spec['climb_m'],
                                tech.figure_of_merit, tech.efficiency)
            seconds = spec['hover_s'] + 2.0 * spec['climb_m'] / CLIMB_M_S + m['km'] * 1e3 / spec['cruise_m_s']
            out[name] = dict(minutes=seconds / 60.0, kwh=float(trip['energy_j']) / 3.6e6)
        elif m['flyer'] == 'personal wing':
            out[name] = dict(minutes=m['km'] * 1e3 / wing_speed_m_s / 60.0 + 2.0, kwh=0.0)
        elif m['flyer'] == 'ferry':
            out[name] = dict(minutes=m['km'] / ferry['commercial_speed_km_h'] * 60.0,
                             kwh=m['km'] * ferry['wh_per_passenger_km'] / 0.7 / 1e3)
    return out


def metropolis(city: dict, modes: dict, per_flight: dict, ferry: dict) -> dict:
    """The planning case's trips by mode at low, base and high travel, and for each flying mode its flights, the
    flyers aloft at the busy hour, its fleet and its energy."""
    out = {}
    for i, case in enumerate(('low', 'base', 'high')):
        daily = city['population'] * city['trips_per_person'][i]
        busy = city['busy_hour_share'][i]
        rows = {}
        for name, m in modes.items():
            t = daily * m['share']
            row = dict(trips_per_day=t, busy_hour_trips=t * busy)
            if name in per_flight:
                f = per_flight[name]
                if m['flyer'] == 'ferry':
                    seats = ferry['passengers'] * 0.7
                    per_day = t * m['km'] / (seats * FERRY_LINE_KM)
                    aloft = t * busy * f['minutes'] / 60.0 / seats
                    row.update(ferry_runs_per_day=per_day, aloft_busy_hour=aloft,
                               fleet_in_service=aloft * (1.0 + FERRY_SPARES))
                else:
                    per_day = t / m['occupancy']
                    aloft = per_day * busy * f['minutes'] / 60.0
                    row.update(flights_per_day=per_day, aloft_busy_hour=aloft,
                               fleet_private=city['population'] * m['owned_per_person'] if m.get('owned_per_person') else None,
                               fleet_shared=per_day / m['flights_per_vehicle_day'] if m.get('flights_per_vehicle_day') else None)
                if m.get('berth_m2'):
                    for key in ('fleet_private', 'fleet_shared'):
                        if row.get(key):
                            row[key.replace('fleet', 'berths_km2')] = row[key] * m['berth_m2'] / 1e6
                row.update(per_km2_busy_hour=aloft / city['area_km2'], minutes_aloft=f['minutes'],
                           energy_gwh_day=t * f['kwh'] / 1e6 if m['flyer'] == 'ferry' else per_day * f['kwh'] / 1e6)
            rows[name] = {k: (round(v, -2) if isinstance(v, float) and v > 1e4 else
                              round(v, 2) if isinstance(v, float) else v) for k, v in row.items()}
        out[case] = rows
    return out


def airspace(aloft: dict, bands: dict, area_km2: float) -> dict:
    """Flyers aloft in each height band at the busy hour: per km2, per km3, and the mean distance to the nearest
    neighbour in a layer of the band's own thickness per layer."""
    out = {}
    for band, b in bands.items():
        n = sum(aloft.get(m, 0.0) * s for m, s in b['modes'].items())
        depth = b['top_m'] - b['bottom_m']
        per_layer = n / area_km2 / b['layers']
        out[band] = dict(bottom_m=b['bottom_m'], top_m=b['top_m'], layers=b['layers'], aloft=round(n, -2),
                         per_km2=round(n / area_km2, 1), per_km3=round(n / area_km2 / (depth / 1e3), 1),
                         nearest_in_layer_m=round(500.0 / np.sqrt(per_layer), 0) if per_layer > 0 else None,
                         speed_m_s=b['speed_m_s'])
    return out


# ---------------------------------------------------------------------------------------------------------------
# The cases

# Earth's reference gliders (sources.json): span and wing area (projected for the paraglider), all-up mass, best glide
# and its speed or the least sink; `fit_polar` reproduces them, and the model's other figures are set beside the
# published ones. The paraglider's projected span is its flat span times Ozone's projected-to-flat ratio (0.79).
EARTH_GLIDERS = {
    'paraglider, EN-B (MAC PARA Eden 7, size 22)': dict(span_m=0.79 * 11.53, area_m2=19.26, mass_kg=67.5, glide=10.0,
                                                        glide_speed_m_s=38.0 / 3.6, max_lift=1.0,
                                                        published=dict(least_sink_m_s=1.05)),
    'sailplane, 18 m (HpH 304S Shark), 600 kg': dict(span_m=18.0, area_m2=11.74, mass_kg=600.0, glide=50.0,
                                                    glide_speed_m_s=125.0 / 3.6, max_lift=1.5,
                                                    published=dict(least_sink_m_s=0.47, at_mass_kg='minimum')),
    'human-powered aircraft (Musculair 2)': dict(span_m=19.5, area_m2=11.7, mass_kg=78.0, glide=37.0,
                                                 least_sink_m_s=0.27, max_lift=1.4,
                                                 published=dict(power_w=250.0, at_speed_m_s=10.0)),
    'human-powered aircraft (MIT Light Eagle)': dict(span_m=34.7, area_m2=31.0, mass_kg=110.0, glide=36.0,
                                                     glide_speed_m_s=7.8, max_lift=1.4,
                                                     published=dict(power_w=200.0, at_speed_m_s=7.8)),
}
PROPELLER_EFFICIENCY = 0.85            # propeller and drive of a human-powered aircraft
# People's sustained mechanical power at the pedals (sources.json): sedentary and recreational adults rode three hours
# at 87% of their first lactate threshold, about 70 W for the women and 100-105 W for the men (Matomaki et al. 2023);
# untrained women's critical power is 1.8 W/kg (about 113 W; Kavaliauskas et al. 2017), untrained men's 155 W (Lei et
# al. 2026); NASA's curve gives young, fit men 215 W for an hour (Webb 1973).
HUMAN_POWER_W = {'an untrained adult for three hours': 70.0, "an untrained woman's critical power": 113.0,
                 "an untrained man's critical power": 155.0, 'a fit young man for an hour': 215.0}
# Parachutes (Knacke 1991): a flat circular canopy's drag coefficient is 0.75-0.80 on its cloth area, and it inflates
# to 0.67-0.70 of its flat diameter.
ROUND_CANOPY_CD0 = 0.78
ROUND_CANOPY_INFLATED = 0.685
# Earth's electric flyer that takes off vertically, for the airframe share (sources.json): the Joby S4's certified
# take-off mass (4,800 lb) and payload (up to 1,000 lb, a pilot and four passengers), its target range (100 miles),
# a cruise of 150 mph under its 200 mph top speed, and rotors at the 2015 design's disc loading (six 9.5 ft rotors on
# 4,000 lb, 46 kg/m2). The design's published hover estimate (560 hp) is the momentum-theory power over a figure of
# merit of 0.58; its trip here hovers 90 s and climbs 450 m.
EVTOLS = {'Joby S4 (JAS4-1)': dict(mtow_kg=4800 * ships.LB, payload_kg=1000 * ships.LB,
                                   disc_m2=4800 * ships.LB / (4000 * ships.LB / (6 * np.pi * (9.5 * ships.FT / 2) ** 2)),
                                   cruise_m_s=150 * 0.44704, range_m=100 * 1609.344, hover_s=90.0, climb_m=450.0)}
# Technology: Joby's 235 Wh/kg pack; a whole-aircraft figure of merit of 0.6 (the 2015 estimate's 0.58); a
# lift-to-drag ratio of 14, between the 11 of Archer's designs and the 16-18 assumed for Beta's and Lilium's; motors,
# controllers and rotors at 2 kW per kg of hover power (motors alone give 3.6-13 kW/kg continuous).
TECH = wi.Electric(airframe_share=0.4, lift_to_drag=14.0, figure_of_merit=0.6, efficiency=0.8, drive_w_kg=2000.0,
                   battery_wh_kg=235.0, usable=0.8, reserve=1.3)

# Lunar designs. The micro glider and the personal wing are proposals at plausible sizes; their drag areas are the
# wing's profile drag (C_D0 on the area) and the pilot's.
MICRO_GLIDER = dict(span_m=9.0, area_m2=4.5, craft_kg=18.0, profile_cd=0.012, pilot_drag_area_m2=0.15,
                    span_efficiency=0.85, max_lift=1.4)
PERSONAL_WING = dict(span_m=11.0, area_m2=7.0, craft_kg=22.0, profile_cd=0.010, pilot_drag_area_m2=0.10,
                     span_efficiency=0.9, max_lift=1.3, efficiency=0.8)
ELECTRIC_CLASSES = {
    'personal flyer (one seat)': dict(payload_kg=100.0, seats=1, cruise_m_s=30.0, range_km=60.0, hover_s=60.0,
                                      climb_m=300.0, rotors=4),
    'sky boat (four seats)': dict(payload_kg=400.0, seats=4, cruise_m_s=40.0, range_km=150.0, hover_s=90.0,
                                  climb_m=500.0, rotors=6),
    'air taxi or sky van (twelve seats, 1.2 t)': dict(payload_kg=1200.0, seats=12, cruise_m_s=45.0, range_km=200.0,
                                                      hover_s=90.0, climb_m=800.0, rotors=8),
    'parcel drone (5 kg)': dict(payload_kg=5.0, seats=None, cruise_m_s=20.0, range_km=20.0, hover_s=60.0,
                                climb_m=150.0, rotors=4),
    'cargo drone (100 kg)': dict(payload_kg=100.0, seats=None, cruise_m_s=25.0, range_km=40.0, hover_s=60.0,
                                 climb_m=300.0, rotors=6),
}
CITY_TRIP_KM = 10.0                       # a trip across part of the metropolis, for the classes' energy
FERRY_LENGTHS_M = (150.0, 200.0, 250.0, 300.0)
FERRY_CRUISE_M_S = 25.0
FERRY_PASSENGER_KG = 150.0                # a person with a bag (90 kg) and a share of the seats, decks and crew
FERRY_ENERGY_HOURS = 2.0                  # fuel between refills
FERRY_HEADWAYS_MIN = (2.0, 3.0, 5.0)
FERRY_STOP_KM, FERRY_DWELL_S, FERRY_MANOEUVRE_S = 4.0, 120.0, 60.0
FERRY_LINE_KM = 30.0                      # a ferry line across the metropolis
FERRY_SPARES = 0.15                       # ferries in service beyond those aloft at the busy hour
CITY_GUST_M_S = 10.7                      # Earth's airship gust; the ring's highest updraft is 7.6 m/s at 10 km
REGIONAL_LENGTHS_M = (245.0, 350.0, 500.0)
REGIONAL_CRUISE_M_S = 35.0
REGIONAL_PASSENGER_KG = 200.0             # seated, for trips of a few hours
REGIONAL_ROUTE_KM = 1500.0
REGIONAL_MEAN_TRIP_KM = 500.0
LONG_HAUL_LENGTHS_M = (500.0, 750.0, 1000.0, 1500.0)
LONG_HAUL_SPEEDS_M_S = (30.0, 40.0, 50.0)
LONG_HAUL_ROUTE_KM = ships.CRUISE_M_S * ships.TRIP_HOURS * 3.6    # the sky-ship study's 48 hours at 30 m/s
# A winged liner that takes off and lands on tilting rotors: 400 kg a seat, as the airships carry a passenger with a
# cabin, so the comparison isolates the physics; a lift-to-drag ratio of 17 and a figure of merit of 0.7.
FAST_LINER = dict(seats=300, kg_per_seat=400.0, wing_loading_pa=3000.0, lift_coefficient=0.5, lift_to_drag=17.0,
                  disc_loading_pa=600.0, figure_of_merit=0.7)
FREIGHT_LENGTHS_M = (500.0, 1000.0, 1500.0)
FREIGHT_CRUISE_M_S = 25.0
PLATFORM_LENGTHS_M = (100.0, 200.0)
BERTH_HOURS = ships.BERTH_HOURS

CITY_W_PER_PERSON = 2000.0              # the metropolis brief's planning figure
CITY = dict(population=100e6, radius_km=32.0, area_km2=np.pi * 32.0 ** 2,
            # Trips per person per day: Tokyo's region (2018) and London's residents (2024/25) make 2.0, the Paris
            # region's residents 3.7 (EGT H2020, walking counted); the base is between them.
            trips_per_person=(2.0, 2.8, 3.7),
            # The busiest hour's share of the day's trips: 8% (the low end of urban roads' K-factor), 10%, and the
            # 13% of Hong Kong's morning peak hour (TCS 2022).
            busy_hour_share=(0.08, 0.10, 0.13))
# The planning case's trips by mode, after the Tokyo 23-ward area (2018 person trip survey: rail 51%, car 8%, bicycle
# 13%, walking 24%), the densest of Earth's transit cities in the sources: sky boats take cars' place at 10%, air
# taxis 2%, personal wings and gliders take part of cycling at 5%, and the backbone's trips split between metro and
# rail and the sky ferries. `km` is a trip's straight-line length (the Paris region's trips average 4.5 km, its
# public transport trips 9 km; a flyer flies straight), `occupancy` the people per flight (cars carry 1.28 in the
# Paris region, 1.40 in Tokyo's), `owned_per_person` the private fleet (Tokyo-to has 222 private cars per 1,000
# people, Paris 260).
MODES = {
    'walking and cycling': dict(share=0.30),
    'metro and rail': dict(share=0.375),
    'sky ferries': dict(share=0.125, km=9.0, flyer='ferry'),
    'sky boats': dict(share=0.10, km=8.0, occupancy=1.3, flyer='sky boat (four seats)',
                      owned_per_person=0.25, flights_per_vehicle_day=20.0, berth_m2=60.0),
    'air taxis': dict(share=0.02, km=10.0, occupancy=3.0, flyer='air taxi or sky van (twelve seats, 1.2 t)',
                      flights_per_vehicle_day=30.0, berth_m2=150.0),
    'personal wings and gliders': dict(share=0.05, km=4.0, occupancy=1.0, flyer='personal wing',
                                       owned_per_person=0.1),
    'other': dict(share=0.03),
}
# Heights above the metropolis's ground where each mode flies, and how many layers of traffic each band holds.
BANDS = {
    'wings, gliders and canopies': dict(bottom_m=30.0, top_m=300.0, layers=3, speed_m_s=7.0,
                                        modes={'personal wings and gliders': 1.0}),
    'sky boats and air taxis': dict(bottom_m=300.0, top_m=1500.0, layers=6, speed_m_s=40.0,
                                    modes={'sky boats': 1.0, 'air taxis': 1.0}),
    'sky ferries': dict(bottom_m=1500.0, top_m=2500.0, layers=2, speed_m_s=FERRY_CRUISE_M_S,
                        modes={'sky ferries': 1.0}),
}
# Spacings for the sky boats' capacity: NASA's urban air mobility routes lie 1,500 ft (457 m) apart, and in its
# unmitigated simulation a horizontal separation of more than 1,800 ft (549 m) with more than 200 ft (61 m)
# vertically kept the collision risk under 10% (Lee et al. 2022); a near mid-air collision is closer than 500 ft
# (152 m; FAA AIM 7-7-3).
SPACINGS = {
    "NASA's urban air mobility analysis": dict(lateral_m=457.0, in_trail_m=549.0, vertical_m=61.0),
    'lanes at the near-miss distance': dict(lateral_m=152.0, in_trail_m=152.0, vertical_m=61.0),
}


def capacity(band: dict, spacing: dict) -> float:
    """Flyers per km2 a band holds as layers of lanes: layers a vertical spacing apart, lanes a lateral spacing apart
    in each, flyers an in-trail spacing apart along each lane."""
    layers = np.floor((band['top_m'] - band['bottom_m']) / spacing['vertical_m'])
    return float(layers * 1e6 / (spacing['lateral_m'] * spacing['in_trail_m']))


def results() -> dict:
    p = load()
    winds, ring, port, sky, form = p['winds'], p['ring'], p['port'], p['ships'], p['form']
    factors = sky['buoyant']['factors']
    band = sky['band']
    summit_km = port['site']['ground_above_sea_m'] / 1e3
    city = air(winds, summit_km + CITY_AIR_KM)
    summit = air(winds, summit_km)
    out = dict(schema=SCHEMA, evidence=EVIDENCE, reading_rule=READING_RULE)
    out['air'] = {name: dict(height_km=round(z, 1), density_kg_m3=round(air(winds, z).density_kg_m3, 3),
                             temperature_c=round(air(winds, z).temperature_k - 273.15, 1))
                  for name, z in dict(city_low=summit_km - CITY_DROP_KM, summit=summit_km,
                                      city_flight=summit_km + CITY_AIR_KM,
                                      regional_docks=summit_km + np.mean(REGIONAL_KM), flight_band=ships.BAND_KM,
                                      fast_cruise=FAST_CRUISE_KM, platforms=PLATFORM_KM).items()}

    # People on wings
    gliders = []
    for name, g in EARTH_GLIDERS.items():
        pol = fit_polar(g['span_m'], g['area_m2'], g['mass_kg'], g['glide'], g.get('glide_speed_m_s'),
                        g.get('least_sink_m_s'), g['max_lift'])
        earth = glider_row(name, pol, g['mass_kg'], G_EARTH, 1.225)
        earth['published'] = g['published']
        if 'power_w' in g['published']:
            earth['model_power_at_published_speed_w'] = round(float(wi.flight_point(
                pol, g['mass_kg'], G_EARTH, 1.225, g['published']['at_speed_m_s'], PROPELLER_EFFICIENCY)['power_w']), 0)
        if g['published'].get('at_mass_kg') == 'minimum':
            light = wi.min_power(pol, 445.0, G_EARTH, 1.225, stall_margin=0.0)
            earth['model_least_sink_at_445_kg_m_s'] = round(float(light['sink_m_s']), 2)
        moon = glider_row(name, pol, g['mass_kg'], G_MOON, summit.density_kg_m3,
                          PROPELLER_EFFICIENCY if 'power_w' in g['published'] else None)
        small = wi.Polar(pol.span_m / np.sqrt(ships.SIMILAR), pol.area_m2 / ships.SIMILAR, pol.drag_area_m2 / ships.SIMILAR,
                         pol.span_efficiency, pol.max_lift_coefficient)
        sixth = glider_row(name + ', a sixth of the wing', small, g['mass_kg'], G_MOON, summit.density_kg_m3)
        gliders.append(dict(name=name, drag_area_m2=round(pol.drag_area_m2, 3), span_efficiency=round(pol.span_efficiency, 2),
                            earth_sea_level=earth, moon_summit=moon, moon_summit_sixth_of_wing=sixth))
    designs = {}
    for name, d in dict(micro_glider=MICRO_GLIDER, personal_wing=PERSONAL_WING).items():
        pol = wi.Polar(d['span_m'], d['area_m2'], d['area_m2'] * d['profile_cd'] + d['pilot_drag_area_m2'],
                       d['span_efficiency'], d['max_lift'])
        mass = PERSON_KG + d['craft_kg']
        designs[name] = dict(moon_summit=glider_row(name, pol, mass, G_MOON, summit.density_kg_m3, d.get('efficiency'),
                                                    HUMAN_POWER_W if d.get('efficiency') else None),
                             earth_sea_level=glider_row(name, pol, mass, G_EARTH, 1.225, d.get('efficiency'),
                                                        HUMAN_POWER_W if d.get('efficiency') else None))
    evac = port['port']['evacuation']
    flat = [float(wi.canopy_diameter_m(PERSON_KG, G_MOON, rho, evac['person']['landing_m_s'], ROUND_CANOPY_CD0))
            for rho in (evac['canopy']['base']['density_kg_m3'], evac['canopy']['top']['density_kg_m3'])]
    earth_flat = float(wi.canopy_diameter_m(PERSON_KG, G_EARTH, 1.225, evac['person']['landing_m_s'], ROUND_CANOPY_CD0))
    canopy = dict(knacke_flat_across_m=[round(d, 1) for d in flat],
                  knacke_inflated_across_m=[round(ROUND_CANOPY_INFLATED * d, 1) for d in flat],
                  knacke_earth_flat_across_m=round(earth_flat, 1),
                  round_canopy_across_m=[evac['canopy']['base']['canopy_across_m'], evac['canopy']['top']['canopy_across_m']],
                  landing_m_s=evac['person']['landing_m_s'],
                  free_fall_m_s=[evac['canopy']['base']['free_fall_m_s'], evac['canopy']['top']['free_fall_m_s']],
                  earth_canopy_across_m=round(float(wi.canopy_diameter_m(PERSON_KG, G_EARTH, 1.225,
                                                                           evac['person']['landing_m_s'],
                                                                           evac['person']['canopy_drag_coefficient'])), 1),
                  descent_minutes=evac['descent_minutes'])
    turning = dict(moon_7_m_s_30_deg_m=round(float(wi.turn_radius_m(7.0, G_MOON, 30.0)), 0),
                   earth_10_m_s_30_deg_m=round(float(wi.turn_radius_m(10.0, G_EARTH, 30.0)), 0))
    out['people_on_wings'] = dict(canopy=canopy, earth_gliders=gliders, lunar_designs=designs, turning=turning,
                                  human_power_w=HUMAN_POWER_W, soaring=soaring(ring))

    # Electric flyers
    share, shares = calibrate_electric(EVTOLS, TECH) if EVTOLS else (TECH.airframe_share, {})
    tech = wi.Electric(**{**TECH.__dict__, 'airframe_share': share})
    out['electric'] = dict(technology=dict(tech.__dict__), calibration=shares,
                           moon=[electric_row(n, s, tech, G_MOON, city.density_kg_m3, CITY_TRIP_KM)
                                 for n, s in ELECTRIC_CLASSES.items()],
                           earth=[electric_row(n, s, tech, G_EARTH, 1.225, CITY_TRIP_KM)
                                  for n, s in ELECTRIC_CLASSES.items()],
                           buoyant_sky_boat=buoyant_boat(ELECTRIC_CLASSES['sky boat (four seats)']['payload_kg'], city))

    # Sky ferries
    full = air(winds, summit_km + FERRY_FULL_KM)
    ferries = []
    for L in FERRY_LENGTHS_M:
        r = ship(city, full, factors, L, FERRY_CRUISE_M_S, FERRY_CRUISE_M_S * FERRY_ENERGY_HOURS * 3.6,
                 FERRY_PASSENGER_KG, CITY_GUST_M_S, hotel_w=100.0)
        stop = FERRY_STOP_KM * 1e3 / FERRY_CRUISE_M_S + FERRY_DWELL_S + FERRY_MANOEUVRE_S
        r['commercial_speed_km_h'] = round(FERRY_STOP_KM / stop * 3600.0, 0)
        r['line_per_hour_per_direction'] = {f'{h:.0f} min': round(r['passengers'] * 60.0 / h, -2) for h in FERRY_HEADWAYS_MIN}
        r['turn_radius_m'] = round(ships.TURN_RADIUS_LENGTHS * L, 0)
        ferries.append(r)
    out['sky_ferries'] = ferries

    # The metropolis at the busy hour
    wing = designs['personal_wing']['moon_summit']
    electric = {r['name']: r for r in out['electric']['moon'] if r.get('closes')}
    ferry = next(f for f in ferries if f['length_m'] == 250.0)
    per_flight = flights(MODES, electric, tech, wing['best_glide_speed_m_s'], ferry, G_MOON, city.density_kg_m3)
    city_case = metropolis(CITY, MODES, per_flight, ferry)
    aloft = {k: v.get('aloft_busy_hour', 0.0) for k, v in city_case['base'].items()}
    out['metropolis'] = dict(city=CITY, modes=MODES, per_flight={k: {kk: round(vv, 3) for kk, vv in v.items()}
                                                              for k, v in per_flight.items()},
                             ferry_length_m=ferry['length_m'], cases=city_case,
                             bands=airspace(aloft, BANDS, CITY['area_km2']),
                             capacity={name: dict(**sp, per_km2=round(capacity(BANDS['sky boats and air taxis'], sp), 0),
                                                  sky_boat_share_held=round(MODES['sky boats']['share'] * capacity(
                                                      BANDS['sky boats and air taxis'], sp) / (
                                                      (aloft['sky boats'] + aloft['air taxis']) / CITY['area_km2']), 3))
                                       for name, sp in SPACINGS.items()},
                             city_power_gw=round(CITY['population'] * CITY_W_PER_PERSON / 1e9, 0))

    # The port: regional and long-haul
    pph = port['port']['occupancy']['passengers_per_hour']
    trips = port['port']['occupancy']['travel_by_trip']
    regional_pph = [x * trips['regional_km2'] / sum(trips.values()) for x in pph]
    long_pph = sky['operations']['long_haul_per_hour']
    reg_fly, reg_full = air(winds, summit_km + np.mean(REGIONAL_KM) + 2.0), air(winds, summit_km + REGIONAL_KM[1] + 4.0)
    reg_gust = ships.ring_row(ring, 30.0)['updraft_max']
    regional = []
    for L in REGIONAL_LENGTHS_M:
        r = ship(reg_fly, reg_full, factors, L, REGIONAL_CRUISE_M_S, REGIONAL_ROUTE_KM, REGIONAL_PASSENGER_KG, reg_gust,
                 hotel_w=200.0)
        r['network'] = network(regional_pph, r['passengers'], REGIONAL_MEAN_TRIP_KM, REGIONAL_CRUISE_M_S, 1.0)
        regional.append(r)
    out['regional'] = dict(passengers_per_hour=[round(x, -2) for x in regional_pph], ships=regional,
                           mean_trip_km=REGIONAL_MEAN_TRIP_KM)

    band_air, band_top = air(winds, ships.BAND_KM), air(winds, ships.FLIGHT_BAND_KM[1])
    gust = band['updraft_max_m_s']
    tower = port['port']['design']
    top_width = 320.0                    # the tower model's width at the top (summit tower study's width profile)
    width = lambda z: top_width + (tower['base_width_m'] - top_width) * max(0.0, 1.0 - z / (tower['height_km'] * 1e3)) ** tower['flare_exponent']
    band_bottom = (ships.FLIGHT_BAND_KM[0] - summit_km) * 1e3
    top = tower['height_km'] * 1e3
    in_band = top - band_bottom
    terminal = form['terminal']
    terminal_roof = terminal['z_m'] + terminal['storeys'] * terminal['storey_m']
    options, by_speed = [], []
    for L in LONG_HAUL_LENGTHS_M:
        r = ship(band_air, band_top, factors, L, ships.CRUISE_M_S, LONG_HAUL_ROUTE_KM, ships.PASSENGER_KG[1], gust)
        r['network'] = network(long_pph, r['passengers'], MEAN_TRIP_KM, ships.CRUISE_M_S, BERTH_HOURS)
        r['mooring'] = mooring(band_air, r, factors, band['design_wind_m_s'], band['service_wind_m_s'])
        berths = int(np.ceil(max(r['network']['berths']) - 1e-9))
        r['crown'] = crown_berths(max(1, berths), r, width(band_bottom), terminal_roof, top)
        r['hydrogen_at_crown_t'] = round(r['crown']['berths'] * r['hydrogen_t'], 0)
        r['trip_hours'] = dict(mean=round(MEAN_TRIP_KM * 1e3 / ships.CRUISE_M_S / 3600.0, 1),
                               longest=round(HALF_CIRCUMFERENCE_KM * 1e3 / ships.CRUISE_M_S / 3600.0, 1))
        options.append(r)
        for v in LONG_HAUL_SPEEDS_M_S:
            s = ship(band_air, band_top, factors, L, v, LONG_HAUL_ROUTE_KM, ships.PASSENGER_KG[1], gust)
            s['mean_trip_hours'] = round(MEAN_TRIP_KM * 1e3 / v / 3600.0, 1)
            by_speed.append(s)
    fast_air = air(winds, FAST_CRUISE_KM)
    fast = winged_liner(fast_air, band_air, **FAST_LINER)
    fast['network'] = network(long_pph, FAST_LINER['seats'], MEAN_TRIP_KM, fast['cruise_m_s'], 1.0)
    out['long_haul'] = dict(passengers_per_hour=long_pph, route_km=round(LONG_HAUL_ROUTE_KM, 0),
                            mean_trip_km=round(MEAN_TRIP_KM, 0), passenger_kg=ships.PASSENGER_KG[1],
                            crown=dict(top_width_m=top_width, band_bottom_above_summit_m=round(band_bottom, 0),
                                       tower_in_band_m=round(in_band, 0), top_above_summit_m=top,
                                       terminal=dict(z_m=terminal['z_m'], roof_m=terminal_roof,
                                                     storeys=terminal['storeys'], floor_m2=terminal['floor_m2']),
                                       width_at_band_bottom_m=round(width(band_bottom), 0),
                                       disks_above_summit_m=[round(d['z_m'], 0) for d in form['disks']],
                                       swing_deg=SWING_DEG, clearance_m=CLEARANCE_M, arm_beyond=ARM_BEYOND),
                            options=options, by_speed=by_speed, fast_winged=fast,
                            wind=dict(design_m_s=band['design_wind_m_s'], service_m_s=band['service_wind_m_s']))

    # Freight and platforms
    out['freighters'] = [freighter(band_air, band_top, factors, L, FREIGHT_CRUISE_M_S, LONG_HAUL_ROUTE_KM, gust)
                         for L in FREIGHT_LENGTHS_M]
    out['platforms'] = [platform(air(winds, PLATFORM_KM), winds, PLATFORM_KM, L, ships.FABRICS[0]) for L in PLATFORM_LENGTHS_M]

    out['producer'] = dict(study='sky_fleet', files={str(Path(f).relative_to(ROOT)): digest(f) for f in
                                                     (__file__, bu.__file__, wi.__file__, ships.__file__)},
                           inputs={str(Path(f).relative_to(ROOT)): digest(f) for f in (SKY_SHIPS, PORT, FORM, GLOBAL_WINDS,
                                                                                            RING)})
    return out


def main() -> int:
    r = results()
    path = HERE / 'results' / 'sky_fleet.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(r, indent=1, default=float) + '\n')
    print(path)
    return 0


if __name__ == '__main__':
    sys.exit(main())
