"""Fire safety in the central port at 0.16 g: how smoke and heat move through its spaces, when detectors and sprinklers
answer, how long people have and how long they need, and what the tower's shafts, water, sealed zone and
firefighting from the air ask.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.port_fire.run   # results/port_fire.json

It reads the summit tower study's port (its zones and their air, floors and programme) and applies the engineering
fire models (engineering/fire): Earth's fire-engineering correlations carried to lunar gravity and the port's air by
Froude scaling, EN 1991-1-2's design fires with a slower lunar bracket, the measured shift of materials' oxygen limits
at lunar gravity, and the hydraulic model of people's movement. Closed-form sums and small time-stepping loops: a
few seconds on one core.
"""
from __future__ import annotations
import hashlib
import json
import math
import sys
from pathlib import Path
import numpy as np

from atmosphere.radiative_convective import thermodynamics as th
from engineering.fire import design_fires as df
from engineering.fire import egress as eg
from engineering.fire import smoke as sm
from research.studies.summit_tower import run as tower
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.port-fire/1'
PORT = ROOT / 'research' / 'studies' / 'summit_tower' / 'results' / 'summit_tower.json'

ATM = 101325.0
INSIDE_K = 293.15                                  # conditioned air inside enclosed floors
DESIGN_AIR = th.earthlike_air(121590.0, 400.0)     # 17.46% oxygen: the design composition at every height
MOON_G = th.MOON.surface_gravity

# Sensors: smoke detectors by the 13 K ceiling-jet surrogate; quick-response sprinklers, RTI 50 (m s)^1/2, 68 C.
SPRINKLER_RTI = 50.0
SPRINKLER_C = 68.0
CLEAR_M = 2.0            # smoke kept above head height
ALARM_S = 30.0           # detection to a voice alarm through the compartment
T_MAX_S = 5400.0

# Floor plates are strips or rings about 40 m deep with open air on both faces, divided into compartments of
# 2,000 m2 (a 50 m length of strip), each with two double doors into its neighbours. The spaces: the room a fire
# starts in, a floor compartment, and a hall through three of a band's four storeys.
PLATE_DEPTH_M = 40.0
SPACES = {
    'room': dict(label='hotel or short-stay room, 30 m2 under a 3.0 m ceiling', area_m2=30.0, ceiling_m=3.0,
                 detector_m=2.0, sprinkler_m=2.5, doors_m=[0.9], travel_m=6.0, alarm_s=0.0),
    'compartment': dict(label='floor compartment, 2,000 m2 under a 3.5 m ceiling', area_m2=2000.0, ceiling_m=3.5,
                        detector_m=6.4, sprinkler_m=3.2, doors_m=[2.4, 2.4], travel_m=45.0, alarm_s=ALARM_S),
    'hall': dict(label='hall through three storeys, 5,000 m2 under a 12 m ceiling', area_m2=5000.0, ceiling_m=12.0,
                 detector_m=6.4, sprinkler_m=3.2, doors_m=[2.4] * 6, travel_m=60.0, alarm_s=ALARM_S),
}
# The room for a fully developed fire: 5 x 6 x 3 m, its door open and its window failed, lined with board.
ROOM_OPENINGS = ((0.9, 2.1), (1.8, 1.5))            # (width, height) m
ROOM_LINING_KW_M2_K = 0.03
# A floor compartment's facade for a fully developed fire: a 50 m run of 2 m-high windows, failed.
COMPARTMENT_OPENING = (50.0 * 2.0, 2.0)             # area m2, height m

# Port uses: the EN 1991-1-2 occupancy that stands for each, the space its fire starts in, and the time people take to
# start moving once alarmed (assumed ranges: awake and familiar, or awake among strangers). Hotel and short-stay fires
# start in a room, judged by room_fire.
USES = {
    'travel': dict(port_use='travel: terminals, docks, customs, lounges, ship servicing',
                   occupancy='transport (public space)', space='hall', premovement_s=(120.0, 300.0)),
    'hotels': dict(port_use='hotels', occupancy='hotel (room)', space='room', premovement_s=None),
    'short_stay': dict(port_use='short stay: compact and transit rooms', occupancy='hotel (room)', space='room',
                       premovement_s=None),
    'commerce': dict(port_use='commerce and trade', occupancy='office', space='compartment', premovement_s=(60.0, 180.0)),
    'shops': dict(port_use='shops, markets, food and drink', occupancy='shopping centre', space='compartment',
                  premovement_s=(120.0, 300.0)),
    'events': dict(port_use='events, conferences, leisure, observation', occupancy='theatre (cinema)', space='hall',
                   premovement_s=(120.0, 300.0)),
    'services': dict(port_use='services and back of house', occupancy='library', space='compartment',
                     premovement_s=(60.0, 180.0)),   # the Annex's heaviest fire load, standing in for stores
}
LUNAR_FLOW = (0.6, 1.0)    # crowd flow at lunar gravity relative to Earth's: assumed, unmeasured

# Clearing a band by air (planning figures): a berth every 60 m of rim, ships of 100 people, one every 5 minutes.
BERTH_SPACING_M, SHIP_CAPACITY, SHIP_CYCLE_S = 60.0, 100, 300.0
# Firefighting from the air: stations every 3 km up the tower, ships at 20 m/s, a minute to launch.
STATION_SPACING_M, RESPONSE_SPEED_M_S, TURNOUT_S = 3000.0, 20.0, 60.0
NFPA1710_TRAVEL_S = 240.0

# Water: NFPA 13's ordinary hazard group 2 (8.1 mm/min over 139 m2) with 946 L/min for hoses, for 90 minutes; pipe
# components rated 1.21 MPa (175 psi), or 2.07 MPa (300 psi) for high-pressure ones.
SPRINKLER_DEMAND_L_MIN = 8.1 * 139.0
HOSE_DEMAND_L_MIN = 946.0
WATER_DURATION_MIN = 90.0
PIPE_RATINGS_PA = dict(standard_175_psi=1.21e6, high_pressure_300_psi=2.07e6)

# Stairs and lift lobbies held at 12.5 Pa above the fire floor (NFPA 92, sprinklered) and no more than a person can
# open a door against (NFPA 101: 133 N; a 0.9 x 2.1 m door with a 45 N closer).
STAIR_MIN_PA = 12.5
DOOR = dict(force_n=133.0, closer_n=45.0, width_m=0.9, height_m=2.1, knob_m=0.075)
EARTH_SHAFT = dict(height_m=400.0, outside_c=-10.0)    # a tall Earth building in winter, for comparison


def digest(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def air_cases(port: dict) -> dict:
    """The air a fire meets: Earth's, and in each zone of the port. Enclosed floors hold 20 C at the outside pressure;
    the open floors take the outside air; the sealed zone is pressurised with ordinary air to the summit's pressure."""
    zones = port['zones']

    def mid(z):
        p = math.sqrt(z['air_at_bottom']['pressure_atm'] * z['air_at_top']['pressure_atm']) * ATM
        t = 0.5 * (z['air_at_bottom']['temperature_c'] + z['air_at_top']['temperature_c']) + 273.15
        return p, t

    ground = zones[0]['air_at_bottom']['pressure_atm'] * ATM
    p0, t0 = mid(zones[0])
    p1, p2 = mid(zones[1])[0], mid(zones[2])[0]
    return {
        'earth': ('Earth: 1 g, 1 atm, 20 C, 20.95% oxygen (the correlations\' own air)', sm.EARTH_REFERENCE),
        'open_0_3km': (f'open floors, 0-3 km: the outside air, {p0 / ATM:.2f} atm and {t0 - 273.15:.0f} C',
                       sm.Air(MOON_G, p0, t0, DESIGN_AIR)),
        'enclosed_3_10km': (f'enclosed floors, 3-10 km: {p1 / ATM:.2f} atm, 20 C', sm.Air(MOON_G, p1, INSIDE_K, DESIGN_AIR)),
        'enclosed_10_15km': (f'enclosed floors, 10-15 km: {p2 / ATM:.2f} atm, 20 C', sm.Air(MOON_G, p2, INSIDE_K, DESIGN_AIR)),
        'sealed_15_24km': (f'sealed zone, 15-24 km, pressurised with ordinary air to the summit\'s {ground / ATM:.2f} atm, 20 C',
                           sm.Air(MOON_G, ground, INSIDE_K, DESIGN_AIR)),
    }


def rounded(x, digits=1):
    return None if x is None or not math.isfinite(x) else round(float(x), digits)


def air_table(cases: dict) -> dict:
    earth = sm.EARTH_REFERENCE
    base_flame = float(sm.flame_height_m(1000.0, 1.0, earth))
    base_plume = float(sm.plume_mass_flow_kg_s(700.0, CLEAR_M, earth))
    out = {}
    for key, (label, air) in cases.items():
        f = sm.froude(air)
        flame = float(sm.flame_height_m(1000.0, 1.0, air))
        plume = float(sm.plume_mass_flow_kg_s(700.0, CLEAR_M, air))
        out[key] = dict(label=label, gravity_m_s2=round(air.gravity_m_s2, 4), pressure_atm=round(air.pressure_pa / ATM, 3),
                        temperature_c=round(air.temperature_k - 273.15, 1), density_kg_m3=round(air.density, 3),
                        oxygen_percent=round(100 * air.oxygen_mole_fraction, 2),
                        heat_per_kg_air_mj=round(air.heat_per_air_kj_kg / 1e3, 2),
                        acts_as_earth_fire_times=round(float(f['heat']), 2), time_factor=round(f['time'], 2),
                        velocity_factor=round(f['velocity'], 3),
                        flame_height_1mw_on_1m_m=round(flame, 2), flame_height_vs_earth=round(flame / base_flame, 2),
                        plume_at_2m_700kw_kg_s=round(plume, 2), plume_vs_earth=round(plume / base_plume, 2),
                        room_flashover_vs_earth=round(float(room_flashover(air) / room_flashover(earth)), 2),
                        ventilation_limit_vs_earth=round(float(room_ventilation(air) / room_ventilation(earth)), 2))
    return out


def room_openings():
    area = sum(w * h for w, h in ROOM_OPENINGS)
    height = sum(w * h * h for w, h in ROOM_OPENINGS) / area   # area-weighted height (EN 1991-1-2 heq)
    surface = 2 * 30.0 + 2 * (5.0 + 6.0) * 3.0 - area
    return area, height, surface


def room_flashover(air):
    area, height, surface = room_openings()
    return sm.flashover_kw(area, height, ROOM_LINING_KW_M2_K, surface, air)


def room_ventilation(air):
    area, height, _ = room_openings()
    return sm.ventilation_limit_kw(area, height, air)


def design_fires(oxygen: float) -> dict:
    out = {}
    for key, u in USES.items():
        occ = df.EN1991_ANNEX_E[u['occupancy']]
        earth, lunar = occ.earth_fire(), df.lunar_estimate(occ.earth_fire(), oxygen)
        out[key] = dict(port_use=u['port_use'], occupancy=occ.name, growth=occ.growth,
                        earth=dict(growth_time_s=earth.growth_time_s, alpha_kw_s2=round(earth.alpha_kw_s2, 5),
                                   peak_kw_m2=earth.peak_kw_m2, fire_load_mj_m2=earth.fire_load_mj_m2),
                        lunar_estimate=dict(growth_time_s=lunar.growth_time_s, alpha_kw_s2=round(lunar.alpha_kw_s2, 5),
                                            peak_kw_m2=round(lunar.peak_kw_m2), fire_load_mj_m2=lunar.fire_load_mj_m2))
    return out


def space_case(space: str, fire: df.DesignFire, air: sm.Air) -> dict:
    """Detection, sprinkler response and the time to untenable smoke for one fire in one space and air."""
    s = SPACES[space]
    alpha, area, ceiling = fire.alpha_kw_s2, s['area_m2'], s['ceiling_m']
    peak = fire.peak_kw_m2 * area
    if space == 'room':
        peak = min(peak, float(room_ventilation(air)))
    detect = float(sm.gas_rise_time_s(sm.SMOKE_DETECTOR_RISE_K, alpha, ceiling, s['detector_m'], air))
    jet = float(sm.ceiling_jet(detect, alpha, ceiling, s['detector_m'], air)[1])
    rise = SPRINKLER_C + 273.15 - air.temperature_k
    sprinkler = sm.activation_time_s(rise, SPRINKLER_RTI, alpha, ceiling, s['sprinkler_m'], air, t_max_s=T_MAX_S)
    q_sprinkler = min(fire.hrr_kw(sprinkler), peak) if math.isfinite(sprinkler) else None
    step = 0.25 if space == 'room' else 1.0
    free = sm.smoke_filling(lambda t: fire.hrr_kw(t, peak), area, ceiling, air, CLEAR_M, t_max_s=T_MAX_S, step_s=step)
    out = dict(smoke_detector_s=rounded(detect), fire_at_detection_kw=rounded(fire.hrr_kw(detect), 0),
               jet_at_detector_m_s=round(jet, 2), sprinkler_s=rounded(sprinkler),
               fire_at_sprinkler_kw=rounded(q_sprinkler, 0),
               untenable_unsprinklered_s=rounded(free['untenable_s']), unsprinklered_governed_by=free['governed_by'],
               untenable_sprinklered_s=None, sprinklered_governed_by=None,
               exhaust_to_hold_m3_s=None, natural_vents_m2=None, flame_height_at_sprinkler_m=None)
    if q_sprinkler is not None:
        held = sm.smoke_filling(lambda t: fire.hrr_kw(t, q_sprinkler), area, ceiling, air, CLEAR_M, t_max_s=T_MAX_S,
                                step_s=step)
        hold = sm.exhaust_to_hold(0.7 * q_sprinkler, CLEAR_M, air)
        diameter = math.sqrt(4.0 * q_sprinkler / (math.pi * fire.peak_kw_m2))
        out.update(untenable_sprinklered_s=rounded(held['untenable_s']), sprinklered_governed_by=held['governed_by'],
                   flame_height_at_sprinkler_m=rounded(float(sm.flame_height_m(q_sprinkler, diameter, air)), 2))
        if space != 'room':
            out.update(exhaust_to_hold_m3_s=rounded(hold['volume_m3_s']),
                       natural_vents_m2=rounded(float(sm.natural_vent_area_m2(hold['mass_kg_s'], ceiling - CLEAR_M,
                                                                              hold['layer_rise_k'], air))))
    return out


def spaces(cases: dict, oxygen: float) -> dict:
    out = {}
    for space, s in SPACES.items():
        fires = sorted({u['occupancy'] for u in USES.values() if u['space'] == space})
        rows = {}
        for name in fires:
            earth_fire = df.EN1991_ANNEX_E[name].earth_fire()
            lunar = df.lunar_estimate(earth_fire, oxygen)
            rows[name] = {key: dict(earth_fire=space_case(space, earth_fire, air),
                                    **({} if key == 'earth' else dict(lunar_estimate=space_case(space, lunar, air))))
                          for key, (_, air) in cases.items()}
        out[space] = dict(label=s['label'], fires=rows)
    return out


def room_fire(space_rows: dict, cases: dict, oxygen: float) -> dict:
    """The room a fire starts in (a hotel or short-stay room): when its own detector and sprinkler answer, when smoke
    reaches the heads of people in it, when it would flash over, and how large and long its fully developed fire is.
    Its occupants are judged, as on Earth, against the room's own alarm, its sprinkler and flashover; the rest of a
    floor is kept apart from it by its walls and door."""
    occ = df.EN1991_ANNEX_E['hotel (room)']
    load = occ.fire_load_80_mj_m2 * SPACES['room']['area_m2']
    out = {}
    for key, (_, air) in cases.items():
        row = {}
        for bracket, fire in (('earth_fire', occ.earth_fire()),
                              ('lunar_estimate', df.lunar_estimate(occ.earth_fire(), oxygen))):
            if key == 'earth' and bracket == 'lunar_estimate':
                continue
            s = space_rows['room']['fires'][occ.name][key][bracket]
            flash, vent = float(room_flashover(air)), float(room_ventilation(air))
            fuel = fire.peak_kw_m2 * SPACES['room']['area_m2']
            size = min(fuel, vent)
            row[bracket] = dict(smoke_detector_s=s['smoke_detector_s'], sprinkler_s=s['sprinkler_s'],
                                fire_at_sprinkler_kw=s['fire_at_sprinkler_kw'],
                                smoke_at_heads_s=s['untenable_unsprinklered_s'],
                                flashover_kw=round(flash), flashover_s=rounded(math.sqrt(flash / fire.alpha_kw_s2)),
                                sprinkler_fire_share_of_flashover=round(s['fire_at_sprinkler_kw'] / flash, 2),
                                fuel_limit_kw=round(fuel), ventilation_limit_kw=round(vent),
                                controlled_by='ventilation' if vent < fuel else 'fuel', fire_kw=round(size),
                                burning_min=round(load * 1e3 / size / 60.0))
        out[key] = row
    return out


def zone_winds(port_results: dict) -> dict:
    """The ring's 3-hourly median and 90th-percentile wind at the middle of each zone (the tower study's climate)."""
    rows = port_results['wind']
    heights = [r['height_above_ground_km'] for r in rows]
    out = {}
    for z in port_results['port']['zones']:
        mid = 0.5 * (z['from_km'] + z['to_km'])
        out[f"{z['from_km']:g}_{z['to_km']:g}km"] = {
            p: round(float(np.interp(mid, heights, [r['ring_3hourly'][p] for r in rows])), 2) for p in ('median', 'p90')}
    return out


ZONE_OF_CASE = {'open_0_3km': '0_3km', 'enclosed_3_10km': '3_10km', 'enclosed_10_15km': '10_15km',
                'sealed_15_24km': '15_24km'}


def compartment_burnout(cases: dict, winds: dict) -> dict:
    """A floor compartment fully involved behind a failed 50 m window run on each face: the heat release its
    openings feed and how long its fire load lasts, in still air (buoyancy alone) and with the zone's median and
    90th-percentile wind blowing through. Earth's design fire; the fuel-controlled peak caps it."""
    area, height = COMPARTMENT_OPENING
    floor = SPACES['compartment']['area_m2']
    out = {}
    for use in ('commerce', 'shops', 'services'):
        occ = df.EN1991_ANNEX_E[USES[use]['occupancy']]
        load, fuel = occ.fire_load_80_mj_m2 * floor, occ.peak_kw_m2 * floor
        row = {}
        for key, (_, air) in cases.items():
            options = dict(still_air=float(sm.ventilation_limit_kw(area, height, air)))
            if key in ZONE_OF_CASE:
                for p, v in winds[ZONE_OF_CASE[key]].items():
                    options[f'{p}_wind'] = max(options['still_air'], float(sm.wind_ventilation_limit_kw(area, v, air)))
            row[key] = {k: dict(fire_mw=round(min(q, fuel) / 1e3), burning_h=round(load * 1e3 / min(q, fuel) / 3600.0, 1),
                                controlled_by='ventilation' if q < fuel else 'fuel') for k, q in options.items()}
        out[use] = dict(occupancy=occ.name, fire_load_gj=round(load / 1e3), fuel_controlled_mw=round(fuel / 1e3),
                        by_air=row)
    return out


def egress(space_rows: dict) -> dict:
    """Leaving the compartment or hall a fire starts in: RSET against the time to untenable smoke, on Earth and in
    the enclosed floors at 3-10 km, with Earth's design fires."""
    out = {}
    for use, u in USES.items():
        if u['space'] == 'room':
            continue
        s = SPACES[u['space']]
        dense = tower.PORT_USES[u['port_use']][1][0]
        people = s['area_m2'] / dense
        row = dict(space=u['space'], people=round(people), premovement_s=list(u['premovement_s']))
        for key, flows in (('earth', (1.0, 1.0)), ('enclosed_3_10km', LUNAR_FLOW)):
            case = space_rows[u['space']]['fires'][u['occupancy']][key]['earth_fire']
            detect = case['smoke_detector_s']
            quick = eg.movement_s(people, s['travel_m'], s['doors_m'], flows[1])
            slow = eg.movement_s(people, s['travel_m'], s['doors_m'], flows[0])
            rset = (eg.rset_s(detect, s['alarm_s'], u['premovement_s'][0], quick['movement_s']),
                    eg.rset_s(detect, s['alarm_s'], u['premovement_s'][1], slow['movement_s']))
            aset = case['untenable_sprinklered_s']
            row[key] = dict(detection_s=detect, movement_s=[round(quick['movement_s']), round(slow['movement_s'])],
                            rset_s=[round(v) for v in rset], untenable_sprinklered_s=aset,
                            untenable_unsprinklered_s=case['untenable_unsprinklered_s'],
                            margin_sprinklered_s=None if aset is None else round(aset - rset[1]),
                            margin_unsprinklered_s=None if case['untenable_unsprinklered_s'] is None
                            else round(case['untenable_unsprinklered_s'] - rset[1]))
        out[use] = row
    return out


def width_at(port: dict, height_m: float) -> float:
    d = port['design']
    top = tower.FRAME.top_width_m
    return top + (d['base_width_m'] - top) * (1.0 - height_m / (d['height_km'] * 1e3)) ** d['flare_exponent']


def band_clearing(port: dict) -> dict:
    """Clearing a whole band by air, at the middle of each zone: storey area, people when every space is full of
    each band type, and the time to lift them off from the rim."""
    floors = tower.PORT_FLOORS
    out = {}
    for z in port['zones']:
        height = 0.5 * (z['from_km'] + z['to_km']) * 1e3
        storey = float(min(floors.area_per_storey_m2, floors.plan_share * width_at(port, height) ** 2))
        rim = floors.storeys * storey / PLATE_DEPTH_M
        people = {use: round(floors.storeys * storey / tower.PORT_USES[u['port_use']][1][0])
                  for use, u in USES.items()}
        worst = max(people.values())
        lift = eg.aerial_evacuation(worst, rim, BERTH_SPACING_M, SHIP_CAPACITY, SHIP_CYCLE_S)
        out[f"{z['from_km']:g}_{z['to_km']:g}km"] = dict(
            storey_m2=round(storey), band_m2=round(floors.storeys * storey), rim_m=round(rim), people_when_full=people,
            berths=lift['berths'], persons_per_s=round(lift['persons_per_s'], 1),
            fullest_band_clear_min=round(lift['clear_s'] / 60.0),
            clear_min_by_use={use: round(n / lift['persons_per_s'] / 60.0, 1) for use, n in people.items()})
    return out


def door_limit_pa() -> float:
    d = DOOR
    return 2.0 * (d['force_n'] - d['closer_n']) * (d['width_m'] - d['knob_m']) / (d['width_m'] * d['width_m'] * d['height_m'])


def systems(port: dict, cases: dict, oxygen: float, winds: dict) -> dict:
    zones = port['zones']
    # Shafts: the stack effect with 20 C inside and the zone's outside air at its middle.
    window = door_limit_pa() - STAIR_MIN_PA
    earth_out = sm.Air(sm.EARTH_REFERENCE.gravity_m_s2, ATM, EARTH_SHAFT['outside_c'] + 273.15)
    earth_dp = float(sm.stack_pressure_pa(EARTH_SHAFT['height_m'], earth_out, INSIDE_K))
    shafts = dict(door_limit_pa=round(door_limit_pa()), stair_minimum_pa=STAIR_MIN_PA,
                  earth_comparison=dict(shaft_m=EARTH_SHAFT['height_m'], outside_c=EARTH_SHAFT['outside_c'],
                                        stack_pa=round(earth_dp),
                                        unbroken_shaft_within_door_window_m=round(window * EARTH_SHAFT['height_m'] / earth_dp)),
                  zones={})
    for z in zones:
        p = math.sqrt(z['air_at_bottom']['pressure_atm'] * z['air_at_top']['pressure_atm']) * ATM
        t = 0.5 * (z['air_at_bottom']['temperature_c'] + z['air_at_top']['temperature_c']) + 273.15
        outside = sm.Air(MOON_G, p, t, DESIGN_AIR)
        per_m = float(sm.stack_pressure_pa(1.0, outside, INSIDE_K))
        span = (z['to_km'] - z['from_km']) * 1e3
        shafts['zones'][f"{z['from_km']:g}_{z['to_km']:g}km"] = dict(
            outside_c=round(t - 273.15, 1), stack_pa_per_100m=round(100 * per_m, 2),
            stack_pa_per_band_spacing=round(tower.PORT_FLOORS.spacing_m * per_m, 1),
            stack_pa_whole_zone=round(span * per_m), unbroken_shaft_within_door_window_m=round(window / per_m, -1))
    # Water.
    water = dict(pressure_kpa_per_100m=dict(moon=round(float(sm.water_pressure_pa(100.0, MOON_G)) / 1e3, 1),
                                            earth=round(float(sm.water_pressure_pa(100.0, sm.EARTH_REFERENCE.gravity_m_s2)) / 1e3, 1)),
                 pressure_zone_m={}, zones_over_tower={})
    height = port['design']['height_km'] * 1e3
    for label, rating in PIPE_RATINGS_PA.items():
        moon = rating / float(sm.water_pressure_pa(1.0, MOON_G))
        earth = rating / float(sm.water_pressure_pa(1.0, sm.EARTH_REFERENCE.gravity_m_s2))
        water['pressure_zone_m'][label] = dict(moon=round(moon), earth=round(earth))
        water['zones_over_tower'][label] = dict(moon=math.ceil(height / moon), earth_gravity=math.ceil(height / earth))
    demand = SPRINKLER_DEMAND_L_MIN + HOSE_DEMAND_L_MIN
    water.update(demand_l_min=round(demand), tank_per_band_m3=round(demand * WATER_DURATION_MIN / 1e3),
                 lift_kwh_per_m3_to_top=round(float(sm.water_pressure_pa(height, MOON_G)) / 3.6e6, 1))
    # The sealed zone: its envelope when pressurised, and the air people meet if it fails.
    ground = zones[0]['air_at_bottom']['pressure_atm']
    sealed = zones[-1]
    sealed_zone = dict(pressurised_to_atm=ground,
                       envelope_kpa=[round((ground - sealed['air_at_bottom']['pressure_atm']) * ATM / 1e3, 1),
                                     round((ground - sealed['air_at_top']['pressure_atm']) * ATM / 1e3, 1)],
                       outside_oxygen_like_earth_at_m=[sealed['air_at_bottom']['oxygen_like_earth_at_m'],
                                                       sealed['air_at_top']['oxygen_like_earth_at_m']],
                       faa_effective_performance_min_at_18000_ft=[20, 30])
    # Topping up by enrichment instead: the oxygen fraction that gives each zone's top the summit's oxygen partial
    # pressure, or Earth's at sea level, without raising the pressure.
    summit_po2 = ground * 100 * oxygen
    enrichment = {f"{z['from_km']:g}_{z['to_km']:g}km": dict(
        pressure_atm=z['air_at_top']['pressure_atm'],
        oxygen_percent_for_summit_breathing=round(summit_po2 / z['air_at_top']['pressure_atm'], 1),
        oxygen_percent_for_earth_sea_level=round(20.95 / z['air_at_top']['pressure_atm'], 1)) for z in zones[1:]}
    # Firefighting from the air.
    reach = STATION_SPACING_M / 2.0
    firefighting = dict(stations=math.ceil(height / STATION_SPACING_M), farthest_m=reach,
                        arrival_s=round(TURNOUT_S + reach / RESPONSE_SPEED_M_S),
                        nfpa1710_first_engine_travel_s=NFPA1710_TRAVEL_S)
    return dict(shafts=shafts, water=water, sealed_zone=sealed_zone, enrichment_instead=enrichment,
                firefighting=firefighting, compartment_burnout=compartment_burnout(cases, winds))


def materials(oxygen: float) -> dict:
    shifts = df.lunar_oxygen_shifts()
    margin = max(max(s.values()) for s in shifts.values())
    return dict(limits_mole_percent=df.LUNAR_OXYGEN_LIMITS, shifts_points=shifts, margin_points=margin,
                design_oxygen_percent=round(100 * oxygen, 2),
                test_oxygen_percent=df.materials_test_oxygen_percent(oxygen, margin),
                burning_rate_factor=round(df.burning_rate_factor(oxygen), 3))


def results() -> dict:
    tower_results = json.loads(PORT.read_text())
    port = tower_results['port']
    oxygen = DESIGN_AIR.fractions['O2']
    cases = air_cases(port)
    space_rows = spaces(cases, oxygen)
    winds = zone_winds(tower_results)
    return dict(
        port=dict(height_km=port['design']['height_km'], floor_km2=port['floor_km2'],
                  floors=dict(storeys_per_band=tower.PORT_FLOORS.storeys, band_spacing_m=tower.PORT_FLOORS.spacing_m,
                              storey_height_m=tower.PORT_FLOORS.storey_height_m, plate_depth_m=PLATE_DEPTH_M,
                              compartment_m2=SPACES['compartment']['area_m2'])),
        air=air_table(cases),
        design_fires=design_fires(oxygen),
        materials=materials(oxygen),
        spaces=space_rows,
        room_fire=room_fire(space_rows, cases, oxygen),
        egress=egress(space_rows),
        band_clearing=band_clearing(port),
        winds_m_s=winds,
        systems=systems(port, cases, oxygen, winds))


def main() -> int:
    r = results()
    files = {p: digest(ROOT / p) for p in ('research/studies/port_fire/run.py', 'engineering/fire/smoke.py',
                                           'engineering/fire/design_fires.py', 'engineering/fire/egress.py',
                                           'research/studies/summit_tower/run.py',
                                           'atmosphere/radiative_convective/thermodynamics.py')}
    product = dict(schema=SCHEMA, producer=dict(study='port_fire', files=files, constants=constants_used(files),
                                                inputs={str(PORT.relative_to(ROOT)): digest(PORT)}),
                   evidence=('First-order fire-engineering arithmetic. Earth correlations (plumes, ceiling jets, '
                             'detector and sprinkler response, smoke filling, flashover, ventilation limits) are '
                             'carried to lunar gravity and the port\'s air by Froude scaling; none has been tested '
                             'at partial gravity beyond small samples. Design fires are EN 1991-1-2\'s, with a slower '
                             'lunar bracket built from Earth oxygen data and short low-gravity tests. Crowd flows at '
                             'lunar gravity, pre-movement times and the aerial evacuation figures are assumptions. '
                             'Nothing here is a fire strategy or a certified design.'),
                   reading_rule=('Times in seconds from ignition. untenable_* is when smoke reaches 2 m above the floor '
                                 'or the layer reaches 200 C; null or absent where a sensor never answers. air.*.'
                                 'acts_as_earth_fire_times: in the same space, a fire of Q behaves as an Earth fire of '
                                 'that many times Q, slowed by time_factor. earth_fire is EN 1991-1-2\'s design fire; '
                                 'lunar_estimate its slower bracket.'), **r)
    out = HERE / 'results' / 'port_fire.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(product, indent=1) + '\n')
    print(out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
