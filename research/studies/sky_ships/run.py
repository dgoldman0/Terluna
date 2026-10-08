"""How large a sky ship the Open Moon's air can reasonably carry: the air ships fly in, the rule of similarity that
lunar gravity sets, how buoyant and winged ships' structure grows with size, the atmosphere's own scale in its
storms, and what else grows with a ship.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.sky_ships.run   # results/sky_ships.json

It couples the climate products (the design run's winds and air over the whole Moon from
climate/gcm/global_winds.py, and the cloud-resolving equatorial ring's winds, vertical motion and storms), the
summit tower study's port (the long-haul traffic at the crown) and the engineering flight models
(engineering/flight). Everything is closed-form sizing over a few thousand candidate ships: a second on one core.
"""
from __future__ import annotations
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

from atmosphere.radiative_convective.thermodynamics import MOON, earthlike_air
from engineering.flight import buoyant as bu
from engineering.flight import winged as wi
from research.studies.summit_tower import run as tower
from research.studies.summit_tower.run import earth_air
from shared.constants import STANDARD_GRAVITY
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.sky-ships/1'
GLOBAL_WINDS = ROOT / 'climate' / 'results' / 'gcm' / 'global_winds_A28_dim5_moon.json'
RING = ROOT / 'climate' / 'results' / 'crm' / 'ring_ring_equator.json'
PORT = ROOT / 'research' / 'studies' / 'summit_tower' / 'results' / 'summit_tower.json'

PRODUCER_FILES = (__file__, bu.__file__, wi.__file__, tower.__file__,
                  ROOT / 'atmosphere' / 'radiative_convective' / 'thermodynamics.py')
G_MOON = MOON.surface_gravity                       # GM/R^2, 1.6242 m/s2
G_EARTH = STANDARD_GRAVITY
AIR = earthlike_air(121590.0, 400.0)                 # the design air's composition (research/decisions.md: 1.2 atm)
AIR_MOLAR_MASS = AIR.molar_mass
EARTH_MOLAR_MASS = 0.0289644
SIMILAR = G_EARTH / G_MOON                           # the length factor that keeps weight-driven shares equal

AIR_HEIGHTS_KM = (0.0, 10.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 55.0, 70.0)
FLIGHT_BAND_KM = (35.0, 45.0)
BAND_KM = 40.0                                       # the middle of the flight band, for sizing
PURITY = 0.98                                        # lifting gas share by volume in the cells
HYDROGEN_HEAT_J_KG = 120e6                           # lower heating value
GUST_FACTOR, RARITY_FACTOR = 1.4, 1.2                # as the summit tower study: 3-second gust, decades
EVIDENCE = ('First-order sizing of buoyant and winged ships with the engineering flight models, fitted to the weight '
            'statements of four rigid airships and the wing groups of ten transports, in the design run\'s air and winds over the whole Moon '
            '(ExoPlaSim, 3-day means) and the cloud-resolving equatorial ring (3-hourly, with storms). No ship has '
            'been designed; loads follow Earth practice (Woodward\'s gust envelope, a deflated gas cell, Pratt\'s gust '
            'formula) carried to lunar gravity.')
READING_RULE = ('Masses in kg or t as labelled, lengths in m, shares of gross lift (buoyant) or of flyer mass '
                '(winged). Heights are above the model\'s sea level. The flight band is 35-45 km; ships in it are '
                'sized in the air at 40 km. "useful share" is gross lift less the empty ship, over gross lift.')


def digest(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def load():
    winds = json.loads(GLOBAL_WINDS.read_text())
    ring = json.loads(RING.read_text())
    port = json.loads(PORT.read_text())
    assert winds['schema'] == 'terluna.climate.gcm-global-winds/1'
    return winds, ring, port


def at(winds: dict, key: str, z_km: float) -> float:
    return float(np.interp(z_km, winds['height_km'], np.asarray(winds[key], dtype=float)))


def moon_air(winds: dict, z_km: float) -> bu.Air:
    return bu.Air(G_MOON, at(winds, 'density_mean_kg_m3', z_km), at(winds, 'temperature_k', z_km), AIR_MOLAR_MASS)


def earth_height_of_density(rho: float) -> float | None:
    """The height (km) in the U.S. Standard Atmosphere with this density, from 2 km below sea level to 32 km."""
    z = np.linspace(-2e3, 32e3, 3401)
    d = earth_air(z)
    if rho > d[0] or rho < d[-1]:
        return None
    return float(np.interp(rho, d[::-1], z[::-1]) / 1e3)


def ring_row(ring: dict, z_km: float) -> dict:
    heights = {float(k.split('_')[0]): k for k in ring['wind_aloft_m_s']}
    near = min(heights, key=lambda h: abs(h - z_km))
    w, v = ring['wind_aloft_m_s'][heights[near]], ring['vertical_wind_m_s'][heights[near]]
    return dict(ring_height_km=near, median=round(w['median'], 1), p99=round(w['p99'], 1), max=round(w['max'], 1),
                vertical_p99=round(v['abs_p99'], 2), updraft_max=round(v['up_max'], 1))


def air_table(winds: dict, ring: dict) -> list:
    rows = []
    for z in AIR_HEIGHTS_KM:
        air = moon_air(winds, z)
        rows.append(dict(
            height_km=z, density_kg_m3=round(air.density_kg_m3, 3),
            density_range_kg_m3=[round(at(winds, 'density_p1_kg_m3', z), 3), round(at(winds, 'density_p99_kg_m3', z), 3)],
            pressure_atm=round(at(winds, 'pressure_pa', z) / 101325.0, 3), temperature_c=round(air.temperature_k - 273.15, 1),
            lift_hydrogen_kg_m3=round(float(bu.net_lift_kg_m3(air, bu.HYDROGEN, PURITY)), 3),
            lift_helium_kg_m3=round(float(bu.net_lift_kg_m3(air, bu.HELIUM, PURITY)), 3),
            earth_height_same_density_km=(None if (e := earth_height_of_density(air.density_kg_m3)) is None
                                          else round(e, 1)),
            gcm_3day=dict(median=round(at(winds, 'speed_p50_m_s', z), 1), p99=round(at(winds, 'speed_p99_m_s', z), 1),
                          max=round(at(winds, 'speed_max_m_s', z), 1),
                          return_50yr=round(at(winds, 'return_50yr_m_s', z), 1)),
            ring_3hourly=ring_row(ring, z)))
    return rows


def band_winds(winds: dict, ring: dict) -> dict:
    """The flight band's winds: the models' highest, and the design and service winds for a moored ship, raised by
    the summit tower study's factors for a 3-second gust and for decades."""
    z = np.array(winds['height_km'])
    inband = (z >= FLIGHT_BAND_KM[0]) & (z <= FLIGHT_BAND_KM[1])
    gcm_max = float(np.nanmax(np.asarray(winds['speed_max_m_s'], dtype=float)[inband]))
    gcm_p99 = float(np.nanmax(np.asarray(winds['speed_p99_m_s'], dtype=float)[inband]))
    r = ring_row(ring, BAND_KM)
    highest = max(gcm_max, r['max'])
    return dict(gcm_3day_max_m_s=round(gcm_max, 1), gcm_3day_p99_m_s=round(gcm_p99, 1), ring_3hourly_max_m_s=r['max'],
                ring_3hourly_p99_m_s=r['p99'], models_highest_m_s=round(highest, 1),
                service_wind_m_s=round(max(gcm_p99, r['p99']) * GUST_FACTOR, 1),
                design_wind_m_s=round(highest * GUST_FACTOR * RARITY_FACTOR, 1),
                updraft_max_m_s=r['updraft_max'], vertical_p99_m_s=r['vertical_p99'],
                storm_cloud_share=dict(land=round(ring['cloud_in_flight_band']['land'], 4),
                                       water=round(ring['cloud_in_flight_band']['water'], 4)))


def storms(ring: dict) -> dict:
    s, t = ring['storms'], ring['storm_tracks']
    pick = lambda d: {k: round(v, 1) for k, v in d.items()}
    return dict(width_km=pick(s['width_km']), cloud_top_km=pick(s['cloud_top_km']),
                updraft_max_m_s=pick(s['updraft_max_m_s']), lifetime_h=t['lifetime_h'],
                drift_m_s=pick(t['drift_m_s']), per_snapshot=round(s['per_snapshot'], 1))


# ---------------------------------------------------------------------------------------------------------------
# Buoyant ships by size

EARTH_SEA = bu.Air(G_EARTH, 1.225, 288.15, EARTH_MOLAR_MASS)
LENGTHS_M = np.geomspace(30.0, 30e3, 241)


def scan(air: bu.Air, design: bu.Design, lift_per_m3: float, lengths=LENGTHS_M, fineness=6.0) -> dict:
    """Rigid hulls of the given lengths flying in `air` with a gross lift per cubic metre of hull."""
    r = bu.rigid(bu.Hull(np.asarray(lengths, dtype=float), fineness), air, lift_per_m3, design)
    return dict(lift_per_m3=lift_per_m3, **r)


def crossing(x, y, level):
    """First x where y falls through (or rises through) a level, by linear interpolation; None if never."""
    y = np.asarray(y, dtype=float) - level
    k = np.nonzero(np.sign(y[:-1]) != np.sign(y[1:]))[0]
    if k.size == 0:
        return None
    i = k[0]
    return float(x[i] - y[i] * (x[i + 1] - x[i]) / (y[i + 1] - y[i]))


def landmarks(r: dict, reference_useful: float) -> dict:
    """Lengths where the ship's structure changes character: where the weight-driven structure overtakes the
    aerodynamic, where the useful lift peaks, where it falls back to a reference ship's share, and where it is gone."""
    L, u = r['length_m'], r['useful_share']
    best = int(np.argmax(u))
    after = slice(best, None)
    return dict(weight_overtakes_gusts_m=crossing(L, r['weight_driven_share'] - r['aerodynamic_share'], 0.0),
                best_length_m=float(L[best]), best_useful_share=float(u[best]),
                back_to_reference_m=crossing(L[after], u[after], reference_useful),
                useful_gone_m=crossing(L[after], u[after], 0.0))


# ---------------------------------------------------------------------------------------------------------------
# Pressure hulls: the largest envelope a fabric holds

def largest_pressure_hull(air: bu.Air, fabric: bu.Fabric, speed_m_s: float, gust_m_s: float, fineness=4.0,
                          safety=4.0) -> dict:
    """The largest non-rigid hull (by diameter) whose crown tension, times the fabric's safety factor, the fabric
    carries at a top speed and gust."""
    lift = float(bu.net_lift_kg_m3(air, bu.HYDROGEN, PURITY))
    design = bu.Design(speed_m_s=speed_m_s, gust_m_s=gust_m_s)
    L = np.geomspace(10.0, 20e3, 400)
    r = bu.pressure_hull(bu.Hull(L, fineness), air, lift, fabric, design, fabric_safety=safety)
    ok = r['fabric_use'] <= 1.0
    if not ok.any():
        return dict(length_m=None)
    i = int(np.nonzero(ok)[0].max())
    h = bu.Hull(L[i], fineness)
    return dict(length_m=round(float(L[i]), -1), diameter_m=round(float(h.diameter_m), 0),
                volume_m3=float(f'{h.volume_m3:.3g}'), gross_lift_t=round(float(lift * h.volume_m3 / 1e3), 0),
                pressure_pa=round(float(r['pressure_pa'][i]), 0),
                dynamic_pressure_pa=round(0.5 * air.density_kg_m3 * speed_m_s ** 2, 0))


# ---------------------------------------------------------------------------------------------------------------
# Winged ships by span

SPANS_M = np.geomspace(5.0, 5000.0, 301)


def wing_scan(gravity: float, density: float, wing_loading_pa: float, gust_true_m_s: float, wing_model: dict,
              spans=SPANS_M, aspect_ratio=9.0, taper=0.3, thickness=0.13, lift_coefficient=0.5,
              manoeuvre=2.5, material=wi.ALUMINIUM) -> dict:
    """The wing's structural share of the flyer's mass against span, for flyers of one wing loading (in pascals)
    cruising at a lift coefficient in air of one density; the load factor is the manoeuvre factor or the gust's,
    whichever is larger, times 1.5."""
    b = np.asarray(spans, dtype=float)
    wing = wi.Wing(b, aspect_ratio, taper, thickness)
    mass = wing_loading_pa * wing.area_m2 / gravity
    per_area = mass / wing.area_m2
    speed = float(wi.speed_for_lift(per_area[0], gravity, density, lift_coefficient))
    dn = wi.gust_load_factor(per_area, wing.mean_chord_m, speed, gust_true_m_s, density, gravity)
    n_ult = 1.5 * np.maximum(manoeuvre, 1.0 + dn)
    share = wi.wing_mass_kg(wing, mass, gravity, n_ult, material, **wing_model) / mass
    return dict(span_m=b, mass_kg=mass, speed_m_s=speed, ultimate_load_factor=n_ult, wing_share=share)


# Wing groups of transports (Roskam, Airplane Design Part V, Appendix A: wing group, flight design gross weight and wing
# area, in lb and ft2; ultimate load factor 3.75 for all) with their spans (the 747's from NASA TP-1625, the C-5's from
# Lockheed Martin, the rest from the aircraft's published specifications; sources.json).
LB, FT2, FT = 0.45359237, 0.09290304, 0.3048
TRANSPORTS = {name: dict(wing_kg=w * LB, mass_kg=g * LB, area_m2=s * FT2, span_m=b) for name, (w, g, s, b) in {
    'Boeing 747-100': (86_402, 710_000, 5_500, 195.7 * FT),
    'Lockheed C-5A': (100_015, 769_000, 6_200, 222.8 * FT),
    'Lockheed C-141B': (35_272, 314_200, 3_228, 48.74),
    'Boeing KC-135': (25_251, 297_000, 2_435, 39.88),
    'Airbus A300 B2': (44_131, 302_000, 2_799, 44.84),
    'Boeing 707-321': (28_647, 302_000, 2_892, 43.41),
    'Douglas DC-8': (27_556, 215_000, 2_773, 43.41),
    'Boeing 727-100': (17_764, 160_000, 1_700, 32.92),
    'Boeing 737-200': (10_613, 115_500, 980, 28.35),
    'Douglas DC-9-10': (9_470, 91_500, 934, 27.25)}.items()}
EARTH_CRUISE_DENSITY = 0.38                # about 10.5 km
EARTH_GUST_EAS = 15.24                     # 50 ft/s derived gust at cruise (the older FAR 25.341), equivalent airspeed
REFERENCE_FLYER = dict(name='An-225', span_m=88.4, area_m2=905.0, mass_kg=640e3)
WING_GUST_BAND_M_S = 20.0                  # Earth's 50 ft/s at the band's density (19.7 m/s true) and the ring's 14.3 m/s x 1.4


def calibrate_wings(transports=TRANSPORTS) -> tuple:
    """Least-squares fit, in relative error, of the non-optimum factor and the areal term to the transports' wing
    groups; returns the model and each transport's modelled wing against its own."""
    rows, caps_kg = [], {}
    for name, t in transports.items():
        wing = wi.Wing(t['span_m'], t['span_m'] ** 2 / t['area_m2'], 0.3, 0.13)
        caps_kg[name] = float(wi.wing_mass_kg(wing, t['mass_kg'], G_EARTH, 3.75, non_optimum=1.0, areal_kg_m2=0.0))
        rows.append([caps_kg[name] / t['wing_kg'], t['area_m2'] / t['wing_kg']])
    (k, w), *_ = np.linalg.lstsq(np.array(rows), np.ones(len(rows)), rcond=None)
    check = {name: dict(wing_t=round(t['wing_kg'] / 1e3, 1), model_t=round((k * caps_kg[name] + w * t['area_m2']) / 1e3, 1),
                        share=round(t['wing_kg'] / t['mass_kg'], 3)) for name, t in transports.items()}
    return dict(non_optimum=float(k), areal_kg_m2=float(w)), check


def winged_results(band_density: float, band_gust: float) -> dict:
    """The wing's share against span on Earth (at a transport's cruise) and on the Moon (in the flight band), for
    flyers of the reference flyer's wing loading, and the spans where the Moon's share matches the reference's."""
    model, check = calibrate_wings()
    ref = REFERENCE_FLYER
    loading = ref['mass_kg'] * G_EARTH / ref['area_m2']
    ar = ref['span_m'] ** 2 / ref['area_m2']
    earth_gust = EARTH_GUST_EAS * (1.225 / EARTH_CRUISE_DENSITY) ** 0.5
    earth = wing_scan(G_EARTH, EARTH_CRUISE_DENSITY, loading, earth_gust, model, aspect_ratio=ar)
    moon = wing_scan(G_MOON, band_density, loading, band_gust, model, aspect_ratio=ar)
    same_air = wing_scan(G_MOON, EARTH_CRUISE_DENSITY, loading, earth_gust, model, aspect_ratio=ar)
    ref_share = float(np.interp(ref['span_m'], earth['span_m'], earth['wing_share']))
    match = lambda r: crossing(r['span_m'], r['wing_share'], ref_share)
    span_band, span_same = match(moon), match(same_air)
    mass_at = lambda b: loading * b ** 2 / ar / G_MOON
    return dict(model={k: round(v, 3) for k, v in model.items()}, calibration=check,
                reference=dict(**ref, wing_loading_pa=round(loading, 0), wing_share=round(ref_share, 3),
                                             cruise_m_s=round(earth['speed_m_s'], 0)),
                moon_same_air=dict(span_m=round(span_same, 0), mass_t=round(mass_at(span_same) / 1e3, -2),
                                   cruise_m_s=round(same_air['speed_m_s'], 0)),
                moon_flight_band=dict(span_m=round(span_band, 0), mass_t=round(mass_at(span_band) / 1e3, -2),
                                      cruise_m_s=round(moon['speed_m_s'], 0),
                                      ultimate_load_factor=round(float(np.interp(span_band, moon['span_m'],
                                                                                 moon['ultimate_load_factor'])), 2)),
                curves=dict(span_m=[round(float(x), 1) for x in earth['span_m']],
                            earth=[round(float(x), 4) for x in earth['wing_share']],
                            moon_flight_band=[round(float(x), 4) for x in moon['wing_share']]))


HULL_LIFT_COEFFICIENT = 0.3              # a lifting hull's lift on its planform at a few degrees' pitch
PLANFORM_SHARE = 0.7                     # planform / (length x diameter) of the streamlined hull


def hull_lift_ratio(air: bu.Air, lift_per_m3: float, speed_m_s: float, diameters_m) -> list:
    """A hybrid's aerodynamic lift at cruise over its buoyant lift, for hulls of fineness 6: q C_L (0.7 L D) over
    (Δρ g 0.475 L D^2), which falls as 1/D."""
    q = 0.5 * air.density_kg_m3 * speed_m_s ** 2
    c = bu.Hull.volume_coefficient
    return [dict(diameter_m=d, ratio=round(q * HULL_LIFT_COEFFICIENT * PLANFORM_SHARE /
                                            (lift_per_m3 * air.gravity_m_s2 * c * d), 2)) for d in diameters_m]


# ---------------------------------------------------------------------------------------------------------------
# What grows with a ship

PASSENGER_KG = (250.0, 400.0, 600.0)     # per passenger: person and baggage 100 kg, cabin, services and crew share
CRUISE_M_S = 30.0
TRIP_HOURS = 48.0                        # the longest trips: half the Moon's circumference at cruise speed is 50 h
FUEL_CELL_EFFICIENCY = 0.5
TANK_SHARE = 0.3                         # hydrogen fuel as a share of fuel and its tanks
BERTH_HOURS = 2.0                        # a call: unload, service, load
BEAM_DRAG = 0.6                          # side force coefficient on the hull's side area (0.8 L D) before it turns


HOTEL_W_PER_PASSENGER = 500.0            # light, heat, water and galleys aboard


def operations(air: bu.Air, design: bu.Design, lengths, design_wind: float, service_wind: float,
               long_haul_per_hour, fineness=6.0, lift_air=None) -> list:
    lift = float(bu.net_lift_kg_m3(lift_air or air, bu.HYDROGEN, PURITY))
    rho_h2 = bu.gas_density(lift_air or air, bu.HYDROGEN, 1.0)     # cells full at the top of the working heights
    rho_he = bu.gas_density(lift_air or air, bu.HELIUM, 1.0)
    hindenburg_h2 = 200_000.0 * bu.gas_density(EARTH_SEA, bu.HYDROGEN, 1.0) * PURITY
    rows = []
    for L in lengths:
        hull = bu.Hull(float(L), fineness)
        r = bu.rigid(hull, air, lift, design)
        cruise = float(bu.power_w(hull, air, CRUISE_M_S, design))
        per_joule = 1.0 / (FUEL_CELL_EFFICIENCY * HYDROGEN_HEAT_J_KG * TANK_SHARE)   # kg of fuel and tanks per J
        fuel = cruise * TRIP_HOURS * 3600.0 * per_joule
        hotel = HOTEL_W_PER_PASSENGER * TRIP_HOURS * 3600.0 * per_joule                 # per passenger
        passengers = [(float(r['useful_kg']) - fuel) / (m + hotel) for m in PASSENGER_KG]
        fuel += hotel * passengers[1]
        payload = float(r['useful_kg']) - fuel
        h2 = hull.volume_m3 * PURITY * rho_h2
        q_design, q_service = 0.5 * air.density_kg_m3 * design_wind ** 2, 0.5 * air.density_kg_m3 * service_wind ** 2
        side = 0.8 * hull.length_m * hull.diameter_m
        n = passengers[1]
        rows.append(dict(
            length_m=float(L), diameter_m=round(hull.diameter_m, 1), volume_m3=float(f'{hull.volume_m3:.3g}'),
            gross_lift_t=round(float(r['gross_lift_kg']) / 1e3, 1), empty_t=round(float(r['empty_kg']) / 1e3, 1),
            useful_share=round(float(r['useful_share']), 3), trip_fuel_and_tanks_t=round(fuel / 1e3, 1),
            payload_t=round(payload / 1e3, 1), passengers=[round(x, -1) for x in passengers[::-1]],
            hydrogen_t=round(h2 / 1e3, 1), hydrogen_heat_tj=round(h2 * HYDROGEN_HEAT_J_KG / 1e12, 1),
            hydrogen_in_hindenburgs=round(h2 / hindenburg_h2, 1), helium_instead_t=round(hull.volume_m3 * PURITY * rho_he / 1e3, 0),
            cruise_power_mw=round(cruise / 1e6, 2),
            energy_per_passenger_km_wh=round(cruise / (n * CRUISE_M_S) / 3.6, 2) if n > 0 else None,
            moored_head_on_kn=round(q_design * design.drag_coefficient * hull.reference_area_m2 / 1e3, 0),
            moored_beam_on_mn=[round(q_service * BEAM_DRAG * side / 1e6, 2), round(q_design * BEAM_DRAG * side / 1e6, 1)],
            berths_for_long_haul=[round(f * BERTH_HOURS / (2.0 * n), 1) for f in long_haul_per_hour] if n > 0 else None))
    return rows


# ---------------------------------------------------------------------------------------------------------------
# Calibration: the rigid airships' weight statements

# Goodyear Aerospace 1975 (NASA CR-137692, vol. III, table 4, from the German weight reports and the Macon's final
# weight statement), in lb. Gross lift is 95 per cent inflation of the gas volume with hydrogen at 0.068 or helium at
# 0.062 lb/ft3. Longitudinals are the first line of the hull's breakdown; the cells' group holds the cells, nettings
# and valves; the power group the power plant, the fuel and oil system and the water recovery. Sizes, top speeds and
# gas volumes from the same report's table 3 and the ships' records (sources.json).
UNIT_LIFT_KG_M3 = dict(hydrogen=0.068 * LB / FT ** 3, helium=0.062 * LB / FT ** 3)
SHIPS = {
    'LZ 129 Hindenburg': dict(gas='hydrogen', length_m=245.0, diameter_m=41.2, gas_ft3=7.06e6, top_speed_kn=72.0,
                              cells=16, gross_lb=456_076, empty_lb=249_000, hull_lb=110_435, longitudinals_lb=28_129,
                              empennage_lb=16_669, cells_lb=25_524 + 1_043 + 2_605, cover_lb=18_753,
                              power_lb=36_465 + 5_625 + 8_336),
    'ZRS-5 Macon': dict(gas='helium', length_m=239.3, diameter_m=132.9 * FT, gas_ft3=6.85e6, top_speed_kn=72.2,
                        cells=12, gross_lb=403_465, empty_lb=236_497, hull_lb=94_320, longitudinals_lb=22_867,
                        empennage_lb=14_116, cells_lb=21_769 + 570 + 3_102, cover_lb=12_606,
                        power_lb=49_759 + 5_546 + 12_922),
    'LZ 127 Graf Zeppelin': dict(gas='hydrogen', length_m=236.6, diameter_m=30.5, gas_ft3=3.9e6, top_speed_kn=128 / 1.852,
                                 cells=16, gross_lb=260_300, empty_lb=137_000, hull_lb=61_117, longitudinals_lb=19_873,
                                 empennage_lb=6_998, cells_lb=12_822 + 715 + 1_047, cover_lb=10_836,
                                 power_lb=21_102 + 4_275 + 5_876),
    'ZR-3 Los Angeles': dict(gas='helium', length_m=200.0, diameter_m=27.64, gas_ft3=2.6e6, top_speed_kn=65.2,
                             cells=14, gross_lb=153_140, empty_lb=89_239, hull_lb=31_299, longitudinals_lb=11_298,
                             empennage_lb=2_750, cells_lb=8_769 + 671 + 781, cover_lb=7_401,
                             power_lb=21_199 + 2_651 + 4_000),
}
REFERENCE_SHIP = 'LZ 129 Hindenburg'
DESIGN_1930 = dict(frame=bu.DURALUMIN, gust_m_s=10.7, drag_coefficient=0.025)
# A modern rigid hull, on the calibrated factors: a carbon-fibre truss (bu.CARBON), the cover and gas cells Goodyear
# proposed in 1975 (Dacron, Tedlar and adhesive at 3.35 oz/yd2; scrim film cells at 2.0 oz/yd2, with the Hindenburg's
# share for nettings and valves), and an electric drive with fuel cells at an assumed 2 kg per kW.
MODERN = dict(frame=bu.CARBON, cover_kg_m2=3.35 * 0.033906, cell_kg_m2=2.0 * 0.033906 * 29_172 / 25_524,
              power_plant_kg_w=2.0e-3)
SPEED_M_S = 36.0                                   # top speed of the ships sized here, near the great rigids' 37 m/s
# Hull laminates: Uretek 3216LV (Vectran), 200 g/m2 and about 995 N/cm; Z2929T-AB (Zylon), 157 g/m2 and 997 N/cm.
FABRICS = (bu.Fabric("today's strongest hull laminates (Vectran, Zylon)", 0.18, 99.6e3),
           bu.Fabric('a laminate twice as strong', 0.18, 199.2e3))


def ship_case(ship: dict):
    """A historical ship as the model sees it: its hull (volume set to its gas volume), the air at low altitude and
    the lift per cubic metre of its gas volume at 95 per cent inflation."""
    v = ship['gas_ft3'] * FT ** 3
    hull = bu.Hull(ship['length_m'], ship['length_m'] / ship['diameter_m'],
                   volume_coefficient=v / (ship['length_m'] * ship['diameter_m'] ** 2))
    return hull, EARTH_SEA, 0.95 * UNIT_LIFT_KG_M3[ship['gas']]


def fit(columns, target):
    """Least squares in relative error: the coefficients c minimising sum((columns @ c / target - 1)^2)."""
    a = np.asarray(columns, dtype=float) / np.asarray(target, dtype=float)[:, None]
    c, *_ = np.linalg.lstsq(a, np.ones(len(target)), rcond=None)
    return c


def calibrate(ships=SHIPS) -> tuple:
    """Fit each group of the rigid-hull model to the ships' weight statements, and return the model's prediction for
    each ship beside its statement."""
    rows = []
    for name, ship in ships.items():
        hull, air, lift = ship_case(ship)
        d = bu.Design(speed_m_s=ship['top_speed_kn'] / 1.943844, cells=ship['cells'], **DESIGN_1930)
        i = bu.ideal_structure(hull, air, lift, d)
        rows.append(dict(name=name, ship=ship, bending=i['gust'] + i['static'], rings=i['rings'], fins=i['fins'],
                         cover_area=hull.area_m2 + 2.0 * i['fin_area_m2'],
                         cell_area=hull.area_m2 + 2.0 * ship['cells'] * np.pi * hull.radius_m ** 2 * 0.75,
                         power=float(bu.power_w(hull, air, d.speed_m_s, d))))
    kg = lambda r, key: r['ship'][key] * LB
    col = lambda f: [f(r) for r in rows]
    (girder,) = fit(np.c_[col(lambda r: r['bending'])], col(lambda r: kg(r, 'longitudinals_lb')))
    frame, frame_gas = fit(np.c_[col(lambda r: r['bending']), col(lambda r: r['rings'])],
                           col(lambda r: kg(r, 'hull_lb') - kg(r, 'longitudinals_lb')))
    (fin,) = fit(np.c_[col(lambda r: r['fins'])], col(lambda r: kg(r, 'empennage_lb')))
    (cover,) = fit(np.c_[col(lambda r: r['cover_area'])], col(lambda r: kg(r, 'cover_lb')))
    (cells,) = fit(np.c_[col(lambda r: r['cell_area'])], col(lambda r: kg(r, 'cells_lb')))
    (power,) = fit(np.c_[col(lambda r: r['power'])], col(lambda r: kg(r, 'power_lb')))
    fixed = float(np.mean(col(lambda r: (kg(r, 'empty_lb') - sum(kg(r, k) for k in (
        'hull_lb', 'empennage_lb', 'cells_lb', 'cover_lb', 'power_lb'))) / kg(r, 'gross_lb'))))
    factors = dict(girder_factor=float(girder), frame_factor=float(frame), frame_gas_factor=float(frame_gas),
                   fin_factor=float(fin), cover_kg_m2=float(cover), cell_kg_m2=float(cells),
                   power_plant_kg_w=float(power), fixed_share=fixed)
    checks = {}
    t = lambda x: round(float(x) / 1e3, 1)
    for r in rows:
        ship = r['ship']
        hull, air, lift = ship_case(ship)
        d = bu.Design(speed_m_s=ship['top_speed_kn'] / 1.943844, cells=ship['cells'], **DESIGN_1930, **factors)
        m = bu.rigid(hull, air, lift, d)
        parts = m['parts_kg']
        checks[r['name']] = dict(
            gross_lift_t=[t(kg(r, 'gross_lb')), t(m['gross_lift_kg'])],
            longitudinals_t=[t(kg(r, 'longitudinals_lb')), t(parts['longitudinals'])],
            frames_t=[t(kg(r, 'hull_lb') - kg(r, 'longitudinals_lb')), t(parts['frames'])],
            empennage_t=[t(kg(r, 'empennage_lb')), t(parts['empennage'])],
            cells_t=[t(kg(r, 'cells_lb')), t(parts['cells'])], cover_t=[t(kg(r, 'cover_lb')), t(parts['cover'])],
            power_t=[t(kg(r, 'power_lb')), t(parts['power_plant'])],
            empty_t=[t(kg(r, 'empty_lb')), t(m['empty_kg'])],
            empty_share=[round(ship['empty_lb'] / ship['gross_lb'], 3), round(float(m['empty_kg'] / m['gross_lift_kg']), 3)],
            weight_driven_share=round(float(m['weight_driven_share']), 3),
            aerodynamic_share=round(float(m['aerodynamic_share']), 3), areal_share=round(float(m['areal_share']), 3))
    return factors, checks


# ---------------------------------------------------------------------------------------------------------------
# The rule of similarity

def equal_density_heights(winds: dict) -> list:
    """Heights on the Moon whose air has the density of Earth's standard atmosphere at a height."""
    z = np.array(winds['height_km'])
    rho = np.array(winds['density_mean_kg_m3'], dtype=float)
    ok = np.isfinite(rho)
    rows = []
    for ze in (0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0):
        d = float(earth_air(ze * 1e3))
        rows.append(dict(earth_km=ze, density_kg_m3=round(d, 3),
                         moon_km=round(float(np.interp(d, rho[ok][::-1], z[ok][::-1])), 1)))
    return rows


# Earth's largest flyers, built or designed (sources.json): the length that sets their size and the mass that
# measures what they carry.
EARTH_FLYERS = (
    dict(name='LZ 129 Hindenburg (flew 1936)', size='length', size_m=245.0, mass='gross lift', mass_t=206.9),
    dict(name='CargoLifter CL160 (designed 2000)', size='length', size_m=260.0, mass='cargo', mass_t=160.0),
    dict(name='SkyCat 1000 hybrid (designed)', size='length', size_m=307.0, mass='payload', mass_t=1000.0),
    dict(name='Antonov An-225 (flew 1988)', size='span', size_m=88.4, mass='take-off mass', mass_t=640.0),
    dict(name='Stratolaunch Roc (flew 2019)', size='span', size_m=117.3, mass='take-off mass', mass_t=590.0),
    dict(name='Lockheed CL-1201 (studied about 1970)', size='span', size_m=340.0, mass='gross mass', mass_t=5375.0),
    dict(name='Mil Mi-26 (flew 1977)', size='rotor', size_m=32.0, mass='take-off mass', mass_t=56.0),
)


def similarity(winds: dict) -> dict:
    sea = at(winds, 'density_mean_kg_m3', 0.0)
    ratio = G_MOON / G_EARTH
    return dict(
        length_factor=round(SIMILAR, 3), mass_factor=round(SIMILAR ** 3, 1),
        earth_flyers=[dict(f, moon_size_m=round(f['size_m'] * SIMILAR, -1),
                           moon_mass_t=float(f"{f['mass_t'] * SIMILAR ** 3:.3g}")) for f in EARTH_FLYERS],
        same_flyer_same_air=dict(weight=round(ratio, 4), hover_power=round(ratio ** 1.5, 4),
                                 energy_per_km=round(ratio, 4), speed_at_a_lift_coefficient=round(ratio ** 0.5, 4)),
        same_flyer_design_air_against_earth_sea_level=dict(
            density_ratio=round(sea / 1.225, 3), hover_power=round(ratio ** 1.5 * (1.225 / sea) ** 0.5, 4),
            speed_at_a_lift_coefficient=round((ratio * 1.225 / sea) ** 0.5, 4)),
        equal_density_heights=equal_density_heights(winds))


# ---------------------------------------------------------------------------------------------------------------
# Results

TABLE_LENGTHS_M = (100.0, 245.0, 500.0, 1000.0, 1500.0, 2000.0, 3000.0, 5000.0, 10000.0, 20000.0)
OPERATION_LENGTHS_M = (245.0, 500.0, 750.0, 1000.0, 1500.0, 2000.0, 3000.0)
TURN_RADIUS_LENGTHS = 4.6          # the Los Angeles turned on a radius of about 3,000 ft, 4.6 lengths (NACA TR-333)


def summarise_scan(r: dict, reference: float) -> dict:
    L = r['length_m']
    rows = []
    for x in TABLE_LENGTHS_M:
        i = int(np.argmin(np.abs(np.log(L / x))))
        rows.append(dict(length_m=x, gross_lift_t=float(f"{r['gross_lift_kg'][i] / 1e3:.3g}"),
                         useful_share=round(float(r['useful_share'][i]), 3),
                         weight_driven_share=round(float(r['weight_driven_share'][i]), 3),
                         aerodynamic_share=round(float(r['aerodynamic_share'][i]), 3),
                         areal_share=round(float(r['areal_share'][i]), 3)))
    lm = {k: (None if v is None else round(v, 3 if k == 'best_useful_share' else -1)) for k, v in
          landmarks(r, reference).items()}
    curves = dict(length_m=[round(float(x), 1) for x in L], useful_share=[round(float(x), 4) for x in r['useful_share']])
    return dict(lift_per_m3=round(r['lift_per_m3'], 3), table=rows, landmarks=lm, **curves)


def results() -> dict:
    winds, ring, port = load()
    band_air, bw = moon_air(winds, BAND_KM), band_winds(winds, ring)
    band_top = moon_air(winds, FLIGHT_BAND_KM[1])
    factors, checks = calibrate()
    reference = 1.0 - SHIPS[REFERENCE_SHIP]['empty_lb'] / SHIPS[REFERENCE_SHIP]['gross_lb']
    moon_gust = bw['updraft_max_m_s']
    base = dict(DESIGN_1930, **factors, speed_m_s=SPEED_M_S)
    d1930, d1930_moon = bu.Design(**base), bu.Design(**{**base, 'gust_m_s': moon_gust})
    modern_earth = bu.Design(**{**base, **MODERN})
    modern_moon = bu.Design(**{**base, **MODERN, 'gust_m_s': moon_gust})
    same_air = bu.Air(G_MOON, 1.225, 288.15, EARTH_MOLAR_MASS)
    earth_lift = 0.95 * UNIT_LIFT_KG_M3['hydrogen']                      # the weight statements' basis
    band_lift = float(bu.net_lift_kg_m3(band_top, bu.HYDROGEN, PURITY))  # cells full at the band's top
    cases = dict(earth_1930s=(EARTH_SEA, d1930, earth_lift),
                 moon_1930s_in_earth_sea_level_air=(same_air, d1930, earth_lift),
                 moon_1930s_flight_band=(band_air, d1930_moon, band_lift),
                 moon_modern_flight_band=(band_air, modern_moon, band_lift),
                 earth_modern_sea_level=(EARTH_SEA, modern_earth, earth_lift))
    buoyant = dict(factors={k: round(v, 5) for k, v in factors.items()}, weight_statements=checks,
                   hindenburg_useful_share=round(reference, 3),
                   cases={name: summarise_scan(scan(air, d, lift), reference) for name, (air, d, lift) in cases.items()})
    pph = port['port']['occupancy']['passengers_per_hour']
    trips = port['port']['occupancy']['travel_by_trip']
    long_haul = [x * trips['long_haul_km2'] / sum(trips.values()) for x in pph]
    ops = operations(band_air, modern_moon, OPERATION_LENGTHS_M, bw['design_wind_m_s'], bw['service_wind_m_s'], long_haul,
                     lift_air=band_top)
    ops_1930 = operations(band_air, d1930_moon, OPERATION_LENGTHS_M, bw['design_wind_m_s'], bw['service_wind_m_s'],
                          long_haul, lift_air=band_top)
    for row in ops:
        row['turn_radius_km'] = round(TURN_RADIUS_LENGTHS * row['length_m'] / 1e3, 1)
        row['half_turn_minutes'] = round(np.pi * TURN_RADIUS_LENGTHS * row['length_m'] / CRUISE_M_S / 60.0, 1)
    pressure = {f.name: {str(v): dict(flight_band=largest_pressure_hull(band_air, f, v, moon_gust),
                                      moon_sea_level=largest_pressure_hull(moon_air(winds, 0.0), f, v, moon_gust),
                                      earth_sea_level=largest_pressure_hull(EARTH_SEA, f, v, 10.7))
                         for v in (25.0, 30.0, 36.0)} for f in FABRICS}
    return dict(
        schema=SCHEMA,
        producer=dict(study='sky_ships', files={str(Path(f).relative_to(ROOT)): digest(f) for f in PRODUCER_FILES},
                      constants=constants_used(PRODUCER_FILES),
                      inputs={str(Path(f).relative_to(ROOT)): digest(f) for f in (GLOBAL_WINDS, RING, PORT)}),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        air=air_table(winds, ring), band=bw, storms=storms(ring), similarity=similarity(winds), buoyant=buoyant,
        pressure_hulls=pressure, winged=winged_results(band_air.density_kg_m3, WING_GUST_BAND_M_S),
        hybrid_hull_lift=dict(flight_band=hull_lift_ratio(band_air, float(bu.net_lift_kg_m3(band_top, bu.HYDROGEN, PURITY)),
                                                          CRUISE_M_S, (40.0, 83.0, 167.0, 250.0)),
                              earth_sea_level=hull_lift_ratio(EARTH_SEA, float(bu.net_lift_kg_m3(EARTH_SEA, bu.HYDROGEN, PURITY)),
                                                              CRUISE_M_S, (40.0, 83.0, 167.0, 250.0))),
        operations=dict(
            passenger_kg=PASSENGER_KG, cruise_m_s=CRUISE_M_S, trip_hours=TRIP_HOURS, berth_hours=BERTH_HOURS,
            long_haul_per_hour=[round(x, -2) for x in long_haul], ships=ops,
            ships_1930s_materials=[dict(length_m=r['length_m'], useful_share=r['useful_share'], payload_t=r['payload_t'],
                                        passengers=r['passengers']) for r in ops_1930]))


def main() -> int:
    r = results()
    out = HERE / 'results' / 'sky_ships.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(r, indent=1) + '\n')
    b = r['buoyant']
    print(f"factors {b['factors']}, LZ 129 useful share {b['hindenburg_useful_share']}")
    for name, c in b['weight_statements'].items():
        print(f"{name:22s} empty {c['empty_t']} t")
    for name, case in b['cases'].items():
        print(f"{name:36s} {case['landmarks']}")
    for row in r['operations']['ships']:
        print(f"{row['length_m']:6.0f} m  gross {row['gross_lift_t']:>9} t  useful {row['useful_share']}  "
              f"passengers {row['passengers']}  hydrogen {row['hydrogen_t']} t")
    print(out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
