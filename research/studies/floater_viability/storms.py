"""Storms, lightning and hydrogen fire for floaters.

Pure functions plus evaluate(). Updrafts, rain and storm counts come from the CM1
equatorial ring; lightning rates from the electrified box and its Moon-wide
translation; the air from the GCM wind product. Ignition and spread figures
are requirements screens with stated inputs.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_GM, MOON_RADIUS

ROOT = Path(__file__).resolve().parents[3]
INPUT_FILES = (
    'climate/results/crm/ring_ring_equator.json',
    'climate/results/gcm/global_winds_A28_dim5_moon.json',
    'research/studies/sky_ships/results/sky_ships.json',
    'research/studies/floater_viability/storms_sources.json',
)
# Lightning on the Moon (research/studies/atmospheric_electricity, joint_synthesis README section 3):
# Moon-mean 0.0043-0.0079 flashes per km² a year, summit port cells 0.011-0.083, polar plain 0-0.0001;
# 106 of the electrified box's 530 flashes (20%) struck the ground; none flashed in the lunar night;
# in-cloud flashes started at a median 35.8 km. These are the products' own numbers, carried as inputs.
FLASH_DENSITY_PER_KM2_YEAR = dict(moon_mean=(0.0043, 0.0079), summit_cells=(0.011, 0.083), polar_plain=(0., 0.0001))
GROUND_STRIKE_SHARE = 106/530
IN_CLOUD_START_MEDIAN_KM = 35.8
H2_LHV_J_KG = 120e6
H2_MOLAR_KG = 0.002016
AIR_MOLAR_KG = 0.02896
GAS_KAPPA = 2/7  # R/c_p of a diatomic ideal gas: the 98% hydrogen lifting mixture


def _positive(**values):
    if any(not math.isfinite(v) or v <= 0 for v in values.values()):
        raise ValueError(f'finite positive inputs required: {tuple(values)}')


def _nonnegative(**values):
    if any(not math.isfinite(v) or v < 0 for v in values.values()):
        raise ValueError(f'finite nonnegative inputs required: {tuple(values)}')


def _fraction(name, value):
    if not math.isfinite(value) or not 0 < value <= 1:
        raise ValueError(f'{name} must lie in (0, 1]')


def _interp(x, xs, ys):
    if not xs[0] <= x <= xs[-1]:
        raise ValueError('outside the stored table')
    for i in range(1, len(xs)):
        if x <= xs[i]:
            f = (x-xs[i-1])/(xs[i]-xs[i-1])
            return ys[i-1] + f*(ys[i]-ys[i-1])
    return ys[-1]


# --- updrafts, downdrafts and rain ------------------------------------------------

def rise_before_full(gas_fraction, start_height_km, heights_km, pressure_pa, temperature_k):
    """Height a soft gas body reaches when lifted until its gas fills the envelope.

    Gas volume grows as T/p at fixed moles; the body is full where
    (T/T0)(p0/p) = 1/f. Below that height a soft body rides the air with no
    restoring force; above it, it must vent or hold superpressure.
    """
    _fraction('gas fraction', gas_fraction)
    p0 = _interp(start_height_km, heights_km, pressure_pa)
    t0 = _interp(start_height_km, heights_km, temperature_k)
    target = 1/gas_fraction
    for i in range(len(heights_km)-1):
        lo, hi = heights_km[i], heights_km[i+1]
        if hi <= start_height_km:
            continue
        lo = max(lo, start_height_km)
        e_lo = _interp(lo, heights_km, temperature_k)/t0*p0/_interp(lo, heights_km, pressure_pa)
        e_hi = _interp(hi, heights_km, temperature_k)/t0*p0/_interp(hi, heights_km, pressure_pa)
        if e_lo <= target <= e_hi:
            return lo + (target-e_lo)/(e_hi-e_lo)*(hi-lo)
    return None


def _log_pressure_height(target_pa, heights_km, pressure_pa):
    """Height where the profile's pressure equals target, interpolating ln p."""
    lp = [math.log(x) for x in pressure_pa]
    t = math.log(target_pa)
    for i in range(1, len(lp)):
        if lp[i-1] >= t >= lp[i]:
            return heights_km[i-1] + (lp[i-1]-t)/(lp[i-1]-lp[i])*(heights_km[i]-heights_km[i-1])
    return None


def _log_pressure_at(height_km, heights_km, pressure_pa):
    return math.exp(_interp(height_km, heights_km, [math.log(x) for x in pressure_pa]))


def rise_before_full_adiabatic(gas_fraction, start_height_km, heights_km, pressure_pa, kappa=GAS_KAPPA):
    """Height where a soft body's gas, cooling on its own adiabat as it rises, fills the envelope.

    Volume grows as (p0/p)^(1-kappa) when the gas exchanges no heat with the air;
    the body is full at p = p0 f^(1/(1-kappa)). Gas that keeps the air's
    temperature fills lower (rise_before_full); real bodies lie between.
    """
    _fraction('gas fraction', gas_fraction)
    p0 = _log_pressure_at(start_height_km, heights_km, pressure_pa)
    return _log_pressure_height(p0*gas_fraction**(1/(1-kappa)), heights_km, pressure_pa)


def vent_share_profile(start_height_km, rise_km, heights_km, pressure_pa, temperature_k, kappa=None):
    """Share of gas a full body of fixed volume vents over a forced rise, from the air's profile.

    Gas at the air's temperature keeps the air's density ratio: 1 - rho(z1)/rho(z0).
    With kappa, the gas cools on its own adiabat: 1 - (p1/p0)^(1-kappa).
    """
    _nonnegative(rise=rise_km)
    z1 = start_height_km + rise_km
    p0 = _log_pressure_at(start_height_km, heights_km, pressure_pa)
    p1 = _log_pressure_at(z1, heights_km, pressure_pa)
    if kappa is not None:
        return 1 - (p1/p0)**(1-kappa)
    t0 = _interp(start_height_km, heights_km, temperature_k)
    t1 = _interp(z1, heights_km, temperature_k)
    return 1 - (p1/p0)*(t0/t1)


def hold_against_updraft(updraft_m_s, air_density_kg_m3, gravity_m_s2, drag_coefficient=.47):
    """Extra weight per projected m² that holds a body at its height in an updraft."""
    _nonnegative(updraft=updraft_m_s)
    _positive(density=air_density_kg_m3, gravity=gravity_m_s2, cd=drag_coefficient)
    return .5*air_density_kg_m3*drag_coefficient*updraft_m_s**2/gravity_m_s2


def rain_brake_minutes(updraft_m_s, rain_mm_h, air_density_kg_m3, gravity_m_s2, drag_coefficient=.47,
                       catch_share=1.):
    """Minutes of kept storm rain that weigh a body down enough to hold against an updraft."""
    _nonnegative(rain=rain_mm_h)
    _fraction('catch share', catch_share)
    need = hold_against_updraft(updraft_m_s, air_density_kg_m3, gravity_m_s2, drag_coefficient)
    return None if rain_mm_h == 0 else 60*need/(rain_mm_h*catch_share)


def vent_share_to_hold(rise_km, scale_height_km):
    """Share of gas a full body vents to stay at constant pressure over a forced rise."""
    _nonnegative(rise=rise_km)
    _positive(scale=scale_height_km)
    return 1-math.exp(-rise_km/scale_height_km)


def rain_surge(rain_mm_h, minutes, drainage_kg_m2_h=0., catch_share=1.):
    """Water a body holds after a burst if it drains at a set rate (1 mm = 1 kg/m²)."""
    _nonnegative(rain=rain_mm_h, minutes=minutes, drainage=drainage_kg_m2_h)
    _fraction('catch share', catch_share)
    return max(0., rain_mm_h*catch_share-drainage_kg_m2_h)*minutes/60


def residual_descent(excess_kg_m2, air_density_kg_m3, gravity_m_s2, drag_coefficient=.47):
    """Drag-balanced sink speed of a body left heavy by excess_kg_m2."""
    _nonnegative(excess=excess_kg_m2)
    _positive(density=air_density_kg_m3, gravity=gravity_m_s2, cd=drag_coefficient)
    return math.sqrt(2*excess_kg_m2*gravity_m_s2/(air_density_kg_m3*drag_coefficient))


def storm_cover(storms_per_snapshot, storm_width_km, ring_length_km):
    """Share of a ring under storm rain at any moment, from counts and widths."""
    _nonnegative(count=storms_per_snapshot)
    _positive(width=storm_width_km, ring=ring_length_km)
    return min(1., storms_per_snapshot*storm_width_km/ring_length_km)


# --- hydrogen mixing and fire -----------------------------------------------------

def mole_to_mass_fraction(mole_fraction, molar_gas=H2_MOLAR_KG, molar_air=AIR_MOLAR_KG):
    _fraction('mole fraction', mole_fraction)
    return mole_fraction*molar_gas/(mole_fraction*molar_gas+(1-mole_fraction)*molar_air)


def flammable_reach(hole_diameter_m, gas_density_kg_m3, air_density_kg_m3, exit_hydrogen_mass_fraction,
                    lower_limit_mole=.04, decay_slope=.208):
    """Distance along a round jet to the lower flammable limit of hydrogen.

    Centreline decay Y_exit/Y = C_Y x/D_ef with D_ef = d sqrt(rho_gas/rho_air)
    (C_Y = 0.208 for a hydrogen round jet; 0.31-0.33 for jets from an orifice in the
    side of a tube). The 4% mole limit is a hydrogen mass fraction near 0.003, so the
    flammable core reaches a few hundred hole diameters.
    """
    _positive(hole=hole_diameter_m, gas=gas_density_kg_m3, air=air_density_kg_m3, slope=decay_slope)
    _fraction('exit fraction', exit_hydrogen_mass_fraction)
    y = mole_to_mass_fraction(lower_limit_mole)
    reach = exit_hydrogen_mass_fraction*math.sqrt(gas_density_kg_m3/air_density_kg_m3)*hole_diameter_m/(decay_slope*y)
    return dict(hole_diameter_m=hole_diameter_m, lower_limit_mass_fraction=y,
                reach_m=reach, reach_hole_diameters=reach/hole_diameter_m)


def momentum_length(hole_diameter_m, overpressure_pa, gas_density_kg_m3, air_density_kg_m3, gravity_m_s2):
    """Distance over which a light-gas jet from a hole is driven by its momentum.

    L_m = M^(3/4)/B^(1/2) with the momentum flux M = (rho_gas/rho_air) u^2 A and the
    buoyancy flux B = g (rho_air - rho_gas)/rho_air u A, u = sqrt(2 dp/rho_gas).
    Beyond L_m the plume rises by buoyancy and dilutes faster than the jet law.
    """
    _positive(hole=hole_diameter_m, dp=overpressure_pa, gas=gas_density_kg_m3, air=air_density_kg_m3,
              gravity=gravity_m_s2)
    if gas_density_kg_m3 >= air_density_kg_m3:
        raise ValueError('a light gas is required')
    u = math.sqrt(2*overpressure_pa/gas_density_kg_m3)
    area = math.pi*hole_diameter_m**2/4
    m = gas_density_kg_m3/air_density_kg_m3*u*u*area
    b = gravity_m_s2*(air_density_kg_m3-gas_density_kg_m3)/air_density_kg_m3*u*area
    return dict(exit_speed_m_s=u, momentum_length_m=m**.75/b**.5)


def upper_limit_in_air(oxygen_mole_fraction_air, limiting_oxygen=.05):
    """Hydrogen fraction above which a hydrogen-air mixture holds too little oxygen to burn."""
    _fraction('oxygen', oxygen_mole_fraction_air)
    _fraction('limiting oxygen', limiting_oxygen)
    return 1-limiting_oxygen/oxygen_mole_fraction_air


def strike_rate(ground_strikes_per_km2_year, body_radius_m, striking_distance_m):
    """Direct strikes a year on a body that ground-strike channels pass by within reach."""
    _nonnegative(density=ground_strikes_per_km2_year, striking=striking_distance_m)
    _positive(radius=body_radius_m)
    area_km2 = math.pi*(body_radius_m+striking_distance_m)**2/1e6
    return dict(capture_area_km2=area_km2, strikes_per_year=ground_strikes_per_km2_year*area_km2)


def burn_energy(hydrogen_kg):
    _nonnegative(hydrogen=hydrogen_kg)
    return hydrogen_kg*H2_LHV_J_KG


def neighbour_dose(hydrogen_kg, burn_seconds, radiant_fraction, distance_m):
    """Radiant energy per m² reaching a neighbour from one burning compartment."""
    _positive(burn=burn_seconds, distance=distance_m)
    _nonnegative(hydrogen=hydrogen_kg, radiant=radiant_fraction)
    power = burn_energy(hydrogen_kg)/burn_seconds
    flux = radiant_fraction*power/(4*math.pi*distance_m**2)
    return dict(fire_power_w=power, flux_w_m2=flux, dose_j_m2=flux*burn_seconds)


def wet_skin_drying_energy(water_kg_m2, latent_j_kg=2.45e6):
    """Heat per m² that must dry a wet skin before it can burn."""
    _nonnegative(water=water_kg_m2)
    return water_kg_m2*latent_j_kg


def oxygen_ingress_years(gas_column_m, contact_area_ratio, oxygen_permeance_mol_m2_s_pa,
                         outside_oxygen_pa, gas_molar_density_mol_m3, limit_fraction=.05):
    """Years for oxygen leaking into an unscrubbed gas volume to reach a fraction.

    Linear ingress at the full outside partial pressure (a short-time upper rate).
    """
    _positive(column=gas_column_m, ratio=contact_area_ratio, permeance=oxygen_permeance_mol_m2_s_pa,
              pressure=outside_oxygen_pa, molar=gas_molar_density_mol_m3)
    _fraction('limit', limit_fraction)
    moles = gas_column_m*gas_molar_density_mol_m3
    flux = oxygen_permeance_mol_m2_s_pa*outside_oxygen_pa*contact_area_ratio
    return limit_fraction*moles/flux/(JULIAN_DAY*JULIAN_YEAR_DAYS)


def scrub_cost_share(oxygen_permeance, hydrogen_permeance, outside_oxygen_pa, inside_hydrogen_pa):
    """Extra hydrogen burned by oxygen-scrubbing microbes over the hydrogen lost by permeation.

    2 H2 + O2 -> 2 H2O consumes two hydrogen per oxygen that enters.
    """
    _positive(o2=oxygen_permeance, h2=hydrogen_permeance, po2=outside_oxygen_pa, ph2=inside_hydrogen_pa)
    return 2*oxygen_permeance*outside_oxygen_pa/(hydrogen_permeance*inside_hydrogen_pa)


def evaluate(root=ROOT, compartment_hydrogen_kg=(1753., 14409., 333840.)):
    """Updraft entrainment, rain brakes, flammable reach, strike rates and fire spread."""
    root = Path(root)
    ring = json.loads((root/INPUT_FILES[0]).read_text())
    winds = json.loads((root/INPUT_FILES[1]).read_text())
    ships = json.loads((root/INPUT_FILES[2]).read_text())
    for data, schema in ((ring, 'terluna.climate.crm-ring/1'), (winds, 'terluna.climate.gcm-global-winds/1'),
                         (ships, 'terluna.research.sky-ships/1')):
        if data.get('schema') != schema:
            raise ValueError('unexpected input schema')
    air = {a['height_km']: a for a in ships['air']}
    g10 = MOON_GM/(MOON_RADIUS+10000.)**2
    rho10 = air[10]['density_kg_m3']
    gas10 = rho10-air[10]['lift_hydrogen_kg_m3']
    h, p, t = winds['height_km'], winds['pressure_pa'], winds['temperature_k']
    entrain = [dict(gas_fraction=f, start_height_km=z0, full_at_km=rise_before_full(f, z0, h, p, t),
                    full_at_km_gas_on_its_adiabat=rise_before_full_adiabatic(f, z0, h, p))
               for f in (.3, .5, .64, .8, .95) for z0 in (10., 20.)]
    up = ring['vertical_wind_m_s']
    updrafts = [dict(height_km=int(k.split('_')[0]), abs_p99_m_s=v['abs_p99'], up_max_m_s=v['up_max'],
                     down_max_m_s=v['down_max']) for k, v in up.items()]
    storms = ring['storms']
    hold = [dict(updraft_m_s=w, drag_coefficient=cd, hold_kg_m2=hold_against_updraft(w, rho10, g10, cd),
                 rain_minutes_at_median_storm=rain_brake_minutes(w, storms['rain_max_mm_h']['median'], rho10, g10, cd),
                 rain_minutes_at_p90_storm=rain_brake_minutes(w, storms['rain_max_mm_h']['p90'], rho10, g10, cd))
            for w in (1., storms['updraft_max_m_s']['median'], storms['updraft_max_m_s']['p90'],
                      up['10_km']['up_max'], storms['updraft_max_m_s']['max']) for cd in (.47, 1.2)]
    venting = [dict(forced_rise_km=r, vent_share=vent_share_profile(10., r, h, p, t),
                    vent_share_gas_on_its_adiabat=vent_share_profile(10., r, h, p, t, GAS_KAPPA))
               for r in (1., 5., 10., 20.)]
    p10 = _log_pressure_at(10., h, p)
    dz = .5
    pressure_scale = 2*dz/math.log(_log_pressure_at(10.-dz, h, p)/_log_pressure_at(10.+dz, h, p))
    density_scale = 2*dz/math.log(_log_pressure_at(10.-dz, h, p)/_interp(10.-dz, h, t)
                                  / (_log_pressure_at(10.+dz, h, p)/_interp(10.+dz, h, t)))
    superpressure = dict(gas_holds_its_temperature_pa_per_km=rho10*g10*1000.,
                         gas_cools_with_the_air_pa_per_km=p10/density_scale)
    surges = [dict(rain_mm_h=r, minutes=m, drainage_kg_m2_h=d, held_kg_m2=rain_surge(r, m, d),
                   sink_m_s=residual_descent(rain_surge(r, m, d), rho10, g10))
              for r in (storms['rain_max_mm_h']['median'], storms['rain_max_mm_h']['p90'], storms['rain_max_mm_h']['max'])
              for m in (6., 60.) for d in (0., 20., 70.)]
    ring_length = 2*math.pi*MOON_RADIUS/1000
    cover = storm_cover(storms['per_snapshot'], storms['width_km']['median'], ring_length)
    exit_fraction = mole_to_mass_fraction(.98)  # hydrogen in the 98% lifting mixture (2% air)
    reach = []
    for d in (.001, .01, .1, 1.):
        lengths = [momentum_length(d, dp, gas10, rho10, g10)['momentum_length_m'] for dp in (100., 300.)]
        for k in (.208, .326):
            r = flammable_reach(d, gas10, rho10, exit_fraction, decay_slope=k)
            reach.append(dict(**r, decay_slope=k, momentum_length_m_at_100_300_pa=lengths,
                              jet_law_holds=r['reach_m'] <= max(lengths)))
    upper = dict(lunar_air_oxygen=.175, earth_air_oxygen=.2095,
                 lunar=upper_limit_in_air(.175), earth=upper_limit_in_air(.2095),
                 lifting_mixture_hydrogen=.98, lifting_mixture_oxygen=.02*.175)
    strikes = []
    for place, (lo, hi) in FLASH_DENSITY_PER_KM2_YEAR.items():
        for density in (lo, hi):
            ground = density*GROUND_STRIKE_SHARE
            for radius in (20., 100., 1000.):
                for reach_m in (50., 200.):
                    s = strike_rate(ground, radius, reach_m)
                    strikes.append(dict(place=place, flashes_per_km2_year=density,
                                        ground_strikes_per_km2_year=ground, body_radius_m=radius,
                                        striking_distance_m=reach_m, **s,
                                        strikes_in_300_years=300*s['strikes_per_year']))
    fires = [dict(compartment_hydrogen_kg=m, burn_seconds=b, radiant_fraction=x, distance_m=d,
                  energy_j=burn_energy(m), **neighbour_dose(m, b, x, d))
             for m in compartment_hydrogen_kg for b in (30., 120.) for x in (.03, .1, .16) for d in (40., 100.)]
    skins = [dict(water_kg_m2=w, drying_j_m2=wet_skin_drying_energy(w)) for w in (.1, .5, 2.)]
    barrer = 3.35e-16
    oxygen = []
    for ratio in (.1, .25, .5):
        # Gas column and gas-contact area (outer skin plus the gas-air partition) per projected m²:
        # a one-layer raft of 40-m spheres on its hexagonal cells, and 200-m and 2-km spheres.
        raft = math.pi/(2*math.sqrt(3))
        for column, contact in ((4/3*20*raft, 5*raft), (400/3, 5.), (4000/3, 5.)):
            perm = ratio*.173*barrer/100e-6
            oxygen.append(dict(o2_to_h2_permeability=ratio, gas_column_m=column, contact_area_ratio=contact,
                               years_to_5_percent=oxygen_ingress_years(column, contact, perm, .175*air[10]['pressure_atm']*101325.,
                                                                      air[10]['pressure_atm']*101325./(8.314462618*(air[10]['temperature_c']+273.15))),
                               scrub_hydrogen_share=scrub_cost_share(ratio, 1., .175, .98)))
    return dict(
        evidence='Storm statistics from the CM1 equatorial ring (one lunar day, 3-hourly), lightning rates from the electrified box and its Moon-wide translation, air from the GCM wind product. Striking distances, radiant fractions, burn times and oxygen permeabilities are stated ranges.',
        reading_rule='Per projected m² unless named. Entrainment heights are where a lifted soft body becomes full, with its gas at the air\'s temperature or cooling on its own adiabat. Venting shares hold a full body\'s volume over a forced rise from 10 km. Flammable reach follows the momentum-jet law; beyond the momentum length the plume rises by buoyancy and dilutes faster, so longer reaches are upper bounds. Strike rates count ground-strike channels passing within reach of a body of that radius. Fire doses are radiant J/m² at a neighbour.',
        updrafts=updrafts, storm_counts=dict(per_snapshot=storms['per_snapshot'], width_km=storms['width_km'],
            rain_max_mm_h=storms['rain_max_mm_h'], updraft_max_m_s=storms['updraft_max_m_s'],
            lifetime_h=ring['storm_tracks']['lifetime_h'], cloud_top_km=storms['cloud_top_km'],
            ring_cover_share=cover),
        entrainment=entrain, hold=hold,
        venting=dict(start_height_km=10., pressure_scale_height_km=pressure_scale,
                     density_scale_height_km=density_scale, cases=venting, superpressure=superpressure),
        rain_surges=surges, flammable_reach=reach, upper_limit=upper,
        lightning=dict(flash_density={k: list(v) for k, v in FLASH_DENSITY_PER_KM2_YEAR.items()},
                       ground_share=GROUND_STRIKE_SHARE,
                       in_cloud_start_median_km=IN_CLOUD_START_MEDIAN_KM, night_flashes=0,
                       strikes=strikes),
        fires=fires, wet_skin=skins, oxygen=oxygen,
        unresolved=['Triggered lightning on a wet floating body in a storm field.',
                    'Flame spread across a wet, compartmented skin, and burn times of real tears.',
                    'Storm extremes over centuries; updrafts sharper than 6-km columns.',
                    'Oxygen permeability of the candidate barrier, measured with its hydrogen permeability.'])
