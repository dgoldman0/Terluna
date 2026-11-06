"""The fault across the ring: what can stop every beam at once, which of those faults can also reach the shade, how
long a restart takes and what storage it needs, and how the beams' control can be kept apart from the tiles' attitude
and station-keeping (ring_fault.md).

A closed-form screen on committed products; it runs no dynamics. It reads:
- the power screen (provisioning/results/power.json): the critical loads and the week's stores, the towns' lifting
  hydrogen, the transmitting tiles and their radio power;
- the attitude schemes (results/attitude_schemes.json): the torques on a radial-facing tile with its balanced filter,
  the momentum stores, and the edge-on turn times, which fix a tile's pitch inertia;
- photon control (results/photon_control.json): the keeping demand of the kept run;
- the ring bundle (results/ring_bundle.json): when free rings first come within 150 m, and the edge ring's light; the
  kept run (results/ring_keeping.json): the keeping's cost per orbit;
- the tile design (results/tile_design.json): a failed tile's times to the 150 m gate and to contact, the stack's
  coverage multiplicity, the service fleet's dispatch rate, the torque allocation and the reference tile's inertia;
- the rolled layout (results/ring_layout_rolled.json): the radius at which each height's ring turns with the Sun, which
  orders the rings in radius at their common line of nodes;
- the swarm's UV transmission (protection/transmission/results/uv_transmission.json) and the loss response's base
  pressure (atmosphere/loss_response/results/loss_response.json).
Values the products give only in their notes (the latch's thresholds, the climate's e-folding, the energy study's local
supplies on domain/provisioning) and the design guesses are carried in DESIGN with their sources
(ring_fault_sources.json).

Sections of the product:
- shade: how fast a tile reaches the 150 m gate and contact once its attitude or keeping control is lost, and how many
  lost tiles open the shade past the latch's thresholds;
- nodes: the speeds at which radius neighbours pass at the stack's common line of nodes, and how often a fragment that
  sweeps the stack's radii meets a tile there;
- meteoroid_storm: what a day of a meteoroid storm holes, against the latch;
- swollen_loss: the cycle times of the losses the latch leads to;
- upper_air: an order-of-magnitude bracket on how fast the air above the base could warm toward the latch;
- recoil: the push a transmitting tile's beams give it, and the step a ring-wide beam stop puts on its keeping;
- trim: how much of the shade rides on the reflectivity trim;
- beams: how fast the critical loads need power back, the stores each fault class needs and the share local supplies
  carry;
- faults: the fault classes, what each reaches with and without the separation, and their restart times;
- separation: what keeping the beams' control apart costs;
- targets: how improbable a fault that opens the shade must be made.

    python -m research.studies.solar_shield_array.ring_fault
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np

from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_GM, MOON_SURFACE_GRAVITY, SPEED_OF_LIGHT
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE/'results/ring_fault.json'
SCHEMA = 'terluna.research.ring-fault/1'
FILES = ['research/studies/solar_shield_array/ring_fault.py']
POWER = ROOT/'provisioning/results/power.json'
ATTITUDE = HERE/'results/attitude_schemes.json'
PHOTON = HERE/'results/photon_control.json'
BUNDLE = HERE/'results/ring_bundle.json'
KEPT = HERE/'results/ring_keeping.json'
TILE = HERE/'results/tile_design.json'
LAYOUT = HERE/'results/ring_layout_rolled.json'
UV = ROOT/'protection/transmission/results/uv_transmission.json'
LOSS = ROOT/'atmosphere/loss_response/results/loss_response.json'
SOURCES = HERE/'ring_fault_sources.json'
INPUTS = (POWER, ATTITUDE, PHOTON, BUNDLE, KEPT, TILE, LAYOUT, UV, LOSS, SOURCES)
HOUR = 3600.
LABELS = ('published', 'recalled', 'product', 'computed', 'design guess')
FUNCTIONS = ('beams', 'attitude_keeping', 'film')

DESIGN = dict(
    # The tile and its shingle (tile_design.md; relative_orbits.md): 10 km squares at an 8.5 km pitch along the ring,
    # shingled 1.5 degrees, rings 1 km apart in radius, the 150 m safety gate and the tile study's 50 m envelope.
    tile_side_m=10000., pitch_m=8500., shingle_deg=1.5, clearance_m=150., envelopes_m=[0., 50.],
    # Areal masses whose inertia is screened: the annulus allocation (constraints review), the mean tile (46 Gt over
    # 28 million tiles, array_industry.md) and the kept run's tiles, against which the turn times are checked.
    areal_masses_kg_m2=dict(annulus_allocation=.013765, mean_tile=.0164, kept_run=None),
    mean_tile_kg=1.64e6, fleet_inventory_Gt=46.,
    design_radius_km=20000.,
    # The latch (product text: atmosphere/loss_response/README.md, 'The runaway is a latch'): the UV transmission
    # through gaps at quiet Sun that tips the upper air into its swollen state, and the heat at that threshold.
    latch_transmission=dict(titania=[.0035, .0082], edge_200_nm=[.0012, .0024]),
    latch_heat_W_m2=dict(titania=[3.9e-5, 8.4e-5], edge_200_nm=[1.9e-5, 3.0e-5]),
    # The upper air's warming (computed bracket, no product gives a response time): the air stays within a few kelvin
    # of the base temperature for three or four e-folds of pressure (product text), so the warmed share of the column
    # above the base is e^-4 to e^-3; nitrogen's heat capacity (recalled); the warming that would carry the exobase
    # to the threshold is an assumption.
    warmed_efolds=[3., 4.], cp_J_kg_K=1040., warming_K=[30., 300.],
    climate_efolding_years=[4.0, 7.1],      # product text: climate/gcm/README.md, 'How long a cold start takes'
    # Fragments from a collision at the nodes (computed geometric screen): the tilt of the fragment's plane, and the
    # relative speed counted as fast.
    fragment_tilts_deg=[0., 12.], fast_m_s=50.,
    # A meteoroid storm (recalled): a stream's flux ten to a thousand times above the sporadic background for a day,
    # holes up to four times the particle's diameter as the swarm's product takes them.
    storm_enhancement=[10., 1000.], storm_days=1., storm_hole_factor=4.,
    # The energy study's local supplies (product: domain/provisioning afe38eda,
    # provisioning/results/energy_options.json, sha256 127c599d87d044e6c923f54071cc71f19d0e04f6e5bd101da93bd19b1cdc42d9,
    # read-only and outside this tree), TW
    # delivered on average: rivers developable, tides, the natural heat flow if all taken, and roofs and pavements with
    # panels stored in hydrogen. The seas' thermal gradients (0.24-0.50 TW) wait on cold deep water and are left out.
    local_tw=dict(rivers=[.2218, .4437], tides=[.0288, .0577], heat_flow=[.05311, .1195], roofs=[.8273, 6.154]),
    # How fast the critical loads need their power back, seconds (published where named, else design guesses).
    tolerance_s=dict(
        close_computing=dict(range=[.01, 1.], kind='design guess', basis='millisecond control rides on its own '
                             'uninterruptible supply'),
        health_care=dict(range=[10., 10.], kind='published', basis='NFPA 110 Level 1, Type 10: emergency power within '
                         '10 s for life safety; Class 96, 96 hours of fuel on site'),
        sky_control=dict(range=[1., 60.], kind='design guess', basis="the towns' flight control rides on their own "
                         'batteries and fuel cells'),
        food_cultures=dict(range=[60., 3600.], kind='design guess', basis='dense aerobic cultures use their dissolved '
                           'oxygen within minutes without aeration (recalled); cooling holds them for about an hour'),
        refrigeration=dict(range=[4*HOUR, 48*HOUR], kind='published', basis='USDA: a closed refrigerator keeps food '
                           'safe 4 hours, a full freezer 48'),
        pumps=dict(range=[2*HOUR, 24*HOUR], kind='design guess', basis='water and sewage pumps; tanks and reservoirs '
                   'buffer them for hours'),
        water_lift=dict(range=[1*JULIAN_DAY, 3*JULIAN_DAY], kind='design guess', basis='the lift fills reservoirs'),
        lamp_farms=dict(range=[2*JULIAN_DAY, 14*JULIAN_DAY], kind='design guess', basis='crops kept dark with air and '
                        "water moving; the field crops already pass the 354-hour lunar night"),
    ),
    # Allowed probability of a shade-opening fault over a horizon (design guess), and the comparison standards.
    accepted_probability=.01, horizons_years=[1e3, 1e4, 1e6, 1e9],
    # The Halloween storms of 2003 (published, NOAA service assessment): 28 satellites reported anomalies and 2 were put
    # out of service; the satellites then operating, about 800-1,000, are recalled.
    halloween_2003=dict(anomalies=28, out_of_service=2, operating=[800, 1000]),
    # A whole grid's black start (published, Iberia 28 April 2025): Portugal's load back in about 12 hours, Spain's in
    # about 24, with their sources intact. Added to every restart that stops the beams.
    grid_restoration_days=[.5, 1.],
    # The atmosphere's mass (product text: atmosphere/loss_response/README.md, the feasibility baseline's column) and
    # the loss just short of the latch behind the titania stack with collisional upper air (same README), for the
    # swollen state's cycle times; the unshielded energy-limited losses come from the loss product.
    atmosphere_kg=3.1e18, loss_short_of_onset_kg_s=[140., 190.],
    aviation_catastrophic_per_flight_hour=1e-9,      # published: FAA AC 25.1309-1
    nuclear_large_release_per_reactor_year=1e-6,     # published: US NRC safety goal
    collateral_tiles=[2, 4],        # neighbours an uncontrolled tile damages before it is contained (design guess)
    margin=10.,                     # on the correlated-loss budget (design guess)
    # The separation's added hardware per tile (design guesses): an independent controller, sensors, a few square
    # metres of photovoltaics and a battery for attitude and keeping alone.
    separate_controller_kg=[10., 100.], separate_controller_W=[10., 100.],
    update_soak_days=[1., 3.], fallback_failure=[1e-3, 1e-4],
    # The trim (product text): devices over 25-89% of each tile change the redirected band's reflection by up to a
    # quarter (array_industry.md); the window stack reflects 4.69% of the sunlight and the annulus film 5.79%
    # (constraints_review.md); each extra percent of sunlight the films absorb adds 1.7-2.9 W/m2 of glow (register,
    # 9 October, from the provisioning screen).
    trim=dict(range=.25, reflected=[.0469, .0579], cover=[.25, .89], glow_W_m2_per_percent=[1.7, 2.9],
              basis='array_industry.md (the shield\'s own loads); constraints_review.md finding 1; decisions.md, '
                    '9 October'),
    # The 1978 reference's transmitter, 10 kg per kW of electricity (product text: energy_options.md), for the
    # transmitting tiles' areal mass; and a thin-film sandwich (design guess).
    transmitter_kg_per_kw=10., transmitter_efficiency=.7, thin_film_transmitter_kg_m2=.2,
)

# The fault classes (proposals and design guesses; precedents in ring_fault_sources.json). 'reach' lists the functions
# a fault takes down directly, without and with the separation; a fault opens the shade where it reaches the film or
# the tiles' attitude and keeping, whose loss brings tiles into contact within hours. 'modes' says how each function
# fails without the separation ('stops', 'holds', 'drifts' or 'unaffected'); 'restart_days' is the time to restore the
# beams; 'spread' and 'detection' are in words.
FAULTS = [
    dict(id='F1', name="the beams' phase and timing reference",
         reach=dict(without=['beams'], with_=['beams']),
         modes=dict(beams='stops', attitude_keeping='unaffected', film='unaffected'),
         spread='at the speed of the reference: every beam that takes it loses coherence within seconds',
         detection='each transmitter sees its pilot lock and its spot efficiency fall',
         restart_days=dict(without=[.1, 7.], with_=[.01, .25]),
         likelihood='occasional: about one constellation-wide timing or reference fault a decade per system',
         precedents=['galileo_2019', 'gps_2016'], kind='design guess'),
    dict(id='F2', name="the receivers' pilot signals",
         reach=dict(without=['beams'], with_=['beams']),
         modes=dict(beams='stops', attitude_keeping='unaffected', film='unaffected'),
         spread='as fast as a common fault in the receivers takes effect; a lost pilot defocuses its beam at once',
         detection='beams douse on their own; the receivers report',
         restart_days=dict(without=[.1, 3.], with_=[.1, 3.]),
         likelihood='occasional: a shared design or update across 0.7 million crowns and up to 48,000 patches',
         precedents=['dickinson_1978', 'doe_nasa_1978'], kind='design guess'),
    dict(id='F3', name='scheduling and command software',
         reach=dict(without=['beams', 'attitude_keeping'], with_=['beams']),
         modes=dict(beams='holds', attitude_keeping='holds', film='unaffected'),
         spread='beams hold their last assignment and fade as their tiles set, over 2-20 hours',
         detection='missed handovers; watchdogs on the schedulers',
         restart_days=dict(without=[.01, 1.], with_=[.01, .5]),
         likelihood='expected: global software services fail for hours several times a decade',
         precedents=['starlink_2025', 'crowdstrike_2024'], kind='design guess'),
    dict(id='F4', name='a bad update, or a latent timer, across the fleet',
         reach=dict(without=['beams', 'attitude_keeping'], with_=['beams']),
         modes=dict(beams='stops', attitude_keeping='drifts', film='unaffected'),
         spread='as fast as the rollout; a latent timer fires everywhere at once',
         detection='canary waves; the second, dissimilar stack disagrees',
         restart_days=dict(without=[.5, 10.], with_=[.1, 3.]),
         likelihood='expected over the fleet life unless no update ever reaches the whole fleet at once',
         precedents=['crowdstrike_2024', 'glonass_2014'], kind='design guess'),
    dict(id='F5', name='an intrusion into command',
         reach=dict(without=['beams', 'attitude_keeping'], with_=['beams']),
         modes=dict(beams='stops', attitude_keeping='drifts', film='unaffected'),
         spread='at the attacker\'s pace once trusted management access is held',
         detection='command authentication, rate limits that refuse fleet-wide moves, anomaly watch',
         restart_days=dict(without=[3., 60.], with_=[1., 14.]),
         likelihood='occasional to expected: constellation-scale intrusions have happened',
         precedents=['viasat_2022'], kind='design guess'),
    dict(id='F6', name='space weather across the fleet',
         reach=dict(without=['beams', 'attitude_keeping'], with_=['beams', 'attitude_keeping']),
         modes=dict(beams='stops', attitude_keeping='holds', film='unaffected'),
         spread='simultaneous: solar protons bring upsets and latch-ups to every tile within hours of a flare and '
                'stay for days; surface and internal charging in Earth\'s magnetotail, which the Moon crosses for a '
                'few days each month (recalled)',
         detection='radiation monitors on the tiles and upstream',
         restart_days=dict(without=[.5, 7.], with_=[.5, 7.]),
         likelihood='expected: Oct-1989-class proton events a few times a solar cycle; the ring sits outside '
                    "Earth's magnetosphere most of each month",
         precedents=['tdrs_1989', 'halloween_2003', 'anik_1994', 'galaxy15_2010'], kind='design guess'),
    dict(id='F7', name='shared navigation and timing for attitude and keeping',
         reach=dict(without=['beams', 'attitude_keeping'], with_=['attitude_keeping']),
         modes=dict(beams='stops', attitude_keeping='drifts', film='unaffected'),
         spread='as fast as the broadcast; every tile that trusts it computes the same wrong target',
         detection='cross-checks against neighbours\' ranging, star trackers and the Sun',
         restart_days=dict(without=[.5, 7.], with_=[0., 0.]),
         likelihood='occasional: GLONASS 2014, GPS 2016 and Galileo 2019 in six years among four systems',
         precedents=['glonass_2014', 'gps_2016', 'galileo_2019'], kind='design guess'),
    dict(id='F8', name='a collision cascade in the rings',
         reach=dict(without=['beams', 'attitude_keeping', 'film'], with_=['beams', 'attitude_keeping', 'film']),
         modes=dict(beams='stops', attitude_keeping='stops', film='stops'),
         spread='fragments return to the common line of nodes every half orbit, 34-39 hours',
         detection='the tiles\' own ranging and the fleet\'s tracking',
         restart_days=dict(without=[180., 3650.], with_=[180., 3650.]),
         likelihood='no precedent at this scale; the concern of crowded low Earth orbit',
         precedents=['kessler_1978'], kind='design guess'),
    dict(id='F9', name='a meteoroid stream',
         reach=dict(without=['film'], with_=['film']),
         modes=dict(beams='unaffected', attitude_keeping='unaffected', film='drifts'),
         spread='hours, over the whole fleet at once',
         detection='forecast from the parent comets',
         restart_days=dict(without=[0., 0.], with_=[0., 0.]),
         likelihood='expected: showers yearly and storms every few decades, adding holes slowly',
         precedents=['grun_1985'], kind='design guess'),
    dict(id='F10', name='an operations error',
         reach=dict(without=['beams', 'attitude_keeping'], with_=['beams']),
         modes=dict(beams='stops', attitude_keeping='holds', film='unaffected'),
         spread='as fast as the command is accepted',
         detection='two-person and two-system authority; refusal of fleet-wide moves',
         restart_days=dict(without=[.01, 2.], with_=[.01, 1.5]),
         likelihood='expected: global routing and configuration errors take networks down for hours',
         precedents=['iberia_2025', 'starlink_2025'], kind='design guess'),
]


def load(path):
    return json.loads(path.read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pair(f):
    return [f(0), f(1)]


def plate_inertia(sigma, side=DESIGN['tile_side_m']):
    """In-plane moment of inertia of a uniform square plate (kg m2)."""
    return sigma*side**4/12.


def turn_check(attitude):
    """The pitch inertia that reproduces the attitude study's edge-on turn times (a rest-to-rest turn of turn_deg at a
    constant torque, half accelerating and half braking), against a uniform plate at each areal mass."""
    theta = math.radians(attitude['design']['turn_deg'])
    rows = {}
    for label, by_torque in attitude['edge_on_pitch']['turn_hours_by_torque'].items():
        sigma = float(label.split('_')[0])/1e3
        inertia = [float(t.split('_')[0])*(h*HOUR/2)**2/theta for t, h in by_torque.items()]
        rows[label] = dict(inertia_kg_m2=inertia, plate_kg_m2=plate_inertia(sigma),
                           over_plate=[i/plate_inertia(sigma) for i in inertia])
    return rows


def shingle_angles(envelope):
    """Pitch error (rad) at which a tile's corner reaches its shingle neighbour's 150 m gate and its surface, the
    conservative bound the tile study uses: the shingle's normal gap, less the envelope on each side, over the half
    diagonal."""
    gap = DESIGN['pitch_m']*math.sin(math.radians(DESIGN['shingle_deg']))
    half_diagonal = DESIGN['tile_side_m']/math.sqrt(2.)
    gate = math.asin(max(gap-envelope-DESIGN['clearance_m'], 0.)/half_diagonal)
    contact = math.asin((gap-envelope)/half_diagonal)
    return dict(shingle_gap_m=gap, gate_rad=gate, contact_rad=contact)


def attitude_clock(attitude, tile):
    """Times for a tile whose attitude control is lost to reach its shingle neighbour's gate and surface: passive drift
    under the light's torque on the balanced filter (with gravity's scaled by inertia to the design radius), the
    momentum store's contents dumped into the body, and an active wrong command at the torque allocation; and the
    rolled law's faster departure."""
    scheme = attitude['schemes']['balanced_400_m']
    run_sigma = attitude['tile_areal_mass_kg_m2']
    radius_ratio = attitude['radius_km']/DESIGN['design_radius_km']
    radiation = [scheme['radiation_torque_N_m']['median'], scheme['radiation_torque_N_m']['p95']]
    gravity_run = [scheme['gravity_torque_N_m']['median'], scheme['gravity_torque_N_m']['p95']]
    stores = dict(chosen=scheme['store_N_m_s']['trim_0.25']['median'],
                  rolled=attitude['schemes']['rolled_80_deg']['store_N_m_s']['trim_0.25']['median'])
    motor = tile['hardware']['initial_torque_allocation_per_axis_N_m']
    reference = tile['hardware']['variants']['annulus_reference']
    sigmas = dict(DESIGN['areal_masses_kg_m2'], kept_run=run_sigma,
                  reference_tile=reference['total_areal_mass_kg_m2'])
    inertias = {k: plate_inertia(v) for k, v in sigmas.items()}
    inertias['reference_tile'] = reference['locked_inertia_kg_m2'][0][0]
    angles = {f'envelope_{int(e)}_m': shingle_angles(e) for e in DESIGN['envelopes_m']}
    run_inertia = plate_inertia(run_sigma)
    rows = {}
    for name, inertia in inertias.items():
        # Gravity-gradient torque goes with the inertia and the cube of the mean motion; the light's does not change.
        torque = [r+g*inertia/run_inertia*radius_ratio**3 for r, g in zip(radiation, gravity_run)]
        row = dict(inertia_kg_m2=inertia, areal_mass_kg_m2=sigmas[name], torque_N_m=torque)
        for env, a in angles.items():
            drift = {k: [math.sqrt(2*a[k]*inertia/t)/HOUR for t in torque] for k in ('gate_rad', 'contact_rad')}
            dump = {s: {k: a[k]/(h/inertia) for k in ('gate_rad', 'contact_rad')} for s, h in stores.items()}
            command = {k: math.sqrt(2*a[k]*inertia/motor)/60. for k in ('gate_rad', 'contact_rad')}
            row[env] = dict(drift_hours_median_p95_torque=drift, store_dump_seconds=dump, wrong_command_minutes=command)
        rows[name] = row
    # The rolled law (the design's night roll, glow.md): to follow its roll the controller applies the attitude torque
    # I dw/dt + w x Iw, which goes with inertia. When control stops, that torque goes missing beside the light's, so a
    # rolled tile departs from its commanded attitude faster; the shingle angles stand in for its clearance, which the
    # attitude study's envelope gives only for attitudes the whole bundle shares.
    rolled = attitude['schemes']['rolled_80_deg']
    rolled_rows = {}
    for name in ('annulus_allocation', 'mean_tile'):
        inertia = inertias[name]
        missing = [rolled['attitude_torque_N_m'][q]*inertia/run_inertia+rolled['radiation_torque_N_m'][q]
                   for q in ('median', 'p95')]
        a = angles['envelope_0_m']
        rolled_rows[name] = dict(missing_torque_N_m=missing,
                                 hours={k: [math.sqrt(2*a[k]*inertia/t)/HOUR for t in missing]
                                        for k in ('gate_rad', 'contact_rad')})
    guarantee = tile['local_failures']['attitude_guarantee_loss']
    return dict(angles=angles, stores_N_m_s=stores, motor_torque_N_m=motor, by_tile=rows, rolled_law=rolled_rows,
                guarantee_loss_check=[dict(rate_deg_s=g['bounded_error_rate_deg_s'], product_s=g['bound_lost_s'],
                                           screen_s=math.degrees(angles['envelope_50_m']['gate_rad'])
                                           / g['bounded_error_rate_deg_s']) for g in guarantee])


def keeping_clock(photon, bundle, tile, kept):
    """Times to the 150 m gate once keeping is lost: free rings (the bundle, 62.7 g/m2 at 15,000 km), scaled to the
    design loading as a steady differential force (time as the inverse square root of the acceleration), and a single
    tile among held neighbours (the tile study, 20,000 km) at the keeping demand scaled to the design loading."""
    sigma_run = photon['tile_areal_mass_kg_m2']
    scale = math.sqrt(DESIGN['areal_masses_kg_m2']['annulus_allocation']/sigma_run)
    free = dict(edge_h=bundle['continuous']['clearance']['first_below_150_m_h'],
                interior_h=bundle['continuous']['interior_clearance']['first_below_150_m_h'],
                minimum_h=bundle['continuous']['clearance']['at_h'])
    free_design = {k: v*scale for k, v in free.items()}
    demand = [photon['demand']['median_m_s2'], photon['demand']['peak_m_s2']]
    demand_design = [d*sigma_run/DESIGN['areal_masses_kg_m2']['annulus_allocation'] for d in demand]
    cases = tile['local_failures']['cases']
    table = {}
    for env in (1, 50, 100):
        acc, gate, contact = [], [], []
        for a in (1e-6, 1e-5, 1e-4):
            c = cases[f'depth_{env}_a_{a:g}'] if f'depth_{env}_a_{a:g}' in cases else cases[f'depth_{env}_a_{a}']
            acc.append(a)
            gate.append(c['first_safety_crossing_s'])
            contact.append(c['first_envelope_contact_s'])
        table[f'envelope_{env}_m'] = dict(residual_m_s2=acc, gate_s=gate, contact_s=contact)

    def at(env, key, a):
        row = table[env]
        return float(np.interp(math.log10(a), np.log10(row['residual_m_s2']), row[key]))
    single = {f'residual_{i}': dict(residual_m_s2=a, gate_min=at('envelope_50_m', 'gate_s', a)/60.,
                                     contact_min=at('envelope_50_m', 'contact_s', a)/60.)
              for i, a in zip(('median', 'peak'), demand_design)}
    edge_light = (1-bundle['starting_shadow_loss']['top_ring'])/(1-bundle['starting_shadow_loss']['interior_rings'])-1
    # Cross-check: the kept run's keeping cost per orbit, spread over the orbit, against the photon study's demand.
    period = kept['period_h']*HOUR
    cost = dict(rings_mean_m_s_per_orbit=kept['keeping']['rings_mean_m_s_per_orbit'],
                top_ring_m_s_per_orbit=kept['keeping']['top_ring_mean_m_s_per_orbit'])
    cost['mean_acceleration_m_s2'] = [cost['rings_mean_m_s_per_orbit']/period, cost['top_ring_m_s_per_orbit']/period]
    cost['over_photon_median_demand'] = [a/demand[0] for a in cost['mean_acceleration_m_s2']]
    return dict(free_rings_product_h=free, free_rings_design_h=free_design, scale_to_design=scale, kept_cost=cost,
                run=dict(areal_mass_kg_m2=sigma_run, radius_step_m=bundle['design']['radius_step_m'],
                         period_h=bundle['period_h']),
                keeping_demand_m_s2=dict(run=demand, design=demand_design), single_tile_table=table,
                single_tile_design=single, edge_ring_light_excess=edge_light)


def shade_threshold(tile, uv, power, transmitters):
    """Tiles lost at once that would carry the UV transmission past the latch, by shield, for losses spread across the
    stack (each lost slot uncovers its singly covered share) and for contiguous losses (everything under them); the
    service fleet's days to replace them; and the transmitting tiles' share."""
    n = tile['fleet_service']['population']
    single = tile['fleet_service']['ideal_service_lattice']['coverage_multiplicity']['1']
    dispatch = tile['fleet_service']['scenarios']['local_spares_1day']['dispatch_capacity_per_day']
    baseline = sorted(uv['levels']['standard']['transmission_by_hole_factor'].values())
    base = [baseline[0], baseline[-1]]
    rows = {}
    for shield, (lo, hi) in DESIGN['latch_transmission'].items():
        margin = [lo-base[1], hi-base[0]]
        spread = [m/single*n for m in margin]
        contiguous = [m*n for m in margin]
        rows[shield] = dict(latch_transmission=[lo, hi], margin_over_standard=margin, tiles_spread=spread,
                            tiles_contiguous=contiguous, share_of_fleet_contiguous=[c/n for c in contiguous],
                            dispatch_days=[c/dispatch for c in contiguous]+[s/dispatch for s in spread])
    lost = {}
    for variant, count in transmitters.items():
        for k in [0]+DESIGN['collateral_tiles']:
            lost[f'{variant}_collateral_{k}'] = [single*c*(1+k)/n for c in count]
    return dict(tiles=n, single_cover_share=single, standard_level=base, dispatch_per_day=dispatch, latch=rows,
                transmitters_all_lost_uncovered=lost, fleet_rebuild_years=n/dispatch/JULIAN_YEAR_DAYS)


def swollen_loss(loss):
    """Cycle times of the losses the latch leads to: just short of its onset behind the titania stack, and the
    unshielded energy-limited bound of the loss product at quiet Sun."""
    unshielded = sorted(loss['unshielded_energy_limited_kg_s']['quiet'].values())
    rates = dict(short_of_onset=DESIGN['loss_short_of_onset_kg_s'], unshielded_quiet=[unshielded[0], unshielded[-1]])
    return {k: dict(kg_s=v, cycle_years=[DESIGN['atmosphere_kg']/r/(JULIAN_YEAR_DAYS*JULIAN_DAY) for r in v])
            for k, v in rates.items()}


def nodes(layout):
    """The stack's rings ordered by the radius at which each turns with the Sun, as the matched layout ranks them; the
    speed at which radius neighbours pass at the common line of nodes, from the difference of their tilts; the strip
    neighbours' speed; and the encounter rate near the nodes for a fragment whose orbit sweeps the stack's radii."""
    stack = layout['stack']
    table = stack['ideal_radius_km']['by_height']
    h_table = np.array([r['height_km'] for r in table])
    r_table = np.array([r['radius_km'] for r in table])
    heights = np.linspace(h_table[0], h_table[-1], stack['rings'])
    radii = np.interp(heights, h_table, r_table)
    tilts = np.arcsin(heights/radii)
    order = np.argsort(radii)
    r_sorted, t_sorted = radii[order], tilts[order]
    speed = np.sqrt(MOON_GM/(1e3*r_sorted))
    v_pair = .5*(speed[1:]+speed[:-1])
    neighbour = 2*v_pair*np.abs(np.sin(.5*np.diff(t_sorted)))
    opposite = np.sign(t_sorted[1:]) != np.sign(t_sorted[:-1])
    v_strip = np.sqrt(MOON_GM/(1e3*.5*(radii[1:]+radii[:-1])))
    strip = 2*v_strip*np.abs(np.sin(.5*np.diff(tilts)))
    n_mean = math.sqrt(MOON_GM/(1e3*DESIGN['design_radius_km'])**3)
    slide = 1.5*n_mean*1e3                      # along-track slide per km of radius step, m/s
    widest = 2*math.sqrt(MOON_GM/(1e3*r_sorted.mean()))*math.sin(.5*(tilts.max()-tilts.min()))
    half_width = DESIGN['tile_side_m']/2
    encounters = {}
    for tilt_deg in DESIGN['fragment_tilts_deg']:
        dtheta = np.abs(math.radians(tilt_deg)-tilts)
        x = np.minimum(1., half_width/np.maximum(1e3*radii*dtheta, 1e-9))
        p = 2/math.pi*np.arcsin(x)
        fast = 2*np.sqrt(MOON_GM/(1e3*radii))*np.sin(.5*dtheta) > DESIGN['fast_m_s']
        encounters[f'tilt_{tilt_deg:g}_deg'] = dict(per_orbit=float(2*p.sum()), fast_per_orbit=float(2*p[fast].sum()),
                                                   fast_share=float(p[fast].sum()/p.sum()))
    return dict(rings=int(stack['rings']),
                tilt_range_deg=[float(np.degrees(tilts.min())), float(np.degrees(tilts.max()))],
                radius_range_km=[float(radii.min()), float(radii.max())],
                radius_neighbours=dict(opposite_side_share=float(opposite.mean()),
                                       speed_m_s=dict(median=float(np.median(neighbour)),
                                                      p90=float(np.percentile(neighbour, 90)),
                                                      max=float(neighbour.max())),
                                       share_over_100_m_s=float((neighbour > 100.).mean())),
                strip_neighbours_speed_m_s=dict(median=float(np.median(strip)), max=float(strip.max())),
                slide_m_s_per_km_step=slide, widest_pair_m_s=widest,
                orbit_hours=float(2*math.pi/n_mean/HOUR), fragment_encounters=encounters)


def upper_air(loss, shade):
    """Order-of-magnitude bracket on the time the air above the base needs to warm toward the latch: the warmed share
    of the column's heat capacity times an assumed warming, over the threshold heat. The products give no response
    time; this bounds its order only."""
    base = loss['assumptions']['base_pressure_pa']
    column = base/MOON_SURFACE_GRAVITY
    shares = [math.exp(-e) for e in DESIGN['warmed_efolds']]
    capacity = [s*column*DESIGN['cp_J_kg_K'] for s in sorted(shares)]
    heat = [min(min(v) for v in DESIGN['latch_heat_W_m2'].values()),
            max(max(v) for v in DESIGN['latch_heat_W_m2'].values())]
    fastest = capacity[0]*DESIGN['warming_K'][0]/heat[1]
    slowest = capacity[1]*DESIGN['warming_K'][1]/heat[0]
    replace = [min(min(r['dispatch_days']) for r in shade['latch'].values()),
               max(max(r['dispatch_days']) for r in shade['latch'].values())]
    return dict(base_pressure_Pa=base, column_kg_m2=column, warmed_capacity_J_m2_K=capacity, threshold_heat_W_m2=heat,
                days=[fastest/JULIAN_DAY, slowest/JULIAN_DAY], climate_efolding_years=DESIGN['climate_efolding_years'],
                service_days_to_replace_threshold=replace)


def recoil(power, attitude, bundle):
    """The radio recoil on a transmitting tile, in use and at its rating, against a film tile's sail force and against
    the push of sunlight a whole-tile collector absorbs; the step a ring-wide beam stop puts on each transmitting tile,
    against the edge ring's extra light that brings free edge rings within 150 m; and the recoil's acceleration by the
    transmitter's areal mass."""
    ring = power['ring']
    tiles = ring['transmitting_tiles']
    in_use = pair(lambda e: ring['rf_in_use_tw'][e]*1e12/tiles[e]/SPEED_OF_LIGHT)
    rated = pair(lambda e: ring['tile']['whole_tile_rf_gw'][e]*1e9/SPEED_OF_LIGHT)
    small = ring['small_collector_variant']
    small_in_use = pair(lambda e: ring['rf_in_use_tw'][e]*1e12/small['transmitting_tiles'][e]/SPEED_OF_LIGHT)
    sail = attitude['schemes']['balanced_400_m']['sail_force_N']
    absorbed = ring['tile']['sunlight_gw']*1e9/SPEED_OF_LIGHT
    area = DESIGN['tile_side_m']**2
    electricity_rated = pair(lambda e: ring['tile']['whole_tile_rf_gw'][e]*1e9/DESIGN['transmitter_efficiency'])
    reference_sigma = pair(lambda e: electricity_rated[e]/1e3*DESIGN['transmitter_kg_per_kw']/area)
    sigmas = dict(film_allocation=DESIGN['areal_masses_kg_m2']['annulus_allocation'],
                  thin_film_sandwich=DESIGN['thin_film_transmitter_kg_m2'], reference_1978=reference_sigma)
    accel = {}
    for name, s in sigmas.items():
        s_pair = s if isinstance(s, list) else [s, s]
        accel[name] = dict(areal_mass_kg_m2=s_pair, in_use_m_s2=[in_use[e]/(s_pair[e]*area) for e in (0, 1)])
    edge = (1-bundle['starting_shadow_loss']['top_ring'])/(1-bundle['starting_shadow_loss']['interior_rings'])-1
    # A whole-tile collector that absorbs all its light shares the rings' sail loading only at the areal mass that
    # gives it a film tile's acceleration (median sail force on the annulus allocation); a crude match of medians.
    film_acceleration = sail['median']/(DESIGN['areal_masses_kg_m2']['annulus_allocation']*area)
    matched = absorbed/film_acceleration/area
    return dict(transmitting_tiles=tiles, in_use_N=in_use, rated_N=rated, tenth_tile_in_use_N=small_in_use,
                film_tile_sail_N=dict(median=sail['median'], p95=sail['p95']),
                in_use_over_film_sail=[f/sail['median'] for f in in_use],
                whole_tile_absorbed_N=absorbed,
                step_share_of_absorbed_push=dict(in_use=[f/absorbed for f in in_use],
                                                 rated=[f/absorbed for f in rated]),
                edge_ring_light_excess=edge, acceleration_by_areal_mass=accel,
                film_tile_acceleration_m_s2=film_acceleration, matched_sandwich_kg_m2=matched,
                reference_over_matched=[s/matched for s in reference_sigma])


def beams(power, faults):
    """The critical loads with how fast each needs power back; for each fault class the days of storage its restart
    asks (the upper end of its restart time), the energy, and what the ground's store must still hold once the local
    flows carry their share (the low end pairs the least load with the most local supply, the high end the most load
    with the least); the share of the ground's critical loads the local supplies carry; the towns' weeks on their
    reserve lift."""
    rf = power['reserves']['ring_fault']
    crit = rf['critical_tw']
    total, ground, towns = rf['critical_total_tw'], rf['ground']['critical_tw'], rf['towns']['critical_tw']
    local = DESIGN['local_tw']
    flows = [sum(v[0] for k, v in local.items() if k != 'roofs'), sum(v[1] for k, v in local.items() if k != 'roofs')]
    with_roofs = [flows[0]+local['roofs'][0], flows[1]+local['roofs'][1]]
    # Conservative pairing: the least local supply against the larger load, the most against the smaller.
    share = dict(flows=[flows[0]/ground[1], flows[1]/ground[0]], flows_and_roofs=[with_roofs[0]/ground[1],
                                                                                 min(1., with_roofs[1]/ground[0])])
    tolerances = {k: dict(v, range_h=[x/HOUR for x in v['range']]) for k, v in DESIGN['tolerance_s'].items()}
    stores = {}
    grid = DESIGN['grid_restoration_days'][1]
    for f in faults:
        for side in ('without', 'with_'):
            restart = f['restart_days'][side][1]
            stopped = 'beams' in f['reach'][side] and restart > 0
            days = restart+grid if stopped else 0.
            stores.setdefault(f['id'], {})[side.rstrip('_')] = dict(
                restart_days=restart, days=days, energy_twh=[t*24*days for t in total],
                ground_store_twh_after_local_flows=[max(ground[e]-flows[0 if e else 1], 0.)*24*days for e in (0, 1)],
                within_a_week=days <= rf['days'])
    weeks_on_reserve_lift = [1/s for s in rf['towns']['share_of_reserve_lift']]
    return dict(critical_tw=crit, critical_total_tw=total, ground_tw=ground, towns_tw=towns, week_twh=rf['energy_twh'],
                local_tw=dict(flows=flows, with_roofs=with_roofs), local_share_of_ground=share,
                towns_weeks_on_reserve_lift=weeks_on_reserve_lift, tolerance=tolerances, stores=stores)


def separation(shade, tile):
    """The separation's added mass and power across the fleet, against the mean tile and the inventory; the largest
    update wave one cause may touch; and how long a staged campaign over the fleet takes."""
    n = tile['fleet_service']['population']
    mass_t = [m*n/1e9 for m in DESIGN['separate_controller_kg']]
    power_gw = [w*n/1e9 for w in DESIGN['separate_controller_W']]
    budget = min(min(r['tiles_contiguous']) for r in shade['latch'].values())
    waves = [budget/(1+k)/DESIGN['margin'] for k in DESIGN['collateral_tiles']]
    campaigns = {f'wave_{int(w)}': [math.ceil(n/w)*d/JULIAN_YEAR_DAYS for d in DESIGN['update_soak_days']]
                 for w in waves}
    # With a second, dissimilar stack that takes control when the updated one fails, a wave leaves a tile
    # uncontrolled only when the fallback fails too, with probability p per tile: the wave may be 1/p times larger.
    fallback = {}
    for p in DESIGN['fallback_failure']:
        wave = [min(w/p, n) for w in waves]
        fallback[f'p_{p:g}'] = dict(wave_tiles=wave, waves=[math.ceil(n/w) for w in wave],
                                    campaign_days=[[math.ceil(n/w)*d for d in DESIGN['update_soak_days']]
                                                   for w in wave])
    return dict(fleet_mass_Mt=mass_t,
                share_of_mean_tile=[m/DESIGN['mean_tile_kg'] for m in DESIGN['separate_controller_kg']],
                share_of_inventory=[m/(DESIGN['fleet_inventory_Gt']*1e3) for m in mass_t], fleet_power_GW=power_gw,
                smallest_latch_budget_tiles=budget, largest_wave_tiles=waves, campaign_years=campaigns,
                with_fallback=fallback)


def meteoroid_storm(uv, shade):
    """The share of the aperture a meteoroid storm holes in a day: the swarm product's yearly holing at the particle's
    own size, times the hole's area at the largest hole factor, the storm's enhancement and its duration; against the
    smallest margin to the latch."""
    rate = uv['holing_rate_per_year_hole_factor_1']*DESIGN['storm_hole_factor']**2
    holed = [rate*e*DESIGN['storm_days']/JULIAN_YEAR_DAYS for e in DESIGN['storm_enhancement']]
    lowest = min(min(r['margin_over_standard']) for r in shade['latch'].values())
    return dict(holing_per_year_at_largest_holes=rate, holed_share=holed,
                over_smallest_margin=[h/lowest for h in holed])


def trim_at_stake():
    """How much of the shade rides on the reflectivity trim: the share of the sunlight the trim could pass to the Moon
    if its devices traded reflection for transmission, and the glow it would add if they traded it for absorption.
    The UV filter is the titania, which absorbs whatever the trim does."""
    t = DESIGN['trim']
    share = [t['range']*t['reflected'][0]*t['cover'][0], t['range']*t['reflected'][1]*t['cover'][1]]
    return dict(sunlight_share_if_transmitted=share,
                glow_W_m2_if_absorbed=[share[0]*100*t['glow_W_m2_per_percent'][0],
                                       share[1]*100*t['glow_W_m2_per_percent'][1]],
                uv_transmission_at_stake=0., basis=t['basis'])


def targets(shade):
    """The annual frequency a shade-opening fault may have for the accepted probability over each horizon, beside the
    aviation and nuclear standards; and the largest number of tiles one cause may leave uncontrolled at once."""
    hours = JULIAN_YEAR_DAYS*24
    n = shade['tiles']
    budget = {shield: [c/(1+k) for k in DESIGN['collateral_tiles'] for c in r['tiles_contiguous'][:1]]
              for shield, r in shade['latch'].items()}
    # A common stress such as a solar proton event reaches every tile: the chance each tile loses control in it may be
    # at most the budget over the fleet, and a tenth of that with the margin, if the losses are independent.
    per_tile = {shield: dict(no_margin=[b/n for b in v], with_margin=[b/n/DESIGN['margin'] for b in v])
                for shield, v in budget.items()}
    record = DESIGN['halloween_2003']
    storm = dict(anomaly_share=[record['anomalies']/record['operating'][1], record['anomalies']/record['operating'][0]],
                 loss_share=[record['out_of_service']/record['operating'][1],
                             record['out_of_service']/record['operating'][0]])
    loosest = max(max(v['no_margin']) for v in per_tile.values())
    storm['loss_share_over_loosest_budget'] = [s/loosest for s in storm['loss_share']]
    return dict(per_year={f'{h:g}_years': DESIGN['accepted_probability']/h for h in DESIGN['horizons_years']},
                aviation_per_year_of_continuous_operation=DESIGN['aviation_catastrophic_per_flight_hour']*hours,
                nuclear_large_release_per_reactor_year=DESIGN['nuclear_large_release_per_reactor_year'],
                uncontrolled_at_once_budget_tiles=budget, per_tile_loss_of_control_per_event=per_tile,
                halloween_2003=storm)


def fault_table():
    rows = []
    for f in FAULTS:
        rows.append(dict(f, reach=dict(without=f['reach']['without'], with_separation=f['reach']['with_']),
                         restart_days=dict(without=f['restart_days']['without'],
                                           with_separation=f['restart_days']['with_']),
                         opens_shade_without=bool({'film', 'attitude_keeping'} & set(f['reach']['without'])),
                         opens_shade_with=bool({'film', 'attitude_keeping'} & set(f['reach']['with_']))))
    return rows


def main():
    power, attitude, photon, bundle, kept, tile, layout, uv, loss = (load(p) for p in INPUTS[:-1])
    transmitters = dict(whole_tile=power['ring']['transmitting_tiles'],
                        tenth_tile=power['ring']['small_collector_variant']['transmitting_tiles'])
    shade = shade_threshold(tile, uv, power, transmitters)
    out = dict(
        schema=SCHEMA,
        producer=dict(files={f: digest(ROOT/f) for f in FILES}, constants=constants_used(FILES)),
        inputs={str(p.relative_to(ROOT)): digest(p) for p in INPUTS},
        evidence=('Closed-form screen on the named products; no dynamics are run. Clocks come from the products\' '
                  'torques, inertias, keeping demand, free-ring and failed-tile runs, scaled to the design loading '
                  'where stated; node speeds from the rolled layout\'s radius at which each height\'s ring turns with '
                  'the Sun; the upper air\'s time is an order-of-magnitude bracket on assumed warming; fault classes, '
                  'restart times, tolerances and the separation\'s hardware are design guesses anchored on the '
                  'precedents in ring_fault_sources.json. Each value carries its kind in ring_fault.md.'),
        reading_rule=('Times in seconds, minutes, hours or days as named; forces N; accelerations m/s2; powers TW; '
                      'energies TWh. Pairs are the low-cost and high-cost ends of the power screen, or the low and '
                      'high ends of a range as named. "with_separation" applies the design rules of ring_fault.md. '
                      'Tile counts are slots of the 28.2 million; "contiguous" losses uncover everything under them, '
                      '"spread" losses only the singly covered share.'),
        design=DESIGN,
        inertia_check=turn_check(attitude),
        attitude=attitude_clock(attitude, tile),
        keeping=keeping_clock(photon, bundle, tile, kept),
        shade=shade,
        nodes=nodes(layout),
        recoil=recoil(power, attitude, bundle),
        trim=trim_at_stake(),
        faults=fault_table(),
    )
    out['meteoroid_storm'] = meteoroid_storm(uv, shade)
    out['swollen_loss'] = swollen_loss(loss)
    out['upper_air'] = upper_air(loss, shade)
    out['beams'] = beams(power, FAULTS)
    out['separation'] = separation(shade, tile)
    out['targets'] = targets(shade)
    ok_turn = all(abs(r-1) < .01 for v in out['inertia_check'].values() for r in v['over_plate'])
    ok_guarantee = all(abs(g['screen_s']/g['product_s']-1) < .01 for g in out['attitude']['guarantee_loss_check'])
    week = power['reserves']['ring_fault']
    ok_week = all(math.isclose(t*24*week['days'], e, rel_tol=2e-3)
                  for t, e in zip(week['critical_total_tw'], week['energy_twh']))
    ok_labels = all(f['kind'] in LABELS and set(f['modes']) == set(FUNCTIONS) for f in FAULTS)
    out['checks'] = dict(turn_times_match_plate_inertia=ok_turn, guarantee_loss_matches_tile_study=ok_guarantee,
                         week_matches_power_screen=ok_week, fault_classes_complete=ok_labels)
    OUT.write_text(json.dumps(out, indent=1)+'\n')
    print(json.dumps(dict(checks=out['checks'], nodes=out['nodes']['radius_neighbours'],
                          upper_air_days=out['upper_air']['days'], recoil_in_use_N=out['recoil']['in_use_N']),
                     indent=1))


if __name__ == '__main__':
    main()
