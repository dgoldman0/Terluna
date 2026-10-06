"""Exhaust isolation of the held screen (step 3 of the integrated comparison).

A screen held upstream of the Moon needs a sunward push that sunlight cannot give, so its reaction mass
leaves into the Moon's hemisphere. This runner takes the zoned design point (26 g/m2 window, 5 g/m2
annulus) on the selected moving trajectory of holding.json and measures three things against the proposed
requirement S7 of research/studies/protection_architecture/requirements.md.

1. The allowance. The September report's energy diagnostic (protection/report.md, section 8) lets a tenth
   of the energy deposited in the escape region do escape work against the Moon's binding energy at three
   lunar radii, so each kg/s of loss budget allows a deposited power, stated here as a share of the jet
   power.
2. The direct path. Each node's propulsive need is split between two thrusters canted from the resultant
   axis, in the plane that keeps both plumes furthest from the Moon and, for comparison, in the plane
   through the Moon. The share of each plume that leaves straight into the protected sphere (four lunar
   radii) is integrated over that sphere's cone from two measured plume profiles: a gridded ion thruster
   (NEXT) and a Hall thruster (BPT-4000), at cants of 45, 60 and 75 degrees.
3. The slow gas. Unionized propellant leaves the thrusters at thermal speed. Test particles leave the
   screen with their node's own velocity, a cosine law about each thruster axis and the effusive speed
   distribution of oxygen molecules or xenon atoms at 500 K, and fly under DE440s point-mass gravity for
   60 days with a fixed-step fourth-order Runge-Kutta integrator. The share that enters the protected
   sphere, by day, gives the share delivered for any neutral lifetime.
4. The magnetosphere. The September magnets hold off charged particles whose gyroradius at the magnetopause is
   small against their stand-off; they cannot hold the fast atoms that charge exchange makes in each plume, about
   a tenth of the beam on the same paths, or the slow gas. With the magnets, the direct path therefore leaves its
   fast-neutral share, and the runner gives the share of the slow gas that must be captured at the thrusters for
   each loss budget.

The solar wind's pickup of the charged exhaust, and the magnetosphere's leak, lie outside every model here and stay
open.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.exhaust_isolation
"""
from __future__ import annotations

import hashlib
import inspect
import json
import resource
import time
from pathlib import Path

import numpy as np

from shared import constants as K
from shared.provenance import constants_used
from protection.dynamics.cycling import CyclingEnvironment
from protection.dynamics.ephemeris import DEFAULT_KERNEL, Ephemeris
from protection.dynamics.model import geometry
from protection.dynamics.optical import length
from .run import CORE_DIVERTED, SCENARIO
from .zoned_aperture import HOLDING, OUT as ZONED, zoned_aperture

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RUN = ROOT/'research/runs/solar_shield_array/exhaust_isolation'
OUT = HERE/'results/exhaust_isolation.json'
FILES = ['research/studies/solar_shield_array/exhaust_isolation.py', 'research/studies/solar_shield_array/zoned_aperture.py',
         'research/studies/solar_shield_array/aperture.py', 'research/studies/solar_shield_array/run.py',
         'protection/dynamics/model.py', 'protection/dynamics/optical.py', 'protection/dynamics/ephemeris.py',
         'protection/dynamics/cycling.py']
DESIGN = dict(case='core_26g_annulus_5g', core_sigma_kg_m2=.026, annulus_sigma_kg_m2=.005, protected_radii=4,
              cants_deg=[45., 60., 75.], budgets_kg_s=[1., 10., 100.], escape_work_fraction=.1, binding_radii=3.,
              unionized_share=[.05, .15], gas_temperature_K=500., molar_mass_kg_mol=dict(O2=.031998, Xe=.131293),
              launch_days=30., launch_every_days=2., flight_days=60., particles_per_thruster=6, step_s=300.,
              check_step_s=150., lifetimes_days=[10., 30., 100.], near_earth_radii=10., escape_m=3e9,
              direct_every_h=24., seed=20261006,
              # Charge exchange in the plume turns about a tenth of the beam into fast atoms that keep the beam's
              # direction (Goebel and Katz, figure 8-11: 0.74 A of charge-exchange ions in space for a 3 kW, 300 V
              # Hall thruster, whose discharge carries 10 A).
              fast_neutral_share_of_beam=.1,
              # The September magnets: stand-off 10 lunar radii (atmosphere/loss_response), and the field that
              # balances 2 nPa of solar-wind pressure at the magnetopause (protection/report.md, section 10).
              standoff_radii=10., magnetopause_field_T=71e-9,
              # The loss response's typical solar wind at 1 au (atmosphere/loss_response/model.py, SOLAR_WIND).
              solar_wind_speed_m_s=4.0e5, budgets_with_magnets_kg_s=[1., 10.])
# Relative energetic-ion intensity per steradian against angle from the thruster axis, read from the published
# figures to about +-30% in each value.
PLUMES = {
    'gridded_ion_NEXT': dict(
        source=('Young, Matlock, Nakles and Crofton, Far Field Plume Distribution and Divergence for NEXT: DART '
                'Mission, AIAA SciTech Forum 2019, figure 5 (retarding potential analyzer above 35 V, grid-edge '
                'angle). Beyond 50 degrees the last decade is extended at its own slope; the paper reports a fall '
                'of over five decades at high angles.'),
        angle_deg=[0., 10., 15., 20., 25., 30., 35., 38., 40., 42., 45., 50.],
        relative=[1., .69, .25, .0625, .031, .0125, .0052, .0021, 6.3e-4, 4.2e-4, 3.1e-4, 1.7e-4]),
    'hall_BPT4000': dict(
        source=('Goebel and Katz, Fundamentals of Electric Propulsion: Ion and Hall Thrusters (JPL, 2008), '
                'figure 8-12 (BPT-4000 ion flux above 50 V at 1 m, laboratory); the flight data of figure 8-6 '
                '(SPT-100 on Express-A) fall alike to about 5% of the axial value at 40 degrees.'),
        angle_deg=[0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90.],
        relative=[.8, 1., .6, .28, .14, .08, .044, .024, .014, .008, .0048, .0028, .0018, .0014, .0012, 8e-4,
                  3.2e-4, 1.6e-4, 8e-5]),
}
SOURCES = [
    'protection/report.md, section 8: exhaust isolation and the energy diagnostic.',
    'atmosphere/loss_response/README.md: returning pickup ions sputter 1-10 molecules each.',
    PLUMES['gridded_ion_NEXT']['source'],
    PLUMES['hall_BPT4000']['source'],
    'Goebel and Katz (2008), chapter 7: a mass utilization efficiency of 95% reported for SPT thrusters, so at '
    'least about 5% of the propellant leaves unionized; 15% is the range assumed for less efficient thrusters.',
    'Heays, Bosman and van Dishoeck, A&A 602, A105 (2017), tables 18-19 scaled to 1 au: O2 photodissociation '
    'about 2.4e-6 per s; photoionisation of O, N2 and O2 about 0.5-1.3e-7 per s. Lifetimes of 10-100 days span '
    'photoionisation with charge exchange and electron impact in the solar wind.',
    'Goebel and Katz (2008), figure 8-11: a space-condition plume model of the BPT-4000 at 3 kW makes 0.74 A of '
    'charge-exchange ions, about a tenth of its beam; each exchange leaves a fast atom on the beam ion\'s path.',
    'protection/report.md, section 10: 71 nT balances 2 nPa of solar-wind pressure; the September magnets '
    '(1.5e21 A m2) stand off at 10 lunar radii (atmosphere/loss_response/README.md; requirements C1). Magnetic '
    'protection may divert charged exhaust and cannot be credited with removing fast neutrals (section 8).',
]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def unit(v):
    """Unit vectors along the last axis; a zero vector stays zero."""
    n = length(v)
    return np.divide(v, n[..., None], out=np.zeros_like(v), where=n[..., None] > 0)


def perpendicular(v):
    """A unit vector perpendicular to each unit vector v."""
    seed = np.where(np.abs(v[..., :1]) < .9, np.array([1., 0., 0.]), np.array([0., 1., 0.]))
    return unit(np.cross(v, seed))


def intensity(profile, theta):
    """Relative energetic-ion intensity per steradian at theta (rad) from the axis; nothing leaves backward."""
    a, r = np.radians(profile['angle_deg']), np.log(profile['relative'])
    slope = (r[-1]-r[-2])/(a[-1]-a[-2])
    log = np.where(theta <= a[-1], np.interp(theta, a, r), r[-1]+slope*(theta-a[-1]))
    return np.where(theta <= np.pi/2, np.exp(log), 0.)


def total_intensity(profile):
    theta = np.linspace(0., np.pi/2, 20001)
    return float(np.trapezoid(intensity(profile, theta)*2*np.pi*np.sin(theta), theta))


def cap_share(profile, axis, centre, half, n_r=12, n_phi=24):
    """Share of a plume along axis (..., 3) leaving into the cone of half-angle half (...) about centre (..., 3)."""
    e1 = perpendicular(centre)
    e2 = np.cross(centre, e1)
    s = half[..., None]*np.sqrt((np.arange(n_r)+.5)/n_r)
    phi = 2*np.pi*(np.arange(n_phi)+.5)/n_phi
    around = (np.cos(phi)[:, None]*e1[..., None, None, :]+np.sin(phi)[:, None]*e2[..., None, None, :])
    direction = np.cos(s)[..., None, None]*centre[..., None, None, :]+np.sin(s)[..., None, None]*around
    theta = np.arccos(np.clip(np.sum(direction*axis[..., None, None, :], axis=-1), -1., 1.))
    solid = 2*np.pi*(1-np.cos(half))
    return intensity(profile, theta).mean(axis=(-1, -2))*solid/total_intensity(profile)


def chunked_share(profile, axis, centre, half, size=20):
    """cap_share over the first (date) axis in chunks, to bound memory."""
    return np.concatenate([cap_share(profile, axis[i:i+size], centre[i:i+size], half[i:i+size])
                           for i in range(0, len(axis), size)])


def plume_axes(force, moon, cant, plane):
    """The two thruster axes of a pair canted from the resultant; plane 'clear' keeps both off the Moon's side."""
    f, m = unit(force), unit(moon)
    u = np.cross(f, m) if plane == 'clear' else m-np.sum(m*f, axis=-1)[..., None]*f
    # Along the Moon's direction either plane is undefined and any perpendicular serves.
    u = np.where(length(u)[..., None] > 1e-9, unit(u), perpendicular(f))
    back = -f*np.cos(cant)
    return back+u*np.sin(cant), back-u*np.sin(cant)


def allowance(propellant_kg_s, jet_power_W):
    binding = K.MOON_GM/(DESIGN['binding_radii']*K.MOON_RADIUS)
    rows = []
    for budget in DESIGN['budgets_kg_s']:
        power = budget*binding/DESIGN['escape_work_fraction']
        rows.append(dict(budget_kg_s=budget, deposited_power_MW=power/1e6, share_of_jet_power=power/jet_power_W,
                         exhaust_kg_s_at_jet_speed=power/(.5*SCENARIO['propulsion']['exhaust_velocity_m_s']**2)))
    return dict(binding_energy_MJ_kg=binding/1e6, jet_power_TW=jet_power_W/1e12, propellant_kg_s=propellant_kg_s,
                per_budget=rows)


def direct_path(fields):
    """Clearance of the plume axes from the protected sphere and the share of measured plumes sent into it."""
    t = fields['t']
    every = int(round(DESIGN['direct_every_h']*3600./(t[1]-t[0])))
    q, res = fields['q'][::every], fields['residual'][::every]
    moon = -q
    half = np.arcsin(np.minimum(DESIGN['protected_radii']*K.MOON_RADIUS/length(q), 1.))
    flow = fields['sigma'][None, :]*length(res)*fields['weights'][None, :]
    out = {}
    for plane in ('clear', 'through_moon'):
        rows = []
        for cant in DESIGN['cants_deg']:
            a, b = plume_axes(res, moon, np.radians(cant), plane)
            gap = [np.degrees(np.arccos(np.clip(np.sum(x*unit(moon), axis=-1), -1., 1.))-half) for x in (a, b)]
            row = dict(cant_deg=cant, propellant_and_power_factor=float(np.cos(np.radians(45.))/np.cos(np.radians(cant))),
                       clearance_min_deg=float(np.minimum(*gap).min()),
                       clearance_flow_median_deg=float(np.median(np.minimum(*gap))))
            for name, profile in PLUMES.items():
                share = .5*(chunked_share(profile, a, unit(moon), half)+chunked_share(profile, b, unit(moon), half))
                row[f'{name}_share'] = float(np.sum(flow*share)/np.sum(flow))
                row[f'{name}_share_max'] = float(share.max())
                # Fast atoms from charge exchange keep the beam's directions and pass any magnetosphere.
                row[f'{name}_fast_neutral_share'] = DESIGN['fast_neutral_share_of_beam']*row[f'{name}_share']
            rows.append(row)
        out[plane] = rows
    out['protected_half_angle_deg'] = [float(np.degrees(half.min())), float(np.degrees(half.max()))]
    out['samples'] = dict(dates=int(len(q)), nodes=int(q.shape[1]))
    out['profile_totals'] = {name: total_intensity(p) for name, p in PLUMES.items()}
    return out


def magnetosphere(direct, implied, allowance_rows, jet_power_W):
    """What the September magnets hold back, and the unionized-gas capture each loss budget then needs.

    Exhaust ions and picked-up ions are held off where their gyroradius at the magnetopause is small against the
    stand-off. Fast neutrals (the charge-exchanged share of each direct-path plume) and the slow gas pass. The
    capture is the share of the slow gas that must be stopped at the thrusters so that fast neutrals and slow gas
    together stay within the allowance, with a 30-day neutral lifetime.
    """
    standoff = DESIGN['standoff_radii']*K.MOON_RADIUS
    gyro = {}
    for name, molar in DESIGN['molar_mass_kg_mol'].items():
        def radius(speed):
            return molar/K.AVOGADRO*speed/(K.ELEMENTARY_CHARGE*DESIGN['magnetopause_field_T'])
        gyro[name] = dict(exhaust_ion_km=radius(SCENARIO['propulsion']['exhaust_velocity_m_s'])/1e3,
                          pickup_ion_km=radius(DESIGN['solar_wind_speed_m_s'])/1e3,
                          pickup_diameter_over_standoff=2*radius(DESIGN['solar_wind_speed_m_s'])/standoff)
    allowed = {row['budget_kg_s']: row['deposited_power_MW'] for row in allowance_rows}
    slow = [mw for species in implied.values() for mw in species['30_days']['deposited_MW']]
    rows = []
    for budget in DESIGN['budgets_with_magnets_kg_s']:
        for row in direct['clear']:
            for name in PLUMES:
                fast = row[f'{name}_fast_neutral_share']*jet_power_W/1e6
                left = allowed[budget]-fast
                capture = [max(0., 1-left/mw) for mw in slow] if left > 0 else None
                rows.append(dict(budget_kg_s=budget, thruster=name, cant_deg=row['cant_deg'], fast_neutral_MW=fast,
                                 allowance_MW=allowed[budget],
                                 slow_gas_capture_needed=[min(capture), max(capture)] if capture else None))
    return dict(standoff_km=standoff/1e3, magnetopause_field_nT=DESIGN['magnetopause_field_T']*1e9,
                gyroradius=gyro, slow_gas_30_day_MW=[min(slow), max(slow)], with_magnets=rows)


def effusive(rng, n, axis, molar_mass):
    """Velocities of n molecules leaving along a unit axis with a cosine law and the effusive speed distribution."""
    a = np.sqrt(2*K.BOLTZMANN*DESIGN['gas_temperature_K']*K.AVOGADRO/molar_mass)
    speed = a*np.sqrt(-np.log(rng.random(n)*rng.random(n)))
    cos_t = np.sqrt(rng.random(n))
    sin_t = np.sqrt(1-cos_t**2)
    phi = 2*np.pi*rng.random(n)
    e1 = perpendicular(axis[None])[0]
    e2 = np.cross(axis, e1)
    d = cos_t[:, None]*axis+sin_t[:, None]*(np.cos(phi)[:, None]*e1+np.sin(phi)[:, None]*e2)
    return speed[:, None]*d


def fly(env, t0, states, days, step):
    """Fixed-step RK4 under point-mass gravity; returns fate (0 flying, 1 protected sphere, 2 near Earth,
    3 escaped), entry time after launch (s) and entry speed (m/s)."""
    n = len(states)
    y, t = states.copy(), t0
    fate = np.zeros(n, int)
    t_in, v_in = np.full(n, np.nan), np.full(n, np.nan)
    protected = DESIGN['protected_radii']*K.MOON_RADIUS

    def f(tt, s):
        return np.hstack([s[:, 3:], env.gravity(tt, s[:, :3])])
    for _ in range(int(round(days*K.JULIAN_DAY/step))):
        idx = np.flatnonzero(fate == 0)
        if not len(idx):
            break
        s = y[idx]
        k1 = f(t, s)
        k2 = f(t+step/2, s+step/2*k1)
        k3 = f(t+step/2, s+step/2*k2)
        k4 = f(t+step, s+step*k3)
        s = s+step/6*(k1+2*k2+2*k3+k4)
        t += step
        y[idx] = s
        r = length(s[:, :3])
        earth = length(s[:, :3]-env.at(t)['positions']['earth'])
        new = np.select([r < protected, earth < DESIGN['near_earth_radii']*K.EARTH_RADIUS, r > DESIGN['escape_m']],
                        [1, 2, 3], 0)
        fate[idx] = new
        hit = idx[new == 1]
        t_in[hit], v_in[hit] = t-t0, length(s[new == 1, 3:])
    return fate, t_in, v_in


def slow_gas(fields, env, run_key):
    """Test particles of unionized propellant from every node, every launch date, both thrusters of each pair."""
    RUN.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(DESIGN['seed'])
    t = fields['t']
    cant = np.radians(DESIGN['cants_deg'][0])
    launches = np.arange(0., DESIGN['launch_days'], DESIGN['launch_every_days'])*K.JULIAN_DAY
    k = DESIGN['particles_per_thruster']
    records = {name: [] for name in DESIGN['molar_mass_kg_mol']}
    cpu = 0.
    for i, tl in enumerate(launches):
        j = int(np.argmin(np.abs(t-tl)))
        q, v, res = fields['q'][j], fields['v'][j], fields['residual'][j]
        a, b = plume_axes(res, -q, cant, 'clear')
        flow = fields['sigma']*length(res)*fields['weights']
        states, weights, species = [], [], []
        for name, molar in DESIGN['molar_mass_kg_mol'].items():
            for node in range(len(q)):
                for axis in (a[node], b[node]):
                    kick = effusive(rng, k, axis, molar)
                    states.append(np.hstack([np.repeat(q[node][None], k, axis=0), v[node]+kick]))
                    weights.append(np.full(k, flow[node]/(2*k)))
                    species.extend([name]*k)
        states, weights, species = np.vstack(states), np.concatenate(weights), np.array(species)
        path = RUN/f'launch_{i}_{run_key}.npz'
        if not path.exists():
            start = time.process_time()
            fate, t_in, v_in = fly(env, t[j], states, DESIGN['flight_days'], DESIGN['step_s'])
            tmp = path.with_name(path.stem+'.tmp.npz')
            np.savez(tmp, fate=fate, t_in=t_in, v_in=v_in, cpu=time.process_time()-start)
            tmp.rename(path)
        z = np.load(path)
        cpu += float(z['cpu'])
        for name in records:
            m = species == name
            records[name].append(dict(weights=weights[m], fate=z['fate'][m], t_in=z['t_in'][m], v_in=z['v_in'][m]))
        if i == 0:
            first = dict(states=states, t0=t[j], fate=z['fate'])
    # The same first launch at half the step checks the integrator.
    check = RUN/f'check_{run_key}.npz'
    if not check.exists():
        fate, _, _ = fly(env, first['t0'], first['states'], DESIGN['flight_days'], DESIGN['check_step_s'])
        np.savez(check, fate=fate)
    fine = np.load(check)['fate']
    out = {}
    for name, rows in records.items():
        w = np.concatenate([r['weights'] for r in rows])
        fate = np.concatenate([r['fate'] for r in rows])
        t_in = np.concatenate([r['t_in'] for r in rows])/K.JULIAN_DAY
        v_in = np.concatenate([r['v_in'] for r in rows])
        hit = fate == 1
        days = np.arange(0, int(DESIGN['flight_days'])+1)
        cumulative = [float(w[hit & (t_in <= d)].sum()/w.sum()) for d in days]
        order = np.argsort(v_in[hit])
        cw = np.cumsum(w[hit][order])
        out[name] = dict(particles=int(len(w)), protected_sphere=float(w[hit].sum()/w.sum()),
                         near_earth=float(w[fate == 2].sum()/w.sum()), escaped=float(w[fate == 3].sum()/w.sum()),
                         still_flying=float(w[fate == 0].sum()/w.sum()), cumulative_by_day=cumulative,
                         delivered_by_lifetime={f'{tau:g}_days': float(np.sum(w[hit]*np.exp(-t_in[hit]/tau))/w.sum())
                                                for tau in DESIGN['lifetimes_days']},
                         entry_speed_median_m_s=float(v_in[hit][order][np.searchsorted(cw, cw[-1]/2)]) if hit.any() else None,
                         entry_energy_mean_MJ_kg=float(np.sum(w[hit]*.5*v_in[hit]**2)/w[hit].sum()/1e6) if hit.any() else None)
    coarse = first['fate']
    out['integrator_check'] = dict(step_s=DESIGN['step_s'], check_step_s=DESIGN['check_step_s'],
                                   particles=int(len(coarse)),
                                   protected_sphere_count=[int(np.sum(coarse == 1)), int(np.sum(fine == 1))],
                                   fate_changed=int(np.sum(coarse != fine)))
    out['propagation_cpu_s'] = cpu
    return out


def run_key():
    """Checkpoint key of the slow gas: its settings, the code that launches and flies it, and its inputs."""
    settings = {k: DESIGN[k] for k in ('core_sigma_kg_m2', 'annulus_sigma_kg_m2', 'protected_radii', 'cants_deg',
                                       'gas_temperature_K', 'molar_mass_kg_mol', 'launch_days', 'launch_every_days',
                                       'flight_days', 'particles_per_thruster', 'step_s', 'check_step_s',
                                       'near_earth_radii', 'escape_m', 'seed')}
    files = ''.join(digest(ROOT/f) for f in FILES if f != 'research/studies/solar_shield_array/exhaust_isolation.py')
    # This runner's analysis can change without touching the flights, so only its flight code enters the key.
    code = ''.join(inspect.getsource(f) for f in (unit, perpendicular, plume_axes, effusive, fly, slow_gas))
    return hashlib.sha256((json.dumps(settings, sort_keys=True)+files+code+digest(HOLDING)).encode()).hexdigest()[:16]


def main():
    wall, cpu0 = time.monotonic(), time.process_time()
    holding = json.loads(HOLDING.read_text())
    coefficients = np.array(holding['trajectory_search']['variable_distance']['coefficients_km'])
    zoned = {c['label']: c for c in json.loads(ZONED.read_text())['cases']}[DESIGN['case']]
    prop = SCENARIO['propulsion']
    jet = prop['efficiency']*zoned['mean_power_TW']*1e12
    # The year of the zoned product, every three hours.
    days = SCENARIO['primary_days']
    eph = Ephemeris(DEFAULT_KERNEL, SCENARIO['epoch_tdb'])
    year = eph.sample(np.linspace(0, days*K.JULIAN_DAY, int(days*8)+1))
    window = eph.sample(np.arange(-3600., (DESIGN['launch_days']+DESIGN['flight_days'])*K.JULIAN_DAY+7200., 1800.))
    eph.close()
    _, fields = zoned_aperture(year, geometry(year), coefficients, DESIGN['core_sigma_kg_m2'],
                               DESIGN['annulus_sigma_kg_m2'], CORE_DIVERTED, fields=True)
    direct = direct_path(fields)
    _, launch = zoned_aperture(window, geometry(window), coefficients, DESIGN['core_sigma_kg_m2'],
                               DESIGN['annulus_sigma_kg_m2'], CORE_DIVERTED, fields=True)
    env = CyclingEnvironment(window, .05, CORE_DIVERTED, 'filter_only')
    key = run_key()
    gas = slow_gas(launch, env, key)
    unionized = [s*zoned['propellant_kg_s'] for s in DESIGN['unionized_share']]
    implied = {}
    for name in DESIGN['molar_mass_kg_mol']:
        g = gas[name]
        implied[name] = {f'{tau}': dict(delivered_kg_s=[u*g['delivered_by_lifetime'][tau] for u in unionized],
                                        deposited_MW=[u*g['delivered_by_lifetime'][tau]*(g['entry_energy_mean_MJ_kg'] or 0.)
                                                      for u in unionized])
                         for tau in g['delivered_by_lifetime']}
    allowed = allowance(zoned['propellant_kg_s'], jet)
    magnets = magnetosphere(direct, implied, allowed['per_budget'], jet)
    out = dict(schema='terluna.research.exhaust-isolation/1',
               producer=dict(files={f: digest(ROOT/f) for f in FILES},
                             constants=constants_used([f for f in FILES if f.endswith('.py')]),
                             inputs={'results/holding.json': digest(HOLDING), 'results/zoned_aperture.json': digest(ZONED),
                                     'de440s.bsp': digest(DEFAULT_KERNEL), 'scenario.json': digest(HERE/'scenario.json')}),
               evidence=('The zoned held screen of zoned_aperture.json (26 g/m2 window, 5 g/m2 annulus) on the selected '
                         'moving trajectory, with DE440s geometry. The direct path samples every node daily over the '
                         'primary year and integrates two measured plume profiles over the cone of the protected sphere; '
                         'the plumes leave straight, with no fields, and a tenth of each plume is charge-exchanged fast '
                         'atoms on the same paths. The slow gas flies test particles from every node every two days '
                         'for a month, each for 60 days under point-mass gravity, without radiation pressure, '
                         'collisions or ionisation; lifetimes enter afterward as weights. The magnetosphere case '
                         'compares ion gyroradii at the September magnets\' magnetopause with their stand-off.'),
               reading_rule=('Compare each delivered share with the allowance share for the same loss budget. Without '
                             'the magnets the whole straight-line share counts and the solar wind\'s pickup of the '
                             'charged exhaust adds to it. With the September magnets the exhaust ions are held off, '
                             'leaving the fast-neutral share, the slow gas and the magnetosphere\'s leak, which is open; '
                             'with_magnets gives the share of the slow gas that must then be captured. Plume values are '
                             'read from published figures to about +-30%; the fast-neutral tenth comes from one Hall '
                             'thruster model.'),
               design=DESIGN, run_key=key, plumes=PLUMES, sources=SOURCES,
               allowance=allowed, direct_path=direct, slow_gas=gas,
               slow_gas_implied=dict(unionized_kg_s=unionized, by_species=implied), magnetosphere=magnets,
               open=('The charged exhaust, 85-95% of the flow, is a heavy-ion source about 10^5 times the solar wind\'s '
                     'own mass flux through the aperture. Without magnets the wind loads, slows and carries it '
                     'downstream over the Moon, energising picked-up oxygen to tens of keV; the loss response credits '
                     'each returning pickup ion with 1-10 sputtered molecules. The September magnets hold off exhaust '
                     'ions and picked-up oxygen, whose gyroradii at the magnetopause are small against the stand-off; '
                     'xenon picked up at solar-wind speed gyrates on a scale near the stand-off, and entry through '
                     'the cusps and by reconnection is open. No model here follows that plasma.'),
               resources=dict(cpu_s=time.process_time()-cpu0, wall_s=time.monotonic()-wall,
                              max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024, numerical_threads=1))
    OUT.write_text(json.dumps(out, indent=1)+'\n')
    print(json.dumps(dict(allowance=out['allowance'], direct=direct, slow_gas={k: {kk: vv for kk, vv in v.items()
                                                                                  if kk != 'cumulative_by_day'}
                                                                              if isinstance(v, dict) else v
                                                                              for k, v in gas.items()},
                          implied=out['slow_gas_implied'], magnetosphere=magnets, resources=out['resources']), indent=1))


if __name__ == '__main__':
    main()
