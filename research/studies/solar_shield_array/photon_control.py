"""Photon control of the kept ring bundle (step 4 of the integrated comparison).

The ring fleet keeps station with photon forces, so it releases no exhaust. This runner asks whether sunlight
can do that work along the twelve kept orbits of ring_keeping_run, read back from its checkpoints.

Translation. The filter's redirected band reflects specularly along each tile's normal (ring_bundle), so
tilting a tile's normal by up to a set angle and trimming the band's reflectivity by up to a set fraction move
its sail force within a bounded set. At each sample the keeping law's demand on each ring's representative is
compared with that set's reach in the demand's own direction (its support), with each tile's lit fraction
behind the bundle's shadows. Demand beyond the reach is velocity the law would have to move to other phases of
the orbit; the reach summed over the orbit bounds what any schedule could deliver. The kept run's demand is
inflated by its fixed keeping frame (relative_orbits.md), so the test is conservative.

Attitude. A flat tile facing radially holds its axis of largest inertia along the radius, the orientation that
gravity gradient makes unstable, and its 1.5 degree shingle tilt holds it off that equilibrium; mutual shadows
move its centre of pressure. Gravity and radiation torques come from tile_torque and the bundle's shadow
moments. Reflectivity trim across the halves of a tile supplies torque in proportion to the light on it;
where trim falls short, near the terminators and in eclipse, a momentum store takes the excess and trim
unloads it later.

Torques and sail forces are evaluated for the kept run's tiles (the seed's 50 g/m2 times its mass ratio);
gravity torque scales with areal mass and the sail force's reach inversely, so both are also stated for
5 and 26 g/m2 tiles.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.photon_control
"""
from __future__ import annotations

import hashlib
import json
import resource
import time
from pathlib import Path

import numpy as np

from shared import constants as K
from shared.provenance import constants_used
from protection.dynamics import ring_bundle as rb
from protection.dynamics.ephemeris import DEFAULT_KERNEL
from protection.dynamics.fleet import sun_points
from protection.dynamics.natural_pattern import retarded_sun
from protection.dynamics.optical import length, unit
from protection.dynamics.ring_keeping import ContinuousKeeping
from protection.dynamics.tile_torque import gravity_torque
from .ring_bundle_run import DESIGN as BUNDLE, SEED, environment, initial
from .ring_keeping_run import DESIGN as KEPT, FILES as KEPT_FILES, OUT as KEPT_OUT, kept

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE/'results/photon_control.json'
FILES = ['research/studies/solar_shield_array/photon_control.py', 'protection/dynamics/tile_torque.py'] + KEPT_FILES
DESIGN = dict(every=2, tilt_levels=[[1., .1], [2., .25], [5., .25]], azimuths=24, trim_levels=[.1, .25],
              unload_s=3600., areal_masses_kg_m2=[.005, .026], suns=8, rotation=.317)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sail(normals, sun_dir, irradiance, lit, fraction):
    """Specular sail force (N) on each tile's clear area along its normal, scaled by its lit fraction."""
    cosine = np.sum(sun_dir*normals, axis=-1)
    pressure = -2*fraction/K.SPEED_OF_LIGHT*irradiance*cosine*np.abs(cosine)
    return (pressure*lit*rb.CLEAR**2)[..., None]*normals


def reach(frames, sun_dir, irradiance, lit, fraction, mass, demand, tilt_deg, trim):
    """Largest acceleration change along each demand that tilt and trim can give (m/s^2), per tile."""
    n0, e1, e2 = frames[:, :, 2], frames[:, :, 0], frames[:, :, 1]
    base = sail(n0, sun_dir, irradiance, lit, fraction)
    size = length(demand)[:, None]
    direction = np.divide(demand, size, out=np.zeros_like(demand), where=size > 0)
    best = np.zeros(len(n0))
    delta = np.radians(tilt_deg)
    psi = 2*np.pi*np.arange(DESIGN['azimuths'])/DESIGN['azimuths']
    candidates = [n0]+[np.cos(delta)*n0+np.sin(delta)*(np.cos(p)*e1+np.sin(p)*e2) for p in psi]
    for normal in candidates:
        for f in (-trim, 0., trim):
            change = (sail(normal, sun_dir, irradiance, lit, fraction*(1+f))-base)/mass[:, None]
            best = np.maximum(best, np.sum(change*direction, axis=1))
    return best


def torques(env, t, everything, count, tilt, mass, fraction):
    """Radiation torque from the shadowed area's moments and gravity-gradient torque on the first count tiles,
    in the bundle's common frame (N m), with each tile's pressure scale at its incidence."""
    sample = env.at(t)
    q = everything[:, :3]
    common = rb.centroid_frame(everything, tilt)
    bodies = rb.eclipse_bodies(q, sample)
    radiation = np.zeros((len(q), 3))
    for source in sun_points(retarded_sun(sample), DESIGN['suns'], DESIGN['rotation']):
        area, mx, my = rb.visible_moments(env, t, q, common, source, sample, bodies).T
        ray = source-q
        cosine = unit(ray)@common[:, 2]
        density = K.SOLAR_CONSTANT*(K.AU/length(ray))**2/DESIGN['suns']*np.abs(cosine)
        pressure = -2*fraction/K.SPEED_OF_LIGHT*density*cosine
        radiation += pressure[:, None]*np.c_[my, -mx, np.zeros(len(q))]
    gravity = gravity_torque(env, t, q[:count], common)*mass[:, None]
    return radiation[:count], gravity


def main():
    wall, cpu0 = time.monotonic(), time.process_time()
    env = environment(KEPT['orbits']*2.0+2)
    states, ring, slot, ratio, built = initial(env, KEPT['radius_step_m'])
    tilt = np.radians(BUNDLE['tilt_deg'])
    reps = states[slot == (BUNDLE['per_ring']-1)//2]
    lags = BUNDLE['pitch_m']/np.linalg.norm(reps[:, 3:], axis=1)
    model = rb.ContinuousRings(env, tilt, np.full(len(reps), ratio), lags, reach=KEPT['copy_reach'])
    origin = unit(env.at(0.)['positions']['sun'])
    t, s, _, _, period, _, key, _ = kept(env, reps, model, origin)
    mid = len(reps)//2
    keeper = ContinuousKeeping(model, reps, mid, origin, KEPT['tau_s'], KEPT['limit_m_s2'])
    mass = model.forces.mass[:len(reps)]
    fraction = env.central_fraction
    edges = np.isin(np.arange(len(reps)), [0, len(reps)-1])
    pick = np.arange(0, len(t), DESIGN['every'])
    dt = np.gradient(t[pick])
    demand, cos_sun, lit_all = [], [], []
    support = {f'{d:g}_deg_{f:g}': [] for d, f in DESIGN['tilt_levels']}
    rad, grav, scale = [], [], []
    sail_force = []
    for i in pick:
        state = s[i]
        everything = model.copies(state)[0]
        force, first, bare = model.forces.light(t[i], everything)
        sail_force.append(length(force[:len(reps)]))
        lit = np.divide(first[:len(reps)], bare[:len(reps)], out=np.zeros(len(reps)), where=bare[:len(reps)] > 0)
        control = keeper.control(t[i], state)
        sun = retarded_sun(env.at(t[i]))
        sun_dir = unit(sun-state[:, :3])
        irradiance = K.SOLAR_CONSTANT*(K.AU/length(sun-state[:, :3]))**2
        frames = rb.tile_frames(state, tilt)
        for d, f in DESIGN['tilt_levels']:
            support[f'{d:g}_deg_{f:g}'].append(reach(frames, sun_dir, irradiance, lit, fraction, mass, control, d, f))
        r, g = torques(env, t[i], everything, len(reps), tilt, mass, fraction)
        cosine = np.sum(sun_dir*frames[:, :, 2], axis=1)
        demand.append(control); cos_sun.append(cosine); lit_all.append(lit); rad.append(r); grav.append(g)
        # Torque that a reflectivity trim of 100% across the tile's halves would give at this incidence.
        scale.append(2*fraction/K.SPEED_OF_LIGHT*irradiance*cosine**2*lit*rb.CLEAR**3/8)
    demand, cos_sun, lit_all = np.array(demand), np.array(cos_sun), np.array(lit_all)
    rad, grav, scale, sail_force = np.array(rad), np.array(grav), np.array(scale), np.array(sail_force)
    need = np.linalg.norm(demand, axis=2)
    keepers = np.arange(len(reps)) != mid

    def share(level):
        reachable = np.array(support[level])
        delivered = np.minimum(need, reachable)
        rows = {}
        for name, m in (('interior', keepers & ~edges), ('edges', edges)):
            total = float(np.sum(need[:, m]*dt[:, None]))
            rows[name] = dict(demand_m_s=total, instant_share=float(np.sum(delivered[:, m]*dt[:, None]))/total,
                              reach_over_demand=float(np.sum(reachable[:, m]*dt[:, None]))/total)
        return rows
    translation = {level: share(level) for level in support}
    # Smallest tilt that reaches each demand, among the tested levels.
    order = list(support)
    reached = np.stack([np.array(support[level]) >= need for level in order])
    first_level = np.where(reached.any(axis=0), reached.argmax(axis=0), len(order))
    level_share = {name: float(np.mean(first_level[:, keepers] == k)) for k, name in enumerate(order+['beyond'])}
    # Attitude: in-plane torques (u, v axes of the common frame) against trim capacity, with a momentum store.
    torque = (rad+grav)[:, :, :2]
    attitude = {}
    for trim in DESIGN['trim_levels']:
        cap = trim*scale
        stored = np.zeros((len(reps), 2))
        peak = np.zeros(len(reps))
        for j in range(len(pick)):
            want = -torque[j]-stored/DESIGN['unload_s']
            applied = np.clip(want, -cap[j][:, None], cap[j][:, None])
            stored = stored+(torque[j]+applied)*dt[j]
            peak = np.maximum(peak, np.linalg.norm(stored, axis=1))
        attitude[f'trim_{trim:g}'] = dict(momentum_store_N_m_s_max=float(peak.max()),
                                          momentum_store_N_m_s_median=float(np.median(peak)))
    lit_samples = scale > 0
    needed_trim = np.max(np.abs(torque), axis=2)/np.where(scale > 0, scale, np.inf)
    # Shadows move the centre of pressure; a moving mass that follows it cancels the radiation torque while lit.
    strong = sail_force > .05*np.max(sail_force)
    offset = np.linalg.norm(rad[:, :, :2], axis=2)[strong]/sail_force[strong]
    orbit = np.minimum((t[pick]/period).astype(int), KEPT['orbits']-1)
    interior = keepers & ~edges
    per_axis = {}
    for axis, name in ((0, 'roll_about_along_track'), (1, 'pitch_about_orbit_normal')):
        mean = np.array([[np.average(rad[orbit == k, j, axis], weights=dt[orbit == k]) for j in range(len(reps))]
                         for k in range(KEPT['orbits'])])
        per_axis[name] = dict(
            radiation_rms_N_m=float(np.sqrt(np.mean(rad[:, interior, axis]**2))),
            radiation_orbit_mean_N_m_median=float(np.median(np.abs(mean[:, interior]))),
            gravity_mean_N_m=float(np.mean(grav[:, interior, axis])))
    sigma = env.sigma*np.asarray(ratio)
    scaled = {}
    for areal in DESIGN['areal_masses_kg_m2']:
        g = grav*areal/sigma
        tt = (rad+g)[:, :, :2]
        trim = np.max(np.abs(tt), axis=2)/np.where(scale > 0, scale, np.inf)
        scaled[f'{areal*1e3:g}_g_m2'] = dict(
            gravity_torque_N_m_median=float(np.median(np.linalg.norm(g, axis=2))),
            trim_needed_lit_median=float(np.median(trim[lit_samples])),
            trim_needed_lit_95=float(np.percentile(trim[lit_samples], 95)),
            sail_reach_factor=float(sigma/areal))
    out = dict(schema='terluna.research.photon-control/1',
               producer=dict(files={f: digest(ROOT/f) for f in FILES},
                             constants=constants_used([f for f in FILES if f.endswith('.py')]),
                             inputs={'results/ring_keeping.json': digest(KEPT_OUT), 'results/joint_seed.json': digest(SEED),
                                     'de440s.bsp': digest(DEFAULT_KERNEL)}),
               evidence=('The twelve kept orbits of ring_keeping_run (run key in kept_run_key), every second '
                         'five-minute sample. Translation: the keeping law\'s acceleration on each representative '
                         'against the reach of a specular sail force on the tile\'s clear area, its lit fraction '
                         'behind the bundle\'s shadows, with the normal tilted up to the stated angle in 24 azimuths '
                         'and the redirected band trimmed by up to the stated fraction. Attitude: radiation torque '
                         'from the visible area\'s first moments (8 solar points) and gravity-gradient torque from a '
                         '3x3 quadrature of the tile, in the bundle\'s common frame; trim torque is the stated '
                         'fraction of the band\'s pressure on half the clear square at a quarter-side lever.'),
               reading_rule=('Translation is feasible where the instant share is near one; where the reach summed over '
                             'the orbit exceeds the demand, a schedule that moves the shortfall to lit phases can '
                             'deliver it. Attitude is held by trim alone where the stored momentum stays small; the '
                             'store is what a wheel or a moving mass must hold. The kept run\'s tiles carry '
                             f"{sigma*1e3:.1f} g/m2; lighter tiles gain sail reach and lose gravity torque in proportion."),
               design=DESIGN, kept_run_key=key, period_h=period/3600., tile_areal_mass_kg_m2=float(sigma),
               samples=int(len(pick)),
               sail_acceleration_normal_incidence_m_s2=float(2*fraction*K.SOLAR_CONSTANT/K.SPEED_OF_LIGHT*rb.CLEAR**2
                                                             /mass[0]),
               demand=dict(peak_m_s2=float(need.max()), median_m_s2=float(np.median(need[:, keepers])),
                           lit_fraction_median=float(np.median(lit_all)),
                           sunlit_share=float(np.mean(lit_all[:, keepers] > 0))),
               translation=translation, first_tilt_level_that_reaches=level_share,
               attitude=dict(radiation_torque_N_m_median=float(np.median(np.linalg.norm(rad, axis=2))),
                             radiation_torque_N_m_95=float(np.percentile(np.linalg.norm(rad, axis=2), 95)),
                             gravity_torque_N_m_median=float(np.median(np.linalg.norm(grav, axis=2))),
                             by_axis=per_axis,
                             pressure_centre_offset_m=dict(median=float(np.median(offset)),
                                                           p95=float(np.percentile(offset, 95)),
                                                           max=float(offset.max())),
                             trim_needed_lit_median=float(np.median(needed_trim[lit_samples])),
                             trim_needed_lit_95=float(np.percentile(needed_trim[lit_samples], 95)),
                             by_trim=attitude, by_areal_mass=scaled),
               resources=dict(cpu_s=time.process_time()-cpu0, wall_s=time.monotonic()-wall,
                              max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024, numerical_threads=1))
    OUT.write_text(json.dumps(out, indent=1)+'\n')
    print(json.dumps({k: v for k, v in out.items() if k not in ('producer', 'evidence', 'reading_rule')}, indent=1))


if __name__ == '__main__':
    main()
