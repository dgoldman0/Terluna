"""Validation of the kept ring bundle with independently propagated tiles and receiver rays (item 4 of the
narrowing).

The kept run (ring_keeping_run) flies one representative per ring and stands the rest of each ring and its
neighbours in as two-body copies, and it measures overlaps as orthographic projections along the Sun's centre. This
runner removes both stand-ins on a patch of the kept bundle: 31 tiles on each of its 19 rings at the kept run's 1 km
radius step, every tile propagated with the full forces (DE440s point masses with finite-square quadrature, the
filter's sail force on its own radial-facing normal, finite-Sun mutual shadows among all 589 tiles, Moon and Earth
eclipses).

Keeping. Each ring's middle tile is kept as the kept run keeps its representative (ContinuousKeeping, thirty-minute
time constant, 0.001 m/s^2 limit). Every other tile of the ring starts on its slot, the middle tile's state carried
along its two-body orbit by the slot's lag, the position the copies assumed, and feels the middle tile's acceleration
and a critically damped pull back to that slot. The patch's two ends and its edge rings get more light than a
complete ring's tiles would; the pull holds them as well, and the checks read the interior.

Checks every five minutes: the smallest surface distance between tiles, over the patch and over its interior (rings
1-17, tiles 6-24); inside the service arc, the orthographic overlaps along each ring and between rings, as
ring_bundle_run measures them. Receiver rays: at every third service sample, rays from the plane through the Moon
normal to the Sun to points spread uniformly over the solar disk, drawn so that each crosses the bundle at least two
rings inside the patch's edges and, as the rings slide apart, a tile pitch inside the ends of its own and both
neighbouring rings' interior tiles (the fifth to the twenty-seventh), are tested against every tile's own square. A
covered ray meets at least one tile; its margin is how far inside the edge of the tile it meets most deeply it
passes. Light time between the bundle and the Moon moves the whole pattern alike and is left out. The keeping's
demand on each tile is set against the sail force of photon_control's tiles.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.bundle_validation
"""
from __future__ import annotations

import hashlib
import inspect
import json
import resource
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from shared import constants as K
from shared.provenance import constants_used
from protection.dynamics import ring_bundle as rb
from protection.dynamics.ephemeris import DEFAULT_KERNEL
from protection.dynamics.fleet import sail_basis
from protection.dynamics.natural_pattern import retarded_sun
from protection.dynamics.optical import unit
from protection.dynamics.ring_keeping import ContinuousKeeping, geometric_elements
from .ring_bundle_run import DESIGN as BUNDLE, FILES as BUNDLE_FILES, SEED, environment, initial, overlaps, sun_angle
from .ring_keeping_run import DESIGN as KEPT

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RUN = ROOT/'research/runs/solar_shield_array/bundle_validation'
OUT = HERE/'results/bundle_validation.json'
PHOTON = HERE/'results/photon_control.json'
FILES = ['research/studies/solar_shield_array/bundle_validation.py', 'protection/dynamics/ring_keeping.py'] + BUNDLE_FILES
DYNAMICS = [f for f in FILES if not f.endswith(('_run.py', 'bundle_validation.py'))]
DESIGN = dict(orbits=4, tau_s=KEPT['tau_s'], limit_m_s2=KEPT['limit_m_s2'], radius_step_m=KEPT['radius_step_m'],
              sample_s=300., rtol=1e-10, atol=1e-3, max_step_s=600., interior_ring_margin=1, interior_slot_margin=6,
              ray_every=3, rays=8192, ray_ring_margin=2, ray_slot_margin=4, ray_chunk=512, seed=20261007)
HOUR = 3600.


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class KeptPatch:
    """Every tile propagated; each ring's middle tile kept to its neighbour ring, the others pulled to their slots."""
    def __init__(self, env, tilt, ratio, states, ring, slot, origin, tau, limit):
        self.forces = rb.RingBundleForces(env, tilt, np.full(len(states), ratio))
        self.ring, self.tau, self.limit = ring, tau, limit
        mid = (BUNDLE['per_ring']-1)//2
        self.middles = np.array([np.nonzero((ring == k) & (slot == mid))[0][0] for k in np.unique(ring)])
        speed = np.linalg.norm(states[self.middles, 3:], axis=1)
        self.lags = (slot-mid)*BUNDLE['pitch_m']/speed[ring]
        self.keeper = ContinuousKeeping(SimpleNamespace(env=env), states[self.middles], len(self.middles)//2, origin,
                                        tau, limit)

    def slots(self, states):
        """Each tile's slot: its ring's middle tile carried along its two-body orbit by the tile's lag."""
        return rb.kepler_shift(states[self.middles][self.ring], self.lags)

    def control(self, t, states):
        target = self.slots(states)
        pull = -2/self.tau*(states[:, 3:]-target[:, 3:])-(states[:, :3]-target[:, :3])/self.tau**2
        pull[self.middles] = 0.
        a = self.keeper.control(t, states[self.middles])[self.ring]+pull
        norm = np.linalg.norm(a, axis=1)
        return a*np.minimum(1., self.limit/np.maximum(norm, 1e-30))[:, None]

    def acceleration(self, t, states):
        return self.forces.acceleration(t, states)+self.control(t, states)


def run_key(states, origin):
    settings = {k: DESIGN[k] for k in ('tau_s', 'limit_m_s2', 'radius_step_m', 'sample_s', 'rtol', 'atol', 'max_step_s')}
    return hashlib.sha256(json.dumps(settings, sort_keys=True).encode()+states.tobytes()+origin.tobytes()+
                          inspect.getsource(KeptPatch).encode()+''.join(digest(ROOT/f) for f in DYNAMICS).encode()
                          ).hexdigest()[:16]


def propagate(model, states, origin, period):
    """Orbit by orbit with checkpoints; sample times and states."""
    RUN.mkdir(parents=True, exist_ok=True)
    key = run_key(states, origin)
    times, samples, current, t0, cpu = [], [], states, 0., 0.
    for k in range(DESIGN['orbits']):
        path = RUN/f'orbit_{k}_{key}.npz'
        t1 = (k+1)*period
        if not path.exists():
            start = time.process_time()
            at, sol = rb.propagate(model.acceleration, current, t0, t1, DESIGN['rtol'], DESIGN['atol'],
                                   DESIGN['max_step_s'])
            grid = np.arange(t0, t1, DESIGN['sample_s'])
            grid = np.r_[grid, t1] if grid[-1] < t1 else grid
            tmp = path.with_name(path.stem+'.tmp.npz')
            np.savez(tmp, t=grid, s=np.stack([at(x) for x in grid]), cpu=time.process_time()-start, nfev=sol.nfev)
            tmp.rename(path)
            print(f'orbit {k}: {time.process_time()-start:.0f} s', flush=True)
        z = np.load(path)
        times.append(z['t'] if k == 0 else z['t'][1:]); samples.append(z['s'] if k == 0 else z['s'][1:])
        cpu += float(z['cpu']); current, t0 = z['s'][-1], t1
    return np.concatenate(times), np.concatenate(samples), cpu, key


def rays(states, tilt, sample, ring, slot, rng, count):
    """Rays (origins on the plane through the Moon normal to the Sun, unit directions toward points spread uniformly
    over the solar disk) that cross the bundle inside the patch's interior.

    Rings slide past one another, so the patch's ends shear apart; a ray in ring k's band is drawn only where rings
    k-1, k and k+1 all still carry interior tiles, a tile pitch inside the nearest of their ends."""
    sun = retarded_sun(sample)
    s_hat = unit(sun)
    mid = (BUNDLE['per_ring']-1)//2
    centre = states[(ring == ring.max()//2) & (slot == mid)][0]
    frame = rb.tile_frames(centre, tilt)
    q0, u, v = centre[:3], frame[:, 0], frame[:, 1]
    rel = states[:, :3]-q0
    along, across = rel@u, rel@v
    lo_slot, hi_slot = DESIGN['ray_slot_margin'], BUNDLE['per_ring']-1-DESIGN['ray_slot_margin']
    level = np.array([np.median(across[ring == k]) for k in range(ring.max()+1)])
    span = np.sort([[along[(ring == k) & (slot == lo_slot)][0], along[(ring == k) & (slot == hi_slot)][0]]
                    for k in range(ring.max()+1)], axis=1)
    bands = []
    for k in range(DESIGN['ray_ring_margin'], ring.max()-DESIGN['ray_ring_margin']+1):
        lo = max(span[j, 0] for j in (k-1, k, k+1))+BUNDLE['pitch_m']
        hi = min(span[j, 1] for j in (k-1, k, k+1))-BUNDLE['pitch_m']
        y0, y1 = sorted(((level[k-1]+level[k])/2, (level[k]+level[k+1])/2))
        if hi > lo:
            bands.append((lo, hi, y0, y1))
    bands = np.array(bands)
    area = (bands[:, 1]-bands[:, 0])*(bands[:, 3]-bands[:, 2])
    pick = bands[rng.choice(len(bands), size=count, p=area/area.sum())]
    x, y = rng.uniform(pick[:, 0], pick[:, 1]), rng.uniform(pick[:, 2], pick[:, 3])
    point = q0+x[:, None]*u+y[:, None]*v
    _, b, c = sail_basis(sun)
    radius = K.SUN_RADIUS*np.sqrt(rng.uniform(0., 1., count))
    angle = rng.uniform(0., 2*np.pi, count)
    source = sun+radius[:, None]*(np.cos(angle)[:, None]*b+np.sin(angle)[:, None]*c)
    d = unit(source-point)
    lam = np.sum(point*s_hat, axis=1)/np.sum(d*s_hat, axis=1)
    return point-lam[:, None]*d, d, len(bands), float(area.sum())


def coverage(origins, directions, states, tilt):
    """Per ray: whether it meets a tile's clear square, how many it meets, and its deepest margin inside an edge."""
    frames = rb.tile_frames(states, tilt)
    q, u, v, n = states[:, :3], frames[:, :, 0], frames[:, :, 1], frames[:, :, 2]
    half = rb.CLEAR/2
    hits, margin = [], []
    for k in range(0, len(origins), DESIGN['ray_chunk']):
        o, d = origins[k:k+DESIGN['ray_chunk']], directions[k:k+DESIGN['ray_chunk']]
        denom = d@n.T
        tau = np.einsum('rtj,tj->rt', q[None]-o[:, None], n)/denom
        point = o[:, None]+tau[..., None]*d[:, None]-q[None]
        x, y = np.einsum('rtj,tj->rt', point, u), np.einsum('rtj,tj->rt', point, v)
        inside = (tau > 0) & (np.abs(x) <= half) & (np.abs(y) <= half)
        depth = np.where(inside, np.minimum(half-np.abs(x), half-np.abs(y)), -np.inf)
        hits.append(inside.sum(axis=1)); margin.append(depth.max(axis=1))
    return np.concatenate(hits), np.concatenate(margin)


def main():
    wall, cpu0 = time.monotonic(), time.process_time()
    env = environment(KEPT['orbits']*2.0+2)
    states, ring, slot, ratio, built = initial(env, DESIGN['radius_step_m'])
    tilt = np.radians(BUNDLE['tilt_deg'])
    origin = unit(env.at(0.)['positions']['sun'])
    model = KeptPatch(env, tilt, ratio, states, ring, slot, origin, DESIGN['tau_s'], DESIGN['limit_m_s2'])
    # Tiles start on their slots, the copies' positions; the middle tiles are the kept run's representatives.
    states = model.slots(states)
    mids = states[model.middles]
    _, a0, n0, _ = geometric_elements(mids[len(mids)//2], mids, origin)
    period = 2*np.pi/n0
    t, s, cpu_run, key = propagate(model, states, origin, period)
    interior = ((ring >= DESIGN['interior_ring_margin']) & (ring <= ring.max()-DESIGN['interior_ring_margin']) &
                (slot >= DESIGN['interior_slot_margin']) & (slot < BUNDLE['per_ring']-DESIGN['interior_slot_margin']))
    rng = np.random.default_rng(DESIGN['seed'])
    clear, inner, service, demand, offset, ray_rows = [], [], [], [], [], []
    for i, (tt, ss) in enumerate(zip(t, s)):
        d, pair = rb.own_attitude_clearance(ss, tilt)
        clear.append((d, int(ring[pair[0]]), int(slot[pair[0]]), int(ring[pair[1]]), int(slot[pair[1]])))
        d, _ = rb.own_attitude_clearance(ss[interior], tilt)
        inner.append(d)
        a = np.linalg.norm(model.control(tt, ss), axis=1)
        demand.append(a)
        offset.append(np.linalg.norm(ss[:, :3]-model.slots(ss)[:, :3], axis=1))
        angle = sun_angle(env, tt, ss)
        if abs(angle) <= np.radians(BUNDLE['service_half_angle_deg']):
            sun = env.at(tt)['positions']['sun']-ss[:, :3].mean(axis=0)
            along, across = overlaps(ss, ring, slot, unit(sun), tilt, interior)
            service.append((tt, along, across))
            if len(service) % DESIGN['ray_every'] == 1:
                o, dirs, bands, area = rays(ss, tilt, env.at(tt), ring, slot, rng, DESIGN['rays'])
                hits, margin = coverage(o, dirs, ss, tilt)
                ray_rows.append(dict(hour=float(tt/HOUR), phase_deg=float(np.degrees(angle)), rays=int(len(hits)),
                                     bands=bands, sampled_area_km2=area/1e6,
                                     uncovered=int(np.count_nonzero(hits == 0)),
                                     margin_min_m=float(margin[hits > 0].min()) if np.any(hits > 0) else None,
                                     layers=[int(np.count_nonzero(hits == j)) for j in range(4)]))
    clear, inner, service = np.array(clear), np.array(inner), np.array(service)
    demand, offset = np.array(demand), np.array(offset)
    photon = json.loads(PHOTON.read_text())
    sail = photon['sail_acceleration_normal_incidence_m_s2']
    kept_median = photon['demand']['median_m_s2']

    def first_below(series):
        below = np.nonzero(series < 150.)[0]
        return float(t[below[0]]/HOUR) if len(below) else None
    k = int(np.argmin(clear[:, 0]))
    middle = np.zeros(len(states), bool); middle[model.middles] = True
    out = dict(schema='terluna.research.bundle-validation/1',
               producer=dict(files={f: digest(ROOT/f) for f in FILES},
                             constants=constants_used([f for f in FILES if f.endswith('.py')]),
                             inputs={'results/joint_seed.json': digest(SEED), 'results/photon_control.json': digest(PHOTON),
                                     'de440s.bsp': digest(DEFAULT_KERNEL)}),
               evidence=('A patch of the kept bundle, 19 rings of 31 tiles at a 1 km radius step, from 2026-10-04 TDB '
                         'for four orbits with every tile propagated (DE440s point masses with finite-square quadrature, '
                         'the filter sail force on each tile\'s radial-facing normal, finite-Sun shadows among all '
                         'tiles in the patch centroid\'s frame, Moon and Earth eclipses), each ring\'s middle tile kept '
                         'as in ring_keeping_run and every other tile pulled to its two-body slot. Clearance every five '
                         'minutes; orthographic overlaps in the service arc; exact ray-square tests of uniformly drawn '
                         'receiver rays toward the finite solar disk through the patch interior.'),
               reading_rule=('The bundle passes with independent tiles when the interior keeps 150 m, the overlaps stay '
                             'positive and no interior ray goes uncovered. The patch ends and edge rings carry light a '
                             'complete ring would not; the whole-patch clearance and the demand at the ends show what '
                             'that costs. Keeping demand is per tile, against the sail acceleration of the kept run\'s '
                             'tiles at normal incidence.'),
               design=DESIGN, run_key=key, period_h=period/HOUR, built=built, propagation_cpu_s=cpu_run,
               clearance=dict(min_m=float(clear[k, 0]), at_h=float(t[k]/HOUR), pair_ring_slot=clear[k, 1:].astype(int).tolist(),
                              first_below_150_m_h=first_below(clear[:, 0]),
                              samples_below_150_m=int(np.count_nonzero(clear[:, 0] < 150.)), samples=len(clear),
                              series_h=(t/HOUR).round(3).tolist()[::6], series_m=clear[:, 0].round(1).tolist()[::6]),
               interior_clearance=dict(min_m=float(inner.min()), at_h=float(t[int(np.argmin(inner))]/HOUR),
                                       first_below_150_m_h=first_below(inner),
                                       samples_below_150_m=int(np.count_nonzero(inner < 150.)),
                                       series_m=inner.round(1).tolist()[::6]),
               service=dict(samples=len(service), along_ring_overlap_min_m=float(service[:, 1].min()),
                            adjacent_ring_overlap_min_m=float(service[:, 2].min())),
               rays=dict(samples=len(ray_rows), tested=int(sum(r['rays'] for r in ray_rows)),
                         uncovered=int(sum(r['uncovered'] for r in ray_rows)),
                         margin_min_m=float(min(r['margin_min_m'] for r in ray_rows if r['margin_min_m'] is not None)),
                         layers_share=(np.sum([r['layers'] for r in ray_rows], axis=0)/
                                       sum(r['rays'] for r in ray_rows)).round(4).tolist(),
                         by_sample=ray_rows),
               keeping=dict(sail_acceleration_normal_incidence_m_s2=sail, kept_run_representative_median_m_s2=kept_median,
                            middle_tiles_m_s2=dict(median=float(np.median(demand[:, middle])),
                                                   p95=float(np.percentile(demand[:, middle], 95))),
                            interior_tiles_m_s2=dict(median=float(np.median(demand[:, interior & ~middle])),
                                                     p95=float(np.percentile(demand[:, interior & ~middle], 95)),
                                                     max=float(demand[:, interior & ~middle].max())),
                            patch_ends_and_edges_m_s2=dict(median=float(np.median(demand[:, ~interior])),
                                                           p95=float(np.percentile(demand[:, ~interior], 95)),
                                                           max=float(demand[:, ~interior].max())),
                            slot_offset_m=dict(interior_max=float(offset[:, interior].max()),
                                               interior_median=float(np.median(offset[:, interior])),
                                               ends_and_edges_max=float(offset[:, ~interior].max()))),
               resources=dict(cpu_s=time.process_time()-cpu0, wall_s=time.monotonic()-wall,
                              max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024, numerical_threads=1))
    OUT.write_text(json.dumps(out, indent=1)+'\n')
    print(json.dumps({k: v for k, v in out.items() if k in ('clearance', 'interior_clearance', 'service', 'keeping')},
                     indent=1)[:3000])
    print(json.dumps({k: v for k, v in out['rays'].items() if k != 'by_sample'}))


if __name__ == '__main__':
    main()
