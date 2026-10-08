"""Self-contained coupled experiment using the committed loaded epoch seed.

Restoration trajectories may continue through contact solely to define an
optimization residual. Only event-stopped, independently audited trajectories
are physical executions. No supply equipment or force model is substituted.
"""
import hashlib
import json
import resource
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.spatial.transform import Rotation, RotationSpline

from shared.provenance import constants_used, constants_changed
from protection.dynamics.service_search import free_command, smooth_step
from protection.dynamics.eclipse_parallel import EclipseParallelPattern, load
from protection.dynamics.attitude_load import rigid_square_load
from protection.dynamics.pattern_return import smooth_arc
from protection.dynamics.square_distance import closest_squares
from protection.dynamics.optical import length
from .fleet_run import environment, digest
from .joint_common import spline, state_at

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
RUN = ROOT / 'research/runs/solar_shield_array/coupled_optimization'
SPECIFIC = 30000 / (1.4 * np.cos(np.pi / 4))
DIMENSION = 8


def inputs(extra=()):
    seed = json.loads((HERE / 'results/joint_seed.json').read_text())
    if constants_changed(seed['producer']['constants']):
        raise ValueError('Seed constants changed')
    for name, sha in seed['producer']['source_hashes'].items():
        if digest(ROOT / name) != sha:
            raise ValueError('Seed source changed: ' + name)
    arrays = {k: np.asarray(v) for k, v in seed['arrays'].items()}
    for name, a in arrays.items():
        if hashlib.sha256(np.asarray(a, dtype='<f8').tobytes()).hexdigest() != seed['float64_little_endian_hashes'][name]:
            raise ValueError('Seed array hash mismatch: ' + name)
    source = dict(seed['producer']['source_hashes'])
    files = ['research/studies/solar_shield_array/coupled_common.py',
             'protection/dynamics/eclipse_parallel.py', 'protection/dynamics/eclipse_parallel.cpp',
             'protection/dynamics/fast_parallel.py', 'protection/dynamics/fast_parallel.cpp',
             'protection/dynamics/service_search.py', 'protection/dynamics/square_distance.py', *extra]
    source.update({p: digest(ROOT / p) for p in files})
    parent = 'natural_free_slow_coarse.json'
    arrival = json.loads((HERE / 'results' / parent).read_text())['arrival_s']
    producer = dict(source_hashes=source, constants=constants_used(source), inputs={
        'joint_seed.json': digest(HERE / 'results/joint_seed.json'),
        parent: digest(HERE / 'results' / parent)})
    return arrays, arrival, producer


def write(name, out, cpu_start, wall_start):
    RUN.mkdir(parents=True, exist_ok=True)
    out['resources'] = dict(cpu_s=time.process_time() - cpu_start,
        wall_s=time.monotonic() - wall_start,
        max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
        numerical_threads=1, new_raw_MiB=sum(p.stat().st_size for p in RUN.rglob('*') if p.is_file()) / 1024**2)
    (HERE / 'results' / name).write_text(json.dumps(out, indent=2, allow_nan=False) + '\n')


def save(name, **arrays):
    RUN.mkdir(parents=True, exist_ok=True)
    path = RUN / (name + '.npz')
    np.savez_compressed(path, **arrays)
    return dict(raw_path=str(path.relative_to(ROOT)), raw_sha256=digest(path))


def read_raw(row):
    path = ROOT / row['raw_path']
    if digest(path) != row['raw_sha256']:
        raise ValueError('Raw hash mismatch')
    return dict(np.load(path))


class Experiment:
    def __init__(self):
        self.arrays, self.arrival, self.producer = inputs()
        self.env = environment(3.)
        self.end = self.arrival + 43200.
        self.ratio = self.arrays['mass_ratio']
        self.mass = 5e6 * self.ratio
        self.capacity = (self.ratio - 1.2) * 5e6 * 300
        rows, cols = np.indices((19, 19))
        patterns = [(-1.)**rows.ravel(), (-1.)**cols.ravel(),
                    (-1.)**(rows + cols).ravel(), (rows.ravel() - 9.) / 9.]
        for i, j in [(325, 347), (338, 341)]:
            v = np.zeros(361); v[i] = 1.; v[j] = -1.; patterns.append(v)
        # Zero mass-weighted mean: these are relative-motion corrections.
        patterns = np.array(patterns).T
        patterns -= np.average(patterns, axis=0, weights=self.mass)
        patterns /= np.max(abs(patterns), axis=0)
        normal = self.command(np.zeros(DIMENSION))(13.8 * 3600).as_matrix()[:, 2]
        self.basis = .1 * patterns[:, None, :] * normal[None, :, None]

    def command(self, x):
        # Common attitude only. The subsequent departure is part of the target.
        duration = 6. + .5 * x[6]
        first = free_command(self.env, self.end, arrival=self.arrival,
                             turn_start=6., turn_hours=duration)
        times = np.unique(np.r_[first.times, self.arrival + 21600., self.end])
        angle = -np.pi / 2 * smooth_step((times - self.arrival - 21600.) / 21600.)
        frames = first(times).as_matrix() @ Rotation.from_rotvec(np.c_[angle, np.zeros((len(times), 2))]).as_matrix()
        return RotationSpline(times, Rotation.from_matrix(frames))

    def controls(self, x):
        # x[:6] are dimensionless; unit coefficient = 0.1 m/s maximum mode.
        u = np.einsum('nic,c->ni', self.basis, x[:6])
        windows = np.array([[21600., (12. + .5 * x[7]) * 3600.]])
        return u[:, None, :], windows

    def thrust(self, x, times):
        u, windows = self.controls(x)
        return np.array([np.einsum('k,nki->ni', smooth_arc(t, windows), u) for t in np.atleast_1d(times)])

    def evaluate(self, x, initial, start, end, fine=False, stop=False, step=120., account=True):
        cpu = time.process_time(); wall = time.monotonic()
        x = np.asarray(x); cmd = self.command(x)
        suns, grid = (16, 32) if fine else (8, 16)
        model = EclipseParallelPattern(self.env, cmd, self.ratio, suns=suns, grid=grid)
        u, windows = self.controls(x)
        def rhs(t, flat):
            y = flat.reshape(361, 6)
            acceleration = model.acceleration(t, y) + np.einsum('k,nki->ni', smooth_arc(t, windows), u)
            return np.c_[y[:, 3:], acceleration].ravel()
        def event(t, flat):
            return closest_squares(flat.reshape(361, 6)[:, :3], cmd(t).as_matrix())[0] - 100.
        event.terminal = stop; event.direction = -1
        sol = solve_ivp(rhs, [start, end], np.asarray(initial).ravel(), method='DOP853',
            rtol=2e-12 if fine else 2e-10, max_step=60 if fine else 120,
            atol=np.tile([1e-5]*3 + [1e-9]*3, 361), dense_output=True, events=event)
        if not sol.success:
            raise RuntimeError(sol.message)
        finish = float(sol.t[-1])
        times = np.unique(np.r_[np.arange(start, finish, step), finish,
                               [t for t in [43200., self.arrival, self.arrival+21600.] if start <= t <= finish]])
        states = sol.sol(times).T.reshape(-1, 361, 6)
        frames = cmd(times).as_matrix()
        distances = np.array([closest_squares(y[:, :3], f)[0] for y, f in zip(states, frames)])
        events = sol.t_events[0].tolist()
        raw = dict(t=times, state=states, frames=frames, x=x, basis=self.basis,
                   controls=u, windows=windows, mass_ratio=self.ratio, distance_m=distances)
        row = dict(start_s=start, end_s=finish, intended_end_s=end, x=x.tolist(),
            event_stopped=stop, clearance_events_s=events, restoration_only=bool(events and not stop),
            minimum_sampled_distance_m=float(distances.min()),
            sampled_dates_below_floor=int(np.count_nonzero(distances < 100 - 1e-7)),
            nfev=sol.nfev, completed_horizon=abs(end-finish)<1e-6,
            suns=suns, area_grid=grid, max_step_s=60 if fine else 120, rtol=2e-12 if fine else 2e-10,
            accepted_return=False, accepted_cycle=False)
        path = spline(raw); mid = (times[:-1]+times[1:])/2
        row['reconstruction_error_m'] = float(length(path(mid)-sol.sol(mid).T.reshape(-1,361,6)[:,:,:3]).max())
        if account:
            nominal = rigid_square_load(cmd, times)['torque_per_mass']
            attitude = []; optical = []; torques = []
            for k, t in enumerate(times):
                light = load(self.env, t, states[k], frames[k], suns=suns, grid=grid, torque=True)
                torque = nominal[k] - light['gravity_torque_per_mass'] - light['radiation_torque_N_m']/self.mass[:,None]
                attitude.append(2*abs(torque).sum(axis=1)/10000)
                optical.append(light['optical']['first_intercept_bolometric_equivalent_W'])
                torques.append(length(torque*self.mass[:,None]))
            acceleration = self.thrust(x, times)
            attitude = np.array(attitude); optical = np.array(optical)
            power = (attitude + length(acceleration))*self.mass*SPECIFIC
            raw.update(power_W=power, attitude_acceleration=attitude, translation=acceleration,
                       first_intercept_W=optical, redirected_W=optical*self.env.central_fraction)
            row.update(electrical_energy_J=float(np.trapezoid(power.sum(axis=1),times)),
                first_intercept_J=float(np.trapezoid(optical.sum(axis=1),times)),
                maximum_individual_power_W=float(power.max()),
                minimum_installed_to_peak=float(np.min(self.capacity/power.max(axis=0))),
                overloaded_members=int(np.count_nonzero(power.max(axis=0)>self.capacity)),
                maximum_torque_N_m=float(np.max(torques)),
                maximum_translation_m_s2=float(length(acceleration).max()),
                maximum_attitude_rate_deg_s=float(np.rad2deg(length(cmd(times,1))).max()),
                actual_generation_W=None, actual_bus_delivery_W=None, actual_storage_state_J=None,
                actual_delivery_losses_J=None, actual_unmet_demand_J=None)
        row['cpu_s'] = time.process_time()-cpu; row['wall_s'] = time.monotonic()-wall
        return row, raw, sol


def limit():
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
