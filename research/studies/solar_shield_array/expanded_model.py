"""Expanded controls around the preserved loaded, all-member coupled model.

The force, eclipse and finite-square kernels are inherited unchanged. This
module exposes integration settings and compact ephemeris inputs, and adds
independent early in-plane and late acquisition controls. A post-contact
trajectory is a restoration diagnostic throughout its remaining duration.
"""
import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.spatial.transform import Rotation, RotationSpline

from protection.dynamics.cycling import CyclingEnvironment
from protection.dynamics.ephemeris import BODY_GM
from protection.dynamics.eclipse_parallel import EclipseParallelPattern, load
from protection.dynamics.attitude_load import rigid_square_load
from protection.dynamics.pattern_return import smooth_arc
from protection.dynamics.pattern_departure import sun_frame
from protection.dynamics.service_search import smooth_step
from protection.dynamics.square_distance import closest_squares
from protection.dynamics.optical import length
from .coupled_common import HERE, ROOT, SPECIFIC, inputs
from .coupled_restart import restore
from .coupled_stable import coalesce
from .coupled_cycle import transplant
from .joint_common import spline, state_at
from .cycling_search import CENTRAL

PACKAGE = HERE / 'local_continuation'
NAMES = ['row_parity_normal', 'column_parity_normal', 'parity_product_normal',
         'row_gradient_normal', 'pair_325_347_normal', 'pair_338_341_normal',
         'row_gradient_in_plane_u', 'column_gradient_in_plane_v',
         'pair_325_347_in_plane_v', 'pair_338_341_in_plane_u',
         'acquire_position_u', 'acquire_position_v', 'acquire_position_n',
         'acquire_velocity_u', 'acquire_velocity_v', 'acquire_velocity_n',
         'arrival_time', 'arrival_turn_duration']
DIM = len(NAMES)


def check_budget():
    """The runner installs an explicit process-clock guard, also usable in tests."""
    return None


def environment_from_compact(path=PACKAGE/'ephemeris.npz'):
    with np.load(path, allow_pickle=False) as z:
        samples = {k: z[k] for k in z.files if not k.startswith(('p_', 'v_'))}
        samples['positions'] = {n: z['p_'+n] for n in BODY_GM}
        samples['velocities'] = {n: z['v_'+n] for n in BODY_GM}
    return CyclingEnvironment(samples, .05, CENTRAL, 'filter_only')


class ExpandedExperiment:
    def __init__(self, package=PACKAGE):
        from .fleet_run import digest
        identity = json.loads((package/'inputs.json').read_text())
        for filename, key in [('ephemeris.npz', 'compact_ephemeris_sha256'), ('basis.npz', 'basis_sha256')]:
            if digest(package/filename) != identity[key]:
                raise ValueError('Compact input hash mismatch: '+filename)
        self.arrays, self.arrival0, self.producer = inputs()
        self.restart, old = restore()
        self.initial = old['initial_optimization_state']
        self.env = environment_from_compact(package/'ephemeris.npz')
        self.ratio = self.arrays['mass_ratio']
        self.mass = 5e6*self.ratio
        self.capacity = (self.ratio-1.2)*5e6*300
        self.early = np.zeros((361, 3, 10))
        self.early[:, :, :6] = old['basis']
        rows, cols = np.indices((19, 19))
        patterns = [(rows.ravel()-9)/9, (cols.ravel()-9)/9]
        for i, j in [(325, 347), (338, 341)]:
            p = np.zeros(361); p[i] = 1.; p[j] = -1.; patterns.append(p)
        frame = self.command(np.zeros(DIM))(13.8*3600).as_matrix()
        for col, (p, axis) in enumerate(zip(patterns, [0, 1, 1, 0]), 6):
            p = p-np.average(p, weights=self.mass); p /= max(abs(p))
            self.early[:, :, col] = .1*p[:, None]*frame[:, axis]
        with np.load(package/'basis.npz', allow_pickle=False) as z:
            self.late = z['late']
            self.reference_centre = z['reference_centre']
        self.basis = self.early

    def dates(self, x):
        arrival = self.arrival0+1800*x[16]
        return np.array([arrival, arrival+21600., arrival+46800.])

    def command(self, x):
        arrival = self.arrival0+1800*x[16]; end = arrival+46800.
        turn = 28800.+3600*x[17]
        times = coalesce(np.r_[np.arange(0., end, 120.), 0., 21600., 43200.,
                               arrival-turn, arrival, arrival+21600., arrival+43200., end])
        times = times[(times >= 0)&(times <= end)]
        first = smooth_step((times-21600.)/21600.)*(1-smooth_step((times-arrival+turn)/turn))
        second = smooth_step((times-arrival-21600.)/21600.)
        angle = -np.pi/2*(first+second)
        frames = np.array([sun_frame(self.env, t) for t in times]) @ Rotation.from_rotvec(
            np.c_[angle, np.zeros((len(times), 2))]).as_matrix()
        return RotationSpline(times, Rotation.from_matrix(frames))

    def controls(self, x):
        arrival = self.dates(x)[0]
        u = np.einsum('nic,c->ni', self.early, x[:10])
        p = np.einsum('nic,c->ni', self.late[:, :, :3], x[10:13])
        v = np.einsum('nic,c->ni', self.late[:, :, 3:], x[13:16])
        windows = np.array([[21600., 43200.], [arrival-28800., arrival-14400.],
                            [arrival-14400., arrival], [arrival+21600., arrival+43200.]])
        return np.stack([u, p, v, u], axis=1), windows

    def onset(self, j):
        if j < 10: return 21600.
        if j < 13: return self.arrival0-28800.
        if j < 16: return self.arrival0-14400.
        # Changing the attitude spline can affect preceding knots. Use a
        # conservative full common prefix through hour 30 for timing columns.
        return 30*3600.

    def evaluate(self, x, initial=None, start=21600., end=None, fine=False,
                 stop=False, max_step=480., output_step=300., account=True):
        check_budget()
        cpu = time.process_time(); x = np.asarray(x)
        cmd = self.command(x); end = self.dates(x)[-1] if end is None else end
        suns, grid = (16, 32) if fine else (8, 16)
        max_step = 60. if fine else max_step
        model = EclipseParallelPattern(self.env, cmd, self.ratio, suns=suns, grid=grid)
        u, windows = self.controls(x)
        def thrust(t): return np.einsum('k,nki->ni', smooth_arc(t, windows), u)
        def rhs(t, flat):
            check_budget()
            y = flat.reshape(361, 6)
            return np.c_[y[:, 3:], model.acceleration(t, y)+thrust(t)].ravel()
        def event(t, flat):
            check_budget()
            return closest_squares(flat.reshape(361, 6)[:, :3], cmd(t).as_matrix())[0]-100.
        event.terminal = stop; event.direction = -1
        y0 = self.initial if initial is None else initial
        sol = solve_ivp(rhs, [start, end], y0.ravel(), method='DOP853',
                        rtol=2e-12 if fine else 2e-10, max_step=max_step,
                        atol=np.tile([1e-5]*3+[1e-9]*3, 361), dense_output=True, events=event)
        if not sol.success: raise RuntimeError(sol.message)
        finish = float(sol.t[-1])
        special = np.r_[43200., 46800., self.dates(x), windows.ravel(),
                        self.arrival0-28800., self.arrival0-14400., 30*3600.]
        times = coalesce(np.r_[np.arange(start, finish, output_step), finish,
                              special[(special >= start)&(special <= finish)]], 1e-6)
        states = sol.sol(times).T.reshape(-1, 361, 6); frames = cmd(times).as_matrix()
        distances = np.array([closest_squares(y[:, :3], f)[0] for y, f in zip(states, frames)])
        raw = dict(t=times, state=states, frames=frames, x=x, controls=u,
                   windows=windows, distance_m=distances, mass_ratio=self.ratio)
        path = spline(raw); mid = (times[:-1]+times[1:])/2
        error = float(length(path(mid)-sol.sol(mid).T.reshape(-1, 361, 6)[:, :, :3]).max())
        row = dict(start_s=start, end_s=finish, intended_end_s=end, x=x.tolist(),
                   completed_horizon=bool(abs(finish-end)<1e-6), event_stopped=stop,
                   clearance_events_s=sol.t_events[0].tolist(),
                   restoration_only=bool(len(sol.t_events[0]) and not stop),
                   nfev=sol.nfev, suns=suns, area_grid=grid, max_step_s=max_step,
                   rtol=2e-12 if fine else 2e-10, reconstruction_error_m=error,
                   minimum_sampled_distance_m=float(distances.min()),
                   accepted_cycle=False)
        if account:
            nominal = rigid_square_load(cmd, times)['torque_per_mass']
            attitude = []; optical = []; torques = []
            for k, t in enumerate(times):
                check_budget()
                light = load(self.env, t, states[k], frames[k], suns=suns, grid=grid, torque=True)
                torque = nominal[k]-light['gravity_torque_per_mass']-light['radiation_torque_N_m']/self.mass[:, None]
                attitude.append(2*abs(torque).sum(axis=1)/10000.)
                optical.append(light['optical']['first_intercept_bolometric_equivalent_W'])
                torques.append(torque*self.mass[:, None])
            att = np.array(attitude); optical = np.array(optical); acc = np.array([thrust(t) for t in times])
            power = (att+length(acc))*self.mass*SPECIFIC
            raw.update(power_W=power, attitude_acceleration=att, translation=acc,
                       first_intercept_W=optical, redirected_W=optical*self.env.central_fraction,
                       torque_N_m=np.array(torques))
            row.update(electrical_energy_J=float(np.trapezoid(power.sum(axis=1), times)),
                       maximum_individual_power_W=float(power.max()),
                       minimum_installed_to_peak=float(np.min(self.capacity/power.max(axis=0))),
                       maximum_translation_m_s2=float(length(acc).max()),
                       maximum_attitude_rate_deg_s=float(np.rad2deg(length(cmd(times, 1))).max()),
                       overloaded_members=int(np.count_nonzero(power.max(axis=0)>self.capacity)))
        row['cpu_s'] = time.process_time()-cpu
        return row, raw, sol

    def terminal(self, raw, x):
        """Actual coupled state against transported repeated phase, no tail STM.

        The centre follows the saved baseline at each matched phase. It is
        retained in the residual: an arbitrary centre displacement is not free.
        """
        cmd = self.command(x); path = spline(raw); dates = self.dates(x)
        initial = self.arrays['initial_epoch_state']; service = self.initial
        reference_states = [initial, service, state_at(path, 46800.)]
        phases = [0., 21600., 46800.]
        targets = np.array([transplant(y, cmd(p).as_matrix(), cmd(p, 1), cmd(t).as_matrix(), cmd(t, 1), c)
                            for y, p, t, c in zip(reference_states, phases, dates, self.reference_centre)])
        residual = state_at(path, dates)-targets
        return residual

    def summary(self, raw, x):
        residual = self.terminal(raw, x)
        return dict(dates_s=self.dates(x).tolist(),
                    stages=['next_service_start', 'next_service_end', 'subsequent_departure'],
                    maximum_position_error_m=length(residual[:, :, :3]).max(axis=1).tolist(),
                    maximum_velocity_error_m_s=length(residual[:, :, 3:]).max(axis=1).tolist(),
                    normalized_rms=float(np.sqrt(np.mean((residual/np.array([1e5]*3+[5.]*3))**2))),
                    commanded_service_attitude_error_rad=0.,
                    recurrence_test_reached=False, accepted_cycle=False)
