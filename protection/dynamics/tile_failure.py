"""Conditional local failure response for finite shingled shield tiles.

The failed tile's displacement from its prescribed healthy slot follows the
linear Hill--Clohessy--Wiltshire equations about a circular lunar orbit. Healthy
neighbours hold prescribed offsets in that rotating RTN frame. This is a local
disturbance experiment, not a propagation of the complete eccentric ring fleet.
The residual acceleration is a specified *difference* from healthy control,
constant in RTN; it is not a claim that sunlight supplies arbitrary thrust.

Parallel squares have an explicit box envelope, including user-specified
structural thickness. Surface clearance is exact for these parallel boxes.
Safety thresholds and the six faces of each collision box are located on the
continuous numerical solution. Face crossings also catch a thin envelope that
would be crossed between output samples. Only the prefix before first contact
has physical meaning; subsequent integration states are discarded.

Optional rescue uses an assumed independent, bidirectional acceleration source
after a specified detection/response delay. Its continuous capped feedback
cancels the differential orbit/force terms and damps position and velocity.
There are no impulses, position resets, automatic repairs or replacement tiles.
Flexible shape, changing mutual shadows, neighbour response, eccentricity,
third-body differential forces, debris and optical coverage are outside this
model. A separate orientation bound is supplied for lost attitude control.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math

import numpy as np
from scipy.integrate import solve_ivp

from shared.constants import MOON_GM, MOON_RADIUS


@dataclass(frozen=True)
class FailureConfig:
    radius_m: float = 20_000_000.0
    side_m: float = 10_000.0
    shingle_deg: float = 1.5
    alongtrack_pitch_m: float = 8_500.0
    ring_step_m: float | None = 1_000.0
    envelope_thickness_m: float = 1.0
    safety_clearance_m: float = 150.0
    duration_s: float = 21_600.0
    residual_acceleration_m_s2: tuple[float, float, float] = (1e-5, 0.0, 0.0)
    initial_state: tuple[float, ...] = (0., 0., 0., 0., 0., 0.)
    response_delay_s: float | None = None
    backup_acceleration_cap_m_s2: float = 1e-4
    recovery_time_scale_s: float = 600.0
    recovery_position_tolerance_m: float = 1.0
    recovery_speed_tolerance_m_s: float = 1e-3
    max_step_s: float = 30.0
    rtol: float = 1e-10
    atol: float = 1e-9

    def validate(self):
        values = [self.radius_m, self.side_m, self.alongtrack_pitch_m,
                  self.duration_s, self.backup_acceleration_cap_m_s2,
                  self.recovery_time_scale_s, self.max_step_s, self.rtol,
                  self.atol, self.recovery_position_tolerance_m,
                  self.recovery_speed_tolerance_m_s]
        if not all(math.isfinite(v) and v > 0 for v in values):
            raise ValueError("radii, sizes, durations, tolerances and caps must be positive finite values")
        if self.radius_m <= MOON_RADIUS:
            raise ValueError("reference orbit must lie above the lunar surface")
        if self.envelope_thickness_m < 0 or self.safety_clearance_m < 0:
            raise ValueError("envelope thickness and safety clearance must be nonnegative")
        if not all(math.isfinite(v) for v in [self.shingle_deg, self.envelope_thickness_m,
                                             self.safety_clearance_m]):
            raise ValueError("geometry values must be finite")
        if self.ring_step_m is not None and (not math.isfinite(self.ring_step_m) or self.ring_step_m <= 0):
            raise ValueError("ring step must be positive or None")
        if self.response_delay_s is not None and (not math.isfinite(self.response_delay_s) or self.response_delay_s < 0):
            raise ValueError("response delay must be nonnegative or None")
        if np.asarray(self.initial_state).shape != (6,) or not np.isfinite(self.initial_state).all():
            raise ValueError("initial_state must contain six finite RTN positions/velocities")
        if np.asarray(self.residual_acceleration_m_s2).shape != (3,) or not np.isfinite(self.residual_acceleration_m_s2).all():
            raise ValueError("residual acceleration must contain three finite RTN components")
        # The experiment is local. This does not independently bound all omitted
        # eccentric-orbit/tide errors; it just prevents a long-duration misuse.
        n = math.sqrt(MOON_GM / self.radius_m**3)
        if n * self.duration_s > math.pi / 2:
            raise ValueError("local failure experiment is limited to one quarter of the reference orbit")


def tile_frame(shingle_deg):
    """Columns: along-square, orbit-normal-square, square normal in RTN."""
    a = math.radians(shingle_deg)
    return np.array([[math.sin(a), 0., math.cos(a)],
                     [math.cos(a), 0., -math.sin(a)],
                     [0., 1., 0.]])


def default_neighbours(config):
    """Two same-ring neighbours and optional aligned inner/outer-ring tiles.

    The latter sample a local aligned encounter, not every possible sliding
    phase of the adjacent rings. Arbitrary extra offsets can be supplied.
    """
    neighbours = [[0., -config.alongtrack_pitch_m, 0.],
                  [0., config.alongtrack_pitch_m, 0.]]
    if config.ring_step_m is not None:
        neighbours += [[-config.ring_step_m, 0., 0.], [config.ring_step_m, 0., 0.]]
    return np.asarray(neighbours, float)


def hcw_acceleration(state, mean_motion):
    """Acceleration in rotating RTN coordinates for lunar circular HCW."""
    x, _, z, xd, yd, _ = np.asarray(state)[:6]
    n = mean_motion
    return np.array([3*n*n*x + 2*n*yd, -2*n*xd, -n*n*z])


def free_hcw_state(initial, times_s, mean_motion):
    """Analytical unforced HCW solution, used as an independent check."""
    x, y, z, xd, yd, zd = np.asarray(initial, float)
    t = np.asarray(times_s, float)
    n = mean_motion
    s, c = np.sin(n*t), np.cos(n*t)
    return np.stack([(4-3*c)*x+s/n*xd+2*(1-c)/n*yd,
                     6*(s-n*t)*x+y-2*(1-c)/n*xd+(4*s-3*n*t)/n*yd,
                     c*z+s/n*zd,
                     3*n*s*x+c*xd+2*s*yd,
                     6*n*(c-1)*x-2*s*xd+(4*c-3)*yd,
                     -n*s*z+c*zd], axis=-1)


def parallel_box_clearances(displacement, neighbours, frame, side_m, thickness_m):
    """Exact surface separation (zero for overlapping parallel box envelopes)."""
    projected = (np.asarray(neighbours)-np.asarray(displacement)) @ frame
    outside = np.maximum(np.abs(projected)-[side_m, side_m, thickness_m], 0.)
    return np.linalg.norm(outside, axis=-1)


def rotation_clearance_lower_bound(initial_clearance_m, side_m, angular_rate_rad_s, time_s,
                                   envelope_thickness_m=0.):
    """Conservative clearance if only one square loses attitude control.

    Every point lies within sqrt(side^2/2 + thickness^2/4) of the centre. An orientation differing
    by at most rate*time moves any point by at most 2*r*sin(angle/2). Subtracting
    that from the held-attitude separation bounds surface clearance. A negative
    bound means unproved separation, not a demonstrated physical collision.
    The rate must bound the actual angular speed; no torque model is assumed.
    """
    if side_m <= 0 or envelope_thickness_m < 0 or angular_rate_rad_s < 0 or np.any(np.asarray(time_s) < 0):
        raise ValueError("size must be positive; rate and times must be nonnegative")
    angle = np.minimum(np.asarray(time_s)*angular_rate_rad_s, math.pi)
    diameter = math.sqrt(2*side_m**2+envelope_thickness_m**2)
    return np.asarray(initial_clearance_m)-diameter*np.sin(angle/2)


def rotation_bound_window_s(initial_clearance_m, clearance_gate_m, side_m, angular_rate_rad_s,
                            envelope_thickness_m=0.):
    """First loss of the conservative rotation-clearance guarantee."""
    if side_m <= 0 or envelope_thickness_m < 0 or angular_rate_rad_s < 0:
        raise ValueError("size must be positive and rate nonnegative")
    margin = initial_clearance_m-clearance_gate_m
    if margin <= 0:
        return 0.
    diameter = math.sqrt(2*side_m**2+envelope_thickness_m**2)
    if angular_rate_rad_s == 0 or margin > diameter:
        return None
    return 2*math.asin(margin/diameter)/angular_rate_rad_s


def simulate_failure(config: FailureConfig, neighbour_offsets_m=None):
    """Return a JSON-safe local response account; no global fleet claims.

    Extrema are sampled on the solver mesh plus 4 subdivisions per accepted
    step; first gate/contact times are continuous-solution roots. Repeat at a
    smaller max_step_s to assess event and sampled-extrema convergence.
    """
    config.validate()
    neighbours = default_neighbours(config) if neighbour_offsets_m is None else np.asarray(neighbour_offsets_m, float)
    if neighbours.ndim != 2 or neighbours.shape[1] != 3 or len(neighbours) == 0 or not np.isfinite(neighbours).all():
        raise ValueError("neighbour offsets must be a nonempty finite (N,3) array")
    frame = tile_frame(config.shingle_deg)
    extent = np.array([config.side_m, config.side_m, config.envelope_thickness_m])
    n = math.sqrt(MOON_GM/config.radius_m**3)
    residual = np.asarray(config.residual_acceleration_m_s2)

    def control(state, active):
        if not active:
            return np.zeros(3)
        tau = config.recovery_time_scale_s
        demand = -residual-hcw_acceleration(state, n)-state[:3]/tau**2-2*state[3:6]/tau
        norm = np.linalg.norm(demand)
        return demand * min(1., config.backup_acceleration_cap_m_s2/max(norm, 1e-300))

    def projected(state):
        return (neighbours-state[:3]) @ frame

    def clearances(state):
        return parallel_box_clearances(state[:3], neighbours, frame, config.side_m, config.envelope_thickness_m)

    events, event_kinds = [], []
    for j in range(len(neighbours)):
        def safety(t, state, j=j):
            return clearances(state)[j]-config.safety_clearance_m
        safety.direction = -1
        events.append(safety)
        event_kinds.append(("safety", j, None, None))
        # Each continuous face coordinate changes sign when entering or leaving
        # that face, even when the whole thin box is traversed in one output step.
        for axis in range(3):
            for sign in (-1., 1.):
                def face(t, state, j=j, axis=axis, sign=sign):
                    return projected(state)[j, axis]-sign*extent[axis]
                events.append(face)
                event_kinds.append(("face", j, axis, sign))

    initial = np.r_[config.initial_state, 0.]
    initial_clearances = clearances(initial)
    safety_hits = [(0., j) for j, d in enumerate(initial_clearances) if d <= config.safety_clearance_m]
    contact_hits = [(0., j) for j, p in enumerate(projected(initial)) if np.all(np.abs(p) <= extent)]
    solutions = []
    boundary = [0., config.duration_s]
    delay = config.response_delay_s
    if delay is not None and 0. < delay < config.duration_s:
        boundary.insert(1, delay)
    state = initial
    for start, end in zip(boundary[:-1], boundary[1:]):
        if contact_hits:
            break
        active = delay is not None and start >= delay

        def rhs(t, y):
            applied = control(y, active)
            return np.r_[y[3:6], hcw_acceleration(y, n)+residual+applied, np.linalg.norm(applied)]

        sol = solve_ivp(rhs, (start, end), state, method="DOP853", dense_output=True,
                        events=events, max_step=config.max_step_s,
                        rtol=config.rtol, atol=config.atol)
        if not sol.success:
            raise RuntimeError(sol.message)
        solutions.append((sol, active))
        for kind, times, states in zip(event_kinds, sol.t_events, sol.y_events):
            for time, hit in zip(times, states):
                label, j, _, _ = kind
                if label == "safety":
                    safety_hits.append((float(time), j))
                elif np.all(np.abs(projected(hit)[j]) <= extent+1e-6):
                    contact_hits.append((float(time), j))
        state = sol.y[:, -1]

    contact = min(contact_hits) if contact_hits else None
    stop = contact[0] if contact else config.duration_s
    safety_hits = [v for v in safety_hits if v[0] <= stop+1e-8]
    safety = min(safety_hits) if safety_hits else None
    sample_times, sample_states, applied_accelerations = [0.], [initial], [np.zeros(3)]
    for sol, active in solutions:
        if sol.t[0] > stop:
            break
        knots = np.unique(np.r_[sol.t[sol.t <= stop], min(stop, sol.t[-1])])
        dense_times = np.unique(np.concatenate([np.linspace(a, b, 5) for a, b in zip(knots[:-1], knots[1:])])) if len(knots) > 1 else knots
        for t, y in zip(dense_times, sol.sol(dense_times).T):
            sample_times.append(float(t))
            sample_states.append(y)
            applied_accelerations.append(control(y, active))
    states = np.asarray(sample_states)
    final = states[-1]
    positions = np.linalg.norm(states[:, :3], axis=1)
    speeds = np.linalg.norm(states[:, 3:6], axis=1)
    all_clearances = np.array([clearances(s) for s in states])
    peak_acceleration = float(np.linalg.norm(applied_accelerations, axis=1).max())
    at_slot = positions[-1] <= config.recovery_position_tolerance_m and speeds[-1] <= config.recovery_speed_tolerance_m_s
    return {
        "schema": "terluna.protection.local-tile-failure/1",
        "evidence": "Local linear circular-orbit disturbance calculation with prescribed healthy neighbours; not fleet collision certification, optical coverage, flexible structure or demonstrated actuator hardware.",
        "config": asdict(config),
        "neighbour_offsets_m": neighbours.tolist(),
        "reference_mean_motion_rad_s": n,
        "initial_minimum_surface_clearance_m": float(initial_clearances.min()),
        "first_safety_crossing_s": None if safety is None else safety[0],
        "first_safety_neighbour_index": None if safety is None else safety[1],
        "first_envelope_contact_s": None if contact is None else contact[0],
        "first_contact_neighbour_index": None if contact is None else contact[1],
        "physical_prefix_end_s": stop,
        "status": "envelope_contact" if contact else "safety_gate_crossed" if safety else "clear_through_horizon",
        "minimum_sampled_surface_clearance_m": float(all_clearances.min()),
        "maximum_sampled_displacement_m": float(positions.max()),
        "maximum_sampled_relative_speed_m_s": float(speeds.max()),
        "maximum_displacement_over_radius": float(positions.max()/config.radius_m),
        "maximum_sampled_backup_acceleration_m_s2": peak_acceleration,
        "backup_integrated_acceleration_m_s": float(final[6]),
        "final_state_rtn_m_m_s": final[:6].tolist(),
        "final_position_error_m": float(positions[-1]),
        "final_speed_error_m_s": float(speeds[-1]),
        "returned_to_slot_by_horizon": bool(contact is None and delay is not None and delay < stop and at_slot),
        "filter_function_restored": False,
        "interpretation": "A recovered position does not repair a failed optical device. Replacement/removal, independent actuator mass and power, adjacent-ring sliding phases and global coverage remain separate. Safety roots are numerical events, not a proof against all possible grazing encounters; repeat with a refined step.",
    }
