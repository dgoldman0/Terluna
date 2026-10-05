"""Finite-pattern return proposals with constrained, reusable shooting maps.

The ephemeris reference has sail-assisted natural dynamics. A gravity
variational model gives a linear control proposal; it is not a coupled fleet
solution. Endpoint assignments and encounter cuts act on every actual member.
"""
import warnings

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import linprog, brentq
from scipy.sparse import eye, kron, csr_matrix, hstack, vstack
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation, RotationSpline

from shared import constants as K
from .active_formation import solar_frame
from .cycling import ray_geometry, visible_sun
from .fleet import safe_normals, square_vertices, beam_margins
from .natural_pattern import transported_roll
from .optical import unit, length


def common_normal(env, t, state, allowance=200000.):
    """A conservative common plane for the small reference neighbourhood."""
    return safe_normals(state[None, :], env.at(t), 4*K.MOON_RADIUS+allowance,
                        2*allowance)[0][0]


def reference_acceleration(env, t, state):
    sample = env.at(t); rays = ray_geometry(state[:3], sample)
    n = common_normal(env, t, state)
    cosine = max(float(n@(-unit(rays['sun']))), 0.)
    visible = visible_sun(rays['sun'], rays['earth'], rays['moon'])
    scale = 2*env.central_fraction*K.SOLAR_CONSTANT/(K.SPEED_OF_LIGHT*env.sigma)
    sail = scale*(K.AU/length(rays['sun']))**2*visible*cosine**2*(9890/10000)**2*n
    return env.gravity(t, state[:3])+sail


def natural_reference(env, initial, start, end):
    def rhs(t, y):
        return np.r_[y[3:], reference_acceleration(env, t, y)]
    solution = solve_ivp(rhs, [start, end], initial, method='DOP853', max_step=300.,
        rtol=2e-12, atol=[1e-5]*3+[1e-9]*3, dense_output=True)
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution


def smooth_arc(t, windows):
    """Sin-squared acceleration shapes, each with unit time integral."""
    values = np.zeros(len(windows))
    for j, (a, b) in enumerate(windows):
        if a <= t <= b:
            values[j] = 2/(b-a)*np.sin(np.pi*(t-a)/(b-a))**2
    return values


def shooting_maps(env, reference, start, end, burn=7200.):
    middle = (start+end)/2
    windows = np.array([[start+900., start+900.+burn],
                        [middle-burn/2, middle+burn/2], [end-burn, end]])
    if windows[0, 1] >= windows[1, 0] or windows[1, 1] >= windows[2, 0]:
        raise ValueError('Control arcs overlap')

    def rhs(t, y):
        phi, control = y[:36].reshape(6, 6), y[36:].reshape(6, 9)
        centre = reference.sol(t)[:3]
        eps = 10.
        force = env.gravity(t, centre+np.r_[np.eye(3), -np.eye(3)]*eps)
        gradient = ((force[:3]-force[3:])/(2*eps)).T
        jacobian = np.block([[np.zeros((3, 3)), np.eye(3)], [gradient, np.zeros((3, 3))]])
        drive = np.zeros((6, 9))
        drive[3:] = np.kron(smooth_arc(t, windows)[None, :], np.eye(3))
        return np.r_[(jacobian@phi).ravel(), (jacobian@control+drive).ravel()]

    solution = solve_ivp(rhs, [start, end], np.r_[np.eye(6).ravel(), np.zeros(54)],
        method='DOP853', max_step=300., rtol=2e-11, atol=1e-9, dense_output=True)
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution, windows


def proposal_state(reference, maps, initial, controls, times):
    times = np.atleast_1d(times)
    matrix = maps.sol(times).T
    phi = matrix[:, :36].reshape(-1, 6, 6)
    response = matrix[:, 36:].reshape(-1, 6, 9)
    relative = initial-reference.sol(maps.t[0])
    return (reference.sol(times).T[:, None, :]+
        np.einsum('tij,nj->tni', phi, relative)+np.einsum('tij,nj->tni', response, controls))


def attitude_command(env, reference, start, end, roll=0.):
    """Transport roll smoothly and recover the service square's orientation."""
    times = np.linspace(start, end, int(np.ceil((end-start)/120))+1)
    orientations = []
    for t in times:
        state = reference.sol(t)
        n = common_normal(env, t, state)
        sun = ray_geometry(state[:3], env.at(t))['sun']
        corners = square_vertices(state[None, :3], n[None, :], sun[None, :], 10000.)[0]
        orientations.append(np.column_stack([(corners[1]-corners[0])/10000,
                            (corners[3]-corners[0])/10000, n]))
    # Match the saved Sun-facing service plane exactly at both boundaries.
    for index, t in [(0, start), (-1, end)]:
        u, b, c = solar_frame(env.at(t)).T
        orientations[index] = np.column_stack([b, -c, -u])
    transported = transported_roll(times, np.array(orientations))
    frames = transported(times).as_matrix()
    desired = solar_frame(env.at(end))[:, 1]
    twist = np.arctan2(frames[-1, :, 1]@desired, frames[-1, :, 0]@desired)
    f = (times-start)/(end-start)
    angle = twist*(3*f*f-2*f*f*f)+roll*np.sin(np.pi*f)**2
    frames = frames@Rotation.from_rotvec(np.c_[np.zeros((len(times), 2)), angle]).as_matrix()
    return RotationSpline(times, Rotation.from_matrix(frames))


def minimum_component_impulse(maps, reference, initial, target, windows, cuts=(), cap=.001):
    """LP minimizes componentwise L1 impulse with exact linear endpoints.

    Report Euclidean impulse separately. A component box conservatively
    enforces the vector acceleration cap for each non-overlapping smooth arc.
    Encounter cuts select one separating axis and are therefore local.
    """
    count = len(initial); size = count*9
    terminal = maps.y[:, -1]
    phi, response = terminal[:36].reshape(6, 6), terminal[36:].reshape(6, 9)
    rhs = target-reference.sol(maps.t[-1])-(initial-reference.sol(maps.t[0]))@phi.T
    scale = np.array([1e4]*3+[1.]*3)
    equal = hstack([kron(eye(count), csr_matrix(response/scale[:, None])), csr_matrix((count*6, size))], format='csr')
    identity = eye(size, format='csr')
    upper = [hstack([identity, -identity]), hstack([-identity, -identity])]
    bounds_rhs = [np.zeros(size*2)]
    if cuts:
        rows, columns, values, limits = [], [], [], []
        initial_relative = initial-reference.sol(maps.t[0])
        for row, cut in enumerate(cuts):
            t, i, j, axis, required = cut
            data = maps.sol(t)
            p, b = data[:36].reshape(6, 6), data[36:].reshape(6, 9)
            coefficient = axis@b[:3]
            baseline = axis@(p[:3]@(initial_relative[i]-initial_relative[j]))
            rows.extend([row]*18)
            columns.extend([*(i*9+np.arange(9)), *(j*9+np.arange(9))])
            values.extend(np.r_[-coefficient, coefficient]/1e4)
            limits.append((baseline-required)/1e4)
        constraint = csr_matrix((values, (rows, columns)), shape=(len(cuts), 2*size))
        upper.append(constraint); bounds_rhs.append(np.array(limits))
    limit = np.repeat(cap*np.diff(windows, axis=1).ravel()/(2*np.sqrt(3)), 3)
    bounds = [(-float(x), float(x)) for _ in range(count) for x in limit]+[(0., None)]*size
    # HiGHS options are forwarded by SciPy; its one-thread setting is explicit.
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='Unrecognized options detected')
        result = linprog(np.r_[np.zeros(size), np.ones(size)], A_ub=vstack(upper, format='csr'),
            b_ub=np.r_[*bounds_rhs], A_eq=equal, b_eq=(rhs/scale).ravel(), bounds=bounds,
            method='highs', options={'threads': 1, 'parallel': False, 'time_limit': 40.})
    output = dict(success=bool(result.success), status=int(result.status), message=result.message,
                  cuts=len(cuts), acceleration_cap_m_s2=cap)
    if result.success:
        control = result.x[:size].reshape(count, 9)
        per_tile = length(control.reshape(count, 3, 3)).sum(axis=1)
        defect = control@response.T-rhs
        output.update(component_L1_mean_m_s=float(result.fun/count),
            mean_delta_v_m_s=float(per_tile.mean()), maximum_delta_v_m_s=float(per_tile.max()),
            endpoint_position_residual_m=float(length(defect[:, :3]).max()),
            endpoint_velocity_residual_m_s=float(length(defect[:, 3:]).max()))
        return control, output
    return None, output


def encounter_scan(times, states, command, side=10000., clearance=100., retain=256):
    """Exact static parallel-square SAT, retaining worst sampled pair witnesses."""
    witnesses = []; minimum = np.inf; count = 0
    for t, state, frame in zip(times, states, command(times).as_matrix()):
        q = state[:, :3]
        pairs = cKDTree(q).query_pairs(np.sqrt(2)*side+clearance, output_type='ndarray')
        if not len(pairs):
            continue
        delta = (q[pairs[:, 0]]-q[pairs[:, 1]])@frame
        values = np.abs(delta)-np.array([side, side, 0.])
        axis = values.argmax(axis=1)
        separation = values[np.arange(len(pairs)), axis]
        minimum = min(minimum, float(separation.min()))
        bad = np.flatnonzero(separation < clearance)
        count += len(bad)
        for k in bad:
            direction = frame[:, axis[k]]*(1 if delta[k, axis[k]] >= 0 else -1)
            extent = side if axis[k] < 2 else 0.
            witnesses.append(dict(time_s=float(t), i=int(pairs[k, 0]), j=int(pairs[k, 1]),
                separating_projection_m=float(separation[k]), axis=direction.tolist(),
                required_projection_m=float(extent+clearance)))
    witnesses.sort(key=lambda w: w['separating_projection_m'])
    return dict(sampled_pair_violations=count, minimum_sampled_separating_projection_m=minimum,
                witnesses=witnesses[:retain])


def beam_and_rate_screen(env, times, states, command):
    margins = []
    for t, state, frame in zip(times, states, command(times).as_matrix()):
        rays = ray_geometry(state[:, :3], env.at(t))
        margins.append(float(beam_margins(state[:, :3], rays['sun'], frame[:, 2],
                                        env.at(t), 4*K.MOON_RADIUS, 10000.).min()))
    return dict(minimum_sampled_beam_margin_deg=float(np.rad2deg(min(margins))),
        maximum_sampled_rate_deg_s=float(np.rad2deg(length(command(times, 1))).max()),
        maximum_sampled_angular_acceleration_deg_s2=float(np.rad2deg(length(command(times, 2))).max()))


def crossing_witnesses(times, states, command, state_at, side=10000., retain=256):
    """Seek exact coplanarity between snapshots, then test finite square edges.

    A shared orientation makes centre coplanarity a scalar root. Candidate
    selection is deliberately generous. This detects missed impacts; absence
    of a detected crossing is not a continuous separation certificate.
    """
    found = []
    frames = command(times).as_matrix()
    for k in range(len(times)-1):
        q0, q1 = states[k, :, :3], states[k+1, :, :3]
        relative_travel = float(2*length((q1-q0)-(q1-q0).mean(axis=0)).max())
        pairs = cKDTree((q0+q1)/2).query_pairs(np.sqrt(2)*side+relative_travel+500., output_type='ndarray')
        if not len(pairs):
            continue
        i, j = pairs.T
        d0 = (q0[i]-q0[j])@frames[k, :, 2]
        d1 = (q1[i]-q1[j])@frames[k+1, :, 2]
        for pair in pairs[d0*d1 < 0]:
            def delta(t):
                y = state_at(t, pair)
                return (y[0, :3]-y[1, :3])@command(t).as_matrix()
            root = brentq(lambda t: delta(t)[2], times[k], times[k+1], xtol=.001)
            projected = delta(root)
            if max(abs(projected[:2])) >= side+100:
                continue
            # At a coplanar crossing, an in-plane escape moves edges apart.
            # Moving the crossing time alone cannot eliminate the encounter.
            axis = int(np.argmax(abs(projected[:2])))
            direction = command(root).as_matrix()[:, axis]*(1 if projected[axis] >= 0 else -1)
            found.append(dict(time_s=float(root), i=int(pair[0]), j=int(pair[1]),
                coplanarity_residual_m=float(abs(projected[2])),
                in_plane_projections_m=projected[:2].tolist(),
                physical_square_intersection=bool(max(abs(projected[:2])) < side),
                separating_projection_m=float(max(abs(projected[:2]))-side),
                axis=direction.tolist(), required_projection_m=side+300.))
    found.sort(key=lambda w: w['separating_projection_m'])
    return dict(crossing_near_encounters=len(found),
        physical_intersections=sum(w['physical_square_intersection'] for w in found),
        witnesses=found[:retain])
