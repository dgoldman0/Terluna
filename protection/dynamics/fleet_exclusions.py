"""Vectorized backend for the reference fleet's conservative swept guards.

Two dense N-by-N arrays replace repeated Python updates of the same pair.
This backend is intended for the few-thousand-member study populations.
The reference implementation in fleet.py defines the guard semantics.
"""
import numpy as np
from scipy.spatial import cKDTree

from .optical import length, unit
from .cycling import sail_basis


def swept_conflicts(times, states, side, acceleration_bound=.15, clearance=100.,
                    solar_angular_radius=.0048, sun_directions=None):
    count = states.shape[1]
    minimum = [np.full((count, count), np.inf), np.full((count, count), np.inf)]
    threshold = np.sqrt(2)*side+clearance
    max_radius = float(length(states[..., :3]).max())
    for k, dt in enumerate(np.diff(times)):
        q0, q1 = states[k, :, :3], states[k+1, :, :3]
        mid = (q0+q1)/2; move = q1-q0
        sag = acceleration_bound*dt*dt/4
        travel = float(length(move).max())
        def collect(coords, delta, limit, matrix, axial=None):
            pairs = cKDTree(coords).query_pairs(limit+travel+sag, output_type='ndarray')
            if not len(pairs):
                return
            i, j = pairs.T
            start = coords[i]-coords[j]-(delta[i]-delta[j])/2
            change = delta[i]-delta[j]
            f = np.clip(-np.sum(start*change, axis=-1)/np.maximum(np.sum(change*change, axis=-1), 1e-100), 0, 1)
            lower = length(start+f[:, None]*change)-sag
            bound = limit
            if axial is not None:
                depth = np.abs(axial[i]-axial[j])+travel
                bound = threshold+solar_angular_radius*depth+3e-7*dt*max_radius
            selected = lower < bound
            i, j, lower = i[selected], j[selected], lower[selected]
            # query_pairs returns each unordered pair exactly once per call.
            matrix[i, j] = np.minimum(matrix[i, j], lower)
        collect(mid, move, threshold, minimum[0])
        if sun_directions is not None:
            u = unit(sun_directions[k]+sun_directions[k+1]); _, b, c = sail_basis(u)
            projected = np.column_stack([mid@b, mid@c])
            delta = np.column_stack([move@b, move@c])
            collect(projected, delta, threshold+2*solar_angular_radius*max_radius,
                    minimum[1], mid@u)
    result = []
    for matrix in minimum:
        i, j = np.nonzero(np.isfinite(matrix))
        result.append({(int(a), int(b)): float(v) for a, b, v in zip(i, j, matrix[i, j])})
    return tuple(result)
