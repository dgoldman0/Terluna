"""Bounded direct-collocation search for specular-sail lunar return loops."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from scipy.sparse import lil_matrix

from shared import constants as K
from protection.dynamics.ephemeris import Ephemeris
from protection.dynamics.cycling import (
    CyclingEnvironment, dro_seed, cr3bp_derivative, to_inertial, from_inertial)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RUN = ROOT/"research/runs/solar_shield_array/cycling"
SCENARIO = json.loads((HERE/"scenario.json").read_text())
CENTRAL = 1-SCENARIO["climate_flux_before_dimming_W_m2"]*SCENARIO["climate_sunlight_scale"]/K.SOLAR_CONSTANT


class ReturnProblem:
    """Hermite–Simpson states, piecewise-linear passive sail angles.

    The boundary returns to the same position/velocity in the *instantaneous*
    Earth–Moon frame. This is one finite ephemeris arc, not a periodic solution
    of the ephemeris. Control regularization selects among feasible solutions;
    its residual must be reported separately from dynamics/closure defects.
    """
    def __init__(self, env, radius, duration, start=0., intervals=64):
        self.env, self.radius, self.duration, self.start = env, radius, duration, start
        self.intervals = intervals
        self.times = np.linspace(start, start+duration, intervals+1)
        self.midtimes = (self.times[:-1]+self.times[1:])/2
        self.L = radius
        self.T = np.sqrt(radius**3/K.MOON_GM)
        self.scale = np.array([self.L]*3+[self.L/self.T]*3)
        self.h = duration/intervals/self.T
        self.frame = env.frames(self.times)
        self.frame_d = env.frames(self.times, 1)
        self.control_weight = 1e-7

    def derivative(self, t, state, angles):
        acceleration = self.env.acceleration(t, state[..., :3]*self.L, angles)
        return np.concatenate([state[..., 3:], acceleration*self.T**2/self.L], axis=-1)

    def unpack(self, vector):
        nodes = vector.reshape(self.intervals+1, 8)
        return nodes[:, :6], nodes[:, 6:]

    def seed(self, local_initial, phase=0.):
        mu = K.MOON_GM/(K.EARTH_GM+K.MOON_GM)
        D = K.EARTH_MOON_DISTANCE
        Tu = np.sqrt(D**3/(K.EARTH_GM+K.MOON_GM))
        seed = local_initial/np.array([D]*3+[D/Tu]*3)
        seed[0] += 1-mu
        sol = solve_ivp(cr3bp_derivative, [0, (1+phase)*self.duration/Tu], seed,
                        args=(mu,), t_eval=(self.times-self.start+phase*self.duration)/Tu,
                        rtol=2e-11, atol=2e-13)
        y = sol.y.T
        y[:, 0] -= 1-mu
        y *= np.array([D]*3+[D/Tu]*3)
        self.anchor = y[0, :3]/self.L
        inertial = to_inertial(y, self.frame, self.frame_d)/self.scale
        angles = np.tile([np.deg2rad(65.), 0.], (len(y), 1))
        return np.column_stack([inertial, angles]).ravel()

    def residual(self, vector, regularize=True):
        state, angles = self.unpack(vector)
        f = self.derivative(self.times, state, angles)
        mid = (state[:-1]+state[1:])/2+self.h/8*(f[:-1]-f[1:])
        fm = self.derivative(self.midtimes, mid, (angles[:-1]+angles[1:])/2)
        defect = state[1:]-state[:-1]-self.h/6*(f[:-1]+4*fm+f[1:])
        local = from_inertial(state[[0, -1]]*self.scale, self.frame[[0, -1]], self.frame_d[[0, -1]])/self.scale
        closure = local[1]-local[0]
        anchor = 1e-5*(local[0, :3]-getattr(self, "anchor", np.array([1., 0., 0.])))
        result = [defect.ravel(), closure, anchor]
        if regularize:
            # Small finite regularization; no claimed exact feasibility from
            # optimizer success alone. Independent propagation is mandatory.
            result += [self.control_weight*np.diff(angles, axis=0).ravel(),
                       self.control_weight*(angles[:, 0]-np.deg2rad(65.))]
        return np.concatenate(result)

    def sparsity(self):
        n = self.intervals
        rows = 6*n+9+2*n+(n+1)
        sp = lil_matrix((rows, 8*(n+1)), dtype=int)
        for i in range(n):
            sp[6*i:6*i+6, 8*i:8*i+16] = 1
        sp[6*n:6*n+6, :6] = 1
        sp[6*n:6*n+6, -8:-2] = 1
        sp[6*n+6:6*n+9, :3] = 1
        base = 6*n+9
        for i in range(n):
            for j in range(2):
                sp[base+2*i+j, 8*i+6+j] = 1
                sp[base+2*i+j, 8*(i+1)+6+j] = 1
        base += 2*n
        for i in range(n+1):
            sp[base+i, 8*i+6] = 1
        return sp.tocsr()

    def solve(self, seed, evaluations=150):
        lower = np.tile([-5]*3+[-5]*3+[0., -4*np.pi], (self.intervals+1, 1)).ravel()
        upper = np.tile([5]*3+[5]*3+[np.pi/2, 4*np.pi], (self.intervals+1, 1)).ravel()
        result = least_squares(self.residual, seed, bounds=(lower, upper),
                               jac_sparsity=self.sparsity(), x_scale="jac",
                               ftol=1e-11, xtol=1e-11, gtol=1e-11,
                               tr_solver="lsmr", tr_options={"atol": 1e-11, "btol": 1e-11},
                               max_nfev=evaluations)
        return result

    def propagate(self, vector, step=1800., rtol=2e-10, perturbation=None, enabled=True):
        state, angles = self.unpack(vector)
        initial = state[0]*self.scale
        if perturbation is not None:
            initial = initial+perturbation
        def control(t):
            return np.array([np.interp(t, self.times, angles[:, j]) for j in range(2)])
        def rhs(t, y):
            return np.r_[y[3:], self.env.acceleration(t, y[:3], control(t), enabled=enabled)]
        # Split at every control knot; the integrator does not step over a
        # piecewise-linear control corner unnoticed.
        states, times = [], []
        for a, b in zip(self.times[:-1], self.times[1:]):
            output = np.linspace(a, b, max(2, int(np.ceil((b-a)/step))+1))
            sol = solve_ivp(rhs, [a, b], initial, method="DOP853", t_eval=output,
                            rtol=rtol, atol=np.array([1e-3]*3+[1e-7]*3), max_step=step)
            if not sol.success:
                raise RuntimeError(sol.message)
            times.extend(output[:-1]); states.extend(sol.y.T[:-1])
            initial = sol.y[:, -1]
        times.append(self.times[-1]); states.append(initial)
        t, y = np.array(times), np.array(states)
        actual_angles = np.stack([control(ti) for ti in t])
        diagnostic = self.env.acceleration(t, y[:, :3], actual_angles, enabled=enabled, diagnostics=True)
        local = from_inertial(y[[0, -1]], self.env.frames(t[[0, -1]]), self.env.frames(t[[0, -1]], 1))
        return {"t": t, "state": y, "angles": actual_angles, "local_closure": local[1]-local[0], **diagnostic}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radius-km", type=float, default=50000.)
    parser.add_argument("--intervals", type=int, default=64)
    parser.add_argument("--evaluations", type=int, default=150)
    parser.add_argument("--start-days", type=float, default=0.)
    parser.add_argument("--phase", type=float, default=0.)
    parser.add_argument("--resume", type=Path, help="Interpolate a saved solution onto this mesh")
    args = parser.parse_args()
    RUN.mkdir(parents=True, exist_ok=True)
    before = time.monotonic()
    local, duration = dro_seed(args.radius_km*1000)
    start = args.start_days*K.JULIAN_DAY
    eph = Ephemeris(epoch=SCENARIO["epoch_tdb"])
    samples = eph.sample(np.arange(start-3600, start+duration+7200, 1800.))
    eph.close()
    env = CyclingEnvironment(samples, SCENARIO["baseline_areal_mass_kg_m2"], CENTRAL)
    problem = ReturnProblem(env, args.radius_km*1000, duration, start, args.intervals)
    seed = problem.seed(local, args.phase)
    if args.resume:
        previous = np.load(args.resume)
        if (float(previous["radius"]) != problem.radius or
                abs(float(previous["duration"])-duration) > 1e-5 or
                float(previous["start"]) != start):
            raise ValueError("Resume geometry does not match requested radius/time")
        nodes = previous["vector"].reshape(-1, 8)
        old_times = np.linspace(start, start+duration, len(nodes))
        seed = np.column_stack([np.interp(problem.times, old_times, nodes[:, j]) for j in range(8)]).ravel()
        problem.anchor = previous["anchor"]
    print(json.dumps({"stage": "seed", "radius_km": args.radius_km, "days": duration/K.JULIAN_DAY,
                      "defect_norm": float(np.linalg.norm(problem.residual(seed, False)))}), flush=True)
    solution = problem.solve(seed, args.evaluations)
    print(json.dumps({"stage": "optimized", "success": solution.success, "message": solution.message,
                      "nfev": solution.nfev, "cost": solution.cost,
                      "max_scaled_defect": float(np.max(np.abs(problem.residual(solution.x, False)))),
                      "elapsed_s": time.monotonic()-before}), flush=True)
    path = RUN/f"candidate_{args.radius_km:g}_{args.start_days:g}_{args.intervals}_{args.phase:g}.npz"
    np.savez_compressed(path, vector=solution.x, radius=args.radius_km*1000, duration=duration,
                        start=start, intervals=args.intervals, phase=args.phase, anchor=problem.anchor)
    trace = problem.propagate(solution.x)
    np.savez_compressed(path.with_name(path.stem+"_trace.npz"), **trace)
    print(json.dumps({"stage": "propagated", "closure_m": float(np.linalg.norm(trace["local_closure"][:3])),
                      "closure_m_s": float(np.linalg.norm(trace["local_closure"][3:])),
                      "range_km": [float(np.min(np.linalg.norm(trace["state"][:, :3], axis=1))/1000),
                                   float(np.max(np.linalg.norm(trace["state"][:, :3], axis=1))/1000)],
                      "elapsed_s": time.monotonic()-before}), flush=True)


if __name__ == "__main__":
    main()
