"""Validation of the garment solver against analytic cases, on the Moon and in an Earth control.

    python -m humanities.dance.garment.validate OUTDIR [--only terminal,cantilever,spin,convergence]

1. Terminal speed. A flat sheet released face-on reaches sqrt(2 m_a g / (rho C_N)); a sheet held at incidence beta
   glides to the pair v_n = sqrt(2 m_a g cos(beta) / (rho C_N)), v_t = sqrt(2 m_a g sin(beta) / (rho C_T)); two
   sheets lying in contact fall at sqrt(2) times the single sheet's speed when shielding is on.
2. Bending length. A strip clamped horizontal with an overhang l = 2 c settles to the heavy-elastica shape; Peirce's
   formula c = l (cos(theta/2) / (8 tan theta))^(1/3) applied to the simulated chord angle theta recovers
   c = (B / (m_a g))^(1/3), which is (g_E / g_M)^(1/3) = 1.82 times longer on the Moon.
3. Steady turn. The circle skirt on the body (arms raised clear) turned at 5, 10, 15 and 20 rpm; hem radius and flare
   angle against a rigid-cone conical pendulum and a hanging-chain (string) model, with and without drag.
4. Convergence: the cantilever in mesh size and timestep, the terminal speed in timestep, the turn in mesh size.
Results go to OUTDIR/validation.json with diagnostic PNGs.
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
import torch

from .body import Body, MotionPose, SpinPose, run
from .cloth import EARTH, MOON, Cloth, Settings, World
from .garments import circle_skirt
from .materials import FABRICS, bending_length
from .patterns import rectangle

TEST_PATTERN = '/media/projectspace/terluna-research/humanities/dance-runs/test_pattern.npz'
WORLDS = (MOON, EARTH)


# ------------------------------------------------------------------------------------------------- 1. terminal speed
def terminal(mat, world, beta_deg=0.0, layers=1, h=0.02, size=0.2, dt=1 / 2400, T=None, shielding=True):
    """Fall of a size x size sheet (or a stack of `layers` sheets 0.8 contact distances apart) held flat at incidence
    beta. Returns the time series of the mean velocity's normal and in-plane components."""
    pc = rectangle(size, size, h)
    V1 = pc.n_vertices
    b = math.radians(beta_deg)
    R = np.array([[math.cos(b), 0, -math.sin(b)], [0, 1, 0], [math.sin(b), 0, math.cos(b)]])
    uv = np.c_[pc.vertex_uv - size / 2, np.zeros(V1)]
    s = Settings(dt=dt, self_collision=layers > 1, shielding=shielding)
    if layers > 1:
        from .garments import merge
        pc = merge([(f'l{k}', rectangle(size, size, h)) for k in range(layers)])
        d = 0.8 * 0.75 * pc.edge_length()
        X0 = np.concatenate([uv + [0, 0, k * d] for k in range(layers)]) @ R.T + [0, 0, 50.0]
    else:
        X0 = uv @ R.T + [0, 0, 50.0]
    cl = Cloth(pc, X0, mat, world, s)
    n = R @ np.array([0, 0, 1.0])
    vt_scale = math.sqrt(2 * mat.areal_density * world.gravity / (world.air_density * s.c_normal))
    if T is None:
        T = 12 * vt_scale / world.gravity + (8.0 if beta_deg else 0.0)
    nsteps = int(T / dt)
    rec = max(1, nsteps // 400)
    ts, vn, vtan = [], [], []
    for k in range(nsteps):
        cl.substep(dt)
        if k % rec == rec - 1:
            v = cl.v.mean(0).cpu().numpy()
            if beta_deg:                                   # measure against the sheet's current mean normal
                _, _, _, nt, At = cl.tri_geometry(cl.x[:cl.V // layers] if layers > 1 else cl.x)
                n = (nt * At[:, None]).sum(0).cpu().numpy()
                n = n / np.linalg.norm(n)
            ts.append((k + 1) * dt)
            vn.append(float(v @ n))
            vtan.append(float(np.linalg.norm(v - (v @ n) * n)))
    return np.array(ts), np.array(vn), np.array(vtan), cl


def terminal_cases(out, fast=False):
    res = []
    mats = ['silk_chiffon', 'silk_satin_heavy'] if fast else list(FABRICS)
    for world in WORLDS:
        for name in mats:
            m = FABRICS[name]
            ts, vn, vt, _ = terminal(m, world)
            an = math.sqrt(2 * m.areal_density * world.gravity / (world.air_density * 1.17))
            res.append(dict(case='face-on', world=world.name, fabric=name, analytic=an, simulated=float(-vn[-20:].mean()),
                            rel_error=float(-vn[-20:].mean() / an - 1), t=ts.tolist()[::4], v=(-vn).tolist()[::4]))
            print(json.dumps({k: v for k, v in res[-1].items() if k not in ('t', 'v')}), flush=True)
    m = FABRICS['silk_chiffon']
    for world in WORLDS:
        beta = 30.0
        rigid = type(m)(**{**m.__dict__, 'name': 'chiffon_mass_rigid', 'bending_rigidity': 1e-3})
        ts, vn, vt, _ = terminal(rigid, world, beta_deg=beta)
        an_n = math.sqrt(2 * m.areal_density * world.gravity * math.cos(math.radians(beta)) / (world.air_density * 1.17))
        an_t = math.sqrt(2 * m.areal_density * world.gravity * math.sin(math.radians(beta)) / (world.air_density * 0.01))
        res.append(dict(case='glide 30 deg, rigid sheet of chiffon areal density', world=world.name,
                        fabric='silk_chiffon (B = 1e-3 N m)', analytic_normal=an_n,
                        simulated_normal=float(abs(vn[-20:].mean())), analytic_tangential=an_t,
                        simulated_tangential=float(vt[-20:].mean()),
                        rel_error_normal=float(abs(vn[-20:].mean()) / an_n - 1),
                        rel_error_tangential=float(vt[-20:].mean() / an_t - 1)))
        print(json.dumps(res[-1]), flush=True)
        ts, vn, vt, cl = terminal(m, world, layers=2)
        an = math.sqrt(2 * 2 * m.areal_density * world.gravity / (world.air_density * 1.17))
        res.append(dict(case='two layers in contact, shielding on', world=world.name, fabric='silk_chiffon',
                        analytic_single_face=an, simulated=float(-vn[-20:].mean()),
                        rel_error=float(-vn[-20:].mean() / an - 1)))
        print(json.dumps(res[-1]), flush=True)
    return res


# ------------------------------------------------------------------------------------------------ 2. bending length
def elastica(l, B, w, n=400):
    """Heavy cantilever, clamped horizontal at s = 0, free at s = l: B theta'' = -w (l - s) cos(theta), theta the
    downward slope. Returns s, x, y (y downward) and the chord angle at the tip."""
    from scipy.integrate import solve_ivp
    from scipy.optimize import brentq

    def shoot(k0):
        return solve_ivp(lambda s, y: [y[1], -(w / B) * (l - s) * math.cos(y[0])], (0, l), [0.0, k0],
                         rtol=1e-10, atol=1e-12, dense_output=True)
    # theta' at s=0 is positive (curving downward); the moment vanishes at the tip: theta'(l) = 0
    g = lambda k0: shoot(k0).y[1, -1]
    k0 = brentq(g, 1e-9, 50.0 / l)
    sol = shoot(k0)
    s = np.linspace(0, l, n)
    th = sol.sol(s)[0]
    x = np.concatenate([[0], np.cumsum(0.5 * (np.cos(th[1:]) + np.cos(th[:-1])) * np.diff(s))])
    y = np.concatenate([[0], np.cumsum(0.5 * (np.sin(th[1:]) + np.sin(th[:-1])) * np.diff(s))])
    return s, x, y, math.atan2(y[-1], x[-1])


def peirce_c(l, theta):
    return l * (math.cos(theta / 2) / (8 * math.tan(theta))) ** (1 / 3)


def cantilever(mat, world, n_el=12, width=0.025, dt_factor=1.0, settle_periods=4.0, start_scale=0.5):
    """Simulated overhang l = 2 c_pred with n_el elements along it; starts from the elastica of a strip with
    start_scale x the bending rigidity (a deliberately wrong guess) and settles with numerical damping."""
    g = world.gravity
    w = mat.areal_density * g
    c = bending_length(mat, g)
    l = 2 * c
    h = l / n_el
    pc = rectangle(width, l + h, h)
    uv = pc.vertex_uv
    s0, x0, y0, _ = elastica(l, start_scale * mat.bending_rigidity, w)
    s = uv[:, 1] - h
    xs = np.where(s > 0, np.interp(s, s0, x0), s)
    zs = np.where(s > 0, -np.interp(s, s0, y0), 0.0)
    X0 = np.c_[xs, uv[:, 0] - width / 2, 1.0 + zs]
    pinned = np.where(uv[:, 1] <= h + 1e-9)[0]
    dt = dt_factor * 0.0141 * h * h / math.sqrt(c ** 3 * g)
    st = Settings(dt=dt, self_collision=False, shielding=False)
    cl = Cloth(pc, X0, mat, world, st, pinned=pinned)
    cl.set_targets(torch.as_tensor(X0[pinned], dtype=torch.float32, device=cl.dev))
    Ts = 4 * math.sqrt(c / g)
    cl.damp.fill_(4.0 / Ts)
    nsteps = int(settle_periods * Ts / dt)
    tip_idx = np.array(pc.rings['end'])
    edge_idx = np.where(np.abs(uv[:, 1] - h) < 1e-9)[0]
    hist = []
    for k in range(nsteps):
        cl.substep(dt)
        if k % max(1, nsteps // 30) == 0 or k == nsteps - 1:
            x = cl.x.cpu().numpy()
            tip, edge = x[tip_idx].mean(0), x[edge_idx].mean(0)
            hist.append(math.atan2(edge[2] - tip[2], tip[0] - edge[0]))
    x = cl.x.cpu().numpy()
    tip, edge = x[tip_idx].mean(0), x[edge_idx].mean(0)
    theta = math.atan2(edge[2] - tip[2], tip[0] - edge[0])
    _, xe, ye, th_exact = elastica(l, mat.bending_rigidity, w)
    mid = np.abs(uv[:, 0] - width / 2) < 0.6 * h
    prof = x[mid]
    prof = prof[np.argsort(prof[:, 0])]
    return dict(world=world.name, fabric=mat.name, c_pred=c, overhang=l, h=h, n_el=n_el, dt=dt, substeps=nsteps,
                theta_sim_deg=math.degrees(theta), theta_elastica_deg=math.degrees(th_exact),
                c_sim=peirce_c(l, theta), c_peirce_from_elastica=peirce_c(l, th_exact),
                rel_error_c=peirce_c(l, theta) / c - 1, settle_drift_deg=math.degrees(hist[-1] - hist[-2]),
                profile_x=(prof[:, 0] - edge[0]).tolist(), profile_z=(prof[:, 2] - edge[2]).tolist(),
                elastica_x=xe[::10].tolist(), elastica_z=(-ye[::10]).tolist())


def cantilever_cases(fast=False):
    res = []
    mats = ['silk_chiffon'] if fast else list(FABRICS)
    for world in WORLDS:
        for name in mats:
            r = cantilever(FABRICS[name], world)
            res.append(r)
            print(json.dumps({k: v for k, v in r.items() if not k.startswith(('profile', 'elastica'))}), flush=True)
    return res


# ------------------------------------------------------------------------------------------------- 3. steady turn
def cone_estimate(omega, g, r_w, L):
    """Rigid straight meridians at flare phi, mass per length proportional to (r_w + s): torque balance about the
    waist. Returns phi (rad)."""
    from scipy.optimize import brentq
    s = np.linspace(0, L, 2001)
    lam = r_w + s
    trap = lambda f: float(np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(s)))

    def f(phi):
        return omega ** 2 * math.cos(phi) * trap((r_w + s * math.sin(phi)) * s * lam) - g * math.sin(phi) * trap(s * lam)
    if omega == 0:
        return 0.0
    return brentq(f, 1e-6, math.pi / 2 - 1e-9) if f(math.pi / 2 - 1e-9) < 0 else math.pi / 2


def chain_estimate(omega, g, rho, r_w, L, m_a, c_n=1.17, c_t=0.01, drag=True, folds=True, n=400, iters=4000,
                   relax=0.2):
    """A meridional string of the skirt (per radian of pattern angle) turning steadily with the waist in still air.
    Loads per unit length: weight, centrifugal force and, with drag, the circumferential air force. With folds the
    excess circumference is taken as a zig-zag whose faces slope at beta (cos(beta) = R / (r_w + s)) and both faces
    meet the flow: 1/2 rho (omega R)^2 (r_w + s) (C_N sin^3 beta + C_T cos^3 beta), an upper bound on fold drag that
    feeds back on itself (more lag, smaller R, deeper folds). Without folds the cloth is a smooth cone with skin
    friction only, a lower bound. The string lies along the sum of the loads below each point."""
    s = (np.arange(n) + 0.5) * L / n
    ds = L / n
    lam = m_a * (r_w + s)
    X = np.c_[r_w + 0.3 * s, np.zeros(n), -s * 0.95]
    res = np.inf
    for it in range(iters):
        Rh = np.hypot(X[:, 0], X[:, 1])
        f = np.zeros((n, 3))
        f[:, 2] = -lam * g
        f[:, 0] = lam * omega ** 2 * X[:, 0]
        f[:, 1] = lam * omega ** 2 * X[:, 1]
        if drag:
            cosb = np.clip(Rh / (r_w + s), 0, 1)
            sinb = np.sqrt(1 - cosb ** 2)
            ceff = (c_n * sinb ** 3 + c_t * cosb ** 3) if folds else np.full(n, c_t)
            vel = omega * np.c_[-X[:, 1], X[:, 0], np.zeros(n)]
            sp = np.linalg.norm(vel, axis=1)
            f += -0.5 * rho * ceff[:, None] * (r_w + s)[:, None] * sp[:, None] * vel
        Fb = np.cumsum((f * ds)[::-1], axis=0)[::-1]
        t = Fb / np.linalg.norm(Fb, axis=1, keepdims=True)
        Xn = np.array([r_w, 0.0, 0.0]) + np.cumsum(t * ds, axis=0) - 0.5 * t * ds
        res = np.abs(Xn - X).max()
        X = (1 - relax) * X + relax * Xn
        if res < 1e-9:
            break
    tip = X[-1] + 0.5 * ds * t[-1]
    R = math.hypot(tip[0], tip[1])
    return dict(hem_radius=R, drop=float(-tip[2]), flare_deg=math.degrees(math.atan2(R - r_w, -tip[2])),
                lag_deg=-math.degrees(math.atan2(tip[1], tip[0])), converged=bool(res < 1e-6), residual=float(res))


def spin_runs(rpms=(5, 10, 15, 20), h=0.02, worlds=WORLDS, fabric='silk_chiffon', settle=5.0, t_up=2.0, hold=5.0,
              n_sub=40, fps=60.0):
    out = []
    body = Body()
    mp = MotionPose(TEST_PATTERN, body, t_raised=1.0, t_lower=1.0, t_hold=1.0)
    p0, q0 = mp(-mp.pre)                              # first frame, arms raised to the side
    mat = FABRICS[fabric]
    for world in worlds:
        G = circle_skirt(body, (p0, q0), h=h)
        st = Settings(dt=1.0 / (fps * n_sub))
        cl = Cloth(G.piece, G.X0, mat, world, st, pinned=G.pinned)
        still = SpinPose(p0, q0, 0.0)
        t0 = time.time()
        run(cl, still, G.attachment, body, 0.0, settle, n_sub, fps, damping=lambda t: 2.0 if t < settle - 1.0 else 0.0)
        x_set, v_set = cl.x.clone(), cl.v.clone()
        hem = np.array(G.piece.rings['outer'])
        inner = np.array(G.piece.rings['inner'])
        uvh = G.piece.vertex_uv[hem]
        theta_p = np.arctan2(uvh[:, 1], uvh[:, 0])
        r_w = G.params['waist_radius_m']
        L = G.params['length_m']
        for rpm in rpms:
            om = rpm * 2 * math.pi / 60
            cl.x.copy_(x_set)
            cl.v.copy_(v_set)
            sp = SpinPose(p0, q0, om, t_up=t_up)
            acc = []

            def cb(t, c):
                if t > t_up + hold - 2.0 + 1e-9:
                    x = c.x.cpu().numpy()
                    zc = x[inner, 2].mean()
                    xy = x[hem, :2]
                    R = np.hypot(xy[:, 0], xy[:, 1])
                    ang = np.arctan2(xy[:, 1], xy[:, 0]) - theta_p - sp.angle(t)
                    lag = math.degrees(math.atan2(np.sin(ang).mean(), np.cos(ang).mean()))
                    acc.append((R.mean(), zc - x[hem, 2].mean(), lag, R.min(), R.max()))
            run(cl, sp, G.attachment, body, 0.0, t_up + hold, n_sub, fps, callback=cb)
            a = np.array(acc)
            Rm, drop, lag = a[:, 0].mean(), a[:, 1].mean(), a[:, 2].mean()
            sim = dict(hem_radius=float(Rm), drop=float(drop), flare_deg=math.degrees(math.atan2(Rm - r_w, drop)),
                       lag_deg=float(-lag), hem_radius_std_in_time=float(a[:, 0].std()),
                       hem_radius_min=float(a[:, 3].mean()), hem_radius_max=float(a[:, 4].mean()))
            phi_c = cone_estimate(om, world.gravity, r_w, L)
            est_cone = dict(flare_deg=math.degrees(phi_c), hem_radius=r_w + L * math.sin(phi_c))
            nd = chain_estimate(om, world.gravity, world.air_density, r_w, L, mat.areal_density, drag=False)
            wd = chain_estimate(om, world.gravity, world.air_density, r_w, L, mat.areal_density, drag=True)
            r = dict(world=world.name, fabric=fabric, rpm=rpm, omega=om, h=h, simulated=sim, rigid_cone_no_drag=est_cone,
                     chain_no_drag=nd, chain_with_drag=wd, wall_s=time.time() - t0)
            out.append(r)
            print(json.dumps(r), flush=True)
    return out


# ------------------------------------------------------------------------------------------------ 4. convergence
def convergence_cases():
    res = dict(cantilever_h=[], cantilever_dt=[], terminal_dt=[])
    m = FABRICS['silk_georgette']
    for world in WORLDS:
        for n_el in ((6, 12, 24) if world is MOON else (6, 12)):
            r = cantilever(m, world, n_el=n_el)
            res['cantilever_h'].append({k: v for k, v in r.items() if not k.startswith(('profile', 'elastica'))})
            print(json.dumps(res['cantilever_h'][-1]), flush=True)
    for f in (2.0, 1.0, 0.5):
        r = cantilever(m, MOON, n_el=6, dt_factor=f)
        res['cantilever_dt'].append({k: v for k, v in r.items() if not k.startswith(('profile', 'elastica'))})
        print(json.dumps(res['cantilever_dt'][-1]), flush=True)
    mc = FABRICS['silk_chiffon']
    for dt in (1 / 1200, 1 / 2400, 1 / 4800):
        ts, vn, vt, _ = terminal(mc, EARTH, dt=dt)
        an = math.sqrt(2 * mc.areal_density * EARTH.gravity / (EARTH.air_density * 1.17))
        res['terminal_dt'].append(dict(dt=dt, rel_error=float(-vn[-20:].mean() / an - 1)))
        print(json.dumps(res['terminal_dt'][-1]), flush=True)
    return res


# ------------------------------------------------------------------------------------------------------ plots
def plots(R, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    col = {'moon': '#3b6fb6', 'earth': '#c8553d'}
    if 'terminal' in R:
        fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
        for r in R['terminal']:
            if r['case'] != 'face-on':
                continue
            k = 0 if r['world'] == 'moon' else 1
            ax[k].plot(r['t'], r['v'], lw=1.4, label=r['fabric'].replace('silk_', ''))
            ax[k].axhline(r['analytic'], color='k', lw=0.6, ls=':')
        for k, wn in enumerate(('moon', 'earth')):
            ax[k].set_title(f'Face-on fall, {wn} (dotted: sqrt(2 m g / rho C_N))', fontsize=9)
            ax[k].set_xlabel('time (s)')
            ax[k].set_ylabel('falling speed (m/s)')
            ax[k].legend(fontsize=7, frameon=False)
        fig.tight_layout()
        fig.savefig(Path(out) / 'terminal_speed.png', dpi=130)
        plt.close(fig)
    if 'cantilever' in R:
        fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
        for r in R['cantilever']:
            c = col[r['world']]
            ax[0].plot(np.array(r['elastica_x']) * 100, np.array(r['elastica_z']) * 100, color=c, lw=0.8, alpha=0.6)
            ax[0].plot(np.array(r['profile_x']) * 100, np.array(r['profile_z']) * 100, 'o', color=c, ms=2.5)
            ax[1].plot(r['c_pred'] * 100, r['c_sim'] * 100, 'o', color=c, label=r['world'])
        lim = [0, max(r['c_pred'] for r in R['cantilever']) * 110]
        ax[1].plot(lim, lim, 'k:', lw=0.7)
        ax[0].set_title('Cantilever l = 2c: simulated (dots) and heavy elastica (lines)', fontsize=9)
        ax[0].set_xlabel('x (cm)')
        ax[0].set_ylabel('z (cm)')
        ax[0].set_aspect('equal')
        ax[1].set_title('Bending length: Peirce from simulation vs (B/w)^(1/3)', fontsize=9)
        ax[1].set_xlabel('predicted c (cm)')
        ax[1].set_ylabel('simulated c (cm)')
        h, l = ax[1].get_legend_handles_labels()
        ax[1].legend(dict(zip(l, h)).values(), dict(zip(l, h)).keys(), fontsize=8, frameon=False)
        fig.tight_layout()
        fig.savefig(Path(out) / 'cantilever.png', dpi=130)
        plt.close(fig)
    if 'spin' in R:
        fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
        for wn in ('moon', 'earth'):
            rs = [r for r in R['spin'] if r['world'] == wn]
            if not rs:
                continue
            rpm = [r['rpm'] for r in rs]
            c = col[wn]
            ax[0].plot(rpm, [r['simulated']['flare_deg'] for r in rs], 'o-', color=c, label=f'{wn}: simulated')
            ax[0].plot(rpm, [r['chain_with_drag']['flare_deg'] for r in rs], '--', color=c, label=f'{wn}: chain, drag')
            ax[0].plot(rpm, [r['chain_no_drag']['flare_deg'] for r in rs], ':', color=c, label=f'{wn}: chain, no drag')
            ax[0].plot(rpm, [r['rigid_cone_no_drag']['flare_deg'] for r in rs], 'x', color=c, label=f'{wn}: rigid cone')
            ax[1].plot(rpm, [r['simulated']['hem_radius'] for r in rs], 'o-', color=c, label=f'{wn}: simulated')
            ax[1].plot(rpm, [r['chain_with_drag']['hem_radius'] for r in rs], '--', color=c, label=f'{wn}: chain, drag')
        ax[0].set_xlabel('turn rate (rpm)')
        ax[0].set_ylabel('flare angle, waist to mean hem (deg)')
        ax[1].set_xlabel('turn rate (rpm)')
        ax[1].set_ylabel('mean hem radius (m)')
        ax[0].legend(fontsize=6.5, frameon=False)
        ax[1].legend(fontsize=6.5, frameon=False)
        ax[0].set_title('Full circle skirt (chiffon) in a steady turn', fontsize=9)
        fig.tight_layout()
        fig.savefig(Path(out) / 'spin_flare.png', dpi=130)
        plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--only', default='terminal,cantilever,spin,convergence')
    ap.add_argument('--spin-h', type=float, default=0.02)
    a = ap.parse_args(argv)
    torch.set_num_threads(3)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    path = out / 'validation.json'
    R = json.loads(path.read_text()) if path.exists() else {}
    for part in a.only.split(','):
        t0 = time.time()
        if part == 'terminal':
            R['terminal'] = terminal_cases(out)
        elif part == 'cantilever':
            R['cantilever'] = cantilever_cases()
        elif part == 'spin':
            R['spin'] = spin_runs(h=a.spin_h)
        elif part == 'convergence':
            R['convergence'] = convergence_cases()
        elif part == 'spin_satin':
            R['spin_satin'] = spin_runs(rpms=(15, 20), worlds=(MOON,), fabric='silk_satin_heavy')
        elif part == 'spin_vacuum':
            R['spin_vacuum'] = spin_runs(rpms=(15, 20), worlds=(World('moon-vacuum', MOON.gravity, 0.0),))
        elif part == 'spin_convergence':
            R['spin_convergence'] = [r for hh in (0.03, 0.015) for r in spin_runs(rpms=(15,), h=hh, worlds=(MOON,))]
        R.setdefault('wall_s', {})[part] = time.time() - t0
        path.write_text(json.dumps(R, indent=1))
        plots(R, out)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
