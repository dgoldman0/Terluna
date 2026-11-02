"""Trial runs of fabrics and cuts in the cloth solver, on the Moon and in an Earth control: the cases the screen
(screen.py) can only estimate.

    python -m humanities.dance.garment.trials OUTDIR [--only drape,turn,stream] [--fabrics a,b] [--fast]

1. Drape (Cusick's test). A disc of fabric 30 cm across, clamped between discs 18 cm across, hangs under its own
   weight. The drape coefficient is the share of the free ring that its shadow covers, (A_shadow - A_disc) /
   (A_specimen - A_disc): 1 for a fabric that stays flat, near 0 for one that hangs straight down. The nodes are the
   folds around the hem. A static test: the air only damps the fall.
2. Turn. Circle skirts (full and half circles, with and without a weighted hem) on the body turned steadily, as in
   validate.py's steady turn: the hanging shape after dressing, then the flare and the hem's lag at each rate of turn.
   A steady state, in still air (no co-rotation), so the drag on folds is an upper bound.
3. Stream. A panel pinned along its short edge to a rod that sets off horizontally and travels at a steady speed: the
   angle the panel trails below the horizontal and how its far edge flutters. The mean angle is a steady quantity; the
   flutter's amplitude and frequency depend on the air's added mass, which the solver does not carry (README,
   shortcomings), so they are indications only.
Results go to OUTDIR/trials.json, one record per case, written as each case finishes.
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
from .cloth import EARTH, MOON, Cloth, Settings
from .garments import skirt, waist_tethers
from .materials import FABRICS
from .patterns import annulus, rectangle
from .validate import TEST_PATTERN

WORLDS = {'moon': MOON, 'earth': EARTH}
CUSICK = dict(specimen_radius_m=0.15, support_radius_m=0.09)


def _save(out, R):
    p = Path(out) / 'trials.json'
    p.write_text(json.dumps(R, indent=1))


def _shape(out, tag, faces, verts):
    """Keep a case's final cloth for figures: faces and vertex positions, float32."""
    d = Path(out) / 'shapes'
    d.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(d / f'{tag}.npz', faces=np.asarray(faces, np.int32), verts=np.asarray(verts, np.float32))


def hem_outline(x, hem, centre):
    """Projected hem: radius against angle about the centre, the shadow's area (shoelace) and whether the outline
    goes once round in order (so that the shoelace area is the shadow's)."""
    xy = x[hem, :2] - centre[:2]
    ang = np.unwrap(np.arctan2(xy[:, 1], xy[:, 0]))
    monotone = bool(np.all(np.diff(ang) > 0) or np.all(np.diff(ang) < 0))
    area = 0.5 * abs(np.dot(xy[:, 0], np.roll(xy[:, 1], -1)) - np.dot(xy[:, 1], np.roll(xy[:, 0], -1)))
    return np.hypot(xy[:, 0], xy[:, 1]), area, monotone


def count_nodes(R, min_prominence):
    """Folds round a closed hem: peaks of its projected radius, counted with wrap-around."""
    from scipy.signal import find_peaks
    n = len(R)
    ext = np.r_[R, R, R]
    pk, _ = find_peaks(ext, prominence=min_prominence)
    return int(((pk >= n) & (pk < 2 * n)).sum())


# --------------------------------------------------------------------------------------------------------- 1. drape
def drape(mat, world, h=0.005, seed=0, T=None, damping_per_s=None, out=None):
    """Cusick's drape test; returns the drape coefficient and nodes with their history in time."""
    rs, Rs = CUSICK['support_radius_m'], CUSICK['specimen_radius_m']
    pc = annulus(rs - 2 * h, Rs, h, seed=seed)
    uv = pc.vertex_uv
    r = np.hypot(uv[:, 0], uv[:, 1])
    pinned = np.where(r <= rs + 1e-9)[0]
    X0 = np.c_[uv, np.full(pc.n_vertices, 1.0)]
    # small steps keep the bending's static error small at this mesh (validate.cantilever's step rule)
    dt = min(1 / 2400, 0.0141 * h * h / math.sqrt(mat.bending_rigidity / mat.areal_density))
    st = Settings(dt=dt, self_collision=True)
    cl = Cloth(pc, X0, mat, world, st, pinned=pinned)
    cl.set_targets(torch.as_tensor(X0[pinned], dtype=torch.float32, device=cl.dev))
    free_len = Rs - rs
    t_unit = math.sqrt(free_len / world.gravity)
    T = T or 25 * t_unit
    cl.damp.fill_(damping_per_s if damping_per_s is not None else 1.0 / t_unit)
    hem = np.array(pc.rings['outer'])
    a_disc = math.pi * rs ** 2
    a_spec = math.pi * Rs ** 2
    n = int(T / dt)
    every = max(1, n // 25)
    hist = []
    t0 = time.time()
    for k in range(n):
        cl.substep(dt)
        if (k + 1) % every == 0 or k == n - 1:
            x = cl.x.cpu().numpy()
            R, area, mono = hem_outline(x, hem, np.zeros(3))
            hist.append(dict(t=round((k + 1) * dt, 3), dc=round((area - a_disc) / (a_spec - a_disc), 4),
                             nodes=count_nodes(R, 0.003), monotone=mono))
    x = cl.x.cpu().numpy()
    R, area, mono = hem_outline(x, hem, np.zeros(3))
    ke = cl.energy()['kinetic']
    if out is not None:
        _shape(out, f'drape_{mat.name}_{world.name}_s{seed}', pc.F, x)
    return dict(fabric=mat.name, world=world.name, h=h, seed=seed, dt=dt, T=T, substeps=n,
                drape_coefficient=round((area - a_disc) / (a_spec - a_disc), 4), nodes=count_nodes(R, 0.003),
                outline_in_order=mono, hem_radius_min_m=round(float(R.min()), 4),
                hem_radius_max_m=round(float(R.max()), 4), hem_drop_m=round(float(1.0 - x[hem, 2].mean()), 4),
                kinetic_J=ke, history=hist, wall_s=round(time.time() - t0, 1))


# ---------------------------------------------------------------------------------------------------------- 2. turn
def turn(mat, world, cut=2 * math.pi, hem_kg_m=0.0, rpms=(15.0, 27.0), h=0.02, settle=None, t_up=2.0, hold=5.0,
         n_sub=40, fps=60.0, out=None, tethered=None):
    """The skirt is dressed and left to hang for settle seconds (5 s on Earth, longer by the gravity clock on the
    Moon, where a skirt laid flat takes that much longer to fall), then turned at each rate from the same hanging
    state."""
    settle = settle if settle is not None else 5.0 * math.sqrt(EARTH.gravity / world.gravity)
    body = Body()
    mp = MotionPose(TEST_PATTERN, body, t_raised=1.0, t_lower=1.0, t_hold=1.0)
    p0, q0 = mp(-mp.pre)
    G = skirt(body, (p0, q0), cut=cut, h=h, hem_weight_kg_m=hem_kg_m)
    # a weighted hem's pull needs tethers to the waistband (the solver's one stretch pass cannot carry it up light cloth)
    tethered = hem_kg_m > 0 if tethered is None else tethered
    cl = Cloth(G.piece, G.X0, mat, world, Settings(dt=1.0 / (fps * n_sub)), pinned=G.pinned, extra_mass=G.extra_mass,
               tethers=waist_tethers(G) if tethered else None)
    still = SpinPose(p0, q0, 0.0)
    t0 = time.time()
    run(cl, still, G.attachment, body, 0.0, settle, n_sub, fps, damping=lambda t: 2.0 if t < settle - 1.0 else 0.0)
    x_set, v_set = cl.x.clone(), cl.v.clone()
    hem = np.array(G.piece.rings['outer'])
    inner = np.array(G.piece.rings['inner'])
    r_w = G.params['waist_radius_m']
    L = G.params['length_m']
    x = cl.x.cpu().numpy()
    c = x[inner].mean(0)
    R, _, _ = hem_outline(x, hem, c)
    rest = dict(hem_radius_m=round(float(R.mean()), 3), hem_drop_m=round(float(c[2] - x[hem, 2].mean()), 3),
                nodes=count_nodes(R, 0.01), flare_deg=round(math.degrees(math.atan2(R.mean() - r_w,
                                                                                    c[2] - x[hem, 2].mean())), 1))
    psi = np.mod(np.arctan2(G.piece.vertex_uv[hem, 1], G.piece.vertex_uv[hem, 0]), 2 * np.pi)
    theta_p = math.pi + psi * 2 * math.pi / cut
    tag = f'{mat.name}_{world.name}_cut{round(math.degrees(cut))}_hem{round(hem_kg_m * 1000)}'
    if out is not None:
        _shape(out, f'turn_{tag}_rest', G.piece.F, x)
    res = dict(fabric=mat.name, world=world.name, cut_deg=round(math.degrees(cut)), hem_kg_m=hem_kg_m, h=h,
               tethered=bool(tethered),
               params={k: v for k, v in G.params.items() if not isinstance(v, np.ndarray)}, at_rest=rest, turns=[])
    for rpm in rpms:
        om = rpm * 2 * math.pi / 60
        cl.x.copy_(x_set)
        cl.v.copy_(v_set)
        sp = SpinPose(p0, q0, om, t_up=t_up)
        acc = []

        def cb(t, c_):
            if t > t_up + hold - 2.0 + 1e-9:
                xx = c_.x.cpu().numpy()
                zc = xx[inner, 2].mean()
                xy = xx[hem, :2]
                Rr = np.hypot(xy[:, 0], xy[:, 1])
                ang = np.arctan2(xy[:, 1], xy[:, 0]) - theta_p - sp.angle(t)
                lag = math.degrees(math.atan2(np.sin(ang).mean(), np.cos(ang).mean()))
                acc.append((Rr.mean(), zc - xx[hem, 2].mean(), lag, Rr.min(), Rr.max()))
        run(cl, sp, G.attachment, body, 0.0, t_up + hold, n_sub, fps, callback=cb)
        a = np.array(acc)
        Rm, drop, lag = a[:, 0].mean(), a[:, 1].mean(), a[:, 2].mean()
        if out is not None:
            _shape(out, f'turn_{tag}_rpm{rpm:g}', G.piece.F, cl.x.cpu().numpy())
        res['turns'].append(dict(rpm=rpm, hem_radius_m=round(float(Rm), 3), drop_m=round(float(drop), 3),
                                 flare_deg=round(math.degrees(math.atan2(Rm - r_w, drop)), 1),
                                 lag_deg=round(float(-lag), 1), hem_radius_spread_m=round(float((a[:, 4] - a[:, 3]).mean()), 3),
                                 hem_radius_std_in_time_m=round(float(a[:, 0].std()), 4)))
    res['wall_s'] = round(time.time() - t0, 1)
    res['strain'] = cl.strain_stats()
    return res


# -------------------------------------------------------------------------------------------------------- 3. stream
def stream(mat, world, U, length=1.2, width=0.4, h=0.02, settle=None, ramp=1.0, hold=5.0, n_sub=40, fps=60.0,
           out=None):
    """A panel hanging from a horizontal rod (its short edge, across the travel) that sets off along +x and holds
    speed U. Returns the mean trailing angle below the horizontal over the last 3 s and the far edge's flutter."""
    pc = rectangle(width, length, h, warp=(0.0, 1.0))
    uv = pc.vertex_uv
    X0 = np.c_[np.zeros(pc.n_vertices), uv[:, 0] - width / 2, 2.0 - uv[:, 1]]       # hanging straight down
    pinned = np.array(pc.rings['start'])
    far = np.array(pc.rings['end'])
    dt = 1.0 / (fps * n_sub)
    cl = Cloth(pc, X0, mat, world, Settings(dt=dt, self_collision=True), pinned=pinned)
    P0 = torch.as_tensor(X0[pinned], dtype=torch.float32, device=cl.dev)
    settle = settle if settle is not None else 6 * math.sqrt(length / world.gravity)

    def xrod(t):                       # smooth start: speed rises as smoothstep over the ramp
        if t <= 0:
            return 0.0
        if t < ramp:
            x = t / ramp
            return U * ramp * (x ** 3 - 0.5 * x ** 4)
        return U * ramp * 0.5 + U * (t - ramp)

    t0w = time.time()
    cl.set_targets(P0)
    cl.damp.fill_(2.0)
    for _ in range(int(settle * fps * n_sub)):
        cl.substep(dt)
    cl.damp.fill_(0.0)
    rec = []
    n_frames = int((ramp + hold) * fps)
    off = torch.zeros(3, dtype=torch.float32, device=cl.dev)
    for f in range(n_frames):
        ta, tb = f / fps, (f + 1) / fps
        T0 = P0 + torch.tensor([xrod(ta), 0, 0], dtype=torch.float32, device=cl.dev)
        T1 = P0 + torch.tensor([xrod(tb), 0, 0], dtype=torch.float32, device=cl.dev)
        cl.set_frame(T0=T0, T1=T1)
        for k in range(n_sub):
            cl.substep(dt, k / n_sub, (k + 1) / n_sub)
        if tb > ramp + hold - 3.0:
            x = cl.x.cpu().numpy()
            rod = np.array([xrod(tb), 0.0, 2.0])
            tip = x[far].mean(0)
            rec.append((tb, rod[0] - tip[0], rod[2] - tip[2], tip[2]))
    if out is not None:
        _shape(out, f'stream_{mat.name}_{world.name}_U{U:g}', pc.F, cl.x.cpu().numpy())
    a = np.array(rec)
    ang = np.degrees(np.arctan2(a[:, 2], a[:, 1]))
    z = a[:, 3] - a[:, 3].mean()
    fr = np.fft.rfftfreq(len(z), 1 / fps)
    spec = np.abs(np.fft.rfft(z))
    f_peak = float(fr[1:][np.argmax(spec[1:])]) if len(z) > 4 else None
    return dict(fabric=mat.name, world=world.name, U_m_s=U, panel_m=[length, width], h=h,
                trail_deg_mean=round(float(ang.mean()), 1), trail_deg_std=round(float(ang.std()), 1),
                far_edge_flutter_rms_m=round(float(z.std()), 4), flutter_peak_hz=None if f_peak is None else round(f_peak, 2),
                strain=cl.strain_stats(), wall_s=round(time.time() - t0w, 1))


# The turn trials: cut (radians), fabric and weighted hem (kg/m) on each world; the Earth controls take two fabrics.
TURN_CASES = {'moon': [(2 * math.pi, 'silk_chiffon', 0.0), (math.pi, 'silk_chiffon', 0.0),
                       (2 * math.pi, 'silk_chiffon', 0.1), (2 * math.pi, 'silk_chiffon', 0.2),
                       (2 * math.pi, 'silk_satin_heavy', 0.0), (math.pi, 'silk_satin_heavy', 0.0),
                       (2 * math.pi, 'wool_plain', 0.0), (math.pi, 'wool_plain', 0.0),
                       (2 * math.pi, 'silk_velvet', 0.0), (2 * math.pi, 'nylon_tulle', 0.0)],
              'earth': [(2 * math.pi, 'silk_chiffon', 0.0), (2 * math.pi, 'wool_plain', 0.0)]}
# Weighted hems, each tethered to the waistband: full circles of chiffon, and half circles whose wind-up index the
# hem brings to 1.1-2.0 (screen.wind_up_index), on the Moon.
HEM_CASES = [(2 * math.pi, 'silk_chiffon', 0.1), (2 * math.pi, 'silk_chiffon', 0.2), (math.pi, 'silk_chiffon', 0.05),
             (math.pi, 'silk_chiffon', 0.1), (math.pi, 'silk_satin_heavy', 0.05)]
STREAM_FABRICS = ('silk_chiffon', 'silk_satin_heavy', 'wool_plain')
# Lee et al. (2011) measured a silk satin of 69 g/m2 at a drape coefficient of 21.9% (sources.json, measured_drape);
# Earth drapes at these bending rigidities place it.
SATIN_LEE = dict(areal_density=0.069, drape=0.219, bending_grid=(0.5e-6, 0.75e-6, 1.0e-6, 2.0e-6, 4.0e-6, 8.0e-6))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--only', default='drape,satin,turn,stream,hem')
    ap.add_argument('--fabrics', default='')
    ap.add_argument('--fast', action='store_true', help='one short case of each part, to check the harness')
    a = ap.parse_args(argv)
    torch.set_num_threads(1)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    path = out / 'trials.json'
    R = json.loads(path.read_text()) if path.exists() else {}
    names = a.fabrics.split(',') if a.fabrics else list(FABRICS)
    keys = ('fabric', 'world', 'cut_deg', 'hem_kg_m', 'U_m_s', 'seed', 'bending_rigidity')
    for part in a.only.split(','):
        R.setdefault(part, [])
        done = {json.dumps({k: r[k] for k in r if k in keys}, sort_keys=True) for r in R[part]}

        def todo(**key):
            return json.dumps(key, sort_keys=True) not in done

        def keep(r):
            R[part].append(r)
            _save(out, R)
            print(json.dumps({k: v for k, v in r.items() if k not in ('history', 'params')}), flush=True)

        if part == 'drape':
            for name in (names[:1] if a.fast else names):
                for w in WORLDS.values():
                    if todo(fabric=name, world=w.name, seed=0):
                        keep(drape(FABRICS[name], w, T=0.3 if a.fast else None, out=out))
        elif part == 'satin':
            base = FABRICS['silk_satin_heavy']
            for B in (SATIN_LEE['bending_grid'][:1] if a.fast else SATIN_LEE['bending_grid']):
                mat = type(base)(**{**base.__dict__, 'name': 'silk_satin_lee', 'areal_density': SATIN_LEE['areal_density'],
                                    'bending_rigidity': B})
                if todo(fabric='silk_satin_lee', world='earth', seed=0, bending_rigidity=B):
                    r = drape(mat, EARTH, T=0.3 if a.fast else None)
                    r['bending_rigidity'] = B
                    keep(r)
        elif part == 'turn':
            for wn, cases in TURN_CASES.items():
                for cut, name, hm in (cases[:1] if a.fast else cases):
                    if name in names and todo(fabric=name, world=wn, cut_deg=round(math.degrees(cut)), hem_kg_m=hm):
                        keep(turn(FABRICS[name], WORLDS[wn], cut=cut, hem_kg_m=hm, out=out,
                                  **(dict(rpms=(27.0,), settle=1.0, hold=2.5, t_up=0.5) if a.fast else {})))
                if a.fast:
                    break
        elif part == 'hem':
            for cut, name, hm in (HEM_CASES[:1] if a.fast else HEM_CASES):
                if name in names and todo(fabric=name, world='moon', cut_deg=round(math.degrees(cut)), hem_kg_m=hm):
                    keep(turn(FABRICS[name], MOON, cut=cut, hem_kg_m=hm, out=out, tethered=True,
                              **(dict(rpms=(27.0,), settle=1.0, hold=2.5, t_up=0.5) if a.fast else {})))
        elif part == 'stream':
            for name in (STREAM_FABRICS[:1] if a.fast else STREAM_FABRICS):
                for w in WORLDS.values():
                    for U in ((2.0,) if a.fast else (1.0, 2.0, 3.0)):
                        if name in names and todo(fabric=name, world=w.name, U_m_s=U):
                            keep(stream(FABRICS[name], w, U, out=out, **(dict(hold=3.2, settle=0.5) if a.fast else {})))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
