"""Dress a garment on a motion and simulate it, writing the vertex sequence at the motion's frame rate.

    python -m humanities.dance.garment.drive skirt moon MOTION.npz OUT.npz [--fabric silk_chiffon] [--h 0.02]
    python -m humanities.dance.garment.drive sleeves earth MOTION.npz OUT.npz

Output (schema terluna.humanities.dance-garment/1): faces (F, 3) int32; verts (N, V, 3) float32 in the motion's world
frame at its fps, frame k at time[k]; time (N,); uv (V, 2) pattern coordinates; pinned (P,) int32; meta, a JSON
string with the garment, fabric, world (gravity, air density), solver settings, dressing sequence, timing, memory and
the evidence statement.
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
import torch

from .body import Body, MotionPose, run
from .cloth import EARTH, MOON, Cloth, Settings
from .garments import circle_skirt, water_sleeves
from .materials import FABRICS

SCHEMA = 'terluna.humanities.dance-garment/1'
WORLDS = {'moon': MOON, 'earth': EARTH}


def dressing_damping(pose: MotionPose, strong=2.0, gentle=0.5):
    """Settling damping (1/s): strong while the arms are raised, gentle while they come down, none in the motion."""
    def f(t):
        if t > 0:
            return 0.0
        tau = t + pose.pre
        return strong if tau < pose.t_raised else gentle
    return f


def simulate(garment: str, world: str, motion, out, fabric=None, h=0.02, substeps=40, t_raised=6.0, t_lower=1.5,
             t_hold=1.5, settings: Settings | None = None, t_end=None, quiet=False):
    torch.set_num_threads(3)
    W = WORLDS[world]
    body = Body()
    pose = MotionPose(str(motion), body, t_raised=t_raised, t_lower=t_lower, t_hold=t_hold)
    p0 = pose(-pose.pre)
    if garment == 'skirt':
        G = circle_skirt(body, p0, h=h)
        mat = FABRICS[fabric or 'silk_chiffon']
    elif garment == 'sleeves':
        G = water_sleeves(body, p0, h=h)
        mat = FABRICS[fabric or 'silk_georgette']
    else:
        raise ValueError(garment)
    s = settings or Settings()
    s.dt = 1.0 / (pose.fps * substeps)
    cloth = Cloth(G.piece, G.X0, mat, W, s, pinned=G.pinned)
    damp = dressing_damping(pose)
    torch.cuda.reset_peak_memory_stats() if torch.cuda.is_available() else None
    t_wall = time.time()
    log = []

    def cb(t, c):
        if abs(t * pose.fps - round(t * pose.fps)) < 1e-6 and int(round(t * pose.fps)) % int(pose.fps) == 0:
            st = c.strain_stats()
            x = c.x
            log.append(dict(t=round(t, 3), z_min=float(x[:, 2].min()), z_mean=float(x[:, 2].mean()),
                            speed_max=float(c.v.norm(dim=1).max()), max_stretch=st["max_stretch"], p99_stretch=st["p99_stretch"], mean_strain=st["mean_abs_strain"],
                            finite=bool(torch.isfinite(x).all())))
            if not quiet:
                print(json.dumps(log[-1]), flush=True)

    run(cloth, pose, G.attachment, body, -pose.pre, 0.0, substeps, pose.fps, damping=damp, callback=cb)
    t_pre = time.time() - t_wall
    first = cloth.x.detach().cpu().numpy().astype(np.float32)
    t1 = pose.duration if t_end is None else min(t_end, pose.duration)
    frames, times = run(cloth, pose, G.attachment, body, 0.0, t1, substeps, pose.fps, record=lambda t: True,
                        damping=damp, callback=cb)
    wall = time.time() - t_wall
    verts = np.concatenate([first[None], frames]).astype(np.float32)
    times = np.r_[0.0, times]
    mem = torch.cuda.max_memory_allocated() / 1e6 if torch.cuda.is_available() else None
    sim_time = pose.pre + t1
    meta = dict(
        schema=SCHEMA, garment=G.name, garment_params=G.params, world=W.name, gravity=W.gravity,
        air_density=W.air_density, fabric=mat.__dict__ if hasattr(mat, '__dict__') else str(mat),
        solver=dict(method='XPBD small steps, Jacobi with mass splitting, CUDA graph', substeps_per_frame=substeps,
                    settings=s.to_dict(), self_contact_distance=cloth.d_self, device=str(cloth.dev),
                    dtype='float32', pairs_max=cloth.stats['pairs_max'], pairs_truncated=cloth.stats['pairs_truncated']),
        dressing=dict(description='first frame held with the arms raised 80 deg to the side while the garment settles '
                                  'with numerical damping, then the arms lowered to the first frame and held',
                      t_raised=pose.t_raised, t_lower=pose.t_lower, t_hold=pose.t_hold,
                      damping_per_s=dict(raised=2.0, lowering_and_hold=0.5, motion=0.0)),
        motion=dict(path=str(motion), fps=pose.fps, frames=len(times), producer=pose.m['meta'].get('producer')),
        timing=dict(wall_s=wall, wall_dressing_s=t_pre, simulated_s=sim_time,
                    motion_seconds_per_compute_second=sim_time / wall, peak_gpu_memory_mb=mem),
        log=log,
        evidence='Simulated with stated fabric values (sources.json) and quasi-steady per-triangle aerodynamics '
                 'validated against analytic cases (report.md); no comparison with filmed fabric. The air carries no '
                 'added mass and no wake; layers in contact share one face to the air.',
        producer='humanities/dance/garment/drive.py')
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out, schema=SCHEMA, faces=G.piece.F.astype(np.int32), verts=verts, time=times,
                        uv=G.piece.vertex_uv.astype(np.float32), pinned=np.asarray(G.pinned, np.int32),
                        meta=json.dumps(meta))
    if not quiet:
        print(f'wrote {out}: {verts.shape}, wall {wall:.0f} s, {sim_time / wall:.3f} s of motion per s, peak {mem} MB')
    return meta


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('garment', choices=['skirt', 'sleeves'])
    ap.add_argument('world', choices=list(WORLDS))
    ap.add_argument('motion')
    ap.add_argument('out')
    ap.add_argument('--fabric', default=None, choices=list(FABRICS))
    ap.add_argument('--h', type=float, default=0.02)
    ap.add_argument('--substeps', type=int, default=40)
    ap.add_argument('--t-end', type=float, default=None)
    a = ap.parse_args(argv)
    simulate(a.garment, a.world, a.motion, a.out, fabric=a.fabric, h=a.h, substeps=a.substeps, t_end=a.t_end)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
