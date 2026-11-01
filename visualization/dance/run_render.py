"""Render a dance motion file to stills and an MP4 with the skinned dancer (Blender 4.0.2, EEVEE).

    python visualization/dance/run_render.py MOTION.npz --out OUT.mp4 [--garment G.npz ...]
           [--preset test|final] [--stills 2.0,4.4,6.5] [--stills-only] [--frames 1:120] [--rebuild]

MOTION.npz follows humanities/dance/motion_io.py (schema terluna.humanities.dance-motion/1). Each --garment file holds
`faces` (F, 3) int and `verts` (N_frames, V, 3) in the dance world frame at the motion's rate (or at its own `fps`).
The rigged dancer (.blend) is built once from the CC0 base mesh and body.json and reused; --rebuild remakes it, which
also rewrites landmark_check.json. Presets: test = 1280x720, 24 samples, soft shadows; video = 1920x1080, 32
samples, filtered (not sampled) soft shadows, for whole sequences; final = 1920x1080, 80 samples, soft shadows,
for stills. Motion blur is off (EEVEE's returned black frames in background renders here).
Blender runs with 3 CPU threads; EEVEE renders on the GPU.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BLENDER = '/usr/bin/blender'
DRIVE = Path('/media/projectspace/terluna-research/humanities')
DANCER = DRIVE / 'dance-runs' / 'renders' / 'dancer_rigged.blend'
REPORT = DRIVE / 'agents' / 'render'


def blender(args, threads=3):
    cmd = [BLENDER, '-b', '--factory-startup', '--threads', str(threads), '--python', str(HERE / 'render.py'), '--']
    return subprocess.run(cmd + [str(a) for a in args], cwd=ROOT).returncode


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('motion', help='motion .npz (motion_io format)')
    ap.add_argument('--out', required=True, help='output .mp4; stills are written beside it as <stem>_tSS.SS.png')
    ap.add_argument('--garment', action='append', default=[], help='cloth sequence .npz (repeatable)')
    ap.add_argument('--preset', choices=['test', 'video', 'final'], default='test')
    ap.add_argument('--stills', default='', help='comma-separated times in seconds')
    ap.add_argument('--stills-only', action='store_true')
    ap.add_argument('--frames', default='', help='render-frame range a:b for the video (1-based, inclusive)')
    ap.add_argument('--cam-az', type=float, default=None,
                    help='camera bearing from the +x axis, degrees (default: 30 for phrases in place, square to the '
                         'travel for travelling phrases)')
    ap.add_argument('--cam-dist', type=float, default=None, help='camera distance, m (default 5.5 in place, 6.8 travelling)')
    ap.add_argument('--dancer', default=str(DANCER), help='rigged dancer .blend (built if missing)')
    ap.add_argument('--report-dir', default=str(REPORT))
    ap.add_argument('--rebuild', action='store_true')
    ap.add_argument('--save-scene', default='', help='also save the assembled scene .blend here')
    a = ap.parse_args(argv)
    if a.rebuild or not Path(a.dancer).exists():
        rc = blender(['build', '--blend', a.dancer, '--report-dir', a.report_dir])
        if rc:
            return rc
    args = ['render', '--blend', a.dancer, '--motion', Path(a.motion).resolve(), '--out', Path(a.out).resolve(),
            '--preset', a.preset]
    if a.cam_az is not None:
        args += ['--cam-az', a.cam_az]
    if a.cam_dist is not None:
        args += ['--cam-dist', a.cam_dist]
    for g in a.garment:
        args += ['--garment', Path(g).resolve()]
    if a.stills:
        args += ['--stills', a.stills]
    if a.stills_only:
        args += ['--stills-only']
    if a.frames:
        args += ['--frames', a.frames]
    if a.save_scene:
        args += ['--save-scene', Path(a.save_scene).resolve()]
    return blender(args)


if __name__ == '__main__':
    sys.exit(main())
