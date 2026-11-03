"""Render a dance motion file to stills and an MP4 with the skinned dancer (Blender 4.0.2, EEVEE).

    python visualization/dance/run_render.py MOTION.npz --out OUT.mp4 [--garment G.npz ...]
           [--preset draft|test|detail|video|final] [--stills 2.0,4.4,6.5] [--stills-preset final] [--stills-only]
           [--frames 1:120] [--cam-az DEG] [--cam-dist M] [--cam-aim M] [--cam-fstop N] [--rebuild] [--bodice R,G,B [--bodice-look matte|satin]] [--yoke R,G,B,A] [--sheer hashed|blend] [--samples N]
           [--bodice-shell STYLE[,NECKLINE] [--trim R,G,B]] [--light stage|soft]

MOTION.npz follows humanities/dance/motion_io.py (schema terluna.humanities.dance-motion/1). Each --garment file holds
`faces` (F, 3) int and `verts` (N_frames, V, 3) in the dance world frame at the motion's rate (or at its own `fps`).
The rigged dancer (.blend) is built once from the CC0 base mesh and body.json and reused; --rebuild remakes it, which
also rewrites landmark_check.json. Presets: draft = 854x480, 24 samples, soft shadows, the default for whole
videos until a complete video is near; test = 1280x720, 24 samples, soft shadows; video = 1920x1080, 32
samples, filtered (not sampled) soft shadows, for whole sequences; final = 1920x1080, 80 samples, soft shadows,
for stills (the default for --stills, which check detail at full size); detail = 1280x720 at final quality. Motion blur is off (EEVEE's returned black frames in background renders here).
Blender runs with 3 CPU threads; EEVEE renders on the GPU.

Costume options for design stills. --bodice R,G,B dresses the leotard zone as a bodice of that sRGB colour, in the
look --bodice-look names: matte (the default), a matte silk with a soft sheen in its own colour to sit with sheer
chiffon, or satin, the glossy satin of the first design stills. --yoke R,G,B,A lays a sheer chiffon yoke of that sRGB
colour over the skin above the bodice, from the bodice's neckline to a jewel neckline at the base of the neck and over
the shoulders to the shoulder point; A (0 to 1) is the share of the skin it veils seen face on, denser toward the
silhouette and at its hems. The yoke is off unless --yoke is given, and it needs a dancer built since 2 November (the
build stores its zone; rebuild an older dancer with --rebuild). Sheer garments render as one surface seen from both
sides; opaque ones are given 1.5 mm of thickness. --bodice-shell STYLE[,NECKLINE] makes the bodice a garment: a shell
over the torso, skinned as the body is, in a structured, ruched or pleated style with a sweetheart, vee or straight
neckline piped in the --trim colour (the README gives its construction). --cam-aim and --cam-fstop frame close stills:
--cam-aim 0.36 --cam-dist 2.1 --cam-fstop 8 takes the standing dancer from the hips to the head.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BLENDER = '/usr/bin/blender'
DRIVE = Path('/media/projectspace/terluna-research/humanities')
DANCER = DRIVE / 'dance-runs' / 'renders' / 'dancer_rigged.blend'
REPORT = DRIVE / 'agents' / 'render'


def gpu_env():
    """The environment for Blender, with OpenGL offloaded to the NVIDIA GPU when its driver is loaded. On a laptop with
    hybrid graphics Blender otherwise renders EEVEE through Mesa on the integrated GPU, which returned black or striped
    frames at 1080p and at 80 samples and took minutes over a still that the NVIDIA GPU renders in seconds."""
    env = dict(os.environ)
    if Path('/proc/driver/nvidia/version').exists():
        env.setdefault('__NV_PRIME_RENDER_OFFLOAD', '1')
        env.setdefault('__GLX_VENDOR_LIBRARY_NAME', 'nvidia')
    return env


def blender(args, threads=3):
    cmd = [BLENDER, '-b', '--factory-startup', '--threads', str(threads), '--python', str(HERE / 'render.py'), '--']
    return subprocess.run(cmd + [str(a) for a in args], cwd=ROOT, env=gpu_env()).returncode


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('motion', help='motion .npz (motion_io format)')
    ap.add_argument('--out', required=True, help='output .mp4; stills are written beside it as <stem>_tSS.SS.png')
    ap.add_argument('--garment', action='append', default=[], help='cloth sequence .npz (repeatable)')
    ap.add_argument('--preset', choices=['draft', 'test', 'detail', 'video', 'final'], default='draft')
    ap.add_argument('--stills-preset', choices=['draft', 'test', 'detail', 'video', 'final'], default='final',
                    help='size and quality of the stills (default final: 1920x1080, 80 samples)')
    ap.add_argument('--stills', default='', help='comma-separated times in seconds')
    ap.add_argument('--stills-only', action='store_true')
    ap.add_argument('--frames', default='', help='render-frame range a:b for the video (1-based, inclusive)')
    ap.add_argument('--cam-az', type=float, default=None,
                    help='camera bearing from the +x axis, degrees (default: 30 for phrases in place, square to the '
                         'travel for travelling phrases)')
    ap.add_argument('--cam-dist', type=float, default=None, help='camera distance, m (default 5.5 in place, 6.8 travelling)')
    ap.add_argument('--cam-aim', type=float, default=None,
                    help='height of the camera\'s aim above the centre of mass, m (default 0.06); the camera rises '
                         'with it, so --cam-aim 0.32 --cam-dist 1.9 frames the bodice level, hips to head')
    ap.add_argument('--cam-fstop', type=float, default=None,
                    help='aperture (default 2.8); a close still wants about 8 to keep the torso in focus')
    ap.add_argument('--dancer', default=str(DANCER), help='rigged dancer .blend (built if missing)')
    ap.add_argument('--report-dir', default=str(REPORT))
    ap.add_argument('--rebuild', action='store_true')
    ap.add_argument('--save-scene', default='', help='also save the assembled scene .blend here')
    ap.add_argument('--bodice', default='', help='sRGB r,g,b: dress the leotard as a bodice of this colour')
    ap.add_argument('--bodice-look', choices=['matte', 'satin'], default='matte',
                    help="the bodice's look: matte silk (default) or glossy satin")
    ap.add_argument('--yoke', default='',
                    help='sRGB r,g,b and strength a (0-1): a sheer chiffon yoke over the skin above the bodice (off '
                         'by default)')
    ap.add_argument('--sheer', choices=['hashed', 'blend'], default='hashed',
                    help='transparency of sheer cloth: hashed (dithered; smooth only at high sample counts) or blend '
                         '(smooth at any count; layers sorted by object)')
    ap.add_argument('--samples', type=int, default=0, help="render samples (0: the preset's)")
    ap.add_argument('--bodice-shell', default='',
                    help='STYLE[,NECKLINE]: make the bodice a shell over the torso (structured, ruched or pleated; '
                         'sweetheart, vee or straight neckline), coloured by --bodice')
    ap.add_argument('--trim', default='0.70,0.52,0.25', help="sRGB r,g,b of the bodice shell's piped edge")
    ap.add_argument('--light', choices=['stage', 'soft'], default='stage',
                    help='the lighting: stage (side light from coloured booms, the default) or soft (an even, stronger front light)')
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
    if a.cam_aim is not None:
        args += ['--cam-aim', a.cam_aim]
    if a.cam_fstop is not None:
        args += ['--cam-fstop', a.cam_fstop]
    for g in a.garment:
        args += ['--garment', Path(g).resolve()]
    if a.stills:
        args += ['--stills', a.stills, '--stills-preset', a.stills_preset]
    if a.stills_only:
        args += ['--stills-only']
    if a.frames:
        args += ['--frames', a.frames]
    if a.save_scene:
        args += ['--save-scene', Path(a.save_scene).resolve()]
    if a.bodice:
        args += ['--bodice', a.bodice, '--bodice-look', a.bodice_look]
    if a.yoke:
        args += ['--yoke', a.yoke]
    if a.sheer != 'hashed':
        args += ['--sheer', a.sheer]
    if a.samples:
        args += ['--samples', a.samples]
    if a.bodice_shell:
        args += ['--bodice-shell', a.bodice_shell, '--trim', a.trim]
    if a.light != 'stage':
        args += ['--light', a.light]
    return blender(args)


if __name__ == '__main__':
    sys.exit(main())
