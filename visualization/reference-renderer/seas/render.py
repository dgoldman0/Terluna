"""Render the sea-appearance study's four coasts through the month: perspective frames from an eye 2 m above
the water, spectral light carried as CIE XYZ.

    NUMBA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 <numba environment>/bin/python \\
        visualization/reference-renderer/seas/render.py [--coasts ...] [--moments ...] [--width 1920] [--spp 8]

Each frame is rendered in blocks of rows, each saved as it completes, so an interrupted run resumes where it
stopped. Finished frames (radiance in cd/m2 as XYZ, the luminance's sampling error and what each pixel met)
go to research/runs/sea_appearance/renderings/; display images are made by frames.py.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import kernels as k  # noqa: E402
import scene  # noqa: E402
import stars  # noqa: E402

OUT = scene.ROOT / "research/runs/sea_appearance/renderings"
ROWS_PER_BLOCK = 30
SOURCES = [HERE / "kernels.py", HERE / "scene.py", HERE / "render.py", HERE / "stars.py"]


def identity(record, spp_side, seed):
    """What fixes a frame: its scene record, sampling and the renderer's source."""
    body = dict(record=record, spp_side=spp_side, seed=seed,
                sources={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in SOURCES})
    return hashlib.sha256(json.dumps(body, sort_keys=True, default=str).encode()).hexdigest()


def frame_name(coast, moment):
    return f"{coast.replace(' ', '_').lower()}__{moment.replace(' ', '_')}"


def render_frame(common, coast, moment, width, height, spp_side, seed, fov, out=OUT):
    name = frame_name(coast, moment)
    final = out / f"{name}.npz"
    started = time.monotonic()
    camera, arrays, record, geometry, light = scene.build(common, coast, moment, width, height, fov, seed)
    ident = identity(record, spp_side, seed)
    if final.exists():
        with np.load(final, allow_pickle=False) as z:
            if json.loads(str(z["record"])).get("identity") == ident:
                print(name, "up to date", flush=True)
                return final
    work = out / "partial" / name
    work.mkdir(parents=True, exist_ok=True)
    (work / "identity.txt").write_text(ident)
    print(name, "scene built in", round(time.monotonic() - started, 1), "s;", json.dumps(record["sea"]["spectrum"]),
          flush=True)
    a = arrays
    for row0 in range(0, height, ROWS_PER_BLOCK):
        row1 = min(height, row0 + ROWS_PER_BLOCK)
        block = work / f"rows_{row0:05d}.npz"
        if block.exists():
            continue
        t0 = time.monotonic()
        mean, error, kind = k.render(width, height, row0, row1, spp_side, seed, camera, a["params"], a["sea_h"],
                                     a["sea_blocks"], a["lean"], a["lean_offsets"], a["lean_sizes"],
                                     a["lean_spacing"], a["lean_levels"], a["residual"], a["ter_h"], a["ter_blocks"],
                                     a["tex"], a["sky"], a["sky_el"], a["lut_s"], a["lut_t"], a["lut_tb"],
                                     a["lut_el"], a["lut_d"], a["sun"], a["earth"], a["water"])
        temporary = block.with_suffix(".tmp.npz")
        np.savez(temporary, mean=mean, error=error, kind=kind)
        temporary.replace(block)
        print(f"  rows {row0}-{row1} in {time.monotonic() - t0:.1f} s", flush=True)
    blocks = sorted(work.glob("rows_*.npz"))
    parts = [np.load(b) for b in blocks]
    xyz = np.concatenate([p["mean"] for p in parts])
    error = np.concatenate([p["error"] for p in parts])
    kind = np.concatenate([p["kind"] for p in parts])
    star_layer, star_record = stars.layer(common, record, geometry, camera, arrays, width, height, fov)
    record.update(identity=ident, spp=spp_side ** 2, seconds=round(time.monotonic() - started, 1),
                  stars=star_record, renderer_sources={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                                      for p in SOURCES},
                  inputs=common.inputs)
    final.parent.mkdir(parents=True, exist_ok=True)
    temporary = final.with_suffix(".tmp.npz")
    np.savez_compressed(temporary, xyz=xyz.astype(np.float32), stars=star_layer.astype(np.float32),
                        error=error.astype(np.float32), kind=kind.astype(np.float16),
                        record=json.dumps(record))
    temporary.replace(final)
    for b in blocks:
        b.unlink()
    (work / "identity.txt").unlink()
    work.rmdir()
    print(name, "done in", record["seconds"], "s", flush=True)
    return final


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--coasts", nargs="*", default=list(scene.CAMERAS))
    parser.add_argument("--moments", nargs="*", default=["noon", "low Sun", "after sunset", "evening",
                                                         "deep evening", "darkest"])
    parser.add_argument("--width", type=int, default=1920)
    parser.add_argument("--height", type=int, default=1080)
    parser.add_argument("--spp-side", type=int, default=4, help="samples per pixel = this squared")
    parser.add_argument("--fov", type=float, default=65.0)
    parser.add_argument("--seed", type=int, default=20380207)
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()
    common = scene.Common()
    for coast in args.coasts:
        for moment in args.moments:
            render_frame(common, coast, moment, args.width, args.height, args.spp_side, args.seed, args.fov,
                         args.out)


if __name__ == "__main__":
    main()
