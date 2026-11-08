"""A trivial cloth sequence for testing the renderer's garment path: a short flared skirt carried by the pelvis.

    python visualization/dance/test_sheet.py MOTION.npz OUT.npz

The skirt is kinematic (a ring of panels hung from the pelvis frame, flaring with the pelvis's turning speed); it has no
cloth physics and exists only to exercise Mesh Cache playback until the garment solver's sequences arrive. Output:
`faces` (F, 3) int, `verts` (N, V, 3) in the dance world frame at the motion's rate, `fps`.
"""
from __future__ import annotations

import sys

import numpy as np


def quat_to_mat(q):
    q = q / np.linalg.norm(q, axis=-1, keepdims=True)
    w, x, y, z = np.moveaxis(q, -1, 0)
    return np.stack([np.stack([1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)], -1),
                     np.stack([2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)], -1),
                     np.stack([2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)], -1)], -2)


def main(src, out, n_around=48, n_down=8, length=0.30):
    with np.load(src, allow_pickle=False) as z:
        names = [str(s) for s in z['segments']]
        pos, quat, fps = z['seg_pos'], z['seg_quat'], float(z['fps'])
    i = names.index('pelvis')
    P, R = pos[:, i], quat_to_mat(quat[:, i])
    yaw = np.unwrap(np.arctan2(R[:, 1, 0], R[:, 0, 0]))
    spin = np.abs(np.gradient(yaw) * fps)                                   # rad/s about the vertical
    flare = np.radians(18) + np.radians(45) * np.clip(spin / 3.0, 0, 1)
    a = np.linspace(0, 2 * np.pi, n_around, endpoint=False)
    s = np.linspace(0, 1, n_down)
    verts = np.empty((len(P), n_around * n_down, 3))
    for k in range(len(P)):
        r = 0.175 + s[:, None] * length * np.sin(flare[k])
        zz = 0.10 - s[:, None] * length * np.cos(flare[k])
        local = np.stack([np.broadcast_to(r * np.cos(a)[None], (n_down, n_around)) - 0.015,
                          np.broadcast_to(r * np.sin(a)[None], (n_down, n_around)),
                          np.broadcast_to(zz, (n_down, n_around))], -1).reshape(-1, 3)
        verts[k] = local @ R[k].T + P[k]
    faces = []
    for d in range(n_down - 1):
        for j in range(n_around):
            v0, v1 = d * n_around + j, d * n_around + (j + 1) % n_around
            faces += [[v0, v1, v1 + n_around], [v0, v1 + n_around, v0 + n_around]]
    np.savez_compressed(out, faces=np.array(faces, int), verts=verts.astype(np.float32), fps=fps,
                        evidence='kinematic test skirt for renderer plumbing; no cloth physics')
    print(f'wrote {out}: {len(verts)} frames, {verts.shape[1]} vertices, {len(faces)} faces')
    return 0


if __name__ == '__main__':
    raise SystemExit(main(sys.argv[1], sys.argv[2]))
