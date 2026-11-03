"""Cloths cut from several fabrics: each triangle carries its own fabric's density, stiffness and air factor and each
vertex its own fabric's bending rigidity and friction. Checked against each fabric's mass, fall and frictional slide,
against each piece simulated alone, through layers of two fabrics meeting in the self-contact, and against a
single-fabric cloth when every triangle is given the same fabric."""
from __future__ import annotations

import dataclasses
import math

import numpy as np
import pytest
import torch

from humanities.dance.garment.cloth import EARTH, MOON, Cloth, Settings, World, merged_tri_material
from humanities.dance.garment.materials import FABRICS
from humanities.dance.garment.patterns import rectangle

torch.set_num_threads(1)

CHIFFON, GEORGETTE, SATIN = FABRICS['silk_chiffon'], FABRICS['silk_georgette'], FABRICS['silk_satin_heavy']
DT = 1 / 2400


def _merge(pieces):
    from humanities.dance.garment.garments import merge
    return merge(pieces)


def _piece_of_vertex(pieces):
    return np.repeat(np.arange(len(pieces)), [p.n_vertices for _, p in pieces])


def test_each_piece_carries_its_own_fabrics_mass():
    """A chiffon panel and a heavy satin panel of different sizes merged into one cloth: each piece's vertex masses sum
    to its own area times its own areal density, and every triangle carries its own fabric's mass."""
    pieces = [('veil', rectangle(0.3, 0.2, 0.025)), ('under', rectangle(0.2, 0.25, 0.025))]
    pc = _merge(pieces)
    cl = Cloth(pc, np.c_[pc.vertex_uv, np.zeros(pc.n_vertices)], [CHIFFON, SATIN], MOON, Settings(),
               tri_material=merged_tri_material(pieces, [0, 1]), device='cpu')
    a, b = pieces[0][1], pieces[1][1]
    m = cl.m.double().numpy()
    assert m[:a.n_vertices].sum() == pytest.approx(CHIFFON.areal_density * a.area(), rel=1e-6)
    assert m[a.n_vertices:].sum() == pytest.approx(SATIN.areal_density * b.area(), rel=1e-6)
    want = np.r_[np.full(len(a.F), CHIFFON.areal_density), np.full(len(b.F), SATIN.areal_density)]
    assert np.allclose(cl.mtri.double().numpy() / cl.area0.double().numpy(), want, rtol=1e-6)
    assert cl.material is CHIFFON and cl.materials == (CHIFFON, SATIN)


def _cantilevers(fabrics, slots, n_sub=400):
    """Strips 4 cm wide and 12 cm long, clamped horizontal by their first two rows of vertices and 10 cm apart, swing
    down for n_sub substeps on Earth. One fabric is given as a single Material, several per triangle."""
    pieces = [(f's{k}', rectangle(0.04, 0.12, 0.01, warp=(0.0, 1.0))) for k in range(len(fabrics))]
    pc = _merge(pieces)
    kv = _piece_of_vertex(pieces)
    uv = pc.vertex_uv
    X0 = np.c_[uv[:, 1], uv[:, 0] + 0.1 * np.asarray(slots)[kv], np.full(pc.n_vertices, 1.0)]
    clamp = np.nonzero(uv[:, 1] < 0.015)[0]
    if len(fabrics) == 1:
        cl = Cloth(pc, X0, fabrics[0], EARTH, Settings(self_collision=False), pinned=clamp, device='cpu')
    else:
        cl = Cloth(pc, X0, list(fabrics), EARTH, Settings(self_collision=False), pinned=clamp,
                   tri_material=merged_tri_material(pieces, range(len(fabrics))), device='cpu')
    cl.set_targets(torch.as_tensor(X0[clamp], dtype=torch.float32))
    for _ in range(n_sub):
        cl.substep(DT)
    return cl.x.double().numpy(), kv


def test_pieces_keep_their_own_physics_in_one_cloth():
    """Strips of chiffon, georgette and heavy satin swinging down from a clamp in one cloth each move as the same strip
    does alone in a single-fabric cloth (within 1e-6 m), which needs every strip's own mass, air factor, stretch,
    shear and bending. Giving the georgette strip chiffon's bending rigidity alone moves it by about a millimetre, so
    the comparison would catch a fabric read from the wrong strip."""
    fabrics = [CHIFFON, GEORGETTE, SATIN]
    x, kv = _cantilevers(fabrics, [0, 1, 2])
    for k, f in enumerate(fabrics):
        alone, _ = _cantilevers([f], [k])
        assert np.abs(x[kv == k] - alone).max() < 1e-6
    swapped, _ = _cantilevers([dataclasses.replace(GEORGETTE, bending_rigidity=CHIFFON.bending_rigidity)], [1])
    assert np.abs(x[kv == 1] - swapped).max() > 3e-4
    tips = [x[kv == k][-5:, 2].mean() for k in range(3)]
    assert tips[2] < tips[1] - 0.01 and tips[1] < tips[0] - 0.003    # air holds back the lighter strips more


def test_sheets_fall_at_their_own_terminal_speeds():
    """Three 10 cm sheets in one cloth fall face-on through Earth's air side by side: chiffon, heavy satin, and an open
    weave of chiffon taking half a closed sheet's normal pressure. Each reaches its own sqrt(2 m_a g / (rho C_N f))."""
    open_weave = dataclasses.replace(CHIFFON, name='open', normal_factor=0.5)
    fabrics = [CHIFFON, SATIN, open_weave]
    pieces = [(f'p{k}', rectangle(0.1, 0.1, 0.025)) for k in range(3)]
    pc = _merge(pieces)
    kv = _piece_of_vertex(pieces)
    X0 = np.c_[pc.vertex_uv[:, 0] + 0.4 * kv, pc.vertex_uv[:, 1], np.full(pc.n_vertices, 5.0)]
    s = Settings(self_collision=False)
    cl = Cloth(pc, X0, fabrics, EARTH, s, tri_material=merged_tri_material(pieces, [0, 1, 2]), device='cpu')
    speed = [math.sqrt(2 * f.areal_density * EARTH.gravity / (EARTH.air_density * s.c_normal * f.normal_factor))
             for f in fabrics]
    n = int(5 * max(speed) / EARTH.gravity / s.dt)            # five time constants of the slowest to settle
    tail, vz = 50, np.zeros(3)
    for k in range(n):
        cl.substep(s.dt)
        if k >= n - tail:
            vz += np.bincount(kv, weights=cl.v[:, 2].double().numpy(), minlength=3) / np.bincount(kv) / tail
    for v, want in zip(-vz, speed):
        assert v == pytest.approx(want, rel=0.01)


def test_vertex_bending_and_friction_come_from_the_vertexs_fabric():
    """A square cut half from chiffon and half from heavy satin: a vertex inside either half takes that fabric's
    bending rigidity and friction, and a vertex on the seam their mean weighted by its triangles' areas, which on this
    even grid is their number. The bending energy of the square rolled into a cylinder uses each vertex's own value."""
    pc = rectangle(0.2, 0.2, 0.025)
    tri = (pc.UV[:, :, 0].mean(1) > 0.1).astype(np.int64)         # triangles right of u = 0.1 are satin
    uv = pc.vertex_uv
    R = 0.15
    X = np.c_[R * np.sin(uv[:, 0] / R), uv[:, 1], R * (1 - np.cos(uv[:, 0] / R))]
    cl = Cloth(pc, X, [CHIFFON, SATIN], MOON, Settings(), tri_material=tri, device='cpu')
    n_all = np.bincount(pc.F.ravel(), minlength=pc.n_vertices)
    n_sat = np.bincount(pc.F[tri == 1].ravel(), minlength=pc.n_vertices)
    seam = (n_sat > 0) & (n_sat < n_all)
    assert np.allclose(uv[seam, 0], 0.1) and seam.sum() == 9
    for name in ('bending_rigidity', 'friction'):
        a, b = getattr(CHIFFON, name), getattr(SATIN, name)
        want = ((n_all - n_sat) * a + n_sat * b) / n_all
        got = getattr(cl, name).double().numpy()
        assert np.allclose(got, want, rtol=1e-6, atol=0)
        assert np.all(got[uv[:, 0] < 0.09] == np.float32(a)) and np.all(got[uv[:, 0] > 0.11] == np.float32(b))
    C = (cl.M @ cl.x).double().numpy()
    Av = cl.Av.double().numpy()
    B = ((n_all - n_sat) * CHIFFON.bending_rigidity + n_sat * SATIN.bending_rigidity) / n_all
    assert cl.energy()['bending'] == pytest.approx((B / (2 * Av) * (C * C).sum(1)).sum(), rel=1e-5)


def test_each_sheet_slides_to_its_own_frictional_stop():
    """A chiffon and a heavy satin sheet in one cloth slide on the floor at 1 m/s with no air: each stops after
    v^2 / (2 mu g), its own Coulomb friction (0.3 and 0.25) acting on its own vertices."""
    vacuum = World('vacuum', EARTH.gravity, 0.0)
    pieces = [(f'p{k}', rectangle(0.1, 0.1, 0.025)) for k in range(2)]
    pc = _merge(pieces)
    kv = _piece_of_vertex(pieces)
    s = Settings(self_collision=False)
    X0 = np.c_[pc.vertex_uv[:, 0], pc.vertex_uv[:, 1] + 0.3 * kv, np.full(pc.n_vertices, s.floor_offset)]
    cl = Cloth(pc, X0, [CHIFFON, SATIN], vacuum, s, tri_material=merged_tri_material(pieces, [0, 1]), device='cpu')
    cl.v[:, 0] = 1.0
    for _ in range(int(0.45 / s.dt)):                          # both stop within 0.41 s
        cl.substep(s.dt)
    x, v = cl.x.double().numpy(), cl.v.double().numpy()
    for k, f in enumerate((CHIFFON, SATIN)):
        assert (x[kv == k, 0] - X0[kv == k, 0]).mean() == pytest.approx(1.0 / (2 * f.friction * EARTH.gravity), rel=0.01)
        assert np.abs(v[kv == k]).max() < 1e-4


def test_chiffon_over_satin_meet_through_the_self_contact():
    """A chiffon sheet laid over a heavy satin sheet, both falling onto a horizontal capsule through Earth's air: the
    satin falls faster and drapes first, and the chiffon settles onto it and stays one contact distance outside it
    over the crest. The stacked vertices take shielded air."""
    pieces = [('under', rectangle(0.2, 0.2, 0.02)), ('veil', rectangle(0.2, 0.2, 0.02))]
    pc = _merge(pieces)
    kv = _piece_of_vertex(pieces)
    d = 0.75 * pc.edge_length()                                  # the self-contact distance
    X0 = np.c_[pc.vertex_uv - 0.1, np.where(kv == 0, 0.09, 0.09 + 1.2 * d)]
    cl = Cloth(pc, X0, [SATIN, CHIFFON], EARTH, Settings(), tri_material=merged_tri_material(pieces, [0, 1]),
               device='cpu')
    assert cl.d_self == pytest.approx(d)
    cl.set_body([0.05])
    a, b = torch.tensor([[0.0, -0.3, 0.03]]), torch.tensor([[0.0, 0.3, 0.03]])
    cl.set_frame(a, a, b, b)
    for _ in range(600):
        cl.substep(DT)
    x = cl.x.double().numpy()
    under, veil = x[kv == 0], x[kv == 1]
    gap = np.linalg.norm(veil[:, None] - under[None], axis=2).min(1)
    assert gap.min() > 0.9 * d
    crest = (under[np.abs(under[:, 0]) < 0.011, 2].mean(), veil[np.abs(veil[:, 0]) < 0.011, 2].mean())
    assert 0.9 * d < crest[1] - crest[0] < 1.2 * d
    assert under[:, 2].min() < veil[:, 2].min()                    # the satin hangs lower beneath the chiffon
    assert int(cl.stacked.sum()) > 0.8 * cl.V


def _drape(material, n_sub=300, **kw):
    """A 0.3 m square of 3 cm cells dropped at 1 m/s onto a horizontal capsule drapes onto the floor (Earth, air,
    body and floor friction and self-contact all acting)."""
    pc = rectangle(0.3, 0.3, 0.03)
    X0 = np.c_[pc.vertex_uv - 0.15, np.full(pc.n_vertices, 0.085)]
    cl = Cloth(pc, X0, material, EARTH, Settings(), device='cpu', **kw)
    cl.set_body([0.05])
    a, b = torch.tensor([[0.0, -0.3, 0.03]]), torch.tensor([[0.0, 0.3, 0.03]])
    cl.set_frame(a, a, b, b)
    cl.v[:] = torch.tensor([0.2, 0.0, -1.0])
    for _ in range(n_sub):
        cl.substep(DT)
    return cl


def test_one_fabric_per_triangle_matches_the_single_material():
    """A single Material keeps the solver's scalar coefficients; the same fabric given per triangle, as one entry or as
    two copies alternating over the triangles, runs on per-element tensors and must agree with it to float32
    rounding: within 1e-6 m after 300 substeps of a drape that reaches the capsule, the floor and the self-contact."""
    one = _drape(SATIN)
    x = one.x.double().numpy()
    assert x[:, 2].min() < 0.0025 and one.stats['pairs_max'] > 0      # the floor and the self-contact were reached
    T = len(one.F)
    for mats, tri in (([SATIN], np.zeros(T, np.int64)), ([SATIN, SATIN], np.arange(T) % 2)):
        per = _drape(mats, tri_material=tri)
        assert isinstance(per.friction, torch.Tensor) and isinstance(one.friction, float)
        assert np.abs(per.x.double().numpy() - x).max() < 1e-6


def test_material_assignment_is_checked():
    pc = rectangle(0.1, 0.1, 0.05)
    X0 = np.c_[pc.vertex_uv, np.zeros(pc.n_vertices)]
    T = len(pc.F)
    with pytest.raises(ValueError):
        Cloth(pc, X0, [CHIFFON, SATIN], MOON, device='cpu')                                  # no assignment
    with pytest.raises(ValueError):
        Cloth(pc, X0, [CHIFFON, SATIN], MOON, tri_material=np.zeros(T - 1, int), device='cpu')
    with pytest.raises(ValueError):
        Cloth(pc, X0, [CHIFFON, SATIN], MOON, tri_material=np.full(T, 2), device='cpu')
    with pytest.raises(ValueError):
        merged_tri_material([('a', pc), ('b', pc)], [0])


def test_outfit_pieces_start_on_their_pins_and_near_their_cut():
    """The gown's sleeves, wrist flounces, gathered shoulder panels and sash start with their pinned edges on their
    attachment targets, pins in vertex order, and their cloth at its cut lengths (the gathered panels within 5%,
    as folded accordions; the full sleeves start straight between their rings, gathered, never stretched)."""
    from humanities.dance.garment.body import Body, MotionPose
    from humanities.dance.garment.garments import arm_sleeves, butterfly_shoulders, cuff_flounces, sash
    from humanities.dance.garment.trials import TEST_PATTERN
    body = Body()
    mp = MotionPose(TEST_PATTERN, body, t_raised=1.0, t_lower=1.0, t_hold=1.0)
    p0 = mp(-mp.pre)
    for G, low in ((arm_sleeves(body, p0, h=0.03), 0.0), (cuff_flounces(body, p0, h=0.03), 0.99),
                   (butterfly_shoulders(body, p0), 0.95), (sash(body, p0, h=0.03), 0.99)):
        assert np.allclose(G.attachment.targets(*p0), G.X0[G.pinned], atol=1e-9)
        assert (np.diff(G.pinned) > 0).all()
        F, X, UV = G.piece.F, G.X0, G.piece.UV
        r = np.concatenate([np.linalg.norm(X[F[:, i]] - X[F[:, j]], axis=1)
                            / np.linalg.norm(UV[:, i] - UV[:, j], axis=1) for i, j in ((0, 1), (1, 2), (2, 0))])
        assert r.max() < 1.05, G.name
        assert r.min() > low, G.name


@pytest.mark.parametrize('slant, top, ring_r', [(0.0, 0.44, 0.055), (0.12, 0.44, 0.055), (0.0, 0.56, 0.058)])
def test_sleeve_tiers_start_on_their_wrist_rings_with_bound_hems(slant, top, ring_r):
    """A pair of sleeve tiers starts with its tops pinned on their wrist rings (pins in vertex order) and its cloth at
    its cut lengths below the gathered top; its bound hem is the last strip of triangles, every one of them touching
    the hem; a slanted hem runs from (1 - slant) to (1 + slant) times the tier's length."""
    from humanities.dance.garment.body import Body, MotionPose
    from humanities.dance.garment.garments import sleeve_tier
    from humanities.dance.garment.trials import TEST_PATTERN
    body = Body()
    mp = MotionPose(TEST_PATTERN, body, t_raised=1.0, t_lower=1.0, t_hold=1.0)
    p0 = mp(-mp.pre)
    G = sleeve_tier(body, p0, length=0.5, top=top, ring_r=ring_r, flare=0.2, slant=slant, h=0.03)
    assert np.isfinite(G.X0).all()                 # a full tier on a small ring once could not widen in time
    assert np.allclose(G.attachment.targets(*p0), G.X0[G.pinned], atol=1e-9)
    assert (np.diff(G.pinned) > 0).all()
    F, X, UV = G.piece.F, G.X0, G.piece.UV
    r = np.concatenate([np.linalg.norm(X[F[:, i]] - X[F[:, j]], axis=1)
                        / np.linalg.norm(UV[:, i] - UV[:, j], axis=1) for i, j in ((0, 1), (1, 2), (2, 0))])
    assert r.max() < 1.05
    assert r.min() > 2 * np.pi * ring_r / top - 0.03                    # gathered into the ring, never stretched
    hem = set(G.piece.rings['p0_outer']) | set(G.piece.rings['p1_outer'])     # the left and right tiers' hems
    bt = np.asarray(G.params['binding_tris'])
    assert len(bt) > 0 and all(set(F[t]) & hem for t in bt)
    uv = G.piece.vertex_uv
    rho = np.hypot(uv[:, 0], uv[:, 1])
    for side in ('p0', 'p1'):                                             # each side in its own pattern
        d = rho[G.piece.rings[f'{side}_outer']] - rho[G.piece.rings[f'{side}_inner']].mean()
        assert d.min() == pytest.approx(0.5 * (1 - slant), abs=0.02)
        assert d.max() == pytest.approx(0.5 * (1 + slant), abs=0.02)


def test_layers_in_one_contact_group_do_not_repel():
    """Layers gathered into one seam may lie closer than the self-contact distance. The chiffon and satin of the
    layered drape above, put in one contact group, settle on the capsule together, the chiffon within a fifth of a
    contact distance of the satin over the crest; with the group numbers 0 they keep one contact distance apart."""
    pieces = [('under', rectangle(0.2, 0.2, 0.02)), ('veil', rectangle(0.2, 0.2, 0.02))]
    pc = _merge(pieces)
    kv = _piece_of_vertex(pieces)
    d = 0.75 * pc.edge_length()
    X0 = np.c_[pc.vertex_uv - 0.1, np.where(kv == 0, 0.09, 0.09 + 1.2 * d)]
    gap = {}
    for label, skip in (('grouped', np.ones(pc.n_vertices, int)), ('apart', np.zeros(pc.n_vertices, int))):
        cl = Cloth(pc, X0, [SATIN, CHIFFON], EARTH, Settings(), tri_material=merged_tri_material(pieces, [0, 1]),
                   device='cpu', contact_skip=skip)
        cl.set_body([0.05])
        a, b = torch.tensor([[0.0, -0.3, 0.03]]), torch.tensor([[0.0, 0.3, 0.03]])
        cl.set_frame(a, a, b, b)
        for _ in range(600):
            cl.substep(DT)
        x = cl.x.double().numpy()
        under, veil = x[kv == 0], x[kv == 1]
        gap[label] = veil[np.abs(veil[:, 0]) < 0.011, 2].mean() - under[np.abs(under[:, 0]) < 0.011, 2].mean()
    assert abs(gap['grouped']) < 0.2 * d
    assert 0.9 * d < gap['apart'] < 1.2 * d
