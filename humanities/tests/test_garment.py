"""Fast checks of the garment solver: the bending operator, the air forces, layer shielding, the bending length's
gravity scaling, the steady-turn estimates and the garment builders. The heavy validation runs live in
humanities/dance/garment/validate.py."""
from __future__ import annotations

import math

import numpy as np
import pytest
import torch

from humanities.dance.garment.cloth import EARTH, MOON, Cloth, Settings, bending_operator
from humanities.dance.garment.materials import FABRICS, bending_length
from humanities.dance.garment.patterns import annulus, rectangle
from humanities.dance.garment import validate as V

torch.set_num_threads(3)


def _cylinder(uv, R, psi):
    a = uv[:, 0] * math.cos(psi) + uv[:, 1] * math.sin(psi)
    q = -uv[:, 0] * math.sin(psi) + uv[:, 1] * math.cos(psi)
    return np.c_[R * np.sin(q / R), a, R * (1 - np.cos(q / R))]


@pytest.mark.parametrize('piece,tol', [(rectangle(0.6, 0.6, 0.02), 0.07), (annulus(0.11, 0.5, 0.02), 0.03)])
def test_bending_operator_flat_and_cylinder(piece, tol):
    M, Av = bending_operator(piece)
    th = 0.6
    Rz = np.array([[math.cos(th), -math.sin(th), 0], [math.sin(th), math.cos(th), 0], [0, 0, 1.0]])
    X = np.c_[piece.vertex_uv, np.zeros(piece.n_vertices)] @ Rz.T + [0.2, -0.1, 1.0]
    assert np.abs(M @ X).max() < 1e-10                     # every isometric flat placement is free of bending
    A = piece.area()
    for psi in (0.0, math.pi / 8, math.pi / 4, math.pi / 2):
        C = M @ _cylinder(piece.vertex_uv, 0.15, psi)
        E = (1 / (2 * Av) * (C * C).sum(1)).sum()
        assert abs(E / (0.5 * A / 0.15 ** 2) - 1) < tol     # (B/2) integral of k^2 dA with B = 1


@pytest.mark.parametrize('world', [MOON, EARTH])
def test_face_on_terminal_speed(world):
    m = FABRICS['silk_chiffon']
    an = math.sqrt(2 * m.areal_density * world.gravity / (world.air_density * 1.17))
    ts, vn, vt, _ = V.terminal(m, world, h=0.025, size=0.1, T=10 * an / world.gravity)
    assert abs(-vn[-10:].mean() / an - 1) < 0.01


def test_two_layers_share_one_face():
    m = FABRICS['silk_chiffon']
    an1 = math.sqrt(2 * m.areal_density * MOON.gravity / (MOON.air_density * 1.17))
    ts, vn, vt, cl = V.terminal(m, MOON, layers=2, h=0.025, size=0.1, T=12 * math.sqrt(2) * an1 / MOON.gravity)
    assert int(cl.stacked.sum()) == cl.V
    assert abs(-vn[-10:].mean() / (math.sqrt(2) * an1) - 1) < 0.03


def test_rigid_glide_normal_and_skin_friction():
    m = FABRICS['silk_chiffon']
    rigid = type(m)(**{**m.__dict__, 'bending_rigidity': 1e-3})
    beta = math.radians(30)
    ts, vn, vt, _ = V.terminal(rigid, EARTH, beta_deg=30.0, h=0.025, size=0.1, dt=1 / 1200, T=7.0)
    an_n = math.sqrt(2 * m.areal_density * EARTH.gravity * math.cos(beta) / (EARTH.air_density * 1.17))
    an_t = math.sqrt(2 * m.areal_density * EARTH.gravity * math.sin(beta) / (EARTH.air_density * 0.01))
    assert abs(abs(vn[-10:].mean()) / an_n - 1) < 0.03
    assert abs(vt[-10:].mean() / an_t - 1) < 0.03


def test_bending_length_scales_with_gravity():
    m = FABRICS['silk_georgette']
    rM = V.cantilever(m, MOON, n_el=6, settle_periods=4)
    rE = V.cantilever(m, EARTH, n_el=6, settle_periods=4)
    ratio = (EARTH.gravity / MOON.gravity) ** (1 / 3)
    assert abs(bending_length(m, MOON.gravity) / bending_length(m, EARTH.gravity) - ratio) < 1e-9
    assert abs(rM['c_sim'] / rE['c_sim'] / ratio - 1) < 0.02
    # at six elements over the overhang the clamp's discretisation shortens c by about a tenth (O(h))
    for r in (rM, rE):
        assert -0.15 < r['rel_error_c'] < 0.0
        assert abs(r['theta_elastica_deg'] - 41.5) < 0.2       # Peirce's 41.5 deg at l = 2c


def test_turn_estimates():
    r_w, L, m_a = 0.114, 0.85, 0.0261
    assert V.cone_estimate(0.0, 1.6242, r_w, L) == 0.0
    still = V.chain_estimate(1e-6, 1.6242, 1.4, r_w, L, m_a, drag=False)
    assert still['flare_deg'] < 0.1
    for om in (1.0, 2.0):
        nd = V.chain_estimate(om, 1.6242, 1.4, r_w, L, m_a, drag=False)
        wd = V.chain_estimate(om, 1.6242, 1.4, r_w, L, m_a, drag=True, folds=False)
        assert wd['converged']
        cone = math.degrees(V.cone_estimate(om, 1.6242, r_w, L))
        assert abs(nd['flare_deg'] - cone) < 6.0                 # a string and a rigid cone nearly agree
        assert abs(nd['lag_deg']) < 1e-6 and wd['lag_deg'] > 0   # drag makes the hem lag the waist
    # lunar gravity flares the same skirt far more at the same rate
    assert V.cone_estimate(1.5, 1.6242, r_w, L) > 2 * V.cone_estimate(1.5, 9.80665, r_w, L)


def test_garment_builders():
    from humanities.dance.garment.body import Body, MotionPose
    from humanities.dance.garment.garments import circle_skirt, water_sleeves
    body = Body()
    pose = MotionPose(V.TEST_PATTERN, body)
    p0 = pose(-pose.pre)
    sk = circle_skirt(body, p0, h=0.03)
    r_w = sk.params['waist_radius_m']
    tgt = sk.attachment.targets(*p0)
    assert np.allclose(np.linalg.norm(tgt[:, :2] - tgt[:, :2].mean(0), axis=1), r_w, atol=1e-6)
    assert np.allclose(sk.X0[sk.pinned], tgt)
    cl = Cloth(sk.piece, sk.X0, FABRICS['silk_chiffon'], MOON, Settings(), pinned=sk.pinned, device='cpu')
    st = cl.strain_stats()
    assert abs(st['max_stretch']) < 1e-5                         # the flat disc is the pattern itself
    sl = water_sleeves(body, p0, h=0.03)
    assert len(sl.pinned) == 2 * len(sl.piece.rings['L_start'])
    cl = Cloth(sl.piece, sl.X0, FABRICS['silk_georgette'], MOON, Settings(), pinned=sl.pinned, device='cpu')
    assert abs(cl.strain_stats()['max_stretch']) < 1e-3         # the straight tube is the rolled pattern
