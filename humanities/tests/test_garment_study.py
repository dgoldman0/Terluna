"""The fabric and cut study: sewn ring sectors and their cones, open weaves and weighted hems in the solver, the skirt
builder, the closed-form screen and its product, and one short drape trial."""
from __future__ import annotations

import json
import math

import numpy as np
import pytest
import torch

from humanities.dance.garment import screen as S
from humanities.dance.garment import validate as V
from humanities.dance.garment.cloth import EARTH, MOON, Cloth, Settings, bending_operator
from humanities.dance.garment.materials import FABRICS
from humanities.dance.garment.patterns import boundary_vertices, sector

torch.set_num_threads(1)


@pytest.mark.parametrize('cut', [math.pi, 1.5 * math.pi, 2 * math.pi, 4 * math.pi])
def test_sector_is_a_sewn_band_with_the_cut_area(cut):
    r_in = 0.8 / cut
    pc = sector(r_in, r_in + 0.4, cut, 0.03)
    assert pc.area() == pytest.approx(0.5 * cut * ((r_in + 0.4) ** 2 - r_in ** 2), rel=2e-3)
    assert set(boundary_vertices(pc.F)) == set(pc.rings['inner']) | set(pc.rings['outer'])     # sewn: no seam edge


@pytest.mark.parametrize('cut', [math.pi, 1.5 * math.pi])
def test_sewn_cone_bends_as_the_cone_it_makes(cut):
    """A partial circle laid smooth is a cone with half-angle asin(cut / 2pi); its bending energy is
    (B/2) cut cot^2(alpha) ln(r_out / r_in) for B = 1, which the seam's turned corners must reproduce."""
    r_in = 0.8 / cut
    pc = sector(r_in, r_in + 0.4, cut, 0.015)
    M, Av = bending_operator(pc)
    sa = cut / (2 * math.pi)
    ca = math.sqrt(1 - sa * sa)
    rho = np.hypot(*pc.vertex_uv.T)
    th = np.mod(np.arctan2(pc.vertex_uv[:, 1], pc.vertex_uv[:, 0]), 2 * math.pi) * 2 * math.pi / cut
    X = np.c_[rho * sa * np.cos(th), rho * sa * np.sin(th), -rho * ca]
    C = M @ X
    E = (1 / (2 * Av) * (C * C).sum(1)).sum()
    assert E == pytest.approx(0.5 * cut * (ca / sa) ** 2 * math.log((r_in + 0.4) / r_in), rel=0.02)


def test_full_and_double_circles_lie_flat_without_bending():
    for cut in (2 * math.pi, 4 * math.pi):
        pc = sector(0.8 / cut, 0.8 / cut + 0.3, cut, 0.03)
        M, _ = bending_operator(pc)
        assert np.abs(M @ np.c_[pc.vertex_uv, np.zeros(pc.n_vertices)]).max() < 1e-12


def test_open_weave_falls_faster_by_the_root_of_its_factor():
    m = FABRICS['silk_chiffon']
    net = type(m)(**{**m.__dict__, 'name': 'open', 'normal_factor': 0.5})
    an = math.sqrt(2 * m.areal_density * MOON.gravity / (MOON.air_density * 1.17 * 0.5))
    ts, vn, vt, _ = V.terminal(net, MOON, h=0.025, size=0.1, T=10 * an / MOON.gravity)
    assert abs(-vn[-10:].mean() / an - 1) < 0.01


def test_weighted_hem_adds_mass_the_air_does_not_see():
    from humanities.dance.garment.body import Body, MotionPose
    from humanities.dance.garment.garments import skirt
    body = Body()
    mp = MotionPose(V.TEST_PATTERN, body, t_raised=1.0, t_lower=1.0, t_hold=1.0)
    for cut in (math.pi, 2 * math.pi):
        G = skirt(body, mp(-mp.pre), cut=cut, hem_weight_kg_m=0.05)
        F, UV = G.piece.F, G.piece.UV
        placed = np.linalg.norm(G.X0[F[:, 1]] - G.X0[F[:, 0]], axis=1)
        cut_len = np.linalg.norm(UV[:, 1] - UV[:, 0], axis=1)
        # it starts as the cone it makes, unstretched but for the chord across each edge's curvature (a flat full
        # circle has none)
        assert np.abs(placed / cut_len - 1).max() < (1e-9 if cut == 2 * math.pi else 2e-3)
        assert G.params['smooth_cone_flare_deg'] == pytest.approx(math.degrees(math.asin(cut / (2 * math.pi))))
        assert G.extra_mass.sum() == pytest.approx(0.05 * G.params['hem_circumference_m'])
        mat = FABRICS['silk_chiffon']
        cl = Cloth(G.piece, G.X0, mat, MOON, Settings(), pinned=G.pinned, extra_mass=G.extra_mass, device='cpu')
        assert float(cl.m.sum()) == pytest.approx(mat.areal_density * G.piece.area() + G.extra_mass.sum(), rel=1e-5)


def test_screen_closed_forms():
    g, rho = MOON.gravity, MOON.air_density
    m = 0.03
    vt = S.face_on_fall(m, g, rho)
    assert S.fall_time(3.0, vt, g) == pytest.approx(3.0 / vt + vt / g * math.log(2), rel=1e-6)
    for U in (1.0, 2.0, 3.0):
        a = math.radians(S.trail_angle(m, g, rho, U))
        assert 0.5 * rho * S.C_N * U ** 2 * math.sin(a) ** 2 == pytest.approx(m * g * math.cos(a), rel=1e-6)
    # a lunar fabric falls and trails as one lighter by (g_E / g_M)(rho_M / rho_E) does on Earth
    G = json.loads(S.OUT.read_text())
    row = next(r for r in G['fabrics'] if r['fabric'] == 'wool_plain')
    mE, (gM, rM), (gE, rE) = row['falls_and_trails_as_on_earth_kg_m2'], *(
        (G['worlds'][w]['gravity_m_s2'], G['worlds'][w]['air_kg_m3']) for w in ('moon', 'earth'))
    assert S.face_on_fall(row['areal_density_kg_m2'], gM, rM) == pytest.approx(S.face_on_fall(mE, gE, rE), rel=1e-3)
    assert S.trail_angle(row['areal_density_kg_m2'], gM, rM, 2.0) == pytest.approx(S.trail_angle(mE, gE, rE, 2.0), rel=1e-3)
    # in the time unit sqrt(L / g) the swing holds no gravity: in one air, the Moon's is Earth's slowed exactly
    sm, se = S.swing(m, MOON.gravity, rho), S.swing(m, EARTH.gravity, rho)
    assert sm['largest_swing_back_deg'] == pytest.approx(se['largest_swing_back_deg'], abs=0.11)
    assert sm['first_through_vertical_s'] / se['first_through_vertical_s'] == pytest.approx(
        math.sqrt(EARTH.gravity / MOON.gravity), rel=5e-3)


def test_meridian_extends_the_validation_string():
    r_w = 0.13
    om = 20 * math.pi / 30
    for g in (MOON.gravity, EARTH.gravity):
        for rho, kw in ((0.0, dict(drag=False)), (1.4, dict(drag=True, folds=False))):
            a = S.meridian(0.0261, g, rho, om, r_w=r_w, folds=False)
            b = V.chain_estimate(om, g, 1.4, r_w, S.SKIRT_LENGTH_M, 0.0261, **kw)
            assert a['converged'] and b['converged']
            assert a['flare_deg'] == pytest.approx(b['flare_deg'], abs=0.3)
            assert a['lag_deg'] == pytest.approx(b['lag_deg'], abs=0.3)
    # a half circle never opens past its smooth cone; a weighted hem opens a light skirt further
    for rho in (0.0, 1.4):
        half = S.meridian(0.0261, MOON.gravity, rho, 6.0, cut=math.pi, r_w=r_w, folds=False)
        assert half['converged'] and half['hem_radius_m'] <= r_w + S.SKIRT_LENGTH_M / 2 + 0.002
    # the wind-up index holds no gravity and no rate of turn; a weighted hem lowers it, and the hem it asks for lands it
    w0 = S.wind_up_index(0.0261, 1.4, r_w=r_w)
    w1 = S.wind_up_index(0.0261, 1.4, hem_kg_m=0.1, r_w=r_w)
    assert w1 < w0 / 5
    assert S.wind_up_index(0.0261, 1.4, hem_kg_m=S.hem_for_index(0.0261, 1.4, 0.5, r_w=r_w), r_w=r_w) == pytest.approx(0.5)


def test_screen_product_regenerates():
    assert json.loads(S.OUT.read_text()) == json.loads(json.dumps(S.results()))


def test_tethers_hold_a_heavy_hem_without_touching_slack_cloth():
    """A strip of chiffon hanging from its top edge with a heavy bar sewn along its foot: the tethers keep the foot
    within their length of the top, and change nothing while the cloth is slack."""
    from humanities.dance.garment.patterns import rectangle
    pc = rectangle(0.06, 0.5, 0.02, warp=(0.0, 1.0))
    uv = pc.vertex_uv
    X0 = np.c_[uv[:, 0], np.zeros(pc.n_vertices), 1.0 - uv[:, 1]]
    top = np.array(pc.rings['start'])
    foot = np.array(pc.rings['end'])
    extra = np.zeros(pc.n_vertices)
    extra[foot] = 0.05 / len(foot)                        # 50 g on a 6 cm strip of 0.8 g
    anchor = top[np.abs(uv[:, None, 0] - uv[top][None, :, 0]).argmin(1)]
    length = 1.02 * np.abs(uv[:, 1] - uv[anchor, 1])
    out = {}
    for tethers in (None, (anchor, length)):
        cl = Cloth(pc, X0, FABRICS['silk_chiffon'], EARTH, Settings(self_collision=False), pinned=top,
                   extra_mass=extra, tethers=tethers, device='cpu')
        cl.set_targets(torch.as_tensor(X0[top], dtype=torch.float32))
        cl.damp.fill_(5.0)
        for _ in range(1500):
            cl.substep(1 / 2400)
        out[tethers is None] = float(1.0 - cl.x[foot, 2].mean())
    assert out[False] <= 0.5 * 1.02 + 1e-3                # tethered: within the tether
    assert out[True] > out[False]                          # untethered, one stretch pass lets the bar pull further


def test_cuff_tethers_hold_each_sleeve_to_its_own_cuff():
    """Water sleeves' tethers anchor every vertex to a vertex of its own sleeve's cuff, at 1.05 times its distance
    from the cuff in the pattern, which for a tube is its distance along the sleeve."""
    from humanities.dance.garment.body import Body, MotionPose
    from humanities.dance.garment.garments import cuff_tethers, water_sleeves
    from humanities.dance.garment.trials import TEST_PATTERN
    body = Body()
    mp = MotionPose(TEST_PATTERN, body, t_raised=1.0, t_lower=1.0, t_hold=1.0)
    G = water_sleeves(body, mp(-mp.pre), h=0.05)
    anchor, length = cuff_tethers(G)
    half = G.piece.n_vertices // 2
    uv = G.piece.vertex_uv
    for side, lo in (('L', 0), ('R', half)):
        ring = set(G.piece.rings[f'{side}_start'])
        idx = np.arange(lo, lo + half)
        assert all(a in ring for a in anchor[idx])
        assert np.allclose(length[idx], 1.05 * uv[idx, 1], atol=1e-9)



def test_tapered_piece_pushes_cloth_to_its_radius_there():
    """A dress form's piece tapers from ra to rb along its axis: a vertex that moved inside the cone in this substep
    goes back out to the radius at its place on the axis (plus the body offset), where a straight capsule of radius
    ra would have left it, and a vertex outside the cone stays put."""
    from humanities.dance.garment.patterns import rectangle
    pc = rectangle(0.1, 0.1, 0.05)
    X0 = np.c_[pc.vertex_uv, np.zeros(pc.n_vertices)] + 5.0
    out = {}
    for rb in (0.1, 0.2):
        c = Cloth(pc, X0, FABRICS['silk_chiffon'], MOON, Settings(self_collision=False), device='cpu')
        c.set_body([0.1], [rb])
        z = torch.zeros((1, 3))
        c.set_frame(z, z, z + torch.tensor([[0.0, 0.0, 1.0]]), z + torch.tensor([[0.0, 0.0, 1.0]]))
        c.frac0.fill_(1.0)
        c.frac1.fill_(1.0)
        c.dt_now = 1.0 / 2400.0
        x = c.x.clone()
        x[0] = torch.tensor([0.12, 0.0, 0.5])            # inside the cone (0.15 there), outside the capsule (0.10)
        x[1] = torch.tensor([0.0, 0.16, 0.5])            # outside both
        disp = torch.zeros_like(x)
        disp[0] = torch.tensor([-0.05, 0.0, 0.0])        # it arrived moving inward in this substep
        dx, _ = c.body_contact(x, disp)
        out[rb] = (x + dx)[:2].numpy()
    off = Settings().body_offset
    assert np.allclose(out[0.1][0], [0.12, 0.0, 0.5], atol=1e-7)          # straight capsule: already outside
    assert np.allclose(out[0.2][0], [0.15 + off, 0.0, 0.5], atol=1e-6)    # cone: out to its radius at mid-length
    assert np.allclose(out[0.2][1], [0.0, 0.16, 0.5], atol=1e-7)


def test_skirt_band_follows_a_dress_forms_waist():
    """With a dress form, a skirt's waistband lies on the form's waist outline plus its clearance, its points spaced
    by arc length as the pattern's inner edge is, and the pattern's inner edge is as long as the band."""
    from humanities.dance.garment.body import Body, MotionPose
    from humanities.dance.garment.garments import BAND_CLEARANCE, WAIST_ABOVE_OMPHALION, skirt
    from humanities.dance.garment.trials import TEST_PATTERN
    body = Body()
    ab = {s['name']: s for s in body.spec['segments']}['abdomen']
    axis_x = ab['origin'][0] + ab['capsules'][0]['a'][0]
    z_w = ab['origin'][2] + WAIST_ABOVE_OMPHALION
    th = (np.arange(72) + 0.5) * 2 * np.pi / 72
    r = 0.13 * 0.09 / np.sqrt((0.09 * np.cos(th - np.pi / 2)) ** 2 + (0.13 * np.sin(th - np.pi / 2)) ** 2)
    body.waist = dict(axis_x=axis_x, z=z_w, theta=th.tolist(), r=r.tolist())     # 26 cm wide, 18 cm deep
    mp = MotionPose(TEST_PATTERN, body, t_raised=1.0, t_lower=1.0, t_hold=1.0)
    cut = np.pi
    G = skirt(body, mp(-mp.pre), cut=cut, length=0.3, h=0.03)
    P = G.attachment.p_rest - np.array([axis_x, 0.0, z_w])
    ang = np.mod(np.arctan2(P[:, 1], P[:, 0]), 2 * np.pi)
    want = np.interp(ang, th, r, period=2 * np.pi) + BAND_CLEARANCE
    assert np.allclose(np.hypot(P[:, 0], P[:, 1]), want, atol=1e-4)
    order = np.argsort(np.mod(ang - np.pi, 2 * np.pi))
    Q = np.vstack([P[order], P[order[:1]]])[:, :2]
    gaps = np.linalg.norm(np.diff(Q, axis=0), axis=1)
    assert gaps.std() / gaps.mean() < 0.02
    assert G.params['waist_circumference_m'] == pytest.approx(gaps.sum(), rel=6e-3)     # chords on a 3 cm grid
    assert G.params['waist_circumference_m'] == pytest.approx(cut * G.params['pattern_inner_radius_m'], rel=1e-9)


def test_dress_form_product_matches_the_body():
    """The dress form that dressed renders collide with: its schema, its pieces on body.json's segments, limbs never
    slimmer than body.json's capsules, and the body.json it was measured against."""
    import hashlib
    from pathlib import Path
    from humanities.dance.motion_io import body_spec
    root = Path(__file__).resolve().parents[2]
    form = json.loads((root / 'visualization' / 'dance' / 'dress_form.json').read_text())
    assert form['schema'] == 'terluna.visualization.dance-dress-form/1'
    body_json = root / 'humanities' / 'dance' / 'body' / 'body.json'
    assert form['producer']['body_json_sha256'] == hashlib.sha256(body_json.read_bytes()).hexdigest()
    spec = {s['name']: s for s in body_spec()['segments']}
    for p in form['pieces']:
        assert p['segment'] in spec
        if p['segment'] not in ('pelvis', 'abdomen', 'thorax'):
            r = min(c['radius'] for c in spec[p['segment']]['capsules'])
            assert min(p['ra'], p['rb']) >= r - 1e-9
    w = form['waist']
    assert len(w['theta']) == len(w['r']) and 0.5 < w['perimeter_m'] < 0.9
