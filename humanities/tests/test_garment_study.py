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


@pytest.mark.parametrize('cut', [math.pi, 1.5 * math.pi, 2 * math.pi])
def test_waist_tethers_reach_the_near_waist_for_every_cut(cut):
    """Each skirt vertex is tethered to the waistband point nearest it in the sewn pattern, across the seam where that
    is shorter, at 1.05 times their pattern distance: on the skirt's starting cone (isometric to its pattern) no
    tether is shorter than the distance it starts at. Cuts between a half and a full circle once tethered cloth to
    the far side of the waist."""
    from humanities.dance.garment.body import Body, MotionPose
    from humanities.dance.garment.garments import skirt, waist_tethers
    from humanities.dance.garment.trials import TEST_PATTERN
    body = Body()
    mp = MotionPose(TEST_PATTERN, body, t_raised=1.0, t_lower=1.0, t_hold=1.0)
    G = skirt(body, mp(-mp.pre), cut=cut, length=0.5, h=0.04)
    anchor, length = waist_tethers(G)
    d0 = np.linalg.norm(G.X0 - G.X0[anchor], axis=1)
    assert (length >= d0 - 1e-9).all()
    # the anchor sits on the vertex's own line down the pattern, within a mesh spacing of it
    rho = np.hypot(*G.piece.vertex_uv.T)
    assert (length <= 1.05 * (rho - rho.min() + 0.04) + 1e-9).all()


def test_sector_hem_cuts_a_slanted_hem():
    """A sector given a hem function ends at that hem: its outer ring lies at r_in plus the hem's distance there, its
    area is the integral of the slanted band, and it stays one sewn band with only its top and hem free."""
    r_in, L, s, cut = 3.0, 0.5, 0.2, 0.15
    hem = lambda u: L * (1 + s * np.cos(2 * np.pi * np.asarray(u)))
    pc = sector(r_in, r_in + L * (1 + s), cut, 0.02, hem=hem)
    uv = pc.vertex_uv
    out = np.array(pc.rings['outer'])
    u = np.mod(np.arctan2(uv[out, 1], uv[out, 0]), 2 * np.pi) / cut
    assert np.allclose(np.hypot(*uv[out].T), r_in + hem(u), atol=1e-12)
    psi = np.linspace(0.0, cut, 20001)
    f = 0.5 * ((r_in + hem(psi / cut)) ** 2 - r_in ** 2)
    want = float(((f[1:] + f[:-1]) * np.diff(psi)).sum() / 2)
    assert pc.area() == pytest.approx(want, rel=2e-3)
    assert set(boundary_vertices(pc.F)) == set(pc.rings['inner']) | set(pc.rings['outer'])


def _oval_waist(body):
    """An oval waist outline 26 cm wide and 18 cm deep about the abdomen capsule's axis, as on a dress form."""
    from humanities.dance.garment.garments import WAIST_ABOVE_OMPHALION
    ab = {s['name']: s for s in body.spec['segments']}['abdomen']
    axis_x = ab['origin'][0] + ab['capsules'][0]['a'][0]
    z_w = ab['origin'][2] + WAIST_ABOVE_OMPHALION
    th = (np.arange(72) + 0.5) * 2 * np.pi / 72
    r = 0.13 * 0.09 / np.sqrt((0.09 * np.cos(th - np.pi / 2)) ** 2 + (0.13 * np.sin(th - np.pi / 2)) ** 2)
    return dict(axis_x=axis_x, z=z_w, theta=th.tolist(), r=r.tolist()), axis_x, z_w


@pytest.mark.parametrize('cut', [math.pi, 2 * math.pi])
def test_skirt_on_a_waist_outline_starts_at_its_pattern_lengths_and_clear_of_the_body(cut):
    """On a waist outline that is not a circle, a skirt starts out along the outline's normal and turns into its
    pattern's own cone about the outline's centre, so the starting cloth keeps its cut lengths everywhere but in that
    blend (laid out from the capsule axis instead, half its triangles started more than 15% off their cut), and no
    free starting point lies inside the body."""
    from humanities.dance.garment.body import Body, MotionPose
    from humanities.dance.garment.garments import skirt
    from humanities.dance.garment.trials import TEST_PATTERN
    body = Body()
    body.waist, _, _ = _oval_waist(body)
    mp = MotionPose(TEST_PATTERN, body, t_raised=1.0, t_lower=1.0, t_hold=1.0)
    p0 = mp(-mp.pre)
    G = skirt(body, p0, cut=cut, length=0.6, h=0.02)
    F, X, UV = G.piece.F, G.X0, G.piece.UV
    e = np.stack([np.linalg.norm(X[F[:, i]] - X[F[:, j]], axis=1) / np.linalg.norm(UV[:, i] - UV[:, j], axis=1)
                  for i, j in ((0, 1), (1, 2), (2, 0))], 1)
    strain = np.abs(e - 1).max(1)
    assert np.median(strain) < 0.01
    assert np.percentile(strain, 95) < 0.15
    free = np.setdiff1d(np.arange(G.piece.n_vertices), G.pinned)
    A, B, r = body.capsules(*p0)
    ab = B - A
    x = X[free]
    t = np.clip(((x[:, None] - A[None]) * ab[None]).sum(2) / (ab * ab).sum(1)[None], 0, 1)
    dist = np.linalg.norm(x[:, None] - (A[None] + t[..., None] * ab[None]), axis=2)
    assert (r[None] + t * (body.cap_rb - body.cap_r)[None] - dist).max() < -0.003


def test_sash_follows_a_dress_forms_waist():
    """With a dress form's waist outline the sash's band lies on the outline plus its clearance, all round."""
    from humanities.dance.garment.body import Body, MotionPose
    from humanities.dance.garment.garments import sash
    from humanities.dance.garment.trials import TEST_PATTERN
    body = Body()
    body.waist, axis_x, z_w = _oval_waist(body)
    mp = MotionPose(TEST_PATTERN, body, t_raised=1.0, t_lower=1.0, t_hold=1.0)
    G = sash(body, mp(-mp.pre), clearance=0.016, h=0.02)
    P = G.attachment.parts[0].p_rest - np.array([axis_x, 0.0, z_w])         # the band, every point of it held
    assert len(P) == G.params['parts'][0]['vertices']
    ang = np.mod(np.arctan2(P[:, 1], P[:, 0]), 2 * np.pi)
    th, r = np.asarray(body.waist['theta']), np.asarray(body.waist['r'])
    assert np.allclose(np.hypot(P[:, 0], P[:, 1]), np.interp(ang, th, r, period=2 * np.pi) + 0.016, atol=2e-4)


def test_develop_keeps_a_grids_columns_and_cuts_a_bulging_band_to_its_longest_rows():
    """A grid on a cylinder develops without strain; a band bulging out at its middle rows, developed to its longest
    rows, starts with every row and column at or under its cut length (compressed, never stretched)."""
    from humanities.dance.garment.patterns import develop, grid
    th = np.linspace(0.0, 1.5, 40)
    z = np.linspace(0.0, 0.1, 9)
    for bulge, common in ((0.0, False), (0.01, True)):
        r = 0.15 + bulge * np.sin(np.pi * z / z[-1]) ** 2
        X = np.stack([r[None, :] * np.cos(th)[:, None], r[None, :] * np.sin(th)[:, None],
                      np.broadcast_to(z, (len(th), len(z)))], -1)
        pc = grid(develop(X, common_rows=common))
        Xv = X.transpose(1, 0, 2).reshape(-1, 3)
        i, j = np.arange(pc.n_vertices) % len(th), np.arange(pc.n_vertices) // len(th)
        va = np.concatenate([pc.F[:, a] for a in (0, 1, 2)])
        vb = np.concatenate([pc.F[:, b] for b in (1, 2, 0)])
        flat = np.concatenate([np.linalg.norm(pc.UV[:, b] - pc.UV[:, a], axis=1) for a, b in ((0, 1), (1, 2), (2, 0))])
        ratio = np.linalg.norm(Xv[vb] - Xv[va], axis=1) / flat
        straight = (i[va] == i[vb]) | (j[va] == j[vb])
        if bulge == 0.0:
            assert np.allclose(ratio[straight], 1.0, atol=1e-6)
            assert np.abs(ratio - 1).max() < 0.01
        else:
            assert ratio[straight].max() < 1.0 + 1e-9
            assert ratio[straight].min() < 0.95                        # the rows inside the bulge start compressed


def test_radial_envelope_reaches_the_outermost_cloth():
    """The envelope of points on a cone is at least each point's radius in its own angle and height."""
    from humanities.dance.garment.garments import radial_envelope
    rng = np.random.default_rng(1)
    th, z = rng.uniform(0, 2 * np.pi, 4000), rng.uniform(0.3, 1.0, 4000)
    r = 0.15 + 0.4 * (1.0 - z)
    env = radial_envelope(np.c_[0.02 + r * np.cos(th), r * np.sin(th), z], 0.02)
    assert (env(th, z) >= r - 1e-12).all()
    assert (env(th, z) <= r + 0.4 * 0.02 + 0.01).all()                  # a bin and its neighbours, not more
    assert env(0.0, 1.5) == 0.0                                         # no cloth there


def test_draped_sash_starts_on_its_pins_with_its_tails_clear_of_the_dressed_cloth():
    """A draped sash on a dress form's waist outline (in the rest pose, where its start is its rest) starts with its
    held edges on their targets, its band's rows and columns at or under their cut lengths and wider toward the knot,
    and its tails flat at their cut, tilted out until every point below the knot lies `margin` outside the cloth
    already dressed."""
    from humanities.dance.garment.body import Body
    from humanities.dance.garment.garments import draped_sash
    body = Body()
    body.waist, axis_x, z_w = _oval_waist(body)
    p0 = (body.origin, np.tile([1.0, 0.0, 0.0, 0.0], (len(body.names), 1)))
    env = lambda th, z: np.where(np.asarray(z) < z_w - 0.08, 0.24, 0.0)          # dressed skirts below the hip
    G = draped_sash(body, p0, h=0.012, h_tail=0.03, envelope=env, margin=0.018)
    assert np.isfinite(G.X0).all()
    assert np.allclose(G.attachment.targets(*p0), G.X0[G.pinned], atol=1e-9)
    offs = G.params['offsets']
    assert [q['name'] for q in G.params['parts']] == ['sash_band', 'sash_knot', 'sash_tail_0', 'sash_tail_1']
    w = G.params['parts'][0]['width_on_body_m']
    assert w[0] > 0.055 and w[1] > 1.4 * w[0]
    F, UV, X = G.piece.F, G.piece.UV, G.X0
    edges = lambda sel: [(F[sel][:, a], F[sel][:, b], np.linalg.norm(UV[sel][:, b] - UV[sel][:, a], axis=1))
                         for a, b in ((0, 1), (1, 2), (2, 0))]
    n_i = len(G.piece.rings['p0_start'])                                  # the band's grid: rows and columns
    for va, vb, flat in edges(np.all((F >= offs[0]) & (F < offs[1]), 1)):
        ka, kb = va - offs[0], vb - offs[0]
        straight = (ka % n_i == kb % n_i) | (ka // n_i == kb // n_i)
        assert (np.linalg.norm(X[vb] - X[va], axis=1) / flat)[straight].max() < 1.0 + 1e-6
    for t in (2, 3):                                                      # the tails: flat at their cut, clear
        for va, vb, flat in edges(np.all((F >= offs[t]) & (F < offs[t + 1]), 1)):
            assert np.allclose(np.linalg.norm(X[vb] - X[va], axis=1) / flat, 1.0, atol=1e-6)
        Xt = X[offs[t]:offs[t + 1]]
        low = Xt[:, 2] < z_w - 0.12
        r = np.hypot(Xt[low, 0] - axis_x, Xt[low, 1])
        assert low.any() and (r >= 0.24 + 0.018 - 1e-9).all()
        assert G.params['parts'][t]['tilt_deg'] > 10.0


def _gathered_tube(compression):
    """A chiffon tube cut from a narrow ring sector (as a sleeve tier is), its top gathered into a ring 0.73 of its
    width and hung from it for 1 s on the Moon; returns the largest edge stretch of each triangle in the row below the
    seam, against the cut."""
    top, L, h, ring_r, flare = 0.30, 0.12, 0.02, 0.035, 0.2
    angle = flare * top / L
    r_in = top / angle
    pc = sector(r_in, r_in + L, angle, h)
    uv = pc.vertex_uv
    rho = np.hypot(uv[:, 0], uv[:, 1])
    th = 2 * np.pi * np.mod(np.arctan2(uv[:, 1], uv[:, 0]), 2 * np.pi) / angle
    d = rho - r_in
    sa = angle / (2 * np.pi)
    r_free = top / (2 * np.pi)
    dg = 2 * (r_free - ring_r)
    slope = (r_free - ring_r) / dg + sa
    R = ring_r + (r_free - ring_r) * np.clip(d / dg, 0, 1) + d * sa
    along = np.where(d < dg, d * math.sqrt(1 - slope ** 2),
                     dg * math.sqrt(1 - slope ** 2) + (d - dg) * math.sqrt(1 - sa * sa))
    X0 = np.c_[R * np.cos(th), R * np.sin(th), 1.0 - along]
    inner = np.array(pc.rings['inner'])
    s = Settings(compression=compression)
    s.dt = 1 / 2400
    cl = Cloth(pc, X0, FABRICS['silk_chiffon'], MOON, s, pinned=inner, device='cpu')
    T = torch.as_tensor(X0[inner], dtype=cl.x.dtype)
    for _ in range(60):
        cl.set_frame(T0=T, T1=T, damping=2.0)
        for k in range(40):
            cl.substep(s.dt, k / 40, (k + 1) / 40)
    x = cl.x.double().numpy()
    F, UV = pc.F, pc.UV
    lp = np.stack([np.linalg.norm(UV[:, b] - UV[:, a], axis=1) for a, b in ((0, 1), (1, 2), (2, 0))], 1)
    lx = np.stack([np.linalg.norm(x[F[:, b]] - x[F[:, a]], axis=1) for a, b in ((0, 1), (1, 2), (2, 0))], 1)
    return (lx / lp).max(1)[d[F].mean(1) < 0.015] - 1


def test_gathered_seam_holds_the_cloth_below_at_its_cut_lengths():
    """Woven cloth buckles under a small fraction of the load it carries in tension, so the solver meets compression
    along warp and weft with 1% of the tensile stiffness. A tube gathered into a ring then keeps the row below the
    seam within a few percent of its cut lengths; with the two-sided stretch used until 3 November (compression 1)
    that row stayed stretched by about half."""
    assert Settings().compression == pytest.approx(0.01)
    assert _gathered_tube(0.01).max() < 0.05
    assert np.median(_gathered_tube(1.0)) > 0.3


def test_skirt_hem_binding_is_its_last_strip():
    """A bound hem is the skirt's last strip of triangles, every one touching the hem, with the binding's fold added
    to the hem's mass."""
    from humanities.dance.garment.body import Body, MotionPose
    from humanities.dance.garment.garments import skirt
    from humanities.dance.garment.trials import TEST_PATTERN
    body = Body()
    mp = MotionPose(TEST_PATTERN, body, t_raised=1.0, t_lower=1.0, t_hold=1.0)
    plain = skirt(body, mp(-mp.pre), cut=math.pi, length=0.4, h=0.04)
    G = skirt(body, mp(-mp.pre), cut=math.pi, length=0.4, h=0.04, binding=True, binding_kg_m=0.002)
    hem = set(G.piece.rings['outer'])
    bt = np.asarray(G.params['binding_tris'])
    assert len(bt) > 0 and all(set(G.piece.F[t]) & hem for t in bt)
    assert plain.params['binding_tris'] == [] and plain.extra_mass is None
    assert G.extra_mass.sum() == pytest.approx(0.002 * G.params['hem_circumference_m'], rel=1e-9)


def test_slits_open_a_skirt_from_its_hem_to_their_tops():
    """Side slits cut a half circle open from the hem up to slit_top down the pattern, at the requested angles: each
    slit's two lips are open edges of the cloth that meet at its closed top, start on one another, end on the hem and
    are bound like it; the pattern itself is unchanged."""
    from humanities.dance.garment.body import Body
    from humanities.dance.garment.garments import skirt
    body = Body()
    p0 = (body.origin, np.tile([1.0, 0.0, 0.0, 0.0], (len(body.names), 1)))
    kw = dict(cut=math.pi, length=0.4, h=0.04, binding=True, binding_kg_m=0.002)
    plain = skirt(body, p0, **kw)
    G = skirt(body, p0, slits=(90.0, -90.0), slit_top=0.2, **kw)
    lips = G.params['slit_lips']
    assert len(lips) == 2
    assert G.piece.n_vertices == plain.piece.n_vertices + sum(len(lo) - 1 for lo, _ in lips)
    assert G.piece.area() == pytest.approx(plain.piece.area(), rel=1e-12)
    assert len(G.piece.rings['outer']) == len(plain.piece.rings['outer']) + 2
    e = np.sort(G.piece.F[:, [0, 1, 1, 2, 2, 0]].reshape(-1, 2), 1)
    u, n = np.unique(e, axis=0, return_counts=True)
    open_edges = {tuple(x) for x in u[n == 1].tolist()}
    uv = G.piece.vertex_uv
    r_in = G.params['pattern_inner_radius_m']
    sp = {s['name']: s for s in body.spec['segments']}
    axis_x = sp['abdomen']['origin'][0] + sp['abdomen']['capsules'][0]['a'][0]
    hem = set(G.piece.rings['outer'])
    for (lo, up), deg in zip(lips, (90.0, -90.0)):
        assert lo[0] == up[0] and len(lo) == len(up) > 3
        for lip in (lo, up):
            assert all((min(a, b), max(a, b)) in open_edges for a, b in zip(lip[:-1], lip[1:]))
            assert lip[-1] in hem
        assert not any(set(t) >= {a, b} for t in G.piece.F.tolist() for a, b in zip(lo[1:], up[1:]))
        assert np.allclose(G.X0[lo], G.X0[up])
        top = np.hypot(*uv[lo[0]]) - r_in
        assert 0.2 - 1e-9 <= top < 0.2 + 0.04
        X = G.X0[lo[0]]
        assert abs(math.degrees(math.atan2(X[1], X[0] - axis_x)) - deg) < 6.0
    lip_v = {v for lo, up in lips for v in lo[1:] + up[1:]}
    bt = set(G.params['binding_tris'])
    assert all(k in bt for k, t in enumerate(G.piece.F.tolist()) if lip_v & set(t))
    assert G.extra_mass.sum() == pytest.approx(plain.extra_mass.sum() + 4 * 0.002 * 0.2, rel=1e-9)


def test_tacked_sash_holds_its_top_edge_and_turns_as_one_piece():
    """hold='top' keeps the band's top edge held all round and tacks its bottom edge about every tack_every; with one
    share, a twist of the torso turns every held point of the sash (band, knot, tails) by the same angle, where the
    default shares turn the top edge further than the bottom."""
    from humanities.dance.garment.body import Body
    from humanities.dance.garment.garments import draped_sash
    body = Body()
    body.waist, axis_x, z_w = _oval_waist(body)
    p0 = (body.origin, np.tile([1.0, 0.0, 0.0, 0.0], (len(body.names), 1)))
    kw = dict(h=0.012, h_tail=0.03)
    E = draped_sash(body, p0, **kw)
    T = draped_sash(body, p0, hold='top', tack_every=0.08, share=0.8, **kw)
    assert T.piece.n_vertices == E.piece.n_vertices
    pins = set(np.asarray(T.pinned).tolist())
    top, bottom = T.piece.rings['p0_start'], T.piece.rings['p0_end']
    assert set(top) <= pins and set(bottom) <= set(np.asarray(E.pinned).tolist())
    tacks = [v for v in bottom if v in pins]
    assert 6 <= len(tacks) <= len(bottom) / 5                                # every 7th point at h = 12 mm
    gaps = np.linalg.norm(np.diff(T.X0[tacks], axis=0), axis=1)
    assert 0.06 < np.median(gaps) < 0.10
    assert np.allclose(T.attachment.targets(*p0), T.X0[T.pinned], atol=1e-9)
    quat = p0[1].copy()
    k = body.index['abdomen']
    quat[k] = [math.cos(math.radians(10.0)), 0.0, 0.0, math.sin(math.radians(10.0))]   # 20 degrees about the vertical
    o = body.origin[k]
    turned = lambda G: np.degrees(np.angle(np.exp(1j * (
        np.arctan2(*(G.attachment.targets(p0[0], quat) - o)[:, [1, 0]].T)
        - np.arctan2(*(G.X0[G.pinned] - o)[:, [1, 0]].T)))))
    t_one, t_rows = turned(T), turned(E)
    assert np.ptp(t_one) < 1e-6 and t_one[0] == pytest.approx(math.degrees(np.angle(0.8 + 0.2 * np.exp(1j * math.radians(20.0)))))
    assert np.ptp(t_rows) > 1.0
