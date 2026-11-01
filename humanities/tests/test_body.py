"""The body model's spec and MuJoCo model regenerate from model.py, hold the stated body, and stand on the floor."""
import json

import pytest

from humanities.dance.body import model as bm


def test_spec_and_mjcf_regenerate():
    assert json.loads(bm.SPEC_OUT.read_text()) == json.loads(json.dumps(bm.spec()))
    assert bm.MJCF_OUT.read_text() == bm.mjcf()


def test_masses_proportions_and_rest_geometry():
    sp = bm.spec()
    assert sum(s['mass_kg'] for s in sp['segments']) == pytest.approx(bm.MASS_KG, abs=0.02)
    by = {s['name']: s for s in sp['segments']}
    assert by['pelvis']['parent'] is None and len(sp['segments']) == 21
    # de Leva's landmarks scaled to 1.64 m: hip, knee and shoulder heights against the usual shares of stature.
    assert 0.49 < bm.MIDH / bm.H < 0.53 and 0.27 < bm.KJC / bm.H < 0.30 and 0.78 < bm.SJC / bm.H < 0.81
    # Every sole point lies on the floor at rest, and every segment's centre of mass is inside the body's span.
    for name, c in sp['contacts'].items():
        world = [o + p for o, p in zip(by[c['segment']]['origin'], c['pos'])]
        if 'palm' not in name:
            assert world[2] == pytest.approx(0.0, abs=1e-9)
    for s in sp['segments']:
        z = s['origin'][2] + s['com'][2]
        assert 0.0 < z < bm.STATURE_M
    for s in sp['segments']:
        for j in s['joints']:
            if j['kind'] == 'hinge':
                lo, hi = j['range_deg']
                assert lo < hi and min(j['torque_nm']) > 0


def test_mujoco_loads_the_model_with_one_motor_per_degree_of_freedom():
    mujoco = pytest.importorskip('mujoco')
    m = mujoco.MjModel.from_xml_path(str(bm.MJCF_OUT))
    assert m.nv == 48 and m.nu == 42
    assert m.body_subtreemass[1] == pytest.approx(bm.MASS_KG, abs=0.02)
