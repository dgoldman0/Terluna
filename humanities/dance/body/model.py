"""The dancer's body: segments, joints, masses, ranges and strength, written once.

    python -m humanities.dance.body.model     # writes body.json (the spec) and dancer.xml (the MuJoCo model)

One description feeds every consumer. The motion engine builds its Pinocchio model from the spec, MuJoCo replays the
optimised motion on dancer.xml, the garment simulation reads the collision capsules and attachment frames, and the
renderer fits its armature to the joint centres.

Frame and rest pose. x points forward, y to the dancer's left, z up; metres, kilograms, seconds. In the rest pose the
dancer stands upright with the feet parallel under the hips, the arms hanging at the sides with the palms toward the
thighs, and the eyes forward. Every segment's frame sits at its proximal joint centre and is aligned with the world in
the rest pose, so a segment's world rotation is also its rotation from rest.

The body. A female dancer of 1.64 m, the stature of the CC0 realistic female base mesh the renderer skins (Blender
Studio's Human Base Meshes v1.4.1), and 55 kg. Segment lengths, masses, centres of mass and radii of gyration follow
de Leva (1996)'s adjustment of Zatsiorsky–Seluyanov's female data (Table 4, as transcribed in sources.json), scaled by
stature: the trunk and head take the rows measured from the cervicale (CERV), whose lengths stack the landmarks, with
their own centre-of-mass shares and radii; a sagittal radius turns about the front-back axis (x), a transverse radius
about the side-to-side axis (y). The joint centres between his landmarks (shoulder, spine and sternoclavicular joints,
the foot's sole points) are anthropometric placements recorded beside each value. Ranges of motion are a trained dancer's, wider than the clinical norms at the hip, spine and ankle.
Joint strength is a set of peak torques for a 55 kg trained woman; the motion engine bounds its torques by them and the
checks report how close each motion comes. Every number here is a modelling choice for one representative body; the
literature record (sources.json) holds what each rests on.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPEC_OUT = HERE / 'body.json'
MJCF_OUT = HERE / 'dancer.xml'
SCHEMA = 'terluna.humanities.dance-body/1'

STATURE_M = 1.64
MASS_KG = 55.0
H = STATURE_M

# --- landmarks in the rest pose (de Leva 1996, female: segment lengths as fractions of his 1.735 m sample's stature)
VERT = H                                   # top of the head
CERV = H * (1 - 243.7 / 1735)              # cervicale, C7; head and neck 243.7 mm
XYPH = CERV - H * 228.0 / 1735             # xiphion; upper trunk 228.0 mm
OMPH = XYPH - H * 205.3 / 1735             # omphalion; middle trunk 205.3 mm
MIDH = OMPH - H * 181.5 / 1735             # mid-point of the hip joint centres; lower trunk 181.5 mm
KJC = MIDH - H * 368.5 / 1735              # knee joint centre; thigh 368.5 mm
LMAL = KJC - H * 432.3 / 1735              # lateral malleolus; shank 432.3 mm
SJC = 1.300                                # shoulder joint centre: 4 cm below the acromion (0.818 H), a placement
EJC = SJC - H * 275.1 / 1735               # elbow; upper arm 275.1 mm
WJC = EJC - H * 264.3 / 1735               # wrist; forearm 264.3 mm
MET3 = WJC - H * 78.0 / 1735               # third metacarpal head; hand 78.0 mm
FOOT_LEN = 0.230                           # heel to toe tip; de Leva's 228.3 mm scaled gives 0.216, the mesh's foot 0.23
HIP_Y = 0.082                              # half the distance between the hip joint centres (placement)
SHOULDER_Y = 0.170                         # half the distance between the shoulder joint centres (placement)
ANKLE_Z = LMAL - 0.002                     # the ankle's axis passes just below the malleolus (placement)


def deg(a):
    return a * math.pi / 180


@dataclass
class Joint:
    """One joint of a segment. A hinge turns about one axis (in the parent's frame at rest, positive in the anatomical
    sense named). A ball joint turns freely; its range is a swing cone of the limb's long axis and a twist range."""
    name: str
    kind: str                               # 'hinge' or 'ball'
    axis: tuple = (0.0, 0.0, 0.0)           # hinge axis
    range_deg: tuple = (0.0, 0.0)           # hinge range
    torque_nm: tuple = (0.0, 0.0)           # hinge torque bounds (negative-direction, positive-direction), magnitudes
    swing_deg: dict = field(default_factory=dict)   # ball: largest swing of the long axis toward each direction
    twist_deg: tuple = (0.0, 0.0)           # ball: twist range about the long axis (internal, external)
    ball_torque_nm: tuple = (0.0, 0.0, 0.0) # ball: torque bound about x, y, z of the segment frame
    long_axis: tuple = (0.0, 0.0, -1.0)     # ball: the limb's long axis at rest


@dataclass
class Segment:
    name: str
    parent: str | None
    origin: tuple                           # proximal joint centre, world, rest pose
    joints: list                            # applied in order from the parent outward
    mass_share: float                       # fraction of body mass
    com: tuple                              # world, rest pose
    gyration: tuple                         # radii of gyration about x, y, z of the segment frame (m)
    capsules: list = field(default_factory=list)   # (a, b, radius), world, rest pose
    tip: tuple | None = None                # distal end for bones that have no child (hand, toes, head)


def side_name(s):
    return 'L' if s > 0 else 'R'


def build_segments():
    segs = []
    m = lambda pct: pct / 100.0

    # Pelvis: de Leva's lower trunk (OMPH-MIDH), centre of mass 49.20% of its length from OMPH.
    lpt = OMPH - MIDH
    segs.append(Segment('pelvis', None, (0.0, 0.0, MIDH), [], m(12.47),
                        (-0.010, 0.0, OMPH - 0.4920 * lpt),
                        (0.433 * lpt, 0.402 * lpt, 0.444 * lpt),
                        capsules=[((-0.015, 0.085, MIDH + 0.035), (-0.015, -0.085, MIDH + 0.035), 0.100),
                                  ((-0.020, 0.0, MIDH - 0.010), (-0.020, 0.0, OMPH - 0.020), 0.115)]))
    # Abdomen: middle trunk (XYPH-OMPH), centre of mass 45.12% from XYPH. Lumbar joint at L3-L4, behind the navel.
    mpt = XYPH - OMPH
    lumbar = [Joint('lumbar_flex', 'hinge', (0, 1, 0), (-35, 60), (100, 120)),
              Joint('lumbar_bend', 'hinge', (-1, 0, 0), (-30, 30), (90, 90)),
              Joint('lumbar_rot', 'hinge', (0, 0, 1), (-20, 20), (60, 60))]
    segs.append(Segment('abdomen', 'pelvis', (-0.040, 0.0, OMPH), lumbar, m(14.65),
                        (-0.010, 0.0, XYPH - 0.4512 * mpt),
                        (0.433 * mpt, 0.354 * mpt, 0.415 * mpt),
                        capsules=[((-0.005, 0.0, OMPH + 0.010), (-0.005, 0.0, XYPH - 0.010), 0.108)]))
    # Thorax: upper trunk (CERV-XYPH) less the shoulder girdles, centre of mass 50.50% of its length below CERV.
    upt = CERV - XYPH
    thoracic = [Joint('thoracic_flex', 'hinge', (0, 1, 0), (-30, 40), (90, 100)),
                Joint('thoracic_bend', 'hinge', (-1, 0, 0), (-25, 25), (80, 80)),
                Joint('thoracic_rot', 'hinge', (0, 0, 1), (-40, 40), (60, 60))]
    segs.append(Segment('thorax', 'abdomen', (-0.050, 0.0, XYPH), thoracic, m(15.45 - 0.80),
                        (-0.020, 0.0, CERV - 0.5050 * upt),
                        (0.466 * upt, 0.314 * upt, 0.449 * upt),
                        capsules=[((-0.010, 0.075, XYPH + 0.115), (-0.010, -0.075, XYPH + 0.115), 0.115),
                                  ((-0.020, 0.0, XYPH + 0.020), (-0.020, 0.0, CERV - 0.060), 0.110)]))
    # Neck: one per cent of the body, from C7-T1 to the skull's base.
    neck = [Joint('neck_flex', 'hinge', (0, 1, 0), (-50, 45), (25, 25)),
            Joint('neck_bend', 'hinge', (-1, 0, 0), (-35, 35), (20, 20)),
            Joint('neck_rot', 'hinge', (0, 0, 1), (-70, 70), (15, 15))]
    head_base = 0.915 * H
    segs.append(Segment('neck', 'thorax', (-0.045, 0.0, CERV - 0.010), neck, m(1.00),
                        (-0.030, 0.0, 0.5 * (CERV + head_base)), (0.035, 0.035, 0.030),
                        capsules=[((-0.030, 0.0, CERV - 0.010), (-0.020, 0.0, head_base), 0.050)]))
    # Head: de Leva's head and neck (VERT-CERV) less the neck; its centre of mass 48.41% of VERT-CERV below the vertex.
    hn = VERT - CERV
    segs.append(Segment('head', 'neck', (-0.015, 0.0, head_base),
                        [Joint('head_nod', 'hinge', (0, 1, 0), (-25, 25), (10, 10))], m(5.68),
                        (0.005, 0.0, VERT - 0.4841 * hn), (0.271 * hn, 0.295 * hn, 0.261 * hn),
                        capsules=[((0.010, 0.0, VERT - 0.095), (0.000, 0.0, VERT - 0.070), 0.092)],
                        tip=(0.0, 0.0, VERT)))

    for s in (1, -1):
        L = side_name(s)
        # Shoulder girdle: sternoclavicular joint to the shoulder joint centre.
        sc = (0.035, s * 0.020, CERV - 0.035)
        clav = [Joint(f'{L}_clav_elev', 'hinge', (s, 0, 0), (-10, 35), (40, 40)),
                Joint(f'{L}_clav_prot', 'hinge', (0, 0, -s), (-15, 20), (40, 40))]
        segs.append(Segment(f'{L}_clavicle', 'thorax', sc, clav, m(0.40),
                            (0.010, s * 0.090, CERV - 0.050), (0.012, 0.040, 0.040),
                            capsules=[(sc, (0.000, s * (SHOULDER_Y - 0.020), SJC + 0.040), 0.030)]))
        # Upper arm, ball joint at the shoulder; centre of mass 57.54% down.
        ua = SJC - EJC
        segs.append(Segment(f'{L}_upper_arm', f'{L}_clavicle', (-0.010, s * SHOULDER_Y, SJC),
                            [Joint(f'{L}_shoulder', 'ball', swing_deg=dict(forward=180, back=60, out=180, across=40),
                                   twist_deg=(90, 90), ball_torque_nm=(45, 50, 30))],
                            m(2.55), (-0.010, s * (SHOULDER_Y + 0.004), SJC - 0.5754 * ua),
                            (0.278 * ua, 0.260 * ua, 0.148 * ua),
                            capsules=[((-0.010, s * SHOULDER_Y, SJC - 0.010), (-0.010, s * 0.178, EJC + 0.020), 0.043)]))
        # Forearm: elbow flexion, then pronation about the long axis. Its radii keep the transverse value on x, which
        # suits the rest pose's palms toward the thighs (de Leva's data hold the palm forward).
        fa = EJC - WJC
        segs.append(Segment(f'{L}_forearm', f'{L}_upper_arm', (-0.010, s * 0.178, EJC),
                            [Joint(f'{L}_elbow', 'hinge', (0, -1, 0), (0, 150), (35, 35)),
                             Joint(f'{L}_pronation', 'hinge', (0, 0, -s), (-80, 80), (8, 8))],
                            m(1.38), (-0.005, s * 0.180, EJC - 0.4559 * fa), (0.257 * fa, 0.261 * fa, 0.094 * fa),
                            capsules=[((-0.010, s * 0.178, EJC - 0.010), (0.003, s * 0.182, WJC + 0.015), 0.034)]))
        # Hand: wrist flexion and radial deviation; centre of mass 74.74% of WJC-MET3.
        hd = WJC - MET3
        segs.append(Segment(f'{L}_hand', f'{L}_forearm', (0.005, s * 0.183, WJC),
                            [Joint(f'{L}_wrist_flex', 'hinge', (-s, 0, 0), (-70, 80), (15, 15)),
                             Joint(f'{L}_wrist_dev', 'hinge', (0, -1, 0), (-30, 20), (12, 12))],
                            m(0.56), (0.008, s * 0.183, WJC - 0.7474 * hd), (0.454 * hd, 0.531 * hd, 0.335 * hd),
                            capsules=[((0.008, s * 0.183, WJC - 0.010), (0.010, s * 0.183, 0.377 * H + 0.025), 0.024)],
                            tip=(0.010, s * 0.183, 0.377 * H)))
        # Thigh, ball joint at the hip; centre of mass 36.12% down.
        th = MIDH - KJC
        segs.append(Segment(f'{L}_thigh', 'pelvis', (0.0, s * HIP_Y, MIDH),
                            [Joint(f'{L}_hip', 'ball', swing_deg=dict(forward=150, back=35, out=110, across=30),
                                   twist_deg=(40, 65), ball_torque_nm=(100, 150, 50))],
                            m(14.78), (0.005, s * HIP_Y, MIDH - 0.3612 * th), (0.369 * th, 0.364 * th, 0.162 * th),
                            capsules=[((0.000, s * (HIP_Y + 0.008), MIDH - 0.030), (0.000, s * HIP_Y, KJC + 0.040), 0.072)]))
        # Shank: knee flexion; centre of mass 44.16% down.
        sh = KJC - LMAL
        segs.append(Segment(f'{L}_shank', f'{L}_thigh', (0.0, s * HIP_Y, KJC),
                            [Joint(f'{L}_knee', 'hinge', (0, 1, 0), (-5, 150), (178, 90))],
                            m(4.81), (-0.010, s * HIP_Y, KJC - 0.4416 * sh), (0.271 * sh, 0.267 * sh, 0.093 * sh),
                            capsules=[((-0.005, s * HIP_Y, KJC - 0.020), (-0.005, s * HIP_Y, ANKLE_Z + 0.030), 0.048)]))
        # Foot: ankle plantarflexion, then subtalar inversion; heel to the metatarsal heads.
        mtp_x = -0.050 + 0.78 * FOOT_LEN
        segs.append(Segment(f'{L}_foot', f'{L}_shank', (0.0, s * HIP_Y, ANKLE_Z),
                            [Joint(f'{L}_ankle', 'hinge', (0, 1, 0), (-25, 85), (29, 120)),
                             Joint(f'{L}_subtalar', 'hinge', (-s, 0, 0), (-20, 30), (30, 30))],
                            m(1.00), (0.007, s * HIP_Y, 0.035),       # with the toes, 40.14% of heel to toe tip
                            (0.139 * FOOT_LEN, 0.279 * FOOT_LEN, 0.299 * FOOT_LEN),
                            capsules=[((-0.030, s * HIP_Y, 0.032), (mtp_x - 0.010, s * HIP_Y, 0.024), 0.030)]))
        # Toes: metatarsophalangeal extension, which takes the dancer onto demi-pointe.
        segs.append(Segment(f'{L}_toes', f'{L}_foot', (mtp_x, s * HIP_Y, 0.022),
                            [Joint(f'{L}_mtp', 'hinge', (0, -1, 0), (-30, 85), (20, 20))],
                            m(0.29), (mtp_x + 0.035, s * HIP_Y, 0.012), (0.012, 0.020, 0.020),
                            capsules=[((mtp_x + 0.005, s * HIP_Y, 0.014), (-0.050 + FOOT_LEN - 0.012, s * HIP_Y, 0.012),
                                       0.013)],
                            tip=(-0.050 + FOOT_LEN, s * HIP_Y, 0.010)))
    return segs


def contact_points():
    """Points of the sole and hand that can take the floor's force, world coordinates at rest. Sole points lie on the
    floor (z = 0) at rest; the heel pair and the two metatarsal heads span the foot's support, the toe point carries
    demi-pointe."""
    pts = {}
    mtp_x = -0.050 + 0.78 * FOOT_LEN
    for s in (1, -1):
        L = side_name(s)
        pts[f'{L}_heel_m'] = (f'{L}_foot', (-0.040, s * (HIP_Y - 0.018), 0.0))
        pts[f'{L}_heel_l'] = (f'{L}_foot', (-0.040, s * (HIP_Y + 0.018), 0.0))
        pts[f'{L}_ball_m'] = (f'{L}_foot', (mtp_x, s * (HIP_Y - 0.032), 0.0))
        pts[f'{L}_ball_l'] = (f'{L}_foot', (mtp_x - 0.025, s * (HIP_Y + 0.036), 0.0))
        pts[f'{L}_toe'] = (f'{L}_toes', (-0.050 + FOOT_LEN - 0.010, s * (HIP_Y - 0.010), 0.0))
        pts[f'{L}_palm'] = (f'{L}_hand', (0.008, s * 0.183, MET3 - 0.020))
    return pts


def segments_by_name(segs=None):
    return {s.name: s for s in (segs or build_segments())}


def total_mass(segs=None):
    return MASS_KG * sum(s.mass_share for s in (segs or build_segments()))


def rel(a, b):
    return tuple(round(x - y, 6) for x, y in zip(a, b))


def spec() -> dict:
    segs = build_segments()
    by = segments_by_name(segs)
    out = dict(schema=SCHEMA,
               frame='x forward, y to the dancer\'s left, z up; metres, kilograms, seconds; rest pose standing, feet '
                     'parallel under the hips, arms at the sides, palms toward the thighs; every segment frame sits at '
                     'its proximal joint centre and is world-aligned at rest',
               stature_m=STATURE_M, mass_kg=MASS_KG,
               evidence='One representative trained female dancer: de Leva (1996) segment parameters scaled by '
                        'stature, anthropometric joint placements, dancer ranges of motion and peak torques chosen for '
                        'a 55 kg trained woman (sources.json). A modelling basis for the motion engine and renderer, '
                        'not a measured individual.',
               segments=[], contacts={})
    for s in segs:
        joints = []
        for j in s.joints:
            if j.kind == 'hinge':
                joints.append(dict(name=j.name, kind='hinge', axis=list(j.axis), range_deg=list(j.range_deg),
                                   torque_nm=list(j.torque_nm)))
            else:
                joints.append(dict(name=j.name, kind='ball', long_axis=list(j.long_axis), swing_deg=j.swing_deg,
                                   twist_deg=list(j.twist_deg), torque_nm=list(j.ball_torque_nm)))
        mass = MASS_KG * s.mass_share
        out['segments'].append(dict(
            name=s.name, parent=s.parent, origin=[round(x, 6) for x in s.origin],
            offset=list(rel(s.origin, by[s.parent].origin)) if s.parent else [round(x, 6) for x in s.origin],
            joints=joints, mass_kg=round(mass, 6), com=list(rel(s.com, s.origin)),
            inertia_kg_m2=[round(mass * r * r, 8) for r in s.gyration],
            capsules=[dict(a=list(rel(a, s.origin)), b=list(rel(b, s.origin)), radius=r) for a, b, r in s.capsules],
            tip=list(rel(s.tip, s.origin)) if s.tip else None))
    for name, (seg, p) in contact_points().items():
        out['contacts'][name] = dict(segment=seg, pos=list(rel(p, by[seg].origin)))
    return out


# --- MuJoCo model -------------------------------------------------------------------------------------------------

def _fmt(v):
    return ' '.join(f'{x:.6g}' for x in v)


def mjcf(gravity: float = 1.6242, density: float = 1.40, viscosity: float = 1.8e-5, friction: float = 0.6,
         timestep: float = 0.001) -> str:
    """The MuJoCo model. Ball joints become MuJoCo ball joints driven by three motors each; hinges keep their ranges;
    torques are bounded by the strength table. Sole and palm points are small spheres, the body's capsules collide
    with the floor only. Fluid forces use MuJoCo's ellipsoid model on every capsule at the given air."""
    sp = spec()
    children = {}
    for s in sp['segments']:
        children.setdefault(s['parent'], []).append(s)
    lines = [f'<mujoco model="terluna_dancer">',
             f'  <option timestep="{timestep}" gravity="0 0 {-gravity:.6g}" density="{density:.6g}" '
             f'viscosity="{viscosity:.6g}" integrator="implicitfast"/>',
             '  <compiler angle="radian" autolimits="true"/>',
             '  <default>',
             f'    <geom contype="0" conaffinity="0" friction="{friction} 0.02 0.002" fluidshape="ellipsoid"/>',
             '    <joint damping="0.5" armature="0.01"/>',
             '  </default>',
             '  <worldbody>',
             f'    <geom name="floor" type="plane" size="50 50 0.1" contype="1" conaffinity="1" friction="{friction} 0.02 0.002" fluidshape="none"/>']
    contacts_by_seg = {}
    for name, c in sp['contacts'].items():
        contacts_by_seg.setdefault(c['segment'], []).append((name, c['pos']))
    actuators = []

    def body(s, indent):
        pad = ' ' * indent
        lines.append(f'{pad}<body name="{s["name"]}" pos="{_fmt(s["offset"])}">')
        if s['parent'] is None:
            lines.append(f'{pad}  <freejoint name="root"/>')
        for j in s['joints']:
            if j['kind'] == 'hinge':
                lo, hi = (deg(x) for x in j['range_deg'])
                lines.append(f'{pad}  <joint name="{j["name"]}" type="hinge" axis="{_fmt(j["axis"])}" range="{lo:.6g} {hi:.6g}"/>')
                actuators.append(f'    <motor name="{j["name"]}" joint="{j["name"]}" ctrlrange="{-j["torque_nm"][0]:.6g} {j["torque_nm"][1]:.6g}"/>')
            else:
                lim = deg(max(j['swing_deg'].values()))
                lines.append(f'{pad}  <joint name="{j["name"]}" type="ball" range="0 {lim:.6g}"/>')
                for k, ax in enumerate(('x', 'y', 'z')):
                    g = ['0'] * 6
                    g[k] = '1'
                    t = j['torque_nm'][k]
                    actuators.append(f'    <motor name="{j["name"]}_{ax}" joint="{j["name"]}" gear="{" ".join(g)}" ctrlrange="{-t:.6g} {t:.6g}"/>')
        I = s['inertia_kg_m2']
        lines.append(f'{pad}  <inertial pos="{_fmt(s["com"])}" mass="{s["mass_kg"]:.6g}" diaginertia="{_fmt(I)}"/>')
        sole = s['name'].endswith(('_foot', '_toes'))        # the sole points carry the feet's contact
        for k, c in enumerate(s['capsules']):
            lines.append(f'{pad}  <geom name="{s["name"]}_cap{k}" type="capsule" fromto="{_fmt(c["a"] + c["b"])}" '
                         f'size="{c["radius"]:.6g}" mass="0" contype="{0 if sole else 2}" conaffinity="{0 if sole else 1}" '
                         f'rgba="0.8 0.7 0.6 1"/>')
        for name, p in contacts_by_seg.get(s['name'], []):
            lines.append(f'{pad}  <geom name="{name}" type="sphere" pos="{_fmt(p[:2] + [p[2] + 0.006])}" size="0.006" '
                         f'mass="0" contype="1" conaffinity="1" fluidshape="none" rgba="0.2 0.8 0.3 1"/>')
            lines.append(f'{pad}  <site name="{name}" pos="{_fmt(p)}" size="0.005"/>')
        for ch in children.get(s['name'], []):
            body(ch, indent + 2)
        lines.append(f'{pad}</body>')

    for root in children[None]:
        body(root, 4)
    lines.append('  </worldbody>')
    lines.append('  <contact>')
    lines.append('  </contact>')
    lines.append('  <actuator>')
    lines.extend(actuators)
    lines.append('  </actuator>')
    lines.append('</mujoco>')
    return '\n'.join(lines) + '\n'


def main() -> int:
    SPEC_OUT.write_text(json.dumps(spec(), indent=1) + '\n')
    MJCF_OUT.write_text(mjcf())
    print(f'wrote {SPEC_OUT.relative_to(HERE.parents[2])} and {MJCF_OUT.relative_to(HERE.parents[2])}; '
          f'{len(build_segments())} segments, {total_mass():.2f} kg')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
