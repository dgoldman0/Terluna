"""The terms the sport study's contests share: the bodies, a player's state, the floor's grips, the hands' reach, and
the contracts between the reach engine (reach.py), the ball windows (windows.py), the athletic vocabulary
(vocabulary.py) and the contests.

    from humanities.sport.state import BODIES, GRIPS, WORLDS, State, Reach

The second round of the sport study (8 November 2026) asks which contests stay open when bodies meet a ball on a floor
at lunar gravity: who can be where, and when, against the time a ball takes to get there. A player is a centre of mass
over a floor with a posture, a horizontal velocity and, when aloft, an arc it cannot change; the hands reach a band
around it. Every figure is per body (BODIES, from the envelope), per world (WORLDS) and per grip (GRIPS: the floor's
friction coefficient, which the contests sweep from the dance's vinyl to a high-grip floor instead of fixing one).

Coordinates are metres: x and y on the floor, z up from it. Headings are radians anticlockwise from +x.

Assumptions, each a value chosen for this screen (no source retrieved this round; sources.json lists them):
- Reaction. A player starts to act 0.20 s after the cue (a release, an opponent's take-off), on either world: the
  nervous system keeps Earth's tempo, like the limbs.
- A quick step. Placing a foot takes 0.25 s, on the muscles' clock.
- Top speed. The legs' swing caps running at 7 m/s on either world. On the Moon the floor's grip binds long before:
  reaching 7 m/s at a grip of 0.8 takes 5.4 s and 19 m.
- The legs' horizontal push. A sprint start's force-velocity profile, F0 (1 - v/V0) with F0 7 N/kg and V0 8.5 m/s,
  the same per kilogram for both bodies (a team-sport athlete's profile from general knowledge of sprint profiling;
  not retrieved). On Earth it binds before the grip; on the Moon the grip binds below about 6-7 m/s.
- Jumps. The take-off follows the vertical force-velocity profile (the dance envelope's) over the push's stroke, and a
  jump sinks only as deep as its take-off speed needs, at least 0.05 m from upright; from the held crouch it needs no
  dip up to the crouch's best take-off.
- The hands. Standing, the fingertips reach 1.33 statures up (2.18 m for the dance body), 0.60 m out from the centre
  of mass, and down to 0.50 m above the floor without a dip; the larger body scales by stature. Aloft, the band rides
  with the centre of mass.
- Catching. A catch takes the ball's momentum over the hand's 0.15 s give, as a throw gives it; a longer give halves
  the push. The arm's own fastest release is 20 m/s, beyond what most stances hold on the Moon.
- Contests. Two players whose hands can reach a ball within 0.30 s of each other contest it; it goes to hands and
  bodies, beyond what these screens decide.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

from humanities.dance import envelope as dance
from humanities.sport import envelope

G_M = envelope.G_M
G_E = envelope.G_E
WORLDS = envelope.WORLDS                  # (('moon', g), ('earth', g))
GRIPS = (0.5, 0.8, 1.2)                   # the floor's grip in every contest: the dance's vinyl, an athletic floor, high grip
MU = envelope.MU                          # the central grip
RHO_MOON = dance.air_density(0.0)         # the design run's air at the sea-level datum
RHO_EARTH = envelope.RHO_EARTH

REACTION_S = 0.20
STEP_S = 0.25
V_MAX_M_S = 7.0
SPRINT_F0_N_KG = 7.0
SPRINT_V0_M_S = 8.5
DIP_M = 0.29                              # a standing jump's countermovement (the solved salto's dip; envelope.jumps)
MIN_DIP_M = 0.05                          # the least a jump from upright sinks first
CONTEST_S = 0.30                          # two players reaching a ball within this of each other contest it
REACH_UP_RATIO = 1.33
REACH_OUT_M = 0.60                        # the dance body; the larger body scales by stature
REACH_LOW_M = 0.50
CATCH_S = envelope.THROW_S
RELEASE_MAX_M_S = 20.0

POSTURES = ('upright', 'ready', 'split', 'moving', 'aloft', 'down')
"""upright: standing, feet under the body (the envelope's 0.13 m beside and 0.10 m ahead or behind).
ready: the wide ready crouch, feet 0.70 m apart, centre of mass 0.75 m up (dance body).
split: the long split stance facing `facing`, the rear foot 0.50 m behind the centre of mass, crouched as ready.
moving: already in a gait, so the feet can be placed where the next push needs them; `next_stance_s` is the time
    until the next foot is down (0 for a shuffle; up to a stride's flight when bounding).
aloft: in an arc it cannot change until it lands; `z` and `vz` give the arc.
down: on the floor after a fall; it must get up before anything else."""

UNITS = (
    # getting set and going
    'ready_crouch', 'lean_start', 'drop_step', 'shuffle', 'bound', 'bounce_cut', 'braking_landing', 'pivot',
    # jumps
    'standing_jump', 'crouch_jump', 'rebound_jump', 'side_hop', 'running_leap',
    # the ball
    'chest_pass', 'split_throw', 'landing_throw', 'air_pass', 'catch', 'scoop', 'kick',
    # contact and the floor
    'shove', 'brace', 'boost', 'fall_and_rise',
)
"""The athletic vocabulary's unit ids (vocabulary/units.json). The reach engine and the contests name the moves they use
by these ids, in order, in Reach.moves and in their results."""


@dataclass(frozen=True)
class Body:
    key: str                              # 'dancer' or 'larger'
    name: str                             # the envelope's name for the body
    mass_kg: float
    stature_m: float
    com_height_m: float                   # standing
    ready_com_m: float
    ready_half_stance_m: float
    side_reach_m: float
    fore_aft_reach_m: float
    split_reach_m: float
    push_off_m: float                     # the full push, from the dip's bottom to take-off
    release_above_com_m: float
    inertia_kg_m2: dict = field(hash=False)
    reach_up_m: float = 0.0
    reach_out_m: float = 0.0
    reach_low_m: float = 0.0

    def com_m(self, posture):
        """The centre of mass's height on the floor in a posture."""
        return self.ready_com_m if posture in ('ready', 'split') else self.com_height_m

    @property
    def dip_m(self):
        """A standing jump's countermovement for this body (the dance body's, scaled by stature)."""
        return DIP_M * self.stature_m / envelope.STATURE_M

    @property
    def crouch_push_m(self):
        """The push from the held ready crouch to the same take-off height as a full jump, without a dip."""
        return self.com_height_m - self.ready_com_m + self.push_off_m - self.dip_m

    def take_off_m_s(self, g, push_m=None):
        """The fastest vertical take-off over a push (Samozino's force-velocity model, the dance envelope's profile)."""
        return dance.samozino_take_off(g, h=self.push_off_m if push_m is None else push_m)

    def hand_band(self, com_z):
        """The heights the hands reach with the centre of mass at com_z: the standing band, carried with the body."""
        return com_z - (self.com_height_m - self.reach_low_m), com_z + (self.reach_up_m - self.com_height_m)

    def push_for(self, vz, g):
        """The shortest push that reaches take-off speed vz at full effort: Samozino's model inverted,
        v^2 = 2 s (F0 (1 - v/(2 V0)) - g) solved for the stroke s (inf where the profile cannot reach vz)."""
        net = dance.F0_N_KG * (1 - vz / (2 * dance.V0_M_S)) - g
        return vz * vz / (2 * net) if net > 0 else math.inf

    def jump_timing(self, posture, vz, g):
        """The quickest take-off at vertical speed vz from a posture: (dip below the posture in m, the dip's seconds,
        the push's seconds), or None where the body cannot reach vz. The legs leave the floor at the take-off height, on
        the toes, push_off - dip above standing. A full-effort push of the stroke s that vz needs (push_for) covers it
        at the mean speed vz/2, and the legs coast at vz through whatever stroke S the posture leaves beyond it, so the
        push takes (S + s)/vz. From upright the body sinks as far as s needs and at least MIN_DIP_M, no faster than g
        (envelope.dip_s). From the ready crouch or the split stance the crouch's stroke is S, and the body sinks further
        only where s exceeds it. Moving, the landing's or the stride's stroke is the countermovement: no dip, and the
        stance lasts a stride's (envelope.STRIDE_S)."""
        s = self.push_for(vz, g)
        if vz <= 0 or s > self.push_off_m + 1e-9:
            return None
        if posture == 'moving':
            return 0.0, 0.0, envelope.STRIDE_S
        if posture in ('ready', 'split'):
            dip = s - self.crouch_push_m if s - self.crouch_push_m > 1e-9 else 0.0
            S = self.crouch_push_m + dip
            return dip, (envelope.dip_s(dip, g) if dip > 0 else 0.0), (S + s) / vz
        lift = self.push_off_m - self.dip_m
        dip = max(MIN_DIP_M, s - lift)
        return dip, envelope.dip_s(dip, g), (dip + lift + s) / vz


def muscle_push(v_along):
    """The legs' horizontal push per kilogram at a speed along the push (the sprint profile, zero beyond V0)."""
    return SPRINT_F0_N_KG * max(0.0, 1 - v_along / SPRINT_V0_M_S)


def _bodies():
    out = {}
    for key, (name, b) in zip(('dancer', 'larger'), envelope.bodies().items()):
        s = b['stature_m'] / envelope.STATURE_M
        out[key] = Body(key=key, name=name, mass_kg=b['mass_kg'], stature_m=b['stature_m'],
                        com_height_m=b['com_height_m'], ready_com_m=b['ready_com_m'],
                        ready_half_stance_m=b['ready_half_stance_m'], side_reach_m=b['side_reach_m'],
                        fore_aft_reach_m=b['fore_aft_reach_m'], split_reach_m=b['split_reach_m'],
                        push_off_m=b['push_off_m'], release_above_com_m=b['release_above_com_m'],
                        inertia_kg_m2=dict(b['inertia_kg_m2']), reach_up_m=round(REACH_UP_RATIO * b['stature_m'], 3),
                        reach_out_m=round(REACH_OUT_M * s, 3), reach_low_m=round(REACH_LOW_M * s, 3))
    return out


BODIES = _bodies()


@dataclass(frozen=True)
class State:
    """A player at one moment."""
    p: tuple = (0.0, 0.0)                 # the centre of mass over the floor, m
    v: tuple = (0.0, 0.0)                 # its horizontal velocity, m/s
    posture: str = 'ready'                # one of POSTURES
    facing: float = 0.0                   # the chest's heading, rad: the split stance's forward, a throw's line
    z: float | None = None                # aloft: the centre of mass's height now (None on the floor)
    vz: float = 0.0                       # aloft: its vertical speed now, up positive
    next_stance_s: float = 0.0            # moving: time until the next foot is down

    def landing_s(self, body, g):
        """Aloft: the time until the feet touch, when the centre of mass is back at the standing height."""
        if self.posture != 'aloft' or self.z is None:
            return 0.0
        dz = self.z - body.com_height_m
        return (self.vz + math.sqrt(max(self.vz * self.vz + 2 * g * dz, 0.0))) / g


@dataclass(frozen=True)
class Reach:
    """What the reach engine answers: whether, when and how the hands can be at a point."""
    ok: bool
    seconds: float                        # reach_time: the earliest time; reach_at: the time asked
    margin_m: float                       # reach_at: distance to spare at that time (negative: short by that much)
    mode: str                             # 'ground', 'jump', 'aloft' or 'landing' (aloft, then on from the landing)
    moves: tuple = ()                     # vocabulary unit ids in order
    committed_s: float = 0.0              # of the time to the point, how long the body was in an arc it could not change
    take_off_s: float | None = None       # when a jump leaves the floor: the moment an opponent can read it

    @property
    def reads_s(self):
        """How long before the point is reached an opponent sees the commitment (0 on the ground)."""
        return self.seconds - self.take_off_s if self.take_off_s is not None else 0.0
