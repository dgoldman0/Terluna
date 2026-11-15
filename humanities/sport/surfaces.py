"""Walls, ledges and bars: the surface moves of a lunar obstacle game on the Moon beside Earth, and routes across an
arena built of obstacles.

    from humanities.sport import surfaces
    from humanities.sport.state import BODIES, G_M, State
    surfaces.tic_tac(BODIES['dancer'], 3.0, 45.0, G_M)                 # a foot's push off a wall
    surfaces.quickest_route(BODIES['dancer'], State(p=surfaces.EXAMPLE_START), surfaces.EXAMPLE_GOAL,
                            surfaces.EXAMPLE, 0.8, G_M)

The sport study serves a paintball game whose arenas are built from obstacles. On the floor a lunar player's changes
of speed and direction are capped by the grip, mu g on average, about 1 m/s per stride, so a sharp cut at speed is a
visible jump (envelope.cut_bounce). A push off a wall or an obstacle is not capped that way: the push presses the foot
or the hand into the surface, so the grip grows with the push. This module times the moves that use surfaces, named by
PARKOUR_UNITS, and routes a player across an arena with them. It adds to the contract without changing state.py or
reach.py: floor travel is the reach engine's, and the moves start and end in its States.

The model and its choices. Values marked (assumed) are this screen's choices, none retrieved; the larger body scales
every length by stature. Gravity enters through each move's own clock; the design air plays no part in any of them.
1. Pushing off a surface. A foot pushes along the surface's normal over the leg's stroke on the vertical jump's
   force-velocity profile with no gravity along the push: Samozino's v^2 = 2 s F0 (1 - v/(2 V0)) at g = 0, which gives
   3.17 m/s over the dance body's 0.36 m (the vertical take-off gives 3.05 on the Moon, 2.38 on Earth). The contact
   first takes up the incoming normal speed, yielding no softer than the push and within the leg's stroke (assumed: up
   to 1.5 F0 while yielding), then pushes out at the push's mean acceleration F0 (1 - v/(2 V0)); its normal impulse
   per kilogram is the incoming plus the outgoing normal speed, and friction gives up to the surface's grip times that
   along the surface, less what holds the body up against gravity through the contact. A one-footed push may reach
   less than the two-footed profile (ONE_LEG_SHARE 0.6 of F0 as the bracket). The floor's cut through one stance
   (floor_cut) changes the velocity by at most mu (2 w + g T), so a sharp cut needs a bounce w.
2. Running up a wall (wall_run): a run meets a wall head on. The first step's normal impulse takes up the approach and
   pushes the body off by delta; each later step pushes off by delta more, from the speed away the body already has
   (the force-velocity profile at the push's mean speed) over the stroke the wall's reach leaves. delta is the leg's
   full push for one step, and for two or three steps the largest equal push that keeps the centre of mass within
   WALL_REACH_M of the wall and the legs able to push. Each step's friction lifts by the wall's grip times its normal
   impulse (delivered at the step's start, assumed), gravity acts throughout, and a step lasts STEP_S or its contact,
   whichever is longer. A flat wall throws the runner off: more steps keep it close but lift it less. A curved wall
   pushes back on a runner along it, v^2 / R, and holds the weight while mu v^2 / R >= g (wall_of_death_m_s).
3. Hands on an edge. One arm-strength assumption serves all three: a trained player's two arms pull ARM_PULL_BW = 1.5
   Earth body weights (assumed; bracket 1.0-2.0), the same on either world, over ARM_STROKE_M. cat_leap: the arms stop
   the descent within their give (ARM_GIVE_M) while v^2/(2 give) + g <= their pull; the feet take the travel into the
   wall over the leg's stroke. climb_up: the centre of mass rises from a straight-arm hang (HANG_BELOW_M under the
   hands) to crouched on the top (CLIMB_ABOVE_M over the edge), pulled then pressed by the arms over two strokes and
   coasting the rest to rest (envelope.rise_and_stop_s); where the arms cannot, a foot comes up (STEP_S) and the legs
   stand the rest. pull_down: the hands on cover pull the body down at g plus the arms' pull, then the legs brake
   (envelope.dip_s's form, against its gravity-limited sink).
4. Over and off obstacles. vault and kong_vault: a running take-off keeps the run's speed; the centre of mass must
   stay VAULT_CLEAR_M above the top over the box's length (assumed: 0.15 m for a speed vault, the hip passing beside
   the hand; 0.25 m for a kong's dive over the hands). The hands on the top hold the body up while it passes, with up
   to their pull (one arm's half of it for a speed vault, both arms for a kong), so the centre of mass falls over the
   box at g less that support; a kong's hands also push forward along the edge of the top's friction cone for
   HAND_PUSH_S. The cheapest take-off reaches the near edge at the clearance height rising just fast enough to stay
   above it to the far edge, and leaves wherever that puts it before the box. drop_roll: the landing speed
   sqrt(2 g h), the legs' braking landing (envelope.braking_landing) and a roll that spreads it over ROLL_M more.
   precision_jump: the quickest arc that lands on an edge and stops on it: the landing sheds min(grip, lean) times its
   vertical impulse, the lean being the centre of pressure's most lead over a narrow edge, (edge/2 + LAND_LEAN_M) / h,
   the feet meeting the edge ahead of the centre of mass (assumed), which must not travel past it in the stance.
   slide: the centre of mass sinks to SLIDE_COM_M no faster than g, keeping the run's speed, then slides at the
   clothing's grip on the floor (SLIDE_GRIP 0.4, assumed).
5. Bars (lache): the body hangs at arm's length as a physical pendulum, I = I_layout + m r^2 about the bar; the arms
   can pull it in by ARM_STROKE_M at the bottom where their pull exceeds the weight and the swing's centripetal need,
   which spins it faster (angular momentum kept); released RELEASE_DEG past the bottom, it flies freely and keeps
   turning at its swing's rate.
6. An arena (Arena of Box and Bar): axis-aligned boxes (bunkers and walls, with heights; their tops are levels) and
   bars. route_time times a sequence of moves from a State, each move reached by floor travel (the reach engine's
   law, around boxes by the shortest path between their corners, BODY_R_M clear), and checks each is feasible; a
   foot's push in a route aims at the next point of the way on, with the full push or only as hard as that aim needs
   (the search tries both; a move without a choice takes whichever reaches that point sooner); a wall is crossed by a
   step that only takes up the run (wall_run with push=False, staying at the wall), a catch of its top, a climb, a
   walk across and a drop and roll;
   quickest_route searches routes of up to three surface moves (a beam over the arena's candidates) against floor
   travel alone and an open floor. Floor travel is the reach engine's: from a standing posture through
   reach.reach_along (its table), and once in a gait by the 'moving' profile's closed form (the grip's parabola, the
   sprint profile's exponential, V_MAX_M_S), which reach.py integrates; a leg is straight in the drift frame, and the
   drift's curve past a corner is not checked against the boxes. A runner takes a corner no faster than the larger of
   a cut through one stance without a bounce and a steady turn on CORNER_R_M at the grip (corner_m_s), braking for it
   at min(grip g, F0): the extra time over arriving at the leg's own speed v is (v - c)^2 / (2 a v).

Source read: Minetti, Luciano, Natalucci and Pavei (2024), Horizontal running inside circular walls of Moon
settlements, Royal Society Open Science 11: 231906 (abstract and the journal page's summary figures only): runners
under 83% body-weight unloading kept horizontal running at 5.4-6.5 m/s inside a 4.7 m radius 'Wall of Death', with a
theoretical least speed of 3.5 m/s.
"""
from __future__ import annotations

import heapq
import math
from dataclasses import dataclass

import numpy as np

from humanities.dance import envelope as dance
from humanities.sport import envelope as env
from humanities.sport import reach, state
from humanities.sport.state import G_E, STEP_S, V_MAX_M_S, State

PARKOUR_UNITS = ('vault', 'kong_vault', 'tic_tac', 'wall_run', 'cat_leap', 'climb_up', 'precision_jump', 'drop_roll',
                 'lache', 'slide', 'pull_down')
WALL_GRIP = 0.8                           # a shoe or a hand on a wall or an obstacle's top (assumed)
ARM_PULL_BW = 1.5                         # a trained player's two-arm pull, Earth body weights (assumed)
ARM_PULL_BRACKET = (1.0, 2.0)
ARM_STROKE_M = 0.60                       # the arms' stroke, hands overhead to the shoulders (dance body)
ARM_GIVE_M = 0.30                         # the arms' give in a catch
HANG_BELOW_M = 1.05                       # the centre of mass below the hands in a straight-arm hang
CLIMB_ABOVE_M = 0.60                      # and above the edge, crouched on the top
HAND_PUSH_S = env.THROW_S                 # a hand's push on a box top
YIELD_SHARE = 1.5                         # the legs take up to this times F0 while yielding (assumed)
ONE_LEG_SHARE = 0.6                       # a one-footed push's force against the two-footed profile (bracket)
WALL_REACH_M = 0.90                       # the centre of mass's farthest distance from a wall a foot still pushes on
WALL_STANDOFF_M = 0.45                    # its distance from the wall as the foot meets it
VAULT_MAX_M = 1.3                         # boxes taller than this are also climbed over (wall_run ... drop_roll)
VAULT_TAKEOFF_M = 0.80                    # the take-off's distance before the box
VAULT_CLEAR_M = dict(vault=0.15, kong_vault=0.25)
VAULT_MIN_M_S = 1.0
TIC_TAC_MIN_M = 1.0                       # faces this high take a foot's push
SLIDE_GRIP = 0.4                          # clothing on the floor (assumed)
SLIDE_COM_M = 0.35                        # the centre of mass in a feet-first slide
ROLL_M = 1.0                              # the path a roll spreads a landing over, beyond the legs' stroke
ROLL_EXIT_M_S = 1.0                       # the travel a roll off a wall's top leaves (assumed)
EDGE_W_M = 0.20
LAND_LEAN_M = 0.30                        # the landing feet's lead on the centre of mass in a precision landing
BODY_R_M = 0.30                           # a body's clearance from boxes on the floor plan
RELEASE_DEG = 30.0
CROSS_M_S = 1.0                           # walking across a wall's top
CORNER_R_M = 1.0                          # the curve a runner may take round a corner (assumed)
SAMPLE_S = 0.005                          # the floor legs' time grid


def _scale(body):
    return body.stature_m / env.STATURE_M


# --- 1. pushing off a surface ----------------------------------------------------------------------------------------

def push_speed(body, f0_share=1.0, stroke_m=None):
    """The speed a push along a surface's normal gives from rest over the leg's stroke, on the jump's force-velocity
    profile with no gravity along the push (dance.samozino_take_off at g = 0)."""
    stroke = body.push_off_m if stroke_m is None else stroke_m
    return dance.samozino_take_off(0.0, f0=f0_share * dance.F0_N_KG, h=stroke)


def push_from(body, u, stroke_m, f0_share=1.0):
    """The speed a push adds over stroke_m to a body already moving away from the surface at u: Samozino's mean force at
    the mean speed, (u + d)^2 - u^2 = 2 s F0 (1 - (u + d/2)/V0), no gravity along it (0 once u reaches V0)."""
    F0, V0 = f0_share * dance.F0_N_KG, dance.V0_M_S
    if u >= V0 or stroke_m <= 0:
        return 0.0
    b = 2 * u + stroke_m * F0 / V0
    return (-b + math.sqrt(b * b + 8 * stroke_m * F0 * (1 - u / V0))) / 2


def contact_s(body, v_in, v_out, f0_share=1.0):
    """A contact that takes up the normal speed v_in and pushes out at v_out: the push's mean acceleration over the
    stroke v_out needs at full effort is F0 (1 - v_out/(2 V0)); the yield is no softer than that and fits within the
    leg's stroke. Returns the contact's seconds and the yield's acceleration."""
    a_push = f0_share * dance.F0_N_KG * (1 - v_out / (2 * dance.V0_M_S))
    a_take = max(a_push, v_in * v_in / (2 * body.push_off_m))
    return (v_in / a_take if v_in > 0 else 0.0) + v_out / a_push, a_take


def floor_cut(speed, turn_deg, grip, g):
    """The floor's cut through one stance (envelope.cut_bounce, unrounded): the bounce w each way that the grip on the
    stance's vertical impulse needs, mu (2 w + g T) = 2 v sin(turn/2), its tell and its commitment (each 2 w / g)."""
    dv = 2 * speed * math.sin(math.radians(turn_deg) / 2)
    w = max(0.0, (dv / grip - g * env.STRIDE_S) / 2)
    return dict(turn_deg=turn_deg, velocity_change_m_s=dv, bounce_m_s=w, bounce_rise_m=w * w / (2 * g),
                tell_s=2 * w / g, commits_s=2 * w / g, total_s=4 * w / g + env.STRIDE_S,
                one_stance_m_s=grip * g * env.STRIDE_S)


def _turn_deg(vt, vn, vt2, vp):
    """The angle between the travel in (vt along the wall, vn into it) and out (vt2 along, vp away)."""
    c = (vt * vt2 - vn * vp) / (math.hypot(vt, vn) * math.hypot(vt2, vp))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def tic_tac(body, speed, angle_deg, g, wall_grip=WALL_GRIP, grip=state.MU, keep_height=True, f0_share=1.0):
    """A foot's push off a wall met at `speed`, travelling angle_deg from the wall's normal (0: head on). The normal
    impulse is the incoming normal speed plus the push's; friction gives up to the wall's grip times that along the
    wall, less the g T that holds the body's height through the contact (where the cone cannot hold even that, the
    body sinks by the rest, sinks_m_s). 'kept': the along-wall speed kept, the push
    turning the travel; 'most': the along-wall speed cut as far as friction allows, the largest turn. Beside it, the
    floor's cut through the same turn at the same speed."""
    a = math.radians(angle_deg)
    vn, vt = speed * math.cos(a), speed * math.sin(a)
    vp = push_speed(body, f0_share)
    J = vn + vp
    T, a_take = contact_s(body, vn, vp, f0_share)
    hold = g * T if keep_height else 0.0
    side = math.sqrt(max((wall_grip * J) ** 2 - hold * hold, 0.0))
    most = vt - side
    kept = dict(speed_out_m_s=math.hypot(vt, vp), turn_deg=_turn_deg(vt, vn, vt, vp), velocity_change_m_s=J)
    turned = dict(speed_out_m_s=math.hypot(most, vp), turn_deg=_turn_deg(vt, vn, most, vp),
                  velocity_change_m_s=math.hypot(J, side))
    return dict(unit='tic_tac', speed_m_s=speed, angle_deg=angle_deg, push_out_m_s=vp, normal_impulse_m_s=J,
                contact_s=T, friction_m_s=wall_grip * J, holds_height_m_s=hold, along_wall_m_s=side,
                sinks_m_s=max(hold - wall_grip * J, 0.0), kept=kept,
                most=turned, yield_m_s2=a_take, yield_ok=a_take <= YIELD_SHARE * f0_share * dance.F0_N_KG,
                floor_cut=floor_cut(speed, kept['turn_deg'], grip, g))


# --- 2. running up a wall --------------------------------------------------------------------------------------------

def wall_run(body, approach, g, wall_grip=WALL_GRIP, steps=1, step_s=STEP_S, push=True):
    """Up a wall met head on at `approach` m/s in 1-3 steps (see the module's note 2): the centre of mass's rise to its
    apex, the hands' top there, the time on the wall, the vertical and outward speeds off it. push=False: the step
    only takes up the approach and stays at the wall (one step: the way up to a ledge)."""
    if not push:
        steps = 1
    vp = push_speed(body) if push else 0.0
    reach_m = WALL_REACH_M * _scale(body)

    def possible(delta):
        for i in range(1, steps):             # step i+1 starts this far from the wall, moving away at i delta
            off = delta * step_s * i * (i + 1) / 2
            if off > reach_m or delta > push_from(body, i * delta, min(body.push_off_m, reach_m - off)) + 1e-12:
                return False
        return True

    delta = vp
    if push and steps > 1:
        lo, hi = 0.0, vp
        for _ in range(50):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if possible(mid) else (lo, mid)
        delta = lo
    first, a_take = contact_s(body, approach, delta) if push else (2 * body.push_off_m / max(approach, 0.1),
                                                                    approach * approach / (2 * body.push_off_m))
    vz = z = t = 0.0
    lifts = []
    for i in range(steps):
        J = (approach if i == 0 else 0.0) + delta
        vz += wall_grip * J
        lifts.append(wall_grip * J)
        dt = max(step_s, first) if i == 0 else step_s
        z += vz * dt - 0.5 * g * dt * dt
        vz -= g * dt
        t += dt
    rise = z + (vz * vz / (2 * g) if vz > 0 else 0.0)
    return dict(unit='wall_run', approach_m_s=approach, steps=steps, push_each_m_s=delta, lifts_m_s=lifts,
                time_on_wall_s=t, rise_on_wall_m=z, vertical_off_m_s=vz, outward_off_m_s=steps * delta, rise_m=rise,
                apex_s=t + max(vz, 0.0) / g, hands_top_m=body.reach_up_m + rise, yield_m_s2=a_take,
                yield_ok=a_take <= YIELD_SHARE * dance.F0_N_KG)


def wall_of_death_m_s(g, radius_m, wall_grip=WALL_GRIP):
    """The least speed that keeps a runner on a curved wall, the wall's push v^2/R holding the weight by friction."""
    return math.sqrt(g * radius_m / wall_grip)


# --- 3. hands on an edge ---------------------------------------------------------------------------------------------

def cat_leap(body, g, v_down=0.0, v_in=0.0, arm_bw=ARM_PULL_BW):
    """The hands catch an edge, the body coming down at v_down and travelling into the wall at v_in. The arms hold
    while v_down^2/(2 give) + g <= their pull; the feet take v_in against the wall over the leg's stroke."""
    give = ARM_GIVE_M * _scale(body)
    pull = arm_bw * G_E
    need = v_down * v_down / (2 * give) + g
    fastest = math.sqrt(2 * give * (pull - g)) if pull > g else 0.0
    legs = v_in * v_in / (2 * body.push_off_m)
    return dict(unit='cat_leap', arms_need_bw=need / G_E, arms_bw=arm_bw, holds=need <= pull + 1e-9,
                fastest_catch_m_s=fastest, drop_above_catch_m=fastest * fastest / (2 * g), legs_m_s2=legs,
                legs_ok=legs <= YIELD_SHARE * dance.F0_N_KG, catch_s=2 * give / v_down if v_down > 0 else 0.0)


def climb_up(body, g, arm_bw=ARM_PULL_BW):
    """From a straight-arm hang to crouched on the top (see the module's note 3)."""
    s = _scale(body)
    H = (HANG_BELOW_M + CLIMB_ABOVE_M) * s
    a = arm_bw * G_E - g
    stroke = 2 * ARM_STROKE_M * s
    if a <= 0:
        return dict(unit='climb_up', ok=False, seconds=math.inf, rise_m=H, how='the arms cannot lift the body')
    if H * g / (a + g) <= stroke:
        return dict(unit='climb_up', ok=True, seconds=env.rise_and_stop_s(H, a, g), rise_m=H, how='arms', arms_m_s2=a)
    v = body.take_off_m_s(g)
    legs = v * v / (2 * body.push_off_m)
    t = env.rise_and_stop_s(stroke, a, g) + STEP_S + env.rise_and_stop_s(H - stroke, legs, g)
    return dict(unit='climb_up', ok=True, seconds=t, rise_m=H, how='arms, then legs', arms_m_s2=a)


def pull_down(body, depth_m, g, arm_bw=ARM_PULL_BW, brake=env.BRAKE_M_S2):
    """Sinking depth_m behind cover with the hands pulling down: the fall at g plus the arms' pull, then the legs'
    brake (envelope.dip_s's form), against envelope.dip_s's sink at g alone."""
    a = g + arm_bw * G_E
    d1 = depth_m * brake / (a + brake)
    t = math.sqrt(2 * d1 / a) + math.sqrt(2 * (depth_m - d1) / brake)
    alone = env.dip_s(depth_m, g)
    return dict(unit='pull_down', depth_m=depth_m, seconds=t, gravity_alone_s=alone, saves_s=alone - t)


# --- 4. over and off obstacles ---------------------------------------------------------------------------------------

def vault(body, height_m, length_m, speed, g, kind='vault', wall_grip=WALL_GRIP, arm_bw=ARM_PULL_BW):
    """Over a box height_m high and length_m long at a run of `speed` (see the module's note 4). Returns the flight
    from the take-off to the landing at the standing height, the speed kept, and what the hands add."""
    z_to = reach.take_off_height(body)
    top = height_m + VAULT_CLEAR_M[kind] * _scale(body)
    F = (0.5 if kind == 'vault' else 1.0) * arm_bw * G_E
    beta = math.atan(wall_grip) if kind == 'kong_vault' else 0.0
    hold = min(F * math.cos(beta), g)
    v = max(speed, 1e-6)

    def arc(g_box):                           # (take-off vz, rise at the near edge) for a fall of g_box over the box
        w = 0.5 * g_box * length_m / v
        return math.sqrt(w * w + 2 * g * max(top - z_to, 0.0)), w

    vz, w = arc(g - hold)
    bare, _ = arc(g)
    t1, t2 = max(vz - w, 0.0) / g, length_m / v
    t3 = (-w + math.sqrt(w * w + 2 * g * (top - body.com_height_m))) / g
    dvx = F * math.sin(beta) * HAND_PUSH_S
    T = t1 + t2 + t3
    x_land = v * T + dvx * (t2 + t3)
    return dict(unit=kind, ok=vz <= body.take_off_m_s(g), height_m=height_m, length_m=length_m, speed_m_s=speed,
                flight_s=T, take_off_before_m=v * t1, lands_from_take_off_m=x_land,
                beyond_box_m=x_land - v * t1 - length_m,
                speed_kept_m_s=speed + dvx, take_off_m_s=vz, without_hands_m_s=bare, hands_hold_share=hold / g,
                hands_forward_m_s=dvx)


def drop_roll(body, level_m, g, speed=0.0, grip=state.MU):
    """A drop from a level (the centre of mass from standing on it to standing below): the landing speed, the legs'
    braking landing over envelope.STROKE_M and a roll that spreads it over ROLL_M more, as Earth body weights."""
    v = math.sqrt(2 * g * level_m)
    roll = ROLL_M * _scale(body)
    legs_bw = (v * v / (2 * env.STROKE_M) + g) / G_E
    roll_bw = (v * v / (2 * (env.STROKE_M + roll)) + g) / G_E
    out = dict(unit='drop_roll', level_m=level_m, fall_s=math.sqrt(2 * level_m / g), landing_m_s=v,
               legs_push_earth_bw=legs_bw, roll_push_earth_bw=roll_bw, roll_s=2 * (env.STROKE_M + roll) / v,
               like_an_earth_drop_m=level_m * g / G_E)
    if speed > 0:
        b = env.braking_landing(speed, v, grip, g)
        out['braking_landing'] = dict(shed_m_s=b['shed_in_landing_m_s'], left_m_s=b['left_m_s'], total_m=b['total_m'])
    return out


def precision_jump(body, gap_m, g, rise_m=0.0, edge_w=EDGE_W_M, edge_grip=WALL_GRIP, n=120):
    """A jump that carries the centre of mass gap_m to an edge rise_m above the take-off's floor and stops on it: the
    quickest arc the take-off allows whose landing sheds its travel (see the module's note 4)."""
    z_to = reach.take_off_height(body)
    z_land = rise_m + body.com_height_m
    room = edge_w / 2 + LAND_LEAN_M * _scale(body)
    lean = room / body.com_height_m
    k = min(edge_grip, lean)
    full = body.take_off_m_s(g)
    best, widest = None, 0.0
    for vz in np.linspace(full / n, full, n):
        disc = vz * vz - 2 * g * (z_land - z_to)
        if disc < 0:
            continue
        vd = math.sqrt(disc)
        T = (vz + vd) / g
        stance = 2 * env.STROKE_M / max(vd, 0.05)
        shed = k * (vd + g * stance)
        vx = gap_m / T
        widest = max(widest, T * min(shed, 2 * room / stance))
        if shed >= vx and vx * stance / 2 <= room and (best is None or T < best['flight_s']):
            best = dict(flight_s=T, take_off_m_s=float(vz), horizontal_m_s=vx, landing_m_s=vd, sheds_m_s=shed,
                        stance_s=stance)
    out = dict(unit='precision_jump', gap_m=gap_m, rise_m=rise_m, edge_w_m=edge_w, ok=best is not None,
               limited_by='balance' if lean < edge_grip else 'grip', shed_share=k, widest_gap_m=widest)
    if best:
        out.update(best)
    return out


def slide(body, speed, g, slide_grip=SLIDE_GRIP, grip=state.MU, cover_m=None):
    """From a run at `speed`, a feet-first slide into cover (see the module's note 4), against running into the cover,
    braking at min(grip g, F0) and then sinking to the same height (envelope.dip_s). cover_m: the cover's distance ahead
    when the choice is made (default the slide's own reach). Returns when each is low and when each is stopped."""
    s = _scale(body)
    sink = body.com_height_m - SLIDE_COM_M * s
    t_sink = math.sqrt(2 * sink / g)
    v1 = max(speed - slide_grip * math.sqrt(2 * g * sink), 0.0)
    slide_m, slide_s = v1 * v1 / (2 * slide_grip * g), v1 / (slide_grip * g)
    reach_m = speed * t_sink + slide_m
    a = min(grip * g, state.SPRINT_F0_N_KG)
    brake_m, brake_s = speed * speed / (2 * a), speed / a
    D = max(reach_m, brake_m) if cover_m is None else cover_m
    run_stopped = (D - brake_m) / speed + brake_s
    out = dict(unit='slide', speed_m_s=speed, sink_s=t_sink, slide_m=slide_m, slide_s=slide_s, reach_m=reach_m,
               cover_m=D, run=dict(low_s=run_stopped + env.dip_s(sink, g), stopped_s=run_stopped))
    if D >= reach_m:
        lead = (D - reach_m) / speed
        out['slide'] = dict(low_s=lead + t_sink, stopped_s=lead + t_sink + slide_s)
    return out


# --- 5. bars ----------------------------------------------------------------------------------------------------------

def lache(body, g, amplitude_deg=None, catch_m_s=None, release_deg=RELEASE_DEG, bar_m=2.0, arm_bw=ARM_PULL_BW,
          pull=True):
    """A swing on a bar and a release (see the module's note 5). The swing starts from rest at amplitude_deg or from the
    bottom at catch_m_s (a runner catching the bar). Returns the period, the speeds, the swing's time from the bottom to
    the release, and the flight to the floor."""
    s = _scale(body)
    m, r, Icm = body.mass_kg, HANG_BELOW_M * s, body.inertia_kg_m2['layout']
    I = Icm + m * r * r
    period = 2 * math.pi * math.sqrt(I / (m * g * r))
    if catch_m_s is not None:
        w0 = catch_m_s / r
    else:
        w0 = math.sqrt(2 * m * g * r * (1 - math.cos(math.radians(amplitude_deg))) / I)
    r1, I1, w1, pulled = r, I, w0, False
    if pull:
        rp = r - ARM_STROKE_M * s
        Ip = Icm + m * rp * rp
        wp = w0 * I / Ip
        if m * (g + wp * wp * rp) <= arm_bw * G_E * m:
            r1, I1, w1, pulled = rp, Ip, wp, True
    work = 0.5 * I1 * w1 * w1 - 0.5 * I * w0 * w0 + m * g * (r - r1)
    th = math.radians(release_deg)
    wr2 = w1 * w1 - 2 * m * g * r1 * (1 - math.cos(th)) / I1
    out = dict(unit='lache', period_s=period, radius_m=r, bottom_m_s=w0 * r, pulled_in=pulled, arms_work_J=work,
               after_pull_m_s=w1 * r1, release_deg=release_deg)
    if wr2 <= 0:
        out.update(ok=False)
        return out
    ths = np.linspace(0.0, th, 401)
    w = np.sqrt(np.maximum(w1 * w1 - 2 * m * g * r1 * (1 - np.cos(ths)) / I1, 1e-12))
    swing = float(np.sum(0.5 * (1 / w[1:] + 1 / w[:-1]) * np.diff(ths)))
    wr = math.sqrt(wr2)
    vr = wr * r1
    vx, vz = vr * math.cos(th), vr * math.sin(th)
    zr = bar_m - r1 * math.cos(th)
    T = (vz + math.sqrt(vz * vz + 2 * g * max(zr - body.com_height_m, 0.0))) / g
    out.update(ok=True, release_m_s=vr, swing_s=swing, release_height_m=zr, flight_s=T,
               lands_beyond_bar_m=r1 * math.sin(th) + vx * T, horizontal_m_s=vx,
               spin_in_flight_deg=math.degrees(wr * T))
    return out


# --- 6. an arena and its routes --------------------------------------------------------------------------------------

@dataclass(frozen=True)
class Box:
    """An axis-aligned box on the floor plan: a bunker or a wall, its top a level."""
    x0: float
    y0: float
    x1: float
    y1: float
    height: float

    @property
    def centre(self):
        return (0.5 * (self.x0 + self.x1), 0.5 * (self.y0 + self.y1))

    def face(self, name):
        """A face's two ends and its outward normal: 'x0' (at x0, facing -x), 'x1', 'y0', 'y1'."""
        return {'x0': ((self.x0, self.y0), (self.x0, self.y1), (-1.0, 0.0)),
                'x1': ((self.x1, self.y0), (self.x1, self.y1), (1.0, 0.0)),
                'y0': ((self.x0, self.y0), (self.x1, self.y0), (0.0, -1.0)),
                'y1': ((self.x0, self.y1), (self.x1, self.y1), (0.0, 1.0))}[name]


@dataclass(frozen=True)
class Bar:
    """A bar between two points on the floor plan at a height."""
    x0: float
    y0: float
    x1: float
    y1: float
    height: float = 2.0


@dataclass(frozen=True)
class Arena:
    boxes: tuple = ()
    bars: tuple = ()


EXAMPLE = Arena(boxes=(Box(1.8, -1.0, 2.8, 1.0, 1.2),      # the start's bunker
                       Box(6.0, -0.5, 7.5, 1.5, 1.0),      # a low box
                       Box(10.0, -2.0, 10.3, 4.0, 3.0),    # a wall
                       Box(15.0, 1.5, 16.0, 3.5, 1.4)),    # the goal's bunker
                bars=(Bar(12.0, -3.0, 12.0, 0.0, 2.0),))
EXAMPLE_START = (1.0, 0.0)
EXAMPLE_GOAL = (14.2, 2.5)
BARRIER = Arena(boxes=(Box(8.0, -6.0, 8.6, 6.0, 1.0),))    # a low barrier across the way, open past its ends
BARRIER_START = (1.0, 0.0)
BARRIER_GOAL = (15.0, 0.0)


def _hits(p, q, b, r):
    """Whether the segment p-q passes through the box grown by r (touching its edge does not count)."""
    lo, hi = (b.x0 - r, b.y0 - r), (b.x1 + r, b.y1 + r)
    t0, t1 = 0.0, 1.0
    for k in (0, 1):
        d = q[k] - p[k]
        if abs(d) < 1e-12:
            if p[k] <= lo[k] + 1e-9 or p[k] >= hi[k] - 1e-9:
                return False
            continue
        ta, tb = (lo[k] - p[k]) / d, (hi[k] - p[k]) / d
        t0, t1 = max(t0, min(ta, tb)), min(t1, max(ta, tb))
        if t0 >= t1 - 1e-9:
            return False
    return True


def clear(arena, p, q, r=BODY_R_M):
    return not any(_hits(p, q, b, r) for b in arena.boxes)


def inside(arena, p, r=BODY_R_M):
    return any(b.x0 - r < p[0] < b.x1 + r and b.y0 - r < p[1] < b.y1 + r for b in arena.boxes)


def floor_path(arena, p, q, r=BODY_R_M):
    """The shortest way on the floor plan from p to q around the boxes (grown by r), through their corners: the points
    after p, or None."""
    p, q = (float(p[0]), float(p[1])), (float(q[0]), float(q[1]))
    if clear(arena, p, q, r):
        return [q]
    e = r + 0.02
    nodes = [p, q] + [c for b in arena.boxes for c in ((b.x0 - e, b.y0 - e), (b.x1 + e, b.y0 - e),
                                                         (b.x0 - e, b.y1 + e), (b.x1 + e, b.y1 + e))
                       if not inside(arena, c, r)]
    dist, prev, todo = {0: 0.0}, {}, [(0.0, 0)]
    while todo:
        d, i = heapq.heappop(todo)
        if i == 1:
            break
        if d > dist.get(i, math.inf):
            continue
        for j in range(1, len(nodes)):
            if j != i and clear(arena, nodes[i], nodes[j], r):
                nd = d + math.dist(nodes[i], nodes[j])
                if nd < dist.get(j, math.inf):
                    dist[j], prev[j] = nd, i
                    heapq.heappush(todo, (nd, j))
    if 1 not in dist:
        return None
    path, i = [], 1
    while i != 0:
        path.append(nodes[i])
        i = prev[i]
    return path[::-1]


def moving_x(t, va, grip, g):
    """The reach engine's 'moving' profile in closed form: x'' = min(mu g, the sprint profile at va + x'), the speed
    along u held at or below V_MAX_M_S - va (reach.py integrates the same law). Returns x and x' (arrays)."""
    t, va = np.broadcast_arrays(np.asarray(t, float), np.asarray(va, float))
    F0, V0 = state.SPRINT_F0_N_KG, state.SPRINT_V0_M_S
    r, a = F0 / V0, grip * g
    cap = np.maximum(V_MAX_M_S - va, 0.0)
    s = np.maximum(t, 0.0)
    u1 = np.clip(V0 * (1 - a / F0) - va, 0.0, cap)
    t1, S0 = u1 / a, va + u1
    capped_by_grip = u1 >= cap
    s2c = np.where(capped_by_grip, 0.0, np.log(np.maximum(V0 - S0, 1e-12) / (V0 - V_MAX_M_S)) / r)

    def muscles(s2):
        e = np.exp(-r * s2)
        return 0.5 * a * t1 * t1 + (V0 - va) * s2 - (V0 - S0) / r * (1 - e), V0 - (V0 - S0) * e - va

    tc = np.where(capped_by_grip, cap / a, t1 + s2c)
    xc = np.where(capped_by_grip, 0.5 * cap * cap / a, muscles(s2c)[0])
    xm, um = muscles(np.maximum(s - t1, 0.0))
    x, u = np.where(s <= t1, 0.5 * a * s * s, xm), np.where(s <= t1, a * s, um)
    return np.where(s > tc, xc + cap * (s - tc), x), np.where(s > tc, cap, u)


def floor_leg(body, st, target, grip, g, arrive='com', t_max=25.0):
    """The quickest floor travel from st to a target on the floor plan, straight in the drift frame: the earliest time
    the centre of mass (arrive='com') or the hands (arrive='hands', reach_out short) can be there on the ground, and
    the State then (the centre of mass on the ray, its velocity v + x' u, 'moving'). From a gait, the closed form;
    from a standing posture, the reach engine (reach.reach_along, its ground mode)."""
    R = body.reach_out_m if arrive == 'hands' else 0.0
    p, v = np.array(st.p, float), np.array(st.v, float)
    q = np.array(target, float)
    if np.hypot(*(q - p)) <= R + 1e-9 and not v.any():
        return 0.0, State(p=tuple(p), v=(0.0, 0.0), posture=st.posture, facing=st.facing)
    if st.posture == 'moving':
        t = np.arange(1, int(t_max / SAMPLE_S) + 1) * SAMPLE_S
        rel = q[None, :] - p[None, :] - v[None, :] * t[:, None]
        d = np.hypot(rel[:, 0], rel[:, 1])
        u = rel / np.maximum(d, 1e-12)[:, None]
        x, _ = moving_x(t - st.next_stance_s, u @ v, grip, g)
        f = x - (d - R)
    else:
        if st.posture not in ('upright', 'ready', 'split', 'down'):
            raise ValueError(f'floor legs start on the floor or in a gait, not {st.posture!r}')
        t = np.arange(0, int(t_max / 0.01) + 1) * 0.01
        low = body.hand_band(body.com_m('ready' if st.posture == 'down' else st.posture))[0]
        z = 0.5 * (low + body.reach_up_m)
        _, m, code = reach.reach_along(body, st, t, np.tile([q[0], q[1], z], (t.size, 1)), grip, g)
        f = np.where(code == 0, m - (body.reach_out_m - R), -1.0)
    hit = np.flatnonzero(f >= 0)
    if hit.size == 0:
        return math.inf, None
    i = int(hit[0])
    ts = float(t[i]) if i == 0 else float(t[i - 1] + (t[i] - t[i - 1]) * f[i - 1] / (f[i - 1] - f[i]))
    rel = q - p - v * ts
    dd = float(np.hypot(*rel))
    u = rel / dd if dd > 1e-12 else np.array([math.cos(st.facing), math.sin(st.facing)])
    va = float(u @ v)
    if st.posture == 'moving':
        x, xd = (float(a[()]) for a in moving_x(ts - st.next_stance_s, va, grip, g))
    else:
        _, xs, xds = reach.profile(body, st.posture, math.atan2(u[1], u[0]), grip, g, v_along=va, facing=st.facing,
                                   t=np.array([ts]), next_stance_s=st.next_stance_s)
        x, xd = float(xs[0]), float(xds[0])
    end = p + v * ts + x * u
    vel = v + xd * u
    return ts, State(p=(float(end[0]), float(end[1])), v=(float(vel[0]), float(vel[1])), posture='moving',
                     facing=math.atan2(vel[1], vel[0]) if np.hypot(*vel) > 1e-9 else st.facing)


def corner_m_s(turn_deg, grip, g):
    """The fastest a runner takes a corner turning turn_deg: the larger of a cut through one stance with no bounce,
    2 c sin(turn/2) <= mu g T, and a steady turn of radius CORNER_R_M at the grip."""
    half = math.sin(math.radians(turn_deg) / 2)
    cut = grip * g * env.STRIDE_S / (2 * half) if half > 1e-9 else math.inf
    return min(V_MAX_M_S, max(cut, math.sqrt(grip * g * CORNER_R_M)))


def _travel(body, st, target, arena, grip, g, arrive='com'):
    """Floor travel to a target around the arena's boxes, slowing for each corner: the seconds, the State at the end
    and the legs."""
    path = floor_path(arena, st.p, target)
    if path is None:
        return None
    total, legs, here = 0.0, [], st.p
    a_brake = min(grip * g, state.SPRINT_F0_N_KG)
    for i, w in enumerate(path):
        if math.dist(st.p, w) < 1e-6:
            continue
        last = i == len(path) - 1
        t, st = floor_leg(body, st, w, grip, g, arrive if last else 'com')
        if not math.isfinite(t):
            return None
        if not last:
            u_in, u_out = _unit((w[0] - here[0], w[1] - here[1])), _unit((path[i + 1][0] - w[0], path[i + 1][1] - w[1]))
            cos_turn = max(-1.0, min(1.0, u_in[0] * u_out[0] + u_in[1] * u_out[1]))
            c = corner_m_s(math.degrees(math.acos(cos_turn)), grip, g)
            v = math.hypot(*st.v)
            if v > c:
                t += (v - c) ** 2 / (2 * a_brake * v)
                st = State(p=st.p, v=(st.v[0] * c / v, st.v[1] * c / v), posture='moving', facing=st.facing)
        here = w
        total += t
        legs.append(dict(unit='floor', to=(round(w[0], 2), round(w[1], 2)), seconds=t))
    return total, st, legs


def _unit(v):
    n = math.hypot(*v)
    return (v[0] / n, v[1] / n) if n > 1e-12 else (0.0, 0.0)


def _step(body, st, mv, arena, grip, g, wall_grip=WALL_GRIP, arm_bw=ARM_PULL_BW):
    """One move of a route from st: (seconds, State after, records) or None where it cannot be done."""
    unit = mv['unit']
    if unit == 'floor':
        return _travel(body, st, mv['to'], arena, grip, g, mv.get('arrive', 'com'))
    s = _scale(body)
    if unit in ('vault', 'kong_vault'):
        b = arena.boxes[mv['box']]
        dx, dy = mv['dir']
        half = 0.5 * ((b.x1 - b.x0) if dx else (b.y1 - b.y0))
        cx, cy = b.centre
        d = VAULT_TAKEOFF_M * s
        for _ in range(3):                    # the take-off's distance follows the speed the run brings to it
            to = (cx - dx * (half + d), cy - dy * (half + d))
            if inside(arena, to, 0.0):
                return None
            trav = _travel(body, st, to, arena, grip, g)
            if trav is None:
                return None
            t0, st1, legs = trav
            v = st1.v[0] * dx + st1.v[1] * dy
            if v < VAULT_MIN_M_S:
                return None
            out = vault(body, b.height, 2 * half, v, g, unit, wall_grip, arm_bw)
            if abs(out['take_off_before_m'] - d) < 0.05:
                break
            d = out['take_off_before_m']
        side = (st1.v[0] - v * dx, st1.v[1] - v * dy)
        land = (to[0] + dx * out['lands_from_take_off_m'] + side[0] * out['flight_s'],
                to[1] + dy * out['lands_from_take_off_m'] + side[1] * out['flight_s'])
        if not out['ok'] or inside(arena, land):
            return None
        vel = (dx * out['speed_kept_m_s'] + side[0], dy * out['speed_kept_m_s'] + side[1])
        return (t0 + out['flight_s'], State(p=land, v=vel, posture='moving'),
                legs + [dict(unit=unit, box=mv['box'], seconds=out['flight_s'], speed_m_s=v)])
    if unit == 'tic_tac':
        b = arena.boxes[mv['box']]
        (ax, ay), (bx, by), (nx, ny) = b.face(mv['face'])
        if b.height < TIC_TAC_MIN_M:
            return None
        at = mv.get('at', 0.5)
        C = (ax + at * (bx - ax) + nx * WALL_STANDOFF_M * s, ay + at * (by - ay) + ny * WALL_STANDOFF_M * s)
        trav = _travel(body, st, C, arena, grip, g)
        if trav is None:
            return None
        t0, st1, legs = trav
        tx, ty = _unit((bx - ax, by - ay))
        vn = max(-(st1.v[0] * nx + st1.v[1] * ny), 0.0)
        vt = st1.v[0] * tx + st1.v[1] * ty
        way = floor_path(arena, C, mv['aim'])
        nxt = way[0] if way else mv['aim']
        ux, uy = _unit((nxt[0] - C[0], nxt[1] - C[1]))
        an, an_t = ux * nx + uy * ny, ux * tx + uy * ty
        if an < 0.05:
            return None
        full = push_speed(body)               # the full push, or only as hard as aiming at the next point needs
        pushes = [full]
        if abs(an_t) > wall_grip * an:
            pushes.append(min(full, max(an * (abs(vt) + wall_grip * vn) / (abs(an_t) - wall_grip * an), 0.05)))
        if mv.get('push') == 'full':
            pushes = pushes[:1]
        elif mv.get('push') == 'aim':
            pushes = pushes[-1:]
        best = None
        for vp in pushes:                     # whichever reaches the next point sooner
            T, _ = contact_s(body, vn, vp)
            side = math.sqrt(max((wall_grip * (vn + vp)) ** 2 - (g * T) ** 2, 0.0))
            vt2 = min(max(vp * an_t / an, vt - side), vt + side)
            vel = (vt2 * tx + vp * nx, vt2 * ty + vp * ny)
            ahead = T + floor_leg(body, State(p=C, v=vel, posture='moving'), nxt, grip, g)[0]
            if best is None or ahead < best[0]:
                best = (ahead, T, vt2, vp, vel)
        _, T, vt2, vp, vel = best
        return (t0 + T, State(p=C, v=vel, posture='moving'),
                legs + [dict(unit='tic_tac', box=mv['box'], face=mv['face'], seconds=T,
                             turn_deg=_turn_deg(vt, vn, vt2, vp) if math.hypot(vt, vn) > 1e-9 else 0.0,
                             speed_out_m_s=math.hypot(*vel))])
    if unit == 'over':                    # a wall crossed: wall_run, cat_leap, climb_up, drop_roll
        b = arena.boxes[mv['box']]
        (ax, ay), (bx, by), (nx, ny) = b.face(mv['face'])
        far = {'x0': 'x1', 'x1': 'x0', 'y0': 'y1', 'y1': 'y0'}[mv['face']]
        (fx0, fy0), (fx1, fy1), (fnx, fny) = b.face(far)
        C = (0.5 * (ax + bx) + nx * WALL_STANDOFF_M * s, 0.5 * (ay + by) + ny * WALL_STANDOFF_M * s)
        trav = _travel(body, st, C, arena, grip, g)
        if trav is None:
            return None
        t0, st1, legs = trav
        vn = max(-(st1.v[0] * nx + st1.v[1] * ny), 0.0)
        run = wall_run(body, vn, g, wall_grip, push=False)    # up the wall, staying at it, to catch the top
        reach_dz = b.height - body.reach_up_m
        z, vz = run['rise_on_wall_m'], run['vertical_off_m_s']
        if run['hands_top_m'] < b.height or WALL_STANDOFF_M * s > body.reach_out_m:
            return None
        to_catch = 0.0 if z >= reach_dz else (vz - math.sqrt(max(vz * vz - 2 * g * (reach_dz - z), 0.0))) / g
        catch = cat_leap(body, g, 0.0, vn, arm_bw)
        climb = climb_up(body, g, arm_bw)
        if not (catch['holds'] and climb['ok']):
            return None
        thick = abs((b.x1 - b.x0) * nx) + abs((b.y1 - b.y0) * ny)
        drop = drop_roll(body, b.height, g)
        T = run['time_on_wall_s'] + to_catch + climb['seconds'] + max(STEP_S, thick / CROSS_M_S) + drop['fall_s'] \
            + drop['roll_s']
        off = (WALL_STANDOFF_M + ROLL_M) * s
        end = (0.5 * (fx0 + fx1) + fnx * off, 0.5 * (fy0 + fy1) + fny * off)
        if inside(arena, end):
            return None
        return (t0 + T, State(p=end, v=(fnx * ROLL_EXIT_M_S, fny * ROLL_EXIT_M_S), posture='moving'),
                legs + [dict(unit=u, box=mv['box'], seconds=x) for u, x in
                        (('wall_run', run['time_on_wall_s'] + to_catch), ('cat_leap', 0.0),
                         ('climb_up', climb['seconds'] + max(STEP_S, thick / CROSS_M_S)),
                         ('drop_roll', drop['fall_s'] + drop['roll_s']))])
    if unit == 'lache':
        bar = arena.bars[mv['bar']]
        if bar.height > body.reach_up_m + 0.05:
            return None
        tx, ty = _unit((bar.x1 - bar.x0, bar.y1 - bar.y0))
        nx, ny = -ty * mv['sign'], tx * mv['sign']
        M = (0.5 * (bar.x0 + bar.x1), 0.5 * (bar.y0 + bar.y1))
        trav = _travel(body, st, M, arena, grip, g)
        if trav is None:
            return None
        t0, st1, legs = trav
        vn = st1.v[0] * nx + st1.v[1] * ny
        if vn < VAULT_MIN_M_S:
            return None
        out = lache(body, g, catch_m_s=vn, bar_m=bar.height, arm_bw=arm_bw)
        if not out['ok']:
            return None
        vt = st1.v[0] * tx + st1.v[1] * ty
        T = out['swing_s'] + out['flight_s']
        beyond = out['lands_beyond_bar_m']
        land = (M[0] + nx * beyond + tx * vt * T, M[1] + ny * beyond + ty * vt * T)
        if inside(arena, land):
            return None
        vel = (nx * out['horizontal_m_s'] + tx * vt, ny * out['horizontal_m_s'] + ty * vt)
        return (t0 + T, State(p=land, v=vel, posture='moving'),
                legs + [dict(unit='lache', bar=mv['bar'], seconds=T, spin_deg=out['spin_in_flight_deg'])])
    if unit == 'slide':
        trav = _travel(body, st, mv['to'], arena, grip, g)
        if trav is None:
            return None
        u = _unit((mv['to'][0] - st.p[0], mv['to'][1] - st.p[1]))
        for _ in range(3):                # the drop point, where the slide at the speed then ends at the cover
            sp = max(trav[1].v[0] * u[0] + trav[1].v[1] * u[1], 0.1)
            sl = slide(body, sp, g, grip=grip)
            drop = (mv['to'][0] - u[0] * sl['reach_m'], mv['to'][1] - u[1] * sl['reach_m'])
            trav = _travel(body, st, drop, arena, grip, g)
            if trav is None:
                return None
        sp = max(trav[1].v[0] * u[0] + trav[1].v[1] * u[1], 0.1)
        sl = slide(body, sp, g, grip=grip)
        T = sl['sink_s'] + sl['slide_s']
        return (trav[0] + T, State(p=tuple(mv['to']), v=(0.0, 0.0), posture='down'),
                trav[2] + [dict(unit='slide', seconds=T, speed_m_s=sp)])
    if unit == 'pull_down':
        out = pull_down(body, mv.get('depth', 0.5), g, arm_bw)
        return out['seconds'], State(p=st.p, v=(0.0, 0.0), posture='ready'), [dict(unit='pull_down',
                                                                                  seconds=out['seconds'])]
    raise ValueError(f'unknown move {unit!r}')


def route_time(body, start, moves, arena, grip, g, wall_grip=WALL_GRIP, arm_bw=ARM_PULL_BW):
    """Times a sequence of moves from a State and checks each is feasible. Moves are dicts: {'unit': 'floor', 'to':
    (x, y), 'arrive': 'com' | 'hands'}; {'unit': 'vault' | 'kong_vault', 'box': i, 'dir': (1, 0) etc.}; {'unit':
    'tic_tac', 'box': i, 'face': 'x0' | 'x1' | 'y0' | 'y1', 'at': 0-1 along it, 'aim': (x, y)}; {'unit': 'over',
    'box': i, 'face': ...} (wall_run, cat_leap, climb_up and drop_roll across a wall); {'unit': 'lache', 'bar': j,
    'sign': 1 | -1}; {'unit': 'slide', 'to': (x, y)}; {'unit': 'pull_down', 'depth': m}. Each surface move is reached
    by floor travel. Returns ok, the seconds, the records with their units, the end State and its speed."""
    st, total, records = start, 0.0, []
    for mv in moves:
        out = _step(body, st, mv, arena, grip, g, wall_grip, arm_bw)
        if out is None:
            return dict(ok=False, seconds=math.inf, failed=mv, records=records, end=st)
        t, st, rec = out
        total += t
        records += rec
    return dict(ok=True, seconds=total, records=records, end=st, end_speed_m_s=math.hypot(*st.v),
                units=[r['unit'] for r in records if r['unit'] != 'floor'])


def candidates(arena):
    """The arena's surface moves for the search: vaults over every box along both axes each way (the take-off decides
    which are possible: on the Moon a running take-off clears a 3 m wall), a foot's push at three places on every face
    of a box tall enough, every box taller than VAULT_MAX_M climbed over from each face, and every bar both ways."""
    out = []
    for i, b in enumerate(arena.boxes):
        out += [dict(unit=k, box=i, dir=d) for k in ('vault', 'kong_vault')
                for d in ((1.0, 0.0), (-1.0, 0.0), (0.0, 1.0), (0.0, -1.0))]
        if b.height >= TIC_TAC_MIN_M:
            out += [dict(unit='tic_tac', box=i, face=f, at=a, push=k) for f in ('x0', 'x1', 'y0', 'y1')
                    for a in (0.25, 0.5, 0.75) for k in ('full', 'aim')]
        if b.height > VAULT_MAX_M:
            out += [dict(unit='over', box=i, face=f) for f in ('x0', 'x1', 'y0', 'y1')]
    out += [dict(unit='lache', bar=j, sign=s) for j in range(len(arena.bars)) for s in (1, -1)]
    return out


def quickest_route(body, start, goal, arena, grip, g, max_moves=3, beam=8, wall_grip=WALL_GRIP):
    """The quickest route from a State to a goal on the floor plan with up to max_moves surface moves (a beam search:
    at each depth the `beam` partial routes with the least time plus the straight distance left at V_MAX_M_S grow by
    every candidate), beside floor travel around the boxes and the open floor's straight leg."""
    goal = (float(goal[0]), float(goal[1]))
    finish = [dict(unit='floor', to=goal)]
    floor = route_time(body, start, finish, arena, grip, g, wall_grip)
    open_floor = route_time(body, start, finish, Arena(), grip, g, wall_grip)
    best = dict(floor, moves=finish)
    frontier = [(0.0, start, [])]
    cands = candidates(arena)
    for _ in range(max_moves):
        grown = []
        for t0, st, moves in frontier:
            for c in cands:
                mv = dict(c, aim=goal) if c['unit'] == 'tic_tac' else c
                out = _step(body, st, mv, arena, grip, g, wall_grip)
                if out is None:
                    continue
                t1, st1, _ = out
                end = _step(body, st1, finish[0], arena, grip, g, wall_grip)
                if end is not None and t0 + t1 + end[0] < best['seconds'] - 1e-9:
                    best = dict(route_time(body, start, moves + [mv] + finish, arena, grip, g, wall_grip),
                                moves=moves + [mv] + finish)
                grown.append((t0 + t1 + math.dist(st1.p, goal) / V_MAX_M_S, t0 + t1, st1, moves + [mv]))
        grown.sort(key=lambda x: x[0])
        frontier = [(t, st, moves) for _, t, st, moves in grown[:beam]]
    return dict(best=best, floor_only=floor, open_floor=open_floor,
                saves_s=floor['seconds'] - best['seconds'] if floor['ok'] else None)
