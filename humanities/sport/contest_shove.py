"""Shoving for position: two players contesting a spot, as in a box-out, on the Moon beside Earth, and what decides it.

    from humanities.sport import contest_shove
    contest_shove.screen()      # the contest's results; contests.py writes them into results/contests.json

Two players meet chest to chest over a spot, one pushing ('shove') and one holding ('brace'). Each pushes or holds
through its hands, HANDS_ABOVE_COM_M above its centre of mass (vocabulary.HANDS_ABOVE_COM_M), so the push's height
above the floor is the posture's centre of mass plus that.

The measures:
1. The steady push a body makes or holds: its floor force is held by its grip, mu m g, and by its lean, since the push
   at the hands' height turns the body back about its rear support, which its weight holds while F hands <= m g lead
   (lead: the support's farthest extent behind the centre of mass along the push, reach.lead0). By stance: the split
   stance facing the push (lead split_reach), the wide crouch squared up to it with the feet across the push line
   (fore_aft_reach: the vocabulary's 'side-on' figure), the wide crouch turned with the feet along the push line
   (ready_half_stance), and upright (fore_aft_reach at the standing height).
2. Moving an opponent: the pusher at its limit against the holder's hold, the two moving together. The pusher's own
   inertia counts in its lean (the d'Alembert force m a at its centre of mass's height) and its grip carries both
   bodies' acceleration; the holder resists with its steady hold throughout. The time to move the holder DISTANCE_M
   from rest, for equal bodies, the larger against the dance body and back, and a pusher leaning DEEPER_M further
   forward over its split stance.
3. Building a lean: the centre of mass moves forward as a linear inverted pendulum about the rear support,
   x'' = (g/h)(x - x_eq), where m g x_eq / hands is the push still held; the time to deepen the lean by DEEPER_M while
   holding none or half of the push, on gravity's clock sqrt(h/g).
4. The opponent yields: the contact vanishes. The pusher's centre of pressure stays where it was for REACTION_S (the
   legs still push), then moves to the front foot; if the capture point x + v/omega lies beyond it, a step (STEP_S)
   lands at the capture point. Reported: the forward travel until the body is still again, the time until it is
   within 1 cm of rest over its support, and whether a step was needed ('fall_and_rise' is the cost if it is not
   taken). Time-stepped.
5. Shoving a jumper: in its dip's fall (envelope.dip_s's free fall) the floor holds nothing, since the body's weight
   has left the floor; in its push the take-off's mean floor force N holds min(mu, lead/hands) N against a shove
   (squared up, a two-footed take-off); aloft at the top of an envelope.APEX_TAKE_OFF_M_S jump, a push F over
   CONTACT_S changes its travel by F t / m, moves the landing by that times the time left, and turns the body by
   F t lever / I by the landing (envelope.attitude_error_deg's form, lever HANDS_ABOVE_COM_M), laid out and tucked.

Assumptions, beside the shared ones in state.py: the contact force is horizontal; a holder resists with its steady
hold while it is moved (sliding or giving ground step by step, at the same force); the push's height and a body's
centre of mass height stay fixed through a lean; the pusher stands on the floor when it shoves a jumper, at its own
steady push; a recovery step can land at any capture point within 1.5 m.
"""
from __future__ import annotations

import math

from humanities.sport import envelope as env
from humanities.sport import state
from humanities.sport.vocabulary import HANDS_ABOVE_COM_M

STANCES = {                                   # posture, the support's extent behind the centre of mass along the push
    'split, facing the push': ('split', 'split_reach_m'),
    'wide crouch, squared up': ('ready', 'fore_aft_reach_m'),
    'wide crouch, turned': ('ready', 'ready_half_stance_m'),
    'upright': ('upright', 'fore_aft_reach_m'),
}
DISTANCE_M = 0.5
DEEPER_M = (0.1, 0.2)
CONTACT_S = 0.2
SETTLE_M = 0.01
STEP_MAX_M = 1.5
DT = 1e-3
WORKERS = 3


def hands_m(body, posture):
    """The push's height above the floor: the posture's centre of mass plus HANDS_ABOVE_COM_M."""
    return body.com_m(posture) + HANDS_ABOVE_COM_M


def steady(body, stance, g, grip, deeper_m=0.0):
    """The steady push a body makes or holds in a stance: the least of its grip and its lean."""
    posture, attr = STANCES[stance]
    lead = getattr(body, attr) + deeper_m
    weight = body.mass_kg * g
    grip_n, lean_n = grip * weight, weight * lead / hands_m(body, posture)
    return dict(push_N=min(grip_n, lean_n), grip_N=grip_n, lean_N=lean_n, share_of_weight=min(grip_n, lean_n) / weight,
                limited_by='grip' if grip_n < lean_n else 'lean', lead_m=lead, hands_m=hands_m(body, posture))


def move(pusher, holder, g, grip, holder_stance='split, facing the push', deeper_m=0.0, distance=DISTANCE_M):
    """The pusher at its limit in the split stance (leaning deeper_m further) against the holder's steady hold. The pair
    moves at the least of the pusher's lean, (m_A g lead/hands - F_B)/(m_B + m_A h/hands), and its grip,
    (mu m_A g - F_B)/(m_A + m_B); the forces then balance each body's momentum."""
    F_B = steady(holder, holder_stance, g, grip)['push_N']
    A = steady(pusher, 'split, facing the push', g, grip, deeper_m)
    mA, mB, h, hands = pusher.mass_kg, holder.mass_kg, pusher.com_m('split'), A['hands_m']
    a_lean = (A['lean_N'] - F_B) / (mB + mA * h / hands)
    a_grip = (A['grip_N'] - F_B) / (mA + mB)
    a = min(a_lean, a_grip)
    if a <= 0:
        return dict(moves=False, holder_hold_N=F_B, pusher_push_N=A['push_N'], acceleration_m_s2=0.0, seconds=None)
    contact = F_B + mB * a
    return dict(moves=True, holder_hold_N=F_B, pusher_push_N=A['push_N'], acceleration_m_s2=a,
                seconds=math.sqrt(2 * distance / a), contact_N=contact, pusher_floor_N=contact + mA * a,
                limited_by='lean' if a_lean < a_grip else 'grip', speed_after_m_s=math.sqrt(2 * distance * a))


def lean_build_s(body, g, deeper_m, held_share):
    """The time to move the centre of mass deeper_m further forward over the split stance's rear foot while holding
    held_share of the push: x(t) = x_eq + (x0 - x_eq) cosh(t / tau), tau = sqrt(h/g)."""
    tau = math.sqrt(body.com_m('split') / g)
    x0 = body.split_reach_m
    x_eq = held_share * x0
    return tau * math.acosh((x0 + deeper_m - x_eq) / (x0 - x_eq))


def yielded(body, g, grip, deeper_m=0.0, dt=DT, t_max=8.0):
    """The pusher's lurch when the opponent gives way, on the linear inverted pendulum x'' = omega^2 (x - p), x from
    the rear foot, p the centre of pressure: p holds the push's place for REACTION_S, then goes to the front foot, and a
    step lands at the capture point if that lies beyond the front foot."""
    h = body.com_m('split')
    w2 = g / h
    w = math.sqrt(w2)
    front = body.split_reach_m + body.fore_aft_reach_m
    x0 = body.split_reach_m + deeper_m
    F = steady(body, 'split, facing the push', g, grip, deeper_m)['push_N']
    p = x0 - F * hands_m(body, 'split') / (body.mass_kg * g)     # where the floor's push sat while it balanced
    x, v, t, peak, step = x0, 0.0, 0.0, 0.0, False
    n_react, n_step = int(round(state.REACTION_S / dt)), int(round(state.STEP_S / dt))
    i = 0
    while True:
        if i == n_react:
            p = front
            step = x + v / w > front
        if step and i == n_react + n_step:
            break
        if not step and i >= n_react:
            break
        v += w2 * (x - p) * dt
        x += v * dt
        peak = max(peak, v)
        i += 1
        t = i * dt
        if t > t_max:
            break
    xi = x + v / w                                   # the capture point: where the support must be to stop the body
    settle = (state.REACTION_S + (state.STEP_S if step else 0.0)
              + max(0.0, math.log(max(abs(x - xi), 1e-12) / SETTLE_M)) / w)
    return dict(step_needed=step, travel_m=xi - x0, recover_s=settle, peak_speed_m_s=peak,
                step_length_m=(xi - front) if step else 0.0, step_feasible=(xi - front) <= STEP_MAX_M if step else True)


def jumper(pusher, target, g, grip, contact_s=CONTACT_S):
    """A shove at a jumper: in the dip's fall, in the push, and aloft at the top of an APEX_TAKE_OFF_M_S jump."""
    F = steady(pusher, 'split, facing the push', g, grip)['push_N']
    m = target.mass_kg
    depth = target.dip_m
    fall_s = math.sqrt(2 * depth * env.BRAKE_M_S2 / (g + env.BRAKE_M_S2) / g)
    v = target.take_off_m_s(g)
    N = m * (v * v / (2 * target.push_off_m) + g)
    holds = min(grip, target.fore_aft_reach_m / hands_m(target, 'ready')) * N
    descent = env.APEX_TAKE_OFF_M_S / g
    t = min(contact_s, descent)
    left = descent - t / 2
    dv = F * t / m
    dL = F * t * HANDS_ABOVE_COM_M
    return dict(shove_N=F, dip_fall_s=fall_s, push_floor_N=N, push_holds_N=holds, push_absorbs=holds >= F,
                aloft=dict(contact_s=t, descent_s=descent, speed_change_m_s=dv, landing_moved_m=dv * left,
                           angular_momentum_N_m_s=dL,
                           turn_by_landing_deg={k: math.degrees(dL / target.inertia_kg_m2[k] * left)
                                                for k in ('layout', 'tuck')}))


def contest(g, grip):
    """Every measure for one world and grip."""
    B = state.BODIES
    out = dict(steady={k: {s: _r(steady(b, s, g, grip)) for s in STANCES} for k, b in B.items()})
    out['move'] = [dict(pusher=a, holder=b, holder_stance=s, deeper_m=d, **_r(move(B[a], B[b], g, grip, s, d)))
                   for a, b in (('dancer', 'dancer'), ('larger', 'dancer'), ('dancer', 'larger'))
                   for s in ('split, facing the push', 'wide crouch, squared up') for d in (0.0,) + DEEPER_M]
    out['lean_build_s'] = {k: {f'{d:g}_m': {f'holding_{s:g}': round(lean_build_s(b, g, d, s), 3) for s in (0.0, 0.5)}
                               for d in DEEPER_M} for k, b in B.items()}
    out['yield'] = {k: {f'{d:g}_m': _r(yielded(b, g, grip, d)) for d in (0.0,) + DEEPER_M} for k, b in B.items()}
    out['jumper'] = {f'{a}_on_{b}': _r(jumper(B[a], B[b], g, grip)) for a, b in (('dancer', 'dancer'),
                                                                                 ('larger', 'dancer'))}
    return out


def _r(d):
    if isinstance(d, dict):
        return {k: _r(v) for k, v in d.items()}
    return round(d, 3) if isinstance(d, float) else d


def screen(worlds=None, grips=None, workers=WORKERS):
    """Closed forms and short pendulum runs: quick enough to run in one process whatever workers says."""
    worlds = worlds or list(state.WORLDS)
    grips = grips or state.GRIPS
    return dict(hands_above_com_m=HANDS_ABOVE_COM_M, distance_m=DISTANCE_M, contact_s=CONTACT_S,
                reaction_s=state.REACTION_S, step_s=state.STEP_S,
                units=dict(push='shove', hold='brace', turn_to_face='pivot', failed_recovery='fall_and_rise'),
                cases=[dict(world=w, grip=grip, **contest(g, grip)) for w, g in worlds for grip in grips])
