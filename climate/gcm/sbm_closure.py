"""Closure (c) and route L in the simplified Betts-Miller scheme, as the Terluna patch puts them into PlaSim
(SBM_CLOSURE_EDITS here, exoplasim_run.SBM_CLOSURE_PATCHES), one column at a time in Python, for its tests (3 November
2026).

Both were built offline in PlaSim's single column under the CM1 boxes' own forcing (trigger_column.py on the research
drive, session-tools/agents/gap_boxes; its notes, notes_convection.md, from "Closure (c) built" on), and every constant
here is that column's as set there.

Closure (c) (rainmod's nclos_bm 1) takes sbm.py's diluted parcel (2e-4 per metre of rise above the lowest level's
condensation level), its level of free convection, its top and the 200 hPa depth test, and relaxes the whole convecting
column, from the parcel's source to its top, over 18 hours toward the parcel's temperature and each layer toward the
humidity Romps's (2014) relation gives for its own detrainment, delta / (delta + gamma). gamma = L Gamma / (R_v T^2) -
g / (R_d T) comes from the column's own lapse; delta is CM1's cloud-body detrainment, 0.3 per km, plus the organized
detrainment -1/2 d ln K / dz where the updraught's kinetic energy K falls, from the IFS's budget of K (Cy36r1, Part IV,
eqs. 6.9-6.11) with the diluted parcel's virtual buoyancy, the IFS's constants and the parcel's dilution as the mixing.
The source keeps the cloud body's humidity, delta_body / (delta_body + gamma). As in the scheme, the rain is the vapour
removed and the heating is shifted by one amount through the relaxed layers so that it releases the rain's latent heat.
A column whose parcel is buoyant but whose cloud is too thin, or whose vapour lies below the reference, goes to
PlaSim's shallow scheme with the shallow cap's top.

Route L (nroutel_bm 1, with closure (c)) starts deep convection aloft where the lowest level's parcel does not convect
deep. Each layer between sigma 0.60 and 0.94 above the boundary layer's top, from which closure (c) launched with the
layer brought to saturation rains, is a candidate; the route starts when any candidate has been moistened over the
previous 24 hours by the processes other than convection, runs on for 12 hours after its last start while candidates
remain, and launches closure (c) from the most humid candidate. The boundary layer's top is the gust patch's depth
(surface_gust.gust_column, gustcol in Fortran) on the state entering the convection with the step's surface fluxes,
the surface pressure in stable or neutral air. The moistening is the column's measure: each level's change of relative
humidity a step, (dq/dt - rh q_s L dT/dt / (R_v T^2)) / q_s, by the dynamics (PlaSim's adiabatic tendencies of
temperature and vapour) and the radiation on the state entering the physics, summed exactly over the last 48 steps of
30 minutes. The order a column takes: closure (c) from the lowest level where its parcel convects deep; otherwise
route L where it is on; otherwise the shallow scheme where the lowest level's parcel is buoyant.

Where PlaSim's arrays depart from the column's (each stated in the patch's record, exoplasim_run.convection_setting):
the moistening counts the GCM's own dynamics and radiation where the column counted its imposed large-scale terms and
radiation; the dynamics' vapour tendency (dqt) is that of q p_s, so it carries the surface pressure's own effect on the
relative humidity at a sigma level; the boundary layer's top takes the winds after the vertical diffusion where the
column held CM1's, and PlaSim's own constants throughout. The change between the state leaving the convection at one
step and the state entering it at the next is not used: in PlaSim's leapfrog the physics advances its provisional
state by two steps while the model's state advances by one, so that difference counts twice the non-convective change
less the step's whole change, convection included (fortran_cl/leapfrog_toy.py on the research drive).

The Fortran leaves the parcel unintegrated above the first level where it is not buoyant over its level of free
convection, and forms route L's candidates only where they can change a decision: only the moistened layers while the
lowest level convects deep, and none while no layer in the range is moistened and the route's run after its last start
is over. Neither changes an output; the tests compare every output with this copy, which forms everything.

Arrays run top first, as PlaSim's levels do.
"""
from __future__ import annotations

import math

import numpy as np

from climate.gcm import massflux as mf
from climate.gcm import sbm
from climate.gcm import surface_gust as sg

# The column's constants as trigger_column.py sets them (LUNAR, COMMON, the closure's and route L's).
TAU_S = 18.0 * 3600.0              # the relaxation time (s)
DILUTION = 2.0e-4                  # the parcel's dilution per metre of rise above its condensation level
DEEP_PA = 20000.0                  # the depth test's threshold (Pa); the shallow cap holds a shallow top within it
MAX_MIX = sbm.MAX_MIX              # the diluted parcel's largest mixing exponent in one step
DELTA_BODY = 0.3e-3                # per metre: CM1's cloud-body detrainment
IFS_KE = dict(virtual_mass=0.5, turbulence=2.0, drag=0.506 * 1.875, w_base=1.0)
KE_SUBSTEP_M = 50.0                # the most metres in one step of the kinetic energy's budget
RV = 461.5                         # J/kg/K: the vapour's gas constant in gamma and in the history
XLV = 2.501e6                      # J/kg: the latent heat in the history's change of relative humidity
LIFT_SIGMA = (0.60, 0.94)          # route L's layers, as fractions of the surface pressure (trigger_design2)
MEMORY_S = 12.0 * 3600.0           # route L runs on this long after its last start
WINDOW_S = 86400.0                 # the span of its history
WINDOW = 48                        # the history's steps, PlaSim's 30 minutes each
STEP_S = 1800.0                    # PlaSim's deltsec
STEP2_S = 3600.0                   # PlaSim's deltsec2, the scheme's shortest relaxation time

CLOSURE_MARKER = 'Terluna: closure (c) and route L'


# ------------------------------------------------------------------------------------------------------- closure (c)

def gradient(f, x) -> np.ndarray:
    """numpy.gradient(f, x) for unevenly spaced x, written out as numpy forms it (second-order differences inside,
    first-order ones at the ends), so that the Fortran can follow it operation by operation."""
    f, x = np.asarray(f, float), np.asarray(x, float)
    n = f.size
    out = np.zeros(n)
    dx = np.diff(x)
    if np.all(dx == dx[0]):                                 # numpy's evenly spaced form, which heights never take
        out[1:-1] = (f[2:] - f[:-2]) / (2.0 * dx[0])
    else:
        for i in range(1, n - 1):
            dx1, dx2 = dx[i - 1], dx[i]
            a = -dx2 / (dx1 * (dx1 + dx2))
            b = (dx2 - dx1) / (dx1 * dx2)
            c = dx1 / (dx2 * (dx1 + dx2))
            out[i] = a * f[i - 1] + b * f[i] + c * f[i + 1]
    out[0] = (f[1] - f[0]) / dx[0]
    out[n - 1] = (f[n - 1] - f[n - 2]) / dx[n - 2]
    return out


def romps_reference(t, q, sigma, ps, tp, top, c, dilution=DILUTION, delta_body=DELTA_BODY) -> dict:
    """Closure (c)'s reference relative humidity on a column whose lowest level is the parcel's source, from the level
    above the source to the parcel's top (trigger_column.romps_reference): Romps's delta / (delta + gamma), 1 where
    gamma is not positive and at most 1. gamma = L Gamma / (R_v T^2) - g / (R_d T) from the column's own lapse Gamma;
    delta is delta_body plus the organized detrainment -1/2 d ln K / dz where the updraught's kinetic energy K falls,
    K from the source's 1/2 w_base^2 by the IFS's budget dK/dz = -2 eps (1 + drag) K + B / (turbulence (1 +
    virtual_mass)), B the diluted parcel's virtual buoyancy, integrated in steps of at most KE_SUBSTEP_M metres with K
    held at least 1e-6. Returns the reference (nan outside the cloud), delta and gamma (per metre), K and the heights
    above the source."""
    t, q, sigma, tp = (np.asarray(a, float) for a in (t, q, sigma, tp))
    n = t.size
    g, rd = c['ga'], c['gascon']
    p = sigma * ps
    eps_v = 1.0 / c['rdbrv'] - 1.0
    tv = t * (1.0 + eps_v * q)
    z = np.zeros(n)
    for k in range(n - 2, -1, -1):
        z[k] = z[k + 1] + rd * 0.5 * (tv[k] + tv[k + 1]) / g * math.log(p[k + 1] / p[k])
    lapse = -gradient(t, z)
    gamma = c['ALV'] * lapse / (RV * t ** 2) - g / (rd * t)
    qs_p = np.array([sbm.saturation(a, b, c) for a, b in zip(tp, p)])
    buoy = g * (tp * (1.0 + eps_v * qs_p) - tv) / tv
    a = 2.0 * dilution * (1.0 + IFS_KE['drag'])
    b = 1.0 / (IFS_KE['turbulence'] * (1.0 + IFS_KE['virtual_mass']))
    kin = np.full(n, np.nan)
    kin[n - 1] = 0.5 * IFS_KE['w_base'] ** 2
    for k in range(n - 2, max(top, 0) - 1, -1):
        steps = max(1, int(math.ceil((z[k] - z[k + 1]) / KE_SUBSTEP_M)))
        h = (z[k] - z[k + 1]) / steps
        kk = kin[k + 1]
        for s in range(steps):
            f = (s + 0.5) / steps
            kk = max(kk + h * (-a * kk + b * ((1.0 - f) * buoy[k + 1] + f * buoy[k])), 1.0e-6)
        kin[k] = kk
    delta = np.full(n, np.nan)
    rh = np.full(n, np.nan)
    for k in range(n - 2, max(top, 0) - 1, -1):
        org = max(0.0, -0.5 * math.log(kin[k] / kin[k + 1]) / (z[k] - z[k + 1]))
        delta[k] = delta_body + org
        rh[k] = min(1.0, delta[k] / (delta[k] + gamma[k])) if gamma[k] > 0.0 else 1.0
    return dict(rh=rh, delta=delta, gamma=gamma, kinetic=kin, z=z)


def closure(t, q, sigma, dsigma, ps, c, akap, step2=STEP2_S, tau=TAU_S) -> dict:
    """Closure (c) on a column whose lowest level is the parcel's source (trigger_column.closure_run with CLOSURE
    'romps'): the diluted parcel (sbm.diluted_parcel at DILUTION), its level of free convection and its top as
    sbm.column finds them, the depth test from the lowest level's condensation level, and where the cloud is deeper
    than DEEP_PA every level from the top to the source relaxed over tau (at least step2) toward the parcel's
    temperature and romps_reference's humidity, the source's the cloud body's. Returns whether it rains (deep), the
    parcel's top (None where it is buoyant nowhere), the top it hands the shallow scheme where it does not rain
    (shallow_top: the parcel's top for a cloud the depth test calls thin, else the shallow cap's), the condensation
    level (base, Pa), the tendencies (K/s, kg/kg/s), the rain formed at each level (kg/m2/s), the reference and the
    parcel's temperature."""
    t, q = np.asarray(t, float), np.asarray(q, float)
    sigma, dsigma = np.asarray(sigma, float), np.asarray(dsigma, float)
    n = t.size
    p, m = ps * sigma, ps * dsigma / c['ga']
    out = dict(deep=False, top=None, shallow_top=None, base=None, dtdt=np.zeros(n), dqdt=np.zeros(n),
               rain=np.zeros(n), rh_ref=np.full(n, np.nan))
    tp = sbm.diluted_parcel(t, q, sigma, ps, DILUTION, akap, c, max_mix=MAX_MIX)
    out['reference_t'] = tp
    buoyant = [k for k in range(n - 2, -1, -1) if tp[k] > t[k]]
    if not buoyant:
        return out
    lfc = buoyant[0]
    top = lfc
    for k in range(lfc - 1, -1, -1):
        if tp[k] <= t[k]:
            break
        top = k
    base = sbm.condensation_level(t[-1], q[-1], p[-1], p[top], akap, c)
    out.update(top=top, shallow_top=top, base=base)
    if base - p[top] <= DEEP_PA:                             # the depth test: a thin cloud goes to the shallow scheme
        return out
    ref = romps_reference(t, q, sigma, ps, tp, top, c)
    rh = ref['rh'].copy()
    g0 = ref['gamma'][n - 1]
    rh[n - 1] = min(1.0, DELTA_BODY / (DELTA_BODY + g0)) if g0 > 0.0 else 1.0
    out['rh_ref'] = rh
    cpk = c['acpd'] * (1.0 + c['ADV'] * q)
    tau = max(tau, step2)
    layers = range(top, n)
    dtr = {k: (tp[k] - t[k]) / tau for k in layers}
    dqr = {k: (rh[k] * sbm.saturation(tp[k], p[k], c) - q[k]) / tau for k in layers}
    pq = -sum(dqr[k] * m[k] for k in layers)
    pt = sum(cpk[k] * dtr[k] * m[k] for k in layers)
    cpm = sum(cpk[k] * m[k] for k in layers)
    if pq <= 0.0:                                            # too dry to rain: the shallow cap's top
        out['shallow_top'] = next((k for k in range(top, n) if p[k] >= base - DEEP_PA), n - 1)
        return out
    corr = (c['ALV'] * pq - pt) / cpm
    cond = sum(max(0.0, -dqr[k]) * m[k] for k in layers)
    for k in layers:
        out['dtdt'][k] = dtr[k] + corr
        out['dqdt'][k] = dqr[k]
        if dqr[k] < 0.0:
            out['rain'][k] = -dqr[k] * pq / cond * m[k]
    out.update(deep=True, shallow_top=None)
    return out


def elevated(t, q, sigma, dsigma, ps, k, c, akap, step2=STEP2_S) -> dict:
    """Closure (c) launched from level k brought to saturation (trigger_column.elevated): the column cut at k, the
    level's vapour set to its saturation, its tendencies and rain padded with zeros below k."""
    t, q = np.asarray(t, float), np.asarray(q, float)
    n = t.size
    qk = q[:k + 1].copy()
    qk[k] = sbm.saturation(t[k], sigma[k] * ps, c)
    o = closure(t[:k + 1], qk, sigma[:k + 1], dsigma[:k + 1], ps, c, akap, step2)
    for key in ('dtdt', 'dqdt', 'rain'):
        o[key] = np.concatenate([o[key], np.zeros(n - k - 1)])
    return o


# ------------------------------------------------------------------------------------------------------- route L

def interp(x: float, xp, fp) -> float:
    """numpy.interp(x, xp, fp) for one x and increasing xp, written out as numpy forms it."""
    if x < xp[0]:
        return float(fp[0])
    if x >= xp[-1]:
        return float(fp[-1])
    j = max(i for i in range(len(xp) - 1) if xp[i] <= x)
    if xp[j] == x:
        return float(fp[j])
    slope = (fp[j + 1] - fp[j]) / (xp[j + 1] - xp[j])
    return float(slope * (x - xp[j]) + fp[j])


def pressure_at(h: float, pf, zf, ps: float) -> float:
    """The pressure (Pa) at height h above the ground, linear in log pressure through the ground and the full levels'
    heights zf (trigger_column.pressure_at)."""
    xs = [0.0] + [float(z) for z in zf[::-1]]
    ys = [math.log(ps)] + [math.log(p) for p in pf[::-1]]
    return math.exp(interp(h, xs, ys))


def surface_fluxes(dshfl: float, devap: float, t_low: float, q_low: float, p_low: float, c: dict) -> tuple:
    """The kinematic surface fluxes of heat (K m/s) and vapour (kg/kg m/s), upward, from PlaSim's dshfl (W/m2) and
    devap (m/s of water), downward, over the lowest level's air density and the dry heat capacity, as the column forms
    them from its sensible and latent heat."""
    rho = p_low / (c['gascon'] * t_low * (1.0 + (1.0 / c['rdbrv'] - 1.0) * q_low))
    return -dshfl / (rho * c['acpd']), -1000.0 * devap / rho


def boundary_layer_pressure(t, q, u, v, sigma, sigmah, ps, wth0, wq0, ust, c) -> float:
    """Route L's boundary-layer top (Pa): the gust's depth (surface_gust.gust_column, the bulk Richardson number with
    the thermals' excess) as a pressure through massflux.heights' full levels, the surface pressure where the surface
    buoyancy flux is not upward."""
    gust = sg.gust_column(t, q, u, v, sigma, sigmah, wth0, wq0, ust, c)
    h = gust['h'] if gust['unstable'] else 0.0
    if not h > 0.0:
        return float(ps)
    pf = ps * np.asarray(sigma, float)
    zf, _ = mf.heights(t, q, pf, np.r_[0.0, ps * np.asarray(sigmah, float)], c['ga'], c)
    return pressure_at(h, pf, zf, ps)


def history_entry(t, q, p, rad, dyn_t, dyn_q, c, step=STEP_S) -> np.ndarray:
    """Each level's change of relative humidity over a step by the processes other than convection (trigger_column's
    history): (dq/dt - rh q_s XLV dT/dt / (RV T^2)) / q_s times the step, on the state t, q at pressures p, with dT/dt
    the radiation's heating rad and the dynamics' dyn_t, dq/dt the dynamics' dyn_q."""
    t, q = np.asarray(t, float), np.asarray(q, float)
    qs = np.array([sbm.saturation(a, b, c) for a, b in zip(t, p)])
    rh = q / qs
    return (np.asarray(dyn_q, float) - rh * qs * XLV * (np.asarray(rad, float) + np.asarray(dyn_t, float))
            / (RV * t ** 2)) / qs * step


class RouteL:
    """Route L's state at one grid point: its history, the change of each level's relative humidity in each of the
    last WINDOW steps (a ring indexed by the step's number, empty where a run starts), and the step of its last start
    (None before the first)."""

    def __init__(self, n: int):
        self.hist = np.zeros((WINDOW, n))
        self.last = None

    def record(self, nstep: int, entry) -> np.ndarray:
        """Put this step's change into its slot and return the history's sum, the change over the previous 24
        hours, added slot by slot as the Fortran adds them."""
        self.hist[nstep % WINDOW] = entry
        total = np.zeros(self.hist.shape[1])
        for j in range(WINDOW):
            total = total + self.hist[j]
        return total


def convect(te, qe, sigma, dsigma, ps, c, akap, step2=STEP2_S, route=None) -> dict:
    """One column's convection by closure (c) on its state entering the convection (te, qe), with route L when route
    is given: a dict of the step's number (nstep), its length (step), the 24-hour change of each level's relative
    humidity (drh24), the boundary layer's top (p_bl, Pa) and the RouteL state (state, updated here). The order of
    trigger_column's scheme 'C+L': closure (c) from the lowest level where it rains; otherwise route L where it is on,
    from its most humid candidate; otherwise the shallow scheme where the lowest level's parcel is buoyant. Returns the
    route ('low', 'L', 'shallow' or None), the tendencies and rain, the shallow scheme's top, the source level and
    route L's candidates (level, relative humidity, 24-hour change), start and run."""
    te, qe = np.asarray(te, float), np.asarray(qe, float)
    n = te.size
    sigma = np.asarray(sigma, float)
    low = closure(te, qe, sigma, dsigma, ps, c, akap, step2)
    zero = np.zeros(n)
    out = dict(route=None, dtdt=zero, dqdt=zero, rain=zero, shallow_top=None, source=None, low=low, candidates=[],
               start=False, on=False)
    best = None
    if route is not None:
        p = ps * sigma
        for k in range(n - 1):
            if not (LIFT_SIGMA[0] <= sigma[k] <= LIFT_SIGMA[1] and p[k] < route['p_bl']):
                continue
            o = elevated(te, qe, sigma, dsigma, ps, k, c, akap, step2)
            if not o['deep']:
                continue
            rh = qe[k] / sbm.saturation(te[k], p[k], c)
            out['candidates'].append((k, rh, float(route['drh24'][k])))
            if best is None or rh > best[1]:
                best = (k, rh, o)
        state = route['state']
        out['start'] = any(d > 0.0 for _, _, d in out['candidates'])
        if out['start']:
            state.last = route['nstep']
        out['on'] = out['start'] or (bool(out['candidates']) and state.last is not None
                                     and (route['nstep'] - state.last) * route['step'] <= MEMORY_S)
    if low['deep']:
        out.update(route='low', dtdt=low['dtdt'], dqdt=low['dqdt'], rain=low['rain'], source=n - 1)
    elif out['on'] and best is not None:
        o = best[2]
        out.update(route='L', dtdt=o['dtdt'], dqdt=o['dqdt'], rain=o['rain'], source=best[0])
    elif low['shallow_top'] is not None:
        out.update(route='shallow', shallow_top=low['shallow_top'])
    return out


# ------------------------------------------------------------------------------------------------------- Fortran

# The namelist's values: (Fortran name, value, what it is), closure (c)'s, then route L's.
CLOSURE_VALUES = (
    ('dil_bm', DILUTION, "closure (c)'s parcel: its dilution per metre of rise above its condensation level"),
    ('dmix_bm', MAX_MIX, "the parcel's largest mixing exponent in one step"),
    ('delb_bm', DELTA_BODY, "the cloud body's detrainment per metre (CM1's)"),
    ('kevm_bm', IFS_KE['virtual_mass'], "the updraught's kinetic energy: virtual mass coefficient (IFS)"),
    ('keturb_bm', IFS_KE['turbulence'], 'its turbulence factor'),
    ('kedrag_bm', IFS_KE['drag'], 'its drag'),
    ('kew_bm', IFS_KE['w_base'], 'its speed at the source (m/s)'),
    ('kesub_bm', KE_SUBSTEP_M, 'the most metres in one step of its budget'),
    ('rv_bm', RV, "the vapour's gas constant in gamma and in route L's history (J/kg/K)"),
)
ROUTE_L_VALUES = (
    ('rlsig1_bm', LIFT_SIGMA[0], "route L's layers: from this share of the surface pressure"),
    ('rlsig2_bm', LIFT_SIGMA[1], 'to this one'),
    ('rlmem_bm', MEMORY_S, 'it runs on this long after its last start (s)'),
    ('rlwin_bm', WINDOW_S, 'the span of its history (s)'),
    ('xlv_bm', XLV, "the latent heat in its history's change of relative humidity (J/kg)"),
)


SOURCE = ('trigger_column.py (CLOSURE romps, L_START closure) and notes_convection.md, '
          'session-tools/agents/gap_boxes on the research drive, 3 November 2026')


def closure_record() -> dict:
    """Closure (c) as a run's configuration records it (exoplasim_run.convection_setting)."""
    return dict(form='c', constants={name: value for name, value, _ in CLOSURE_VALUES},
                parcel='the lowest level\'s air lifted undiluted to its condensation level, then diluted per metre of '
                       'rise (sbm.diluted_parcel), its level of free convection and top as the scheme finds them',
                depth='deep only where the cloud from the condensation level to the top is thicker than the depth '
                      'test\'s threshold; a thinner one, or one too dry to rain, takes the shallow scheme with the '
                      'shallow cap\'s top',
                reference='every level from the source to the top toward the parcel\'s temperature and Romps\'s (2014) '
                          'delta / (delta + gamma), delta the cloud body\'s detrainment plus -1/2 d ln K / dz from the '
                          'IFS\'s kinetic-energy budget (Cy36r1 Part IV, eqs. 6.9-6.11) with the parcel\'s virtual '
                          'buoyancy; the source the cloud body\'s delta_body / (delta_body + gamma)',
                relaxation='over tau_s; the rain the vapour removed, the heating shifted by one amount to release its '
                           'latent heat; the rain through PlaSim\'s handling of the scheme\'s convective rain',
                replaces='the reference humidity rh (rh_bm is not read)', source=SOURCE,
                code='climate/gcm/sbm_closure.py (closure; sbmclos, sbmdilp and sbmromps in Fortran)')


def route_l_record() -> dict:
    """Route L as a run's configuration records it, with the ways PlaSim's arrays depart from the single column."""
    return dict(constants={name: value for name, value, _ in ROUTE_L_VALUES}, history_steps=WINDOW,
                order='closure (c) from the lowest level where it rains; else route L where it is on; else the shallow '
                      'scheme where the lowest level\'s parcel is buoyant',
                candidates='layers from rlsig1_bm to rlsig2_bm of the surface pressure above the boundary layer\'s top '
                           'from which closure (c), launched with the layer brought to saturation, rains',
                start='a candidate whose relative humidity the processes other than convection raised over the last '
                      'history_steps steps (their sum above 0)',
                memory='on for rlmem_bm seconds after its last start while candidates remain',
                launch='closure (c) from the most humid candidate',
                boundary_layer='the gust patch\'s depth (gustcol) on the state entering the convection, its winds '
                               'after the vertical diffusion, with the step\'s surface fluxes over the lowest '
                               'level\'s air density; the surface pressure in stable or neutral air',
                moistening='(dq/dt - rh qs xlv_bm dT/dt / (rv_bm T^2)) / qs over each step, dT/dt the dynamics\' '
                           '(stt) and the radiation\'s (dtdtswr + dtdtlwr) heating, dq/dt the dynamics\' vapour (dqt, '
                           'of q p_s), on the state entering the physics',
                restart='its history and last starts in the restart (dsbmrlhist, dsbmrllast); empty from a restart '
                        'without them',
                departures=['the moistening counts the GCM\'s own dynamics and radiation in place of the column\'s '
                            'imposed large-scale terms and radiation; dqt, the tendency of q p_s, adds the surface '
                            'pressure\'s own effect on the relative humidity at a sigma level',
                            'the change between the state leaving the convection and the state entering it at the '
                            'next step is not used: in PlaSim\'s leapfrog it counts twice the non-convective change '
                            'less the step\'s whole change (fortran_cl/leapfrog_toy.py on the research drive)',
                            'the boundary layer\'s top takes the winds after the vertical diffusion where the column '
                            'held CM1\'s, and PlaSim\'s constants throughout'],
                source=SOURCE, code='climate/gcm/sbm_closure.py (convect, RouteL; sbmclosc, sbmpbl and sbmrlhist in '
                                    'Fortran)')


def fortran_literal(value: float) -> str:
    """A Fortran real literal of the value, exact in double precision (the runs build with -fdefault-real-8)."""
    text = repr(float(value))
    return text if ('e' in text or '.' in text) else text + '.'


def namelist_values() -> dict:
    """The rainmod namelist entries (name, text) for closure (c) and route L's constants, as the runner writes them."""
    return {name.upper(): repr(float(value)) for name, value, _ in CLOSURE_VALUES + ROUTE_L_VALUES}


def _declarations() -> str:
    lines = [f'      integer :: nclos_bm = 0 ! {CLOSURE_MARKER}: 1 takes closure (c) in the scheme\n']
    lines += [f'      real :: {name} = {fortran_literal(value)} ! {what}\n' for name, value, what in CLOSURE_VALUES]
    lines += ['      integer :: nroutel_bm = 0 ! 1 adds route L to closure (c)\n']
    lines += [f'      real :: {name} = {fortran_literal(value)} ! {what}\n' for name, value, what in ROUTE_L_VALUES]
    lines += [f'      integer, parameter :: NRLWIN = {WINDOW} ! route L\'s history: its steps\n',
              '      real, allocatable :: drlhist(:,:,:) ! its change of relative humidity by level and step\n',
              '      real, allocatable :: drllast(:) ! the step of its last start at each grid point\n',
              '      integer :: nsbmroute(NHOR) = 0 ! the convection at each grid point: 1 from the lowest level, '
              '2 by route L, 3 shallow\n',
              '      integer :: nsbmsrc(NHOR) = 0 ! the level its deep convection rises from\n']
    return ''.join(lines)


NAMELIST_NAMES = (['nclos_bm'] + [name for name, _, _ in CLOSURE_VALUES] + ['nroutel_bm']
                  + [name for name, _, _ in ROUTE_L_VALUES])
INTEGER_NAMES = ('nclos_bm', 'nroutel_bm')

_CAP_LINE = ("      integer :: nshcap_bm = 0 ! Terluna: the shallow cap: 1 holds a shallow top within dpdeep_bm of its "
             "base\n")
_CAP_NAMELIST = "      namelist/rainmod_nl/nshcap_bm   ! Terluna: the shallow cap\n"
_CAP_BROADCAST = "      call mpbci(nshcap_bm)\n"
_DEPTH_DECLARATION = ("      real :: zsub, zstep, zlift, zsig0, zsig1, zw, zpbase   ! Terluna: the cloud-depth test "
                      "and the mixing within layers\n")
_BEFORE_LOOP = ("      zqe(:,:)=dq(:,1:NLEV)+dqdt(:,1:NLEV)*deltsec2\n!\n      do jhor=1,NHOR\n"
                "       zp(:)=dp(jhor)*sigma(:)\n")
_PARCEL = "!      the parcel, from the lowest level up\n"
_LCL_END = "      end function sbmlcl\n"
_MKSBM_END = "      end subroutine mksbm\n"
_RAININI_END = "      return\n      end subroutine rainini\n"
_RAINSTOP = "      subroutine rainstop\n      use rainmod\n!\n"

# mksbm's own routines (after sbmlcl, so they see its variables): one column by closure (c) with route L, the diluted
# parcel, closure (c) on the column's top kn levels, its reference, the boundary layer's top and route L's history.
ROUTINES = """!
      subroutine sbmclosc(kh)
!     Terluna: closure (c) and route L at grid point kh (nclos_bm 1; climate/gcm/sbm_closure.py, convect): closure (c)
!     from the lowest level where it rains; otherwise, with route L (nroutel_bm 1), closure (c) from the most humid of
!     the layers between rlsig1_bm and rlsig2_bm of the surface pressure above the boundary layer's top from which,
!     brought to saturation, it rains, where the route is on: started this step by such a layer moistened over the
!     history's span (zdrh24 > 0), or running on within rlmem_bm of its last start while such layers remain;
!     otherwise PlaSim's shallow scheme where the lowest level's parcel is buoyant, with the top closure (c) gives it
      integer, intent(in) :: kh
      real :: zt(NLEV),zq(NLEV),zdtl(NLEV),zdql(NLEV),zrnl(NLEV),zrhl(NLEV)
      real :: zt2(NLEV),zq2(NLEV),zdte(NLEV),zdqe(NLEV),zrne(NLEV),zrhe(NLEV)
      real :: zdtb(NLEV),zdqb(NLEV),zrnb(NLEV)
      real :: zpbl,zrhc,zrhbest
      logical :: ldeep,ldeepe,lstart,lon,lmoist,lrun
      integer :: jl,jtopl,jshl,jtope,jshe,jtopb,jbest,icand
!
      zt(:)=zte(kh,:)
      zq(:)=zqe(kh,:)
      call sbmclos(NLEV,zt,zq,dp(kh),ldeep,jtopl,jshl,zdtl,zdql,zrnl,zrhl)
      nsbmroute(kh)=0
      nsbmsrc(kh)=0
      jbest=0
      jtopb=0
      lon=.false.
      if(nroutel_bm == 1) then
       zpbl=sbmpbl(kh,zt,zq)
       icand=0
       lstart=.false.
       zrhbest=0.
!      the candidates where they can count: only the moistened layers, for the start, while the lowest level
!      convects deep; none while no layer in the range is moistened and the route's run after its last start is over
       lmoist=.false.
       do jl=1,NLEM
        if(sigma(jl) >= rlsig1_bm .and. sigma(jl) <= rlsig2_bm .and. zp(jl) < zpbl .and. zdrh24(kh,jl) > 0.)    &
     &   lmoist=.true.
       enddo
       lrun=(real(nstep)-drllast(kh))*deltsec <= rlmem_bm
       do jl=1,NLEM
        if(sigma(jl) < rlsig1_bm .or. sigma(jl) > rlsig2_bm .or. zp(jl) >= zpbl) cycle
        if(ldeep .and. (lstart .or. .not. zdrh24(kh,jl) > 0.)) cycle
        if(.not. (ldeep .or. lmoist .or. lrun)) cycle
        zt2(1:jl)=zt(1:jl)
        zq2(1:jl)=zq(1:jl)
        zq2(jl)=sbmqs(zt(jl),zp(jl))
        call sbmclos(jl,zt2(1:jl),zq2(1:jl),dp(kh),ldeepe,jtope,jshe,zdte(1:jl),zdqe(1:jl),zrne(1:jl),zrhe(1:jl))
        if(.not. ldeepe) cycle
        icand=icand+1
        if(zdrh24(kh,jl) > 0.) lstart=.true.
        zrhc=zq(jl)/sbmqs(zt(jl),zp(jl))
        if(icand == 1 .or. zrhc > zrhbest) then
         zrhbest=zrhc
         jbest=jl
         jtopb=jtope
         zdtb(:)=0.
         zdqb(:)=0.
         zrnb(:)=0.
         zdtb(1:jl)=zdte(1:jl)
         zdqb(1:jl)=zdqe(1:jl)
         zrnb(1:jl)=zrne(1:jl)
        endif
       enddo
       if(lstart) drllast(kh)=real(nstep)
       lon=lstart .or. (icand > 0 .and. (real(nstep)-drllast(kh))*deltsec <= rlmem_bm)
      endif
      if(ldeep) then
       call sbmapply(kh,jtopl,NLEV,zdtl,zdql,zrnl)
       nsbmroute(kh)=1
       nsbmsrc(kh)=NLEV
      elseif(lon .and. jbest > 0) then
       call sbmapply(kh,jtopb,jbest,zdtb,zdqb,zrnb)
       nsbmroute(kh)=2
       nsbmsrc(kh)=jbest
      elseif(jshl > 0) then
       if(nshallow == 1) then
        ishallow(kh)=1
        ilift(kh)=NLEP
        itop(kh)=jshl
       endif
       nsbmroute(kh)=3
      endif
      end subroutine sbmclosc
!
      subroutine sbmapply(kh,ktop,kn,pdtr,pdqr,prain)
!     Terluna: closure (c)'s tendencies and rain at grid point kh, from the top ktop to the source kn, as the scheme's
!     deep convection adds them
      integer, intent(in) :: kh,ktop,kn
      real, intent(in) :: pdtr(NLEV),pdqr(NLEV),prain(NLEV)
      integer :: jl
      do jl=ktop,kn
       dtdt(kh,jl)=dtdt(kh,jl)+pdtr(jl)
       dqdt(kh,jl)=dqdt(kh,jl)+pdqr(jl)
       zdtc(kh,jl)=pdtr(jl)
       if(pdqr(jl) < 0.) then
        dprcl(kh,jl)=prain(jl)
        icclev(kh,jl)=1
        icctot(kh)=icctot(kh)+1
        if(sigma(jl) < 0.4 .and. jl == ktop) icclev(kh,jl)=2
       endif
      enddo
      end subroutine sbmapply
!
      subroutine sbmclos(kn,pt,pq,pps,ldeep,ktop,kshal,pdtr,pdqr,prain,prh)
!     Terluna: closure (c) on the levels 1 to kn, level kn the parcel's source (sbm_closure.closure): the diluted
!     parcel (sbmdilp), its level of free convection and top, the depth test from its condensation level, and where the
!     cloud is deeper than dpdeep_bm every level from the top to the source relaxed over tau_bm toward the parcel's
!     temperature and the reference humidity of sbmromps, the heating shifted by one amount so that it releases the
!     rain's latent heat. ldeep where it rains; ktop the parcel's top and kshal the top for the shallow scheme where it
!     does not (the parcel's top for a thin cloud, the shallow cap's for one too dry to rain; 0 where the parcel is
!     buoyant nowhere); the tendencies (K/s, kg/kg/s), the rain formed in dprcl's units and the reference by level
      integer, intent(in) :: kn
      real, intent(in) :: pt(kn),pq(kn),pps
      logical, intent(out) :: ldeep
      integer, intent(out) :: ktop,kshal
      real, intent(out) :: pdtr(kn),pdqr(kn),prain(kn),prh(kn)
      real :: ztp(NLEV),zcpk(NLEV),zbase,ztauc,zsum,zpq,zpt,zcpm,zcorr,zcond
      integer :: jl,jlfc
!
      ldeep=.false.
      ktop=0
      kshal=0
      pdtr(:)=0.
      pdqr(:)=0.
      prain(:)=0.
      prh(:)=0.
      call sbmdilp(kn,pt,pq,pps,ztp)
      jlfc=0
      do jl=kn-1,1,-1
       if(ztp(jl) > pt(jl)) then
        jlfc=jl
        exit
       endif
      enddo
      if(jlfc == 0) return
      ktop=jlfc
      do jl=jlfc-1,1,-1
       if(ztp(jl) <= pt(jl)) exit
       ktop=jl
      enddo
      kshal=ktop
      zbase=sbmlcl(pt(kn),pq(kn),zp(kn),zp(ktop))
      if(zbase-zp(ktop) <= dpdeep_bm) return
      call sbmromps(kn,pt,pq,pps,ztp,ktop,prh)
      ztauc=max(tau_bm,deltsec2)
      zcpk(1:kn)=acpd*(1.+ADV*pq(1:kn))
      zsum=0.
      zpt=0.
      zcpm=0.
      do jl=ktop,kn
       pdtr(jl)=(ztp(jl)-pt(jl))/ztauc
       pdqr(jl)=(prh(jl)*sbmqs(ztp(jl),zp(jl))-pq(jl))/ztauc
       zsum=zsum+pdqr(jl)*zm(jl)
       zpt=zpt+zcpk(jl)*pdtr(jl)*zm(jl)
       zcpm=zcpm+zcpk(jl)*zm(jl)
      enddo
      zpq=-zsum
      if(zpq <= 0.) then
       do jl=ktop,kn
        if(zp(jl) >= zbase-dpdeep_bm) exit
       enddo
       kshal=min(jl,kn)
       pdtr(:)=0.
       pdqr(:)=0.
       return
      endif
      zcorr=(ALV*zpq-zpt)/zcpm
      zcond=0.
      do jl=ktop,kn
       zcond=zcond+max(0.,-pdqr(jl))*zm(jl)
      enddo
      do jl=ktop,kn
       pdtr(jl)=pdtr(jl)+zcorr
       if(pdqr(jl) < 0.) prain(jl)=-pdqr(jl)*zpq/zcond*dsigma(jl)
      enddo
      ldeep=.true.
      kshal=0
      end subroutine sbmclos
!
      subroutine sbmdilp(kn,pt,pq,pps,ptp)
!     Terluna: closure (c)'s parcel from level kn (sbm.diluted_parcel): lifted dry adiabatically with its vapour kept,
!     undiluted, to the condensation level of its air, and above it taking in the surrounding air at dil_bm per metre
!     of rise, each layer crossed in equal steps of log pressure of at most dmix_bm of mixing exponent, each lifting it
!     dry adiabatically, mixing in the air interpolated there in log pressure and condensing the mixture; the metres
!     of a step from the surrounding air's virtual temperature; its temperature by level. Above the first level where
!     it is not buoyant over its level of free convection, which closure (c) does not read, it is held at the
!     surrounding air's temperature
      integer, intent(in) :: kn
      real, intent(in) :: pt(kn),pq(kn),pps
      real, intent(out) :: ptp(kn)
      real :: ztv(NLEV),zbase,zlb,zscale,ztn,zqn,zx1,zx0,zabove,zdzl,zwidth,zxa,zxb,zmixed,zw,zwm,zdz,zmx
      integer :: jl,js,isteps
      logical :: lbuoy
!
      zbase=sbmlcl(pt(kn),pq(kn),zp(kn),zp(1))
      zlb=log(zbase/pps)
      ztv(1:kn)=pt(1:kn)*(1.+(1./rdbrv-1.)*pq(1:kn))
      zscale=gascon/ga
      ztn=pt(kn)
      zqn=pq(kn)
      call sbmsat(ztn,zqn,zp(kn))
      ptp(kn)=ztn
      lbuoy=.false.
      do jl=kn-1,1,-1
       zx1=log(sigma(jl+1))
       zx0=log(sigma(jl))
       zabove=min(zx1,zlb)-zx0
       if(zabove <= 0.) then
        ztn=ztn*(sigma(jl)/sigma(jl+1))**akap
        call sbmsat(ztn,zqn,zp(jl))
       else
        zdzl=zscale*0.5*(ztv(jl)+ztv(jl+1))*zabove
        isteps=max(1,ceiling(dil_bm*zdzl/dmix_bm))
        zwidth=(zx1-zx0)/real(isteps)
        zxa=zx1
        do js=1,isteps
         if(js == isteps) then
          zxb=zx0
         else
          zxb=zx1-real(js)*zwidth
         endif
         ztn=ztn*exp(akap*(zxb-zxa))
         zmixed=min(zxa,zlb)-zxb
         if(zmixed > 0. .and. dil_bm > 0.) then
          zw=(zx1-zxb)/(zx1-zx0)
          zwm=(zx1-0.5*(min(zxa,zlb)+zxb))/(zx1-zx0)
          zdz=zscale*((1.-zwm)*ztv(jl+1)+zwm*ztv(jl))*zmixed
          zmx=1.-exp(-dil_bm*zdz)
          ztn=(1.-zmx)*ztn+zmx*((1.-zw)*pt(jl+1)+zw*pt(jl))
          zqn=(1.-zmx)*zqn+zmx*((1.-zw)*pq(jl+1)+zw*pq(jl))
         endif
         call sbmsat(ztn,zqn,pps*exp(zxb))
         zxa=zxb
        enddo
       endif
       ptp(jl)=ztn
       if(ptp(jl) > pt(jl)) then
        lbuoy=.true.
       elseif(lbuoy) then
        ptp(1:jl-1)=pt(1:jl-1)
        return
       endif
      enddo
      end subroutine sbmdilp
!
      subroutine sbmromps(kn,pt,pq,pps,ptp,ktop,prh)
!     Terluna: closure (c)'s reference relative humidity on the levels 1 to kn (sbm_closure.romps_reference): from
!     the level above the source kn to the parcel's top ktop, Romps's (2014) delta / (delta + gamma), 1 where gamma is
!     not positive and at most 1, gamma = L Gamma / (R_v T^2) - g / (R_d T) from the column's own lapse as
!     numpy.gradient forms it, delta the cloud body's delb_bm plus the organized detrainment -1/2 d ln K / dz where the
!     updraught's kinetic energy K falls, K by the IFS's budget (Cy36r1, Part IV, eqs. 6.9-6.11) with the diluted
!     parcel's virtual buoyancy in steps of at most kesub_bm metres; at the source delb_bm / (delb_bm + gamma)
      integer, intent(in) :: kn,ktop
      real, intent(in) :: pt(kn),pq(kn),pps,ptp(kn)
      real, intent(out) :: prh(kn)
      real :: ztv(NLEV),zz(NLEV),zlapse(NLEV),zgam(NLEV),zbuoy(NLEV),zkin(NLEV)
      real :: zepsv,za,zb,zdx1,zdx2,zkk,zh,zf,zorg,zdel
      integer :: jl,js,isteps
!
      zepsv=1./rdbrv-1.
      ztv(1:kn)=pt(1:kn)*(1.+zepsv*pq(1:kn))
      zz(kn)=0.
      do jl=kn-1,1,-1
       zz(jl)=zz(jl+1)+gascon*0.5*(ztv(jl)+ztv(jl+1))/ga*log(zp(jl+1)/zp(jl))
      enddo
      zlapse(1)=-((pt(2)-pt(1))/(zz(2)-zz(1)))
      zlapse(kn)=-((pt(kn)-pt(kn-1))/(zz(kn)-zz(kn-1)))
      do jl=2,kn-1
       zdx1=zz(jl)-zz(jl-1)
       zdx2=zz(jl+1)-zz(jl)
       zlapse(jl)=-(-zdx2/(zdx1*(zdx1+zdx2))*pt(jl-1)+(zdx2-zdx1)/(zdx1*zdx2)*pt(jl)+zdx1/(zdx2*(zdx1+zdx2))*pt(jl+1))
      enddo
      do jl=1,kn
       zgam(jl)=ALV*zlapse(jl)/(rv_bm*pt(jl)**2)-ga/(gascon*pt(jl))
       zbuoy(jl)=ga*(ptp(jl)*(1.+zepsv*sbmqs(ptp(jl),zp(jl)))-ztv(jl))/ztv(jl)
      enddo
      za=2.*dil_bm*(1.+kedrag_bm)
      zb=1./(keturb_bm*(1.+kevm_bm))
      zkin(kn)=0.5*kew_bm**2
      do jl=kn-1,max(ktop,1),-1
       isteps=max(1,ceiling((zz(jl)-zz(jl+1))/kesub_bm))
       zh=(zz(jl)-zz(jl+1))/real(isteps)
       zkk=zkin(jl+1)
       do js=0,isteps-1
        zf=(real(js)+0.5)/real(isteps)
        zkk=max(zkk+zh*(-za*zkk+zb*((1.-zf)*zbuoy(jl+1)+zf*zbuoy(jl))),1.e-6)
       enddo
       zkin(jl)=zkk
      enddo
      prh(:)=0.
      do jl=kn-1,max(ktop,1),-1
       zorg=max(0.,-0.5*log(zkin(jl)/zkin(jl+1))/(zz(jl)-zz(jl+1)))
       zdel=delb_bm+zorg
       if(zgam(jl) > 0.) then
        prh(jl)=min(1.,zdel/(zdel+zgam(jl)))
       else
        prh(jl)=1.
       endif
      enddo
      if(zgam(kn) > 0.) then
       prh(kn)=min(1.,delb_bm/(delb_bm+zgam(kn)))
      else
       prh(kn)=1.
      endif
      end subroutine sbmromps
!
      real function sbmpbl(kh,pt,pq)
!     Terluna: route L's boundary-layer top (Pa) at grid point kh (sbm_closure.boundary_layer_pressure): the gust
!     patch's depth (gustcol: the bulk Richardson number with the thermals' excess) on the state entering the
!     convection, its winds after the vertical diffusion, with this step's surface fluxes over the lowest level's air
!     density, as a pressure linear in log pressure through the ground and the full levels' heights (massflux.heights);
!     the surface pressure where the surface buoyancy flux is not upward
      integer, intent(in) :: kh
      real, intent(in) :: pt(NLEV),pq(NLEV)
      real :: zu(NLEV),zv(NLEV),zzf(NLEV),ztv(NLEV)
      real :: zrho,zwth0,zwq0,zust,zh,zwst,zzh,zr,ztop,zx0,zx1,zy0,zy1
      integer :: jl
!
      zu(:)=du(kh,1:NLEV)+dudt(kh,1:NLEV)*deltsec2
      zv(:)=dv(kh,1:NLEV)+dvdt(kh,1:NLEV)*deltsec2
      zrho=zp(NLEV)/(gascon*pt(NLEV)*(1.+(1./rdbrv-1.)*pq(NLEV)))
      zwth0=-dshfl(kh)/(zrho*acpd)
      zwq0=-1000.*devap(kh)/zrho
      zust=dust3(kh)**(1./3.)
      call gustcol(NLEV,pt,pq,zu,zv,sigma,sigmah,ga,gascon,akap,rdbrv,zwth0,zwq0,zust,zh,zwst)
      sbmpbl=dp(kh)
      if(.not. zh > 0.) return
      ztv(:)=pt(:)*(1.+(1./rdbrv-1.)*pq(:))
      zzh=0.
      do jl=NLEV,1,-1
       zr=gascon*ztv(jl)/ga
       zzf(jl)=zzh+zr*log(dp(kh)*sigmah(jl)/zp(jl))
       if(jl > 1) then
        ztop=dp(kh)*sigmah(jl-1)
       else
        ztop=zp(1)*0.5
       endif
       zzh=zzh+zr*log(dp(kh)*sigmah(jl)/ztop)
      enddo
      zx0=0.
      zy0=log(dp(kh))
      do jl=NLEV,1,-1
       zx1=zzf(jl)
       zy1=log(zp(jl))
       if(zh < zx1) then
        if(zh == zx0) then
         sbmpbl=exp(zy0)
        else
         sbmpbl=exp((zy1-zy0)/(zx1-zx0)*(zh-zx0)+zy0)
        endif
        return
       endif
       zx0=zx1
       zy0=zy1
      enddo
      sbmpbl=exp(zy0)
      end function sbmpbl
!
      subroutine sbmrlhist
!     Terluna: route L's history (sbm_closure.history_entry, RouteL.record): each level's change of relative
!     humidity this step by the dynamics (PlaSim's adiabatic tendencies of temperature, stt at the grid points, and of
!     vapour, dqt) and the radiation (dtdtswr and dtdtlwr) on the state entering the physics, (dq/dt - rh qs xlv_bm
!     dT/dt / (rv_bm T^2)) / qs over deltsec, into the step's slot of a ring of NRLWIN steps; zdrh24 the ring's sum,
!     the change over the previous rlwin_bm seconds
      real :: zdtt(NHOR,NLEV),zqs,zrh
      integer :: jh,jl,jw,jslot
!
      if(nint(rlwin_bm/deltsec) /= NRLWIN) call mpabort('Terluna route L: rlwin_bm must span NRLWIN time steps')
      call sbmdtt(zdtt)
      jslot=mod(nstep,NRLWIN)+1
      do jl=1,NLEV
       do jh=1,NHOR
        zqs=sbmqs(dt(jh,jl),dp(jh)*sigma(jl))
        zrh=dq(jh,jl)/zqs
        drlhist(jh,jl,jslot)=(dqt(jh,jl)-zrh*zqs*xlv_bm*(dtdtswr(jh,jl)+dtdtlwr(jh,jl)+zdtt(jh,jl))              &
     &                       /(rv_bm*dt(jh,jl)**2))/zqs*deltsec
       enddo
      enddo
      zdrh24(:,:)=0.
      do jw=1,NRLWIN
       zdrh24(:,:)=zdrh24(:,:)+drlhist(:,:,jw)
      enddo
      end subroutine sbmrlhist
"""

# The dynamics' heating at the grid points, a routine of its own so that rainmod.f90 compiles on PlaSim's module alone.
DYNAMICS_HEATING = """
!     =================
!     SUBROUTINE SBMDTT
!     =================

      subroutine sbmdtt(pdtt)
      use rainmod
!
!     Terluna: route L's dynamics' heating (K/s) at the grid points: the spectral tendency stt as the step's spectral
!     part leaves it, transformed as mkdqtgp transforms sqt, with the scaling PlaSim's energy diagnostics give it
!
      real, intent(out) :: pdtt(NHOR,NLEV)
      real :: zstt(NESP,NLEV)
      call mpgallsp(zstt,stt,NLEV)
      zstt(:,:)=zstt(:,:)*ct*ww
      call sp2fl(zstt,pdtt,NLEV)
      call fc2gp(pdtt,NLON,NLPP*NLEV)
      return
      end subroutine sbmdtt
"""

RESTART_GET = f"""!     {CLOSURE_MARKER}: route L's history, empty where a restart does not hold it
      if(nclos_bm == 1 .and. nroutel_bm == 1) then
       allocate(drlhist(NHOR,NLEV,NRLWIN))
       allocate(drllast(NHOR))
       drlhist(:,:,:)=0.
       drllast(:)=-1.e30
       irlread=0
       if(mypid == NROOT) call restart_has_array('dsbmrlhist',irlread)
       call mpbci(irlread)
       if(irlread == 1) then
        call mpgetgp('dsbmrlhist',drlhist,NHOR,NLEV*NRLWIN)
        call mpgetgp('dsbmrllast',drllast,NHOR,1)
       endif
      endif
"""
RESTART_PUT = f"""      if(nclos_bm == 1 .and. nroutel_bm == 1) then   ! {CLOSURE_MARKER}: route L's history
       call mpputgp('dsbmrlhist',drlhist,NHOR,NLEV*NRLWIN)
       call mpputgp('dsbmrllast',drllast,NHOR,1)
      endif
"""

SBM_CLOSURE_EDITS = [
    (_CAP_LINE, _CAP_LINE + _declarations()),
    (_CAP_NAMELIST, _CAP_NAMELIST
     + f"      namelist/rainmod_nl/{','.join(NAMELIST_NAMES[:6])}   ! {CLOSURE_MARKER}\n"
     + f"      namelist/rainmod_nl/{','.join(NAMELIST_NAMES[6:len(CLOSURE_VALUES) + 1])}\n"
     + f"      namelist/rainmod_nl/{','.join(NAMELIST_NAMES[len(CLOSURE_VALUES) + 1:])}\n"),
    (_CAP_BROADCAST, _CAP_BROADCAST + ''.join(f"      call mpbc{'i' if n in INTEGER_NAMES else 'r'}({n})"
                                              + (f"   ! {CLOSURE_MARKER}\n" if i == 0 else "\n")
                                              for i, n in enumerate(NAMELIST_NAMES))),
    (_DEPTH_DECLARATION, _DEPTH_DECLARATION
     + f"      real :: zdrh24(NHOR,NLEV)   ! {CLOSURE_MARKER}: route L's 24-hour change of relative humidity\n"),
    (_BEFORE_LOOP, "      zqe(:,:)=dq(:,1:NLEV)+dqdt(:,1:NLEV)*deltsec2\n"
     f"      if(nclos_bm == 1 .and. nroutel_bm == 1) call sbmrlhist   ! {CLOSURE_MARKER}: route L's history\n"
     "!\n      do jhor=1,NHOR\n       zp(:)=dp(jhor)*sigma(:)\n"),
    (_PARCEL, f"!      {CLOSURE_MARKER}: the column by closure (c) instead, with route L where it is set\n"
              "       if(nclos_bm == 1) then\n        call sbmclosc(jhor)\n        cycle\n       endif\n" + _PARCEL),
    (_LCL_END, _LCL_END + ROUTINES),
    (_MKSBM_END, _MKSBM_END + DYNAMICS_HEATING),
    (_RAININI_END, RESTART_GET + _RAININI_END),
    (_RAINSTOP, _RAINSTOP + RESTART_PUT),
]
