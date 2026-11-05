"""The second redesign's cloud layers only (4 November 2026): a switch on which levels DESIGN2 relaxes, nsubcl_bm in
rainmod's namelist, for one GCM test (SBM_SUBCLOUD_EDITS here, exoplasim_run.SBM_SUBCLOUD_PATCHES over the second
redesign's and the heat rule's patches, rainmod's nsubcl_bm 1 under nclos_bm 3).

DESIGN2 relaxes every level from the trigger's top to the source toward 0.70 of its own saturation, Frierson's (2007)
form inherited from the design scheme, so it rains the boundary layer's vapour above 0.70 and marks those layers as
convective cloud. In R2 that drains the land's sub-cloud layer: with the falling rain's evaporation counted, the
convection takes 1.42 mm/day from it over the lunar day and 2.95 in the afternoon, against the design years' 0.28.
CM1's own convection moistens its land sub-cloud layer by 0.94 mm/day over the lunar day and leaves it nearly
unchanged while it rains: the vapour its resolved transport takes out of the layer, 3.1 mm/day, the falling rain's
evaporation returns, and it forms no condensate below its condensation level
(session-tools/agents/gap_boxes/followups/subcloud/ on the research drive).

With nsubcl_bm 1 ('cloud'), the levels below the condensation level of the lowest level's air (the one the trigger and
the plume take) are left wholly alone: their reference vapour is their own vapour and their reference warming nil, so
the scheme neither dries, moistens nor heats them, and the boundary-layer scheme, the surface and the falling rain's
evaporation carry their budget. The cloud layers relax as in DESIGN2. A column with no cloud layer above its reference
vapour, or none below its reference temperature, forms no rain and goes to the shallow scheme with the cap's top, as
DESIGN2's 'dry' columns do. Where DESIGN2 changes no level below the condensation level the column is DESIGN2's, bit for
bit. Moist enthalpy and water are kept. With nsubcl_bm 0, the default, the routine is DESIGN2's.

Arrays run top first, as PlaSim's levels do.
"""
from __future__ import annotations

import numpy as np

from climate.gcm import sbm
from climate.gcm import sbm_redesign2 as sr2
from climate.gcm import sbm_redesign2_heat as srh

RULES = ('all', 'cloud')                            # nsubcl_bm 0 and 1
SUBCLOUD_MARKER = "Terluna: the second redesign's cloud layers only"
DRIVE = ('session-tools/agents/gap_boxes/followups/subcloud/subcloud_rule.py on the research drive '
         '(rules_subcloud.txt)')


def column(t, q, sigma, dsigma, ps, c, akap, step2=sr2.STEP2_S, tau=sr2.TAU_S, rule='cloud') -> dict:
    """One column's convection by the second redesign with the levels named: 'all' is sbm_redesign2.column itself;
    'cloud' leaves the levels below the condensation level wholly alone and relaxes the rest as DESIGN2 does."""
    if rule not in RULES:
        raise ValueError(f'the sub-cloud rule is one of {RULES}, not {rule!r}')
    if rule == 'all':
        return sr2.column(t, q, sigma, dsigma, ps, c, akap, step2=step2, tau=tau)
    t, q = np.asarray(t, float), np.asarray(q, float)
    sigma, dsigma = np.asarray(sigma, float), np.asarray(dsigma, float)
    n = t.size
    p, m = ps * sigma, ps * dsigma / c['ga']
    out = dict(decision='none', deep=False, dtdt=np.zeros(n), dqdt=np.zeros(n), rain=np.zeros(n), lfc=None, top=None,
               lcl=None, base=None, shallow_top=None, t_ref=None, q_ref=np.full(n, np.nan), trigger_t=None,
               excess=0.0, deficit=0.0, pt=0.0, moist=0.0, pq=0.0, limited=False)
    tp = sbm.parcel(t, q, sigma, p, 0.0, akap, c)
    lcl = sbm.condensation_level(t[-1], q[-1], p[-1], p[0], akap, c)
    out.update(trigger_t=tp, lcl=lcl)
    lfc, top = sr2.cloud_layer(t, tp, p, lcl)
    if top is None:
        return out
    base = sbm.condensation_level(t[-1], q[-1], p[-1], p[top], akap, c)
    out.update(lfc=lfc, top=top, base=base, shallow_top=top)
    if base - p[top] <= sr2.DEEP_PA:
        out['decision'] = 'thin'
        return out
    t_ref = sr2.zbp_reference(t[-1], q[-1], sigma, ps, c, akap, lcl)
    layers = range(top, n)
    q_ref = np.full(n, np.nan)
    for k in layers:
        q_ref[k] = sr2.RH_REF * sbm.saturation(t[k], p[k], c)
        if p[k] > lcl:                                       # below the condensation level: left alone
            q_ref[k] = q[k]
    out.update(t_ref=t_ref, q_ref=q_ref)
    tau = max(tau, step2)
    cpk = c['acpd'] * (1.0 + c['ADV'] * q)
    dtr = {k: (max(t_ref[k] - t[k], 0.0) / tau if p[k] <= lcl else 0.0) for k in layers}
    excess = sum(max(0.0, q[k] - q_ref[k]) * m[k] for k in layers)
    deficit = sum(max(0.0, q_ref[k] - q[k]) * m[k] for k in layers)
    pt = sum(cpk[k] * dtr[k] * m[k] for k in layers)
    out.update(excess=float(excess), deficit=float(deficit), pt=float(pt))
    if excess <= 0.0 or pt <= 0.0:
        out['decision'] = 'dry'
        out['shallow_top'] = next((k for k in range(top, n) if p[k] >= base - sr2.DEEP_PA), n - 1)
        return out
    moist = min(deficit, sr2.F_MOIST * excess)
    rain_amount = excess - moist
    limit = (tau / step2) * pt * tau / c['ALV']
    rain_amount = min(rain_amount, limit)
    rain_total = rain_amount / tau
    dry_share = (rain_amount + moist) / excess
    dqr = np.zeros(n)
    for k in layers:
        d = q_ref[k] - q[k]
        if d < 0.0:
            dqr[k] = d / tau * dry_share
        elif d > 0.0 and deficit > 0.0:
            dqr[k] = d / tau * (moist / deficit)
    scale = c['ALV'] * rain_total / pt
    for k in layers:
        out['dtdt'][k] = dtr[k] * scale
        out['dqdt'][k] = dqr[k]
        if dqr[k] < 0.0:
            out['rain'][k] = -dqr[k] * m[k] * (rain_amount / (rain_amount + moist))
    out.update(decision='deep', deep=True, shallow_top=None, pq=float(rain_total), moist=float(moist),
               limited=bool(limit < excess - moist))
    return out


# ------------------------------------------------------------------------------------------------------- Fortran

VALUES = (('nsubcl_bm', 1, "the second redesign's levels: 0 every level from the top to the source (DESIGN2's), 1 the "
                           'cloud layers only, the levels below the condensation level left alone'),)

_DECLARATION_ANCHOR = srh.declarations()
_DECLARATION = f"      integer :: nsubcl_bm = 0 ! {SUBCLOUD_MARKER} (nclos_bm 3): 0 every level, 1 cloud layers\n"
_NAMELIST = f"      namelist/rainmod_nl/nsubcl_bm   ! {SUBCLOUD_MARKER}\n"
_BROADCAST = f"      call mpbci(nsubcl_bm)   ! {SUBCLOUD_MARKER}\n"
_LOOP_OLD = ('       zqr(jl)=rhref_bm*sbmqs(zt(jl),zp(jl))\n'
             '       zdt(jl)=max(ztr(jl)-zt(jl),0.)/ztau\n')
_LOOP_NEW = _LOOP_OLD + (f'!      {SUBCLOUD_MARKER}: with nsubcl_bm 1 a level below zlcl is left alone\n'
                         '       if(nsubcl_bm == 1 .and. zp(jl) > zlcl) then\n'
                         '        zqr(jl)=zq(jl)\n'
                         '        zdt(jl)=0.\n'
                         '       endif\n')

SBM_SUBCLOUD_EDITS = [
    (_DECLARATION_ANCHOR, _DECLARATION_ANCHOR + _DECLARATION),
    (srh._NAMELIST, srh._NAMELIST + _NAMELIST),
    (srh._BROADCAST, srh._BROADCAST + _BROADCAST),
    (_LOOP_OLD, _LOOP_NEW),
]

DEPARTURES = [
    "the switch is rainmod's nsubcl_bm (0 every level, DESIGN2's; 1 the cloud layers only), read from the namelist, "
    'which the runner sets',
    "the levels left alone are those whose pressure exceeds zlcl, the condensation level the trigger and the plume "
    'take, as the copy takes them',
    'sbmapply marks convective cloud only where the column removes vapour, so no sub-cloud level is marked',
]


def declarations() -> str:
    """The switch's declaration in rainmod's module, at DESIGN2's value 0."""
    return _DECLARATION


def record(rule: str = 'cloud') -> dict:
    """The sub-cloud rule as a run's configuration records it (exoplasim_run.convection_setting) when the switch is
    on."""
    if rule != 'cloud':
        raise ValueError("only the cloud layers are recorded; DESIGN2's own rule leaves the configuration as it was")
    return dict(rule='cloud', switch='nsubcl_bm 1', constants={'nsubcl_bm': 1},
                description='the levels below the condensation level of the lowest level\'s air are left wholly alone '
                            '(their reference vapour their own, their reference warming nil); the cloud layers relax '
                            'as in DESIGN2, and a column with nothing to rain there takes the shallow scheme',
                departures=list(DEPARTURES), source=DRIVE,
                code='climate/gcm/sbm_redesign2_subcloud.py (column; the edit to sbmred2 in Fortran)')
