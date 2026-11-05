"""The second redesign's heat beyond its reference (4 November 2026): a switch on DESIGN2's energy rule, nheat_bm in
rainmod's namelist, for one GCM test (SBM_HEAT_EDITS here, exoplasim_run.SBM_HEAT_PATCHES over the second redesign's
patch, rainmod's nheat_bm 1 under nclos_bm 3).

DESIGN2's rule 5 (sbm_redesign2.column) lets the rain's latent heat warm only the layers standing below their reference
temperature, each in proportion to its shortfall and none past its reference within deltsec2, and holds the rain to what
those layers can take: rain <= (tau / deltsec2) pt tau / L. On the GCM's columns that limit binds on two thirds or more
of the deep columns at the three gap-test boxes, and on land it holds the afternoon's rain at once to about a quarter of
what the vapour's relaxation gives (session-tools/agents/gap_boxes/followups/report.md and heat_limit/ on the research
drive). CM1's own convection, raining, puts 0.40-0.66 of its gross heating (microphysics, resolved transport and
turbulence) on land and at the dry box into layers standing at or above DESIGN2's reference, where its net tendency is
small because the imposed circulation cools them (heat_limit/cm1_layers.json): the rule matches CM1's net tendency
there, which in a GCM the dynamics supply, and not the gross heating the scheme itself must supply.

With nheat_bm 1 ('spread'), where the limit would bind the rain is not held. The layers below their reference are warmed
at the limit's own rate, each reaching its reference in deltsec2, and the rest of the rain's latent heat is spread over
the cloud layers, the relaxed levels at or above the condensation level of the lowest level's air, by one temperature
increment a second. Where the limit does not bind the column is DESIGN2's, bit for bit, and the decisions are DESIGN2's
everywhere: a column with no layer below its reference still forms no rain. Moist enthalpy and water are kept. With
nheat_bm 0, the default, the routine is DESIGN2's.

Arrays run top first, as PlaSim's levels do.
"""
from __future__ import annotations

import numpy as np

from climate.gcm import sbm_redesign2 as sr2

RULES = ('reference', 'spread')                     # nheat_bm 0 and 1
HEAT_MARKER = "Terluna: the second redesign's heat beyond its reference"
DRIVE = 'session-tools/agents/gap_boxes/followups/heat_limit/heat_rule.py on the research drive (rules_heat_limit.txt)'


def column(t, q, sigma, dsigma, ps, c, akap, step2=sr2.STEP2_S, tau=sr2.TAU_S, rule='spread') -> dict:
    """One column's convection by the second redesign with the energy rule named: 'reference' is sbm_redesign2.column
    itself; 'spread' leaves its decisions and, where its limit on the rain does not bind, everything as it is, and where
    the limit binds rains the vapour's relaxation unheld, warms the layers below their reference at the limit's rate and
    spreads the rest of the latent heat over the cloud layers by one increment (K/s). Adds 'spread', the increment."""
    if rule not in RULES:
        raise ValueError(f'the heat rule is one of {RULES}, not {rule!r}')
    o = sr2.column(t, q, sigma, dsigma, ps, c, akap, step2=step2, tau=tau)
    o['spread'] = 0.0
    if rule == 'reference' or o['decision'] != 'deep' or not o['limited']:
        return o
    t, q = np.asarray(t, float), np.asarray(q, float)
    sigma, dsigma = np.asarray(sigma, float), np.asarray(dsigma, float)
    n = t.size
    p, m = ps * sigma, ps * dsigma / c['ga']
    tau = max(tau, step2)
    layers = range(o['top'], n)
    cpk = c['acpd'] * (1.0 + c['ADV'] * q)
    t_ref, q_ref = o['t_ref'], o['q_ref']
    excess, deficit, pt = o['excess'], o['deficit'], o['pt']
    moist = min(deficit, sr2.F_MOIST * excess)
    rain_amount = excess - moist                              # kg/m2 over tau: the vapour's relaxation, not held
    rain_total = rain_amount / tau
    dry_share = (rain_amount + moist) / excess
    scale = tau / step2                                       # the toward part at the limit's own rate
    cm = 0.0
    for k in layers:
        if p[k] <= o['lcl']:
            cm += cpk[k] * m[k]
    shift = (c['ALV'] * rain_total - pt * scale) / cm         # K/s over the cloud layers
    dtdt, dqdt, rain = np.zeros(n), np.zeros(n), np.zeros(n)
    for k in layers:
        d = q_ref[k] - q[k]
        if d < 0.0:
            dqdt[k] = d / tau * dry_share
        elif d > 0.0 and deficit > 0.0:
            dqdt[k] = d / tau * (moist / deficit)
        dtdt[k] = max(t_ref[k] - t[k], 0.0) / tau * scale
        if p[k] <= o['lcl']:
            dtdt[k] = dtdt[k] + shift
        if dqdt[k] < 0.0:
            rain[k] = -dqdt[k] * m[k] * (rain_amount / (rain_amount + moist))
    o.update(dtdt=dtdt, dqdt=dqdt, rain=rain, pq=float(rain_total), moist=float(moist), spread=float(shift),
             limited=False)
    return o


# ------------------------------------------------------------------------------------------------------- Fortran

# The namelist's switch: (Fortran name, value of the run that takes it, what it is).
VALUES = (('nheat_bm', 1, "the second redesign's heat rule: 0 DESIGN2's, heat only toward the reference with the rain "
                          "held; 1 the rest of the heat spread over the cloud layers, the rain not held"),)

_DECLARATION_ANCHOR = sr2.declarations().splitlines(keepends=True)[-1]
_DECLARATION = f"      integer :: nheat_bm = 0 ! {HEAT_MARKER} (nclos_bm 3): 0 DESIGN2's, 1 spread\n"
_NAMELIST = f"      namelist/rainmod_nl/nheat_bm   ! {HEAT_MARKER}\n"
_BROADCAST = f"      call mpbci(nheat_bm)   ! {HEAT_MARKER}\n"
_LOCALS_OLD = '      real :: zlcl,zbase,zex,zdef,zpt,zmst,zram,zlim,zrtot,zdry,zscl,zfr,zd\n'
_LOCALS_NEW = '      real :: zlcl,zbase,zex,zdef,zpt,zmst,zram,zlim,zrtot,zdry,zscl,zfr,zd\n      real :: zspr,zcm\n'
_HOLD_OLD = '      zram=min(zram,zlim)\n'
_HOLD_NEW = (f"!     {HEAT_MARKER}: with nheat_bm 1 the rain is not held where the limit binds\n"
             "      zspr=0.\n"
             "      if(nheat_bm /= 1 .or. zram <= zlim) zram=min(zram,zlim)\n")
_SCALE_OLD = '      zscl=ALV*zrtot/zpt\n'
_SCALE_NEW = ('      zscl=ALV*zrtot/zpt\n'
              '      if(nheat_bm == 1 .and. zram > zlim) then\n'
              '       zscl=ztau/deltsec2\n'
              '       zcm=0.\n'
              '       do jl=jtop,NLEV\n'
              '        if(zp(jl) <= zlcl) zcm=zcm+zcp(jl)*zm(jl)\n'
              '       enddo\n'
              '       zspr=(ALV*zrtot-zpt*zscl)/zcm\n'
              '      endif\n')
_HEAT_OLD = '       zdt(jl)=zdt(jl)*zscl\n'
_HEAT_NEW = '       zdt(jl)=zdt(jl)*zscl\n       if(nheat_bm == 1 .and. zp(jl) <= zlcl) zdt(jl)=zdt(jl)+zspr\n'

SBM_HEAT_EDITS = [
    (_DECLARATION_ANCHOR, _DECLARATION_ANCHOR + _DECLARATION),
    (sr2._NAMELIST, sr2._NAMELIST + _NAMELIST),
    (sr2._BROADCASTS, sr2._BROADCASTS + _BROADCAST),
    (_LOCALS_OLD, _LOCALS_NEW),
    (_HOLD_OLD, _HOLD_NEW),
    (_SCALE_OLD, _SCALE_NEW),
    (_HEAT_OLD, _HEAT_NEW),
]

DEPARTURES = [
    "the switch is rainmod's nheat_bm (0 DESIGN2's rule, 1 the spread), read from the namelist, which the runner sets",
    "the cloud layers are the levels from the trigger's top to the lowest whose pressure is at or below the "
    "condensation level the trigger and the plume take (zlcl), as the copy takes them",
    "the spread enters the convection's heating diagnostic (zdtc) with the rest of the heating, through sbmapply",
]


def declarations() -> str:
    """The switch's declaration in rainmod's module, at DESIGN2's value 0."""
    return _DECLARATION


def record(rule: str = 'spread') -> dict:
    """The heat rule as a run's configuration records it (exoplasim_run.convection_setting) when the switch is on."""
    if rule != 'spread':
        raise ValueError("only the spread is recorded; DESIGN2's own rule leaves the configuration as it was")
    return dict(rule='spread', switch='nheat_bm 1', constants={'nheat_bm': 1},
                description="where DESIGN2's limit on the rain would bind, the rain is the vapour's relaxation "
                            'unheld; the layers below their reference are warmed at the limit\'s rate, each reaching '
                            'it in deltsec2, and the rest of the latent heat is spread over the cloud layers (at or '
                            'above the condensation level) by one increment; elsewhere and in every decision the '
                            'column is DESIGN2\'s',
                departures=list(DEPARTURES), source=DRIVE,
                code='climate/gcm/sbm_redesign2_heat.py (column; the edits to sbmred2 in Fortran)')
