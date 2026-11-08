"""The redesigned closure of the simplified Betts-Miller scheme, as the Terluna patch puts it into PlaSim
(SBM_REDESIGN_EDITS here, exoplasim_run.SBM_REDESIGN_PATCHES, rainmod's nclos_bm 2), one column at a time in Python,
for its tests (4 November 2026).

The scheme was designed and tested offline in PlaSim's single column under the CM1 boxes' own forcing (sbm_redesign.py
and test_sbm_redesign.py in session-tools/agents/gap_boxes/redesign on the research drive, its form DESIGN, the module's
sha256 DRIVE_SHA256; report.md there and the sections "The redesign, phase A" and "The drop-evaporation check" of
notes_convection.md beside it). This copy carries that form alone, and every constant is closure (c)'s (sbm_closure.py).

1. The trigger is the design scheme's parcel: the lowest level's air lifted undiluted, dry and then moist
   adiabatically, one step a layer (sbm.parcel). Its level of free convection is the lowest level above the source
   where it is buoyant and its top the last level of that buoyant layer; the column convects deep where the cloud from
   the lowest level's condensation level to that top is thicker than DEEP_PA, and a thinner cloud goes to PlaSim's
   shallow scheme with the parcel's top.
2. The reference, over every level from the source to the trigger's top: the temperature of the parcel diluted at
   DILUTION per metre above its condensation level (sbm.diluted_parcel); the vapour, in the buoyant cloud (the levels
   at or above both the condensation level and the trigger's level of free convection), Romps's (2014)
   delta / (delta + gamma) of saturation at that temperature (sbm_closure.romps_reference), and beneath it, in the
   forced cloud, the sub-cloud layer and at the source, the air's own vapour.
3. A layer is warmed toward the reference only where the reference stands above it, and the rest of the rain's latent
   heat is spread by one amount through the relaxed layers, so that the column's moist enthalpy is kept: a free
   troposphere warmer than the diluted parcel is left to the radiation and the circulation.
4. Where the relaxed layers hold more vapour than their reference, the excess rains over TAU_S and the heating releases
   its latent heat; where they hold less, the column moistens without rain by Frierson's (2007) qref form, the
   reference vapour scaled by one factor so that no rain forms and the temperature shifted by one amount so that no
   heat is added, in place of the shallow scheme's call for such a column.

The Fortran (sbmred, called from PlaSim's mksbm under nclos_bm 2 over closure (c)'s patch) takes the trigger parcel
from the design scheme's own lines in mksbm, which with eps_bm 0 and nsub_bm 1, as the runner sets them, are
sbm.parcel's one step a layer; it reuses closure (c)'s sbmlcl, sbmdilp, sbmromps and sbmapply. Closure (c)'s sbmdilp
leaves its parcel unintegrated above the first level where it is not buoyant over its level of free convection, the
levels closure (c) never reads; the redesign reads the diluted parcel up to the undiluted parcel's top, so under
nclos_bm 2 the patch lets sbmdilp rise through every level, as sbm.diluted_parcel does. The rain goes into dprcl, and
the convective cloud layers (icclev, icctot) are marked where the raining column removes vapour, as the design scheme
marks them; a column that moistens forms no rain and marks none. Its heating enters the convection's heating
diagnostic as the raining column's does. nsbmroute holds the decision at each grid point (1 rains, 3 shallow, 4
moistens, 0 none). The tests compare every output with this copy.

Arrays run top first, as PlaSim's levels do.
"""
from __future__ import annotations

import numpy as np

from climate.gcm import sbm
from climate.gcm import sbm_closure as scl

# Closure (c)'s constants, as the redesign takes them.
TAU_S = scl.TAU_S                  # the relaxation time (s)
DILUTION = scl.DILUTION            # the reference parcel's dilution per metre above its condensation level
DEEP_PA = scl.DEEP_PA              # the depth test's threshold (Pa)
MAX_MIX = scl.MAX_MIX              # the diluted parcel's largest mixing exponent in one step
STEP2_S = scl.STEP2_S              # PlaSim's deltsec2, the shortest relaxation time
DECISIONS = ('none', 'thin', 'moisten', 'deep')
DECISION_CODES = dict(none=0, deep=1, thin=3, moisten=4)   # nsbmroute in the Fortran

DRIVE_SHA256 = '031d414bb9d2311fa6cb7aaf2c62abb5f5909cb8fd2521a6c2aa341dc37f6380'
REDESIGN_MARKER = 'Terluna: the redesigned closure'


def cloud_layer(t, tp):
    """A parcel's level of free convection (the lowest level above the source where it is buoyant) and its top (the
    last level of that buoyant layer), or (None, None)."""
    n = t.size
    lfc = next((k for k in range(n - 2, -1, -1) if tp[k] > t[k]), None)
    if lfc is None:
        return None, None
    top = lfc
    for k in range(lfc - 1, -1, -1):
        if tp[k] <= t[k]:
            break
        top = k
    return lfc, top


def reference_vapour(t, q, sigma, ps, c, t_ref, top, base, lfc):
    """The reference vapour (kg/kg) and relative humidity (of saturation at the reference temperature, nan where the
    reference is the air's own vapour) from the top to the source: Romps's delta / (delta + gamma) at the levels at or
    above both the condensation level base (Pa) and the level of free convection lfc, the air's own vapour beneath
    them and at the source."""
    n = t.size
    p = ps * sigma
    q_ref, rh = np.full(n, np.nan), np.full(n, np.nan)
    ref = scl.romps_reference(t, q, sigma, ps, t_ref, top, c)
    for k in range(top, n):
        if k < n - 1 and p[k] <= base and k <= lfc:
            rh[k] = ref['rh'][k]
            q_ref[k] = rh[k] * sbm.saturation(t_ref[k], p[k], c)
        else:
            q_ref[k] = q[k]
    return q_ref, rh


def column(t, q, sigma, dsigma, ps, c, akap, step2=STEP2_S, tau=TAU_S) -> dict:
    """One column's convection by the redesign: the decision ('none': the trigger buoyant nowhere; 'thin': its cloud no
    deeper than DEEP_PA, left to the shallow scheme with the parcel's top; 'moisten': the no-rain branch; 'deep': it
    rains), the tendencies (K/s, kg/kg/s) and the rain formed at each level (kg/m2/s), with the trigger's level of free
    convection and top, the condensation level (Pa), the references, the vapour above the reference over the relaxed
    layers (pq, kg/m2/s once over tau), the heating toward the reference (pt, W/m2), the shift (K/s) and the no-rain
    branch's factor."""
    t, q = np.asarray(t, float), np.asarray(q, float)
    sigma, dsigma = np.asarray(sigma, float), np.asarray(dsigma, float)
    n = t.size
    p, m = ps * sigma, ps * dsigma / c['ga']
    out = dict(decision='none', deep=False, dtdt=np.zeros(n), dqdt=np.zeros(n), rain=np.zeros(n), lfc=None, top=None,
               base=None, shallow_top=None, t_ref=None, q_ref=np.full(n, np.nan), rh_ref=np.full(n, np.nan),
               trigger_t=None, pq=0.0, pt=0.0, shift=0.0, factor=None)
    tp = sbm.parcel(t, q, sigma, ps * sigma, 0.0, akap, c)
    out['trigger_t'] = tp
    lfc, top = cloud_layer(t, tp)
    if top is None:
        return out
    base = sbm.condensation_level(t[-1], q[-1], p[-1], p[top], akap, c)
    out.update(lfc=lfc, top=top, base=base, shallow_top=top)
    if base - p[top] <= DEEP_PA:                             # the depth test: a thin cloud goes to the shallow scheme
        out['decision'] = 'thin'
        return out
    t_ref = sbm.diluted_parcel(t, q, sigma, ps, DILUTION, akap, c, max_mix=MAX_MIX)
    q_ref, rh = reference_vapour(t, q, sigma, ps, c, t_ref, top, base, lfc)
    out.update(t_ref=t_ref, q_ref=q_ref, rh_ref=rh)
    tau = max(tau, step2)
    cpk = c['acpd'] * (1.0 + c['ADV'] * q)
    layers = range(top, n)
    dtr = {k: max(t_ref[k] - t[k], 0.0) / tau for k in layers}   # warmed only toward a reference above the air
    pq = sum((q[k] - q_ref[k]) * m[k] for k in layers) / tau
    pt = sum(cpk[k] * dtr[k] * m[k] for k in layers)
    cpm = sum(cpk[k] * m[k] for k in layers)
    out.update(pq=float(pq), pt=float(pt))
    if pq > 0.0:                                             # deep: the vapour above the reference rains
        dqr = {k: (q_ref[k] - q[k]) / tau for k in layers}
        pq_sum = -sum(dqr[k] * m[k] for k in layers)
        corr = (c['ALV'] * pq_sum - pt) / cpm
        cond = sum(max(0.0, -dqr[k]) * m[k] for k in layers)
        for k in layers:
            out['dtdt'][k] = dtr[k] + corr
            out['dqdt'][k] = dqr[k]
            if dqr[k] < 0.0:
                out['rain'][k] = -dqr[k] * pq_sum / cond * m[k]
        out.update(decision='deep', deep=True, shift=float(corr), shallow_top=None, pq=float(pq_sum))
        return out
    factor = sum(q[k] * m[k] for k in layers) / sum(q_ref[k] * m[k] for k in layers)
    corr = -pt / cpm
    for k in layers:                                         # Frierson's qref form: no net rain, no net heating
        out['dqdt'][k] = (factor * q_ref[k] - q[k]) / tau
        out['dtdt'][k] = dtr[k] + corr
    out.update(decision='moisten', shift=float(corr), factor=float(factor), shallow_top=None)
    return out


def budgets(o, q, sigma, dsigma, ps, c) -> dict:
    """A column's convective budgets: the rain formed (kg/m2/s), the vapour change (kg/m2/s), the heating (W/m2) and the
    moist enthalpy change (W/m2, nil to rounding for every decision)."""
    m = ps * np.asarray(dsigma, float) / c['ga']
    cpk = c['acpd'] * (1.0 + c['ADV'] * np.asarray(q, float))
    heat = float(np.sum(cpk * o['dtdt'] * m))
    vap = float(np.sum(o['dqdt'] * m))
    return dict(rain=float(np.sum(o['rain'])), vapour=vap, heating=heat, enthalpy=heat + c['ALV'] * vap)


def enthalpy_flux(o, q, sigma, dsigma, ps, c) -> np.ndarray:
    """The upward flux of moist enthalpy cp T + L q (W/m2) through the half level beneath each full level implied by
    the tendencies: the enthalpy the layers above it gain (the rain's latent heat counted where the scheme heats, its
    water leaving without enthalpy). Positive is upward; beneath the lowest level it is nil to rounding."""
    m = ps * np.asarray(dsigma, float) / c['ga']
    cpk = c['acpd'] * (1.0 + c['ADV'] * np.asarray(q, float))
    return np.cumsum((cpk * o['dtdt'] + c['ALV'] * o['dqdt']) * m)


# ------------------------------------------------------------------------------------------------------- Fortran

SOURCE = ('sbm_redesign.py (DESIGN), report.md and notes_convection.md ("The redesign, phase A", "The '
          'drop-evaporation check"), session-tools/agents/gap_boxes/redesign on the research drive, 4 November 2026')

# Departures of the Fortran from the Python copy, each recorded with the run.
DEPARTURES = [
    "the trigger parcel is the design scheme's own in mksbm, which with eps_bm 0 and nsub_bm 1 (the runner refuses "
    "the others) is sbm.parcel's one step a layer",
    "the depth test's threshold, the relaxation time and closure (c)'s constants come from rainmod's namelist "
    "(dpdeep_bm, tau_bm and closure (c)'s values), which the runner sets to the copy's",
    "closure (c)'s sbmdilp rises through every level under nclos_bm 2, where closure (c) leaves it unintegrated above "
    "the first level where it is not buoyant over its level of free convection, so that the redesign reads the "
    "diluted parcel up to the undiluted parcel's top as sbm.diluted_parcel gives it",
    "the convective cloud layers (icclev, icctot) are marked where the raining column removes vapour, as the design "
    "scheme marks them; a column that moistens forms no rain and marks none, and the copy has no cloud layers",
    "the moistening column's heating enters the convection's heating diagnostic (zdtc) as the raining column's does",
]


def redesign_record() -> dict:
    """The redesigned closure as a run's configuration records it (exoplasim_run.convection_setting)."""
    return dict(form='redesign', switch='nclos_bm 2', constants={name: value for name, value, _ in scl.CLOSURE_VALUES},
                trigger='the lowest level\'s air lifted undiluted one step a layer, the design scheme\'s parcel; its '
                        'level of free convection and top as the scheme finds them',
                depth='deep only where the cloud from the lowest level\'s condensation level to the trigger\'s top is '
                      'thicker than the depth test\'s threshold; a thinner one takes the shallow scheme with the '
                      'parcel\'s top',
                reference='from the source to the trigger\'s top, the temperature of the parcel diluted at dil_bm per '
                          'metre above its condensation level (sbmdilp); the vapour, at or above both the '
                          'condensation level and the trigger\'s level of free convection, Romps\'s (2014) delta / '
                          '(delta + gamma) of saturation at that temperature as closure (c) forms it (sbmromps), and '
                          'beneath them and at the source the air\'s own',
                energy='a layer warmed only toward a reference above it; the rest of the rain\'s latent heat spread '
                       'by one amount through the relaxed layers, keeping the column\'s moist enthalpy',
                rain='where the relaxed layers hold more vapour than their reference, the excess over tau_s, through '
                     'PlaSim\'s handling of the scheme\'s convective rain',
                moistening='where they hold less, no rain: the reference vapour scaled by one factor so that none '
                           'forms and the temperature shifted by one amount so that no heat is added (Frierson 2007, '
                           'qref), in place of the shallow scheme\'s call',
                relaxation='over tau_s', replaces='the reference humidity rh (rh_bm is not read); no route L',
                departures=list(DEPARTURES), source=SOURCE, drive_module_sha256=DRIVE_SHA256,
                code='climate/gcm/sbm_redesign.py (column; sbmred in Fortran, with closure (c)\'s sbmdilp, sbmromps '
                     'and sbmapply)')


_TRIGGER_FOUND = ("!      its level of free convection, the lowest above the start where it is buoyant, and its top, "
                  "the last level\n")
_HOOK = (f"!      {REDESIGN_MARKER}: under nclos_bm 2 the column takes the redesign on this parcel\n"
         "       if(nclos_bm == 2) then\n        call sbmred(jhor)\n        cycle\n       endif\n")
_DILP_HOLD = "       elseif(lbuoy) then\n        ptp(1:jl-1)=pt(1:jl-1)\n        return\n"
_DILP_FULL = (f"!      {REDESIGN_MARKER}: the redesign (nclos_bm 2) reads this parcel up to the undiluted "
              "parcel's top,\n!      so for it the parcel rises through every level\n"
              "       elseif(lbuoy .and. nclos_bm == 1) then\n        ptp(1:jl-1)=pt(1:jl-1)\n        return\n")
_RLHIST_END = "      end subroutine sbmrlhist\n"

# mksbm's routine for the redesign, after closure (c)'s routines, so that it sees mksbm's variables and calls them.
ROUTINES = """!
      subroutine sbmred(kh)
!     Terluna: the redesigned closure at grid point kh (nclos_bm 2; climate/gcm/sbm_redesign.py, column) on the design
!     scheme's parcel ztp, the lowest level's air lifted undiluted one step a layer: its level of free convection and
!     top, and the depth test from the condensation level of the lowest level's air, a cloud no thicker than dpdeep_bm
!     taking the shallow scheme with the parcel's top. Over a deeper cloud every level from the source to the top
!     relaxes over tau_bm toward the parcel diluted at dil_bm per metre above its condensation level (sbmdilp), a layer
!     warmed only toward a reference above it, and toward the vapour of sbmromps's relative humidity at that
!     temperature in the buoyant cloud (at or above both the condensation level and the level of free convection),
!     the air's own vapour beneath it and at the source. Where the relaxed layers hold more vapour than their
!     reference the excess rains, the heating shifted by one amount to release its latent heat (sbmapply); where they
!     hold less, the reference vapour is scaled by one factor so that no rain forms and the heating shifted so that
!     none is added (Frierson 2007's qref form). nsbmroute: 1 rains, 3 shallow, 4 moistens, 0 none
      integer, intent(in) :: kh
      real :: zt(NLEV),zq(NLEV),ztr(NLEV),zrh(NLEV),zqr(NLEV),zdt(NLEV),zdq(NLEV),zrn(NLEV)
      real :: zbase,zpq,zpt,zcpm,zsum,zcorr,zcond,zqm,zqrm,zfac
      integer :: jl,jlfc,jtop
!
      nsbmroute(kh)=0
      zt(:)=zte(kh,:)
      zq(:)=zqe(kh,:)
      jlfc=0
      do jl=NLEM,1,-1
       if(ztp(jl) > zt(jl)) then
        jlfc=jl
        exit
       endif
      enddo
      if(jlfc == 0) return
      jtop=jlfc
      do jl=jlfc-1,1,-1
       if(ztp(jl) <= zt(jl)) exit
       jtop=jl
      enddo
      zbase=sbmlcl(zt(NLEV),zq(NLEV),zp(NLEV),zp(jtop))
      if(zbase-zp(jtop) <= dpdeep_bm) then
       if(nshallow == 1) then
        ishallow(kh)=1
        ilift(kh)=NLEP
        itop(kh)=jtop
       endif
       nsbmroute(kh)=3
       return
      endif
      call sbmdilp(NLEV,zt,zq,dp(kh),ztr)
      call sbmromps(NLEV,zt,zq,dp(kh),ztr,jtop,zrh)
      zpq=0.
      zpt=0.
      zcpm=0.
      do jl=jtop,NLEV
       if(jl < NLEV .and. zp(jl) <= zbase .and. jl <= jlfc) then
        zqr(jl)=zrh(jl)*sbmqs(ztr(jl),zp(jl))
       else
        zqr(jl)=zq(jl)
       endif
       zdt(jl)=max(ztr(jl)-zt(jl),0.)/ztau
       zpq=zpq+(zq(jl)-zqr(jl))*zm(jl)
       zpt=zpt+zcp(jl)*zdt(jl)*zm(jl)
       zcpm=zcpm+zcp(jl)*zm(jl)
      enddo
      zpq=zpq/ztau
      if(zpq > 0.) then
       zsum=0.
       zcond=0.
       do jl=jtop,NLEV
        zdq(jl)=(zqr(jl)-zq(jl))/ztau
        zsum=zsum+zdq(jl)*zm(jl)
        zcond=zcond+max(0.,-zdq(jl))*zm(jl)
       enddo
       zpq=-zsum
       zcorr=(ALV*zpq-zpt)/zcpm
       zrn(:)=0.
       do jl=jtop,NLEV
        zdt(jl)=zdt(jl)+zcorr
        if(zdq(jl) < 0.) zrn(jl)=-zdq(jl)*zpq/zcond*dsigma(jl)
       enddo
       call sbmapply(kh,jtop,NLEV,zdt,zdq,zrn)
       nsbmroute(kh)=1
       return
      endif
      zqm=0.
      zqrm=0.
      do jl=jtop,NLEV
       zqm=zqm+zq(jl)*zm(jl)
       zqrm=zqrm+zqr(jl)*zm(jl)
      enddo
      zfac=zqm/zqrm
      zcorr=-zpt/zcpm
      do jl=jtop,NLEV
       zdq(jl)=(zfac*zqr(jl)-zq(jl))/ztau
       zdt(jl)=zdt(jl)+zcorr
       dtdt(kh,jl)=dtdt(kh,jl)+zdt(jl)
       dqdt(kh,jl)=dqdt(kh,jl)+zdq(jl)
       zdtc(kh,jl)=zdt(jl)
      enddo
      nsbmroute(kh)=4
      end subroutine sbmred
"""

SBM_REDESIGN_EDITS = [
    (_TRIGGER_FOUND, _HOOK + _TRIGGER_FOUND),
    (_DILP_HOLD, _DILP_FULL),
    (_RLHIST_END, _RLHIST_END + ROUTINES),
]
