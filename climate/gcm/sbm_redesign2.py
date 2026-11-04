"""The second redesign of the simplified Betts-Miller scheme's closure, as the Terluna patch puts it into PlaSim
(SBM_REDESIGN2_EDITS here, exoplasim_run.SBM_REDESIGN2_PATCHES, rainmod's nclos_bm 3), one column at a time in Python,
for its tests (4 November 2026).

The scheme was designed and gated offline in PlaSim's single column under the CM1 boxes' own forcing (sbm_redesign2.py,
its form DESIGN2, and the trigger of post_trigger.py, in session-tools/agents/gap_boxes/redesign2 on the research drive,
their sha256 DRIVE_SHA256 and DRIVE_TRIGGER_SHA256; report.md there and the section "The second redesign, phase A" of
notes_convection.md beside it). This copy carries that form with that trigger alone, and every constant is the drive
module's.

1. The trigger is the design scheme's parcel: the lowest level's air lifted undiluted, dry and then moist
   adiabatically, one step a layer (sbm.parcel). Its level of free convection is the first level above the source, at
   or above the condensation level of the lowest level's air, where the parcel is warmer than the air, and its top the
   last level of that buoyant layer, so a buoyant layer next to the source, below the condensation level, starts no
   cloud. The column convects deep where the cloud from the condensation level to that top is thicker than DEEP_PA; a
   thinner cloud goes to PlaSim's shallow scheme with the parcel's top.
2. The reference temperature is the zero-buoyancy state of a plume from the lowest level's air (zbp_reference),
   integrated in NSUB steps a layer with PlaSim's own condensation and entraining EPS_Z / z per metre (z above the
   ground) of air at RH_REF of saturation. Nothing above the lowest level enters it.
3. The reference vapour is RH_REF of the air's own saturation at every relaxed level, from the top to the source.
4. The vapour the layers above their reference hold over it (E) supplies the rain, and at most F_MOIST of it moistens
   the layers below theirs (D): the rain is E - min(D, F_MOIST E) over TAU_S.
5. The rain's latent heat warms only the layers standing below their reference temperature, each in proportion to its
   shortfall and none past its reference within PlaSim's deltsec2, and the rain is held to what those layers can take.
   Warmer air gets no convective heating, and the column's moist enthalpy and water are kept.
6. A column with no relaxed layer above its reference vapour, or none below its reference temperature, forms no rain
   and goes to the shallow scheme with the shallow cap's top. A column whose parcel is buoyant nowhere at or above the
   condensation level takes neither scheme.

The Fortran (sbmred2 with the plume's sbmzbp, called from PlaSim's mksbm under nclos_bm 3 over the first redesign's and
closure (c)'s patches) takes the trigger parcel from the design scheme's own lines in mksbm, which with eps_bm 0 and
nsub_bm 1, as the runner sets them, are sbm.parcel's one step a layer. It reuses the depth patch's sbmlcl and closure
(c)'s sbmapply, which adds the tendencies, puts the rain into dprcl and marks the convective cloud layers (icclev,
icctot) where the column removes vapour, as the design scheme marks them. nsbmroute holds the decision at each grid
point (DECISION_CODES). The tests compare every output with this copy.

Arrays run top first, as PlaSim's levels do.
"""
from __future__ import annotations

import math

import numpy as np

from climate.gcm import sbm
from climate.gcm import sbm_redesign as srd

TAU_S = 12.0 * 3600.0              # the relaxation time (s)
EPS_Z = 0.5                        # the plume's entrainment EPS_Z / z per metre, z its height above the ground
RH_REF = 0.70                      # the reference relative humidity of the air's own saturation; the entrained air's
NSUB = 16                          # the plume's steps a layer above the condensation level
F_MOIST = 0.1                      # the largest share of the vapour removed that moistens the layers below theirs
DEEP_PA = 20000.0                  # the depth test's threshold (Pa)
STEP2_S = 3600.0                   # PlaSim's deltsec2 at the runs' 30-minute step, the shortest relaxation time
DECISIONS = ('none', 'thin', 'dry', 'deep')
DECISION_CODES = dict(none=0, deep=1, thin=3, dry=5)       # nsbmroute in the Fortran

DRIVE_SHA256 = '393a09f5e09c186a1200265f3005e4d55c26f80e9cfb843e18ab5e0da66fbdc1'
DRIVE_TRIGGER_SHA256 = '7b975fabc877972d9c2804f067848fa05ef1aa57fe4548b3096af166738b3e26'
REDESIGN2_MARKER = 'Terluna: the second redesign'


def zbp_reference(t_src, q_src, sigma, ps, c, akap, base, eps=EPS_Z, rh=RH_REF, nsub=NSUB) -> np.ndarray:
    """The zero-buoyancy state of a plume from the lowest level's air (temperature t_src, vapour q_src at sigma[-1]),
    its temperature by level, top first. Below base, the condensation level of that air (Pa), the plume rises dry
    adiabatically with its vapour kept and undiluted, and a layer wholly below it is crossed in one step, as sbm.parcel
    crosses it. Every other layer is crossed in nsub equal steps of log pressure, each lifting the plume dry
    adiabatically to the step's end and condensing it (sbm.condense), then mixing in a share 1 - exp(-e dz) of air at
    the plume's own temperature and rh of saturation: dz the step's metres above base, e = eps / z with z the plume's
    height above the ground at the step's middle, both from the plume's own virtual temperature. Mixing at equal
    temperature lowers the plume's vapour and leaves its temperature; the next step's lift and condensation carry the
    deficit. With eps 0 this is the undiluted parcel integrated in nsub steps a layer (sbm.diluted_parcel at dilution 0
    with min_steps nsub). The result reads nothing of the column above its lowest level."""
    sigma = np.asarray(sigma, float)
    n = sigma.size
    p = ps * sigma
    lb = math.log(base / ps)
    scale = c['gascon'] / c['ga']
    epsv = 1.0 / c['rdbrv'] - 1.0
    tp = np.zeros(n)
    tn, qn = sbm.condense(t_src, q_src, p[-1], c)
    tp[-1] = tn
    z = scale * tn * (1.0 + epsv * qn) * math.log(1.0 / sigma[-1])     # the source's height above the ground
    for k in range(n - 2, -1, -1):
        x1, x0 = math.log(sigma[k + 1]), math.log(sigma[k])
        if min(x1, lb) - x0 <= 0.0:                     # wholly below the base: one dry step, as sbm.parcel
            t_old = tn
            tn *= (sigma[k] / sigma[k + 1]) ** akap
            tn, qn = sbm.condense(tn, qn, p[k], c)
            z += scale * 0.5 * (t_old + tn) * (1.0 + epsv * qn) * (x1 - x0)
            tp[k] = tn
            continue
        width = (x1 - x0) / nsub
        xa = x1
        for s in range(1, nsub + 1):
            xb = x0 if s == nsub else x1 - s * width
            t_old = tn
            tn *= math.exp(akap * (xb - xa))
            tn, qn = sbm.condense(tn, qn, ps * math.exp(xb), c)
            dz_step = scale * 0.5 * (t_old + tn) * (1.0 + epsv * qn) * (xa - xb)
            mixed = min(xa, lb) - xb
            if mixed > 0.0 and eps > 0.0:
                dz = dz_step * mixed / (xa - xb)
                e = eps / (z + 0.5 * dz_step)
                mix = 1.0 - math.exp(-e * dz)
                qn = (1.0 - mix) * qn + mix * rh * sbm.saturation(tn, ps * math.exp(xb), c)
            z += dz_step
            xa = xb
        tp[k] = tn
    return tp


def cloud_layer(t, tp, p, base):
    """The trigger's level of free convection, the first level above the source at or above the condensation level
    base (Pa) where the parcel tp is warmer than the air t, and its top, the last level of that buoyant layer; or
    (None, None)."""
    n = t.size
    lfc = next((k for k in range(n - 2, -1, -1) if tp[k] > t[k] and p[k] <= base), None)
    if lfc is None:
        return None, None
    top = lfc
    for k in range(lfc - 1, -1, -1):
        if tp[k] <= t[k]:
            break
        top = k
    return lfc, top


def column(t, q, sigma, dsigma, ps, c, akap, step2=STEP2_S, tau=TAU_S) -> dict:
    """One column's convection by the second redesign: the decision ('none': the parcel buoyant nowhere at or above the
    condensation level; 'thin': its cloud no deeper than DEEP_PA, left to the shallow scheme with the parcel's top;
    'dry': no relaxed layer above its reference vapour or none below its reference temperature, left to the shallow
    scheme with the shallow cap's top; 'deep': it rains), the tendencies (K/s, kg/kg/s) and the rain formed at each
    level (kg/m2/s), with the trigger's level of free convection and top, the condensation levels (Pa; lcl the one the
    trigger and the plume take, bracketed by the top level's pressure, base the depth test's, bracketed by the
    trigger's top), the references, the excess E and deficit D of vapour over the relaxed layers (kg/m2), the heating
    toward the reference pt (W/m2), the moistening (kg/m2 over tau), the rain (pq, kg/m2/s) and whether the heat's
    limit held the rain."""
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
    lfc, top = cloud_layer(t, tp, p, lcl)
    if top is None:
        return out
    base = sbm.condensation_level(t[-1], q[-1], p[-1], p[top], akap, c)
    out.update(lfc=lfc, top=top, base=base, shallow_top=top)
    if base - p[top] <= DEEP_PA:                             # the depth test: a thin cloud goes to the shallow scheme
        out['decision'] = 'thin'
        return out
    t_ref = zbp_reference(t[-1], q[-1], sigma, ps, c, akap, lcl)
    layers = range(top, n)
    q_ref = np.full(n, np.nan)
    for k in layers:
        q_ref[k] = RH_REF * sbm.saturation(t[k], p[k], c)
    out.update(t_ref=t_ref, q_ref=q_ref)
    tau = max(tau, step2)
    cpk = c['acpd'] * (1.0 + c['ADV'] * q)
    dtr = {k: max(t_ref[k] - t[k], 0.0) / tau for k in layers}              # warmed only toward a reference above
    excess = sum(max(0.0, q[k] - q_ref[k]) * m[k] for k in layers)          # kg/m2
    deficit = sum(max(0.0, q_ref[k] - q[k]) * m[k] for k in layers)
    pt = sum(cpk[k] * dtr[k] * m[k] for k in layers)
    out.update(excess=float(excess), deficit=float(deficit), pt=float(pt))
    if excess <= 0.0 or pt <= 0.0:                           # nothing to rain, or no layer the plume could warm
        out['decision'] = 'dry'
        out['shallow_top'] = next((k for k in range(top, n) if p[k] >= base - DEEP_PA), n - 1)
        return out
    moist = min(deficit, F_MOIST * excess)                   # kg/m2 moved into the layers below their reference
    rain_amount = excess - moist                             # kg/m2 over tau, at least (1 - F_MOIST) E here
    limit = (tau / step2) * pt * tau / c['ALV']              # what the layers below their reference can take
    rain_amount = min(rain_amount, limit)                    # within step2, with none carried past its reference
    rain_total = rain_amount / tau                           # kg/m2/s, positive since E and pt are
    dry_share = (rain_amount + moist) / excess
    dqr = np.zeros(n)
    for k in layers:
        d = q_ref[k] - q[k]
        if d < 0.0:
            dqr[k] = d / tau * dry_share                     # the layers above their reference give rain and moistening
        elif d > 0.0 and deficit > 0.0:
            dqr[k] = d / tau * (moist / deficit)             # those below get their share of the moistening
    scale = c['ALV'] * rain_total / pt
    for k in layers:
        out['dtdt'][k] = dtr[k] * scale
        out['dqdt'][k] = dqr[k]
        if dqr[k] < 0.0:
            out['rain'][k] = -dqr[k] * m[k] * (rain_amount / (rain_amount + moist))
    out.update(decision='deep', deep=True, shallow_top=None, pq=float(rain_total), moist=float(moist),
               limited=bool(limit < excess - moist))
    return out


budgets = srd.budgets
enthalpy_flux = srd.enthalpy_flux


# ------------------------------------------------------------------------------------------------------- Fortran

SOURCE = ('sbm_redesign2.py (DESIGN2) with the trigger of post_trigger.py, report.md and notes_convection.md ("The '
          'second redesign, phase A"), session-tools/agents/gap_boxes/redesign2 on the research drive, 4 November 2026')

# The namelist's values: (Fortran name, value, what it is). The relaxation time and the depth test's threshold are the
# scheme's own tau_bm and dpdeep_bm, which the runner sets to TAU_S and DEEP_PA.
VALUES = (
    ('nzbp_bm', NSUB, "its reference plume's steps a layer"),
    ('epsz_bm', EPS_Z, "the plume's entrainment epsz_bm / z per metre, z its height above the ground"),
    ('rhref_bm', RH_REF, "the reference relative humidity, of the air's own saturation, and the entrained air's"),
    ('fmoist_bm', F_MOIST, 'the largest share of the vapour removed that moistens the layers below their reference'),
)
INTEGER_NAMES = ('nzbp_bm',)

# Departures of the Fortran from the Python copy, each recorded with the run.
DEPARTURES = [
    "the trigger parcel is the design scheme's own in mksbm, which with eps_bm 0 and nsub_bm 1 (the runner refuses "
    "the others) is sbm.parcel's one step a layer",
    "the relaxation time, the depth test's threshold and the form's constants come from rainmod's namelist (tau_bm, "
    "dpdeep_bm, nzbp_bm, epsz_bm, rhref_bm and fmoist_bm), which the runner sets to the copy's",
    "the rain goes into dprcl in mkrain's units, the vapour removed times the rain's share of it times dsigma, where "
    "the copy gives kg/m2/s from the layer's mass; the two agree to rounding",
    "the convective cloud layers (icclev, icctot) are marked where the raining column removes vapour, as the design "
    "scheme marks them, the source and sub-cloud layers included; the copy has no cloud layers",
    "a thin or dry column goes to PlaSim's shallow scheme with the top the copy names, and a column with no level of "
    "free convection to neither scheme",
    "the raining column's heating enters the convection's heating diagnostic (zdtc)",
]


def record() -> dict:
    """The second redesign as a run's configuration records it (exoplasim_run.convection_setting)."""
    return dict(form='redesign2', switch='nclos_bm 3', constants={name: value for name, value, _ in VALUES},
                trigger='the lowest level\'s air lifted undiluted one step a layer, the design scheme\'s parcel; its '
                        'level of free convection the first level above the source, at or above the condensation '
                        'level of that air, where it is warmer than the air, and its top the last level of that '
                        'buoyant layer',
                depth='deep only where the cloud from the condensation level to the trigger\'s top is thicker than '
                      'the depth test\'s threshold; a thinner one takes the shallow scheme with the parcel\'s top',
                reference='every level from the source to the top toward the zero-buoyancy state of a plume from the '
                          'lowest level\'s air, entraining epsz_bm / z per metre (z above the ground) of air at '
                          'rhref_bm of saturation, integrated in nzbp_bm steps a layer with PlaSim\'s condensation '
                          '(sbmzbp); nothing above the lowest level enters it',
                humidity='rhref_bm of the air\'s own saturation at every relaxed level',
                rain='the vapour the layers above their reference hold over it, less at most fmoist_bm of it that '
                     'moistens the layers below theirs, over tau_s, through PlaSim\'s handling of the scheme\'s '
                     'convective rain',
                energy='the rain\'s latent heat warms only the layers below their reference temperature, in '
                       'proportion to the shortfall and none past it within deltsec2, the rain held to what they can '
                       'take; warmer air gets no convective heating',
                routing='a column with no layer above its reference vapour or none below its reference temperature '
                        'takes the shallow scheme with the shallow cap\'s top; one with no level of free convection '
                        'neither scheme',
                relaxation='over tau_s', replaces='the reference humidity rh (rh_bm is not read); no route L',
                departures=list(DEPARTURES), source=SOURCE, drive_module_sha256=DRIVE_SHA256,
                drive_trigger_sha256=DRIVE_TRIGGER_SHA256,
                code='climate/gcm/sbm_redesign2.py (column; sbmred2 and sbmzbp in Fortran, with the depth patch\'s '
                     'sbmlcl and closure (c)\'s sbmapply)')


def declarations() -> str:
    """The form's namelist variables in rainmod's module, at the copy's values."""
    lines = []
    for i, (name, value, what) in enumerate(VALUES):
        kind = 'integer' if name in INTEGER_NAMES else 'real'
        literal = str(int(value)) if name in INTEGER_NAMES else repr(float(value))
        note = f'{REDESIGN2_MARKER} (nclos_bm 3): {what}' if i == 0 else what
        lines.append(f'      {kind} :: {name} = {literal} ! {note}\n')
    return ''.join(lines)


def namelist_values() -> dict:
    """The rainmod namelist entries (name, text) for the form's constants, as the runner writes them."""
    return {name.upper(): str(int(value)) if name in INTEGER_NAMES else repr(float(value)) for name, value, _ in VALUES}


_SOURCE_LINE = '      integer :: nsbmsrc(NHOR) = 0 ! the level its deep convection rises from\n'
_ROUTE_L_NAMELIST = '      namelist/rainmod_nl/nroutel_bm,rlsig1_bm,rlsig2_bm,rlmem_bm,rlwin_bm,xlv_bm\n'
_LAST_BROADCAST = '      call mpbcr(xlv_bm)\n'
_NAMELIST = f"      namelist/rainmod_nl/{','.join(name for name, _, _ in VALUES)}   ! {REDESIGN2_MARKER}\n"
_BROADCASTS = ''.join(f"      call mpbc{'i' if name in INTEGER_NAMES else 'r'}({name})"
                      + (f'   ! {REDESIGN2_MARKER}\n' if i == 0 else '\n') for i, (name, _, _) in enumerate(VALUES))
_HOOK = (f"!      {REDESIGN2_MARKER}: under nclos_bm 3 the column takes the second redesign on this parcel\n"
         "       if(nclos_bm == 3) then\n        call sbmred2(jhor)\n        cycle\n       endif\n")
_SBMRED_END = '      end subroutine sbmred\n'

# mksbm's routines for the second redesign, after the first redesign's, so that they see mksbm's variables and call
# its routines.
ROUTINES = f"""!
      subroutine sbmred2(kh)
!     {REDESIGN2_MARKER} of the closure at grid point kh (nclos_bm 3; climate/gcm/sbm_redesign2.py, column), on
!     the design scheme's parcel ztp, the lowest level's air lifted undiluted one step a layer. Its level of free
!     convection is the first level above the source, at or above the condensation level of the lowest level's air,
!     where it is warmer than the air, and its top the last level of that buoyant layer; a cloud from the condensation
!     level to the top no thicker than dpdeep_bm takes the shallow scheme with the parcel's top. Over a deeper cloud
!     every level from the source to the top relaxes over tau_bm, the temperature toward the zero-buoyancy state of a
!     plume from the lowest level's air (sbmzbp) and the vapour toward rhref_bm of the air's own saturation. The vapour
!     the layers above their reference hold over it rains, less at most fmoist_bm of it that moistens the layers below
!     theirs; its latent heat warms only the layers below their reference temperature, in proportion to the shortfall,
!     and the rain is held to what they can take within deltsec2 (sbmapply adds it all). A column with no layer above
!     its reference vapour or none below its reference temperature takes the shallow scheme with the shallow cap's top.
!     nsbmroute: 1 rains, 3 shallow (thin), 5 shallow (dry), 0 none
      integer, intent(in) :: kh
      real :: zt(NLEV),zq(NLEV),ztr(NLEV),zqr(NLEV),zdt(NLEV),zdq(NLEV),zrn(NLEV)
      real :: zlcl,zbase,zex,zdef,zpt,zmst,zram,zlim,zrtot,zdry,zscl,zfr,zd
      integer :: jl,jlfc,jtop
!
      nsbmroute(kh)=0
      zt(:)=zte(kh,:)
      zq(:)=zqe(kh,:)
      zlcl=sbmlcl(zt(NLEV),zq(NLEV),zp(NLEV),zp(1))
      jlfc=0
      do jl=NLEM,1,-1
       if(ztp(jl) > zt(jl) .and. zp(jl) <= zlcl) then
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
      call sbmzbp(zt(NLEV),zq(NLEV),dp(kh),zlcl,ztr)
      zex=0.
      zdef=0.
      zpt=0.
      do jl=jtop,NLEV
       zqr(jl)=rhref_bm*sbmqs(zt(jl),zp(jl))
       zdt(jl)=max(ztr(jl)-zt(jl),0.)/ztau
       zex=zex+max(0.,zq(jl)-zqr(jl))*zm(jl)
       zdef=zdef+max(0.,zqr(jl)-zq(jl))*zm(jl)
       zpt=zpt+zcp(jl)*zdt(jl)*zm(jl)
      enddo
      if(zex <= 0. .or. zpt <= 0.) then
       if(nshallow == 1) then
        ishallow(kh)=1
        ilift(kh)=NLEP
        do jl=jtop,NLEV
         if(zp(jl) >= zbase-dpdeep_bm) exit
        enddo
        itop(kh)=min(jl,NLEV)
       endif
       nsbmroute(kh)=5
       return
      endif
      zmst=min(zdef,fmoist_bm*zex)
      zram=zex-zmst
      zlim=(ztau/deltsec2)*zpt*ztau/ALV
      zram=min(zram,zlim)
      zrtot=zram/ztau
      zdry=(zram+zmst)/zex
      zscl=ALV*zrtot/zpt
      zfr=zram/(zram+zmst)
      zrn(:)=0.
      do jl=jtop,NLEV
       zd=zqr(jl)-zq(jl)
       zdq(jl)=0.
       if(zd < 0.) then
        zdq(jl)=zd/ztau*zdry
       elseif(zd > 0. .and. zdef > 0.) then
        zdq(jl)=zd/ztau*(zmst/zdef)
       endif
       zdt(jl)=zdt(jl)*zscl
       if(zdq(jl) < 0.) zrn(jl)=-zdq(jl)*zfr*dsigma(jl)
      enddo
      call sbmapply(kh,jtop,NLEV,zdt,zdq,zrn)
      nsbmroute(kh)=1
      end subroutine sbmred2
!
      subroutine sbmzbp(pt0,pq0,pps,plcl,ptr)
!     {REDESIGN2_MARKER}'s reference temperature by level (sbm_redesign2.zbp_reference): the zero-buoyancy
!     state of a plume from the lowest level's air (pt0, pq0). Below its condensation level plcl it rises dry
!     adiabatically with its vapour kept, undiluted, and a layer wholly below plcl is crossed in one step, as the
!     scheme's parcel crosses it. Every other layer is crossed in nzbp_bm equal steps of log pressure, each lifting the
!     plume dry adiabatically to its end and condensing it (sbmsat), then mixing in a share 1 - exp(-e dz) of air at
!     the plume's own temperature and rhref_bm of saturation, dz the step's metres above plcl and e = epsz_bm / z, z the
!     plume's height above the ground at the step's middle, the metres from the plume's own virtual temperature. It
!     reads nothing of the column above its lowest level
      real, intent(in) :: pt0,pq0,pps,plcl
      real, intent(out) :: ptr(NLEV)
      real :: zlb,zscale,zepsv,ztn,zqn,zz,zx1,zx0,zwidth,zxa,zxb,zto,zdzs,zmixed,zdz,ze,zmx
      integer :: jl,js
!
      zlb=log(plcl/pps)
      zscale=gascon/ga
      zepsv=1./rdbrv-1.
      ztn=pt0
      zqn=pq0
      call sbmsat(ztn,zqn,zp(NLEV))
      ptr(NLEV)=ztn
      zz=zscale*ztn*(1.+zepsv*zqn)*log(1./sigma(NLEV))
      do jl=NLEM,1,-1
       zx1=log(sigma(jl+1))
       zx0=log(sigma(jl))
       if(min(zx1,zlb)-zx0 <= 0.) then
        zto=ztn
        ztn=ztn*(sigma(jl)/sigma(jl+1))**akap
        call sbmsat(ztn,zqn,zp(jl))
        zz=zz+zscale*0.5*(zto+ztn)*(1.+zepsv*zqn)*(zx1-zx0)
        ptr(jl)=ztn
        cycle
       endif
       zwidth=(zx1-zx0)/real(nzbp_bm)
       zxa=zx1
       do js=1,nzbp_bm
        if(js == nzbp_bm) then
         zxb=zx0
        else
         zxb=zx1-real(js)*zwidth
        endif
        zto=ztn
        ztn=ztn*exp(akap*(zxb-zxa))
        call sbmsat(ztn,zqn,pps*exp(zxb))
        zdzs=zscale*0.5*(zto+ztn)*(1.+zepsv*zqn)*(zxa-zxb)
        zmixed=min(zxa,zlb)-zxb
        if(zmixed > 0. .and. epsz_bm > 0.) then
         zdz=zdzs*zmixed/(zxa-zxb)
         ze=epsz_bm/(zz+0.5*zdzs)
         zmx=1.-exp(-ze*zdz)
         zqn=(1.-zmx)*zqn+zmx*rhref_bm*sbmqs(ztn,pps*exp(zxb))
        endif
        zz=zz+zdzs
        zxa=zxb
       enddo
       ptr(jl)=ztn
      enddo
      end subroutine sbmzbp
"""

SBM_REDESIGN2_EDITS = [
    (_SOURCE_LINE, _SOURCE_LINE + declarations()),
    (_ROUTE_L_NAMELIST, _ROUTE_L_NAMELIST + _NAMELIST),
    (_LAST_BROADCAST, _LAST_BROADCAST + _BROADCASTS),
    (srd._HOOK, srd._HOOK + _HOOK),
    (_SBMRED_END, _SBMRED_END + ROUTINES),
]
