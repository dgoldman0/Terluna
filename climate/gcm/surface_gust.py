"""Beljaars's gust in PlaSim's surface transfer, in place of its free-convection term over water and its least wind of
1 m/s (28 October 2026; the record is the section "Beljaars's gust in PlaSim" in README.md beside this file).

PlaSim's surflx (fluxmod.f90) gives the surface its transfer of heat and vapour |U| C_N F_h and of momentum |U| C_N F_m:
|U| the lowest level's wind, at least zumin (1 m/s); C_N = (kappa / ln(z_l / z0 + 1))^2 at that level's height z_l; F the
Louis-ECHAM factors of the bulk Richardson number. In unstable air over water F_h is instead Miller, Beljaars and
Palmer's (1992) interpolation toward the free convection of tank experiments, (1 + (cfree dthv^(1/3) / (|U| C_N))^1.25)^0.8,
whose constant 0.0016 m s-1 K-1/3 carries Earth's gravity and air. Beljaars (1995), a co-author of that interpolation,
took the large eddies of the convective boundary layer instead: they sweep the surface air at about their convective
velocity w*, which joins the mean wind in the bulk formulas. CM1's surface layer takes that form over land and, in the
Terluna builds since 28 October, over water; COARE's ocean algorithm takes it too (Fairall and others 2003). With the
gust (fluxmod's ngust = 1) the wind in surflx's bulk formulas is

    |U|_eff = max( sqrt(|U|^2 + (beta w*)^2), u_min ),   w* = (g / theta_v,n  h  (w'theta_v')_0)^(1/3)

with w* in unstable air and none in stable or neutral air, beta 1 and u_min 0.1 m/s as in CM1's sfclay (its vconvc and
its least wind). The bulk Richardson number takes |U|_eff, and every surface takes Louis's unstable factors, the
free-convection term over water dropped. The surface buoyancy flux (w'theta_v')_0 comes from the surface fluxes of the
step before (dshfl and devap, as pblnl forms them), and the boundary layer's depth h from the nonlocal scheme's rule
(pbl_scheme.nonlocal_column): the bulk Richardson number from the lowest level up, B u*^2 in its shear, zero at h for a
parcel carrying the thermals' excess b (w'theta_v')_0 / w_m. zumin keeps its other use, the least wind shear in vdiff's
diffusivities.

`surface_transfer` is surflx's transfer on one of PlaSim's columns, with or without the gust; `gust_column` the depth and
convective velocity of a column. GUST_PATCHES put the gust into fluxmod.f90 behind ngust, over the nonlocal scheme's
patch (pbl_fortran.PBL_PATCHES), whose constants the depth reads. Arrays run top first, as PlaSim's levels do.
"""
from __future__ import annotations

import numpy as np

from climate.gcm import pbl_fortran as pbf
from climate.gcm import pbl_scheme as ps

GUST_MARKER = "Terluna: Beljaars's gust in the surface transfer"
GUST_BETA = 1.0          # the gust's factor on w*: CM1's vconvc (sfclay.F), Beljaars (1995)
GUST_UMIN = 0.1          # m/s: the least wind in CM1's bulk formulas (sfclay.F)
EXOPLASIM_ZUMIN = 1.0    # fluxmod's zumin, compared there with the squared wind (m2/s2), so 1 m/s
EXOPLASIM_CFREE = 0.0016  # m s-1 K-1/3: Miller, Beljaars and Palmer's constant, fluxmod's cfree


def bulk_depth(u, v, thv, zf, thref: float, ricr: float, ustar: float, g: float, b_ustar: float) -> float:
    """The nonlocal scheme's depth (CAM's pblintd, as pbl_scheme.nonlocal_column takes it): the bulk Richardson number
    between the lowest level and each level above, g (thv_k - thref) (z_k - z_n) / (thv_n (|V_k - V_n|^2 + b u*^2)), from
    the lowest level up; h where it reaches ricr, linear in the number between levels, or the top level's height."""
    nb = len(thv) - 1
    ri_below = 0.0
    for k in range(nb - 1, -1, -1):
        vvk = max((u[k] - u[nb]) ** 2 + (v[k] - v[nb]) ** 2 + b_ustar * ustar * ustar, ps.TINY)
        ri = g * (thv[k] - thref) * (zf[k] - zf[nb]) / (thv[nb] * vvk)
        if ri >= ricr:
            if ri - ri_below <= 0.0:
                return float(zf[k + 1])
            return float(zf[k + 1] + (ricr - ri_below) / (ri - ri_below) * (zf[k] - zf[k + 1]))
        ri_below = ri
    return float(zf[0])


def gust_column(t, q, u, v, sigma, sigmah, wth0: float, wq0: float, ustar: float, c: dict, prm=ps.PARAMETERS) -> dict:
    """One column's boundary-layer depth h (m) and convective velocity w* = (g / thv_n * (w'thv')0 * h)^(1/3) (m/s) in
    unstable air, from its state (top first), its kinematic surface fluxes of potential temperature (K m/s) and vapour
    (kg/kg m/s), upward, and its friction velocity; h and w* are 0 in stable or neutral air. The depth is the nonlocal
    scheme's in unstable air: a first depth from the lowest level's own thv, w_m = (u*^3 + beta_m eps kappa B h)^(1/3)
    from it, and the depth again for a parcel carrying the excess b (w'thv')0 / w_m."""
    c = ps.plasim_constants(c)
    g, gascon, akap, rdbrv = c['ga'], c['gascon'], c['akap'], c['rdbrv']
    t, q, u, v = (np.asarray(a, float) for a in (t, q, u, v))
    sigma, sigmah = np.asarray(sigma, float), np.asarray(sigmah, float)
    nb = t.size - 1
    eps_v = 1.0 / rdbrv - 1.0
    theta = t / sigma ** akap
    thv = theta * (1.0 + eps_v * q)
    zf, _ = ps.heights(t, q, sigma, sigmah, g, gascon, rdbrv)
    wthv0 = wth0 * (1.0 + eps_v * q[nb]) + eps_v * theta[nb] * wq0
    if not wthv0 > 0.0:
        return dict(unstable=False, wthv0=float(wthv0), h=0.0, wstar=0.0)
    buoy = g / thv[nb] * wthv0
    depth = lambda thref: bulk_depth(u, v, thv, zf, thref, prm['ricr_unstable'], ustar, g, prm['b_ustar'])
    h0 = depth(thv[nb])
    wm0 = (ustar ** 3 + prm['beta_m'] * prm['eps_surface'] * ps.KARMAN * buoy * h0) ** (1.0 / 3.0)
    excess = prm['b_excess'] * wthv0 / wm0
    h = depth(thv[nb] + excess)
    return dict(unstable=True, wthv0=float(wthv0), buoyancy=float(buoy), h_first=h0, excess=float(excess), h=h,
                wstar=float((buoy * h) ** (1.0 / 3.0)))


def surface_transfer(col: dict, sigma, sigmah, c: dict, ngust: int = 0, cfree: float = EXOPLASIM_CFREE,
                     zumin: float = EXOPLASIM_ZUMIN, b: float = ps.STABILITY, cc: float = ps.STABILITY,
                     d: float = ps.STABILITY, beta: float = GUST_BETA, umin: float = GUST_UMIN,
                     prm=ps.PARAMETERS) -> dict:
    """surflx's transfer on one column (ntsa 2), step for step. col holds PlaSim's dt and dq (the levels top first, then
    the surface's), du and dv (the levels), dz0 (m), dls (the land share: below 1 the column counts as water), dp (Pa),
    and the surface fluxes of the step before: dshfl (W/m2) and devap (m/s of water), downward, and dust3 (u*^3).
    Returns the squared wind in the bulk formulas (m2/s2), the lowest level's height, the bulk Richardson number, the
    factors F_h and F_m, the transfer coefficients dtransh and dtransm (m/s) and, with the gust, its h and w*."""
    c = ps.plasim_constants(c)
    ga, gascon, akap, rdbrv = c['ga'], c['gascon'], c['akap'], c['rdbrv']
    sigma, sigmah = np.asarray(sigma, float), np.asarray(sigmah, float)
    n = sigma.size
    dt, dq = np.asarray(col['dt'], float), np.asarray(col['dq'], float)
    du, dv = np.asarray(col['du'], float), np.asarray(col['dv'], float)
    eps_v = 1.0 / rdbrv - 1.0
    dtsa = dt[n - 1] * sigma[n - 1] ** (-akap) * (1.0 + eps_v * dq[n - 1])
    wind2 = du[n - 1] * du[n - 1] + dv[n - 1] * dv[n - 1]
    out = dict(h=0.0, wstar=0.0)
    if ngust == 1:
        s = pbf.surface_from_plasim(col['dshfl'], col['devap'], col['dust3'], dt[n], col['dp'], dq[n - 1], c)
        gust = gust_column(dt[:n], dq[:n], du, dv, sigma, sigmah, s['wth0'], s['wq0'], s['ustar'], c, prm)
        out.update(h=gust['h'], wstar=gust['wstar'])
        absu2 = max(umin * umin, wind2 + (beta * gust['wstar']) ** 2)
    else:
        absu2 = max(zumin, wind2)
    znl = -gascon * 0.5 * (dt[n - 1] + dtsa) * np.log(sigma[n - 1]) / ga
    bz0 = znl / col['dz0']
    ri = ga * znl / (absu2 * dtsa) * (dtsa - dt[n] * (1.0 + eps_v * dq[n]))
    kblnz2 = (ps.KARMAN / np.log(bz0 + 1.0)) ** 2
    if ri <= 0.0:
        denom = 1.0 + 3.0 * cc * b * np.sqrt(-1.0 * ri * (bz0 + 1.0)) * kblnz2
        rifm = 1.0 - 2.0 * b * ri / denom
        if col['dls'] < 1.0 and ngust != 1:
            dth = -(dtsa - dt[n] * (1.0 + eps_v * dq[n]))
            rifh = (1.0 + (cfree * dth ** (1.0 / 3.0) / np.sqrt(absu2) / kblnz2) ** 1.25) ** 0.8
        else:
            rifh = 1.0 - 3.0 * b * ri / denom
    else:
        denom = np.sqrt(1 + d * ri)
        rifh = 1.0 / (1.0 + 3.0 * b * ri * denom)
        rifm = 1.0 / (1.0 + 2.0 * b * ri / denom)
    trans = np.sqrt(absu2) * kblnz2
    out.update(absu2=float(absu2), znl=float(znl), ri=float(ri), rifh=float(rifh), rifm=float(rifm),
               dtransh=float(trans * rifh), dtransm=float(trans * rifm), dtsa=float(dtsa))
    return out


# ------------------------------------------------------------------------------------------------------- Fortran

def declarations() -> str:
    return (f'      integer :: ngust = 0 ! {GUST_MARKER}: 1 adds it (gustnl)\n'
            f"      real :: gust_beta = {GUST_BETA!r} ! the gust's factor on w* (CM1's vconvc; Beljaars 1995)\n"
            f'      real :: gust_umin = {GUST_UMIN!r} ! the least wind in the bulk formulas with the gust (m/s; CM1)\n')


NAMELIST_NAMES = ('ngust', 'gust_beta', 'gust_umin')

# The routines, appended after surflx. gustnl: the columns, PlaSim's arrays to gustcol's arguments; gustcol: one column.
GUST_CODE = """
!     =================
!     SUBROUTINE GUSTNL
!     =================

      subroutine gustnl(pgust2)
      use fluxmod
!
!     Terluna: Beljaars's (1995) gust in the surface transfer, called by surflx with ngust = 1: each column's squared
!     gust (gust_beta w*)**2 from its state, its kinematic surface fluxes of the step before (dshfl and devap over the
!     surface air's density ps/(R Ts) and the moist heat capacity mkshfl uses, as pblnl forms them) and its friction
!     velocity (dust3 = u*^3), through gustcol; none in stable or neutral air.
!
      real pgust2(NHOR)
      real zt(NLEV),zq(NLEV),zu(NLEV),zv(NLEV)
      real zrho,zwth0,zwq0,zust,zh,zwst
      integer jhor
!
      do jhor=1,NHOR
       zrho=dp(jhor)/(gascon*dt(jhor,NLEP))
       zwth0=-dshfl(jhor)/(zrho*acpd*AMAX1(1.,1.+ADV*dq(jhor,NLEV)))
       zwq0=-1000.*devap(jhor)/zrho
       zust=dust3(jhor)**(1./3.)
       zt(:)=dt(jhor,1:NLEV)
       zq(:)=dq(jhor,1:NLEV)
       zu(:)=du(jhor,1:NLEV)
       zv(:)=dv(jhor,1:NLEV)
       call gustcol(NLEV,zt,zq,zu,zv,sigma,sigmah,ga,gascon,akap,rdbrv,zwth0,zwq0,zust,zh,zwst)
       pgust2(jhor)=(gust_beta*zwst)**2
      enddo
      return
      end subroutine gustnl

!     ==================
!     SUBROUTINE GUSTCOL
!     ==================

      subroutine gustcol(kn,pt,pq,pu,pv,psig,psigh,pg,pgas,pkap,prdbrv,pwth0,pwq0,pust,ph,pwst)
      use fluxmod
!
!     Terluna: one column's boundary-layer depth ph (m) and convective velocity pwst = (g/thv_n (w'thv')0 h)**(1/3)
!     (m/s) in unstable air, both 0 in stable or neutral air (surface_gust.gust_column step for step). The depth is the
!     nonlocal scheme's (pblcol): the bulk Richardson number from the lowest level up, pbl_bu u*^2 in its shear, at
!     pbl_ricru for a parcel carrying the thermals' excess pbl_b (w'thv')0/wm, wm from a first depth without it.
!
      integer kn
      real pt(kn),pq(kn),pu(kn),pv(kn),psig(kn),psigh(kn)
      real pg,pgas,pkap,prdbrv,pwth0,pwq0,pust,ph,pwst
!
      real zth(kn),zthv(kn),ztv(kn),zzf(kn),zzh(kn)
      real zepsv,zthvn,zwthv0,zbuoy,zust3,zh0,zwm0,zexc
      integer jk
!
      zepsv=1./prdbrv-1.
      zth(:)=pt(:)/psig(:)**pkap
      zthv(:)=zth(:)*(1.+zepsv*pq(:))
      ztv(:)=pt(:)*(1.+zepsv*pq(:))
      zzh(kn)=0.
      do jk=kn,2,-1
       zzh(jk-1)=zzh(jk)+pgas*ztv(jk)/pg*LOG(psigh(jk)/psigh(jk-1))
      enddo
      zzf(:)=zzh(:)+pgas*ztv(:)/pg*LOG(psigh(:)/psig(:))
      zthvn=zthv(kn)
      zwthv0=pwth0*(1.+zepsv*pq(kn))+zepsv*zth(kn)*pwq0
      ph=0.
      pwst=0.
      if(zwthv0 > 0.) then
       zbuoy=pg/zthvn*zwthv0
       zust3=pust**3
       zh0=gustdepth(zthvn)
       zwm0=(zust3+pbl_betam*pbl_eps*vonkarman*zbuoy*zh0)**(1./3.)
       zexc=pbl_b*zwthv0/zwm0
       ph=gustdepth(zthvn+zexc)
       pwst=(zbuoy*ph)**(1./3.)
      endif
      return
!
      contains
!
      real function gustdepth(pthref)
!     the bulk Richardson number from the lowest level up; h linear in the number between levels (CAM's pblintd)
      real pthref
      real zrib,zri,zvvk
      integer jl
      zrib=0.
      do jl=kn-1,1,-1
       zvvk=MAX((pu(jl)-pu(kn))**2+(pv(jl)-pv(kn))**2+pbl_bu*pust*pust,1.e-36)
       zri=pg*(zthv(jl)-pthref)*(zzf(jl)-zzf(kn))/(zthvn*zvvk)
       if(zri >= pbl_ricru) then
        if(zri-zrib <= 0.) then
         gustdepth=zzf(jl+1)
        else
         gustdepth=zzf(jl+1)+(pbl_ricru-zrib)/(zri-zrib)*(zzf(jl)-zzf(jl+1))
        endif
        return
       endif
       zrib=zri
      enddo
      gustdepth=zzf(1)
      end function gustdepth
!
      end subroutine gustcol
"""

CFREE_LINE = pbf.module_edits()[0][0]                          # the free-convection patch's line, the PBL patch's anchor
SURFLX_EDITS = [
    ("      real zabsu2(NHOR)    ! squared wind speed (m2/s2)\n",
     "      real zabsu2(NHOR)    ! squared wind speed (m2/s2)\n"
     "      real zgust2(NHOR)    ! Terluna: the gust's square (m2/s2), with ngust = 1\n"),
    ("      zabsu2(:)=AMAX1(zumin,du(:,NLEV)*du(:,NLEV)+dv(:,NLEV)*dv(:,NLEV))\n",
     f"      if(ngust == 1) then   ! {GUST_MARKER} (gustnl), at least gust_umin\n"
     "       call gustnl(zgust2)\n"
     "       zabsu2(:)=AMAX1(gust_umin*gust_umin,du(:,NLEV)*du(:,NLEV)+dv(:,NLEV)*dv(:,NLEV)+zgust2(:))\n"
     "      else\n"
     "       zabsu2(:)=AMAX1(zumin,du(:,NLEV)*du(:,NLEV)+dv(:,NLEV)*dv(:,NLEV))\n"
     "      endif\n"),
    ("        if(dls(jhor) < 1.) then\n         zdth=-(dtsa(jhor)",
     "        if(dls(jhor) < 1. .and. ngust /= 1) then   ! Terluna: with the gust, Louis's factor as over land\n"
     "         zdth=-(dtsa(jhor)"),
    ("      end subroutine surflx\n", "      end subroutine surflx\n" + GUST_CODE),
]


def module_edits() -> list:
    """The module variables after cfree's line, the namelist's line after the nonlocal scheme's, their broadcasts."""
    return [
        (CFREE_LINE, CFREE_LINE + declarations()),
        (pbf.namelist_text(), pbf.namelist_text()[:-1] + ' &\n     &                ,' + ','.join(NAMELIST_NAMES) + '\n'),
        ("      call mpbcr(cfree)\n",
         "      call mpbcr(cfree)\n      call mpbci(ngust)\n      call mpbcr(gust_beta)\n      call mpbcr(gust_umin)\n"),
    ]


GUST_EDITS = module_edits() + SURFLX_EDITS
GUST_PATCHES = [('fluxmod.f90', GUST_MARKER, GUST_EDITS)]
