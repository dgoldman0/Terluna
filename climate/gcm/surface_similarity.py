"""CM1's surface-layer similarity over water in PlaSim's surface transfer, beside Beljaars's gust (29 October 2026; the
record is the section "Similarity over water in PlaSim" in README.md beside this file).

PlaSim's surflx (fluxmod.f90) gives each surface its transfer of heat and vapour |U| C_N F_h and of momentum |U| C_N F_m,
C_N = (kappa / ln(z_l / z0 + 1))^2 at its lowest level's height z_l and F the Louis-ECHAM factors of the bulk Richardson
number; with Beljaars's gust (surface_gust.py, ngust 1) |U| carries the boundary layer's convective velocity. Set beside
Monin-Obukhov similarity integrated to the same height on the same air and wind, that transfer evaporated 0.72 of the
similarity's over the Near-side Sea box's air and 0.82 over the GCM's own (rules and outcome in session-tools/evaporation
on the research drive, 29 October), most of the difference in light wind over warm water, where similarity's
free-convection functions and its roughness lengths for heat and vapour carry more than Louis's factors and one
Charnock roughness. With nsimw 1 (fluxmod), over ice-free water the transfer comes from CM1's own surface layer
(sfclay.F of the Terluna builds, sfcmodel 2 with isftcflx 0), carried to the lowest level's height:

- the bulk Richardson number Ri_b = g/theta z_l (thv_l - thv_s) / U^2 sets the regime: psi = -10 ln(z_l/z0) at least
  -10 for Ri_b >= 0.2, -5 Ri_b ln(z_l/z0) / (1.1 - 5 Ri_b) at least -10 below it, nil at Ri_b = 0, and in unstable air
  Paulson's functions of z_l/L from CM1's tables at steps of 0.01, at most 0.9 ln(z_l/z0), z_l/L from the step before's
  u* and temperature scale; air that was unstable the step before stays out of the stable regimes;
- Charnock's roughness 0.0185 u*^2/g at the planet's gravity with the smooth-flow term 0.11 x 1.5e-5/u*, at most
  2.85 mm, and Fairall's roughness for heat and vapour 5.5e-5 (u* z0/nu)^-0.6 within 2e-9 and 1e-4 m,
  psi_q = max(ln((z_l + z0q)/z0q) - psi_h, 2);
- u* = (u*_old + kappa U/(ln(z_l/z0) - psi_m)) / 2 as CM1 averages it, the temperature scale kappa (theta_l - theta_s)/psi_q,
  and the roughness from the new u*, each carried to the next step as CM1 carries UST, MOL and ZNT;
- the transfer velocities kappa u*/psi_q for heat and vapour and u*^2/U for momentum, U the gust's wind at least 0.1 m/s.

These enter surflx as its stability factors against the neutral |U| C_N, so mkstress, mkshfl, mkevap, the 10 m wind and
the debug output take them unchanged. Land, sea ice and partly frozen sea keep PlaSim's transfer. CM1's subgrid wind
(Mahrt and Sun's, 0.19 m/s at its 6 km grid) is left out. `similarity_column` is one step of one column (simwcol in
Fortran); `surface_transfer` is surflx's transfer with it; SIMW_PATCHES put it into fluxmod.f90 behind nsimw, over the
gust's patch.
"""
from __future__ import annotations

import math

import numpy as np

from climate.gcm import pbl_scheme as ps
from climate.gcm import surface_gust as sg

SIMW_MARKER = 'Terluna: similarity over water'
KARMAN = 0.4                     # CM1's karman and PlaSim's vonkarman
CHARNOCK = 0.0185                # sfclay's CZO
SMOOTH = 0.11 * 1.5e-5           # sfclay's smooth-flow term (m2/s), over u*
Z0_MAX = 2.85e-3                 # m: sfclay's cap on the roughness over water
Z0Q_COEFF, Z0Q_EXPONENT = 5.5e-5, -0.60    # Fairall's roughness for heat and vapour (m), on the roughness Reynolds number
Z0Q_MIN, Z0Q_MAX = 2.0e-9, 1.0e-4
PSI_FLOOR = 2.0                  # sfclay's least psi_q
ZOL_MIN = -9.9999                # its range of z/L in unstable air
WIND_MIN = 0.1                   # m/s: its least wind
START = dict(ust=0.1, mol=0.0, z0=1.0e-4)   # the state before a column's first step


def tables():
    """CM1's sfclayinit: Paulson's functions for momentum and heat at z/L = -0.01 n, n = 0..1000."""
    zoln = -np.arange(1001, dtype=float) * 0.01
    x = (1.0 - 16.0 * zoln) ** 0.25
    psim = 2.0 * np.log(0.5 * (1.0 + x)) + np.log(0.5 * (1.0 + x * x)) - 2.0 * np.arctan(x) + 2.0 * np.arctan(1.0)
    psih = 2.0 * np.log(0.5 * (1.0 + (1.0 - 16.0 * zoln) ** 0.5))
    return psim, psih


PSIM_TABLE, PSIH_TABLE = tables()


def similarity_column(za: float, tha: float, thg: float, thva: float, thvg: float, ta: float, absu: float, g: float,
                      state: dict = START) -> dict:
    """One step of CM1's surface layer over water on one column (simwcol, line for line): the lowest level's height za
    (m), its potential and virtual potential temperature tha and thva against the surface's thg and thvg (K, on one
    reference pressure), its temperature ta (K), the wind in the bulk formulas absu (m/s) and gravity g, from the state
    of the step before (u*, the temperature scale and the roughness). Returns the transfer velocities of heat and vapour
    (transh) and of momentum (transm), the parts, and the state for the next step."""
    ust, mol, znt = float(state['ust']), float(state['mol']), float(state['z0'])
    gz1oz0 = math.log(za / znt)
    wspd = max(absu, WIND_MIN)
    br = g / tha * za * (thva - thvg) / (wspd * wspd)
    if mol < 0.0:
        br = min(br, 0.0)
    if br >= 0.2:
        regime = 1
        psim = max(-10.0 * gz1oz0, -10.0)
        psih = psim
    elif br > 0.0:
        regime = 2
        psim = max(-5.0 * br * gz1oz0 / (1.1 - 5.0 * br), -10.0)
        psih = psim
    elif br == 0.0:
        regime = 3
        psim = 0.0
        psih = 0.0
    else:
        regime = 4
        zol = br * gz1oz0 if ust < 0.01 else KARMAN * g / tha * za * mol / (ust * ust)
        zol = max(min(zol, 0.0), ZOL_MIN)
        n = int(-zol * 100.0)
        r = -zol * 100.0 - n
        psim = PSIM_TABLE[n] + r * (PSIM_TABLE[n + 1] - PSIM_TABLE[n])
        psih = PSIH_TABLE[n] + r * (PSIH_TABLE[n + 1] - PSIH_TABLE[n])
        psih = min(psih, 0.9 * gz1oz0)
        psim = min(psim, 0.9 * gz1oz0)
    psix = gz1oz0 - psim
    visc = (1.32 + 0.009 * (ta - 273.15)) * 1.0e-5
    restar = ust * znt / visc
    z0q = min(max(Z0Q_COEFF * restar ** Z0Q_EXPONENT, Z0Q_MIN), Z0Q_MAX)
    psiq = max(math.log((za + z0q) / z0q) - psih, PSI_FLOOR)
    ust = 0.5 * ust + 0.5 * KARMAN * wspd / psix
    mol = KARMAN * (tha - thg) / psiq
    znt = min(CHARNOCK * ust * ust / g + SMOOTH / ust, Z0_MAX)
    return dict(transh=ust * KARMAN / psiq, transm=ust * ust / wspd, regime=regime, br=br, psim=float(psim),
                psih=float(psih), psix=psix, psiq=psiq, z0q=z0q, wspd=wspd, state=dict(ust=ust, mol=mol, z0=znt))


def applies(col: dict) -> bool:
    """Ice-free water: PlaSim's land share below 1 and no ice cover."""
    return col['dls'] < 1.0 and col.get('dicec', 0.0) <= 0.0


def surface_transfer(col: dict, sigma, sigmah, c: dict, state: dict = START, nsimw: int = 1, **gust) -> dict:
    """surflx's transfer with the gust (surface_gust.surface_transfer, ngust 1) and, with nsimw 1 over ice-free water,
    CM1's similarity in its place, entered as the stability factors against PlaSim's neutral transfer; col as
    surface_gust takes it, with dicec its ice cover (0 if absent)."""
    out = sg.surface_transfer(col, sigma, sigmah, c, ngust=1, **gust)
    if nsimw != 1 or not applies(col):
        return out
    cc = ps.plasim_constants(c)
    sigma = np.asarray(sigma, float)
    n = sigma.size
    dt, dq = np.asarray(col['dt'], float), np.asarray(col['dq'], float)
    eps_v = 1.0 / cc['rdbrv'] - 1.0
    tha = dt[n - 1] * sigma[n - 1] ** (-cc['akap'])
    s = similarity_column(out['znl'], tha, dt[n], tha * (1.0 + eps_v * dq[n - 1]), dt[n] * (1.0 + eps_v * dq[n]),
                          dt[n - 1], math.sqrt(out['absu2']), cc['ga'], state)
    ztrans = math.sqrt(out['absu2']) * (KARMAN / math.log(out['znl'] / col['dz0'] + 1.0)) ** 2
    rifh, rifm = s['transh'] / ztrans, s['transm'] / ztrans
    out.update(rifh=rifh, rifm=rifm, dtransh=ztrans * rifh, dtransm=ztrans * rifm, similarity=s)
    return out


# ------------------------------------------------------------------------------------------------------- Fortran

def declarations() -> str:
    return (f'      integer :: nsimw = 0 ! {SIMW_MARKER}: 1 takes CM1\'s surface layer over ice-free water (simwcol)\n'
            f"      real :: dsimust(NHOR) = {START['ust']!r} ! its u* (m/s), carried as CM1 carries UST\n"
            f"      real :: dsimmol(NHOR) = {START['mol']!r} ! its temperature scale (K), CM1's MOL\n"
            f"      real :: dsimz0(NHOR) = {START['z0']!r} ! its roughness (m), CM1's ZNT\n"
            '      real :: simw_psim(0:1000), simw_psih(0:1000) ! Paulson\'s functions at z/L = -0.01 n (sfclayinit)\n'
            '      integer :: nsimwinit = 0\n')


NAMELIST_NAMES = ('nsimw',)

SIMW_CODE = f"""
!     ===================
!     SUBROUTINE SIMWCOL
!     ===================

      subroutine simwcol(kh,pza,pabsu,ptrans,prifh,prifm)
      use fluxmod
!
!     Terluna: CM1's surface layer over water (sfclay.F, sfcmodel 2, isftcflx 0) on column kh at the lowest level's
!     height pza, the wind in the bulk formulas pabsu (surface_similarity.similarity_column line for line): from the
!     step before's u*, temperature scale and roughness (dsimust, dsimmol, dsimz0), the transfer velocities of heat and
!     vapour (u* k/psi_q) and of momentum (u*^2/U), returned as the factors prifh, prifm against PlaSim's neutral
!     transfer ptrans, and the state carried to the next step.
!
      integer kh,jn
      real pza,pabsu,ptrans,prifh,prifm
      real zust,zmol,znt,zgz,zws,zbr,zpsim,zpsih,zzol,zr,zpsix,zvisc,zre,zz0q,zpsiq
      real ztha,zthg,zthva,zthvg,zepsv
!
      if(nsimwinit == 0) then
       do jn=0,1000
        zr=-REAL(jn)*0.01
        zzol=(1.-16.*zr)**0.25
        simw_psim(jn)=2.*LOG(0.5*(1.+zzol))+LOG(0.5*(1.+zzol*zzol))-2.*ATAN(zzol)+2.*ATAN(1.)
        simw_psih(jn)=2.*LOG(0.5*(1.+(1.-16.*zr)**0.5))
       enddo
       nsimwinit=1
      endif
      zepsv=1./rdbrv-1.
      ztha=dt(kh,NLEV)*sigma(NLEV)**(-akap)
      zthg=dt(kh,NLEP)
      zthva=ztha*(1.+zepsv*dq(kh,NLEV))
      zthvg=zthg*(1.+zepsv*dq(kh,NLEP))
      zust=dsimust(kh)
      zmol=dsimmol(kh)
      znt=dsimz0(kh)
      zgz=LOG(pza/znt)
      zws=MAX(pabsu,{WIND_MIN!r})
      zbr=ga/ztha*pza*(zthva-zthvg)/(zws*zws)
      if(zmol < 0.) zbr=MIN(zbr,0.)
      if(zbr >= 0.2) then
       zpsim=MAX(-10.*zgz,-10.)
       zpsih=zpsim
      elseif(zbr > 0.) then
       zpsim=MAX(-5.*zbr*zgz/(1.1-5.*zbr),-10.)
       zpsih=zpsim
      elseif(zbr == 0.) then
       zpsim=0.
       zpsih=0.
      else
       if(zust < 0.01) then
        zzol=zbr*zgz
       else
        zzol={KARMAN!r}*ga/ztha*pza*zmol/(zust*zust)
       endif
       zzol=MAX(MIN(zzol,0.),{ZOL_MIN!r})
       jn=INT(-zzol*100.)
       zr=-zzol*100.-REAL(jn)
       zpsim=simw_psim(jn)+zr*(simw_psim(jn+1)-simw_psim(jn))
       zpsih=simw_psih(jn)+zr*(simw_psih(jn+1)-simw_psih(jn))
       zpsih=MIN(zpsih,0.9*zgz)
       zpsim=MIN(zpsim,0.9*zgz)
      endif
      zpsix=zgz-zpsim
      zvisc=(1.32+0.009*(dt(kh,NLEV)-273.15))*1.e-5
      zre=zust*znt/zvisc
      zz0q=MIN(MAX({Z0Q_COEFF!r}*zre**({Z0Q_EXPONENT!r}),{Z0Q_MIN!r}),{Z0Q_MAX!r})
      zpsiq=MAX(LOG((pza+zz0q)/zz0q)-zpsih,{PSI_FLOOR!r})
      zust=0.5*zust+0.5*{KARMAN!r}*zws/zpsix
      zmol={KARMAN!r}*(ztha-zthg)/zpsiq
      znt=MIN({CHARNOCK!r}*zust*zust/ga+{SMOOTH!r}/zust,{Z0_MAX!r})
      dsimust(kh)=zust
      dsimmol(kh)=zmol
      dsimz0(kh)=znt
      prifh=zust*{KARMAN!r}/zpsiq/ptrans
      prifm=zust*zust/zws/ptrans
      return
      end subroutine simwcol
"""

GUST_LAST_DECLARATION = sg.declarations().splitlines(keepends=True)[-1]
GUST_NAMELIST_END = ',' + ','.join(sg.NAMELIST_NAMES) + '\n'
SIMW_EDITS = [
    (GUST_LAST_DECLARATION, GUST_LAST_DECLARATION + declarations()),
    (GUST_NAMELIST_END, GUST_NAMELIST_END[:-1] + ',' + ','.join(NAMELIST_NAMES) + '\n'),
    ("      call mpbcr(gust_umin)\n", "      call mpbcr(gust_umin)\n      call mpbci(nsimw)\n"),
    ("       dtransh(jhor)=ztrans*zrifh(jhor)\n",
     f"       if(nsimw == 1 .and. dls(jhor) < 1. .and. dicec(jhor) <= 0.) then   ! {SIMW_MARKER} (simwcol)\n"
     "        call simwcol(jhor,znl(jhor),SQRT(zabsu2(jhor)),ztrans,zrifh(jhor),zrifm(jhor))\n"
     "       endif\n"
     "       dtransh(jhor)=ztrans*zrifh(jhor)\n"),
    ("      end subroutine gustcol\n", "      end subroutine gustcol\n" + SIMW_CODE),
]
SIMW_PATCHES = [('fluxmod.f90', SIMW_MARKER, SIMW_EDITS)]
