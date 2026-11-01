"""Cloud water for PlaSim's large-scale condensation: the condensate held in each grid point's layer above freezing,
evaporated where the air is unsaturated, and turned to rain by autoconversion and accretion.

    climate/gcm/.venv/bin/python -m climate.gcm.condensate        # the fall-speed fit from the lunar fall law

PlaSim's large-scale condensation (rainmod.f90, mklsp) turns each layer's supersaturation into precipitation in the same
step, mkrain evaporates part of it on the way down, and the model holds no cloud water between steps; the cloud water
its radiation sees is diagnosed (mkclouds). Over the Near-side Sea after midnight CM1's sea box condenses 1.14 mm/day
at 5.2-11.8 km and re-evaporates 0.85 of it from cloud in the same layers, where PlaSim with its deep convection held off
rains all it condenses there (night_water.py). With the switch (ncwat 1, in PlaSim's own namelist so that the restart
can carry the water, PLASIM_SWITCH_PATCHES in exoplasim_run.py):
- mklsp hands what condenses in a layer at or above freezing, by the layer's temperature after condensing, to the
  layer's cloud water instead of to the precipitation;
- mkcwat, after mklsp and before mkrain, goes down each column. Cloud water in a layer below freezing goes to the
  precipitation at once. Where the layer's air is unsaturated the cloud water evaporates toward saturation, with mklsp's
  correction for the latent heat, at most the water there. What remains turns to rain by Khairoutdinov and Kogan's
  (2000) autoconversion, 1350 qc^2.47 Nc^-1.79, and accretion, 67 (qc qr)^1.15 (kg/kg/s; qc and qr in kg/kg, Nc per
  cm3), as CM1's Morrison scheme takes them (module_mp_morr_two_moment, PRC and PRA), with the droplet number of the
  CM1 boxes, 100 per cm3; qr is the rain formed above the layer over the air's density and the rain's fall speed;
- what the cloud hands on joins mklsp's precipitation for mkrain, which carries it down, evaporates part of it below
  and turns it to snow where the air is below freezing.
The water is neither advected nor mixed: it stays in the layer and grid point where it formed. Radiation keeps PlaSim's
diagnosed clouds. The restart carries the water (dqcw), and a run that starts from a restart without it starts with none.

The rain's fall speed. Accretion needs the rain's mixing ratio, qr = F/(rho V) for a flux F. The lunar fall law
(atmosphere/rain/results/rain_fall_law.json: the single-drop solver's mass-weighted mean speed of an exponential spectrum
in the Moon's air, CM1's rain law) with Marshall and Palmer's intercept, 8e6 per m4, gives V = 0.79, 0.99 and 1.26 m/s
at 1, 3 and 10 mm/day. A power law V = V0 (F/F0)^b, fitted by least squares in logs from 0.3 to 30 mm/day, gives
V0 = 0.99 m/s at F0 = 3 mm/day and b = 0.20, within 1.4% (fall_speed_fit). The intercept is Earth's, an assumption.

Time stepping, as mklsp's: each rate brings the step's predicted state (the leapfrog's state two steps ahead, dq +
deltsec2 dqdt) toward saturation, and the cloud water, a grid-point field of one time level, changes by each rate times
one step, so that over the leapfrog's alternate steps the vapour and the cloud water exchange the same water.

The Python copy (mklsp, mkcwat) follows the Fortran line by line on columns (rows) and levels (top first), for the
tests; the stochastic rain option (nstorain) is left out.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FALL_LAW = ROOT / 'atmosphere' / 'rain' / 'results' / 'rain_fall_law.json'
MARKER = 'Terluna: cloud water'
DAY = 86400.0
N0_MARSHALL_PALMER = 8.0e6        # per m4: Marshall and Palmer's (1948) intercept, Earth's
FIT_RANGE_MM_DAY = (0.3, 30.0)
FR0_MM_DAY = 3.0
PARAMETERS = dict(ncw_cm3=100.0,                       # the CM1 boxes' droplet number (case droplets_cm3)
                  vtr0_cw=0.99,                        # m/s, the rain's mass-weighted fall speed at fr0_cw
                  fr0_cw=FR0_MM_DAY / DAY,             # kg/m2/s, 3 mm/day
                  vtrexp_cw=0.20)                      # the fall speed's exponent on the flux
AUTO = (1350.0, 2.47, -1.79)      # Khairoutdinov and Kogan (2000), as CM1's Morrison scheme takes them
ACCRETION = (67.0, 1.15)
QC_AUTO_MIN, QC_ACC_MIN = 1.0e-6, 1.0e-8               # Morrison's thresholds for the two (kg/kg)


# ------------------------------------------------------------------------------------------------ the fall speed

def fall_speed_fit(path: Path = FALL_LAW, n0: float = N0_MARSHALL_PALMER, low_high=FIT_RANGE_MM_DAY,
                   f0_mm_day: float = FR0_MM_DAY) -> dict:
    """The rain's mass-weighted fall speed against its flux from the lunar fall law's table, for exponential spectra
    with intercept n0: the flux of a spectrum with slope lambda is rho_w pi n0 / lambda^4 times its mass-weighted
    speed, both from the table; fitted as V0 (F/F0)^b in logs over the fluxes low_high (mm/day)."""
    d = json.loads(Path(path).read_text())
    t = d['table']
    lnl = t['ln_lambda_first'] + t['ln_lambda_step'] * np.arange(t['points'])
    lam = np.exp(lnl)
    v = np.exp(np.asarray(t['ln_mass_weighted'], float))
    flux = 1000.0 * math.pi * n0 / lam ** 4 * v                          # kg/m2/s
    mm = flux * DAY
    keep = (mm >= low_high[0]) & (mm <= low_high[1])
    b, a = np.polyfit(np.log(mm[keep] / f0_mm_day), np.log(v[keep]), 1)
    fit = np.exp(a) * (mm[keep] / f0_mm_day) ** b
    at = lambda m: float(np.exp(np.interp(math.log(m), np.log(mm[::-1]), np.log(v[::-1]))))
    return dict(v0_m_s=float(np.exp(a)), f0_mm_day=f0_mm_day, exponent=float(b), n0_per_m4=n0,
                range_mm_day=list(low_high), largest_relative_error=float(np.max(np.abs(fit / v[keep] - 1.0))),
                speeds_m_s={f'{m:g} mm/day': at(m) for m in (1.0, 3.0, 10.0)},
                source=str(Path(path).relative_to(ROOT)) if Path(path).is_relative_to(ROOT) else str(path))


# ------------------------------------------------------------------------------------------- the Python copy

def saturation(t, ps, sigma, c):
    """mklsp's saturation specific humidity at temperature t (K) and pressure ps sigma."""
    q = c['rdbrv'] * c['ra1'] * np.exp(c['ra2'] * (t - c['TMELT']) / (t - c['ra4'])) / (ps * sigma)
    q = np.minimum(q, c['rdbrv'])
    cor = 1.0 / (1.0 - (1.0 / c['rdbrv'] - 1.0) * q)
    return q * cor, cor


def mklsp(s: dict, c: dict, ncwat: int = 1) -> None:
    """PlaSim's mklsp on the state s (columns by levels, top first; updated in place): the condensation of each layer's
    predicted supersaturation in two iterations, its heating, and its precipitation dprll, which with ncwat 1 goes to
    the cloud water dqcw where the layer stands at or above freezing after condensing (dcwcond the rate it hands on)."""
    nlev = s['dq'].shape[1]
    s['dprll'][:] = 0.0
    s['dprsll'][:] = 0.0
    if ncwat == 1:
        s['dcwcond'][:] = 0.0
    for k in range(nlev):
        zq = s['dq'][:, k] + s['dqdt'][:, k] * c['deltsec2']
        zt = s['dt'][:, k] + s['dtdt'][:, k] * c['deltsec2']
        zlcpe = c['ALV'] / (c['acpd'] * (1.0 + c['ADV'] * s['dq'][:, k]))
        zqsat, zcor = saturation(zt, s['dp'], s['sigma'][k], c)
        zqsat = np.where(zqsat < 0.0, 0.0, zqsat)
        m = zq > zqsat
        if m.any():
            zdqdt = (zqsat - zq) / (1.0 + zlcpe * c['ra2'] * (c['TMELT'] - c['ra4']) * zqsat * zcor / (zt - c['ra4']) ** 2)
            ztn = zt - zdqdt * zlcpe
            zqn = zq + zdqdt
            zqsat2, zcor2 = saturation(ztn, s['dp'], s['sigma'][k], c)
            zdqdt = (zqsat2 - zqn) / (1.0 + zlcpe * c['ra2'] * (c['TMELT'] - c['ra4']) * zqsat2 * zcor2
                                      / (ztn - c['ra4']) ** 2)
            ztn = ztn - zdqdt * zlcpe
            zqn = zqn + zdqdt
            zdqdt = (zqn - zq) / c['deltsec2']
            zdtdt = (ztn - zt) / c['deltsec2']
            s['dqdt'][:, k] = np.where(m, s['dqdt'][:, k] + zdqdt, s['dqdt'][:, k])
            s['dtdt'][:, k] = np.where(m, s['dtdt'][:, k] + zdtdt, s['dtdt'][:, k])
            s['dprll'][:, k] = np.where(m, -np.minimum(zdqdt, 0.0) * s['dsigma'][k], s['dprll'][:, k])
        if ncwat == 1:
            warm = (s['dprll'][:, k] > 0.0) & (s['dt'][:, k] + s['dtdt'][:, k] * c['deltsec2'] >= c['TMELT'])
            s['dcwcond'][:, k] = np.where(warm, s['dprll'][:, k] / s['dsigma'][k], s['dcwcond'][:, k])
            s['dqcw'][:, k] = np.where(warm, s['dqcw'][:, k] + s['dcwcond'][:, k] * c['deltsec'], s['dqcw'][:, k])
            s['dprll'][:, k] = np.where(warm, 0.0, s['dprll'][:, k])


def mkcwat(s: dict, c: dict, prm: dict = PARAMETERS) -> dict:
    """The cloud water's step after mklsp (state s updated in place): by level from the top, the water below freezing
    to the precipitation, the evaporation where unsaturated, then autoconversion and accretion; returns each process's
    rate by level (kg/kg/s), the evaporation and the water handed to the precipitation."""
    n, nlev = s['dq'].shape
    zflux = np.zeros(n)
    out = dict(evaporation=np.zeros((n, nlev)), to_precipitation=np.zeros((n, nlev)))
    for k in range(nlev):
        zq = s['dq'][:, k] + s['dqdt'][:, k] * c['deltsec2']
        zt = s['dt'][:, k] + s['dtdt'][:, k] * c['deltsec2']
        cold = (s['dqcw'][:, k] > 0.0) & (zt < c['TMELT'])
        zconv = np.where(cold, s['dqcw'][:, k] / c['deltsec'], 0.0)
        s['dqcw'][:, k] = np.where(cold, 0.0, s['dqcw'][:, k])
        zlcpe = c['ALV'] / (c['acpd'] * (1.0 + c['ADV'] * s['dq'][:, k]))
        zqsat, zcor = saturation(zt, s['dp'], s['sigma'][k], c)
        dry = (s['dqcw'][:, k] > 0.0) & (zq < zqsat)
        need = (zqsat - zq) / c['deltsec2'] / (1.0 + zlcpe * c['ra2'] * (c['TMELT'] - c['ra4']) * zqsat * zcor
                                                 / (zt - c['ra4']) ** 2)
        zevap = np.where(dry, np.minimum(need, s['dqcw'][:, k] / c['deltsec']), 0.0)
        s['dqcw'][:, k] = np.where(dry, np.maximum(s['dqcw'][:, k] - zevap * c['deltsec'], 0.0), s['dqcw'][:, k])
        s['dqdt'][:, k] = np.where(dry, s['dqdt'][:, k] + zevap, s['dqdt'][:, k])
        s['dtdt'][:, k] = np.where(dry, s['dtdt'][:, k] - zevap * zlcpe, s['dtdt'][:, k])
        zrho = s['dp'] * s['sigma'][k] / (c['gascon'] * zt)
        qc = s['dqcw'][:, k]
        zauto = np.where(qc >= QC_AUTO_MIN, AUTO[0] * np.maximum(qc, QC_AUTO_MIN) ** AUTO[1] * prm['ncw_cm3'] ** AUTO[2],
                         0.0)
        wet = (qc >= QC_ACC_MIN) & (zflux > 0.0)
        safe = np.where(wet, zflux, prm['fr0_cw'])
        zqr = safe / (zrho * prm['vtr0_cw'] * (safe / prm['fr0_cw']) ** prm['vtrexp_cw'])
        zacc = np.where(wet, ACCRETION[0] * (np.maximum(qc, QC_ACC_MIN) * zqr) ** ACCRETION[1], 0.0)
        zx = np.minimum(zauto + zacc, qc / c['deltsec'])
        zconv = zconv + zx
        s['dqcw'][:, k] = np.maximum(qc - zx * c['deltsec'], 0.0)
        s['dprll'][:, k] = s['dprll'][:, k] + zconv * s['dsigma'][k]
        zflux = zflux + (s['dprll'][:, k] - s['dprsll'][:, k] + s['dprcl'][:, k] - s['dprscl'][:, k]) * s['dp'] / c['ga']
        out['evaporation'][:, k] = zevap
        out['to_precipitation'][:, k] = zconv
    return out


# ------------------------------------------------------------------------------------------------ the Fortran

MODULE_EDITS = [
    ("      real :: gamma   = 0.01  ! tuning parameter for evaporation of precip.\n",
     "      real :: gamma   = 0.01  ! tuning parameter for evaporation of precip.\n"
     f"      real :: ncw_cm3   = 100.   ! {MARKER} (pumamod's ncwat 1): the droplets' number (per cm3)\n"
     "      real :: vtr0_cw   = 0.99   ! the rain's mass-weighted fall speed (m/s) at the flux fr0_cw\n"
     "      real :: fr0_cw    = 3.4722222222222224e-05 ! that flux (kg/m2/s; 3 mm/day)\n"
     "      real :: vtrexp_cw = 0.20   ! the fall speed's exponent on the flux\n"),
    ("      real :: dprsll(NHOR,NLEV)     ! large scale snow for each level\n",
     "      real :: dprsll(NHOR,NLEV)     ! large scale snow for each level\n"
     f"      real :: dcwcond(NHOR,NLEV) = 0. ! {MARKER}: what mklsp hands the cloud (kg/kg/s)\n"),
    ("      namelist/rainmod_nl/kbeta,nprl,nprc,ndca,ncsurf,nmoment,nshallow  &\n",
     f"      namelist/rainmod_nl/ncw_cm3,vtr0_cw,fr0_cw,vtrexp_cw   ! {MARKER}\n"
     "      namelist/rainmod_nl/kbeta,nprl,nprc,ndca,ncsurf,nmoment,nshallow  &\n"),
    ("      call mpbcr(gamma)\n",
     "      call mpbcr(gamma)\n      call mpbcr(ncw_cm3)\n      call mpbcr(vtr0_cw)\n      call mpbcr(fr0_cw)\n"
     "      call mpbcr(vtrexp_cw)\n"),
    ("      if(nprl == 1) call mklsp\n",
     "      if(nprl == 1) call mklsp\n"
     f"      if(nprl == 1 .and. ncwat == 1) call mkcwat   ! {MARKER}\n"),
]
MKLSP_EDITS = [
    ("      dprll(:,:)=0.\n      dprsll(:,:)=0.\n",
     "      dprll(:,:)=0.\n      dprsll(:,:)=0.\n"
     f"      if(ncwat == 1) dcwcond(:,:)=0.   ! {MARKER}\n"),
    ("        dprll(:,jlev)=-AMIN1(zdqdt(:),0.0)*dsigma(jlev)\n\n       endwhere ! (zq(:) > zqsto(:))\n",
     "        dprll(:,jlev)=-AMIN1(zdqdt(:),0.0)*dsigma(jlev)\n\n       endwhere ! (zq(:) > zqsto(:))\n"
     f"!      {MARKER} (ncwat 1): what condenses at or above freezing, by the layer's temperature after\n"
     "!      condensing, stays as the layer's cloud water\n"
     "       if(ncwat == 1) then\n"
     "        where(dprll(:,jlev) > 0. .and. dt(:,jlev)+dtdt(:,jlev)*deltsec2 >= TMELT)\n"
     "         dcwcond(:,jlev)=dprll(:,jlev)/dsigma(jlev)\n"
     "         dqcw(:,jlev)=dqcw(:,jlev)+dcwcond(:,jlev)*deltsec\n"
     "         dprll(:,jlev)=0.\n"
     "        endwhere\n"
     "       endif\n"),
]
MKCWAT_CODE = f"""
!     ==================
!     SUBROUTINE MKCWAT
!     ==================
!
!     {MARKER} (ncwat 1), after mklsp, down each column: the cloud water of a layer below freezing goes to the
!     precipitation at once; where the layer's air is unsaturated the cloud water evaporates toward saturation, with
!     mklsp's latent-heat correction, at most the water there; the rest turns to rain by Khairoutdinov and Kogan's
!     (2000) autoconversion, 1350 qc**2.47 Nc**(-1.79), and accretion, 67 (qc qr)**1.15 (kg/kg/s; qc, qr in kg/kg,
!     Nc per cm3), as CM1's Morrison scheme takes them, qr the rain formed above over the air's density and the
!     fall speed vtr0_cw (F/fr0_cw)**vtrexp_cw; what the cloud hands on joins dprll for mkrain
!     (climate/gcm/condensate.py).
!
      subroutine mkcwat
      use rainmod
      real zq(NHOR),zt(NHOR),zqsat(NHOR),zcor(NHOR),zlcpe(NHOR),zrho(NHOR)
      real zevap(NHOR),zconv(NHOR),zauto(NHOR),zacc(NHOR),zqr(NHOR)
      real zflux(NHOR),zx(NHOR),zsafe(NHOR),zqc(NHOR)
      zflux(:)=0.
      do jlev=1,NLEV
       zq(:)=dq(:,jlev)+dqdt(:,jlev)*deltsec2
       zt(:)=dt(:,jlev)+dtdt(:,jlev)*deltsec2
       zconv(:)=0.
       zevap(:)=0.
       where(dqcw(:,jlev) > 0. .and. zt(:) < TMELT)
        zconv(:)=dqcw(:,jlev)/deltsec
        dqcw(:,jlev)=0.
       endwhere
       zlcpe(:)=ALV/(acpd*(1.+ADV*dq(:,jlev)))
       zqsat(:)=rdbrv*ra1*exp(ra2*(zt(:)-TMELT)/(zt(:)-ra4))/(dp(:)*sigma(jlev))
       zqsat(:)=AMIN1(zqsat(:),rdbrv)
       zcor(:)=1./(1.-(1./rdbrv-1.)*zqsat(:))
       zqsat(:)=zqsat(:)*zcor(:)
       where(dqcw(:,jlev) > 0. .and. zq(:) < zqsat(:))
        zevap(:)=AMIN1((zqsat(:)-zq(:))/deltsec2                           &
     &       /(1.0+zlcpe(:)*ra2*(TMELT-ra4)*zqsat(:)*zcor(:)/(zt(:)-ra4)**2)  &
     &       ,dqcw(:,jlev)/deltsec)
        dqcw(:,jlev)=AMAX1(dqcw(:,jlev)-zevap(:)*deltsec,0.)
        dqdt(:,jlev)=dqdt(:,jlev)+zevap(:)
        dtdt(:,jlev)=dtdt(:,jlev)-zevap(:)*zlcpe(:)
       endwhere
       zrho(:)=dp(:)*sigma(jlev)/(gascon*zt(:))
       zqc(:)=dqcw(:,jlev)
       zauto(:)=0.
       zacc(:)=0.
       where(zqc(:) >= 1.E-6)
        zauto(:)=1350.*zqc(:)**2.47*ncw_cm3**(-1.79)
       endwhere
       zsafe(:)=fr0_cw
       where(zqc(:) >= 1.E-8 .and. zflux(:) > 0.) zsafe(:)=zflux(:)
       zqr(:)=zsafe(:)/(zrho(:)*vtr0_cw*(zsafe(:)/fr0_cw)**vtrexp_cw)
       where(zqc(:) >= 1.E-8 .and. zflux(:) > 0.)
        zacc(:)=67.*(zqc(:)*zqr(:))**1.15
       endwhere
       zx(:)=AMIN1(zauto(:)+zacc(:),zqc(:)/deltsec)
       zconv(:)=zconv(:)+zx(:)
       dqcw(:,jlev)=AMAX1(zqc(:)-zx(:)*deltsec,0.)
       dprll(:,jlev)=dprll(:,jlev)+zconv(:)*dsigma(jlev)
       zflux(:)=zflux(:)+(dprll(:,jlev)-dprsll(:,jlev)+dprcl(:,jlev)       &
     &                   -dprscl(:,jlev))*dp(:)/ga
       if(ndiaggp==1 .and. ndiaggp3d >= 25) then
        dgp3d(:,jlev,22)=dqcw(:,jlev)
        dgp3d(:,jlev,23)=dcwcond(:,jlev)*dp(:)*dsigma(jlev)/ga/1000.
        dgp3d(:,jlev,24)=zevap(:)*dp(:)*dsigma(jlev)/ga/1000.
        dgp3d(:,jlev,25)=zconv(:)*dp(:)*dsigma(jlev)/ga/1000.
       endif
      enddo
      return
      end subroutine mkcwat
"""
ROUTINE_EDITS = [("      end subroutine mklsp\n", "      end subroutine mklsp\n" + MKCWAT_CODE)]
CLOUD_WATER_EDITS = MODULE_EDITS + MKLSP_EDITS + ROUTINE_EDITS
CLOUD_WATER_PATCHES = [('rainmod.f90', MARKER, CLOUD_WATER_EDITS)]
# The cloud water's diagnostics, with ndiaggp 1 and four slots beyond PlaSim's 21 (NDIAGGP3D 4 in plasim_namelist):
# the water each layer's cloud holds (a snapshot at each output's last step), the water mklsp hands it, the vapour its
# evaporation gives back and the water it hands to the precipitation, the last three water fluxes as codes 31-37 are.
DIAGNOSTICS = {42: ('qcw', 'cloud_water_held_by_the_layer', 'kg kg-1'),
               43: ('dcw_cond', 'water_the_cloud_gains_from_large_scale_condensation', 'm s-1'),
               44: ('dq_cwevap', 'water_vapour_gained_in_the_layer_by_the_evaporation_of_cloud_water', 'm s-1'),
               45: ('dcw_rain', 'water_the_cloud_hands_to_the_precipitation', 'm s-1')}
EXTRA_SLOTS = 4


def namelist_values(prm: dict = PARAMETERS) -> dict:
    """rainmod_namelist's entries for the cloud water's constants."""
    return {k.upper(): repr(float(v)) for k, v in prm.items()}


def setting(prm: dict = PARAMETERS) -> dict:
    """The cloud water as a run's configuration records it."""
    fit = fall_speed_fit()
    return dict(scheme='cloud', parameters={k: float(v) for k, v in prm.items()},
                autoconversion='Khairoutdinov and Kogan (2000): 1350 qc^2.47 Nc^-1.79 kg/kg/s, Nc per cm3, as CM1\'s '
                               'Morrison scheme (PRC) takes it, for qc of at least 1e-6 kg/kg',
                accretion='67 (qc qr)^1.15 kg/kg/s, as CM1\'s Morrison scheme (PRA) takes it, for qc of at least 1e-8 '
                          'kg/kg, qr the rain formed above over the air\'s density and the fall speed',
                droplets='the CM1 boxes\' 100 per cm3 (their case\'s droplets_cm3)',
                fall_speed=dict(fit, note='Marshall and Palmer\'s intercept, Earth\'s, on the lunar fall law'),
                held='at or above freezing (the layer\'s temperature after condensing); below it the condensate '
                     'falls at once, as in PlaSim',
                transport='none: the water stays in its layer and grid point',
                radiation="PlaSim's diagnosed clouds, unchanged",
                patch='climate/gcm/condensate.py CLOUD_WATER_PATCHES (rainmod ncw_cm3, vtr0_cw, fr0_cw, vtrexp_cw; '
                      'pumamod ncwat and dqcw through exoplasim_run.PLASIM_SWITCH_PATCHES)')


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.parse_args(argv)
    fit = fall_speed_fit()
    print(json.dumps(fit, indent=1))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
