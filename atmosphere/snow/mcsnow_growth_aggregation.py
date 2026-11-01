"""Lunar snowflakes from crystals that grow and aggregate together, in the lunar column and along storm paths.

The 24 October aggregation runs (mcsnow_aggregation.py) injected crystals at the
size they reach at the freezing level.  Here McSnow grows them from 4 um spheres
with the frozen improved component while they settle and aggregate, so growth
and aggregation act together.  Two settings carry a storm's water supply:

* the reference lunar column (the 301 K ladder from -15 C, 16.28 km) at water
  saturation, either sharing a steady supply among the crystals present in each
  cell in proportion to the growth each requests, or without a supply limit;
* the 354 archived storm paths (mcsnow_weather.py, earth_error_band.py), each
  carrying its own column whose air follows the storm at the path's position, in
  the archived vapour or held at water saturation, with a puff of crystals
  released once in the path's release cell.

Patches to the pinned McSnow source, each declared here and off by default:
* the three of mcsnow_aggregation.py (monomer laws in the namelist, a fall-speed
  multiplier, the lunar atmosphere constants);
* the improved component's one change, the extra oblate mean-ventilation term
  omitted (plate_base), with the error band's capacitance and axis-growth
  multipliers (band_mass, band_axes);
* deposition substeps limited to a fractional mass change (depo_max_fraction),
  removal below a sublimation residual (depo_residual) and a single-crystal size
  censor (depo_censor), as the storm probes ended;
* a cell water budget: a supply Q spread evenly over the column depth and shared
  by demand within each cell (water_supply).  For a column of identical crystals
  it equals the per-crystal rule of the single-crystal column program;
* cells driven by their own forcing tables (box_array), grouped into sub-columns
  within which particles settle (box_sub, box_settle), seeded once (ncl_until),
  with arrivals at 0 C recorded and the unused bulk scheme skipped (skip_bulk);
* binary superparticle snapshots (dt_snap).
With every new option at its default the build reproduces the 24 October runs;
parity() checks that byte for byte.
"""
from __future__ import annotations

import argparse
import difflib
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np

from shared.constants import MOON_SURFACE_GRAVITY
from atmosphere.snow import mcsnow_aggregation as agg
from atmosphere.snow import welss_aggregation_build as builder
from atmosphere.snow import welss_aggregation as author

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CACHE = ROOT / 'research/runs/snow_agent_2/growagg'
SOURCE = agg.SOURCE
PRODUCT = HERE / 'results/mcsnow_growth_aggregation.json'
ICE_DENSITY = 919.          # McSnow's rhoi, as the column program uses it
SEED_DIAMETER = 4e-6
SEED_MASS = ICE_DENSITY * math.pi / 6 * SEED_DIAMETER ** 3
RESIDUAL_MASS = ICE_DENSITY * math.pi / 6 * 2e-7 ** 3   # the storm probes' end: a dense 0.2 um sphere
CENSOR_M = 1.                                            # the storm probes' size censor for single crystals
REFERENCE_FLUX = agg.LUNAR_REFERENCE_FLUX
STORM_RUNS = ROOT / 'research/runs/snow_agent/earth_error_band/weather_runs'
EARTH_BAND = HERE / 'results/earth_error_band.json'
SIZE_RANGE = HERE / 'results/lunar_size_range.json'
_rep = agg._replace_once


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# --- source patches ----------------------------------------------------------------

NML_VARIABLES = ('water_supply, depo_max_fraction, band_mass, band_axes, ncl_until, dt_snap, plate_base, box_array, &\n'
                 '            box_sub, box_settle, skip_bulk, depo_residual, depo_censor')


def _namelist_source(text):
    text = _rep(text, '    coll_kern_type    ! 0: A-kernel, 1: D-kernel\n',
                '    coll_kern_type, & ! 0: A-kernel, 1: D-kernel\n'
                '    box_sub           ! cells in each forced sub-column (Terluna path columns)\n')
    text = _rep(text, '    mono_alf, mono_bet, mono_sig, mono_gam  ! monomer m-D and A-D laws (Terluna)\n  LOGICAL     ::   &\n',
                '    mono_alf, mono_bet, mono_sig, mono_gam, &  ! monomer m-D and A-D laws (Terluna)\n'
                '    water_supply,     &  ! [kg/m2/s] supply spread over the column, shared by demand per cell; <0 off (Terluna)\n'
                '    depo_max_fraction, & ! [1] largest fractional mass change per deposition substep; 0 off (Terluna)\n'
                '    band_mass, band_axes, & ! [1] capacitance and axis-growth multipliers (Terluna error band)\n'
                '    ncl_until,        &  ! [s] nucleation only up to this time; <0 always (Terluna)\n'
                '    dt_snap,          &  ! [s] interval of binary superparticle snapshots; 0 off (Terluna)\n'
                '    depo_residual,    &  ! [kg] a sublimating particle below this mass is removed; 0 off (Terluna)\n'
                '    depo_censor          ! [m] a single crystal growing to this span is removed and recorded; 0 off (Terluna)\n'
                '  LOGICAL     ::   &\n'
                '    plate_base,    &  ! omit the extra oblate mean-ventilation term (Terluna improved component)\n'
                '    box_array,     &  ! independent forced cells without settling (Terluna parcels)\n'
                '    box_settle,    &  ! forced cells settle within sub-columns of box_sub cells (Terluna)\n'
                '    skip_bulk,     &  ! skip the two-moment bulk scheme run alongside, which no particle uses (Terluna)\n')
    text = _rep(text, 'coll_kern_type, sp_kern_sig, &\n            prolate_intp, band_speed, mono_alf, mono_bet, mono_sig, mono_gam\n',
                'coll_kern_type, sp_kern_sig, &\n            prolate_intp, band_speed, mono_alf, mono_bet, mono_sig, mono_gam, &\n'
                f'            {NML_VARIABLES}\n')
    text = _rep(text, '  sp_frag, coll_kern_type, sp_kern_sig, prolate_intp, band_speed, mono_alf, mono_bet, mono_sig, mono_gam\n',
                '  sp_frag, coll_kern_type, sp_kern_sig, prolate_intp, band_speed, mono_alf, mono_bet, mono_sig, mono_gam, &\n'
                f'  {NML_VARIABLES}\n')
    text = _rep(text, '  mono_gam       = 2.00_wp\n',
                '  mono_gam       = 2.00_wp\n  water_supply   = -1.0_wp\n  depo_max_fraction = 0.0_wp\n'
                '  band_mass      = 1.0_wp\n  band_axes      = 1.0_wp\n  ncl_until      = -1.0_wp\n  dt_snap        = 0.0_wp\n'
                '  plate_base     = .FALSE.\n  box_array      = .FALSE.\n'
                '  box_sub        = 1\n  box_settle     = .FALSE.\n  skip_bulk      = .FALSE.\n  depo_residual  = 0.0_wp\n'
                '  depo_censor    = 0.0_wp\n')
    return text


DEPOSITION_ROUTINES = '''!
! Terluna: deposition in substeps whose fractional mass change stays below depo_max_fraction.
! With depo_max_fraction <= 0 it is one call of sp_depo, as before.
!
SUBROUTINE sp_depo_limited( sp, atmo, dt, q_heat, l_delnode, dep_rate, dep_mono, heat_rate )
  TYPE(t_sp),   INTENT(inout)  :: sp
  TYPE(t_atmo), INTENT(in)     :: atmo
  REAL(wp),     INTENT(in)     :: dt
  REAL(wp),     INTENT(inout)  :: q_heat
  LOGICAL,      INTENT(out)    :: l_delnode
  REAL(wp),     INTENT(inout)  :: dep_rate, dep_mono, heat_rate

  TYPE(t_sp) :: trial
  REAL(wp)   :: remaining, h, qt, dr, dmo, hr, m0, m1
  LOGICAL    :: gone
  INTEGER    :: halvings

  depo_censored = .FALSE.
  IF ( depo_max_fraction <= 0._wp ) THEN
    CALL sp_depo( sp, atmo, dt, q_heat, l_delnode, dep_rate, dep_mono, heat_rate )
    RETURN
  END IF
  l_delnode = .FALSE.
  remaining = dt
  h = dt
  DO WHILE ( remaining > 1.e-9_wp*dt )
    h = MIN( h, remaining )
    halvings = 0
    DO
      trial = sp
      qt = q_heat
      dr = 0._wp
      dmo = 0._wp
      hr = 0._wp
      CALL sp_depo( trial, atmo, h, qt, gone, dr, dmo, hr )
      IF ( gone ) EXIT
      m0 = sp%m_i + sp%m_f + sp%m_r + sp%m_w
      m1 = trial%m_i + trial%m_f + trial%m_r + trial%m_w
      IF ( .NOT. (m1 == m1) .OR. ABS(m1) > HUGE(m1) ) THEN
        WRITE(*,*) 'Terluna: non-finite deposition step; particle m_i, V_i, phi, d, mm, xi, z, m_w:', &
          & sp%m_i, sp%V_i, sp%phi, sp%d, sp%mm, sp%xi, sp%z, sp%m_w
        WRITE(*,*) 'Terluna: air T, p, pv, ssat, rho; step:', atmo%T, atmo%p, atmo%pv, atmo%ssat, atmo%rho, h
        ERROR STOP 'Terluna: non-finite deposition step'
      END IF
      IF ( m1 < depo_residual ) THEN   ! sublimated to the residual: removed, as the storm probes ended
        gone = .TRUE.
        EXIT
      END IF
      IF ( ABS(m1-m0) <= depo_max_fraction*m0 .OR. halvings >= 40 ) EXIT
      h = 0.5_wp*h
      halvings = halvings + 1
    END DO
    sp = trial
    q_heat = qt
    dep_rate  = dep_rate  + dr*h/dt
    dep_mono  = dep_mono  + dmo*h/dt
    heat_rate = heat_rate + hr
    remaining = remaining - h
    IF ( gone ) THEN
      l_delnode = .TRUE.
      RETURN
    END IF
    IF ( depo_censor > 0._wp .AND. sp%mm < 1.01_wp ) THEN   ! the storm probes' size censor, for single crystals
      CALL m2d_particle( sp )
      IF ( sp%d >= depo_censor ) THEN
        depo_censored = .TRUE.
        l_delnode = .TRUE.
        RETURN
      END IF
    END IF
    IF ( halvings == 0 ) h = 2._wp*h
  END DO
END SUBROUTINE sp_depo_limited

!
! Terluna: the deposition rate a particle requests in its present state, before any supply limit.
!
REAL(wp) FUNCTION depo_demand( sp, atmo, dt ) RESULT( rate )
  TYPE(t_sp),   INTENT(in) :: sp
  TYPE(t_atmo), INTENT(in) :: atmo
  REAL(wp),     INTENT(in) :: dt

  TYPE(t_sp) :: trial
  REAL(wp)   :: vt, mdiff, qdiff, lhc, V_new, phi_tend, saved

  rate = 0._wp
  IF ( sp%m_w >= TINY(sp%m_w) ) RETURN
  trial = sp
  CALL m2d_particle( trial )
  vt = vterm_snow( atmo, trial )
  saved = supply_factor
  supply_factor = 1._wp
  CALL vapourDiffusion( trial, vt, dt, atmo, mdiff, qdiff, lhc, V_new, phi_tend )
  supply_factor = saved
  rate = MAX( mdiff, 0._wp )
END FUNCTION depo_demand

'''


def _deposition_source(text):
    text = _rep(text, '  USE mo_sp_nml,       ONLY: sp_habit\n',
                '  USE mo_sp_nml,       ONLY: sp_habit, plate_base, band_mass, band_axes, depo_max_fraction, depo_residual, &\n'
                '                              depo_censor\n')
    text = _rep(text, '  INTEGER, PARAMETER :: subcycles = 1\n',
                '  INTEGER, PARAMETER :: subcycles = 1\n  REAL(wp), PUBLIC :: supply_factor = 1._wp, requested_rate = 0._wp\n'
                '  LOGICAL,  PUBLIC :: depo_censored = .FALSE.\n')
    text = _rep(text, '  PUBLIC  :: sp_depo, sp_depo_iteration, ventilation_coef_IGF\n',
                '  PUBLIC  :: sp_depo, sp_depo_iteration, ventilation_coef_IGF, sp_depo_limited, depo_demand\n')
    start = text.index('SUBROUTINE vapourDiffusion( sp,')
    end = text.index('END SUBROUTINE vapourDiffusion', start)
    part = text[start:end]
    part = _rep(part, '  C = capacitance( sp, Re )\n', '  C = capacitance( sp, Re ) * band_mass\n')
    anchor = '  dMdt =  4.0_wp * pi * C *   Dv * (atmo%pv - es)                         / (Rv * atmo%T * (1._wp+Hi))\n'
    part = _rep(part, anchor, anchor + '  requested_rate = dMdt\n'
                '  IF ( supply_factor < 1._wp .AND. dMdt > 0._wp ) dMdt = dMdt*supply_factor\n')
    text = text[:start] + part + text[end:]
    text = _rep(text, '  IGF_vent = IGF * ventilation_coef_IGF( ReSc, r_c, V_tot**z1_3, sp%phi ) /                       &\n'
                      '                   ventilation_coef_IGF( ReSc, r_a, V_tot**z1_3, sp%phi ) \n',
                '  IGF_vent = IGF * ventilation_coef_IGF( ReSc, r_c, V_tot**z1_3, sp%phi ) /                       &\n'
                '                   ventilation_coef_IGF( ReSc, r_a, V_tot**z1_3, sp%phi ) \n'
                '  IGF_vent = IGF_vent * band_axes\n')
    text = _rep(text, '    fv      = fv + c2 / phi * X_equiv**(3./2.)\n',
                '    IF ( .NOT. plate_base ) fv      = fv + c2 / phi * X_equiv**(3./2.)\n')
    text = _rep(text, '!\n! implicit euler for particle temperature', DEPOSITION_ROUTINES + '!\n! implicit euler for particle temperature')
    return text


BOX_ATMOSPHERE = '''
!
! Terluna: forcing tables for independent cells (box_array). Each cell k follows its own
! temperature, pressure and vapour pressure in time; a negative vapour pressure means
! water saturation at the cell's temperature.  box_end(k) ends the cell's evolution.
!
SUBROUTINE load_box_forcing()
  INTEGER :: k, it, ncell, unit, nsub
  LOGICAL :: binary
  REAL(wp) :: header(4)
  REAL(4), ALLOCATABLE :: tmp(:,:)
  INQUIRE(FILE='box_forcing.bin', EXIST=binary)
  IF ( binary ) THEN
    ! float64 ncell, ntime, dt, nsub; float64 box_end(ncell), box_release(nsub); float32 T, p, pv (ntime, ncell)
    OPEN(NEWUNIT=unit, FILE='box_forcing.bin', ACCESS='STREAM', FORM='UNFORMATTED', STATUS='old', ACTION='read')
    READ(unit) header
    ncell = NINT(header(1)); box_nt = NINT(header(2)); box_dt = header(3); nsub = NINT(header(4))
    IF ( ncell /= grid%nz ) STOP 'box forcing: cell count differs from nz'
    ALLOCATE( box_T(box_nt,ncell), box_p(box_nt,ncell), box_pv(box_nt,ncell), box_end(ncell), box_release(nsub) )
    ALLOCATE( tmp(box_nt,ncell) )
    READ(unit) box_end
    READ(unit) box_release
    READ(unit) tmp
    box_T = REAL(tmp,wp)
    READ(unit) tmp
    box_p = REAL(tmp,wp)
    READ(unit) tmp
    box_pv = REAL(tmp,wp)
    CLOSE(unit)
    DEALLOCATE( tmp )
    RETURN
  END IF
  ALLOCATE( box_release(1) )
  box_release = 1._wp
  OPEN(NEWUNIT=unit, FILE='box_forcing.txt', STATUS='old', ACTION='read')
  READ(unit,*) ncell, box_nt, box_dt
  IF ( ncell /= grid%nz ) STOP 'box forcing: cell count differs from nz'
  ALLOCATE( box_T(box_nt,ncell), box_p(box_nt,ncell), box_pv(box_nt,ncell), box_end(ncell) )
  DO k = 1, ncell
    READ(unit,*) box_end(k)
    DO it = 1, box_nt
      READ(unit,*) box_T(it,k), box_p(it,k), box_pv(it,k)
    END DO
  END DO
  CLOSE(unit)
END SUBROUTINE load_box_forcing

SUBROUTINE update_box_forcing( time )
  REAL(wp), INTENT(in) :: time
  INTEGER  :: k, i0, i1
  REAL(wp) :: s, w, T, p, pv
  s  = time / box_dt
  i0 = MIN( INT(s) + 1, box_nt )
  i1 = MIN( i0 + 1, box_nt )
  w  = MERGE( s - REAL(i0-1,wp), 0._wp, i0 < box_nt )
  DO k = 1, grid%nz
    T  = (1._wp-w)*box_T(i0,k)  + w*box_T(i1,k)
    p  = (1._wp-w)*box_p(i0,k)  + w*box_p(i1,k)
    pv = (1._wp-w)*box_pv(i0,k) + w*box_pv(i1,k)
    patmo(k)%T     = T
    patmo(k)%p     = p
    patmo(k)%psatw = psat_water( T )
    patmo(k)%psati = psat_ice( T )
    IF ( box_pv(i0,k) < 0._wp ) pv = patmo(k)%psatw
    patmo(k)%pv    = pv
    patmo(k)%rh    = 100._wp * pv / patmo(k)%psatw
    patmo(k)%ssat  = pv / patmo(k)%psati - 1._wp
    patmo(k)%eta   = dynamic_viscosity( T )
    patmo(k)%etai  = 1._wp / patmo(k)%eta
    patmo(k)%rho   = (p - pv)/(Rd*T) + pv/(Rv*T)
    patmo(k)%qv    = (Rd/Rv)*pv/(p - (1._wp-Rd/Rv)*pv)
    patmo(k)%qc    = 0._wp
  END DO
END SUBROUTINE update_box_forcing
'''


def _atmosphere_source(text):
    text = _rep(text, '  TYPE(t_diag), ALLOCATABLE :: diag(:)         ! 1D profile with diagnostic fields\n',
                '  TYPE(t_diag), ALLOCATABLE :: diag(:)         ! 1D profile with diagnostic fields\n'
                '  REAL(wp), ALLOCATABLE, PUBLIC :: box_T(:,:), box_p(:,:), box_pv(:,:), box_end(:)  ! Terluna parcels\n'
                '  REAL(wp), ALLOCATABLE, PUBLIC :: box_release(:)  ! release cell of each sub-column (Terluna)\n'
                '  INTEGER,  PUBLIC :: box_nt = 0\n  REAL(wp), PUBLIC :: box_dt = 600._wp\n')
    text = _rep(text, '  PUBLIC  ::            &\n    diag,               &\n',
                '  PUBLIC  ::            &\n    load_box_forcing,   &\n    update_box_forcing, &\n    diag,               &\n')
    text = _rep(text, '\nEND MODULE mo_atmo', BOX_ATMOSPHERE + '\nEND MODULE mo_atmo')
    text = _rep(text, '  USE mo_sp_nml,       ONLY: lwc0, interpol_atmos\n',
                '  USE mo_sp_nml,       ONLY: lwc0, interpol_atmos, box_array\n')
    text = _rep(text, '  atmo%T = T_3 - tlr*z\n', '  atmo%T = T_3 - tlr*z\n  IF ( box_array ) atmo%T = T_3   ! forced cells: replaced by their tables\n')
    return text


def _nucleation_source(text):
    text = _rep(text, '  USE mo_sp_nml,       ONLY: agg_hist, xi0, nclmass, nclrate, nu_gam\n',
                '  USE mo_sp_nml,       ONLY: agg_hist, xi0, nclmass, nclrate, nu_gam, box_array, box_sub\n')
    text = _rep(text, '  USE mo_grid,         ONLY: t_grid, coord2index\n', '  USE mo_grid,         ONLY: t_grid, coord2index, grid\n')
    text = _rep(text, '  USE mo_atmo,         ONLY: t_atmo, get_atmo\n', '  USE mo_atmo,         ONLY: t_atmo, get_atmo, box_release\n')
    text = _rep(text, '  INTEGER :: i, ncl_in_dbox\n', '  INTEGER :: i, ncl_in_dbox, per, s, j, cell\n')
    text = _rep(text, '    sp%z = nh1 + (nh2-nh1)*rand\n',
                '    sp%z = nh1 + (nh2-nh1)*rand\n'
                '    IF ( box_array ) sp%z = nh1 + (nh2-nh1)*(REAL(i,wp)-0.5_wp)/REAL(ncl_in_dbox,wp)\n'
                '    IF ( box_array .AND. box_sub > 1 ) THEN   ! a puff in each sub-column\'s release cell\n'
                '      per  = MAX(1, ncl_in_dbox/SIZE(box_release))\n'
                '      s    = MIN((i-1)/per, SIZE(box_release)-1)\n'
                '      j    = MOD(i-1, per)\n'
                '      cell = s*box_sub + NINT(box_release(s+1))\n'
                '      sp%z = grid%dom_top - (REAL(cell-1,wp) + (REAL(j,wp)+0.5_wp)/REAL(per,wp))*grid%dz\n'
                '    END IF\n')
    return text


ARRIVAL_ROUTINE = '''    !
    ! Terluna: append a superparticle leaving a forced sub-column (reason 1: its cell reached 0 C;
    ! 2: it settled below the sub-column; 3: a single crystal reached the size censor) to arrivals.bin:
    ! time, sub-column, height, mass, maximum dimension, monomers, multiplicity, reason, as 8-byte reals.
    !
    SUBROUTINE record_arrival( sp, itc, ic, reason )
      TYPE(t_sp), INTENT(in) :: sp
      INTEGER,    INTENT(in) :: itc, ic, reason
      TYPE(t_sp) :: copy
      INTEGER    :: unit
      copy = sp
      CALL m2d_particle( copy )
      OPEN(NEWUNIT=unit, FILE='arrivals.bin', ACCESS='STREAM', FORM='UNFORMATTED', POSITION='APPEND', STATUS='UNKNOWN')
      WRITE(unit) REAL(itc,wp)*dtc, REAL((ic-1)/box_sub+1,wp), copy%z, copy%m_i+copy%m_f+copy%m_r+copy%m_w, copy%d, &
                & REAL(copy%mm,wp), REAL(copy%xi,wp), REAL(reason,wp)
      CLOSE(unit)
    END SUBROUTINE record_arrival

'''


SNAPSHOT_ROUTINE = '''    !
    ! Terluna: append every active superparticle (height, mass, maximum dimension, monomers,
    ! multiplicity, last fall speed) to snapshots.bin as a stream of 8-byte reals.
    !
    SUBROUTINE write_snapshot( cont, itc )
      TYPE(t_cont), INTENT(in) :: cont
      INTEGER,      INTENT(in) :: itc
      TYPE(t_sp) :: copy
      INTEGER    :: icell, j, n, unit
      REAL(wp), ALLOCATABLE :: buf(:,:)
      n = 0
      DO icell = 1, grid%nz
        DO j = 1, cont%cellSParray(icell,1,1)%noParts
          IF ( cont%sp_array(cont%cellSParray(icell,1,1)%startPos+j-1)%cellId >= 1 ) n = n + 1
        END DO
      END DO
      ALLOCATE( buf(6,MAX(n,1)) )
      n = 0
      DO icell = 1, grid%nz
        DO j = 1, cont%cellSParray(icell,1,1)%noParts
          copy = cont%sp_array(cont%cellSParray(icell,1,1)%startPos+j-1)
          IF ( copy%cellId < 1 ) CYCLE
          n = n + 1
          CALL m2d_particle( copy )
          buf(1,n) = copy%z
          buf(2,n) = copy%m_i + copy%m_f + copy%m_r + copy%m_w
          buf(3,n) = copy%d
          buf(4,n) = REAL(copy%mm,wp)
          buf(5,n) = REAL(copy%xi,wp)
          buf(6,n) = copy%v
        END DO
      END DO
      OPEN(NEWUNIT=unit, FILE='snapshots.bin', ACCESS='STREAM', FORM='UNFORMATTED', POSITION='APPEND', STATUS='UNKNOWN')
      WRITE(unit) REAL(itc,wp)*dtc, REAL(n,wp)
      IF ( n > 0 ) WRITE(unit) buf(:,1:n)
      CLOSE(unit)
      DEALLOCATE( buf )
    END SUBROUTINE write_snapshot

'''


def _loop_source(text):
    text = _rep(text, '                                & nucl_type, sp_kernel_id\n',
                '                                & nucl_type, sp_kernel_id, water_supply, box_array, ncl_until, dt_snap, &\n'
                '                                & box_sub, box_settle, skip_bulk\n')
    text = _rep(text, '  USE mo_atmo,              ONLY: t_atmo, get_atmo, diag\n',
                '  USE mo_atmo,              ONLY: t_atmo, get_atmo, diag, load_box_forcing, update_box_forcing, box_end\n')
    text = _rep(text, '  USE mo_depo,              ONLY: sp_depo, sp_depo_iteration\n',
                '  USE mo_depo,              ONLY: sp_depo, sp_depo_iteration, sp_depo_limited, depo_demand, supply_factor, &\n'
                '                                & depo_censored\n'
                '  USE mo_mass2diam,         ONLY: m2d_particle\n')
    text = _rep(text, '    REAL(wp)              :: m_shed_tot, dtg, q_heat\n',
                '    REAL(wp)              :: m_shed_tot, dtg, q_heat\n    REAL(wp)              :: demand, supply_cell\n')
    text = _rep(text, '      CALL two_moment_mcrph_init_mem(grid)\n',
                '      CALL two_moment_mcrph_init_mem(grid)\n      IF ( box_array ) CALL load_box_forcing()\n')
    text = _rep(text, '        diag%heat_rate = 0.0_wp\n',
                '        diag%heat_rate = 0.0_wp\n        IF ( box_array ) CALL update_box_forcing( REAL(itc-1,wp)*dtc )\n')
    text = _rep(text, '          nsp      = cont%cellSParray(jk,jc,jb)%noParts\n'
                      '          startPos = cont%cellSParray(jk,jc,jb)%startPos\n\n          IF( nsp > 0 ) THEN\n#ifndef NDEBUG\n',
                '          nsp      = cont%cellSParray(jk,jc,jb)%noParts\n'
                '          startPos = cont%cellSParray(jk,jc,jb)%startPos\n'
                '          IF ( box_array ) THEN\n'
                '            IF ( REAL(itc-1,wp)*dtc >= box_end(ic) ) CYCLE\n'
                '          END IF\n\n          IF( nsp > 0 ) THEN\n#ifndef NDEBUG\n')
    text = _rep(text, '          n_splinter_tot = 0_i8\n          m_shed_tot = 0.0_wp\n          multiplication_n_colls = 0\n',
                '#ifndef WDROPS\n'
                '          supply_factor = 1.0_wp\n'
                '          IF ( water_supply >= 0.0_wp .AND. diffusion_type == 1 ) THEN\n'
                '            demand = 0.0_wp\n'
                '            DO isp = 1, nsp\n'
                '              sp => cont%sp_array(startPos+isp-1)\n'
                '              IF( sp%cellId < 1 ) CYCLE\n'
                '              CALL get_atmo( sp%z, atmo )\n'
                '              demand = demand + REAL(sp%xi,wp)*depo_demand( sp, atmo, dtc/REAL(idtgm,wp) )\n'
                '            END DO\n'
                '            supply_cell = water_supply/(grid%dom_top-grid%dom_bottom)*grid%dbox\n'
                '            IF ( demand > supply_cell ) supply_factor = supply_cell/demand\n'
                '          END IF\n'
                '#endif\n'
                '          n_splinter_tot = 0_i8\n          m_shed_tot = 0.0_wp\n          multiplication_n_colls = 0\n')
    text = _rep(text, '                  CALL sp_depo( sp, atmo, dtg, q_heat, l_delnode, diag(ic)%dep_rate, diag(ic)%dep_mono, diag(ic)%heat_rate )\n',
                '                  CALL sp_depo_limited( sp, atmo, dtg, q_heat, l_delnode, diag(ic)%dep_rate, diag(ic)%dep_mono, &\n'
                '                                      & diag(ic)%heat_rate )\n'
                '                  IF ( l_delnode .AND. depo_censored ) CALL record_arrival( sp, itc, ic, 3 )\n')
    text = _rep(text, '            CALL sp_sedi( sp, atmo, dtc )\n',
                '            IF ( .NOT. box_array .OR. box_settle ) CALL sp_sedi( sp, atmo, dtc )\n'
                '            IF ( box_array .AND. box_settle ) THEN\n'
                '              IF ( sp%z < grid%dom_top - REAL((ic-1)/box_sub+1,wp)*REAL(box_sub,wp)*grid%dz ) THEN\n'
                '                CALL record_arrival( sp, itc, ic, 2 )\n'
                '                CALL delParticle( cont%cellSParray(jk,jc,jb), sp )\n'
                '                CYCLE partloop\n'
                '              END IF\n'
                '            END IF\n')
    text = _rep(text, '            CALL get_atmo( sp%z, atmo )\n            checkForHydroBreakup = (hydrobreakup_type > 0)\n',
                '            CALL get_atmo( sp%z, atmo )\n'
                '            IF ( box_array .AND. box_settle ) THEN\n'
                '              IF ( atmo%T >= 273.15_wp ) THEN\n'
                '                CALL record_arrival( sp, itc, ic, 1 )\n'
                '                CALL delParticle( cont%cellSParray(jk,jc,jb), sp )\n'
                '                CYCLE partloop\n'
                '              END IF\n'
                '            END IF\n'
                '            checkForHydroBreakup = (hydrobreakup_type > 0)\n')
    text = _rep(text, '        CALL two_moment_mcrph_init_values( grid, bnd_data%number_flux, bnd_data%mass_flux )\n',
                '        IF ( .NOT. skip_bulk ) THEN\n'
                '        CALL two_moment_mcrph_init_values( grid, bnd_data%number_flux, bnd_data%mass_flux )\n')
    text = _rep(text, '                        l_cv=.TRUE.          )\n', '                        l_cv=.TRUE.          )\n        END IF\n')
    text = _rep(text, '        IF (nucl_type == 1) THEN\n',
                '        IF (nucl_type == 1 .AND. (ncl_until < 0.0_wp .OR. REAL(itc,wp)*dtc <= ncl_until)) THEN\n')
    text = _rep(text, '        CALL write_solution( cont, itc, grid, twomom_data%d )\n',
                '        CALL write_solution( cont, itc, grid, twomom_data%d )\n'
                '        IF ( dt_snap > 0.0_wp ) THEN\n'
                '          IF ( MOD(itc, MAX(1, NINT(dt_snap/dtc))) == 0 ) CALL write_snapshot( cont, itc )\n'
                '        END IF\n')
    text = _rep(text, '  CONTAINS\n    \n    SUBROUTINE tic( t1 )\n',
                '  CONTAINS\n    \n' + ARRIVAL_ROUTINE + SNAPSHOT_ROUTINE + '    SUBROUTINE tic( t1 )\n')
    return text


def patched_sources(world):
    """The aggregation build's three patches plus this module's growth, budget, parcel and snapshot patches."""
    src = SOURCE / 'src'
    out = agg.patched_sources(world)
    out['mo_sp_nml.f90'] = _namelist_source(out['mo_sp_nml.f90'])
    out['mo_depo.f90'] = _deposition_source((src / 'mo_depo.f90').read_text())
    out['mo_atmo.f90'] = _atmosphere_source((src / 'mo_atmo.f90').read_text())
    out['mo_nucleation.f90'] = _nucleation_source((src / 'mo_nucleation.f90').read_text())
    out['mo_loop_array.f90'] = _loop_source((src / 'mo_loop_array.f90').read_text())
    return out


def build_path(world):
    """Each set of patched sources gets its own build, named by the hash of those sources."""
    sources = patched_sources(world)
    digest = hashlib.sha256(''.join(name + sources[name] for name in sorted(sources)).encode()).hexdigest()
    return CACHE / f'build_{world}_{digest[:12]}'


def build(world):
    path = build_path(world)
    if path.exists():
        return verify_build(world)
    manifest = json.loads(builder.MANIFEST.read_text())
    for item in manifest['files']:
        builder.verify((SOURCE / item['path']).read_bytes(), item)
    path.mkdir(parents=True)
    tree = path / 'source'
    shutil.copytree(SOURCE / 'src', tree / 'src')
    diff = ''
    for name, body in sorted(patched_sources(world).items()):
        original = (SOURCE / 'src' / name).read_text()
        (tree / 'src' / name).write_text(body)
        diff += ''.join(difflib.unified_diff(original.splitlines(True), body.splitlines(True),
                                             fromfile='pinned/' + name, tofile=f'{world}/' + name))
    main = (SOURCE / 'src/main.f90').read_text()
    (path / 'main.f90').write_text(builder.deterministic_main(main))
    (path / 'growth_aggregation.patch').write_text(diff)
    compiler = str(Path(shutil.which('gfortran')).resolve())
    files = sorted(p for p in (tree / 'src').rglob('*') if p.suffix.lower() == '.f90' and p.name != 'main.f90') + [path / 'main.f90']
    order = builder.compilation_order(files, compiler, tree / 'src')
    objects, commands = [], []
    with (path / 'build.log').open('w') as log:
        for i, p in enumerate(order):
            obj = path / f'{i:02d}_{p.stem}.o'
            objects.append(str(obj))
            cmd = [compiler, *builder.FLAGS, '-I' + str(tree / 'src'), '-I' + str(path), '-J' + str(path), '-c', str(p), '-o', str(obj)]
            commands.append(cmd)
            subprocess.run(cmd, cwd=path, stdout=log, stderr=subprocess.STDOUT, check=True)
        cmd = [compiler, *builder.FLAGS, *objects, '-o', str(path / 'MCSNOW')]
        commands.append(cmd)
        subprocess.run(cmd, cwd=path, stdout=log, stderr=subprocess.STDOUT, check=True)
    record = dict(schema='terluna.atmosphere.snow-mcsnow-growth-aggregation-build/1', world=world,
                  revision=manifest['revision'], source_dir=str(SOURCE), source_manifest_sha256=sha(builder.MANIFEST),
                  patch_sha256=sha(path / 'growth_aggregation.patch'),
                  compiler=subprocess.check_output([compiler, '--version'], text=True).splitlines()[0],
                  compiler_sha256=sha(compiler), flags=builder.FLAGS, commands=commands,
                  executable_sha256=sha(path / 'MCSNOW'))
    (path / 'build.json').write_text(json.dumps(record, indent=1) + '\n')
    return record


def verify_build(world):
    path = build_path(world)
    record = json.loads((path / 'build.json').read_text())
    if sha(path / 'MCSNOW') != record['executable_sha256']:
        raise ValueError('Executable changed')
    for name, body in patched_sources(world).items():
        if (path / 'source/src' / name).read_text() != body:
            raise ValueError('Patched source differs: ' + name)
    return record


# --- runs ----------------------------------------------------------------------------

TRIMMED = ('distribution', 'hei2massdens', 'hprofiles', 'twomom_d', 'atmo.dat')


def run(name, world, config, *, seed=20261025, forcing=None, forcing_link=None, keep_all=False):
    """Run McSnow once.  Profile outputs are written only at the start and end (dt_1dprof = time_end);
    keep_all=False then deletes those bulky profile files, keeping the final snapshot, the binary
    snapshots and the log, with every deleted file's hash recorded."""
    identity = verify_build(world)
    path = CACHE / 'runs' / name
    if (path / 'run.json').exists():
        record = json.loads((path / 'run.json').read_text())
        if record['configuration'] != config:
            raise ValueError('Existing run has different inputs: ' + name)
        return record
    path.mkdir(parents=True, exist_ok=True)
    (path / 'input.nml').write_text(author.namelist(config))
    for item in json.loads(builder.MANIFEST.read_text())['files']:
        if item['path'].startswith('input/'):
            shutil.copyfile(SOURCE / item['path'], path / Path(item['path']).name)
    if forcing is not None:
        (path / 'box_forcing.txt').write_text(forcing)
    if forcing_link is not None:
        (path / 'box_forcing.bin').symlink_to(Path(forcing_link).resolve())
    start = time.monotonic()
    with (path / 'run.log').open('w') as log:
        proc = subprocess.run(['bash', '-c', f'ulimit -v 8000000; exec {build_path(world) / "MCSNOW"} {seed}'],
                              cwd=path, stdout=log, stderr=subprocess.STDOUT)
    removed = {}
    if not keep_all:
        for p in sorted(path.iterdir()):
            if p.is_file() and p.name.startswith(TRIMMED):
                removed[p.name] = dict(sha256=sha(p), bytes=p.stat().st_size)
                p.unlink()
    record = dict(name=name, world=world, configuration=config, seed=seed, returncode=proc.returncode,
                  wall_s=time.monotonic() - start, executable_sha256=identity['executable_sha256'],
                  forcing_sha256=sha(path / 'box_forcing.txt') if forcing is not None else
                  (sha(path / 'box_forcing.bin') if forcing_link is not None else None),
                  removed_outputs=removed,
                  kept_outputs={p.name: sha(p) for p in sorted(path.iterdir())
                                if p.is_file() and not p.is_symlink() and p.suffix in ('.dat', '.bin')})
    (path / 'run.json').write_text(json.dumps(record, indent=1) + '\n')
    return record


def column_config(water_mm_h, flux_factor, *, saturated=False, stick=2, band=(1., 1., 1.), sp_per_step=2.,
                  area=5., duration_h=72., dz=20., dtc=10., source_fraction=.1, depo_fraction=.01, snap_h=12.,
                  collisions=True, nu_gam=10., top_C=-15.):
    """The reference lunar column at water saturation with 4 um seeds entering the top tenth.  Water-limited
    unless saturated: the supply water_mm_h is spread over the column depth and shared by demand in each cell.
    top_C below -15 deepens the column on the same lapse rate, up to the habit table's -30 C."""
    lunar = agg.lunar_column()
    depth = lunar['depth_15c'] if top_C == -15. else -top_C / lunar['tlr']
    if top_C < -30.:
        raise ValueError('The habit table ends at -30 C')
    number_flux = REFERENCE_FLUX * flux_factor
    xi = max(1, int(round(number_flux * area * dtc / sp_per_step)))
    config = agg.configuration('moon', monomer_mass=SEED_MASS, water_mm_h=number_flux * SEED_MASS * 3600.,
                               monomer=agg.LUNAR_MONOMER, stick=stick, band_speed=band[2], source_fraction=source_fraction,
                               duration=duration_h * 3600., depth=depth, dz=dz, dtc=dtc, xi=xi, nu_gam=nu_gam,
                               max_sps=400000)
    depth = config['grid_nml']['dom_top']
    config['grid_nml']['box_area'] = float(area)
    belt = source_fraction * depth
    per_step = max(1, int(round(number_flux * area * dtc / xi)))
    config['sp_nml']['nclrate'] = float((per_step + 1e-9) / (dtc * (belt - 1.) * area))
    end = config['general_nml']['time_end']
    config['general_nml'].update(dt_1dprof=float(end), dt_2dprof=float(end))
    config['sp_nml'].update(sp_habit=1, sp_kernel_id=2 if collisions else 0,
                            water_supply=-1. if saturated else water_mm_h / 3600., depo_max_fraction=float(depo_fraction),
                            band_mass=float(band[0]), band_axes=float(band[1]), ncl_until=-1., dt_snap=float(snap_h * 3600.),
                            plate_base=True, box_array=False)
    # Water saturation from the freezing level to the top: the 'LWC region' of McSnow's ideal profile.
    config['atmo_nml'].update(IGF_id=2, ssat=0., h0=0., h1=0., h2=float(depth))
    return config


# --- parity with the 24 October aggregation runs ---------------------------------------------

def parity(tag=None):
    """With every new option at its default, the build reproduces two 24 October runs: an Earth run in all of
    its output files, and a lunar run in its final superparticle state with the profile output trimmed.  Each
    build is tagged by its patch hash, so a changed build is checked afresh."""
    tag = tag or build_path('moon').name.split('_')[-1]
    out = {}
    for name, world, trim in (('earth_r1_s2', 'earth', False), ('moon_q3_f0.01_s2', 'moon', True)):
        reference = agg.CACHE / 'runs' / name
        old = json.loads((reference / 'run.json').read_text())
        config = json.loads(json.dumps(old['configuration']))
        if trim:
            end = config['general_nml']['time_end']
            config['general_nml'].update(dt_1dprof=float(end), dt_2dprof=float(end))
        label = f'parity_{tag}_' + name
        record = run(label, world, config, seed=old['seed'], keep_all=not trim)
        path = CACHE / 'runs' / label
        names = ['mass2fr.dat'] if trim else sorted(p.name for p in reference.glob('*.dat'))
        differ = [n for n in names if not (path / n).exists() or sha(path / n) != sha(reference / n)]
        out[name] = dict(world=world, build=build_path(world).name, compared_files=len(names), differing_files=differ,
                         trimmed_profiles=trim, identical=not differ and record['returncode'] == 0, wall_s=record['wall_s'])
        if not trim:
            for p in path.iterdir():
                if p.is_file() and p.name.startswith(TRIMMED):
                    p.unlink()
    # The build also reproduces one of this module's own columns made with the first build, in McSnow's final
    # superparticle state (mass2fr.dat, which recomputes every fall speed).  The snapshots differ only in the
    # unset fall speed of the two seeds created in each snapshot's step.
    first = json.loads((CACHE / 'runs/col_q3_f0.01_base/run.json').read_text())
    label = f'repro_{tag}_col_q3_f0.01_base'
    record = run(label, 'moon', first['configuration'], seed=first['seed'])
    same = sha(CACHE / 'runs' / label / 'mass2fr.dat') == sha(CACHE / 'runs/col_q3_f0.01_base/mass2fr.dat')
    mine, theirs = read_snapshots(label), read_snapshots('col_q3_f0.01_base')
    other = sum(int(np.count_nonzero(a[:, :5] != b[:, :5])) for (_, a), (_, b) in zip(mine, theirs)
                if a.shape == b.shape) if all(a.shape == b.shape for (_, a), (_, b) in zip(mine, theirs)) else None
    out['col_q3_f0.01_base'] = dict(world='moon', build=build_path('moon').name, earlier_executable=first['executable_sha256'],
                                    compared_files=1, identical=same and record['returncode'] == 0, wall_s=record['wall_s'],
                                    snapshot_entries_differing_outside_fall_speed=other)
    (CACHE / 'parity.json').write_text(json.dumps(out, indent=1) + '\n')
    return out


# --- parcels along the storm paths ---------------------------------------------------------

ROW_FIELDS = ('id', 'window_s', 'age_s', 'x_m', 'y_m', 'height_m', 'mass_kg', 'ice_volume_m3', 'aspect_ratio',
              'diameter_m', 'equatorial_dimension_m', 'polar_dimension_m', 'density_kg_m3', 'terminal_speed_m_s',
              'temperature_K', 'pressure_Pa', 'vapour_pressure_Pa', 'ice_rh', 'status')
BOX_DT = 600.          # s, the storm archive's cadence


def reference_concentration():
    """Mean crystal concentration in the water-limited reference column (3 mm/h shared by the reference flux):
    flux times descent time over depth, from the single-crystal column of plan step 3."""
    sizes = json.loads(SIZE_RANGE.read_text())
    ref = next(r for r in sizes['grid'] if r['water_mm_h'] == 3. and r['flux_factor'] == 1.)
    depth = sizes['column']['depth_m']
    return ref['number_flux_m2_s'] * ref['baseline']['descent_h'] * 3600. / depth


# --- columns that follow the storm paths -------------------------------------------------------------

COLUMN_TOP_M, COLUMN_BOTTOM_M, COLUMN_DZ = 46000., 24000., 500.     # m, in the archived storm's heights
N_SUB = int(round((COLUMN_TOP_M - COLUMN_BOTTOM_M) / COLUMN_DZ))
PATH_AREA = 20.                                                    # m2: the release cell holds 10^4 m3


def path_column_forcing(label, saturated):
    """Air for one sub-column per storm path: 500 m cells from 46 km down to 24 km, below every freezing level of
    the archive, at the probe's horizontal position every ten minutes while the probe lives and then carried by
    the archived wind at the probe's last height, until the archive ends.  Returns the binary table the build
    reads, the release cells and a summary."""
    from atmosphere.snow import mcsnow_weather as weather
    from atmosphere.snow import moist_history as history
    record = weather.verify_weather()
    data = np.load(weather.CACHE / 'weather/fields.npy', mmap_mode='r')
    z, dx = np.asarray(record['z_m']), record['dx_m']
    raw = np.load(STORM_RUNS / label / 'raw.npz')
    rows, col = raw['rows'], {k: i for i, k in enumerate(ROW_FIELDS)}
    n = len(raw['age_s'])
    window_end = (len(data) - 1) * 600.
    nt = int(round(window_end / BOX_DT)) + 1
    centres = COLUMN_TOP_M - COLUMN_DZ * (np.arange(N_SUB) + .5)
    T = np.empty((n, N_SUB, nt), np.float32)
    P = np.empty_like(T)
    V = np.empty_like(T)
    release = raw['release']
    span = window_end - release[:, 0]
    life = raw['age_s']
    tracks = []
    for k in range(n):
        r = rows[rows[:, col['id']] == k]
        r = r[np.argsort(r[:, col['age_s']], kind='stable')]
        age = r[:, col['age_s']]
        keep = np.concatenate([[True], np.diff(age) > 0])
        tracks.append((age[keep], r[keep, col['x_m']], r[keep, col['y_m']], r[keep, col['height_m']]))
    x = np.array([tr[1][0] for tr in tracks])
    y = np.array([tr[2][0] for tr in tracks])
    for i in range(nt):
        a = i * BOX_DT
        for k, (ages, xs, ys, zs) in enumerate(tracks):
            if a <= life[k]:
                x[k], y[k] = np.interp(a, ages, xs), np.interp(a, ages, ys)
        for start in np.unique(release[:, 0]):
            g = np.flatnonzero(release[:, 0] == start)
            seconds = min(start + a, window_end)
            pos = np.column_stack([np.repeat(x[g], N_SUB), np.repeat(y[g], N_SUB), np.tile(centres, len(g))])
            air = weather.weather_air(history.sample(data, pos, z, dx, seconds, 1))
            T[g, :, i] = air[0].reshape(len(g), N_SUB)
            P[g, :, i] = air[1].reshape(len(g), N_SUB)
            V[g, :, i] = -1. if saturated else air[2].reshape(len(g), N_SUB)
            later = g[a >= life[g]]
            if later.size:
                top = np.array([np.clip(tracks[k][3][-1], z[0], z[-1]) for k in later])
                wind = history.sample(data, np.column_stack([x[later], y[later], top]), z, dx, seconds, 1)
                x[later] += wind[:, 6] * BOX_DT
                y[later] += wind[:, 7] * BOX_DT
    release_cell = 1 + np.floor((COLUMN_TOP_M - release[:, 3]) / COLUMN_DZ)
    header = np.array([n * N_SUB, nt, BOX_DT, n], float)
    ends = np.repeat(span, N_SUB)
    blob = (header.tobytes() + ends.astype('<f8').tobytes() + release_cell.astype('<f8').tobytes()
            + T.reshape(n * N_SUB, nt).astype('<f4').tobytes() + P.reshape(n * N_SUB, nt).astype('<f4').tobytes()
            + V.reshape(n * N_SUB, nt).astype('<f4').tobytes())
    freezing = []
    for k in range(n):
        warm = T[k, :, 0] >= 273.15
        freezing.append(float(centres[np.argmax(warm)]) if warm.any() else None)
    summary = dict(paths=n, cells_per_column=N_SUB, cell_m=COLUMN_DZ, top_m=COLUMN_TOP_M, bottom_m=COLUMN_BOTTOM_M,
                   times=nt, bottom_cell_warmest_K=float(T[:, -1, :].min()),
                   release_cells=release_cell.tolist(), first_warm_cell_height_at_release_m=freezing)
    return blob, summary


def storm_condensate(label):
    """The storm's own cloud ice and snow along each path's sub-column, as kg per square metre, at release and as
    the mean over the archive window: the water a parcel population can be set against."""
    from atmosphere.snow import mcsnow_weather as weather
    from atmosphere.snow import moist_history as history
    from shared.constants import CM1_DRY_AIR_GAS_CONSTANT as Rd
    record = weather.verify_weather()
    data = np.load(weather.CACHE / 'weather/fields.npy', mmap_mode='r')
    z, dx = np.asarray(record['z_m']), record['dx_m']
    raw = np.load(STORM_RUNS / label / 'raw.npz')
    release = raw['release']
    centres = COLUMN_TOP_M - COLUMN_DZ * (np.arange(N_SUB) + .5)
    out = []
    for start in np.unique(release[:, 0]):
        g = np.flatnonzero(release[:, 0] == start)
        pos = np.column_stack([np.repeat(release[g, 1], N_SUB), np.repeat(release[g, 2], N_SUB), np.tile(centres, len(g))])
        values = []
        for seconds in np.arange(start, (len(data) - 1) * 600. + 1., 3600.):
            s = history.sample(data, pos, z, dx, float(seconds), 1)
            T = s[:, 0] * (s[:, 1] / history.t.CM1_P00) ** (Rd / history.t.CM1_CP)
            rho = s[:, 1] / (Rd * T)
            values.append(((s[:, 4] + s[:, 5]) * rho * COLUMN_DZ).reshape(len(g), N_SUB).sum(axis=1))
        values = np.array(values)
        out += [(int(k), float(values[0, j]), float(values[:, j].mean())) for j, k in enumerate(g)]
    out.sort()
    return np.array([[a, b] for _, a, b in out])      # per path: at release, mean over the window (kg/m2)


def path_column_config(n_paths, spans, concentration, *, stick=2, band=(1., 1., 1.), per_release=64, dtc=10.,
                       depo_fraction=.01, snap_s=1800.):
    depth = n_paths * N_SUB * COLUMN_DZ
    end = math.ceil((max(spans) + dtc) / 600.) * 600.
    xi = max(1, int(round(concentration * COLUMN_DZ * PATH_AREA / per_release)))
    config = agg.configuration('moon', monomer_mass=SEED_MASS, water_mm_h=1., monomer=agg.LUNAR_MONOMER, stick=stick,
                               band_speed=band[2], source_fraction=1., duration=end, dz=COLUMN_DZ, dtc=dtc, xi=xi,
                               max_sps=max(400000, 2 * n_paths * per_release))
    config['grid_nml'].update(dom_top=float(depth), box_area=float(PATH_AREA), nz=n_paths * N_SUB)
    config['general_nml'].update(time_end=float(end), dt_1dprof=float(end), dt_2dprof=float(end), interpol_atmos=False)
    # nrp0 only sizes McSnow's particle array in this mode (2 nrp0 V / xi0 entries); 1 leaves it at maxTotalNoSps.
    config['sp_nml'].update(nrp0=1, nclrate=float((n_paths * per_release + .5) / (dtc * depth * PATH_AREA)), sp_habit=1,
                            water_supply=-1., depo_max_fraction=float(depo_fraction), band_mass=float(band[0]),
                            band_axes=float(band[1]), ncl_until=float(dtc), dt_snap=float(snap_s), plate_base=True,
                            box_array=True, box_sub=N_SUB, box_settle=True, skip_bulk=True, depo_residual=RESIDUAL_MASS,
                            depo_censor=CENSOR_M)
    config['atmo_nml'].update(IGF_id=2, ssat=0., h0=0., h1=0., h2=float(depth), nh1=0., nh2=float(depth))
    return config


PATH_CASES = {f'{storm}_c{factor:g}': (storm, factor) for storm in ('archived', 'saturated') for factor in (1., .1, .01)}
PATH_VARIANTS = {'base': 64, 'sp4': 256}


def forcing_file(storm, edge):
    label = f'{storm}__{edge or "baseline"}'
    path = CACHE / 'forcing' / f'{label}.bin'
    if not path.exists():
        blob, summary = path_column_forcing(label, storm == 'saturated')
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(blob)
        summary.update(paths_sha256=sha(STORM_RUNS / label / 'raw.npz'), forcing_sha256=sha(path))
        path.with_suffix('.json').write_text(json.dumps(summary) + '\n')
    return path, json.loads(path.with_suffix('.json').read_text())


def path_case(case, variant='base', edge=None):
    storm, factor = PATH_CASES[case]
    forcing, summary = forcing_file(storm, edge)
    raw = np.load(STORM_RUNS / f'{storm}__{edge or "baseline"}' / 'raw.npz')
    spans = ((summary['times'] - 1) * BOX_DT - raw['release'][:, 0]).tolist()
    name = f'path_{case}_{variant}' + (f'_{edge}' if edge else '')
    config = path_column_config(summary['paths'], spans, reference_concentration() * factor, band=band_factors(edge),
                                per_release=PATH_VARIANTS[variant])
    record = run(name, 'moon', config, forcing_link=forcing)
    if record['returncode'] != 0:
        raise RuntimeError('McSnow failed: ' + name)
    return name, record


def read_arrivals(name):
    path = CACHE / 'runs' / name / 'arrivals.bin'
    if not path.exists():
        return np.zeros((0, 8))
    return np.fromfile(path, dtype='<f8').reshape(-1, 8)   # time, column, height, mass, span, monomers, multiplicity, reason


def path_summary(name):
    """Flakes reaching the freezing level from each path's column, and the largest flakes each column held."""
    run_record = json.loads((CACHE / 'runs' / name / 'run.json').read_text())
    sp = run_record['configuration']['sp_nml']
    n = run_record['configuration']['grid_nml']['nz'] // N_SUB
    arr = read_arrivals(name)
    t, k, zz, m, d, mono, xi, reason = arr.T if len(arr) else (np.zeros(0),) * 8
    warm = reason == 1
    w = xi[warm] * m[warm]
    pooled = dict(arrivals=int(warm.sum()), below_column=int((reason == 2).sum()))
    if warm.sum() >= 5:
        pooled.update(d50_mm=weighted_quantile(d[warm], w, .5) * 1e3, d90_mm=weighted_quantile(d[warm], w, .9) * 1e3,
                      d99_mm=weighted_quantile(d[warm], w, .99) * 1e3, largest_mm=float(d[warm].max() * 1e3),
                      mean_monomers=float((w * mono[warm]).sum() / w.sum()),
                      arrival_time_h=dict(p10=float(np.percentile(t[warm], 10) / 3600), p50=float(np.percentile(t[warm], 50) / 3600),
                                          p90=float(np.percentile(t[warm], 90) / 3600)))
    per_path = []
    for c in range(1, n + 1):
        s = warm & (k == c)
        if s.sum() >= 3:
            ws = xi[s] * m[s]
            per_path.append([c - 1, int(s.sum()), float(ws.sum()), weighted_quantile(d[s], ws, .5) * 1e3,
                             weighted_quantile(d[s], ws, .9) * 1e3, float(d[s].max() * 1e3)])
    peak = np.zeros(n)
    top = run_record['configuration']['grid_nml']['dom_top']
    for _, a in read_snapshots(name):
        z, mm, dd, mo, x, v = a.T
        cell = np.clip(((top - z) // COLUMN_DZ).astype(int), 0, n * N_SUB - 1)
        column = cell // N_SUB
        for c in np.unique(column):
            s = column == c
            if s.sum() >= 3:
                peak[c] = max(peak[c], weighted_quantile(dd[s], x[s] * mm[s], .9))
    reached = [r for r in per_path]
    def frac(limit):
        return float(np.mean([r[4] > limit for r in reached])) if reached else None
    return dict(pooled=pooled, columns_with_arrivals=len(reached),
                arriving_d90_mm=dict(p10=float(np.percentile([r[4] for r in reached], 10)) if reached else None,
                                     p50=float(np.percentile([r[4] for r in reached], 50)) if reached else None,
                                     p90=float(np.percentile([r[4] for r in reached], 90)) if reached else None),
                fraction_arriving_d90_over={f'{x * 100:g}_cm': frac(x * 1e3) for x in THRESHOLDS_M},
                column_peak_d90_mm=dict(p50=float(np.percentile(peak * 1e3, 50)), p90=float(np.percentile(peak * 1e3, 90)),
                                        max=float(peak.max() * 1e3)),
                per_path=[[round(v, 4) if isinstance(v, float) else v for v in r] for r in per_path])


# --- reading the runs ----------------------------------------------------------------------------

def read_snapshots(name):
    raw = np.fromfile(CACHE / 'runs' / name / 'snapshots.bin', dtype='<f8')
    out, i = [], 0
    while i < raw.size:
        t, n = raw[i], int(raw[i + 1])
        i += 2
        out.append((float(t), raw[i:i + 6 * n].reshape(n, 6).copy()))
        i += 6 * n
    return out          # columns: height, mass, maximum dimension, monomers, multiplicity, fall speed


def weighted_quantile(values, weights, q):
    return agg.weighted_quantile(np.asarray(values), np.asarray(weights), q)


THRESHOLDS_M = (.10, .14)       # the load analogy's survivable span, lunar equivalents of 4-6 cm Earth flakes
POOL_H = 12.                    # h: steady-state snapshots pooled for each column's statistics


def column_summary(name, *, bin_m=1000., arriving_m=2000., snapshot=-1, pool_h=0.):
    """Sizes by depth and the snow arriving at the freezing level, weighted by mass flux (multiplicity x mass x
    fall speed), from one binary snapshot or from every snapshot in the last pool_h hours up to it, pooled with
    equal weight; and the depth below the top where the 90th-percentile span first exceeds each threshold."""
    run_record = json.loads((CACHE / 'runs' / name / 'run.json').read_text())
    grid = run_record['configuration']['grid_nml']
    top, area = grid['dom_top'], grid['box_area']
    snaps = read_snapshots(name)
    t = snaps[snapshot][0]
    chosen = [a for s, a in snaps if t - pool_h * 3600. - 1. <= s <= t]
    a = np.concatenate(chosen)
    z, m, d, mono, xi, v = a.T
    xi = xi / len(chosen)
    # Superparticles created in the step a snapshot is written have no fall speed yet (the field is unset);
    # they are 4 um seeds and carry no weight.
    v = np.where(np.isfinite(v) & (v > 1e-6) & (v < 50.), v, 0.)
    w = xi * m * v
    bins = []
    for upper in np.arange(top, 0, -bin_m):
        s = (z < upper) & (z >= upper - bin_m)
        if s.sum() >= 3 and w[s].sum() > 0:
            bins.append(dict(depth_below_top_km=float((top - upper + bin_m / 2) / 1000), superparticles=int(s.sum()),
                             d50_mm=weighted_quantile(d[s], w[s], .5) * 1e3, d90_mm=weighted_quantile(d[s], w[s], .9) * 1e3,
                             d99_mm=weighted_quantile(d[s], w[s], .99) * 1e3, largest_mm=float(d[s].max() * 1e3),
                             mean_monomers=float((w[s] * mono[s]).sum() / w[s].sum()),
                             ice_kg_m3=float((xi[s] * m[s]).sum() / (bin_m * area))))
    onset = {f'{x * 100:g}_cm': next((b['depth_below_top_km'] for b in bins if b['d90_mm'] > x * 1e3), None)
             for x in THRESHOLDS_M}
    low = z < arriving_m
    arriving = dict(layer_m=arriving_m, superparticles=int(low.sum()))
    if low.sum() >= 5 and w[low].sum() > 0:
        melted = (6 * m[low] / (math.pi * 1000.)) ** (1 / 3)
        arriving.update(water_mm_h=float(w[low].sum() / (arriving_m * area) * 3600.),
                        d50_mm=weighted_quantile(d[low], w[low], .5) * 1e3,
                        d90_mm=weighted_quantile(d[low], w[low], .9) * 1e3,
                        d99_mm=weighted_quantile(d[low], w[low], .99) * 1e3,
                        largest_mm=float(d[low].max() * 1e3),
                        melted_d50_mm=weighted_quantile(melted, w[low], .5) * 1e3,
                        mean_monomers=float((w[low] * mono[low]).sum() / w[low].sum()),
                        aggregate_mass_fraction=float(w[low][mono[low] > 1.01].sum() / w[low].sum()))
    return dict(time_h=t / 3600., pooled_snapshots=len(chosen), column_superparticles=int(len(z) / len(chosen)),
                exhausted_fraction=float(np.mean(xi * len(chosen) <= 2)), bins=bins, onset_depth_km=onset,
                arriving=arriving)




def growth_parity():
    """Crystals growing without collisions in the water-limited reference column (3 mm/h shared by the reference
    flux) against the single-crystal column program of plan step 3, which gives 3.14 mm after 54.2 hours."""
    config = column_config(3., 1., collisions=False, source_fraction=.01, nu_gam=1e4, dtc=60., sp_per_step=1.,
                           duration_h=72., snap_h=6.)
    record = run('growth_parity_q3_f1', 'moon', config)
    snaps = read_snapshots('growth_parity_q3_f1')
    t, a = snaps[-1]
    z, m, d, mono, xi, v = a.T
    low = z < 300.
    sizes = json.loads(SIZE_RANGE.read_text())
    ref = next(r for r in sizes['grid'] if r['water_mm_h'] == 3. and r['flux_factor'] == 1.)['baseline']
    w = xi[low] * m[low] * v[low]
    out = dict(returncode=record['returncode'], wall_s=record['wall_s'], time_h=t / 3600., superparticles=int(len(z)),
               arriving_superparticles=int(low.sum()),
               arriving_d50_mm=weighted_quantile(d[low], w, .5) * 1e3, arriving_d10_mm=weighted_quantile(d[low], w, .1) * 1e3,
               arriving_d90_mm=weighted_quantile(d[low], w, .9) * 1e3,
               arriving_mass_median_mg=weighted_quantile(m[low], w, .5) * 1e6,
               column_program_diameter_mm=ref['diameter_m'] * 1e3, column_program_mass_mg=ref['mass_kg'] * 1e6,
               column_program_descent_h=ref['descent_h'],
               arriving_water_mm_h=float(w.sum() / (300. * config['grid_nml']['box_area']) * 3600.))
    (CACHE / 'growth_parity.json').write_text(json.dumps(out, indent=1) + '\n')
    return out


# --- the cases ----------------------------------------------------------------------------------

COLUMN_CASES = {
    'q3_f1': dict(water_mm_h=3., flux_factor=1., duration_h=72.),
    'q3_f0.1': dict(water_mm_h=3., flux_factor=.1, duration_h=48.),
    'q3_f0.01': dict(water_mm_h=3., flux_factor=.01, duration_h=36.),
    'q1_f1': dict(water_mm_h=1., flux_factor=1., duration_h=96.),
    'q0.3_f1': dict(water_mm_h=.3, flux_factor=1., duration_h=132.),
    'q0.3_f0.1': dict(water_mm_h=.3, flux_factor=.1, duration_h=72.),
    'sat_f1': dict(water_mm_h=0., flux_factor=1., saturated=True, duration_h=36.),
    'sat_f0.1': dict(water_mm_h=0., flux_factor=.1, saturated=True, duration_h=36.),
    'sat_f0.01': dict(water_mm_h=0., flux_factor=.01, saturated=True, duration_h=36.),
    # The snow zone of the design-climate storm boxes reaches far above -15 C (snow_zone_depth.py): the same
    # light snow from the habit table's coldest level, -30 C, 32.6 km above the freezing level.
    'deep_q0.3_f1': dict(water_mm_h=.3, flux_factor=1., duration_h=168., top_C=-30.),
    'deep_q0.3_f0.1': dict(water_mm_h=.3, flux_factor=.1, duration_h=96., top_C=-30.),
}
# Numerical variants: the sampled column area and superparticles per step together (same multiplicity), or
# superparticles alone (smaller multiplicity).
VARIANTS = {'base': dict(area=5., sp_per_step=2.), 'volume4': dict(area=20., sp_per_step=8.),
            'sp2': dict(area=5., sp_per_step=4.), 'sp4': dict(area=5., sp_per_step=8.),
            'half': dict(area=5., sp_per_step=1.)}
BAND_EDGES = ('combined_small', 'combined_large')


def band_factors(edge):
    if edge is None:
        return (1., 1., 1.)
    return tuple(json.loads(EARTH_BAND.read_text())['column']['allowance_F'][edge]['factors'])


def column_case(case, variant='base', edge=None):
    options = dict(COLUMN_CASES[case], **VARIANTS[variant])
    name = f'col_{case}_{variant}' + (f'_{edge}' if edge else '')
    config = column_config(band=band_factors(edge), snap_h=6., **options)
    record = run(name, 'moon', config)
    if record['returncode'] != 0:
        raise RuntimeError('McSnow failed: ' + name)
    return name, record




def path_ice_peak(name):
    """Each path column's largest ice load through the run, kg per square metre."""
    record = json.loads((CACHE / 'runs' / name / 'run.json').read_text())
    grid = record['configuration']['grid_nml']
    n = grid['nz'] // N_SUB
    peak = np.zeros(n)
    for _, a in read_snapshots(name):
        z, m, d, mono, xi, v = a.T
        column = ((grid['dom_top'] - z) // COLUMN_DZ).astype(int) // N_SUB
        peak = np.maximum(peak, np.bincount(np.clip(column, 0, n - 1), weights=xi * m, minlength=n) / grid['box_area'])
    return peak


def convergence(base, refined):
    """Base against a refined variant: 90th-percentile span by depth, onset depths and the arriving snow."""
    a, b = column_summary(base, pool_h=POOL_H), column_summary(refined, pool_h=POOL_H)
    rows = []
    for x in a['bins']:
        y = next((r for r in b['bins'] if r['depth_below_top_km'] == x['depth_below_top_km']), None)
        if y:
            rows.append([x['depth_below_top_km'], round(x['d90_mm'], 2), round(y['d90_mm'], 2),
                         round(y['d90_mm'] / x['d90_mm'] - 1, 4)])
    keys = ('water_mm_h', 'd50_mm', 'd90_mm', 'd99_mm', 'largest_mm')
    arriving = {k: [a['arriving'].get(k), b['arriving'].get(k)] for k in keys}
    # the deepest bin above the first two successive bins whose 90th percentiles differ by more than 25%
    agree = rows[-1][0] if rows else 0.
    for i in range(len(rows) - 1):
        if abs(rows[i][3]) > .25 and abs(rows[i + 1][3]) > .25:
            agree = rows[i - 1][0] if i > 0 else 0.
            break
    return dict(base=base, refined=refined, superparticles=[a['column_superparticles'], b['column_superparticles']],
                d90_by_depth=rows, onset_depth_km=[a['onset_depth_km'], b['onset_depth_km']], arriving=arriving,
                d90_agrees_to_km=agree, agreement_rule='successive 1 km bins differing by more than 25% end the agreement')


def assemble():
    """Compact product from every completed run."""
    runs = sorted(p.name for p in (CACHE / 'runs').iterdir() if (p / 'run.json').exists())
    columns, parcels = {}, {}
    for name in runs:
        record = json.loads((CACHE / 'runs' / name / 'run.json').read_text())
        sp = record['configuration']['sp_nml']
        if name.startswith('col_'):
            snaps = read_snapshots(name)
            final = column_summary(name, pool_h=POOL_H)
            previous = column_summary(name, snapshot=-3, pool_h=POOL_H) if len(snaps) > 4 else None
            columns[name] = dict(wall_s=record['wall_s'], executable_sha256=record['executable_sha256'], water_supply_mm_h=sp['water_supply'] * 3600. if sp['water_supply'] >= 0 else None,
                                 multiplicity=sp['xi0'], area_m2=record['configuration']['grid_nml']['box_area'],
                                 band=[sp['band_mass'], sp['band_axes'], sp['band_speed']], final=final,
                                 previous=dict(time_h=previous['time_h'], onset_depth_km=previous['onset_depth_km'],
                                               arriving=previous['arriving']) if previous else None)
        elif name.startswith('path_'):
            storm = 'saturated' if name.startswith('path_saturated') else 'archived'
            edge = next((e for e in BAND_EDGES if name.endswith(e)), None)
            label = f'{storm}__{edge or "baseline"}'
            summary = path_summary(name)
            condensate = storm_condensate(label)
            peak = path_ice_peak(name)
            n_release = round(math.floor(sp['nclrate'] * sp['dtc'] * record['configuration']['grid_nml']['dom_top']
                                         * record['configuration']['grid_nml']['box_area']) / len(peak))
            parcels[name] = dict(wall_s=record['wall_s'], executable_sha256=record['executable_sha256'],
                                 superparticles_per_release=n_release,
                                 concentration_m3=sp['xi0'] * n_release / (COLUMN_DZ * PATH_AREA), multiplicity=sp['xi0'],
                                 band=[sp['band_mass'], sp['band_axes'], sp['band_speed']],
                                 forcing_sha256=record['forcing_sha256'],
                                 summary={k: v for k, v in summary.items() if k != 'per_path'},
                                 ice_peak_kg_m2=dict(p50=float(np.median(peak)), p90=float(np.percentile(peak, 90))),
                                 storm_ice_and_snow_kg_m2=dict(release_p50=float(np.median(condensate[:, 0])),
                                                               window_mean_p50=float(np.median(condensate[:, 1]))))
    product = dict(
        schema='terluna.atmosphere.snow-mcsnow-growth-aggregation/1',
        evidence=('Numerical McSnow columns and storm-path columns in which crystals grow by the frozen improved '
                  'component and aggregate together, without breakup; Earth-checked aggregation (mcsnow_aggregation) '
                  'and a growth parity with the single-crystal column program; no empirical lunar validation.'),
        reading_rule=('Sizes are maximum dimensions weighted by mass flux (reference columns) or mass (arrivals). A '
                      'size below the top few kilometres is a result only where its refined variant agrees '
                      '(convergence); columns pool their last 12 hours.'),
        thresholds_m=list(THRESHOLDS_M), reference_concentration_m3=reference_concentration(), seed_diameter_m=SEED_DIAMETER,
        parity=json.loads((CACHE / 'parity.json').read_text()),
        growth_parity=json.loads((CACHE / 'growth_parity.json').read_text()),
        builds={w: dict(executable_sha256=b['executable_sha256'], patch_sha256=b['patch_sha256'], revision=b['revision'])
                for w, b in ((w, verify_build(w)) for w in ('earth', 'moon'))},
        columns=columns, parcels=parcels,
        convergence=[convergence(f'col_{c}_base', f'col_{c}_{v}') for c in COLUMN_CASES for v in ('sp4', 'sp2', 'volume4', 'half')
                     if (CACHE / 'runs' / f'col_{c}_{v}' / 'run.json').exists()],
        inputs=dict(earth_error_band_sha256=sha(EARTH_BAND), lunar_size_range_sha256=sha(SIZE_RANGE)),
        producer_sha256=sha(Path(__file__)), raw_directory=str(CACHE.relative_to(ROOT)))
    PRODUCT.write_text(json.dumps(product, separators=(',', ':'), allow_nan=False) + '\n')
    return product


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('build', 'parity', 'growth_parity', 'column', 'paths', 'assemble'))
    parser.add_argument('--case')
    parser.add_argument('--variant', default='base')
    parser.add_argument('--edge', default=None)
    args = parser.parse_args()
    if args.action == 'build':
        for world in ('earth', 'moon'):
            print(world, build(world)['executable_sha256'])
    elif args.action == 'parity':
        print(json.dumps(parity(), indent=1))
    elif args.action == 'growth_parity':
        print(json.dumps(growth_parity(), indent=1))
    elif args.action == 'column':
        name, record = column_case(args.case, args.variant, args.edge)
        summary = column_summary(name)
        print(name, 'wall', round(record['wall_s']), 's; onset', summary['onset_depth_km'], '; arriving',
              {k: round(v, 2) if isinstance(v, float) else v for k, v in summary['arriving'].items()}, flush=True)
    elif args.action == 'paths':
        name, record = path_case(args.case, args.variant, args.edge)
        summary = path_summary(name)
        print(name, 'wall', round(record['wall_s']), 's; paths delivering snow', summary['columns_with_arrivals'],
              '; their arriving 90th percentiles', summary['arriving_d90_mm'], flush=True)
    elif args.action == 'assemble':
        product = assemble()
        print(len(product['columns']), 'columns;', len(product['parcels']), 'path-column runs')


if __name__ == '__main__':
    main()
