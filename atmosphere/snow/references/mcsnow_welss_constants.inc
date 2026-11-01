!
! Copyright (C) 2004-2024, DWD
! See ./AUTHORS.txt for a list of authors
! See ./LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
!

MODULE mo_constants

  IMPLICIT NONE

  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6,37) !< single precision
  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12,307) !< double precision
  INTEGER, PARAMETER :: i2 = SELECTED_INT_KIND(4)   !< at least 2 byte integer
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)   !< at least 4 byte integer
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)   !< at least 8 byte integer
  !INTEGER, PARAMETER :: i8 = i4
  INTEGER, PARAMETER :: mp = wp ! enable mixed precision by switching to mp = sp

  REAL(wp), PARAMETER :: pi = ATAN(1._wp) * 4._wp

  REAL(wp), PARAMETER :: z1_6   = 1.0_wp/6.0_wp
  REAL(wp), PARAMETER :: z1_3   = 1.0_wp/3.0_wp
  REAL(wp), PARAMETER :: z2_3   = 2.0_wp/3.0_wp

  REAL(wp), PARAMETER :: pi_6   = pi/6.0_wp
  REAL(wp), PARAMETER :: pi_4   = pi/4.0_wp
  REAL(wp), PARAMETER :: pi_2   = pi/2.0_wp 
  
  REAL(wp), PARAMETER :: z4_pi  = 4.0_wp/pi
  REAL(wp), PARAMETER :: z6_pi  = 6.0_wp/pi

  REAL(wp), PARAMETER :: z4pi_3 = 4._wp*pi/3._wp
  REAL(wp), PARAMETER :: z3_4pi = 3._wp/(4._wp*pi)

  REAL(wp), PARAMETER :: zs1_2pi = 1/SQRT(2.0_wp*pi)


END MODULE mo_constants


!
! Copyright (C) 2004-2024, DWD
! See ./AUTHORS.txt for a list of authors
! See ./LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
!

MODULE mo_atmo_types
  USE mo_constants,    ONLY: wp, pi

  IMPLICIT NONE

  TYPE :: t_atmo
    REAL(wp) :: &
      p,        &  ! [Pa]    pressure
      pv,       &  ! [Pa]    vapor pressure
      qv,       &  ! [kg/kg] specific humidity of vapor
      qc,       &  ! [kg/kg] specific humidity of cloud liquid water
      rho,      &  ! [kg/m3] density
      T,        &  ! [K]     temperature
      eta,      &  ! [kg/ms] dynamic viscosity
      etai,     &  ! [ms/kg] inverse of dynamic viscosity
      ssat,     &  ! [1]     supersaturation over ice (ssat=pv/e_sat,i-1)
      rh,       &  ! [%]     rel. humidity
      psatw,    &
      psati
  END TYPE t_atmo

  TYPE :: t_diag
     REAL(wp) :: &
          rime_rate,  & ! riming growth rate (total)
          dep_rate,   & ! depositional growth rate (total)
          dep_mono,   & ! depositional growth rate of mono-crystals
          heat_rate     ! warming cooling due to latent heat
  END type t_diag

  ! constants
  !
  REAL(wp), PARAMETER ::    &
    ps    = 101325._wp,     &  ! [Pa]    surface pressure
    grav  = 9.81_wp,        &  ! [m/s2]  gravitational acceleration
    tlr   = 0.0062_wp,      &  ! [K/m]   temp. lapse rate
    eta0  = 1.718e-5_wp,    &  ! [kg/ms] dynamic viscosity of dry air @ T=273 K
    Rd    = 287.102_wp,     &  ! [J/kgK] gas constant of dry air
    Rv    = 461.401_wp,     &  ! [J/kgK] gas constant of cloud water vapor
    c_pa  = 1005_wp,        &  ! [J/kgK] specific heat capacity of dry air
    c_pi  = 2.108e3_wp,     &  ! [J/kgK] specific heat capacity of ice
    c_pw  = 4.205e3_wp,     &  ! [J/kgK] specific heat capacity of water
    rhlr  = 0.05_wp,        &  ! [%/m]   rel. humidity lapse rate
    rhoi  = 919_wp,         &  ! [kg/m3] density of ice (916.2@0°C, 918.9@-10°C,
                               !         919.4@-20°C, 920@-30°C)
    rhol  = 1e3_wp
  

  REAL(wp), PARAMETER ::    &
    A_i  = 2.18745584e1_wp, &  !..constant in saturation pressure - ice
    A_w  = 1.72693882e1_wp, &  !..constant in saturation pressure - liquid
    B_i  = 7.66000000e0_wp, &  !..constant in saturation pressure - ice
    B_w  = 3.58600000e1_wp, &  !..constant in saturation pressure - liquid
    T_3  = 273.15_wp,       &  ! [K]     melting temperature of ice/snow
    e_3  = 6.10780000e2_wp     !..saturation pressure at T = T_3


  ! derived constants
  REAL(wp), PARAMETER ::       &
    g_Rtlr = grav/(Rd*tlr),    &
    rhoii  = 1._wp/rhoi,       & 
    T_3i   = 1._wp/T_3,         &
    rholi  = 1._wp/rhol,       &
    m2r    = rholi*0.75_wp/pi, &
    r2m    = 1._wp/m2r

END MODULE mo_atmo_types

module geometry_constants
use mo_constants
use mo_atmo_types
implicit none
  ! m-D and A-D constants from 
  !   Morrison and Grabowski 'A Novel Approach...' (2007)
  !
  ! NOTE: coeffs in Table 1, p. 1532 are given in CGS, 
  REAL(wp), PARAMETER ::                       &
!    These are 'crystal with sector-like branches P1b' of M96
!    unr_bet = 2.02_wp,                         &
!    unr_alf = 0.00142_wp*10_wp**(2_wp*unr_bet-3_wp), & ! 10_wp**.. to scale to MKS
!    unr_gam = 1.97_wp,                         &
!    unr_sig = 0.55_wp*10_wp**(2_wp*unr_gam-4_wp),  &
!   Jussi's aggregates of dendrites
    ! unr_alf = 0.01243_wp,          &
    ! unr_bet = 2.00000_wp,          &
    ! unr_sig = 0.05625_wp,          &
    ! unr_gam = 1.81000_wp,          &
!   Aggregates of side planes, columns and bullets of M96
    unr_bet = 2.1_wp,                         &
    unr_alf = 0.0028_wp*10_wp**(2_wp*unr_bet-3_wp), & ! 10_wp**.. to scale to MKS
    unr_gam = 1.88_wp,                         &
    unr_sig = 0.2285_wp*10_wp**(2_wp*unr_gam-4_wp),  &
    sph_gam = 2.0_wp,                          &
    sph_sig = pi_4,                            &
    gra_gam = 2.0_wp,                          &
    gra_sig = pi_4,                            &
    ! derived m-D and A-D constants
    unr_alf_i   = 1/unr_alf,           &
    unr_bet_i   = 1/unr_bet,           &
    betm3       = unr_bet - 3._wp,     &
    betm3_i     = 1/betm3,             &
    Dth         = (rhoi * pi_6 * unr_alf_i)**betm3_i, &   ! Eq. (8)
    mth         = pi_6 * rhoi * Dth**3, &
    z1_6alf     = z1_6*unr_alf_i,      &
    pi_6alf     = pi*z1_6*unr_alf_i,   &
    z6alf_pi    = 6*unr_alf/pi,        &
    z3_bet      = 3*unr_bet_i

  ! m-D and A-D constants for m2d_dda
  !   by A. Seifert and J. Leinonen
  REAL(wp), PARAMETER ::               &
    unr_alf_dda = 0.01243_wp,          &
    unr_bet_dda = 2.00000_wp,          &
    unr_sig_dda = 0.05625_wp,          &
    unr_gam_dda = 1.81000_wp,          &
    a_m = 0.72_wp, &
    p_m = 0.52_wp, &
    ! derived m-D and A-D constants
    unr_alf_dda_i = 1._wp/unr_alf_dda, &
    unr_bet_dda_i = 1._wp/unr_bet_dda, &
    unr_sig_dda_i = 1._wp/unr_sig_dda, &
    unr_gam_dda_i = 1._wp/unr_gam_dda, &
    p_m_i = 1._wp / p_m, &
    a_m_i = 1._wp / a_m
  ! properties of plate in Jussis aggregation-model   
    REAL(wp), PARAMETER  :: &
    mono_bet_plates = 2.47,                            &
    mono_alf_plates = 0.7453,                          &
    mono_gam_plates = 2.00,                            &
    mono_sig_plates = 0.6495,                          &
    mono_epsilon_plates = 0.474-1.0,                   &
    mono_zeta_plates    = 0.00125058

end module geometry_constants
