!
! Copyright (C) 2004-2024, DWD
! See ./AUTHORS.txt for a list of authors
! See ./LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
!

REAL(wp) FUNCTION boehmVel( atmo, sp, partInfo ) RESULT( vterm_bohm )
  TYPE(t_atmo),                     INTENT(IN)    :: atmo
  TYPE(t_sp),                       INTENT(IN)    :: sp
  TYPE(t_CE_boehm), OPTIONAL, INTENT(OUT)   :: partInfo

  !REAL(wp) :: ALPHA 
  !REAL(wp) :: ZAEHLER, NENNER, XXX, CDP, CDPRIM, BETA, NRE0, CDO, GAM, NRE, DEL0, X0
  !REAL(wp) :: alf, beta, logm
  REAL(wp) :: fx
  REAL(wp) :: q, k, gama_big, N_Re, DEL0, d, mtot, X, X_0, V_tot, r_a, r_c
  REAL(wp) :: k_pro, gama_big_pro, N_Re_pro, vterm_bohm_pro, X_pro
  
  mtot = sp%m_w + sp%m_r + sp%m_f + sp%m_i

  IF( BTEST(sp%statusb,sp_stat%liquidOut)  ) THEN
  !IF( BTEST(sp%statusb,sp_stat%liquidSurf) ) THEN   ! smooth
    X_0 = 6.7D6
  ELSE                                              ! rough
    X_0 = 2.8D6
  ENDIF

  ! characteristic length is defined regarding the radius/a-axis, transformation needed for prolates
  IF (sp%phi .GT. 1._wp) THEN   
    IF (sp%mm < 1.01_wp) THEN
      q = z4_pi
      IF ( (sp_habit == 1 .AND. sp%mm < 1.01_wp) .OR. (sp_habit == 2) ) THEN
        CALL spheroid_geometry(sp, r_a, r_c, V_tot)
      ELSE
        V_tot = pi_6 * sp%d**3 ! TODO what about melting particles, they might not be spherical
      ENDIF
      d = (4.0_wp * V_tot / pi / sp%phi)**(1.0_wp/3.0_wp)
    ELSE
      q = sp%p_area / (pi_4 * sp%d**2 / sp%phi)
      d = sp%d / sp%phi
    END IF                                  
  ELSE
    q = sp%p_area / (pi_4 * sp%d**2)
    d = sp%d
  END IF

  ! cylindrical assumption q = 4/pi
  N_Re       = Reynolds_Boehm(atmo, sp%phi, mtot, q, X_0, k, gama_big, X)
  vterm_bohm = N_Re * atmo%eta / (d * atmo%rho)         ! I (22)
  ! possible interpolation between prolate and cylindrical geometry.
  ! TRAIL data (McCorquodale & Westbrook 2020) show transition instead of sudden jump
  IF ( (sp%phi .GT. 1._wp) .AND. (prolate_intp) ) THEN
    ! prolate assumption q = A / A_ce 
    N_Re_pro       = Reynolds_Boehm(atmo, sp%phi, mtot, sp%p_area / (pi_4 * sp%d**2 / sp%phi),    &
                              X_0, k_pro, gama_big_pro, X_pro)
    vterm_bohm_pro = N_Re_pro * atmo%eta / ( sp%d / sp%phi * atmo%rho)

    !logm        = LOG(mtot)
    ! alf         = 8.60e-2_wp * logm + 1.722 
    ! beta        = 3.08e-2_wp * logm + 1.691  

    ! fx          = MIN(1._wp, beta * exp(-alf * sp%phi))         ! limiter to prevent overestimation
    fx          =  exp( -0.3_wp * (sp%phi-1.0))                   ! [0.16, 0.0094] complex fit
    vterm_bohm  = fx * vterm_bohm_pro + (1._wp - fx) * vterm_bohm
  END IF


  IF( PRESENT( partInfo ) ) THEN

    IF( BTEST(sp%statusb,sp_stat%liquidSurf) ) THEN   ! smooth
      partInfo%SIG       = 0.0_wp
    ELSE                                              ! rough
      partInfo%SIG       = 0.25_wp
    ENDIF
    IF( sp%phi > 1.0_wp ) THEN     ! flow is not axisymmetric
      partInfo%j = 1
    ELSE
      partInfo%j = 2
    ENDIF

    partInfo%r       = 0.5_wp * d                                         ! 92c, p.106 "radius of drop or cylinder"
    IF (sp%phi .GT. 1._wp) THEN
      partInfo%r_e     = partInfo%r * SQRT(z4_pi * sp%phi )               ! B99 p. 170 definition in text
    ELSE
      partInfo%r_e     = SQRT( sp%p_area / pi )                           ! 99 eq. 3
    END IF
    partInfo%v       = vterm_bohm
    partInfo%m       = mtot
    partInfo%alpha   = sp%phi
    partInfo%alpha_e = MIN(sp%phi, 1.0_wp) * partInfo%r/partInfo%r_e      ! 99 eq. 4
    partInfo%q_e     = q
    IF( q > 1.0_wp ) THEN                                                 ! 99 eq. 5
      partInfo%q_e    = q * pi / 4.0_wp
    END IF
    partInfo%k       = k

    ! Berechnung weiterer Groessen fuer die Bestimmung der Effizienz
    IF( sp%phi .LE. 1.0_wp ) THEN ! 99 chapter 2.2
      DEL0=3.60_wp
    ELSE
      DEL0=4.54_wp
    ENDIF
    partInfo%DELS   = DEL0*partInfo%r/SQRT(N_Re*gama_big)                 ! 99 eq. 1
    partInfo%CDI    = X/N_Re**2 - 24.0_Wp*partInfo%k/N_Re                 ! 99 eq. 18

  ENDIF
END FUNCTION boehmVel

! inversion of the Re-Cd relation assumed by Boehm '92 including turbulence correction and low
! Reynolds number correction
REAL(wp) FUNCTION Reynolds_Boehm(atmo, alpha, mtot, q, X_0, k, gama_big, X) RESULT(N_Re)
  TYPE(t_atmo), INTENT(IN)    :: atmo
  REAL(wp), INTENT(IN)        :: alpha      ! aspect ratio, original naming of Boehm
  REAL(wp), INTENT(IN)        :: mtot
  REAL(wp), INTENT(IN)        :: q
  REAL(wp), INTENT(IN)        :: X_0
  REAL(wp), INTENT(OUT)       :: k
  REAL(wp), INTENT(OUT)       :: gama_big
  REAL(wp), INTENT(OUT)       :: X

  REAL(wp) :: C_DP, C_DP_prim, beta, gama_small, C_DO, N_Re0

  CALL hydrodynamicValuesBoehm( alpha, q, k, gama_big, C_DP )

  ! 99 (7) or I (17)/(18)
  X = 8.0_wp*mtot*grav*atmo%rho/(pi*(atmo%eta**2)*MAX(alpha,1.0_wp)*q**(1.0_wp/4.0_wp))                 ! porosity correction also for q > 1
  !X = 8.0_wp*mtot*grav*atmo%rho/(pi*(atmo%eta**2)*MAX(alpha,1.0_wp)*MAX(q**(1.0_wp/4.0_wp),q))         ! original Version, Boehm Paper Version
  !X = 8.0_wp*mtot*grav*atmo%rho/(pi*(atmo%eta**2)*MAX(alpha,1.0_wp)*MIN(q**(1.0_wp/4.0_wp),1.0_wp))    ! Best fit to B92Fig3

  ! I (23) turbulence correction
  C_DP_prim = C_DP*((1.0_wp+1.6_wp*(X/X_0)**2)/(1.0_wp+(X/X_0)**2))

  ! I (21)
  beta = SQRT(1.0_wp+C_DP_prim/6.0_wp/k*SQRT(X/C_DP_prim))-1
  ! I (20)
  N_Re0 = 6.0_wp*k/C_DP_prim*beta**2

  ! low Reynolds number
  ! I (16)
  C_DO = 4.5_wp*k**2*MAX(alpha,1.0_wp)
  ! I (26)
  gama_small = (C_DO - C_DP) / (4.0_wp * C_DP)
  ! I (27)
  N_Re  = N_Re0*(1.0_wp + (2.0_wp*beta*EXP(-beta*gama_small))/((2.0_wp+beta)*(1.0_wp+beta)) )

END FUNCTION Reynolds_Boehm

!
! helper function for Boehm efficiency
!
SUBROUTINE hydrodynamicValuesBoehm( phi, q, k, gama_big, C_DP )
  REAL(wp), INTENT(in)  :: phi
  REAL(wp), INTENT(in)  :: q
  REAL(wp), INTENT(out) :: k
  REAL(wp), INTENT(out) :: gama_big
  REAL(wp), INTENT(out) :: C_DP

  ! 99 (9) attention sqrt(phi) should be just alpha, see I (11)/(12)/(13)/(14) | alpha = phi
  k = min(max(0.82_wp+0.18_wp*phi,0.85_wp),0.37_wp+0.63_wp/phi, &
    & 1.33_wp/(max(log(phi),0.0_wp)+1.19_wp))
  ! 99 (11) or I (5)
  gama_big = max(1.0_wp, min(1.98_wp,3.76_wp-8.41_wp*phi+9.18*phi**2-3.53*phi**3))
  ! 99 (10) or I (7)/(8)
  C_DP = max(0.292_wp*k*gama_big,0.492_wp-0.2_wp/sqrt(phi))
  C_DP = max(1.0_wp,q*(1.46_wp*q-0.46_wp))*C_DP
END SUBROUTINE hydrodynamicValuesBoehm

!
! Copyright (C) 2004-2024, DWD
! See ./AUTHORS.txt for a list of authors
! See ./LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
!

FUNCTION Ecbohm( atmo, partInfoX, partInfoY ) RESULT( Ecbohm_v )
  TYPE(t_atmo), INTENT(IN)           :: atmo
  TYPE(t_CE_boehm), TARGET, INTENT(IN) :: partInfoX, partInfoY
  REAL(wp)                           :: Ecbohm_v


  TYPE(t_CE_boehm), POINTER :: partInfo1, partInfo2
  REAL(wp) :: rhoa, nu
  REAL(wp) :: jstar, rr1, re1, v1, m1, alf1, qe1, & 
       &      dels1, cdi1, vi1, rr2, re2, V2, m2, &
       &      alf2, qe2, kk2, dels2, sig1, sig2
  REAL(wp) :: vmin, mstar, rstar, astar, phi1, arg1, phi2, arg2, vstar, delstar, alfstar, qstar,  &
       &      kstar, gamstar, nrestar, cdpstar, cdostar, betaexp, gammaexp, cdistar,              &
       &      gg, ff, hh, ffhh, bb, k2x, k2y, cx, cy, delx, dely, tdel, ep

  nu = atmo%eta
  rhoa = atmo%rho
  vmin = 1.D-10

  ! the collector particle falls at greater speed (or should it be the greater one?)
  IF( partInfoX%v > partInfoY%v ) THEN
    partInfo1 => partInfoX
    partInfo2 => partInfoY
  ELSE
    partInfo1 => partInfoY
    partInfo2 => partInfoX
  ENDIF

  jstar  = partInfo1%j
  vi1    = partInfo1%v*(2.D0+partInfo1%j)/4.D0 ! Boehm99 below Eq (23)
  vi1    = VI1*MIN(partInfo1%alpha,1.0_wp)     ! Boehm92c Eq (6), but not in Boehm99 ! TODO remove or not?
  cdi1   = partInfo1%CDI

  rr1    = partInfo1%r
  re1    = partInfo1%r_e
  V1     = partInfo1%v
  m1     = partInfo1%m
  alf1   = partInfo1%alpha_e
  qe1    = partInfo1%q_e
  dels1  = partInfo1%DELS
  sig1   = partInfo1%SIG

  rr2    = partInfo2%r
  re2    = partInfo2%r_e
  V2     = partInfo2%v
  m2     = partInfo2%m
  alf2   = partInfo2%alpha_e
  qe2    = partInfo2%q_e
  dels2  = partInfo2%DELS
  sig2   = partInfo2%SIG

  kk2    = partInfo2%k
  
  ! reduced mass of two-body system denoted by (*), Boehm99 eq. 13
  mstar = m1*m2/(m1+m2)         
  ! corresponding radius, Boehm99 eq. 14
  rstar = rr1*rr2/(re1+re2)
  ! area of intersection, Boehm99 eq. 15
  astar = pi*rstar**2
  
  ! Assume that terminal fall speed is a random variable and follows a log-normal distribution
  ! (Appendix of Boehm 1992c)
  IF (sig1.lt.1.d-10) THEN
     phi1 = SIGN(1.D0,LOG(V1/V2))
  ELSE
     ! Boehm92c, argument in eq. 18
     arg1 = (LOG(v1/v2) + sig1**2)/(2.d0*sig1)
     ! Boehm92c, approximation of error function as given by (20)
     phi1 = TANH( 2.d0/SQRT(pi)*arg1 + .1012d0*arg1**3 )
  ENDIF
  IF (sig2.lt.1.d-10) THEN
     phi2 = SIGN(1.d0,LOG(v2/v1))
  ELSE
     ! Boehm92c eq. 18
     arg2 = (LOG(v2/v1) + sig2**2)/(2.d0*sig2)
     ! Boehm92c eq. 20
     phi2 = TANH( 2.d0/SQRT(pi)*arg2 + .1012d0*arg2**3 )
  ENDIF

  ! expected value of difference of fall speeds, Boehm92c eq. 17
  vstar = MAX(ABS( v1*phi1 + v2*phi2 ),vmin)
  
  ! boundary layer thickness, Boehm99 eq. 17
  delstar = dels1*SQRT(v1/vstar)+dels2*SQRT(v2/vstar)

  ! Boehm99 eq. 6
  alfstar = rstar*(alf1/re1+alf2/re2)
  qstar   = rstar*( qe1/re1+ qe2/re2)
  
  CALL hydrodynamicValuesBoehm( alfstar, qstar, kstar, gamstar, cdpstar )

  ! reynolds number, Boehm99 eq. 19
  nrestar = 2.D0*rstar*vstar*rhoa/nu
  
  ! Oseen drag coefficient, Boehm92a eq. 16
  cdostar = 4.5D0*kstar**2*MAX(alfstar,1.D0)

  ! Boehm99 eq. 20
  betaexp = SQRT(nrestar*cdpstar/(6.d0*kstar))

  ! Boehm92a eq. 26
  gammaexp=(cdostar-cdpstar)/(4.d0*cdpstar)
  
  ! inertial drag coefficient, Boehm92a eq. 25
  cdistar = cdpstar*(1.d0+4.d0/betaexp*(1.d0-exp(-gammaexp*betaexp)))

  ! The G of Boehm99 below eq. 23
  gg = 6.d0*kstar*nu/(rhoa*rstar*vstar*cdistar)
  
  ! The F of Boehm99 below eq. 23
  ff = SQRT(gg**2+cdi1*vi1**2/cdistar/vstar**2)
  
  ! The H Boehm99 below eq. 23
  hh = 2.d0*mstar/(rhoa*astar*cdistar*delstar)

  ! The term F/H occurs often
  ffhh=ff/hh
  
  if (ffhh.GE.10) THEN
     ! Boehm99 eq. 27
     Ecbohm_v = hh*log((ff+gg)/(2.d0*ff))+ff-gg
  ELSE
     ! Boehm99 eq. 23
     Ecbohm_v = hh*log(cosh(ffhh)+(1.d0+gg)/ff*sinh(ffhh))-gg
  ENDIF

  ! Boehm99 eq. 23
  ecbohm_v = ecbohm_v**(2.d0/jstar)

  ! Boehm94 eq. 3
  bb = 3.d0*v1/(re1+dels1)
  IF(alf2.LT.1) THEN
     ! Boehm94 eq. 8
     k2x = 0.57d0+.43d0*alf2
  ELSE
     ! Boehm94 eq. 9
     k2x = kk2**1.15d0
  ENDIF
  k2y = kk2
  
  ! Boehm94 eq. 7 
  cx = m2/(3.d0*pi*re2*nu)/k2x
  cy = m2/(3.d0*pi*re2*nu)/k2y
  
  IF (bb*cy.ge.0.5d0) THEN
     
     ! Boehm99 eq. 25 and 26
     delx = sqrt(2.d0*bb*cx/jstar+1.d0) 
     dely = sqrt(2.d0*bb*cy-1.d0)

     ! Boehm94 eq. 14
     IF(bb*cy.ge.1.d0) THEN
        ! Boehm94 eq. 14 and Boehm99 eq. 24
        tdel = cy/dely*atan(dely/(bb*cy-1.d0))
        ep = ((cosh(delx*tdel/cx)+sinh(delx*tdel/cx)/delx)*exp(-tdel/cx))**(-jstar)
     ELSE
        ! Boehm94 eq. 14 and Boehm99 eq. 30
        tdel=cy/dely*(pi-atan(dely/(1.d0-bb*cy)))
        ep=(2.d0*delx/(1.d0+delx)*exp(-(delx-1.d0)*tdel/cx))**jstar
     ENDIF
  ELSE
     ep=0.d0
  ENDIF
  ecbohm_v=ep*ecbohm_v
  
END FUNCTION Ecbohm
