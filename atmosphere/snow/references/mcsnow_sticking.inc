!
! Copyright (C) 2004-2024, DWD
! See ./AUTHORS.txt for a list of authors
! See ./LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
!

PURE FUNCTION stick_effi_tempOnly( T ) RESULT( se )
  USE mo_sp_nml, ONLY: istick=>sp_stick_effi
  REAL(wp), INTENT(IN) :: T
  REAL(wp)             :: se
  REAL(wp)             :: Tc
  !INTEGER, PARAMETER   :: istick = 3

  Tc = T - T_3
  IF     (istick.eq.0) THEN
        se = 1.0_wp
  ELSEIF (istick.eq.1) THEN
     ! after Pruppacher and Klett, p. 691, originally Mitchell 1988
     IF( Tc >= 0.0_wp ) THEN
        se = 1.0_wp
     ELSE IF ( Tc >= -4._wp ) THEN
        se = 0.1_wp
     ELSE IF ( Tc >= -6._wp ) THEN
        se = 0.6_wp
     ELSE IF ( Tc >= -9._wp ) THEN
        se = 0.1_wp
     ELSE IF ( Tc >= -12.5_wp ) THEN
        se = 0.4_wp
     ELSE IF ( Tc >= -17_wp ) THEN
        se = 1.0_wp
     ELSE IF ( Tc >= -20_wp ) THEN
        se = 0.4_wp
     ELSE
        se = 0.25_wp
     END IF
  ELSEIF (istick.eq.2) THEN
     ! piecewise linear sticking efficiency with maximum at -15 C,
     ! inspired by Figure 14 of Connolly et al. ACP 2012, doi:10.5194/acp-12-2055-2012
     ! Value at -40 C is based on Kajikawa and Heymsfield as cited by Philips et al. (2015, JAS)
     IF ( Tc >= 0_wp ) THEN
        se = 0.14
     ELSEIF ( Tc >= -10_wp ) THEN
        se = -0.01_wp*(Tc+10_wp)+0.24_wp
     ELSEIF ( Tc >= -15_wp ) THEN
        se = -0.08_wp*(Tc+15_wp)+0.64_wp
     ELSEIF ( Tc >= -20_wp ) THEN
        se =  0.10_wp*(Tc+20_wp)+0.14_wp
     ELSEIF ( Tc >= -40_wp ) THEN
        se = 0.005_wp*(Tc+40_wp)+0.04_wp
     ELSE
        se =  0.04_wp
     END IF
  ELSEIF (istick.eq.3) THEN
      ! as option 1, but with factor 0.5, i.e., the lower range of Figure 14
     IF ( Tc >= 0_wp ) THEN
        se = 0.14
     ELSEIF ( Tc >= -10_wp ) THEN
        se = -0.01_wp*(Tc+10_wp)+0.24_wp
     ELSEIF ( Tc >= -15_wp ) THEN
        se = -0.08_wp*(Tc+15_wp)+0.64_wp
     ELSEIF ( Tc >= -20_wp ) THEN
        se =  0.10_wp*(Tc+20_wp)+0.14_wp
     ELSEIF ( Tc >= -40_wp ) THEN
        se = 0.005_wp*(Tc+40_wp)+0.04_wp
     ELSE
        se =  0.04_wp
     END IF
     se = 0.5 * se
  ELSEIF (istick.eq.4) THEN
     ! inspired by PK (option 1) and Connolly (option 2)
     IF ( Tc >= -1 ) THEN
       se = 1.0
     ELSEIF (  Tc >= -6 ) THEN
       se = 0.19*(Tc+6)+0.24
     ELSEIF (  Tc >= -10 ) THEN
       se = 0.24
     ELSEIF (  Tc >= -12.5 ) THEN
       se = -0.304*(Tc+12.5)+1.0
     ELSEIF (  Tc >= -17 ) THEN
       se = 1.0
     ELSEIF (  Tc >= -20 ) THEN
       se =  0.286666*(Tc+20)+0.14
     ELSEIF (  Tc >= -40 ) THEN
       se = 0.005*(Tc+40)+0.04
     ELSE
       se =  0.04
     END IF
  END IF

END FUNCTION stick_effi_tempOnly
