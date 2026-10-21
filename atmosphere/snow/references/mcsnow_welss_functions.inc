!
! Copyright (C) 2004-2024, DWD
! See ./AUTHORS.txt for a list of authors
! See ./LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
!

SUBROUTINE m2d_particle( sp )
  TYPE(t_sp),   INTENT(inout) :: sp

  REAL(wp) :: d_solid, a_solid

  CALL m2d_particle_solid(  sp, d_solid, a_solid )
  IF ( (sp_habit == 0) .OR. ((sp_habit==1) .AND. (sp%mm>1)) ) call iceGeo_init(sp) ! If geometry not prognostic, reset V_i to meaningful value
  CALL m2d_particle_liquid( sp, d_solid, a_solid )
END SUBROUTINE m2d_particle

!
!  calculate the geometry of unrimed aggregates dependent on the monomer number 
!
PURE SUBROUTINE agg_monomer_depend( mm, mono_alf, mono_bet, mono_sig, mono_gam, unr_alf_agg, unr_bet_agg, unr_sig_agg, unr_gam_agg)
  REAL(wp), INTENT(in)  :: mm
  REAL(wp), INTENT(in)  :: mono_alf, mono_bet, mono_sig, mono_gam
  REAL(wp), INTENT(out) :: unr_alf_agg, unr_bet_agg, unr_sig_agg, unr_gam_agg
  REAL(wp), PARAMETER :: &
   alf_mod1 = -0.846, & !first parameter in rational function of a in m=aD**b
   alf_mod2 = 0.547, &  !second parameter in rational function of a in m=aD**b
   bet_mod1 = -0.146, & !first parameter in rational function of b in m=aD**b
   bet_mod2 = 0.338, &  !second parameter in rational function of b in m=aD**b
   sig_mod1 = -0.460, & !first parameter in rational function of c in A=cD**d
   sig_mod2 = 0.287, &  !second parameter in rational function of c in A=cD**d
   gam_mod1 = -0.022, & !first parameter in rational function of d in A=cD**d
   gam_mod2 = -0.165    !second parameter in rational function of d in A=cD**d

  REAL(wp) :: log10Ofmm

  log10Ofmm = log10(mm)
  
  !modify the mass
  unr_alf_agg = 10**((alf_mod1 *log10Ofmm)/(1+alf_mod2*log10Ofmm))*mono_alf
  unr_bet_agg = (bet_mod1*log10Ofmm)/(1+bet_mod2*log10Ofmm)+mono_bet
  !modify the area
  unr_sig_agg = 10**((sig_mod1 *log10Ofmm)/(1+sig_mod2*log10Ofmm))*mono_sig
  unr_gam_agg = (gam_mod1*log10Ofmm)/(1+gam_mod2*log10Ofmm)+mono_gam

END SUBROUTINE agg_monomer_depend

!
! Find max. dimension and area of a non-spherical crystal
!
PURE SUBROUTINE habit_geometry(sp,d_solid,a_solid)
  TYPE(t_sp), INTENT(in)   :: sp
  REAL(wp),   INTENT(out)  :: d_solid, a_solid

  REAL(wp) :: r_a, r_c, xi, V_tot

  IF( (sp%V_i + sp%V_r + (sp%m_f * rhoii)) < TINY(sp%m_i)) THEN
    d_solid = 0.0_wp
    a_solid = 0.0_wp
    RETURN
  ENDIF

  CALL spheroid_geometry(sp, r_a, r_c, V_tot)

  d_solid = 2._wp * MAX(r_a,r_c)
  
  ! A = xi * A_ellipsis 
  xi      = porosity_xi(sp)
  a_solid = xi * pi * r_a * max(r_a, r_c)
END SUBROUTINE habit_geometry

PURE SUBROUTINE spheroid_geometry(sp, r_a, r_c, V_tot)
  TYPE(t_sp),         INTENT(in)  :: sp
  REAL(wp),           INTENT(out) :: r_a, r_c, V_tot

  IF (sp%mm > 1._wp .AND. sp_habit > 1) THEN ! first fill rime into agg spheroid structure
    IF( (sp%V_i - sp%m_i*rhoii) > (sp%m_f*rhoii+sp%V_r) ) THEN ! rime fits inside
      V_tot = sp%V_i
    ELSE                                                       ! rime (partially) outside
      V_tot = sp%V_r + (sp%m_i+sp%m_f)*rhoii
    ENDIF
  ELSE
    V_tot = sp%V_i + sp%V_r + (sp%m_f * rhoii)
  ENDIF
  r_a   = ( z3_4pi * V_tot / sp%phi )**z1_3 
  r_c   = r_a * sp%phi
END SUBROUTINE spheroid_geometry

!
! Consideration of porosity for area following JH15 for monomers and Shima 2020 for aggregates
! TO DO: Think about porosity, when rime is included, particle could get more dense
PURE ELEMENTAL FUNCTION porosity_xi( sp ) result ( xi )
  TYPE(t_sp), INTENT(in) :: sp

  REAL(wp) :: xi, rho_app, r_a, r_c, V_tot, V_crit, Fr

  CALL spheroid_geometry(sp, r_a, r_c, V_tot)
  rho_app = (sp%m_i + sp%m_r + sp%m_f)/V_tot

  ! Aggregates need special treatment because of their low apparent density even at sphericity
  IF (sp%mm > 1._wp .AND. sp_habit > 1) THEN
    !xi = (rho_app * rhoii) ** exp(-sp%phi)        ! Shima 2020, Eq. 5 + 6 
    V_crit  = sp%V_i - sp%m_i*rhoii
    Fr = MERGE( (sp%m_f*rhoii+sp%V_r)/V_crit, 1.0_wp, v_crit>0.0_wp )      ! v_crit = 0 is possible
    IF( Fr < 1.0_wp ) THEN ! rime fits inside, interpolate
      xi = (1.0_wp-Fr) * (rho_app * rhoii)**(2.0_wp/3.0_wp) + Fr * 1.0_wp
    ELSE
      xi = 1.0_wp
    ENDIF
  ELSE IF ( (sp_habit > 0) .AND. sp%phi < 1._wp) THEN
    ! include branching: xi = (1 - phi) * (rho/rho_i) + phi (J&H 2015) ;
    ! Alternative: xi = ( (sp%m_i + sp%m_r) / (sp%V_i + sp%V_r) ) * rhoii 
    xi = (1._wp - sp%phi) * (rho_app * rhoii) + sp%phi
  ELSE
    xi = 1._wp          ! no reduction of effective area for hollowing or non-habit particles
  END IF 

END FUNCTION porosity_xi

PURE SUBROUTINE unrimed_geometry(sp,dunr,aunr)
   ! select the geometry of the unrimed particles (either pristine or aggregates)
   ! based on the aggeotype namelist parameter
   TYPE(t_sp),   INTENT(in)     :: sp
   REAL(wp),     INTENT(out)    :: dunr, aunr
   REAL(wp) :: unr_alf_agg, unr_bet_agg, unr_sig_agg, unr_gam_agg,dice,aice

   ! ice, distinguish between spherical small ice and power-law geometry, m_i = 0 is safe
   dice = (z6_pi * sp%m_i * rhoii)**z1_3
   aice = sph_sig * dice**sph_gam

  IF ( agggeo_type .eq. 3 ) THEN !   monomer number dependent fit to Jussi's plates
     call agg_monomer_depend( sp%mm, mono_alf_plates, mono_bet_plates , mono_sig_plates, mono_gam_plates, &
                            &        unr_alf_agg,     unr_bet_agg ,     unr_sig_agg,     unr_gam_agg)
     
     dunr = (sp%m_i / unr_alf_agg)**(1./unr_bet_agg)
     aunr = unr_sig_agg * dunr**unr_gam_agg
     
     IF (dunr<dice) THEN
        dunr = dice
        !aunr = aice
     END IF
  ELSE
    IF( agggeo_type .eq. 2 .AND. (sp%mm < 1.01) ) THEN !different geometry for monomers
      dunr = (sp%m_i / mono_alf_plates)**(1./mono_bet_plates)
      aunr = mono_sig_plates*dunr**mono_gam_plates
      IF (dunr<dice) THEN
        dunr = dice
        !aunr = aice
      END IF
    ELSE ! aggregate M96 in case of aggeotype<=2 and sp%mm>1
      IF( sp%m_i < mth) THEN
        ! small ice
        dunr = dice
        aunr = aice
      ELSE
        ! unrimed ice
        dunr = (sp%m_i * unr_alf_i)**unr_bet_i
        aunr = unr_sig * dunr**unr_gam 
        aunr = max(aunr,aice)
      END IF
    END IF
  END IF
END SUBROUTINE unrimed_geometry

!
! calculate d and p_area of the solid (ice+rime) part of the particle
!
PURE SUBROUTINE m2d_particle_solid( sp, d_solid, a_solid)
  TYPE(t_sp),   INTENT(in)        :: sp
  REAL(wp),     INTENT(out)       :: d_solid, a_solid

  IF ( ( sp_habit == 1 .AND. sp%mm < 1.01 ) .OR. (sp_habit == 2) ) THEN
    CALL habit_geometry( sp, d_solid, a_solid )
  ELSE 
    IF ( snowgeo_type .eq. 3 ) THEN
      CALL m2d_similarity( sp, d_solid, a_solid )
    ELSE
      CALL m2d_fillin(     sp, d_solid, a_solid )  ! snowgeo_type.eq.2.or.snowgeo_type.eq.22
    ENDIF
  END IF
END SUBROUTINE m2d_particle_solid


!
! compute diameter from m_i, m_r, v_r using the fill-in model
! compute proj. area using linear interpolation as in
!   based on Morrison and Grabowski 'A Novel Approach...' 2007
! 
PURE SUBROUTINE m2d_fillin( sp, d_solid, a_solid )
  TYPE(t_sp),   INTENT(in)    :: sp
  REAL(wp),     INTENT(out)   :: d_solid, a_solid

  REAL(wp) :: dunr, aunr
  REAL(wp) :: Fr, v_crit

  
  !get geometry of unrimed particles
  CALL unrimed_geometry(sp,dunr,aunr)

  ! rime, distinguish between fill-in and spherical graupel, m_i = 0 and/or m_r = 0 is safe
  v_crit  = pi_6*dunr**3 - sp%m_i*rhoii
  Fr = MERGE( (sp%m_f*rhoii+sp%v_r)/v_crit, 1.0_wp, v_crit>0.0_wp )      ! v_crit = 0 for small ice or no ice
  IF( Fr < 1.0_wp ) THEN ! rime inside (or no rime at all)
    IF(snowgeo_type.eq.22) Fr = sqrt(Fr)                  ! sqrt interpolation of cross sectional area
    d_solid = dunr
    a_solid = (1.0_wp-Fr) * aunr + Fr * gra_sig*dunr**gra_gam
  ELSE                   ! rime outside => spherical rime particle, i.e. graupel
    d_solid = (z6_pi*(sp%v_r + (sp%m_i+sp%m_f)*rhoii))**z1_3
    a_solid = gra_sig*d_solid**gra_gam
  END IF
END SUBROUTINE m2d_fillin

SUBROUTINE m2d_particle_liquid(sp, d_solid, a_solid, smoothFraction)
  TYPE(t_sp), INTENT(inout)         :: sp
  REAL(wp),   INTENT(in)            :: d_solid, a_solid
  REAL(wp),   OPTIONAL, INTENT(out) :: smoothFraction

  REAL(wp), PARAMETER :: wr = 33 ! Seifert, A. Geosci. Model Dev., 7, 463-478, 2014
  REAL(wp) :: Fr, v_crit, v_w_out, v_w_torus, d_torus, v_total, d_volume, d_max_w, part, smoothPart
  REAL(wp), PARAMETER :: lambda = 4.7_wp/1000.0_wp

  sp%statusb = IBCLR( sp%statusb, sp_stat%liquidOut)
  sp%statusb = IBCLR( sp%statusb, sp_stat%liquidFillin)
  sp%statusb = IBCLR( sp%statusb, sp_stat%liquidSurf)

  ! water, distinguish between water at surface and water inside rime/no water, m_w=0 is safe, m_r=m_i=0 is safe
  v_w_out = sp%m_r*rhoii + sp%m_w*rholi - sp%v_r
  IF( v_w_out > 0.0_wp ) THEN ! water outside rime
    sp%statusb = IBSET( sp%statusb, sp_stat%liquidOut )

    v_total    = (sp%m_i+sp%m_f)*rhoii + sp%v_r + v_w_out
    d_volume   = (v_total*z6_pi)**z1_3          
    d_max_w    = d_volume*EXP(wr*d_volume)                 ! max drop diameter, Seifert et al. GeoSci. Model Dev. 7,463,2014 eq.3
    ! distinguish between solid core and pure water drop
    IF( d_solid > 0.0_wp ) THEN   ! solid core
      v_crit  = pi_6*d_solid**3 - (sp%m_i+sp%m_f)*rhoii - sp%v_r
      Fr = MERGE( v_w_out/v_crit, 1.0_wp, v_crit>0.0_wp )  ! v_crit = 0 for small ice or graupel 
      ! distinguish between spherical fill-in and water torus
      IF( Fr < 1.0_wp ) THEN          ! water fits inside solid maximum dimension sphere
        sp%statusb = IBSET( sp%statusb, sp_stat%liquidFillin )
        IF(snowgeo_type.eq.22) Fr = sqrt(Fr)               ! sqrt interpolation of cross sectional area
          sp%d       = d_solid
          sp%p_area  = (1.0_wp-Fr) * a_solid + Fr * pi_4*d_solid**2
          smoothPart = Fr
      ELSE                            ! build up water torus
        v_w_torus  = v_w_out - v_crit                       ! Fr >=1.0 => v_w_torus >= 0.0
        d_torus    = SQRT(8.0_wp*v_w_torus/(pi**2*d_solid)) ! assume water to form half-torus with inner radius equal to d_solid/2
        IF( d_torus < d_solid ) THEN    ! max diameter is the sum of the solid diameter plus the torus diameter
          sp%d       = (d_solid + d_torus)
          ! if v_crit>0 and Fr > 1 => smooth at whole surface
          ! if v_crit=0, i.e. small Ice or graupel, smoothPart = A_halftorus/(A_halftorus - 2*A_sphericalcap)
          !smoothPart = MERGE( 1.0_wp, pi_2*d_torus/(d_solid+d_torus*(pi_2-1.0_wp)), v_crit>0.0_wp )
          smoothPart = MERGE( pi_2*d_torus/(d_solid+d_torus*(pi_2-1.0_wp)), 1.0_wp, v_crit < EPSILON(v_crit) )
        ELSE                            ! interpolate between sphere+torus and max dimension of pure water drop
          part       = (2.0_wp*d_solid/(d_solid+d_torus))**3
          sp%d       = (d_solid + d_torus)*part + (1.0_wp-part)*d_max_w
          sp%statusb = IBSET( sp%statusb, sp_stat%liquidSurf )
          smoothPart = 1.0_wp
        ENDIF  
        sp%p_area  = pi_4*sp%d**2
      END IF
    ELSE                          ! no ice/rime, pure water drop
          sp%d       = d_max_w
          sp%p_area  = pi_4*sp%d**2
          ! pure water, flatten rain drop
          sp%phi     = exp(-sp%d/lambda) + (1.0_wp - exp(-sp%d/lambda)) / (1.0_wp+sp%d/lambda)  
          sp%statusb = IBSET( sp%statusb, sp_stat%liquidSurf )
          smoothPart = 1.0_wp
    ENDIF
  ELSE                        ! water inside rime ( or no water at all )
          sp%d       = d_solid
          sp%p_area  = a_solid
          smoothPart = 0.0_wp
  ENDIF
  IF( PRESENT(smoothFraction) ) smoothFraction = smoothPart
END SUBROUTINE m2d_particle_liquid

PURE SUBROUTINE m2d_similarity( sp, d_solid, a_solid )
  TYPE(t_sp),   INTENT(in)    :: sp
  REAL(wp),     INTENT(out)   :: d_solid, a_solid

  REAL(wp)  :: Fr, vtan, s, kappa, f, df, ds, &
       &       dunr, aunr, dgrp, agrp, coeff, v_r

  IF( (sp%m_i + sp%m_f + sp%m_r) < TINY(sp%m_i) ) THEN
    d_solid = 0.0_wp
    a_solid = 0.0_wp
    RETURN
  ENDIF

  Fr  = sp%m_r / (sp%m_i + sp%m_f + sp%m_r)
  v_r = sp%v_r + sp%m_f*rhoii               ! zero order insertation of m_f ! TODO howto similarity with m_f?

  !get geometry of unrimed particles
  CALL unrimed_geometry(sp,dunr,aunr)

  dgrp = (z6_pi*v_r)**z1_3
  agrp = pi/4.*dgrp**2

  d_solid  = dunr
  s  = 1.0
  ds = 1.0
  IF (Fr < 1e-11_wp) THEN
     ! unrimed ice
     a_solid = aunr     
  ELSEIF (Fr < 0.9999_wp) THEN
     kappa = z6_pi*v_r/d_solid**3     
     IF (kappa < 1) THEN
        s = 0.8*(ATAN(kappa)*a_m_i)**p_m_i
     ELSE
        s = kappa**z1_3+1._wp
     END IF     
     DO WHILE (ds/s > 1e-4)
        vtan = TANH( a_m*s**p_m )
        f  = vtan*(s+1._wp)**3 - kappa
        df = (s+1._wp)**2*( (1._wp-vtan*vtan)*a_m*p_m * s**(p_m-1)*(s+1) + 3.*vtan)
        ds = f/df
        s  = s-ds
     END DO     
     d_solid = d_solid*(s+1._wp)
     kappa = (6.*v_r/pi)/d_solid**3
     coeff = 0.32*tanh(7.36*kappa)+0.52*kappa   ! rho_rime = 270/700, dendrites large
     agrp  = pi/4.*d_solid**2
     a_solid  = coeff * (agrp-aunr) + aunr
  ELSE
     ! graupel
     d_solid  = dgrp
     a_solid  = agrp
  END IF
END SUBROUTINE m2d_similarity

SUBROUTINE iceGeo_init(sp)
  TYPE(t_sp),   INTENT(inout)     :: sp

  REAL(wp)            :: aunr, dunr, V_i
    
  CALL unrimed_geometry(sp,dunr,aunr)
  IF( (agggeo_type > 1) .AND. (sp%mm < 1.01) ) THEN
    !sp%phi = MIN(mono_zeta_plates * dunr**mono_epsilon_plates, 1.0_wp)
    !sp%V_i = 3.0_wp*sqrt(3.0_wp)/8.0_wp *dunr**3 * sp%phi ! volume of hexagonal prism
    sp%V_i = sp%m_i * rhoii ! solid ice plate
    sp%phi = MIN( (sp%V_i + sp%V_r + (sp%m_f * rhoii)) / (pi_6 * dunr**3), 1.0_wp ) ! fill-in of rime
    ! for only ice this corresponds to (alpha/rhoi)/(pi/6) d**(beta-3), which seems close enough
    ! 0.00125058*d**-0.526 vs 0.0015488772*d**-0.53
  ELSE
    sp%V_i = pi_6 * dunr**3
    sp%phi = 1.0_wp
  ENDIF
END SUBROUTINE iceGeo_init

!
! Copyright (C) 2004-2024, DWD
! See ./AUTHORS.txt for a list of authors
! See ./LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
!

ELEMENTAL FUNCTION dynamic_viscosity( T ) result(eta)
    REAL(wp), INTENT(in) :: T
    REAL(wp)             :: eta
    ! REAL(wp) :: Tc
    REAL(wp), PARAMETER  :: A_S = 1.458e-6, T_S = 110.4
    ! Tc = T - T_3
    ! Where is this from?
    ! eta = eta0 * (1._wp + Tc*(0.00285215366705471478_wp - Tc*6.9848661233993015133e-6_wp))  

    ! Sutherland law
    eta = A_S * T**(3./2.) / ( T + T_S )
END FUNCTION dynamic_viscosity

!
! saturation vapor pressure over ice 
!
PURE ELEMENTAL FUNCTION psat_ice( T ) result( psati )
    REAL(wp), INTENT(in) :: T
    REAL(wp)             :: psati
    ! Beheng, lecture scrpt
    ! psati = 610.64_wp * EXP( 21.88_wp * Tc / (Tc + 265.5_wp) )
    psati = e_3 * EXP( A_i * (T - T_3) / (T - B_i) )
    ! alternative Goff-Gratch equation
    ! tt0 = 273.16_wp / atmo%T
    ! psati = 611.73_wp * 10**(-9.09718_wp*(tt0-1._wp) - 3.56654_wp*LOG10(tt0) + &
    !        0.876793_wp*(1._wp-1._wp/tt0))
END FUNCTION psat_ice

!
! saturation vapor pressure over water
!
PURE ELEMENTAL FUNCTION psat_water( T ) result ( psatw )
    REAL(wp), INTENT(in) :: T
    REAL(wp)             :: psatw
    ! Beheng, lecture scrpt
    ! psatw = 610.7_wp * EXP( 17.15_wp * Tc / (Tc + 234.9_wp) )
    psatw = e_3 * EXP( A_w * (T - T_3) / (T - B_w) )
    ! alternative Goff-Gratch equation
    ! tts = 373.15_wp / atmo%T
    ! psatw = 101325_wp * 10**(-7.90298_wp*(tts-1._wp) + 5.02808_wp*LOG10(tts) - &
    !  1.3816e-7_wp*(10**(11.344*(1._wp-1._wp/tts))-1._wp) + &
    !  8.1328e-3_wp*(10**(3.49149_wp*(1._wp-tts)-1._wp)))
END FUNCTION psat_water

!
! [J/msK] thermal conductivity of dry air (PK '97, p. 508)
!
PURE ELEMENTAL FUNCTION thermalCond_dry( T ) result ( Ka )
    REAL(wp), INTENT(in) :: T
    REAL(wp)             :: Ka

    Ka = 418.4e-5_wp*( 5.69_wp + 0.017_wp*(T - T_3) )
END FUNCTION thermalCond_dry

!
! [m2/s] diffusivity of water vapor (PK '97, p. 503)
! NOTE: Ghosh et al. (2007) proposed a more accurate expression
!
PURE ELEMENTAL FUNCTION diffusivity_vapour( T, p ) result( Dv )
    REAL(wp), INTENT(in) :: T, p
    REAL(wp)             :: Dv

    Dv = 0.211e-4_wp * (T * T_3i)**1.94_wp * ps / p
END FUNCTION diffusivity_vapour

!
! [J/kg] latent heat release of deposition/sublimation
!
PURE ELEMENTAL FUNCTION latentHeat_depo( T ) result ( Ld )
    REAL(wp), INTENT(in) :: T
    REAL(wp)             :: Ld
    REAL(wp) :: Tc
    Tc = T - T_3
    Ld = 2834100_wp - Tc*(290_wp + Tc*4_wp) 
END FUNCTION latentHeat_depo

!
! [J/kg] latent heat release of evaporation/condensation
!
PURE ELEMENTAL FUNCTION latentHeat_evap( T ) result ( Le )
    REAL(wp), INTENT(in) :: T
    REAL(wp)             :: Le
    REAL(wp) :: Tc
    Tc = T - T_3
    Le = 2500800_wp + Tc*(-2360_wp + Tc*(1.6_wp - 0.06_wp*Tc)) ! Slavko Wikipedia, more efficient, no power involved
    !Le = 2.5008e6_wp * (T_3/T)**(0.167_wp + 3.67e-4_wp* T) ! RH87a
END FUNCTION latentHeat_evap

!
! Copyright (C) 2004-2024, DWD
! See ./AUTHORS.txt for a list of authors
! See ./LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
!

SUBROUTINE vapourDiffusion( sp, vt, dtg_f, atmo, dMdt, dQdt, latentHeatCoeff, V_new, phi_tend, factor )
  TYPE(t_sp),    INTENT(in)    :: sp
  REAL(wp),      INTENT(in)    :: vt
  REAL(wp),      INTENT(in)    :: dtg_f
  TYPE(t_atmo),  INTENT(in)    :: atmo
  REAL(wp),      INTENT(out)   :: dMdt
  REAL(wp),      INTENT(out)   :: dQdt
  REAL(wp),      INTENT(out)   :: latentHeatCoeff
  REAL(wp),      INTENT(out)   :: V_new
  REAL(wp),      INTENT(out)   :: phi_tend
  REAL(wp), OPTIONAL, INTENT(out) :: factor

  REAL(wp) :: es, L, Dv, Ka, C, sqrtRe, cubrPr, cubrSc, Hi, Re

  !CALL m2d_particle( sp ) ! update particle dimensions before velocity and capacitance computation (already done in sp_depo)
  !vt = vterm_snow( atmo, sp )  ! modify heat and vapour transfer by ventilation
  Re = atmo%rho * vt * sp%d / atmo%eta              ! we can use d_max since we re-scale in ventilation for habit particles
  sqrtRe = SQRT( Re )
  Dv = diffusivity_vapour( atmo%T, atmo%p )
  Ka = thermalCond_dry( atmo%T )
  cubrPr = ( atmo%eta*c_pa/Ka )**z1_3
  cubrSc = ( atmo%eta/(atmo%rho*Dv) )**z1_3
  C = capacitance( sp, Re )

  IF( (sp%m_i + sp%m_f + sp%m_r) < TINY(sp%m_i ) ) THEN ! only water -> condensation/evaporation
    Dv = Dv*ventilation_coef_w( sqrtRe*cubrSc )
    Ka = Ka*ventilation_coef_w( sqrtRe*cubrPr )

    es = atmo%psatw
    L  = latentHeat_evap( atmo%T )
  ELSE                                         ! no water outside -> deposition/sublimation
    Dv = Dv*ventilation_coef_i( sqrtRe*cubrSc, sp%phi )  ! CHRS for old bugged version (reproduce paper) comment these two lines
    Ka = Ka*ventilation_coef_i( sqrtRe*cubrPr, sp%phi )

    es = atmo%psati
    L  = latentHeat_depo( atmo%T )
  ENDIF
  Hi = (L/(Rv*atmo%T) - 1.0_wp)*(L*Dv*es)/(Ka*Rv*atmo%T**2)

  ! neglecting curvature effects relevant for small particles
  dMdt =  4.0_wp * pi * C *   Dv * (atmo%pv - es)                         / (Rv * atmo%T * (1._wp+Hi))
  ! sum heat for possible melting at surface temperature T_3 ( psat_ice(T_3) = psat_water(T_3), and assumed L ~ const )
  dQdt = -4.0_wp * pi * C * ( Dv * (atmo%pv - (atmo%T/T_3)*psat_ice(T_3)) / (Rv * atmo%T) * L + Ka * (atmo%T - T_3) ) ! sum heat for possible melting
  latentHeatCoeff = L
  IF( PRESENT(factor) ) factor = L/(4.0_wp*pi*C*Ka)

  ! habit prediction only when ice exists
  IF ( dMdt > 0._wp ) THEN
    IF ( ( sp_habit == 1 .AND. sp%mm < 1.01 ) .OR. (sp_habit == 2) ) THEN
      CALL habit_predict(sp, dMdt, dtg_f, atmo%T, Dv, vt, sqrtRe*cubrSc, V_new, phi_tend)
    END IF
  END IF
END SUBROUTINE vapourDiffusion

PURE ELEMENTAL FUNCTION ventilation_coef_i( X, phi ) result ( fv )
  REAL(wp), INTENT(in) :: X
  REAL(wp), INTENT(in) :: phi
  REAL(wp)             :: fv, X_equiv
  REAL(wp), PARAMETER  :: c1 = 0.0280_wp
  REAL(wp), PARAMETER  :: c2 = 0.0028_wp

  IF (phi >= 1._wp) THEN
    X_equiv = X * phi**(-1./3.)
  ELSE IF (phi < 1._wp) THEN
    X_equiv = X * phi**(1./6.)
  END IF

! spheres:  fv = 1.00 + 0.14*X^2,         for X < 1   
!           fv = 0.86 + 0.28*X,           for X > 1  
  IF( X_equiv < 1.0_wp ) THEN
    fv = 1.0_wp  + 0.14_wp*X_equiv**2
  ELSE
    fv = 0.86_wp + 0.28_wp*X_equiv
  END IF

! prolates: fv = fv + c1 * phi * X,     for phi > 1
! oblates:  fv = fv + c2 / phi * X^3/2, for phi < 1 
  IF (phi > 1._wp) THEN
  ! sqrt(Re_equiv) = sqrt(Re_c / phi * phi^(1/3) )= sqrt(Re_c * phi^(-2/3)) = sqrt(Re_c) * phi^(-1/3) 
    fv      = fv + c1 * phi * X_equiv
  ELSE IF (phi < 1._wp) THEN
  ! sqrt(Re_Dequiv) = sqrt(Re_a * phi^(1/3)) = X_a * phi^(1/6)
    fv      = fv + c2 / phi * X_equiv**(3./2.)
  END IF
END FUNCTION ventilation_coef_i


!
! ventilation coefficient of water 
! could be extended by a shape-dependent term but effect should be small because oblate flattening
!
PURE ELEMENTAL FUNCTION ventilation_coef_w( X ) result ( fv )
  REAL(wp), INTENT(in) :: X
  REAL(wp)             :: fv

  IF( X < 1.4_wp ) THEN
    fv = 1.0_wp  + 0.108*X**2
  ELSE
    fv = 0.78_wp + 0.308*X
  END IF
END FUNCTION ventilation_coef_w

PURE FUNCTION capacitance( sp, Re ) result ( C )
  TYPE(t_sp),   INTENT(in)    :: sp
  REAL(wp),     INTENT(in)    :: Re
  REAL(wp)                    :: C

  REAL(wp), PARAMETER :: ReLimit(2) = (/ 3000.0_wp, 6000.0_wp /)

  REAL(wp) :: v_crit, Fr, C_unr, d, d_solid, a_solid, part, eps, r_a, r_c, V_tot
  REAL(wp), PARAMETER ::  C_sph = 0.5_wp ! graupel (spherical ice)

  ! capacitance of explicitly evolving particles/spheroids, following Chen&Lamb '94 (p.9)
  ! TODO: there is no dependency of capacitance on porosity. This might lead to an overestimation
  ! TODO: of capacitance (cf. Westbrook et al. 2008).
  IF ( ( (sp_habit == 1 .AND. sp%mm < 1.01_wp) .OR. (sp_habit ==2 ) )                             &
    .AND. ( sp%m_i + sp%m_r + sp%m_f > 0.0_wp) ) THEN
    CALL spheroid_geometry(sp, r_a, r_c, V_tot)

    IF (sp%phi < 1._wp) THEN                ! Chen&Lamb '94 Eq. 39
      eps = SQRT(1._wp - sp%phi**2)
      C   = r_a * eps / ASIN(eps)
    ELSE IF (sp%phi > 1._wp) THEN           ! Chen&Lamb '94 Eq. 40
      eps = SQRT( 1._wp - sp%phi**(-2) )
      C   = r_c * eps / LOG( (1._wp + eps) * sp%phi )
    ELSE
      C   = (z3_4pi * V_tot)**z1_3          ! spherical particle C = r
    END IF
  ELSE
    C_unr = 0.25_wp                        ! aggregates ( eq. A9 of Soelch and Kaercher 2010) after Westbrook 2008
    IF(sp%mm < 1.01_wp) C_unr = 1.0_wp / pi     ! spherical plate PK 13-77
    v_crit  = pi_6 * sp%d**3 - sp%m_i * rhoii                                                 ! v_crit based on sp%d, hence Fr <=1
    Fr = MERGE( (sp%m_f*rhoii + MAX(sp%v_r, sp%m_r*rhoii + sp%m_w * rholi)) / v_crit,              &
                1.0_wp, v_crit > 0.0_wp )  ! v_crit = 0 for small ice 
    d  = sp%d
    IF( BTEST(sp%statusb, sp_stat%liquidOut) .AND. (sp%m_i+sp%m_f+sp%m_r > 0.0_wp)                 &
        &.AND. (.NOT. BTEST(sp%statusb, sp_stat%liquidFillin)) ) THEN ! 8 Torus, 9 Torus to rain
      part = MIN( 1.0_wp, (Re - ReLimit(1)) / (ReLimit(2) - ReLimit(1)) )
      IF( part > 0.0_wp ) THEN
        CALL m2d_particle_solid( sp, d_solid, a_solid )
        d = (1.0_wp - part) * sp%d + part * d_solid
      END IF
    ENDIF
    C  = d * ( (1.0_wp - Fr) * C_unr + Fr * C_sph )
  END IF
END FUNCTION capacitance
