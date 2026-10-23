! Diagnostic caller: original McSnow geometry and drag, with no growth steps.
! Controls alter a copy after geometry; they are not adopted particle states.
program mcsnow_drag_driver
  use mo_constants, only: wp, mp, pi
  use mo_sp_types, only: t_sp
  use mo_atmo_types, only: t_atmo, Rd, Rv, rhoi, grav
  use mo_atmo, only: psat_water, dynamic_viscosity
  use mo_sp_nml, only: sp_habit, sp_vterm_id, snowgeo_type, agggeo_type, prolate_intp, rm, Vm
  use mo_mass2diam, only: m2d_particle, spheroid_geometry, porosity_xi
  use mo_velocity, only: hydrodynamicValuesBoehm, Reynolds_Boehm, boehmVel, vterm_snow
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  type(t_sp) :: original, particle
  type(t_atmo) :: air
  real(wp) :: temperature, pressure, mass, volume, phi, r_a, r_c, v_total, xi, q, k, gamma, cd, best, re
  real(wp) :: via_re, via_boehm, via_selector, values(26)
  ! The rough-particle value in the pinned boehmVel, also used by the separate Reynolds call.
  real(wp), parameter :: rough_best_scale=2.8e6_wp
  integer :: unit, ios, case_id, control
  character(len=1024) :: input_path

  call get_command_argument(1,input_path)
  if (len_trim(input_path)==0) error stop 'Input file required'
  open(newunit=unit,file=trim(input_path),status='old',action='read')
  sp_habit=1; sp_vterm_id=3; snowgeo_type=2; agggeo_type=1; prolate_intp=.true.
  rm=1e-5_wp; Vm=0._wp
  write(*,'(a)') 'case_id,control_id,temperature_K,pressure_Pa,vapor_pressure_Pa,air_density_kg_m3,'// &
    'dynamic_viscosity_Pa_s,ice_density_kg_m3,gravity_m_s2,input_mass_kg,ice_volume_m3,aspect_ratio,'// &
    'equatorial_dimension_m,polar_dimension_m,diameter_m,native_projected_area_m2,native_porosity,'// &
    'drag_mass_kg,projected_area_m2,area_ratio,k,gamma_big,pressure_drag_coefficient,best_number,'// &
    'reynolds_number,speed_from_reynolds_m_s,speed_boehm_m_s,terminal_speed_m_s'
  do
    read(unit,*,iostat=ios) case_id,temperature,pressure,mass,volume,phi
    if (ios<0) exit
    if (ios/=0) error stop 'Malformed input row'
    if (.not.all(ieee_is_finite([temperature,pressure,mass,volume,phi]))) error stop 'Nonfinite input'
    if (min(temperature,pressure,mass,volume,phi)<=0.or.phi>1) error stop 'Requires positive oblate state'
    air%T=temperature; air%p=pressure; air%pv=psat_water(temperature)
    air%rho=(pressure-air%pv)/(Rd*temperature)+air%pv/(Rv*temperature)
    air%eta=dynamic_viscosity(temperature)
    original%statusb=0; original%mm=1._mp
    original%m_w=0._wp; original%m_r=0._wp; original%m_f=0._wp; original%V_r=0._wp
    original%m_i=mass; original%V_i=volume; original%phi=phi
    call m2d_particle(original)
    call spheroid_geometry(original,r_a,r_c,v_total)
    xi=porosity_xi(original)
    do control=0,3
      particle=original
      ! 0: archived inputs. 1: area only. 2: mass only. 3: solid-envelope mass and area.
      if (control==1.or.control==3) particle%p_area=pi/4._wp*particle%d**2
      if (control==2.or.control==3) particle%m_i=rhoi*volume
      q=particle%p_area/(pi/4._wp*particle%d**2)
      call hydrodynamicValuesBoehm(phi,q,k,gamma,cd)
      re=Reynolds_Boehm(air,phi,particle%m_i,q,rough_best_scale,k,gamma,best)
      via_re=re*air%eta/(particle%d*air%rho)
      via_boehm=boehmVel(air,particle)
      via_selector=vterm_snow(air,particle)
      values=[temperature,pressure,air%pv,air%rho,air%eta,rhoi,grav,mass,volume,phi, &
        2*r_a,2*r_c,particle%d,original%p_area,xi,particle%m_i,particle%p_area,q,k,gamma,cd,best,re, &
        via_re,via_boehm,via_selector]
      if (.not.all(ieee_is_finite(values))) error stop 'Nonfinite trace'
      write(*,'(i0,a,i0,a,*(es25.17e3,:,","))') case_id,',',control,',',values
    end do
  end do
  close(unit)
end program mcsnow_drag_driver
