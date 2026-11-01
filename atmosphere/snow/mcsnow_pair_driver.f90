! Prescribed pure-ice pairs; original geometry, drag and collection efficiency.
! The outer ckt0 encounter area follows mo_coll_utils; no growth or random speed.
program mcsnow_pair_driver
  use mo_constants, only: wp, mp, pi, i8
  use mo_sp_types, only: t_sp
  use mo_atmo_types, only: t_atmo, Rd, Rv, rhoi
  use mo_atmo, only: psat_water, psat_ice, dynamic_viscosity
  use mo_sp_nml, only: sp_habit, sp_vterm_id, snowgeo_type, agggeo_type, prolate_intp, rm, Vm, sp_stick_effi
  use mo_mass2diam, only: m2d_particle
  use mo_velocity, only: vterm_snow
  use mo_colleffi, only: Ecbohm_spheroid, stick_effi_tempOnly
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  type(t_sp) :: left, right
  type(t_atmo) :: air
  real(wp) :: t, pressure, m1, v1, phi1, m2, v2, phi2, speed1, speed2, radius1, radius2
  real(wp) :: efficiency, reverse_efficiency, stick, area, speed, kernel, values(22)
  integer :: unit, ios, case_id, habit, sticking
  character(len=1024) :: input_path

  call get_command_argument(1,input_path)
  if (len_trim(input_path)==0) error stop 'Input file required'
  open(newunit=unit,file=trim(input_path),status='old',action='read')
  sp_vterm_id=3; snowgeo_type=2; agggeo_type=1; prolate_intp=.true.
  rm=1e-5_wp; Vm=0._wp
  write(*,'(a)') 'case_id,temperature_K,pressure_Pa,mass_left_kg,mass_right_kg,'// &
    'volume_left_m3,volume_right_m3,phi_left,phi_right,diameter_left_m,diameter_right_m,'// &
    'area_left_m2,area_right_m2,speed_left_m_s,speed_right_m_s,encounter_area_m2,'// &
    'relative_speed_m_s,sticking,collection_efficiency,reverse_collection_efficiency,kernel_m3_s,'// &
    'air_density_kg_m3,dynamic_viscosity_Pa_s'
  do
    read(unit,*,iostat=ios) case_id,t,pressure,habit,sticking,m1,v1,phi1,m2,v2,phi2
    if (ios<0) exit
    if (ios/=0) error stop 'Malformed pair'
    if (.not.all(ieee_is_finite([t,pressure,m1,v1,phi1,m2,v2,phi2]))) error stop 'Nonfinite input'
    if (min(t,pressure,m1,v1,phi1,m2,v2,phi2)<=0) error stop 'Nonpositive input'
    if (habit<0.or.habit>1.or.sticking<0.or.sticking>1) error stop 'Unsupported switch'
    if (habit==1.and.(m1>rhoi*v1*(1+1e-12_wp).or.m2>rhoi*v2*(1+1e-12_wp))) error stop 'Invalid density'
    sp_habit=habit; sp_stick_effi=sticking
    air%T=t; air%p=pressure; air%pv=psat_water(t)
    air%psatw=air%pv; air%psati=psat_ice(t)
    air%rh=100._wp; air%ssat=air%pv/air%psati-1._wp
    air%rho=(pressure-air%pv)/(Rd*t)+air%pv/(Rv*t)
    air%eta=dynamic_viscosity(t); air%etai=1._wp/air%eta
    air%qv=(Rd/Rv)*air%pv/(pressure-(1._wp-Rd/Rv)*air%pv); air%qc=0._wp
    call particle(left,m1,v1,phi1)
    call particle(right,m2,v2,phi2)
    speed1=vterm_snow(air,left); speed2=vterm_snow(air,right)
    radius1=.5_wp*left%d*min(1._wp,1._wp/sqrt(left%phi))
    radius2=.5_wp*right%d*min(1._wp,1._wp/sqrt(right%phi))
    efficiency=Ecbohm_spheroid(air,left,right)
    reverse_efficiency=Ecbohm_spheroid(air,right,left)
    stick=stick_effi_tempOnly(t)
    area=pi*(radius1+radius2)**2
    speed=abs(speed1-speed2)
    kernel=efficiency*area*speed
    values=[t,pressure,m1,m2,v1,v2,real(left%phi,wp),real(right%phi,wp), &
      real(left%d,wp),real(right%d,wp),real(left%p_area,wp),real(right%p_area,wp), &
      speed1,speed2,area,speed,stick,efficiency,reverse_efficiency,kernel,air%rho,air%eta]
    if (.not.all(ieee_is_finite(values))) error stop 'Nonfinite pair output'
    write(*,'(i0,a,*(es25.17e3,:,","))') case_id,',',values
  end do
  close(unit)
contains
  subroutine particle(sp,mass,volume,phi)
    type(t_sp), intent(out) :: sp
    real(wp), intent(in) :: mass, volume, phi
    sp%statusb=0; sp%xi=1_i8; sp%mm=1._mp
    sp%m_w=0._wp; sp%m_r=0._wp; sp%m_f=0._wp; sp%V_r=0._wp
    sp%m_i=mass; sp%V_i=volume; sp%phi=phi
    sp%T=t; sp%init_temp=t; sp%init_ssat=air%ssat
    sp%z=0._wp; sp%cellId=1; sp%start_pe=0; sp%local_idx=1
    call m2d_particle(sp)
  end subroutine particle
end program mcsnow_pair_driver
