! Fixed-state caller. Original geometry, air and velocity modules are unchanged.
program mcsnow_joint_drag_driver
  use mo_constants, only: wp, mp, pi
  use mo_sp_types, only: t_sp
  use mo_atmo_types, only: t_atmo, Rd, Rv, grav
  use mo_atmo, only: psat_water, dynamic_viscosity
  use mo_sp_nml, only: sp_habit, sp_vterm_id, snowgeo_type, agggeo_type, prolate_intp, rm, Vm
  use mo_mass2diam, only: m2d_particle, spheroid_geometry
  use mo_velocity, only: vterm_snow
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  type(t_sp) :: sp
  type(t_atmo) :: air
  real(wp) :: mass, a, c, volume, phi, values(11)
  integer :: unit, ios, id
  character(len=1024) :: input_path
  call get_command_argument(1,input_path)
  open(newunit=unit,file=trim(input_path),status='old',action='read')
  sp_habit=1; sp_vterm_id=3; snowgeo_type=2; agggeo_type=1; prolate_intp=.true.
  rm=1e-5_wp; Vm=0._wp
  do
    read(unit,*,iostat=ios) id,air%T,air%p,mass,a,c
    if (ios<0) exit
    if (ios/=0) error stop 'Malformed input'
    if (.not.all(ieee_is_finite([air%T,air%p,mass,a,c]))) error stop 'Nonfinite input'
    if (min(air%T,air%p,mass,a,c)<=0) error stop 'Nonpositive input'
    air%pv=psat_water(air%T)
    air%rho=(air%p-air%pv)/(Rd*air%T)+air%pv/(Rv*air%T)
    air%eta=dynamic_viscosity(air%T)
    sp%statusb=0; sp%mm=1._mp
    sp%m_w=0._wp; sp%m_r=0._wp; sp%m_f=0._wp; sp%V_r=0._wp
    sp%m_i=mass; sp%V_i=pi/6._wp*a*a*c; sp%phi=c/a
    call m2d_particle(sp)
    call spheroid_geometry(sp,a,c,volume)
    values=[air%rho,air%eta,grav,mass,2*a,2*c,volume,sp%phi,sp%d,sp%p_area,vterm_snow(air,sp)]
    if (.not.all(ieee_is_finite(values))) error stop 'Nonfinite output'
    write(*,'(i0,1x,*(es25.17e3,1x))') id,values
  end do
  close(unit)
end program mcsnow_joint_drag_driver
