! Repeated prescribed-air trajectories; growth, geometry and drag use native McSnow.
! stdin: six diagnostic controls. stdout: one row per requested state, then DONE.
program mcsnow_tune_driver
  use mo_constants, only: wp, mp, pi, i8
  use mo_sp_types, only: t_sp
  use mo_atmo_types, only: t_atmo, Rd, Rv, rhoi
  use mo_atmo, only: psat_water, psat_ice, dynamic_viscosity
  use mo_atmo_nml, only: IGF_id, read_IGF
  use mo_sp_nml, only: sp_habit, sp_vterm_id, snowgeo_type, agggeo_type, prolate_intp, rm, Vm
  use mo_mass2diam, only: m2d_particle, spheroid_geometry
  use mo_velocity, only: vterm_snow
  use mo_depo, only: sp_depo, tune, tune_cap_count
  use, intrinsic :: iso_fortran_env, only: input_unit, output_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  type(t_sp) :: particle
  type(t_atmo) :: air
  real(wp), allocatable :: settings(:,:), times(:,:)
  integer, allocatable :: counts(:)
  real(wp) :: dt, q_heat, heat_rate, dep_rate, dep_mono, deposited, initial_mass, a, c, volume, speed
  integer :: ncase, i, j, k, last_step, next_step, unit, ios
  logical :: deleted
  character(len=1024) :: input_path
  call get_command_argument(1,input_path)
  open(newunit=unit,file=trim(input_path),status='old',action='read')
  read(unit,*) ncase
  allocate(settings(7,ncase),times(8,ncase),counts(ncase))
  do i=1,ncase
    ! temperature, pressure, water saturation ratio, mass, full a, full c, dt, nout
    read(unit,*) settings(:,i),counts(i)
    if (counts(i)<1.or.counts(i)>8) error stop 'Invalid output count'
    read(unit,*) times(1:counts(i),i)
  end do
  close(unit)
  IGF_id=2
  call read_IGF()
  sp_habit=1; sp_vterm_id=3; snowgeo_type=2; agggeo_type=1; prolate_intp=.true.
  rm=1e-5_wp; Vm=0._wp
  write(output_unit,'(a)') 'READY'
  flush(output_unit)
  do
    read(input_unit,*,iostat=ios) tune
    if (ios<0) exit
    if (ios/=0.or.any(.not.ieee_is_finite(tune))) error stop 'Invalid controls'
    do i=1,ncase
      air%T=settings(1,i); air%p=settings(2,i)
      air%psatw=psat_water(air%T); air%psati=psat_ice(air%T)
      air%pv=settings(3,i)*air%psatw
      air%rh=100._wp*air%pv/air%psatw; air%ssat=air%pv/air%psati-1._wp
      air%eta=dynamic_viscosity(air%T); air%etai=1._wp/air%eta
      air%rho=(air%p-air%pv)/(Rd*air%T)+air%pv/(Rv*air%T)
      air%qv=(Rd/Rv)*air%pv/(air%p-(1._wp-Rd/Rv)*air%pv); air%qc=0._wp
      ! Reset all evolved fields before every case and every parameter vector.
      particle%m_w=0._wp; particle%m_r=0._wp; particle%m_f=0._wp; particle%v_r=0._wp
      particle%statusb=0; particle%xi=1_i8; particle%mm=1._mp
      particle%z=0._wp; particle%cellId=1; particle%start_pe=0; particle%local_idx=1
      particle%T=air%T; particle%init_temp=air%T; particle%init_ssat=air%ssat
      particle%V_i=pi/6._wp*settings(5,i)**2*settings(6,i)
      particle%phi=settings(6,i)/settings(5,i); particle%m_i=settings(4,i)
      if (particle%m_i<0._wp) particle%m_i=particle%V_i*rhoi
      if (particle%m_i<=0.or.particle%m_i>particle%V_i*rhoi*(1._wp+1e-12_wp)) error stop 'Invalid density'
      initial_mass=particle%m_i; dt=settings(7,i)
      q_heat=0._wp; heat_rate=0._wp; deposited=0._wp; tune_cap_count=0
      last_step=0
      do j=1,counts(i)
        next_step=nint(times(j,i)/dt)
        if (abs(next_step*dt-times(j,i))>1e-8_wp.or.next_step<last_step) error stop 'Invalid time'
        do k=last_step+1,next_step
          dep_rate=0._wp; dep_mono=0._wp
          call sp_depo(particle,air,dt,q_heat,deleted,dep_rate,dep_mono,heat_rate)
          if (deleted) error stop 'Particle deleted'
          deposited=deposited+dep_rate*dt
        end do
        last_step=next_step
        call m2d_particle(particle)
        call spheroid_geometry(particle,a,c,volume)
        speed=vterm_snow(air,particle)
        if (.not.all(ieee_is_finite([particle%m_i,a,c,volume,speed]))) error stop 'Nonfinite state'
        if (min(particle%m_i,a,c,volume)<=0.or.speed<0) error stop 'Invalid state'
        write(output_unit,'(*(es25.17e3,1x))') particle%m_i,2*a,2*c,speed, &
          (particle%m_i-initial_mass-deposited)/particle%m_i, &
          volume/(pi/6._wp*(2*a)**2*(2*c))-1._wp,real(tune_cap_count,wp), &
          particle%m_i/volume/rhoi
      end do
    end do
    write(output_unit,'(a)') 'DONE'
    flush(output_unit)
  end do
end program mcsnow_tune_driver
