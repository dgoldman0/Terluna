! Batch interface to native McSnow for diagnostic crystals in prescribed air.
module mcsnow_weather_driver
  use iso_c_binding, only: c_double, c_int
  use mo_constants, only: wp, mp, pi, i8
  use mo_sp_types, only: t_sp
  use mo_atmo_types, only: t_atmo, Rd, Rv, rhoi, grav
  use mo_atmo, only: psat_water, psat_ice, dynamic_viscosity
  use mo_atmo_nml, only: IGF_id, read_IGF
  use mo_sp_nml, only: sp_habit, sp_vterm_id, snowgeo_type, agggeo_type, prolate_intp, rm, Vm
  use mo_mass2diam, only: m2d_particle, spheroid_geometry
  use mo_velocity, only: vterm_snow
  use mo_depo, only: sp_depo, supply_rate
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
contains
  subroutine weather_init() bind(C)
    IGF_id=2
    call read_IGF()
    sp_habit=1; sp_vterm_id=3; snowgeo_type=2; agggeo_type=1; prolate_intp=.true.
    rm=1e-5_wp; Vm=0._wp; supply_rate=-1._wp
  end subroutine

  subroutine weather_constants(values) bind(C)
    real(c_double),intent(out) :: values(4)
    values=[rhoi,Rd,Rv,grav]
  end subroutine

  subroutine weather_evolve(n,state,forcing,duration,saturated,max_step,max_fraction,min_mass,max_diameter,output,status) bind(C)
    integer(c_int),value :: n,saturated
    real(c_double),intent(inout) :: state(5,n)
    real(c_double),intent(in) :: forcing(3,n),duration(n)
    real(c_double),value :: max_step,max_fraction,min_mass,max_diameter
    real(c_double),intent(out) :: output(14,n)
    integer(c_int),intent(out) :: status(n)
    type(t_sp) :: sp,trial
    type(t_atmo) :: air
    real(wp) :: dt,elapsed,dm,dep,dep_mono,q,heat,total_heat,total_dep,a,c,vol,peak,m0
    logical :: removed
    integer :: i,steps
    do i=1,n
      output(:,i)=0._wp;status(i)=0;elapsed=0._wp;total_heat=0._wp;total_dep=0._wp;steps=0
      if (minval(state(1:3,i))<=0.or.minval(forcing(:,i))<=0.or.duration(i)<0) then
        status(i)=-1;cycle
      end if
      air%T=forcing(1,i);air%p=forcing(2,i);air%pv=forcing(3,i)
      air%psatw=psat_water(air%T);air%psati=psat_ice(air%T)
      if (saturated==1) air%pv=air%psatw
      air%rh=100._wp*air%pv/air%psatw;air%ssat=air%pv/air%psati-1._wp
      air%eta=dynamic_viscosity(air%T);air%etai=1._wp/air%eta
      air%rho=(air%p-air%pv)/(Rd*air%T)+air%pv/(Rv*air%T)
      air%qv=(Rd/Rv)*air%pv/(air%p-(1._wp-Rd/Rv)*air%pv);air%qc=0._wp
      sp%m_w=0._wp;sp%m_r=0._wp;sp%v_r=0._wp;sp%m_f=0._wp;sp%statusb=0
      sp%xi=1_i8;sp%mm=1._mp;sp%z=0._wp;sp%cellId=1;sp%start_pe=0;sp%local_idx=1
      sp%m_i=state(1,i);sp%V_i=state(2,i);sp%phi=state(3,i)
      sp%T=air%T;sp%init_temp=state(4,i);sp%init_ssat=state(5,i)
      call m2d_particle(sp)
      peak=sp%d;m0=sp%m_i
      do while(elapsed<duration(i))
        if (sp%m_i<=min_mass) then
          status(i)=1;exit
        end if
        if (sp%d>=max_diameter) then
          status(i)=2;exit
        end if
        dt=min(max_step,duration(i)-elapsed)
        do
          trial=sp;q=0._wp;heat=0._wp;dep=0._wp;dep_mono=0._wp
          call sp_depo(trial,air,dt,q,removed,dep,dep_mono,heat)
          dm=trial%m_i-sp%m_i
          if (.not.removed.and.abs(dm)<=max_fraction*sp%m_i) exit
          dt=dt/2
          if (dt<1e-12_wp) then
            status(i)=-2;exit
          end if
        end do
        if (status(i)<0) exit
        sp=trial;elapsed=elapsed+dt;steps=steps+1
        total_dep=total_dep+dm;total_heat=total_heat+heat
        call m2d_particle(sp)
        peak=max(peak,real(sp%d,wp))
        if (steps>10000000.or..not.all(ieee_is_finite([sp%m_i,sp%V_i,sp%phi,real(sp%d,wp)]))) then
          status(i)=-3;exit
        end if
      end do
      if (sp%m_i<=min_mass) status(i)=1
      if (sp%d>=max_diameter) status(i)=2
      call spheroid_geometry(sp,a,c,vol)
      state(:,i)=[sp%m_i,sp%V_i,real(sp%phi,wp),real(sp%init_temp,wp),real(sp%init_ssat,wp)]
      output(:,i)=[real(sp%d,wp),2*a,2*c,sp%m_i/vol,vterm_snow(air,sp), &
        total_dep,total_heat,peak,elapsed,real(steps,wp),air%psatw,air%psati,air%pv,sp%m_i-m0-total_dep]
    end do
  end subroutine
end module mcsnow_weather_driver
