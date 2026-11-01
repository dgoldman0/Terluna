! Prescribed-air and steady monodisperse-column adapter to pinned McSnow physics.
program mcsnow_lunar_driver
  use mo_constants, only: wp, mp, pi, i8
  use mo_sp_types, only: t_sp
  use mo_atmo_types, only: t_atmo, Rd, Rv, rhoi, grav, T_3
  use mo_atmo, only: psat_water, psat_ice, dynamic_viscosity
  use mo_atmo_nml, only: IGF_id, read_IGF, get_IGF
  use mo_sp_nml, only: sp_habit, sp_vterm_id, snowgeo_type, agggeo_type, prolate_intp, rm, Vm
  use mo_mass2diam, only: m2d_particle, spheroid_geometry
  use mo_velocity, only: vterm_snow
  use mo_depo, only: sp_depo, growth_trace, supply_rate, requested_rate
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  type(t_sp) :: particle, trial
  type(t_atmo) :: air
  real(wp) :: temperature_K, pressure_Pa, vapor_pressure_Pa, dt_s, duration_s, output_interval_s
  real(wp) :: initial_mass_kg, initial_equatorial_m, initial_polar_m
  real(wp) :: cloud_depth_m, freezing_pressure_Pa, pressure_exponent, lapse_K_m, rh_water
  real(wp) :: water_supply_kg_m2_s, number_flux_m2_s, max_mass_fraction, diameter_stop_m
  real(wp) :: elapsed, height, step, next_output, last_output, initial_mass, deposited_mass
  real(wp) :: q_heat, heat_rate, q_trial, heat_trial, dep_rate, dep_mono
  real(wp) :: speed, r_a, r_c, volume, dm, fraction, maximum_fraction, peak_diameter
  real(wp) :: evaluation_time, evaluation_height, evaluation_T, evaluation_p, accepted_trace(22)
  real(wp) :: last_requested, last_admitted, budget
  integer :: growth_mode, column_mode, unit, trace_unit, info_unit, steps, rejected, limited, status
  logical :: prolate_interpolation, native_saturation, growth_enabled, delete_particle
  character(len=1024) :: input_path
  namelist /crystal/ temperature_K, pressure_Pa, vapor_pressure_Pa, dt_s, duration_s, &
    output_interval_s, initial_mass_kg, initial_equatorial_m, initial_polar_m, growth_mode, &
    prolate_interpolation, native_saturation, column_mode, cloud_depth_m, freezing_pressure_Pa, &
    pressure_exponent, lapse_K_m, rh_water, water_supply_kg_m2_s, number_flux_m2_s, &
    max_mass_fraction, diameter_stop_m, growth_enabled

  column_mode=0; cloud_depth_m=0._wp; freezing_pressure_Pa=0._wp; pressure_exponent=0._wp
  lapse_K_m=0._wp; rh_water=1._wp; water_supply_kg_m2_s=-1._wp; number_flux_m2_s=1._wp
  max_mass_fraction=0._wp; diameter_stop_m=1._wp; growth_enabled=.true.; native_saturation=.true.
  call get_command_argument(1,input_path)
  open(newunit=unit,file=trim(input_path),status='old',action='read')
  read(unit,nml=crystal)
  close(unit)
  if (min(dt_s,duration_s,output_interval_s,initial_equatorial_m,initial_polar_m)<=0) error stop 'Invalid input'
  if (column_mode<0.or.column_mode>1.or.growth_mode/=2) error stop 'Unsupported configuration'
  if (column_mode==1.and.min(cloud_depth_m,freezing_pressure_Pa,number_flux_m2_s)<=0) error stop 'Invalid column'
  if (max_mass_fraction<0.or.max_mass_fraction>1) error stop 'Invalid refinement control'
  if (column_mode==0.and.water_supply_kg_m2_s>=0) error stop 'Supply requires a column'
  IGF_id=growth_mode
  call read_IGF()
  sp_habit=1; sp_vterm_id=3; snowgeo_type=2; agggeo_type=1; prolate_intp=prolate_interpolation
  rm=1e-5_wp; Vm=0._wp
  height=cloud_depth_m
  call set_air()
  particle%xi=1_i8; particle%mm=1._mp
  particle%z=0._wp; particle%cellId=1; particle%start_pe=0; particle%local_idx=1
  particle%T=air%T; particle%init_temp=air%T; particle%init_ssat=air%ssat
  particle%V_i=pi/6._wp*initial_equatorial_m**2*initial_polar_m
  particle%phi=initial_polar_m/initial_equatorial_m
  particle%m_i=initial_mass_kg
  if (initial_mass_kg<0) particle%m_i=particle%V_i*rhoi
  if (particle%m_i<=0.or.particle%m_i>particle%V_i*rhoi*(1._wp+1e-12_wp)) error stop 'Invalid density'
  initial_mass=particle%m_i
  elapsed=0._wp; q_heat=0._wp; heat_rate=0._wp; deposited_mass=0._wp
  steps=0; rejected=0; limited=0; status=0; maximum_fraction=0._wp; peak_diameter=0._wp
  next_output=output_interval_s; last_output=-1._wp
  open(newunit=unit,file='trajectory.csv',status='new',action='write')
  write(unit,'(a)') 'elapsed_s,height_above_freezing_m,temperature_K,pressure_Pa,vapor_pressure_Pa,'// &
    'mass_kg,equatorial_dimension_m,polar_dimension_m,diameter_m,density_kg_m3,ice_volume_m3,'// &
    'aspect_ratio,projected_area_m2,terminal_speed_m_s,integrated_deposition_kg,'// &
    'ice_mass_balance_residual_kg,latent_heat_J,allocated_water_kg'
  open(newunit=trace_unit,file='step_trace.csv',status='new',action='write')
  write(trace_unit,'(a)') 'evaluation_elapsed_s,state_elapsed_s,step_s,height_above_freezing_m,'// &
    'temperature_K,pressure_Pa,fall_speed_m_s,requested_rate_kg_s,admitted_rate_kg_s,supply_rate_kg_s,'// &
    'ventilation_vapor,effective_growth_factor,deposition_density_kg_m3'
  open(newunit=info_unit,file='air.csv',status='new',action='write')
  write(info_unit,'(a)') 'temperature_K,pressure_Pa,vapor_pressure_Pa,rh_percent,density_kg_m3,'// &
    'native_gravity_m_s2,native_ice_density_kg_m3,inherent_growth_factor'
  write(info_unit,'(*(es25.17e3,:,","))') air%T,air%p,air%pv,air%rh,air%rho,grav,rhoi,get_IGF(air%T)
  close(info_unit)
  call write_state()
  do while (elapsed<duration_s)
    call set_air()
    call m2d_particle(particle)
    speed=vterm_snow(air,particle)
    evaluation_time=elapsed; evaluation_height=height; evaluation_T=air%T; evaluation_p=air%p
    step=min(dt_s,duration_s-elapsed,next_output-elapsed)
    if (column_mode==1) step=min(step,height/speed)
    if (step<=0._wp) error stop 'Nonpositive integration step'
    supply_rate=-1._wp
    if (water_supply_kg_m2_s>=0) supply_rate=water_supply_kg_m2_s/cloud_depth_m*speed/number_flux_m2_s
    do
      trial=particle; q_trial=q_heat; heat_trial=heat_rate; dep_rate=0._wp; dep_mono=0._wp
      requested_rate=0._wp; growth_trace=0._wp
      if (growth_enabled) then
        call sp_depo(trial,air,step,q_trial,delete_particle,dep_rate,dep_mono,heat_trial)
        if (delete_particle) error stop 'Particle lost'
      end if
      dm=trial%m_i-particle%m_i
      fraction=abs(dm)/particle%m_i
      if (max_mass_fraction==0.or.fraction<=max_mass_fraction) exit
      step=step*0.5_wp; rejected=rejected+1
      if (step<1e-10_wp) error stop 'Step refinement exhausted'
    end do
    accepted_trace=growth_trace; last_requested=requested_rate; last_admitted=dm/step
    if (supply_rate>=0.and.requested_rate>supply_rate) limited=limited+1
    particle=trial; q_heat=q_trial; heat_rate=heat_trial
    deposited_mass=deposited_mass+dep_rate*step
    maximum_fraction=max(maximum_fraction,fraction)
    elapsed=elapsed+step; steps=steps+1
    if (steps>10000000) error stop 'Step-count resource limit'
    if (column_mode==1) height=max(0._wp,height-speed*step)
    call m2d_particle(particle)
    peak_diameter=max(peak_diameter,particle%d)
    if (.not.all(ieee_is_finite([particle%m_i,particle%d,particle%V_i,particle%phi]))) error stop 'Nonfinite state'
    if (column_mode==1.and.height<=1e-9_wp) status=1
    if (diameter_stop_m>0.and.particle%d>=diameter_stop_m) status=2
    if (elapsed>=next_output.or.elapsed>=duration_s.or.status/=0.or.steps==1) then
      write(trace_unit,'(*(es25.17e3,:,","))') evaluation_time,elapsed,step,evaluation_height, &
        evaluation_T,evaluation_p,speed,last_requested,last_admitted,supply_rate, &
        accepted_trace(11),accepted_trace(14),accepted_trace(15)
    end if
    if (elapsed>=next_output.or.elapsed>=duration_s.or.status/=0) then
      call set_air()
      call write_state()
      next_output=next_output+output_interval_s
    end if
    if (status/=0) exit
  end do
  close(unit); close(trace_unit)
  open(newunit=info_unit,file='summary.csv',status='new',action='write')
  write(info_unit,'(a)') 'elapsed_s,status_code,accepted_steps,rejected_steps,water_limited_steps,'// &
    'maximum_step_mass_fraction,peak_diameter_m,initial_mass_kg,integrated_deposition_kg,'// &
    'latent_heat_J,height_above_freezing_m,allocated_water_kg'
  budget=-1._wp
  if (water_supply_kg_m2_s>=0) budget=water_supply_kg_m2_s/number_flux_m2_s*(cloud_depth_m-height)/cloud_depth_m
  write(info_unit,'(*(es25.17e3,:,","))') elapsed,real(status,wp),real(steps,wp),real(rejected,wp), &
    real(limited,wp),maximum_fraction,peak_diameter,initial_mass,deposited_mass,heat_rate,height,budget
  close(info_unit)
contains
  subroutine set_air()
    air%T=temperature_K; air%p=pressure_Pa; air%pv=vapor_pressure_Pa
    if (column_mode==1) then
      air%T=T_3-lapse_K_m*height
      air%p=freezing_pressure_Pa*(air%T/T_3)**pressure_exponent
    end if
    air%psatw=psat_water(air%T); air%psati=psat_ice(air%T)
    if (native_saturation.or.column_mode==1) air%pv=rh_water*air%psatw
    air%rh=100._wp*air%pv/air%psatw; air%ssat=air%pv/air%psati-1._wp
    air%eta=dynamic_viscosity(air%T); air%etai=1._wp/air%eta
    air%rho=(air%p-air%pv)/(Rd*air%T)+air%pv/(Rv*air%T)
    air%qv=(Rd/Rv)*air%pv/(air%p-(1._wp-Rd/Rv)*air%pv); air%qc=0._wp
  end subroutine set_air
  subroutine write_state()
    call m2d_particle(particle)
    call spheroid_geometry(particle,r_a,r_c,volume)
    budget=-1._wp
    if (water_supply_kg_m2_s>=0) budget=water_supply_kg_m2_s/number_flux_m2_s*(cloud_depth_m-height)/cloud_depth_m
    write(unit,'(*(es25.17e3,:,","))') elapsed,height,air%T,air%p,air%pv,particle%m_i,2*r_a,2*r_c, &
      particle%d,particle%m_i/volume,particle%V_i,particle%phi,particle%p_area,vterm_snow(air,particle), &
      deposited_mass,particle%m_i-initial_mass-deposited_mass,heat_rate,budget
    flush(unit)
    last_output=elapsed
  end subroutine write_state
end program mcsnow_lunar_driver
