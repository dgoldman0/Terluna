! Thin prescribed-air driver; all growth, habit, geometry and drag stay native.
program mcsnow_crystal_driver
  use mo_constants, only: wp, mp, pi, i8
  use mo_sp_types, only: t_sp
  use mo_atmo_types, only: t_atmo, Rd, Rv, rhoi, grav
  use mo_atmo, only: psat_water, psat_ice, dynamic_viscosity, diffusivity_vapour
  use mo_atmo_nml, only: IGF_id, read_IGF, get_IGF
  use mo_sp_nml, only: sp_habit, sp_vterm_id, snowgeo_type, agggeo_type, prolate_intp, rm, Vm
  use mo_mass2diam, only: m2d_particle, spheroid_geometry
  use mo_velocity, only: vterm_snow
  use mo_depo, only: sp_depo
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  type(t_sp) :: particle
  type(t_atmo) :: air
  real(wp) :: temperature_K, pressure_Pa, vapor_pressure_Pa, dt_s, duration_s, output_interval_s
  real(wp) :: initial_mass_kg, initial_equatorial_m, initial_polar_m
  real(wp) :: q_heat, dep_rate, dep_mono, heat_rate, deposited_mass, initial_mass, elapsed
  real(wp) :: r_a, r_c, volume, fall_speed, old_mass, maximum_step_mass_fraction
  real(wp) :: branching_ratio
  integer :: branch_unit
  logical :: previous_branch, current_branch
  integer :: growth_mode, nsteps, output_steps, k, unit
  logical :: delete_particle, prolate_interpolation
  character(len=1024) :: input_path
  namelist /crystal/ temperature_K, pressure_Pa, vapor_pressure_Pa, dt_s, duration_s, &
    output_interval_s, initial_mass_kg, initial_equatorial_m, initial_polar_m, growth_mode, prolate_interpolation

  call get_command_argument(1,input_path)
  if (len_trim(input_path)==0) error stop 'Input namelist required'
  open(newunit=unit,file=trim(input_path),status='old',action='read')
  read(unit,nml=crystal)
  close(unit)
  if (dt_s<=0.or.duration_s<=0.or.output_interval_s<=0) error stop 'Nonpositive time'
  nsteps=nint(duration_s/dt_s); output_steps=nint(output_interval_s/dt_s)
  if (abs(nsteps*dt_s-duration_s)>1e-8_wp.or.output_steps<1) error stop 'Invalid time division'
  if (abs(output_steps*dt_s-output_interval_s)>1e-8_wp) error stop 'Invalid output division'
  if (initial_equatorial_m<=0.or.initial_polar_m<=0) error stop 'Invalid axes'
  if (growth_mode<1.or.growth_mode>2) error stop 'Unsupported growth mode'
  IGF_id=growth_mode
  call read_IGF()
  sp_habit=1; sp_vterm_id=3; snowgeo_type=2; agggeo_type=1
  prolate_intp=prolate_interpolation
  rm=1e-5_wp; Vm=0._wp

  air%T=temperature_K; air%p=pressure_Pa; air%pv=vapor_pressure_Pa
  air%psatw=psat_water(air%T); air%psati=psat_ice(air%T)
  air%rh=100._wp*air%pv/air%psatw; air%ssat=air%pv/air%psati-1._wp
  air%eta=dynamic_viscosity(air%T); air%etai=1._wp/air%eta
  air%rho=(air%p-air%pv)/(Rd*air%T)+air%pv/(Rv*air%T)
  air%qv=(Rd/Rv)*air%pv/(air%p-(1._wp-Rd/Rv)*air%pv); air%qc=0._wp

  particle%xi=1_i8; particle%mm=1._mp
  particle%z=0._wp; particle%cellId=1; particle%start_pe=0; particle%local_idx=1
  particle%T=air%T; particle%init_temp=air%T; particle%init_ssat=air%ssat
  particle%V_i=pi/6._wp*initial_equatorial_m**2*initial_polar_m
  particle%phi=initial_polar_m/initial_equatorial_m
  particle%m_i=initial_mass_kg
  if (initial_mass_kg<0._wp) particle%m_i=particle%V_i*rhoi
  if (particle%m_i<=0.or.particle%m_i>particle%V_i*rhoi*(1._wp+1e-12_wp)) error stop 'Invalid density'
  initial_mass=particle%m_i
  q_heat=0._wp; heat_rate=0._wp; deposited_mass=0._wp; dep_rate=0._wp; dep_mono=0._wp
  open(newunit=unit,file='air.csv',status='new',action='write')
  write(unit,'(a)') 'temperature_K,pressure_Pa,vapor_pressure_Pa,rh_percent,ice_supersaturation,density_kg_m3,'// &
    'dynamic_viscosity_Pa_s,native_ice_density_kg_m3,native_gravity_m_s2,inherent_growth_factor'
  write(unit,'(*(es25.17e3,:,","))') air%T,air%p,air%pv,air%rh,air%ssat,air%rho,air%eta,rhoi,grav,get_IGF(air%T)
  close(unit)
  open(newunit=unit,file='trajectory.csv',status='new',action='write')
  write(unit,'(a)') 'elapsed_s,mass_kg,equatorial_dimension_m,polar_dimension_m,diameter_m,density_kg_m3,'// &
    'ice_volume_m3,aspect_ratio,projected_area_m2,terminal_speed_m_s,last_deposition_rate_kg_s,'// &
    'integrated_deposition_kg,ice_mass_balance_residual_kg,integrated_heat_J,particle_heat_diagnostic_J'
  elapsed=0._wp; maximum_step_mass_fraction=0._wp
  open(newunit=branch_unit,file='branch_diagnostic.csv',status='new',action='write')
  write(branch_unit,'(a)') 'elapsed_s,branch_threshold_ratio,branch_criterion_met'
  call write_state()
  call branch_diagnostic(.true.)
  do k=1,nsteps
    dep_rate=0._wp; dep_mono=0._wp
    old_mass=particle%m_i
    call sp_depo(particle,air,dt_s,q_heat,delete_particle,dep_rate,dep_mono,heat_rate)
    if (delete_particle) error stop 'Particle deleted'
    deposited_mass=deposited_mass+dep_rate*dt_s
    maximum_step_mass_fraction=max(maximum_step_mass_fraction,abs(particle%m_i-old_mass)/old_mass)
    elapsed=k*dt_s
    call branch_diagnostic(.false.)
    if (mod(k,output_steps)==0.or.k==nsteps) call write_state()
  end do
  close(unit)
  close(branch_unit)
  open(newunit=unit,file='step_diagnostic.csv',status='new',action='write')
  write(unit,'(a)') 'maximum_step_mass_fraction,steps,wp_decimal_precision,mp_decimal_precision'
  write(unit,'(es25.17e3,3(a,i0))') maximum_step_mass_fraction,',',nsteps,',',precision(1._wp),',',precision(1._mp)
  close(unit)
contains
  subroutine branch_diagnostic(first)
    ! Diagnostic only: repeats the source branch condition using native gas/geometry.
    ! IGF1 native Dv includes mean ventilation, which is private; report only IGF2+ threshold.
    logical,intent(in) :: first
    if (growth_mode/=2) return
    call m2d_particle(particle)
    branching_ratio=0.5_wp*particle%d/200e-6_wp
    current_branch=branching_ratio>=1._wp
    if (first) then
      write(branch_unit,'(es25.17e3,a,es25.17e3,a,i0)') elapsed,',',branching_ratio,',',merge(1,0,current_branch)
      flush(branch_unit)
    else if (current_branch.neqv.previous_branch) then
      write(branch_unit,'(es25.17e3,a,es25.17e3,a,i0)') elapsed,',',branching_ratio,',',merge(1,0,current_branch)
      flush(branch_unit)
    end if
    previous_branch=current_branch
  end subroutine branch_diagnostic
  subroutine write_state()
    call m2d_particle(particle)
    call spheroid_geometry(particle,r_a,r_c,volume)
    fall_speed=vterm_snow(air,particle)
    if (.not.all(ieee_is_finite([particle%m_i,r_a,r_c,volume,fall_speed,particle%phi]))) error stop 'Nonfinite state'
    if (min(particle%m_i,r_a,r_c,volume,particle%phi)<=0.or.fall_speed<0) error stop 'Invalid state'
    write(unit,'(*(es25.17e3,:,","))') elapsed,particle%m_i,2*r_a,2*r_c,particle%d, &
      particle%m_i/volume,particle%V_i,particle%phi,particle%p_area,fall_speed,dep_rate, &
      deposited_mass,particle%m_i-initial_mass-deposited_mass,heat_rate,q_heat
    flush(unit)
  end subroutine write_state
end program mcsnow_crystal_driver
