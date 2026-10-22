! Guarded C interface to the separately cached Harrington et al. (2021) routines.
! This file contains no copied polynomial fits. Source air-property formulas are
! repeated below at the native REAL precision so a caller can audit coupling.
module facet_attachment_interface
  use, intrinsic :: iso_c_binding, only: c_double, c_int
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: ieee_exceptions
  implicit none
contains
  subroutine facet_attachment_eval(inputs, outputs, status) bind(C)
    real(c_double), intent(in) :: inputs(6)
    real(c_double), intent(out) :: outputs(20)
    integer(c_int), intent(out) :: status
    real :: si, tk, pressure, a, c, m, aa, ac, at, req
    real :: slope, ar, cr, ur, nu1, nu2, nd1, nd2
    real :: cap, capdelta, dv, ka, esi, lv, lf, ls, speed, density
    real :: mfp, delta, thermal_length, sa, sc, pmb
    real, external :: capacitance, ei
    real, parameter :: pi=3.14159, mw=.01802, runiv=8.314, rv=461.5, rd=287., cp=1004.
    type(ieee_status_type) :: saved_status
    logical :: bad_invalid, bad_zero, bad_overflow

    outputs=0._c_double
    status=1_c_int
    if (.not.all(ieee_is_finite(inputs))) return
    if (inputs(1)<=0.or.inputs(2)<273.15_c_double-70._c_double.or.inputs(2)>273.15_c_double-1._c_double) return
    if (minval(inputs(3:5))<=0.or.inputs(6)<1.or.inputs(6)>15) return
    call ieee_get_status(saved_status)
    call ieee_set_halting_mode(ieee_invalid,.false.)
    call ieee_set_halting_mode(ieee_divide_by_zero,.false.)
    call ieee_set_halting_mode(ieee_overflow,.false.)
    call ieee_set_flag(ieee_all,.false.)
    si=real(inputs(1)); tk=real(inputs(2)); pressure=real(inputs(3))
    a=real(inputs(4)); c=real(inputs(5)); m=real(inputs(6))
    if (.not.all(ieee_is_finite([si,tk,pressure,a,c,m]))) then
      call ieee_set_status(saved_status)
      return
    endif
    if (min(si,pressure,a,c)<=0.) then
      call ieee_set_status(saved_status)
      return
    endif
    req=(a**2*c)**(1./3.)
    call getfittingparams(m,slope,ar,cr,ur,nu1,nu2,nd1,nd2)
    call getdepcoeffparam(si,tk,pressure,a,c,m,slope,ar,cr,ur,nu1,nu2,nd1,nd2,0,req,aa,ac,at)

    ! These repeat the source formulas, including its distinct mean-free-path
    ! constants. They expose conventions; they do not change coefficient inputs.
    pmb=pressure/100.
    dv=2.11e-5*(tk/273.15)**1.94*1013.25/pmb
    ka=2.3823e-2+7.1177e-5*(tk-273.15)
    lv=2.5008e6*(273.15/tk)**(0.167+3.67e-4*tk)
    lf=0.3337e6+(2031.0-10.467*(tk-273.15))*(tk-273.15)
    ls=lv+lf
    esi=ei(tk)
    speed=sqrt(8.0*runiv*tk/(pi*mw))
    density=pressure/(rd*tk)
    mfp=6.6e-8*(1013.25/pmb)*(tk/273.15)
    delta=6.6e-6*(1013.15/pmb)*(tk/293.15)/100.
    thermal_length=ka/(at*density*cp)*sqrt(2.*pi/(rd*tk))
    cap=capacitance(a,c)
    capdelta=capacitance(a+delta,c+delta)
    call getscrit(tk,sa,sc)
    outputs=real([aa,ac,at,cap,capdelta,dv,ka,esi,ls,speed,density,mfp,delta, &
                  thermal_length,sa,sc,rv,rd,cp,real(precision(1.))],c_double)
    call ieee_get_flag(ieee_invalid,bad_invalid)
    call ieee_get_flag(ieee_divide_by_zero,bad_zero)
    call ieee_get_flag(ieee_overflow,bad_overflow)
    status=0_c_int
    if (bad_invalid.or.bad_zero.or.bad_overflow.or..not.all(ieee_is_finite(outputs))) status=2_c_int
    if (min(aa,ac,at)<=0.or.max(aa,ac,at)>1.or.min(sa,sc)<=0.) status=2_c_int
    if (min(cap,capdelta,dv,ka,esi,ls,speed,density)<=0.) status=2_c_int
    ! Domain monotonicity: capacitance lies between inscribed and enclosing
    ! sphere radii, and grows when both semiaxes increase. Detect native
    ! near-sphere cancellation; allow only machine-rounding slack.
    if (cap<min(a,c)*(1.-8.*epsilon(a)).or.cap>max(a,c)*(1.+8.*epsilon(a))) status=2_c_int
    if (capdelta<min(a+delta,c+delta)*(1.-8.*epsilon(a))) status=2_c_int
    if (capdelta>max(a+delta,c+delta)*(1.+8.*epsilon(a))) status=2_c_int
    if (capdelta<cap*(1.-8.*epsilon(a))) status=2_c_int
    call ieee_set_status(saved_status)
  end subroutine facet_attachment_eval

  subroutine facet_attachment_demo_state(index, inputs) bind(C)
    ! Match original demonstration input arithmetic, including single rounding.
    integer(c_int), value :: index
    real(c_double), intent(out) :: inputs(6)
    real :: a, c, phi, t0, temp, si, pressure, m
    m=10.; phi=1.; t0=273.15; temp=-15.+t0
    pressure=10000.; si=5./100
    a=float(index)*1.e-6; c=a*phi
    inputs=real([si,temp,pressure,a,c,m],c_double)
  end subroutine facet_attachment_demo_state
end module facet_attachment_interface
