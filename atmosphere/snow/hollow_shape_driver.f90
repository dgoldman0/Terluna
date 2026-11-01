! Local C interface to the pinned Harrington--Pokrifka 2026 routines.
! No published fitting coefficients are copied into this interface.
subroutine hollow_shape_eval(phi_in, habit, coordinates, result) bind(C)
  use iso_c_binding
  implicit none
  real(c_double), value :: phi_in
  integer(c_int), value :: habit, coordinates
  real(c_double), intent(out) :: result(8)
  real :: phi, hafun, hcfun, hc10fun, hc01fun, ha10fun, ha01fun
  real :: hacapfun, hccapfun
  integer :: ihab, iwd
  phi=real(phi_in); ihab=habit; iwd=coordinates
  result(1)=hafun(phi,ihab,iwd)
  result(2)=hcfun(phi,ihab,iwd)
  result(3)=ha10fun(phi,ihab,iwd)
  result(4)=hc10fun(phi,ihab,iwd)
  result(5)=ha01fun(phi,ihab,iwd)
  result(6)=hc01fun(phi,ihab,iwd)
  result(7)=hacapfun(phi,ihab)
  result(8)=hccapfun(phi,ihab)
end subroutine hollow_shape_eval
