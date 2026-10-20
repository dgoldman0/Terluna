! CM1's RRTMG longwave (rrtmg_lw, as the lunar build compiles it) on columns read from standard input, clear sky.
! Input: ncol nlay cp grav; then per column the ground temperature, the nlay+1 interface pressures (hPa) and
! temperatures (K), ground first, and the nlay layer pressures, temperatures and H2O, CO2 and O2 mole fractions.
! Output: per column the heating of each layer (K/day), then the upward and downward fluxes (W/m2) at the nlay+1
! interfaces, all clear sky.
program lw_columns
  use parkind, only: im => kind_im, rb => kind_rb
  use parrrtm, only: nbndlw, ngptlw
  use rrtmg_lw_init, only: rrtmg_lw_ini
  use rrtmg_lw_rad, only: rrtmg_lw
  use module_ra_rrtmg_lw, only: rrtmg_lwlookuptable
  implicit none
  integer(im) :: ncol, nlay, icld, i, k
  real(rb) :: cp, grav
  real(rb), allocatable :: play(:,:), plev(:,:), tlay(:,:), tlev(:,:), tsfc(:), h2o(:,:), co2(:,:), o2(:,:)
  real(rb), allocatable :: zero(:,:), emis(:,:), re(:,:), cld(:,:,:), tauaer(:,:,:)
  real(rb), allocatable :: uflx(:,:), dflx(:,:), hr(:,:), uflxc(:,:), dflxc(:,:), hrc(:,:), ucln(:,:), dcln(:,:)
  read(*,*) ncol, nlay, cp, grav
  allocate(play(ncol,nlay), plev(ncol,nlay+1), tlay(ncol,nlay), tlev(ncol,nlay+1), tsfc(ncol))
  allocate(h2o(ncol,nlay), co2(ncol,nlay), o2(ncol,nlay), zero(ncol,nlay), emis(ncol,nbndlw), re(ncol,nlay))
  allocate(cld(ngptlw,ncol,nlay), tauaer(ncol,nlay,nbndlw))
  allocate(uflx(ncol,nlay+1), dflx(ncol,nlay+1), hr(ncol,nlay), uflxc(ncol,nlay+1), dflxc(ncol,nlay+1), hrc(ncol,nlay))
  allocate(ucln(ncol,nlay+1), dcln(ncol,nlay+1))
  do i = 1, ncol
    read(*,*) tsfc(i)
    read(*,*) (plev(i,k), k = 1, nlay+1)
    read(*,*) (tlev(i,k), k = 1, nlay+1)
    read(*,*) (play(i,k), k = 1, nlay)
    read(*,*) (tlay(i,k), k = 1, nlay)
    read(*,*) (h2o(i,k), k = 1, nlay)
    read(*,*) (co2(i,k), k = 1, nlay)
    read(*,*) (o2(i,k), k = 1, nlay)
  end do
  zero = 0.0_rb
  emis = 1.0_rb
  re = 10.0_rb
  cld = 0.0_rb
  tauaer = 0.0_rb
  icld = 0
  call rrtmg_lwlookuptable                      ! the tables from RRTMG_LW_DATA, as CM1's rrtmg_lwinit reads them first
  call rrtmg_lw_ini(cp, grav)
  call rrtmg_lw(ncol, nlay, icld, play, plev, tlay, tlev, tsfc, h2o, zero, co2, zero, zero, o2, &
                zero, zero, zero, zero, emis, 0, 0, 0, cld, cld, cld, cld, cld, re, re, re, tauaer, &
                uflx, dflx, hr, uflxc, dflxc, hrc, ucln, dcln, 0)
  do i = 1, ncol
    write(*,'(*(es16.8))') (hr(i,k), k = 1, nlay)          ! without cloud (icld 0) the total sky is the clear sky,
    write(*,'(*(es16.8))') (uflx(i,k), k = 1, nlay+1)      ! which RRTMG fills; its separate clear-sky arrays are
    write(*,'(*(es16.8))') (dflx(i,k), k = 1, nlay+1)      ! filled only with cloud on
  end do
end program lw_columns
