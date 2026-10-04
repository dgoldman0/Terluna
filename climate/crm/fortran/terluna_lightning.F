!-----------------------------------------------------------------------------------------------------------------------
! Terluna: the electric field of a CM1 storm's charge and WRF-ELEC's cylindrical lightning (climate/crm/cm1_elec.py)
!
! The field: Poisson's equation for the potential of the net charge density q (C/m3), periodic in x and y, zero at the
! ground (the lowest w level) and at the model top, on CM1's stretched levels. Each horizontal Fourier mode, with the
! eigenvalue of the second-order horizontal Laplacian, leaves a tridiagonal system in height; the ghost level below the
! ground and above the top mirrors the potential with the opposite sign. The field is minus the centred gradient.
!
! Lightning: one call of the cylindrical scheme of WRF-ELEC's light1d (module_discharge.F, MicroTed/wrf4-elec e43041b;
! Ziegler and MacGorman 1994, Fierro et al. 2013), with its newer rule for the share of charge removed. Every column within
! the discharge radius of a point where the field exceeds the breakdown field joins the discharge if it holds charge
! density above chgthr anywhere below the top level; in those columns the charge above chgthr is reduced, by shares set
! so that the discharge carries comparable positive and negative charge. The removed charge goes to the small ions with
! the opposite sign, as in WRF-ELEC. Unlike WRF-ELEC's version, x and y spacings may differ and the domain wraps.
!
! Apart from CM1's FFT (module singleton) the module stands alone, so climate/tests/test_cm1_elec_field.py checks it
! against numpy.
!-----------------------------------------------------------------------------------------------------------------------
module terluna_lightning
  implicit none
  private
  public :: potential, field, energy, breakdown_field, discharge

  real, parameter, public :: eps0 = 8.8541878128e-12   ! F/m
  real, parameter, public :: chgthr = 0.1e-9           ! C/m3: charge density discharges act on (light1d's chgthr)
  real, parameter, public :: disfrac = 0.3             ! nominal share of charge a discharge removes (nssl_disfrac)

contains

  subroutine potential(ni, nj, nk, dx, dy, zh1, zf1, q, phi)
    use singleton, only : fftn, fftkind
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: dx, dy
    real, dimension(nk), intent(in) :: zh1                ! scalar levels (m)
    real, dimension(nk+1), intent(in) :: zf1              ! w levels (m); zf1(1) the ground, zf1(nk+1) the top
    real, dimension(ni,nj,nk), intent(in) :: q            ! C/m3
    real, dimension(ni,nj,nk), intent(out) :: phi         ! V
    complex(fftkind), dimension(:,:,:), allocatable :: a
    complex(fftkind), dimension(ni*nj) :: work
    complex(fftkind), dimension(nk) :: cp, dp
    real(fftkind), dimension(nk) :: lo, up, dg
    real(fftkind) :: lam, pi, den
    integer :: i, j, k

    allocate( a(ni,nj,nk) )
    pi = 4.0_fftkind*atan(1.0_fftkind)
    do k = 1,nk
      if( k .eq. 1 )then
        lo(k) = 1.0_fftkind/((zf1(2) - zf1(1))*2.0_fftkind*(zh1(1) - zf1(1)))
      else
        lo(k) = 1.0_fftkind/((zf1(k+1) - zf1(k))*(zh1(k) - zh1(k-1)))
      endif
      if( k .eq. nk )then
        up(k) = 1.0_fftkind/((zf1(nk+1) - zf1(nk))*2.0_fftkind*(zf1(nk+1) - zh1(nk)))
      else
        up(k) = 1.0_fftkind/((zf1(k+1) - zf1(k))*(zh1(k+1) - zh1(k)))
      endif
      dg(k) = -lo(k) - up(k)
    enddo
    dg(1) = dg(1) - lo(1)                                 ! ghost below the ground: -phi(1)
    dg(nk) = dg(nk) - up(nk)                              ! ghost above the top: -phi(nk)

    !$omp parallel do default(shared) private(k,work)
    do k = 1,nk
      work = reshape( cmplx( -q(:,:,k)/eps0, 0.0, kind=fftkind ), (/ ni*nj /) )
      call fftn( work, (/ ni, nj /) )
      a(:,:,k) = reshape( work, (/ ni, nj /) )
    enddo

    !$omp parallel do default(shared) private(i,j,k,lam,cp,dp,den)
    do j = 1,nj
      do i = 1,ni
        lam = (2.0_fftkind - 2.0_fftkind*cos(2.0_fftkind*pi*(i-1)/ni))/dx**2
        if( nj .gt. 1 ) lam = lam + (2.0_fftkind - 2.0_fftkind*cos(2.0_fftkind*pi*(j-1)/nj))/dy**2
        ! Thomas algorithm for lo(k) a(k-1) + (dg(k) - lam) a(k) + up(k) a(k+1) = rhs(k)
        den = dg(1) - lam
        cp(1) = up(1)/den
        dp(1) = a(i,j,1)/den
        do k = 2,nk
          den = dg(k) - lam - lo(k)*cp(k-1)
          cp(k) = up(k)/den
          dp(k) = (a(i,j,k) - lo(k)*dp(k-1))/den
        enddo
        a(i,j,nk) = dp(nk)
        do k = nk-1,1,-1
          a(i,j,k) = dp(k) - cp(k)*a(i,j,k+1)
        enddo
      enddo
    enddo

    !$omp parallel do default(shared) private(k,work)
    do k = 1,nk
      work = reshape( a(:,:,k), (/ ni*nj /) )
      call fftn( work, (/ ni, nj /), inv=.true. )
      phi(:,:,k) = real( reshape( work, (/ ni, nj /) ) )
    enddo
  end subroutine potential

  subroutine field(ni, nj, nk, dx, dy, zh1, zf1, phi, ex, ey, ez, emag)
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: dx, dy
    real, dimension(nk), intent(in) :: zh1
    real, dimension(nk+1), intent(in) :: zf1
    real, dimension(ni,nj,nk), intent(in) :: phi
    real, dimension(ni,nj,nk), intent(out) :: ex, ey, ez, emag      ! V/m
    integer :: i, j, k, ip, im, jp, jm
    real :: pu, pl, zu, zl
    !$omp parallel do default(shared) private(i,j,k,ip,im,jp,jm,pu,pl,zu,zl)
    do k = 1,nk
      if( k .eq. 1 )then
        zl = 2.0*zf1(1) - zh1(1)
      else
        zl = zh1(k-1)
      endif
      if( k .eq. nk )then
        zu = 2.0*zf1(nk+1) - zh1(nk)
      else
        zu = zh1(k+1)
      endif
      do j = 1,nj
        jp = j + 1
        if( jp .gt. nj ) jp = 1
        jm = j - 1
        if( jm .lt. 1 ) jm = nj
        do i = 1,ni
          ip = i + 1
          if( ip .gt. ni ) ip = 1
          im = i - 1
          if( im .lt. 1 ) im = ni
          ex(i,j,k) = -(phi(ip,j,k) - phi(im,j,k))/(2.0*dx)
          ey(i,j,k) = 0.0
          if( nj .gt. 1 ) ey(i,j,k) = -(phi(i,jp,k) - phi(i,jm,k))/(2.0*dy)
          if( k .eq. 1 )then
            pl = -phi(i,j,1)
          else
            pl = phi(i,j,k-1)
          endif
          if( k .eq. nk )then
            pu = -phi(i,j,nk)
          else
            pu = phi(i,j,k+1)
          endif
          ez(i,j,k) = -(pu - pl)/(zu - zl)
          emag(i,j,k) = sqrt( ex(i,j,k)**2 + ey(i,j,k)**2 + ez(i,j,k)**2 )
        enddo
      enddo
    enddo
  end subroutine field

  double precision function energy(ni, nj, nk, dx, dy, zf1, q, phi)
    ! electrostatic energy (J): half the integral of charge density times potential
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: dx, dy
    real, dimension(nk+1), intent(in) :: zf1
    real, dimension(ni,nj,nk), intent(in) :: q, phi
    integer :: i, j, k
    energy = 0.0d0
    do k = 1,nk
      do j = 1,nj
        do i = 1,ni
          energy = energy + 0.5d0*dble(q(i,j,k))*dble(phi(i,j,k))*dx*dy*(zf1(k+1) - zf1(k))
        enddo
      enddo
    enddo
  end function energy

  subroutine breakdown_field(ni, nj, nk, rho, option, ebrk)
    ! option 1: WRF-ELEC's nssl_ibrkd = 4, 284 kV/m at 1.225 kg/m3 scaled by air density and held within 50-180 kV/m;
    ! option 2: the same scaling without the bounds, the runaway-breakdown field's proportion to air density
    integer, intent(in) :: ni, nj, nk, option
    real, dimension(ni,nj,nk), intent(in) :: rho
    real, dimension(ni,nj,nk), intent(out) :: ebrk
    if( option .eq. 2 )then
      ebrk = 284.0e3*rho/1.225
    else
      ebrk = max( 50.0e3, min( 284.0e3*rho/1.225, 180.0e3 ) )
    endif
  end subroutine breakdown_field

  subroutine discharge(ni, nj, nk, dx, dy, zf1, rho, emag, ebrk, radius, q, dep, ninit, ncol, nregion, qpos, qneg)
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: dx, dy, radius
    real, dimension(nk+1), intent(in) :: zf1
    real, dimension(ni,nj,nk), intent(in) :: rho, emag, ebrk
    real, dimension(ni,nj,nk), intent(inout) :: q        ! C/m3, updated
    real, dimension(ni,nj,nk), intent(out) :: dep        ! C/m3 handed to the small ions
    integer, intent(out) :: ninit, ncol, nregion
    double precision, intent(out) :: qpos, qneg          ! positive and negative charge removed (C)
    integer, dimension(ni,nj) :: col, label
    integer, dimension(2,ni*nj) :: stack
    integer, dimension(4), parameter :: nbi = (/ 1, -1, 0, 0 /), nbj = (/ 0, 0, 1, -1 /)
    integer :: i, j, k, ii, jj, ia, ja, di, dj, nri, nrj, top, n, nb
    double precision :: chgpos, chgneg, totpos, totneg, fracp, fracn, x, y, dv, chg
    real :: qmax

    ! columns within the radius of a point above the breakdown field
    col = 0
    ninit = 0
    nri = int(radius/dx)
    nrj = 0
    if( nj .gt. 1 ) nrj = int(radius/dy)
    do j = 1,nj
      do i = 1,ni
        n = count( emag(i,j,:) .gt. ebrk(i,j,:) )
        ninit = ninit + n
        if( n .gt. 0 )then
          do dj = -nrj,nrj
            do di = -nri,nri
              if( (di*dx)**2 + (dj*dy)**2 .le. radius**2 )then
                col( modulo(i - 1 + di, ni) + 1, modulo(j - 1 + dj, nj) + 1 ) = 1
              endif
            enddo
          enddo
        endif
      enddo
    enddo

    ! the charge available in those columns that hold charge above chgthr
    chgpos = 0.0d0
    chgneg = 0.0d0
    totpos = 0.0d0
    totneg = 0.0d0
    ncol = 0
    do j = 1,nj
      do i = 1,ni
        if( col(i,j) .eq. 1 )then
          qmax = maxval( abs( q(i,j,1:nk-1) ) )
          if( qmax .le. chgthr )then
            col(i,j) = 0
          else
            ncol = ncol + 1
            do k = 1,nk-1
              dv = dble(dx)*dy*(zf1(k+1) - zf1(k))
              if( q(i,j,k) .gt. chgthr )then
                chgpos = chgpos + (q(i,j,k) - chgthr)*dv
              elseif( q(i,j,k) .lt. -chgthr )then
                chgneg = chgneg + (q(i,j,k) + chgthr)*dv
              endif
              totneg = totneg + min(0.0, q(i,j,k))*dv
              totpos = totpos + max(0.0, q(i,j,k))*dv
            enddo
          endif
        endif
      enddo
    enddo

    ! the shares removed (light1d: the newer rule, falling back to the original when either share left is below 0.1)
    if( chgpos .lt. 1.0d0 ) chgpos = 1.0d0
    if( chgneg .gt. -1.0d0 ) chgneg = -1.0d0
    fracn = disfrac*0.5d0*(abs(chgneg) + abs(totneg))
    fracp = disfrac*0.5d0*(chgpos + totpos)
    if( sign(1.0d0, totneg + totpos) .ne. sign(1.0d0, chgneg + totneg + chgpos + totpos) )then
      fracn = min(fracn, fracp)
      fracp = fracn
    endif
    x = 0.0d0
    if( abs(totneg) .gt. 0.0d0 ) x = (chgneg + fracn)/chgneg
    y = 0.0d0
    if( abs(totpos) .gt. 0.0d0 ) y = (chgpos - fracp)/abs(chgpos)
    if( x .lt. 0.1d0 .or. y .lt. 0.1d0 )then
      chg = disfrac*max(chgpos, abs(chgneg))
      chg = min(chg, chgpos, abs(chgneg))
      x = (abs(chgneg) - chg)/abs(chgneg)
      y = (chgpos - chg)/abs(chgpos)
    endif

    dep = 0.0
    qpos = 0.0d0
    qneg = 0.0d0
    do k = 1,nk
      dv = dble(dx)*dy*(zf1(k+1) - zf1(k))
      do j = 1,nj
        do i = 1,ni
          if( col(i,j) .eq. 1 )then
            if( q(i,j,k) .lt. -chgthr )then
              dep(i,j,k) = -(q(i,j,k) + chgthr)*(1.0 - x)
              qneg = qneg + dep(i,j,k)*dv
            elseif( q(i,j,k) .gt. chgthr )then
              dep(i,j,k) = -(q(i,j,k) - chgthr)*(1.0 - y)
              qpos = qpos - dep(i,j,k)*dv
            endif
            q(i,j,k) = q(i,j,k) + dep(i,j,k)
          endif
        enddo
      enddo
    enddo

    ! separate discharge regions: connected sets of discharge columns (sides touching, the domain wrapping)
    label = 0
    nregion = 0
    do j = 1,nj
      do i = 1,ni
        if( col(i,j) .eq. 1 .and. label(i,j) .eq. 0 )then
          nregion = nregion + 1
          label(i,j) = nregion
          top = 1
          stack(:,1) = (/ i, j /)
          do while( top .gt. 0 )
            ii = stack(1,top)
            jj = stack(2,top)
            top = top - 1
            do nb = 1,4
              ia = modulo(ii - 1 + nbi(nb), ni) + 1
              ja = modulo(jj - 1 + nbj(nb), nj) + 1
              if( col(ia,ja) .eq. 1 .and. label(ia,ja) .eq. 0 )then
                label(ia,ja) = nregion
                top = top + 1
                stack(:,top) = (/ ia, ja /)
              endif
            enddo
          enddo
        endif
      enddo
    enddo
  end subroutine discharge

end module terluna_lightning
