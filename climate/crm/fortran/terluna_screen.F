!-----------------------------------------------------------------------------------------------------------------------
! Terluna: screening layers at cloud edges (climate/crm/cm1_elec.py), after WRF-ELEC's screen (elec/module_screen.F,
! MicroTed/wrf4-elec e43041b; Ziegler, MacGorman, Dye and Ray 1991), on one process and a domain that wraps around
!
! Clear air conducts better than cloud, so where cloud meets clear air the conduction current the field drives into
! the edge does not match the current leaving it, and charge gathers there. A cloud point (cloud water and ice above
! 0.01 g/kg) is an edge along x, y or z where the cloud carries on to one side (above 0.015 and 0.0225 g/kg at the next
! two points) and clear air lies beyond it on the other (below 0.01 g/kg at the next point and 0.001 g/kg at the two
! after, and in z also at the fourth and fifth). The edge's normal comes from the gradient of cloud content, the field
! normal to it inside and outside from one-sided differences of the potential (zero when the two point opposite
! ways), and the charge the edge gathers in a step of dt is
!   -2 eps0 / ds * (E_air sigma_air - E_cloud sigma_cloud) / (sigma_air + sigma_cloud)
!     * (1 - exp(-dt (sigma_air + sigma_cloud) / (2 eps0)))
! with ds the cube root of the cell's volume. It goes to the small ions, capped so that the net charge there reaches at
! most 0.25 nC/m3 of either sign. The vertical differences take the distances between CM1's scalar levels; WRF-ELEC
! takes the depth of the layer below for all of them, the same on an even grid.
!-----------------------------------------------------------------------------------------------------------------------
module terluna_screen
  implicit none
  private
  public :: screen, earth_air_conductivity

  real, parameter :: qctol = 1.0e-5, qctol2 = 1.0e-6, fin1 = 1.5, fin2 = 2.25
  real, parameter :: epascl = 5.0, sc_thresh = 0.25e-9, eps = 8.8592e-12   ! WRF-ELEC's eperao

contains

  real function earth_air_conductivity(z)
    ! WRF-ELEC's clear-air (polar) conductivity (S/m) from 0 to 10 km in 500-m steps, carried on beyond on its last
    ! step's slope
    real, intent(in) :: z
    real, dimension(21), parameter :: table = (/ 5.5e-14, 6.0e-14, 7.0e-14, 7.8e-14, 8.8e-14, 1.0e-13, 1.1e-13,       &
                                                 1.3e-13, 1.5e-13, 1.7e-13, 1.9e-13, 2.1e-13, 2.3e-13, 2.6e-13,       &
                                                 2.8e-13, 3.1e-13, 3.4e-13, 3.8e-13, 4.2e-13, 4.5e-13, 4.9e-13 /)
    integer :: k
    k = min( max( int(z/500.0) + 1, 1 ), 20 )
    earth_air_conductivity = table(k) + (table(k+1) - table(k))*(z - 500.0*(k-1))/500.0
  end function earth_air_conductivity

  subroutine screen(ni, nj, nk, dx, dy, zh1, zf1, dt, phi, cloud, q, sig_air, sig_cloud, dep, qpos, qneg)
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: dx, dy, dt
    real, dimension(nk), intent(in) :: zh1                      ! scalar levels (m)
    real, dimension(nk+1), intent(in) :: zf1                    ! w levels (m)
    real, dimension(ni,nj,nk), intent(in) :: phi                ! V
    real, dimension(ni,nj,nk), intent(in) :: cloud              ! cloud water and ice (kg/kg)
    real, dimension(ni,nj,nk), intent(in) :: q                  ! net charge (C/m3)
    real, dimension(nk), intent(in) :: sig_air, sig_cloud       ! conductivity by level (S/m)
    real, dimension(ni,nj,nk), intent(out) :: dep               ! charge handed to the small ions (C/m3)
    double precision, intent(out) :: qpos, qneg                 ! positive and negative charge added (C)
    integer :: i, j, k, kp4, kp5, km4, km5
    real :: ascl, dccx, dccy, dccz, dcpx, dcpy, dcpz, dapx, dapy, dapz, dccldn, ecn, ean, dz, ds, t8
    real, allocatable, dimension(:,:,:) :: gx, gy, gz, c

    allocate( gx(ni,nj,nk), gy(ni,nj,nk), gz(ni,nj,nk), c(ni,nj,nk) )
    c = max(cloud, 0.0)
    gx = 0.0
    gy = 0.0
    gz = 0.0
    do k = 2,nk-1
      dz = zf1(k+1) - zf1(k)
      do j = 1,nj
        do i = 1,ni
          gx(i,j,k) = (c(ii(i+1),j,k) - c(ii(i-1),j,k))*0.5/dx
          gy(i,j,k) = (c(i,jj(j+1),k) - c(i,jj(j-1),k))*0.5/dy
          gz(i,j,k) = (c(i,j,k+1) - c(i,j,k-1))*0.5/dz
        enddo
      enddo
    enddo

    dep = 0.0
    qpos = 0.0d0
    qneg = 0.0d0
    !$omp parallel do default(shared) reduction(+:qpos,qneg)                                                         &
    !$omp private(i,j,k,kp4,kp5,km4,km5,ascl,dccx,dccy,dccz,dcpx,dcpy,dcpz,dapx,dapy,dapz,dccldn,ecn,ean, &
    !$omp dz,ds,t8)
    do k = 2,nk-1
      dz = zf1(k+1) - zf1(k)
      ds = (dx*dy*dz)**(1.0/3.0)
      kp4 = min(nk-1, k+4)
      kp5 = min(nk-1, k+5)
      km4 = max(1, k-4)
      km5 = max(1, k-5)
      do j = 1,nj
        do i = 1,ni
          if( c(i,j,k) .le. qctol ) cycle
          ascl = 6.0
          dccx = 0.0
          dccy = 0.0
          dccz = 0.0
          dcpx = 0.0
          dcpy = 0.0
          dcpz = 0.0
          dapx = 0.0
          dapy = 0.0
          dapz = 0.0
          ! x
          if( c(ii(i+1),j,k) .lt. qctol .and. c(ii(i+2),j,k) .lt. qctol2 .and. c(ii(i+3),j,k) .lt. qctol2 .and.     &
              c(ii(i-1),j,k) .gt. qctol*fin1 .and. c(ii(i-2),j,k) .gt. qctol*fin2 )then
            ascl = ascl - 1.0
            dccx = 0.5*(gx(i,j,k) + gx(ii(i-1),j,k))
            dcpx = -(phi(i,j,k) - phi(ii(i-1),j,k))/dx
            dapx = -(phi(ii(i+2),j,k) - phi(ii(i+1),j,k))/dx
          endif
          if( c(ii(i-1),j,k) .lt. qctol .and. c(ii(i-2),j,k) .lt. qctol2 .and. c(ii(i-3),j,k) .lt. qctol2 .and.     &
              c(ii(i+1),j,k) .gt. qctol*fin1 .and. c(ii(i+2),j,k) .gt. qctol*fin2 )then
            ascl = ascl - 1.0
            dccx = 0.5*(gx(ii(i+1),j,k) + gx(i,j,k))
            dcpx = -(phi(ii(i+1),j,k) - phi(i,j,k))/dx
            dapx = -(phi(ii(i-1),j,k) - phi(ii(i-2),j,k))/dx
          endif
          if( sign(1.0, dcpx) .ne. sign(1.0, dapx) )then
            dcpx = 0.0
            dapx = 0.0
          endif
          ! y
          if( nj .gt. 1 )then
            if( c(i,jj(j+1),k) .lt. qctol .and. c(i,jj(j+2),k) .lt. qctol2 .and. c(i,jj(j+3),k) .lt. qctol2 .and.   &
                c(i,jj(j-1),k) .gt. qctol*fin1 .and. c(i,jj(j-2),k) .gt. qctol*fin2 )then
              ascl = ascl - 1.0
              dccy = 0.5*(gy(i,j,k) + gy(i,jj(j-1),k))
              dcpy = -(phi(i,j,k) - phi(i,jj(j-1),k))/dy
              dapy = -(phi(i,jj(j+2),k) - phi(i,jj(j+1),k))/dy
            endif
            if( c(i,jj(j-1),k) .lt. qctol .and. c(i,jj(j-2),k) .lt. qctol2 .and. c(i,jj(j-3),k) .lt. qctol2 .and.   &
                c(i,jj(j+1),k) .gt. qctol*fin1 .and. c(i,jj(j+2),k) .gt. qctol*fin2 )then
              ascl = ascl - 1.0
              dccy = 0.5*(gy(i,jj(j+1),k) + gy(i,j,k))
              dcpy = -(phi(i,jj(j+1),k) - phi(i,j,k))/dy
              dapy = -(phi(i,jj(j-1),k) - phi(i,jj(j-2),k))/dy
            endif
            if( sign(1.0, dcpy) .ne. sign(1.0, dapy) )then
              dcpy = 0.0
              dapy = 0.0
            endif
          endif
          ! z
          if( k .gt. 3 .and. k .lt. nk-3 )then
            if( c(i,j,k+1) .lt. qctol .and. c(i,j,k+2) .lt. qctol2 .and. c(i,j,k+3) .lt. qctol2 .and.               &
                c(i,j,kp4) .lt. qctol2 .and. c(i,j,kp5) .lt. qctol2 .and.                                         &
                c(i,j,k-1) .gt. qctol*fin1 .and. c(i,j,k-2) .gt. qctol*fin2 )then
              ascl = ascl - 1.0
              dccz = 0.5*(gz(i,j,k) + gz(i,j,k-1))
              dcpz = -(phi(i,j,k) - phi(i,j,k-1))/(zh1(k) - zh1(k-1))
              dapz = -(phi(i,j,k+2) - phi(i,j,k+1))/(zh1(k+2) - zh1(k+1))
            endif
            if( c(i,j,k-1) .lt. qctol .and. c(i,j,k-2) .lt. qctol2 .and. c(i,j,k-3) .lt. qctol2 .and.               &
                c(i,j,km4) .lt. qctol2 .and. c(i,j,km5) .lt. qctol2 .and.                                         &
                c(i,j,k+1) .gt. qctol*fin1 .and. c(i,j,k+2) .gt. qctol*fin2 )then
              ascl = ascl - 1.0
              dccz = 0.5*(gz(i,j,k+1) + gz(i,j,k))
              dcpz = -(phi(i,j,k+1) - phi(i,j,k))/(zh1(k+1) - zh1(k))
              dapz = -(phi(i,j,k-1) - phi(i,j,k-2))/(zh1(k-1) - zh1(k-2))
            endif
            if( sign(1.0, dcpz) .ne. sign(1.0, dapz) )then
              dcpz = 0.0
              dapz = 0.0
            endif
          endif
          if( ascl .gt. epascl ) cycle                          ! not an edge in any direction
          dccldn = sqrt(dccx**2 + dccy**2 + dccz**2)
          if( dccldn .lt. 1.0e-10 ) cycle
          ecn = -(dccx*dcpx + dccy*dcpy + dccz*dcpz)/dccldn      ! the field along the outward normal, in cloud
          ean = -(dccx*dapx + dccy*dapy + dccz*dapz)/dccldn      ! and in the clear air beyond
          t8 = -2.0*eps/ds*(ean*sig_air(k) - ecn*sig_cloud(k))/(sig_air(k) + sig_cloud(k))                          &
               *(1.0 - exp(-dt*(sig_air(k) + sig_cloud(k))/(2.0*eps)))
          ! the net charge there reaches at most 0.25 nC/m3 of either sign
          if( t8 .gt. 0.0 )then
            if( q(i,j,k) .gt. sc_thresh )then
              t8 = 0.0
            elseif( q(i,j,k) + t8 .gt. sc_thresh )then
              t8 = sc_thresh - q(i,j,k)
            endif
          elseif( t8 .lt. 0.0 )then
            if( q(i,j,k) .lt. -sc_thresh )then
              t8 = 0.0
            elseif( q(i,j,k) + t8 .lt. -sc_thresh )then
              t8 = -sc_thresh - q(i,j,k)
            endif
          endif
          dep(i,j,k) = t8
          if( t8 .gt. 0.0 ) qpos = qpos + dble(t8)*dx*dy*dz
          if( t8 .lt. 0.0 ) qneg = qneg + dble(t8)*dx*dy*dz
        enddo
      enddo
    enddo
    deallocate( gx, gy, gz, c )

  contains

    integer function ii(i)
      integer, intent(in) :: i
      ii = modulo(i - 1, ni) + 1
    end function ii

    integer function jj(j)
      integer, intent(in) :: j
      jj = modulo(j - 1, nj) + 1
    end function jj

  end subroutine screen

end module terluna_screen
