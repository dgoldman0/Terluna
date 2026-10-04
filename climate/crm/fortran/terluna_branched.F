!-----------------------------------------------------------------------------------------------------------------------
! Terluna: one flash of WRF-ELEC's branched lightning on a domain that wraps around (climate/crm/cm1_elec.py)
!
! WRF-ELEC's lightmsz (module_discharge_msz.F, MicroTed/wrf4-elec e43041b; MacGorman, Straka and Ziegler 2001) starts a
! flash at a point where the field reaches 0.9 of the breakdown field, runs a channel up and down from it until the
! field falls below 0.8 of breakdown, and spreads each end through the connected cells whose potential differs from the
! channel's by more than its 200 V/m internal field allows. It neutralizes the charge there, handing the opposite charge
! to the small ions; it counts a ground strike when a downward channel reaches below zgrnd (or, with zgrnd <= 0, air
! warmer than tgrnd); and it makes nitrogen oxides along the channel at Wang et al.'s yield per metre for the channel's
! pressure.
!
! Its channels stop at the edge of its arrays, so here each flash's arrays are recentred on the flash's starting point
! across the wrapping domain. The starting point is drawn as lightmsz draws it, uniformly among the points of levels 2
! to nk-1 above 0.9 of breakdown, colder than 275 K and not yet used in the step, with the same Park and Miller generator
! seeded from the step (so a restarted run repeats its flashes); lightmsz is then given that point alone.
!-----------------------------------------------------------------------------------------------------------------------
module terluna_branched
  implicit none
  private
  public :: branched_flash, flash_record, park_miller

  type flash_record
    integer :: kind = 0                     ! -1 no starting point left; 0 none; 1 in cloud; 2 negative and 3 positive
                                            ! charge to ground
    integer :: i0 = 0, j0 = 0, k0 = 0       ! starting point
    real :: e0 = 0.0, ebrk0 = 0.0           ! field and breakdown field there (V/m)
    double precision :: qpos = 0.0d0        ! positive charge neutralized (C)
    double precision :: qneg = 0.0d0        ! negative charge neutralized (C)
    double precision :: nox = 0.0d0         ! nitrogen oxides made (mol)
    integer :: npos = 0, nneg = 0           ! channel points of the positive leader (in negative charge) and of the
                                            ! negative leader (in positive charge): lightmsz's flshp and -flshn
    real :: area = 0.0                      ! area of the columns the channels reach (m2)
    real :: zlo = -1.0, zhi = -1.0          ! lowest and highest level the channels reach (m)
  end type flash_record

  real, parameter :: start_frac = 0.9      ! lightmsz's elgtfdel
  real, parameter :: start_tmax = 275.0    ! lightmsz's starting points are colder than this (K)
  real, parameter :: stop_frac = 0.8       ! WRF-ELEC's nssl_elgtfestop
  real, parameter :: cloud_min = 1.0e-5    ! cloud water and ice (kg/kg) that set the cloud top, as WRF-ELEC's driver

contains

  real function park_miller(seed)
    ! Park and Miller's minimal standard generator, as lightmsz's ran0mpi (Numerical Recipes' ran0)
    integer, intent(inout) :: seed
    integer, parameter :: ia = 16807, im = 2147483647, iq = 127773, ir = 2836, mask = 123459876
    real, parameter :: am = 1.0/im
    integer :: k
    seed = ieor(seed, mask)
    k = seed/iq
    seed = ia*(seed - k*iq) - ir*k
    if( seed .lt. 0 ) seed = seed + im
    park_miller = am*seed
    seed = ieor(seed, mask)
  end function park_miller

  subroutine branched_flash(ni, nj, nk, dx, dy, dz, rho, prs, tk, cloud, phi, ex, ey, ez, emag, q, ebrk, used, seed,  &
                            nstep, dtp, eint, tgrnd, zgrnd, lightrad, iunit, dep, rec)
    use module_discharge_msz, only : lightmsz
    integer, intent(in) :: ni, nj, nk, nstep, iunit
    real, intent(in) :: dx, dy, dtp, eint, tgrnd, zgrnd, lightrad
    real, dimension(nk), intent(in) :: dz                        ! layer depths (m)
    real, dimension(ni,nj,nk), intent(in) :: rho, prs, tk, cloud ! kg/m3, Pa, K, kg/kg
    real, dimension(ni,nj,nk), intent(in) :: phi, ex, ey, ez, emag, q, ebrk   ! V, V/m, C/m3, V/m
    integer, dimension(ni,nj,nk), intent(inout) :: used          ! starting points used in this step
    integer, intent(inout) :: seed
    real, dimension(ni,nj,nk), intent(out) :: dep                ! charge handed to the small ions (C/m3)
    type(flash_record), intent(out) :: rec

    real, allocatable, dimension(:,:,:), save :: sctot, pot, dzz, db, pn, temg, sciona, flshi, flshn, flshp, clnox
    real, allocatable, dimension(:,:,:,:), save :: elec
    real, allocatable, dimension(:,:), save :: light, lightdens, lightdis, fic, ficp, ficn, fcg, fcgp, fcgn
    integer, allocatable, dimension(:,:), save :: lightmask, lightinit, kcldtop
    integer, allocatable, dimension(:,:,:), save :: initmask
    integer :: i, j, k, n, ichoose, i0, j0, k0, si, sj, ir, jr, lgtstp, loccur, flashcount
    double precision :: cnoxtot, etot, dv
    real :: z
    logical, dimension(ni,nj) :: reached

    rec = flash_record()
    dep = 0.0

    ! the starting point
    n = count( emag(:,:,2:nk-1) .ge. ebrk(:,:,2:nk-1)*start_frac .and. tk(:,:,2:nk-1) .lt. start_tmax .and.       &
               used(:,:,2:nk-1) .eq. 0 )
    if( n .eq. 0 )then
      rec%kind = -1
      return
    endif
    ichoose = int( real(n)*park_miller(seed) ) + 1
    ichoose = min( ichoose, n )
    n = 0
    search: do k = 2,nk-1
      do j = 1,nj
        do i = 1,ni
          if( emag(i,j,k) .ge. ebrk(i,j,k)*start_frac .and. tk(i,j,k) .lt. start_tmax .and. used(i,j,k) .eq. 0 )then
            n = n + 1
            if( n .eq. ichoose )then
              i0 = i
              j0 = j
              k0 = k
              exit search
            endif
          endif
        enddo
      enddo
    enddo search
    used(i0,j0,k0) = 1
    rec%i0 = i0
    rec%j0 = j0
    rec%k0 = k0
    rec%e0 = emag(i0,j0,k0)
    rec%ebrk0 = ebrk(i0,j0,k0)

    ! the arrays, recentred on the starting point, with a halo of one column that wraps
    if( .not. allocated(sctot) )then
      allocate( sctot(0:ni+1,nk,0:nj+1), pot(0:ni+1,nk,0:nj+1), dzz(0:ni+1,nk,0:nj+1), db(0:ni+1,nk,0:nj+1),     &
                pn(0:ni+1,nk,0:nj+1), temg(0:ni+1,nk,0:nj+1), sciona(0:ni+1,nk,0:nj+1), flshi(0:ni+1,nk,0:nj+1),  &
                flshn(0:ni+1,nk,0:nj+1), flshp(0:ni+1,nk,0:nj+1), clnox(0:ni+1,nk,0:nj+1), elec(0:ni+1,nk,0:nj+1,4) )
      allocate( light(0:ni+1,0:nj+1), lightdens(0:ni+1,0:nj+1), lightdis(0:ni+1,0:nj+1), fic(0:ni+1,0:nj+1),      &
                ficp(0:ni+1,0:nj+1), ficn(0:ni+1,0:nj+1), fcg(0:ni+1,0:nj+1), fcgp(0:ni+1,0:nj+1),                 &
                fcgn(0:ni+1,0:nj+1), lightmask(0:ni+1,0:nj+1), lightinit(0:ni+1,0:nj+1), kcldtop(0:ni+1,0:nj+1),    &
                initmask(0:ni+1,0:nj+1,nk) )
    endif
    si = max(1, ni/2) - i0
    sj = max(1, nj/2) - j0
    kcldtop = 1
    do jr = 0,nj+1
      j = modulo(jr - 1 - sj, nj) + 1
      do k = 1,nk
        do ir = 0,ni+1
          i = modulo(ir - 1 - si, ni) + 1
          sctot(ir,k,jr) = q(i,j,k)
          pot(ir,k,jr) = phi(i,j,k)
          elec(ir,k,jr,1) = ex(i,j,k)
          elec(ir,k,jr,2) = ey(i,j,k)
          elec(ir,k,jr,3) = ez(i,j,k)
          elec(ir,k,jr,4) = emag(i,j,k)
          dzz(ir,k,jr) = dz(k)
          db(ir,k,jr) = rho(i,j,k)
          pn(ir,k,jr) = prs(i,j,k)
          temg(ir,k,jr) = tk(i,j,k)
          if( cloud(i,j,k) .ge. cloud_min ) kcldtop(ir,jr) = max(k, kcldtop(ir,jr))
        enddo
      enddo
    enddo
    sciona = 0.0
    flshi = 0.0
    flshn = 0.0
    flshp = 0.0
    clnox = 0.0
    light = 0.0
    lightdens = 0.0
    lightdis = 0.0
    fic = 0.0
    ficp = 0.0
    ficn = 0.0
    fcg = 0.0
    fcgp = 0.0
    fcgn = 0.0
    lightmask = 0
    lightinit = 0
    initmask = 1                                              ! every point but the chosen one is closed to lightmsz
    initmask(modulo(i0 - 1 + si, ni) + 1, modulo(j0 - 1 + sj, nj) + 1, k0) = 0

    lgtstp = 0
    loccur = 1
    flashcount = 0
    cnoxtot = 0.0d0
    etot = 0.0d0
    call lightmsz(1, lgtstp, loccur, dtp, nstep, 2, stop_frac, cnoxtot,                                         &
                  0, ni, nj, nk, 1, ni, 1, nj, 1, nk, 1, 1,                                                    &
                  0, ni+1, 0, nj+1, 1, nk, temg,                                                               &
                  dx, dy, lightrad, 120.0e3, 0.3, dzz, sctot, pot,                                             &
                  elec, sciona, 0, 0, 0, 0, 0, 0, 0, 0,                                                        &
                  iunit, db, pn, light, lightdens, etot, lightdis,                                             &
                  fic, ficp, ficn, fcg, fcgp, fcgn,                                                            &
                  flshi, flshn, flshp, lightmask, lightinit, initmask,                                         &
                  eint, tgrnd, zgrnd, .true., 1, 1, flashcount, clnox, kcldtop)

    ! back to the domain: the ions' charge (C/kg there, C/m3 here), the channels and the nitrogen oxides
    reached = .false.
    do jr = 0,nj+1
      j = modulo(jr - 1 - sj, nj) + 1
      do k = 1,nk
        do ir = 0,ni+1
          i = modulo(ir - 1 - si, ni) + 1
          if( ir .ge. 1 .and. ir .le. ni .and. jr .ge. 1 .and. jr .le. nj )then
            dep(i,j,k) = dep(i,j,k) + sciona(ir,k,jr)*rho(i,j,k)
          elseif( sciona(ir,k,jr) .ne. 0.0 )then              ! a halo cell is its column's twin across the edge
            dep(i,j,k) = dep(i,j,k) + sciona(ir,k,jr)*rho(i,j,k)
          endif
          if( ir .ge. 1 .and. ir .le. ni .and. jr .ge. 1 .and. jr .le. nj )then
            if( flshn(ir,k,jr) .ne. 0.0 .or. flshp(ir,k,jr) .ne. 0.0 )then
              reached(i,j) = .true.
              rec%nneg = rec%nneg + nint(abs(flshn(ir,k,jr)))   ! lightmsz counts these down from zero
              rec%npos = rec%npos + nint(flshp(ir,k,jr))
              z = sum(dz(1:k)) - 0.5*dz(k)
              if( rec%zlo .lt. 0.0 .or. z .lt. rec%zlo ) rec%zlo = z
              rec%zhi = max(rec%zhi, z)
            endif
          endif
        enddo
      enddo
    enddo
    do k = 1,nk
      dv = dble(dx)*dy*dz(k)
      do j = 1,nj
        do i = 1,ni
          if( dep(i,j,k) .lt. 0.0 ) rec%qpos = rec%qpos - dep(i,j,k)*dv
          if( dep(i,j,k) .gt. 0.0 ) rec%qneg = rec%qneg + dep(i,j,k)*dv
        enddo
      enddo
    enddo
    rec%area = count(reached)*dx*dy
    rec%nox = cnoxtot
    rec%kind = loccur
    if( loccur .lt. 1 .or. loccur .gt. 3 ) rec%kind = 0
  end subroutine branched_flash

end module terluna_branched
