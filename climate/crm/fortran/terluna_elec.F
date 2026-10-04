!-----------------------------------------------------------------------------------------------------------------------
! Terluna: the electricity of CM1 storms (climate/crm/cm1_elec.py; research/studies/atmospheric_electricity)
!
! terluna_nssl_elec takes CM1's calls of its NSSL driver and passes them to WRF-ELEC's (MicroTed/wrf4-elec e43041b;
! Mansell et al. 2005, 2010; Fierro et al. 2013), whose arrays run (i,k,j) where CM1's run (i,j,k). It carries the charge
! of cloud water, rain, cloud ice, snow, graupel, hail and the small ions' net charge in CM1's passive tracers 1-6 (7 with
! hail: scw, scr, sci, scs, sch, sciona, schl; C/kg), which CM1 advects without its positivity limiter (pdtra = 0).
!
! terluna_elec_step follows the microphysics each time step with what WRF-ELEC does outside it: the electric field of the
! net charge (terluna_lightning), lightning, and, where asked, leakage through the air's conductivity. Its vertical
! component is kept for the next step's inductive charging.
!
! Settings from CM1's namelist (param8):
!   var6  ipelec: 0 off, 2 non-inductive charging, 3 with inductive charging
!   var7  isaund, the non-inductive charging law: 12 Saunders and Peck with Brooks's critical rime accretion rate, as
!         WRF-ELEC's default; 1 Takahashi (with nonigrd = -1 in &nssl_mp_params and takahashi.txt in the run folder)
!   var8  lightning: 0 none; 1 WRF-ELEC's cylinders with its breakdown field (nssl_ibrkd = 4); 2 the same with the
!         breakdown field 284 kV/m times air density over 1.225 kg/m3, unbounded
!   var9  leakage: 0 none; 1 net charge relaxes at the rate sigma/eps0, sigma by height from terluna_conductivity.txt
!         (cloud where cloud water and ice reach 0.01 g/kg), its opposite going to the small ions
!   var10 discharge radius (m); 12 km if zero
! Each discharge is logged to terluna_flashes.txt, and every ten steps terluna_field.txt takes the largest field, the
! domain's charge and electrostatic energy and the largest charging rates.
!-----------------------------------------------------------------------------------------------------------------------
module terluna_elec_module
  implicit none
  private
  public :: terluna_nssl_elec, terluna_elec_step

  real, allocatable, save :: ezsave(:,:,:)        ! vertical field (V/m) at scalar points, WRF's (i,k,j) order
  real, allocatable, save :: winduc(:,:,:), wnoninduc(:,:,:)   ! the step's charging rates (C/m3/s), WRF's order
  integer, save :: nsig = 0
  real, allocatable, save :: zsig(:), sigclear(:), sigcloud(:)
  logical, save :: sigread = .false.
  integer, parameter :: maxdisch = 20             ! discharges per time step at most

contains

!-----------------------------------------------------------------------------------------------------------------------
! WRF-ELEC's NSSL driver with CM1's arguments; arrays not passed by CM1 (or not in use) stay absent
!-----------------------------------------------------------------------------------------------------------------------
  subroutine terluna_nssl_elec(th, qv, qc, qr, qi, qs, qh, qhl, cn, cna, cni, ccw, crw, cci, csw, chw, chl,            &
                               zrw, zhw, zhl, vhw, vhl, pii, p, w, dn, dz, dtp, itimestep, RAINNC, RAINNCV, dbz,       &
                               dx, dy, re_cloud, re_ice, re_snow, has_reqc, has_reqi, has_reqs, diagflag,               &
                               ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte,  &
                               pt3d)
    use input, only : ibp, iep, jbp, jep, kbp, kep, npt, var6, nssl_3moment
    use index_module, only : lhl, lvh, lvhl, lccna, lcina
    use module_mp_nssl_2mom, only : nssl_2mom_driver
    integer, intent(in) :: ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte
    integer, intent(in) :: itimestep
    real, intent(in) :: dtp
    real, intent(in), optional :: dx, dy
    integer, intent(in), optional :: has_reqc, has_reqi, has_reqs
    logical, intent(in), optional :: diagflag
    real, dimension(ims:ime,jms:jme,kms:kme), intent(inout) :: qv, qc, qr, qs, qh
    real, dimension(ims:ime,jms:jme,kms:kme), intent(inout), optional :: th, qi, qhl, cn, cna, cni, ccw, crw, cci,    &
                                                   csw, chw, chl, zrw, zhw, zhl, vhw, vhl, dbz, re_cloud, re_ice, re_snow
    real, dimension(ims:ime,jms:jme,kms:kme), intent(in) :: pii, p, w, dn, dz
    real, dimension(ims:ime,jms:jme) :: RAINNC, RAINNCV
    real, dimension(ibp:iep,jbp:jep,kbp:kep,npt), intent(inout) :: pt3d

    real, allocatable, dimension(:,:,:), save :: wth, wqv, wqc, wqr, wqi, wqs, wqh, wqhl, wcn, wccw, wcrw, wcci,     &
                                                 wcsw, wchw, wchl, wvhw, wvhl, wpii, wp, ww, wdn, wdz, wdbz, wrec,     &
                                                 wrei, wres, wscw, wscr, wsci, wscs, wsch, wsciona, wschl, wsctot
    integer :: ipelec, nchg
    logical :: elec, hail, diag
    logical, save :: checked = .false.

    ipelec = nint(var6)
    elec = ipelec .gt. 0
    hail = present(qhl) .and. lhl .gt. 1
    nchg = 6
    if( hail ) nchg = 7
    diag = .false.
    if( present(diagflag) ) diag = diagflag
    if( .not. checked )then
      checked = .true.
      if( nssl_3moment ) stop 'Terluna: the electrified NSSL build carries no three-moment arrays'
      if( lccna .gt. 1 .or. lcina .gt. 1 ) stop 'Terluna: the electrified NSSL build carries no activated CCN or IN'
      if( elec .and. .not. present(ccw) ) stop 'Terluna: WRF-ELEC charges only the two-moment NSSL scheme'
      if( elec .and. npt .lt. nchg ) stop 'Terluna: electrified NSSL needs npt = 6 passive tracers (7 with hail)'
      if( elec .and. ( ibp .ne. ims .or. iep .ne. ime .or. jbp .ne. jms .or. jep .ne. jme .or. kbp .ne. kms .or.    &
                       kep .ne. kme ) ) stop 'Terluna: the passive tracers must span the model arrays'
      print *,'Terluna: WRF-ELEC NSSL driver, ipelec = ',ipelec,', charge tracers = ',nchg
    endif

    call load(th, wth, .true.)
    call load(qv, wqv, .true.)
    call load(qc, wqc, .true.)
    call load(qr, wqr, .true.)
    call load(qi, wqi, .true.)
    call load(qs, wqs, .true.)
    call load(qh, wqh, .true.)
    call load(qhl, wqhl, hail)
    call load(cn, wcn, .true.)
    call load(ccw, wccw, .true.)
    call load(crw, wcrw, .true.)
    call load(cci, wcci, .true.)
    call load(csw, wcsw, .true.)
    call load(chw, wchw, .true.)
    call load(chl, wchl, hail)
    call load(vhw, wvhw, lvh .gt. 1)                      ! CM1 passes q3d(...,lvh) even when lvh = 0
    call load(vhl, wvhl, hail .and. lvhl .gt. 1)
    call load(pii, wpii, .true.)
    call load(p, wp, .true.)
    call load(w, ww, .true.)
    call load(dn, wdn, .true.)
    call load(dz, wdz, .true.)
    call load(dbz, wdbz, diag)
    call load(re_cloud, wrec, flag(has_reqc))
    call load(re_ice, wrei, flag(has_reqi))
    call load(re_snow, wres, flag(has_reqs))
    call tracer(1, wscw, elec)
    call tracer(2, wscr, elec)
    call tracer(3, wsci, elec)
    call tracer(4, wscs, elec)
    call tracer(5, wsch, elec)
    call tracer(6, wsciona, elec)
    call tracer(7, wschl, elec .and. hail)
    if( elec )then
      ! the charging rates (C/m3/s) are only written here; given them, the driver also prints its charging totals
      if( .not. allocated(wsctot) )then
        allocate( wsctot(ims:ime,kms:kme,jms:jme), winduc(ims:ime,kms:kme,jms:jme), wnoninduc(ims:ime,kms:kme,jms:jme) )
        winduc = 0.0                                      ! the driver writes the interior only
        wnoninduc = 0.0
      endif
      if( .not. allocated(ezsave) )then
        allocate( ezsave(ims:ime,kms:kme,jms:jme) )
        ezsave = 0.0
      endif
    endif

    call nssl_2mom_driver(th=wth, qv=wqv, qc=wqc, qr=wqr, qi=wqi, qs=wqs, qh=wqh, qhl=wqhl, cn=wcn,                    &
                          ccw=wccw, crw=wcrw, cci=wcci, csw=wcsw, chw=wchw, chl=wchl, vhw=wvhw, vhl=wvhl,                &
                          pii=wpii, p=wp, w=ww, dn=wdn, dz=wdz, dtp=dtp, itimestep=itimestep,                           &
                          RAINNC=RAINNC, RAINNCV=RAINNCV, dbz=wdbz, dx=dx, dy=dy,                                       &
                          re_cloud=wrec, re_ice=wrei, re_snow=wres, has_reqc=has_reqc, has_reqi=has_reqi,               &
                          has_reqs=has_reqs, diagflag=diagflag,                                                         &
                          scw=wscw, scr=wscr, sci=wsci, scs=wscs, sch=wsch, schl=wschl, sciona=wsciona, sctot=wsctot,   &
                          f_sciona=elec, f_scion=.false., elecz=ezsave, ipelectmp=ipelec, induc=winduc,                &
                          noninduc=wnoninduc,                                                                           &
                          ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme,     &
                          kms=kms, kme=kme, its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)

    call store(th, wth)
    call store(qv, wqv)
    call store(qc, wqc)
    call store(qr, wqr)
    call store(qi, wqi)
    call store(qs, wqs)
    call store(qh, wqh)
    call store(qhl, wqhl)
    call store(cn, wcn)
    call store(ccw, wccw)
    call store(crw, wcrw)
    call store(cci, wcci)
    call store(csw, wcsw)
    call store(chw, wchw)
    call store(chl, wchl)
    call store(vhw, wvhw)
    call store(vhl, wvhl)
    call store(dbz, wdbz)
    call store(re_cloud, wrec)
    call store(re_ice, wrei)
    call store(re_snow, wres)
    if( elec )then
      call untracer(1, wscw)
      call untracer(2, wscr)
      call untracer(3, wsci)
      call untracer(4, wscs)
      call untracer(5, wsch)
      call untracer(6, wsciona)
      if( hail ) call untracer(7, wschl)
    endif

  contains

    logical function flag(has)
      integer, intent(in), optional :: has
      flag = .false.
      if( present(has) ) flag = has .eq. 1
    end function flag

    subroutine load(src, dst, use)
      real, dimension(ims:ime,jms:jme,kms:kme), intent(in), optional :: src
      real, allocatable, dimension(:,:,:), intent(inout) :: dst
      logical, intent(in) :: use
      integer :: i, j, k
      if( present(src) .and. use )then
        if( .not. allocated(dst) ) allocate( dst(ims:ime,kms:kme,jms:jme) )
        !$omp parallel do default(shared) private(i,j,k)
        do j = jms,jme
          do k = kms,kme
            do i = ims,ime
              dst(i,k,j) = src(i,j,k)
            enddo
          enddo
        enddo
      elseif( allocated(dst) )then
        deallocate( dst )
      endif
    end subroutine load

    subroutine store(dst, src)
      real, dimension(ims:ime,jms:jme,kms:kme), intent(inout), optional :: dst
      real, allocatable, dimension(:,:,:), intent(in) :: src
      integer :: i, j, k
      if( present(dst) .and. allocated(src) )then
        !$omp parallel do default(shared) private(i,j,k)
        do j = jts,jte
          do k = kts,kte
            do i = its,ite
              dst(i,j,k) = src(i,k,j)
            enddo
          enddo
        enddo
      endif
    end subroutine store

    subroutine tracer(n, dst, use)
      integer, intent(in) :: n
      real, allocatable, dimension(:,:,:), intent(inout) :: dst
      logical, intent(in) :: use
      integer :: i, j, k
      if( use )then
        if( .not. allocated(dst) ) allocate( dst(ims:ime,kms:kme,jms:jme) )
        !$omp parallel do default(shared) private(i,j,k)
        do j = jms,jme
          do k = kms,kme
            do i = ims,ime
              dst(i,k,j) = pt3d(i,j,k,n)
            enddo
          enddo
        enddo
      elseif( allocated(dst) )then
        deallocate( dst )
      endif
    end subroutine tracer

    subroutine untracer(n, src)
      integer, intent(in) :: n
      real, allocatable, dimension(:,:,:), intent(in) :: src
      integer :: i, j, k
      !$omp parallel do default(shared) private(i,j,k)
      do j = jts,jte
        do k = kts,kte
          do i = its,ite
            pt3d(i,j,k,n) = src(i,k,j)
          enddo
        enddo
      enddo
    end subroutine untracer

  end subroutine terluna_nssl_elec

!-----------------------------------------------------------------------------------------------------------------------
! The field, lightning and leakage after the microphysics
!-----------------------------------------------------------------------------------------------------------------------
  subroutine terluna_elec_step(nstep, mtime, dt, zh, zf, rho, q3d, pt3d)
    use input, only : ni, nj, nk, ib, ie, jb, je, kb, ke, ibm, iem, jbm, jem, kbm, kem, numq, ibp, iep, jbp, jep,     &
                      kbp, kep, npt, dx, dy, var6, var8, var9, var10, terrain_flag
    use index_module, only : lc, li, lhl
    use terluna_lightning, only : potential, field, energy, breakdown_field, discharge, eps0
    integer, intent(in) :: nstep
    double precision, intent(in) :: mtime
    real, intent(in) :: dt
    real, dimension(ib:ie,jb:je,kb:ke), intent(in) :: zh, rho
    real, dimension(ib:ie,jb:je,kb:ke+1), intent(in) :: zf
    real, dimension(ibm:iem,jbm:jem,kbm:kem,numq), intent(in) :: q3d
    real, dimension(ibp:iep,jbp:jep,kbp:kep,npt), intent(inout) :: pt3d

    real, allocatable, dimension(:,:,:), save :: q, phi, ex, ey, ez, emag, ebrk, dep, rhoc
    real, dimension(nk) :: zh1
    real, dimension(nk+1) :: zf1
    integer :: i, j, k, n, nchg, ilight, ileak, iter, ninit, ncol, nregion, loc(3), flashes
    real :: radius, sig, frac, time, e0, rates(2)
    double precision :: w0, w1, qpos, qneg, totpos, totneg
    logical, save :: started = .false.

    if( nint(var6) .le. 0 ) return
    if( .not. started )then                              ! a run's first step: mark where it takes the logs up
      started = .true.
      call mark_start(real(mtime))
    endif
    if( terrain_flag ) stop 'Terluna: the electric field solver takes flat ground'
    ilight = nint(var8)
    ileak = nint(var9)
    radius = var10
    if( radius .le. 0.0 ) radius = 12000.0
    nchg = 6
    if( lhl .gt. 1 ) nchg = 7
    time = real(mtime) + dt
    if( .not. allocated(q) )then
      allocate( q(ni,nj,nk), phi(ni,nj,nk), ex(ni,nj,nk), ey(ni,nj,nk), ez(ni,nj,nk), emag(ni,nj,nk) )
      allocate( ebrk(ni,nj,nk), dep(ni,nj,nk), rhoc(ni,nj,nk) )
    endif
    if( .not. allocated(ezsave) )then
      allocate( ezsave(ib:ie,kb:ke,jb:je) )
      ezsave = 0.0
    endif
    do k = 1,nk
      zh1(k) = zh(1,1,k)
      zf1(k) = zf(1,1,k)
    enddo
    zf1(nk+1) = zf(1,1,nk+1)

    !$omp parallel do default(shared) private(i,j,k,n)
    do k = 1,nk
      do j = 1,nj
        do i = 1,ni
          rhoc(i,j,k) = rho(i,j,k)
          q(i,j,k) = 0.0
          do n = 1,nchg
            q(i,j,k) = q(i,j,k) + pt3d(i,j,k,n)
          enddo
          q(i,j,k) = q(i,j,k)*rho(i,j,k)
        enddo
      enddo
    enddo

    if( ileak .eq. 1 )then
      call read_conductivity()
      !$omp parallel do default(shared) private(i,j,k,sig,frac)
      do k = 1,nk
        do j = 1,nj
          do i = 1,ni
            if( q3d(i,j,k,lc) + q3d(i,j,k,li) .ge. 1.0e-5 )then
              sig = interp(zh1(k), sigcloud)
            else
              sig = interp(zh1(k), sigclear)
            endif
            frac = 1.0 - exp( -dt*sig/eps0 )
            pt3d(i,j,k,6) = pt3d(i,j,k,6) - frac*q(i,j,k)/rho(i,j,k)
            q(i,j,k) = q(i,j,k)*(1.0 - frac)
          enddo
        enddo
      enddo
    endif

    call potential(ni, nj, nk, dx, dy, zh1, zf1, q, phi)
    call field(ni, nj, nk, dx, dy, zh1, zf1, phi, ex, ey, ez, emag)

    flashes = 0
    if( ilight .ge. 1 )then
      call breakdown_field(ni, nj, nk, rhoc, ilight, ebrk)
      do iter = 1,maxdisch
        if( .not. any( emag .gt. ebrk ) ) exit
        loc = maxloc( emag/ebrk )
        e0 = emag(loc(1),loc(2),loc(3))
        w0 = energy(ni, nj, nk, dx, dy, zf1, q, phi)
        call discharge(ni, nj, nk, dx, dy, zf1, rhoc, emag, ebrk, radius, q, dep, ninit, ncol, nregion, qpos, qneg)
        if( ncol .eq. 0 ) exit                            ! above breakdown where no column holds charge to remove
        !$omp parallel do default(shared) private(i,j,k)
        do k = 1,nk
          do j = 1,nj
            do i = 1,ni
              pt3d(i,j,k,6) = pt3d(i,j,k,6) + dep(i,j,k)/rho(i,j,k)
            enddo
          enddo
        enddo
        call potential(ni, nj, nk, dx, dy, zh1, zf1, q, phi)
        call field(ni, nj, nk, dx, dy, zh1, zf1, phi, ex, ey, ez, emag)
        w1 = energy(ni, nj, nk, dx, dy, zf1, q, phi)
        call log_flash(time, nstep, iter, ninit, e0, ebrk(loc(1),loc(2),loc(3)),                                 &
                       (loc(1)-0.5)*dx, (loc(2)-0.5)*dy, zh1(loc(3)), ncol, nregion, qpos, qneg, w0, w1,           &
                       ni, nj, nk, dep, zh1)
        flashes = flashes + 1
      enddo
      if( flashes .eq. maxdisch ) print *,'Terluna: ',maxdisch,' discharges in step ',nstep
    endif

    !$omp parallel do default(shared) private(i,j,k)
    do j = 1,nj
      do k = 1,nk
        do i = 1,ni
          ezsave(i,k,j) = ez(i,j,k)
        enddo
      enddo
    enddo

    if( mod(nstep, 10) .eq. 0 .or. flashes .gt. 0 )then
      totpos = 0.0d0
      totneg = 0.0d0
      do k = 1,nk
        do j = 1,nj
          do i = 1,ni
            if( q(i,j,k) .gt. 0.0 )then
              totpos = totpos + dble(q(i,j,k))*dx*dy*(zf1(k+1) - zf1(k))
            else
              totneg = totneg + dble(q(i,j,k))*dx*dy*(zf1(k+1) - zf1(k))
            endif
          enddo
        enddo
      enddo
      loc = maxloc( emag )
      rates = 0.0
      if( allocated(wnoninduc) ) rates = (/ maxval(abs(wnoninduc(1:ni,1:nk,1:nj))), maxval(abs(winduc(1:ni,1:nk,1:nj))) /)
      call log_field(time, nstep, emag(loc(1),loc(2),loc(3)), (loc(1)-0.5)*dx, (loc(2)-0.5)*dy, zh1(loc(3)),          &
                     totpos, totneg, energy(ni, nj, nk, dx, dy, zf1, q, phi), flashes, rates)
    endif

  contains

    real function interp(z, vs)
      ! log-linear in the conductivity, held at the ends of the table
      real, intent(in) :: z
      real, dimension(nsig), intent(in) :: vs
      integer :: m
      interp = vs(nsig)
      if( z .le. zsig(1) )then
        interp = vs(1)
        return
      endif
      do m = 1,nsig-1
        if( z .le. zsig(m+1) )then
          interp = exp( log(vs(m)) + (log(vs(m+1)) - log(vs(m)))*(z - zsig(m))/(zsig(m+1) - zsig(m)) )
          return
        endif
      enddo
    end function interp

  end subroutine terluna_elec_step

  subroutine log_flash(time, nstep, iter, ninit, emax, ebrk, x, y, z, ncol, nregion, qpos, qneg, w0, w1,             &
                       ni, nj, nk, dep, zh1)
    real, intent(in) :: time, emax, ebrk, x, y, z
    integer, intent(in) :: nstep, iter, ninit, ncol, nregion, ni, nj, nk
    double precision, intent(in) :: qpos, qneg, w0, w1
    real, dimension(ni,nj,nk), intent(in) :: dep
    real, dimension(nk), intent(in) :: zh1
    real :: zlo, zhi
    integer :: k, unit
    zlo = -1.0
    zhi = -1.0
    do k = 1,nk
      if( any( dep(:,:,k) .ne. 0.0 ) )then
        if( zlo .lt. 0.0 ) zlo = zh1(k)
        zhi = zh1(k)
      endif
    enddo
    call open_log('terluna_flashes.txt', unit)
    write(unit,'(f13.1,i9,i4,i8,2es11.3,3f11.1,2i6,4es12.4,2f11.1)') time, nstep, iter, ninit, emax, ebrk, x, y, z,   &
          ncol, nregion, qpos, qneg, w0, w1, zlo, zhi
    close(unit)
  end subroutine log_flash

  subroutine log_field(time, nstep, emax, x, y, z, totpos, totneg, w, flashes, rates)
    real, intent(in) :: time, emax, x, y, z, rates(2)
    integer, intent(in) :: nstep, flashes
    double precision, intent(in) :: totpos, totneg, w
    integer :: unit
    call open_log('terluna_field.txt', unit)
    write(unit,'(f13.1,i9,es11.3,3f11.1,3es12.4,i4,2es11.3)') time, nstep, emax, x, y, z, totpos, totneg, w, flashes, &
          rates
    close(unit)
  end subroutine log_field

  subroutine open_log(name, unit)
    ! a log opened for appending, its column names written when it is new
    character(len=*), intent(in) :: name
    integer, intent(out) :: unit
    logical :: exists
    inquire( file=name, exist=exists )
    open( newunit=unit, file=name, status='unknown', position='append', action='write' )
    if( exists ) return
    if( name .eq. 'terluna_flashes.txt' )then
      write(unit,'(a)') '# time_s step discharge points_above e_max_V_m e_break_V_m x_m y_m z_m columns regions '//  &
                        'positive_C negative_C energy_before_J energy_after_J z_low_m z_high_m'
    else
      write(unit,'(a)') '# time_s step e_max_V_m x_m y_m z_m positive_C negative_C energy_J discharges '//         &
                        'noninductive_max_C_m3_s inductive_max_C_m3_s'
    endif
  end subroutine open_log

  subroutine mark_start(time)
    ! rows a restart replaces (those after the time the run takes up from) are dropped when the logs are read
    real, intent(in) :: time
    integer :: unit
    call open_log('terluna_flashes.txt', unit)
    write(unit,'(a,f13.1)') '# run from ', time
    close(unit)
    call open_log('terluna_field.txt', unit)
    write(unit,'(a,f13.1)') '# run from ', time
    close(unit)
  end subroutine mark_start

  subroutine read_conductivity()
    ! terluna_conductivity.txt: a count, then height (m), clear-air and cloud conductivity (S/m) by line
    integer :: unit, n
    logical :: exists
    if( sigread ) return
    sigread = .true.
    inquire( file='terluna_conductivity.txt', exist=exists )
    if( .not. exists ) stop 'Terluna: leakage needs terluna_conductivity.txt'
    open( newunit=unit, file='terluna_conductivity.txt', status='old', action='read' )
    read(unit,*) nsig
    allocate( zsig(nsig), sigclear(nsig), sigcloud(nsig) )
    do n = 1,nsig
      read(unit,*) zsig(n), sigclear(n), sigcloud(n)
    enddo
    close(unit)
    print *,'Terluna: conductivity at ',nsig,' heights; clear air ',sigclear(1),' to ',sigclear(nsig),' S/m'
  end subroutine read_conductivity

end module terluna_elec_module
