!-----------------------------------------------------------------------------------------------------------------------
! Terluna: the electricity of CM1 storms (climate/crm/cm1_elec.py; research/studies/atmospheric_electricity)
!
! terluna_nssl_elec takes CM1's calls of its NSSL driver and passes them to WRF-ELEC's (MicroTed/wrf4-elec e43041b;
! Mansell et al. 2005, 2010; Fierro et al. 2013), whose arrays run (i,k,j) where CM1's run (i,j,k). It carries the charge
! of cloud water, rain, cloud ice, snow, graupel, hail and the small ions' net charge in CM1's passive tracers 1-6 (7 with
! hail: scw, scr, sci, scs, sch, sciona, schl; C/kg), which CM1 advects without its positivity limiter (pdtra = 0).
!
! With charging on it follows WRF-ELEC's own driver (include_microphysics_driver_elec.F): the microphysics runs without
! sedimentation, and the step's sedimentation then runs in sub-steps (WRF-ELEC's 0.75 s, or var5), each followed by
! the electric field of the net charge (terluna_lightning) and lightning, so a discharge acts as the field reaches
! breakdown. Screening layers, leakage and point discharge from the ground, where asked, follow once a step, and the
! last field's vertical component is kept for the next step's inductive charging.
!
! Settings from CM1's namelist (param8):
!   var2  point discharge from the ground (WRF-ELEC has none): its onset field (V/m) at 1.225 kg/m3, scaled by the
!         air's density at the ground; 0 none. Standler and Winn (1979) measured onsets of 3 kV/m over dense vegetation
!         and 5 kV/m on a barren ridge
!   var3  screening layers at cloud edges (WRF-ELEC's iscreen, terluna_screen.F): 0 none; 1 with WRF-ELEC's Earth
!         clear-air conductivity and a tenth of it in cloud; 2 with both by height from terluna_conductivity.txt
!   var4  for branched lightning, the height (m) below which a downward channel counts as reaching the ground; zero or
!         less for WRF-ELEC's rule, air warmer than 266.16 K
!   var5  the sub-step (s) for sedimentation and lightning; WRF-ELEC's 0.75 s if zero
!   var6  ipelec: 0 off, 2 non-inductive charging, 3 with inductive charging
!   var7  isaund, the non-inductive charging law: 12 Saunders and Peck with Brooks's critical rime accretion rate, as
!         WRF-ELEC's default; 11 the same law with a smooth critical rate and charging below -32.5 C, as WRF-ELEC's
!         test case; 1 Takahashi (with nonigrd = -1 in &nssl_mp_params and takahashi.txt in the run folder)
!   var8  lightning: 0 none; 1 WRF-ELEC's cylinders (light1d) and 3 its branched scheme (lightmsz), with its breakdown
!         field (nssl_ibrkd = 4: 284 kV/m times air density over 1.225 kg/m3, within 50-180 kV/m); 2 and 4 the same
!         without the 180-kV/m cap, which in lunar air applies below about 36 km (the 50-kV/m floor stays: above
!         about 100 km in lunar air, 14 km in Earth's)
!   var9  leakage: 0 none; 1 net charge relaxes at the rate sigma/eps0, sigma by height from terluna_conductivity.txt
!         (cloud where cloud water and ice reach 0.01 g/kg), its opposite going to the small ions
!   var10 cylinder radius (m); 12 km if zero
! Each discharge or flash is logged to terluna_flashes.txt and lightmsz's own report on the first 2000 flashes of a
! run to terluna_msz.log; every ten steps terluna_field.txt takes the largest field, the domain's charge and
! electrostatic energy and the largest charging rates, and terluna_ground.txt the field at the ground and the charge
! the ground's points give off.
! The cloud ice size the radiation reads is capped at 140 microns, where RRTMG's ice optics end, and
! terluna_ice_optics.txt logs the ice the cap cuts.
!-----------------------------------------------------------------------------------------------------------------------
module terluna_elec_module
  implicit none
  private
  public :: terluna_nssl_elec, terluna_elec_restart_write

  real, allocatable, save :: ezsave(:,:,:)        ! vertical field (V/m) at scalar points, WRF's (i,k,j) order
  real, allocatable, save :: winduc(:,:,:), wnoninduc(:,:,:)   ! the step's charging rates (C/m3/s), WRF's order
  integer, save :: nsig = 0
  real, allocatable, save :: zsig(:), sigclear(:), sigcloud(:)
  logical, save :: sigread = .false.
  integer, save :: mszunit = 0                    ! newunit's numbers are negative
  logical, save :: mszopen = .false.
  integer, save :: mszcount = 0                   ! flashes lightmsz has reported in this run
  integer, parameter :: mszmax = 2000             ! lightmsz's report (about 2.5 kB a flash) stops after these
  logical, save :: started = .false.
  integer, parameter :: maxdisch = 20             ! cylinder discharges per sub-step at most
  integer, parameter :: maxflash = 120            ! branched flashes per sub-step at most (WRF-ELEC's maxiter)
  real, parameter :: dts_wrf = 0.75               ! WRF-ELEC's sub-step (s)
  real, parameter :: tgrnd = 266.16               ! WRF-ELEC's nssl_tgrnd (K)
  real, parameter :: eint = 200.0                 ! WRF-ELEC's nssl_einternal (V/m)
  real, allocatable, save :: landx(:)              ! 1 over land, 0 over water, by column along x
  double precision, save :: corpos = 0.0d0, corneg = 0.0d0   ! charge the ground's points gave off since the last log
  real, save :: coremax = 0.0, coredry = 0.0      ! the largest field at the ground (V/m) since then, before discharge,
                                                  ! over all land and where no particles reach the lowest level

contains

!-----------------------------------------------------------------------------------------------------------------------
! WRF-ELEC's NSSL driver with CM1's arguments; arrays not passed by CM1 (or not in use) stay absent
!-----------------------------------------------------------------------------------------------------------------------
  subroutine terluna_nssl_elec(th, qv, qc, qr, qi, qs, qh, qhl, cn, cna, cni, ccw, crw, cci, csw, chw, chl,            &
                               zrw, zhw, zhl, vhw, vhl, pii, p, w, dn, dz, dtp, itimestep, RAINNC, RAINNCV, dbz,       &
                               dx, dy, re_cloud, re_ice, re_snow, has_reqc, has_reqi, has_reqs, diagflag,               &
                               ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte,  &
                               pt3d, mtime)
    use input, only : ibp, iep, jbp, jep, kbp, kep, npt, var2, var3, var4, var5, var6, var8, var9, var10, nssl_3moment
    use index_module, only : lhl, lvh, lvhl, lccna, lcina
    use module_mp_nssl_2mom, only : nssl_2mom_driver, nssl_2mom_sed_driver
    use module_boxmgsetup, only : terluna_ebrk_lo, terluna_ebrk_hi
    use terluna_lightning, only : potential, field, energy, breakdown_field, discharge, eps0
    use terluna_branched, only : branched_flash, flash_record
    use terluna_screen, only : screen, earth_air_conductivity, ground_discharge
    integer, intent(in) :: ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte
    integer, intent(in) :: itimestep
    real, intent(in) :: dtp
    real, intent(in), optional :: dx, dy
    integer, intent(in), optional :: has_reqc, has_reqi, has_reqs
    logical, intent(in), optional :: diagflag
    double precision, intent(in), optional :: mtime
    real, dimension(ims:ime,jms:jme,kms:kme), intent(inout) :: qv, qc, qr, qs, qh
    real, dimension(ims:ime,jms:jme,kms:kme), intent(inout), optional :: th, qi, qhl, cn, cna, cni, ccw, crw, cci,    &
                                                   csw, chw, chl, zrw, zhw, zhl, vhw, vhl, dbz, re_cloud, re_ice, re_snow
    real, dimension(ims:ime,jms:jme,kms:kme), intent(in) :: pii, p, w, dn, dz
    real, dimension(ims:ime,jms:jme) :: RAINNC, RAINNCV
    real, dimension(ibp:iep,jbp:jep,kbp:kep,npt), intent(inout) :: pt3d

    real, allocatable, dimension(:,:,:), save :: wth, wqv, wqc, wqr, wqi, wqs, wqh, wqhl, wcn, wccw, wcrw, wcci,     &
                                                 wcsw, wchw, wchl, wvhw, wvhl, wpii, wp, ww, wdn, wdz, wdbz, wrec,     &
                                                 wrei, wres, wscw, wscr, wsci, wscs, wsch, wsciona, wschl, wsctot
    ! the field on CM1's (i,j,k) order, for the solver and the lightning schemes
    real, allocatable, dimension(:,:,:), save :: q, phi, ex, ey, ez, emag, ebrk, dep, rho3, prs3, tk3, cld3
    integer, allocatable, dimension(:,:,:), save :: used
    real, allocatable, dimension(:), save :: zh1, zf1, dz1
    integer :: ipelec, nchg, ilight, ileak, iscreen, nsub, isub, ni, nj, nk, seed, k
    real :: dts, radius, time0
    logical :: elec, hail, diag
    logical, save :: checked = .false.

    ipelec = nint(var6)
    elec = ipelec .gt. 0
    hail = present(qhl) .and. lhl .gt. 1
    nchg = 6
    if( hail ) nchg = 7
    diag = .false.
    if( present(diagflag) ) diag = diagflag
    ni = ite - its + 1
    nj = jte - jts + 1
    nk = kte - kts + 1
    ilight = nint(var8)
    ileak = nint(var9)
    iscreen = nint(var3)
    radius = var10
    if( radius .le. 0.0 ) radius = 12000.0
    if( .not. checked )then
      checked = .true.
      if( nssl_3moment ) stop 'Terluna: the electrified NSSL build carries no three-moment arrays'
      if( lccna .gt. 1 .or. lcina .gt. 1 ) stop 'Terluna: the electrified NSSL build carries no activated CCN or IN'
      if( elec .and. .not. present(ccw) ) stop 'Terluna: WRF-ELEC charges only the two-moment NSSL scheme'
      if( elec .and. npt .lt. nchg ) stop 'Terluna: electrified NSSL needs npt = 6 passive tracers (7 with hail)'
      if( elec .and. ( ibp .ne. ims .or. iep .ne. ime .or. jbp .ne. jms .or. jep .ne. jme .or. kbp .ne. kms .or.    &
                       kep .ne. kme ) ) stop 'Terluna: the passive tracers must span the model arrays'
      if( elec .and. ( ilight .lt. 0 .or. ilight .gt. 4 ) ) stop 'Terluna: var8 (lightning) takes 0 to 4'
      if( elec .and. ( iscreen .lt. 0 .or. iscreen .gt. 2 ) ) stop 'Terluna: var3 (screening layers) takes 0 to 2'
      if( ilight .eq. 2 .or. ilight .eq. 4 ) terluna_ebrk_hi = 1.0e30
      print *,'Terluna: WRF-ELEC NSSL driver, ipelec = ',ipelec,', charge tracers = ',nchg,', lightning = ',ilight,  &
              ', screening = ',iscreen,', point discharge onset (V/m at 1.225 kg/m3) = ',var2
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
        wsctot = 0.0
      endif
      if( .not. allocated(ezsave) )then
        allocate( ezsave(ims:ime,kms:kme,jms:jme) )
        ezsave = 0.0
        if( present(mtime) )then                          ! a restarted run takes up the field its restart kept
          if( mtime .gt. 0.0d0 ) call read_ez(mtime)
        endif
      endif
    endif

    ! the microphysics; with charging on, without sedimentation, which follows in sub-steps as in WRF-ELEC
    call nssl_2mom_driver(th=wth, qv=wqv, qc=wqc, qr=wqr, qi=wqi, qs=wqs, qh=wqh, qhl=wqhl, cn=wcn,                    &
                          ccw=wccw, crw=wcrw, cci=wcci, csw=wcsw, chw=wchw, chl=wchl, vhw=wvhw, vhl=wvhl,                &
                          pii=wpii, p=wp, w=ww, dn=wdn, dz=wdz, dtp=dtp, itimestep=itimestep,                           &
                          RAINNC=RAINNC, RAINNCV=RAINNCV, dbz=wdbz, dx=dx, dy=dy,                                       &
                          re_cloud=wrec, re_ice=wrei, re_snow=wres, has_reqc=has_reqc, has_reqi=has_reqi,               &
                          has_reqs=has_reqs, diagflag=diagflag,                                                         &
                          scw=wscw, scr=wscr, sci=wsci, scs=wscs, sch=wsch, schl=wschl, sciona=wsciona, sctot=wsctot,   &
                          f_sciona=elec, f_scion=.false., elecz=ezsave, ipelectmp=ipelec, induc=winduc,                &
                          noninduc=wnoninduc, isedonly_in=merge(1, 0, elec),                                            &
                          ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme,     &
                          kms=kms, kme=kme, its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)

    if( elec )then
      if( .not. started )then                             ! a run's first step: mark where it takes the logs up
        started = .true.
        if( present(mtime) ) call mark_start(real(mtime))
      endif
      if( .not. allocated(q) )then
        allocate( q(ni,nj,nk), phi(ni,nj,nk), ex(ni,nj,nk), ey(ni,nj,nk), ez(ni,nj,nk), emag(ni,nj,nk),           &
                  ebrk(ni,nj,nk), dep(ni,nj,nk), rho3(ni,nj,nk), prs3(ni,nj,nk), tk3(ni,nj,nk), cld3(ni,nj,nk),      &
                  used(ni,nj,nk), zh1(nk), zf1(nk+1), dz1(nk) )
      endif
      ! the levels: CM1's scalar levels lie midway between its w levels, the ground flat at zero
      zf1(1) = 0.0
      do k = 1,nk
        dz1(k) = wdz(its,k,jts)
        zf1(k+1) = zf1(k) + dz1(k)
        zh1(k) = zf1(k) + 0.5*dz1(k)
      enddo
      time0 = 0.0
      if( present(mtime) ) time0 = real(mtime)
      used = 0
      seed = 1 + modulo(itimestep, 2147483000)
      dts = dts_wrf
      if( var5 .gt. 0.0 ) dts = var5
      nsub = max(1, nint(dtp/dts))
      dts = dtp/nsub
      do isub = 1,nsub
        call nssl_2mom_sed_driver(th=wth, qv=wqv, qc=wqc, qr=wqr, qi=wqi, qs=wqs, qh=wqh, qhl=wqhl, cn=wcn,          &
                          ccw=wccw, crw=wcrw, cci=wcci, csw=wcsw, chw=wchw, chl=wchl, vhw=wvhw, vhl=wvhl,                &
                          pii=wpii, p=wp, w=ww, dn=wdn, dz=wdz, dtp=dts, itimestep=itimestep,                           &
                          lastloop=( isub .eq. nsub ), RAINNC=RAINNC, RAINNCV=RAINNCV, dbz=wdbz, dx=dx, dy=dy,          &
                          re_cloud=wrec, re_ice=wrei, re_snow=wres, has_reqc=has_reqc, has_reqi=has_reqi,               &
                          has_reqs=has_reqs, diagflag=diagflag,                                                         &
                          scw=wscw, scr=wscr, sci=wsci, scs=wscs, sch=wsch, schl=wschl, sciona=wsciona, sctot=wsctot,   &
                          f_sciona=.true., f_scion=.false., elecz=ezsave, ipelectmp=ipelec, isedonly_in=0,              &
                          ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme,     &
                          kms=kms, kme=kme, its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)
        call fields()
        if( ilight .ge. 1 ) call lightning(time0 + isub*dts, isub)
      enddo
      RAINNCV(its:ite,jts:jte) = RAINNC(its:ite,jts:jte)   ! the step's precipitation (CM1 zeroes RAINNC before the call)
      if( iscreen .ge. 1 ) call screening(dtp)              ! once a step, after the sub-steps, as in WRF-ELEC
      if( ileak .eq. 1 )then
        call leakage(dtp)
        call fields()
      endif
      call point_discharge(dtp)                         ! and, without it (var2 0), the field at the ground's log
      if( var2 .gt. 0.0 ) call fields()
      call save_ez()
      if( mod(itimestep, 10) .eq. 0 ) call log_field(time0 + dtp, itimestep)
      if( allocated(wrei) ) call cap_ice_size(time0 + dtp, itimestep)
    endif

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

    subroutine fields()
      ! the net charge (C/m3), the air, and the potential and field of the charge, on CM1's order
      integer :: i, j, k, ii, jj
      !$omp parallel do default(shared) private(i,j,k,ii,jj)
      do k = 1,nk
        do j = 1,nj
          jj = jts - 1 + j
          do i = 1,ni
            ii = its - 1 + i
            q(i,j,k) = wscw(ii,k,jj) + wscr(ii,k,jj) + wsci(ii,k,jj) + wscs(ii,k,jj) + wsch(ii,k,jj) + wsciona(ii,k,jj)
            if( hail ) q(i,j,k) = q(i,j,k) + wschl(ii,k,jj)
            q(i,j,k) = q(i,j,k)*wdn(ii,k,jj)
            rho3(i,j,k) = wdn(ii,k,jj)
            prs3(i,j,k) = wp(ii,k,jj)
            tk3(i,j,k) = wth(ii,k,jj)*wpii(ii,k,jj)
            cld3(i,j,k) = wqc(ii,k,jj) + wqi(ii,k,jj)
          enddo
        enddo
      enddo
      call potential(ni, nj, nk, dx, dy, zh1, zf1, q, phi)
      call field(ni, nj, nk, dx, dy, zh1, zf1, phi, ex, ey, ez, emag)
    end subroutine fields

    subroutine cap_ice_size(time, nstep)
      ! RRTMG's ice optics end at 140 microns (CM1's Morrison caps its cloud ice there too), so the cloud ice size the
      ! radiation reads is capped there. Lunar ice grows past it, and every ten steps terluna_ice_optics.txt logs the
      ! ice the cap cuts, which the radiation treats as optically thicker than its size gives
      real, intent(in) :: time
      integer, intent(in) :: nstep
      real, parameter :: rmax = 140.0e-6
      double precision :: dm, mass, cut, cutr
      real :: rbig
      integer :: i, j, k, n, unit
      mass = 0.0d0
      cut = 0.0d0
      cutr = 0.0d0
      rbig = 0.0
      n = 0
      do j = jts,jte
        do k = kts,kte
          do i = its,ite
            if( wqi(i,k,j) .gt. 0.0 )then
              dm = dble(wqi(i,k,j))*wdn(i,k,j)*dx*dy*wdz(i,k,j)
              mass = mass + dm
              if( wrei(i,k,j) .gt. rmax )then
                cut = cut + dm
                cutr = cutr + dm*wrei(i,k,j)
                n = n + 1
                if( wqi(i,k,j) .ge. 1.0e-6 ) rbig = max(rbig, wrei(i,k,j))
              endif
            endif
            wrei(i,k,j) = min(wrei(i,k,j), rmax)
          enddo
        enddo
      enddo
      if( mod(nstep, 10) .ne. 0 ) return
      call open_log('terluna_ice_optics.txt', unit)
      write(unit,'(f13.1,i9,2es12.4,2f9.1,i8)') time, nstep, mass, cut, 1.0e6*cutr/max(cut, 1.0d-30), 1.0e6*rbig, n
      close(unit)
    end subroutine cap_ice_size

    subroutine screening(dt)
      ! WRF-ELEC's screening layers from the field of the step's last solve: Earth's clear-air conductivity and a tenth
      ! of it in cloud (var3 = 1), or both by height from terluna_conductivity.txt (var3 = 2)
      real, intent(in) :: dt
      real, dimension(nk) :: sa, sc
      double precision :: qp, qn
      integer :: k
      if( iscreen .eq. 1 )then
        do k = 1,nk
          sa(k) = earth_air_conductivity(zh1(k))
          sc(k) = 0.1*sa(k)
        enddo
      else
        call read_conductivity()
        do k = 1,nk
          sa(k) = interp(zh1(k), sigclear)
          sc(k) = interp(zh1(k), sigcloud)
        enddo
      endif
      call screen(ni, nj, nk, dx, dy, zh1, zf1, dt, phi, cld3, q, sa, sc, dep, qp, qn)
      call to_ions(dep)
      if( mod(itimestep, 10) .eq. 0 ) print '(a,i9,2es12.4)', ' Terluna: screening charge this step (C):',       &
                                             itimestep, qp, qn
    end subroutine screening

    subroutine to_ions(add)
      ! charge handed to the small ions (C/m3 on CM1's order) into their tracer (C/kg) and into the net charge
      real, dimension(ni,nj,nk), intent(in) :: add
      integer :: i, j, k, ii, jj
      !$omp parallel do default(shared) private(i,j,k,ii,jj)
      do k = 1,nk
        do j = 1,nj
          jj = jts - 1 + j
          do i = 1,ni
            ii = its - 1 + i
            wsciona(ii,k,jj) = wsciona(ii,k,jj) + add(i,j,k)/wdn(ii,k,jj)
            q(i,j,k) = q(i,j,k) + add(i,j,k)
          enddo
        enddo
      enddo
    end subroutine to_ions

    subroutine lightning(time, isub)
      ! WRF-ELEC's discharges for one sub-step: cylinders until no point is above the breakdown field (maxdisch at most),
      ! or branched flashes one at a time, the field solved again after each, until no starting point is left
      real, intent(in) :: time
      integer, intent(in) :: isub
      type(flash_record) :: rec
      integer :: iter, ninit, ncol, nregion, loc(3)
      double precision :: w0, w1, qpos, qneg
      real :: e0, over
      call breakdown_field(ni, nj, nk, rho3, ebrk)
      if( ilight .le. 2 )then
        do iter = 1,maxdisch
          if( .not. any( emag .gt. ebrk ) ) exit
          loc = maxloc( emag/ebrk )
          e0 = emag(loc(1),loc(2),loc(3))
          w0 = energy(ni, nj, nk, dx, dy, zf1, q, phi)
          call discharge(ni, nj, nk, dx, dy, zf1, rho3, emag, ebrk, radius, q, dep, ninit, ncol, nregion, qpos, qneg)
          if( ncol .eq. 0 ) exit                          ! above breakdown where no column holds charge to remove
          q = q - dep                                     ! discharge has already taken it off q
          call to_ions(dep)
          call potential(ni, nj, nk, dx, dy, zh1, zf1, q, phi)
          call field(ni, nj, nk, dx, dy, zh1, zf1, phi, ex, ey, ez, emag)
          w1 = energy(ni, nj, nk, dx, dy, zf1, q, phi)
          call log_flash(time, itimestep, isub, 9, loc(1), loc(2), loc(3), e0, ebrk(loc(1),loc(2),loc(3)),          &
                         e0/ebrk(loc(1),loc(2),loc(3)), ncol*dx*dy, nregion, qpos, qneg, w0, w1, depth_lo(dep),       &
                         depth_hi(dep), ninit, 0.0d0)
          if( iter .eq. maxdisch ) print *,'Terluna: ',maxdisch,' discharges in sub-step ',isub,' of step ',itimestep
        enddo
      else
        if( .not. mszopen )then
          open( newunit=mszunit, file='terluna_msz.log', status='unknown', position='append', action='write' )
          mszopen = .true.
        endif
        do iter = 1,maxflash
          if( mszcount .eq. mszmax )then
            write(mszunit,'(a,i0,a)') '# lightmsz''s report stops after ',mszmax,' flashes; terluna_flashes.txt '//   &
                                      'has every flash'
            close(mszunit)
            open( newunit=mszunit, file='/dev/null', status='old', action='write' )
            mszcount = mszmax + 1
          endif
          w0 = energy(ni, nj, nk, dx, dy, zf1, q, phi)
          over = maxval( emag/ebrk )
          call branched_flash(ni, nj, nk, dx, dy, dz1, rho3, prs3, tk3, cld3, phi, ex, ey, ez, emag, q, ebrk, used,     &
                              seed, itimestep, dts, eint, tgrnd, var4, radius, mszunit, dep, rec)
          if( rec%kind .ge. 0 .and. mszcount .lt. mszmax ) mszcount = mszcount + 1
          if( rec%kind .le. 0 ) exit                      ! no starting point left, or no flash from the one chosen
          call to_ions(dep)
          call potential(ni, nj, nk, dx, dy, zh1, zf1, q, phi)
          call field(ni, nj, nk, dx, dy, zh1, zf1, phi, ex, ey, ez, emag)
          w1 = energy(ni, nj, nk, dx, dy, zf1, q, phi)
          call log_flash(time, itimestep, isub, rec%kind, rec%i0, rec%j0, rec%k0, rec%e0, rec%ebrk0, over, rec%area,  &
                         1, rec%qpos, rec%qneg, w0, w1, rec%zlo, rec%zhi, rec%npos + rec%nneg, rec%nox)
        enddo
        if( iter .gt. maxflash )then
          ! WRF-ELEC's driver then takes 5 % of the net charge everywhere (to the small ions) and moves on
          dep = -0.05*q
          call to_ions(dep)
          call potential(ni, nj, nk, dx, dy, zh1, zf1, q, phi)
          call field(ni, nj, nk, dx, dy, zh1, zf1, phi, ex, ey, ez, emag)
          print *,'Terluna: ',maxflash,' flashes in sub-step ',isub,' of step ',itimestep,'; 5 % of all charge removed'
        endif
      endif
    end subroutine lightning

    real function depth_lo(d)
      real, dimension(ni,nj,nk), intent(in) :: d
      integer :: k
      depth_lo = -1.0
      do k = 1,nk
        if( any( d(:,:,k) .ne. 0.0 ) )then
          depth_lo = zh1(k)
          return
        endif
      enddo
    end function depth_lo

    real function depth_hi(d)
      real, dimension(ni,nj,nk), intent(in) :: d
      integer :: k
      depth_hi = -1.0
      do k = nk,1,-1
        if( any( d(:,:,k) .ne. 0.0 ) )then
          depth_hi = zh1(k)
          return
        endif
      enddo
    end function depth_hi

    subroutine leakage(dt)
      ! net charge relaxes through the air's conductivity, its opposite going to the small ions
      real, intent(in) :: dt
      integer :: i, j, k
      real :: sig, frac
      call read_conductivity()
      !$omp parallel do default(shared) private(i,j,k,sig,frac)
      do k = 1,nk
        do j = 1,nj
          do i = 1,ni
            if( cld3(i,j,k) .ge. 1.0e-5 )then
              sig = interp(zh1(k), sigcloud)
            else
              sig = interp(zh1(k), sigclear)
            endif
            frac = 1.0 - exp( -dt*sig/eps0 )
            dep(i,j,k) = -frac*q(i,j,k)
          enddo
        enddo
      enddo
      call to_ions(dep)
    end subroutine leakage

    subroutine point_discharge(dt)
      ! point discharge from the ground (terluna_screen's ground_discharge) from the field at the lowest level, its ions
      ! joining the small ions' charge there, which CM1 then carries with the air (the small ions do not drift in the
      ! field here). WRF-ELEC's scheme hands the small ions' charge at once to any particles in a cell, so where cloud or
      ! rain reaches the lowest level the ions go to the particles. Without point discharge (var2 0) it only logs the
      ! field at the ground
      real, intent(in) :: dt
      double precision :: qp, qn
      real :: emax, onset
      integer :: i, j, ii, jj, n, nwet, unit
      logical :: wet
      call land_columns()
      onset = var2
      if( var2 .le. 0.0 ) onset = 1.0e30
      dep = 0.0
      call ground_discharge(ni, nj, dx, dy, dz1(1), dt, onset, landx(its:ite), ez(:,:,1), rho3(:,:,1), dep(:,:,1),   &
                            qp, qn, emax, n)
      nwet = 0
      do j = 1,nj
        jj = jts - 1 + j
        do i = 1,ni
          ii = its - 1 + i
          if( landx(ii) .lt. 0.5 ) cycle
          wet = wqc(ii,1,jj) + wqr(ii,1,jj) + wqi(ii,1,jj) + wqs(ii,1,jj) + wqh(ii,1,jj) .gt. 1.0e-9
          if( hail ) wet = wet .or. wqhl(ii,1,jj) .gt. 1.0e-9
          if( wet .and. dep(i,j,1) .ne. 0.0 ) nwet = nwet + 1
          if( .not. wet ) coredry = max(coredry, abs(ez(i,j,1)))
        enddo
      enddo
      if( var2 .gt. 0.0 ) call to_ions(dep)
      corpos = corpos + qp
      corneg = corneg + qn
      coremax = max(coremax, emax)
      if( mod(itimestep, 10) .eq. 0 )then
        call open_log('terluna_ground.txt', unit)
        write(unit,'(f13.1,i9,2es12.4,2f11.1,2i8)') time0 + dt, itimestep, corpos, corneg, coremax, coredry, n, nwet
        close(unit)
        corpos = 0.0d0
        corneg = 0.0d0
        coremax = 0.0
        coredry = 0.0
      endif
    end subroutine point_discharge

    subroutine land_columns()
      ! land and open water along x from terluna_surface.txt, as CM1 reads it (initsfc = 9); land throughout without it
      real :: x0, x1, xl, lu, t1, t2
      integer :: unit, nseg, m, i
      logical :: exists
      if( allocated(landx) ) return
      allocate( landx(ims:ime) )
      landx = 1.0
      inquire( file='terluna_surface.txt', exist=exists )
      if( .not. exists ) return
      open( newunit=unit, file='terluna_surface.txt', status='old', action='read' )
      read(unit,*) nseg
      do m = 1,nseg
        read(unit,*) x0, x1, xl, lu, t1, t2
        do i = ims,ime
          if( (i - 0.5)*dx .ge. x0 .and. (i - 0.5)*dx .lt. x1 ) landx(i) = merge(1.0, 0.0, nint(xl) .eq. 1)
        enddo
      enddo
      close(unit)
      print *,'Terluna: land under ',count(landx(its:ite) .gt. 0.5),' of ',ite-its+1,' columns along x'
    end subroutine land_columns

    subroutine save_ez()
      integer :: i, j, k
      !$omp parallel do default(shared) private(i,j,k)
      do j = 1,nj
        do k = 1,nk
          do i = 1,ni
            ezsave(its-1+i,k,jts-1+j) = ez(i,j,k)
          enddo
        enddo
      enddo
    end subroutine save_ez

    subroutine log_field(time, nstep)
      real, intent(in) :: time
      integer, intent(in) :: nstep
      double precision :: totpos, totneg
      real :: rates(2)
      integer :: i, j, k, unit, loc(3)
      totpos = 0.0d0
      totneg = 0.0d0
      do k = 1,nk
        do j = 1,nj
          do i = 1,ni
            if( q(i,j,k) .gt. 0.0 )then
              totpos = totpos + dble(q(i,j,k))*dx*dy*dz1(k)
            else
              totneg = totneg + dble(q(i,j,k))*dx*dy*dz1(k)
            endif
          enddo
        enddo
      enddo
      rates = (/ maxval(abs(wnoninduc(its:ite,kts:kte,jts:jte))), maxval(abs(winduc(its:ite,kts:kte,jts:jte))) /)
      loc = maxloc( emag )
      call open_log('terluna_field.txt', unit)
      write(unit,'(f13.1,i9,es11.3,3f11.1,3es12.4,2es11.3)') time, nstep, emag(loc(1),loc(2),loc(3)),                 &
            (loc(1)-0.5)*dx, (loc(2)-0.5)*dy, zh1(loc(3)), totpos, totneg, energy(ni, nj, nk, dx, dy, zf1, q, phi), rates
      close(unit)
    end subroutine log_field

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

    subroutine log_flash(time, nstep, isub, kind, i0, j0, k0, e0, ebrk0, over, extent, regions, qpos, qneg, w0, w1,  &
                         zlo, zhi, points, nox)
      real, intent(in) :: time, e0, ebrk0, over, extent, zlo, zhi
      integer, intent(in) :: nstep, isub, kind, i0, j0, k0, regions, points
      double precision, intent(in) :: qpos, qneg, w0, w1, nox
      integer :: unit
      call open_log('terluna_flashes.txt', unit)
      write(unit,'(f13.1,i9,i5,i3,3f11.1,2es11.3,f9.3,es11.3,i5,4es12.4,2f11.1,i8,es11.3)') time, nstep, isub, kind, &
            (i0-0.5)*dx, (j0-0.5)*dy, zh1(k0), e0, ebrk0, over, extent, regions, qpos, qneg, w0, w1, zlo, zhi,       &
            points, nox
      close(unit)
    end subroutine log_flash

  end subroutine terluna_nssl_elec

!-----------------------------------------------------------------------------------------------------------------------
! The vertical field the inductive charging takes from the step before, kept with each restart as WRF-ELEC keeps elecz:
! cm1.F calls this after writing a restart, and the file's name carries the time the run takes up from (whole seconds)
!-----------------------------------------------------------------------------------------------------------------------
  subroutine terluna_elec_restart_write(mtime)
    double precision, intent(in) :: mtime
    integer :: unit
    if( .not. allocated(ezsave) ) return
    open( newunit=unit, file=ez_name(mtime), access='stream', form='unformatted', status='replace', action='write' )
    write(unit) lbound(ezsave), ubound(ezsave)
    write(unit) ezsave
    close(unit)
  end subroutine terluna_elec_restart_write

  subroutine read_ez(mtime)
    double precision, intent(in) :: mtime
    integer :: unit, b(6)
    logical :: exists
    inquire( file=ez_name(mtime), exist=exists )
    if( .not. exists )then
      print *,'Terluna: no ',trim(ez_name(mtime)),'; the first step''s inductive charging takes no field'
      return
    endif
    open( newunit=unit, file=ez_name(mtime), access='stream', form='unformatted', status='old', action='read' )
    read(unit) b
    if( any( b .ne. (/ lbound(ezsave), ubound(ezsave) /) ) ) stop 'Terluna: the kept vertical field does not fit the grid'
    read(unit) ezsave
    close(unit)
    print *,'Terluna: vertical field taken up from ',trim(ez_name(mtime))
  end subroutine read_ez

  function ez_name(mtime) result(name)
    double precision, intent(in) :: mtime
    character(len=40) :: name
    write(name,'(a,i0,a)') 'terluna_ez_', nint(mtime), '.bin'
  end function ez_name

  subroutine open_log(name, unit)
    ! a log opened for appending, its column names written when it is new
    character(len=*), intent(in) :: name
    integer, intent(out) :: unit
    logical :: exists
    inquire( file=name, exist=exists )
    open( newunit=unit, file=name, status='unknown', position='append', action='write' )
    if( exists ) return
    if( name .eq. 'terluna_flashes.txt' )then
      write(unit,'(a)') '# time_s step substep kind x_m y_m z_m e_V_m e_break_V_m e_over_break_max extent_m2 '//    &
                        'regions positive_C negative_C energy_before_J energy_after_J z_low_m z_high_m points nox_mol'
      write(unit,'(a)') '# kind: 1 in cloud, 2 negative and 3 positive charge to ground (branched); 9 cylinders'
    elseif( name .eq. 'terluna_ground.txt' )then
      write(unit,'(a)') '# The field at the ground (the lowest level) over land and, where var2 sets its onset (V/m at '//  &
                        '1.225 kg/m3), point discharge (Standler and Winn 1979; J = a (E**2 - E0**2), 1 nA/m2 at 8 kV/m)'
      write(unit,'(a)') '# Every ten steps: time_s step positive_C negative_C (given off since the line before) '//      &
                        'e_ground_max_V_m e_ground_dry_max_V_m (largest field before discharge since then, over all '//  &
                        'land and where no particles reach the lowest level) columns columns_wet (discharging this step)'
    elseif( name .eq. 'terluna_ice_optics.txt' )then
      write(unit,'(a)') '# Cloud ice reaches the radiation (RRTMG) with its effective radius capped at 140 microns, '//  &
                        'where RRTMG''s ice optics end, as CM1''s Morrison caps it.'
      write(unit,'(a)') '# For the radiation and cloud-optics work: the ice above the cap is treated as optically '//    &
                        'thicker than its size gives, by its radius over 140 microns.'
      write(unit,'(a)') '# Every ten steps: time_s step ice_kg ice_above_140um_kg mean_radius_above_um '//                &
                        'max_radius_um (where ice is 0.001 g/kg or more) cells_above'
    else
      write(unit,'(a)') '# time_s step e_max_V_m x_m y_m z_m positive_C negative_C energy_J '//                     &
                        'noninductive_max_C_m3_s inductive_max_C_m3_s'
    endif
  end subroutine open_log

  subroutine mark_start(time)
    ! rows a restart replaces (those after the time the run takes up from) are dropped when the logs are read
    real, intent(in) :: time
    integer :: unit
    call open_log('terluna_ground.txt', unit)
    write(unit,'(a,f13.1)') '# run from ', time
    close(unit)
    call open_log('terluna_flashes.txt', unit)
    write(unit,'(a,f13.1)') '# run from ', time
    close(unit)
    call open_log('terluna_field.txt', unit)
    write(unit,'(a,f13.1)') '# run from ', time
    close(unit)
    call open_log('terluna_ice_optics.txt', unit)
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
