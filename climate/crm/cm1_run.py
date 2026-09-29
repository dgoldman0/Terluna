"""CM1, a cloud-resolving model, run at lunar gravity: fetch, patch, build, set up and run cases.

    climate/gcm/.venv/bin/python -m climate.crm.cm1_run build                   # Moon and Earth-gravity executables
    climate/gcm/.venv/bin/python -m climate.crm.cm1_run setup ring              # write a case's inputs
    climate/gcm/.venv/bin/python -m climate.crm.cm1_run run ring --hours 6      # run or resume, for at most 6 wall hours
    climate/gcm/.venv/bin/python -m climate.crm.cm1_run status ring
    climate/gcm/.venv/bin/python -m climate.crm.cm1_run stop ring               # stop cleanly at the next restart file

CM1 (George Bryan, NCAR; MIT-style licence) is downloaded at a pinned hash and built outside the
repository, in TERLUNA_CM1_HOME (default /media/projectspace/terluna-research/cm1). The runner patches
a fresh copy of its source for every build, the way climate/gcm/exoplasim_run.py patches PlaSim:

- gravity is a build-time constant (TERLUNA_G) in the dynamics, the RRTMG radiation (which converts
  pressure to mass paths with it) and the CAPE diagnostic;
- the Morrison scheme's fall speeds follow gravity: a particle regime V = A D^B has a Reynolds number
  growing as the Best number X^((B+1)/3), and X is proportional to g, so V scales as g^((B+1)/3) at a
  fixed size (Stokes droplets, B = 2, fall g times slower; hail, B = 0.5, as the square root of g);
- the Sun keeps a solar day of var19 seconds at var17 W/m2 above the air, starting at hour angle
  var18 degrees, and on a domain that wraps the Moon (var16 = its length in m) the hour angle grows
  eastward along x, so the terminator crosses the domain;
- the air holds the design's carbon dioxide and oxygen (TERLUNA_CO2, TERLUNA_O2), no ozone, methane,
  nitrous oxide or halocarbons;
- the soil's five layers are 5-80 cm thick, reaching below the month-long day's thermal wave;
- land and water can be set along x from a file (initsfc = 9);
- large-scale nudging of potential temperature and vapour can be confined above var14 m, reaching full
  strength at var15 m;
- with var13 = 1 a large-scale vertical wind against height is read from terluna_wls.txt into CM1's own
  large-scale vertical advection (its dolsw option), which carries temperature, vapour, condensate and
  wind with it.

Runs live in climate/crm/runs (a link to the external drive), one folder per case, and restart from
their latest restart file.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import time
import urllib.request

HERE = Path(__file__).resolve().parent
RUNS = HERE / 'runs'
CM1_HOME = Path(os.environ.get('TERLUNA_CM1_HOME', '/media/projectspace/terluna-research/cm1'))
SOURCE = dict(version='cm1r22.0', url='https://www2.mmm.ucar.edu/people/bryan/cm1/cm1r22.0.tar.gz',
              sha256='05f990bb79055873c49499c92d0ff447524a2795016a474447a7e87dd5a86226',
              licence='MIT-style (see getcode.html on the CM1 site); not redistributed here')
EARTH_G = 9.81                    # CM1's own value, kept for the Earth-gravity control

DEFAULTS = "#ifndef TERLUNA_G\n#define TERLUNA_G 9.81\n#endif\n"
GAS_DEFAULTS = ("#ifndef TERLUNA_CO2\n#define TERLUNA_CO2 400.e-6\n#endif\n"
                "#ifndef TERLUNA_O2\n#define TERLUNA_O2 0.209488\n#endif\n")

MORRISON_BLOCK = """         G = 9.81
         ! Terluna: fall speeds at the host gravity. A regime V = A D^B has Re ~ X^((B+1)/3) with the
         ! Best number X proportional to g, so V ~ g^((B+1)/3) at a fixed size; the Stokes law for
         ! cloud droplets takes G itself. The caps on mean fall speeds scale the same way.
         G = terluna_g_host
         AI = AI*(G/9.81)**((BI+1.)/3.)
         AS = AS*(G/9.81)**((BS+1.)/3.)
         AR = AR*(G/9.81)**((BR+1.)/3.)
         AG = AG*(G/9.81)**((BG+1.)/3.)
         VCAPR = 9.1*(G/9.81)**((BR+1.)/3.)
         VCAPS = 1.2*(G/9.81)**((BS+1.)/3.)
         VCAPG = 20.*(G/9.81)**((BG+1.)/3.)
         VCAPI = 1.2*(G/9.81)**((BI+1.)/3.)
"""

SUN_BLOCK = """!!!          albd = 0.07
        ENDIF

        IF( var19.gt.1.0 )THEN
          ! Terluna: a solar day of var19 s and var17 W/m2 of sunlight above the air. The hour angle
          ! starts at var18 degrees (0 = noon, -90 = sunrise); on a domain that wraps the Moon
          ! (var16 = its length in m) it also grows eastward along x: on the circle of latitude ctrlat
          ! by 360 degrees per var16, and on a great circle tilted var12 degrees to the equator, with x
          ! from its northward equator crossing, by the longitude each column reaches.
          solcon = var17
          do j=1,nj
          do i=1,ni
            hrang(i,j) = ( var18 + 360.0*real((mtime+0.5*dtrad)/var19) )*pi/180.0
            terluna_lat = ctrlat*pi/180.0
            if( var16.gt.1.0 )then
              if( var12.gt.0.0 )then
                terluna_s = 2.0*pi*xh(i)/var16
                terluna_lat = asin( sin(var12*pi/180.0)*sin(terluna_s) )
                hrang(i,j) = hrang(i,j) + atan2( cos(var12*pi/180.0)*sin(terluna_s) , cos(terluna_s) )
              else
                hrang(i,j) = hrang(i,j) + 2.0*pi*xh(i)/var16
              endif
            endif
            coszen(i,j) = min( 1.0 , max( -1.0 , cos(terluna_lat)*cos(hrang(i,j)) ) )
          enddo
          enddo
        ENDIF
"""

SURFACE_BLOCK = """      ELSEIF( initsfc.eq.9 )THEN

        ! Terluna: land and water along x from terluna_surface.txt: a segment count, then one line
        ! per segment: x_west (m), x_east (m), xland (1 land, 2 water), land-use index, tsk (K), tmn (K)
        open(unit=97,file='terluna_surface.txt',status='old',action='read')
        read(97,*) nseg
        do l=1,nseg
          read(97,*) sx1,sx2,sxl,slu,stsk,stmn
          do j=jb,je
          do i=ib,ie
            if( xh(i).ge.sx1 .and. xh(i).lt.sx2 )then
              xland(i,j) = sxl
              tsk(i,j) = stsk
            endif
          enddo
          enddo
          do j=jbl,jel
          do i=ibl,iel
            if( xh(i).ge.sx1 .and. xh(i).lt.sx2 )then
              lu_index(i,j) = slu
              tmn(i,j) = stmn
            endif
          enddo
          enddo
        enddo
        close(unit=97)

      ELSEIF( initsfc.ne.1 .and. initsfc.ne.2 )THEN
"""

LSW_BLOCK = """    IF( var13.gt.0.5 .and. var13.lt.1.5 )THEN
      ! Terluna: the large-scale vertical wind from terluna_wls.txt: a level count, then one line per level,
      ! height (m) and vertical wind (m/s), heights increasing; interpolated to the w levels, zero outside
      open(unit=97,file='terluna_wls.txt',status='old',action='read')
      read(97,*) terluna_n
      allocate( terluna_z(terluna_n) , terluna_w(terluna_n) )
      do terluna_l=1,terluna_n
        read(97,*) terluna_z(terluna_l),terluna_w(terluna_l)
      enddo
      close(unit=97)
      do k=1,nk+1
        wprof(k) = 0.0
        do terluna_l=1,terluna_n-1
          if( zf(1,1,k).ge.terluna_z(terluna_l) .and. zf(1,1,k).le.terluna_z(terluna_l+1) )then
            wprof(k) = terluna_w(terluna_l) + (terluna_w(terluna_l+1)-terluna_w(terluna_l))     &
                     *(zf(1,1,k)-terluna_z(terluna_l))/(terluna_z(terluna_l+1)-terluna_z(terluna_l))
          endif
        enddo
        if( myid.eq.0 ) print *,'  Terluna large-scale w (k, z, w): ',k,zf(1,1,k),wprof(k)
      enddo
      deallocate( terluna_z , terluna_w )
    ENDIF

    IF( var13.gt.1.5 )THEN
      ! Terluna: the large-scale vertical wind varying along x, from terluna_wls2d.txt: the column and level
      ! counts, the level heights (m), then for each column, evenly spaced from x = 0 around the ring, its
      ! vertical wind (m/s) at those heights; interpolated to the w levels, zero outside the table
      open(unit=97,file='terluna_wls2d.txt',status='old',action='read')
      read(97,*) terluna_nc,terluna_n
      allocate( terluna_z(terluna_n) , terluna_wc(terluna_nc,terluna_n) )
      read(97,*) (terluna_z(terluna_l),terluna_l=1,terluna_n)
      do terluna_c=1,terluna_nc
        read(97,*) (terluna_wc(terluna_c,terluna_l),terluna_l=1,terluna_n)
      enddo
      close(unit=97)
      if( .not. allocated(terluna_w2) ) allocate( terluna_w2(ib:ie,kb:ke+1) )
      terluna_w2 = 0.0
      do i=ib,ie
        terluna_x = ( real(i+myi1-1) - 0.5 )*dx
        terluna_c = 1 + modulo( int(floor(terluna_x/var16*real(terluna_nc))) , terluna_nc )
        do k=2,nk
          do terluna_l=1,terluna_n-1
            if( zf(1,1,k).ge.terluna_z(terluna_l) .and. zf(1,1,k).le.terluna_z(terluna_l+1) )then
              terluna_w2(i,k) = terluna_wc(terluna_c,terluna_l)                                          &
                  + (terluna_wc(terluna_c,terluna_l+1)-terluna_wc(terluna_c,terluna_l))                  &
                  *(zf(1,1,k)-terluna_z(terluna_l))/(terluna_z(terluna_l+1)-terluna_z(terluna_l))
            endif
          enddo
        enddo
      enddo
      wprof = 0.0
      do k=2,nk
        wprof(k) = sum( terluna_w2(1:ni,k) )/real(ni)
      enddo
      if( myid.eq.0 ) print *,'  Terluna large-scale w varies along x: columns, levels ',terluna_nc,terluna_n
      deallocate( terluna_z , terluna_wc )
    ENDIF

"""

WSUB_BLOCK = """    end subroutine wsub


    subroutine terluna_wsub(ix  ,jy  ,kz  ,a  ,wprof,c1,c2,mh,rr0,rf0,weps,dumz,subs,ustag)
    ! Terluna: large-scale vertical advection as in wsub, with each column's own vertical wind (terluna_w2,
    ! on w levels, averaged to u points when ustag = 1) once one has been read; otherwise wsub itself.
    use input, only : ib,ie,jb,je,kb,ke,ngxy,ngz,nk,rdz,terluna_w2
    implicit none

    integer, intent(in) :: ix,jy,kz,ustag
    real, intent(in), dimension(1-ngxy:ix+ngxy,1-ngxy:jy+ngxy,1-ngz:kz+ngz)   :: a
    real, intent(in), dimension(kb:ke) :: wprof
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2,mh,rr0,rf0
    double precision, intent(in) :: weps
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumz,subs

    integer :: i,j,k
    real :: div,wk,wk1

    if( .not. allocated(terluna_w2) )then
      call wsub(ix,jy,kz,a,wprof,c1,c2,mh,rr0,rf0,weps,dumz,subs)
      return
    endif

    DO j=1,jy
      do k=2,(nk-1)
        do i=1,ix
          wk = terluna_w2(i,k)
          if( ustag.eq.1 ) wk = 0.5*( terluna_w2(i-1,k) + terluna_w2(i,k) )
          if( wk.le.0.0 )then
            dumz(i,j,k) = rf0(1,1,k)*wk*upstrpd(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
          elseif( k.eq.2 )then
            dumz(i,j,k) = rf0(1,1,k)*wk*(c1(i,j,k)*a(i,j,k-1)+c2(i,j,k)*a(i,j,k))
          else
            dumz(i,j,k) = rf0(1,1,k)*wk*upstrpd(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),weps)
          endif
        enddo
      enddo
      do i=1,ix
        wk = terluna_w2(i,nk)
        if( ustag.eq.1 ) wk = 0.5*( terluna_w2(i-1,nk) + terluna_w2(i,nk) )
        dumz(i,j,nk) = rf0(1,1,nk)*wk*(c1(i,j,nk)*a(i,j,nk-1)+c2(i,j,nk)*a(i,j,nk))
        dumz(i,j,1) = 0.0
        dumz(i,j,nk+1) = 0.0
      enddo
      do k=1,nk
        do i=1,ix
          wk = terluna_w2(i,k)
          wk1 = terluna_w2(i,k+1)
          if( ustag.eq.1 )then
            wk = 0.5*( terluna_w2(i-1,k) + terluna_w2(i,k) )
            wk1 = 0.5*( terluna_w2(i-1,k+1) + terluna_w2(i,k+1) )
          endif
          div = (rf0(1,1,k+1)*wk1-rf0(1,1,k)*wk)*rdz*mh(1,1,k)
          subs(i,j,k) = ( -(dumz(i,j,k+1)-dumz(i,j,k))*rdz*mh(1,1,k)   &
                          +a(i,j,k)*div )*rr0(1,1,k)
        enddo
      enddo
    ENDDO

    end subroutine terluna_wsub
"""

CORIOLIS_BLOCK = """      if( var12.gt.0.0 )then
        ! Terluna: the Coriolis parameter along a great circle tilted var12 degrees to the equator, x from its
        ! northward equator crossing and var16 m round: fcor (twice the rotation rate) times sin(latitude)
        do j=jb,je
        do i=ib,ie
          f2d(i,j) = fcor*sin(var12*pi/180.0)*sin(2.0*pi*xh(i)/var16)
        enddo
        enddo
        if( myid.eq.0 ) print *,'  Terluna tilted ring: f2d at columns 1 and ni/4: ',f2d(1,1),f2d(max(1,ni/4),1)
      endif

"""

BETA_OLD = """              ! beta plane:
              uten(i,j,k)=uten(i,j,k)+0.125*(f2d(i,j)+f2d(i-1,j))           &
                                           *( (v3d(i  ,j,k)+v3d(i  ,j+1,k)) &
                                             +(v3d(i-1,j,k)+v3d(i-1,j+1,k)) )
              vten(i,j,k)=vten(i,j,k)-0.125*(f2d(i,j)+f2d(i,j-1))           &
                                           *( (u3d(i,j  ,k)+u3d(i+1,j  ,k)) &
                                             +(u3d(i,j-1,k)+u3d(i+1,j-1,k)) )
"""
BETA_NEW = """              ! beta plane:
              ! Terluna: acting on departures from the base-state wind, whose balancing pressure gradient is implied
              uten(i,j,k)=uten(i,j,k)+0.125*(f2d(i,j)+f2d(i-1,j))           &
                         *( (v3d(i  ,j,k)-v0(i  ,j,k)+v3d(i  ,j+1,k)-v0(i  ,j+1,k)) &
                           +(v3d(i-1,j,k)-v0(i-1,j,k)+v3d(i-1,j+1,k)-v0(i-1,j+1,k)) )
              vten(i,j,k)=vten(i,j,k)-0.125*(f2d(i,j)+f2d(i,j-1))           &
                         *( (u3d(i,j  ,k)-u0(i,j  ,k)+u3d(i+1,j  ,k)-u0(i+1,j  ,k)) &
                           +(u3d(i,j-1,k)-u0(i,j-1,k)+u3d(i+1,j-1,k)-u0(i+1,j-1,k)) )
"""

NUDGE_WEIGHT = ("          terluna_w = 1.0   ! Terluna: above var14 m only, fully above var15 m\n"
                "          if( var15.gt.var14 ) terluna_w = min(1.0,max(0.0,(zh(1,1,k)-var14)/(var15-var14)))\n")

ADVECTION_DECLARATIONS = """      integer, save :: terluna_nb = 0, terluna_nl = 0   ! Terluna: large-scale horizontal advection by local time
      real, dimension(:), allocatable, save :: terluna_az
      real, dimension(:,:), allocatable, save :: terluna_ath,terluna_aqv
      real :: terluna_h,terluna_f,terluna_th,terluna_qv
      integer :: terluna_b0,terluna_b1,terluna_l,terluna_m
"""

ADVECTION_BLOCK = """!--------------------------------------------------------------------
!  Terluna: large-scale horizontal advection by local time (var11 = 1), the heating and moistening the
!  planet-wide day-night circulation brings to a site. terluna_lsadv.txt holds the local-time bin count
!  and the level count, the level heights (m), then for each bin, evenly spaced in hour angle from -180
!  degrees, the theta tendency (K/s) and the vapour tendency (kg/kg/s) at each level. They are applied
!  the same across the domain, linear in local time between bin centres (round the day) and in height,
!  and zero outside the levels. The hour angle is the radiation's, var18 + 360 mtime / var19 degrees.

    IF( var11.gt.0.5 )THEN
      if( terluna_nb.eq.0 )then
        open(unit=96,file='terluna_lsadv.txt',status='old',action='read')
        read(96,*) terluna_nb,terluna_nl
        allocate( terluna_az(terluna_nl) , terluna_ath(terluna_nl,terluna_nb) , terluna_aqv(terluna_nl,terluna_nb) )
        read(96,*) terluna_az
        do terluna_m=1,terluna_nb
          read(96,*) terluna_ath(:,terluna_m)
          read(96,*) terluna_aqv(:,terluna_m)
        enddo
        close(unit=96)
        if( myid.eq.0 ) print *,'  Terluna: large-scale horizontal advection by local time; bins, levels = ',terluna_nb,terluna_nl
      endif
      terluna_h = modulo( var18 + 360.0*real(mtime/var19) + 180.0 , 360.0 )*real(terluna_nb)/360.0 - 0.5
      terluna_m = int(floor(terluna_h))
      terluna_f = terluna_h - real(terluna_m)
      terluna_b0 = modulo(terluna_m,terluna_nb) + 1
      terluna_b1 = modulo(terluna_m+1,terluna_nb) + 1
      do k=1,nk
        terluna_th = 0.0
        terluna_qv = 0.0
        if( zh(1,1,k).ge.terluna_az(1) .and. zh(1,1,k).le.terluna_az(terluna_nl) )then
          terluna_l = 1
          do terluna_m=1,terluna_nl-1
            if( zh(1,1,k).ge.terluna_az(terluna_m) ) terluna_l = terluna_m
          enddo
          tem = (zh(1,1,k)-terluna_az(terluna_l))/(terluna_az(terluna_l+1)-terluna_az(terluna_l))
          terluna_th = (1.0-terluna_f)*((1.0-tem)*terluna_ath(terluna_l,terluna_b0)+tem*terluna_ath(terluna_l+1,terluna_b0))  &
                      +terluna_f*((1.0-tem)*terluna_ath(terluna_l,terluna_b1)+tem*terluna_ath(terluna_l+1,terluna_b1))
          terluna_qv = (1.0-terluna_f)*((1.0-tem)*terluna_aqv(terluna_l,terluna_b0)+tem*terluna_aqv(terluna_l+1,terluna_b0))  &
                      +terluna_f*((1.0-tem)*terluna_aqv(terluna_l,terluna_b1)+tem*terluna_aqv(terluna_l+1,terluna_b1))
        endif
        do j=1,nj
        do i=1,ni
          thten1(i,j,k) = thten1(i,j,k) + terluna_th
        enddo
        enddo
        if( imoist.eq.1 )then
          do j=1,nj
          do i=1,ni
            qten(i,j,k,nqv) = qten(i,j,k,nqv) + terluna_qv
          enddo
          enddo
        endif
      enddo
    ENDIF

"""

# (file, marker, edits) in the order applied. An edit is (old, new, count): old must occur exactly count
# times; an empty old prepends new. A marker appears only in text its patch adds.
CPP_RULE = ('$(CPP) $(DM) $(OMP) $(DP) $(ADV) $(OUTPUTOPT) $*.F > $*.f90',
            '$(CPP) $(DM) $(OMP) $(DP) $(ADV) $(OUTPUTOPT) $(TERLUNA) $*.F > $*.f90', 1)
# The Makefile's GNU sections: MPI for 3-D runs, OpenMP for 2-D ones (CM1 needs ny >= 3 under MPI).
MAKEFILE = {
    'mpi': ('Makefile', 'Terluna: GNU compiler with MPI', [
        ('#FC = mpif90\n#OPTS = -ffree-form -ffree-line-length-none -O2 -finline-functions --param=max-vartrack-size=0 '
         '-fallow-argument-mismatch\n#CPP  = cpp -C -P -traditional -Wno-invalid-pp-token -ffreestanding\n#DM = -DMPI\n',
         '# Terluna: GNU compiler with MPI\nFC = mpif90\nOPTS = -ffree-form -ffree-line-length-none -O2 -finline-functions '
         '--param=max-vartrack-size=0 -fallow-argument-mismatch\nCPP  = cpp -C -P -traditional -Wno-invalid-pp-token '
         '-ffreestanding\nDM = -DMPI\n', 1), CPP_RULE]),
    'omp': ('Makefile', 'Terluna: GNU compiler with OpenMP', [
        ('#FC   = gfortran\n#OPTS = -ffree-form -ffree-line-length-none -O2 -finline-functions --param=max-vartrack-size=0 '
         '-fopenmp\n#CPP  = cpp -C -P -traditional -Wno-invalid-pp-token -ffreestanding\n#OMP  = -DOPENMP\n',
         '# Terluna: GNU compiler with OpenMP\nFC   = gfortran\nOPTS = -ffree-form -ffree-line-length-none -O2 '
         '-finline-functions --param=max-vartrack-size=0 -fopenmp\nCPP  = cpp -C -P -traditional -Wno-invalid-pp-token '
         '-ffreestanding\nOMP  = -DOPENMP\n', 1), CPP_RULE]),
}
PATCHES = [
    ('constants.F', 'Terluna: gravity set at build time', [
        ('', DEFAULTS, 1),
        ('        g      = 9.81\n        rd     = 287.04\n',
         '        g      = TERLUNA_G   ! Terluna: gravity set at build time\n        rd     = 287.04\n', 1)]),
    ('module_ra_etc.F', 'Terluna: gravity for RRTMG', [
        ('', DEFAULTS, 1),
        ('   REAL    , PARAMETER :: g = 9.81  ! acceleration due to gravity (m {s}^-2)',
         '   REAL    , PARAMETER :: g = TERLUNA_G  ! Terluna: gravity for RRTMG, set at build time', 1)]),
    ('getcape.F', 'Terluna: gravity for CAPE', [
        ('', DEFAULTS, 1),
        ('    real, parameter :: g     = 9.81\n', '    real, parameter :: g     = TERLUNA_G   ! Terluna: gravity for CAPE\n', 1)]),
    ('morrison.F', 'Terluna: fall speeds at the host gravity', [
        ("     REAL, PRIVATE ::      BI,BC,BS,BR,BG ! 'B' PARAMETER IN FALLSPEED-DIAM RELATIONSHIP\n",
         "     REAL, PRIVATE ::      BI,BC,BS,BR,BG ! 'B' PARAMETER IN FALLSPEED-DIAM RELATIONSHIP\n"
         "     REAL, PRIVATE ::      VCAPR,VCAPS,VCAPG,VCAPI ! Terluna: fall-speed caps\n", 1),
        ('SUBROUTINE GRAUPEL_INIT(cm1hail,cm1inum,cm1ndcnst,cm1db)\n',
         'SUBROUTINE GRAUPEL_INIT(cm1hail,cm1inum,cm1ndcnst,cm1db)\n      use constants, only : terluna_g_host => g\n', 1),
        ('         G = 9.81\n', MORRISON_BLOCK, 1),
        ('UMS=MIN(UMS,1.2*dum)', 'UMS=MIN(UMS,VCAPS*dum)', 3),
        ('UNS=MIN(UNS,1.2*dum)', 'UNS=MIN(UNS,VCAPS*dum)', 3),
        ('UMR=MIN(UMR,9.1*dum)', 'UMR=MIN(UMR,VCAPR*dum)', 5),
        ('UNR=MIN(UNR,9.1*dum)', 'UNR=MIN(UNR,VCAPR*dum)', 5),
        ('UMG=MIN(UMG,20.*dum)', 'UMG=MIN(UMG,VCAPG*dum)', 3),
        ('UNG=MIN(UNG,20.*dum)', 'UNG=MIN(UNG,VCAPG*dum)', 3),
        ('UMI=MIN(UMI,1.2*(rhosu/rho(k))**0.35)', 'UMI=MIN(UMI,VCAPI*(rhosu/rho(k))**0.35)', 1),
        ('UNI=MIN(UNI,1.2*(rhosu/rho(k))**0.35)', 'UNI=MIN(UNI,VCAPI*(rhosu/rho(k))**0.35)', 1)]),
    ('radiation_driver.F', 'Terluna: a solar day of var19 s', [
        ('      use constants, only : pi,g,cp,cpl,cpi,cv,cvv,rd,cvdcp,degdpi\n',
         '      use constants, only : pi,g,cp,cpl,cpi,cv,cvv,rd,cvdcp,degdpi\n'
         '      use input, only : var12,var16,var17,var18,var19   ! Terluna: the lunar day\n', 1),
        ('      real :: saltitude,sazimuth,zen,rtime\n',
         '      real :: saltitude,sazimuth,zen,rtime\n      real :: terluna_lat,terluna_s   ! Terluna: tilted rings\n', 1),
        ('!!!          albd = 0.07\n        ENDIF\n', SUN_BLOCK, 1)]),
    ('module_ra_rrtmg_lw.F', 'Terluna: the design air, longwave', [
        ('', GAS_DEFAULTS + '! Terluna: the design air, longwave\n', 1),
        ('      co2 = (280. + 90.*exp(0.02*(yr-2000)))*1.e-6\n', '      co2 = TERLUNA_CO2\n', 1),
        ('    data ch4 / 1774.e-9 /', '    data ch4 / 0.0 /', 1),
        ('    data n2o / 319.e-9 /', '    data n2o / 0.0 /', 1),
        ('    data cfc11 / 0.251e-9 /', '    data cfc11 / 0.0 /', 1),
        ('    data cfc12 / 0.538e-9 /', '    data cfc12 / 0.0 /', 1),
        ('    data cfc22 / 0.169e-9 /', '    data cfc22 / 0.0 /', 1),
        ('    data ccl4 / 0.093e-9 /', '    data ccl4 / 0.0 /', 1),
        ('    data o2 / 0.209488 /', '    data o2 / TERLUNA_O2 /', 1),
        ('o3vmr(ncol,k) = o3mmr(k) * amdo', 'o3vmr(ncol,k) = 0.0*o3mmr(k) * amdo', 2)]),
    ('module_ra_rrtmg_sw.F', 'Terluna: the design air, shortwave', [
        ('', GAS_DEFAULTS + '! Terluna: the design air, shortwave\n', 1),
        ('      co2 = (280. + 90.*exp(0.02*(yr-2000)))*1.e-6\n', '      co2 = TERLUNA_CO2\n', 1),
        ('    data ch4 / 1774.e-9 /', '    data ch4 / 0.0 /', 1),
        ('    data n2o / 319.e-9 /', '    data n2o / 0.0 /', 1),
        ('    data o2 / 0.209488 /', '    data o2 / TERLUNA_O2 /', 1),
        ('o3vmr(ncol,k) = o3mmr(k) * amdo', 'o3vmr(ncol,k) = 0.0*o3mmr(k) * amdo', 2)]),
    ('init_surface.F', 'Terluna: land and water along x', [
        ('      integer :: i,j,k,l\n      real :: x1,x2,xcoast\n',
         '      integer :: i,j,k,l\n      real :: x1,x2,xcoast\n      integer :: nseg,slu\n      real :: sx1,sx2,sxl,stsk,stmn\n', 1),
        ('      ELSEIF( initsfc.ne.1 .and. initsfc.ne.2 )THEN\n', SURFACE_BLOCK, 1),
        ('      slab_dzs(1)=.01\n', '      slab_dzs(1)=.05   ! Terluna: layers of 5-80 cm\n', 1)]),
    ('solve1.F', 'Terluna: above var14 m only', [
        ('      real :: dttmp,rtime,rdt,tem,tem0,tem1,tem2,thrad,prad\n',
         '      real :: dttmp,rtime,rdt,tem,tem0,tem1,tem2,thrad,prad,terluna_w\n', 1),
        ('          tem1 = -lsnudgefac*( thavg(k)-lsnudge_th(k,lsnudge_count) )/(lsnudge_tau)\n',
         NUDGE_WEIGHT + '          tem1 = -lsnudgefac*terluna_w*( thavg(k)-lsnudge_th(k,lsnudge_count) )/(lsnudge_tau)\n', 1),
        ('          tem1 = -lsnudgefac*( qavg(k,nqv)-lsnudge_qv(k,lsnudge_count) )/(lsnudge_tau)\n',
         NUDGE_WEIGHT + '          tem1 = -lsnudgefac*terluna_w*( qavg(k,nqv)-lsnudge_qv(k,lsnudge_count) )/(lsnudge_tau)\n', 1)]),
    ('param.F', 'Terluna: large-scale vertical wind from terluna_wls.txt', [
        ('      IF( testcase.eq.7 ) dolsw = .true.\n',
         '      IF( testcase.eq.7 ) dolsw = .true.\n'
         '      IF( var13.gt.0.5 ) dolsw = .true.   ! Terluna: large-scale vertical wind from terluna_wls.txt\n', 1)]),
    ('base.F', 'Terluna: the large-scale vertical wind from terluna_wls.txt', [
        ('      real :: z1,z2,z3,z4,z5,z6,z7,z8,z9\n',
         '      real :: z1,z2,z3,z4,z5,z6,z7,z8,z9\n      integer :: terluna_n,terluna_l,terluna_nc,terluna_c\n'
         '      real :: terluna_x\n      real, dimension(:), allocatable :: terluna_z,terluna_w\n'
         '      real, dimension(:,:), allocatable :: terluna_wc\n', 1),
        ('    ! boundary conditions:\n    wprof(1) = 0.0\n    wprof(nk+1) = 0.0\n',
         LSW_BLOCK + '    ! boundary conditions:\n    wprof(1) = 0.0\n    wprof(nk+1) = 0.0\n', 1)]),
    ('input.F', 'Terluna: large-scale vertical wind by column', [
        ('  MODULE input\n\n  implicit none\n\n  public\n',
         '  MODULE input\n\n  implicit none\n\n  public\n\n      real, dimension(:,:), allocatable :: terluna_w2'
         '   ! Terluna: large-scale vertical wind by column and w level\n', 1)]),
    ('adv_routines.F', 'Terluna: large-scale vertical advection as in wsub', [
        ('  public :: movesfc,wsub,zsgrad\n', '  public :: movesfc,wsub,zsgrad,terluna_wsub\n', 1),
        ('    end subroutine wsub\n', WSUB_BLOCK, 1)]),
    ('adv.F', 'Terluna: w by column', [
        ('use adv_routines, only : advsaxi,wsub,', 'use adv_routines, only : advsaxi,wsub,terluna_wsub,', 1),
        ('use adv_routines, only : advuaxi,wsub,', 'use adv_routines, only : advuaxi,wsub,terluna_wsub,', 1),
        ('use adv_routines, only : advvaxi,wsub,', 'use adv_routines, only : advvaxi,wsub,terluna_wsub,', 1),
        ('call     wsub(ni  ,nj  ,nk  ,s  ,wprof,c1,c2,mh,rr0,rf0,weps,dumz,subs)',
         'call terluna_wsub(ni  ,nj  ,nk  ,s  ,wprof,c1,c2,mh,rr0,rf0,weps,dumz,subs,0)   ! Terluna: w by column', 1),
        ('call     wsub(ni+1,nj  ,nk  ,u3d,wprof,c1,c2,mh,rr0,rf0,weps,dumz,subs)',
         'call terluna_wsub(ni+1,nj  ,nk  ,u3d,wprof,c1,c2,mh,rr0,rf0,weps,dumz,subs,1)   ! Terluna: w by column', 1),
        ('call     wsub(ni  ,nj+1,nk  ,v3d,wprof,c1,c2,mh,rr0,rf0,weps,dumz,subs)',
         'call terluna_wsub(ni  ,nj+1,nk  ,v3d,wprof,c1,c2,mh,rr0,rf0,weps,dumz,subs,0)   ! Terluna: w by column', 1)]),
    ('param.F', 'Terluna: the Coriolis parameter along a great circle', [
        ('      if( axisymm.eq.1 .or. ny.le.3 )  betaplane = 0\n',
         '      if( axisymm.eq.1 .or. ( ny.le.3 .and. var12.le.0.0 ) )  betaplane = 0'
         '   ! Terluna: kept on a tilted ring, where f varies along x\n', 1),
        ('    if( lspgrad.eq.3 )then\n      ! for lspgrad = 3\n',
         CORIOLIS_BLOCK + '    if( lspgrad.eq.3 )then\n      ! for lspgrad = 3\n', 1)]),
    ('solve2.F', 'Terluna: acting on departures from the base-state wind', [(BETA_OLD, BETA_NEW, 1)]),
    ('solve1.F', 'Terluna: large-scale horizontal advection by local time', [
        ('      real :: dttmp,rtime,rdt,tem,tem0,tem1,tem2,thrad,prad,terluna_w\n',
         '      real :: dttmp,rtime,rdt,tem,tem0,tem1,tem2,thrad,prad,terluna_w\n' + ADVECTION_DECLARATIONS, 1),
        ('!-------------------------------------------------------------------\n!    NOTE:  now ok to change dum7,dum8\n',
         ADVECTION_BLOCK + '!-------------------------------------------------------------------\n'
         '!    NOTE:  now ok to change dum7,dum8\n', 1)]),
]


def apply_patch(text: str, name: str, marker: str, edits) -> str:
    """The text with one patch applied; unchanged if its marker is already present."""
    if marker in text:
        return text
    for old, new, count in edits:
        if not old:
            text = new + text
            continue
        found = text.count(old)
        if found != count:
            raise RuntimeError(f'CM1 {name}: expected {count} of {old.strip()[:60]!r}, found {found}')
        text = text.replace(old, new)
    if marker not in text:
        raise RuntimeError(f'CM1 {name}: the patch did not add its marker {marker!r}')
    return text


def fetch() -> Path:
    """The pristine CM1 tree, downloaded and unpacked if missing, its tarball checked against the pin."""
    CM1_HOME.mkdir(parents=True, exist_ok=True)
    tarball = CM1_HOME / f"{SOURCE['version']}.tar.gz"
    if not tarball.exists():
        part = tarball.with_name(tarball.name + '.part')
        with urllib.request.urlopen(SOURCE['url']) as response, open(part, 'wb') as out:
            shutil.copyfileobj(response, out)
        part.replace(tarball)
    digest = hashlib.sha256(tarball.read_bytes()).hexdigest()
    if digest != SOURCE['sha256']:
        raise RuntimeError(f'{tarball} has sha256 {digest}, expected {SOURCE["sha256"]}; not using it')
    tree = CM1_HOME / SOURCE['version']
    if not (tree / 'src' / 'Makefile').exists():
        with tarfile.open(tarball) as archive:
            archive.extractall(CM1_HOME, filter='data')
    return tree


def planet():
    return json.loads((HERE.parent / 'gcm' / 'products' / 'moon_gcm_configuration.json').read_text())['planet']


def design_air():
    """Carbon dioxide and oxygen mixing ratios of the design case's air, from its GCM run."""
    progress = json.loads((HERE.parent / 'gcm' / 'runs' / GCM_RUN / 'progress.json').read_text())
    cfg = progress['configuration']
    gases = cfg['gases_bar']
    total = sum(gases.values())
    return dict(co2=gases['pCO2'] / total, o2=gases['pO2'] / total, surface_pa=cfg['pressure_pa'],
                sunlight_w_m2=progress['inputs']['flux_w_m2'], run=GCM_RUN)


# label: (gravity, parallel mode)
BUILDS = {'moon': (lambda: planet()['gravity_m_s2'], 'mpi'), 'earth_g': (lambda: EARTH_G, 'mpi'),
          'moon_omp': (lambda: planet()['gravity_m_s2'], 'omp'), 'earth_g_omp': (lambda: EARTH_G, 'omp')}


def build(label: str, jobs: int = 8) -> Path:
    """Compile CM1 with the Terluna patches at the build's gravity; returns the executable."""
    tree = fetch()
    air = design_air()
    gravity, mode = BUILDS[label][0](), BUILDS[label][1]
    folder = CM1_HOME / 'build' / label
    exe = folder / 'cm1.exe'
    flags = f"-DTERLUNA_G={gravity!r} -DTERLUNA_CO2={air['co2']:.6e} -DTERLUNA_O2={air['o2']:.6f}"
    src = folder / 'src'
    if src.exists():
        shutil.rmtree(src)
    shutil.copytree(tree / 'src', src)
    (folder / 'run').mkdir(parents=True, exist_ok=True)
    hashes = {}
    for name, marker, edits in [MAKEFILE[mode], *PATCHES]:
        path = src / name
        path.write_text(apply_patch(path.read_text(), name, marker, edits))
    for name in sorted({p[0] for p in [MAKEFILE[mode], *PATCHES]}):
        hashes[name] = hashlib.sha256((src / name).read_bytes()).hexdigest()[:16]
    log = folder / 'build.log'
    with open(log, 'w') as out:
        done = subprocess.run(['make', f'-j{jobs}', f'TERLUNA={flags}'], cwd=src, stdout=out, stderr=subprocess.STDOUT)
    if done.returncode != 0 or not (folder / 'run' / 'cm1.exe').exists():
        raise RuntimeError(f'CM1 build {label} failed; see {log}')
    (folder / 'run' / 'cm1.exe').replace(exe)
    record = dict(label=label, gravity_m_s2=gravity, parallel=mode, flags=flags, source=SOURCE, patched=hashes,
                  executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest()[:16], air=air)
    (folder / 'build.json').write_text(json.dumps(record, indent=1) + '\n')
    return exe


# --------------------------------------------------------------------------------------------------
# Cases

# The chosen design (5% dimmer shield) with the GCM's Sun, tilt and sunlight split corrected, settled years. Every
# case set up before 2026-09-29 took its forcing from A28_dim5, years 15-24.
GCM_RUN, GCM_YEARS = 'A28_dim5_moon', (20, 29)
EQUATOR_BAND_DEG = 10.0                              # GCM rows averaged for the equatorial profile
ROW_BAND_DEG = 3.0                                   # half the GCM's row spacing: the one row nearest a latitude
PLACEHOLDER_LAND = (28, "28,      20.,   .50,   .95,   10.,    4.,  2.00, 25.0e5,'Terluna placeholder land'\n")
# PlaSim's bucket (landmod.f90): soil holds up to wsmax m of water and land evaporates at the rate of open
# water times min(1, water / (drhsfull * wsmax)), the same factor as CM1's moisture availability.
PLASIM_SOIL = dict(wsmax=0.5, drhsfull=0.4)

CASES = {
    'ring': dict(
        purpose='the equator as a 2-D ring 10,917 km round, the Sun crossing it once a lunar day; seas and lakes '
                'of the 28% scenario at sea level, the 5% dimmer shield, the GCM design case as reference',
        build='moon_omp', dx_target_m=6000.0, ranks=8, z_grid=dict(dz0=100.0, z_fine=1000.0, stretch=1.08, dz_max=2000.0,
                                                                ztop=150000.0),
        start_hour_angle_deg=-90.0, days=59.06, output_s=10800.0, restart_s=43200.0, segment_s=86400.0,
        nudge=dict(tau_s=259200.0, z_start_m=8000.0, z_full_m=16000.0, ramp_s=86400.0),
        land_moisture=0.50, droplets_cm3=100.0),
    # The gravity pair: the same air in pressure terms over the same sea, at lunar and at Earth gravity.
    # Lengths and times of the Earth case are the Moon's times g_moon/g_earth, as the dynamics scale, so any
    # difference comes from what does not scale: how fast drops and ice fall and how fast rain forms. No
    # radiation. The air starts as the ring's own air over its seas, and above 8-16 km (lunar) its mean is
    # held to the reference, which takes away the heat and vapour convection brings up, as in
    # weak-temperature-gradient experiments; below, the boundary layer and convection are free. (Holding the
    # air at every height kept the boundary layer at the reference's mean humidity, too dry to convect.)
    'pair_moon': dict(kind='pair', build='moon_omp', nx=256, dx_moon_m=6000.0,
                      z_grid=dict(dz0=100.0, z_fine=1000.0, stretch=1.08, dz_max=2000.0, ztop=150000.0),
                      days=20.0, output_s=10800.0, restart_s=86400.0, segment_s=432000.0, initial='ring_sea',
                      nudge=dict(tau_s=259200.0, z_start_m=8000.0, z_full_m=16000.0, ramp_s=86400.0), droplets_cm3=100.0),
    'pair_earth': dict(kind='pair', build='earth_g_omp', nx=256, dx_moon_m=6000.0,
                       z_grid=dict(dz0=100.0, z_fine=1000.0, stretch=1.08, dz_max=2000.0, ztop=150000.0),
                       days=20.0, output_s=10800.0, restart_s=86400.0, segment_s=432000.0, initial='ring_sea',
                       nudge=dict(tau_s=259200.0, z_start_m=8000.0, z_full_m=16000.0, ramp_s=86400.0), droplets_cm3=100.0),
}
# Rings along circles of latitude near the poles, set up as the equatorial ring except that the Sun never
# climbs above 90 degrees less the latitude, the Coriolis force of the latitude turns departures from the
# reference wind (whose own pressure gradient CM1 supplies, lspgrad = 1), the reference is the GCM row
# nearest the latitude, and the land evaporates as readily as the GCM's land on that row. The Sun keeps the
# equinox: the Moon's 1.5-degree tilt raises and lowers the daily sunlight at 80 degrees by about a quarter
# over the year, alternately in the two hemispheres.
for _lat, _tag in ((80.0, '80n'), (-80.0, '80s')):
    CASES[f'ring_{_tag}'] = dict(
        CASES['ring'], latitude_deg=_lat, band_deg=ROW_BAND_DEG, land_moisture='gcm',
        purpose=f'the circle of latitude {abs(_lat):.0f} {"N" if _lat > 0 else "S"} as a 2-D ring, the Sun crossing it '
                'once a lunar day, seas and lakes of the 28% scenario at sea level, land as wet as the GCM\'s there, '
                'the 5% dimmer shield, the GCM design case\'s row nearest the latitude as reference')
# The same with the GCM's mean vertical wind at the latitude imposed as large-scale vertical advection: a closed
# ring cannot rise or sink on average, and the GCM's air sinks over the polar caps at about 2-7 mm/s through
# most of its depth (the descending branch of its overturning between the equator and the poles).
for _lat, _tag in ((80.0, '80n'), (-80.0, '80s'), (70.0, '70n'), (-70.0, '70s'), (45.0, '45n'), (-45.0, '45s')):
    CASES[f'ring_{_tag}_lsw'] = dict(
        CASES['ring'], latitude_deg=_lat, band_deg=ROW_BAND_DEG, land_moisture='gcm', large_scale_w='gcm',
        purpose=f'the circle of latitude {abs(_lat):.0f} {"N" if _lat > 0 else "S"} as a 2-D ring, the Sun crossing it '
                'once a lunar day, seas and lakes of the 28% scenario at sea level, land as wet as the GCM\'s there, '
                'the 5% dimmer shield, the GCM design case\'s row nearest the latitude as reference, and its mean '
                'vertical wind there imposed as large-scale vertical advection')


def band_vertical_wind(lat, sigma_edges, flux, rows, radius, gravity):
    """Mean vertical wind (m/s, up positive) at the GCM's layer interfaces over the band covered by `rows`
    (a mask on lat, rows ordered north to south), from its mass budget. `flux` is the zonal and time mean of
    surface pressure times northward wind by layer and row (Pa m/s). Its column mean is removed first, since
    no net mass crosses a latitude circle in a steady state; then what converges on the band above an
    interface sinks through it. The band's edges lie midway to the next rows out."""
    import numpy as np
    lat, sigma_edges, flux = (np.asarray(a, dtype=float) for a in (lat, sigma_edges, flux))
    dsig = np.diff(sigma_edges)
    flux = flux - (flux * dsig[:, None]).sum(axis=0) / dsig.sum()
    above = np.vstack([np.zeros(lat.size), np.cumsum(flux * dsig[:, None], axis=0)]) / gravity   # kg/m/s
    idx = np.flatnonzero(rows)
    jn, js = idx.min(), idx.max()
    edge_n, edge_s = 0.5 * (lat[jn] + lat[jn - 1]), 0.5 * (lat[js] + lat[js + 1])
    across = lambda edge, j0, j1: 2 * np.pi * radius * np.cos(np.radians(edge)) * 0.5 * (above[:, j0] + above[:, j1])
    into = across(edge_s, js, js + 1) - across(edge_n, jn, jn - 1)                 # kg/s converging above each interface
    area = 2 * np.pi * radius ** 2 * (np.sin(np.radians(edge_n)) - np.sin(np.radians(edge_s)))
    return -into / area, (float(edge_s), float(edge_n))                             # kg/m2/s, up positive


def plasim_half_levels(sigma):
    """PlaSim's layer interfaces from its full levels: the top at 0 and each full level midway between the
    interfaces around it (plasim.f90). The output's levp holds midpoints between full levels instead."""
    import numpy as np
    edges = [0.0]
    for s in np.asarray(sigma, dtype=float):
        edges.append(2.0 * s - edges[-1])
    return np.array(edges)


def gcm_vertical_wind(latitude, band, run=GCM_RUN, years=GCM_YEARS) -> dict:
    """The GCM's mean vertical wind over the band of rows within `band` degrees of a latitude, against height
    above the ground (from the rows' mean virtual temperatures), ready for terluna_wls.txt."""
    import numpy as np
    import netCDF4
    gravity, rd = planet()['gravity_m_s2'], 287.04
    folder = HERE.parent / 'gcm' / 'runs' / run / 'model'
    flux, temp, vap, surface = [], [], [], []
    for year in range(years[0], years[1] + 1):
        with netCDF4.Dataset(folder / f'MOST.{year:05d}.nc') as d:
            ps = np.asarray(d['ps'][:], float) * 100.0                              # hPa in the output
            flux.append((ps[:, None] * np.asarray(d['va'][:], float)).mean(axis=(0, 3)))
            temp.append(np.asarray(d['ta'][:], float).mean(axis=(0, 3)))
            vap.append(np.asarray(d['hus'][:], float).mean(axis=(0, 3)))
            surface.append(ps.mean(axis=(0, 2)))
            lat = np.asarray(d['lat'][:], float)
            sigma = np.asarray(d['lev'][:], float)
    flux, temp, vap, surface = (np.mean(a, axis=0) for a in (flux, temp, vap, surface))
    edges = plasim_half_levels(sigma)
    rows = np.abs(lat - latitude) < band
    mass, span = band_vertical_wind(lat, edges, flux, rows, planet()['radius_m'], gravity)
    t, q, p = temp[:, rows].mean(axis=1), vap[:, rows].mean(axis=1), surface[rows].mean()
    z = np.zeros(edges.size)                                                        # interface heights, top first
    for k in range(sigma.size - 1, -1, -1):
        if edges[k] > 0:
            z[k] = z[k + 1] + rd * t[k] * (1 + 0.608 * q[k]) / gravity * np.log(edges[k + 1] / edges[k])
    rho = p * edges / (rd * np.interp(edges, sigma, t))
    w = np.where(edges > 0, mass / np.maximum(rho, 1e-9), 0.0)
    keep = edges > 0                                                                # the top interface has no height
    return dict(z_m=[float(v) for v in z[keep][::-1]], w_m_s=[float(v) for v in w[keep][::-1]], band_edges_deg=list(span),
                rows_deg=[float(v) for v in lat[rows]], run=run, years=list(years))


def write_vertical_wind(path: Path, profile: dict, ztop: float) -> None:
    """terluna_wls.txt: the level count, then height (m) and vertical wind (m/s) from the ground up, ending at
    zero at the model top."""
    z, w = list(profile['z_m']), list(profile['w_m_s'])
    if z[-1] < ztop:
        z, w = z + [ztop], w + [0.0]
    path.write_text(f'{len(z)}\n' + ''.join(f'{a:.3f} {b:.6e}\n' for a, b in zip(z, w)))

# Rings along great circles tilted 45 degrees to the equator, in two mirrored pairs: A heads north across the
# equator at 0 E and A-prime at 180 E, so the two cross there and at 180 E, one at 45 N where the other is at
# 45 S; B and B-prime are the pair turned 90 degrees. Every column takes its latitude, longitude, Sun, Coriolis
# parameter, surface and the GCM's mean vertical wind from its place on the path; the upper air is held to one
# profile, the mean along the path (the GCM's air is nearly uniform with latitude).
for _name, _node in (('ring_a', 0.0), ('ring_a_prime', 180.0), ('ring_b', 90.0), ('ring_b_prime', 270.0)):
    CASES[_name] = dict(
        CASES['ring'], kind='tilted', tilt_deg=45.0, node_deg=_node, land_moisture='gcm', large_scale_w='gcm',
        purpose=f'a great circle tilted 45 degrees to the equator, heading north across it at {_node:.0f} E, as a 2-D '
                'ring the Sun crosses once a lunar day, seas and lakes of the 28% scenario at sea level, land as wet '
                'as the GCM\'s beneath it, the 5% dimmer shield, and the GCM design case\'s mean vertical wind at each '
                'column\'s latitude')
# Rings forced from the corrected GCM (set up on 2026-09-29). Two great circles tilted 70 degrees, heading north across
# the equator at 45 and 135 E: together they reach 70 degrees on either side with 3,200-4,600 km of land in each band
# of latitude, and they cross each other nearly at right angles (83 degrees) at 62.8 N, 180 E and 62.8 S, 0 E, where
# comfort is in question. With them the equatorial ring again, set up as a great circle of no tilt so that, as on the
# tilted rings, each column takes its own surface, land wetness and ground temperature and the GCM's vertical wind;
# it crosses each steep ring on the equator.
for _name, _tilt, _node in (('ring_70_45e', 70.0, 45.0), ('ring_70_135e', 70.0, 135.0), ('ring_equator', 0.0, 0.0)):
    CASES[_name] = dict(
        CASES['ring'], kind='tilted', tilt_deg=_tilt, node_deg=_node, land_moisture='gcm', large_scale_w='gcm',
        purpose=(f'a great circle tilted {_tilt:.0f} degrees to the equator, heading north across it at {_node:.0f} E,'
                 if _tilt > 0.0 else 'the equator, a great circle of no tilt,')
                + ' as a 2-D ring the Sun crosses once a lunar day, seas and lakes of the 28% scenario at sea level, land as '
                'wet as the GCM\'s beneath it, the 5% dimmer shield, and the corrected GCM design case\'s mean vertical '
                'wind at each column\'s latitude')
WETNESS_CLASSES =((20, 0.05), (21, 0.10), (22, 0.15), (23, 0.20), (25, 0.30), (26, 0.40), (27, 0.50), (28, 0.60),
                   (29, 0.75), (30, 0.90))                    # land-use rows given over to the placeholder land at these wetnesses
# A three-dimensional box where rings A and A-prime cross (0 N, 0 E): ring A's patch there given a second horizontal
# dimension, to show what the rings' two-dimensionality does. It keeps ring A's column spacing, levels, upper-air
# reference and wind along the ring (its x axis runs along ring A, 45 degrees east of north there), and takes the
# land wetness, ground temperature and GCM vertical wind of ring A's columns within 100 km of the crossing, all
# land. The Sun is the site's, the same across the box; on the equator there is no Coriolis force. The box cannot
# make the planet-wide day-night circulation that the rings carry along their length.
# The box also takes the heating and moistening the rings' own day-night circulation brings there: the advection by
# the flow of each ring's longest waves (wavelengths of 2,180 km and more), along the ring and up or down, the mean of
# rings A and A-prime over their second lunar day by local time, refitted with the day's first four harmonics, up to
# 16 km.
CASES['box_0e'] = dict(
    CASES['ring'], kind='box', nx=64, ny=64, site=dict(ring='ring_a', lat_deg=0.0, lon_deg=0.0, radius_m=100.0e3),
    land_moisture='gcm', large_scale_w='gcm',
    day_night=dict(rings=('ring_a', 'ring_a_prime'), from_day=29.5, wavenumbers=5, bins=36, harmonics=4, top_m=16000.0),
    purpose='ring A\'s patch where it crosses A-prime (0 N, 0 E) as a 3-D box 385 km square, all land as wet as the '
            'GCM\'s beneath the patch, the Sun crossing it once a lunar day, the 5% dimmer shield, ring A\'s upper-air '
            'reference and the GCM design case\'s mean vertical wind there')
# Tests of the land surface every CM1 case shares, at the same site and under the same forcing in boxes a quarter the
# size (32 x 32 columns, 192 km square): the placeholder land as in box_0e, then the GCM's land one property at a time
# and all together. PlaSim's land (landmod.f90 as run) has a roughness length of 2 m everywhere, a soil of thermal
# inertia sqrt(1.8 W/m/K x 2.4e6 J/m3/K) = 2,080 SI (CM1's slab soil, whose diffusivity is fixed, takes it as THERIN
# 4.97 against the placeholder's 4) and the placeholder's albedo of 0.20; at 0 E its bucket stays full (wetness 1).
# These boxes also write the radiation at the ground.
GCM_LAND = dict(roughness_cm=200.0, therin=4.97, wetness=1.0)
for _tag, _land in (('', None), ('_rough', dict(roughness_cm=GCM_LAND['roughness_cm'])),
                    ('_soil', dict(therin=GCM_LAND['therin'])), ('_gcm_land', GCM_LAND)):
    CASES[f'box_0e_small{_tag}'] = dict(
        CASES['box_0e'], nx=32, ny=32, land=_land, surface_output=True,
        purpose='box_0e at a quarter of its size (192 km square), ' + (
            'with the placeholder land' if _land is None else 'with the placeholder land given the GCM\'s '
            + ', '.join({'roughness_cm': 'roughness (2 m)', 'therin': 'soil thermal inertia',
                         'wetness': 'full wetness'}[k] for k in _land)))


def gcm_soil(folder: Path) -> dict:
    """PlaSim's soil capacity and free-evaporation share for a GCM run: the defaults unless its
    landmod_namelist sets them."""
    import re
    values = dict(PLASIM_SOIL)
    path = folder / 'landmod_namelist'
    text = path.read_text().lower() if path.exists() else ''
    for key in values:
        found = re.search(rf'\b{key}\s*=\s*([-+0-9.e]+)', text)
        if found:
            values[key] = float(found.group(1))
    return values


def gcm_equator_profile(run=GCM_RUN, years=GCM_YEARS, band=EQUATOR_BAND_DEG, latitude=0.0):
    """The GCM's reference for a ring on the equator, or on another circle of latitude: mean potential
    temperature, vapour and eastward wind of the sea columns within `band` degrees of the latitude against
    height above the sea (their surface), the sea surface and land ground temperatures, the air near the
    ground over the sea, and how readily the land there evaporates."""
    import numpy as np
    import netCDF4
    from climate.gcm.compare import area_weights
    gravity, rd, kappa = planet()['gravity_m_s2'], 287.04, 287.04 / 1005.7
    folder = HERE.parent / 'gcm' / 'runs' / run / 'model'
    fields = {k: [] for k in ('ta', 'hus', 'ua', 'ps', 'ts', 'tas', 'mrso')}
    for year in range(years[0], years[1] + 1):
        with netCDF4.Dataset(folder / f'MOST.{year:05d}.nc') as d:
            for k in fields:
                fields[k].append(np.asarray(d[k][:], dtype=float))
            lat = np.asarray(d['lat'][:], dtype=float)
            sigma = np.asarray(d['lev'][:], dtype=float)
            lsm = np.asarray(d['lsm'][:], dtype=float)[0] > 0.5
    f = {k: np.concatenate(v) for k, v in fields.items()}
    height = np.loadtxt(folder.parent / 'inputs' / 'moon_topography.sra', skiprows=1).ravel().reshape(lsm.shape) / gravity
    rows = np.abs(lat - latitude) < band
    sea = rows[:, None] & ~lsm & (height < 1.0)                           # seas at sea level; lakes stand higher
    land = rows[:, None] & lsm
    w = area_weights(lat, 64)
    mean = lambda x, m: float((x * w)[m].sum() / w[m].sum())
    ps = f['ps'] * 100.0                                                  # hPa in the output
    tv = f['ta'] * (1 + 0.608 * f['hus'])
    z = np.zeros_like(f['ta'])                                            # heights above the ground
    z[:, -1] = rd * tv[:, -1] / gravity * np.log(1 / sigma[-1])
    for k in range(len(sigma) - 2, -1, -1):
        z[:, k] = z[:, k + 1] + rd * 0.5 * (tv[:, k] + tv[:, k + 1]) / gravity * np.log(sigma[k + 1] / sigma[k])
    theta = f['ta'] * (1e5 / (sigma[None, :, None, None] * ps[:, None])) ** kappa
    tm = lambda x: x.mean(axis=0)                                         # time mean, (lev, lat, lon)
    profile = dict(z_m=[mean(tm(z)[k], sea) for k in range(len(sigma))],
                   theta_k=[mean(tm(theta)[k], sea) for k in range(len(sigma))],
                   qv_kg_kg=[mean(tm(f['hus'])[k], sea) for k in range(len(sigma))],
                   u_m_s=[mean(tm(f['ua'])[k], rows[:, None] & np.ones_like(lsm)) for k in range(len(sigma))])
    # Land stands higher than the ring's flat ground at sea level: carry its ground temperature down at
    # the lapse rate of the GCM's lowest kilometres.
    lapse = (profile['theta_k'][-2] - profile['theta_k'][-1]) / (profile['z_m'][-2] - profile['z_m'][-1])
    lapse = gravity / 1005.7 - lapse * (1e5 / mean(tm(ps), sea)) ** -kappa
    land_sea_level = tm(f['ts']) + lapse * height
    soil = gcm_soil(folder)
    wetness = np.minimum(1.0, f['mrso'] / (soil['drhsfull'] * soil['wsmax']))
    return dict(profile=profile, sea_surface_k=mean(tm(f['ts']), sea), land_ground_k=mean(land_sea_level, land),
                land_ground_in_place_k=mean(tm(f['ts']), land), land_height_m=mean(height, land), lapse_k_m=lapse,
                air_over_sea_k=mean(tm(f['tas']), sea), surface_pa_sea=mean(tm(ps), sea), sea_columns=int(sea.sum()),
                land_wetness=mean(tm(wetness), land), soil=soil,
                run=run, years=list(years), band_deg=band, latitude_deg=latitude, rows_deg=[float(v) for v in lat[rows]])


def vertical_grid(dz0, z_fine, stretch, dz_max, ztop):
    """Heights (m) of the w levels: dz0 up to z_fine, growing by `stretch` per level to dz_max, then even."""
    levels, z, dz = [0.0], 0.0, dz0
    while z < ztop - 1e-6:
        if z >= z_fine:
            dz = min(dz * stretch, dz_max)
        z = min(z + dz, ztop)
        levels.append(z)
    if levels[-1] - levels[-2] < 0.5 * dz_max:                            # no sliver at the top
        levels.pop(-2)
    return levels


def ring_surface(nx, dx, sea_k, land_k, band=1.0, latitude=0.0):
    """Land and water segments along a ring on a circle of latitude (x eastward from 0 E): water where the
    28% scenario's seas and rain-fed lakes cover at least half of the band within `band` degrees of the
    latitude. A ring column narrower than the product's 0.25-degree columns takes the one under its centre;
    a wider one, the mean of those whose centres it holds."""
    import numpy as np
    from climate.gcm import boundary
    product = boundary.lakes_product(*boundary.LAKES)
    lat, lon, water = product['lat_deg'], product['lon_deg'], product['water_fraction']
    along = water[np.abs(lat - latitude) < band].mean(axis=0)             # 0.25-degree columns
    if nx >= lon.size:
        centres = (np.arange(nx) + 0.5) * 360.0 / nx
        share = along[np.minimum((centres / (360.0 / lon.size)).astype(int), lon.size - 1)]
    else:
        holder = np.minimum((lon / (360.0 / nx)).astype(int), nx - 1)
        share = np.bincount(holder, weights=along, minlength=nx) / np.bincount(holder, minlength=nx)
    wet = share >= 0.5
    segments, start = [], 0
    for i in range(1, nx + 1):
        if i == nx or wet[i] != wet[start]:
            x0 = start * dx if start > 0 else -1e9
            x1 = i * dx if i < nx else 1e9
            segments.append((x0, x1, 2 if wet[start] else 1, 16 if wet[start] else PLACEHOLDER_LAND[0],
                             sea_k if wet[start] else land_k, sea_k if wet[start] else land_k))
            start = i
    return segments, float(wet.mean()), product['sha256']


def tilted_path(nx: int, tilt_deg: float, node_deg: float) -> dict:
    """Latitude, longitude and the longitude gained since x = 0 (0 to 360, degrees) at the column centres of a
    ring along a great circle tilted tilt_deg to the equator, x eastward from its northward crossing at node_deg E."""
    import numpy as np
    s = 2.0 * np.pi * (np.arange(nx) + 0.5) / nx
    tilt = np.radians(tilt_deg)
    lat = np.degrees(np.arcsin(np.sin(tilt) * np.sin(s)))
    gained = np.degrees(np.unwrap(np.arctan2(np.cos(tilt) * np.sin(s), np.cos(s))))
    return dict(lat=lat, lon=(node_deg + gained) % 360.0, dlon=gained)


def wetness_class(wetness) -> int:
    """The land-use row whose wetness is nearest."""
    return min(WETNESS_CLASSES, key=lambda c: abs(c[1] - float(wetness)))[0]


def landuse_rows(classes, land=None) -> dict:
    """Placeholder-land rows of the land-use table for (index, wetness) pairs. `land` may override the albedo
    (albedo_percent), roughness length (roughness_cm), soil thermal inertia (therin, as the table's THERIN) and
    wetness (for every row)."""
    base = PLACEHOLDER_LAND[1]
    rows = {index: base.replace('28,', f'{index},', 1).replace('.50', f'{w:.2f}'.lstrip('0') if w < 1 else f'{w:.2f}', 1)
            .replace("'Terluna placeholder land'", f"'Terluna placeholder land, wetness {w:.2f}'") for index, w in classes}
    if not land:
        return rows
    columns = dict(albedo_percent=1, wetness=2, roughness_cm=4, therin=5)       # after the index: ALBD SLMO SFEM SFZ0 THERIN
    labels = dict(albedo_percent='albedo {:g}%', roughness_cm='roughness {:g} cm', therin='THERIN {:g}')
    for index, w in classes:
        fields = rows[index].split(',', 7)
        for key, value in land.items():
            fields[columns[key]] = f'   {value:.2f}' if key in ('wetness', 'therin') else f'   {value:.0f}.'
        extra = ''.join(', ' + labels[k].format(v) for k, v in land.items() if k in labels)
        fields[7] = fields[7].split(',')[0] + f",'Terluna placeholder land, wetness {land.get('wetness', w):.2f}{extra}'\n"
        rows[index] = ','.join(fields)
    return rows


def gcm_path_reference(lat_path, lon_path, tilt_deg, zw, run=GCM_RUN, years=GCM_YEARS) -> dict:
    """The GCM's reference along a tilted ring. One upper-air profile for the whole path: the mean over its
    columns of the sea-column profiles of the GCM rows nearest them. The wind along the ring: each row's
    zonal-mean eastward wind times the cosine between east and the ring's heading, cos(tilt)/cos(latitude),
    averaged over the columns. For each column: the sea surface temperature of its latitude (rows interpolated),
    the ground temperature and wetness of the nearest GCM land cell (ground carried to sea level at its row's lapse
    rate), and the GCM's mean vertical wind at its latitude on the model's w levels zw."""
    import numpy as np
    import netCDF4
    gravity, rd, kappa = planet()['gravity_m_s2'], 287.04, 287.04 / 1005.7
    folder = HERE.parent / 'gcm' / 'runs' / run / 'model'
    soil = gcm_soil(folder)
    acc, count, flux = {}, 0, []
    for year in range(years[0], years[1] + 1):
        with netCDF4.Dataset(folder / f'MOST.{year:05d}.nc') as d:
            f = {k: np.asarray(d[k][:], dtype=float) for k in ('ta', 'hus', 'ua', 'va', 'ps', 'ts', 'tas', 'mrso')}
            lat = np.asarray(d['lat'][:], float)
            lon = np.asarray(d['lon'][:], float)
            sigma = np.asarray(d['lev'][:], float)
            lsm = np.asarray(d['lsm'][:], float)[0] > 0.5
        ps = f['ps'] * 100.0
        pieces = dict(theta=f['ta'] * (1e5 / (sigma[None, :, None, None] * ps[:, None])) ** kappa, tv=f['ta'] * (1 + 0.608 * f['hus']),
                      qv=f['hus'], u=f['ua'], ps=ps, ts=f['ts'], tas=f['tas'],
                      wet=np.minimum(1.0, f['mrso'] / (soil['drhsfull'] * soil['wsmax'])))
        for k, v in pieces.items():
            acc[k] = acc.get(k, 0.0) + v.sum(axis=0)
        count += f['ta'].shape[0]
        flux.append((ps[:, None] * f['va']).mean(axis=(0, 3)))
    m = {k: v / count for k, v in acc.items()}
    flux = np.mean(flux, axis=0)
    height = np.loadtxt(folder.parent / 'inputs' / 'moon_topography.sra', skiprows=1).ravel().reshape(lsm.shape) / gravity
    sea = ~lsm & (height < 1.0)
    edges = plasim_half_levels(sigma)
    rows = {}
    for j in range(1, lat.size - 1):
        if abs(lat[j]) > tilt_deg + 8.0:
            continue
        s_cells = sea[j] if sea[j].any() else np.ones(lon.size, bool)
        tv = m['tv'][:, j, s_cells].mean(axis=1)
        z = np.zeros(sigma.size)                                              # layer heights above the sea, bottom last
        z[-1] = rd * tv[-1] / gravity * np.log(1 / sigma[-1])
        for k in range(sigma.size - 2, -1, -1):
            z[k] = z[k + 1] + rd * 0.5 * (tv[k] + tv[k + 1]) / gravity * np.log(sigma[k + 1] / sigma[k])
        theta = m['theta'][:, j, s_cells].mean(axis=1)
        lapse = gravity / 1005.7 - (theta[-2] - theta[-1]) / (z[-2] - z[-1]) * (1e5 / m['ps'][j, s_cells].mean()) ** -kappa
        mass, _ = band_vertical_wind(lat, edges, flux, lat == lat[j], planet()['radius_m'], gravity)
        tvi = m['tv'][:, j].mean(axis=1)
        zi = np.zeros(edges.size)
        for k in range(sigma.size - 1, -1, -1):
            if edges[k] > 0:
                zi[k] = zi[k + 1] + rd * tvi[k] / gravity * np.log(edges[k + 1] / edges[k])
        rho = m['ps'][j].mean() * edges / (rd * np.interp(edges, sigma, m['tv'][:, j].mean(axis=1)))
        keep = edges > 0
        w = np.where(keep, mass / np.maximum(rho, 1e-9), 0.0)
        rows[float(lat[j])] = dict(z=z, theta=theta, qv=m['qv'][:, j, s_cells].mean(axis=1), u=m['u'][:, j].mean(axis=1),
                                   sea_k=float(m['ts'][j, sea[j]].mean()) if sea[j].any() else np.nan,
                                   air_k=float(m['tas'][j, s_cells].mean()), ps=float(m['ps'][j, s_cells].mean()), lapse=lapse,
                                   w_z=zi[keep][::-1], w=w[keep][::-1])
    lats = np.array(sorted(rows))
    pick = lambda key: np.array([rows[la][key] for la in lats])
    lat_path, lon_path = np.asarray(lat_path, float), np.asarray(lon_path, float)
    nearest = lats[np.abs(lat_path[:, None] - lats[None, :]).argmin(axis=1)]
    weights = {la: float((nearest == la).mean()) for la in lats}
    mean_row = lambda key: sum(weights[la] * rows[la][key] for la in lats)
    heading = np.cos(np.radians(tilt_deg)) / np.cos(np.radians(lat_path))           # cosine between east and the ring
    u_rows = pick('u')                                                                # (rows, levels)
    u_along = np.mean([np.array([np.interp(la, lats, u_rows[:, k]) for k in range(sigma.size)]) * h
                       for la, h in zip(lat_path, heading)], axis=0)
    sea_rows = pick('sea_k')
    ok = np.isfinite(sea_rows)
    sea_k = np.interp(lat_path, lats[ok], sea_rows[ok])
    # nearest GCM land cell for ground and wetness (great-circle nearest among land cells)
    glat, glon = np.meshgrid(lat, lon, indexing='ij')
    land_lat, land_lon = glat[lsm], glon[lsm]
    lapse_cells = np.interp(glat, lats, pick('lapse'))
    land_ground = (m['ts'] + lapse_cells * height)[lsm]
    land_wet = m['wet'][lsm]
    def nearest_land(la, lo):
        a, b = np.radians(la), np.radians(land_lat)
        cosd = np.sin(a) * np.sin(b) + np.cos(a) * np.cos(b) * np.cos(np.radians(lo - land_lon))
        return int(np.argmax(cosd))
    near = np.array([nearest_land(la, lo) for la, lo in zip(lat_path, lon_path)])
    zw = np.asarray(zw, float)
    w_rows = np.array([np.interp(zw, rows[la]['w_z'], rows[la]['w'], left=0.0, right=0.0) for la in lats])
    w_columns = np.array([[np.interp(la, lats, w_rows[:, k]) for k in range(zw.size)] for la in lat_path])
    profile = dict(z_m=list(mean_row('z')), theta_k=list(mean_row('theta')), qv_kg_kg=list(mean_row('qv')), u_m_s=list(u_along))
    return dict(profile=profile, sea_surface_k=float(sea_k.mean()), land_ground_k=float(land_ground[near].mean()),
                air_over_sea_k=float(mean_row('air_k')), surface_pa_sea=float(mean_row('ps')),
                sea_k=sea_k, land_k=land_ground[near], wetness=land_wet[near], w_columns=w_columns,
                land_wetness=float(land_wet[near].mean()), soil=soil, rows_deg=[float(v) for v in lats],
                row_weights=weights, run=run, years=list(years), tilt_deg=tilt_deg)


def path_surface(lat, lon, dx, sea_k, land_k, wetness, half_width_deg=0.5):
    """Land and water column by column along a path: water where the 28% scenario's seas and rain-fed lakes cover
    at least half of a box half_width_deg of latitude (and as far east and west) around the column's centre. One
    segment per column, with its own sea or land temperature and, on land, the land-use row of its wetness."""
    import numpy as np
    from climate.gcm import boundary
    product = boundary.lakes_product(*boundary.LAKES)
    plat, plon, water = product['lat_deg'], product['lon_deg'], product['water_fraction']
    nx = len(lat)
    share = np.zeros(nx)
    for i, (la, lo) in enumerate(zip(lat, lon)):
        r = np.abs(plat - la) <= half_width_deg
        c = np.abs((plon - lo + 180.0) % 360.0 - 180.0) * np.cos(np.radians(la)) <= half_width_deg
        share[i] = water[np.ix_(r, c)].mean()
    wet = share >= 0.5
    segments = []
    for i in range(nx):
        x0 = i * dx if i > 0 else -1e9
        x1 = (i + 1) * dx if i < nx - 1 else 1e9
        if wet[i]:
            segments.append((x0, x1, 2, 16, float(sea_k[i]), float(sea_k[i])))
        else:
            segments.append((x0, x1, 1, wetness_class(wetness[i]), float(land_k[i]), float(land_k[i])))
    return segments, float(wet.mean()), product['sha256']


def box_site(cfg: dict, zw, nx_ring: int) -> dict:
    """What a box takes from the ring it extends: the ring's GCM reference (upper air and wind along the ring), and
    over the ring's columns within the site's radius of it, the mean land wetness and ground temperature of the land
    columns and the mean GCM vertical wind on the w levels zw."""
    import numpy as np
    site = cfg['site']
    ring = CASES[site['ring']]
    track = tilted_path(nx_ring, ring['tilt_deg'], ring['node_deg'])
    ref = gcm_path_reference(track['lat'], track['lon'], ring['tilt_deg'], zw)
    a, b = np.radians(site['lat_deg']), np.radians(track['lat'])
    cosd = np.sin(a) * np.sin(b) + np.cos(a) * np.cos(b) * np.cos(np.radians(track['lon'] - site['lon_deg']))
    centre = int(np.argmax(cosd))
    length = 2.0 * math.pi * planet()['radius_m']
    reach = int(site['radius_m'] // (length / nx_ring))
    cols = (centre + np.arange(-reach, reach + 1)) % nx_ring
    segments, _, surface_sha = path_surface(track['lat'][cols], track['lon'][cols], 1.0, ref['sea_k'][cols],
                                            ref['land_k'][cols], ref['wetness'][cols])
    land = np.array([s[2] == 1 for s in segments])
    if not land.any():
        raise ValueError(f'no land within {site["radius_m"] / 1000:.0f} km of the site on {site["ring"]}')
    return dict(ref=ref, columns=[int(c) for c in cols], land_share=float(land.mean()), surface_product=surface_sha,
                wetness=float(ref['wetness'][cols][land].mean()), ground_k=float(ref['land_k'][cols][land].mean()),
                w_m_s=[float(v) for v in np.asarray(ref['w_columns'])[cols].mean(axis=0)])


def write_day_night(path: Path, forcing: dict) -> None:
    """terluna_lsadv.txt: the local-time bin and level counts, the level heights (m), then per bin (evenly spaced in
    hour angle from -180 degrees) a line of theta tendencies (K/s) and a line of vapour tendencies (kg/kg/s)."""
    import numpy as np
    z, th, qv = (np.asarray(forcing[k], float) for k in ('z_m', 'theta_k_s', 'qv_kg_kg_s'))
    path.write_text(f'{th.shape[0]} {z.size}\n' + ' '.join(f'{v:.3f}' for v in z) + '\n'
                    + ''.join(' '.join(f'{v:.6e}' for v in a) + '\n' + ' '.join(f'{v:.6e}' for v in b) + '\n'
                              for a, b in zip(th, qv)))


def write_vertical_wind_2d(path: Path, zw, w_columns) -> None:
    """terluna_wls2d.txt: the column and level counts, the level heights (m), then one line per column of its
    vertical wind (m/s), columns evenly spaced from x = 0."""
    import numpy as np
    w_columns = np.asarray(w_columns, float)
    path.write_text(f'{w_columns.shape[0]} {len(zw)}\n' + ' '.join(f'{z:.3f}' for z in zw) + '\n'
                    + ''.join(' '.join(f'{v:.6e}' for v in row) + '\n' for row in w_columns))


def set_namelist(text: str, section: str, key: str, value) -> str:
    """Set one entry of a Fortran namelist held as text; the entry must exist in its section."""
    import re
    start = text.index(f'&{section}\n')
    end = text.index('\n /', start)
    block = text[start:end]
    pattern = re.compile(rf'^(\s*{re.escape(key)}\s*=\s*)([^,\n]*)(,?)', re.M)
    if not pattern.search(block):
        raise KeyError(f'{key} not in &{section}')
    if isinstance(value, bool):
        value = '.true.' if value else '.false.'
    block = pattern.sub(lambda m: f'{m.group(1)}{value}{m.group(3)}', block, count=1)
    return text[:start] + block + text[end:]


def time_scale(cfg, gravity) -> float:
    """Factor on the case's lengths and times: 1 at lunar gravity, g_moon/g for a pair case at gravity g."""
    return planet()['gravity_m_s2'] / gravity if cfg.get('kind') == 'pair' else 1.0


def case_settings(cfg, nx, dx, nz, ztop, air, gravity):
    """Namelist entries of a case, by section, on top of CM1's RCE template."""
    length = nx * dx
    s = time_scale(cfg, gravity)
    pair = cfg.get('kind') == 'pair'
    box = cfg.get('kind') == 'box'                                       # 3-D, one site: the Sun the same across it
    latitude = cfg['site']['lat_deg'] if box else cfg.get('latitude_deg', 0.0)
    tilted = cfg.get('kind') == 'tilted'
    turning = turning_rate(cfg, latitude)
    rotating = turning != 0.0
    return {
        'param0': dict(nx=nx, ny=cfg['ny'] if box else 1, nz=nz, ppnode=cfg.get('ranks', 8), timeformat=3, timestats=1),
        'param1': dict(dx=round(dx, 3), dy=round(dx, 3), dz=round(ztop / nz, 1), dtl=round(40.0 * s, 3), cfl_limit=1.0,
                       timax=round(cfg['days'] * 86400.0 * s), run_time=-999.9, tapfrq=round(cfg['output_s'] * s, 3),
                       rstfrq=round(cfg['restart_s'] * s, 3), statfrq=round(3600.0 * s, 3), prclfrq=1.0e9),
        'param2': dict(cm1setup=2, testcase=0, adapt_dt=1, irst=0, rstnum=1, ipbl=2, sgsmodel=0, tconfig=2,
                       horizturb=0, irdamp=2, psolver=3, ptype=5, ihail=0, icor=int(rotating), betaplane=int(tilted and rotating),
                       lspgrad=int(rotating and not tilted),
                       eqtset=2, idiss=1, wbc=1, ebc=1, sbc=1, nbc=1, bbc=3, tbc=1, isnd=7, iwnd=0, itern=0, iinit=0,
                       irandp=1, iorigin=1, apmasscon=1),
        # The Rayleigh layer's time and the PBL scheme's asymptotic length scale are Earth's values stretched
        # by the ratio of gravities, as the dynamics stretch.
        'param3': dict(fcor=float(f'{turning:.6e}'), rdalpha=round(gravity / EARTH_G / 300.0, 8), zd=0.77 * ztop,
                       l_inf=round(75.0 * EARTH_G / gravity, 1), ndcnst=cfg['droplets_cm3']),
        'param6': dict(stretch_z=4, ztop=ztop),
        'param8': dict(**(dict(var11=1.0 if cfg.get('day_night') else 0.0) if box else {}),
                       var12=cfg['tilt_deg'] if tilted else 0.0,
                       var13=(2.0 if tilted else 1.0) if cfg.get('large_scale_w') else 0.0,
                       var14=cfg['nudge']['z_start_m'] * s, var15=cfg['nudge']['z_full_m'] * s,
                       var16=0.0 if pair or box else round(length, 3), var17=0.0 if pair else round(air['sunlight_w_m2'], 3),
                       var18=0.0 if pair else cfg['start_hour_angle_deg'] + (cfg['node_deg'] if tilted else
                                                                             cfg['site']['lon_deg'] if box else 0.0),
                       var19=0.0 if pair else round(solar_day_s(), 1)),
        'param9': dict(output_format=1, output_filetype=2, output_sfcparams=int(bool(cfg.get('surface_output'))),
                       output_tke=0, output_km=0,
                       output_kh=0, output_uinterp=1, output_vinterp=int(rotating or box), output_v=0, output_winterp=1,
                       output_radten=int(not box), output_cape=1, output_cin=1, output_lcl=1, output_lfc=1, output_pwat=1,
                       output_lwp=1, **(dict(output_u=0, output_w=0, output_dbz=0) if box else {})),
        'param11': dict(radopt=0 if pair else 2, dtrad=1800.0, ctrlat=latitude, ctrlon=0.0, year=2014),
        'param12': dict(isfcflx=1, sfcmodel=2, oceanmodel=1, initsfc=1 if pair else 9, season=1),
        'param14': dict(dodomaindiag=True, diagfrq=round(cfg['output_s'] * s, 3)),
        'param16': dict(restart_format=1, restart_filetype=2, restart_reset_frqtim=True),
        'param19': dict(do_lsnudge=True, do_lsnudge_u=True, do_lsnudge_v=rotating or box, do_lsnudge_th=True,
                        do_lsnudge_qv=True, lsnudge_tau=round(cfg['nudge']['tau_s'] * s, 3), lsnudge_start=1.0,
                        lsnudge_end=1.0e12, lsnudge_ramp_time=round(cfg['nudge']['ramp_s'] * s, 3)),
    }


def coriolis(latitude_deg: float) -> float:
    """The Coriolis parameter (1/s) at a latitude, from the Moon's sidereal rotation."""
    return 2.0 * planet()['rotation_rate_rad_s'] * math.sin(math.radians(latitude_deg))


def turning_rate(cfg: dict, latitude: float) -> float:
    """CM1's fcor for a case: the Coriolis parameter of its latitude; on a tilted ring twice the rotation rate,
    which CM1 multiplies by sin(latitude) column by column; none on a great circle of no tilt, the equator."""
    if cfg.get('kind') == 'tilted':
        return coriolis(90.0) if cfg['tilt_deg'] > 0.0 else 0.0
    return coriolis(latitude)


def ring_grid(latitude_deg: float, dx_target_m: float, ranks: int):
    """Columns, their width and the length of a ring around the Moon on a circle of latitude; the column
    count a multiple of the ranks."""
    length = 2.0 * math.pi * planet()['radius_m'] * math.cos(math.radians(latitude_deg))
    nx = int(round(length / dx_target_m / ranks)) * ranks
    return nx, length / nx, length


def land_moisture(cfg: dict, ref: dict) -> float:
    """The land's moisture availability: the case's own value, or the GCM's for 'gcm'."""
    value = cfg.get('land_moisture', 0.5)
    return float(ref['land_wetness']) if value == 'gcm' else float(value)


def solar_day_s():
    return planet()['solar_day_days'] * 86400.0


def water_mean_profile(snapshots, land) -> dict:
    """Mean potential temperature and vapour by level over the water columns of some snapshots, with the
    mean surface pressure, 2 m temperature and 2 m vapour there."""
    import numpy as np
    sea = ~np.asarray(land, bool)
    total, count = {}, 0
    for d in snapshots:
        for key, value in (('theta_k', d['th'][:, sea].mean(axis=1)), ('qv_kg_kg', d['qv'][:, sea].mean(axis=1)),
                           ('surface_pa', d['psfc'][sea].mean()), ('air_2m_k', d['t2'][sea].mean()),
                           ('qv_2m_kg_kg', d['q2'][sea].mean())):
            total[key] = total.get(key, 0.0) + np.asarray(value, dtype=float)
        count += 1
    return {key: value / count for key, value in total.items()}


def ring_sea_profile(days: float = 3.0, name: str = 'ring') -> dict:
    """The ring's own air over its seas, averaged over its last `days`: where the gravity pair starts."""
    from climate.crm import ring_analysis as ra
    case = RUNS / name
    geo = ra.case_geometry(case)
    tap = geo['record']['configuration']['output_s']
    outputs = sorted(int(p.name[8:14]) for p in case.glob('cm1out_t*_s.dat'))
    last = (outputs[-1] - 1) * tap
    use = [n for n in outputs if (n - 1) * tap >= last - days * 86400.0]
    profile = water_mean_profile((ra.read_snapshot(case, n) for n in use), geo['land'])
    profile.update(z_m=geo['zh'], span_days=[(use[0] - 1) * tap / 86400.0, last / 86400.0], case=name)
    return profile


def setup(name: str) -> Path:
    """Write a case's inputs into climate/crm/runs/<name>: namelist, sounding, grid, nudging profile,
    surface segments, land-use table and links to the executable and radiation tables."""
    import numpy as np
    cfg = CASES[name]
    tree = fetch()
    exe = CM1_HOME / 'build' / cfg['build'] / 'cm1.exe'
    if not exe.exists():
        raise RuntimeError(f'build {cfg["build"]} first')
    build_record = json.loads((exe.parent / 'build.json').read_text())
    air, gravity = build_record['air'], build_record['gravity_m_s2']
    case = RUNS / name
    if (case / 'progress.json').exists():
        raise RuntimeError(f'{case} has started; remove it to set up afresh')
    case.mkdir(parents=True, exist_ok=True)
    box = cfg.get('kind') == 'box'
    latitude = cfg['site']['lat_deg'] if box else cfg.get('latitude_deg', 0.0)
    tilted = cfg.get('kind') == 'tilted'
    s = time_scale(cfg, gravity)
    pair = cfg.get('kind') == 'pair'
    if pair:
        nx, dx = cfg['nx'], cfg['dx_moon_m'] * s
        length = nx * dx
    elif box:                                                             # the spacing of the ring it extends
        nx_ring, dx, _ = ring_grid(0.0, cfg['dx_target_m'], cfg['ranks'])
        nx, length = cfg['nx'], cfg['nx'] * dx
    else:
        nx, dx, length = ring_grid(latitude, cfg['dx_target_m'], cfg['ranks'])
    zw = [z * s for z in vertical_grid(**cfg['z_grid'])]
    nz = len(zw) - 1
    track = tilted_path(nx, cfg['tilt_deg'], cfg['node_deg']) if tilted else None
    site = box_site(cfg, zw, nx_ring) if box else None
    ref = (gcm_path_reference(track['lat'], track['lon'], cfg['tilt_deg'], zw) if tilted else site['ref'] if box
           else gcm_equator_profile(band=cfg.get('band_deg', EQUATOR_BAND_DEG), latitude=latitude))
    zh = 0.5 * (np.array(zw[1:]) + np.array(zw[:-1]))
    prof = ref['profile']
    zg = np.array(prof['z_m'][::-1]) * s                                  # GCM levels, bottom first
    th = np.array(prof['theta_k'][::-1])
    qv = np.array(prof['qv_kg_kg'][::-1])
    uu = np.array(prof['u_m_s'][::-1]) * (0.0 if pair else 1.0)
    top_slope = (th[-1] - th[-2]) / (zg[-1] - zg[-2])
    theta_at = lambda z: np.where(z <= zg[-1], np.interp(z, zg, th), th[-1] + top_slope * (z - zg[-1]))
    qv_at = lambda z: np.exp(np.interp(z, zg, np.log(np.maximum(qv, 1e-9))))
    u_at = lambda z: np.interp(z, zg, uu)
    kappa = 287.04 / 1005.7
    theta_sfc = ref['air_over_sea_k'] * (1e5 / ref['surface_pa_sea']) ** kappa
    lines = [f"{ref['surface_pa_sea'] / 100:12.4f} {theta_sfc:12.4f} {qv[0] * 1000:12.5f}"]
    start_theta, start_qv, initial = theta_at, qv_at, dict(source='GCM reference')
    if cfg.get('initial') == 'ring_sea':
        start = ring_sea_profile()
        zr = np.asarray(start['z_m']) * s
        start_theta = lambda z: np.interp(z, zr, start['theta_k'])
        start_qv = lambda z: np.exp(np.interp(z, zr, np.log(np.maximum(start['qv_kg_kg'], 1e-12))))
        theta_2m = float(start['air_2m_k']) * (1e5 / float(start['surface_pa'])) ** kappa
        lines = [f"{float(start['surface_pa']) / 100:12.4f} {theta_2m:12.4f} {float(start['qv_2m_kg_kg']) * 1000:12.5f}"]
        initial = dict(source='the ring over its seas', case=start['case'], span_days=start['span_days'])
    # the sounding ends at the model top, above the highest scalar level, so rounding cannot leave it short
    lines += [f'{z:12.3f} {float(start_theta(z)):12.4f} {float(start_qv(z)) * 1000:12.5f} {float(u_at(z)):8.3f} {0.0:8.3f}'
              for z in [*zh, zw[-1]]]
    (case / 'input_sounding').write_text('\n'.join(lines) + '\n')
    (case / 'input_grid_z').write_text(''.join(f'{z:.3f}\n' for z in zw))
    nudge = [' *  Header:   lsnudge_time1 (s)   lsnudge_time2 (s)', '                     0.0                1.0e30 ',
             ' *  Profile.   Note: values will be ignored if nudging for that variable is off.',
             ' *    z (m)   theta (K)    qv (g/kg)  u (m/s, grnd-reltv) v (m/s, grnd-reltv)']
    nudge += [f'  {z:12.4f} {float(theta_at(z)):12.4f} {float(qv_at(z)) * 1000:12.6f} {float(u_at(z)):10.4f} {0.0:10.4f}'
              for z in [0.0, *zh, zw[-1]]]
    (case / 'lsnudge_0001.dat').write_text('\n'.join(nudge) + '\n')
    if pair:
        segments, water_share, surface_sha = [(-1e9, 1e9, 2, 16, ref['sea_surface_k'], ref['sea_surface_k'])], 1.0, None
    elif box:
        segments = [(-1e9, 1e9, 1, wetness_class(site['wetness']), site['ground_k'], site['ground_k'])]
        water_share, surface_sha = 0.0, site['surface_product']
    elif tilted:
        segments, water_share, surface_sha = path_surface(track['lat'], track['lon'], dx, ref['sea_k'], ref['land_k'], ref['wetness'])
    else:
        segments, water_share, surface_sha = ring_surface(nx, dx, ref['sea_surface_k'], ref['land_ground_k'],
                                                          latitude=latitude)
    (case / 'terluna_surface.txt').write_text(f'{len(segments)}\n' + ''.join(
        f'{x0:.1f} {x1:.1f} {xl:.1f} {lu:d} {tsk:.3f} {tmn:.3f}\n' for x0, x1, xl, lu, tsk, tmn in segments))
    table = (tree / 'run' / 'LANDUSE.TBL').read_text().splitlines(keepends=True)
    index, row = PLACEHOLDER_LAND
    table = [row.replace('.50', f"{land_moisture(cfg, ref):.2f}".lstrip('0')) if line.startswith(f'{index},') else line
             for line in table]
    if tilted or box:                                                     # land-use rows for the wetness classes
        rows = landuse_rows(WETNESS_CLASSES, cfg.get('land'))
        table = [next((r for i, r in rows.items() if line.startswith(f'{i},')), line) for line in table]
    (case / 'LANDUSE.TBL').write_text(''.join(table))
    vertical_wind = None
    if box and cfg.get('large_scale_w') == 'gcm':
        vertical_wind = dict(z_m=[float(z) for z in zw], w_m_s=site['w_m_s'], ring=cfg['site']['ring'],
                             ring_columns=site['columns'], rows_deg=ref['rows_deg'])
        write_vertical_wind(case / 'terluna_wls.txt', vertical_wind, zw[-1])
    elif tilted and cfg.get('large_scale_w') == 'gcm':
        write_vertical_wind_2d(case / 'terluna_wls2d.txt', zw, ref['w_columns'])
        wc = np.asarray(ref['w_columns'])
        vertical_wind = dict(z_m=[float(z) for z in zw], path_mean_m_s=[float(v) for v in wc.mean(axis=0)],
                             min_m_s=float(wc.min()), max_m_s=float(wc.max()), rows_deg=ref['rows_deg'])
    elif cfg.get('large_scale_w') == 'gcm':
        vertical_wind = gcm_vertical_wind(latitude, cfg.get('band_deg', EQUATOR_BAND_DEG))
        write_vertical_wind(case / 'terluna_wls.txt', vertical_wind, zw[-1])
    day_night = None
    if box and cfg.get('day_night'):                                      # the rings' day-night circulation there
        from climate.crm.ring_crossings import day_night_advection
        dn = cfg['day_night']
        forcing = day_night_advection(dn['rings'], cfg['site']['lat_deg'], cfg['site']['lon_deg'], dn['from_day'],
                                      dn['wavenumbers'], dn['bins'], dn['harmonics'], dn['top_m'], cfg['site']['radius_m'])
        write_day_night(case / 'terluna_lsadv.txt', forcing)
        low = forcing['z_m'] < 3000.0
        spread = lambda a, scale: dict(mean=float(a[:, low].mean() * scale), min=float(a[:, low].mean(axis=1).min() * scale),
                                       max=float(a[:, low].mean(axis=1).max() * scale))
        day_night = dict({k: forcing[k] for k in ('rings', 'wavenumbers', 'shortest_wave_km', 'harmonics', 'top_m', 'from_day')},
                         lowest_3_km=dict(theta_k_day=spread(forcing['theta_k_s'], 86400.0),
                                          qv_g_kg_day=spread(forcing['qv_kg_kg_s'], 8.64e7)))
    for link, target in (('cm1.exe', exe), ('RRTMG_LW_DATA', tree / 'run' / 'RRTMG_LW_DATA'),
                         ('RRTMG_SW_DATA', tree / 'run' / 'RRTMG_SW_DATA')):
        path = case / link
        if path.is_symlink() or path.exists():
            path.unlink()
        path.symlink_to(target)
    text = (tree / 'run' / 'config_files' / 'cpm_RadConvEquil' / 'namelist.input').read_text()
    settings = case_settings(cfg, nx, dx, nz, zw[-1], air, gravity)
    settings['param12'].update(tsk0=round(ref['sea_surface_k'], 3), tmn0=round(ref['land_ground_k'], 3), xland0=2.0, lu0=16)
    for section, entries in settings.items():
        for key, value in entries.items():
            text = set_namelist(text, section, key, value)
    (case / 'namelist.template').write_text(text)
    digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]
    if tilted or box:                                                     # per-column arrays as rounded lists
        ref = {k: ([round(float(x), 4) for x in v] if isinstance(v, np.ndarray) else v) for k, v in ref.items() if k != 'w_columns'}
    record = dict(case=name, configuration=cfg, time_scale=s,
                  grid=dict(nx=nx, dx_m=dx, length_m=length, nz=nz, ztop_m=zw[-1], **(dict(ny=cfg['ny']) if box else {})),
                  site=dict(cfg['site'], ring_columns=site['columns'], land_share=site['land_share'], wetness=site['wetness'],
                            wetness_class=wetness_class(site['wetness']), ground_k=site['ground_k']) if box else None,
                  build=build_record, reference=ref, initial=initial, water_share=water_share, surface_product=surface_sha,
                  path=dict(tilt_deg=cfg['tilt_deg'], node_deg=cfg['node_deg'], lat=[round(float(v), 4) for v in track['lat']],
                            lon=[round(float(v), 4) for v in track['lon']], dlon=[round(float(v), 4) for v in track['dlon']])
                  if tilted else None,
                  latitude_deg=latitude, coriolis_1_s=turning_rate(cfg, latitude),
                  land_moisture=None if pair else site['wetness'] if box else land_moisture(cfg, ref),
                  large_scale_w=vertical_wind,
                  day_night=day_night,
                  runner=digest(__file__),
                  inputs={p: digest(case / p) for p in ('input_sounding', 'input_grid_z', 'lsnudge_0001.dat',
                                                         'terluna_surface.txt', 'LANDUSE.TBL', 'namelist.template',
                                                         'terluna_wls.txt', 'terluna_wls2d.txt', 'terluna_lsadv.txt')
                          if (case / p).exists()})
    (case / 'case.json').write_text(json.dumps(record, indent=1) + '\n')
    return case


def latest_restart(case: Path):
    """Index of the newest complete restart set (cm1rst_tNNNNNN_{i,s,u,v,w,x}.dat), or 0."""
    import re
    found = {int(m.group(1)) for p in case.glob('cm1rst_t*_s.dat') if (m := re.match(r'cm1rst_t(\d+)_s\.dat$', p.name))}
    complete = [n for n in found if all((case / f'cm1rst_t{n:06d}_{part}.dat').exists() for part in 'isuvwx')]
    return max(complete, default=0)


def stats(case: Path) -> dict:
    """CM1's statistics file (cm1out_stats.dat) as arrays, one entry per statistics time."""
    import numpy as np
    lines = (case / 'cm1out_stats.ctl').read_text().splitlines()
    i = next(k for k, line in enumerate(lines) if line.lower().startswith('vars'))
    names = [line.split()[0] for line in lines[i + 1:i + 1 + int(lines[i].split()[1])]]
    raw = np.fromfile(case / 'cm1out_stats.dat', dtype='<f4')
    data = raw[:raw.size // len(names) * len(names)].reshape(-1, len(names))
    return {name: data[:, k] for k, name in enumerate(names)}


def _progress(case: Path) -> dict:
    path = case / 'progress.json'
    return json.loads(path.read_text()) if path.exists() else dict(segments=[])


def _save_progress(case: Path, progress: dict):
    tmp = case / 'progress.json.part'
    tmp.write_text(json.dumps(progress, indent=1) + '\n')
    tmp.replace(case / 'progress.json')


def _unlimited_stack():
    import resource
    resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))


def case_length_s(cfg: dict, scale: float) -> float:
    """Model seconds a case runs: its length in days, cut to the last whole restart interval, the last state a
    later run can resume from (two lunar days at half-day restarts end at day 59.0)."""
    restart_s = cfg['restart_s'] * scale
    return math.floor(cfg['days'] * 86400.0 * scale / restart_s + 1e-6) * restart_s


def run(name: str, hours: float, threads: int = 8) -> dict:
    """Run a case in segments of segment_s model seconds, each from the newest restart file, until the
    case's length, the wall-clock budget or a STOP file ends it."""
    cfg = CASES[name]
    case = RUNS / name
    record = json.loads((case / 'case.json').read_text())
    scale = record.get('time_scale', 1.0)                                # model seconds per lunar-equivalent second
    restart_s, segment_s = cfg['restart_s'] * scale, cfg['segment_s'] * scale
    total = case_length_s(cfg, scale)
    progress = _progress(case)
    (case / 'run.lock').write_text(str(os.getpid()))
    deadline = time.time() + hours * 3600.0
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_STACKSIZE='512M')
    try:
        while True:
            done = progress['segments'][-1]['model_s'] if progress['segments'] else 0.0
            if done >= total - 1.0:
                break
            if (case / 'STOP').exists():
                (case / 'STOP').unlink()
                break
            if progress['segments'] and time.time() + progress['segments'][-1]['wall_s'] > deadline:
                break
            rst = latest_restart(case)
            if rst and abs(rst * restart_s - done) > 1.0:
                raise RuntimeError(f'restart {rst} does not match the recorded model time {done} s')
            text = (case / 'namelist.template').read_text()
            text = set_namelist(text, 'param2', 'irst', 1 if rst else 0)
            text = set_namelist(text, 'param2', 'rstnum', max(rst, 1))
            text = set_namelist(text, 'param1', 'run_time', round(min(segment_s, total - done), 3))
            (case / 'namelist.input').write_text(text)
            log = case / f'cm1_segment_{len(progress["segments"]) + 1:03d}.log'
            start = time.time()
            with open(log, 'w') as out:
                done_run = subprocess.run(['nice', '-n', '5', './cm1.exe'], cwd=case, env=env, stdout=out,
                                          stderr=subprocess.STDOUT, preexec_fn=_unlimited_stack)
            if done_run.returncode != 0 or 'Program terminated normally' not in log.read_text()[-4000:]:
                raise RuntimeError(f'CM1 stopped abnormally; see {log}')
            if latest_restart(case) <= rst:
                raise RuntimeError(f'CM1 wrote no new restart file in {case}')
            reached = latest_restart(case) * restart_s
            progress['segments'].append(dict(segment=len(progress['segments']) + 1, model_s=reached,
                                             wall_s=round(time.time() - start, 1), log=log.name,
                                             finished=time.strftime('%Y-%m-%d %H:%M:%S')))
            progress['case'] = record['case']
            _save_progress(case, progress)
    finally:
        (case / 'run.lock').unlink(missing_ok=True)
    return progress


def status(name: str) -> str:
    case = RUNS / name
    progress = _progress(case)
    cfg = CASES[name]
    scale = json.loads((case / 'case.json').read_text()).get('time_scale', 1.0)
    done = (progress['segments'][-1]['model_s'] if progress['segments'] else 0.0) / scale   # lunar-equivalent
    total = case_length_s(cfg, scale) / scale
    wall = sum(s['wall_s'] for s in progress['segments'])
    rate = done * scale / wall if wall else float('nan')
    running = (case / 'run.lock').exists()
    return (f"{name}: {done / 86400:.2f} of {total / 86400:.2f} days ({100 * done / total:.0f}%), "
            f"{wall / 3600:.1f} wall hours, {rate:.0f} model s per wall s; {'running' if running else 'idle'}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest='command', required=True)
    b = sub.add_parser('build')
    b.add_argument('labels', nargs='*', default=sorted(BUILDS))
    b.add_argument('--jobs', type=int, default=8)
    s = sub.add_parser('setup')
    s.add_argument('case', choices=sorted(CASES))
    r = sub.add_parser('run')
    r.add_argument('case', choices=sorted(CASES))
    r.add_argument('--hours', type=float, default=4.0, help='wall-clock budget')
    r.add_argument('--threads', type=int, default=8)
    for command in ('status', 'stop'):
        sub.add_parser(command).add_argument('case', choices=sorted(CASES))
    args = parser.parse_args(argv)
    if args.command == 'run':
        run(args.case, args.hours, args.threads)
        print(status(args.case))
    elif args.command == 'status':
        print(status(args.case))
    elif args.command == 'stop':
        (RUNS / args.case / 'STOP').write_text('stop after the current segment\n')
        print(f'{args.case}: will stop after the current segment')
    if args.command == 'build':
        for label in args.labels:
            exe = build(label, args.jobs)
            print(f'{label}: {exe}')
    elif args.command == 'setup':
        case = setup(args.case)
        record = json.loads((case / 'case.json').read_text())
        g = record['grid']
        print(f"{case}: {g['nx']}{' x ' + str(g['ny']) if 'ny' in g else ''} x {g['nz']} points, dx {g['dx_m']:.1f} m, top {g['ztop_m'] / 1000:.0f} km, "
              f"water {record['water_share']:.2f} of the ring")
    return 0


if __name__ == '__main__':
    sys.exit(main())
