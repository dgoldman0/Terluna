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

from climate.crm import cm1_elec
from shared.constants import (CM1_GRAVITY, CM1_RAIN_BREAKUP_CUTOFF, CM1_RAIN_BREAKUP_RATE, CM1_RAIN_FALL_CAP,
                              CM1_RAIN_FALL_EXPONENT, CM1_RAIN_MEAN_DIAMETER_CAP)

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
# Terrain (terrain_flag, itern = 4, heights from perts.dat). CM1 refuses its large-scale nudging over terrain, where
# one model level lies at different heights in different columns. The terrain build nudges instead the domain mean of
# each column's departure from the reference at the column's own height (the reference given on the first column's
# levels, at sea level on the rings); without terrain it nudges as before. CM1 also switches its Rayleigh layer to
# damp toward the domain mean whenever it nudges, which it then refuses over terrain; the terrain build keeps the
# layer damping toward the initial state there. It refuses a large-scale vertical wind over terrain too; the terrain
# build allows it, advecting on the first column's metrics, which understate vertical gradients over ground at
# height z by the factor 1 - z/ztop (under 6% up to 8.3 km of a 150-km model).
TERRAIN_NUDGE = """        ! Terluna: nudging over terrain, each column at its own height: at each level, the mean of the columns'
        ! departures from the reference at their own heights, weighted by how fully each column's height there is
        ! held (above var14 m, fully above var15 m), applied to each column with its weight. The reference is
        ! given on the first column's levels; where each column's levels fall on it, and their weights, are found
        ! once and kept. Without terrain, the domain mean against the reference.
        do k=1,nk
          terluna_dev(k) = {mean}(k{q})-lsnudge_{v}(k,lsnudge_count)
        enddo
        if( terrain_flag )then
          if( .not. allocated(terluna_nk) )then
            allocate( terluna_nk(ni,nj,nk) , terluna_nf(ni,nj,nk) , terluna_nw(ni,nj,nk) , terluna_sw(nk) )
            terluna_sw = 0.0
            do j=1,nj
            do i=1,ni
              terluna_kz = 1
              do k=1,nk
                do while( terluna_kz.lt.nk-1 .and. zh(1,1,terluna_kz+1).lt.zh(i,j,k) )
                  terluna_kz = terluna_kz+1
                enddo
                terluna_nk(i,j,k) = terluna_kz
                terluna_nf(i,j,k) = max(0.0,min(1.0,(zh(i,j,k)-zh(1,1,terluna_kz))/(zh(1,1,terluna_kz+1)-zh(1,1,terluna_kz))))
                terluna_nw(i,j,k) = 1.0
                if( var15.gt.var14 ) terluna_nw(i,j,k) = min(1.0,max(0.0,(zh(i,j,k)-var14)/(var15-var14)))
                terluna_sw(k) = terluna_sw(k) + terluna_nw(i,j,k)
              enddo
            enddo
            enddo
          endif
          !$omp parallel do default(shared) private(i,j,k,terluna_kz)
          do k=1,nk
            terluna_dev(k) = 0.0
            if( terluna_sw(k).gt.0.0 )then
              do j=1,nj
              do i=1,ni
                terluna_kz = terluna_nk(i,j,k)
                terluna_dev(k) = terluna_dev(k) + terluna_nw(i,j,k)*( {field}                                  &
                    - ( lsnudge_{v}(terluna_kz,lsnudge_count)                                                   &
                      + terluna_nf(i,j,k)*(lsnudge_{v}(terluna_kz+1,lsnudge_count)-lsnudge_{v}(terluna_kz,lsnudge_count)) ) )
              enddo
              enddo
              terluna_dev(k) = terluna_dev(k)/terluna_sw(k)
              do j=1,nj
              do i=1,ni
                {ten} = {ten} - lsnudgefac*terluna_nw(i,j,k)*terluna_dev(k)/(lsnudge_tau)
              enddo
              enddo
            endif
          enddo
          terluna_dev = 0.0   ! applied above, column by column; the domain-wide step below adds nothing
        endif
"""
NUDGE_TH = "        if( myid.eq.0 ) print *,'  applying large-scale nudging of theta, lsnudge_count = ',lsnudge_count\n        rdt = 1.0/dt\n"
NUDGE_QV = "        if( myid.eq.0 ) print *,'  applying large-scale nudging of qv, lsnudge_count = ',lsnudge_count\n        rdt = 1.0/dt\n"
SURFACE_2D_BLOCK = """
        ! Terluna: every column's ground from terluna_surface2d.txt when it exists (boxes over terrain): the
        ! column counts in x and y, then one line per column, x fastest: tsk (K), tmn (K), land-use index.
        ! Columns are found by position, so the halo takes the neighbours across the box's periodic edges.
        inquire(file='terluna_surface2d.txt',exist=terluna_two)
        if( terluna_two )then
          open(unit=97,file='terluna_surface2d.txt',status='old',action='read')
          read(97,*) terluna_nx,terluna_ny
          allocate( terluna_tsk(terluna_nx,terluna_ny) , terluna_tmn(terluna_nx,terluna_ny) , terluna_lu(terluna_nx,terluna_ny) )
          do terluna_j=1,terluna_ny
          do terluna_i=1,terluna_nx
            read(97,*) terluna_tsk(terluna_i,terluna_j),terluna_tmn(terluna_i,terluna_j),terluna_lu(terluna_i,terluna_j)
          enddo
          enddo
          close(unit=97)
          do j=jb,je
          do i=ib,ie
            tsk(i,j) = terluna_tsk(modulo(floor(xh(i)/dx),terluna_nx)+1,modulo(floor(yh(j)/dy),terluna_ny)+1)
          enddo
          enddo
          do j=jbl,jel
          do i=ibl,iel
            tmn(i,j) = terluna_tmn(modulo(floor(xh(i)/dx),terluna_nx)+1,modulo(floor(yh(j)/dy),terluna_ny)+1)
            lu_index(i,j) = terluna_lu(modulo(floor(xh(i)/dx),terluna_nx)+1,modulo(floor(yh(j)/dy),terluna_ny)+1)
          enddo
          enddo
          deallocate( terluna_tsk , terluna_tmn , terluna_lu )
          if( myid.eq.0 ) print *,'  Terluna: ground temperature and land use for every column from terluna_surface2d.txt'
        endif
"""
ADVECTION_OLD = """      do k=1,nk
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
"""
ADVECTION_COLUMNS = """      ! Terluna: over terrain each column takes the forcing at its own height above sea level, its levels
      ! found once and kept; the corner, the lowest and the highest column's lowest forced level are printed
      if( .not. allocated(terluna_al) )then
        allocate( terluna_al(ni,nj,nk) , terluna_at(ni,nj,nk) )
        terluna_al = 0
        terluna_at = 0.0
        do j=1,nj
        do i=1,ni
        do k=1,nk
          if( zh(i,j,k).ge.terluna_az(1) .and. zh(i,j,k).le.terluna_az(terluna_nl) )then
            terluna_l = 1
            do terluna_m=1,terluna_nl-1
              if( zh(i,j,k).ge.terluna_az(terluna_m) ) terluna_l = terluna_m
            enddo
            terluna_al(i,j,k) = terluna_l
            terluna_at(i,j,k) = (zh(i,j,k)-terluna_az(terluna_l))/(terluna_az(terluna_l+1)-terluna_az(terluna_l))
          endif
        enddo
        enddo
        enddo
        terluna_ii = 1
        terluna_jj = 1
        terluna_iq = 1
        terluna_jq = 1
        do j=1,nj
        do i=1,ni
          if( zh(i,j,1).lt.zh(terluna_ii,terluna_jj,1) )then
            terluna_ii = i
            terluna_jj = j
          endif
          if( zh(i,j,1).gt.zh(terluna_iq,terluna_jq,1) )then
            terluna_iq = i
            terluna_jq = j
          endif
        enddo
        enddo
        if( myid.eq.0 )then
          print *,'  Terluna: day-night forcing at each column''s own height; levels from ',terluna_az(1),' to ',  &
                  terluna_az(terluna_nl),' m'
          print *,'  Terluna: corner column, lowest level at ',zh(1,1,1),' m, forced level index ',terluna_al(1,1,1)
          print *,'  Terluna: lowest column ',terluna_ii,terluna_jj,', lowest level at ',zh(terluna_ii,terluna_jj,1),   &
                  ' m, forced level index ',terluna_al(terluna_ii,terluna_jj,1)
          print *,'  Terluna: highest column ',terluna_iq,terluna_jq,', lowest level at ',zh(terluna_iq,terluna_jq,1),  &
                  ' m, forced level index ',terluna_al(terluna_iq,terluna_jq,1)
        endif
      endif
      !$omp parallel do default(shared) private(i,j,k,terluna_l,tem)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        terluna_l = terluna_al(i,j,k)
        if( terluna_l.gt.0 )then
          tem = terluna_at(i,j,k)
          thten1(i,j,k) = thten1(i,j,k)                                                                          &
              +(1.0-terluna_f)*((1.0-tem)*terluna_ath(terluna_l,terluna_b0)+tem*terluna_ath(terluna_l+1,terluna_b0))  &
              +terluna_f*((1.0-tem)*terluna_ath(terluna_l,terluna_b1)+tem*terluna_ath(terluna_l+1,terluna_b1))
          if( imoist.eq.1 ) qten(i,j,k,nqv) = qten(i,j,k,nqv)                                                     &
              +(1.0-terluna_f)*((1.0-tem)*terluna_aqv(terluna_l,terluna_b0)+tem*terluna_aqv(terluna_l+1,terluna_b0))  &
              +terluna_f*((1.0-tem)*terluna_aqv(terluna_l,terluna_b1)+tem*terluna_aqv(terluna_l+1,terluna_b1))
        endif
      enddo
      enddo
      enddo
"""
W3_BLOCK = """    IF( var13.gt.2.5 )THEN
      ! Terluna: the large-scale vertical wind against height from terluna_wls.txt at each column's own height
      ! (boxes over terrain): interpolated to every column's w levels, zero outside the table. The corner, the
      ! lowest and the highest column are printed with the level the nudging starts from (var14 m).
      open(unit=97,file='terluna_wls.txt',status='old',action='read')
      read(97,*) terluna_n
      allocate( terluna_z(terluna_n) , terluna_w(terluna_n) )
      do terluna_l=1,terluna_n
        read(97,*) terluna_z(terluna_l),terluna_w(terluna_l)
      enddo
      close(unit=97)
      if( .not. allocated(terluna_w3) ) allocate( terluna_w3(ib:ie,jb:je,kb:ke+1) )
      terluna_w3 = 0.0
      do j=jb,je
      do i=ib,ie
      do k=2,nk
        do terluna_l=1,terluna_n-1
          if( zf(i,j,k).ge.terluna_z(terluna_l) .and. zf(i,j,k).le.terluna_z(terluna_l+1) )then
            terluna_w3(i,j,k) = terluna_w(terluna_l) + (terluna_w(terluna_l+1)-terluna_w(terluna_l))   &
                *(zf(i,j,k)-terluna_z(terluna_l))/(terluna_z(terluna_l+1)-terluna_z(terluna_l))
          endif
        enddo
      enddo
      enddo
      enddo
      wprof = 0.0
      do k=2,nk
        wprof(k) = sum( terluna_w3(1:ni,1:nj,k) )/real(ni*nj)
      enddo
      terluna_ci = 1
      terluna_cj = 1
      do j=1,nj
      do i=1,ni
        if( zf(i,j,1).lt.zf(terluna_ci(2),terluna_cj(2),1) )then
          terluna_ci(2) = i
          terluna_cj(2) = j
        endif
        if( zf(i,j,1).gt.zf(terluna_ci(3),terluna_cj(3),1) )then
          terluna_ci(3) = i
          terluna_cj(3) = j
        endif
      enddo
      enddo
      if( myid.eq.0 )then
        do terluna_q=1,3
          i = terluna_ci(terluna_q)
          j = terluna_cj(terluna_q)
          terluna_kk = nk
          do k=nk,1,-1
            if( zh(i,j,k).gt.var14 ) terluna_kk = k
          enddo
          print *,'  Terluna column ',i,j,': ground at ',zf(i,j,1),' m; nudged from level ',terluna_kk,' at ',   &
                  zh(i,j,terluna_kk),' m; large-scale w at its lowest w levels ',terluna_w3(i,j,2),terluna_w3(i,j,3)
        enddo
      endif
      deallocate( terluna_z , terluna_w )
    ENDIF

"""
TERRAIN_PATCHES = [
    ('param.F', 'Terluna: over terrain the Rayleigh layer damps toward the initial state', [
        ('      IF( do_lsnudge .and. (do_lsnudge_u .or. do_lsnudge_v .or. do_lsnudge_th) .and. irdamp.eq.1 )THEN\n',
         '      IF( do_lsnudge .and. (do_lsnudge_u .or. do_lsnudge_v .or. do_lsnudge_th) .and. irdamp.eq.1 &\n'
         '          .and. .not. terrain_flag )THEN   ! Terluna: over terrain the Rayleigh layer damps toward the initial state\n', 1),
        ('        if( terrain_flag  .or.  axisymm.eq.1  )then\n          print *\n'
         "          print *,'  Cannot use large-scale vertical velocity with terrain or axisymmetric model '\n",
         '        if( axisymm.eq.1 )then   ! Terluna: large-scale vertical wind allowed over terrain\n          print *\n'
         "          print *,'  Cannot use large-scale vertical velocity with terrain or axisymmetric model '\n", 1)]),
    ('lsnudge.F', 'Terluna: nudged column by column over terrain', [
        ("      if( terrain_flag )then\n        if(myid.eq.0)then\n        print *\n"
         "        print *,'  cannot use lsnudge with terrain (for now)  '\n",
         "      if( terrain_flag .and. .false. )then   ! Terluna: nudged column by column over terrain (solve1.F)\n"
         "        if(myid.eq.0)then\n        print *\n        print *,'  cannot use lsnudge with terrain (for now)  '\n", 1)]),
    ('solve1.F', 'Terluna: nudging over terrain', [
        ('      real :: dttmp,rtime,rdt,tem,tem0,tem1,tem2,thrad,prad,terluna_w\n',
         '      real :: dttmp,rtime,rdt,tem,tem0,tem1,tem2,thrad,prad,terluna_w\n'
         '      integer :: terluna_kz   ! Terluna: nudging over terrain\n      real :: terluna_fz\n'
         '      real, dimension(nk) :: terluna_dev\n      real, dimension(:), allocatable, save :: terluna_sw\n'
         '      integer, dimension(:,:,:), allocatable, save :: terluna_nk\n'
         '      real, dimension(:,:,:), allocatable, save :: terluna_nf,terluna_nw\n', 1),
        (NUDGE_TH, NUDGE_TH + TERRAIN_NUDGE.format(mean='thavg', q='', v='th', field='th0(i,j,k)+tha(i,j,k)',
                                                   ten='thten1(i,j,k)'), 1),
        ('          tem1 = -lsnudgefac*terluna_w*( thavg(k)-lsnudge_th(k,lsnudge_count) )/(lsnudge_tau)\n',
         '          tem1 = -lsnudgefac*terluna_w*terluna_dev(k)/(lsnudge_tau)\n', 1),
        (NUDGE_QV, NUDGE_QV + TERRAIN_NUDGE.format(mean='qavg', q=',nqv', v='qv', field='qa(i,j,k,nqv)',
                                                   ten='qten(i,j,k,nqv)'), 1),
        ('          tem1 = -lsnudgefac*terluna_w*( qavg(k,nqv)-lsnudge_qv(k,lsnudge_count) )/(lsnudge_tau)\n',
         '          tem1 = -lsnudgefac*terluna_w*terluna_dev(k)/(lsnudge_tau)\n', 1)]),
    ('init_surface.F', "Terluna: every column's ground", [
        ('      real :: sx1,sx2,sxl,stsk,stmn\n',
         '      real :: sx1,sx2,sxl,stsk,stmn\n'
         "      logical :: terluna_two   ! Terluna: every column's ground from terluna_surface2d.txt\n"
         '      integer :: terluna_nx,terluna_ny,terluna_i,terluna_j\n'
         '      real, dimension(:,:), allocatable :: terluna_tsk,terluna_tmn\n'
         '      integer, dimension(:,:), allocatable :: terluna_lu\n', 1),
        ('        close(unit=97)\n\n      ELSEIF( initsfc.ne.1 .and. initsfc.ne.2 )THEN\n',
         '        close(unit=97)\n' + SURFACE_2D_BLOCK + '\n      ELSEIF( initsfc.ne.1 .and. initsfc.ne.2 )THEN\n', 1)]),
    ('solve1.F', "Terluna: over terrain each column takes the forcing at its own height", [
        ('      integer :: terluna_b0,terluna_b1,terluna_l,terluna_m\n',
         '      integer :: terluna_b0,terluna_b1,terluna_l,terluna_m\n'
         "      integer, dimension(:,:,:), allocatable, save :: terluna_al   ! Terluna: each column's own height\n"
         '      real, dimension(:,:,:), allocatable, save :: terluna_at\n'
         '      integer :: terluna_ii,terluna_jj,terluna_iq,terluna_jq\n', 1),
        (ADVECTION_OLD, ADVECTION_COLUMNS, 1)]),
    ('input.F', "Terluna: large-scale vertical wind at each column's own height", [
        ('   ! Terluna: large-scale vertical wind by column and w level\n',
         '   ! Terluna: large-scale vertical wind by column and w level\n'
         "      real, dimension(:,:,:), allocatable :: terluna_w3   ! Terluna: large-scale vertical wind at each column's own height\n", 1)]),
    ('base.F', "Terluna: the large-scale vertical wind against height from terluna_wls.txt at each column's own height", [
        ('      real, dimension(:,:), allocatable :: terluna_wc\n',
         '      real, dimension(:,:), allocatable :: terluna_wc\n'
         '      integer :: terluna_ci(3),terluna_cj(3),terluna_q,terluna_kk   ! Terluna: sample columns\n', 1),
        ('    IF( var13.gt.1.5 )THEN\n', '    IF( var13.gt.1.5 .and. var13.lt.2.5 )THEN\n', 1),
        ('    ! boundary conditions:\n    wprof(1) = 0.0\n    wprof(nk+1) = 0.0\n',
         W3_BLOCK + '    ! boundary conditions:\n    wprof(1) = 0.0\n    wprof(nk+1) = 0.0\n', 1)]),
    ('adv_routines.F', "Terluna: each column's own wind, staggered once", [
        ('    use input, only : ib,ie,jb,je,kb,ke,ngxy,ngz,nk,rdz,terluna_w2\n',
         '    use input, only : ib,ie,jb,je,kb,ke,ngxy,ngz,nk,rdz,terluna_w2,terluna_w3\n', 1),
        ('    integer :: i,j,k\n    real :: div,wk,wk1\n',
         '    integer :: i,j,k\n    real :: div,wk,wk1\n'
         "    real, dimension(:,:,:,:), allocatable, save :: terluna_ws   ! Terluna: each column's own wind, staggered once\n", 1),
        ('    if( .not. allocated(terluna_w2) )then\n      call wsub(ix,jy,kz,a,wprof,c1,c2,mh,rr0,rf0,weps,dumz,subs)\n'
         '      return\n    endif\n',
         '    if( .not. allocated(terluna_w2) .and. .not. allocated(terluna_w3) )then\n'
         '      call wsub(ix,jy,kz,a,wprof,c1,c2,mh,rr0,rf0,weps,dumz,subs)\n      return\n    endif\n'
         '    if( .not. allocated(terluna_ws) )then   ! the winds never change: at scalar (0), u (1) and v (2) points\n'
         '      allocate( terluna_ws(ib:ie,jb:je,kb:ke+1,0:2) )\n'
         '      terluna_ws = 0.0\n'
         '      do k=kb,ke+1\n      do j=jb+1,je\n      do i=ib+1,ie\n'
         '        if( allocated(terluna_w3) )then\n'
         '          terluna_ws(i,j,k,0) = terluna_w3(i,j,k)\n'
         '          terluna_ws(i,j,k,1) = 0.5*( terluna_w3(i-1,j,k) + terluna_w3(i,j,k) )\n'
         '          terluna_ws(i,j,k,2) = 0.5*( terluna_w3(i,j-1,k) + terluna_w3(i,j,k) )\n'
         '        else\n'
         '          terluna_ws(i,j,k,0) = terluna_w2(i,k)\n'
         '          terluna_ws(i,j,k,1) = 0.5*( terluna_w2(i-1,k) + terluna_w2(i,k) )\n'
         '          terluna_ws(i,j,k,2) = terluna_w2(i,k)\n'
         '        endif\n'
         '      enddo\n      enddo\n      enddo\n'
         '    endif\n', 1),
        ('          wk = terluna_w2(i,k)\n          if( ustag.eq.1 ) wk = 0.5*( terluna_w2(i-1,k) + terluna_w2(i,k) )\n',
         '          wk = terluna_ws(i,j,k,ustag)\n', 1),
        ('        wk = terluna_w2(i,nk)\n        if( ustag.eq.1 ) wk = 0.5*( terluna_w2(i-1,nk) + terluna_w2(i,nk) )\n',
         '        wk = terluna_ws(i,j,nk,ustag)\n', 1),
        ('          wk = terluna_w2(i,k)\n          wk1 = terluna_w2(i,k+1)\n          if( ustag.eq.1 )then\n'
         '            wk = 0.5*( terluna_w2(i-1,k) + terluna_w2(i,k) )\n'
         '            wk1 = 0.5*( terluna_w2(i-1,k+1) + terluna_w2(i,k+1) )\n          endif\n',
         '          wk = terluna_ws(i,j,k,ustag)\n          wk1 = terluna_ws(i,j,k+1,ustag)\n', 1)]),
    ('adv.F', "Terluna: v points take each column's own wind", [
        ('call terluna_wsub(ni  ,nj+1,nk  ,v3d,wprof,c1,c2,mh,rr0,rf0,weps,dumz,subs,0)   ! Terluna: w by column',
         "call terluna_wsub(ni  ,nj+1,nk  ,v3d,wprof,c1,c2,mh,rr0,rf0,weps,dumz,subs,2)   ! Terluna: v points take each column's own wind", 1)]),
]

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


def design_air(run=None):
    """Carbon dioxide and oxygen mixing ratios of the design case's air and its sunlight, from its GCM run (or
    another run's)."""
    run = run or GCM_RUN
    progress = json.loads((HERE.parent / 'gcm' / 'runs' / run / 'progress.json').read_text())
    cfg = progress['configuration']
    gases = cfg['gases_bar']
    total = sum(gases.values())
    return dict(co2=gases['pCO2'] / total, o2=gases['pO2'] / total, surface_pa=cfg['pressure_pa'],
                sunlight_w_m2=progress['inputs']['flux_w_m2'], run=run)


# The curved air's light (climate/gcm/twilight.py): RRTMG lights a column only while its Sun is up and lets each
# column's infrared out through flat layers. With var20 = 1 the twilight build adds, after RRTMG, what the Moon's tall
# curved air takes beyond a flat column of the same air, by each column's solar elevation (below the horizon too) from
# terluna_twilight.txt: each shell's share of sunlight as heating of the levels it overlaps, and the ground's share to
# the net shortwave at the ground in the column's clear share; and in every column at every hour, the infrared the same
# geometry lets out beyond the flat column, each shell's at its design temperature scaled with the levels' own to the
# fourth power, and the ground's taken from the infrared reaching it. Shells above the model's top are left out with
# the air they stand for.
TWILIGHT_TABLE = 'terluna_twilight.txt'
TWILIGHT_PRODUCT = 'twilight.json'                 # the twilight ledger a case's table comes from, unless it names one
TWILIGHT_MODES = ('curved', 'full')          # as climate/gcm/exoplasim_run.py: 'full' adds the shells' widening
TWILIGHT_PATCH = ('radiation_driver.F', "Terluna: the curved air's twilight", [
    ('  public :: radiation_driver\n',
     '  public :: radiation_driver\n'
     "  ! Terluna: the curved air's twilight table, read once from terluna_twilight.txt\n"
     '  integer, save :: terluna_ntw = -1, terluna_nsh = 0\n'
     '  real, save :: terluna_lwg = 0.0, terluna_tg = 1.0\n'
     '  real, allocatable, save :: terluna_tel(:), terluna_tgr(:), terluna_tz(:), terluna_tsh(:,:), terluna_lw(:), terluna_tt(:)\n', 1),
    ('      use input, only : var12,var16,var17,var18,var19   ! Terluna: the lunar day\n',
     '      use input, only : var12,var16,var17,var18,var19,var20   ! Terluna: the lunar day, the twilight\n', 1),
    ('      real :: terluna_lat,terluna_s   ! Terluna: tilted rings\n',
     '      real :: terluna_lat,terluna_s   ! Terluna: tilted rings\n'
     '      integer :: terluna_l,terluna_m,terluna_lo,terluna_hi   ! Terluna: the twilight\n'
     '      real :: terluna_el,terluna_u,terluna_clear,terluna_e,terluna_f,terluna_ov\n'
     '      real, allocatable :: terluna_s1(:)\n', 1),
    ("                     xcoszen=coszen,yr=yr,julian=julian,mp_physics=mp_physics,co2out=co2out ) ! jararias 2013/08/14\n"
     "\n            if( nstep.le.1 .and. myid.eq.0 ) print *,'  co2 = ',co2out\n",
     "                     xcoszen=coszen,yr=yr,julian=julian,mp_physics=mp_physics,co2out=co2out ) ! jararias 2013/08/14\n"
     "\n            if( nstep.le.1 .and. myid.eq.0 ) print *,'  co2 = ',co2out\n"
     """
            IF( var20.gt.0.5 )THEN
              ! Terluna: the curved air's twilight. What the Moon's tall curved air adds to a flat column of the
              ! same air (climate/gcm/twilight.py): by each column's solar elevation, below the horizon too, each
              ! shell's share of sunlight as heating of the levels it overlaps and the ground's share in the column's
              ! clear share; at every hour, each shell's extra infrared escape as cooling, scaled with the levels'
              ! temperature (dum3) to the fourth power, and the ground's from the infrared reaching it. RRTMG's
              ! tendencies are temperature tendencies here, converted below with them.
              IF( terluna_ntw.lt.0 )THEN
                open(unit=96,file='terluna_twilight.txt',status='old',action='read')
                read(96,*) terluna_ntw,terluna_nsh
                allocate( terluna_tel(terluna_ntw), terluna_tgr(terluna_ntw), terluna_tz(terluna_nsh+1),  &
                          terluna_tsh(terluna_nsh,terluna_ntw), terluna_lw(terluna_nsh), terluna_tt(terluna_nsh) )
                read(96,*) terluna_tz
                do terluna_l=1,terluna_ntw
                  read(96,*) terluna_tel(terluna_l),terluna_tgr(terluna_l),terluna_tsh(:,terluna_l)
                enddo
                read(96,*) terluna_lw,terluna_lwg
                read(96,*) terluna_tt,terluna_tg
                close(unit=96)
                if( myid.eq.0 ) print *,'  twilight table: ',terluna_ntw,' elevations, ',terluna_nsh,' shells'
              ENDIF
              allocate( terluna_s1(terluna_nsh) )
              do j=1,nj
              do i=1,ni
                terluna_el = asin( max(-1.0,min(1.0,coszen(i,j))) )*180.0/pi
                terluna_s1(:) = 0.0
                if( terluna_el.ge.terluna_tel(1) .and. terluna_el.le.terluna_tel(terluna_ntw) )then
                  terluna_lo = 1
                  terluna_hi = terluna_ntw
                  do while( terluna_hi-terluna_lo.gt.1 )
                    terluna_m = (terluna_lo+terluna_hi)/2
                    if( terluna_tel(terluna_m).gt.terluna_el )then
                      terluna_hi = terluna_m
                    else
                      terluna_lo = terluna_m
                    endif
                  enddo
                  terluna_u = (terluna_el-terluna_tel(terluna_lo))/(terluna_tel(terluna_hi)-terluna_tel(terluna_lo))
                  terluna_s1(:) = solcon*( (1.0-terluna_u)*terluna_tsh(:,terluna_lo) + terluna_u*terluna_tsh(:,terluna_hi) )
                  terluna_clear = 1.0
                  do k=1,nk
                    terluna_clear = terluna_clear*(1.0-max(0.0,min(1.0,cldfra(i,j,k))))
                  enddo
                  gsw(i,j) = gsw(i,j) + terluna_clear*solcon*( (1.0-terluna_u)*terluna_tgr(terluna_lo)  &
                                                            + terluna_u*terluna_tgr(terluna_hi) )
                endif
                glw(i,j) = glw(i,j) - terluna_lwg*(tsk(i,j)/terluna_tg)**4
                do k=1,nk
                  terluna_e = 0.0
                  terluna_f = 0.0
                  do terluna_l=1,terluna_nsh
                    terluna_ov = min(zf(i,j,k+1),terluna_tz(terluna_l+1)) - max(zf(i,j,k),terluna_tz(terluna_l))
                    if( terluna_ov.gt.0.0 )then
                      terluna_ov = terluna_ov/(terluna_tz(terluna_l+1)-terluna_tz(terluna_l))
                      terluna_e = terluna_e + terluna_s1(terluna_l)*terluna_ov
                      terluna_f = terluna_f + terluna_lw(terluna_l)*terluna_ov*(dum3(i,j,k)/terluna_tt(terluna_l))**4
                    endif
                  enddo
                  swten(i,j,k) = swten(i,j,k) + terluna_e/(rho(i,j,k)*cp*(zf(i,j,k+1)-zf(i,j,k)))
                  lwten(i,j,k) = lwten(i,j,k) - terluna_f/(rho(i,j,k)*cp*(zf(i,j,k+1)-zf(i,j,k)))
                enddo
              enddo
              enddo
              deallocate( terluna_s1 )
            ENDIF
""", 1)])

# label: (gravity, parallel mode)
BUILDS = {'moon': (lambda: planet()['gravity_m_s2'], 'mpi'), 'earth_g': (lambda: EARTH_G, 'mpi'),
          'moon_omp': (lambda: planet()['gravity_m_s2'], 'omp'), 'earth_g_omp': (lambda: EARTH_G, 'omp'),
          'moon_omp_terrain': (lambda: planet()['gravity_m_s2'], 'omp')}
BUILD_PATCHES = {'moon_omp_terrain': TERRAIN_PATCHES}               # patches beyond PATCHES, by build
# Faster builds (the author's decision of 2026-09-29, in the decisions register). Two rebuilds of moon_omp at -O3 for
# this laptop's processor were tested on 2026-09-30 (build_check.py): moon_omp_o3 lets gfortran fuse multiplies and
# adds and call glibc's vector exp and log in the loops it vectorizes; moon_omp_o3_strict does neither (no contraction,
# and -nostdinc drops the vector-maths declarations gfortran pre-includes, the compiler's own module folder named again
# for omp_lib and the IEEE modules). Neither matched moon_omp bit for bit, both stayed within a last-bit nudge of it,
# and moon_omp_o3 ran 11-15% faster, so every OpenMP build compiles as it does (OMP_DEFAULT). moon_omp_o2 keeps CM1's
# own -O2 for comparisons; the MPI builds keep -O2, untested.
OMP_OPTS = '\nOPTS = -ffree-form -ffree-line-length-none -O2 -finline-functions --param=max-vartrack-size=0 -fopenmp\n'
OMP_DEFAULT = '-O3 -march=native'
OPTIMIZATION = {'moon_omp_o2': '-O2', 'moon_omp_o3': '-O3 -march=native',
                'moon_omp_o3_strict': '-O3 -march=native -ffp-contract=off -nostdinc'}
for _label in OPTIMIZATION:
    BUILDS[_label] = BUILDS['moon_omp']


# The rain test (the author's go-ahead of 2026-09-30): lunar gravity everywhere except in how fast rain, snow, graupel,
# ice and cloud droplets fall, which take Earth's speeds, so a box run with it shows what the slow fall at lunar
# gravity does to the rain reaching the ground and to the humidity of the lower air. Droplet activation keeps the host
# gravity.
EARTH_FALL_PATCH = ('morrison.F', 'Terluna: Earth fall speeds for the rain test', [
    ("     REAL, PRIVATE ::      VCAPR,VCAPS,VCAPG,VCAPI ! Terluna: fall-speed caps\n",
     "     REAL, PRIVATE ::      VCAPR,VCAPS,VCAPG,VCAPI ! Terluna: fall-speed caps\n"
     "     REAL, PARAMETER, PRIVATE :: GFALL = 9.81 ! Terluna: Earth fall speeds for the rain test\n", 1),
    ('(G/9.81)', '(GFALL/9.81)', 8),
    ('ACN(K) = G*RHOW/(18.*MU(K))', 'ACN(K) = GFALL*RHOW/(18.*MU(K))', 1)])
BUILDS['moon_omp_earth_fall'] = BUILDS['moon_omp']
BUILD_PATCHES['moon_omp_earth_fall'] = [EARTH_FALL_PATCH]
BUILDS['moon_omp_twilight'] = BUILDS['moon_omp']
BUILD_PATCHES['moon_omp_twilight'] = [TWILIGHT_PATCH]
# Boxes over terrain with sea around it (2026-10-15, Jura Island): the terrain build's two-dimensional surface sets each
# column's ground and land use but keeps the domain's land flag, so a water column would run CM1's soil with water's
# properties. This build makes its water columns (land use 16) sea, whose surface temperature CM1 holds.
SEA_COLUMNS_PATCH = ('init_surface.F', 'Terluna: the two-dimensional surface water columns are sea', [
    ('            tsk(i,j) = terluna_tsk(modulo(floor(xh(i)/dx),terluna_nx)+1,modulo(floor(yh(j)/dy),terluna_ny)+1)\n',
     '            tsk(i,j) = terluna_tsk(modulo(floor(xh(i)/dx),terluna_nx)+1,modulo(floor(yh(j)/dy),terluna_ny)+1)\n'
     '            xland(i,j) = merge(2.0,1.0,terluna_lu(modulo(floor(xh(i)/dx),terluna_nx)+1,'
     'modulo(floor(yh(j)/dy),terluna_ny)+1).eq.16)   ! Terluna: the two-dimensional surface water columns are sea\n', 1)])
BUILDS['moon_omp_terrain_sea'] = BUILDS['moon_omp_terrain']
BUILD_PATCHES['moon_omp_terrain_sea'] = TERRAIN_PATCHES + [SEA_COLUMNS_PATCH]
# The same with the curved air's twilight hook, which places each shell's share by its column's own face heights, so
# over terrain a column takes the shells above its ground.
BUILDS['moon_omp_terrain_sea_twilight'] = BUILDS['moon_omp_terrain']
BUILD_PATCHES['moon_omp_terrain_sea_twilight'] = TERRAIN_PATCHES + [SEA_COLUMNS_PATCH, TWILIGHT_PATCH]

# The rain limits (the author's go-ahead of 2026-10-19). Morrison's limits on rain's size and speed are Earth's, and in
# the lunar runs the rain's mean diameter sits at Earth's breakup equilibrium (rain_sizes.py). The patch scales them
# from the namelist and keeps CM1's own at its default; its header in morrison.F says what it changes, its default and
# its source. It applies after the core patches above, which give the build's gravity G and the fall-speed cap VCAPR.
RAIN_LIMITS_MARKER = 'Terluna rain limits'
RAIN_LIMITS_HEADER = """! --------------------------------------------------------------------------------------------------------------
! Terluna rain-limits patch: Terluna (github.com/dgoldman0/Terluna), climate/crm/cm1_run.py RAIN_LIMITS_PATCHES,
! 2026-10-19. The derivation is atmosphere/rain/breakup.py, its product atmosphere/rain/results/rain_breakup.json.
! What it changes. Morrison's scheme limits rain's exponential spectrum, of mean diameter Dm = 1/lambda, with Earth's
!   values. Self-collection is offset by breakup above Dm = D_b (namelist brkup_db, 300 um when 0) by the factor
!   2 - exp(2300/m (Dm - D_b)) (Verlinde and Cotton 1993), which balances at D_b + ln2/2300 m = 0.60 mm; Dm is capped
!   at 2.8 mm (LAMMINR); rain's mean fall speeds are capped at VCAPR (9.1 m/s in the reference air on Earth, scaled to
!   the build's gravity by Terluna's core patch).
! Default. brkup_gscale = 0 (namelist param3) leaves all of them as they are: D_b from brkup_db, the rate 2300/m, the
!   2.8 mm cap and the build's VCAPR.
! With brkup_gscale = 1 they follow gravity g, the build's, through the capillary length sqrt(sigma/(rho_w g)), which
!   sets the largest drop that holds together and, at equal air density, the sizes at which colliding drops reach a
!   given Weber number. With s = sqrt(9.81/g), the breakup's two lengths, D_b and 1/2300 m, scale by
!   s*brkup_rhor**brkup_rhoexp. brkup_rhor is the surface air's density over Earth's sea-level air at the same
!   temperature, which in CM1, whose air has Earth's gas constant, is the surface pressure over 1013.25 hPa, and
!   brkup_rhoexp says how much denser air, slowing the drops, softens their collisions (1.0 and 0.0 by default). The
!   cap on Dm scales by s, as the largest stable drop does, and VCAPR becomes 9.1*(g/9.81)**0.25 m/s, as that drop's
!   terminal speed does.
! Source. Colliding raindrops coalesce or break by their Weber number (Straub et al. 2010). With a single-drop solver's
!   terminal speeds (Loftus and Wordsworth 2021), the breakup-weighted sizes at which pairs reach equal Weber numbers
!   follow s to within 0.2% at Earth's air density for g from the Moon's to Earth's, and grow as
!   (rho_air/rho_air,Earth)**0.69 at 1.2 times Earth's density.
! --------------------------------------------------------------------------------------------------------------
"""
RAIN_LIMITS_BLOCK = """         D_b = cm1db  ! cutoff diameter for breakup
         ! Terluna rain limits (the header of this file): CM1's own unless brkup_gscale = 1
         KB_r = 2300.
         IF (brkup_gscale.EQ.1) THEN
           IF (brkup_rhor.LE.0.) STOP 'Terluna rain limits: brkup_rhor must be positive'
           terluna_scap = SQRT(9.81/G)
           terluna_scol = terluna_scap*brkup_rhor**brkup_rhoexp
           D_b = D_b*terluna_scol
           KB_r = KB_r/terluna_scol
           LAMMINR = LAMMINR/terluna_scap
           VCAPR = 9.1*(G/9.81)**0.25
         END IF
         if(dowr) write(outfile,*) '  Terluna rain limits: brkup_gscale',brkup_gscale,' D_b',D_b,' rate',KB_r,  &
                                   ' cap on 1/lambda',1./LAMMINR,' VCAPR',VCAPR
"""
RAIN_LIMITS_PATCHES = [
    ('input.F', RAIN_LIMITS_MARKER, [
        ('!  Begin namelists:\n',
         '      integer :: brkup_gscale            ! Terluna rain limits: 0 keeps CM1\'s own, 1 scales them with gravity\n'
         '      real :: brkup_rhor,brkup_rhoexp    ! Terluna rain limits: surface air density over Earth\'s, its exponent\n'
         '!  Begin namelists:\n', 1),
        ('l_inf,brkup_db,ndcnst,nt_c, &', 'l_inf,brkup_db,brkup_gscale,brkup_rhor,brkup_rhoexp,ndcnst,nt_c, &', 1)]),
    ('cm1.F', RAIN_LIMITS_MARKER, [
        ('      brkup_db = 0.0\n',
         '      brkup_db = 0.0\n'
         '      brkup_gscale = 0     ! Terluna rain limits: CM1\'s own unless the namelist sets 1\n'
         '      brkup_rhor = 1.0\n      brkup_rhoexp = 0.0\n', 1)]),
    ('param.F', RAIN_LIMITS_MARKER, [
        ('      call MPI_BCAST(brkup_db,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)\n',
         '      call MPI_BCAST(brkup_db,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)\n'
         '      call MPI_BCAST(brkup_gscale,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)   ! Terluna rain limits\n'
         '      call MPI_BCAST(brkup_rhor,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)\n'
         '      call MPI_BCAST(brkup_rhoexp,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)\n', 1),
        ("      if(dowr) write(outfile,*) 'brkup_db  =',brkup_db\n",
         "      if(dowr) write(outfile,*) 'brkup_db  =',brkup_db\n"
         "      if(dowr) write(outfile,*) 'brkup_gscale =',brkup_gscale   ! Terluna rain limits\n"
         "      if(dowr) write(outfile,*) 'brkup_rhor =',brkup_rhor\n"
         "      if(dowr) write(outfile,*) 'brkup_rhoexp =',brkup_rhoexp\n", 1)]),
    ('morrison.F', RAIN_LIMITS_MARKER, [
        ('', RAIN_LIMITS_HEADER, 1),
        ('     REAL, PRIVATE ::      D_b     ! cutoff diameter for breakup\n',
         '     REAL, PRIVATE ::      D_b     ! cutoff diameter for breakup\n'
         '     REAL, PRIVATE ::      KB_r    ! Terluna rain limits: the breakup law\'s rate (1/m), CM1\'s 2300 by default\n',
         1),
        ('SUBROUTINE GRAUPEL_INIT(cm1hail,cm1inum,cm1ndcnst,cm1db)\n',
         'SUBROUTINE GRAUPEL_INIT(cm1hail,cm1inum,cm1ndcnst,cm1db)\n'
         '      use input, only : brkup_gscale,brkup_rhor,brkup_rhoexp,dowr,outfile   ! Terluna rain limits\n', 1),
        ('      real,    intent(inout) :: cm1db\n',
         '      real,    intent(inout) :: cm1db\n      real :: terluna_scap,terluna_scol   ! Terluna rain limits\n', 1),
        ('         D_b = cm1db  ! cutoff diameter for breakup\n', RAIN_LIMITS_BLOCK, 1),
        ('dum=2.-exp(2300.*(1./lamr(k)-dum1))',
         'dum=2.-exp(KB_r*(1./lamr(k)-dum1))   ! Terluna rain limits: the rate, 2300/m by default', 2)]),
]
BUILDS['moon_omp_rain_limits'] = BUILDS['moon_omp']
BUILD_PATCHES['moon_omp_rain_limits'] = RAIN_LIMITS_PATCHES
# Terluna's own settings: a case's rain_limits is 'earth' (CM1's own limits, written out) or 'lunar' (the values of
# atmosphere/rain's breakup product, brkup_rhor from the case's own surface pressure).
RAIN_LIMITS_PRODUCT = 'atmosphere/rain/results/rain_breakup.json'
RAIN_LIMITS_SCHEMA = 'terluna.atmosphere.rain-breakup/1'
RAIN_LIMITS_KEYS = ('brkup_gscale', 'brkup_rhor', 'brkup_rhoexp')
SEA_LEVEL_PA = 101325.0                    # Earth's sea-level air, the reference of brkup_rhor


def rain_limits(gravity: float, settings: dict) -> dict:
    """The rain limits the rain-limits patch sets at a build's gravity from its namelist entries, as RAIN_LIMITS_BLOCK
    computes them, beside the fall-speed cap the core patch sets: the breakup cutoff (m), the breakup law's rate (1/m)
    and the mean diameter where it balances (m), the cap on the mean diameter (m) and the cap on the mean fall speeds in
    the reference air (m/s)."""
    cutoff = settings.get('brkup_db') or CM1_RAIN_BREAKUP_CUTOFF
    rate, cap = CM1_RAIN_BREAKUP_RATE, CM1_RAIN_MEAN_DIAMETER_CAP
    fall = CM1_RAIN_FALL_CAP * (gravity / CM1_GRAVITY)**((CM1_RAIN_FALL_EXPONENT + 1.0) / 3.0)
    if settings.get('brkup_gscale', 0) == 1:
        if not settings['brkup_rhor'] > 0:
            raise ValueError('brkup_rhor must be positive')
        capillary = math.sqrt(CM1_GRAVITY / gravity)
        collisional = capillary * settings['brkup_rhor']**settings['brkup_rhoexp']
        cutoff, rate, cap = cutoff * collisional, rate / collisional, cap * capillary
        fall = CM1_RAIN_FALL_CAP * (gravity / CM1_GRAVITY)**0.25
    return dict(breakup_cutoff_m=cutoff, breakup_rate_per_m=rate, equilibrium_m=cutoff + math.log(2.0) / rate,
                mean_diameter_cap_m=cap, fall_speed_cap_m_s=fall)


def rain_limit_settings(kind: str, surface_pa: float) -> tuple:
    """The &param3 entries of a case's rain limits, and the product they come from: 'earth' writes CM1's own values
    out, 'lunar' the breakup product's exponent with brkup_rhor the case's surface pressure over 1013.25 hPa."""
    earth = dict(brkup_db=CM1_RAIN_BREAKUP_CUTOFF, brkup_gscale=0, brkup_rhor=1.0, brkup_rhoexp=0.0)
    if kind == 'earth':
        return earth, None
    if kind != 'lunar':
        raise ValueError(f'rain_limits is earth or lunar, not {kind!r}')
    raw = (HERE.parents[1] / RAIN_LIMITS_PRODUCT).read_bytes()
    product = json.loads(raw)
    if product.get('schema') != RAIN_LIMITS_SCHEMA:
        raise ValueError(f'{RAIN_LIMITS_PRODUCT}: expected {RAIN_LIMITS_SCHEMA}, found {product.get("schema")}')
    namelist = product['lunar_limits']['namelist']
    settings = dict(earth, brkup_gscale=namelist['brkup_gscale'], brkup_rhoexp=namelist['brkup_rhoexp'],
                    brkup_rhor=round(surface_pa / SEA_LEVEL_PA, 4))
    return settings, {RAIN_LIMITS_PRODUCT: hashlib.sha256(raw).hexdigest()}


def with_rain_limits(text: str, settings: dict) -> str:
    """A namelist with the rain-limits entries set in &param3, added after brkup_db where the namelist lacks them."""
    import re
    start = text.index('&param3\n')
    end = text.index('\n /', start)
    block = text[start:end]
    for key in reversed(RAIN_LIMITS_KEYS):
        if not re.search(rf'^\s*{key}\s*=', block, re.M):
            line = re.search(r'^(\s*)brkup_db\s*=.*$', block, re.M)
            block = block[:line.end()] + f'\n{line.group(1)}{key} = 0,' + block[line.end():]
    text = text[:start] + block + text[end:]
    for key, value in settings.items():
        text = set_namelist(text, 'param3', key, value)
    return text


def apply_rain_limits(case: Path, cfg: dict, gravity: float):
    """Write a case's rain limits into its namelist template, from the surface pressure of its sounding; returns the
    record for case.json, or None for a case that does not set them."""
    if not cfg.get('rain_limits'):
        return None
    surface_pa = float((case / 'input_sounding').read_text().split()[0]) * 100.0
    settings, product = rain_limit_settings(cfg['rain_limits'], surface_pa)
    path = case / 'namelist.template'
    path.write_text(with_rain_limits(path.read_text(), settings))
    return dict(kind=cfg['rain_limits'], namelist=settings, limits=rain_limits(gravity, settings), product=product)


# The rain law (plan item 11, the author's go-ahead of 2026-10-21). Morrison's rain falls by one power law, which the
# lunar builds slow at every size by one factor; the single-drop solver slows each size by its own, and no power law
# follows it through both regimes. The scheme reads its law only through three moments of rain's exponential spectrum,
# so the patch takes them from a table of the solver's (atmosphere/rain/fall_law.py), with rain's self-collection at
# its lunar rate beside it. It applies after the rain-limits patch, whose limits it needs, and keeps CM1's own law and
# rate at its default.
RAIN_LAW_MARKER = 'Terluna rain law'
RAIN_LAW_TABLE = 'terluna_rain_law.txt'
RAIN_LAW_HEADER = """! --------------------------------------------------------------------------------------------------------------
! Terluna rain-law patch: Terluna (github.com/dgoldman0/Terluna), climate/crm/cm1_run.py RAIN_LAW_PATCHES, 2026-10-21.
! The derivation is atmosphere/rain/fall_law.py, its product atmosphere/rain/results/rain_fall_law.json.
! What it changes. Morrison's rain falls at V = AR (RHOSU/RHO)**0.54 D**BR, which the scheme reads only through three
!   moments of rain's exponential spectrum of slope lambda: the mass-weighted and number-weighted mean fall speeds UMR
!   and UNR (its sedimentation, and its collection of snow and graupel) and the ventilation integral in its
!   evaporation EPSR. With rain_fall_law = 1 (namelist param3) all three come from the table terluna_rain_law.txt in
!   the run directory: on a uniform grid in ln(lambda), the natural logs of the two speeds (m/s) and of the integral
!   of D**1.5 V**0.5 exp(-lambda D) dD (m3 s**-0.5) in the reference air, interpolated linearly in ln(lambda) and held
!   at the table's ends; the density factor multiplies the speeds and its square root the integral. Rain's
!   self-collection rate, 5.78 m3/kg/s (Seifert and Beheng 2001), becomes rain_selfcoll where that is positive.
! Default. rain_fall_law = 0 and rain_selfcoll = 0 leave the scheme as it is: the build's power law and 5.78. The cap
!   on rain's mean fall speeds (VCAPR) stays, and rain's collection of cloud ice keeps the power law.
! Source. The table integrates a single-drop solver's terminal speeds (Loftus and Wordsworth 2021) in the Moon's air
!   over Morrison's spectra; the self-collection rate carries Earth's by the ratio of geometric collision kernels.
! --------------------------------------------------------------------------------------------------------------
"""
RAIN_LAW_FUNCTION = """!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REAL FUNCTION TERLUNA_RAIN(LAM,M)
      ! Terluna rain law: rain's mass-weighted (M = 1) or number-weighted (M = 2) mean fall speed (m/s), or the
      ! ventilation integral of its evaporation (M = 3, m3 s**-0.5), in the reference air, for an exponential
      ! spectrum of slope LAM (1/m): the table of terluna_rain_law.txt, linear in ln(LAM), held at its ends
      IMPLICIT NONE
      REAL, INTENT(IN) :: LAM
      INTEGER, INTENT(IN) :: M
      REAL :: X,W
      INTEGER :: I
      X = (LOG(LAM)-TERLUNA_RX0)/TERLUNA_RDX
      I = MIN(MAX(INT(X),0),TERLUNA_RN-2)
      W = MIN(MAX(X-REAL(I),0.),1.)
      TERLUNA_RAIN = EXP((1.-W)*TERLUNA_RT(M,I+1)+W*TERLUNA_RT(M,I+2))
      END FUNCTION TERLUNA_RAIN

"""
RAIN_LAW_BLOCK = """         ! Terluna rain law (the header of this file): CM1's own unless rain_fall_law = 1 or rain_selfcoll > 0
         KSC_r = 5.78
         IF (rain_selfcoll.GT.0.) KSC_r = rain_selfcoll
         IF (rain_fall_law.EQ.1 .AND. .NOT. TERLUNA_RLAW) THEN
           open(unit=95,file='terluna_rain_law.txt',status='old',action='read')
           read(95,*) TERLUNA_RN,TERLUNA_RX0,TERLUNA_RDX
           allocate( TERLUNA_RT(3,TERLUNA_RN) )
           do i=1,TERLUNA_RN
             read(95,*) TERLUNA_RT(1,i),TERLUNA_RT(2,i),TERLUNA_RT(3,i)
           enddo
           close(unit=95)
           TERLUNA_RLAW = .TRUE.
         ELSE IF (rain_fall_law.NE.0 .AND. rain_fall_law.NE.1) THEN
           STOP 'Terluna rain law: rain_fall_law must be 0 or 1'
         END IF
         if(dowr) write(outfile,*) '  Terluna rain law: rain_fall_law',rain_fall_law,' table rows',TERLUNA_RN,  &
                                   ' self-collection',KSC_r
         IF (TERLUNA_RLAW .and. dowr) write(outfile,*) '  Terluna rain law at lambda 1000/m: UMR',TERLUNA_RAIN(1000.,1), &
                                   ' UNR',TERLUNA_RAIN(1000.,2),' ventilation',TERLUNA_RAIN(1000.,3)
"""
RAIN_LAW_EPSR = ("        EPSR = 2.*PI*N0RR(K)*RHO(K)*DV(K)*                           &\n"
                 "                   (F1R/(LAMR(K)*LAMR(K))+                       &\n"
                 "                    F2R*(ARN(K)*RHO(K)/MU(K))**0.5*                      &\n"
                 "                    SC(K)**(1./3.)*CONS9/                   &\n"
                 "                (LAMR(K)**CONS34))\n")
RAIN_LAW_SEDIMENT = '      UNR = ARN(K)*CONS6/DLAMR**BR\n      UMR = ARN(K)*CONS4/(DLAMR**BR)\n'
RAIN_LAW_PATCHES = [
    ('input.F', RAIN_LAW_MARKER, [
        ("      real :: brkup_rhor,brkup_rhoexp    ! Terluna rain limits: surface air density over Earth's, its exponent\n",
         "      real :: brkup_rhor,brkup_rhoexp    ! Terluna rain limits: surface air density over Earth's, its exponent\n"
         '      integer :: rain_fall_law           ! Terluna rain law: 0 keeps the power law, 1 reads terluna_rain_law.txt\n'
         "      real :: rain_selfcoll              ! Terluna rain law: rain's self-collection rate (m3/kg/s), 5.78 when 0\n",
         1),
        ('brkup_rhor,brkup_rhoexp,ndcnst,nt_c, &', 'brkup_rhor,brkup_rhoexp,rain_fall_law,rain_selfcoll,ndcnst,nt_c, &', 1)]),
    ('cm1.F', RAIN_LAW_MARKER, [
        ('      brkup_rhoexp = 0.0\n',
         '      brkup_rhoexp = 0.0\n'
         '      rain_fall_law = 0    ! Terluna rain law: the power law unless the namelist sets 1\n'
         '      rain_selfcoll = 0.0\n', 1)]),
    ('param.F', RAIN_LAW_MARKER, [
        ('      call MPI_BCAST(brkup_rhoexp,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)\n',
         '      call MPI_BCAST(brkup_rhoexp,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)\n'
         '      call MPI_BCAST(rain_fall_law,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)   ! Terluna rain law\n'
         '      call MPI_BCAST(rain_selfcoll,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)\n', 1),
        ("      if(dowr) write(outfile,*) 'brkup_rhoexp =',brkup_rhoexp\n",
         "      if(dowr) write(outfile,*) 'brkup_rhoexp =',brkup_rhoexp\n"
         "      if(dowr) write(outfile,*) 'rain_fall_law =',rain_fall_law   ! Terluna rain law\n"
         "      if(dowr) write(outfile,*) 'rain_selfcoll =',rain_selfcoll\n", 1)]),
    ('morrison.F', RAIN_LAW_MARKER, [
        ('', RAIN_LAW_HEADER, 1),
        ("     REAL, PRIVATE ::      KB_r    ! Terluna rain limits: the breakup law's rate (1/m), CM1's 2300 by default\n",
         "     REAL, PRIVATE ::      KB_r    ! Terluna rain limits: the breakup law's rate (1/m), CM1's 2300 by default\n"
         "     REAL, PRIVATE ::      KSC_r   ! Terluna rain law: rain's self-collection rate (m3/kg/s), CM1's 5.78 by default\n"
         "     LOGICAL, PRIVATE ::   TERLUNA_RLAW = .FALSE.   ! Terluna rain law: rain's moments from terluna_rain_law.txt\n"
         '     INTEGER, PRIVATE ::   TERLUNA_RN = 0\n'
         '     REAL, PRIVATE ::      TERLUNA_RX0,TERLUNA_RDX\n'
         '     REAL, ALLOCATABLE, PRIVATE :: TERLUNA_RT(:,:)\n', 1),
        ('SUBROUTINE GRAUPEL_INIT(cm1hail,cm1inum,cm1ndcnst,cm1db)\n',
         RAIN_LAW_FUNCTION + 'SUBROUTINE GRAUPEL_INIT(cm1hail,cm1inum,cm1ndcnst,cm1db)\n'
         '      use input, only : rain_fall_law,rain_selfcoll   ! Terluna rain law\n', 1),
        ("                                   ' cap on 1/lambda',1./LAMMINR,' VCAPR',VCAPR\n",
         "                                   ' cap on 1/lambda',1./LAMMINR,' VCAPR',VCAPR\n" + RAIN_LAW_BLOCK, 1),
        ('            UMR = ARN(K)*CONS4/(LAMR(K)**BR)\n',
         '            IF (TERLUNA_RLAW) THEN   ! Terluna rain law: the table\'s moments\n'
         '            UMR = ARN(K)/AR*TERLUNA_RAIN(LAMR(K),1)\n            ELSE\n'
         '            UMR = ARN(K)*CONS4/(LAMR(K)**BR)\n            END IF\n', 4),
        ('            UNR = ARN(K)*CONS6/LAMR(K)**BR\n',
         '            IF (TERLUNA_RLAW) THEN   ! Terluna rain law: the table\'s moments\n'
         '            UNR = ARN(K)/AR*TERLUNA_RAIN(LAMR(K),2)\n            ELSE\n'
         '            UNR = ARN(K)*CONS6/LAMR(K)**BR\n            END IF\n', 4),
        (RAIN_LAW_SEDIMENT,
         '      IF (TERLUNA_RLAW) THEN   ! Terluna rain law: the table\'s moments\n'
         '      UNR = ARN(K)/AR*TERLUNA_RAIN(DLAMR,2)\n      UMR = ARN(K)/AR*TERLUNA_RAIN(DLAMR,1)\n      ELSE\n'
         + RAIN_LAW_SEDIMENT + '      END IF\n', 1),
        (RAIN_LAW_EPSR,
         "        IF (TERLUNA_RLAW) THEN   ! Terluna rain law: the table's ventilation integral\n"
         "        EPSR = 2.*PI*N0RR(K)*RHO(K)*DV(K)*                           &\n"
         "                   (F1R/(LAMR(K)*LAMR(K))+                       &\n"
         "                    F2R*(ARN(K)/AR*RHO(K)/MU(K))**0.5*                      &\n"
         "                    SC(K)**(1./3.)*TERLUNA_RAIN(LAMR(K),3))\n"
         "        ELSE\n" + RAIN_LAW_EPSR + "        END IF\n", 2),
        ('nragg(k) = -5.78*dum*qr3d(k)*nr3d(k)*rho(k)',
         'nragg(k) = -KSC_r*dum*qr3d(k)*nr3d(k)*rho(k)   ! Terluna rain law: the self-collection rate, 5.78 by default', 2)]),
]
BUILDS['moon_omp_rain_law'] = BUILDS['moon_omp']
BUILD_PATCHES['moon_omp_rain_law'] = RAIN_LIMITS_PATCHES + RAIN_LAW_PATCHES + [TWILIGHT_PATCH]
# CM1's simple boundary layer (ipbl = 2, Bryan and Rotunno 2009), found 2026-10-21 in the Near-side Sea box's water
# budget and fixed with the author's go-ahead of that day. Its driver (turb.F, block dobr09) takes its tendencies from
# turbsz, turbuz and turbvz, which return the tendency of the density-weighted field (rho times the change per
# second); the ordinary turbulence path divides the density out (turbtend.F: times rr, or 1/rho at u and v points),
# but this one hands them to solve1.F, which adds thpten, qvpten, qcpten, qipten, upten and vpten as tendencies of
# theta, the mixing ratios and the wind. Every run with ipbl = 2 therefore gave its air rho times the heat, water and
# drag its surface passed: 1.32-1.35 times the evaporation and 1.46-1.51 times the sensible heat in the boxes, where the
# air at the ground weighs 1.38 kg/m3. The same code stands in the cm1r22.0 release and NCAR's main branch. The patch
# divides the density out where the driver's block ends; solve1.F's averaging of upten and vpten from scalar points,
# which this scheme fills at u and v points, is left as CM1 has it.
PBL_DENSITY_MARKER = 'Terluna: the simple PBL tendencies per unit mass'
PBL_DENSITY_PATCH = ('turb.F', PBL_DENSITY_MARKER, [
    ('      ENDIF\n\n      endif  dobr09\n', """      ENDIF

      ! Terluna: the simple PBL tendencies per unit mass. turbsz, turbuz and turbvz return rho times the change per
      ! second, as turbtend.F's own path expects before it divides by the density; solve1.F adds these as tendencies
      ! of theta, the mixing ratios and the wind, so the density is divided out here (at u and v points for the wind).
      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=1,nk
        do j=1,nj
        do i=1,ni
          thpten(i,j,k) = thpten(i,j,k)*rr(i,j,k)
          qvpten(i,j,k) = qvpten(i,j,k)*rr(i,j,k)
          qcpten(i,j,k) = qcpten(i,j,k)*rr(i,j,k)
          qipten(i,j,k) = qipten(i,j,k)*rr(i,j,k)
        enddo
        enddo
        do j=1,nj
        do i=1,ni+1
          upten(i,j,k) = upten(i,j,k)/(0.5*(rho(i-1,j,k)+rho(i,j,k)))
        enddo
        enddo
        do j=1,nj+1
        do i=1,ni
          vpten(i,j,k) = vpten(i,j,k)/(0.5*(rho(i,j-1,k)+rho(i,j,k)))
        enddo
        enddo
      enddo

      endif  dobr09
""", 1)])
BUILDS['moon_omp_rain_law_pbl'] = BUILDS['moon_omp']
BUILD_PATCHES['moon_omp_rain_law_pbl'] = RAIN_LIMITS_PATCHES + RAIN_LAW_PATCHES + [TWILIGHT_PATCH, PBL_DENSITY_PATCH]
RAIN_LAW_PRODUCT = 'atmosphere/rain/results/rain_fall_law.json'
RAIN_LAW_SCHEMA = 'terluna.atmosphere.rain-fall-law/1'
RAIN_LAW_KEYS = ('rain_fall_law', 'rain_selfcoll')


def rain_law_settings(kind: str) -> tuple:
    """The &param3 entries of a case's rain law, the text of the table it reads and the product they come from:
    'earth' writes CM1's own out (the power law and 5.78, no table), 'lunar' the fall-law product's."""
    earth = dict(rain_fall_law=0, rain_selfcoll=0.0)
    if kind == 'earth':
        return earth, None, None
    if kind != 'lunar':
        raise ValueError(f'rain_law is earth or lunar, not {kind!r}')
    raw = (HERE.parents[1] / RAIN_LAW_PRODUCT).read_bytes()
    product = json.loads(raw)
    if product.get('schema') != RAIN_LAW_SCHEMA:
        raise ValueError(f'{RAIN_LAW_PRODUCT}: expected {RAIN_LAW_SCHEMA}, found {product.get("schema")}')
    t = product['table']
    rows = zip(t['ln_mass_weighted'], t['ln_number_weighted'], t['ln_ventilation'])
    text = (f"{t['points']} {t['ln_lambda_first']:.9e} {t['ln_lambda_step']:.9e}\n"
            + ''.join(f'{a:.9e} {b:.9e} {c:.9e}\n' for a, b, c in rows))
    settings = dict(rain_fall_law=int(product['namelist']['rain_fall_law']),
                    rain_selfcoll=float(product['namelist']['rain_selfcoll']))
    return settings, text, {RAIN_LAW_PRODUCT: hashlib.sha256(raw).hexdigest()}


def with_rain_law(text: str, settings: dict) -> str:
    """A namelist with the rain-law entries set in &param3, added after brkup_rhoexp where the namelist lacks them."""
    import re
    start = text.index('&param3\n')
    end = text.index('\n /', start)
    block = text[start:end]
    for key in reversed(RAIN_LAW_KEYS):
        if not re.search(rf'^\s*{key}\s*=', block, re.M):
            line = re.search(r'^(\s*)brkup_rhoexp\s*=.*$', block, re.M)
            block = block[:line.end()] + f'\n{line.group(1)}{key} = 0,' + block[line.end():]
    text = text[:start] + block + text[end:]
    for key, value in settings.items():
        text = set_namelist(text, 'param3', key, value)
    return text


def apply_rain_law(case: Path, cfg: dict):
    """Write a case's rain law into its namelist template and its table beside it; returns the record for case.json,
    or None for a case that does not set it. The lunar law needs the lunar rain limits, whose cap and breakup it
    assumes."""
    if not cfg.get('rain_law'):
        return None
    if cfg['rain_law'] == 'lunar' and cfg.get('rain_limits') != 'lunar':
        raise ValueError('the lunar rain law needs the lunar rain limits')
    settings, table, product = rain_law_settings(cfg['rain_law'])
    path = case / 'namelist.template'
    path.write_text(with_rain_law(path.read_text(), settings))
    if table:
        (case / RAIN_LAW_TABLE).write_text(table)
    return dict(kind=cfg['rain_law'], namelist=settings, table=RAIN_LAW_TABLE if table else None, product=product)

# Electrified storms (stage 2 of the atmospheric-electricity study): WRF-ELEC's NSSL microphysics with the field,
# lightning and leakage of climate/crm/cm1_elec.py, at the Moon's gravity and at Earth's for the benchmark storm.
BUILD_FILES = {}                                                     # whole source files put in before patching, by build
for _label, _host in (('moon_omp_elec', 'moon_omp'), ('earth_g_omp_elec', 'earth_g_omp')):
    BUILDS[_label] = BUILDS[_host]
    BUILD_PATCHES[_label] = cm1_elec.PATCHES
    BUILD_FILES[_label] = cm1_elec.sources(CM1_HOME)


def optimization_for(label: str, mode: str):
    """The optimisation flags a build compiles with: its own if it has any, OMP_DEFAULT for the other OpenMP builds, and
    the Makefile's for the MPI builds (None)."""
    return OPTIMIZATION.get(label, OMP_DEFAULT if mode == 'omp' else None)


def optimization_patch(label: str, flags: str):
    """The Makefile patch that compiles an OpenMP build with flags in place of -O2."""
    return ('Makefile', f'Terluna: optimisation for {label}',
            [(OMP_OPTS, OMP_OPTS.replace(' -O2 ', f' {flags} ').rstrip('\n') + f'   # Terluna: optimisation for {label}\n', 1)])


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
    for name, text in BUILD_FILES.get(label, {}).items():
        (src / name).write_text(text(), encoding='latin-1')
        hashes[name] = hashlib.sha256((src / name).read_bytes()).hexdigest()[:16]
    patches = [MAKEFILE[mode], *PATCHES, *BUILD_PATCHES.get(label, [])]
    optimization = optimization_for(label, mode)
    if optimization and '-nostdinc' in optimization:     # omp_lib and the IEEE modules live in the compiler's own folder
        finclude = subprocess.run(['gfortran', '-print-file-name=finclude'], capture_output=True, text=True,
                                  check=True).stdout.strip()
        optimization += f' -I{finclude} -fintrinsic-modules-path {finclude}'
    if optimization:
        patches.append(optimization_patch(label, optimization))
    for name, marker, edits in patches:
        path = src / name
        path.write_text(apply_patch(path.read_text(encoding='latin-1'), name, marker, edits), encoding='latin-1')
    for name in sorted({p[0] for p in patches}):
        hashes[name] = hashlib.sha256((src / name).read_bytes()).hexdigest()[:16]
    log = folder / 'build.log'
    with open(log, 'w') as out:
        done = subprocess.run(['make', f'-j{jobs}', f'TERLUNA={flags}'], cwd=src, stdout=out, stderr=subprocess.STDOUT)
    if done.returncode != 0 or not (folder / 'run' / 'cm1.exe').exists():
        raise RuntimeError(f'CM1 build {label} failed; see {log}')
    (folder / 'run' / 'cm1.exe').replace(exe)
    record = dict(label=label, gravity_m_s2=gravity, parallel=mode, flags=flags, optimization=optimization or '-O2',
                  source=SOURCE, patched=hashes, executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest()[:16], air=air)
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
# The 45-degree rings again, on the design case GCM_RUN names and its land's wetness (their first runs took the run
# before it, whose land evaporated ten times as much), with the curved air's twilight and side light added to RRTMG.
for _tag in ('45n', '45s'):
    CASES[f'ring_{_tag}_lsw_twilight'] = dict(
        CASES[f'ring_{_tag}_lsw'], build='moon_omp_twilight', twilight='full',
        purpose=CASES[f'ring_{_tag}_lsw']['purpose'] + ', with the curved air\'s twilight and side light '
                '(climate/results/gcm/twilight.json)')
# The 45-degree rings on the design case of the GCM's corrected radiation (the curved air's light and infrared, CO2
# and O2 in its sunlight; 2026-10-14), at the sunlight scale of 0.970 that holds the design temperature, with and
# without the curved air's terms in CM1: a case's 'gcm' run gives its reference, vertical wind, land and sunlight.
CORRECTED_GCM = ('A28_corrected_970', (10, 19))
for _tag in ('45n', '45s'):
    CASES[f'ring_{_tag}_corrected'] = dict(
        CASES[f'ring_{_tag}_lsw'], gcm=CORRECTED_GCM,
        purpose=CASES[f'ring_{_tag}_lsw']['purpose'].replace('the 5% dimmer shield', 'the 3% dimmer shield')
        .replace("the GCM design case's", "the corrected radiation's design case's"))
    CASES[f'ring_{_tag}_corrected_twilight'] = dict(
        CASES[f'ring_{_tag}_corrected'], build='moon_omp_twilight', twilight='full',
        purpose=CASES[f'ring_{_tag}_corrected']['purpose'] + ', with the curved air\'s twilight and side light '
                '(climate/results/gcm/twilight.json)')


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
# ring_70_45e again with its ground and raised lakes at their real heights, in the build that nudges the upper air
# column by column at each column's own height. Only the terrain and what it requires differ from the flat
# ring_70_45e, its twin: the ground and lake temperatures carried to their heights, and the Rayleigh layer at the top
# (above 115 km) damping toward the initial state, as CM1 requires with terrain, where the flat ring damps toward the
# domain mean.
CASES['ring_70_45e_terrain'] = dict(CASES['ring_70_45e'], build='moon_omp_terrain', terrain=dict(passes=2),
                                    purpose=CASES['ring_70_45e']['purpose'] + ', with the ground and lakes at their real heights')
# Land-use rows given over to the placeholder land at these wetnesses. Rows 31-33 (the USGS table's residential and
# commercial land, which no case uses) hold the driest since 2026-09-30: the corrected GCM's ground poleward of 45
# degrees is drier than 0.05 almost everywhere (median 0.005-0.013), and cases set up before then, the three flat rings
# among them, gave such land 0.05.
WETNESS_CLASSES = ((31, 0.003), (32, 0.01), (33, 0.02), (20, 0.05), (21, 0.10), (22, 0.15), (23, 0.20), (25, 0.30),
                   (26, 0.40), (27, 0.50), (28, 0.60), (29, 0.75), (30, 0.90))
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
# The land test on dry ground: those boxes' land evaporated 101-105 W/m2, so their land surface mattered little to the
# air. The 45 N ring on the corrected radiation's design case is a desert whose land gives the air its heat, and by
# day its air near the ground runs far warmer than the GCM's. This ring is ring_45n_corrected with the GCM's roughness
# and soil, keeping the GCM's own wetness there as the control does, so the pair differs in those two alone.
CASES['ring_45n_corrected_gcm_land'] = dict(
    CASES['ring_45n_corrected'], land={k: GCM_LAND[k] for k in ('roughness_cm', 'therin')},
    purpose=CASES['ring_45n_corrected']['purpose'] + ', its land given the GCM\'s roughness (2 m) and soil thermal '
            'inertia')
# The confirming test of the heat budget (the author's go-ahead of 2026-10-16): ring_45n_corrected with the heat and
# vapour the GCM's circulation carries sideways at 45 N imposed by height, the same at every hour through CM1's
# large-scale advection, so the pair differs in that alone.
CASES['ring_45n_corrected_sideways'] = dict(
    CASES['ring_45n_corrected'], sideways=True,
    purpose=CASES['ring_45n_corrected']['purpose'] + ', with the heat and vapour the GCM\'s circulation carries '
            'sideways at that latitude imposed by height')
# The rain test's box: box_0e_small's inputs as they were written, run with Earth's fall speeds, so the pair differs in
# how fast rain, snow, graupel, ice and cloud droplets fall and in nothing else. Both keep the forcing of the GCM before
# its correction, with which box_0e_small was set up.
CASES['box_0e_small_earth_fall'] = dict(
    CASES['box_0e_small'], build='moon_omp_earth_fall', inputs_from='box_0e_small',
    purpose='box_0e_small\'s inputs as written, with rain, snow, graupel, ice and cloud droplets falling at Earth\'s speeds')
# The rain-limits test (proposed on 2026-10-19, after the author's go-ahead to patch and test the limits that day; not
# set up or run): box_0e_small's inputs as they were written, on the rain-limits build, once with CM1's own rain limits
# written out and once with the lunar ones, so the pair shares one executable and differs in its namelist's rain limits
# alone. Both write the budgets of potential temperature and vapour and the radiative heating, for the column's heating
# by process.
RAIN_TEST_OUTPUT = dict(param9=dict(output_thbudget=1, output_qvbudget=1, output_radten=1))
for _tag, _kind in (('rain_control', 'earth'), ('lunar_rain', 'lunar')):
    CASES[f'box_0e_small_{_tag}'] = dict(
        CASES['box_0e_small'], build='moon_omp_rain_limits', inputs_from='box_0e_small', rain_limits=_kind,
        namelist=RAIN_TEST_OUTPUT,
        purpose='box_0e_small\'s inputs as written, on the rain-limits build with ' + (
            'CM1\'s own rain limits written out' if _kind == 'earth' else
            'the lunar rain limits of atmosphere/rain/results/rain_breakup.json') + ', writing the heat and vapour budgets')
# The rain-limits test in its cheaper form, which the author chose to run on the night of 19 October: both arms run the
# small box's second lunar day again from its day-29.5 restart, on the one rain-limits executable, differing in their
# rain limits alone. The drops' sizes settle within an hour of a change, so the lunar day is ample.
for _tag in ('rain_control', 'lunar_rain'):
    CASES[f'box_0e_small_{_tag}_day29'] = dict(
        CASES[f'box_0e_small_{_tag}'], inputs_from=None, restart_from=dict(case='box_0e_small', day=29.5), days=59.0,
        purpose=CASES[f'box_0e_small_{_tag}']['purpose'] + "; the second lunar day again from box_0e_small's day-29.5 "
                'restart')
# The pair of boxes on the design climate (the author's go-ahead of 2026-10-21): CM1's next runs, with the lunar rain
# limits and the rain law, where the GCM's rain evaporation is put to the test. One stands over the Near-side Sea at
# 5 S, 45 W in the southern Procellarum region, open water 120 km from the nearest land; its twin over the land at 0 E,
# where the small box stood. No ring passes either on the design climate, so each takes its forcing from the GCM's
# column at its site (gcm_site.py): the GCM's mean air there as the reference, its mean vertical wind, and its own
# budget of heat and vapour by local time, from the cold-cloud run's branch that writes PlaSim's moistening by process.
# The sea is held at the GCM's temperature there; the land is the placeholder land at the GCM's wetness, standing at
# the GCM's ground there, 0.7-1.3 km up, with heights above it. Both have the small box's size and spacing, start at
# local sunrise and run two lunar days, the second analysed, with the curved air's twilight from the design climate's
# ledger, and write the heat and vapour budgets for the share of falling rain that evaporates. The first pair took the
# forcing to 16 km and held the air to the GCM's from 8 km, fully from 16 km; over the sea the nudging stood in for the
# GCM's deep convection and the box barely rained (gcm_site.py). With the author's go-ahead of 2026-10-21 the forcing
# reaches the top, the air is held only above the top of the GCM's convection (nudge held_from 'convection_top': from
# the GCM's first level above it, fully from the next), the wind still nudged throughout, and the build divides the
# density out of CM1's simple boundary layer's tendencies (PBL_DENSITY_PATCH).
GCM_SITE_RUN, GCM_SITE_YEARS = 'A28R_301k_cold_moist', (0, 5)
SITE_BOX = dict(CASES['ring'], kind='box', nx=32, ny=32, build='moon_omp_rain_law_pbl', gcm=(GCM_SITE_RUN, GCM_SITE_YEARS),
                large_scale_w='gcm', day_night='gcm', land_moisture='gcm', land=None, surface_output=True,
                rain_limits='lunar', rain_law='lunar', twilight='full', twilight_product='twilight_301k.json',
                namelist=RAIN_TEST_OUTPUT, wind_by_local_time=True,
                nudge=dict(CASES['ring']['nudge'], ramp_s=0.0, held_from='convection_top'))
for _name, _site, _where in (('box_near_side_sea', dict(lat_deg=-5.0, lon_deg=-45.0, surface='sea'),
                              'over the Near-side Sea at 5 S, 45 W'),
                             ('box_0e_design', dict(lat_deg=0.0, lon_deg=0.0, surface='land'), 'over the land at 0 N, 0 E')):
    CASES[_name] = dict(
        SITE_BOX, site=dict(gcm=True, half_lat_deg=5.6, half_lon_deg=8.5, **_site),
        start_hour_angle_deg=-90.0 - _site['lon_deg'],
        purpose=f'a 3-D box 192 km square {_where} on the design climate (the cold-cloud run at 298.5 K), forced by the '
                "GCM's column there through its whole depth: its mean air held only above the top of the GCM's "
                'convection, its mean vertical wind and its own heating and moistening by local time; the lunar rain '
                "limits and rain law, CM1's boundary layer per unit mass, the twilight, two lunar days from sunrise")
# The sea box's third lunar day (the plan of 2026-10-24): box_near_side_sea run on from its day-59 restart through days
# 59-88.5 on the same build, inputs and settings, its wind nudged by local time as before (setup_restart carries the
# windows on by whole lunar days), to see whether its upper air settles and its rain by time of day holds.
CASES['box_near_side_sea_day3'] = dict(
    CASES['box_near_side_sea'], restart_from=dict(case='box_near_side_sea', day=59.0), days=88.5,
    purpose=CASES['box_near_side_sea']['purpose'] + '; its third lunar day, days 59-88.5, from its day-59 restart')
# The resolution test of the sea box's quiet nights (the author's go-ahead of 25 October 2026). After midnight the sea
# box barely rains where the GCM with its new convection rains 1.45 mm/day, and its night boundary layer, 6.1 km deep,
# is as deep as its columns are wide, so its eddies go unresolved; coarse grids delay convection's onset (Petch and
# others 2002). Three boxes start from the settled box's own night: box_near_side_sea_day3's air at day 84.5 (its
# output 677, 40 degrees of hour angle past midnight on its third lunar day and 14 minutes into a window of its wind by
# local time) averaged over its columns, with CM1's random perturbations of 0.25 K. Each runs two Earth days under that
# hour's forcing held fixed (setup_night): the clock that times the Sun and the day-night table (var19) stopped at that
# hour, the wind nudged toward that window's for the whole run, and the source's hold on the air above the GCM's
# convection, mean vertical wind, sea, twilight, rain law and limits, build and physics. They differ in their grid
# alone: 6 km on 32 x 32 columns as the box (192 km, the control), 3 km on 64 x 64 (192 km) and 2 km on 64 x 64
# (128 km), the first time step divided as the spacing. To write little, the 3-D output comes every 6 hours without the
# budget fields (the domain diagnostics keep their means; CM1 writes the radiation's whenever it runs radiation), the
# domain diagnostics every half hour, the statistics every 15 minutes and restarts every 12 hours.
NIGHT_SOURCE = dict(case='box_near_side_sea_day3', day=84.5)
NIGHT_OUTPUT = dict(param9=dict(output_thbudget=0, output_qvbudget=0))
STOPPED_CLOCK_S = 1.0e15            # var19 for a held hour: the Sun and the table move 6e-8 degrees in two days
for _tag, _refine, _columns in (('6km', 1, 32), ('3km', 2, 64), ('2km', 3, 64)):
    CASES[f'box_near_side_sea_night_{_tag}'] = dict(
        CASES['box_near_side_sea'], night_from=dict(NIGHT_SOURCE, refine=_refine), nx=_columns, ny=_columns, days=2.0,
        output_s=21600.0, diag_s=1800.0, stat_s=900.0, restart_s=43200.0, segment_s=86400.0, namelist=NIGHT_OUTPUT,
        purpose=f"the sea box's settled night at {6 // _refine} km spacing on {_columns} x {_columns} columns: "
                "box_near_side_sea_day3's air at day 84.5, 40 degrees past midnight, averaged over its columns, and "
                "the forcing of that hour held for two Earth days")
# The sea box's third lunar day with Earth's boundary-layer length (the author's go-ahead of 27 October 2026). CM1's
# simple boundary layer (ipbl 2) mixes with K = l^2 sqrt(max(S^2 - N^2/Pr, small)), whose asymptotic length l_inf
# case_settings stretches by the gravities from Earth's 75 m to 453 m, a scaling no experiment here has tested; at the
# top of the sea's night boundary layer that scheme carries most of the heat and vapour exchanged (box_lid.py). The case
# runs box_near_side_sea_day3 again from the same day-59 restart, sunrise to sunrise, on the same build, inputs, forcing
# and physics with l_inf at 75 m, a namelist value (restart_namelist, applied over the source's namelist). Its domain
# diagnostics keep the day-3 run's 3-hour cadence, so the two compare output for output; its 3-D output comes every 6
# hours without the budgets' fields, as the held nights' does; and it keeps only its newest two restarts besides the one
# it starts from (keep_restarts). The held night at 6 km with the same length is the quick look before it.
CASES['box_near_side_sea_day3_linf75'] = dict(
    CASES['box_near_side_sea_day3'], output_s=21600.0, keep_restarts=2,
    restart_namelist=dict(NIGHT_OUTPUT, param3=dict(l_inf=75.0)),
    purpose=CASES['box_near_side_sea_day3']['purpose'] + "; again with Earth's boundary-layer length, l_inf 75 m")
CASES['box_near_side_sea_night_6km_linf75'] = dict(
    CASES['box_near_side_sea_night_6km'], namelist=dict(NIGHT_OUTPUT, param3=dict(l_inf=75.0)),
    purpose=CASES['box_near_side_sea_night_6km']['purpose'] + "; with Earth's boundary-layer length, l_inf 75 m")
# The sea box forced by the fifteen-level GCM (register, 27 October): whether CM1's quiet nights survive the
# circulation of the GCM whose night rain is in question, which differs from the Kuo run's that forced the box. The
# case is box_near_side_sea with one change, the run whose column forces it: the fifteen-level Betts-Miller baseline
# rerun on today's shared source, its years 0-1 bit for bit the baseline's, and run on to four years. With it come that
# run's reference air, sea, mean vertical wind, heating and moistening by local time, wind by local time and the height
# the air is held from (its convection's top). Two lunar days from sunrise, then the third from the day-59 restart, as
# the Kuo-forced box ran.
GCM_L15_RUN, GCM_L15_YEARS = 'A28R_301k_cold_sbm_rh07_deep200_l15cap_cwsource', (0, 3)
CASES['box_near_side_sea_l15'] = dict(
    CASES['box_near_side_sea'], gcm=(GCM_L15_RUN, GCM_L15_YEARS),
    purpose='a 3-D box 192 km square over the Near-side Sea at 5 S, 45 W, forced by the fifteen-level Betts-Miller '
            "GCM's column there (its years 0-3) through its whole depth: its mean air held only above the top of the "
            "GCM's convection, its mean vertical wind and its own heating and moistening by local time; otherwise as "
            'box_near_side_sea, two lunar days from sunrise')
CASES['box_near_side_sea_l15_day3'] = dict(
    CASES['box_near_side_sea_l15'], restart_from=dict(case='box_near_side_sea_l15', day=59.0), days=88.5,
    purpose=CASES['box_near_side_sea_l15']['purpose'] + '; its third lunar day, days 59-88.5, from its day-59 restart')
# Its fourth lunar day: the third left the box still changing, its night's rain up from 0.50 to 1.62 mm/day after
# midnight and its middle air cooling by up to 1.4 K over the day, so it runs on from the third day's restart.
CASES['box_near_side_sea_l15_day4'] = dict(
    CASES['box_near_side_sea_l15'], restart_from=dict(case='box_near_side_sea_l15_day3', day=88.5), days=118.0,
    purpose=CASES['box_near_side_sea_l15']['purpose'] + "; its fourth lunar day, days 88.5-118, from the third day's "
            'day-88.5 restart')
# The late evening on finer grids (register, 28 October). From 135 to 180 degrees of hour angle the box forced by the
# fifteen-level GCM rains 0.1-0.2 mm/day under falling inhibition while that GCM rains 2.3-3.6, and PlaSim's scheme
# applied to the box's air there does not convect either. Three boxes test whether the lull is the 6 km grid's, as the
# night was tested on 25 October: box_near_side_sea_l15_day4's air at day 107.625 (its output 862, 142 degrees of hour
# angle), averaged over its columns, with that hour's forcing held for two Earth days (setup_night), on 6, 3 and 2 km.
EVENING_SOURCE = dict(case='box_near_side_sea_l15_day4', day=107.625)
for _tag, _refine, _columns in (('6km', 1, 32), ('3km', 2, 64), ('2km', 3, 64)):
    CASES[f'box_near_side_sea_l15_evening_{_tag}'] = dict(
        CASES['box_near_side_sea_l15'], night_from=dict(EVENING_SOURCE, refine=_refine), nx=_columns, ny=_columns,
        days=2.0, output_s=21600.0, diag_s=1800.0, stat_s=900.0, restart_s=43200.0, segment_s=86400.0,
        namelist=NIGHT_OUTPUT,
        purpose=f"the late evening of the sea box forced by the fifteen-level GCM at {6 // _refine} km spacing on "
                f"{_columns} x {_columns} columns: box_near_side_sea_l15_day4's air at day 107.625, 142 degrees of "
                "hour angle, averaged over its columns, and the forcing of that hour held for two Earth days")
# The highland box (the author's go-ahead of 2026-10-01; the plan of 2026-09-29 in the decisions register): box_0e's
# size (64 x 64 columns, 385 km square) over high ground on ring_70_45e's path at 44.7 S, 246.1 E on the far side, its x
# axis along the ring. The ground stands at the atlas's heights (each column the mean of its 0.25-degree neighbourhood,
# smoothed by two passes of a 1-2-1 filter and blended toward the edge columns' mean over 50 km at each side, so it
# joins where the box wraps round), all of it land at the site's wetness, the ground temperature of every column carried
# from the ring's sea-level ground at the GCM's lapse rate; it runs in the terrain build. The upper air, the vertical
# wind and the day-night circulation come from ring_70_45e, forced from the corrected GCM.
CASES['box_highland'] = dict(
    CASES['box_0e'], site=dict(ring='ring_70_45e', lat_deg=-44.7, lon_deg=246.1, radius_m=100.0e3),
    day_night=dict(rings=('ring_70_45e',), from_day=29.5, wavenumbers=5, bins=36, harmonics=4, top_m=16000.0),
    build='moon_omp_terrain', terrain=dict(half_width_deg=0.25, passes=2, taper_km=50.0), surface_output=True,
    purpose='box_0e\'s size over high ground on ring_70_45e\'s path (44.7 S, 246.1 E), the ground at its heights, forced '
            'by that ring and the corrected GCM')
# box_highland again with the terrain build of 2026-10-02, which holds and forces each column at its own height above
# sea level: the hold to the GCM's air from 8 km above sea level over every column, the ring's day-night heating and
# moistening and the GCM's vertical wind at each column's own height. The first run placed all three by the levels of
# its corner column at 2.9 km (README, "A box over the highlands").
CASES['box_highland_own_height'] = dict(
    CASES['box_highland'],
    purpose=CASES['box_highland']['purpose'] + ', the forcing held and applied at each column\'s own height')
# Jura Island (2026-10-15): the 2.5 km arc of Montes Jura round the flooded Sinus Iridum, 58,000 km2 at 38.4-51.6 N
# and 25.4-46.9 W in the 28% scenario's atlas, as a 3-D box about the island's centre at the 45 N ring's spacing,
# 599 by 539 km, so 70-100 km of sea lie beyond its outer coasts: the island at its heights, the sea and islets round
# it from the atlas, the corrected radiation's design case as reference, and the day-night circulation of the 45 N
# twilight ring on the same design case, with the curved air's twilight added to RRTMG, so the box is complete,
# although the 45 N pair put the twilight's local effect at hundredths of a kelvin.
CASES['box_jura'] = dict(
    CASES['box_0e'], nx=100, ny=90, site=dict(ring='ring_45n_corrected', lat_deg=45.0, lon_deg=323.875, radius_m=100.0e3),
    day_night=dict(rings=('ring_45n_corrected_twilight',), from_day=29.5, wavenumbers=5, bins=36, harmonics=4,
                   top_m=16000.0),
    build='moon_omp_terrain_sea_twilight', twilight='full',
    terrain=dict(half_width_deg=0.125, passes=2, taper_km=50.0, sea=True), gcm=CORRECTED_GCM, surface_output=True,
    purpose='Jura Island, the arc of Montes Jura round the flooded Sinus Iridum (box centre 45.0 N, 36.1 W), as a 3-D '
            'box 599 by 539 km with the island at its heights and the atlas\'s sea and islets round it, the 3% dimmer '
            'shield, the corrected radiation\'s design case as reference, the 45 N twilight ring\'s day-night '
            'circulation and the curved air\'s twilight and side light (climate/results/gcm/twilight.json)')
# The benchmark storm of the electrified build (stage 2 of the atmospheric-electricity study): CM1's own supercell
# (Weisman and Rotunno 2000: the Weisman and Klemp sounding, the quarter-circle hodograph, a warm bubble, 1 km spacing,
# 120 km square, two hours) at Earth's gravity, with WRF-ELEC's NSSL microphysics with hail, its default charging, its
# branched lightning on its 0.75-s sub-steps, and output every 5 minutes; supercell_elec_cylinders has WRF-ELEC's
# cylindrical lightning in its place.
# The equatorial box electrified (stage 2): box_0e's inputs as they were written (its site, surface, grid, starting
# air and forcing, from the GCM before its correction), with the NSSL microphysics with hail in place of Morrison,
# given lunar fall speeds, its CCN at box_0e's 100 droplets per cm3, WRF-ELEC's charging and its branched lightning on
# sub-steps scaled to the box's layers and fall speeds, a downward channel striking within 5 km of the ground
# (cm1_elec.LUNAR); its restarts let windows of its storms run again with output every few minutes, other charging
# laws, leakage or the unbounded breakdown field.
# It and the windows below that ran before leakage and the leader's crossing became the lunar defaults (2026-10-05) keep
# the settings they ran with: box_0e_elec the cap until day 18, when its namelist took lightning 4 (README), and a
# channel within 5 km of the ground; none of them leakage, until the leakage windows, or the crossing.
CASES['box_0e_elec'] = dict(
    CASES['box_0e'], build='moon_omp_elec',
    elec=dict(cm1_elec.LUNAR, lightning=3, ground_m=5000.0, leakage=0, leader_v_m=0.0),
    inputs_from='box_0e',
    purpose=CASES['box_0e']['purpose'] + ', with the NSSL microphysics and WRF-ELEC\'s charging and branched lightning')
# The electrified box again with every correction of 2026-10-04 to 2026-10-06, stage 2's main run: box_0e_elec's inputs
# and two lunar days with the lunar defaults from the start (cm1_elec.LUNAR: the breakdown field's cap lifted,
# WRF-ELEC's ground rule with the leader's crossing at 1 kV/m, leakage) and the build's rain and cloud-ice ventilation at
# lunar fall speeds (2026-10-06), three-hourly output and twelve-hourly restarts,
# from whose restarts windows of its stormy days run again with output every few minutes.
CASES['box_0e_elec_corrected'] = dict(
    CASES['box_0e_elec'], elec=dict(cm1_elec.LUNAR),
    purpose=CASES['box_0e']['purpose'] + ', with the NSSL microphysics and WRF-ELEC\'s charging and branched lightning '
            'under the lunar rules of 2026-10-05: the breakdown field\'s cap lifted, the leader\'s crossing to the '
            'ground and leakage')
# Windows of the main run's stormy days with output every ten minutes, to follow each storm's life: the first lunar
# day's busiest storms (days 11.5-13.0, 129 of its 199 flashes) and the second's (days 40.5-42.0, the days the earlier
# windows ran), from box_0e_elec_corrected's twelve-hourly restarts under its own settings. CM1 on several threads
# makes each a new realization of those days.
CASES['box_0e_elec_corrected_storms_first'] = dict(
    CASES['box_0e_elec_corrected'], inputs_from=None, restart_from=dict(case='box_0e_elec_corrected', day=11.5),
    days=13.0, output_s=600.0,
    purpose=CASES['box_0e_elec_corrected']['purpose'] + '; its first lunar day\'s busiest storms again from its day-11.5 '
            'restart with output every ten minutes')
CASES['box_0e_elec_corrected_storms_second'] = dict(
    CASES['box_0e_elec_corrected'], inputs_from=None, restart_from=dict(case='box_0e_elec_corrected', day=40.5),
    days=42.0, output_s=600.0,
    purpose=CASES['box_0e_elec_corrected']['purpose'] + '; its second lunar day\'s busiest storms again from its '
            'day-40.5 restart with output every ten minutes')
# The fine box (stage 2): box_0e_elec's site at a third of its spacing (2.0 km) over a box a third as wide (128 km, the
# same 64 by 64 columns), started from box_0e_elec's air averaged over its columns and its mean skin temperature at a
# day of its run a few hours before a stormy window (day 10.75, four hours before box_0e_elec's first flash), under
# the same Sun, forcing, land and electricity, with output every 15 minutes. It resolves the storm cells and draws the
# flashes on a grid three times finer. It ran two model days without a flash and was extended toward box_0e_elec's
# busiest lightning (days 12.7-14.3; 2026-10-04, at the author's direction); its deep convection collapsed after
# coarse day 13.25 and the author stopped it at day 3.0 of its run.
CASES['box_0e_elec_fine'] = dict(
    CASES['box_0e_elec'], inputs_from=None, fine_from=dict(case='box_0e_elec', day=10.75, refine=3), days=4.0,
    output_s=900.0, restart_s=10800.0, segment_s=43200.0,
    purpose=CASES['box_0e_elec']['purpose'] + '; a box a third as wide at a third of the spacing, started from that '
            'run\'s averaged air a few hours before a stormy window')
# box_0e_elec's first lunar day's storms again (days 10.5-14.5, from its day-10.5 restart, before its first flash at day
# 10.92), with the breakdown field's 180-kV/m cap lifted (the author's decision, 2026-10-04), beside the first run under
# the cap; and the same with point discharge from the ground, which the author left to the agent's judgement, its onset
# Standler and Winn's 3 kV/m over dense vegetation scaled by the density at the ground.
CASES['box_0e_elec_uncapped'] = dict(
    CASES['box_0e_elec'], inputs_from=None, restart_from=dict(case='box_0e_elec', day=10.5), days=14.5,
    elec=dict(cm1_elec.LUNAR, ground_m=5000.0, leakage=0, leader_v_m=0.0),
    purpose=CASES['box_0e_elec']['purpose'] + '; its first lunar day\'s storms again from its day-10.5 restart, with '
            'the breakdown field\'s 180-kV/m cap lifted')
CASES['box_0e_elec_uncapped_corona'] = dict(
    CASES['box_0e_elec_uncapped'], elec=dict(cm1_elec.LUNAR, ground_m=5000.0, leakage=0, leader_v_m=0.0, corona_v_m=3000.0),
    purpose=CASES['box_0e_elec_uncapped']['purpose'] + ', and point discharge from the ground')
# Ground strikes by WRF-ELEC's own rule (ground_m -1, its default: a downward channel reaching air warmer than -7 C in
# charge of the matching sign). The Earth benchmarks ran with ground_m 0, which lightmsz reads as only the two lowest
# levels, so supercell_elec_ground_rule runs the benchmark storm again with the rule. In the lunar storms the main
# negative charge sits at -8 to -32 C with -7 C just beneath it, as on Earth, and box_0e_elec_ground_rule runs the coarse
# run's busiest storms (days 40.5-42, the cap lifted) under the same rule: an upper bound, since it takes a channel that
# reaches -7 C, 34 km up, to bridge the rest of the way to the ground.
CASES['box_0e_elec_ground_rule'] = dict(
    CASES['box_0e_elec'], inputs_from=None, restart_from=dict(case='box_0e_elec', day=40.5), days=42.0,
    elec=dict(cm1_elec.LUNAR, ground_m=-1.0, leakage=0, leader_v_m=0.0),
    purpose=CASES['box_0e_elec']['purpose'] + '; its busiest storms again from its day-40.5 restart, with WRF-ELEC\'s '
            'own ground-strike rule')
# The charging law. Saunders and Peck's law, which every run so far uses, nearly stops for slow lunar graupel, which
# rimes below its threshold; Takahashi's grows linearly with impact speed and, in stage 1's column estimate, charged the
# storms two hundred to five hundred times faster. box_0e_elec_takahashi runs box_0e_elec_ground_rule's window under
# Takahashi's law (WRF-ELEC's isaund 1, with its size and speed factor, on its table, inside whose 0 to -30 C and cloud
# water range the lunar charging zone lies), and box_0e_elec_takahashi_first box_0e_elec_uncapped's first 1.5 days
# (from day 10.5, before the first flash, with that run's 5-km ground rule); each differs from its twin only in the law.
# Both ran on 2026-10-05: as many flashes and as much charge as under Saunders and Peck, whose charging NSSL's hail does
# most of, and an inverted dipole, since Takahashi's table charges rimed ice positively at the storms' low cloud water
# (README, "The charging law").
CASES['box_0e_elec_takahashi'] = dict(
    CASES['box_0e_elec_ground_rule'], elec=dict(cm1_elec.LUNAR, isaund=1, leakage=0, leader_v_m=0.0),
    purpose=CASES['box_0e_elec']['purpose'] + '; its busiest storms again from its day-40.5 restart, under '
            'Takahashi\'s charging law')
CASES['box_0e_elec_takahashi_first'] = dict(
    CASES['box_0e_elec_takahashi'], restart_from=dict(case='box_0e_elec', day=10.5), days=12.0,
    elec=dict(cm1_elec.LUNAR, isaund=1, ground_m=5000.0, leakage=0, leader_v_m=0.0),
    purpose=CASES['box_0e_elec']['purpose'] + '; its first lunar day\'s storms again from its day-10.5 restart, under '
            'Takahashi\'s charging law')
# Leakage. With nothing to conduct it, the charge evaporating cloud leaves on the small ions stays for weeks: 2,000-3,600
# C of each sign at 20-40 km through the lunar night, and 11.5 kC net at day 40.5 (7.0 kC on snow, 4.1 kC on the ions),
# whose negative partner fell to the ground on precipitation; WRF-ELEC hands the ions' charge to each new cloud. The
# lunar clear air would relax it in about ten minutes (stage 1's conductivity). box_0e_elec_leakage and
# box_0e_elec_leakage_first run the charging-law windows under Saunders and Peck with leakage (var9 1: net charge relaxing
# at sigma/eps0, the clear air's and the cloud's conductivity by height), each differing from its twin only in that.
# Both ran on 2026-10-05: the leftover charge goes within the hour, the field at the ground away from storms falls from
# about 20 kV/m to a median 0.2-1.6, and the lightning stays as it was (README, "Leakage").
CASES['box_0e_elec_leakage'] = dict(
    CASES['box_0e_elec_ground_rule'], elec=dict(cm1_elec.LUNAR, leakage=1, leader_v_m=0.0),
    purpose=CASES['box_0e_elec']['purpose'] + '; its busiest storms again from its day-40.5 restart, with the charge '
            'leaking through the air\'s conductivity')
CASES['box_0e_elec_leakage_first'] = dict(
    CASES['box_0e_elec_leakage'], restart_from=dict(case='box_0e_elec', day=10.5), days=12.0,
    elec=dict(cm1_elec.LUNAR, leakage=1, ground_m=5000.0, leader_v_m=0.0),
    purpose=CASES['box_0e_elec']['purpose'] + '; its first lunar day\'s storms again from its day-10.5 restart, with '
            'the charge leaking through the air\'s conductivity')
# The leader's crossing (README, "The leader's crossing"). WRF-ELEC's rule counts a downward channel at -7 C, 34 km up in
# the lunar air, as reaching the ground; a leader stops when the potential difference between its tip and the air ahead
# runs out (Lalande et al. 2002; Mazur and Ruhnke), which over that path turns on the channel's internal field: about 1
# kV/m at Earth's sea-level density for a thermalized leader (Mansell 2000; Boggs et al. 2018), 1-10 kV/m in Lalande et
# al.'s range. box_0e_elec_leader_1kv and box_0e_elec_leader_10kv run box_0e_elec_leakage's window with a ground strike
# needing the crossing at those fields (scaled by density), each differing from it only in that.
CASES['box_0e_elec_leader_1kv'] = dict(
    CASES['box_0e_elec_leakage'], elec=dict(cm1_elec.LUNAR, leakage=1, leader_v_m=1000.0),
    purpose=CASES['box_0e_elec']['purpose'] + '; its busiest storms again from its day-40.5 restart, a ground strike '
            'needing a leader with a 1-kV/m internal field to cross to the ground')
CASES['box_0e_elec_leader_10kv'] = dict(
    CASES['box_0e_elec_leader_1kv'], elec=dict(cm1_elec.LUNAR, leakage=1, leader_v_m=10000.0),
    purpose=CASES['box_0e_elec']['purpose'] + '; its busiest storms again from its day-40.5 restart, a ground strike '
            'needing a leader with a 10-kV/m internal field to cross to the ground')
# The same storm with CM1's own copy of the NSSL scheme and no electricity, to check that WRF-ELEC's copy, run through
# terluna_elec.F, makes the same storm.
CASES['supercell_nssl'] = dict(
    kind='sample', sample='supercell', build='earth_g_omp', days=7200.0 / 86400.0, output_s=300.0, restart_s=1800.0,
    segment_s=7200.0, namelist=dict(param2=dict(ptype=27)),
    purpose='CM1\'s supercell at Earth\'s gravity with CM1\'s own NSSL microphysics, the electrified benchmark\'s twin')
# The benchmarks below ran with ground_m 0 (only the two lowest levels count), before WRF-ELEC's own rule became the
# default on 2026-10-05, and keep it; supercell_elec_ground_rule runs the storm with the rule.
CASES['supercell_elec'] = dict(
    kind='sample', sample='supercell', build='earth_g_omp_elec', days=7200.0 / 86400.0, output_s=300.0,
    restart_s=1800.0, segment_s=7200.0, elec=dict(ground_m=0.0),
    purpose='CM1\'s supercell at Earth\'s gravity with WRF-ELEC\'s NSSL microphysics, charging and branched lightning, '
            'the electrified build\'s benchmark')
CASES['supercell_elec_ground_rule'] = dict(
    CASES['supercell_elec'], elec=dict(ground_m=-1.0),
    purpose='CM1\'s supercell at Earth\'s gravity with WRF-ELEC\'s NSSL microphysics, charging and branched lightning '
            'and its own ground-strike rule, the benchmark again with ground strikes possible')
CASES['supercell_elec_leader_1kv'] = dict(
    CASES['supercell_elec_ground_rule'], elec=dict(ground_m=-1.0, leader_v_m=1000.0),
    purpose='CM1\'s supercell at Earth\'s gravity with WRF-ELEC\'s NSSL microphysics, charging and branched lightning '
            'and its own ground-strike rule, a ground strike needing a leader with a 1-kV/m internal field to cross to '
            'the ground: the leader\'s crossing on Earth')
CASES['supercell_elec_cylinders'] = dict(
    CASES['supercell_elec'], elec=dict(ground_m=0.0, lightning=1),
    purpose='CM1\'s supercell at Earth\'s gravity with WRF-ELEC\'s NSSL microphysics, charging and cylindrical '
            'lightning')
# The benchmark at the fine box's 2-km spacing: supercell_elec on a grid twice as coarse over the same 120 km, and the
# storm in WRF-ELEC's own test settings (its em_quarter_ss namelist: 2 km over 84 km, 800 CCN per cm3, Saunders and
# Peck's law as isaund 11 sets it, screening layers and branched lightning), each for two hours (WRF-ELEC's test runs
# one).
CASES['supercell_elec_2km'] = dict(
    CASES['supercell_elec'], namelist=dict(param0=dict(nx=60, ny=60), param1=dict(dx=2000.0, dy=2000.0)),
    purpose='CM1\'s supercell at Earth\'s gravity with WRF-ELEC\'s NSSL microphysics, charging and branched lightning '
            'at 2-km spacing, the benchmark at the fine box\'s spacing')
CASES['supercell_elec_wrf'] = dict(
    CASES['supercell_elec'], elec=dict(ground_m=0.0, isaund=11, screen=1),
    namelist=dict(param0=dict(nx=42, ny=42), param1=dict(dx=2000.0, dy=2000.0), nssl2mom_params=dict(ccn=0.8e9)),
    purpose='CM1\'s supercell at Earth\'s gravity in WRF-ELEC\'s own test settings: 2-km spacing over 84 km, 800 CCN '
            'per cm3, Saunders and Peck\'s law as isaund 11 sets it, screening layers and branched lightning')


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
    text = lambda w: f'{w:.2f}' if abs(w - round(w, 2)) < 1e-9 else f'{w:.3f}'          # 0.05, 0.90, 0.003
    rows = {index: base.replace('28,', f'{index},', 1).replace('.50', text(w).lstrip('0') if w < 1 else text(w), 1)
            .replace("'Terluna placeholder land'", f"'Terluna placeholder land, wetness {text(w)}'") for index, w in classes}
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
                lapse_k_m=np.interp(lat_path, lats, pick('lapse')),
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


def atlas_heights(lat, lon, half_width_deg=0.5):
    """Height above sea level (m) of the surface under each point, from the 28% scenario's atlas, and whether it is
    water: on water (at least half the box within half_width_deg wet, as in path_surface) the water surface, at sea
    level on the seas and at their own levels on raised lakes; on land the mean height of the box's dry cells; nothing
    below sea level."""
    import numpy as np
    from climate.gcm import boundary
    product = boundary.lakes_product(*boundary.LAKES)
    plat, plon, level = product['lat_deg'], product['lon_deg'], product['meta']['level_m']
    water, height, surface = product['water_fraction'], product['mean_height_m'], product['water_surface_m']
    raw, wet_column = np.zeros(len(lat)), np.zeros(len(lat), bool)
    for i, (la, lo) in enumerate(zip(lat, lon)):
        r = np.abs(plat - la) <= half_width_deg
        c = np.abs((plon - lo + 180.0) % 360.0 - 180.0) * np.cos(np.radians(la)) <= half_width_deg
        w, h, s = (a[np.ix_(r, c)] for a in (water, height, surface))
        wet_column[i] = w.mean() >= 0.5
        raw[i] = (s[w >= 0.5].mean() if wet_column[i] else h[w < 0.5].mean()) - level
    return np.maximum(raw, 0.0), wet_column


def path_heights(lat, lon, half_width_deg=0.5, passes=2):
    """Height above sea level (m) of the surface under each column of a path (atlas_heights), the land then smoothed
    by `passes` of a 1-2-1 filter along the ring, water held at its level."""
    import numpy as np
    raw, wet_column = atlas_heights(lat, lon, half_width_deg)
    out = raw.copy()
    for _ in range(passes):
        out = np.where(wet_column, raw, 0.25 * np.roll(out, 1) + 0.5 * out + 0.25 * np.roll(out, -1))
    return out


def box_columns(lat0: float, lon0: float, heading_deg: float, nx: int, ny: int, dx: float) -> tuple:
    """Latitude and longitude (degrees, arrays ny by nx) of a box's columns about a site at its centre: x along
    heading_deg (clockwise from north), y 90 degrees to its left, azimuthal-equidistant about the site."""
    import numpy as np
    radius = planet()['radius_m']
    x = (np.arange(nx) + 0.5) * dx - 0.5 * nx * dx
    y = (np.arange(ny) + 0.5) * dx - 0.5 * ny * dx
    xx, yy = np.meshgrid(x, y)
    a = np.radians(heading_deg)
    east, north = xx * np.sin(a) - yy * np.cos(a), xx * np.cos(a) + yy * np.sin(a)
    delta, bearing = np.hypot(east, north) / radius, np.arctan2(east, north)
    la0, lo0 = np.radians(lat0), np.radians(lon0)
    lat = np.arcsin(np.sin(la0) * np.cos(delta) + np.cos(la0) * np.sin(delta) * np.cos(bearing))
    lon = lo0 + np.arctan2(np.sin(bearing) * np.sin(delta) * np.cos(la0), np.cos(delta) - np.sin(la0) * np.sin(lat))
    return np.degrees(lat), np.degrees(lon) % 360.0


def box_heights(raw, dx: float, passes: int = 2, taper_m: float = 50.0e3, wet=None) -> dict:
    """A box's ground heights (m, ny by nx) from the atlas's (atlas_heights), smoothed by `passes` of a 1-2-1 filter in
    x and in y and blended toward the mean of the box's edge columns over taper_m at each side (a cosine taper), so the
    ground joins where the box wraps round; the largest slope and the heights' spread with them. Columns that wet marks
    as water keep their level."""
    import numpy as np
    raw = np.asarray(raw, dtype=float)
    held = np.zeros(raw.shape, bool) if wet is None else np.asarray(wet, bool)
    h = raw.copy()
    for _ in range(passes):
        for axis in (0, 1):
            h = np.where(held, raw, 0.25 * np.roll(h, 1, axis) + 0.5 * h + 0.25 * np.roll(h, -1, axis))
    edge = float(np.concatenate([h[0], h[-1], h[1:-1, 0], h[1:-1, -1]]).mean())
    ny, nx = h.shape

    def weight(n):                                                   # 0 at the edge, 1 beyond taper_m from it
        d = np.minimum(np.arange(n) + 0.5, n - np.arange(n) - 0.5) * dx
        return np.where(d < taper_m, 0.5 * (1.0 - np.cos(np.pi * d / taper_m)), 1.0)
    w = weight(ny)[:, None] * weight(nx)[None, :]
    h = np.where(held, raw, np.maximum(edge + (h - edge) * w, 0.0))
    slope = np.hypot((np.roll(h, -1, 1) - h) / dx, (np.roll(h, -1, 0) - h) / dx)
    return dict(heights=h, edge_m=edge, max_slope=float(slope.max()), p99_slope=float(np.percentile(slope, 99)),
                min_m=float(h.min()), max_m=float(h.max()), mean_m=float(h.mean()))


def box_site(cfg: dict, zw, nx_ring: int, run=GCM_RUN, years=GCM_YEARS) -> dict:
    """What a box takes from the ring it extends: the ring's GCM reference (upper air and wind along the ring), and
    over the ring's columns within the site's radius of it, the mean land wetness and ground temperature of the land
    columns and the mean GCM vertical wind on the w levels zw. A ring along a circle of latitude gives its row's
    reference (parallel_box_site), and a box no ring passes takes everything from the GCM (gcm_box_site)."""
    import numpy as np
    site = cfg['site']
    if site.get('gcm'):
        return gcm_box_site(cfg, zw, run, years)
    ring = CASES[site['ring']]
    if ring.get('kind') != 'tilted':
        return parallel_box_site(cfg, zw, nx_ring, ring, run, years)
    track = tilted_path(nx_ring, ring['tilt_deg'], ring['node_deg'])
    ref = gcm_path_reference(track['lat'], track['lon'], ring['tilt_deg'], zw, run, years)
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
    a0, a1 = np.radians(track['lat'][(centre - 1) % nx_ring]), np.radians(track['lat'][(centre + 1) % nx_ring])
    dl = np.radians(track['lon'][(centre + 1) % nx_ring] - track['lon'][(centre - 1) % nx_ring])
    heading = np.degrees(np.arctan2(np.sin(dl) * np.cos(a1), np.cos(a0) * np.sin(a1) - np.sin(a0) * np.cos(a1) * np.cos(dl)))
    return dict(ref=ref, columns=[int(c) for c in cols], land_share=float(land.mean()), surface_product=surface_sha,
                wetness=float(ref['wetness'][cols][land].mean()), ground_k=float(ref['land_k'][cols][land].mean()),
                w_m_s=[float(v) for v in np.asarray(ref['w_columns'])[cols].mean(axis=0)],
                heading_deg=float(heading % 360.0), lapse_k_m=float(np.asarray(ref['lapse_k_m'])[cols].mean()))


def gcm_box_site(cfg: dict, zw, run, years) -> dict:
    """box_site for a box no ring passes: from the GCM's column at the site (gcm_site.py), its mean air as the
    reference, its surface, its mean vertical wind on the w levels zw and the table of its own heating and moistening
    by local time; the box's x axis runs east, and it stands at the site's own surface, its heights above it."""
    from climate.crm import gcm_site
    site = cfg['site']
    f = gcm_site.site_forcing(run, years, {k: site[k] for k in ('lat_deg', 'lon_deg', 'half_lat_deg', 'half_lon_deg',
                                                                 'surface')})
    cm1 = gcm_site.to_cm1(f, zw)
    s, prof = f['surface'], f['profile']
    ref = dict(profile={k: prof[k] for k in ('z_m', 'theta_k', 'qv_kg_kg', 'u_m_s', 'v_m_s')},
               sea_surface_k=s['ground_k'], land_ground_k=s['ground_k'], land_wetness=s['wetness'],
               air_over_sea_k=s['air_2m_k'], surface_pa_sea=s['pressure_pa'], lapse_k_m=0.0,
               rows_deg=sorted({c[0] for c in f['cells']}), cells=f['cells'], run=run, years=list(years))
    return dict(ref=ref, columns=f['cells'], land_share=0.0 if site['surface'] == 'sea' else 1.0, surface_product=None,
                wetness=s['wetness'], ground_k=s['ground_k'], sea_k=s['ground_k'], w_m_s=cm1['w_m_s'],
                heading_deg=90.0, lapse_k_m=0.0, day_night=cm1['day_night'], forcing=f)


def held_heights(nudge: dict, convection: dict, ztop: float, scale: float) -> dict:
    """A GCM-fed box's nudging with its air held to the GCM's only above the GCM's convection (gcm_site's
    convection_top): from the GCM's first level above the convection's top, fully from the next or from the model's
    top, the heights divided by the case's scale as the namelist multiplies them back."""
    if convection['held_from_m'] is None:
        raise ValueError("the GCM's convection reaches its top level; there is no air above it to hold")
    full = convection['held_fully_from_m'] if convection['held_fully_from_m'] is not None else ztop
    return dict(nudge, z_start_m=convection['held_from_m'] / scale, z_full_m=full / scale)


def gcm_site_record(f: dict) -> dict:
    """What case.json keeps of a box's GCM site: its cells and outputs, the budgets behind its table, its top and where
    the GCM forms its rain."""
    return dict({k: f[k] for k in ('run', 'years', 'site', 'cells', 'cell_heights_m', 'outputs', 'phases', 'harmonics',
                                   'window_deg', 'surface', 'budget', 'convection')},
                top_m=f['day_night']['top_m'])


def parallel_box_site(cfg: dict, zw, nx_ring: int, ring: dict, run=GCM_RUN, years=GCM_YEARS) -> dict:
    """box_site for a ring along a circle of latitude: the GCM's reference on the ring's row as the ring takes it
    (gcm_equator_profile, its sea surface, its land's wetness and ground carried to sea level, its lapse), the row's
    mean vertical wind on the w levels zw, and the ring's columns within the site's radius; the box's x axis runs east
    along the parallel."""
    import numpy as np
    site = cfg['site']
    latitude, band = ring['latitude_deg'], ring.get('band_deg', EQUATOR_BAND_DEG)
    if abs(site['lat_deg'] - latitude) > 2.0 * band:
        raise ValueError(f"the site at {site['lat_deg']} degrees lies off {site['ring']} at {latitude} degrees")
    ref = gcm_equator_profile(run, years, band=band, latitude=latitude)
    wind = gcm_vertical_wind(latitude, band, run, years)
    w = np.interp(zw, wind['z_m'], wind['w_m_s'], left=0.0, right=0.0)
    length = 2.0 * math.pi * planet()['radius_m'] * math.cos(math.radians(latitude))
    dx = length / nx_ring
    gap = ((np.arange(nx_ring) + 0.5) * dx - (site['lon_deg'] % 360.0) / 360.0 * length + 0.5 * length) % length - 0.5 * length
    cols = np.flatnonzero(np.abs(gap) <= site['radius_m'])
    return dict(ref=dict(ref, sea_k=ref['sea_surface_k']), columns=[int(c) for c in cols], land_share=None,
                surface_product=None, wetness=float(ref['land_wetness']), ground_k=float(ref['land_ground_k']),
                w_m_s=[float(v) for v in w], heading_deg=90.0, lapse_k_m=float(ref['lapse_k_m']),
                sea_k=float(ref['sea_surface_k']))


def write_day_night(path: Path, forcing: dict) -> None:
    """terluna_lsadv.txt: the local-time bin and level counts, the level heights (m), then per bin (evenly spaced in
    hour angle from -180 degrees) a line of theta tendencies (K/s) and a line of vapour tendencies (kg/kg/s)."""
    import numpy as np
    z, th, qv = (np.asarray(forcing[k], float) for k in ('z_m', 'theta_k_s', 'qv_kg_kg_s'))
    path.write_text(f'{th.shape[0]} {z.size}\n' + ' '.join(f'{v:.3f}' for v in z) + '\n'
                    + ''.join(' '.join(f'{v:.6e}' for v in a) + '\n' + ' '.join(f'{v:.6e}' for v in b) + '\n'
                              for a, b in zip(th, qv)))


def sideways_forcing(run, years, latitude, band, ref) -> dict:
    """The heat and vapour the GCM's circulation carries sideways at a ring's latitude, as terluna_lsadv.txt takes
    them in one local-time bin, so the same at every hour: the zonal-mean meridional circulation's advection and the
    resolved eddies' flux convergence (ring_heat_budget.gcm_sideways), on the GCM's levels at its sea columns' heights
    with a level at the ground holding the lowest level's values; theta tendencies (K/s), the temperature tendencies
    over each level's Exner function, and vapour tendencies (kg/kg/s)."""
    import numpy as np
    from climate.crm.ring_heat_budget import gcm_sideways
    g = gcm_sideways(run, years, latitude, band)
    sigma = np.asarray(g['sigma'], float)
    heat = (np.asarray(g['mean_meridional_k_day']) + np.asarray(g['eddies_k_day'])) / 86400.0      # K/s, top first
    vapour = (np.asarray(g['vapour']['mean_meridional_g_kg_day'])
              + np.asarray(g['vapour']['eddies_g_kg_day'])) / 8.64e7                              # kg/kg/s
    exner = (sigma * ref['surface_pa_sea'] / 1.0e5) ** (287.04 / 1005.7)
    z, th, qv = np.asarray(ref['profile']['z_m'], float)[::-1], (heat / exner)[::-1], vapour[::-1]   # ground first
    return dict(z_m=np.r_[0.0, z], theta_k_s=np.r_[th[0], th][None, :], qv_kg_kg_s=np.r_[qv[0], qv][None, :],
                rows_deg=g['rows_deg'])


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
    settings = {
        'param0': dict(nx=nx, ny=cfg['ny'] if box else 1, nz=nz, ppnode=cfg.get('ranks', 8), timeformat=3, timestats=1,
                       terrain_flag=bool(cfg.get('terrain'))),
        'param1': dict(dx=round(dx, 3), dy=round(dx, 3), dz=round(ztop / nz, 1), dtl=round(40.0 * s, 3), cfl_limit=1.0,
                       timax=round(cfg['days'] * 86400.0 * s), run_time=-999.9, tapfrq=round(cfg['output_s'] * s, 3),
                       rstfrq=round(cfg['restart_s'] * s, 3), statfrq=round(3600.0 * s, 3), prclfrq=1.0e9),
        'param2': dict(cm1setup=2, testcase=0, adapt_dt=1, irst=0, rstnum=1, ipbl=2, sgsmodel=0, tconfig=2,
                       horizturb=0, irdamp=1 if cfg.get('terrain') else 2, psolver=3, ptype=5, ihail=0, icor=int(rotating),
                       betaplane=int(tilted and rotating),
                       lspgrad=int(rotating and not tilted),
                       eqtset=2, idiss=1, wbc=1, ebc=1, sbc=1, nbc=1, bbc=3, tbc=1, isnd=7, iwnd=0,
                       itern=4 if cfg.get('terrain') else 0, iinit=0,
                       irandp=1, iorigin=1, apmasscon=1),
        # The Rayleigh layer's time and the PBL scheme's asymptotic length scale are Earth's values stretched
        # by the ratio of gravities, as the dynamics stretch.
        'param3': dict(fcor=float(f'{turning:.6e}'), rdalpha=round(gravity / EARTH_G / 300.0, 8), zd=0.77 * ztop,
                       l_inf=round(75.0 * EARTH_G / gravity, 1), ndcnst=cfg['droplets_cm3']),
        'param6': dict(stretch_z=4, ztop=ztop),
        'param8': dict(**(dict(var11=1.0 if cfg.get('day_night') else 0.0) if box
                          else dict(var11=1.0) if cfg.get('sideways') else {}),
                       var12=cfg['tilt_deg'] if tilted else 0.0,
                       var13=(3.0 if box and cfg.get('terrain') else 2.0 if tilted else 1.0)
                       if cfg.get('large_scale_w') else 0.0,          # 3: at each column's own height
                       var14=cfg['nudge']['z_start_m'] * s, var15=cfg['nudge']['z_full_m'] * s,
                       var16=0.0 if pair or box else round(length, 3), var17=0.0 if pair else round(air['sunlight_w_m2'], 3),
                       var18=0.0 if pair else cfg['start_hour_angle_deg'] + (cfg['node_deg'] if tilted else
                                                                             cfg['site']['lon_deg'] if box else 0.0),
                       var19=0.0 if pair else round(solar_day_s(), 1),
                       **(dict(var20=1.0) if cfg.get('twilight') else {})),
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
    if cfg.get('elec') is not None:                                       # NSSL with charge, its CCN the case's droplets
        for section, entries in cm1_elec.namelist_settings(cfg['elec']).items():
            settings.setdefault(section, {}).update(entries)
        settings['nssl2mom_params'] = dict(ccn=cfg['droplets_cm3'] * 1.0e6)
    return settings


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


def setup_from(name: str, cfg: dict) -> Path:
    """Set up a case with another case's inputs as they were written and only its executable changed, so the pair
    differs in the build alone, whatever the setup has learned since. An electrified case also changes the namelist to
    the NSSL microphysics with the charge tracers and its electricity, and writes the files that needs."""
    source = RUNS / cfg['inputs_from']
    exe = CM1_HOME / 'build' / cfg['build'] / 'cm1.exe'
    if not exe.exists():
        raise RuntimeError(f'build {cfg["build"]} first')
    case = RUNS / name
    if (case / 'progress.json').exists():
        raise RuntimeError(f'{case} has started; remove it to set up afresh')
    case.mkdir(parents=True, exist_ok=True)
    for item in sorted(source.iterdir()):
        if item.name.startswith(('cm1out_', 'cm1rst_', 'cm1_segment_')) or item.name in (
                'progress.json', 'run.lock', 'STOP', 'namelist.input', 'cm1.exe', 'case.json'):
            continue
        target = case / item.name
        if target.is_symlink() or target.exists():
            target.unlink()
        if item.is_symlink():
            target.symlink_to(item.resolve())
        elif item.is_file():
            shutil.copy2(item, target)
    if (case / 'cm1.exe').is_symlink() or (case / 'cm1.exe').exists():
        (case / 'cm1.exe').unlink()
    (case / 'cm1.exe').symlink_to(exe)
    record = json.loads((source / 'case.json').read_text())
    build_record = json.loads((exe.parent / 'build.json').read_text())
    record.update(case=name, purpose=cfg['purpose'], build=build_record,
                  inputs_from=dict(case=cfg['inputs_from'], build=record['build']['label'],
                                   executable_sha256=record['build']['executable_sha256']))
    if cfg.get('namelist') or cfg.get('rain_limits'):  # output switches, and the rain limits, on the inputs as written
        text = (case / 'namelist.template').read_text()
        for section, entries in cfg.get('namelist', {}).items():
            for key, value in entries.items():
                text = set_namelist(text, section, key, value)
        (case / 'namelist.template').write_text(text)
        record.update(configuration=cfg, rain_limits=apply_rain_limits(case, cfg, build_record.get('gravity_m_s2')))
        record.setdefault('inputs', {})['namelist.template'] = hashlib.sha256(
            (case / 'namelist.template').read_bytes()).hexdigest()[:16]
    if cfg.get('elec') is not None:
        text = (case / 'namelist.template').read_text()
        settings = cm1_elec.namelist_settings(cfg['elec'])
        settings['nssl2mom_params'] = dict(ccn=cfg['droplets_cm3'] * 1.0e6)
        for section, entries in settings.items():
            for key, value in entries.items():
                text = set_namelist(text, section, key, value)
        text, files = cm1_elec.run_files(case, cfg['elec'], text)
        (case / 'namelist.template').write_text(text)
        record.update(configuration=cfg, electricity=dict(cm1_elec.SETTINGS, **cfg['elec'], files=files))
        record.setdefault('inputs', {})['namelist.template'] = hashlib.sha256(
            (case / 'namelist.template').read_bytes()).hexdigest()[:16]
    (case / 'case.json').write_text(json.dumps(record, indent=1) + '\n')
    return case


def fine_start(snapshot: dict) -> dict:
    """A finer box's starting air from one output of a coarser box: the means over its columns of potential
    temperature (K), vapour (kg/kg) and wind (m/s) by level, and at the surface its pressure (Pa), the potential
    temperature and vapour at 2 m and the skin temperature (K)."""
    import numpy as np
    mean = lambda name: np.asarray(snapshot[name], float).mean(axis=-1)
    kappa = 287.04 / 1005.7
    theta_2m = np.asarray(snapshot['t2'], float) * (1.0e5 / np.asarray(snapshot['psfc'], float)) ** kappa
    return dict(theta_k=mean('th'), qv_kg_kg=mean('qv'), u_m_s=mean('uinterp'), v_m_s=mean('vinterp'),
                surface_pa=float(mean('psfc')), theta_2m_k=float(theta_2m.mean()), qv_2m_kg_kg=float(mean('q2')),
                skin_k=float(mean('tsk')))


def fine_sounding(start: dict, zw) -> str:
    """CM1's input_sounding (isnd = 7) from a finer box's starting air: the surface line, then height (m), potential
    temperature (K), vapour (g/kg) and wind (m/s) at the scalar levels and at the model top, where the potential
    temperature carries on the slope of the top two levels."""
    import numpy as np
    zw = np.asarray(zw, float)
    zh = 0.5 * (zw[1:] + zw[:-1])
    th, qv, u, v = (np.asarray(start[k], float) for k in ('theta_k', 'qv_kg_kg', 'u_m_s', 'v_m_s'))
    if th.size != zh.size:
        raise ValueError(f'{th.size} levels of air for {zh.size} scalar levels')
    top = th[-1] + (th[-1] - th[-2]) / (zh[-1] - zh[-2]) * (zw[-1] - zh[-1])
    rows = [*zip(zh, th, qv, u, v), (zw[-1], top, qv[-1], u[-1], v[-1])]
    lines = [f"{start['surface_pa'] / 100:12.4f} {start['theta_2m_k']:12.4f} {start['qv_2m_kg_kg'] * 1000:12.5f}"]
    lines += [f'{z:12.3f} {t:12.4f} {q * 1000:12.5f} {uu:8.3f} {vv:8.3f}' for z, t, q, uu, vv in rows]
    return '\n'.join(lines) + '\n'


def fine_namelist(text: str, refine: int, day: float, cfg: dict) -> tuple:
    """A coarser box's namelist for a finer box: the spacing and the first time step divided by refine, the finer
    case's length, output and restarts, and the Sun's hour angle (var18, which also times the day-night forcing) at the
    coarser run's day it starts from. Returns the namelist and the finer box's spacing and starting hour angle."""
    dx = float(namelist_value(text, 'param1', 'dx')) / refine
    dtl = float(namelist_value(text, 'param1', 'dtl')) / refine
    h0, day_s = float(namelist_value(text, 'param8', 'var18')), float(namelist_value(text, 'param8', 'var19'))
    hour = (h0 + 360.0 * day * 86400.0 / day_s + 180.0) % 360.0 - 180.0
    for section, key, value in (('param1', 'dx', round(dx, 3)), ('param1', 'dy', round(dx, 3)),
                                ('param1', 'dtl', round(dtl, 3)), ('param1', 'timax', round(cfg['days'] * 86400.0)),
                                ('param1', 'run_time', -999.9), ('param1', 'tapfrq', float(cfg['output_s'])),
                                ('param1', 'rstfrq', float(cfg['restart_s'])), ('param2', 'irst', 0),
                                ('param2', 'rstnum', 1), ('param8', 'var18', round(hour, 4)),
                                ('param14', 'diagfrq', float(cfg['output_s']))):
        text = set_namelist(text, section, key, value)
    return text, dict(dx_m=dx, dtl_s=dtl, hour_angle_deg=hour)


def setup_fine(name: str, cfg: dict) -> Path:
    """Set up a finer box inside a coarser box's run (cfg['fine_from']: the coarser case, the day of its run to start
    from and the factor its spacing is divided by). The finer box keeps the coarser one's column count, vertical grid,
    nudging, vertical wind, day-night forcing, land and electricity, and starts from the coarser run's output at that
    day: its air averaged over the columns, by level, and its mean skin temperature over the land, under the Sun and
    the day-night forcing of that hour. CM1 cannot move a run onto a finer grid, so the finer box builds its own clouds
    from there."""
    import numpy as np
    from climate.crm import ring_analysis as ra
    fine = cfg['fine_from']
    if fine.get('day') is None:
        raise RuntimeError(f"{name}: choose the day of {fine['case']} to start from (fine_from['day']) first")
    source = RUNS / fine['case']
    src = json.loads((source / 'case.json').read_text())
    if src['configuration'].get('terrain'):
        raise RuntimeError(f"{fine['case']} has terrain, which a finer box does not carry")
    exe = CM1_HOME / 'build' / cfg['build'] / 'cm1.exe'
    if not exe.exists():
        raise RuntimeError(f'build {cfg["build"]} first')
    case = RUNS / name
    if (case / 'progress.json').exists():
        raise RuntimeError(f'{case} has started; remove it to set up afresh')
    output_s = src['configuration']['output_s']
    n = round(fine['day'] * 86400.0 / output_s) + 1
    if abs((n - 1) * output_s - fine['day'] * 86400.0) > 1.0 or not (source / f'cm1out_t{n:06d}_s.dat').exists():
        raise RuntimeError(f"{fine['case']} has no output at day {fine['day']}")
    case.mkdir(parents=True, exist_ok=True)
    start = fine_start(ra.read_snapshot(source, n))
    (case / 'input_sounding').write_text(fine_sounding(start, np.loadtxt(source / 'input_grid_z')))
    for item in ('input_grid_z', 'lsnudge_0001.dat', 'terluna_wls.txt', 'terluna_lsadv.txt', 'LANDUSE.TBL'):
        if (source / item).exists():
            shutil.copy2(source / item, case / item)
    segments = np.loadtxt(source / 'terluna_surface.txt', skiprows=1, ndmin=2)
    (case / 'terluna_surface.txt').write_text(f'{len(segments)}\n' + ''.join(
        f'{x0:.1f} {x1:.1f} {xl:.1f} {int(lu):d} {start["skin_k"] if xl == 1 else tsk:.3f} {tmn:.3f}\n'
        for x0, x1, xl, lu, tsk, tmn in segments))
    text, grid = fine_namelist((source / 'namelist.template').read_text(), fine['refine'], fine['day'], cfg)
    electricity = None
    if cfg.get('elec') is not None:
        for section, entries in cm1_elec.namelist_settings(cfg['elec']).items():
            for key, value in entries.items():
                text = set_namelist(text, section, key, value)
        text, files = cm1_elec.run_files(case, cfg['elec'], text)
        electricity = dict(cm1_elec.SETTINGS, **cfg['elec'], files=files)
    (case / 'namelist.template').write_text(text)
    tree = fetch()
    for link, target in (('cm1.exe', exe), ('RRTMG_LW_DATA', tree / 'run' / 'RRTMG_LW_DATA'),
                         ('RRTMG_SW_DATA', tree / 'run' / 'RRTMG_SW_DATA')):
        path = case / link
        if path.is_symlink() or path.exists():
            path.unlink()
        path.symlink_to(target)
    digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]
    lon = src['configuration'].get('site', {}).get('lon_deg', 0.0)
    record = dict(src, case=name, purpose=cfg['purpose'], inputs_from=None,
                  configuration=dict(cfg, start_hour_angle_deg=round(grid['hour_angle_deg'] - lon, 4)),
                  grid=dict(src['grid'], dx_m=grid['dx_m'], length_m=src['grid']['nx'] * grid['dx_m']),
                  build=json.loads((exe.parent / 'build.json').read_text()), electricity=electricity,
                  initial=dict(source=f"{fine['case']}'s output {n} (day {fine['day']}), averaged over its columns",
                               case=fine['case'], day=fine['day'], snapshot=n,
                               hour_angle_deg=round(grid['hour_angle_deg'], 4),
                               **{k: round(start[k], 4) for k in ('surface_pa', 'theta_2m_k', 'skin_k')},
                               qv_2m_kg_kg=round(start['qv_2m_kg_kg'], 7)),
                  runner=digest(__file__),
                  inputs={p: digest(case / p) for p in ('input_sounding', 'input_grid_z', 'lsnudge_0001.dat',
                                                         'terluna_surface.txt', 'LANDUSE.TBL', 'namelist.template',
                                                         'terluna_wls.txt', 'terluna_lsadv.txt')
                          if (case / p).exists()})
    (case / 'case.json').write_text(json.dumps(record, indent=1) + '\n')
    return case


RESTART_INPUTS = ('input_sounding', 'input_grid_z', 'lsnudge_0001.dat', 'terluna_wls.txt', 'terluna_lsadv.txt', TWILIGHT_TABLE,
                  'terluna_surface.txt', 'LANDUSE.TBL', 'cm1out_s.ctl', 'cm1out_stats.ctl', 'cm1out_metadata.ctl',
                  'cm1out_metadata.dat', RAIN_LAW_TABLE)


def restart_windows(source: Path, case: Path, start_s: float, end_s: float, bins: int = 36):
    """The wind-nudging windows of a case run on from the restart of a source that nudged its wind by local time
    (lsnudge_0001.dat onward, each window a bins-th of the solar day from the source's start, as setup writes them):
    the windows from the one holding start_s to the last starting before end_s, numbered from 1 as CM1 reads them,
    each with its own times and the profile of the source's window at the same phase of the day, the last
    open-ended. The source's windows at one phase must hold the same profile. Returns the record for case.json, or
    None when the source nudged toward one profile."""
    found = sorted(source.glob('lsnudge_0*.dat'))
    if len(found) < 2:
        return None
    step = solar_day_s() / bins
    rows = [p.read_text().splitlines() for p in found]
    for j, r in enumerate(rows):                                          # every source window at one phase alike
        if r[:1] + r[2:] != rows[j % bins][:1] + rows[j % bins][2:]:
            raise RuntimeError(f'{found[j].name} differs from {found[j % bins].name}, a lunar day earlier')
    first = int(math.floor(start_s / step + 1e-9))                        # 0-based index of the window holding start_s
    last = int(math.ceil(end_s / step - 1e-9)) - 1
    if last - first + 1 > 100:
        raise ValueError(f'{last - first + 1} nudging windows; CM1 reads at most 100')
    for k, i in enumerate(range(first, last + 1), 1):
        t1, t2 = i * step, (1.0e30 if i == last else (i + 1) * step)
        r = rows[i % bins]
        (case / f'lsnudge_{k:04d}.dat').write_text('\n'.join(r[:1] + [f'  {t1:22.1f} {t2:22.6e}'] + r[2:]) + '\n')
    return dict(windows=last - first + 1, window_s=step, first_window_of_source=first + 1,
                source=f"{source.name}'s windows by local time, carried on by whole lunar days (restart_windows)")


def window_root(source: Path) -> Path:
    """The case whose wind-nudging windows a run on from source's restart carries: source itself, or, when source was
    itself run on from another case's restart, the case the chain began from. A run-on case numbers its windows from
    its own restart, so restart_windows, which counts phases from the run's start, would take them a few windows out
    of phase (one window, 0.82 days, for a box's third lunar day run on to a fourth)."""
    while True:
        earlier = json.loads((source / 'case.json').read_text()).get('restart_from')
        if not earlier:
            return source
        source = source.parent / earlier['case']


def setup_restart(name: str, cfg: dict) -> Path:
    """Set up a case that runs on from another case's restart (cfg['restart_from']: the case and the day), so that a
    window of a run can run again under other settings, or run on: the other case's inputs, restart files and kept
    vertical field, its rain law's table and its wind-nudging windows carried over the case's span (restart_windows),
    and its namelist with this case's electricity. Its progress starts at the restart, and the runner takes up from
    there; output and restarts keep the other case's numbering."""
    fr = cfg['restart_from']
    source = RUNS / fr['case']
    src = json.loads((source / 'case.json').read_text())
    scale = src.get('time_scale', 1.0)
    restart_s = cfg['restart_s'] * scale
    start_s = fr['day'] * 86400.0 * scale
    n = round(start_s / restart_s)
    files = [source / f'cm1rst_t{n:06d}_{part}.dat' for part in 'isuvwx']
    if abs(n * restart_s - start_s) > 1.0 or not all(f.exists() for f in files):
        raise RuntimeError(f"{fr['case']} has no restart at day {fr['day']}")
    ez = source / f'terluna_ez_{round(n * restart_s)}.bin'
    if cfg.get('elec') is not None and not ez.exists():
        raise RuntimeError(f'{ez} is missing: the restart\'s vertical field for the inductive charging')
    exe = CM1_HOME / 'build' / cfg['build'] / 'cm1.exe'
    if not exe.exists():
        raise RuntimeError(f'build {cfg["build"]} first')
    case = RUNS / name
    if (case / 'progress.json').exists():
        raise RuntimeError(f'{case} has been set up; remove it to set up afresh')
    case.mkdir(parents=True, exist_ok=True)
    for item in RESTART_INPUTS:
        if (source / item).exists():
            shutil.copy2(source / item, case / item)
    for f in files + ([ez] if ez.exists() else []):
        shutil.copy2(f, case / f.name)
    windows = restart_windows(window_root(source), case, n * restart_s, case_length_s(cfg, scale))
    if windows and abs(windows['window_s'] - src['nudging_wind']['window_s']) > 1e-3:
        raise RuntimeError(f"{fr['case']}'s windows are {src['nudging_wind']['window_s']} s long, not {windows['window_s']}")
    text = (source / 'namelist.template').read_text()
    if abs(cfg['output_s'] - src['configuration']['output_s']) > 0.5:   # a window may write output more often;
        text = set_namelist(text, 'param1', 'tapfrq', cfg['output_s'] * scale)   # CM1 resets its schedule on restart
    for section, entries in (cfg.get('restart_namelist') or {}).items():    # the case's own changes to the source's
        for key, value in entries.items():                                   # namelist (no earlier case sets any)
            text = set_namelist(text, section, key, value)
    electricity = None
    if cfg.get('elec') is not None:
        for section, entries in cm1_elec.namelist_settings(cfg['elec']).items():
            for key, value in entries.items():
                text = set_namelist(text, section, key, value)
        text, efiles = cm1_elec.run_files(case, cfg['elec'], text)
        electricity = dict(cm1_elec.SETTINGS, **cfg['elec'], files=efiles)
    (case / 'namelist.template').write_text(text)
    build_record = json.loads((exe.parent / 'build.json').read_text())
    rain = apply_rain_limits(case, cfg, build_record.get('gravity_m_s2'))
    tree = fetch()
    for link, target in (('cm1.exe', exe), ('RRTMG_LW_DATA', tree / 'run' / 'RRTMG_LW_DATA'),
                         ('RRTMG_SW_DATA', tree / 'run' / 'RRTMG_SW_DATA')):
        path = case / link
        if path.is_symlink() or path.exists():
            path.unlink()
        path.symlink_to(target)
    made = next((s for s in _progress(source)['segments'] if s['model_s'] >= n * restart_s - 1.0), {})
    digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]
    windows_written = tuple(f'lsnudge_{k:04d}.dat' for k in range(2, windows['windows'] + 1)) if windows else ()
    record = dict(src, case=name, purpose=cfg['purpose'], inputs_from=None, configuration=dict(cfg),
                  build=build_record, electricity=electricity, rain_limits=rain,
                  restart_from=dict(case=fr['case'], day=fr['day'], restart=n, model_s=n * restart_s,
                                    made_by=made.get('executable')),
                  **({'nudging_wind': dict(src.get('nudging_wind') or {}, **windows)} if windows else {}),
                  runner=digest(__file__),
                  inputs={p: digest(case / p) for p in RESTART_INPUTS + ('namelist.template',) + windows_written
                          if (case / p).exists()})
    (case / 'case.json').write_text(json.dumps(record, indent=1) + '\n')
    _save_progress(case, dict(segments=[], case=name, start_s=n * restart_s,
                              start=f"{fr['case']}'s restart {n} (day {fr['day']})"))
    return case


NIGHT_INPUTS = ('input_grid_z', 'terluna_wls.txt', 'terluna_lsadv.txt', 'terluna_surface.txt', 'LANDUSE.TBL',
                TWILIGHT_TABLE, RAIN_LAW_TABLE)


def hour_angle_at(h0_deg: float, t_s: float, day_s: float) -> float:
    """The hour angle (degrees, -180 to 180, 0 = noon) at model time t_s of a run that starts at h0_deg."""
    return (h0_deg + 360.0 * t_s / day_s + 180.0) % 360.0 - 180.0


def sun_elevation_deg(latitude_deg: float, hour_angle_deg: float) -> float:
    """The Sun's elevation as the Terluna Sun puts it, at the equinox: sin(elevation) = cos(latitude) cos(hour angle)."""
    return math.degrees(math.asin(math.cos(math.radians(latitude_deg)) * math.cos(math.radians(hour_angle_deg))))


def held_wind_window(source: Path, t_s: float) -> tuple:
    """The source's wind-nudging window that holds model time t_s, as CM1 steps through them (on to the next once the
    time passes a window's end): its file and its lines."""
    for path in sorted(source.glob('lsnudge_0*.dat')):
        rows = path.read_text().splitlines()
        t1, t2 = (float(v) for v in rows[1].split())
        if t1 <= t_s <= t2:
            return path, rows
    raise RuntimeError(f'{source.name} has no wind-nudging window holding model time {t_s:.0f} s')


def setup_night(name: str, cfg: dict) -> Path:
    """Set up a box that starts from another box's settled air at one hour and holds that hour's forcing
    (cfg['night_from']: the case, the day of its run and the factor its spacing is divided by). The air is the
    source's output at that day averaged over its columns, by level (fine_start), with CM1's random perturbations; the
    clock that times the Sun and the day-night table (var19) stops at that hour (var18), so both hold their values;
    the wind is nudged toward the source's window at that hour for the whole run; and every other input is the
    source's as written: the hold on the air above the GCM's convection, the mean vertical wind, the surface, the
    twilight, the rain law and limits. The box takes the case's columns, the source's spacing and first time step
    divided by refine, and the case's length, output, diagnostics, statistics and restarts."""
    import numpy as np
    from climate.crm import ring_analysis as ra
    night = cfg['night_from']
    source = RUNS / night['case']
    src = json.loads((source / 'case.json').read_text())
    if src['configuration'].get('terrain'):
        raise RuntimeError(f"{night['case']} has terrain, which a box started from its mean air does not carry")
    exe = CM1_HOME / 'build' / cfg['build'] / 'cm1.exe'
    if not exe.exists():
        raise RuntimeError(f'build {cfg["build"]} first')
    case = RUNS / name
    if (case / 'progress.json').exists():
        raise RuntimeError(f'{case} has started; remove it to set up afresh')
    output_s = src['configuration']['output_s']
    n = round(night['day'] * 86400.0 / output_s) + 1
    t_s = (n - 1) * output_s
    if abs(t_s - night['day'] * 86400.0) > 1.0 or not (source / f'cm1out_t{n:06d}_s.dat').exists():
        raise RuntimeError(f"{night['case']} has no output at day {night['day']}")
    window, rows = held_wind_window(source, t_s)
    case.mkdir(parents=True, exist_ok=True)
    start = fine_start(ra.read_snapshot(source, n))
    (case / 'input_sounding').write_text(fine_sounding(start, np.loadtxt(source / 'input_grid_z')))
    for item in NIGHT_INPUTS:
        if (source / item).exists():
            shutil.copy2(source / item, case / item)
    segments = np.loadtxt(source / 'terluna_surface.txt', skiprows=1, ndmin=2)
    if (segments[:, 2] == 1).any():                                       # land takes the mean skin, as a finer box's
        (case / 'terluna_surface.txt').write_text(f'{len(segments)}\n' + ''.join(
            f'{x0:.1f} {x1:.1f} {xl:.1f} {int(lu):d} {start["skin_k"] if xl == 1 else tsk:.3f} {tmn:.3f}\n'
            for x0, x1, xl, lu, tsk, tmn in segments))
    (case / 'lsnudge_0001.dat').write_text('\n'.join(rows[:1] + [f'  {0.0:22.1f} {1.0e30:22.6e}'] + rows[2:]) + '\n')
    source_text = (source / 'namelist.template').read_text()
    text, grid = fine_namelist(source_text, night['refine'], night['day'], cfg)
    for section, key, value in (('param0', 'nx', cfg['nx']), ('param0', 'ny', cfg['ny']),
                                ('param1', 'statfrq', float(cfg['stat_s'])), ('param14', 'diagfrq', float(cfg['diag_s'])),
                                ('param8', 'var19', STOPPED_CLOCK_S)):
        text = set_namelist(text, section, key, value)
    for section, entries in (cfg.get('namelist') or {}).items():          # the case's output switches
        for key, value in entries.items():
            text = set_namelist(text, section, key, value)
    (case / 'namelist.template').write_text(text)
    tree = fetch()
    for link, target in (('cm1.exe', exe), ('RRTMG_LW_DATA', tree / 'run' / 'RRTMG_LW_DATA'),
                         ('RRTMG_SW_DATA', tree / 'run' / 'RRTMG_SW_DATA')):
        path = case / link
        if path.is_symlink() or path.exists():
            path.unlink()
        path.symlink_to(target)
    digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]
    site = src['configuration']['site']
    hour = grid['hour_angle_deg']
    h0, day_s = float(namelist_value(source_text, 'param8', 'var18')), float(namelist_value(source_text, 'param8', 'var19'))
    t1, t2 = (float(v) for v in rows[1].split())
    record = dict(src, case=name, purpose=cfg['purpose'], inputs_from=None, restart_from=None,
                  configuration=dict(cfg, start_hour_angle_deg=round(hour - site['lon_deg'], 4)),
                  grid=dict(src['grid'], nx=cfg['nx'], ny=cfg['ny'], dx_m=grid['dx_m'], length_m=cfg['nx'] * grid['dx_m']),
                  build=json.loads((exe.parent / 'build.json').read_text()),
                  initial=dict(source=f"{night['case']}'s output {n} (day {night['day']}), averaged over its columns, "
                                      "with CM1's random potential-temperature perturbations of up to 0.25 K (irandp 1)",
                               case=night['case'], day=night['day'], snapshot=n,
                               **{k: round(start[k], 4) for k in ('surface_pa', 'theta_2m_k', 'skin_k')},
                               qv_2m_kg_kg=round(start['qv_2m_kg_kg'], 7)),
                  held=dict(hour_angle_deg=round(hour, 4), sun_elevation_deg=round(sun_elevation_deg(site['lat_deg'], hour), 2),
                            clock_s=STOPPED_CLOCK_S, day_night_table=f"{night['case']}'s, at that hour",
                            wind_window=dict(source_file=window.name, source_time_s=[t1, t2],
                                             hour_angle_deg=[round(hour_angle_at(h0, t1, day_s), 3),
                                                             round(hour_angle_at(h0, t2, day_s), 3)])),
                  nudging_wind=dict(windows=1, source=f"{night['case']}'s {window.name}, held for the whole run"),
                  runner=digest(__file__),
                  inputs={p: digest(case / p) for p in ('input_sounding', 'lsnudge_0001.dat', 'namelist.template',
                                                         *NIGHT_INPUTS) if (case / p).exists()})
    (case / 'case.json').write_text(json.dumps(record, indent=1) + '\n')
    return case


def namelist_value(text: str, section: str, key: str) -> str:
    """One entry of a Fortran namelist held as text, as written."""
    import re
    start = text.index(f'&{section}\n')
    block = text[start:text.index('\n /', start)]
    return re.search(rf'^\s*{re.escape(key)}\s*=\s*([^,\n!]*)', block, re.M).group(1).strip()


def setup_sample(name: str, cfg: dict) -> Path:
    """Set up one of CM1's own sample cases (run/config_files/<sample>) with the case's changes to its namelist and,
    when it is electrified, its electricity."""
    tree = fetch()
    exe = CM1_HOME / 'build' / cfg['build'] / 'cm1.exe'
    if not exe.exists():
        raise RuntimeError(f'build {cfg["build"]} first')
    case = RUNS / name
    if (case / 'progress.json').exists():
        raise RuntimeError(f'{case} has started; remove it to set up afresh')
    case.mkdir(parents=True, exist_ok=True)
    text = (tree / 'run' / 'config_files' / cfg['sample'] / 'namelist.input').read_text()
    settings = {'param1': dict(timax=round(cfg['days'] * 86400.0, 3), run_time=-999.9, tapfrq=cfg['output_s'],
                               rstfrq=cfg['restart_s']),
                'param16': dict(restart_format=1, restart_filetype=2, restart_reset_frqtim=True)}
    for section, entries in cfg.get('namelist', {}).items():
        settings.setdefault(section, {}).update(entries)
    if cfg.get('elec') is not None:
        for section, entries in cm1_elec.namelist_settings(cfg['elec']).items():
            settings.setdefault(section, {}).update(entries)
    for section, entries in settings.items():
        for key, value in entries.items():
            text = set_namelist(text, section, key, value)
    electricity = None
    if cfg.get('elec') is not None:
        text, files = cm1_elec.run_files(case, cfg['elec'], text)
        electricity = dict(cm1_elec.SETTINGS, **cfg['elec'], files=files)
    (case / 'namelist.template').write_text(text)
    if (case / 'cm1.exe').is_symlink() or (case / 'cm1.exe').exists():
        (case / 'cm1.exe').unlink()
    (case / 'cm1.exe').symlink_to(exe)
    value = lambda section, key: float(namelist_value(text, section, key))
    grid = dict(nx=int(value('param0', 'nx')), ny=int(value('param0', 'ny')), nz=int(value('param0', 'nz')),
                dx_m=value('param1', 'dx'), dy_m=value('param1', 'dy'), dz_m=value('param1', 'dz'))
    record = dict(case=name, configuration=cfg, time_scale=1.0, sample=cfg['sample'], grid=grid,
                  build=json.loads((exe.parent / 'build.json').read_text()), electricity=electricity,
                  runner=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16],
                  inputs={'namelist.template': hashlib.sha256((case / 'namelist.template').read_bytes()).hexdigest()[:16]})
    (case / 'case.json').write_text(json.dumps(record, indent=1) + '\n')
    return case


def write_twilight(path: Path, mode: str, product: str = None) -> dict:
    """terluna_twilight.txt from the twilight product's shells: the elevation and shell counts, the shells' height
    edges (m), then for each solar elevation (degrees) what the curved air adds to the ground and to each shell, per
    unit sunlight at the beam; then the extra infrared each shell and the ground lose to space (W/m2) and the
    temperatures it is given at. 'full' adds the shells' widening to both; with the product's identity."""
    if mode not in TWILIGHT_MODES:
        raise ValueError(f'twilight must be one of {TWILIGHT_MODES}')
    product = HERE.parent / 'results' / 'gcm' / (product or TWILIGHT_PRODUCT)
    data = json.loads(product.read_text())
    term, cone, lw = data['cm1_term'], data['cm1_term_cone'], data['longwave']
    full = mode == 'full'
    edges = [1e3 * z for z in term['height_edges_km']]
    rows = [f"{len(term['elevation_deg'])} {len(edges) - 1}", ' '.join(f'{z:.3f}' for z in edges)]
    for i, el in enumerate(term['elevation_deg']):
        shells = [a[i] + (b[i] if full else 0.0) for a, b in zip(term['shells'], cone['shells'])]
        rows.append(' '.join(f'{x:.6e}' for x in [el, term['ground'][i]] + shells))
    escape = [a + (b if full else 0.0) for a, b in zip(lw['shells']['curved'], lw['shells']['cone'])]
    rows.append(' '.join(f'{x:.6e}' for x in escape + [lw['ground_curved_w_m2']]))
    rows.append(' '.join(f'{x:.6e}' for x in lw['shells']['temperature_k'] + [lw['ground_temperature_k']]))
    path.write_text('\n'.join(rows) + '\n')
    return dict(product=str(product.relative_to(HERE.parents[1])), sha256=hashlib.sha256(product.read_bytes()).hexdigest()[:16],
                mode=mode, elevations=len(term['elevation_deg']), shells=len(edges) - 1)


def held_share_above(zg, theta_at, qv_below, surface_pa, gravity, ztop):
    """Vapour (kg/kg) at any height: the GCM's own up to its top level zg[-1], and above it whatever keeps that level's
    share of ice saturation as the temperature falls. The pressure comes from the potential temperature profile
    hydrostatically, at the build's gravity over the case's heights, as CM1 takes it from the sounding."""
    import numpy as np
    from atmosphere.radiative_convective.thermodynamics import saturation_pressure
    cp, rd = 1005.7, 287.04
    z = np.linspace(0.0, max(ztop, zg[-1]) * 1.01, 4001)
    q = qv_below(z)
    inv_tv = 1.0 / (theta_at(z) * (1.0 + 0.61 * q))
    exner = (surface_pa / 1e5) ** (rd / cp) - gravity / cp * np.concatenate(
        [[0.0], np.cumsum(0.5 * (inv_tv[1:] + inv_tv[:-1]) * np.diff(z))])
    p, t = 1e5 * exner ** (cp / rd), theta_at(z) * exner
    e = q * p / (0.622 + q)
    share = float(np.interp(zg[-1], z, e / saturation_pressure(t)))
    e_held = share * saturation_pressure(t)
    log_held = np.log(0.622 * e_held / (p - e_held))
    return lambda zz: np.where(np.asarray(zz) <= zg[-1], qv_below(zz), np.exp(np.interp(zz, z, log_held)))


def setup(name: str) -> Path:
    """Write a case's inputs into climate/crm/runs/<name>: namelist, sounding, grid, nudging profile,
    surface segments, land-use table and links to the executable and radiation tables."""
    import numpy as np
    cfg = CASES[name]
    if cfg.get('night_from'):
        return setup_night(name, cfg)
    if cfg.get('restart_from'):
        return setup_restart(name, cfg)
    if cfg.get('fine_from'):
        return setup_fine(name, cfg)
    if cfg.get('inputs_from'):
        return setup_from(name, cfg)
    if cfg.get('kind') == 'sample':
        return setup_sample(name, cfg)
    tree = fetch()
    exe = CM1_HOME / 'build' / cfg['build'] / 'cm1.exe'
    if not exe.exists():
        raise RuntimeError(f'build {cfg["build"]} first')
    build_record = json.loads((exe.parent / 'build.json').read_text())
    air, gravity = build_record['air'], build_record['gravity_m_s2']
    run, years = cfg.get('gcm', (GCM_RUN, GCM_YEARS))
    if run != air.get('run'):                     # the run's sunlight, in a build compiled with the same air
        own = design_air(run)
        if abs(own['co2'] - air['co2']) > 1e-12 or abs(own['o2'] - air['o2']) > 1e-12:
            raise RuntimeError(f"build {cfg['build']} holds the air of {air.get('run')}, not of {run}")
        air = dict(air, sunlight_w_m2=own['sunlight_w_m2'], run=run)
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
        host = {} if cfg['site'].get('gcm') else CASES[cfg['site']['ring']]   # a tilted ring's, or none: the equator's
        nx_ring, dx, _ = ring_grid(0.0 if host.get('kind') == 'tilted' else host.get('latitude_deg', 0.0),
                                   cfg['dx_target_m'], cfg['ranks'])
        nx, length = cfg['nx'], cfg['nx'] * dx
    else:
        nx, dx, length = ring_grid(latitude, cfg['dx_target_m'], cfg['ranks'])
    zw = [z * s for z in vertical_grid(**cfg['z_grid'])]
    nz = len(zw) - 1
    track = tilted_path(nx, cfg['tilt_deg'], cfg['node_deg']) if tilted else None
    site = box_site(cfg, zw, nx_ring, run, years) if box else None
    if box and cfg['nudge'].get('held_from') == 'convection_top':        # the air held above the GCM's convection only
        cfg = dict(cfg, nudge=held_heights(cfg['nudge'], site['forcing']['convection'], zw[-1], s))
    ref = (gcm_path_reference(track['lat'], track['lon'], cfg['tilt_deg'], zw, run, years) if tilted else site['ref']
           if box else gcm_equator_profile(run, years, band=cfg.get('band_deg', EQUATOR_BAND_DEG), latitude=latitude))
    zh = 0.5 * (np.array(zw[1:]) + np.array(zw[:-1]))
    prof = ref['profile']
    zg = np.array(prof['z_m'][::-1]) * s                                  # GCM levels, bottom first
    th = np.array(prof['theta_k'][::-1])
    qv = np.array(prof['qv_kg_kg'][::-1])
    uu = np.array(prof['u_m_s'][::-1]) * (0.0 if pair else 1.0)
    vv = np.array(prof.get('v_m_s', [0.0] * len(prof['u_m_s']))[::-1]) * (0.0 if pair else 1.0)
    top_slope = (th[-1] - th[-2]) / (zg[-1] - zg[-2])
    theta_at = lambda z: np.where(z <= zg[-1], np.interp(z, zg, th), th[-1] + top_slope * (z - zg[-1]))
    qv_below = lambda z: np.exp(np.interp(z, zg, np.log(np.maximum(qv, 1e-9))))
    # Above the GCM's top level the vapour falls with the temperature, keeping that level's share of ice saturation.
    # Holding the vapour itself left the top of ring_45n_corrected_sideways 6-15 times saturated at the start, and the
    # ice that formed there at once stayed at 1.1 g/m2 for the whole run (2026-10-17).
    qv_at = held_share_above(zg, theta_at, qv_below, ref['surface_pa_sea'], gravity, zw[-1])
    u_at = lambda z: np.interp(z, zg, uu)
    v_at = lambda z: np.interp(z, zg, vv)
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
    lines += [f'{z:12.3f} {float(start_theta(z)):12.4f} {float(start_qv(z)) * 1000:12.5f} {float(u_at(z)):8.3f} '
              f'{float(v_at(z)):8.3f}'
              for z in [*zh, zw[-1]]]
    (case / 'input_sounding').write_text('\n'.join(lines) + '\n')
    (case / 'input_grid_z').write_text(''.join(f'{z:.3f}\n' for z in zw))
    nudge = [' *  Header:   lsnudge_time1 (s)   lsnudge_time2 (s)', '                     0.0                1.0e30 ',
             ' *  Profile.   Note: values will be ignored if nudging for that variable is off.',
             ' *    z (m)   theta (K)    qv (g/kg)  u (m/s, grnd-reltv) v (m/s, grnd-reltv)']
    nudge += [f'  {z:12.4f} {float(theta_at(z)):12.4f} {float(qv_at(z)) * 1000:12.6f} {float(u_at(z)):10.4f} '
              f'{float(v_at(z)):10.4f}'
              for z in [0.0, *zh, zw[-1]]]
    nudging_wind = None
    if box and cfg['site'].get('gcm') and cfg.get('wind_by_local_time'):   # the GCM's day-night wind, window by window
        from climate.crm import gcm_site
        heights = np.array([0.0, *zh, zw[-1]])
        windows = gcm_site.wind_windows(site['forcing'], heights, cfg['start_hour_angle_deg'] + cfg['site']['lon_deg'],
                                        solar_day_s(), case_length_s(cfg, s))
        if len(windows) > 100:
            raise ValueError(f'{len(windows)} nudging windows; CM1 reads at most 100')
        for i, window in enumerate(windows, 1):
            t1, t2 = window['time_s']
            rows = nudge[:1] + [f'  {t1:22.1f} {1.0e30 if i == len(windows) else t2:22.6e}'] + nudge[2:4]
            rows += [f'  {z:12.4f} {float(theta_at(z)):12.4f} {float(qv_at(z)) * 1000:12.6f} {u:10.4f} {v:10.4f}'
                     for z, u, v in zip(heights, window['u_m_s'], window['v_m_s'])]
            (case / f'lsnudge_{i:04d}.dat').write_text('\n'.join(rows) + '\n')
        nudging_wind = dict(windows=len(windows), window_s=windows[0]['time_s'][1] - windows[0]['time_s'][0],
                            source="the GCM's wind by local time at the site (gcm_site.wind_windows)",
                            lowest_level_speed_m_s=dict(min=float(min(np.hypot(w['u_m_s'][1], w['v_m_s'][1])
                                                                      for w in windows)),
                                                        max=float(max(np.hypot(w['u_m_s'][1], w['v_m_s'][1])
                                                                      for w in windows))))
    else:
        (case / 'lsnudge_0001.dat').write_text('\n'.join(nudge) + '\n')
    if pair:
        segments, water_share, surface_sha = [(-1e9, 1e9, 2, 16, ref['sea_surface_k'], ref['sea_surface_k'])], 1.0, None
    elif box:
        segments = [(-1e9, 1e9, 1, wetness_class(site['wetness']), site['ground_k'], site['ground_k'])]
        water_share, surface_sha = 0.0, site['surface_product']
        if cfg['site'].get('surface') == 'sea':                              # all of it sea, held at the GCM's
            segments, water_share = [(-1e9, 1e9, 2, 16, site['sea_k'], site['sea_k'])], 1.0
        if cfg.get('terrain'):                    # the ground at its heights: all of it land, or with the atlas's sea
            t = cfg['terrain']
            blat, blon = box_columns(cfg['site']['lat_deg'], cfg['site']['lon_deg'], site['heading_deg'], nx, cfg['ny'], dx)
            raw, wet = atlas_heights(blat.ravel(), blon.ravel(), t['half_width_deg'])
            sea = bool(t.get('sea'))
            wet = wet.reshape(blat.shape) & sea
            shape = box_heights(raw.reshape(blat.shape), dx, t['passes'], t['taper_km'] * 1000.0, wet if sea else None)
            heights = shape.pop('heights')
            np.asarray(heights, dtype='<f4').tofile(case / 'perts.dat')   # CM1's itern = 4: nx by ny reals, x fastest
            ground = site['ground_k'] - site['lapse_k_m'] * heights       # carried from sea level to each column's height
            water = (site['sea_k'] - site['lapse_k_m'] * heights) if sea else ground   # the sea, and raised lakes
            row = wetness_class(site['wetness'])
            (case / 'terluna_surface2d.txt').write_text(f'{nx} {cfg["ny"]}\n' + ''.join(
                f'{s:.3f} {s:.3f} 16\n' if w else f'{g:.3f} {g:.3f} {row:d}\n'
                for g, s, w in zip(ground.ravel(), water.ravel(), wet.ravel())))
            water_share = float(wet.mean())
            terrain = dict(t, source='the 28% scenario atlas, 0.25-degree cells', heading_deg=site['heading_deg'],
                           water_share=float(wet.mean()), lapse_k_m=site['lapse_k_m'], **shape)
    elif tilted:
        sea_k, land_k = ref['sea_k'], ref['land_k']                     # carried to sea level
        if cfg.get('terrain'):                                            # the ground and raised lakes at their heights
            heights = path_heights(track['lat'], track['lon'], passes=cfg['terrain']['passes'])
            sea_k, land_k = sea_k - ref['lapse_k_m'] * heights, land_k - ref['lapse_k_m'] * heights
            np.asarray(heights, dtype='<f4').tofile(case / 'perts.dat')   # CM1's itern = 4: nx by ny reals
            slope = np.abs(np.diff(np.append(heights, heights[0]))) / dx
            terrain = dict(cfg['terrain'], source='the 28% scenario atlas, 0.5-degree boxes', max_m=float(heights.max()),
                           mean_m=float(heights.mean()), max_slope=float(slope.max()), p99_slope=float(np.percentile(slope, 99)))
        segments, water_share, surface_sha = path_surface(track['lat'], track['lon'], dx, sea_k, land_k, ref['wetness'])
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
    elif cfg.get('land'):                         # a circle of latitude's one land row, its wetness as written above
        row = landuse_rows(((index, round(land_moisture(cfg, ref), 2)),), cfg['land'])[index]
        table = [row if line.startswith(f'{index},') else line for line in table]
    (case / 'LANDUSE.TBL').write_text(''.join(table))
    vertical_wind = None
    if box and cfg.get('large_scale_w') == 'gcm':
        vertical_wind = dict(z_m=[float(z) for z in zw], w_m_s=site['w_m_s'], ring=cfg['site'].get('ring'),
                             ring_columns=site['columns'], rows_deg=ref['rows_deg'])
        write_vertical_wind(case / 'terluna_wls.txt', vertical_wind, zw[-1])
    elif tilted and cfg.get('large_scale_w') == 'gcm':
        write_vertical_wind_2d(case / 'terluna_wls2d.txt', zw, ref['w_columns'])
        wc = np.asarray(ref['w_columns'])
        vertical_wind = dict(z_m=[float(z) for z in zw], path_mean_m_s=[float(v) for v in wc.mean(axis=0)],
                             min_m_s=float(wc.min()), max_m_s=float(wc.max()), rows_deg=ref['rows_deg'])
    elif cfg.get('large_scale_w') == 'gcm':
        vertical_wind = gcm_vertical_wind(latitude, cfg.get('band_deg', EQUATOR_BAND_DEG), run, years)
        write_vertical_wind(case / 'terluna_wls.txt', vertical_wind, zw[-1])
    sideways = None
    if cfg.get('sideways') and not (box or tilted):                      # the GCM's transport between latitudes
        forcing = sideways_forcing(run, years, latitude, cfg.get('band_deg', EQUATOR_BAND_DEG), ref)
        write_day_night(case / 'terluna_lsadv.txt', forcing)
        low = forcing['z_m'] < 2000.0
        sideways = dict(rows_deg=forcing['rows_deg'], lowest_2_km=dict(
            theta_k_day=float(forcing['theta_k_s'][0, low].mean() * 86400.0),
            qv_g_kg_day=float(forcing['qv_kg_kg_s'][0, low].mean() * 8.64e7)))
    day_night = None
    if box and cfg.get('day_night') == 'gcm':                            # the GCM's own budget at the site
        from climate.crm import gcm_site
        forcing = gcm_site.to_cm1(site['forcing'], zw, theta_at, qv_at)['day_night']
        write_day_night(case / 'terluna_lsadv.txt', forcing)
        low = forcing['z_m'] < 3000.0
        spread = lambda a, scale: dict(mean=float(a[:, low].mean() * scale), min=float(a[:, low].mean(axis=1).min() * scale),
                                       max=float(a[:, low].mean(axis=1).max() * scale))
        f = site['forcing']
        day_night = dict(source="the GCM's own budget at the site (gcm_site.py)", run=f['run'], years=f['years'],
                         harmonics=f['harmonics'], top_m=f['day_night']['top_m'],
                         lowest_3_km=dict(theta_k_day=spread(forcing['theta_k_s'], 86400.0),
                                          qv_g_kg_day=spread(forcing['qv_kg_kg_s'], 8.64e7)))
    elif box and cfg.get('day_night'):                                    # the rings' day-night circulation there
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
    for section, entries in (cfg.get('namelist') or {}).items():          # a case's own output switches
        for key, value in entries.items():
            text = set_namelist(text, section, key, value)
    electricity = None
    if cfg.get('elec') is not None:
        text, electricity = cm1_elec.run_files(case, cfg['elec'], text)
        electricity = dict(cm1_elec.SETTINGS, **cfg['elec'], files=electricity)
    (case / 'namelist.template').write_text(text)
    rain = apply_rain_limits(case, cfg, gravity)
    rain_law = apply_rain_law(case, cfg)
    twilight = (write_twilight(case / TWILIGHT_TABLE, cfg['twilight'], cfg.get('twilight_product', TWILIGHT_PRODUCT))
                if cfg.get('twilight') else None)
    digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]
    if tilted or box:                                                     # per-column arrays as rounded lists
        ref = {k: ([round(float(x), 6 if k == 'lapse_k_m' else 4) for x in v] if isinstance(v, np.ndarray) else v)
               for k, v in ref.items() if k != 'w_columns'}
    record = dict(case=name, configuration=cfg, time_scale=s,
                  grid=dict(nx=nx, dx_m=dx, length_m=length, nz=nz, ztop_m=zw[-1], **(dict(ny=cfg['ny']) if box else {})),
                  site=dict(cfg['site'], ring_columns=site['columns'], land_share=site['land_share'], wetness=site['wetness'],
                            wetness_class=wetness_class(site['wetness']), ground_k=site['ground_k']) if box else None,
                  build=build_record, gcm=dict(run=run, years=list(years)), air=air,
                  reference=ref, initial=initial, water_share=water_share, surface_product=surface_sha,
                  path=dict(tilt_deg=cfg['tilt_deg'], node_deg=cfg['node_deg'], lat=[round(float(v), 4) for v in track['lat']],
                            lon=[round(float(v), 4) for v in track['lon']], dlon=[round(float(v), 4) for v in track['dlon']])
                  if tilted else None,
                  latitude_deg=latitude, coriolis_1_s=turning_rate(cfg, latitude),
                  land_moisture=None if pair else site['wetness'] if box else land_moisture(cfg, ref),
                  large_scale_w=vertical_wind,
                  day_night=day_night, sideways=sideways,
                  terrain=terrain if cfg.get('terrain') else None,
                  electricity=electricity,
                  twilight=twilight,
                  rain_limits=rain, rain_law=rain_law, nudging_wind=nudging_wind,
                  gcm_site=gcm_site_record(site['forcing']) if box and site.get('forcing') else None,
                  runner=digest(__file__),
                  inputs={p: digest(case / p) for p in ('input_sounding', 'input_grid_z', 'lsnudge_0001.dat',
                                                         'terluna_surface.txt', 'LANDUSE.TBL', 'namelist.template',
                                                         'terluna_wls.txt', 'terluna_wls2d.txt', 'terluna_lsadv.txt',
                                                         'perts.dat', 'terluna_surface2d.txt', TWILIGHT_TABLE,
                                                         RAIN_LAW_TABLE, *sorted(q.name for q in case.glob('lsnudge_0*.dat'))[1:])
                          if (case / p).exists()})
    (case / 'case.json').write_text(json.dumps(record, indent=1) + '\n')
    return case


def latest_restart(case: Path):
    """Index of the newest complete restart set (cm1rst_tNNNNNN_{i,s,u,v,w,x}.dat), or 0."""
    import re
    found = {int(m.group(1)) for p in case.glob('cm1rst_t*_s.dat') if (m := re.match(r'cm1rst_t(\d+)_s\.dat$', p.name))}
    complete = [n for n in found if all((case / f'cm1rst_t{n:06d}_{part}.dat').exists() for part in 'isuvwx')]
    return max(complete, default=0)


def prune_restarts(case: Path, keep: int, start: int = 0) -> list:
    """Remove the complete restart sets a case has written (numbers above start, the set it started from) beyond its
    newest keep; returns the numbers removed. Incomplete sets, the starting set and the newest keep stay."""
    import re
    found = {int(m.group(1)) for p in case.glob('cm1rst_t*_s.dat') if (m := re.match(r'cm1rst_t(\d+)_s\.dat$', p.name))}
    complete = sorted(n for n in found if all((case / f'cm1rst_t{n:06d}_{part}.dat').exists() for part in 'isuvwx'))
    newest = set(complete[-keep:]) if keep > 0 else set()
    removed = [n for n in complete if n > start and n not in newest]
    for n in removed:
        for part in 'isuvwx':
            (case / f'cm1rst_t{n:06d}_{part}.dat').unlink()
    return removed


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
    executable = hashlib.sha256((case / 'cm1.exe').resolve().read_bytes()).hexdigest()[:16]   # the build each segment ran
    try:
        while True:
            done = progress['segments'][-1]['model_s'] if progress['segments'] else progress.get('start_s', 0.0)
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
                                             finished=time.strftime('%Y-%m-%d %H:%M:%S'), executable=executable,
                                             threads=threads))
            progress['case'] = record['case']
            _save_progress(case, progress)
            if cfg.get('keep_restarts'):                                  # a case that keeps only its newest restarts
                prune_restarts(case, cfg['keep_restarts'], round(progress.get('start_s', 0.0) / restart_s))
    finally:
        (case / 'run.lock').unlink(missing_ok=True)
    return progress


def status(name: str) -> str:
    case = RUNS / name
    progress = _progress(case)
    cfg = CASES[name]
    scale = json.loads((case / 'case.json').read_text()).get('time_scale', 1.0)
    start = progress.get('start_s', 0.0)                         # a case run on from another's restart starts there
    done = (progress['segments'][-1]['model_s'] if progress['segments'] else start) / scale   # lunar-equivalent
    total = case_length_s(cfg, scale) / scale
    wall = sum(s['wall_s'] for s in progress['segments'])
    rate = (done * scale - start) / wall if wall else float('nan')
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
        if record.get('sample'):
            print(f"{case}: CM1's {record['sample']} case, {g['nx']} x {g['ny']} x {g['nz']} points, dx {g['dx_m']:.0f} m")
        else:
            print(f"{case}: {g['nx']}{' x ' + str(g['ny']) if 'ny' in g else ''} x {g['nz']} points, dx {g['dx_m']:.1f} m, "
                  f"top {g['ztop_m'] / 1000:.0f} km, water {record['water_share']:.2f} of the ring")
    return 0


if __name__ == '__main__':
    sys.exit(main())
