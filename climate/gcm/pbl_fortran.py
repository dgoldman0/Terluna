"""The nonlocal boundary layer of pbl_scheme.py as Fortran for PlaSim's vertical diffusion (fluxmod.f90, vdiff), behind
the namelist switch nlpbl, and a completion of PlaSim's diagnostic of the convection's vapour (rainmod.f90): the routines'
source and the edits that put them into ExoPlaSim's source, as exoplasim_run.FREE_CONVECTION_EDITS and
massflux_fortran.MASSFLUX_EDITS hold theirs. The runner applies both patches for a run that sets the boundary-layer
options (exoplasim_run.ensure_patched); the vdiff patch goes over the free-convection constant's patch.

The vdiff patch (PBL_PATCHES):
- fluxmod's module variables gain nlpbl (0, the default, leaves vdiff as it was; 1 runs the scheme) and the scheme's
  constants (PARAMETER_NAMES, pbl_scheme.PARAMETERS under Fortran names), each in fluxmod's namelist and broadcast.
- vdiff, after its local diffusivities and with nlpbl 1, keeps the state it diffuses (zt0, zq0) and calls pblnl, which
  for each column takes PlaSim's surface fluxes from surflx (dshfl and devap, downward, over the surface air's density
  ps / (R T_s) and PlaSim's moist heat capacity for the sensible heat, as mkshfl and mkevap form them) and the friction
  velocity (dust3, u*^3 from mkstress), and hands them to pblcol. pblcol is pbl_scheme.nonlocal_column step for step: it
  raises the diffusivities inside the boundary layer to the profile's and adds the countergradient and entrainment
  increments to the state. vdiff then solves as before, its tendencies taken against the state it kept.
- Runs without the switch compute every number as before: the new statements sit behind nlpbl.

The diagnostic patch (DIAGNOSTIC_PATCHES): PlaSim's gridpoint diagnostic 16, the vapour each layer gains by convection
(code 36, a water flux, m/s), is preset to minus the vapour tendency so far by the simplified Betts-Miller scheme
(mksbm) and the mass-flux scheme (mkmassflux), as Kuo's routine presets it, but only Kuo's adds the tendency afterwards.
The patch adds it in both, so the slot holds the scheme's vapour by layer (with PlaSim's shallow scheme's under mksbm).
The slot is a diagnostic only, written only when a run asks for code 36.
"""
from __future__ import annotations

from climate.gcm import pbl_scheme as ps

PBL_MARKER = 'Terluna: the nonlocal boundary layer'
DIAGNOSTIC_MARKER = "Terluna: the convection scheme's vapour by layer"
# (Fortran name, pbl_scheme.PARAMETERS key, what the variable is)
PARAMETER_NAMES = (
    ('pbl_b', 'b_excess', "thermals' temperature excess constant b (Holtslag and Boville 1993)"),
    ('pbl_a', 'a_countergradient', 'countergradient constant a'),
    ('pbl_eps', 'eps_surface', "surface layer's share of the boundary layer"),
    ('pbl_betam', 'beta_m', 'unstable wind-gradient constant (Dyer 1974)'),
    ('pbl_betah', 'beta_h', 'unstable temperature-gradient constant'),
    ('pbl_betas', 'beta_s', 'stable gradient constant'),
    ('pbl_ricrs', 'ricr_stable', 'critical bulk Richardson number, stable or neutral (CAM3)'),
    ('pbl_ricru', 'ricr_unstable', 'critical bulk Richardson number, unstable (Hong et al. 2006)'),
    ('pbl_bu', 'b_ustar', 'u*^2 factor in the bulk number (Vogelezang and Holtslag 1996)'),
    ('pbl_ent', 'a_entrain', 'entrainment flux at h against the surface buoyancy flux (Hong et al. 2006)'),
    ('pbl_entu', 'c_entrain_ustar', "u*^3 factor in the entrainment's velocity scale"),
)


def fortran_literal(value) -> str:
    return repr(float(value))


def fortran_values(prm=None) -> dict:
    """Each Fortran constant's value: pbl_scheme.PARAMETERS, or prm."""
    prm = ps.PARAMETERS if prm is None else prm
    return {name: float(prm[key]) for name, key, _ in PARAMETER_NAMES}


def declarations(values=None) -> str:
    """The module variables, one a line: the switch, then the constants at their defaults."""
    values = fortran_values() if values is None else values
    lines = [f'      integer :: nlpbl = 0 ! {PBL_MARKER}: 1 runs it in vdiff (pblnl)\n']
    for name, _, what in PARAMETER_NAMES:
        lines.append(f'      real :: {name} = {fortran_literal(values[name])} ! {what}\n')
    return ''.join(lines)


def namelist_names() -> list:
    return ['nlpbl'] + [name for name, _, _ in PARAMETER_NAMES]


def namelist_text() -> str:
    """fluxmod's namelist line ending in cfree, continued with this patch's names on one more line."""
    return ('     &                ,zumin,vdiff_lamm,vdiff_b,vdiff_c,vdiff_d,cfree &\n'
            '     &                ,' + ','.join(namelist_names()) + '\n')


def broadcast_text() -> str:
    lines = ['      call mpbci(nlpbl)\n'] + [f'      call mpbcr({name})\n' for name, _, _ in PARAMETER_NAMES]
    return ''.join(lines)


# The routines, appended after vdiff. pblnl: the columns, PlaSim's arrays to pblcol's arguments; pblcol: one column.
PBL_CODE = """
!     ================
!     SUBROUTINE PBLNL
!     ================

      subroutine pblnl(pt,pq,pu,pv,pkm,pkh)
      use fluxmod
!
!     Terluna: the nonlocal boundary layer (Holtslag and Boville 1993, as CAM3 implements it, with the explicit
!     entrainment flux of Hong, Noh and Dudhia 2006), called by vdiff with nlpbl = 1 after its local diffusivities:
!     each column's kinematic surface fluxes, from surflx's dshfl and devap over the surface air's density ps/(R Ts)
!     and the moist heat capacity mkshfl uses, and its friction velocity (dust3 = u*^3), go to pblcol, which raises
!     the diffusivities pkm and pkh (sigma form) inside the boundary layer and adds its explicit fluxes to pt and pq.
!
      real pt(NHOR,NLEV),pq(NHOR,NLEV),pu(NHOR,NLEV),pv(NHOR,NLEV)
      real pkm(NHOR,NLEV),pkh(NHOR,NLEV)
      real zt(NLEV),zq(NLEV),zu(NLEV),zv(NLEV),zkm(NLEM),zkh(NLEM),zinfo(10)
      real zrho,zwth0,zwq0,zust
      integer jhor
!
      do jhor=1,NHOR
       zrho=dp(jhor)/(gascon*dt(jhor,NLEP))
       zwth0=-dshfl(jhor)/(zrho*acpd*AMAX1(1.,1.+ADV*dq(jhor,NLEV)))
       zwq0=-1000.*devap(jhor)/zrho
       zust=dust3(jhor)**(1./3.)
       zt(:)=pt(jhor,:)
       zq(:)=pq(jhor,:)
       zu(:)=pu(jhor,:)
       zv(:)=pv(jhor,:)
       zkm(:)=pkm(jhor,1:NLEM)
       zkh(:)=pkh(jhor,1:NLEM)
       call pblcol(NLEV,zt,zq,zu,zv,sigma,sigmah,dsigma,ga,gascon,akap,rdbrv,deltsec2,zwth0,zwq0,zust,zkm,zkh,zinfo)
       pt(jhor,:)=zt(:)
       pq(jhor,:)=zq(:)
       pkm(jhor,1:NLEM)=zkm(:)
       pkh(jhor,1:NLEM)=zkh(:)
      enddo
      return
      end subroutine pblnl

!     =================
!     SUBROUTINE PBLCOL
!     =================

      subroutine pblcol(kn,pt,pq,pu,pv,psig,psigh,pdsig,pg,pgas,pkap,prdbrv,pdt,pwth0,pwq0,pust,pkm,pkh,pinfo)
      use fluxmod
!
!     Terluna: the nonlocal boundary layer on one column (top first), pbl_scheme.nonlocal_column step for step. The
!     depth h comes from the bulk Richardson number between the lowest level and the levels above (Vogelezang and
!     Holtslag's form, critical pbl_ricru in unstable air with the thermals' excess pbl_b (w'thv')0/wm, pbl_ricrs in
!     stable air); inside h the diffusivity is k w z (1-z/h)**2 with w = wm = (u*^3 + betam eps k B h)^(1/3) for
!     momentum and wm/Pr for heat and vapour, Pr from similarity at eps h, or u*/phi_h(z/L) in stable air; heat and
!     vapour gain the countergradient flux K gamma, gamma = a w* (w'c')0/(wm^2 h), and an entrainment flux at h,
!     -pbl_ent (thv/g) (w*^3 + pbl_entu u*^3)/h for thv, carried by each across the jump between the levels bracketing
!     h and falling as (z/h)^3 below. pkm and pkh (sigma form) become the larger of their own and the profile's inside
!     h; the explicit fluxes' increments join pt and pq, the vapour's only where none leaves it negative.
!     pinfo: h, w*, wm, Pr, excess, entrainment velocity, unstable (1), vapour kept (1), inversion interface, B.
!
      integer kn
      real pt(kn),pq(kn),pu(kn),pv(kn),psig(kn),psigh(kn),pdsig(kn)
      real pg,pgas,pkap,prdbrv,pdt,pwth0,pwq0,pust
      real pkm(kn-1),pkh(kn-1),pinfo(10)
!
      real zth(kn),zthv(kn),ztv(kn),zzf(kn),zzh(kn)
      real ztoz(kn-1),zfac(kn-1),zkmp(kn-1),zkhp(kn-1)
      real zfcgt(kn-1),zfcgq(kn-1),zfent(kn-1),zfenq(kn-1),zgt(kn-1),zgq(kn-1)
      real zdtinc(kn),zdqinc(kn)
      real zepsv,zthvn,zwthv0,zbuoy,zust3,zeps,zh0,zwm0,zexc,zh,zwst3,zwst,zphim,zphih,zwm,zpr
      real zgamt,zgamq,zz,zshape,zbneg,zwt,zrthl,zkused,zwent3,zwent,zfthv,zdthv,zwe,zft,zfq,zr3
      integer jk,jb,ja
      logical lunst,lin(kn-1),lkept
!
      zepsv=1./prdbrv-1.
      zth(:)=pt(:)/psig(:)**pkap
      zthv(:)=zth(:)*(1.+zepsv*pq(:))
      ztv(:)=pt(:)*(1.+zepsv*pq(:))
      zzh(kn)=0.
      do jk=kn,2,-1
       zzh(jk-1)=zzh(jk)+pgas*ztv(jk)/pg*LOG(psigh(jk)/psigh(jk-1))
      enddo
      zzf(:)=zzh(:)+pgas*ztv(:)/pg*LOG(psigh(:)/psig(:))
      do jk=1,kn-1
       zrthl=(pdsig(jk+1)+pdsig(jk))/(pt(jk)*pdsig(jk+1)+pt(jk+1)*pdsig(jk))
       ztoz(jk)=pg*psigh(jk)*zrthl/pgas
      enddo
      zthvn=zthv(kn)
      zwthv0=pwth0*(1.+zepsv*pq(kn))+zepsv*zth(kn)*pwq0
      zbuoy=pg/zthvn*zwthv0
      lunst=(zwthv0 > 0.)
      zust3=pust**3
      zeps=pbl_eps
!
!     the depth, the velocity scales and the countergradient terms
!
      if(lunst) then
       zh0=pbldepth(zthvn,pbl_ricru)
       zwm0=(zust3+pbl_betam*zeps*vonkarman*zbuoy*zh0)**(1./3.)
       zexc=pbl_b*zwthv0/zwm0
       zh=pbldepth(zthvn+zexc,pbl_ricru)
       zwst3=zbuoy*zh
       zwst=zwst3**(1./3.)
       zphim=(zust3+pbl_betam*zeps*vonkarman*zbuoy*zh)**(1./3.)
       zphih=SQRT(zust3+pbl_betah*zeps*vonkarman*zbuoy*zh)
       zwm=zphim
       zpr=SQRT(pust)*zphim/zphih+zeps*vonkarman*pbl_a*zwst/zwm
       zgamt=pbl_a*zwst*pwth0/(zwm*zwm*zh)
       zgamq=pbl_a*zwst*pwq0/(zwm*zwm*zh)
      else
       zh=pbldepth(zthvn,pbl_ricrs)
       zexc=0.
       zwst3=0.
       zwst=0.
       zwm=0.
       zpr=1.
       zgamt=0.
       zgamq=0.
      endif
!
!     the profile's diffusivities inside h and the larger of those and the local ones (sigma form)
!
      zkmp(:)=0.
      zkhp(:)=0.
      do jk=1,kn-1
       lin(jk)=(zzh(jk) < zh)
       zfac(jk)=pdt*ztoz(jk)*ztoz(jk)/(psig(jk+1)-psig(jk))
       if(lin(jk)) then
        zz=zzh(jk)
        zshape=zz*(1.-zz/zh)**2
        if(lunst) then
         zkmp(jk)=vonkarman*zwm*zshape
         zkhp(jk)=zkmp(jk)/zpr
        else
         zbneg=-zbuoy
         if(pust <= 0.) then
          zwt=0.
         elseif(zz*vonkarman*zbneg <= zust3) then
          zwt=pust*zust3/(zust3+pbl_betas*zz*vonkarman*zbneg)
         else
          zwt=pust*zust3/(pbl_betas*zust3+zz*vonkarman*zbneg)
         endif
         zkmp(jk)=vonkarman*zwt*zshape
         zkhp(jk)=zkmp(jk)
        endif
        pkm(jk)=MAX(pkm(jk),zkmp(jk)*zfac(jk))
        pkh(jk)=MAX(pkh(jk),zkhp(jk)*zfac(jk))
       endif
      enddo
!
!     the countergradient fluxes with the diffusivity used, and the entrainment flux (kinematic, upward)
!
      zfcgt(:)=0.
      zfcgq(:)=0.
      zfent(:)=0.
      zfenq(:)=0.
      if(lunst) then
       do jk=1,kn-1
        if(lin(jk)) then
         zkused=pkh(jk)/zfac(jk)
         zfcgt(jk)=zkused*zgamt
         zfcgq(jk)=zkused*zgamq
        endif
       enddo
      endif
      jb=kn
      do jk=1,kn
       if(zzf(jk) <= zh) then
        jb=jk
        exit
       endif
      enddo
      zwe=0.
      ja=0
      if(lunst .and. jb > 1) then
       ja=jb-1
       zwent3=zwst3+pbl_entu*zust3
       zwent=zwent3**(1./3.)
       zfthv=-pbl_ent*zthvn/pg*zwent3/zh
       zdthv=zthv(ja)-zthv(jb)
       if(zdthv <= 0.) then
        zwe=zwent
       else
        zwe=MIN(-zfthv/zdthv,zwent)
       endif
       zft=-zwe*(zth(ja)-zth(jb))
       zfq=-zwe*(pq(ja)-pq(jb))
       zfent(ja)=zft
       zfenq(ja)=zfq
       do jk=ja+1,kn-1
        zr3=(zzh(jk)/zh)**3
        zfent(jk)=zft*zr3
        zfenq(jk)=zfq*zr3
       enddo
      endif
!
!     the explicit increments: each interface's flux out of the layer above into the layer below it
!
      zgt(:)=psigh(1:kn-1)**pkap*pdt*ztoz(:)*(zfcgt(:)+zfent(:))
      zgq(:)=pdt*ztoz(:)*(zfcgq(:)+zfenq(:))
      zdtinc(1)=zgt(1)/pdsig(1)
      zdqinc(1)=zgq(1)/pdsig(1)
      do jk=2,kn-1
       zdtinc(jk)=(zgt(jk)-zgt(jk-1))/pdsig(jk)
       zdqinc(jk)=(zgq(jk)-zgq(jk-1))/pdsig(jk)
      enddo
      zdtinc(kn)=(0.-zgt(kn-1))/pdsig(kn)
      zdqinc(kn)=(0.-zgq(kn-1))/pdsig(kn)
      lkept=ANY(pq(:)+zdqinc(:) < 0.)
      if(lkept) zdqinc(:)=0.
      pt(:)=pt(:)+zdtinc(:)
      pq(:)=pq(:)+zdqinc(:)
      pinfo(:)=0.
      pinfo(1)=zh
      pinfo(2)=zwst
      pinfo(3)=zwm
      pinfo(4)=zpr
      pinfo(5)=zexc
      pinfo(6)=zwe
      if(lunst) pinfo(7)=1.
      if(lkept) pinfo(8)=1.
      pinfo(9)=real(ja)
      pinfo(10)=zbuoy
      return
!
      contains
!
      real function pbldepth(pthref,pricr)
!     the bulk Richardson number from the lowest level up; h linear in the number between levels (CAM's pblintd)
      real pthref,pricr
      real zrib,zri,zvvk
      integer jl
      zrib=0.
      do jl=kn-1,1,-1
       zvvk=MAX((pu(jl)-pu(kn))**2+(pv(jl)-pv(kn))**2+pbl_bu*pust*pust,1.e-36)
       zri=pg*(zthv(jl)-pthref)*(zzf(jl)-zzf(kn))/(zthvn*zvvk)
       if(zri >= pricr) then
        if(zri-zrib <= 0.) then
         pbldepth=zzf(jl+1)
        else
         pbldepth=zzf(jl+1)+(pricr-zrib)/(zri-zrib)*(zzf(jl)-zzf(jl+1))
        endif
        return
       endif
       zrib=zri
      enddo
      pbldepth=zzf(1)
      end function pbldepth
!
      end subroutine pblcol
"""

VDIFF_EDITS = [
    ("      real zke(NHOR,NLEV),zken(NHOR,NLEV)\n",
     "      real zke(NHOR,NLEV),zken(NHOR,NLEV)\n"
     f"      real zt0(NHOR,NLEV),zq0(NHOR,NLEV) ! {PBL_MARKER}: the state before its fluxes\n"),
    ("        zkdiffm(jhor,jlev)=zfac*zmixm2*zrifm*zdvds*zrthl*zrthl*zrthl\n       enddo\n      enddo\n",
     "        zkdiffm(jhor,jlev)=zfac*zmixm2*zrifm*zdvds*zrthl*zrthl*zrthl\n       enddo\n      enddo\n"
     "!\n"
     f"!     {PBL_MARKER} (pblnl), with nlpbl = 1: its diffusivities inside h and its explicit fluxes\n"
     "!\n"
     "      if(nlpbl == 1) then\n"
     "       zt0(:,:)=zt(:,:)\n"
     "       zq0(:,:)=zq(:,:)\n"
     "       call pblnl(zt,zq,zu,zv,zkdiffm,zkdiffh)\n"
     "      endif\n"),
    ("      zdqdt(:,1:NLEV)=(zqn(:,1:NLEV)-zq(:,1:NLEV))/deltsec2\n",
     "      zdqdt(:,1:NLEV)=(zqn(:,1:NLEV)-zq(:,1:NLEV))/deltsec2\n"
     "      if(nlpbl == 1) zdqdt(:,1:NLEV)=(zqn(:,1:NLEV)-zq0(:,1:NLEV))/deltsec2\n"),
    ("      zdtdt(:,1:NLEV)=(ztn(:,1:NLEV)-zt(:,1:NLEV))/deltsec2\n",
     "      zdtdt(:,1:NLEV)=(ztn(:,1:NLEV)-zt(:,1:NLEV))/deltsec2\n"
     "      if(nlpbl == 1) zdtdt(:,1:NLEV)=(ztn(:,1:NLEV)-zt0(:,1:NLEV))/deltsec2\n"),
    ("      end subroutine vdiff\n", "      end subroutine vdiff\n" + PBL_CODE),
]


def module_edits() -> list:
    return [
        ("      real :: cfree      = 0.0016 ! Terluna: the free-convection constant over water (m s-1 K-1/3; Miller et al. "
         "1992)\n",
         "      real :: cfree      = 0.0016 ! Terluna: the free-convection constant over water (m s-1 K-1/3; Miller et al. "
         "1992)\n" + declarations()),
        ("     &                ,zumin,vdiff_lamm,vdiff_b,vdiff_c,vdiff_d,cfree\n", namelist_text()),
        ("      call mpbcr(cfree)\n", "      call mpbcr(cfree)\n" + broadcast_text()),
    ]


PBL_EDITS = module_edits() + VDIFF_EDITS
PBL_PATCHES = [('fluxmod.f90', PBL_MARKER, PBL_EDITS)]

DIAGNOSTIC_EDITS = [
    ("        dgp3d(:,jlev,8)=zdtcv(:,jlev)\n",
     "        dgp3d(:,jlev,8)=zdtcv(:,jlev)\n"
     f"!       {DIAGNOSTIC_MARKER}: the tendency added to the slot's preset\n"
     "        dgp3d(:,jlev,16)=dgp3d(:,jlev,16)+dqdt(:,jlev)*dp(:)*dsigma(jlev)/ga/1000.\n"),
    ("        dgp3d(:,jlev,8)=zdtc(:,jlev)+zdtdts(:,jlev)\n",
     "        dgp3d(:,jlev,8)=zdtc(:,jlev)+zdtdts(:,jlev)\n"
     f"!       {DIAGNOSTIC_MARKER}: the tendency added to the slot's preset\n"
     "        dgp3d(:,jlev,16)=dgp3d(:,jlev,16)+dqdt(:,jlev)*dp(:)*dsigma(jlev)/ga/1000.\n"),
]
DIAGNOSTIC_PATCHES = [('rainmod.f90', DIAGNOSTIC_MARKER, DIAGNOSTIC_EDITS)]


def surface_from_plasim(dshfl: float, devap: float, dust3: float, ts: float, ps_pa: float, q_low: float, c: dict):
    """pblnl's arithmetic: the kinematic surface fluxes of potential temperature (K m/s) and vapour (kg/kg m/s), upward,
    and the friction velocity, from surflx's dshfl (W/m2, downward), devap (m/s of water, downward) and dust3 (u*^3),
    the surface temperature, the surface pressure and the lowest level's vapour before this step."""
    rho = ps_pa / (c['gascon'] * ts)
    return dict(wth0=-dshfl / (rho * c['acpd'] * max(1.0, 1.0 + c['ADV'] * q_low)), wq0=-1000.0 * devap / rho,
                ustar=dust3 ** (1.0 / 3.0))
