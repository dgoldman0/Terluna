"""The mass-flux convection scheme of massflux.py as Fortran for PlaSim (rainmod.f90), behind PlaSim's switch nprc = 3,
with Kuo's scheme (1) and the simplified Betts-Miller scheme (2) kept as they are: the routine's source and the edits
that put it into ExoPlaSim's source, as exoplasim_run.SBM_CODE and SBM_EDITS hold that scheme's. The patch needs the
Betts-Miller patch under it, which every run's source has. It is not among the patches exoplasim_run.ensure_patched
applies yet.

The source (MASSFLUX_CODE) has three routines:
- mkmassflux, called in rainstep in Kuo's place: each column's state after the other processes' tendencies so far, its
  surface fluxes (PlaSim's dshfl and dlhfl, upward), its boundary-layer forcing (those tendencies with the dynamics'
  heating and moistening), the large-scale vertical velocity dw and the wind speed go to mfcolumn. Its tendencies join
  dtdt and dqdt; its rain and snow by layer, net of its downdraught's evaporation, go to mkrain as dprcl and dprscl;
  the layers its updraught crosses become mkclouds' convective cloud, as Kuo's rain levels do. PlaSim's diffusive
  shallow scheme is not called, since the scheme convects shallow itself.
- mfcolumn, the scheme on one column: massflux.column with melt 'plasim', step for step, its arrays top first.
- mkdttgp, the dynamics' temperature tendency at the grid points: the spectral tendency stt transformed as mkdqtgp
  transforms sqt, with PlaSim's own scaling of it (ct ww), for the closure's boundary-layer term.

The edits (MASSFLUX_EDITS) add the scheme's parameters to rainmod's module variables, its namelist and its broadcasts,
the call in rainstep and the routines, and one rule in mkrain: under nprc 3, the convective water mkrain evaporates in
a layer is at most the least running total of the water handed on from that layer down, less what it has evaporated
above, so the downdraught's evaporation, already taken from the handed water, never leaves a running total negative
(massflux.evaporate is its copy). Under nprc 1 and 2 that cap is 1e30, so mkrain's arithmetic is as before.

The parameters (PARAMETER_NAMES) are massflux.PARAMETERS under Fortran names ending in _mf; fortran_values gives
their values for the source's defaults and a run's namelist. Six options exist only in the Python copy, for its
offline variants, and must keep their defaults for the Fortran (PYTHON_ONLY), the entrainment that follows the
environment (massflux_entrainment.py) among them.

A second patch (START_PATCHES, 26 October) goes over this one for runs that set the surface parcel's starting speed from
the boundary layer (massflux.surface_start, the runner's --massflux-start): mfstart and its two parameters
(START_PARAMETER_NAMES), off by default.
"""
from __future__ import annotations

from climate.gcm import massflux as mf

MASSFLUX_MARKER = 'Terluna: the mass-flux scheme'
# (Fortran name, massflux.PARAMETERS key, Fortran type, what the variable is)
PARAMETER_NAMES = (
    ('epsbg_mf', 'eps_bg', 'real', "updraught entrainment aloft (/m)"),
    ('epsnb_mf', 'eps_nb', 'real', "its added part at the cloud base (/m)"),
    ('lnb_mf', 'l_nb', 'real', "the height over which that part decays (m)"),
    ('shfac_mf', 'shallow_factor', 'real', "shallow clouds' entrainment against deep ones'"),
    ('gvm_mf', 'gamma_vm', 'real', "kinetic energy: virtual mass coefficient"),
    ('fturb_mf', 'f_turb', 'real', "kinetic energy: turbulence factor"),
    ('betacd_mf', 'beta_cd', 'real', "kinetic energy: drag"),
    ('w_mf', 'w_base_deep', 'real', "deep and mid-level test parcels' starting speed (m/s)"),
    ('dt_mf', 'dt_deep', 'real', "test parcels' temperature excess (K)"),
    ('dq_mf', 'dq_deep', 'real', "test parcels' vapour excess (kg/kg)"),
    ('deepdp_mf', 'deep_dp', 'real', "a cloud thicker than this convects deep (Pa)"),
    ('departure_mf', 'departure_dp', 'real', "deep test parcels depart up to this far above the ground (Pa)"),
    ('subf_mf', 'sub_factor', 'real', "test parcels' sub-cloud mixing against the deep base rate"),
    ('zsrc_mf', 'source_height_earth', 'real', "the surface parcel's source height at Earth's gravity (m)"),
    ('ustar_mf', 'ustar', 'real', "friction velocity of the surface excess and similarity (m/s)"),
    ('c0_mf', 'c0', 'real', "rain formation rate (/s)"),
    ('lcrit_mf', 'l_crit', 'real', "its critical condensate (kg/kg)"),
    ('lmin_mf', 'l_min', 'real', "no rain below this condensate (kg/kg)"),
    ('wmaxp_mf', 'w_max_precip', 'real', "updraught speed's cap in the rain rate (m/s)"),
    ('tbf_mf', 't_bf', 'real', "Bergeron-Findeisen enhancement below this temperature (K)"),
    ('tice_mf', 't_ice', 'real', "condensate all frozen below this temperature (K)"),
    ('etad_mf', 'eta_down', 'real', "downdraught's mass flux at its start against the base flux"),
    ('epsd_mf', 'eps_down', 'real', "downdraught's mixing (/m)"),
    ('botdp_mf', 'down_bottom_dp', 'real', "downdraught spreads over the lowest this (Pa)"),
    ('tau_mf', 'tau_factor', 'real', "deep adjustment time in turnover times"),
    ('taumin_mf', 'tau_min', 'real', "its floor (s)"),
    ('tstar_mf', 't_star', 'real', "PCAPEbl's T* (K; 0: g H / cp, H the cloud's depth)"),
    ('ublmin_mf', 'u_bl_min', 'real', "sub-cloud wind's floor over water (m/s)"),
    ('midrh_mf', 'midlevel_rh', 'real', "mid-level departures at least this moist"),
    ('midtop_mf', 'midlevel_top', 'real', "mid-level departures below this sigma"),
    ('cfl_mf', 'cfl', 'real', "mass flux's limit in layer masses a step"),
    ('nelev_mf', 'elevated_ascent', 'int', "1: elevated deep departures only under large-scale ascent"),
    ('nsub_mf', 'substeps', 'int', "the plume's steps across each half layer"),
)
# Options of the Python copy that the Fortran does not carry, at the values the Fortran's arithmetic takes.
PYTHON_ONLY = dict(rh_dependence=False, delta_form='ratio', source_dp=0.0, excess_scale='source',
                   stability='free_convection', entrainment_form='fixed')


def fortran_values(prm=None) -> dict:
    """The Fortran parameters' values for the Python copy's parameters prm (massflux.PARAMETERS by default): t_star
    'cloud' as 0, the elevated-ascent switch as 1 or 0. Refuses parameters the Fortran does not carry."""
    prm = dict(mf.PARAMETERS) if prm is None else prm
    for key, value in PYTHON_ONLY.items():
        if prm[key] != value:
            raise ValueError(f'the Fortran routine carries {key} = {value!r} only, not {prm[key]!r}')
    out = {}
    for name, key, kind, _ in PARAMETER_NAMES:
        value = prm[key]
        if key == 't_star':
            value = 0.0 if value == 'cloud' else float(value)
        if kind == 'int':
            out[name] = int(value)
        else:
            out[name] = float(value)
    return out


def fortran_literal(value) -> str:
    """A value as Fortran reads it back to the same double (with -fdefault-real-8)."""
    if isinstance(value, int):
        return str(value)
    text = repr(float(value))
    return text if ('.' in text or 'e' in text) else text + '.'


def declarations(values=None) -> str:
    values = fortran_values() if values is None else values
    lines = []
    for name, _, kind, what in PARAMETER_NAMES:
        decl = f"      {'integer' if kind == 'int' else 'real'} :: {name} = {fortran_literal(values[name])}"
        lines.append(f'{decl:42s} ! {MASSFLUX_MARKER} (nprc = 3): {what}\n' if not lines else f'{decl:42s} ! {what}\n')
    return ''.join(lines)


def _wrapped(prefix: str, items, width: int = 120) -> str:
    """Fortran statements prefix + comma-separated items, as many statements as the width needs."""
    out, line = [], []
    for item in items:
        if line and len(prefix) + len(','.join(line + [item])) > width:
            out.append(prefix + ','.join(line))
            line = []
        line.append(item)
    if line:
        out.append(prefix + ','.join(line))
    return ''.join(f'{text}\n' for text in out)


def namelist_text() -> str:
    return _wrapped('      namelist/rainmod_nl/', [name for name, _, _, _ in PARAMETER_NAMES])


def broadcast_text() -> str:
    return ''.join(f"      call {'mpbci' if kind == 'int' else 'mpbcr'}({name})\n"
                   for name, _, kind, _ in PARAMETER_NAMES)


WRAPPER_CODE = """
!     ======================
!     SUBROUTINE MKMASSFLUX
!     ======================

      subroutine mkmassflux
      use rainmod
      implicit none
!
!     Terluna: the mass-flux convection scheme (climate/gcm/massflux.py, its record massflux.md), in place of Kuo's
!     with nprc = 3. Each column's state after the other processes' tendencies so far goes to mfcolumn with the
!     surface's fluxes (upward), the boundary layer's forcing (those tendencies with the dynamics' heating and
!     moistening), the large-scale vertical velocity and the wind speed. Its tendencies join dtdt and dqdt; its rain
!     and snow by layer, net of its downdraught's evaporation, go to mkrain as dprcl and dprscl, which evaporates
!     them on the way down under its floor and turns them to rain or snow by temperature; the layers its updraught
!     crosses are mkclouds' convective cloud. PlaSim's diffusive shallow scheme is not called: the scheme convects
!     shallow itself.
!
      integer :: jhor, jlev, iland, iinfo(4)
      real :: zt(NLEV), zq(NLEV), zpf(NLEV), zph(0:NLEV), zbt(NLEV), zbq(NLEV), zom(NLEV), zwsp(NLEV), ztph(NLEV)
      real :: zdtc(NLEV), zdqc(NLEV), zrain(NLEV), zsnow(NLEV), zmu(0:NLEV), zmd(0:NLEV), zinfo(10)
      real :: zdtt(NHOR,NLEV), zdtcv(NHOR,NLEV)
!
      icctot(:)=0
      icclev(:,:)=0
      dprcl(:,:)=0.
      dprscl(:,:)=0.
      zdtcv(:,:)=0.
      if(ndiaggp==1) then
       do jlev=1,NLEV
        dgp3d(:,jlev,16)=-dqdt(:,jlev)*dp(:)*dsigma(jlev)/ga/1000.
        dgp3d(:,jlev,18)=dqt(:,jlev)*dp(:)*dsigma(jlev)/ga/1000.
       enddo
      end if
      call mkdttgp(zdtt)
      do jhor=1,NHOR
       zt(:)=dt(jhor,1:NLEV)+dtdt(jhor,1:NLEV)*deltsec2
       zq(:)=dq(jhor,1:NLEV)+dqdt(jhor,1:NLEV)*deltsec2
       zpf(:)=dp(jhor)*sigma(:)
       zph(0)=0.
       zph(1:NLEV)=dp(jhor)*sigmah(1:NLEV)
       zbt(:)=dtdt(jhor,1:NLEV)+zdtt(jhor,:)
       zbq(:)=dqdt(jhor,1:NLEV)+dqt(jhor,1:NLEV)
       zom(:)=dw(jhor,1:NLEV)
       zwsp(:)=sqrt(du(jhor,1:NLEV)**2+dv(jhor,1:NLEV)**2)
       ztph(:)=dt(jhor,1:NLEV)
       iland=0
       if(dls(jhor) > 0.5) iland=1
       call mfcolumn(NLEV,zt,zq,zpf,zph,ga,deltsec2,-dshfl(jhor),-dlhfl(jhor),zbt,zbq,zom,zwsp,iland,ztph,         &
     &               zdtc,zdqc,zrain,zsnow,zmu,zmd,iinfo,zinfo)
       if(iinfo(1) == 0) cycle
       dtdt(jhor,1:NLEV)=dtdt(jhor,1:NLEV)+zdtc(:)
       dqdt(jhor,1:NLEV)=dqdt(jhor,1:NLEV)+zdqc(:)
       zdtcv(jhor,:)=zdtc(:)
       dprcl(jhor,:)=(zrain(:)+zsnow(:))*ga/dp(jhor)
       dprscl(jhor,:)=zsnow(:)*ga/dp(jhor)
!      the convective cloud: the layers from the one holding the updraught's top to the one above its base
       do jlev=iinfo(3),iinfo(2)
        icclev(jhor,jlev)=1
        icctot(jhor)=icctot(jhor)+1
       enddo
       if(sigma(iinfo(3)) < 0.4) icclev(jhor,iinfo(3))=2
      enddo
      if(ndiaggp==1) then
       do jlev=1,NLEV
        dgp3d(:,jlev,8)=zdtcv(:,jlev)
       enddo
      end if
      return
      end subroutine mkmassflux

"""
DTT_CODE = """
!     ===================
!     SUBROUTINE MKDTTGP
!     ===================

      subroutine mkdttgp(pdtt)
      use rainmod
      implicit none
!
!     Terluna: the dynamics' temperature tendency (K/s) at the grid points, for the mass-flux scheme's boundary-layer
!     forcing: the spectral tendency stt as the step's spectral part leaves it, transformed as mkdqtgp transforms
!     sqt, with the scaling PlaSim's energy diagnostics give it
!
      real, intent(out) :: pdtt(NHOR,NLEV)
      real :: zstt(NESP,NLEV)
      call mpgallsp(zstt,stt,NLEV)
      zstt(:,:)=zstt(:,:)*ct*ww
      call sp2fl(zstt,pdtt,NLEV)
      call fc2gp(pdtt,NLON,NLPP*NLEV)
      return
      end subroutine mkdttgp
"""
COLUMN_CODE = """
!     =====================
!     SUBROUTINE MFCOLUMN
!     =====================

      subroutine mfcolumn(kn,pt,pq,ppf,pph,pg,pdt,pshf,plhf,pbt,pbq,pom,pwsp,kland,ptph,                       &
     &                    pdtdt,pdqdt,prain,psnow,pmu,pmd,kinfo,pinfo)
      use rainmod
      implicit none
!
!     Terluna: the mass-flux scheme on one column (climate/gcm/massflux.py, its column with PlaSim's phase rule),
!     levels top first: full levels 1..kn, half levels 0 (the top) .. kn (the ground), layer k between half levels
!     k-1 and k. In: temperature (K), specific humidity, full and half level pressures (Pa), gravity, the step over
!     which the tendencies act (s), the surface's sensible and latent heat fluxes (W/m2, upward), the tendencies of
!     temperature and vapour from processes other than convection (K/s, kg/kg/s), the large-scale vertical velocity
!     (Pa/s), the wind speed (m/s), land (1) or water (0), and the temperature mkrain takes the falling water's phase
!     from (K). Out: the tendencies (K/s, kg/kg/s), the rain and snow handed to mkrain by layer, net of the
!     downdraught's evaporation (kg/m2/s), the updraught's and downdraught's mass fluxes on half levels (kg/m2/s), and
!     kinfo (type: 0 none, 1 deep, 2 shallow, 3 mid-level; base and top half levels; 1 if rooted at the ground) and
!     pinfo (base mass flux, base and top heights, base and top pressures, PCAPE, PCAPEbl, tau, 1 if limited by the
!     layers' mass, the updraught's mean speed).
!
      integer, intent(in) :: kn, kland
      real, intent(in) :: pt(kn), pq(kn), ppf(kn), pph(0:kn), pg, pdt, pshf, plhf
      real, intent(in) :: pbt(kn), pbq(kn), pom(kn), pwsp(kn), ptph(kn)
      real, intent(out) :: pdtdt(kn), pdqdt(kn), prain(kn), psnow(kn), pmu(0:kn), pmd(0:kn)
      integer, intent(out) :: kinfo(4)
      real, intent(out) :: pinfo(10)
!
!     the column as the scheme reads it
      real :: zzf(kn), zzh(0:kn), zdp(kn), zmass(kn), zs(kn), ztv(kn), zqs(kn), zrh(kn), zth(0:kn), zqh(0:kn)
!     a path (full and half levels alternately from a departure point up) and a plume's points along it
      integer :: npath, ipk(2*kn+2), ipi(2*kn+2)
      real :: zpz(2*kn+2), zpp(2*kn+2), zpte(2*kn+2), zpqe(2*kn+2)
      integer :: npt, iptk(4*kn+8), iptl(4*kn+8), iptx(4*kn+8)
      real :: zptz(4*kn+8), zptp(4*kn+8), zptt(4*kn+8), zptqv(4*kn+8), zptl(4*kn+8), zptsl(4*kn+8), zptqt(4*kn+8)
      real :: zptm(4*kn+8), zptkin(4*kn+8), zptb(4*kn+8), zptrain(4*kn+8), zptsnow(4*kn+8)
!     the updraught, downdraught and their fluxes
      real :: zmu(0:kn), zslu(0:kn), zqtu(0:kn), zrainu(kn), zsnowu(kn)
      real :: zmd(0:kn), zsd(0:kn), zqd(0:kn), zevd(kn)
      real :: zfs(0:kn), zfq(0:kn), zdsdt(kn), zdqdt(kn), znr(kn), zns(kn), zmuu(0:kn)
      integer :: itype, ib, it, ibase, k, ksh, kdep, i, jj
      logical :: lok, lshok, lrooted, lbuoy(kn)
      real :: zcp, zlv, zls, zeps1, zg, zps, zdts, zdqs, zws, zdtl, zdql, zt0, zq0
      real :: zdepth, ztopz, zbsl, zbqt, zbkin, zzcond, zshdep
      real :: zwmean, zpc, zpbl, ztau, ztauc, ztaubl, zrate, zmb, zlimit, zmnet, zsupply, zhup, zhenv, zexc
      real :: zx, zub, zsum
!
      zg = pg
      zcp = acpd
      zlv = ALV
      zls = ALS
      zeps1 = 1./rdbrv-1.
      zps = pph(kn)
      pdtdt(:) = 0.
      pdqdt(:) = 0.
      prain(:) = 0.
      psnow(:) = 0.
      pmu(:) = 0.
      pmd(:) = 0.
      kinfo(:) = 0
      pinfo(:) = 0.
!
!     heights from the virtual temperatures, the full level at its pressure within the layer
      zzh(kn) = 0.
      do k = kn, 1, -1
       ztv(k) = pt(k)*(1.+zeps1*pq(k))
       zx = gascon*ztv(k)/zg
       zzf(k) = zzh(k)+zx*log(pph(k)/ppf(k))
       if (pph(k-1) > 0.) then
        zzh(k-1) = zzh(k)+zx*log(pph(k)/pph(k-1))
       else
        zzh(k-1) = zzh(k)+zx*log(pph(k)/(ppf(k)*0.5))
       endif
      enddo
      do k = 1, kn
       zdp(k) = pph(k)-pph(k-1)
       zmass(k) = zdp(k)/zg
       zs(k) = zcp*pt(k)+zg*zzf(k)
       zqs(k) = mfqsat(pt(k),ppf(k))
       zrh(k) = pq(k)/zqs(k)
      enddo
!     the surroundings at half levels, linear in log pressure between the full levels on each side
      do i = 1, kn-1
       zx = log(pph(i)/ppf(i))/log(ppf(i+1)/ppf(i))
       zth(i) = (1.-zx)*pt(i)+zx*pt(i+1)
       zqh(i) = (1.-zx)*pq(i)+zx*pq(i+1)
      enddo
      zth(kn) = pt(kn)
      zqh(kn) = pq(kn)
      zth(0) = pt(1)
      zqh(0) = pq(1)
!
!     the test parcels and the type (classify)
      call mfexcess(zdts,zdqs,zws)
      call mfsurface(zdtl,zdql)
      zt0 = pt(kn)+zdtl+max(dt_mf,zdts)
      zq0 = pq(kn)+zdql+max(dq_mf,zdqs)
      call mftrial(kn,zt0,zq0,zws,.true.,.false.,lshok,zshdep)
      itype = 0
      lrooted = .false.
      kdep = 0
      do k = kn, 1, -1
       if (zps-ppf(k) > departure_mf) exit
       if (k < kn .and. nelev_mf == 1 .and. pom(k) >= 0.) cycle
       if (k == kn) then
        zx = max(w_mf,zws)
        call mftrial(k,pt(k)+zdtl+max(dt_mf,zdts),pq(k)+zdql+max(dq_mf,zdqs),zx,.false.,.true.,lok,zdepth)
       else
        call mftrial(k,pt(k)+dt_mf,pq(k)+dq_mf,w_mf,.false.,.true.,lok,zdepth)
       endif
       if (lok .and. zdepth > deepdp_mf) then
        itype = 1
        kdep = k
        lrooted = (k == kn)
        exit
       endif
      enddo
      if (itype == 1) then
       call mfkeep(kdep)
      else
       if (lshok .and. zshdep <= deepdp_mf) then
        itype = 2
        call mftrial(kn,zt0,zq0,zws,.true.,.false.,lok,zdepth)
        call mfkeep(kn)
        lrooted = .true.
       else
        do k = kn, 1, -1
         if (ppf(k) < midtop_mf*zps) exit
         if (pom(k) < 0. .and. zrh(k) >= midrh_mf) then
          call mftrial(k,pt(k)+dt_mf,pq(k)+dq_mf,w_mf,.false.,.true.,lok,zdepth)
          if (lok) then
           itype = 3
           call mfkeep(k)
           exit
          endif
         endif
        enddo
       endif
      endif
      kinfo(1) = itype
      if (itype == 0) return
      if (lrooted) kinfo(4) = 1
!
!     the updraught from the test parcel's cloud base, at unit mass flux there
      call mfupdraught(itype == 2)
      kinfo(2) = ib
      kinfo(3) = it
!     the downdraught
      call mfdowndraught(itype == 2)
!     the unit tendencies, PCAPE and the closures
      call mftend(1.0)
      zmuu(:) = zmu(:)
      do i = 1, kn-1
       if (i > ib) zmuu(i) = zmu(ib)*(zps-pph(i))/(zps-pph(ib))
      enddo
      zpc = 0.
      do jj = 2, npt
       zpc = zpc+0.5*(max(zptb(jj-1)/zg,0.)+max(zptb(jj)/zg,0.))*(zptp(jj-1)-zptp(jj))
      enddo
      pinfo(2) = zzh(ib)
      pinfo(3) = ztopz
      pinfo(4) = pph(ib)
      pinfo(5) = exp(mfinterp(ztopz))
      pinfo(6) = zpc
      pinfo(10) = zwmean
      zmb = 0.
      if (itype == 1) then
       ztauc = max(ztopz-zzh(ib),1.)/max(zwmean,0.1)
       ztau = max(tau_mf*ztauc,taumin_mf,pdt)
       zpbl = 0.
       if (lrooted) then
        if (kland == 1) then
         ztaubl = ztauc
        else
         zsum = 0.
         zx = 0.
         do k = kn, 1, -1
          if (pph(k) <= pph(ib)) exit
          zsum = zsum+pwsp(k)*(pph(k)-max(pph(k-1),pph(ib)))
          zx = zx+(pph(k)-max(pph(k-1),pph(ib)))
         enddo
         zub = pwsp(kn)
         if (zx > 0.) zub = zsum/zx
         ztaubl = max(zzh(ib),1.)/max(zub,ublmin_mf)
        endif
        zsum = 0.
        do k = kn, 1, -1
         if (pph(k) <= pph(ib)) exit
         zsum = zsum+(pbt(k)*(1.+zeps1*pq(k))+zeps1*pt(k)*pbq(k))*(pph(k)-max(pph(k-1),pph(ib)))
        enddo
        zx = tstar_mf
        if (zx <= 0.) zx = zg*max(ztopz-zzh(ib),1.)/zcp
        zpbl = ztaubl/zx*zsum
       endif
!      the rate at which the unit tendencies reduce PCAPE, over the levels where the updraught is buoyant
       lbuoy(:) = .false.
       do jj = 1, npt
        if (iptk(jj) == 1 .and. zptb(jj) > 0.) lbuoy(iptx(jj)) = .true.
       enddo
       zrate = 0.
       do k = 1, kn
        if (lbuoy(k)) zrate = zrate+(zdsdt(k)*(1.+zeps1*pq(k))+zeps1*pt(k)*zdqdt(k))/ztv(k)*zdp(k)
       enddo
       if (zrate > 0.) zmb = max(0.,(zpc-zpbl)/ztau/zrate)
       pinfo(7) = zpbl
       pinfo(8) = ztau
      elseif (itype == 2) then
       zsupply = 0.
       do k = kn, 1, -1
        if (pph(k) <= pph(ib)) exit
        zsupply = zsupply+(zcp*pbt(k)+zlv*pbq(k))*(pph(k)-max(pph(k-1),pph(ib)))/zg
       enddo
       zhup = zslu(ib)+zlv*zqtu(ib)
       zhenv = zcp*zth(ib)+zg*zzh(ib)+zlv*zqh(ib)
       zexc = zhup-zhenv
       if (zexc > 0.) zmb = max(0.,zsupply/zexc)
      else
       zmb = max(0.,-pom(ksh)/zg)
      endif
!     the mass flux's limit: its net flux may cross no layer faster than cfl_mf times that layer's mass a step
      zlimit = huge(1.)
      do i = 1, kn-1
       zmnet = abs(zmuu(i)+zmd(i))
       if (zmnet > 0.) zlimit = min(zlimit,cfl_mf*min(zmass(i),zmass(i+1))/pdt/zmnet)
      enddo
      if (zmb > zlimit) pinfo(9) = 1.
      zmb = min(zmb,zlimit)
      pinfo(1) = zmb
      call mftend(zmb)
      pdtdt(:) = zdsdt(:)
      pdqdt(:) = zdqdt(:)
      prain(:) = znr(:)
      psnow(:) = zns(:)
      pmu(:) = zmuu(:)*zmb
      pmd(:) = zmd(:)*zmb
      return
!
      contains
!
      real function mfqsat(pt_,pp_)
!     PlaSim's saturation specific humidity
      real, intent(in) :: pt_, pp_
      real :: zq_
      zq_ = rdbrv*ra1*exp(ra2*(pt_-TMELT)/(pt_-ra4))/pp_
      zq_ = min(zq_,rdbrv)
      mfqsat = zq_/(1.-(1./rdbrv-1.)*zq_)
      end function mfqsat
!
      real function mfdqsdt(pt_,pp_)
      real, intent(in) :: pt_, pp_
      real :: zq_, zc_
      zq_ = rdbrv*ra1*exp(ra2*(pt_-TMELT)/(pt_-ra4))/pp_
      zq_ = min(zq_,rdbrv)
      zc_ = 1./(1.-(1./rdbrv-1.)*zq_)
      mfdqsdt = zq_*zc_*zc_*ra2*(TMELT-ra4)/(pt_-ra4)**2
      end function mfdqsdt
!
      real function mfice(pt_)
!     the frozen share of condensate: nil above freezing, all below tice_mf, quadratic between (IFS)
      real, intent(in) :: pt_
      real :: za_
      if (pt_ >= TMELT) then
       mfice = 0.
      elseif (pt_ <= tice_mf) then
       mfice = 1.
      else
       za_ = (pt_-tice_mf)/(TMELT-tice_mf)
       mfice = 1.-za_*za_
      endif
      end function mfice
!
      real function mflat(pt_)
      real, intent(in) :: pt_
      mflat = ALV+mfice(pt_)*(ALS-ALV)
      end function mflat
!
      real function mfvirt(pt_,pq_)
      real, intent(in) :: pt_, pq_
      mfvirt = pt_*(1.+zeps1*pq_)
      end function mfvirt
!
      subroutine mfadjust(psl_,pqt_,pz_,pp_,pt_,pqv_,pl_)
!     temperature, vapour and condensate of air with liquid-ice static energy psl_ and total water pqt_: Newton on
!     cp T + g z - L(T) (qt - qs(T)) - sl = 0 when saturated
      real, intent(in) :: psl_, pqt_, pz_, pp_
      real, intent(out) :: pt_, pqv_, pl_
      real :: ztd_, zqs_, zl_, zdl_, za_, zf_, zdf_, zst_
      integer :: jit_
      ztd_ = (psl_-zg*pz_)/zcp
      if (pqt_ <= mfqsat(ztd_,pp_)) then
       pt_ = ztd_
       pqv_ = pqt_
       pl_ = 0.
       return
      endif
      pt_ = ztd_
      do jit_ = 1, 12
       zqs_ = mfqsat(pt_,pp_)
       zl_ = mflat(pt_)
       zdl_ = 0.
       if (tice_mf < pt_ .and. pt_ < TMELT) then
        za_ = (pt_-tice_mf)/(TMELT-tice_mf)
        zdl_ = -2.*za_/(TMELT-tice_mf)*(ALS-ALV)
       endif
       zf_ = zcp*pt_+zg*pz_-zl_*(pqt_-zqs_)-psl_
       zdf_ = zcp-zdl_*(pqt_-zqs_)+zl_*mfdqsdt(pt_,pp_)
       zst_ = zf_/zdf_
       pt_ = pt_-zst_
       if (abs(zst_) < 1.e-7) exit
      enddo
      zqs_ = mfqsat(pt_,pp_)
      pqv_ = zqs_
      pl_ = max(pqt_-zqs_,0.)
      end subroutine mfadjust
!
      real function mfmix(pphi_,pea_,peb_,peps_,pdz_)
!     a plume's value after rising pdz_ while mixing in surrounding air at peps_ per metre, the surrounding value
!     linear from pea_ to peb_: the exact solution of d phi/dz = eps (phi_env - phi)
      real, intent(in) :: pphi_, pea_, peb_, peps_, pdz_
      real :: zgr_, zdec_
      if (peps_ <= 0. .or. pdz_ <= 0.) then
       mfmix = pphi_
       return
      endif
      zgr_ = (peb_-pea_)/pdz_
      zdec_ = exp(-peps_*pdz_)
      mfmix = peb_-zgr_/peps_+(pphi_-pea_+zgr_/peps_)*zdec_
      end function mfmix
!
      real function mfconv(pl_,pt_,pw_)
!     the rate per metre at which condensate turns to rain or snow (IFS (6.35)-(6.37))
      real, intent(in) :: pl_, pt_, pw_
      real :: zcbf_, zlc_, zw_
      if (pl_ <= lmin_mf) then
       mfconv = 0.
       return
      endif
      zcbf_ = 1.
      if (pt_ < tbf_mf) zcbf_ = 1.+0.5*sqrt(min(tbf_mf-pt_,tbf_mf-tice_mf))
      zlc_ = lcrit_mf/zcbf_
      zw_ = min(max(pw_,0.05),wmaxp_mf)
      mfconv = c0_mf*zcbf_/(0.75*zw_)*(1.-exp(-(pl_/zlc_)**2))
      end function mfconv
!
      subroutine mfwetbulb(pt_,pq_,pp_,ptw_,pqw_)
!     air brought to saturation by evaporating water into it at constant moist static energy
      real, intent(in) :: pt_, pq_, pp_
      real, intent(out) :: ptw_, pqw_
      real :: zh_, zf_, zst_
      integer :: jit_
      zh_ = zcp*pt_+zlv*pq_
      ptw_ = pt_
      do jit_ = 1, 12
       zf_ = zcp*ptw_+zlv*mfqsat(ptw_,pp_)-zh_
       zst_ = zf_/(zcp+zlv*mfdqsdt(ptw_,pp_))
       ptw_ = ptw_-zst_
       if (abs(zst_) < 1.e-7) exit
      enddo
      pqw_ = mfqsat(ptw_,pp_)
      end subroutine mfwetbulb
!
      subroutine mfexcess(pdt_,pdq_,pws_)
!     the surface test parcels' excess from the surface fluxes, 1.5 times each over rho sigma_w (IFS (6.19)), with
!     sigma_w (6.20) at the parcel's source height (or the lowest level's, if that is lower)
      real, intent(out) :: pdt_, pdq_, pws_
      real :: zrho_, zfl_, zz_
      zrho_ = pph(kn)/(gascon*mfvirt(pt(kn),pq(kn)))
      zfl_ = (pshf/zcp+zeps1*pt(kn)*plhf/zlv)/zrho_
      zz_ = zzf(kn)
      if (zsrc_mf > 0.) zz_ = min(zsrc_mf*9.80665/zg,zz_)
      pws_ = (1.2*ustar_mf**3+1.5*0.4*zg/pt(kn)*zz_*max(zfl_,0.))**(1./3.)
      pdt_ = 1.5*max(pshf,0.)/(zrho_*zcp*pws_)
      pdq_ = 1.5*max(plhf,0.)/(zrho_*zlv*pws_)
      end subroutine mfexcess
!
      real function mfpsi(pz_)
!     the integrated stability function for heat in unstable air, phi_h = (1 - 16 zeta)^(-1/3)
      real, intent(in) :: pz_
      real :: zy_
      zy_ = (1.-16.*pz_)**(1./3.)
      mfpsi = 1.5*log((zy_*zy_+zy_+1.)/3.)-sqrt(3.)*atan((2.*zy_+1.)/sqrt(3.))+4.*atan(1.)/sqrt(3.)
      end function mfpsi
!
      subroutine mfsurface(pdt_,pdq_)
!     the air at the source height against the lowest level's, by Monin-Obukhov similarity in a convective
!     surface layer
      real, intent(out) :: pdt_, pdq_
      real :: zzs_, zrho_, zwt_, zwq_, zex_, zthe_, zwtv_, zob_, zf_
      pdt_ = 0.
      pdq_ = 0.
      zzs_ = zsrc_mf*9.80665/zg
      if (zzs_ <= 0. .or. zzs_ >= zzf(kn)) return
      zrho_ = pph(kn)/(gascon*mfvirt(pt(kn),pq(kn)))
      zwt_ = pshf/(zrho_*zcp)
      zwq_ = plhf/(zrho_*zlv)
      zex_ = (ppf(kn)/zps)**(gascon/zcp)
      zthe_ = pt(kn)/zex_
      zwtv_ = zwt_*(1.+zeps1*pq(kn))+zeps1*zthe_*zwq_
      if (zwtv_ <= 0.) return
      zob_ = -ustar_mf**3*zthe_*(1.+zeps1*pq(kn))/(0.4*zg*zwtv_)
      zf_ = (log(zzf(kn)/zzs_)-mfpsi(zzf(kn)/zob_)+mfpsi(zzs_/zob_))/0.4
      pdt_ = zwt_/ustar_mf*zf_*zex_
      pdq_ = zwq_/ustar_mf*zf_
      end subroutine mfsurface
!
      subroutine mfpath(kstart_,lhalf_)
!     the plume's path: full and half levels alternately from full level kstart_ (or half level kstart_ with
!     lhalf_) up to the top full level, with the surroundings there
      integer, intent(in) :: kstart_
      logical, intent(in) :: lhalf_
      integer :: k_
      npath = 0
      if (lhalf_) then
       npath = npath+1
       ipk(npath) = 2
       ipi(npath) = kstart_
       zpz(npath) = zzh(kstart_)
       zpp(npath) = pph(kstart_)
       zpte(npath) = zth(kstart_)
       zpqe(npath) = zqh(kstart_)
       k_ = kstart_+1
      else
       k_ = kstart_+1
      endif
      do while (k_ > 1)
       k_ = k_-1
       npath = npath+1
       ipk(npath) = 1
       ipi(npath) = k_
       zpz(npath) = zzf(k_)
       zpp(npath) = ppf(k_)
       zpte(npath) = pt(k_)
       zpqe(npath) = pq(k_)
       if (k_ > 1) then
        npath = npath+1
        ipk(npath) = 2
        ipi(npath) = k_-1
        zpz(npath) = zzh(k_-1)
        zpp(npath) = pph(k_-1)
        zpte(npath) = zth(k_-1)
        zpqe(npath) = zqh(k_-1)
       endif
      enddo
      end subroutine mfpath
!
      subroutine mfascend(psl_,pqt_,pkin_,pzb_,lshal_,lrain_,lsub_,kbase_,ptop_)
!     one plume along the path from its first point: dry to its condensation level (mixing at the sub-cloud rate
!     with lsub_), then mixing, detraining and raining (with lrain_), until its kinetic energy runs out
      real, intent(in) :: psl_, pqt_, pkin_, pzb_
      logical, intent(in) :: lshal_, lrain_, lsub_
      integer, intent(out) :: kbase_
      real, intent(out) :: ptop_
      real :: zsl_, zqt_, zm_, zkin_, zta_, zqva_, zla_, ztv0_, zzb_, zfac_
      real :: zza_, zzbb_, zpa_, zpb_, ztea_, zteb_, zqea_, zqeb_, zdz_, zrhm_, zep_, zde_, zlp_
      real :: ztb_, zqvb_, zlb_, zform_, zca_, zl0_, zdec_, zkept_, ztvb_, zbb_, zmu_, zbm_, zk2_, zkinb_
      real :: zmb_, zshare_, zmm_, zf0_, zf1_
      integer :: j_, js_, nsub_
      logical :: lstop_, lcloudy_
      zsl_ = psl_
      zqt_ = pqt_
      zm_ = 1.
      zkin_ = pkin_
      call mfadjust(zsl_,zqt_,zpz(1),zpp(1),zta_,zqva_,zla_)
      ztv0_ = mfvirt(zpte(1),zpqe(1))
      npt = 1
      iptk(1) = ipk(1)
      iptx(1) = ipi(1)
      iptl(1) = 0
      zptz(1) = zpz(1)
      zptp(1) = zpp(1)
      zptt(1) = zta_
      zptqv(1) = zqva_
      zptl(1) = zla_
      zptsl(1) = zsl_
      zptqt(1) = zqt_
      zptm(1) = zm_
      zptkin(1) = zkin_
      zptb(1) = zg*(mfvirt(zta_,zqva_)-zla_*zta_-ztv0_)/ztv0_
      zptrain(1) = 0.
      zptsnow(1) = 0.
      kbase_ = 0
      if (zla_ > 0.) kbase_ = 1
      zzb_ = pzb_
      ptop_ = -1.
      nsub_ = max(1,nsub_mf)
      zfac_ = 1./(fturb_mf*(1.+gvm_mf))
      lstop_ = .false.
      do j_ = 2, npath
       do js_ = 0, nsub_-1
        zf0_ = real(js_)/real(nsub_)
        zf1_ = real(js_+1)/real(nsub_)
        zza_ = zpz(j_-1)+zf0_*(zpz(j_)-zpz(j_-1))
        zzbb_ = zpz(j_-1)+zf1_*(zpz(j_)-zpz(j_-1))
        zpa_ = exp(log(zpp(j_-1))+zf0_*log(zpp(j_)/zpp(j_-1)))
        zpb_ = exp(log(zpp(j_-1))+zf1_*log(zpp(j_)/zpp(j_-1)))
        ztea_ = zpte(j_-1)+zf0_*(zpte(j_)-zpte(j_-1))
        zteb_ = zpte(j_-1)+zf1_*(zpte(j_)-zpte(j_-1))
        zqea_ = zpqe(j_-1)+zf0_*(zpqe(j_)-zpqe(j_-1))
        zqeb_ = zpqe(j_-1)+zf1_*(zpqe(j_)-zpqe(j_-1))
        zdz_ = zzbb_-zza_
        zrhm_ = 0.5*(zqea_+zqeb_)/mfqsat(0.5*(ztea_+zteb_),sqrt(zpa_*zpb_))
        lcloudy_ = (kbase_ > 0)
        if (lcloudy_) then
         zep_ = epsbg_mf+epsnb_mf*exp(-max(0.5*(zza_+zzbb_)-zzb_,0.)/lnb_mf)
         if (lshal_) zep_ = zep_*shfac_mf
         zde_ = zep_*(1.6-min(zrhm_,1.))
        elseif (lsub_) then
         zep_ = subf_mf*(epsbg_mf+epsnb_mf)
         zde_ = zep_
        else
         zep_ = 0.
         zde_ = 0.
        endif
        zlp_ = zptl(npt)
        zsl_ = mfmix(zsl_,zcp*ztea_+zg*zza_,zcp*zteb_+zg*zzbb_,zep_,zdz_)
        zqt_ = mfmix(zqt_,zqea_,zqeb_,zep_,zdz_)
        call mfadjust(zsl_,zqt_,zzbb_,zpb_,ztb_,zqvb_,zlb_)
        zform_ = 0.
        if (zlb_ > 0. .and. lrain_) then
         zca_ = mfconv(zlb_,ztb_,sqrt(2.*max(zkin_,0.)))
         if (zca_ > 0.) then
          zl0_ = zlp_*exp(-zep_*zdz_)
          zdec_ = exp(-zca_*zdz_)
          zkept_ = zl0_*zdec_+(zlb_-zl0_)/(zca_*zdz_)*(1.-zdec_)
          zform_ = max(zlb_-max(zkept_,0.),0.)
         endif
         if (zform_ > 0.) then
          zqt_ = zqt_-zform_
          zsl_ = zsl_+mflat(ztb_)*zform_
          call mfadjust(zsl_,zqt_,zzbb_,zpb_,ztb_,zqvb_,zlb_)
         endif
        endif
        ztvb_ = mfvirt(zteb_,zqeb_)
        zbb_ = zg*(mfvirt(ztb_,zqvb_)-zlb_*ztb_-ztvb_)/ztvb_
!       kinetic energy, with drag from mixing at the larger of the rates (IFS (6.10))
        zmu_ = max(zep_,zde_)
        zbm_ = 0.5*(zptb(npt)+zbb_)
        if (zmu_ > 0.) then
         zk2_ = 2.*(1.+betacd_mf)*zmu_
         zdec_ = exp(-zk2_*zdz_)
         zkinb_ = zkin_*zdec_+zfac_*zbm_*(1.-zdec_)/zk2_
        else
         zkinb_ = zkin_+zfac_*zbm_*zdz_
        endif
        zmb_ = zm_*exp((zep_-zde_)*zdz_)
        if (lcloudy_ .and. zbb_ < 0. .and. zkinb_ > 0.) zmb_ = min(zmb_,zm_*sqrt(min(zkinb_/zkin_,1.)))
        if (zkinb_ <= 0.) then
         if (zkin_ > zkinb_) then
          ptop_ = zza_+zkin_/(zkin_-zkinb_)*zdz_
         else
          ptop_ = zza_
         endif
         lstop_ = .true.
         zmb_ = 0.
        endif
        zshare_ = mfice(ztb_)
        zmm_ = 0.5*(zm_+zmb_)
        npt = npt+1
        if (js_ == nsub_-1) then
         iptk(npt) = ipk(j_)
        else
         iptk(npt) = 3
        endif
        iptx(npt) = ipi(j_)
        if (ipk(j_) == 1) then
         iptl(npt) = ipi(j_)
        else
         iptl(npt) = ipi(j_)+1
        endif
        zptz(npt) = zzbb_
        zptp(npt) = zpb_
        zptt(npt) = ztb_
        zptqv(npt) = zqvb_
        zptl(npt) = zlb_
        zptsl(npt) = zsl_
        zptqt(npt) = zqt_
        zptm(npt) = zmb_
        zptkin(npt) = max(zkinb_,0.)
        zptb(npt) = zbb_
        zptrain(npt) = zmm_*zform_*(1.-zshare_)
        zptsnow(npt) = zmm_*zform_*zshare_
        zm_ = zmb_
        zkin_ = zkinb_
        if (kbase_ == 0 .and. zlb_ > 0.) then
         kbase_ = npt
         zzb_ = zzbb_
        endif
        if (lstop_) exit
       enddo
       if (lstop_) exit
      enddo
      if (ptop_ < 0.) ptop_ = zptz(npt)
      end subroutine mfascend
!
      real function mfinterp(pz_)
!     log pressure at height pz_, linear in height between the plume's points (the last one's beyond them)
      real, intent(in) :: pz_
      integer :: j_
      if (pz_ <= zptz(1)) then
       mfinterp = log(zptp(1))
       return
      endif
      do j_ = 2, npt
       if (pz_ <= zptz(j_)) then
        mfinterp = log(zptp(j_-1))+(pz_-zptz(j_-1))/(zptz(j_)-zptz(j_-1))*(log(zptp(j_))-log(zptp(j_-1)))
        return
       endif
      enddo
      mfinterp = log(zptp(npt))
      end function mfinterp
!
      subroutine mftrial(kd_,pt0_,pq0_,pw0_,lshal_,lrain_,lok_,pdepth_)
!     one test parcel from full level kd_: whether it saturates at a half level with kinetic energy to spare, and the
!     cloud's depth in pressure from that half level to its top
      integer, intent(in) :: kd_
      real, intent(in) :: pt0_, pq0_, pw0_
      logical, intent(in) :: lshal_, lrain_
      logical, intent(out) :: lok_
      real, intent(out) :: pdepth_
      integer :: kb_, jb_, j_
      real :: ztop_
      call mfpath(kd_,.false.)
      call mfascend(zcp*pt0_+zg*zzf(kd_),pq0_,0.5*pw0_*pw0_,zzf(kd_),lshal_,lrain_,.true.,kb_,ztop_)
      lok_ = .false.
      pdepth_ = 0.
      if (kb_ == 0) return
      jb_ = 0
      do j_ = kb_, npt
       if (iptk(j_) == 2) then
        jb_ = j_
        exit
       endif
      enddo
      if (jb_ == 0) return
      if (zptkin(jb_) <= 0.) return
      lok_ = .true.
      pdepth_ = zptp(jb_)-exp(mfinterp(ztop_))
!     what the updraught starts from, should this parcel decide the type
      zbsl = zptsl(jb_)
      zbqt = zptqt(jb_)
      zbkin = zptkin(jb_)
      zzcond = zptz(kb_)
      ibase = iptx(jb_)
      end subroutine mftrial
!
      subroutine mfkeep(kd_)
!     the deciding parcel's departure level (its cloud base was kept by its trial)
      integer, intent(in) :: kd_
      ksh = kd_
      end subroutine mfkeep
!
      subroutine mfupdraught(lshal_)
!     the updraught from the deciding parcel's base half level at unit mass flux: its half-level mass flux,
!     liquid-ice static energy and total water, its rain and snow by layer, its top and mean speed
      logical, intent(in) :: lshal_
      integer :: kb_, j_, n_
      real :: zsum_
      ib = ibase
      call mfpath(ib,.true.)
      call mfascend(zbsl,zbqt,zbkin,zzcond,lshal_,.true.,.false.,kb_,ztopz)
      zmu(:) = 0.
      zslu(:) = 0.
      zqtu(:) = 0.
      zrainu(:) = 0.
      zsnowu(:) = 0.
      it = ib
      do j_ = 1, npt
       if (iptk(j_) == 2 .and. zptm(j_) > 0.) then
        zmu(iptx(j_)) = zptm(j_)
        zslu(iptx(j_)) = zptsl(j_)
        zqtu(iptx(j_)) = zptqt(j_)
        it = iptx(j_)
       endif
      enddo
      do j_ = 2, npt
       zrainu(iptl(j_)) = zrainu(iptl(j_))+zptrain(j_)
       zsnowu(iptl(j_)) = zsnowu(iptl(j_))+zptsnow(j_)
      enddo
      zsum_ = 0.
      n_ = 0
      do j_ = 1, npt
       if (zptz(j_) <= ztopz+1.e-6) then
        zsum_ = zsum_+sqrt(2.*zptkin(j_))
        n_ = n_+1
       endif
      enddo
      zwmean = 0.
      if (n_ > 0) zwmean = zsum_/real(n_)
      end subroutine mfupdraught
!
      subroutine mfdowndraught(lshal_)
!     the saturated downdraught at unit updraught base flux, from its level of free sinking
      logical, intent(in) :: lshal_
      integer :: kmin_, i_, k_, lfs_
      real :: zhmin_, zh_, ztw_, zqw_, ztm_, zqm_, zlm_, zmdn_, zfall_, zbot_, zdz_, ztd_, zte_, zqe_
      real :: zneed_, ze_, zf_, zsdn_, zqdn_, zpre_(kn), zabove_
      zmd(:) = 0.
      zsd(:) = 0.
      zqd(:) = 0.
      zevd(:) = 0.
      zpre_(:) = zrainu(:)+zsnowu(:)
      if (sum(zpre_) <= 0. .or. lshal_) return
!     the level of least moist static energy below the cloud top
      kmin_ = 0
      zhmin_ = huge(1.)
      do k_ = 1, kn
       if (ppf(k_) > pph(it)) then
        zh_ = zs(k_)+zlv*pq(k_)
        if (zh_ < zhmin_) then
         zhmin_ = zh_
         kmin_ = k_
        endif
       endif
      enddo
      lfs_ = 0
      do i_ = it+1, ib-1
       if (pph(i_) < ppf(kmin_) .or. zmu(i_) <= 0.) cycle
       call mfwetbulb(zth(i_),zqh(i_),pph(i_),ztw_,zqw_)
       call mfadjust(0.5*(zslu(i_)+zcp*ztw_+zg*zzh(i_)),0.5*(zqtu(i_)+zqw_),zzh(i_),pph(i_),ztm_,zqm_,zlm_)
       zabove_ = 0.
       do k_ = 1, i_
        zabove_ = zabove_+zpre_(k_)
       enddo
       if (mfvirt(ztm_,zqm_)-zlm_*ztm_ < mfvirt(zth(i_),zqh(i_)) .and. zabove_ > 0.) then
        lfs_ = i_
        zsdn_ = zcp*ztm_+zg*zzh(i_)
        zqdn_ = zqm_
        exit
       endif
      enddo
      if (lfs_ == 0) return
      zmdn_ = -etad_mf
      zfall_ = 0.
      do k_ = 1, lfs_
       zfall_ = zfall_+zpre_(k_)
      enddo
      zmd(lfs_) = zmdn_
      zsd(lfs_) = zsdn_
      zqd(lfs_) = zqdn_
      zbot_ = zps-botdp_mf
      do i_ = lfs_+1, kn
       k_ = i_
       zdz_ = zzh(i_-1)-zzh(i_)
       if (i_ < kn) then
        zte_ = zth(i_)
        zqe_ = zqh(i_)
       else
        zte_ = pt(kn)
        zqe_ = pq(kn)
       endif
       zsdn_ = mfmix(zsdn_,zcp*zth(i_-1)+zg*zzh(i_-1),zcp*zte_+zg*zzh(i_),epsd_mf,zdz_)
       zqdn_ = mfmix(zqdn_,zqh(i_-1),zqe_,epsd_mf,zdz_)
       ztd_ = (zsdn_-zg*zzh(i_))/zcp
       call mfwetbulb(ztd_,zqdn_,pph(i_),ztw_,zqw_)
       zneed_ = max(zqw_-zqdn_,0.)*(-zmdn_)
       zfall_ = zfall_+zpre_(k_)
       ze_ = min(zneed_,zfall_)
       if (zneed_ > 0.) then
        zf_ = ze_/zneed_
        ztd_ = ztd_+zf_*(ztw_-ztd_)
        zqdn_ = zqdn_+zf_*(zqw_-zqdn_)
        zsdn_ = zcp*ztd_+zg*zzh(i_)
       endif
       zfall_ = zfall_-ze_
       zevd(k_) = zevd(k_)+ze_
       if (i_ == kn) exit
       if (mfvirt(ztd_,zqdn_) > mfvirt(zte_,zqe_)) exit
       if (pph(i_) > zbot_) zmdn_ = -etad_mf*max(zps-pph(i_),0.)/botdp_mf
       zmd(i_) = zmdn_
       zsd(i_) = zsdn_
       zqd(i_) = zqdn_
      enddo
      end subroutine mfdowndraught
!
      subroutine mftend(pscale_)
!     the column's tendencies for the updraught and downdraught at a base flux of pscale_: the convergence of their
!     fluxes against the surroundings' compensating motion taken upwind, the fluxes falling linearly in pressure
!     below the base, the rain and snow formed as sources and the downdraught's evaporation, the falling water's
!     phase following mkrain (all rain above freezing in ptph, all snow at or below it); heating over PlaSim's
!     moist heat capacity
      real, intent(in) :: pscale_
      real :: zmu_(0:kn), zmd_(0:kn), zshb_, zqhb_, zfus_, zfuq_, zsh_, zqh_, zzf_, zr_(kn), zsn_(kn), zev_(kn)
      real :: zpr_, zpsn_, zheat_, ze_, zer_, zes_
      integer :: i_, k_
      zmu_(:) = zmu(:)*pscale_
      zmd_(:) = zmd(:)*pscale_
      zfs(:) = 0.
      zfq(:) = 0.
      call mfup(ib,zmu_(ib)+zmd_(ib),zshb_,zqhb_)
      zfus_ = zmu_(ib)*(zslu(ib)-zshb_)
      zfuq_ = zmu_(ib)*(zqtu(ib)-zqhb_)
      do i_ = 1, kn-1
       if (i_ <= ib) then
        call mfup(i_,zmu_(i_)+zmd_(i_),zsh_,zqh_)
        zfs(i_) = zmu_(i_)*(zslu(i_)-zsh_)+zmd_(i_)*(zsd(i_)-zsh_)
        zfq(i_) = zmu_(i_)*(zqtu(i_)-zqh_)+zmd_(i_)*(zqd(i_)-zqh_)
       else
        zzf_ = (zps-pph(i_))/(zps-pph(ib))
        zmu_(i_) = zmu_(ib)*zzf_
        call mfup(i_,zmu_(i_)+zmd_(i_),zsh_,zqh_)
        zfs(i_) = zzf_*zfus_+zmd_(i_)*(zsd(i_)-zsh_)
        zfq(i_) = zzf_*zfuq_+zmd_(i_)*(zqd(i_)-zqh_)
       endif
      enddo
      zr_(:) = zrainu(:)*pscale_
      zsn_(:) = zsnowu(:)*pscale_
      zev_(:) = zevd(:)*pscale_
      zpr_ = 0.
      zpsn_ = 0.
      do k_ = 1, kn
       zheat_ = zlv*zr_(k_)+zls*zsn_(k_)
       zpr_ = zpr_+zr_(k_)
       zpsn_ = zpsn_+zsn_(k_)
       if (ptph(k_) > TMELT) then
        zpr_ = zpr_+zpsn_
        zpsn_ = 0.
       else
        zpsn_ = zpr_+zpsn_
        zpr_ = 0.
       endif
       ze_ = min(zev_(k_),zpr_+zpsn_)
       if (ptph(k_) > TMELT) then
        zer_ = ze_
        zes_ = 0.
       else
        zer_ = 0.
        zes_ = ze_
       endif
       zpr_ = zpr_-zer_
       zpsn_ = zpsn_-zes_
       zheat_ = zheat_-(zlv*zer_+zls*zes_)
       znr(k_) = zr_(k_)-zer_
       zns(k_) = zsn_(k_)-zes_
       zdsdt(k_) = zg*(zfs(k_)-zfs(k_-1))/zdp(k_)+zg*zheat_/zdp(k_)
       zdqdt(k_) = zg*(zfq(k_)-zfq(k_-1))/zdp(k_)-zg*(zr_(k_)+zsn_(k_)-zer_-zes_)/zdp(k_)
       zdsdt(k_) = zdsdt(k_)/(zcp*(1.+ADV*pq(k_)))
      enddo
      end subroutine mftend
!
      subroutine mfup(ki_,pmnet_,psh_,pqh_)
!     the surrounding air crossing half level ki_, upwind: from the layer above when the net cloud flux rises
      integer, intent(in) :: ki_
      real, intent(in) :: pmnet_
      real, intent(out) :: psh_, pqh_
      if (pmnet_ > 0.) then
       psh_ = zs(ki_)
       pqh_ = pq(ki_)
      else
       psh_ = zs(ki_+1)
       pqh_ = pq(ki_+1)
      endif
      end subroutine mfup
!
      end subroutine mfcolumn
"""
MASSFLUX_CODE = WRAPPER_CODE + DTT_CODE + COLUMN_CODE
# mkrain's floor under nprc 3 (massflux.evaporate is its copy).
MKRAIN_FLOOR = f"""!
!*    {MASSFLUX_MARKER}'s floor (nprc = 3): in each layer mkrain may evaporate of the convective water at most the
!     least running total of the water handed on from that layer down, less what it has evaporated above, so the
!     downdraught's evaporation, already taken from the handed water, never leaves a running total negative
!
      zfloor(:,:)=1.E30
      ztaken(:)=0.
      if(nprc == 3) then
       zrem(:)=0.
       do jlev=1,NLEV
        zrem(:)=zrem(:)+dprcl(:,jlev)
        zfloor(:,jlev)=zrem(:)
       enddo
       do jlev=NLEM,1,-1
        zfloor(:,jlev)=MIN(zfloor(:,jlev),zfloor(:,jlev+1))
       enddo
      endif
"""
MASSFLUX_EDITS = [
    ("      real :: cin_bm  = 0.    ! its limit on the parcel's inhibition (J/kg; 0: none)\n",
     "      real :: cin_bm  = 0.    ! its limit on the parcel's inhibition (J/kg; 0: none)\n" + declarations()),
    ("      namelist/rainmod_nl/tau_bm,rh_bm,eps_bm,cin_bm   ! Terluna: the simplified Betts-Miller scheme\n",
     "      namelist/rainmod_nl/tau_bm,rh_bm,eps_bm,cin_bm   ! Terluna: the simplified Betts-Miller scheme\n"
     + namelist_text()),
    ("      call mpbcr(cin_bm)\n", "      call mpbcr(cin_bm)\n" + broadcast_text()),
    ("      if(nprc == 2) call mksbm   ! Terluna: the simplified Betts-Miller scheme\n",
     "      if(nprc == 2) call mksbm   ! Terluna: the simplified Betts-Miller scheme\n"
     f"      if(nprc == 3) call mkmassflux   ! {MASSFLUX_MARKER}\n"),
    ("      end subroutine mkdqtgp\n", "      end subroutine mkdqtgp\n" + MASSFLUX_CODE),
    ("      real :: zlcp(NHOR)\n",
     "      real :: zlcp(NHOR)\n"
     f"      real :: zfloor(NHOR,NLEV), ztaken(NHOR), zrem(NHOR) ! {MASSFLUX_MARKER}'s floor (nprc = 3)\n"),
    ("      zdqdtl(:,:)=0.\n", "      zdqdtl(:,:)=0.\n" + MKRAIN_FLOOR),
    ("       if(nevapprec.eq.1) then\n",
     "       if(nevapprec.eq.1) then\n         zrem(:)=MAX(zfloor(:,jlev)-ztaken(:),0.)\n"),
    ("     &                   ,zprc(:))\n           zprc(:)=zprc(:)-zdqdt(:)\n",
     "     &                   ,zprc(:),zrem(:))\n           zprc(:)=zprc(:)-zdqdt(:)\n"
     "           zrem(:)=zrem(:)-zdqdt(:)\n           ztaken(:)=ztaken(:)+zdqdt(:)\n"),
    ("     &                   ,zprsc(:))\n           zprsc(:)=zprsc(:)-zdqdt(:)\n",
     "     &                   ,zprsc(:),zrem(:))\n           zprsc(:)=zprsc(:)-zdqdt(:)\n"
     "           zrem(:)=zrem(:)-zdqdt(:)\n           ztaken(:)=ztaken(:)+zdqdt(:)\n"),
]
MASSFLUX_PATCHES = [('rainmod.f90', MASSFLUX_MARKER, MASSFLUX_EDITS)]

# The surface deep test parcel's starting speed (2026-10-26; massflux.md, the third round): a second patch, over the
# mass-flux patch, applied for a run that sets it and kept in the source after, as the Betts-Miller scheme's options
# are. mfstart gives the boundary layer's convective velocity w* (massflux.boundary_layer_velocity), and with nstart_mf
# 1 the surface deep test parcel starts at cstart_mf w* in place of the larger of w_mf and the excess's velocity scale
# (massflux.surface_start). A parcel may then start with no kinetic energy, so the organized detrainment's ratio of
# kinetic energies is taken only where the energy below is positive, as it always is with nstart_mf 0, the default,
# under which mfcolumn computes as before.
START_MARKER = "Terluna: the surface parcel's starting speed (mass-flux scheme)"
# (Fortran name, massflux.PARAMETERS key, Fortran type, what the variable is)
START_PARAMETER_NAMES = (
    ('nstart_mf', 'start', 'int', "1: the surface deep test parcel starts at cstart_mf w*"),
    ('cstart_mf', 'c_start', 'real', "its starting speed in units of the boundary layer's w*"),
)
START_CHOICES = {'ifs': 0, 'boundary_layer': 1}


def start_values(prm=None) -> dict:
    """The starting-speed parameters' Fortran values for the Python copy's parameters prm (massflux.PARAMETERS by
    default): start 'ifs' as 0 and 'boundary_layer' as 1."""
    prm = dict(mf.PARAMETERS) if prm is None else prm
    if prm['start'] not in START_CHOICES:
        raise ValueError(f"the starting speed is one of {tuple(START_CHOICES)}, not {prm['start']!r}")
    return dict(nstart_mf=START_CHOICES[prm['start']], cstart_mf=float(prm['c_start']))


def start_declarations(values=None) -> str:
    """The starting-speed parameters' declarations, the source's defaults (start 'ifs') unless values are given."""
    values = start_values() if values is None else values
    return (f'!     {START_MARKER} (nprc = 3)\n'
            + ''.join(f"{'      ' + ('integer' if kind == 'int' else 'real') + ' :: ' + name + ' = ' + fortran_literal(values[name]):42s}"
                      f' ! {what}\n' for name, _, kind, what in START_PARAMETER_NAMES))


START_CODE = """
!     =====================
!     SUBROUTINE MFSTART
!     =====================
!
      subroutine mfstart(pdtl_,pdql_,pwst_)
!     Terluna: the boundary layer's convective velocity w* = (g / theta_v (w'theta_v')_s z_i)^(1/3) (Deardorff 1970),
!     z_i the depth of the dry convective layer the surface parcel's source air rises through: where the surroundings'
!     virtual potential temperature first exceeds its own, or where it saturates lifted dry, whichever is lower, each
!     linear in height between full levels; nil without a rising surface buoyancy flux (massflux.boundary_layer_velocity)
      real, intent(in) :: pdtl_, pdql_
      real, intent(out) :: pwst_
      real :: zrho_, zfl_, zkap_, zt0_, zq0_, zthv0_, zthva_, zthvb_, zda_, zdb_, zzt_, zzl_
      integer :: k_
      pwst_ = 0.
      zrho_ = pph(kn)/(gascon*mfvirt(pt(kn),pq(kn)))
      zfl_ = (pshf/zcp+zeps1*pt(kn)*plhf/zlv)/zrho_
      if (zfl_ <= 0.) return
      zkap_ = gascon/zcp
      zt0_ = pt(kn)+pdtl_
      zq0_ = pq(kn)+pdql_
      zthv0_ = mfvirt(zt0_,zq0_)*(zps/ppf(kn))**zkap_
      zthvb_ = mfvirt(pt(kn),pq(kn))*(zps/ppf(kn))**zkap_
      if (zthvb_ >= zthv0_) then
       zzt_ = zzf(kn)
      else
       zzt_ = zzf(1)
       do k_ = kn-1, 1, -1
        zthva_ = mfvirt(pt(k_),pq(k_))*(zps/ppf(k_))**zkap_
        if (zthva_ > zthv0_) then
         zzt_ = zzf(k_+1)+(zthv0_-zthvb_)/(zthva_-zthvb_)*(zzf(k_)-zzf(k_+1))
         exit
        endif
        zthvb_ = zthva_
       enddo
      endif
      zdb_ = mfqsat(zt0_,ppf(kn))-zq0_
      if (zdb_ <= 0.) then
       zzl_ = zzf(kn)
      else
       zzl_ = zzf(1)
       do k_ = kn-1, 1, -1
        zda_ = mfqsat(zt0_*(ppf(k_)/ppf(kn))**zkap_,ppf(k_))-zq0_
        if (zda_ <= 0.) then
         zzl_ = zzf(k_+1)+zdb_/(zdb_-zda_)*(zzf(k_)-zzf(k_+1))
         exit
        endif
        zdb_ = zda_
       enddo
      endif
      pwst_ = (zg/zthv0_*zfl_*min(zzt_,zzl_))**(1./3.)
      end subroutine mfstart
!
"""
# Edits to mfcolumn, which the column's tests apply to COLUMN_CODE (start_column_code), and to rainmod's module.
START_CODE_EDITS = [
    ("      real :: zx, zub, zsum\n", "      real :: zx, zub, zsum, zwst\n"),
    ("      call mfsurface(zdtl,zdql)\n",
     "      call mfsurface(zdtl,zdql)\n      zwst = 0.\n      if (nstart_mf == 1) call mfstart(zdtl,zdql,zwst)\n"),
    ("        zx = max(w_mf,zws)\n", "        zx = max(w_mf,zws)\n        if (nstart_mf == 1) zx = cstart_mf*zwst\n"),
    ("        if (lcloudy_ .and. zbb_ < 0. .and. zkinb_ > 0.) zmb_ = min(zmb_,zm_*sqrt(min(zkinb_/zkin_,1.)))\n",
     "        if (lcloudy_ .and. zbb_ < 0. .and. zkinb_ > 0. .and. zkin_ > 0.) "
     "zmb_ = min(zmb_,zm_*sqrt(min(zkinb_/zkin_,1.)))\n"),
    ("      subroutine mfpath(kstart_,lhalf_)\n", START_CODE.lstrip('\n') + "      subroutine mfpath(kstart_,lhalf_)\n"),
]
START_MODULE_EDITS = [
    (declarations().splitlines(True)[-1], declarations().splitlines(True)[-1] + start_declarations()),
    (namelist_text().splitlines(True)[-1],
     namelist_text().splitlines(True)[-1] + '      namelist/rainmod_nl/'
     + ','.join(name for name, _, _, _ in START_PARAMETER_NAMES) + '\n'),
    ("      call mpbci(nsub_mf)\n", "      call mpbci(nsub_mf)\n" + ''.join(
        f"      call {'mpbci' if kind == 'int' else 'mpbcr'}({name})\n" for name, _, kind, _ in START_PARAMETER_NAMES)),
]
START_EDITS = START_MODULE_EDITS + START_CODE_EDITS
START_PATCHES = [('rainmod.f90', START_MARKER, START_EDITS)]


def start_column_code() -> str:
    """mfcolumn with the starting-speed edits, as the patched source holds it (for the column's tests)."""
    text = COLUMN_CODE
    for old, new in START_CODE_EDITS:
        if text.count(old) != 1:
            raise RuntimeError(f'the column code does not hold {old.strip()!r} once')
        text = text.replace(old, new)
    return text
