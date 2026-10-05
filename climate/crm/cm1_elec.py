"""CM1 with electrified storms: WRF-ELEC's NSSL two-moment microphysics in place of CM1's own copy, the electric field,
lightning and leakage, and the patches that fit them to CM1 and to lunar gravity. Stage 2 of the atmospheric-electricity
study (research/studies/atmospheric_electricity).

WRF-ELEC (MicroTed/wrf4-elec; Mansell et al. 2005, 2010; Fierro et al. 2013) carries charge on every particle type
through every microphysical process, with non-inductive and inductive charging and small ions. Its NSSL module is the
same scheme as CM1 r22's (Mansell maintains both) with the electrification kept and CM1's driver options dropped. The
build fetches it and WRF-ELEC's branched lightning (module_discharge_msz.F) at a pinned commit (in the public domain
under the WRF notice), drops the module in for CM1's module_mp_nssl_2mom.F, and adds three Terluna files
(climate/crm/fortran): terluna_elec.F passes CM1's calls to WRF-ELEC's driver with the arrays reordered from CM1's
(i,j,k) to WRF's (i,k,j), keeps the charges in CM1's passive tracers, and, as WRF-ELEC's own driver does, runs the
sedimentation in sub-steps, each followed by the field and lightning; terluna_lightning.F solves Poisson's equation by
FFT and carries WRF-ELEC's cylindrical discharge (light1d); terluna_branched.F runs one flash of the branched scheme
(lightmsz) on a domain that wraps around; terluna_screen.F ports WRF-ELEC's screening layers at cloud edges. Settings
are in the header of terluna_elec.F; CM1 passes var6 and var7 to the module's set-up as the charging switch and law.

Not carried: CM1's own water budget from the NSSL scheme (the condensation, evaporation and rain totals CM1's copy adds
to qbudget), three-moment arrays, activated CCN and IN, and terrain under the field solver.
"""
from __future__ import annotations
import hashlib
import shutil
import urllib.request
from pathlib import Path

ELEC_SOURCE = dict(
    repository='https://github.com/MicroTed/wrf4-elec', commit='e43041b000ebfd77c4a528f2ca3847746e839aa8',
    url='https://raw.githubusercontent.com/MicroTed/wrf4-elec/e43041b000ebfd77c4a528f2ca3847746e839aa8/elec/'
        'module_mp_nssl_2mom_elec.F',
    sha256='9203693d76b841f176bb51fa1f444dce990212c6590ace2ce204748f1ea0ffc3', bytes=1118061, licence='WRF public-domain notice (LICENSE.txt)')


def fetch_module(home: Path) -> Path:
    """WRF-ELEC's NSSL module at the pinned commit, downloaded once and checked against the pin."""
    folder = home / 'wrf4-elec' / ELEC_SOURCE['commit'][:12]
    path = folder / 'module_mp_nssl_2mom_elec.F'
    if not path.exists():
        folder.mkdir(parents=True, exist_ok=True)
        part = path.with_name(path.name + '.part')
        with urllib.request.urlopen(ELEC_SOURCE['url']) as response, open(part, 'wb') as out:
            shutil.copyfileobj(response, out)
        part.replace(path)
    data = path.read_bytes()
    if len(data) != ELEC_SOURCE['bytes'] or hashlib.sha256(data).hexdigest() != ELEC_SOURCE['sha256']:
        raise RuntimeError(f'{path} does not match the pinned WRF-ELEC module; not using it')
    return path


# WRF-ELEC's branched lightning (lightmsz; MacGorman, Straka and Ziegler 2001), at the same commit
MSZ_SOURCE = dict(
    url='https://raw.githubusercontent.com/MicroTed/wrf4-elec/e43041b000ebfd77c4a528f2ca3847746e839aa8/elec/'
        'module_discharge_msz.F',
    sha256='1cde4951dfc5e01f72bccb3437be6e6c03dd0e239b49b8f6b6485d9217a7c283', bytes=165455)


def fetch_msz(home: Path) -> Path:
    """WRF-ELEC's branched lightning at the pinned commit, downloaded once and checked against the pin."""
    folder = home / 'wrf4-elec' / ELEC_SOURCE['commit'][:12]
    path = folder / 'module_discharge_msz.F'
    if not path.exists():
        folder.mkdir(parents=True, exist_ok=True)
        part = path.with_name(path.name + '.part')
        with urllib.request.urlopen(MSZ_SOURCE['url']) as response, open(part, 'wb') as out:
            shutil.copyfileobj(response, out)
        part.replace(path)
    data = path.read_bytes()
    if len(data) != MSZ_SOURCE['bytes'] or hashlib.sha256(data).hexdigest() != MSZ_SOURCE['sha256']:
        raise RuntimeError(f'{path} does not match the pinned WRF-ELEC branched lightning; not using it')
    return path


def quiet_msz(text: str) -> str:
    """lightmsz's report to standard output and error goes to the unit it is given (terluna_msz.log): every write to
    units 6, 0 or * in the routine after it takes up that unit."""
    import re
    start = text.index('       iunit = iunit0')
    end = text.index('      end subroutine lightmsz')
    body = re.sub(r'(?i)\bwrite *\( *(6|0|\*) *,', 'write(iunit,', text[start:end])
    return text[:start] + body + text[end:]


# WRF-ELEC's trilinear interpolation (mlint2), which lightmsz calls, stands outside the multigrid module in this file
BOXMG_SOURCE = dict(
    url='https://raw.githubusercontent.com/MicroTed/wrf4-elec/e43041b000ebfd77c4a528f2ca3847746e839aa8/elec/'
        'module_boxmgsetup.F',
    sha256='eca6882798713e9b940219ec64aad168fc22ce900149dc33e1754451cdfa11ac', bytes=264072)


def fetch_boxmg(home: Path) -> Path:
    """WRF-ELEC's multigrid set-up at the pinned commit, downloaded once and checked against the pin."""
    folder = home / 'wrf4-elec' / ELEC_SOURCE['commit'][:12]
    path = folder / 'module_boxmgsetup.F'
    if not path.exists():
        folder.mkdir(parents=True, exist_ok=True)
        part = path.with_name(path.name + '.part')
        with urllib.request.urlopen(BOXMG_SOURCE['url']) as response, open(part, 'wb') as out:
            shutil.copyfileobj(response, out)
        part.replace(path)
    data = path.read_bytes()
    if len(data) != BOXMG_SOURCE['bytes'] or hashlib.sha256(data).hexdigest() != BOXMG_SOURCE['sha256']:
        raise RuntimeError(f'{path} does not match the pinned WRF-ELEC multigrid set-up; not using it')
    return path


def mlint2(text: str) -> str:
    """WRF-ELEC's mlint2 as its own file: the subroutine after the multigrid module, to its closing end."""
    start = text.index('      subroutine mlint2')
    end = text.index('\n      return\n      end\n', start) + len('\n      return\n      end\n')
    return ('! WRF-ELEC (MicroTed/wrf4-elec e43041b, elec/module_boxmgsetup.F): mlint2, trilinear interpolation for '
            'lightmsz;\n! taken out by climate/crm/cm1_elec.py\n' + text[start:end])


FORTRAN = Path(__file__).resolve().parent / 'fortran'
REPO = Path(__file__).resolve().parents[2]
CONDUCTIVITY = REPO / 'atmosphere' / 'electricity' / 'results' / 'conductivity_moon.json'
TAKAHASHI = REPO / 'atmosphere' / 'electricity' / 'inputs' / 'takahashi.txt'

# A case's electricity (its 'elec' entry overrides these): WRF-ELEC's defaults of non-inductive (Saunders and Peck with
# Brooks's critical rime accretion rate) and inductive charging, its branched lightning (lightning 3; 1 for its
# cylinders, 2 and 4 for those two without the breakdown field's 180-kV/m cap), its breakdown field and its 0.75-s sub-step
# (substep_s 0), its ground-strike rule (ground_m -1, WRF-ELEC's nssl_zgrnd default: a downward channel reaching air
# warmer than -7 C in charge of the matching sign; lightmsz takes it only below zero, and at zero counts only a channel
# reaching the two lowest levels, as the Earth benchmarks of 2026-10-03/04 ran), no leakage, and
# the NSSL scheme with hail, and no screening layers (screen 1 for WRF-ELEC's, with Earth's conductivity; 2 with the
# conductivity below). Leakage and screening take the Moon's conductivity at solar minimum, for 100 aerosol particles
# per cm3 in clear air and 0.1 g/m3 of cloud water in cloud (atmosphere/electricity/conductivity.py). Point discharge
# from the ground, which WRF-ELEC does not have, is off (corona_v_m 0); otherwise it is its onset field at 1.225 kg/m3
# (Standler and Winn 1979: 3000 over dense vegetation, 5000 on a barren ridge). A downward channel that meets the
# ground-strike rule strikes as WRF-ELEC has it (leader_v_m 0), or only if a leader with that internal field (V/m at
# 1.225 kg/m3, scaled by density) keeps the potential its streamer zone needs all the way to the ground.
SETTINGS = dict(ipelec=3, isaund=12, lightning=3, leakage=0, radius_m=12000.0, hail=True, substep_s=0.0, ground_m=-1.0,
                screen=0, conductivity=dict(sun='solar_minimum', clear_air='100_per_cm3', cloud='0.1_g_m3'),
                corona_v_m=0.0, leader_v_m=0.0)
# The lunar boxes. WRF-ELEC's 0.75-s sub-step lets graupel settle through about 1 % of the 500-m layers of its supercell
# in each; in the boxes' charging zone (25-35 km) the layers are 2 km deep and graupel falls at 0.44 of Earth's speed,
# so 0.75 x 4 / 0.44 = 6.8 s keeps that share. The boxes take WRF-ELEC's own ground-strike rule (ground_m -1: a downward
# channel reaching air warmer than -7 C in matching charge), since the lunar storms hold their charge as Earth's do, with
# -7 C just beneath the main negative charge (the author's decision, 2026-10-05). It counts a channel at -7 C, about
# 34 km up, as reaching the ground, so its ground strikes are an upper bound. The height rule the boxes first took, a
# channel within 5 km of the ground, kept the demand for matching charge where the storms have none and allowed no
# ground strike. WRF-ELEC caps the breakdown field at 180 kV/m, which on Earth applies only below about 4.5 km, beneath
# where flashes start; the denser lunar air would put it below about 36 km, where lunar flashes start, so the boxes
# lift it (lightning 4; the author's decision, 2026-10-04). Charge leaks through stage 1's conductivity (leakage 1; the
# author's decision, 2026-10-05): without it the charge evaporating cloud leaves on the small ions stays for weeks and
# held the field at the ground near 20 kV/m all over the box, where leakage clears it within the hour and leaves the
# lightning as it was.
LUNAR = dict(substep_s=6.8, ground_m=-1.0, lightning=4, leakage=1)


def sources(home: Path) -> dict:
    """Whole files the electrified builds put in CM1's source before patching: name -> function giving its text."""
    return {'module_mp_nssl_2mom.F': lambda: fetch_module(home).read_text(encoding='latin-1'),
            'module_discharge_msz.F': lambda: quiet_msz(fetch_msz(home).read_text(encoding='latin-1')),
            'terluna_mlint2.F': lambda: mlint2(fetch_boxmg(home).read_text(encoding='latin-1')),
            'terluna_lightning.F': lambda: (FORTRAN / 'terluna_lightning.F').read_text(),
            'terluna_branched.F': lambda: (FORTRAN / 'terluna_branched.F').read_text(),
            'terluna_screen.F': lambda: (FORTRAN / 'terluna_screen.F').read_text(),
            'terluna_elec.F': lambda: (FORTRAN / 'terluna_elec.F').read_text()}


def namelist_settings(elec: dict) -> dict:
    """Namelist entries of an electrified case, by section: the NSSL scheme (with hail or graupel only), the charges in
    CM1's passive tracers without its positivity limiter, and the settings terluna_elec.F reads."""
    e = dict(SETTINGS, **elec)
    return {'param2': dict(ptype=27 if e['hail'] else 26, iptra=1, npt=7 if e['hail'] else 6, pdtra=0),
            'param8': dict(var1=float(e['leader_v_m']), var2=float(e['corona_v_m']), var3=float(e['screen']), var4=float(e['ground_m']),
                           var5=float(e['substep_s']), var6=float(e['ipelec']),
                           var7=float(e['isaund']), var8=float(e['lightning']), var9=float(e['leakage']),
                           var10=float(e['radius_m']))}


def run_files(case: Path, elec: dict, namelist: str) -> tuple:
    """Write the files an electrified case reads into its folder and return its namelist with the NSSL module's own
    section where needed, and a record of what was written: Takahashi's table (isaund 1 or 3, with nonigrd = -1) and
    the conductivity by height for leakage."""
    e = dict(SETTINGS, **elec)
    record = {}
    if e['isaund'] in (1, 3):
        if not TAKAHASHI.exists():
            raise RuntimeError(f'{TAKAHASHI} is missing; restore it with atmosphere/electricity/fetch_inputs.py')
        shutil.copyfile(TAKAHASHI, case / 'takahashi.txt')
        namelist = namelist.rstrip('\n') + '\n\n &nssl_mp_params\n nonigrd = -1,\n /\n'
        record['takahashi'] = hashlib.sha256(TAKAHASHI.read_bytes()).hexdigest()[:16]
    if e['leakage'] or e['screen'] == 2:
        import json
        product = json.loads(CONDUCTIVITY.read_text())
        choice = e['conductivity']
        column = product[choice['sun']]
        z = [1000.0 * v for v in column['z_km']]
        clear, cloud = column['clear_air_s_m'][choice['clear_air']], column['cloud_s_m'][choice['cloud']]
        (case / 'terluna_conductivity.txt').write_text(f'{len(z)}\n' + ''.join(
            f'{zz:.1f} {a:.6e} {b:.6e}\n' for zz, a, b in zip(z, clear, cloud)))
        record['conductivity'] = dict(choice, product=str(CONDUCTIVITY.relative_to(REPO)),
                                      product_sha256=hashlib.sha256(CONDUCTIVITY.read_bytes()).hexdigest()[:16])
    return namelist, record


# CM1 reads these of the module's settings (param.F, init3d.F); WRF-ELEC's copy keeps them private.
PUBLIC_LINES = [
    '  real, private :: rho_qr = 1000., cnor = 8.0e5  ! cnor is set in namelist!!  rain params',
    '  real, private :: rho_qs =  100., cnos = 3.0e6  ! set in namelist!!  snow params',
    '  real, private :: rho_qh =  500., cnoh = 4.0e5  ! set in namelist!!  graupel params',
    '  real, private :: rho_qhl=  800., cnohl = 4.0e4 ! set in namelist!!  hail params',
    '  real, private :: hldnmn = 500.0  ! minimum hail density (for variable density hail)',
    '  real   , private :: ccn            = 0.6e+09   ! set in namelist!! Central plains CCN value',
    '  integer, private :: infall = 4   ! 0 -> uses number-wgt for N; NO correction applied (results in excessive size sorting)',
    '  integer, private :: icdx = 6 ! (graupel) 0=Ferrier; 1=leave drag coef. cd fixed; 2=vary by density, 4=set by user with cdxmin,cdxmax,etc.',
    '  integer, private :: icdxhl = 6 ! (hail) 0=Ferrier; 1=leave drag coef. cd fixed; 2=vary by density, 4=set by user with cdxmin,cdxmax,etc.',
    '  integer, private :: inucopt = 0',
    '  integer, private :: irenuc = 5      ! =1 to always allow renucleation of droplets within the cloud (do no use, obsolete)',
    '  integer, private :: iehw = 1            ! 0 -> ehw=ehw0; 1 -> old ehw; 2 -> test ehw with Mason table data',
    '  integer, private :: iehlw = 1           ! 0 -> ehlw=ehlw0; 1 -> old ehlw; 2 -> test ehlw with Mason table data',
    '  real   , private :: ehw0 = 0.9 ! 0.5          ! constant or max assumed graupel-droplet collection efficiency',
    '  real   , private :: ehlw0 = 0.9 ! 0.75        ! constant or max assumed hail-droplet collection efficiency',
    '  real   , private :: dfrz = 0.15e-3 ! 0.25e-3  ! minimum diameter of frozen drops from Bigg freezing (used for vfrz) for iacr > 1',
    '  integer, private ::  ihlcnh = -1  ! which graupel -> hail conversion to use',
    '  integer, private :: imurain = 1 ! 3 for gamma-volume, 1 for gamma-diameter DSD for rain.',
    '  integer, private :: iferwisventr = 2 ! =1 for Ferrier rwvent, =2 for Wisner rwvent (imurain=1)',
    '  integer, private :: dmrauto       = 0 ',
    '  integer, private :: ioldlimiter = 0 ! test switch for new(=0) or old(=1) size limiter at the end of GS for 3-moment categories',
    '  integer, private  :: lccn = 9 ! 0 or 9, other indices adjusted accordingly',
    '  real, private    :: alphah  = 0.0 ! set in namelist!! shape parameter for ZIEG graupel',
    '  real, private    :: alphahl = 1.0 ! set in namelist!! shape parameter for ZIEG hail',
]

OMP_ANCHOR = '          ancuten(its:ite,1,kts:kte,:) = 0.0\n          thproclocal(:,:) = 0.0\n\n\n     DO jy = jts,jye\n'
OMP_DIRECTIVE = (
    '#if defined(OPENMP)\n'
    '! Terluna: the slabs in parallel, as in CM1\'s copy of the module\n'
    '!$OMP PARALLEL DO DEFAULT(SHARED) &\n'
    '!$OMP PRIVATE(ix,jy,kz,il,loopcnt,xfall,axtra2d,an,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t00,t77,dbz2d,vzf2d,dn1,pn,wn, &\n'
    '!$OMP dz2d,dz2dinv,ltemq,ssival,dp1,elec2,sciona2d,t8s,t9s,ssat,qvapor,ssifac,rainprod2d,evapprod2d,kediagloc, &\n'
    '!$OMP alpha2d,hailmax1d,hailmaxk1,tmp,tmpchg,dv,dv1,c1,temp1,refl,flag) &\n'
    '!$OMP FIRSTPRIVATE(ancuten,thproclocal,tke2d) &\n'
    '!$OMP REDUCTION(+:chgiona1,chgiona2,chgiona3,chgneg1,chgpos1,chgneg2,chgpos2,chgneg3,chgpos3,sctot3, &\n'
    '!$OMP scwtot,scrtot,scitot,scstot,schtot,schltot,cwmass1,cwmass2,rwmass1,rwmass2,icemass1,icemass2, &\n'
    '!$OMP swmass1,swmass2,grmass1,grmass2,hlmass1,hlmass2,wvol5,wvol10,ctswin,ctswip,ctswwn,ctswwp, &\n'
    '!$OMP ctghsn,ctghsp,ctghin,ctghip,ctghwn,ctghwp,ctgsn,ctgsp,ctgin,ctgip,cthsn,cthsp,cthin,cthip, &\n'
    '!$OMP timesed1,timesed2,timesed3,timesetvt,timevtcalc) &\n'
    '!$OMP REDUCTION(max:scwmax,scrmax,scimax,scsmax,schmax,schlmax,sctotmax,scionmax,zmaxsed) &\n'
    '!$OMP REDUCTION(min:scwmin,scrmin,scimin,scsmin,schmin,schlmin,sctotmin,scionmin)\n'
    '#endif\n')

SED_OMP_ANCHOR = ('          ancuten(its:ite,1,kts:kte,:) = 0.0\n          thproclocal(:,:) = 0.0\n\n'
                  '     DO jy = jts,jye\n     \n     xfall(:,:,:) = 0.0\n')
SED_OMP_DIRECTIVE = (
    '#if defined(OPENMP)\n'
    '! Terluna: the slabs in parallel, as in the microphysics driver\n'
    '!$OMP PARALLEL DO DEFAULT(SHARED) &\n'
    '!$OMP PRIVATE(ix,jy,kz,il,xfall,axtra2d,an,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t00,t77,dbz2d,vzf2d,dn1,pn,wn, &\n'
    '!$OMP dz2d,dz2dinv,ltemq,elec2,sciona2d,kediagloc,alpha2d,hailmax1d,hailmaxk1,tmp,tmpchg,dv,dv1,refl) &\n'
    '!$OMP FIRSTPRIVATE(ancuten,thproclocal) &\n'
    '!$OMP REDUCTION(+:chgiona1,chgiona2,chgiona3,chgneg1,chgpos1,chgneg2,chgpos2,chgneg3,chgpos3,sctot3, &\n'
    '!$OMP scwtot,scrtot,scitot,scstot,schtot,schltot,cwmass1,cwmass2,rwmass1,rwmass2,icemass1,icemass2, &\n'
    '!$OMP swmass1,swmass2,grmass1,grmass2,hlmass1,hlmass2,wvol5,wvol10,ctswin,ctswip,ctswwn,ctswwp, &\n'
    '!$OMP ctghsn,ctghsp,ctghin,ctghip,ctghwn,ctghwp,ctgsn,ctgsp,ctgin,ctgip,cthsn,cthsp,cthin,cthip, &\n'
    '!$OMP timesed,timesed1,timesed2,timesed3,timesetvt,timevtcalc,timegs,timenucond) &\n'
    '!$OMP REDUCTION(max:scwmax,scrmax,scimax,scsmax,schmax,schlmax,sctotmax,scionmax,zmaxsed) &\n'
    '!$OMP REDUCTION(min:scwmin,scrmin,scimin,scsmin,schmin,schlmin,sctotmin,scionmin)\n'
    '#endif\n')
SED_PRINT = ("#else\n         IF ( present(scw) .and. ipelec > 0 ) THEN !{\n#endif\n           iunit = 6\n"
             "         IF ( mytask == 0 ) THEN\n           write(iunit,'(a,3(2x,1pe13.5))' ) 'pre-sed: pos/neg/tot = '")

MODULE_PATCH = ('module_mp_nssl_2mom.F', 'Terluna: WRF-ELEC NSSL module in CM1', [
    ('', '! Terluna: WRF-ELEC NSSL module in CM1 (MicroTed/wrf4-elec e43041b), patched by climate/crm/cm1_elec.py\n', 1),
    *[(line + '\n', line.replace('private', 'public', 1) + '\n', 1) for line in PUBLIC_LINES],
    # CM1 sets the module's equation set to its own (eqtset), as CM1's copy of the module allows, and the charging switch
    # and law, which in WRF come with nssl_params; the switch also sets up the charge arrays and reads Takahashi's table
    ('     & myrank, mpiroot &\n     )\n', '     & myrank, mpiroot, eqtset_tmp, ipelec_tmp, isaund_tmp &\n     )\n', 1),
    ('   integer, intent(in) :: ipctmp,mixphase\n',
     '   integer, intent(in) :: ipctmp,mixphase\n'
     '   integer, optional, intent(in) :: eqtset_tmp, ipelec_tmp, isaund_tmp ! Terluna: from CM1\n', 1),
    ('    ipconc = ipctmp\n', '    ipconc = ipctmp\n'
     "    IF ( present( eqtset_tmp ) ) eqtset = eqtset_tmp ! Terluna: as CM1's copy\n"
     '    IF ( present( ipelec_tmp ) ) ipelec = ipelec_tmp ! Terluna: charging on or off\n'
     '    IF ( present( isaund_tmp ) ) isaund = isaund_tmp ! Terluna: the non-inductive charging law\n', 1),
    # WRF writes namelist.output before the module adds its own settings to it; CM1 does not
    ("          open(15,file='namelist.output',status='old',action='readwrite', position='append',form='formatted')\n",
     "          open(15,file='namelist.output',status='unknown',action='readwrite', position='append',form='formatted') "
     "! Terluna: CM1 writes none\n", 1),
    # CM1's copy runs its driver's slabs in parallel; WRF-ELEC's leaves that to WRF's tiles. What each slab sets is
    # private to it, the charge and timing totals are summed across threads, and what the driver sets before the loop
    # (infdo, the fall-speed moments) stays shared.
    (OMP_ANCHOR, OMP_ANCHOR.replace('     DO jy = jts,jye\n', OMP_DIRECTIVE + '     DO jy = jts,jye\n'), 1),
    # the sedimentation driver, which the sub-steps call, likewise; its totals are printed at the step's last sub-step
    (SED_OMP_ANCHOR, SED_OMP_ANCHOR.replace('     DO jy = jts,jye\n', SED_OMP_DIRECTIVE + '     DO jy = jts,jye\n'),
     1),
    (SED_PRINT, SED_PRINT.replace('ipelec > 0 ) THEN !{', 'ipelec > 0 .and. lastlooptmp ) THEN !{ ! Terluna'), 1),
    # the sedimentation driver sets the effective radii radiation reads once the sub-steps run; it hands cloud ice's on
    # whole (its 200-micron bound lifted), and terluna_elec.F caps it at RRTMG's 140 microns and logs what it cuts
    ('             re_ice(ix,kz,jy)   = MAX(10.01E-6, MIN(t2(ix,1,kz), 200.E-6))\n',
     '             re_ice(ix,kz,jy)   = MAX(10.01E-6, MIN(t2(ix,1,kz), 1.E-2)) ! Terluna: capped in terluna_elec.F\n', 1),
    ('             IF ( .not. present(qi) ) re_ice(ix,kz,jy)  = MAX(10.E-6, MIN(t3(ix,1,kz), 200.E-6))\n',
     '             IF ( .not. present(qi) ) re_ice(ix,kz,jy)  = MAX(10.E-6, MIN(t3(ix,1,kz), 1.E-2)) ! Terluna\n', 1),
])

# WRF-ELEC's branched lightning on one process: WRF-ELEC runs it only under MPI, where a few lines outside its MPI
# blocks read the tile count. Its breakdown field takes the bounds the build sets (module_boxmgsetup) and is computed
# afresh at each call, since terluna_branched.F recentres the arrays on each flash; its debugging output is off. Where a
# case sets a leader's internal field (var1), a downward channel that meets the ground-strike rule strikes only if the
# leader can cross the rest of the way (module_boxmgsetup's terluna_crosses).
MSZ_PATCH = ('module_discharge_msz.F', 'Terluna: branched lightning on one process', [
    ('', '! Terluna: branched lightning on one process (MicroTed/wrf4-elec e43041b), patched by climate/crm/cm1_elec.py\n',
     1),
    ('      use module_boxmgsetup, only: igslg0, jgslg0, kgslg0\n',
     '      use module_boxmgsetup, only: igslg0, jgslg0, kgslg0, terluna_ebrk_lo, terluna_ebrk_hi, &\n'
     '     &   terluna_leader_ei, terluna_vapour, terluna_crosses ! Terluna\n', 2),
    ('      integer :: idownward\n', '      integer :: idownward\n      logical :: tcross ! Terluna\n', 1),
    ('      if (  idownward == 1 .and. &\n'
     '          (trje(in,it,llz) .lt. zgrnd',
     '      tcross = .true. ! Terluna: a downward leader\'s crossing to the ground (module_boxmgsetup)\n'
     '      IF ( terluna_leader_ei > 0.0 .and. idownward == 1 .and. &\n'
     '          (trje(in,it,llz) .lt. zgrnd .or. (zgrnd < 0. .and. trje(in,it,ltemg) > tgrnd) .or. kc <= 2 ) ) THEN\n'
     '        tcross = terluna_crosses(trje(in,1,lpot), trje(in,1,llz), trje(in,it,llz), nz, pot(ic,1:nz,jc), &\n'
     '                                 db(ic,1:nz,jc), dzz(ic,1:nz,jc), terluna_vapour(ic,1:nz,jc))\n'
     '      ENDIF\n'
     '      if (  idownward == 1 .and. tcross .and. &\n'
     '          (trje(in,it,llz) .lt. zgrnd', 1),
    ('      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION, MPI_SUM, local_communicator, '
     'mpi_error_code)\n',
     '      mpitotoutdp(1:2) = mpitotindp(1:2) ! Terluna: one process\n', 1),
    ('      IF ( ntasks == 0 ) THEN\n', '      IF ( .false. ) THEN ! Terluna: one process; WRF always has a task\n', 1),
    ('      IF ( .false. .and. ntasks == 1 ) THEN ! only do horizontal if single processor -- for now.\n',
     '      IF ( .false. ) THEN ! Terluna: one process (the streamline stays vertical, as in WRF-ELEC)\n', 1),
    ('       IF ( mpi_setup_flag < 0 ) THEN\n         CALL TASK_PROC_MAP()\n       ENDIF\n',
     '       ! Terluna: one process, no MPI task map\n', 1),
    ('      integer, parameter :: ndebug = 1\n', '      integer, parameter :: ndebug = 0 ! Terluna\n', 1),
    ('      IF ( .not. allocated( ebrkd ) ) THEN ! {\n',
     '      IF ( allocated( ebrkd ) ) deallocate( ebrkd, ebrkdp, zlev, cghgt, t2, t4 ) ! Terluna: afresh each call\n'
     '      IF ( .not. allocated( ebrkd ) ) THEN ! {\n', 1),
    ('        ebrkd(ix,jy,kz)= Max( 50.e3, Min( ebrkdp(ix,jy,kz), 180.0e3 ) )\n',
     '        ebrkd(ix,jy,kz)= Max( terluna_ebrk_lo, Min( ebrkdp(ix,jy,kz), terluna_ebrk_hi ) ) ! Terluna\n', 1),
])

# CM1's calls go to terluna_nssl_elec with the charge tracers. WRF-ELEC's driver reports rain as WRF does (RAINNC
# accumulated, RAINNCV this step, mm) and has none of CM1's own rain, budget, map-factor, radiation-flag and diagnostic
# arguments; CM1 converts the step's rain as it does for P3 (mm to cm, and a rate per second).
DRIVER_PATCH = ('mp_driver.F', 'Terluna: WRF-ELEC NSSL driver', [
    ('                              RAIN = rain,                   &\n'
     '                              hail = hail,                   &\n'
     '                              nrain = nrain,                 &\n'
     '                              prate = prate,                 &\n',
     '                              RAINNC = dum5(ib,jb,1), RAINNCV = dum5(ib,jb,2), &\n', 3),
    ('                              ruh = ruh, rvh = rvh, rmh = rmh, &\n', '', 3),
    ('                              tcond = qbudget(1),            &\n'
     '                              tevac = qbudget(2),            &\n'
     '                              tevar = qbudget(5),            &\n'
     '                              train = qbudget(6),            &\n'
     '                              rr    = dum3,                  &\n', '', 3),
    ('                              dorad = dorad,                 &\n', '', 2),
    ('                              dorad = .false.,                 &\n', '', 1),
    ('                  ib3d=ib3d,ie3d=ie3d,jb3d=jb3d,je3d=je3d,kb3d=kb3d,ke3d=ke3d, &\n'
     '                  nout3d=nout3d,out3d=out3d,                                 &\n', '', 3),
    # WRF's domain extents only mark a charging-free strip along the domain's edges; CM1's boxes and rings are periodic
    ('                              ims = ib ,ime = ie , jms = jb ,jme = je, kms = kb,kme = ke,  &\n',
     '                              ids = -1000000, ide = 1000000, jds = -1000000, jde = 1000000, kds = 1, kde = nk+1, &\n'
     '                              ims = ib ,ime = ie , jms = jb ,jme = je, kms = kb,kme = ke,  &\n', 3),
    ('          IF ( ptype .eq. 26 ) THEN  ! graupel only\n',
     '          ! Terluna: WRF-ELEC NSSL driver: this step\'s rain (mm) comes back in dum5(:,:,2)\n'
     '          dum5(:,:,1:2) = 0.0\n'
     '          IF ( ptype .eq. 26 ) THEN  ! graupel only\n', 1),
    ('          ENDIF\n\n        IF(eqtset.eq.2)THEN\n          ! for mass conservation:\n',
     '          ENDIF\n\n'
     '          do j=1,nj\n          do i=1,ni\n'
     '            prate(i,j) = dum5(i,j,2)/dt\n'
     '            do n=1,nrain\n              rain(i,j,n) = rain(i,j,n) + 0.1*dum5(i,j,2)\n            enddo\n'
     '          enddo\n          enddo\n\n'
     '        IF(eqtset.eq.2)THEN\n          ! for mass conservation:\n', 1),
    ('    use module_mp_nssl_2mom, only : nssl_2mom_driver\n',
     "    use terluna_elec_module, only : terluna_nssl_elec ! Terluna: WRF-ELEC's NSSL driver, reordered\n", 1),
    ('             call nssl_2mom_driver(                          &\n',
     '             call terluna_nssl_elec(                         &\n', 3),
    ('                              its = 1 ,ite = ni, jts = 1,jte = nj, kts = 1,kte = nk)\n',
     '                              its = 1 ,ite = ni, jts = 1,jte = nj, kts = 1,kte = nk, pt3d = pt3d, mtime = mtime)\n',
     3),
    ('                         getdbz,getvt,getsed,getqdiags,dotbud,doqbud)\n',
     '                         getdbz,getvt,getsed,getqdiags,dotbud,doqbud,pt3d,mtime)\n', 1),
    ('    logical, intent(in) :: dotbud,doqbud\n',
     '    logical, intent(in) :: dotbud,doqbud\n'
     '    real, intent(inout), dimension(ibp:iep,jbp:jep,kbp:kep,npt) :: pt3d   ! Terluna: the charge tracers\n'
     '    double precision, intent(in) :: mtime   ! Terluna: the time at the step\'s start, for the lightning logs\n', 1),
])

# CM1 passes the charge tracers and the time to the microphysics, and keeps the vertical field the inductive charging
# takes from the step before with each restart (as WRF-ELEC keeps elecz).
CM1_PATCH = ('cm1.F', 'Terluna: charge tracers to the microphysics', [
    ('                         getdbz,getvt,getsed,getqdiags,dotbud,doqbud)\n      endif\n',
     '                         getdbz,getvt,getsed,getqdiags,dotbud,doqbud,pt3d,mtime)   ! Terluna: charge tracers to the microphysics\n'
     '      endif\n', 1),
    ('      use mp_driver_module, only : mp_driver\n',
     '      use mp_driver_module, only : mp_driver\n'
     '      use terluna_elec_module, only : terluna_elec_restart_write   ! Terluna: the field kept with each restart\n', 1),
    ('          ! end_restart_write\n',
     '          ! end_restart_write\n'
     '        call terluna_elec_restart_write(mtime)   ! Terluna: the field kept with each restart\n', 1),
])

# The module's set-up takes the charging switch (var6) and law (var7).
PARAM_PATCH = ('param.F', 'Terluna: charging switch and law', [
    ('                     eqtset_tmp=eqtset,myrank=myid,mpiroot=0)\n',
     '                     eqtset_tmp=eqtset,myrank=myid,mpiroot=0, &\n'
     '                     ipelec_tmp=nint(var6),isaund_tmp=nint(var7))   ! Terluna: charging switch and law\n', 2),
])

# The passive tracers carry charge and start without any.
INIT3D_PATCH = ('init3d.F', 'Terluna: the tracers carry charge', [
    ('pta(i,j,k,n)=0.001', 'pta(i,j,k,n)=0.0   ! Terluna: the tracers carry charge', 3),
])

MAKEFILE_PATCH = ('Makefile', 'Terluna: electrified NSSL', [
    ('\tmodule_mp_nssl_2mom.F \\\n',
     '\tmodule_mp_nssl_2mom.F \\\n\tterluna_lightning.F \\\n\tterluna_mlint2.F \\\n\tmodule_discharge_msz.F \\\n'
     '\tterluna_branched.F \\\n\tterluna_screen.F \\\n\tterluna_elec.F \\\n', 1),
    ('mp_driver.o: constants.o input.o misclibs.o', 'mp_driver.o: terluna_elec.o constants.o input.o misclibs.o', 1),
    ('poiss.o: input.o singleton.o\n',
     'poiss.o: input.o singleton.o\n'
     '# Terluna: electrified NSSL\n'
     'terluna_lightning.o: singleton.o\n'
     'module_discharge_msz.o: terluna_lightning.o\n'
     'terluna_branched.o: module_discharge_msz.o terluna_lightning.o\n'
     'terluna_elec.o: input.o constants.o module_mp_nssl_2mom.o terluna_lightning.o terluna_branched.o terluna_screen.o\n'
     'cm1.o: terluna_elec.o\n', 1),
])

# Lunar gravity. Graupel and hail under the drag laws (icdx 1-5) and cloud droplets (Stokes) fall with gr itself. Rain,
# snow and cloud ice follow fitted laws V = a D^b, as do graupel and hail under WRF-ELEC's default icdx = 6 (the
# Milbrandt and Morrison 2013 table, by particle density) or icdx < 1, so their speeds and the graupel and hail
# coefficient axx (also used for collection and ventilation) scale by (g/9.81)^((b+1)/3), as the Morrison patch scales
# its laws: rain with b = 0.8 as in the Morrison runs, snow with Ferrier's 0.42 (isnowfall = 2), cloud ice by its option
# (icefallopt 1, 2, 3: 1.415, 0.6635, 0.55), graupel and hail with each particle's tabulated exponent bxx. The scheme
# sets a particle's speeds and coefficients only where it is present, so the scaling acts only there; at Earth's 9.81
# m/s2 it leaves them as they are.
GRAVITY_TEXT = (
    '!     Terluna: fall speeds at the host gravity (climate/crm/cm1_elec.py), where each particle type is present\n'
    '      tlg = TERLUNA_G/9.81\n'
    '      IF ( tlg /= 1.0 ) THEN\n'
    '      tlbice = 0.55\n'
    '      IF ( icefallopt == 1 ) tlbice = 1.415\n'
    '      IF ( icefallopt == 2 ) tlbice = 0.6635\n'
    '      DO mgs = 1,ngscnt\n'
    '        IF ( lr > 1 .and. ( ildo == 0 .or. ildo == lr ) ) THEN\n'
    '          IF ( qx(mgs,lr) > qxmin(lr) ) vtxbar(mgs,lr,:) = vtxbar(mgs,lr,:)*tlg**(1.8/3.0)\n'
    '        ENDIF\n'
    '        IF ( ls > 1 .and. ( ildo == 0 .or. ildo == ls ) ) THEN\n'
    '          IF ( qx(mgs,ls) > qxmin(ls) ) vtxbar(mgs,ls,:) = vtxbar(mgs,ls,:)*tlg**(1.42/3.0)\n'
    '        ENDIF\n'
    '        IF ( li > 1 .and. ( ildo == 0 .or. ildo == li ) ) THEN\n'
    '          IF ( qx(mgs,li) > qxmin(li) ) vtxbar(mgs,li,:) = vtxbar(mgs,li,:)*tlg**((tlbice+1.0)/3.0)\n'
    '        ENDIF\n'
    '        IF ( lh > 1 .and. ( ildo == 0 .or. ildo == lh ) .and. ( icdx == 6 .or. icdx <= 0 ) ) THEN\n'
    '          IF ( qx(mgs,lh) > qxmin(lh) ) THEN\n'
    '            tlfac = tlg**((bxx(mgs,lh)+1.0)/3.0)\n'
    '            vtxbar(mgs,lh,:) = vtxbar(mgs,lh,:)*tlfac\n'
    '            axx(mgs,lh) = axx(mgs,lh)*tlfac\n'
    '          ENDIF\n'
    '        ENDIF\n'
    '        IF ( lhl > 1 .and. ( ildo == 0 .or. ildo == lhl ) .and. ( icdxhl == 6 .or. icdxhl <= 0 ) ) THEN\n'
    '          IF ( qx(mgs,lhl) > qxmin(lhl) ) THEN\n'
    '            tlfac = tlg**((bxx(mgs,lhl)+1.0)/3.0)\n'
    '            vtxbar(mgs,lhl,:) = vtxbar(mgs,lhl,:)*tlfac\n'
    '            axx(mgs,lhl) = axx(mgs,lhl)*tlfac\n'
    '          ENDIF\n'
    '        ENDIF\n'
    '      ENDDO\n'
    '      ENDIF\n\n')
GRAVITY_PATCH = ('module_mp_nssl_2mom.F', 'Terluna: fall speeds at the host gravity', [
    ('', '#ifndef TERLUNA_G\n#define TERLUNA_G 9.81\n#endif\n', 1),
    ('      real, parameter :: gr = 9.8\n', '      real, parameter :: gr = TERLUNA_G ! Terluna: host gravity\n', 1),
    ('      real gr\n', '      real gr\n      real tlg, tlfac, tlbice ! Terluna\n', 1),
    ('      gr = 9.8\n', '      gr = TERLUNA_G ! Terluna: host gravity\n', 1),
    ("      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: END OF ROUTINE'\n",
     GRAVITY_TEXT + "      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: END OF ROUTINE'\n", 1),
])

PATCHES = [MODULE_PATCH, GRAVITY_PATCH, MSZ_PATCH, DRIVER_PATCH, CM1_PATCH, PARAM_PATCH, INIT3D_PATCH,
           MAKEFILE_PATCH]
