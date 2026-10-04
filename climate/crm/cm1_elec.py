"""CM1 with electrified storms: WRF-ELEC's NSSL two-moment microphysics in place of CM1's own copy, the electric field,
lightning and leakage after it, and the patches that fit them to CM1 and to lunar gravity. Stage 2 of the
atmospheric-electricity study (research/studies/atmospheric_electricity).

WRF-ELEC (MicroTed/wrf4-elec; Mansell et al. 2005, 2010; Fierro et al. 2013) carries charge on every particle type
through every microphysical process, with non-inductive and inductive charging and small ions. Its NSSL module is the
same scheme as CM1 r22's (Mansell maintains both) with the electrification kept and CM1's driver options dropped. The
build fetches it at a pinned commit (in the public domain under the WRF notice), drops it in for CM1's
module_mp_nssl_2mom.F, and adds two Terluna files (climate/crm/fortran): terluna_elec.F passes CM1's calls to
WRF-ELEC's driver with the arrays reordered from CM1's (i,j,k) to WRF's (i,k,j) and keeps the charges in CM1's passive
tracers; after the microphysics it solves for the field and runs lightning and leakage with terluna_lightning.F, which
solves Poisson's equation by FFT and carries WRF-ELEC's cylindrical discharge scheme (light1d). Settings are in the
header of terluna_elec.F; CM1 passes var6 and var7 to the module's set-up as the charging switch and law.

Not carried: CM1's own water budget from the NSSL scheme (the condensation, evaporation and rain totals CM1's copy adds
to qbudget), three-moment arrays, activated CCN and IN, terrain under the field solver, and WRF-ELEC's screening-layer
and 3-D discharge options.
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


FORTRAN = Path(__file__).resolve().parent / 'fortran'
REPO = Path(__file__).resolve().parents[2]
CONDUCTIVITY = REPO / 'atmosphere' / 'electricity' / 'results' / 'conductivity_moon.json'
TAKAHASHI = REPO / 'atmosphere' / 'electricity' / 'inputs' / 'takahashi.txt'

# A case's electricity (its 'elec' entry overrides these): WRF-ELEC's defaults of non-inductive (Saunders and Peck with
# Brooks's critical rime accretion rate) and inductive charging, its cylindrical lightning and breakdown field, no
# leakage, and the NSSL scheme with hail. Leakage takes the Moon's conductivity at solar minimum, for 100 aerosol
# particles per cm3 in clear air and 0.1 g/m3 of cloud water in cloud (atmosphere/electricity/conductivity.py).
SETTINGS = dict(ipelec=3, isaund=12, lightning=1, leakage=0, radius_m=12000.0, hail=True,
                conductivity=dict(sun='solar_minimum', clear_air='100_per_cm3', cloud='0.1_g_m3'))


def sources(home: Path) -> dict:
    """Whole files the electrified builds put in CM1's source before patching: name -> function giving its text."""
    return {'module_mp_nssl_2mom.F': lambda: fetch_module(home).read_text(encoding='latin-1'),
            'terluna_lightning.F': lambda: (FORTRAN / 'terluna_lightning.F').read_text(),
            'terluna_elec.F': lambda: (FORTRAN / 'terluna_elec.F').read_text()}


def namelist_settings(elec: dict) -> dict:
    """Namelist entries of an electrified case, by section: the NSSL scheme (with hail or graupel only), the charges in
    CM1's passive tracers without its positivity limiter, and the settings terluna_elec.F reads."""
    e = dict(SETTINGS, **elec)
    return {'param2': dict(ptype=27 if e['hail'] else 26, iptra=1, npt=7 if e['hail'] else 6, pdtra=0),
            'param8': dict(var6=float(e['ipelec']), var7=float(e['isaund']), var8=float(e['lightning']),
                           var9=float(e['leakage']), var10=float(e['radius_m']))}


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
    if e['leakage']:
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
     '                              its = 1 ,ite = ni, jts = 1,jte = nj, kts = 1,kte = nk, pt3d = pt3d)\n', 3),
    ('                         getdbz,getvt,getsed,getqdiags,dotbud,doqbud)\n',
     '                         getdbz,getvt,getsed,getqdiags,dotbud,doqbud,pt3d)\n', 1),
    ('    logical, intent(in) :: dotbud,doqbud\n',
     '    logical, intent(in) :: dotbud,doqbud\n'
     '    real, intent(inout), dimension(ibp:iep,jbp:jep,kbp:kep,npt) :: pt3d   ! Terluna: the charge tracers\n', 1),
])

# The field, lightning and leakage after each step's microphysics, while the step's new state is in the 3d arrays.
CM1_PATCH = ('cm1.F', 'Terluna: field, lightning and leakage', [
    ('      use mp_driver_module, only : mp_driver\n',
     '      use mp_driver_module, only : mp_driver\n'
     '      use terluna_elec_module, only : terluna_elec_step   ! Terluna: field, lightning and leakage\n', 1),
    ('                         getdbz,getvt,getsed,getqdiags,dotbud,doqbud)\n      endif\n',
     '                         getdbz,getvt,getsed,getqdiags,dotbud,doqbud,pt3d)\n'
     '        call terluna_elec_step(nstep,mtime,dt,zh,zf,rho,q3d,pt3d)\n      endif\n', 1),
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
     '\tmodule_mp_nssl_2mom.F \\\n\tterluna_lightning.F \\\n\tterluna_elec.F \\\n', 1),
    ('mp_driver.o: constants.o input.o misclibs.o', 'mp_driver.o: terluna_elec.o constants.o input.o misclibs.o', 1),
    ('cm1.o: constants.o input.o param.o', 'cm1.o: terluna_elec.o constants.o input.o param.o', 1),
    ('poiss.o: input.o singleton.o\n',
     'poiss.o: input.o singleton.o\n'
     '# Terluna: electrified NSSL\n'
     'terluna_lightning.o: singleton.o\n'
     'terluna_elec.o: input.o constants.o module_mp_nssl_2mom.o terluna_lightning.o\n', 1),
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

PATCHES = [MODULE_PATCH, GRAVITY_PATCH, DRIVER_PATCH, CM1_PATCH, PARAM_PATCH, INIT3D_PATCH, MAKEFILE_PATCH]
