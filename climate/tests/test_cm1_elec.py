"""Checks of the electrified CM1's patches and case files (climate/crm/cm1_elec.py): every patch marks its file, the
patches fit CM1 r22 and the pinned WRF-ELEC files, and an electrified case gets the NSSL scheme, the charge tracers,
the settings terluna_elec.F reads and the files leakage and Takahashi's charging need."""
import re
import json

import pytest

from climate.crm import cm1_elec as e
from climate.crm import cm1_run as c


def test_every_electrified_patch_adds_its_marker():
    for name, marker, edits in e.PATCHES:
        assert any(marker in new for _, new, _ in edits), name


def patched(name):
    """A CM1 file, or one of WRF-ELEC's, with the patches every OpenMP build takes and then the electrified ones."""
    tree = c.CM1_HOME / c.SOURCE['version'] / 'src'
    folder = c.CM1_HOME / 'wrf4-elec' / e.ELEC_SOURCE['commit'][:12]
    pinned = {'module_mp_nssl_2mom.F': folder / 'module_mp_nssl_2mom_elec.F',
              'module_discharge_msz.F': folder / 'module_discharge_msz.F'}
    if not tree.exists() or not pinned.get(name, tree).exists():
        pytest.skip('the CM1 source or the WRF-ELEC file is missing')
    if name in pinned:
        text = e.sources(c.CM1_HOME)[name]()
    else:
        text = (tree / name).read_text(encoding='latin-1')
        for pname, marker, edits in [c.MAKEFILE['omp'], *c.PATCHES]:
            if pname == name:
                text = c.apply_patch(text, name, marker, edits)
    for pname, marker, edits in e.PATCHES:
        if pname == name:
            text = c.apply_patch(text, name, marker, edits)
    return text


def directive_names(directive: str, clause: str = r'(?:FIRSTPRIVATE|PRIVATE|REDUCTION)') -> set:
    flat = directive.replace('&\n!$OMP', '')
    clauses = re.findall(rf'\b{clause}\((?:[+a-z]+:)?([^)]*)\)', flat)
    return {name.strip() for c_ in clauses for name in c_.split(',')}


def check_directive(routine: str, directive: str):
    """Every variable the directive names is one of the routine's own, and nothing private is set before the loop (a
    private copy starts unset; firstprivate copies and sums start from the value set there)."""
    for name in directive_names(directive):
        assert re.search(rf'\b{name}\b', routine, re.I), name
    body = routine[routine.index('implicit none'):]
    for name in directive_names(directive, r'PRIVATE'):
        for line in body.splitlines():
            code = line.split('!')[0]
            if '::' in code or re.match(r'\s*do\s', code, re.I):
                continue
            assert not re.search(rf'(^|[^a-z0-9_]){name}\s*(\(.*\))?\s*=[^=]', code, re.I), (name, line)


def test_cm1_calls_the_reordering_driver_with_the_charge_tracers_and_the_time():
    driver = patched('mp_driver.F')
    assert driver.count('call terluna_nssl_elec(') == 3 and 'call nssl_2mom_driver(' not in driver
    assert driver.count('kts = 1,kte = nk, pt3d = pt3d, mtime = mtime)') == 3
    assert driver.count('getdbz,getvt,getsed,getqdiags,dotbud,doqbud,pt3d,mtime)') == 1
    main = patched('cm1.F')
    assert 'doqbud,pt3d,mtime)' in main and 'terluna_elec_step' not in main
    # the vertical field the inductive charging takes from the step before is kept with each restart
    write = main.index('call     restart_write(')
    assert write < main.index('call terluna_elec_restart_write(mtime)') < main.index("next rsttim = ")
    assert patched('param.F').count('ipelec_tmp=nint(var6),isaund_tmp=nint(var7)') == 2
    assert 'pta(i,j,k,n)=0.001' not in patched('init3d.F')
    make = patched('Makefile')
    assert ('\tterluna_lightning.F \\\n\tterluna_mlint2.F \\\n\tmodule_discharge_msz.F \\\n\tterluna_branched.F \\\n'
            '\tterluna_screen.F \\\n\tterluna_elec.F \\\n') in make
    assert 'mp_driver.o: terluna_elec.o' in make and 'terluna_branched.o: module_discharge_msz.o' in make
    assert 'cm1.o: terluna_elec.o' in make


def test_the_module_takes_cm1_s_settings_and_runs_both_drivers_slabs_in_parallel():
    module = patched('module_mp_nssl_2mom.F')
    assert module.count('!$OMP PARALLEL DO') == 2
    assert 'IF ( present( isaund_tmp ) ) isaund = isaund_tmp' in module
    assert 'real, parameter :: gr = TERLUNA_G' in module
    assert {'an', 'elec2', 'sciona2d', 'ctghin', 'zmaxsed', 'scwmin', 'tke2d'} <= directive_names(e.OMP_DIRECTIVE)
    main = module[module.index('SUBROUTINE nssl_2mom_driver('):]
    check_directive(main[:main.index('!$OMP PARALLEL DO')], e.OMP_DIRECTIVE)
    sed = module[module.index('SUBROUTINE nssl_2mom_sed_driver('):]
    check_directive(sed[:sed.index('!$OMP PARALLEL DO')], e.SED_OMP_DIRECTIVE)
    # the sedimentation driver prints its charge totals once a step, at the last sub-step
    assert "ipelec > 0 .and. lastlooptmp ) THEN !{ ! Terluna" in sed
    # and hands cloud ice's effective radius on whole; terluna_elec.F caps what the radiation reads at RRTMG's 140
    # microns and logs the ice the cap cuts
    assert 'MIN(t2(ix,1,kz), 200.E-6)' not in sed and sed.count('MIN(t2(ix,1,kz), 1.E-2)') == 1
    driver = (e.FORTRAN / 'terluna_elec.F').read_text()
    assert 'real, parameter :: rmax = 140.0e-6' in driver and 'wrei(i,k,j) = min(wrei(i,k,j), rmax)' in driver
    assert "call open_log('terluna_ice_optics.txt', unit)" in driver


def test_branched_lightning_runs_on_one_process_and_reports_to_its_unit():
    msz = patched('module_discharge_msz.F')
    routine = msz[msz.index('      subroutine lightmsz('):msz.index('      end subroutine lightmsz')]
    assert not re.search(r'(?i)\bwrite *\( *(6|0|\*) *,', routine[routine.index('       iunit = iunit0'):])
    assert 'mpitotoutdp(1:2) = mpitotindp(1:2) ! Terluna: one process' in routine     # outside its MPI blocks
    assert 'CALL TASK_PROC_MAP()' not in routine and 'ndebug = 0 ! Terluna' in msz
    assert msz.count('use module_boxmgsetup, only: igslg0, jgslg0, kgslg0, terluna_ebrk_lo, terluna_ebrk_hi, &') == 2
    # a ground strike by the rule needs the leader's crossing where a case sets its internal field
    assert routine.index('tcross = terluna_crosses(trje(in,1,lpot)') < routine.index('if (  idownward == 1 .and. tcross')
    assert 'Max( terluna_ebrk_lo, Min( ebrkdp(ix,jy,kz), terluna_ebrk_hi ) )' in routine
    assert routine.index('deallocate( ebrkd, ebrkdp, zlev, cghgt, t2, t4 )') < routine.index('allocate( ebrkd(')
    boxmg = c.CM1_HOME / 'wrf4-elec' / e.ELEC_SOURCE['commit'][:12] / 'module_boxmgsetup.F'
    if boxmg.exists():
        mlint2 = e.sources(c.CM1_HOME)['terluna_mlint2.F']()
        assert '      subroutine mlint2' in mlint2 and mlint2.rstrip().endswith('end')
        assert 'module ' not in mlint2.lower().replace('module_boxmgsetup.f', '')


def test_fall_speeds_scale_only_where_each_particle_is_present_and_not_at_earth_gravity():
    text = e.GRAVITY_TEXT
    assert 'tlg = TERLUNA_G/9.81' in text and 'IF ( tlg /= 1.0 ) THEN' in text
    for species in ('lr', 'ls', 'li', 'lh', 'lhl'):
        guard = text.index(f'IF ( qx(mgs,{species}) > qxmin({species}) )')
        assert guard < text.index(f'vtxbar(mgs,{species},:) = vtxbar(mgs,{species},:)*')
    for species in ('lh', 'lhl'):                                    # the exponent exists only where the particle does
        assert text.index(f'IF ( qx(mgs,{species}) > qxmin({species}) ) THEN') < text.index(f'bxx(mgs,{species})')
    module = patched('module_mp_nssl_2mom.F')
    setvtz = module[module.index('SUBROUTINE setvtz'):]
    assert setvtz.index("      if ( qx(mgs,lh) .gt. qxmin(lh) ) then") < setvtz.index('bxx(mgs,lh) = mmgraupvt(indxr,3)')


def test_an_electrified_case_takes_nssl_and_the_charge_tracers():
    hail = e.namelist_settings({})
    assert hail['param2'] == dict(ptype=27, iptra=1, npt=7, pdtra=0)
    assert hail['param8'] == dict(var1=0.0, var2=0.0, var3=0.0, var4=-1.0, var5=0.0, var6=3.0, var7=12.0, var8=3.0, var9=0.0,
                                  var10=12000.0)
    graupel = e.namelist_settings(dict(hail=False, ipelec=2, lightning=2, leakage=1, radius_m=6000.0, corona_v_m=3000.0,
                                       leader_v_m=1000.0))
    assert graupel['param2'] == dict(ptype=26, iptra=1, npt=6, pdtra=0)
    assert graupel['param8'] == dict(var1=1000.0, var2=3000.0, var3=0.0, var4=-1.0, var5=0.0, var6=2.0, var7=12.0, var8=2.0, var9=1.0,
                                     var10=6000.0)
    assert e.namelist_settings(dict(screen=1, isaund=11))['param8']['var3'] == 1.0
    moon = e.namelist_settings(e.LUNAR)
    assert moon['param8']['var4'] == -1.0 and moon['param8']['var5'] == 6.8           # WRF-ELEC's own ground rule
    assert moon['param8']['var8'] == 4.0                     # branched lightning without the 180-kV/m cap
    assert moon['param8']['var9'] == 1.0                     # leakage through the air's conductivity
    assert moon['param8']['var1'] == 1000.0                  # a ground strike needs the leader to cross


def test_lightning_2_and_4_lift_only_the_breakdown_field_s_cap():
    text = (e.FORTRAN / 'terluna_elec.F').read_text()
    assert '      if( ilight .eq. 2 .or. ilight .eq. 4 ) terluna_ebrk_hi = 1.0e30\n' in text
    assert 'terluna_ebrk_lo =' not in text


def test_leakage_and_lunar_screening_write_the_conductivity_by_height(tmp_path):
    if not e.CONDUCTIVITY.exists():
        pytest.skip('the conductivity product is missing')
    assert 'conductivity' in e.run_files(tmp_path, dict(screen=2), ' &param0\n /\n')[1]
    assert 'conductivity' not in e.run_files(tmp_path, dict(screen=1), ' &param0\n /\n')[1]
    text, record = e.run_files(tmp_path, dict(leakage=1), ' &param0\n /\n')
    assert text == ' &param0\n /\n' and record['conductivity']['cloud'] == '0.1_g_m3'
    lines = (tmp_path / 'terluna_conductivity.txt').read_text().splitlines()
    column = json.loads(e.CONDUCTIVITY.read_text())['solar_minimum']
    assert int(lines[0]) == len(column['z_km']) == len(lines) - 1
    z, clear, cloud = map(float, lines[5].split())
    assert z == pytest.approx(1000.0 * column['z_km'][4])
    assert clear == pytest.approx(column['clear_air_s_m']['100_per_cm3'][4], rel=1e-6)
    assert cloud == pytest.approx(column['cloud_s_m']['0.1_g_m3'][4], rel=1e-6)


def test_takahashi_charging_brings_its_table_and_the_module_s_switch(tmp_path):
    if not e.TAKAHASHI.exists():
        pytest.skip('takahashi.txt is missing (atmosphere/electricity/fetch_inputs.py)')
    text, record = e.run_files(tmp_path, dict(isaund=1), ' &param0\n /\n')
    assert text.endswith(' &nssl_mp_params\n nonigrd = -1,\n /\n') and 'takahashi' in record
    assert (tmp_path / 'takahashi.txt').read_bytes() == e.TAKAHASHI.read_bytes()
