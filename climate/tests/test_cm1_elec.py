"""Checks of the electrified CM1's patches and case files (climate/crm/cm1_elec.py): every patch marks its file, the
patches fit CM1 r22 and the pinned WRF-ELEC module, and an electrified case gets the NSSL scheme, the charge tracers,
the settings terluna_elec.F reads and the files leakage and Takahashi's charging need."""
import json

import pytest

from climate.crm import cm1_elec as e
from climate.crm import cm1_run as c


def test_every_electrified_patch_adds_its_marker():
    for name, marker, edits in e.PATCHES:
        assert any(marker in new for _, new, _ in edits), name


def patched(name):
    """A CM1 file, or WRF-ELEC's module, with the patches every OpenMP build takes and then the electrified ones."""
    tree = c.CM1_HOME / c.SOURCE['version'] / 'src'
    module = c.CM1_HOME / 'wrf4-elec' / e.ELEC_SOURCE['commit'][:12] / 'module_mp_nssl_2mom_elec.F'
    if not tree.exists() or not module.exists():
        pytest.skip('the CM1 source or the WRF-ELEC module is missing')
    if name == 'module_mp_nssl_2mom.F':
        text = module.read_text(encoding='latin-1')
    else:
        text = (tree / name).read_text(encoding='latin-1')
        for pname, marker, edits in [c.MAKEFILE['omp'], *c.PATCHES]:
            if pname == name:
                text = c.apply_patch(text, name, marker, edits)
    for pname, marker, edits in e.PATCHES:
        if pname == name:
            text = c.apply_patch(text, name, marker, edits)
    return text


def test_cm1_calls_the_reordering_driver_and_follows_it_with_the_field():
    driver = patched('mp_driver.F')
    assert driver.count('call terluna_nssl_elec(') == 3 and 'call nssl_2mom_driver(' not in driver
    assert driver.count('kts = 1,kte = nk, pt3d = pt3d)') == 3
    assert driver.count('getdbz,getvt,getsed,getqdiags,dotbud,doqbud,pt3d)') == 1
    main = patched('cm1.F')
    assert main.index('doqbud,pt3d)') < main.index('call terluna_elec_step(nstep,mtime,dt,zh,zf,rho,q3d,pt3d)')
    assert patched('param.F').count('ipelec_tmp=nint(var6),isaund_tmp=nint(var7)') == 2
    assert 'pta(i,j,k,n)=0.001' not in patched('init3d.F')
    make = patched('Makefile')
    assert '\tterluna_lightning.F \\\n\tterluna_elec.F \\\n' in make and 'mp_driver.o: terluna_elec.o' in make


def test_the_module_takes_cm1_s_settings_and_runs_its_slabs_in_parallel():
    module = patched('module_mp_nssl_2mom.F')
    assert module.count('!$OMP PARALLEL DO') == 1
    assert 'IF ( present( isaund_tmp ) ) isaund = isaund_tmp' in module
    assert 'real, parameter :: gr = TERLUNA_G' in module
    # every variable the directive names is one of the driver's own
    import re
    declarations = module[module.index('SUBROUTINE nssl_2mom_driver('):module.index('!$OMP PARALLEL DO')]
    clauses = re.findall(r'(?:PRIVATE|REDUCTION)\((?:[+a-z]+:)?([^)]*)\)', e.OMP_DIRECTIVE.replace('&\n!$OMP', ''))
    listed = {name.strip() for clause in clauses for name in clause.split(',')}
    assert {'an', 'elec2', 'sciona2d', 'ctghin', 'zmaxsed', 'scwmin', 'tke2d'} <= listed
    for name in listed:
        assert re.search(rf'\b{name}\b', declarations, re.I), name
    # a private copy starts unset, so nothing private may be set before the loop (firstprivate copies start from the
    # value set there, and sums start from it as well)
    flat = e.OMP_DIRECTIVE.replace('&\n!$OMP', '')
    private = {name.strip() for name in re.search(r'\bPRIVATE\(([^)]*)\)', flat).group(1).split(',')}
    body = declarations[declarations.index('implicit none'):]
    for name in private:
        for line in body.splitlines():
            code = line.split('!')[0]
            if '::' in code or re.match(r'\s*do\s', code, re.I):
                continue
            assert not re.search(rf'(^|[^a-z0-9_]){name}\s*(\(.*\))?\s*=[^=]', code, re.I), (name, line)


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
    assert hail['param8'] == dict(var6=3.0, var7=12.0, var8=1.0, var9=0.0, var10=12000.0)
    graupel = e.namelist_settings(dict(hail=False, ipelec=2, lightning=2, leakage=1, radius_m=6000.0))
    assert graupel['param2'] == dict(ptype=26, iptra=1, npt=6, pdtra=0)
    assert graupel['param8'] == dict(var6=2.0, var7=12.0, var8=2.0, var9=1.0, var10=6000.0)


def test_leakage_writes_the_conductivity_by_height(tmp_path):
    if not e.CONDUCTIVITY.exists():
        pytest.skip('the conductivity product is missing')
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
