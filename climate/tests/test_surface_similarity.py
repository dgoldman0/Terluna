"""CM1's similarity over water in PlaSim's surface transfer (climate/gcm/surface_similarity.py): its tables are Paulson's
and its neutral limit the log law; carried from step to step it settles to a fixed point, and over the warm lunar sea in
light wind it transfers more heat and vapour than Louis's factors with the gust, while land and sea ice keep PlaSim's
transfer; PlaSim's surflx with the patch, compiled with ExoPlaSim's flags against a stub of PlaSim's module, gives the
Python copy's transfer and carried state to 1e-10 on lunar and Earth-like columns from fresh and settled states; with
the switch off it prints every digit the gust's surflx prints; the patch applies once over the gust's; the patched
fluxmod.f90 compiles with PlaSim's modules; and the runner records the similarity and writes its namelist."""
import math
import re
import shutil
import subprocess

import numpy as np
import pytest

from climate.gcm import exoplasim_run as er
from climate.gcm import pbl_fortran as pbf
from climate.gcm import pbl_scheme as ps
from climate.gcm import surface_gust as sg
from climate.gcm import surface_similarity as ssf
from climate.tests.test_pbl_fortran import FORTRAN, GCM_FLAGS, apply, num, twilight_fluxmod, undo
from climate.tests.test_surface_gust import STUB, cases, surflx, with_gust


def without_similarity(text):
    return undo(text, ssf.SIMW_EDITS, ssf.SIMW_MARKER)


def with_similarity(text):
    return apply(with_gust(without_similarity(text)), ssf.SIMW_EDITS, ssf.SIMW_MARKER)


def settle(x, sigma, half, c, steps=200):
    """The column's transfer after the similarity's state has been carried to its fixed point on fixed air."""
    state = dict(ssf.START)
    for _ in range(steps):
        out = ssf.surface_transfer(x, sigma, half, c, state=state)
        state = out['similarity']['state']
    return out, state


# ------------------------------------------------------------------------------------------------------- Python

def test_the_tables_are_paulsons_and_the_neutral_limit_is_the_log_law():
    assert ssf.PSIM_TABLE[0] == pytest.approx(0.0, abs=1e-15) and ssf.PSIH_TABLE[0] == pytest.approx(0.0, abs=1e-15)
    x = 17.0 ** 0.25                                                        # z/L = -1
    assert ssf.PSIM_TABLE[100] == pytest.approx(2 * math.log((1 + x) / 2) + math.log((1 + x * x) / 2)
                                                - 2 * math.atan(x) + math.pi / 2, rel=1e-12)
    assert ssf.PSIH_TABLE[100] == pytest.approx(2 * math.log((1 + math.sqrt(17.0)) / 2), rel=1e-12)
    state = dict(ssf.START)
    for _ in range(200):                                                    # neutral air: thv the same at both ends
        s = ssf.similarity_column(300.0, 295.0, 295.0, 296.0, 296.0, 290.0, 4.0, 1.62, state)
        state = s['state']
    assert s['regime'] == 3
    assert state['ust'] == pytest.approx(ssf.KARMAN * 4.0 / math.log(300.0 / state['z0']), rel=1e-9)
    assert s['transh'] == pytest.approx(state['ust'] * ssf.KARMAN / math.log((300.0 + s['z0q']) / s['z0q']), rel=1e-9)
    assert state['z0'] == pytest.approx(ssf.CHARNOCK * state['ust'] ** 2 / 1.62 + ssf.SMOOTH / state['ust'], rel=1e-9)
    assert ssf.Z0Q_MIN <= s['z0q'] <= ssf.Z0Q_MAX


def test_it_settles_and_parts_from_louiss_factors_with_the_wind():
    """Calm and unstable over the warm lunar sea, the similarity carries more heat and vapour than Louis's factors with
    the gust; at 9 m/s it carries less, its roughness for heat and vapour (Fairall's) far below the momentum's."""
    c, half, cols = cases()[0]                                   # lunar, fifteen levels
    sigma = ps.full_sigma(half)
    for j in range(5):                                            # the sea's five columns
        x = cols[j]
        out, state = settle(x, sigma, half, c)
        again = ssf.surface_transfer(x, sigma, half, c, state=state)
        assert again['dtransh'] == pytest.approx(out['dtransh'], rel=1e-9)
        assert again['similarity']['state']['ust'] == pytest.approx(state['ust'], rel=1e-9)
        assert out['absu2'] == sg.surface_transfer(x, sigma, half, c, ngust=1)['absu2']   # the same wind, the gust's
    calm, windy = (settle(cols[j], sigma, half, c)[0] for j in (1, 2))
    louis_calm, louis_windy = (sg.surface_transfer(cols[j], sigma, half, c, ngust=1) for j in (1, 2))
    assert calm['similarity']['regime'] == 4 and calm['dtransh'] > 1.3 * louis_calm['dtransh']
    assert windy['dtransh'] < 0.85 * louis_windy['dtransh']
    assert windy['similarity']['z0q'] < 0.1 * windy['similarity']['state']['z0']
    for j in (5, 6, 7):                                           # land keeps PlaSim's transfer
        x = cols[j]
        out = ssf.surface_transfer(x, sigma, half, c)
        louis = sg.surface_transfer(x, sigma, half, c, ngust=1)
        assert 'similarity' not in out and out['dtransh'] == louis['dtransh'] and out['dtransm'] == louis['dtransm']
    ice = dict(cols[0], dicec=0.3)                                # so does sea ice
    assert 'similarity' not in ssf.surface_transfer(ice, sigma, half, c)
    off = ssf.surface_transfer(cols[0], sigma, half, c, nsimw=0)
    assert off['dtransh'] == sg.surface_transfer(cols[0], sigma, half, c, ngust=1)['dtransh']


# ------------------------------------------------------------------------------------------------------- Fortran

SIM_STUB = STUB.replace('{declarations}', '      real :: dicec(NHOR) = 0.\n{declarations}')
HARNESS = """
      program harness
      use fluxmod
      implicit none
      integer :: jhor, jlev
      read(*,*) ga, gascon, akap, rdbrv, acpd, adv
      read(*,*) sigma, sigmah
      read(*,*) zumin, cfree, ngust, gust_beta, gust_umin, nsimw
      do jhor=1,NHOR
       read(*,*) dp(jhor), dshfl(jhor), devap(jhor), dust3(jhor), dz0(jhor), dls(jhor), dicec(jhor)
       read(*,*) dsimust(jhor), dsimmol(jhor), dsimz0(jhor)
       read(*,*) (dt(jhor,jlev), jlev=1,NLEP)
       read(*,*) (dq(jhor,jlev), jlev=1,NLEP)
       read(*,*) (du(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dv(jhor,jlev), jlev=1,NLEV)
      enddo
      call surflx
      do jhor=1,NHOR
       write(*,'(*(es25.16e3))') dtransh(jhor), dtransm(jhor), dtsa(jhor), dsimust(jhor), dsimmol(jhor), dsimz0(jhor)
      enddo
      end program harness
"""


def compile_similarity(folder, code, nhor, nlev):
    folder.mkdir(parents=True, exist_ok=True)
    decl = sg.declarations() + ssf.declarations() + pbf.declarations()
    files = {'stub.f90': SIM_STUB.format(nhor=nhor, nlev=nlev, declarations=decl), 'surflx.f90': code,
             'harness.f90': HARNESS}
    for name, text in files.items():
        (folder / name).write_text(text)
    subprocess.run(['gfortran', *GCM_FLAGS, '-o', 'harness', *files], cwd=folder, check=True, capture_output=True)
    return folder / 'harness'


def harness_input(c, half, cols, states, nsimw):
    sigma = ps.full_sigma(half)
    akap = c['gascon'] / c['acpd']
    lines = [num([c['ga'], c['gascon'], akap, c['rdbrv'], c['acpd'], c['ADV']]), num(np.r_[sigma, half]),
             f'{sg.EXOPLASIM_ZUMIN!r} {sg.EXOPLASIM_CFREE!r} 1 {sg.GUST_BETA!r} {sg.GUST_UMIN!r} {int(nsimw)}']
    for x, s in zip(cols, states):
        lines.append(num([x['dp'], x['dshfl'], x['devap'], x['dust3'], x['dz0'], x['dls'], x.get('dicec', 0.0)]))
        lines.append(num([s['ust'], s['mol'], s['z0']]))
        lines += [num(x[k]) for k in ('dt', 'dq', 'du', 'dv')]
    return '\n'.join(lines) + '\n'


def run(exe, text, ncols):
    rows = subprocess.run([str(exe)], input=text, capture_output=True, text=True, check=True).stdout.strip('\n')
    rows = rows.split('\n')
    assert len(rows) == ncols
    return rows


def columns_and_states(c, half, cols):
    """Each case's columns with a sea-ice twin of the first, from the fresh state and from a settled one."""
    sigma = ps.full_sigma(half)
    cols = cols + [dict(cols[0], dicec=0.4)]
    fresh = [dict(ssf.START) for _ in cols]
    settled = [settle(x, sigma, half, c, steps=40)[1] if ssf.applies(x) else dict(ssf.START) for x in cols]
    return cols, (fresh, settled)


@FORTRAN
def test_the_fortran_surflx_gives_the_python_copys_transfer_and_state(tmp_path):
    source = twilight_fluxmod()
    if source is None:
        pytest.skip('no ExoPlaSim source at hand')
    code = surflx(with_similarity(source.read_text()), last='simwcol')
    for i, (c, half, cols) in enumerate(cases()):
        cols, state_sets = columns_and_states(c, half, cols)
        exe = compile_similarity(tmp_path / f'case{i}', code, len(cols), half.size)
        for states in state_sets:
            rows = run(exe, harness_input(c, half, cols, states, 1), len(cols))
            for x, s, row in zip(cols, states, rows):
                h, m, tsa, ust, mol, z0 = (float(v) for v in row.split())
                py = ssf.surface_transfer(x, ps.full_sigma(half), half, c, state=s)
                assert h == pytest.approx(py['dtransh'], rel=1e-10) and m == pytest.approx(py['dtransm'], rel=1e-10)
                assert tsa == pytest.approx(py['dtsa'], rel=1e-12)
                new = py['similarity']['state'] if 'similarity' in py else s
                assert ust == pytest.approx(new['ust'], rel=1e-10) and z0 == pytest.approx(new['z0'], rel=1e-10)
                assert mol == pytest.approx(new['mol'], rel=1e-10, abs=1e-14)


@FORTRAN
def test_with_the_switch_off_it_prints_every_digit_the_gusts_surflx_prints(tmp_path):
    source = twilight_fluxmod()
    if source is None:
        pytest.skip('no ExoPlaSim source at hand')
    text = source.read_text()
    old = surflx(with_gust(without_similarity(text)), last='gustcol')
    new = surflx(with_similarity(text), last='simwcol')
    assert old != new
    for i, (c, half, cols) in enumerate(cases()):
        cols, (fresh, _) = columns_and_states(c, half, cols)
        exe_old = compile_similarity(tmp_path / f'old{i}', old, len(cols), half.size)
        exe_new = compile_similarity(tmp_path / f'new{i}', new, len(cols), half.size)
        a = run(exe_old, harness_input(c, half, cols, fresh, 0), len(cols))
        b = run(exe_new, harness_input(c, half, cols, fresh, 0), len(cols))
        assert a == b


def test_the_patch_applies_once_over_the_gusts():
    assert not any(ssf.SIMW_MARKER in old for old, _ in ssf.SIMW_EDITS)
    assert all(len(line) <= 132 for _, new in ssf.SIMW_EDITS for line in new.splitlines())
    # it stands on the gust's last declaration, the gust's namelist names and broadcast, and the gust's last routine
    gust_text = ''.join(new for _, new in sg.GUST_EDITS)
    assert all(old in gust_text for old, _ in ssf.SIMW_EDITS if old != "       dtransh(jhor)=ztrans*zrifh(jhor)\n")
    source = twilight_fluxmod()
    if source is None:
        pytest.skip('no ExoPlaSim source at hand')
    text = with_similarity(source.read_text())
    assert apply(text, ssf.SIMW_EDITS, ssf.SIMW_MARKER) == text and without_similarity(text) != text
    assert with_gust(without_similarity(text)) == without_similarity(text)
    assert text.count('        call simwcol(jhor,') == 1 and text.count('      subroutine simwcol(') == 1
    for name in ssf.NAMELIST_NAMES:
        assert len(re.findall(rf'integer :: {name} = ', text)) == 1
        assert len(re.findall(rf'call mpbci\({name}\)', text)) == 1
        assert len(re.findall(rf',gust_umin,{name}\b', text)) == 1


@FORTRAN
def test_the_patched_fluxmod_compiles_with_plasims_modules(tmp_path):
    source = twilight_fluxmod()
    if source is None:
        pytest.skip('no ExoPlaSim source at hand')
    src = source.parent
    (tmp_path / 'resmod.f90').write_text('      module resmod\n      parameter(NLAT_ATM = 32)\n'
                                         '      parameter(NLEV_ATM = 15)\n      parameter(NPRO_ATM = 1)\n'
                                         '      end module resmod\n')
    shutil.copy(src / 'plasimmod.f90', tmp_path / 'plasimmod.f90')
    (tmp_path / 'fluxmod.f90').write_text(with_similarity(source.read_text()))
    for name in ('resmod.f90', 'plasimmod.f90', 'fluxmod.f90'):
        subprocess.run(['gfortran', *GCM_FLAGS, '-c', name], cwd=tmp_path, check=True, capture_output=True)


def test_the_runner_records_the_similarity_and_writes_its_namelist():
    s = er.surface_gust_setting('beljaars_similarity')
    assert s['beta'] == sg.GUST_BETA and s['similarity_over_water']['constants']['charnock'] == ssf.CHARNOCK
    assert s['similarity_over_water']['start'] == ssf.START
    assert 'similarity_over_water' not in er.surface_gust_setting('beljaars')
    text = (er.HERE / 'exoplasim_run.py').read_text()
    assert "'fluxmod_namelist', 'NSIMW', '1'" in text
    assert "the similarity over water goes over the gust's patch" in text
    assert '+ (ssf.SIMW_PATCHES if similarity else [])' in text
