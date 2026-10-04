"""Checks of one flash of WRF-ELEC's branched lightning as the electrified CM1 runs it (climate/crm/fortran/
terluna_branched.F with lightmsz, patched for one process): compiled with CM1's FFT and the field solver, on a small
grid with charge regions stacked over one point. Skipped where gfortran, the CM1 source or WRF-ELEC's files are missing."""
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from climate.crm import cm1_elec as e
from climate.crm import cm1_run as c

HERE = Path(__file__).resolve().parents[1] / 'crm' / 'fortran'
SINGLETON = c.CM1_HOME / c.SOURCE['version'] / 'src' / 'singleton.F'
FOLDER = c.CM1_HOME / 'wrf4-elec' / e.ELEC_SOURCE['commit'][:12]

DRIVER = """
program check
  use module_boxmgsetup
  use terluna_lightning
  use terluna_branched
  implicit none
  integer :: ni, nj, nk, seed, u, nflash, k
  real :: dx, dy, zgrnd
  real, allocatable :: zh1(:), zf1(:), dz(:), q(:,:,:), rho(:,:,:), prs(:,:,:), tk(:,:,:), cloud(:,:,:)
  real, allocatable :: phi(:,:,:), ex(:,:,:), ey(:,:,:), ez(:,:,:), emag(:,:,:), ebrk(:,:,:), dep(:,:,:), total(:,:,:)
  integer, allocatable :: used(:,:,:)
  type(flash_record) :: rec
  open(newunit=u, file='in.bin', access='stream', form='unformatted', status='old')
  read(u) ni, nj, nk, nflash, dx, dy, zgrnd
  allocate( zh1(nk), zf1(nk+1), dz(nk), q(ni,nj,nk), rho(ni,nj,nk), prs(ni,nj,nk), tk(ni,nj,nk), cloud(ni,nj,nk) )
  allocate( phi(ni,nj,nk), ex(ni,nj,nk), ey(ni,nj,nk), ez(ni,nj,nk), emag(ni,nj,nk), ebrk(ni,nj,nk), dep(ni,nj,nk) )
  allocate( total(ni,nj,nk), used(ni,nj,nk) )
  read(u) zh1, zf1, q, rho, prs, tk, cloud
  close(u)
  dz = zf1(2:nk+1) - zf1(1:nk)
  used = 0
  seed = 7
  total = 0.0
  call potential(ni, nj, nk, dx, dy, zh1, zf1, q, phi)
  call field(ni, nj, nk, dx, dy, zh1, zf1, phi, ex, ey, ez, emag)
  call breakdown_field(ni, nj, nk, rho, ebrk)
  open(newunit=u, file='out.txt', status='replace')
  do k = 1,nflash
    call branched_flash(ni, nj, nk, dx, dy, dz, rho, prs, tk, cloud, phi, ex, ey, ez, emag, q, ebrk, used, seed,  &
                        1, 1.0, 200.0, 266.16, zgrnd, 12000.0, 90, dep, rec)
    write(u,*) rec%kind, rec%i0, rec%j0, rec%k0, rec%e0, rec%ebrk0, rec%qpos, rec%qneg, rec%nox, rec%npos, rec%nneg,  &
               rec%area, rec%zlo, rec%zhi
    if( rec%kind .le. 0 ) exit
    total = total + dep
    q = q + dep
    call potential(ni, nj, nk, dx, dy, zh1, zf1, q, phi)
    call field(ni, nj, nk, dx, dy, zh1, zf1, phi, ex, ey, ez, emag)
  enddo
  close(u)
  open(newunit=u, file='dep.bin', access='stream', form='unformatted', status='replace')
  write(u) total
  close(u)
end program check
"""


@pytest.fixture(scope='module')
def program(tmp_path_factory):
    needed = [FOLDER / 'module_discharge_msz.F', FOLDER / 'module_boxmgsetup.F', SINGLETON]
    if shutil.which('gfortran') is None or not all(p.exists() for p in needed):
        pytest.skip('gfortran, the CM1 source or WRF-ELEC\'s files are missing')
    folder = tmp_path_factory.mktemp('branched')
    files = e.sources(c.CM1_HOME)
    msz = files['module_discharge_msz.F']()
    for name, marker, edits in e.PATCHES:
        if name == 'module_discharge_msz.F':
            msz = c.apply_patch(msz, name, marker, edits)
    (folder / 'msz.F').write_text(msz, encoding='latin-1')
    (folder / 'mlint2.F').write_text(files['terluna_mlint2.F'](), encoding='latin-1')
    (folder / 'check.f90').write_text(DRIVER)
    subprocess.run(['gfortran', '-O2', '-ffree-form', '-ffree-line-length-none', '-cpp', str(SINGLETON),
                    str(HERE / 'terluna_lightning.F'), 'mlint2.F', 'msz.F', str(HERE / 'terluna_branched.F'),
                    'check.f90', '-o', 'check'], cwd=folder, check=True, capture_output=True)
    return folder


def column(nk=40, top=20000.0):
    zf1 = np.linspace(0.0, top, nk + 1)
    zh1 = 0.5 * (zf1[:-1] + zf1[1:])
    return zh1, zf1


def air(shape, zh1):
    """A standard-like atmosphere: 300 K at the ground falling 6.5 K/km, pressure and density with a 8-km scale height,
    cloud from 2 to 12 km."""
    ones = np.ones(shape[:2] + (1,))
    tk = (300.0 - 6.5e-3 * zh1)[None, None, :] * ones
    prs = (1.0e5 * np.exp(-zh1 / 8000.0))[None, None, :] * ones
    rho = prs / (287.0 * tk)
    cloud = np.where((zh1 > 2000.0) & (zh1 < 12000.0), 1.0e-3, 0.0)[None, None, :] * ones
    return rho, prs, tk, cloud


def blob(q, i, j, zh1, zc, strength, dx, width=2000.0, depth=1500.0):
    ni, nj, _ = q.shape
    x = (np.arange(ni) - i) * dx
    x = (x + ni * dx / 2) % (ni * dx) - ni * dx / 2                     # distances across the wrapping edge
    y = (np.arange(nj) - j) * dx
    y = (y + nj * dx / 2) % (nj * dx) - nj * dx / 2
    r2 = x[:, None] ** 2 + y[None, :] ** 2
    q += strength * np.exp(-r2 / width ** 2)[:, :, None] * np.exp(-((zh1 - zc) / depth) ** 2)[None, None, :]


def run(folder, q, zh1, zf1, dx, nflash=1, zgrnd=-1.0):
    ni, nj, nk = q.shape
    rho, prs, tk, cloud = air(q.shape, zh1)
    with open(folder / 'in.bin', 'wb') as out:
        np.array([ni, nj, nk, nflash], np.int32).tofile(out)
        np.array([dx, dx, zgrnd], np.float32).tofile(out)
        for a in (zh1, zf1, q, rho, prs, tk, cloud):
            np.asarray(a, np.float32).ravel(order='F').tofile(out)
    done = subprocess.run([str(folder / 'check')], cwd=folder, capture_output=True, text=True)
    assert done.returncode == 0, done.stderr[-2000:]
    rows = [line.split() for line in (folder / 'out.txt').read_text().splitlines()]
    names = ('kind', 'i0', 'j0', 'k0', 'e0', 'ebrk0', 'qpos', 'qneg', 'nox', 'npos', 'nneg', 'area', 'zlo', 'zhi')
    flashes = [dict(zip(names, (float(v) for v in r))) for r in rows]
    dep = np.fromfile(folder / 'dep.bin', np.float32).reshape(nk, nj, ni).transpose(2, 1, 0)
    return flashes, dep


def layout(regions, i=16, j=16):
    """Charge regions (height m, peak C/m3) stacked over one point of a 32 x 32 km domain 20 km deep."""
    zh1, zf1 = column()
    q = np.zeros((32, 32, 40))
    for zc, strength in regions:
        blob(q, i, j, zh1, zc, strength, 1000.0)
    return q, zh1, zf1


DIPOLE = [(9000.0, 3.0e-9), (6000.0, -3.0e-9)]
LOWER = [(11000.0, 1.0e-9), (7500.0, -6.0e-9), (4500.0, 4.0e-9)]     # a strong lower positive charge region


def test_a_flash_between_two_charge_regions_neutralizes_both_and_makes_nitrogen_oxides(program):
    q, zh1, zf1 = layout(DIPOLE)
    flashes, dep = run(program, q, zh1, zf1, 1000.0)
    f = flashes[0]
    assert f['kind'] == 1                                                # in cloud
    assert f['e0'] >= 0.9 * f['ebrk0']
    assert f['qpos'] > 0.0 and f['qneg'] == pytest.approx(f['qpos'], rel=1e-4)   # equal and opposite
    assert np.sum(dep[q > 0]) < 0.0 and np.sum(dep[q < 0]) > 0.0
    assert f['nox'] > 0.0 and f['npos'] > 0 and f['nneg'] > 0
    assert 4000.0 < f['zlo'] < f['zhi'] < 12000.0
    far = np.abs(np.arange(32) - 16)[:, None] + np.abs(np.arange(32) - 16)[None, :] > 12
    assert not np.any(dep[far])                                         # nothing far from the storm


def test_a_flash_across_the_domain_s_edge_reaches_both_sides(program):
    q, zh1, zf1 = layout(DIPOLE, i=0)
    flashes, dep = run(program, q, zh1, zf1, 1000.0)
    assert flashes[0]['kind'] == 1
    changed = np.any(dep != 0.0, axis=2)
    assert changed[:3, 14:19].any() and changed[-3:, 14:19].any()


def test_flashes_continue_until_no_starting_point_is_left(program):
    q, zh1, zf1 = layout(DIPOLE)
    flashes, dep = run(program, q, zh1, zf1, 1000.0, nflash=60)
    assert flashes[-1]['kind'] <= 0 and len(flashes) > 1
    assert all(f['kind'] in (1, 2, 3) for f in flashes[:-1])


def test_a_channel_down_through_a_lower_positive_region_strikes_the_ground_with_negative_charge(program):
    q, zh1, zf1 = layout(LOWER)
    flashes, dep = run(program, q, zh1, zf1, 1000.0)                    # ground below the -7 C level (5.2 km)
    f = flashes[0]
    assert f['kind'] == 2
    assert f['qneg'] > 0.0 and f['qpos'] == 0.0                         # the negative charge goes to the ground
    assert f['npos'] > 0 and f['nneg'] == 0
    assert np.all(dep >= 0.0)


def test_the_mirror_image_storm_strikes_with_positive_charge(program):
    q, zh1, zf1 = layout(LOWER)
    normal, _ = run(program, q, zh1, zf1, 1000.0)
    flashes, dep = run(program, -q, zh1, zf1, 1000.0)
    f, g = flashes[0], normal[0]
    assert f['kind'] == 3
    assert f['qpos'] == pytest.approx(g['qneg'], rel=1e-4) and f['qneg'] == 0.0
    assert (f['npos'], f['nneg']) == (g['nneg'], g['npos'])
    assert np.all(dep <= 0.0)


def test_the_ground_height_decides_whether_the_same_channel_strikes(program):
    q, zh1, zf1 = layout(LOWER)
    high, _ = run(program, q, zh1, zf1, 1000.0, zgrnd=5000.0)          # the lunar setting
    low, _ = run(program, q, zh1, zf1, 1000.0, zgrnd=3000.0)
    assert high[0]['kind'] == 2 and low[0]['kind'] == 1
