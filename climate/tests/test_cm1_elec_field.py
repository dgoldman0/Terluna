"""Checks of the electrified CM1's field solver and lightning (climate/crm/fortran/terluna_lightning.F), compiled with
CM1's FFT and run on small grids: the potential against a direct solve of the same difference equations, against the
exact potential of a charged layer between grounded plates, and the cylindrical discharge against WRF-ELEC's rule.
Skipped where gfortran or the CM1 source is missing."""
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from climate.crm import cm1_run

EPS0 = 8.8541878128e-12
HERE = Path(__file__).resolve().parents[1] / 'crm' / 'fortran'
SINGLETON = cm1_run.CM1_HOME / cm1_run.SOURCE['version'] / 'src' / 'singleton.F'

DRIVER = """
program check
  use module_boxmgsetup
  use terluna_lightning
  implicit none
  integer :: ni, nj, nk, option, ninit, ncol, nregion, u
  real :: dx, dy, radius
  real, allocatable :: zh1(:), zf1(:), q(:,:,:), rho(:,:,:), phi(:,:,:), ex(:,:,:), ey(:,:,:), ez(:,:,:), emag(:,:,:)
  real, allocatable :: ebrk(:,:,:), dep(:,:,:)
  double precision :: w, qpos, qneg
  open(newunit=u, file='in.bin', access='stream', form='unformatted', status='old')
  read(u) ni, nj, nk, option, dx, dy, radius
  allocate( zh1(nk), zf1(nk+1), q(ni,nj,nk), rho(ni,nj,nk), phi(ni,nj,nk), ex(ni,nj,nk), ey(ni,nj,nk), ez(ni,nj,nk) )
  allocate( emag(ni,nj,nk), ebrk(ni,nj,nk), dep(ni,nj,nk) )
  read(u) zh1, zf1, q, rho
  close(u)
  call potential(ni, nj, nk, dx, dy, zh1, zf1, q, phi)
  call field(ni, nj, nk, dx, dy, zh1, zf1, phi, ex, ey, ez, emag)
  w = energy(ni, nj, nk, dx, dy, zf1, q, phi)
  if( option .eq. 2 )then                     ! unbounded
    terluna_ebrk_lo = 0.0
    terluna_ebrk_hi = 1.0e30
  endif
  call breakdown_field(ni, nj, nk, rho, ebrk)
  call discharge(ni, nj, nk, dx, dy, zf1, rho, emag, ebrk, radius, q, dep, ninit, ncol, nregion, qpos, qneg)
  open(newunit=u, file='out.bin', access='stream', form='unformatted', status='replace')
  write(u) phi, ex, ey, ez, emag, ebrk, q, dep, w, qpos, qneg, ninit, ncol, nregion
  close(u)
end program check
"""


@pytest.fixture(scope='module')
def program(tmp_path_factory):
    if shutil.which('gfortran') is None or not SINGLETON.exists():
        pytest.skip('gfortran or the CM1 source (singleton.F) is missing')
    folder = tmp_path_factory.mktemp('lightning')
    (folder / 'check.f90').write_text(DRIVER)
    subprocess.run(['gfortran', '-O2', '-fopenmp', '-ffree-form', '-ffree-line-length-none', '-cpp', str(SINGLETON),
                    str(HERE / 'terluna_lightning.F'), 'check.f90', '-o', 'check'], cwd=folder, check=True,
                   capture_output=True)
    return folder


def run(folder, q, zh1, zf1, dx, dy, rho=None, radius=12000.0, option=1):
    ni, nj, nk = q.shape
    rho = np.full(q.shape, 1.0) if rho is None else rho
    with open(folder / 'in.bin', 'wb') as out:
        np.array([ni, nj, nk, option], np.int32).tofile(out)
        np.array([dx, dy, radius], np.float32).tofile(out)
        for a in (zh1, zf1, q, rho):
            np.asarray(a, np.float32).ravel(order='F').tofile(out)
    subprocess.run([str(folder / 'check')], cwd=folder, check=True, capture_output=True,
                   env=dict(os.environ, OMP_NUM_THREADS='2'))
    raw = (folder / 'out.bin').read_bytes()
    n = ni * nj * nk
    arrays = np.frombuffer(raw[:8 * n * 4], np.float32).reshape(8, nk, nj, ni).transpose(0, 3, 2, 1)
    w, qpos, qneg = np.frombuffer(raw[8 * n * 4:8 * n * 4 + 24], np.float64)
    ninit, ncol, nregion = np.frombuffer(raw[8 * n * 4 + 24:], np.int32)
    names = ('phi', 'ex', 'ey', 'ez', 'emag', 'ebrk', 'q', 'dep')
    return dict(zip(names, arrays), energy=w, qpos=qpos, qneg=qneg, ninit=ninit, ncol=ncol, nregion=nregion)


def stretched(nk, dz0, stretch, top):
    zf = [0.0]
    dz = dz0
    while len(zf) <= nk:
        zf.append(zf[-1] + dz)
        dz *= stretch
    zf = np.array(zf) * top / zf[-1]
    return 0.5 * (zf[:-1] + zf[1:]), zf


def direct_potential(q, zh1, zf1, dx, dy):
    """The same difference equations, assembled and solved directly."""
    ni, nj, nk = q.shape
    idx = lambda i, j, k: (i % ni) + ni * ((j % nj) + nj * k)
    a = np.zeros((q.size, q.size))
    for k in range(nk):
        dzk = zf1[k + 1] - zf1[k]
        lo = 1.0 / (dzk * (2.0 * (zh1[0] - zf1[0]) if k == 0 else zh1[k] - zh1[k - 1]))
        up = 1.0 / (dzk * (2.0 * (zf1[-1] - zh1[-1]) if k == nk - 1 else zh1[k + 1] - zh1[k]))
        for j in range(nj):
            for i in range(ni):
                r = idx(i, j, k)
                a[r, r] -= lo + up + 2.0 / dx ** 2 + (2.0 / dy ** 2 if nj > 1 else 0.0)
                a[r, idx(i + 1, j, k)] += 1.0 / dx ** 2
                a[r, idx(i - 1, j, k)] += 1.0 / dx ** 2
                if nj > 1:
                    a[r, idx(i, j + 1, k)] += 1.0 / dy ** 2
                    a[r, idx(i, j - 1, k)] += 1.0 / dy ** 2
                if k == 0:
                    a[r, r] -= lo                                         # ghost below the ground holds -phi
                else:
                    a[r, idx(i, j, k - 1)] += lo
                if k == nk - 1:
                    a[r, r] -= up                                         # ghost above the top holds -phi
                else:
                    a[r, idx(i, j, k + 1)] += up
    rhs = -np.array([q[i, j, k] for k in range(nk) for j in range(nj) for i in range(ni)]) / EPS0
    sol = np.linalg.solve(a, rhs)
    return np.array(sol).reshape(nk, nj, ni).transpose(2, 1, 0)


def test_the_potential_solves_the_difference_equations_on_stretched_levels(program):
    rng = np.random.default_rng(3)
    zh1, zf1 = stretched(10, 200.0, 1.15, 20000.0)
    q = rng.normal(0.0, 1.0e-9, (8, 6, 10))
    out = run(program, q, zh1, zf1, 2000.0, 3000.0)
    ref = direct_potential(q, zh1, zf1, 2000.0, 3000.0)
    assert np.abs(out['phi'] - ref).max() < 1.0e-5 * np.abs(ref).max()
    # the field is minus the centred gradient, with the mirrored ghosts at the ground and the top
    ext = np.concatenate([-ref[:, :, :1], ref, -ref[:, :, -1:]], axis=2)
    zext = np.concatenate([[2 * zf1[0] - zh1[0]], zh1, [2 * zf1[-1] - zh1[-1]]])
    ez = -(ext[:, :, 2:] - ext[:, :, :-2]) / (zext[2:] - zext[:-2])
    ex = -(np.roll(ref, -1, 0) - np.roll(ref, 1, 0)) / 4000.0
    ey = -(np.roll(ref, -1, 1) - np.roll(ref, 1, 1)) / 6000.0
    scale = np.abs(ez).max()
    assert np.abs(out['ez'] - ez).max() < 1.0e-4 * scale
    assert np.abs(out['ex'] - ex).max() < 1.0e-4 * scale and np.abs(out['ey'] - ey).max() < 1.0e-4 * scale
    assert out['energy'] == pytest.approx(0.5 * np.sum(q * ref * 2000.0 * 3000.0 * np.diff(zf1)), rel=1.0e-4)


def test_a_charged_layer_between_grounded_plates_has_its_exact_potential(program):
    nk, top = 200, 10000.0
    zf1 = np.linspace(0.0, top, nk + 1)
    zh1 = 0.5 * (zf1[:-1] + zf1[1:])
    q0, z1, z2 = 1.0e-9, 4000.0, 6000.0
    q = np.where((zh1 > z1) & (zh1 < z2), q0, 0.0)[None, None, :] * np.ones((4, 4, 1))
    out = run(program, q, zh1, zf1, 1000.0, 1000.0)
    # phi'' = -q/eps0, phi(0) = phi(top) = 0: the layer's charge q0 (z2 - z1) splits between the plates
    s = q0 * (z2 - z1) / EPS0
    zc = 0.5 * (z1 + z2)
    exact = np.where(zh1 < z1, s * (top - zc) / top * zh1,
                     np.where(zh1 > z2, s * zc / top * (top - zh1),
                              s * (top - zc) / top * zh1 - q0 / EPS0 * (zh1 - z1) ** 2 / 2.0))
    assert np.abs(out['phi'][2, 2] - exact).max() < 1.0e-3 * exact.max()
    assert np.abs(out['ex']).max() < 1.0e-6 * np.abs(out['ez']).max()


def light1d_shares(q, columns, dv, thr=0.1e-9, frac=0.3):
    """WRF-ELEC's light1d rule for the shares of charge left (x negative, y positive), written out again."""
    inside = q[columns][:, :-1]
    v = dv[None, :-1]
    chgpos = np.sum(np.where(inside > thr, inside - thr, 0.0) * v)
    chgneg = np.sum(np.where(inside < -thr, inside + thr, 0.0) * v)
    totpos = np.sum(np.maximum(inside, 0.0) * v)
    totneg = np.sum(np.minimum(inside, 0.0) * v)
    chgpos, chgneg = max(chgpos, 1.0), min(chgneg, -1.0)
    fn, fp = frac * 0.5 * (abs(chgneg) + abs(totneg)), frac * 0.5 * (chgpos + totpos)
    if np.sign(totneg + totpos) != np.sign(chgneg + totneg + chgpos + totpos):
        fn = fp = min(fn, fp)
    x = (chgneg + fn) / chgneg if abs(totneg) > 0 else 0.0
    y = (chgpos - fp) / abs(chgpos) if abs(totpos) > 0 else 0.0
    if x < 0.1 or y < 0.1:
        chg = min(frac * max(chgpos, abs(chgneg)), chgpos, abs(chgneg))
        x, y = (abs(chgneg) - chg) / abs(chgneg), (chgpos - chg) / abs(chgpos)
    return x, y


def dipole(q, i, j, zh1, strength=5.0e-9):          # 5 nC/m3 in one 2-km column: about 275 kV/m
    q[i, j] += strength * (np.exp(-((zh1 - 9000.0) / 1000.0) ** 2) - np.exp(-((zh1 - 6000.0) / 1000.0) ** 2))


def test_a_discharge_removes_charge_as_light1d_does_and_counts_separate_regions(program):
    nk, dx = 30, 2000.0
    zf1 = np.linspace(0.0, 15000.0, nk + 1)
    zh1 = 0.5 * (zf1[:-1] + zf1[1:])
    q = np.zeros((24, 20, nk))
    dipole(q, 5, 5, zh1)
    dipole(q, 16, 14, zh1, 1.0e-10)                                      # weak, far away: below breakdown
    out = run(program, q, zh1, zf1, dx, dx, radius=4000.0)
    assert out['ninit'] > 0 and out['nregion'] == 1
    changed = np.any(out['dep'] != 0.0, axis=2)
    assert changed[5, 5] and not changed[16, 14]
    assert np.allclose(out['q'], q + out['dep'], atol=1.0e-15)
    above = out['emag'] > out['ebrk']
    near = np.zeros(changed.shape, bool)
    for i, j in zip(*np.nonzero(above.any(axis=2))):
        for di in range(-2, 3):
            for dj in range(-2, 3):
                if (di * dx) ** 2 + (dj * dx) ** 2 <= 4000.0 ** 2:
                    near[(i + di) % 24, (j + dj) % 20] = True
    assert not np.any(changed & ~near)                                    # only columns within the radius
    columns = near & (np.abs(q[:, :, :-1]).max(axis=2) > 0.1e-9)
    x, y = light1d_shares(q, columns, dx * dx * np.diff(zf1))
    dv = dx * dx * np.diff(zf1)
    qn = q[columns]
    expect_neg = np.sum(np.where(qn < -0.1e-9, -(qn + 0.1e-9) * (1 - x), 0.0) * dv)
    expect_pos = np.sum(np.where(qn > 0.1e-9, (qn - 0.1e-9) * (1 - y), 0.0) * dv)
    assert out['qneg'] == pytest.approx(expect_neg, rel=1.0e-4)
    assert out['qpos'] == pytest.approx(expect_pos, rel=1.0e-4)
    # two strong storms make two regions; one across the domain's edge stays one
    dipole(q, 16, 14, zh1)
    assert run(program, q, zh1, zf1, dx, dx, radius=4000.0)['nregion'] == 2
    q = np.zeros((24, 20, nk))
    dipole(q, 0, 10, zh1)
    dipole(q, 23, 10, zh1)
    out = run(program, q, zh1, zf1, dx, dx, radius=4000.0)
    assert out['nregion'] == 1 and np.any(out['dep'][23, 10] != 0.0) and np.any(out['dep'][0, 10] != 0.0)
