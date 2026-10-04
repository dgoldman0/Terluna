"""Checks of the screening layers the electrified CM1 gathers at cloud edges (climate/crm/fortran/terluna_screen.F, after
WRF-ELEC's screen), compiled and run on a block of cloud in a uniform field. Skipped where gfortran is missing."""
import math
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parents[1] / 'crm' / 'fortran'
EPS = 8.8592e-12

DRIVER = """
program check
  use terluna_screen
  implicit none
  integer :: ni, nj, nk, u, k
  real :: dx, dy, dt
  real, allocatable :: zh1(:), zf1(:), phi(:,:,:), cloud(:,:,:), q(:,:,:), sa(:), sc(:), dep(:,:,:)
  double precision :: qpos, qneg
  open(newunit=u, file='in.bin', access='stream', form='unformatted', status='old')
  read(u) ni, nj, nk, dx, dy, dt
  allocate( zh1(nk), zf1(nk+1), phi(ni,nj,nk), cloud(ni,nj,nk), q(ni,nj,nk), sa(nk), sc(nk), dep(ni,nj,nk) )
  read(u) zh1, zf1, phi, cloud, q
  close(u)
  do k = 1,nk
    sa(k) = earth_air_conductivity(zh1(k))
    sc(k) = 0.1*sa(k)
  enddo
  call screen(ni, nj, nk, dx, dy, zh1, zf1, dt, phi, cloud, q, sa, sc, dep, qpos, qneg)
  open(newunit=u, file='out.bin', access='stream', form='unformatted', status='replace')
  write(u) dep, sa, qpos, qneg
  close(u)
end program check
"""


@pytest.fixture(scope='module')
def program(tmp_path_factory):
    if shutil.which('gfortran') is None:
        pytest.skip('gfortran is missing')
    folder = tmp_path_factory.mktemp('screen')
    (folder / 'check.f90').write_text(DRIVER)
    subprocess.run(['gfortran', '-O2', '-fopenmp', '-ffree-form', '-ffree-line-length-none', '-cpp',
                    str(HERE / 'terluna_screen.F'), 'check.f90', '-o', 'check'], cwd=folder, check=True,
                   capture_output=True)
    return folder


NI, NJ, NK, DX, DZ, DT = 16, 16, 30, 1000.0, 500.0, 6.0
ZF = np.arange(NK + 1) * DZ
ZH = 0.5 * (ZF[1:] + ZF[:-1])


def run(folder, phi, cloud, q):
    with open(folder / 'in.bin', 'wb') as out:
        np.array([NI, NJ, NK], np.int32).tofile(out)
        np.array([DX, DX, DT], np.float32).tofile(out)
        for a in (ZH, ZF, phi, cloud, q):
            np.asarray(a, np.float32).ravel(order='F').tofile(out)
    subprocess.run([str(folder / 'check')], cwd=folder, check=True, capture_output=True,
                   env=dict(os.environ, OMP_NUM_THREADS='2'))
    raw = (folder / 'out.bin').read_bytes()
    n = NI * NJ * NK
    dep = np.frombuffer(raw[:4 * n], np.float32).reshape(NK, NJ, NI).transpose(2, 1, 0)
    sa = np.frombuffer(raw[4 * n:4 * n + 4 * NK], np.float32)
    qpos, qneg = np.frombuffer(raw[4 * n + 4 * NK:], np.float64)
    return dep, sa, qpos, qneg


def block(i0, i1, k0, k1, j0=5, j1=12):
    cloud = np.zeros((NI, NJ, NK))
    for i in range(i0, i1 + 1):
        cloud[i % NI, j0:j1 + 1, k0:k1 + 1] = 1.0e-3
    return cloud


def expected(e, sa):
    """The charge an edge gathers in a step when the field along its outward normal is e inside and out."""
    sc = 0.1 * sa
    ds = (DX * DX * DZ) ** (1.0 / 3.0)
    return -2.0 * EPS / ds * e * (sa - sc) / (sa + sc) * (1.0 - math.exp(-DT * (sa + sc) / (2.0 * EPS)))


def test_earth_s_conductivity_table_and_its_slope_beyond(program):
    _, sa, _, _ = run(program, np.zeros((NI, NJ, NK)), np.zeros((NI, NJ, NK)), np.zeros((NI, NJ, NK)))
    assert sa[0] == pytest.approx(5.5e-14 + (6.0e-14 - 5.5e-14) * 0.5, rel=1e-5)        # 250 m
    assert sa[19] == pytest.approx(4.5e-13 + 0.4e-13 * 0.5, rel=1e-5)                    # 9750 m
    assert sa[23] == pytest.approx(4.5e-13 + 0.4e-13 * 11750.0 / 500.0 - 0.4e-13 * 19, rel=1e-5)   # 11750 m


def test_an_upward_field_charges_the_cloud_top_negative_and_its_base_positive(program):
    cloud = block(5, 12, 7, 19)                                  # 3.5 to 10 km
    phi = -100.0 * ZH[None, None, :] * np.ones((NI, NJ, 1))      # 100 V/m upward
    dep, sa, qpos, qneg = run(program, phi, cloud, np.zeros_like(cloud))
    top, base = dep[6:12, 6:12, 19], dep[6:12, 6:12, 7]        # away from the block's sides
    assert np.allclose(top, expected(100.0, sa[19]), rtol=1e-4) and np.all(top < 0.0)
    assert np.allclose(base, expected(-100.0, sa[7]), rtol=1e-4) and np.all(base > 0.0)
    inside = dep[6:12, 6:12, 9:18]
    assert not np.any(inside) and not np.any(dep[cloud == 0.0])
    sides = dep[5, 6:12, 9:18]                                   # an edge along x, the field along it
    assert not np.any(sides)
    assert qpos > 0.0 and qneg < 0.0


def test_the_charge_there_stays_within_a_quarter_nc_per_cubic_metre(program):
    cloud = block(5, 12, 7, 19)
    phi = -100.0 * ZH[None, None, :] * np.ones((NI, NJ, 1))
    q = np.zeros_like(cloud)
    q[8, 8, 19] = -0.25e-9 + 1.0e-13                             # top: room for 1e-13 C/m3 more
    q[8, 8, 7] = 0.3e-9                                          # base: already beyond the cap
    dep, *_ = run(program, phi, cloud, q)
    assert dep[8, 8, 19] == pytest.approx(-1.0e-13, rel=1e-3) and dep[8, 8, 7] == 0.0
    assert dep[9, 9, 7] > 0.0


def test_edges_along_x_gather_charge_across_the_domain_s_wrap(program):
    cloud = block(13, 19, 2, 27, j0=0, j1=NJ - 1)               # columns 13-15 and 0-3, every row, 1 to 13.5 km
    x = (np.arange(NI) + 0.5) * DX
    phi = -50.0 * x[:, None, None] * np.ones((1, NJ, NK))       # 50 V/m toward +x
    dep, sa, *_ = run(program, phi, cloud, np.zeros_like(cloud))
    # at the +x edge (column 3) the outward normal runs with the field, at the -x edge (column 13) against it
    assert np.allclose(dep[3, :, 10], expected(50.0, sa[10]), rtol=1e-4)
    assert np.allclose(dep[13, :, 10], expected(-50.0, sa[10]), rtol=1e-4)
    assert not np.any(dep[[14, 15, 0, 1, 2], :, 10])
