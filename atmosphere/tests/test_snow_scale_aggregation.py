"""Controls for the external-solver adapter; no compilation or MPI in pytest."""
import hashlib
from pathlib import Path
import struct

import numpy as np
import pytest

from atmosphere.snow import scale_aggregation as a
from atmosphere.snow import scale_collision_observer as observer


@pytest.mark.parametrize('kind', ['ice','golovin'])
def test_seeded_population_has_declared_water_in_every_collision_cell(kind):
    p=a.initial_population(128,kind,42)
    n=p['n'].reshape(a.NCELL,128)
    m=p['mass'].reshape(a.NCELL,128)
    water=a.ICE_WATER_CONTENT if kind=='ice' else a.LIQUID_WATER_CONTENT
    np.testing.assert_allclose((n*m).sum(axis=1)/a.CELL_VOLUME,water,rtol=2e-15)
    if kind=='ice':
        np.testing.assert_allclose(p['mass'],4*np.pi/3*p['re']**2*p['rp']*p['rho'],rtol=3e-15)
    # Spatial seed cannot silently change the prescribed size distribution.
    q=a.initial_population(128,kind,43)
    np.testing.assert_array_equal(q['mass'],p['mass'])
    assert not np.array_equal(q['z'],p['z'])
    assert np.all((p['z']>0)&(p['z']<a.GRID[0]*a.DX))


def test_native_restart_layout_has_correct_integer_widths_and_record_boundaries(tmp_path):
    from scipy.io import FortranFile, FortranEOFError
    p=a.initial_population(4,'ice',0)
    path=tmp_path/'particles'
    a.write_restart(path,p)
    with path.open('rb') as f:
        assert struct.unpack('>i',f.read(4))[0]==28
        assert struct.unpack('>d5i',f.read(28))==(0.,8,8,128,1,1)
    with FortranFile(path,'r',header_dtype='>u4') as f:
        f.read_record('u1')
        np.testing.assert_array_equal(f.read_record('>i8'),p['n'])
        for key in ('rk','x','y','z','r'):
            np.testing.assert_array_equal(f.read_record('>f8'),p[key])
        for _ in range(4): assert np.all(f.read_record('>f8')==0.)
        np.testing.assert_array_equal(f.read_record('>i2'),p['phase'])
        for key in ('re','rp','rho'):
            np.testing.assert_array_equal(f.read_record('>f8'),p[key])
        f.read_record('>f8');f.read_record('>f8')
        np.testing.assert_array_equal(f.read_record('>i4'),p['monomers'])
        with pytest.raises(FortranEOFError): f.read_record('u1')
    with pytest.raises(FileExistsError): a.write_restart(path,p)


def test_golovin_distribution_normalization_moments_and_quadrature_refinement():
    # Integrate the mass measure independently to recover N and M2. This also
    # catches mistakes in mass/volume units, Bessel scaling, and tail truncation.
    number=a.LIQUID_NUMBER;water=a.LIQUID_WATER_CONTENT;mean=water/number
    mass=mean*np.geomspace(1e-7,1e8,60001)
    middle=np.sqrt(mass[:-1]*mass[1:])
    for t in (0.,1200.,3600.):
        cdf,norm=a.golovin_mass_cdf(mass,number,water,t)
        refined,_=a.golovin_mass_cdf(mass,number,water,t,points=48000)
        assert abs(norm-1)<2e-6
        assert np.max(np.abs(cdf-refined))<2e-6
        assert np.all(np.diff(cdf)>=-1e-15)
        weight=np.diff(cdf)
        inferred=[water*np.sum(weight/middle),water,water*np.sum(weight*middle)]
        exact=a.golovin_moments([number,water,2*number*mean**2],t)
        np.testing.assert_allclose(inferred,exact,rtol=2e-5)


def test_cdf_distance_includes_both_sides_of_atoms():
    # A single atom versus a continuous CDF which reaches 0.4 at the atom:
    # the larger error is immediately after the jump, 1 - 0.4 = 0.6.
    assert a.distribution_distance(np.array([2.,2.]),np.array([1.,3.]),
                                   lambda x:np.full_like(x,.4))==pytest.approx(.6)


@pytest.mark.parametrize('kwargs',[{'dt':0.},{'duration':-300.},{'dt':float('nan')},
                                  {'sd':0},{'kind':'other'},{'lunar':True,'kind':'golovin'}])
def test_invalid_configuration_does_not_create_a_run(tmp_path,kwargs):
    with pytest.raises(ValueError): a.configure(tmp_path,'invalid',**kwargs)
    assert not list(tmp_path.iterdir())


def test_observer_rejects_unknown_source():
    with pytest.raises(ValueError,match='Unidentified'):
        observer.patch('unknown implementation')


def test_observer_patch_preserves_every_original_source_line_when_cache_available():
    path=a.CACHE/'work'/observer.SOURCE_REL
    if not path.exists(): pytest.skip('Optional pinned native source not cached')
    original=path.read_text()
    assert hashlib.sha256(original.encode()).hexdigest()==observer.SOURCE_SHA256
    modified=observer.patch(original)
    # All changes are additions. Exact native output parity is a separate run
    # check, not inferred from this textual control.
    from difflib import SequenceMatcher
    operations=SequenceMatcher(None,original.splitlines(),modified.splitlines()).get_opcodes()
    assert all(tag in ('equal','insert') for tag,*_ in operations)
