"""Unit controls for source identity, water accounting and particle weighting."""
import json
import tarfile
from pathlib import Path

import numpy as np
import pytest

from atmosphere.snow import scale_build, scale_column
from shared.constants import MOON_SURFACE_GRAVITY


def test_water_ledger_includes_changed_gas_density_and_accumulated_precipitation():
    # Two 2 m tall boxes, per square metre. Vapor 0.20 -> 0.10 kg,
    # suspended water 0.04 -> 0.08 kg, fallen water 0 -> 0.06 kg.
    rho=np.array([1.,.95]).reshape(2,1,1,1)
    qv=np.array([.1,.05/.95]).reshape(2,1,1,1)
    qh=np.array([.02,.04/.95]).reshape(2,1,1,1)
    accumulation=np.array([0.,.00006]).reshape(2,1,1)
    vapor,hydro,precip,total=scale_column.water_history(rho,qv,qh,accumulation,2.,1.)
    np.testing.assert_allclose(vapor,[.2,.1],atol=1e-15)
    np.testing.assert_allclose(hydro,[.04,.08],atol=1e-15)
    np.testing.assert_allclose(precip,[0.,.06],atol=1e-15)
    np.testing.assert_allclose(total,[.24,.24],atol=1e-15)
    assert abs((vapor+hydro)[1]-(vapor+hydro)[0])==pytest.approx(.06)


@pytest.mark.parametrize('bad',[-.1,float('nan'),float('inf')])
def test_bad_native_water_is_an_error(bad):
    rho=np.ones((2,1,1,1));qv=np.full_like(rho,.1);qh=np.zeros_like(rho)
    qv[-1]=bad
    with pytest.raises(ValueError):
        scale_column.water_history(rho,qv,qh,np.zeros((2,1,1)),1.,1.)


def test_particle_quantiles_use_physical_multiplicity_and_distinguish_mass():
    size=np.array([1.,10.]);multiplicity=np.array([99.,1.])
    assert scale_column.weighted_quantile(size,multiplicity,[.5])==[1.]
    assert scale_column.weighted_quantile(size,multiplicity*size**3,[.5])==[10.]
    assert scale_column.weighted_quantile(np.array([]),np.array([]))==[None]*4
    with pytest.raises(ValueError):
        scale_column.weighted_quantile(size,np.array([-1.,2.]))


def test_lunar_configuration_uses_shared_gravity_and_keeps_closed_water_supply(tmp_path):
    folder=scale_column.configure(tmp_path,'lunar',lunar=True,move=True)
    settings=json.loads((folder/'settings.json').read_text())
    assert settings['gravity_m_s2']==MOON_SURFACE_GRAVITY
    text=(folder/'run.conf').read_text()
    assert "ATMOS_DYN_TYPE = 'OFF'" in text
    assert 'HEAT_RATE = 0.000000000000000d+00' in text
    assert 'sdm_aslset = -1' in text
    assert 'HISTORY_DEFAULT_DATATYPE = \'REAL8\'' in text
    before=(folder/'run.conf').read_bytes()
    with pytest.raises(FileExistsError):
        scale_column.configure(tmp_path,'lunar',lunar=True)
    assert (folder/'run.conf').read_bytes()==before


@pytest.mark.parametrize('name,kwargs',[('../main',{}),('bad',{'dt':0.}),('bad',{'sd':2.5}),('bad',{'dt':.7})])
def test_bad_case_fails_before_creating_output(tmp_path,name,kwargs):
    with pytest.raises(ValueError):
        scale_column.configure(tmp_path,name,**kwargs)
    assert not list(tmp_path.iterdir())


def test_unidentified_upstream_source_is_never_patched():
    for name in scale_build.PATCH_IDENTITIES:
        with pytest.raises(ValueError,match='identity mismatch'):
            scale_build.patch_text(name,'unidentified source\n')


def test_wrong_archive_does_not_extract_or_build(tmp_path):
    archive=tmp_path/'bad.tar.gz';archive.write_bytes(b'wrong archive')
    with pytest.raises(ValueError,match='SHA256 mismatch'):
        scale_build.prepare(tmp_path/'cache',archive)
    assert not (tmp_path/'cache/work').exists()


def test_reviewed_patches_reproduce_actual_pinned_source_when_cache_available():
    # No network, compile or numerical run in the routine test suite.
    archive=Path('/tmp/terluna-scale-feasibility/source.tar.gz')
    if not archive.exists():
        pytest.skip('Optional pinned upstream archive is not cached')
    assert scale_build.sha256(archive)==scale_build.ARCHIVE_SHA256
    with tarfile.open(archive) as source:
        prefix=f'SCALE-SDM_mixed-phase_Shima2019-{scale_build.REVISION}/'
        for name,(_,expected) in scale_build.PATCH_IDENTITIES.items():
            text=source.extractfile(prefix+name).read().decode()
            patched=scale_build.patch_text(name,text)
            import hashlib
            assert hashlib.sha256(patched.encode()).hexdigest()==expected
