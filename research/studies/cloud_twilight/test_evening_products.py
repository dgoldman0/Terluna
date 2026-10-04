"""Evidence and input contracts for the regional cloud appearance products."""
import hashlib
import json
from pathlib import Path
import numpy as np
import pytest
from shared.constants import SYNODIC_MONTH_DAYS

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(name):
    return json.loads((HERE/'results'/name).read_text())


def test_regional_sample_contract_excludes_first_cycle_and_retains_denominators():
    p=read('regional_evenings.json')
    assert p['schema']=='terluna.research.regional-cloud-evenings/1'
    assert len(p['climate_inputs'])==4
    for c in p['climate_inputs']:
        assert SYNODIC_MONTH_DAYS<=c['span_days'][0]<=c['span_days'][1]<=2*SYNODIC_MONTH_DAYS
        assert c['snapshots']==236
    for r in p['regional_bins']:
        assert r['samples']>0
        assert 0<=r['cloud_fraction']<=1
    assert digest(ROOT/p['producer']['file'])==p['producer']['sha256']
    for file,expected in p['producer']['inputs'].items():
        path=ROOT/file
        if path.exists():
            assert digest(path)==expected


def test_geometric_cloud_limits_and_numerical_twilight_threshold_order():
    pytest.importorskip('numba')
    from illumination.cloud_light.evening_context import last_sunlit_depression
    p=json.loads((ROOT/'illumination/cloud_light/results/evening_context.json').read_text())
    thresholds=p['thresholds']
    assert np.all(np.diff([r['level_lux'] for r in thresholds])<0)
    assert np.all(np.diff([r['sun_depression_deg'] for r in thresholds])>0)
    assert last_sunlit_depression(80000,5)>last_sunlit_depression(40000,5)
    assert last_sunlit_depression(40000,0)>last_sunlit_depression(40000,5)


@pytest.mark.parametrize('mode',['selected','history','sensitivity','deep'])
def test_scene_sources_and_completed_transport(mode):
    p=read(f'evening_scenes_{mode}.json')
    assert p['scenes']
    for scene in p['scenes']:
        assert SYNODIC_MONTH_DAYS<=scene['time_day']<=2*SYNODIC_MONTH_DAYS
        assert scene['sampled_truncated_paths']==0
        assert scene['negative_Y_pixels']==0
        assert scene['surface_lux']>=0
        assert scene['recess_lux']>=0
        for file,expected in scene['producer'].items():
            assert digest(ROOT/file)==expected
        archive=ROOT/scene['archive']['path']
        if archive.exists():
            assert digest(archive)==scene['archive']['sha256']
            with np.load(archive) as f:
                np.testing.assert_allclose(f['xyz'],f['group_xyz'].mean(axis=2),rtol=1e-12,atol=1e-10)


def test_history_preserves_one_observer_and_view_while_weather_changes():
    scenes=read('evening_scenes_history.json')['scenes']
    assert len(scenes)==7
    assert len({s['observer_longitude_deg'] for s in scenes})==1
    assert len({s['observer_latitude_deg'] for s in scenes})==1
    assert len({s['camera_azimuth_in_ring_basis_deg'] for s in scenes})==1
    assert len({s['source']['cloud_sha256'] for s in scenes})==7
    assert np.all(np.diff([s['snapshot'] for s in scenes])>0)
    assert np.all(np.diff([s['observer_sun_elevation_deg'] for s in scenes])<0)
