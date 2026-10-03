"""Check coupling conventions that materially change twilight interpretation."""
import json

import numpy as np
import pytest

from illumination.cloud_light.model import MolecularColumn, atlas_proxy, photopic_weights
from research.studies.cloud_twilight.run import HERE, ROOT, digest, path_row, read_clouds
from shared.constants import MOON_RADIUS


def test_corrected_cloud_product_and_occurrence():
    summary, section, context = read_clouds()
    assert summary['gcm']['run'] == 'A28_dim5_moon'
    assert summary['storms']['cloud_top_km']['median'] == pytest.approx(22.0693135)
    assert context['peak_rain_storm_hour_angle_deg'] < 90
    assert np.asarray(section['fraction']).shape == (36, 81)


def test_sunward_offset_and_scattering_convention():
    vacuum = MolecularColumn(MOON_RADIUS, [0., 100000.], np.zeros((1, 3)), [450., 550., 650.])
    r = path_row(vacuum, np.ones(3), 40., 100., 5.)
    assert r['cloud_local_sun_elevation_deg'] > -5
    assert r['cloud_hour_angle_deg'] < 95
    assert r['scattering_angle_deg'] == pytest.approx(r['apparent_elevation_deg'] + 5)
    assert r['two_path_source_fraction'] == pytest.approx(1.)
    assert path_row(vacuum, np.ones(3), 40., 400., 5.)['line_of_sight_blocked']
    assert path_row(vacuum, np.ones(3), 40., 400., 5.)['two_path_source_fraction'] == 0


def test_two_paths_are_combined_spectrally():
    r = path_row(atlas_proxy(), photopic_weights(), 40., 100., 5.)
    product_of_means = r['cloud_solar_transmission'] * r['observer_path_source_weighted_transmission']
    assert r['two_path_source_fraction'] > product_of_means * 1.05
    assert r['two_path_source_fraction'] < min(r['cloud_solar_transmission'],
                                             r['observer_path_source_weighted_transmission'])


def test_committed_result_is_current_and_bounded():
    product = json.loads((HERE / 'results/cloud_twilight.json').read_text())
    for path, expected in product['producer']['files'].items():
        assert digest(ROOT / path) == expected, path
    for path, record in product['producer']['inputs'].items():
        assert digest(ROOT / path) == record['sha256'], path
    for row in product['overhead'] + product['views']:
        assert 0 <= row['two_path_source_fraction'] <= min(
            row['cloud_solar_transmission'], row['observer_path_source_weighted_transmission']) + 1e-14
        if row['line_of_sight_blocked']:
            assert row['two_path_source_fraction'] == 0
