"""Stored evidence and cross-lane contracts for the consistent sky study."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent


def test_stored_spherical_products_match_their_producers_and_inputs():
    sky = json.loads((ROOT/'illumination/sky/results/solved_sky.json').read_text())
    study = json.loads((HERE/'results/spherical_sky.json').read_text())
    assert sky['schema']=='terluna.illumination.solved-spherical-sky/1'
    assert study['schema']=='terluna.research.optical-comfort-spherical-sky/1'
    producers = [sky['producer'],study['producer']]
    producers += [world['optical_column']['producer'] for world in sky['worlds'].values()]
    for producer in producers:
        for group in ('files','inputs'):
            for name,expected in producer.get(group,{}).items():
                assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected,name


def test_stored_sky_convergence_and_independent_checks_are_explicit():
    sky = json.loads((ROOT/'illumination/sky/results/solved_sky.json').read_text())
    for world in sky['worlds'].values():
        assert world['settings']['order_converged']
        assert world['settings']['albedo']==.1
        assert world['spectral_check']['maximum_absolute_relative_lux_error']<.015
        for row in world['samples']:
            assert row['total_horizontal_lux']>=0
            assert np.isclose(row['total_horizontal_lux'],row['direct_horizontal_lux']+row['diffuse_horizontal_lux'])
    checks = sky['additional_checks']['numerical_checks.json']
    assert checks['all_monte_carlo_paths_completed']
    for row in checks['spectral_monte_carlo']:
        if row['mean_diffuse_lux']==0:
            assert row['relative_difference'] is None
            assert row['percent_comparison_resolved'] is False


def test_every_existing_gaze_and_scene_received_the_consistent_sky():
    new = json.loads((HERE/'results/spherical_sky.json').read_text())
    old = json.loads((HERE/'results/native_sky.json').read_text())
    assert len(new['angular_views'])==len(old['angular_views'])==630
    assert len(new['scenes'])==len(old['scenes'])==24
    assert all(r['flux_quality']=='solved_column_spherical_multiple_scattering' for r in new['angular_views'])
    for a,b in zip(new['angular_views'],old['angular_views'],strict=True):
        assert all(a[k]==b[k] for k in ('world','sun_deg','gaze_deg','azimuth_from_sun_deg'))


def test_available_atlases_supply_the_same_surface_flux_as_the_study():
    from .angular import AngularSky
    sky = json.loads((ROOT/'illumination/sky/results/solved_sky.json').read_text())
    if any(not (ROOT/w['atlas']['path']).is_file() for w in sky['worlds'].values()):
        pytest.skip('Local spherical sky archives are required for the angular-flux check')
    for world in sky['worlds'].values():
        with np.load(ROOT/world['atlas']['path'],allow_pickle=False) as atlas:
            for sun in (1.,10.,90.):
                i = np.flatnonzero(atlas['suns']==sun)[0]
                angular = AngularSky(np.sin(np.radians(atlas['elevations'])),np.radians(atlas['azimuths']),atlas['xyz'][i,:,:,1])
                row = next(r for r in world['samples'] if r['sun_deg']==sun)
                assert angular.horizontal==pytest.approx(row['diffuse_horizontal_lux'],rel=1e-13)
