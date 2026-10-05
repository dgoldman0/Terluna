"""Independent geometry and physical-unit checks for the supplemental view."""
import importlib.util
import json
import math
import os
from pathlib import Path

import numpy as np
import pytest


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def module(name):
    spec = importlib.util.spec_from_file_location(name, HERE / (name + '.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_portrait_keeps_earths_angular_size_and_elevation():
    from research.studies.sea_appearance import scenes

    data = json.loads((ROOT / 'research/studies/sea_appearance/results/nubium-midnight.json').read_text())
    scene = data['scene']
    camera = module('nubium_guides').Camera(scene['camera'])
    earth = scene['earth']
    centre = scenes.unit(earth['elevation_deg'], earth['azimuth_deg'])
    x, y = camera.project(centre)
    # The description rounds azimuth to 0.1 degree before this reconstruction.
    assert float(x[0]) == pytest.approx(earth['x'], abs=0.7)
    assert float(y[0]) == pytest.approx(earth['y'], abs=0.2)
    # Construct two limb rays in angular space, then project and recover them.
    tangent = camera.r - np.dot(camera.r, centre) * centre
    tangent /= np.linalg.norm(tangent)
    radius = math.radians(earth['diameter_deg'] / 2)
    limb = np.stack([centre * math.cos(radius) + sign * tangent * math.sin(radius)
                     for sign in [-1, 1]])
    px, py = camera.project(limb)
    recovered = camera.direction(px, py)
    separation = math.degrees(math.acos(np.dot(*recovered)))
    assert separation == pytest.approx(earth['diameter_deg'], abs=1e-10)
    assert px[1] - px[0] == pytest.approx(earth['width_px'], abs=0.1)


def test_midnight_is_minimum_sun_not_minimum_ground_light():
    folder = ROOT / 'research/studies/sea_appearance/results'
    product = json.loads((folder / 'nubium-midnight.json').read_text())
    calendar = json.loads((folder / 'lighting_calendar.json').read_text())
    station = next(s for s in calendar['sites'] if s['short'] == 'S Nubium')
    series = station['series']
    hours = np.array(calendar['time_hours'])
    selected = product['selection']
    assert selected['hour'] == hours[np.argmin(series['sun_elevation_deg'])]
    light = np.array(series['ground_lux_from_sun']) + series['ground_lux_from_earth']
    assert selected['hour'] != hours[np.argmin(light)]
    assert product['scene']['earth']['lit_fraction'] > .98
    assert product['scene']['sun']['elevation_deg'] < -60


@pytest.mark.skipif(not os.environ.get('TERLUNA_RADIANCE_BIN'), reason='Optional official Radiance executables not configured')
def test_human_operator_keeps_absolute_xyz_units(tmp_path, monkeypatch):
    bins = Path(os.environ['TERLUNA_RADIANCE_BIN'])
    monkeypatch.syspath_prepend(str(HERE))
    condition = module('nubium_human_view').condition
    y = np.geomspace(1e-4, 1e4, 128)[None, :] * np.linspace(.5, 1.5, 96)[:, None]
    xyz = y[..., None] * np.array([.95047, 1., 1.08883])
    record = condition(xyz, {'horizontal_fov_deg': 65, 'vertical_fov_deg': 50},
                       tmp_path / 'luminance-ramp', bins)
    assert record['input_xyz_luminance_roundtrip_max_relative'] < .01
    assert record['input_xyz_roundtrip_error_per_largest_component'] <= 1/256 + 1e-7
    assert record['display_peak_cd_m2'] == 100
    assert (tmp_path / 'luminance-ramp.png').is_file()
    curve = np.loadtxt(tmp_path / 'luminance-ramp.curve.txt')
    assert np.all(np.diff(curve[:, 0]) > 0)
    assert np.all(np.diff(curve[:, 1]) >= -1e-6)
    assert curve[:, 1].min() >= .099
    assert curve[:, 1].max() <= 100.01
