"""The Jura Island flyover keeps its camera and panels, and reads its water level and terrain from the products.

    python3 -m pytest visualization/jura-island/tests

None of these tests needs OpenGL. The tests that build the terrain skip when geography's LOLA and GRAIL inputs are
not fetched (python -m geography.fetch_inputs --download). The original film's numbers come from its record,
repository_import.json.
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pytest

TOOL = Path(__file__).resolve().parents[1]
ROOT = TOOL.parents[1]
sys.path[:0] = [str(TOOL), str(ROOT)]
import film  # noqa: E402
import terrain  # noqa: E402
from geography import fetch_inputs  # noqa: E402

ORIGINAL = json.loads((TOOL / 'repository_import.json').read_text())['source_record']
FETCHED = all((fetch_inputs.input_dir() / name).is_file() for name in terrain.INPUTS)
needs_inputs = pytest.mark.skipif(not FETCHED, reason='LOLA and GRAIL inputs are not fetched')


@pytest.fixture(scope='module')
def ground():
    return terrain.build()


def products_copy(tmp_path):
    island, atlas = json.loads(terrain.ISLAND.read_text()), json.loads(terrain.ATLAS.read_text())
    return island, atlas, tmp_path / 'shores.json', tmp_path / 'atlas.json'


# The camera

def test_camera_path_is_finite_and_keeps_the_original_bounds():
    checks, original = film.camera_checks(), ORIGINAL['checks']
    assert checks['frames'] == film.FPS * film.DURATION == original['frames']
    assert checks['finite']
    assert checks['minimum_altitude_km'] == pytest.approx(original['minimum_camera_altitude_above_water_km'], abs=1e-6)
    assert checks['maximum_rotation_deg_per_frame'] == pytest.approx(
        original['maximum_camera_rotation_degrees_per_frame'], rel=1e-4)
    # the camera never looks straight down, where a heading would be undefined
    low, high = checks['pitch_range_deg']
    assert -80 < low < high < -15
    assert checks['maximum_heading_step_deg'] < 0.5


def test_camera_moves_without_jumps_and_starts_and_ends_at_rest():
    def steps(times):
        p = film.path(times)
        eye = film.xyz(p[:, 0], p[:, 1], p[:, 2])
        angles = film.orientation(times)
        return np.linalg.norm(np.diff(eye, axis=0), axis=1), np.abs(np.diff(angles, axis=0)).max(axis=1)

    coarse, fine = film.frame_times(), np.linspace(0, film.DURATION, 10 * film.FPS * film.DURATION + 1)
    move, turn = steps(coarse)
    fine_move, fine_turn = steps(fine)
    # a jump would not shrink with the step; a smooth path's steps shrink tenfold
    assert fine_move.max() < 0.11 * move.max() and fine_turn.max() < 0.11 * turn.max()
    assert np.all(np.abs(np.diff(film.KEY_ANGLES[:, 0])) < np.pi)          # headings unwrapped at the keys
    # the clamped ends: at 0 s and 44 s the camera neither moves nor turns
    dt = 1e-4
    for t0 in (0.0, film.DURATION - dt):
        speed, spin = (x * film.FPS * dt for x in (move.max(), turn.max()))
        a, b = steps(np.array([t0, t0 + dt]))
        assert a[0] < 1e-3 * speed and b[0] < 1e-3 * spin


def test_camera_flies_inside_the_regional_terrain():
    p = film.path(film.frame_times())
    assert film.MESH['west'] + 5 < p[:, 0].min() and p[:, 0].max() < film.MESH['east'] - 5
    assert film.MESH['south'] + 3 < p[:, 1].min() and p[:, 1].max() < film.MESH['north'] - 5
    assert 16.9 < p[:, 2].min() and p[:, 2].max() < 500    # the climb settles on 480 km from about 489


@needs_inputs
def test_camera_stays_well_above_the_ground(ground):
    p = film.path(film.frame_times())
    clearance = p[:, 2] - np.maximum(ground.heights(p[:, 0], p[:, 1]), 0) / 1000
    assert clearance.min() == pytest.approx(ORIGINAL['checks']['minimum_camera_height_above_terrain_km'], abs=1e-3)


# The panels

def test_panels_change_four_times_and_hold_between_fades():
    times = film.frame_times()
    states = [film.chapter_state(t) for t in times]
    index, opacity = np.array([s[0] for s in states]), np.array([s[1] for s in states])
    changes = np.flatnonzero(np.diff(index))
    assert len(changes) == ORIGINAL['checks']['reference_panel_changes'] and np.all(np.diff(index) >= 0)
    assert np.allclose(times[changes + 1], [c[0] for c in film.CHAPTERS[1:]])
    holds = [np.sum((index == k) & (opacity == 1.0)) / film.FPS for k in range(len(film.CHAPTERS))]
    assert min(holds) >= ORIGINAL['checks']['minimum_full_opacity_hold_seconds']
    # each fade is a smoothstep over the fade's length, so no frame steps the panel's opacity by more than its slope
    assert film.FADE_SECONDS == ORIGINAL['fade_seconds']
    assert np.abs(np.diff(opacity)).max() <= 1.5 / film.FADE_SECONDS / film.FPS + 1e-9
    assert np.all(opacity[changes] < 0.01) and np.all(opacity[changes + 1] == 0.0)


def test_panels_show_the_original_sections():
    for (start, end, title, points), section in zip(film.CHAPTERS, ORIGINAL['sections'], strict=True):
        assert (start, end, title, list(points)) == (section['start_seconds'], section['end_seconds'],
                                                     section['title'], section['points'])
        assert len(points) == ORIGINAL['checks']['reference_points_per_section']


def test_fades_follow_a_smoothstep():
    end = film.CHAPTERS[0][1]
    for s in np.linspace(0, 1, 11):
        t = end - s * film.FADE_SECONDS
        assert film.chapter_state(t - 1e-12)[1] == pytest.approx(s * s * (3 - 2 * s), abs=1e-9)


# The inputs

def test_the_island_study_and_geography_pin_the_inputs_the_original_used():
    manifest = {f['name']: f['sha256'] for f in fetch_inputs.manifest()['files']}
    recorded = json.loads(terrain.ISLAND.read_text())['producer']['inputs']
    for name in terrain.STUDY_INPUTS:
        assert recorded[f'geography/inputs/{name}'] == manifest[name]
    for name in terrain.INPUTS:
        assert ORIGINAL['inputs'][name] == manifest[name]


@needs_inputs
def test_changed_input_bytes_stop_the_build(monkeypatch):
    real = fetch_inputs.manifest()
    recorded = json.loads(terrain.ISLAND.read_text())['producer']['inputs']
    assert set(terrain.check_inputs(recorded)) == {f'geography/inputs/{n}' for n in terrain.INPUTS}

    def changed():
        m = json.loads(json.dumps(real))
        next(f for f in m['files'] if f['name'] == 'ldem_4.img')['sha256'] = '0' * 64
        return m

    monkeypatch.setattr(terrain.fetch_inputs, 'manifest', changed)
    with pytest.raises(ValueError, match='ldem_4.img'):
        terrain.check_inputs(recorded)
    monkeypatch.undo()
    with pytest.raises(ValueError, match='island study'):
        terrain.check_inputs(dict(recorded, **{'geography/inputs/ldem_16.img': '0' * 64}))


# The water level and the water share

def test_water_level_and_share_come_from_the_products():
    p = terrain.products()
    island, atlas = json.loads(terrain.ISLAND.read_text()), json.loads(terrain.ATLAS.read_text())
    assert p['sea_level_m'] == island['sea_level_m'] == atlas['sea_level_m']
    assert p['water_share'] == atlas['water_share']
    captions = film.captions(p['sea_level_m'], p['water_share'], p['pixels_per_degree'])
    assert captions['level'] == f'Common water level {island["sea_level_m"]:,.1f} m relative to the geoid'.replace(
        '-', '−')
    assert captions['share'].startswith(f'{round(100 * atlas["water_share"])}%')
    assert 'Common-water-level view' in captions['view'] and 'illustrative lighting' in captions['view']


def test_changed_products_carry_through_to_the_texts(tmp_path):
    island, atlas, island_path, atlas_path = products_copy(tmp_path)
    island['sea_level_m'] = atlas['sea_level_m'] = -1600.0
    atlas['water_share'] = 0.3
    atlas_path.write_text(json.dumps(atlas))
    island['producer']['inputs']['geography/results/atlas.json'] = terrain.digest(atlas_path)
    island_path.write_text(json.dumps(island))
    p = terrain.products(island_path, atlas_path)
    captions = film.captions(p['sea_level_m'], p['water_share'], p['pixels_per_degree'])
    assert captions['level'] == 'Common water level −1,600.0 m relative to the geoid'
    assert captions['share'] == '30% of the Moon under water'


def test_products_that_disagree_stop_the_build(tmp_path):
    island, atlas, island_path, atlas_path = products_copy(tmp_path)
    atlas['water_share'] = 0.3                     # an atlas other than the one the island study measured against
    atlas_path.write_text(json.dumps(atlas))
    island_path.write_text(json.dumps(island))
    with pytest.raises(ValueError, match='atlas differs'):
        terrain.products(island_path, atlas_path)
    island['producer']['inputs']['geography/results/atlas.json'] = terrain.digest(atlas_path)
    island['sea_level_m'] = -1600.0
    island_path.write_text(json.dumps(island))
    with pytest.raises(ValueError, match='different sea levels'):
        terrain.products(island_path, atlas_path)
    island['schema'] = 'terluna.research.jura-shores/2'
    island_path.write_text(json.dumps(island))
    with pytest.raises(ValueError, match='Unsupported island product'):
        terrain.products(island_path, atlas_path)


@needs_inputs
def test_heights_stand_on_the_products_water_level(ground):
    assert ground.sea_level_m == json.loads(terrain.ISLAND.read_text())['sea_level_m']
    # a pixel's height plus the level gives back its height above the geoid, LOLA's radius less the geoid
    raw = np.memmap(fetch_inputs.path('ldem_16.img'), dtype='<i2', mode='r', shape=(2880, 5760))
    i, j = 200, 300
    row, col = round((90 - ground.lat[i]) * 16 - 0.5), round((ground.lon[j] + 360) * 16 - 0.5)
    lat, lon = ground.lat[i:i + 1], ground.lon[j:j + 1] + 360
    above_geoid = raw[row, col] * 0.5 - terrain.geoid(terrain.Grid(lat, lon, 16), max_degree=200)[0, 0]
    assert ground.z[i, j] + ground.sea_level_m == pytest.approx(above_geoid, abs=1e-6)


# The terrain against the island study

@needs_inputs
def test_terrain_measures_as_the_island_study(ground):
    island = json.loads(terrain.ISLAND.read_text())['island']
    assert round(float(ground.cell_area_km2()[ground.island].sum())) == island['area_km2']
    peak = np.unravel_index(np.argmax(np.where(ground.island, ground.z, -np.inf)), ground.z.shape)
    assert round(float(ground.z[peak])) == island['summit_m'] == round(ground.summit[2])
    assert ground.lat[peak[0]] == pytest.approx(island['summit_lat_deg'], abs=1e-3)
    assert ground.lon[peak[1]] + 360 == pytest.approx(island['summit_lon_deg'], abs=1e-3)
    rows, cols = np.nonzero(ground.island)
    assert [round(float(ground.lat[rows].min()), 2), round(float(ground.lat[rows].max()), 2)] == island['lat_range_deg']
    assert [round(float(ground.lon[cols].min()) + 360, 2),
            round(float(ground.lon[cols].max()) + 360, 2)] == island['lon_range_deg']
    edges = np.concatenate([ground.island[0], ground.island[-1], ground.island[:, 0], ground.island[:, -1]])
    assert not edges.any()                          # the island closes inside the crop


@needs_inputs
def test_terrain_reproduces_the_original_island_and_soundings(ground):
    original = ORIGINAL['terrain']
    # the same pixels: their count, and the sha256 of the packed mask the original saved
    assert list(ground.island.shape) == original['island_mask_shape']
    assert int(ground.island.sum()) == original['island_pixels']
    assert hashlib.sha256(np.packbits(ground.island).tobytes()).hexdigest() == original['island_mask_sha256']
    assert ground.cell_area_km2()[ground.island].sum() == pytest.approx(original['island_area_km2'], abs=0.01)
    # the original interpolated a 4 px/deg geoid and set its level 0.037 m higher: across the crop 19 pixels move by
    # up to 0.527 m, on the island the largest move is 0.368 m, the summit moves 0.089 m and a sounding 0.142 m
    assert ground.z[ground.island].max() == pytest.approx(original['summit_m'], abs=0.1)
    recorded = {s['name']: s for s in ORIGINAL['soundings']}
    for s in terrain.soundings(ground):
        r = recorded[s['name']]
        assert (s['longitude_deg'], s['latitude_deg']) == (r['longitude_deg'], r['latitude_deg'])
        assert s['depth_m'] == pytest.approx(r['depth_m'], abs=0.15)
        assert s['display_depth_m'] == r['display_depth_m']


# The renderer's arithmetic

def made_up_terrain():
    """A Terrain on the flyover's grids with made-up heights, for the arithmetic that needs no inputs."""
    grid = terrain.regional_grid(16)
    lat4, lon4 = 90.0 - (np.arange(720) + 0.5) / 4, (np.arange(1440) + 0.5) / 4 - 360.0
    shape = (grid.lat_deg.size, grid.lon_deg.size)
    return terrain.Terrain(lat=grid.lat_deg, lon=grid.lon_deg - 360.0, z=np.zeros(shape), island=np.zeros(shape, bool),
                           lat4=lat4, lon4=lon4, z4=np.add.outer(lat4, 1000 * lon4), sea_level_m=0.0, water_share=0.28,
                           pixels_per_degree=16, summit=(0.0, 0.0, 0.0), hashes={})


def test_depth_texels_sit_on_the_pixel_centres():
    # OpenGL samples texel i of n at its centre, (i + 0.5) / n: every pixel centre must fall there
    t = made_up_terrain()
    bounds = t.texture_bounds()
    u = terrain.texture_uv(t.lon, np.full(t.lon.shape, t.lat[0]), bounds)[:, 0]
    v = terrain.texture_uv(np.full(t.lat.shape, t.lon[0]), t.lat, bounds)[:, 1]
    np.testing.assert_allclose(u * t.lon.size - 0.5, np.arange(t.lon.size), atol=1e-9)
    np.testing.assert_allclose(v * t.lat.size - 0.5, np.arange(t.lat.size), atol=1e-9)


def test_the_surround_closes_its_seam():
    t = made_up_terrain()
    lon, lat, z = t.surround()
    assert lon[-1] - lon[0] == pytest.approx(360.0) and np.allclose(np.diff(lon), 0.5)
    assert np.array_equal(lat, t.lat4[::2]) and z.shape == (lat.size, lon.size)
    assert np.array_equal(z[:, -1], z[:, 0]) and np.array_equal(z[:, :-1], t.z4[::2, ::2])
