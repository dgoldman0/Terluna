"""Build the explorer page's display assets from checked products; output is ignored.

Each asset is a compact display encoding of one product. The page reads them as
pictures and lookup tables; the light values themselves come from the calendar
products through the shared evaluator.

    sky-<world>.bin   Solved spherical sky atlas (illumination/sky), log luminance and
                      chromaticity on the atlas's own Sun, view-elevation and azimuth grid
    earth.jpg         NASA Blue Marble surface and cloud composite (photographic prior)
    moon.jpg          Seas by depth and land by height with relief shading (geography atlas)
    moon-height.png   Height above sea level on a signed square-root scale (geography atlas)
    places.json       Calendar coasts, heritage sites, seas and islands with their settings
    stars.json        Yale Bright Star Catalogue stars to magnitude 5.5

Missing sources are reported with the command that restores them. A hash
mismatch stops the build; no other source is substituted.
"""
from __future__ import annotations

import hashlib
import io
import json
from pathlib import Path
import sys
import urllib.request

import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from shared.constants import MOON_RADIUS  # noqa: E402

SKY_PRODUCT = ROOT / 'illumination/sky/results/solved_sky.json'
EARTH_MANIFEST = ROOT / 'visualization/sea-appearance/nubium-earth-inputs.json'
EARTH_DIR = ROOT / 'visualization/sea-appearance/results/earth-inputs'
GEO_GRID = ROOT / 'geography/products/atlas_28pct_4ppd.npz'
GEO_ATLAS = ROOT / 'geography/results/atlas.json'
HERITAGE = ROOT / 'research/studies/conservation/results/heritage_register.json'
STARS = ROOT / 'illumination/stars/bright_stars.json'
CALENDAR = ROOT / 'illumination/calendar/results/astronomy.json'
WORLDS = {'moon': 'moon_1.2atm', 'earth': 'earth_control'}
LOG_RANGE = (-16.0, 6.0)          # log10 cd m^-2 covered by the uint16 luminance code


class MissingSource(FileNotFoundError):
    pass


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def sky_atlas(world: str) -> tuple[bytes, dict, dict]:
    """Encode one solved-sky atlas; the recorded atlas hash must match."""
    product = json.loads(SKY_PRODUCT.read_text())
    record = product['worlds'][WORLDS[world]]['atlas']
    path = ROOT / record['path']
    if not path.exists():
        raise MissingSource(f'{rel(path)} (regenerate with the optical comfort study: '
                            'research/studies/optical_comfort/README.md)')
    if digest(path) != record['sha256']:
        raise ValueError(f'Sky atlas {rel(path)} differs from solved_sky.json')
    data = np.load(path)
    xyz = np.asarray(data['xyz'], dtype=np.float64)              # sun, elevation, azimuth, XYZ
    total = np.maximum(xyz.sum(axis=-1), 1e-300)
    y = np.clip(xyz[..., 1], 10 ** LOG_RANGE[0], 10 ** LOG_RANGE[1])
    code_y = np.round((np.log10(y) - LOG_RANGE[0]) / (LOG_RANGE[1] - LOG_RANGE[0]) * 65535)
    code_x = np.round(np.clip(xyz[..., 0] / total, 0, 1) * 65535)
    code_c = np.round(np.clip(xyz[..., 1] / total, 0, 1) * 65535)
    packed = np.stack([code_y, code_x, code_c], axis=-1).astype('<u2')
    flux = np.asarray(data['direct_horizontal_xyz'])[:, 1] + np.asarray(data['diffuse_horizontal_xyz'])[:, 1]
    beam = np.asarray(data['direct_normal_xyz'], dtype=np.float64)
    beam_xy = (beam[:, :2] / np.maximum(beam.sum(axis=1), 1e-300)[:, None]).tolist()
    header = dict(world=WORLDS[world], suns=data['suns'].tolist(), elevations=data['elevations'].tolist(),
                  azimuths=data['azimuths'].tolist(), log_range=list(LOG_RANGE),
                  horizontal_lux=flux.tolist(),
                  diffuse_lux=np.asarray(data['diffuse_horizontal_xyz'])[:, 1].tolist(),
                  direct_normal_lux=beam[:, 1].tolist(), direct_xy=beam_xy,
                  layout='uint16 little-endian [sun][elevation][azimuth][log10 Y, x, y]',
                  reading_rule='Y is luminance in cd m^-2; azimuth is measured from the Sun, symmetric; '
                               'elevations above the horizon only.')
    return packed.tobytes(), header, {rel(path): record['sha256']}


def ensure_earth_inputs() -> tuple[list[Path], dict]:
    manifest = json.loads(EARTH_MANIFEST.read_text())
    paths, hashes = [], {}
    for item in manifest['files'][:2]:
        path = EARTH_DIR / item['filename']
        if not path.exists():
            EARTH_DIR.mkdir(parents=True, exist_ok=True)
            print('Fetching', item['url'])
            with urllib.request.urlopen(item['url'], timeout=120) as response:
                payload = response.read()
            if hashlib.sha256(payload).hexdigest() != item['sha256']:
                raise ValueError(f'Downloaded {item["filename"]} does not match its recorded hash')
            path.write_bytes(payload)
        if digest(path) != item['sha256']:
            raise ValueError(f'{rel(path)} differs from {rel(EARTH_MANIFEST)}')
        paths.append(path)
        hashes[rel(path)] = item['sha256']
    credits = [dict(credit=i['credit'], source=i['source_page'], reading_rule=i['reading_rule'])
               for i in manifest['files'][:2]]
    return paths, dict(hashes=hashes, credits=credits)


def earth_texture() -> tuple[bytes, dict]:
    (surface_path, cloud_path), record = ensure_earth_inputs()
    surface = np.asarray(Image.open(surface_path).convert('RGB'), dtype=np.float64) / 255
    cloud = np.asarray(Image.open(cloud_path).convert('L'), dtype=np.float64) / 255
    # Linear-light composite: clouds as a white layer. The map's faint haze floor is
    # removed so open ocean keeps its colour; thick cloud reaches full opacity.
    linear = surface ** 2.2
    alpha = (np.clip((cloud - 0.14) / 0.66, 0, 1) ** 1.15)[..., None]
    mixed = linear * (1 - alpha) + 0.92 * alpha
    image = Image.fromarray(np.round(np.clip(mixed, 0, 1) ** (1 / 2.2) * 255).astype(np.uint8))
    image = image.resize((1024, 512), Image.LANCZOS)
    out = io.BytesIO()
    image.save(out, 'JPEG', quality=86, optimize=True, progressive=True)
    return out.getvalue(), record


def load_geography():
    if not GEO_GRID.exists():
        raise MissingSource(f'{rel(GEO_GRID)} (regenerate with: python -m geography.atlas)')
    grid = np.load(GEO_GRID)
    meta = json.loads(str(grid['metadata']))
    atlas = json.loads(GEO_ATLAS.read_text())
    if meta['producer'] != atlas['producer'] or abs(meta['sea_level_m'] - atlas['sea_level_m']) > 0.05:
        raise ValueError('The geography grid product and results/atlas.json come from different runs')
    return grid, meta, atlas


def colour_ramp(value, stops):
    xs = np.array([s[0] for s in stops])
    cs = np.array([s[1] for s in stops], dtype=np.float64)
    return np.stack([np.interp(value, xs, cs[:, k]) for k in range(3)], axis=-1)


def moon_textures(grid, meta) -> tuple[bytes, bytes, dict]:
    height = np.asarray(grid['height_m'], dtype=np.float64)
    lat = np.asarray(grid['lat_deg'])
    lon = np.asarray(grid['lon_deg'])
    level = float(meta['sea_level_m'])
    water = np.asarray(grid['water_label']) > 0
    if lat[0] < lat[-1]:                                            # images run north to south
        height, water, lat = height[::-1], water[::-1], lat[::-1]
    shift = int(np.argmin(np.abs(((lon + 180) % 360) - 180 - (-179.875))))
    height, water = np.roll(height, -shift, axis=1), np.roll(water, -shift, axis=1)
    rel_h = height - level
    ppd = 4
    dy = MOON_RADIUS * np.pi / 180 / ppd
    dx = dy * np.cos(np.radians(lat))[:, None]
    gy, gx = np.gradient(np.where(water, level, height))
    nx, ny = -gx / dx * 2.2, gy / dy * 2.2                           # relief exaggerated for display
    norm = np.sqrt(nx ** 2 + ny ** 2 + 1)
    light = np.array([-0.5, 0.5, 0.707])                            # from the north-west, 45 degrees up
    shade = (nx * light[0] + ny * light[1] + light[2]) / norm
    shade = np.clip(shade / light[2], 0.35, 1.45)
    depth = np.where(water, -rel_h, 0)
    sea = colour_ramp(np.sqrt(np.maximum(depth, 0)), [
        (0, (0.40, 0.84, 0.86)), (7, (0.20, 0.66, 0.80)), (16, (0.11, 0.45, 0.71)),
        (32, (0.07, 0.29, 0.58)), (55, (0.05, 0.17, 0.42)), (100, (0.03, 0.10, 0.30))])
    land = colour_ramp(np.maximum(rel_h, 0), [
        (0, (0.78, 0.71, 0.55)), (800, (0.71, 0.62, 0.47)), (2500, (0.60, 0.53, 0.43)),
        (5000, (0.62, 0.59, 0.55)), (8000, (0.86, 0.84, 0.81)), (12000, (0.95, 0.94, 0.92))])
    land = land * shade[..., None]
    sea = sea * (0.82 + 0.18 * np.clip(shade, 0.6, 1.3))[..., None]
    rgb = np.where(water[..., None], sea, land)
    image = Image.fromarray(np.round(np.clip(rgb, 0, 1) * 255).astype(np.uint8))
    out = io.BytesIO()
    image.save(out, 'JPEG', quality=88, optimize=True, progressive=True)
    # Height at 2 px/deg, 128 = sea level, signed square root of height / 10 km.
    coarse = rel_h.reshape(360, 2, 720, 2).mean(axis=(1, 3))
    code = np.clip(np.round(128 + np.sign(coarse) * np.sqrt(np.abs(coarse) / 10000) * 127), 0, 255)
    wet = water.reshape(360, 2, 720, 2).mean(axis=(1, 3)) >= 0.5
    code = np.where(wet, np.minimum(code, 127), np.maximum(code, 129))
    height_png = io.BytesIO()
    Image.fromarray(code.astype(np.uint8)).save(height_png, 'PNG', optimize=True)
    info = dict(sea_level_m=level, share=meta.get('share'), texture='1440x720, longitude -180..180 east, '
                'latitude 90..-90', height_code='value = 128 + sign(h) * sqrt(|h| / 10 km) * 127; '
                'below 128 is sea, above 128 is land')
    return out.getvalue(), height_png.getvalue(), info


def places(atlas) -> dict:
    astronomy = json.loads(CALENDAR.read_text())
    heritage = json.loads(HERITAGE.read_text())
    if abs(heritage['sea_level_m'] - atlas['sea_level_m']) > 0.5:
        raise ValueError('The heritage register uses another sea level than the geography atlas')
    coasts = [dict(id=f'coast-{s["id"]}', name=s['name'], description=s['description'],
                   latitude=s['latitude'], longitude=s['longitude'], group='Study coasts')
              for s in astronomy['sites']]
    sites = []
    for s in heritage['sites']:
        sites.append(dict(id=s['id'], name=s['name'], latitude=s['lat'], longitude=s['lon'],
                          group='Landing sites and heritage', year=s.get('year'), kind=s.get('kind'),
                          status=s.get('status'), depth_m=s.get('depth_m'),
                          water_body=s.get('water_body'), marks=s.get('marks')))
    seas = []
    for b in atlas['bodies']:
        if not b.get('water_names'):
            continue
        lat_c, lon_c = b['centroid'][1], b['centroid'][0]
        seas.append(dict(id=f'sea-{b["label"]}', name=b['water_names'][0], latitude=lat_c,
                         longitude=((lon_c + 180) % 360) - 180, group='Seas',
                         area_km2=b['area_km2'], mean_depth_m=b['mean_depth_m'],
                         max_depth_m=b['max_depth_m'], also=b['water_names'][1:6]))
    islands = []
    for i in sorted(atlas['islands'], key=lambda x: -x['area_km2'])[:14]:
        if not i.get('features'):
            continue
        islands.append(dict(id=f'island-{i["label"]}', name=i['features'][0], latitude=i['summit'][1],
                            longitude=((i['summit'][0] + 180) % 360) - 180, group='Island summits',
                            area_km2=i['area_km2'], summit_m=i['summit_m'], also=i['features'][1:4]))
    return dict(schema='terluna.visualization.explorer-places/1',
                reading_rule='Coasts are the calendar sites; heritage settings, depths and water bodies '
                             'come from the conservation register at the atlas sea level; sea centroids '
                             'and island summits from the geography atlas. Coordinates in degrees north '
                             'and east.',
                groups=[coasts, sites, seas, islands])


def stars() -> dict:
    catalogue = json.loads(STARS.read_text())
    keep = [s for s in catalogue['stars'] if s[2] <= 5.5]
    return dict(schema='terluna.visualization.explorer-stars/1', columns=catalogue['columns'],
                source=catalogue['source'], notes=catalogue['notes'],
                stars=[[round(v, 5) for v in s] for s in keep])


def build(destination: Path) -> dict:
    """Write every asset into destination/assets and return their source record."""
    out = destination / 'assets'
    out.mkdir(parents=True, exist_ok=True)
    record = dict(sources={}, assets={})

    def write(name, payload, sources):
        (out / name).write_bytes(payload)
        record['assets'][name] = hashlib.sha256(payload).hexdigest()
        record['sources'].update(sources)

    headers = {}
    for world in WORLDS:
        payload, header, sources = sky_atlas(world)
        headers[world] = header
        write(f'sky-{world}.bin', payload, sources)
    record['sources'][rel(SKY_PRODUCT)] = digest(SKY_PRODUCT)
    earth, earth_record = earth_texture()
    write('earth.jpg', earth, earth_record['hashes'])
    grid, meta, atlas = load_geography()
    texture, height_png, moon_info = moon_textures(grid, meta)
    write('moon.jpg', texture, {rel(GEO_GRID): digest(GEO_GRID), rel(GEO_ATLAS): digest(GEO_ATLAS)})
    write('moon-height.png', height_png, {})
    write('places.json', (json.dumps(places(atlas), separators=(',', ':')) + '\n').encode(),
          {rel(HERITAGE): digest(HERITAGE), rel(CALENDAR): digest(CALENDAR)})
    write('stars.json', (json.dumps(stars(), separators=(',', ':')) + '\n').encode(), {rel(STARS): digest(STARS)})
    manifest = dict(schema='terluna.visualization.explorer-assets/1', sky=headers, moon=moon_info,
                    earth_credits=earth_record['credits'],
                    evidence='Display encodings of the named products; light values come from the '
                             'calendar products through evaluator.mjs.')
    write('assets.json', (json.dumps(manifest, separators=(',', ':')) + '\n').encode(), {})
    return record


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--out', type=Path, default=HERE / 'dist')
    result = build(parser.parse_args().out)
    for name, value in result['assets'].items():
        print(f'{name:18s} {value[:12]}')
