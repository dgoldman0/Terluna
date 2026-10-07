"""Build the public sky page's display assets from checked products; output is ignored.

Each asset is a compact display encoding of one product. The page reads them as
pictures and lookup tables; the light values themselves come from the calendar
products through the shared evaluator.

    sky-moon.bin   Solved spherical sky atlas of the Open Moon (illumination/sky): log
                   luminance and chromaticity on the atlas's Sun, view-elevation and azimuth grid
    earth.jpg      NASA Blue Marble surface and cloud composite (photographic prior)
    places.json    Six places with their names and settings (calendar coasts, the
                   conservation register at the geography atlas sea level)
    stars.json     Yale Bright Star Catalogue stars to magnitude 5.5
    assets.json    The atlas grid and the light of the Earth control sky for comparisons

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


# The page's places: names and one-line settings written for visitors; the
# coordinates, depths and heights come from the products named in `source`.
PLACES = [
    dict(id='storms', source=('coast', '0'), name='Ocean of Storms', latin='Oceanus Procellarum',
         line='A western shore where Earth hangs low over the sea.'),
    dict(id='showers', source=('coast', '1'), name='Sea of Showers', latin='Mare Imbrium',
         line='A southern shore below the Carpathian Mountains, with Earth high overhead.'),
    dict(id='tranquility', source=('heritage', 'apollo-11'), name='Tranquility Base', latin='Apollo 11, 1969',
         line='Where people first walked on the Moon, now under {depth} m of sea.'),
    dict(id='descartes', source=('heritage', 'apollo-16'), name='Descartes Highlands', latin='Apollo 16, 1972',
         line='Dry highlands {height} m above the sea, where Apollo 16 landed.'),
    dict(id='smyth', source=('coast', '4'), name="Smyth's Sea", latin='Mare Smythii',
         line="On the Moon's eastern edge, where Earth rises and sets on the sea horizon."),
    dict(id='cleverness', source=('coast', '5'), name='Sea of Cleverness', latin='Mare Ingenii',
         line='On the far side, where Earth never rises.'),
]


def places(atlas) -> dict:
    astronomy = json.loads(CALENDAR.read_text())
    heritage = json.loads(HERITAGE.read_text())
    if abs(heritage['sea_level_m'] - atlas['sea_level_m']) > 0.5:
        raise ValueError('The heritage register uses another sea level than the geography atlas')
    coasts = {s['id']: s for s in astronomy['sites']}
    sites = {s['id']: s for s in heritage['sites']}
    out = []
    for p in PLACES:
        kind, key = p['source']
        if kind == 'coast':
            site = coasts[key]
            lon, lat, setting = site['longitude'], site['latitude'], 'coast'
            line = p['line']
        else:
            site = sites[key]
            lon, lat = site['lon'], site['lat']
            setting = 'sea' if site['status'] == 'submerged' else 'land'
            height = round(site['height_above_geoid_m'] - atlas['sea_level_m'], -1)
            line = p['line'].format(depth=f"{site.get('depth_m') or 0:,}", height=f'{height:,.0f}')
        out.append(dict(id=p['id'], name=p['name'], latin=p['latin'], line=line, latitude=lat,
                        longitude=lon, setting=setting))
    return dict(schema='terluna.visualization.sky-places/2',
                reading_rule='Coordinates in degrees north and east. `setting` is coast, sea or land; '
                             'depths and heights come from the conservation register at the atlas sea level.',
                places=out)


def stars() -> dict:
    catalogue = json.loads(STARS.read_text())
    keep = [s for s in catalogue['stars'] if s[2] <= 5.5]
    return dict(schema='terluna.visualization.explorer-stars/1', columns=catalogue['columns'],
                source=catalogue['source'], notes=catalogue['notes'],
                stars=[[round(v, 5) for v in s] for s in keep])


def build(destination: Path) -> dict:
    """Write every asset into destination/assets and return their source record."""
    out = destination / 'assets'
    if out.exists():
        for old in out.iterdir():
            old.unlink()
    out.mkdir(parents=True, exist_ok=True)
    record = dict(sources={}, assets={})

    def write(name, payload, sources):
        (out / name).write_bytes(payload)
        record['assets'][name] = hashlib.sha256(payload).hexdigest()
        record['sources'].update(sources)

    payload, moon, sources = sky_atlas('moon')
    write('sky-moon.bin', payload, sources)
    _, earth_sky, sources = sky_atlas('earth')
    record['sources'].update(sources)
    record['sources'][rel(SKY_PRODUCT)] = digest(SKY_PRODUCT)
    earth, earth_record = earth_texture()
    write('earth.jpg', earth, earth_record['hashes'])
    atlas = json.loads(GEO_ATLAS.read_text())
    write('places.json', (json.dumps(places(atlas), separators=(',', ':')) + '\n').encode(),
          {rel(GEO_ATLAS): digest(GEO_ATLAS), rel(HERITAGE): digest(HERITAGE), rel(CALENDAR): digest(CALENDAR)})
    write('stars.json', (json.dumps(stars(), separators=(',', ':')) + '\n').encode(), {rel(STARS): digest(STARS)})
    manifest = dict(schema='terluna.visualization.sky-assets/2', sky=moon,
                    earth_control=dict(suns=earth_sky['suns'], horizontal_lux=earth_sky['horizontal_lux'],
                                       reading_rule='Horizontal light on the ground under the solver\'s Earth '
                                                    'control sky, by Sun elevation; used for comparisons.'),
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
