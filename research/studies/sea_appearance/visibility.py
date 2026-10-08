"""Visibility screening of near-full Earth and a terrestrial full-Moon control.

Consumes the sea study's physical light, an observed EPIC image and quantitative
LROC reflectance. Runs unmodified HDR-VDP-3.0.7 externally. Predictions concern
achromatic detail detection under specified conditions, not exact appearance,
colour recognition, an individual's eyesight, or a validated outdoor scene.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import subprocess
import zipfile
from pathlib import Path

import numpy as np
from scipy.io import savemat
from scipy.ndimage import map_coordinates, gaussian_filter1d


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        while data := f.read(1024 * 1024):
            h.update(data)
    return h.hexdigest()


def verify_model(archive, directory):
    """Reject modified external model files instead of trusting a version label."""
    checked = {}
    with zipfile.ZipFile(archive) as z:
        for info in z.infolist():
            if info.is_dir():
                continue
            parts = Path(info.filename).parts
            relative = Path(*parts[1:])
            path = directory / relative
            expected = hashlib.sha256(z.read(info)).hexdigest()
            if digest(path) != expected:
                raise ValueError('External model file mismatch: '+str(relative))
            checked[str(relative)] = expected
    return dict(file_count=len(checked),
                tree_sha256=hashlib.sha256(json.dumps(checked,sort_keys=True).encode()).hexdigest())


def grid(ppd, field_deg):
    """Small rectilinear angular window, with exact ray solid-angle weights."""
    n = int(round(ppd * field_deg))
    focal = ppd * 180 / np.pi
    q = (np.arange(n) + .5 - n / 2) / focal
    x, y = np.meshgrid(q, -q)
    ct = 1 / np.sqrt(1 + x*x + y*y)
    return x, y, ct**3 / focal**2, ct


def normalize_beam(direct, omega, cosine, illuminance):
    integral = float(np.sum(direct * omega * cosine))
    if integral <= 0 or not np.isfinite(integral):
        raise ValueError('Invalid direct-light integral')
    return direct * (illuminance / integral)


def featureless_inner(reference, radius_fraction, omega, cosine):
    """Remove azimuthal interior texture, preserving limb and direct integral.

    An annular mean supplies a smooth radial profile. A cosine taper restricts
    changes to r<0.9, with full removal at r<0.7. A weighted offset preserves the
    beam integral. Thus detection cannot be explained by changing disk size,
    limb, or total incident light. The result is an explicit counterfactual.
    """
    r = radius_fraction
    nbin = max(16, round(np.sqrt(np.count_nonzero(r < 1)) / 2))
    k = np.clip((r * nbin).astype(int), 0, nbin-1)
    inside = r < 1
    count = np.bincount(k[inside], minlength=nbin)
    total = np.bincount(k[inside], weights=reference[inside], minlength=nbin)
    good = count > 0
    profile = np.interp(np.arange(nbin), np.flatnonzero(good), total[good]/count[good])
    profile = gaussian_filter1d(profile, 1)
    smooth = np.interp(r * nbin - .5, np.arange(nbin), profile)
    window = .5 * (1 + np.cos(np.pi * np.clip((r-.7)/.2, 0, 1)))
    residual = reference-smooth
    weights = window * omega * cosine
    correction = float(np.sum(residual*weights)/weights.sum())
    counter = reference-window*(residual-correction)
    if np.min(counter) <= 0:
        raise ValueError('Counterfactual is nonpositive; do not silently clip its energy')
    return counter, r < .7


def epic_pattern(path, band, x, y, radius):
    """Observed within-band structure; no interpretation of early absolute units."""
    import h5py
    with h5py.File(path) as f:
        a = f[f'Band{band}nm'][:]
        lat = f['Geolocation/Earth/Latitude'][:]
        mask = np.isfinite(lat) & (np.abs(lat) <= 90) & np.isfinite(a)
        yy, xx = np.where(mask)
        cx, cy = (xx.min()+xx.max())/2, (yy.min()+yy.max())/2
        rx, ry = (xx.max()-xx.min()+1)/2, (yy.max()-yy.min()+1)/2
        a = np.where(mask, np.maximum(a, 0), 0)
        qx = cx + x/math.tan(radius)*rx
        qy = cy - y/math.tan(radius)*ry
        pattern = map_coordinates(a, [qy, qx], order=1, mode='constant', cval=0)
        record = dict(filename=Path(path).name, band_nm=band,
                      time=str(f[f'Band{band}nm'].attrs['time'].decode()),
                      centre_lon_lat_deg=[float(f.attrs['centroid_mean_longitude'][0]),
                                          float(f.attrs['centroid_mean_latitude'][0])],
                      source_disk_bbox=[int(xx.min()), int(yy.min()), int(xx.max()), int(yy.max())],
                      interpretation='Early public EPIC test product; use observed relative within-band structure only. Absolute scale is supplied by the sea study. Historical viewing geometry/clouds, not the 2038 surface map.')
    return pattern, record


def pds_tile(path):
    with Path(path).open('rb') as f:
        header = f.read(10000).decode('ascii', errors='replace').split('\r\nEND\r\n')[0]
    def value(key):
        m = re.search(r'^\s*'+re.escape(key)+r'\s*=\s*([^\r\n<]+)', header, re.M)
        if not m:
            raise ValueError('Missing PDS field: '+key)
        return m[1].strip()
    if value('SAMPLE_TYPE') != 'PC_REAL' or int(value('SAMPLE_BITS')) != 32:
        raise ValueError('Expected little-endian 32-bit PDS reflectance')
    rows, cols = int(value('LINES')), int(value('LINE_SAMPLES'))
    offset = (int(value('^IMAGE'))-1)*int(value('RECORD_BYTES'))
    a = np.memmap(path, dtype='<f4', mode='r', offset=offset, shape=(rows, cols))
    return a, dict(ppd=float(value('MAP_RESOLUTION')),
                   line_offset=float(value('LINE_PROJECTION_OFFSET')),
                   sample_offset=float(value('SAMPLE_PROJECTION_OFFSET')),
                   west=float(value('WESTERNMOST_LONGITUDE')),
                   product_version=value('PRODUCT_VERSION_ID'))


def moon_pattern(inputs, x, y, radius):
    """Orthographic near-side reflectance control at 0 longitude/latitude.

    The normalized 566-nm map is not a full-phase BRDF. Its spatial pattern is
    anchored to measured full-Moon mean luminance. Unmapped polar caps are
    boundary-extended and excluded from the interior test region.
    """
    u, v = x/math.tan(radius), y/math.tan(radius)
    inside = u*u+v*v <= 1
    z = np.sqrt(np.maximum(1-u*u-v*v, 0))
    lon = np.degrees(np.arctan2(u, z)) % 360
    lat = np.degrees(np.arcsin(np.clip(v, -1, 1)))
    output = np.zeros_like(x)
    extended = inside & (np.abs(lat) > 69.99)
    lat = np.clip(lat, -69.99, 69.99)
    records = []
    for side in 'NS':
        for centre in ['0450', '3150']:
            filename = f'WAC_HAPKE_566NM_E350{side}{centre}.IMG'
            a, meta = pds_tile(inputs/filename)
            select = inside & ((lat >= 0) if side == 'N' else (lat < 0))
            select &= (lon < 90) if centre == '0450' else (lon >= 270)
            col = meta['sample_offset'] + lon[select]*meta['ppd']
            row = meta['line_offset'] - lat[select]*meta['ppd']
            values = map_coordinates(a, [row, col], order=1, mode='nearest')
            if np.any(~np.isfinite(values) | (values <= 0)):
                raise ValueError('Invalid or missing LROC samples in control; no synthetic substitution')
            output[select] = values
            records.append(dict(filename=filename, **meta))
    return output, dict(tiles=records, polar_extension_pixel_fraction=float(extended.sum()/inside.sum()),
                        interpretation='LROC 566-nm relative spatial pattern at its reference g=i=60 degrees, e=0; not a full-phase photometric prediction. Disk mean independently anchored to the measured bright full-Moon range.')


def physical_light(root):
    from research.studies.sea_appearance import lighting, regimes
    from shared.constants import EARTH_RADIUS, MOON_RADIUS, EARTH_MOON_DISTANCE
    product = json.loads((root/'research/studies/sea_appearance/results/nubium-midnight.json').read_text())
    s = product['scene']
    g = lighting.geometry(np.array([product['jd_tt']]), *s['camera']['lon_lat_deg'])
    sky = lighting.Sky()
    weights = lighting.Earthlight(sky).weights(g['earth_phase_angle'], g['earth_distance_m'], g['earth_sunlight_factor'])[0]
    el = float(g['earth_elevation'][0])
    beam = regimes.normal_beam(sky, el)
    fields = regimes.Fields()
    background = regimes.lookup(fields.field(el, weights), el, 0)
    sun_el, sun_az = float(g['sun_elevation'][0]), float(g['sun_azimuth'][0])
    background += regimes.lookup(fields.field(sun_el), el, abs((float(g['earth_azimuth'][0])-sun_az+180)%360-180))
    floor = lighting.FLOOR_LUX/np.pi
    background += floor*np.array([regimes.D65[0]/regimes.D65[1], 1, (1-sum(regimes.D65))/regimes.D65[1]])
    above, below = weights@fields.xyz, (weights*beam)@fields.xyz
    radius = math.asin(EARTH_RADIUS/float(g['earth_distance_m'][0]))
    return dict(earth_diameter_deg=math.degrees(2*radius), earth_normal_xyz=below.tolist(),
                earth_above_air_xyz=above.tolist(), earth_direct_transmission=float(below[1]/above[1]),
                earth_mean_full_disk_cd_m2=float(below[1]/(np.pi*np.sin(radius)**2)),
                sky_xyz=background.tolist(), moon_diameter_deg=math.degrees(2*math.asin(MOON_RADIUS/EARTH_MOON_DISTANCE)),
                atmosphere='Study clear molecular 1.2-atm case; aerosols, lunar clouds and refraction omitted',
                units='Y of XYZ is photopic cd/m2 for radiance or lux for irradiance')


def prepare(root, inputs, out):
    out.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((Path(__file__).with_name('visibility_inputs.json')).read_text())
    for row in manifest['files']:
        if digest(inputs/row['filename']) != row['sha256']:
            raise ValueError('Source hash mismatch: '+row['filename'])
    light = physical_light(root)
    model_paths = ['research/studies/sea_appearance/results/nubium-midnight.json',
                   'research/studies/sea_appearance/lighting.py',
                   'research/studies/sea_appearance/regimes.py',
                   'illumination/earthlight/model.py', 'illumination/earthlight/fetch_inputs.py',
                   'geography/lunar_ephemeris.py', 'shared/constants.json',
                   'illumination/sky/results/solved_sky.json',
                   'research/runs/optical_comfort/spherical/moon_1.2atm_standard.npz']
    model_hashes = {str(path):digest(root/path) for path in model_paths}
    cases, source_records = [], {}
    setups = [
        ('earth_551_age24', 'earth', 551, 24, 120, 6, 1),
        ('earth_551_age50', 'earth', 551, 50, 120, 6, 1),
        ('earth_551_age70', 'earth', 551, 70, 120, 6, 1),
        ('earth_quarter_contrast_age70', 'earth', 551, 70, 120, 6, .25),
        ('earth_443_age24', 'earth', 443, 24, 120, 6, 1),
        ('earth_680_age24', 'earth', 680, 24, 120, 6, 1),
        ('earth_60ppd', 'earth', 551, 24, 60, 6, 1),
        ('earth_wider_surround', 'earth', 551, 24, 60, 12, 1),
        ('moon_age24', 'moon', 566, 24, 120, 6, 1),
        ('moon_age70', 'moon', 566, 70, 120, 6, 1),
        ('moon_quarter_contrast_age70', 'moon', 566, 70, 120, 6, .25),
    ]
    cache = {}
    for name, body, band, age, ppd, field, strength in setups:
        key = body, band, ppd, field
        if key not in cache:
            x, y, omega, ct = grid(ppd, field)
            radius = math.radians(light[body+'_diameter_deg']/2)
            r = np.hypot(x, y)/math.tan(radius)
            if body == 'earth':
                pattern, src = epic_pattern(inputs/'epic_1b_20151117002712_00.h5', band, x, y, radius)
                normal = light['earth_normal_xyz'][1]
                background = light['sky_xyz'][1]
            else:
                pattern, src = moon_pattern(inputs, x, y, radius)
                normal = 5000*np.pi*np.sin(radius)**2
                background = .01  # declared control surround, not a measured coastal sky
            direct = normalize_beam(np.where(r<1, pattern, 0), omega, ct, normal)
            ref = background+direct
            counter, roi = featureless_inner(ref, r, omega, ct)
            source_records[body+str(band)] = src
            cache[key] = ref, counter, roi, r, omega, ct, background, normal
        ref, counter, roi, r, omega, ct, background, normal = cache[key]
        reference = counter + strength*(ref-counter)
        path = out/(name+'.mat')
        savemat(path, dict(reference=reference, comparison=counter, roi=roi.astype(np.uint8),
                          radius_fraction=r, ppd=float(ppd), background=float(background)), do_compression=True)
        integral = float(np.sum((reference-background)*omega*ct))
        counter_integral = float(np.sum((counter-background)*omega*ct))
        cases.append(dict(name=name, file=str(path), age=age, ppd=ppd, field_deg=field, identical=False,
                          contrast_fraction=strength, source=body+str(band), normal_lux=normal,
                          direct_integral_relative_error=abs(integral/normal-1),
                          counter_integral_relative_error=abs(counter_integral/normal-1),
                          interior_luminance_p10_p50_p90=np.percentile(reference[roi], [10,50,90]).tolist()))
    # An identical pair must return zero detection, independently of scene brightness.
    first = cases[0]
    cases.append({**first, 'name':'identical_control', 'identical':True})
    record = dict(schema='terluna.research.earth-visibility-inputs/1', producer_sha256=digest(__file__),
                  model_sha256=hashlib.sha256(json.dumps(model_hashes,sort_keys=True).encode()).hexdigest(),
                  model_inputs=model_hashes,
                  reading_rule='Absolute photopic luminance stimuli; compare interior structure at unchanged limb and normal illuminance. HDR-VDP detects differences, not scene recognition or exact appearance.',
                  physical_light=light, cases=cases, sources=source_records,
                  inputs={row['filename']:row['sha256'] for row in manifest['files']},
                  evidence='Computed stimuli and external vision-model screening; no new observer experiment.',
                  limitations=['Within-band Earth contrast used as an achromatic luminance proxy; three observed bands tested separately.',
                               'Early EPIC public test granule has inconsistent L1A/L1B metadata and unspecified radiance unit scale. No absolute EPIC calibration claimed.',
                               'Earth historical phase/view and clouds differ from the 2038 scene; the test concerns plausible observed spatial structure.',
                               'Moon map phase dependence and opposition contrast are not solved; mean luminance anchor is empirical, spatial pattern is a control.',
                               'Moon surround 0.01 cd/m2 is prescribed; there is no calibrated full coastal photograph or human field validation.',
                               'Static attended achromatic detection; colour recognition, gaze dynamics, eye disease and individual perception are not certified.'])
    (out/'stimuli.json').write_text(json.dumps(record, indent=2)+'\n')
    return record


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--inputs', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--octave', type=Path)
    p.add_argument('--hdrvdp', type=Path)
    p.add_argument('--image-package-list', default='')
    p.add_argument('--statistics-package-list', default='')
    args = p.parse_args()
    record = prepare(args.root, args.inputs, args.out)
    if args.octave:
        if not args.hdrvdp: p.error('--octave requires --hdrvdp')
        verified=verify_model(args.inputs/'hdrvdp-3.0.7.zip',args.hdrvdp)
        (args.out/'external-model.json').write_text(json.dumps(verified,indent=2)+'\n')
        config = dict(hdrvdp=str(args.hdrvdp), out=str(args.out), cases=record['cases'],
                      image_package_list=args.image_package_list, statistics_package_list=args.statistics_package_list)
        cp = args.out/'backend.json'; cp.write_text(json.dumps(config))
        subprocess.run([str(args.octave),'--no-history','--no-init-file','--quiet',
                        str(Path(__file__).with_name('visibility_backend.m')),str(cp)], check=True)


if __name__ == '__main__':
    main()
