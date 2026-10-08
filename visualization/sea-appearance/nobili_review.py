"""Image-space fidelity screens for Nobili; not physical or perceptual validation.

Fit the illuminated lower limb jointly with its centre: a fixed-centre radial
screen confounds small placement shifts with size changes. The partial-limb fit
is itself approximate at the generated image's resolution. Brightness metrics
are linearized display values, not recovered scene radiance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.optimize import least_squares
from scipy.stats import spearmanr


def luminance(path):
    rgb = np.asarray(Image.open(path).convert('RGB'), dtype=float) / 255
    linear = np.where(rgb <= .04045, rgb / 12.92, ((rgb + .055) / 1.055) ** 2.4)
    return linear @ np.array([.2126729, .7151522, .072175])


def measure(path, display):
    y = luminance(path)
    height, width = y.shape
    rw, rh = display['camera']['frame']
    projection = display['earth_projection']
    cx, cy = np.array(projection['centre']) * [width / rw, height / rh]
    rx, ry = np.array(projection['disk_dimensions_px']) / 2 * [width / rw, height / rh]
    phi = np.linspace(.15, np.pi - .15, 100)
    radius = np.linspace(.7, 1.35, 181)
    x = cx + rx * np.cos(phi)[:, None] * radius
    yy = cy + ry * np.sin(phi)[:, None] * radius
    sample = map_coordinates(gaussian_filter(y, .5), [yy, x], order=1)
    gradient = -np.gradient(sample, radius, axis=1)
    edge = radius[np.argmax(gradient, axis=1)]
    ex = cx + rx * np.cos(phi) * edge
    ey = cy + ry * np.sin(phi) * edge
    fit = least_squares(
        lambda q: np.hypot((ex - cx - q[0]) / rx, (ey - cy - q[1]) / ry) - q[2],
        [0, 0, 1], loss='soft_l1', f_scale=.015,
        bounds=([-rx * .25, -ry * .25, .7], [rx * .25, ry * .25, 1.3]))
    lo, hi = int(.265 * height), int(.37 * height)
    skyline = (lo + np.argmax(-np.diff(gaussian_filter(y, (1, 0))[lo:hi], axis=0), axis=0)) / height
    # Compare the lower illuminated interior after registering the fitted limb.
    u, v = np.meshgrid(np.linspace(-.8, .8, 81), np.linspace(.03, .8, 41))
    interior = u * u + v * v < .64
    texture = map_coordinates(y, [cy + fit.x[1] + ry * fit.x[2] * v[interior],
                                  cx + fit.x[0] + rx * fit.x[2] * u[interior]], order=1)
    record = dict(frame=[width, height], sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        fitted_lower_limb=dict(centre_shift_px=fit.x[:2].tolist(), radius_scale=float(fit.x[2]),
                              edge_rms_radial_fraction=float(np.sqrt(np.mean(fit.fun ** 2)))),
        fixed_centre_radial_screen_p10_median_p90=np.percentile(edge, [10, 50, 90]).tolist(),
        rim_minmax=[float(skyline.min()), float(skyline.max())],
        median_linear_display_Y=dict(sky=float(np.median(y[:int(.12 * height)])),
            rim=float(np.median(y[int(.39 * height):int(.52 * height)])),
            water=float(np.median(y[int(.64 * height):]))))
    return record, skyline, texture


def main():
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference', type=Path, default=here / 'illustrations/nobili-night-reference.png')
    p.add_argument('--display', type=Path, default=here / 'illustrations/nobili-night-display.json')
    p.add_argument('--image', type=Path, action='append')
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    display = json.loads(a.display.read_text())
    ref, skyline, texture = measure(a.reference, display)
    records = {}
    for path in a.image or [here / 'illustrations/nobili-night.png']:
        record, line, values = measure(path, display)
        expected = np.interp(np.linspace(0, 1, len(line)), np.linspace(0, 1, len(skyline)), skyline)
        record['skyline_rms_fraction_of_height'] = float(np.sqrt(np.mean((line - expected) ** 2)))
        record['registered_earth_interior_rank_correlation'] = float(spearmanr(texture, values).statistic)
        records[path.name] = record
    report = dict(schema='terluna.visualization.nobili-image-checks/1', reference=ref, images=records,
        method='Lower illuminated limb: gradient samples with joint robust centre/radius fit. Skyline: strongest vertical edge in its known band. Earth texture: rank correlation in the registered lower illuminated interior.',
        limits=['Image-space comparison only; no recovered radiometry, wave-height measurement, human-vision validation or proof of accurate geography.',
                'A fixed-centre radial screen confounds disk size with placement; the fitted radius is the size estimate.',
                'The partial-limb fit is approximate and exposure-dependent; it is not an astrometric measurement.'],
        producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
