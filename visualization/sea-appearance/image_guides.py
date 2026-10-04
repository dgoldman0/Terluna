"""Image-generation guides from the sea study's stored frame and coastal terrain.

This is a geometry and display-colour guide, not the spectral path tracer.
Sky and sea interpolate the scene's sRGB readouts; the sea has no waves here.
Terrain is projected from the geography product. Its illustrative Lambertian
shading conveys the stored Sun direction, without cast shadows or interreflection.
The AI-generated final needs a separate visual review against this guide.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from PIL import Image
from scipy.ndimage import map_coordinates


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rgb(value):
    return np.array([int(value[i:i + 2], 16) / 255 for i in (1, 3, 5)])


def linear(value):
    return np.where(value <= 0.04045, value / 12.92, ((value + 0.055) / 1.055) ** 2.4)


def encoded(value):
    value = np.maximum(value, 0)
    return np.where(value <= 0.0031308, 12.92 * value, 1.055 * value ** (1 / 2.4) - 0.055)


def field(points, width, height):
    """Interpolate the stored display colours down each sampled column, then across."""
    columns = []
    xs = []
    for side in ("left edge", "centre", "right edge"):
        selected = sorted((p for p in points if p["label"].endswith(side)), key=lambda p: p["y"])
        xs.append(selected[0]["x"])
        ys = [p["y"] for p in selected]
        colours = linear(np.array([rgb(p["srgb"]) for p in selected]))
        columns.append(np.stack([np.interp(np.arange(height), ys, colours[:, c]) for c in range(3)], -1))
    result = np.empty((height, width, 3))
    for y in range(height):
        for c in range(3):
            result[y, :, c] = np.interp(np.arange(width), xs, [a[y, c] for a in columns])
    return result


def water_structure(root, scene, camera, pixels, sea, region, out):
    """Projected domain wave realization; Fresnel contrast is a display approximation.

    Geometry comes from the study's explicitly borrowed Nectaris spectrum.
    The background sea colours remain the study's values. This quick guide does
    not solve radiance, facet shadowing or reflected terrain.
    """
    from illumination.water_surface.hotfile import read_hotfile
    from illumination.water_surface.model import wavenumber
    from illumination.water_surface import realization, reflection
    from research.studies.sea_appearance.scenes import unit
    from shared.constants import MOON_RADIUS, MOON_SURFACE_GRAVITY

    width, height = scene["camera"]["frame"]
    waves = scene["waves"]
    source = root / waves["file"]
    assert waves["source"].startswith("borrowed: the Nectaris restart spectrum at hour 816,")
    month_path = root / "climate/waves/results/shore_month.json"
    month = json.loads(month_path.read_text())
    node = next(s["lon_lat"] for s in month["shores"] if s["short"] == "Nectaris")
    spectrum = read_hotfile(source, nodes=[node])
    frequency, directions, variance = spectrum["frequency"], spectrum["direction"], spectrum["variance"][0]
    rad = np.radians(directions)
    weights = variance.sum(axis=0)
    original_direction = math.degrees(math.atan2(np.sum(weights * np.sin(rad)), np.sum(weights * np.cos(rad))))
    toward = (90 - waves["wind_toward_compass_deg"]) % 360
    turned = (directions + toward - original_direction) % 360
    order = np.argsort(turned)
    depth, g = waves["depth_m"], MOON_SURFACE_GRAVITY
    fp = float(frequency[np.argmax(variance.sum(axis=1))])
    kp = float(wavenumber(np.array([fp]), depth, g)[0])
    k_from = float(wavenumber(np.array([frequency[-1]]), depth, g)[0])
    densities = [realization.swan_density(frequency, turned[order], variance[:, order], depth, g),
                 realization.short_density(waves["friction_velocity_m_s"], toward, g, k_from, 2 * math.pi * fp / kp)]
    # Smaller tiles suffice for a structure guide; no full spectral render runs.
    tiles = realization.realize(densities, depth, g, 20380207,
                                cascades=((512.0, 512), (31.7, 512), (1.9, 512)))
    pyramids = [realization.lean_pyramid(t.slope) for t in tiles]
    eye = scene["camera"]["eye_height_m"]
    horizon = scene["camera"]["horizon_row_at_centre"]
    focal = width / (2 * camera.th)
    distance = np.geomspace(2.0, 3000.0, 6000)
    footprint = np.maximum(distance / focal, distance**2 / (eye * focal))
    ocean = pixels.copy()
    index = float(reflection.refractive_index(550.0))

    for x in range(width):
        ray = camera.direction(np.array([x + 0.5]), np.array([horizon]))[0]
        azimuth = math.atan2(ray[0], ray[1])
        xx, yy = distance * math.sin(azimuth), distance * math.cos(azimuth)
        tile = tiles[0]
        z = map_coordinates(tile.height, [yy / tile.spacing_m, xx / tile.spacing_m], order=1, mode="grid-wrap")
        sx, sy = np.zeros_like(distance), np.zeros_like(distance)
        for tile, pyramid in zip(tiles, pyramids):
            level = np.clip(np.floor(np.log2(np.maximum(1, footprint / tile.spacing_m))).astype(int), 0, len(pyramid) - 1)
            for li in np.unique(level):
                selected = level == li
                cell = tile.spacing_m * 2**int(li)
                coords = [yy[selected] / cell, xx[selected] / cell]
                sx[selected] += map_coordinates(pyramid[li][0], coords, order=1, mode="grid-wrap")
                sy[selected] += map_coordinates(pyramid[li][1], coords, order=1, mode="grid-wrap")
        elevation = np.degrees(np.arctan2(z - eye - distance**2 / (2 * MOON_RADIUS), distance))
        _, rows = camera.project(unit(elevation, np.full_like(elevation, math.degrees(azimuth))))
        length = np.sqrt(distance**2 + (eye - z)**2)
        vx, vy, vz = -xx / length, -yy / length, (eye - z) / length
        cosine = np.maximum(0.0, (-sx * vx - sy * vy + vz) / np.sqrt(1 + sx**2 + sy**2))
        # A facet contrast guide only: clipped to avoid implying unmodelled
        # direct solar glints. The source's ensemble sea colour sets the mean.
        contrast = np.clip(reflection.fresnel(cosine, index) / np.maximum(reflection.fresnel(vz, index), 1e-6), 0.3, 2.0)
        top = height
        for j in range(len(distance)):
            row = max(math.ceil(horizon), int(math.floor(rows[j])))
            if row < top:
                ocean[row:top, x] = sea[row:top, x] * contrast[j]
                top = row
        # Geometry above the flat horizon is intentionally left to a full renderer.
        ocean[region[:, x] != 2, x] = pixels[region[:, x] != 2, x]

    Image.fromarray(np.uint8(np.clip(np.rint(encoded(ocean) * 255), 0, 255))).save(out / "ic-1-wave-structure-guide.png")
    record = {
        "schema":"terluna.visualization.sea-wave-guide/1",
        "inputs":{str(p.relative_to(root)):digest(p) for p in (source, month_path)},
        "producer_sha256":digest(__file__),
        "evidence":"Domain wave realization from the scene's borrowed Nectaris spectrum, turned to the stored local wind. Projected height and filtered slopes guide water structure; the Fresnel contrast treatment is approximate and is not a radiance solution.",
        "reading_rule":"Use wave scales and perspective as an image-generation reference, with the separate calculated colour guide. No validation of a lunar sea; this scene borrows an Earth-law wave model's lunar spectrum from another sea.",
        "seed":20380207,
        "source_hs_m":waves["significant_height_m"],
        "guide_geometry_hs_m":float(4 * np.std(tiles[0].height)),
        "source_peak_wavelength_m":waves["peak_wavelength_m"],
        "guide_peak_wavelength_m":float(2 * math.pi / kp),
        "source_mss":waves["mean_square_slope"]["total"],
        "drawn_mss":float(sum(np.trace(realization.sample_covariance(t)) for t in tiles)),
        "cascades":[[t.size_m,t.cells] for t in tiles],
        "omissions":["full radiance transport","statistical residual slope below the last tile","facet shadowing","terrain reflections","height variations of short-wave tiles","wave visibility above the flat horizon"],
    }
    (out / "ic-1-wave-structure-guide.json").write_text(json.dumps(record,indent=2)+"\n")
    print(json.dumps({k:record[k] for k in ("source_hs_m","guide_geometry_hs_m","source_peak_wavelength_m","guide_peak_wavelength_m","source_mss","drawn_mss")},indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--out", type=Path)
    parser.add_argument("--water-structure", action="store_true")
    args = parser.parse_args()
    root = args.root.resolve()
    sys.path.insert(0, str(root))
    from research.studies.sea_appearance.scenes import Frame, unit
    from shared.constants import MOON_RADIUS

    source = root / "research/studies/sea_appearance/results/scenes.json"
    terrain = root / "research/runs/sea_appearance/terrain/ingenii_coast.npz"
    scene = next(s for s in json.loads(source.read_text())["scenes"]
                 if s["coast"] == "Ingenii coast" and s["moment"] == "noon")
    out = args.out or root / "visualization/sea-appearance/results/illustrations"
    out.mkdir(parents=True, exist_ok=True)
    width, height = scene["camera"]["frame"]
    horizon = scene["camera"]["horizon_row_at_centre"]
    camera = Frame(scene["camera"]["view_azimuth_deg"], scene["camera"]["pitch_deg"])
    eye = scene["camera"]["eye_height_m"]
    sky = field(scene["sky"], width, height)
    sea = field(scene["sea"], width, height)
    pixels = sky.copy()
    pixels[math.ceil(horizon):] = sea[math.ceil(horizon):]
    region = np.zeros((height, width), np.uint8)
    region[math.ceil(horizon):] = 2
    with np.load(terrain, allow_pickle=False) as z:
        heights = z["height_m"].astype(float) - scene["tide_m"]
        axis = z["x_m"]
    spacing = float(axis[1] - axis[0])
    gy, gx = np.gradient(heights, spacing)
    sun = unit(scene["sun"]["elevation_deg"], scene["sun"]["azimuth_deg"])
    direct_horizontal = scene["sun"]["normal_illuminance_lux"] * sun[2]
    direct_share = min(1.0, direct_horizontal / scene["light"]["ground_lux"])
    apparent_land = linear(rgb(scene["land"][0]["apparent_srgb"]))
    horizon_point = next(p for p in scene["sky"] if p["label"] == "just above the horizon, centre")
    haze = linear(rgb(horizon_point["srgb"]))
    land_reference = scene["land"][0]
    reference_distance = land_reference["highest_point_distance_km"] * 1000
    reference_trans = land_reference["haze_transmittance_550nm"]
    intrinsic = np.maximum(0, (apparent_land - (1 - reference_trans) * haze) / reference_trans)
    d = np.arange(50.0, 105000.0, 50.0)
    skyline_y = np.full(width, horizon)

    # One distance profile per screen column: a lightweight terrain projection.
    # Pitch changes the azimuth slightly up the column; this approximation is
    # measured below against the scene's stored skyline samples.
    for x in range(width):
        ray = camera.direction(np.array([x + 0.5]), np.array([horizon]))[0]
        azimuth = math.degrees(math.atan2(ray[0], ray[1])) % 360
        a = math.radians(azimuth)
        xx, yy = d * math.sin(a), d * math.cos(a)
        ix, iy = (xx - axis[0]) / spacing, (yy - axis[0]) / spacing
        inside = (ix >= 0) & (iy >= 0) & (ix <= len(axis) - 1) & (iy <= len(axis) - 1)
        dd = d[inside]
        coords = [iy[inside], ix[inside]]
        z = map_coordinates(heights, coords, order=1, mode="nearest")
        slopes_x = map_coordinates(gx, coords, order=1, mode="nearest")
        slopes_y = map_coordinates(gy, coords, order=1, mode="nearest")
        elevations = np.degrees(np.arctan2(z - eye - dd * dd / (2 * MOON_RADIUS), dd))
        _, screen_y = camera.project(unit(elevations, np.full_like(elevations, azimuth)))
        cosine = (-slopes_x * sun[0] - slopes_y * sun[1] + sun[2]) / np.sqrt(1 + slopes_x**2 + slopes_y**2)
        # Isotropic diffuse fill and normal-dependent direct illumination, used
        # only to show orientation. This is not a solved local irradiance field.
        shade = (1 - direct_share) + direct_share * np.maximum(cosine, 0) / sun[2]
        transmit = reference_trans ** (dd / reference_distance)
        colours = intrinsic[None] * shade[:, None] * transmit[:, None] + haze[None] * (1 - transmit[:, None])
        top = math.ceil(horizon)
        for j in range(len(dd)):
            row = max(0, int(math.floor(screen_y[j])))
            if z[j] > 0 and row < top:
                pixels[row:top, x] = colours[j]
                region[row:top, x] = 1
                top = row
        skyline_y[x] = top

    image = np.uint8(np.clip(np.rint(encoded(pixels) * 255), 0, 255))
    Image.fromarray(image).save(out / "ic-1-calculated-guide.png")
    # Component guide: the sky and true terrain over the untextured sea's colour field.
    Image.fromarray(np.where(region == 1, 255, 0).astype(np.uint8)).save(out / "ic-1-land-mask.png")
    errors = []
    for p in scene["land_profile"]:
        x = int(p["x"])
        expected = horizon - p["px_above_horizon"]
        errors.append({"x":x,"expected_y":round(expected, 2),"guide_y":float(skyline_y[x]),
                       "difference_px":round(float(skyline_y[x] - expected), 2)})
    record = {
        "schema":"terluna.visualization.sea-image-guide/1",
        "scene_id":"IC-1",
        "producer":{"script":"visualization/sea-appearance/image_guides.py", "sha256":digest(__file__)},
        "inputs":{str(p.relative_to(root)):digest(p) for p in (source, terrain)},
        "evidence":"Projection of the stored LOLA terrain and interpolation of the scene's display colours. Terrain illumination is an illustrative normal-based approximation, not the spectral path tracer.",
        "reading_rule":"Use as geometry, colour and illumination-direction guidance for image generation; retain this provenance outside the final image. Individual waves, cast shadows, interreflection and terrain reflections are absent.",
        "width":width,"height":height,"horizon_y":horizon,
        "region_fractions":{"sky":float(np.mean(region==0)),"land":float(np.mean(region==1)),"sea":float(np.mean(region==2))},
        "skyline_comparison":errors,
        "maximum_skyline_difference_px":max(abs(p["difference_px"]) for p in errors),
        "sun_elevation_deg":scene["sun"]["elevation_deg"],"sun_azimuth_deg":scene["sun"]["azimuth_deg"],
        "direct_horizontal_light_share":float(direct_share),
    }
    (out / "ic-1-calculated-guide.json").write_text(json.dumps(record,indent=2)+"\n")
    print(json.dumps({k:record[k] for k in ("region_fractions","maximum_skyline_difference_px","direct_horizontal_light_share")},indent=2))
    if args.water_structure:
        water_structure(root, scene, camera, pixels, sea, region, out)


if __name__ == "__main__":
    main()
