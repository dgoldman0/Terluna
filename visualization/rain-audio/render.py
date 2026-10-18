"""Render the rain's sound in the three target climates beside Earth's, and the thunder, with a manifest.

    python3 visualization/rain-audio/render.py --output <folder>    # about four minutes on one thread

For each surface it renders one Earth texture, mono and as a stereo pair. For each target climate it carries the
texture to the Moon at two of the GCM's rain rates on that target's land (scenes.gcm_rates): a median rainy 6.5 hours
and a heavy one. Each uses the impact-force transfer from the Earth scene, 5 mm/h in Earth's air of the target's
temperature and humidity, to the Moon scene (scenes.py, transfer.py), and reaches the microphones through the target's
air (propagation.py). From those renders it writes, per target:

    <surface>-earth-moon-matched.wav           mono: Earth, the Moon's median rain, its heavy rain, 12 s each, each
                                               at -35.4 LUFS
    <surface>-earth-moon-shared-gain.wav       mono: the whole renders with one gain, the model's relative levels
    <surface>-stereo-earth-moon-matched.wav    the same in stereo
    <surface>-stereo-earth-moon-shared-gain.wav
    concrete-large-drops-shared-gain.wav       the Moon's heavy rain, then the same with 1% of its water in drops
                                               of 8-12 mm
    grass-leaves-distance.wav                  one 1 m2 patch of heavy rain near, middle and far, Earth then Moon

and, across the targets, <surface>-targets-shared-gain.wav (the Earth texture, then the Moon's median rain in the hot
world, the resort and the lock, then their heavy rain, 12 s each with one gain) and thunder-earth-moon-shared-gain.wav
and -matched.wav (four events below their flashes, with one gain, and each at -23 LUFS). A gain that would lift a true
peak above -1 dBFS stops there. Raw renders are float WAV at the texture model's digital level. Every level is
digital: none is a sound pressure in the world except the thunder's, whose calibration the manifest records. Rendered
audio stays out of Git; the research drive keeps it.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for path in (str(HERE), str(ROOT)):
    if path not in sys.path:
        sys.path.insert(0, path)

import levels  # noqa: E402
import propagation  # noqa: E402
import scenes  # noqa: E402
import texture  # noqa: E402
import thunder_layer  # noqa: E402
import transfer  # noqa: E402

SECONDS = 30.0
PASSAGE_SECONDS = 12.0
DISTANCE_SECONDS = 6.0
RATES = ('median', 'heavy')
MONO_SEEDS = {'concrete': 4101, 'puddles': 5101, 'grass-leaves': 6101}       # the rain work's lunar renders
STEREO_SEEDS = {'concrete': 7201, 'puddles': 8201, 'grass-leaves': 9201}     # and its stereo renders
TAIL_FRACTION = 0.01
TAIL_AREA_M2 = 4.0
# One rain patch of 1 m2 at three places: near-left, middle-front and far-right (x right, y forward, z up, m).
DISTANCE_POINTS = (('near-left', (-1.0, 1.0, 0.0)), ('middle-front', (0.0, 5.0, 0.0)),
                   ('far-right', (12.0, 16.0, 0.0)))
THUNDER_ORDER = ('earth-cloud-flash', 'earth-ground-strike', 'moon-cloud-flash', 'moon-ground-strike')
THUNDER_AZIMUTH_DEG = 30.0
THUNDER_MATCHED_LUFS = -23.0
CODE = ('render.py', 'texture.py', 'scenes.py', 'transfer.py', 'propagation.py', 'thunder_layer.py', 'levels.py')


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class Writer:
    """Writes files under one folder and records each one."""

    def __init__(self, root: Path):
        self.root = Path(root)
        self.files = []

    def wav(self, name, x, bits=24, **record):
        path = levels.write_wav(self.root / name, x, bits=bits)
        self.files.append(dict(file=name, sha256=sha256(path), bits=bits, **levels.describe(x), **record))
        return path


def matched(passages, target=levels.TARGET_LUFS, ceiling_dbfs=-1.0):
    """Each passage faded and brought to the target loudness on its own, unless that would lift its true peak above
    the ceiling, when its gain stops there; returns the audio, gains and spans."""
    faded = [levels.fade(p) for p in passages]
    gains = [min(levels.gain_to(p, target), ceiling_dbfs - 20 * np.log10(levels.true_peak(p))) for p in faded]
    audio, spans = levels.sequence([p * 10**(g / 20) for p, g in zip(faded, gains)])
    return audio, gains, spans


def shared(passages, reference=0, target=levels.TARGET_LUFS, ceiling_dbfs=-1.0):
    """Every passage faded and scaled by one gain: the gain that brings passage `reference` to the target loudness,
    or less if that would lift any passage's true peak above the ceiling."""
    faded = [levels.fade(p) for p in passages]
    gain = min(levels.gain_to(faded[reference], target),
               ceiling_dbfs - 20 * np.log10(max(levels.true_peak(p) for p in faded)))
    audio, spans = levels.sequence([p * 10**(gain / 20) for p in faded])
    return audio, gain, spans


def compare(earth, moon):
    """The Moon's level against Earth's, by root-mean-square and by loudness (dB)."""
    return dict(rms=levels.rms_dbfs(moon) - levels.rms_dbfs(earth),
                loudness=levels.loudness(moon) - levels.loudness(earth))


def render_surface(surface, product, rates, writer, targets, seconds, spectrum, earth_mm_h, record):
    model = texture.load(surface)
    guarded, mono_source = texture.synthesize_guarded(model, seconds, MONO_SEEDS[surface])
    pair, stereo_source = texture.carrier_pair(model, seconds, STEREO_SEEDS[surface])
    earth_mono = texture.crop(guarded, seconds)
    writer.wav(f'raw/{surface}-earth-mono.wav', earth_mono, bits=32, kind='raw')
    entry = dict(model_sha256=sha256(texture.model_path(surface)), mono_source=mono_source,
                 stereo_source=stereo_source, targets={})
    moons = {}
    grid = propagation.FIELD_FREQUENCIES
    n = round(PASSAGE_SECONDS * levels.SAMPLE_RATE)
    names = ['earth'] + [f'moon, {r} rain' for r in RATES]
    for target in targets:
        rate = scenes.gcm_rates(rates, target)
        earth = scenes.scene(product, target, 'earth', surface, earth_mm_h, spectrum)
        cov_earth = propagation.field_covariance(grid, earth.air)
        # The source is calibrated once, to Earth's field at the nominal geometry, and kept for the Moon, so the field
        # keeps the texture's spectrum at each microphone on Earth and the Moon's differs by its physics alone.
        calibration = 1 / np.sqrt((cov_earth['power_left'] + cov_earth['power_right']) / 2)
        earth_stereo = texture.crop(propagation.apply_matrix(
            pair, grid, propagation.covariance_factor(cov_earth, calibration)), seconds)
        writer.wav(f'raw/{target}-{surface}-earth-stereo.wav', earth_stereo, bits=32, kind='raw')
        result = dict(earth=earth.describe(), rates_mm_h=rate, by_rate={},
                      field=dict(earth=cov_earth['metadata']))
        monos, stereos, magnitudes, moon_scenes = [earth_mono], [earth_stereo], {}, {}
        cov_moon = None
        for name in RATES:
            moon = scenes.scene(product, target, 'moon', surface, rate[name], spectrum)
            magnitude, transfer_meta = transfer.make_transfer(earth, moon)
            filtered, filter_meta = transfer.apply_transfer(guarded, transfer.FREQUENCIES, magnitude)
            moon_mono = texture.crop(filtered, seconds)
            if cov_moon is None:
                cov_moon = propagation.field_covariance(grid, moon.air)
                result['field']['moon'] = cov_moon['metadata']
            lunar_pair = np.column_stack([transfer.apply_transfer(pair[:, c], transfer.FREQUENCIES, magnitude)[0]
                                          for c in range(2)])
            moon_stereo = texture.crop(propagation.apply_matrix(
                lunar_pair, grid, propagation.covariance_factor(cov_moon, calibration)), seconds)
            writer.wav(f'raw/{target}-{surface}-moon-{name}-mono.wav', moon_mono, bits=32, kind='raw')
            writer.wav(f'raw/{target}-{surface}-moon-{name}-stereo.wav', moon_stereo, bits=32, kind='raw')
            monos.append(moon_mono)
            stereos.append(moon_stereo)
            magnitudes[name], moon_scenes[name] = magnitude, moon
            moons[(target, name)] = moon_mono
            result['by_rate'][name] = dict(rain_mm_h=rate[name], moon=moon.describe(), transfer=transfer_meta,
                                           filter=filter_meta, moon_minus_earth_db=dict(
                                               mono=compare(earth_mono, moon_mono),
                                               stereo=compare(earth_stereo, moon_stereo)))
        folder = f'{target}/'
        audio, gains, spans = matched([m[:n] for m in monos])
        writer.wav(folder + f'{surface}-earth-moon-matched.wav', audio, kind='audition', passages=names,
                   gains_db=gains, spans_s=spans, rates_mm_h=rate)
        audio, gain, spans = shared(monos)
        writer.wav(folder + f'{surface}-earth-moon-shared-gain.wav', audio, kind='audition', passages=names,
                   gain_db=gain, spans_s=spans, rates_mm_h=rate)
        audio, gains, spans = matched([m[:n] for m in stereos])
        writer.wav(folder + f'{surface}-stereo-earth-moon-matched.wav', audio, kind='audition', passages=names,
                   gains_db=gains, spans_s=spans, rates_mm_h=rate)
        audio, gain, spans = shared(stereos)
        writer.wav(folder + f'{surface}-stereo-earth-moon-shared-gain.wav', audio, kind='audition', passages=names,
                   gain_db=gain, spans_s=spans, rates_mm_h=rate)
        if surface == 'concrete':
            heavy = scenes.scene(product, target, 'moon', surface, rate['heavy'], spectrum, tail_fraction=TAIL_FRACTION)
            heavy_magnitude, _ = transfer.make_transfer(earth, heavy)
            bed = texture.crop(transfer.apply_transfer(guarded, transfer.FREQUENCIES, heavy_magnitude)[0], seconds)
            events, tail_meta = transfer.large_drop_events(earth_mono, earth, heavy, seconds,
                                                           MONO_SEEDS[surface] + 6000, TAIL_AREA_M2)
            audio, gain, spans = shared([monos[RATES.index('heavy') + 1], bed + events])
            writer.wav(folder + 'concrete-large-drops-shared-gain.wav', audio, kind='audition',
                       passages=['moon, heavy rain', 'moon, heavy rain with 1% of the water in 8-12 mm drops'],
                       gain_db=gain, spans_s=spans, rain_mm_h=rate['heavy'])
            result['large_drops'] = dict(scene=heavy.describe(), events=tail_meta)
        if surface == 'grass-leaves':
            source = propagation.apply_matrix(pair, grid, np.eye(2)[None, :, :] * calibration[:, None, None])[:, 0]
            lunar_source = transfer.apply_transfer(source, transfer.FREQUENCIES, magnitudes['heavy'])[0]
            passages, labels, demo = [], [], []
            d = round(DISTANCE_SECONDS * levels.SAMPLE_RATE)
            for world, src, air in (('earth', source, earth.air), ('moon', lunar_source, moon_scenes['heavy'].air)):
                for name, xyz in DISTANCE_POINTS:
                    sound, meta = propagation.propagate_mono(src, xyz, air, source_area_m2=1.0)
                    passages.append(texture.crop(sound, seconds)[:d])
                    labels.append(f'{world} {name}')
                    demo.append(dict(world=world, point=name, distances_m=meta['distances_m'],
                                     absorption_db_at_1_4_10_khz=[f['absorption_db_at_1_4_10_khz']
                                                                  for f in meta['absorption_filters']]))
            faded = [levels.fade(p) for p in passages]
            gains = []
            for k in (0, 3):                      # one gain per world, from its near passage
                gains += [levels.gain_to(faded[k])] * 3
            audio, spans = levels.sequence([p * 10**(g / 20) for p, g in zip(faded, gains)])
            writer.wav(folder + 'grass-leaves-distance.wav', audio, kind='audition', passages=labels,
                       gains_db=gains, spans_s=spans, moon_rain_mm_h=rate['heavy'])
            result['distance'] = demo
        entry['targets'][target] = result
    # Each target pairs the Moon with Earth's air at the same temperature and humidity, as the rain product does, so
    # the texture stands for a slightly different Earth rain in each; this records by how much.
    probes = np.array([100.0, 1000.0, 4000.0, 10000.0])
    earth_q = {t: transfer.source_spectrum(scenes.scene(product, t, 'earth', surface, earth_mm_h, spectrum), probes)
               for t in targets}
    first = earth_q[targets[0]]
    entry['earth_reference_db_against_first_target'] = dict(
        frequency_hz=probes.tolist(), first_target=targets[0],
        **{t: (10 * np.log10(q / first)).tolist() for t, q in earth_q.items()})
    if len(targets) > 1:
        order = [(t, r) for r in RATES for t in targets]
        audio, gain, spans = shared([earth_mono[:n]] + [moons[k][:n] for k in order])
        writer.wav(f'{surface}-targets-shared-gain.wav', audio, kind='audition',
                   passages=['earth'] + [f'moon, {t}, {r} rain' for t, r in order], gain_db=gain, spans_s=spans)
    record[surface] = entry


def render_thunder(writer, sound_speed_m_s):
    product = thunder_layer.load()
    pressures, records = [], []
    for name in THUNDER_ORDER:
        p, meta = thunder_layer.render(product, name)
        pressures.append(thunder_layer.stereo(p, THUNDER_AZIMUTH_DEG, sound_speed_m_s))
        records.append(meta)
    # One gain for all four, bringing the loudest true peak to -1 dBFS: the digital level stands for pressure.
    peak = max(levels.true_peak(p) for p in pressures)
    pascal_per_full_scale = peak / 10**(-1 / 20)
    faded = [levels.fade(p / pascal_per_full_scale) for p in pressures]
    audio, spans = levels.sequence(faded, gap_seconds=1.5)
    writer.wav('thunder-earth-moon-shared-gain.wav', audio, kind='audition', passages=list(THUNDER_ORDER),
               spans_s=spans, full_scale_pa=pascal_per_full_scale,
               full_scale_db_spl=20 * np.log10(pascal_per_full_scale / thunder_layer.P_REF))
    audio, gains, spans = matched([p / pascal_per_full_scale for p in pressures], THUNDER_MATCHED_LUFS)
    writer.wav('thunder-earth-moon-matched.wav', audio, kind='audition', passages=list(THUNDER_ORDER),
               gains_db=gains, spans_s=spans)
    return dict(product=str(thunder_layer.PRODUCT.relative_to(ROOT)), product_sha256=thunder_layer.sha256(),
                azimuth_deg=THUNDER_AZIMUTH_DEG, full_scale_pa=pascal_per_full_scale, events=records)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--output', type=Path, default=HERE / 'renders')
    parser.add_argument('--seconds', type=float, default=SECONDS)
    parser.add_argument('--earth-rain-mm-h', type=float, default=scenes.EARTH_REFERENCE_MM_H,
                        help='the Earth rain the texture stands for')
    parser.add_argument('--spectrum', type=float, nargs=2, metavar=('DM_MM', 'SHAPE'))
    parser.add_argument('--targets', nargs='+', default=list(scenes.TARGETS), choices=scenes.TARGETS)
    parser.add_argument('--surfaces', nargs='+', default=list(texture.SURFACES), choices=texture.SURFACES)
    parser.add_argument('--no-thunder', action='store_true')
    args = parser.parse_args(argv)
    if not PASSAGE_SECONDS <= args.seconds <= 120:
        parser.error('Use 12-120 seconds')
    start = time.monotonic()
    product, rates = scenes.load_arrivals(), scenes.load_rates()
    writer = Writer(args.output)
    record = {}
    for surface in args.surfaces:
        render_surface(surface, product, rates, writer, args.targets, args.seconds,
                       tuple(args.spectrum) if args.spectrum else None, args.earth_rain_mm_h, record)
        print(f'{surface}: {time.monotonic() - start:.0f} s', flush=True)
    thunder_record = None if args.no_thunder else render_thunder(
        writer, scenes.air(product, 'resort', 'moon')['sound_speed_m_s'])
    manifest = dict(
        schema='terluna.visualization.rain-audio-renders/2',
        evidence=('Rendered auditions: the Earth texture models of CC0 recordings, carried to the Moon at the GCM\'s '
                  'rain rates by an assumed impact-force transfer and propagated by direct paths. Levels are digital '
                  'and relative except the thunder\'s; no lunar rain sound has been measured, and listening judges '
                  'these renders.'),
        settings=dict(seconds=args.seconds, passage_seconds=PASSAGE_SECONDS, earth_rain_mm_h=args.earth_rain_mm_h,
                      moon_rain_mm_h={t: scenes.gcm_rates(rates, t) for t in args.targets},
                      moon_rates=dict(median='the median wet 6.5-hour mean on the land',
                                      heavy='the 99th percentile of the wet 6.5-hour means on the land'),
                      spectrum=args.spectrum or [product['central_spectrum']['dm_mm'],
                                                 product['central_spectrum']['shape_mu']],
                      targets=args.targets, surfaces=args.surfaces, target_lufs=levels.TARGET_LUFS,
                      spacing_m=propagation.SPACING_M, listener_height_m=propagation.LISTENER_HEIGHT_M,
                      radius_edges_m=list(propagation.RADIUS_EDGES_M), tail_fraction=TAIL_FRACTION,
                      tail_area_m2=TAIL_AREA_M2),
        inputs={str(scenes.ARRIVALS.relative_to(ROOT)): sha256(scenes.ARRIVALS),
                str(scenes.RATES.relative_to(ROOT)): sha256(scenes.RATES),
                **({str(thunder_layer.PRODUCT.relative_to(ROOT)): thunder_layer.sha256()} if thunder_record else {}),
                **{f'visualization/rain-audio/models/{s}.npz': sha256(texture.model_path(s)) for s in args.surfaces}},
        code={f'visualization/rain-audio/{c}': sha256(HERE / c) for c in CODE},
        environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
        surfaces=record, thunder=thunder_record, files=writer.files, elapsed_s=time.monotonic() - start)
    (Path(args.output) / 'manifest.json').write_text(json.dumps(manifest, indent=1, allow_nan=False) + '\n')
    print(f'wrote {len(writer.files)} files and manifest.json to {args.output} in {time.monotonic() - start:.0f} s')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
