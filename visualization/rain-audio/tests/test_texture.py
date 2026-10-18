"""The texture models are the recovered ones, carry their CC0 sources, reproduce their fitted statistics, and render
the rain work's recorded samples bit for bit.

    python3 -m pytest visualization/rain-audio/tests
"""
import hashlib
import json

import numpy as np
import pytest

import texture

MANIFEST = json.loads((texture.MODELS / 'models.json').read_text())
INPUTS = json.loads((texture.HERE / 'inputs.json').read_text())
ROWS = {row['surface']: row for row in MANIFEST['models']}


def samples_sha256(x):
    return hashlib.sha256(np.ascontiguousarray(x, dtype='<f4').tobytes()).hexdigest()


def test_the_models_are_the_recorded_files_with_cc0_sources():
    assert MANIFEST['schema'] == 'terluna.visualization.rain-texture-models/1'
    assert set(ROWS) == set(texture.SURFACES)
    recordings = {f['name']: f for f in INPUTS['files']}
    for surface, row in ROWS.items():
        path = texture.model_path(surface)
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row['sha256'] and path.stat().st_size == row['bytes']
        model = texture.load(surface)
        assert {k: list(np.asarray(v).shape) for k, v in model.items()} == row['arrays']
        source = row['source']
        assert str(model['source_sha256']) == source['sha256'] == recordings[source['file']]['sha256']
        assert source['licence'] == recordings[source['file']]['licence'] == 'CC0 1.0'
        assert [float(model['training_start_seconds']),
                float(model['training_start_seconds'] + model['training_duration_seconds'])] == row['fit']['interval_s']
        assert int(model['renderer_version']) == 3 and int(model['bands']) == 64


@pytest.mark.parametrize('surface', texture.SURFACES)
def test_generated_envelopes_reproduce_the_fitted_quantiles(surface):
    model = texture.load(surface)
    env = texture.generate(model, 60.0, 7)
    knots = np.asarray(model['latent_knots'], dtype=float)
    quantiles = np.asarray(model['amplitude_quantiles'], dtype=float)
    for z, p in ((-1.2816, 0.1), (0.0, 0.5), (1.2816, 0.9)):
        ratio = np.quantile(env, p, axis=1) / quantiles[:, int(np.argmin(abs(knots - z)))]
        assert np.all((ratio > 0.9) & (ratio < 1.1)), (surface, p)
        assert abs(np.median(ratio) - 1) < 0.03


@pytest.mark.parametrize('surface', texture.SURFACES)
def test_a_short_render_reproduces_the_recovered_renderers_samples(surface):
    row = ROWS[surface]['recovered_renderer_short_render']
    model = texture.load(surface)
    guarded, _ = texture.synthesize_guarded(model, row['seconds'], row['seed'])
    assert samples_sha256(texture.crop(guarded, row['seconds'])) == row['samples_sha256']
    if surface == 'concrete':
        # Another seed makes another rain, and the stereo pair's first channel is the same render before the clip
        # is cut from its guards.
        assert not np.array_equal(texture.synthesize(model, row['seconds'], row['seed'] + 1),
                                  texture.crop(guarded, row['seconds']))
        pair, _ = texture.carrier_pair(model, row['seconds'], row['seed'])
        assert np.array_equal(pair[:, 0], guarded) and not np.array_equal(pair[:, 1], guarded)


def test_the_concrete_master_reproduces_the_archive_bit_for_bit():
    row = ROWS['concrete']
    master = row['archive_master']
    model = texture.load('concrete')
    x = texture.master(model, master['seconds'], master['seed'])
    assert samples_sha256(x) == master['samples_sha256']
    assert float(np.max(np.abs(x))) == master['peak']
    assert 20 * np.log10(np.sqrt(np.mean(x.astype(np.float64)**2))) == pytest.approx(master['rms_dbfs'], abs=1e-9)
    # A render's level follows the recording's level over the fitted interval, within 1.5 dB.
    assert abs(master['rms_dbfs'] - 20 * np.log10(float(model['reference_rms']))) < 1.5


def test_refitting_from_the_recording_returns_the_model():
    recording = texture.HERE / 'inputs' / ROWS['concrete']['source']['file']
    if not recording.is_file():
        pytest.skip('the recordings are not restored; run fetch_inputs.py')
    fitted = texture.train(recording, *[ROWS['concrete']['fit']['interval_s'][0],
                                        np.diff(ROWS['concrete']['fit']['interval_s'])[0]])
    model = texture.load('concrete')
    for key in ('amplitude_quantiles', 'modulation_cholesky', 'radial_reference_quantiles', 'radial_pilot_quantiles',
                'mean_amplitude', 'reference_rms'):
        assert np.array_equal(np.asarray(fitted[key]), np.asarray(model[key])), key
