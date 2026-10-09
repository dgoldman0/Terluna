"""Spectral boundaries and inherited-transfer checks, not biological validation."""
import importlib
import json

import numpy as np
import pytest

from shared.provenance import constants_changed, constants_used

model = importlib.import_module('research.studies.uv_biology.run')


def test_uv_window_has_absolute_transmission_and_preserves_protected_wavelengths():
    wavelengths = np.arange(202., 501.)
    baseline = np.full(wavelengths.shape, 0.001)
    unfiltered = np.full(wavelengths.shape, 2.)
    case = model.window_spectrum(wavelengths, baseline, unfiltered, (295, 315, 0.01))
    band = (wavelengths >= 295) & (wavelengths < 315)
    np.testing.assert_array_equal(case[~band], baseline[~band])
    np.testing.assert_allclose(case[band], 0.02)
    np.testing.assert_array_equal(baseline, np.full(wavelengths.shape, 0.001))


def test_reference_window_and_zero_flux():
    wavelengths = np.arange(202., 501.)
    baseline = np.linspace(0., 1., len(wavelengths))
    reference = model.window_spectrum(wavelengths, baseline, baseline * 100, None)
    np.testing.assert_array_equal(reference, baseline)
    metrics = model.metrics(wavelengths, np.zeros_like(wavelengths))
    assert metrics['uv_index'] == 0
    assert metrics['uv_energy_w_m2'] == 0
    assert all(value == 0 for value in metrics['band_photons_umol_m2_s'].values())
    with pytest.raises(ValueError):
        model.metrics(wavelengths[::2], baseline[::2])


def test_metric_scaling_and_nonoverlapping_uv_bands():
    wavelengths = np.arange(202., 501.)
    unit = model.metrics(wavelengths, np.ones_like(wavelengths))
    double = model.metrics(wavelengths, 2 * np.ones_like(wavelengths))
    assert unit['uv_energy_w_m2'] == 120
    assert double['uv_index'] == pytest.approx(2 * unit['uv_index'])
    assert double['bee_uv_relative_catch'] == pytest.approx(2 * unit['bee_uv_relative_catch'])


def test_runner_keeps_every_shortwave_bin_and_records_no_adoption(tmp_path, monkeypatch):
    monkeypatch.setattr(model, 'OUT', tmp_path / 'uv.json')
    product = model.run()
    saved = json.loads(model.OUT.read_text())
    assert saved['schema'] == product['schema']
    assert product['adoption']['current_spectrum_retained'] is True
    assert product['adoption']['added_uv_adopted'] is False
    for row in product['by_sun_degrees'].values():
        assert all(case['unchanged_below_295nm'] for case in row['cases'].values())
        narrow = row['cases']['long_uvb_310_315_10pct']['uv_index']
        wide = row['cases']['uvb_295_315_10pct']['uv_index']
        assert narrow < wide
        assert row['cases']['reference']['added_uv_energy_w_m2'] == 0
    for path, digest in product['producer']['inputs'].items():
        assert model.digest(model.ROOT / path) == digest


def test_saved_product_pins_current_dependencies_and_reproduces(tmp_path, monkeypatch):
    saved = json.loads((model.HERE / 'results/uv_biology.json').read_text())
    provenance = saved['producer']
    for group in ('files', 'inputs'):
        for path, recorded_digest in provenance[group].items():
            assert model.digest(model.ROOT / path) == recorded_digest, path
    assert not constants_changed(provenance['constants'])
    assert provenance['constants'] == constants_used([
        model.ROOT / path for path in provenance['files']
    ])
    monkeypatch.setattr(model, 'OUT', tmp_path / 'regenerated_uv.json')
    model.run()
    assert json.loads(model.OUT.read_text()) == saved
