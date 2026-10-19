"""The thunder layer renders the thunder product's levels, span and audible bands."""
import numpy as np
import pytest
from scipy import signal

import thunder_layer

PRODUCT = thunder_layer.load()


@pytest.mark.parametrize('name', thunder_layer.EVENTS)
def test_each_band_has_the_products_loudest_second(name):
    pressure, meta = thunder_layer.render(PRODUCT, name)
    levels = thunder_layer.event(PRODUCT, name)
    for band, record in meta['bands'].items():
        assert record['rendered_db'] == pytest.approx(record['product_db'], abs=0.01)
    rendered = {float(b) for b in meta['bands']}
    for band, level in levels['peak_db_by_band'].items():
        threshold = dict(zip(thunder_layer.thunder.BANDS, thunder_layer.thunder.THRESHOLD_DB))[float(band)]
        assert (float(band) in rendered) == (np.isnan(threshold) or level >= threshold)
    # Filtered out of the whole event, each strong band keeps its level within half a decibel.
    for band in (31.5, 63.0, 125.0):
        sos = signal.butter(4, [band / np.sqrt(2), band * np.sqrt(2)], 'bandpass', fs=thunder_layer.SAMPLE_RATE,
                            output='sos')
        level = 20 * np.log10(thunder_layer.loudest_second(signal.sosfiltfilt(sos, pressure)) / thunder_layer.P_REF)
        assert level == pytest.approx(levels['peak_db_by_band'][f'{band:.1f}'], abs=0.5)


def test_the_envelope_holds_the_products_span():
    rng = np.random.default_rng(3)
    for span in (14.0, 34.0):
        env = thunder_layer.envelope(rng, span)
        cumulative = np.cumsum(env**2)
        t = np.arange(len(env)) / thunder_layer.SAMPLE_RATE
        q5, q95 = np.interp(np.array([.05, .95]) * cumulative[-1], cumulative, t)
        assert q95 - q5 == pytest.approx(span, rel=0.01)


def test_the_lunar_cloud_flash_is_deep_and_the_ground_strike_broad():
    moon = thunder_layer.event(PRODUCT, 'moon-cloud-flash')['peak_db_by_band']
    strike = thunder_layer.event(PRODUCT, 'moon-ground-strike')['peak_db_by_band']
    assert max(moon, key=moon.get) in ('16.0', '31.5') and moon['1000.0'] < 0
    assert strike['2000.0'] > 60
