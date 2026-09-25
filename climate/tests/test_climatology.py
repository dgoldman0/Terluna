"""The GCM climatology's cloud decks and its composite of cloud against the Sun."""
import numpy as np

from climate.gcm import climatology as cg


def test_decks_split_layers_at_sigma():
    sigma = np.array([0.1, 0.3, 0.6, 0.9])
    cl = np.zeros((2, 4, 1, 1))
    cl[:, 0], cl[:, 1], cl[:, 2], cl[:, 3] = 0.2, 0.5, 0.1, 0.4
    low, high = cg.decks(cl, sigma)
    assert np.allclose(low, 0.4) and np.allclose(high, 0.5)


def test_sun_composite_recovers_an_afternoon_peak():
    lat = np.array([10.0, -10.0])
    lon = np.arange(64) * 5.625
    steps = 40
    subsolar = (np.arange(steps) * 36.0) % 360.0          # the Sun moves a tenth of a day per output
    czen = np.cos(np.radians((lon[None, :] - subsolar[:, None] + 180.0) % 360.0 - 180.0)).clip(0.0)
    czen = np.repeat(czen[:, None, :], lat.size, axis=1)
    hour = (lon[None, :] - subsolar[:, None] + 180.0) % 360.0 - 180.0
    clt = np.repeat(np.exp(-((hour - 50.0) / 25.0) ** 2)[:, None, :], lat.size, axis=1)
    assert np.allclose(cg.subsolar_longitudes(czen, lat, lon), subsolar % 360.0, atol=5.7)
    composite, centres = cg.sun_composite(clt, czen, lat, lon)
    assert composite.shape == (lat.size, cg.HOUR_BINS)
    assert abs(centres[np.argmax(composite[0])] - 50.0) <= 10.0
    assert composite[0].min() < 0.05
