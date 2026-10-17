"""Checks of the copy of PlaSim's clear-sky shortwave: its levels, its energy balance and its limits."""
import numpy as np

from climate.gcm import plasim_shortwave as p

SIGMA = np.array([0.0383, 0.1191, 0.21085, 0.31685, 0.4368, 0.5668, 0.69935, 0.82335, 0.9241, 0.9833])  # T21 L10


def moist_column(**over):
    q = 0.018 * SIGMA ** 3                                                # moist near the ground, dry aloft
    return dict(dict(sigma=SIGMA, ps=121590.0, t=200.0 + 95.0 * SIGMA, q=q, albedo=0.2, rcoeff=0.3857,
                     solar1=0.502, gravity=1.6242), **over)


def test_half_levels_are_plasims():
    h = p.half_levels(SIGMA)
    assert np.allclose(h[[0, 4, 8, 9]], [0.0766, 0.5, 0.9666, 1.0], atol=1e-4)


def test_energy_balances_and_a_clear_column_reflects_its_ground():
    for mu in (0.1, 0.5, 1.0):
        s = p.clear_sky(moist_column(), mu)
        assert np.isclose(sum(s.values()), 1.0)
        assert 0.0 < s['air'] < 0.5 and 0.0 < s['reflected'] < 0.5
    bare = p.clear_sky(moist_column(q=np.zeros(10), rcoeff=0.0), 0.7)
    assert np.isclose(bare['reflected'], 0.2) and np.isclose(bare['air'], 0.0, atol=1e-12)
    rayleigh = p.clear_sky(moist_column(q=np.zeros(10)), 0.7)
    assert np.isclose(rayleigh['air'], 0.0, atol=1e-12) and rayleigh['reflected'] > 0.2   # scattering only


def test_every_form_raises_the_air_absorption_from_plasims():
    col = moist_column()
    base = p.daily_split(col, 0.0)['air']
    for form, k in (('path', 5.0), ('scale', 1.3), ('linear', 0.002)):
        raised = p.daily_split(col, 0.0, form, k)['air']
        assert raised > base + 0.02
    assert np.isclose(p.daily_split(col, 0.0, 'path', 1.0)['air'], base)  # the neutral settings are PlaSim's own
    assert np.isclose(p.daily_split(col, 0.0, 'scale', 1.0)['air'], base)
    # A grazing Sun through a very wet column: the scaled band keeps a sliver of transmission, as the patch does.
    assert p.band_absorptivity(np.array([1e4]), 0.498, 'scale', 1.356)[0] == 0.999
    wet = p.clear_sky(moist_column(q=np.full(10, 0.05)), 0.02, 'scale', 1.356)
    assert np.isclose(sum(wet.values()), 1.0) and min(wet.values()) >= 0.0


FIT = [3.0, 2000.0, 0.7, 7.0, 0.8]                                       # a CO2-like fit (a, b, c, d, n)


def test_a_gas_band_absorbs_nothing_at_no_path_and_never_all():
    a = p.gas_absorptivity(np.array([0.0, 1e-4, 1e-2, 1.0, 1e3, 1e9]), FIT)
    assert a[0] == 0.0 and np.all(np.diff(a) > 0.0) and a[-1] <= 0.999


def test_layer_transmissivities_telescope_to_the_whole_path():
    path = np.linspace(0.1, 1.0, 10) * 1e-3
    a = lambda u: p.gas_absorptivity(u, FIT)
    down, up = p.path_transmissivities(path, 2.0, a)
    assert np.isclose(np.prod(down), 1.0 - a(2.0 * path.sum()))
    # The light the ground returns crosses the whole slant path, then BETA times the vertical path back up.
    assert np.isclose(np.prod(down) * np.prod(up), 1.0 - a((2.0 + p.BETA) * path.sum()))


def test_gas_bands_are_plasims_own_at_no_gas_and_take_sunlight_from_the_ground():
    gases = dict(fits={name: FIT for name, _, _ in p.GAS_BANDS}, fractions=dict(CO2=0.0, O2=0.0),
                 magnification=[0.49, 0.42, 0.069])
    lunar = [0.184, 0.170, 0.0136]
    for mu in (0.05, 0.5, 1.0):
        base = p.clear_sky(moist_column(magnification=lunar), mu)
        none = p.clear_sky(moist_column(magnification=lunar, gases=gases), mu)
        assert all(np.isclose(base[k], none[k], rtol=0.0, atol=1e-14) for k in base)
        moon = p.clear_sky(moist_column(magnification=lunar, gases=dict(gases, fractions=dict(CO2=6.1e-4, O2=0.194))), mu)
        assert np.isclose(sum(moon.values()), 1.0)
        assert moon['air'] > base['air'] and moon['ground'] < base['ground'] and moon['reflected'] < base['reflected']


def test_the_water_exponent_at_one_is_plasims_and_lower_weights_the_upper_air_more():
    col = moist_column()
    top = lambda net: net[0] - net[3]                                   # the top three layers' absorption
    for mu in (0.2, 1.0):
        base = p.net_flux(col, mu, 'scale', 1.3)
        assert np.array_equal(base, p.net_flux(dict(col, water_exponent=1.0), mu, 'scale', 1.3))
        assert top(p.net_flux(dict(col, water_exponent=0.5), mu, 'scale', 1.3)) > top(base)
