"""The ecology screen binds its inputs, reproduces the surface-light product's UV index and keeps its physics checkable."""
import hashlib
import json
from pathlib import Path

import numpy as np

from biosphere.ecology import screen
from shared.provenance import constants_changed

ROOT = Path(__file__).resolve().parents[2]
PRODUCT = json.loads(screen.OUT.read_text())


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def test_product_binds_producer_inputs_and_constants():
    assert PRODUCT['schema'] == screen.SCHEMA
    for path, h in {**PRODUCT['producer']['files'], **PRODUCT['producer']['inputs']}.items():
        assert digest(ROOT / path) == h, path
    assert not constants_changed(PRODUCT['producer']['constants'])


def test_uv_index_matches_the_surface_light_product():
    surface = json.loads(screen.SURFACE.read_text())
    for case, key in (('moon_1.2atm/design', 'moon'), ('earth_control/unfiltered', 'earth')):
        stored = surface['cases'][case]['summaries']['90.0']['0.1']['uv_index']
        assert abs(PRODUCT['ultraviolet']['by_sun']['90.0']['uv_index'][key] - stored) / stored < 0.01


def test_no_uv_b_reaches_the_ground():
    for row in PRODUCT['ultraviolet']['by_sun'].values():
        assert row['uv_b_photons_umol_m2_s']['moon'] == 0.0
        assert row['previtamin_d_band_photons_ratio'] < 1e-4


def test_pigment_template_peaks_at_its_maximum():
    w = np.arange(300.0, 700.0, 0.5)
    for peak in screen.RECEPTORS_NM.values():
        s = screen.govardovskii(w, peak)
        assert abs(w[np.argmax(s)] - peak) <= 3.0
        assert 0.95 < s.max() < 1.05


def test_erythema_spectrum_follows_the_cie_definition():
    assert screen.erythema(np.array([290.0]))[0] == 1.0
    assert np.isclose(screen.erythema(np.array([318.0]))[0], 10 ** (0.094 * -20))
    assert np.isclose(screen.erythema(np.array([350.0]))[0], 10 ** (0.015 * -210))


def test_stokes_fall_on_earth_matches_the_textbook():
    # A 10 um sphere of unit density falls at about 3.0-3.1 mm/s in sea-level air (slip included).
    v, _ = screen.fall_speed(10e-6, 1000.0, screen.STANDARD_GRAVITY, *screen.EARTH_AIR, screen.AIR_MOLAR_MASS)
    assert 2.9e-3 < v < 3.2e-3


def test_lunar_fall_scales_with_gravity():
    for row in PRODUCT['settling']['particles'].values():
        ratio = row['fall_m_s']['moon'] / row['fall_m_s']['earth']
        assert 0.14 < ratio < 0.18          # g ratio 0.166, with the warmer air's viscosity and slip


def test_trace_gas_inventories_are_consistent():
    gases = PRODUCT['trace_gases']
    assert np.isclose(gases['methane']['tg_per_ppm'], gases['air_mol'] * 1e-6 * 16.04 / 1e12, rtol=1e-3)
    lo, hi = gases['methane']['lifetime_years'][0], gases['methane']['lifetime_years'][-1]
    assert 100 < lo < hi < 5000
    assert all(r > 0 for r in gases['nitrous_oxide']['rise_ppb_per_year'])
