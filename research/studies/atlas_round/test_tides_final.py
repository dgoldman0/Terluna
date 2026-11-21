"""Checks of the design surface's tide: with today's frame the turned responses are tides.py's own, a turn about the
pole carries the forcing's monomials with the ground, and once computed the record holds together (the seven seas with
finite ranges, the nodes the moisture map's sea pixels, the frame the climate's, the product matching its hash)."""
import json

import numpy as np
import pytest

from geography import tides as T
from research.studies.atlas_round import final_drainage as FD
from research.studies.atlas_round import tides_final as TF


def small_seas():
    lat = 90.0 - 5.0 * (np.arange(36) + 0.5)
    lon = 5.0 * (np.arange(72) + 0.5)
    labels = np.zeros((36, 72), int)
    labels[10:16, 5:20] = 1
    labels[20:26, 40:60] = 2
    return dict(labels=labels, depth=np.where(labels > 0, 100.0, 0.0), lat=lat, lon=lon)


def test_the_turned_responses_are_tides_own_in_todays_frame():
    response = T.Response(small_seas(), (1, 2), nmax=8, water_density=1000.0)
    ours = TF.basis_responses(response, np.eye(3))
    theirs = T.basis_responses(response)
    for a, b in zip(ours[:3], theirs[:3]):
        np.testing.assert_allclose(a, b, rtol=0, atol=1e-12)


def test_a_turn_about_the_pole_carries_the_forcing_with_the_ground():
    response = T.Response(small_seas(), (1, 2), nmax=8, water_density=1000.0)
    quarter = np.array([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])   # the turned x axis is today's y
    _, _, values, _ = TF.basis_responses(response, quarter)
    la, lo = np.meshgrid(np.radians(response.grid['lat']), np.radians(response.grid['lon']), indexing='ij')
    today = T.monomial_values(T.unit_points(la, lo))
    np.testing.assert_allclose(values[..., T.MONOMIALS.index('x')], today[..., T.MONOMIALS.index('y')], atol=1e-12)
    np.testing.assert_allclose(values[..., T.MONOMIALS.index('z')], today[..., T.MONOMIALS.index('z')], atol=1e-12)


@pytest.fixture(scope='module')
def record():
    if not TF.OUT.exists():
        pytest.skip('run python -m research.studies.atlas_round.tides_final first')
    return json.loads(TF.OUT.read_text())


def test_the_record_holds_together(record):
    assert sorted(int(k) for k in record['seas']) == list(range(1, 8))
    for sea in record['seas'].values():
        m = sea['monthly_range_median_m']
        assert 0.0 <= m['area_median'] <= m['area_p95'] <= m['maximum'] < 50.0
        assert sea['nodes'] > 0
    frame_path, _ = FD.frame_record()
    assert record['settings']['frame']['record'] == str(frame_path.relative_to(TF.ROOT))
    assert record['settings']['frame']['sha256'] == TF.sha256(frame_path)
    assert len(record['dynamics']) == len(TF.DYNAMIC_SEAS)
    if TF.PRODUCT.exists():
        assert TF.sha256(TF.PRODUCT) == record['product']['sha256']
        with np.load(TF.PRODUCT, allow_pickle=False) as z:
            assert z['rows'].size == sum(s['nodes'] for s in record['seas'].values())
            assert np.isfinite(z['monthly_range_median_m']).all()
