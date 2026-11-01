"""The Earth error band: declared source patches, edge logic and the committed product."""
import json

import pytest

from atmosphere.snow import earth_error_band as band


def baseline_available():
    return (band.model.BASELINE / 'source/src/mo_depo.f90').exists()


def test_each_band_factor_enters_the_deposition_source_once():
    if not baseline_available():
        pytest.skip('pinned McSnow source is not installed')
    _, sources = band.depo_source('moon')
    body = sources['mo_depo.f90']
    assert body.count('REAL(wp), PUBLIC :: band_factor(3)=1._wp') == 1
    assert body.count('capacitance( sp, Re ) * band_factor(1)') == 1
    assert body.count('IGF_vent=IGF_vent*band_factor(2)') == 1
    speed = 'vt = vterm_snow( atmo, sp ) * band_factor(3)'
    assert body.count(speed) == 1
    start = body.index('SUBROUTINE sp_depo( sp,')
    assert start < body.index(speed) < body.index('END SUBROUTINE sp_depo', start)


def test_drivers_scale_settling_and_report_the_scaled_speed():
    column = band.column_driver()
    assert 'speed=vterm_snow(air,particle)*band_factor(3)' in column
    assert 'vterm_snow(air,particle)*band_factor(3), &' in column
    assert 'band_mass=1._wp; band_axes=1._wp; band_speed=1._wp' in column
    assert 'vterm_snow(air,sp)*band_factor(3)' in band.weather_driver()


def test_edges_change_one_factor_at_a_time():
    sizing = {f: dict(low=dict(factor=.5 + i), high=dict(factor=2. + i)) for i, f in enumerate(band.FACTORS)}
    edges = band.edges(sizing)
    assert edges['mass_low'] == (.5, 1., 1.) and edges['speed_high'] == (1., 1., 4.)
    assert len(edges) == 6


def test_combined_edges_follow_the_size_each_edge_produced():
    sizing = {f: dict(low=dict(factor=.8), high=dict(factor=1.2)) for f in band.FACTORS}
    sizes = {'mass_low': 1., 'mass_high': 3., 'axes_low': 5., 'axes_high': 2., 'speed_low': 4., 'speed_high': 4.}
    small, large = band.combined(band.edges(sizing), sizing, sizes.get)
    assert small == (.8, 1.2, .8) and large == (1.2, .8, 1.2)


@pytest.fixture(scope='module')
def product():
    if not band.PRODUCT.exists():
        pytest.skip('Earth error band product not assembled')
    return json.loads(band.PRODUCT.read_text())


def test_committed_band_reproduces_the_frozen_component_at_factor_one(product):
    assert all(row['identical'] for row in product['parity']['column'])
    assert all(product['parity']['weather']['identical'].values())
    for factor in band.FACTORS:
        assert product['band'][factor]['low']['factor'] < 1 < product['band'][factor]['high']['factor']


def test_committed_spreads_summarize_their_tables(product):
    for case, table in product['column'].items():
        sizes = [v['diameter_m'] for k, v in table.items() if not k.endswith('dt0.5')]
        assert product['spread']['band'][case]['ratio'] == pytest.approx(max(sizes) / min(sizes), rel=1e-12)
    baselines = [t['baseline']['diameter_m'] for t in product['column'].values()]
    assert product['spread']['moisture']['column']['ratio'] == pytest.approx(max(baselines) / min(baselines), rel=1e-12)
