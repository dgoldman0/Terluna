"""The GCM's rain evaporation: the fit to CM1's share and the product's two measures."""
import json

import pytest

from climate.gcm import rain_evaporation as re_


def test_the_fit_interpolates_gamma_in_its_log():
    runs = [dict(gamma=g, shares=dict(cm1_box=dict(evaporated_share=s))) for g, s in ((0.01, 0.3), (0.04, 0.6))]
    assert re_.fit(runs, 0.45) == pytest.approx(0.02, rel=1e-9)
    assert re_.fit(runs, 0.9) is None                                    # outside the runs given


def test_the_product_matches_cm1_between_its_runs_and_its_two_measures_agree():
    if not re_.OUT.exists():
        pytest.skip('the rain evaporation has not been measured')
    out = json.loads(re_.OUT.read_text())
    assert out['cm1']['lunar_limits'] == pytest.approx(0.48)
    gammas = sorted(m['gamma'] for m in out['measured'])
    assert gammas[0] < out['gamma_for_cm1_lunar_limits'] < gammas[-1]
    for m in out['measured']:
        g = m['shares']['global']
        assert abs(g['evaporated_share'] - g['evaporated_share_by_residual']) < 0.02
