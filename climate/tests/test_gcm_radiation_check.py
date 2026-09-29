"""Checks of the GCM radiation check's pure pieces: which columns it takes and how it builds them."""
import numpy as np

from climate.gcm import radiation_check as g


def fields():
    lat = np.array([55.0, 45.0, 25.0, 2.8, -2.8, -25.0, -45.0, -65.0])
    lon = np.arange(8) * 45.0
    lsm = np.zeros((8, 8), bool)
    lsm[:, :4] = True                                                     # land in the west, sea in the east
    prw = np.add.outer(np.linspace(150, 200, 8), np.arange(8.0))
    height = np.where(lsm, 1000.0, 0.0)
    height[3, 2] = 8000.0                                                 # a tropical highland
    return dict(lat=lat, lon=lon, lsm=lsm, prw=prw, height_m=height)


def test_columns_span_dry_to_humid_land_a_highland_the_site_and_seas_by_latitude():
    f = fields()
    picks = g.choose_columns(f, site=(0.0, 0.0))
    labels = [p['label'] for p in picks]
    assert labels[0] == 'site land' and 'tropical highland' in labels
    high = next(p for p in picks if p['label'] == 'tropical highland')
    assert (high['j'], high['i']) == (3, 2)
    for p in picks:
        assert f['lsm'][p['j'], p['i']] == p['label'].startswith(('site', 'land', 'tropical'))
    seas = [p for p in picks if p['label'].startswith('sea')]
    assert [abs(f['lat'][p['j']]) for p in seas] == sorted(abs(f['lat'][p['j']]) for p in seas)
    assert len({(p['j'], p['i']) for p in picks}) == len(picks)          # no column twice


def test_a_column_profile_runs_from_the_ground_to_50_pa():
    sigma = np.array([0.1, 0.5, 0.9])
    f = dict(sigma=sigma, ps=np.full((1, 1), 1000.0), ta=np.array([220.0, 260.0, 290.0])[:, None, None],
             hus=np.array([1e-5, 1e-3, 1e-2])[:, None, None], tas=np.full((1, 1), 295.0))
    prof = g.column_profile(f, 0, 0)
    assert prof['ps'] == 1.0e5 and prof['t2'] == 295.0
    assert np.allclose(prof['p'][:3], sigma * 1.0e5) and np.isclose(prof['p'][-1], 50.0)
    assert np.all(prof['t'][3:] == 220.0) and np.all(prof['qv'][3:] == 1e-7)
