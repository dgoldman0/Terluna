"""The page builder's lists of regions, districts and coastal seas: the links each unit's card names."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import render as R                                     # noqa: E402


def unit(code, **more):
    return dict(code=code, name=code.lower(), tongue='t', at=dict(turned=[0.0, 0.0]), area_km2=1.0, description='',
                **more)


def test_links_hold_every_unit_a_card_names_once_and_in_order():
    # a coastal sea that ten districts feed, two of which also face it, and a land district with nine coastal seas
    lands = [unit(f'L1.{k}', region='L1', coasts=['S1.1']) for k in range(1, 11)]
    lands[0]['coasts'] = [f'S1.{k}' for k in range(1, 10)]
    seas = [unit(f'S1.{k}', region='L1') for k in range(1, 10)]
    seas[0].update(faces=['L1.2', 'L1.1'], fed_by=[f'L1.{k}' for k in range(1, 11)] + ['L9.9'])
    record = dict(regions=[unit('L1', kind='land')], districts=lands, coastal_seas=seas,
                  tongues=dict(t=dict(label='Tongue')), summary={})
    data = R.districts_data(record)
    by_code = {u['code']: u for u in data['districts']}
    index = {u['code']: u['index'] for u in data['districts']}
    assert [u['index'] for u in data['districts']] == list(range(1, 20))          # as the districts texture holds them
    # the shore districts first, then the others that feed it, none twice, and none the record does not hold
    assert by_code['S1.1']['links'] == [index[c] for c in ['L1.2', 'L1.1'] + [f'L1.{k}' for k in range(3, 11)]]
    assert by_code['L1.1']['links'] == [index[f'S1.{k}'] for k in range(1, 10)]   # all nine: no cut to seven
    assert by_code['L1.5']['links'] == [index['S1.1']]
    assert data['regions'][0]['links'] == [] and data['regions'][0]['index'] == 0
