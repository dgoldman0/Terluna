"""The joint synthesis's products bind the files and inputs that made them."""
import hashlib
import json
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PRODUCTS = ('night_sky', 'ledgers', 'magnets_map', 'places', 'upper_air_no', 'air_chemistry')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


@pytest.mark.parametrize('name', PRODUCTS)
def test_product_binds_producer_and_inputs(name):
    d = json.loads((HERE / 'results' / f'{name}.json').read_text())
    assert d['schema'].startswith('terluna.research.joint-synthesis')
    for path, h in d['producer']['files'].items():
        assert digest(ROOT / path) == h, path
    for path, h in d['producer'].get('inputs', {}).items():
        assert digest(ROOT / path) == h, path


def test_loss_rates_carry_cycle_times():
    air = json.loads((HERE / 'results' / 'ledgers.json').read_text())['air']
    for row in air['losses'].values():
        assert row['cycle_years']
