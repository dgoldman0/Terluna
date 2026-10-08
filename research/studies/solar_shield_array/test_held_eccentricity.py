import json

import pytest

from . import held_eccentricity as he
from . import ring_layout as rl


@pytest.mark.skipif(not he.OUT.exists(), reason='held_eccentricity.json not written')
def test_product_is_current_and_reproduces_the_layout():
    product = json.loads(he.OUT.read_text())
    assert product['schema'] == he.SCHEMA
    for name, value in product['producer']['files'].items():
        assert he.digest(he.ROOT/name) == value, name
    assert product['inputs']['results/ring_layout.json'] == he.digest(rl.OUT)
    assert product['checks']['reproduces_layout']
    for row in product['by_clearance'].values():
        share = row['share']
        assert 0 <= share['mean'] <= share['max'] < 1
        # Holding takes more at the wider clearance, where the radius steps leave less room for the eccentricity to
        # rise from ring to ring.
    assert product['by_clearance']['1000']['share']['max'] > product['by_clearance']['600']['share']['max']
