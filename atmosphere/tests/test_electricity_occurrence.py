"""The Moon-wide flash climatology reproduces the box at its own cells and binds its inputs."""
import hashlib
import json
from pathlib import Path

import numpy as np

from atmosphere.electricity import occurrence

PRODUCT = json.loads(occurrence.OUT.read_text())


def test_the_box_cells_take_the_box_rate():
    conv = json.loads(occurrence.CONVECTION.read_text())
    lat, lon = np.asarray(conv['lat_deg']), np.asarray(conv['lon_deg'])
    field = np.asarray(PRODUCT['central_flash_density'])
    at_box = np.mean([occurrence.bilinear(c[0], c[1], field, lat, lon) for c in PRODUCT['box']['gcm_cells']])
    assert abs(at_box / PRODUCT['box']['flash_density'] - 1) < 0.01


def test_sea_factor_lowers_only_the_seas():
    for by_sea in PRODUCT['cases'].values():
        assert abs(by_sea['1']['land_mean'] - by_sea['0.1']['land_mean']) < 1e-9
        assert by_sea['0.1']['sea_mean'] < by_sea['1']['sea_mean']


def test_product_binds_its_inputs():
    digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]
    for name, h in PRODUCT['producer']['inputs'].items():
        assert digest(occurrence.ROOT / name) == h, name
    assert digest(occurrence.HERE / 'occurrence.py') == PRODUCT['producer']['files']['electricity/occurrence.py']
