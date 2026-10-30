"""Junction loads in a falling snowflake: the contact network's checks and its product."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from atmosphere.snow import breakup_network as bn


def test_contacts_equilibrium_and_tree_sums_hold_on_a_small_aggregate():
    check = bn.selftest(np.random.default_rng(bn.SEED + 2))
    assert check['passed']
    assert check['edges'] == 127                          # one contact per merge: a spanning tree of 128 plates
    assert check['contact_distance_error'] < 1e-9
    assert check['closest_pair'] > 1.0 - 1e-9
    assert check['tree_sum_vs_brute_force'] < 1e-8
    assert check['uniform_drag_max_moment'] < 1e-9        # drag in proportion to weight loads no junction


def test_product_was_written_by_this_module_and_its_ratios_follow_its_fits():
    product = json.loads(bn.PRODUCT.read_text())
    assert product['schema'] == bn.SCHEMA
    assert product['producer_sha256'] == hashlib.sha256(Path(bn.__file__).read_bytes()).hexdigest()
    for name, digest in product['inputs']['imported_files_sha256'].items():
        assert hashlib.sha256((bn.HERE / name).read_bytes()).hexdigest() == digest, name
    assert product['selftest']['passed']
    assert product['max_equilibrium_residual'] < 1e-8
    g_ratio = product['inputs']['gravity']['ratio']
    for key, rows in product['spans'].items():
        fit = product['fits'][key]
        for row in rows.values():
            b = row['mass_exponent']
            assert row['load_analogy_ratio'] == pytest.approx(g_ratio ** (1 / b))
            assert row['force_ratio'] == pytest.approx(g_ratio ** (1 / (b + fit['force_norm_exponent'])))
            assert row['moment_ratio'] == pytest.approx(g_ratio ** (1 / (b + 1 + fit['moment_norm_exponent'])))
