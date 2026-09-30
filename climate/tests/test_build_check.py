"""Checks of the build comparison's pure pieces: the one-bit nudge and the verdict."""
import numpy as np
import pytest

from climate.crm import build_check as bc


def test_a_nudge_flips_the_last_bit_of_every_field_value_and_spares_zeros_and_whole_numbers(tmp_path):
    values = np.array([0.0, 1.0, 2.0, 1.5e-3, 91.97725, -3.25e-7, 294.23636, 7.0], dtype='<f4')
    for which, flipped, count in ((0, [0, 0, 0, 1, 1, 1, 1, 0], 4), (1, [0, 0, 0, 0, 1, 0, 1, 0], 2)):
        path = tmp_path / f'cm1rst_t00000{which}_s.dat'
        values.tofile(path)
        done = bc.nudge(path, which)
        after = np.fromfile(path, dtype='<f4')
        assert list(after.view('<u4') ^ values.view('<u4')) == flipped
        assert done == dict(file=path.name, values=count, of=8, every_other=bool(which))


def test_a_build_passes_when_identical_or_within_the_largest_nudge():
    def run(**rms):
        return dict(identical=False, differences=dict(fields={k: dict(rms=[v, v]) for k, v in rms.items()}))
    nudged = [run(t2=0.4, qv=1e-4), run(t2=0.6, qv=2e-4)]
    assert bc.verdict(dict(identical=True), nudged) == dict(identical=True, within_nudge=True)
    assert bc.verdict(run(t2=0.5, qv=2e-4), nudged)['within_nudge']
    fails = bc.verdict(run(t2=0.7, qv=1e-4), nudged)
    assert not fails['within_nudge'] and fails['fields']['t2']['ratio'] == pytest.approx(0.7 / 0.6)
