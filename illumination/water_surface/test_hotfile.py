"""SWAN restart files read as variance spectra, checked against a built file and the nearside run's printed spectra."""
import struct
from pathlib import Path

import numpy as np
import pytest

from .hotfile import read_hotfile
from .model import iter_swan

ROOT = Path(__file__).resolve().parents[2]


def fortran(*records):
    return b''.join(struct.pack('<i', len(r)) + r + struct.pack('<i', len(r)) for r in records)


def test_hotfile_spectra_are_read_as_variance_density(tmp_path):
    frequency, direction = np.array([.1, .2]), np.array([45., 135., 225., 315.])
    lonlat = [(10.0, 20.0), (10.0, 21.0)]
    action = np.arange(16, dtype='<f4').reshape(2, 2, 4) / 7
    records = [b'41.51AB'.ljust(20), b'Test sea'.ljust(20), *(struct.pack('<i', 1),) * 3,
               struct.pack('<3i', 2, 1, 2), *(struct.pack('<2f', *xy) for xy in lonlat),
               struct.pack('<i', 2), *(struct.pack('<f', f) for f in frequency),
               struct.pack('<i', 4), *(struct.pack('<f', d) for d in direction), b'20000101.030000'.ljust(20)]
    for p in range(2):
        records += [struct.pack('<i', p + 2), action[p].tobytes()]
    path = tmp_path / 'end.hot'
    path.write_bytes(fortran(*records))
    spectra = read_hotfile(path, nodes=[(10.0, 21.2)])
    assert spectra['index'].tolist() == [1] and spectra['time'].hour == 3
    expected = action[1] * (2 * np.pi * frequency)[:, None] * 2 * np.pi * np.pi / 180
    assert spectra['variance'][0] == pytest.approx(expected.astype(float), rel=1e-6)
    path.write_bytes(fortran(*records[:-1]))
    with pytest.raises(ValueError):
        read_hotfile(path)


def test_a_nearside_restart_file_matches_the_printed_spectra_at_the_same_time():
    segment = ROOT / 'research/runs/waves/nearside/dt300/s0720_0768'
    if not (segment / 'end.hot').exists():
        pytest.skip('The nearside wave run lives on the research drive')
    hot = read_hotfile(segment / 'end.hot', nodes=np.loadtxt(segment / 'refs.xy'))
    printed = [r for r in iter_swan(segment / 'refs2d.spc') if r['time'] == hot['time']]
    assert len(printed) == 1
    printed = printed[0]
    assert hot['lonlat'] == pytest.approx(printed['locations'])
    assert np.allclose(hot['direction'], printed['direction'], atol=1e-3)
    scale = printed['variance'].max(axis=(1, 2))[:, None, None]
    assert np.max(np.abs(hot['variance'] - printed['variance']) / scale) < 1e-4
