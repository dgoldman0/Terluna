"""Checks of the CM1 file store on synthetic cases: compressing, reading back bit for bit, taking up an interrupted
compression, the free-space floor and the listings."""
import hashlib
import json
import os

import numpy as np
import pytest

from climate.crm import cm1_files as cf


class Crash(BaseException):
    """A stop part-way through a step, as a killed process or a power cut leaves it."""


@pytest.fixture(autouse=True)
def quiet_machine(monkeypatch):
    """No CM1 process works in the synthetic cases, and the pieces are small enough to need many to a file."""
    monkeypatch.setattr(cf, 'cm1_processes', lambda: [])
    monkeypatch.setattr(cf, 'CHUNK', 1 << 16)


def synthetic_case(folder, outputs=(1, 2, 3), restarts=(1, 2)):
    """A case folder with 3-D outputs (float32 fields with structure and noise, about 0.4 MB each), complete restart
    sets of six files of several sizes, and files that stay raw."""
    folder.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(7)
    data = {}
    for n in outputs:
        field = (np.linspace(250.0, 320.0, 100_000) + rng.normal(0.0, 0.5, 100_000)).astype('<f4')
        path = folder / f'cm1out_t{n:06d}_s.dat'
        field.tofile(path)
        data[path.name] = path.read_bytes()
    for n in restarts:
        for k, part in enumerate(cf.RESTART_PARTS):
            values = rng.integers(0, 2 ** 16, 1000 * (k + 1), dtype='<u4')
            path = folder / f'cm1rst_t{n:06d}_{part}.dat'
            values.tofile(path)
            data[path.name] = path.read_bytes()
    for name, text in (('cm1out_stats.dat', b'\x00' * 400), ('cm1out_diag_000001_s.dat', b'\x01' * 400),
                       ('cm1out_t000001_w.dat', b'\x02' * 400), ('namelist.input', b' &param1\n /\n')):
        (folder / name).write_bytes(text)
    return data


def test_the_round_trip_is_bit_for_bit(tmp_path):
    case = tmp_path / 'case'
    data = synthetic_case(case)
    before = sorted(p.name for p in case.iterdir())
    result = cf.compress(case, floor_bytes=0)
    assert sorted(result['compressed']) == sorted(data) and not result['failed'] and not result['left']
    for name in data:
        assert not (case / name).exists() and cf.compressed(case / name).exists()
    for name in ('cm1out_stats.dat', 'cm1out_diag_000001_s.dat', 'cm1out_t000001_w.dat', 'namelist.input'):
        assert (case / name).exists() and not cf.compressed(case / name).exists()      # these stay raw
    lines = cf.manifest(case)
    assert sorted(lines) == sorted(data)
    for name, raw in data.items():
        assert lines[name]['raw_sha256'] == hashlib.sha256(raw).hexdigest() and lines[name]['raw_bytes'] == len(raw)
        assert lines[name]['zst_bytes'] == cf.compressed(case / name).stat().st_size
    assert result['raw_bytes'] == sum(len(v) for v in data.values()) == result['raw_bytes_after']
    assert result['disk_bytes_after'] < result['disk_bytes_before']
    assert cf.verify(case)['ok'] == sorted(data)
    assert sorted(cf.decompress(case, floor_bytes=0)) == sorted(data)
    for name, raw in data.items():                                     # back raw, the .zst kept beside it
        assert (case / name).read_bytes() == raw and cf.compressed(case / name).exists()
    assert sorted(p.name for p in case.iterdir() if not p.name.endswith(cf.SUFFIX)) == sorted(before + [cf.MANIFEST])


def test_a_compressed_file_keeps_the_raw_file_s_time(tmp_path):
    case = tmp_path / 'case'
    synthetic_case(case, outputs=(1,), restarts=())
    path = case / 'cm1out_t000001_s.dat'
    os.utime(path, ns=(1_700_000_000_000_000_000, 1_700_000_000_000_000_000))
    cf.compress(case, floor_bytes=0)
    assert cf.compressed(path).stat().st_mtime_ns == 1_700_000_000_000_000_000
    cf.decompress(case, [path.name], floor_bytes=0)
    assert path.stat().st_mtime_ns == 1_700_000_000_000_000_000


def test_reading_a_compressed_file_equals_reading_the_raw_one(tmp_path):
    raw_case, packed = tmp_path / 'raw', tmp_path / 'packed'
    synthetic_case(raw_case)
    synthetic_case(packed)
    cf.compress(packed, floor_bytes=0)
    for name in ('cm1out_t000002_s.dat', 'cm1rst_t000001_s.dat', 'cm1rst_t000002_x.dat'):
        a, b = raw_case / name, packed / name
        assert not b.exists() and cf.exists(b)
        for dtype, kw in (('<f4', {}), ('<u4', {}), ('u1', {}), ('<f4', dict(count=7)), ('<f4', dict(offset=8)),
                          ('<f8', dict(count=3, offset=4)), ('<f4', dict(count=10 ** 9))):
            x, y = np.fromfile(a, dtype=dtype, **kw), cf.fromfile(b, dtype=dtype, **kw)
            assert x.dtype == y.dtype and x.shape == y.shape and x.tobytes() == y.tobytes()
            assert y.flags.writeable                                    # as np.fromfile's
            assert cf.fromfile(a, dtype=dtype, **kw).tobytes() == x.tobytes()   # a raw file is np.fromfile's
        m, n = np.memmap(a, dtype='<f4', mode='r'), cf.memmap(b, dtype='<f4', mode='r')
        assert m.shape == n.shape and m.tobytes() == n.tobytes() and not n.flags.writeable
        m, n = (np.memmap(a, dtype='<f4', mode='r', offset=4, shape=(5,)),
                cf.memmap(b, dtype='<f4', mode='r', offset=4, shape=(5,)))
        assert m.tobytes() == n.tobytes()
        assert isinstance(cf.memmap(a, dtype='<f4', mode='r'), np.memmap)
        assert cf.sha256(b) == cf.sha256(a) == hashlib.sha256(a.read_bytes()).hexdigest()
        assert cf.size(b) == cf.size(a) == a.stat().st_size
    with pytest.raises(ValueError, match='only reads'):
        cf.memmap(packed / 'cm1out_t000002_s.dat', dtype='<f4', mode='r+')
    with pytest.raises(FileNotFoundError):
        cf.fromfile(packed / 'cm1out_t000009_s.dat', dtype='<f4')
    with pytest.raises(ValueError):
        cf.fromfile(packed / 'cm1out_t000002_s.dat', dtype='<f4', offset=10 ** 9)
    assert cf.written_bytes(packed) == cf.written_bytes(raw_case) == sum(
        p.stat().st_size for p in raw_case.iterdir() if p.is_file() and not p.is_symlink())


def test_the_manifest_line_stands_before_the_raw_file_goes(tmp_path, monkeypatch):
    case = tmp_path / 'case'
    data = synthetic_case(case, outputs=(1, 2), restarts=(1,))
    seen = []
    remove = cf._remove_raw

    def watched(path):
        lines = [json.loads(line) for line in (case / cf.MANIFEST).read_text().splitlines()]
        assert path.name in {line['name'] for line in lines}            # the line is written ...
        assert cf.compressed(path).exists() and path.exists()           # ... beside its .zst, the raw file still there
        assert not path.with_name(path.name + cf.SUFFIX + cf.PART).exists()
        seen.append(path.name)
        remove(path)

    monkeypatch.setattr(cf, '_remove_raw', watched)
    cf.compress(case, floor_bytes=0)
    assert sorted(seen) == sorted(data)


def test_a_corrupted_compressed_file_is_caught_and_the_raw_file_kept(tmp_path, monkeypatch):
    case = tmp_path / 'case'
    data = synthetic_case(case, outputs=(1, 2), restarts=())
    differs = cf._differs

    def damaged(part, digest, count):                                   # a byte flipped on the way to the drive
        if part.name.startswith('cm1out_t000001'):
            with part.open('r+b') as f:
                f.seek(part.stat().st_size // 2)
                byte = f.read(1)
                f.seek(-1, os.SEEK_CUR)
                f.write(bytes([byte[0] ^ 0xFF]))
        return differs(part, digest, count)

    monkeypatch.setattr(cf, '_differs', damaged)
    result = cf.compress(case, floor_bytes=0)
    assert [f['name'] for f in result['failed']] == ['cm1out_t000001_s.dat']
    assert result['compressed'] == ['cm1out_t000002_s.dat']                        # and the work went on
    failed = case / 'cm1out_t000001_s.dat'
    assert failed.read_bytes() == data[failed.name]
    assert not cf.compressed(failed).exists() and not list(case.glob('*.part'))
    assert failed.name not in cf.manifest(case)
    monkeypatch.setattr(cf, '_differs', differs)                       # a later run takes it up
    assert cf.compress(case, floor_bytes=0)['compressed'] == [failed.name]
    good = cf.compressed(case / 'cm1out_t000002_s.dat')                # damage after the commit: verify and the
    raw = bytearray(good.read_bytes())                                 # readers catch it
    raw[len(raw) // 2] ^= 0xFF
    good.write_bytes(bytes(raw))
    report = cf.verify(case)
    assert [m['name'] for m in report['mismatched']] == ['cm1out_t000002_s.dat'] and report['ok'] == [failed.name]
    with pytest.raises(RuntimeError, match='could not decompress'):
        cf.fromfile(case / 'cm1out_t000002_s.dat', dtype='<f4')
    with pytest.raises(RuntimeError):
        cf.decompress(case, ['cm1out_t000002_s.dat'], floor_bytes=0)
    assert not (case / 'cm1out_t000002_s.dat').exists() and not list(case.glob('*.part'))


@pytest.mark.parametrize('stop', ['compressing', 'checking', 'renaming', 'writing_line', 'cutting_line', 'removing'])
def test_each_interrupted_state_recovers(tmp_path, monkeypatch, stop):
    """A stop at each step of one file leaves a state the next compress settles, ending as a clean compression."""
    case = tmp_path / 'case'
    data = synthetic_case(case, outputs=(1, 2), restarts=())
    first = case / 'cm1out_t000001_s.dat'
    real = dict(differs=cf._differs, rename=os.rename, append=cf._append_line, remove=cf._remove_raw)

    def crash(*args, **kwargs):
        raise Crash(stop)

    if stop == 'compressing':                                           # mid-stream: a short .part stands
        def short(part, digest, count):
            part.write_bytes(part.read_bytes()[:100])
            raise Crash(stop)
        monkeypatch.setattr(cf, '_differs', short)
    elif stop == 'checking':
        monkeypatch.setattr(cf, '_differs', crash)
    elif stop == 'renaming':
        monkeypatch.setattr(cf.os, 'rename', crash)
    elif stop == 'writing_line':                                        # renamed, no line
        monkeypatch.setattr(cf, '_append_line', crash)
    elif stop == 'cutting_line':                                        # half a line written
        def half(folder, entry):
            with (folder / cf.MANIFEST).open('ab') as out:
                out.write(json.dumps(entry).encode()[:40])
            raise Crash(stop)
        monkeypatch.setattr(cf, '_append_line', half)
    else:                                                               # line written, raw file still there
        monkeypatch.setattr(cf, '_remove_raw', crash)
    with pytest.raises(Crash):
        cf.compress(case, floor_bytes=0)
    assert first.read_bytes() == data[first.name]                      # whatever stopped, the raw file stands
    monkeypatch.setattr(cf, '_differs', real['differs'])
    monkeypatch.setattr(cf.os, 'rename', real['rename'])
    monkeypatch.setattr(cf, '_append_line', real['append'])
    monkeypatch.setattr(cf, '_remove_raw', real['remove'])
    result = cf.compress(case, floor_bytes=0)
    if stop in ('compressing', 'checking', 'renaming'):
        assert result['parts_removed'] == [first.name + cf.SUFFIX + cf.PART]
    if stop == 'removing':
        assert result['finished'] == [first.name] and result['compressed'] == ['cm1out_t000002_s.dat']
    else:
        assert result['compressed'] == [first.name, 'cm1out_t000002_s.dat']
    assert not result['failed'] and not result['left'] and not list(case.glob('*.part'))
    assert all(not (case / n).exists() for n in data)
    assert cf.verify(case)['ok'] == sorted(data)
    text = (case / cf.MANIFEST).read_text()
    assert all(json.loads(line)['name'] in data for line in text.splitlines()[-2:])
    again = cf.compress(case, floor_bytes=0)                            # a third run finds everything done
    assert sorted(again['done']) == sorted(data) and not again['compressed']


def test_a_raw_file_that_differs_from_its_line_is_left_with_its_zst(tmp_path):
    case = tmp_path / 'case'
    data = synthetic_case(case, outputs=(1,), restarts=())
    cf.compress(case, floor_bytes=0)
    path = case / 'cm1out_t000001_s.dat'
    path.write_bytes(b'another run wrote this')                         # beside the .zst of the first
    result = cf.compress(case, floor_bytes=0)
    assert [x['name'] for x in result['left']] == [path.name] and not result['compressed']
    assert path.read_bytes() == b'another run wrote this' and cf.compressed(path).exists()
    assert cf.manifest(case)[path.name]['raw_sha256'] == hashlib.sha256(data[path.name]).hexdigest()


def test_a_stopped_decompression_leaves_no_raw_file_and_its_part_goes(tmp_path):
    case = tmp_path / 'case'
    synthetic_case(case, outputs=(1,), restarts=())
    cf.compress(case, floor_bytes=0)
    (case / 'cm1out_t000001_s.dat.part').write_bytes(b'half a file')   # a restore stopped part-way
    assert cf.glob(case, 'cm1out_t*_s.dat') == [case / 'cm1out_t000001_s.dat']
    result = cf.compress(case, floor_bytes=0)
    assert result['parts_removed'] == ['cm1out_t000001_s.dat.part'] and result['done'] == ['cm1out_t000001_s.dat']


def test_the_free_space_floor_stops_the_work(tmp_path, monkeypatch):
    case = tmp_path / 'case'
    data = synthetic_case(case, outputs=(1, 2, 3), restarts=())
    size = len(data['cm1out_t000001_s.dat'])
    free = iter([10 ** 12, 30 * 10 ** 9 + size - 1])                   # room for one file, then below the floor
    monkeypatch.setattr(cf, 'free_bytes', lambda folder: next(free))
    with pytest.raises(cf.SpaceFloor) as stopped:
        cf.compress(case)
    assert stopped.value.result['compressed'] == ['cm1out_t000001_s.dat']
    assert (case / 'cm1out_t000002_s.dat').read_bytes() == data['cm1out_t000002_s.dat']
    assert not list(case.glob('*.part')) and not cf.compressed(case / 'cm1out_t000002_s.dat').exists()
    monkeypatch.setattr(cf, 'free_bytes', lambda folder: 30 * 10 ** 9 + size)    # exactly the floor plus the file
    assert cf.compress(case)['compressed'] == ['cm1out_t000002_s.dat', 'cm1out_t000003_s.dat']
    monkeypatch.setattr(cf, 'free_bytes', lambda folder: 0)
    with pytest.raises(cf.SpaceFloor):                                  # writing back raw keeps the floor too
        cf.decompress(case, ['cm1out_t000001_s.dat'])
    assert not (case / 'cm1out_t000001_s.dat').exists()


def test_the_listings_find_raw_and_compressed_files_alike(tmp_path):
    case = tmp_path / 'case'
    synthetic_case(case, outputs=(1, 2, 3), restarts=(4, 5, 6))
    for name in ('cm1out_t000002_s.dat', 'cm1rst_t000004_s.dat', 'cm1rst_t000004_x.dat', 'cm1rst_t000005_i.dat',
                 'cm1rst_t000005_s.dat', 'cm1rst_t000005_u.dat', 'cm1rst_t000005_v.dat', 'cm1rst_t000005_w.dat',
                 'cm1rst_t000005_x.dat'):
        assert cf.compress_file(case / name, floor_bytes=0)['status'] == 'compressed'
    (case / 'cm1rst_t000006_w.dat').unlink()                            # set 6 incomplete
    (case / 'cm1out_t000004_s.dat.zst.part').write_bytes(b'')           # being written
    assert cf.glob(case, 'cm1out_t*_s.dat') == [case / f'cm1out_t{n:06d}_s.dat' for n in (1, 2, 3)]
    assert [p.name for p in case.glob('cm1out_t*_s.dat')] != [p.name for p in cf.glob(case, 'cm1out_t*_s.dat')]
    assert cf.restart_numbers(case) == [4, 5]
    assert [p.name for p in cf.restart_files(case, 5)] == [f'cm1rst_t000005_{p}.dat' for p in 'isuvwx']
    assert cf.kept_names(case) == sorted(p.name for p in cf.glob(case, 'cm1*.dat') if cf.KEPT.fullmatch(p.name))
    assert cf.glob(tmp_path / 'absent', 'cm1out_t*_s.dat') == []
    assert cf.exists(case / 'cm1out_t000002_s.dat') and not cf.exists(case / 'cm1out_t000004_s.dat')
    assert cf.remove(case / 'cm1rst_t000004_s.dat') == ['cm1rst_t000004_s.dat.zst'] and cf.restart_numbers(case) == [5]


def test_a_case_in_use_is_left_alone(tmp_path, monkeypatch):
    case = tmp_path / 'case'
    data = synthetic_case(case, outputs=(1,), restarts=())
    (case / 'run.lock').write_text('1234')
    with pytest.raises(cf.InUse, match='run.lock'):
        cf.compress(case, floor_bytes=0)
    (case / 'run.lock').unlink()
    monkeypatch.setattr(cf, 'cm1_processes', lambda: [os.getpid()])     # a CM1 process working in the case
    monkeypatch.chdir(case)
    with pytest.raises(cf.InUse, match='runs in it'):
        cf.compress(case, floor_bytes=0)
    assert (case / 'cm1out_t000001_s.dat').read_bytes() == data['cm1out_t000001_s.dat']
    monkeypatch.chdir(tmp_path)
    assert cf.compress(case, floor_bytes=0)['compressed'] == ['cm1out_t000001_s.dat']


def test_the_command_line_compresses_verifies_and_decompresses(tmp_path, capsys):
    case = tmp_path / 'case'
    data = synthetic_case(case, outputs=(1, 2), restarts=(1,))
    assert cf.main(['compress', str(case), '--floor-gb', '0']) == 0
    assert '8 compressed' in capsys.readouterr().out
    assert cf.main(['verify', str(case)]) == 0 and '8 match' in capsys.readouterr().out
    assert cf.main(['decompress', str(case), 'cm1rst_t000001_s.dat', '--floor-gb', '0']) == 0
    assert (case / 'cm1rst_t000001_s.dat').read_bytes() == data['cm1rst_t000001_s.dat']
    (case / 'run.lock').write_text('1')
    assert cf.main(['compress', str(case), '--floor-gb', '0']) == 1 and 'skipped' in capsys.readouterr().out


def snapshot_case(folder, numbers=(1, 2)):
    """A case of four columns and three levels with GrADS control files, as ring_analysis reads them."""
    folder.mkdir(parents=True)
    names = [('th', 3), ('qv', 3), ('psfc', 0), ('t2', 0)]
    (folder / 'cm1out_s.ctl').write_text('xdef 4 linear 0 1\nzdef 3 levels\n1\n2\n3\n'
                                         f'vars {len(names)}\n' + ''.join(f'{n} {k} 99 {n}\n' for n, k in names) +
                                         'endvars\n')
    for n in numbers:
        values = np.arange(3 * 4 * 2 + 4 * 2, dtype='<f4') * (1.0 + n / 7.0)
        values.tofile(folder / f'cm1out_t{n:06d}_s.dat')


def test_ring_analysis_reads_a_compressed_snapshot_as_the_raw_one(tmp_path):
    from climate.crm import ring_analysis as ra
    raw_case, packed = tmp_path / 'raw', tmp_path / 'packed'
    snapshot_case(raw_case)
    snapshot_case(packed)
    cf.compress(packed, floor_bytes=0)
    assert not (packed / 'cm1out_t000002_s.dat').exists()
    for n in (1, 2):
        a, b = ra.read_snapshot(raw_case, n), ra.read_snapshot(packed, n)
        assert sorted(a) == sorted(b) == ['psfc', 'qv', 't2', 'th']
        assert all(a[k].shape == b[k].shape and a[k].tobytes() == b[k].tobytes() for k in a)


def test_the_wave_exporter_reads_a_compressed_snapshot_s_bytes_as_the_raw_one(tmp_path):
    from climate.crm import wave_forcing as wf
    layout = [('th', 3), *((k, 1) for k in wf.FIELDS), ('qv', 2)]
    nx = 5
    raw_case, packed = tmp_path / 'raw', tmp_path / 'packed'
    for folder in (raw_case, packed):
        folder.mkdir()
        values = np.concatenate([np.full(n * nx, 280.0 + k, '<f4') + np.arange(n * nx, dtype='<f4') / 100
                                 for k, (_, n) in enumerate(layout)])
        values.tofile(folder / 'cm1out_t000003_s.dat')
    cf.compress(packed, floor_bytes=0)
    a = wf.read_surface(raw_case / 'cm1out_t000003_s.dat', nx, layout)
    b = wf.read_surface(packed / 'cm1out_t000003_s.dat', nx, layout)
    assert a[1] == b[1] and all(a[0][k].tobytes() == b[0][k].tobytes() for k in wf.FIELDS)
    with pytest.raises(ValueError, match='inconsistent'):
        wf.read_surface(packed / 'cm1out_t000003_s.dat', nx + 1, layout)
