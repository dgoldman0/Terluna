"""The ledger of the Terluna patches in an ExoPlaSim source folder (climate/gcm/patch_ledger.py): a patch applied is
recorded with its version, a file carrying another version of a patch or a marker with no record is refused, a folder
patched before the ledger is recorded only where each file is exactly the untouched file with the current patches, a
rebuild writes that file, and ensure_patched keeps the ledger and refuses a revised patch."""
from __future__ import annotations

import shutil
import sys
import types

import pytest

from climate.gcm import exoplasim_run as er
from climate.gcm import patch_ledger as plg

UNTOUCHED = "      program demo\n      x = 1\n      y = 2\n      end program demo\n"
A1 = ('a.f90', 'Terluna: patch A', [('      x = 1\n', '      x = 1   ! Terluna: patch A\n      x = x + 1\n')])
A2 = ('a.f90', 'Terluna: patch A', [('      x = 1\n', '      x = 1   ! Terluna: patch A\n      x = x + 2\n')])
B = ('a.f90', 'Terluna: patch B', [('      y = 2\n', '      y = 2   ! Terluna: patch B\n      y = y * 3\n')])


def folders(tmp_path, text):
    """An untouched folder holding a.f90 as UNTOUCHED and a source folder holding a.f90 as text."""
    untouched, src = tmp_path / 'untouched', tmp_path / 'src'
    untouched.mkdir()
    src.mkdir()
    (untouched / 'a.f90').write_text(UNTOUCHED)
    (src / 'a.f90').write_text(text)
    return untouched, src


def test_a_patch_recorded_as_applied_passes_the_check():
    ledger = {}
    plg.note(ledger, *A1)
    plg.check(ledger, '/nowhere', *A1)
    assert ledger == {'a.f90': {'Terluna: patch A': plg.fingerprint(A1[2])}}


def test_a_revised_patch_is_refused():
    ledger = {}
    plg.note(ledger, *A1)
    with pytest.raises(RuntimeError, match='another version'):
        plg.check(ledger, '/nowhere', *A2)


def test_a_marker_with_no_record_is_refused():
    with pytest.raises(RuntimeError, match='no record of its version'):
        plg.check({}, '/nowhere', *A1)


def test_the_fingerprint_follows_every_edit():
    assert plg.fingerprint(A1[2]) != plg.fingerprint(A2[2])
    assert plg.fingerprint(A1[2]) == plg.fingerprint([tuple(e) for e in A1[2]])


def test_record_takes_a_file_that_is_the_untouched_one_with_the_current_patches(tmp_path):
    text = plg.apply(plg.apply(UNTOUCHED, A2[2], 'a.f90'), B[2], 'a.f90')
    untouched, src = folders(tmp_path, text)
    result = plg.record(src, untouched, [A2, B])
    assert result == {'a.f90': dict(ok=True, markers=['Terluna: patch A', 'Terluna: patch B'], reason=None)}
    ledger = plg.read(src)
    plg.check(ledger, src, *A2)
    plg.check(ledger, src, *B)


def test_record_leaves_a_file_with_an_earlier_version_unrecorded(tmp_path):
    untouched, src = folders(tmp_path, plg.apply(UNTOUCHED, A1[2], 'a.f90'))
    result = plg.record(src, untouched, [A2, B])
    assert not result['a.f90']['ok'] and result['a.f90']['markers'] == ['Terluna: patch A']
    with pytest.raises(RuntimeError, match='no record of its version'):
        plg.check(plg.read(src), src, *A2)


def test_rebuild_writes_the_untouched_file_with_the_current_patches(tmp_path):
    untouched, src = folders(tmp_path, plg.apply(plg.apply(UNTOUCHED, A1[2], 'a.f90'), B[2], 'a.f90'))
    plg.rebuild(src, untouched, 'a.f90', [A2, B])
    assert (src / 'a.f90').read_text() == plg.apply(plg.apply(UNTOUCHED, A2[2], 'a.f90'), B[2], 'a.f90')
    plg.check(plg.read(src), src, *A2)
    assert plg.verify(src, untouched, [A2, B])['a.f90']['ok']


def test_a_file_carrying_no_marker_needs_no_record(tmp_path):
    untouched, src = folders(tmp_path, UNTOUCHED)
    assert plg.verify(src, untouched, [A2, B]) == {}


def test_the_catalogue_holds_every_patch_a_run_applies_in_the_same_order_within_each_file():
    catalogue = [(n, m) for n, m, _ in er.patch_catalogue()]
    assert len(set(catalogue)) == len(catalogue)
    for switches in (dict(), dict(twilight=True, land_store=True), dict(twilight=True, gust=True, similarity=True),
                     dict(twilight=True, night=True, cloud_water=True, rain_drops=True, humidity_mean=True),
                     dict(twilight=True, depth=True, cap=True, gust=True, similarity=True, rain_drops=True,
                          humidity_mean=True, land_store=True)):
        chain = [(n, m) for n, m, _ in er.patch_chain(**switches)]
        assert set(chain) <= set(catalogue)
        for name in {n for n, _ in chain}:
            places = [catalogue.index(e) for e in dict.fromkeys(e for e in chain if e[0] == name)]
            assert places == sorted(places)


@pytest.mark.skipif(not (er.UNTOUCHED_SOURCE / 'radmod.f90').exists(), reason='no untouched ExoPlaSim sources at hand')
def test_ensure_patched_keeps_the_ledger_and_refuses_a_revised_patch(tmp_path, monkeypatch):
    package = tmp_path / 'exoplasim'
    src = package / 'plasim' / 'src'
    src.mkdir(parents=True)
    (package / 'plasim' / 'run').mkdir()
    for name in ('radmod.f90', 'rainmod.f90'):
        shutil.copy(er.UNTOUCHED_SOURCE / name, src / name)
    monkeypatch.setitem(sys.modules, 'exoplasim', types.SimpleNamespace(__file__=str(package / '__init__.py')))
    first = er.ensure_patched()
    ledger = plg.read(src)
    assert sorted(ledger) == ['radmod.f90', 'rainmod.f90']
    assert sum(len(v) for v in ledger.values()) == len(er.PATCHES)
    assert er.ensure_patched() == first                       # the second run finds every patch current
    name, marker, edits = er.PATCHES[0]
    old, new = edits[0]
    revised = [(name, marker, [(old, new + '!     revised\n')] + list(edits[1:]))] + er.PATCHES[1:]
    monkeypatch.setattr(er, 'PATCHES', revised)
    with pytest.raises(RuntimeError, match='another version'):
        er.ensure_patched()
    plg.rebuild(src, er.UNTOUCHED_SOURCE, name, er.patch_catalogue())
    assert '!     revised\n' in (src / name).read_text()
    er.ensure_patched()
