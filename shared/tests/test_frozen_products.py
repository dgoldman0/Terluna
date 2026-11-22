"""The frozen products' manifest against Git and the files at hand, and the rule that skips a test for an absent one."""
import hashlib
import subprocess

import pytest

from shared import frozen_products as fp

ENTRIES = fp.entries()


def git(*args):
    return subprocess.run(['git', '-C', str(fp.ROOT), *args], capture_output=True, text=True)


needs_git = pytest.mark.skipif(git('rev-parse', '--git-dir').returncode != 0, reason='not a Git checkout')


def test_the_manifest_lists_each_product_once_with_its_size_digest_and_blob():
    paths = [e['path'] for e in ENTRIES]
    assert paths and len(set(paths)) == len(paths) and paths == sorted(paths)
    for e in ENTRIES:
        assert set(e) == {'path', 'bytes', 'sha256', 'git_blob'}
        assert e['bytes'] > 0 and len(e['sha256']) == 64 and len(e['git_blob']) == 40
        assert not e['path'].startswith('/') and '..' not in e['path'].split('/')


@needs_git
def test_no_frozen_product_is_tracked_and_every_one_is_ignored():
    paths = [e['path'] for e in ENTRIES]
    assert git('ls-files', '--', *paths).stdout.split() == []
    ignored = git('check-ignore', '--no-index', '--', *paths).stdout.split()
    assert sorted(ignored) == paths


@needs_git
def test_history_holds_each_product_at_its_recorded_size():
    for e in ENTRIES:
        size = git('cat-file', '-s', e['git_blob'])
        if size.returncode:
            pytest.skip("this clone's history lacks the frozen products' blobs")
        assert int(size.stdout) == e['bytes'], e['path']


@needs_git
def test_the_smallest_product_comes_out_of_history_with_its_digest():
    e = min(ENTRIES, key=lambda entry: entry['bytes'])
    if git('cat-file', '-e', e['git_blob']).returncode:
        pytest.skip("this clone's history lacks the frozen products' blobs")
    assert hashlib.sha256(fp.from_history(e)).hexdigest() == e['sha256']


def test_every_product_at_hand_has_its_recorded_bytes():
    states = {e['path']: fp.state(e) for e in ENTRIES}
    assert [p for p, s in states.items() if s == 'differs'] == []


def test_an_absent_product_is_recognised_from_the_error_that_opening_it_raises(tmp_path):
    e = ENTRIES[0]
    try:
        open(tmp_path / e['path'], 'rb')
    except FileNotFoundError as error:
        direct = error
    assert fp.absent_frozen_product(direct, root=tmp_path) == e['path']
    try:
        try:
            open(tmp_path / e['path'], 'rb')
        except FileNotFoundError as inner:
            raise RuntimeError('could not load the returns') from inner
    except RuntimeError as error:
        wrapped = error
    assert fp.absent_frozen_product(wrapped, root=tmp_path) == e['path']


def test_other_missing_files_and_other_errors_are_left_as_failures(tmp_path):
    try:
        open(tmp_path / 'some/other/file.json', 'rb')
    except FileNotFoundError as error:
        other = error
    assert fp.absent_frozen_product(other, root=tmp_path) is None
    assert fp.absent_frozen_product(ValueError('a wrong value'), root=tmp_path) is None
    present = tmp_path / ENTRIES[0]['path']
    present.parent.mkdir(parents=True)
    present.write_bytes(b'{}')
    try:
        open(tmp_path / ENTRIES[0]['path'] / 'inside', 'rb')
    except (FileNotFoundError, NotADirectoryError) as error:
        beside = error
    assert fp.absent_frozen_product(beside, root=tmp_path) is None


def test_a_link_to_nothing_counts_as_absent_and_a_changed_file_as_differing(tmp_path):
    e = dict(ENTRIES[0])
    path = tmp_path / e['path']
    path.parent.mkdir(parents=True)
    path.symlink_to(tmp_path / 'nowhere')
    assert fp.state(e, root=tmp_path) == 'absent'
    assert e['path'] in fp.absent_paths(root=tmp_path)
    path.unlink()
    path.write_bytes(b'not the product')
    assert fp.state(e, root=tmp_path) == 'differs'
