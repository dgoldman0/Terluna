"""Temporary byte-verified import; final tree excludes this helper and workflow."""
from __future__ import annotations
import base64
import hashlib
import json
import lzma
import os
from pathlib import Path
import subprocess
import sys
import tempfile

BASE = '79a60a1a51b4552ac9ae8e8bb295b567d9254995'
PAYLOAD_SHA256 = 'b836b4cc52028e13aa8aa49174d83a90ab5ba42251d9fe5adb423464778a0588'
DESTINATION = 'environment-import-result'
ROOT = Path.cwd().resolve()

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def sha(data):
    return hashlib.sha256(data).hexdigest()

def git(*args, data=None):
    return subprocess.run(['git', *args], input=data, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, check=True).stdout

def safe_path(relative):
    result = (ROOT / relative).resolve()
    require(result.is_relative_to(ROOT) and '.git' not in Path(relative).parts,
            'Unsafe path: ' + relative)
    return result

encoded = ''.join((ROOT / '.environment-import' / f'part-{i:02d}').read_text() for i in range(6))
require(sha(encoded.encode()) == PAYLOAD_SHA256, 'Transport digest mismatch')
payload = json.loads(lzma.decompress(base64.b64decode(encoded, validate=True)))
manifest = payload['manifest']
require(manifest['remote_base_commit'] == BASE, 'Unexpected base')
expected = manifest['files']
require(len(expected) == 37, 'Unexpected file count')
require(set(payload['files']).isdisjoint(payload['generate']), 'Duplicate transport/generated paths')
require(set(payload['files']) | set(payload['generate']) == set(expected), 'Incomplete manifest')
base_paths = set(git('ls-tree', '-r', '--name-only', BASE).decode().splitlines())
for path in expected:
    safe_path(path)
    if path in manifest['preimages']:
        data = git('show', f'{BASE}:{path}')
        old = manifest['preimages'][path]
        require(sha(data) == old['sha256'], 'Preimage hash mismatch: ' + path)
        require(hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()
                == old['git_blob_sha'], 'Preimage Git hash mismatch: ' + path)
    else:
        require(path not in base_paths, 'Unexpected existing path: ' + path)
for path, text in payload['files'].items():
    data = text.encode('utf-8')
    require(sha(data) == expected[path], 'Supplied file mismatch: ' + path)
    destination = safe_path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(data)

subprocess.run([sys.executable, '-m', 'unittest', 'discover', '-s', 'tests', '-v'], check=True)
with tempfile.TemporaryDirectory(prefix='terluna-generated-') as temp:
    generated = Path(temp)
    subprocess.run([sys.executable, 'research/run_environment_screens.py', '--out', str(generated)], check=True)
    mismatches = []
    for path in payload['generate']:
        override = ROOT / '.environment-import' / 'overrides' / path
        source = override if override.is_file() else generated / Path(path).name
        data = source.read_bytes()
        if sha(data) != expected[path]:
            mismatches.append({'path': path, 'expected': expected[path], 'actual': sha(data)})
        else:
            destination = safe_path(path)
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(data)
    if mismatches:
        print('GENERATED_HASH_MISMATCHES=' + json.dumps(mismatches, indent=2), flush=True)
        raise RuntimeError('Generated outputs differ; exact original bytes are required')

for path, digest in expected.items():
    require(sha(safe_path(path).read_bytes()) == digest, 'Final file mismatch: ' + path)
print('EXACT_FILES_VERIFIED=37', flush=True)
with tempfile.TemporaryDirectory(prefix='terluna-checks-') as temp:
    subprocess.run([sys.executable, 'research/check.py', '--output', str(Path(temp) / 'checks.json')], check=True)
subprocess.run([sys.executable, 'ensemble/tools/validate_setup.py'], check=True)

# Reset only the index to the original tree and stage exactly the intended files.
git('read-tree', BASE)
git('add', '--', *sorted(expected))
changed = set(git('diff', '--cached', '--name-only', BASE).decode().splitlines())
require(changed == set(expected), 'Final tree changed an unexpected path')
tree = git('write-tree').decode().strip()
message = ('Add bounded atmosphere, climate and long-night research models\n\n'
           'Import the prepared 37-file research change set with exact SHA-256 verification.\n'
           'Add numerical models, conditional results, proofs and 39 tests; rewrite the README\n'
           'and deprecate synthetic training as research evidence while preserving its history.\n'
           'Retain original baselines, accepted seeds and locked editorial controls unchanged.\n'
           'Original run/source records remain historical; complete physical validation is pending.\n')
commit = git('-c', 'user.name=github-actions[bot]', '-c',
             'user.email=41898282+github-actions[bot]@users.noreply.github.com',
             'commit-tree', tree, '-p', BASE, data=message.encode()).decode().strip()
require(not subprocess.run(['git', 'ls-remote', '--exit-code', '--heads', 'origin', DESTINATION],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0,
        'Result branch already exists; refusing to overwrite')
git('push', 'origin', f'{commit}:refs/heads/{DESTINATION}')
print('IMPORT_TREE=' + tree, flush=True)
print('IMPORT_COMMIT=' + commit, flush=True)
print('Destination branch only; main has not been updated.', flush=True)
