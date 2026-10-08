#!/usr/bin/env python3
"""Recheck this analysis package without changing Git refs, index or worktrees.

Only temporary Git objects are written, outside the repository. --out optionally
writes a new audit record to an explicit path. No models or test suites run.
"""
from __future__ import annotations
import argparse
import datetime
import itertools
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, help='Write a new audit JSON here; default prints summary only')
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    root = here.parents[1]
    snapshot = json.loads((here / 'evidence/snapshot.json').read_text())
    local = json.loads((here / 'evidence/local-state.json').read_text())
    failures = []
    checks = []
    def git(*argv, env=None):
        return subprocess.run(['git', '-C', str(root), '-c', 'gc.auto=0', *argv],
                              text=True, capture_output=True, env=env)
    def value(*argv):
        result = git(*argv)
        if result.returncode:
            raise RuntimeError(result.stderr)
        return result.stdout.strip()
    def check(name, ok, detail):
        checks.append(dict(check=name, passed=bool(ok), detail=detail))
        if not ok:
            failures.append(name)
    branches = snapshot['branches']
    for branch in branches:
        current = value('rev-parse', branch['ref'])
        check('ref '+branch['ref'], current == branch['sha'], dict(expected=branch['sha'], actual=current))
        base = value('merge-base', branches[0]['sha'], branch['sha'])
        counts = [int(x) for x in value('rev-list', '--left-right', '--count', branches[0]['sha']+'...'+branch['sha']).split()]
        paths = value('diff', '--name-only', base, branch['sha']).splitlines()
        check('ancestry '+branch['ref'], base == branch['main_merge_base'] and counts == branch['main_only_branch_only'] and paths == branch['changed_paths'], dict(base=base, counts=counts, changed_paths=len(paths)))
    for path, groups in local['preserved_trees'].items():
        if path == 'shared/provenance.py':
            continue  # Deliberately absent in the two older branches.
        for expected, refs in groups.items():
            for ref in refs:
                sha = next(b['sha'] for b in branches if b['ref'] == ref)
                actual = value('rev-parse', sha+':'+path)
                check('preserved '+ref+':'+path, actual == expected, actual)
    old = snapshot['excluded_merged_branch']
    check('excluded branch contained in main', git('merge-base', '--is-ancestor', old['sha'], branches[0]['sha']).returncode == 0, old['ref'])
    for worktree in snapshot['worktrees']:
        if worktree['branch'] == 'refs/heads/research/branch-integration-analysis':
            continue  # This package is intentionally uncommitted there.
        result = subprocess.run(['git', '-C', worktree['path'], 'status', '--porcelain=v1'], text=True, capture_output=True)
        check('source worktree '+worktree['branch'], result.returncode == 0 and result.stdout == worktree['status'], result.stdout or result.stderr or 'clean')
    objects = value('rev-parse', '--path-format=absolute', '--git-path', 'objects')
    with tempfile.TemporaryDirectory(prefix='terluna-merge-analysis-') as directory:
        env = os.environ.copy()
        env.update(GIT_OBJECT_DIRECTORY=directory, GIT_ALTERNATE_OBJECT_DIRECTORIES=objects, GIT_OPTIONAL_LOCKS='0')
        for a, b in itertools.combinations(branches, 2):
            expected = next(p for p in snapshot['pairwise_merges'] if p['left'] == a['ref'] and p['right'] == b['ref'])
            result = git('merge-tree', '--write-tree', '--name-only', a['sha'], b['sha'], env=env)
            lines = result.stdout.splitlines()
            end = lines.index('') if '' in lines else len(lines)
            conflicts = lines[1:end] if result.returncode == 1 else []
            base = value('merge-base', a['sha'], b['sha'])
            ap = set(value('diff', '--name-only', base, a['sha']).splitlines())
            bp = set(value('diff', '--name-only', base, b['sha']).splitlines())
            ok = result.returncode == expected['exit_code'] and conflicts == expected['conflict_paths'] and sorted(ap & bp) == expected['overlap_paths']
            check('pair '+a['ref']+' + '+b['ref'], ok, dict(exit_code=result.returncode, conflicts=conflicts, overlaps=len(ap & bp), stderr=result.stderr))
    for path in sorted((here / 'evidence').glob('*.json')):
        json.loads(path.read_text())
        check('JSON '+path.name, True, 'parsed')
    link_count = 0
    for path in sorted(here.glob('*.md')):
        for target in re.findall(r'\[[^\]]*\]\(([^\s)]+)\)', path.read_text()):
            if '://' in target or target.startswith('#'):
                continue
            target = target.split('#', 1)[0]
            if not target:
                continue
            link_count += 1
            check('link '+path.name+' -> '+target, (path.parent / target).exists(), target)
    result = dict(schema='terluna.research.integration-audit/1', checked_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  scope='Pinned sources and review artifacts only; no cumulative integration or scientific validation',
                  checks=checks, failures=failures, local_links_checked=link_count)
    if args.out:
        args.out.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(result='PASS' if not failures else 'FAIL', checks=len(checks), failures=failures, local_links=link_count), indent=2))
    return 1 if failures else 0

if __name__ == '__main__':
    raise SystemExit(main())
