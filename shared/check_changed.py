#!/usr/bin/env python3
"""Select checks for the working tree, including consumers of recorded data products.

No scientific modules are imported while planning. Imports, literal file paths
and product metadata supply dependency edges; unfamiliar changes fall back to
their owner's tests or the full suite. This is a conservative static selector,
not proof of coverage for dynamically constructed or external dependencies.
"""
from __future__ import annotations

import argparse
import ast
from collections import defaultdict, deque
from dataclasses import asdict, dataclass, field
import fnmatch
from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import posixpath
import re
import shlex
import subprocess
import sys
import tempfile
import tomllib

ROOT = Path(__file__).resolve().parents[1]
DOMAINS = {
    "atmosphere", "biosphere", "climate", "engineering", "geography", "geology",
    "habitation", "humanities", "illumination", "protection", "provisioning", "resources",
}
JS_GROUPS = {
    "test-immersion": ("immersion/",),
    "test-column": ("atmosphere/column/",),
    "test-light-transport": ("illumination/references/",),
    "test-column-viewer": ("visualization/atmospheric-columns/",),
    "test-calendar": ("illumination/calendar/", "visualization/light-calendar/"),
}
CONTROL = {"Makefile", "shared/check_changed.py"}
GLOBAL = {"pyproject.toml", "conftest.py", "research/requirements.txt"}
# These are planning/check records, not executable inputs by default. Actual
# readers still participate in the graph; this only disables the owner fallback.
RECORDS = {
    "research/status.json", "research/checks.json", "research/tooling_checks.json",
    "habitation/crater_gardens/design.json", "habitation/crater_gardens/buildout.json",
}
CODE_SUFFIXES = {".py", ".js", ".mjs", ".cjs"}
PROSE_SUFFIXES = {".md", ".rst", ".txt", ".bib"}
SKIP = ("archive/", "immersion/docs/history/")
JS_STRINGS = re.compile(r'''["'`]([^"'`\n]+)["'`]''')


def git(root: Path, *args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=root).decode("utf-8", "surrogateescape")


def paths(output: str) -> set[str]:
    return set(filter(None, output.split("\0")))


def inventory(root: Path, base: str) -> tuple[set[str], set[str]]:
    current = paths(git(root, "ls-files", "-z", "--cached", "--others", "--exclude-standard"))
    previous = paths(git(root, "ls-tree", "-r", "--name-only", "-z", base))
    changed = paths(git(root, "diff", "--name-only", "--no-renames", "-z", base, "--"))
    changed |= paths(git(root, "ls-files", "--others", "--exclude-standard", "-z"))
    return current | previous, changed


def owner(path: str) -> str | None:
    parts = PurePosixPath(path).parts
    if parts[:2] == ("research", "studies") and len(parts) > 2:
        return "/".join(parts[:3])
    if parts and parts[0] in DOMAINS | {"shared", "immersion"}:
        return parts[0]
    if parts and parts[0] == "visualization" and len(parts) > 1:
        return "/".join(parts[:2])
    return None


class References:
    def __init__(self, files: set[str]):
        self.files = files
        self.extensions = {PurePosixPath(p).suffix for p in files} - {""}
        self.suffixes: dict[str, set[str]] = defaultdict(set)
        for path in files:
            parts = path.split("/")
            for n in range(len(parts)):
                self.suffixes["/".join(parts[n:])].add(path)

    @lru_cache(maxsize=32768)
    def resolve(self, value: str, source: str) -> set[str]:
        """Prefer exact/local paths; ambiguous suffixes include every candidate."""
        value = value.strip()
        if not value or len(value) > 512 or "://" in value or "\n" in value:
            return set()
        if PurePosixPath(value).suffix not in self.extensions and value not in self.suffixes:
            return set()
        rooted = value.startswith("./") or "/" in value
        value = posixpath.normpath(value)
        if rooted and value in self.files:
            return {value}
        parent = posixpath.dirname(source)
        while parent:
            local = posixpath.normpath(posixpath.join(parent, value))
            if local in self.files:
                return {local}
            parent = posixpath.dirname(parent)
        if ("*" in value or "?" in value) and PurePosixPath(value).suffix:
            # A bare '*.json' in a runtime directory does not read every JSON
            # file in its domain. Only resolve globs with a known directory.
            literal_prefix = re.split(r"[*?]", value, maxsplit=1)[0]
            if "/" not in literal_prefix or literal_prefix in {"./", "../"}:
                return set()
            local = posixpath.normpath(posixpath.join(posixpath.dirname(source), value))
            return {p for p in self.files if fnmatch.fnmatch(p, value) or fnmatch.fnmatch(p, local)}
        if not PurePosixPath(value).suffix:
            return set()  # e.g. an external model's Makefile is not this repository's Makefile.
        return self.suffixes.get(value.lstrip("./"), set())

    def strings(self, value: str, source: str) -> set[str]:
        found = set(self.resolve(value, source))
        # importlib.import_module and Python -m entry points.
        if re.fullmatch(r"[A-Za-z_]\w*(?:\.[A-Za-z_]\w*)+", value):
            found.update(self.module(value))
        return found

    @lru_cache(maxsize=16384)
    def module(self, name: str) -> set[str]:
        base = name.replace(".", "/")
        candidates = {base + ".py", base + "/__init__.py"}
        # Package initializers also execute when importing a child module.
        parts = base.split("/")
        candidates.update("/".join(parts[:n]) + "/__init__.py" for n in range(1, len(parts)))
        return candidates & self.files


def expression(node: ast.AST, names: dict[str, str]) -> str | None:
    """Resolve common Path expressions without evaluating repository code."""
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    if isinstance(node, ast.Name):
        return names.get(node.id)
    if isinstance(node, ast.JoinedStr):
        return "".join(n.value if isinstance(n, ast.Constant) else "*" for n in node.values)
    if isinstance(node, ast.BinOp):
        left, right = expression(node.left, names), expression(node.right, names)
        if isinstance(node.op, ast.Div) and right is not None:
            return posixpath.join(left or "", right)
        if isinstance(node.op, ast.Add) and left is not None and right is not None:
            return left + right
    if isinstance(node, ast.Call):
        if isinstance(node.func, ast.Name) and node.func.id in {"Path", "PurePath", "PurePosixPath"} and node.args:
            value = expression(node.args[0], names)
            if value and isinstance(node.args[0], ast.Constant) and "/" not in value:
                value = "./" + value
            return value
        if isinstance(node.func, ast.Attribute) and node.func.attr in {"resolve", "absolute"}:
            return expression(node.func.value, names)
    if isinstance(node, ast.Attribute) and node.attr == "parent":
        base = expression(node.value, names)
        return (posixpath.dirname(base) or ".") if base is not None else None
    if (isinstance(node, ast.Subscript) and isinstance(node.value, ast.Attribute)
            and node.value.attr == "parents" and isinstance(node.slice, ast.Constant)
            and isinstance(node.slice.value, int) and node.slice.value >= 0):
        base = expression(node.value.value, names)
        if base is not None:
            for _ in range(node.slice.value + 1):
                base = posixpath.dirname(base)
            return base or "."
    return None


def python_edges(source: str, text: str, refs: References, is_test: bool) -> dict[str, set[str]]:
    tree = ast.parse(text, filename=source)
    names = {"__file__": source}
    for node in tree.body:
        if isinstance(node, (ast.Assign, ast.AnnAssign)) and node.value is not None:
            value = expression(node.value, names)
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            if value is not None:
                for target in targets:
                    if isinstance(target, ast.Name):
                        names[target.id] = value
    result: dict[str, set[str]] = defaultdict(set)

    def visit(node: ast.AST, consumer: str, classes: tuple[str, ...] = (), reading=False) -> None:
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant) and isinstance(node.value.value, str):
            return  # Docstrings do not read the files they describe.
        if isinstance(node, ast.ClassDef):
            classes += (node.name,)
        if is_test and isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name.startswith("test_"):
            consumer = "::".join((source, *classes, node.name))
        if isinstance(node, ast.Call):
            function = node.func.attr if isinstance(node.func, ast.Attribute) else getattr(node.func, "id", "")
            reading = reading or function in {"read_text", "read_bytes", "open", "read", "load", "loadtxt", "read_csv",
                                              "exists", "is_file", "is_dir", "stat"}
            reading = reading or any(word in function.lower() for word in ("digest", "sha256", "hash_file"))
            if reading and isinstance(node.func, ast.Attribute):
                base = expression(node.func.value, names)
                if base is not None:
                    result[consumer].update(refs.strings(base, source))
            if function in {"glob", "rglob"} and isinstance(node.func, ast.Attribute) and node.args:
                directory = expression(node.func.value, names)
                pattern = expression(node.args[0], names)
                if directory and pattern:
                    result[consumer].update(refs.resolve(posixpath.join(directory, pattern), source))
        if isinstance(node, ast.Import):
            for alias in node.names:
                result[consumer].update(refs.module(alias.name))
        if isinstance(node, ast.ImportFrom):
            if node.level:
                parent = source.split("/")[:-node.level]
                name = ".".join(parent + ([node.module] if node.module else []))
            else:
                name = node.module or ""
            result[consumer].update(refs.module(name))
            for alias in node.names:
                result[consumer].update(refs.module(name + "." + alias.name))
        value = expression(node, names)
        if value is not None:
            result[consumer].update(p for p in refs.strings(value, source)
                                    if reading or PurePosixPath(p).suffix not in PROSE_SUFFIXES)
            if isinstance(node, (ast.BinOp, ast.JoinedStr, ast.Call)):
                return  # The full path supersedes its bare filename fragments.
        for child in ast.iter_child_nodes(node):
            visit(child, consumer, classes, reading)

    visit(tree, source)
    return result


def json_references(text: str, source: str, refs: References) -> set[str]:
    """Read file identities, not paths mentioned incidentally in descriptive prose.

Markdown citations are dependencies only when included in an input/producer
record or accompanied by a digest. Editing a cited decision must not invalidate
an otherwise unchanged scenario or historical result.
"""
    result: set[str] = set()
    identity_keys = {"inputs", "files", "dependencies", "bound_inputs", "hashes", "code"}
    citations = {"cited_files", "citations", "references", "source_notes", "evidence", "description", "notes", "material"}

    def visit(value, identity=False):
        if isinstance(value, dict):
            pinned = identity or (bool(value.keys() & {"sha256", "sha", "hash", "digest"})
                                  and bool(value.keys() & {"path", "file", "filename", "path_at_integration"}))
            for key, child in value.items():
                if key in citations:
                    continue
                context = pinned or key in identity_keys
                if context:
                    add(key, context)
                visit(child, context)
        elif isinstance(value, list):
            for child in value:
                visit(child, identity)
        elif isinstance(value, str):
            add(value, identity)

    def add(value, identity):
        # Exact references can contain spaces in filenames; narrative text cannot
        # become an edge merely because a regex finds a filename inside it.
        for path in refs.resolve(value, source):
            if identity or PurePosixPath(path).suffix not in PROSE_SUFFIXES:
                result.add(path)

    visit(json.loads(text))
    return result


def test_files(root: Path, files: set[str]) -> set[str]:
    config = tomllib.loads((root / "pyproject.toml").read_text())["tool"]["pytest"]["ini_options"]
    roots = config["testpaths"]
    excluded = config.get("norecursedirs", [])
    tests = set()
    for path in files:
        if not (PurePosixPath(path).name.startswith("test_") and path.endswith(".py")
                and (root / path).is_file()):
            continue
        for start in roots:
            if path.startswith(start.rstrip("/") + "/"):
                below = path[len(start.rstrip("/")) + 1:].split("/")[:-1]
                if not any(fnmatch.fnmatch(part, pattern) for part in below for pattern in excluded):
                    tests.add(path)
    return tests


def graph(root: Path, files: set[str], tests: set[str]) -> dict[str, set[str]]:
    refs = References(files)
    reverse: dict[str, set[str]] = defaultdict(set)
    # Cache dependency extraction only, never test success. Every source's bytes
    # are hashed again, and inventory/selector changes invalidate resolved paths.
    cache_path = root / ".pytest_cache" / "check_changed" / "dependencies.json"
    identity = hashlib.sha256(Path(__file__).read_bytes() +
                              json.dumps([sorted(files), sorted(tests)]).encode()).hexdigest()
    try:
        cache = json.loads(cache_path.read_text())
        records = cache["records"] if cache["identity"] == identity else {}
        if not isinstance(records, dict):
            records = {}
    except (OSError, ValueError, KeyError, TypeError):
        records = {}
    updated = False
    for path in sorted(files):
        if path.startswith(SKIP) or not (root / path).is_file():
            continue
        suffix = PurePosixPath(path).suffix
        if suffix not in CODE_SUFFIXES | {".json"}:
            continue
        raw = (root / path).read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        cached = records.get(path, {})
        if not isinstance(cached, dict):
            cached = {}
        valid_edges = (isinstance(cached.get("edges"), dict)
                       and all(isinstance(p, str) and isinstance(ds, list) and all(isinstance(d, str) for d in ds)
                               for p, ds in cached["edges"].items()))
        if cached.get("sha256") == digest and valid_edges:
            edges = cached["edges"]
        elif path == "shared/check_changed.py":
            # The selector's path tables describe routing; they do not consume
            # the named scientific inputs. Its own tests are selected explicitly.
            edges = {}
        elif suffix == ".py":
            text = raw.decode(errors="replace")
            edges = python_edges(path, text, refs, path in tests)
        else:
            text = raw.decode(errors="replace")
            dependencies: set[str] = set()
            if suffix == ".json":
                dependencies = json_references(text, path, refs)
            else:
                for match in JS_STRINGS.finditer(text):
                    dependencies.update(refs.strings(match.group(1), path))
            edges = {path: dependencies}
        if cached.get("sha256") != digest or not valid_edges:
            records[path] = {"sha256": digest, "edges": {p: sorted(ds) for p, ds in edges.items()}}
            updated = True
        for consumer, dependencies in edges.items():
            for dependency in dependencies:
                if dependency != consumer:
                    reverse[dependency].add(consumer)
    if updated:
        temporary = None
        try:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            with tempfile.NamedTemporaryFile(mode="w", dir=cache_path.parent, delete=False) as handle:
                temporary = Path(handle.name)
                json.dump({"identity": identity, "records": records}, handle)
            temporary.replace(cache_path)
        except OSError:
            pass  # A read-only cache must not prevent verification.
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
    refs.resolve.cache_clear()
    refs.module.cache_clear()
    for path in files:
        if PurePosixPath(path).name == "conftest.py":
            parent = posixpath.dirname(path)
            reverse[path].update(t for t in tests if not parent or t.startswith(parent + "/"))
    return reverse


@dataclass
class Plan:
    changed: list[str]
    python: list[str] = field(default_factory=list)
    targets: list[str] = field(default_factory=list)
    reasons: dict[str, list[str]] = field(default_factory=dict)
    uncovered: list[str] = field(default_factory=list)
    full: bool = False


def plan(root: Path, files: set[str], changed: set[str]) -> Plan:
    result = Plan(sorted(changed))
    if not changed:
        return result
    global_changes = {
        p for p in changed if p in GLOBAL or PurePosixPath(p).name.startswith("requirements") or
        (p.startswith("shared/") and p not in CONTROL and not p.startswith("shared/tests/")
         and PurePosixPath(p).suffix not in PROSE_SUFFIXES and p != "shared/check_layers.py")
    }
    if global_changes:
        result.full = True
        result.reasons["check-full"] = sorted(global_changes)
        return result
    tests = test_files(root, files)
    reverse = graph(root, files, tests)
    queue = deque(sorted(changed))
    parents: dict[str, str | None] = dict.fromkeys(changed)
    while queue:
        dependency = queue.popleft()
        for consumer in sorted(reverse.get(dependency, ())):
            if consumer not in parents:
                parents[consumer] = dependency
                queue.append(consumer)

    def chain(path: str) -> list[str]:
        result = [path]
        while parents.get(result[-1]) is not None:
            result.append(parents[result[-1]])
        return list(reversed(result))

    selected = {p for p in parents if p.split("::")[0] in tests}
    for target, prefixes in JS_GROUPS.items():
        affected = sorted(p for p in parents if p.startswith(prefixes)
                          and PurePosixPath(p).suffix not in PROSE_SUFFIXES)
        if affected:
            result.targets.append(target)
            result.reasons[target] = chain(affected[0])
    if any(p.startswith("ensemble/") for p in changed) or "ensemble/tools/validate_setup.py" in parents:
        result.targets.append("check-ensemble")
    if ("research/check.py" in parents or "research/provenance.json" in changed
            or any(p.startswith("research/baselines/") for p in changed)):
        result.targets.append("check-provenance")
    if any(PurePosixPath(p).suffix in CODE_SUFFIXES for p in changed):
        result.targets.append("check-layers")
    if changed & CONTROL:
        selected.add("shared/tests/test_check_changed.py")
        result.reasons["shared/tests/test_check_changed.py"] = sorted(changed & CONTROL)
    if "shared/check_layers.py" in changed:
        selected.add("shared/tests/test_check_layers.py")
    # When no known consumer reaches a changed executable/data file, do not
    # interpret that as evidence that it has no effect. Use its owner's suite.
    forward: dict[str, set[str]] = defaultdict(set)
    for dependency, consumers in reverse.items():
        for consumer in consumers:
            forward[consumer].add(dependency)
    reached = set(tests) | {p for p in forward if p.split("::")[0] in tests}
    reached.update(p for p in files if any(p.startswith(prefixes) for prefixes in JS_GROUPS.values()))
    reached.update({"research/check.py", "ensemble/tools/validate_setup.py"})
    queue = deque(reached)
    while queue:
        for dependency in forward.get(queue.popleft(), ()):
            if dependency not in reached:
                reached.add(dependency)
                queue.append(dependency)
    for path in sorted(changed):
        if (path in reached or path in CONTROL | RECORDS or path.startswith(("ensemble/", "archive/"))
                or PurePosixPath(path).suffix in PROSE_SUFFIXES or path == ".gitignore"
                or any(path.startswith(prefixes) for prefixes in JS_GROUPS.values())):
            continue
        scope = owner(path)
        local = {t for t in tests if scope and t.startswith(scope + "/")}
        if local:
            selected.update(local)
            for test in local:
                result.reasons.setdefault(test, [path, "owner fallback: " + scope, test])
        elif scope:
            result.uncovered.append(path)
        else:
            result.full = True
            result.reasons.setdefault("check-full", []).append(path)
    # A module selection subsumes any individually selected tests in it.
    result.python = sorted(p for p in selected if "::" not in p or p.split("::")[0] not in selected)
    for test in result.python:
        result.reasons.setdefault(test, chain(test))
    result.targets = sorted(set(result.targets))
    return result


def validate_changes(root: Path, base: str, changed: list[str]) -> None:
    subprocess.run(["git", "diff", "--check", base, "--"], cwd=root, check=True)
    for path in changed:
        file = root / path
        if not file.is_file():
            continue
        if file.suffix == ".json":
            json.loads(file.read_text())
        elif file.suffix == ".toml":
            tomllib.loads(file.read_text())
        elif file.suffix == ".py":
            ast.parse(file.read_text(), filename=path)


def commands(result: Plan, python: str) -> list[list[str]]:
    if result.full:
        return [["make", "check-full", "PYTHON=" + python]]
    jobs = [["make", target, "PYTHON=" + python] for target in result.targets]
    if result.python:
        jobs.append([python, "-m", "pytest", *result.python])
    return jobs


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", default="HEAD", help="Compare this revision with the working tree (default: HEAD)")
    parser.add_argument("--files", nargs="+", help="Explicit input paths, including ignored inputs if needed")
    parser.add_argument("--include", nargs="+", default=[], help="Additional affected inputs alongside the Git changes")
    parser.add_argument("--plan", action="store_true", help="Explain the selection without running checks")
    parser.add_argument("--json", action="store_true", help="Print the plan as JSON without running checks")
    args = parser.parse_args(argv)
    try:
        files, changed = inventory(ROOT, args.base)
        if args.files:
            changed = {posixpath.normpath(p) for p in args.files}
        changed.update(posixpath.normpath(p) for p in args.include)
        if any(p.startswith(("/", "../")) or p == ".." for p in changed):
            raise ValueError("Input paths must be relative to the repository")
        files |= changed
        selected = plan(ROOT, files, changed)
        if args.json:
            print(json.dumps(asdict(selected), indent=2))
            return 0
        print(f"Changed paths: {len(selected.changed)}; selected Python files/tests: {len(selected.python)}", flush=True)
        for job in commands(selected, sys.executable):
            print("  " + shlex.join(job), flush=True)
        for path in selected.uncovered:
            print(f"No configured automated coverage found: {path}; inspect the artifact and run its documented checks.", flush=True)
        if args.plan:
            for target, why in sorted(selected.reasons.items()):
                print(target + ": " + " -> ".join(why))
            return 0
        validate_changes(ROOT, args.base, selected.changed)
        env = os.environ.copy()
        for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
            env.setdefault(key, "1")
        for job in commands(selected, sys.executable):
            subprocess.run(job, cwd=ROOT, env=env, check=True)
        print("Selected checks passed. Full repository verification: " + ("performed." if selected.full else "not requested."))
        return 0
    except (OSError, ValueError, SyntaxError, subprocess.CalledProcessError) as exc:
        print(f"check-changed failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
