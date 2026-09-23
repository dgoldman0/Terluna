#!/usr/bin/env python3
"""Check that code imports only across the lanes the repository allows.

Lanes, by top-level folder:
  shared         the foundation: constants and conventions every lane may use
  domain         atmosphere, climate, biosphere, protection, engineering,
                 geography, habitation, illumination
  research       the hub and its cross-domain studies
  visualization  scientific and engineering rendering of domain results
  immersion-engine  immersion/engine: rendering and runtime systems for any world
  immersion-world   immersion/world: world definitions (landscape, sites, weather, flora)
  immersion         the experiences and tests, which choose a world and run the engine
  immersion-bake    immersion/bake: turns domain products into runtime assets
  ensemble, archive

Rules (source lane -> lanes it may import):
  shared -> shared
  domain -> shared, domain
  research -> shared, domain, research
  visualization -> anything except archive (it may measure the experience)
  immersion-engine -> shared, immersion-engine   (the engine never imports a world)
  immersion-world -> shared, immersion-engine, immersion-world
  immersion -> shared and every immersion part except bake (reads baked assets)
  immersion-bake -> shared, every immersion part, domain
  ensemble -> ensemble;  archive -> nothing;  nothing imports archive (archived code does
  not run and is not scanned; live code is checked for imports into it)

Python `import`/`from` statements naming a repository top-level folder, and
JavaScript import/require/new URL specifiers with relative paths, are checked.
Historical harnesses in immersion/docs/history are skipped.

Usage: python shared/check_layers.py   (exit 1 on any violation)
"""
from __future__ import annotations

import ast
import re
import subprocess
import sys
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
DOMAINS = {"atmosphere", "climate", "biosphere", "protection", "engineering",
           "geography", "habitation", "illumination"}
ALLOWED = {
    "shared": {"shared"},
    "domain": {"shared", "domain"},
    "research": {"shared", "domain", "research"},
    "visualization": {"shared", "domain", "research", "visualization", "immersion", "immersion-engine",
                      "immersion-world", "immersion-bake"},
    "immersion-engine": {"shared", "immersion-engine"},
    "immersion-world": {"shared", "immersion-engine", "immersion-world"},
    "immersion": {"shared", "immersion-engine", "immersion-world", "immersion"},
    "immersion-bake": {"shared", "immersion-engine", "immersion-world", "immersion", "immersion-bake", "domain"},
    "ensemble": {"ensemble"},
    "archive": set(),
}
SKIP_PREFIXES = ("immersion/docs/history/", "immersion/node_modules/", "immersion/dist/", "archive/")
JS_SPECIFIER = re.compile(
    r"""(?:\bimport\s+(?:[^'";]*?\bfrom\s*)?|\bimport\s*\(\s*|\brequire\s*\(\s*|\bnew\s+URL\s*\(\s*)(['"])([^'"]+)\1""")


def lane(path: str) -> str | None:
    parts = PurePosixPath(path).parts
    if not parts:
        return None
    top = parts[0]
    if top == "immersion":
        sub = parts[1] if len(parts) > 1 else ""
        return {"bake": "immersion-bake", "engine": "immersion-engine", "world": "immersion-world"}.get(sub, "immersion")
    if top in DOMAINS:
        return "domain"
    return top if top in ALLOWED else None


def allowed(source_path: str, target_path: str) -> bool:
    """Whether code at source_path may import target_path (repository-relative)."""
    source, target = lane(source_path), lane(target_path)
    return source is None or target is None or target in ALLOWED[source]


def tracked_sources() -> list[str]:
    out = subprocess.run(["git", "ls-files", "*.py", "*.js", "*.mjs", "*.cjs"], cwd=ROOT,
                         capture_output=True, text=True, check=True).stdout.split()
    return [p for p in out if not p.startswith(SKIP_PREFIXES)]


def python_targets(path: str) -> list[tuple[int, str]]:
    try:
        tree = ast.parse((ROOT / path).read_text(), filename=path)
    except SyntaxError:
        return []
    found = []
    for node in ast.walk(tree):
        names = []
        if isinstance(node, ast.Import):
            names = [a.name for a in node.names]
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            names = [node.module]
        for name in names:
            top = name.split(".")[0]
            if (ROOT / top).is_dir() and lane(top):
                found.append((node.lineno, top))
    return found


def js_targets(path: str) -> list[tuple[int, str]]:
    text = (ROOT / path).read_text(errors="ignore")
    found = []
    for match in JS_SPECIFIER.finditer(text):
        spec = match.group(2)
        if not spec.startswith("."):
            continue  # packages and node: builtins
        target = (ROOT / path).parent.joinpath(spec).resolve()
        try:
            rel = target.relative_to(ROOT).as_posix()
        except ValueError:
            found.append((text.count("\n", 0, match.start()) + 1, "<outside repository>"))
            continue
        found.append((text.count("\n", 0, match.start()) + 1, rel))
    return found


def main() -> int:
    violations = []
    checked = 0
    for path in tracked_sources():
        source = lane(path)
        if source is None:
            continue
        checked += 1
        targets = python_targets(path) if path.endswith(".py") else js_targets(path)
        for line, target in targets:
            target_lane = lane(target) if target != "<outside repository>" else "outside"
            if target_lane is None:
                continue
            if target_lane not in ALLOWED[source]:
                violations.append(f"{path}:{line}: {source} imports {target} ({target_lane})")
    for v in violations:
        print("LAYER VIOLATION", v)
    print(f"check_layers: {checked} files checked, {len(violations)} violations")
    return 1 if violations else 0


if __name__ == "__main__":
    sys.exit(main())
