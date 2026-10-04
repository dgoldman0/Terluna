"""The shared constants a data product depends on, recorded by name and value.

A product that pins the whole of constants.json goes stale whenever any lane
adds a constant. A producer records instead the constants its pinned Python
files read, found from their imports of shared.constants, with their values.
The product's tests compare those values with the current ones, so it stays
current through unrelated additions and shows any change to a value it used.

    from shared.provenance import constants_used, constants_changed
    producer = dict(files={...}, constants=constants_used(files))
    assert not constants_changed(product['producer']['constants'])
"""
from __future__ import annotations

import ast
import json
from pathlib import Path

from shared import constants

ROOT = Path(__file__).resolve().parents[1]


def constants_read(files):
    """Names from shared.constants that the given Python files import or read as attributes."""
    names = set()
    for file in files:
        path = Path(file) if Path(file).is_absolute() else ROOT / file
        if path.suffix != ".py":
            continue
        tree = ast.parse(path.read_text(), filename=str(path))
        aliases, dotted = set(), False
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.module == "shared.constants":
                names.update(a.name for a in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module == "shared":
                aliases.update(a.asname or a.name for a in node.names if a.name == "constants")
            elif isinstance(node, ast.Import):
                for a in node.names:
                    if a.name == "shared.constants":
                        if a.asname:
                            aliases.add(a.asname)
                        else:
                            dotted = True
        for node in ast.walk(tree):
            if not isinstance(node, ast.Attribute):
                continue
            base = node.value
            if isinstance(base, ast.Name) and base.id in aliases:
                names.add(node.attr)
            elif (dotted and isinstance(base, ast.Attribute) and base.attr == "constants"
                  and isinstance(base.value, ast.Name) and base.value.id == "shared"):
                names.add(node.attr)
    if names & {"DATA", "*"}:
        raise ValueError("Read named constants; DATA or a star import would pin the whole table")
    unknown = sorted(n for n in names if not hasattr(constants, n))
    if unknown:
        raise ValueError(f"Not in shared.constants: {unknown}")
    return sorted(names)


def _plain(value):
    return json.loads(json.dumps(value))


def constants_used(files):
    """{name: value} of every shared constant the given Python files read."""
    return {name: _plain(getattr(constants, name)) for name in constants_read(files)}


def constants_changed(recorded):
    """{name: (recorded, current)} for each recorded constant whose current value differs or is gone."""
    changed = {}
    for name, value in (recorded or {}).items():
        current = _plain(getattr(constants, name)) if hasattr(constants, name) else None
        if current != value:
            changed[name] = (value, current)
    return changed
