#!/usr/bin/env python3
"""Check the immutable planning import and workspace; never grant manuscript clearance."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import runpy
import sys
import uuid


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def local_path(root: Path, relative: str) -> Path:
    candidate = (root / relative).resolve()
    require(candidate.is_relative_to(root.resolve()), f"Path escapes workspace: {relative}")
    return candidate


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate(root: Path) -> dict:
    root = root.resolve()
    snapshot = root / "planning"
    setup = read_json(root / "repository_setup.json")
    require(setup["snapshot_directory"] == "planning", "Unexpected snapshot location")
    checksums = read_json(snapshot / "checksums.json")
    require(isinstance(checksums, dict) and bool(checksums), "Missing checksum map")
    for relative, expected in checksums.items():
        require(sha256(local_path(snapshot, relative)) == expected,
                f"Snapshot checksum mismatch: {relative}")
    actual_files = {str(p.relative_to(snapshot)) for p in snapshot.rglob("*") if p.is_file()}
    expected_files = set(checksums) | {"checksums.json"}
    require(actual_files == expected_files, "Unexpected or missing planning snapshot files")
    require(len(actual_files) == setup["snapshot_file_count"], "Snapshot file count changed")

    lock = read_json(snapshot / "control_lock.json")
    require(sha256(local_path(snapshot, lock["charter"])) == lock["charter_sha256"],
            "Charter lock mismatch")
    for source in lock["source_controls"]:
        require(sha256(local_path(snapshot, source["path"])) == source["sha256"],
                f"Source-control lock mismatch: {source['name']}")

    manifest = read_json(snapshot / "paper_manifest.json")
    papers = manifest["papers"]
    require(manifest["paper_count"] == len(papers) == 5, "Expected five paper identities")
    require({p["role"] for p in papers} == {"root", "physical", "biosphere", "human", "engineering"},
            "Paper roles changed")
    ids = [p["paper_id"] for p in papers]
    require(len(set(ids)) == 5 and ids == lock["paper_ids"], "Duplicate or changed paper IDs")
    require(set(setup["workspaces"]) == set(ids), "Workspace mapping differs from paper register")
    outline = (snapshot / "Open_Moon_Ensemble_Outline_plan-01.md").read_text(encoding="utf-8")
    for paper in papers:
        identifier = uuid.UUID(paper["full_uuid"])
        require(identifier.version == 4 and str(identifier)[:13] == paper["paper_id"],
                f"Invalid UUID identity: {paper['paper_id']}")
        require(paper["paper_id"] in outline, f"Paper absent from outline: {paper['paper_id']}")
        workspace = local_path(root, setup["workspaces"][paper["paper_id"]])
        text = (workspace / "README.md").read_text(encoding="utf-8")
        require(paper["paper_id"] in text and paper["full_uuid"] in text and paper["title"] in text,
                f"Workspace identity mismatch: {paper['paper_id']}")
        for target in ("../../planning/Open_Moon_Ensemble_Outline_plan-01.md",
                       "../../planning/Open_Moon_Editorial_Charter_control-01.md"):
            require(target in text and (workspace / target).is_file(),
                    f"Workspace planning link missing: {paper['paper_id']}")

    diagnostic = runpy.run_path(str(snapshot / "twilight_diagnostic.py"))["diagnostic"]
    actual = diagnostic()
    expected = read_json(snapshot / "twilight_diagnostic.json")
    require(actual.keys() == expected.keys(), "Diagnostic output fields changed")
    for key, value in expected.items():
        require(math.isclose(actual[key], value, rel_tol=1e-12, abs_tol=1e-12),
                f"Diagnostic mismatch: {key}")
    return {"result": "PASS", "scope": "planning snapshot and workspace integrity only",
            "snapshot_files": len(actual_files), "paper_workspaces": len(papers),
            "scientific_or_manuscript_clearance": False}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ensemble-root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    try:
        result = validate(args.ensemble_root)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
