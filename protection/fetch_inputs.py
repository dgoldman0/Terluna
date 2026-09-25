#!/usr/bin/env python3
"""Restore exact external protection inputs; never change the recorded hashes."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parent


def restore(archive: Path | None = None, download: bool = False) -> dict:
    manifest = json.loads((ROOT / "inputs.json").read_text(encoding="utf-8"))
    target = ROOT / manifest["input_directory"]
    pending: list[tuple[Path, bytes]] = []
    checked = []
    with zipfile.ZipFile(archive) if archive is not None else _NoArchive() as source:
        for item in manifest["files"]:
            name = item["name"]
            if Path(name).name != name or name in (".", ".."):
                raise ValueError(f"Unsafe input name: {name}")
            path = target / name
            if path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest() == item["sha256"]:
                checked.append(name)
                continue
            if source is not None and "archive_member" in item:
                info = source.getinfo(item["archive_member"])
                if info.file_size != item["bytes"]:
                    raise ValueError(f"Archive size differs for {name}")
                data = source.read(info)
            elif download:
                req = urllib.request.Request(item["url"], headers={"User-Agent": "Terluna-input-restorer/1"})
                with urllib.request.urlopen(req, timeout=30) as response:
                    data = response.read(item["bytes"] + 1)
            else:
                raise FileNotFoundError(f"Missing or changed input: {path}. Use --archive or --download.")
            if len(data) != item["bytes"] or hashlib.sha256(data).hexdigest() != item["sha256"]:
                raise ValueError(f"Input checksum/size mismatch: {name}; upstream or archive data need review.")
            pending.append((path, data))
            checked.append(name)
    # Validate every input before replacing any destination.
    target.mkdir(parents=True, exist_ok=True)
    for path, data in pending:
        with tempfile.NamedTemporaryFile(dir=target, delete=False) as temporary:
            temporary.write(data)
            temp_path = Path(temporary.name)
        try:
            temp_path.replace(path)
        finally:
            temp_path.unlink(missing_ok=True)
    return {"status": "PASS", "scope": "exact input bytes only", "checked": checked,
            "restored": len(pending), "scientific_source_admission": False}


class _NoArchive:
    def __enter__(self):
        return None
    def __exit__(self, *args):
        return False


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--archive", type=Path)
    group.add_argument("--download", action="store_true")
    args = parser.parse_args()
    try:
        print(json.dumps(restore(args.archive, args.download), indent=2))
        return 0
    except (OSError, ValueError, KeyError, zipfile.BadZipFile) as exc:
        print(f"INPUTS BLOCKED: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
