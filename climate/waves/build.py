"""Build stock SWAN 41.51 and an AGROW numerical-similarity sensitivity.

Run from the repository root with gfortran, gcc, make and Perl on PATH. Source
bytes stay in the explicitly chosen build directory; existing builds are never
overwritten. The patch removes one dimensional numerical default and does not
validate the wind-input parameterization for lunar conditions.
"""
from __future__ import annotations

import argparse
import difflib
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import tarfile
import time
import urllib.request

from shared.constants import STANDARD_GRAVITY

SOURCE_URL = "https://swanmodel.sourceforge.io/download/zip/swan4151.tar.gz"
SOURCE_SHA256 = "325c229f3dde0b812db4ec9ae6fa1cd4086d2109562c389576590ce8f94b26bb"
AGROW_SOURCE_SHA256 = "b164c7b9ac8f94672b637e114dbacad503b624c8c1c8ae50339b80c0083558af"
MAX_ARCHIVE_BYTES = 4 * 1024**2
MAX_EXTRACTED_BYTES = 64 * 1024**2


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def extract_source(archive: Path, destination: Path) -> Path:
    """Verify the pinned archive and extract regular files into a fresh directory."""
    if sha256(archive) != SOURCE_SHA256:
        raise ValueError("SWAN archive SHA256 differs from the pinned source")
    with tarfile.open(archive, "r:gz") as stream:
        members = stream.getmembers()
        names = set()
        total = 0
        for member in members:
            name = PurePosixPath(member.name)
            if (name.is_absolute() or ".." in name.parts or not name.parts
                    or name.parts[0] != "swan4151"
                    or not (member.isdir() or member.isfile())
                    or name in names):
                raise ValueError(f"Unsafe or duplicate SWAN archive member: {member.name}")
            names.add(name)
            total += member.size
        if total > MAX_EXTRACTED_BYTES:
            raise ValueError("SWAN archive exceeds the extraction size limit")
        destination.mkdir(parents=True, exist_ok=False)
        for member in members:
            target = destination.joinpath(*PurePosixPath(member.name).parts)
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                with stream.extractfile(member) as source, target.open("xb") as output:
                    shutil.copyfileobj(source, output)
    return destination / "swan4151"


def patch_agrow(source: Path) -> dict:
    """Keep the arbitrary AGROW frequency knee at fixed frequency/gravity."""
    path = source / "swancom3.ftn"
    before = path.read_bytes()
    if hashlib.sha256(before).hexdigest() != AGROW_SOURCE_SHA256:
        raise ValueError("swancom3.ftn is not the pinned original SWAN 41.51 file")
    text = before.decode("ascii")
    old = "            FREQ1 = SIGMA/PI2                                             40.88"
    if text.count(old) != 1:
        raise ValueError("SWIND0 FREQ1 anchor differs from original SWAN 41.51")
    new = (
        "!     Terluna sensitivity: keep the arbitrary AGROW knee at fixed f/g.\n"
        "!     This is a numerical choice, not lunar empirical validation.\n"
        f"            FREQ1 = (SIGMA/PI2) * ({STANDARD_GRAVITY:.9g}/GRAV)"
    )
    after = text.replace(old, new)
    difference = "".join(difflib.unified_diff(
        text.splitlines(keepends=True), after.splitlines(keepends=True),
        fromfile="a/swancom3.ftn", tofile="b/swancom3.ftn"))
    path.write_bytes(after.encode("ascii"))
    return dict(path=path.name, before_sha256=AGROW_SOURCE_SHA256,
                after_sha256=sha256(path), unified_diff=difference,
                reference_gravity_m_s2=STANDARD_GRAVITY,
                evidence="numerical similarity sensitivity; no lunar empirical validation")


def compile_source(source: Path, environment: dict) -> dict:
    """Use SWAN's compiler configuration and serial build targets."""
    start = time.monotonic()
    commands = [["make", "config"], ["make", "-j1", "ser"]]
    with (source / "build.log").open("x") as log:
        for command in commands:
            subprocess.run(command, cwd=source, env=environment, check=True,
                           stdout=log, stderr=subprocess.STDOUT, timeout=600)
    executable = source / "swan.exe"
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise RuntimeError(f"SWAN build did not produce an executable in {source}")
    macros = (source / "macros.inc").read_text()
    return dict(executable=str(executable), executable_sha256=sha256(executable),
                elapsed_wall_s=time.monotonic() - start, commands=commands,
                macros_inc=macros, macros_sha256=sha256(source / "macros.inc"),
                build_log_sha256=sha256(source / "build.log"))


def build(build_root: Path, *, archive: Path | None = None, download: bool = False) -> dict:
    """Prepare two independent source trees and record their completed builds."""
    if (archive is None) == (not download):
        raise ValueError("Choose exactly one of an archive path or explicit download")
    build_root = build_root.absolute()
    if build_root.exists() or build_root.is_symlink():
        raise FileExistsError(f"Build directory already exists: {build_root}")
    programs = {name: shutil.which(name) for name in ("gfortran", "gcc", "make", "perl")}
    missing = [name for name, path in programs.items() if path is None]
    if missing:
        raise RuntimeError("Missing build programs on PATH: " + ", ".join(missing))
    compiler_version = subprocess.check_output(
        [programs["gfortran"], "--version"], text=True, timeout=15).strip()
    if archive is not None:
        archive = archive.resolve(strict=True)
        if archive.stat().st_size > MAX_ARCHIVE_BYTES or sha256(archive) != SOURCE_SHA256:
            raise ValueError("SWAN archive differs from the pinned source")
    build_root.mkdir(parents=True, exist_ok=False)
    pinned_archive = build_root / "swan4151.tar.gz"
    if download:
        with urllib.request.urlopen(SOURCE_URL, timeout=60) as response:
            data = response.read(MAX_ARCHIVE_BYTES + 1)
        if len(data) > MAX_ARCHIVE_BYTES or hashlib.sha256(data).hexdigest() != SOURCE_SHA256:
            raise ValueError("Downloaded SWAN source differs from the pinned archive")
        pinned_archive.write_bytes(data)
    else:
        shutil.copyfile(archive, pinned_archive)
    pristine = extract_source(pinned_archive, build_root / "source")
    stock = build_root / "stock"
    patched = build_root / "patched_agrow"
    shutil.copytree(pristine, stock)
    shutil.copytree(pristine, patched)
    patch = patch_agrow(patched)
    (build_root / "agrow.patch").write_text(patch["unified_diff"])
    environment = dict(os.environ, FC="gfortran", OMP_NUM_THREADS="1",
                       OPENBLAS_NUM_THREADS="1")
    for name in ("MAKEFLAGS", "MFLAGS", "GNUMAKEFLAGS"):
        environment.pop(name, None)
    stock_build = compile_source(stock, environment)
    patched_build = compile_source(patched, environment)
    record = dict(schema="terluna.climate.swan-build/1", software="SWAN", version="41.51",
                  source_archive=dict(url=SOURCE_URL, sha256=sha256(pinned_archive),
                                      size_bytes=pinned_archive.stat().st_size),
                  producer_sha256=sha256(Path(__file__)),
                  compiler_version=compiler_version, programs=programs,
                  build_mode="serial", make_jobs=1,
                  stock_executable=dict(path="stock/swan.exe", sha256=sha256(stock / "swan.exe")),
                  patched_executable=dict(path="patched_agrow/swan.exe",
                                          sha256=sha256(patched / "swan.exe")),
                  patch_manifest=dict(reference_gravity_m_s2=STANDARD_GRAVITY,
                                      files=[{key: patch[key] for key in
                                              ("path", "before_sha256", "after_sha256")}],
                                      unified_diff=patch["unified_diff"],
                                      evidence=patch["evidence"]),
                  builds=dict(stock=stock_build, patched_agrow=patched_build))
    (build_root / "build.json").write_text(json.dumps(record, indent=2) + "\n")
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--archive", type=Path, help="Previously downloaded pinned SWAN archive")
    source.add_argument("--download", action="store_true", help="Download the pinned upstream archive")
    parser.add_argument("--build-root", type=Path, required=True, help="New build directory")
    args = parser.parse_args()
    record = build(args.build_root, archive=args.archive, download=args.download)
    print(json.dumps({name: info["executable"] for name, info in record["builds"].items()}, indent=2))


if __name__ == "__main__":
    main()
