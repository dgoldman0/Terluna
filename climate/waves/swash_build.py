"""Build the pinned official SWASH release in the external wave-data directory."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tarfile
import time
import urllib.request

from climate.waves.model import ROOT, sha256

URL = "https://swash.sourceforge.io/download/zip/swash-12.01.tar.gz"
DIGEST = "c8a2e149e4c1583f3f5a10fa6e9fa110e6db7d1eb41978420a4a387877ac397f"
RUNS = ROOT / "research/runs/waves/swash_source"


def build(directory=RUNS):
    directory = Path(directory).absolute()
    directory.mkdir(parents=True, exist_ok=True)
    archive = directory / "swash-12.01.tar.gz"
    if not archive.exists():
        with urllib.request.urlopen(URL, timeout=60) as response:
            content = response.read(4*1024**2)
        archive.write_bytes(content)
    if sha256(archive) != DIGEST:
        raise ValueError("SWASH archive differs from its pinned bytes")
    source = directory / "source/swash-12.01"
    if not source.exists():
        with tarfile.open(archive) as stream:
            stream.extractall(directory/"source", filter="data")
    manifest = directory / "build.json"
    executable = source / "swash.exe"
    with tarfile.open(archive) as stream:
        original={Path(member.name).name:hashlib.sha256(stream.extractfile(member).read()).hexdigest()
                  for member in stream if member.isfile()}
    # make config rewrites the compiler configuration; the scientific source
    # stays byte-identical to the pinned archive before and after compilation.
    for name,digest in original.items():
        if manifest.exists() and name=='macros.inc':
            continue
        if sha256(source/name)!=digest:
            raise ValueError(f'SWASH source differs from the pinned archive: {name}')
    if manifest.exists():
        record = json.loads(manifest.read_text())
        if (record['source_sha256']!=DIGEST or record['source_files']!=original or
                sha256(executable) != record["executable_sha256"]):
            raise ValueError("SWASH executable differs from its build record")
        return executable
    source_hashes = {p.name: sha256(p) for p in source.iterdir() if p.is_file()}
    started = time.monotonic()
    commands = [["make", "config"], ["make", "-j1", "ser", "FLAGS_OPT=-O3"]]
    environment = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    for key in ("MAKEFLAGS", "MFLAGS", "GNUMAKEFLAGS"):
        environment.pop(key, None)
    with (directory/"build.log").open("x") as log:
        for command in commands:
            subprocess.run(command, cwd=source, env=environment, stdout=log,
                           stderr=subprocess.STDOUT, check=True, timeout=900)
    record = dict(schema="terluna.climate.swash-build/1", version="12.01",
                  source_url=URL, source_sha256=DIGEST, source_files=source_hashes,
                  license="GNU GPL v3 or later; Delft University of Technology",
                  source_modifications=[], commands=commands,
                  compiler=subprocess.check_output(["gfortran","--version"], text=True).splitlines()[0],
                  executable_sha256=sha256(executable), log_sha256=sha256(directory/"build.log"),
                  elapsed_wall_s=time.monotonic()-started)
    manifest.write_text(json.dumps(record,indent=2)+"\n")
    print(json.dumps({"executable":str(executable),"elapsed_wall_s":record["elapsed_wall_s"]}),flush=True)
    return executable


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=RUNS)
    build(parser.parse_args().directory)
