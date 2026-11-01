"""Build the pinned general SCALE-SDM morphology release in an external cache.

The two source changes are standard Fortran logical syntax and configurable
liquid-settling gravity. Upstream code remains in the cache, with its BSD license.
See references/scale_bsd_2_clause.txt for the Team SCALE copyright and license.
"""
from __future__ import annotations

import argparse
import difflib
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tarfile
import time
import urllib.request

REVISION = "4eaf4f74aa03d091d9778eff373b816f178a962f"
ARCHIVE_SHA256 = "5cff5e454b8d7b22d7049f785b8afe4f584e5fa25580f915115d878b746779ff"
SOURCE_URL = f"https://codeload.github.com/Shima-Lab/SCALE-SDM_mixed-phase_Shima2019/tar.gz/{REVISION}"
CASE_REL = Path("scale-les/test/case/warmbubble/2D_Khain04_sdm")
PATCH_IDENTITIES = {
    "contrib/SDM/scale_atmos_phy_mp_sdm.F90": (
        "f1c5aa0b3c9f7522ff7283ccad4a311d7979abf31e2c6cbfe82ca83f5d493c89",
        "7dd8d44b34cbccba201ed75021f39435bba07045fdb91ce13db2003f73fc1605"),
    "contrib/SDM/sdm_motion.f90": (
        "7e4161af4c5bce443ae33fbb3243e9696f8401118d80cf694abb7dcfcf86bbbd",
        "1fdf72e7c01a8749c5075e9814673c649515a7fb0406608cf1f109522531c4a1"),
}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def patch_text(relative_path, original):
    """Reject any unidentified source before applying the reviewed change."""
    before_hash, after_hash = PATCH_IDENTITIES[relative_path]
    if hashlib.sha256(original.encode()).hexdigest() != before_hash:
        raise ValueError(f"Source identity mismatch: {relative_path}")
    if relative_path.endswith("scale_atmos_phy_mp_sdm.F90"):
        result, count = re.subn(r"==\s*(\.(?:true|false)\.)", r".eqv. \1", original, flags=re.I)
        assert count == 5
    else:
        start = original.index("  subroutine sdm_getvz_liq(")
        end = original.index("  end subroutine sdm_getvz_liq")
        chunk = original[start:end].replace(
            "    real(RP), parameter :: grav = 9.80665E+2_RP       ! gravity acceleration [cm/s^2]",
            "    real(RP) :: grav                               ! configured gravity [cm/s^2]")
        chunk = chunk.replace("    if( .not. doprecipitation ) then",
                              "    grav = grav_mks * m2cm\n\n    if( .not. doprecipitation ) then")
        result = original[:start] + chunk + original[end:]
    if hashlib.sha256(result.encode()).hexdigest() != after_hash:
        raise ValueError(f"Unexpected patch result: {relative_path}")
    return result


def prepare(cache, archive=None, gravity_patch=True):
    """Create a fresh cache; an archive mismatch is an error, never a fallback."""
    cache = Path(cache).resolve()
    cache.mkdir(parents=True, exist_ok=True)
    if (cache / "work").exists() or (cache / "pristine").exists():
        raise FileExistsError("Use a fresh cache to preserve existing builds and runs")
    if archive is None:
        archive = cache / "source.tar.gz"
        with urllib.request.urlopen(SOURCE_URL, timeout=60) as response:
            archive.write_bytes(response.read())
    archive = Path(archive).resolve()
    if sha256(archive) != ARCHIVE_SHA256:
        raise ValueError("SCALE source archive SHA256 mismatch")
    pristine = cache / "pristine"
    with tarfile.open(archive) as src:
        src.extractall(pristine, filter="data")
    source = pristine / f"SCALE-SDM_mixed-phase_Shima2019-{REVISION}"
    work = cache / "work"
    shutil.copytree(source, work)
    patches = []
    for rel, hashes in PATCH_IDENTITIES.items():
        if not gravity_patch and rel.endswith("sdm_motion.f90"):
            continue
        path = work / rel
        before = path.read_text()
        after = patch_text(rel, before)
        path.write_text(after)
        diff = "".join(difflib.unified_diff(before.splitlines(True), after.splitlines(True),
                                           fromfile="a/"+rel, tofile="b/"+rel))
        patch = cache / (path.name + ".patch")
        patch.write_text(diff)
        patches.append(dict(path=rel, original_sha256=hashes[0], patched_sha256=hashes[1],
                            patch=patch.name, patch_sha256=sha256(patch)))
    record = dict(revision=REVISION, url=SOURCE_URL, archive_sha256=sha256(archive),
                  bytes=archive.stat().st_size, patches=patches,
                  license_sha256=sha256(work / "LICENSE"), gravity_patch=gravity_patch)
    (cache / "source.json").write_text(json.dumps(record, indent=2)+"\n")
    return work


def build(work, netcdf_include, netcdf_libs, workers=2):
    if workers not in (1, 2):
        raise ValueError("This feasibility runner permits one or two build workers")
    work = Path(work).resolve()
    command = ["make", f"-j{workers}", "build", "ENABLE_SDM=T", "SCALE_SYS=Linux64-gnu-ompi",
               "GNUCPPOPT=-cpp", "MODDIROPT=-J",
               "FFLAGS_FAST=-O2 -g -cpp -m64 -fconvert=big-endian -frecord-marker=4 "
               "-fallow-argument-mismatch -ffree-line-length-none -fimplicit-none "
               "-fmodule-private -fbacktrace -fcheck=bounds -ffpe-trap=invalid,zero,overflow",
               "FFLAGS_CORE=-m64 -fconvert=big-endian -frecord-marker=4 -fallow-argument-mismatch",
               "NETCDF_INCLUDE="+netcdf_include, "NETCDF_LIBS="+netcdf_libs]
    env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
    start = time.monotonic()
    logpath = work.parent / "build.log"
    with logpath.open("x") as log:
        proc = subprocess.run(command, cwd=work / CASE_REL, env=env,
                              stdout=log, stderr=subprocess.STDOUT)
    record = dict(command=command, cwd=str(work / CASE_REL), exit_code=proc.returncode,
                  elapsed_s=time.monotonic()-start, log_sha256=sha256(logpath),
                  compiler=subprocess.check_output(["gfortran", "--version"], text=True).splitlines()[0],
                  environment={k: env[k] for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS")})
    if proc.returncode == 0:
        record["executable_sha256"] = {name: sha256(work / CASE_REL / name)
                                        for name in ("scale-les_init", "scale-les")}
    (work.parent / "build.json").write_text(json.dumps(record, indent=2)+"\n")
    if proc.returncode:
        raise RuntimeError(f"SCALE build failed; see {logpath}")
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--archive", type=Path)
    parser.add_argument("--netcdf-include", required=True, help="e.g. -I/path/to/include")
    parser.add_argument("--netcdf-libs", required=True, help="link flags, including any rpath")
    parser.add_argument("--workers", type=int, choices=(1, 2), default=2)
    parser.add_argument("--earth-portability-only", action="store_true")
    args = parser.parse_args()
    work = prepare(args.cache, args.archive, not args.earth_portability_only)
    print(json.dumps(build(work, args.netcdf_include, args.netcdf_libs, args.workers), indent=2))


if __name__ == "__main__":
    main()
