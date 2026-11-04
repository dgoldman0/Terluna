"""Build the ocean models from their pinned sources on the research drive.

MITgcm builds unchanged, once for each experiment's compile-time code (its
grid size and options), since its lunar settings are run-time inputs. SCHISM
builds once with schism/lunar.patch applied: its gravity and rotation are
compiled constants, and the Moon's tide replaces Earth's constituents. Every
build checks that the checkout sits at the commit pinned in sources.json and
that the NetCDF packages match their hashes, and SCHISM's build keeps a ledger
of the patch and the files it changed beside the binary.

The drive holds src/ (checkouts), deps/ (NetCDF), builds/ and runs/; set
TERLUNA_OCEAN_DRIVE to move it. `python -m climate.ocean.builds schism` builds
SCHISM (about ten minutes on three cores, ParMETIS included).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCES = json.loads((HERE / "sources.json").read_text())
PATCH = HERE / "schism" / "lunar.patch"
DRIVE = Path(os.environ.get("TERLUNA_OCEAN_DRIVE", "/media/projectspace/terluna-research/ocean"))
SYSTEM_NETCDF = Path("/usr/lib/x86_64-linux-gnu/libnetcdf.so.19")
MPI_INCLUDE = Path("/usr/lib/x86_64-linux-gnu/openmpi/include")
CMAKE = Path(os.environ.get("TERLUNA_OCEAN_VENV", "/media/projectspace/terluna-research/venvs/ocean")) / "bin" / "cmake"
CHECKOUT = {"mitgcm": "MITgcm", "schism": "schism"}


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def files_digest(paths) -> str:
    """One digest over named files: each relative name and its content."""
    digest = hashlib.sha256()
    for name, path in sorted(paths.items()):
        digest.update(name.encode() + b"\0" + Path(path).read_bytes() + b"\0")
    return digest.hexdigest()


def run(command, cwd=None, log=None, env=None):
    """Run a command, appending its output to `log`; raise with the log's tail on failure."""
    with open(log, "a") if log else open(os.devnull, "w") as out:
        out.write(f"\n$ {' '.join(map(str, command))}\n")
        out.flush()
        done = subprocess.run([str(c) for c in command], cwd=cwd, stdout=out, stderr=subprocess.STDOUT, env=env)
    if done.returncode != 0:
        tail = Path(log).read_text().splitlines()[-25:] if log else []
        raise RuntimeError(f"{command[0]} failed ({done.returncode}):\n" + "\n".join(tail))


def checkout(name: str) -> Path:
    """The pinned checkout of `name`, cloned at its tag if absent; refuses any other commit."""
    pin = SOURCES["models"][name]
    path = DRIVE / "src" / CHECKOUT[name]
    if not (path / ".git").exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        run(["git", "clone", "--depth", "1", "--branch", pin["tag"], pin["url"], path], log=DRIVE / "logs" / f"clone_{name}.log")
    head = subprocess.run(["git", "-C", str(path), "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
    if head != pin["commit"]:
        raise RuntimeError(f"{path} is at {head}, not the pinned {pin['tag']} {pin['commit']}")
    return path


def netcdf_prefix() -> Path:
    """Headers and libraries for NetCDF-Fortran: the pinned packages unpacked, linked to the system's C runtime."""
    deps = DRIVE / "deps"
    prefix = deps / "prefix"
    if (prefix / "lib" / "libnetcdff.so").exists() and (prefix / "lib" / "libnetcdf.so").exists():
        return prefix
    if not SYSTEM_NETCDF.exists():
        raise RuntimeError(f"{SYSTEM_NETCDF} (the NetCDF C runtime) is missing")
    deps.mkdir(parents=True, exist_ok=True)
    for package in SOURCES["libraries"]["packages"]:
        deb = deps / package["url"].rsplit("/", 1)[1]
        if not deb.exists():
            urllib.request.urlretrieve(package["url"], deb)
        if sha256(deb) != package["sha256"]:
            raise RuntimeError(f"{deb.name}: sha256 differs from sources.json")
        run(["dpkg-deb", "-x", deb, deps / "root"])
    lib = prefix / "lib"
    lib.mkdir(parents=True, exist_ok=True)
    multiarch = deps / "root" / "usr" / "lib" / "x86_64-linux-gnu"
    for link, target in (("libnetcdff.so", multiarch / "libnetcdff.so.7.1.0"),
                         ("libnetcdff.so.7", multiarch / "libnetcdff.so.7.1.0"),
                         ("libnetcdf.so", SYSTEM_NETCDF)):
        (lib / link).unlink(missing_ok=True)
        (lib / link).symlink_to(target)
    include = prefix / "include"
    if not include.exists():
        include.symlink_to(deps / "root" / "usr" / "include")
    return prefix


def schism_build_dir() -> Path:
    return DRIVE / "builds" / f"schism-{SOURCES['models']['schism']['tag']}-lunar-{sha256(PATCH)[:12]}"


def schism(jobs: int = 3) -> Path:
    """The lunar SCHISM binary, built if absent: pinned source, lunar.patch, NetCDF and ParMETIS."""
    source = checkout("schism")
    build = schism_build_dir()
    ledger_path = build / "terluna_ledger.json"
    found = sorted((build / "build" / "bin").glob("pschism*")) if (build / "build" / "bin").exists() else []
    if found and ledger_path.exists():
        ledger = json.loads(ledger_path.read_text())
        if ledger["patch_sha256"] == sha256(PATCH) and ledger["upstream_commit"] == SOURCES["models"]["schism"]["commit"]:
            return found[0]
    prefix = netcdf_prefix()
    if build.exists():
        shutil.rmtree(build)
    build.mkdir(parents=True)
    log = build / "build.log"
    for part in ("src", "cmake"):
        shutil.copytree(source / part, build / part, symlinks=True)
    with open(PATCH) as patch:
        done = subprocess.run(["patch", "-p1", "-d", str(build)], stdin=patch, capture_output=True, text=True)
    if done.returncode != 0:
        raise RuntimeError("lunar.patch did not apply:\n" + done.stdout + done.stderr)
    changed = sorted({line.split()[1][2:] for line in PATCH.read_text().splitlines() if line.startswith("+++ ")})
    ledger = dict(schema="terluna.climate.ocean-schism-build/1", upstream_tag=SOURCES["models"]["schism"]["tag"],
                  upstream_commit=SOURCES["models"]["schism"]["commit"], patch="climate/ocean/schism/lunar.patch",
                  patch_sha256=sha256(PATCH),
                  patched_files={name: dict(upstream_sha256=sha256(source / name), built_sha256=sha256(build / name))
                                 for name in changed})
    ledger_path.write_text(json.dumps(ledger, indent=1) + "\n")
    (build / "build").mkdir()
    env = dict(os.environ, NETCDF_C_DIR=str(prefix), NETCDF_FORTRAN_DIR=str(prefix))
    run([CMAKE if CMAKE.exists() else shutil.which("cmake"),
         "-C", build / "cmake/SCHISM.local.build", "-C", build / "cmake/SCHISM.local.ubuntu",
         f"-DNetCDF_C_DIR={prefix}", f"-DNetCDF_FORTRAN_DIR={prefix}", f"-DNetCDF_ROOT={prefix}",
         f"-DNetCDF_C_ROOT={prefix}", f"-DNetCDF_Fortran_ROOT={prefix}",
         "-DCMAKE_POLICY_VERSION_MINIMUM=3.5", build / "src"],   # ParMETIS 4.0.3 predates CMake 4's minimum
        cwd=build / "build", log=log, env=env)
    run(["make", f"-j{jobs}", "pschism"], cwd=build / "build", log=log, env=env)
    found = sorted((build / "build" / "bin").glob("pschism*"))
    if not found:
        raise RuntimeError(f"no pschism binary in {build / 'build' / 'bin'}")
    ledger["binary"] = dict(name=found[0].name, sha256=sha256(found[0]))
    ledger_path.write_text(json.dumps(ledger, indent=1) + "\n")
    return found[0]


def mitgcm(code: dict[str, str], *, mpi: bool = True, jobs: int = 3) -> Path:
    """An MITgcm binary for one experiment's compile-time files (`code`: file name to text), built if absent."""
    source = checkout("mitgcm")
    digest = hashlib.sha256(json.dumps(dict(code=code, mpi=mpi, commit=SOURCES["models"]["mitgcm"]["commit"]),
                                       sort_keys=True).encode()).hexdigest()[:12]
    build = DRIVE / "builds" / f"mitgcm-{SOURCES['models']['mitgcm']['tag']}-{digest}"
    exe = build / "mitgcmuv"
    if exe.exists():
        return exe
    if build.exists():
        shutil.rmtree(build)
    (build / "code").mkdir(parents=True)
    for name, text in code.items():
        (build / "code" / name).write_text(text)
    log = build / "build.log"
    env = dict(os.environ, MPI_INC_DIR=str(MPI_INCLUDE)) if mpi else dict(os.environ)
    command = [source / "tools" / "genmake2", "-rootdir", source, "-mods", build / "code",
               "-optfile", source / "tools" / "build_options" / "linux_amd64_gfortran"] + (["-mpi"] if mpi else [])
    run(command, cwd=build, log=log, env=env)
    run(["make", "depend"], cwd=build, log=log, env=env)
    run(["make", f"-j{jobs}"], cwd=build, log=log, env=env)
    if not exe.exists():
        raise RuntimeError(f"no mitgcmuv in {build}")
    (build / "terluna_build.json").write_text(json.dumps(dict(
        schema="terluna.climate.ocean-mitgcm-build/1", upstream_tag=SOURCES["models"]["mitgcm"]["tag"],
        upstream_commit=SOURCES["models"]["mitgcm"]["commit"], mpi=mpi, code_digest=digest,
        binary_sha256=sha256(exe)), indent=1) + "\n")
    return exe


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("model", choices=["schism", "netcdf"])
    parser.add_argument("--jobs", type=int, default=3)
    args = parser.parse_args()
    print(schism(args.jobs) if args.model == "schism" else netcdf_prefix())


if __name__ == "__main__":
    main()
