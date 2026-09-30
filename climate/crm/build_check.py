"""Whether a faster CM1 build reproduces moon_omp: one model day from copies of a finished case's newest restart files,
run with each build side by side on the same threads, and with moon_omp from restart files nudged in their last bit,
whose spread is the tolerance (the author's decision of 2026-09-29, in the decisions register).

    climate/gcm/.venv/bin/python -m climate.crm.build_check ring_equator moon_omp_o3 moon_omp_o3_strict --threads 2

writes ../results/crm/build_check_<case>.json; the copies run in runs/build_check/.

A build whose restart files match moon_omp's byte for byte after the day gives the same results. One whose files
differ is held to the nudge: its snapshots may differ from moon_omp's no more than the nudged runs' do, field by field.
A compiler's changes touch the last bit of every operation everywhere at once, so the nudge does the same to the state:
it flips the last bit of every value of the scalar restart file that is neither zero nor a whole number (so no mask or
count changes), and the second nudge of every other such value. A last-bit change moves individual storms within days,
as the runs' random starting perturbations do, and the analyses compare lunar-day statistics, so a build within that
spread changes nothing they report. A nudge of a single value can land on a field the model recomputes and may not
reach the whole ring within the day.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import numpy as np
from climate.crm import cm1_run as c
from climate.crm import ring_analysis as ra

CHECK = c.RUNS / 'build_check'
REFERENCE = 'moon_omp'
INPUTS = ('case.json', 'input_grid_z', 'input_sounding', 'LANDUSE.TBL', 'lsnudge_0001.dat', 'namelist.template',
          'RRTMG_LW_DATA', 'RRTMG_SW_DATA', 'terluna_surface.txt', 'terluna_wls2d.txt', 'perts.dat', 'cm1out_s.ctl',
          'cm1out_w.ctl')
SURFACE = ('t2', 'q2', 'psfc', 'prate', 's10', 'tsk')
SECTIONS = ('th', 'qv', 'uinterp', 'winterp')
NUDGES = 2                                  # nudged copies of the reference, each a different value's last bit


def nudge(path: Path, which: int) -> dict:
    """Flip in place the last bit of every value of a restart file that is neither zero nor a whole number (which 0),
    or of every other such value (which 1)."""
    raw = np.fromfile(path, dtype='<u4')
    values = raw.view('<f4')
    fit = np.isfinite(values) & (values != 0.0) & (values != np.round(values))
    if which:
        fit &= np.arange(values.size) % 2 == 0
    raw[fit] ^= np.uint32(1)
    raw.tofile(path)
    return dict(file=path.name, values=int(fit.sum()), of=int(values.size), every_other=bool(which))


def copy_case(case: Path, name: str, label: str, restart: int, days: float, which: int | None = None) -> dict:
    """A copy of case that runs days more from restart set restart with build label, nudged when which is given."""
    out = CHECK / name
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    for item in INPUTS:
        src = case / item
        if src.is_symlink():
            (out / item).symlink_to(src.resolve())
        elif src.exists():
            shutil.copy2(src, out / item)
    for part in 'isuvwx':
        shutil.copy2(case / f'cm1rst_t{restart:06d}_{part}.dat', out / f'cm1rst_t{restart:06d}_{part}.dat')
    (out / 'cm1.exe').symlink_to(c.CM1_HOME / 'build' / label / 'cm1.exe')
    text = (out / 'namelist.template').read_text()
    text = c.set_namelist(text, 'param2', 'irst', 1)
    text = c.set_namelist(text, 'param2', 'rstnum', restart)
    text = c.set_namelist(text, 'param1', 'run_time', round(days * 86400.0, 3))
    (out / 'namelist.input').write_text(text)
    return dict(name=name, label=label, folder=out,
                nudge=None if which is None else nudge(out / f'cm1rst_t{restart:06d}_s.dat', which))


def cm1_seconds(log: Path) -> float:
    """The run's own time as CM1 reports it at the end of its log."""
    line = [l for l in log.read_text(errors='replace').splitlines() if 'Total time' in l][-1]
    return float(line.split(':')[1].split()[0])


def run_side_by_side(copies: list, threads: int) -> None:
    """Run every copy at once on threads each, noting the time CM1 reports for each; any abnormal end stops the
    check."""
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_STACKSIZE='512M')
    running = []
    for copy in copies:
        log = open(copy['folder'] / 'cm1.log', 'w')
        running.append((copy, log, subprocess.Popen(['nice', '-n', '5', './cm1.exe'], cwd=copy['folder'], env=env,
                                                     stdout=log, stderr=subprocess.STDOUT, preexec_fn=c._unlimited_stack)))
    for copy, log, process in running:
        process.wait()
        log.close()
        if process.returncode != 0 or 'Program terminated normally' not in (copy['folder'] / 'cm1.log').read_text()[-4000:]:
            raise RuntimeError(f"CM1 stopped abnormally in {copy['folder']}")
        copy['wall_s'] = cm1_seconds(copy['folder'] / 'cm1.log')


def rms(a, b) -> float:
    return float(np.sqrt(np.mean((np.asarray(a, float) - np.asarray(b, float)) ** 2)))


def differences(copy: Path, reference: Path) -> dict:
    """Root-mean-square and domain-mean differences from the reference, field by field, over the snapshots both wrote."""
    numbers = sorted(int(p.name[8:14]) for p in reference.glob('cm1out_t*_s.dat'))
    numbers = [n for n in numbers if (copy / f'cm1out_t{n:06d}_s.dat').exists()]
    out = {}
    for key in (*SURFACE, *SECTIONS):
        pairs = [(ra.read_snapshot(copy, n)[key], ra.read_snapshot(reference, n)[key]) for n in numbers]
        out[key] = dict(rms=[rms(a, b) for a, b in pairs], mean=[float(np.mean(a) - np.mean(b)) for a, b in pairs])
    return dict(snapshots=numbers, fields=out)


def restart_bytes_match(copy: Path, reference: Path, restart: int) -> bool:
    return all(hashlib.sha256((copy / f'cm1rst_t{restart:06d}_{part}.dat').read_bytes()).digest()
               == hashlib.sha256((reference / f'cm1rst_t{restart:06d}_{part}.dat').read_bytes()).digest()
               for part in 'isuvwx')


def verdict(build: dict, nudged: list) -> dict:
    """Identical, or for each field whether the build's root-mean-square difference, averaged over the snapshots, is
    no larger than the largest of the nudged runs'."""
    if build['identical']:
        return dict(identical=True, within_nudge=True)
    fields = {}
    for key, value in build['differences']['fields'].items():
        own = float(np.mean(value['rms']))
        spread = max(float(np.mean(n['differences']['fields'][key]['rms'])) for n in nudged)
        fields[key] = dict(rms=own, nudge_rms=spread, ratio=own / spread if spread > 0 else float('inf'))
    return dict(identical=False, within_nudge=all(f['ratio'] <= 1.0 for f in fields.values()), fields=fields)


def check(case_name: str, labels, threads: int = 2, days: float = 1.0, nudges: int = NUDGES) -> dict:
    case = c.RUNS / case_name
    restart = c.latest_restart(case)
    if not restart:
        raise RuntimeError(f'{case_name} has no restart files')
    copies = [copy_case(case, f'{case_name}_{REFERENCE}', REFERENCE, restart, days)]
    copies += [copy_case(case, f'{case_name}_{label}', label, restart, days) for label in labels]
    copies += [copy_case(case, f'{case_name}_{REFERENCE}_nudge{k + 1}', REFERENCE, restart, days, which=k) for k in range(nudges)]
    run_side_by_side(copies, threads)
    reference = copies[0]['folder']
    restart_s = json.loads((case / 'case.json').read_text())['configuration']['restart_s']
    last = restart + int(round(days * 86400.0 / restart_s))
    builds = {}
    for copy in copies[1:]:
        builds[copy['name']] = dict(label=copy['label'], nudge=copy['nudge'], wall_s=copy['wall_s'],
                                    speed_vs_reference=copies[0]['wall_s'] / copy['wall_s'],
                                    identical=restart_bytes_match(copy['folder'], reference, last),
                                    differences=differences(copy['folder'], reference))
    nudged = [b for b in builds.values() if b['nudge']]
    for b in builds.values():
        if not b['nudge']:
            b['verdict'] = verdict(b, nudged)
    executables = {label: json.loads((c.CM1_HOME / 'build' / label / 'build.json').read_text())
                   for label in (REFERENCE, *labels)}
    return dict(case=case_name, restart=restart, days=days, threads=threads, runs_at_once=len(copies),
                reference=dict(label=REFERENCE, wall_s=copies[0]['wall_s']),
                builds_used={label: dict(optimization=e.get('optimization', '-O2'), executable_sha256=e['executable_sha256'])
                             for label, e in executables.items()},
                runs=builds)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('case')
    parser.add_argument('labels', nargs='+')
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--days', type=float, default=1.0)
    parser.add_argument('--nudges', type=int, default=NUDGES)
    args = parser.parse_args(argv)
    result = check(args.case, args.labels, args.threads, args.days, args.nudges)
    result['producer'] = dict(domain='climate', files={'crm/build_check.py': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]})
    path = ra.RESULTS / f'build_check_{args.case}.json'
    path.write_text(json.dumps(result, indent=1) + '\n')
    for name, run in result['runs'].items():
        v = run.get('verdict', {})
        print(f"{name}: {run['wall_s']:.0f} s ({run['speed_vs_reference']:.2f} x the reference's speed)"
              + ('' if run['nudge'] else f"; identical {v.get('identical')}, within the nudge {v.get('within_nudge')}"))
    return 0


if __name__ == '__main__':
    sys.exit(main())
