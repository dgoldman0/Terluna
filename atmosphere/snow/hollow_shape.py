"""Local interface to archived 2026 crystal shape functions, not a growth law.

Original coefficients, data and native executables stay in the ignored cache.
Both archived demonstration conventions are reproduced before interpreting the
component. The public wrapper rejects the source's silent aspect-ratio clipping
and the converted coordinate's excursion beyond the original fitting interval.
"""
from __future__ import annotations
import argparse
import ctypes
import hashlib
import io
import json
import math
from pathlib import Path
import shutil
import subprocess
import zipfile
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CACHE = ROOT / 'research/runs/snow_scale/hollow_growth'
STEM = 'harrington-pokrifka-revisiting-theories-for-the-growth-of-single-crystalline-ice-2026'
ARCHIVE = CACHE / (STEM + '.zip')
ARCHIVE_SHA = '3fa016d36ae11dad221b2c9b300a5fe928ed253ac92dd8acdb2887291f32bc36'
DRIVER = HERE / 'hollow_shape_driver.f90'
BUILD = CACHE / 'native_build'
SOURCE = 'hfunctions-jas-2026.f90'
FLAGS = ['-O0', '-std=legacy', '-fno-fast-math', '-fcheck=all',
         '-finit-real=snan', '-ffpe-trap=invalid,zero,overflow', '-ffree-line-length-none']
FIELDS = ('prism_edge', 'basal_edge', 'prism_middle', 'basal_at_prism_middle',
          'prism_at_basal_middle', 'basal_middle', 'prism_capacitance', 'basal_capacitance')
TABLES = ('hfunction-eqns-hex-jas2026.dat', 'hfunction-eqns-hexconvQ-jas2026.dat')
# Native output column positions, including phi in column zero.
COLUMNS = {'cylinder': (2, 1, 4, 3, 6, 5, 15, 16),
           'hexagonal_prism': (11, 12, 10, 14, 9, 13, 8, 7)}
NATIVE_COSINE = .866  # Preserve the published coordinate conversion exactly.


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def producer():
    return {str(p.relative_to(ROOT)): sha(p) for p in (Path(__file__), DRIVER)}


def archive_members(archive=ARCHIVE):
    if sha(archive) != ARCHIVE_SHA:
        raise ValueError('Hollow-growth archive hash mismatch')
    with zipfile.ZipFile(archive) as z:
        return {Path(n).name: z.read(n) for n in z.namelist() if not n.endswith('/')}


def build(output=BUILD, *, compiler='gfortran'):
    """Compile unchanged source demo and a one-parameter coordinate variant."""
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    members = archive_members()
    original = members[SOURCE].decode()
    if original.count('parameter(iwd = 0)') != 1 or original.count('end program hfunctions') != 1:
        raise ValueError('Unexpected original program boundaries')
    compiler_path = shutil.which(compiler)
    if compiler_path is None:
        raise FileNotFoundError(compiler)
    compiler_path = str(Path(compiler_path).resolve())
    identity = producer()
    output.mkdir(parents=True)
    (output / 'native_routines.f90').write_text(original.split('end program hfunctions', 1)[1])
    shutil.copyfile(DRIVER, output / DRIVER.name)
    commands = [[compiler_path, *FLAGS, '-shared', '-fPIC', 'native_routines.f90',
                 DRIVER.name, '-o', 'libhollow_shape.so']]
    for iwd in (0, 1):
        source = original if iwd == 0 else original.replace('parameter(iwd = 0)', 'parameter(iwd = 1)')
        (output / f'native_demo_{iwd}.f90').write_text(source)
        commands.append([compiler_path, *FLAGS, f'native_demo_{iwd}.f90', '-o', f'native_demo_{iwd}'])
    with (output / 'build.log').open('w') as log:
        for command in commands:
            subprocess.run(command, cwd=output, stdout=log, stderr=subprocess.STDOUT, check=True)
        for iwd in (0, 1):
            subprocess.run([str(output / f'native_demo_{iwd}')], cwd=output, stdout=log,
                           stderr=subprocess.STDOUT, check=True)
    if producer() != identity:
        raise RuntimeError('Wrapper changed during build')
    record = dict(schema='terluna.atmosphere.hollow-shape-build/1', producer=identity,
        archive_sha256=ARCHIVE_SHA, source_sha256=hashlib.sha256(members[SOURCE]).hexdigest(),
        commands=commands, compiler_path=compiler_path, compiler_sha256=sha(compiler_path),
        compiler_version=subprocess.check_output([compiler_path, '--version'], text=True).splitlines()[0],
        files={p.name: sha(p) for p in sorted(output.iterdir()) if p.is_file()},
        changes=['No numerical changes to original fitting routines.',
                 'Extract routines after native main program for local C interface.',
                 'Second demo changes only the published iwd parameter from zero to one.'],
        evidence='Numerical geometry component only; no attachment, thermal or hollow-rim evolution solver.')
    (output / 'build.json').write_text(json.dumps(record, indent=2) + '\n')
    return record


def validate_inputs(phi, habit, coordinates):
    if habit not in COLUMNS or coordinates not in ('wood', 'native'):
        raise ValueError('Explicit cylinder/hexagonal_prism and wood/native coordinates required')
    if not math.isfinite(phi) or not .01 <= phi <= 100.:
        raise ValueError('Shape fit requires 0.01 <= aspect ratio <= 100; no clipping')
    if habit == 'hexagonal_prism' and coordinates == 'native' and phi > 100 * NATIVE_COSINE:
        raise ValueError('Converted Wood aspect ratio exceeds 100')


class ShapeLibrary:
    def __init__(self, build_dir=BUILD):
        self.build_dir = Path(build_dir)
        self.identity = json.loads((self.build_dir / 'build.json').read_text())
        if self.identity['producer'] != producer() or self.identity['archive_sha256'] != ARCHIVE_SHA:
            raise ValueError('Shape wrapper/build identity changed')
        for name, digest in self.identity['files'].items():
            if sha(self.build_dir / name) != digest:
                raise ValueError('Shape build artifact changed: ' + name)
        self.lib = ctypes.CDLL(str((self.build_dir / 'libhollow_shape.so').resolve()))
        self.lib.hollow_shape_eval.argtypes = [ctypes.c_double, ctypes.c_int, ctypes.c_int,
                                             ctypes.POINTER(ctypes.c_double)]
        self.lib.hollow_shape_eval.restype = None

    def evaluate(self, phi, habit, coordinates='native'):
        phi = float(phi)
        validate_inputs(phi, habit, coordinates)
        values = (ctypes.c_double * 8)()
        self.lib.hollow_shape_eval(phi, 1 if habit == 'cylinder' else 2,
                                   0 if coordinates == 'wood' else 1, values)
        if any(not math.isfinite(x) or x <= 0 for x in values):
            raise FloatingPointError('Nonpositive or nonfinite native shape function')
        return dict(zip(FIELDS, map(float, values)))


def audit(library):
    members = archive_members()
    rows = []
    for iwd, name in enumerate(TABLES):
        archived = np.loadtxt(io.BytesIO(members[name]))
        native = np.loadtxt(library.build_dir / name)
        if archived.shape != (500, 17) or native.shape != archived.shape:
            raise ValueError('Unexpected archived shape-table dimensions')
        delta = np.abs(native / archived - 1)
        maximum = float(delta.max())
        if maximum > 5e-6:
            raise AssertionError('Native demo differs materially from archived table')
        wrapper_max = 0.; evaluated = 0; rejected = 0
        for row in native:
            for habit, columns in COLUMNS.items():
                # The original float grid's first printed value is very slightly
                # below 0.01; restore only that known single-precision endpoint.
                phi = float(row[0])
                if phi < .01 and abs(phi / .01 - 1) < 1e-7:
                    phi = .01
                try:
                    answer = library.evaluate(phi, habit, 'wood' if iwd == 0 else 'native')
                except ValueError:
                    rejected += 1
                    continue
                wrapper_max = max(wrapper_max, max(abs(answer[f] / row[c] - 1)
                                                   for f, c in zip(FIELDS, columns)))
                evaluated += 1
        if wrapper_max > 5e-6:
            raise AssertionError('Shape wrapper differs from native demonstration')
        rows.append(dict(coordinates='wood' if iwd == 0 else 'native', table=name,
            archived_sha256=hashlib.sha256(members[name]).hexdigest(), native_sha256=sha(library.build_dir / name),
            numerical_entries=8500, maximum_relative_archive_difference=maximum,
            maximum_relative_wrapper_difference=wrapper_max,
            wrapper_evaluations=evaluated, out_of_fit_hexagonal_requests=rejected))
    return dict(comparisons=rows, build=library.identity,
        limits=['Archive comparison establishes numerical reproduction, not measured growth accuracy.',
                'Native conversion uses 0.866 and fits are guarded in the source Wood-coordinate interval.',
                'Capacitance functions do not respond to iwd in the original code; the Wood demonstration is a source audit, not a complete consistent physical state.',
                'A hollow-rim transport/step-source growth model is not provided by these shape functions alone.'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', type=Path)
    args = parser.parse_args()
    if args.build:
        print(json.dumps(build(args.build), indent=2))
    else:
        print(json.dumps(audit(ShapeLibrary()), indent=2))
