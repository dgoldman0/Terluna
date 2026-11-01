"""Pinned local interface to Harrington et al. (2021) facet coefficients.

The original Fortran stays in the ignored research cache. A build removes two
unused assignments that read uninitialized variables, splits off the published
demo, and links the routines to our C interface. An explicit optional variant
uses equivalent stable spheroid-capacitance arithmetic. This is a
component audit, not a validated crystal-growth or ventilation implementation.
"""
from __future__ import annotations

import argparse
import ctypes
import difflib
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DRIVER = HERE / 'facet_attachment_driver.f90'
SOURCE = ROOT / 'research/runs/snow_scale/growth_physics_sources/harrington2021_code/depcoeffparamsubs.f90'
SOURCE_SHA256 = 'f6a05f4ca899a2e08851b8cc3c65f449ab819104c9f79b868be75fa779637ef1'
SOURCE_DOI = 'https://doi.org/10.26208/s7de-et44'
SOURCE_URL = 'https://www.datacommons.psu.edu/download/meteorology/harrington/Harrington_CodeCalculateDepositionCoefficients2021.zip'
ARCHIVE_SHA256 = '34cbab3aea0bde80794c4edd9c994d1d406484b265f047cb486e93c86d7e8259'
FLAGS = ['-O0', '-std=legacy', '-fno-fast-math', '-fcheck=all', '-finit-real=snan',
         '-ffpe-trap=invalid,zero,overflow', '-ffree-line-length-none']
FIELDS = ('alpha_a', 'alpha_c', 'alpha_t', 'capacitance_m', 'boundary_capacitance_m',
          'diffusivity_m2_s', 'conductivity_W_m_K', 'ice_saturation_pressure_Pa',
          'latent_sublimation_J_kg', 'vapor_molecular_speed_m_s', 'dry_air_density_kg_m3',
          'vapor_mean_free_path_m', 'boundary_offset_m', 'thermal_jump_length_m',
          'characteristic_supersaturation_a', 'characteristic_supersaturation_c',
          'vapor_gas_constant_J_kg_K', 'dry_air_gas_constant_J_kg_K',
          'dry_air_heat_capacity_J_kg_K', 'decimal_precision')
INPUT_FIELDS = ('ice_supersaturation', 'temperature_K', 'pressure_Pa',
                'equatorial_radius_m', 'polar_radius_m', 'm')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def producer():
    return {Path(__file__).name: sha(__file__), DRIVER.name: sha(DRIVER)}


def prepare_source(original, *, stable_capacitance=False):
    if hashlib.sha256(original.encode()).hexdigest() != SOURCE_SHA256:
        raise ValueError('Unrecognized Harrington source bytes')
    patched = original
    for axis in ('a', 'c'):
        line = f'      scratio_{axis} = min(max(slocal_diff_{axis}/scrit_{axis},scrat_low),scrat_hi)'
        # The c-axis statement has three additional indentation spaces.
        matches = [row for row in patched.splitlines(keepends=True) if row.strip() == line.strip()]
        if len(matches) != 1:
            raise ValueError('Unexpected dead-assignment patch anchor')
        patched = patched.replace(matches[0], '', 1)
    if stable_capacitance:
        begin = patched.index('function capacitance(a,c)')
        end = patched.index('end function capacitance', begin) + len('end function capacitance')
        patched = patched[:begin] + '''function capacitance(a,c)
  ! Numerically stable, analytically equivalent spheroid capacitance.
  ! Differences of squares are factored; numerator and inverse-function
  ! argument use the same dimensionless eccentricity variable.
  implicit none
  real capacitance,a,c,x
  if (c<a) then
    x=sqrt((a-c)*(a+c))/a
    capacitance=a*x/asin(x)
  elseif (c>a) then
    x=sqrt((c-a)*(c+a))/a
    capacitance=a*x/asinh(x)
  else
    capacitance=a
  endif
end function capacitance''' + patched[end:]
    marker = 'end program depcoefftest'
    if patched.count(marker) != 1:
        raise ValueError('Unexpected demonstration boundary')
    routines = patched.split(marker, 1)[1]
    return patched, routines


def build(output, *, source=SOURCE, precision='single', compiler='gfortran', stable_capacitance=False):
    """Build local library and original demo; never download or redistribute code."""
    output = Path(output).resolve()
    source = Path(source).resolve()
    if output.exists():
        raise FileExistsError(output)
    if precision not in ('single', 'double'):
        raise ValueError('precision must be single or double')
    if sha(source) != SOURCE_SHA256:
        raise ValueError('Source hash mismatch')
    original = source.read_text()
    patched, routines = prepare_source(original, stable_capacitance=stable_capacitance)
    identity = producer()
    compiler_path = shutil.which(compiler)
    if compiler_path is None:
        raise FileNotFoundError(compiler)
    compiler_path = str(Path(compiler_path).resolve())
    flags = [*FLAGS]
    if precision == 'double':
        flags += ['-fdefault-real-8', '-fdefault-double-8']
    output.mkdir(parents=True)
    (output / 'native_demo.f90').write_text(patched)
    (output / 'native_routines.f90').write_text(routines)
    (output / 'source.patch').write_text(''.join(difflib.unified_diff(
        original.splitlines(keepends=True), patched.splitlines(keepends=True),
        fromfile='depcoeffparamsubs.f90.original', tofile='depcoeffparamsubs.f90.audited')))
    shutil.copyfile(DRIVER, output / DRIVER.name)
    commands = [
        [compiler_path, *flags, 'native_demo.f90', '-o', 'native_demo'],
        [compiler_path, *flags, '-fPIC', '-shared', 'native_routines.f90', DRIVER.name,
         '-o', 'libfacet_attachment.so'],
    ]
    with (output / 'build.log').open('w') as log:
        for command in commands:
            subprocess.run(command, cwd=output, stdout=log, stderr=subprocess.STDOUT, check=True)
    with (output / 'native_demo.stdout').open('w') as log:
        subprocess.run([str(output / 'native_demo')], cwd=output, stdout=log,
                       stderr=subprocess.STDOUT, check=True)
    if identity != producer() or sha(source) != SOURCE_SHA256:
        raise RuntimeError('Producer or source changed during build')
    files = {name: sha(output / name) for name in (
        'native_demo.f90', 'native_routines.f90', 'source.patch', DRIVER.name,
        'native_demo', 'libfacet_attachment.so', 'build.log',
        'native_demo.stdout', 'depcoeffparamtest.dat')}
    record = dict(schema='terluna.atmosphere.facet-attachment-build/1', producer=identity,
        source_path=str(source), source_sha256=SOURCE_SHA256, source_doi=SOURCE_DOI,
        source_url=SOURCE_URL, archive_sha256=ARCHIVE_SHA256, precision=precision,
        stable_capacitance=stable_capacitance,
        flags=flags, commands=commands, compiler_path=compiler_path,
        compiler_sha256=sha(compiler_path),
        compiler_version=subprocess.check_output([compiler_path, '--version'], text=True).splitlines()[0],
        files=files,
        patches=['Remove exactly two unused scratio_a/c assignments reading uninitialized scrat_low/hi.',
                 'Extract routines after the published main program for the shared library.'] +
                (['Replace only spheroid capacitance arithmetic with equivalent factored eccentricity/asin and asinh forms.']
                 if stable_capacitance else []),
        wrapper='C-double interface converts inputs to native REAL; original coefficient operations unchanged.',
        source_limits=['Positive ice supersaturation; no sublimation.',
                       'Characteristic-supersaturation fits guarded to -70 through -1 C; m in [1,15].',
                       'Native coefficients have an existing floor of 1e-6.',
                       'No ventilation input; assumes source air transport properties.',
                       'Source metadata requires annual notification of customization; no notification sent.'],
        evidence='Compiled published component; native demo is a numerical reference, not empirical validation.')
    (output / 'build.json').write_text(json.dumps(record, indent=2, allow_nan=False) + '\n')
    return record


def verify_build(build_dir):
    build_dir = Path(build_dir).resolve()
    record = json.loads((build_dir / 'build.json').read_text())
    for name, digest in record['files'].items():
        if sha(build_dir / name) != digest:
            raise ValueError('Build artifact changed: ' + name)
    if record['source_sha256'] != SOURCE_SHA256 or record['producer'] != producer():
        raise ValueError('Source or wrapper identity differs')
    return record


class AttachmentLibrary:
    """Fast, serial ctypes calls; native fitting routines are not thread audited."""

    def __init__(self, build_dir):
        self.build_dir = Path(build_dir).resolve()
        self.identity = verify_build(self.build_dir)
        self.library = ctypes.CDLL(str(self.build_dir / 'libfacet_attachment.so'))
        self.library.facet_attachment_eval.argtypes = [ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int)]
        self.library.facet_attachment_eval.restype = None
        self.library.facet_attachment_demo_state.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_double)]
        self.library.facet_attachment_demo_state.restype = None

    def evaluate(self, ice_supersaturation, temperature_K, pressure_Pa,
                 equatorial_radius_m, polar_radius_m, m=10.):
        values = tuple(map(float, (ice_supersaturation, temperature_K, pressure_Pa,
                                  equatorial_radius_m, polar_radius_m, m)))
        if not all(math.isfinite(x) for x in values):
            raise ValueError('All facet inputs must be finite')
        if values[0] <= 0 or not 273.15 - 70. <= values[1] <= 273.15 - 1.:
            raise ValueError('Requires positive ice supersaturation and -70 <= T_C <= -1')
        if min(values[2:5]) <= 0 or not 1 <= values[5] <= 15:
            raise ValueError('Requires positive pressure and semiaxes, and 1 <= m <= 15')
        inputs = (ctypes.c_double * 6)(*values)
        outputs = (ctypes.c_double * len(FIELDS))()
        status = ctypes.c_int()
        self.library.facet_attachment_eval(inputs, outputs, ctypes.byref(status))
        if status.value:
            raise FloatingPointError(f'Native facet evaluation status {status.value}; inputs={values}')
        return dict(zip(FIELDS, map(float, outputs)))

    def demo_state(self, index):
        if not 10 <= index <= 1000:
            raise ValueError('Native demonstration indices are 10 through 1000')
        inputs = (ctypes.c_double * 6)()
        self.library.facet_attachment_demo_state(index, inputs)
        return dict(zip(INPUT_FIELDS, map(float, inputs)))


def audit(single_dir, double_dir):
    """Compare both wrappers with their native demos, then quantify precision drift."""
    libraries = [AttachmentLibrary(single_dir), AttachmentLibrary(double_dir)]
    runs = []
    all_outputs = []
    for library in libraries:
        if library.identity['precision'] != ('single' if not runs else 'double'):
            raise ValueError('Audit requires single then double builds')
        rows = []
        for line in (library.build_dir / 'native_demo.stdout').read_text().splitlines():
            tokens = line.split()
            if tokens and tokens[0] == 'test':
                rows.append(list(map(float, tokens[1:])))
        if len(rows) != 991:
            raise ValueError('Native demonstration did not produce 991 rows')
        maximum_relative_error = 0.
        outputs = []
        for index, row in zip(range(10, 1001), rows):
            answer = library.evaluate(**library.demo_state(index))
            outputs.append(answer)
            for field, expected in zip(('alpha_a', 'alpha_c'), row[2:]):
                maximum_relative_error = max(maximum_relative_error, abs(answer[field]/expected - 1.))
        tolerance = 6e-7 if library.identity['precision'] == 'single' else 6e-15
        if maximum_relative_error > tolerance:
            raise AssertionError('Wrapper differs from native demonstration')
        runs.append(dict(precision=library.identity['precision'], rows=len(rows),
                         maximum_relative_difference_from_printed_native=maximum_relative_error,
                         tolerance=tolerance, build_sha256=sha(library.build_dir / 'build.json')))
        all_outputs.append(outputs)
    differences = {field: max(abs(single[field]/double[field]-1.)
        for single, double in zip(*all_outputs)) for field in FIELDS if field != 'decimal_precision'}
    guards = []
    base = dict(ice_supersaturation=.05, temperature_K=258.15, pressure_Pa=10000.,
                equatorial_radius_m=1e-4, polar_radius_m=1e-4, m=10.)
    for field, bad in [('ice_supersaturation',0.), ('temperature_K',273.),
                       ('temperature_K',200.), ('pressure_Pa',0.),
                       ('equatorial_radius_m',0.), ('polar_radius_m',float('nan')), ('m',16.)]:
        try:
            libraries[0].evaluate(**{**base, field: bad})
        except ValueError:
            guards.append(field)
        else:
            raise AssertionError('Invalid input accepted: ' + field)
    return dict(schema='terluna.atmosphere.facet-attachment-audit/1', native_demo_comparisons=runs,
        maximum_relative_single_double_differences=differences, invalid_input_guards_checked=guards,
        evidence='Numerical component parity and precision sensitivity; no empirical growth validation.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    command = sub.add_parser('build')
    command.add_argument('output', type=Path)
    command.add_argument('--source', type=Path, default=SOURCE)
    command.add_argument('--precision', choices=('single', 'double'), default='single')
    command.add_argument('--stable-capacitance', action='store_true')
    command = sub.add_parser('audit')
    command.add_argument('single', type=Path)
    command.add_argument('double', type=Path)
    command.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'build':
        result = build(args.output, source=args.source, precision=args.precision,
                       stable_capacitance=args.stable_capacitance)
        print(json.dumps({'precision': result['precision'], 'files': result['files']}, indent=2))
    else:
        result = audit(args.single, args.double)
        args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
        print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
