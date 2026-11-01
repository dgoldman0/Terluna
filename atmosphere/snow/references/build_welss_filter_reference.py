"""Compile an independent NCL 6.6.2 five-point filter fixture.

Fetch the two URLs in SOURCES into an external directory, then run:
    python build_welss_filter_reference.py --source-dir DIR --output PATH

External NCL Fortran bytes remain outside this repository. This driver imports
no Python filter implementation and executes the pinned official routines.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import struct
import subprocess
import tempfile

SOURCES = {
    'filtrx.f': {
        'sha256': 'f4215266fdc27778f6297fb30b4ec7cf1b2a868c12ef797a182b1e9bbba171bf',
        'url': 'https://raw.githubusercontent.com/NCAR/ncl/6.6.2/ni/src/lib/nfpfort/filtrx.f',
    },
    'wrunave_dp.f': {
        'sha256': 'c21c13343dfe1060d5b25d4df7913fd375c73c8c92d929ad2e70b1285a264e02',
        'url': 'https://raw.githubusercontent.com/NCAR/ncl/6.6.2/ni/src/lib/nfpfort/wrunave_dp.f',
    },
}
MISSING = -9999999.
DRIVER = '''program reference
implicit none
integer, parameter :: nwt=5, npts=13
integer :: ier, ios, lwork
double precision :: wt(nwt), rounded(nwt), resp(2*nwt-1), freq(2*nwt-1)
double precision :: x(npts), work(npts+4)
call DFILTRQ(nwt,dble(.2),-999.d0,1.d0,0,wt,resp,freq,ier)
if (ier /= 0) error stop 'DFILTRQ failed'
rounded=dble(real(wt))
write(*,'(5(es26.17e3,1x))') wt
write(*,'(5(es26.17e3,1x))') rounded
do
  read(*,*,iostat=ios) x
  if (ios /= 0) exit
  lwork=npts+4
  call DWGTRUNAVE(x,npts,rounded,nwt,0,-9999999.d0,work,lwork,ier)
  if (ier /= 0) error stop 'DWGTRUNAVE failed'
  write(*,'(13(es26.17e3,1x))') x
enddo
end program
'''


def float32(value):
    return struct.unpack('f', struct.pack('f', value))[0]


def cases():
    impulse = [0.] * 13
    impulse[6] = 1.
    missing = [7.] * 13
    missing[6] = None
    return [
        {'name': 'constant', 'input': [7.] * 13},
        {'name': 'impulse', 'input': impulse},
        {'name': 'missing_center', 'input': missing},
        {'name': 'float32_ramp', 'input': [float32((i-4)/7) for i in range(13)]},
    ]


def generate(source_dir):
    source_dir = Path(source_dir)
    source_bytes = {}
    for name, record in SOURCES.items():
        data = (source_dir/name).read_bytes()
        if hashlib.sha256(data).hexdigest() != record['sha256']:
            raise ValueError(f'{name}: pinned NCL 6.6.2 source hash mismatch')
        source_bytes[name] = data
    inputs = cases()
    stdin = '\n'.join(' '.join(format(MISSING if x is None else x, '.17g')
                              for x in row['input']) for row in inputs) + '\n'
    flags = ['-O0', '-fno-fast-math', '-ffixed-line-length-none']
    with tempfile.TemporaryDirectory(prefix='welss-ncl-filter-') as directory:
        work = Path(directory)
        for name, data in source_bytes.items():
            (work/name).write_bytes(data)
        (work/'reference.f90').write_text(DRIVER)
        command = ['gfortran', *flags, 'filtrx.f', 'wrunave_dp.f',
                   'reference.f90', '-o', 'reference']
        subprocess.run(command, cwd=work, check=True, capture_output=True, text=True)
        run = subprocess.run([str(work/'reference')], cwd=work, input=stdin,
                             check=True, capture_output=True, text=True)
    rows = [[float(x) for x in line.split()] for line in run.stdout.splitlines()]
    if len(rows) != len(inputs)+2 or any(len(row) != 13 for row in rows[2:]):
        raise RuntimeError('unexpected NCL reference output shape')
    for row, output in zip(inputs, rows[2:]):
        row['expected'] = [None if x == MISSING else x for x in output]
    return {
        'schema': 'terluna.atmosphere.welss-filter-reference/1',
        'evidence': 'Independently compiled source parity fixture using synthetic arrays; not an archived population run or empirical validation.',
        'source_version': 'NCL 6.6.2',
        'source_credit': 'NCAR Command Language (NCL), UCAR/NCAR. Official source files are retrieved separately and are not redistributed here.',
        'external_source_files': SOURCES,
        'generation_script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'driver_sha256': hashlib.sha256(DRIVER.encode()).hexdigest(),
        'compiler': subprocess.run(['gfortran', '--version'], check=True,
                                   capture_output=True, text=True).stdout.splitlines()[0],
        'compiler_flags': flags,
        'settings': {'nwt': 5, 'fca': float32(.2), 'fcb': -999., 'nsigma': 1.,
                     'ihp': 0, 'kopt': 0, 'npts': 13, 'missing_sentinel': MISSING},
        'precision': 'DFILTRQ receives float32(.2) promoted to double; resulting double weights are rounded to default REAL (binary32) then promoted for DWGTRUNAVE. Synthetic inputs and arithmetic in DWGTRUNAVE are double; float32_ramp inputs are rounded explicitly. Output is serialized from double; missing values become null.',
        'scope': 'One official weight-generation and running-average pass, including conditional normalization, missing-input propagation, and two masked endpoints on each side. The archive producer and NCL runtime version remain unverified.',
        'weights_double': rows[0],
        'weights_float32_promoted': rows[1],
        'cases': inputs,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(generate(args.source_dir), indent=2,
                                     allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
