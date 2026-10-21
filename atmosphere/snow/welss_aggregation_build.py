"""Build the declared pinned McSnow author solver without changing its physics.

External source and executable bytes live in the user-specified scratch directory.
Every source is checked against a recorded Git blob and SHA256. The sole patch
replaces clock-based RNG seeding with a required integer command-line argument.
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
import time
import urllib.request

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / 'welss_aggregation_sources.json'
FLAGS = ['-O3', '-DNDEBUG', '-cpp']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify(data, item):
    blob = hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
    if len(data) != item['bytes'] or blob != item['git_blob_sha1'] or hashlib.sha256(data).hexdigest() != item['sha256']:
        raise ValueError('source identity mismatch: '+item['path'])


def deterministic_main(original):
    """Keep all physics unchanged; seed-vector entries share the supplied integer."""
    declaration = '  INTEGER :: values(1:8), k\n'
    assign = '  seed(:) = values(8)\n'
    if original.count(declaration) != 1 or original.count(assign) != 1:
        raise ValueError('RNG patch context differs')
    return original.replace(declaration, declaration+'  CHARACTER(len=32) :: deterministic_seed\n').replace(assign,
        '  CALL get_command_argument(1, deterministic_seed)\n'
        '  READ(deterministic_seed, *) values(8)\n'+assign)


def compilation_order(files, compiler, src):
    provided, needed = {}, {}
    for path in files:
        body = subprocess.check_output([compiler, *FLAGS, '-E', '-I'+str(src), str(path)], text=True)
        names = re.findall(r'^\s*module\s+(?!procedure\b|function\b|subroutine\b)(\w+)', body, re.I|re.M)
        for name in names:
            if name.lower() in provided:
                raise ValueError('duplicate module: '+name)
            provided[name.lower()] = path
        needed[path] = {name.lower() for name in re.findall(r'^\s*use\s+(\w+)', body, re.I|re.M)}
    dependencies = {p: {provided[n] for n in names if n in provided and provided[n] != p} for p,names in needed.items()}
    result = []
    while dependencies:
        ready = sorted(p for p, need in dependencies.items() if not need)
        if not ready:
            raise ValueError('cyclic or unparsed Fortran module dependencies')
        result.extend(ready)
        for p in ready:
            dependencies.pop(p)
        for need in dependencies.values():
            need.difference_update(ready)
    return result


def build(source, output, *, fetch=False, compiler='gfortran'):
    source, output = Path(source).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError('build directory already exists: '+str(output))
    identity = {str(p):sha(p) for p in (Path(__file__),MANIFEST)}
    manifest = json.loads(MANIFEST.read_text())
    for item in manifest['files']:
        p = source/item['path']
        if not p.exists():
            if not fetch:
                raise FileNotFoundError(p)
            data = urllib.request.urlopen(item['url'],timeout=60).read()
            verify(data,item)
            p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(data)
        verify(p.read_bytes(),item)
    output.mkdir(parents=True)
    patched = output/'main.f90'
    original = (source/'src/main.f90').read_text()
    changed = deterministic_main(original);patched.write_text(changed)
    patch = ''.join(difflib.unified_diff(original.splitlines(True),changed.splitlines(True),fromfile='pinned/src/main.f90',tofile='deterministic/main.f90'))
    (output/'rng.patch').write_text(patch)
    compiler_path = str(Path(shutil.which(compiler)).resolve())
    version = subprocess.check_output([compiler_path,'--version'],text=True).splitlines()[0]
    files = sorted(p for p in (source/'src').rglob('*') if p.suffix.lower() == '.f90' and p.name != 'main.f90')+[patched]
    order = compilation_order(files,compiler_path,source/'src')
    commands, objects = [], []
    started = time.monotonic()
    with (output/'build.log').open('w') as log:
        for index,p in enumerate(order):
            obj = output/f'{index:02d}_{p.stem}.o';objects.append(str(obj))
            cmd = [compiler_path,*FLAGS,'-I'+str(source/'src'),'-I'+str(output),'-J'+str(output),'-c',str(p),'-o',str(obj)]
            commands.append(cmd);subprocess.run(cmd,cwd=output,stdout=log,stderr=subprocess.STDOUT,check=True)
        cmd=[compiler_path,*FLAGS,*objects,'-o',str(output/'MCSNOW')]
        commands.append(cmd);subprocess.run(cmd,cwd=output,stdout=log,stderr=subprocess.STDOUT,check=True)
    for item in manifest['files']:
        verify((source/item['path']).read_bytes(),item)
    if identity != {str(p):sha(p) for p in (Path(__file__),MANIFEST)}:
        raise RuntimeError('builder or manifest changed during build')
    product=dict(schema='terluna.atmosphere.welss-author-build/1',revision=manifest['revision'],
        producer={p.name:identity[str(p)] for p in (Path(__file__),MANIFEST)},compiler=version,
        compiler_path=compiler_path,compiler_sha256=sha(compiler_path),flags=FLAGS,openmp=False,
        source_dir=str(source),source_files={i['path']:i['sha256'] for i in manifest['files']},
        rng_intervention={'original_sha256':sha(source/'src/main.f90'),'patched_sha256':sha(patched),
                          'unified_diff':patch,'seed_rule':'All GNU Fortran seed-vector entries equal required first integer command-line argument.'},
        commands=commands,executable_sha256=sha(output/'MCSNOW'),build_log_sha256=sha(output/'build.log'),elapsed_s=time.monotonic()-started)
    (output/'build.json').write_text(json.dumps(product,indent=2,allow_nan=False)+'\n')
    return product


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',required=True,type=Path);p.add_argument('--output',required=True,type=Path)
    p.add_argument('--fetch',action='store_true');p.add_argument('--compiler',default='gfortran')
    a=p.parse_args();r=build(a.source,a.output,fetch=a.fetch,compiler=a.compiler)
    print(json.dumps({k:r[k] for k in ['revision','compiler','elapsed_s','executable_sha256']},indent=2))

if __name__=='__main__':
    main()
