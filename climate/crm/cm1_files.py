"""CM1's 3-D outputs and restart files kept raw or compressed with zstd, read bit for bit the same either way.

    climate/gcm/.venv/bin/python -m climate.crm.cm1_files compress box_0e      # every 3-D output and restart, one at a time
    climate/gcm/.venv/bin/python -m climate.crm.cm1_files verify box_0e        # every compressed file against its hash
    climate/gcm/.venv/bin/python -m climate.crm.cm1_files decompress box_0e cm1rst_t000085_s.dat

A case's 3-D outputs (cm1out_tNNNNNN_s.dat, 20-80 MB each) and restart sets (cm1rst_tNNNNNN_{i,s,u,v,w,x}.dat) hold
nearly everything CM1 writes. zstd at level 3 shrinks a sea box's output about 1.9 times and its largest restart file
about 3.6 times, losslessly. compress takes a case's files one at a time through the system's zstd, streaming each, so
no file is held whole in memory or copied anywhere else:
- it stops (SpaceFloor) when the drive's free space is below FLOOR_BYTES plus the file's size;
- it reads the file in CHUNK pieces into a SHA-256 and into zstd -3 -T1, which writes <file>.zst.part beside it;
- it decompresses the .part through zstd -dc into a second SHA-256, after fsyncing it and asking the kernel to drop it
  from the page cache, so the check reads it back from the drive; if the two hashes differ it removes the .part, keeps
  the raw file and goes on;
- it fsyncs the .part, renames it <file>.zst with the raw file's times, appends the file's line (name, raw size, raw
  SHA-256, compressed size, time) to the case's compression_manifest.jsonl and fsyncs that, and only then removes the
  raw file.
A compress taken up again removes any leftover .part and redoes that file. Where a raw file and its .zst both stand,
the raw file goes only when the manifest holds the file's line and the raw file still matches it; with no line the
.zst goes and the file is redone; a line the raw file no longer matches leaves both for a person to look at. A file
whose raw is gone and whose line stands is done. decompress writes raw files back by the same streamed, hash-checked
steps and keeps the .zst; verify checks every compressed file against its line. A case holding run.lock, or one a CM1
process is working in, is left alone.

The readers here take the raw file's path and read <path>.zst when the raw file is absent: fromfile is np.fromfile,
bit for bit; memmap stands in for np.memmap, read only (a compressed file is decompressed whole, so a raw file's memory
map stays the lighter reader); sha256 and size are the raw file's. The listings (glob, restart_numbers) find a file
kept either way. Everything else in a case folder (statistics, diagnostics, metadata, control files, the other output
kinds, namelists, logs) stays raw and is read as before.
"""
from __future__ import annotations
import argparse
import fnmatch
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
RUNS = HERE / 'runs'
ZSTD = '/usr/bin/zstd'                # the system's zstd; the climate venv carries no zstandard module
LEVEL = 3                             # level 19 gains about 4% more at many times the time
CHUNK = 4 << 20                       # bytes read, hashed and piped at a time
FLOOR_BYTES = 30 * 10**9              # free space kept on the drive beyond the file being written
MANIFEST = 'compression_manifest.jsonl'
SUFFIX = '.zst'
PART = '.part'                        # a file being written: <file>.zst.part compressing, <file>.part restoring
KEPT = re.compile(r'(cm1out_t\d{6}_s|cm1rst_t\d{6}_[isuvwx])\.dat')    # what compress takes; everything else stays raw
RESTART_PARTS = 'isuvwx'
PROCESS_CHECK_S = 30.0                # how often compress looks again for CM1 processes in its case


class SpaceFloor(RuntimeError):
    """The drive's free space is below the floor plus the file about to be written."""


class InUse(RuntimeError):
    """A case holds run.lock or a CM1 process is working in it."""


def compressed(path) -> Path:
    """The compressed file kept for path: path.zst."""
    path = Path(path)
    return path.with_name(path.name + SUFFIX)


def exists(path) -> bool:
    """Whether path stands raw or compressed."""
    return Path(path).exists() or compressed(path).exists()


def _decompressed(zpath, chunk: int | None = None):
    """zstd -dc's output for a compressed file, in pieces (CHUNK bytes unless chunk is given); raises when zstd fails,
    as it does on a damaged or cut file (its frames carry their own checksum)."""
    chunk = chunk or CHUNK
    process = subprocess.Popen([ZSTD, '-d', '-c', '-q', '--', str(zpath)], stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    finished = False
    try:
        while block := process.stdout.read(chunk):
            yield block
        finished = True
    finally:
        if not finished:                                               # the reader stopped early: end this zstd
            process.kill()
        process.stdout.close()
        error = process.stderr.read().decode(errors='replace').strip()
        process.stderr.close()
        code = process.wait()
    if code != 0:
        raise RuntimeError(f'zstd could not decompress {zpath}: {error or f"exit {code}"}')


def _chunks(path, chunk: int | None = None):
    """The raw bytes of path in pieces: read from the file, or decompressed from path.zst when the file is absent."""
    path, chunk = Path(path), chunk or CHUNK
    try:
        stream = path.open('rb')
    except FileNotFoundError:
        if not compressed(path).exists():
            raise
        yield from _decompressed(compressed(path), chunk)
        return
    with stream:
        while block := stream.read(chunk):
            yield block


def _whole(zpath) -> bytearray:
    """A compressed file decompressed into one writable buffer, as a reader that needs the whole file takes it."""
    out = bytearray()
    for block in _decompressed(zpath):
        out += block
    return out


def fromfile(path, dtype=float, count: int = -1, offset: int = 0):
    """np.fromfile(path, dtype, count, offset=offset) of a file kept raw or compressed: path itself, or path.zst
    decompressed when path is absent, bit for bit the same array (writable, as np.fromfile's is)."""
    import numpy as np
    path = Path(path)
    try:
        return np.fromfile(path, dtype=dtype, count=count, offset=offset)
    except FileNotFoundError:
        if not compressed(path).exists():
            raise
    data = _whole(compressed(path))
    itemsize = np.dtype(dtype).itemsize
    available = (len(data) - offset) // itemsize                        # whole items after offset, as np.fromfile reads
    if available < 0:
        raise ValueError('negative dimensions are not allowed')
    items = available if count < 0 else min(count, available)
    if items == 0:
        return np.empty(0, dtype=dtype)
    return np.frombuffer(data, dtype=dtype, count=items, offset=offset)


def memmap(path, dtype='uint8', mode: str = 'r', offset: int = 0, shape=None):
    """np.memmap(path, dtype, mode='r', offset=offset, shape=shape) of a file kept raw or compressed, for reading: the
    memory map itself when the raw file stands, or else the decompressed file as a read-only array of the same values
    and shape."""
    import numpy as np
    if mode != 'r':
        raise ValueError(f"mode {mode!r}: cm1_files.memmap only reads (mode 'r')")
    path = Path(path)
    try:
        return np.memmap(path, dtype=dtype, mode='r', offset=offset, shape=shape)
    except FileNotFoundError:
        if not compressed(path).exists():
            raise
    data = _whole(compressed(path))
    itemsize = np.dtype(dtype).itemsize
    if shape is None:
        if (len(data) - offset) % itemsize:
            raise ValueError('Size of available data is not a multiple of the data-type size.')
        shape = ((len(data) - offset) // itemsize,)
    shape = (shape,) if isinstance(shape, int) else tuple(shape)
    items = 1
    for n in shape:
        items *= n
    if offset + items * itemsize > len(data):
        raise ValueError('mmap length is greater than file size')
    out = np.frombuffer(data, dtype=dtype, count=items, offset=offset).reshape(shape)
    out.flags.writeable = False
    return out


def sha256(path) -> str:
    """The SHA-256 (hex) of a file's raw bytes, read from the file or from the decompression of path.zst."""
    digest = hashlib.sha256()
    for block in _chunks(path):
        digest.update(block)
    return digest.hexdigest()


_MANIFESTS = {}


def manifest(folder) -> dict:
    """A case's compression_manifest.jsonl: the latest whole line for each file name, a line a crash cut short passed
    over. Kept parsed until the file changes, for readers that ask once a snapshot."""
    path = Path(folder) / MANIFEST
    try:
        st = path.stat()
    except FileNotFoundError:
        return {}
    key = (st.st_size, st.st_mtime_ns, st.st_ino)
    held = _MANIFESTS.get(path)
    if held is not None and held[0] == key:
        return held[1]
    try:
        text = path.read_text()
    except FileNotFoundError:
        return {}
    lines = {}
    for line in text.splitlines():
        try:
            entry = json.loads(line)
        except json.JSONDecodeError:
            continue
        if isinstance(entry, dict) and {'name', 'raw_bytes', 'raw_sha256', 'zst_bytes'} <= entry.keys():
            lines[entry['name']] = entry
    _MANIFESTS[path] = (key, lines)
    return lines


def size(path) -> int:
    """A file's raw size in bytes: the raw file's, or for a compressed one its manifest line's (counted from the
    decompression when the line is missing)."""
    path = Path(path)
    try:
        return path.stat().st_size
    except FileNotFoundError:
        if not compressed(path).exists():
            raise
    line = manifest(path.parent).get(path.name)
    return line['raw_bytes'] if line else sum(len(block) for block in _decompressed(compressed(path)))


def glob(folder, pattern: str) -> list:
    """folder.glob(pattern) over files kept raw or compressed: the raw path of every file matching pattern whose raw
    file or .zst stands, sorted; files being written (.part) are left out."""
    folder = Path(folder)
    match = re.compile(fnmatch.translate(pattern)).match
    names = set()
    try:
        entries = os.scandir(folder)
    except (FileNotFoundError, NotADirectoryError):
        return []
    with entries:
        for entry in entries:
            name = entry.name
            if name.endswith(PART):
                continue
            if name.endswith(SUFFIX):
                name = name[:-len(SUFFIX)]
            if match(name):
                names.add(name)
    return [folder / name for name in sorted(names)]


def restart_files(case, n: int) -> list:
    """The raw paths of restart set n's six files."""
    return [Path(case) / f'cm1rst_t{n:06d}_{part}.dat' for part in RESTART_PARTS]


def restart_numbers(case) -> list:
    """The numbers of a case's complete restart sets, every one of the six files raw or compressed, in order."""
    names = {p.name for p in glob(case, 'cm1rst_t*.dat')}
    numbers = {int(m.group(1)) for name in names if (m := re.fullmatch(r'cm1rst_t(\d+)_s\.dat', name))}
    return sorted(n for n in numbers if all(p.name in names for p in restart_files(case, n)))


def remove(path) -> list:
    """Remove a file kept raw or compressed, both if both stand; returns the names removed. Its manifest line stays,
    as the record of what was compressed."""
    removed = []
    for p in (Path(path), compressed(path)):
        try:
            p.unlink()
            removed.append(p.name)
        except FileNotFoundError:
            pass
    return removed


def written_bytes(case) -> int:
    """The bytes a case's files hold raw, as CM1 wrote them: each raw file's size and each compressed file's raw size
    (once where both stand); the manifest, files being written and links left out."""
    case = Path(case)
    total, held = 0, manifest(case)
    with os.scandir(case) as entries:
        files = {e.name: e for e in entries if e.is_file(follow_symlinks=False)}
    for name, entry in files.items():
        if name == MANIFEST or name.endswith(PART):
            continue
        if name.endswith(SUFFIX) and KEPT.fullmatch(name[:-len(SUFFIX)]):
            raw = name[:-len(SUFFIX)]
            if raw not in files:
                total += held[raw]['raw_bytes'] if raw in held else size(case / raw)
            continue
        total += entry.stat(follow_symlinks=False).st_size
    return total


def free_bytes(folder) -> int:
    """Free space on the drive holding folder, in bytes, as an ordinary user may fill it (df's Avail)."""
    st = os.statvfs(folder)
    return st.f_bavail * st.f_frsize


def _check_space(folder, needed: int, floor_bytes: float):
    free = free_bytes(folder)
    if free < floor_bytes + needed:
        raise SpaceFloor(f'{free / 1e9:.1f} GB free on the drive of {folder}, below the floor of {floor_bytes / 1e9:.0f} '
                         f'GB plus {needed / 1e9:.2f} GB')


def _drop_cache(fd: int):
    """Let the kernel drop a file's pages from its page cache, so streaming hundreds of gigabytes does not push other
    work's files out of memory (and so a check after fsync reads the drive)."""
    try:
        os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
    except (AttributeError, OSError):
        pass


def _fsync_folder(folder):
    fd = os.open(folder, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _append_line(folder, entry: dict):
    """Append one line to a case's manifest and fsync it, starting a new line if a crash left the last one cut."""
    path = Path(folder) / MANIFEST
    new = not path.exists()
    try:
        with path.open('rb') as stream:
            stream.seek(-1, os.SEEK_END)
            whole = stream.read(1) == b'\n'
    except (FileNotFoundError, OSError):                                # absent, or empty
        whole = True
    with path.open('ab') as out:
        out.write((b'' if whole else b'\n') + json.dumps(entry).encode() + b'\n')
        out.flush()
        os.fsync(out.fileno())
    if new:
        _fsync_folder(folder)


def _remove_raw(path: Path):
    """The last step of compressing a file, once its .zst and manifest line are on the drive."""
    path.unlink()


_VERSION = []


def zstd_version() -> str:
    if not _VERSION:
        text = subprocess.run([ZSTD, '--version'], capture_output=True, text=True, check=True).stdout
        m = re.search(r'v(\d+\.\d+\.\d+)', text)
        _VERSION.append(m.group(1) if m else text.strip())
    return _VERSION[0]


def _hash_decompressed(zpath) -> tuple:
    """The SHA-256 (hex) and byte count of a compressed file's decompression, streamed."""
    digest, count = hashlib.sha256(), 0
    for block in _decompressed(zpath):
        digest.update(block)
        count += len(block)
    return digest.hexdigest(), count


def _differs(part: Path, digest: str, count: int):
    """Why a .part does not decompress to the raw file's bytes, or None when it does."""
    try:
        check, n = _hash_decompressed(part)
    except RuntimeError as e:
        return str(e)
    if (check, n) != (digest, count):
        return f'the .part decompresses to {n} bytes ({check[:16]}), the raw file holds {count} ({digest[:16]})'
    return None


def compress_file(path, floor_bytes: float = FLOOR_BYTES) -> dict:
    """One file through compress's steps (the module's docstring). Returns its manifest line with status 'compressed',
    or status 'failed' and the reason, with the raw file kept and no .part left; raises SpaceFloor before touching
    anything when the drive is too full."""
    path = Path(path)
    folder, part, target = path.parent, path.with_name(path.name + SUFFIX + PART), compressed(path)
    if target.exists() or part.exists():
        raise RuntimeError(f'{target.name} or its .part already stands; compress() settles those first')
    before = path.stat()
    _check_space(folder, before.st_size, floor_bytes)
    raw_digest, read = hashlib.sha256(), 0                              # one pass: hash and compress together
    with part.open('xb') as out:
        process = subprocess.Popen([ZSTD, f'-{LEVEL}', '-T1', '-q', '-c'], stdin=subprocess.PIPE, stdout=out,
                                   stderr=subprocess.PIPE)
        try:
            with path.open('rb') as raw:
                while block := raw.read(CHUNK):
                    raw_digest.update(block)
                    read += len(block)
                    process.stdin.write(block)
                _drop_cache(raw.fileno())
        except BrokenPipeError:                                         # zstd stopped; its exit code says why
            pass
        finally:                                                        # zstd ends on its input's end, whatever happened
            try:
                process.stdin.close()
            except BrokenPipeError:
                pass
            error = process.stderr.read().decode(errors='replace').strip()
            process.stderr.close()
            code = process.wait()
        out.flush()
        os.fsync(out.fileno())                                          # on the drive, and out of the page cache,
        _drop_cache(out.fileno())                                       # before the check reads it back
    after = path.stat()
    if code != 0:
        reason = f'zstd stopped with exit {code}: {error}'
    elif read != before.st_size or (after.st_size, after.st_mtime_ns) != (before.st_size, before.st_mtime_ns):
        reason = 'the raw file changed while it was read'
    else:
        reason = _differs(part, raw_digest.hexdigest(), read)          # the check, streamed
    if reason is not None:
        part.unlink()
        return dict(name=path.name, status='failed', reason=reason)
    fd = os.open(part, os.O_RDONLY)                                     # commit: fsync, rename, line, then removal
    try:
        os.fsync(fd)
        _drop_cache(fd)
    finally:
        os.close(fd)
    os.utime(part, ns=(before.st_atime_ns, before.st_mtime_ns))
    os.rename(part, target)
    _fsync_folder(folder)
    entry = dict(name=path.name, raw_bytes=read, raw_sha256=raw_digest.hexdigest(), zst_bytes=target.stat().st_size,
                 time=time.strftime('%Y-%m-%d %H:%M:%S'), raw_mtime_ns=before.st_mtime_ns,
                 zstd=f'{zstd_version()} -{LEVEL} -T1')
    _append_line(folder, entry)
    _remove_raw(path)
    return dict(entry, status='compressed')


def cm1_processes() -> list:
    """The process ids whose command line runs a CM1 executable (cm1.exe or a build's cm1_<hash>.exe), from ps."""
    text = subprocess.run(['ps', '-eo', 'pid=,args='], capture_output=True, text=True, check=True).stdout
    own = os.getpid()
    out = []
    for line in text.splitlines():
        pid, _, args = line.strip().partition(' ')
        if re.search(r'(^|[/\s])cm1[^/\s]*\.exe(\s|$)', args) and pid.isdigit() and int(pid) != own:
            out.append(int(pid))
    return out


def in_use(case):
    """Why a case may not be compressed now: its runner's run.lock, or a CM1 process working in it (or one whose
    folder cannot be read); None when neither."""
    case = Path(case)
    if (case / 'run.lock').exists():
        return 'it holds run.lock'
    here = case.resolve()
    for pid in cm1_processes():
        try:
            cwd = Path(os.readlink(f'/proc/{pid}/cwd')).resolve()
        except OSError:
            if Path(f'/proc/{pid}').exists():
                return f'CM1 process {pid} runs in a folder that cannot be read'
            continue                                                    # ended since ps listed it
        if cwd == here or here in cwd.parents:
            return f'CM1 process {pid} runs in it'
    return None


def _exclusive(case) -> int:
    """Lock the case's folder (flock), so two compressions of one case cannot run at once; returns the descriptor to
    close."""
    import fcntl
    fd = os.open(case, os.O_RDONLY)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        os.close(fd)
        raise InUse(f'{Path(case).name}: another compression holds it')
    return fd


def kept_names(case) -> list:
    """The names (raw) of the 3-D outputs and restart files a case holds raw or compressed, the files compress takes."""
    return [p.name for p in glob(case, 'cm1*.dat') if KEPT.fullmatch(p.name)]


def disk_bytes(case, names) -> tuple:
    """For the named files of a case: their raw bytes, and the bytes they take on the drive as kept (raw, or .zst)."""
    case, held = Path(case), manifest(case)
    raw = disk = 0
    for name in names:
        path = case / name
        if path.exists():
            n = path.stat().st_size
            raw, disk = raw + n, disk + n
        elif compressed(path).exists():
            raw += held[name]['raw_bytes'] if name in held else size(path)
            disk += compressed(path).stat().st_size
    return raw, disk


def compress(case, floor_bytes: float = FLOOR_BYTES, report=None) -> dict:
    """Compress every 3-D output and restart file of a case, one at a time (compress_file), taking up an interrupted
    compression where it stopped (the module's docstring); report, when given, takes a line for each file. Returns
    what it did: the files compressed, finished (an earlier commit's raw file removed), done before, failed (raw kept)
    and left alone, with the files' raw bytes and their bytes on the drive before and after. Raises InUse for a case
    in use and SpaceFloor when the drive is too full to go on, either with the counts so far as .result."""
    case = Path(case)
    say = report or (lambda line: None)
    reason = in_use(case)
    if reason:
        raise InUse(f'{case.name}: {reason}')
    lock = _exclusive(case)
    result = dict(case=case.name, compressed=[], finished=[], done=[], failed=[], left=[], parts_removed=[])
    try:
        for part in sorted(p for p in case.iterdir() if p.name.endswith(PART)):    # leftovers of a stopped run
            stem = part.name[:-len(PART)]
            if KEPT.fullmatch(stem[:-len(SUFFIX)] if stem.endswith(SUFFIX) else stem):
                part.unlink()
                result['parts_removed'].append(part.name)
                say(f'{part.name}: left by an earlier run, removed')
        names = kept_names(case)
        result['raw_bytes'], result['disk_bytes_before'] = disk_bytes(case, names)
        lines = manifest(case)
        looked = time.monotonic()
        for name in names:
            raw, target = case / name, compressed(case / name)
            line = lines.get(name)
            if raw.is_symlink() or target.is_symlink():
                result['left'].append(dict(name=name, reason='a link, not a file of this case'))
                continue
            if raw.exists() and target.exists():                       # stopped between the rename and the removal
                if line is None:
                    target.unlink()
                    say(f'{name}: a .zst without its manifest line, removed to redo')
                elif (line['raw_bytes'] == raw.stat().st_size and line['zst_bytes'] == target.stat().st_size
                      and sha256(raw) == line['raw_sha256']):
                    _remove_raw(raw)
                    result['finished'].append(name)
                    say(f'{name}: its .zst and line stood, raw file removed')
                    continue
                else:
                    result['left'].append(dict(name=name, reason='the raw file and .zst both stand and the raw file '
                                                                 'differs from its manifest line'))
                    say(f'{name}: left as it is, the raw file differs from its manifest line')
                    continue
            if not raw.exists():
                if line is not None:
                    result['done'].append(name)
                else:
                    result['left'].append(dict(name=name, reason='a .zst with no raw file and no manifest line'))
                    say(f'{name}: left as it is, a .zst with no manifest line')
                continue
            if time.monotonic() - looked > PROCESS_CHECK_S:             # a CM1 process now, or run.lock
                reason, looked = in_use(case), time.monotonic()
            else:
                reason = 'it holds run.lock' if (case / 'run.lock').exists() else None
            if reason:
                raise InUse(f'{case.name}: {reason} (stopped before {name})')
            entry = compress_file(raw, floor_bytes)
            if entry['status'] == 'compressed':
                result['compressed'].append(name)
                say(f"{name}: {entry['raw_bytes'] / 1e6:.1f} MB to {entry['zst_bytes'] / 1e6:.1f} MB")
            else:
                result['failed'].append(dict(name=name, reason=entry['reason']))
                say(f"{name}: FAILED, raw file kept: {entry['reason']}")
    except (SpaceFloor, InUse) as e:
        result['raw_bytes_after'], result['disk_bytes_after'] = disk_bytes(case, kept_names(case))
        e.result = result
        raise
    finally:
        os.close(lock)
    result['raw_bytes_after'], result['disk_bytes_after'] = disk_bytes(case, kept_names(case))
    return result


def _restore(zpath: Path, target: Path, line: dict, floor_bytes: float):
    """Write a compressed file's raw bytes to target by compress's steps: through target.part, hashed as it streams,
    checked against the manifest line, fsynced and renamed, with the raw file's own time."""
    _check_space(target.parent, line['raw_bytes'], floor_bytes)
    part = target.with_name(target.name + PART)
    digest, count = hashlib.sha256(), 0
    try:
        with part.open('wb') as out:
            for block in _decompressed(zpath):
                digest.update(block)
                count += len(block)
                out.write(block)
            out.flush()
            os.fsync(out.fileno())
            _drop_cache(out.fileno())
        if digest.hexdigest() != line['raw_sha256'] or count != line['raw_bytes']:
            raise RuntimeError(f'{zpath} decompresses to {count} bytes ({digest.hexdigest()[:16]}), not the '
                               f"{line['raw_bytes']} bytes ({line['raw_sha256'][:16]}) of its manifest line")
    except BaseException:
        part.unlink(missing_ok=True)
        raise
    if 'raw_mtime_ns' in line:
        os.utime(part, ns=(line['raw_mtime_ns'], line['raw_mtime_ns']))
    os.rename(part, target)
    _fsync_folder(target.parent)


def copy_raw(source, target, floor_bytes: float = FLOOR_BYTES) -> bool:
    """Copy a file kept raw or compressed to target as a raw file: shutil.copy2 of the raw file, or the .zst
    decompressed and checked against its manifest line (_restore). Returns whether it decompressed."""
    source, target = Path(source), Path(target)
    if source.exists():
        shutil.copy2(source, target)
        return False
    zpath = compressed(source)
    if not zpath.exists():
        raise FileNotFoundError(source)
    line = manifest(source.parent).get(source.name)
    if line is None:
        raise RuntimeError(f'{zpath} has no line in {source.parent / MANIFEST} to check it against')
    _restore(zpath, target, line, floor_bytes)
    return True


def decompress(case, names=None, floor_bytes: float = FLOOR_BYTES) -> list:
    """Write a case's named files (all its compressed ones by default) back raw beside their .zst, each checked
    against its manifest line; a file already raw is passed over. Returns the names written."""
    case = Path(case)
    if names is None:
        names = [p.name for p in glob(case, 'cm1*.dat') if compressed(p).exists()]
    lines, written = manifest(case), []
    for name in names:
        name = name[:-len(SUFFIX)] if name.endswith(SUFFIX) else name
        path = case / name
        if path.exists():
            continue
        if not compressed(path).exists():
            raise FileNotFoundError(f'{path}: neither it nor its .zst stands')
        if name not in lines:
            raise RuntimeError(f'{name}: no line in {case / MANIFEST} to check it against')
        _restore(compressed(path), path, lines[name], floor_bytes)
        written.append(name)
    return written


def restore(case, names, floor_bytes: float = FLOOR_BYTES) -> list:
    """The files a CM1 run is about to read, written back raw where they stand compressed (decompress); returns the
    names written."""
    return decompress(case, [n for n in names if not (Path(case) / n).exists()], floor_bytes)


def verify(case) -> dict:
    """Every compressed file of a case decompressed and checked against its manifest line: the names that match,
    differ or fail to decompress, the .zst files with no line, and the lines whose files are gone."""
    case = Path(case)
    lines = manifest(case)
    out = dict(case=case.name, ok=[], mismatched=[], unrecorded=[], gone=[], raw=[])
    for p in glob(case, 'cm1*.dat'):
        zpath = compressed(p)
        if not zpath.exists():
            continue
        line = lines.get(p.name)
        if line is None:
            out['unrecorded'].append(p.name)
            continue
        try:
            digest, count = _hash_decompressed(zpath)
        except RuntimeError as e:
            out['mismatched'].append(dict(name=p.name, reason=str(e)))
            continue
        if digest == line['raw_sha256'] and count == line['raw_bytes']:
            out['ok'].append(p.name)
        else:
            out['mismatched'].append(dict(name=p.name, reason=f'{count} bytes, {digest[:16]}; the line holds '
                                                              f"{line['raw_bytes']} bytes, {line['raw_sha256'][:16]}"))
    for name in sorted(lines):
        if not compressed(case / name).exists():
            out['raw' if (case / name).exists() else 'gone'].append(name)
    return out


def case_folder(run: str) -> Path:
    """A case named on the command line: its folder under runs/, or the path given."""
    return Path(run) if os.sep in run else RUNS / run


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest='command', required=True)
    c = sub.add_parser('compress', help="a case's 3-D outputs and restarts, one at a time")
    c.add_argument('runs', nargs='+')
    c.add_argument('--floor-gb', type=float, default=FLOOR_BYTES / 1e9, help='free space to keep beyond each file')
    d = sub.add_parser('decompress', help='files back raw beside their .zst (all of the case when none is named)')
    d.add_argument('run')
    d.add_argument('files', nargs='*')
    d.add_argument('--floor-gb', type=float, default=FLOOR_BYTES / 1e9)
    v = sub.add_parser('verify', help="every compressed file of a case against its manifest line")
    v.add_argument('run')
    args = parser.parse_args(argv)
    if args.command == 'compress':
        status = 0
        for run in args.runs:
            try:
                r = compress(case_folder(run), args.floor_gb * 1e9, report=lambda line: print(line, flush=True))
            except InUse as e:                                          # skipped; the other cases go on
                print(f'{run}: skipped, {e}', flush=True)
                status = 1
                continue
            except SpaceFloor as e:
                print(f'{run}: stopped, {e}', flush=True)
                return 2
            print(f"{r['case']}: {len(r['compressed'])} compressed, {len(r['finished'])} finished, {len(r['done'])} done "
                  f"before, {len(r['failed'])} failed, {len(r['left'])} left alone; {r['disk_bytes_before'] / 1e9:.2f} GB "
                  f"on the drive before, {r['disk_bytes_after'] / 1e9:.2f} GB after", flush=True)
            status = status or int(bool(r['failed'] or r['left']))
        return status
    if args.command == 'decompress':
        written = decompress(case_folder(args.run), args.files or None, args.floor_gb * 1e9)
        print(f'{args.run}: {len(written)} files written raw')
        return 0
    r = verify(case_folder(args.run))
    print(f"{r['case']}: {len(r['ok'])} match, {len(r['mismatched'])} differ, {len(r['unrecorded'])} without a line, "
          f"{len(r['gone'])} lines whose files are gone, {len(r['raw'])} lines whose files stand raw")
    for item in r['mismatched']:
        print(f"  {item['name']}: {item['reason']}")
    for name in r['unrecorded']:
        print(f'  {name}: no manifest line')
    return 1 if r['mismatched'] or r['unrecorded'] else 0


if __name__ == '__main__':
    sys.exit(main())
