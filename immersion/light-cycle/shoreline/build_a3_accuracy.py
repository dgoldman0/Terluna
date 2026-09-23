"""Build the offline A3 data inspector. All data are local and hash-checked."""
from pathlib import Path
import hashlib,json,sys
ROOT=Path(__file__).resolve().parent

def build(output=None):
 output=output or ROOT/'Open_Moon_A3_Accuracy_Lab.html'
 bench=ROOT/'data/a3/generated/benchmark.json';valid=ROOT/'A3_VALIDATION.json';display=ROOT/'data/a3/generated/display.json'
 for p in [bench,valid,display]:
  if not p.is_file():raise FileNotFoundError(p)
 v=json.loads(valid.read_text())
 if hashlib.sha256(bench.read_bytes()).hexdigest()!=v['benchmark_sha256']:raise ValueError('Benchmark bytes differ from validation')
 text=(ROOT/'index.a3-accuracy-lab.html').read_text()
 for tag,p in [('__A3_BENCHMARK__',bench),('__A3_VALIDATION__',valid),('__A3_DISPLAY__',display),('__A3_SCRIPT__',ROOT/'src/a3-accuracy-lab.js')]:
  s=p.read_text()
  if tag!='__A3_SCRIPT__' and '</script' in s.lower():raise ValueError('Unsafe inline data delimiter')
  text=text.replace(tag,s)
 output.write_text(text)
 print(json.dumps(dict(output=str(output),bytes=output.stat().st_size,sha256=hashlib.sha256(output.read_bytes()).hexdigest())))
if __name__=='__main__':build()
