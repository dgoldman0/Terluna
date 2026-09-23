"""Reproduce B1 without network calls or repository writes.

Existing frame records are reused only when the renderer's complete input/source
identity matches. Scientific integration is deliberately expensive.
"""
from __future__ import annotations
import argparse,os,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]

def main():
 p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=4);p.add_argument('--viewer-only',action='store_true');a=p.parse_args()
 if a.workers<1:p.error('workers must be positive')
 env={**os.environ,'OPENBLAS_NUM_THREADS':'1','NUMBA_NUM_THREADS':'1'}
 def run(*args):
  print('RUN',*args,flush=True);subprocess.run([sys.executable,*map(str,args)],cwd=ROOT,env=env,check=True)
 if not a.viewer_only:
  run('tools/b1/test_b1.py');run('tools/b1/check_spectral.py')
  groups=[(['clear','fair'],[60,30,6,0,-2,-6],16),(['clear','fair'],[12],64),(['high'],[6,-2],16),(['thin','thick'],[6],16)]
  for cases,suns,spp in groups:
   run('tools/b1/render_scene.py','--scenarios',*cases,'--suns',*suns,'--width',160,'--height',96,'--spp',spp,'--workers',a.workers,'--output','data/b1/generated/final')
  groups=[(['clear','fair'],[60,30,12,6,0,-2,-6]),(['high'],[6,-2]),(['thin','thick'],[6])]
  for cases,suns in groups:
   run('tools/b1/render_scene.py','--scenarios',*cases,'--suns',*suns,'--projection','hemisphere','--width',80,'--height',40,'--spp',12,'--workers',a.workers,'--output','data/b1/generated/final')
   run('tools/b1/diagnostics.py','--scenarios',*cases,'--suns',*suns,'--size',28,'--spp',12,'--workers',a.workers)
  run('tools/b1/render_scene.py','--scenarios','fair','--suns',12,'--black','--width',160,'--height',96,'--spp',16,'--workers',a.workers,'--output','data/b1/generated/controls')
 run('build_b1_painter.py');run('tools/b1/check_browser.py');run('tools/b1/validate_b1.py')
if __name__=='__main__':main()
