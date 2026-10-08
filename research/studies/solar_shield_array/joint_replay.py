"""Source-bound tighter replay, with a 90-second maximum step.

The common runner is preserved byte-for-byte with its coarse execution. This
entry point supplies 16 rather than 8 solar sources, 32 rather than 16 area
nodes per axis at partial eclipses, 90 rather than 120 second maximum steps,
and rtol 2e-12 rather than 2e-10. It records this wrapper in the producer hash.
"""
from . import joint_execute as runner
from .joint_common import HERE,digest
import json

_original_setup=runner.setup
_original_solve=runner.solve_ivp

def setup(extra):
    return _original_setup([*extra,'research/studies/solar_shield_array/joint_replay.py'])

def solve(*args,**kwargs):
    kwargs['max_step']=90.
    return _original_solve(*args,**kwargs)

def main():
    import sys
    if '--fine' not in sys.argv:sys.argv.append('--fine')
    runner.setup=setup;runner.solve_ivp=solve
    runner.main()
    path=HERE/'results/joint_corridor_108_fine.json'
    out=json.loads(path.read_text())
    out['replay_settings']=dict(suns=16,body_area_grid_per_axis=32,max_step_s=90.,rtol=2e-12,
        control_source='Unchanged coarse inverse-dynamics commands; replay begins at epoch zero with no state resets.')
    path.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
if __name__=='__main__':main()
