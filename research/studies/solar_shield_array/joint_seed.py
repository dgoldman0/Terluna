"""Persist exact canonical states for the next bounded local experiment."""
import time,hashlib
import numpy as np
from .joint_common import *


def main():
    before=time.monotonic();cpu=time.process_time();p,producer=setup(['research/studies/solar_shield_array/joint_seed.py'])
    entries=p['closure_validation.json']['runs'];service=load_raw(entries[0]);departure=load_raw(entries[1])
    bridge=json.loads((HERE/'results/joint_bridge.json').read_text())['cases'][0]
    producer['inputs']['joint_bridge.json']=digest(HERE/'results/joint_bridge.json')
    arrays=dict(initial_epoch_state=service['state'][0],first_service_end_state=service['state'][-1],
        original_twelve_hour_end_state=departure['state'][-1],mass_ratio=service['mass_ratio'])
    out=dict(schema='terluna.research.joint-seed/1',producer=producer,completed=True,
        units=dict(position='m, Moon-relative ICRF',velocity='m/s',time='s from 2026-10-04 TDB'),
        arrays={k:v.tolist() for k,v in arrays.items()},
        float64_little_endian_hashes={k:hashlib.sha256(np.asarray(v,dtype='<f8').tobytes()).hexdigest() for k,v in arrays.items()},
        selected_arrival_s=bridge['arrival_s'],selected_spec=bridge['spec'],
        reading_rule='Exact round-trip initial states and unchanged hardware ratios, retained independently of ignored trajectories. joint_bridge.build reconstructs the corridor from first_service_end_state, mass_ratio, selected_arrival_s and selected_spec. This seed alone does not constitute a replay or provide operational acceptance.')
    write('joint_seed.json',out,before,cpu)
if __name__=='__main__':main()
