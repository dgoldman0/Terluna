"""Export/restore compact optimizer state independently of ignored run files."""
import argparse
import base64
import hashlib
import io
import json
import time
import numpy as np
from .coupled_common import HERE, ROOT, inputs, read_raw, write, digest


def restore(path=None):
    data=json.loads((path or HERE/'results/coupled_restart.json').read_text())
    for name,expected in data['producer']['source_hashes'].items():
        if digest(ROOT/name)!=expected:raise ValueError('Restart source changed: '+name)
    payload=base64.b64decode(data['arrays_npz_base64'],validate=True)
    if hashlib.sha256(payload).hexdigest()!=data['payload_sha256']:
        raise ValueError('Restart payload changed')
    with np.load(io.BytesIO(payload),allow_pickle=False) as raw:
        arrays={k:raw[k] for k in raw.files}
    for name,array in arrays.items():
        h=hashlib.sha256(np.asarray(array,dtype='<f8').tobytes()).hexdigest()
        if h!=data['float64_hashes'][name]:raise ValueError('Restart array changed: '+name)
    return data,arrays


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--inspect',action='store_true');args=ap.parse_args()
    if args.inspect:
        data,arrays=restore();print(data['continuation']);print({k:v.shape for k,v in arrays.items()});return
    cpu=time.process_time();wall=time.monotonic()
    opt=json.loads((HERE/'results/coupled_optimization.json').read_text())
    bench=json.loads((HERE/'results/coupled_benchmark.json').read_text())
    selected=read_raw(opt['selected']);service=read_raw(bench['service'])
    model=read_raw(opt['checkpoint'])
    _,_,producer=inputs(['research/studies/solar_shield_array/coupled_restart.py',
        'research/studies/solar_shield_array/coupled_cycle.py','protection/dynamics/coupled_control.py',
        'research/studies/solar_shield_array/coupled_stable.py','research/studies/solar_shield_array/coupled_refine.py',
        'research/studies/solar_shield_array/coupled_contacts.py',
        'research/studies/solar_shield_array/coupled_optimize.py'])
    producer['inputs'].update({n:digest(HERE/'results'/n) for n in ['coupled_benchmark.json','coupled_optimization.json']})
    last=opt['iterations'][opt['checkpoint']['last_iteration']]
    arrays=dict(initial_optimization_state=service['state'][-1],selected_terminal_state=selected['state'][-1],
        basis=selected['basis'],x=np.array(opt['selected_x']),
        linearization_anchor_x=np.array(last['base_x']),cuts=model['cuts'],gaps=model['gaps'],
        gap_jac=model['gap_jac'],terminal_jac=model['terminal_jac'],
        endpoint_state_jac=model['state_jac'][-1],endpoint_frame_jac=model['frame_jac'][-1],
        last_replayed_step=np.array(last['step']))
    packed=io.BytesIO();np.savez_compressed(packed,**arrays);payload=packed.getvalue()
    out=dict(schema='terluna.research.coupled-restart/1',producer=producer,completed=True,
        arrays_npz_base64=base64.b64encode(payload).decode('ascii'),payload_sha256=hashlib.sha256(payload).hexdigest(),
        array_shapes={k:list(v.shape) for k,v in arrays.items()},
        float64_hashes={k:hashlib.sha256(np.asarray(v,dtype='<f8').tobytes()).hexdigest() for k,v in arrays.items()},
        continuation=opt['checkpoint'],selected_terminal_time_s=float(selected['t'][-1]),
        selected_terminal_restoration_only=bool(opt['selected'].get('clearance_events_s')),
        cycle_debt=opt['final_cycle_debt'],accepted_return=False,accepted_cycle=False,
        reading_rule='Exact float64 restart controls, coupled starting/terminal arrays, basis and trust radius. The saved terminal trace can follow failed contacts and is an optimization initializer, not an executed physical state. Retained local half-spaces/contact Jacobian and terminal Jacobian use linearization_anchor_x. Endpoint response includes the latest measured secant. Call restore() without ignored runs to recover them. Repropagate selected x from the actual initial_optimization_state at hour 6 and rebuild date-dependent responses when extending the horizon. Full saved trajectories/Jacobians remain source-bound reproducible raw products, not prerequisites for reading this checkpoint.')
    write('coupled_restart.json',out,cpu,wall);restore();print('restart bytes',(HERE/'results/coupled_restart.json').stat().st_size)


if __name__=='__main__':main()
