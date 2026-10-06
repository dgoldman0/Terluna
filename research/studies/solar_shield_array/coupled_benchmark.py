"""Reproduce the loaded slow turn and measure exact coupled derivatives."""
import signal
import time
import numpy as np
from .coupled_common import Experiment, DIMENSION, limit, inputs, write, save, length


def main():
    limit(); signal.alarm(550)
    cpu = time.process_time(); wall = time.monotonic()
    ex = Experiment(); _, _, producer = inputs(['research/studies/solar_shield_array/coupled_benchmark.py'])
    out = dict(schema='terluna.research.coupled-benchmark/1', producer=producer,
        completed=False, evaluations=[], accepted_return=False, accepted_cycle=False,
        scope='Finite perturbations replay all 361 coupled tiles. Post-contact continuation is restoration-only.')
    x = np.zeros(DIMENSION)
    service, raw, sol = ex.evaluate(x, ex.arrays['initial_epoch_state'], 0., 21600.)
    service.update(save('baseline_service', **raw)); out['service']=service
    initial = raw['state'][-1]
    write('coupled_benchmark.json', out, cpu, wall)
    print('service',service['cpu_s'],flush=True)
    end = 14.5*3600
    values = []
    for label, value in [('base', 0.), ('repeat', 0.), ('direction_h', .2), ('direction_half_h', .1)]:
        x = np.zeros(DIMENSION); x[0] = value
        row, raw, sol = ex.evaluate(x, initial, 21600., end)
        row['name']=label; row.update(save('benchmark_'+label, **raw))
        out['evaluations'].append(row); values.append(raw)
        write('coupled_benchmark.json', out, cpu, wall)
        print(label, row['cpu_s'], row['clearance_events_s'][:1], flush=True)
    base, repeat, full, half = values
    a=(full['state']-base['state'])/.2; b=(half['state']-base['state'])/.1
    out['repeatability']=dict(maximum_position_m=float(length(repeat['state'][:,:,:3]-base['state'][:,:,:3]).max()),
        maximum_velocity_m_s=float(length(repeat['state'][:,:,3:]-base['state'][:,:,3:]).max()))
    out['directional_convergence']=dict(parameter='row_parity_normal_impulse',
        physical_steps_m_s=[.02,.01], maximum_position_derivative_disagreement_m=float(length(a[:,:,:3]-b[:,:,:3]).max()),
        maximum_response_m=float(length(full['state'][:,:,:3]-base['state'][:,:,:3]).max()),
        half_step_prediction_error_m=float(length(half['state'][:,:,:3]-base['state'][:,:,:3]-.1*a[:,:,:3]).max()))
    out['completed']=True; write('coupled_benchmark.json',out,cpu,wall)
    print(out['repeatability'],out['directional_convergence'],flush=True)


if __name__=='__main__':
    main()
