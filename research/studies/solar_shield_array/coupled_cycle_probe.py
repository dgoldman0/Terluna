"""Budget-gated coupled check of the cycle-tail initializer, never acceptance."""
import json
import signal
import time
import numpy as np
from protection.dynamics.natural_pattern import local_coverage
from .coupled_common import HERE, limit, inputs, write, save, read_raw, spline, state_at, digest, length
from .coupled_cycle import CycleDebt, transplant
from .coupled_stable import StableExperiment
from .coupled_refine import truncate


def main():
    limit();cpu=time.process_time();wall=time.monotonic()
    products={n:json.loads((HERE/'results'/('coupled_'+n+'.json')).read_text())
              for n in ['benchmark','optimization','validation']}
    recovered=products['validation'].get('recovery',{})
    recovery_charge=sum(recovered.get(k,0.) for k in ['completed_execution_cpu_s','additional_failed_export_allowance_s','additional_finalize_attempt_cpu_s'])
    spent=products['optimization']['previous_numerical_charge_s']+sum(products[n]['resources']['cpu_s'] for n in ['optimization','validation'])+recovery_charge
    remaining=1800-spent
    # Preserve at least 45 s for final scientific extraction/report products.
    if remaining<280:raise RuntimeError('Cycle-tail probe exceeds remaining allocation')
    signal.alarm(int(remaining-45))
    ex=StableExperiment();bench=products['benchmark'];opt=products['optimization']
    x=np.array(products['validation']['selected_x'])
    selected=opt['selected']
    if products['validation']['validated_case']=='last_rejected_contact_rebuild':
        selected=next(r['replay'] for r in reversed(opt['iterations']) if 'replay' in r)
    prefix=read_raw(selected)
    _,_,producer=inputs(['research/studies/solar_shield_array/coupled_cycle.py',
                        'research/studies/solar_shield_array/coupled_stable.py',
                        'research/studies/solar_shield_array/coupled_refine.py',
                        'research/studies/solar_shield_array/coupled_cycle_probe.py'])
    producer['inputs'].update({'coupled_'+n+'.json':digest(HERE/'results'/('coupled_'+n+'.json')) for n in products})
    pending=dict(schema='terluna.research.coupled-cycle-probe/1',producer=producer,completed=False,
        accepted_return=False,accepted_next_service=False,accepted_cycle=False,accepted_fleet=False,
        remaining_cpu_s_before_probe=remaining,selected_x=x.tolist(),
        scope='Diagnostic tail through failed contacts; no physical service or cycle acceptance.')
    def expired(signum,frame):
        pending['stopped_reason']='Reserved final reporting allocation; tail probe exceeded remaining stage budget'
        write('coupled_cycle_probe.json',pending,cpu,wall)
        raise SystemExit(2)
    signal.signal(signal.SIGALRM,expired)
    write('coupled_cycle_probe.json',pending,cpu,wall)
    row,raw,_=ex.evaluate(x,prefix['state'][-1],float(prefix['t'][-1]),ex.end,step=240.)
    row.update(save('selected_cycle_tail_diagnostic',**raw));row['restoration_only']=True
    service=read_raw(bench['service']);full_base=read_raw(bench['evaluations'][0])
    base=truncate(full_base,opt['horizon_s']);base_tail=read_raw(opt['tail'])
    base_tail={k:(np.r_[full_base[k][full_base['t']>=opt['horizon_s']],v[1:]] if k in
        ['t','state','frames','power_W','attitude_acceleration','first_intercept_W','redirected_W','translation','distance_m'] else v)
        for k,v in base_tail.items()}
    debt=CycleDebt(ex,base,base_tail,service)
    predicted_residual=debt.residual(prefix,x)[1]
    actual=state_at(spline(raw),debt.dates);targets=[];cmd=ex.command(x)
    for date,phase,reference in zip(debt.dates,debt.phase_dates,debt.base_tail):
        state=state_at(spline(service if phase<=21600 else prefix),phase)
        targets.append(transplant(state,cmd(phase).as_matrix(),cmd(phase,1),
            cmd(date).as_matrix(),cmd(date,1),reference.mean(axis=0)))
    residual=actual-np.array(targets);error=residual-predicted_residual
    out=dict(schema='terluna.research.coupled-cycle-probe/1',producer=producer,completed=False,
        accepted_return=False,accepted_next_service=False,accepted_cycle=False,accepted_fleet=False,
        remaining_cpu_s_before_probe=remaining,tail=row,selected_x=x.tolist(),
        scope='Diagnostic propagation through already failed contacts; no physical execution or service credit after the first clearance stop. Tests the terminal initializer, not a feasible cycle.',
        stages=['next_service_start','next_service_end','subsequent_departure'],dates_s=debt.dates.tolist(),
        actual_maximum_position_debt_m=length(residual[:,:,:3]).max(axis=1).tolist(),
        actual_maximum_velocity_debt_m_s=length(residual[:,:,3:]).max(axis=1).tolist(),
        terminal_proxy_maximum_position_error_m=length(error[:,:,:3]).max(axis=1).tolist(),
        terminal_proxy_maximum_velocity_error_m_s=length(error[:,:,3:]).max(axis=1).tolist(),
        actual_normalized_cycle_residual_rms=float(np.sqrt(np.mean((residual/debt.scale)**2))),
        predicted_normalized_cycle_residual_rms=debt.summary(prefix,x)['normalized_rms'],
        hypothetical_control_energy_through_subsequent_departure_J=
            bench['service']['electrical_energy_J']+selected['electrical_energy_J']+row['electrical_energy_J'],
        actual_executed_energy_J=products['validation']['account']['electrical_energy_J'])
    write('coupled_cycle_probe.json',out,cpu,wall)
    coverage=[]
    for t in np.linspace(ex.arrival,ex.arrival+21600.,13):
        state=state_at(spline(raw),t)
        coverage.append(local_coverage(state,ex.env.at(t),state.mean(axis=0),limb_count=16,interior=8,rotation=.137))
    out['diagnostic_next_service']=dict(dates=13,sources_per_date=25,
        minimum_source_coverage=min(r['minimum_source_coverage'] for r in coverage),
        credited_service_s=0.,physical_passage_completed=False)
    out['completed']=True;write('coupled_cycle_probe.json',out,cpu,wall)
    print('cycle probe',out['actual_maximum_position_debt_m'],out['terminal_proxy_maximum_position_error_m'],
          out['diagnostic_next_service'],out['resources']['cpu_s'],flush=True)


if __name__=='__main__':main()
