"""One successive correction using measured coupled displacement defects."""
import json
import time
import signal
import resource

import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.spatial import cKDTree

from shared.provenance import constants_used, constants_changed
from protection.dynamics.energy_coast import coast_command
from protection.dynamics.pattern_departure import expansion_maps, expansion_state
from protection.dynamics.departure_control import control_modes, mode_impulses, choose_cut, minimize_impulse
from protection.dynamics.coast_control import signed_impulse
from protection.dynamics.pattern_return import encounter_scan, crossing_witnesses, beam_and_rate_screen
from .energy_screen import RUN
from .fleet_run import environment, digest
from .cycling_search import ROOT, HERE


def main():
    signal.alarm(180);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    before=time.monotonic()
    p=json.loads((HERE/'results/energy_screen.json').read_text())
    tested=json.loads((HERE/'results/energy.json').read_text())
    sources={**tested['producer']['source_hashes'],
        'protection/dynamics/coast_control.py':digest(ROOT/'protection/dynamics/coast_control.py'),
        'research/studies/solar_shield_array/energy_correct.py':digest(__file__)}
    for path,h in sources.items():
        if digest(ROOT/path)!=h:raise ValueError('Source changed: '+path)
    if constants_changed(tested['producer']['constants']):raise ValueError('Constants changed')
    entry=tested['runs'][0];case=p['cases'][entry['case']]
    if digest(ROOT/entry['raw_path'])!=entry['raw_sha256']:raise ValueError('Raw changed')
    raw=np.load(ROOT/entry['raw_path']);initial=raw['state'][0]
    start,end,burn=p['start_s'],p['end_s'],p['burn_s'];env=environment(1.)
    maps=expansion_maps(env,initial,start,end,burn);modes=control_modes(env,start)
    command=coast_command(env,start,end,**case['schedule'])
    original=np.array(case['parameters_m_s']);original_impulses=mode_impulses(modes,original)
    proposal=expansion_state(maps,initial,original_impulses,raw['t'])
    difference=raw['state']-proposal
    defect=CubicHermiteSpline(raw['t'],difference[:,:,:3],difference[:,:,3:])
    last=float(raw['t'][-1]);inverse=np.linalg.inv(maps.sol(last)[6:42].reshape(6,6))

    def corrected_state(parameters,times):
        times=np.atleast_1d(times)
        result=expansion_state(maps,initial,mode_impulses(modes,parameters),times)
        for k,t in enumerate(times):
            if t<=last:
                result[k]+=np.c_[defect(t),defect(t,1)]
            else:
                phi=maps.sol(t)[6:42].reshape(6,6)
                result[k]+=difference[-1]@(phi@inverse).T
        return result

    parameters=original.copy();cuts=[];seen=set();records=[]
    times=np.arange(start,end+1,120.);fine_t=np.arange(start,end+1,30.)
    for iteration in range(4):
        states=corrected_state(parameters,times)
        for t,state,frame in zip(times,states,command(times).as_matrix()):
            pairs=cKDTree(state[:,:3]).query_pairs(14500.,output_type='ndarray')
            if not len(pairs):continue
            projected=(state[pairs[:,0],:3]-state[pairs[:,1],:3])@frame
            scores=(abs(projected)-[10000.,10000.,0.]).max(axis=1)
            response=maps.sol(t)[42:].reshape(6,3)[:3]
            for k in np.argsort(scores)[:12]:
                i,j=pairs[k];gradient=frame.T@response@(modes[i]-modes[j])
                cut=choose_cut(t,i,j,projected[k],gradient,parameters,400.)
                if cut is None:continue
                key=(round(t,1),int(i),int(j),cut['axis'],cut['sign'])
                if key not in seen:cuts.append(cut);seen.add(key)
        previous=parameters.copy()
        parameters,fit=minimize_impulse(modes,parameters,cuts,burn)
        row=dict(iteration=iteration,nonnegative_optimization=fit);records.append(row)
        if not fit['success']:
            parameters,fit=signed_impulse(modes,previous,cuts,burn)
        row['optimization']=fit
        if not fit['success']:break
        fine=corrected_state(parameters,fine_t)
        scan=encounter_scan(fine_t,fine,command,clearance=350.,retain=8)
        roots=crossing_witnesses(fine_t,fine,command,lambda t,ids:corrected_state(parameters,[t])[0,ids],retain=8)
        beam=beam_and_rate_screen(env,fine_t,fine,command)
        row.update(encounters=scan,crossings=roots,beam_and_rate=beam)
        if not scan['sampled_pair_violations'] and not roots['crossing_near_encounters'] and beam['minimum_sampled_beam_margin_deg']>0:break
    passed=bool(fit['success'] and not records[-1].get('encounters',{}).get('sampled_pair_violations',1)
        and not records[-1]['crossings']['crossing_near_encounters'] and records[-1]['beam_and_rate']['minimum_sampled_beam_margin_deg']>0)
    path=RUN/'corrected_proposal.npz'
    np.savez_compressed(path,initial=initial,t=fine_t,state=corrected_state(parameters,fine_t),
        impulses=mode_impulses(modes,parameters),frames=command(fine_t).as_matrix())
    out=dict(schema='terluna.research.energy-coast-correction/1',
        producer=dict(source_hashes=sources,constants=constants_used(sources),
            inputs={n:digest(HERE/'results'/n) for n in ['energy_screen.json','energy.json']}),
        start_s=start,end_s=end,burn_s=burn,case=case['index'],schedule=case['schedule'],
        original_parameters_m_s=original.tolist(),parameters_m_s=parameters.tolist(),
        iterations=records,cuts=len(cuts),proposal_pass=passed,accepted_departure=False,accepted_fleet=False,
        nominal_attitude_m_s=case['nominal_attitude_m_s'],translation_mean_m_s=fit.get('mean_delta_v_m_s'),
        defect_scope='Measured coupled-minus-variational defect through the clearance stop. Remaining 0.51 h uses gravity transport of that endpoint defect. Control sensitivities still use gravity variational response; a complete coupled retry is mandatory.',
        measured_defect_end_s=last,maximum_measured_position_defect_m=float(np.linalg.norm(difference[:,:,:3],axis=-1).max()),
        raw_path=str(path.relative_to(ROOT)),raw_sha256=digest(path),elapsed_s=time.monotonic()-before,
        max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
    (HERE/'results/energy_correction.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(stage='correction',passed=passed,parameters=out['parameters_m_s'],
        mean_impulse_m_s=out['translation_mean_m_s'],elapsed_s=out['elapsed_s'])),flush=True)


if __name__=='__main__':main()
