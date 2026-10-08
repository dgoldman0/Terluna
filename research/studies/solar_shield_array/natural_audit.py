"""Numerical, prefix-clearance, optical and inventory audit of bounded results."""
import time,resource,signal
import numpy as np
from protection.dynamics.service_search import free_command
from protection.dynamics.square_distance import distance_guard
from protection.dynamics.service_margin import service_margin
from protection.dynamics.optical import length
from .joint_common import spline,state_at
from .natural_common import *


def read(name,producer):
    p=HERE/'results'/name;producer['inputs'][name]=digest(p);return json.loads(p.read_text())


def raw(entry):
    p=ROOT/entry['raw_path'];assert digest(p)==entry['raw_sha256'];return np.load(p)


def main():
    signal.alarm(90);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    before=time.monotonic();cpu=time.process_time()
    arrays,producer=setup(['research/studies/solar_shield_array/natural_audit.py'])
    fine=read('natural_free_slow_fine.json',producer);coarse=read('natural_free_slow_coarse.json',producer)
    account=read('natural_account.json',producer);scout=read('natural_fast_scout.json',producer);original=read('natural_scout.json',producer)
    env=environment(3.);cmd=free_command(env,fine['arrival_s']+21600,arrival=fine['arrival_s'],turn_start=6.,turn_hours=6.)
    entry=fine['runs'][1];path=spline(raw(entry));dates=np.unique(np.r_[np.arange(21600.,43200.,60.),43200.]);states=state_at(path,dates)
    coarse_path=spline(raw(coarse['runs'][1]));error=float(length(states[:,:,:3]-coarse_path(dates)).max())+entry['reconstruction_error_m']+coarse['runs'][1]['reconstruction_error_m']
    guard=distance_guard(dates,states,cmd,point_error=error)
    first=spline(raw(fine['runs'][0]));reference=read('joint_regression.json',producer);ref=spline(raw(reference['runs'][0]));worst=None;count=0
    for t in np.linspace(0,21600,25):
        y=state_at(first,t);centre=state_at(ref,t).mean(axis=0)
        for angle in np.arange(12)*2*np.pi/12+.319:
            result=service_margin(y,cmd(t).as_matrix(),env.at(t),centre,[np.cos(angle),np.sin(angle)])
            count+=1
            if worst is None or result['margin_m']<worst['margin_m']:worst=dict(time_s=float(t),angle_rad=float(angle),**result)
    # Quantify only the reduced centre-force approximation, not missing shadows.
    a=raw(original['cases'][0]);b=raw(scout['cases'][0]);bp=spline(b)
    reduced_error=[]
    for stop in [43200.,172800.,864000.]:
        mask=a['t']<=stop;ts=a['t'][mask]
        reduced_error.append(dict(through_s=stop,maximum_position_difference_m=float(length(a['state'][mask,:,:3]-bp(ts)).max())))
    face=read('natural_face_refined_coarse.json',producer);proposal=read('natural_refine.json',producer);proposed=spline(raw(proposal));max_difference=0.
    for run in face['runs']:
        values=raw(run);max_difference=max(max_difference,float(length(values['state'][:,:,:3]-proposed(values['t'])).max()))
    measured=next(x for x in account['cases'] if x['case']=='free_slow')
    out=dict(schema='terluna.research.natural-audit/1',producer=producer,completed=True,accepted_return=False,accepted_cycle=False,accepted_fleet=False,
        first_twelve_hours=dict(replayed_position_difference_m=error,clearance=guard,electrical_energy_J=measured['first_twelve_hours_energy_J'],loaded_baseline_J=51.859e12,
            energy_reduction_fraction=1-measured['first_twelve_hours_energy_J']/51.859e12,zero_electric_translation=measured['translation_mean_m_s']==0,
            passed_local_prefix=error<50 and guard['minimum_conditional_all_pair_clearance_m']>100 and worst['margin_m']>0 and measured['overloaded_members']==0),
        initial_service_adversary=dict(probes=count,worst=worst,continuous_certificate=False),finite_area_vs_centre_scout=reduced_error,
        cheap_ray_candidate_maximum_coupled_vs_proposal_position_difference_m=max_difference,
        replacement_service_required_from_s=21600.,next_service_executed=False,handovers_executed=0,
        inventory_scope='Each listed scout cadence supplies a necessary time-capacity inventory lower bound for the same six-hour service credit. Spatial matching and handover are unproved. No failed-prefix energy is extrapolated to operating fleet demand.')
    write('natural_audit.json',out,before,cpu);print(out['first_twelve_hours'],out['resources'],flush=True)
if __name__=='__main__':main()
