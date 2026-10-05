"""Ten-day all-member independent scout; no terminal lattice or held formation."""
import signal,resource,time
import numpy as np
from scipy.optimize import brentq
from protection.dynamics.service_search import free_command,natural_opportunities,cadence_account
from protection.dynamics.cycle_control import independent_path
from protection.dynamics.square_distance import closest_squares
from protection.dynamics.natural_pattern import local_coverage
from protection.dynamics.optical import length
from .natural_common import *


def main():
    signal.alarm(400);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    before=time.monotonic();cpu=time.process_time()
    arrays,producer=setup(['protection/dynamics/cycle_control.py','research/studies/solar_shield_array/natural_scout.py'])
    env=environment(10.5);initial=arrays['first_service_end_state'];ratio=arrays['mass_ratio'];end=10*86400.
    out=dict(schema='terluna.research.natural-scout/1',producer=producer,cases=[],completed=False,accepted_return=False,accepted_cycle=False,
        scope='Independent finite-source 361-member proposal trajectories. Trajectories beyond a sampled contact are mathematical diagnostics, not executed physical sequences. Hypothetical Sun-facing coverage requires an executed turn for feathered families.')
    for spec in [dict(name='face',tilt=0.,axis=0),dict(name='edge_a',tilt=90.,axis=0),dict(name='edge_b',tilt=90.,axis=1)]:
        cmd=free_command(env,end,tilt=spec['tilt'],axis=spec['axis'])
        sol=independent_path(env,initial,21600,end,cmd,mass_ratio=ratio)
        def state(t):return sol.sol(t).reshape(361,6)
        times=np.unique(np.r_[np.arange(21600,64800,300.),np.arange(64800,end,1800.),end]);ys=sol.sol(times).T.reshape(-1,361,6)
        first=None;minimum=np.inf
        for k,(t,y) in enumerate(zip(times,ys)):
            distance,pair,projection=closest_squares(y[:,:3],cmd(t).as_matrix());minimum=min(minimum,distance)
            if first is None and distance<100:
                lo=times[max(0,k-1)]
                def gap(x):return closest_squares(state(x)[:,:3],cmd(x).as_matrix())[0]-100
                hit=brentq(gap,lo,t,xtol=.001) if gap(lo)>0 else t
                dd,pp,pr=closest_squares(state(hit)[:,:3],cmd(hit).as_matrix())
                first=dict(time_s=float(hit),pair=pp.tolist(),distance_m=float(dd),projection_m=pr.tolist())
        opportunities=natural_opportunities(env,lambda t:state(t).mean(axis=0),86400,end-21600)
        for item in opportunities:
            t=item['midpass_s'];y=state(t);c=y.mean(axis=0)
            coverage=local_coverage(y,env.at(t),c,limb_count=8,interior=0)
            item.update(hypothetical_Sun_facing_midpass_coverage=coverage,cadence=cadence_account(item['arrival_s'],mass_kg=float((5e6*ratio).sum())))
        row=dict(spec=spec,first_clearance_failure=first,minimum_sampled_distance_m=float(minimum),opportunities=opportunities,
            minimum_lunar_radius_m=float(length(ys[:,:,:3]).min()),nfev=sol.nfev,
            **save('scout_'+spec['name'],t=times,state=ys,mass_ratio=ratio))
        out['cases'].append(row);write('natural_scout.json',out,before,cpu)
        print(spec['name'],'first',first,'opportunities',[(round(x['arrival_s']/3600,2),round(x['hypothetical_Sun_facing_midpass_coverage']['minimum_source_coverage'],3)) for x in opportunities],time.monotonic()-before,flush=True)
    out['completed']=True;write('natural_scout.json',out,before,cpu)
if __name__=='__main__':main()
