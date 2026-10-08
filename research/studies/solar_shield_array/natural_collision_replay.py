"""Check natural-release collision witnesses without interpolated positions."""
import json
import numpy as np
from scipy.integrate import solve_ivp

from shared import constants as K
from shared.provenance import constants_used, constants_changed
from protection.dynamics.natural_traffic import NaturalTraffic
from protection.dynamics.active_global import square_pair_separation
from protection.dynamics.optical import unit
from .natural_traffic_run import identities
from .fleet_run import environment, digest
from .cycling_search import ROOT, HERE
from .cycling_analysis import compact_series_json


class ExactMembers(NaturalTraffic):
    def path(self,t,plane,slot=None,phase=None,parity=None):
        if slot is None or t!=self.replay_time:
            raise ValueError('Exact witness states are available only at their replay date')
        indices=np.array([self.lookup[(int(j),int(k))] for j,k in zip(plane,slot)])
        state=self.exact_states[indices]
        return state,self.env.gravity(t,state[:,:3]),self.exact_tangent[indices]


def main():
    src={**identities(),str(__file__).split(str(ROOT)+'/')[-1]:digest(__file__)}
    path=HERE/'results/natural_traffic.json';p=json.loads(path.read_text())
    if p['producer']['source_hashes']!=identities() or constants_changed(p['producer']['constants']):
        raise ValueError('Natural diagnostic identity mismatch')
    env=environment(30);rows=[];epsilon=2e-6
    for row in p['nominal_return_subset_intersections']:
        if not row['example_ids']:continue
        ids=np.unique(np.array(row['example_ids']).reshape(-1,2),axis=0)
        f=ExactMembers(env);j,k=ids.T
        phase=2*np.pi*k/f.count[j]+f.phase[j]
        # Two neighbouring initial phases supply the actual advected tangent.
        # All three states are independently integrated at the same epoch.
        initial=f.initial_states(j[:,None],phase[:,None]+epsilon*np.array([-1,0,1]),(k%2)[:,None])
        def rhs(t,y):
            y=y.reshape(-1,6)
            return np.column_stack([y[:,3:],env.gravity(t,y[:,:3])]).ravel()
        t=row['time_days']*K.JULIAN_DAY
        sol=solve_ivp(rhs,[0,t],initial.ravel(),method='DOP853',t_eval=[t],
            max_step=300.,rtol=2e-12,atol=np.tile([1e-5]*3+[1e-9]*3,initial.size//6))
        if not sol.success:raise RuntimeError(sol.message)
        truth=sol.y[:,-1].reshape(-1,3,6)
        f.exact_states=truth[:,1];f.exact_tangent=unit(truth[:,2,:3]-truth[:,0,:3]);f.replay_time=t
        f.lookup={tuple(map(int,x)):i for i,x in enumerate(ids)}
        pairs=np.array([[f.lookup[tuple(a)],f.lookup[tuple(b)]] for a,b in row['example_ids']])
        d=f.facets(t,j,k)
        separation=square_pair_separation(d['state'][:,:3],d['a'],d['b'],pairs,f.layout.side)
        rows.append(dict(day=row['time_days'],integer_members=len(ids),example_ids=row['example_ids'],
            maximum_separating_axis_values_m=separation.tolist(),
            confirmed_square_intersections=int(np.count_nonzero(separation<-1e-5)),
            minimum_beam_margins_deg=np.rad2deg(d['margins'].min(axis=0)).tolist(),nfev=sol.nfev))
    out=dict(schema='terluna.research.natural-collision-replay/1',
        producer=dict(source_hashes=src,constants=constants_used(src),input_product_sha256=digest(path)),
        max_step_s=300.,rtol=2e-12,phase_difference_rad=epsilon,series=rows,
        evidence='Independent common-epoch IVPs of actual integer collision witnesses and neighbouring phases. Finite-square intersections use the independently propagated centres and advected tangents.',
        scope='Confirms individual intersections in the point-gravity, ideal-attitude model; no global collision census, avoidance manoeuvre or swept-path certificate is claimed.')
    if {name:digest(ROOT/name) for name in src}!=src:raise ValueError('Sources changed during replay')
    (HERE/'results/natural_collision_replay.json').write_text(compact_series_json(out))
    print(json.dumps(out['series']),flush=True)


if __name__=='__main__':main()
