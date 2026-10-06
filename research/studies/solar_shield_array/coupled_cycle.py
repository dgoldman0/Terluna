"""A repeated command schedule and explicitly approximate cycle-tail debt.

The baseline tail is all-member coupled propagation. Transport of a changed
prefix to that tail uses a gravity STM only and receives no feasibility credit.
This terminal initializer includes next service and the following departure;
accepted cycle claims require their coupled replay and optical tests.
"""
import numpy as np
from scipy.spatial.transform import Rotation, RotationSpline
from protection.dynamics.service_search import free_command, smooth_step
from protection.dynamics.packing_control import response_maps, matrices
from .coupled_common import Experiment, spline, state_at


class CycleExperiment(Experiment):
    def __init__(self):
        super().__init__()
        self.end = self.arrival+46800.

    def command(self,x):
        duration=(6.+.5*x[6])*3600
        first=free_command(self.env,self.end,arrival=self.arrival,turn_start=6.,turn_hours=duration/3600)
        knots=[self.arrival+21600.,self.arrival+21600.+duration,self.end]
        times=np.unique(np.r_[first.times,knots]);times=times[times<=self.end]
        angle=-np.pi/2*smooth_step((times-self.arrival-21600.)/duration)
        frames=first(times).as_matrix()@Rotation.from_rotvec(np.c_[angle,np.zeros((len(times),2))]).as_matrix()
        return RotationSpline(times,Rotation.from_matrix(frames))

    def controls(self,x):
        controls,windows=super().controls(x)
        return np.repeat(controls,2,axis=1),np.r_[windows,windows+self.arrival]


def transplant(state, old_frame, old_rate, new_frame, new_rate, centre):
    relative=state-state.mean(axis=0)
    q=relative[:,:3]@old_frame
    v=relative[:,3:]@old_frame-np.cross(old_rate,q)
    return centre+np.c_[q@new_frame.T,(v+np.cross(new_rate,q))@new_frame.T]


class CycleDebt:
    def __init__(self,ex,prefix,tail,service):
        self.ex=ex;self.prefix=spline(prefix);self.tail=spline(tail);self.service=spline(service)
        self.horizon=float(prefix['t'][-1]);self.reference_end=prefix['state'][-1]
        self.cmd=ex.command(np.zeros(8))
        self.dates=np.array([ex.arrival,ex.arrival+21600.,ex.end])
        self.phase_dates=np.array([0.,21600.,ex.end-ex.arrival])
        self.base_tail=state_at(self.tail,self.dates)
        # A repeated departure pulse acts on the later departure too.
        windows=np.array([[ex.arrival+21600.,ex.arrival+43200.]])
        maps=response_maps(ex.env,lambda t:state_at(self.tail,t).mean(axis=0),
                           self.horizon,ex.end,windows,np.eye(3))
        self.phi=matrices(maps,self.dates,1)[:,:,:6]
        self.response=matrices(maps,self.dates,1)[:,:,6:]
        self.scale=np.array([1e5]*3+[5.]*3)

    def residual(self,raw,x):
        command=self.ex.command(x);prefix=spline(raw)
        change=raw['state'][-1]-self.reference_end
        predicted=self.base_tail+np.einsum('tij,nj->tni',self.phi,change)
        u=self.ex.controls(x)[0][:,0]
        predicted+=np.einsum('tij,nj->tni',self.response,u)
        targets=[]
        for date,phase,reference in zip(self.dates,self.phase_dates,self.base_tail):
            state=state_at(self.service if phase<=21600 else prefix,phase)
            targets.append(transplant(state,command(phase).as_matrix(),command(phase,1),
                command(date).as_matrix(),command(date,1),reference.mean(axis=0)))
        residual=predicted-np.array(targets)
        return (residual/self.scale).ravel(),residual

    def summary(self,raw,x):
        vector,residual=self.residual(raw,x)
        return dict(normalized_rms=float(np.sqrt(np.mean(vector**2))),
            dates_s=self.dates.tolist(),stages=['next_service_start','next_service_end','subsequent_departure'],
            maximum_position_debt_m=np.linalg.norm(residual[:,:,:3],axis=-1).max(axis=1).tolist(),
            maximum_velocity_debt_m_s=np.linalg.norm(residual[:,:,3:],axis=-1).max(axis=1).tolist(),
            model='Coupled baseline tail; gravity-only transport of changed prefix and repeated burn. Diagnostic terminal initializer, not a feasible return or service certificate.')
