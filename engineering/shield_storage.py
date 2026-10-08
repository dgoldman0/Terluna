"""Per-member chronological local-bus storage, including losses and power caps."""
import numpy as np


def required_capacity(times,load_W,generation_W,charge=.95,discharge=.95,depth=.8):
    """Nominal energy needed with a full initial battery and surplus recharge.

    The maximum drawdown of the cumulative net-energy path gives each isolated
    bus's capacity; it grants no transfer from other buses. Generation remains
    a caller-specified scenario. A separate rating limits charge/discharge.
    """
    load=np.asarray(load_W);gen=np.asarray(generation_W)
    net=(load[:-1]+load[1:]-gen[:-1]-gen[1:])/2
    delta=np.maximum(net,0)/discharge-np.maximum(-net,0)*charge
    draw=np.zeros(load.shape[1]);maximum=draw.copy()
    for d,dt in zip(delta,np.diff(times)):
        draw=np.maximum(0,draw+d*dt);maximum=np.maximum(maximum,draw)
    return maximum/depth


def dispatch(times,load_W,generation_W,nominal_J,limit_W,charge=.95,discharge=.95,depth=.8):
    times=np.asarray(times);load=np.asarray(load_W);gen=np.asarray(generation_W)
    if load.shape!=gen.shape or load.ndim!=2 or np.any(np.diff(times)<=0):raise ValueError('Invalid chronological arrays')
    if np.any(load<0) or np.any(gen<0):raise ValueError('Negative power')
    capacity=np.asarray(nominal_J)*depth;state=capacity.copy();trace=[state.copy()]
    missing=np.zeros_like(state);curtailed=missing.copy();loss=missing.copy();first=None
    supplied=np.zeros_like(state);generated=np.zeros_like(state);demand=np.zeros_like(state)
    for k,dt in enumerate(np.diff(times)):
        l=(load[k]+load[k+1])/2;g=(gen[k]+gen[k+1])/2;net=g-l
        charging=np.minimum(np.maximum(net,0),limit_W)
        added=np.minimum(charging*dt*charge,capacity-state)
        requested=np.minimum(np.maximum(-net,0),limit_W)
        drawn=np.minimum(requested*dt/discharge,state)
        delivered=drawn*discharge;state+=added-drawn
        miss=np.maximum(-net,0)*dt-delivered
        if first is None and np.any(miss>1e-5):first=float(times[k])
        missing+=miss;curtailed+=np.maximum(net,0)*dt-added/charge
        loss+=added*(1/charge-1)+drawn*(1-discharge)
        supplied+=l*dt-miss;generated+=g*dt;demand+=l*dt
        trace.append(state.copy())
    # Initial energy + generated bus = final stored + supplied + losses + spill.
    residual=capacity+generated-state-supplied-loss-curtailed
    return dict(state_J=np.array(trace),unserved_J=missing,curtailed_J=curtailed,
        storage_loss_J=loss,first_unserved_s=first,balance_residual_J=residual,
        supplied_J=supplied,generation_J=generated,demand_J=demand)
