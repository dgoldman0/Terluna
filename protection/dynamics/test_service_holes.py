from types import SimpleNamespace
import numpy as np
from .service_holes import uncovered_rows
from .cycling import sail_basis


def test_continuous_receiver_audit_finds_an_interior_hole():
    sun=np.array([1.5e11,0.,0.]);zero=np.zeros(3)
    sample=dict(positions=dict(sun=sun,earth=np.array([4e8,0.,0.])),sun_v=zero,sun_a=zero,earth_v=zero,earth_a=zero)
    env=SimpleNamespace(at=lambda t:sample)
    _,b,c=sail_basis(sun);ix,iy=np.meshgrid(np.arange(-9,10),np.arange(-9,10))
    q=np.array([15e6,0.,0.])+8500*(ix.ravel()[:,None]*b+iy.ravel()[:,None]*c)
    y=np.c_[q,np.zeros_like(q)]
    base=lambda ts:np.broadcast_to(y,(len(np.atleast_1d(ts)),361,6)).copy()
    maps=SimpleNamespace(sol=lambda ts:np.zeros((54,len(np.atleast_1d(ts)))))
    constraints,full=uncovered_rows(env,base,base,maps,0.,1)
    assert full['minimum_source_coverage']>1-1e-9 and not constraints
    def missing(ts):
        result=base(ts);result[:,180,:3]+=1e7*b
        return result
    constraints,holes=uncovered_rows(env,base,missing,maps,0.,1)
    assert holes['minimum_source_coverage']<.95 and constraints
