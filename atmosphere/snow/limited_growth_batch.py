"""Batch the frozen limited-growth reconstruction without changing its arithmetic.

Each row is one layer with the same mass grid. Array reconstruction is evaluated
for all layers together, but scalar vapour caps and each row's dot-product
reduction retain the frozen implementation's operation order. This is an
execution optimization, not a new transport or physical approximation.
"""
from __future__ import annotations

import numpy as np

from atmosphere.snow.limited_growth import LimitedGrowth


class BatchLimitedGrowth:
    """Evaluate existing LimitedGrowth operators on one layer-by-mass array.

    Returns (number_transfers_2d, deposited_rates_per_layer, minimum_CFL_limit).
    Callers seeking bitwise parity with a previous sequential ledger should sum
    deposited rates in the same left-to-right order, not a vector reduction.
    """
    def __init__(self, operators):
        if not operators or not all(isinstance(o, LimitedGrowth) for o in operators):
            raise ValueError('a nonempty sequence of LimitedGrowth operators is required')
        first=operators[0]
        if not all(np.array_equal(o.masses,first.masses) for o in operators):
            raise ValueError('batched operators must share exactly the same mass grid')
        self.masses=first.masses.copy()
        self.widths=first.widths.copy()
        self.centroid_spacing=first.centroid_spacing.copy()
        self.face_offset=(first.face_masses-first.centroids).copy()
        self.spacing=first.spacing.copy()
        self.face_growth=np.stack([o.face_growth for o in operators])
        self.cfl=np.array([o.cfl for o in operators])
        self.layers=len(operators)

    def evaluate(self, number, production):
        n=np.asarray(number,dtype=float)
        if n.shape!=self.face_growth.shape or not np.isfinite(n).all() or (n<0).any():
            raise ValueError('population must be finite, nonnegative and match layer/mass shape')
        production=np.asarray(production,dtype=float)
        if production.ndim==0:
            production=np.full(self.layers,float(production))
        if production.shape!=(self.layers,) or not np.isfinite(production).all() or (production<0).any():
            raise ValueError('production must be finite nonnegative scalar or per-layer vector')
        density=n/self.widths
        slope=np.zeros_like(density)
        if self.masses.size>2:
            left=(density[:,1:-1]-density[:,:-2])/self.centroid_spacing[:-1]
            right=(density[:,2:]-density[:,1:-1])/self.centroid_spacing[1:]
            h_left,h_right=self.centroid_spacing[:-1],self.centroid_spacing[1:]
            centered=(h_right*left+h_left*right)/(h_left+h_right)
            slope[:,1:-1]=np.where(left*right>0,np.sign(centered)*np.minimum(np.abs(centered),
                                             np.minimum(2*np.abs(left),2*np.abs(right))),0.)
        face=density+slope*self.face_offset
        face[:,:-1]=np.minimum(np.maximum(face[:,:-1],np.minimum(density[:,:-1],density[:,1:])),
                              np.maximum(density[:,:-1],density[:,1:]))
        transfers=self.face_growth*face
        # Keep the scalar implementation's vector length, layout and dot order.
        demand=np.array([float(np.dot(self.spacing,row)) for row in transfers])
        factor=np.array([min(1.,float(p)/float(d)) if d>0 else 1. for p,d in zip(production,demand)])
        transfers*=factor[:,None]
        deposited=np.array([float(np.dot(self.spacing,row)) for row in transfers])
        ratios=np.divide(n,transfers,out=np.full_like(n,np.inf),where=transfers>0)
        limits=self.cfl*np.min(ratios,axis=1)
        return transfers,deposited,float(limits.min())
