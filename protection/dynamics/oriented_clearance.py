"""All-pair oriented-square clearance with conservative projection pruning."""
import numpy as np
from scipy.spatial import cKDTree
from .oriented_tiles import square_distance
from .optical import length


def projection_lower_bounds(delta,fa,fb,side=10000.):
    a=np.swapaxes(fa,1,2);b=np.swapaxes(fb,1,2)
    axes=np.concatenate([a,b,np.cross(a[:,:,None,:],b[:,None,:,:]).reshape(-1,9,3)],axis=1)
    norm=length(axes);axes=axes/np.maximum(norm,1e-100)[:,:,None]
    radius=side/2*(np.abs(np.einsum('nkj,nja->nka',axes,fa[:,:,:2])).sum(axis=2)+np.abs(np.einsum('nkj,nja->nka',axes,fb[:,:,:2])).sum(axis=2))
    gap=np.abs(np.einsum('nkj,nj->nk',axes,delta))-radius
    return np.maximum(np.max(np.where(norm>1e-12,gap,0.),axis=1),0.)


def closest_oriented(q,frames,side=10000.):
    tree=cKDTree(q);d,j=tree.query(q,k=2);i=int(np.argmin(d[:,1]));j=int(j[i,1])
    best=square_distance(q[i],frames[i],q[j],frames[j],side);pair=[i,j]
    pairs=tree.query_pairs(np.sqrt(2)*side+best+1e-6,output_type='ndarray')
    lower=projection_lower_bounds(q[pairs[:,0]]-q[pairs[:,1]],frames[pairs[:,0]],frames[pairs[:,1]],side)
    for k in np.argsort(lower):
        if lower[k]>=best:break
        i,j=pairs[k];value=square_distance(q[i],frames[i],q[j],frames[j],side)
        if value<best:best=value;pair=[int(i),int(j)]
        if best<1e-8:break
    return float(best),pair
