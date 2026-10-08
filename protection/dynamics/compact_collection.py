"""The existing first-intercept ledger using exact compact shadow unions."""
import numpy as np

from shared import constants as K
from .compact_shadow import compact_illumination
from .natural_pattern import retarded_sun
from .collection import require_uneclipsed_tiles
from .cycling import ray_geometry
from .fleet import sun_points
from .optical import length, unit


def compact_collection(env,t,state,frame,suns=64,rotation=.317,clear_side=9890.):
    q=state[:,:3];sample=env.at(t)
    require_uneclipsed_tiles(ray_geometry(q,sample),10000.)
    first=np.zeros(len(q));front=first.copy();back=first.copy();unshadowed=first.copy()
    force=np.zeros_like(q)
    for source in sun_points(retarded_sun(sample),suns,rotation):
        ray=source-q;cosine=unit(ray)@frame[:,2]
        power=K.SOLAR_CONSTANT*(K.AU/length(ray))**2*clear_side**2*abs(cosine)/suns
        unique=power*compact_illumination(q,frame,source,clear_side)
        unshadowed+=power;first+=unique
        front+=np.where(cosine>=0,unique,0.);back+=np.where(cosine<0,unique,0.)
        force-=(2*env.central_fraction/K.SPEED_OF_LIGHT*unique*cosine)[:,None]*frame[:,2]
    return dict(first_intercept_bolometric_equivalent_W=first,front_first_intercept_W=front,
        back_first_intercept_W=back,unshadowed_aperture_W=unshadowed,
        redirected_band_optical_W=env.central_fraction*first,
        mutual_shadow_equivalent_loss_W=unshadowed-first,reflected_force_N=force)
