"""Measure actual holes in the saved giant-square probe and maneuver scales.

Distances are to simultaneous finite-square footprints at the same ephemeris
date. Dilating those footprints is only an optimistic reach envelope: it adds
effective area and ignores assignments, collisions and velocity changes. No
dilated result is labelled achieved coverage or a feasible transfer.
"""
import json
from pathlib import Path

import numpy as np
import shapely

from shared import constants as K
from shared.provenance import constants_used
from protection.dynamics.fleet import fleet_acceleration,square_vertices,projected_polygons,sun_points
from protection.dynamics.cycling import sail_basis
from protection.dynamics.optical import unit,arriving_ray_vectors
from .fleet_run import read_raw,environment,identities,digest
from .cycling_search import HERE,ROOT,SCENARIO,CENTRAL
from .cycling_analysis import compact_series_json


def main():
    name='dense_double_1000km'; raw,meta=read_raw(name)
    states=raw['state']; t=raw['t']; env=environment(meta['config']['days'])
    radius=4*K.MOON_RADIUS; side=1e6
    axis=np.linspace(-radius,radius,129); x,y=np.meshgrid(axis,axis)
    xy=np.column_stack([x.ravel(),y.ravel()]); xy=xy[np.sum(xy**2,axis=1)<=radius**2]
    boundary=np.arange(1024)*2*np.pi/1024
    witness_xy=np.concatenate([xy,radius*np.column_stack([np.cos(boundary),np.sin(boundary)])])
    points=shapely.points(witness_xy)
    budgets=np.array([0,1000,10000,100000,500000,1000000.])
    dates=[0.,7.5,15.,22.5,29.25,30.]; records=[]
    for day in dates:
        i=int(np.argmin(abs(t-day*K.JULIAN_DAY))); state=states[i]; sample=env.at(t[i])
        diag=fleet_acceleration(env,t[i],state,radius,side,diagnostics=True)
        sun,_=arriving_ray_vectors(sample['positions']['sun'],sample['positions']['earth'],
            sample['sun_v'],sample['earth_v'],sample['sun_a'],sample['earth_a'])
        u=unit(sun); _,b,c=sail_basis(sun)
        tau=np.maximum(state[:,:3]@u,0)/K.SPEED_OF_LIGHT
        q=state[:,:3]-(state[:,3:]+sample['moon_v'])*tau[:,None]
        keep=q@u>0
        vertices=square_vertices(q[keep],diag['normal'][keep],np.broadcast_to(sun,q[keep].shape),side-110.)
        distances=[]; max_gap=-1; witness=None
        for k,source in enumerate(sun_points(sun,16,.217)):
            polygons=projected_polygons(vertices,source,u,b,c)
            polygons=polygons[shapely.is_valid(polygons)&(shapely.area(polygons)>1e-8)]
            matches,distance=shapely.STRtree(polygons).query_nearest(points,return_distance=True,all_matches=False)
            ordered=np.empty(len(points)); ordered[matches[0]]=distance
            distances.append(ordered[:len(xy)])
            j=int(np.argmax(ordered))
            if ordered[j]>max_gap:
                max_gap=float(ordered[j]); witness=dict(point_xy_km=(witness_xy[j]/1000).tolist(),solar_direction_index=k)
        distances=np.concatenate(distances)
        record=dict(day=day,mean_unintercepted_fraction_on_spatial_grid=float(np.mean(distances>1e-6)),
            maximum_identified_distance_to_a_footprint_m=max_gap,witness=witness,
            area_solar_ray_distance_percentiles_m=dict(zip(['50','90','95','99'],map(float,np.percentile(distances,[50,90,95,99])))),
            optimistic_dilated_ray_reach_fraction={str(d):float(np.mean(distances<=d+1e-6)) for d in budgets})
        records.append(record); print(json.dumps(record),flush=True)
    prop=SCENARIO['propulsion']; ve=prop['exhaust_velocity_m_s']; eta=prop['efficiency']
    cant=np.cos(np.deg2rad(prop['cant_deg'])); a_sail=2*CENTRAL*K.SOLAR_CONSTANT/(K.SPEED_OF_LIGHT*.05)*2/(3*np.sqrt(3))
    distances=[100.,10000.,100000.,500000.,max(r['maximum_identified_distance_to_a_footprint_m'] for r in records)]
    scales=[]
    for distance in distances:
        for hours in [6.,24.,72.]:
            duration=hours*3600; acc=4*distance/duration**2; dv=4*distance/duration
            scales.append(dict(distance_m=distance,hours=hours,symmetric_rest_to_rest_acceleration_m_s2=acc,
                maneuver_delta_v_m_s=dv,mean_power_per_10km_tile_W=5e6*acc*ve/(2*eta*cant),
                propellant_per_10km_tile_kg=5e6*dv/(ve*cant),
                fraction_of_maximum_ideal_transverse_photon_acceleration=acc/a_sail))
    sources={**identities(),'research/studies/solar_shield_array/active_gap_scale.py':digest(__file__),
        'research/studies/solar_shield_array/cycling_analysis.py':digest(HERE/'cycling_analysis.py')}
    result=dict(schema='terluna.research.active-gap-scale/1',
        producer=dict(source_hashes=sources,constants=constants_used(sources),raw_sha256=meta['sha256']),
        evidence='Spatial distance witnesses on actual simultaneous finite-Sun footprints of the earlier giant-square geometry probe; ideal maneuver scales only',
        parameters=dict(case=name,receiver_grid_width=129,receiver_interior_points=len(xy),boundary_points=1024,solar_directions=16),
        records=records,ideal_maneuver_scales=scales,propulsion_assumptions=prop,
        maximum_transverse_photon_acceleration_m_s2=a_sail,
        reading_rule='The largest identified hole is a sampled geometric witness, not a global maximum certificate. Dilation is an optimistic reach envelope with invented extra area, not achieved coverage. Rest-to-rest scales omit orbital gravity, initial/final velocity constraints, finite-extent force errors, target assignment, collisions and return legs; they are neither an optimized solution nor a rigorous thrust lower bound. Long lead times can exploit orbital phase changes.',
        global_coverage_demonstrated=False,delivered_power_W=None)
    (HERE/'results/active_gap_scale.json').write_text(compact_series_json(result))


if __name__=='__main__':main()
