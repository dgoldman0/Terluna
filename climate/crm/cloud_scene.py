"""Export a geographical CM1 ring section for optical scene calculations."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from climate.crm.cloud_columns import dry_density, read_cm1_constants,sha256,SPECIES,NUMBERS
from climate.crm.ring_analysis import case_geometry,read_snapshot,hour_angle
from climate.crm.cm1_run import RUNS
from shared.constants import SYNODIC_MONTH_DAYS,MOON_RADIUS

SCHEMA='terluna.climate.cloud-section/1'


def unit_vectors(lat,lon):
    a,b=np.radians(lat),np.radians(lon)
    return np.stack((np.cos(a)*np.cos(b),np.cos(a)*np.sin(b),np.sin(a)),axis=-1)


def export(case_name,snapshot,observer_column,expected_hash,output,half_columns=100):
    case=RUNS/case_name
    geo=case_geometry(case);record=geo['record']
    if record['reference']['run']!='A28_dim5_moon' or not record.get('path'):
        raise ValueError('A corrected flat great-circle ring is required')
    scalar=case/f'cm1out_t{snapshot:06d}_s.dat'
    if sha256(scalar)!=expected_hash:
        raise ValueError('Selected cloud snapshot differs from regional input admission')
    d=read_snapshot(case,snapshot)
    constants,source_hashes=read_cm1_constants(record)
    rho,temperature=dry_density(d['prs'],d['th'],d['qv'],constants)
    lat,lon=[np.array(record['path'][k]) for k in ('lat','lon')]
    vectors=unit_vectors(lat,lon)
    n=len(lat);i=observer_column
    offsets=np.arange(-half_columns,half_columns+1)
    indices=(i+offsets)%n
    up=vectors[i]
    tangent=vectors[(i+1)%n]-vectors[(i-1)%n]
    tangent-=np.dot(tangent,up)*up;tangent/=np.linalg.norm(tangent)
    cross=np.cross(up,tangent)
    basis=np.stack((tangent,cross,up))
    time=(snapshot-1)*record['configuration']['output_s']/86400
    h=hour_angle(geo,time*86400,SYNODIC_MONTH_DAYS*86400)[i]
    sun_lon=np.radians(lon[i]-h)
    sun=basis@np.array([np.cos(sun_lon),np.sin(sun_lon),0.])
    metadata=dict(schema=SCHEMA,case=case_name,snapshot=snapshot,time_day=time,observer_column=i,
                  observer_latitude_deg=float(lat[i]),observer_longitude_deg=float(lon[i]),
                  observer_hour_angle_deg=float(h),observer_sun_elevation_deg=float(np.degrees(np.arcsin(sun[2]))),
                  dx_m=geo['dx'],radius_m=MOON_RADIUS,droplets_cm3=record['configuration']['droplets_cm3'],
                  source=dict(file=str(scalar),sha256=expected_hash,case_sha256=sha256(case/'case.json')),
                  producer=dict(file='climate/crm/cloud_scene.py',sha256=sha256(__file__)),
                  thermodynamic_source_hashes=source_hashes,
                  evidence='Saved evolving CM1 section in the second lunar cycle; the simulation resolves the along-ring '
                           'and vertical dimensions. Cross-ring extent is supplied by the optical scenario.',
                  reading_rule='Local Cartesian basis: x along the ring, y across it, z radially outward. '
                               'Each column is bounded by great-circle planes, each vertical cell by concentric shells. '
                               'The observer is at the centre column, simulated sea level.')
    if not SYNODIC_MONTH_DAYS<=time<=2*SYNODIC_MONTH_DAYS:
        raise ValueError('Second-cycle scene required')
    arrays={k:d[k][:,indices] for k in (*SPECIES,*NUMBERS,'prs','th','qv')}
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output,metadata=json.dumps(metadata),height_edges_m=geo['zw'],
                        column_index=indices,column_offset=offsets,latitude_deg=lat[indices],longitude_deg=lon[indices],
                        density_kg_m3=rho[:,indices].astype('f4'),temperature_k=temperature[:,indices].astype('f4'),
                        sun_local=sun,basis=basis,**arrays)
    return dict(path=str(output),sha256=sha256(output),schema=SCHEMA,metadata=metadata)
