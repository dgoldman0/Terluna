"""Render the long-range candidate and resolve its narrow low-horizon view."""
from pathlib import Path
import json
from climate.crm.cloud_scene import export
from climate.crm.cloud_columns import sha256
from illumination.cloud_light.render import render
from illumination.cloud_light.detail import low_horizon

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
RUNS=ROOT/'research/runs/optical_comfort/evening'


def main():
    source=HERE/'results/deep_evenings.json';p=json.loads(source.read_text())
    rows=[]
    for pick in p['selected_scenes']:
        name=f"deep_{pick['case']}_{pick['observer_column']}_{pick['snapshot']}"
        section=RUNS/'sections'/f'{name}.npz'
        admission=export(pick['case'],pick['snapshot'],pick['observer_column'],pick['source_sha256'],section,half_columns=160)
        # This selected equatorial target lies west, along decreasing ring index.
        if pick['case']!='ring_equator' or pick['cloud_longitude_deg']>=pick['observer_longitude_deg']:
            raise ValueError('Verify the camera direction for this new long-range candidate')
        parent=render(section,RUNS/'scenes'/f'{name}_overview.npz',180.,photons=256,shape=(8,16),surface_photons=65536)
        parent.update(selection=pick['selection'],initial_pick=pick,delta_hours=0,section=admission)
        (RUNS/'scenes'/f'{name}_overview.json').write_text(json.dumps(parent,indent=2)+'\n')
        print('Detailed horizon view',name,flush=True)
        detail=low_horizon(parent,RUNS/'scenes'/f'{name}_detail.npz')
        (RUNS/'scenes'/f'{name}_detail.json').write_text(json.dumps(detail,indent=2,allow_nan=False)+'\n')
        rows.append(detail)
        print('Finished',name,'open lux',detail['surface_lux'],'low opening',detail['recess_lux'],flush=True)
    product=dict(schema='terluna.research.cloud-evening-scenes/1',mode='deep',
        producer=dict(file=str(Path(__file__).relative_to(ROOT)),sha256=sha256(__file__),deep_input_sha256=sha256(source)),
        evidence='Full spectral cloud transport for the expanded long-range search. A 0.5–15 degree view resolves the distant cloud; a separate 0.5–5 degree opening probes dark foreground viewing.',scenes=rows)
    (HERE/'results/evening_scenes_deep.json').write_text(json.dumps(product,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
