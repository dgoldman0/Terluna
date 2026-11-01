"""Run and retain independent numerical controls for the Welss Earth comparison."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from atmosphere.snow import welss_column as w

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def cases():
    base = dict(mass_nodes=16,height_nodes_per_cell=2,bin_nodes=3,rtol=1e-7,max_step_m=20.)
    return dict(base=base,mass=dict(base,mass_nodes=32),height=dict(base,height_nodes_per_cell=4),
                bins=dict(base,bin_nodes=5),ode=dict(base,rtol=1e-9,max_step_m=10.))


def identity():
    return dict(suite_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                model_files=w.source_identity())


def compare_cases(rows):
    """Relative integrated profile changes; no automatic convergence verdict."""
    base=rows['base'];z=np.asarray(base['height_m']);result={}
    for name,row in rows.items():
        if name=='base':
            continue
        if row['producer'] != base['producer'] or row['reference'] != base['reference']:
            raise ValueError('comparison cases have different producers or archive inputs')
        result[name]={}
        for mode in base['profiles']:
            result[name][mode]={}
            for kind in ('raw_binned','filtered_binned'):
                metrics={}
                for key in w.MOMENTS:
                    a=np.asarray(base['profiles'][mode][kind][key],dtype=float)
                    b=np.asarray(row['profiles'][mode][kind][key],dtype=float)
                    valid=np.isfinite(a)&np.isfinite(b)&(z<w.SOURCE_TOP)
                    denominator=np.sum(np.abs(b[valid]))
                    metrics[key]=dict(relative_L1=float(np.sum(np.abs(a[valid]-b[valid]))/denominator),
                                      relative_bias=float(np.sum(a[valid]-b[valid])/denominator))
                result[name][mode][kind]=metrics
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference',type=Path,required=True)
    parser.add_argument('--output',type=Path,default=HERE/'results/welss_column.json')
    parser.add_argument('--cache',type=Path,default=ROOT/'research/runs/welss_column')
    parser.add_argument('--wall-limit-s',type=float,default=600.)
    args=parser.parse_args()
    before=identity();reference=w.read_archive(args.reference)
    args.cache.mkdir(parents=True,exist_ok=True)
    rows={}
    for name,settings in cases().items():
        print(f'Welss numerical control: {name}',flush=True)
        result=w.integrate(**settings,wall_limit_s=args.wall_limit_s,progress=True)
        if identity()!=before or hashlib.sha256(args.reference.read_bytes()).hexdigest()!=reference['sha256']:
            raise RuntimeError('suite producer or archive input changed during integration')
        row=w.product(result,reference,before['model_files'])
        (args.cache/(name+'.json')).write_text(json.dumps(row,indent=2,allow_nan=False)+'\n')
        rows[name]=row
        print(json.dumps(dict(case=name,diagnostics=row['diagnostics'],comparison=row['profiles']['last5h']['comparison'])),flush=True)
    product=dict(schema='terluna.atmosphere.welss-comparison/1',producer=before,
                 evidence='Independent declared-source Earth deposition/settling comparison to archived author simulations, with separate numerical refinements; not empirical validation.',
                 reading_rule='Compare each independently refined case with base. Relative L1 uses the refined profile absolute integral below 3000 m. A single refinement per variable bounds tested sensitivity, not a formal convergence order.',
                 cases=rows,refinement_changes=compare_cases(rows))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(product,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
