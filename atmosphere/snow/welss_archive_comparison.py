"""Compare frozen Welss characteristics to the archive's unsmoothed profiles.

The physical integrations and initial one-sided-filter diagnostics remain intact.
This separate product records the archive evidence and direct raw comparisons.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from atmosphere.snow import welss_column as w
from atmosphere.snow import welss_comparison as suite

HERE=Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def mean_lattice(values,samples):
    """Compatibility with a mean of three-significant-digit ASCII quantities.

    The local decimal quantum is approximate when samples cross a power of ten.
    Compatibility identifies neither unique sample count nor sampling times.
    """
    a=np.asarray(values,dtype=float)
    a=a[np.isfinite(a)&(a>0)]
    quantum=10.**(np.floor(np.log10(a))-2.)
    scaled=a*samples/quantum
    distance=np.abs(scaled-np.round(scaled))
    tolerance=2.*np.spacing(a.astype(np.float32)).astype(float)*samples/quantum
    return dict(compatible=int(np.sum(distance<=tolerance)),count=int(a.size),
                compatible_fraction=float(np.mean(distance<=tolerance)))


def profile_metrics(model,reference,height):
    m,r,z=map(lambda x:np.asarray(x,dtype=float),(model,reference,height))
    valid=np.isfinite(m)&np.isfinite(r)&(r>0)&(z<w.SOURCE_TOP)
    rows={}
    for label,mask in [('all_populated',valid),('below_nucleation',valid&(z<w.SOURCE_BOTTOM)),
                       ('nucleation_belt',valid&(z>=w.SOURCE_BOTTOM))]:
        d=m[mask]-r[mask];relative=d/r[mask];peak=int(np.argmax(np.abs(relative)))
        rows[label]=dict(count=int(mask.sum()),relative_L1=float(np.abs(d).sum()/r[mask].sum()),
                         relative_bias=float(d.sum()/r[mask].sum()),
                         relative_RMSE=float(np.sqrt(np.mean(d*d))/np.mean(r[mask])),
                         maximum_absolute_relative_difference=float(np.abs(relative[peak])),
                         maximum_relative_difference_height_m=float(z[mask][peak]))
    return rows


def make_product(physical_path,archive_path):
    physical_path,archive_path=map(Path,(physical_path,archive_path))
    paths=[physical_path,archive_path,Path(__file__),HERE/'welss_archive_sources.json']
    before={str(p):sha(p) for p in paths}
    physical=json.loads(physical_path.read_text())
    if physical['producer']!=suite.identity():
        raise ValueError('frozen physical-result producer identity differs from available source')
    reference=w.read_archive(archive_path)
    data=reference['arrays'];z=data['height'];rows={}
    for name,case in physical['cases'].items():
        if case['producer']['files']!=w.source_identity() or case['reference']['sha256']!=reference['sha256']:
            raise ValueError('case producer or archive identity differs')
        if not np.array_equal(case['height_m'],z):
            raise ValueError('archive and model height bins differ')
        rows[name]={mode:{key:profile_metrics(profile['raw_binned'][key],data[key],z)
                          for key in w.MOMENTS[:4]}
                    for mode,profile in case['profiles'].items()}
    evidence=dict(boundary_bins_m=[0.,20.,2960.,2980.],
                  all_totals_finite_at_boundary_bins=bool(all(np.isfinite(data[key][[0,1,148,149]]).all() for key in w.MOMENTS[:4])),
                  interpretation='Finite domain-end and source-end bins exclude the verified five-point opt0 pass with source zero-to-missing handling. Direct raw-bin comparison is used.',
                  arithmetic={key:{str(n):mean_lattice(data[key],n) for n in (4,30,31,32,61,62,93)} for key in w.MOMENTS[:4]},
                  arithmetic_limit='The 31-snapshot lattice also admits multiples such as62 and93. Agreement supports, but does not identify uniquely, a31-snapshot five-hour window at600s cadence; raw time series and executable remain unavailable.')
    detailed_changes={}
    base=physical['cases']['base']
    for name,case in physical['cases'].items():
        if name=='base':
            continue
        detailed_changes[name]={mode:{key:profile_metrics(base['profiles'][mode]['raw_binned'][key],profile['raw_binned'][key],z)
                                      for key in w.MOMENTS}
                                for mode,profile in case['profiles'].items()}
    raw_changes={name:{mode:changes['raw_binned'] for mode,changes in values.items()}
                 for name,values in physical['refinement_changes'].items()}
    result=dict(schema='terluna.atmosphere.welss-archive-comparison/1',
                producer=dict(file='atmosphere/snow/welss_archive_comparison.py',sha256=before[str(Path(__file__))],
                              source_ledger_sha256=before[str(HERE/'welss_archive_sources.json')]),
                physical_product=dict(file='atmosphere/snow/results/welss_column.json',sha256=before[str(physical_path)]),
                reference={k:v for k,v in reference.items() if k!='arrays'},
                evidence='Comparison between deterministic declared-source calculations and archived author numerical profiles, not empirical validation.',
                reading_rule='Use these direct raw-profile metrics for archive interpretation. The earlier one-sided-filter metrics in the frozen physical product are processing diagnostics, not matched archive comparisons.',
                archive_processing_evidence=evidence,comparison=rows,raw_refinement_changes=raw_changes,
                regional_refinement_changes=detailed_changes,
                averaging='All sampled physical trajectories reach ground before5h; continuous last5h, four final snapshots and final-time deterministic profiles coincide. An unknown archived window still changes stochastic sampling.',
                remaining=['Exact archived executable and preprocessing configuration.',
                           'Origin of residual number-density differences inside the nucleation belt.',
                           'Aggregation population comparison and empirical constraints.'])
    if before!={str(p):sha(p) for p in paths} or physical['producer']!=suite.identity():
        raise RuntimeError('source or input changed while comparing archive')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference',type=Path,required=True)
    parser.add_argument('--physical-result',type=Path,default=HERE/'results/welss_column.json')
    parser.add_argument('--output',type=Path,default=HERE/'results/welss_archive_comparison.json')
    args=parser.parse_args()
    result=make_product(args.physical_result,args.reference)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result['comparison']['base']['last5h'],indent=2))


if __name__=='__main__':
    main()
