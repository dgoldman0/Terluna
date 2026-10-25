"""Reproduce the isolated capacitance case's trial rejection and startup checks."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from atmosphere.snow import welss_column as w
from atmosphere.snow import welss_compact_capacitance as c
from atmosphere.snow import welss_archive_comparison as archive

HERE=Path(__file__).resolve().parent


def main():
    bound=[Path(__file__),HERE/'welss_compact_capacitance.py',HERE/'welss_archive_comparison.py']
    identity={p.name:archive.sha(p) for p in bound};source=w.source_identity()
    masses,_,heights,_=w.seed_quadrature();z0=float(heights[16])
    original_geom,original_solver=w.geometry,w.solve_ivp
    cases={};changes={}
    try:
        w.geometry=lambda mass:c.compact_geometry(mass,original_geom)
        for rtol in [1e-7,1e-9]:
            ends=[]
            for factor in [1.,.5]:
                counter=dict(invalid_rk_stages=0)
                w.solve_ivp=lambda fun,span,y,**kw:c.startup_limited_solver(original_solver,fun,span,y,counter=counter,first_step_factor=factor,**kw)
                sol=w.trajectory(z0,masses,rtol=rtol,max_step_m=20.)
                end=sol.sol(0.);ends.append(end)
                cases[f'{rtol:g}_{factor:g}']=dict(rtol=rtol,first_step_factor=factor,evaluations=sol.nfev,
                     **counter,ground_mass_kg=np.exp(end[:len(masses)]).tolist(),ground_age_s=end[len(masses):].tolist())
            changes[f'{rtol:g}']=dict(maximum_ground_mass_relative_difference=float(np.max(np.abs(np.exp(ends[0][:len(masses)]-ends[1][:len(masses)])-1))),
                                    maximum_ground_age_relative_difference=float(np.max(np.abs(ends[0][len(masses):]/ends[1][len(masses):]-1))),
                                    maximum_ground_age_difference_s=float(np.max(np.abs(ends[0][len(masses):]-ends[1][len(masses):]))))
    finally:
        w.geometry,w.solve_ivp=original_geom,original_solver
    if identity!={p.name:archive.sha(p) for p in bound} or source!=w.source_identity():
        raise RuntimeError('startup probe identity changed')
    failures=[]
    for filename in ['welss_compact_capacitance_initial_failure.log','welss_compact_capacitance_bound_failure.log']:
        # The /tmp originals were copied to the project-space drive on 24 October (results/tmp_relocation.json).
        path=HERE.parents[1]/'research/runs/snow_tmp_copy/lunar_snow_earth_population'/filename
        if path.exists():
            failures.append(dict(file=filename,sha256=archive.sha(path),status='Rejected attempt; invalid internal trial mass; no output product admitted'))
    product=dict(schema='terluna.atmosphere.welss-capacitance-startup/1',producer=identity,frozen_source=source,
                 source_height_m=z0,seed_masses_kg=masses.tolist(),cases=cases,halved_first_step_changes=changes,
                 failed_attempts=failures,scope='One formerly failing source height at two ODE tolerances; endpoint startup sensitivity, not full intervention convergence.')
    (HERE/'results/welss_capacitance_startup.json').write_text(json.dumps(product,indent=2,allow_nan=False)+'\n')
    print(json.dumps(changes,indent=2))


if __name__=='__main__':
    main()
