"""Positive limited finite-volume vapour growth on the existing mass pivots.

Number density is reconstructed linearly in mass cells. Slopes are limited on
true cell centroids, which differ from pivots on nonuniform grids. Face values
stay between neighboring cell-average densities. This limits reconstructed
fluxes; particle populations are never clipped or renormalized.

The discrete deposition demand is sum (m_next-m_i)*F_i, including the absorbing
upper ghost. A common vapour-supply factor caps those fluxes. Thus number and
the existing discrete first-mass moment have a closed ledger. This is explicitly
a different finite-grid growth quadrature from sum n_i*G(m_i); both approximate
the same continuous growth law as the grid is refined. Collisions and all
physical coefficients are outside this module.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

import numpy as np

from atmosphere.snow import aggregation as ag, flakes as sf, sectional, growth_transport as gt
from atmosphere.snow.habits import GEOMETRIES
from shared.provenance import constants_used

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'atmosphere/snow/results/limited_growth.json'
SCHEMA = 'terluna.atmosphere.snow-limited-growth/1'
FILES = ('atmosphere/snow/limited_growth.py',) + gt.FILES


def prepare(masses):
    """Geometric boundary ghosts, midpoint faces and true cell centroids."""
    masses = sectional._grid(masses).copy()
    ghost = masses[-1]**2/masses[-2]
    next_mass = np.r_[masses[1:], ghost]
    faces = .5*(masses+next_mass)
    lower = .5*(masses[0]+masses[0]**2/masses[1])
    edges = np.r_[lower, faces]
    widths = np.diff(edges)
    centers = .5*(edges[:-1]+edges[1:])
    return dict(masses=masses, face_masses=faces, upper_ghost_mass=ghost,
                spacing=next_mass-masses, edges=edges, widths=widths,
                centroids=centers, centroid_spacing=np.diff(centers))


class LimitedGrowth:
    """Reusable face-law evaluation for one fixed atmosphere layer.

    evaluate(number, production) returns number transfers across outflow faces,
    the actual discrete deposited mass per volume/time, and an outgoing-flux CFL
    limit. Combine adjacent transfers to form the number tendency. The last
    transfer exits at upper_ghost_mass; it is never counted as snowfall.
    """
    def __init__(self, masses, face_growth, *, cfl=.3):
        self.grid = prepare(masses)
        for key, value in self.grid.items():
            setattr(self, key, value)
        self.face_growth = np.asarray(face_growth, dtype=float).copy()
        if (self.face_growth.shape != self.masses.shape or not np.isfinite(self.face_growth).all()
                or (self.face_growth<0).any()):
            raise ValueError('face growth must be finite, nonnegative and match masses')
        if not math.isfinite(cfl) or not 0<cfl<=1:
            raise ValueError('CFL must lie in (0,1]')
        self.cfl = cfl

    def face_density(self, number):
        n = np.asarray(number, dtype=float)
        if n.shape != self.masses.shape or not np.isfinite(n).all() or (n<0).any():
            raise ValueError('population must be finite, nonnegative and match masses')
        density = n/self.widths
        slope = np.zeros_like(density)
        if density.size>2:
            left = (density[1:-1]-density[:-2])/self.centroid_spacing[:-1]
            right = (density[2:]-density[1:-1])/self.centroid_spacing[1:]
            h_left, h_right = self.centroid_spacing[:-1], self.centroid_spacing[1:]
            centered = (h_right*left+h_left*right)/(h_left+h_right)
            same = (left*right>0)
            slope[1:-1] = np.where(same, np.sign(centered)*np.minimum(np.abs(centered),
                                                    np.minimum(2*np.abs(left), 2*np.abs(right))), 0.)
        face = density+slope*(self.face_masses-self.centroids)
        # Bound only reconstructed face values, never the population state.
        face[:-1] = np.minimum(np.maximum(face[:-1], np.minimum(density[:-1],density[1:])),
                               np.maximum(density[:-1],density[1:]))
        return face

    def evaluate(self, number, production):
        if not math.isfinite(production) or production<0:
            raise ValueError('vapour production must be finite and nonnegative')
        n = np.asarray(number, dtype=float)
        transfers = self.face_growth*self.face_density(n)
        demand = float(np.dot(self.spacing, transfers))
        if demand>0:
            transfers *= min(1., production/demand)
        deposited = float(np.dot(self.spacing, transfers))
        outgoing = transfers>0
        limit = float(self.cfl*np.min(n[outgoing]/transfers[outgoing])) if outgoing.any() else math.inf
        return transfers, deposited, limit


def integrate(operator, numbers, times, production, *, max_step_s=30., deadline=None):
    """Growth-only box using the unchanged positive SSPRK2 integrator."""
    masses, initial, times = gt._validate(operator.masses, numbers, times, production)
    if not math.isfinite(max_step_s) or max_step_s<=0:
        raise ValueError('maximum step must be positive')
    def evaluate(n):
        gt._deadline(deadline)
        transfers, deposited, limit = operator.evaluate(n, production)
        dn = -transfers
        dn[1:] += transfers[:-1]
        return dn, np.array([deposited, transfers[-1]*operator.upper_ghost_mass, transfers[-1]]), limit
    start=time.monotonic()
    state=initial.copy(); snapshots=[state.copy()]; ledgers=[np.zeros(3)]; ledger=np.zeros(3)
    t=0.; steps=0; rejected=0; smallest=math.inf; largest=0.
    for end in times[1:]:
        while t<end:
            state,dt,rates,_,retry=sectional._advance(state,evaluate,min(max_step_s,end-t),operator.cfl)
            ledger+=rates*dt; t+=dt; steps+=1; rejected+=retry
            smallest,largest=min(smallest,dt),max(largest,dt)
        snapshots.append(state.copy()); ledgers.append(ledger.copy())
    return dict(masses=masses,numbers=np.asarray(snapshots),times=times,
                deposited_mass=np.asarray(ledgers)[:,0],upper_exit_mass=np.asarray(ledgers)[:,1],
                upper_exit_number=np.asarray(ledgers)[:,2],max_step_s=max_step_s,
                diagnostics=dict(steps=steps,rejected_steps=rejected,minimum_step_s=smallest,
                                 maximum_step_s=largest,wall_seconds=time.monotonic()-start))


def identity(inputs):
    return dict(files={p:gt.sha(ROOT/p) for p in FILES},inputs={p:gt.sha(ROOT/p) for p in inputs},
                constants=constants_used(FILES))


def analyse(*, budget_s=180.):
    start=time.monotonic(); deadline=start+budget_s
    baseline_path='atmosphere/snow/results/growth_transport.json'
    baseline=json.loads((ROOT/baseline_path).read_text())
    for kind in ('files','inputs'):
        for path,expected in baseline['producer'][kind].items():
            if gt.sha(ROOT/path)!=expected:
                raise RuntimeError('Stale growth reference: '+path)
    world,binding=sf.moon_column(); settings=baseline['settings']; geometry=GEOMETRIES[settings['geometry']]
    state=world.state(settings['height_m']); production=settings['vapour_production_kg_m3_s']
    def growth(m):
        d=geometry.size(m)
        velocity,re=ag.speeds(d,m,world.gravity,state,geometry.area_ratio)
        return np.maximum(ag.deposition(d,re,state,geometry.capacitance_ratio),0.)
    inputs=dict(binding['inputs']);inputs[baseline_path]=gt.sha(ROOT/baseline_path)
    before=identity(inputs)
    if before['inputs']!=inputs:
        raise RuntimeError('Column input identity changed')
    runs=[]
    for bins in (96,144,216):
        reference=next(r for r in baseline['runs'] if r['id']==f'm{bins}_dt30')
        masses=np.asarray(reference['mass_pivots_kg']); initial=np.asarray(reference['initial_numbers_m3'])
        times=np.array([r['time_s'] for r in reference['snapshots']])
        grid=prepare(masses); operator=LimitedGrowth(masses,growth(grid['face_masses']))
        for step in ((30.,15.,7.5) if bins==144 else (30.,)):
            result=integrate(operator,initial,times,production,max_step_s=step,deadline=deadline)
            rows=[]; initial_mass=float(np.dot(initial,masses))
            for k,(t,n) in enumerate(zip(times,result['numbers'])):
                old=reference['snapshots'][k]
                row=gt.moments(masses,n,geometry)
                support=old['once_projected_support_maximum_m']
                above=geometry.size(masses)>support
                dep,exited=float(result['deposited_mass'][k]),float(result['upper_exit_mass'][k])
                row.update(time_s=float(t),characteristic=old['characteristic'],
                           once_projected_characteristic=old['once_projected_characteristic'],
                           once_projected_support_maximum_m=support,
                           mass_share_beyond_once_projected_support=float(np.dot(masses[above],n[above])/row['mass_kg_m3']),
                           excess_variance_over_once_projection_kg2=row['variance_mass_kg2']-old['once_projected_characteristic']['variance_mass_kg2'],
                           water_ledger=dict(initial=initial_mass,deposited=dep,remaining=row['mass_kg_m3'],upper_exit=exited,
                                             residual=initial_mass+dep-exited-row['mass_kg_m3']),
                           number_residual=float(initial.sum()-n.sum()-result['upper_exit_number'][k]))
                rows.append(row)
            previous=reference['snapshots'][-1]; current=rows[-1]; exact=current['once_projected_characteristic']
            old_excess=previous['relative_mass_variance']-exact['relative_mass_variance']
            new_excess=current['relative_mass_variance']-exact['relative_mass_variance']
            comparison=dict(excess_relative_variance_ratio=abs(new_excess/old_excess),
                            signed_excess_relative_variance=new_excess,
                            mass_above_5mm_ratio=current['mass_share_above_diameter']['0.005']/previous['mass_share_above_diameter']['0.005'],
                            reference_id=reference['id'])
            runs.append(dict(id=f'm{bins}_dt{step:g}',mass_bins=bins,max_step_s=step,diagnostics=result['diagnostics'],
                             snapshots=rows,comparison_with_upwind=comparison,mass_pivots_kg=masses.tolist(),
                             final_numbers_m3=result['numbers'][-1].tolist()))
            print(f'limited m{bins} dt{step:g}: variance={current["relative_mass_variance"]:.6g}, '
                  f'mass>5mm={current["mass_share_above_diameter"]["0.005"]:.6g}, '
                  f'{result["diagnostics"]["wall_seconds"]:.2f}s',flush=True)
    byid={r['id']:r for r in runs}
    coarse,half,fine=[byid[id]['snapshots'][-1] for id in ('m144_dt30','m144_dt15','m144_dt7.5')]
    time_controls={}
    for key in ('relative_mass_variance','mass_share_beyond_once_projected_support'):
        spatial=abs(fine[key]-(fine['once_projected_characteristic'][key] if key=='relative_mass_variance' else 0.))
        last=abs(half[key]-fine[key])
        time_controls[key]=dict(coarse_to_half=abs(coarse[key]-half[key]),half_to_quarter=last,
                               residual_spatial_difference=spatial,passed=last<.05*spatial if spatial else last==0)
    gains=[r['comparison_with_upwind'] for r in runs if r['max_step_s']==30.]
    grid_errors=[abs(byid[f'm{b}_dt30']['snapshots'][-1]['relative_mass_variance']
                     -byid[f'm{b}_dt30']['snapshots'][-1]['characteristic']['relative_mass_variance'])
                 for b in (96,144,216)]
    acceptance=dict(grid_error_reduces=grid_errors[2]<grid_errors[1]<grid_errors[0],
        improves_variance_and_fixed_tail_on_every_grid=all(
        c['excess_relative_variance_ratio']<.5 and c['mass_above_5mm_ratio']<.5 for c in gains),
        time_controls_pass=all(c['passed'] for c in time_controls.values()))
    if identity(inputs)!=before:
        raise RuntimeError('Producer/input identity changed during limited-growth comparison')
    return dict(schema=SCHEMA,producer=before,settings=settings,runs=runs,time_controls=time_controls,
                acceptance=acceptance,grid_relative_variance_errors=grid_errors,wall_seconds=time.monotonic()-start,
                evidence='Limited finite-volume growth comparison against independently integrated continuous masses '
                         'and the frozen first-order growth product. No change to physical rate coefficients; '
                         'the finite-grid growth quadrature is changed explicitly. Closed box, no coupled-column '
                         'convergence or empirical validation.',
                reading_rule='The deposited-water ledger uses sum of pivot mass gaps times limited number flux, '
                             'including the ghost exit. Passing the improvement gates is not vanishing transport '
                             'error: retain mass/time controls and remaining reference discrepancies.')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--budget-s',type=float,default=180.)
    args=parser.parse_args();product=analyse(budget_s=args.budget_s)
    OUT.write_text(json.dumps(product,indent=1,default=float)+'\n')
    print(f'Wrote {OUT.relative_to(ROOT)} in {product["wall_seconds"]:.2f}s; gates={product["acceptance"]}',flush=True)


if __name__=='__main__': main()
