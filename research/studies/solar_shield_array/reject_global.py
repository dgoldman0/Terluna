"""Preserve the rejected rigid schedule without claiming its aborted validation."""
import json
import numpy as np

from shared import constants as K
from shared.provenance import constants_used
from protection.dynamics.active_global import Traffic,Layout,pole_collisions
from .active_global_run import sources,budget
from .fleet_run import environment,digest
from .cycling_search import ROOT,HERE
from .cycling_analysis import compact_series_json


def main():
    env=environment(30);f=Traffic(env)
    log=ROOT/'research/runs/solar_shield_array/active/global_run.log'
    partial=[]
    for line in log.read_text().splitlines():
        if line.startswith('{'):
            row=json.loads(line)
            if row.get('stage')=='coverage':partial.append({k:v for k,v in row.items() if k!='stage'})
    raw=ROOT/'research/runs/solar_shield_array/active/global_validate_geometry.json'
    geometry=json.loads(raw.read_text())
    src={**sources(),'research/studies/solar_shield_array/active_global_validate.py':digest(HERE/'active_global_validate.py')}
    if geometry['producer']['source_hashes']!=src:raise ValueError('Geometry source changed')
    src['research/studies/solar_shield_array/reject_global.py']=digest(__file__)
    costs=budget(f,np.linspace(0,30*K.JULIAN_DAY,9),plane_stride=32,phase_count=64)
    compact=Traffic(env,layout=Layout(radial_step=1000.))
    held=json.loads((HERE/'results/holding.json').read_text())['aperture_quadrature']['variable_distance_finer']
    periods=2*np.pi/np.abs(f.omega)
    result=dict(schema='terluna.research.active-global-rejected/1',
        producer=dict(source_hashes=src,constants=constants_used(src),
            input_products={'holding.json':digest(HERE/'results/holding.json')}),
        status='REJECTED_HIGH_RECURRING_COST',inventory=f.inventory,
        optical_base_mass_kg=f.inventory*f.layout.side**2*env.sigma,
        preliminary_budget=costs,separation=f.separation_certificate(),
        nominal_period_days=[float(periods.min()/K.JULIAN_DAY),float(periods.max()/K.JULIAN_DAY)],
        front_passages_per_day=float(np.sum(f.count/periods)*K.JULIAN_DAY),
        nominal_handover_interval_s=[float(np.min(periods/f.count)),float(np.max(periods/f.count))],
        interrupted_global_coverage=dict(completed=False,dates_completed=len(partial),
            last_day=partial[-1]['day'],solar_directions=16,area_receivers=8192,rim_receivers=1024,
            ray_tests=len(partial)*16*(8192+1024),series=partial,
            log_path=str(log.relative_to(ROOT)),log_sha256=digest(log)),
        completed_local_geometry={k:v for k,v in geometry.items() if k!='producer'},
        completed_local_geometry_raw=dict(path=str(raw.relative_to(ROOT)),sha256=digest(raw)),
        compact_rejected=dict(inventory=compact.inventory,return_pole=pole_collisions(compact,0)),
        held_comparison=dict(mean_power_TW=held['mean_power_TW'],total_mass_kg=held['total_mass_kg'],
            base_optical_mass_kg=held['base_optical_mass_kg'],
            mean_residual_acceleration_m_s2=held['area_mean_residual_mm_s2']/1000),
        acceptance=dict(global_four_lunar_radius_coverage_demonstrated=False,
            low_energy_active_gap_repair_demonstrated=False,delivered_electrical_power_demonstrated_W=None),
        evidence='Rejected prescribed-circle comparison. Preliminary force quadrature and completed optical tests are retained; the 30-day coverage/census run was interrupted after the user challenged its cost.',
        unperformed=['Full 30-day spatial survey','All-integer force census','Force quadrature convergence',
            'Planned eight-day rigid-target feedback replay','Bounded-slew optical attitude realization'],
        reason='Lower acceleration per kilogram is outweighed by large circulating inventory and continuous enforcement of Sun-following circles. This does not measure minimum-intervention gap repair on natural trajectories.',
        redirection='Propagate the same starting fleet gravitationally, cancel only its photon force, and measure the coverage and collision failures before assigning gap-repair costs.')
    (HERE/'results/active_global_rejected.json').write_text(compact_series_json(result))
    print(json.dumps(dict(saved='active_global_rejected.json',power_bounds_TW=(np.array(costs['mean_power_bounds_W'])/1e12).tolist())),flush=True)


if __name__=='__main__':main()
