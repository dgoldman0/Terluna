"""Price the inherited first Earth-route fuel stock against the actual seed fleet.

This transparent seed-only addendum changes no annual flow or equipment result.
The initial fuel is bought on Earth and delivered to LEO under the inherited
external launch contract. Later source-supplied fuel remains in physical flows.
"""
from __future__ import annotations
from copy import deepcopy
import argparse
import hashlib
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def apply(product,engineering):
    result=deepcopy(product);old=engineering['original_architecture'];v=old['peak_flows']['Earth_vehicle']
    hull_per_ship=old['assets']['Earth_hull_tanks']['mass_kg']/old['engineering']['counts']['Earth_route_ships']
    fraction=engineering['transport']['depots']['seed_earth_fleet_fraction'];I=engineering['normalization']['USD2026_to_USD2041'];t=engineering['industrial']['transport']
    for case in result['cases']:
        if not case['annual'] or case.get('seed') is None:continue
        seed=case['seed']
        if 'priced_initial_route_fuel' in seed:raise ValueError('Initial route fuel already priced')
        ships=math.ceil(case['annual'][0]['asset_stocks_kg']['Earth_hull_tanks']/hull_per_ship-1e-10)
        fuel=ships*fraction*v['burns']['ascent_kg'];purchase=fuel*t['earth_propellant_USD2026_kg']*I;launch=fuel*t['LEO_USD2026_kg']*I
        seed['priced_initial_route_fuel']=dict(Earth_route_ship_count=ships,seed_fleet_fraction=fraction,
            fuel_to_LEO_kg=fuel,external_propellant_purchase_USD2041=purchase,Earth_launch_USD2041=launch,total_USD2041=purchase+launch,
            formula='Actual opening Earth-route hull count × inherited starting fleet fraction × inherited ascent propellant per ship. Earth propellant purchase and external Earth-to-LEO launch are paid once; no local/source feed is credited or charged for this Earth-provided stock.')
        seed['cost']['initial_route_fuel_USD2041']=purchase+launch;seed['cost']['total_USD2041']+=purchase+launch
        seed['remaining_unpriced_requirements']=[
            'A qualified sequence for the first source/depot deployment, including any fuel beyond the priced initial Earth-route stock, storage filling and commissioning losses.',
            'Site and source access, passenger access and preproduction staff accommodation beyond the average material and payroll screen.',
            'Civil and service-facility details not established by the inherited per-resident, per-staff and per-bed proxies.'
        ]
    result['rules'].append('The bootstrap addendum prices the inherited initial Earth-route fuel formula against the actual opening fleet. It adds Earth procurement and launch once to seed cost, leaves annual physics and source feed unchanged, and records the remaining unpriced bootstrap requirements.')
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('source',type=Path);ap.add_argument('output',type=Path);args=ap.parse_args()
    engineering_path=HERE/'inputs/engineering.json';result=apply(json.loads(args.source.read_text()),json.loads(engineering_path.read_text()))
    result['postprocessing']=dict(path='provisioning/enclosed_industry/bootstrap.py',sha256=sha(__file__),source_product_sha256=sha(args.source),engineering_sha256=sha(engineering_path),scope='Seed initial external route-fuel procurement only; annual physics unchanged.')
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(args.output)

if __name__=='__main__':main()
