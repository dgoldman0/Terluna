#!/usr/bin/env python3
"""Reproduce this bounded research pass; no automatic whole-world PASS status.

Run from repository root: python research/run_environment_screens.py
Optional original ZIP: --protection-archive /path/to/Lunar_Protection_Model.zip
Outputs default to research/results/environment_screens (small JSON/CSV tables).
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,sys,platform,zipfile
from dataclasses import asdict,replace
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from atmosphere.thermal_column import ColumnConfig,solve_column,R,AREA
from atmosphere.spectral_interface import audit_optical_files,inventory_audit,historical_solar_bands,gap_energy_bounds,deposited_heat
from biosphere.long_night import carbon_trace,periodic_storage_requirement,simulate_store,oxygen_periodic,PERIOD_DAYS
from climate.cycle import ClimateConfig,solve_climate


def plain(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x).__name__)

def write_json(path,obj):
    path.write_text(json.dumps(obj,indent=2,default=plain,allow_nan=False)+'\n')

def write_csv(path,rows):
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',newline='') as f:
        out=csv.DictWriter(f,fieldnames=keys,lineterminator="\n");out.writeheader()
        for row in rows:
            out.writerow({k:json.dumps(v,default=plain) if isinstance(v,(dict,list)) else v for k,v in row.items()})


def run(out,archive=None):
    out.mkdir(parents=True,exist_ok=True)
    import scipy
    configurations=[]
    for tb in [150.,180.,210.]:
        for pb in [.001,.1,10.]:
            for shape in ['low','middle','high']:
                configurations.append(ColumnConfig(lower_temperature_k=tb,lower_pressure_pa=pb,heating_shape=shape))
    for multiplier in [.5,2.]:configurations.append(ColumnConfig(conductivity_multiplier=multiplier))
    heats=[0.,1e-8,1e-7,3e-7,1e-6,3e-6,1e-5,3e-5]
    rows=[];profiles=[]
    for config in configurations:
        previous=None
        for q in heats:
            try:
                row,pf,sol=solve_column(q,config,previous=previous)
                row['numerical_status']='converged';previous=sol
                if config==ColumnConfig():profiles.append({'heat_W_m2':q,'profile':{k:v[::10] for k,v in pf.items()}})
            except (ValueError,RuntimeError,FloatingPointError) as exc:
                row=dict(lower_temperature_k=config.lower_temperature_k,lower_pressure_pa=config.lower_pressure_pa,
                    heating_shape=config.heating_shape,conductivity_multiplier=config.conductivity_multiplier,
                    deposited_heat_w_m2=q,numerical_status='failed',error=str(exc))
            rows.append(row)
    write_csv(out/'molecular_columns.csv',rows)
    write_json(out/'molecular_profiles.json',profiles)

    carbon=[]
    for assimilation in [1.,2.,3.]:
        for night_temp in [268.,278.,288.]:
            for demand in [.1,.25,.5,1.]:
                net,dt,_=carbon_trace(mean_light_assimilation=assimilation,night_temp_k=night_temp,night_demand_fraction=demand)
                r=periodic_storage_requirement(net,dt)
                r.update(day_usable_assimilation=assimilation,night_temperature_K=night_temp,
                         night_demand_fraction=demand,units='day-reference-maintenance equivalents')
                for cap in [1.,2.,4.,8.,16.]:
                    s=simulate_store(net,dt,cap,cycles=4)
                    carbon.append({**r,'capacity':cap,'last_cycle_unmet':s['cycles'][-1]['unmet'],
                        'last_cycle_spill':s['cycles'][-1]['spill'],'ledger_residual':s['ledger_residual']})
    write_csv(out/'carbon_budget.csv',carbon)
    oxygen=[]
    for resp in [.25,.5,1.,2.]:
        for k in [.01,.05,.1,.2,.5]:
            r=oxygen_periodic(respiration_g_m3_day=resp,mean_day_production_g_m3_day=3*resp,exchange_per_day=k)
            r.pop('trace_g_m3');oxygen.append(r)
    write_csv(out/'oxygen_budget.csv',oxygen)

    climates=[];coupled=[]
    for geography in ['dry','concentrated','distributed']:
        for participation in [.02,.1,1.]:
            for olr in [180.,200.,220.]:
                for transport in [.3,1.]:
                    cfg=ClimateConfig(layout=geography,participating_atmosphere=participation,
                                      olr_at_273_w_m2=olr,heat_transport_w_m2_k=transport)
                    r,trace=solve_climate(cfg);climates.append(r)
                    if olr==200. and transport==1.:
                        for latitude_index in [0,2]:
                            cell=latitude_index*cfg.nlon+cfg.nlon//4
                            net,dt,_=carbon_trace(steps=cfg.steps,mean_light_assimilation=3.,
                                night_demand_fraction=.25,light_trace=trace['incident_W_m2'][:,cell],
                                temperature_trace=trace['temperature_K'][:,cell],light_reference=600.)
                            req=periodic_storage_requirement(net,dt)
                            coupled.append({**req,'layout':geography,'participating_atmosphere':participation,
                                'latitude_deg':float(np.degrees(np.arcsin(trace['sin_latitude'][latitude_index]))),
                                'cell_min_K':float(trace['temperature_K'][:,cell].min()),
                                'cell_max_K':float(trace['temperature_K'][:,cell].max()),
                                'light_reference_W_m2':600.,'hypothetical_trait_assumptions':True,
                                'cold_or_heat_injury_not_modelled':True})
    write_csv(out/'climate_budget.csv',climates)
    write_csv(out/'climate_carbon_interface.csv',coupled)

    audit=inventory_audit();audit['actual_inputs_inspected']=False
    existing=ROOT/'protection/sources'
    if archive:
        with zipfile.ZipFile(archive) as z:
            solar=z.read('Lunar_Protection_Model/sources/TSIS1_HSRS_stride100.csv').decode()
            titania=z.read('Lunar_Protection_Model/sources/TiO2_Siefke.yml').decode()
        audit=audit_optical_files(solar,titania)
        audit['source_archive_sha256']=hashlib.sha256(Path(archive).read_bytes()).hexdigest()
    elif (existing/'TSIS1_HSRS_stride100.csv').exists() and (existing/'TiO2_Siefke.yml').exists():
        audit=audit_optical_files((existing/'TSIS1_HSRS_stride100.csv').read_text(),(existing/'TiO2_Siefke.yml').read_text())
    historic=historical_solar_bands()
    audit['historical_band_gap_bounds']=gap_energy_bounds(historic,audit['uncovered_material_interval_nm'])
    audit['historical_band_reference']='Ribas et al. Table 4, Sun column; 1993-based composite, not present-day resolved irradiance'
    write_json(out/'spectral_audit.json',audit)
    write_json(out/'historical_solar_bands.json',[asdict(b) for b in historic])
    band_screen=[]
    # All response quantities below are explicit hypothetical scenarios, not measured filter results.
    for transmission in [1e-5,1e-4,1e-3]:
        for eta in [.03,.10,.30]:
            for radius in [2.,3.,5.]:
                supplied=[replace(b,transmission=transmission,atmospheric_absorptance=1.,
                    heat_fraction=eta,absorption_radius_R=radius,
                    provenance=b.provenance+'; imposed uniform response parameters') for b in historic]
                q,ledger=deposited_heat(supplied,(.1,118.))
                band_screen.append(dict(transmission=transmission,heat_fraction=eta,
                    absorption_radius_R=radius,atmospheric_absorptance=1.,
                    deposited_heat_W_m2=q,uniform_parameters_are_assumptions=True,
                    excludes_wavelengths_longer_than_118nm=True,
                    excludes_chemical_escape_and_radiative_cooling=True))
    write_csv(out/'band_heat_scenarios.csv',band_screen)

    # Existing shield power reproduced from its published base-case massflow and ve.
    thrust_power=125313596024069.12
    contamination=[dict(intercepted_fraction=f,thermalization_fraction=eta,
        deposited_heat_W_m2=thrust_power*f*eta/AREA,
        origin='Inherited candidate shield exhaust kinetic power; imposed intercepted/thermalized fractions')
        for f in [1e-8,1e-7,1e-6,1e-5] for eta in [.1,1.]]
    write_csv(out/'plume_heat_interface.csv',contamination)
    summary=dict(base_commit='79a60a1a51b4552ac9ae8e8bb295b567d9254995',
        atmosphere_cases=len(rows),atmosphere_numerically_converged=sum(r['numerical_status']=='converged' for r in rows),
        atmosphere_with_kinetic_flags=sum('KINETIC_ESCAPE_SENSITIVITY_REQUIRED' in r.get('domain_flags',[]) for r in rows),
        carbon_scenarios=len(carbon),oxygen_scenarios=len(oxygen),climate_scenarios=len(climates),
        climate_carbon_interfaces=len(coupled),historical_band_response_scenarios=len(band_screen),
        full_spectrum_chemistry_closure='BLOCKED_INPUTS_AND_UNIMPLEMENTED_PROCESSES',
        biological_viability='UNESTABLISHED_REQUIREMENTS_ONLY',
        climate_calibration='UNESTABLISHED_PARAMETER_SCREEN',
        integrated_habitability='UNESTABLISHED',python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
    write_json(out/'summary.json',summary)
    write_json(out/'configurations.json',dict(column_base=asdict(ColumnConfig()),column_heats=heats,
        column_configurations=[asdict(c) for c in configurations],climate_base=asdict(ClimateConfig()),
        biological_parameters_are_hypothetical=True,period_days=PERIOD_DAYS))
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=ROOT/'research/results/environment_screens')
    parser.add_argument('--protection-archive',type=Path)
    args=parser.parse_args();run(args.out,args.protection_archive)
