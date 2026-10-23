"""Identifiability of spiral-step kinetics under unmeasured vapour depletion.

This tests a candidate explanation, not a correction to either experiment.
A supplied chamber humidity ratio can change both inferred kinetic coefficient
and step energy while preserving the measured speed and spacing. The actual
three-dimensional chamber field and other growing crystals are not supplied.
"""
from __future__ import annotations
import json
import math
from pathlib import Path
from shared.constants import BOLTZMANN
from atmosphere.snow import facet_attachment, facet_growth, flakes
from atmosphere.snow.crystal_parameter_audit import energy_comparison, PRISM_MOLECULAR_AREA
from atmosphere.snow.surface_transfer import prism_spiral_rate
from atmosphere.snow.surface_transfer_assessment import identity

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
RESULT=HERE/'results/step_humidity_audit.json'


def infer_under_depletion(apparent_beta, apparent_energy, apparent_sigma, fraction):
    """Exact algebraic inversion for supplied sigma_true=fraction*sigma_app.

    V=beta*sigma and L=19*Omega*kappa/(kT*ln(1+sigma)). The small-sigma
    limit is beta_true=beta_app/f and kappa_true=f*kappa_app.
    """
    if any(not math.isfinite(x) or x<=0 for x in (apparent_beta,apparent_energy,apparent_sigma,fraction)) or fraction>1:
        raise ValueError('Positive finite parameters and humidity fraction in (0,1] required')
    actual=fraction*apparent_sigma
    return dict(actual_surface_sigma=actual,
        inferred_beta_m_s=apparent_beta/fraction,
        inferred_energy_J_m=apparent_energy*math.log1p(actual)/math.log1p(apparent_sigma))


def mean_attachment_for_depletion(kinetic_length_m,chamber_length_m,fraction):
    """Ideal planar chamber: f=alpha_diff/(alpha_bar+alpha_diff).

    This is a supplied geometry sensitivity, not the measured chamber.
    Values above one reject the planar configuration under alpha<=1.
    """
    if any(not math.isfinite(x) or x<=0 for x in (kinetic_length_m,chamber_length_m,fraction)) or fraction>1:
        raise ValueError('Positive lengths and humidity fraction in (0,1] required')
    return kinetic_length_m/chamber_length_m*(1/fraction-1)


def diffusion_path_for_growth(kinetic_length, kinetic_speed, reservoir_sigma, normal_rate, fraction):
    """Planar path required if the observed spiral face is the only sink."""
    if any(not math.isfinite(x) or x <= 0 for x in
           (kinetic_length, kinetic_speed, reservoir_sigma, normal_rate, fraction)) or fraction > 1:
        raise ValueError('Positive transport and growth scales, fraction in (0,1] required')
    return kinetic_length*kinetic_speed*reservoir_sigma*(1-fraction)/normal_rate


def run():
    library=facet_attachment.AttachmentLibrary(facet_growth.BUILD)
    matched=next(r for r in energy_comparison()['cases'] if r['temperature_matched_reported_peak'])
    T=flakes.T0+matched['temperature_C']
    app_energy=matched['microscopy_ledge_energy_J_m']
    target_energy=matched['broad_implied_ledge_energy_J_m']
    fraction=target_energy/app_energy
    air=library.evaluate(.01,T,101000.,20e-6,20e-6,10.)
    vth=air['vapor_molecular_speed_m_s']/4
    rhoeq=air['ice_saturation_pressure_Pa']/(air['vapor_gas_constant_J_kg_K']*T)
    vkin=rhoeq/facet_growth.ICE_DENSITY*vth
    X0=air['diffusivity_m2_s']/vth
    rows=[]
    for sigma in (.001,.005,.01,.02,.05):
        app=prism_spiral_rate(T,sigma,app_energy)
        corrected=infer_under_depletion(app['step_coefficient_m_s'],app_energy,sigma,fraction)
        spacing=19*PRISM_MOLECULAR_AREA*corrected['inferred_energy_J_m']/(BOLTZMANN*T*math.log1p(corrected['actual_surface_sigma']))
        velocity=corrected['inferred_beta_m_s']*corrected['actual_surface_sigma']
        # Libbrecht 2015 reanalysis equation (3), one collecting terrace and
        # alpha near a step set to one. This is an effective collection length.
        collection=.39e-9*corrected['inferred_beta_m_s']/vkin
        rows.append(dict(apparent_sigma=sigma,apparent=app,**corrected,
            reconstituted_spacing_m=spacing,reconstituted_step_velocity_m_s=velocity,
            inferred_energy_to_broad_ratio=corrected['inferred_energy_J_m']/target_energy,
            effective_one_terrace_collection_length_m=collection,
            spiral_mean_attachment=app['normal_rate_m_s']/(vkin*corrected['actual_surface_sigma']),
            all_spiral_plane_path_needed_m=diffusion_path_for_growth(X0,vkin,sigma,app['normal_rate_m_s'],fraction),
            spacing_relative_residual=abs(spacing/app['spacing_m']-1),
            velocity_relative_residual=abs(velocity/app['lateral_step_velocity_m_s']-1)))
    planes=[]
    for H in (.2e-3,2e-3,16e-3):
        alpha=mean_attachment_for_depletion(X0,H,fraction)
        mixtures=[]
        for row in rows:
            spiral=row['spiral_mean_attachment']
            rough=(alpha-spiral)/(1-spiral)
            mixtures.append(dict(apparent_sigma=row['apparent_sigma'],spiral_mean_attachment=spiral,
                unit_attachment_area_fraction=rough,mixture_feasible=0<=rough<=1))
        planes.append(dict(supplied_chamber_length_m=H,mean_attachment=alpha,
            compatible_with_alpha_at_most_one=alpha<=1,spiral_and_unit_attachment_mixtures=mixtures))
    return dict(schema='terluna.atmosphere.step-humidity-audit/1',
        evidence='conditional inverse-problem consistency; no measured chamber correction or adopted kinetic law',
        inputs=[identity(p) for p in [Path(__file__).resolve(),HERE/'crystal_morphology_sources.json',
            HERE/'crystal_parameter_audit.py',HERE/'spatial_facet_growth.py',HERE/'surface_transfer.py',
            HERE/'facet_attachment.py',HERE/'facet_growth.py',ROOT/'shared/constants.json']],
        native_build_identity=library.identity,temperature_K=T,air=air,kinetic_speed_m_s=vkin,
        kinetic_length_m=X0,matched_energy_comparison=matched,required_low_sigma_humidity_fraction=fraction,
        cases=rows,planar_chamber_sensitivity=planes,limits=[
            'The humidity ratio is chosen to reconcile two energy inferences, not independently measured or calibrated against a chamber field.',
            'At nonzero supersaturation, exactly preserving the measured spacing requires the logarithmic correction shown. A single low-sigma energy differs slightly from that exact reconstruction.',
            'Measured lateral velocity and spacing, and hence normal growth at that measurement, are unchanged by the inversion.',
            'The effective migration length assumes one collecting terrace and unit attachment near a step. It is not a measurement of prism surface diffusion.',
            'The required path and two-surface mixture test whether the measured spiral face alone could sustain the hypothesized depletion. Unit-attachment patches are a supplied upper-efficiency sink, not observed rough-area measurements.',
            'The ideal planar chamber lengths are supplied alternatives; the 2022 chamber geometry, active ice area and local humidity remain unresolved.',
            'The 2015 reanalysis concerns basal data and estimates its own migration length only within roughly a factor of three. It motivates a test, not a correction transferred to the 2022 prism experiment.',
            'No corrected parameters are inserted into the Earth or lunar growth calculations.'])


if __name__=='__main__':
    result=run()
    RESULT.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(fraction=result['required_low_sigma_humidity_fraction'],
        collection_length_m=result['cases'][0]['effective_one_terrace_collection_length_m'],
        planes=result['planar_chamber_sensitivity']),indent=2))
