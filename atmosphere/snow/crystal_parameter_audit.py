"""Conditional tests of crystal radius conventions and step-energy compatibility.

The geometry mappings preserve either maximum span, hexagonal basal area, or
apothem. None reproduces the three-dimensional field of a hexagonal sheath.
The step-energy comparison asks whether two published inferences can share one
classical equilibrium parameter; it does not select either as a growth law.
"""
from __future__ import annotations
import copy
import json
import math
from pathlib import Path
import numpy as np
import scipy
from shared.constants import BOLTZMANN
from atmosphere.snow import facet_attachment, facet_growth, flakes
from atmosphere.snow.spatial_facet_growth import libbrecht_parameters
from atmosphere.snow.surface_transfer_assessment import observed_state, calculate, identity

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULT = HERE/'results/crystal_parameter_audit.json'
PRISM_MOLECULAR_AREA = 8.30e-20  # Miyamoto 2022, equation (4), m^2/molecule.
MAPPINGS = dict(circumscribed=1., equal_basal_area=math.sqrt(3*math.sqrt(3)/(2*math.pi)),
                inscribed=math.sqrt(3)/2)


def mapped_state(state, convention):
    """Map the radius and its target derivative together; retain measured flow."""
    if convention not in MAPPINGS:
        raise ValueError('Declare circumscribed, equal_basal_area, or inscribed')
    factor = MAPPINGS[convention]
    out = copy.deepcopy(state)
    out['radius_m'] *= factor
    out['fitted_rates_a_c_mass'][0] *= factor
    out['radius_convention'] = convention
    out['radius_factor'] = factor
    return out


def nucleation_barrier(temperature_K, ledge_energy_J_m, molecular_area_m2=PRISM_MOLECULAR_AREA, *, polynuclear=True):
    """Classical circular-nucleus sigma0; fractional, not percent.

    G(r)=2*pi*r*kappa-pi*r^2*Delta_mu/Omega. Its maximum gives
    pi*kappa^2*Omega/(kT*Delta_mu). Polynuclear normal growth introduces
    the factor 1/3 printed by Libbrecht and Rickerby (2011, p. 4).
    """
    if any(not math.isfinite(x) or x <= 0 for x in
           (temperature_K, ledge_energy_J_m, molecular_area_m2)):
        raise ValueError('Positive finite temperature, energy, and molecular area required')
    return math.pi*ledge_energy_J_m**2*molecular_area_m2/((3 if polynuclear else 1)*(BOLTZMANN*temperature_K)**2)


def energy_comparison():
    rows = []
    for tc in (-25., -20., -15.3, -10., -5.3, -2.6):
        T = flakes.T0+tc
        broad = libbrecht_parameters(tc, narrow=False)['prism_critical']
        implied = math.sqrt(3*broad*(BOLTZMANN*T)**2/(math.pi*PRISM_MOLECULAR_AREA))
        for kappa in (1.5e-10, 8.5e-10):
            sigma0 = nucleation_barrier(T, kappa)
            rows.append(dict(temperature_C=tc, broad_prism_sigma0=broad,
                broad_implied_ledge_energy_J_m=implied, microscopy_ledge_energy_J_m=kappa,
                microscopy_energy_to_broad_ratio=kappa/implied,
                equivalent_polynuclear_sigma0=sigma0,
                equivalent_barrier_to_broad_ratio=sigma0/broad,
                temperature_matched_reported_peak=tc == -15.3 and kappa == 8.5e-10))
    return dict(cases=rows, molecular_area_m2=PRISM_MOLECULAR_AREA,
        formula='sigma0 = pi*kappa^2*Omega/(3*k_B^2*T^2)',
        limits=[
            'Compare broad-facet parameters; the narrow-facet SDAK barrier is not an equilibrium step-energy estimate.',
            'The 2011 experiment measured basal growth. Only its classical relation is used here with the separately sourced prism molecular area.',
            'Microscopy kappa endpoints span temperatures. Only the reported -15.3 C peak is explicitly temperature matched here.',
            'Different step sources, substrates, pressure, surface states and inferred local humidity remain possible reasons for incompatibility.',
            'Mononuclear instead of polynuclear convention multiplies the equivalent barrier by three and does not remove this discrepancy.',
            'No author-confirmed error or replacement parameter is inferred; neither source is silently rescaled.'])


def run():
    library = facet_attachment.AttachmentLibrary(facet_growth.BUILD)
    states = [observed_state(library, age) for age in (4.5, 8., 27.)]
    cases = []
    for state in states:
        for convention in MAPPINGS:
            mapped = mapped_state(state, convention)
            for xs in (0., 14e-9, 5e-6, 10e-6):
                for delta in (0., 1e-9):
                    row = calculate(mapped, .1, xs, 1., n=32, capillary_length=delta)
                    row.update(radius_convention=convention, radius_factor=mapped['radius_factor'])
                    cases.append(row)
            print(json.dumps(dict(age=state['age_min'], convention=convention,
                completed=len(cases), failures=sum(not x['complete'] for x in cases))), flush=True)
    refinement = []
    for convention in MAPPINGS:
        for xs in (0., 5e-6):
            for delta in (0., 1e-9):
                row = calculate(mapped_state(states[1], convention), .1, xs, 1., n=64, capillary_length=delta)
                row.update(radius_convention=convention, radius_factor=MAPPINGS[convention])
                refinement.append(row)
    rows = cases+refinement
    return dict(schema='terluna.atmosphere.crystal-parameter-audit/1',
        evidence='conditional geometry sensitivity and cross-source parameter compatibility; no Earth or lunar growth validation',
        numpy_version=np.__version__, scipy_version=scipy.__version__,
        inputs=[identity(p) for p in [Path(__file__).resolve(), HERE/'crystal_parameter_sources.json',
            HERE/'surface_transfer_assessment.py', HERE/'surface_transfer.py', HERE/'spatial_crystal_growth.py',
            HERE/'spatial_diffusion.py', HERE/'spatial_facet_growth.py', HERE/'facet_attachment.py',
            HERE/'facet_growth.py', HERE/'takahashi_growth.py', HERE/'crystal_transport_replay.py',
            ROOT/'shared/constants.json']], native_build_identity=library.identity,
        states=states, radius_factors=MAPPINGS, cases=cases, refinement=refinement,
        step_energy=energy_comparison(), completed_cases=sum(r['complete'] for r in rows),
        failed_cases=sum(not r['complete'] for r in rows), limits=[
            'Cylinder mappings preserve different geometric quantities. The equal-area radius preserves a regular hexagon basal area and its area derivative, not its exterior diffusion field.',
            'Rim fraction 0.1 and cavity depth half the half-length remain supplied. No radius or migration distance is fitted.',
            'Measured maximum-span flow and ventilation are held fixed across mappings. Gas properties are independent of the mapped radius.',
            'The target radial derivative uses the same radius conversion as the geometry; axial and mass targets are unchanged.',
            'Half-width, half-length and mass histories are separate population fits, not a tracked specimen or a supplied joint-error distribution.',
            'Mean facet advance is not a prediction of maximum span; separate three-dimensional sheaths remain unrepresented.'])


if __name__ == '__main__':
    result = run()
    RESULT.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('completed_cases', 'failed_cases')}, indent=2))
