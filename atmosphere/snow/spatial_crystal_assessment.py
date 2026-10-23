"""Moving-interface numerical controls and a conditional warm Earth restart.

This own event-grid implementation tests source-inspired mechanisms. It is not
an exact reproduction of the published cellular automata, and a fitted seed
mass with unmeasured cavity shape is not a measured initial crystal profile.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import time
import numpy as np
import scipy
from atmosphere.snow import facet_attachment, facet_growth, flakes
from atmosphere.snow import takahashi_growth as lab, crystal_transport_replay as replay
from atmosphere.snow.spatial_crystal_growth import GrowingCrystal, grid
from atmosphere.snow.spatial_diffusion import Attachment
from atmosphere.snow.spatial_facet_growth import libbrecht_laws

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULT = HERE / 'results/spatial_crystal_assessment.json'


def file_id(path):
    return dict(path=str(path.relative_to(ROOT)), sha256=facet_growth.sha(path))


def profile(crystal):
    """Filled-cell envelope and cross-sectional ice areas, including partial cells."""
    rows = []
    for j in np.where(np.any(crystal.fraction > 0, axis=0))[0]:
        solid = np.where(crystal.solid[:, j])[0]
        occupied = np.where(crystal.fraction[:, j] > 0)[0]
        rows.append(dict(z_lower_m=float(crystal.zf[j]*crystal.dx),
            z_upper_m=float(crystal.zf[j+1]*crystal.dx),
            filled_outer_radius_m=float(crystal.rf[max(solid)+1]*crystal.dx) if len(solid) else 0.,
            occupied_outer_radius_m=float(crystal.rf[max(occupied)+1]*crystal.dx),
            ice_area_m2=float(math.pi*crystal.dx**2*np.dot(np.diff(crystal.rf**2), crystal.fraction[:, j]))))
    return rows


def evolve(crystal, forcing, end, interval, wall_limit):
    start = time.monotonic()
    trace = [crystal.state()]
    residual = 0.
    max_nonlinear = 0.
    error = None
    try:
        while crystal.time < end:
            target = min(end, interval*(math.floor((crystal.time+1e-10)/interval)+1))
            out = crystal.advance(until=target, **forcing(crystal.time))
            residual = max(residual, out.get('mass_residual', 0.))
            max_nonlinear = max(max_nonlinear, out.get('nonlinear_residual', 0.))
            if abs(crystal.time/interval-round(crystal.time/interval)) < 1e-9:
                trace.append(crystal.state())
            if time.monotonic()-start > wall_limit:
                raise RuntimeError('Declared per-case wall-time limit reached')
    except (ValueError, ArithmeticError, RuntimeError) as exc:
        error = str(exc)
    return dict(complete=crystal.time >= end and error is None, error=error,
        target_elapsed_s=end, trace=trace, final=crystal.state(),
        final_profile=profile(crystal), wall_elapsed_s=time.monotonic()-start,
        max_water_residual=residual, max_nonlinear_residual=max_nonlinear)


def cylinder_controls():
    # Infinite cylinder with a finite cylindrical humidity boundary. Radial
    # diffusion, finite attachment and isothermal heat have an exact time integral.
    air = dict(ice_saturation_pressure_Pa=100., vapor_gas_constant_J_kg_K=460.,
        vapor_molecular_speed_m_s=500., latent_sublimation_J_kg=2.8e6,
        diffusivity_m2_s=2e-5, conductivity_W_m_K=.02)
    T, si, alpha = 240., .1, .2
    R0, R1, Rout = 5e-6, 10e-6, 100e-6
    rhoeq = 100/(460*T)
    D = air['diffusivity_m2_s']
    X0 = D/125
    F = lambda r: .5*r*r*math.log(Rout/r)+.25*r*r
    rows = []
    for heat in (False, True):
        L = air['latent_sublimation_J_kg']
        chi = (L/(460*T)-1)/T*L*rhoeq*D/air['conductivity_W_m_K'] if heat else 0.
        exact = facet_growth.ICE_DENSITY/(rhoeq*D*si)*((1+chi)*(F(R1)-F(R0))+X0/alpha*(R1-R0))
        for dx in (1.25e-6, .625e-6, .3125e-6, .15625e-6):
            core = math.ceil(20e-6/dx)
            g = GrowingCrystal(dx, R0, dx, core, 1, extent=20., adaptive=False)
            g.rf = grid(core, Rout/dx)
            g.zf = np.array([0., 1.])
            g.far = Rout/dx
            g.volume = math.pi*np.diff(g.rf**2)[:, None]
            g.fraction = ((g.rf[:-1]+g.rf[1:])[:, None]/2*dx < R0).astype(float)
            g.solid = g.fraction > 0
            g.reflect_top = True
            g.initial_mass = g.mass()
            residual = 0.
            while g.state()['a_m'] < R1-dx/10:
                out = g.advance(air, T, si, [Attachment('constant', alpha)]*2,
                    1e10, kinks=False, heating=heat, tolerance=1e-11)
                residual = max(residual, out['mass_residual'])
            rows.append(dict(dx_m=dx, heating=heat, exact_time_s=exact,
                calculated_time_s=g.time, relative_time_error=g.time/exact-1,
                final=g.state(), max_water_residual=residual))
    return dict(air=air, temperature_K=T, ice_supersaturation=si, alpha=alpha,
        initial_radius_m=R0, final_radius_m=R1, outer_radius_m=Rout, cases=rows)


def capillary_controls(library, wall_limit):
    T = flakes.T0-15
    air = library.evaluate(.01, T, 100000., 2e-6, 5e-6, 10.)
    rows = []
    for extent in (40., 80.):
        for dx in (.25e-6, .125e-6, .0625e-6):
            for delta in (0., .3e-9, 1e-9):
                g = GrowingCrystal(dx, 2e-6, 5e-6, math.ceil(25e-6/dx), math.ceil(10e-6/dx), extent=extent)
                def forcing(_):
                    return dict(air=air, T=T, si=.01,
                        laws=[Attachment('constant', 1.), Attachment('nucleation', 1., .03)],
                        heating=False, capillary_length=delta)
                row = dict(dx_m=dx, capillary_length_m=delta, extent=extent,
                    **evolve(g, forcing, 6., 1., wall_limit))
                rows.append(row)
                print(json.dumps(dict(section='capillary', dx_m=dx, delta_m=delta,
                    extent=extent, complete=row['complete'], error=row['error'], final=row['final'])), flush=True)
    return dict(source='Libbrecht 2013 section 5.1 and Figure 5', air=air,
        temperature_K=T, pressure_Pa=100000., ice_supersaturation=.01,
        cases=rows, limits=[
            'Source-inspired qualitative control, not an exact reproduction: finite-volume/event solver and air convention differ.',
            'No latent heating, as in the source example; finite outer boundary is independently varied.',
            'Outer-terrace curvature proxy suppresses growth only; it does not implement sublimation or a general curvature law.',
            'Grid spacing resolves the approximate 0.15 micrometre kinetic length at the two finest levels.'])


def warm_restarts(library, wall_limit):
    tc = -5.3
    T = flakes.T0+tc
    initial = lab.observed(tc, 3.)
    a, c = initial['equatorial_dimension_m']/2, initial['polar_dimension_m']/2
    air = library.evaluate(.1, T, 101000., a, c, 10.)
    si = flakes.liquid_saturation_pressure(T)/air['ice_saturation_pressure_Pa']-1
    flow_air = replay.air_state(tc)
    nu = flow_air['viscosity_Pa_s']/flow_air['density_kg_m3']
    rows = []
    for dx in (1e-6, .5e-6):
        for label, sd in (('local', 0.), ('terraces', 14e-9)):
            g = GrowingCrystal(dx, a, c, math.ceil(150e-6/dx), math.ceil(1000e-6/dx), extent=80.)
            seed = g.measured_seed(a, c, initial['mass_kg'], 2e-6)
            def forcing(seconds):
                age = 3.+seconds/60
                fv = replay.jw_ventilation(lab.fall_speed(tc, age)*lab.observed(tc, age)['equatorial_dimension_m']/nu,
                    .63, replay.branch_habit(tc, age))
                return dict(air=air, T=T, si=si, laws=libbrecht_laws(tc),
                    surface_diffusion_length=sd, kinks=label == 'terraces',
                    capillary_length=1e-9, ventilation=fv)
            row = dict(dx_m=dx, label=label, surface_diffusion_length_m=sd,
                capillary_length_m=1e-9, seed=seed, **evolve(g, forcing, 300., 60., wall_limit))
            observed = lab.observed(tc, 3.+g.time/60)
            row['observed_fit_at_actual_final_age'] = observed
            row['final_to_observed_ratios_a_c_mass'] = [
                2*row['final']['a_m']/observed['equatorial_dimension_m'],
                2*row['final']['c_m']/observed['polar_dimension_m'],
                row['final']['mass_kg']/observed['mass_kg']]
            rows.append(row)
            print(json.dumps(dict(section='warm', dx_m=dx, label=label,
                complete=row['complete'], error=row['error'], final=row['final'],
                ratios=row['final_to_observed_ratios_a_c_mass'])), flush=True)
    return dict(temperature_K=T, pressure_Pa=101000., ice_supersaturation=si,
        air=air, ventilation_air=flow_air, observed_initial=initial, cases=rows,
        limits=[
            'The observed initial mass and dimensions constrain an unmeasured hollow profile; the two-micrometre rim is supplied.',
            'Fit ages 3 to 8 minutes; ventilation uses prescribed observed size and fall speed, not predicted drag.',
            '0.5 and 1 micrometre grids do not resolve the roughly 0.15 micrometre kinetic/tip length; these are coarse diagnostics.',
            'The 2015 terrace rule uses n_SD=x_s/molecular_step pixels; its physical collection distance changes with grid spacing.',
            'Local case disables rough-kink promotion; terrace case includes both rough kinks and terrace-distance enhancement.',
            'No adjoining-facet surface transfer, stochastic step sources, substrate, three-dimensional branching, riming or aggregation.',
            'Whole-fit dimensions and mass have no supplied joint uncertainty. No trajectory is adopted as an Earth or lunar growth law.'])


def run(section, wall_limit):
    library = facet_attachment.AttachmentLibrary(facet_growth.BUILD)
    result = dict(schema='terluna.atmosphere.spatial-crystal-assessment/1',
        evidence='numerical moving-interface controls and conditional coarse Earth restart; no morphology validation',
        numpy_version=np.__version__, scipy_version=scipy.__version__,
        per_case_wall_limit_s=wall_limit, inputs=[file_id(p) for p in [
            Path(__file__).resolve(), HERE/'spatial_crystal_growth.py', HERE/'spatial_diffusion.py',
            HERE/'spatial_facet_growth.py', HERE/'facet_attachment.py', HERE/'facet_growth.py',
            HERE/'takahashi_growth.py', HERE/'crystal_transport_replay.py',
            HERE/'spatial_growth_sources.json', ROOT/'shared/constants.json']])
    if section in ('all', 'controls'):
        result['cylinder_controls'] = cylinder_controls()
        result['capillary_controls'] = capillary_controls(library, wall_limit)
    if section in ('all', 'earth'):
        result['warm_restarts'] = warm_restarts(library, wall_limit)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--section', choices=['all', 'controls', 'earth'], default='all')
    parser.add_argument('--wall-limit', type=float, default=1200.)
    parser.add_argument('--output', type=Path, default=RESULT)
    args = parser.parse_args()
    if not math.isfinite(args.wall_limit) or args.wall_limit <= 0:
        parser.error('Positive finite wall limit required')
    result = run(args.section, args.wall_limit)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
