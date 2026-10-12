"""Mechanical requirements for hypothetical living aerophytes.

These pure functions are screening models, not a solved balloon shape or a
material design. Stress allowances are assembled, wet, aged design values;
measured coupon breaking stress must not be substituted without qualification.
"""
from __future__ import annotations

import math
from shared.constants import MOON_GM, MOON_RADIUS


def _positive(**values):
    if any(not math.isfinite(v) or v <= 0 for v in values.values()):
        raise ValueError(f'Expected finite positive parameters: {values}')


def _nonnegative(**values):
    if any(not math.isfinite(v) or v < 0 for v in values.values()):
        raise ValueError(f'Expected finite nonnegative parameters: {values}')


def gravity_at_height(height_m):
    """Lunar GM/r²; heights refer to the same spherical sea-level datum as air."""
    _nonnegative(height_m=height_m)
    return MOON_GM / (MOON_RADIUS + height_m)**2


def pressure_loads(radius_m, air_density_kg_m3, gas_density_kg_m3, gravity_m_s2,
                   bottom_overpressure_pa=0., gust_m_s=0., gust_cp=1.):
    """Full-gas vertical head and additive local aerodynamic load sensitivity.

    d(p_inside-p_outside)/dz=(rho_air-rho_gas)*g. Bottom pressure is
    nonnegative. Crown pressure is used uniformly when sizing later screens;
    this overbounds local pressure, but does not solve wrinkling/load paths.
    A co-moving organism can still meet gust/shear across its finite body.
    """
    _positive(radius_m=radius_m, air_density=air_density_kg_m3,
              gas_density=gas_density_kg_m3, gravity=gravity_m_s2)
    _nonnegative(bottom_pressure=bottom_overpressure_pa, gust=gust_m_s, cp=gust_cp)
    if gas_density_kg_m3 >= air_density_kg_m3:
        raise ValueError('Gas must be lighter than surrounding air')
    contrast = air_density_kg_m3 - gas_density_kg_m3
    head = 2 * radius_m * contrast * gravity_m_s2
    q = .5 * air_density_kg_m3 * gust_m_s**2
    return dict(hydrostatic_head_pa=head, gust_dynamic_pressure_pa=q,
                gust_pressure_allowance_pa=gust_cp*q,
                bottom_overpressure_pa=bottom_overpressure_pa,
                crown_static_pressure_pa=bottom_overpressure_pa+head,
                design_differential_pressure_pa=bottom_overpressure_pa+head+gust_cp*q,
                pressure_gradient_pa_m=contrast*gravity_m_s2,
                lift_density_kg_m3=contrast)


def isothermal_hydrostatic_head(radius_m, centre_height_m, centre_pressure_pa,
                                 centre_air_density_kg_m3, centre_gas_density_kg_m3,
                                 bottom_overpressure_pa=0.):
    """Independent large-height check; both gases isothermal, g varies as GM/r².

    Gas density is its value at the same reference pressure as ambient air.
    This is an analytic consistency check, not a resolved Terluna temperature
    profile. The actual air table/thermal model determines applicability.
    """
    _positive(radius_m=radius_m, centre_pressure=centre_pressure_pa,
              air_density=centre_air_density_kg_m3, gas_density=centre_gas_density_kg_m3)
    _nonnegative(centre_height=centre_height_m, bottom_pressure=bottom_overpressure_pa)
    if radius_m > centre_height_m:
        raise ValueError('Envelope would intersect the reference surface')
    if centre_gas_density_kg_m3 >= centre_air_density_kg_m3:
        raise ValueError('Gas must be lighter than surrounding air')
    rc = MOON_RADIUS + centre_height_m
    rb, rt = rc-radius_m, rc+radius_m
    phi_bottom = MOON_GM*(1/rc-1/rb)
    phi_top = MOON_GM*(1/rc-1/rt)
    ra = centre_air_density_kg_m3 / centre_pressure_pa
    rg = centre_gas_density_kg_m3 / centre_pressure_pa
    pa_bottom = centre_pressure_pa*math.exp(-ra*phi_bottom)
    pa_top = centre_pressure_pa*math.exp(-ra*phi_top)
    pg_top = (pa_bottom+bottom_overpressure_pa)*math.exp(-rg*(phi_top-phi_bottom))
    return dict(ambient_bottom_pa=pa_bottom, ambient_top_pa=pa_top,
                gas_top_pa=pg_top, crown_pressure_pa=pg_top-pa_top,
                local_density_scale_height_m=centre_pressure_pa/(centre_air_density_kg_m3*gravity_at_height(centre_height_m)),
                body_height_over_density_scale_height=2*radius_m*centre_air_density_kg_m3*gravity_at_height(centre_height_m)/centre_pressure_pa)


def spherical_skin(radius_m, differential_pressure_pa, allowable_stress_pa,
                   material_density_kg_m3, minimum_thickness_m=0.,
                   barrier_areal_mass_kg_m2=0.):
    """Uniform sphere under the supplied bounding pressure: sigma=pR/(2t).

    Structural thickness and an additional barrier are separate. If the same
    layer serves both roles, use its thickness floor and zero additional barrier.
    Density must include the structural layer's hydration in the stated thickness.
    """
    _positive(radius_m=radius_m, allowable=allowable_stress_pa, density=material_density_kg_m3)
    _nonnegative(pressure=differential_pressure_pa, thickness=minimum_thickness_m,
                 barrier=barrier_areal_mass_kg_m2)
    stress_thickness = differential_pressure_pa*radius_m/(2*allowable_stress_pa)
    thickness = max(minimum_thickness_m, stress_thickness)
    area = 4*math.pi*radius_m**2
    film_mass = area*thickness*material_density_kg_m3
    barrier_mass = area*barrier_areal_mass_kg_m2
    return dict(area_m2=area, structural_thickness_m=thickness,
                stress_required_thickness_m=stress_thickness,
                structural_areal_mass_kg_m2=thickness*material_density_kg_m3,
                structural_mass_kg=film_mass, additional_barrier_mass_kg=barrier_mass,
                total_structure_mass_kg=film_mass+barrier_mass,
                membrane_resultant_n_m=differential_pressure_pa*radius_m/2,
                thickness_over_radius=thickness/radius_m,
                thin_shell_screen_applicable=thickness/radius_m <= .01)


def lobed_skin_tendons(radius_m, differential_pressure_pa, lobe_radius_m,
                       film_allowable_pa, film_density_kg_m3,
                       tendon_allowable_pa, tendon_density_kg_m3,
                       minimum_film_thickness_m=0., barrier_areal_mass_kg_m2=0.,
                       area_multiplier=1.1, additional_supported_force_n=0.):
    """Ideal lobed film plus meridional tendons; equivalent spherical volume.

    Film hoop resultant p*r_lobe (cylindrical panel approximation). N equal
    meridional tendons, each length pi*R, share p*pi*R²+extra_force at the equator.
    No anchorage, junction, tendon sleeve, internal partition, loss-of-tendon,
    creep incompatibility or shape-stability mass is hidden in these estimates.
    The additional force is an optional independently specified load-path allowance,
    not an instruction to double-count payload already represented by pressure.
    """
    _positive(radius_m=radius_m, lobe_radius=lobe_radius_m, film_allowable=film_allowable_pa,
              film_density=film_density_kg_m3, tendon_allowable=tendon_allowable_pa,
              tendon_density=tendon_density_kg_m3, area_multiplier=area_multiplier)
    _nonnegative(pressure=differential_pressure_pa, thickness=minimum_film_thickness_m,
                 barrier=barrier_areal_mass_kg_m2, extra_force=additional_supported_force_n)
    if lobe_radius_m > radius_m or area_multiplier < 1:
        raise ValueError('Lobe radius must not exceed body radius; area multiplier must be >=1')
    count = max(4, math.ceil(math.pi*radius_m/lobe_radius_m))
    film_thickness = max(minimum_film_thickness_m,
                         differential_pressure_pa*lobe_radius_m/film_allowable_pa)
    area = area_multiplier*4*math.pi*radius_m**2
    force = differential_pressure_pa*math.pi*radius_m**2+additional_supported_force_n
    tendon_force = force/count
    tendon_cross_section = tendon_force/tendon_allowable_pa
    length = math.pi*radius_m
    tendon_mass = count*length*tendon_cross_section*tendon_density_kg_m3
    film_mass = area*film_thickness*film_density_kg_m3
    barrier_mass = area*barrier_areal_mass_kg_m2
    return dict(film_area_m2=area, structural_thickness_m=film_thickness,
                structural_areal_mass_kg_m2=film_thickness*film_density_kg_m3,
                film_mass_kg=film_mass, additional_barrier_mass_kg=barrier_mass,
                tendon_count=count, tendon_length_m=length,
                tendon_force_n=tendon_force, total_equatorial_force_n=force,
                tendon_cross_section_m2=tendon_cross_section,
                tendon_mass_kg=tendon_mass,
                total_structure_mass_kg=film_mass+barrier_mass+tendon_mass,
                film_thickness_over_local_radius=film_thickness/lobe_radius_m,
                thin_panel_screen_applicable=film_thickness/lobe_radius_m <= .01,
                load_path_and_shape_validated=False)


def supported_payload(radius_m, lift_density_kg_m3, structure_mass_kg,
                      loaded_lift_fraction=1., additional_mass_kg=0.):
    """Remaining wet biomass allowance. Negative values reject this mass screen."""
    _positive(radius_m=radius_m, lift_density=lift_density_kg_m3)
    _nonnegative(structure=structure_mass_kg, additional_mass=additional_mass_kg)
    if not 0 < loaded_lift_fraction <= 1:
        raise ValueError('Loaded lift fraction must be in (0,1]')
    volume = 4*math.pi*radius_m**3/3
    gross = lift_density_kg_m3*volume
    remaining = loaded_lift_fraction*gross-structure_mass_kg-additional_mass_kg
    return dict(full_gross_supported_mass_kg=gross, loaded_lift_fraction=loaded_lift_fraction,
                unused_lift_capacity_kg=(1-loaded_lift_fraction)*gross,
                remaining_wet_payload_kg=remaining,
                remaining_wet_payload_kg_per_projected_m2=remaining/(math.pi*radius_m**2),
                mass_screen_pass=remaining >= 0)


def spherical_payload_window(lift_density_kg_m3, gravity_m_s2, allowable_stress_pa,
                             material_density_kg_m3, target_wet_payload_kg_m2=0.,
                             loaded_lift_fraction=1., extra_pressure_pa=0.,
                             barrier_areal_mass_kg_m2=0.):
    """Analytic size window for crown-sized, stress-limited homogeneous spheres.

    B(R)=aR-bR²-c, per projected area. Excludes a minimum-thickness floor,
    non-spherical shapes and nonstructural masses; it is a screen-specific bound,
    never a universal maximum size. Required film thickness must later be checked.
    """
    _positive(lift=lift_density_kg_m3, gravity=gravity_m_s2,
              allowable=allowable_stress_pa, density=material_density_kg_m3)
    _nonnegative(target=target_wet_payload_kg_m2, extra_pressure=extra_pressure_pa,
                 barrier=barrier_areal_mass_kg_m2)
    if not 0 < loaded_lift_fraction <= 1:
        raise ValueError('Loaded lift fraction must be in (0,1]')
    a = 4*loaded_lift_fraction*lift_density_kg_m3/3 - 2*material_density_kg_m3*extra_pressure_pa/allowable_stress_pa
    b = 4*material_density_kg_m3*lift_density_kg_m3*gravity_m_s2/allowable_stress_pa
    c = 4*barrier_areal_mass_kg_m2
    optimum = max(0.,a/(2*b))
    peak = a*optimum-b*optimum**2-c
    disc = a*a-4*b*(c+target_wet_payload_kg_m2)
    feasible = a > 0 and disc >= 0
    window = [(a-math.sqrt(disc))/(2*b),(a+math.sqrt(disc))/(2*b)] if feasible else None
    return dict(optimum_radius_m=optimum, maximum_wet_payload_kg_per_projected_m2=peak,
                target_wet_payload_kg_per_projected_m2=target_wet_payload_kg_m2,
                feasible_radius_window_m=window, homogeneous_sphere_screen_feasible=feasible,
                polynomial_coefficients=dict(a=a,b=b,c=c), universal_size_bound=False)


def neutral_trim(radius_m, air_density_kg_m3, lifting_mix_density_kg_m3,
                 supported_mass_kg, hydrogen_mole_fraction=.98):
    """Exact constant-density trim via a mixed-gas bladder and ambient-air ballonet.

    supported_mass includes every actual solid/wet mass and emergency inventory,
    including partitions. Unused lift capacity is excluded from actual mass.
    Hydrogen mole/volume fraction is within the lifting mixture; the other
    fraction is ambient air, consistently with the sky_ships98%-purity product.
    Leakage barrier area bounds are geometry only, not gas-transfer predictions.
    """
    _positive(radius_m=radius_m, air_density=air_density_kg_m3, gas_density=lifting_mix_density_kg_m3)
    _nonnegative(supported_mass=supported_mass_kg)
    if not 0 < hydrogen_mole_fraction <= 1 or lifting_mix_density_kg_m3 >= air_density_kg_m3:
        raise ValueError('Invalid lifting mixture')
    species_density = lifting_mix_density_kg_m3-(1-hydrogen_mole_fraction)*air_density_kg_m3
    if species_density < 0:
        raise ValueError('Mixture density incompatible with stated air fraction')
    volume = 4*math.pi*radius_m**3/3
    area = 4*math.pi*radius_m**2
    gross = (air_density_kg_m3-lifting_mix_density_kg_m3)*volume
    fraction = supported_mass_kg/gross
    out = dict(full_gross_supported_mass_kg=gross,
               required_lifting_mix_volume_fraction=fraction,
               required_hydrogen_volume_fraction=fraction*hydrogen_mole_fraction,
               neutral_trim_possible=fraction <= 1,
               full_mixture_mass_kg=lifting_mix_density_kg_m3*volume,
               full_hydrogen_mass_kg=species_density*volume,
               full_displaced_air_mass_kg=air_density_kg_m3*volume,
               unused_lift_capacity_kg=gross-supported_mass_kg)
    if fraction > 1:
        return dict(**out, actual_mixture_mass_kg=None, actual_hydrogen_mass_kg=None,
                    ballonet_air_mass_kg=None, total_mass_closure_kg=None,
                    gas_barrier_area_spherical_minimum_m2=None, gas_barrier_area_full_exterior_m2=area)
    gas_mass = lifting_mix_density_kg_m3*volume*fraction
    ballast = air_density_kg_m3*volume*(1-fraction)
    return dict(**out,actual_mixture_mass_kg=gas_mass,
                actual_hydrogen_mass_kg=species_density*volume*fraction,
                ballonet_air_mass_kg=ballast,
                total_mass_closure_kg=supported_mass_kg+gas_mass+ballast-air_density_kg_m3*volume,
                gas_barrier_area_spherical_minimum_m2=area*fraction**(2/3),
                gas_barrier_area_full_exterior_m2=area)


def partition_mass(radius_m, partition_areal_mass_kg_m2, projected_area_multiples=1.):
    """Explicit planar partition sensitivity; one great-circle sheet=piR².

    Shape/attachment/stress and whether this also supplies a gas barrier are
    unsolved; users must not add an identical layer again elsewhere.
    """
    _positive(radius_m=radius_m)
    _nonnegative(areal_mass=partition_areal_mass_kg_m2, multiples=projected_area_multiples)
    return math.pi*radius_m**2*projected_area_multiples*partition_areal_mass_kg_m2


def initial_hole_leak(hole_area_m2, differential_pressure_pa, gas_density_kg_m3,
                      ambient_pressure_pa, discharge_coefficient=.6):
    """Initial continuum subsonic orifice outflow; excludes tear growth/diffusion.

    Q=Cd*A*sqrt(2*dp/rho). Valid only for dp/pambient<=0.1; outside this
    bounded use a compressible flow model. A zero pressure result means no
    pressure-driven flow, not zero diffusive leakage through an open hole.
    """
    _positive(gas_density=gas_density_kg_m3, ambient_pressure=ambient_pressure_pa)
    _nonnegative(area=hole_area_m2, differential_pressure=differential_pressure_pa)
    if not 0 < discharge_coefficient <= 1:
        raise ValueError('Invalid discharge coefficient')
    ratio = differential_pressure_pa/ambient_pressure_pa
    if ratio > .1:
        raise ValueError('Pressure ratio outside low-compressibility screen')
    volume_s = discharge_coefficient*hole_area_m2*math.sqrt(2*differential_pressure_pa/gas_density_kg_m3)
    return dict(initial_volume_m3_s=volume_s,
                initial_gas_mass_kg_s=volume_s*gas_density_kg_m3,
                overpressure_fraction_of_ambient=ratio,
                tear_propagation_modelled=False, diffusion_modelled=False)


def residual_descent_speed(excess_mass_kg, gravity_m_s2, air_density_kg_m3,
                           projected_area_m2, drag_coefficient=.47):
    """Drag-balanced descent from residual negative buoyancy after any active trim.

    Still-air, steady-speed diagnostic only: excludes added mass, acceleration,
    deformation, turbulent up/down drafts and impacts. Rain that is drained or
    compensated by gas/air trim is not counted again as residual excess mass.
    """
    _nonnegative(excess_mass=excess_mass_kg)
    _positive(gravity=gravity_m_s2, air_density=air_density_kg_m3,
              area=projected_area_m2, cd=drag_coefficient)
    return math.sqrt(2*excess_mass_kg*gravity_m_s2 /
                     (air_density_kg_m3*drag_coefficient*projected_area_m2))
