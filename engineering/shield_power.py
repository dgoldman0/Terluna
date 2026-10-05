"""Explicit optical-to-electric and material sensitivities for shield accounts."""


def routed_power(redirected_W, capture, conversion, delivery, propulsion_W=0., other_W=0.):
    """Operations draw from the collector bus before the final delivery link.

    Capture excludes missed beams. Unconverted captured optical power is an
    ideal collector heat account; optical rejection or thermal routing needs
    device-specific treatment. No hardware feasibility is implied.
    """
    if any(not 0 <= e <= 1 for e in (capture, conversion, delivery)):
        raise ValueError('Efficiencies must lie in [0, 1]')
    if min(redirected_W, propulsion_W, other_W) < 0:
        raise ValueError('Powers must be nonnegative')
    captured = redirected_W*capture
    gross = captured*conversion
    bus = gross-propulsion_W-other_W
    return dict(redirected_optical_W=redirected_W, captured_optical_W=captured,
        uncaptured_optical_W=redirected_W-captured, gross_electric_W=gross,
        ideal_collector_waste_heat_W=captured-gross, propulsion_W=propulsion_W,
        other_operations_W=other_W, net_bus_W=bus,
        delivered_W=max(bus,0.)*delivery,
        final_link_loss_W=max(bus,0.)*(1-delivery), deficit_W=max(-bus,0.))


def optical_requirement(delivered_W, capture, conversion, delivery, propulsion_W=0., other_W=0.):
    if min(capture,conversion,delivery) <= 0 or max(capture,conversion,delivery)>1:
        raise ValueError('Positive efficiencies at most one required')
    if min(delivered_W,propulsion_W,other_W)<0:
        raise ValueError('Powers must be nonnegative')
    return (propulsion_W+other_W+delivered_W/delivery)/(capture*conversion)


def exhaust_rate(propulsion_W, exhaust_m_s, efficiency):
    """Consumed exhaust from bus-to-jet power; cant already enters the load."""
    if propulsion_W<0 or exhaust_m_s<=0 or not 0<efficiency<=1:
        raise ValueError('Invalid propulsion account')
    return 2*efficiency*propulsion_W/exhaust_m_s**2
