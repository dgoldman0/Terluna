"""How much warming inert fluorinated trace gases could add to the Moon's clear-sky greenhouse.

    python -m atmosphere.radiative_convective.trace_gases     # writes results/trace_gases.json

For each gas the calculation adds a well-mixed absorber that is grey across the
gas's strongest infrared band, at a range of band optical depths, to the Earth
control and to the lunar columns, and records the clear-sky drop in outgoing
longwave radiation (the instantaneous radiative forcing, no stratospheric
adjustment) line by line. Two numbers come from the literature: the band
limits and the IPCC AR6 radiative efficiency (Forster et al. 2021, Table
7.SM.7; W m^-2 per ppb, all-sky, stratosphere-adjusted, for Earth). The AR6
value fixes the band's mean cross-section through the Earth control's thin-limit
forcing, which then converts optical depth to mixing ratio on the Moon.

Why a band model: the Moon holds about six times more air above each square
metre per unit of surface pressure, so a part per billion is a six times thicker
absorber and the band centres saturate at a few ppb; the forcing then stops
growing in proportion to the amount. Real bands have wings and hot bands that
keep adding weakly past this point, so the grey-band curve is a lower bound for
large amounts, and the conversion inherits the difference between AR6's all-sky
adjusted value and this clear-sky instantaneous one (stated as a factor range).

Measured bands replace the grey band where cross-sections are at hand: SF6
and NF3 from the PNNL infrared database (Sharpe et al. 2004; 278, 298 and 323
K, 560 or 600 to 3000 cm^-1 on this grid) and CF4's strong band from NCAR
(Massie et al. 1991; 203-293 K, 1255-1290 cm^-1), as copied into AER's
cross-section database for LBLRTM (inputs.json). Each layer takes the
cross-sections interpolated to its temperature, held at the nearest measured
set outside their range, and the gas is added well mixed at 0.01-1000 ppb.
The Moon's curves need no conversion from optical depth, and their wings,
hot bands and weaker bands keep adding past the grey band's saturation.
All-sky estimates scale them by the factor that brings the Earth control's
thin-limit clear-sky forcing onto AR6's all-sky adjusted value. That factor
absorbs clouds, stratospheric adjustment and the overlapping gases these
columns lack (N2O, CH4 and O3). CF4's band lies under N2O and CH4, so its
factor is the smallest, and its lunar values hold only where those gases
are absent.
"""
from __future__ import annotations
import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys
import numpy as np

from atmosphere.radiative_convective import climate as cl, thermodynamics as th, optics, longwave as lw, fetch_inputs

HERE = Path(__file__).resolve().parent
# Strongest band (cm^-1) and AR6 radiative efficiency (W m^-2 ppb^-1). Lifetimes (years) from AR6 Table 7.SM.7.
GASES = {
    'SF6': dict(band=(925.0, 955.0), re_ar6=0.567, lifetime_yr=3200),
    'CF4': dict(band=(1270.0, 1290.0), re_ar6=0.099, lifetime_yr=50000),
    'NF3': dict(band=(890.0, 920.0), re_ar6=0.204, lifetime_yr=569),
}
TAUS = (1e-3, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0)
CLEAR_TO_ALLSKY = (0.6, 0.9)     # plausible all-sky adjusted / clear-sky instantaneous ratio for window absorbers
# Measured cross-section sets (HITRAN format, from AER's LBLRTM database; see inputs.json).
MEASURED = {'SF6': ('aer_xs_SF6AT1', 'aer_xs_SF6AT2', 'aer_xs_SF6AT3'),
            'NF3': ('aer_xs_NF3AT1', 'aer_xs_NF3AT2', 'aer_xs_NF3AT3'),
            'CF4': tuple(f'aer_xs_F14AT{i}' for i in range(1, 7))}
PPB = (0.01, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0)


def base_optics(scenario, ts, cfg=optics.LONGWAVE):
    col = scenario.column(ts, scenario.longwave_levels)
    nu, tau = optics.optical_depth(col, cfg)
    return col, nu, tau


def olr(col, nu, tau):
    up, _ = lw.fluxes(nu, tau, col['t_k'], col['t_k'][0])
    return float(lw.integrate(nu, up)[-1])


def forcing_curve(col, nu, tau, band):
    """OLR drop for a well-mixed grey absorber of total optical depth t across `band` (per unit mass column)."""
    inside = (nu >= band[0]) & (nu <= band[1])
    share = col['layer_mass_kg_m2'] / col['layer_mass_kg_m2'].sum()
    ref = olr(col, nu, tau)
    out = []
    for t in TAUS:
        extra = np.zeros_like(tau)
        extra[:, inside] = t * share[:, None]
        out.append(ref - olr(col, nu, tau + extra))
    return ref, np.array(out)


def read_cross_sections(name):
    """One HITRAN-format set: temperature (K), wavenumbers (cm^-1) and cross-sections (cm^2 per molecule).
    Values follow the header in fixed ten-character fields, ten to a line; anything after the stated
    number of points (padding, a DOS end-of-file byte) is ignored."""
    lines = fetch_inputs.path(name).read_text().splitlines()
    head = lines[0].split()
    nu_min, nu_max, npts, t = float(head[1]), float(head[2]), int(head[3]), float(head[4])
    values = []
    for line in lines[1:]:
        values.extend(float(line[i:i + 10]) for i in range(0, len(line), 10) if line[i:i + 10].strip())
        if len(values) >= npts:
            break
    if len(values) < npts:
        raise ValueError(f'{name}: {len(values)} values for {npts} points')
    return t, np.linspace(nu_min, nu_max, npts), np.maximum(np.array(values[:npts]), 0.0)


def measured_tau_per_ppb(col, nu, gas):
    """Optical depth (layers x wavenumbers) of 1 ppb of `gas`, well mixed in dry air, with each layer's
    cross-sections interpolated linearly in temperature between the measured sets."""
    sets = sorted((read_cross_sections(n) for n in MEASURED[gas]), key=lambda s: s[0])
    temps = np.array([s[0] for s in sets])
    on_grid = [np.interp(nu, s[1], s[2], left=0.0, right=0.0) for s in sets]
    dry = col['layer_column_cm2'] * (1 - col['layer_x_h2o'])
    tau = np.empty((dry.size, nu.size))
    for k, t in enumerate(col['layer_t_k']):
        t = min(max(float(t), temps[0]), temps[-1])
        j = min(int(np.searchsorted(temps, t)), temps.size - 1)
        i = max(j - 1, 0)
        w = 0.0 if i == j else (t - temps[i]) / (temps[j] - temps[i])
        tau[k] = ((1 - w) * on_grid[i] + w * on_grid[j]) * dry[k] * 1e-9
    return tau


def measured_forcing(col, nu, tau, gas):
    """Clear-sky OLR drop (W/m^2) at each of PPB for the measured bands of `gas`."""
    per_ppb = measured_tau_per_ppb(col, nu, gas)
    ref = olr(col, nu, tau)
    return ref, np.array([ref - olr(col, nu, tau + x * per_ppb) for x in PPB])


def _ppb_at(curve, target):
    """Mixing ratio (ppb) at which a measured forcing curve reaches `target` W/m^2, or None beyond 1000 ppb."""
    if curve[-1] < target:
        return None
    return float(np.exp(np.interp(target, curve, np.log(PPB))))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--ts', type=float, default=288.0)
    parser.add_argument('--out', type=Path, default=HERE / 'results' / 'trace_gases.json')
    parser.add_argument('--processes', type=int, default=0)
    args = parser.parse_args(argv)
    cfg = replace(optics.LONGWAVE, processes=args.processes)
    scenarios = dict(earth_1atm=cl.Scenario('earth_1.0atm', th.EARTH, 101325.0, 400.0),
                     moon_1atm=cl.Scenario('moon_1.0atm', th.MOON, 101325.0, 400.0),
                     moon_12atm=cl.Scenario('moon_1.2atm', th.MOON, 121590.0, 400.0))
    columns = {name: base_optics(sc, args.ts, cfg) for name, sc in scenarios.items()}
    result = dict(schema='terluna.atmosphere.trace-gas-forcing/1', surface_temperature_k=args.ts,
                  band_optical_depths=list(TAUS), evidence=__doc__.split('\n\n')[1].replace('\n', ' '),
                  gases={})
    for gas, info in GASES.items():
        entry = dict(band_cm1=info['band'], re_ar6_w_m2_ppb=info['re_ar6'], lifetime_yr=info['lifetime_yr'],
                     columns={})
        curves = {}
        for name, (col, nu, tau) in columns.items():
            ref, curve = forcing_curve(col, nu, tau, info['band'])
            molecules = float(np.sum(col['layer_column_cm2'] * (1 - col['layer_x_h2o'])))
            curves[name] = (ref, curve, molecules)
        # Thin-limit clear-sky forcing per unit band optical depth on the Earth control.
        e_ref, e_curve, e_molecules = curves['earth_1atm']
        slope_earth = e_curve[0] / TAUS[0]
        for name, (ref, curve, molecules) in curves.items():
            # AR6 all-sky adjusted RE = ratio x clear-sky RE, so the band optical depth per ppb on Earth is
            # tau_ppb_earth = RE_AR6 / (ratio x slope_earth); it scales with the air column.
            per_ppb = {f'{r:g}': info['re_ar6'] / (r * slope_earth) * molecules / e_molecules for r in CLEAR_TO_ALLSKY}
            # The same ratio converts this column's clear-sky forcing to an all-sky adjusted estimate.
            entry['columns'][name] = dict(
                clear_sky_olr=ref, clear_sky_forcing_w_m2=curve.round(4).tolist(),
                thin_limit_w_m2_per_unit_tau=float(curve[0] / TAUS[0]),
                band_tau_per_ppb=per_ppb,
                thin_limit_allsky_w_m2_per_ppb={f'{r:g}': float(r * curve[0] / TAUS[0] * per_ppb[f'{r:g}'])
                                                for r in CLEAR_TO_ALLSKY},
                ppb_for_allsky_forcing={f'{r:g}': {f'{f:g}': _ppb_for(r * curve, per_ppb[f'{r:g}'], f)
                                                   for f in (1, 2, 5, 10, 20)} for r in CLEAR_TO_ALLSKY},
                saturated_band_clear_sky_w_m2=float(curve[-1]))
        result['gases'][gas] = entry
        print(gas, {n: (round(v['thin_limit_allsky_w_m2_per_ppb']['0.9'], 3), round(v['saturated_band_clear_sky_w_m2'], 2),
                        v['ppb_for_allsky_forcing']['0.9']) for n, v in entry['columns'].items()}, flush=True)
    result['measured'] = dict(mixing_ratios_ppb=list(PPB), gases={})
    for gas in MEASURED:
        curves = {name: measured_forcing(col, nu, tau, gas)[1] for name, (col, nu, tau) in columns.items()}
        # AR6's all-sky adjusted efficiency over this model's clear-sky Earth value converts every curve.
        ratio = GASES[gas]['re_ar6'] / float(curves['earth_1atm'][0] / PPB[0])
        entry = dict(sets_k=sorted(read_cross_sections(n)[0] for n in MEASURED[gas]), ar6_over_clear_sky_earth=ratio,
                     columns={})
        for name, curve in curves.items():
            entry['columns'][name] = dict(
                clear_sky_forcing_w_m2=curve.round(4).tolist(),
                thin_limit_clear_sky_w_m2_per_ppb=float(curve[0] / PPB[0]),
                allsky_forcing_w_m2=(ratio * curve).round(4).tolist(),
                ppb_for_allsky_forcing={f'{f:g}': _ppb_at(ratio * curve, f) for f in (1, 2, 5, 10, 20)})
        result['measured']['gases'][gas] = entry
        print(gas, 'measured', {n: (round(v['thin_limit_clear_sky_w_m2_per_ppb'], 3), v['clear_sky_forcing_w_m2'][-1])
                                for n, v in entry['columns'].items()},
              'AR6/clear-sky on Earth', round(ratio, 2), flush=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=1) + '\n')
    return 0


def _ppb_for(curve, tau_per_ppb, target):
    """Mixing ratio (ppb) whose grey band gives `target` W/m^2, or None beyond the saturated band."""
    if curve[-1] < target:
        return None
    tau = float(np.exp(np.interp(target, curve, np.log(TAUS))))
    return tau / tau_per_ppb


if __name__ == '__main__':
    sys.exit(main())
