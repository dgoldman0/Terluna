"""Explicit band energy accounting; refuse unprovided spectral information.

The existing optical archive leaves a ~24.8--120.18 nm material-data gap, and
its solar table starts at 202 nm. This interface never interpolates across gaps.
Heating, chemistry and deposited radius must be supplied with provenance.
"""
from __future__ import annotations
from dataclasses import dataclass
import math

@dataclass(frozen=True)
class Band:
    lower_nm: float
    upper_nm: float
    irradiance_w_m2: float | None
    transmission: float | None
    atmospheric_absorptance: float | None
    heat_fraction: float | None
    absorption_radius_R: float | None
    provenance: str


def deposited_heat(bands,required_interval_nm):
    if not bands: raise ValueError('No spectral data: atmospheric result BLOCKED')
    lo,hi=required_interval_nm
    ordered=sorted(bands,key=lambda x:x.lower_nm)
    cursor=lo;total=0.;ledger=[]
    for b in ordered:
        if abs(b.lower_nm-cursor)>1e-9 or b.upper_nm<=b.lower_nm:
            raise ValueError('Spectral gap or overlap: atmospheric result BLOCKED')
        fields=(b.irradiance_w_m2,b.transmission,b.atmospheric_absorptance,b.heat_fraction,b.absorption_radius_R)
        if any(v is None for v in fields) or not b.provenance:
            raise ValueError('Unspecified spectral/chemical input: atmospheric result BLOCKED')
        if not all(math.isfinite(v) for v in fields):raise ValueError('Nonfinite band data')
        f,tr,absorbed,eta,r=fields
        if f<0 or r<1 or any(not 0<=v<=1 for v in (tr,absorbed,eta)):
            raise ValueError('Unphysical band input')
        incident=f*r*r/4  # cross-section / lunar surface area
        absorbed_power=incident*tr*absorbed
        heat=absorbed_power*eta
        ledger.append(dict(lower_nm=b.lower_nm,upper_nm=b.upper_nm,
            incoming_per_surface_W_m2=incident,absorbed_W_m2=absorbed_power,
            heat_W_m2=heat,other_energy_W_m2=absorbed_power-heat,provenance=b.provenance))
        total+=heat;cursor=b.upper_nm
    if abs(cursor-hi)>1e-9:raise ValueError('Incomplete requested wavelength range')
    return total,ledger


def inventory_audit():
    return dict(status='BLOCKED_SPECTRUM_CHEMISTRY_CLOSURE',
        solar_irradiance_nm=[202.0,2729.9],
        titania_index_nm=[120.181141,125122.7623],
        xray_model_min_energy_ev=50.0,
        xray_longest_evaluated_nm=1239.841984/50,
        uncovered_material_interval_nm=[1239.841984/50,120.181141],
        missing=['Resolved solar EUV/FUV irradiance below 202 nm; historical integrated bands available separately',
                 'Optical response in the gap, with relevant material state',
                 'Species absorption/photoelectron heating/photochemical rates',
                 'Atomic transport and infrared cooling',
                 'Kinetic/momentum closure when molecular Jeans approximation fails'],
        consequence='Visible transmission and Lyman-alpha opacity cannot determine total XUV heating or escape.')


def audit_optical_files(solar_text,titania_text):
    """Measure the supplied input ranges, retaining the inherited model's 50 eV X-ray limit."""
    import csv,io,hashlib
    solar=list(csv.reader(io.StringIO(solar_text)))
    wave=[float(r[0]) for r in solar[1:] if r]
    nk=[]
    for line in titania_text.splitlines():
        row=line.split()
        if len(row)==3:
            try:nk.append([float(v) for v in row])
            except ValueError:pass
    if not wave or not nk:raise ValueError('Unreadable spectral inputs')
    result=inventory_audit()
    result['solar_irradiance_nm']=[min(wave),max(wave)]
    result['titania_index_nm']=[min(r[0] for r in nk)*1000,max(r[0] for r in nk)*1000]
    result['uncovered_material_interval_nm'][1]=result['titania_index_nm'][0]
    expected={'solar':'0ae235c6912f26773c0e9b71ea3924aaf8856c74a511dc07e86b60b3c1b88172',
              'titania':'9b03947949e54cd98f58a989bd031a5fddb1a36a9d9f7fd64d5ff2bd68c2885c'}
    actual={'solar':hashlib.sha256(solar_text.encode()).hexdigest(),
            'titania':hashlib.sha256(titania_text.encode()).hexdigest()}
    if actual!=expected:raise ValueError('Historical optical input hash mismatch; audit BLOCKED')
    result['actual_inputs_inspected']=True
    result['input_sha256']={'solar':hashlib.sha256(solar_text.encode()).hexdigest(),
                           'titania':hashlib.sha256(titania_text.encode()).hexdigest()}
    return result


def historical_solar_bands():
    """Ribas et al. Table 4, Sun column; integrated energy, not resolved spectra.

    The solar reference uses 1993 mid-cycle data plus complementary/inferred
    bands (section 3.5). Angstrom -> nm /10; erg/s/cm2 -> W/m2 *0.001.
    Transmission, absorption, heating efficiency and radius remain unspecified.
    """
    values=[(.1,2.,.15),(2.,10.,.70),(10.,36.,2.05),(36.,92.,1.00),(92.,118.,.74)]
    return [Band(a,b,f*.001,None,None,None,None,
        'Ribas et al. arXiv:astro-ph/0412253v1, Table 4 Sun column; historical composite')
        for a,b,f in values]


def gap_energy_bounds(bands, gap):
    """Bounds without assuming where irradiance lies inside a reported band."""
    lower=upper=total=0.
    for band in bands:
        f=band.irradiance_w_m2
        if f is None or not math.isfinite(f) or f<0:raise ValueError('Valid band energy required')
        total+=f
        if band.lower_nm>=gap[0] and band.upper_nm<=gap[1]:lower+=f
        if band.upper_nm>gap[0] and band.lower_nm<gap[1]:upper+=f
    if total<=0:raise ValueError('Positive total energy required')
    return dict(total_irradiance_W_m2=total,gap_energy_lower_bound_W_m2=lower,
        gap_energy_upper_bound_W_m2=upper,lower_fraction=lower/total,upper_fraction=upper/total,
        within_band_distribution_assumed=False)
