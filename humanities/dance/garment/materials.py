"""Fabrics for the dance garments. Every number and its basis is recorded in sources.json beside this file.

Areal density and bending rigidity come in measured pairs where a measurement was found: silk chiffon and silk
georgette (Fuss 2024, cantilever bending to ASTM D1388 on one fabric each), and a plain-weave wool and a plain-weave
cotton (Sun et al. 2017, KES-F2; the solver's isotropic bending takes the mean of warp and weft). A silk satin has a
measured weight and drape coefficient (Lee et al. 2011), and its bending rigidity is the one that makes the solver's
drape test give that drape. Crepe de chine and
the heavier satin rest on stated brackets; their central values are modelling choices inside those brackets. Organza,
habotai, taffeta, velvet, tulle and viscose jersey have weights from supplier specifications and no bending
measurement found: their bending rigidity is a modelling choice, set from a chosen bending length on Earth (a crisp
sheer, a soft plain weave, a crisp taffeta, a heavy pile, a sized net, a limp knit) and bracketed in sources.json, so
any result that rests on their stiffness carries that choice. Tulle's open weave takes a share of a closed sheet's
normal pressure (normal_factor) from porous-plate drag. Shear stiffness is a bracket from the usual KES-F range for
light woven fabrics (G about 0.2-1 gf/cm/deg); tensile stiffness only has to make woven cloth effectively inextensible
at the loads met here (strains below 1e-3), which it does; the knit and the net are given softer values.
"""
from __future__ import annotations

import json
from pathlib import Path

from .cloth import Material

GF_CM_PER_DEG = 0.980665 * 180.0 / 3.141592653589793      # 1 gf/cm per degree of shear in N/m per radian

FABRICS = {
    # measured pair: 26.1 g/m2, G_P = 0.502 uN m (Fuss 2024, Table 1)
    'silk_chiffon': Material('silk_chiffon', areal_density=0.0261, bending_rigidity=0.502e-6,
                             shear_stiffness=0.3 * GF_CM_PER_DEG, tensile_stiffness=5000.0, friction=0.3),
    # measured pair: 34.3 g/m2, G_P = 0.602 uN m (Fuss 2024, Table 1)
    'silk_georgette': Material('silk_georgette', areal_density=0.0343, bending_rigidity=0.602e-6,
                               shear_stiffness=0.3 * GF_CM_PER_DEG, tensile_stiffness=5000.0, friction=0.3),
    # bracket: 0.06-0.07 kg/m2 (12-16 momme), B 1-3 uN m; central values chosen
    'silk_crepe_de_chine': Material('silk_crepe_de_chine', areal_density=0.065, bending_rigidity=1.5e-6,
                                    shear_stiffness=0.4 * GF_CM_PER_DEG, tensile_stiffness=5000.0, friction=0.3),
    # bracket: 0.08-0.10 kg/m2 (19-23 momme charmeuse to light duchess), B 2-6 uN m; central values chosen
    'silk_satin_heavy': Material('silk_satin_heavy', areal_density=0.090, bending_rigidity=3.0e-6,
                                 shear_stiffness=0.5 * GF_CM_PER_DEG, tensile_stiffness=5000.0, friction=0.25),
    # measured pair: 199.5 g/m2, KES-F2 B 14.3 (warp) and 9.00 (weft) uN m (Sun et al. 2017, sample 14)
    'wool_plain': Material('wool_plain', areal_density=0.1995, bending_rigidity=11.65e-6,
                           shear_stiffness=0.8 * GF_CM_PER_DEG, tensile_stiffness=5000.0, friction=0.35),
    # measured pair: 91.9 g/m2, KES-F2 B 7.30 (warp) and 3.83 (weft) uN m (Sun et al. 2017, sample 1)
    'cotton_plain': Material('cotton_plain', areal_density=0.0919, bending_rigidity=5.565e-6,
                             shear_stiffness=0.8 * GF_CM_PER_DEG, tensile_stiffness=5000.0, friction=0.35),
    # measured weight and drape coefficient (Lee et al. 2011: 69 g/m2, 21.9%); bending inferred by matching the
    # solver's drape test to that drape (trials.py, satin grid): 0.80 uN m
    'silk_satin': Material('silk_satin', areal_density=0.069, bending_rigidity=0.80e-6,
                           shear_stiffness=0.4 * GF_CM_PER_DEG, tensile_stiffness=5000.0, friction=0.25),
    # supplier weight (6 momme); bending chosen: a crisp sized sheer, bending length 3.0 cm on Earth
    'silk_organza': Material('silk_organza', areal_density=0.026, bending_rigidity=6.9e-6,
                             shear_stiffness=0.5 * GF_CM_PER_DEG, tensile_stiffness=5000.0, friction=0.3),
    # supplier weight (8 momme); bending chosen: a soft plain weave, bending length 1.3 cm on Earth
    'silk_habotai': Material('silk_habotai', areal_density=0.035, bending_rigidity=0.75e-6,
                             shear_stiffness=0.3 * GF_CM_PER_DEG, tensile_stiffness=5000.0, friction=0.3),
    # supplier bracket 14-19 momme; bending chosen: a crisp taffeta, bending length 3.5 cm on Earth
    'silk_taffeta': Material('silk_taffeta', areal_density=0.070, bending_rigidity=29.4e-6,
                             shear_stiffness=0.8 * GF_CM_PER_DEG, tensile_stiffness=5000.0, friction=0.3),
    # supplier weight; bending chosen: a heavy pile that drapes, bending length 1.7 cm on Earth
    'silk_velvet': Material('silk_velvet', areal_density=0.170, bending_rigidity=8.2e-6,
                            shear_stiffness=0.6 * GF_CM_PER_DEG, tensile_stiffness=5000.0, friction=0.4),
    # supplier weights 14-19 g/m2; bending chosen: a sized net, bending length 2.0 cm on Earth; open weave takes 0.6
    # of a closed sheet's normal pressure (assumed open area about half; porous-plate drag, Steiros & Hultmark 2018)
    'nylon_tulle': Material('nylon_tulle', areal_density=0.015, bending_rigidity=1.18e-6,
                            shear_stiffness=0.2 * GF_CM_PER_DEG, tensile_stiffness=1000.0, friction=0.3,
                            normal_factor=0.6),
    # supplier bracket 140-180 g/m2; bending chosen: a limp knit, bending length 1.3 cm on Earth; a knit stretches
    'viscose_jersey': Material('viscose_jersey', areal_density=0.160, bending_rigidity=3.45e-6,
                               shear_stiffness=0.3 * GF_CM_PER_DEG, tensile_stiffness=200.0, friction=0.35),
}

# How each fabric's bending rigidity stands: a measurement, inferred from a measured drape coefficient through the
# solver's drape test, a stated bracket, or a modelling choice (sources.json).
BENDING_STATUS = {'silk_chiffon': 'measured', 'silk_georgette': 'measured', 'wool_plain': 'measured',
                  'cotton_plain': 'measured', 'silk_satin': 'inferred', 'silk_crepe_de_chine': 'bracket',
                  'silk_satin_heavy': 'bracket',
                  'silk_organza': 'chosen', 'silk_habotai': 'chosen', 'silk_taffeta': 'chosen',
                  'silk_velvet': 'chosen', 'nylon_tulle': 'chosen', 'viscose_jersey': 'chosen'}


def fabric(name: str) -> Material:
    return FABRICS[name]


def bending_length(m: Material, g: float) -> float:
    """Peirce's bending length c = (B / w)^(1/3), w the fabric's weight per unit area."""
    return (m.bending_rigidity / (m.areal_density * g)) ** (1.0 / 3.0)


def sources() -> dict:
    return json.loads((Path(__file__).resolve().parent / 'sources.json').read_text())
