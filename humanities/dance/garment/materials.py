"""Fabrics for the dance garments. Every number and its basis is recorded in sources.json beside this file.

Areal density and bending rigidity come in measured pairs where a measurement was found (Fuss 2024: silk chiffon
and silk georgette, cantilever bending to ASTM D1388 on one fabric each). Crepe de chine and the heavier satin rest on
stated brackets; their central values are modelling choices inside those brackets. Shear stiffness is a bracket from
the usual KES-F range for light silks (G about 0.2-0.6 gf/cm/deg); tensile stiffness only has to make the cloth
effectively inextensible at the loads met here (strains below 1e-3), which it does.
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
}


def fabric(name: str) -> Material:
    return FABRICS[name]


def bending_length(m: Material, g: float) -> float:
    """Peirce's bending length c = (B / w)^(1/3), w the fabric's weight per unit area."""
    return (m.bending_rigidity / (m.areal_density * g)) ** (1.0 / 3.0)


def sources() -> dict:
    return json.loads((Path(__file__).resolve().parent / 'sources.json').read_text())
