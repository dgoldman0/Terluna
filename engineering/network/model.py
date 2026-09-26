"""The industrial network model: the structure of the source-to-use ledger.

A skeleton. It fixes the vocabulary the engineering domain's network model uses
(node kinds, commodities, link properties and node constraints) and checks that
a network built from them is consistent, including the author's traffic-safety
rule that passive failures miss planets. The vocabulary follows section 10 of
the reconstructed 19 September architecture
(engineering/reference/industrial_architecture/RECOVERED_INDUSTRIAL_ARCHITECTURE.md,
a reconstruction written on 2026-09-26); it supplies names and relationships,
and no numbers are taken from it. Demand enters from data products such as the
protection module catalogue (ledger.py).
"""
from __future__ import annotations
from dataclasses import dataclass, field

COMMODITIES = {
    'nitrogen_feedstock': 'nitrogen-bearing feedstock and product',
    'water_oxygen_hydrogen': 'water, oxygen and hydrogen',
    'carbon_nutrients': 'carbon and nutrients',
    'structural_materials': 'structural metals, glass, silicates and oxides such as titania',
    'conductors': 'conductor and superconductor material',
    'reactor_fuels': 'reactor fuels and breeding materials',
    'propellant': 'propellant and reaction mass',
    'precision_components': 'precision components, electronics and power hardware',
    'replacement_hardware': 'finished replacement units',
    'waste_coproducts': 'rejected waste and coproducts',
}

NODE_KINDS = {
    'source_extraction': 'source mines and refineries',
    'source_power': 'source power plants and radiators',
    'factory': 'factory complexes',
    'launch': 'launch and injection systems',
    'relay_depot': 'relay depots and tug bases',
    'cislunar_capture': 'cislunar capture and braking complexes',
    'lunar_orbital_depot': 'lunar orbital depots',
    'construction': 'surface and orbital construction industry',
    'protection_service': 'protection-system factories and service depots',
    'storage_release': 'environmental storage and release systems',
    'recycling': 'recycling and reclamation plants',
}

LINK_PROPERTIES = (
    'mass_throughput_kg_s', 'transit_time_s', 'packet_kg', 'packets_per_s', 'container_mass_fraction',
    'propulsion_energy_j_kg', 'propellant_kg_per_kg', 'loss_fraction', 'capture_energy_j_kg',
    'momentum_destination', 'safety_corridor', 'passive_failure_outcome', 'buffer_kg')

NODE_CONSTRAINTS = (
    'power_w', 'conversion_efficiency', 'heat_rejection_temperature_k', 'radiator_area_m2',
    'manufacturing_capacity_kg_s', 'critical_imports', 'maintenance_demand_kg_s', 'fault_recovery_s',
    'discharge_limits')

# research/decisions.md, transport and safety: passive failures miss planets, and nothing
# reaches an inhabited body before verified capture.
SAFE_PASSIVE_FAILURE_OUTCOMES = ('misses_planets', 'returns_to_source', 'held_at_depot')


@dataclass
class Node:
    id: str
    kind: str
    location: str
    capacities_kg_s: dict = field(default_factory=dict)   # commodity -> kg/s
    constraints: dict = field(default_factory=dict)       # keys from NODE_CONSTRAINTS


@dataclass
class Link:
    source: str
    target: str
    commodity: str
    properties: dict = field(default_factory=dict)        # keys from LINK_PROPERTIES


@dataclass
class Network:
    nodes: list
    links: list

    def validate(self):
        """Raise ValueError on the first inconsistency; return self."""
        ids = [n.id for n in self.nodes]
        if len(ids) != len(set(ids)):
            raise ValueError('Node ids must be unique')
        for n in self.nodes:
            if n.kind not in NODE_KINDS:
                raise ValueError(f'Unknown node kind: {n.kind}')
            for c in n.capacities_kg_s:
                if c not in COMMODITIES:
                    raise ValueError(f'Unknown commodity at {n.id}: {c}')
            for k in n.constraints:
                if k not in NODE_CONSTRAINTS:
                    raise ValueError(f'Unknown constraint at {n.id}: {k}')
        known = set(ids)
        for link in self.links:
            if link.source not in known or link.target not in known:
                raise ValueError(f'Link endpoints must be nodes: {link.source} -> {link.target}')
            if link.commodity not in COMMODITIES:
                raise ValueError(f'Unknown commodity on link: {link.commodity}')
            for k in link.properties:
                if k not in LINK_PROPERTIES:
                    raise ValueError(f'Unknown link property: {k}')
            outcome = link.properties.get('passive_failure_outcome')
            if outcome is not None and outcome not in SAFE_PASSIVE_FAILURE_OUTCOMES:
                raise ValueError(f'Passive failure must miss planets ({link.source} -> {link.target}: {outcome})')
        return self
