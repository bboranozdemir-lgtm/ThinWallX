"""Sectalix — Numerical analysis tool for arbitrary thin-walled structural sections."""

from __future__ import annotations

from sectalix.cells import Cell, CellTopology, extract_cell_topology
from sectalix.closed_section import ClosedSection
from sectalix.closed_shear_flow import (
    ClosedSegmentShearFlow,
    ClosedShearFlowResult,
    calculate_closed_shear_flow,
)
from sectalix.closed_torsion import (
    ClosedTorsionWarpingResult,
    compute_closed_shear_center,
    compute_closed_torsion_warping,
)
from sectalix.exceptions import (
    GeometryError,
    SingularSectionError,
    SectalixError,
    TopologyError,
    ValidationError,
)
from sectalix.mixed_section import MixedSection
from sectalix.mixed_shear_flow import (
    MixedSegmentShearFlow,
    MixedShearFlowResult,
    calculate_mixed_shear_flow,
)
from sectalix.mixed_topology import MixedTopology, extract_mixed_topology
from sectalix.mixed_torsion import (
    MixedTorsionWarpingResult,
    compute_mixed_shear_center,
    compute_mixed_torsion_warping,
)
from sectalix.primitives import Node, Segment
from sectalix.properties import SectionProperties, compute_properties
from sectalix.section import Section
from sectalix.shear_center import (
    ShearCenterResult,
    compute_shear_center,
    compute_shear_flow_torque,
)
from sectalix.shear_flow import (
    SegmentShearFlow,
    ShearFlowResult,
    calculate_shear_flow,
)
from sectalix.shear_load import ShearLoad
from sectalix.torsion import (
    SegmentWarping,
    TorsionWarpingResult,
    compute_torsion_warping,
)
from sectalix.validation import (
    cluster_nodes,
    validate_section_geometry_and_topology,
)

from sectalix.stress import AppliedLoads, SegmentStressProfile, StressRecoveryResult, calculate_stresses

__version__ = "1.0.1"

__all__ = [
    "AppliedLoads",
    "SegmentStressProfile",
    "StressRecoveryResult",
    "calculate_stresses",
    "Cell",
    "CellTopology",
    "ClosedSection",
    "ClosedSegmentShearFlow",
    "ClosedShearFlowResult",
    "ClosedTorsionWarpingResult",
    "GeometryError",
    "MixedSection",
    "MixedSegmentShearFlow",
    "MixedShearFlowResult",
    "MixedTopology",
    "MixedTorsionWarpingResult",
    "Node",
    "Section",
    "SectionProperties",
    "Segment",
    "SegmentShearFlow",
    "SegmentWarping",
    "ShearCenterResult",
    "ShearFlowResult",
    "ShearLoad",
    "SingularSectionError",
    "SectalixError",
    "TopologyError",
    "TorsionWarpingResult",
    "ValidationError",
    "calculate_closed_shear_flow",
    "calculate_mixed_shear_flow",
    "calculate_shear_flow",
    "cluster_nodes",
    "compute_closed_shear_center",
    "compute_closed_torsion_warping",
    "compute_mixed_shear_center",
    "compute_mixed_torsion_warping",
    "compute_properties",
    "compute_shear_center",
    "compute_shear_flow_torque",
    "compute_torsion_warping",
    "extract_cell_topology",
    "extract_mixed_topology",
    "validate_section_geometry_and_topology",
]
