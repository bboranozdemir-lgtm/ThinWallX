"""ThinWallX — Numerical analysis tool for arbitrary thin-walled structural sections."""

from __future__ import annotations

from thinwallx.cells import Cell, CellTopology, extract_cell_topology
from thinwallx.closed_section import ClosedSection
from thinwallx.closed_shear_flow import (
    ClosedSegmentShearFlow,
    ClosedShearFlowResult,
    calculate_closed_shear_flow,
)
from thinwallx.closed_torsion import (
    ClosedTorsionWarpingResult,
    compute_closed_shear_center,
    compute_closed_torsion_warping,
)
from thinwallx.exceptions import (
    GeometryError,
    SingularSectionError,
    ThinWallXError,
    TopologyError,
    ValidationError,
)
from thinwallx.mixed_section import MixedSection
from thinwallx.mixed_shear_flow import (
    MixedSegmentShearFlow,
    MixedShearFlowResult,
    calculate_mixed_shear_flow,
)
from thinwallx.mixed_topology import MixedTopology, extract_mixed_topology
from thinwallx.mixed_torsion import (
    MixedTorsionWarpingResult,
    compute_mixed_shear_center,
    compute_mixed_torsion_warping,
)
from thinwallx.primitives import Node, Segment
from thinwallx.properties import SectionProperties, compute_properties
from thinwallx.section import Section
from thinwallx.shear_center import (
    ShearCenterResult,
    compute_shear_center,
    compute_shear_flow_torque,
)
from thinwallx.shear_flow import (
    SegmentShearFlow,
    ShearFlowResult,
    calculate_shear_flow,
)
from thinwallx.shear_load import ShearLoad
from thinwallx.torsion import (
    SegmentWarping,
    TorsionWarpingResult,
    compute_torsion_warping,
)
from thinwallx.validation import (
    cluster_nodes,
    validate_section_geometry_and_topology,
)

from thinwallx.stress import AppliedLoads, SegmentStressProfile, StressRecoveryResult, calculate_stresses

__version__ = "1.0.0"

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
    "ThinWallXError",
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
