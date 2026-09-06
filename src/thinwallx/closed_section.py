"""ClosedSection class representing closed single-cell and multi-cell thin-walled cross sections."""

from __future__ import annotations

import math
from typing import Sequence

import numpy as np

from thinwallx.cells import Cell, CellTopology, extract_cell_topology
from thinwallx.closed_shear_flow import ClosedShearFlowResult, calculate_closed_shear_flow
from thinwallx.closed_torsion import (
    ClosedTorsionWarpingResult,
    compute_closed_shear_center,
    compute_closed_torsion_warping,
)
from thinwallx.primitives import Node, Segment
from thinwallx.properties import SectionProperties, compute_properties
from thinwallx.section import Section
from thinwallx.shear_center import ShearCenterResult
from thinwallx.shear_load import ShearLoad


class ClosedSection(Section):
    """A closed thin-walled cross section composed of one or more bounded cells.

    Attributes:
        segments: Tuple of straight centerline segments.
        canonical_nodes: Tuple of clustered representative Node objects.
        properties: Computed SectionProperties.
        cell_topology: Extracted planar CellTopology.
    """

    def __init__(
        self,
        segments: Sequence[Segment],
        validate: bool = True,
        node_tolerance: float = 1e-9,
        safety_factor: float = 1e4,
    ) -> None:
        """Initialize and validate a ClosedSection.

        Args:
            segments: Sequence of Segment objects.
            validate: Whether to run planar cell topology validation (default True).
            node_tolerance: Absolute Euclidean coordinate tolerance for merging nodes.
            safety_factor: Safety multiplier for positive-definiteness checks.
        """
        seg_tuple = tuple(segments)
        self._segments = seg_tuple
        self._node_tolerance = node_tolerance
        self._safety_factor = safety_factor
        self._is_validated = validate

        if validate:
            self._cell_topology = extract_cell_topology(
                seg_tuple,
                node_tolerance=node_tolerance,
                safety_factor=safety_factor,
            )
            self._canonical_nodes = self._cell_topology.canonical_nodes
        else:
            self._cell_topology = None  # type: ignore[assignment]
            all_nodes = [node for seg in seg_tuple for node in (seg.p1, seg.p2)]
            from thinwallx.validation import cluster_nodes
            canonical_nodes, _ = cluster_nodes(all_nodes, tol=node_tolerance)
            self._canonical_nodes = tuple(canonical_nodes)

        self._properties: SectionProperties = compute_properties(self._segments)

    def validate(self) -> None:
        """Run planar cell topology validation on this closed section.

        Raises:
            GeometryError: If any geometric property is invalid or mixed open/closed topology.
            TopologyError: If topology is empty, disconnected, or non-cellular.
            SingularSectionError: If cell flexibility matrix H is singular.
        """
        self._cell_topology = extract_cell_topology(
            self._segments,
            node_tolerance=self._node_tolerance,
            safety_factor=self._safety_factor,
        )
        self._canonical_nodes = self._cell_topology.canonical_nodes
        self._is_validated = True

    @property
    def is_closed(self) -> bool:
        """True for ClosedSection."""
        return True

    @property
    def is_open(self) -> bool:
        """False for ClosedSection."""
        return False

    @property
    def is_mixed(self) -> bool:
        """False for ClosedSection."""
        return False

    @property
    def topology(self) -> str:
        """Topology classification: 'closed'."""
        return "closed"

    @property
    def topology_type(self) -> str:
        """Topology classification: 'closed'."""
        return "closed"

    @property
    def cell_topology(self) -> CellTopology:
        """Planar cell topology container."""
        if self._cell_topology is None:
            self.validate()
        assert self._cell_topology is not None
        return self._cell_topology

    @property
    def cells(self) -> tuple[Cell, ...]:
        """Extracted bounded counter-clockwise cells."""
        return self.cell_topology.cells

    @property
    def cell_count(self) -> int:
        """Number of bounded cells n_c."""
        return self.cell_topology.cell_count

    @property
    def B(self) -> np.ndarray:
        """Oriented cell-wall incidence matrix of shape (n_c, E)."""
        return self.cell_topology.B

    @property
    def H(self) -> np.ndarray:
        """Cell flexibility matrix H = B diag(L/t) B^T of shape (n_c, n_c).

        Raises:
            OverflowError: If H elements exceed IEEE-754 float64 range.
        """
        return self.cell_topology.H

    def calculate_shear_flow(
        self,
        vx: float | ShearLoad,
        vy: float | None = None,
        cut_param: float = 0.5,
        safety_factor: float = 1e4,
        spanning_tree_edges: Sequence[int] | None = None,
    ) -> ClosedShearFlowResult:
        """Compute exact closed-section transverse shear flow.

        Args:
            vx: Transverse shear force in x, or ShearLoad object.
            vy: Transverse shear force in y.
            cut_param: Parameter xi_cut in (0, 1) for chord cuts (default 0.5).
            safety_factor: Multiplier for positive-definiteness checks.
            spanning_tree_edges: Optional sequence of segment indices for spanning tree.

        Returns:
            ClosedShearFlowResult object.
        """
        return calculate_closed_shear_flow(
            self,
            vx=vx,
            vy=vy,
            cut_param=cut_param,
            safety_factor=safety_factor,
            spanning_tree_edges=spanning_tree_edges,
        )

    def compute_shear_center(
        self,
        safety_factor: float = 1e4,
    ) -> ShearCenterResult:
        """Compute exact shear-center location for this closed section.

        Args:
            safety_factor: Multiplier for positive-definiteness checks.

        Returns:
            ShearCenterResult with offsets (ex, ey) and coordinates (x_s, y_s).
        """
        return compute_closed_shear_center(self, safety_factor=safety_factor)

    @property
    def shear_center(self) -> tuple[float, float]:
        """Absolute coordinates of the shear center (xs, ys)."""
        res = self.compute_shear_center()
        return (res.x, res.y)

    @property
    def sc_offset(self) -> tuple[float, float]:
        """Centroid-relative offset of the shear center (ex, ey)."""
        res = self.compute_shear_center()
        return (res.ex, res.ey)

    def torsion_properties(
        self,
        root_node_idx: int = 0,
        safety_factor: float = 1e4,
    ) -> ClosedTorsionWarpingResult:
        """Compute exact closed-section Bredt-Batho torsion and warping properties.

        Args:
            root_node_idx: Canonical node index used as root for raw omega propagation.
            safety_factor: Multiplier for positive-definiteness of H.

        Returns:
            ClosedTorsionWarpingResult containing J_BB, C_w, phi, F, and warping fields.
        """
        return compute_closed_torsion_warping(
            self, root_node_idx=root_node_idx, safety_factor=safety_factor
        )

    def compute_torsion_warping(
        self,
        root_node_idx: int = 0,
        safety_factor: float = 1e4,
    ) -> ClosedTorsionWarpingResult:
        """Alias for torsion_properties."""
        return self.torsion_properties(
            root_node_idx=root_node_idx, safety_factor=safety_factor
        )

    @property
    def J(self) -> float:
        """Bredt-Batho Saint-Venant torsion constant J_BB."""
        return self.torsion_properties().J

    @property
    def J_total(self) -> float:
        """Total Saint-Venant torsion constant J_total = J_BB for closed section."""
        return self.J

    @property
    def J_BB(self) -> float:
        """Bredt-Batho closed-cell torsion constant J_BB = J for closed section."""
        return self.J

    @property
    def J_open(self) -> float:
        """Open-strip torsion constant (0.0 for closed section)."""
        return 0.0

    @property
    def Cw(self) -> float:
        """Closed-section warping constant C_w."""
        return self.torsion_properties().Cw

    @property
    def sectorial_coordinates(self) -> tuple[float, ...]:
        """Normalized principal sectorial coordinates omega at section canonical nodes."""
        return self.torsion_properties().node_omega

    def translated(self, dx: float, dy: float) -> ClosedSection:
        """Return a new ClosedSection translated by (dx, dy)."""
        translated_segments = [
            Segment(
                p1=Node(x=seg.p1.x + dx, y=seg.p1.y + dy, id=seg.p1.id),
                p2=Node(x=seg.p2.x + dx, y=seg.p2.y + dy, id=seg.p2.id),
                t=seg.t,
                id=seg.id,
            )
            for seg in self._segments
        ]
        return ClosedSection(
            translated_segments,
            validate=False,
            node_tolerance=self._node_tolerance,
            safety_factor=self._safety_factor,
        )

    def rotated(
        self, angle_rad: float, origin: tuple[float, float] = (0.0, 0.0)
    ) -> ClosedSection:
        """Return a new ClosedSection rotated counter-clockwise by angle_rad around origin."""
        x0, y0 = origin
        cos_a = math.cos(angle_rad)
        sin_a = math.sin(angle_rad)

        def rotate_point(x: float, y: float) -> tuple[float, float]:
            rx = x - x0
            ry = y - y0
            return (x0 + rx * cos_a - ry * sin_a, y0 + rx * sin_a + ry * cos_a)

        rotated_segments = []
        for seg in self._segments:
            p1_rot = rotate_point(seg.p1.x, seg.p1.y)
            p2_rot = rotate_point(seg.p2.x, seg.p2.y)
            rotated_segments.append(
                Segment(
                    p1=Node(x=p1_rot[0], y=p1_rot[1], id=seg.p1.id),
                    p2=Node(x=p2_rot[0], y=p2_rot[1], id=seg.p2.id),
                    t=seg.t,
                    id=seg.id,
                )
            )
        return ClosedSection(
            rotated_segments,
            validate=False,
            node_tolerance=self._node_tolerance,
            safety_factor=self._safety_factor,
        )

    @classmethod
    def from_segments(
        cls,
        segments: Sequence[Segment],
        validate: bool = True,
        node_tolerance: float = 1e-9,
        safety_factor: float = 1e4,
    ) -> ClosedSection:
        """Build ClosedSection from a sequence of Segment objects."""
        return cls(
            segments,
            validate=validate,
            node_tolerance=node_tolerance,
            safety_factor=safety_factor,
        )

    @classmethod
    def from_tuples(
        cls,
        segment_tuples: Sequence[tuple[tuple[float, float], tuple[float, float], float]],
        validate: bool = True,
        node_tolerance: float = 1e-9,
        safety_factor: float = 1e4,
    ) -> ClosedSection:
        """Build ClosedSection from a sequence of ((x1, y1), (x2, y2), thickness) tuples."""
        segments = [
            Segment(
                p1=Node(x=pt1[0], y=pt1[1]),
                p2=Node(x=pt2[0], y=pt2[1]),
                t=t,
                id=i,
            )
            for i, (pt1, pt2, t) in enumerate(segment_tuples)
        ]
        return cls(
            segments,
            validate=validate,
            node_tolerance=node_tolerance,
            safety_factor=safety_factor,
        )
