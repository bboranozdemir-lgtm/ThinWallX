"""MixedSection class representing mixed open-closed thin-walled cross sections (Sectalix v0.6)."""

from __future__ import annotations

import math
from typing import Sequence

import numpy as np

from sectalix.cells import Cell
from sectalix.mixed_shear_flow import MixedShearFlowResult, calculate_mixed_shear_flow
from sectalix.mixed_topology import MixedTopology, extract_mixed_topology
from sectalix.mixed_torsion import (
    MixedTorsionWarpingResult,
    compute_mixed_shear_center,
    compute_mixed_torsion_warping,
)
from sectalix.primitives import Node, Segment
from sectalix.properties import SectionProperties, compute_properties
from sectalix.section import Section
from sectalix.shear_center import ShearCenterResult
from sectalix.shear_load import ShearLoad


class MixedSection(Section):
    """A mixed open-closed thin-walled cross section composed of closed cells and open branches.

    Attributes:
        segments: Tuple of straight centerline segments.
        canonical_nodes: Tuple of clustered representative Node objects.
        properties: Computed SectionProperties.
        mixed_topology: Extracted MixedTopology container.
    """

    def __init__(
        self,
        segments: Sequence[Segment],
        validate: bool = True,
        node_tolerance: float = 1e-9,
        safety_factor: float = 1e4,
    ) -> None:
        """Initialize and validate a MixedSection.

        Args:
            segments: Sequence of Segment objects.
            validate: Whether to run mixed topology validation (default True).
            node_tolerance: Absolute Euclidean coordinate tolerance for merging nodes.
            safety_factor: Safety multiplier for positive-definiteness checks.
        """
        seg_tuple = tuple(segments)
        self._segments = seg_tuple
        self._node_tolerance = node_tolerance
        self._safety_factor = safety_factor
        self._is_validated = validate

        if validate:
            self._mixed_topology = extract_mixed_topology(
                seg_tuple,
                node_tolerance=node_tolerance,
                safety_factor=safety_factor,
            )
            self._canonical_nodes = self._mixed_topology.canonical_nodes
        else:
            self._mixed_topology = None  # type: ignore[assignment]
            all_nodes = [node for seg in seg_tuple for node in (seg.p1, seg.p2)]
            from sectalix.validation import cluster_nodes
            canonical_nodes, _ = cluster_nodes(all_nodes, tol=node_tolerance)
            self._canonical_nodes = tuple(canonical_nodes)

        self._properties: SectionProperties = compute_properties(self._segments)
        self._torsion_warping_result: MixedTorsionWarpingResult | None = None
        self._torsion_constants_result: tuple[float, float, float] | None = None

    def validate(self) -> None:
        """Run mixed topology validation on this section.

        Raises:
            GeometryError: If any geometric property is invalid.
            TopologyError: If topology is empty, disconnected, or invalid.
            SingularSectionError: If cell flexibility matrix H is singular.
        """
        self._mixed_topology = extract_mixed_topology(
            self._segments,
            node_tolerance=self._node_tolerance,
            safety_factor=self._safety_factor,
        )
        self._canonical_nodes = self._mixed_topology.canonical_nodes
        self._is_validated = True

    @property
    def mixed_topology(self) -> MixedTopology:
        """Mixed topology container."""
        if self._mixed_topology is None:
            self.validate()
        assert self._mixed_topology is not None
        return self._mixed_topology

    @property
    def cell_topology(self) -> MixedTopology:
        """Alias for mixed_topology to support cellular interfaces."""
        return self.mixed_topology

    @property
    def topology(self) -> str:
        """Topology classification: 'open', 'closed', or 'mixed'."""
        return self.mixed_topology.topology_type

    @property
    def topology_type(self) -> str:
        """Topology classification: 'open', 'closed', or 'mixed'."""
        return self.mixed_topology.topology_type

    @property
    def is_mixed(self) -> bool:
        """True if section contains both closed cells and open branches."""
        return self.mixed_topology.is_mixed

    @property
    def is_closed(self) -> bool:
        """True if section contains closed cells and zero open branches."""
        return self.mixed_topology.is_closed

    @property
    def is_open(self) -> bool:
        """True if section contains zero closed cells."""
        return self.mixed_topology.is_open

    @property
    def cells(self) -> tuple[Cell, ...]:
        """Extracted bounded counter-clockwise cells."""
        return self.mixed_topology.cells

    @property
    def cell_count(self) -> int:
        """Number of bounded cells n_c."""
        return self.mixed_topology.cell_count

    @property
    def B(self) -> np.ndarray:
        """Oriented cell-wall incidence matrix of shape (n_c, E)."""
        return self.mixed_topology.B

    @property
    def H(self) -> np.ndarray:
        """Physical cell flexibility matrix H of shape (n_c, n_c).

        Raises:
            OverflowError: If H elements exceed IEEE-754 float64 range.
        """
        return self.mixed_topology.H

    @property
    def H_scaled(self) -> np.ndarray:
        """Dimensionless cell flexibility matrix of shape (n_c, n_c)."""
        return self.mixed_topology.H_scaled

    @property
    def e_max_H(self) -> int:
        """Binary scaling exponent of H."""
        return self.mixed_topology.e_max_H

    @property
    def block_cell_indices(self) -> tuple[tuple[int, ...], ...]:
        """Global cell indices grouped by independent cyclic block."""
        return self.mixed_topology.block_cell_indices

    @property
    def H_block_scaled(self) -> tuple[np.ndarray, ...]:
        """Independently scaled flexibility submatrix for each cyclic block."""
        return self.mixed_topology.H_block_scaled

    @property
    def e_max_H_block(self) -> tuple[int, ...]:
        """Binary exponent paired with each independently scaled H block."""
        return self.mixed_topology.e_max_H_block

    @property
    def closed_edges(self) -> tuple[int, ...]:
        """Segment indices of closed cyclic walls E_c."""
        return self.mixed_topology.closed_edges

    @property
    def open_edges(self) -> tuple[int, ...]:
        """Segment indices of open bridge branches E_o."""
        return self.mixed_topology.open_edges

    @property
    def junction_nodes(self) -> tuple[int, ...]:
        """Canonical node indices connecting closed cells and open branches."""
        return self.mixed_topology.junction_nodes

    @property
    def articulation_nodes(self) -> tuple[int, ...]:
        """Canonical node indices that are cut vertices."""
        return self.mixed_topology.articulation_nodes

    @property
    def free_tip_nodes(self) -> tuple[int, ...]:
        """Canonical node indices of degree 1 in the full graph."""
        return self.mixed_topology.free_tip_nodes

    @property
    def cyclic_components(self) -> tuple[tuple[int, ...], ...]:
        """Connected components of the cyclic subgraph G_c."""
        return self.mixed_topology.cyclic_components

    @property
    def cyclic_blocks(self) -> tuple[tuple[int, ...], ...]:
        """Biconnected blocks of the cyclic subgraph G_c."""
        return self.mixed_topology.cyclic_blocks

    @property
    def open_components(self) -> tuple[tuple[int, ...], ...]:
        """Connected tree components of the open bridge forest E_o."""
        return self.mixed_topology.open_components

    @property
    def cut_edges(self) -> tuple[int, ...]:
        """Selected chord segment indices strictly from E_c."""
        return self.mixed_topology.cut_edges

    def calculate_shear_flow(
        self,
        vx: float | ShearLoad,
        vy: float | None = None,
        cut_param: float = 0.5,
        safety_factor: float = 1e4,
        spanning_tree_edges: Sequence[int] | None = None,
    ) -> MixedShearFlowResult:
        """Compute exact mixed-section transverse shear flow.

        Args:
            vx: Transverse shear force in x, or ShearLoad object.
            vy: Transverse shear force in y.
            cut_param: Parameter xi_cut in (0, 1) for chord cuts (default 0.5).
            safety_factor: Multiplier for positive-definiteness checks.
            spanning_tree_edges: Optional sequence of segment indices for spanning tree.

        Returns:
            MixedShearFlowResult object.
        """
        return calculate_mixed_shear_flow(
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
        """Compute exact shear-center location for this mixed section.

        Args:
            safety_factor: Multiplier for positive-definiteness checks.

        Returns:
            ShearCenterResult with offsets (ex, ey) and coordinates (x_s, y_s).
        """
        return compute_mixed_shear_center(self, safety_factor=safety_factor)

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
    ) -> MixedTorsionWarpingResult:
        """Compute exact mixed-section hybrid torsion and warping properties.

        Args:
            root_node_idx: Canonical node index used as root for raw omega propagation.
            safety_factor: Multiplier for positive-definiteness of H.

        Returns:
            MixedTorsionWarpingResult containing J_total, J_BB, J_open, C_w, and warping fields.
        """
        res = compute_mixed_torsion_warping(
            self, root_node_idx=root_node_idx, safety_factor=safety_factor
        )
        self._torsion_warping_result = res
        self._torsion_constants_result = (res.J_total, res.J_BB, res.J_open)
        return res

    def compute_torsion_warping(
        self,
        root_node_idx: int = 0,
        safety_factor: float = 1e4,
    ) -> MixedTorsionWarpingResult:
        """Alias for torsion_properties."""
        return self.torsion_properties(
            root_node_idx=root_node_idx, safety_factor=safety_factor
        )

    def torsion_constants(
        self,
        safety_factor: float = 1e4,
    ) -> tuple[float, float, float]:
        """Compute Saint-Venant torsion constants (J_total, J_BB, J_open) without shear center or warping.

        Returns:
            tuple of (J_total, J_BB, J_open)
        """
        if self._torsion_constants_result is None:
            from sectalix.mixed_torsion import compute_mixed_torsion_constant
            self._torsion_constants_result = compute_mixed_torsion_constant(
                self, safety_factor=safety_factor
            )
        return self._torsion_constants_result

    @property
    def J(self) -> float:
        """Total hybrid Saint-Venant torsion constant J_total = J_BB + J_open."""
        if self._torsion_warping_result is not None:
            return self._torsion_warping_result.J
        return self.torsion_constants()[0]

    @property
    def J_total(self) -> float:
        """Total hybrid Saint-Venant torsion constant J_total = J_BB + J_open."""
        if self._torsion_warping_result is not None:
            return self._torsion_warping_result.J_total
        return self.torsion_constants()[0]

    @property
    def J_BB(self) -> float:
        """Bredt-Batho closed-cell torsion constant J_BB."""
        if self._torsion_warping_result is not None:
            return self._torsion_warping_result.J_BB
        return self.torsion_constants()[1]

    @property
    def J_open(self) -> float:
        """Open-branch strip torsion constant J_open = sum_{E_o} (L*t^3/3)."""
        if self._torsion_warping_result is not None:
            return self._torsion_warping_result.J_open
        return self.torsion_constants()[2]

    @property
    def Cw(self) -> float:
        """Mixed-section warping constant C_w."""
        return self.torsion_properties().Cw

    @property
    def sectorial_coordinates(self) -> tuple[float, ...]:
        """Normalized principal sectorial coordinates omega at section canonical nodes."""
        return self.torsion_properties().node_omega

    def translated(self, dx: float, dy: float) -> MixedSection:
        """Return a new MixedSection translated by (dx, dy)."""
        translated_segments = [
            Segment(
                p1=Node(x=seg.p1.x + dx, y=seg.p1.y + dy, id=seg.p1.id),
                p2=Node(x=seg.p2.x + dx, y=seg.p2.y + dy, id=seg.p2.id),
                t=seg.t,
                id=seg.id,
            )
            for seg in self._segments
        ]
        return MixedSection(
            translated_segments,
            validate=False,
            node_tolerance=self._node_tolerance,
            safety_factor=self._safety_factor,
        )

    def rotated(
        self, angle_rad: float, origin: tuple[float, float] = (0.0, 0.0)
    ) -> MixedSection:
        """Return a new MixedSection rotated counter-clockwise by angle_rad around origin."""
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
        return MixedSection(
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
    ) -> MixedSection:
        """Build MixedSection from a sequence of Segment objects."""
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
    ) -> MixedSection:
        """Build MixedSection from a sequence of ((x1, y1), (x2, y2), thickness) tuples."""
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
