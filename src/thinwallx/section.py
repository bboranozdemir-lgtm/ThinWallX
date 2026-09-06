"""Section class representing an open thin-walled cross section."""

from __future__ import annotations

import math
from typing import Sequence, TYPE_CHECKING

import numpy as np

from thinwallx.primitives import Node, Segment
from thinwallx.properties import SectionProperties, compute_properties
from thinwallx.shear_center import ShearCenterResult, compute_shear_center
from thinwallx.shear_flow import ShearFlowResult, calculate_shear_flow
from thinwallx.shear_load import ShearLoad
from thinwallx.torsion import (
    SegmentWarping,
    TorsionWarpingResult,
    compute_torsion_warping,
)
from thinwallx.validation import validate_section_geometry_and_topology

if TYPE_CHECKING:
    from thinwallx.stress import AppliedLoads, StressRecoveryResult



class Section:
    """An open thin-walled cross section composed of centerline segments.

    Attributes:
        segments: Tuple of straight centerline segments.
        canonical_nodes: Tuple of unique clustered nodes.
        properties: Computed SectionProperties.
    """

    def calculate_stresses(self, loads: AppliedLoads) -> StressRecoveryResult:
        """Recover v0.7 stresses; inherited by ClosedSection and MixedSection."""
        from thinwallx.stress import calculate_stresses
        return calculate_stresses(self, loads)

    def stresses(self, loads: AppliedLoads) -> StressRecoveryResult:
        """Alias for calculate_stresses."""
        return self.calculate_stresses(loads)

    def __init__(
        self,
        segments: Sequence[Segment],
        validate: bool = True,
        node_tolerance: float = 1e-9,
    ) -> None:
        """Initialize and validate a Section.

        Args:
            segments: Sequence of Segment objects.
            validate: Whether to run geometry and topology validation (default True).
            node_tolerance: Absolute Euclidean tolerance, in coordinate units, for
                merging shared node coordinates.
        """
        seg_tuple = tuple(segments)
        if validate:
            canonical_nodes, _ = validate_section_geometry_and_topology(
                seg_tuple, node_tolerance=node_tolerance
            )
            self._canonical_nodes = tuple(canonical_nodes)
        else:
            self._canonical_nodes = tuple(
                {seg.p1 for seg in seg_tuple} | {seg.p2 for seg in seg_tuple}
            )

        self._segments = seg_tuple
        self._node_tolerance = node_tolerance
        self._is_validated = validate
        self._properties: SectionProperties = compute_properties(self._segments)

    def validate(self) -> None:
        """Run geometry and topology validation on this section.

        Raises:
            GeometryError: If any geometric property is invalid.
            TopologyError: If section topology is invalid (empty, disconnected, cyclical, duplicate).
        """
        canonical_nodes, _ = validate_section_geometry_and_topology(
            self._segments, node_tolerance=self._node_tolerance
        )
        self._canonical_nodes = tuple(canonical_nodes)
        self._is_validated = True

    @property
    def is_closed(self) -> bool:
        """False for open Section."""
        return False

    @property
    def is_open(self) -> bool:
        """True for open Section."""
        return True

    @property
    def is_mixed(self) -> bool:
        """False for open Section."""
        return False

    @property
    def topology(self) -> str:
        """Topology classification: 'open'."""
        return "open"

    @property
    def topology_type(self) -> str:
        """Topology classification: 'open'."""
        return "open"

    @property
    def segments(self) -> tuple[Segment, ...]:
        """Centerline segments comprising the section."""
        return self._segments

    @property
    def nodes(self) -> tuple[Node, ...]:
        """Canonical nodes of the section."""
        return self._canonical_nodes

    @property
    def properties(self) -> SectionProperties:
        """Full computed section properties."""
        return self._properties

    @property
    def area(self) -> float:
        """Total cross-sectional area A."""
        return self._properties.area

    @property
    def centroid(self) -> tuple[float, float]:
        """Centroid coordinates (cx, cy)."""
        return self._properties.centroid

    @property
    def cx(self) -> float:
        """Centroid horizontal coordinate."""
        return self._properties.cx

    @property
    def cy(self) -> float:
        """Centroid vertical coordinate."""
        return self._properties.cy

    @property
    def ix_raw(self) -> float:
        """Raw second moment of area about origin x-axis."""
        return self._properties.ix_raw

    @property
    def iy_raw(self) -> float:
        """Raw second moment of area about origin y-axis."""
        return self._properties.iy_raw

    @property
    def ixy_raw(self) -> float:
        """Raw product moment of area about origin."""
        return self._properties.ixy_raw

    @property
    def Ix(self) -> float:
        """Centroidal second moment of area about horizontal axis."""
        return self._properties.ix

    @property
    def Iy(self) -> float:
        """Centroidal second moment of area about vertical axis."""
        return self._properties.iy

    @property
    def Ixy(self) -> float:
        """Centroidal product moment of area."""
        return self._properties.ixy

    @property
    def I1(self) -> float:
        """Major principal second moment of area (I1 >= I2)."""
        return self._properties.i1

    @property
    def I2(self) -> float:
        """Minor principal second moment of area."""
        return self._properties.i2

    @property
    def theta_p(self) -> float:
        """Principal-axis angle in radians for I1, in (-pi/2, pi/2]."""
        return self._properties.theta_p

    @property
    def theta_p_deg(self) -> float:
        """Principal-axis angle in degrees."""
        return self._properties.theta_p_deg

    @property
    def is_degenerate(self) -> bool:
        """True if the section is isotropic (all in-plane axes are principal)."""
        return self._properties.is_degenerate

    @property
    def inertia_matrix(self) -> np.ndarray:
        """2D inertia tensor matrix: [[Ix, -Ixy], [-Ixy, Iy]]."""
        return self._properties.inertia_matrix

    def translated(self, dx: float, dy: float) -> Section:
        """Return a new Section translated by (dx, dy).

        Args:
            dx: Horizontal translation offset.
            dy: Vertical translation offset.

        Returns:
            New translated Section.
        """
        translated_segments = [
            Segment(
                p1=Node(x=seg.p1.x + dx, y=seg.p1.y + dy, id=seg.p1.id),
                p2=Node(x=seg.p2.x + dx, y=seg.p2.y + dy, id=seg.p2.id),
                t=seg.t,
                id=seg.id,
            )
            for seg in self._segments
        ]
        return Section(translated_segments, validate=False, node_tolerance=self._node_tolerance)

    def rotated(
        self, angle_rad: float, origin: tuple[float, float] = (0.0, 0.0)
    ) -> Section:
        """Return a new Section rotated counter-clockwise by angle_rad around origin.

        Transformation:
            x' = x0 + (x - x0)*cos(phi) - (y - y0)*sin(phi)
            y' = y0 + (x - x0)*sin(phi) + (y - y0)*cos(phi)

        Args:
            angle_rad: Rotation angle in radians (counter-clockwise).
            origin: Center of rotation (x0, y0). Default (0, 0).

        Returns:
            New rotated Section.
        """
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
        return Section(rotated_segments, validate=False, node_tolerance=self._node_tolerance)

    @classmethod
    def from_segments(
        cls,
        segments: Sequence[Segment],
        validate: bool = True,
        node_tolerance: float = 1e-9,
    ) -> Section:
        """Build Section from a sequence of Segment objects."""
        return cls(segments, validate=validate, node_tolerance=node_tolerance)

    @classmethod
    def from_tuples(
        cls,
        segment_tuples: Sequence[tuple[tuple[float, float], tuple[float, float], float]],
        validate: bool = True,
        node_tolerance: float = 1e-9,
    ) -> Section:
        """Build Section from a sequence of ((x1, y1), (x2, y2), thickness) tuples.

        Args:
            segment_tuples: Sequence of ((x1, y1), (x2, y2), t).
            validate: Whether to validate geometry and topology.
            node_tolerance: Absolute tolerance, in coordinate units, for merging
                shared nodes.
        """
        segments = [
            Segment(
                p1=Node(x=pt1[0], y=pt1[1]),
                p2=Node(x=pt2[0], y=pt2[1]),
                t=t,
                id=i,
            )
            for i, (pt1, pt2, t) in enumerate(segment_tuples)
        ]
        return cls(segments, validate=validate, node_tolerance=node_tolerance)

    def calculate_shear_flow(
        self,
        vx: float | ShearLoad,
        vy: float | None = None,
        safety_factor: float = 1e4,
    ) -> ShearFlowResult:
        """Compute exact thin-walled shear-flow distribution for this open section.

        Args:
            vx: Transverse shear force in x direction, or a ShearLoad object.
            vy: Transverse shear force in y direction.
            safety_factor: Multiplier for numerical positive-definiteness check.

        Returns:
            ShearFlowResult object containing exact segment shear-flow fields.
        """
        return calculate_shear_flow(self, vx=vx, vy=vy, safety_factor=safety_factor)

    def compute_shear_center(
        self,
        safety_factor: float = 1e4,
    ) -> ShearCenterResult:
        """Compute exact shear-center location and centroid-relative offsets.

        Args:
            safety_factor: Multiplier for numerical positive-definiteness check.

        Returns:
            ShearCenterResult object with offsets (ex, ey) and coordinates (x_s, y_s).
        """
        return compute_shear_center(self, safety_factor=safety_factor)

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
    ) -> TorsionWarpingResult:
        """Compute open-section Saint-Venant torsion and warping properties.

        Args:
            root_node_idx: Canonical node index used as root for raw propagation (default 0).
            safety_factor: Multiplier for numerical positive-definiteness in shear-center solve.

        Returns:
            TorsionWarpingResult containing J, C_w, normalized omega field, and segment warpings.
        """
        return compute_torsion_warping(
            self, root_node_idx=root_node_idx, safety_factor=safety_factor
        )

    def compute_torsion_warping(
        self,
        root_node_idx: int = 0,
        safety_factor: float = 1e4,
    ) -> TorsionWarpingResult:
        """Alias for torsion_properties."""
        return self.torsion_properties(
            root_node_idx=root_node_idx, safety_factor=safety_factor
        )

    @property
    def J(self) -> float:
        """Saint-Venant open-section torsion constant J."""
        return self.torsion_properties().J

    @property
    def J_total(self) -> float:
        """Total Saint-Venant torsion constant J_total = J_open for open section."""
        return self.J

    @property
    def J_BB(self) -> float:
        """Bredt-Batho closed-cell torsion constant (0.0 for open section)."""
        return 0.0

    @property
    def J_open(self) -> float:
        """Open-strip torsion constant J_open = J for open section."""
        return self.J

    @property
    def Cw(self) -> float:
        """Warping constant C_w."""
        return self.torsion_properties().Cw

    @property
    def sectorial_coordinates(self) -> tuple[float, ...]:
        """Normalized principal sectorial coordinates omega at section canonical nodes."""
        return self.torsion_properties().node_omega

