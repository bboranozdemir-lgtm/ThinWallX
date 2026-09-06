"""Geometric primitives: Node and Segment dataclasses.

All line integrals implemented herein follow the exact closed-form straight-segment
formulation specified in ThinWallX ACTIVE_PHASE.md.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

from thinwallx.exceptions import GeometryError


@dataclass(frozen=True)
class Node:
    """A 2D point representing a centerline node or vertex.

    Attributes:
        x: Horizontal coordinate.
        y: Vertical coordinate.
        id: Optional identifier or label for the node.
    """

    x: float
    y: float
    id: Any = None

    def __post_init__(self) -> None:
        """Validate coordinates for finiteness."""
        if not math.isfinite(self.x) or not math.isfinite(self.y):
            raise GeometryError(
                f"Node coordinates must be finite. Got x={self.x}, y={self.y} (id={self.id})."
            )

    @property
    def coords(self) -> tuple[float, float]:
        """Return (x, y) coordinates as a tuple."""
        return (self.x, self.y)

    def distance_to(self, other: Node) -> float:
        """Compute Euclidean distance to another Node."""
        return math.hypot(self.x - other.x, self.y - other.y)

    def is_close(self, other: Node, tol: float = 1e-9) -> bool:
        """Check if two nodes coincide within Euclidean tolerance."""
        return self.distance_to(other) <= tol


@dataclass(frozen=True)
class Segment:
    """A straight thin-walled centerline segment with constant thickness.

    Attributes:
        p1: Start node of the centerline segment.
        p2: End node of the centerline segment.
        t: Wall thickness (> 0).
        id: Optional identifier or label.
    """

    p1: Node
    p2: Node
    t: float
    id: Any = None

    def __post_init__(self) -> None:
        """Validate thickness and length."""
        if not math.isfinite(self.t):
            raise GeometryError(f"Segment thickness must be finite. Got t={self.t} (id={self.id}).")
        if self.t <= 0.0:
            raise GeometryError(
                f"Segment thickness must be strictly positive (t > 0). Got t={self.t} (id={self.id})."
            )
        # Reject only an exactly collapsed float representation.  A dimensional
        # absolute cutoff would make otherwise valid geometry depend on the
        # caller's choice of units; topological near-coincidence is handled by
        # Section.node_tolerance instead.
        length = self.length
        if not math.isfinite(length):
            raise GeometryError(
                f"Segment length must be finite. Got length={length} between "
                f"p1=({self.p1.x}, {self.p1.y}) and p2=({self.p2.x}, {self.p2.y}) "
                f"(id={self.id})."
            )
        if length == 0.0:
            raise GeometryError(
                f"Segment length must be strictly positive. Got length={length:.6e} between "
                f"p1=({self.p1.x}, {self.p1.y}) and p2=({self.p2.x}, {self.p2.y}) (id={self.id})."
            )

    @property
    def length(self) -> float:
        """Centerline length: L = sqrt((x2 - x1)^2 + (y2 - y1)^2)."""
        return self.p1.distance_to(self.p2)

    @property
    def area(self) -> float:
        """Wall area contribution: A_i = t * L."""
        return self.t * self.length

    def int_x(self) -> float:
        """Exact line integral: \\int_L x ds = L * (x1 + x2) / 2.

        Derivation:
            Parameterize line s in [0, L]: x(s) = x1 + (x2 - x1)*s/L.
            \\int_0^L x(s) ds = x1*L + (x2 - x1)*L/2 = L*(x1 + x2)/2.
            Source: ThinWallX ACTIVE_PHASE.md, Exact Segment Integrals.
        """
        return self.length * (self.p1.x + self.p2.x) / 2.0

    def int_y(self) -> float:
        """Exact line integral: \\int_L y ds = L * (y1 + y2) / 2.

        Derivation:
            Parameterize line s in [0, L]: y(s) = y1 + (y2 - y1)*s/L.
            \\int_0^L y(s) ds = y1*L + (y2 - y1)*L/2 = L*(y1 + y2)/2.
            Source: ThinWallX ACTIVE_PHASE.md, Exact Segment Integrals.
        """
        return self.length * (self.p1.y + self.p2.y) / 2.0

    def int_x2(self) -> float:
        """Exact line integral: \\int_L x^2 ds = (L / 3) * (x1^2 + x1*x2 + x2^2).

        Derivation:
            With x(s) = x1*(1 - s/L) + x2*(s/L):
            \\int_0^L (1 - s/L)^2 ds = L/3,
            \\int_0^L 2*(1 - s/L)*(s/L) ds = L/3,
            \\int_0^L (s/L)^2 ds = L/3.
            Expanding (x1*(1-s/L) + x2*(s/L))^2 yields:
            (L/3)*(x1^2 + x1*x2 + x2^2).
            Source: ThinWallX ACTIVE_PHASE.md, Exact Segment Integrals.
        """
        x1, x2 = self.p1.x, self.p2.x
        return (self.length / 3.0) * (x1 * x1 + x1 * x2 + x2 * x2)

    def int_y2(self) -> float:
        """Exact line integral: \\int_L y^2 ds = (L / 3) * (y1^2 + y1*y2 + y2^2).

        Derivation:
            Analogous to int_x2 with y coordinates.
            Source: ThinWallX ACTIVE_PHASE.md, Exact Segment Integrals.
        """
        y1, y2 = self.p1.y, self.p2.y
        return (self.length / 3.0) * (y1 * y1 + y1 * y2 + y2 * y2)

    def int_xy(self) -> float:
        """Exact line integral: \\int_L xy ds = (L / 6) * (2*x1*y1 + x1*y2 + x2*y1 + 2*x2*y2).

        Derivation:
            With x(s) = x1*(1 - s/L) + x2*(s/L) and y(s) = y1*(1 - s/L) + y2*(s/L):
            \\int_0^L (1 - s/L)^2 ds = L/3,
            \\int_0^L (1 - s/L)*(s/L) ds = L/6,
            \\int_0^L (s/L)^2 ds = L/3.
            Combining terms:
            x1*y1*(L/3) + (x1*y2 + x2*y1)*(L/6) + x2*y2*(L/3)
            = (L/6) * (2*x1*y1 + x1*y2 + x2*y1 + 2*x2*y2).
            Source: ThinWallX ACTIVE_PHASE.md, Exact Segment Integrals.
        """
        x1, y1 = self.p1.x, self.p1.y
        x2, y2 = self.p2.x, self.p2.y
        return (self.length / 6.0) * (
            2.0 * x1 * y1 + x1 * y2 + x2 * y1 + 2.0 * x2 * y2
        )

    def reversed(self) -> Segment:
        """Return a copy of the segment with reversed endpoint orientation."""
        return Segment(p1=self.p2, p2=self.p1, t=self.t, id=self.id)
