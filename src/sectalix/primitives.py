"""Geometric primitives for straight thin-wall centerline segments.

The line integrals in this module use the closed-form straight-segment identities
documented in ``docs/THEORY_AND_CONVENTIONS.md``.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

from sectalix.exceptions import GeometryError


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
        # Reject only an exactly collapsed float representation. A dimensional
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
        """Closed-form line integral: \\int_L x ds = L * (x1 + x2) / 2.

        Parameterizing the straight segment with s in [0, L], x(s) varies
        linearly between x1 and x2, giving the stated expression directly.
        """
        return self.length * (self.p1.x + self.p2.x) / 2.0

    def int_y(self) -> float:
        """Closed-form line integral: \\int_L y ds = L * (y1 + y2) / 2.

        Parameterizing the straight segment with s in [0, L], y(s) varies
        linearly between y1 and y2, giving the stated expression directly.
        """
        return self.length * (self.p1.y + self.p2.y) / 2.0

    def int_x2(self) -> float:
        """Closed-form line integral: \\int_L x^2 ds = (L / 3) * (x1^2 + x1*x2 + x2^2).

        This follows by squaring the linear endpoint interpolation for x(s)
        and integrating the resulting quadratic polynomial over [0, L].
        """
        x1, x2 = self.p1.x, self.p2.x
        return (self.length / 3.0) * (x1 * x1 + x1 * x2 + x2 * x2)

    def int_y2(self) -> float:
        """Closed-form line integral: \\int_L y^2 ds = (L / 3) * (y1^2 + y1*y2 + y2^2)."""
        y1, y2 = self.p1.y, self.p2.y
        return (self.length / 3.0) * (y1 * y1 + y1 * y2 + y2 * y2)

    def int_xy(self) -> float:
        """Closed-form line integral: \\int_L xy ds.

        For linearly interpolated x(s) and y(s),

        \\int_L xy ds = (L / 6) *
        (2*x1*y1 + x1*y2 + x2*y1 + 2*x2*y2).
        """
        x1, y1 = self.p1.x, self.p1.y
        x2, y2 = self.p2.x, self.p2.y
        return (self.length / 6.0) * (
            2.0 * x1 * y1 + x1 * y2 + x2 * y1 + 2.0 * x2 * y2
        )

    def reversed(self) -> Segment:
        """Return a copy of the segment with reversed endpoint orientation."""
        return Segment(p1=self.p2, p2=self.p1, t=self.t, id=self.id)
