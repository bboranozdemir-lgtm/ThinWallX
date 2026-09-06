"""Section properties calculation for ThinWallX v0.1.

Implements exact closed-form straight-segment formulas for area, centroid,
second moments of area, and principal properties strictly according to ACTIVE_PHASE.md.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import sys
from typing import Iterable, Sequence

import numpy as np

from thinwallx.primitives import Segment
from thinwallx.exceptions import GeometryError, TopologyError


@dataclass(frozen=True)
class SectionProperties:
    """Computed section properties for an open thin-walled section.

    Attributes:
        area: Total cross-sectional area A = \\sum t_i * L_i.
        cx: Centroid horizontal coordinate \\bar{x}.
        cy: Centroid vertical coordinate \\bar{y}.
        ix_raw: Raw second moment of area about origin x-axis: I_{x,0} = \\int y^2 dA.
        iy_raw: Raw second moment of area about origin y-axis: I_{y,0} = \\int x^2 dA.
        ixy_raw: Raw product moment of area about origin: I_{xy,0} = \\int xy dA.
        ix: Centroidal second moment of area: I_x = I_{x,0} - A * \\bar{y}^2.
        iy: Centroidal second moment of area: I_y = I_{y,0} - A * \\bar{x}^2.
        ixy: Centroidal product moment of area: I_{xy} = I_{xy,0} - A * \\bar{x} * \\bar{y}.
        i1: Major principal second moment (I_1 >= I_2).
        i2: Minor principal second moment.
        theta_p: Principal-axis angle in radians for I_1 relative to x-axis, in (-pi/2, pi/2].
        is_degenerate: True if section inertia is isotropic (I_x == I_y and I_xy == 0).
    """

    area: float
    cx: float
    cy: float
    ix_raw: float
    iy_raw: float
    ixy_raw: float
    ix: float
    iy: float
    ixy: float
    i1: float
    i2: float
    theta_p: float
    is_degenerate: bool

    @property
    def centroid(self) -> tuple[float, float]:
        """Return (cx, cy) centroid tuple."""
        return (self.cx, self.cy)

    @property
    def theta_p_deg(self) -> float:
        """Principal-axis angle in degrees."""
        return math.degrees(self.theta_p)

    @property
    def inertia_matrix(self) -> np.ndarray:
        """2D inertia tensor matrix: [[I_x, -I_{xy}], [-I_{xy}, I_y]]."""
        return np.array([[self.ix, -self.ixy], [-self.ixy, self.iy]], dtype=float)


def _finite_fsum(values: Iterable[float], quantity: str) -> float:
    """Accurately sum floats and fail explicitly if the result is unrepresentable."""
    try:
        result = math.fsum(values)
    except (OverflowError, ValueError) as exc:
        raise GeometryError(
            f"{quantity} exceeds the supported floating-point range."
        ) from exc
    if not math.isfinite(result):
        raise GeometryError(f"{quantity} exceeds the supported floating-point range.")
    return result


def compute_properties(
    segments: Sequence[Segment], degeneracy_tol: float = 1e-12
) -> SectionProperties:
    """Compute exact centerline section properties from a sequence of segments.

    Engineering Formulas & Derivation References:
    1. Area:
       A = \\sum_i t_i * L_i
       (ACTIVE_PHASE.md, Area)

    2. Centroid:
       \\bar{x} = (1 / A) * \\sum_i t_i * \\int_{L_i} x ds = (1 / A) * \\sum_i t_i * L_i * (x_1 + x_2) / 2
       \\bar{y} = (1 / A) * \\sum_i t_i * \\int_{L_i} y ds = (1 / A) * \\sum_i t_i * L_i * (y_1 + y_2) / 2
       (ACTIVE_PHASE.md, Centroid)

    3. Raw Second Moments (about origin):
       I_{x,0} = \\sum_i t_i * \\int_{L_i} y^2 ds = \\sum_i t_i * (L_i / 3) * (y_1^2 + y_1*y_2 + y_2^2)
       I_{y,0} = \\sum_i t_i * \\int_{L_i} x^2 ds = \\sum_i t_i * (L_i / 3) * (x_1^2 + x_1*x_2 + x_2^2)
       I_{xy,0} = \\sum_i t_i * \\int_{L_i} xy ds = \\sum_i t_i * (L_i / 6) * (2*x_1*y_1 + x_1*y_2 + x_2*y_1 + 2*x_2*y_2)
       (ACTIVE_PHASE.md, Raw Second Moments)

    4. Centroidal Second Moments (Parallel Axis Theorem):
       I_x = I_{x,0} - A * \\bar{y}^2
       I_y = I_{y,0} - A * \\bar{x}^2
       I_{xy} = I_{xy,0} - A * \\bar{x} * \\bar{y}
       (ACTIVE_PHASE.md, Centroidal Second Moments)

    5. Principal Second Moments & Tensor:
       I = [[I_x, -I_{xy}], [-I_{xy}, I_y]]
       Eigenvalues ordered I_1 >= I_2:
       I_mean = (I_x + I_y) / 2
       R = sqrt(((I_x - I_y) / 2)^2 + I_{xy}^2)
       I_1 = I_mean + R
       I_2 = I_mean - R
       (ACTIVE_PHASE.md, Principal Properties)

    6. Principal-Axis Angle:
       theta_p = 0.5 * atan2(-2 * I_{xy}, I_x - I_y)
       Quadrant-safe formulation. If |I_x - I_y| <= tol and |I_{xy}| <= tol,
       section is degenerate (isotropic), theta_p = 0.0 by convention.
       (ACTIVE_PHASE.md, Principal Properties)

    Args:
        segments: Sequence of valid straight Segment objects.
        degeneracy_tol: Relative numerical tolerance used to determine isotropic
            degeneracy from the inertia-tensor scale.

    Returns:
        SectionProperties dataclass with all computed quantities.
    """
    if len(segments) == 0:
        raise TopologyError("Section properties require at least one segment.")
    if not math.isfinite(degeneracy_tol) or degeneracy_tol < 0.0:
        raise GeometryError(
            "The degeneracy tolerance must be finite and non-negative. "
            f"Got {degeneracy_tol}."
        )

    areas = [seg.area for seg in segments]
    if any(not math.isfinite(area) or area <= 0.0 for area in areas):
        raise GeometryError(
            "Every segment area t*L must be finite and strictly positive in floating-point."
        )
    total_area = _finite_fsum(areas, "Section area")

    # Use a deterministic endpoint as a local origin.  Relative coordinates
    # avoid forming first moments of order A*translation and make results
    # invariant to both segment ordering and large rigid translations.
    reference_x, reference_y = min(
        (node.x, node.y)
        for seg in segments
        for node in (seg.p1, seg.p2)
    )
    midpoint_offsets_x = [
        (seg.p1.x - reference_x) / 2.0 + (seg.p2.x - reference_x) / 2.0
        for seg in segments
    ]
    midpoint_offsets_y = [
        (seg.p1.y - reference_y) / 2.0 + (seg.p2.y - reference_y) / 2.0
        for seg in segments
    ]
    cx_offset = _finite_fsum(
        (
            (area / total_area) * midpoint_x
            for area, midpoint_x in zip(areas, midpoint_offsets_x, strict=True)
        ),
        "Centroid x-offset",
    )
    cy_offset = _finite_fsum(
        (
            (area / total_area) * midpoint_y
            for area, midpoint_y in zip(areas, midpoint_offsets_y, strict=True)
        ),
        "Centroid y-offset",
    )
    cx = reference_x + cx_offset
    cy = reference_y + cy_offset
    if not math.isfinite(cx) or not math.isfinite(cy):
        raise GeometryError("Section centroid exceeds the supported floating-point range.")

    # Preserve the specified raw origin moments for the public API.  fsum makes
    # their accumulation insensitive to segment order at normal float scale.
    ix_raw = _finite_fsum(
        (seg.t * seg.int_y2() for seg in segments), "Raw I_x"
    )
    iy_raw = _finite_fsum(
        (seg.t * seg.int_x2() for seg in segments), "Raw I_y"
    )
    ixy_raw = _finite_fsum(
        (seg.t * seg.int_xy() for seg in segments), "Raw I_xy"
    )

    # Evaluate the parallel-axis shift in its algebraically equivalent exact
    # centerline form, using coordinates relative to the centroid:
    #   I_x = sum(t * integral((y-cy)^2) ds), and analogously for I_y/I_xy.
    # Directly subtracting I_0 - A*c^2 loses all meaningful digits after a
    # large translation.  For a linear segment, midpoint/delta integration is
    # exactly L*(m^2 + delta^2/12), and the product integral is exactly
    # L*(mx*my + dx*dy/12).
    ix_terms: list[float] = []
    iy_terms: list[float] = []
    ixy_terms: list[float] = []
    for seg, area in zip(segments, areas, strict=True):
        x1 = (seg.p1.x - reference_x) - cx_offset
        x2 = (seg.p2.x - reference_x) - cx_offset
        y1 = (seg.p1.y - reference_y) - cy_offset
        y2 = (seg.p2.y - reference_y) - cy_offset
        midpoint_x = x1 / 2.0 + x2 / 2.0
        midpoint_y = y1 / 2.0 + y2 / 2.0
        dx = x2 - x1
        dy = y2 - y1

        ix_terms.append(area * (midpoint_y * midpoint_y + dy * dy / 12.0))
        iy_terms.append(area * (midpoint_x * midpoint_x + dx * dx / 12.0))
        ixy_terms.append(area * (midpoint_x * midpoint_y + dx * dy / 12.0))

    ix = _finite_fsum(ix_terms, "Centroidal I_x")
    iy = _finite_fsum(iy_terms, "Centroidal I_y")
    ixy = _finite_fsum(ixy_terms, "Centroidal I_xy")
    # Symmetric geometries can leave a cancellation residual of a few ulps in
    # Ixy.  Use a purely relative roundoff bound (no dimensional absolute floor).
    moment_scale = max(abs(ix), abs(iy))
    if abs(ixy) <= 32.0 * sys.float_info.epsilon * moment_scale:
        ixy = 0.0

    # Principal moments calculation
    diff = ix - iy
    r = math.hypot(diff / 2.0, ixy)
    i_mean = ix / 2.0 + iy / 2.0
    i1 = i_mean + r
    i2 = i_mean - r
    if not math.isfinite(i1) or not math.isfinite(i2):
        raise GeometryError(
            "Principal moments exceed the supported floating-point range."
        )

    # A covariance/inertia tensor is positive semidefinite.  Clamp only a
    # negative eigenvalue attributable to roundoff; never erase a small but
    # physically positive minor moment.
    eigenvalue_scale = max(abs(i1), abs(ix), abs(iy), abs(ixy))
    negative_roundoff_tol = 32.0 * sys.float_info.epsilon * eigenvalue_scale
    if i2 < 0.0:
        if abs(i2) <= negative_roundoff_tol:
            i2 = 0.0
        else:
            raise GeometryError(
                "Computed inertia tensor is not positive semidefinite within floating-point "
                f"roundoff (I2={i2})."
            )

    # Degeneracy check
    scale = max(abs(ix), abs(iy), abs(ixy))
    is_degenerate = scale == 0.0 or (
        abs(ix - iy) <= degeneracy_tol * scale
        and abs(ixy) <= degeneracy_tol * scale
    )

    if is_degenerate:
        theta_p = 0.0
    else:
        # Quadrant-safe atan2 convention: 0.5 * atan2(-2*Ixy, Ix - Iy)
        theta_p = 0.5 * math.atan2(-2.0 * ixy, ix - iy)
        # Principal axes are unoriented lines modulo pi.  Normalize the one
        # possible boundary value to the documented interval (-pi/2, pi/2].
        if theta_p <= -math.pi / 2.0:
            theta_p += math.pi

    return SectionProperties(
        area=total_area,
        cx=cx,
        cy=cy,
        ix_raw=ix_raw,
        iy_raw=iy_raw,
        ixy_raw=ixy_raw,
        ix=ix,
        iy=iy,
        ixy=ixy,
        i1=i1,
        i2=i2,
        theta_p=theta_p,
        is_degenerate=is_degenerate,
    )
