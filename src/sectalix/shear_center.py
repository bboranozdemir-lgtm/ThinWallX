"""Thin-wall shear-center analysis for open sections.

The shear center is obtained from centroidal torque integration of the solved
transverse shear-flow fields. Governing conventions are documented in
``docs/THEORY_AND_CONVENTIONS.md``.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import TYPE_CHECKING

import numpy as np

from sectalix.exceptions import GeometryError
from sectalix.shear_flow import ShearFlowResult, calculate_shear_flow
from sectalix.shear_load import ShearLoad

if TYPE_CHECKING:
    from sectalix.section import Section


def compute_shear_flow_torque(result: ShearFlowResult) -> float:
    """Compute the torque of the physical shear-flow field about the centroid.

    For each straight segment, the centroid-relative cross-product factor is
    constant along the segment, so the torque contribution reduces to that
    factor multiplied by the analytical integral of the scalar shear flow.

    A local reference shift is used for centroid-relative distances to reduce
    cancellation when the entire geometry is translated far from the origin.

    Args:
        result: A solved ShearFlowResult from calculate_shear_flow.

    Returns:
        Total scalar torque about the section centroid, positive about +z.

    Raises:
        GeometryError: If any computed torque contribution or total is non-finite.
    """
    if hasattr(result, "torque") and not isinstance(result, ShearFlowResult):
        return float(result.torque)

    if not isinstance(result, ShearFlowResult):
        raise GeometryError(
            f"Expected a ShearFlowResult instance, got {type(result).__name__}."
        )

    section = result.section
    segments = section.segments
    all_nodes = [node for seg in segments for node in (seg.p1, seg.p2)]
    ref_x, ref_y = min((n.x, n.y) for n in all_nodes)
    total_area = section.area

    # Centroid relative to local reference
    cx_loc = sum(
        seg.area * ((seg.p1.x - ref_x) + (seg.p2.x - ref_x)) / 2.0
        for seg in segments
    ) / total_area
    cy_loc = sum(
        seg.area * ((seg.p1.y - ref_y) + (seg.p2.y - ref_y)) / 2.0
        for seg in segments
    ) / total_area

    total_torque = 0.0
    for sf in result.segment_flows:
        seg = sf.segment
        xc1 = (seg.p1.x - ref_x) - cx_loc
        yc1 = (seg.p1.y - ref_y) - cy_loc
        tx = sf.tangent[0]
        ty = sf.tangent[1]

        # [r_c(0) x t]_z = xc1 * ty - yc1 * tx
        cross_z = xc1 * ty - yc1 * tx
        q_integral = sf.scalar_integral()
        seg_torque = cross_z * q_integral

        if not math.isfinite(seg_torque):
            raise GeometryError(
                f"Non-finite segment torque contribution on segment {seg.id}: "
                f"cross_z={cross_z}, q_integral={q_integral}."
            )
        total_torque += seg_torque

    if not math.isfinite(total_torque):
        raise GeometryError(
            f"Non-finite total shear-flow torque computed: total_torque={total_torque}."
        )

    return total_torque


@dataclass(frozen=True)
class ShearCenterResult:
    """Container holding the results of a thin-walled shear-center analysis.

    Attributes:
        section: The analyzed Section.
        ex: Centroid-relative shear-center offset in x (e_x = x_s - cx).
        ey: Centroid-relative shear-center offset in y (e_y = y_s - cy).
        x: Absolute horizontal shear-center coordinate x_s.
        y: Absolute vertical shear-center coordinate y_s.
        tz_x: Centroidal torque for unit +x transverse shear load [1, 0]^T.
        tz_y: Centroidal torque for unit +y transverse shear load [0, 1]^T.
        centroid: Centroid coordinates (cx, cy).
        flow_x: Solved ShearFlowResult for unit +x basis load.
        flow_y: Solved ShearFlowResult for unit +y basis load.
    """

    section: Section
    ex: float
    ey: float
    x: float
    y: float
    tz_x: float
    tz_y: float
    centroid: tuple[float, float]
    flow_x: ShearFlowResult
    flow_y: ShearFlowResult

    def __post_init__(self) -> None:
        """Validate shear center results for finiteness."""
        for name, val in [
            ("ex", self.ex),
            ("ey", self.ey),
            ("x", self.x),
            ("y", self.y),
            ("tz_x", self.tz_x),
            ("tz_y", self.tz_y),
            ("centroid.cx", self.centroid[0]),
            ("centroid.cy", self.centroid[1]),
        ]:
            if not math.isfinite(val):
                raise GeometryError(f"ShearCenterResult property '{name}' must be finite. Got {val}.")

    @property
    def offset(self) -> tuple[float, float]:
        """Centroid-relative offset (ex, ey)."""
        return (self.ex, self.ey)

    @property
    def coordinates(self) -> tuple[float, float]:
        """Absolute coordinates (x_s, y_s)."""
        return (self.x, self.y)

    def residual_torque(
        self,
        vx: float | ShearLoad,
        vy: float | None = None,
        safety_factor: float = 1e4,
    ) -> float:
        """Evaluate residual torque for a transverse shear load V = [Vx, Vy]^T.

        The residual is

            T_z(V) - (e_x * V_y - e_y * V_x).

        It should be zero within the numerical tolerance of the solved model.

        Args:
            vx: Force in x direction, or a ShearLoad object.
            vy: Force in y direction (ignored if vx is ShearLoad).
            safety_factor: Multiplier for positive definiteness check.

        Returns:
            Residual torque scalar.
        """
        if isinstance(vx, ShearLoad):
            load = vx
        else:
            if vy is None:
                raise GeometryError("vy must be provided when vx is a scalar force.")
            load = ShearLoad(vx=float(vx), vy=float(vy))

        flow = calculate_shear_flow(self.section, load, safety_factor=safety_factor)
        tz_actual = compute_shear_flow_torque(flow)
        tz_predicted = self.ex * load.vy - self.ey * load.vx
        return tz_actual - tz_predicted


def compute_shear_center(
    section: Section,
    safety_factor: float = 1e4,
) -> ShearCenterResult:
    """Compute the shear-center location for an open thin-walled section.

    Unit transverse shear loads in the x and y directions are solved, their
    centroidal torques are integrated, and the offsets follow from

        T_z(V) = e_x * V_y - e_y * V_x.

    Args:
        section: Validated open Section object.
        safety_factor: Safety multiplier on machine epsilon for positive definiteness.

    Returns:
        ShearCenterResult object containing offsets, absolute coordinates, and basis torques.

    Raises:
        GeometryError: If section is invalid or values are non-finite.
        SingularSectionError: If the section inertia matrix is singular/rank-deficient.
    """
    from sectalix.section import Section

    if not isinstance(section, Section):
        raise GeometryError(f"Expected a Section instance, got {type(section).__name__}.")
    if hasattr(section, "is_closed") and section.is_closed:
        from sectalix.closed_torsion import compute_closed_shear_center
        return compute_closed_shear_center(section, safety_factor=safety_factor)

    section.validate()

    if not math.isfinite(safety_factor) or safety_factor <= 0.0:
        raise GeometryError(
            f"safety_factor must be finite and strictly positive. Got safety_factor={safety_factor}."
        )

    # Unit basis load cases:
    # 1. V^{(x)} = [1, 0]^T
    flow_x = calculate_shear_flow(
        section, vx=1.0, vy=0.0, safety_factor=safety_factor
    )
    tz_x = compute_shear_flow_torque(flow_x)

    # 2. V^{(y)} = [0, 1]^T
    flow_y = calculate_shear_flow(
        section, vx=0.0, vy=1.0, safety_factor=safety_factor
    )
    tz_y = compute_shear_flow_torque(flow_y)

    # From T_z(V) = e_x * V_y - e_y * V_x:
    # V = [1, 0]^T ==> T_z^{(x)} = -e_y ==> e_y = -T_z^{(x)}
    # V = [0, 1]^T ==> T_z^{(y)} =  e_x ==> e_x =  T_z^{(y)}
    ex = tz_y
    ey = -tz_x

    cx, cy = section.centroid
    xs = cx + ex
    ys = cy + ey

    if not (
        math.isfinite(ex)
        and math.isfinite(ey)
        and math.isfinite(xs)
        and math.isfinite(ys)
    ):
        raise GeometryError(
            f"Non-finite shear center calculated: ex={ex}, ey={ey}, xs={xs}, ys={ys}."
        )

    return ShearCenterResult(
        section=section,
        ex=ex,
        ey=ey,
        x=xs,
        y=ys,
        tz_x=tz_x,
        tz_y=tz_y,
        centroid=(cx, cy),
        flow_x=flow_x,
        flow_y=flow_y,
    )
