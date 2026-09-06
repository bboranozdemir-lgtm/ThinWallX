"""Validation and error handling tests for Sectalix v0.3 Shear Center Analysis.

Verifies:
1. Singular/rank-deficient sections raise SingularSectionError.
2. Invalid safety factors raise GeometryError.
3. Non-finite values or results raise GeometryError.
4. Input type checks raise GeometryError.
5. Residual torque argument validation.
"""

from __future__ import annotations

import math
import pytest

from sectalix.exceptions import GeometryError, SingularSectionError
from sectalix.primitives import Node, Segment
from sectalix.section import Section
from sectalix.shear_center import (
    ShearCenterResult,
    compute_shear_center,
    compute_shear_flow_torque,
)
from sectalix.shear_flow import calculate_shear_flow


def make_single_strip() -> Section:
    """A single flat strip, which has zero moment of area about one in-plane axis."""
    return Section.from_tuples([((0.0, 0.0), (100.0, 0.0), 2.0)])


def make_valid_angle() -> Section:
    """A valid non-singular angle section."""
    return Section.from_tuples(
        [((0.0, 0.0), (50.0, 0.0), 2.0), ((0.0, 0.0), (0.0, 50.0), 2.0)]
    )


class TestShearCenterValidation:
    """Validation checks for shear center analysis."""

    def test_singular_section_raises_singular_section_error(self) -> None:
        """Criterion 21: Singular / rank-deficient section fails explicitly."""
        strip = make_single_strip()
        with pytest.raises(SingularSectionError, match="singular or effectively rank-deficient"):
            compute_shear_center(strip)

    def test_invalid_safety_factor_raises_geometry_error(self) -> None:
        """Negative, zero, or non-finite safety factor raises GeometryError."""
        sec = make_valid_angle()
        for bad_sf in [0.0, -1.0, float("nan"), float("inf"), -1e5]:
            with pytest.raises(GeometryError, match="safety_factor must be finite and strictly positive"):
                compute_shear_center(sec, safety_factor=bad_sf)

    def test_non_section_raises_geometry_error(self) -> None:
        """Passing an object other than Section raises GeometryError."""
        with pytest.raises(GeometryError, match="Expected a Section instance"):
            compute_shear_center("not_a_section")  # type: ignore[arg-type]

    def test_torque_on_non_shear_flow_result_raises_geometry_error(self) -> None:
        """Passing non-ShearFlowResult to compute_shear_flow_torque raises GeometryError."""
        with pytest.raises(GeometryError, match="Expected a ShearFlowResult instance"):
            compute_shear_flow_torque(123.45)  # type: ignore[arg-type]

    def test_result_finiteness_validation(self) -> None:
        """ShearCenterResult with non-finite values raises GeometryError."""
        sec = make_valid_angle()
        flow_x = calculate_shear_flow(sec, 1.0, 0.0)
        flow_y = calculate_shear_flow(sec, 0.0, 1.0)

        with pytest.raises(GeometryError, match="must be finite"):
            ShearCenterResult(
                section=sec,
                ex=float("nan"),
                ey=0.0,
                x=0.0,
                y=0.0,
                tz_x=0.0,
                tz_y=0.0,
                centroid=(0.0, 0.0),
                flow_x=flow_x,
                flow_y=flow_y,
            )

    def test_residual_torque_missing_vy(self) -> None:
        """Calling residual_torque with scalar vx but missing vy raises GeometryError."""
        sec = make_valid_angle()
        res = sec.compute_shear_center()
        with pytest.raises(GeometryError, match="vy must be provided"):
            res.residual_torque(100.0)
