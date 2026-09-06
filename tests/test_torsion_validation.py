"""Validation and error handling tests for ThinWallX v0.4 Torsion and Warping.

Verifies:
1. Singular/rank-deficient sections raise SingularSectionError.
2. Invalid root node indices raise GeometryError.
3. Non-Section input raises GeometryError.
4. Out-of-bounds or non-finite evaluation arguments raise GeometryError.
5. Dataclass validation rejects invalid values.
"""

from __future__ import annotations

import math
import pytest

from thinwallx.exceptions import GeometryError, SingularSectionError
from thinwallx.primitives import Node, Segment
from thinwallx.section import Section
from thinwallx.torsion import (
    SegmentWarping,
    TorsionWarpingResult,
    compute_torsion_warping,
)


def make_single_strip() -> Section:
    """A single flat strip with singular 2D inertia matrix."""
    return Section.from_tuples([((0.0, 0.0), (100.0, 0.0), 2.0)])


def make_valid_angle() -> Section:
    """A valid non-singular angle section."""
    return Section.from_tuples(
        [((0.0, 0.0), (50.0, 0.0), 2.0), ((0.0, 0.0), (0.0, 50.0), 2.0)]
    )


class TestTorsionValidation:
    """Validation tests for torsion and warping analysis."""

    def test_singular_section_raises_singular_section_error(self) -> None:
        """Singular section fails when computing shear center pole."""
        strip = make_single_strip()
        with pytest.raises(SingularSectionError):
            compute_torsion_warping(strip)

    def test_invalid_root_node_index_raises_geometry_error(self) -> None:
        """Root node index outside valid range raises GeometryError."""
        sec = make_valid_angle()
        with pytest.raises(GeometryError, match="root_node_idx must be in"):
            compute_torsion_warping(sec, root_node_idx=999)
        with pytest.raises(GeometryError, match="root_node_idx must be in"):
            compute_torsion_warping(sec, root_node_idx=-1)

    def test_non_section_input_raises_geometry_error(self) -> None:
        """Passing non-Section object raises GeometryError."""
        with pytest.raises(GeometryError, match="Expected a Section instance"):
            compute_torsion_warping("invalid")  # type: ignore[arg-type]

    def test_xi_out_of_bounds_raises_geometry_error(self) -> None:
        """Evaluating omega_at_xi outside [0, 1] raises GeometryError."""
        sec = make_valid_angle()
        res = sec.torsion_properties()
        sw = res[0]

        with pytest.raises(GeometryError, match="xi must be in"):
            sw.omega_at_xi(-0.5)
        with pytest.raises(GeometryError, match="xi must be in"):
            sw.omega_at_xi(1.5)
        with pytest.raises(GeometryError, match="xi must be finite"):
            sw.omega_at_xi(float("nan"))

    def test_s_non_finite_raises_geometry_error(self) -> None:
        """Evaluating omega_at with non-finite s raises GeometryError."""
        sec = make_valid_angle()
        res = sec.torsion_properties()
        sw = res[0]

        with pytest.raises(GeometryError, match="s must be finite"):
            sw.omega_at(float("inf"))

    def test_result_finiteness_and_positivity_validation(self) -> None:
        """TorsionWarpingResult validates J > 0 and Cw >= 0."""
        sec = make_valid_angle()

        # J <= 0 rejected
        with pytest.raises(GeometryError, match="J must be strictly positive"):
            TorsionWarpingResult(
                section=sec,
                J=-1.0,
                Cw=0.0,
                omega_mean=0.0,
                shear_center=(0.0, 0.0),
                node_omega=(0.0, 0.0, 0.0),
                segment_warpings=(),
                raw_node_omega=(0.0, 0.0, 0.0),
            )

        # Cw < 0 rejected
        with pytest.raises(GeometryError, match="cannot be negative"):
            TorsionWarpingResult(
                section=sec,
                J=10.0,
                Cw=-5.0,
                omega_mean=0.0,
                shear_center=(0.0, 0.0),
                node_omega=(0.0, 0.0, 0.0),
                segment_warpings=(),
                raw_node_omega=(0.0, 0.0, 0.0),
            )
