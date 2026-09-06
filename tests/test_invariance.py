"""Tests for physical and mathematical invariance properties.

Verifies acceptance criteria 4, 5, 6, and 7:
- Criterion 4: Invariance to reversing segment directions.
- Criterion 5: Invariance to reordering segment definitions.
- Criterion 6: Translation transforms centroid but leaves centroidal Ix, Iy, Ixy unchanged.
- Criterion 7: Rotation transforms inertia tensor according to tensor transformation law
  and preserves principal invariants I1, I2.
"""

import itertools
import math
import numpy as np
import pytest

from thinwallx.primitives import Node, Segment
from thinwallx.section import Section


def make_channel_section() -> Section:
    """Create a standard asymmetric/symmetric channel section for testing."""
    # Channel: web height 100, top flange 50, bottom flange 50, thickness 4
    # Endpoints: (50, 50) -> (0, 50) -> (0, -50) -> (50, -50)
    s1 = Segment(Node(50.0, 50.0), Node(0.0, 50.0), t=4.0)
    s2 = Segment(Node(0.0, 50.0), Node(0.0, -50.0), t=4.0)
    s3 = Segment(Node(0.0, -50.0), Node(50.0, -50.0), t=4.0)
    return Section([s1, s2, s3])


def make_asymmetric_section() -> Section:
    """Create an asymmetric multi-segment open section."""
    s1 = Segment(Node(0.0, 120.0), Node(60.0, 120.0), t=3.0)
    s2 = Segment(Node(0.0, 0.0), Node(0.0, 120.0), t=4.0)
    s3 = Segment(Node(-40.0, 0.0), Node(0.0, 0.0), t=5.0)
    s4 = Segment(Node(-40.0, 0.0), Node(-40.0, 25.0), t=3.0)
    return Section([s1, s2, s3, s4])


class TestReversalInvariance:
    """Criterion 4: Invariance to reversing segment direction."""

    @pytest.mark.parametrize("sec_factory", [make_channel_section, make_asymmetric_section])
    def test_reversing_all_segments(self, sec_factory) -> None:
        base_sec = sec_factory()
        rev_segments = [s.reversed() for s in base_sec.segments]
        rev_sec = Section(rev_segments)

        assert math.isclose(base_sec.area, rev_sec.area, rel_tol=1e-12)
        assert math.isclose(base_sec.cx, rev_sec.cx, rel_tol=1e-12, abs_tol=1e-12)
        assert math.isclose(base_sec.cy, rev_sec.cy, rel_tol=1e-12, abs_tol=1e-12)
        assert math.isclose(base_sec.Ix, rev_sec.Ix, rel_tol=1e-12)
        assert math.isclose(base_sec.Iy, rev_sec.Iy, rel_tol=1e-12)
        assert math.isclose(base_sec.Ixy, rev_sec.Ixy, rel_tol=1e-12, abs_tol=1e-12)
        assert math.isclose(base_sec.I1, rev_sec.I1, rel_tol=1e-12)
        assert math.isclose(base_sec.I2, rev_sec.I2, rel_tol=1e-12)
        assert math.isclose(base_sec.theta_p, rev_sec.theta_p, rel_tol=1e-12, abs_tol=1e-12)

    def test_reversing_subset_of_segments(self) -> None:
        base_sec = make_channel_section()
        # Reverse only segment 0 and 2
        mixed_segments = [
            base_sec.segments[0].reversed(),
            base_sec.segments[1],
            base_sec.segments[2].reversed(),
        ]
        mixed_sec = Section(mixed_segments)

        assert math.isclose(base_sec.area, mixed_sec.area, rel_tol=1e-12)
        assert math.isclose(base_sec.cx, mixed_sec.cx, rel_tol=1e-12)
        assert math.isclose(base_sec.cy, mixed_sec.cy, rel_tol=1e-12, abs_tol=1e-12)
        assert math.isclose(base_sec.Ix, mixed_sec.Ix, rel_tol=1e-12)
        assert math.isclose(base_sec.Iy, mixed_sec.Iy, rel_tol=1e-12)
        assert math.isclose(base_sec.Ixy, mixed_sec.Ixy, rel_tol=1e-12, abs_tol=1e-12)
        assert math.isclose(base_sec.I1, mixed_sec.I1, rel_tol=1e-12)
        assert math.isclose(base_sec.I2, mixed_sec.I2, rel_tol=1e-12)


class TestReorderingInvariance:
    """Criterion 5: Invariance to segment input order."""

    def test_all_permutations_of_segments(self) -> None:
        base_sec = make_channel_section()
        segs = list(base_sec.segments)

        for perm in itertools.permutations(segs):
            perm_sec = Section(perm)
            assert math.isclose(base_sec.area, perm_sec.area, rel_tol=1e-12)
            assert math.isclose(base_sec.cx, perm_sec.cx, rel_tol=1e-12)
            assert math.isclose(base_sec.cy, perm_sec.cy, rel_tol=1e-12, abs_tol=1e-12)
            assert math.isclose(base_sec.Ix, perm_sec.Ix, rel_tol=1e-12)
            assert math.isclose(base_sec.Iy, perm_sec.Iy, rel_tol=1e-12)
            assert math.isclose(base_sec.Ixy, perm_sec.Ixy, rel_tol=1e-12, abs_tol=1e-12)
            assert math.isclose(base_sec.I1, perm_sec.I1, rel_tol=1e-12)
            assert math.isclose(base_sec.I2, perm_sec.I2, rel_tol=1e-12)
            assert math.isclose(base_sec.theta_p, perm_sec.theta_p, rel_tol=1e-12, abs_tol=1e-12)


class TestTranslationInvariance:
    """Criterion 6: Translation changes centroid coordinates but preserves centroidal moments."""

    @pytest.mark.parametrize("dx, dy", [(100.0, -50.0), (-345.67, 890.12), (0.0, 0.0)])
    def test_translation(self, dx: float, dy: float) -> None:
        base_sec = make_asymmetric_section()
        trans_sec = base_sec.translated(dx, dy)

        # Area invariant
        assert math.isclose(trans_sec.area, base_sec.area, rel_tol=1e-12)

        # Centroid shifts exactly by (dx, dy)
        assert math.isclose(trans_sec.cx, base_sec.cx + dx, rel_tol=1e-12)
        assert math.isclose(trans_sec.cy, base_sec.cy + dy, rel_tol=1e-12)

        # Centroidal second moments are identical
        assert math.isclose(trans_sec.Ix, base_sec.Ix, rel_tol=1e-10)
        assert math.isclose(trans_sec.Iy, base_sec.Iy, rel_tol=1e-10)
        assert math.isclose(trans_sec.Ixy, base_sec.Ixy, rel_tol=1e-10, abs_tol=1e-10)

        # Principal moments and angle are identical
        assert math.isclose(trans_sec.I1, base_sec.I1, rel_tol=1e-10)
        assert math.isclose(trans_sec.I2, base_sec.I2, rel_tol=1e-10)
        assert math.isclose(trans_sec.theta_p, base_sec.theta_p, rel_tol=1e-10, abs_tol=1e-10)


class TestRotationCovariance:
    """Criterion 7: Rotation transforms inertia tensor correctly and preserves eigenvalues."""

    @pytest.mark.parametrize(
        "angle_rad",
        [
            math.pi / 6.0,
            math.pi / 4.0,
            math.pi / 3.0,
            math.pi / 2.0,
            math.pi,
            -math.pi / 5.0,
            1.234567,
        ],
    )
    def test_rotation_transform(self, angle_rad: float) -> None:
        base_sec = make_asymmetric_section()
        rot_sec = base_sec.rotated(angle_rad)

        # 1. Area is invariant
        assert math.isclose(rot_sec.area, base_sec.area, rel_tol=1e-12)

        # 2. Centroid rotates according to standard 2D rotation matrix:
        # [cx', cy']^T = R * [cx, cy]^T
        c = math.cos(angle_rad)
        s = math.sin(angle_rad)
        R = np.array([[c, -s], [s, c]], dtype=float)

        expected_centroid = R @ np.array([base_sec.cx, base_sec.cy])
        assert math.isclose(rot_sec.cx, expected_centroid[0], rel_tol=1e-10, abs_tol=1e-10)
        assert math.isclose(rot_sec.cy, expected_centroid[1], rel_tol=1e-10, abs_tol=1e-10)

        # 3. Inertia tensor transforms as I_rot = R * I_base * R^T
        I_base = base_sec.inertia_matrix
        expected_I_rot = R @ I_base @ R.T
        actual_I_rot = rot_sec.inertia_matrix

        np.testing.assert_allclose(actual_I_rot, expected_I_rot, rtol=1e-9, atol=1e-9)

        # 4. Principal invariants (I1, I2) are identical
        assert math.isclose(rot_sec.I1, base_sec.I1, rel_tol=1e-9)
        assert math.isclose(rot_sec.I2, base_sec.I2, rel_tol=1e-9)

        # 5. Principal axis angle rotates by angle_rad modulo pi
        # Both theta_p describe lines, so they are defined mod pi in (-pi/2, pi/2]
        def normalize_angle(theta: float) -> float:
            # Wrap to (-pi/2, pi/2]
            return (theta + math.pi / 2.0) % math.pi - math.pi / 2.0

        expected_theta_p = normalize_angle(base_sec.theta_p + angle_rad)
        actual_theta_p = normalize_angle(rot_sec.theta_p)

        # Difference should be near 0 or pi
        diff = abs(expected_theta_p - actual_theta_p)
        assert math.isclose(min(diff, abs(diff - math.pi)), 0.0, abs_tol=1e-9)
