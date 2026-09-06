"""Invariance and covariance tests for Sectalix v0.3 Shear Center Analysis.

Verifies:
1. Rigid translation covariance (both moderate and 1e12 large translation).
2. Rigid rotation covariance (nontrivial angle, e.g. 37 degrees).
3. Segment endpoint reversal invariance.
4. Segment reordering invariance.
5. Basis-load magnitude scaling invariance.
6. Combined load residual-torque cancellation for arbitrary (Vx, Vy).
"""

from __future__ import annotations

import math
import numpy as np
import pytest

from sectalix.primitives import Node, Segment
from sectalix.section import Section
from sectalix.shear_center import compute_shear_center, compute_shear_flow_torque
from sectalix.shear_flow import calculate_shear_flow
from sectalix.shear_load import ShearLoad


def make_channel_section() -> Section:
    """Standard C-channel section for invariance tests."""
    return Section.from_tuples(
        [
            ((0.0, -50.0), (0.0, 50.0), 3.0),
            ((0.0, 50.0), (40.0, 50.0), 3.0),
            ((0.0, -50.0), (40.0, -50.0), 3.0),
        ]
    )


def make_l_section() -> Section:
    """Standard unsymmetric L-angle section for invariance tests."""
    return Section.from_tuples(
        [
            ((0.0, 0.0), (60.0, 0.0), 2.5),
            ((0.0, 0.0), (0.0, 40.0), 3.5),
        ]
    )


class TestTranslationInvariance:
    """Test translation covariance under moderate and large coordinate offsets."""

    def test_moderate_translation_covariance(self) -> None:
        """Verify e_s is unchanged and S translates by a."""
        base = make_channel_section()
        dx, dy = 125.0, -350.0
        translated = base.translated(dx, dy)

        res_base = base.compute_shear_center()
        res_trans = translated.compute_shear_center()

        # Offsets e_x, e_y are strictly invariant
        assert res_trans.ex == pytest.approx(res_base.ex, abs=1e-11)
        assert res_trans.ey == pytest.approx(res_base.ey, abs=1e-11)

        # Absolute coordinates translate by (dx, dy)
        assert res_trans.x == pytest.approx(res_base.x + dx, abs=1e-11)
        assert res_trans.y == pytest.approx(res_base.y + dy, abs=1e-11)

    def test_large_translation_stability(self) -> None:
        """Verify 1e12 coordinate offset causes no precision loss or cancellation."""
        base = make_channel_section()
        large_offset = 1e12
        translated = base.translated(large_offset, -large_offset)

        res_base = base.compute_shear_center()
        res_trans = translated.compute_shear_center()

        # Offsets remain accurate to machine precision
        assert res_trans.ex == pytest.approx(res_base.ex, rel=1e-10, abs=1e-10)
        assert res_trans.ey == pytest.approx(res_base.ey, rel=1e-10, abs=1e-10)

        # Basis torques are preserved
        assert res_trans.tz_x == pytest.approx(res_base.tz_x, rel=1e-10, abs=1e-10)
        assert res_trans.tz_y == pytest.approx(res_base.tz_y, rel=1e-10, abs=1e-10)

        # Absolute coordinates translate consistently
        assert res_trans.x == pytest.approx(res_base.x + large_offset, rel=1e-11)
        assert res_trans.y == pytest.approx(res_base.y - large_offset, rel=1e-11)


class TestRotationCovariance:
    """Test rotation covariance under non-trivial angles."""

    @pytest.mark.parametrize("angle_deg", [30.0, 37.0, 45.0, 71.5, 120.0, -42.3])
    def test_rotation_covariance(self, angle_deg: float) -> None:
        """Verify e_s' = R * e_s for arbitrary non-trivial rotation angles."""
        base = make_channel_section()
        theta = math.radians(angle_deg)
        c, s = math.cos(theta), math.sin(theta)
        r_mat = np.array([[c, -s], [s, c]])

        rotated = base.rotated(theta, origin=(0.0, 0.0))

        res_base = base.compute_shear_center()
        res_rot = rotated.compute_shear_center()

        # Expected rotated offset vector: e_s' = R * e_s
        expected_offset = r_mat @ np.array([res_base.ex, res_base.ey])

        assert res_rot.ex == pytest.approx(expected_offset[0], rel=1e-10, abs=1e-10)
        assert res_rot.ey == pytest.approx(expected_offset[1], rel=1e-10, abs=1e-10)

        # Expected rotated absolute coordinates: S' = R * S
        expected_s = r_mat @ np.array([res_base.x, res_base.y])
        assert res_rot.x == pytest.approx(expected_s[0], rel=1e-10, abs=1e-10)
        assert res_rot.y == pytest.approx(expected_s[1], rel=1e-10, abs=1e-10)


class TestEndpointReversalAndReorderingInvariance:
    """Test invariance to segment endpoint reversal and list permutation."""

    def test_segment_endpoint_reversal_invariance(self) -> None:
        """Reversing segment start/end nodes must not change the physical shear center."""
        base = make_channel_section()
        # Reverse 2 of the 3 segments
        reversed_segments = [
            base.segments[0].reversed(),
            base.segments[1],
            base.segments[2].reversed(),
        ]
        rev_section = Section(reversed_segments)

        res_base = base.compute_shear_center()
        res_rev = rev_section.compute_shear_center()

        assert res_rev.ex == pytest.approx(res_base.ex, abs=1e-12)
        assert res_rev.ey == pytest.approx(res_base.ey, abs=1e-12)
        assert res_rev.x == pytest.approx(res_base.x, abs=1e-12)
        assert res_rev.y == pytest.approx(res_base.y, abs=1e-12)

    def test_segment_reordering_invariance(self) -> None:
        """Reordering segment sequence in Section definition must not change shear center."""
        base = make_channel_section()
        # Permute segments: [2, 0, 1]
        permuted_segments = [base.segments[2], base.segments[0], base.segments[1]]
        perm_section = Section(permuted_segments)

        res_base = base.compute_shear_center()
        res_perm = perm_section.compute_shear_center()

        assert res_perm.ex == pytest.approx(res_base.ex, abs=1e-12)
        assert res_perm.ey == pytest.approx(res_base.ey, abs=1e-12)
        assert res_perm.x == pytest.approx(res_base.x, abs=1e-12)
        assert res_perm.y == pytest.approx(res_base.y, abs=1e-12)


class TestLoadScalingAndResidualTorque:
    """Test load scaling linearity and general combined load residual torque."""

    @pytest.mark.parametrize("scale", [0.01, 2.5, 100.0, -5.0])
    def test_scaled_basis_loads_yield_identical_shear_center(self, scale: float) -> None:
        """Verify T_z(c * V) / c yields identical shear-center offset."""
        sec = make_channel_section()
        res = sec.compute_shear_center()

        # Scale Vx basis load
        flow_scaled_x = calculate_shear_flow(sec, vx=scale, vy=0.0)
        tz_scaled_x = compute_shear_flow_torque(flow_scaled_x)
        ey_from_scaled = - (tz_scaled_x / scale)

        # Scale Vy basis load
        flow_scaled_y = calculate_shear_flow(sec, vx=0.0, vy=scale)
        tz_scaled_y = compute_shear_flow_torque(flow_scaled_y)
        ex_from_scaled = tz_scaled_y / scale

        assert ex_from_scaled == pytest.approx(res.ex, rel=1e-10, abs=1e-10)
        assert ey_from_scaled == pytest.approx(res.ey, rel=1e-10, abs=1e-10)

    @pytest.mark.parametrize("vx,vy", [
        (100.0, 200.0),
        (-500.0, 1500.0),
        (2500.0, -3200.0),
        (-9999.0, -8888.0),
    ])
    def test_general_combined_load_residual_torque_is_zero(
        self, vx: float, vy: float
    ) -> None:
        """Criterion 13: For arbitrary (Vx, Vy), T_z - (ex * Vy - ey * Vx) == 0."""
        # Exercise both length-unit changes and load magnitudes. The tolerance
        # has torque units and no absolute floor. 32 epsilon allows for the
        # short segment integration/solve chain and final cancellation, while
        # requiring relative accuracy better than 8e-15 in every case.
        for length_scale in (1e-3, 1.0, 1e3):
            sec = Section.from_tuples([
                ((0.0, 0.0), (60.0 * length_scale, 0.0), 2.5 * length_scale),
                ((0.0, 0.0), (0.0, 40.0 * length_scale), 3.5 * length_scale),
            ])
            res = sec.compute_shear_center()
            for load_scale in (1e-6, 1.0, 1e6):
                fx, fy = vx * load_scale, vy * load_scale
                residual = res.residual_torque(vx=fx, vy=fy)
                torque_scale = abs(res.ex * fy) + abs(res.ey * fx)
                assert math.isfinite(residual) and torque_scale > 0.0
                assert abs(residual) <= 32 * np.finfo(float).eps * torque_scale, (
                    length_scale, load_scale, residual, torque_scale
                )
