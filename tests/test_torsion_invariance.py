"""Invariance and physical symmetry tests for ThinWallX v0.4 Torsion and Warping.

Verifies:
1. Root independence of normalized sectorial field and C_w.
2. Segment endpoint reversal invariance.
3. Segment definition reordering invariance.
4. Translation invariance (moderate and 1e12 large-coordinate offset).
5. Rigid proper rotation invariance (37 degrees).
6. Reflection / mirror symmetry behavior (omega reverses sign).
7. Zero-mean normalization integral int_A omega dA = 0.
"""

from __future__ import annotations

import math
import numpy as np
import pytest

from thinwallx.primitives import Node, Segment
from thinwallx.section import Section
from thinwallx.torsion import compute_torsion_warping


def make_channel_section() -> Section:
    """Standard C-channel section for invariance tests."""
    return Section.from_tuples(
        [
            ((0.0, -50.0), (0.0, 50.0), 3.0),
            ((0.0, 50.0), (40.0, 50.0), 3.0),
            ((0.0, -50.0), (40.0, -50.0), 3.0),
        ]
    )


class TestRootIndependence:
    """Invariant 2: Changing tree root must not affect normalized omega or C_w."""

    def test_channel_root_independence(self) -> None:
        """Verify different root nodes produce identical normalized omega and C_w."""
        sec = make_channel_section()
        v_count = len(sec.nodes)

        results = [
            compute_torsion_warping(sec, root_node_idx=r)
            for r in range(v_count)
        ]

        base = results[0]
        for idx, res in enumerate(results[1:], start=1):
            # J and Cw must be identical
            assert res.J == pytest.approx(base.J, abs=1e-12)
            assert res.Cw == pytest.approx(base.Cw, abs=1e-12)

            # Normalized node_omega must be identical
            for w_base, w_res in zip(base.node_omega, res.node_omega):
                assert w_res == pytest.approx(w_base, abs=1e-12)

            # Raw values must differ only by a uniform constant shift
            diffs = [
                w_raw_res - w_raw_base
                for w_raw_res, w_raw_base in zip(res.raw_node_omega, base.raw_node_omega)
            ]
            first_diff = diffs[0]
            for d in diffs[1:]:
                assert d == pytest.approx(first_diff, abs=1e-12)

    def test_branched_t_section_root_independence(self) -> None:
        """Verify root independence on a branched T-section."""
        sec = Section([
            Segment(p1=Node(0.0, 0.0), p2=Node(0.0, 100.0), t=4.0, id=0),
            Segment(p1=Node(-40.0, 100.0), p2=Node(0.0, 100.0), t=6.0, id=1),
            Segment(p1=Node(0.0, 100.0), p2=Node(40.0, 100.0), t=6.0, id=2),
        ])
        v_count = len(sec.nodes)
        results = [
            compute_torsion_warping(sec, root_node_idx=r)
            for r in range(v_count)
        ]
        base = results[0]
        for res in results[1:]:
            assert res.J == pytest.approx(base.J, abs=1e-12)
            assert res.Cw == pytest.approx(base.Cw, abs=1e-12)
            for w_base, w_res in zip(base.node_omega, res.node_omega):
                assert w_res == pytest.approx(w_base, abs=1e-12)

    def test_asymmetric_branched_section_root_independence(self) -> None:
        """Verify root independence on an asymmetric branched tree section."""
        sec = Section([
            Segment(p1=Node(0.0, 0.0), p2=Node(0.0, 80.0), t=4.0, id=0),
            Segment(p1=Node(-50.0, 80.0), p2=Node(0.0, 80.0), t=3.0, id=1),
            Segment(p1=Node(0.0, 80.0), p2=Node(30.0, 80.0), t=5.0, id=2),
            Segment(p1=Node(0.0, 0.0), p2=Node(40.0, 0.0), t=3.5, id=3),
        ])
        v_count = len(sec.nodes)
        results = [
            compute_torsion_warping(sec, root_node_idx=r)
            for r in range(v_count)
        ]
        base = results[0]
        for res in results[1:]:
            assert res.J == pytest.approx(base.J, rel=1e-12)
            assert res.Cw == pytest.approx(base.Cw, rel=1e-12)
            for w_base, w_res in zip(base.node_omega, res.node_omega):
                assert w_res == pytest.approx(w_base, rel=1e-12, abs=1e-12)

            diffs = [
                w_raw_res - w_raw_base
                for w_raw_res, w_raw_base in zip(res.raw_node_omega, base.raw_node_omega)
            ]
            first_diff = diffs[0]
            for d in diffs[1:]:
                assert d == pytest.approx(first_diff, rel=1e-12, abs=1e-12)


class TestReversalAndReorderingInvariance:
    """Invariants 3 & 4: Reversal of endpoints and reordering of segments."""

    def test_endpoint_reversal_invariance(self) -> None:
        """Reversing segment endpoints does not change J, C_w, or physical omega."""
        base = make_channel_section()
        rev_segments = [
            base.segments[0].reversed(),
            base.segments[1],
            base.segments[2].reversed(),
        ]
        rev_sec = Section(rev_segments)

        res_base = base.torsion_properties()
        res_rev = rev_sec.torsion_properties()

        assert res_rev.J == pytest.approx(res_base.J, abs=1e-12)
        assert res_rev.Cw == pytest.approx(res_base.Cw, rel=1e-12)

        # Compare normalized sectorial coordinates at shared nodes by position
        for n_base, w_base in zip(base.nodes, res_base.node_omega):
            for n_rev, w_rev in zip(rev_sec.nodes, res_rev.node_omega):
                if n_base.distance_to(n_rev) < 1e-8:
                    assert w_rev == pytest.approx(w_base, abs=1e-12)

    def test_segment_reordering_invariance(self) -> None:
        """Permuting the segment definitions must not change J, C_w, or omega."""
        base = make_channel_section()
        perm_segments = [base.segments[2], base.segments[0], base.segments[1]]
        perm_sec = Section(perm_segments)

        res_base = base.torsion_properties()
        res_perm = perm_sec.torsion_properties()

        assert res_perm.J == pytest.approx(res_base.J, abs=1e-12)
        assert res_perm.Cw == pytest.approx(res_base.Cw, rel=1e-12)

        for n_base, w_base in zip(base.nodes, res_base.node_omega):
            for n_perm, w_perm in zip(perm_sec.nodes, res_perm.node_omega):
                if n_base.distance_to(n_perm) < 1e-8:
                    assert w_perm == pytest.approx(w_base, abs=1e-12)


class TestTranslationInvariance:
    """Invariant 5: Rigid translation invariance, including 1e12 large translation."""

    def test_moderate_translation_invariance(self) -> None:
        """Verify J, C_w, and omega are invariant under rigid shift."""
        base = make_channel_section()
        trans = base.translated(150.0, -350.0)

        res_base = base.torsion_properties()
        res_trans = trans.torsion_properties()

        assert res_trans.J == pytest.approx(res_base.J, abs=1e-12)
        assert res_trans.Cw == pytest.approx(res_base.Cw, abs=1e-12)
        for w_base, w_trans in zip(res_base.node_omega, res_trans.node_omega):
            assert w_trans == pytest.approx(w_base, abs=1e-12)

    def test_large_translation_stability(self) -> None:
        """Verify 1e12 coordinate offset causes no precision loss or cancellation."""
        base = make_channel_section()
        offset = 1e12
        trans = base.translated(offset, -offset)

        res_base = base.torsion_properties()
        res_trans = trans.torsion_properties()

        assert res_trans.J == pytest.approx(res_base.J, rel=1e-11)
        assert res_trans.Cw == pytest.approx(res_base.Cw, rel=1e-11)
        for w_base, w_trans in zip(res_base.node_omega, res_trans.node_omega):
            assert w_trans == pytest.approx(w_base, rel=1e-10, abs=1e-10)


class TestRotationAndReflectionInvariance:
    """Invariants 6 & 7: Rigid proper rotation and reflection behavior."""

    def test_rotation_invariance(self) -> None:
        """Proper rotation preserves J, C_w, and omega at corresponding nodes."""
        base = make_channel_section()
        theta = math.radians(37.0)
        c, s = math.cos(theta), math.sin(theta)
        r_mat = np.array([[c, -s], [s, c]])

        rot = base.rotated(theta, origin=(0.0, 0.0))

        res_base = base.torsion_properties()
        res_rot = rot.torsion_properties()

        assert res_rot.J == pytest.approx(res_base.J, rel=1e-11)
        assert res_rot.Cw == pytest.approx(res_base.Cw, rel=1e-11)

        # Match rotated nodes: r' = R * r
        for n_base, w_base in zip(base.nodes, res_base.node_omega):
            expected_pos = r_mat @ np.array([n_base.x, n_base.y])
            for n_rot, w_rot in zip(rot.nodes, res_rot.node_omega):
                dist = math.hypot(n_rot.x - expected_pos[0], n_rot.y - expected_pos[1])
                if dist < 1e-8:
                    assert w_rot == pytest.approx(w_base, rel=1e-10, abs=1e-10)

    def test_reflection_reverses_sectorial_sign(self) -> None:
        """Reflection reverses the sign of omega while preserving J and C_w."""
        base = make_channel_section()
        # Reflect across x-axis: y -> -y
        refl_segments = [
            Segment(
                p1=Node(seg.p1.x, -seg.p1.y),
                p2=Node(seg.p2.x, -seg.p2.y),
                t=seg.t,
                id=seg.id,
            )
            for seg in base.segments
        ]
        refl_sec = Section(refl_segments)

        res_base = base.torsion_properties()
        res_refl = refl_sec.torsion_properties()

        assert res_refl.J == pytest.approx(res_base.J, rel=1e-11)
        assert res_refl.Cw == pytest.approx(res_base.Cw, rel=1e-11)

        for n_base, w_base in zip(base.nodes, res_base.node_omega):
            for n_refl, w_refl in zip(refl_sec.nodes, res_refl.node_omega):
                if abs(n_refl.x - n_base.x) < 1e-8 and abs(n_refl.y - (-n_base.y)) < 1e-8:
                    assert w_refl == pytest.approx(-w_base, rel=1e-10, abs=1e-10)


class TestZeroMeanNormalization:
    """Criterion 8: Normalized sectorial coordinate satisfies int_A omega dA = 0."""

    def test_zero_mean_property_across_sections(self) -> None:
        """Verify int_A omega dA is zero within floating-point tolerance."""
        sec = make_channel_section()
        res = sec.torsion_properties()
        assert abs(res.integral_omega_da) < 1e-11
