"""Tests for physical and mathematical invariance of shear-flow fields (ThinWallX v0.2).

Verifies acceptance criteria 7, 8, 9, 10, 11, 12, 13:
- Criterion 7: Free-edge zero flow (q = 0 at every geometric free edge).
- Criterion 8: Resultant recovery (sum of integrated shear-flow vectors recovers V).
- Criterion 9: Load linearity and superposition.
- Criterion 10: Segment endpoint reversal preserves the physical shear-flow vector field.
- Criterion 11: Segment reordering preserves physical results.
- Criterion 12: Translation invariance.
- Criterion 13: Rotation covariance.
"""

import itertools
import math
import numpy as np
import pytest

from thinwallx.primitives import Node, Segment
from thinwallx.section import Section
from thinwallx.shear_load import ShearLoad


def make_channel_section() -> Section:
    """Channel section: web 100, flanges 50, thickness 4.0."""
    s1 = Segment(Node(50.0, 50.0), Node(0.0, 50.0), t=4.0)
    s2 = Segment(Node(0.0, 50.0), Node(0.0, -50.0), t=4.0)
    s3 = Segment(Node(0.0, -50.0), Node(50.0, -50.0), t=4.0)
    return Section([s1, s2, s3])


def make_t_section() -> Section:
    """Branched T-section: flange width 120, tf=6; web height 100, tw=4."""
    s1 = Segment(Node(-60.0, 0.0), Node(0.0, 0.0), t=6.0)
    s2 = Segment(Node(60.0, 0.0), Node(0.0, 0.0), t=6.0)
    s3 = Segment(Node(0.0, 0.0), Node(0.0, -100.0), t=4.0)
    return Section([s1, s2, s3])


def make_asymmetric_section() -> Section:
    """Asymmetric 4-segment open tree section."""
    s1 = Segment(Node(0.0, 120.0), Node(60.0, 120.0), t=3.0)
    s2 = Segment(Node(0.0, 0.0), Node(0.0, 120.0), t=4.0)
    s3 = Segment(Node(-40.0, 0.0), Node(0.0, 0.0), t=5.0)
    s4 = Segment(Node(-40.0, 0.0), Node(-40.0, 25.0), t=3.0)
    return Section([s1, s2, s3, s4])


class TestFreeEdgeCondition:
    """Criterion 7: Every geometric free edge satisfies q = 0."""

    @pytest.mark.parametrize("sec_factory", [make_channel_section, make_t_section, make_asymmetric_section])
    @pytest.mark.parametrize("vx, vy", [(0.0, 500.0), (300.0, 0.0), (-250.0, 400.0)])
    def test_free_edge_zero_flow(self, sec_factory, vx: float, vy: float) -> None:
        sec = sec_factory()
        res = sec.calculate_shear_flow(vx=vx, vy=vy)

        # Count occurrences of each canonical node across all segment endpoints
        # A node with degree 1 is a geometric free edge
        all_nodes = [node for s in sec.segments for node in (s.p1, s.p2)]
        from thinwallx.validation import cluster_nodes
        _, mapping = cluster_nodes(all_nodes, tol=sec._node_tolerance)

        deg: dict[int, int] = {}
        for s_idx in range(len(sec.segments)):
            u = mapping[2 * s_idx]
            v = mapping[2 * s_idx + 1]
            deg[u] = deg.get(u, 0) + 1
            deg[v] = deg.get(v, 0) + 1

        for s_idx, sf in enumerate(res.segment_flows):
            u = mapping[2 * s_idx]
            v = mapping[2 * s_idx + 1]
            # If node1 is a leaf, q at xi=0 must be 0
            if deg[u] == 1:
                assert math.isclose(sf.q0, 0.0, abs_tol=1e-10), (
                    f"Segment {s_idx} start is a free edge but q0={sf.q0}"
                )
            # If node2 is a leaf, q at xi=1 must be 0
            if deg[v] == 1:
                assert math.isclose(sf.q_end, 0.0, abs_tol=1e-10), (
                    f"Segment {s_idx} end is a free edge but q_end={sf.q_end}"
                )


class TestResultantRecovery:
    """Criterion 8: Integrated physical shear-flow vectors recover the applied resultant."""

    @pytest.mark.parametrize("sec_factory", [make_channel_section, make_t_section, make_asymmetric_section])
    @pytest.mark.parametrize(
        "vx, vy",
        [(0.0, 1000.0), (1000.0, 0.0), (450.0, -780.0), (-1234.56, 789.01)],
    )
    def test_resultant_recovery(self, sec_factory, vx: float, vy: float) -> None:
        sec = sec_factory()
        res = sec.calculate_shear_flow(vx=vx, vy=vy)

        recovered = res.recovered_resultant
        expected = np.array([vx, vy], dtype=float)

        np.testing.assert_allclose(recovered, expected, rtol=1e-10, atol=1e-10)
        assert res.resultant_error < 1e-9 * max(np.linalg.norm(expected), 1.0)


class TestSuperposition:
    """Criterion 9: Load linearity and superposition."""

    @pytest.mark.parametrize("sec_factory", [make_channel_section, make_asymmetric_section])
    def test_linearity_and_superposition(self, sec_factory) -> None:
        sec = sec_factory()

        v1 = ShearLoad(200.0, -150.0)
        v2 = ShearLoad(-50.0, 300.0)
        a, b = 2.5, -1.7

        res1 = sec.calculate_shear_flow(v1)
        res2 = sec.calculate_shear_flow(v2)
        res_comb = sec.calculate_shear_flow(a * v1 + b * v2)

        xi_test = [0.0, 0.25, 0.5, 0.75, 1.0]
        for s_idx in range(len(sec.segments)):
            for xi in xi_test:
                q1 = res1.flow_at(s_idx, xi)
                q2 = res2.flow_at(s_idx, xi)
                q_comb = res_comb.flow_at(s_idx, xi)

                expected_q = a * q1 + b * q2
                assert math.isclose(q_comb, expected_q, rel_tol=1e-11, abs_tol=1e-11)


class TestReversalInvariance:
    """Criterion 10: Segment endpoint reversal preserves the physical shear-flow vector field."""

    @pytest.mark.parametrize("sec_factory", [make_channel_section, make_t_section, make_asymmetric_section])
    def test_reversing_segments_preserves_vector_field(self, sec_factory) -> None:
        base_sec = sec_factory()
        # Reverse all segments
        rev_segments = [s.reversed() for s in base_sec.segments]
        rev_sec = Section(rev_segments)

        load = ShearLoad(350.0, -600.0)
        res_base = base_sec.calculate_shear_flow(load)
        res_rev = rev_sec.calculate_shear_flow(load)

        # For reversed segment, normalized position xi corresponds to 1 - xi on base segment,
        # but the tangent is reversed, so q_rev(xi) = -q_base(1 - xi),
        # making the physical vector q_rev * t_rev = (-q_base) * (-t_base) = q_base * t_base!
        xi_samples = [0.0, 0.2, 0.5, 0.8, 1.0]
        for i in range(len(base_sec.segments)):
            for xi in xi_samples:
                vec_base = res_base.vector_at(i, xi)
                vec_rev = res_rev.vector_at(i, 1.0 - xi)

                np.testing.assert_allclose(vec_base, vec_rev, rtol=1e-10, atol=1e-10)

        # Resultant recovery is identical
        np.testing.assert_allclose(
            res_base.recovered_resultant, res_rev.recovered_resultant, rtol=1e-11, atol=1e-11
        )


class TestReorderingInvariance:
    """Criterion 11: Segment reordering preserves physical results."""

    def test_all_segment_permutations(self) -> None:
        base_sec = make_channel_section()
        segs = list(base_sec.segments)
        load = ShearLoad(400.0, -700.0)
        res_base = base_sec.calculate_shear_flow(load)

        for perm in itertools.permutations(segs):
            perm_sec = Section(perm)
            res_perm = perm_sec.calculate_shear_flow(load)

            # 1. Recovered resultant is invariant
            np.testing.assert_allclose(
                res_perm.recovered_resultant, res_base.recovered_resultant, rtol=1e-11, atol=1e-11
            )

            # 2. Physical vector flow field and scalar flow on each segment are invariant
            for perm_idx, perm_seg in enumerate(perm_sec.segments):
                # Find matching segment in base_sec
                base_idx = next(
                    b_i for b_i, b_s in enumerate(base_sec.segments)
                    if b_s.p1.is_close(perm_seg.p1) and b_s.p2.is_close(perm_seg.p2)
                )

                for xi in [0.0, 0.25, 0.5, 0.75, 1.0]:
                    vec_perm = res_perm.vector_at(perm_idx, xi)
                    vec_base = res_base.vector_at(base_idx, xi)
                    np.testing.assert_allclose(vec_perm, vec_base, rtol=1e-11, atol=1e-11)

                    q_perm = res_perm.flow_at(perm_idx, xi)
                    q_base = res_base.flow_at(base_idx, xi)
                    assert math.isclose(q_perm, q_base, rel_tol=1e-11, abs_tol=1e-11)


class TestTranslationInvariance:
    """Criterion 12: Translation invariance."""

    @pytest.mark.parametrize("dx, dy", [(150.0, -220.0), (-500.0, 1000.0)])
    def test_translation(self, dx: float, dy: float) -> None:
        base_sec = make_asymmetric_section()
        trans_sec = base_sec.translated(dx, dy)

        load = ShearLoad(-250.0, 800.0)
        res_base = base_sec.calculate_shear_flow(load)
        res_trans = trans_sec.calculate_shear_flow(load)

        # Scalar shear flow along each segment must be completely invariant
        xi_samples = [0.0, 0.25, 0.5, 0.75, 1.0]
        for s_idx in range(len(base_sec.segments)):
            for xi in xi_samples:
                q_base = res_base.flow_at(s_idx, xi)
                q_trans = res_trans.flow_at(s_idx, xi)
                assert math.isclose(q_base, q_trans, rel_tol=1e-10, abs_tol=1e-10)

    def test_large_translation_invariance(self) -> None:
        """Centroid-relative shear flow must be invariant under huge translations (1e12)."""
        base_sec = make_asymmetric_section()
        trans_sec = base_sec.translated(1.0e12, -1.0e12)

        load = ShearLoad(500.0, -800.0)
        res_base = base_sec.calculate_shear_flow(load)
        res_trans = trans_sec.calculate_shear_flow(load)

        # Verify scalar shear flow and physical vector flow match to within 1e-10 tolerance
        xi_samples = [0.0, 0.25, 0.5, 0.75, 1.0]
        for s_idx in range(len(base_sec.segments)):
            for xi in xi_samples:
                q_base = res_base.flow_at(s_idx, xi)
                q_trans = res_trans.flow_at(s_idx, xi)
                assert math.isclose(q_base, q_trans, rel_tol=1e-10, abs_tol=1e-10)

                vec_base = res_base.vector_at(s_idx, xi)
                vec_trans = res_trans.vector_at(s_idx, xi)
                np.testing.assert_allclose(vec_trans, vec_base, rtol=1e-10, atol=1e-10)

        # Recovered resultant must match load
        np.testing.assert_allclose(res_trans.recovered_resultant, load.vector, rtol=1e-10)



class TestRotationCovariance:
    """Criterion 13: Rotation covariance."""

    @pytest.mark.parametrize("angle_rad", [math.pi / 6.0, math.pi / 4.0, math.pi / 2.0, 1.2345])
    def test_rotation_covariance(self, angle_rad: float) -> None:
        base_sec = make_asymmetric_section()
        rot_sec = base_sec.rotated(angle_rad)

        c = math.cos(angle_rad)
        s = math.sin(angle_rad)
        R = np.array([[c, -s], [s, c]], dtype=float)

        base_load = ShearLoad(300.0, -500.0)
        # Rotated load vector V_rot = R * V_base
        v_rot = R @ base_load.vector
        rot_load = ShearLoad(float(v_rot[0]), float(v_rot[1]))

        res_base = base_sec.calculate_shear_flow(base_load)
        res_rot = rot_sec.calculate_shear_flow(rot_load)

        # At corresponding positions on each segment:
        # Physical shear-flow vector must rotate: q_rot(xi) * t_rot == R * (q_base(xi) * t_base)
        xi_samples = [0.0, 0.25, 0.5, 0.75, 1.0]
        for s_idx in range(len(base_sec.segments)):
            for xi in xi_samples:
                vec_base = res_base.vector_at(s_idx, xi)
                expected_vec_rot = R @ vec_base
                actual_vec_rot = res_rot.vector_at(s_idx, xi)

                np.testing.assert_allclose(
                    actual_vec_rot, expected_vec_rot, rtol=1e-9, atol=1e-9
                )
