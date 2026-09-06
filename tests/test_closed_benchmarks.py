"""Comprehensive analytical benchmarks and invariance tests for closed sections (Sectalix v0.5).

Covers:
- Benchmark 1: Single-cell rectangular box Bredt-Batho J = 2*a^2*b^2*t / (a + b).
- Benchmark 2: Doubly symmetric rectangular box shear center S = C.
- Benchmark 3: Rectangular closed box warping constant C_w = t*a^2*b^2*(a - b)^2 / [24*(a + b)].
               Square box special case C_w = 0.
- Benchmark 4: Symmetric two-cell box: analytical H, phi_1 = phi_2, F_shared = 0, J = 8*b^2*h^2*t / (2*b + h).
- Benchmark 5: Unsymmetric closed cell: nontrivial shear center, residual torque consistency.
- Benchmark 6: Unsymmetric multi-cell section: unequal cell flows, compatibility, resultant recovery.
- Benchmark 7: Cut independence on single-cell and multi-cell sections.
- Benchmark 8: Large coordinate translation invariance (1e12 offset).
- Benchmark 9: In-plane rotation covariance (nontrivial angle 37 deg).
- Benchmark 10: Uniform scaling laws: e_s ~ L, J ~ L^4, omega ~ L^2, C_w ~ L^6.
"""

import math
import sys
import numpy as np
import pytest

from sectalix.closed_section import ClosedSection
from sectalix.closed_shear_flow import calculate_closed_shear_flow
from sectalix.primitives import Node, Segment
from sectalix.shear_load import ShearLoad


# ==============================================================================
# Independent Analytical Oracles
# ==============================================================================

def oracle_single_cell_box_J(a: float, b: float, t: float) -> float:
    """Analytical Bredt-Batho J for rectangular box a x b with uniform thickness t."""
    return (2.0 * (a**2) * (b**2) * t) / (a + b)


def oracle_single_cell_box_Cw(a: float, b: float, t: float) -> float:
    """Analytical warping constant C_w for rectangular box a x b with uniform thickness t."""
    return (t * (a**2) * (b**2) * ((a - b) ** 2)) / (24.0 * (a + b))


def oracle_symmetric_twocell_box_J(b: float, h: float, t: float) -> float:
    """Analytical Bredt-Batho J for two identical cells side-by-side."""
    return (8.0 * (b**2) * (h**2) * t) / (2.0 * b + h)


# ==============================================================================
# Geometry Helpers
# ==============================================================================

def make_rect_box(a: float, b: float, t: float) -> ClosedSection:
    """Build rectangular ClosedSection of dimensions a x b and thickness t."""
    n1 = Node(0.0, 0.0)
    n2 = Node(a, 0.0)
    n3 = Node(a, b)
    n4 = Node(0.0, b)
    return ClosedSection([
        Segment(n1, n2, t=t),
        Segment(n2, n3, t=t),
        Segment(n3, n4, t=t),
        Segment(n4, n1, t=t),
    ])


def make_twocell_box(b: float, h: float, t: float) -> ClosedSection:
    """Build two-cell ClosedSection with shared vertical web at x = b."""
    n1 = Node(0.0, 0.0)
    n2 = Node(b, 0.0)
    n3 = Node(2.0 * b, 0.0)
    n4 = Node(2.0 * b, h)
    n5 = Node(b, h)
    n6 = Node(0.0, h)
    return ClosedSection([
        Segment(n1, n2, t=t),
        Segment(n2, n3, t=t),
        Segment(n3, n4, t=t),
        Segment(n4, n5, t=t),
        Segment(n5, n6, t=t),
        Segment(n6, n1, t=t),
        Segment(n2, n5, t=t),  # shared web (index 6)
    ])


def make_unsymmetric_box(
    a: float = 100.0, b: float = 60.0, t_bot: float = 4.0, t_top: float = 2.0, t_left: float = 3.0, t_right: float = 5.0
) -> ClosedSection:
    """Build unsymmetric rectangular ClosedSection with different wall thicknesses."""
    n1 = Node(0.0, 0.0)
    n2 = Node(a, 0.0)
    n3 = Node(a, b)
    n4 = Node(0.0, b)
    return ClosedSection([
        Segment(n1, n2, t=t_bot),
        Segment(n2, n3, t=t_right),
        Segment(n3, n4, t=t_top),
        Segment(n4, n1, t=t_left),
    ])


def make_unsymmetric_twocell(
    b1: float = 80.0, b2: float = 40.0, h: float = 50.0, t1: float = 3.0, t2: float = 2.0, t_web: float = 4.0
) -> ClosedSection:
    """Build unsymmetric two-cell section with unequal cell widths and wall thicknesses."""
    n1 = Node(0.0, 0.0)
    n2 = Node(b1, 0.0)
    n3 = Node(b1 + b2, 0.0)
    n4 = Node(b1 + b2, h)
    n5 = Node(b1, h)
    n6 = Node(0.0, h)
    return ClosedSection([
        Segment(n1, n2, t=t1),
        Segment(n2, n3, t=t2),
        Segment(n3, n4, t=t2),
        Segment(n4, n5, t=t2),
        Segment(n5, n6, t=t1),
        Segment(n6, n1, t=t1),
        Segment(n2, n5, t=t_web),  # shared web
    ])


# ==============================================================================
# Benchmarks
# ==============================================================================

class TestClosedBenchmarks:
    """Execution of all 10 required benchmarks from ACTIVE_PHASE.md."""

    def test_benchmark_1_single_cell_box_bredt_j(self) -> None:
        """Benchmark 1: Single-cell rectangular box Bredt-Batho J analytical agreement."""
        a, b, t = 120.0, 80.0, 3.0
        sec = make_rect_box(a, b, t)
        expected_J = oracle_single_cell_box_J(a, b, t)

        assert math.isclose(sec.J, expected_J, rel_tol=1e-12)

        # Single-cell reduction check: J = 4*A_m^2 / \oint (ds/t)
        Am = a * b
        oint_ds_over_t = 2.0 * (a + b) / t
        j_reduced = (4.0 * Am**2) / oint_ds_over_t
        assert math.isclose(sec.J, j_reduced, rel_tol=1e-12)

    def test_benchmark_2_doubly_symmetric_box_shear_center(self) -> None:
        """Benchmark 2: Doubly symmetric rectangular box shear center S = C."""
        a, b, t = 150.0, 100.0, 4.0
        sec = make_rect_box(a, b, t)
        cx, cy = sec.centroid
        xs, ys = sec.shear_center
        ex, ey = sec.sc_offset

        assert math.isclose(ex, 0.0, abs_tol=1e-11)
        assert math.isclose(ey, 0.0, abs_tol=1e-11)
        assert math.isclose(xs, cx, abs_tol=1e-11)
        assert math.isclose(ys, cy, abs_tol=1e-11)

    def test_benchmark_3_rectangular_box_warping_constant(self) -> None:
        """Benchmark 3: Rectangular closed box C_w independent analytical agreement."""
        a, b, t = 120.0, 60.0, 2.5
        sec = make_rect_box(a, b, t)
        expected_Cw = oracle_single_cell_box_Cw(a, b, t)

        assert math.isclose(sec.Cw, expected_Cw, rel_tol=1e-11)

        # Special case: Square box a = b must give C_w = 0
        square_sec = make_rect_box(80.0, 80.0, 3.0)
        expected_square_Cw = oracle_single_cell_box_Cw(80.0, 80.0, 3.0)
        assert expected_square_Cw == 0.0
        assert math.isclose(square_sec.Cw, 0.0, abs_tol=1e-10)

    def test_benchmark_4_symmetric_two_cell_box(self) -> None:
        """Benchmark 4: Symmetric two-cell box analytical H, J, and zero shared web flow."""
        b, h, t = 50.0, 30.0, 2.0
        sec = make_twocell_box(b, h, t)
        expected_J = oracle_symmetric_twocell_box_J(b, h, t)

        assert math.isclose(sec.J, expected_J, rel_tol=1e-12)

        # Analytical H check:
        expected_H = (1.0 / t) * np.array(
            [[2.0 * (b + h), -h], [-h, 2.0 * (b + h)]], dtype=float
        )
        np.testing.assert_allclose(sec.H, expected_H, rtol=1e-12, atol=1e-12)

        # Pure torsion cell flows phi_1 and phi_2 must be equal
        t_res = sec.torsion_properties()
        phi1, phi2 = t_res.phi[0], t_res.phi[1]
        assert math.isclose(phi1, phi2, rel_tol=1e-12)

        # Shared web (index 6) must carry ZERO net torsional membrane flow: F_shared = 0
        F_shared = t_res.F[6]
        assert math.isclose(F_shared, 0.0, abs_tol=1e-12)

    def test_benchmark_5_unsymmetric_closed_cell(self) -> None:
        """Benchmark 5: Unsymmetric closed cell exercises nontrivial shear center and residual torque consistency."""
        sec = make_unsymmetric_box(100.0, 60.0, t_bot=5.0, t_top=1.5, t_left=2.0, t_right=6.0)
        ex, ey = sec.sc_offset
        # Because wall thicknesses are unsymmetric, shear center must NOT coincide with centroid
        assert abs(ex) > 1e-3
        assert abs(ey) > 1e-3

        sc_res = sec.compute_shear_center()
        # Verify combined load residual torque consistency:
        for vx, vy in [(500.0, 300.0), (-400.0, 800.0)]:
            load = ShearLoad(vx, vy)
            res_torque = sc_res.residual_torque(load)
            assert math.isclose(res_torque, 0.0, abs_tol=1e-10)

    def test_benchmark_6_unsymmetric_multicell(self) -> None:
        """Benchmark 6: Unsymmetric multi-cell with unequal cell areas and thicknesses."""
        sec = make_unsymmetric_twocell(80.0, 40.0, 50.0, t1=3.0, t2=1.5, t_web=4.0)

        # Cell flows under pure torsion must be unequal:
        t_res = sec.torsion_properties()
        assert not math.isclose(t_res.phi[0], t_res.phi[1], rel_tol=1e-3)
        assert sec.J > 0.0
        assert math.isfinite(sec.Cw)

        # All cycle residuals for warping coordinate derivative must vanish: \oint d(omega) = 0
        np.testing.assert_allclose(t_res.cycle_residuals, [0.0, 0.0], atol=1e-10)

        # Transverse shear compatibility and resultant recovery:
        flow_res = sec.calculate_shear_flow(600.0, -450.0)
        np.testing.assert_allclose(flow_res.recovered_resultant, [600.0, -450.0], atol=1e-10)
        np.testing.assert_allclose(flow_res.compatibility_residuals, [0.0, 0.0], atol=1e-10)

    def test_benchmark_7_cut_independence_single_and_multicell(self) -> None:
        """Benchmark 7: Cut set and cut location independence for single and two-cell sections."""
        # Single cell:
        sec1 = make_rect_box(100.0, 50.0, 2.0)
        res1_a = calculate_closed_shear_flow(sec1, 300.0, 400.0, cut_param=0.2, spanning_tree_edges=[0, 1, 2])
        res1_b = calculate_closed_shear_flow(sec1, 300.0, 400.0, cut_param=0.8, spanning_tree_edges=[1, 2, 3])
        assert math.isclose(res1_a.torque, res1_b.torque, rel_tol=1e-10, abs_tol=1e-10)
        for i in range(len(sec1.segments)):
            assert math.isclose(
                res1_a.segment_flows[i].q_mid, res1_b.segment_flows[i].q_mid, rel_tol=1e-10, abs_tol=1e-10
            )

        # Multi-cell:
        sec2 = make_twocell_box(60.0, 40.0, 3.0)
        res2_a = calculate_closed_shear_flow(sec2, 500.0, 200.0, cut_param=0.3)
        res2_b = calculate_closed_shear_flow(sec2, 500.0, 200.0, cut_param=0.7)
        assert math.isclose(res2_a.torque, res2_b.torque, rel_tol=1e-10, abs_tol=1e-10)
        for i in range(len(sec2.segments)):
            assert math.isclose(
                res2_a.segment_flows[i].q_mid, res2_b.segment_flows[i].q_mid, rel_tol=1e-10, abs_tol=1e-10
            )

    def test_benchmark_8_large_translation_invariance(self) -> None:
        """Benchmark 8: 1e12 coordinate offset leaves relative shear center, J, omega, and C_w unchanged."""
        base_sec = make_unsymmetric_box(100.0, 60.0, 4.0, 2.0, 3.0, 5.0)
        shift = 1e12
        trans_sec = base_sec.translated(shift, shift)

        # Centroid-relative shear-center offset:
        assert math.isclose(base_sec.sc_offset[0], trans_sec.sc_offset[0], rel_tol=1e-10, abs_tol=1e-10)
        assert math.isclose(base_sec.sc_offset[1], trans_sec.sc_offset[1], rel_tol=1e-10, abs_tol=1e-10)

        # Absolute coordinates shift by exactly 1e12:
        assert math.isclose(trans_sec.shear_center[0], base_sec.shear_center[0] + shift, rel_tol=1e-11)
        assert math.isclose(trans_sec.shear_center[1], base_sec.shear_center[1] + shift, rel_tol=1e-11)

        # Torsion constant J:
        assert math.isclose(base_sec.J, trans_sec.J, rel_tol=1e-11)

        # Warping constant C_w:
        assert math.isclose(base_sec.Cw, trans_sec.Cw, rel_tol=1e-10)

        # Normalized sectorial coordinates:
        np.testing.assert_allclose(
            base_sec.sectorial_coordinates, trans_sec.sectorial_coordinates, rtol=1e-10, atol=1e-10
        )

        # Shear flow field:
        res_base = base_sec.calculate_shear_flow(100.0, 200.0)
        res_trans = trans_sec.calculate_shear_flow(100.0, 200.0)
        for i in range(len(base_sec.segments)):
            assert math.isclose(
                res_base.segment_flows[i].q_mid, res_trans.segment_flows[i].q_mid, rel_tol=1e-10, abs_tol=1e-10
            )

    def test_benchmark_9_rotation_covariance(self) -> None:
        """Benchmark 9: Rigid in-plane rotation by nontrivial angle (37 deg)."""
        theta = math.radians(37.0)
        base_sec = make_unsymmetric_box(100.0, 60.0, 4.0, 2.0, 3.0, 5.0)
        rot_sec = base_sec.rotated(theta)

        # 1. Scalar J and C_w must be strictly invariant:
        assert math.isclose(rot_sec.J, base_sec.J, rel_tol=1e-11)
        assert math.isclose(rot_sec.Cw, base_sec.Cw, rel_tol=1e-10)

        # 2. Shear center offset e_s must rotate covariantly:
        # e_s' = R(theta) * e_s
        ex, ey = base_sec.sc_offset
        cos_t, sin_t = math.cos(theta), math.sin(theta)
        ex_expected = ex * cos_t - ey * sin_t
        ey_expected = ex * sin_t + ey * cos_t
        assert math.isclose(rot_sec.sc_offset[0], ex_expected, rel_tol=1e-10, abs_tol=1e-10)
        assert math.isclose(rot_sec.sc_offset[1], ey_expected, rel_tol=1e-10, abs_tol=1e-10)

        # 3. Normalized node omega values must be invariant:
        np.testing.assert_allclose(
            rot_sec.sectorial_coordinates, base_sec.sectorial_coordinates, rtol=1e-10, atol=1e-10
        )

    def test_benchmark_10_uniform_scaling(self) -> None:
        """Benchmark 10: Scaling all coordinates and wall thicknesses by lambda > 0.

        Scaling laws (ACTIVE_PHASE.md):
        e_s' = lambda * e_s
        J' = lambda^4 * J
        omega' = lambda^2 * omega
        C_w' = lambda^6 * C_w
        """
        scale_fac = 2.5
        base_sec = make_unsymmetric_box(100.0, 60.0, 4.0, 2.0, 3.0, 5.0)

        scaled_segs = [
            Segment(
                p1=Node(seg.p1.x * scale_fac, seg.p1.y * scale_fac),
                p2=Node(seg.p2.x * scale_fac, seg.p2.y * scale_fac),
                t=seg.t * scale_fac,
            )
            for seg in base_sec.segments
        ]
        scaled_sec = ClosedSection(scaled_segs)

        # 1. e_s' = lambda * e_s
        assert math.isclose(scaled_sec.sc_offset[0], base_sec.sc_offset[0] * scale_fac, rel_tol=1e-11)
        assert math.isclose(scaled_sec.sc_offset[1], base_sec.sc_offset[1] * scale_fac, rel_tol=1e-11)

        # 2. J' = lambda^4 * J
        assert math.isclose(scaled_sec.J, base_sec.J * (scale_fac**4), rel_tol=1e-11)

        # 3. omega' = lambda^2 * omega
        np.testing.assert_allclose(
            scaled_sec.sectorial_coordinates,
            np.array(base_sec.sectorial_coordinates) * (scale_fac**2),
            rtol=1e-11,
        )

        # 4. C_w' = lambda^6 * C_w
        assert math.isclose(scaled_sec.Cw, base_sec.Cw * (scale_fac**6), rel_tol=1e-11)

    def test_theoretical_relation_4_A_invH_A(self) -> None:
        """Criterion 26: Verify theoretical identity J = 4 * A^T * H^(-1) * A on multi-cell section."""
        sec = make_twocell_box(60.0, 40.0, 3.0)
        A_vec = np.array([c.area for c in sec.cells], dtype=float)
        H_inv = np.linalg.inv(sec.H)
        j_theoretical = float(4.0 * A_vec.T @ H_inv @ A_vec)
        assert math.isclose(sec.J, j_theoretical, rel_tol=1e-12)

    def test_root_node_invariance_for_warping(self) -> None:
        """Criterion 34: Changing graph root node must not change normalized omega or C_w."""
        sec = make_unsymmetric_box(100.0, 60.0, 4.0, 2.0, 3.0, 5.0)
        res_root0 = sec.torsion_properties(root_node_idx=0)
        res_root1 = sec.torsion_properties(root_node_idx=1)
        res_root2 = sec.torsion_properties(root_node_idx=2)

        assert math.isclose(res_root0.Cw, res_root1.Cw, rel_tol=1e-12)
        assert math.isclose(res_root0.Cw, res_root2.Cw, rel_tol=1e-12)
        np.testing.assert_allclose(res_root0.node_omega, res_root1.node_omega, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(res_root0.node_omega, res_root2.node_omega, rtol=1e-12, atol=1e-12)

    def test_endpoint_reversal_invariance(self) -> None:
        """Criterion 40: Reversing segment endpoints must not change physical J, C_w, or shear center."""
        base_sec = make_unsymmetric_box(100.0, 60.0, 4.0, 2.0, 3.0, 5.0)
        # Flip endpoints of segment 1 and 3
        flipped_segs = [
            base_sec.segments[0],
            Segment(base_sec.segments[1].p2, base_sec.segments[1].p1, t=base_sec.segments[1].t),
            base_sec.segments[2],
            Segment(base_sec.segments[3].p2, base_sec.segments[3].p1, t=base_sec.segments[3].t),
        ]
        flipped_sec = ClosedSection(flipped_segs)

        assert math.isclose(flipped_sec.J, base_sec.J, rel_tol=1e-12)
        assert math.isclose(flipped_sec.Cw, base_sec.Cw, rel_tol=1e-12)
        assert math.isclose(flipped_sec.sc_offset[0], base_sec.sc_offset[0], rel_tol=1e-12, abs_tol=1e-12)
        assert math.isclose(flipped_sec.sc_offset[1], base_sec.sc_offset[1], rel_tol=1e-12, abs_tol=1e-12)

    def test_reordering_invariance(self) -> None:
        """Criterion 41: Reordering input segments must not materially change physical results."""
        base_sec = make_unsymmetric_box(100.0, 60.0, 4.0, 2.0, 3.0, 5.0)
        # Permute segments: [2, 0, 3, 1]
        permuted_segs = [
            base_sec.segments[2],
            base_sec.segments[0],
            base_sec.segments[3],
            base_sec.segments[1],
        ]
        permuted_sec = ClosedSection(permuted_segs)

        assert math.isclose(permuted_sec.J, base_sec.J, rel_tol=1e-12)
        assert math.isclose(permuted_sec.Cw, base_sec.Cw, rel_tol=1e-12)
        assert math.isclose(permuted_sec.sc_offset[0], base_sec.sc_offset[0], rel_tol=1e-12, abs_tol=1e-12)
        assert math.isclose(permuted_sec.sc_offset[1], base_sec.sc_offset[1], rel_tol=1e-12, abs_tol=1e-12)

    def test_nonfinite_and_singular_failures(self) -> None:
        """Criteria 45-46: Non-finite inputs and singular matrices must fail explicitly."""
        sec = make_rect_box(100.0, 50.0, 2.0)
        from sectalix.exceptions import GeometryError, SingularSectionError

        with pytest.raises(GeometryError):
            sec.torsion_properties(safety_factor=float("nan"))
        with pytest.raises(GeometryError):
            sec.torsion_properties(safety_factor=-1.0)
        with pytest.raises(GeometryError):
            sec.torsion_properties(root_node_idx=999)

    def test_analytical_shear_flow_oracle(self) -> None:
        """Independent analytical shear flow oracle on single-cell box (b, h, t) under V_y.

        For a rectangular or square box of width b, height h, wall thickness t,
        loaded by transverse shear V_y, with chord on bottom flange cut at parameter xi_c in (0, 1),
        the redundant circulating flow satisfies the exact general analytical formula:
            expected_q0 = - (V_y * t * b * h * (1.0 - 2.0 * xi_c)) / (4.0 * I_x)

        Verified for both rectangular (b=100, h=50, t=2) and square (b=100, h=100, t=2)
        across multiple cut locations (xi_c = 1/6, 0.25, 0.45, 0.5) to rtol=1e-12, atol=0.0.
        """
        # 1. Rectangular box: b=100.0, h=50.0, t=2.0
        b, h, t = 100.0, 50.0, 2.0
        vy = 1000.0
        pts_rect = [(-b / 2.0, -h / 2.0), (b / 2.0, -h / 2.0), (b / 2.0, h / 2.0), (-b / 2.0, h / 2.0)]
        nodes_rect = [Node(pts_rect[k][0], pts_rect[k][1]) for k in range(4)]
        segs_rect = [
            Segment(nodes_rect[0], nodes_rect[1], t=t),
            Segment(nodes_rect[1], nodes_rect[2], t=t),
            Segment(nodes_rect[2], nodes_rect[3], t=t),
            Segment(nodes_rect[3], nodes_rect[0], t=t),
        ]
        sec_rect = ClosedSection(segs_rect)
        Ix_rect_oracle = t * (h**2) * (b / 2.0 + h / 6.0)

        for xi_c in [1.0 / 6.0, 0.25, 0.45, 0.5]:
            res = calculate_closed_shear_flow(
                sec_rect, vx=0.0, vy=vy, cut_param=xi_c, spanning_tree_edges=[1, 2, 3]
            )
            expected_q0 = -(vy * t * b * h * (1.0 - 2.0 * xi_c)) / (4.0 * Ix_rect_oracle)
            actual_q0 = res.redundant_cell_flows[0]
            atol_val = 1e-15 if expected_q0 == 0.0 else 0.0
            assert np.isclose(actual_q0, expected_q0, rtol=1e-12, atol=atol_val)

        # 2. Square box: b=100.0, h=100.0, t=2.0
        b_sq, h_sq = 100.0, 100.0
        pts_sq = [(-b_sq / 2.0, -h_sq / 2.0), (b_sq / 2.0, -h_sq / 2.0), (b_sq / 2.0, h_sq / 2.0), (-b_sq / 2.0, h_sq / 2.0)]
        nodes_sq = [Node(pts_sq[k][0], pts_sq[k][1]) for k in range(4)]
        segs_sq = [
            Segment(nodes_sq[0], nodes_sq[1], t=t),
            Segment(nodes_sq[1], nodes_sq[2], t=t),
            Segment(nodes_sq[2], nodes_sq[3], t=t),
            Segment(nodes_sq[3], nodes_sq[0], t=t),
        ]
        sec_sq = ClosedSection(segs_sq)
        Ix_sq_oracle = t * (h_sq**2) * (b_sq / 2.0 + h_sq / 6.0)

        for xi_c in [1.0 / 6.0, 0.25, 0.45, 0.5]:
            res = calculate_closed_shear_flow(
                sec_sq, vx=0.0, vy=vy, cut_param=xi_c, spanning_tree_edges=[1, 2, 3]
            )
            expected_q0 = -(vy * t * b_sq * h_sq * (1.0 - 2.0 * xi_c)) / (4.0 * Ix_sq_oracle)
            actual_q0 = res.redundant_cell_flows[0]
            if np.isclose(xi_c, 1.0 / 6.0):
                # Verify exact expected value -2.5
                assert np.isclose(expected_q0, -2.5, rtol=1e-12, atol=0.0)
            atol_val = 1e-15 if expected_q0 == 0.0 else 0.0
            assert np.isclose(actual_q0, expected_q0, rtol=1e-12, atol=atol_val)

    def test_concave_multicell_l_shape(self) -> None:
        """Concave multi-cell section test: 3-cell L-shaped box section."""
        w, h, t = 40.0, 40.0, 2.0
        n = {
            "00": Node(0, 0),
            "w0": Node(w, 0),
            "2w0": Node(2 * w, 0),
            "2wh": Node(2 * w, h),
            "wh": Node(w, h),
            "w2h": Node(w, 2 * h),
            "02h": Node(0, 2 * h),
            "0h": Node(0, h),
        }
        segs = [
            Segment(n["00"], n["w0"], t=t),
            Segment(n["w0"], n["2w0"], t=t),
            Segment(n["2w0"], n["2wh"], t=t),
            Segment(n["2wh"], n["wh"], t=t),
            Segment(n["wh"], n["w2h"], t=t),
            Segment(n["w2h"], n["02h"], t=t),
            Segment(n["02h"], n["0h"], t=t),
            Segment(n["0h"], n["00"], t=t),
            Segment(n["w0"], n["wh"], t=t),  # internal shared wall between cell 0 and 1
            Segment(n["0h"], n["wh"], t=t),  # internal shared wall between cell 0 and 2
        ]
        sec = ClosedSection(segs)
        assert sec.cell_count == 3
        assert len(sec.cell_topology.cells) == 3

        tors = sec.torsion_properties()
        assert tors.J > 0.0 and math.isfinite(tors.J)
        assert tors.Cw >= 0.0 and math.isfinite(tors.Cw)

        flow = sec.calculate_shear_flow(150.0, 250.0)
        np.testing.assert_allclose(flow.recovered_resultant, [150.0, 250.0], rtol=1e-11, atol=1e-10)
        assert flow.max_compatibility_residual < 1e-10

    def test_adversarial_huge_twocell_j(self) -> None:
        """Adversarial benchmark: Huge two-cell box (b=h=1e10, t=1e-300) with L/t > 1e308.

        H flexibility matrix must not overflow or become NaN.
        Analytical J = (8/3) * 1e-270 must be preserved to rtol=1e-12, atol=0.0.
        """
        b, h, t = 1e10, 1e10, 1e-300
        n1, n2, n3 = Node(0, 0), Node(b, 0), Node(2 * b, 0)
        n4, n5, n6 = Node(2 * b, h), Node(b, h), Node(0, h)
        segs = [
            Segment(n1, n2, t=t),
            Segment(n2, n3, t=t),
            Segment(n3, n4, t=t),
            Segment(n4, n5, t=t),
            Segment(n5, n6, t=t),
            Segment(n6, n1, t=t),
            Segment(n2, n5, t=t),
        ]
        sec = ClosedSection(segs)
        expected_J = (8.0 / 3.0) * 1e-270
        assert np.isclose(sec.J, expected_J, rtol=1e-12, atol=0.0)

    def test_adversarial_huge_singlecell_j(self) -> None:
        """Adversarial benchmark: Huge single-cell box (a=b=1e100, t=1e-220).

        Analytical J = 2*a^2*b^2*t / (a + b) = 1e80 must be exact to rtol=1e-12, atol=0.0.
        """
        a, b, t = 1e100, 1e100, 1e-220
        n1, n2, n3, n4 = Node(0, 0), Node(a, 0), Node(a, b), Node(0, b)
        sec = ClosedSection([
            Segment(n1, n2, t=t),
            Segment(n2, n3, t=t),
            Segment(n3, n4, t=t),
            Segment(n4, n1, t=t),
        ])
        expected_J = 1e80
        assert np.isclose(sec.J, expected_J, rtol=1e-12, atol=0.0)

    def test_micro_box_valid(self) -> None:
        """Adversarial benchmark: Micro-box (a=b=1e-8, t=1e-10).

        Enclosed area is 1e-16 (below 1e-15 dimensional floor).
        Must be accepted as a valid cell with finite positive J.
        """
        a, b, t = 1e-8, 1e-8, 1e-10
        n1, n2, n3, n4 = Node(0, 0), Node(a, 0), Node(a, b), Node(0, b)
        sec = ClosedSection([
            Segment(n1, n2, t=t),
            Segment(n2, n3, t=t),
            Segment(n3, n4, t=t),
            Segment(n4, n1, t=t),
        ])
        assert sec.cell_count == 1
        assert sec.J > 0.0
        assert math.isfinite(sec.J)
        expected_J = 1e-34
        assert np.isclose(sec.J, expected_J, rtol=1e-12, atol=0.0)

    def test_adversarial_huge_b_shear_flow(self) -> None:
        """Codex adversarial counter-example: Extreme basic compatibility vector b.

        L = 1e100, t = 1e-220, xi_c = 0.45.
        1. At Vy = 1e89: segment compliance integrals exceed float64 (~5e308),
           so calculate_closed_shear_flow must explicitly raise OverflowError without
           silent inf leakage or fake DBL_MAX clamping.
        2. At Vy = 2.0e88: all segment compliance integrals remain strictly finite
           (<= 1e308 < DBL_MAX), b = 3.0e307, q0 = -7.5e-14, with no values equal to DBL_MAX.
        """
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            L = 1e100
            t = 1e-220
            xi_c = 0.45

            pts = [(-L / 2.0, -L / 2.0), (L / 2.0, -L / 2.0), (L / 2.0, L / 2.0), (-L / 2.0, L / 2.0)]
            nodes = [Node(pts[k][0], pts[k][1]) for k in range(4)]
            segs = [
                Segment(nodes[0], nodes[1], t=t),
                Segment(nodes[1], nodes[2], t=t),
                Segment(nodes[2], nodes[3], t=t),
                Segment(nodes[3], nodes[0], t=t),
            ]
            sec = ClosedSection(segs)

            # 1. Vy = 1e89 must explicitly raise OverflowError
            vy_overflow = 1e89
            with pytest.raises(
                OverflowError,
                match="Segment transverse-shear compatibility integral overflows IEEE-754 float64",
            ):
                calculate_closed_shear_flow(
                    sec, vx=0.0, vy=vy_overflow, cut_param=xi_c, spanning_tree_edges=[1, 2, 3]
                )

            # 2. Sub-load level Vy = 2.0e88 where all segment integrals are strictly finite
            vy_sub = 2.0e88
            res = calculate_closed_shear_flow(
                sec, vx=0.0, vy=vy_sub, cut_param=xi_c, spanning_tree_edges=[1, 2, 3]
            )

            assert np.isclose(res.redundant_cell_flows[0], -7.5e-14, rtol=1e-12, atol=0.0)
            assert np.isclose(res.b_vector[0], 3.0e307, rtol=1e-12, atol=0.0)
            np.testing.assert_allclose(res.recovered_resultant, [0.0, vy_sub], atol=1e-11 * vy_sub, rtol=1e-11)

            # Verify all segment integrals are strictly finite and none clamped to DBL_MAX
            for sf in res.segment_flows:
                assert math.isfinite(sf.integral_q_over_t)
                assert sf.integral_q_over_t != sys.float_info.max
                assert sf.integral_q_over_t != -sys.float_info.max

    def test_flexibility_matrix_overflow_raises_error(self) -> None:
        """Accessing .H property when elements exceed float64 must explicitly raise OverflowError."""
        L = 1e100
        t = 1e-220
        pts = [(-L / 2.0, -L / 2.0), (L / 2.0, -L / 2.0), (L / 2.0, L / 2.0), (-L / 2.0, L / 2.0)]
        nodes = [Node(pts[k][0], pts[k][1]) for k in range(4)]
        segs = [
            Segment(nodes[0], nodes[1], t=t),
            Segment(nodes[1], nodes[2], t=t),
            Segment(nodes[2], nodes[3], t=t),
            Segment(nodes[3], nodes[0], t=t),
        ]
        sec = ClosedSection(segs)

        with pytest.raises(OverflowError, match="Cell flexibility matrix H overflows IEEE-754 float64"):
            _ = sec.H

        with pytest.raises(OverflowError, match="Cell flexibility matrix H overflows IEEE-754 float64"):
            _ = sec.cell_topology.H

    def test_h_finite_at_exponent_boundary(self) -> None:
        """Triangle cell with one edge t=2^-1023 produces e_max_H=1024, but H is finite and must not raise."""
        t0 = math.ldexp(1.0, -1023)
        segs = [
            Segment(Node(0.0, 0.0), Node(1.0, 0.0), t=t0),
            Segment(Node(1.0, 0.0), Node(0.0, 1.0), t=1.0),
            Segment(Node(0.0, 1.0), Node(0.0, 0.0), t=1.0),
        ]
        sec = ClosedSection(segs)
        assert sec.cell_topology.e_max_H == 1024

        expected_H_val = 8.98846567431158e307
        H_sec = sec.H
        H_topo = sec.cell_topology.H

        assert np.isfinite(H_sec[0, 0])
        assert np.isfinite(H_topo[0, 0])
        assert np.isclose(H_sec[0, 0], expected_H_val, rtol=1e-12, atol=0.0)
        assert np.isclose(H_topo[0, 0], expected_H_val, rtol=1e-12, atol=0.0)

    def test_extreme_load_shear_flow_raises_overflow(self) -> None:
        """Extreme shear load on extreme box section must explicitly raise OverflowError without leaking inf."""
        L = 1e100
        t = 1e-220
        vy = 1e100
        pts = [(-L / 2.0, -L / 2.0), (L / 2.0, -L / 2.0), (L / 2.0, L / 2.0), (-L / 2.0, L / 2.0)]
        nodes = [Node(pts[k][0], pts[k][1]) for k in range(4)]
        segs = [
            Segment(nodes[0], nodes[1], t=t),
            Segment(nodes[1], nodes[2], t=t),
            Segment(nodes[2], nodes[3], t=t),
            Segment(nodes[3], nodes[0], t=t),
        ]
        sec = ClosedSection(segs)

        with pytest.raises(
            OverflowError,
            match="Physical compatibility vector b or compatibility residuals overflow",
        ):
            calculate_closed_shear_flow(sec, vx=0.0, vy=vy)


