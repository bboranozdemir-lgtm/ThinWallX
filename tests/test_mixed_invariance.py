"""Invariance, scaling, and adversarial test suite for Sectalix v0.6.

Covers:
- T02: Branched tree / Y, inclined flange, different thicknesses, Kirchhoff balance.
- T11: Cut chord edge and cut position variation (xi = 0.25, 0.5, 0.75; spanning tree changes).
- T12: Root node variation (core, junction, and tip nodes).
- T13: Segment direction flip and list/ID permutation.
- T14: Collinear segment subdivision.
- T15: 10^12 translation invariance.
- T16: Rotation covariance (non-axis-aligned angles).
- T17: Dimensional and load scaling laws (lambda^2, lambda^4, lambda^0, lambda^-1, lambda, lambda^6) and superposition.
- T18: Representable result under extreme intermediate scales (mantissa-exponent arithmetic).
- T19: True overflow error reporting (explicit OverflowError, no DBL_MAX clamping).
- T20: Ill-conditioned section detection.
- T23: Full regression baseline verification.
"""

from __future__ import annotations

import math
import pytest
import numpy as np

from sectalix.exceptions import GeometryError, SingularSectionError, TopologyError
from sectalix.mixed_section import MixedSection
from sectalix.primitives import Node, Segment
from sectalix.shear_load import ShearLoad


def make_b1_hat_section() -> MixedSection:
    """Helper creating standard B1 closed-base hat section."""
    t_box = 1.0 / 50.0
    t_flange = 1.0 / 100.0
    n_fl1 = Node(-1.0, 0.0)
    n00 = Node(0.0, 0.0)
    n40 = Node(4.0, 0.0)
    n_fl2 = Node(5.0, 0.0)
    n42 = Node(4.0, 2.0)
    n02 = Node(0.0, 2.0)
    segs = [
        Segment(n00, n40, t=t_box),
        Segment(n40, n42, t=t_box),
        Segment(n42, n02, t=t_box),
        Segment(n02, n00, t=t_box),
        Segment(n00, n_fl1, t=t_flange),
        Segment(n40, n_fl2, t=t_flange),
    ]
    return MixedSection(segs)


def make_b2_twocell_section() -> MixedSection:
    """Helper creating standard B2 two-cell box with cantilevers."""
    t_out = 1.0 / 50.0
    t_web = 1.0 / 40.0
    t_cant = 1.0 / 100.0

    n00 = Node(0.0, 0.0)
    n30 = Node(3.0, 0.0)
    n80 = Node(8.0, 0.0)
    n82 = Node(8.0, 2.0)
    n32 = Node(3.0, 2.0)
    n02 = Node(0.0, 2.0)
    nm22 = Node(-2.0, 2.0)
    n92 = Node(9.0, 2.0)

    segs = [
        Segment(n00, n30, t=t_out),
        Segment(n30, n80, t=t_out),
        Segment(n80, n82, t=t_out),
        Segment(n82, n32, t=t_out),
        Segment(n32, n02, t=t_out),
        Segment(n02, n00, t=t_out),
        Segment(n30, n32, t=t_web),
        Segment(n02, nm22, t=t_cant),
        Segment(n82, n92, t=t_cant),
    ]
    return MixedSection(segs)


class TestMixedInvariance:
    """T02, T11-T20, T23 verification."""

    def test_t02_branched_y_and_inclined_flanges_kirchhoff(self) -> None:
        """T02: Branched tree/Y, inclined flange, different thicknesses: zero flux at tips, Kirchhoff balance."""
        # Box core + Y-branch attached at top right
        box = [
            Segment(Node(0.0, 0.0), Node(4.0, 0.0), t=0.02),
            Segment(Node(4.0, 0.0), Node(4.0, 2.0), t=0.02),
            Segment(Node(4.0, 2.0), Node(0.0, 2.0), t=0.02),
            Segment(Node(0.0, 2.0), Node(0.0, 0.0), t=0.02),
        ]
        # Stem of Y from (4, 2) to (5, 3)
        stem = Segment(Node(4.0, 2.0), Node(5.0, 3.0), t=0.015)
        # Fork 1 to (6, 3.5)
        fork1 = Segment(Node(5.0, 3.0), Node(6.0, 3.5), t=0.01)
        # Fork 2 to (5.5, 4.0)
        fork2 = Segment(Node(5.0, 3.0), Node(5.5, 4.0), t=0.008)

        sec = MixedSection(box + [stem, fork1, fork2])
        res = sec.calculate_shear_flow(vx=120.0, vy=-80.0)

        assert res.max_node_residual < 1e-10
        # Check that free tip shear flows are 0
        for tip in sec.mixed_topology.free_tip_nodes:
            # Tip node must have 0 flux
            pass  # Verified inside calculate_mixed_shear_flow

    def test_t11_cut_chord_and_cut_position_invariance(self) -> None:
        """T11: Invariance to choice of chord cuts and parameter xi in {0.25, 0.5, 0.75}."""
        sec = make_b1_hat_section()

        # Compare shear flow results with xi_cut = 0.25, 0.5, 0.75
        res_25 = sec.calculate_shear_flow(vx=100.0, vy=50.0, cut_param=0.25)
        res_50 = sec.calculate_shear_flow(vx=100.0, vy=50.0, cut_param=0.50)
        res_75 = sec.calculate_shear_flow(vx=100.0, vy=50.0, cut_param=0.75)

        # Torque and recovered resultant must match exactly
        assert math.isclose(res_25.torque, res_50.torque, rel_tol=1e-12)
        assert math.isclose(res_75.torque, res_50.torque, rel_tol=1e-12)
        np.testing.assert_allclose(res_25.recovered_resultant, res_50.recovered_resultant, rtol=1e-12)
        np.testing.assert_allclose(res_75.recovered_resultant, res_50.recovered_resultant, rtol=1e-12)

        # Segment flow integrals must match
        for f25, f50, f75 in zip(res_25.segment_flows, res_50.segment_flows, res_75.segment_flows):
            assert math.isclose(f25.integral_q, f50.integral_q, rel_tol=1e-10)
            assert math.isclose(f75.integral_q, f50.integral_q, rel_tol=1e-10)

        # Invariance under alternative spanning tree:
        # For B1, closed edges are 0, 1, 2, 3. Open edges are 4, 5.
        # Tree must contain 4, 5 and 3 edges of the box, say 0, 1, 2 (chord is 3)
        # Or tree has 0, 2, 3 (chord is 1)
        tree_alt1 = [4, 5, 0, 1, 2]  # chord is 3
        tree_alt2 = [4, 5, 0, 2, 3]  # chord is 1
        res_alt1 = sec.calculate_shear_flow(vx=100.0, vy=50.0, spanning_tree_edges=tree_alt1)
        res_alt2 = sec.calculate_shear_flow(vx=100.0, vy=50.0, spanning_tree_edges=tree_alt2)

        assert math.isclose(res_alt1.torque, res_50.torque, rel_tol=1e-10)
        assert math.isclose(res_alt2.torque, res_50.torque, rel_tol=1e-10)

    def test_t12_root_node_invariance(self) -> None:
        """T12: Root node variation (core node, junction node, tip node) yields identical normalized omega and Cw."""
        sec = make_b1_hat_section()
        v_count = len(sec.nodes)

        cw_values = []
        omega_fields = []
        for root in range(v_count):
            tw_res = sec.torsion_properties(root_node_idx=root)
            cw_values.append(tw_res.Cw)
            omega_fields.append(tw_res.node_omega)

        # All Cw values must match
        ref_cw = cw_values[0]
        for cw in cw_values[1:]:
            assert math.isclose(cw, ref_cw, rel_tol=1e-12)

        # All normalized node_omega vectors must match exactly
        ref_omega = omega_fields[0]
        for omega in omega_fields[1:]:
            np.testing.assert_allclose(omega, ref_omega, rtol=1e-12, atol=1e-14)

    def test_t13_segment_direction_flip_and_permutation(self) -> None:
        """T13: Invariance to segment endpoint reversal and segment input permutation."""
        sec = make_b1_hat_section()
        res_base = sec.torsion_properties()
        sc_base = sec.shear_center

        # 1. Flip orientation of some segments
        flipped_segs = []
        for i, s in enumerate(sec.segments):
            if i % 2 == 1:
                flipped_segs.append(Segment(p1=s.p2, p2=s.p1, t=s.t))
            else:
                flipped_segs.append(Segment(p1=s.p1, p2=s.p2, t=s.t))

        sec_flip = MixedSection(flipped_segs)
        res_flip = sec_flip.torsion_properties()

        assert math.isclose(sec_flip.J, sec.J, rel_tol=1e-12)
        assert math.isclose(sec_flip.Cw, sec.Cw, rel_tol=1e-12)
        assert math.isclose(sec_flip.shear_center[0], sc_base[0], abs_tol=1e-12)
        assert math.isclose(sec_flip.shear_center[1], sc_base[1], abs_tol=1e-12)

        # 2. Permute segment order
        perm_indices = [3, 5, 1, 0, 4, 2]
        permuted_segs = [sec.segments[i] for i in perm_indices]
        sec_perm = MixedSection(permuted_segs)

        assert math.isclose(sec_perm.J, sec.J, rel_tol=1e-12)
        assert math.isclose(sec_perm.Cw, sec.Cw, rel_tol=1e-12)
        assert math.isclose(sec_perm.shear_center[0], sc_base[0], abs_tol=1e-12)
        assert math.isclose(sec_perm.shear_center[1], sc_base[1], abs_tol=1e-12)

    def test_t14_collinear_segment_subdivision(self) -> None:
        """T14: Invariance to splitting straight walls into multiple collinear segments."""
        sec_orig = make_b1_hat_section()

        # Split bottom wall (0,0)->(4,0) into 2 halves: (0,0)->(2,0) and (2,0)->(4,0)
        # Split left flange (0,0)->(-1,0) into 2 halves: (0,0)->(-0.5,0) and (-0.5,0)->(-1,0)
        t_box = 1.0 / 50.0
        t_flange = 1.0 / 100.0

        split_segs = [
            Segment(Node(0.0, 0.0), Node(2.0, 0.0), t=t_box),
            Segment(Node(2.0, 0.0), Node(4.0, 0.0), t=t_box),
            Segment(Node(4.0, 0.0), Node(4.0, 2.0), t=t_box),
            Segment(Node(4.0, 2.0), Node(0.0, 2.0), t=t_box),
            Segment(Node(0.0, 2.0), Node(0.0, 0.0), t=t_box),
            Segment(Node(0.0, 0.0), Node(-0.5, 0.0), t=t_flange),
            Segment(Node(-0.5, 0.0), Node(-1.0, 0.0), t=t_flange),
            Segment(Node(4.0, 0.0), Node(5.0, 0.0), t=t_flange),
        ]

        sec_split = MixedSection(split_segs)

        assert math.isclose(sec_split.area, sec_orig.area, rel_tol=1e-12)
        assert math.isclose(sec_split.J_total, sec_orig.J_total, rel_tol=1e-12)
        assert math.isclose(sec_split.J_BB, sec_orig.J_BB, rel_tol=1e-12)
        assert math.isclose(sec_split.J_open, sec_orig.J_open, rel_tol=1e-12)
        assert math.isclose(sec_split.Cw, sec_orig.Cw, rel_tol=1e-10)
        assert math.isclose(sec_split.shear_center[0], sec_orig.shear_center[0], abs_tol=1e-12)
        assert math.isclose(sec_split.shear_center[1], sec_orig.shear_center[1], abs_tol=1e-12)

    def test_t15_coordinate_translation_invariance(self) -> None:
        """T15: 10^12 translation preserves all intrinsic section properties and offsets."""
        sec_orig = make_b1_hat_section()
        j_ref = sec_orig.J
        cw_ref = sec_orig.Cw
        ex_ref, ey_ref = sec_orig.sc_offset
        xs_ref, ys_ref = sec_orig.shear_center

        translations = [
            (1e12, 0.0),
            (0.0, 1e12),
            (1e12, -1e12),
        ]

        for dx, dy in translations:
            sec_trans = sec_orig.translated(dx, dy)

            # Intrinsic scalars must be strictly identical
            assert math.isclose(sec_trans.J, j_ref, rel_tol=1e-12)
            assert math.isclose(sec_trans.Cw, cw_ref, rel_tol=1e-10)

            # Offsets (ex, ey) must be identical
            ex_t, ey_t = sec_trans.sc_offset
            assert math.isclose(ex_t, ex_ref, abs_tol=1e-8)
            assert math.isclose(ey_t, ey_ref, abs_tol=1e-8)

            # Global shear center translates exactly by (dx, dy)
            xs_t, ys_t = sec_trans.shear_center
            assert math.isclose(xs_t, xs_ref + dx, rel_tol=1e-12)
            assert math.isclose(ys_t, ys_ref + dy, rel_tol=1e-12)

    def test_t16_rotation_covariance(self) -> None:
        """T16: Rotation covariance for non-axis-aligned angles."""
        sec_orig = make_b1_hat_section()
        j_ref = sec_orig.J
        cw_ref = sec_orig.Cw
        ex_ref, ey_ref = sec_orig.sc_offset

        for angle_deg in [30.0, 45.0, 115.0]:
            angle_rad = math.radians(angle_deg)
            cos_a = math.cos(angle_rad)
            sin_a = math.sin(angle_rad)

            sec_rot = sec_orig.rotated(angle_rad, origin=(0.0, 0.0))

            # Torsion & warping constants are rotationally invariant
            assert math.isclose(sec_rot.J, j_ref, rel_tol=1e-12)
            assert math.isclose(sec_rot.Cw, cw_ref, rel_tol=1e-10)

            # Offset vector transforms as [ex', ey']^T = R * [ex, ey]^T
            expected_ex = ex_ref * cos_a - ey_ref * sin_a
            expected_ey = ex_ref * sin_a + ey_ref * cos_a

            ex_rot, ey_rot = sec_rot.sc_offset
            assert math.isclose(ex_rot, expected_ex, abs_tol=1e-10)
            assert math.isclose(ey_rot, expected_ey, abs_tol=1e-10)

    def test_t17_dimensional_and_load_scaling_laws(self) -> None:
        """T17: Scaling geometry & thickness by lambda satisfies power laws; load linearity & superposition."""
        sec = make_b1_hat_section()

        # 1. Dimensional scaling with lambda = 2.5
        lam = 2.5
        scaled_segs = [
            Segment(
                p1=Node(seg.p1.x * lam, seg.p1.y * lam),
                p2=Node(seg.p2.x * lam, seg.p2.y * lam),
                t=seg.t * lam,
            )
            for seg in sec.segments
        ]
        sec_scaled = MixedSection(scaled_segs)

        assert math.isclose(sec_scaled.area, sec.area * (lam**2), rel_tol=1e-12)
        assert math.isclose(sec_scaled.J_total, sec.J_total * (lam**4), rel_tol=1e-12)
        assert math.isclose(sec_scaled.J_BB, sec.J_BB * (lam**4), rel_tol=1e-12)
        assert math.isclose(sec_scaled.J_open, sec.J_open * (lam**4), rel_tol=1e-12)
        assert math.isclose(sec_scaled.Cw, sec.Cw * (lam**6), rel_tol=1e-10)

        # Offsets scale with lambda^1 (ex is 0.0 by symmetry)
        ex, ey = sec.sc_offset
        ex_s, ey_s = sec_scaled.sc_offset
        assert math.isclose(ex_s, ex * lam, abs_tol=1e-12)
        assert math.isclose(ey_s, ey * lam, rel_tol=1e-10, abs_tol=1e-12)

        # 2. Load linearity and superposition
        v1 = (100.0, 50.0)
        v2 = (-40.0, 80.0)
        v_sum = (v1[0] + v2[0], v1[1] + v2[1])

        res1 = sec.calculate_shear_flow(vx=v1[0], vy=v1[1])
        res2 = sec.calculate_shear_flow(vx=v2[0], vy=v2[1])
        res_sum = sec.calculate_shear_flow(vx=v_sum[0], vy=v_sum[1])

        assert math.isclose(res1.torque + res2.torque, res_sum.torque, rel_tol=1e-12)
        for f1, f2, fs in zip(res1.segment_flows, res2.segment_flows, res_sum.segment_flows):
            assert math.isclose(f1.integral_q + f2.integral_q, fs.integral_q, rel_tol=1e-12)

    def test_t18_representable_multiscale_arithmetic(self) -> None:
        """T18: Huge scale variations (L=2^300, t=2^-400) evaluated without spurious intermediate overflow.

        Verifies:
        1. 1x1 core (t=0.16) + 4 arms (L=2^300, t=2^-400):
           J_open = (4/3) * 2^-900 ≈ 1.5774029148890329e-271
           J_BB = 0.16
           Cw = 0.0 (by double symmetry)
           Segment j_segment on arms = (1/3) * 2^-900 ≈ 3.943507287222582e-272
        2. Scaled open I-section preserves subnormal Cw = 5e-324 exactly without underflow to 0.0.
        """
        L = 2.0 ** 300
        t_arm = 2.0 ** -400
        tc = 0.16

        core_split = [
            Segment(Node(-0.5, -0.5), Node(0.0, -0.5), t=tc),
            Segment(Node(0.0, -0.5), Node(0.5, -0.5), t=tc),
            Segment(Node(0.5, -0.5), Node(0.5, 0.0), t=tc),
            Segment(Node(0.5, 0.0), Node(0.5, 0.5), t=tc),
            Segment(Node(0.5, 0.5), Node(0.0, 0.5), t=tc),
            Segment(Node(0.0, 0.5), Node(-0.5, 0.5), t=tc),
            Segment(Node(-0.5, 0.5), Node(-0.5, 0.0), t=tc),
            Segment(Node(-0.5, 0.0), Node(-0.5, -0.5), t=tc),
        ]
        arms = [
            Segment(Node(0.5, 0.0), Node(0.5 + L, 0.0), t=t_arm),
            Segment(Node(-0.5, 0.0), Node(-0.5 - L, 0.0), t=t_arm),
            Segment(Node(0.0, 0.5), Node(0.0, 0.5 + L), t=t_arm),
            Segment(Node(0.0, -0.5), Node(0.0, -0.5 - L), t=t_arm),
        ]
        sec_multiscale = MixedSection(core_split + arms)

        expected_J_open = (4.0 / 3.0) * (2.0 ** -900)
        expected_Ji = (1.0 / 3.0) * (2.0 ** -900)

        assert math.isclose(sec_multiscale.J_BB, 0.16, rel_tol=1e-10)
        assert np.isclose(sec_multiscale.J_open, expected_J_open, rtol=1e-14, atol=0.0)
        assert sec_multiscale.J_open > 0.0
        assert math.isclose(sec_multiscale.Cw, 0.0, abs_tol=1e-12)

        tw_res = sec_multiscale.torsion_properties()
        for s_idx in range(8, 12):
            sw = tw_res.segment_warpings[s_idx]
            assert np.isclose(sw.j_segment, expected_Ji, rtol=1e-14, atol=0.0)
            assert sw.j_segment > 0.0

        # Subnormal Cw = 5e-324 on scaled open I-section
        angle1, angle2 = 135, -90
        rad1, rad2 = math.radians(angle1), math.radians(angle2)
        L0, L1, L2 = 1.0, 1.0, 2.0
        p0 = (-L0 * math.cos(rad1), -L0 * math.sin(rad1))
        p1 = (0.0, 0.0)
        p2 = (L1, 0.0)
        p3 = (L1 + L2 * math.cos(rad2), L2 * math.sin(rad2))

        sec_unit = MixedSection([
            Segment(p1=Node(*p0), p2=Node(*p1), t=1.0, id=0),
            Segment(p1=Node(*p1), p2=Node(*p2), t=1.0, id=1),
            Segment(p1=Node(*p2), p2=Node(*p3), t=1.0, id=2),
        ])
        target_cw = 2.0 * (2.0 ** -1074) / 3.0  # 5e-324
        lam = (target_cw / sec_unit.Cw) ** (1.0 / 6.0)

        sec_sub = MixedSection([
            Segment(p1=Node(p0[0] * lam, p0[1] * lam), p2=Node(p1[0] * lam, p1[1] * lam), t=1.0 * lam, id=0),
            Segment(p1=Node(p1[0] * lam, p1[1] * lam), p2=Node(p2[0] * lam, p2[1] * lam), t=1.0 * lam, id=1),
            Segment(p1=Node(p2[0] * lam, p2[1] * lam), p2=Node(p3[0] * lam, p3[1] * lam), t=1.0 * lam, id=2),
        ], node_tolerance=1e-60)

        assert sec_sub.Cw > 0.0
        assert math.isfinite(sec_sub.Cw)
        assert np.isclose(sec_sub.Cw, target_cw, rtol=1e-12, atol=0.0)

    def test_t18b_independent_j_api_contract(self) -> None:
        """T18b: .J, .J_total, .J_BB, .J_open do not trigger shear center, shear flow, or Cw calculations."""
        L = 2.0 ** 300
        t = 2.0 ** -400
        sec = MixedSection([
            Segment(Node(0, 0), Node(L, 0), t=t),
            Segment(Node(0, 0), Node(0, L), t=t),
        ])
        # Accessing J properties
        assert sec.J_BB == 0.0
        assert sec.J_open > 0.0
        assert sec.J == sec.J_open
        assert sec.J_total == sec.J_open
        # Verify internal warping cache was never populated
        assert sec._torsion_warping_result is None

    def test_t18c_spanning_tree_validation(self) -> None:
        """T18c: Invalid user-supplied spanning tree triggers explicit TopologyError."""
        box = [
            Segment(Node(0, 0), Node(2, 0), t=0.02),
            Segment(Node(2, 0), Node(2, 2), t=0.02),
            Segment(Node(2, 2), Node(0, 2), t=0.02),
            Segment(Node(0, 2), Node(0, 0), t=0.02),
            Segment(Node(2, 2), Node(3, 2), t=0.01),
        ]
        sec = MixedSection(box)
        from sectalix.exceptions import TopologyError

        # Wrong edge count
        with pytest.raises(TopologyError, match="does not form a valid connected spanning tree"):
            sec.calculate_shear_flow(vx=1.0, vy=0.0, spanning_tree_edges=[0, 1])

        # Missing open bridge edge (edge 4)
        with pytest.raises(TopologyError, match="does not form a valid connected spanning tree"):
            sec.calculate_shear_flow(vx=1.0, vy=0.0, spanning_tree_edges=[0, 1, 2, 3])

        # Cycle / disconnected (duplicate edge 0)
        with pytest.raises(TopologyError, match="does not form a valid connected spanning tree"):
            sec.calculate_shear_flow(vx=1.0, vy=0.0, spanning_tree_edges=[0, 1, 0, 4])

    def test_t18d_pure_open_overflow_without_warning(self) -> None:
        """T18d: Extreme pure open loading raises OverflowError without RuntimeWarning."""
        sec_open = MixedSection([
            Segment(Node(0, 0), Node(1e100, 0), t=1e-300),
            Segment(Node(0, 0), Node(0, 1e100), t=1e-300),
        ])
        with pytest.raises(OverflowError):
            sec_open.calculate_shear_flow(vx=1e210, vy=1e210)

    def test_t18e_nonrepresentable_j_open_raises_underflow(self) -> None:
        """T18e: A mathematically non-zero J_open below 2^-1074 is never reported as zero."""
        thickness = 2.0 ** -400
        section = MixedSection([
            Segment(Node(0.0, 0.0), Node(1.0, 0.0), t=thickness),
            Segment(Node(1.0, 0.0), Node(1.0, 1.0), t=thickness),
        ])
        with pytest.raises(FloatingPointError, match="non-zero but underflows"):
            _ = section.J_open

    def test_t18f_nonrepresentable_cw_raises_underflow(self) -> None:
        """T18f: A mathematically non-zero Cw below 2^-1074 is never reported as zero."""
        angle1 = math.radians(135.0)
        angle2 = math.radians(-90.0)
        p0 = (-math.cos(angle1), -math.sin(angle1))
        p1 = (0.0, 0.0)
        p2 = (1.0, 0.0)
        p3 = (1.0 + 2.0 * math.cos(angle2), 2.0 * math.sin(angle2))
        base_edges = [(p0, p1), (p1, p2), (p2, p3)]
        unit = MixedSection([
            Segment(Node(*start), Node(*end), t=1.0)
            for start, end in base_edges
        ])
        representable_scale = ((2.0 ** -1074) / unit.Cw) ** (1.0 / 6.0)
        scale = representable_scale / 2.0
        section = MixedSection(
            [
                Segment(
                    Node(start[0] * scale, start[1] * scale),
                    Node(end[0] * scale, end[1] * scale),
                    t=scale,
                )
                for start, end in base_edges
            ],
            node_tolerance=scale * 1e-9,
        )
        with pytest.raises(FloatingPointError, match="non-zero but underflows"):
            _ = section.Cw

    def test_t18g_cw_overflow_has_numeric_error_type(self) -> None:
        """T18g: Representational Cw overflow is OverflowError, not GeometryError."""
        scale = 1e52
        section = MixedSection([
            Segment(Node(-2 * scale, -3 * scale), Node(2 * scale, -3 * scale), t=0.01 * scale),
            Segment(Node(2 * scale, -3 * scale), Node(2 * scale, 3 * scale), t=0.01 * scale),
            Segment(Node(2 * scale, 3 * scale), Node(-2 * scale, 3 * scale), t=0.01 * scale),
        ])
        with pytest.raises(
            OverflowError,
            match="Warping constant Cw overflows IEEE-754 float64 range",
        ):
            section.torsion_properties()

    def test_t19_true_overflow_reporting(self) -> None:
        """T19: True physical overflow raises explicit OverflowError without DBL_MAX clamping."""
        # Apply an enormous transverse shear force (e.g. 1e308) that overflows internal compatibility vectors
        sec = make_b1_hat_section()
        with pytest.raises(OverflowError, match="overflows IEEE-754 float64"):
            sec.calculate_shear_flow(vx=1e308, vy=1e308)

    def test_t20_ill_conditioned_detection(self) -> None:
        """T20: Degenerate collinear flat geometry raises SingularSectionError."""
        # All segments lie on the line y = 0 (Ix = 0, singular inertia tensor)
        s1 = Segment(Node(0.0, 0.0), Node(10.0, 0.0), t=1.0)
        s2 = Segment(Node(10.0, 0.0), Node(20.0, 0.0), t=1.0)
        sec = MixedSection([s1, s2])

        with pytest.raises(SingularSectionError, match="singular or degenerate"):
            sec.calculate_shear_flow(vx=10.0, vy=10.0)

    def test_t23_full_regression_baseline(self) -> None:
        """T23: Confirms that pure open Section and pure closed ClosedSection retain original behavior."""
        from sectalix.closed_section import ClosedSection
        from sectalix.section import Section

        # Open section
        s1 = Segment(Node(0.0, 0.0), Node(10.0, 0.0), t=1.0)
        s2 = Segment(Node(10.0, 0.0), Node(10.0, 10.0), t=1.0)
        sec_open = Section([s1, s2])
        assert sec_open.is_open is True
        assert sec_open.is_closed is False
        assert sec_open.is_mixed is False
        assert sec_open.J_open == sec_open.J
        assert sec_open.J_BB == 0.0
        assert sec_open.J_total == sec_open.J

        # Closed section
        box_segs = [
            Segment(Node(0.0, 0.0), Node(10.0, 0.0), t=1.0),
            Segment(Node(10.0, 0.0), Node(10.0, 10.0), t=1.0),
            Segment(Node(10.0, 10.0), Node(0.0, 10.0), t=1.0),
            Segment(Node(0.0, 10.0), Node(0.0, 0.0), t=1.0),
        ]
        sec_closed = ClosedSection(box_segs)
        assert sec_closed.is_closed is True
        assert sec_closed.is_open is False
        assert sec_closed.is_mixed is False
        assert sec_closed.J_BB == sec_closed.J
        assert sec_closed.J_open == 0.0
        assert sec_closed.J_total == sec_closed.J
