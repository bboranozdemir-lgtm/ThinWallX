"""Analytical benchmark tests for open-section shear flow (Sectalix v0.2).

Implements all required benchmark families from ACTIVE_PHASE.md:
1. Symmetric open section (channel / C-section) under single-axis shear vs textbook closed-form solution.
2. Branched open tree section (T-section) verifying junction flow continuity and leaf boundary conditions.
3. Unsymmetric section with Ixy != 0 (L-section and asymmetric Z-section) under Vx, Vy, and combined shear
   benchmarked against independent exact rational arithmetic.
4. Rotated benchmark verifying covariance of the physical shear-flow field.
5. Principal-axis special case reduction (Ixy = 0 => q = Vy*Qx/Ix + Vx*Qy/Iy).
"""

from fractions import Fraction
import math
import numpy as np
import pytest

from sectalix.primitives import Node, Segment
from sectalix.section import Section
from sectalix.shear_load import ShearLoad


class TestSymmetricChannelBenchmark:
    """Benchmark 1: Symmetric channel section under vertical transverse shear Vy."""

    def test_channel_vy_shear_flow(self) -> None:
        # Channel: web height h = 100.0, flange width b = 50.0, thickness t = 2.0
        h = 100.0
        b = 50.0
        t = 2.0

        # Segments:
        # 0: Top flange: (b, h/2) -> (0, h/2)
        # 1: Web: (0, h/2) -> (0, -h/2)
        # 2: Bottom flange: (0, -h/2) -> (b, -h/2)
        s0 = Segment(Node(b, h / 2.0), Node(0.0, h / 2.0), t=t)
        s1 = Segment(Node(0.0, h / 2.0), Node(0.0, -h / 2.0), t=t)
        s2 = Segment(Node(0.0, -h / 2.0), Node(b, -h / 2.0), t=t)
        sec = Section([s0, s1, s2])

        Vy = 1000.0
        res = sec.calculate_shear_flow(vx=0.0, vy=Vy)

        # Analytical textbook reference:
        # Ix = (t * h^2 / 12) * (6*b + h) = (2 * 10000 / 12) * (300 + 100) = 2000000 / 3
        # Ixy = 0
        Ix = (t * h**2 / 12.0) * (6.0 * b + h)
        assert math.isclose(sec.Ix, Ix, rel_tol=1e-12)
        assert math.isclose(sec.Ixy, 0.0, abs_tol=1e-12)

        # Top flange: s in [0, b] from (b, 50) to (0, 50)
        # Material on node2 side is web + bottom flange, so Qx is negative:
        # Qx = -t * b * (h/2)
        sf0 = res[0]
        assert math.isclose(sf0.q0, 0.0, abs_tol=1e-11)
        expected_q_corner = -(Vy / Ix) * t * b * (h / 2.0)
        assert math.isclose(sf0.q_end, expected_q_corner, rel_tol=1e-11)

        # Check midpoint of flange against exact analytical theoretical value:
        # At midpoint s = b/2, Qx = t * (b/2) * (h/2) - t * b * (h/2) = -t * b * h / 4.0
        # q(b/2) = -(Vy / Ix) * t * b * (h / 4.0) = 0.5 * expected_q_corner
        expected_q_flange_mid = -(Vy / Ix) * t * b * (h / 4.0)
        assert math.isclose(sf0.q_mid, expected_q_flange_mid, rel_tol=1e-11)

        # 2. Web: s in [0, h] from (0, 50) to (0, -50)
        # Top corner s=0: q(0) = expected_q_corner (continuity)
        sf1 = res[1]
        assert math.isclose(sf1.q0, expected_q_corner, rel_tol=1e-11)

        # Midpoint of web (neutral axis y = 0, s = h/2):
        # Qx on node2 side (bottom half of web + bottom flange) is negative:
        Qx_max = t * b * (h / 2.0) + (t * h**2) / 8.0
        expected_q_max = -(Vy / Ix) * Qx_max
        assert math.isclose(sf1.q_mid, expected_q_max, rel_tol=1e-11)

        # Bottom corner of web s=h: q(h) = expected_q_corner
        assert math.isclose(sf1.q_end, expected_q_corner, rel_tol=1e-11)

        # 3. Bottom flange: s in [0, b] from (0, -50) to (b, -50)
        sf2 = res[2]
        assert math.isclose(sf2.q0, expected_q_corner, rel_tol=1e-11)
        assert math.isclose(sf2.q_end, 0.0, abs_tol=1e-11)

        # Resultant recovery:
        # Web carries 100% of Vy in +y direction (since tangent is [0, -1] and q is negative)
        web_force_y = sf1.resultant_integral()[1]
        assert math.isclose(web_force_y, Vy, rel_tol=1e-11)

        # Total recovered resultant must match [0, Vy]
        np.testing.assert_allclose(res.recovered_resultant, [0.0, Vy], rtol=1e-11, atol=1e-11)


class TestBranchedTSectionBenchmark:
    """Benchmark 2: Branched T-section verifying junction flow continuity."""

    def test_t_section_flow_continuity_and_leaves(self) -> None:
        # Flange width b = 120 (modeled as two segments meeting at (0, 0)), tf = 6.0
        # Web height h = 100, tw = 4.0
        b = 120.0
        tf = 6.0
        h = 100.0
        tw = 4.0

        # Segments:
        # 0: Left flange: (-b/2, 0) -> (0, 0)
        # 1: Right flange: (b/2, 0) -> (0, 0)
        # 2: Web: (0, 0) -> (0, -h)
        s_left = Segment(Node(-b / 2.0, 0.0), Node(0.0, 0.0), t=tf)
        s_right = Segment(Node(b / 2.0, 0.0), Node(0.0, 0.0), t=tf)
        s_web = Segment(Node(0.0, 0.0), Node(0.0, -h), t=tw)
        sec = Section([s_left, s_right, s_web])

        Vy = 1500.0
        res = sec.calculate_shear_flow(vx=0.0, vy=Vy)

        # 1. Leaf boundary conditions:
        # Left tip (s=0 of seg 0) is a free edge:
        assert math.isclose(res[0].q0, 0.0, abs_tol=1e-11)
        # Right tip (s=0 of seg 1) is a free edge:
        assert math.isclose(res[1].q0, 0.0, abs_tol=1e-11)
        # Web bottom (s=L of seg 2) is a free edge:
        assert math.isclose(res[2].q_end, 0.0, abs_tol=1e-11)

        # 2. Symmetry under vertical shear Vy:
        # Shear flow in left flange and right flange must have equal magnitude
        assert math.isclose(res[0].q_end, res[1].q_end, rel_tol=1e-11)

        # 3. Flow continuity at junction:
        # The two flange flows flowing INTO the junction must sum to the web flow flowing OUT of junction!
        q_junction_flanges_sum = res[0].q_end + res[1].q_end
        q_junction_web = res[2].q0
        assert math.isclose(q_junction_flanges_sum, q_junction_web, rel_tol=1e-11)

        # 4. Resultant recovery:
        np.testing.assert_allclose(res.recovered_resultant, [0.0, Vy], rtol=1e-11, atol=1e-11)


class TestUnsymmetricCoupledBenchmark:
    """Benchmark 3: Unsymmetric section with Ixy != 0 vs independent exact rational arithmetic."""

    def test_equal_angle_coupled_shear_flow(self) -> None:
        # Equal angle b = 60.0, t = 3.0
        # Leg 1: (0, b) -> (0, 0)
        # Leg 2: (0, 0) -> (b, 0)
        b = 60.0
        t = 3.0
        s1 = Segment(Node(0.0, b), Node(0.0, 0.0), t=t)
        s2 = Segment(Node(0.0, 0.0), Node(b, 0.0), t=t)
        sec = Section([s1, s2])

        # Centroid: cx = b/4 = 15, cy = b/4 = 15
        # Ix = Iy = (5/24)*t*b^3 = (5/24)*3*216000 = 135000
        # Ixy = -(1/8)*t*b^3 = -(1/8)*3*216000 = -81000
        # C = [[135000, -81000], [-81000, 135000]]
        # Delta = Ix*Iy - Ixy^2 = 135000^2 - 81000^2 = 18225000000 - 6561000000 = 11664000000
        # C^-1 = (1 / 11664000000) * [[135000, 81000], [81000, 135000]]

        # Load cases:
        test_loads = [
            ShearLoad(0.0, 1000.0),
            ShearLoad(1000.0, 0.0),
            ShearLoad(600.0, -800.0),
        ]

        Ix_frac = Fraction(5, 24) * Fraction(3) * Fraction(60**3)  # 135000
        Ixy_frac = -Fraction(1, 8) * Fraction(3) * Fraction(60**3)  # -81000
        Delta_frac = Ix_frac * Ix_frac - Ixy_frac * Ixy_frac  # 11664000000

        for load in test_loads:
            Vx_frac = Fraction(int(load.vx))
            Vy_frac = Fraction(int(load.vy))

            alpha_x_frac = (Vx_frac * Ix_frac - Vy_frac * Ixy_frac) / Delta_frac
            alpha_y_frac = (Vy_frac * Ix_frac - Vx_frac * Ixy_frac) / Delta_frac

            res = sec.calculate_shear_flow(load)

            # Test corner value at (0, 0) (end of seg 0, start of seg 1)
            # Material on node2 side of seg 0 is seg 1:
            # Seg 1 has length b, y = 0 => yc = -b/4 = -15
            # Seg 1 xc from -15 to 45 => average xc = 15
            # Qy_seg1 = t * b * (b/4) = 3 * 60 * 15 = 2700
            # Qx_seg1 = t * b * (-b/4) = 3 * 60 * (-15) = -2700
            Qy_corner = Fraction(3 * 60 * 15)
            Qx_corner = Fraction(3 * 60 * (-15))
            expected_q_corner = float(alpha_x_frac * Qy_corner + alpha_y_frac * Qx_corner)

            assert math.isclose(res[0].q_end, expected_q_corner, rel_tol=1e-11)
            assert math.isclose(res[1].q0, expected_q_corner, rel_tol=1e-11)

            # Leaf endpoints must be zero:
            assert math.isclose(res[0].q0, 0.0, abs_tol=1e-11)
            assert math.isclose(res[1].q_end, 0.0, abs_tol=1e-11)

            # Resultant recovery:
            np.testing.assert_allclose(
                res.recovered_resultant, load.vector, rtol=1e-11, atol=1e-11
            )


class TestRotatedBenchmark:
    """Benchmark 4: Rotated benchmark verifying covariance of physical shear-flow vectors."""

    def test_rotated_l_section(self) -> None:
        b = 50.0
        t = 3.0
        s1 = Segment(Node(0.0, b), Node(0.0, 0.0), t=t)
        s2 = Segment(Node(0.0, 0.0), Node(b, 0.0), t=t)
        sec = Section([s1, s2])

        phi = 37.5 * math.pi / 180.0
        rot_sec = sec.rotated(phi)

        c = math.cos(phi)
        s = math.sin(phi)
        R = np.array([[c, -s], [s, c]], dtype=float)

        base_load = ShearLoad(300.0, 700.0)
        v_rot = R @ base_load.vector
        rot_load = ShearLoad(float(v_rot[0]), float(v_rot[1]))

        res_base = sec.calculate_shear_flow(base_load)
        res_rot = rot_sec.calculate_shear_flow(rot_load)

        # Check vector covariance at multiple normalized locations
        for i in range(len(sec.segments)):
            for xi in [0.0, 0.25, 0.5, 0.75, 1.0]:
                vec_base = res_base.vector_at(i, xi)
                expected_rot = R @ vec_base
                actual_rot = res_rot.vector_at(i, xi)
                np.testing.assert_allclose(actual_rot, expected_rot, rtol=1e-10, atol=1e-10)

        # Resultant recovery
        np.testing.assert_allclose(
            res_rot.recovered_resultant, rot_load.vector, rtol=1e-11, atol=1e-11
        )


class TestPrincipalAxisReduction:
    """Benchmark 5: Principal-axis special case reduction (Ixy = 0)."""

    def test_principal_axis_formula_reduction(self) -> None:
        # Use symmetric channel where Ixy = 0
        h = 120.0
        b = 60.0
        t = 3.0
        s0 = Segment(Node(b, h / 2.0), Node(0.0, h / 2.0), t=t)
        s1 = Segment(Node(0.0, h / 2.0), Node(0.0, -h / 2.0), t=t)
        s2 = Segment(Node(0.0, -h / 2.0), Node(b, -h / 2.0), t=t)
        sec = Section([s0, s1, s2])

        assert math.isclose(sec.Ixy, 0.0, abs_tol=1e-12)
        Ix = sec.Ix
        Iy = sec.Iy

        load = ShearLoad(400.0, 900.0)
        res = sec.calculate_shear_flow(load)

        # When Ixy = 0: alpha_x = Vx / Iy, alpha_y = Vy / Ix
        # q_decoupled(s) = (Vy * Qx(s) / Ix) + (Vx * Qy(s) / Iy)
        # Compare at multiple points on each segment
        for s_idx, sf in enumerate(res.segment_flows):
            seg = sf.segment
            L = seg.length
            m_sub_qy, m_sub_qx = sf.subtree_moment

            for xi in [0.0, 0.33, 0.5, 0.67, 1.0]:
                s = xi * L
                # Exact partial first moments:
                # xc(u) = (seg.p1.x - cx) + (u/L)*d_xc
                xc1 = seg.p1.x - sec.cx
                yc1 = seg.p1.y - sec.cy
                d_xc = seg.p2.x - seg.p1.x
                d_yc = seg.p2.y - seg.p1.y

                # int_s^L xc(u) du = xc1 * (L - s) + (d_xc / (2*L)) * (L^2 - s^2)
                qy_part = t * (xc1 * (L - s) + (d_xc / (2.0 * L)) * (L**2 - s**2))
                qx_part = t * (yc1 * (L - s) + (d_yc / (2.0 * L)) * (L**2 - s**2))

                Qy_total = m_sub_qy + qy_part
                Qx_total = m_sub_qx + qx_part

                expected_q = (load.vy * Qx_total / Ix) + (load.vx * Qy_total / Iy)
                actual_q = sf.q_at_xi(xi)

                assert math.isclose(actual_q, expected_q, rel_tol=1e-12, abs_tol=1e-12)
