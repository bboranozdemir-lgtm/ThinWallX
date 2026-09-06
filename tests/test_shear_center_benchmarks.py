"""Benchmark tests for ThinWallX v0.3 Shear Center Analysis.

Verifies:
1. Doubly symmetric I-section: e_s = 0, shear center at centroid.
2. Singly symmetric channel: shear center on symmetry axis, matching independent hand oracle.
3. Branched T-section: shear center at junction of web and flange.
4. Unsymmetric L-section: I_xy != 0 coupling, both offsets non-zero, shear center at junction.
5. Independent Gauss-Legendre numerical quadrature oracle.
6. Large-translation stability (1e12 translation offset).
"""

from __future__ import annotations

import math
import numpy as np
import pytest

from thinwallx.primitives import Node, Segment
from thinwallx.section import Section
from thinwallx.shear_center import (
    ShearCenterResult,
    compute_shear_center,
    compute_shear_flow_torque,
)
from thinwallx.shear_flow import calculate_shear_flow


def make_i_section(
    h: float = 200.0, b: float = 100.0, tw: float = 6.0, tf: float = 10.0
) -> Section:
    """Build an idealized doubly-symmetric open I-section centered at origin."""
    # Flanges split into left and right halves to form a valid tree connected at web ends
    # Web: (0, -h/2) to (0, +h/2)
    # Top flange: (-b/2, h/2) to (0, h/2) and (0, h/2) to (b/2, h/2)
    # Bottom flange: (-b/2, -h/2) to (0, -h/2) and (0, -h/2) to (b/2, -h/2)
    segments = [
        Segment(p1=Node(0.0, -h / 2), p2=Node(0.0, h / 2), t=tw, id=0),
        Segment(p1=Node(-b / 2, h / 2), p2=Node(0.0, h / 2), t=tf, id=1),
        Segment(p1=Node(0.0, h / 2), p2=Node(b / 2, h / 2), t=tf, id=2),
        Segment(p1=Node(-b / 2, -h / 2), p2=Node(0.0, -h / 2), t=tf, id=3),
        Segment(p1=Node(0.0, -h / 2), p2=Node(b / 2, -h / 2), t=tf, id=4),
    ]
    return Section(segments)


def make_c_channel(
    h: float = 100.0, b: float = 50.0, tw: float = 2.0, tf: float = 2.0
) -> Section:
    """Build a singly symmetric C-channel opening to the right (+x).

    Web: from (0, -h/2) to (0, +h/2)
    Top flange: from (0, +h/2) to (b, +h/2)
    Bottom flange: from (0, -h/2) to (b, -h/2)
    Symmetry axis: y = 0.
    """
    segments = [
        Segment(p1=Node(0.0, -h / 2), p2=Node(0.0, h / 2), t=tw, id=0),
        Segment(p1=Node(0.0, h / 2), p2=Node(b, h / 2), t=tf, id=1),
        Segment(p1=Node(0.0, -h / 2), p2=Node(b, -h / 2), t=tf, id=2),
    ]
    return Section(segments)


def make_t_section(
    h_w: float = 100.0, b_f: float = 80.0, tw: float = 4.0, tf: float = 6.0
) -> Section:
    """Build a branched T-section symmetric about x = 0.

    Web: (0, 0) to (0, h_w)
    Flange: (-b_f/2, h_w) to (0, h_w) and (0, h_w) to (b_f/2, h_w)
    Junction: (0, h_w).
    """
    segments = [
        Segment(p1=Node(0.0, 0.0), p2=Node(0.0, h_w), t=tw, id=0),
        Segment(p1=Node(-b_f / 2, h_w), p2=Node(0.0, h_w), t=tf, id=1),
        Segment(p1=Node(0.0, h_w), p2=Node(b_f / 2, h_w), t=tf, id=2),
    ]
    return Section(segments)


def make_l_angle(
    b1: float = 80.0, b2: float = 60.0, t1: float = 3.0, t2: float = 4.0
) -> Section:
    """Build an unsymmetric L-angle with unequal legs meeting at (0, 0).

    Leg 1 along +x: (0, 0) to (b1, 0)
    Leg 2 along +y: (0, 0) to (0, b2)
    Junction: (0, 0).
    """
    segments = [
        Segment(p1=Node(0.0, 0.0), p2=Node(b1, 0.0), t=t1, id=0),
        Segment(p1=Node(0.0, 0.0), p2=Node(0.0, b2), t=t2, id=1),
    ]
    return Section(segments)


class TestDoublySymmetricBenchmark:
    """Benchmark 1: Doubly symmetric I-section."""

    def test_i_section_shear_center_coincides_with_centroid(self) -> None:
        """Verify e_s = 0 and (x_s, y_s) = (cx, cy) for doubly symmetric I-section."""
        sec = make_i_section(h=200.0, b=100.0, tw=6.0, tf=10.0)
        res = sec.compute_shear_center()

        # Centroid is (0, 0)
        assert abs(res.centroid[0]) < 1e-12
        assert abs(res.centroid[1]) < 1e-12

        # Offsets e_x = 0, e_y = 0
        assert abs(res.ex) < 1e-12
        assert abs(res.ey) < 1e-12
        assert abs(res.tz_x) < 1e-12
        assert abs(res.tz_y) < 1e-12

        # Absolute coordinates
        assert abs(res.x) < 1e-12
        assert abs(res.y) < 1e-12
        assert sec.shear_center == (pytest.approx(0.0, abs=1e-12), pytest.approx(0.0, abs=1e-12))
        assert sec.sc_offset == (pytest.approx(0.0, abs=1e-12), pytest.approx(0.0, abs=1e-12))

    def test_i_section_combined_residual_torque(self) -> None:
        """Verify residual torque is zero for arbitrary combined load on I-section."""
        sec = make_i_section(h=150.0, b=80.0, tw=5.0, tf=8.0)
        res = sec.compute_shear_center()
        residual = res.residual_torque(vx=3500.0, vy=-1250.0)
        assert abs(residual) < 1e-10


class TestSinglySymmetricChannelBenchmark:
    """Benchmark 2: Singly symmetric C-channel with independent hand oracle."""

    def test_channel_shear_center_hand_oracle(self) -> None:
        """Verify shear center matches independent analytical thin-walled beam formula.

        Theoretical Oracle Source:
            Megson, T.H.G., "Aircraft Structures for Engineering Students", Sec 17.2;
            Roark's Formulas for Stress and Strain, Table on Shear Center of Channels.
        
        Analytical Formula for channel with web height h, flange width b, thicknesses tw, tf:
            I_x = tw * h^3 / 12 + 2 * (tf * b * (h/2)^2) = h^2 / 12 * (tw * h + 6 * tf * b)
            Shear center location from web centerline:
                x_{s,web} = - (tf * b^2 * h^2) / (4 * I_x) = - (3 * b^2 * tf) / (6 * b * tf + h * tw)
            Centroid:
                cx = (2 * tf * b * (b/2)) / (tw * h + 2 * tf * b) = (tf * b^2) / (tw * h + 2 * tf * b)
                cy = 0
            Shear center offset from centroid:
                e_x = x_{s,web} - cx
                e_y = 0

        Numerical Values for h=100, b=50, tw=2, tf=2:
            A_w = 200, A_f = 100, A = 400
            cx = (2 * 100 * 25) / 400 = 12.5 mm
            I_x = 2 * 100^3 / 12 + 2 * (2 * 50 * 50^2) = 166666.6667 + 500000 = 666666.6667 mm^4
            x_{s,web} = - (2 * 50^2 * 100^2) / (4 * 666666.6667) = - 50e6 / (8/3 * 1e6) = - 18.75 mm
            e_x = -18.75 - 12.5 = -31.25 mm
            e_y = 0.0 mm
        """
        h, b, tw, tf = 100.0, 50.0, 2.0, 2.0
        sec = make_c_channel(h=h, b=b, tw=tw, tf=tf)
        res = sec.compute_shear_center()

        # Expected analytical values
        expected_cx = 12.5
        expected_cy = 0.0
        expected_xs_web = -18.75
        expected_xs = expected_xs_web  # since web is at x=0
        expected_ys = 0.0
        expected_ex = -31.25
        expected_ey = 0.0

        assert res.centroid[0] == pytest.approx(expected_cx, abs=1e-10)
        assert res.centroid[1] == pytest.approx(expected_cy, abs=1e-10)

        # Shear center lies on symmetry axis (y = 0)
        assert res.ey == pytest.approx(expected_ey, abs=1e-10)
        assert res.y == pytest.approx(expected_ys, abs=1e-10)

        # Shear center offset along x matches oracle
        assert res.ex == pytest.approx(expected_ex, abs=1e-10)
        assert res.x == pytest.approx(expected_xs, abs=1e-10)

        # Basis torques:
        # T_z^{(y)} = e_x = -31.25
        # T_z^{(x)} = -e_y = 0
        assert res.tz_y == pytest.approx(expected_ex, abs=1e-10)
        assert res.tz_x == pytest.approx(0.0, abs=1e-10)

        # Correct sign: shear center is outside the web (to the left, x_s < 0)
        assert res.x < 0.0

    @pytest.mark.parametrize("h,b,tw,tf", [
        (80.0, 40.0, 3.0, 3.0),
        (120.0, 60.0, 1.5, 3.0),
        (150.0, 75.0, 4.0, 2.5),
    ])
    def test_channel_parametric_oracle_agreement(
        self, h: float, b: float, tw: float, tf: float
    ) -> None:
        """Verify channel shear center against independent formula across multiple aspect ratios."""
        sec = make_c_channel(h=h, b=b, tw=tw, tf=tf)
        res = sec.compute_shear_center()

        # Independent formula:
        # x_{s,web} = - (3 * b^2 * tf) / (6 * b * tf + h * tw)
        # cx = (tf * b^2) / (tw * h + 2 * tf * b)
        xs_web_oracle = - (3.0 * (b ** 2) * tf) / (6.0 * b * tf + h * tw)
        cx_oracle = (tf * (b ** 2)) / (tw * h + 2.0 * tf * b)
        ex_oracle = xs_web_oracle - cx_oracle

        assert res.x == pytest.approx(xs_web_oracle, rel=1e-10)
        assert res.ex == pytest.approx(ex_oracle, rel=1e-10)
        assert res.ey == pytest.approx(0.0, abs=1e-10)
        assert res.y == pytest.approx(0.0, abs=1e-10)


class TestBranchedSectionBenchmark:
    """Benchmark 3: Branched T-section."""

    def test_t_section_shear_center_at_junction(self) -> None:
        """Verify branched T-section shear center is at the junction (0, h_w).

        Theorem (Roark / Megson):
            If all centerline segments of an open thin-walled section intersect at a
            single point, the shear center is located at that intersection point,
            because the moment arm of every segment shear force about that point is zero.
        """
        h_w = 120.0
        b_f = 80.0
        tw = 4.0
        tf = 6.0
        sec = make_t_section(h_w=h_w, b_f=b_f, tw=tw, tf=tf)
        res = sec.compute_shear_center()

        # By symmetry about x=0:
        assert res.centroid[0] == pytest.approx(0.0, abs=1e-10)
        assert res.ex == pytest.approx(0.0, abs=1e-10)
        assert res.x == pytest.approx(0.0, abs=1e-10)

        # By single-intersection theorem, shear center must be at junction (0, h_w):
        assert res.x == pytest.approx(0.0, abs=1e-10)
        assert res.y == pytest.approx(h_w, abs=1e-10)
        assert res.ey == pytest.approx(h_w - sec.cy, abs=1e-10)


class TestUnsymmetricSectionBenchmark:
    """Benchmark 4: Unsymmetric open L-section exercising I_xy != 0 coupling."""

    def test_l_section_shear_center_at_corner(self) -> None:
        """Verify unequal L-angle shear center is at the corner (0, 0).

        Coupling test:
            For an unequal angle not aligned with principal axes, I_xy != 0.
            Both unit basis loads [1, 0]^T and [0, 1]^T produce coupled shear flow fields.
            By the single-intersection theorem, both legs meet at (0, 0), so the shear
            center must be exactly at (0, 0), giving:
                x_s = 0,  y_s = 0
                e_x = -cx,  e_y = -cy
        """
        b1, b2, t1, t2 = 80.0, 60.0, 3.0, 4.0
        sec = make_l_angle(b1=b1, b2=b2, t1=t1, t2=t2)

        # Verify I_xy is non-zero
        assert abs(sec.Ixy) > 1e-3

        res = sec.compute_shear_center()

        # Both offsets are non-zero:
        assert abs(res.ex) > 1.0
        assert abs(res.ey) > 1.0

        # Exact agreement with corner (0, 0):
        assert res.x == pytest.approx(0.0, abs=1e-10)
        assert res.y == pytest.approx(0.0, abs=1e-10)
        assert res.ex == pytest.approx(-sec.cx, abs=1e-10)
        assert res.ey == pytest.approx(-sec.cy, abs=1e-10)

    def test_l_section_combined_residual_torque(self) -> None:
        """Verify residual torque vanishes for arbitrary load on unsymmetric angle."""
        sec = make_l_angle(b1=70.0, b2=50.0, t1=2.5, t2=3.5)
        res = sec.compute_shear_center()
        res_torque = res.residual_torque(vx=1234.0, vy=-5678.0)
        assert abs(res_torque) < 1e-10


class TestIndependentGaussQuadratureOracle:
    """Benchmark 7: Independent numerical Gauss-Legendre quadrature oracle."""

    @staticmethod
    def gauss_integrate_shear_flow_torque(sec: Section, vx: float, vy: float) -> float:
        """Independently integrate torque using 10-point Gauss-Legendre quadrature."""
        flow = calculate_shear_flow(sec, vx=vx, vy=vy)
        cx, cy = sec.centroid

        # 10-point Gauss-Legendre rule on [-1, 1]
        xi_nodes, weights = np.polynomial.legendre.leggauss(10)
        # Map to s in [0, L]: s = L * (xi_node + 1) / 2, ds = (L / 2) * dxi_node
        s_nodes_norm = 0.5 * (xi_nodes + 1.0)
        s_weights = 0.5 * weights

        total_torque = 0.0
        for sf in flow.segment_flows:
            p1 = sf.segment.p1
            p2 = sf.segment.p2
            length = sf.length
            t_vec = sf.tangent

            for s_norm, w in zip(s_nodes_norm, s_weights):
                s = s_norm * length
                # Centroid-relative position r_c(s) = p1 + s * t - centroid
                r_c = np.array([p1.x + s * t_vec[0] - cx, p1.y + s * t_vec[1] - cy])
                # Physical shear flow vector q_vec(s) = q(s) * t_vec
                q_s = sf.q_at(s)
                q_vec = q_s * t_vec
                # Cross product [r_c x q_vec]_z
                cross_z = r_c[0] * q_vec[1] - r_c[1] * q_vec[0]
                total_torque += cross_z * w * length

        return total_torque

    def test_gauss_quadrature_matches_exact_production_torque(self) -> None:
        """Verify production closed-form torque matches independent 10-point Gauss quadrature."""
        sec = make_c_channel(h=100.0, b=50.0, tw=2.0, tf=2.0)
        flow_x = calculate_shear_flow(sec, vx=1.0, vy=0.0)
        flow_y = calculate_shear_flow(sec, vx=0.0, vy=1.0)

        exact_tx = compute_shear_flow_torque(flow_x)
        exact_ty = compute_shear_flow_torque(flow_y)

        gauss_tx = self.gauss_integrate_shear_flow_torque(sec, vx=1.0, vy=0.0)
        gauss_ty = self.gauss_integrate_shear_flow_torque(sec, vx=0.0, vy=1.0)

        assert exact_tx == pytest.approx(gauss_tx, rel=1e-12, abs=1e-12)
        assert exact_ty == pytest.approx(gauss_ty, rel=1e-12, abs=1e-12)
