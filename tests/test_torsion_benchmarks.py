"""Benchmark tests for ThinWallX v0.4 Torsion and Warping Analysis.

Verifies:
1. Doubly symmetric I-section: analytical J and Cw = tf * b^3 * h^2 / 24, web omega = 0.
2. Singly symmetric channel: independent Gauss-Legendre quadrature oracle and textbook Cw.
3. Branched T-section: single-intersection theorem omega = 0, Cw = 0, exact J.
4. Unsymmetric L-angle: single-intersection theorem omega = 0, Cw = 0, exact J.
5. Uniform geometric scaling law: J ~ lambda^4, omega ~ lambda^2, Cw ~ lambda^6.
"""

from __future__ import annotations

import math
import sys
import numpy as np
import pytest

from thinwallx.primitives import Node, Segment
from thinwallx.section import Section
from thinwallx.torsion import (
    SegmentWarping,
    TorsionWarpingResult,
    compute_torsion_warping,
)


def make_i_section(
    h: float = 200.0, b: float = 100.0, tw: float = 6.0, tf: float = 8.0
) -> Section:
    """Idealized doubly-symmetric open I-section centered at origin."""
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
    """Singly symmetric C-channel opening to +x."""
    segments = [
        Segment(p1=Node(0.0, -h / 2), p2=Node(0.0, h / 2), t=tw, id=0),
        Segment(p1=Node(0.0, h / 2), p2=Node(b, h / 2), t=tf, id=1),
        Segment(p1=Node(0.0, -h / 2), p2=Node(b, -h / 2), t=tf, id=2),
    ]
    return Section(segments)


def make_t_section(
    h_w: float = 100.0, b_f: float = 80.0, tw: float = 4.0, tf: float = 6.0
) -> Section:
    """Branched T-section symmetric about x = 0."""
    segments = [
        Segment(p1=Node(0.0, 0.0), p2=Node(0.0, h_w), t=tw, id=0),
        Segment(p1=Node(-b_f / 2, h_w), p2=Node(0.0, h_w), t=tf, id=1),
        Segment(p1=Node(0.0, h_w), p2=Node(b_f / 2, h_w), t=tf, id=2),
    ]
    return Section(segments)


def make_l_angle(
    b1: float = 80.0, b2: float = 60.0, t1: float = 3.0, t2: float = 4.0
) -> Section:
    """Unsymmetric L-angle meeting at (0, 0)."""
    segments = [
        Segment(p1=Node(0.0, 0.0), p2=Node(b1, 0.0), t=t1, id=0),
        Segment(p1=Node(0.0, 0.0), p2=Node(0.0, b2), t=t2, id=1),
    ]
    return Section(segments)


class TestISectionBenchmark:
    """Benchmark 1: Doubly symmetric I-section analytical J and Cw."""

    def test_i_section_torsion_and_warping_exact(self) -> None:
        """Verify J and C_w match the exact analytical thin-walled I-section formulas.

        Analytical Formulas:
            J = (2 * b * tf^3 + h * tw^3) / 3
            Cw = (tf * b^3 * h^2) / 24
            omega(web) = 0
            omega(flange tips) = +/- b * h / 4
        """
        h, b, tw, tf = 200.0, 100.0, 6.0, 8.0
        sec = make_i_section(h=h, b=b, tw=tw, tf=tf)
        res = sec.torsion_properties()

        # Analytical J
        expected_J = (2.0 * b * (tf ** 3) + h * (tw ** 3)) / 3.0
        assert res.J == pytest.approx(expected_J, rel=1e-12)
        assert sec.J == pytest.approx(expected_J, rel=1e-12)

        # Analytical Cw
        expected_Cw = (tf * (b ** 3) * (h ** 2)) / 24.0
        assert res.Cw == pytest.approx(expected_Cw, rel=1e-12)
        assert sec.Cw == pytest.approx(expected_Cw, rel=1e-12)

        # Zero-mean integral
        assert abs(res.integral_omega_da) < 1e-10

        # Web sectorial coordinate must be zero everywhere
        # Segment 0 is the web
        web_warping = res[0]
        assert abs(web_warping.omega1) < 1e-12
        assert abs(web_warping.omega2) < 1e-12
        assert abs(web_warping.omega_at_xi(0.5)) < 1e-12
        assert abs(web_warping.cw_segment) < 1e-12

        # Flange tip values: magnitude is b * h / 4 = 100 * 200 / 4 = 5000.0
        expected_tip_mag = b * h / 4.0

        # Exact signed analytical sectorial coordinates at all 4 tips and junctions:
        # Segment 1: (-b/2, h/2) -> (0, h/2)  ==> omega1 is left-top (+bh/4), omega2 is top junction (0)
        assert res[1].omega1 == pytest.approx(+expected_tip_mag, abs=1e-10)
        assert abs(res[1].omega2) < 1e-10

        # Segment 2: (0, h/2) -> (b/2, h/2)   ==> omega1 is top junction (0), omega2 is right-top (-bh/4)
        assert abs(res[2].omega1) < 1e-10
        assert res[2].omega2 == pytest.approx(-expected_tip_mag, abs=1e-10)

        # Segment 3: (-b/2, -h/2) -> (0, -h/2) ==> omega1 is left-bottom (-bh/4), omega2 is bottom junction (0)
        assert res[3].omega1 == pytest.approx(-expected_tip_mag, abs=1e-10)
        assert abs(res[3].omega2) < 1e-10

        # Segment 4: (0, -h/2) -> (b/2, -h/2)  ==> omega1 is bottom junction (0), omega2 is right-bottom (+bh/4)
        assert abs(res[4].omega1) < 1e-10
        assert res[4].omega2 == pytest.approx(+expected_tip_mag, abs=1e-10)

        # Antisymmetric pattern across web centerline:
        # Top-left and top-right must have opposite signs
        assert res[1].omega1 == pytest.approx(-res[2].omega2, abs=1e-10)
        # Bottom-left and bottom-right must have opposite signs
        assert res[3].omega1 == pytest.approx(-res[4].omega2, abs=1e-10)


class TestBendingWarpingOrthogonality:
    """Verify bending-warping orthogonality int_A omega * xc dA = 0 and int_A omega * yc dA = 0."""

    @staticmethod
    def compute_sectorial_product_moments(sec: Section) -> tuple[float, float, float]:
        """Compute exact integrals int_A omega * xc dA, int_A omega * yc dA and scale."""
        res = sec.torsion_properties()
        cx, cy = sec.centroid

        i_wx_terms: list[float] = []
        i_wy_terms: list[float] = []
        for sw in res.segment_warpings:
            seg = sw.segment
            xc1 = seg.p1.x - cx
            yc1 = seg.p1.y - cy
            xc2 = seg.p2.x - cx
            yc2 = seg.p2.y - cy
            w1 = sw.omega1
            w2 = sw.omega2

            # Exact integral of product of two linear functions over [0, L]:
            # \int_0^L f(s) g(s) ds = (L / 6) * (2*f1*g1 + f1*g2 + f2*g1 + 2*f2*g2)
            factor = seg.t * seg.length / 6.0
            term_x = factor * (2.0 * w1 * xc1 + w1 * xc2 + w2 * xc1 + 2.0 * w2 * xc2)
            term_y = factor * (2.0 * w1 * yc1 + w1 * yc2 + w2 * yc1 + 2.0 * w2 * yc2)

            i_wx_terms.append(term_x)
            i_wy_terms.append(term_y)

        iwx = math.fsum(i_wx_terms)
        iwy = math.fsum(i_wy_terms)

        max_len = max((seg.length for seg in sec.segments), default=1.0)
        # Characteristic dimensional scale for I_wx, I_wy (dimension L^5): Area * L_max^3
        scale = sec.area * (max_len ** 3)
        return iwx, iwy, scale

    def test_i_section_orthogonality(self) -> None:
        """I-section bending-warping product moments are zero to machine precision."""
        sec = make_i_section(h=200.0, b=100.0, tw=6.0, tf=8.0)
        iwx, iwy, scale = self.compute_sectorial_product_moments(sec)
        assert abs(iwx) / scale < 1e-12
        assert abs(iwy) / scale < 1e-12

    def test_channel_orthogonality(self) -> None:
        """Channel section bending-warping product moments are zero to machine precision."""
        sec = make_c_channel(h=100.0, b=50.0, tw=2.0, tf=2.0)
        iwx, iwy, scale = self.compute_sectorial_product_moments(sec)
        assert abs(iwx) / scale < 1e-12
        assert abs(iwy) / scale < 1e-12

    def test_unit_scale_absolute_orthogonality(self) -> None:
        """Unit-scale sections satisfy absolute orthogonality < 1e-12 directly."""
        sec_unit = Section.from_tuples([
            ((0.0, -0.5), (0.0, 0.5), 0.02),
            ((0.0, 0.5), (0.5, 0.5), 0.02),
            ((0.0, -0.5), (0.5, -0.5), 0.02),
        ])
        iwx, iwy, _ = self.compute_sectorial_product_moments(sec_unit)
        assert abs(iwx) < 1e-12
        assert abs(iwy) < 1e-12

    def test_branched_and_unsymmetric_orthogonality(self) -> None:
        """Branched T-section and unsymmetric L-angle satisfy orthogonality."""
        sec_t = make_t_section(h_w=100.0, b_f=80.0, tw=4.0, tf=6.0)
        iwx_t, iwy_t, scale_t = self.compute_sectorial_product_moments(sec_t)
        assert abs(iwx_t) / scale_t < 1e-12
        assert abs(iwy_t) / scale_t < 1e-12

        sec_l = make_l_angle(b1=80.0, b2=60.0, t1=3.0, t2=4.0)
        iwx_l, iwy_l, scale_l = self.compute_sectorial_product_moments(sec_l)
        assert abs(iwx_l) / scale_l < 1e-12
        assert abs(iwy_l) / scale_l < 1e-12


class TestChannelBenchmark:
    """Benchmark 2: Singly symmetric channel section."""

    def test_channel_analytical_cw_and_j(self) -> None:
        """Verify C-channel J and Cw against independent analytical formula.

        Textbook Reference (Roark / Galambos / Murray):
            For a channel opening to the right with web height h, flange width b,
            and uniform thicknesses tw, tf:
                Cw = (h^2 * b^3 * tf / 12) * (2*h*tw + 3*b*tf) / (h*tw + 6*b*tf)
        """
        h, b, tw, tf = 100.0, 50.0, 2.0, 2.0
        sec = make_c_channel(h=h, b=b, tw=tw, tf=tf)
        res = sec.torsion_properties()

        # Expected J
        expected_J = (2.0 * b * (tf ** 3) + h * (tw ** 3)) / 3.0
        assert res.J == pytest.approx(expected_J, rel=1e-12)

        # Expected Cw
        factor = (2.0 * h * tw + 3.0 * b * tf) / (h * tw + 6.0 * b * tf)
        expected_Cw = (h ** 2 * b ** 3 * tf / 12.0) * factor
        assert res.Cw == pytest.approx(expected_Cw, rel=1e-11)

        # Zero-mean condition
        assert abs(res.integral_omega_da) < 1e-10

    def test_channel_independent_gauss_quadrature_oracle(self) -> None:
        """Verify Cw matches independent 20-point Gauss-Legendre quadrature."""
        h, b, tw, tf = 120.0, 45.0, 3.0, 2.5
        sec = make_c_channel(h=h, b=b, tw=tw, tf=tf)
        res = sec.torsion_properties()

        # Independent Gauss-Legendre integration of omega(s)^2 along segments
        xi_nodes, weights = np.polynomial.legendre.leggauss(20)
        s_norm = 0.5 * (xi_nodes + 1.0)
        w_scaled = 0.5 * weights

        gauss_cw = 0.0
        for sw in res.segment_warpings:
            seg_integral = 0.0
            for xi, w in zip(s_norm, w_scaled):
                val = sw.omega_at_xi(xi)
                seg_integral += (val ** 2) * w
            gauss_cw += sw.segment.t * sw.segment.length * seg_integral

        assert res.Cw == pytest.approx(gauss_cw, rel=1e-12)


class TestBranchedAndUnsymmetricBenchmarks:
    """Benchmark 3 & 4: Sections where all lines of action meet at one point."""

    def test_branched_t_section_has_zero_warping(self) -> None:
        """Theorem: Centerline segments intersecting at a point have Cw = 0."""
        sec = make_t_section(h_w=100.0, b_f=80.0, tw=4.0, tf=6.0)
        res = sec.torsion_properties()

        # Exact J
        expected_J = (80.0 * (6.0 ** 3) + 100.0 * (4.0 ** 3)) / 3.0
        assert res.J == pytest.approx(expected_J, rel=1e-12)

        # Cw must be zero
        assert res.Cw == pytest.approx(0.0, abs=1e-12)

        # Normalized omega must be zero everywhere
        for node_val in res.node_omega:
            assert abs(node_val) < 1e-12

    def test_unsymmetric_l_angle_has_zero_warping(self) -> None:
        """Theorem: L-angle legs intersect at corner, so Cw = 0."""
        sec = make_l_angle(b1=80.0, b2=60.0, t1=3.0, t2=4.0)
        res = sec.torsion_properties()

        # Exact J
        expected_J = (80.0 * (3.0 ** 3) + 60.0 * (4.0 ** 3)) / 3.0
        assert res.J == pytest.approx(expected_J, rel=1e-12)

        # Cw must be zero
        assert res.Cw == pytest.approx(0.0, abs=1e-12)

        # Normalized omega must be zero everywhere
        for node_val in res.node_omega:
            assert abs(node_val) < 1e-12


class TestUniformGeometricScaling:
    """Benchmark 7: Uniform geometric scaling law."""

    @pytest.mark.parametrize("lam", [0.5, 1.8, 2.5, 10.0])
    def test_geometric_scaling_law(self, lam: float) -> None:
        """Criterion 24: J' = lam^4 * J, omega' = lam^2 * omega, Cw' = lam^6 * Cw."""
        base = make_c_channel(h=100.0, b=50.0, tw=2.0, tf=2.0)
        res_base = base.torsion_properties()

        scaled_segments = [
            Segment(
                p1=Node(seg.p1.x * lam, seg.p1.y * lam),
                p2=Node(seg.p2.x * lam, seg.p2.y * lam),
                t=seg.t * lam,
                id=seg.id,
            )
            for seg in base.segments
        ]
        scaled = Section(scaled_segments)
        res_scaled = scaled.torsion_properties()

        # J scales as lam^4
        assert res_scaled.J == pytest.approx(res_base.J * (lam ** 4), rel=1e-11)

        # Cw scales as lam^6
        assert res_scaled.Cw == pytest.approx(res_base.Cw * (lam ** 6), rel=1e-11)

        # Sectorial coordinates scale as lam^2
        for w_base, w_scaled in zip(res_base.node_omega, res_scaled.node_omega):
            assert w_scaled == pytest.approx(w_base * (lam ** 2), rel=1e-11, abs=1e-11)


class TestExtremeScaleRobustness:
    """Verify numeric stability against intermediate underflow in J and overflow in Cw."""

    def test_extreme_j_intermediate_underflow_prevention(self) -> None:
        """L = 2^300, t = 2^-400 yields J = 2 * (2^-900) / 3 > 0 without underflowing to 0.0."""
        L = 2.0 ** 300
        t = 2.0 ** -400
        sec = Section([
            Segment(p1=Node(0.0, 0.0), p2=Node(L, 0.0), t=t, id=0),
            Segment(p1=Node(0.0, 0.0), p2=Node(0.0, L), t=t, id=1),
        ])
        res = sec.torsion_properties()

        # Two segments of length L, thickness t: J = 2 * ((L * t / 3) * t) * t
        expected_J = 2.0 * ((2.0 ** -900) / 3.0)
        assert res.J > 0.0
        assert math.isfinite(res.J)
        assert res.J == pytest.approx(expected_J, rel=1e-12, abs=0.0)

    def test_extreme_cw_intermediate_squaring_overflow_prevention(self) -> None:
        """L = 2^270, t = 2^-340 where omega ~ 2^538 does not trigger squaring overflow."""
        L = 2.0 ** 270
        t = 2.0 ** -340
        h, b, tw, tf = L, L, t, t

        sec = Section([
            Segment(p1=Node(0.0, -h / 2), p2=Node(0.0, h / 2), t=tw, id=0),
            Segment(p1=Node(-b / 2, h / 2), p2=Node(0.0, h / 2), t=tf, id=1),
            Segment(p1=Node(0.0, h / 2), p2=Node(b / 2, h / 2), t=tf, id=2),
            Segment(p1=Node(-b / 2, -h / 2), p2=Node(0.0, -h / 2), t=tf, id=3),
            Segment(p1=Node(0.0, -h / 2), p2=Node(b / 2, -h / 2), t=tf, id=4),
        ])
        res = sec.torsion_properties()

        # Analytical Cw = tf * b^3 * h^2 / 24 = 2^-340 * 2^810 * 2^540 / 24 = 2^1010 / 24
        expected_Cw = (2.0 ** 1010) / 24.0
        assert math.isfinite(res.Cw)
        assert res.Cw > 0.0
        assert res.Cw == pytest.approx(expected_Cw, rel=1e-11)

    def test_extreme_j_characteristic_scaling_subnormal_prevention(self) -> None:
        """L = 2^300, t = 2^-458 yields J = (2/3) * 2^-1074 ~= 5e-324 without flushing to 0.0."""
        L = 2.0 ** 300
        t = 2.0 ** -458
        sec = Section([
            Segment(p1=Node(0.0, 0.0), p2=Node(L, 0.0), t=t, id=0),
            Segment(p1=Node(0.0, 0.0), p2=Node(0.0, L), t=t, id=1),
        ])
        res = sec.torsion_properties()

        expected_J = 2.0 * (2.0 ** -1074) / 3.0
        assert res.J > 0.0
        assert math.isfinite(res.J)
        assert res.J == pytest.approx(expected_J, rel=1e-12, abs=0.0)

    def test_extreme_cw_sum_of_squares_opposite_signs_adversarial(self) -> None:
        """3-segment chain where one endpoint has a^2 >= 1.05 * DBL_MAX and naive squaring yields inf.

        Verifies that:
        1. On segment 0, a and b have opposite signs and a >= sqrt(1.05 * DBL_MAX) ~= 1.37e154.
        2. Naive polynomial evaluation (a*a + a*b + b*b) overflows to inf.
        3. ThinWallX complete squares with characteristic scaling evaluates Cw to a finite, positive value < DBL_MAX.
        """
        angle1 = 135
        angle2 = -90
        rad1 = math.radians(angle1)
        rad2 = math.radians(angle2)
        L0, L1, L2 = 1.0, 1.0, 2.0
        p0 = (-L0 * math.cos(rad1), -L0 * math.sin(rad1))
        p1 = (0.0, 0.0)
        p2 = (L1, 0.0)
        p3 = (L1 + L2 * math.cos(rad2), L2 * math.sin(rad2))

        sec_unit = Section([
            Segment(p1=Node(*p0), p2=Node(*p1), t=1.0, id=0),
            Segment(p1=Node(*p1), p2=Node(*p2), t=1.0, id=1),
            Segment(p1=Node(*p2), p2=Node(*p3), t=1.0, id=2),
        ])
        res_unit = sec_unit.torsion_properties()
        sw0_unit = res_unit.segment_warpings[0]
        a0_unit = sw0_unit.omega1 * math.sqrt(sw0_unit.segment.t * sw0_unit.segment.length / 3.0)

        # Scale factor lambda so that a^2 >= 1.05 * DBL_MAX (a >= 1.37e154)
        target_a = math.sqrt(1.051) * math.sqrt(sys.float_info.max)
        lam = (target_a / a0_unit) ** (1.0 / 3.0)

        sec = Section([
            Segment(p1=Node(p0[0] * lam, p0[1] * lam), p2=Node(p1[0] * lam, p1[1] * lam), t=1.0 * lam, id=0),
            Segment(p1=Node(p1[0] * lam, p1[1] * lam), p2=Node(p2[0] * lam, p2[1] * lam), t=1.0 * lam, id=1),
            Segment(p1=Node(p2[0] * lam, p2[1] * lam), p2=Node(p3[0] * lam, p3[1] * lam), t=1.0 * lam, id=2),
        ])
        res = sec.torsion_properties()

        # Check segment 0 endpoint values
        sw0 = res.segment_warpings[0]
        scale_0 = math.sqrt(sw0.segment.t * sw0.segment.length / 3.0)
        a = sw0.omega1 * scale_0
        b = sw0.omega2 * scale_0

        # 1. Endpoints have opposite signs
        assert a * b < 0.0
        # 2. a^2 >= 1.05 * DBL_MAX (a >= 1.37e154)
        assert a >= math.sqrt(1.05) * math.sqrt(sys.float_info.max)
        # 3. Naive squaring overflows to inf in float64
        assert not math.isfinite(a * a)
        assert not math.isfinite((a * a) + (a * b) + (b * b))

        # 4. ThinWallX result is strictly finite, positive, and bounded by DBL_MAX
        assert math.isfinite(res.Cw)
        assert 0.0 < res.Cw < sys.float_info.max
        assert np.isclose(res.Cw, 1.7781306565303429e308, rtol=1e-14, atol=0.0)

    def test_extreme_cw_subnormal_prevention(self) -> None:
        """Cw = 5e-324 subnormal regression test."""
        angle1, angle2 = 135, -90
        rad1, rad2 = math.radians(angle1), math.radians(angle2)
        L0, L1, L2 = 1.0, 1.0, 2.0
        p0 = (-L0 * math.cos(rad1), -L0 * math.sin(rad1))
        p1 = (0.0, 0.0)
        p2 = (L1, 0.0)
        p3 = (L1 + L2 * math.cos(rad2), L2 * math.sin(rad2))

        sec_unit = Section([
            Segment(p1=Node(*p0), p2=Node(*p1), t=1.0, id=0),
            Segment(p1=Node(*p1), p2=Node(*p2), t=1.0, id=1),
            Segment(p1=Node(*p2), p2=Node(*p3), t=1.0, id=2),
        ])
        res_unit = sec_unit.torsion_properties()
        target_cw = 2.0 * (2.0 ** -1074) / 3.0  # subnormal 5e-324
        lam = (target_cw / res_unit.Cw) ** (1.0 / 6.0)

        sec_sub = Section([
            Segment(p1=Node(p0[0] * lam, p0[1] * lam), p2=Node(p1[0] * lam, p1[1] * lam), t=1.0 * lam, id=0),
            Segment(p1=Node(p1[0] * lam, p1[1] * lam), p2=Node(p2[0] * lam, p2[1] * lam), t=1.0 * lam, id=1),
            Segment(p1=Node(p2[0] * lam, p2[1] * lam), p2=Node(p3[0] * lam, p3[1] * lam), t=1.0 * lam, id=2),
        ], node_tolerance=1e-60)
        res_sub = sec_sub.torsion_properties()

        assert res_sub.Cw > 0.0
        assert math.isfinite(res_sub.Cw)
        assert res_sub.Cw == pytest.approx(target_cw, rel=1e-12, abs=0.0)

    def test_heterogeneous_multiscale_cw_four_branch_tree(self) -> None:
        """Four-branch H=2^300, t=1 and tips l=2^-200, ts=2^-700 open tree: Expected Cw = 2.534788755060213e-211."""
        H = 2.0 ** 300
        t_main = 1.0
        l_tip = 2.0 ** -200
        ts_tip = 2.0 ** -700

        segments = [
            Segment(p1=Node(0.0, 0.0), p2=Node(H, 0.0), t=t_main, id=0),
            Segment(p1=Node(H, 0.0), p2=Node(H, l_tip), t=ts_tip, id=1),
            Segment(p1=Node(0.0, 0.0), p2=Node(-H, 0.0), t=t_main, id=2),
            Segment(p1=Node(-H, 0.0), p2=Node(-H, -l_tip), t=ts_tip, id=3),
            Segment(p1=Node(0.0, 0.0), p2=Node(0.0, H), t=t_main, id=4),
            Segment(p1=Node(0.0, H), p2=Node(-l_tip, H), t=ts_tip, id=5),
            Segment(p1=Node(0.0, 0.0), p2=Node(0.0, -H), t=t_main, id=6),
            Segment(p1=Node(0.0, -H), p2=Node(l_tip, -H), t=ts_tip, id=7),
        ]
        sec = Section(segments, node_tolerance=2.0 ** -250)
        res = sec.torsion_properties()

        expected_cw = 2.534788755060213e-211
        assert math.isfinite(res.Cw)
        assert res.Cw > 0.0
        assert np.isclose(res.Cw, expected_cw, rtol=1e-14, atol=0.0)

    def test_heterogeneous_multiscale_ji_long_and_short_branches(self) -> None:
        """Two orthogonal branches L=2^300, t=2^-400 and third branch L=sqrt(2), t=2^-40: Long branches Ji = 3.943507287222582e-272."""
        L = 2.0 ** 300
        t_long = 2.0 ** -400
        t_short = 2.0 ** -40

        sec = Section([
            Segment(p1=Node(0.0, 0.0), p2=Node(L, 0.0), t=t_long, id=0),
            Segment(p1=Node(0.0, 0.0), p2=Node(0.0, L), t=t_long, id=1),
            Segment(p1=Node(0.0, 0.0), p2=Node(1.0, 1.0), t=t_short, id=2),
        ])
        res = sec.torsion_properties()

        expected_Ji = 3.943507287222582e-272
        sw0 = res.segment_warpings[0]
        sw1 = res.segment_warpings[1]

        assert math.isfinite(sw0.j_segment)
        assert sw0.j_segment > 0.0
        assert np.isclose(sw0.j_segment, expected_Ji, rtol=1e-14, atol=0.0)

        assert math.isfinite(sw1.j_segment)
        assert sw1.j_segment > 0.0
        assert np.isclose(sw1.j_segment, expected_Ji, rtol=1e-14, atol=0.0)

        assert math.isfinite(res.J)
        assert res.J > 0.0

