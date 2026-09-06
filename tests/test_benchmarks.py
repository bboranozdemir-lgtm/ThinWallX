"""Analytical benchmark tests for thin-walled open sections.

Implements all 5 required benchmark families specified in ACTIVE_PHASE.md:
1. Single straight centerline wall (vertical & inclined).
2. L-shaped section (equal and unequal leg angle).
3. Channel / C-shaped section.
4. T-shaped section (branched-node topology).
5. Asymmetric multi-segment open section.
Plus the isotropic degenerate section case.

Every test validates against independent closed-form analytical derivations.
"""

from fractions import Fraction
import math
import numpy as np
import pytest

from sectalix.primitives import Node, Segment
from sectalix.section import Section


class TestSingleWallBenchmark:
    """Benchmark Family 1: Single straight centerline wall."""

    def test_vertical_wall(self) -> None:
        # Height H = 120.0, thickness t = 3.0, centered at origin
        H = 120.0
        t = 3.0
        p1 = Node(0.0, -H / 2.0)
        p2 = Node(0.0, H / 2.0)
        sec = Section([Segment(p1, p2, t=t)])

        # Analytical reference values (thin-wall line integration model)
        expected_area = t * H
        expected_cx = 0.0
        expected_cy = 0.0
        expected_ix = (t * H**3) / 12.0
        expected_iy = 0.0
        expected_ixy = 0.0
        expected_i1 = expected_ix
        expected_i2 = 0.0
        expected_theta_p = 0.0

        assert math.isclose(sec.area, expected_area, rel_tol=1e-12)
        assert math.isclose(sec.cx, expected_cx, abs_tol=1e-12)
        assert math.isclose(sec.cy, expected_cy, abs_tol=1e-12)
        assert math.isclose(sec.Ix, expected_ix, rel_tol=1e-12)
        assert math.isclose(sec.Iy, expected_iy, abs_tol=1e-12)
        assert math.isclose(sec.Ixy, expected_ixy, abs_tol=1e-12)
        assert math.isclose(sec.I1, expected_i1, rel_tol=1e-12)
        assert math.isclose(sec.I2, expected_i2, abs_tol=1e-12)
        assert math.isclose(sec.theta_p, expected_theta_p, abs_tol=1e-12)

    def test_inclined_wall_30_deg(self) -> None:
        # Length L = 100.0, thickness t = 2.0, inclined at alpha = 30 deg (pi/6)
        L = 100.0
        t = 2.0
        alpha = math.pi / 6.0
        c = math.cos(alpha)
        s = math.sin(alpha)

        p1 = Node(-0.5 * L * c, -0.5 * L * s)
        p2 = Node(0.5 * L * c, 0.5 * L * s)
        sec = Section([Segment(p1, p2, t=t)])

        expected_area = t * L
        # Integral s^2 ds from -L/2 to L/2 is L^3 / 12
        I_line = (t * L**3) / 12.0
        expected_ix = I_line * s**2
        expected_iy = I_line * c**2
        expected_ixy = I_line * s * c
        expected_i1 = I_line
        expected_i2 = 0.0

        assert math.isclose(sec.area, expected_area, rel_tol=1e-12)
        assert math.isclose(sec.cx, 0.0, abs_tol=1e-12)
        assert math.isclose(sec.cy, 0.0, abs_tol=1e-12)
        assert math.isclose(sec.Ix, expected_ix, rel_tol=1e-12)
        assert math.isclose(sec.Iy, expected_iy, rel_tol=1e-12)
        assert math.isclose(sec.Ixy, expected_ixy, rel_tol=1e-12)
        assert math.isclose(sec.I1, expected_i1, rel_tol=1e-12)
        assert math.isclose(sec.I2, expected_i2, abs_tol=1e-12)


class TestLSectionBenchmark:
    """Benchmark Family 2: L-shaped thin-walled section (equal and unequal leg angle)."""

    def test_equal_leg_angle(self) -> None:
        # Equal angle b = 80.0, t = 4.0
        # Leg 1: (0, b) -> (0, 0)
        # Leg 2: (0, 0) -> (b, 0)
        b = 80.0
        t = 4.0

        s1 = Segment(Node(0.0, b), Node(0.0, 0.0), t=t)
        s2 = Segment(Node(0.0, 0.0), Node(b, 0.0), t=t)
        sec = Section([s1, s2])

        # Analytical derivations:
        # A = 2*b*t
        # cx = b/4, cy = b/4
        # Ix = (5/24)*t*b^3
        # Iy = (5/24)*t*b^3
        # Ixy = -(1/8)*t*b^3 = -(3/24)*t*b^3
        # I1 = (1/3)*t*b^3
        # I2 = (1/12)*t*b^3
        # theta_p = 45 deg (pi/4)
        expected_area = 2.0 * b * t
        expected_c = b / 4.0
        expected_ix = (5.0 / 24.0) * t * (b**3)
        expected_iy = (5.0 / 24.0) * t * (b**3)
        expected_ixy = -(1.0 / 8.0) * t * (b**3)
        expected_i1 = (1.0 / 3.0) * t * (b**3)
        expected_i2 = (1.0 / 12.0) * t * (b**3)
        expected_theta_p = math.pi / 4.0

        assert math.isclose(sec.area, expected_area, rel_tol=1e-12)
        assert math.isclose(sec.cx, expected_c, rel_tol=1e-12)
        assert math.isclose(sec.cy, expected_c, rel_tol=1e-12)
        assert math.isclose(sec.Ix, expected_ix, rel_tol=1e-12)
        assert math.isclose(sec.Iy, expected_iy, rel_tol=1e-12)
        assert math.isclose(sec.Ixy, expected_ixy, rel_tol=1e-12)
        assert math.isclose(sec.I1, expected_i1, rel_tol=1e-12)
        assert math.isclose(sec.I2, expected_i2, rel_tol=1e-12)
        assert math.isclose(sec.theta_p, expected_theta_p, rel_tol=1e-12)
        assert math.isclose(sec.theta_p_deg, 45.0, rel_tol=1e-12)

    def test_unequal_leg_angle(self) -> None:
        # Unequal angle: vertical leg b2 = 100, horizontal leg b1 = 50, t = 5.0
        b1 = 50.0
        b2 = 100.0
        t = 5.0

        s1 = Segment(Node(0.0, b2), Node(0.0, 0.0), t=t)
        s2 = Segment(Node(0.0, 0.0), Node(b1, 0.0), t=t)
        sec = Section([s1, s2])

        expected_area = t * (b1 + b2)
        expected_cx = (b1**2) / (2.0 * (b1 + b2))
        expected_cy = (b2**2) / (2.0 * (b1 + b2))
        expected_ix = (1.0 / 3.0) * t * (b2**3) - expected_area * (expected_cy**2)
        expected_iy = (1.0 / 3.0) * t * (b1**3) - expected_area * (expected_cx**2)
        expected_ixy = -expected_area * expected_cx * expected_cy

        diff = expected_ix - expected_iy
        r = math.hypot(diff / 2.0, expected_ixy)
        expected_i1 = (expected_ix + expected_iy) / 2.0 + r
        expected_i2 = (expected_ix + expected_iy) / 2.0 - r
        expected_theta_p = 0.5 * math.atan2(-2.0 * expected_ixy, diff)

        assert math.isclose(sec.area, expected_area, rel_tol=1e-12)
        assert math.isclose(sec.cx, expected_cx, rel_tol=1e-12)
        assert math.isclose(sec.cy, expected_cy, rel_tol=1e-12)
        assert math.isclose(sec.Ix, expected_ix, rel_tol=1e-12)
        assert math.isclose(sec.Iy, expected_iy, rel_tol=1e-12)
        assert math.isclose(sec.Ixy, expected_ixy, rel_tol=1e-12)
        assert math.isclose(sec.I1, expected_i1, rel_tol=1e-12)
        assert math.isclose(sec.I2, expected_i2, rel_tol=1e-12)
        assert math.isclose(sec.theta_p, expected_theta_p, rel_tol=1e-12)


class TestChannelSectionBenchmark:
    """Benchmark Family 3: Channel / C-shaped thin-walled section."""

    def test_channel_section(self) -> None:
        # Channel with web height h = 150, flange width b = 75, thickness t = 4.5
        h = 150.0
        b = 75.0
        t = 4.5

        # Flange 1: (b, h/2) -> (0, h/2)
        # Web: (0, h/2) -> (0, -h/2)
        # Flange 2: (0, -h/2) -> (b, -h/2)
        s1 = Segment(Node(b, h / 2.0), Node(0.0, h / 2.0), t=t)
        s2 = Segment(Node(0.0, h / 2.0), Node(0.0, -h / 2.0), t=t)
        s3 = Segment(Node(0.0, -h / 2.0), Node(b, -h / 2.0), t=t)
        sec = Section([s1, s2, s3])

        # Analytical derivations:
        expected_area = t * (2.0 * b + h)
        expected_cx = (b**2) / (2.0 * b + h)
        expected_cy = 0.0
        expected_ix = (t * h**2 / 12.0) * (6.0 * b + h)
        expected_iy = (2.0 / 3.0) * t * (b**3) - (t * b**4) / (2.0 * b + h)
        expected_ixy = 0.0
        expected_i1 = max(expected_ix, expected_iy)
        expected_i2 = min(expected_ix, expected_iy)
        expected_theta_p = 0.0 if expected_ix >= expected_iy else math.pi / 2.0

        assert math.isclose(sec.area, expected_area, rel_tol=1e-12)
        assert math.isclose(sec.cx, expected_cx, rel_tol=1e-12)
        assert math.isclose(sec.cy, expected_cy, abs_tol=1e-12)
        assert math.isclose(sec.Ix, expected_ix, rel_tol=1e-12)
        assert math.isclose(sec.Iy, expected_iy, rel_tol=1e-12)
        assert math.isclose(sec.Ixy, expected_ixy, abs_tol=1e-12)
        assert math.isclose(sec.I1, expected_i1, rel_tol=1e-12)
        assert math.isclose(sec.I2, expected_i2, rel_tol=1e-12)
        assert math.isclose(sec.theta_p, expected_theta_p, abs_tol=1e-12)


class TestTSectionBenchmark:
    """Benchmark Family 4: T-shaped thin-walled section (validates branched-node topology)."""

    def test_t_section_branched(self) -> None:
        # Flange width b = 120, tf = 6.0
        # Web height h = 100, tw = 4.0
        b = 120.0
        tf = 6.0
        h = 100.0
        tw = 4.0

        # Three segments meeting at junction (0, 0)
        s_left = Segment(Node(-b / 2.0, 0.0), Node(0.0, 0.0), t=tf)
        s_right = Segment(Node(b / 2.0, 0.0), Node(0.0, 0.0), t=tf)
        s_web = Segment(Node(0.0, 0.0), Node(0.0, -h), t=tw)
        sec = Section([s_left, s_right, s_web])

        # Analytical derivations:
        expected_area = b * tf + h * tw
        expected_cx = 0.0
        expected_cy = -(tw * h**2) / (2.0 * expected_area)
        expected_iy = (tf * b**3) / 12.0
        expected_ix0 = (tw * h**3) / 3.0
        expected_ix = expected_ix0 - expected_area * (expected_cy**2)
        expected_ixy = 0.0
        expected_i1 = max(expected_ix, expected_iy)
        expected_i2 = min(expected_ix, expected_iy)

        assert math.isclose(sec.area, expected_area, rel_tol=1e-12)
        assert math.isclose(sec.cx, expected_cx, abs_tol=1e-12)
        assert math.isclose(sec.cy, expected_cy, rel_tol=1e-12)
        assert math.isclose(sec.Ix, expected_ix, rel_tol=1e-12)
        assert math.isclose(sec.Iy, expected_iy, rel_tol=1e-12)
        assert math.isclose(sec.Ixy, expected_ixy, abs_tol=1e-12)
        assert math.isclose(sec.I1, expected_i1, rel_tol=1e-12)
        assert math.isclose(sec.I2, expected_i2, rel_tol=1e-12)


class TestAsymmetricSectionBenchmark:
    """Benchmark Family 5: Asymmetric multi-segment open section."""

    def test_asymmetric_z_section(self) -> None:
        # Asymmetric Z-section with independent exact rational arithmetic
        # Bottom flange: (-30, 0) -> (0, 0), t = 2.5
        # Web: (0, 0) -> (0, 80), t = 2.5
        # Top flange: (0, 80) -> (60, 80), t = 2.5
        t = 2.5
        s1 = Segment(Node(-30.0, 0.0), Node(0.0, 0.0), t=t)
        s2 = Segment(Node(0.0, 0.0), Node(0.0, 80.0), t=t)
        s3 = Segment(Node(0.0, 80.0), Node(60.0, 80.0), t=t)
        sec = Section([s1, s2, s3])

        # Exact rational values:
        # Area = 2.5 * (30 + 80 + 60) = 425
        # Qy = 2.5 * (30*(-15) + 80*0 + 60*30) = 3375 -> cx = 3375/425 = 135/17
        # Qx = 2.5 * (30*0 + 80*40 + 60*80) = 20000 -> cy = 20000/425 = 800/17
        # Ix,0 = 2.5 * (80^3 / 3 + 60*80^2) = 4160000 / 3
        # Ix = 4160000/3 - 425*(800/17)^2 = 22720000 / 51
        # Iy,0 = 2.5 * (30^3 / 3 + 60^3 / 3) = 202500
        # Iy = 202500 - 425*(135/17)^2 = 2986875 / 17
        # Ixy,0 = 2.5 * 80 * (60^2 / 2) = 360000
        # Ixy = 360000 - 425*(135/17)*(800/17) = 360000 - 2700000 / 17 = 3420000 / 17
        expected_area = 425.0
        expected_cx = float(Fraction(135, 17))
        expected_cy = float(Fraction(800, 17))
        expected_ix = float(Fraction(22720000, 51))
        expected_iy = float(Fraction(2986875, 17))
        expected_ixy = float(Fraction(3420000, 17))

        diff = expected_ix - expected_iy
        r = math.hypot(diff / 2.0, expected_ixy)
        expected_i1 = (expected_ix + expected_iy) / 2.0 + r
        expected_i2 = (expected_ix + expected_iy) / 2.0 - r
        expected_theta_p = 0.5 * math.atan2(-2.0 * expected_ixy, diff)

        assert math.isclose(sec.area, expected_area, rel_tol=1e-12)
        assert math.isclose(sec.cx, expected_cx, rel_tol=1e-12)
        assert math.isclose(sec.cy, expected_cy, rel_tol=1e-12)
        assert math.isclose(sec.Ix, expected_ix, rel_tol=1e-12)
        assert math.isclose(sec.Iy, expected_iy, rel_tol=1e-12)
        assert math.isclose(sec.Ixy, expected_ixy, rel_tol=1e-12)
        assert math.isclose(sec.I1, expected_i1, rel_tol=1e-12)
        assert math.isclose(sec.I2, expected_i2, rel_tol=1e-12)
        assert math.isclose(sec.theta_p, expected_theta_p, rel_tol=1e-12)


class TestDegenerateIsotropicBenchmark:
    """Degenerate isotropic section test (equal principal moments)."""

    def test_symmetric_cross_section(self) -> None:
        # Cross with 4 equal arms of length b = 50.0, thickness t = 2.0
        b = 50.0
        t = 2.0
        s1 = Segment(Node(0.0, 0.0), Node(b, 0.0), t=t)
        s2 = Segment(Node(0.0, 0.0), Node(0.0, b), t=t)
        s3 = Segment(Node(0.0, 0.0), Node(-b, 0.0), t=t)
        s4 = Segment(Node(0.0, 0.0), Node(0.0, -b), t=t)
        sec = Section([s1, s2, s3, s4])

        # Centroid at (0, 0)
        # Ix = 2 * (t * b^3 / 3) = (2/3)*t*b^3
        # Iy = 2 * (t * b^3 / 3) = (2/3)*t*b^3
        # Ixy = 0
        expected_ix = (2.0 / 3.0) * t * (b**3)
        assert math.isclose(sec.Ix, expected_ix, rel_tol=1e-12)
        assert math.isclose(sec.Iy, expected_ix, rel_tol=1e-12)
        assert math.isclose(sec.Ixy, 0.0, abs_tol=1e-12)
        assert sec.is_degenerate is True
        assert sec.theta_p == 0.0
        assert math.isclose(sec.I1, expected_ix, rel_tol=1e-12)
        assert math.isclose(sec.I2, expected_ix, rel_tol=1e-12)
