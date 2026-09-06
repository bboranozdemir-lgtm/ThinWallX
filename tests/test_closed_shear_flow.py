"""Tests for closed-section transverse shear flow (Sectalix v0.5)."""

import math
import numpy as np
import pytest

from sectalix.closed_section import ClosedSection
from sectalix.closed_shear_flow import calculate_closed_shear_flow
from sectalix.exceptions import GeometryError, SingularSectionError
from sectalix.primitives import Node, Segment
from sectalix.shear_load import ShearLoad


def make_box_section(a: float = 100.0, b: float = 50.0, t: float = 2.0) -> ClosedSection:
    """Helper returning a ClosedSection for a rectangular box a x b with thickness t."""
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


def make_twocell_box_section(
    b: float = 60.0, h: float = 40.0, t: float = 3.0
) -> ClosedSection:
    """Helper returning a ClosedSection for two side-by-side cells."""
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
        Segment(n2, n5, t=t),  # shared web
    ])


class TestClosedTransverseShearFlow:
    """Tests for closed-section transverse shear flow calculation."""

    def test_single_cell_resultant_recovery(self) -> None:
        """Final compatible shear flow must recover applied resultant Vx, Vy exactly."""
        sec = make_box_section(100.0, 60.0, 2.5)
        for vx, vy in [(1000.0, 0.0), (0.0, 500.0), (750.0, -400.0)]:
            res = sec.calculate_shear_flow(vx, vy)
            np.testing.assert_allclose(
                res.recovered_resultant, [vx, vy], rtol=1e-11, atol=1e-10
            )
            assert res.resultant_error < 1e-9

    def test_single_cell_zero_twist_compatibility(self) -> None:
        """Zero-twist compatibility oint (q/t) ds = 0 must hold for all load cases."""
        sec = make_box_section(80.0, 50.0, 3.0)
        for vx, vy in [(500.0, 0.0), (0.0, 1000.0), (1200.0, 800.0)]:
            res = sec.calculate_shear_flow(vx, vy)
            assert res.max_compatibility_residual < 1e-10

    def test_cut_continuity(self) -> None:
        """Physical flow must be continuous across the virtual cut on chord segments."""
        sec = make_box_section(120.0, 60.0, 4.0)
        res = sec.calculate_shear_flow(0.0, 1000.0, cut_param=0.35)

        for flow in res.segment_flows:
            if flow.is_chord:
                # Check value from left and right of cut_param
                eps = 1e-12
                q_left = flow.q_at_xi(0.35 - eps)
                q_right = flow.q_at_xi(0.35 + eps)
                q_at_cut = flow.q_at_xi(0.35)
                assert math.isclose(q_left, q_at_cut, rel_tol=1e-10, abs_tol=1e-10)
                assert math.isclose(q_right, q_at_cut, rel_tol=1e-10, abs_tol=1e-10)

    def test_cut_location_independence(self) -> None:
        """Varying the interior cut parameter location must leave final flow unchanged."""
        sec = make_box_section(100.0, 50.0, 2.0)
        res_mid = sec.calculate_shear_flow(400.0, 600.0, cut_param=0.5)
        res_quarter = sec.calculate_shear_flow(400.0, 600.0, cut_param=0.25)
        res_threequarter = sec.calculate_shear_flow(400.0, 600.0, cut_param=0.75)

        # Check evaluations at multiple points on all segments
        for i in range(len(sec.segments)):
            for xi in [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]:
                q_mid = res_mid.segment_flows[i].q_at_xi(xi)
                q_qtr = res_quarter.segment_flows[i].q_at_xi(xi)
                q_tq = res_threequarter.segment_flows[i].q_at_xi(xi)
                assert math.isclose(q_mid, q_qtr, rel_tol=1e-10, abs_tol=1e-10)
                assert math.isclose(q_mid, q_tq, rel_tol=1e-10, abs_tol=1e-10)

    def test_spanning_tree_independence(self) -> None:
        """Choosing different valid spanning trees (different chord cut) must yield identical physical flow."""
        sec = make_box_section(100.0, 50.0, 2.0)
        # Tree 1: default spanning tree (edges 0, 1, 2 as tree, edge 3 as chord)
        res_tree1 = calculate_closed_shear_flow(sec, 500.0, 300.0, spanning_tree_edges=[0, 1, 2])
        # Tree 2: edges 1, 2, 3 as tree, edge 0 as chord
        res_tree2 = calculate_closed_shear_flow(sec, 500.0, 300.0, spanning_tree_edges=[1, 2, 3])

        for i in range(len(sec.segments)):
            for xi in [0.0, 0.2, 0.5, 0.8, 1.0]:
                q1 = res_tree1.segment_flows[i].q_at_xi(xi)
                q2 = res_tree2.segment_flows[i].q_at_xi(xi)
                assert math.isclose(q1, q2, rel_tol=1e-10, abs_tol=1e-10)

        # Torque and shear center must also match exactly
        assert math.isclose(res_tree1.torque, res_tree2.torque, rel_tol=1e-10, abs_tol=1e-10)

    def test_multi_cell_resultant_and_compatibility(self) -> None:
        """Two-cell section under combined shear must satisfy both cell compatibilities and recover resultant."""
        sec = make_twocell_box_section(50.0, 30.0, 2.0)
        vx, vy = 800.0, -600.0
        res = sec.calculate_shear_flow(vx, vy)

        np.testing.assert_allclose(res.recovered_resultant, [vx, vy], rtol=1e-11, atol=1e-10)
        assert len(res.compatibility_residuals) == 2
        np.testing.assert_allclose(res.compatibility_residuals, [0.0, 0.0], atol=1e-10)

    def test_doubly_symmetric_box_shear_center(self) -> None:
        """For doubly symmetric box, shear center coincides with centroid: S = C."""
        sec = make_box_section(100.0, 60.0, 3.0)
        cx, cy = sec.centroid
        xs, ys = sec.shear_center
        ex, ey = sec.sc_offset

        assert math.isclose(ex, 0.0, abs_tol=1e-10)
        assert math.isclose(ey, 0.0, abs_tol=1e-10)
        assert math.isclose(xs, cx, abs_tol=1e-10)
        assert math.isclose(ys, cy, abs_tol=1e-10)

    def test_invalid_shear_load_rejected(self) -> None:
        sec = make_box_section()
        with pytest.raises(GeometryError, match="must be finite"):
            sec.calculate_shear_flow(float("nan"), 0.0)
        with pytest.raises(GeometryError, match="must be finite"):
            sec.calculate_shear_flow(0.0, float("inf"))

    def test_exact_evaluation_api(self) -> None:
        """Section 8: Exact evaluation API q_i(xi) and q_i(s) on physical walls."""
        sec = make_box_section(100.0, 50.0, 2.0)
        res = sec.calculate_shear_flow(200.0, 400.0)
        for flow in res.segment_flows:
            L = flow.length
            for xi in [0.0, 0.25, 0.5, 0.75, 1.0]:
                q_xi = flow.q_at_xi(xi)
                q_s = flow.q_at_s(xi * L)
                assert math.isclose(q_xi, q_s, rel_tol=1e-12)

        # Also test out of range xi fails explicitly
        with pytest.raises(GeometryError, match="must lie in"):
            res.segment_flows[0].q_at_xi(-0.5)
        with pytest.raises(GeometryError, match="must lie in"):
            res.segment_flows[0].q_at_xi(1.5)
        with pytest.raises(GeometryError, match="must be finite"):
            res.segment_flows[0].q_at_xi(float("nan"))

