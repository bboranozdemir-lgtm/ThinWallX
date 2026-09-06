"""Unit tests for Node and Segment primitives."""

import math
import pytest

from sectalix.exceptions import GeometryError
from sectalix.primitives import Node, Segment


class TestNode:
    """Tests for Node primitive."""

    def test_valid_node(self) -> None:
        n = Node(1.5, -2.5, id="N1")
        assert n.x == 1.5
        assert n.y == -2.5
        assert n.id == "N1"
        assert n.coords == (1.5, -2.5)

    def test_node_rejects_nan_and_inf(self) -> None:
        with pytest.raises(GeometryError, match="must be finite"):
            Node(float("nan"), 0.0)
        with pytest.raises(GeometryError, match="must be finite"):
            Node(0.0, float("inf"))
        with pytest.raises(GeometryError, match="must be finite"):
            Node(float("-inf"), float("nan"))

    def test_distance_and_is_close(self) -> None:
        n1 = Node(0.0, 0.0)
        n2 = Node(3.0, 4.0)
        assert math.isclose(n1.distance_to(n2), 5.0)
        assert not n1.is_close(n2, tol=1e-3)

        n3 = Node(1e-10, -1e-10)
        assert n1.is_close(n3, tol=1e-9)


class TestSegment:
    """Tests for Segment primitive and its exact closed-form line integrals."""

    def test_valid_segment(self) -> None:
        p1 = Node(0.0, 0.0)
        p2 = Node(3.0, 4.0)
        seg = Segment(p1, p2, t=2.0, id="S1")
        assert seg.length == 5.0
        assert seg.area == 10.0
        assert seg.t == 2.0
        assert seg.id == "S1"

    def test_segment_rejects_non_positive_thickness(self) -> None:
        p1 = Node(0.0, 0.0)
        p2 = Node(1.0, 1.0)
        with pytest.raises(GeometryError, match="strictly positive"):
            Segment(p1, p2, t=0.0)
        with pytest.raises(GeometryError, match="strictly positive"):
            Segment(p1, p2, t=-1.5)

    def test_segment_rejects_non_finite_thickness(self) -> None:
        p1 = Node(0.0, 0.0)
        p2 = Node(1.0, 1.0)
        with pytest.raises(GeometryError, match="must be finite"):
            Segment(p1, p2, t=float("nan"))
        with pytest.raises(GeometryError, match="must be finite"):
            Segment(p1, p2, t=float("inf"))

    def test_segment_rejects_zero_length(self) -> None:
        p1 = Node(2.0, 3.0)
        p2 = Node(2.0, 3.0)
        with pytest.raises(GeometryError, match="strictly positive"):
            Segment(p1, p2, t=1.0)

    def test_exact_integrals_horizontal_segment(self) -> None:
        # Segment from (0, 2) to (4, 2) of length L=4
        # \int x ds = 4 * (0 + 4)/2 = 8
        # \int y ds = 4 * (2 + 2)/2 = 8
        # \int x^2 ds = (4/3) * (0 + 0 + 16) = 64/3
        # \int y^2 ds = (4/3) * (4 + 4 + 4) = 16
        # \int xy ds = (4/6) * (0 + 0 + 8 + 2*4*2) = (4/6) * 24 = 16.0
        seg = Segment(Node(0.0, 2.0), Node(4.0, 2.0), t=1.0)
        assert math.isclose(seg.length, 4.0)
        assert math.isclose(seg.int_x(), 8.0)
        assert math.isclose(seg.int_y(), 8.0)
        assert math.isclose(seg.int_x2(), 64.0 / 3.0)
        assert math.isclose(seg.int_y2(), 16.0)
        assert math.isclose(seg.int_xy(), 16.0)

    def test_exact_integrals_inclined_segment(self) -> None:
        # Segment from (1, 2) to (4, 6): dx=3, dy=4, L=5
        # \int x ds = 5 * (1 + 4)/2 = 12.5
        # \int y ds = 5 * (2 + 6)/2 = 20.0
        # \int x^2 ds = (5/3) * (1 + 4 + 16) = 5 * 21 / 3 = 35.0
        # \int y^2 ds = (5/3) * (4 + 12 + 36) = 5 * 52 / 3 = 260/3
        # \int xy ds = (5/6) * (2*1*2 + 1*6 + 4*2 + 2*4*6) = (5/6) * (4 + 6 + 8 + 48) = 5 * 66 / 6 = 55.0
        seg = Segment(Node(1.0, 2.0), Node(4.0, 6.0), t=1.5)
        assert math.isclose(seg.length, 5.0)
        assert math.isclose(seg.int_x(), 12.5)
        assert math.isclose(seg.int_y(), 20.0)
        assert math.isclose(seg.int_x2(), 35.0)
        assert math.isclose(seg.int_y2(), 260.0 / 3.0)
        assert math.isclose(seg.int_xy(), 55.0)

    def test_reversed_segment_integral_invariance(self) -> None:
        seg = Segment(Node(1.0, 2.0), Node(4.0, 6.0), t=1.5)
        rev = seg.reversed()
        assert math.isclose(seg.int_x(), rev.int_x())
        assert math.isclose(seg.int_y(), rev.int_y())
        assert math.isclose(seg.int_x2(), rev.int_x2())
        assert math.isclose(seg.int_y2(), rev.int_y2())
        assert math.isclose(seg.int_xy(), rev.int_xy())
