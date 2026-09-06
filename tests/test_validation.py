"""Unit tests for geometry and topology validation.

Verifies that invalid geometry and unsupported topology fail explicitly with
testable exceptions (GeometryError, TopologyError).
"""

import pytest

from thinwallx.exceptions import GeometryError, TopologyError
from thinwallx.primitives import Node, Segment
from thinwallx.section import Section


class TestGeometryValidation:
    """Tests for geometric input validation."""

    def test_empty_section_rejected(self) -> None:
        with pytest.raises(TopologyError, match="at least one segment"):
            Section([])

    def test_non_finite_coordinates_rejected(self) -> None:
        with pytest.raises(GeometryError, match="must be finite"):
            Node(float("nan"), 0.0)
        with pytest.raises(GeometryError, match="must be finite"):
            Node(1.0, float("inf"))

    def test_non_positive_thickness_rejected(self) -> None:
        n1 = Node(0.0, 0.0)
        n2 = Node(1.0, 0.0)
        with pytest.raises(GeometryError, match="strictly positive"):
            Segment(n1, n2, t=0.0)
        with pytest.raises(GeometryError, match="strictly positive"):
            Segment(n1, n2, t=-2.0)

    def test_non_finite_thickness_rejected(self) -> None:
        n1 = Node(0.0, 0.0)
        n2 = Node(1.0, 0.0)
        with pytest.raises(GeometryError, match="must be finite"):
            Segment(n1, n2, t=float("nan"))
        with pytest.raises(GeometryError, match="must be finite"):
            Segment(n1, n2, t=float("inf"))

    def test_zero_length_segment_rejected(self) -> None:
        n1 = Node(2.0, 3.0)
        n2 = Node(2.0, 3.0)
        with pytest.raises(GeometryError, match="strictly positive"):
            Segment(n1, n2, t=1.0)

    def test_coincident_endpoints_within_tolerance_rejected(self) -> None:
        n1 = Node(0.0, 0.0)
        n2 = Node(1e-10, 0.0)  # Length < 1e-9 tolerance
        # Creating segment alone: 1e-10 > 1e-12 so segment might pass segment post_init,
        # but Section validation with node_tolerance=1e-9 clusters them as coincident!
        with pytest.raises(GeometryError, match="coincident endpoints"):
            Section([Segment(n1, n2, t=1.0)], node_tolerance=1e-9)


class TestTopologyValidation:
    """Tests for topology validation (duplicate segments, disconnected geometry, closed loops)."""

    def test_duplicate_segments_same_direction_rejected(self) -> None:
        n1 = Node(0.0, 0.0)
        n2 = Node(10.0, 0.0)
        seg1 = Segment(n1, n2, t=2.0)
        seg2 = Segment(n1, n2, t=3.0)
        with pytest.raises(TopologyError, match="Duplicate segment detected"):
            Section([seg1, seg2])

    def test_duplicate_segments_reversed_direction_rejected(self) -> None:
        n1 = Node(0.0, 0.0)
        n2 = Node(10.0, 0.0)
        seg1 = Segment(n1, n2, t=2.0)
        seg2 = Segment(n2, n1, t=2.0)
        with pytest.raises(TopologyError, match="Duplicate segment detected"):
            Section([seg1, seg2])

    def test_disconnected_segments_rejected(self) -> None:
        # Two parallel independent segments
        seg1 = Segment(Node(0.0, 0.0), Node(10.0, 0.0), t=1.0)
        seg2 = Segment(Node(0.0, 5.0), Node(10.0, 5.0), t=1.0)
        with pytest.raises(TopologyError, match="Disconnected geometry detected"):
            Section([seg1, seg2])

    def test_disconnected_subtrees_rejected(self) -> None:
        # One T-section and one separate detached segment
        t1 = Segment(Node(-5.0, 10.0), Node(0.0, 10.0), t=1.0)
        t2 = Segment(Node(0.0, 10.0), Node(5.0, 10.0), t=1.0)
        t3 = Segment(Node(0.0, 10.0), Node(0.0, 0.0), t=1.0)
        isolated = Segment(Node(100.0, 0.0), Node(100.0, 10.0), t=1.0)
        with pytest.raises(TopologyError, match="Disconnected geometry detected"):
            Section([t1, t2, t3, isolated])

    def test_closed_loop_triangle_rejected(self) -> None:
        n1 = Node(0.0, 0.0)
        n2 = Node(10.0, 0.0)
        n3 = Node(5.0, 8.66)
        seg1 = Segment(n1, n2, t=1.0)
        seg2 = Segment(n2, n3, t=1.0)
        seg3 = Segment(n3, n1, t=1.0)
        with pytest.raises(TopologyError, match="Closed-loop geometry detected"):
            Section([seg1, seg2, seg3])

    def test_closed_loop_box_section_rejected(self) -> None:
        n1 = Node(0.0, 0.0)
        n2 = Node(10.0, 0.0)
        n3 = Node(10.0, 10.0)
        n4 = Node(0.0, 10.0)
        seg1 = Segment(n1, n2, t=1.0)
        seg2 = Segment(n2, n3, t=1.0)
        seg3 = Segment(n3, n4, t=1.0)
        seg4 = Segment(n4, n1, t=1.0)
        with pytest.raises(TopologyError, match="Closed-loop geometry detected"):
            Section([seg1, seg2, seg3, seg4])

    def test_valid_branched_topology_accepted(self) -> None:
        # A valid branched T-section must pass validation without error
        n_left = Node(-10.0, 20.0)
        n_mid = Node(0.0, 20.0)
        n_right = Node(10.0, 20.0)
        n_bot = Node(0.0, 0.0)
        sec = Section([
            Segment(n_left, n_mid, t=2.0),
            Segment(n_right, n_mid, t=2.0),
            Segment(n_mid, n_bot, t=2.0),
        ])
        assert len(sec.segments) == 3
        assert len(sec.nodes) == 4
