"""Tests for planar cell topology extraction and validation (ThinWallX v0.5)."""

import math
import numpy as np
import pytest

from thinwallx.cells import Cell, CellTopology, extract_cell_topology
from thinwallx.closed_section import ClosedSection
from thinwallx.exceptions import GeometryError, SingularSectionError, TopologyError
from thinwallx.primitives import Node, Segment


def make_box_segments(a: float, b: float, t: float) -> list[Segment]:
    """Helper creating a 4-segment rectangular box of width a, height b."""
    n1 = Node(0.0, 0.0)
    n2 = Node(a, 0.0)
    n3 = Node(a, b)
    n4 = Node(0.0, b)
    return [
        Segment(n1, n2, t=t),  # bottom (seg 0)
        Segment(n2, n3, t=t),  # right (seg 1)
        Segment(n3, n4, t=t),  # top (seg 2)
        Segment(n4, n1, t=t),  # left (seg 3)
    ]


def make_twocell_box_segments(b: float, h: float, t: float) -> list[Segment]:
    """Helper creating two rectangular cells side-by-side with shared internal web.

    Cell 1: [0, b] x [0, h]
    Cell 2: [b, 2b] x [0, h]
    Shared web: x=b from (b, 0) to (b, h).
    """
    n1 = Node(0.0, 0.0)
    n2 = Node(b, 0.0)
    n3 = Node(2.0 * b, 0.0)
    n4 = Node(2.0 * b, h)
    n5 = Node(b, h)
    n6 = Node(0.0, h)

    # 7 segments total:
    # seg 0: (0, 0) -> (b, 0) bottom-left
    # seg 1: (b, 0) -> (2b, 0) bottom-right
    # seg 2: (2b, 0) -> (2b, h) right outer
    # seg 3: (2b, h) -> (b, h) top-right
    # seg 4: (b, h) -> (0, h) top-left
    # seg 5: (0, h) -> (0, 0) left outer
    # seg 6: (b, 0) -> (b, h) shared internal web
    return [
        Segment(n1, n2, t=t),
        Segment(n2, n3, t=t),
        Segment(n3, n4, t=t),
        Segment(n4, n5, t=t),
        Segment(n5, n6, t=t),
        Segment(n6, n1, t=t),
        Segment(n2, n5, t=t),  # shared web
    ]


class TestCellTopologyExtraction:
    """Tests for planar cell extraction and validation."""

    def test_single_cell_box_topology(self) -> None:
        a, b, t = 100.0, 50.0, 2.0
        segs = make_box_segments(a, b, t)
        topo = extract_cell_topology(segs)

        assert topo.cell_count == 1
        assert topo.segment_count == 4
        assert topo.node_count == 4

        cell = topo.cells[0]
        assert cell.id == 0
        assert math.isclose(cell.area, a * b, rel_tol=1e-12)
        assert math.isclose(cell.perimeter, 2.0 * (a + b), rel_tol=1e-12)

        # Incidence matrix B: shape (1, 4), values in {-1, 1}
        assert topo.B.shape == (1, 4)
        assert np.all(np.abs(topo.B) == 1.0)

        # Cell flexibility matrix H: scalar = 2*(a + b)/t
        expected_H = 2.0 * (a + b) / t
        assert topo.H.shape == (1, 1)
        assert math.isclose(float(topo.H[0, 0]), expected_H, rel_tol=1e-12)

    def test_two_cell_box_topology_and_shared_web(self) -> None:
        b, h, t = 60.0, 40.0, 3.0
        segs = make_twocell_box_segments(b, h, t)
        topo = extract_cell_topology(segs)

        # 6 nodes, 7 segments => m = 7 - 6 + 1 = 2
        assert topo.cell_count == 2
        assert topo.segment_count == 7
        assert topo.node_count == 6

        # Both cells have area b * h
        for cell in topo.cells:
            assert math.isclose(cell.area, b * h, rel_tol=1e-12)

        # Shared web is segment index 6
        # Adjacent cells must have opposite signs on shared web
        b_web = topo.B[:, 6]
        assert np.count_nonzero(b_web) == 2
        c1, c2 = np.where(b_web != 0.0)[0]
        assert topo.B[c1, 6] == -topo.B[c2, 6]

        # Analytical H check per ACTIVE_PHASE.md:
        # H = (1/t) * [[2*(b+h), -h], [-h, 2*(b+h)]]
        expected_H = (1.0 / t) * np.array(
            [[2.0 * (b + h), -h], [-h, 2.0 * (b + h)]], dtype=float
        )
        np.testing.assert_allclose(topo.H, expected_H, rtol=1e-12, atol=1e-12)

    def test_segment_order_invariance(self) -> None:
        """Reordering segments must yield the same physical cells and areas."""
        a, b, t = 100.0, 50.0, 2.0
        segs = make_box_segments(a, b, t)
        topo1 = extract_cell_topology(segs)

        # Reverse order of segments
        segs_rev = list(reversed(segs))
        topo2 = extract_cell_topology(segs_rev)

        assert topo1.cell_count == topo2.cell_count == 1
        assert math.isclose(topo1.cells[0].area, topo2.cells[0].area, rel_tol=1e-12)
        assert math.isclose(topo1.cells[0].perimeter, topo2.cells[0].perimeter, rel_tol=1e-12)
        assert math.isclose(float(topo1.H[0, 0]), float(topo2.H[0, 0]), rel_tol=1e-12)

    def test_segment_direction_reversal(self) -> None:
        """Reversing individual segment endpoints must preserve cell area and H."""
        a, b, t = 100.0, 50.0, 2.0
        segs = make_box_segments(a, b, t)
        # Reverse p1 and p2 of segments 0 and 2
        segs_flipped = [
            Segment(segs[0].p2, segs[0].p1, t=t),
            segs[1],
            Segment(segs[2].p2, segs[2].p1, t=t),
            segs[3],
        ]
        topo = extract_cell_topology(segs_flipped)
        assert topo.cell_count == 1
        assert math.isclose(topo.cells[0].area, a * b, rel_tol=1e-12)
        assert math.isclose(float(topo.H[0, 0]), 2.0 * (a + b) / t, rel_tol=1e-12)

    def test_mixed_open_closed_topology_explicitly_rejected(self) -> None:
        """Section containing closed cell + dangling open branch must raise GeometryError."""
        a, b, t = 80.0, 40.0, 2.0
        box = make_box_segments(a, b, t)
        # Add a dangling antenna branch attached to node (a, b)
        antenna = Segment(Node(a, b), Node(a + 20.0, b + 20.0), t=t)
        mixed_segs = box + [antenna]

        with pytest.raises(GeometryError, match="Mixed open/closed topology detected"):
            extract_cell_topology(mixed_segs)

        with pytest.raises(GeometryError, match="Mixed open/closed topology detected"):
            ClosedSection(mixed_segs)

    def test_pure_open_tree_rejected_by_closed_topology(self) -> None:
        """Pure open section (m = 0) must be rejected by extract_cell_topology."""
        s1 = Segment(Node(0.0, 0.0), Node(10.0, 0.0), t=2.0)
        s2 = Segment(Node(10.0, 0.0), Node(10.0, 10.0), t=2.0)
        with pytest.raises(TopologyError, match="No closed cells detected"):
            extract_cell_topology([s1, s2])

        with pytest.raises(TopologyError, match="No closed cells detected"):
            ClosedSection([s1, s2])

    def test_collapsed_cell_rejected(self) -> None:
        """Cell with near-zero enclosed area relative to perimeter^2 must be rejected."""
        L = 1e5
        h = 1e-7  # > node_tolerance 1e-9, but Area = 0.01 <= 1e-12 * (2e5)^2 = 0.04
        n1 = Node(0.0, 0.0)
        n2 = Node(L, 0.0)
        n3 = Node(L, h)
        n4 = Node(0.0, h)
        collapsed_segs = [
            Segment(n1, n2, t=1.0),
            Segment(n2, n3, t=1.0),
            Segment(n3, n4, t=1.0),
            Segment(n4, n1, t=1.0),
        ]
        with pytest.raises(GeometryError, match="collapsed/zero area"):
            extract_cell_topology(collapsed_segs)

    def test_closed_section_wrapper_properties(self) -> None:
        """Verify ClosedSection accessors."""
        a, b, t = 120.0, 60.0, 4.0
        sec = ClosedSection.from_segments(make_box_segments(a, b, t))
        assert sec.is_closed is True
        assert sec.cell_count == 1
        assert len(sec.cells) == 1
        assert math.isclose(sec.cells[0].area, a * b)
        assert math.isclose(sec.area, 2.0 * (a + b) * t)
