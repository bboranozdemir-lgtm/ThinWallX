"""Tests for mixed open-closed topology decomposition (ThinWallX v0.6).

Covers:
- T01: Single-cell + branch; multi-cell + multi-branch (E_c, E_o, junctions, rank).
- T03: Bridge connecting two closed loops (k_c = 2, bridge in E_o, no free tip).
- T04: Single articulation node connecting two closed loops (figure-8 / bowtie).
- T05: Additional path connecting two nodes of a cell (becomes cycle, not open).
- T06: Pure open and pure closed reductions in MixedSection.
- T21: Disconnected geometry and invalid section error handling.
- T22: API immutability and caching consistency.
"""

from __future__ import annotations

import math
import pytest
import numpy as np

from thinwallx.exceptions import GeometryError, SingularSectionError, TopologyError
from thinwallx.mixed_section import MixedSection
from thinwallx.mixed_topology import MixedTopology, extract_mixed_topology
from thinwallx.primitives import Node, Segment


def make_box_segments(x0: float, y0: float, w: float, h: float, t: float) -> list[Segment]:
    """Helper to create a rectangular box starting at (x0, y0)."""
    p1 = Node(x0, y0)
    p2 = Node(x0 + w, y0)
    p3 = Node(x0 + w, y0 + h)
    p4 = Node(x0, y0 + h)
    return [
        Segment(p1, p2, t=t),
        Segment(p2, p3, t=t),
        Segment(p3, p4, t=t),
        Segment(p4, p1, t=t),
    ]


class TestMixedTopology:
    """T01, T03, T04, T05, T06, T21, T22 verification."""

    def test_t01_single_cell_plus_single_antenna(self) -> None:
        """T01: Single closed box with an open antenna branch."""
        box = make_box_segments(0.0, 0.0, 4.0, 2.0, 0.02)
        # Antenna from node (4, 2) to (6, 2)
        antenna = Segment(Node(4.0, 2.0), Node(6.0, 2.0), t=0.01)
        segs = box + [antenna]

        sec = MixedSection(segs)
        topo = sec.mixed_topology

        assert topo.is_mixed is True
        assert topo.is_closed is False
        assert topo.is_open is False
        assert topo.topology_type == "mixed"
        assert topo.cycle_rank == 1
        assert topo.cell_count == 1
        assert len(topo.closed_edges) == 4
        assert len(topo.open_edges) == 1
        assert topo.open_edges == (4,)
        assert len(topo.junction_nodes) == 1
        assert len(topo.free_tip_nodes) == 1

        # Check B matrix: (1, 5) with zero on antenna column
        assert topo.B.shape == (1, 5)
        assert topo.B[0, 4] == 0.0
        assert math.isclose(topo.cell_areas[0], 8.0, rel_tol=1e-12)

    def test_t01_two_cell_box_with_two_cantilevers(self) -> None:
        """T01: Two-cell box with two open cantilever flanges (B2 topology)."""
        # Outer box (0,0)-(8,2) with web at x=3
        # Left cell [0, 3] x [0, 2], Right cell [3, 8] x [0, 2]
        n00 = Node(0.0, 0.0)
        n30 = Node(3.0, 0.0)
        n80 = Node(8.0, 0.0)
        n82 = Node(8.0, 2.0)
        n32 = Node(3.0, 2.0)
        n02 = Node(0.0, 2.0)

        # Cantilever nodes
        nm22 = Node(-2.0, 2.0)
        n92 = Node(9.0, 2.0)

        t_out = 0.02
        t_web = 0.025
        t_cant = 0.01

        segs = [
            Segment(n00, n30, t=t_out),  # 0
            Segment(n30, n80, t=t_out),  # 1
            Segment(n80, n82, t=t_out),  # 2
            Segment(n82, n32, t=t_out),  # 3
            Segment(n32, n02, t=t_out),  # 4
            Segment(n02, n00, t=t_out),  # 5
            Segment(n30, n32, t=t_web),  # 6 (common web)
            Segment(n02, nm22, t=t_cant),  # 7 (cantilever left)
            Segment(n82, n92, t=t_cant),  # 8 (cantilever right)
        ]

        sec = MixedSection(segs)
        topo = sec.mixed_topology

        assert topo.cycle_rank == 2
        assert topo.cell_count == 2
        assert len(topo.closed_edges) == 7
        assert len(topo.open_edges) == 2
        assert set(topo.open_edges) == {7, 8}
        assert len(topo.junction_nodes) == 2
        assert len(topo.free_tip_nodes) == 2

        # Check enclosed cell areas
        areas = sorted(topo.cell_areas)
        assert math.isclose(areas[0], 6.0, rel_tol=1e-12)
        assert math.isclose(areas[1], 10.0, rel_tol=1e-12)

        # Check B shape
        assert topo.B.shape == (2, 9)
        assert topo.B[0, 7] == 0.0 and topo.B[0, 8] == 0.0
        assert topo.B[1, 7] == 0.0 and topo.B[1, 8] == 0.0

    def test_t03_bridge_connecting_two_separate_closed_loops(self) -> None:
        """T03: Two closed loops connected by an open bridge (k_c = 2)."""
        box1 = make_box_segments(0.0, 0.0, 2.0, 2.0, 0.02)  # edges 0-3
        box2 = make_box_segments(5.0, 0.0, 2.0, 2.0, 0.02)  # edges 4-7
        # Bridge connecting (2, 1) on box1 to (5, 1) on box2
        # Note: box1 has vertices (0,0), (2,0), (2,2), (0,2). To connect at (2,1), split right wall
        # Or connect directly at vertex (2, 2) of box1 to (5, 2) of box2!
        bridge = Segment(Node(2.0, 2.0), Node(5.0, 2.0), t=0.01)
        segs = box1 + box2 + [bridge]

        sec = MixedSection(segs)
        topo = sec.mixed_topology

        assert topo.is_mixed is True
        assert len(topo.open_edges) == 1
        assert topo.open_edges == (8,)  # bridge edge is edge 8
        assert len(topo.closed_edges) == 8
        assert len(topo.cyclic_components) == 2  # k_c = 2
        assert topo.cycle_rank == 2  # 1 from each loop
        assert topo.cell_count == 2
        assert len(topo.junction_nodes) == 2
        assert len(topo.free_tip_nodes) == 0  # no free tip: bridge connects two closed loops!

        # Global sectorial coordinate continuity and zero-mean verification
        tw_res = sec.torsion_properties()
        assert math.isfinite(tw_res.J_total)
        assert math.isfinite(tw_res.Cw)
        assert tw_res.Cw >= 0.0

        # Zero-mean sectorial coordinate verification: int_A omega dA = 0
        total_omega_integral = sum(
            seg.t * seg.length * 0.5 * (tw_res.node_omega[u] + tw_res.node_omega[v])
            for (u, v), seg in zip(topo.canonical_edges, sec.segments)
        )
        assert math.isclose(total_omega_integral, 0.0, abs_tol=1e-12)

        # Global continuity across all segments and junction nodes
        for (u, v), sw in zip(topo.canonical_edges, tw_res.segment_warpings):
            d_omega = tw_res.node_omega[v] - tw_res.node_omega[u]
            expected_d_omega = sw.p * sw.segment.length
            assert math.isclose(d_omega, expected_d_omega, rel_tol=1e-10, abs_tol=1e-12)

    def test_t03_multiscale_barbell_blocks_conditioning(self) -> None:
        """T03b: Barbell with multiscale blocks (L/t=1 vs L/t=2^300) does not falsely trigger SingularSectionError."""
        t2 = 2.0 ** -300
        box1_split = [
            Segment(Node(0, 0), Node(1, 0), t=1.0),
            Segment(Node(1, 0), Node(1, 0.5), t=1.0),
            Segment(Node(1, 0.5), Node(1, 1), t=1.0),
            Segment(Node(1, 1), Node(0, 1), t=1.0),
            Segment(Node(0, 1), Node(0, 0), t=1.0),
        ]
        bridge = [
            Segment(Node(1, 0.5), Node(6, 0.5), t=0.1)
        ]
        box2_split = [
            Segment(Node(6, 0), Node(7, 0), t=t2),
            Segment(Node(7, 0), Node(7, 1), t=t2),
            Segment(Node(7, 1), Node(6, 1), t=t2),
            Segment(Node(6, 1), Node(6, 0.5), t=t2),
            Segment(Node(6, 0.5), Node(6, 0), t=t2),
        ]
        sec_barbell = MixedSection(box1_split + bridge + box2_split)
        topo = sec_barbell.mixed_topology
        assert topo.cycle_rank == 2
        assert topo.cell_count == 2
        assert len(topo.cyclic_components) == 2
        assert len(topo.open_edges) == 1

    def test_t03c_extreme_barbell_uses_independent_block_solves(self) -> None:
        """T03c: H blocks differing by more than 1074 binary exponents remain independent."""
        length = 2.0 ** 100
        thin = 2.0 ** -1000
        large_cell = [
            Segment(Node(0.0, 0.0), Node(length, 0.0), t=thin),
            Segment(Node(length, 0.0), Node(length, length), t=thin),
            Segment(Node(length, length), Node(0.0, length), t=thin),
            Segment(Node(0.0, length), Node(0.0, 0.0), t=thin),
        ]
        small_cell = [
            Segment(Node(1.0, -1.0), Node(2.0, -1.0), t=0.1),
            Segment(Node(2.0, -1.0), Node(2.0, -2.0), t=0.1),
            Segment(Node(2.0, -2.0), Node(1.0, -2.0), t=0.1),
            Segment(Node(1.0, -2.0), Node(1.0, -1.0), t=0.1),
        ]
        bridge = [Segment(Node(0.0, 0.0), Node(1.0, -1.0), t=0.01)]
        section = MixedSection(large_cell + small_cell + bridge)
        topology = section.mixed_topology

        assert topology.block_cell_indices == ((0,), (1,))
        assert topology.e_max_H_block == (1101, 4)
        assert all(matrix.shape == (1, 1) for matrix in topology.H_block_scaled)

        expected_j_open = math.sqrt(2.0) * 0.01**3 / 3.0
        expected_j_total = 0.1 + expected_j_open
        assert section.J_BB == pytest.approx(0.1, rel=1e-14, abs=0.0)
        assert section.J_open == pytest.approx(expected_j_open, rel=1e-14, abs=0.0)
        assert section.J_total == pytest.approx(expected_j_total, rel=1e-14, abs=0.0)

        # A representable extreme-load slice also exercises independent q0 block solves.
        result = section.calculate_shear_flow(vx=2.0 ** -700, vy=0.0)
        assert np.all(np.isfinite(result.redundant_cell_flows))
        assert np.all(result.redundant_cell_flows != 0.0)
        assert result.relative_residual <= 1e-10

    def test_t04_two_loops_meeting_at_single_articulation_node(self) -> None:
        """T04: Two closed loops meeting at a single articulation vertex (bowtie / figure-8)."""
        # Loop 1: (0,0)-(2,0)-(2,2)-(0,2)
        # Loop 2: (2,2)-(4,2)-(4,4)-(2,4)
        # Meeting at vertex (2, 2)
        l1 = make_box_segments(0.0, 0.0, 2.0, 2.0, 0.02)
        l2 = make_box_segments(2.0, 2.0, 2.0, 2.0, 0.02)
        segs = l1 + l2

        sec = MixedSection(segs)
        topo = sec.mixed_topology

        assert topo.cycle_rank == 2
        assert topo.cell_count == 2
        assert len(topo.closed_edges) == 8
        assert len(topo.open_edges) == 0
        assert topo.is_closed is True
        assert len(topo.cyclic_components) == 1  # connected via articulation vertex
        assert len(topo.cyclic_blocks) == 2  # two biconnected blocks!
        assert len(topo.articulation_nodes) == 1

    def test_t05_additional_path_forms_new_cycle_not_open(self) -> None:
        """T05: Additional path connecting two nodes of a cell forms a new cycle."""
        box = make_box_segments(0.0, 0.0, 2.0, 2.0, 0.02)
        # Path connecting (0,0) and (2,2) externally via (0, -1), (3, -1), (3, 2)
        p1 = Segment(Node(0.0, 0.0), Node(0.0, -1.0), t=0.02)
        p2 = Segment(Node(0.0, -1.0), Node(3.0, -1.0), t=0.02)
        p3 = Segment(Node(3.0, -1.0), Node(3.0, 2.0), t=0.02)
        p4 = Segment(Node(3.0, 2.0), Node(2.0, 2.0), t=0.02)
        segs = box + [p1, p2, p3, p4]

        sec = MixedSection(segs)
        topo = sec.mixed_topology

        # Entire graph is 2-edge connected; 0 open edges!
        assert len(topo.open_edges) == 0
        assert len(topo.closed_edges) == 8
        assert topo.cycle_rank == 2
        assert topo.cell_count == 2

    def test_t06_pure_open_reduction(self) -> None:
        """T06: MixedSection on a pure open tree reduces to open section behavior."""
        # Simple I-section or channel
        s1 = Segment(Node(0.0, 0.0), Node(10.0, 0.0), t=2.0)
        s2 = Segment(Node(10.0, 0.0), Node(10.0, 20.0), t=2.0)
        s3 = Segment(Node(10.0, 20.0), Node(0.0, 20.0), t=2.0)
        segs = [s1, s2, s3]

        sec = MixedSection(segs)
        topo = sec.mixed_topology

        assert topo.is_open is True
        assert topo.is_closed is False
        assert topo.is_mixed is False
        assert topo.topology_type == "open"
        assert topo.cycle_rank == 0
        assert topo.cell_count == 0
        assert len(topo.closed_edges) == 0
        assert len(topo.open_edges) == 3

        # Hybrid J should equal J_open exactly
        assert sec.J_BB == 0.0
        assert math.isclose(sec.J_open, sec.J, rel_tol=1e-12)
        assert math.isclose(sec.J_total, sec.J, rel_tol=1e-12)

    def test_t06_pure_closed_reduction(self) -> None:
        """T06: MixedSection on a pure closed box reduces to closed section behavior."""
        box = make_box_segments(0.0, 0.0, 4.0, 2.0, 0.02)
        sec = MixedSection(box)
        topo = sec.mixed_topology

        assert topo.is_closed is True
        assert topo.is_open is False
        assert topo.is_mixed is False
        assert topo.topology_type == "closed"
        assert topo.cycle_rank == 1
        assert len(topo.open_edges) == 0

        # J_open should be 0.0, J_BB == J_total
        assert sec.J_open == 0.0
        assert math.isclose(sec.J_BB, sec.J, rel_tol=1e-12)
        assert math.isclose(sec.J_total, sec.J, rel_tol=1e-12)

    def test_t21_invalid_geometry_rejection(self) -> None:
        """T21: Disconnected, non-finite, and degenerate sections must raise explicit errors."""
        # Empty
        with pytest.raises(TopologyError, match="Section must contain at least one segment"):
            MixedSection([])

        # Disconnected
        s1 = Segment(Node(0.0, 0.0), Node(1.0, 0.0), t=1.0)
        s2 = Segment(Node(5.0, 5.0), Node(6.0, 5.0), t=1.0)
        with pytest.raises(TopologyError, match="Disconnected geometry detected"):
            MixedSection([s1, s2])

        # Non-positive thickness
        with pytest.raises(GeometryError, match="strictly positive"):
            MixedSection([Segment(Node(0.0, 0.0), Node(1.0, 0.0), t=0.0)])

        # Non-positive length
        with pytest.raises(GeometryError, match="strictly positive"):
            MixedSection([Segment(Node(0.0, 0.0), Node(0.0, 0.0), t=1.0)])


    def test_t22_immutability_and_caching(self) -> None:
        """T22: Section properties are cached and immutable."""
        box = make_box_segments(0.0, 0.0, 4.0, 2.0, 0.02)
        ant = Segment(Node(4.0, 2.0), Node(5.0, 2.0), t=0.01)
        sec = MixedSection(box + [ant])

        res1 = sec.torsion_properties()
        res2 = sec.torsion_properties()
        assert res1.J == res2.J
        assert res1.Cw == res2.Cw
        assert res1.shear_center == res2.shear_center

        # Verify translation produces valid new section
        sec_trans = sec.translated(10.0, 20.0)
        assert math.isclose(sec_trans.area, sec.area, rel_tol=1e-12)
        assert math.isclose(sec_trans.J, sec.J, rel_tol=1e-12)
        assert math.isclose(sec_trans.Cw, sec.Cw, rel_tol=1e-12)
