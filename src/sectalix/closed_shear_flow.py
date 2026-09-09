"""Closed-section transverse shear-flow analysis.

Implements the documented thin-wall closed-section formulation:
1. Virtual cut-open topology preserving 100% of physical material.
2. Particular basic shear flow q_b(s) via the open-tree shear-flow solver.
3. Compatibility vector b_c = sum_i B_ci \\int (q_{b,i}/t_i) ds.
4. Direct solve of H * q_0 = -b for redundant cell circulations q_0.
5. Assembly of continuous physical shear flow q_i(s) = q_{b,i}(s) + sum_c B_ci * q_{0,c}.
6. Verification of zero-twist cell compatibility: \\oint_c (q/t) ds = 0.
7. Verification of transverse resultant recovery: sum_i \\int q_i t_hat_i ds = V.
8. Exact centroidal torque T_z for shear-center computation.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import sys
from typing import TYPE_CHECKING, Sequence

import numpy as np

from sectalix.exceptions import GeometryError, SingularSectionError
from sectalix.primitives import Node, Segment
from sectalix.shear_flow import SegmentShearFlow, _compute_open_tree_flows
from sectalix.shear_load import ShearLoad

if TYPE_CHECKING:
    from sectalix.cells import CellTopology
    from sectalix.section import Section


@dataclass(frozen=True)
class ClosedSegmentShearFlow:
    """Exact shear-flow field along a single physical segment of a closed section.

    Supports exact evaluation q(xi) for xi in [0, 1] (s = xi * L).
    For tree segments, q(xi) is a single quadratic polynomial.
    For chord segments containing a virtual cut, q(xi) is evaluated via the two
    continuous analysis pieces shifted by the compatible redundant cell circulation.

    Attributes:
        segment: The underlying physical Segment.
        is_chord: True if this segment contained a virtual analysis cut.
        q0: Value at s=0 (xi=0).
        q_mid: Value at s=L/2 (xi=0.5).
        q_end: Value at s=L (xi=1.0).
        integral_q: Exact integral \\int_0^L q(s) ds.
        integral_q_over_t: Exact integral \\int_0^L (q(s)/t) ds.
        q0_circulation: Constant redundant cell circulation on this segment sum_c B_ci * q_{0,c}.
    """

    segment: Segment
    is_chord: bool
    integral_q: float
    integral_q_over_t: float
    q0_circulation: float
    # Internal polynomial representations:
    # If not is_chord: poly = (a0, a1, a2)
    # If is_chord: poly_a = (a0, a1, a2), poly_b = (a0, a1, a2), cut_param = xi_cut
    _poly_tree: tuple[float, float, float] | None = None
    _poly_a: tuple[float, float, float] | None = None
    _poly_b: tuple[float, float, float] | None = None
    _cut_param: float = 0.5

    @property
    def tangent(self) -> np.ndarray:
        """Unit tangent vector [tx, ty]^T pointing from p1 to p2."""
        dx = self.segment.p2.x - self.segment.p1.x
        dy = self.segment.p2.y - self.segment.p1.y
        length = self.segment.length
        return np.array([dx / length, dy / length], dtype=float)

    @property
    def length(self) -> float:
        """Segment length L."""
        return self.segment.length

    @property
    def q0(self) -> float:
        """Scalar shear flow at start (s=0, xi=0)."""
        return self.q_at_xi(0.0)

    @property
    def q_mid(self) -> float:
        """Scalar shear flow at midpoint (s=L/2, xi=0.5)."""
        return self.q_at_xi(0.5)

    @property
    def q_end(self) -> float:
        """Scalar shear flow at end (s=L, xi=1.0)."""
        return self.q_at_xi(1.0)

    def q_at_xi(self, xi: float) -> float:
        """Evaluate scalar shear flow q at normalized position xi in [0, 1]."""
        if not math.isfinite(xi):
            raise GeometryError(f"Normalized coordinate xi must be finite. Got xi={xi}.")
        if xi < -1e-7 or xi > 1.0 + 1e-7:
            raise GeometryError(
                f"Normalized coordinate xi must lie in [0, 1]. Got xi={xi}."
            )
        xi_clamped = max(0.0, min(1.0, xi))

        if not self.is_chord:
            assert self._poly_tree is not None
            a0, a1, a2 = self._poly_tree
            return a0 + a1 * xi_clamped + a2 * (xi_clamped**2)
        else:
            assert self._poly_a is not None and self._poly_b is not None
            xi_cut = self._cut_param
            if xi_clamped <= xi_cut:
                # Piece a: parameter in [0, 1] along piece a
                xi_a = xi_clamped / xi_cut if xi_cut > 0.0 else 0.0
                a0, a1, a2 = self._poly_a
                return (a0 + a1 * xi_a + a2 * (xi_a**2)) + self.q0_circulation
            else:
                # Piece b: parameter in [0, 1] along piece b
                xi_b = (xi_clamped - xi_cut) / (1.0 - xi_cut) if xi_cut < 1.0 else 1.0
                a0, a1, a2 = self._poly_b
                return (a0 + a1 * xi_b + a2 * (xi_b**2)) + self.q0_circulation

    def q_at_s(self, s: float) -> float:
        """Evaluate scalar shear flow q at physical distance s in [0, L]."""
        return self.q_at_xi(s / self.segment.length)


@dataclass(frozen=True)
class ClosedShearFlowResult:
    """Complete results of transverse shear-flow analysis for a closed section.

    Attributes:
        section: The analysed Section (or ClosedSection).
        load: Applied transverse ShearLoad (vx, vy).
        segment_flows: Tuple of ClosedSegmentShearFlow for each physical segment.
        redundant_cell_flows: Solved redundant cell circulations q_0 of shape (n_c,).
        b_vector: Basic compatibility vector b of shape (n_c,).
        recovered_resultant: Recovered resultant force [Vx, Vy]^T.
        compatibility_residuals: Vector of \\oint_c (q/t) ds for each cell.
        torque: Centroidal torque T_z of the final compatible flow.
        alpha: Bending gradient vector alpha from C * alpha = V.
        c_matrix: Centroidal coordinate second-moment matrix C.
    """

    section: Section
    load: ShearLoad
    segment_flows: tuple[ClosedSegmentShearFlow, ...]
    redundant_cell_flows: np.ndarray
    b_vector: np.ndarray
    recovered_resultant: np.ndarray
    compatibility_residuals: np.ndarray
    torque: float
    alpha: np.ndarray
    c_matrix: np.ndarray

    @property
    def resultant_error(self) -> float:
        """Euclidean norm of force recovery error: ||V_recovered - V||."""
        return float(np.linalg.norm(self.recovered_resultant - self.load.vector))

    @property
    def max_compatibility_residual(self) -> float:
        """Maximum absolute zero-twist compatibility residual across all cells."""
        return float(np.max(np.abs(self.compatibility_residuals)))


def _find_spanning_tree(
    v_count: int,
    canonical_edges: Sequence[tuple[int, int]],
    forced_edges: Sequence[int] | None = None,
) -> tuple[list[int], list[int]]:
    """Determine spanning tree and chord edges using Kruskal algorithm.

    Args:
        v_count: Number of canonical nodes.
        canonical_edges: List of (u, v) pairs for each segment.
        forced_edges: Optional list of edge indices to prioritize for the tree.

    Returns:
        (tree_edge_indices, chord_edge_indices)
    """
    parent = list(range(v_count))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i: int, j: int) -> bool:
        root_i = find(i)
        root_j = find(j)
        if root_i == root_j:
            return False
        parent[root_j] = root_i
        return True

    tree_edges: list[int] = []
    chord_edges: list[int] = []

    # If forced edges provided, add them first
    if forced_edges is not None:
        forced_set = set(forced_edges)
        for e_idx in forced_edges:
            u, v = canonical_edges[e_idx]
            if union(u, v):
                tree_edges.append(e_idx)
            else:
                chord_edges.append(e_idx)
        for e_idx in range(len(canonical_edges)):
            if e_idx not in forced_set:
                u, v = canonical_edges[e_idx]
                if union(u, v):
                    tree_edges.append(e_idx)
                else:
                    chord_edges.append(e_idx)
    else:
        for e_idx, (u, v) in enumerate(canonical_edges):
            if union(u, v):
                tree_edges.append(e_idx)
            else:
                chord_edges.append(e_idx)

    return tree_edges, chord_edges


def calculate_closed_shear_flow(
    section: Section,
    vx: float | ShearLoad,
    vy: float | None = None,
    cut_param: float = 0.5,
    safety_factor: float = 1e4,
    spanning_tree_edges: Sequence[int] | None = None,
) -> ClosedShearFlowResult:
    """Compute exact closed-section transverse shear flow under load (Vx, Vy).

    Args:
        section: Validated ClosedSection (or Section with cellular topology).
        vx: Transverse shear force in x direction, or ShearLoad object.
        vy: Transverse shear force in y direction.
        cut_param: Parameter xi_cut in (0, 1) where chords are virtually cut. Default 0.5.
        safety_factor: Multiplier for numerical positive-definiteness check on C and H.
        spanning_tree_edges: Optional sequence of segment indices to use as spanning tree.

    Returns:
        ClosedShearFlowResult containing exact physical segment flows and verification.
    """
    from sectalix.cells import CellTopology, extract_cell_topology

    if not math.isfinite(cut_param) or cut_param <= 0.01 or cut_param >= 0.99:
        raise GeometryError(
            f"cut_param must be strictly between 0.01 and 0.99. Got {cut_param}."
        )
    if not math.isfinite(safety_factor) or safety_factor <= 0.0:
        raise GeometryError(
            f"safety_factor must be finite and strictly positive. Got {safety_factor}."
        )

    # 1. Parse load
    if isinstance(vx, ShearLoad):
        load = vx
    else:
        if vy is None:
            raise GeometryError("vy must be provided when vx is a scalar force.")
        load = ShearLoad(vx=float(vx), vy=float(vy))

    # 2. Section geometry and planar cell topology
    if hasattr(section, "cell_topology") and isinstance(
        section.cell_topology, CellTopology
    ):
        cell_topo: CellTopology = section.cell_topology
    else:
        cell_topo = extract_cell_topology(
            section.segments,
            node_tolerance=section._node_tolerance,
            safety_factor=safety_factor,
        )

    canonical_nodes = cell_topo.canonical_nodes
    canonical_edges = cell_topo.canonical_edges
    B = cell_topo.B
    n_c = cell_topo.cell_count
    e_count = len(section.segments)
    v_count = len(canonical_nodes)

    # 3. Construct centroidal coordinate second-moment matrix C
    ix = section.Ix
    iy = section.Iy
    ixy = section.Ixy
    c_matrix = np.array([[iy, ixy], [ixy, ix]], dtype=float)

    eigvals = np.linalg.eigvalsh(c_matrix)
    lambda_min = float(eigvals[0])
    lambda_max = float(eigvals[1])
    eps_mach = float(np.finfo(float).eps)
    c_threshold = safety_factor * eps_mach * lambda_max

    if lambda_min <= c_threshold or lambda_min <= 0.0:
        raise SingularSectionError(
            f"Section inertia matrix C is singular or rank-deficient: "
            f"lambda_min={lambda_min:.3e}, threshold={c_threshold:.3e}."
        )

    # Solve C * alpha = V directly
    alpha = np.linalg.solve(c_matrix, load.vector)
    if not np.all(np.isfinite(alpha)):
        raise GeometryError(
            f"Non-finite solution vector alpha in closed shear-flow solve: alpha={alpha}."
        )

    # 4. Spanning tree and chord cuts
    tree_indices, chord_indices = _find_spanning_tree(
        v_count=v_count,
        canonical_edges=canonical_edges,
        forced_edges=spanning_tree_edges,
    )
    assert len(chord_indices) == n_c, f"Expected {n_c} chords, got {len(chord_indices)}"

    # Build the virtual cut-open tree segments and topology
    # Physical material is preserved 100% exactly.
    cut_segments: list[Segment] = []
    cut_canonical_nodes: list[Node] = list(canonical_nodes)
    cut_canonical_edges: list[tuple[int, int]] = []

    # Map tree segments directly
    tree_to_cut_map: dict[int, int] = {}
    for s_idx in tree_indices:
        tree_to_cut_map[s_idx] = len(cut_segments)
        cut_segments.append(section.segments[s_idx])
        cut_canonical_edges.append(canonical_edges[s_idx])

    # For each chord edge, split into piece a (u -> cut_node1) and piece b (cut_node2 -> v)
    chord_to_cut_pieces: dict[int, tuple[int, int]] = {}
    for k_chord, s_idx in enumerate(chord_indices):
        seg = section.segments[s_idx]
        u, v = canonical_edges[s_idx]
        p1 = seg.p1
        p2 = seg.p2
        t = seg.t

        x_cut = p1.x + cut_param * (p2.x - p1.x)
        y_cut = p1.y + cut_param * (p2.y - p1.y)

        cut_node1 = Node(x=x_cut, y=y_cut)
        cut_node2 = Node(x=x_cut, y=y_cut)

        idx_cut1 = len(cut_canonical_nodes)
        cut_canonical_nodes.append(cut_node1)
        idx_cut2 = len(cut_canonical_nodes)
        cut_canonical_nodes.append(cut_node2)

        # Piece a: from u to cut_node1
        seg_a = Segment(p1=p1, p2=cut_node1, t=t)
        piece_a_idx = len(cut_segments)
        cut_segments.append(seg_a)
        cut_canonical_edges.append((u, idx_cut1))

        # Piece b: from cut_node2 to v
        seg_b = Segment(p1=cut_node2, p2=p2, t=t)
        piece_b_idx = len(cut_segments)
        cut_segments.append(seg_b)
        cut_canonical_edges.append((idx_cut2, v))

        chord_to_cut_pieces[s_idx] = (piece_a_idx, piece_b_idx)

    # 5. Solve basic cut-open shear flow q_b using exact mechanics
    ref_x = min(n.x for n in canonical_nodes)
    ref_y = min(n.y for n in canonical_nodes)

    cut_flows = _compute_open_tree_flows(
        segments=cut_segments,
        canonical_nodes=cut_canonical_nodes,
        canonical_edges=cut_canonical_edges,
        total_area=section.area,
        alpha=alpha,
        ref_origin=(ref_x, ref_y),
    )

    # 6. Integrate basic shear flow on each original physical segment
    # I_qb[i] = \int_0^{L_i} q_{b,i}(s) ds
    # For tree segments: L * (a0 + a1/2 + a2/3)
    # For chord segments: sum of integrals of piece a and piece b
    I_qb = np.zeros(e_count, dtype=float)

    for s_idx in tree_indices:
        flow = cut_flows[tree_to_cut_map[s_idx]]
        L = flow.length
        I_qb[s_idx] = L * (flow.a0 + 0.5 * flow.a1 + (1.0 / 3.0) * flow.a2)

    for s_idx in chord_indices:
        idx_a, idx_b = chord_to_cut_pieces[s_idx]
        flow_a = cut_flows[idx_a]
        flow_b = cut_flows[idx_b]
        int_a = flow_a.length * (flow_a.a0 + 0.5 * flow_a.a1 + (1.0 / 3.0) * flow_a.a2)
        int_b = flow_b.length * (flow_b.a0 + 0.5 * flow_b.a1 + (1.0 / 3.0) * flow_b.a2)
        I_qb[s_idx] = int_a + int_b

    # 7. Assemble dimensionless basic-flow compatibility vector b_scaled: shape (n_c,)
    # w_i = (\int q_b ds) / t_i decomposed using frexp to avoid intermediate float overflow
    m_w: list[float] = []
    e_w: list[int] = []
    for i, seg in enumerate(section.segments):
        iq = float(I_qb[i])
        if iq == 0.0:
            m_w.append(0.0)
            e_w.append(0)
        else:
            m_iq, e_iq = math.frexp(iq)
            m_t, e_t = math.frexp(seg.t)
            m_div, e_div = math.frexp(m_iq / m_t)
            m_w.append(m_div)
            e_w.append((e_iq - e_t) + e_div)

    nonzeros_w = [e_w[i] for i in range(e_count) if m_w[i] != 0.0]
    if not nonzeros_w:
        e_max_w = 0
        w_scaled = np.zeros(e_count, dtype=float)
    else:
        e_max_w = max(nonzeros_w)
        w_scaled = np.zeros(e_count, dtype=float)
        for i in range(e_count):
            if m_w[i] != 0.0:
                delta_e = e_w[i] - e_max_w
                if delta_e >= -1100:
                    w_scaled[i] = math.ldexp(m_w[i], delta_e)

    # Dimensionless compatibility vector b_scaled is strictly finite (entries in [-E, E])
    b_scaled = B @ w_scaled

    # 8. Solve compatibility system H_scaled @ q0_scaled = -b_scaled
    H_scaled = cell_topo.H_scaled
    e_max_H = cell_topo.e_max_H

    q0_scaled = np.linalg.solve(H_scaled, -b_scaled)
    if not np.all(np.isfinite(q0_scaled)):
        raise GeometryError(
            f"Non-finite scaled redundant cell flow vector q0_scaled encountered: {q0_scaled}."
        )

    # Reconstruct physical redundant cell circulation q_0 = ldexp(q0_scaled, e_max_w - e_max_H)
    delta_q = e_max_w - e_max_H
    q_0 = np.zeros(n_c, dtype=float)
    for c in range(n_c):
        try:
            q_0[c] = math.ldexp(float(q0_scaled[c]), delta_q)
        except OverflowError:
            q_0[c] = math.copysign(math.inf, q0_scaled[c])

    # Construct physical b vector: ldexp(b_scaled, e_max_w)
    b_vec = np.zeros(n_c, dtype=float)
    for c in range(n_c):
        val = float(b_scaled[c])
        if val != 0.0:
            try:
                b_vec[c] = math.ldexp(val, e_max_w)
            except OverflowError:
                b_vec[c] = math.copysign(math.inf, val)

    if not np.all(np.isfinite(b_vec)):
        raise OverflowError(
            "Physical compatibility vector b or compatibility residuals overflow IEEE-754 float64 "
            "under the applied shear load. Use dimensionless formulation or reduce load magnitude."
        )

    # 9. Segment constant circulations: q_{0,seg} = B^T * q_0
    q0_seg = B.T @ q_0

    # 10. Assemble physical ClosedSegmentShearFlow for every segment
    physical_segment_flows: list[ClosedSegmentShearFlow] = []

    for s_idx in range(e_count):
        seg = section.segments[s_idx]
        L = seg.length
        t = seg.t
        circ = float(q0_seg[s_idx])
        int_total = float(I_qb[s_idx] + circ * L)
        try:
            int_over_t = float(int_total / t)
        except OverflowError:
            int_over_t = math.copysign(math.inf, int_total)
        if not math.isfinite(int_over_t):
            raise OverflowError(
                "Segment transverse-shear compatibility integral overflows IEEE-754 float64. "
                "Use dimensionless formulation or reduce load magnitude."
            )

        if s_idx in tree_indices:
            flow = cut_flows[tree_to_cut_map[s_idx]]
            # Shift a0 by circ: q(xi) = (a0 + circ) + a1 * xi + a2 * xi^2
            poly_tree = (flow.a0 + circ, flow.a1, flow.a2)
            seg_flow = ClosedSegmentShearFlow(
                segment=seg,
                is_chord=False,
                integral_q=int_total,
                integral_q_over_t=int_over_t,
                q0_circulation=circ,
                _poly_tree=poly_tree,
            )
        else:
            idx_a, idx_b = chord_to_cut_pieces[s_idx]
            flow_a = cut_flows[idx_a]
            flow_b = cut_flows[idx_b]
            poly_a = (flow_a.a0, flow_a.a1, flow_a.a2)
            poly_b = (flow_b.a0, flow_b.a1, flow_b.a2)
            seg_flow = ClosedSegmentShearFlow(
                segment=seg,
                is_chord=True,
                integral_q=int_total,
                integral_q_over_t=int_over_t,
                q0_circulation=circ,
                _poly_a=poly_a,
                _poly_b=poly_b,
                _cut_param=cut_param,
            )
        physical_segment_flows.append(seg_flow)

    # 11. Verification: Resultant recovery
    # V_rec = sum_i \int q_i * \hat{t}_i ds
    v_rec = np.zeros(2, dtype=float)
    for flow in physical_segment_flows:
        v_rec += flow.integral_q * flow.tangent

    if not np.all(np.isfinite(v_rec)):
        raise GeometryError(
            f"Non-finite recovered resultant in closed shear-flow: {v_rec}."
        )

    # 12. Verification: Zero-twist compatibility
    # res_c = sum_i B_ci * \int (q_i / t_i) ds = b_c + (H * q_0)_c
    # In scaled form: res_scaled = b_scaled + H_scaled @ q0_scaled to eliminate inf - inf cancellation
    res_scaled = b_scaled + H_scaled @ q0_scaled
    compat_residuals = np.zeros(n_c, dtype=float)
    for c in range(n_c):
        val = float(res_scaled[c])
        if val == 0.0:
            compat_residuals[c] = 0.0
        else:
            try:
                compat_residuals[c] = math.ldexp(val, e_max_w)
            except OverflowError:
                compat_residuals[c] = math.copysign(math.inf, val)

    # 13. Centroidal torque T_z using local reference shift to avoid precision loss under large translations
    # For each segment: lever arm h_i = x_{c,1} * t_y - y_{c,1} * t_x is constant along straight segment.
    # T_z = sum_i h_i * \int q_i ds
    cx_loc = sum(
        seg.area * ((seg.p1.x - ref_x) + (seg.p2.x - ref_x)) / 2.0
        for seg in section.segments
    ) / section.area
    cy_loc = sum(
        seg.area * ((seg.p1.y - ref_y) + (seg.p2.y - ref_y)) / 2.0
        for seg in section.segments
    ) / section.area

    torque_terms: list[float] = []
    for flow in physical_segment_flows:
        seg = flow.segment
        xc1 = (seg.p1.x - ref_x) - cx_loc
        yc1 = (seg.p1.y - ref_y) - cy_loc
        tx = flow.tangent[0]
        ty = flow.tangent[1]
        lever_arm = xc1 * ty - yc1 * tx
        torque_terms.append(lever_arm * flow.integral_q)
    torque_z = math.fsum(torque_terms)

    # 14. Guard against silent float64 overflow in physical quantities
    if not np.all(np.isfinite(b_vec)) or not np.all(
        np.isfinite(compat_residuals)
    ):
        raise OverflowError(
            "Physical compatibility vector b or compatibility residuals overflow IEEE-754 float64 "
            "under the applied shear load. Use dimensionless formulation or reduce load magnitude."
        )
    for seg_flow in physical_segment_flows:
        if not math.isfinite(seg_flow.integral_q_over_t):
            raise OverflowError(
                "Segment transverse-shear compatibility integral overflows IEEE-754 float64."
            )

    return ClosedShearFlowResult(
        section=section,
        load=load,
        segment_flows=tuple(physical_segment_flows),
        redundant_cell_flows=q_0,
        b_vector=b_vec,
        recovered_resultant=v_rec,
        compatibility_residuals=compat_residuals,
        torque=float(torque_z),
        alpha=alpha,
        c_matrix=c_matrix,
    )
