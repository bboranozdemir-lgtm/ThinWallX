"""Mixed open-closed transverse shear-flow analysis.

Implements the documented thin-wall mixed open-closed formulation:
1. Spanning tree containing all open bridges E_o; exactly n_c chords chosen strictly from E_c.
2. Zero virtual cuts on open branches; virtual cuts placed only within closed chord walls.
3. Particular basic shear flow q_b(s) via the open-tree shear-flow solver with q_b=0 at real free tips.
4. Dimensionless compatibility vector b_c = sum_i B_ci \\int (q_{b,i}/t_i) ds via frexp/ldexp scaling.
5. Direct solve of H_scaled * q0_scaled = -b_scaled for closed cell circulations.
6. Assembly of continuous physical shear flow q_e(s) = q_{b,e}(s) + (B^T q_0)_e.
7. Verification of Kirchhoff balance at all nodes (including junctions), zero-twist compatibility,
   and exact transverse resultant recovery.
8. Centroidal torque T_z for mixed shear-center calculation.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import TYPE_CHECKING, Sequence

import numpy as np

from sectalix.exceptions import GeometryError, SingularSectionError, TopologyError
from sectalix.primitives import Node, Segment
from sectalix.shear_flow import _compute_open_tree_flows
from sectalix.shear_load import ShearLoad

if TYPE_CHECKING:
    from sectalix.mixed_topology import MixedTopology
    from sectalix.section import Section


@dataclass(frozen=True)
class MixedSegmentShearFlow:
    """Exact shear-flow field along a single physical segment of a mixed section.

    Attributes:
        segment: The underlying physical Segment.
        is_chord: True if this segment contained a virtual analysis cut.
        is_open: True if this segment is an open bridge / branch (in E_o).
        integral_q: Exact integral \\int_0^L q(s) ds.
        integral_q_over_t: Exact integral \\int_0^L (q(s)/t) ds.
        q0_circulation: Constant redundant cell circulation on this segment (B^T q_0)_e.
    """

    segment: Segment
    is_chord: bool
    is_open: bool
    integral_q: float
    integral_q_over_t: float
    q0_circulation: float
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
                xi_a = xi_clamped / xi_cut if xi_cut > 0.0 else 0.0
                a0, a1, a2 = self._poly_a
                return (a0 + a1 * xi_a + a2 * (xi_a**2)) + self.q0_circulation
            else:
                xi_b = (xi_clamped - xi_cut) / (1.0 - xi_cut) if xi_cut < 1.0 else 1.0
                a0, a1, a2 = self._poly_b
                return (a0 + a1 * xi_b + a2 * (xi_b**2)) + self.q0_circulation

    def q_at_s(self, s: float) -> float:
        """Evaluate scalar shear flow q at physical distance s in [0, L]."""
        return self.q_at_xi(s / self.segment.length)


@dataclass(frozen=True)
class MixedShearFlowResult:
    """Complete results of transverse shear-flow analysis for a mixed section.

    Attributes:
        section: The analysed Section or MixedSection.
        load: Applied transverse ShearLoad (vx, vy).
        segment_flows: Tuple of MixedSegmentShearFlow for each physical segment.
        redundant_cell_flows: Solved redundant cell circulations q_0 of shape (n_c,).
        b_vector: Basic compatibility vector b of shape (n_c,).
        recovered_resultant: Recovered resultant force [Vx, Vy]^T.
        compatibility_residuals: Vector of \\oint_c (q/t) ds for each cell.
        torque: Centroidal torque T_z of the final compatible flow.
        alpha: Bending gradient vector alpha from C * alpha = V.
        c_matrix: Centroidal coordinate second-moment matrix C.
        relative_residual: Dimensionless solver backward error eta.
        condition_number_estimate: Condition number estimate of H_scaled.
        max_node_residual: Maximum Kirchhoff node equilibrium residual.
    """

    section: Section
    load: ShearLoad
    segment_flows: tuple[MixedSegmentShearFlow, ...]
    redundant_cell_flows: np.ndarray
    b_vector: np.ndarray
    recovered_resultant: np.ndarray
    compatibility_residuals: np.ndarray
    torque: float
    alpha: np.ndarray
    c_matrix: np.ndarray
    relative_residual: float = 0.0
    condition_number_estimate: float = 1.0
    max_node_residual: float = 0.0


def calculate_mixed_shear_flow(
    section: Section,
    vx: float | ShearLoad,
    vy: float | None = None,
    cut_param: float = 0.5,
    safety_factor: float = 1e4,
    spanning_tree_edges: Sequence[int] | None = None,
) -> MixedShearFlowResult:
    """Compute exact transverse shear flow for an arbitrary connected mixed open-closed section.

    Args:
        section: Validated Section or MixedSection.
        vx: ShearLoad object or transverse shear force in x.
        vy: Transverse shear force in y (required if vx is float).
        cut_param: Parameter xi_cut in (0, 1) for chord cuts (default 0.5).
        safety_factor: Multiplier for numerical positive-definiteness check.
        spanning_tree_edges: Optional sequence of segment indices to use as spanning tree.

    Returns:
        MixedShearFlowResult containing exact segment shear-flow fields.

    Raises:
        GeometryError: For non-finite inputs, invalid cut parameters, or singular inertia matrix.
        TopologyError: If spanning tree cuts an open bridge or does not span all nodes.
        OverflowError: If compatibility quantities overflow IEEE-754 binary64 range.
    """
    if isinstance(vx, ShearLoad):
        load = vx
    else:
        if vy is None:
            raise GeometryError("vy must be provided when vx is a scalar force.")
        load = ShearLoad(vx=float(vx), vy=float(vy))

    if not math.isfinite(cut_param) or cut_param <= 0.0 or cut_param >= 1.0:
        raise GeometryError(
            f"cut_param must lie strictly inside (0, 1). Got {cut_param}."
        )

    from sectalix.mixed_topology import MixedTopology, extract_mixed_topology

    if hasattr(section, "mixed_topology") and isinstance(
        section.mixed_topology, MixedTopology
    ):
        mixed_topo: MixedTopology = section.mixed_topology
    else:
        mixed_topo = extract_mixed_topology(
            section.segments,
            node_tolerance=getattr(section, "_node_tolerance", 1e-9),
            safety_factor=safety_factor,
        )

    canonical_nodes = mixed_topo.canonical_nodes
    canonical_edges = mixed_topo.canonical_edges
    v_count = len(canonical_nodes)
    e_count = len(section.segments)
    n_c = mixed_topo.cycle_rank
    B = mixed_topo.B

    # 1. Inertia tensor and bending gradient solve
    ix = section.Ix
    iy = section.Iy
    ixy = section.Ixy

    c_mat = np.array([[iy, ixy], [ixy, ix]], dtype=float)
    eigvals_c = np.linalg.eigvalsh(c_mat)
    min_eig_c = float(eigvals_c[0])
    max_eig_c = float(eigvals_c[-1])
    eps_mach = float(np.finfo(float).eps)
    c_threshold = safety_factor * eps_mach * max_eig_c

    if min_eig_c <= c_threshold or min_eig_c <= 0.0:
        raise SingularSectionError(
            f"Section second-moment matrix C is singular or degenerate (min_eig={min_eig_c:.6e}, "
            f"threshold={c_threshold:.6e}). Cannot solve transverse shear-flow equations."
        )

    v_vector = np.array([load.vx, load.vy], dtype=float)
    alpha = np.linalg.solve(c_mat, v_vector)
    if not np.all(np.isfinite(alpha)):
        raise OverflowError(
            "Bending gradient vector alpha overflows IEEE-754 float64 under the applied shear load. "
            "Use dimensionless formulation or reduce load magnitude."
        )

    # 2. Spanning tree and chord partition
    if spanning_tree_edges is not None:
        if len(spanning_tree_edges) != v_count - 1:
            raise TopologyError(
                "Provided spanning_tree_edges does not form a valid connected spanning tree."
            )
        tree_set = set(spanning_tree_edges)
        if len(tree_set) != v_count - 1:
            raise TopologyError(
                "Provided spanning_tree_edges does not form a valid connected spanning tree."
            )
        for s_idx in mixed_topo.open_edges:
            if s_idx not in tree_set:
                raise TopologyError(
                    "Provided spanning_tree_edges does not form a valid connected spanning tree."
                )

        parent = list(range(v_count))

        def find(i: int) -> int:
            path = []
            while parent[i] != i:
                path.append(i)
                i = parent[i]
            for node in path:
                parent[node] = i
            return i

        for e_idx in tree_set:
            if e_idx < 0 or e_idx >= e_count:
                raise TopologyError(
                    "Provided spanning_tree_edges does not form a valid connected spanning tree."
                )
            u, v = canonical_edges[e_idx]
            root_u = find(u)
            root_v = find(v)
            if root_u == root_v:
                raise TopologyError(
                    "Provided spanning_tree_edges does not form a valid connected spanning tree."
                )
            parent[root_u] = root_v

        root_0 = find(0)
        for i in range(1, v_count):
            if find(i) != root_0:
                raise TopologyError(
                    "Provided spanning_tree_edges does not form a valid connected spanning tree."
                )

        chord_indices = [i for i in range(e_count) if i not in tree_set]
        tree_indices = [i for i in range(e_count) if i in tree_set]
        for ch in chord_indices:
            if ch not in mixed_topo.closed_edges:
                raise TopologyError(
                    f"Chord segment {ch} is not in E_c. Only closed-cell chords may be virtually cut."
                )
    else:
        chord_indices = list(mixed_topo.cut_edges)
        chord_set = set(chord_indices)
        tree_indices = [i for i in range(e_count) if i not in chord_set]

    # Handle pure open section boundary case (n_c == 0)
    if n_c == 0:
        ref_x = min(n.x for n in canonical_nodes)
        ref_y = min(n.y for n in canonical_nodes)

        with np.errstate(over="raise"):
            try:
                open_flows = _compute_open_tree_flows(
                    segments=section.segments,
                    canonical_nodes=canonical_nodes,
                    canonical_edges=canonical_edges,
                    total_area=section.area,
                    alpha=alpha,
                    ref_origin=(ref_x, ref_y),
                )

                cx_loc = sum(
                    s.area * ((s.p1.x - ref_x) + (s.p2.x - ref_x)) / 2.0
                    for s in section.segments
                ) / section.area
                cy_loc = sum(
                    s.area * ((s.p1.y - ref_y) + (s.p2.y - ref_y)) / 2.0
                    for s in section.segments
                ) / section.area

                segment_flows: list[MixedSegmentShearFlow] = []
                v_rec = np.zeros(2, dtype=float)
                tz_total = 0.0

                for s_idx in range(e_count):
                    seg = section.segments[s_idx]
                    flow = open_flows[s_idx]
                    L = seg.length
                    t = seg.t
                    int_q = float(L * (flow.a0 + 0.5 * flow.a1 + (1.0 / 3.0) * flow.a2))
                    int_q_over_t = float(int_q / t)

                    if not (math.isfinite(int_q) and math.isfinite(int_q_over_t)):
                        raise OverflowError(
                            f"Shear flow integral overflows IEEE-754 float64 on segment {s_idx}."
                        )

                    seg_flow = MixedSegmentShearFlow(
                        segment=seg,
                        is_chord=False,
                        is_open=True,
                        integral_q=int_q,
                        integral_q_over_t=int_q_over_t,
                        q0_circulation=0.0,
                        _poly_tree=(flow.a0, flow.a1, flow.a2),
                    )
                    segment_flows.append(seg_flow)
                    v_rec += int_q * seg_flow.tangent

                    # Local torque contribution anchored to ref_origin
                    x1_c = (seg.p1.x - ref_x) - cx_loc
                    y1_c = (seg.p1.y - ref_y) - cy_loc
                    tx = (seg.p2.x - seg.p1.x) / L
                    ty = (seg.p2.y - seg.p1.y) / L
                    pc = x1_c * ty - y1_c * tx
                    tz_total += pc * int_q
            except (FloatingPointError, OverflowError):
                raise OverflowError(
                    "Pure open shear flow resultant, torque, or flow integrals overflow IEEE-754 float64."
                )

        if not (np.all(np.isfinite(v_rec)) and math.isfinite(tz_total)):
            raise OverflowError(
                "Recovered resultant or torque overflows IEEE-754 float64 under applied shear load."
            )

        # Verify Kirchhoff node balance across all nodes
        node_res = np.zeros(v_count, dtype=float)
        node_flow_sum = np.zeros(v_count, dtype=float)
        for s_idx, (u, v) in enumerate(canonical_edges):
            sf = segment_flows[s_idx]
            node_res[u] -= sf.q0
            node_res[v] += sf.q_end
            node_flow_sum[u] += abs(sf.q0)
            node_flow_sum[v] += abs(sf.q_end)
        max_node_residual = float(np.max(np.abs(node_res)))
        flow_scale = max(
            max(
                max(abs(sf.q0), abs(sf.q_end)) for sf in segment_flows
            ) if segment_flows else 0.0,
            float(np.linalg.norm(v_vector)),
        )
        for node_idx in range(v_count):
            violates_balance = (
                abs(node_res[node_idx]) > 0.0
                if flow_scale == 0.0
                else abs(node_res[node_idx]) / flow_scale > 1e-10
            )
            if violates_balance:
                raise GeometryError(
                    f"Kirchhoff shear flow balance violated at node {node_idx}: "
                    f"residual={node_res[node_idx]:.6e}, scale={flow_scale:.6e}."
                )

        return MixedShearFlowResult(
            section=section,
            load=load,
            segment_flows=tuple(segment_flows),
            redundant_cell_flows=np.zeros(0, dtype=float),
            b_vector=np.zeros(0, dtype=float),
            recovered_resultant=v_rec,
            compatibility_residuals=np.zeros(0, dtype=float),
            torque=float(tz_total),
            alpha=alpha,
            c_matrix=c_mat,
            relative_residual=0.0,
            condition_number_estimate=1.0,
            max_node_residual=max_node_residual,
        )

    # 3. Construct cut-open tree graph by splitting chord edges
    cut_canonical_nodes = list(canonical_nodes)
    cut_segments: list[Segment] = []
    cut_canonical_edges: list[tuple[int, int]] = []

    tree_to_cut_map: dict[int, int] = {}
    for s_idx in tree_indices:
        tree_to_cut_map[s_idx] = len(cut_segments)
        cut_segments.append(section.segments[s_idx])
        cut_canonical_edges.append(canonical_edges[s_idx])

    chord_to_cut_pieces: dict[int, tuple[int, int]] = {}
    for s_idx in chord_indices:
        seg = section.segments[s_idx]
        u, v = canonical_edges[s_idx]
        p1 = seg.p1
        p2 = seg.p2
        t = seg.t

        # Interior cut point
        xc = p1.x + cut_param * (p2.x - p1.x)
        yc = p1.y + cut_param * (p2.y - p1.y)

        idx_cut1 = len(cut_canonical_nodes)
        cut_node1 = Node(x=xc, y=yc, id=idx_cut1)
        cut_canonical_nodes.append(cut_node1)

        idx_cut2 = len(cut_canonical_nodes)
        cut_node2 = Node(x=xc, y=yc, id=idx_cut2)
        cut_canonical_nodes.append(cut_node2)

        # Piece a: u -> cut_node1
        seg_a = Segment(p1=p1, p2=cut_node1, t=t)
        piece_a_idx = len(cut_segments)
        cut_segments.append(seg_a)
        cut_canonical_edges.append((u, idx_cut1))

        # Piece b: cut_node2 -> v
        seg_b = Segment(p1=cut_node2, p2=p2, t=t)
        piece_b_idx = len(cut_segments)
        cut_segments.append(seg_b)
        cut_canonical_edges.append((idx_cut2, v))

        chord_to_cut_pieces[s_idx] = (piece_a_idx, piece_b_idx)

    # 4. Solve basic cut-open shear flow q_b using exact mechanics
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

    # 5. Integrate basic shear flow on each original physical segment
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

    # 6. Decompose basic compatibility integrals before any physical division.
    # w_i = (\int q_b ds) / t_i decomposed using frexp to avoid intermediate float overflow
    m_w: list[float] = []
    e_w: list[int] = []
    for i in range(e_count):
        seg = section.segments[i]
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

    # 7. Solve H_block_scaled @ q0_block_scaled = -b_block_scaled for each
    # independent cyclic block. A barbell's uncoupled cells must never share an
    # exponent merely because an open bridge connects their parent graph.
    q_0 = np.zeros(n_c, dtype=float)
    b_vec = np.zeros(n_c, dtype=float)
    compat_residuals = np.zeros(n_c, dtype=float)
    q0_seg = np.zeros(e_count, dtype=float)
    relative_residual = 0.0
    cond_h = 1.0

    for cell_ids, H_block, e_max_H_block in zip(
        mixed_topo.block_cell_indices,
        mixed_topo.H_block_scaled,
        mixed_topo.e_max_H_block,
    ):
        block_indices = np.array(cell_ids, dtype=int)
        B_block = B[block_indices, :]
        block_edges = np.flatnonzero(np.any(B_block != 0.0, axis=0))
        nonzero_exponents = [e_w[i] for i in block_edges if m_w[i] != 0.0]
        e_max_w_block = max(nonzero_exponents) if nonzero_exponents else 0
        w_scaled_block = np.zeros(e_count, dtype=float)
        for edge_idx in block_edges:
            if m_w[edge_idx] == 0.0:
                continue
            exponent_shift = e_w[edge_idx] - e_max_w_block
            if exponent_shift >= -1100:
                w_scaled_block[edge_idx] = math.ldexp(
                    m_w[edge_idx], exponent_shift
                )

        b_scaled_block = B_block @ w_scaled_block
        q0_scaled_block = np.linalg.solve(H_block, -b_scaled_block)
        if not np.all(np.isfinite(q0_scaled_block)):
            raise GeometryError(
                "Non-finite scaled redundant cell flow vector encountered in an "
                f"independent cyclic block: {q0_scaled_block}."
            )

        scaled_residual = H_block @ q0_scaled_block + b_scaled_block
        norm_res = float(np.linalg.norm(scaled_residual, ord=np.inf))
        norm_H = float(np.linalg.norm(H_block, ord=np.inf))
        norm_q0 = float(np.linalg.norm(q0_scaled_block, ord=np.inf))
        norm_b = float(np.linalg.norm(b_scaled_block, ord=np.inf))
        denominator = norm_H * norm_q0 + norm_b
        block_relative_residual = norm_res / denominator if denominator > 0.0 else 0.0
        relative_residual = max(relative_residual, block_relative_residual)

        eigvals_h = np.linalg.eigvalsh(H_block)
        block_condition = (
            float(eigvals_h[-1] / eigvals_h[0])
            if eigvals_h[0] > 0.0
            else math.inf
        )
        cond_h = max(cond_h, block_condition)

        q_exponent = e_max_w_block - e_max_H_block
        for local_idx, global_idx in enumerate(cell_ids):
            scaled_q = float(q0_scaled_block[local_idx])
            scaled_b = float(b_scaled_block[local_idx])
            scaled_r = float(scaled_residual[local_idx])
            try:
                q_0[global_idx] = math.ldexp(scaled_q, q_exponent)
                b_vec[global_idx] = math.ldexp(scaled_b, e_max_w_block)
                compat_residuals[global_idx] = math.ldexp(
                    scaled_r, e_max_w_block
                )
            except OverflowError as exc:
                raise OverflowError(
                    "Physical redundant flow, compatibility vector, or compatibility "
                    "residual overflows IEEE-754 float64 under applied load."
                ) from exc

        q0_segment_scaled = B_block.T @ q0_scaled_block
        for edge_idx in block_edges:
            scaled_value = float(q0_segment_scaled[edge_idx])
            try:
                q0_seg[edge_idx] = math.ldexp(scaled_value, q_exponent)
            except OverflowError as exc:
                raise OverflowError(
                    "Segment redundant circulation overflows IEEE-754 float64 under applied load."
                ) from exc

    if not (
        np.all(np.isfinite(q_0))
        and np.all(np.isfinite(b_vec))
        and np.all(np.isfinite(compat_residuals))
        and np.all(np.isfinite(q0_seg))
    ):
        raise OverflowError(
            "Physical compatibility quantities overflow IEEE-754 float64 under applied load."
        )

    # 8. Assemble physical MixedSegmentShearFlow for every segment
    physical_segment_flows: list[MixedSegmentShearFlow] = []
    for s_idx in range(e_count):
        seg = section.segments[s_idx]
        L = seg.length
        t = seg.t
        is_op = s_idx in mixed_topo.open_edges
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
            poly_tree = (flow.a0 + circ, flow.a1, flow.a2)
            seg_flow = MixedSegmentShearFlow(
                segment=seg,
                is_chord=False,
                is_open=is_op,
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
            seg_flow = MixedSegmentShearFlow(
                segment=seg,
                is_chord=True,
                is_open=is_op,
                integral_q=int_total,
                integral_q_over_t=int_over_t,
                q0_circulation=circ,
                _poly_a=poly_a,
                _poly_b=poly_b,
                _cut_param=cut_param,
            )
        physical_segment_flows.append(seg_flow)

    # 9. Resultant recovery verification
    v_rec = np.zeros(2, dtype=float)
    for flow in physical_segment_flows:
        v_rec += flow.integral_q * flow.tangent

    diff_v = np.linalg.norm(v_rec - v_vector)
    resultant_scale = float(np.linalg.norm(v_vector))
    violates_resultant = (
        diff_v > 0.0
        if resultant_scale == 0.0
        else diff_v / resultant_scale > 1e-10
    )
    if violates_resultant:
        raise GeometryError(
            f"Transverse resultant recovery error: recovered {v_rec}, applied {v_vector}, diff={diff_v:.6e}."
        )

    # 10. Kirchhoff balance at all nodes
    node_res = np.zeros(v_count, dtype=float)
    node_flow_sum = np.zeros(v_count, dtype=float)
    for s_idx, (u, v) in enumerate(canonical_edges):
        sf = physical_segment_flows[s_idx]
        node_res[u] -= sf.q0
        node_res[v] += sf.q_end
        node_flow_sum[u] += abs(sf.q0)
        node_flow_sum[v] += abs(sf.q_end)
    max_node_residual = float(np.max(np.abs(node_res)))
    flow_scale = max(
        max(
            max(abs(sf.q0), abs(sf.q_end)) for sf in physical_segment_flows
        ) if physical_segment_flows else 0.0,
        float(np.linalg.norm(v_vector)),
    )
    for node_idx in range(v_count):
        violates_balance = (
            abs(node_res[node_idx]) > 0.0
            if flow_scale == 0.0
            else abs(node_res[node_idx]) / flow_scale > 1e-10
        )
        if violates_balance:
            raise GeometryError(
                f"Kirchhoff shear flow balance violated at node {node_idx}: "
                f"residual={node_res[node_idx]:.6e}, scale={flow_scale:.6e}."
            )

    # 11. Centroidal torque T_z using local reference shift to eliminate cancellation under large translations
    cx_loc = sum(
        seg.area * ((seg.p1.x - ref_x) + (seg.p2.x - ref_x)) / 2.0
        for seg in section.segments
    ) / section.area
    cy_loc = sum(
        seg.area * ((seg.p1.y - ref_y) + (seg.p2.y - ref_y)) / 2.0
        for seg in section.segments
    ) / section.area

    tz_total = 0.0
    for flow in physical_segment_flows:
        seg = flow.segment
        L = seg.length
        tx = (seg.p2.x - seg.p1.x) / L
        ty = (seg.p2.y - seg.p1.y) / L
        x1_c = (seg.p1.x - ref_x) - cx_loc
        y1_c = (seg.p1.y - ref_y) - cy_loc
        pc = x1_c * ty - y1_c * tx
        tz_total += pc * flow.integral_q

    if not math.isfinite(tz_total):
        raise OverflowError("Computed torque overflows IEEE-754 float64 under applied shear load.")

    return MixedShearFlowResult(
        section=section,
        load=load,
        segment_flows=tuple(physical_segment_flows),
        redundant_cell_flows=q_0,
        b_vector=b_vec,
        recovered_resultant=v_rec,
        compatibility_residuals=compat_residuals,
        torque=float(tz_total),
        alpha=alpha,
        c_matrix=c_mat,
        relative_residual=relative_residual,
        condition_number_estimate=cond_h,
        max_node_residual=max_node_residual,
    )
