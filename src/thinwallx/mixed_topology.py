"""Planar graph decomposition and topology extraction for mixed open-closed sections (ThinWallX v0.6).

Implements deterministic graph decomposition:
1. Bridge detection via Tarjan DFS: partitions edges into cyclic edges E_c and open bridges E_o.
2. Connected cyclic components G_c,j and biconnected cyclic blocks.
3. Generalized Euler-Poincaré cycle-rank: n_c = sum(|E_c,j| - |V_c,j| + 1) = |E_c| - |V_c| + k_c.
4. Planar half-edge cyclic face-walk on each cyclic component to extract bounded CCW cells.
5. Topological role identification: junction_nodes, articulation_nodes, free_tip_nodes.
6. Dimensionless cell flexibility matrix H_scaled = B diag(rho_scaled) B^T and positive-definiteness check.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
import math
from typing import Sequence

import numpy as np

from thinwallx.cells import Cell
from thinwallx.exceptions import GeometryError, SingularSectionError, TopologyError
from thinwallx.primitives import Node, Segment
from thinwallx.validation import (
    _validate_node_tolerance,
    _validate_segment_intersections,
    cluster_nodes,
)


@dataclass(frozen=True)
class MixedTopology:
    """Container for the complete decomposed topology of a mixed open-closed section.

    Attributes:
        canonical_nodes: Tuple of clustered representative Node objects.
        canonical_edges: Tuple of (u, v) canonical node index pairs for each segment.
        closed_edges: Tuple of segment indices belonging to at least one cycle (E_c).
        open_edges: Tuple of segment indices that are open bridges / branches (E_o).
        cells: Tuple of extracted bounded counter-clockwise Cell objects.
        cell_incidence: Oriented cell-wall incidence matrix of shape (n_c, E).
        cell_areas: Tuple of enclosed median-line areas for each cell.
        cyclic_components: Tuple of edge-index tuples for each connected component of G_c.
        cyclic_blocks: Tuple of edge-index tuples for each biconnected block of G_c.
        open_components: Tuple of edge-index tuples for each connected tree component of E_o.
        junction_nodes: Canonical node indices incident to both closed and open edges.
        articulation_nodes: Canonical node indices that are cut vertices of the full graph.
        free_tip_nodes: Canonical node indices of degree 1 in the full graph.
        cycle_rank: Total cycle rank n_c = |E_c| - |V_c| + k_c.
        cut_edges: Selected chord segment indices strictly from E_c (length n_c).
        H_scaled: Dimensionless cell flexibility matrix of shape (n_c, n_c).
        e_max_H: Binary scaling exponent of the flexibility matrix.
        block_cell_indices: Global cell indices for each independent cyclic block.
        H_block_scaled: Independently scaled flexibility submatrix for each cyclic block.
        e_max_H_block: Binary scaling exponent paired with each H_block_scaled matrix.
        topology_type: Classification string: 'open', 'closed', or 'mixed'.
    """

    canonical_nodes: tuple[Node, ...]
    canonical_edges: tuple[tuple[int, int], ...]
    closed_edges: tuple[int, ...]
    open_edges: tuple[int, ...]
    cells: tuple[Cell, ...]
    cell_incidence: np.ndarray
    cell_areas: tuple[float, ...]
    cyclic_components: tuple[tuple[int, ...], ...]
    cyclic_blocks: tuple[tuple[int, ...], ...]
    open_components: tuple[tuple[int, ...], ...]
    junction_nodes: tuple[int, ...]
    articulation_nodes: tuple[int, ...]
    free_tip_nodes: tuple[int, ...]
    cycle_rank: int
    cut_edges: tuple[int, ...]
    H_scaled: np.ndarray
    e_max_H: int = 0
    block_cell_indices: tuple[tuple[int, ...], ...] = ()
    H_block_scaled: tuple[np.ndarray, ...] = ()
    e_max_H_block: tuple[int, ...] = ()
    topology_type: str = "mixed"
    _H: np.ndarray | None = None

    @property
    def B(self) -> np.ndarray:
        """Oriented cell-wall incidence matrix of shape (n_c, E)."""
        return self.cell_incidence

    @property
    def H(self) -> np.ndarray:
        """Physical cell flexibility matrix H of shape (n_c, n_c).

        Raises:
            OverflowError: If H elements exceed IEEE-754 float64 range.
        """
        if self._H is not None:
            return self._H
        if self.cycle_rank == 0:
            return np.zeros((0, 0), dtype=np.float64)
        try:
            H_physical = np.zeros_like(self.H_scaled, dtype=np.float64)
            for cell_ids, H_block, e_max_block in zip(
                self.block_cell_indices,
                self.H_block_scaled,
                self.e_max_H_block,
            ):
                for local_r, global_r in enumerate(cell_ids):
                    for local_c, global_c in enumerate(cell_ids):
                        val = float(H_block[local_r, local_c])
                        if val == 0.0:
                            H_physical[global_r, global_c] = 0.0
                        else:
                            m, e = math.frexp(val)
                            H_physical[global_r, global_c] = math.ldexp(
                                m, e + e_max_block
                            )
            return H_physical
        except OverflowError:
            raise OverflowError(
                "Cell flexibility matrix H overflows IEEE-754 float64. "
                "Use mixed_topology.H_scaled and e_max_H for dimensionless analysis."
            )

    @property
    def cell_count(self) -> int:
        """Number of bounded cells n_c."""
        return len(self.cells)

    @property
    def segment_count(self) -> int:
        """Total number of physical segments E."""
        return len(self.canonical_edges)

    @property
    def node_count(self) -> int:
        """Total number of canonical nodes V."""
        return len(self.canonical_nodes)

    @property
    def is_mixed(self) -> bool:
        """True if section contains both closed cells and open branches."""
        return len(self.closed_edges) > 0 and len(self.open_edges) > 0

    @property
    def is_closed(self) -> bool:
        """True if section contains closed cells and zero open branches."""
        return len(self.closed_edges) > 0 and len(self.open_edges) == 0

    @property
    def is_open(self) -> bool:
        """True if section contains zero closed cells."""
        return len(self.closed_edges) == 0


def extract_mixed_topology(
    segments: Sequence[Segment],
    node_tolerance: float = 1e-9,
    safety_factor: float = 1e4,
) -> MixedTopology:
    """Extract, validate, and decompose planar topology for an arbitrary connected section.

    Partitions edges into cyclic cells E_c and open bridge branches E_o, verifies generalized
    Euler-Poincaré cycle rank across cyclic components, identifies junction nodes, and
    constructs dimensionless flexibility matrix H_scaled.

    Args:
        segments: Sequence of physical Segment objects.
        node_tolerance: Absolute Euclidean coordinate tolerance for node clustering.
        safety_factor: Multiplier on machine epsilon for positive-definiteness of H.

    Returns:
        MixedTopology instance containing complete decomposition.

    Raises:
        GeometryError: For non-finite geometry, collapsed cells, or self-intersecting segments.
        TopologyError: For empty sections, disconnected geometry, or cycle-rank mismatch.
        SingularSectionError: If cell flexibility matrix H is singular or rank-deficient.
    """
    _validate_node_tolerance(node_tolerance)

    if len(segments) == 0:
        raise TopologyError("Section must contain at least one segment.")

    if not math.isfinite(safety_factor) or safety_factor <= 0.0:
        raise GeometryError(
            f"safety_factor must be finite and strictly positive. Got {safety_factor}."
        )

    # 1. Validate individual segments
    all_nodes: list[Node] = []
    for i, seg in enumerate(segments):
        if not math.isfinite(seg.t) or seg.t <= 0.0:
            raise GeometryError(
                f"Segment {i} has invalid thickness t={seg.t}. Thickness must be finite and > 0."
            )
        if not math.isfinite(seg.length) or seg.length <= 0.0:
            raise GeometryError(
                f"Segment {i} has invalid length (L={seg.length:.6e}); length must be finite and > 0."
            )
        all_nodes.append(seg.p1)
        all_nodes.append(seg.p2)

    # 2. Cluster endpoints into canonical nodes
    canonical_nodes, mapping = cluster_nodes(all_nodes, tol=node_tolerance)
    v_count = len(canonical_nodes)
    e_count = len(segments)

    # 3. Build edge list in canonical node indices
    canonical_edges: list[tuple[int, int]] = []
    seen_edges: dict[tuple[int, int], int] = {}
    for seg_idx in range(e_count):
        u = mapping[2 * seg_idx]
        v = mapping[2 * seg_idx + 1]
        if u == v:
            raise GeometryError(
                f"Segment {seg_idx} has coincident endpoints within tolerance {node_tolerance:.1e}."
            )
        edge_key = (min(u, v), max(u, v))
        if edge_key in seen_edges:
            prior_seg = seen_edges[edge_key]
            raise TopologyError(
                f"Duplicate segment detected: segment {seg_idx} connects canonical nodes "
                f"({u}, {v}), which is already connected by segment {prior_seg}."
            )
        seen_edges[edge_key] = seg_idx
        canonical_edges.append((u, v))

    # Validate non-endpoint segment intersections and collinear overlaps
    _validate_segment_intersections(segments, canonical_edges, node_tolerance)

    # 4. Check global connectivity
    adj_graph: list[list[tuple[int, int]]] = [[] for _ in range(v_count)]
    for s_idx, (u, v) in enumerate(canonical_edges):
        adj_graph[u].append((v, s_idx))
        adj_graph[v].append((u, s_idx))

    visited_conn: set[int] = set()
    queue_conn: deque[int] = deque([0])
    visited_conn.add(0)
    while queue_conn:
        curr = queue_conn.popleft()
        for nbr, _ in adj_graph[curr]:
            if nbr not in visited_conn:
                visited_conn.add(nbr)
                queue_conn.append(nbr)

    if len(visited_conn) < v_count:
        raise TopologyError(
            f"Disconnected geometry detected: visited {len(visited_conn)} of {v_count} nodes."
        )

    # 5. Bridge detection & Articulation nodes via Tarjan DFS
    tin = [-1] * v_count
    low = [-1] * v_count
    timer = 0
    is_bridge = [False] * e_count
    articulation_set: set[int] = set()

    def dfs_bridges(u: int, p_edge: int = -1) -> None:
        nonlocal timer
        tin[u] = low[u] = timer
        timer += 1
        children = 0
        for v, e_idx in adj_graph[u]:
            if e_idx == p_edge:
                continue
            if tin[v] != -1:
                # Back-edge
                if tin[v] < low[u]:
                    low[u] = tin[v]
            else:
                # Tree-edge
                children += 1
                dfs_bridges(v, e_idx)
                if low[v] < low[u]:
                    low[u] = low[v]
                if low[v] > tin[u]:
                    is_bridge[e_idx] = True
                if p_edge != -1 and low[v] >= tin[u]:
                    articulation_set.add(u)
        if p_edge == -1 and children > 1:
            articulation_set.add(u)

    dfs_bridges(0, -1)

    closed_edges_list: list[int] = []
    open_edges_list: list[int] = []
    for s_idx in range(e_count):
        if is_bridge[s_idx]:
            open_edges_list.append(s_idx)
        else:
            closed_edges_list.append(s_idx)

    closed_edges = tuple(closed_edges_list)
    open_edges = tuple(open_edges_list)

    # Free tip nodes: degree 1 in full graph
    free_tips: list[int] = [u for u in range(v_count) if len(adj_graph[u]) == 1]
    free_tip_nodes = tuple(free_tips)
    articulation_nodes = tuple(sorted(articulation_set))

    # Open components: connected components formed by open edges (forest of bridges)
    open_components_list: list[tuple[int, ...]] = []
    if open_edges:
        open_adj: dict[int, list[tuple[int, int]]] = {}
        for s_idx in open_edges:
            u, v = canonical_edges[s_idx]
            open_adj.setdefault(u, []).append((v, s_idx))
            open_adj.setdefault(v, []).append((u, s_idx))

        visited_open_edges: set[int] = set()
        for s_idx in open_edges:
            if s_idx in visited_open_edges:
                continue
            comp_edges: list[int] = []
            queue_oe: deque[int] = deque([s_idx])
            visited_open_edges.add(s_idx)
            while queue_oe:
                curr_s = queue_oe.popleft()
                comp_edges.append(curr_s)
                u_curr, v_curr = canonical_edges[curr_s]
                for node in (u_curr, v_curr):
                    for nbr, next_s in open_adj.get(node, []):
                        if next_s not in visited_open_edges:
                            visited_open_edges.add(next_s)
                            queue_oe.append(next_s)
            open_components_list.append(tuple(sorted(comp_edges)))

    open_components = tuple(open_components_list)

    # Junction nodes: incident to at least one closed edge and at least one open edge
    closed_incident_nodes: set[int] = set()
    for s_idx in closed_edges:
        u, v = canonical_edges[s_idx]
        closed_incident_nodes.add(u)
        closed_incident_nodes.add(v)

    open_incident_nodes: set[int] = set()
    for s_idx in open_edges:
        u, v = canonical_edges[s_idx]
        open_incident_nodes.add(u)
        open_incident_nodes.add(v)

    junction_nodes = tuple(sorted(closed_incident_nodes & open_incident_nodes))

    # 6. Cyclic components & Euler-Poincaré cycle rank
    # G_c = (V_c, E_c)
    bounded_cells: list[Cell] = []
    cyclic_components_list: list[tuple[int, ...]] = []
    cyclic_blocks_list: list[tuple[int, ...]] = []
    total_cycle_rank = 0

    if not closed_edges:
        # Pure open section
        k_c = 0
        total_cycle_rank = 0
        cell_incidence = np.zeros((0, e_count), dtype=float)
        cell_areas = ()
        cut_edges = ()
        H_scaled = np.zeros((0, 0), dtype=float)
        e_max_H = 0
        block_cell_indices_tuple: tuple[tuple[int, ...], ...] = ()
        H_block_scaled_tuple: tuple[np.ndarray, ...] = ()
        e_max_H_block_tuple: tuple[int, ...] = ()
        _H = np.zeros((0, 0), dtype=float)
        topology_type = "open"
    else:
        # Build adjacency for G_c
        closed_adj: dict[int, list[tuple[int, int]]] = {node: [] for node in closed_incident_nodes}
        for s_idx in closed_edges:
            u, v = canonical_edges[s_idx]
            closed_adj[u].append((v, s_idx))
            closed_adj[v].append((u, s_idx))

        # Find connected components of G_c
        visited_closed_nodes: set[int] = set()
        for start_node in sorted(closed_incident_nodes):
            if start_node in visited_closed_nodes:
                continue
            comp_nodes: list[int] = []
            comp_edges: list[int] = []
            queue_cc = deque([start_node])
            visited_closed_nodes.add(start_node)
            while queue_cc:
                curr = queue_cc.popleft()
                comp_nodes.append(curr)
                for nbr, s_idx in closed_adj[curr]:
                    if s_idx not in comp_edges:
                        comp_edges.append(s_idx)
                    if nbr not in visited_closed_nodes:
                        visited_closed_nodes.add(nbr)
                        queue_cc.append(nbr)

            comp_edges_tuple = tuple(sorted(comp_edges))
            cyclic_components_list.append(comp_edges_tuple)

            # Component cycle rank: n_c,j = |E_c,j| - |V_c,j| + 1
            v_cj = len(comp_nodes)
            e_cj = len(comp_edges)
            n_cj = e_cj - v_cj + 1
            if n_cj < 1:
                raise TopologyError(
                    f"Cyclic component has invalid cycle rank ({n_cj} < 1)."
                )
            total_cycle_rank += n_cj

            # Extract biconnected blocks for this cyclic component
            # Using standard DFS edge-stack
            tin_b: dict[int, int] = {}
            low_b: dict[int, int] = {}
            timer_b = 0
            edge_stack: list[int] = []

            def dfs_blocks(u: int, p_edge: int = -1) -> None:
                nonlocal timer_b
                tin_b[u] = low_b[u] = timer_b
                timer_b += 1
                for v, s_idx in closed_adj[u]:
                    if s_idx == p_edge:
                        continue
                    if s_idx not in edge_stack:
                        edge_stack.append(s_idx)
                    if v in tin_b:
                        if tin_b[v] < low_b[u]:
                            low_b[u] = tin_b[v]
                    else:
                        dfs_blocks(v, s_idx)
                        if low_b[v] < low_b[u]:
                            low_b[u] = low_b[v]
                        if low_b[v] >= tin_b[u]:
                            # Pop block edges
                            block: list[int] = []
                            while True:
                                popped = edge_stack.pop()
                                block.append(popped)
                                if popped == s_idx:
                                    break
                            cyclic_blocks_list.append(tuple(sorted(block)))

            dfs_blocks(start_node, -1)

            # Extract planar faces for this cyclic component
            # Outgoing half-edges for component nodes
            outgoing_he: dict[int, list[tuple[int, int, int, int]]] = {u: [] for u in comp_nodes}
            for s_idx in comp_edges:
                u, v = canonical_edges[s_idx]
                outgoing_he[u].append((u, v, s_idx, +1))
                outgoing_he[v].append((v, u, s_idx, -1))

            # Sort counter-clockwise by polar angle
            for u in comp_nodes:
                u_pt = canonical_nodes[u]

                def polar_angle(he: tuple[int, int, int, int]) -> float:
                    tgt = canonical_nodes[he[1]]
                    return math.atan2(tgt.y - u_pt.y, tgt.x - u_pt.x)

                outgoing_he[u].sort(key=polar_angle)

            # Cyclic successor map
            next_he: dict[tuple[int, int, int], tuple[int, int, int, int]] = {}
            for u in comp_nodes:
                nbrs = outgoing_he[u]
                k = len(nbrs)
                for idx, he in enumerate(nbrs):
                    prev_he = nbrs[(idx - 1) % k]
                    incoming_key = (he[1], u, he[2])
                    next_he[incoming_key] = prev_he

            visited_he: set[tuple[int, int, int]] = set()
            raw_faces: list[list[tuple[int, int, int, int]]] = []

            for u in comp_nodes:
                for he in outgoing_he[u]:
                    he_key = (he[0], he[1], he[2])
                    if he_key in visited_he:
                        continue
                    face_walk: list[tuple[int, int, int, int]] = []
                    curr_he = he
                    while True:
                        curr_key = (curr_he[0], curr_he[1], curr_he[2])
                        if curr_key in visited_he:
                            break
                        visited_he.add(curr_key)
                        face_walk.append(curr_he)
                        inc_key = (curr_he[0], curr_he[1], curr_he[2])
                        curr_he = next_he[inc_key]

                    if face_walk:
                        raw_faces.append(face_walk)

            # Classify faces by shoelace area
            comp_bounded_cells: list[Cell] = []
            for face in raw_faces:
                poly_pts = [canonical_nodes[he[0]] for he in face]
                num_pts = len(poly_pts)
                x0, y0 = poly_pts[0].x, poly_pts[0].y
                signed_area = 0.5 * sum(
                    (poly_pts[k].x - x0) * (poly_pts[(k + 1) % num_pts].y - y0)
                    - (poly_pts[(k + 1) % num_pts].x - x0) * (poly_pts[k].y - y0)
                    for k in range(num_pts)
                )

                if signed_area == 0.0:
                    raise GeometryError(
                        "Degenerate zero-area face detected during planar cell extraction."
                    )

                if signed_area > 0.0:
                    # CCW bounded cell
                    c_id = len(bounded_cells) + len(comp_bounded_cells)
                    node_indices = tuple(he[0] for he in face)
                    seg_indices = tuple(he[2] for he in face)
                    seg_orientations = tuple(he[3] for he in face)
                    perim = sum(segments[he[2]].length for he in face)

                    # Scale-aware collapsed check: |A| / perim^2 <= 1e-12
                    m_A, e_A = math.frexp(signed_area)
                    m_p, e_p = math.frexp(perim)
                    rel_m = m_A / (m_p * m_p)
                    rel_exp = e_A - 2 * e_p
                    try:
                        is_collapsed = (rel_exp < -45) or (math.ldexp(rel_m, rel_exp) <= 1e-12)
                    except OverflowError:
                        is_collapsed = False

                    if is_collapsed:
                        raise GeometryError(
                            f"Cell {c_id} has collapsed/zero area (A={signed_area:.6e}, perimeter={perim:.6e})."
                        )

                    cell_obj = Cell(
                        id=c_id,
                        node_indices=node_indices,
                        segment_indices=seg_indices,
                        segment_orientations=seg_orientations,
                        area=float(signed_area),
                        perimeter=float(perim),
                    )
                    comp_bounded_cells.append(cell_obj)

            if len(comp_bounded_cells) != n_cj:
                raise TopologyError(
                    f"Extracted cell count ({len(comp_bounded_cells)}) does not match cycle rank "
                    f"({n_cj} = E - V + 1) for cyclic component."
                )

            bounded_cells.extend(comp_bounded_cells)

        k_c = len(cyclic_components_list)
        # Verify generalized Euler-Poincaré: n_c = |E_c| - |V_c| + k_c
        expected_n_c = len(closed_edges) - len(closed_incident_nodes) + k_c
        if total_cycle_rank != expected_n_c or len(bounded_cells) != expected_n_c:
            raise TopologyError(
                f"Global cycle rank mismatch: extracted {len(bounded_cells)} cells, "
                f"expected {expected_n_c} (|E_c|={len(closed_edges)}, |V_c|={len(closed_incident_nodes)}, k_c={k_c})."
            )

        # Assemble cell-wall incidence matrix B: shape (n_c, E)
        n_c = len(bounded_cells)
        cell_incidence = np.zeros((n_c, e_count), dtype=float)
        for cell in bounded_cells:
            for s_idx, ori in zip(cell.segment_indices, cell.segment_orientations):
                cell_incidence[cell.id, s_idx] = float(ori)

        # Shared wall condition verification
        for s_idx in closed_edges:
            incident_cells = np.where(cell_incidence[:, s_idx] != 0.0)[0]
            if len(incident_cells) == 2:
                c1, c2 = incident_cells
                if cell_incidence[c1, s_idx] * cell_incidence[c2, s_idx] != -1.0:
                    raise TopologyError(
                        f"Shared wall {s_idx} does not have opposite incidence signs in cells {c1} and {c2}."
                    )
            elif len(incident_cells) > 2:
                raise TopologyError(
                    f"Segment {s_idx} belongs to {len(incident_cells)} cells (maximum 2 allowed)."
                )

        # Dimensionless cell flexibility matrix H_scaled = B @ diag(rho_scaled) @ B^T
        m_rhos = []
        e_rhos = []
        for s_idx in range(e_count):
            seg = segments[s_idx]
            mL, eL = math.frexp(seg.length)
            mt, et = math.frexp(seg.t)
            m_r, e_norm = math.frexp(mL / mt)
            m_rhos.append(m_r)
            e_rhos.append((eL - et) + e_norm)

        closed_e_rhos = [e_rhos[i] for i in closed_edges]
        e_max_H = max(closed_e_rhos)
        rho_scaled = np.zeros(e_count, dtype=float)
        for i in closed_edges:
            delta_e = e_rhos[i] - e_max_H
            if delta_e >= -1100:
                rho_scaled[i] = math.ldexp(m_rhos[i], delta_e)

        H_scaled = cell_incidence @ np.diag(rho_scaled) @ cell_incidence.T

        # Positive-definiteness check analyzed independently per cyclic block / component
        # This prevents heterogeneous scales between decoupled cells (e.g. barbell sections)
        # from falsely triggering SingularSectionError.
        eps_mach = float(np.finfo(float).eps)

        # Map cells to blocks and retain each block's independent numerical scale.
        # Independent biconnected cyclic blocks have exactly zero H coupling, so a
        # single global exponent would erase a sufficiently smaller block.
        block_cell_indices_list: list[tuple[int, ...]] = []
        H_block_scaled_list: list[np.ndarray] = []
        e_max_H_block_list: list[int] = []
        assigned_cells: set[int] = set()
        for block_idx, block in enumerate(cyclic_blocks_list):
            block_set = set(block)
            block_cell_indices = tuple(
                c.id for c in bounded_cells
                if set(c.segment_indices).issubset(block_set)
            )
            if not block_cell_indices:
                continue
            duplicate_cells = assigned_cells.intersection(block_cell_indices)
            if duplicate_cells:
                raise TopologyError(
                    "A bounded cell was assigned to more than one independent cyclic block: "
                    f"{sorted(duplicate_cells)}."
                )
            assigned_cells.update(block_cell_indices)

            block_e_rhos = [e_rhos[i] for i in block]
            e_max_block = max(block_e_rhos)
            rho_block = np.zeros(e_count, dtype=float)
            for i in block:
                delta_e = e_rhos[i] - e_max_block
                if delta_e >= -1100:
                    rho_block[i] = math.ldexp(m_rhos[i], delta_e)

            B_block = cell_incidence[block_cell_indices, :]
            H_block = B_block @ np.diag(rho_block) @ B_block.T
            eigvals_b = np.linalg.eigvalsh(H_block)
            l_min_b = float(eigvals_b[0])
            l_max_b = float(eigvals_b[-1])
            thresh_b = safety_factor * eps_mach * l_max_b

            if l_min_b <= thresh_b or l_min_b <= 0.0:
                raise SingularSectionError(
                    f"Cyclic block {block_idx} flexibility matrix H is singular or rank-deficient: "
                    f"lambda_min={l_min_b:.6e}, threshold={thresh_b:.6e}."
                )

            block_cell_indices_list.append(block_cell_indices)
            H_block_scaled_list.append(H_block)
            e_max_H_block_list.append(e_max_block)

        if assigned_cells != set(range(n_c)):
            missing = sorted(set(range(n_c)) - assigned_cells)
            raise TopologyError(
                "Not every bounded cell was assigned to an independent cyclic block; "
                f"missing cell indices: {missing}."
            )

        block_cell_indices_tuple = tuple(block_cell_indices_list)
        H_block_scaled_tuple = tuple(H_block_scaled_list)
        e_max_H_block_tuple = tuple(e_max_H_block_list)

        # Reconstruct physical H if strictly within IEEE-754 float64 range
        _H = None
        try:
            H_physical = np.zeros_like(H_scaled, dtype=np.float64)
            for cell_ids, H_block, e_max_block in zip(
                block_cell_indices_tuple,
                H_block_scaled_tuple,
                e_max_H_block_tuple,
            ):
                for local_r, global_r in enumerate(cell_ids):
                    for local_c, global_c in enumerate(cell_ids):
                        val = float(H_block[local_r, local_c])
                        if val == 0.0:
                            H_physical[global_r, global_c] = 0.0
                        else:
                            m, e = math.frexp(val)
                            H_physical[global_r, global_c] = math.ldexp(
                                m, e + e_max_block
                            )
            _H = H_physical
        except OverflowError:
            _H = None

        cell_areas = tuple(c.area for c in bounded_cells)
        topology_type = "mixed" if open_edges else "closed"

    # Select deterministic cut chords strictly from E_c (none from E_o)
    # Spanning forest of E_o is built first, then extended into E_c
    parent_uf = list(range(v_count))

    def find_uf(i: int) -> int:
        path = []
        while parent_uf[i] != i:
            path.append(i)
            i = parent_uf[i]
        for node in path:
            parent_uf[node] = i
        return i

    def union_uf(i: int, j: int) -> bool:
        root_i = find_uf(i)
        root_j = find_uf(j)
        if root_i == root_j:
            return False
        parent_uf[root_i] = root_j
        return True

    # 1. Add all open bridge edges to spanning tree (never forms cycles)
    for s_idx in open_edges:
        u, v = canonical_edges[s_idx]
        union_uf(u, v)

    # 2. Add closed edges to complete spanning tree; chords are non-tree closed edges
    tree_closed: list[int] = []
    chord_closed: list[int] = []
    for s_idx in closed_edges:
        u, v = canonical_edges[s_idx]
        if union_uf(u, v):
            tree_closed.append(s_idx)
        else:
            chord_closed.append(s_idx)

    cut_edges = tuple(chord_closed)
    if len(cut_edges) != total_cycle_rank:
        raise TopologyError(
            f"Number of chord cuts ({len(cut_edges)}) does not match cycle rank ({total_cycle_rank})."
        )

    return MixedTopology(
        canonical_nodes=tuple(canonical_nodes),
        canonical_edges=tuple(canonical_edges),
        closed_edges=closed_edges,
        open_edges=open_edges,
        cells=tuple(bounded_cells),
        cell_incidence=cell_incidence,
        cell_areas=cell_areas,
        cyclic_components=tuple(cyclic_components_list),
        cyclic_blocks=tuple(cyclic_blocks_list),
        open_components=open_components,
        junction_nodes=junction_nodes,
        articulation_nodes=articulation_nodes,
        free_tip_nodes=free_tip_nodes,
        cycle_rank=total_cycle_rank,
        cut_edges=cut_edges,
        H_scaled=H_scaled,
        e_max_H=e_max_H,
        block_cell_indices=block_cell_indices_tuple,
        H_block_scaled=H_block_scaled_tuple,
        e_max_H_block=e_max_H_block_tuple,
        topology_type=topology_type,
        _H=_H,
    )
