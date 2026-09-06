"""Planar cell topology extraction and validation for closed sections (Sectalix v0.5).

Implements deterministic planar face-walk extraction for straight-segment networks:
- Node clustering and canonical representation
- Planar half-edge cyclic ordering by polar angle
- Bounded counter-clockwise face extraction
- Verification of Euler cycle rank: n_c = E - V + 1
- Explicit rejection of mixed open/closed topology (dangling branches)
- Exact polygon shoelace cell area
- Oriented cell-wall incidence matrix B in {-1, 0, +1}^(n_c x E)
- Cell flexibility matrix H = B diag(L/t) B^T and scale-aware positive-definiteness check
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
import math
from typing import Sequence

import numpy as np

from sectalix.exceptions import GeometryError, SingularSectionError, TopologyError
from sectalix.primitives import Node, Segment
from sectalix.validation import (
    _validate_node_tolerance,
    _validate_segment_intersections,
    cluster_nodes,
)


@dataclass(frozen=True)
class Cell:
    """A single bounded counter-clockwise cell in a closed thin-walled section.

    Attributes:
        id: Zero-based cell index.
        node_indices: Ordered tuple of canonical node indices defining the cell polygon (CCW).
        segment_indices: Ordered tuple of physical segment indices forming the cell boundary.
        segment_orientations: Tuple of +1 / -1 indicating if segment i is traversed
            from p1 to p2 (+1) or p2 to p1 (-1) along the CCW cell boundary.
        area: Exact median-line enclosed area A_c > 0.
        perimeter: Total perimeter length L_c = sum(L_i for i in cell).
    """

    id: int
    node_indices: tuple[int, ...]
    segment_indices: tuple[int, ...]
    segment_orientations: tuple[int, ...]
    area: float
    perimeter: float

    def __post_init__(self) -> None:
        if not math.isfinite(self.area) or self.area <= 0.0:
            raise GeometryError(
                f"Cell {self.id} must have a finite, strictly positive area. Got {self.area}."
            )
        if not math.isfinite(self.perimeter) or self.perimeter <= 0.0:
            raise GeometryError(
                f"Cell {self.id} must have a finite, strictly positive perimeter. Got {self.perimeter}."
            )


@dataclass(frozen=True)
class CellTopology:
    """Container for the complete extracted planar cell topology of a closed section.

    Attributes:
        canonical_nodes: Tuple of clustered representative Node objects.
        canonical_edges: Tuple of (u, v) canonical node index pairs for each segment.
        cells: Tuple of extracted CCW bounded Cell objects.
        B: Cell-wall incidence matrix of shape (n_c, E), values in {-1, 0, +1}.
        H_scaled: Dimensionless cell flexibility matrix of shape (n_c, n_c).
        e_max_H: Binary scaling exponent of the flexibility matrix.
    """

    canonical_nodes: tuple[Node, ...]
    canonical_edges: tuple[tuple[int, int], ...]
    cells: tuple[Cell, ...]
    B: np.ndarray
    H_scaled: np.ndarray
    e_max_H: int = 0
    _H: np.ndarray | None = None

    @property
    def H(self) -> np.ndarray:
        """Cell flexibility matrix H = B @ diag(L/t) @ B^T of shape (n_c, n_c).

        Raises:
            OverflowError: If H elements exceed IEEE-754 float64 range.
        """
        if self._H is not None:
            return self._H
        try:
            H_physical = np.empty_like(self.H_scaled, dtype=np.float64)
            for r in range(self.H_scaled.shape[0]):
                for c in range(self.H_scaled.shape[1]):
                    val = float(self.H_scaled[r, c])
                    if val == 0.0:
                        H_physical[r, c] = 0.0
                    else:
                        m, e = math.frexp(val)
                        H_physical[r, c] = math.ldexp(m, e + self.e_max_H)
            return H_physical
        except OverflowError:
            raise OverflowError(
                "Cell flexibility matrix H overflows IEEE-754 float64. "
                "Use cell_topology.H_scaled and e_max_H for dimensionless analysis."
            )

    @property
    def cell_count(self) -> int:
        """Number of bounded cells n_c."""
        return len(self.cells)

    @property
    def segment_count(self) -> int:
        """Number of physical segments E."""
        return len(self.canonical_edges)

    @property
    def node_count(self) -> int:
        """Number of canonical nodes V."""
        return len(self.canonical_nodes)


def extract_cell_topology(
    segments: Sequence[Segment],
    node_tolerance: float = 1e-9,
    safety_factor: float = 1e4,
) -> CellTopology:
    """Extract, validate, and construct planar cell topology for a closed cellular section.

    Args:
        segments: Sequence of physical Segment objects.
        node_tolerance: Absolute Euclidean coordinate tolerance for node clustering.
        safety_factor: Multiplier on machine epsilon for positive-definiteness of H.

    Returns:
        CellTopology instance containing canonical nodes, cells, B, and H.

    Raises:
        GeometryError: For non-finite geometry, collapsed cells, or mixed open/closed topology.
        TopologyError: For empty sections, disconnected geometry, or invalid cell counts.
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

    # 4. Check connectivity
    adj_graph: list[list[int]] = [[] for _ in range(v_count)]
    for u, v in canonical_edges:
        adj_graph[u].append(v)
        adj_graph[v].append(u)

    visited_conn: set[int] = set()
    queue_conn: deque[int] = deque([0])
    visited_conn.add(0)
    while queue_conn:
        curr = queue_conn.popleft()
        for nbr in adj_graph[curr]:
            if nbr not in visited_conn:
                visited_conn.add(nbr)
                queue_conn.append(nbr)

    if len(visited_conn) < v_count:
        raise TopologyError(
            f"Disconnected geometry detected: visited {len(visited_conn)} of {v_count} nodes."
        )

    # 5. Cycle rank check: for planar connected graph, m = E - V + 1
    cycle_rank = e_count - v_count + 1
    if cycle_rank < 1:
        raise TopologyError(
            f"No closed cells detected (cycle rank m = {cycle_rank} < 1). "
            "For open sections, use Section."
        )

    # 6. Planar half-edge embedding
    # Directed half-edge: (from_node, to_node, seg_idx, orientation)
    # where orientation = +1 if from_node == canonical_edges[seg_idx][0], else -1
    outgoing_half_edges: list[list[tuple[int, int, int, int]]] = [[] for _ in range(v_count)]
    for seg_idx, (u, v) in enumerate(canonical_edges):
        outgoing_half_edges[u].append((u, v, seg_idx, +1))
        outgoing_half_edges[v].append((v, u, seg_idx, -1))

    # At each canonical node, sort outgoing half-edges counter-clockwise by polar angle
    for u in range(v_count):
        u_node = canonical_nodes[u]

        def polar_angle(he: tuple[int, int, int, int]) -> float:
            target_node = canonical_nodes[he[1]]
            dx = target_node.x - u_node.x
            dy = target_node.y - u_node.y
            return math.atan2(dy, dx)

        outgoing_half_edges[u].sort(key=polar_angle)

    # Precompute cyclic successor:
    # When arriving at node v via half-edge (u, v, s, ori), its reverse outgoing half-edge
    # at v is (v, u, s, -ori).
    # In the CCW cyclic order at node v, the half-edge that turns tightest left into the face
    # is the one immediately preceding (v, u, s, -ori).
    next_half_edge: dict[tuple[int, int, int], tuple[int, int, int, int]] = {}
    for u in range(v_count):
        nbrs = outgoing_half_edges[u]
        k = len(nbrs)
        for idx, he in enumerate(nbrs):
            # he is outgoing from u: (u, v, s, ori)
            # The previous half-edge in cyclic CCW order at u:
            prev_he = nbrs[(idx - 1) % k]
            # Key: incoming half-edge to u, which is reverse of he: (v, u, s)
            incoming_key = (he[1], u, he[2])
            next_half_edge[incoming_key] = prev_he

    # 7. Face walk
    visited_half_edges: set[tuple[int, int, int]] = set()
    raw_faces: list[list[tuple[int, int, int, int]]] = []

    for u in range(v_count):
        for he in outgoing_half_edges[u]:
            he_key = (he[0], he[1], he[2])
            if he_key in visited_half_edges:
                continue

            # Trace face boundary
            face_walk: list[tuple[int, int, int, int]] = []
            curr_he = he
            while True:
                curr_key = (curr_he[0], curr_he[1], curr_he[2])
                if curr_key in visited_half_edges:
                    break
                visited_half_edges.add(curr_key)
                face_walk.append(curr_he)

                inc_key = (curr_he[0], curr_he[1], curr_he[2])
                curr_he = next_half_edge[inc_key]

            if face_walk:
                raw_faces.append(face_walk)

    # 8. Classify faces by shoelace signed area
    bounded_cells: list[Cell] = []
    exterior_faces: list[list[tuple[int, int, int, int]]] = []

    # Map each half-edge to the face type it belongs to:
    # face_id >= 0 for bounded cell id, -1 for exterior face
    he_to_cell_id: dict[tuple[int, int, int], int] = {}

    for face in raw_faces:
        # Polygon vertices
        poly_nodes = [canonical_nodes[he[0]] for he in face]
        num_pts = len(poly_nodes)
        # Translation-invariant shoelace formula anchored at the first vertex:
        # Prevents catastrophic cancellation when coordinates are large (e.g. 1e12)
        x0 = poly_nodes[0].x
        y0 = poly_nodes[0].y
        signed_area = 0.5 * sum(
            (poly_nodes[k].x - x0) * (poly_nodes[(k + 1) % num_pts].y - y0)
            - (poly_nodes[(k + 1) % num_pts].x - x0) * (poly_nodes[k].y - y0)
            for k in range(num_pts)
        )

        if signed_area == 0.0:
            raise GeometryError(
                "Degenerate zero-area face detected during planar cell extraction."
            )

        if signed_area > 0.0:
            # Counter-clockwise bounded cell!
            cell_id = len(bounded_cells)
            node_idx_tuple = tuple(he[0] for he in face)
            seg_idx_tuple = tuple(he[2] for he in face)
            ori_tuple = tuple(he[3] for he in face)
            perim = sum(segments[he[2]].length for he in face)

            # Scale-aware cell area validation: A_c must be strictly positive and not collapsed.
            # Relative check: |A| / perim^2 <= 1e-12 using frexp to avoid overflow/underflow
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
                    f"Cell {cell_id} has collapsed/zero area (A={signed_area:.6e}, perimeter={perim:.6e})."
                )

            cell_obj = Cell(
                id=cell_id,
                node_indices=node_idx_tuple,
                segment_indices=seg_idx_tuple,
                segment_orientations=ori_tuple,
                area=float(signed_area),
                perimeter=float(perim),
            )
            bounded_cells.append(cell_obj)

            for he in face:
                he_to_cell_id[(he[0], he[1], he[2])] = cell_id
        else:
            # Exterior unbounded face
            exterior_faces.append(face)
            for he in face:
                he_to_cell_id[(he[0], he[1], he[2])] = -1

    n_c = len(bounded_cells)

    # 9. Validate Euler cycle rank
    if n_c != cycle_rank:
        raise TopologyError(
            f"Extracted bounded cell count ({n_c}) does not match cycle rank ({cycle_rank} = E - V + 1). "
            "Section geometry is non-planar or invalid."
        )

    # 10. Check that every physical segment belongs to 1 or 2 bounded cells,
    # and EXPLICITLY REJECT mixed open/closed topology (dangling open branches)
    # A segment whose two half-edges both belong to the exterior face (-1) is a dangling branch!
    for s_idx, (u, v) in enumerate(canonical_edges):
        cell_fwd = he_to_cell_id.get((u, v, s_idx), -1)
        cell_rev = he_to_cell_id.get((v, u, s_idx), -1)

        # Count how many bounded cells this segment belongs to
        bounded_cell_count = (1 if cell_fwd >= 0 else 0) + (1 if cell_rev >= 0 else 0)

        if bounded_cell_count == 0:
            raise GeometryError(
                f"Mixed open/closed topology detected: segment {s_idx} connects nodes ({u}, {v}) "
                "outside of any bounded cell. A section containing both closed cells and "
                "dangling open branches is explicitly not supported in v0.5."
            )

        if cell_fwd >= 0 and cell_rev >= 0 and cell_fwd == cell_rev:
            raise GeometryError(
                f"Segment {s_idx} is a bridge/antenna inside cell {cell_fwd}. "
                "Such geometry is unsupported in v0.5."
            )

    # 11. Assemble oriented cell-wall incidence matrix B: shape (n_c, E)
    # B_{ci} = +1 if segment i follows CCW boundary of cell c,
    #          -1 if segment i opposes CCW boundary of cell c,
    #           0 if segment i is not in cell c.
    B = np.zeros((n_c, e_count), dtype=float)
    for cell in bounded_cells:
        for s_idx, ori in zip(cell.segment_indices, cell.segment_orientations):
            B[cell.id, s_idx] = float(ori)

    # Verify shared wall condition: for any internal common wall, adjacent cells have opposite signs
    for s_idx in range(e_count):
        incident_cells = np.where(B[:, s_idx] != 0.0)[0]
        if len(incident_cells) == 2:
            c1, c2 = incident_cells
            if B[c1, s_idx] * B[c2, s_idx] != -1.0:
                raise TopologyError(
                    f"Shared wall {s_idx} does not have opposite incidence signs in cells {c1} and {c2}."
                )
        elif len(incident_cells) > 2:
            raise TopologyError(
                f"Segment {s_idx} belongs to {len(incident_cells)} cells (maximum 2 allowed)."
            )

    # 12. Assemble dimensionless cell flexibility matrix H_scaled = B @ diag(rho_scaled) @ B^T
    # Decompose L_i / t_i using frexp to avoid intermediate float overflow when L_i / t_i > 1.79e308
    m_rhos = []
    e_rhos = []
    for seg in segments:
        mL, eL = math.frexp(seg.length)
        mt, et = math.frexp(seg.t)
        m_r, e_norm = math.frexp(mL / mt)
        m_rhos.append(m_r)
        e_rhos.append((eL - et) + e_norm)

    e_max_H = max(e_rhos)
    rho_scaled = np.zeros(e_count, dtype=float)
    for i in range(e_count):
        delta_e = e_rhos[i] - e_max_H
        if delta_e >= -1100:
            rho_scaled[i] = math.ldexp(m_rhos[i], delta_e)
        else:
            rho_scaled[i] = 0.0

    H_scaled = B @ np.diag(rho_scaled) @ B.T

    # 13. Scale-aware positive-definiteness validation on strictly finite dimensionless H_scaled
    eigvals = np.linalg.eigvalsh(H_scaled)
    lambda_min = float(eigvals[0])
    lambda_max = float(eigvals[-1])
    eps_mach = float(np.finfo(float).eps)
    h_threshold = safety_factor * eps_mach * lambda_max

    if lambda_min <= h_threshold or lambda_min <= 0.0:
        raise SingularSectionError(
            f"Cell flexibility matrix H is singular or rank-deficient: lambda_min={lambda_min:.6e}, "
            f"threshold={h_threshold:.6e}."
        )

    # Reconstruct physical H if strictly within IEEE-754 float64 range
    _H: np.ndarray | None = None
    try:
        H_physical = np.empty_like(H_scaled, dtype=np.float64)
        for r in range(H_scaled.shape[0]):
            for c in range(H_scaled.shape[1]):
                val = float(H_scaled[r, c])
                if val == 0.0:
                    H_physical[r, c] = 0.0
                else:
                    m, e = math.frexp(val)
                    H_physical[r, c] = math.ldexp(m, e + e_max_H)
        _H = H_physical
    except OverflowError:
        _H = None

    return CellTopology(
        canonical_nodes=tuple(canonical_nodes),
        canonical_edges=tuple(canonical_edges),
        cells=tuple(bounded_cells),
        B=B,
        H_scaled=H_scaled,
        e_max_H=e_max_H,
        _H=_H,
    )
