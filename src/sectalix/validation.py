"""Input and topology validation for Sectalix v0.1."""

from __future__ import annotations

from collections import deque
import math
from typing import Sequence

from sectalix.exceptions import GeometryError, TopologyError
from sectalix.primitives import Node, Segment


def _validate_node_tolerance(tol: float) -> None:
    """Require a finite, non-negative absolute coordinate tolerance."""
    if not math.isfinite(tol) or tol < 0.0:
        raise GeometryError(
            f"The node tolerance must be finite and non-negative. Got {tol}."
        )


def cluster_nodes(
    nodes: Sequence[Node], tol: float = 1e-9
) -> tuple[list[Node], dict[int, int]]:
    """Cluster nodes by transitive spatial proximity within tolerance.

    Clusters are the connected components of the undirected proximity graph in
    which two nodes are adjacent when their Euclidean separation is at most
    ``tol``.  This single-linkage definition makes cluster membership invariant
    to node/segment input order.  ``tol`` is an absolute distance expressed in
    the same units as the coordinates.

    Args:
        nodes: Sequence of Node objects to cluster.
        tol: Euclidean distance tolerance for merging nodes.

    Returns:
        A tuple of:
            - canonical_nodes: List of unique representative Node objects.
            - mapping: Dict mapping original node index (in `nodes`) to canonical node index.
    """
    _validate_node_tolerance(tol)

    parent = list(range(len(nodes)))
    rank = [0] * len(nodes)

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(first: int, second: int) -> None:
        first_root = find(first)
        second_root = find(second)
        if first_root == second_root:
            return
        if rank[first_root] < rank[second_root]:
            first_root, second_root = second_root, first_root
        parent[second_root] = first_root
        if rank[first_root] == rank[second_root]:
            rank[first_root] += 1

    # This is intentionally all-pairs: v0.1 sections are small, and the simple
    # proximity graph is easier to audit than a spatial-index implementation.
    for first in range(len(nodes)):
        for second in range(first + 1, len(nodes)):
            if nodes[first].is_close(nodes[second], tol=tol):
                union(first, second)

    canonical_nodes: list[Node] = []
    mapping: dict[int, int] = {}
    root_to_canonical: dict[int, int] = {}
    for index, node in enumerate(nodes):
        root = find(index)
        if root not in root_to_canonical:
            root_to_canonical[root] = len(canonical_nodes)
            canonical_nodes.append(node)
        mapping[index] = root_to_canonical[root]

    return canonical_nodes, mapping


def _cross(ax: float, ay: float, bx: float, by: float) -> float:
    """Return the scalar 2D cross product of vectors a and b."""
    return ax * by - ay * bx


def _point_lies_on_segment(point: Node, segment: Segment, tol: float) -> bool:
    """Test point/segment incidence using perpendicular distance and projection."""
    dx = segment.p2.x - segment.p1.x
    dy = segment.p2.y - segment.p1.y
    px = point.x - segment.p1.x
    py = point.y - segment.p1.y
    length = segment.length
    perpendicular_distance = abs(_cross(dx, dy, px, py)) / length
    if perpendicular_distance > tol:
        return False
    projection = (px * dx + py * dy) / length
    return -tol <= projection <= length + tol


def _segments_intersect(first: Segment, second: Segment, tol: float) -> bool:
    """Return whether two finite segments intersect within an absolute tolerance.

    The proper-intersection branch is the standard orientation test: endpoints
    of each segment must lie on opposite sides of the other supporting line.
    Point-on-segment checks cover endpoint contact and collinear cases.
    """
    r_x = first.p2.x - first.p1.x
    r_y = first.p2.y - first.p1.y
    s_x = second.p2.x - second.p1.x
    s_y = second.p2.y - second.p1.y

    q1_x = second.p1.x - first.p1.x
    q1_y = second.p1.y - first.p1.y
    q2_x = second.p2.x - first.p1.x
    q2_y = second.p2.y - first.p1.y
    p1_x = first.p1.x - second.p1.x
    p1_y = first.p1.y - second.p1.y
    p2_x = first.p2.x - second.p1.x
    p2_y = first.p2.y - second.p1.y

    first_side_1 = _cross(r_x, r_y, q1_x, q1_y)
    first_side_2 = _cross(r_x, r_y, q2_x, q2_y)
    second_side_1 = _cross(s_x, s_y, p1_x, p1_y)
    second_side_2 = _cross(s_x, s_y, p2_x, p2_y)
    first_cross_tol = tol * first.length
    second_cross_tol = tol * second.length

    crosses_first = (first_side_1 > first_cross_tol and first_side_2 < -first_cross_tol) or (
        first_side_1 < -first_cross_tol and first_side_2 > first_cross_tol
    )
    crosses_second = (
        second_side_1 > second_cross_tol and second_side_2 < -second_cross_tol
    ) or (second_side_1 < -second_cross_tol and second_side_2 > second_cross_tol)
    if crosses_first and crosses_second:
        return True

    return any(
        (
            _point_lies_on_segment(first.p1, second, tol),
            _point_lies_on_segment(first.p2, second, tol),
            _point_lies_on_segment(second.p1, first, tol),
            _point_lies_on_segment(second.p2, first, tol),
        )
    )


def _collinear_overlap_length(first: Segment, second: Segment, tol: float) -> float:
    """Return projected overlap when segments are collinear within tolerance."""
    if not (
        _point_lies_on_segment(second.p1, first, tol)
        or _point_lies_on_segment(second.p2, first, tol)
        or _point_lies_on_segment(first.p1, second, tol)
        or _point_lies_on_segment(first.p2, second, tol)
    ):
        return 0.0

    dx = first.p2.x - first.p1.x
    dy = first.p2.y - first.p1.y
    length = first.length
    # Both endpoints of the second segment must be close to the first segment's
    # supporting line; otherwise this is a transverse intersection, not overlap.
    line_distance_1 = abs(
        _cross(dx, dy, second.p1.x - first.p1.x, second.p1.y - first.p1.y)
    ) / length
    line_distance_2 = abs(
        _cross(dx, dy, second.p2.x - first.p1.x, second.p2.y - first.p1.y)
    ) / length
    if line_distance_1 > tol or line_distance_2 > tol:
        return 0.0

    unit_x = dx / length
    unit_y = dy / length
    second_start = (
        (second.p1.x - first.p1.x) * unit_x
        + (second.p1.y - first.p1.y) * unit_y
    )
    second_end = (
        (second.p2.x - first.p1.x) * unit_x
        + (second.p2.y - first.p1.y) * unit_y
    )
    overlap_start = max(0.0, min(second_start, second_end))
    overlap_end = min(length, max(second_start, second_end))
    return max(0.0, overlap_end - overlap_start)


def _validate_segment_intersections(
    segments: Sequence[Segment], canonical_edges: Sequence[tuple[int, int]], tol: float
) -> None:
    """Reject overlapping walls and intersections not represented by endpoints."""
    for first_index in range(len(segments)):
        for second_index in range(first_index + 1, len(segments)):
            first = segments[first_index]
            second = segments[second_index]
            shared_nodes = set(canonical_edges[first_index]) & set(
                canonical_edges[second_index]
            )

            if _collinear_overlap_length(first, second, tol) > tol:
                raise TopologyError(
                    f"Segments {first_index} and {second_index} overlap along their "
                    "centerlines; overlapping wall length would be counted twice."
                )

            if not shared_nodes and _segments_intersect(first, second, tol):
                raise TopologyError(
                    f"Segments {first_index} and {second_index} intersect without a shared "
                    "endpoint node. Split all walls at intersections before v0.1 analysis."
                )


def validate_section_geometry_and_topology(
    segments: Sequence[Segment], node_tolerance: float = 1e-9
) -> tuple[list[Node], list[tuple[int, int]]]:
    """Validate geometry and topology for v0.1 open thin-walled sections.

    Enforces all acceptance criteria:
    - Non-empty section
    - All coordinates and thicknesses finite
    - All thicknesses strictly positive (t > 0)
    - All segments strictly positive length (L > 0)
    - No duplicate segments
    - Connected graph (single component)
    - Open section only (no closed loops / cycles; trees and branched trees allowed)

    Args:
        segments: Sequence of Segment objects.
        node_tolerance: Absolute coordinate tolerance, in the geometry's units,
            for identifying shared nodes.

    Returns:
        A tuple of:
            - canonical_nodes: List of unique canonical nodes.
            - canonical_edges: List of (u, v) canonical node index pairs for each segment.

    Raises:
        GeometryError: If any geometric property (coordinates, thickness, length) is invalid.
        TopologyError: If topology is empty, disconnected, cyclical, or has duplicate segments.
    """
    _validate_node_tolerance(node_tolerance)

    if len(segments) == 0:
        raise TopologyError("Section must contain at least one segment.")

    # 1. Individual segment validation (already checked in dataclass __post_init__,
    #    re-verified here for completeness)
    all_nodes: list[Node] = []
    for i, seg in enumerate(segments):
        if not math.isfinite(seg.t) or seg.t <= 0.0:
            raise GeometryError(
                f"Segment {i} has invalid thickness t={seg.t}. Thickness must be finite and > 0."
            )
        if not math.isfinite(seg.length) or seg.length == 0.0:
            raise GeometryError(
                f"Segment {i} has invalid length (L={seg.length:.6e}); length must be "
                "finite and strictly positive."
            )
        all_nodes.append(seg.p1)
        all_nodes.append(seg.p2)

    # 2. Cluster endpoints into canonical nodes
    canonical_nodes, mapping = cluster_nodes(all_nodes, tol=node_tolerance)

    # 3. Build edge list in terms of canonical node indices
    canonical_edges: list[tuple[int, int]] = []
    seen_edges: dict[tuple[int, int], int] = {}

    for seg_idx in range(len(segments)):
        u = mapping[2 * seg_idx]
        v = mapping[2 * seg_idx + 1]

        if u == v:
            raise GeometryError(
                f"Segment {seg_idx} has coincident endpoints within tolerance {node_tolerance:.1e}: "
                f"p1=({segments[seg_idx].p1.x}, {segments[seg_idx].p1.y}), "
                f"p2=({segments[seg_idx].p2.x}, {segments[seg_idx].p2.y})."
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

    # Graph topology alone misses crossings through segment interiors and
    # collinear overlaps.  These must be explicit endpoint nodes in v0.1.
    _validate_segment_intersections(segments, canonical_edges, node_tolerance)

    v_count = len(canonical_nodes)
    e_count = len(canonical_edges)

    # 4. Build adjacency list
    adj: list[list[int]] = [[] for _ in range(v_count)]
    for u, v in canonical_edges:
        adj[u].append(v)
        adj[v].append(u)

    # 5. Check connectivity via BFS/DFS
    visited: set[int] = set()
    queue: deque[int] = deque([0])
    visited.add(0)

    while queue:
        curr = queue.popleft()
        for nbr in adj[curr]:
            if nbr not in visited:
                visited.add(nbr)
                queue.append(nbr)

    if len(visited) < v_count:
        num_components = 1
        unvisited = set(range(v_count)) - visited
        while unvisited:
            num_components += 1
            start = next(iter(unvisited))
            q = deque([start])
            unvisited.remove(start)
            while q:
                c = q.popleft()
                for nb in adj[c]:
                    if nb in unvisited:
                        unvisited.remove(nb)
                        q.append(nb)
        raise TopologyError(
            f"Disconnected geometry detected: section has {num_components} disconnected components "
            f"(visited {len(visited)} of {v_count} nodes)."
        )

    # 6. Check for closed loops (cycles)
    # For a connected undirected graph, it is a tree (open section) if and only if E = V - 1.
    # If E >= V, there is at least one cycle.
    if e_count >= v_count:
        # Perform DFS cycle tracing to provide informative diagnostics
        cycle_parent: dict[int, int | None] = {0: None}
        cycle_visited: set[int] = set([0])
        cycle_found: list[int] = []

        def dfs(node: int, parent: int | None) -> bool:
            for neighbor in adj[node]:
                if neighbor == parent:
                    continue
                if neighbor in cycle_visited:
                    # Cycle found! Backtrack path
                    cycle_found.append(neighbor)
                    curr: int | None = node
                    while curr is not None and curr != neighbor:
                        cycle_found.append(curr)
                        curr = cycle_parent.get(curr)
                    cycle_found.append(neighbor)
                    return True
                cycle_visited.add(neighbor)
                cycle_parent[neighbor] = node
                if dfs(neighbor, node):
                    return True
            return False

        dfs(0, None)
        cycle_desc = " -> ".join(str(n) for n in cycle_found) if cycle_found else "cycle detected"
        raise TopologyError(
            f"Closed-loop geometry detected: unsupported in v0.1 open-section analysis. "
            f"Found loop involving canonical nodes: {cycle_desc} (E={e_count} >= V={v_count})."
        )

    return canonical_nodes, canonical_edges
