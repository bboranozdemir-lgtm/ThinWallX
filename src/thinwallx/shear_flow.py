"""Exact thin-walled shear-flow analysis for open sections (ThinWallX v0.2).

Solves C * alpha = V for arbitrary transverse shear load V = [Vx, Vy]^T and computes
exact closed-form shear flow q(s) along straight centerline segments for open tree topologies.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
import math
from typing import TYPE_CHECKING, Sequence

import numpy as np

from thinwallx.exceptions import GeometryError, SingularSectionError
from thinwallx.primitives import Segment
from thinwallx.shear_load import ShearLoad
from thinwallx.validation import cluster_nodes

if TYPE_CHECKING:
    from thinwallx.section import Section


@dataclass(frozen=True)
class SegmentShearFlow:
    """Exact shear-flow distribution along a single straight centerline segment.

    Parameterization:
        s in [0, L] measured from node1 to node2.
        xi = s / L in [0, 1].

    Scalar shear flow:
        q(xi) = a0 + a1 * xi + a2 * xi^2
        q(s)  = a0 + (a1 / L) * s + (a2 / L^2) * s^2

    Physical shear-flow vector:
        q_vec(s) = q(s) * tangent

    Attributes:
        segment: The underlying straight Segment.
        a0: Constant term of quadratic polynomial in xi = s/L.
        a1: Linear coefficient in xi = s/L.
        a2: Quadratic coefficient in xi = s/L.
        subtree_moment: Centroid-relative first-moment vector [Qy, Qx] of the
            node2-side subtree attached to node2.
    """

    segment: Segment
    a0: float
    a1: float
    a2: float
    subtree_moment: tuple[float, float]

    def __post_init__(self) -> None:
        """Validate polynomial coefficients and moments for finiteness."""
        if not (math.isfinite(self.a0) and math.isfinite(self.a1) and math.isfinite(self.a2)):
            raise GeometryError(
                f"SegmentShearFlow polynomial coefficients must be finite. "
                f"Got a0={self.a0}, a1={self.a1}, a2={self.a2}."
            )
        if not (math.isfinite(self.subtree_moment[0]) and math.isfinite(self.subtree_moment[1])):
            raise GeometryError(
                f"SegmentShearFlow subtree moment must be finite. Got {self.subtree_moment}."
            )

    @property
    def tangent(self) -> np.ndarray:
        """Unit tangent vector pointing from node1 to node2: [tx, ty]^T."""
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
        return self.a0

    @property
    def q_mid(self) -> float:
        """Scalar shear flow at midpoint (s=L/2, xi=0.5)."""
        return self.a0 + 0.5 * self.a1 + 0.25 * self.a2

    @property
    def q_end(self) -> float:
        """Scalar shear flow at end (s=L, xi=1)."""
        return self.a0 + self.a1 + self.a2

    def q_at_xi(self, xi: float) -> float:
        """Evaluate scalar shear flow q at normalized position xi in [0, 1]."""
        if not math.isfinite(xi):
            raise GeometryError(f"Normalized coordinate xi must be finite. Got xi={xi}.")
        if xi < -1e-7 or xi > 1.0 + 1e-7:
            raise GeometryError(
                f"Normalized coordinate xi must be in [0, 1]. Got xi={xi}."
            )
        xi_clamped = max(0.0, min(1.0, xi))
        return self.a0 + self.a1 * xi_clamped + self.a2 * (xi_clamped * xi_clamped)

    def q_at(self, s: float) -> float:
        """Evaluate scalar shear flow q at position s in [0, L]."""
        if not math.isfinite(s):
            raise GeometryError(f"Position s must be finite. Got s={s}.")
        return self.q_at_xi(s / self.length)

    def vector_at_xi(self, xi: float) -> np.ndarray:
        """Evaluate physical shear-flow vector q(xi) * tangent at xi in [0, 1]."""
        return self.q_at_xi(xi) * self.tangent

    def vector_at(self, s: float) -> np.ndarray:
        """Evaluate physical shear-flow vector q(s) * tangent at s in [0, L]."""
        return self.q_at(s) * self.tangent

    def scalar_integral(self) -> float:
        """Exact integral of scalar shear flow along the segment: \\int_0^L q(s) ds."""
        # \\int_0^1 (a0 + a1*xi + a2*xi^2) dxi = a0 + a1/2 + a2/3
        return self.length * (self.a0 + 0.5 * self.a1 + (1.0 / 3.0) * self.a2)

    def resultant_integral(self) -> np.ndarray:
        """Exact vector integral: \\int_0^L q(s)*tangent ds = (\\int_0^L q(s) ds) * tangent."""
        return self.scalar_integral() * self.tangent


@dataclass(frozen=True)
class ShearFlowResult:
    """Container holding the results of a thin-walled shear-flow analysis.

    Attributes:
        section: The analyzed Section.
        load: Applied transverse shear load.
        segment_flows: Tuple of SegmentShearFlow objects corresponding to section.segments.
        alpha: Solution vector alpha = [alpha_x, alpha_y]^T from C * alpha = V.
        c_matrix: Centroidal coordinate second-moment matrix C = [[Iy, Ixy], [Ixy, Ix]].
    """

    section: Section
    load: ShearLoad
    segment_flows: tuple[SegmentShearFlow, ...]
    alpha: np.ndarray
    c_matrix: np.ndarray

    def __getitem__(self, idx: int) -> SegmentShearFlow:
        """Get SegmentShearFlow by segment index."""
        return self.segment_flows[idx]

    def __len__(self) -> int:
        """Number of segments."""
        return len(self.segment_flows)

    def flow_at(self, segment_idx: int, xi: float) -> float:
        """Evaluate scalar shear flow on segment `segment_idx` at `xi` in [0, 1]."""
        return self.segment_flows[segment_idx].q_at_xi(xi)

    def vector_at(self, segment_idx: int, xi: float) -> np.ndarray:
        """Evaluate physical shear-flow vector on segment `segment_idx` at `xi` in [0, 1]."""
        return self.segment_flows[segment_idx].vector_at_xi(xi)

    @property
    def recovered_resultant(self) -> np.ndarray:
        """Sum of integrated physical shear flow vectors: \\sum_i \\int_0^L q_i(s) t_i ds."""
        total = np.zeros(2, dtype=float)
        for sf in self.segment_flows:
            total += sf.resultant_integral()
        return total

    @property
    def resultant_error(self) -> float:
        """Euclidean error norm between recovered resultant and applied load vector."""
        return float(np.linalg.norm(self.recovered_resultant - self.load.vector))


def _compute_open_tree_flows(
    segments: Sequence[Segment],
    canonical_nodes: Sequence[Node],
    canonical_edges: Sequence[tuple[int, int]],
    total_area: float,
    alpha: np.ndarray,
    ref_origin: tuple[float, float] | None = None,
) -> list[SegmentShearFlow]:
    """Compute exact shear flow fields for an open tree network.

    Args:
        segments: Sequence of Segment objects in the tree.
        canonical_nodes: Sequence of canonical Node objects.
        canonical_edges: Canonical edge index pairs (u, v) for each segment.
        total_area: Total physical section area A.
        alpha: Bending gradient vector from C * alpha = V.
        ref_origin: Optional (ref_x, ref_y) coordinate shift origin.

    Returns:
        List of SegmentShearFlow objects for each segment.
    """
    v_count = len(canonical_nodes)
    if ref_origin is None:
        ref_x = min(n.x for n in canonical_nodes)
        ref_y = min(n.y for n in canonical_nodes)
    else:
        ref_x, ref_y = ref_origin

    # Centroid relative to local reference:
    cx_loc = sum(
        seg.area * ((seg.p1.x - ref_x) + (seg.p2.x - ref_x)) / 2.0
        for seg in segments
    ) / total_area
    cy_loc = sum(
        seg.area * ((seg.p1.y - ref_y) + (seg.p2.y - ref_y)) / 2.0
        for seg in segments
    ) / total_area

    # Adjacency list: adj[u] -> list of (v, seg_index)
    adj: list[list[tuple[int, int]]] = [[] for _ in range(v_count)]
    for s_idx, (u, v) in enumerate(canonical_edges):
        adj[u].append((v, s_idx))
        adj[v].append((u, s_idx))

    # Precompute full centroid-relative first-moment vector for each segment
    seg_full_moments: list[np.ndarray] = []
    for seg in segments:
        xc1 = (seg.p1.x - ref_x) - cx_loc
        yc1 = (seg.p1.y - ref_y) - cy_loc
        xc2 = (seg.p2.x - ref_x) - cx_loc
        yc2 = (seg.p2.y - ref_y) - cy_loc
        length = seg.length
        t = seg.t
        qy = t * length * (xc1 + xc2) / 2.0
        qx = t * length * (yc1 + yc2) / 2.0
        seg_full_moments.append(np.array([qy, qx], dtype=float))

    segment_flows: list[SegmentShearFlow] = []
    for i, seg in enumerate(segments):
        u_node, v_node = canonical_edges[i]
        # v_node is the node2 endpoint of segment i
        # Find all segments in the node2-side subtree by BFS from v_node without crossing segment i
        visited_nodes: set[int] = {v_node}
        queue: deque[int] = deque([v_node])
        subtree_seg_indices: list[int] = []

        while queue:
            curr = queue.popleft()
            for neighbor, s_idx in adj[curr]:
                if s_idx == i:
                    continue
                if neighbor not in visited_nodes:
                    visited_nodes.add(neighbor)
                    queue.append(neighbor)
                    subtree_seg_indices.append(s_idx)

        m_subtree = np.zeros(2, dtype=float)
        for s_idx in subtree_seg_indices:
            m_subtree += seg_full_moments[s_idx]

        if not np.all(np.isfinite(m_subtree)):
            raise GeometryError(
                f"Non-finite subtree moment encountered for segment {i}: m_subtree={m_subtree}."
            )

        xc1 = (seg.p1.x - ref_x) - cx_loc
        yc1 = (seg.p1.y - ref_y) - cy_loc
        xc2 = (seg.p2.x - ref_x) - cx_loc
        yc2 = (seg.p2.y - ref_y) - cy_loc
        d_xc = xc2 - xc1
        d_yc = yc2 - yc1
        t = seg.t
        length = seg.length

        rc1_q = np.array([xc1, yc1], dtype=float)
        d_rc_q = np.array([d_xc, d_yc], dtype=float)

        with np.errstate(over="raise"):
            try:
                m_total_at_start = m_subtree + seg_full_moments[i]
                a0 = float(np.dot(alpha, m_total_at_start))
                a1 = -float(t * length * np.dot(alpha, rc1_q))
                a2 = -float(0.5 * t * length * np.dot(alpha, d_rc_q))
            except (FloatingPointError, OverflowError):
                raise OverflowError(
                    f"Shear flow polynomial coefficients overflow IEEE-754 float64 on segment {i} under applied load."
                )

        if not (math.isfinite(a0) and math.isfinite(a1) and math.isfinite(a2)):
            raise OverflowError(
                f"Non-finite polynomial coefficients encountered on segment {i}: "
                f"a0={a0}, a1={a1}, a2={a2}."
            )

        segment_flows.append(
            SegmentShearFlow(
                segment=seg,
                a0=a0,
                a1=a1,
                a2=a2,
                subtree_moment=(float(m_subtree[0]), float(m_subtree[1])),
            )
        )

    return segment_flows


def calculate_shear_flow(
    section: Section,
    vx: float | ShearLoad,
    vy: float | None = None,
    safety_factor: float = 1e4,
) -> ShearFlowResult:
    """Compute exact thin-walled shear-flow distribution for an open section.

    Solves the exact thin-walled bending-shear equilibrium equation:
        dq/ds = -alpha_x * (t * x_c(s)) - alpha_y * (t * y_c(s))

    Algorithm:
    1. Construct centroidal coordinate second-moment matrix C:
       C = [[Iy,  Ixy],
            [Ixy, Ix ]]
    2. Numerical Positive-Definiteness check:
       Eigenvalues of C: lambda_min > safety_factor * eps_mach * lambda_max
       If singular or rank-deficient, raises SingularSectionError.
    3. Solve:
       C * alpha = V for alpha = [alpha_x, alpha_y]^T
    4. Cut-side first moment for straight segment i with s in [0, L_i]:
       m_i(s) = m_{node2-side subtree} + m_{partial}(s)
       where m_{partial}(s) = t_i * \\int_s^{L_i} r_c(u) du
       m_i(s) = [Q_y(s), Q_x(s)]^T
    5. Scalar shear flow:
       q_i(s) = alpha^T * m_i(s) = alpha_x * Q_y(s) + alpha_y * Q_x(s)
       Polynomial in xi = s/L_i: q_i(xi) = a0 + a1*xi + a2*xi^2.

    Args:
        section: Validated open Section object.
        vx: Transverse shear force in x direction, or a ShearLoad object.
        vy: Transverse shear force in y direction (ignored if vx is ShearLoad).
        safety_factor: Safety multiplier on machine epsilon for positive definiteness.

    Returns:
        ShearFlowResult object containing segment flows and resultant verification.

    Raises:
        GeometryError: If shear loads are non-finite.
        SingularSectionError: If the section inertia matrix C is singular/rank-deficient.
    """
    from thinwallx.section import Section

    if not isinstance(section, Section):
        raise GeometryError(f"Expected a Section instance, got {type(section).__name__}.")
    if hasattr(section, "is_closed") and section.is_closed:
        from thinwallx.closed_shear_flow import calculate_closed_shear_flow
        return calculate_closed_shear_flow(section, vx=vx, vy=vy, safety_factor=safety_factor)  # type: ignore[return-value]

    section.validate()

    if not math.isfinite(safety_factor) or safety_factor <= 0.0:
        raise GeometryError(
            f"safety_factor must be finite and strictly positive. Got safety_factor={safety_factor}."
        )

    if isinstance(vx, ShearLoad):
        load = vx
    else:
        if vy is None:
            raise GeometryError("vy must be provided when vx is a scalar force.")
        load = ShearLoad(vx=float(vx), vy=float(vy))

    # 1. Construct centroidal coordinate second-moment matrix C
    ix = section.Ix
    iy = section.Iy
    ixy = section.Ixy

    c_matrix = np.array([[iy, ixy], [ixy, ix]], dtype=float)

    # 2. Check numerical positive definiteness of C
    # Eigenvalues of 2x2 symmetric matrix
    eigvals = np.linalg.eigvalsh(c_matrix)  # Sorted ascending: [lambda_min, lambda_max]
    lambda_min = float(eigvals[0])
    lambda_max = float(eigvals[1])

    eps_mach = float(np.finfo(float).eps)
    threshold = safety_factor * eps_mach * lambda_max

    if lambda_min <= threshold or lambda_min <= 0.0:
        raise SingularSectionError(
            f"Centroidal coordinate second-moment matrix C is singular or effectively rank-deficient. "
            f"lambda_min={lambda_min:.3e}, lambda_max={lambda_max:.3e}, threshold={threshold:.3e}. "
            f"Open-section 2D shear-flow analysis requires a positive-definite C matrix."
        )

    # 3. Solve C * alpha = V directly
    alpha = np.linalg.solve(c_matrix, load.vector)
    if not np.all(np.isfinite(alpha)):
        raise GeometryError(
            f"Non-finite solution vector alpha encountered in shear-flow solve: alpha={alpha}."
        )

    # 4. Analyze open tree topology to partition cut-side components
    segments = section.segments
    all_nodes = [node for seg in segments for node in (seg.p1, seg.p2)]
    canonical_nodes, mapping = cluster_nodes(all_nodes, tol=section._node_tolerance)

    seg_count = len(segments)

    # Local reference shift:
    # Use deterministic reference origin to preserve full floating-point precision
    # when section is located at large translation offsets (e.g. 1e12):
    ref_x, ref_y = min((n.x, n.y) for n in all_nodes)
    total_area = section.area

    # Build canonical edges:
    canonical_edges = [
        (mapping[2 * s_idx], mapping[2 * s_idx + 1]) for s_idx in range(seg_count)
    ]

    segment_flows = _compute_open_tree_flows(
        segments=segments,
        canonical_nodes=canonical_nodes,
        canonical_edges=canonical_edges,
        total_area=total_area,
        alpha=alpha,
        ref_origin=(ref_x, ref_y),
    )

    result = ShearFlowResult(
        section=section,
        load=load,
        segment_flows=tuple(segment_flows),
        alpha=alpha,
        c_matrix=c_matrix,
    )
    if not np.all(np.isfinite(result.recovered_resultant)):
        raise GeometryError(
            f"Non-finite recovered resultant encountered: {result.recovered_resultant}."
        )
    return result
