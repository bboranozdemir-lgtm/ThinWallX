"""Exact open-section Saint-Venant torsion and warping analysis (Sectalix v0.4).

Computes:
1. Saint-Venant open-section torsion constant J = sum_i (1/3) * L_i * t_i^3
2. Raw and normalized principal sectorial coordinates omega(s) with shear center as pole.
3. Warping constant C_w = int_A omega^2 dA = sum_i (t_i * L_i / 3) * (omega_{i,1}^2 + omega_{i,1} * omega_{i,2} + omega_{i,2}^2)

Governing Equations & Conventions (ACTIVE_PHASE.md):
- Sectorial differential with shear center S = (x_s, y_s) as pole:
      d(omega_raw) = [r_S x dr]_z
- Along straight segment i with constant tangent t_i:
      p_i = [r_S(0) x t_i]_z = x_{S,1} * t_y - y_{S,1} * t_x = const
      omega_raw,i(s) = omega_raw,i(0) + p_i * s
      Delta omega_i = p_i * L_i
- Tree propagation:
      Exact propagation from arbitrary root node: omega_raw(root) = 0
- Area-weighted mean sectorial coordinate:
      bar{omega} = (1 / A) * sum_i t_i * (L_i / 2) * (omega_raw,i,1 + omega_raw,i,2)
- Normalized principal sectorial coordinate:
      omega(s) = omega_raw(s) - bar{omega}  ==>  int_A omega dA = 0
- Warping constant:
      C_w = sum_i (t_i * L_i / 3) * (omega_{i,1}^2 + omega_{i,1} * omega_{i,2} + omega_{i,2}^2)
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
import math
import sys
from typing import TYPE_CHECKING, Sequence

import numpy as np

from sectalix.exceptions import GeometryError
from sectalix.primitives import Node, Segment
from sectalix.validation import cluster_nodes

if TYPE_CHECKING:
    from sectalix.section import Section


@dataclass(frozen=True)
class SegmentWarping:
    """Exact warping properties along a single straight centerline segment.

    Parameters
    ----------
    segment : Segment
        The underlying straight centerline segment.
    omega1 : float
        Normalized sectorial coordinate at start node (s = 0).
    omega2 : float
        Normalized sectorial coordinate at end node (s = L).
    p : float
        Sectorial rate / perpendicular moment arm p_i = [r_S(0) x t_i]_z.
    cw_segment : float
        Segment contribution to warping constant:
        int_0^L omega(s)^2 t_i ds = (t_i * L_i / 3) * (omega1^2 + omega1 * omega2 + omega2^2).
    j_segment : float
        Segment contribution to Saint-Venant torsion constant: (1/3) * L_i * t_i^3.
    """

    segment: Segment
    omega1: float
    omega2: float
    p: float
    cw_segment: float
    j_segment: float

    def __post_init__(self) -> None:
        """Validate properties for finiteness and valid values."""
        for name, val in [
            ("omega1", self.omega1),
            ("omega2", self.omega2),
            ("p", self.p),
            ("cw_segment", self.cw_segment),
            ("j_segment", self.j_segment),
        ]:
            if not math.isfinite(val):
                raise GeometryError(
                    f"SegmentWarping property '{name}' must be finite. Got {val}."
                )
        if self.j_segment < 0.0:
            raise GeometryError(
                f"SegmentWarping j_segment cannot be negative. Got {self.j_segment}."
            )
        if self.cw_segment < 0.0:
            scale = max(abs(self.omega1), abs(self.omega2)) * math.sqrt(self.segment.t) * math.sqrt(self.length)
            m, e = math.frexp(scale)
            m_tol, e_tol = math.frexp(1e-12 * m * m)
            try:
                tol = math.ldexp(m_tol, e_tol + 2 * e)
            except OverflowError:
                tol = 0.0
            if self.cw_segment < -tol:
                raise GeometryError(
                    f"SegmentWarping cw_segment cannot be negative. Got {self.cw_segment}."
                )

    @property
    def length(self) -> float:
        """Segment length L."""
        return self.segment.length

    def omega_at_xi(self, xi: float) -> float:
        """Evaluate normalized sectorial coordinate at normalized position xi in [0, 1]."""
        if not math.isfinite(xi):
            raise GeometryError(f"Normalized coordinate xi must be finite. Got xi={xi}.")
        if xi < -1e-7 or xi > 1.0 + 1e-7:
            raise GeometryError(
                f"Normalized coordinate xi must be in [0, 1]. Got xi={xi}."
            )
        xi_clamped = max(0.0, min(1.0, xi))
        return (1.0 - xi_clamped) * self.omega1 + xi_clamped * self.omega2

    def omega_at(self, s: float) -> float:
        """Evaluate normalized sectorial coordinate at position s in [0, L]."""
        if not math.isfinite(s):
            raise GeometryError(f"Position s must be finite. Got s={s}.")
        return self.omega_at_xi(s / self.length)


@dataclass(frozen=True)
class TorsionWarpingResult:
    """Container holding results of open-section torsion and warping analysis.

    Attributes:
        section: The analyzed Section.
        J: Saint-Venant open-section torsion constant sum_i (1/3) * L_i * t_i^3.
        Cw: Warping constant int_A omega^2 dA.
        omega_mean: Area-weighted mean raw sectorial coordinate bar{omega}.
        shear_center: Coordinates of shear center (x_s, y_s) used as sectorial pole.
        node_omega: Tuple of normalized principal sectorial coordinates at section.nodes.
        segment_warpings: Tuple of SegmentWarping objects corresponding to section.segments.
        raw_node_omega: Tuple of raw sectorial coordinates at section.nodes before normalization.
    """

    section: Section
    J: float
    Cw: float
    omega_mean: float
    shear_center: tuple[float, float]
    node_omega: tuple[float, ...]
    segment_warpings: tuple[SegmentWarping, ...]
    raw_node_omega: tuple[float, ...]

    def __post_init__(self) -> None:
        """Validate torsion and warping results."""
        for name, val in [
            ("J", self.J),
            ("Cw", self.Cw),
            ("omega_mean", self.omega_mean),
            ("shear_center.x", self.shear_center[0]),
            ("shear_center.y", self.shear_center[1]),
        ]:
            if not math.isfinite(val):
                raise GeometryError(f"TorsionWarpingResult property '{name}' must be finite. Got {val}.")

        if self.J <= 0.0:
            raise GeometryError(f"Saint-Venant torsion constant J must be strictly positive. Got {self.J}.")

        if self.Cw < 0.0:
            raise GeometryError(f"Warping constant C_w cannot be negative. Got {self.Cw}.")

        for val in self.node_omega:
            if not math.isfinite(val):
                raise GeometryError(f"Non-finite node sectorial coordinate encountered: {val}.")

    def __getitem__(self, idx: int) -> SegmentWarping:
        """Get SegmentWarping by segment index."""
        return self.segment_warpings[idx]

    def __len__(self) -> int:
        """Number of segments."""
        return len(self.segment_warpings)

    def omega_at(self, segment_idx: int, xi: float) -> float:
        """Evaluate normalized sectorial coordinate on segment `segment_idx` at `xi` in [0, 1]."""
        return self.segment_warpings[segment_idx].omega_at_xi(xi)

    @property
    def integral_omega_da(self) -> float:
        """Exact integral of normalized sectorial coordinate over the section: int_A omega dA.

        By definition of zero-mean principal normalization, this must be zero within floating-point tolerance.
        """
        contributions = [
            sw.segment.area * 0.5 * (sw.omega1 + sw.omega2)
            for sw in self.segment_warpings
        ]
        return math.fsum(contributions)


def compute_torsion_warping(
    section: Section,
    root_node_idx: int = 0,
    safety_factor: float = 1e4,
) -> TorsionWarpingResult:
    """Compute open-section Saint-Venant torsion constant J, sectorial coordinates, and warping constant C_w.

    Governing Mathematical Formulation (ACTIVE_PHASE.md):
    1. Saint-Venant open-section torsion constant:
           J = sum_i (1/3) * L_i * t_i^3
       using compensated summation (math.fsum).
    2. Obtain accepted v0.3 shear-center location S = (x_s, y_s) and centroid-relative offset e_s = [e_x, e_y]^T.
    3. Position relative to shear center:
           r_S,k = r_{c,k} - e_s
       using local reference coordinates to maintain precision under large translations (1e12).
    4. Segment sectorial rate:
           p_i = [r_S(0) x t_i]_z = x_{S,1} * t_{y,i} - y_{S,1} * t_{x,i} = const
    5. Propagate raw sectorial coordinates along connected tree topology from root_node_idx:
           omega_raw(root) = 0
           omega_raw(v) - omega_raw(u) = p_i * L_i
    6. Area-weighted mean:
           bar{omega} = (1 / A) * sum_i t_i * (L_i / 2) * (omega_raw,i,1 + omega_raw,i,2)
    7. Normalized principal sectorial coordinate:
           omega = omega_raw - bar{omega}
    8. Exact warping constant:
           C_w = sum_i (t_i * L_i / 3) * (omega_{i,1}^2 + omega_{i,1} * omega_{i,2} + omega_{i,2}^2)

    Args:
        section: Validated open Section object.
        root_node_idx: Index of canonical node to initialize omega_raw = 0 (default 0).
        safety_factor: Safety multiplier on machine epsilon for positive definiteness check in shear-center solver.

    Returns:
        TorsionWarpingResult containing J, C_w, normalized node and segment fields.

    Raises:
        GeometryError: If section is invalid or intermediate/final values are non-finite.
        SingularSectionError: If the section inertia matrix is singular/rank-deficient.
    """
    from sectalix.section import Section

    if not isinstance(section, Section):
        raise GeometryError(f"Expected a Section instance, got {type(section).__name__}.")
    if hasattr(section, "is_closed") and section.is_closed:
        from sectalix.closed_torsion import compute_closed_torsion_warping
        return compute_closed_torsion_warping(
            section, root_node_idx=root_node_idx, safety_factor=safety_factor
        )

    section.validate()

    segments = section.segments
    if len(segments) == 0:
        raise GeometryError("Section must contain at least one segment.")

    # 1. Saint-Venant torsion constant J and segment contributions J_i
    # Decompose into mantissa and integer exponent using math.frexp and math.ldexp
    # to eliminate intermediate underflow/overflow risk under heterogeneous scales.
    m_list: list[float] = []
    e_list: list[int] = []
    for seg in segments:
        mL, eL = math.frexp(seg.length)
        mt, et = math.frexp(seg.t)
        # J_i = (1/3) * L * t^3 = ((mL * mt^3) / 3) * 2^(eL + 3*et)
        m = (mL * (mt ** 3)) / 3.0
        e = eL + 3 * et
        m_norm, e_norm = math.frexp(m)
        m_list.append(m_norm)
        e_list.append(e + e_norm)

    e_max = max(e_list)
    scaled_mantissas: list[float] = []
    for m, e in zip(m_list, e_list):
        shift = e - e_max
        if shift < -1100:
            scaled_mantissas.append(0.0)
        else:
            scaled_mantissas.append(math.ldexp(m, shift))

    total_mantissa = math.fsum(scaled_mantissas)
    try:
        J = math.ldexp(total_mantissa, e_max)
    except OverflowError:
        J = float("inf")

    if not math.isfinite(J) or J <= 0.0:
        raise GeometryError(
            f"Computed Saint-Venant torsion constant J must be finite and positive. Got {J}."
        )

    j_contributions: list[float] = []
    for m, e in zip(m_list, e_list):
        try:
            val = math.ldexp(m, e)
        except OverflowError:
            val = float("inf")
        j_contributions.append(val)

    # 2. Accepted v0.3 shear center pole
    sc_result = section.compute_shear_center(safety_factor=safety_factor)
    ex, ey = sc_result.ex, sc_result.ey
    xs, ys = sc_result.x, sc_result.y

    # Canonical nodes and mapping
    all_nodes = [node for seg in segments for node in (seg.p1, seg.p2)]
    canonical_nodes, mapping = cluster_nodes(all_nodes, tol=section._node_tolerance)
    v_count = len(canonical_nodes)

    if not (0 <= root_node_idx < v_count):
        raise GeometryError(
            f"root_node_idx must be in [0, {v_count - 1}]. Got {root_node_idx}."
        )

    # Local reference shift for coordinates
    ref_x, ref_y = min((n.x, n.y) for n in all_nodes)
    total_area = section.area

    cx_loc = sum(
        seg.area * ((seg.p1.x - ref_x) + (seg.p2.x - ref_x)) / 2.0
        for seg in segments
    ) / total_area
    cy_loc = sum(
        seg.area * ((seg.p1.y - ref_y) + (seg.p2.y - ref_y)) / 2.0
        for seg in segments
    ) / total_area

    # Graph adjacency for tree propagation
    # adj[u] -> list of (v, seg_idx, is_forward)
    adj: list[list[tuple[int, int, bool]]] = [[] for _ in range(v_count)]
    seg_rates: list[float] = []
    canonical_edges: list[tuple[int, int]] = []

    for s_idx, seg in enumerate(segments):
        u = mapping[2 * s_idx]
        v = mapping[2 * s_idx + 1]
        canonical_edges.append((u, v))

        # Position of segment.p1 relative to shear center:
        # r_{S,1} = r_{c,1} - e_s, computed with local reference shift directly
        # from segment's own coordinates (not cluster_nodes representative):
        xc1 = (seg.p1.x - ref_x) - cx_loc
        yc1 = (seg.p1.y - ref_y) - cy_loc
        rx1 = xc1 - ex
        ry1 = yc1 - ey

        # Unit tangent from p1 to p2:
        dx = seg.p2.x - seg.p1.x
        dy = seg.p2.y - seg.p1.y
        length = seg.length
        tx = dx / length
        ty = dy / length

        # Moment arm / sectorial rate p_i = [r_S(0) x t]_z
        p_i = rx1 * ty - ry1 * tx

        if not math.isfinite(p_i):
            raise GeometryError(
                f"Non-finite sectorial rate p_i computed on segment {seg.id}: p_i={p_i}."
            )

        seg_rates.append(p_i)
        adj[u].append((v, s_idx, True))
        adj[v].append((u, s_idx, False))

    # 3. Propagate raw sectorial coordinate from root_node_idx
    raw_node_omega = [0.0] * v_count
    visited = [False] * v_count
    queue = deque([root_node_idx])
    visited[root_node_idx] = True

    while queue:
        curr = queue.popleft()
        curr_val = raw_node_omega[curr]
        for neighbor, s_idx, is_forward in adj[curr]:
            if not visited[neighbor]:
                visited[neighbor] = True
                delta = seg_rates[s_idx] * segments[s_idx].length
                if is_forward:
                    raw_node_omega[neighbor] = curr_val + delta
                else:
                    raw_node_omega[neighbor] = curr_val - delta
                queue.append(neighbor)

    # 4. Area-weighted mean sectorial coordinate bar{omega}
    mean_numerator_terms: list[float] = []
    for s_idx, seg in enumerate(segments):
        u, v = canonical_edges[s_idx]
        w_u = raw_node_omega[u]
        w_v = raw_node_omega[v]
        mean_numerator_terms.append(seg.area * 0.5 * (w_u + w_v))

    omega_mean = math.fsum(mean_numerator_terms) / total_area
    if not math.isfinite(omega_mean):
        raise GeometryError(f"Non-finite mean sectorial coordinate bar{{omega}} computed: {omega_mean}.")

    # 5. Normalized principal sectorial coordinates
    norm_node_omega = tuple(w - omega_mean for w in raw_node_omega)

    # 6. Segment-level warping and C_w calculation
    # Reconstruct segment C_w integration using local segment maximum w_max_i
    # and math.frexp / math.ldexp decomposition to eliminate intermediate underflow/overflow
    # even under heterogeneous multiscale geometry.
    m_cw_list: list[float] = []
    e_cw_list: list[int] = []

    for s_idx, seg in enumerate(segments):
        node_u, node_v = canonical_edges[s_idx]
        w1 = norm_node_omega[node_u]
        w2 = norm_node_omega[node_v]
        w_max_i = max(abs(w1), abs(w2))

        if w_max_i == 0.0:
            m_cw_list.append(0.0)
            e_cw_list.append(0)
            continue

        w1_norm = w1 / w_max_i
        w2_norm = w2 / w_max_i

        # Complete squares formulation on normalized endpoint values [-1, 1]:
        # w1_norm^2 + w1_norm * w2_norm + w2_norm^2 = (w1_norm + 0.5 * w2_norm)^2 + 0.75 * w2_norm^2
        u_c = w1_norm + 0.5 * w2_norm
        v_c = (math.sqrt(3.0) / 2.0) * w2_norm
        K_hat = u_c * u_c + v_c * v_c

        # C_{w,i} = (t_i * L_i * w_max_i^2 / 3.0) * K_hat
        mt, et = math.frexp(seg.t)
        mL, eL = math.frexp(seg.length)
        mw, ew = math.frexp(w_max_i)
        mK, eK = math.frexp(K_hat / 3.0)

        m_i = mt * mL * (mw * mw) * mK
        e_i = et + eL + 2 * ew + eK

        m_norm, e_norm = math.frexp(m_i)
        m_cw_list.append(m_norm)
        e_cw_list.append(e_i + e_norm)

    non_zero_e = [e for m, e in zip(m_cw_list, e_cw_list) if m > 0.0]
    if not non_zero_e:
        Cw = 0.0
        cw_contributions = [0.0] * len(segments)
    else:
        e_cw_max = max(non_zero_e)
        scaled_cw_mantissas: list[float] = []
        for m, e in zip(m_cw_list, e_cw_list):
            if m == 0.0:
                scaled_cw_mantissas.append(0.0)
            else:
                shift = e - e_cw_max
                if shift < -1100:
                    scaled_cw_mantissas.append(0.0)
                else:
                    scaled_cw_mantissas.append(math.ldexp(m, shift))

        total_cw_mantissa = math.fsum(scaled_cw_mantissas)
        try:
            Cw = math.ldexp(total_cw_mantissa, e_cw_max)
        except OverflowError:
            Cw = float("inf")

        cw_contributions = []
        for m, e in zip(m_cw_list, e_cw_list):
            if m == 0.0:
                cw_contributions.append(0.0)
            else:
                try:
                    val = math.ldexp(m, e)
                except OverflowError:
                    val = float("inf")
                cw_contributions.append(val)

    segment_warpings: list[SegmentWarping] = []
    for s_idx, seg in enumerate(segments):
        node_u, node_v = canonical_edges[s_idx]
        w1 = norm_node_omega[node_u]
        w2 = norm_node_omega[node_v]
        segment_warpings.append(
            SegmentWarping(
                segment=seg,
                omega1=float(w1),
                omega2=float(w2),
                p=seg_rates[s_idx],
                cw_segment=cw_contributions[s_idx],
                j_segment=j_contributions[s_idx],
            )
        )

    if not math.isfinite(Cw):
        raise GeometryError(f"Non-finite warping constant C_w computed: {Cw}.")
    if Cw < 0.0:
        raise GeometryError(f"Warping constant C_w cannot be negative. Got {Cw}.")

    return TorsionWarpingResult(
        section=section,
        J=J,
        Cw=Cw,
        omega_mean=omega_mean,
        shear_center=(xs, ys),
        node_omega=norm_node_omega,
        segment_warpings=tuple(segment_warpings),
        raw_node_omega=tuple(raw_node_omega),
    )
