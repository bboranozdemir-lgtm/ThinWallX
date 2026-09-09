"""Mixed open-closed Saint-Venant torsion and warping mechanics.

Implements:
1. Hybrid Saint-Venant torsion constant:
       J_total = J_BB + J_open
   where J_BB = 4 A^T H^{-1} A for cyclic cells E_c, and
   J_open = sum_{e in E_o} (L_e * t_e^3 / 3) for open bridge branches only.
   (Closed cell walls are never double-counted in J_open).
2. Mixed shear center S from exact centroidal torques of final compatible transverse flows:
       e_x = T_z^(y),  e_y = -T_z^(x),  S = (cx + e_x, cy + e_y)
3. Continuous sectorial coordinate rate across junctions and bridges:
       d(omega)/ds = p_e - F_e / t_e  (e in E_c)
       d(omega)/ds = p_e               (e in E_o)
   with pole at mixed shear center S.
4. Cycle consistency: \\oint_c d(omega) = 2 A_c - (H f)_c = 0 for all cells.
5. Global graph propagation, zero-mean normalization int_A omega dA = 0, and stable C_w integration.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
import math
from typing import TYPE_CHECKING, Sequence

import numpy as np

from sectalix.exceptions import GeometryError
from sectalix.shear_center import ShearCenterResult
from sectalix.torsion import SegmentWarping, TorsionWarpingResult

if TYPE_CHECKING:
    from sectalix.mixed_topology import MixedTopology
    from sectalix.section import Section


@dataclass(frozen=True)
class MixedTorsionWarpingResult(TorsionWarpingResult):
    """Complete results of mixed open-closed torsion and warping analysis.

    Attributes:
        J_total: Total hybrid Saint-Venant torsion constant J_BB + J_open.
        J_BB: Bredt-Batho closed-cell torsion constant.
        J_open: Open-branch strip torsion constant sum_{E_o} (L*t^3/3).
        phi: Solved cell circulation vector per unit G*theta' of shape (n_c,).
        F: Wall membrane shear-flow coefficient vector of shape (E,).
        cell_areas: Enclosed median-line cell areas of shape (n_c,).
        cycle_residuals: Residuals of \\oint_c d(omega) for each cell.
        relative_residual: Dimensionless solver backward error.
        condition_number_estimate: Condition number estimate of H_scaled.
    """

    J_total: float = 0.0
    J_BB: float = 0.0
    J_open: float = 0.0
    phi: np.ndarray = None  # type: ignore[assignment]
    F: np.ndarray = None  # type: ignore[assignment]
    cell_areas: np.ndarray = None  # type: ignore[assignment]
    cycle_residuals: np.ndarray = None  # type: ignore[assignment]
    relative_residual: float = 0.0
    condition_number_estimate: float = 1.0


def _compute_reference_shift(coords: Sequence[float]) -> float:
    min_val = min(coords)
    max_val = max(coords)
    if min_val >= 0.0:
        return min_val
    if max_val <= 0.0:
        return max_val
    return 0.0


_NONZERO_UNDERFLOW_MESSAGE = (
    "Calculated value is non-zero but underflows IEEE-754 float64 "
    "subnormal range (< 5e-324)."
)


def _checked_ldexp(
    mantissa: float,
    exponent: int,
    *,
    overflow_message: str,
) -> float:
    """Reconstruct a non-zero scaled value without silent IEEE-754 loss."""
    if mantissa == 0.0:
        return 0.0
    try:
        value = math.ldexp(mantissa, exponent)
    except OverflowError as exc:
        raise OverflowError(overflow_message) from exc
    if not math.isfinite(value):
        raise OverflowError(overflow_message)
    if value == 0.0:
        raise FloatingPointError(_NONZERO_UNDERFLOW_MESSAGE)
    return value


def _sum_mantissa_exponent_terms(
    mantissas: Sequence[float],
    exponents: Sequence[int],
    *,
    overflow_message: str,
) -> float:
    """Sum non-negative mantissa/exponent terms with one checked reconstruction."""
    nonzero = [(m, e) for m, e in zip(mantissas, exponents) if m != 0.0]
    if not nonzero:
        return 0.0
    e_max = max(e for _, e in nonzero)
    aligned = [
        math.ldexp(m, e - e_max) if e - e_max >= -1100 else 0.0
        for m, e in nonzero
    ]
    total_mantissa = math.fsum(aligned)
    return _checked_ldexp(
        total_mantissa,
        e_max,
        overflow_message=overflow_message,
    )


def _compute_open_torsion_terms(
    section: Section,
    mixed_topo: MixedTopology,
) -> tuple[float, list[float]]:
    """Return J_open and representable per-segment strip contributions."""
    e_count = len(section.segments)
    contributions = [0.0] * e_count
    if not mixed_topo.open_edges:
        return 0.0, contributions

    mantissas: list[float] = []
    exponents: list[int] = []
    for s_idx in mixed_topo.open_edges:
        seg = section.segments[s_idx]
        m_length, e_length = math.frexp(seg.length)
        m_thickness, e_thickness = math.frexp(seg.t)
        m_term, e_norm = math.frexp((m_length * m_thickness**3) / 3.0)
        e_term = e_length + 3 * e_thickness + e_norm
        mantissas.append(m_term)
        exponents.append(e_term)
        try:
            contributions[s_idx] = math.ldexp(m_term, e_term)
        except OverflowError:
            contributions[s_idx] = math.inf

    total = _sum_mantissa_exponent_terms(
        mantissas,
        exponents,
        overflow_message="Open-section torsion constant J_open overflows IEEE-754 float64 range.",
    )
    return total, contributions


def _solve_bredt_batho_blocks(
    mixed_topo: MixedTopology,
    *,
    require_physical_fields: bool,
) -> tuple[
    float,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    float,
    float,
]:
    """Solve independent Bredt-Batho systems using each block's own exponent."""
    n_c = mixed_topo.cycle_rank
    e_count = mixed_topo.segment_count
    if n_c == 0:
        return (
            0.0,
            np.zeros(0, dtype=float),
            np.zeros(e_count, dtype=float),
            np.zeros(e_count, dtype=float),
            np.zeros(e_count, dtype=np.int64),
            0.0,
            1.0,
        )

    phi = np.zeros(n_c, dtype=float)
    F = np.zeros(e_count, dtype=float)
    F_scaled = np.zeros(e_count, dtype=float)
    F_exponents = np.zeros(e_count, dtype=np.int64)
    j_mantissas: list[float] = []
    j_exponents: list[int] = []
    max_relative_residual = 0.0
    max_condition = 1.0

    for cell_ids, H_block, e_max_H_block in zip(
        mixed_topo.block_cell_indices,
        mixed_topo.H_block_scaled,
        mixed_topo.e_max_H_block,
    ):
        m_areas: list[float] = []
        e_areas: list[int] = []
        for cell_idx in cell_ids:
            m_area, e_area = math.frexp(mixed_topo.cells[cell_idx].area)
            m_areas.append(m_area)
            e_areas.append(e_area)
        e_max_area = max(e_areas)
        A_scaled = np.array(
            [math.ldexp(m, e - e_max_area) for m, e in zip(m_areas, e_areas)],
            dtype=float,
        )

        phi_scaled = np.linalg.solve(H_block, 2.0 * A_scaled)
        if not np.all(np.isfinite(phi_scaled)):
            raise GeometryError(
                "Non-finite scaled circulation vector in a Bredt-Batho block solve: "
                f"{phi_scaled}."
            )

        residual = H_block @ phi_scaled - 2.0 * A_scaled
        norm_res = float(np.linalg.norm(residual, ord=np.inf))
        norm_H = float(np.linalg.norm(H_block, ord=np.inf))
        norm_phi = float(np.linalg.norm(phi_scaled, ord=np.inf))
        norm_rhs = float(np.linalg.norm(2.0 * A_scaled, ord=np.inf))
        denominator = norm_H * norm_phi + norm_rhs
        block_relative_residual = norm_res / denominator if denominator > 0.0 else 0.0
        max_relative_residual = max(max_relative_residual, block_relative_residual)

        eigvals = np.linalg.eigvalsh(H_block)
        block_condition = (
            float(eigvals[-1] / eigvals[0]) if eigvals[0] > 0.0 else math.inf
        )
        max_condition = max(max_condition, block_condition)

        scalar = float(np.dot(A_scaled, phi_scaled))
        if scalar <= 0.0 or not math.isfinite(scalar):
            raise GeometryError(
                f"Scaled torsion scalar must be positive finite. Got {scalar}."
            )
        m_scalar, e_scalar = math.frexp(scalar)
        j_mantissas.append(m_scalar)
        j_exponents.append(1 + 2 * e_max_area - e_max_H_block + e_scalar)

        phi_exponent = e_max_area - e_max_H_block
        block_indices = np.array(cell_ids, dtype=int)
        for local_idx, global_idx in enumerate(cell_ids):
            try:
                physical_phi = math.ldexp(float(phi_scaled[local_idx]), phi_exponent)
            except OverflowError as exc:
                if require_physical_fields:
                    raise OverflowError(
                        "Cell circulation phi overflows IEEE-754 float64 range."
                    ) from exc
                physical_phi = math.copysign(math.inf, float(phi_scaled[local_idx]))
            if require_physical_fields and not math.isfinite(physical_phi):
                raise OverflowError(
                    "Cell circulation phi overflows IEEE-754 float64 range."
                )
            phi[global_idx] = physical_phi

        B_block = mixed_topo.B[block_indices, :]
        F_block_scaled = B_block.T @ phi_scaled
        for edge_idx, scaled_value in enumerate(F_block_scaled):
            if scaled_value == 0.0:
                continue
            if F_scaled[edge_idx] != 0.0:
                raise GeometryError(
                    f"Closed edge {edge_idx} was assigned to more than one cyclic block."
                )
            F_scaled[edge_idx] = float(scaled_value)
            F_exponents[edge_idx] = phi_exponent
            try:
                physical_F = math.ldexp(float(scaled_value), phi_exponent)
            except OverflowError as exc:
                if require_physical_fields:
                    raise OverflowError(
                        "Wall membrane flow coefficient F overflows IEEE-754 float64 range."
                    ) from exc
                physical_F = math.copysign(math.inf, float(scaled_value))
            if require_physical_fields and not math.isfinite(physical_F):
                raise OverflowError(
                    "Wall membrane flow coefficient F overflows IEEE-754 float64 range."
                )
            F[edge_idx] = physical_F

    J_bb = _sum_mantissa_exponent_terms(
        j_mantissas,
        j_exponents,
        overflow_message="Closed-cell torsion constant J_BB overflows IEEE-754 float64 range.",
    )
    if J_bb <= 0.0:
        raise GeometryError(
            f"Torsion constant J_BB must be strictly positive. Got {J_bb}."
        )
    return (
        J_bb,
        phi,
        F,
        F_scaled,
        F_exponents,
        max_relative_residual,
        max_condition,
    )


def compute_mixed_shear_center(
    section: Section,
    safety_factor: float = 1e4,
) -> ShearCenterResult:
    """Compute exact shear center for an arbitrary connected mixed open-closed section.

    Solves unit +x and +y transverse shear cases with zero-twist cell compatibility,
    and integrates the exact centroidal torque T_z of the final compatible flow:
        e_y = -T_z^(x)
        e_x =  T_z^(y)
        x_s = cx + e_x,  y_s = cy + e_y

    Args:
        section: Validated Section or MixedSection.
        safety_factor: Multiplier for positive-definiteness check.

    Returns:
        ShearCenterResult containing offsets (ex, ey) and coordinates (x_s, y_s).
    """
    from sectalix.mixed_shear_flow import calculate_mixed_shear_flow

    res_x = calculate_mixed_shear_flow(
        section=section,
        vx=1.0,
        vy=0.0,
        safety_factor=safety_factor,
    )
    tz_x = res_x.torque
    ey = -tz_x

    res_y = calculate_mixed_shear_flow(
        section=section,
        vx=0.0,
        vy=1.0,
        safety_factor=safety_factor,
    )
    tz_y = res_y.torque
    ex = tz_y

    cx, cy = section.centroid
    xs = cx + ex
    ys = cy + ey

    return ShearCenterResult(
        section=section,
        ex=float(ex),
        ey=float(ey),
        x=float(xs),
        y=float(ys),
        tz_x=float(tz_x),
        tz_y=float(tz_y),
        centroid=(float(cx), float(cy)),
        flow_x=res_x,  # type: ignore[arg-type]
        flow_y=res_y,  # type: ignore[arg-type]
    )


def compute_mixed_torsion_constant(
    section: Section,
    safety_factor: float = 1e4,
) -> tuple[float, float, float]:
    """Compute hybrid Saint-Venant torsion constants independently of shear center / warping.

    Returns:
        tuple (J_total, J_BB, J_open)
    """
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

    J_open, _ = _compute_open_torsion_terms(section, mixed_topo)
    J_bb, *_ = _solve_bredt_batho_blocks(
        mixed_topo,
        require_physical_fields=False,
    )
    J_total = float(J_bb + J_open)
    if not math.isfinite(J_total):
        raise OverflowError(
            "Total torsion constant J_total overflows IEEE-754 float64 range."
        )
    return (float(J_total), float(J_bb), float(J_open))


def compute_mixed_torsion_warping(
    section: Section,
    root_node_idx: int = 0,
    safety_factor: float = 1e4,
) -> MixedTorsionWarpingResult:
    """Compute exact mixed-section hybrid torsion constant J and warping constant C_w.

    Args:
        section: Validated Section or MixedSection.
        root_node_idx: Canonical node index used as root for raw omega propagation.
        safety_factor: Multiplier for numerical positive-definiteness of H.

    Returns:
        MixedTorsionWarpingResult containing J_total, J_BB, J_open, C_w, normalized omega, and F.
    """
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
    n_c = mixed_topo.cycle_rank
    v_count = len(canonical_nodes)
    e_count = len(section.segments)

    if root_node_idx < 0 or root_node_idx >= v_count:
        raise GeometryError(
            f"root_node_idx must be in [0, {v_count - 1}]. Got {root_node_idx}."
        )

    # 1. Hybrid Saint-Venant torsion. Open strips and each independent cyclic
    # block retain their own binary scale through their final reconstruction.
    J_open, j_seg_list = _compute_open_torsion_terms(section, mixed_topo)
    (
        J_bb,
        phi,
        F_vec,
        F_scaled,
        F_exponents,
        relative_residual,
        cond_h,
    ) = _solve_bredt_batho_blocks(
        mixed_topo,
        require_physical_fields=True,
    )
    A_vec = np.array([c.area for c in mixed_topo.cells], dtype=float)
    cycle_residuals = np.zeros(n_c, dtype=float)

    J_total = float(J_bb + J_open)
    if not math.isfinite(J_total):
        raise OverflowError(
            "Total torsion constant J_total overflows IEEE-754 float64 range."
        )

    # 3. Shear center pole S = (x_s, y_s) and centroid offset (ex, ey)
    sc_res = compute_mixed_shear_center(section=section, safety_factor=safety_factor)
    ex, ey = sc_res.ex, sc_res.ey
    xs, ys = sc_res.x, sc_res.y

    # Local reference shift to eliminate cancellation under large translations (e.g. 1e12)
    # Uses 0.0 if section straddles origin (e.g. 1x1 core with 2^300 arms) to avoid destroying core coordinates
    ref_x = _compute_reference_shift([n.x for n in canonical_nodes])
    ref_y = _compute_reference_shift([n.y for n in canonical_nodes])
    cx_loc = sum(
        seg.area * ((seg.p1.x - ref_x) + (seg.p2.x - ref_x)) / 2.0
        for seg in section.segments
    ) / section.area
    cy_loc = sum(
        seg.area * ((seg.p1.y - ref_y) + (seg.p2.y - ref_y)) / 2.0
        for seg in section.segments
    ) / section.area
    rS_loc = (cx_loc + ex, cy_loc + ey)

    # 4. Segment sectorial rates g_i
    # g_i = p_i - F_i / t_i  (for closed edges E_c)
    # g_i = p_i               (for open edges E_o, since F_i = 0)
    p_rates = np.zeros(e_count, dtype=float)
    g_rates = np.zeros(e_count, dtype=float)

    for i, seg in enumerate(section.segments):
        L = seg.length
        t = seg.t
        tx = (seg.p2.x - seg.p1.x) / L
        ty = (seg.p2.y - seg.p1.y) / L

        r1_loc_x = seg.p1.x - ref_x
        r1_loc_y = seg.p1.y - ref_y
        rx = r1_loc_x - rS_loc[0]
        ry = r1_loc_y - rS_loc[1]
        pi = rx * ty - ry * tx
        p_rates[i] = pi

        if i in mixed_topo.closed_edges:
            val_fs = float(F_scaled[i])
            if val_fs == 0.0:
                f_over_t = 0.0
            else:
                mt, et = math.frexp(t)
                m_div, e_div = math.frexp(val_fs / mt)
                exp_tot = (int(F_exponents[i]) - et) + e_div
                try:
                    f_over_t = math.ldexp(m_div, exp_tot)
                except OverflowError:
                    f_over_t = math.copysign(math.inf, m_div)
            g_rates[i] = pi - f_over_t
        else:
            g_rates[i] = pi

    # 5. Cycle consistency verification on all cells
    # Evaluated with local cell reference corner to remain completely immune to distant arms
    cycle_residuals = np.zeros(n_c, dtype=float)
    if n_c > 0:
        for c_idx, cell in enumerate(mixed_topo.cells):
            cell_nodes = [canonical_nodes[u] for u in cell.node_indices]
            ref_cx = min(n.x for n in cell_nodes)
            ref_cy = min(n.y for n in cell_nodes)

            cell_p_sum = 0.0
            cell_f_sum = 0.0
            for s_idx, ori in zip(cell.segment_indices, cell.segment_orientations):
                seg = section.segments[s_idx]
                L = seg.length
                tx = (seg.p2.x - seg.p1.x) / L
                ty = (seg.p2.y - seg.p1.y) / L
                r1x = seg.p1.x - ref_cx
                r1y = seg.p1.y - ref_cy
                p_local = r1x * ty - r1y * tx
                cell_p_sum += ori * p_local * L

                val_fs = float(F_scaled[s_idx])
                if val_fs != 0.0:
                    mt, et = math.frexp(seg.t)
                    m_div, e_div = math.frexp(val_fs / mt)
                    exp_tot = (int(F_exponents[s_idx]) - et) + e_div
                    try:
                        f_over_t = math.ldexp(m_div, exp_tot)
                    except OverflowError:
                        f_over_t = math.copysign(math.inf, m_div)
                else:
                    f_over_t = 0.0
                cell_f_sum += ori * f_over_t * L

            res_c = cell_p_sum - cell_f_sum
            cycle_residuals[c_idx] = res_c
            if res_c == 0.0:
                rel_res = 0.0
            else:
                m_residual, e_residual = math.frexp(abs(res_c))
                m_area, e_area = math.frexp(cell.area)
                try:
                    rel_res = math.ldexp(
                        m_residual / (2.0 * m_area),
                        e_residual - e_area,
                    )
                except OverflowError:
                    rel_res = math.inf
            if rel_res > 1e-10:
                raise GeometryError(
                    f"Cycle warping consistency residual too large for cell {c_idx}: {rel_res:.6e}."
                )

    # 6. Global graph propagation of raw sectorial coordinates
    # Undirected adjacency over the full connected graph G
    adj: list[list[tuple[int, int, float]]] = [[] for _ in range(v_count)]
    for s_idx, (u, v) in enumerate(canonical_edges):
        adj[u].append((v, s_idx, 1.0))
        adj[v].append((u, s_idx, -1.0))

    omega_raw = np.full(v_count, np.nan, dtype=float)
    omega_raw[root_node_idx] = 0.0

    queue = deque([root_node_idx])
    visited = {root_node_idx}

    while queue:
        curr = queue.popleft()
        for neighbor, s_idx, sign in adj[curr]:
            if neighbor not in visited:
                visited.add(neighbor)
                d_omega = sign * g_rates[s_idx] * section.segments[s_idx].length
                omega_raw[neighbor] = omega_raw[curr] + d_omega
                queue.append(neighbor)

    if not np.all(np.isfinite(omega_raw)):
        raise GeometryError("Non-finite raw sectorial coordinate encountered during propagation.")

    # 7. Zero-mean normalization: int_A omega dA = 0
    total_area = section.area
    weighted_omega_sum = 0.0
    for s_idx, (u, v) in enumerate(canonical_edges):
        seg = section.segments[s_idx]
        w1 = omega_raw[u]
        w2 = omega_raw[v]
        weighted_omega_sum += seg.t * seg.length * 0.5 * (w1 + w2)

    omega_mean = weighted_omega_sum / total_area
    omega_norm = omega_raw - omega_mean

    # 8. Exact analytical C_w integration with exponent-aligned mantissa addition
    # Eliminates intermediate underflow to preserve microscopic / subnormal values down to 5e-324
    m_cw_list: list[float] = []
    e_cw_list: list[int] = []

    for s_idx, (u, v) in enumerate(canonical_edges):
        seg = section.segments[s_idx]
        w1 = float(omega_norm[u])
        w2 = float(omega_norm[v])
        w_max_i = max(abs(w1), abs(w2))

        if w_max_i == 0.0:
            m_cw_list.append(0.0)
            e_cw_list.append(0)
            continue

        w1_norm = w1 / w_max_i
        w2_norm = w2 / w_max_i

        u_c = w1_norm + 0.5 * w2_norm
        v_c = (math.sqrt(3.0) / 2.0) * w2_norm
        K_hat = u_c * u_c + v_c * v_c

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
        Cw_total = 0.0
        cw_contributions = [0.0] * e_count
    else:
        Cw_total = _sum_mantissa_exponent_terms(
            m_cw_list,
            e_cw_list,
            overflow_message="Warping constant Cw overflows IEEE-754 float64 range.",
        )

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

    if not math.isfinite(Cw_total):
        raise OverflowError(
            "Warping constant Cw overflows IEEE-754 float64 range."
        )
    if Cw_total < 0.0:
        raise GeometryError(
            f"Warping constant C_w must be non-negative. Got {Cw_total}."
        )

    segment_warpings: list[SegmentWarping] = []
    for s_idx, (u, v) in enumerate(canonical_edges):
        seg = section.segments[s_idx]
        w1 = float(omega_norm[u])
        w2 = float(omega_norm[v])
        segment_warpings.append(
            SegmentWarping(
                segment=seg,
                omega1=w1,
                omega2=w2,
                p=float(g_rates[s_idx]),
                cw_segment=float(cw_contributions[s_idx]),
                j_segment=float(j_seg_list[s_idx]),
            )
        )

    return MixedTorsionWarpingResult(
        section=section,
        J=float(J_total),
        Cw=float(Cw_total),
        omega_mean=float(omega_mean),
        shear_center=(float(xs), float(ys)),
        node_omega=tuple(float(w) for w in omega_norm),
        segment_warpings=tuple(segment_warpings),
        raw_node_omega=tuple(float(w) for w in omega_raw),
        J_total=float(J_total),
        J_BB=float(J_bb),
        J_open=float(J_open),
        phi=phi,
        F=F_vec,
        cell_areas=A_vec,
        cycle_residuals=cycle_residuals,
        relative_residual=relative_residual,
        condition_number_estimate=cond_h,
    )
