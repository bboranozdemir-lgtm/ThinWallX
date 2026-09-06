"""Closed-section Bredt-Batho Saint-Venant torsion and warping mechanics (ThinWallX v0.5).

Implements:
1. Direct solve of Bredt-Batho multi-cell torsion equations:
       H * phi = 2 * A
   where phi = q_cell / (G * theta') are cell circulations per unit G * theta'.
2. Membrane shear-flow coefficient vector:
       F = B^T * phi  ==>  q_i^(T) = G * theta' * F_i
3. Bredt-Batho torsion constant:
       J_BB = 2 * A^T * phi
   (No open-strip term sum_i (L_i * t_i^3 / 3) is added to closed sections).
4. Closed-section shear center from exact centroidal torques under basis transverse shears:
       e_x = T_z^(y),  e_y = -T_z^(x)
5. Sectorial coordinate rate including membrane-flow correction:
       g_i = d(omega)/ds = p_i - F_i / t_i
   with shear center as pole: p_i = [(r_1 - r_S) x t_i]_z.
6. Cycle consistency identity:
       \\oint_c d(omega) = 2 * A_c - (H * phi)_c = 0
7. Graph propagation, zero-mean normalization int_A omega dA = 0, and stable C_w integration.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
import math
from typing import TYPE_CHECKING, Sequence

import numpy as np

from thinwallx.exceptions import GeometryError
from thinwallx.primitives import Segment
from thinwallx.shear_center import ShearCenterResult
from thinwallx.torsion import SegmentWarping, TorsionWarpingResult

if TYPE_CHECKING:
    from thinwallx.cells import CellTopology
    from thinwallx.section import Section


@dataclass(frozen=True)
class ClosedTorsionWarpingResult(TorsionWarpingResult):
    """Complete results of closed-section torsion and warping analysis.

    Extends TorsionWarpingResult with Bredt-Batho multi-cell quantities:
        phi: Solved cell circulation vector per unit G*theta' of shape (n_c,).
        F: Wall membrane shear-flow coefficient vector of shape (E,).
        cell_areas: Enclosed median-line cell areas of shape (n_c,).
        cycle_residuals: Residuals of \\oint_c d(omega) for each cell.
    """

    phi: np.ndarray = None  # type: ignore[assignment]
    F: np.ndarray = None  # type: ignore[assignment]
    cell_areas: np.ndarray = None  # type: ignore[assignment]
    cycle_residuals: np.ndarray = None  # type: ignore[assignment]


def compute_closed_shear_center(
    section: Section,
    safety_factor: float = 1e4,
) -> ShearCenterResult:
    """Compute exact shear center for a closed thin-walled section.

    Solves unit +x and +y transverse shear cases with zero-twist cell compatibility,
    and integrates the exact centroidal torque T_z of the final compatible flow:
        e_y = -T_z^(x)
        e_x =  T_z^(y)
        x_s = cx + e_x,  y_s = cy + e_y

    Args:
        section: Validated ClosedSection (or cellular Section).
        safety_factor: Multiplier for positive-definiteness check.

    Returns:
        ShearCenterResult containing offsets (ex, ey) and coordinates (x_s, y_s).
    """
    from thinwallx.closed_shear_flow import calculate_closed_shear_flow

    # Basis load 1: unit shear in +x direction [1, 0]^T
    res_x = calculate_closed_shear_flow(
        section=section,
        vx=1.0,
        vy=0.0,
        safety_factor=safety_factor,
    )
    tz_x = res_x.torque
    ey = -tz_x

    # Basis load 2: unit shear in +y direction [0, 1]^T
    res_y = calculate_closed_shear_flow(
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


def compute_closed_torsion_warping(
    section: Section,
    root_node_idx: int = 0,
    safety_factor: float = 1e4,
) -> ClosedTorsionWarpingResult:
    """Compute exact closed-section Bredt-Batho torsion constant J and warping constant C_w.

    Args:
        section: Validated ClosedSection (or cellular Section).
        root_node_idx: Canonical node index used as root for raw omega propagation.
        safety_factor: Multiplier for numerical positive-definiteness of H.

    Returns:
        ClosedTorsionWarpingResult containing J_BB, C_w, normalized omega field, and Bredt quantities.
    """
    from thinwallx.cells import CellTopology, extract_cell_topology

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
    v_count = len(canonical_nodes)
    e_count = len(section.segments)

    if root_node_idx < 0 or root_node_idx >= v_count:
        raise GeometryError(
            f"root_node_idx must be in [0, {v_count - 1}]. Got {root_node_idx}."
        )

    # 1. Bredt-Batho cell compatibility: H * phi = 2 * A
    # Solved in dimensionless scaled coordinates: H_scaled * phi_scaled = 2 * A_scaled
    # to eliminate intermediate overflow/underflow across disparate scales
    H_scaled = cell_topo.H_scaled
    e_max_H = cell_topo.e_max_H

    A_vec = np.array([c.area for c in cell_topo.cells], dtype=float)
    mA_list = []
    eA_list = []
    for c in cell_topo.cells:
        mA, eA = math.frexp(c.area)
        mA_list.append(mA)
        eA_list.append(eA)
    e_max_A = max(eA_list)

    A_scaled = np.zeros(n_c, dtype=float)
    for c in range(n_c):
        delta_e = eA_list[c] - e_max_A
        if delta_e >= -1100:
            A_scaled[c] = math.ldexp(mA_list[c], delta_e)
        else:
            A_scaled[c] = 0.0

    phi_scaled = np.linalg.solve(H_scaled, 2.0 * A_scaled)
    if not np.all(np.isfinite(phi_scaled)):
        raise GeometryError(
            f"Non-finite scaled circulation vector phi_scaled in Bredt-Batho solve: {phi_scaled}."
        )

    # Recombine phi = 2^(e_max_A - e_max_H) * phi_scaled
    phi = np.zeros(n_c, dtype=float)
    for c in range(n_c):
        try:
            phi[c] = math.ldexp(float(phi_scaled[c]), e_max_A - e_max_H)
        except OverflowError:
            phi[c] = math.copysign(math.inf, phi_scaled[c])

    # 2. Bredt-Batho torsion constant: J_BB = 2 * A^T * phi
    # J_BB = 2 * (2^e_max_A * A_scaled)^T * (2^(e_max_A - e_max_H) * phi_scaled)
    #      = 2^(1 + 2*e_max_A - e_max_H) * (A_scaled^T * phi_scaled)
    s = float(np.dot(A_scaled, phi_scaled))
    if s <= 0.0 or not math.isfinite(s):
        raise GeometryError(f"Scaled torsion scalar must be positive finite. Got {s}.")
    ms, es = math.frexp(s)
    exp_J = (1 + 2 * e_max_A - e_max_H) + es
    try:
        J_bb = math.ldexp(ms, exp_J)
    except OverflowError:
        J_bb = math.inf

    if not math.isfinite(J_bb) or J_bb <= 0.0:
        raise GeometryError(
            f"Torsion constant J_BB must be finite and strictly positive. Got {J_bb}."
        )

    # 3. Membrane shear-flow coefficient vector: F = B^T * phi
    F_vec = B.T @ phi
    F_scaled = B.T @ phi_scaled

    # 4. Shear center pole S = (x_s, y_s) and centroid offset (ex, ey)
    sc_res = compute_closed_shear_center(section=section, safety_factor=safety_factor)
    ex, ey = sc_res.ex, sc_res.ey
    xs, ys = sc_res.x, sc_res.y

    # Local reference shift to avoid loss of precision under large translations (e.g. 1e12):
    ref_x = min(n.x for n in canonical_nodes)
    ref_y = min(n.y for n in canonical_nodes)
    cx_loc = sum(
        seg.area * ((seg.p1.x - ref_x) + (seg.p2.x - ref_x)) / 2.0
        for seg in section.segments
    ) / section.area
    cy_loc = sum(
        seg.area * ((seg.p1.y - ref_y) + (seg.p2.y - ref_y)) / 2.0
        for seg in section.segments
    ) / section.area
    rS_loc = (cx_loc + ex, cy_loc + ey)

    # 5. Segment sectorial rate: g_i = p_i - F_i / t_i
    # Lever arm p_i = [(r_{1,i} - r_S) x t_i]_z
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

        # Robust F_i / t_i using mantissa/exponent scaling:
        val_fs = float(F_scaled[i])
        if val_fs == 0.0:
            f_over_t = 0.0
        else:
            mt, et = math.frexp(t)
            m_div, e_div = math.frexp(val_fs / mt)
            exp_tot = (e_max_A - e_max_H - et) + e_div
            try:
                f_over_t = math.ldexp(m_div, exp_tot)
            except OverflowError:
                f_over_t = math.copysign(math.inf, m_div)

        g_rates[i] = pi - f_over_t

    # 6. Verify cycle consistency: \oint_c d(omega) = sum_i B_ci * g_i * L_i = 0
    L_vec = np.array([seg.length for seg in section.segments], dtype=float)
    cycle_residuals = B @ (g_rates * L_vec)

    # 7. Graph propagation of raw sectorial coordinates
    # Build undirected adjacency: adj[u] -> list of (v, seg_idx, direction_sign)
    adj: list[list[tuple[int, int, float]]] = [[] for _ in range(v_count)]
    for s_idx, (u, v) in enumerate(canonical_edges):
        # Traversing u -> v follows segment direction (+1)
        # Traversing v -> u opposes segment direction (-1)
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

    # 8. Zero-mean normalization: int_A omega dA = 0
    # Segment endpoints:
    total_area = section.area
    weighted_omega_sum = 0.0
    for s_idx, (u, v) in enumerate(canonical_edges):
        seg = section.segments[s_idx]
        w1 = omega_raw[u]
        w2 = omega_raw[v]
        # Average omega on straight segment is (w1 + w2) / 2
        weighted_omega_sum += seg.t * seg.length * 0.5 * (w1 + w2)

    omega_mean = weighted_omega_sum / total_area
    omega_norm = omega_raw - omega_mean

    # 9. Compute segment contributions and total warping constant C_w
    # C_w = sum_i (t_i * L_i / 3) * (w1^2 + w1*w2 + w2^2)
    # Using stable square completion:
    # a = w1 * sqrt(t*L/3), b = w2 * sqrt(t*L/3)
    # a^2 + ab + b^2 = (a + 0.5*b)^2 + 0.75*b^2
    cw_contributions: list[float] = []
    segment_warpings: list[SegmentWarping] = []

    for s_idx, (u, v) in enumerate(canonical_edges):
        seg = section.segments[s_idx]
        w1 = float(omega_norm[u])
        w2 = float(omega_norm[v])
        t = seg.t
        L = seg.length

        scale = math.sqrt(t) * math.sqrt(L / 3.0)
        a = w1 * scale
        b = w2 * scale
        u_sq = a + 0.5 * b
        v_sq = (math.sqrt(3.0) / 2.0) * b
        cw_seg = u_sq * u_sq + v_sq * v_sq

        cw_contributions.append(cw_seg)
        segment_warpings.append(
            SegmentWarping(
                segment=seg,
                omega1=w1,
                omega2=w2,
                p=float(g_rates[s_idx]),
                cw_segment=float(cw_seg),
                j_segment=0.0,  # Open-strip term is 0 for closed section Bredt-Batho mechanics
            )
        )

    Cw_total = math.fsum(cw_contributions)
    if not math.isfinite(Cw_total) or Cw_total < 0.0:
        raise GeometryError(
            f"Warping constant C_w must be finite and non-negative. Got {Cw_total}."
        )

    return ClosedTorsionWarpingResult(
        section=section,
        J=J_bb,
        Cw=float(Cw_total),
        omega_mean=float(omega_mean),
        shear_center=(float(xs), float(ys)),
        node_omega=tuple(float(w) for w in omega_norm),
        segment_warpings=tuple(segment_warpings),
        raw_node_omega=tuple(float(w) for w in omega_raw),
        phi=phi,
        F=F_vec,
        cell_areas=A_vec,
        cycle_residuals=cycle_residuals,
    )
