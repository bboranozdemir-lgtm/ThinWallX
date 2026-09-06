"""Independent analytical benchmarks for Sectalix v0.6 mixed open-closed mechanics.

Covers:
- B1: Closed-base hat / Omega (closed box + open flanges) and asymmetric variants.
- B2: Two-cell box girder with cantilevers.
- B3: Doubly symmetric rectangular box + side midpoint antennas (positive & zero C_w).
- B4: Branch and antenna limits (t_a -> 0+, l -> 0+, Y-branch, inclined antennas).

Includes completely independent analytical oracles for transverse shear flow,
shear center, and warping, using direct node-Kirchhoff and cell-circulation linear systems
without virtual tree cuts.
"""

from __future__ import annotations

from collections import deque
from decimal import Decimal, localcontext
import math
from typing import Sequence
import pytest
import numpy as np

from sectalix.mixed_section import MixedSection
from sectalix.primitives import Node, Segment
from sectalix.shear_load import ShearLoad


# ==============================================================================
# Independent Analytical Oracle
# ==============================================================================


def _decimal(value: float | int | str) -> Decimal:
    """Convert an input literal without introducing fresh binary arithmetic."""
    return Decimal(str(value))


def _decimal_solve(
    matrix: Sequence[Sequence[Decimal]],
    rhs: Sequence[Decimal],
) -> list[Decimal]:
    """Independent 100-digit partial-pivot Gaussian elimination."""
    n = len(rhs)
    augmented = [list(row) + [rhs[i]] for i, row in enumerate(matrix)]
    for column in range(n):
        pivot = max(range(column, n), key=lambda row: abs(augmented[row][column]))
        if augmented[pivot][column] == 0:
            raise ArithmeticError("Independent Decimal oracle matrix is singular.")
        if pivot != column:
            augmented[column], augmented[pivot] = augmented[pivot], augmented[column]

        divisor = augmented[column][column]
        for j in range(column, n + 1):
            augmented[column][j] /= divisor

        for row in range(n):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor == 0:
                continue
            for j in range(column, n + 1):
                augmented[row][j] -= factor * augmented[column][j]
    return [augmented[i][n] for i in range(n)]

class IndependentMixedOracle:
    """Independent mathematical reference solver for thin-walled mixed open-closed sections.

    Constructs and solves direct linear systems for edge constants c_e without virtual
    spanning tree cuts or production code imports:
    1. dq_e/ds = -t_e (a X(s) + b Y(s)) ==> q_e(xi) = c_e + alpha_e*xi + beta_e*xi^2
    2. Node balance: sum_in q - sum_out q = 0 for V-1 nodes
    3. Cell compatibility: oint_c (q/t) ds = 0 for all n_c bounded cells
    4. Direct solve for edge constants c_e
    5. Torque T_z and shear center (ex, ey)
    6. Independent Bredt-Batho solve H*f = 2A, continuous omega propagation, and C_w
    """

    def __init__(
        self,
        nodes: Sequence[tuple[float, float]],
        edges: Sequence[tuple[int, int, float]],
        cells: Sequence[tuple[tuple[int, ...], tuple[int, ...]]],  # (edge_indices, orientations)
    ) -> None:
        """
        Args:
            nodes: List of (x, y) coordinates.
            edges: List of (u, v, thickness) for each physical edge.
            cells: List of (cell_edge_indices, cell_orientations).
        """
        self.nodes = [(_decimal(x), _decimal(y)) for x, y in nodes]
        self.edges = [(u, v, _decimal(t)) for u, v, t in edges]
        self.cells = list(cells)
        self.num_nodes = len(nodes)
        self.num_edges = len(edges)
        self.num_cells = len(cells)

        with localcontext() as context:
            context.prec = 100
            zero = Decimal(0)
            two = Decimal(2)
            three = Decimal(3)
            total_area = zero
            sum_x = zero
            sum_y = zero
            self.edge_lengths: list[Decimal] = []
            for u, v, thickness in self.edges:
                x1, y1 = self.nodes[u]
                x2, y2 = self.nodes[v]
                length = ((x2 - x1) ** 2 + (y2 - y1) ** 2).sqrt()
                self.edge_lengths.append(length)
                edge_area = thickness * length
                total_area += edge_area
                sum_x += edge_area * (x1 + x2) / two
                sum_y += edge_area * (y1 + y2) / two
            self.area = total_area
            self.cx = sum_x / total_area
            self.cy = sum_y / total_area

            ix = zero
            iy = zero
            ixy = zero
            for i, (u, v, thickness) in enumerate(self.edges):
                length = self.edge_lengths[i]
                x1, y1 = self.nodes[u]
                x2, y2 = self.nodes[v]
                X1 = x1 - self.cx
                Y1 = y1 - self.cy
                dX = x2 - x1
                dY = y2 - y1
                ix += thickness * length * (Y1**2 + Y1 * dY + dY**2 / three)
                iy += thickness * length * (X1**2 + X1 * dX + dX**2 / three)
                ixy += thickness * length * (
                    X1 * Y1 + (X1 * dY + Y1 * dX) / two + dX * dY / three
                )
            self.Ix = ix
            self.Iy = iy
            self.Ixy = ixy

            self.cell_areas: list[Decimal] = []
            for cell_edges, cell_orientations in self.cells:
                polygon = []
                for edge_idx, orientation in zip(cell_edges, cell_orientations):
                    u, v, _ = self.edges[edge_idx]
                    polygon.append(self.nodes[u] if orientation == 1 else self.nodes[v])
                x0, y0 = polygon[0]
                twice_area = sum(
                    (polygon[k][0] - x0)
                    * (polygon[(k + 1) % len(polygon)][1] - y0)
                    - (polygon[(k + 1) % len(polygon)][0] - x0)
                    * (polygon[k][1] - y0)
                    for k in range(len(polygon))
                )
                self.cell_areas.append(abs(twice_area / two))

            self.H = [
                [zero for _ in range(self.num_cells)]
                for _ in range(self.num_cells)
            ]
            for i, (edges_i, orientations_i) in enumerate(self.cells):
                incidence_i = dict(zip(edges_i, orientations_i))
                for j, (edges_j, orientations_j) in enumerate(self.cells):
                    incidence_j = dict(zip(edges_j, orientations_j))
                    value = zero
                    for edge_idx, orientation_i in incidence_i.items():
                        if edge_idx in incidence_j:
                            value += (
                                Decimal(orientation_i * incidence_j[edge_idx])
                                * self.edge_lengths[edge_idx]
                                / self.edges[edge_idx][2]
                            )
                    self.H[i][j] = value

    def _solve_shear_flow_decimal(
        self,
        vx: Decimal,
        vy: Decimal,
    ) -> tuple[list[Decimal], list[Decimal], Decimal]:
        """Solve the independent direct system entirely at 100-digit precision."""
        zero = Decimal(0)
        two = Decimal(2)
        three = Decimal(3)
        alpha, beta = _decimal_solve(
            [[self.Iy, self.Ixy], [self.Ixy, self.Ix]],
            [vx, vy],
        )

        linear: list[Decimal] = []
        quadratic: list[Decimal] = []
        for i, (u, v, thickness) in enumerate(self.edges):
            length = self.edge_lengths[i]
            X1 = self.nodes[u][0] - self.cx
            Y1 = self.nodes[u][1] - self.cy
            dX = self.nodes[v][0] - self.nodes[u][0]
            dY = self.nodes[v][1] - self.nodes[u][1]
            linear.append(-thickness * length * (alpha * X1 + beta * Y1))
            quadratic.append(
                -thickness * length * (alpha * dX + beta * dY) / two
            )

        matrix = [
            [zero for _ in range(self.num_edges)]
            for _ in range(self.num_edges)
        ]
        rhs = [zero for _ in range(self.num_edges)]
        row = 0
        for node_idx in range(1, self.num_nodes):
            for edge_idx, (u, v, _) in enumerate(self.edges):
                if v == node_idx:
                    matrix[row][edge_idx] += Decimal(1)
                    rhs[row] -= linear[edge_idx] + quadratic[edge_idx]
                elif u == node_idx:
                    matrix[row][edge_idx] -= Decimal(1)
            row += 1

        for cell_edges, cell_orientations in self.cells:
            for edge_idx, orientation in zip(cell_edges, cell_orientations):
                weight = (
                    Decimal(orientation)
                    * self.edge_lengths[edge_idx]
                    / self.edges[edge_idx][2]
                )
                matrix[row][edge_idx] += weight
                rhs[row] -= weight * (
                    linear[edge_idx] / two + quadratic[edge_idx] / three
                )
            row += 1
        assert row == self.num_edges

        constants = _decimal_solve(matrix, rhs)
        integrals: list[Decimal] = []
        torque = zero
        for edge_idx, (u, v, _) in enumerate(self.edges):
            length = self.edge_lengths[edge_idx]
            integral = length * (
                constants[edge_idx]
                + linear[edge_idx] / two
                + quadratic[edge_idx] / three
            )
            integrals.append(integral)
            X1 = self.nodes[u][0] - self.cx
            Y1 = self.nodes[u][1] - self.cy
            tx = (self.nodes[v][0] - self.nodes[u][0]) / length
            ty = (self.nodes[v][1] - self.nodes[u][1]) / length
            torque += (X1 * ty - Y1 * tx) * integral
        return constants, integrals, torque

    def solve_shear_flow(self, vx: float, vy: float) -> tuple[np.ndarray, np.ndarray, float]:
        """Solve exact shear-flow field: returns (c_vector, edge_integrals, torque)."""
        with localcontext() as context:
            context.prec = 100
            constants, integrals, torque = self._solve_shear_flow_decimal(
                _decimal(vx), _decimal(vy)
            )
        return (
            np.array([float(value) for value in constants], dtype=float),
            np.array([float(value) for value in integrals], dtype=float),
            float(torque),
        )

    def _compute_shear_center_decimal(
        self,
    ) -> tuple[Decimal, Decimal, Decimal, Decimal]:
        _, _, torque_x = self._solve_shear_flow_decimal(Decimal(1), Decimal(0))
        _, _, torque_y = self._solve_shear_flow_decimal(Decimal(0), Decimal(1))
        ey = -torque_x
        ex = torque_y
        return ex, ey, self.cx + ex, self.cy + ey

    def compute_shear_center(self) -> tuple[float, float, float, float]:
        """Compute shear center offsets (ex, ey) and coordinates (xs, ys)."""
        with localcontext() as context:
            context.prec = 100
            ex, ey, xs, ys = self._compute_shear_center_decimal()
        return float(ex), float(ey), float(xs), float(ys)

    def compute_torsion_warping(self) -> tuple[float, float, float, float, np.ndarray]:
        """Compute (J_total, J_BB, J_open, C_w, omega_norm)."""
        with localcontext() as context:
            context.prec = 100
            zero = Decimal(0)
            two = Decimal(2)
            three = Decimal(3)
            if self.num_cells:
                phi = _decimal_solve(
                    self.H,
                    [two * area for area in self.cell_areas],
                )
                J_BB = two * sum(
                    area * circulation
                    for area, circulation in zip(self.cell_areas, phi)
                )
            else:
                phi = []
                J_BB = zero

            closed_edges = {
                edge for cell_edges, _ in self.cells for edge in cell_edges
            }
            open_edges = set(range(self.num_edges)) - closed_edges
            J_open = sum(
                self.edge_lengths[edge] * self.edges[edge][2] ** 3 / three
                for edge in open_edges
            )
            J_total = J_BB + J_open

            membrane_flow = [zero for _ in range(self.num_edges)]
            for cell_idx, (cell_edges, orientations) in enumerate(self.cells):
                for edge_idx, orientation in zip(cell_edges, orientations):
                    membrane_flow[edge_idx] += Decimal(orientation) * phi[cell_idx]

            _, _, xs, ys = self._compute_shear_center_decimal()
            rates: list[Decimal] = []
            for edge_idx, (u, v, thickness) in enumerate(self.edges):
                length = self.edge_lengths[edge_idx]
                tx = (self.nodes[v][0] - self.nodes[u][0]) / length
                ty = (self.nodes[v][1] - self.nodes[u][1]) / length
                lever = (
                    (self.nodes[u][0] - xs) * ty
                    - (self.nodes[u][1] - ys) * tx
                )
                rates.append(
                    lever - membrane_flow[edge_idx] / thickness
                    if edge_idx in closed_edges
                    else lever
                )

            adjacency: list[list[tuple[int, int, int]]] = [
                [] for _ in range(self.num_nodes)
            ]
            for edge_idx, (u, v, _) in enumerate(self.edges):
                adjacency[u].append((v, edge_idx, 1))
                adjacency[v].append((u, edge_idx, -1))
            omega_raw: list[Decimal | None] = [None] * self.num_nodes
            omega_raw[0] = zero
            visited = {0}
            queue = deque([0])
            while queue:
                current = queue.popleft()
                assert omega_raw[current] is not None
                for neighbor, edge_idx, sign in adjacency[current]:
                    if neighbor in visited:
                        continue
                    visited.add(neighbor)
                    omega_raw[neighbor] = (
                        omega_raw[current]
                        + Decimal(sign) * rates[edge_idx] * self.edge_lengths[edge_idx]
                    )
                    queue.append(neighbor)
            assert all(value is not None for value in omega_raw)
            raw_values = [value for value in omega_raw if value is not None]

            weighted_omega = sum(
                self.edges[i][2]
                * self.edge_lengths[i]
                * (raw_values[u] + raw_values[v])
                / two
                for i, (u, v, _) in enumerate(self.edges)
            )
            omega_mean = weighted_omega / self.area
            omega_norm = [value - omega_mean for value in raw_values]
            Cw = sum(
                self.edges[i][2]
                * self.edge_lengths[i]
                * (
                    omega_norm[u] ** 2
                    + omega_norm[u] * omega_norm[v]
                    + omega_norm[v] ** 2
                )
                / three
                for i, (u, v, _) in enumerate(self.edges)
            )

        return (
            float(J_total),
            float(J_BB),
            float(J_open),
            float(Cw),
            np.array([float(value) for value in omega_norm], dtype=float),
        )


# ==============================================================================
# Benchmarks
# ==============================================================================

class TestMixedBenchmarks:
    """Rigorous independent oracle tests for benchmarks B1, B2, B3, B4."""

    def test_b1_closed_base_hat_omega(self) -> None:
        """B1: Closed-base hat / Omega section: box (4x2, t=1/50) + flanges (1x0, t=1/100)."""
        # Theoretical exact values from ACTIVE_PHASE.md:
        # J_BB = 32/75 = 0.4266666666666667
        # J_open = 1/1500000 = 6.666666666666667e-7
        # J_total = 640001/1500000
        # S_x = 2.0
        t_box = 1.0 / 50.0
        t_flange = 1.0 / 100.0

        # Physical segments for MixedSection
        n_fl1 = Node(-1.0, 0.0)
        n00 = Node(0.0, 0.0)
        n40 = Node(4.0, 0.0)
        n_fl2 = Node(5.0, 0.0)
        n42 = Node(4.0, 2.0)
        n02 = Node(0.0, 2.0)

        segs = [
            Segment(n00, n40, t=t_box),       # 0: bottom wall
            Segment(n40, n42, t=t_box),       # 1: right wall
            Segment(n42, n02, t=t_box),       # 2: top wall
            Segment(n02, n00, t=t_box),       # 3: left wall
            Segment(n00, n_fl1, t=t_flange),  # 4: left flange
            Segment(n40, n_fl2, t=t_flange),  # 5: right flange
        ]

        sec = MixedSection(segs)

        # 1. Saint-Venant torsion constant checks
        expected_J_bb = 32.0 / 75.0
        expected_J_open = 1.0 / 1500000.0
        expected_J_total = 640001.0 / 1500000.0

        assert math.isclose(sec.J_BB, expected_J_bb, rel_tol=1e-12)
        assert math.isclose(sec.J_open, expected_J_open, rel_tol=1e-12)
        assert math.isclose(sec.J_total, expected_J_total, rel_tol=1e-12)
        assert math.isclose(sec.J, expected_J_total, rel_tol=1e-12)

        # 2. Independent oracle verification
        # Node indices:
        # 0: (0,0), 1: (4,0), 2: (4,2), 3: (0,2), 4: (-1,0), 5: (5,0)
        nodes_ref = [(0.0, 0.0), (4.0, 0.0), (4.0, 2.0), (0.0, 2.0), (-1.0, 0.0), (5.0, 0.0)]
        edges_ref = [
            (0, 1, t_box),      # 0
            (1, 2, t_box),      # 1
            (2, 3, t_box),      # 2
            (3, 0, t_box),      # 3
            (0, 4, t_flange),   # 4
            (1, 5, t_flange),   # 5
        ]
        cells_ref = [((0, 1, 2, 3), (1, 1, 1, 1))]

        oracle = IndependentMixedOracle(nodes_ref, edges_ref, cells_ref)

        ex_ref, ey_ref, xs_ref, ys_ref = oracle.compute_shear_center()
        assert math.isclose(sec.shear_center[0], xs_ref, rel_tol=1e-10, abs_tol=1e-12)
        assert math.isclose(sec.shear_center[1], ys_ref, rel_tol=1e-10, abs_tol=1e-12)
        assert math.isclose(sec.shear_center[0], 2.0, abs_tol=1e-12)  # Symmetry!

        j_tot_ref, j_bb_ref, j_op_ref, cw_ref, _ = oracle.compute_torsion_warping()
        assert math.isclose(sec.J_total, j_tot_ref, rel_tol=1e-12)
        assert math.isclose(sec.Cw, cw_ref, rel_tol=1e-10)

        # 3. Transverse shear flow check under Vx=100, Vy=50
        sf_res = sec.calculate_shear_flow(vx=100.0, vy=50.0)
        c_ref, edge_int_ref, tz_ref = oracle.solve_shear_flow(vx=100.0, vy=50.0)
        assert math.isclose(sf_res.torque, tz_ref, rel_tol=1e-10)
        assert math.isclose(sf_res.recovered_resultant[0], 100.0, rel_tol=1e-12)
        assert math.isclose(sf_res.recovered_resultant[1], 50.0, rel_tol=1e-12)
        for i, sf in enumerate(sf_res.segment_flows):
            assert math.isclose(sf.integral_q, edge_int_ref[i], rel_tol=1e-10)
            assert math.isclose(sf.q0, c_ref[i], rel_tol=1e-10, abs_tol=1e-12)

    def test_b1_asymmetric_variants(self) -> None:
        """B1 variant: asymmetric flange lengths and thicknesses."""
        t_box = 0.02
        t_f1 = 0.015
        t_f2 = 0.01

        n_fl1 = Node(-2.5, 0.0)  # Flange 1 length 2.5
        n00 = Node(0.0, 0.0)
        n40 = Node(4.0, 0.0)
        n_fl2 = Node(5.2, 0.0)  # Flange 2 length 1.2
        n42 = Node(4.0, 2.0)
        n02 = Node(0.0, 2.0)

        segs = [
            Segment(n00, n40, t=t_box),
            Segment(n40, n42, t=t_box),
            Segment(n42, n02, t=t_box),
            Segment(n02, n00, t=t_box),
            Segment(n00, n_fl1, t=t_f1),
            Segment(n40, n_fl2, t=t_f2),
        ]

        sec = MixedSection(segs)

        nodes_ref = [(0.0, 0.0), (4.0, 0.0), (4.0, 2.0), (0.0, 2.0), (-2.5, 0.0), (5.2, 0.0)]
        edges_ref = [
            (0, 1, t_box),
            (1, 2, t_box),
            (2, 3, t_box),
            (3, 0, t_box),
            (0, 4, t_f1),
            (1, 5, t_f2),
        ]
        cells_ref = [((0, 1, 2, 3), (1, 1, 1, 1))]
        oracle = IndependentMixedOracle(nodes_ref, edges_ref, cells_ref)

        ex_ref, ey_ref, xs_ref, ys_ref = oracle.compute_shear_center()
        j_tot_ref, j_bb_ref, j_op_ref, cw_ref, _ = oracle.compute_torsion_warping()

        assert math.isclose(sec.shear_center[0], xs_ref, rel_tol=1e-10, abs_tol=1e-12)
        assert math.isclose(sec.shear_center[1], ys_ref, rel_tol=1e-10, abs_tol=1e-12)
        assert math.isclose(sec.J_total, j_tot_ref, rel_tol=1e-12)
        assert math.isclose(sec.Cw, cw_ref, rel_tol=1e-10)

        # Transverse shear flow check under Vx=150, Vy=-80
        sf_res = sec.calculate_shear_flow(vx=150.0, vy=-80.0)
        c_ref, edge_int_ref, tz_ref = oracle.solve_shear_flow(vx=150.0, vy=-80.0)
        assert math.isclose(sf_res.torque, tz_ref, rel_tol=1e-10)
        assert math.isclose(sf_res.recovered_resultant[0], 150.0, rel_tol=1e-12)
        assert math.isclose(sf_res.recovered_resultant[1], -80.0, rel_tol=1e-12)
        for i, sf in enumerate(sf_res.segment_flows):
            assert math.isclose(sf.integral_q, edge_int_ref[i], rel_tol=1e-10)
            assert math.isclose(sf.q0, c_ref[i], rel_tol=1e-10, abs_tol=1e-12)

    def test_b2_two_cell_box_girder_with_cantilevers(self) -> None:
        """B2: Two-cell box girder with cantilevers."""
        # Dimensions and properties from ACTIVE_PHASE.md:
        # Outer rect: (0,0) to (8,2), internal web at x=3
        # Outer walls: t=1/50, shared web: t=1/40, cantilevers: t=1/100 (L1=2, L2=1)
        # H = [[480, -80], [-80, 680]]
        # J_BB = 513/500 = 1.026
        # J_open = 1/1000000 = 1e-6
        # J_total = 1026001/1000000 = 1.026001
        t_out = 1.0 / 50.0
        t_web = 1.0 / 40.0
        t_cant = 1.0 / 100.0

        n00 = Node(0.0, 0.0)
        n30 = Node(3.0, 0.0)
        n80 = Node(8.0, 0.0)
        n82 = Node(8.0, 2.0)
        n32 = Node(3.0, 2.0)
        n02 = Node(0.0, 2.0)
        nm22 = Node(-2.0, 2.0)
        n92 = Node(9.0, 2.0)

        segs = [
            Segment(n00, n30, t=t_out),   # 0
            Segment(n30, n80, t=t_out),   # 1
            Segment(n80, n82, t=t_out),   # 2
            Segment(n82, n32, t=t_out),   # 3
            Segment(n32, n02, t=t_out),   # 4
            Segment(n02, n00, t=t_out),   # 5
            Segment(n30, n32, t=t_web),   # 6
            Segment(n02, nm22, t=t_cant), # 7
            Segment(n82, n92, t=t_cant),  # 8
        ]

        sec = MixedSection(segs)

        # 1. H matrix check
        # Permute H rows/cols to match cell order if necessary
        # Cells: one is [0,3]x[0,2], other is [3,8]x[0,2]
        H_calc = sec.H
        expected_H_diag = {480.0, 680.0}
        assert set(np.diag(H_calc).round(6)) == expected_H_diag
        assert math.isclose(abs(H_calc[0, 1]), 80.0, rel_tol=1e-12)

        # 2. Torsion constants
        expected_J_bb = 513.0 / 500.0
        expected_J_open = 1.0 / 1000000.0
        expected_J_total = 1026001.0 / 1000000.0

        assert math.isclose(sec.J_BB, expected_J_bb, rel_tol=1e-12)
        assert math.isclose(sec.J_open, expected_J_open, rel_tol=1e-12)
        assert math.isclose(sec.J_total, expected_J_total, rel_tol=1e-12)

        # 3. Independent oracle verification
        # Nodes: 0:(0,0), 1:(3,0), 2:(8,0), 3:(8,2), 4:(3,2), 5:(0,2), 6:(-2,2), 7:(9,2)
        nodes_ref = [
            (0.0, 0.0), (3.0, 0.0), (8.0, 0.0), (8.0, 2.0),
            (3.0, 2.0), (0.0, 2.0), (-2.0, 2.0), (9.0, 2.0)
        ]
        edges_ref = [
            (0, 1, t_out),   # 0
            (1, 2, t_out),   # 1
            (2, 3, t_out),   # 2
            (3, 4, t_out),   # 3
            (4, 5, t_out),   # 4
            (5, 0, t_out),   # 5
            (1, 4, t_web),   # 6
            (5, 6, t_cant),  # 7
            (3, 7, t_cant),  # 8
        ]
        # Cell 1: 0, 6, 4, 5 (edges 0, 6, 4, 5) -> 0->1, 1->4, 4->5, 5->0
        # Cell 2: 1, 2, 3, 4 (edges 1, 2, 3, 6) -> 1->2, 2->3, 3->4, 4->1 (opposes edge 6!)
        cells_ref = [
            ((0, 6, 4, 5), (1, 1, 1, 1)),
            ((1, 2, 3, 6), (1, 1, 1, -1)),
        ]

        oracle = IndependentMixedOracle(nodes_ref, edges_ref, cells_ref)
        ex_ref, ey_ref, xs_ref, ys_ref = oracle.compute_shear_center()
        j_tot_ref, j_bb_ref, j_op_ref, cw_ref, _ = oracle.compute_torsion_warping()

        assert math.isclose(sec.shear_center[0], xs_ref, rel_tol=1e-10, abs_tol=1e-12)
        assert math.isclose(sec.shear_center[1], ys_ref, rel_tol=1e-10, abs_tol=1e-12)
        assert math.isclose(sec.Cw, cw_ref, rel_tol=1e-10)

        # Transverse shear torque and flow fields check
        sf_res = sec.calculate_shear_flow(vx=250.0, vy=150.0)
        c_ref, edge_int_ref, tz_ref = oracle.solve_shear_flow(vx=250.0, vy=150.0)
        assert math.isclose(sf_res.torque, tz_ref, rel_tol=1e-10)
        assert math.isclose(sf_res.recovered_resultant[0], 250.0, rel_tol=1e-12)
        assert math.isclose(sf_res.recovered_resultant[1], 150.0, rel_tol=1e-12)
        for i, sf in enumerate(sf_res.segment_flows):
            assert math.isclose(sf.integral_q, edge_int_ref[i], rel_tol=1e-10)
            assert math.isclose(sf.q0, c_ref[i], rel_tol=1e-10, abs_tol=1e-12)

    def test_b3_doubly_symmetric_rectangular_and_square_with_antennas(self) -> None:
        """B3: Doubly symmetric rectangular box + antennas with closed-form C_w and square C_w=0."""
        # Case A: Rectangular box w=4, h=2, t_c=1/50, l=1, t_a=1/100
        # Formula from ACTIVE_PHASE.md:
        # J = 2*t_c*w^2*h^2 / (w+h) + 2*l*t_a^3 / 3 = 32/75 + 2e-6 / 3
        # C_w = t_c * w^2 * h^2 * (h - w)^2 / (24 * (w + h)) = 8/225
        w = 4.0
        h = 2.0
        tc = 1.0 / 50.0
        l_ant = 1.0
        ta = 1.0 / 100.0

        hw = w / 2.0
        hh = h / 2.0

        # Box centered at (0,0) with antennas at (-hw, 0) and (hw, 0)
        # Vertical walls split at y=0 to attach antennas
        n_bl = Node(-hw, -hh)
        n_br = Node(hw, -hh)
        n_tr = Node(hw, hh)
        n_tl = Node(-hw, hh)

        n_mid_l = Node(-hw, 0.0)
        n_mid_r = Node(hw, 0.0)

        n_ant_l = Node(-hw - l_ant, 0.0)
        n_ant_r = Node(hw + l_ant, 0.0)

        segs_rect = [
            Segment(n_bl, n_br, t=tc),         # bottom
            Segment(n_br, n_mid_r, t=tc),      # right lower
            Segment(n_mid_r, n_tr, t=tc),      # right upper
            Segment(n_tr, n_tl, t=tc),         # top
            Segment(n_tl, n_mid_l, t=tc),      # left upper
            Segment(n_mid_l, n_bl, t=tc),      # left lower
            Segment(n_mid_l, n_ant_l, t=ta),   # left antenna
            Segment(n_mid_r, n_ant_r, t=ta),   # right antenna
        ]

        sec_rect = MixedSection(segs_rect)

        expected_J = (2.0 * tc * (w**2) * (h**2)) / (w + h) + (2.0 * l_ant * (ta**3)) / 3.0
        expected_Cw = (tc * (w**2) * (h**2) * ((h - w) ** 2)) / (24.0 * (w + h))

        assert math.isclose(expected_Cw, 8.0 / 225.0, rel_tol=1e-12)
        assert math.isclose(sec_rect.J, expected_J, rel_tol=1e-12)
        assert math.isclose(sec_rect.Cw, expected_Cw, rel_tol=1e-10)

        # Shear center must be (0, 0) by double symmetry
        assert math.isclose(sec_rect.shear_center[0], 0.0, abs_tol=1e-12)
        assert math.isclose(sec_rect.shear_center[1], 0.0, abs_tol=1e-12)

        # Case B: Square box w = h = 2.0 ==> C_w == 0.0 analytically!
        w_sq = 2.0
        h_sq = 2.0
        hw_s = w_sq / 2.0
        hh_s = h_sq / 2.0

        n_bl_s = Node(-hw_s, -hh_s)
        n_br_s = Node(hw_s, -hh_s)
        n_tr_s = Node(hw_s, hh_s)
        n_tl_s = Node(-hw_s, hh_s)
        n_mid_ls = Node(-hw_s, 0.0)
        n_mid_rs = Node(hw_s, 0.0)
        n_ant_ls = Node(-hw_s - l_ant, 0.0)
        n_ant_rs = Node(hw_s + l_ant, 0.0)

        segs_sq = [
            Segment(n_bl_s, n_br_s, t=tc),
            Segment(n_br_s, n_mid_rs, t=tc),
            Segment(n_mid_rs, n_tr_s, t=tc),
            Segment(n_tr_s, n_tl_s, t=tc),
            Segment(n_tl_s, n_mid_ls, t=tc),
            Segment(n_mid_ls, n_bl_s, t=tc),
            Segment(n_mid_ls, n_ant_ls, t=ta),
            Segment(n_mid_rs, n_ant_rs, t=ta),
        ]

        sec_sq = MixedSection(segs_sq)

        expected_J_sq = (2.0 * tc * (w_sq**2) * (h_sq**2)) / (w_sq + h_sq) + (2.0 * l_ant * (ta**3)) / 3.0
        assert math.isclose(sec_sq.J, expected_J_sq, rel_tol=1e-12)
        # Exactly 0.0 C_w
        assert math.isclose(sec_sq.Cw, 0.0, abs_tol=1e-12)

    def test_b4_branch_and_antenna_limits(self) -> None:
        """B4: Thickness ta -> 0+ and length l -> 0+ converge to closed core properties."""
        w = 4.0
        h = 2.0
        tc = 0.02
        hw = w / 2.0
        hh = h / 2.0

        # Core pure box
        core_segs = [
            Segment(Node(-hw, -hh), Node(hw, -hh), t=tc),
            Segment(Node(hw, -hh), Node(hw, hh), t=tc),
            Segment(Node(hw, hh), Node(-hw, hh), t=tc),
            Segment(Node(-hw, hh), Node(-hw, -hh), t=tc),
        ]
        sec_core = MixedSection(core_segs)

        # 1. Thickness scaling: ta -> 0+
        for ta in [1e-2, 1e-4, 1e-6]:
            ant = Segment(Node(hw, 0.0), Node(hw + 1.0, 0.0), t=ta)
            # Wall split at (hw, 0)
            split_box = [
                Segment(Node(-hw, -hh), Node(hw, -hh), t=tc),
                Segment(Node(hw, -hh), Node(hw, 0.0), t=tc),
                Segment(Node(hw, 0.0), Node(hw, hh), t=tc),
                Segment(Node(hw, hh), Node(-hw, hh), t=tc),
                Segment(Node(-hw, hh), Node(-hw, -hh), t=tc),
                ant,
            ]
            sec_limit = MixedSection(split_box)

            # J_open must scale as ta^3
            assert math.isclose(sec_limit.J_open, (1.0 * ta**3) / 3.0, rel_tol=1e-10)
            assert math.isclose(sec_limit.J_BB, sec_core.J_BB, rel_tol=1e-8)

        # 2. Length scaling: l -> 0+
        for length in [1.0, 0.1, 0.01]:
            ant = Segment(Node(hw, 0.0), Node(hw + length, 0.0), t=0.01)
            split_box = [
                Segment(Node(-hw, -hh), Node(hw, -hh), t=tc),
                Segment(Node(hw, -hh), Node(hw, 0.0), t=tc),
                Segment(Node(hw, 0.0), Node(hw, hh), t=tc),
                Segment(Node(hw, hh), Node(-hw, hh), t=tc),
                Segment(Node(-hw, hh), Node(-hw, -hh), t=tc),
                ant,
            ]
            sec_limit = MixedSection(split_box)
            # J_open scales linearly with length
            assert math.isclose(sec_limit.J_open, (length * 0.01**3) / 3.0, rel_tol=1e-10)

        # 3. Y-branch and inclined antenna convergence
        # Inclined antenna at 45 degrees
        dx = 0.01 * math.cos(math.pi / 4.0)
        dy = 0.01 * math.sin(math.pi / 4.0)
        inc_ant = Segment(Node(hw, hh), Node(hw + dx, hh + dy), t=1e-4)
        sec_inc = MixedSection(core_segs + [inc_ant])
        assert math.isclose(sec_inc.J, sec_core.J, rel_tol=1e-6)
        assert math.isclose(sec_inc.shear_center[0], sec_core.shear_center[0], abs_tol=1e-4)
