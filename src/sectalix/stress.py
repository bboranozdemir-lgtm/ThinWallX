"""Linear-elastic stress recovery for supported thin-wall section models.

Polynomials are assembled in xi=s/L and only converted to dimensional coefficients
for the public API. Normal resultants follow the specified section-face convention:
N=int(sigma dA), Mx=-int(Y sigma dA), My=int(X sigma dA), B=int(omega sigma dA).
Secondary flow uses dq/ds=-M_omega*t*omega/Cw, Kirchhoff balance and zero cell twist.
No sampling, quadrature, or polynomial fitting is used by this module.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any, TYPE_CHECKING

import numpy as np

from sectalix.exceptions import GeometryError, SingularSectionError

if TYPE_CHECKING:
    from sectalix.section import Section

_OVERFLOW = "Recovered stress overflows IEEE-754 float64 range."
_UNDERFLOW = "Calculated value is non-zero but underflows IEEE-754 float64 subnormal range (< 5e-324)."
_ZERO = (0.0, 0)
Scaled = tuple[float, int]


def _scaled(value: float) -> Scaled:
    return math.frexp(float(value))


def _mul(*values: Scaled) -> Scaled:
    m, e = 1.0, 0
    for mantissa, exponent in values:
        if mantissa == 0.0:
            return _ZERO
        m, shift = math.frexp(m * mantissa)
        e += exponent + shift
    return m, e


def _div(value: Scaled, denominator: Scaled) -> Scaled:
    if denominator[0] == 0.0:
        raise SingularSectionError("Stress recovery denominator is zero.")
    m, shift = math.frexp(value[0] / denominator[0])
    return m, value[1] - denominator[1] + shift


def _sum(*values: Scaled) -> Scaled:
    terms = [v for v in values if v[0] != 0.0]
    if not terms:
        return _ZERO
    exponent = max(v[1] for v in terms)
    m = math.fsum(math.ldexp(v[0], v[1] - exponent) for v in terms)
    m, shift = math.frexp(m)
    return m, exponent + shift


def _number(value: Scaled) -> float:
    if value[0] == 0.0:
        return 0.0
    try:
        result = math.ldexp(*value)
    except OverflowError as exc:
        raise OverflowError(_OVERFLOW) from exc
    if not math.isfinite(result):
        raise OverflowError(_OVERFLOW)
    if result == 0.0:
        raise FloatingPointError(_UNDERFLOW)
    return result


def _product(*values: float) -> Scaled:
    return _mul(*(_scaled(v) for v in values))


def _neg(value: Scaled) -> Scaled:
    return -value[0], value[1]


def _unresolved_warping_resistance(section: Section, cw: float) -> bool:
    """Stress inversion requires Cw/(t_max*L_max**5) > 1e-12.

    This is a numerical resolution gate, not a change to the geometric Cw.
    The audit's relative threshold rejects roundoff-only warping in concurrent
    straight-leg sections. Compare normalized mantissa/exponent pairs so neither
    L**5 nor the threshold needs to be representable as a physical float.
    A positive subnormal is not automatically singular: an absolute float floor
    would change this decision on changing units and discard resolved geometry.
    """
    if cw == 0.0:
        return True
    length = max(s.length for s in section.segments)
    thickness = max(s.t for s in section.segments)
    threshold = _product(1e-12, thickness, length, length, length, length, length)
    mantissa, exponent = _scaled(cw)
    return (exponent, mantissa) <= (threshold[1], threshold[0])


def _poly_value(coefficients: tuple[float, ...], x: float) -> float:
    # Ascending powers, with exponent arithmetic through the complete sum.
    terms = []
    power = _scaled(1.0)
    for coefficient in coefficients:
        terms.append(_mul(_scaled(coefficient), power))
        power = _mul(power, _scaled(x))
    return _number(_sum(*terms))


def _real_roots(coefficients: list[float]) -> list[float]:
    """Closed-form real roots of degree <=3; coefficients in ascending order.

    Quadratics use the cancellation-resistant q formula; cubics use the
    algebraic companion-matrix eigenvalue formulation. Unlike Cardano's raw
    discriminant this also handles nearly quadratic coefficients without
    cubing enormous coefficient ratios. No spatial search or fitting is used.
    """
    while coefficients and coefficients[-1] == 0.0:
        coefficients = coefficients[:-1]
    degree = len(coefficients) - 1
    if degree <= 0:
        return []
    scale = max(abs(c) for c in coefficients)
    c = [v / scale for v in coefficients]
    if degree == 1:
        return [-c[0] / c[1]]
    if degree == 2:
        a, b, d = c[2], c[1], c[0]
        discriminant = math.fsum([b*b, -4*a*d])
        if discriminant < 0.0:
            return []
        q = -0.5 * (b + math.copysign(math.sqrt(discriminant), b))
        return [-b / (2*a)] if q == 0.0 else [q/a, d/q]
    roots = np.roots(c[::-1])
    return [float(r.real) for r in roots if abs(r.imag) < 1e-12]


@dataclass(frozen=True)
class AppliedLoads:
    """Consistent units: forces F, moments F L, B F L², yield stress F/L²."""

    N: float = 0.0
    Vx: float = 0.0
    Vy: float = 0.0
    Mx: float = 0.0
    My: float = 0.0
    Tsv: float = 0.0
    B: float = 0.0
    M_omega: float = 0.0
    sigma_yield: float | None = None

    def __post_init__(self) -> None:
        for name in ("N", "Vx", "Vy", "Mx", "My", "Tsv", "B", "M_omega"):
            if not math.isfinite(getattr(self, name)):
                raise GeometryError(f"Applied load {name} must be finite.")
        if self.sigma_yield is not None and (
            not math.isfinite(self.sigma_yield) or self.sigma_yield <= 0.0
        ):
            raise GeometryError("sigma_yield must be finite and strictly positive.")


@dataclass(frozen=True)
class SegmentStressProfile:
    """Dimensional descending coefficients; s follows the original segment direction.

    peak_sigma_vm and s_peak are computed on construction. On an open strip,
    tau_surface is the envelope of the two opposing surface stresses.
    """

    segment_id: Any
    length: float
    thickness: float
    sigma_zz_coeffs: tuple[float, float]
    tau_membrane_coeffs: tuple[float, float, float]
    tau_sv_surface: float
    peak_sigma_vm: float = field(init=False)
    s_peak: float = field(init=False)
    _sigma_xi: tuple[float, float] = field(init=False, repr=False)
    _tau_xi: tuple[float, float, float] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if not all(math.isfinite(v) for v in (
            self.length, self.thickness, self.tau_sv_surface,
            *self.sigma_zz_coeffs, *self.tau_membrane_coeffs,
        )):
            raise GeometryError("Stress profile inputs must be finite.")
        if self.length <= 0.0 or self.thickness <= 0.0 or self.tau_sv_surface < 0.0:
            raise GeometryError("Stress profile length/thickness must be positive and surface amplitude non-negative.")
        a, b = self.sigma_zz_coeffs
        u, v, w = self.tau_membrane_coeffs
        sigma = (b, _number(_product(a, self.length)))
        tau = (w, _number(_product(v, self.length)), _number(_product(u, self.length, self.length)))
        object.__setattr__(self, "_sigma_xi", sigma)
        object.__setattr__(self, "_tau_xi", tau)
        scale = max(abs(z) for z in (*sigma, *tau, self.tau_sv_surface))
        candidates = [0.0, 1.0]
        if scale != 0.0:
            s0, s1 = (z/scale for z in sigma)
            t0, t1, t2 = (z/scale for z in tau)
            surface = self.tau_sv_surface/scale
            candidates += _real_roots([t0, t1, t2])
            # Both smooth envelope branches are valid polynomials on their
            # corresponding sign intervals. Extraneous candidates are harmless:
            # the true envelope is always evaluated at every candidate.
            for sign in (-1.0, 1.0):
                a0, a1, a2 = sign*t0+surface, sign*t1, sign*t2
                derivative = [
                    2*s0*s1+6*a0*a1,
                    2*s1*s1+6*(a1*a1+2*a0*a2),
                    18*a1*a2, 12*a2*a2,
                ]
                # A convex branch has no interior maximum. Certifying this
                # avoids ill-conditioned, irrelevant roots far outside [0,1].
                vertex = min(1.0, max(0.0, -derivative[2]/(3*derivative[3]))) if derivative[3] else 0.0
                curvature = derivative[1]+2*derivative[2]*vertex+3*derivative[3]*vertex*vertex
                if curvature < 0.0:
                    candidates += _real_roots(derivative)
        best, best_x = -1.0, 0.0
        for x in sorted(set(candidates)):
            if math.isfinite(x) and 0.0 <= x <= 1.0:
                value = self._vm_at_xi(x)
                if value > best:
                    best, best_x = value, x
        object.__setattr__(self, "peak_sigma_vm", best)
        object.__setattr__(self, "s_peak", best_x*self.length)

    def _xi(self, s: float) -> float:
        if not math.isfinite(s) or s < 0.0 or s > self.length:
            raise GeometryError("Stress evaluation coordinate s must lie in [0, L].")
        return s/self.length

    def eval_sigma_zz(self, s: float) -> float:
        return _poly_value(self._sigma_xi, self._xi(s))

    def eval_tau_membrane(self, s: float) -> float:
        return _poly_value(self._tau_xi, self._xi(s))

    def eval_tau_surface(self, s: float) -> float:
        return _number(_sum(_scaled(abs(self.eval_tau_membrane(s))), _scaled(self.tau_sv_surface)))

    def _vm_at_xi(self, x: float) -> float:
        sigma = _poly_value(self._sigma_xi, x)
        tau = _number(_sum(_scaled(abs(_poly_value(self._tau_xi, x))), _scaled(self.tau_sv_surface)))
        scale = max(abs(sigma), tau)
        if scale == 0.0:
            return 0.0
        return _number(_product(scale, math.hypot(sigma/scale, math.sqrt(3.0)*(tau/scale))))

    def eval_sigma_vm(self, s: float) -> float:
        return self._vm_at_xi(self._xi(s))


@dataclass(frozen=True)
class StressRecoveryResult:
    segment_profiles: dict[Any, SegmentStressProfile]
    max_sigma_vm: float
    peak_segment_id: Any
    peak_s: float
    peak_location_xy: tuple[float, float]
    load_factor: float | None
    resultant_N: float
    resultant_Mx: float
    resultant_My: float
    resultant_B: float


def _secondary_static_moments(section: Section, topology: Any, warp: Any) -> list[list[Scaled]]:
    """Integrate t*omega, balance a spanning tree, then apply H circulations.

    S=c+a1*xi+a2*xi². At every node D*c=-sum(incoming(a1+a2)).
    Chord constants initially vanish; closed-block corrections enforce int S/t=0.
    This is an endpoint-cut particular solution; final fields are cut independent.
    """
    edges = topology.canonical_edges
    count = topology.node_count
    coefficients: list[list[Scaled]] = []
    rhs = [_ZERO for _ in range(count)]
    for (u, v), seg, sw in zip(edges, section.segments, warp.segment_warpings):
        a1 = _product(seg.t, seg.length, sw.omega1)
        a2 = _mul(_product(seg.t, seg.length, 0.5), _sum(_scaled(sw.omega2), _scaled(-sw.omega1)))
        coefficients.append([_ZERO, a1, a2])
        rhs[v] = _sum(rhs[v], _neg(a1), _neg(a2))
    adjacency: list[list[tuple[int, int]]] = [[] for _ in range(count)]
    cuts = set(topology.cut_edges)
    for i, (u, v) in enumerate(edges):
        if i not in cuts:
            adjacency[u].append((v, i))
            adjacency[v].append((u, i))
    parents = {0: (-1, -1)}
    order = [0]
    for node in order:
        for neighbor, edge in adjacency[node]:
            if neighbor not in parents:
                parents[neighbor] = (node, edge)
                order.append(neighbor)
    for node in reversed(order[1:]):
        parent, edge = parents[node]
        coefficients[edge][0] = rhs[node] if edges[edge][1] == node else _neg(rhs[node])
        rhs[parent] = _sum(rhs[parent], rhs[node])

    for ids, H, eH in zip(topology.block_cell_indices, topology.H_block_scaled, topology.e_max_H_block):
        B = topology.B[np.array(ids), :]
        integrals = []
        for seg, (c0, c1, c2) in zip(section.segments, coefficients):
            mean = _sum(c0, _mul(c1, _scaled(0.5)), _div(c2, _scaled(3.0)))
            integrals.append(_div(_mul(mean, _scaled(seg.length)), _scaled(seg.t)))
        b = [_sum(*(_mul(_scaled(float(sign)), value) for sign, value in zip(row, integrals))) for row in B]
        nonzero = [e for m, e in b if m != 0.0]
        eb = max(nonzero) if nonzero else 0
        scaled_b = np.array([math.ldexp(m, e-eb) if m else 0.0 for m, e in b])
        solution = np.linalg.solve(H, -scaled_b)
        for edge, value in enumerate(B.T @ solution):
            m, e = math.frexp(float(value))
            coefficients[edge][0] = _sum(coefficients[edge][0], (m, e+eb-eH))
    return coefficients


def calculate_stresses(section: Section, loads: AppliedLoads) -> StressRecoveryResult:
    """Recover section-face stresses and exact polynomial peak envelopes.

    Missing segment IDs use their input indices. Duplicate/colliding IDs are
    rejected so the profile dictionary cannot silently discard a physical wall.
    Zero load with a yield stress has infinite load factor (no finite yielding load).
    """
    from sectalix.section import Section
    from sectalix.mixed_topology import extract_mixed_topology
    from sectalix.mixed_torsion import _compute_open_torsion_terms, _solve_bredt_batho_blocks
    if not isinstance(section, Section) or not isinstance(loads, AppliedLoads):
        raise TypeError("calculate_stresses expects a Section and AppliedLoads.")
    section.validate()
    identifiers = [seg.id if seg.id is not None else i for i, seg in enumerate(section.segments)]
    try:
        if len(set(identifiers)) != len(identifiers):
            raise GeometryError("Stress profiles require unique segment IDs (including fallback indices).")
    except TypeError as exc:
        raise GeometryError("Stress profile segment IDs must be hashable.") from exc

    topology = extract_mixed_topology(section.segments, node_tolerance=section._node_tolerance)
    # Use each class's accepted warping implementation; no warping calculation
    # is needed for ordinary axial/bending/transverse recovery.
    warp = section.torsion_properties() if (loads.B or loads.M_omega) else None
    cw = warp.Cw if warp is not None else 0.0
    unresolved_cw = _unresolved_warping_resistance(section, cw) if warp is not None else False
    if unresolved_cw and loads.B != 0.0:
        raise SingularSectionError("Section has zero warping resistance (Cw=0). Bimoment B cannot be sustained.")
    if unresolved_cw and loads.M_omega != 0.0:
        raise SingularSectionError("Section has zero warping resistance (Cw=0). Warping moment M_omega cannot be sustained.")
    torsion_constant = 0.0
    membrane_torsion: list[Scaled] = [_ZERO]*len(section.segments)
    if loads.Tsv:
        j_open, _ = _compute_open_torsion_terms(section, topology)
        j_bb, _, _, f_scaled, f_exponents, _, _ = _solve_bredt_batho_blocks(
            topology, require_physical_fields=False
        )
        torsion_constant = _number(_sum(_scaled(j_open), _scaled(j_bb)))
        membrane_torsion = [(float(m), int(e)) for m,e in zip(f_scaled,f_exponents)]

    refx, refy = section.segments[0].p1.coords
    cx = _div(_sum(*(_product(s.t, s.length, 0.5, (s.p1.x-refx)+(s.p2.x-refx)) for s in section.segments)), _scaled(section.area))
    cy = _div(_sum(*(_product(s.t, s.length, 0.5, (s.p1.y-refy)+(s.p2.y-refy)) for s in section.segments)), _scaled(section.area))
    xcenter, ycenter = _number(cx), _number(cy)
    bx = by = _ZERO
    if loads.Mx or loads.My:
        inertia_scale = max(section.Ix, section.Iy, abs(section.Ixy))
        if inertia_scale == 0.0:
            raise SingularSectionError("Bending inertia is singular.")
        mat = np.array([[section.Iy/inertia_scale, section.Ixy/inertia_scale], [section.Ixy/inertia_scale, section.Ix/inertia_scale]])
        ev = np.linalg.eigvalsh(mat)
        if ev[0] <= 1e4*np.finfo(float).eps*ev[-1]:
            raise SingularSectionError("Bending inertia is singular or ill-conditioned.")
        determinant = _sum(_product(section.Ix, section.Iy), _neg(_product(section.Ixy, section.Ixy)))
        bx = _div(_sum(_product(loads.My, section.Ix), _product(loads.Mx, section.Ixy)), determinant)
        by = _neg(_div(_sum(_product(loads.Mx, section.Iy), _product(loads.My, section.Ixy)), determinant))
    flow = None
    load_exponent = 0
    if loads.Vx or loads.Vy:
        load_exponent = math.frexp(max(abs(loads.Vx), abs(loads.Vy)))[1]
        try:
            with np.errstate(over="raise", invalid="raise", divide="raise"):
                flow = section.calculate_shear_flow(math.ldexp(loads.Vx, -load_exponent), math.ldexp(loads.Vy, -load_exponent))
        except (OverflowError, FloatingPointError) as exc:
            raise OverflowError(_OVERFLOW) from exc
    secondary = _secondary_static_moments(section, topology, warp) if loads.M_omega else None
    profiles = {}
    normal_integrals: list[list[Scaled]] = [[], [], [], []]
    for i, (key, seg) in enumerate(zip(identifiers, section.segments)):
        L, t = seg.length, seg.t
        X, Y = (seg.p1.x-refx)-xcenter, (seg.p1.y-refy)-ycenter
        dx, dy = seg.p2.x-seg.p1.x, seg.p2.y-seg.p1.y
        w0 = warp.segment_warpings[i].omega1 if warp else 0.0
        dw = (warp.segment_warpings[i].omega2-w0) if warp else 0.0
        b0 = _div(_product(loads.B, w0), _scaled(cw)) if loads.B else _ZERO
        b1 = _div(_product(loads.B, dw), _scaled(cw)) if loads.B else _ZERO
        sigma = [
            _sum(_div(_scaled(loads.N), _scaled(section.area)), _mul(bx, _scaled(X)), _mul(by, _scaled(Y)), b0),
            _sum(_mul(bx, _scaled(dx)), _mul(by, _scaled(dy)), b1),
        ]
        tau = [_ZERO, _ZERO, _ZERO]
        if flow is not None:
            tau[0] = _div(_mul(_scaled(flow.segment_flows[i].q0), (0.5, load_exponent+1)), _scaled(t))
            ax, ay = map(float, flow.alpha)
            tau[1] = _mul(_sum(_product(-L, ax, X), _product(-L, ay, Y)), (0.5, load_exponent+1))
            tau[2] = _mul(_sum(_product(-0.5, L, ax, dx), _product(-0.5, L, ay, dy)), (0.5, load_exponent+1))
        if secondary is not None:
            factor = _neg(_div(_scaled(loads.M_omega), _product(cw, t)))
            tau = [_sum(v, _mul(factor, moment)) for v, moment in zip(tau, secondary[i])]
        surface = 0.0
        if loads.Tsv:
            if i in topology.closed_edges:
                tau[0] = _sum(tau[0], _div(_mul(_scaled(loads.Tsv), membrane_torsion[i]), _product(torsion_constant, t)))
            else:
                surface = _number(_div(_product(abs(loads.Tsv), t), _scaled(torsion_constant)))
        profile = SegmentStressProfile(
            key, L, t,
            (_number(_div(sigma[1], _scaled(L))), _number(sigma[0])),
            (_number(_div(tau[2], _product(L, L))), _number(_div(tau[1], _scaled(L))), _number(tau[0])),
            surface,
        )
        profiles[key] = profile
        # int_0^1 (s0+s1 xi)*(v0+v1 xi) dxi, evaluated without
        # materializing possibly overflowing intermediate products.
        for output, (v0, v1) in zip(normal_integrals, [(1.0, 0.0), (-Y, -dy), (X, dx), (w0, dw)]):
            mean = _sum(_mul(sigma[0], _scaled(v0)), _mul(sigma[0], _product(v1, 0.5)), _mul(sigma[1], _product(v0, 0.5)), _div(_mul(sigma[1], _scaled(v1)), _scaled(3.0)))
            output.append(_mul(_product(t, L), mean))
    peak = max(profiles.values(), key=lambda p: p.peak_sigma_vm)
    index = identifiers.index(peak.segment_id)
    seg = section.segments[index]
    xi = peak.s_peak/seg.length
    location = (seg.p1.x+xi*(seg.p2.x-seg.p1.x), seg.p1.y+xi*(seg.p2.y-seg.p1.y))
    factor = None if loads.sigma_yield is None else (
        _number(_div(_scaled(loads.sigma_yield), _scaled(peak.peak_sigma_vm)))
        if peak.peak_sigma_vm else math.inf
    )
    return StressRecoveryResult(profiles, peak.peak_sigma_vm, peak.segment_id, peak.s_peak, location, factor, *(_number(_sum(*v)) for v in normal_integrals))
