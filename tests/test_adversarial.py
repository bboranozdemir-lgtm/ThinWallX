"""Adversarial regression tests found during the independent v0.1 audit."""

import itertools
import math

import numpy as np
import pytest

from thinwallx.exceptions import GeometryError, TopologyError
from thinwallx.primitives import Node, Segment
from thinwallx.section import Section
from thinwallx.validation import cluster_nodes


def _asymmetric_section() -> Section:
    return Section(
        [
            Segment(Node(0.0, 120.0), Node(60.0, 120.0), 3.0),
            Segment(Node(0.0, 0.0), Node(0.0, 120.0), 4.0),
            Segment(Node(-40.0, 0.0), Node(0.0, 0.0), 5.0),
            Segment(Node(-40.0, 0.0), Node(-40.0, 25.0), 3.0),
        ]
    )


def test_large_translation_preserves_centroidal_inertia() -> None:
    """Raw-moment cancellation must not corrupt centroidal properties."""
    base = _asymmetric_section()
    moved = base.translated(1.0e12, -1.0e12)

    np.testing.assert_allclose(
        moved.inertia_matrix,
        base.inertia_matrix,
        rtol=2.0e-13,
        atol=0.0,
    )
    assert math.isclose(moved.I1, base.I1, rel_tol=2.0e-13)
    assert math.isclose(moved.I2, base.I2, rel_tol=2.0e-13)
    assert moved.I2 >= 0.0


def test_small_scale_inclined_wall_is_not_zeroed_or_marked_isotropic() -> None:
    """Moment and degeneracy tolerances must scale with the actual tensor."""
    length = 1.0e-4
    thickness = 1.0e-4
    alpha = math.pi / 6.0
    c = math.cos(alpha)
    s = math.sin(alpha)
    section = Section(
        [
            Segment(
                Node(-0.5 * length * c, -0.5 * length * s),
                Node(0.5 * length * c, 0.5 * length * s),
                thickness,
            )
        ]
    )

    line_inertia = thickness * length**3 / 12.0
    assert math.isclose(section.Ix, line_inertia * s**2, rel_tol=2.0e-15)
    assert math.isclose(section.Iy, line_inertia * c**2, rel_tol=2.0e-15)
    assert math.isclose(section.Ixy, line_inertia * s * c, rel_tol=2.0e-15)
    assert math.isclose(section.I1, line_inertia, rel_tol=2.0e-15)
    assert section.I2 == 0.0
    assert section.is_degenerate is False
    assert math.isclose(section.theta_p, alpha - math.pi / 2.0, abs_tol=2.0e-15)


def test_small_positive_minor_principal_moment_is_not_clamped_to_zero() -> None:
    """A highly anisotropic section is not the same as a rank-one section."""
    tiny_thickness = 1.0e-15
    section = Section(
        [
            Segment(Node(-1.0, 0.0), Node(0.0, 0.0), 1.0),
            Segment(Node(0.0, 0.0), Node(1.0, 0.0), 1.0),
            Segment(Node(0.0, -1.0), Node(0.0, 0.0), tiny_thickness),
            Segment(Node(0.0, 0.0), Node(0.0, 1.0), tiny_thickness),
        ]
    )

    expected_i2 = 2.0 * tiny_thickness / 3.0
    assert section.I2 > 0.0
    assert math.isclose(section.I2, expected_i2, rel_tol=2.0e-3)
    assert section.is_degenerate is False


def test_tiny_nonzero_segment_can_be_used_with_zero_clustering_tolerance() -> None:
    """Exact zero length, rather than an absolute dimensional cutoff, is invalid."""
    length = 1.0e-13
    section = Section(
        [Segment(Node(0.0, 0.0), Node(length, 0.0), 1.0)],
        node_tolerance=0.0,
    )

    assert section.area == length
    assert section.Iy > 0.0


@pytest.mark.parametrize("bad_tolerance", [-1.0, math.nan, math.inf])
def test_invalid_node_tolerance_is_rejected(bad_tolerance: float) -> None:
    segment = Segment(Node(0.0, 0.0), Node(1.0, 0.0), 1.0)
    with pytest.raises(GeometryError, match="node tolerance"):
        Section([segment], node_tolerance=bad_tolerance)


def test_proximity_clustering_is_transitive_and_order_independent() -> None:
    """A~B and B~C must define one stable topological joint."""
    nodes = [Node(0.0, 0.0), Node(0.75, 0.0), Node(1.5, 0.0)]

    for permutation in itertools.permutations(nodes):
        canonical_nodes, mapping = cluster_nodes(permutation, tol=1.0)
        assert len(canonical_nodes) == 1
        assert len(set(mapping.values())) == 1


@pytest.mark.parametrize("scale", [1.0e-6, 1.0, 1.0e6])
def test_absolute_node_tolerance_scales_with_coordinate_units(scale: float) -> None:
    """Scaling coordinates and the absolute tolerance preserves connectivity."""
    tolerance = 1.0e-6 * scale
    gap = 0.5 * tolerance
    section = Section(
        [
            Segment(Node(0.0, 0.0), Node(1.0 * scale, 0.0), scale),
            Segment(
                Node(1.0 * scale + gap, 0.0),
                Node(2.0 * scale, 0.0),
                scale,
            ),
        ],
        node_tolerance=tolerance,
    )
    assert len(section.nodes) == 3


def test_duplicate_segments_with_nearly_coincident_endpoints_are_rejected() -> None:
    tolerance = 1.0e-8
    segments = [
        Segment(Node(0.0, 0.0), Node(10.0, 0.0), 1.0),
        Segment(Node(10.0 + 0.25 * tolerance, 0.0), Node(0.25 * tolerance, 0.0), 2.0),
    ]

    with pytest.raises(TopologyError, match="Duplicate segment"):
        Section(segments, node_tolerance=tolerance)


def test_self_crossing_chain_is_rejected_as_a_geometric_loop() -> None:
    """A graph-theoretic tree can still contain a loop at an interior crossing."""
    segments = [
        Segment(Node(0.0, 0.0), Node(10.0, 10.0), 1.0),
        Segment(Node(10.0, 10.0), Node(0.0, 10.0), 1.0),
        Segment(Node(0.0, 10.0), Node(10.0, 0.0), 1.0),
    ]

    with pytest.raises(TopologyError, match="intersect"):
        Section(segments)


def test_unsplit_interior_t_junction_is_rejected_explicitly() -> None:
    """Branched joints are valid, but every junction must be an endpoint node."""
    segments = [
        Segment(Node(0.0, 0.0), Node(10.0, 0.0), 1.0),
        Segment(Node(5.0, 0.0), Node(5.0, 5.0), 1.0),
    ]

    with pytest.raises(TopologyError, match="Split all walls at intersections"):
        Section(segments)


def test_partially_overlapping_segments_are_rejected() -> None:
    """Collinear wall overlap must not be silently counted twice."""
    segments = [
        Segment(Node(0.0, 0.0), Node(10.0, 0.0), 1.0),
        Segment(Node(0.0, 0.0), Node(5.0, 0.0), 1.0),
    ]

    with pytest.raises(TopologyError, match="overlap"):
        Section(segments)


def test_generic_section_against_independent_exact_gauss_reference() -> None:
    """Two Gauss points exactly integrate the linear/quadratic line fields."""
    segments = [
        Segment(Node(-3.0, 2.0), Node(1.0, 5.0), 0.7),
        Segment(Node(1.0, 5.0), Node(4.0, -1.0), 1.3),
        Segment(Node(4.0, -1.0), Node(8.0, 3.0), 0.4),
    ]
    section = Section(segments)

    # This oracle does not call any production integral.  The two-point
    # Gauss-Legendre rule is exact for polynomials through degree three, while
    # x, y, x^2, y^2, and xy restricted to a straight segment have degree <= 2.
    gauss_points = (-1.0 / math.sqrt(3.0), 1.0 / math.sqrt(3.0))
    weighted_points: list[tuple[float, float, float]] = []
    for segment in segments:
        length = math.hypot(
            segment.p2.x - segment.p1.x,
            segment.p2.y - segment.p1.y,
        )
        point_weight = segment.t * length / 2.0
        for coordinate in gauss_points:
            x = (
                (1.0 - coordinate) * segment.p1.x
                + (1.0 + coordinate) * segment.p2.x
            ) / 2.0
            y = (
                (1.0 - coordinate) * segment.p1.y
                + (1.0 + coordinate) * segment.p2.y
            ) / 2.0
            weighted_points.append((point_weight, x, y))

    area = math.fsum(weight for weight, _, _ in weighted_points)
    cx = math.fsum(weight * x for weight, x, _ in weighted_points) / area
    cy = math.fsum(weight * y for weight, _, y in weighted_points) / area
    ix_raw = math.fsum(weight * y * y for weight, _, y in weighted_points)
    iy_raw = math.fsum(weight * x * x for weight, x, _ in weighted_points)
    ixy_raw = math.fsum(weight * x * y for weight, x, y in weighted_points)
    ix = math.fsum(weight * (y - cy) ** 2 for weight, _, y in weighted_points)
    iy = math.fsum(weight * (x - cx) ** 2 for weight, x, _ in weighted_points)
    ixy = math.fsum(
        weight * (x - cx) * (y - cy) for weight, x, y in weighted_points
    )

    assert math.isclose(section.area, area, rel_tol=2.0e-15)
    assert math.isclose(section.cx, cx, rel_tol=2.0e-15)
    assert math.isclose(section.cy, cy, rel_tol=2.0e-15)
    assert math.isclose(section.ix_raw, ix_raw, rel_tol=2.0e-15)
    assert math.isclose(section.iy_raw, iy_raw, rel_tol=2.0e-15)
    assert math.isclose(section.ixy_raw, ixy_raw, rel_tol=2.0e-15)
    assert math.isclose(section.Ix, ix, rel_tol=2.0e-15)
    assert math.isclose(section.Iy, iy, rel_tol=2.0e-15)
    assert math.isclose(section.Ixy, ixy, rel_tol=2.0e-15)

    reference_tensor = np.array([[ix, -ixy], [-ixy, iy]])
    reference_moments = np.linalg.eigvalsh(reference_tensor)
    assert math.isclose(section.I1, reference_moments[1], rel_tol=2.0e-15)
    assert math.isclose(section.I2, reference_moments[0], rel_tol=2.0e-15)

    principal_direction = np.array(
        [math.cos(section.theta_p), math.sin(section.theta_p)]
    )
    np.testing.assert_allclose(
        section.inertia_matrix @ principal_direction,
        section.I1 * principal_direction,
        rtol=2.0e-15,
        atol=2.0e-15,
    )
