"""Unit tests for Section class API, builders, and property accessors."""

import math
import numpy as np
import pytest

from sectalix.primitives import Node, Segment
from sectalix.section import Section


def test_section_builder_from_tuples() -> None:
    # Build channel using from_tuples
    tuples = [
        ((50.0, 50.0), (0.0, 50.0), 4.0),
        ((0.0, 50.0), (0.0, -50.0), 4.0),
        ((0.0, -50.0), (50.0, -50.0), 4.0),
    ]
    sec = Section.from_tuples(tuples)
    assert len(sec.segments) == 3
    assert len(sec.nodes) == 4
    assert math.isclose(sec.area, 800.0)
    assert math.isclose(sec.cy, 0.0, abs_tol=1e-12)


def test_section_builder_from_segments() -> None:
    s1 = Segment(Node(0.0, 0.0), Node(10.0, 0.0), t=2.0)
    sec = Section.from_segments([s1])
    assert len(sec.segments) == 1
    assert math.isclose(sec.area, 20.0)


def test_section_property_accessors() -> None:
    s1 = Segment(Node(0.0, 0.0), Node(20.0, 0.0), t=2.0)
    s2 = Segment(Node(20.0, 0.0), Node(20.0, 10.0), t=2.0)
    sec = Section([s1, s2])

    assert sec.area == sec.properties.area
    assert sec.cx == sec.properties.cx
    assert sec.cy == sec.properties.cy
    assert sec.centroid == (sec.cx, sec.cy)
    assert sec.ix_raw == sec.properties.ix_raw
    assert sec.iy_raw == sec.properties.iy_raw
    assert sec.ixy_raw == sec.properties.ixy_raw
    assert sec.Ix == sec.properties.ix
    assert sec.Iy == sec.properties.iy
    assert sec.Ixy == sec.properties.ixy
    assert sec.I1 == sec.properties.i1
    assert sec.I2 == sec.properties.i2
    assert sec.theta_p == sec.properties.theta_p
    assert sec.theta_p_deg == math.degrees(sec.theta_p)
    assert sec.is_degenerate == sec.properties.is_degenerate
    np.testing.assert_array_equal(sec.inertia_matrix, sec.properties.inertia_matrix)
