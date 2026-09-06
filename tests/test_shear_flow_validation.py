"""Unit tests for shear-flow validation and error handling (Sectalix v0.2).

Verifies acceptance criteria 16 and 17:
- Criterion 16: Singular / effectively rank-deficient inertia matrices fail explicitly
  without pseudoinverse or silent clamping (SingularSectionError).
- Criterion 17: Non-finite shear loads fail explicitly (GeometryError).
"""

import math
import numpy as np
import pytest


from sectalix.exceptions import GeometryError, SingularSectionError
from sectalix.primitives import Node, Segment
from sectalix.section import Section
from sectalix.shear_flow import calculate_shear_flow
from sectalix.shear_load import ShearLoad


class TestShearLoadValidation:
    """Criterion 17: Non-finite shear loads fail explicitly."""

    def test_finite_shear_loads(self) -> None:
        load = ShearLoad(100.0, -250.0)
        assert load.vx == 100.0
        assert load.vy == -250.0
        assert load.vector.tolist() == [100.0, -250.0]

    def test_rejects_nan_and_inf(self) -> None:
        with pytest.raises(GeometryError, match="must be finite"):
            ShearLoad(float("nan"), 10.0)
        with pytest.raises(GeometryError, match="must be finite"):
            ShearLoad(10.0, float("inf"))
        with pytest.raises(GeometryError, match="must be finite"):
            ShearLoad(float("-inf"), 0.0)

    def test_section_calculate_shear_flow_rejects_nan_and_inf(self) -> None:
        # Build valid angle section
        s1 = Segment(Node(0.0, 50.0), Node(0.0, 0.0), t=3.0)
        s2 = Segment(Node(0.0, 0.0), Node(50.0, 0.0), t=3.0)
        sec = Section([s1, s2])

        with pytest.raises(GeometryError, match="must be finite"):
            sec.calculate_shear_flow(float("nan"), 100.0)
        with pytest.raises(GeometryError, match="must be finite"):
            sec.calculate_shear_flow(100.0, float("inf"))
        with pytest.raises(GeometryError, match="vy must be provided"):
            sec.calculate_shear_flow(100.0)

    def test_non_section_input_rejected(self) -> None:
        with pytest.raises(GeometryError, match="Expected a Section"):
            calculate_shear_flow("not_a_section", 100.0, 100.0)  # type: ignore

    def test_unvalidated_invalid_section_rejected(self) -> None:
        from sectalix.exceptions import TopologyError

        # Build closed loop with validate=False
        n1 = Node(0.0, 0.0)
        n2 = Node(10.0, 0.0)
        n3 = Node(0.0, 10.0)
        loop = [Segment(n1, n2, t=1.0), Segment(n2, n3, t=1.0), Segment(n3, n1, t=1.0)]
        unvalidated_sec = Section(loop, validate=False)
        with pytest.raises(TopologyError, match="Closed-loop geometry detected"):
            calculate_shear_flow(unvalidated_sec, 0.0, 100.0)

    def test_section_validate_method(self) -> None:
        s1 = Segment(Node(0.0, 50.0), Node(0.0, 0.0), t=3.0)
        s2 = Segment(Node(0.0, 0.0), Node(50.0, 0.0), t=3.0)
        sec = Section([s1, s2], validate=False)
        assert sec._is_validated is False
        sec.validate()
        assert sec._is_validated is True



class TestSingularSectionValidation:
    """Criterion 16: Singular or rank-deficient inertia matrices fail explicitly."""

    def test_single_straight_wall_fails_shear_flow(self) -> None:
        # A single straight centerline wall has Iy = 0 (or Ix = 0 if horizontal)
        # and lambda_min = 0, so C is singular in 2D.
        s = Segment(Node(0.0, 0.0), Node(0.0, 100.0), t=4.0)
        sec = Section([s])

        with pytest.raises(SingularSectionError, match="singular or effectively rank-deficient"):
            sec.calculate_shear_flow(vx=100.0, vy=0.0)

        with pytest.raises(SingularSectionError, match="singular or effectively rank-deficient"):
            calculate_shear_flow(sec, vx=0.0, vy=100.0)

    def test_two_collinear_segments_fail_shear_flow(self) -> None:
        # Two collinear vertical segments still form a single straight line
        s1 = Segment(Node(0.0, 0.0), Node(0.0, 50.0), t=3.0)
        s2 = Segment(Node(0.0, 50.0), Node(0.0, 100.0), t=3.0)
        sec = Section([s1, s2])

        with pytest.raises(SingularSectionError, match="singular or effectively rank-deficient"):
            sec.calculate_shear_flow(vx=50.0, vy=50.0)


class TestCoordinateBoundsValidation:
    """Validation of normalized position xi and arc-length s."""

    def test_xi_out_of_bounds_rejected(self) -> None:
        s1 = Segment(Node(0.0, 50.0), Node(0.0, 0.0), t=3.0)
        s2 = Segment(Node(0.0, 0.0), Node(50.0, 0.0), t=3.0)
        sec = Section([s1, s2])
        res = sec.calculate_shear_flow(0.0, 1000.0)

        with pytest.raises(GeometryError, match=r"\[0, 1\]"):
            res[0].q_at_xi(-0.1)
        with pytest.raises(GeometryError, match=r"\[0, 1\]"):
            res[0].q_at_xi(1.1)
        with pytest.raises(GeometryError, match="must be finite"):
            res[0].q_at_xi(float("nan"))

    def test_s_finite_validation(self) -> None:
        s1 = Segment(Node(0.0, 50.0), Node(0.0, 0.0), t=3.0)
        s2 = Segment(Node(0.0, 0.0), Node(50.0, 0.0), t=3.0)
        sec = Section([s1, s2])
        res = sec.calculate_shear_flow(0.0, 1000.0)

        with pytest.raises(GeometryError, match="must be finite"):
            res[0].q_at(float("inf"))


class TestScaleIndependenceAndOverflow:
    """Scale independence of singularity check and overflow protection."""

    def test_meter_scale_small_inertia_section_succeeds(self) -> None:
        # Thin-walled channel defined in meters:
        # Dimensions: web height 0.001 m (1 mm), flanges 0.0005 m (0.5 mm), thickness 0.0001 m (0.1 mm)
        # Section properties have Ix, Iy ~ 10^-14 to 10^-15 m^4
        # Scale-independent singularity check must NOT raise SingularSectionError!
        h = 0.001
        b = 0.0005
        t = 0.0001
        s0 = Segment(Node(b, h / 2.0), Node(0.0, h / 2.0), t=t)
        s1 = Segment(Node(0.0, h / 2.0), Node(0.0, -h / 2.0), t=t)
        s2 = Segment(Node(0.0, -h / 2.0), Node(b, -h / 2.0), t=t)
        sec = Section([s0, s1, s2], node_tolerance=1e-12)

        assert sec.I1 < 1e-12
        assert sec.I2 < 1e-12
        # Must compute successfully without SingularSectionError:
        res = sec.calculate_shear_flow(vx=0.0, vy=1.0)
        np.testing.assert_allclose(res.recovered_resultant, [0.0, 1.0], rtol=1e-10)

    def test_extreme_load_overflow_raises_error(self) -> None:
        # Finite extreme load of 1e308 causes floating-point overflow during solve
        s1 = Segment(Node(0.0, 50.0), Node(0.0, 0.0), t=3.0)
        s2 = Segment(Node(0.0, 0.0), Node(50.0, 0.0), t=3.0)
        sec = Section([s1, s2])

        with pytest.raises(GeometryError, match="Non-finite"):
            sec.calculate_shear_flow(vx=1e308, vy=1e308)

    @pytest.mark.parametrize("bad_factor", [0.0, -10.0, float("nan"), float("inf")])
    def test_invalid_safety_factor_rejected(self, bad_factor: float) -> None:
        s1 = Segment(Node(0.0, 50.0), Node(0.0, 0.0), t=3.0)
        s2 = Segment(Node(0.0, 0.0), Node(50.0, 0.0), t=3.0)
        sec = Section([s1, s2])

        with pytest.raises(GeometryError, match="safety_factor"):
            sec.calculate_shear_flow(vx=0.0, vy=100.0, safety_factor=bad_factor)

