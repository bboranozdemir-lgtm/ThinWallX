"""Sampling is exclusively an independent test comparison, never production."""
import math
import numpy as np
import pytest
from sectalix import SegmentStressProfile


@pytest.mark.parametrize("surface",[0.0,0.3,4.0])
def test_10000_point_grid_and_exact_interior_peak(surface):
    # Interior extremum lies exactly on one of 10000 grid nodes. This allows
    # a <1e-11 comparison without claiming a general grid has that accuracy.
    r=4000/9999
    p=SegmentStressProfile("p",1,.02,(0,.7),(-1,2*r,1-r*r),surface)
    xs=np.linspace(0,1,10000)
    grid=np.sqrt(.7**2+3*(np.abs(-xs*xs+2*r*xs+1-r*r)+surface)**2)
    assert abs(p.peak_sigma_vm-float(max(grid)))<1e-11
    assert p.s_peak==pytest.approx(r,abs=1e-13)


def test_cubic_stationary_point():
    # P=(k*x)^2+3*(1-x²)^2, interior stationary roots ±sqrt(1-k²/6).
    # Those are minima; evaluating only stationary roots would miss the ends.
    p=SegmentStressProfile(1,1,.1,(1,0),(-1,0,1),0)
    assert p.peak_sigma_vm==pytest.approx(math.sqrt(3),rel=1e-14)
    assert p.s_peak==0


def test_nontrivial_cubic_global_peak():
    # Exact rational stationary point x=1/5, sigma=x/5+143/25,
    # tau=1-x². The other stationary point is a minimum.
    p=SegmentStressProfile(1,1,.1,(.2,5.72),(-1,0,1),0)
    assert p.s_peak==pytest.approx(.2,abs=1e-14)
    assert p.peak_sigma_vm==pytest.approx(math.sqrt(35.9424),rel=1e-14)


@pytest.mark.parametrize("quadratic",[1e-20,1e-100,1e-150])
def test_nearly_linear_shear(quadratic):
    p=SegmentStressProfile(1,1,.1,(.2,.3),(quadratic,-1,.5),.1)
    assert p.s_peak==1
    assert p.peak_sigma_vm==pytest.approx(math.hypot(.5,math.sqrt(3)*.6),rel=1e-14)


@pytest.mark.parametrize("quadratic",[1e-150,1e-160,1e-300])
def test_convex_tiny_quadratic_skips_companion_solve(monkeypatch,quadratic):
    def forbidden(*args,**kwargs):
        raise AssertionError("A convex envelope branch needs no cubic roots.")
    monkeypatch.setattr(np,"roots",forbidden)
    p=SegmentStressProfile(1,1,.1,(.2,.3),(quadratic,-1,.5),.1)
    assert p.s_peak==1
    assert p.eval_tau_surface(.5)==pytest.approx(.1,abs=1e-15)
    assert p.peak_sigma_vm==pytest.approx(math.hypot(.5,math.sqrt(3)*.6),rel=1e-14)


def test_zero_constant_linear_and_sign_change():
    z=SegmentStressProfile(0,2,.1,(0,0),(0,0,0),0)
    assert z.peak_sigma_vm==0 and z.s_peak==0
    p=SegmentStressProfile(1,2,.1,(2,-1),(0,1,-1),2)
    assert p.eval_tau_surface(1)==2
    assert p.peak_sigma_vm==pytest.approx(math.sqrt(36))
    assert p.s_peak==2


def test_random_quartics_dominate_grid_with_resolution_bound():
    rng=np.random.default_rng(709)
    xs=np.linspace(0,1,10000)
    for _ in range(150):
        s=rng.uniform(-2,2,2); t=rng.uniform(-2,2,3); a=float(rng.uniform(0,2))
        p=SegmentStressProfile(0,1,.1,tuple(s),tuple(t),a)
        g=np.sqrt((s[0]*xs+s[1])**2+3*(np.abs(t[0]*xs**2+t[1]*xs+t[2])+a)**2)
        assert p.peak_sigma_vm>=float(max(g))-1e-12
        # Independent grid resolution bound, not a relaxation of root accuracy.
        assert p.peak_sigma_vm-float(max(g))<1e-6
        assert p.eval_sigma_vm(p.s_peak)==p.peak_sigma_vm
