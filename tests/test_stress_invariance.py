"""v0.7 covariance, validation and IEEE-754 stress recovery tests."""
from dataclasses import replace
import math
import sys

import pytest

from thinwallx import (AppliedLoads, ClosedSection, GeometryError, MixedSection,
                      Node, Section, Segment, SegmentStressProfile,
                      SingularSectionError, calculate_stresses)
from tests.test_stress_benchmarks import fixtures


LOADS=AppliedLoads(N=3,Vx=2,Vy=-5,Mx=7,My=-11,Tsv=13,B=17,M_omega=19)


def transformed(section, angle=0.0, shift=(0.0,0.0), reverse=False):
    c,s=math.cos(angle),math.sin(angle)
    def point(n):
        return Node(c*n.x-s*n.y+shift[0],s*n.x+c*n.y+shift[1])
    edges=[Segment(point(e.p1),point(e.p2),e.t,e.id) for e in section.segments]
    if reverse:
        edges=[e.reversed() for e in reversed(edges)]
    return type(section)(edges)


@pytest.mark.parametrize("kind",["I","box","hat"])
@pytest.mark.parametrize("shift",[(1e12,0),(0,1e12),(1e12,-1e12)])
def test_large_translation(kind,shift):
    section,*_=fixtures(kind)
    base=section.stresses(LOADS)
    moved=transformed(section,shift=shift).stresses(LOADS)
    assert moved.max_sigma_vm==pytest.approx(base.max_sigma_vm,rel=1e-10,abs=0)
    for key,p in base.segment_profiles.items():
        q=moved.segment_profiles[key]
        for xi in (0,.23,.5,1):
            assert q.eval_sigma_zz(xi*q.length)==pytest.approx(p.eval_sigma_zz(xi*p.length),rel=1e-10,abs=1e-11)
            assert q.eval_tau_membrane(xi*q.length)==pytest.approx(p.eval_tau_membrane(xi*p.length),rel=1e-10,abs=1e-11)
            assert q.eval_sigma_vm(xi*q.length)==pytest.approx(p.eval_sigma_vm(xi*p.length),rel=1e-10,abs=0)


@pytest.mark.parametrize("kind",["I","box","hat"])
@pytest.mark.parametrize("angle",[.37,1.13])
def test_rotation_covariance(kind,angle):
    sec,*_=fixtures(kind)
    c,s=math.cos(angle),math.sin(angle)
    load=replace(LOADS,Vx=c*LOADS.Vx-s*LOADS.Vy,Vy=s*LOADS.Vx+c*LOADS.Vy,
                 Mx=c*LOADS.Mx-s*LOADS.My,My=s*LOADS.Mx+c*LOADS.My)
    a,b=sec.stresses(LOADS),transformed(sec,angle=angle).stresses(load)
    assert a.max_sigma_vm==pytest.approx(b.max_sigma_vm,rel=1e-10)
    for key,p in a.segment_profiles.items():
        q=b.segment_profiles[key]
        for xi in (.0,.3,.8,1):
            assert p.eval_sigma_vm(xi*p.length)==pytest.approx(q.eval_sigma_vm(xi*q.length),rel=1e-10)


@pytest.mark.parametrize("kind",["I","box","hat"])
def test_direction_and_order_invariance(kind):
    sec,*_=fixtures(kind)
    a,b=sec.stresses(LOADS),transformed(sec,reverse=True).stresses(LOADS)
    assert a.max_sigma_vm==pytest.approx(b.max_sigma_vm,rel=1e-11)
    for key,p in a.segment_profiles.items():
        q=b.segment_profiles[key]
        for xi in (0,.21,.5,1):
            assert p.eval_sigma_zz(xi*p.length)==pytest.approx(q.eval_sigma_zz((1-xi)*q.length),rel=1e-10,abs=1e-11)
            assert p.eval_tau_membrane(xi*p.length)==pytest.approx(-q.eval_tau_membrane((1-xi)*q.length),rel=1e-10,abs=1e-11)
            assert p.eval_sigma_vm(xi*p.length)==pytest.approx(q.eval_sigma_vm((1-xi)*q.length),rel=1e-10)


def square(cls=ClosedSection):
    points=[(-1,-1),(1,-1),(1,1),(-1,1)]
    return cls([Segment(Node(*points[i]),Node(*points[(i+1)%4]),.125,i) for i in range(4)])


@pytest.mark.parametrize("cls",[ClosedSection,MixedSection])
def test_zero_warping_and_zero_load(cls):
    sec=square(cls)
    with pytest.raises(SingularSectionError,match="Bimoment B cannot be sustained"):
        sec.stresses(AppliedLoads(B=1))
    with pytest.raises(SingularSectionError,match="Warping moment"):
        sec.stresses(AppliedLoads(M_omega=1))
    zero=sec.stresses(AppliedLoads(sigma_yield=250))
    assert zero.max_sigma_vm==0
    assert math.isinf(zero.load_factor)
    assert all(p.eval_sigma_vm(0)==0 for p in zero.segment_profiles.values())
    assert sec.stresses(AppliedLoads(Tsv=1)).max_sigma_vm>0


@pytest.mark.parametrize("name",["N","Vx","Vy","Mx","My","B","M_omega","Tsv"])
def test_nonfinite_loads_rejected(name):
    with pytest.raises(GeometryError,match="finite"):
        AppliedLoads(**{name:math.inf})


@pytest.mark.parametrize("yield_stress",[0,-1,math.inf,math.nan])
def test_invalid_yield(yield_stress):
    with pytest.raises(GeometryError):
        AppliedLoads(sigma_yield=yield_stress)


def test_type_ids_and_coordinate_validation():
    sec,*_=fixtures("I")
    with pytest.raises(TypeError):
        calculate_stresses(sec,object())
    bad=Section([Segment(e.p1,e.p2,e.t,'duplicate') for e in sec.segments])
    with pytest.raises(GeometryError,match="unique"):
        bad.stresses(AppliedLoads())
    p=sec.stresses(AppliedLoads(N=1)).segment_profiles[0]
    for s in (-1,math.inf,math.nan,p.length*1.01):
        with pytest.raises(GeometryError):
            p.eval_sigma_vm(s)


@pytest.mark.parametrize("loads",[AppliedLoads(N=sys.float_info.max),AppliedLoads(B=sys.float_info.max),AppliedLoads(Vy=sys.float_info.max),AppliedLoads(Tsv=sys.float_info.max),AppliedLoads(M_omega=sys.float_info.max)])
def test_physical_stress_overflow(loads):
    sec,*_=fixtures("hat")
    with pytest.raises(OverflowError,match="Recovered stress overflows IEEE-754 float64 range"):
        sec.stresses(loads)


def test_subnormal_preserved_and_below_range_rejected():
    sec=square()  # area=1, hence sigma=N exactly.
    value=math.ldexp(1.0,-1074)
    res=sec.stresses(AppliedLoads(N=value))
    assert res.max_sigma_vm==value
    for p in res.segment_profiles.values():
        assert p.eval_sigma_zz(0)==value
    big=ClosedSection([Segment(e.p1,e.p2,.5,e.id) for e in sec.segments])
    with pytest.raises(FloatingPointError,match="non-zero but underflows"):
        big.stresses(AppliedLoads(N=value))


def test_vm_avoids_squaring_overflow():
    p=SegmentStressProfile('huge',1,.1,(0,1e200),(0,0,1e200),0)
    assert p.peak_sigma_vm==pytest.approx(2e200,rel=1e-14)


def test_intermediate_bimoment_product_overflow_is_avoided():
    # omega~L² and Cw~L^6; B*omega overflows but recovered stress is finite.
    base,*_=fixtures('I')
    scale=1e10
    sec=Section([Segment(Node(e.p1.x*scale,e.p1.y*scale),Node(e.p2.x*scale,e.p2.y*scale),e.t*scale,e.id) for e in base.segments])
    loads=AppliedLoads(B=1e300)
    result=sec.stresses(loads)
    assert math.isfinite(result.max_sigma_vm)
    assert result.max_sigma_vm==pytest.approx(6e300/(.04*4**3*6**2/24)/scale**4,rel=1e-12)


def test_pure_torsion_does_not_request_overflowing_cw():
    # J~L^4 is finite here while Cw~L^6 is not. Pure Tsv needs only J.
    base,*_=fixtures('I')
    scale=1e52
    section=Section([Segment(Node(e.p1.x*scale,e.p1.y*scale),Node(e.p2.x*scale,e.p2.y*scale),e.t*scale,e.id) for e in base.segments])
    result=section.stresses(AppliedLoads(Tsv=1e150))
    expected_j=(6*.02**3+8*.04**3)/3
    assert result.max_sigma_vm==pytest.approx(math.sqrt(3)*.04*1e150/expected_j/scale**3,rel=1e-12)


@pytest.mark.parametrize("factor",[1e-8,1e8])
def test_geometric_scaling_at_fixed_stress(factor):
    section,*_=fixtures('hat')
    scaled=MixedSection([Segment(Node(e.p1.x*factor,e.p1.y*factor),Node(e.p2.x*factor,e.p2.y*factor),e.t*factor,e.id) for e in section.segments],node_tolerance=factor*1e-9)
    loads=replace(LOADS,N=LOADS.N*factor**2,Vx=LOADS.Vx*factor**2,Vy=LOADS.Vy*factor**2,Mx=LOADS.Mx*factor**3,My=LOADS.My*factor**3,Tsv=LOADS.Tsv*factor**3,M_omega=LOADS.M_omega*factor**3,B=LOADS.B*factor**4)
    ref,actual=section.stresses(LOADS),scaled.stresses(loads)
    assert actual.max_sigma_vm==pytest.approx(ref.max_sigma_vm,rel=1e-10)
    for key,p in actual.segment_profiles.items():
        q=ref.segment_profiles[key]
        assert p.eval_sigma_vm(.37*p.length)==pytest.approx(q.eval_sigma_vm(.37*q.length),rel=1e-10)


@pytest.mark.parametrize("kind",['I','box','hat'])
def test_secondary_flow_node_and_cycle_balance(kind):
    section,*_=fixtures(kind)
    result=section.stresses(AppliedLoads(M_omega=2))
    from thinwallx.mixed_topology import extract_mixed_topology
    topology=extract_mixed_topology(section.segments)
    residual=[0.0]*topology.node_count
    scale=0.0
    integrals=[]
    for i,(u,v) in enumerate(topology.canonical_edges):
        p=result.segment_profiles[i]
        start,end=p.eval_tau_membrane(0)*p.thickness,p.eval_tau_membrane(p.length)*p.thickness
        residual[u]-=start;residual[v]+=end
        scale=max(scale,abs(start),abs(end))
        a,b,c=p.tau_membrane_coeffs
        integrals.append(p.length*(a*p.length**2/3+b*p.length/2+c))
    assert max(abs(v) for v in residual)<=scale*1e-10
    for row in topology.B:
        integral=sum(float(s)*v for s,v in zip(row,integrals))
        natural=sum(abs(float(s)*v) for s,v in zip(row,integrals))
        assert abs(integral)<=natural*1e-10


@pytest.mark.parametrize("cls",[Section,MixedSection])
@pytest.mark.parametrize("kind",["angle","tee","cross"])
@pytest.mark.parametrize("scale",[1e-8,1.0,1e8])
@pytest.mark.parametrize("warping_load",["B","M_omega"])
def test_roundoff_only_warping_is_rejected(cls,kind,scale,warping_load):
    # All supporting straight lines meet at (0,0): omega=Cw=0 analytically.
    tips=[(100,0),(0,60)]
    if kind in ("tee","cross"):
        tips.append((-40,0))
    if kind == "cross":
        tips.append((0,-70))
    section=cls([
        Segment(Node(0,0),Node(x*scale,y*scale),2*scale,i)
        for i,(x,y) in enumerate(tips)
    ],node_tolerance=scale*1e-9)
    loads=AppliedLoads(N=12345,Mx=54321,My=-67890,**{warping_load:98765})
    with pytest.raises(SingularSectionError,match="cannot be sustained"):
        section.stresses(loads)
    # Ordinary bending remains meaningful even when warping resistance is zero.
    # Keep N : M/L fixed when changing length units; otherwise this check
    # measures cancellation between widely separated load components instead.
    result=section.stresses(AppliedLoads(N=12345,Mx=54321*scale,My=-67890*scale))
    assert result.resultant_N==pytest.approx(12345,rel=1e-10)
    assert result.resultant_Mx==pytest.approx(54321*scale,rel=1e-10,abs=0)
    assert result.resultant_My==pytest.approx(-67890*scale,rel=1e-10,abs=0)


def test_resolution_threshold_has_no_intermediate_range_loss():
    from fractions import Fraction
    from types import SimpleNamespace
    from thinwallx.stress import _unresolved_warping_resistance
    for length,thickness in [(1e70,1e-70),(1e-70,1e50)]:
        section=SimpleNamespace(segments=[SimpleNamespace(length=length,t=thickness)])
        exact=Fraction.from_float(1e-12)*Fraction.from_float(thickness)*Fraction.from_float(length)**5
        below=float(exact/2)
        above=float(exact*2)
        assert below>0 and math.isfinite(above)
        assert _unresolved_warping_resistance(section,below)
        assert not _unresolved_warping_resistance(section,above)


def test_resolved_subnormal_cw_is_not_a_singularity():
    # Uniform scaling preserves the dimensionless warping resistance. Cw is
    # subnormal here; its quantization allows ~1e-7 relative precision, not 1e-10.
    base,*_=fixtures('I')
    scale=math.ldexp(1.0,-175)
    section=Section([
        Segment(Node(e.p1.x*scale,e.p1.y*scale),Node(e.p2.x*scale,e.p2.y*scale),e.t*scale,e.id)
        for e in base.segments
    ],node_tolerance=scale*1e-9)
    assert 0<section.Cw<sys.float_info.min
    result=section.stresses(AppliedLoads(B=math.ldexp(1.0,-700)))
    from fractions import Fraction as Q
    expected=float(Q(6)/(Q(1,25)*4**3*6**2/24))
    assert result.max_sigma_vm==pytest.approx(expected,rel=1e-7,abs=0)
