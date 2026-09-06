"""Independent Fraction oracles from v0.8 §5.2, not production properties."""
from fractions import Fraction as Q
from dataclasses import replace
import math

import pytest
from thinwallx import AppliedLoads
from thinwallx.dxf import loads_dxf
from tests.test_dxf_import import document, line, lw, poly, options


@pytest.mark.parametrize("builder,version",[(poly,"AC1009"),(lw,"AC1015")])
def test_t42_rectangle(builder,version):
    b,h,t=Q(4),Q(2),Q(1,100)
    text=document(builder([(-2,-1),(2,-1),(2,1),(-2,1)],True),version)
    sec=loads_dxf(text,options=options()).section
    area=2*t*(b+h)
    ix=t*h*h*(b/2+h/6)
    iy=t*b*b*(h/2+b/6)
    j=2*t*b*b*h*h/(b+h)
    for got,expected in [(sec.area,area),(sec.Ix,ix),(sec.Iy,iy),(sec.J,j)]:
        assert got==pytest.approx(float(expected),rel=1e-12,abs=0)
    result=sec.calculate_stresses(AppliedLoads(N=3,Mx=2,My=1,Tsv=1))
    for seg in sec.segments:
        profile=result.segment_profiles[seg.id]
        for u in (Q(0),Q(1,4),Q(1)):
            x=Q(seg.p1.x)+(Q(seg.p2.x)-Q(seg.p1.x))*u
            y=Q(seg.p1.y)+(Q(seg.p2.y)-Q(seg.p1.y))*u
            sigma=3/area+x/iy-2*y/ix
            tau=1/(2*b*h*t)
            expected=math.sqrt(float(sigma*sigma+3*tau*tau))
            assert profile.eval_sigma_vm(float(u)*seg.length)==pytest.approx(expected,rel=1e-12,abs=0)


def test_t16_open_l_inertia_and_stress():
    b,h,t=Q(2),Q(1),Q(1,100)
    sec=loads_dxf(document(line((0,0),(2,0))+line((0,0),(0,1))),options=options()).section
    a=t*(b+h)
    cx,cy=b*b/(2*(b+h)),h*h/(2*(b+h))
    ix,iy=t*h**3/3-a*cy**2,t*b**3/3-a*cx**2
    ixy=-a*cx*cy
    det=ix*iy-ixy*ixy
    assert sec.Ixy==pytest.approx(float(ixy),rel=1e-12,abs=0)
    result=sec.calculate_stresses(AppliedLoads(N=1,Mx=2,My=3))
    for seg in sec.segments:
        p=result.segment_profiles[seg.id]
        x,y=Q(seg.p1.x)-cx,Q(seg.p1.y)-cy
        sigma=1/a+(3*ix+2*ixy)/det*x-(2*iy+3*ixy)/det*y
        assert p.eval_sigma_zz(0)==pytest.approx(float(sigma),rel=1e-12,abs=0)


def test_t43_open_q_oracle():
    sec=loads_dxf(document(line((0,0),(1,0))+line((0,0),(0,1))),options=options()).section
    result=sec.calculate_shear_flow(vx=0,vy=1)
    for seg,flow in zip(sec.segments,result.segment_flows):
        for s in (Q(0),Q(1,4),Q(1,2),Q(1)):
            expected=(-Q(3,4)+3*s-Q(9,4)*s*s) if seg.p2.x else (Q(3,4)+3*s-Q(15,4)*s*s)
            if expected:
                assert flow.q_at_xi(float(s))==pytest.approx(float(expected),rel=1e-12,abs=0)
            else:
                assert abs(flow.q_at_xi(float(s)))<=1e-12


def test_t37_barbell():
    entities=lw([(0,0),(1,0),(1,1),(0,1)],True)
    entities+=lw([(3,0),(4,0),(4,1),(3,1)],True)+line((1,0),(3,0))
    sec=loads_dxf(document(entities),options=options()).section
    t=Q(1,100)
    expected=2*t+Q(2,3)*t**3
    assert sec.J==pytest.approx(float(expected),rel=1e-12,abs=0)


@pytest.mark.parametrize("cut",[Q(1,6),Q(1,4),Q(9,20),Q(1,2)])
def test_t43_closed_cut_oracle(cut):
    from thinwallx import calculate_closed_shear_flow
    b,h,t,vy=Q(4),Q(2),Q(1,100),Q(1)
    sec=loads_dxf(document(lw([(-2,-1),(2,-1),(2,1),(-2,1)],True)),options=options()).section
    bottom=next(i for i,s in enumerate(sec.segments) if s.p1.y==s.p2.y==-1)
    result=calculate_closed_shear_flow(sec,vx=0.,vy=float(vy),cut_param=float(cut),
                                     spanning_tree_edges=[i for i in range(4) if i!=bottom])
    ix=t*h*h*(b/2+h/6)
    expected=-(vy*t*b*h*(1-2*cut))/(4*ix)
    got=float(result.redundant_cell_flows[0])
    if expected:
        assert got==pytest.approx(float(expected),rel=1e-12,abs=0)
    else:
        assert abs(got)<=1e-12
    sf=result.segment_flows[bottom]
    left=sf.q_at_xi(math.nextafter(float(cut),0.))
    right=sf.q_at_xi(math.nextafter(float(cut),1.))
    assert abs(left-right)<=1e-10*float(vy/h)
    assert result.recovered_resultant[1]==pytest.approx(1.,rel=1e-12,abs=0)
