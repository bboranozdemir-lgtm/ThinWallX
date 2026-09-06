"""Independent 100-digit direct-equilibrium stress oracles (no production solves).

The v0.6 test oracle assembles geometry from literal nodes/edges/cells and uses
Decimal Gaussian elimination, without virtual cuts or production property calls.
Here all warping, secondary-flow and stress reference arithmetic stays Decimal.
"""

from decimal import Decimal as D, localcontext
import math

import pytest

from sectalix import AppliedLoads, ClosedSection, MixedSection, Node, Section, Segment
from tests.test_mixed_benchmarks import IndependentMixedOracle, _decimal_solve


def fixtures(kind):
    if kind == "I":
        nodes = [(0,-3),(0,3),(-2,3),(2,3),(-2,-3),(2,-3)]
        edges = [(0,1,.02),(2,1,.04),(1,3,.04),(4,0,.04),(0,5,.04)]
        cells, cls = [], Section
    elif kind == "asymmetric_channel":
        nodes = [(0,0),(0,6),(4,6),(2,0)]
        edges = [(0,1,.02),(1,2,.04),(0,3,.03)]
        cells, cls = [], Section
    else:
        nodes = [(0,0),(4,0),(4,2),(0,2)]
        edges = [(0,1,.02),(1,2,.02),(2,3,.02),(3,0,.02)]
        cells = [((0,1,2,3),(1,1,1,1))]
        cls = ClosedSection
        if kind in ("hat", "asymmetric_hat"):
            nodes += [(-1,0),(5,0)]
            edges += [(0,4,.01),(1,5,.01)]
            cls = MixedSection
            if kind == "asymmetric_hat":
                nodes[-1] = (6,1)
                edges[-1] = (1,5,.03)
    section = cls([Segment(Node(*nodes[u]),Node(*nodes[v]),t,id=i) for i,(u,v,t) in enumerate(edges)])
    return section, nodes, edges, cells


def reference(nodes, edges, cells, loads):
    with localcontext() as ctx:
        ctx.prec = 100
        o = IndependentMixedOracle(nodes, edges, cells)
        cv, _, _ = o._solve_shear_flow_decimal(D(str(loads.Vx)),D(str(loads.Vy)))
        ax, ay = _decimal_solve([[o.Iy,o.Ixy],[o.Ixy,o.Ix]],[D(str(loads.Vx)),D(str(loads.Vy))])
        phi = _decimal_solve(o.H,[2*a for a in o.cell_areas]) if cells else []
        F = [D(0)]*len(edges)
        closed = set()
        for k,(ce,signs) in enumerate(cells):
            for i,sign in zip(ce,signs):
                F[i] += sign*phi[k]
                closed.add(i)
        J = 2*sum(a*p for a,p in zip(o.cell_areas,phi)) + sum(o.edge_lengths[i]*o.edges[i][2]**3/3 for i in range(len(edges)) if i not in closed)
        _,_,xs,ys = o._compute_shear_center_decimal()
        w = {0:D(0)}
        while len(w)<len(nodes):
            for i,(u,v,t) in enumerate(o.edges):
                x,y = o.nodes[u]
                dx,dy = o.nodes[v][0]-x,o.nodes[v][1]-y
                delta = (x-xs)*dy-(y-ys)*dx-F[i]*o.edge_lengths[i]/t
                if u in w and v not in w:
                    w[v]=w[u]+delta
                elif v in w and u not in w:
                    w[u]=w[v]-delta
        mean = sum(t*o.edge_lengths[i]*(w[u]+w[v])/2 for i,(u,v,t) in enumerate(o.edges))/o.area
        w = {i:value-mean for i,value in w.items()}
        Cw = sum(t*o.edge_lengths[i]*(w[u]**2+w[u]*w[v]+w[v]**2)/3 for i,(u,v,t) in enumerate(o.edges))
        # Direct physical-edge secondary flow, independent of the tree algorithm.
        a1 = [-D(str(loads.M_omega))*t*o.edge_lengths[i]*w[u]/Cw for i,(u,v,t) in enumerate(o.edges)]
        a2 = [-D(str(loads.M_omega))*t*o.edge_lengths[i]*(w[v]-w[u])/(2*Cw) for i,(u,v,t) in enumerate(o.edges)]
        matrix=[]; rhs=[]
        for node in range(1,len(nodes)):
            matrix.append([D(int(v==node)-int(u==node)) for u,v,t in edges])
            rhs.append(-sum(a1[i]+a2[i] for i,(u,v,t) in enumerate(edges) if v==node))
        for ce,signs in cells:
            row=[D(0)]*len(edges); b=D(0)
            for i,sign in zip(ce,signs):
                row[i]=sign*o.edge_lengths[i]/o.edges[i][2]
                b-=row[i]*(a1[i]/2+a2[i]/3)
            matrix.append(row);rhs.append(b)
        constants=_decimal_solve(matrix,rhs)
        result=[]
        for i,(u,v,t) in enumerate(o.edges):
            X,Y=o.nodes[u][0]-o.cx,o.nodes[u][1]-o.cy
            dx,dy=o.nodes[v][0]-o.nodes[u][0],o.nodes[v][1]-o.nodes[u][1]
            L=o.edge_lengths[i]
            det=o.Ix*o.Iy-o.Ixy**2
            bx=(D(str(loads.My))*o.Ix+D(str(loads.Mx))*o.Ixy)/det
            by=-(D(str(loads.Mx))*o.Iy+D(str(loads.My))*o.Ixy)/det
            surface=abs(D(str(loads.Tsv)))*t/J if i not in closed else D(0)
            samples=[]
            for xi in (D(0),D('0.23'),D('0.5'),D(1)):
                sigma=D(str(loads.N))/o.area+bx*(X+dx*xi)+by*(Y+dy*xi)+D(str(loads.B))*(w[u]+(w[v]-w[u])*xi)/Cw
                qv=cv[i]-t*L*((ax*X+ay*Y)*xi+(ax*dx+ay*dy)*xi**2/2)
                qw=constants[i]+a1[i]*xi+a2[i]*xi**2
                tau=(qv+qw+D(str(loads.Tsv))*F[i]/J)/t
                vm=(sigma**2+3*(abs(tau)+surface)**2).sqrt()
                samples.append((float(xi),float(sigma),float(tau),float(abs(tau)+surface),float(vm)))
            result.append(samples)
        return result


@pytest.mark.parametrize("kind",["I","box","hat","asymmetric_channel","asymmetric_hat"])
@pytest.mark.parametrize("loads",[
    AppliedLoads(B=7), AppliedLoads(Vx=3,Vy=-4), AppliedLoads(Tsv=-5),
    AppliedLoads(M_omega=2),
    AppliedLoads(N=2,Mx=3,My=-4,B=5,Vx=-6,Vy=7,Tsv=-8,M_omega=9,sigma_yield=250),
])
def test_independent_stress_oracle(kind,loads):
    section,nodes,edges,cells=fixtures(kind)
    if kind.startswith("asymmetric"):
        oracle = IndependentMixedOracle(nodes, edges, cells)
        assert abs(oracle.Ixy) > D("0.01")
        assert abs(section.Ixy) > 0.01
    expected=reference(nodes,edges,cells,loads)
    actual=section.stresses(loads)
    scale=max(abs(value) for samples in expected for row in samples for value in row[1:])
    for i,samples in enumerate(expected):
        p=actual.segment_profiles[i]
        for xi,sigma,tau,surface,vm in samples:
            values=(p.eval_sigma_zz(xi*p.length),p.eval_tau_membrane(xi*p.length),p.eval_tau_surface(xi*p.length),p.eval_sigma_vm(xi*p.length))
            assert values == pytest.approx((sigma,tau,surface,vm),rel=1e-10,abs=scale*1e-12)
    assert actual.resultant_N == pytest.approx(loads.N,abs=1e-11)
    assert actual.resultant_Mx == pytest.approx(loads.Mx,abs=1e-11)
    assert actual.resultant_My == pytest.approx(loads.My,abs=1e-11)
    assert actual.resultant_B == pytest.approx(loads.B,abs=1e-11)
    if loads.sigma_yield:
        assert actual.load_factor == pytest.approx(loads.sigma_yield/actual.max_sigma_vm,rel=1e-14)


def test_i_tip_bimoment_analytical_fraction():
    from fractions import Fraction as Q
    section,*_=fixtures("I")
    cw=Q(1,25)*4**3*6**2/24
    omega_tip=Q(4*6,4)
    expected=float(7*omega_tip/cw)
    r=section.stresses(AppliedLoads(B=7))
    assert r.max_sigma_vm == pytest.approx(expected,rel=1e-13)
    assert r.peak_location_xy[0] in (-2,2)
    assert r.peak_location_xy[1] in (-3,3)


def test_box_pure_torque_fraction():
    from fractions import Fraction as Q
    section,*_=fixtures("box")
    r=section.stresses(AppliedLoads(Tsv=5))
    expected=float(Q(5)/(2*8*Q(1,50)))
    for p in r.segment_profiles.values():
        assert p.eval_tau_membrane(p.length/2)==pytest.approx(expected,rel=1e-13)
        assert p.tau_sv_surface==0
    assert r.max_sigma_vm==pytest.approx(math.sqrt(3)*expected,rel=1e-13)
