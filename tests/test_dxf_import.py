"""T16--T41: restricted DXF syntax and explicit topology normalization."""
from dataclasses import replace
import math
from pathlib import Path

import pytest

from sectalix import Section, ClosedSection, MixedSection, GeometryError, TopologyError
from sectalix.dxf import DxfImportOptions, ThicknessMap, loads_dxf, read_dxf


def options(**kw):
    return DxfImportOptions("mm", "mm", ThicknessMap(default=.01), **kw)


def document(entities, version="AC1015", extra_header=()):
    pairs = [(0,"SECTION"),(2,"HEADER"),(9,"$ACADVER"),(1,version), *extra_header,
             (0,"ENDSEC"),(0,"SECTION"),(2,"ENTITIES"),*entities,(0,"ENDSEC"),(0,"EOF")]
    return "".join(f"{c}\n{v}\n" for c,v in pairs)


def line(a, b, extra=()):
    return [(0,"LINE"),(8,"0"),(10,a[0]),(20,a[1]),(11,b[0]),(21,b[1]),*extra]


def lw(points, closed=False, extra=()):
    return [(0,"LWPOLYLINE"),(8,"0"),(90,len(points)),(70,int(closed)),*extra,
            *(item for p in points for item in [(10,p[0]),(20,p[1])])]


def poly(points, closed=False):
    return [(0,"POLYLINE"),(8,"0"),(66,1),(70,int(closed)),
            *(item for p in points for item in [(0,"VERTEX"),(10,p[0]),(20,p[1])]),
            (0,"SEQEND")]


def signature(section):
    return [(s.p1.coords, s.p2.coords, s.t, s.id) for s in section.segments]


@pytest.mark.parametrize("kind,version", [("line","AC1009"),("poly","AC1009"),("lw","AC1015")])
def test_t16_t18_entities(kind,version):
    points=[(0,0),(2,0),(2,1)]
    entities = line(points[0],points[1])+line(points[1],points[2]) if kind=="line" else poly(points) if kind=="poly" else lw(points)
    result=loads_dxf(document(entities,version), options=options())
    assert type(result.section) is Section
    assert result.report.edge_count==2
    assert result.report.cycle_rank==0


@pytest.mark.parametrize("version", ["AC1009","AC1015"])
def test_t17_closed(version):
    result=loads_dxf(document(poly([(0,0),(2,0),(2,1),(0,1)],True),version),options=options())
    assert type(result.section) is ClosedSection
    assert result.report.cycle_rank==1


@pytest.mark.parametrize("text", ["AutoCAD Binary DXF\r\n",document(line((0,0),(1,0)),"AC1027"),
                                 document(lw([(0,0),(1,0)]),"AC1009"),
                                 document(line((0,0),(1,0))).replace("0\nEOF\n",""),
                                 document([(0,"VERTEX"),(10,0),(20,0)]),
                                 document([(0,"POLYLINE"),(0,"VERTEX"),(10,0),(20,0)]),
                                 document(line((0,0),(1,0)))+"junk\n"])
def test_t19_t20_invalid_syntax(text):
    with pytest.raises((ValueError,GeometryError)):
        loads_dxf(text,options=options())


@pytest.mark.parametrize("extra", [[(42,"1e-300")],[(43,1)],[(39,1)],[(210,1)],[(230,-1)],[(38,1)]])
def test_t21_t24_reject_geometry(extra):
    with pytest.raises(GeometryError):
        loads_dxf(document(lw([(0,0),(1,0)],extra=extra)),options=options())


@pytest.mark.parametrize("kind", ["ARC","SPLINE","INSERT","3DFACE"])
def test_t21_unsupported(kind):
    with pytest.raises(GeometryError):
        loads_dxf(document([(0,kind),(8,"0")]),options=options())


def test_t22_plane():
    text=document(line((0,0),(1,0),[(30,2),(31,2)]))
    assert loads_dxf(text,options=options(plane_z=2)).report.max_plane_displacement==0
    with pytest.raises(GeometryError):
        loads_dxf(text,options=options())


def test_t23_thickness():
    text=document(line((0,0),(1,0),[(5,"A")])).replace("8\n0\n","8\nTHICK_5.0\n")
    assert loads_dxf(text,options=options()).section.segments[0].t==5
    opts=replace(options(),thickness=ThicknessMap(by_handle={"a":2.}))
    assert loads_dxf(text,options=opts).section.segments[0].t==2
    with pytest.raises(ValueError,match="Unused"):
        loads_dxf(text,options=replace(options(),thickness=ThicknessMap(by_layer={"absent":2.})))
    with pytest.raises((ValueError,GeometryError)):
        loads_dxf(text.replace("THICK_5.0","THICK_-2"),options=options())


def test_t24_lineweight_is_not_t():
    result=loads_dxf(document(line((0,0),(1,0),[(370,100)])),options=options())
    assert result.section.segments[0].t==.01


def test_t25_units():
    text=document(line((0,0),(1000,0)),extra_header=[(9,"$INSUNITS"),(70,4)])
    out=loads_dxf(text,options=replace(options(),target_length_unit="m",node_tolerance=1e-12))
    assert out.section.segments[0].length==1
    assert out.section.segments[0].t==pytest.approx(1e-5,rel=1e-15,abs=0)
    with pytest.raises(ValueError):
        loads_dxf(text,options=replace(options(),source_length_unit="m"))


def test_t26_t28_invariance():
    edges=[((0,0),(2,0)),((2,0),(2,1))]
    a=loads_dxf(document(sum((line(p,q) for p,q in edges),[])),options=options()).section
    b=loads_dxf(document(sum((line(q,p) for p,q in reversed(edges)),[])),options=options()).section
    assert signature(a)==signature(b)
    offset=1e12
    c=loads_dxf(document(sum((line((p[0]+offset,p[1]+offset),(q[0]+offset,q[1]+offset)) for p,q in edges),[])),options=options()).section
    assert c.Ix==pytest.approx(a.Ix,rel=1e-10,abs=0)


@pytest.mark.parametrize("token,error", [("1e-999",FloatingPointError),("1e999",OverflowError),("nan",ValueError)])
def test_t29_conversion_range(token,error):
    with pytest.raises(error):
        loads_dxf(document(line((0,0),(token,0))),options=options(node_tolerance=0))


@pytest.mark.parametrize("gap", [0.,.5e-6,1e-6])
def test_t30_cluster(gap):
    text=document(line((-1,0),(0,0))+line((gap,0),(0,1)))
    result=loads_dxf(text,options=options(node_tolerance=1e-6))
    assert result.report.node_count==3
    assert result.report.max_endpoint_displacement<=1e-6


def test_t31_chain():
    text=document(line((-2,0),(0,0))+line((.75,0),(.75,3))+line((1.5,0),(4,0)))
    with pytest.raises(TopologyError,match="transitive"):
        loads_dxf(text,options=options(node_tolerance=1))


@pytest.mark.parametrize("second", [line((1,0),(0,0)),line((.5,0),(2,0))])
def test_t32_duplicate_overlap(second):
    with pytest.raises(TopologyError):
        loads_dxf(document(line((0,0),(1,0))+second),options=options())


def test_t33_t34_t_split():
    text=document(line((-1,0),(1,0))+line((0,0),(0,1)))
    result=loads_dxf(text,options=options())
    assert result.report.edge_count==3
    assert any(c["operation"]=="split" for c in result.report.changes)
    assert result.section.area==pytest.approx(.03,rel=1e-15,abs=0)
    with pytest.raises(TopologyError):
        loads_dxf(text,options=options(junction_policy="reject"))
    near=document(line((-1,0),(1,0))+line((0,1e-10),(0,1)))
    assert loads_dxf(near,options=options()).report.edge_count==3


def test_t35_x():
    with pytest.raises(TopologyError,match="X crossing"):
        loads_dxf(document(line((-1,0),(1,0))+line((0,-1),(0,1))),options=options())


def test_t36_mixed_and_explicit_closed():
    text=document(lw([(0,0),(1,0),(1,1),(0,1)],True)+line((1,1),(2,1)))
    assert type(loads_dxf(text,options=options()).section) is MixedSection
    with pytest.raises((GeometryError,TopologyError)):
        loads_dxf(text,options=options(section_kind="closed"))


def test_t38_closing_repeat():
    result=loads_dxf(document(lw([(0,0),(1,0),(1,1),(0,1),(0,0)],True)),options=options())
    assert result.report.removed_closing_vertices==1
    assert result.report.edge_count==4


def test_t39_layer_filter():
    text=document(line((0,0),(1,0))+[(0,"INSERT"),(8,"IGNORE")])
    out=loads_dxf(text,options=options(layers=frozenset({"0"})))
    assert out.report.excluded_entities=={"INSERT":1}


def test_t40_duplicate_handle():
    with pytest.raises(ValueError,match="handle"):
        loads_dxf(document(line((0,0),(1,0),[(5,"A")])+line((1,0),(1,1),[(5,"a")])),options=options())


def test_t40_read_file(tmp_path):
    p=tmp_path/"simple.dxf"
    p.write_text(document(line((0,0),(1,0))),encoding="ascii")
    assert read_dxf(p,options=options()).report.edge_count==1


def test_t41_limits():
    with pytest.raises(ValueError,match="budget"):
        loads_dxf(document(line((0,0),(1,0))),options=options(max_pair_checks=1))


def test_t20_duplicate_coordinate():
    with pytest.raises(ValueError, match="duplicate"):
        loads_dxf(document(line((0,0),(1,0),[(10,2)])),options=options())


def test_t21_subnormal_bulge_not_zero():
    with pytest.raises(GeometryError):
        loads_dxf(document(lw([(0,0),(1,0)],extra=[(42,"1e-999")])),options=options())


def test_t22_polyline_3d_flag():
    text=document(poly([(0,0),(1,0)])).replace("70\n0\n","70\n8\n")
    with pytest.raises(GeometryError):
        loads_dxf(text,options=options())


def test_t30_nextafter_outside_tolerance():
    text=document(line((-1,0),(0,0))+line((math.nextafter(1e-6,math.inf),0),(1,1)))
    with pytest.raises(TopologyError):
        loads_dxf(text,options=options(node_tolerance=1e-6))


def test_t34_ambiguous_projection():
    text=document(line((-1,0),(1,0))+line((2e-6,-1),(2e-6,1))+line((1e-6,1e-6),(1,2)))
    with pytest.raises(TopologyError,match="incompatible"):
        loads_dxf(text,options=options(node_tolerance=1.1e-6))


def test_t36_disconnected():
    with pytest.raises(TopologyError):
        loads_dxf(document(line((0,0),(1,0))+line((3,0),(4,0))),options=options())


def test_t37_star_junction():
    text=document(lw([(0,0),(1,0),(1,1),(0,1)],True)+
                  line((1,1),(2,1))+line((1,1),(1,2))+line((1,1),(2,2)))
    result=loads_dxf(text,options=options())
    assert result.report.edge_count==7
    assert result.report.cycle_rank==1
    assert type(result.section) is MixedSection


def test_t39_paperspace_report_and_conflict():
    text=document(line((0,0),(1,0))+[(0,"CIRCLE"),(67,1),(410,"Sheet")])
    assert loads_dxf(text,options=options()).report.excluded_entities=={"CIRCLE":1}
    with pytest.raises(ValueError):
        loads_dxf(text.replace("410\nSheet","410\nModel"),options=options())


def test_t40_legacy_encoding(tmp_path):
    text=document(line((0,0),(1,0)),extra_header=[(9,"$DWGCODEPAGE"),(3,"ANSI_1254")])
    text=text.replace("8\n0\n","8\nÇELİK\n")
    p=tmp_path/"legacy.dxf"
    p.write_bytes(text.encode("cp1254"))
    assert read_dxf(p,options=options(),encoding="cp1254").report.edge_count==1
    with pytest.raises(UnicodeDecodeError):
        read_dxf(p,options=options())
