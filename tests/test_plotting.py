"""T44--T59 headless rendering, numerical data, determinism and ownership."""
import hashlib
import os
from pathlib import Path
import struct
import subprocess
import sys
from dataclasses import replace
from xml.etree import ElementTree

import numpy as np
import pytest

from sectalix import AppliedLoads, Section, Segment, Node
from sectalix.plotting import (
    PlotStyle, GeometryPlotOptions, ShearPlotOptions, StressPlotOptions,
    plot_geometry, plot_shear_flow, plot_stresses, _frame, _sample_points, _flow_candidates,
)
from tests.test_serialization import shape

STYLE=PlotStyle(figsize=(6.,4.),dpi=60)


@pytest.mark.parametrize("fmt",["png","svg"])
@pytest.mark.parametrize("which",["geometry","shear","stress"])
def test_t45_t54_t55_outputs(tmp_path,fmt,which):
    sec=shape()
    def render(p):
        if which=="geometry":
            return plot_geometry(sec,p,options=GeometryPlotOptions(style=STYLE,show_node_ids=True,show_segment_ids=True,show_shear_center=True,show_principal_axes=True))
        if which=="shear":
            return plot_shear_flow(sec,sec.calculate_shear_flow(vx=0,vy=1),p,options=ShearPlotOptions(style=STYLE))
        return plot_stresses(sec,sec.calculate_stresses(AppliedLoads(N=2)),p,options=StressPlotOptions(style=STYLE))
    a,b=tmp_path/f"a.{fmt}",tmp_path/f"b.{fmt}"
    report=render(a)
    render(b)
    assert a.read_bytes()==b.read_bytes()
    assert report.segment_count==2
    if fmt=="png":
        raw=a.read_bytes()
        assert raw[:8]==b"\x89PNG\r\n\x1a\n"
        assert struct.unpack(">II",raw[16:24])==(360,240)
    else:
        root=ElementTree.parse(a).getroot()
        assert root.tag.endswith("svg")
        assert "viewBox" in root.attrib


def test_t46_t52_local_coordinates():
    sec=shape()
    translated=sec.translated(1e12,1e12)
    for i in range(2):
        assert np.array_equal(_sample_points(sec,_frame(sec),i,[0.,1.]),
                              _sample_points(translated,_frame(translated),i,[0.,1.]))


def test_t47_isotropic(tmp_path):
    from tests.test_serialization import shape
    from sectalix import ClosedSection
    p=[Node(0,0),Node(1,0),Node(1,1),Node(0,1)]
    sec=ClosedSection([Segment(p[i],p[(i+1)%4],.01) for i in range(4)])
    report=plot_geometry(sec,tmp_path/"iso.svg",options=GeometryPlotOptions(style=STYLE,show_principal_axes=True))
    assert any("isotropic" in s for s in report.skipped_annotations)


def test_t48_reverse_flow():
    sec=shape()
    rev=Section([s.reversed() for s in sec.segments])
    a,b=sec.calculate_shear_flow(vx=1,vy=2),rev.calculate_shear_flow(vx=1,vy=2)
    for x,y in zip(a.segment_flows,b.segment_flows):
        np.testing.assert_allclose(x.vector_at_xi(.25),y.vector_at_xi(.75),rtol=1e-12,atol=0)
        assert _flow_candidates(x)[0]==0
        assert _flow_candidates(x)[-1]==1


def test_t49_closed_plot(tmp_path):
    from sectalix import ClosedSection
    sec=shape(ClosedSection)
    plot_shear_flow(sec,sec.calculate_shear_flow(vx=0,vy=1),tmp_path/"closed.svg",options=ShearPlotOptions(style=STYLE))


@pytest.mark.parametrize("samples",[17,65,257])
def test_t50_t51_t53_peak_and_profiles(tmp_path,samples):
    sec=shape()
    result=sec.calculate_stresses(AppliedLoads(N=2,Mx=1))
    maximum=result.max_sigma_vm
    report=plot_stresses(sec,result,tmp_path/f"{samples}.svg",
                         options=StressPlotOptions(style=STYLE,samples_per_segment=samples,abscissa="concatenated"))
    assert result.max_sigma_vm==maximum
    assert report.field_scale==float(maximum).hex()
    assert "NOT continuous perimeter" in report.path.read_text(encoding="utf-8")


def test_t57_rcparams_and_other_figures(tmp_path):
    import matplotlib
    from matplotlib import pyplot as plt
    fig=plt.figure()
    before=dict(matplotlib.rcParams)
    plot_geometry(shape(),tmp_path/"style.png",options=GeometryPlotOptions(style=STYLE))
    assert dict(matplotlib.rcParams)==before
    assert plt.fignum_exists(fig.number)
    plt.close(fig)


def test_t56_repeated_cleanup(tmp_path):
    from matplotlib._pylab_helpers import Gcf
    before=len(Gcf.get_all_fig_managers())
    for i in range(100):
        plot_geometry(shape(),tmp_path/"repeat.svg",options=GeometryPlotOptions(style=STYLE),overwrite=True)
    assert len(Gcf.get_all_fig_managers())==before


def test_t58_title_and_t59_overwrite(tmp_path):
    sec=shape()
    p=tmp_path/"safe.svg"
    plot_geometry(sec,p,options=GeometryPlotOptions(style=replace(STYLE,title='<script>& $x$')))
    text=p.read_text(encoding="utf-8")
    assert "<script>" not in text
    assert "&lt;script&gt;" in text
    with pytest.raises(FileExistsError):
        plot_geometry(sec,p)
    with pytest.raises(ValueError):
        plot_geometry(sec,tmp_path/"bad.pdf")


def test_t44_no_matplotlib_imported():
    code="import sys; import sectalix; import sectalix.serialization; import sectalix.dxf; import sectalix.plotting; assert 'matplotlib' not in sys.modules"
    env=dict(os.environ,PYTHONPATH=str(Path("src").absolute()))
    subprocess.run([sys.executable,"-W","error","-c",code],env=env,check=True,capture_output=True)


@pytest.mark.parametrize("fmt", ["svg","png"])
def test_t54_clean_processes(tmp_path,fmt):
    code="""from sectalix import Section,Segment,Node
from sectalix.plotting import plot_geometry
import sys
s=Section([Segment(Node(0,0),Node(1,1),.01)])
plot_geometry(s,sys.argv[1])
"""
    env=dict(os.environ,PYTHONPATH=str(Path("src").absolute()),MPLBACKEND="Agg")
    for name in (f"a.{fmt}",f"b.{fmt}"):
        subprocess.run([sys.executable,"-W","error","-c",code,str(tmp_path/name)],env=env,check=True,capture_output=True)
    assert (tmp_path/f"a.{fmt}").read_bytes()==(tmp_path/f"b.{fmt}").read_bytes()


def test_t59_failed_save_cleans(tmp_path,monkeypatch):
    from matplotlib.figure import Figure
    from matplotlib._pylab_helpers import Gcf
    before=len(Gcf.get_all_fig_managers())
    def fail(*args,**kwargs):
        raise OSError("simulated full disk")
    monkeypatch.setattr(Figure,"savefig",fail)
    with pytest.raises(OSError):
        plot_geometry(shape(),tmp_path/"fail.png",options=GeometryPlotOptions(style=STYLE))
    assert not (tmp_path/"fail.png").exists()
    assert not list(tmp_path.glob(".sectalix-*"))
    assert len(Gcf.get_all_fig_managers())==before


@pytest.mark.parametrize("value", [5e-324, sys.float_info.max])
def test_t52_extreme_field(tmp_path,value):
    from sectalix import SegmentStressProfile, StressRecoveryResult
    sec=shape()
    profiles={s.id:SegmentStressProfile(s.id,s.length,s.t,(0.,value),(0.,0.,0.),0.) for s in sec.segments}
    result=StressRecoveryResult(profiles,value,0,0.,sec.segments[0].p1.coords,None,0.,0.,0.,0.)
    report=plot_stresses(sec,result,tmp_path/"extreme.svg",options=StressPlotOptions(style=STYLE))
    assert report.field_scale==value.hex()


@pytest.mark.parametrize("exponent", [-200,200])
def test_t52_extreme_geometry(tmp_path,exponent):
    scale=2.**exponent
    sec=Section([Segment(Node(0.,0.),Node(scale,0.),.01*scale),
                 Segment(Node(scale,0.),Node(scale,scale),.01*scale)],node_tolerance=0.)
    report=plot_geometry(sec,tmp_path/"geometry.svg",options=GeometryPlotOptions(style=STYLE))
    assert abs(report.display_length_exponent-exponent)<=2


def test_t44_missing_matplotlib(tmp_path):
    code="""import sys, importlib.abc
class Block(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'matplotlib' or fullname.startswith('matplotlib.'):
            raise ModuleNotFoundError('blocked matplotlib')
sys.meta_path.insert(0, Block())
from sectalix import Section, Segment, Node
from sectalix.serialization import to_json
from sectalix.plotting import plot_geometry
s=Section([Segment(Node(0,0),Node(1,1),.01)])
assert to_json(s)
try:
    plot_geometry(s,sys.argv[1])
except ImportError as exc:
    assert 'optional matplotlib' in str(exc)
else:
    raise AssertionError('plot unexpectedly succeeded')
"""
    env=dict(os.environ,PYTHONPATH=str(Path("src").absolute()))
    subprocess.run([sys.executable,"-W","error","-c",code,str(tmp_path/"missing.png")],
                   env=env,check=True,capture_output=True)
    assert not (tmp_path/"missing.png").exists()
