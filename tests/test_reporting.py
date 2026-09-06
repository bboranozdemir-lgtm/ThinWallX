"""v1 T21--T30: calculation sheet structure, values, safe links and writes."""
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

import pytest

from sectalix import AppliedLoads, ClosedSection, MixedSection
from sectalix.reporting import calculation_report, write_report, inspect_section
from sectalix.serialization import UnitSystem
from tests.test_serialization import shape


@pytest.mark.parametrize("cls",[None,ClosedSection,MixedSection])
def test_t21_t26_t29_t30_report(cls):
    sec=shape() if cls is None else shape(cls)
    loads=AppliedLoads(N=1,Mx=.1,sigma_yield=250.)
    result=sec.calculate_stresses(loads)
    text=calculation_report(sec,loads,result,units=UnitSystem("mm","N"),
                            input_path="section.json",generated_at=datetime(2026,1,1,tzinfo=timezone.utc))
    for heading in ("Geometry and model","Section characteristics","Applied loads","Critical stress",
                    "Recovered resultants","Graphical appendices"):
        assert "## "+heading in text
    assert "2026-01-01T00:00:00+00:00" in text
    assert "section.json" in text and "force=N" in text
    for value in (sec.area,sec.Ixy,sec.theta_p,result.max_sigma_vm,result.load_factor,result.resultant_N):
        assert format(float(value),".17g") in text
    assert "not necessarily uniform" in text
    assert "not a design-code safety factor" in text


def test_t27_relative_links_and_escaping(tmp_path):
    sec=shape()
    loads=AppliedLoads(N=1)
    image=tmp_path/"a (b).png"
    image.write_bytes(b"image-fixture")
    report=tmp_path/"report.md"
    text=calculation_report(sec,loads,sec.calculate_stresses(loads),report_path=report,
                            input_path="<script>|bad\n# title",plots={"Safe":image})
    assert "](./a%20%28b%29.png)" in text
    assert "<script>" not in text and "\n# title" not in text
    assert "&#124;" in text


def test_t28_atomic(tmp_path,monkeypatch):
    sec=shape()
    loads=AppliedLoads(N=1)
    result=sec.calculate_stresses(loads)
    path=tmp_path/"report.md"
    write_report(sec,loads,result,path)
    old=path.read_bytes()
    with pytest.raises(FileExistsError):
        write_report(sec,loads,result,path)
    assert path.read_bytes()==old
    assert not list(tmp_path.glob(".sectalix-*"))


def test_zero_stress_unbounded_factor():
    sec=shape()
    loads=AppliedLoads(sigma_yield=250)
    text=calculation_report(sec,loads,sec.calculate_stresses(loads))
    assert "Unbounded (zero recovered stress" in text
    assert "not a general capacity guarantee" in text


def test_t23_independent_l_properties():
    sec=shape()
    info=inspect_section(sec)
    # shape is (0,0)->(2,0)->(2,1), t=1/100; local line integrals.
    area=Fraction(3,100)
    cx=Fraction(4,3)
    cy=Fraction(1,6)
    ix=Fraction(1,300)-area*cy*cy
    assert info.properties["A"]==pytest.approx(float(area),rel=1e-12,abs=0)
    assert info.properties["cx"]==pytest.approx(float(cx),rel=1e-12,abs=0)
    assert info.properties["Ix"]==pytest.approx(float(ix),rel=1e-12,abs=0)


def test_report_missing_image_fails_before_write(tmp_path):
    sec=shape()
    loads=AppliedLoads(N=1)
    with pytest.raises(ValueError):
        write_report(sec,loads,sec.calculate_stresses(loads),tmp_path/"r.md",plots={"x":tmp_path/"missing.png"})
    assert not (tmp_path/"r.md").exists()


def test_t28_publication_failure_cleans_temporary(tmp_path, monkeypatch):
    import sectalix.serialization as codec
    sec = shape()
    loads = AppliedLoads(N=1)
    def fail(source, destination):
        raise OSError("simulated publication failure")
    monkeypatch.setattr(codec.os, "link", fail)
    with pytest.raises(OSError, match="publication failure"):
        write_report(sec, loads, sec.calculate_stresses(loads), tmp_path / "r.md")
    assert list(tmp_path.iterdir()) == []
