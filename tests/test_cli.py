"""v1 T01--T20: real module processes and application error categories."""
from __future__ import annotations
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from sectalix import AppliedLoads
from sectalix.cli import main
from sectalix.serialization import from_json, to_json, write_json, UnitSystem
from tests.test_serialization import shape

ROOT = Path(__file__).resolve().parents[1]
DXF = ROOT/"tests"/"fixtures"/"dxf"


def run(*args: str, stdin: str | None = None) -> subprocess.CompletedProcess[str]:
    env = dict(os.environ, PYTHONPATH=str(ROOT/"src"), PYTHONIOENCODING="utf-8")
    return subprocess.run([sys.executable, "-W", "error", "-m", "sectalix", *map(str,args)],
                          input=stdin, capture_output=True, text=True, encoding="utf-8", env=env)


@pytest.mark.parametrize("flag",["--help","--version"])
def test_t01(flag):
    p=run(flag)
    assert p.returncode==0 and not p.stderr
    assert ("1.0.1" if flag=="--version" else "convert-dxf") in p.stdout


@pytest.mark.parametrize("args",[[],["bad"],["inspect"],["inspect","x.json","--bad"]])
def test_t02_usage(args):
    p=run(*args)
    assert p.returncode==1 and p.stderr and not p.stdout


def test_t03_missing():
    p=run("inspect","absent.json")
    assert p.returncode==2 and not p.stdout


@pytest.mark.parametrize("name,kind",[("open_l_r12","open"),("rectangle_r12","closed"),("barbell_ac1015","mixed")])
def test_t04_t06_inspection(name,kind):
    p=run("inspect",DXF/(name+".dxf"),"--json")
    assert p.returncode==0, p.stderr
    record=json.loads(p.stdout)
    assert record["topology"]==kind
    assert set(("A","Ix","Iy","J","Cw","sx","sy")) <= set(record["properties"])
    restored=from_json(json.dumps(record["section"]))
    assert restored.units.length=="mm"


def test_t07_geometry(tmp_path):
    d=json.loads(to_json(shape()))
    d["data"]["segments"][0]["t"]["f64"]="-0x1.0000000000000p+0"
    p=run("inspect","-","--json",stdin=json.dumps(d))
    assert p.returncode==3 and not p.stdout and "Traceback" not in p.stderr


def test_t08_t12_load_file_and_override(tmp_path):
    sec=tmp_path/"sec.json"
    loads=tmp_path/"loads.json"
    write_json(shape(),sec,units=UnitSystem("mm","N"))
    write_json(AppliedLoads(N=2),loads,units=UnitSystem("mm","N"))
    p=run("analyze",sec,"--loads",loads,"--N","3","--stdout")
    assert p.returncode==0,p.stderr
    result=from_json(p.stdout).value
    assert result.resultant_N==pytest.approx(3,rel=1e-12,abs=0)


def test_t09_negative_scientific():
    p=run("analyze",DXF/"open_l_r12.dxf","--N","1","--Mx","-1e-2")
    assert p.returncode==0,p.stderr
    assert from_json(p.stdout).value.resultant_Mx==pytest.approx(-.01,rel=1e-12,abs=0)


def test_t10_t12_artifacts(tmp_path):
    p=run("analyze",DXF/"rectangle_r12.dxf","--N","1","--plots-dir",tmp_path/"plots",
          "--report",tmp_path/"report.md","--output",tmp_path/"results.json")
    assert p.returncode==0,p.stderr
    assert not p.stdout
    assert {p.name for p in (tmp_path/"plots").iterdir()}=={"geometry.png","shear_flow.png","stress_vm.png"}
    assert (tmp_path/"report.md").read_text(encoding="utf-8").startswith("# Sectalix")
    assert from_json((tmp_path/"results.json").read_text()).value.max_sigma_vm>0


def test_t13_singular():
    p=run("analyze",DXF/"open_l_r12.dxf","--B","1")
    assert p.returncode==5 and not p.stdout


@pytest.mark.parametrize("value",["1e309","1e-999"])
def test_t14_cli_range(value):
    p=run("analyze",DXF/"open_l_r12.dxf","--N",value)
    assert p.returncode==6 and not p.stdout


def test_t14_mechanical_overflow():
    p=run("analyze",DXF/"open_l_r12.dxf","--N","1e308")
    assert p.returncode==6 and not p.stdout


def test_t15_conversion(tmp_path):
    out=tmp_path/"converted.json"
    p=run("convert-dxf",DXF/"rectangle_r12.dxf","--output",out)
    assert p.returncode==0,p.stderr
    assert len(from_json(out.read_text()).value.segments)==4


@pytest.mark.parametrize("fmt",["svg","png"])
def test_t16_plot(tmp_path,fmt):
    p=run("plot",DXF/"open_l_r12.dxf","--show-ids","--show-axes","--output",tmp_path/("geom."+fmt))
    assert p.returncode==0,p.stderr


def test_t17_t18_stdin():
    p=run("inspect","-","--json",stdin=to_json(shape()))
    assert p.returncode==0 and not p.stderr
    assert json.loads(p.stdout)["node_count"]==3
    p=run("analyze","-","--N","1","--stdout",stdin=to_json(shape()))
    assert p.returncode==0 and not p.stderr
    assert from_json(p.stdout).value.resultant_N==pytest.approx(1,rel=1e-12,abs=0)


def test_t19_missing_matplotlib(tmp_path):
    code="""import sys, importlib.abc
class Block(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'matplotlib' or fullname.startswith('matplotlib.'):
            raise ModuleNotFoundError('matplotlib intentionally unavailable')
sys.meta_path.insert(0, Block())
from sectalix.cli import main
raise SystemExit(main(sys.argv[1:]))
"""
    env=dict(os.environ,PYTHONPATH=str(ROOT/"src"))
    for args,expected in [
        (["inspect",str(DXF/"open_l_r12.dxf")],0),
        (["analyze",str(DXF/"open_l_r12.dxf"),"--N","1"],0),
        (["plot",str(DXF/"open_l_r12.dxf"),"--output",str(tmp_path/"x.png")],2),
    ]:
        p=subprocess.run([sys.executable,"-W","error","-c",code,*args],env=env,capture_output=True,text=True)
        assert p.returncode==expected,p.stderr


def test_t20_no_clobber(tmp_path):
    target=tmp_path/"out.json"
    target.write_text("keep",encoding="utf-8")
    p=run("convert-dxf",DXF/"open_l_r12.dxf","--output",target)
    assert p.returncode==2
    assert target.read_text()=="keep"


def test_topology_code():
    raw=to_json(shape())
    data=json.loads(raw)
    data["data"]["nodes"][2]["x"]["f64"]="0x1.0000000000000p+4"
    data["data"]["nodes"][3]["x"]["f64"]="0x1.0000000000000p+4"
    p=run("inspect","-",stdin=json.dumps(data))
    assert p.returncode==4,p.stderr


def test_unit_mismatch(tmp_path):
    sec,loads=tmp_path/"s.json",tmp_path/"l.json"
    write_json(shape(),sec,units=UnitSystem("mm","N"))
    write_json(AppliedLoads(N=1),loads,units=UnitSystem("m","N"))
    p=run("analyze",sec,"--loads",loads)
    assert p.returncode==2 and "units differ" in p.stderr
