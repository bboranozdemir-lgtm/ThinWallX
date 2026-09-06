"""T01--T15: independent wire records and adversarial codec tests."""
import copy
import json
import math
import struct
import sys
from dataclasses import replace

import pytest

from thinwallx import Node, Segment, Section, ClosedSection, MixedSection, AppliedLoads, GeometryError
from thinwallx.serialization import (
    DecodeLimits, UnitSystem, from_dict, from_json, to_dict, to_json, read_json, write_json,
)


def shape(cls=Section):
    points = [Node(0., 0., ("origin", 1)), Node(2., 0., "B"), Node(2., 1., "C"), Node(0., 1., "D")]
    edges = [(0, 1), (1, 2)] if cls is Section else [(0, 1), (1, 2), (2, 3), (3, 0)]
    return cls([Segment(points[a], points[b], .01, i) for i, (a, b) in enumerate(edges)])


@pytest.mark.parametrize("cls", [Section, ClosedSection, MixedSection])
def test_t01_t03_section_roundtrip(cls):
    source = shape(cls)
    wire = to_dict(source, units=UnitSystem("mm", "N"))
    output = from_json(to_json(source, units=UnitSystem("mm", "N")))
    assert type(output.value) is cls
    assert output.value.segments == source.segments
    assert to_dict(output) == wire
    wire["data"]["nodes"][0]["x"]["f64"] = "0x1.0000000000000p+0"
    assert source.segments[0].p1.x == 0.


def test_t03_pure_open_mixed():
    source = MixedSection(shape().segments)
    assert type(from_json(to_json(source)).value) is MixedSection


@pytest.mark.parametrize("x", [0., -0., 5e-324, -5e-324, sys.float_info.max, math.nextafter(1., 2.)])
def test_t04_t06_bits(x):
    result = from_json(to_json(AppliedLoads(N=x))).value
    assert struct.pack(">d", result.N) == struct.pack(">d", x)


def test_t04_all_loads():
    loads = AppliedLoads(1., -2., 3., -4., 5., -6., 7., -8., 250.)
    assert from_json(to_json(loads)).value == loads


@pytest.mark.parametrize("loads", [AppliedLoads(N=2.), AppliedLoads(sigma_yield=250.)])
def test_t05_results(loads):
    source = shape().calculate_stresses(loads)
    assert to_dict(from_json(to_json(source))) == to_dict(source)


@pytest.mark.parametrize("identifier", [None, "1", 1, True, 1.0, ("A", 2, -0.)])
def test_t07_ids(identifier):
    source = Section([Segment(Node(0, 0, identifier), Node(1, 1), .01, identifier)])
    out = from_json(to_json(source)).value
    assert type(out.segments[0].id) is type(identifier)
    assert out.segments[0].id == identifier


def test_t07_unsupported():
    with pytest.raises(TypeError):
        to_json(Section([Segment(Node(0, 0), Node(1, 0), .01, object())]))


@pytest.mark.parametrize("change", ["extra", "missing", "badref", "unused", "negative", "validate", "duplicate_ref"])
def test_t08_t11_invalid(change):
    d = to_dict(shape())
    if change == "extra": d["data"]["extra"] = 1
    if change == "missing": del d["data"]["node_tolerance"]
    if change == "badref": d["data"]["segments"][0]["p1"] = "absent"
    if change == "unused": d["data"]["nodes"].append({"ref": "extra", "id": {"type": "none"}, "x": {"f64":"0x0.0p+0"}, "y":{"f64":"0x0.0p+0"}})
    if change == "negative": d["data"]["segments"][0]["t"]["f64"] = "-0x1.0000000000000p+0"
    if change == "validate": d["data"]["validate"] = False
    if change == "duplicate_ref": d["data"]["nodes"][1]["ref"] = "n0"
    with pytest.raises(GeometryError if change == "negative" else ValueError):
        from_dict(d)


@pytest.mark.parametrize("token", ["nan", "inf", "0x1p+0", "0x1.0000000000001p-9999", "0x1.0p+9999"])
def test_t09_bad_float(token):
    d = to_dict(AppliedLoads())
    d["data"]["N"] = {"f64": token}
    with pytest.raises(ValueError):
        from_dict(d)


def test_t08_duplicate_key():
    with pytest.raises(ValueError, match="duplicate"):
        from_json('{"format":"thinwallx","format":"other"}')


def test_t09_literal_infinity():
    with pytest.raises(ValueError):
        from_json('{"bad":Infinity}')


def test_t12_units():
    d = from_json(to_json(shape(), units=UnitSystem("m", "kN")))
    assert from_json(to_json(d)).units == d.units
    with pytest.raises(ValueError):
        to_json(d, units=UnitSystem("mm", "N"))


def test_t13_tampered_peak():
    d = to_dict(shape().calculate_stresses(AppliedLoads(N=2)))
    d["data"]["segment_profiles"][0]["profile"]["peak_sigma_vm"] = {"f64": "0x0.0p+0"}
    with pytest.raises(ValueError):
        from_dict(d)


def test_t13_colliding_keys():
    d = to_dict(shape().calculate_stresses(AppliedLoads(N=2)))
    d["data"]["segment_profiles"][1] = copy.deepcopy(d["data"]["segment_profiles"][0])
    with pytest.raises(ValueError, match="colliding"):
        from_dict(d)


def test_t14_limits():
    with pytest.raises(ValueError, match="nesting"):
        from_json("[" * 80)
    with pytest.raises(ValueError, match="byte"):
        from_json(to_json(shape()), limits=DecodeLimits(max_bytes=10))
    with pytest.raises(ValueError, match="array"):
        from_dict(to_dict(shape()), limits=DecodeLimits(max_nodes=2))
    with pytest.raises(ValueError):
        DecodeLimits(max_depth=True)


def test_t15_atomic_and_deterministic(tmp_path):
    source = shape()
    a, b = tmp_path / "a.json", tmp_path / "b.json"
    write_json(source, a)
    write_json(source, b)
    assert a.read_bytes() == b.read_bytes()
    assert read_json(a).value.segments == source.segments
    with pytest.raises(FileExistsError):
        write_json(source, a)
    write_json(source, a, overwrite=True)
    assert not list(tmp_path.glob(".thinwallx-*"))


def test_t06_handwritten_wire():
    # Independent IEEE-754 minimum-subnormal token, not produced by the codec.
    fields = ("N", "Vx", "Vy", "Mx", "My", "Tsv", "B", "M_omega")
    d = {"format":"thinwallx", "schema_version":"0.8", "number_encoding":"float64-hex",
         "kind":"applied_loads", "units":{"length":"mm","force":"N"},
         "data":{k:{"f64":"0x0.0p+0"} for k in fields}}
    d["data"].update(N={"f64":"0x0.0000000000001p-1022"}, sigma_yield=None)
    assert struct.pack(">d", from_dict(d).value.N) == bytes.fromhex("0000000000000001")


def test_t01_near_coincident_raw_endpoints():
    sec = Section([Segment(Node(-1., 0.), Node(-0., 0., "a"), .01, "first"),
                   Segment(Node(1e-10, 0., "b"), Node(0., 1.), .01, "second")])
    out = from_json(to_json(sec)).value
    for before, after in zip(sec.segments, out.segments):
        for a, b in zip((before.p1, before.p2), (after.p1, after.p2)):
            assert struct.pack(">dd", a.x, a.y) == struct.pack(">dd", b.x, b.y)
            assert a.id == b.id


def test_t11_unvalidated_cycle_rejected():
    from thinwallx import TopologyError
    sec = Section(shape(ClosedSection).segments, validate=False)
    with pytest.raises(TopologyError):
        to_dict(sec)


def test_t14_surrogate_and_wrong_units():
    with pytest.raises(ValueError):
        from_json('"\ud800"')
    with pytest.raises(ValueError):
        UnitSystem(length=[])


def test_t15_failed_atomic_write(tmp_path, monkeypatch):
    from thinwallx.serialization import _atomic_write
    def fail(path):
        path.write_text("partial", encoding="utf-8")
        raise OSError("simulated write failure")
    with pytest.raises(OSError):
        _atomic_write(tmp_path/"target.json", fail, False)
    assert not list(tmp_path.iterdir())
