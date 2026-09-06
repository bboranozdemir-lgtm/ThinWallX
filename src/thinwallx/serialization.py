"""Strict, lossless v0.8 document codec. See ACTIVE_PHASE_v0.8.md §2.

Float hex strings preserve binary64 bits; typed IDs never become JSON keys.
No geometry repair or executable object deserialization is performed.
"""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, fields
import json
import math
import os
from pathlib import Path
import re
import tempfile
from typing import Any, Callable

from thinwallx.closed_section import ClosedSection
from thinwallx.mixed_section import MixedSection
from thinwallx.primitives import Node, Segment
from thinwallx.section import Section
from thinwallx.stress import AppliedLoads, SegmentStressProfile, StressRecoveryResult

Serializable = Section | ClosedSection | MixedSection | AppliedLoads | StressRecoveryResult
_LENGTH_UNITS = frozenset(("m", "mm", "cm", "in", "ft", "unspecified"))
_FORCE_UNITS = frozenset(("N", "kN", "lbf", "unspecified"))
_LOADS = ("N", "Vx", "Vy", "Mx", "My", "Tsv", "B", "M_omega", "sigma_yield")
_RESULT = ("max_sigma_vm", "peak_s", "resultant_N", "resultant_Mx", "resultant_My", "resultant_B")
_PROFILE = ("length", "thickness", "tau_sv_surface", "peak_sigma_vm", "s_peak")


@dataclass(frozen=True)
class UnitSystem:
    length: str = "unspecified"
    force: str = "unspecified"

    def __post_init__(self) -> None:
        if not isinstance(self.length, str) or not isinstance(self.force, str) or self.length not in _LENGTH_UNITS or self.force not in _FORCE_UNITS:
            raise ValueError("units: unsupported length or force unit")


@dataclass(frozen=True)
class DecodeLimits:
    max_bytes: int = 16777216
    max_depth: int = 64
    max_nodes: int = 10000
    max_segments: int = 5000
    max_string_chars: int = 4096
    max_profiles: int = 5000

    def __post_init__(self) -> None:
        for f in fields(self):
            v = getattr(self, f.name)
            if type(v) is not int or v <= 0:
                raise ValueError(f"limits.{f.name}: expected positive integer")


@dataclass(frozen=True)
class DecodedDocument:
    value: Serializable
    units: UnitSystem


def _record(data: object, keys: tuple[str, ...], path: str) -> Mapping[str, Any]:
    if not isinstance(data, Mapping) or set(data) != set(keys):
        raise ValueError(f"{path}: expected exactly fields {keys}")
    return data


def _string(value: object, path: str) -> str:
    if not isinstance(value, str):
        raise ValueError(f"{path}: expected string")
    try:
        value.encode("utf-8")
    except UnicodeEncodeError as exc:
        raise ValueError(f"{path}: invalid Unicode surrogate") from exc
    return value


def _f64(value: float) -> dict[str, str]:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError("physical scalar must be a real binary64-compatible number")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError("physical scalar must be finite")
    return {"f64": value.hex()}


def _unfloat(data: object, path: str) -> float:
    rec = _record(data, ("f64",), path)
    token = _string(rec["f64"], path + ".f64")
    try:
        value = float.fromhex(token)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{path}: invalid canonical F64") from exc
    if not math.isfinite(value) or value.hex() != token:
        raise ValueError(f"{path}: expected canonical finite F64")
    return value


def _id(value: Any, depth: int = 0) -> dict[str, Any]:
    if depth > 64:
        raise ValueError("ID nesting exceeds 64")
    if value is None:
        return {"type": "none"}
    if type(value) is str:
        return {"type": "str", "value": _string(value, "id")}
    if type(value) is bool:
        return {"type": "bool", "value": value}
    if type(value) is int:
        return {"type": "int", "value": str(value)}
    if type(value) is float:
        if not math.isfinite(value):
            raise TypeError("Portable float IDs must be finite")
        return {"type": "float", "value": _f64(value)}
    if type(value) is tuple:
        return {"type": "tuple", "items": [_id(v, depth + 1) for v in value]}
    raise TypeError(f"Unsupported portable ID type: {type(value).__name__}")


def _unid(data: object, path: str) -> Any:
    if not isinstance(data, Mapping) or "type" not in data:
        raise ValueError(f"{path}: expected typed ID")
    tag = data["type"]
    keys = ("type",) if tag == "none" else ("type", "items") if tag == "tuple" else ("type", "value")
    rec = _record(data, keys, path)
    if tag == "none":
        return None
    if tag == "str":
        return _string(rec["value"], path)
    if tag == "bool" and type(rec["value"]) is bool:
        return rec["value"]
    if tag == "int":
        token = _string(rec["value"], path)
        if re.fullmatch(r"0|-?[1-9][0-9]*", token):
            return int(token)
    if tag == "float":
        return _unfloat(rec["value"], path)
    if tag == "tuple" and isinstance(rec["items"], list):
        return tuple(_unid(v, f"{path}.items[{i}]") for i, v in enumerate(rec["items"]))
    raise ValueError(f"{path}: invalid portable ID")


def _walk(value: object, limits: DecodeLimits, path: str = "$", depth: int = 0) -> None:
    if depth > limits.max_depth:
        raise ValueError(f"{path}: nesting limit exceeded")
    if isinstance(value, str):
        _string(value, path)
        if len(value) > limits.max_string_chars:
            raise ValueError(f"{path}: string limit exceeded")
    elif isinstance(value, Mapping):
        for k, v in value.items():
            _string(k, path)
            _walk(k, limits, path, depth + 1)
            _walk(v, limits, f"{path}.{k}", depth + 1)
    elif isinstance(value, list):
        # Every format array has at most the largest declared item budget.
        if len(value) > max(limits.max_nodes, limits.max_segments, limits.max_profiles):
            raise ValueError(f"{path}: array limit exceeded")
        for i, v in enumerate(value):
            _walk(v, limits, f"{path}[{i}]", depth + 1)
    elif value is not None and type(value) not in (bool, int, float):
        raise ValueError(f"{path}: not a JSON value")


def _array(value: object, path: str, minimum: int, maximum: int) -> list[Any]:
    if not isinstance(value, list) or not minimum <= len(value) <= maximum:
        raise ValueError(f"{path}: invalid array size")
    return value


def _decode_section(data: object, limits: DecodeLimits) -> Section:
    d = _record(data, ("class", "node_tolerance", "safety_factor", "nodes", "segments"), "data")
    classes = {"Section": Section, "ClosedSection": ClosedSection, "MixedSection": MixedSection}
    name = _string(d["class"], "data.class")
    if name not in classes:
        raise ValueError("data.class: unsupported class")
    tol = _unfloat(d["node_tolerance"], "data.node_tolerance")
    if tol < 0:
        raise ValueError("data.node_tolerance: must be nonnegative")
    nodes: dict[str, Node] = {}
    for i, raw in enumerate(_array(d["nodes"], "data.nodes", 2, limits.max_nodes)):
        path = f"data.nodes[{i}]"
        n = _record(raw, ("ref", "id", "x", "y"), path)
        ref = _string(n["ref"], path + ".ref")
        if not ref or ref in nodes:
            raise ValueError(f"{path}.ref: empty or duplicate reference")
        nodes[ref] = Node(_unfloat(n["x"], path + ".x"), _unfloat(n["y"], path + ".y"), _unid(n["id"], path + ".id"))
    segments: list[Segment] = []
    used: set[str] = set()
    for i, raw in enumerate(_array(d["segments"], "data.segments", 1, limits.max_segments)):
        path = f"data.segments[{i}]"
        s = _record(raw, ("id", "p1", "p2", "t"), path)
        a, b = (_string(s[k], path + "." + k) for k in ("p1", "p2"))
        if a not in nodes or b not in nodes:
            raise ValueError(f"{path}: undefined endpoint")
        used.update((a, b))
        segments.append(Segment(nodes[a], nodes[b], _unfloat(s["t"], path + ".t"), _unid(s["id"], path + ".id")))
    if used != set(nodes):
        raise ValueError("data.nodes: unused reference")
    if name == "Section":
        if d["safety_factor"] is not None:
            raise ValueError("data.safety_factor: Section requires null")
        return Section(segments, validate=True, node_tolerance=tol)
    sf = _unfloat(d["safety_factor"], "data.safety_factor")
    return classes[name](segments, validate=True, node_tolerance=tol, safety_factor=sf)


def _decode_result(data: object, limits: DecodeLimits) -> StressRecoveryResult:
    d = _record(data, (*_RESULT, "segment_profiles", "peak_segment_id", "peak_location_xy", "load_factor"), "data")
    profiles: dict[Any, SegmentStressProfile] = {}
    for i, raw in enumerate(_array(d["segment_profiles"], "data.segment_profiles", 1, limits.max_profiles)):
        path = f"data.segment_profiles[{i}]"
        entry = _record(raw, ("key", "profile"), path)
        key = _unid(entry["key"], path + ".key")
        p = _record(entry["profile"], (*_PROFILE, "segment_id", "sigma_zz_coeffs", "tau_membrane_coeffs"), path + ".profile")
        sid = _unid(p["segment_id"], path + ".profile.segment_id")
        if _id(key) != _id(sid) or key in profiles:
            raise ValueError(f"{path}: mismatched or colliding profile key")
        vals = {k: _unfloat(p[k], path + "." + k) for k in _PROFILE}
        sigma = tuple(_unfloat(v, path + ".sigma_zz_coeffs") for v in _array(p["sigma_zz_coeffs"], path, 2, 2))
        tau = tuple(_unfloat(v, path + ".tau_membrane_coeffs") for v in _array(p["tau_membrane_coeffs"], path, 3, 3))
        profile = SegmentStressProfile(sid, vals["length"], vals["thickness"], sigma, tau, vals["tau_sv_surface"])
        for k in ("s_peak", "peak_sigma_vm"):
            if float(getattr(profile, k)).hex() != vals[k].hex():
                raise ValueError(f"{path}.{k}: archived peak differs from computed peak")
        profiles[key] = profile
    values = {k: _unfloat(d[k], "data." + k) for k in _RESULT}
    peak_id = _unid(d["peak_segment_id"], "data.peak_segment_id")
    if peak_id not in profiles or _id(peak_id) != _id(profiles[peak_id].segment_id):
        raise ValueError("data.peak_segment_id: undefined or mismatched ID")
    maximum = max(p.peak_sigma_vm for p in profiles.values())
    if values["max_sigma_vm"] != maximum or profiles[peak_id].peak_sigma_vm != maximum:
        raise ValueError("data.max_sigma_vm: inconsistent maximum")
    if values["peak_s"].hex() != float(profiles[peak_id].s_peak).hex():
        raise ValueError("data.peak_s: inconsistent location")
    loc = tuple(_unfloat(v, "data.peak_location_xy") for v in _array(d["peak_location_xy"], "data.peak_location_xy", 2, 2))
    factor: float | None = None
    if d["load_factor"] == {"special": "positive_infinity"}:
        if maximum != 0:
            raise ValueError("data.load_factor: infinite factor requires zero stress")
        factor = math.inf
    elif d["load_factor"] is not None:
        factor = _unfloat(d["load_factor"], "data.load_factor")
        if factor <= 0:
            raise ValueError("data.load_factor: must be positive")
    return StressRecoveryResult(profiles, peak_segment_id=peak_id, peak_location_xy=loc, load_factor=factor, **values)


def from_dict(data: Mapping[str, object], *, limits: DecodeLimits | None = None) -> DecodedDocument:
    limits = limits or DecodeLimits()
    _walk(data, limits)
    d = _record(data, ("format", "schema_version", "number_encoding", "kind", "units", "data"), "$")
    if (d["format"], d["schema_version"], d["number_encoding"]) != ("thinwallx", "0.8", "float64-hex"):
        raise ValueError("$: unsupported format/version/encoding")
    units = UnitSystem(**_record(d["units"], ("length", "force"), "units"))
    if d["kind"] == "section":
        value = _decode_section(d["data"], limits)
    elif d["kind"] == "applied_loads":
        rec = _record(d["data"], _LOADS, "data")
        values = {k: None if k == "sigma_yield" and rec[k] is None else _unfloat(rec[k], "data." + k) for k in _LOADS}
        value = AppliedLoads(**values)
    elif d["kind"] == "stress_result":
        value = _decode_result(d["data"], limits)
    else:
        raise ValueError("kind: unsupported kind")
    return DecodedDocument(value, units)


def to_dict(value: Serializable | DecodedDocument, *, units: UnitSystem | None = None) -> dict[str, object]:
    if isinstance(value, DecodedDocument):
        if units is not None and units != value.units:
            raise ValueError("units: relabeling a document is forbidden")
        units, value = value.units, value.value
    units = units or UnitSystem()
    if not isinstance(units, UnitSystem):
        raise TypeError("units: expected UnitSystem")
    data: dict[str, Any]
    if type(value) in (Section, ClosedSection, MixedSection):
        # Validate raw geometry without recomputing section properties/caches.
        if type(value) is Section:
            from thinwallx.validation import validate_section_geometry_and_topology
            validate_section_geometry_and_topology(value.segments, node_tolerance=value._node_tolerance)
        elif type(value) is ClosedSection:
            from thinwallx.cells import extract_cell_topology
            extract_cell_topology(value.segments, node_tolerance=value._node_tolerance, safety_factor=value._safety_factor)
        else:
            from thinwallx.mixed_topology import extract_mixed_topology
            extract_mixed_topology(value.segments, node_tolerance=value._node_tolerance, safety_factor=value._safety_factor)
        nodes, segments = [], []
        for i, s in enumerate(value.segments):
            for j, n in enumerate((s.p1, s.p2)):
                nodes.append({"ref": f"n{2*i+j}", "id": _id(n.id), "x": _f64(n.x), "y": _f64(n.y)})
            segments.append({"id": _id(s.id), "p1": f"n{2*i}", "p2": f"n{2*i+1}", "t": _f64(s.t)})
        data = {"class": type(value).__name__, "node_tolerance": _f64(value._node_tolerance),
                "safety_factor": None if type(value) is Section else _f64(value._safety_factor),
                "nodes": nodes, "segments": segments}
        kind = "section"
    elif type(value) is AppliedLoads:
        data = {k: None if getattr(value, k) is None else _f64(getattr(value, k)) for k in _LOADS}
        kind = "applied_loads"
    elif type(value) is StressRecoveryResult:
        entries = []
        for key, p in value.segment_profiles.items():
            rec = {k: _f64(getattr(p, k)) for k in _PROFILE}
            rec.update(segment_id=_id(p.segment_id), sigma_zz_coeffs=[_f64(v) for v in p.sigma_zz_coeffs],
                       tau_membrane_coeffs=[_f64(v) for v in p.tau_membrane_coeffs])
            entries.append({"key": _id(key), "profile": rec})
        data = {k: _f64(getattr(value, k)) for k in _RESULT}
        data.update(segment_profiles=entries, peak_segment_id=_id(value.peak_segment_id),
                    peak_location_xy=[_f64(v) for v in value.peak_location_xy],
                    load_factor=None if value.load_factor is None else
                    {"special": "positive_infinity"} if value.load_factor == math.inf else _f64(value.load_factor))
        _decode_result(data, DecodeLimits(max_profiles=max(5000, len(entries))))
        kind = "stress_result"
    else:
        raise TypeError(f"Unsupported serializable type: {type(value).__name__}")
    return {"format": "thinwallx", "schema_version": "0.8", "number_encoding": "float64-hex",
            "kind": kind, "units": {"length": units.length, "force": units.force}, "data": data}


def to_json(value: Serializable | DecodedDocument, *, units: UnitSystem | None = None, indent: int | None = None) -> str:
    if indent is not None and (type(indent) is not int or indent < 0):
        raise ValueError("indent must be a nonnegative integer or None")
    return json.dumps(to_dict(value, units=units), ensure_ascii=False, allow_nan=False,
                      sort_keys=True, separators=(",", ":") if indent is None else None, indent=indent)


def _pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for k, v in pairs:
        if k in result:
            raise ValueError(f"duplicate JSON key: {k}")
        result[k] = v
    return result


def _constant(token: str) -> None:
    raise ValueError(f"forbidden JSON constant: {token}")


def from_json(text: str, *, limits: DecodeLimits | None = None) -> DecodedDocument:
    limits = limits or DecodeLimits()
    _string(text, "$")
    if len(text.encode("utf-8")) > limits.max_bytes:
        raise ValueError("$: byte limit exceeded")
    depth, quoted, escaped = 0, False, False
    for c in text:
        if quoted:
            if escaped:
                escaped = False
            elif c == "\\":
                escaped = True
            elif c == '"':
                quoted = False
        elif c == '"':
            quoted = True
        elif c in "[{":
            depth += 1
            if depth > limits.max_depth:
                raise ValueError("$: nesting limit exceeded")
        elif c in "]}":
            depth -= 1
    data = json.loads(text, object_pairs_hook=_pairs, parse_constant=_constant)
    return from_dict(data, limits=limits)


def _atomic_write(path: str | Path, writer: Callable[[Path], None], overwrite: bool) -> Path:
    """Publish only a complete file; exclusive hard-link avoids no-clobber races."""
    target = Path(path).absolute()
    if not overwrite and target.exists():
        raise FileExistsError(target)
    fd, name = tempfile.mkstemp(prefix=".thinwallx-", suffix=target.suffix, dir=target.parent)
    os.close(fd)
    temporary = Path(name)
    try:
        writer(temporary)
        if overwrite:
            os.replace(temporary, target)
        else:
            os.link(temporary, target)
    finally:
        temporary.unlink(missing_ok=True)
    return target


def write_json(value: Serializable | DecodedDocument, path: str | Path, *, units: UnitSystem | None = None, overwrite: bool = False) -> Path:
    text = to_json(value, units=units) + "\n"
    return _atomic_write(path, lambda p: p.write_text(text, encoding="utf-8", newline="\n"), overwrite)


def read_json(path: str | Path, *, limits: DecodeLimits | None = None) -> DecodedDocument:
    limits = limits or DecodeLimits()
    with Path(path).open("rb") as stream:
        raw = stream.read(limits.max_bytes + 1)
    if len(raw) > limits.max_bytes:
        raise ValueError("$: byte limit exceeded")
    return from_json(raw.decode("utf-8"), limits=limits)
