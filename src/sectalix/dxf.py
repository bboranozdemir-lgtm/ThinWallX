"""Restricted ASCII DXF centerline import for AC1009/AC1015 files.

Entity interpretation follows Autodesk DXF LINE/POLYLINE/LWPOLYLINE references.
Exact rational projection predicates avoid overflow in tolerance squares and
orientation products. Snapping is an explicit, reported importer operation.
"""
from __future__ import annotations

from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass, field
from fractions import Fraction
import codecs
import math
from pathlib import Path
import re
from types import MappingProxyType
from typing import Any

from sectalix import Node, Segment, Section, ClosedSection, MixedSection
from sectalix.exceptions import GeometryError, TopologyError
from sectalix.mixed_topology import extract_mixed_topology
from sectalix.serialization import _unfloat, _record
from sectalix.validation import cluster_nodes

Point = tuple[float, float]
Pair = tuple[int, str, int]
_UNITS = {"m": Fraction(1), "mm": Fraction(1, 1000), "cm": Fraction(1, 100),
          "in": Fraction(127, 5000), "ft": Fraction(381, 1250), "unspecified": Fraction(1)}
_DECIMAL = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z")
_META = {5, 6, 8, 48, 60, 62, 67, 100, 102, 330, 370, 410, 420, 430, 440, 999}


def _finite(value: float, name: str, *, positive: bool = False) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{name}: expected finite real number")
    if positive and value <= 0:
        raise ValueError(f"{name}: expected positive value")


@dataclass(frozen=True)
class ThicknessMap:
    by_handle: Mapping[str, float] = field(default_factory=dict)
    by_layer: Mapping[str, float] = field(default_factory=dict)
    use_layer_name: bool = True
    default: float | None = None

    def __post_init__(self) -> None:
        if type(self.use_layer_name) is not bool:
            raise ValueError("use_layer_name: expected bool")
        handle: dict[str, float] = {}
        for key, value in self.by_handle.items():
            if not isinstance(key, str) or not re.fullmatch(r"[0-9a-fA-F]{1,16}", key):
                raise ValueError("by_handle: invalid handle")
            if key.upper() in handle:
                raise ValueError("by_handle: duplicate canonical handle")
            _finite(value, "by_handle", positive=True)
            handle[key.upper()] = value
        layer: dict[str, float] = {}
        for key, value in self.by_layer.items():
            if not isinstance(key, str) or not key:
                raise ValueError("by_layer: invalid layer")
            _finite(value, "by_layer", positive=True)
            layer[key] = value
        if self.default is not None:
            _finite(self.default, "default thickness", positive=True)
        object.__setattr__(self, "by_handle", MappingProxyType(handle))
        object.__setattr__(self, "by_layer", MappingProxyType(layer))

    @classmethod
    def from_dict(cls, data: Mapping[str, object]) -> ThicknessMap:
        d = _record(data, ("by_handle", "by_layer", "use_layer_name", "default"), "thickness")
        maps: dict[str, dict[str, float]] = {}
        for name in ("by_handle", "by_layer"):
            if not isinstance(d[name], Mapping):
                raise ValueError(f"thickness.{name}: expected mapping")
            maps[name] = {k: _unfloat(v, f"thickness.{name}.{k}") for k, v in d[name].items()}
        return cls(**maps, use_layer_name=d["use_layer_name"],
                   default=None if d["default"] is None else _unfloat(d["default"], "thickness.default"))


@dataclass(frozen=True)
class DxfImportOptions:
    source_length_unit: str
    target_length_unit: str
    thickness: ThicknessMap
    section_kind: str = "auto"
    layers: frozenset[str] | None = None
    node_tolerance: float = 1e-9
    plane_z: float = 0.0
    plane_tolerance: float = 0.0
    junction_policy: str = "split_endpoint"
    safety_factor: float = 1e4
    max_bytes: int = 16777216
    max_entities: int = 10000
    max_output_segments: int = 5000
    max_pair_checks: int = 2000000

    def __post_init__(self) -> None:
        if self.source_length_unit not in _UNITS or self.target_length_unit not in _UNITS:
            raise ValueError("Unsupported length unit")
        if (self.source_length_unit == "unspecified") != (self.target_length_unit == "unspecified"):
            raise ValueError("Cannot convert specified and unspecified units")
        if not isinstance(self.thickness, ThicknessMap):
            raise TypeError("thickness must be ThicknessMap")
        if self.section_kind not in ("auto", "open", "closed", "mixed"):
            raise ValueError("Invalid section_kind")
        if self.junction_policy not in ("split_endpoint", "reject"):
            raise ValueError("Invalid junction_policy")
        if self.layers is not None and (not isinstance(self.layers, frozenset) or
                                        any(not isinstance(x, str) or not x for x in self.layers)):
            raise ValueError("layers must be a frozenset of nonempty strings")
        for name in ("node_tolerance", "plane_z", "plane_tolerance", "safety_factor"):
            _finite(getattr(self, name), name)
        if self.node_tolerance < 0 or self.plane_tolerance < 0 or self.safety_factor <= 0:
            raise ValueError("Invalid tolerance or safety_factor")
        for name in ("max_bytes", "max_entities", "max_output_segments", "max_pair_checks"):
            if type(getattr(self, name)) is not int or getattr(self, name) <= 0:
                raise ValueError(f"{name}: expected positive integer")


@dataclass(frozen=True)
class DxfImportReport:
    version: str
    source_length_unit: str
    target_length_unit: str
    conversion_ratio: str
    selected_entities: Mapping[str, int]
    excluded_entities: Mapping[str, int]
    provenance: tuple[Mapping[str, Any], ...]
    changes: tuple[Mapping[str, Any], ...]
    max_plane_displacement: float
    max_endpoint_displacement: float
    removed_closing_vertices: int
    node_count: int
    edge_count: int
    component_count: int
    cycle_rank: int
    section_class: str
    thickness_sources: tuple[Mapping[str, Any], ...]


@dataclass(frozen=True)
class DxfImportResult:
    section: Section | ClosedSection | MixedSection
    report: DxfImportReport


@dataclass(frozen=True)
class _Entity:
    kind: str
    pairs: tuple[Pair, ...]
    line: int

    def one(self, code: int, default: str | None = None) -> str:
        matches = [v for c, v, _ in self.pairs if c == code]
        if len(matches) > 1:
            raise ValueError(f"line {self.line} {self.kind}: duplicate code {code}")
        if not matches and default is None:
            raise ValueError(f"line {self.line} {self.kind}: missing code {code}")
        return matches[0] if matches else default  # type: ignore[return-value]


def _rational(token: str) -> Fraction:
    if len(token) > 4096 or not _DECIMAL.fullmatch(token):
        raise ValueError("Invalid DXF decimal token")
    # Bound exponent construction before allocating enormous Python integers.
    exponent = re.search(r"[eE]([+-]?[0-9]+)$", token)
    if exponent and (len(exponent[1].lstrip("+-")) > 5 or abs(int(exponent[1])) > 10000):
        raise ValueError("DXF decimal exponent exceeds resource limit")
    return Fraction(token)


def _float(value: Fraction) -> float:
    out = float(value)
    if not math.isfinite(out):
        raise OverflowError("DXF converted value exceeds binary64")
    if value and out == 0:
        raise FloatingPointError("DXF nonzero value underflows binary64")
    return out


def _integer(token: str) -> int:
    if not re.fullmatch(r"[+-]?[0-9]{1,10}", token):
        raise ValueError("Invalid DXF integer")
    return int(token)


def _parse(text: str, options: DxfImportOptions) -> tuple[dict[str, str], list[_Entity]]:
    if not isinstance(text, str):
        raise TypeError("DXF text must be str")
    if len(text.encode("utf-8")) > options.max_bytes:
        raise ValueError("DXF byte limit exceeded")
    if text.startswith("AutoCAD Binary DXF") or "\x00" in text:
        raise ValueError("Binary DXF is unsupported")
    lines = text.splitlines()
    # Only blank lines AFTER EOF may be discarded (empty values are significant).
    while lines and not lines[-1].strip():
        lines.pop()
    if len(lines) % 2:
        raise ValueError("DXF has an unpaired group code")
    pairs: list[Pair] = []
    for i in range(0, len(lines), 2):
        code = _integer(lines[i].strip())
        if not 0 <= code <= 1071:
            raise ValueError(f"line {i+1}: invalid group code")
        pairs.append((code, lines[i+1].strip(), i+1))
    sections: dict[str, list[Pair]] = {}
    i = 0
    while i < len(pairs) and pairs[i][:2] != (0, "EOF"):
        if pairs[i][:2] != (0, "SECTION") or i+1 >= len(pairs) or pairs[i+1][0] != 2:
            raise ValueError(f"line {pairs[i][2]}: expected SECTION/name")
        name = pairs[i+1][1]
        if name not in ("HEADER", "ENTITIES", "TABLES", "BLOCKS", "CLASSES", "OBJECTS") or name in sections:
            raise ValueError("Unsupported or duplicate DXF section")
        i += 2
        body: list[Pair] = []
        while i < len(pairs) and pairs[i][:2] != (0, "ENDSEC"):
            if pairs[i][:2] in ((0, "SECTION"), (0, "EOF")):
                raise ValueError("Unterminated DXF section")
            body.append(pairs[i])
            i += 1
        if i == len(pairs):
            raise ValueError("Missing ENDSEC")
        sections[name] = body
        i += 1
    if i != len(pairs)-1 or not pairs or pairs[i][:2] != (0, "EOF"):
        raise ValueError("Missing EOF or trailing data")
    if "HEADER" not in sections or "ENTITIES" not in sections:
        raise ValueError("HEADER and ENTITIES are required")
    header: dict[str, str] = {}
    body = sections["HEADER"]
    j = 0
    seen: set[str] = set()
    while j < len(body):
        if body[j][0] == 999:
            j += 1
            continue
        if body[j][0] != 9:
            raise ValueError("Expected HEADER variable")
        name = body[j][1]
        if name in seen:
            raise ValueError("Duplicate HEADER variable")
        seen.add(name)
        j += 1
        vals: list[Pair] = []
        while j < len(body) and body[j][0] != 9:
            if body[j][0] != 999:
                vals.append(body[j])
            j += 1
        expected = {"$ACADVER": 1, "$INSUNITS": 70, "$DWGCODEPAGE": 3}
        if name in expected:
            if len(vals) != 1 or vals[0][0] != expected[name]:
                raise ValueError(f"Invalid HEADER {name}")
            header[name] = vals[0][1]
    if header.get("$ACADVER") not in ("AC1009", "AC1015"):
        raise ValueError("Unsupported or missing ACADVER")
    entities: list[_Entity] = []
    body = sections["ENTITIES"]
    j = 0
    handles: set[str] = set()
    while j < len(body):
        if body[j][0] == 999:
            j += 1
            continue
        if body[j][0] != 0:
            raise ValueError(f"line {body[j][2]}: expected entity")
        kind, line = body[j][1], body[j][2]
        j += 1
        values: list[Pair] = []
        while j < len(body) and body[j][0] != 0:
            values.append(body[j])
            j += 1
        entity = _Entity(kind, tuple(values), line)
        handle = entity.one(5, "")
        if handle:
            if not re.fullmatch(r"[0-9A-Fa-f]{1,16}", handle) or handle.upper() in handles:
                raise ValueError(f"line {line}: duplicate/invalid handle")
            handles.add(handle.upper())
        entities.append(entity)
        if len(entities) > options.max_entities:
            raise ValueError("DXF entity limit exceeded")
    return header, entities


def _codes(entity: _Entity, allowed: set[int], repeated: set[int] = frozenset()) -> None:
    counts: Counter[int] = Counter()
    control_depth = 0
    for c, v, line in entity.pairs:
        if c not in allowed | _META:
            raise GeometryError(f"line {line} {entity.kind}: unsupported group {c}")
        counts[c] += 1
        if counts[c] > 1 and c not in repeated | {100, 102, 330, 999}:
            raise ValueError(f"line {line} {entity.kind}: duplicate group {c}")
        if c == 102:
            if v.startswith("{"):
                control_depth += 1
            elif v == "}":
                control_depth -= 1
            else:
                raise ValueError("Invalid control group")
            if control_depth < 0:
                raise ValueError("Unbalanced control group")
    if control_depth:
        raise ValueError("Unbalanced control group")


def _flat(entity: _Entity) -> None:
    for c in (39, 40, 41, 43, 42):
        for code, value, line in entity.pairs:
            if code == c and _rational(value) != 0:
                raise GeometryError(f"line {line}: nonzero bulge/width/extrusion thickness")
    for c, default in ((210, "0"), (220, "0"), (230, "1")):
        if _rational(entity.one(c, default)) != Fraction(default):
            raise GeometryError(f"line {entity.line}: unsupported extrusion normal")


def _selected(entity: _Entity, options: DxfImportOptions) -> bool:
    space = _integer(entity.one(67, "0"))
    layout = entity.one(410, "")
    if space not in (0, 1) or layout and ((layout == "Model") != (space == 0)):
        raise ValueError(f"line {entity.line}: contradictory model/paper-space")
    return space == 0 and (options.layers is None or entity.one(8, "0") in options.layers)


def _normalize(raw: list[tuple[Point, Point, float, dict[str, Any]]], options: DxfImportOptions
               ) -> tuple[list[Segment], list[dict[str, Any]], list[dict[str, Any]], float]:
    tol = options.node_tolerance
    n = 2*len(raw)
    # Includes clustering, diameter, projections, intersections and core validation.
    if 10*n*n > options.max_pair_checks:
        raise ValueError("DXF pair-check budget exceeded")
    endpoints = [p for a, b, _, _ in raw for p in (a, b)]
    order = sorted(range(n), key=lambda i: (*endpoints[i], *(float(v).hex() for v in endpoints[i])))
    nodes, mapping = cluster_nodes([Node(*endpoints[i]) for i in order], tol)
    index = {original: mapping[j] for j, original in enumerate(order)}
    members = {k: [endpoints[i] for i in range(n) if index[i] == k] for k in range(len(nodes))}
    for group in members.values():
        if any(math.dist(a, b) > tol for a in group for b in group):
            raise TopologyError("Ambiguous transitive node cluster exceeds tolerance diameter")
    points = [node.coords for node in nodes]
    edges = [(index[2*i], index[2*i+1]) for i in range(len(raw))]
    if any(a == b for a, b in edges):
        raise GeometryError("Clustering collapsed a segment")
    cuts: dict[int, list[tuple[Fraction, int]]] = {}
    moves: dict[int, Point] = {}
    tol2 = Fraction(tol)**2
    for v, p in enumerate(points):
        for e, (a, b) in enumerate(edges):
            if v in (a, b):
                continue
            pa, pb = points[a], points[b]
            dx, dy = Fraction(pb[0])-Fraction(pa[0]), Fraction(pb[1])-Fraction(pa[1])
            px, py = Fraction(p[0])-Fraction(pa[0]), Fraction(p[1])-Fraction(pa[1])
            l2 = dx*dx+dy*dy
            u = (px*dx+py*dy)/l2
            if not 0 < u < 1 or (px*dy-py*dx)**2 > tol2*l2:
                continue
            if u*u*l2 <= tol2 or (1-u)**2*l2 <= tol2:
                raise TopologyError("Ambiguous endpoint-near-interior contact")
            if options.junction_policy == "reject":
                raise TopologyError("Endpoint/interior contact requires splitting")
            target = (_float(Fraction(pa[0])+u*dx), _float(Fraction(pa[1])+u*dy))
            if v in moves and moves[v] != target:
                raise TopologyError("Multiple incompatible T projections")
            moves[v] = target
            cuts.setdefault(e, []).append((u, v))
    for v, target in moves.items():
        if any(math.dist(old, target) > tol for old in members[v]):
            raise TopologyError("T snapping exceeds total endpoint movement budget")
        points[v] = target
    changes: list[dict[str, Any]] = []
    maximum = 0.0
    for i, old in enumerate(endpoints):
        new = points[index[i]]
        displacement = math.dist(old, new)
        maximum = max(maximum, displacement)
        if old != new:
            changes.append({"operation": "snap", "endpoint": i, "before": old, "after": new})
    pieces: list[tuple[Point, Point, float, dict[str, Any]]] = []
    for e, (a, b) in enumerate(edges):
        sequence = [a]+[v for _, v in sorted(cuts.get(e, []))]+[b]
        if len(sequence) > 2:
            changes.append({"operation": "split", "source": raw[e][3], "points": tuple(points[v] for v in sequence)})
        for u, v in zip(sequence, sequence[1:]):
            p, q = sorted((points[u], points[v]))
            if p == q:
                raise GeometryError("Normalization produced zero length")
            pieces.append((p, q, raw[e][2], raw[e][3]))
    if len(pieces) > options.max_output_segments or 40*len(pieces)**2 > options.max_pair_checks:
        raise ValueError("DXF normalized geometry resource limit exceeded")
    pieces.sort(key=lambda s: (s[0], s[1], s[2]))
    if len({(p, q) for p, q, _, _ in pieces}) != len(pieces):
        raise TopologyError("Duplicate normalized segment")
    # Exact orientations reject interior crossings before floating core checks.
    def cross(a: Point, b: Point, p: Point) -> Fraction:
        return ((Fraction(b[0])-Fraction(a[0]))*(Fraction(p[1])-Fraction(a[1])) -
                (Fraction(b[1])-Fraction(a[1]))*(Fraction(p[0])-Fraction(a[0])))
    for i, (a, b, _, _) in enumerate(pieces):
        for c, d, _, _ in pieces[i+1:]:
            c1, c2, c3, c4 = cross(a, b, c), cross(a, b, d), cross(c, d, a), cross(c, d, b)
            if c1*c2 < 0 and c3*c4 < 0:
                raise TopologyError("Interior X crossing is unsupported")
            if c1 == c2 == c3 == c4 == 0:
                axis = 0 if a[0] != b[0] else 1
                lo = max(min(a[axis], b[axis]), min(c[axis], d[axis]))
                hi = min(max(a[axis], b[axis]), max(c[axis], d[axis]))
                if lo < hi:
                    raise TopologyError("Collinear overlap is unsupported")
    canonical = {p: Node(*p, f"n{i}") for i, p in enumerate(sorted({p for a, b, _, _ in pieces for p in (a, b)}))}
    segments = [Segment(canonical[a], canonical[b], t, f"e{i}") for i, (a, b, t, _) in enumerate(pieces)]
    provenance = [dict(source, segment_id=f"e{i}") for i, (_, _, _, source) in enumerate(pieces)]
    return segments, provenance, changes, maximum


def loads_dxf(text: str, *, options: DxfImportOptions) -> DxfImportResult:
    if not isinstance(options, DxfImportOptions):
        raise TypeError("options must be DxfImportOptions")
    header, entities = _parse(text, options)
    version = header["$ACADVER"]
    if "$INSUNITS" in header:
        unit = {0: "unspecified", 1: "in", 2: "ft", 4: "mm", 5: "cm", 6: "m"}.get(_integer(header["$INSUNITS"]))
        if unit != options.source_length_unit:
            raise ValueError("INSUNITS contradicts source_length_unit")
    ratio = _UNITS[options.source_length_unit]/_UNITS[options.target_length_unit]
    plane_z = Fraction(options.plane_z)*ratio
    plane_tol = Fraction(options.plane_tolerance)
    selected: Counter[str] = Counter()
    excluded: Counter[str] = Counter()
    sources: list[dict[str, Any]] = []
    raw: list[tuple[Point, Point, float, dict[str, Any]]] = []
    used_h: set[str] = set()
    used_l: set[str] = set()
    max_z = Fraction(0)
    removed = 0

    def point(entity: _Entity, xc: int, yc: int, z: Fraction) -> Point:
        nonlocal max_z
        dz = abs(z*ratio-plane_z)
        if dz > plane_tol:
            raise GeometryError(f"line {entity.line}: out-of-plane coordinate")
        max_z = max(max_z, dz)
        return (_float(_rational(entity.one(xc))*ratio), _float(_rational(entity.one(yc))*ratio))

    i = 0
    while i < len(entities):
        entity = entities[i]
        i += 1
        vertices: list[_Entity] = []
        if entity.kind == "POLYLINE":
            while i < len(entities) and entities[i].kind == "VERTEX":
                vertices.append(entities[i])
                i += 1
            if i >= len(entities) or entities[i].kind != "SEQEND":
                raise ValueError(f"line {entity.line}: POLYLINE missing SEQEND")
            _codes(entities[i], set())
            i += 1
        elif entity.kind in ("VERTEX", "SEQEND"):
            raise ValueError(f"line {entity.line}: orphan {entity.kind}")
        if not _selected(entity, options):
            excluded[entity.kind] += 1
            continue
        if entity.kind not in ("LINE", "POLYLINE", "LWPOLYLINE") or entity.kind == "LWPOLYLINE" and version == "AC1009":
            raise GeometryError(f"line {entity.line}: unsupported {entity.kind} for {version}")
        selected[entity.kind] += 1
        layer, handle = entity.one(8, "0"), entity.one(5, "").upper()
        tm = options.thickness
        if handle in tm.by_handle:
            thickness, origin = Fraction(tm.by_handle[handle]), "handle"
            used_h.add(handle)
        elif layer in tm.by_layer:
            thickness, origin = Fraction(tm.by_layer[layer]), "layer"
            used_l.add(layer)
        elif tm.use_layer_name and layer.startswith("THICK_"):
            thickness, origin = _rational(layer[6:]), "layer_name"
        elif tm.default is not None:
            thickness, origin = Fraction(tm.default), "default"
        else:
            raise GeometryError(f"line {entity.line}: missing thickness for {layer}")
        t = _float(thickness*ratio)
        if t <= 0:
            raise GeometryError(f"line {entity.line}: thickness must be positive")
        source = {"line": entity.line, "handle": handle, "layer": layer, "entity": entity.kind}
        sources.append(dict(source, rule=origin, thickness=t))
        points: list[Point] = []
        closed = False
        if entity.kind == "LINE":
            _codes(entity, {10, 20, 30, 11, 21, 31, 39, 210, 220, 230})
            _flat(entity)
            points = [point(entity, 10, 20, _rational(entity.one(30, "0"))),
                      point(entity, 11, 21, _rational(entity.one(31, "0")))]
        elif entity.kind == "POLYLINE":
            _codes(entity, {10, 20, 30, 39, 40, 41, 66, 70, 210, 220, 230})
            _flat(entity)
            flags = _integer(entity.one(70, "0"))
            if flags < 0 or flags & ~129:
                raise GeometryError("Unsupported POLYLINE flags")
            if _rational(entity.one(10, "0")) or _rational(entity.one(20, "0")):
                raise GeometryError("POLYLINE dummy coordinates must be zero")
            closed = bool(flags & 1)
            elevation = _rational(entity.one(30, "0"))
            for v in vertices:
                _codes(v, {10, 20, 30, 40, 41, 42, 70})
                _flat(v)
                if _integer(v.one(70, "0")) or _rational(v.one(30, "0")):
                    raise GeometryError("Unsupported VERTEX flags or z")
                if v.one(8, "0") not in ("0", layer):
                    raise GeometryError("VERTEX layer contradicts POLYLINE")
                points.append(point(v, 10, 20, elevation))
        else:
            _codes(entity, {10, 20, 38, 39, 40, 41, 42, 43, 70, 90, 91, 210, 220, 230}, {10, 20, 40, 41, 42, 91})
            # Width/bulge checked per vertex and globally without singleton lookup.
            _flat(entity)
            flags = _integer(entity.one(70, "0"))
            if flags < 0 or flags & ~129:
                raise GeometryError("Unsupported LWPOLYLINE flags")
            closed = bool(flags & 1)
            elevation = _rational(entity.one(38, "0"))
            current: list[Pair] = []
            vertex_records: list[_Entity] = []
            for pair in entity.pairs:
                if pair[0] == 10:
                    if current:
                        vertex_records.append(_Entity("LWVERTEX", tuple(current), current[0][2]))
                    current = [pair]
                elif pair[0] in (20, 40, 41, 42, 91):
                    if not current:
                        raise ValueError("LWPOLYLINE vertex field before x")
                    current.append(pair)
            if current:
                vertex_records.append(_Entity("LWVERTEX", tuple(current), current[0][2]))
            if len(vertex_records) != _integer(entity.one(90)):
                raise ValueError("LWPOLYLINE vertex count mismatch")
            for v in vertex_records:
                _codes(v, {10, 20, 40, 41, 42, 91})
                points.append(point(v, 10, 20, elevation))
        if closed and len(points) > 1 and points[0] == points[-1]:
            points.pop()
            removed += 1
        if len(points) < (3 if closed else 2):
            raise GeometryError("Insufficient polyline vertices")
        pairs = list(zip(points, points[1:]))
        if closed:
            pairs.append((points[-1], points[0]))
        for edge, (a, b) in enumerate(pairs):
            if a == b:
                raise GeometryError("Zero-length DXF edge")
            raw.append((a, b, t, dict(source, edge=edge)))
            if len(raw) > options.max_output_segments:
                raise ValueError("DXF segment limit exceeded")
    if set(options.thickness.by_handle) != used_h or set(options.thickness.by_layer) != used_l:
        raise ValueError("Unused thickness mapping")
    if not raw:
        raise GeometryError("Empty DXF selection")
    segments, provenance, changes, max_move = _normalize(raw, options)
    topology = extract_mixed_topology(segments, node_tolerance=options.node_tolerance, safety_factor=options.safety_factor)
    kind = topology.topology_type if options.section_kind == "auto" else options.section_kind
    cls = {"open": Section, "closed": ClosedSection, "mixed": MixedSection}[kind]
    kwargs: dict[str, Any] = {"node_tolerance": options.node_tolerance, "validate": True}
    if cls is not Section:
        kwargs["safety_factor"] = options.safety_factor
    section = cls(segments, **kwargs)
    report = DxfImportReport(version, options.source_length_unit, options.target_length_unit, str(ratio),
                             MappingProxyType(dict(selected)), MappingProxyType(dict(excluded)),
                             tuple(MappingProxyType(v) for v in provenance),
                             tuple(MappingProxyType(v) for v in changes), _float(max_z), max_move, removed,
                             len(topology.canonical_nodes), len(segments), 1, topology.cycle_rank, cls.__name__,
                             tuple(MappingProxyType(v) for v in sources))
    return DxfImportResult(section, report)


def read_dxf(path: str | Path, *, options: DxfImportOptions, encoding: str = "ascii") -> DxfImportResult:
    with Path(path).open("rb") as stream:
        raw = stream.read(options.max_bytes+1)
    if len(raw) > options.max_bytes:
        raise ValueError("DXF byte limit exceeded")
    text = raw.decode(encoding, errors="strict")
    header, _ = _parse(text, options)
    if "$DWGCODEPAGE" in header:
        page = header["$DWGCODEPAGE"]
        name = "cp" + page[5:] if page.startswith("ANSI_") else page
        if codecs.lookup(name).name != codecs.lookup(encoding).name:
            raise ValueError("DWGCODEPAGE contradicts explicit encoding")
    return loads_dxf(text, options=options)

