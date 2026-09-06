"""Calculation-sheet presentation; mechanics are delegated to frozen APIs."""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from fractions import Fraction
import html
import json
import math
import os
from pathlib import Path
from typing import Mapping
from urllib.parse import quote

from sectalix import Section, ClosedSection, MixedSection
from sectalix.serialization import UnitSystem, _atomic_write, _f64, to_dict
from sectalix.stress import AppliedLoads, StressRecoveryResult


@dataclass(frozen=True)
class SectionInspection:
    topology: str
    node_count: int
    segment_count: int
    properties: Mapping[str, float]


def _finite(value: float) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise OverflowError("Reported physical quantity exceeds float64")
    return result


def _rational_float(value: Fraction) -> float:
    result = _finite(float(value))
    if value and result == 0:
        raise FloatingPointError("Nonzero reported quantity underflows float64")
    return result


def inspect_section(section: Section) -> SectionInspection:
    """Obtain the required properties from frozen public APIs, once."""
    if type(section) not in (Section, ClosedSection, MixedSection):
        raise TypeError("Expected a Sectalix section")
    section.validate()
    if isinstance(section, MixedSection):
        topology = section.mixed_topology.topology_type
    else:
        topology = "closed" if isinstance(section, ClosedSection) else "open"
    warp = section.torsion_properties()
    props = {
        "A": section.area, "cx": section.cx, "cy": section.cy,
        "Ix": section.Ix, "Iy": section.Iy, "Ixy": section.Ixy,
        "I1": section.I1, "I2": section.I2, "theta_p": section.theta_p,
        "J": warp.J, "sx": warp.shear_center[0], "sy": warp.shear_center[1],
        "Cw": warp.Cw,
        "total_length": _rational_float(sum((Fraction(s.length) for s in section.segments), Fraction())),
        "t_min": min(s.t for s in section.segments), "t_max": max(s.t for s in section.segments),
    }
    return SectionInspection(topology, len(section.nodes), len(section.segments),
                             {k: _finite(v) for k, v in props.items()})


def inspection_json(section: Section, units: UnitSystem, inspection: SectionInspection | None = None) -> str:
    """v1 observation envelope, not an additional v0.8 serializable kind."""
    info = inspection or inspect_section(section)
    record = {
        "format": "sectalix-inspection", "schema_version": "1.0",
        "number_encoding": "float64-hex", "topology": info.topology,
        "node_count": info.node_count, "segment_count": info.segment_count,
        "units": {"length": units.length, "force": units.force},
        "properties": {k: _f64(v) for k, v in info.properties.items()},
        "section": to_dict(section, units=units),
    }
    return json.dumps(record, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _escape(value: object) -> str:
    """Prevent supplied paths and IDs from introducing Markdown/HTML structure."""
    text = str(value).replace("\r", " ").replace("\n", " ")
    special = "\\*_{}[]()|!#$" + chr(96)
    return "".join(f"&#{ord(c)};" if c in special else html.escape(c, quote=True) for c in text)


def _number(value: float) -> str:
    return format(_finite(value), ".17g")


def _table(rows: list[tuple[str, str, str]]) -> str:
    return "\n".join(["| Quantity | Value | Unit |", "|---|---:|---|"] +
                     [f"| {_escape(k)} | {v} | {_escape(u)} |" for k, v, u in rows])


def property_rows(info: SectionInspection, units: UnitSystem) -> list[tuple[str, str, str]]:
    length = units.length
    powers = {"A": 2, "Ix": 4, "Iy": 4, "Ixy": 4, "I1": 4, "I2": 4, "J": 4, "Cw": 6}
    return [(k, _number(v), "rad" if k == "theta_p" else
             f"{length}^{powers[k]}" if k in powers else length) for k, v in info.properties.items()]


def inspection_text(info: SectionInspection, units: UnitSystem) -> str:
    return (f"Sectalix section: {info.topology}\nNodes: {info.node_count}\nSegments: {info.segment_count}\n"
            + "\n".join(f"{k}: {v} [{u}]" for k, v, u in property_rows(info, units)) + "\n")


def calculation_report(
    section: Section, loads: AppliedLoads, result: StressRecoveryResult, *,
    units: UnitSystem | None = None, input_path: str = "<memory>",
    report_path: str | Path = "report.md", plots: Mapping[str, str | Path] | None = None,
    generated_at: datetime | None = None, inspection: SectionInspection | None = None,
) -> str:
    """Render a GFM sheet with 17-digit values and explicit model limitations."""
    unit = units or UnitSystem()
    if not isinstance(unit, UnitSystem) or not isinstance(loads, AppliedLoads):
        raise TypeError("Expected UnitSystem and AppliedLoads")
    to_dict(result)
    ids = [s.id if s.id is not None else i for i, s in enumerate(section.segments)]
    if len(set(ids)) != len(ids) or set(ids) != set(result.segment_profiles):
        raise ValueError("Report result segment IDs do not match section")
    for key, segment in zip(ids, section.segments):
        profile = result.segment_profiles[key]
        if profile.length != segment.length or profile.thickness != segment.t:
            raise ValueError("Report result profile geometry mismatch")
    if (loads.sigma_yield is None) != (result.load_factor is None):
        raise ValueError("Report load factor does not match yield-stress availability")
    info = inspection or inspect_section(section)
    stamp = generated_at or datetime.now(timezone.utc)
    if stamp.tzinfo is None or stamp.utcoffset() is None:
        raise ValueError("Report timestamp must be timezone-aware")
    lf = "Not supplied" if result.load_factor is None else (
        "Unbounded (zero recovered stress; not a general capacity guarantee)"
        if result.load_factor == math.inf else _number(result.load_factor))
    length, force = unit.length, unit.force
    load_rows = []
    for name in ("N", "Vx", "Vy", "Mx", "My", "Tsv", "B", "M_omega", "sigma_yield"):
        value = getattr(loads, name)
        dim = force if name in ("N", "Vx", "Vy") else (
            f"{force}/{length}^2" if name == "sigma_yield" else
            f"{force} {length}^2" if name == "B" else f"{force} {length}")
        load_rows.append((name, "Not supplied" if value is None else _number(value), dim))
    peak_rows = [
        ("sigma_vm_max", _number(result.max_sigma_vm), f"{force}/{length}^2"),
        ("peak_segment_id", _escape(repr(result.peak_segment_id)), "ID"),
        ("s_peak", _number(result.peak_s), length),
        ("x_peak", _number(result.peak_location_xy[0]), length),
        ("y_peak", _number(result.peak_location_xy[1]), length),
        ("load_factor", lf, "dimensionless"),
    ]
    balance = ["| Resultant | Applied | Recovered | Recovered - applied | Unit |",
               "|---|---:|---:|---:|---|"]
    for name in ("N", "Mx", "My", "B"):
        applied, recovered = getattr(loads, name), getattr(result, "resultant_"+name)
        delta = _rational_float(Fraction(recovered)-Fraction(applied))
        dim = force if name == "N" else f"{force} {length}^2" if name == "B" else f"{force} {length}"
        balance.append(f"| {name} | {_number(applied)} | {_number(recovered)} | {_number(delta)} | {_escape(dim)} |")
    geometry_rows = [
        ("segments", str(info.segment_count), "count"), ("nodes", str(info.node_count), "count"),
        ("total_length", _number(info.properties["total_length"]), length),
        ("t_min", _number(info.properties["t_min"]), length),
        ("t_max", _number(info.properties["t_max"]), length),
    ]
    images = []
    for label, path in (plots or {}).items():
        image_path = Path(path).absolute()
        if not image_path.is_file() or image_path.suffix.lower() not in (".png", ".svg"):
            raise ValueError("Report image must be an existing PNG or SVG")
        relative = os.path.relpath(image_path, Path(report_path).absolute().parent).replace(os.sep, "/")
        images.append(f"![{_escape(label)}](./{quote(relative, safe='/._-')})")
    sections = [
        "# Sectalix v1.0.1 Calculation Report",
        f"Generated: {stamp.isoformat()}\n\nInput: {_escape(input_path)}\n\n"
        f"Units: length={_escape(length)}, force={_escape(force)}",
        "## Geometry and model assumptions\n\n"
        f"Topology: {_escape(info.topology)}\n\n"+_table(geometry_rows)+"\n\n"
        "Linear-elastic thin-wall centerline model. Thickness is constant per segment, "
        "not necessarily uniform across the section. Corner overlaps and root radii are neglected. "
        "No buckling, plasticity, fatigue, connection or code-compliance verification is performed.",
        "## Section characteristics\n\n"+_table(property_rows(info, unit)),
        "## Applied loads\n\n"+_table(load_rows),
        "## Critical stress and elastic yield factor\n\n"+_table(peak_rows)+"\n\n"
        "The load factor is the proportional elastic first-yield multiplier. It is not a "
        "design-code safety factor and does not establish regulatory approval.",
        "## Recovered resultants and equilibrium differences\n\n"+"\n".join(balance),
        "## Graphical appendices\n\n"+("\n\n".join(images) if images else "No plots requested."),
    ]
    return "\n\n".join(sections)+"\n"


def write_report(
    section: Section, loads: AppliedLoads, result: StressRecoveryResult, path: str | Path,
    *, units: UnitSystem | None = None, input_path: str = "<memory>",
    plots: Mapping[str, str | Path] | None = None, generated_at: datetime | None = None,
    overwrite: bool = False, inspection: SectionInspection | None = None,
) -> Path:
    text = calculation_report(section, loads, result, units=units, input_path=input_path,
                              report_path=path, plots=plots, generated_at=generated_at, inspection=inspection)
    return _atomic_write(path, lambda p: p.write_text(text, encoding="utf-8", newline="\n"), overwrite)
