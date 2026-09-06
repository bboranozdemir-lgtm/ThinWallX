"""Sectalix argparse CLI. Stdout is data only; errors use stderr."""
from __future__ import annotations

import argparse
from decimal import Decimal, InvalidOperation
import math
import re
from pathlib import Path
import sys
from typing import Sequence

from sectalix import __version__, Section, AppliedLoads
from sectalix.dxf import DxfImportOptions, ThicknessMap, read_dxf
from sectalix.exceptions import GeometryError, TopologyError, SingularSectionError
from sectalix.reporting import inspect_section, inspection_json, inspection_text, write_report
from sectalix.serialization import (
    DecodeLimits, DecodedDocument, UnitSystem, read_json, from_json, to_json, _atomic_write,
)

_LOAD_FIELDS = ("N", "Vx", "Vy", "Mx", "My", "Tsv", "B", "M_omega", "sigma_yield")
_UNITS = ("m", "mm", "cm", "in", "ft", "unspecified")


class UsageError(Exception):
    """CLI option error, distinct from a malformed document."""


class _Parser(argparse.ArgumentParser):
    def __init__(self, *args: object, **kwargs: object) -> None:
        super().__init__(*args, **kwargs)
        self._negative_number_matcher = re.compile(r"^-(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?$")

    def error(self, message: str) -> None:
        raise UsageError(message)


def _number(token: str) -> float:
    try:
        value = Decimal(token)
    except InvalidOperation as exc:
        raise argparse.ArgumentTypeError("expected a real number") from exc
    if not value.is_finite():
        raise GeometryError("CLI numeric inputs must be finite")
    result = float(value)
    if not math.isfinite(result):
        raise OverflowError("CLI scalar exceeds float64")
    if value and result == 0:
        raise FloatingPointError("Nonzero CLI scalar underflows float64")
    return result


def parser() -> argparse.ArgumentParser:
    root = _Parser(prog="sectalix", description="Thin-wall inspection, analysis and reports.")
    root.add_argument("--version", action="version", version=f"Sectalix {__version__}")
    commands = root.add_subparsers(dest="command", required=True)
    for name in ("inspect", "analyze", "plot", "convert-dxf"):
        p = commands.add_parser(name)
        p.add_argument("input_file", help="Section JSON/DXF path; - for JSON stdin")
        p.add_argument("--output", required=name in ("plot", "convert-dxf"))
        p.add_argument("--overwrite", action="store_true")
        p.add_argument("--thickness", type=_number, help="Fallback DXF thickness in source units")
        if name == "convert-dxf":
            p.add_argument("--source-unit", choices=_UNITS, default="mm")
            p.add_argument("--target-unit", choices=_UNITS, default="mm")
            p.add_argument("--section-kind", choices=("auto", "open", "closed", "mixed"), default="auto")
        else:
            p.add_argument("--dxf-unit", choices=_UNITS, default="mm")
        if name == "inspect":
            p.add_argument("--json", action="store_true")
        elif name == "analyze":
            p.add_argument("--loads", help="AppliedLoads JSON; scalar flags override named fields")
            for f in _LOAD_FIELDS:
                p.add_argument("--"+f.replace("_", "-"), dest=f, type=_number, default=None)
            p.add_argument("--plots-dir")
            p.add_argument("--report")
            p.add_argument("--stdout", action="store_true")
        elif name == "plot":
            p.add_argument("--show-thickness", action=argparse.BooleanOptionalAction, default=True)
            p.add_argument("--show-axes", action="store_true")
            p.add_argument("--show-ids", action="store_true")
            p.add_argument("--show-shear-center", action="store_true")
    return root


def _stdin_document() -> DecodedDocument:
    limit = DecodeLimits().max_bytes
    if hasattr(sys.stdin, "buffer"):
        raw = sys.stdin.buffer.read(limit+1)
        if len(raw) > limit:
            raise ValueError("stdin exceeds byte budget")
        text = raw.decode("utf-8")
    else:
        text = sys.stdin.read(limit+1)
    return from_json(text)


def _section(args: argparse.Namespace) -> DecodedDocument:
    name = args.input_file
    if name == "-":
        if args.command == "convert-dxf":
            raise UsageError("convert-dxf requires a DXF file")
        document = _stdin_document()
    elif Path(name).suffix.lower() == ".json":
        if args.command == "convert-dxf":
            raise UsageError("convert-dxf requires a DXF file")
        document = read_json(name)
    elif Path(name).suffix.lower() == ".dxf":
        source = args.source_unit if args.command == "convert-dxf" else args.dxf_unit
        target = args.target_unit if args.command == "convert-dxf" else source
        if args.thickness is not None and args.thickness <= 0:
            raise GeometryError("DXF thickness must be positive")
        imported = read_dxf(name, options=DxfImportOptions(
            source, target, ThicknessMap(default=args.thickness),
            section_kind=args.section_kind if args.command == "convert-dxf" else "auto",
        ))
        document = DecodedDocument(imported.section, UnitSystem(target, "unspecified"))
    else:
        raise ValueError("Input extension must be .json or .dxf")
    if not isinstance(document.value, Section):
        raise ValueError("Input document must contain a section")
    return document


def _loads(args: argparse.Namespace, units: UnitSystem) -> tuple[AppliedLoads, UnitSystem]:
    if args.loads == "-" and args.input_file == "-":
        raise UsageError("Section and loads cannot both consume stdin")
    values = {f: None if f == "sigma_yield" else 0.0 for f in _LOAD_FIELDS}
    if args.loads:
        doc = _stdin_document() if args.loads == "-" else read_json(args.loads)
        if not isinstance(doc.value, AppliedLoads):
            raise ValueError("--loads must contain AppliedLoads")
        if doc.units.length != units.length:
            raise ValueError("Section and loads length units differ; automatic conversion is forbidden")
        if units.force != "unspecified" and doc.units.force != units.force:
            raise ValueError("Section and loads force units differ")
        units = UnitSystem(units.length, doc.units.force)
        values = {f: getattr(doc.value, f) for f in _LOAD_FIELDS}
    explicit = {f: getattr(args, f) for f in _LOAD_FIELDS if getattr(args, f) is not None}
    if not args.loads and not explicit:
        raise UsageError("analyze requires --loads or at least one load flag")
    values.update(explicit)
    return AppliedLoads(**values), units


def _output(text: str, path: str | None, overwrite: bool, also_stdout: bool = False) -> None:
    if path and path != "-":
        _atomic_write(path, lambda p: p.write_text(text, encoding="utf-8", newline="\n"), overwrite)
    if not path or path == "-" or also_stdout:
        sys.stdout.write(text)


def _preflight(args: argparse.Namespace) -> None:
    """Files are individually atomic; reject collisions before publishing any."""
    outputs: list[Path] = []
    if args.output and args.output != "-":
        outputs.append(Path(args.output))
    if getattr(args, "report", None):
        outputs.append(Path(args.report))
    plots = Path(args.plots_dir) if getattr(args, "plots_dir", None) else None
    if plots:
        outputs.extend(plots/name for name in ("geometry.png", "shear_flow.png", "stress_vm.png"))
    resolved = [p.resolve() for p in outputs]
    if len(set(resolved)) != len(resolved):
        raise UsageError("Output paths must be distinct")
    inputs = [Path(n).resolve() for n in (args.input_file, getattr(args, "loads", None)) if n and n != "-"]
    if any(p in inputs for p in resolved):
        raise UsageError("Output must not replace an input file")
    for p in outputs:
        if p.exists() and (not args.overwrite or not p.is_file()):
            raise FileExistsError(p)
        if not p.parent.is_dir() and not (plots and p.parent == plots):
            raise FileNotFoundError(f"Output parent does not exist: {p.parent}")


def _execute(args: argparse.Namespace) -> None:
    if args.command == "plot" and args.output == "-":
        raise UsageError("plot requires a PNG/SVG file; binary stdout is unsupported")
    _preflight(args)
    document = _section(args)
    section, units = document.value, document.units
    if args.command == "inspect":
        info = inspect_section(section)
        text = inspection_json(section, units, info)+"\n" if args.json else inspection_text(info, units)
        _output(text, args.output, args.overwrite)
    elif args.command == "convert-dxf":
        _output(to_json(document)+"\n", args.output, args.overwrite)
    elif args.command == "plot":
        from sectalix.plotting import plot_geometry, GeometryPlotOptions
        plot_geometry(section, args.output, units=units, overwrite=args.overwrite,
                      options=GeometryPlotOptions(show_thickness=args.show_thickness,
                                                  show_principal_axes=args.show_axes,
                                                  show_node_ids=args.show_ids, show_segment_ids=args.show_ids,
                                                  show_shear_center=args.show_shear_center))
    else:
        loads, units = _loads(args, units)
        result = section.calculate_stresses(loads)
        output = to_json(result, units=units)+"\n"
        info = inspect_section(section) if args.report else None
        images: dict[str, Path] = {}
        if args.plots_dir:
            from sectalix.plotting import plot_geometry, plot_shear_flow, plot_stresses
            folder = Path(args.plots_dir)
            flow = section.calculate_shear_flow(vx=loads.Vx, vy=loads.Vy)
            folder.mkdir(parents=True, exist_ok=True)
            images = {"Geometry": folder/"geometry.png", "Transverse shear flow": folder/"shear_flow.png",
                      "Von Mises stress": folder/"stress_vm.png"}
            plot_geometry(section, images["Geometry"], units=units, overwrite=args.overwrite)
            plot_shear_flow(section, flow, images["Transverse shear flow"], units=units, overwrite=args.overwrite)
            plot_stresses(section, result, images["Von Mises stress"], units=units, overwrite=args.overwrite)
        if args.report:
            write_report(section, loads, result, args.report, units=units, input_path=args.input_file,
                         plots=images, overwrite=args.overwrite, inspection=info)
        _output(output, args.output, args.overwrite, args.stdout)


def main(argv: Sequence[str] | None = None) -> int:
    """Return the application exit code; argparse help/version exit successfully."""
    try:
        args = parser().parse_args(argv)
        _execute(args)
        return 0
    except UsageError as exc:
        code, message = 1, str(exc)
    except GeometryError as exc:
        code, message = 3, str(exc)
    except TopologyError as exc:
        code, message = 4, str(exc)
    except SingularSectionError as exc:
        code, message = 5, str(exc)
    except (OverflowError, FloatingPointError) as exc:
        code, message = 6, str(exc)
    except (OSError, ValueError, TypeError, ImportError) as exc:
        code, message = 2, str(exc)
    print(f"sectalix: {message}", file=sys.stderr)
    return code
