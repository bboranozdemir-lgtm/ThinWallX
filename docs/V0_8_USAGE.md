# Sectalix v0.8 interchange format

This document describes the v0.8 JSON interchange format, the supported ASCII DXF import path, and the optional plotting interface.

The wire-format identifier `"format": "thinwallx"` is intentionally retained for backward compatibility with files created before the project was renamed to Sectalix. It is a serialization identifier, not the current package name.

## JSON interchange

```python
from sectalix import Node, Segment, Section
from sectalix.serialization import UnitSystem, to_json, from_json

section = Section([
    Segment(Node(0.0, 0.0), Node(2.0, 0.0), 0.01, id="flange"),
    Segment(Node(0.0, 0.0), Node(0.0, 1.0), 0.01, id="web"),
])

units = UnitSystem(length="mm", force="N")
wire = to_json(section, units=units)
document = from_json(wire)

restored_section = document.value
assert document.units == units
assert restored_section.segments == section.segments
```

`from_dict`, `from_json`, and `read_json` return a `DecodedDocument` containing both the decoded value and its unit metadata. Unit labels are descriptive and do not rescale numerical values.

Physical floating-point values are serialized using canonical `float.hex` strings so that supported IEEE-754 binary64 values can be reconstructed exactly. Ordinary decimal JSON numbers are not accepted in fields defined by the v0.8 schema as encoded floating-point values.

The JSON schema is available at [`schemas/sectalix-0.8.schema.json`](../schemas/sectalix-0.8.schema.json). Additional cross-reference, topology, and result-consistency checks are performed by the Python decoder.

Supported object identifiers are limited to serializable scalar or tuple-like values defined by the codec. Arbitrary Python objects are not serialized or reconstructed.

`AppliedLoads` and `StressRecoveryResult` can use the same interchange system. A stress-result document does not automatically embed the source geometry and load document; applications that require a complete audit record should store those inputs alongside the result.

`write_json` does not overwrite an existing file unless `overwrite=True` is requested explicitly. The decoder also applies configurable limits to input size, nesting, and record counts.

## ASCII DXF import

```python
from sectalix.dxf import DxfImportOptions, ThicknessMap, read_dxf

options = DxfImportOptions(
    source_length_unit="mm",
    target_length_unit="mm",
    thickness=ThicknessMap(by_layer={"WEB": 2.0}, default=1.0),
    node_tolerance=1e-9,
    section_kind="auto",
)

imported = read_dxf("section.dxf", options=options)
section = imported.section
report = imported.report
```

Thickness can be assigned by entity handle, layer mapping, a `THICK_<value>` layer name, or an explicit fallback thickness. More specific mappings take precedence. Thickness mapping is interpreted in the source DXF length unit, while node tolerance is interpreted in the target unit.

The importer intentionally supports a restricted subset of ASCII DXF suitable for straight centerline geometry:

- AC1009: `LINE` and straight 2D `POLYLINE`
- AC1015: the above plus straight `LWPOLYLINE`

Curved bulge segments, unsupported widths or extrusion settings, 3D polylines, `INSERT`, arcs, and other unsupported entities are rejected rather than approximated silently.

The importer is not a general-purpose CAD kernel. Its purpose is to extract supported straight centerline geometry with explicit thickness information and pass that geometry through Sectalix validation.

Endpoint snapping and supported open-junction handling are controlled by the import options. Any geometric changes made by the importer are recorded in the returned import report. Intersections or overlaps that cannot be handled unambiguously are rejected.

Automatic topology classification returns an open `Section`, a closed `ClosedSection`, or a `MixedSection` according to the validated graph. A caller may request a specific section kind; incompatible geometry then raises an error.

## Plotting

Plotting is optional:

```console
python -m pip install -e ".[plots]"
```

Example:

```python
from sectalix import AppliedLoads
from sectalix.plotting import (
    plot_geometry,
    plot_shear_flow,
    plot_stresses,
    GeometryPlotOptions,
    StressPlotOptions,
)

plot_geometry(
    section,
    "geometry.svg",
    options=GeometryPlotOptions(
        show_shear_center=True,
        show_principal_axes=True,
    ),
)

flow = section.calculate_shear_flow(vx=0.0, vy=1.0)
plot_shear_flow(section, flow, "shear.png")

stress = section.calculate_stresses(AppliedLoads(N=1.0, Mx=0.1))
plot_stresses(
    section,
    stress,
    "stress.svg",
    options=StressPlotOptions(quantity="sigma_vm", abscissa="segment"),
)
```

Plotting functions create and close their own figures and do not require an interactive GUI. PNG and SVG output are supported.

Thickness bands are a visualization aid for the centerline model; they should not be interpreted as a finite-width solid-section mesh. Reported mechanical extrema come from the analysis result, not from visual plot sampling.

Byte-identical image output is only expected within a consistent rendering environment. Font, Matplotlib, and platform differences can change the rendered bytes without changing the underlying mechanical result.

## Examples and tests

The source repository includes:

- `examples/sample_sections/` — representative open, closed, and mixed DXF inputs
- `examples/sample_calculation/` — selected calculation output and plots
- `examples/quickstart_api.py` — a compact Python API example

From a source checkout:

```console
python examples/quickstart_api.py
sectalix inspect examples/sample_sections/rectangle.dxf
python -m pytest -W error
```

Passing tests provide regression and benchmark evidence for the documented implementation. They do not establish correctness for every possible geometry or replace independent engineering validation for a safety-critical application.
