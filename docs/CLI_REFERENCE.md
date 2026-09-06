# ThinWallX 1.0 CLI reference

## Installation

Python 3.10 or newer is required. From a source checkout run `python -m pip install .`; for plots use `python -m pip install ".[plots]"`. NumPy is the sole required runtime dependency. Matplotlib is optional and lazily loaded. Tests additionally need pytest. Both `thinwallx` and `python -m thinwallx` invoke the same entry point. `thinwallx --version` prints `ThinWallX 1.0.0`; `--help` and each subcommand's `--help` describe accepted options.

## Input and output rules

Section input is a v0.8 ThinWallX JSON document or supported ASCII DXF, selected by case-insensitive filename extension. A positional `-` means UTF-8 JSON stdin, not DXF. See [the frozen I/O guide](V0_8_USAGE.md) for schema, supported entities and thickness mapping. Coordinates and loads must use a coherent unit system; JSON metadata is not a conversion request.

All four commands accept `--output PATH`, `--overwrite` and `--thickness VALUE`. Output is required for plot and convert-dxf. Thickness is a fallback for DXF entities without an assigned thickness, in source length units; entity/layer resolution follows the existing importer. It must be positive. Normal section-consuming commands accept `--dxf-unit UNIT` (default mm), used as both source and target unit. convert-dxf instead has independent source and target options. Units are m, mm, cm, in, ft or unspecified; incompatible unspecified conversions are rejected by the importer.

Existing targets are protected unless --overwrite is explicit. Outputs may not replace the input section or loads file even with --overwrite. Duplicate output paths are rejected before publication. Individual files use atomic publication; the entire multi-file analysis is not a transaction. Completed plots can remain if a later output fails. Ordinary output parents must exist; --plots-dir may create its directory tree.

Stdout contains only requested data, never progress messages. Diagnostics go to stderr. Numeric flags support negative scientific notation, including `--Mx -1e-2`. Nonfinite input and out-of-range conversion are rejected. No interactive window or GUI loop is started.

## inspect

`thinwallx inspect INPUT [--json] [--output PATH] [--dxf-unit mm] [--thickness VALUE] [--overwrite]`

The default human-readable output lists topology, node and segment counts, A, cx, cy, Ix, Iy, Ixy, I1, I2, theta_p (radians), J, sx, sy, Cw, total length and minimum/maximum segment thickness. Without --output, it goes to stdout. --output - also selects stdout.

--json produces a **v1.0 observation envelope**, not a new frozen v0.8 codec kind. Its format is thinwallx-inspection, schema_version is 1.0, number_encoding is float64-hex, and properties contain F64-encoded values. It embeds an unmodified v0.8 document under `section`. Extract that member for subsequent section input; the whole inspection envelope is deliberately not accepted by the v0.8 decoder.

Inspection requests all listed properties. A section whose shear center or warping calculation is singular or unrepresentable fails explicitly even if its area alone could be computed. The CLI does not invent absent properties or hide a failure as null.

## analyze

`thinwallx analyze INPUT [--loads LOADS.json] [load flags] [--output RESULT.json] [--stdout] [--plots-dir DIR] [--report REPORT.md] [--dxf-unit mm] [--thickness VALUE] [--overwrite]`

Load flags are `--N`, `--Vx`, `--Vy`, `--Mx`, `--My`, `--Tsv`, `--B`, `--M-omega` and `--sigma-yield`. An explicit flag overrides only its corresponding loads-file field. Without a file, unspecified force/moment fields are zero and yield stress is absent. At least a file or one explicit flag is required; --N 0 is a valid zero-load request.

The loads document must be the AppliedLoads kind. Its length metadata must match the section; force metadata must match if the section specifies it. An unspecified section force label adopts the load document's label. No numerical rescaling occurs. Scalar flags use that same unit system. Section and loads cannot both consume stdin; --loads - is supported with a section filename.

Output is the unchanged v0.8 StressRecoveryResult document, including exact profile coefficients and resultant fields. Without --output, JSON goes to stdout. --stdout additionally echoes JSON when a file output is chosen. This is not a second computation.

--plots-dir saves geometry.png, shear_flow.png and stress_vm.png. The shear-flow image shows **transverse Vx/Vy flow**, not the sum of torsion and secondary warping effects. The stress image uses the complete recovered stress result. --report writes a Markdown calculation sheet and links the requested images relatively. Report creation requests full section properties and can therefore fail on an unavailable Cw even when a more limited calculation could succeed.

The report includes the timestamp with timezone, input path, units, assumptions, geometry, properties, all loads, maximum von Mises stress and location, yield multiplier, and recovered-minus-applied N/Mx/My/B. The yield multiplier is not regulatory certification. Reports differ by generation time; direct reporting API callers can supply an aware datetime for reproducibility.

## plot

`thinwallx plot INPUT --output FIGURE.png [--show-thickness | --no-show-thickness] [--show-axes] [--show-ids] [--show-shear-center] [--dxf-unit mm] [--thickness VALUE] [--overwrite]`

Only PNG and SVG file output is supported; binary stdout is not. Thickness bands default on. IDs show both nodes and segments. Axes show principal directions. Shear-center calculation is requested only by its option. Equal aspect ratio and frozen headless plotting behavior apply. Without Matplotlib, nonplot commands still work; a requested plot fails with code 2.

## convert-dxf

`thinwallx convert-dxf INPUT.dxf --output SECTION.json [--source-unit mm] [--target-unit mm] [--thickness VALUE] [--section-kind auto|open|closed|mixed] [--overwrite]`

The importer performs supported length conversion, validation and topology classification. Defaults are mm/mm/auto. --output - sends JSON to stdout. Input must be a DXF file. No CAD dependency, GUI, curve tessellation or unsupported entity approximation is introduced.

For a minimal input, create an ASCII DXF ENTITIES section containing a LINE on layer THICK_0.01, with group codes 10/20 for its start and 11/21 for its end. Valid mechanical inspection also requires nonsingular geometry: use the complete [open L fixture](../tests/fixtures/dxf/open_l_r12.dxf), not a single line. Convert it to obtain a complete canonical JSON example:

`thinwallx convert-dxf tests/fixtures/dxf/open_l_r12.dxf --output section.json`

The resulting document can be passed unchanged to inspect/analyze. Canonical JSON represents float64 values as F64 hexadecimal strings; replacing them by ordinary JSON numbers is not valid. [The JSON schema](../schemas/thinwallx-0.8.schema.json) specifies all required keys and rejection rules.

## Exit status

| Code | Meaning | Typical cause |
|---|---|---|
| 0 | Success | Help, version or completed output |
| 1 | Invocation | Unknown flag, missing load specification, colliding paths |
| 2 | I/O or format | Missing file, malformed JSON, unit mismatch, missing Matplotlib |
| 3 | Geometry | Negative thickness, invalid finite geometry |
| 4 | Topology | Disconnected or invalid graph |
| 5 | Singular mechanics | Unsupported nonzero bimoment on zero warping resistance |
| 6 | Numerical range | Overflow or nonzero underflow |

Exception category follows the frozen API; this interface does not reclassify unsupported geometry silently. Unexpected programming errors are not swallowed.

## Automation examples

Bash, with explicit failure propagation:

```sh
set -e
thinwallx convert-dxf shape.dxf --output section.json --thickness 0.01
cat section.json | thinwallx analyze - --N 100 --Mx -1e-2 > result.json
thinwallx analyze section.json --N 100 --plots-dir figures --report report.md --output result2.json
```

PowerShell:

```powershell
$env:PYTHONIOENCODING = "utf-8"
Get-Content -Raw section.json | python -m thinwallx inspect -
if ($LASTEXITCODE -ne 0) { throw "ThinWallX inspection failed" }
python -m thinwallx analyze section.json --N 100 --output result.json
if ($LASTEXITCODE -ne 0) { throw "ThinWallX analysis failed" }
```

Use --output for encoding-independent UTF-8 disk publication, especially with legacy Windows shell redirection. No command uploads data or publishes packages.
