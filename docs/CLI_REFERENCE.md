# Sectalix 1.0.1 CLI reference

## Installation

Python 3.10 or newer is required. Install the release with:

```bash
python -m pip install sectalix
```

From a source checkout:

```bash
python -m pip install -e .
```

For plotting support:

```bash
python -m pip install "sectalix[plots]"
```

NumPy is the only required runtime dependency. Matplotlib is optional and loaded only when plotting is requested. Both `sectalix` and `python -m sectalix` use the same command-line entry point.

`sectalix --version` prints the installed Sectalix version. `--help` and each subcommand's `--help` show the accepted options.

## Input and output rules

Section input is either a supported v0.8 JSON document or a supported ASCII DXF file, selected by filename extension. A positional `-` means UTF-8 JSON from standard input, not DXF.

See the [v0.8 interchange guide](V0_8_USAGE.md) for the JSON schema, supported DXF entities, and thickness mapping.

Coordinates and loads must use a consistent unit system. JSON unit metadata is descriptive and does not request automatic numerical conversion.

The section-consuming commands support `--output PATH`, `--overwrite`, and `--thickness VALUE` where applicable. Plot and DXF-conversion commands require an explicit output target.

Thickness is a fallback for DXF entities that do not receive thickness from a more specific mapping. It must be positive. Standard section-consuming commands accept `--dxf-unit UNIT`; `convert-dxf` provides separate source and target unit options.

Existing targets are protected unless `--overwrite` is explicitly supplied. Outputs may not overwrite the input section or load file. Duplicate output paths are rejected before writing.

Individual output files are written atomically. A multi-file analysis is not a single transaction, so plots already written may remain if a later requested output fails.

Standard output contains requested data only. Diagnostics are written to standard error. Numeric options accept scientific notation, including negative values such as `--Mx -1e-2`. Non-finite or out-of-range inputs are rejected.

## `inspect`

```text
sectalix inspect INPUT [--json] [--output PATH] [--dxf-unit mm] [--thickness VALUE] [--overwrite]
```

The default human-readable output reports:

- topology
- node and segment counts
- area
- centroid
- `Ix`, `Iy`, `Ixy`
- principal moments and principal-axis angle
- torsion constant `J`
- shear-center coordinates
- warping constant `Cw`
- total centerline length
- minimum and maximum segment thickness

Without `--output`, the result is written to standard output. `--output -` also selects standard output.

`--json` produces a v1.0 inspection envelope. The envelope includes encoded property values and embeds the original v0.8 section document under `section`. The whole inspection envelope is not itself accepted as a section input; extract its `section` member for reuse.

Inspection requests the complete listed property set. If a requested derived quantity is singular or numerically unsupported, inspection fails explicitly instead of replacing that value with `null`.

## `analyze`

```text
sectalix analyze INPUT [--loads LOADS.json] [load flags] [--output RESULT.json] [--stdout] [--plots-dir DIR] [--report REPORT.md] [--dxf-unit mm] [--thickness VALUE] [--overwrite]
```

Load options are:

- `--N`
- `--Vx`
- `--Vy`
- `--Mx`
- `--My`
- `--Tsv`
- `--B`
- `--M-omega`
- `--sigma-yield`

An explicit command-line load value overrides only the corresponding field from a load document. Without a load file, unspecified force and moment fields default to zero and yield stress remains unspecified.

At least a load file or one explicit load option must be supplied. An explicitly requested zero load, such as `--N 0`, is valid.

A load document must decode as `AppliedLoads`. Its length metadata must be compatible with the section metadata, and force metadata must also be compatible when both are specified. No automatic numerical rescaling is performed.

Section input and load input cannot both use standard input simultaneously. `--loads -` is supported when the section comes from a file.

The primary result is a v0.8 `StressRecoveryResult` document containing the analytical segment profiles and recovered resultant fields. Without `--output`, JSON is written to standard output. `--stdout` additionally echoes JSON when a file output is also requested.

`--plots-dir` writes:

- `geometry.png`
- `shear_flow.png`
- `stress_vm.png`

The shear-flow image represents transverse `Vx/Vy` flow. The stress image uses the full recovered stress result.

`--report` writes a Markdown calculation sheet with section properties, loads, peak von Mises stress, elastic first-yield multiplier, resultant recovery information, and relative links to any requested plots.

The elastic first-yield multiplier is not a code-based safety factor or regulatory certification.

## `plot`

```text
sectalix plot INPUT --output FIGURE.png [--show-thickness | --no-show-thickness] [--show-axes] [--show-ids] [--show-shear-center] [--dxf-unit mm] [--thickness VALUE] [--overwrite]
```

PNG and SVG output are supported. Binary image output is not written to standard output.

Thickness bands are enabled by default. IDs can be shown for nodes and segments. Principal directions and the shear center can be requested explicitly.

Plots use equal aspect ratio and a headless rendering path. If Matplotlib is not installed, non-plotting commands remain available while a requested plot exits with an I/O/format error code.

## `convert-dxf`

```text
sectalix convert-dxf INPUT.dxf --output SECTION.json [--source-unit mm] [--target-unit mm] [--thickness VALUE] [--section-kind auto|open|closed|mixed] [--overwrite]
```

The importer performs supported length conversion, centerline validation, and topology classification. Defaults are `mm` for source and target units and `auto` for topology classification.

`--output -` writes JSON to standard output. Input must be a supported ASCII DXF file.

The importer deliberately does not approximate unsupported CAD entities or act as a general CAD kernel. See the [interchange guide](V0_8_USAGE.md) for the supported subset.

From a source checkout, a complete example is:

```bash
sectalix convert-dxf examples/sample_sections/open_l.dxf --output section.json
```

The resulting JSON can be supplied directly to `inspect` or `analyze`.

Canonical v0.8 JSON stores supported float64 values as hexadecimal strings. Replacing those encoded fields with ordinary decimal JSON numbers does not satisfy the v0.8 format. The [JSON schema](../schemas/sectalix-0.8.schema.json) documents the required structure.

## Exit status

| Code | Meaning | Typical cause |
|---|---|---|
| 0 | Success | Help, version, or completed output |
| 1 | Invocation | Unknown option, missing load specification, colliding paths |
| 2 | I/O or format | Missing file, malformed JSON, unit mismatch, missing Matplotlib |
| 3 | Geometry | Invalid thickness or finite geometry |
| 4 | Topology | Disconnected or unsupported graph |
| 5 | Singular mechanics | Singular section quantity for the requested calculation |
| 6 | Numerical range | Overflow or nonzero underflow |

Expected Sectalix exceptions are mapped to these categories. Unexpected programming errors are not silently swallowed or reclassified as valid engineering output.

## Automation examples

Bash:

```bash
set -e
sectalix convert-dxf shape.dxf --output section.json --thickness 0.01
cat section.json | sectalix analyze - --N 100 --Mx -1e-2 > result.json
sectalix analyze section.json --N 100 --plots-dir figures --report report.md --output result2.json
```

PowerShell:

```powershell
$env:PYTHONIOENCODING = "utf-8"
Get-Content -Raw section.json | python -m sectalix inspect -
if ($LASTEXITCODE -ne 0) { throw "Sectalix inspection failed" }
python -m sectalix analyze section.json --N 100 --output result.json
if ($LASTEXITCODE -ne 0) { throw "Sectalix analysis failed" }
```

For reliable UTF-8 file output across shells, prefer `--output` when practical. Sectalix CLI commands do not upload analysis data or publish packages.
