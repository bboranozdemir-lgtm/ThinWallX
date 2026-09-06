# ThinWallX 1.0

Python thin-wall centerline section analysis: open, closed and mixed topology,
section properties, shear flow, shear center, torsion, warping and combined stress.
JSON/DXF interchange, headless technical plots and Markdown calculation reports.

## Install and run

Python >=3.10. Install with `python -m pip install .`, or
`python -m pip install ".[plots]"` for optional Matplotlib output.

```sh
thinwallx --version
thinwallx inspect tests/fixtures/dxf/rectangle_r12.dxf
thinwallx analyze tests/fixtures/dxf/rectangle_r12.dxf --N 100 --output result.json
```

`python -m thinwallx` is equivalent to the installed console command.
Existing files are protected by default. JSON uses exact float64 encoding.

## Documentation

- [CLI reference](docs/CLI_REFERENCE.md)
- [Theory and sign conventions](docs/THEORY_AND_CONVENTIONS.md)
- [Verification benchmarks](docs/VERIFICATION_BENCHMARKS.md)
- [Python I/O and plotting guide](docs/V0_8_USAGE.md)

## Verification and limits

Run `python -m pytest -W error`. Packaging tests build wheel/sdist and install
the wheel into a temporary environment. No package publication is performed.
This is a linear-elastic thin-wall model, not a buckling, plasticity or design-code
verification tool. A reported elastic yield multiplier is not safety certification.
