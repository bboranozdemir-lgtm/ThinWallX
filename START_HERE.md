# ThinWallX — Developer & Contributor Guide

Welcome to ThinWallX, a production-quality numerical analysis library for arbitrary thin-walled structural cross-sections.

## Project Status

ThinWallX has completed its initial v0.1–v1.0 roadmap. The numerical core (v0.1–v0.8) and CLI/reporting modules (v1.0) are formally **FROZEN** with 559 regression tests passing at 100% with zero warnings (`-W error`).

## Quick Developer Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/bboranozdemir-lgtm/ThinWallX.git
   cd ThinWallX
   ```

2. Create and activate a virtual environment:
   ```bash
   python -m venv .venv
   # Windows:
   .venv\Scripts\activate
   # Linux/macOS:
   source .venv/bin/activate
   ```

3. Install in editable development mode with test and plotting extras:
   ```bash
   python -m pip install -e ".[test,plots]"
   python -m pip install "setuptools>=77.0.3" wheel
   ```

4. Run the full test suite with strict error checking:
   ```bash
   python -m pytest -W error
   ```

## Repository Structure

- `src/thinwallx/`: Core numerical analysis library, CLI, and reporting modules.
- `tests/`: 559 comprehensive unit, benchmark, and regression tests.
- `docs/`: Technical reference documentation:
  - `docs/THEORY_AND_CONVENTIONS.md`: Formulations, sign conventions, and coordinate rules.
  - `docs/CLI_REFERENCE.md`: CLI command options, subcommands, and UNIX piping.
  - `docs/VERIFICATION_BENCHMARKS.md`: Analytical and independent oracle verification.
  - `docs/archive/`: Historical phase specifications and implementation audit reports.
- `examples/`: Sample section models, analysis outputs, and scripts.

## Contribution Guidelines

1. **Frozen Mechanics:** Core analytical models must not be altered without explicit validation against analytical benchmarks.
2. **Strict Quality Gates:** All pull requests must pass the complete 12-job cross-platform CI matrix (`pytest -W error` on Ubuntu, Windows, macOS across Python 3.10–3.13).
3. **No Heavy Dependencies:** Runtime dependencies remain strictly minimal (`numpy` and standard library, optional `matplotlib` for plotting). No heavy CAD, mesh, or symbolic algebra packages.
