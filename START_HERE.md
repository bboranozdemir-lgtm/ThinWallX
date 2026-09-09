# Sectalix — Developer Guide

This guide describes the local development setup and repository layout for Sectalix.

## Setup

1. Clone the repository:

   ```bash
   git clone https://github.com/berkeboranozdemir/Sectalix.git
   cd Sectalix
   ```

2. Create and activate a virtual environment:

   ```bash
   python -m venv .venv
   ```

   Windows:

   ```bash
   .venv\Scripts\activate
   ```

   Linux/macOS:

   ```bash
   source .venv/bin/activate
   ```

3. Install the package with test and plotting dependencies:

   ```bash
   python -m pip install -e ".[test,plots]"
   python -m pip install "setuptools>=77.0.3" wheel
   ```

4. Run the test suite:

   ```bash
   python -m pytest -W error
   ```

## Repository Structure

- `src/sectalix/` — numerical models, topology handling, stress recovery, I/O, CLI, and reporting
- `tests/` — unit, benchmark, invariance, validation, packaging, and regression tests
- `docs/THEORY_AND_CONVENTIONS.md` — mathematical model and sign conventions
- `docs/VERIFICATION_BENCHMARKS.md` — benchmark definitions and verification notes
- `docs/CLI_REFERENCE.md` — command-line interface reference
- `docs/V0_8_USAGE.md` — interchange-format details retained for compatibility
- `schemas/` — JSON schema definitions
- `examples/` — sample DXF sections and generated calculation outputs

## Development Principles

Changes to the numerical mechanics should be accompanied by a clear derivation or reference and an independent verification case whenever practical.

Please keep the following points in mind:

- preserve explicit units and sign conventions
- avoid silent numerical fallbacks
- add regression tests for corrected numerical behavior
- distinguish exact straight-segment integration from thin-wall engineering approximations
- document new modeling assumptions and limitations
- keep optional visualization separate from the numerical core

## Continuous Integration

The GitHub Actions test matrix runs on Ubuntu, Windows, and macOS using Python 3.10 through 3.13. Warnings are treated as errors during the test run.

Cross-platform tests reduce the risk of implementation regressions, but they do not replace engineering validation against independent reference problems.
