# AGENTS.md — Codex Instructions

## Required Reading Order
Before making changes, read only:
1. `PROJECT_SPEC.md`
2. `ACTIVE_PHASE.md`
3. the relevant existing source/test files for the active phase

Read `ROADMAP.md` only when the active phase requires context about planned interfaces.
Do not repeatedly re-read unrelated documentation.

## Scope Control
- Implement ONLY the version named in `ACTIVE_PHASE.md`.
- Never implement future-roadmap features proactively.
- Never advance `ACTIVE_PHASE.md`.
- Never reinterpret scope to include "helpful" extra features.
- If a requested change conflicts with `PROJECT_SPEC.md` or `ACTIVE_PHASE.md`, stop and report the conflict.

## Engineering & Python Implementation Rules
- Use clean, modern Python (>= 3.10) with explicit type hinting.
- Use lightweight `dataclasses` for geometric primitives such as `Node` and `Segment`.
- Do not invent engineering formulas.
- Implement the exact closed-form straight-segment formulas specified in the specifications.
- Do not replace exact segment formulas with numerical sampling or quadrature.
- Record the source, derivation, or benchmark basis for every nontrivial formula.
- Strictly adhere to the sign conventions and inertia-matrix structure.
- Use the quadrant-safe `atan2` principal-angle convention.
- Keep dependencies minimal: use `numpy` for linear algebra/eigenvalues and `pytest` for testing.
- Add tests with every implemented calculation.
- Prefer analytical benchmarks when available.
- Numerical tolerances in tests must be justified.
- Do not hide invalid geometry, zero-length segments, non-positive thicknesses, disconnected geometry, closed loops, non-finite values, or unsupported cases.
- Validation failures must be explicit and testable.

## Architecture Guidance
The v1.0 production architecture is complete and frozen:
- geometric primitives and section containers (`primitives.py`, `section.py`, `closed_section.py`, `mixed_section.py`)
- cell detection and graph topology (`cells.py`, `mixed_topology.py`, `validation.py`)
- exact analytical mechanics (`properties.py`, `shear_flow.py`, `shear_center.py`, `torsion.py`, `closed_shear_flow.py`, `closed_torsion.py`, `mixed_shear_flow.py`, `mixed_torsion.py`, `stress.py`)
- I/O and visualization (`serialization.py`, `dxf.py`, `plotting.py`)
- user interfaces and reporting (`cli.py`, `reporting.py`)
- regression tests (`tests/`, 559 tests)

## Token / Context Efficiency
- Inspect only files relevant to the active task.
- Do not summarize the whole repository unless explicitly asked.
- Do not refactor unrelated files.
- Do not generate large design documents unless explicitly requested.
- Do not duplicate information already present in project markdown files.
- Make minimal targeted edits.
- Prefer running targeted tests first, then the full test suite at phase completion.

## Phase Completion Procedure
When all `ACTIVE_PHASE.md` acceptance criteria are satisfied:
1. run the required targeted tests,
2. run the full test suite with `pytest`,
3. report:
   - files changed
   - formulas/methods implemented
   - benchmark results
   - validation behavior
   - known mathematical/numerical risks
   - remaining items, if any
4. STOP.

Do not start the next roadmap version.
Only the user may authorize a new phase.
