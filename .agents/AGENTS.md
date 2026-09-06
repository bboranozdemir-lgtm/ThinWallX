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
- Implement the exact closed-form straight-segment formulas specified in `ACTIVE_PHASE.md`.
- Do not replace exact segment formulas with numerical sampling or quadrature in v0.1.
- Record the source, derivation, or benchmark basis for every nontrivial formula.
- Strictly adhere to the sign conventions and inertia-matrix structure in `ACTIVE_PHASE.md`.
- Use the quadrant-safe `atan2` principal-angle convention specified in `ACTIVE_PHASE.md`.
- Keep dependencies minimal: use `numpy` for linear algebra/eigenvalues and `pytest` for testing.
- Add tests with every implemented calculation.
- Prefer analytical benchmarks when available.
- Numerical tolerances in tests must be justified.
- Do not hide invalid geometry, zero-length segments, non-positive thicknesses, disconnected geometry, closed loops, non-finite values, or unsupported cases.
- Validation failures must be explicit and testable.

## Architecture Guidance
Keep the v0.1 architecture minimal.
A reasonable decomposition is:
- geometric primitives
- section/topology container
- validation
- section-property calculations
- tests

Do not create abstractions for future shear flow, torsion, warping, closed-cell, FEM, DXF, GUI, or optimization features during v0.1.

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
