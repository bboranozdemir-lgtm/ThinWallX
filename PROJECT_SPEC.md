# ThinWallX — Project Specification

## Purpose
ThinWallX is a numerical analysis tool for arbitrary thin-walled structural sections.

The project focuses on engineering quantities that are meaningful for thin-walled steel and box-type sections, while keeping the implementation transparent, testable, and mathematically traceable.

## Implementation Stack
- **Language:** Python (>= 3.10)
- **Core Numerical Dependencies:** `numpy`
- **Testing Framework:** `pytest`
- **Visualization:** matplotlib (minimal, technical plots and geometry verification, zero interactive/GUI runtimes)
- **Constraint:** Zero heavy symbolic (e.g., `sympy`) or CAD/mesh runtime dependencies.

## Core Principles
- Numerical correctness is more important than feature count.
- Do not implement features before their roadmap phase becomes active.
- Every implemented engineering formula must have a documented source or derivation note.
- Every major numerical result must have at least one independent benchmark.
- No silent numerical fallbacks.
- Units and sign conventions must be explicit.
- Prefer small, testable modules over large abstractions.
- Do not add GUI, web, mobile, AI/chatbot, optimization, or unrelated features.

## Modeling Assumptions
- Pure thin-wall centerline integration ($t \ll L$).
- Corner overlaps at junctions and root radii are explicitly neglected in the centerline formulation.
- Through-thickness local contributions such as the exact rectangular-strip $t^3/12$ term are neglected unless a later phase explicitly adds them.
- v0.1 is limited to connected open-section topologies. Closed loops are intentionally unsupported until the roadmap explicitly activates closed-section analysis.

## Coordinate Convention
Use a right-handed 2D section coordinate system:
- $x$: horizontal
- $y$: vertical
- Positive moments/products follow standard structural mechanics conventions documented in the active phase.
- All geometry calculations must be invariant to segment input order and segment endpoint direction where physically appropriate.

## Units
The computational core is unit-agnostic but requires internally consistent units.
Examples and tests shall use mm, mm², mm⁴ unless explicitly stated otherwise.

## Planned Final Scope
By v1.0, the project should support:
- arbitrary thin-walled centerline geometry
- open and closed sections
- multi-cell sections
- section area and centroid
- Ix, Iy, Ixy
- principal axes/properties
- shear-flow distributions
- shear-center location
- Saint-Venant torsion constant J
- sectorial coordinate / warping function
- warping constant Cw
- technical plots
- JSON input
- optional DXF import
- benchmark and regression test suite

## Out of Scope Until Explicitly Added
- full finite-element structural analysis
- nonlinear material behavior
- buckling analysis
- composite materials
- 3D solid/shell meshing
- GUI-first development
