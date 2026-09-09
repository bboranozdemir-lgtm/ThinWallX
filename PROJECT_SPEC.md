# Sectalix — Project Scope

## Purpose

Sectalix is a numerical analysis library for thin-walled structural cross-sections represented by piecewise-straight centerline segments.

The project focuses on section-level quantities that are useful in structural and mechanical engineering studies while keeping the mathematical assumptions, sign conventions, and numerical limitations explicit.

## Computational Model

- **Language:** Python 3.10+
- **Core numerical dependency:** `numpy`
- **Testing:** `pytest`
- **Optional plotting:** `matplotlib`
- **Geometry model:** straight centerline segments with constant thickness per segment
- **Material model:** homogeneous, linear-elastic thin wall

The computational core is unit-agnostic. Inputs must use a consistent unit system.

## Implemented Capabilities

Sectalix supports the following section-level calculations within the documented thin-wall assumptions:

- area and centroid
- `Ix`, `Iy`, and `Ixy`
- principal moments and principal-axis orientation
- open-section transverse shear flow
- shear-center location
- open-section Saint-Venant thin-strip torsion constant
- closed and multi-cell Bredt-Batho torsion
- sectorial coordinates and warping constant `Cw`
- mixed open/closed topologies
- linear-elastic stress recovery under combined section resultants
- analytical peak search for von Mises stress along straight segments
- JSON interchange
- ASCII DXF import for supported entities
- technical plotting and Markdown calculation reports
- command-line and Python interfaces

## Geometry Scope

The geometry may be non-standard, asymmetric, branched, closed, multi-cell, or mixed, provided that it can be represented by a valid piecewise-straight thin-wall centerline network supported by the relevant solver.

This does **not** mean that Sectalix represents an arbitrary continuum cross-section. Curved walls must first be represented by straight segments, and finite-width solid geometry is not modeled directly.

## Main Assumptions

- Thin-wall centerline idealization is applicable to the section being studied.
- Segment thickness is positive and constant along each segment.
- Corner fillets, root radii, and local overlap volumes at wall intersections are neglected.
- The material response is linear elastic.
- The section topology must satisfy the validity conditions documented for the selected solver.
- Numerical calculations use IEEE-754 double precision and explicitly reject unsupported singular or severely ill-conditioned cases where required.

## Outside Scope

Sectalix does not perform:

- full 2D or 3D finite-element analysis
- shell or solid meshing
- nonlinear material analysis
- local or global buckling checks
- plastic-section or plastic-hinge analysis
- fatigue or fracture assessment
- connection design
- structural-system analysis
- code-based member design or safety certification

## Verification Approach

The implementation is checked using a combination of:

- independently derived analytical benchmark cases
- high-precision or independently assembled numerical references where appropriate
- rigid translation, rotation, segment-order, and orientation invariance tests
- topology and invalid-input tests
- extreme-scale and floating-point edge cases
- cross-platform continuous integration

The verification suite supports confidence in the documented model and tested cases, but it is not a mathematical proof for every possible geometry or an external engineering certification.

## Engineering Direction

Sectalix is intended as a transparent section-analysis framework that can support research, education, and engineering studies involving non-standard thin-walled profiles. Further engineering validation should compare representative practical sections against trusted analytical references, published data, or independent established software before use in safety-critical design workflows.
