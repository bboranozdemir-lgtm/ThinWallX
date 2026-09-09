# Sectalix

[![CI](https://github.com/berkeboranozdemir/Sectalix/actions/workflows/ci.yml/badge.svg)](https://github.com/berkeboranozdemir/Sectalix/actions/workflows/ci.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Sectalix** is a Python library and command-line tool for the analysis of thin-walled structural cross-sections represented by piecewise-straight centerline segments with constant thickness on each segment.

The library combines closed-form straight-segment integration with classical thin-wall formulations for open, closed multi-cell, and mixed open/closed topologies. The aim is to keep the mathematical model transparent, testable, and useful for engineering studies without requiring finite-element meshing for section-level calculations.

![Sectalix Von Mises Stress Analysis](examples/sample_calculation/stress_vm.png)

---

## Capabilities

- **Closed-form segment integration:** Area, centroid, second moments of area (`Ix`, `Iy`, `Ixy`), principal moments, and principal-axis orientation are evaluated directly from the straight centerline segments.
- **Open-section shear flow and shear center:** Transverse shear flow is obtained from equilibrium on open tree topologies, followed by shear-center evaluation from the resulting torque.
- **Closed and multi-cell sections:** Closed-cell compatibility and torsional response are handled with the Bredt-Batho thin-wall formulation.
- **Mixed open/closed topology:** Sections containing both closed cells and open branches are decomposed and analyzed within the documented thin-wall model.
- **Warping quantities:** Sectorial coordinates and the warping constant `Cw` are evaluated using Vlasov-type thin-wall kinematics and the conventions documented in the theory notes.
- **Stress recovery:** Linear-elastic normal and shear stresses can be recovered for combined section resultants, including analytical search for the peak von Mises stress along each straight segment.
- **Engineering interface:** A command-line interface supports section inspection, analysis, plotting, calculation reports, JSON interchange, and ASCII DXF import.
- **Minimal runtime dependencies:** Core calculations use the Python standard library and `numpy`; `matplotlib` is optional for plotting.

---

## Installation

### From PyPI

```bash
python -m pip install sectalix
```

### From source

```bash
git clone https://github.com/berkeboranozdemir/Sectalix.git
cd Sectalix
pip install -e .
```

For plotting:

```bash
pip install -e ".[plots]"
```

---

## Quickstart

### Command line

```bash
# Inspect section properties
sectalix inspect examples/sample_sections/rectangle.dxf

# Recover stresses for an applied load set
sectalix analyze examples/sample_sections/rectangle.dxf \
  --N 10000 --Vy 5000 --Mx 250000 --sigma-yield 355 \
  --report calculation_report.md \
  --plots-dir ./output_plots

# Plot the imported centerline geometry
sectalix plot examples/sample_sections/rectangle.dxf --output geometry.png
```

The rectangle example is intentionally small and is included to demonstrate the interface rather than represent a practical design case. Plot-producing commands require the optional `plots` dependency.

### Python API

```python
from sectalix import Node, Segment, Section, AppliedLoads

# Open channel section: 50 x 100 x 50 mm, t = 2 mm
n0 = Node(x=50.0, y=100.0)
n1 = Node(x=0.0, y=100.0)
n2 = Node(x=0.0, y=0.0)
n3 = Node(x=50.0, y=0.0)

section = Section([
    Segment(p1=n0, p2=n1, t=2.0),
    Segment(p1=n1, p2=n2, t=2.0),
    Segment(p1=n2, p2=n3, t=2.0),
])

print(f"Area: A = {section.area:.2f} mm^2")
print(f"Centroid: C = ({section.cx:.2f}, {section.cy:.2f}) mm")
print(f"Inertia: Ix = {section.Ix:.2f}, Iy = {section.Iy:.2f} mm^4")
print(f"Saint-Venant torsion constant: J = {section.J:.2f} mm^4")
print(f"Warping constant: Cw = {section.Cw:.2f} mm^6")

loads = AppliedLoads(
    N=10000.0,
    Vy=5000.0,
    Mx=250000.0,
    sigma_yield=355.0,
)
results = section.calculate_stresses(loads)

print(f"Peak von Mises stress: {results.max_sigma_vm:.2f} MPa")
print(f"Elastic first-yield load factor: {results.load_factor:.3f}")
```

---

## Supported section models

| Section model | Shear flow | Torsion | Warping | Stress recovery |
|---|---|---|---|---|
| **Open branched sections** (I, C, L, T, Z, hat) | Open-tree equilibrium | Thin-strip Saint-Venant approximation | Sectorial-coordinate formulation | Yes |
| **Closed multi-cell sections** (boxes, tubes, box girders) | Closed-cell compatibility | Bredt-Batho | Compatible sectorial field | Yes |
| **Mixed sections** (closed cells with open branches) | Combined open/closed formulation | Closed-cell + open-branch model | Continuous sectorial field | Yes |

---

## Model assumptions and limitations

Sectalix is a **thin-wall centerline model**, not a general solid-section or shell finite-element solver.

- Geometry is represented by straight wall-centerline segments.
- Each segment has a constant, positive thickness.
- The formulation assumes a homogeneous, linear-elastic thin wall.
- Corner fillets, root radii, local corner overlap volumes, and detailed through-thickness geometry are not represented.
- The implemented torsion, shear-flow, and warping relations are thin-wall engineering formulations; they are not full finite-thickness elasticity solutions.
- Local or global buckling, nonlinear material response, fatigue, connection behavior, and structural design-code checks are outside the scope of the library.
- The reported elastic load factor is a first-yield indicator, not a code-based safety factor or design certification.

Within these assumptions, the geometry may be non-standard, asymmetric, branched, multi-cell, or mixed, provided that it can be represented by a valid piecewise-straight centerline network supported by the relevant section class.

---

## Verification

The numerical implementation is checked with analytical benchmark cases, high-precision or independently assembled reference calculations where appropriate, rigid-body and orientation invariance tests, adversarial numerical cases, and regression tests.

Continuous integration runs the test suite on Ubuntu, Windows, and macOS with Python 3.10 through 3.13. Passing tests provide evidence for the documented model and cases; they are not a proof of correctness for every possible geometry or engineering application.

See **[Verification Benchmarks](docs/VERIFICATION_BENCHMARKS.md)** for the benchmark definitions and limitations.

---

## Documentation

- **[Theory and Conventions](docs/THEORY_AND_CONVENTIONS.md)** — mathematical model, coordinate system, signs, and governing equations.
- **[Verification Benchmarks](docs/VERIFICATION_BENCHMARKS.md)** — analytical and independent verification cases.
- **[CLI Reference](docs/CLI_REFERENCE.md)** — command-line usage, inputs, outputs, and exit codes.
- **[Developer Guide](START_HERE.md)** — local setup, repository structure, and testing.
- **[Project Scope](PROJECT_SPEC.md)** — implemented scope and explicit limitations.

---

## License

Sectalix is released under the [MIT License](LICENSE).
