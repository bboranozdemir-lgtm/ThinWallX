# Verification benchmarks

## Verification philosophy

Sectalix is verified with a combination of analytical benchmark cases, independently assembled numerical references, invariance tests, topology and invalid-input tests, and floating-point stress cases.

The test suite is evidence for the documented model and tested cases. It is not a proof for every possible geometry, and it is not an external engineering certification.

The main theoretical bases are:

- closed-form integration along straight thin-wall centerline segments
- equilibrium of open-section shear flow
- Bredt-Batho closed-cell compatibility and torsion
- Vlasov-type sectorial-coordinate and warping relations

The governing conventions are summarized in [Theory and Conventions](THEORY_AND_CONVENTIONS.md).

## Basic section properties

For a thin-walled L-section with horizontal leg `b`, vertical leg `h`, and uniform thickness `t`, joined at the origin,

$$A=t(b+h),$$

$$c_x=\frac{b^2}{2(b+h)},\qquad
c_y=\frac{h^2}{2(b+h)},$$

$$I_x=\frac{th^3}{3}-Ac_y^2,$$

$$I_y=\frac{tb^3}{3}-Ac_x^2,$$

$$I_{xy}=-Ac_xc_y.$$

The raw product moment is zero because each leg lies on a coordinate axis. This case checks centroid shifting and the sign of `Ixy` independently of the production property routine.

Additional symmetric sections verify zero centroid offsets, zero product moment where required, principal moments, and principal-axis reconstruction.

Rigid translations, rotations, segment-order changes, and segment reversal are used as metamorphic checks where the underlying physical quantity should be invariant or transform predictably.

## Open-section shear flow

Open tree benchmarks check:

- zero shear flow at free edges
- signed equilibrium at junctions
- recovery of the applied transverse resultant
- the quadratic straight-segment shear-flow field and its analytical antiderivative
- consistency under segment reversal and section transformations

Unsymmetric sections exercise both components of the section inertia matrix used in the shear-flow solution.

## Shear center

Simple analytical cases are used to verify shear-center sign conventions and offsets.

For a symmetric channel with flange width `b`, web height `h`, and uniform thickness, the classical thin-wall offset from the web is

$$e=\frac{3b^2}{6b+h}$$

on the side opposite the flanges under the documented orientation convention.

L-, C-, U-, and related orientation variants are used to check coordinate and sign handling.

## Open-section torsion and warping

For a thin-walled open section, the implemented Saint-Venant thin-strip approximation is

$$J=\frac13\sum_i L_it_i^3.$$

For a symmetric thin-wall I-section with flange width `b`, web height `h`, and uniform thickness,

$$J=\frac{(2b+h)t^3}{3},$$

and the sectorial-coordinate benchmark gives

$$C_w=\frac{tb^3h^2}{24}.$$

Warping tests also check:

- zero area-weighted mean after normalization
- consistency of the sectorial field under root selection
- segment-order and reversal invariance
- expected symmetry properties

## Closed and multi-cell sections

For a rectangular centerline box of width `b`, height `h`, and uniform thickness `t`, the Bredt-Batho torsion constant is

$$J_{BB}=\frac{4b^2h^2}{2(b+h)/t}
=\frac{2tb^2h^2}{b+h}.$$

The corresponding thin-wall second moment about the horizontal centroidal axis is

$$I_x=th^2\left(\frac b2+\frac h6\right).$$

Closed-section tests verify:

- bounded-cell extraction and signed cell areas
- cell-edge incidence matrices
- Bredt-Batho compatibility
- independence of the final physical shear flow from the chosen virtual cut
- multi-cell shared-wall coupling
- concave bounded faces
- scale-aware solution of very thin-wall systems

Small multi-cell systems are compared against independently assembled analytical or high-precision linear systems.

## Mixed open/closed sections

Mixed-section benchmarks include sections with closed cells connected by open bridges and sections with closed cells plus open branches.

The tests verify:

- bridge and cyclic-region decomposition
- correct inclusion of open branches in

  $$J_{open}=\frac13\sum L_et_e^3$$

- exclusion of closed walls from the open-strip contribution
- compatible shear-flow and sectorial-coordinate fields across junctions
- independent behavior of disconnected cyclic blocks coupled only through the open graph
- invariance under geometry and segment-definition changes

## Stress recovery

Independent normal-stress checks use axial stress and unsymmetric bending relations under the documented sign convention.

The stress tests cover:

- recovery of `N`, `Mx`, `My`, and warping bimoment resultants from the computed stress field
- transverse and torsional shear-stress contributions
- combined loading
- zero-load behavior
- first-yield multiplier behavior when a positive yield stress is supplied

The peak von Mises search is verified with cases whose extrema occur at endpoints, interior stationary points, and membrane shear sign changes. Plot samples are not used as the reference for the reported maximum.

## Interchange, DXF, plotting, and packaging

I/O tests cover:

- JSON round trips and malformed-schema rejection
- exact preservation of supported float64 values in the hexadecimal interchange representation
- representative open, closed, and mixed DXF fixtures
- entity ordering and coordinate/unit metadata handling
- PNG and SVG output generation
- command-line exit behavior
- calculation-report generation
- wheel and source-distribution packaging
- installed CLI invocation from a temporary environment

Visualization tests verify production and deterministic numerical inputs; pixel-identical rendering across all font or Matplotlib versions is not treated as a mechanical requirement.

## Numerical robustness

The suite includes cases designed to expose floating-point problems such as:

- large rigid translations
- cancellation in centroidal moment calculations
- near-overflow and near-underflow intermediate values
- positive subnormal results
- ill-conditioned section matrices
- scale differences between independent closed-cell blocks

Where appropriate, expected values are assembled with higher-precision arithmetic, exponent-aware calculations, or formulas that are independent of the production implementation.

A finite result alone is not considered sufficient evidence if the governing system is too ill conditioned for the documented numerical contract.

## Continuous integration

The repository CI runs the test suite on:

- Ubuntu
- Windows
- macOS

with Python 3.10, 3.11, 3.12, and 3.13. Warnings are treated as errors.

Cross-platform CI is intended to catch implementation and packaging regressions. Engineering validation for a practical application should additionally compare representative sections against trusted published references, measurements where available, or independent established analysis software.
