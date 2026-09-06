# Verification benchmarks

## Evidence and independence

v1.0 preserves the previous 520 tests and adds CLI, reporting and real packaging checks. Passing tests are regression evidence, not a proof for all finite geometries. The mathematical bases below are closed-form segment integration, equilibrium, Bredt-Batho compatibility and Vlasov warping. A production-derived expected value verifies transport/presentation only; it is not labelled an independent mechanical oracle.

The phase specifications record the frozen conventions. [Theory](THEORY_AND_CONVENTIONS.md) collects their equations and [CLI reference](CLI_REFERENCE.md) explains observable errors. Classical background is Vlasov's *Thin-Walled Elastic Beams* and the Bredt-Batho thin-wall closed-cell relation. Those names identify the theoretical framework; no unverified page quotation or external numerical dataset is claimed.

## v0.1: exact line-area properties

Straight-strip integrals are derived by expanding linear endpoint interpolation. For an L with horizontal leg b and vertical leg h, uniform t, joined at the origin:

$$A=t(b+h),\quad c_x=\frac{b^2}{2(b+h)},\quad c_y=\frac{h^2}{2(b+h)},$$
$$I_x=\frac{th^3}{3}-Ac_y^2,\quad I_y=\frac{tb^3}{3}-Ac_x^2,\quad I_{xy}=-Ac_xc_y.$$

The raw product moment is zero because each leg lies on an axis. This independently checks sign and centroid shifting. Symmetric I shapes check zero centroid offsets and product moment. Principal eigenvalues and reconstructed tensors test quadrant handling, rotations and isotropy. Translation and segment reversal are metamorphic checks, not new analytical references.

## v0.2–v0.4: open flows, center and warping

Tree benchmarks check free-edge q=0, junction balance and the exact quadratic flow antiderivative. Unsymmetric shapes exercise both components of C alpha=V. Reversal tests compare physical vectors rather than blindly equating directed scalar q.

Intersecting straight-leg L/T profiles have their idealized shear center at the intersection. For a symmetric channel of flange width b, web height h and uniform thickness, the external offset from the web is 3b²/(6b+h), on the side opposite the flanges. This must not be confused with the centroid-relative offset. C and U orientation variants test signs.

For a symmetric thin-wall I section, flange width b, web height h and uniform t, J=(2b+h)t³/3 and Cw=t b³ h²/24. Signed flange warping follows the pole convention, not merely its squared integral. Independent tests additionally integrate omega, omega X and omega Y. Root choice, edge order and reversed definitions test graph integration.

Extreme-scale benchmarks address intermediate products, cancellation, near-DBL_MAX values and positive subnormal totals. Expected values use exponent arithmetic or high-precision arithmetic, rather than reusing the production integrator.

## v0.5: closed cells

A rectangular centerline box b by h with uniform t has

$$J_{BB}=\frac{4b^2h^2}{2(b+h)/t}
=\frac{2tb^2h^2}{b+h},\qquad I_x=th^2(b/2+h/6).$$

These are independent of any Section property call. For the specified reference cut at xi_c in the transverse oracle, q0=-(Vy t b h(1-2xi_c))/(4Ix). q0 is cut-dependent; total physical wall flow is not. Tests vary xi_c, including the zero-circulation symmetry case, rather than treating one arbitrary cut as universal.

Multi-cell tests assemble shared-wall compatibility and compare small independent analytical/high-precision systems. Concave bounded faces test signed area extraction. Very thin walls test scaled H reconstruction and solve behavior separately: a finite solution does not guarantee a representable physical H accessor.

## v0.6: mixed graphs

Barbell sections have distinct cell blocks joined by a bridge. Tests verify zero inter-block coupling, bridge inclusion in the spanning tree, continuous sectorial fields and global area-weighted normalization. Hat profiles exercise coupled closed walls and open branches. Jopen includes bridges only; adding closed wall strip torsion would double count the frozen model.

Heterogeneous block sizes and thicknesses test independent exponent scales. Extreme open contributions, subnormal Cw, root/reversal invariance, invalid user trees and explicit range exceptions are regression cases. They deliberately distinguish exact zero, small numerical residuals and unrepresentable nonzero values.

## v0.7: stresses

Independent normal-stress formulas use N/A and the unsymmetric bending inverse in the documented sign convention. Decimal-based benchmarks validate analytical profiles and recovered N/Mx/My/B. Open-strip surface torsion is distinguished from closed-cell membrane torsion. Combined load cases check the surface envelope rather than adding magnitudes of unrelated normal stresses.

Quartic peak tests include endpoints, interior stationary points and membrane sign changes. A dense plot cannot substitute for these roots. Zero loads with yield stress yield an unbounded elastic multiplier. Nonzero B or M_omega on unsupported zero-Cw states must fail explicitly.

## v0.8: interchange and visualization

JSON round trips compare float64 representations and object fields, including extreme scales, signed zero and malformed-schema rejection. DXF fixtures cover open L, closed rectangle and mixed barbell imports, entity order and supported coordinate/unit transformations. Mechanical import checks use analytical references, while codec round trips check preservation.

Headless PNG/SVG tests check file production, deterministic behavior in a fixed rendering environment, no overwrite and cleanup. Pixels and font layout need not be identical across different Matplotlib/font versions; mathematical results do not depend on rendering samples.

## v1.0 traceability

| Requirement | Test evidence |
|---|---|
| T01–T20 | test_cli.py: subprocess help, codes, inspect, loads, plots, report, pipes, optional dependency and no-clobber |
| T21–T30 | test_reporting.py: report fields, independent L properties, relative links, atomic failure, open/closed/mixed |
| T31–T36 | test_packaging.py: version, entry point, markers, documents, real wheel/sdist and installed invocation |
| T37–T40 | Entire unchanged v0.1–v0.8 test suite |

A temporary environment installs the built wheel with no dependency downloads and checks that imports resolve there, not to the checkout. Existing NumPy is shared for this test; this is not an offline wheel bundle or a test of every supported Python/OS combination.

## Tolerances and limitations

Dimensionless invariance/compatibility criteria typically use 1e-10; analytical ordinary-scale comparisons often use 1e-12 or tighter. Individual frozen tests state their own justified limits. Tiny nonzero oracles require zero absolute tolerance; default np.isclose atol can otherwise accept zero incorrectly. True zero checks require a dimensional local scale or an exact representation comparison as appropriate.

Subnormal results cannot promise normal-range relative accuracy: representable spacing is approximately 4.94e-324. Large rigid translations can destroy input detail before the solver sees it. Ill-conditioning is rejected at the frozen eigenvalue threshold rather than masked by arbitrary regularization. No finite test suite establishes correctness outside the documented centerline assumptions.
