# Theory and conventions

## Model and applicability

Sectalix represents a thin-walled cross-section by straight wall-centerline segments. Each segment has a constant positive thickness. The area measure is therefore `t ds`, rather than the union of finite-width wall polygons.

The model assumes a homogeneous, linear-elastic thin wall. Corner overlaps, fillets, root radii, local plate bending, buckling, yielding redistribution, fatigue, connections, and design-code checks are outside the present formulation. The quantities `J` and `Cw` are geometric section constants within the documented thin-wall models; compatible material properties are required to form stiffnesses such as `GJ` and `E Cw`.

Sectalix uses IEEE-754 double-precision arithmetic. The implementation includes scale-aware calculations and explicit numerical checks, but it does not provide arbitrary-precision arithmetic.

The v0.8 interchange identifier `"format": "thinwallx"` is retained for backward compatibility with files created before the project was renamed to Sectalix.

## Axes, loads, and resultants

The section coordinate system is right-handed: `x` points right, `y` points up, and `z` points out of the section plane. Let

$$X=x-c_x,\qquad Y=y-c_y$$

be centroid-relative coordinates.

Positive `N` is tensile. The normal-stress resultant convention is

$$N=\int_A\sigma\,dA,$$

$$M_x=-\int_A Y\sigma\,dA,$$

$$M_y=\int_A X\sigma\,dA,$$

$$B=\int_A\omega^\ast\sigma\,dA.$$

`Vx` and `Vy` are the integrals of the transverse shear-flow vector components. Positive `Tsv` is torque about `+z`. `M_omega` is the secondary warping moment used by the stress-recovery formulation.

Typical dimensions are:

- `N`, `Vx`, `Vy`: force
- `Mx`, `My`, `Tsv`, `M_omega`: force × length
- `B`: force × length²
- stress: force / length²

Sectalix is unit-agnostic and does not automatically convert values between unit systems.

## Closed-form straight-segment integration

For a straight segment parameterized by

$$\mathbf r(s)=\mathbf r_1+s\mathbf u,\qquad 0\le s\le L,$$

where `u` is the unit tangent, any field that varies linearly between endpoint values can be integrated directly.

For a linear field `f` and two linear fields `f`, `g`,

$$\int_0^L f\,ds=\frac{L}{2}(f_1+f_2),$$

$$\int_0^L fg\,ds=\frac{L}{6}(2f_1g_1+f_1g_2+f_2g_1+2f_2g_2).$$

Multiplication by thickness gives the corresponding thin-wall area integrals. For segment `i`,

$$A_i=t_iL_i.$$

The total area and first moments give

$$A=\sum_i A_i,$$

$$c_x=\frac{1}{A}\sum_i A_i\frac{x_{i1}+x_{i2}}{2},\qquad
c_y=\frac{1}{A}\sum_i A_i\frac{y_{i1}+y_{i2}}{2}.$$

Raw second moments follow from the same straight-segment identities. The centroidal quantities are equivalent to

$$I_x=I_{x,0}-Ac_y^2,$$

$$I_y=I_{y,0}-Ac_x^2,$$

$$I_{xy}=I_{xy,0}-Ac_xc_y.$$

The implementation evaluates centroidal moments using local reference shifts where useful to reduce cancellation under large rigid translations.

## Principal properties

Sectalix uses the inertia tensor

$$\mathbf I=
\begin{pmatrix}
I_x & -I_{xy}\\
-I_{xy} & I_y
\end{pmatrix}.$$

Let

$$d=\operatorname{hypot}(I_x-I_y,2I_{xy}).$$

The ordered principal moments are

$$I_{1,2}=\frac{I_x+I_y\pm d}{2},\qquad I_1\ge I_2.$$

The principal-axis angle associated with `I1` is

$$\theta_p=\frac12\operatorname{atan2}(-2I_{xy},I_x-I_y),$$

mapped to the documented interval. For an isotropic inertia tensor, the principal direction is not unique; the implementation returns zero by convention.

## Open-section transverse shear flow

For an open tree topology, define

$$\mathbf C=
\begin{pmatrix}
I_y&I_{xy}\\
I_{xy}&I_x
\end{pmatrix},\qquad
\mathbf C\boldsymbol\alpha=
\begin{pmatrix}V_x\\V_y\end{pmatrix}.$$

For a directed segment with normalized coordinate `xi = s/L`, the partial first-moment vector is

$$\mathbf m_{\rm partial}(\xi)=tL\left[
\xi(X_1,Y_1)^T+
\frac{\xi^2}{2}(\Delta x,\Delta y)^T
\right].$$

The shear flow along that segment is a quadratic polynomial,

$$q(\xi)=a_0-tL(\alpha_xX_1+\alpha_yY_1)\xi
-\frac{tL}{2}(\alpha_x\Delta x+\alpha_y\Delta y)\xi^2.$$

The integration constant `a0` follows from subtree equilibrium. Free ends have zero flow, and signed flow balance is enforced at junctions. Integrals of `q` and its moments are evaluated from polynomial antiderivatives rather than sampled plot values.

The inertia matrix used in the shear solution must be positive definite and sufficiently well conditioned for the requested calculation.

## Shear center

For a reference point `r_ref`, the torque contribution of a straight segment is obtained from

$$T_{\rm ref}=\int (\mathbf r-\mathbf r_{\rm ref})\times(q\mathbf u)\,ds.$$

Because the cross product factor is constant along a straight segment, the segment torque reduces to that factor multiplied by `∫q ds`.

Evaluating the centroidal torque for unit basis shear loads gives the centroid-relative shear-center offsets. With the implemented sign convention,

$$e_x=T_c(V_x=0,V_y=1),$$

$$e_y=-T_c(V_x=1,V_y=0),$$

and

$$S=(c_x+e_x,c_y+e_y).$$

## Closed cells and Bredt-Batho torsion

For a closed or multi-cell section, bounded faces are identified from a planar half-edge representation. Oriented cell areas are evaluated with the shoelace relation using local coordinates.

Let `B` be the signed cell-edge incidence matrix and

$$\rho_e=\frac{L_e}{t_e}.$$

The cell flexibility matrix is

$$\mathbf H=\mathbf B\,\operatorname{diag}(\rho_e)\mathbf B^T.$$

For Saint-Venant torsion in the Bredt-Batho thin-wall model,

$$\mathbf H\boldsymbol\phi=2\mathbf A_c,$$

$$\mathbf F=\mathbf B^T\boldsymbol\phi,$$

$$J_{BB}=2\mathbf A_c^T\boldsymbol\phi.$$

This is a thin-wall closed-cell formulation, not a full finite-thickness torsion solution.

For transverse shear, the section is represented by a compatible basic flow plus constant cell circulations. The circulations are obtained from the same cell compatibility structure so that the final physical flow is independent of the virtual cut used to form the open tree.

## Mixed open/closed sections

Mixed sections are decomposed into cyclic regions and open bridge branches. Closed-cell mechanics are applied to the cyclic regions, while open branches contribute their thin-strip Saint-Venant term.

The implemented torsion model uses

$$J_{\rm open}=\frac13\sum_{e\in E_{\rm open}}L_et_e^3,$$

$$J_{\rm total}=J_{BB}+J_{\rm open}.$$

Closed walls are not included again in the open-strip sum.

## Warping quantities

The shear center is used as the pole for sectorial-coordinate calculations. On an open edge,

$$\frac{d\omega}{ds}=p,$$

while on a closed wall the compatible closed-cell correction gives

$$\frac{d\omega}{ds}=p-\frac{F_e}{t_e}.$$

The sectorial field is propagated through the section graph and shifted to zero area-weighted mean,

$$\omega^\ast=\omega-rac{\int_A\omega\,dA}{A}.$$

The warping constant is then

$$C_w=\sum_e\frac{t_eL_e}{3}
\left[(\omega_1^\ast)^2+\omega_1^\ast\omega_2^\ast+(\omega_2^\ast)^2\right].$$

Small numerical residuals near theoretical zero are treated according to the documented floating-point tolerances; a near-zero computed value is not automatically interpreted as an exact mathematical zero.

## Stress recovery

With

$$D=I_xI_y-I_{xy}^2,$$

the normal stress is

$$\sigma_{zz}=\frac{N}{A}
+\frac{M_yI_x+M_xI_{xy}}{D}X
-\frac{M_xI_y+M_yI_{xy}}{D}Y
+\frac{B}{C_w}\omega^\ast.$$

The membrane shear stress is obtained from the combined transverse, secondary-warping, and closed-wall torsional shear flows divided by thickness. Open-wall Saint-Venant torsion is represented by the opposing surface stresses defined by the thin-strip model.

The reported surface von Mises stress is

$$\sigma_{vm}=\sqrt{\sigma_{zz}^2+3\tau_{\rm surface}^2}.$$

Along each straight segment, the normal stress is linear and the membrane shear stress is polynomial. Sectalix evaluates the candidate extrema analytically from the resulting polynomial conditions instead of using plot discretization to determine the reported peak.

If a positive yield stress is supplied, the proportional elastic first-yield multiplier is

$$\lambda=\frac{\sigma_{yield}}{\max\sigma_{vm}}.$$

This is an elastic first-yield indicator only. It is not a buckling check, code-based resistance, or safety certification.

## Numerical and input contract

Valid calculations require finite coordinates, positive segment lengths and thicknesses, and a connected planar topology compatible with the selected solver.

The implementation uses local origins, compensated or aligned summation where appropriate, scaled linear solves, and exponent-aware arithmetic in extreme-scale paths to reduce avoidable overflow, underflow, and cancellation. These techniques improve robustness but do not extend double precision into arbitrary precision.

Severely ill-conditioned systems are rejected rather than silently regularized. Values that cannot be represented within the supported numerical range raise explicit numerical or geometry errors.

JSON serialization can preserve the exact stored float64 representation using hexadecimal strings. Unit metadata is descriptive; it does not perform automatic conversion.
