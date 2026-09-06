# Theory and conventions

## Model and applicability

Sectalix uses straight centerline segments with constant positive thickness per segment. It assumes a homogeneous, linear-elastic thin wall. Area is t ds, not a union of finite-width polygons: corner overlaps, fillets and local plate bending are omitted. Buckling, yielding redistribution, fatigue, connections and design-code checks are outside scope. J and Cw are geometric constants; stiffnesses are GJ and ECw when compatible material properties are supplied outside this section library.

The frozen source is the executable convention. Sectalix application versions and the
v0.8 interchange schema version are versioned independently. The rename preserves the
legacy `format: "thinwallx"` wire identifier so existing files remain readable. See
[verification](VERIFICATION_BENCHMARKS.md) and [CLI reference](CLI_REFERENCE.md).

## Axes, loads and resultants

x points right, y up and z out of the section, forming a right-handed basis. X=x-cx and Y=y-cy are centroid-relative. Positive N is tensile. The exact normal-resultant convention is

$$N=\int_A\sigma\,dA,\quad M_x=-\int_A Y\sigma\,dA,\quad
M_y=\int_A X\sigma\,dA,\quad B=\int_A\omega^\ast\sigma\,dA.$$

Vx and Vy are the integrals of the transverse flow vector components. Positive Tsv is torque about +z. M_omega is the secondary warping moment with the signed secondary-flow convention below. Units are F for N,V; FL for Mx,My,Tsv,M_omega; FL² for B; F/L² for stress. Unit metadata never converts JSON quantities.

## Exact straight-segment integration

Let r(s)=r1+s u, 0≤s≤L, u=(r2-r1)/L. For any two linear fields f,g with endpoint values f1,f2,g1,g2,

$$\int_0^L f\,ds=\frac{L}{2}(f_1+f_2),$$
$$\int_0^L fg\,ds=\frac{L}{6}(2f_1g_1+f_1g_2+f_2g_1+2f_2g_2).$$

Multiplication by t gives area integrals. Thus Ai=tL, first moments are Ai(x1+x2)/2 and Ai(y1+y2)/2, and A=sum Ai. The centroid is the first-moment sum divided by A. Raw Ix integrates y², raw Iy integrates x² and raw Ixy integrates xy. Their centroidal forms are Ix=Ix_raw-A cy², Iy=Iy_raw-A cx² and Ixy=Ixy_raw-A cx cy. Production uses local reference shifts to avoid subtracting unnecessarily large raw global moments. These identities do not authorize a numerically unstable raw-global implementation.

Ixy is positive for positive correlation of X and Y. The inertia tensor for axis rotation is

$$I=\begin{pmatrix}I_x&-I_{xy}\\-I_{xy}&I_y\end{pmatrix}.$$

Its ordered eigenvalues I1≥I2 are the principal moments. With d=hypot(Ix-Iy,2Ixy), I1,2=(Ix+Iy±d)/2. The I1 axis angle is theta=atan2(-2Ixy,Ix-Iy)/2, mapped to (-pi/2,pi/2]. Isotropic axes are not uniquely defined; the frozen convention chooses zero. Do not confuse this inertia tensor with the positive-off-diagonal matrix C used below.

## Shear flow and tree equilibrium

$$C=\begin{pmatrix}I_y&I_{xy}\\I_{xy}&I_x\end{pmatrix},\qquad
C\alpha=(V_x,V_y)^T.$$

On a directed segment, with xi=s/L,

$$m_{\rm partial}(\xi)=tL\left[\xi(X_1,Y_1)^T+
\frac{\xi^2}{2}(\Delta x,\Delta y)^T\right],$$
$$q(\xi)=a_0-tL(\alpha_xX_1+\alpha_yY_1)\xi
-\frac{tL}{2}(\alpha_x\Delta x+\alpha_y\Delta y)\xi^2.$$

a0 is set by the directed subtree balance. Each physical edge is counted once; reversal reverses the scalar flow convention while preserving q u. Free ends have zero flow and junctions satisfy signed Kirchhoff balance. Integrals of q and its first moments use polynomial antiderivatives, not plot samples. C must be positive definite and sufficiently conditioned: lambda_min must exceed 10^4 epsilon_machine lambda_max.

## Torque and shear center

For a reference point rref, p=(x1-xref)uy-(y1-yref)ux is constant on a straight segment. Its +z torque is p integral(q ds). Summing gives Tref. About the centroid, ex=Tc(Vx=0,Vy=1), ey=-Tc(Vx=1,Vy=0); S=C+(ex,ey). This is consistent with ex Vy-ey Vx. Internal calculations preserve local offsets; adding a small offset to a huge global coordinate is still limited by the spacing of representable float64 coordinates.

## Closed cells and mixed topology

A planar half-edge embedding identifies bounded faces. Oriented cell areas follow the shoelace formula using a local reference. The signed cell-edge incidence B assigns opposite signs to a shared wall. With rho_e=L_e/t_e,

$$H=B\,\operatorname{diag}(\rho_e)B^T,\qquad H\phi=2A_c,$$
$$F=B^T\phi,\qquad J_{BB}=2A_c^T\phi.$$

The transverse solution uses a virtual-cut tree flow qb and constant cell circulations q0: H q0=-b, b=B integral(qb/t ds), q=qb+B^Tq0. Cut choices must not change the physical solution. These equations are the thin-wall Bredt-Batho compatibility and torque relations, not a finite-thickness torsion solution.

Tarjan bridge detection partitions mixed sections into open bridges and cyclic components. The cyclic subgraph obeys nc=|Ec|-|Vc|+kc. Independent cyclic blocks have no H coupling and are scaled and solved separately; a remote cell must not erase another through a global scale choice.

$$J_{\rm open}=\frac{1}{3}\sum_{e\in E_{\rm open}}L_et_e^3,\qquad
J_{\rm total}=J_{BB}+J_{\rm open}.$$

Closed walls are excluded from the open sum. This mixed model is the frozen thin-wall approximation, not an additional full-solid torsion correction.

## Warping

The pole is the shear center. For an open edge d omega/ds=p; on a closed edge d omega/ds=p-F_e/t_e. Signed traversal maintains continuous node values and closed-cycle compatibility. Subtract the area-weighted mean over every edge:

$$\omega^\ast=\omega-\frac{\int_A\omega\,dA}{A},\qquad
C_w=\sum_e\frac{t_eL_e}{3}
\left((\omega_1^\ast)^2+\omega_1^\ast\omega_2^\ast+(\omega_2^\ast)^2\right).$$

Mean removal makes root choice irrelevant, subject to rounding. Shear-center warping has zero bending cross moments in the underlying model. Near-zero computed values are not automatically exact mathematical zero. A zero-resistance section cannot carry nonzero B or M_omega.

## Stress recovery and peak search

With D=Ix Iy-Ixy²,

$$\sigma_{zz}=\frac{N}{A}+
\frac{M_yI_x+M_xI_{xy}}{D}X-
\frac{M_xI_y+M_yI_{xy}}{D}Y+
\frac{B}{C_w}\omega^\ast.$$

The secondary membrane flow is -M_omega/Cw times the directed warping static moment, with closed-cell compatibility corrections. Membrane stress is the combined transverse, secondary and closed-wall torsional flow divided by t. Closed torsional flow is Tsv F/Jtotal. On an open wall the opposing Saint-Venant surface stresses have magnitude |Tsv|t/Jtotal. The conservative surface envelope is |tau_membrane|+|tau_sv_surface|.

$$\sigma_{vm}=\sqrt{\sigma_{zz}^2+3\tau_{\rm surface}^2}.$$

Normal stress is linear and membrane stress quadratic on a segment. Squared von Mises stress is quartic on each fixed-sign membrane interval. Endpoints, membrane sign-change roots and admissible stationary roots of the quartic derivative are candidates. Plot discretization does not determine the reported maximum. With supplied positive yield stress, the proportional elastic first-yield multiplier is sigma_yield/max_sigma_vm; a zero-stress state has an unbounded multiplier, not an infinite engineering capacity.

## Numerical and input contract

Finite endpoints, positive length/thickness and valid connected planar topology are required. Node clustering uses the frozen tolerance rules; geometry smaller than that tolerance is not recovered by CLI inference. Input coordinates already rounded to the same float cannot be separated by local shifting.

Mantissa/exponent products (frexp/ldexp), aligned sums, scaled solves and local origins reduce intermediate overflow and underflow. They do not provide arbitrary precision. Unrepresentable required physical outputs raise numerical exceptions instead of DBL_MAX clamping. A genuinely nonzero result below the subnormal range is not silently advertised as exact zero. Subnormal relative accuracy remains limited by quantization. Ill-conditioned matrices are rejected, not regularized without consent.

JSON uses exact float64 hexadecimal strings and validates all fields through the frozen codec. Reports print 17 significant digits. Decimal CLI input is checked for nonzero-to-zero conversion and overflow. See [v0.8 format rules](V0_8_USAGE.md) for exact interchange and DXF limits.
