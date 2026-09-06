"""Headless deterministic technical plots; mechanics remain in frozen evaluators.

All dimensional quantities are normalized before passing to matplotlib.
Rational local differences preserve represented coordinates at large origins.
"""
from __future__ import annotations

from dataclasses import dataclass, field, fields
from fractions import Fraction
import math
from pathlib import Path
import platform
from typing import Any, Callable

import numpy as np

from thinwallx import Section, ClosedSection, MixedSection
from thinwallx.closed_shear_flow import ClosedShearFlowResult, ClosedSegmentShearFlow
from thinwallx.mixed_shear_flow import MixedShearFlowResult, MixedSegmentShearFlow
from thinwallx.shear_flow import ShearFlowResult, SegmentShearFlow
from thinwallx.stress import StressRecoveryResult
from thinwallx.serialization import UnitSystem, _atomic_write, to_dict

Flow = ShearFlowResult | ClosedShearFlowResult | MixedShearFlowResult
SegmentFlow = SegmentShearFlow | ClosedSegmentShearFlow | MixedSegmentShearFlow


def _positive(value: object, name: str, integer: bool = False) -> None:
    if integer:
        valid = type(value) is int and value > 0
    else:
        valid = type(value) in (int, float) and math.isfinite(value) and value > 0
    if not valid:
        raise ValueError(f"{name}: expected positive finite {'integer' if integer else 'number'}")


def _check_inputs(section: object, options: object, expected: type) -> None:
    if type(section) not in (Section, ClosedSection, MixedSection):
        raise TypeError("Expected Section, ClosedSection or MixedSection")
    if not isinstance(options, expected):
        raise TypeError(f"options must be {expected.__name__}")


@dataclass(frozen=True)
class PlotStyle:
    figsize: tuple[float, float] = (8.0, 6.0)
    dpi: int = 150
    font_family: str = "DejaVu Sans"
    font_size: float = 9.0
    background: str = "white"
    line_width: float = 1.2
    grid: bool = False
    max_samples: int = 200000
    svg_hashsalt: str = "thinwallx-v0.8"
    title: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.figsize, tuple) or len(self.figsize) != 2:
            raise ValueError("figsize must contain two values")
        for v in self.figsize:
            _positive(v, "figsize")
        for k in ("dpi", "max_samples"):
            _positive(getattr(self, k), k, True)
        for k in ("font_size", "line_width"):
            _positive(getattr(self, k), k)
        if type(self.grid) is not bool or self.title is not None and not isinstance(self.title, str):
            raise ValueError("Invalid grid/title")
        for k in ("font_family", "background", "svg_hashsalt"):
            if not isinstance(getattr(self, k), str) or not getattr(self, k):
                raise ValueError(f"Invalid {k}")
        if any(v*self.dpi > 20000 for v in self.figsize) or self.figsize[0]*self.figsize[1]*self.dpi**2 > 100000000:
            raise ValueError("Plot pixel budget exceeded")


def _options(value: Any) -> None:
    if not isinstance(value.style, PlotStyle):
        raise TypeError("style must be PlotStyle")
    for f in fields(value):
        if f.name.startswith("show_") and type(getattr(value, f.name)) is not bool:
            raise ValueError(f"{f.name}: expected bool")


@dataclass(frozen=True)
class GeometryPlotOptions:
    style: PlotStyle = field(default_factory=PlotStyle)
    show_thickness: bool = True
    show_node_ids: bool = False
    show_segment_ids: bool = False
    show_centroid: bool = True
    show_shear_center: bool = False
    show_principal_axes: bool = False
    thickness_display_factor: float = 1.0

    def __post_init__(self) -> None:
        _options(self)
        _positive(self.thickness_display_factor, "thickness_display_factor")


@dataclass(frozen=True)
class ShearPlotOptions:
    style: PlotStyle = field(default_factory=PlotStyle)
    samples_per_segment: int = 33
    arrows_per_segment: int = 3
    show_colorbar: bool = True
    show_centroid: bool = True
    arrow_mode: str = "normalized"

    def __post_init__(self) -> None:
        _options(self)
        _positive(self.samples_per_segment, "samples_per_segment", True)
        _positive(self.arrows_per_segment, "arrows_per_segment", True)
        if self.samples_per_segment < 2 or self.arrow_mode != "normalized":
            raise ValueError("Invalid shear sampling/arrow mode")


@dataclass(frozen=True)
class StressPlotOptions:
    style: PlotStyle = field(default_factory=PlotStyle)
    quantity: str = "sigma_vm"
    samples_per_segment: int = 65
    show_peak: bool = True
    abscissa: str = "segment"
    segment_order: tuple[int, ...] | None = None

    def __post_init__(self) -> None:
        _options(self)
        _positive(self.samples_per_segment, "samples_per_segment", True)
        if self.samples_per_segment < 2 or self.quantity not in ("sigma_zz", "sigma_vm", "tau_membrane", "tau_surface"):
            raise ValueError("Invalid stress quantity/sampling")
        if self.abscissa not in ("segment", "concatenated"):
            raise ValueError("Invalid abscissa")


@dataclass(frozen=True)
class PlotReport:
    path: Path
    format: str
    pixel_size: tuple[int, int] | None
    segment_count: int
    sample_count: int
    display_origin: tuple[float, float]
    display_length_exponent: int
    field_scale: str | None
    skipped_annotations: tuple[str, ...]
    renderer_versions: dict[str, str]


def _power(k: int) -> Fraction:
    return Fraction(2**k) if k >= 0 else Fraction(1, 2**(-k))


@dataclass(frozen=True)
class _Frame:
    origin: tuple[float, float]
    exponent: int

    def point(self, x: float, y: float) -> tuple[float, float]:
        divisor = _power(self.exponent)
        return (float((Fraction(x)-Fraction(self.origin[0]))/divisor),
                float((Fraction(y)-Fraction(self.origin[1]))/divisor))

    def length(self, value: float) -> float:
        return float(Fraction(value)/_power(self.exponent))


def _frame(section: Section, extras: tuple[tuple[float, float], ...] = (), band: float = 0.) -> _Frame:
    points = [n.coords for s in section.segments for n in (s.p1, s.p2)]
    origin = min(points)
    span = max(abs(Fraction(v)-Fraction(origin[j])) for p in points+list(extras) for j, v in enumerate(p))
    if band:
        span = max(span, max(Fraction(s.t)*Fraction(band) for s in section.segments))
    if not span:
        raise ValueError("Cannot plot collapsed geometry")
    exponent = span.numerator.bit_length()-span.denominator.bit_length()+1
    return _Frame(origin, exponent)


def _sample_points(section: Section, frame: _Frame, i: int, xi: list[float]) -> np.ndarray:
    s = section.segments[i]
    a, b = np.array(frame.point(*s.p1.coords)), np.array(frame.point(*s.p2.coords))
    return np.array([(1-u)*a+u*b for u in xi])


def _signature(section: Section) -> list[tuple[object, ...]]:
    return [(s.p1.coords, s.p2.coords, s.t, s.id) for s in section.segments]


def _backend() -> tuple[Any, Any, Any]:
    try:
        import matplotlib
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg
    except ImportError as exc:
        raise ImportError("Technical plots require the optional matplotlib dependency (thinwallx[plots]).") from exc
    return matplotlib, Figure, FigureCanvasAgg


def _draw_base(ax: Any, section: Section, frame: _Frame, units: UnitSystem, centroid: bool) -> None:
    for s in section.segments:
        a, b = frame.point(*s.p1.coords), frame.point(*s.p2.coords)
        ax.plot([a[0], b[0]], [a[1], b[1]], color="#444444", zorder=1)
    if centroid:
        ax.plot(*frame.point(section.cx, section.cy), "o", color="black", label="C", zorder=5)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(f"(x - {frame.origin[0]:.17g}) / 2^{frame.exponent} [{units.length}]")
    ax.set_ylabel(f"(y - {frame.origin[1]:.17g}) / 2^{frame.exponent} [{units.length}]")
    ax.margins(.12)


def _render(section: Section, path: str | Path, style: PlotStyle, units: UnitSystem | None,
            frame: _Frame, callback: Callable[..., tuple[int, float | None, list[str]]],
            overwrite: bool, panels: int = 1) -> PlotReport:
    if type(section) not in (Section, ClosedSection, MixedSection):
        raise TypeError("Expected Section, ClosedSection or MixedSection")
    # Fresh validation avoids mutating caller caches.
    kwargs = {"node_tolerance": section._node_tolerance}
    if type(section) is not Section:
        kwargs["safety_factor"] = section._safety_factor
    type(section)(section.segments, validate=True, **kwargs)
    units = units or UnitSystem()
    if not isinstance(units, UnitSystem):
        raise TypeError("units must be UnitSystem")
    fmt = Path(path).suffix.lower()[1:]
    if fmt not in ("png", "svg"):
        raise ValueError("Plot output must be .png or .svg")
    if Path(path).exists() and not overwrite:
        raise FileExistsError(path)
    mpl, Figure, Canvas = _backend()
    # Reset all render settings locally, without touching the caller's backend.
    rc = {k: v for k, v in mpl.rcParamsDefault.items() if k not in ("backend", "backend_fallback")}
    rc.update({"font.family": style.font_family, "font.size": style.font_size,
               "figure.facecolor": style.background, "axes.facecolor": style.background,
               "lines.linewidth": style.line_width, "axes.grid": style.grid,
               "text.usetex": False, "text.parse_math": False, "svg.hashsalt": style.svg_hashsalt,
               "svg.fonttype": "none", "savefig.bbox": None})
    with mpl.rc_context(rc):
        fig = Figure(figsize=style.figsize, dpi=style.dpi)
        Canvas(fig)
        try:
            axes = [fig.add_subplot(1, panels, i+1) for i in range(panels)]
            count, scale, notes = callback(fig, axes, units, mpl)
            if count > style.max_samples:
                raise ValueError("Plot sample budget exceeded")
            if scale is not None and not math.isfinite(scale):
                raise OverflowError("Nonfinite plot field scale")
            if style.title is not None:
                fig.suptitle(style.title)
            # Reserve space for a one-panel colorbar's ticks AND dimensional label.
            right = .84 if panels == 1 and scale is not None else .94
            fig.subplots_adjust(left=.13, bottom=.18, right=right, top=.87, wspace=.55)
            metadata = {"Software": "ThinWallX v0.8"} if fmt == "png" else {"Date": None, "Creator": "ThinWallX v0.8"}
            target = _atomic_write(path, lambda p: fig.savefig(p, format=fmt, dpi=style.dpi, metadata=metadata), overwrite)
            for s in section.segments:
                if frame.length(s.length)*min(style.figsize)*style.dpi < 1:
                    notes.append("Some segments are below pixel resolution; they are not physical zeros.")
                    break
            return PlotReport(target, fmt, tuple(int(v*style.dpi) for v in style.figsize) if fmt == "png" else None,
                              len(section.segments), count, frame.origin, frame.exponent,
                              None if scale is None else float(scale).hex(), tuple(notes),
                              {"python": platform.python_version(), "numpy": np.__version__, "matplotlib": mpl.__version__})
        finally:
            fig.clear()


def plot_geometry(section: Section | ClosedSection | MixedSection, path: str | Path, *,
                  options: GeometryPlotOptions | None = None, units: UnitSystem | None = None,
                  overwrite: bool = False) -> PlotReport:
    o = GeometryPlotOptions() if options is None else options
    _check_inputs(section, o, GeometryPlotOptions)
    center = tuple(section.shear_center) if o.show_shear_center else None
    frame = _frame(section, (center,) if center is not None else (), o.thickness_display_factor if o.show_thickness else 0.)
    if 2*len(section.segments) > o.style.max_samples:
        raise ValueError("Plot sample budget exceeded")

    def draw(fig: Any, axes: list[Any], unit: UnitSystem, mpl: Any) -> tuple[int, None, list[str]]:
        from matplotlib.patches import Polygon
        ax = axes[0]
        notes: list[str] = []
        _draw_base(ax, section, frame, unit, o.show_centroid)
        for i, s in enumerate(section.segments):
            a, b = np.array(frame.point(*s.p1.coords)), np.array(frame.point(*s.p2.coords))
            if o.show_thickness:
                # Display coordinates can coincide for sub-pixel segments.
                # Physical tangent must not be obtained from that collapsed pair.
                tangent = np.array([(s.p2.x-s.p1.x)/s.length, (s.p2.y-s.p1.y)/s.length])
                normal = np.array([-tangent[1], tangent[0]])
                width = float(Fraction(s.t)*Fraction(o.thickness_display_factor)/_power(frame.exponent)/2)
                vertices = np.array([a+width*normal, b+width*normal, b-width*normal, a-width*normal])
                if not np.isfinite(vertices).all():
                    raise OverflowError("Thickness band is not representable")
                ax.add_patch(Polygon(vertices, facecolor="#cccccc", edgecolor="none", alpha=.6, zorder=0))
            if o.show_segment_ids:
                ax.text(*((a+b)/2), f"{i}: {s.id!r}")
        if o.show_node_ids:
            seen = set()
            for s in section.segments:
                for n in (s.p1, s.p2):
                    if n.coords not in seen:
                        ax.text(*frame.point(*n.coords), repr(n.id), va="bottom")
                        seen.add(n.coords)
        if center is not None:
            ax.plot(*frame.point(*center), "x", color="purple", label="S", zorder=6)
        if o.show_principal_axes:
            if section.is_degenerate:
                notes.append("Principal direction undefined (isotropic inertia).")
                ax.text(.02,.97,notes[-1],transform=ax.transAxes,va="top")
            else:
                c = np.array(frame.point(section.cx, section.cy))
                for j, theta in enumerate((section.theta_p, section.theta_p+math.pi/2)):
                    direction = np.array([math.cos(theta), math.sin(theta)])*.25
                    ax.plot([c[0]-direction[0],c[0]+direction[0]], [c[1]-direction[1],c[1]+direction[1]],
                            "--",label=f"I{j+1}", color=("#0072B2", "#D55E00")[j])
        if o.thickness_display_factor != 1:
            ax.set_title(f"Thickness displayed x{o.thickness_display_factor:g}; centerline model")
        else:
            ax.set_title("Centerline model; junction overlaps neglected")
        if ax.get_legend_handles_labels()[0]:
            ax.legend()
        return 2*len(section.segments), None, notes
    return _render(section,path,o.style,units,frame,draw,overwrite)


def _flow_candidates(flow: SegmentFlow) -> list[float]:
    """Quadratic extrema: derivative a1+2*a2*u=0 on each analysis piece."""
    pieces: list[tuple[tuple[float, float, float], float, float]]
    if isinstance(flow, SegmentShearFlow):
        pieces = [((flow.a0,flow.a1,flow.a2),0.,1.)]
    elif flow.is_chord:
        pieces = [(flow._poly_a,0.,flow._cut_param),(flow._poly_b,flow._cut_param,1.)]
    else:
        pieces = [(flow._poly_tree,0.,1.)]
    values = {0.,1.}
    for poly,a,b in pieces:
        if poly is None or not all(math.isfinite(v) for v in poly):
            raise ValueError("Invalid shear polynomial")
        values.update((a,b))
        if poly[2]:
            root = -Fraction(float(poly[1]))/(2*Fraction(float(poly[2])))
            if 0 < root < 1:
                values.add(float(Fraction(a)+(Fraction(b)-Fraction(a))*root))
    return sorted(values)


def _finite_field(value: float) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise OverflowError("Nonfinite plot field")
    return value


def _contour(fig: Any, ax: Any, positions: list[np.ndarray], values: list[list[float]],
             scale: float, signed: bool, label: str, show_bar: bool = True) -> None:
    from matplotlib.collections import LineCollection
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    norm = Normalize(-1 if signed else 0, 1)
    cmap = "coolwarm" if signed else "viridis"
    for p, v in zip(positions,values):
        normalized = [x/scale if scale else 0. for x in v]
        segments = np.stack((p[:-1],p[1:]),axis=1)
        # Average already-normalized colors; not physical stress integration.
        colors = [(a+b)/2 for a,b in zip(normalized,normalized[1:])]
        ax.add_collection(LineCollection(segments,array=np.array(colors),cmap=cmap,norm=norm,linewidths=4,zorder=3))
    if show_bar:
        bar = fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),ax=ax,fraction=.045,pad=.04)
        if scale:
            bar.set_label(f"{label} / {scale:.8g}")
        else:
            bar.set_ticks([0])
            bar.set_ticklabels(["0"])
            bar.set_label(label+" = 0")


def plot_shear_flow(section: Section | ClosedSection | MixedSection, flow: Flow, path: str | Path, *,
                    options: ShearPlotOptions | None = None, units: UnitSystem | None = None,
                    overwrite: bool = False) -> PlotReport:
    o = ShearPlotOptions() if options is None else options
    _check_inputs(section, o, ShearPlotOptions)
    if not isinstance(flow,(ShearFlowResult,ClosedShearFlowResult,MixedShearFlowResult)):
        raise TypeError("Expected shear-flow result")
    if _signature(section) != _signature(flow.section) or len(flow.segment_flows)!=len(section.segments):
        raise ValueError("Shear result does not match section")
    if len(section.segments)*(o.samples_per_segment+6+o.arrows_per_segment)>o.style.max_samples:
        raise ValueError("Plot sample budget exceeded")
    frame=_frame(section)
    positions,values,arrows=[],[],[]
    for i,sf in enumerate(flow.segment_flows):
        if sf.segment != section.segments[i]:
            raise ValueError("Shear segment does not match section")
        xi=sorted(set(np.linspace(0.,1.,o.samples_per_segment).tolist()+_flow_candidates(sf)))
        values.append([_finite_field(sf.q_at_xi(u)) for u in xi])
        positions.append(_sample_points(section,frame,i,xi))
        us=[(j+1)/(o.arrows_per_segment+1) for j in range(o.arrows_per_segment)]
        arrows.append((us,[_finite_field(sf.q_at_xi(u)) for u in us]))
    scale=max(abs(v) for vals in values for v in vals)

    def draw(fig: Any, axes: list[Any], unit: UnitSystem, mpl: Any) -> tuple[int,float,list[str]]:
        ax=axes[0]
        _draw_base(ax,section,frame,unit,o.show_centroid)
        _contour(fig,ax,positions,values,scale,True,f"q [{unit.force}/{unit.length}]",o.show_colorbar)
        for i,(us,vals) in enumerate(arrows):
            pts=_sample_points(section,frame,i,us)
            s=section.segments[i]
            dx,dy=s.p2.x-s.p1.x,s.p2.y-s.p1.y
            tangent=(dx/s.length,dy/s.length)
            if scale:
                for p,v in zip(pts,vals):
                    if v:
                        delta=np.array(tangent)*(.08*(v/scale))
                        ax.annotate("",xy=p+delta,xytext=p,arrowprops={"arrowstyle":"->","color":"black"})
        ax.set_title("Physical q * tangent; arrow lengths normalized" if scale else "q = 0")
        return sum(map(len,values))+sum(len(u) for u,_ in arrows),scale,[]
    return _render(section,path,o.style,units,frame,draw,overwrite)


def plot_stresses(section: Section | ClosedSection | MixedSection, result: StressRecoveryResult,
                  path: str | Path, *, options: StressPlotOptions | None = None,
                  units: UnitSystem | None = None, overwrite: bool = False) -> PlotReport:
    o=StressPlotOptions() if options is None else options
    _check_inputs(section, o, StressPlotOptions)
    to_dict(result)  # validates archived coefficients/peaks without mutating them
    keys=[s.id if s.id is not None else i for i,s in enumerate(section.segments)]
    if len(set(keys))!=len(keys) or set(keys)!=set(result.segment_profiles):
        raise ValueError("Stress result IDs do not match section")
    order=o.segment_order if o.segment_order is not None else tuple(range(len(keys)))
    if not isinstance(order,tuple) or any(type(i) is not int for i in order) or sorted(order)!=list(range(len(keys))):
        raise ValueError("segment_order must be a complete permutation")
    if len(keys)*(o.samples_per_segment+5)>o.style.max_samples:
        raise ValueError("Plot sample budget exceeded")
    frame=_frame(section)
    positions,values,abscissae=[],[],[]
    for i in order:
        s,p=section.segments[i],result.segment_profiles[keys[i]]
        if s.length!=p.length or s.t!=p.thickness:
            raise ValueError("Stress profile geometry mismatch")
        ss=sorted(set([float(u)*p.length for u in np.linspace(0.,1.,o.samples_per_segment)]+[0.,p.length,p.s_peak]))
        # Include the membrane quadratic's stationary point for exact color bounds.
        a,b,_=p.tau_membrane_coeffs
        if a:
            root=-Fraction(float(b))/(2*Fraction(float(a)))
            if 0<root<Fraction(p.length):
                ss=sorted(set(ss+[float(root)]))
        vals=[_finite_field(getattr(p,"eval_"+o.quantity)(u)) for u in ss]
        values.append(vals)
        positions.append(_sample_points(section,frame,i,[u/p.length for u in ss]))
        abscissae.append([frame.length(u) for u in ss])
    peak_i=keys.index(result.peak_segment_id)
    peak_seg=section.segments[peak_i]
    u=Fraction(result.peak_s)/Fraction(peak_seg.length)
    expected=tuple(float((1-u)*Fraction(a)+u*Fraction(b)) for a,b in zip(peak_seg.p1.coords,peak_seg.p2.coords))
    if any(abs(Fraction(g)-Fraction(e)) > Fraction(peak_seg.length)*Fraction(1e-10)
           for g,e in zip(result.peak_location_xy,expected)):
        raise ValueError("Stress peak coordinates do not match section")
    scale=result.max_sigma_vm if o.quantity=="sigma_vm" else max(abs(v) for vs in values for v in vs)
    signed=o.quantity in ("sigma_zz","tau_membrane")

    def draw(fig: Any, axes: list[Any], unit: UnitSystem, mpl: Any) -> tuple[int,float,list[str]]:
        ax,profile_ax=axes
        _draw_base(ax,section,frame,unit,True)
        label=f"{o.quantity} [{unit.force}/{unit.length}^2]"
        _contour(fig,ax,positions,values,scale,signed,label)
        start=0.
        for i,ss,vs in zip(order,abscissae,values):
            xs=[v+start for v in ss] if o.abscissa=="concatenated" else ss
            profile_ax.plot(xs,[v/scale if scale else 0. for v in vs],label=f"{i}: {keys[i]!r}")
            if o.abscissa=="concatenated":
                profile_ax.axvline(start,color="gray",linewidth=.5)
                start+=frame.length(section.segments[i].length)
        profile_ax.set_xlabel(("Segment-concatenated coordinate\nNOT continuous perimeter" if o.abscissa=="concatenated"
                               else "Local segment s")+f" / 2^{frame.exponent} [{unit.length}]")
        profile_ax.set_ylabel(label+(f" / {scale:.8g}" if scale else " = 0"))
        profile_ax.legend()
        ax.set_title("Centerline fiber contour\n(not a through-thickness field)")
        if o.show_peak:
            ax.plot(*frame.point(*result.peak_location_xy),"*",color="red",markersize=10,zorder=6)
            fig.text(.05,.03,f"Analytical max sigma_vm = {result.max_sigma_vm:.17g}; segment {result.peak_segment_id!r}")
        return sum(map(len,values)),scale,[]
    return _render(section,path,o.style,units,frame,draw,overwrite,panels=2)
