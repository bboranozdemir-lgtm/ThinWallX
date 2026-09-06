"""Reproducible v0.8 I/O and headless plot examples (no GUI).

Run from the repository with src on PYTHONPATH:
    python examples/v08_io_and_plots.py --output examples/v08_output
"""
from __future__ import annotations

import argparse
from pathlib import Path

from thinwallx import AppliedLoads
from thinwallx.dxf import DxfImportOptions, ThicknessMap, read_dxf
from thinwallx.serialization import UnitSystem, write_json, read_json
from thinwallx.plotting import (
    plot_geometry, plot_shear_flow, plot_stresses,
    PlotStyle, GeometryPlotOptions, ShearPlotOptions, StressPlotOptions,
)


def generate(output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    fixtures = Path(__file__).resolve().parents[1] / "tests" / "fixtures" / "dxf"
    units = UnitSystem("mm", "N")
    for name in ("open_l_r12", "rectangle_r12", "barbell_ac1015"):
        imported = read_dxf(fixtures / (name + ".dxf"),
                            options=DxfImportOptions("mm", "mm", ThicknessMap()))
        section = imported.section
        write_json(section, output / (name + ".json"), units=units, overwrite=True)
        restored = read_json(output / (name + ".json"))
        assert restored.value.segments == section.segments
        loads = AppliedLoads(N=1.0, Mx=0.1, My=0.05, sigma_yield=250.0)
        stress = section.calculate_stresses(loads)
        write_json(loads, output / (name + "_loads.json"), units=units, overwrite=True)
        write_json(stress, output / (name + "_stress.json"), units=units, overwrite=True)
        flow = section.calculate_shear_flow(vx=0., vy=1.)
        for extension in ("png", "svg"):
            plot_geometry(section, output / (name + "_geometry." + extension), units=units, overwrite=True,
                          options=GeometryPlotOptions(show_segment_ids=True, show_shear_center=True,
                                                      show_principal_axes=True))
            plot_shear_flow(section, flow, output / (name + "_shear." + extension), units=units,
                            overwrite=True, options=ShearPlotOptions())
            plot_stresses(section, stress, output / (name + "_stress." + extension), units=units,
                          overwrite=True, options=StressPlotOptions(style=PlotStyle(figsize=(12., 6.))))
        print(f"{name}: {type(section).__name__}, {len(section.segments)} edges, A={section.area:.17g}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    generate(parser.parse_args().output)

