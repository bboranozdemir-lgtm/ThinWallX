"""ThinWallX — Python API Quickstart Example."""
from thinwallx import (
    Node,
    Segment,
    Section,
    AppliedLoads,
)

# 1. Define an open channel section (C-Channel: 50x100x50 mm, t=2 mm)
n0 = Node(x=50.0, y=100.0, id=0)
n1 = Node(x=0.0, y=100.0, id=1)
n2 = Node(x=0.0, y=0.0, id=2)
n3 = Node(x=50.0, y=0.0, id=3)

segments = [
    Segment(p1=n0, p2=n1, t=2.0, id="flange_top"),
    Segment(p1=n1, p2=n2, t=2.0, id="web"),
    Segment(p1=n2, p2=n3, t=2.0, id="flange_bot"),
]
section = Section(segments)

# 2. Access cross-section characteristics
print("--- Cross-Section Properties ---")
print(f"Area: A = {section.area:.2f} mm^2")
print(f"Centroid: C = ({section.cx:.2f}, {section.cy:.2f}) mm")
print(f"Moments of Inertia: Ix = {section.Ix:.2f}, Iy = {section.Iy:.2f} mm^4")
print(f"Saint-Venant Torsion: J = {section.J:.2f} mm^4")
print(f"Warping Constant: Cw = {section.Cw:.2f} mm^6")

# 3. Perform stress recovery under combined external loads
loads = AppliedLoads(N=10000.0, Vy=5000.0, Mx=250000.0, sigma_yield=355.0)
results = section.calculate_stresses(loads)
print("\n--- Stress Recovery Analysis ---")
print(f"Peak Von Mises Stress: {results.max_sigma_vm:.2f} MPa")
print(f"Elastic Load Factor: lambda = {results.load_factor:.3f}")
