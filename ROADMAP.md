# Sectalix — Roadmap

## Release seal — v1.0.1 (historical core v1.0.0 preserved)

v0.1–v1.0 teslimat kapsamı **COMPLETED & FROZEN** olarak mühürlenmiştir.
559 testlik v1.0 çekirdeği ve testleri değiştirilemez; v1.0.1 yeniden adlandırma
doğrulaması bu çekirdeğe ek olarak bir uyumluluk testi içerir. Yol A yalnızca
dağıtım ve CI/CD hazırlığıdır.

| Phase | Delivered scope | Status |
|---|---|---|
| v0.1 | Geometry, centroid, inertia and principal axes | COMPLETED & FROZEN |
| v0.2 | Open shear flow | COMPLETED & FROZEN |
| v0.3 | Shear center | COMPLETED & FROZEN |
| v0.4 | Open torsion and warping | COMPLETED & FROZEN |
| v0.5 | Closed and multicell mechanics | COMPLETED & FROZEN |
| v0.6 | Mixed open/closed mechanics | COMPLETED & FROZEN |
| v0.7 | Combined stress recovery | COMPLETED & FROZEN |
| v0.8 | JSON, DXF and technical plots | COMPLETED & FROZEN |
| v0.9 | Release preparation incorporated in v1.0; no separate versioned implementation record | COMPLETED & FROZEN within v1.0 |
| v1.0 | CLI, calculation reports, documentation and packaging | COMPLETED & FROZEN |

The v0.9 row does not assert a separate historical release or test run.
Original scope summaries below are retained for traceability, not future-work authorization.

## v0.1 — Geometry and Basic Section Properties
- centerline segment representation
- constant thickness per segment
- connected open-section topology, including branched open sections
- area
- centroid
- Ix, Iy, Ixy
- principal section properties
- analytical benchmarks
- explicit geometry validation

## v0.2 — Open-Section Shear Flow
- open thin-walled sections only
- first-moment / shear-flow formulation
- shear-flow distribution
- benchmark cases

## v0.3 — Shear Center
- shear-center computation
- sign-convention validation
- technical visualization
- benchmark cases

## v0.4 — Torsion and Warping
- Saint-Venant torsion constant J
- sectorial coordinate / warping function
- warping constant Cw
- benchmark cases

## v0.5 — Closed and Multi-Cell Sections
- closed-section topology support
- closed-cell shear flow
- compatibility equations
- multi-cell systems
- variable segment thickness support where appropriate

## v1.0 — Production-Quality Technical Release
- stable open/closed/multi-cell workflow
- JSON input
- optional DXF import
- complete validation suite
- technical plots and reports
- documentation of conventions, formulas, assumptions, and limitations

## Rule
The roadmap is descriptive only.
Codex must never advance to another version unless `ACTIVE_PHASE.md` explicitly names that version.
