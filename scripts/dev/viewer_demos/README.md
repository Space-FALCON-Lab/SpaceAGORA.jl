# Viewer demo drivers

Scripts that produce the interactive-viewer pages used to exercise and show
the visualizer. Each writes a run and its page under
`output/viewer_demos/<case>/` (override the root with
`SPACEAGORA_VIEWER_DEMO_OUT`) and, where a case has one, the artifact
variant `artifact.html` built by `build_cdn_page.py`.

| Driver | Case | Assets needed |
|---|---|---|
| `earth_4day.jl` | `examples/AGORA_Earth.jl` propagated for four days | none |
| `iss_demo.jl` | ISS in a 408 km orbit drawn with NASA's ISS (B) model | none |
| `iss_hypr.jl` | HyPR/PSO rendezvous path to the ISS, LQ-MPC tracked, planned paths overlaid | none |
| `arm_demo.jl` | robot-arm cloth dynamics | none |
| `odyssey_two_orbits.jl` | two Odyssey aerobraking orbits with the energy-depletion guidance | none |
| `magellan_aerobraking.jl` | Magellan aerobraking at Venus, two orbits from 1993-05-26, SPICE ghost | GRAM (Venus), Magellan `AEROBRAK.BSP` (auto-downloaded) |
| `cassini_titan_flyby.jl [TA\|T5\|all]` | Cassini's TA and T5 Titan flybys, Titan-centred, SPICE ghost | GRAM (Titan), Cassini SCPSE kernels (auto-downloaded) |
| `apollo11_lunar_orbit.jl` | Apollo 11 lunar module in the parking orbit, two orbits | SPICE (Moon frames) |

The mission cases share `common.jl`: kernel download into
`SPICE/spk/missions/`, `spkezr` sampling relative to the central body in
J2000 (the integrator's frame), apoapsis and closest-approach searches for
the initial epoch, `spice_reference` (the ghost table for
`export_visualization(...; references=...)`), and a separation report
between the run and the kernel over the saved times. Rerun a case with
`SPACEAGORA_DEMO_FORCE=1` to redo the simulation instead of reusing the
results on disk.

`build_cdn_page.py <viewer.html> <out.html> [title heading orbit span foot]`
turns an exported page into the variant the claude.ai artifact host can
show: three.js from jsdelivr and the viewer modules concatenated into one
script, because that host blocks `data:` script URLs. Local pages need no
such step; open `simulation_results_viewer.html` from disk.
