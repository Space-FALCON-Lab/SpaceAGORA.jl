# Interactive visualization ownership

The viewer reads saved trajectories and a scene sidecar. It is opt-in and
does not change the propagation equations, atmosphere model or controller.

| Owner | Responsibility |
| --- | --- |
| `src/analysis/visualization/scene/` | Scene metadata, geometry conversion, body rotation samples, atmosphere overlays, saved-row export and HTML assembly |
| `src/vehicle/structure/mesh_geometry.jl` | Shared mesh reading, bounds, articulation and point sampling |
| `src/simulation/callbacks/save_fields.jl` | Optional viewer records sampled alongside the existing results |
| `src/simulation/engine/persistence.jl` | Write the sidecar after a flagged run |
| `src/simulation/campaigns/monte_carlo_visualization.jl` | Per-sample output directories and ensemble export over the existing campaign runner |
| `viewer/` | Browser rendering, playback and standalone CSV/JSON import |

With visualization disabled, existing default save fields and run behavior
remain unchanged. Enabling it adds pose, density, Sun and reported thruster
fields when available. Explicit save-field lists receive missing visual
fields once. Exporting an already saved bundle does not run a simulation.

This first integration excludes the GRAM epoch reconstruction, aerodynamic
mesh force model, lunar landing control, touchdown handling and plume-surface
physics developed alongside the viewer. Their source remains in PR121.

Adapted from Evan Yu's PR121 at `c03d1dce`; see the user guide and
`viewer/README.md` for usage and third-party renderer attribution.
