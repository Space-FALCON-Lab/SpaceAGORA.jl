# Documentation figures

These scripts regenerate the figures in `docs/src/assets/figures/`. The
documentation build does not run them: when the code or page behind a figure
changes, rerun its script and commit the new files. Run every command from the
repository root.

| Figure | Used on | Script | Needs |
| --- | --- | --- | --- |
| `configuration_anatomy.svg` | Configuring a Simulation | `diagrams.jl` | Julia only |
| `environments_map.svg` | Installation & Environment | `diagrams.jl` | Julia only |
| `mission_workflow_timeline.svg` | Mission Workflows | `diagrams.jl` | Julia only |
| `quickstart_altitude_speed.png`, `quickstart_inertial_trajectory.png` | Quickstart | `quickstart_plots.jl` | the examples environment |
| `odyssey_panel_cap_comparison.png` | Odyssey Guidance with a Frozen Atmosphere | `odyssey_figures.jl` | a run of the Odyssey example and the examples environment |
| `aerobraking_pass_anatomy.png` | Energy Depletion Guidance | `odyssey_figures.jl` | a run of the Odyssey example and the examples environment |
| `viewer_overview.jpg`, `viewer_heating_closeup.jpg`, `viewer_face_inspection.jpg` | Interactive Visualization (the close-up also heads the repository README) | `viewer_screenshots.jl` | the viewer environment, Google Chrome, Node.js 22 or newer |

The run-flow diagram on the home page and on How a Simulation Fits Together,
`docs/src/assets/architecture_run_flow.svg`, does not come from these scripts;
`svg.jl` only borrows its colours.

## Commands

The diagrams are drawn by `diagrams.jl` with the small SVG writer in `svg.jl`:

```text
julia --startup-file=no scripts/docs_figures/diagrams.jl
```

The quickstart plots come from a new quickstart run in a temporary folder:

```text
julia --project=examples scripts/docs_figures/quickstart_plots.jl
```

The Odyssey figures read the tutorial's default comparison of the 90 and 30
degree panel caps. Set up the Odyssey example environment as the tutorial
describes, then:

```text
julia --project=examples/odyssey_surrogate_env examples/odyssey_surrogate.jl --output=output/odyssey_figures
julia --project=examples scripts/docs_figures/odyssey_figures.jl output/odyssey_figures
```

The viewer screenshots run the Earth example from Interactive Visualization and
a Mars aerobraking pass drawn with the NASA Mars Odyssey model from
`data/scene/models/`, then capture the pages in headless Chrome through
`capture_viewer.mjs`. Set `CHROME` if the browser is not in its default place:

```text
julia --project=. scripts/setup_visualizer.jl .visualizer-env
julia --project=.visualizer-env scripts/docs_figures/viewer_screenshots.jl
```

The screenshots depend on the camera reaching the same framing, so check them
by eye after a rerun. The Mars Odyssey model is a NASA 3D Resources asset; keep
its credit with any page that shows it, as `data/scene/models/README.md`
describes.
