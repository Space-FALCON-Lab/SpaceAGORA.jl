# Interactive Visualization

The viewer turns saved results into an HTML page: orbit the central body,
scrub through time, select a spacecraft and inspect its trajectory. The page
opens from disk and can be shared with someone who does not have Julia.

## Make your first page

After completing the [Quickstart](quickstart.md), run this from the repository
root to put the no-GRAM example and its viewer in a dedicated output directory:

```julia
using SpaceAGORA
SpaceAGORA.run_cli(["run", "--example=AGORA_Earth_NoGRAM.jl", "--visualize",
                   "--output-dir=output/viewer-earth"])
```

Open `output/viewer-earth/simulation_results_viewer.html` under the repository.
The same directory contains the normal results and
`simulation_results_scene.json`. Repeating the command overwrites those files.
Choose another `--output-dir` to keep an earlier run. See
[Simulation Outputs](outputs.md) for the full output behavior.

For a configuration you already built:

```julia
run_simulation(args; visualization=true)
```

For a completed run that already has its scene sidecar:

```julia
export_visualization("output/viewer-earth/simulation_results"; max_frames=2000)
SpaceAGORA.run_cli(["visualize", "--run=output/viewer-earth"])
```

`with_visualization_scene(args, true)` sets the sidecar flag when preparing
a configuration. This saves the extra fields and sidecar; use
`export_visualization` afterwards to build the HTML page. Existing runs
without a sidecar need their configuration to create one, or use the
standalone data-file viewer below.

## Read the page

- Drag to orbit, scroll to zoom, and use the timeline to play or scrub.
- Select a spacecraft and press F to follow it. Escape clears the selection.
- Switch between inertial and planet-fixed views. Ground tracks follow the
  body's rotation sampled from the run's own frame model.
- Nearby spacecraft show link boxes and articulated parts. A run without
  attitude data uses an illustrative velocity-aligned orientation.
- Quaternions are scalar-last. Instantaneous altitude is geodetic;
  `sc1_periapsis_altitude` and the default ensemble color are spherical
  periapsis altitude, as defined on the outputs page.

The saved heat-rate history comes from the simulation. Surface heating colors
are an illustrative viewer overlay, not an additional thermal solution.
Viewer-derived dynamic pressure uses its available velocity approximation;
use the verified study diagnostics for quantitative atmosphere-relative
pressure comparisons.

## Size and optional detail

```julia
export_visualization("output/simulation_results";
    max_frames=2000, data_budget_mb=150.0,
    texture_resolution="4k", trail_orbits=3, frame=:inertial)
```

The exporter decimates saved rows to the frame and data limits. The page
reports its resulting cadence. Smaller textures or `textures=false` reduce
the page size; no texture changes the recorded trajectory.

STL, OBJ and uncompressed GLB models can replace link boxes:

```julia
export_visualization("output/simulation_results";
    models=Dict(1 => "data/models/iss_nasa_3d_resources_b.glb"),
    model_scale=2.4, model_rotation_deg=Dict(1 => (-90, 0, -90)))
```

Match the spacecraft ID, scale and orientation to your model. The model is
display geometry; it does not replace the aerodynamic or structural model.
`data/models/README.md` and `data/textures/manifest.toml` record asset sources.

Atmosphere displays use density values already saved along the trajectory.
The sidecar also adds an altitude profile for the built-in exponential and
piecewise exponential models, whose evaluations have no model state.
For GRAM, surrogate, empirical and user-defined models, export adds the
entry-interface shell without making extra atmosphere calls. It does not
reinitialize the atmosphere, change its epoch or advance its sampler.

## Ensembles

```julia
result, page = run_monte_carlo_visualization(
    seed -> build_args(seed), 1:20, "output/campaign"; threads=4)
```

Supply your existing `build_args(seed)` configuration builder. Each sample
writes to its own results and checkpoint directories, including when the
builder supplies a shared checkpoint root. The page compares the trajectories and colors
them by final spherical periapsis altitude. Supply `scalar` and `scalar_name`
for another quantity. An existing campaign can be exported with
`export_ensemble_visualization("output/campaign")`.

## Standalone files

```text
python3 viewer/build_standalone.py --out output/viewer.html
```

Open the page and choose a CSV or JSON file, body and epoch. CSV accepts
SpaceAGORA result columns or `time,x,y,z` with optional velocity and
quaternion columns. Check units and the inertial-frame convention before
interpreting the view. The standalone rotation is an approximation; it does
not replace the frame metadata saved by a simulation.

The default page embeds its renderer and assets and makes no network
requests. Imported models with external resource URLs can request them. For
a host that refuses embedded `data:` scripts, `viewer/build_cdn_page.py`
rewrites an exported page to load the pinned three.js and mp4-muxer from
`cdn.jsdelivr.net` with the renderer modules inlined; the built page lists
its network hosts and third-party notices in its footer, and adds fonts
from Google Fonts only with `--fonts`. See `viewer/README.md`.

See `viewer/README.md` for detailed controls, renderer development and
licensing. Terrain payload export, landing dynamics and plume-surface
physics are separate work and are not part of this viewer integration.
