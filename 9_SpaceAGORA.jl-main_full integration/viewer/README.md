# The viewer has moved

The renderer and its build tools now live in the independent
[SpaceAGORAVisualizer.jl repository](https://github.com/Space-FALCON-Lab/SpaceAGORAVisualizer.jl).
That package retains the original renderer attribution and third-party notices.

SpaceAGORA owns saved simulation data, scene metadata and portable viewer JSON.
The optional package turns that data into a self-contained HTML page.

From the SpaceAGORA repository root, install the tested pair in a separate project:

```sh
julia --project=. scripts/setup_visualizer.jl .visualizer-env
julia --project=.visualizer-env
```

```julia
using SpaceAGORA, SpaceAGORAVisualizer
payload = export_viewer_payload("output/my-run/simulation_results"; textures=false)
SpaceAGORAVisualizer.visualize(payload)
```

For a new mission, see the [mission workflow](../docs/src/user/mission_workflows.md).
For existing exports and display options, see the [visualization guide](../docs/src/user/visualization.md).
