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

The page is written beside the results as
`<results_directory>/simulation_results_viewer.html`, using
`args.simulation_settings.results_directory`.

For a completed run that already has its scene sidecar:

```julia
export_visualization("output/viewer-earth/simulation_results"; max_frames=2000)
SpaceAGORA.run_cli(["visualize", "--run=output/viewer-earth"])
```

`export_visualization` takes the results prefix:
`<results_directory>/simulation_results`, with no file extension. The CLI
`visualize --run=` accepts either that prefix or the results directory.

`with_visualization_scene(args, true)` sets the sidecar flag when preparing
a configuration. This saves the extra fields and sidecar; use
`export_visualization` afterwards to build the HTML page. Existing runs
without a sidecar need their configuration to create one, or use the
standalone data-file viewer below.

## Read the page

- Drag to orbit, scroll to zoom, and use the timeline to play or scrub.
- Click a spacecraft marker or visible model to select it. Its selection
  panel opens at the top right. Press F to follow it; Escape clears the selection.
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

## Plot another saved value

Click a quantity's name or value in the selection panel to plot its history.
Altitude, speed, mass, density, heat rate, drag and wind already have clickable
rows when the corresponding data is available.

To add the saved `sc1_periapsis_altitude` column, pass its suffix without the
spacecraft prefix. It contains osculating spherical periapsis altitude in metres:

```julia
export_visualization("output/viewer-earth/simulation_results";
    channels=[(column="periapsis_altitude", label="Spherical periapsis altitude",
               unit="m", digits=1, log=false)])
```

Reopen the generated page, select a spacecraft and click its **Spherical
periapsis altitude** row in the selection panel. Scrubbing the timeline updates
the displayed value. Channel names come from the saved results headers; see
[Simulation Outputs](outputs.md) for names, units and definitions.

Each spacecraft must have the corresponding column (`sc1_periapsis_altitude`,
`sc2_periapsis_altitude`, and so on). If any column is absent, the channel is
omitted without a message. If its row does not appear, check the spelling and
the headers for every spacecraft in the results table.

Missing values appear as gaps, and interpolation does not bridge a gap. Some
plots automatically use a logarithmic axis when their positive values span a
wide range. Zeros also appear as gaps on that axis; clear the **log** checkbox
to show them on a linear scale. Negative values disable logarithmic scaling.

If saved rows repeat a spacecraft's state after impact or touchdown while the
run continues, those held quantities plot as flat segments. The viewer uses the
saved rows and does not automatically end a history at deactivation.

Channel values use the same saved-row decimation as the trajectory. Choose
enough `max_frames` to retain the features you want to inspect.

## Size and optional detail

```julia
export_visualization("output/viewer-earth/simulation_results";
    max_frames=2000, data_budget_mb=150.0,
    texture_resolution="4k", trail_orbits=3, frame=:inertial)
```

The exporter decimates saved rows to the frame and data limits. The page
reports its resulting cadence. Smaller textures or `textures=false` reduce
the page size; no texture changes the recorded trajectory. Export defaults to
`texture_resolution="4k"`. Earth, Mars and Moon also include 8192 x 4096
textures: opt in with `texture_resolution="8k"` or select the largest available
tier with `:best`. Venus and Titan fall back to 4k. An 8k image has four times
the decoded pixels of 4k, and its embedded JPEG also makes the HTML larger.
The standalone builder likewise defaults to 4k; use `--textures 8k` for
a page containing the higher tiers.

STL, OBJ and uncompressed GLB models can replace link boxes:

```julia
export_visualization("output/viewer-earth/simulation_results";
    models=Dict(1 => "data/models/iss_nasa_3d_resources_b.glb"),
    model_scale=2.4, model_rotation_deg=Dict(1 => (-90, 0, -90)))
```

Match the spacecraft ID, scale and orientation to your model. The model is
display geometry; it does not replace the aerodynamic or structural model.
`data/models/README.md` and `data/textures/manifest.toml` record asset sources.
Bundled mission display meshes include Apollo Lunar Module, Cassini with and
without Huygens, CYGNSS, Magellan and Mars Odyssey, in addition to ISS. Supply
the corresponding GLB path through `models`; models are embedded only when
selected. Their model coordinates need mission-specific scale and rotation.
These six models and three higher-resolution textures add about 35.6 MB to
the checkout. No asset download or native atmosphere library is needed.

Atmosphere displays use density values already saved along the trajectory.
The sidecar automatically adds an altitude profile for the exponential and
piecewise exponential models and for the native-free `GRAMGridAtmosphereModel`.
These evaluations read the model without advancing a sampler or changing its epoch.

Fixed-grid profiles use the equator and longitude zero, over the overlap between
the grid's altitude coverage and zero to 1.5 times the entry-interface altitude.
A grid that does not cover the equator has no equatorial profile. The viewer
draws no density shells outside the sampled altitude range. A global density map
is included only when the grid covers every sampled latitude; regional grids
are not stretched over the globe. The default map height is 0.6 times the entry
interface, bounded by the grid's altitude range, and the sidecar records that
height. These displays describe the stored snapshot, not a newly generated
atmosphere at the simulation epoch.

Native GRAM, hybrid native/surrogate, empirical and custom models keep the
entry-interface shell and saved trajectory density without extra export-time
queries. Standalone `atmosphere_spec(args; sample_model=true, density_params=...)`
remains an explicit advanced operation that can change such a model's state;
scene export never opts into it.

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
licensing. Landing dynamics and plume-surface physics are separate from the viewer.

## Regional terrain in exported pages

For a completed Moon run with a scene sidecar, add a local terrain bundle:

```julia
export_visualization("output/moon-run/simulation_results";
    terrain="output/terrain/site.json", terrain_max_grid=512)
```

Reopen `output/moon-run/simulation_results_viewer.html`. The site bundle
contains `site.json`, the named DEM metadata and `.f32` samples, and optional
JPEG tiles. The Moon download tools and legacy imagery converter live in
`scripts/dev/terrain/`; `docs/reference/lunar_imagery_sources.md` explains their
Python requirements, commands, data sources and output layout. Downloading
assets is a separate step. Existing local bundles export without a network
connection or native GRAM.

The local site bundle follows `load_site_terrain`: regional DEMs in priority
order, explicit common reference radius, planetocentric latitude, east-positive
longitude and little-endian Float32 metre heights at cell centres. Export uses
the canonical loader, keeps the original outer edges, and samples reduced grids
uniformly at their new cell centres. `terrain_max_grid` must be an integer in
1:2048. The shared reduction factor can produce a singleton row or column. The
result approximates the source surface; `site.height_m` records the source query,
while markers and radar queries use the exported approximation.

The first covering grid supplies the height, with the canonical zero fallback
outside all grids. Longitude wraps without iteration. Radar clearance is radial:
`norm(position_body_m) - reference_radius_m - terrain_height_m`. The terrain
radius comes from the DEM, including for dust placement; scene planet metadata
is preserved. If the scene globe is an ellipsoid or uses another radius, the
terrain boundary can differ from that globe. The run and terrain must refer to
the same body. Use a consistent reference sphere to align their datums. Even
with matching radii, a visible rim is expected where the regional terrain meets
the reference globe at a different height. Export performs no datum conversion.

Select a spacecraft and read its **height above terrain** in the selection
panel. With a single spacecraft, **Follow** also selects it. This value is a
radial clearance based on the exported grid, not the saved geodetic altitude.
Zoom toward the site to inspect the local imagery mesh.

Imagery is optional. With no imagery, height queries and the height-above-terrain
readout remain available, while the normal globe is drawn. With a quadtree,
every declared local JPEG is embedded. Missing declared indices or tiles are
errors. Each tile has `level`, `x`, `y`, `file`, and positive `m_per_px`; the
root tile is required and sparse deeper levels inherit available ancestors.
The index uses a regional `root`, integer `tile_px`, and matching `max_level`.
Source, attribution and resolution metadata accompany the imagery.

To bound export size, the site accepts at most 64 DEMs and 128 MiB of DEM samples,
with 1 MiB per JSON document. Imagery accepts at most 8192 tiles through level
20, 4096 pixels per side, 8 MiB per JPEG, and 128 MiB total. Tile dimensions are
checked against their JPEG header. Paths, including symlinks, must resolve
inside the site or imagery directory. Export neither downloads terrain nor
changes propagation, touchdown or atmosphere state.
