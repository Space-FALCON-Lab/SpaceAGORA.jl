# Concepts

This page gives the mental model behind the main workflows: which operating
mode to start in, which entry point to use, and where outputs go. It assumes
you can already run commands from the repository root; the contributor
documentation in the repository is not needed.

Once the examples environment is prepared
(`julia --startup-file=no examples/setup.jl`), the shortest complete run is:

```text
julia --project=. src/interfaces/cli/main.jl run --example=AGORA_Basic_Quickstart.jl --output-dir=output/cli_run
```

What to read next:

- [Your First Script](first_script.md)
- [Assets & Modes](../assets.md)
- [Public API](../generated/public_api.md)

## Two operating modes

### No-GRAM mode

No-GRAM mode needs nothing beyond the public checkout: no GRAM installation,
no SPICE kernels and no licensed data. Use it to learn the package, for smoke
runs and lightweight studies, and for the first run on a new machine.

### Higher-fidelity mode

Switch to the higher-fidelity path when your models need GRAM atmosphere and
winds, SPICE-backed geometry and frames, or other machine-local or licensed
assets. Each of these needs its own preparation; [Assets & Modes](../assets.md)
lists them and the command that checks which ones your machine has.

## Three common entry surfaces

### Example scripts

The scripts in `examples/` are complete scenarios owned by the repository and
the quickest way to run one. They run in the examples environment.

### CLI

`src/interfaces/cli/main.jl` is the packaged command surface. Its commands run
examples (`run`), the telemetry verification study (`telemetry`), the benchmark
studies (`benchmark`), the asset inventory (`assets`) and the viewer
(`visualize`); see [CLI](../cli.md).

### Root package API

Scripts and extensions should depend on what `SpaceAGORA` itself exports. The
generated Public API page is the list of supported exports.

## Outputs stay local

Reports, docs builds and similar generated files go to paths that Git ignores:
`output/`, `docs/build/`, `docs/site/` and `docs/src/generated/`.

## Performance knob: `isolate_state`

By default, `run_simulation(config; isolate_state=true)` works on a deep copy
of the configuration, so repeated or concurrent runs never share mutable
model state.

`isolate_state=false` skips that copy, and the run then advances your own
configuration objects. Use it as a speed setting only after you have measured
that the copy matters for your setup, and only when no other run uses the same
configuration object at the same time. Some models document it as the way to
read back what they recorded during a run, such as a command log. The RPO
planner lifecycle refuses `isolate_state=false` at startup.

## Query a regional terrain map

A terrain model returns surface height at a latitude and longitude. Queries
are independent of propagation: adding a terrain model does not turn on
landing guidance, a touchdown event, a radar altimeter or terrain in the
viewer.

Coordinates are planetocentric latitude and east-positive longitude, in
degrees. Heights are in meters above a stated reference sphere, not geodetic
heights above an ellipsoid, and every grid in a model uses that one reference
radius.

```julia
using SpaceAGORA

# Rows run north to south; columns run west to east.
heights = Float32[100 110 120 130; 200 210 220 230; 300 310 320 330]
radius_m = 3_396_200.0
grid = DEMGrid(heights, 0.0, 3.0, 10.0, 14.0;
               reference_radius_m=radius_m)
terrain = DEMTerrainModel([grid]; reference_radius_m=radius_m,
                          fallback_height_m=-5.0)

terrain_height(terrain, 1.5, 12.0)       # 215.0 m above the sphere
terrain_radius(terrain, 1.5, 12.0)       # 3_396_415.0 m from the body center
dem_grid_covers(grid, 1.5, 12.0)        # true
terrain_height(terrain, 5.0, 12.0)       # -5.0, the selected fallback
```

The grid bounds are the outer cell edges, and each sample sits at a cell
center; in the example the first sample is at latitude 2.5° and longitude
10.5°. A query interpolates bilinearly between centers; in the half cell along
an outer edge it clamps to the outermost row or column of centers. A point
that no grid covers gets the model's `fallback_height_m` (zero by default).
Where grids overlap, the first grid that covers the point wins, so list the
preferred grid first.
`DEMGrid` copies its height matrix and `DEMTerrainModel` copies its grids, so
changing the input arrays afterwards does not change the model.

A regional grid may cross the longitude seam: give increasing bounds such as
`350.0, 370.0`, and queries at −5° and 355° then return the same height.
Longitudes are compared modulo 360°. Latitudes must lie in −90° to 90°, and a
grid's longitude span must be greater than zero and less than 360°; full-globe
grids are not supported. Samples, bounds and query coordinates must be finite,
so missing-height values such as NaN are refused, and each grid needs at least
two rows and two columns.

An in-memory grid may leave out `reference_radius_m` and take its datum from
the model. A grid that does record a radius must match the model's radius
exactly, and `terrain_radius(terrain, lat, lon, radius_m)` refuses a radius
that differs from the model's. For a bare sphere with no terrain, use
`terrain_radius(NoTerrainModel(), lat, lon, radius_m)`. A model type of your
own can subtype `AbstractTerrainModel` and add a `terrain_height` method; the
[Terrain API](../generated/terrain_api.md) lists the functions.

### Load local grid files

`load_dem_grid("region.json")` reads a JSON header and the binary file beside
it, `region.f32`. For the matrix above, the header is:

```json
{
  "rows": 3,
  "cols": 4,
  "lat_min": 0.0,
  "lat_max": 3.0,
  "lon_min": 10.0,
  "lon_max": 14.0,
  "reference_radius_m": 3396200.0,
  "units": "m",
  "layout": "row-major, north to south, west to east, little-endian float32"
}
```

`rows`, `cols`, the four bounds and `reference_radius_m` are required; `units`,
`layout` and `format` are optional, but when present they must agree with the
format. The binary holds exactly `rows × cols` little-endian Float32 samples,
one row after another from north to south, each row ordered west to east. This
is not Julia's column-major matrix order. The loader refuses a file with too
few or too many bytes, invalid dimensions, nonfinite heights or no reference
radius. It reads local files only and never downloads maps.

For several grids and a named query point, `load_site_terrain("site.json")`
returns `(terrain, site)`. A minimal site file is:

```json
{
  "site": {"name": "Example", "lat_deg": 1.5, "lon_deg": 12.0},
  "dem": [{"name": "region", "reference_radius_m": 3396200.0}]
}
```

The `dem` entries are in priority order, and each names a grid header in the
same folder (here `region.json`). Every entry and every header must state the
same radius; no default planetary radius is assumed. The returned `site` holds
the point's `lat_deg`, `lon_deg`, `name`, its terrain `height_m` and the
`reference_radius_m`.
