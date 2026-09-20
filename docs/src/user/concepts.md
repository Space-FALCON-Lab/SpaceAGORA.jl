# Concepts

Use this page when you want the mental model behind the main workflows without
digging into maintainer-only architecture documents.

This page is for users who can already run commands and now want to understand
which route, mode, or surface is appropriate for their work.

Shortest successful command:

```text
julia --project=. src/cli/main.jl run --example=AGORA_Basic_Quickstart.jl --output-dir=output/cli_run
```

What to read next:

- [First Simulation](first_simulation.md)
- [Assets & Modes](../assets.md)
- [Public API](../generated/public_api.md)

## Two operating modes

### No-GRAM mode

Use no-GRAM mode when you want:

- onboarding
- smoke runs
- lightweight studies
- the shortest route to a successful run on a new machine

### Higher-fidelity mode

Use the higher-fidelity path when you need:

- mission-quality atmosphere and winds
- SPICE-backed geometry and frame data
- other machine-local or licensed assets

## Three common entry surfaces

### Example scripts

Example scripts are the quickest way to run a repository-owned scenario.

### CLI

The CLI is the stable packaged command surface for examples, verification,
benchmarks, and asset inspection.

### Root package API

For scripts and extension work, depend on the exported `SpaceAGORA` surface.
Treat the generated Public API page as the supported list of exports.

## Outputs stay local

SpaceAGORA regenerates reports, docs builds, and similar artifacts into ignored
paths such as `output/`, `docs/build/`, `docs/site/`, and `docs/src/generated/`.

## Performance knob: `isolate_state`

`SpaceAGORA.run_simulation(...; isolate_state=true)` deep-copies the
configuration by default so repeated or concurrent runs do not alias mutable
state.

Treat `isolate_state=false` as an expert-only throughput lever after you have
measured the setup cost and confirmed that the configuration instance is not
shared across overlapping or state-mutating runs.

## Query a regional terrain map

Use a terrain model to query surface height at a latitude and longitude. These
queries are independent of propagation: adding a terrain model does not enable
landing guidance, a touchdown event, a radar altimeter or viewer terrain.

Coordinates are **planetocentric latitude and east-positive longitude in
degrees**. Heights are metres above a stated reference sphere, not geodetic
heights above an ellipsoid. Keep the same reference radius for every grid in a
model.

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
terrain_radius(terrain, 1.5, 12.0)       # 3_396_415.0 m from the body centre
dem_grid_covers(grid, 1.5, 12.0)        # true
terrain_height(terrain, 5.0, 12.0)       # -5.0, the selected fallback
```

Grid bounds describe the **outer cell edges**, and samples sit at cell centres.
The query interpolates between centres and clamps to the nearest centre at an
outer edge. The example's first sample is at latitude 2.5°, longitude 10.5°.
When grids overlap, the first matching grid wins; put the preferred grid first.
`DEMGrid` copies its height matrix, and `DEMTerrainModel` makes independent copies
of its grids. Changing the input arrays later does not change the model.

Regional grids may cross the longitude seam: use increasing bounds such as
`350.0, 370.0`; queries at −5° and 355° then agree. Longitudes are compared modulo
360°. Latitude must be within −90° to 90°, and a grid's longitude span must be
greater than zero and less than 360°. Full-globe grids are not supported.
Samples, bounds and queries must be finite; missing-height values such as NaN
are rejected. Each grid needs at least two rows and two columns.

An in-memory grid may omit `reference_radius_m` when its datum is supplied by
the model. If the grid records a radius, it must match the model's explicit
positive radius. `terrain_radius(terrain, lat, lon, radius_m)` also rejects a
conflicting radius. For a flat sphere, use
`terrain_radius(NoTerrainModel(), lat, lon, radius_m)`.

### Load local grid files

`load_dem_grid("region.json")` reads a JSON header and the neighbouring
`region.f32`. For the matrix above, the header is:

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

The binary contains exactly `rows × cols` little-endian Float32 samples, one row
after another from north to south, each row ordered west to east. This is not
Julia's default matrix storage order. The loader rejects missing samples, extra
bytes, invalid dimensions, nonfinite heights and missing reference radii. It
reads local files only; it does not download maps.

For ordered grids and a named query point, `load_site_terrain("site.json")`
returns `(terrain, site)`. A minimal site file is:

```json
{
  "site": {"name": "Example", "lat_deg": 1.5, "lon_deg": 12.0},
  "dem": [{"name": "region", "reference_radius_m": 3396200.0}]
}
```

Every `dem` entry names a neighbouring grid header. Every grid declares the same
radius; no default planetary radius is assumed. The returned `site` contains the
point's coordinates, name, height and reference radius.
