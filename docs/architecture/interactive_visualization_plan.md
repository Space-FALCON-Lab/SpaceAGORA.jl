# Interactive Visualization Plan

Status: design record, agreed 2026-09-14. All five phases plus phase 6
(atmosphere) implemented 2026-09-14: scene layer, sidecar and `link_pose` save field (1); viewer MVP
under `viewer/`, single-file bundler, textures for the five launch bodies,
`run_simulation(args; visualization=true)` (2); close-up assemblies with
attitude and link poses, thruster/facet glyphs, follow camera, picking and
info panel, STL override (3); ensembles from Monte Carlo campaigns and
constellation members with a common time axis, scalar colouring, sample
selector and full-history view (4); `spaceagora visualize`, `run --visualize`,
the docs page, `assets check` reporting the textures, and 8k texture tiers
for Earth, Mars and the Moon (5). The hillshade fallback from topography
harmonics was not built: every launch body has a texture, and an
unregistered body draws a flat colour with a graticule.

## Goal

An interactive, post-hoc 3D view of a SpaceAGORA run: the central body with a
detailed surface map, every spacecraft in the run, and, when the camera is
close enough, a 3D model of each spacecraft assembled from the link dimensions
in the simulation configuration. It never runs inside the integrator loop.

## Decisions

| Question | Decision |
|---|---|
| Viewer technology | Standalone three.js viewer, one self-contained HTML file per run |
| Scale | Single spacecraft, constellations (thousands of members), and Monte Carlo ensembles |
| Bodies | Earth, Mars, Venus, Titan, Moon at launch; adding a body is a texture plus one manifest entry |
| Default frame | Inertial (J2000) camera with the globe rotating underneath |
| Articulated links | Joint deflections are shown; a new save field records per-link poses |
| Attitude without orientation simulation | Model is aligned with the velocity vector |
| Texture assets | Committed to the repo under `data/textures/` |
| When it runs | Only when explicitly flagged; the default run writes nothing extra |

## What the repo already provides

- `Link.dims` is a box (x, y, z) and every non-root link carries `r` and `q`
  relative to the root bus. `Joint` records attachment points on both links.
  No new geometry schema is needed.
- Runs write an Arrow file (plus optional CSV) with per-spacecraft columns
  `sc{i}_pos_{1..3}`, `sc{i}_vel_{1..3}`, altitude, latitude, longitude, mass,
  and `sc{i}_q_{1..4}` when `orientation_sim` is on. This is the trajectory
  feed; nothing is duplicated.
- `r_intor_p!` and `_j2000_to_body_fixed_state` give the J2000 to body-fixed
  rotation for both the SPICE and simple ephemerides paths. The viewer never
  needs SPICE; the exporter samples the rotation at playback cadence.
- The RPO station loader already parses STL, so CAD overrides are cheap.
- Topography harmonics exist for Earth, Mars, and Venus and can drive a
  hillshade fallback for a body with no texture.

Non-root link poses are not part of the integrated state on the standard
path. They live in the mutable `Link` objects and are changed by control code
through `rotate_link` (solar-panel articulation, aerobraking attitude
executors). Only the cloth/robot-arm coupled path integrates link states
(`arm_r`, `arm_q` in the state shape). Showing deflections therefore means
snapshotting link poses at each saved step, not reading the state vector.

## Architecture

Two layers with a documented JSON contract between them.

### Layer 1: scene layer (Julia, no heavy dependencies)

Module `SceneVisualization` under `src/analysis/visualization/scene/`,
included from `core/simulation_model.jl` before the callbacks so the engine
and the save fields can reach it:

| File | Owns |
|---|---|
| `scene_visualization.jl` | Module, `VisualizationScene`, `PlanetSpec`, `SpacecraftGeometry`, `LinkBox` and glyph types |
| `spacecraft_geometry.jl` | `SpacecraftModel` to list of placed boxes, thruster/facet/joint glyphs, `link_pose` layout, velocity-aligned attitude; optional STL reference |
| `planet_spec.jl` | Planet radius/flattening/spin, texture lookup from the texture manifest, sampled body-fixed quaternion table |
| `scene_export.jl` | Scene assembly from a configuration, `scene.json` writer and reader, frame budget arithmetic |
| `viewer_bundle.jl` (phase 2) | Assembles the single-file HTML: viewer JS, textures, and binary trajectory data |

Public API on root `SpaceAGORA`:

- `export_visualization(prefix; kwargs...)`: reads `<prefix>_scene.json` and
  the Arrow file it names, writes `<prefix>_viewer.html`.
- `export_visualization(args::SimulationConfiguration; kwargs...)`: same, for
  the results bundle of `args`.
- `with_visualization_scene(args, flag)`: copy of a configuration with the
  sidecar flag set.
- `write_viewer_dev_payload(prefix, out_js)`: payload for `viewer/dev.html`.
- `export_monte_carlo_visualization(campaign_dir; kwargs...)`: ensemble mode
  (phase 4).

Keyword arguments: `max_frames`, `data_budget_mb`, `trail_s`, `frame`
(`:inertial` or `:planet_fixed` initial camera), `speed`, `title`,
`textures`. Texture resolution is fixed at 4096 x 2048 by the shipped assets.

Serialization uses `JSON` (already in the Manifest through PlotlyJS, now a
direct dependency). `ci_stale_deps_gate` requires it to be used, which it is.

### Layer 2: viewer (three.js, top-level `viewer/`)

| Path | Owns |
|---|---|
| `viewer/vendor/` | Pinned `three.module.js`, `OrbitControls.js` (MIT, license file alongside) |
| `viewer/src/` | Plain ES modules: `globe.js`, `spacecraft.js`, `lod.js`, `timeline.js`, `ensemble.js`, `ui.js`, `main.js` |
| `viewer/template.html` | Page skeleton with placeholders the Julia bundler fills |

No Node toolchain. The Julia bundler (`viewer_bundle.jl`) puts the vendored
modules and every viewer module into the page's import map as `data:` URLs
(bare specifiers `three`, `three/addons/...`, `viewer/<name>.js`), so the
modules keep real `import` statements and `viewer/dev.html` can map the same
specifiers to local files for development. Textures are base64 data URIs.
Trajectory data is Float32 little-endian base64 decoded into typed arrays on
load. Phase 2 ships `data.js`, `globe.js`, `spacecraft.js`, `timeline.js`,
`ui.js`, `main.js`; `lod.js` and `ensemble.js` come with phases 3 and 4.

## Data contract

### `simulation_results_scene.json` (written by the run when flagged)

```
{
  "schema": 1,
  "epoch": {"et_start_s": ..., "utc": "2014-05-27T05:00:00.000Z"},
  "planet": {
    "name": "Mars", "equatorial_radius_m": ..., "polar_radius_m": ...,
    "spin_rad_s": [0, 0, 7.088e-5], "inertial_frame": "J2000", "texture": "mars",
    "rotation": {"t_s": [...], "q_pi": [[x, y, z, w], ...]}
  },
  "spacecraft": [
    {"id": 1, "name": "sc1", "stl_path": null, "bounding_radius_m": ...,
     "links": [
       {"name": "root",  "root": true,  "dims_m": [x, y, z], "r_m": [0, 0, 0], "q": [0, 0, 0, 1], "mass_kg": ...},
       {"name": "link1", "root": false, "dims_m": [..], "r_m": [..], "q": [..], "mass_kg": ...}
     ],
     "thrusters": [{"link": 1, "location_m": [..], "direction": [..], "max_thrust_n": ..}],
     "facets":    [{"link": 1, "name": "..", "normal": [..], "cp_m": [..], "area_m2": ..}],
     "joints":    [{"link1": 1, "link2": 2, "p1_m": [..], "p2_m": [..]}]
    }
  ],
  "orientation_sim": false,
  "results": {
    "feather": "simulation_results.feather",
    "link_pose": {"field": "link_pose", "stride": 7,
                  "layout": ["rx","ry","rz","qx","qy","qz","qw"],
                  "columns": "sc{i}_link_pose_{1..7n}"}
  }
}
```

Link indices in `thrusters`, `facets` and `joints` are 1-based into `links`
(root is 1). The texture key is `lowercase(planet.name)` and is resolved
against the texture manifest by the viewer bundler (phase 2).

The rotation table is sampled at `data_rate` from the same transform the
lat/lon save fields use, so ground tracks in the viewer match the CSV.

### New save field `link_pose`

Per spacecraft, the non-root links concatenated: 7 floats each (`r` relative
to root, `q` relative to root, scalar-last), columns
`sc{i}_link_pose_{1..7n}` in the link order the sidecar records
(`SpacecraftGeometry.links[2:end]`). Added by `default_save_fields` only when
the visualization flag is on and at least one spacecraft has non-root links.
The cloth/robot-arm coupled path integrates a separate arm chain
(`arm_r`/`arm_q`) that is not a spacecraft `Link`; rendering that chain and
cloth nodes is out of scope for the first release.

### Trajectory encoding in the HTML

`window.SPACEAGORA_VIEWER = {scene, frames, textures, models, options}`.
`frames` holds `t_s`, `pos_km` (row-major frame, spacecraft, xyz; Float64
when `pos_dtype` is `"f64"`, which the bundler chooses for runs of at most
64 spacecraft so metre-scale models do not jitter, Float32 otherwise),
optional `vel_kms`, `q` and `mass_kg`, and optional `link_pose` (`stride`,
per-spacecraft `counts` and `offsets`, `data`), each a base64 little-endian
block, plus `count`, `sats`, `source_rows`, `stride_rows`. `models` maps a
spacecraft id to an embedded STL (`url`, `scale`). Positions in km relative
to planet center. The exporter decimates to
`max_frames` (default 2000) and then to `data_budget_mb` (default 150 MB),
whichever is smaller, with linear position interpolation and quaternion slerp
in the viewer. At 32768 members a 150 MB budget yields roughly 380 frames,
enough for playback of a constellation but not for close-up attitude work;
the viewer shows the effective cadence in its info panel.

### Monte Carlo ensembles

`run_monte_carlo_visualization(build_args, seeds, campaign_dir)` runs one
sample per seed through `run_monte_carlo`, each writing its own bundle and
sidecar under `<campaign_dir>/sample_NNNN`, then writes
`ensemble_manifest.json` (index, seed, success, scalar, label, directory)
with a scalar per sample from the caller's `scalar(df)` (default: final
periapsis altitude in km) and builds `ensemble_viewer.html`.
`export_ensemble_visualization(dir)` merges any `sample_NNNN` or
`sat_<i>_id_<id>` subdirectories (the latter is what
`run_constellation_ensemble` writes), resampling every sample onto one time
axis with `NaN` outside its own span; sample `i`, spacecraft `k` becomes
pseudo-spacecraft `(i-1)*per_sample + k` and the page gains
`payload.ensemble`. The viewer colours samples by the scalar (viridis),
draws every history as a faint line, ringing the selected sample and
dimming the rest, and the markers at the scrubbed time are the dispersion
cloud.

## Flagging

- `SimulationSettings.save_visualization_scene::Bool = false` next to the
  existing `save_csv`. When true the run writes
  `<results_directory>/simulation_results_scene.json` after the results
  bundle and enables the `link_pose` save field. It does not build the HTML.
- `run_simulation(args; visualization=true)` sets the flag and builds the
  HTML at the end of the run; `SPACEAGORA_VISUALIZATION=1` in the environment
  is the same switch for unmodified example scripts.
- CLI: `spaceagora run --example=<file> --visualize` and
  `spaceagora visualize --run=<prefix> [--max-frames=N] [--frame=planet_fixed]`.

## Rendering behavior

- Globe: oblate sphere from `Rp_e`/`Rp_p`, equirectangular texture, rotated by
  the sampled body-fixed quaternion. Toggle to planet-fixed camera. Optional
  hillshade from topography harmonics when no texture exists.
- Far view: one point cloud for all spacecraft markers (a small shader with a
  per-point hidden/selected state), so constellations of tens of thousands
  render in one draw call. Trails as fading line segments for runs with at
  most 64 spacecraft, sized in orbits (default three) from the period
  estimated off periapsis passages. Labels are sprites, same threshold.
- Near view: when camera distance to the selected spacecraft drops below a
  multiple of its bounding radius, the marker is replaced by a group of
  `BoxGeometry` meshes, one per link, placed from `link_pose` when saved and
  from the configured pose otherwise. Thrusters render as small cones and
  facets as translucent quads, both toggleable. A per-spacecraft STL
  override replaces the boxes when present.
- Attitude: from `sc{i}_q` when orientation simulation is on. Otherwise body
  +x is aligned with the inertial velocity and body +z with the nadir
  direction projected orthogonal to velocity. That matches the flow-along-+x
  convention the aerodynamic reference areas assume.
- Time: play/pause, speed, scrub, wall-clock and elapsed readouts, frame
  interpolation.
- Camera: free orbit about the planet, follow-spacecraft, planet-fixed.
- Picking: click a marker for id, altitude, latitude, longitude, speed, mass.

## Textures

`data/textures/<body>_<source>_<res>.jpg` with `data/textures/manifest.toml`
recording body name, resolution tier, source, license, and the east
longitude of the image's left edge (`lon_left_deg`; Titan's USGS mosaic is 0
to 360 E, the others are 180 W to 180 E). `scripts/dev/build_viewer_textures.py`
rebuilds any tier from the sources. Earth, Mars and the Moon also ship 8192 x
4096 tiers (4 to 7 MB), picked by default; the viewer downscales on load
when the GPU limit is smaller. All five launch bodies ship as 4096 x 2048
JPEGs, 1.1 to 2.3 MB each: Blue Marble Next Generation (Earth), Viking MDIM 2.1 via
NASA Trek tiles (Mars), Magellan C3-MDIR colorized topography (Venus),
Cassini ISS P19658 (Titan), LRO WAC via the NASA SVS CGI Moon Kit (Moon).
`.gitignore` allowlists the directory and `data/assets_manifest.toml` carries
a `visualization_textures` asset; `assets check` reporting is phase 5. Adding
a body means dropping a file in and adding one manifest entry.

## Contract and gate updates

- `docs/architecture/canonical_topology_contract.md` and
  `src_completeness_contract.md`: add `src/analysis/visualization/scene/*`
  as a canonical owner and top-level `viewer/*` as the viewer root, with
  `src/viewer/*` forbidden.
- `ci_path_policy_gate.jl` and `ci_src_completeness_contract_gate.jl`: extend
  the inventories to match.
- `ci_public_api_surface_gate.jl` and `ci_docs_contract_gate.jl`: the new
  exports need docstrings and a docs page.
- `ci_no_artifact_files_gate.jl` already scans tracked HTML; the template
  must stay free of the tokens it forbids.
- Project.toml: `JSON3` added under deps and compat.

## Tests

- `test/unit/visualization/`: link tree to box list (positions, orientations,
  bounding radius), planet rotation table equals `r_intor_p!` at sampled
  times for both ephemerides paths, decimation and budget arithmetic, JSON
  round trip, velocity-aligned attitude frame orthonormality, HTML bundle
  contains every expected block.
- Integration (`test/integration/`): run the quickstart with the flag, assert
  the sidecar and viewer HTML exist and the embedded data length matches the
  Arrow row count after decimation. No headless browser in CI; the viewer JS
  is exercised manually and by a syntax-only check if a Node binary is found.
- Existing suites: flag off by default, so no golden output changes.

## Phases

1. Scene layer and sidecar: types, geometry conversion, planet spec with
   rotation sampling, JSON writer, `link_pose` save field, settings flag,
   unit tests, contract and gate updates.
2. Viewer MVP: globe with texture and rotation, marker cloud, trails,
   timeline, free camera with an inertial/planet-fixed toggle, single-file
   bundler, end-to-end test (done; markers are a `Points` cloud rather than
   an instanced mesh, one draw call either way).
3. Level of detail: box assemblies with attitude and `link_pose`, thruster and
   facet glyphs, follow camera, picking and info panel, STL override (done;
   visibility switches on the projected size of the bounding radius with
   hysteresis, follow mode uses a floating origin, trails are measured in
   orbits from the estimated period and fade with age).
4. Constellations and ensembles: data budget decimation, trail thresholds,
   Monte Carlo campaign wrapper and manifest, ensemble overlay and selector
   (done).
5. Polish: CLI `visualize` subcommand and `run --visualize`, docs page,
   `assets check` integration, 8k texture tiers with a GPU-limit fallback
   (done; hillshade fallback dropped, see status).

## Phase 6: atmosphere (added 2026-09-14)

The sidecar gains an optional `atmosphere` block: model type name, the
entry-interface altitude, a density profile from the surface to 1.5 x EI
sampled through the run's own `getDensity` (the engine passes its integrator
parameters, so the 7-argument GRAM and NRLMSISE-00 forms work), and, for
models that vary horizontally, a latitude/longitude density map at 0.6 x EI
on a 5 degree grid. A `density` save field joins `link_pose` under the flag;
the bundler embeds `density_kg_m3`, `heat_rate_w_m2`, `drag_n` (magnitude)
and `wind_ms` blocks. The viewer draws a Fresnel limb glow at the EI,
density shells with opacity from the profile, and the draped map; trails
colour by heat rate (default when saved), dynamic pressure, density,
altitude or speed on an inferno scale, and the selection panel reads the
same quantities. A wind arrow was left out because the saved wind vector's
frame is not recorded; the panel shows wind speed only.

## Robot arms (added 2026-09-14)

The cloth robot-arm chain, left out of phase 1, is now rendered. Geometry:
`SpacecraftGeometry.arm` (`ArmGeometry`: links with vector, COM offset,
radius, mass; mount offset; reach) from the `RobotArmPlan` a
`RobotArmControlEffector` carries for that spacecraft, found with the same
lookup the engine uses to couple the state. Poses: the `arm_pose` save field
writes `sc{i}_arm_pose_{1..7n}`, each arm link's COM relative to the
spacecraft (inertial, metres) and inertial quaternion from `arm_r`/`arm_q`,
and the bundler embeds it as an `arm_pose` block beside `link_pose`. The
viewer draws cylinders per link in an unrotated sibling of the assembly.
Cloth nodes remain out of scope.

## Out of scope for the first release

Real-time streaming from a running integrator, cloth node rendering,
atmosphere or density overlays, video export, and multi-body scenes with more
than one central body.
