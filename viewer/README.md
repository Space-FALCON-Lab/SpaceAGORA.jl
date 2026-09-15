# SpaceAGORA viewer

Browser viewer for a finished run: the central body with its surface texture,
every spacecraft as a marker with a trail, a playback timeline, and an
inertial or planet-fixed camera. It never runs inside the integrator; it reads
the results bundle and the scene sidecar written when
`simulation_settings.save_visualization_scene = true`.

Design record: `docs/architecture/interactive_visualization_plan.md`.

## Producing a page

```julia
using SpaceAGORA
run_simulation(args; visualization=true)            # sidecar + viewer page in one go
export_visualization("output/simulation_results")   # or later, from an existing bundle
```

`export_visualization` writes `output/simulation_results_viewer.html`, a
self-contained file: three.js, these modules, the texture and the trajectory
are all embedded, so it opens from disk with no server and no network.
Keywords: `max_frames`, `data_budget_mb`, `trail_orbits` or `trail_s`, `frame`
(`:inertial` or `:planet_fixed`), `speed`, `title`, `textures`, `stl`, `stl_scale`.

## Layout

| Path | Owns |
|---|---|
| `vendor/three.module.js`, `vendor/OrbitControls.js`, `vendor/STLLoader.js`, `vendor/OBJLoader.js`, `vendor/GLTFLoader.js`, `vendor/BufferGeometryUtils.js` | three.js r160 (MIT, `LICENSE.three`); GLTFLoader's one relative import is rewritten to the bare `three/addons/utils/...` specifier so it resolves from a `data:` URL |
| `src/data.js` | payload decoding, frame interpolation, rotation table lookup |
| `src/globe.js` | textured oblate sphere, graticule, rotation with the body |
| `src/spacecraft.js` | marker cloud, fading trails, labels, orbit period estimate |
| `src/lod.js` | close-up box assemblies, thruster/facet glyphs, STL override, visibility by projected size |
| `src/ensemble.js` | ensemble coloring, spaghetti histories, sample selector |
| `src/paths.js` | reference polylines in inertial, RTN or body frames |
| `src/references.js` | reference ghosts: translucent copies of a spacecraft driven by an external state table |
| `src/video.js` | MP4 export through WebCodecs and mp4-muxer (`vendor/mp4-muxer.mjs`, MIT), WebM fallback |
| `src/terrain.js` | landing-site terrain: DEM-displaced nested patches draped with imagery that sharpens toward the site, and the hole they cut in the globe |
| `src/dust.js` | regolith blown off the surface by a descent engine: a GPU particle sheet driven by the `frames.plume` block, a ground haze disk and a scour mark |
| `src/plumes.js` | thruster plumes: one additive, flickering cone pair per thruster, driven by the recorded firing levels, with the idle cone glyphs dimmed |
| `src/plots.js` | time-history plot panel: any quantity in the panels opens its history over the run, with the playback cursor and click-to-seek |
| `src/atmosphere.js` | limb glow, density shells, density map |
| `src/colormaps.js` | inferno and viridis |
| `src/timeline.js` | playback clock and time formatting |
| `src/ui.js` | overlay controls and info panel |
| `src/main.js` | scene setup and the render loop |
| `template.html` | page skeleton the Julia bundler fills |
| `dev.html` | development harness (see below) |

The Julia side is `src/analysis/visualization/scene/viewer_bundle.jl`.
Modules import each other through bare specifiers (`viewer/globe.js`,
`three`) that the page's import map resolves, to `data:` URLs in the bundled
file and to local paths in `dev.html`. No build tool is involved.

## Developing the viewer

1. Produce a payload from any flagged run:
   `write_viewer_dev_payload("output/simulation_results", "viewer/dev_data.js")`.
2. Serve this directory: `python3 -m http.server 8000` from `viewer/`.
3. Open `http://localhost:8000/dev.html`, edit files under `src/`, reload.

`dev_data.js` is ignored by git.

## Controls

Space plays and pauses, the arrow keys step by one speed unit, the slider
scrubs. Drag to orbit, wheel to zoom, right-drag to pan. "Planet-fixed"
counter-rotates the scene so the body stands still and orbits precess past it.

Click a marker to select a spacecraft: the panel at the top right shows its
altitude, geodetic latitude and longitude, radius, speed, mass (when the run
saved it) and whether it is currently drawn as a marker or a 3D model.
"Follow" (or F) keeps the selected spacecraft at the center and moves the
camera close enough to see the assembly; Esc deselects. The trail selector
sets the history drawn behind each spacecraft in orbits (the period is
estimated from periapsis passages in the saved trajectory, so a run shorter
than one orbit shows "longer than run" and the trail covers the whole run).

## Close-up models

When a spacecraft's bounding radius covers a few pixels the marker gives way
to its assembly: one box per link from `Link.dims`, placed by the saved
attitude quaternion or, for runs without orientation state, an attitude with
body +x along the velocity and body +z toward nadir. Non-root links follow
the recorded `link_pose` columns when the run had them, otherwise their
configured pose. Thrusters are orange cones with the apex at the thruster
location pointing along the thrust direction, facets are translucent squares
of the facet area with their normal drawn, and the body axes are red, green
and blue for x, y, z. Each of these has a toggle. A 3D model passed through
`export_visualization(...; models=Dict(id => path), model_scale=..., model_rotation_deg=...)`
(STL, OBJ, glTF or GLB, parsed in the browser) replaces the boxes for that
spacecraft; the glyphs stay.

Runs with at most 64 spacecraft embed positions as Float64 so the meter-scale
models do not jitter at planetary distances; follow mode also re-centers the
scene on the followed spacecraft for the same reason.

## Heating overlay

With density and velocity in the frames, `src/lod.js` can swap every heatable
mesh (link boxes, model meshes) to a shader that colors faces by
½ρV³ cos θ (inferno, log scale, three decades below the run's peak), θ from
the face normal and the airspeed (inertial velocity minus ω × r). Toggle
"heating"; the legend reads W/cm².

## Thruster plumes

`frames.thruster_level` (Float32, frame-major then spacecraft then that
spacecraft's thrusters, 0 to 1) with `frames.thruster_counts` drives
`src/plumes.js`. A plume hangs on the thruster's own link group, runs along the
thruster's `direction`, and takes its reach and color from the rated thrust:
long and orange-white for a main engine, short and blue-white for an attitude
jet. Reach and brightness follow the level through a fractional exponent so a
jet firing a fraction of a percent of its rating still reads; level 0 draws
nothing and dims the static cone glyph. Toggle "plumes"; the selection panel
lists a `thruster k level` row per thruster.

## Ensembles

`payload.ensemble.nominal` (0-based sample index) marks the unperturbed
sample: white line, own label, never dimmed; `src/ensemble.js` computes the
per-time mean and RTN dispersion of the samples and draws a translucent
3σ tube (radial by cross-track) around the nominal, with outline rings and
a readout of the three 3σ values in the panel.

### Samples

A page built by `export_ensemble_visualization` (or
`run_monte_carlo_visualization`) carries `payload.ensemble`: the samples
become pseudo-spacecraft on one time axis, with `NaN` positions where a
sample has not started or has already ended. `src/ensemble.js` colors them
by the manifest scalar (viridis), draws every sample's full history as a
faint line, and adds a sample selector; the selected sample's marker is
ringed and the others dimmed. Assemblies are built for at most 256
spacecraft.

## Texture resolution

`data/textures/manifest.toml` may hold several tiers per body. The exporter
embeds the largest by default (`texture_resolution=:best`); the globe
downscales on load when the GPU's maximum texture size is smaller than the
image, and the info panel shows the tier and the GPU limit.

## Planned paths

`payload.paths` carries reference polylines (`src/paths.js`): points in km,
a frame (`inertial`, `rtn` of a target spacecraft, or that spacecraft's
`body` frame), color and dash style. RTN and body paths are rebuilt every
frame from the target's position, velocity or attitude, and the group sits
at the floating origin like the markers.

## Video export

`src/video.js` drives `viewer.renderAt(t)` (main.js exposes it along with
`setRecording` and `resizeTo`) frame by frame, encodes with `VideoEncoder`
(H.264) and muxes with mp4-muxer in memory, then triggers a download. The
render loop stands still while recording. Import map key `mp4-muxer`.

## Standalone page

`standalone.html` plus `src/standalone.js` and `build_standalone.py` make a
page that builds the payload in the browser from a CSV or JSON table, a body
name and an epoch: IAU rotation table, box spacecraft, model override,
reference ghost, and ensembles from several files. `data.js` accepts typed
arrays as well as base64 blocks for that reason. See the user guide.

## Reference ghosts

`payload.references` carries state tables (`src/references.js`): Float64
times and positions on their own grid, optional velocities and attitude
quaternions, the index of the spacecraft whose geometry the ghost copies,
color and opacity. Each ghost is a translucent copy of that spacecraft's
model override (through `loadModelObject`, shared with `lod.js`) or of its
link boxes, with a ring marker and label, a line over the whole reference
span, and the same pixel-size switch as the assemblies. The selection panel
of the flown spacecraft shows the separation from each ghost that refers to
it. Times outside the table hide the ghost.

## Robot arm

When `payload.scene.spacecraft[i].arm` is present and the frames carry an
`arm_pose` block, `src/lod.js` builds one cylinder per arm link (along the
link's own vector, joint sphere at its origin, a red tip at the end
effector) in a sibling group of the box assembly. Arm poses are inertial and
relative to the spacecraft in meters, so the group sits at the spacecraft
without the body rotation; each frame the joint origin is recovered from the
saved center of mass and quaternion.

## Atmosphere

`src/atmosphere.js` builds, from `payload.scene.atmosphere`, a limb-glow
shell at the entry interface (a Fresnel shader, additive, double-sided so it
reads from inside during a pass), translucent density shells whose opacity
follows the sampled profile, and a density map draped at one altitude when
the model varies horizontally. Trails can be colored by heat rate, dynamic
pressure, density, altitude or speed through inferno (`src/colormaps.js`);
the physical ones use a log scale over the run's range.

## Conventions

- Scene units are kilometers, J2000 axes.
- Quaternions are scalar-last. A sidecar or state quaternion `q` is the
  active body-to-inertial rotation, so it is applied to a three.js object
  directly. (In Julia `rot(q)` is the passive inertial-to-body matrix.)
- The sphere's texture seam is placed so that longitude 0 is on body +x; a
  texture whose left edge is not 180 W is shifted by its manifest
  `lon_left_deg`.
