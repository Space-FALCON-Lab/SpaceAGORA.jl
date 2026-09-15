# Interactive Visualization

Use this page to turn a finished run into a 3D page you can orbit, scrub and
share: the central body with its surface texture, every spacecraft, its
trail, and, up close, the spacecraft's links, thrusters and panels as they
moved. Nothing here runs inside the integrator; the viewer reads the results
bundle after the fact.

Shortest successful command:

```text
SPACEAGORA_VISUALIZATION=1 julia --project=. examples/AGORA_Earth.jl
# then open output/simulation_results_viewer.html
```

What to read next:

- [Simulation Outputs](outputs.md) for the sidecar and column layout
- `viewer/README.md` for the controls and how to develop the viewer

## Turning it on

The viewer is opt-in. Any of these produces the scene sidecar next to the
results bundle and builds the page:

```julia
run_simulation(args; visualization=true)                 # in a script
export_visualization("output/simulation_results")        # later, from an existing flagged run
```

```text
spaceagora run --example=AGORA_Earth.jl --visualize      # unmodified example
spaceagora visualize --run=output                        # page from an existing flagged run
SPACEAGORA_VISUALIZATION=1 julia --project=. my_script.jl
```

`with_visualization_scene(args, true)` returns a configuration with the
flag set when you build configurations by hand. Runs that pass their own
`save_fields` still get the link poses appended.

## What the page shows

- **Body.** An oblate sphere with the registered texture (Earth, Mars, Venus,
  Titan and the Moon ship at 4k; Earth, Mars and the Moon also at 8k, which
  the exporter picks by default). The globe rotates with the body's frame,
  so the ground track is where the latitude and longitude columns say.
- **Spacecraft.** Markers in the far view, with fading trails measured in
  orbits (default three; the period is estimated from periapsis passages)
  and labels for runs with at most 64 spacecraft. Click one to read its
  altitude, latitude, longitude, speed and mass; press F to follow it.
- **Close-up models.** Near a spacecraft the marker becomes one box per link
  from `Link.dims`, placed by the saved attitude quaternion or, for runs
  without orientation state, an attitude with body x along the velocity and
  z toward nadir. Panels follow the recorded `link_pose` columns, so an
  aerobraking pass shows the panels feathering. Thrusters are cones, facets
  translucent squares. `stl=Dict(id => "bus.stl")` swaps in a CAD mesh.
- **Frames.** Inertial by default; "Planet-fixed" holds the body still.

## CAD and other 3D models

Any spacecraft can be drawn from a 3D model instead of its link boxes:

```julia
export_visualization("output/simulation_results";
    models=Dict(1 => "data/models/iss_nasa_3d_resources_b.glb"),
    model_scale=2.4,                            # metres per model unit (or a Dict per id)
    model_rotation_deg=Dict(1 => (-90, 0, -90))) # XYZ Euler angles in the body frame
```

```text
spaceagora visualize --run=output --model=1=data/models/iss_nasa_3d_resources_b.glb --model-scale=2.4
```

STL, OBJ and glTF/GLB are supported and embedded in the page; a `.gltf` must
carry its buffers inline, and glTF files must not require Draco, meshopt or
KTX2 (the exporter refuses them and names the conversion command). By
default the model's bounding-box centre is placed on the spacecraft
(`model_center=false` keeps the file's own origin). The model sits in the
body frame at the root link, so attitude, link poses and glyphs still
apply, and the selection panel's "3D model" row reports the parsed mesh
count or the reason a model failed. `data/models/README.md` lists the
shipped models; the NASA ISS model there is public domain.

The same file can feed a planner. `sample_model_pointcloud(path; n_points,
scale, rotation_deg)` returns surface samples after the viewer's scale and
rotation, centred the same way, so an RPO station built from it (for
example through the `station_points` keyword of the CubeSat MPC demo
builder) is exactly the geometry the page draws:

```julia
points = sample_model_pointcloud("data/models/iss_nasa_3d_resources_b.glb";
    n_points=12000, scale=2.4, rotation_deg=(-90, 0, -90))
demo = build_rpo_cubesat_mpc_demo(; station_points=points, station_keepout_radius_m=3.0,
    start_rtn=SVector(-170.0, -100.0, 60.0), goal_rtn=SVector(70.0, 20.0, 0.0),
    search_margin_m=80.0, sample_ds_m=0.5, mission_time=900.0)
run_simulation(demo.args; visualization=true)
export_visualization(prefix; models=Dict(201 => iss), model_scale=2.4, model_rotation_deg=Dict(201 => (-90, 0, -90)),
    paths=[(name="HyPR plan", points_m=Matrix(demo.plan_result.path), frame=:rtn, target=2)])
```

`paths` overlays reference polylines: a planned path in a target's
radial/transverse/normal frame (`frame=:rtn`, re-expressed from the target's
state every frame), in a spacecraft's body frame (`:body`) or in inertial
axes, with a colour and dashed or solid style, toggled on the page.

## Reference ghosts

`references` draws a second, translucent copy of a spacecraft that follows a
state table the integrator did not produce: a SPICE reconstruction of the
real mission, a telemetry record, or a plan. The ghost copies the geometry
(and the 3D model) of the spacecraft it refers to, carries its own line over
the whole reference span and a ring marker, and the selection panel of the
flown spacecraft reports the separation between the two at the current time.
Times are seconds from the run epoch; positions are inertial (J2000,
planet-centred) metres, embedded as Float64.

```julia
# Sample the mission SPK on the run's saved times, relative to the planet centre.
et0 = str2et(scene_epoch_utc)
states = [spkezr("MAGELLAN", et0 + t, "J2000", "NONE", "VENUS")[1] for t in times_s]
pos_m = 1e3 .* reduce(hcat, [s[1:3] for s in states])
vel_mps = 1e3 .* reduce(hcat, [s[4:6] for s in states])
export_visualization(prefix; models=Dict(1 => "data/models/magellan_nasa_3d_resources.glb"),
    references=[(name="Magellan (SPICE)", t_s=times_s, pos_m=pos_m, vel_mps=vel_mps, target=1, color="#ff8c69")])
```

Attitude comes from a `q` column (4 x N, scalar-last) when given and is
velocity-aligned otherwise. `scripts/dev/viewer_demos/` has drivers that
build such ghosts for Magellan at Venus and Cassini at Titan.

## Robot arms

A spacecraft whose control model carries a `RobotArmControlEffector` with a
`RobotArmPlan` runs the coupled cloth-arm dynamics, and the viewer draws that
arm: the sidecar records the arm links (vector, radius, mass, mount offset)
from the plan's `ClothArmModel`, and the `arm_pose` save field records every
arm link's centre of mass relative to the spacecraft and its inertial
quaternion at each step, straight from the integrated `arm_r`/`arm_q` state.
Up close the arm appears as cylinders with joint spheres and a marked tip,
moving through the planned motion while the bus reacts. The bounding radius
used for the close-up switch grows to the arm's reach.

## Options

```julia
export_visualization(prefix;
    max_frames=2000,          # rows embedded (decimated evenly)
    data_budget_mb=150.0,     # hard cap on the embedded trajectory
    trail_orbits=3,           # or trail_s=...
    frame=:inertial,          # or :planet_fixed
    texture_resolution=:best, # or "4k" to keep the page small
    stl=Dict(1 => "bus.stl"), stl_scale=1.0,
    title="My run")
```

The page is self-contained (three.js, texture and data are embedded) and
opens from disk. An 8k Earth adds about 5 MB; pass `texture_resolution="4k"`
or `textures=false` to shrink it.

## Standalone viewer, no simulation needed

The viewer also runs on its own from a data file, with no Julia and no
SpaceAGORA run behind it. `viewer/build_standalone.py` assembles a single
HTML page that opens from disk:

```bash
python3 viewer/build_standalone.py --out output/spaceagora_viewer.html          # a form: pick the file, body and epoch in the page
python3 viewer/build_standalone.py --out page.html --data run.csv --planet venus --epoch 1993-05-26T00:00:07Z \
    --model data/models/magellan_nasa_3d_resources.glb --model-rotation 0,0,90    # opens straight into the run
python3 viewer/build_standalone.py --out page.html --data nominal.csv s1.csv s2.csv s3.csv --planet mars  # an ensemble
```

The data is a CSV with a `time` column (seconds) and either SpaceAGORA's own
result columns (`sc1_pos_1..3`, `sc1_vel_1..3`, `sc1_q_1..4`, drag, density,
link poses; the `simulation_results.csv` a run writes opens as is, several
spacecraft included) or plain `x,y,z[,vx,vy,vz][,qx,qy,qz,qw]`, in metres or
kilometres, positions inertial (J2000) about the body's centre; a JSON form
with `time` and `spacecraft: [{pos, vel, q}]` is accepted too. The page
builds everything the Julia bundler would: the body's rotation from the IAU
2009 pole and prime meridian at the given epoch (so ground tracks agree with
the SPICE-driven pages to a fraction of a degree), a box spacecraft of the
given size, the model override (centred from the parsed geometry), and an
optional reference ghost from a second table. Several files open as an
ensemble with the first as the nominal. `--cdn` builds the variant the
claude.ai artifact host can show (three.js from a CDN); the default page is
fully offline. The page keeps a "Load other data" button, so one built page
serves any number of files.

## Ensembles

Monte Carlo samples and constellation members merge into one page:

```julia
result, page = run_monte_carlo_visualization(seed -> build_args(seed), 1:50, "output/campaign";
    threads=4, scalar=df -> df[end, "sc1_periapsis_altitude"] / 1000, scalar_name="final periapsis (km)")
```

Each sample writes its bundle to `output/campaign/sample_0001` and so on,
the campaign writes `ensemble_manifest.json`, and `ensemble_viewer.html`
overlays every sample coloured by the scalar, with a sample selector and a
faint full-history view of all of them. `export_ensemble_visualization(dir)`
does the same for any directory of `sample_NNNN` or `sat_<i>_id_<id>`
subdirectories, including the per-member output of
`run_constellation_ensemble`. From the CLI:
`spaceagora visualize --run=output/campaign --ensemble`.

### Nominal sample and the 3σ tube

`run_monte_carlo_visualization(...; nominal=seed)` runs one more sample with
that seed (the unperturbed configuration) labelled "nominal", records it in
the manifest, and the page draws it as a bright white line with its own
label while the Monte Carlo traces stay faint and coloured by the sample
scalar. Around the nominal (or the sample mean when there is none) a
translucent tube shows the 3σ dispersion of the samples at every time, its
cross-section the radial and cross-track standard deviations in the
nominal's RTN frame, and the ensemble panel reports the 3σ radial,
along-track and cross-track values at the current time together with how
many samples are present. The tube and the traces can be switched off and
the σ multiple set to 1, 2 or 3 in the panel. The standalone page does the
same for several files.

## Atmosphere

Runs with an atmosphere model get three more layers, each with a toggle:

- **Limb glow** at the entry-interface altitude, so the pass visibly dips
  into the atmosphere and the EI reads as a boundary.
- **Density shells** between 0.3 and 1.0 of the EI, with opacity from a
  density profile the sidecar samples from the run's own model (the engine
  samples GRAM and NRLMSISE-00 too, using its integrator state).
- **Density map** at 0.6 of the EI for models that vary with latitude and
  longitude (GRAM, NRLMSISE-00, tabulated), draped on the body as an inferno
  colour scale with its range in the info panel. Altitude-only models skip it.

The pass itself is coloured by heat rate by default; the trail colour
selector also offers dynamic pressure (from the saved density and speed),
density, altitude and speed, with a legend. The selection panel reads out
density, dynamic pressure, heat rate, drag and wind speed at the scrubbed
time. These come from the `density` save field, added alongside `link_pose`
when the flag is on, and the existing heat-rate, drag and wind columns.

## Textures

Textures live in `data/textures/` with `manifest.toml` recording body,
resolution tier, source, licence and the longitude of the image's left
edge. `scripts/dev/build_viewer_textures.py --tier 8k --body earth` rebuilds a
tier from its public-domain source; the 16k Earth tier is possible but not
committed. `assets check` reports the directory. Adding a body is one JPEG
and one manifest entry.
