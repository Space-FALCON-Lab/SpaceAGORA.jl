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
    model_scale=2.4,                            # meters per model unit (or a Dict per id)
    model_rotation_deg=Dict(1 => (-90, 0, -90))) # XYZ Euler angles in the body frame
```

```text
spaceagora visualize --run=output --model=1=data/models/iss_nasa_3d_resources_b.glb --model-scale=2.4
```

STL, OBJ and glTF/GLB are supported and embedded in the page; a `.gltf` must
carry its buffers inline, and glTF files must not require Draco, meshopt or
KTX2 (the exporter refuses them and names the conversion command). By
default the model's bounding-box center is placed on the spacecraft
(`model_center=false` keeps the file's own origin). The model sits in the
body frame at the root link, so attitude, link poses and glyphs still
apply, and the selection panel's "3D model" row reports the parsed mesh
count or the reason a model failed. `data/models/README.md` lists the
shipped models; the NASA ISS model there is public domain.

The same file can feed a planner. `sample_model_pointcloud(path; n_points,
scale, rotation_deg)` returns surface samples after the viewer's scale and
rotation, centerd the same way, so an RPO station built from it (for
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
axes, with a color and dashed or solid style, toggled on the page.

### Posing parts of a model

`model_articulations=Dict(id => [(region=..., axis=..., angle_deg=...), ...])`
rotates parts of a model before the scale and rotation: every vertex inside
`region` (an axis-aligned box in model units, `x_min`/`x_max`/... , missing
bounds unbounded) turns by `angle_deg` about `axis` through `pivot`
(default the region's bounding-box center). The same articulations go to
`mesh_aero_panels` for the aerodynamics, so the picture and the drag agree.
The Magellan demo turns the NASA model's cruise-canted wings broadside this
way.

## Site terrain

`export_visualization(prefix; terrain="data/terrain/moon/apollo11/site.json")`
embeds a landing site's digital elevation grids and imagery levels (see
[Lunar Landing](lunar_landing.md) for the fetch script). The page drapes each
imagery level over a patch displaced by the DEM, nested from a 4° window at
83 m/px down to a few hundred meters at 0.65 m/px and drawn all the time, so
the ground sharpens by itself as the camera closes in; the globe is cut open
under the outermost patch, a ring marks the site, and the selection panel
reports the height above the terrain.

What is draped is the albedo-normalized copy of each level when the site
directory carries one -- `fetch_moon_site.py` writes `level_k_albedo.jpg` beside
every `level_k.jpg`, and `terrain_payload` embeds both as `url` and
`albedo_url` -- because the page lights the ground itself and a mosaic carries
the illumination of whenever it was taken. The normalization is a flat field:
each level is divided by its own Gaussian low-pass (a radius of an eighth of the
patch) and rescaled to its own mean, which removes every brightness gradient
broader than the kernel. Shading inside a crater is at the scale of the crater
and survives, and slow genuine albedo variation is flattened along with the
illumination; the fetch script's `derive_albedo` documents both limits. A site
without the albedo images is drawn from `url` as before.

Deriving a second copy of a site (new albedo images, say) needs no network:
`--reuse DIR` copies the DEMs and the imagery originals of an existing site
directory, and the landing demo draws whichever copy `--site-json PATH` or
`SPACEAGORA_TERRAIN_SITE` names.

```
python3 scripts/dev/terrain/fetch_moon_site.py --site 0.67416 23.47314 \
    --out output/terrain/moon/apollo11 --reuse data/terrain/moon/apollo11
SPACEAGORA_TERRAIN_SITE=$PWD/output/terrain/moon/apollo11/site.json \
    julia --project=. scripts/dev/viewer_demos/apollo11_landing.jl
```

## Lunar reflectance and terrain shadows

Regolith does not scatter like a Lambert surface: the full Moon is almost
uniformly bright across its disk, and the surface brightens sharply as the phase
angle closes. When the run carries `frames.sun_dir` the page shades the terrain
patches and the globe with a lunar reflectance instead -- the lunar-Lambert law,
Lommel-Seeliger scattering blended with Lambert by the phase-dependent
coefficient of McEwen (1991), times a Henyey-Greenstein particle phase function
and Hapke's shadow-hiding opposition surge `B0 / (1 + tan(g/2)/h)`. The phase
term is normalized to 1 at a phase angle of 45°, so an ordinary view keeps the
brightness the Lambert page had while a view down-sun is about 2.7 times
brighter and one at 150° about 4 times darker, which is the phase curve the Moon
has. The globe keeps its Lambert shading for bodies that are not airless and
dark (only the Moon, Luna and Mercury take the regolith law), and a run without
a sun direction never leaves the Lambert branch at all.

The relief is shaded from the DEM, not from the drawn mesh: on load the page
turns each height grid into a normal map at half the grid's cell, so slopes the
mesh cannot resolve still shade correctly under a moving sun. The finest grid
travels at its native resolution for this reason -- `terrain_payload` keeps it
to `finest_max_grid` (2048) samples a side against `max_grid` (512) for the
rest, because a NAC digital terrain model at its own 4 m carries slope that the
same grid at 8 m does not: 10.9 degrees RMS against 8.0 at Tranquility Base.
Heights ride as Int16 steps about a base (`heights_i16`, `height_base_m`,
`height_scale_m`), which resolves the grid to a centimeter in half the bytes
Float32 would take; `quantize_heights=false` restores the Float32 `heights`
block.

Below the DEM's own cell the page adds micro-relief: a tiled detail normal map,
built once at load from multi-octave value noise and a few hundred small
craters, scaled to 0.12 RMS slope (about 7 degrees at a meter, which continues
the site grid's own slope-versus-baseline trend) and sampled by the levels at
least 2 m/px sharp, one tile every 8 m of ground. It is a texture rather than a
claim about the site, and under a 10 degree sun it is most of what the eye reads
as ground; `options.microRelief = false` turns it off. Its own mipmaps fade it
out as the camera pulls away.

The ground also shadows itself. The height grids are uploaded as textures and
the terrain's fragment shader marches 64 samples toward the Sun, from 8 m to
6 km in a geometric progression, testing the height of the ray above the
reference sphere against the ground under it; the first solid hit stops the
march, and the edge is softened over the Sun's own angular radius, which makes a
penumbra that widens with the distance to the occluder. That horizon test
multiplies the vehicle's shadow map, so a crater rim throws a shadow hundreds of
meters long while the lander's legs still throw theirs inside the 100 m box the
shadow map covers. It costs one texture tap a sample and nothing at all where
the Sun is below the local horizon.

## Dust

When the run carried a `PlumeSurfaceInteractionModel` (see
[Lunar Landing](lunar_landing.md)) the results table gains the seven
`sc{i}_plume_*` columns, the bundler adds a `frames.plume` block, and the page
draws the regolith the descent engine blows off the surface. Nothing is drawn
where the erosion rate is zero, which is everywhere above the erosion onset
height, so the sheet appears by itself in the last tens of meters of a descent.

What it is drawn to look like is NASA Langley's plume-surface interaction
tests, where an engine fires into a bin of lunar simulant in a vacuum chamber:
seen from above, a continuous translucent veil streams out of the impingement
point in billowing filaments that fade with radius; seen from the side it stays
flat against the ground, a couple of meters thick, while a fine haze builds up
over the following seconds until the hardware behind it is barely readable. So
the page draws five things, all of them driven by the effector's numbers:

- **the sheet**, three very flat disks stacked over the impingement point whose
  fragment shader is fractional Brownian motion in polar coordinates, scrolling
  outward at the ejecta speed with about thirty angular cells, which is what
  makes the filaments radial. Its radius grows from a couple of plume
  footprints at the moment erosion starts to a few tens of meters, and its
  optical depth follows the erosion rate. Each disk is treated as a slab: a
  thin layer seen edge-on is a long way through, so the sheet reads as a low
  wall from a grazing camera and as a veil from overhead, and the disks have no
  visible rim at a shallow angle;
- **the ejecta**, a sparse few hundred points on the fast tail of the speed
  distribution, leaving at one to three degrees above the local slope and
  falling back under the body's gravity;
- **the scour crater**, a floor swept smooth and slightly brighter inside a
  ring of darker piled-up soil, growing with the mass already moved;
- **the haze**, a bank of sixteen view-aligned puffs -- flattened ellipsoids of
  fines, wide and low, each turned to face the camera about the local vertical,
  so the bank veils the scene from a camera standing in it as well as from one
  above it. Its opacity is an exponentially decaying integral of the erosion
  rate, so it builds up over seconds and thins out after the engine stops
  (`hazeOpacity` and `hazeDecayS` are the knobs; their defaults are the vacuum
  case, thinner and shorter-lived than a chamber's, which has residual gas to
  keep the finest grains aloft);
- **the sheet's shadow**, a soft ellipse on the ground offset and stretched
  away from the sun by the height of the dust over it -- at the 10.6 degrees
  Apollo 11 landed under it runs about fifty meters down-sun -- so the lit dust
  reads as a volume and not as a decal.

Everything the sun touches is scaled by a Henyey-Greenstein phase function with
g = 0.65: fine dust scatters strongly forward, which is why the sheet is
blinding looking down-sun through it and almost clear looking the other way,
the way the Apollo films and the Langley footage show it. The sun direction,
its irradiance and the camera's exposure come from the lighting handle below;
the shaders tone map and encode their output exactly as every lit material
does, so one exposure applies to the whole frame. A page with no lighting
handle -- or one that never raised its exposure -- keeps the fixed brightness
the sheet had before that module existed. Positions inside the sheet, the puffs
and the points are derived in the shaders from the clock, so a frame update
writes uniforms and nothing else.

`THREE.Points` does render into three r160's shadow map with a custom depth
material, and the module carries one, but it is off by default
(`pointShadows`): the particles are sized for the perspective camera, so in the
shadow camera's orthographic pass each lands on a texel or two and changes
almost nothing, and sized up a few hundred opaque disks would speckle the
ground where a dust cloud casts one soft darkening. The decal does that.

The "dust" toggle in the toolbar turns all of it off, and the selection panel
gains the seven plume quantities — engine height (m), plume shear and pressure
(Pa), erosion rate (kg/s), eroded mass (kg), ejecta speed (m/s) and ground
effect (N) — each with a time history like every other row.

## Sun lighting and path tracing

When the run's ephemerides can resolve the Sun -- SPICE with the body's
kernels furnished, or `SimpleEphemeridesModel` at Earth -- `default_save_fields`
adds `sun_dir_1..3`, the unit vector from the planet's center to the Sun in the
inertial frame of the saved positions, and the page lights the scene from it.
The directional sun follows that vector over the timeline, rotated into the
scene the same way positions are, so the terminator stands where it stood and
a landing at dawn is lit like one. The shadow map is fitted to a few tens of
meters around the followed vehicle each frame, so the vehicle shadows itself
and drops its shadow on the ground: Apollo 11 landed with the Sun 10.6 degrees
up, and the LM's shadow stretches five times its own height across Tranquility
Base. A run without the columns keeps the fixed light the viewer always had.

### Earthshine

Away from Earth, and when the ephemerides resolve it, `default_save_fields`
also writes `earth_dir_1..3`, the unit vector from the planet's center to
Earth, and the page adds a second directional light along it. Its irradiance is
the sunlight the Earth reflects back: 0.15 W/m² at the Moon with a full Earth,
scaled by the Lambert-sphere phase law for the Sun-Earth-body angle and tinted
toward the blue of ocean and cloud. Apollo 11 landed under a 70 degree Earth
phase, half lit, so the Earth stood 59 degrees up in the west at 0.076 W/m² --
one ten-thousandth of sunlight, thirteen stops down, which is exactly why the
earthshine in the Apollo photographs needed a long exposure and why it takes
one here too. The hemisphere guess it replaces survives only as the fallback
for a run without the columns; a run that has them gets, instead, the one fill
term a vacuum scene really has, the sunlight the ground bounces back up
(irradiance `albedo x E x cos(theta)` on a downward-facing surface, with the
sky side of the hemisphere light set to black).

### The physical camera

Exposure is a camera setting in the photographic sense rather than a fudge
factor. The sun's irradiance comes from the body's mean distance from the Sun
(1361 W/m² at 1 au, so 1361 at the Moon and 587 at Mars), the render buffer is
kept in units of "a white Lambertian surface facing the Sun reads 1.0", and one
scale factor turns a buffer value into an absolute radiance. The display
mapping is then an exposure value at ISO 100: `EV = log2(L x S / K)` with
`K = 12.5` and the radiance converted to luminance at 98 lm/W, the luminous
efficacy of unattenuated sunlight. A 0.12-albedo surface facing the Sun meters
EV +15.3 at 1 au -- the "sunny 16" exposure a photographer would have set,
arrived at from the radiometry rather than assumed -- and EV +14.1 at Mars.
`export_visualization(...; ev=11.0)` fixes the page's EV instead of letting it
adapt.

One EV cannot serve the whole page, though: four stops separate an overview of
the sunlit disc from a close-up of a landing site under a 10.6 degree sun, and
a setting that holds one clips the other. So the camera meters and adapts, the
way an eye does. Every few frames the scene is drawn into a 64-pixel render
target, which three leaves linear and untone-mapped, and the log-average
luminance of the pixels carrying any light at all -- the black sky is excluded,
or an overview of a small disc on a large dark frame would meter itself white
-- gives the EV that would put that average on middle gray. The camera walks
toward it with a half-second time constant, between EV +6 (an earthshine-lit
shadow side) and EV +17 (full sun on a bright surface); a frame with almost
nothing lit in it leaves the exposure alone. The cost is one small extra pass
with the interface meshes hidden and the shadow map frozen.

`[` and `]` are then exposure compensation, a third of a stop at a time, and
the adaptation carries them: the info panel reads `EV +12.4 (auto, comp +0.7)`.
The **auto exposure** checkbox pins the camera at whatever it is reading, and
`export_visualization(...; ev=...)` starts a page fixed at a given value. A
fixed frame of Tranquility Base at the metered EV +15.3 is genuinely dark --
the Sun was 10.6 degrees up, so the ground returns a sixth of what that EV is
set for, and the earthshine on the shadow side is thirteen stops below it --
which is what the adaptation and the compensation are for.

The **lighting** selector offers "path traced when paused". In that mode the
page draws in real time while the timeline runs or the camera moves; once both
have been still for about 300 ms it hands the scene to
[three-gpu-pathtracer](https://github.com/gkjohnson/three-gpu-pathtracer) and
accumulates progressively refined samples into the same canvas, counting them
in the info panel, until the next change resets it. The Sun becomes a circular
light of its true angular diameter (0.5334 degrees), so shadow edges carry the
penumbra they have in reality, against a black sky with no ambient term --
which on the Moon is the whole of it. Markers, labels, trails, thruster cones
and facet glyphs are left out of the traced scene; so is the globe sphere when
site terrain has cut a hole in it, since the hole is a shader `discard` the
path tracer does not run. The library is loaded from the CDN on demand: the
page built by `export_visualization` is self-contained and offers real-time
lighting only, while the CDN pages
(`scripts/dev/viewer_demos/build_cdn_page.py`, `build_standalone.py --cdn`)
carry the import-map entries for it. It is handed the same irradiances (the sun
disc's radiance is calibrated to the directional sun's) and the EV the meter has
reached at hand-off, so switching modes changes the sampling and not the
brightness. The meter stands down while an image is accumulating, and changing
the exposure does not discard one: it is applied on the way out of the renderer,
so an image already gathered is still the right image.

### Thermal foils

The 3D models NASA publishes carry their thermal foils as plain diffuse colors,
which renders a spacecraft that looks like painted cardboard. On a glTF the page
recognizes them -- by material name when the model gives a useful one, otherwise
by base color: a saturated warm color is gold-coated Kapton, a bright neutral one
aluminized Mylar or bare aluminum -- and rebuilds them as `MeshPhysicalMaterial`
with the metal's measured normal-incidence reflectance as its base color (gold
1.00/0.77/0.34, aluminum 0.91/0.92/0.92), metalness near one, a roughness of
about 0.4 for crinkled foil, and a clearcoat over the Kapton for the film above
the metal. Textured, transparent and mid-to-dark gray materials -- the black
blankets, the painted structure, the decals -- are left exactly as they were, and
the model status line in the info panel says how many materials were swapped.

A metal has no diffuse color of its own, so it is only as interesting as what it
reflects. The page renders a 64-pixel cube map of the lit scene around the
vehicle with the vehicle itself stepped aside, prefilters it into a PMREM, and
gives it to the foil materials as their environment map; it is rebuilt when the
Sun has turned by more than three degrees or the vehicle has moved half a
kilometer, so a descent re-probes a few tens of times and a cruise once. The
result is a lander whose skirt carries the ground's color underneath and the
black of space above it.

## Heating on the spacecraft

When the run saved density (`save_visualization_scene` records it) the
close-up shades every face of the links and of the 3D model by the
free-molecular incident energy flux, ½ρV³ cos θ with θ the angle between the
face's outward normal and the airspeed (inertial velocity minus the
co-rotating atmosphere), on an inferno log scale spanning three decades
below the run's peak, with the legend in W/cm². Faces turned away from the
flow stay cold, so a Magellan pass shows the leading dish and the wings
lighting up while the lee side stays dark. It is a geometric incidence
map with full accommodation, not the thermal model's heat rate (which is a
spacecraft-level number in the selection panel); the "heating" toggle turns
it off and restores the plain materials.

## Thruster plumes

A control effector that drives named thrusters reports each one's firing level
(0 to 1) through the `control_thruster_levels(effector, i)` hook. When any
effector reports levels for a spacecraft, the run saves them as
`sc{i}_thruster_level_{k}` columns, one per thruster in the scene's order (the
spacecraft's links in order, each link's `thrusters` in order), and the page
draws a plume on every firing thruster.

`ApolloDescentControlModel` computes them inside `calcControlEffect!`: the
descent engine's level is its actual thrust over the thruster's rating, and the
attitude jets' come from a least-norm allocation of the commanded body torque
over the jets' torque arms about the spacecraft reference point, each clipped
into 0 to 1. A maneuver burn (`BaseThrusterModel`) reports its commanded
throttle, and the RPO controller reports its six-axis allocation. An effector
that drives no thruster returns `nothing` and the columns are not written.

Each plume is an additive, depth-write-free cone pair (a bright core inside a
translucent shroud) attached to the thruster's link, running along the
thruster's `direction` vector, with a shader that fades it toward the tip and
flickers it in time. The rated thrust sets the look: a 45 kN descent engine
burns long and orange-white, a 445 N attitude jet puffs short and blue-white.
Reach and brightness follow the level through a fractional exponent rather
than linearly, because a jet holding an attitude asks for well under a percent
of its rating and a linear mapping would draw nothing at all; a thruster at
level 0 is not drawn, and its static cone glyph is dimmed while it is idle.

How much of that look survives depends on the air around the nozzle. An exhaust
is luminous because it has ambient gas to shock, mix with and burn against: a
descent engine at sea level draws a bright sooty column, by about 70 km
(1e-5 kg/m³) it is a thin diffuse glow, and in vacuum there is nothing left to
excite -- the Apollo films of the powered descent record a bare engine bell
against the ground, not a flame. So the page fades the shroud out and leaves a
faint blue-white core as the ambient density falls, following the logarithm of
the density between those two. The density comes from the run (`frames.density`,
saved whenever the visualization scene is); a run that saved none falls back to
the scene, where a missing atmosphere block means the run flew with
`NoAtmosphereModel` and any other model keeps the full atmospheric look.
`createPlumes` also takes `vacuum: true | false` to force either. The plume's
brightness is deliberately left out of the camera's exposure: it is emissive,
and the ten-fold lift a grazing lunar page applies to its ground would lift the
vacuum core back into a flame.
The "plumes" toggle turns them off, and the selection panel gains a
`thruster k level` row per thruster whose history plots like any other
quantity.

## Time histories and face inspection

Every number in the selection panel is a quantity with a history: altitude,
latitude, longitude, radius, speed, airspeed (the inertial velocity minus the
co-rotating atmosphere), mass, density, dynamic pressure, heat rate, drag,
wind and the separation from each reference ghost. Click one and a plot panel
opens with that quantity over the whole run, the playback cursor on it and
the current value in the header; click or drag in the plot to seek, tick
"log" for a logarithmic axis (chosen automatically for densities and heating
that span decades), and click the quantity again to close. The Monte Carlo
panel's 3σ readouts open their histories the same way.

The spacecraft body is interactive in the close-up: click a face of a link
box or of the 3D model and a face panel opens under the selection panel, with
a marker and outward-normal arrow on the picked face. It reports the face
(which link, box or model surface, the hit point in the link frame), the
incidence angle between the face's outward normal and the airspeed, the local
heating ½ρV³cosθ the overlay shades and the Newtonian ram pressure ρV²cos²θ,
the airspeed and its direction in body axes, the face normal in body axes
(which follows the link's recorded pose), the inertial position and velocity,
the attitude quaternion (the saved one, or the velocity-aligned attitude the
assembly is drawn with), and the picked link's pose when the run recorded
link poses. Each of those is clickable for its history too; vector
quantities plot one line per component. Esc clears the face first, then the
selection.

## Saving a video

The toolbar's "Save video…" renders the animation frame by frame between two
elapsed times at a chosen simulated-seconds-per-video-second rate, frames per
second and size (canvas, 720p, 1080p, 1440p), encodes it with the browser's
H.264 encoder and saves an MP4 (mp4-muxer, vendored). The render is driven
deterministically, so the file does not depend on the machine's frame rate;
a 1080p minute takes a few seconds on a laptop. Browsers without WebCodecs
H.264 (some Firefox builds) get a "Record WebM" fallback that captures ten
seconds of live playback. Note that the claude.ai artifact host blocks
downloads started by the page, so record from the offline page or the
standalone page.

## Reference ghosts

`references` draws a second, translucent copy of a spacecraft that follows a
state table the integrator did not produce: a SPICE reconstruction of the
real mission, a telemetry record, or a plan. The ghost copies the geometry
(and the 3D model) of the spacecraft it refers to, carries its own line over
the whole reference span and a ring marker, and the selection panel of the
flown spacecraft reports the separation between the two at the current time.
Times are seconds from the run epoch; positions are inertial (J2000,
planet-centerd) meters, embedded as Float64.

```julia
# Sample the mission SPK on the run's saved times, relative to the planet center.
et0 = str2et(scene_epoch_utc)
states = [spkezr("MAGELLAN", et0 + t, "J2000", "NONE", "VENUS")[1] for t in times_s]
pos_m = 1e3 .* reduce(hcat, [s[1:3] for s in states])
vel_mps = 1e3 .* reduce(hcat, [s[4:6] for s in states])
export_visualization(prefix; models=Dict(1 => "data/models/magellan_nasa_3d_resources.glb"),
    references=[(name="Magellan (SPICE)", t_s=times_s, pos_m=pos_m, vel_mps=vel_mps, target=1, color="#ff8c69")])
```

Attitude comes from a `q` column (4 x N, scalar-last) when given and is
velocity-aligned otherwise. `scripts/dev/viewer_demos/` has drivers that
build such ghosts for Magellan at Venus, Odyssey at Mars and Cassini at Titan.

## Robot arms

A spacecraft whose control model carries a `RobotArmControlEffector` with a
`RobotArmPlan` runs the coupled cloth-arm dynamics, and the viewer draws that
arm: the sidecar records the arm links (vector, radius, mass, mount offset)
from the plan's `ClothArmModel`, and the `arm_pose` save field records every
arm link's center of mass relative to the spacecraft and its inertial
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
spacecraft included) or plain `x,y,z[,vx,vy,vz][,qx,qy,qz,qw]`, in meters or
kilometers, positions inertial (J2000) about the body's center; a JSON form
with `time` and `spacecraft: [{pos, vel, q}]` is accepted too. The page
builds everything the Julia bundler would: the body's rotation from the IAU
2009 pole and prime meridian at the given epoch (so ground tracks agree with
the SPICE-driven pages to a fraction of a degree), a box spacecraft of the
given size, the model override (centerd from the parsed geometry), and an
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
overlays every sample colored by the scalar, with a sample selector and a
faint full-history view of all of them. `export_ensemble_visualization(dir)`
does the same for any directory of `sample_NNNN` or `sat_<i>_id_<id>`
subdirectories, including the per-member output of
`run_constellation_ensemble`. From the CLI:
`spaceagora visualize --run=output/campaign --ensemble`.

### Nominal sample and the 3σ tube

`run_monte_carlo_visualization(...; nominal=seed)` runs one more sample with
that seed (the unperturbed configuration) labeled "nominal", records it in
the manifest, and the page draws it as a bright white line with its own
label while the Monte Carlo traces stay faint and colored by the sample
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
  color scale with its range in the info panel. Altitude-only models skip it.

The pass itself is colored by heat rate by default; the trail color
selector also offers dynamic pressure (from the saved density and speed),
density, altitude and speed, with a legend. The selection panel reads out
density, dynamic pressure, heat rate, drag and wind speed at the scrubbed
time. These come from the `density` save field, added alongside `link_pose`
when the flag is on, and the existing heat-rate, drag and wind columns.

## Textures

Textures live in `data/textures/` with `manifest.toml` recording body,
resolution tier, source, license and the longitude of the image's left
edge. `scripts/dev/build_viewer_textures.py --tier 8k --body earth` rebuilds a
tier from its public-domain source; the 16k Earth tier is possible but not
committed. `assets check` reports the directory. Adding a body is one JPEG
and one manifest entry.
