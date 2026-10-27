# Lunar Landing

SpaceAGORA can fly a powered descent from orbit to the surface in six degrees
of freedom. An Apollo-style quadratic guidance law commands the descent
engine's thrust and the attitude that points it, an RCS attitude controller
tracks that attitude, and a local terrain bundle acts as the radar altimeter,
ends the run at touchdown and gives the viewer its surface.

`scripts/dev/viewer_demos/apollo11_landing.jl` flies Apollo 11 from powered
descent initiation (PDI, 1969-07-20 20:05:05 UTC, 15.24 km above Tranquility
Base and 480 km uprange on the descent orbit) to touchdown, with NASA's lunar
module display model and a terrain bundle of the site. It demonstrates the
descent, terrain-contact and plume models on one documented case. It does not
reconstruct the flown trajectory or establish flight accuracy: the targets are
only modeled on Apollo 11's, the vehicle is a box with a stand-in inertia and
the LM's mass, thrust and RCS quad layout, and the terrain is whichever bundle
you point it at.

The descent types live in `SpaceAGORA.SimulationModel`, written `SM` below
(`SM = SpaceAGORA.SimulationModel`). The terrain loaders and the plume model
are exported from the package root.

The demo starts at exactly `1969-07-20T20:05:05.000` UTC on the engine's
millisecond clock. When the SPICE kernels are installed, its unit test checks
that epoch against SPICE, including the calendar conversion for a date before
1972.

## Terrain models

`DEMTerrainModel` holds one or more `DEMGrid`s in priority order, usually
finest first. Each grid is a regular grid of heights in meters above one
explicit reference sphere, indexed by planetocentric latitude and
east-positive longitude. `terrain_height(model, lat_deg, lon_deg)`
interpolates bilinearly in the first grid that covers the point and returns
the model's fallback height (zero by default) everywhere else.
`NoTerrainModel()` is the bare sphere.

`load_site_terrain("data/terrain/moon/apollo11/site.json")` reads a site
bundle. `site.json` names the site and lists its DEMs in priority order, each
with its reference radius, and all of them must share that radius. Every DEM
has a JSON header and a `.f32` file of little-endian Float32 heights in meters
at cell centers. The loader returns the terrain model and a site record with
the site's coordinates, its name, the terrain height there and the reference
radius.

`scripts/dev/terrain/fetch_moon_site.py` writes such a bundle for a lunar
site: a window of the LOLA global grid around the site, a coarser LOLA window
over the whole descent corridor, the LROC NAC digital terrain model where one
exists (Apollo 11 has one), and an imagery quadtree along the corridor for the
viewer. [Building lunar terrain and imagery](lunar_imagery_sources.md)
documents its sources, options, registration corrections, resolution
measurements and Python requirements. Downloading is a separate step from
running the simulation:

```sh
python3 scripts/dev/terrain/fetch_moon_site.py --site 0.67416 23.47314 --nac-half-deg 0.03 --out data/terrain/moon/apollo11
```

For Apollo 11 the NAC window is narrowed to ±0.03°. The registered NAC
terrain model there is a narrow north–south strip with the site near its
eastern edge, so the fully valid square around the site is only about
±0.036°. The tool refuses a window that contains no-data cells rather than
invent terrain, so its default of ±0.06° fails at this site.

The bundle lives under `data/terrain/`, which Git ignores.
`SPACEAGORA_TERRAIN_SITE` points the demo at a `site.json` somewhere else,
such as a regenerated bundle under `output/terrain`. Any DEM in the same
layout can be used, including a synthetic one: the unit test
`test/unit/gnc/apollo11_landing_tests.jl` writes a flat two-grid bundle and
builds the descent and plume models on it without network access.

Use the bundle's reference radius as the explicit datum wherever the descent
meets the ground. The guidance configuration, the touchdown event and the
plume model all take it, and they reject a radius that differs from the DEM's.

## Powered descent guidance

`SM.ApolloDescentGuidanceModel` implements the quadratic law of the LM
guidance computer (Klumpp, *Apollo Lunar Descent Guidance*, 1974) in a
site-fixed frame: x uprange, y crossrange and z up, with velocities relative
to the rotating body.

```math
a_c = a_T - \frac{6\,(v_T + v)}{T} + \frac{12\,(r_T - r)}{T^2}
```

The time-to-go `T` is solved again every cycle from the uprange cubic that
includes a jerk target. Two quadratic phases follow the Apollo programs: P63
braking, whose targets lie beyond the high gate so that the phase ends at the
high-gate altitude, and P64 approach, which aims at a hover point above the
site. A P66 rate-of-descent phase then nulls the horizontal velocity and
descends at a configured rate, slower below a configured altitude, until the
engine cutoff altitude; the defaults are 1 m/s, 0.5 m/s below 10 m, and
cutoff at 1.7 m. `SM.apollo11_descent_targets()` returns targets modeled on
Apollo 11's, and `SM.DescentPhaseTargets` lets you set your own.

The thrust command is `m |a_c - g|`, with the rotating-frame terms included.
The descent engine's throttle follows the DPS envelope: it throttles between
`throttle_min` and `throttle_max` (10 and 60 percent by default), and a demand
above that band commands full thrust. For the first `ignition_trim_s` (26 s by
default) after the first guidance update, the engine holds `trim_throttle`
(10% by default). The attitude command points body `-z`, the engine axis,
along the thrust and turns the windows (body `+x`) toward the projection of up
plus the flight direction, so the windows face up during braking and face the
site once the vehicle pitches back. `SM.descent_attitude_command` returns a
scalar-last quaternion in the engine's initial-condition convention: `rot(q)`
maps inertial vectors into the body frame, and its transpose maps body axes
into the inertial frame.

`SM.ApolloDescentControlModel` slews the engine thrust toward the command at
a limited rate, runs a rate-limited attitude loop whose torques are bounded by
the RCS authority, applies the thrust along the vehicle's actual engine axis,
and drains propellant through both the engine and the RCS. The two effectors
share an `SM.ApolloDescentState` (phase, time-to-go, thrust, commanded
attitude, radar altitude, site-frame state and touchdown record). Build both
from the same state with the same `spacecraft_indices`, because the control
flies only the commands its guidance writes there. When a configuration has
both models, `run_simulation` refuses it before the run starts if a spacecraft
is guided or controlled without its counterpart on that state. A guidance
model of another type may also write the state, so the check skips spacecraft
without Apollo guidance, and likewise for control. To
read the state after the run, run with `isolate_state=false`, as the demo
does, or pass `return_results=true` and read it from the configuration the
results return. The demo also saves thrust, throttle, time-to-go, radar
altitude, phase and attitude error as result columns.

The fragment below assumes `using SpaceAGORA` and the Apollo 11 bundle in
place; the commented lines show where the models go in the configuration:

```julia
SM = SpaceAGORA.SimulationModel
terrain, site = load_site_terrain("data/terrain/moon/apollo11/site.json")
braking, approach = SM.apollo11_descent_targets()
gcfg = SM.ApolloDescentConfig(reference_radius_m=site.reference_radius_m,
    site_lat_deg=site.lat_deg, site_lon_deg=site.lon_deg, site_height_m=site.height_m,
    approach_azimuth_deg=270.0, braking=braking, approach=approach)
state = SM.ApolloDescentState(1)
guidance = SM.ApolloDescentGuidanceModel(gcfg, state, terrain)
control = SM.ApolloDescentControlModel(SM.ApolloDescentControlConfig(touchdown_height_m=3.7), gcfg, state, terrain)
# ... guidance_model=SM.GuidanceModel(guidance_effectors=(guidance,), guidance_rates=[2.0]),
#     control_model=SM.ControlModel(control_effectors=(control,), control_rates=[0.05]), orientation_sim=true
```

The control effector reports a landing through `touchdown_spec` (see
[Terrain contact for a landing controller](stop_conditions.md#Terrain-contact-for-a-landing-controller)).
For that spacecraft the engine replaces the usual impact stop with a
touchdown event on the terrain: the run ends when the vehicle's reference
point comes down to `touchdown_height_m` above the ground, and the touchdown
time, ground-relative velocity and miss distance are recorded in the shared
state. `touchdown_height_m` is the height of the vehicle's reference point
above the footpads. The default is 2.5 m; the demo uses 3.7 m for the display
model at its scale.

The control effector also reports what its thrusters are doing through
`control_thruster_levels`. The descent engine's firing level is its actual
thrust over its rating. The RCS jets' levels come from a least-norm
allocation of the commanded body torque over their torque arms about the
spacecraft reference point, clipped to the range 0 to 1 and refreshed every
control cycle. The demo's lunar module carries one descent engine and sixteen
jets in four quads (two lateral, one up and one down per quad), so the run
writes `sc1_thruster_level_1` to `sc1_thruster_level_17` and the viewer draws
a plume on each firing thruster: the engine burns throughout the descent, and
the jets fire briefly at the phase changes. In the default `actuation=:ideal`,
the thruster layout only feeds that allocation; the thrusters on the link are
display geometry, and the control effector applies the forces itself.

With `actuation=:mounted`, the mounted thrusters apply the forces instead,
and the descent engine can carry gimbals (`ThrusterGimbal`). The gimbals hold
their angles unless `ApolloDescentControlConfig(engine_trim=true)`. Then each
control cycle commands them to the angles that put the engine's thrust line
through the center of mass (the body origin), within the gimbals' own angle
and rate limits, as the LM's gimbal drive trimmed the descent engine. Trim
matters when the engine is mounted off the axis through the center of mass.
`engine_trim=true` without `actuation=:mounted` is refused.

## Plume-surface interaction

[`PlumeSurfaceInteractionModel`](@ref) is a dynamic effector for what the
descent engine's exhaust does to the ground: the pressure and shear stress it
applies, the regolith it erodes, the speed at which the grains leave, and the
small thrust augmentation the reflected plume gives the vehicle within two
nozzle diameters of the ground (the default `ground_effect_cutoff`).
[Plume Interaction](plume_interaction.md) describes the closure, its
parameters and its limits. The model is built from the descent control
effector, whose actuator state carries the engine's actual thrust, together
with the terrain and the same explicit datum, so it always sees the throttle
the controller flew and the ground the altimeter measured. Continuing the
fragment above:

```julia
plume = PlumeSurfaceInteractionModel(control, terrain; reference_radius_m=site.reference_radius_m)
# ... dynamic_effectors = (gravity..., plume)
```

When the effector is in the run, the default save fields add seven result
columns per spacecraft: `sc{i}_plume_height_m`, `sc{i}_plume_pressure_pa`,
`sc{i}_plume_shear_pa`, `sc{i}_plume_erosion_kg_s`, `sc{i}_plume_eroded_kg`,
`sc{i}_plume_ejecta_mps` and `sc{i}_plume_ground_effect_n`.
`plume_erosion_onset_height(config, thrust_n)` gives the height at which
erosion starts for a given thrust, which helps size a scenario without
running it. The demo prints it next to the first erosion it saw, the peak
values and the total eroded mass. Those values are what the phenomenological
closure gives with its default parameters on this case; they are not Apollo
flight measurements.

## Surface in the viewer

With `SpaceAGORAVisualizer` loaded,
`export_visualization(prefix; terrain="data/terrain/moon/apollo11/site.json")`
embeds the bundle's grids and, when the bundle has one, its imagery quadtree.
The page drapes each imagery node it needs over geometry displaced by the
DEM, so the surface sharpens as the camera closes in on the site and stays
covered out to the horizon while the lander is still uprange. The selection
panel gains a height-above-terrain quantity, and the info panel names the
terrain and its finest resolution.
[Regional terrain in exported pages](visualization.md#Regional-terrain-in-exported-pages)
describes the export options and the bundle contract.

## Running the demo

The demo loads the optional `SpaceAGORAVisualizer` renderer and then
continues in the examples environment, so prepare both environments once.
From the repository root:

```sh
julia --startup-file=no scripts/setup_visualizer.jl /path/to/visualizer-env
julia --startup-file=no examples/setup.jl
julia --project=/path/to/visualizer-env scripts/dev/viewer_demos/apollo11_landing.jl --output-dir output/apollo11-landing
```

The demo also needs the site bundle (or `SPACEAGORA_TERRAIN_SITE`), the
tracked lunar module display model in `data/scene/models/`, the LP165P
gravity coefficients in `data/planets/gravity/`, and the SPICE kernels from
the `data/GRAMSuite.jl` submodule (or `SPACEAGORA_SPICE_PATH`) for the lunar
orientation. It needs no atmosphere and no native GRAM. `--duration-s` is only
a cap, because the touchdown event ends the run, and an existing nonempty
output directory is refused. The run prints the PDI state, the touchdown
time, the vertical and horizontal contact speeds, the miss distance and the
propellant used, the phase start times, the plume summary and the number of
thruster-level columns, and then writes the standalone viewer page. Setting
`SPACEAGORA_DEMO_CDN=1` also runs the separate CDN page tool when it is
installed.

The output directory contains `simulation_results.csv`,
`simulation_results.feather`, `simulation_results_scene.json` and the
self-contained `simulation_results_viewer.html`. Open the HTML file, click the
spacecraft marker or label, and press **F** to follow the lander. Pause and
move the time slider to inspect braking, approach and the final vertical
descent. The selection panel shows height above terrain; click a value to
plot its history. Switch between the inertial and planet-fixed views to see
the same landing in either frame. Near the surface the viewer shows the site
imagery, the engine and RCS plume glyphs, and modeled dust.

The selection panel recomputes terrain clearance at each saved frame. The
recorded radar-altitude column updates only on the guidance cycle, so it can
lag that value.

Results are saved every 0.25 seconds. The touchdown event can fall between
two samples, so the last saved row and viewer frame can still show thrust
shortly before contact. The printed touchdown time and the returned control
state refer to the event itself, where engine thrust, torque and thruster
levels are zeroed. The first radar-altitude entries can be missing before the
first guidance update; these gaps are not negative terrain clearance.
