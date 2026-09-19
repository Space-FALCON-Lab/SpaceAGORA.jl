# Lunar Landing

A powered descent from orbit to the surface, flown six-degree-of-freedom: an
Apollo-style quadratic guidance commands the descent engine's thrust and the
attitude that points it, an RCS attitude controller tracks that attitude, and a
local terrain bundle serves as the radar altimeter, ends the run at touchdown
and shapes the surface the viewer draws.

`scripts/dev/viewer_demos/apollo11_landing.jl` flies Apollo 11 from powered
descent initiation (PDI, 1969-07-20 20:05:05 UTC, 15.24 km above Tranquility
Base and 480 km uprange on the descent orbit) to touchdown, with NASA's lunar
module display model and a terrain bundle of the site. It demonstrates the
descent, terrain-contact and plume models on one documented case. It is not a
reconstruction of the flown trajectory and it does not establish flight
accuracy: the targets are in the spirit of Apollo 11's, the vehicle is a
box-and-inertia stand-in with the LM's mass, thrust and quad layout, and the
terrain is whichever bundle you point it at.

The descent types live in `SpaceAGORA.SimulationModel`; the examples below use
`SM = SpaceAGORA.SimulationModel`. The terrain loaders and the plume model are
exported at the package root.

The requested start remains exactly `1969-07-20T20:05:05.000` UTC on the
engine's millisecond clock. The demo regression checks this against SPICE,
including the pre-1972 calendar conversion.

## Terrain models

`DEMTerrainModel` holds one or more `DEMGrid`s, regular planetocentric
latitude and longitude grids of heights above one explicit reference sphere,
finest first, and answers `terrain_height(model, lat_deg, lon_deg)` by
bilinear interpolation in the first grid that covers the point, falling back
to its configured fallback height elsewhere; `NoTerrainModel()` is the flat
sphere. `load_site_terrain("data/terrain/moon/apollo11/site.json")` reads a site
bundle: `site.json` names the site and lists the DEMs in priority order, each
with a JSON header, a little-endian Float32 `.f32` array of metre heights at
cell centres, and the reference radius every grid must share. The loader
returns the model and a site record with the coordinates, the name, the
terrain height at the site and that radius.

`scripts/dev/terrain/fetch_moon_site.py` writes such a bundle for a lunar
site: the LOLA global grid in a window around the site and a coarser LOLA
window over the whole descent corridor, the LROC NAC digital terrain model
where one exists (Apollo 11 does), and an imagery quadtree along the descent
corridor for the viewer. Its sources, options, registration corrections,
resolution measurements and Python requirements are documented in the
repository's `docs/reference/lunar_imagery_sources.md`. Downloading is a
separate step from running the simulation:

```
python3 scripts/dev/terrain/fetch_moon_site.py --site 0.67416 23.47314 --nac-half-deg 0.03 --out data/terrain/moon/apollo11
```

The NAC window is narrowed to ±0.03° for Apollo 11: the registered NAC digital
terrain model is a narrow north–south strip and the site sits near its eastern
edge, so the fully valid square around the site is about ±0.036°; the producer
refuses a window that contains no-data cells rather than invent terrain, and the
tool's default of ±0.06° therefore fails at this site.

The bundle stays under `data/` (not tracked); `SPACEAGORA_TERRAIN_SITE` points
the demo at a `site.json` somewhere else, such as a regenerated bundle under
`output/terrain`. Any DEM in the same layout can be dropped in, including a
synthetic one for tests: the unit test beside the demo writes a flat two-grid
bundle and runs the descent models against it without network access.

Use the bundle's reference radius as the explicit datum everywhere the descent
touches the ground: the guidance configuration, the touchdown event and the
plume model all take it, and they reject a radius that differs from the DEM's.

## Powered descent guidance

`SM.ApolloDescentGuidanceModel` implements the LM guidance computer's quadratic
law (Klumpp, *Apollo Lunar Descent Guidance*, 1974) in a site-fixed frame
(x uprange, y crossrange, z up; velocities relative to the rotating body):

```math
a_c = a_T - \frac{6\,(v_T + v)}{T} + \frac{12\,(r_T - r)}{T^2}
```

with the time-to-go `T` re-solved every cycle from the uprange cubic that adds
a jerk target. Two quadratic phases follow the Apollo programs, P63 braking
(targets beyond the high gate, the phase ends at the high-gate altitude) and
P64 approach (a hover point above the site), then a P66 rate-of-descent phase
nulls horizontal velocity and descends at a configured rate, slower below a
configured altitude, until the engine cutoff altitude.
`SM.apollo11_descent_targets()` returns targets in the spirit of Apollo 11's;
`SM.DescentPhaseTargets` lets you set your own.

The thrust command is `m |a_c - g|` with the rotating-frame terms included;
the descent engine's throttle follows the DPS envelope (`throttle_min` to
`throttle_max`, 10 to 60 percent by default, and full thrust above it). The
attitude command points body `-z` (the engine axis) along the thrust with the
windows (body `+x`) along the projection of up plus the flight direction:
windows up during braking, facing the site once the vehicle pitches back.
`SM.descent_attitude_command` returns a scalar-last quaternion in the engine's
initial-condition convention: `rot(q)` maps inertial vectors into the body
frame, and its transpose maps body axes into the inertial frame.

`SM.ApolloDescentControlModel` throttles the engine toward the command with a
slew limit, runs a rate-limited attitude loop whose torques are bounded by the
RCS authority, applies the thrust along the vehicle's actual engine axis, and
drains propellant through both. The two effectors share an
`SM.ApolloDescentState` (phase, time-to-go, thrust, commanded attitude, radar
altitude, site-frame state, touchdown record). Run the simulation with
`isolate_state=false` when you want to read that state afterwards; the demo
also saves thrust, throttle, time-to-go, radar altitude, phase and attitude
error as result columns.

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
[Stopping on a Condition](stop_conditions.md)): the engine replaces its impact
event with a touchdown event on the terrain, the run ends when the vehicle's
reference point reaches `touchdown_height_m` above the ground, and the
touchdown time, ground-relative velocity and miss distance are recorded in the
shared state. The `touchdown_height_m` is the height of the vehicle's
reference point above the footpads: 3.7 m for the display model at the demo's
scale.

The control effector also reports what its thrusters are doing, through
`control_thruster_levels`: the descent engine's firing level is its actual
thrust over the thruster's rating, and the RCS jets' levels come from a
least-norm allocation of the commanded body torque over their torque arms
about the spacecraft reference point, clipped into 0 to 1 and refreshed every
control cycle. The demo's lunar module carries one descent engine and sixteen
jets in four quads (two lateral, one up and one down per quad), so the run
writes `sc1_thruster_level_1..17` and the viewer draws a plume on each firing
thruster: the engine burns throughout the descent and the jets puff at the
phase handovers. The layout is what the allocation sees; the thrusters on the
link are display geometry and the control effector applies the forces.

## Plume-surface interaction

[`PlumeSurfaceInteractionModel`](@ref) is a dynamic effector for what the
descent engine's exhaust does to the ground: the pressure and shear stress it
lays down, the regolith it erodes, the speed the grains leave at, and the small
thrust augmentation the reflected plume gives the vehicle in the last couple of
nozzle diameters. [Plume Interaction](plume_interaction.md) describes the
closure, its parameters and its limits. It is constructed with the descent
control effector (whose actuator state carries the engine's actual thrust), the
terrain and the same explicit datum, so it always sees the throttle the
controller flew and the ground the altimeter measured:

```julia
plume = PlumeSurfaceInteractionModel(control, terrain; reference_radius_m=site.reference_radius_m)
# ... dynamic_effectors = (gravity..., plume)
```

When the effector is in the run, the default save fields publish its state
as seven result columns per spacecraft: `sc{i}_plume_height_m`,
`sc{i}_plume_pressure_pa`, `sc{i}_plume_shear_pa`, `sc{i}_plume_erosion_kg_s`,
`sc{i}_plume_eroded_kg`, `sc{i}_plume_ejecta_mps` and
`sc{i}_plume_ground_effect_n`. `plume_erosion_onset_height(config, thrust_n)`
gives the height at which erosion starts for a given thrust, for sizing a
scenario without running one; the demo prints it beside the first erosion it
saw, the peaks and the total eroded mass. Those values are what the
phenomenological closure produces with its default parameters on this case;
they are not Apollo flight measurements.

## Surface in the viewer

`export_visualization(prefix; terrain="data/terrain/moon/apollo11/site.json")`
embeds the bundle's grids and, when the bundle has one, its imagery quadtree.
The page drapes each node it needs over geometry displaced by the DEM, so the
surface sharpens as the camera closes in on the site and stays covered out to
the horizon while the lander is still uprange. The selection panel gains a
height-above-terrain quantity and the info panel names the terrain and its
finest resolution. [Interactive Visualization](visualization.md) describes
the export options and the bundle contract.

## Running the demo

```
julia --project=. scripts/dev/viewer_demos/apollo11_landing.jl --output-dir output/apollo11-landing
```

The demo needs the site bundle (or `SPACEAGORA_TERRAIN_SITE`), the tracked
lunar module display model in `data/models/`, the LP165P gravity coefficients
in `data/Gravity_harmonics_data/` and the SPICE starter assets
(`SPACEAGORA_SPICE_PATH`) for the lunar orientation; it needs no atmosphere and
no native GRAM. `--duration-s` is only a cap: the touchdown event ends the run.
An existing nonempty output directory is refused. The run prints the PDI state,
the touchdown time, vertical and horizontal contact speed, miss distance and
propellant used, the phase start times, the plume summary and the number of
thruster-level columns, then writes the standalone viewer page. Optional
`SPACEAGORA_DEMO_CDN=1` invokes the separate CDN page tool when it is installed.


The output directory contains `simulation_results.csv`,
`simulation_results.feather`, `simulation_results_scene.json` and the
self-contained `simulation_results_viewer.html`. Open the HTML file, click the
spacecraft marker or label, and press **F** to follow the lander. Pause and move
the time slider to inspect braking, approach and the final vertical descent.
The selection panel shows height above terrain; click a value to plot its
history. Switch between inertial and planet-fixed views to inspect the same
landing in either frame. Near the surface, the viewer shows the site imagery,
engine and RCS plume glyphs, and modeled dust.

Results are saved every 0.25 seconds. The touchdown event can fall between those
samples, so the last saved row and viewer frame can still show thrust shortly
before contact. The printed touchdown time and the returned control state refer
to the actual event, where engine thrust, torque and thruster levels are zeroed.
Early radar-altitude entries can be unavailable before the first guidance
update. These gaps are not negative terrain clearance.
