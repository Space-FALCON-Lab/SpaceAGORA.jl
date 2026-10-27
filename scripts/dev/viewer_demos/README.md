# Local viewer demonstrations

Each script here sets up a mission, simulates it, and exports the saved results
as a self-contained HTML page. SpaceAGORA owns the mission configurations and
the export helpers in this folder; the separate SpaceAGORAVisualizer package
draws the page. Research reconstructions belong with their studies, not here.

## Setup

The HTML renderer comes from the optional `SpaceAGORAVisualizer` companion.
Prepare a separate environment for it once, from the repository root:

```bash
julia --startup-file=no scripts/setup_visualizer.jl /path/to/visualizer-env
```

Setup installs the commit pinned in `packages/VisualizerSource.toml` and does
not add the renderer to SpaceAGORA's own environment. Every demo loads the
companion itself, so start each one with `--project=/path/to/visualizer-env`.
For demos launched through `scripts/dev/run.jl`, set
`SPACEAGORA_VISUALIZER_ENV=/path/to/visualizer-env` so the launcher can find
that environment. The browser and page-building tools are under
`SpaceAGORAVisualizer.viewer_directory()` in the companion's installation.

The ISS rendezvous and the mission demonstrations below also continue in the
examples environment through `examples/common.jl`, so for those prepare it too
with `julia --startup-file=no examples/setup.jl`. The four basic demonstrations
do not need it.

## Basic demonstrations

These four use analytic ephemerides and either no atmosphere or an exponential
model. They need no native GRAM, SPICE kernels or private telemetry. The ISS
demonstration draws the tracked NASA model from `data/scene/models/`; that is
display geometry, not the mass or aerodynamic mesh.

| Script | Default simulated time | What it shows |
| --- | --- | --- |
| `iss_demo.jl` | 11,120 s | A J2 orbit with the NASA ISS display model |
| `arm_demo.jl` | 60 s | Planned, controlled robot-arm motion on an orbiting bus |
| `earth_4day.jl` | Four days | A long, eccentric Earth trajectory |
| `odyssey_two_orbits.jl` | Two initial Keplerian periods | Mars energy-depletion control with an exponential atmosphere |

For a short first run:

```sh
julia --project=/path/to/visualizer-env scripts/dev/viewer_demos/iss_demo.jl --duration-s 60 --output-dir output/iss-first-run
```

All four accept `--duration-s` and `--output-dir`, with the value after a
space as above; they do not accept the `--key=value` form of the main CLI.
Without `--output-dir` they write to `output/viewer_demos/<script-name>`, and
`SPACEAGORA_VIEWER_DEMO_OUT` changes that parent folder. A run refuses an
output folder that already has files in it, so it can never overwrite earlier
work; choose a new folder for each run.

Each run writes the simulation results, the scene description and the HTML page
into its folder, then prints the page's path. Open it once the simulation has
finished; the viewer replays the saved results. A short run checks the setup
and the export, but shows no full orbit, aerobraking pass or completed arm
maneuver. The Mars case shows the analytic model and the controller at work;
it is not fitted to flight data, and the minimum altitude it prints is the
lowest among the saved samples.

The Odyssey demonstration has a 0.1-second controller, so it raises the
solver's step limit to cover that many control stops over the requested
duration, with a margin for rejected steps. Only the step limit changes, not
the tolerances.

## Rendezvous around the ISS

`iss_hypr.jl` plans and flies a CubeSat approach around a point cloud sampled
from the tracked ISS display model. It needs the starter SPICE kernels in the
GRAMSuite asset folder, but not native GRAM or private telemetry. The station
model and the controller are a demonstration, not an ISS flight model.

Copy the prepared visualizer environment into a fresh HYPR-enabled environment,
leaving the original unchanged, and prepare the examples environment if needed:

```text
julia --startup-file=no scripts/setup_hypr.jl --migrate-from /path/to/visualizer-env /path/to/hypr-visualizer-env
julia --startup-file=no examples/setup.jl
```

Start with the short hop:

```sh
SPACEAGORA_DEMO_SMOKE=1 julia --project=/path/to/hypr-visualizer-env scripts/dev/viewer_demos/iss_hypr.jl
```

Leave out `SPACEAGORA_DEMO_SMOKE` for the complete approach:

```sh
julia --project=/path/to/hypr-visualizer-env scripts/dev/viewer_demos/iss_hypr.jl
```

The script prints the HTML path under
`output/viewer_demos/iss_hypr_<inputs digest>/`; `SPACEAGORA_VIEWER_DEMO_OUT`
changes the parent folder. The saved plan and its provenance sit beside the
results. Running again with the same inputs reuses the recorded simulation and
plan only if the plan, scene, Feather results and control log all match the
hashes in the provenance. A missing, empty or changed file starts a new run,
and so does an older output without a recorded Feather hash. A reused run
shows the plan that produced its trajectory. The short hop only checks the
pipeline; it does not show the complete approach.

Open the printed HTML, click the marker beside the station label and press
**F** to follow the chaser. The blue dashed line is the reference the chaser
followed, in the station's local orbital frame, and the yellow line joins the
planner's waypoints. Playback shows saved simulation states.

A seeded plan with a fixed iteration budget comes out the same at any thread
count; a wall-clock stopping budget can stop at a different iteration. Older
versions drew different random streams for the particles, so check tracking
and clearance against the plan recorded with each run.

The reference speed is capped at 0.25 m/s, and the retimer also limits the
acceleration along the path, from the chaser's starting speed down to rest at
the goal. The run ends a fixed hold time after the reference does (60 s, or
10 s for the short hop) and never lasts more than 1800 s; a plan whose
reference would end later is refused. The LQ-MPC controller has no collision
constraints, so when you change the scenario check the simulated path and its
clearance, not only the planned path.

The station rotates with its circular orbit so that its body-frame point cloud
stays aligned with the planner's RTN frame and the saved 3D attitude. This
relies on the example's circular equatorial orbit and principal-axis rotation;
it does not model rendezvous with a tumbling station.

## Mission demonstrations and SPICE comparisons

`apollo11_lunar_orbit.jl`, `magellan_aerobraking.jl`, `odyssey_aerobraking.jl`
and `cassini_titan_flyby.jl [TA|T5|all]` are the mission cases from the
original viewer work. The first three default to two modeled orbits. Each
Cassini flyby defaults to five hours starting 2.5 hours before the closest
approach found by search, and `TA` is the default flyby. All four accept
`--duration-s <positive seconds>` and `--output-dir <new folder>` for a bounded
run, except that `all` refuses `--output-dir` because each flyby needs its own
default folder. An existing output folder with files in it is refused, never
reused or cleared.

These drivers use the tracked mission display models in `data/scene/models/`
and need existing SPICE kernels (`SPACEAGORA_SPICE_PATH`). Magellan, Odyssey and
Cassini also need native GRAM installed; each builds its GRAM model explicitly
at the simulation's initial epoch. That does not validate GRAM's native time
system, coordinates, datum or climatology. Mission navigation SPKs are fetched
from NAIF only when missing, into `SPACEAGORA_MISSION_SPK_DIR` or the
`spk/missions` folder of the configured SPICE tree. Native GRAM itself is never
downloaded.

Odyssey asks for TES mapping year 2 with `mars_map_year=2`. In the native
wrapper this was tested with, that request alone does not apply the year to the
native model; a later parameter-setting call is needed. The example keeps its
existing numerical setup, so do not read it as a validated year-2 atmosphere.
An explicitly applied profile, and its effect on the trajectory, would need
their own validation.

The references are central-body-relative, geometric J2000 SPICE states in SI
units at the saved simulation times. Times without SPK coverage are reported
and left out. The target is given by its one-based spacecraft index, not its
public id. The separation is a comparison between sampled models, not a measure
of mission reconstruction accuracy. Apollo uses an illustrative parking orbit
with SPICE orientation, not a flown LM reference. Without integrated attitude
these examples use the viewer's velocity-aligned pose, so a NASA model on
screen is not a reconstructed attitude history.

Titan's static 5x5 field is zero-tide, with fully normalized coefficients, a
reference radius of 2,575,000 m and GM = 8,978,126,919,238.97 m^3/s^2 from NASA
PGDA product 91; no time-varying tide is modeled. Setting
`SPACEAGORA_DEMO_ODYSSEY_ATMOSPHERE=accelerometer` swaps GRAM for the tracked
accelerometer density table
`data/telemetry/Odyssey/odyssey_accelerometer_density.feather`; the demo never
makes that swap on its own. `SPACEAGORA_DEMO_MAGELLAN_ANTENNA=forward|aft` sets
the model's pose. `SPACEAGORA_DEMO_MESH_AERO=1` turns on SpaceAGORA's mesh
aerodynamic model and fits the selected geometry from scratch rather than
reusing an earlier fit.

The output is the standalone `simulation_results_viewer.html`. With
`SPACEAGORA_DEMO_CDN=1` a run also calls the separate `build_cdn_page.py` tool
from the companion's viewer folder, and stops with an error if that tool is not
installed (the standalone page is already written by then).
`check_odyssey.jl [results_directory]` reads the saved articulation and altitude
of the basic Odyssey demonstration without running it again; with no argument
it reads the `odyssey_two_orbits` results under the viewer-demo output folder.

### Apollo 11 powered descent

`apollo11_landing.jl [--duration-s S] [--output-dir DIR]` flies the lunar module
from powered descent initiation to touchdown over a local terrain bundle:
quadratic descent guidance, RCS attitude control, the terrain-contact touchdown
event and the descent-engine plume model. The seventeen thruster levels and the
plume diagnostics are saved as columns, and the terrain is embedded in the
page. It needs the site bundle that `scripts/dev/terrain/fetch_moon_site.py`
writes (the header of `apollo11_landing.jl` gives the Apollo 11 fetch command,
with `--nac-half-deg 0.03`; `SPACEAGORA_TERRAIN_SITE` points at another site), the
tracked lunar module model, the LP165P coefficients and the SPICE starter
assets. It uses no atmosphere and no native GRAM. The touchdown event ends the
run, so `--duration-s` is only a cap. See the
[Lunar Landing](../../../docs/src/user/lunar_landing.md) guide. This
demonstrates the descent models on a documented case; it is not a
reconstruction of the flown trajectory or a claim of flight accuracy.

### Mission time and reference alignment

The event searches keep their estimated SPICE epoch. Each run converts its
start to the simulator's millisecond clock before it samples the initial state,
and the saved scene and SPICE reference use that same epoch. Differences later
in the run therefore compare the simulation with the navigation reconstruction,
not a sub-millisecond timestamp mismatch at the start. This does not validate
the atmosphere's native time or coordinate conventions. Apollo's lunar orbit
stays a nominal example, not a reconstruction from an Apollo flight kernel.

Full-length runs check that these examples work and that their clocks, initial
states and reference tables agree. They do not establish flight accuracy: the
saved separation from the navigation reference shows how far each run drifts,
and a quantitative reconstruction would need further model validation.

### Leap-second start times

The mission helpers refuse a start inside a UTC leap second, or close enough to
round into one at the simulator's millisecond resolution, with an explicit
error; choose a start outside it. Ordinary minute, day and calendar rollovers
work. The check is in the mission helpers only and adds no leap-second fields
to the general `InitialTime` API.

## CYGNSS research tools

The CYGNSS reconstruction study (in the SpaceAGORAMissionReconstructions
repository) owns its data pipeline, research drivers and tested helpers. Its
guide separates the Earthdata inputs from the private slew telemetry and
explains the limits of interpretation. The study reuses this folder's export
helpers rather than copying the renderer.
