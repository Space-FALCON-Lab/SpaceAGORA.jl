# Local viewer demonstrations

Run these scripts from the repository with its Julia environment instantiated.
The four basic demonstrations below use analytic ephemerides and either no
atmosphere or an exponential model.
Native GRAM, SPICE kernels and private telemetry are not required. The ISS display
uses the tracked NASA model in `data/models/`; it is display geometry, not the
physical mass or aerodynamic mesh.

| Script | Default simulated time | What it demonstrates |
| --- | --- | --- |
| `iss_demo.jl` | 11,120 s | A J2 orbit with the NASA ISS display model |
| `arm_demo.jl` | 60 s | Planned, controlled robot-arm motion on an orbiting bus |
| `earth_4day.jl` | Four days | A long, eccentric Earth trajectory |
| `odyssey_two_orbits.jl` | Two initial Keplerian periods | Mars energy-depletion control with an exponential atmosphere |

For a short first run:

```sh
julia --project=. scripts/dev/viewer_demos/iss_demo.jl --duration-s 60 --output-dir output/iss-first-run
```

All four scripts accept `--duration-s` and `--output-dir` with space-separated
values, as in the command above. These scripts do not accept the `--key=value`
spelling used by the main CLI. Without an explicit
output directory, they use `output/viewer_demos/<script-name>`; the environment
variable `SPACEAGORA_VIEWER_DEMO_OUT` changes that parent directory. An existing
nonempty output directory is rejected, so a new run cannot erase earlier work.
Choose a new directory for another run.

Each run writes simulation results, its scene description and a self-contained HTML
page into that directory, then prints the HTML path. Open the HTML after the
simulation finishes. The viewer replays saved results. A short run checks the
setup and export, but does not demonstrate a full orbit, aerobraking passage or
completed robot-arm maneuver. The Mars scenario is a demonstration of the
analytic model and controller, not a fit to flight telemetry; its printed minimum
altitude is the minimum among saved samples.

The Odyssey demonstration explicitly uses Tsit5 and budgets solver steps for its
0.1-second controller over the requested duration, with a margin for rejected
steps. This changes the iteration ceiling, not the integration tolerances.

## Mission demonstrations and SPICE comparisons

`apollo11_lunar_orbit.jl`, `magellan_aerobraking.jl`, `odyssey_aerobraking.jl`, and
`cassini_titan_flyby.jl [TA|T5|all]` retain the PR121 mission cases. The first
three default to two modeled orbits; each Cassini case defaults to five hours
starting 2.5 hours before the searched closest approach. All accept
`--duration-s <positive seconds>` and `--output-dir <fresh directory>` for a
bounded run. Existing nonempty outputs are refused, never reused or removed.
`SPACEAGORA_DEMO_FORCE` no longer overrides that protection.

These drivers use the tracked mission display models in `data/models/` and
require existing SPICE kernels (`SPACEAGORA_SPICE_PATH`). Magellan, Odyssey and
Cassini require installed native GRAM support; each wrapper is constructed explicitly
at the simulation's initial epoch. This does not validate native time-system,
coordinate, datum or climatology assumptions. Mission navigation SPKs are fetched
from NAIF only when missing, into `SPACEAGORA_MISSION_SPK_DIR` or the configured
SPICE tree's `spk/missions` directory. Native GRAM itself is never downloaded.

Odyssey requests TES mapping year 2 with `mars_map_year=2`. In the tested native
wrapper, that request alone does not apply the year selection to the native
model; a later parameter-setting call is needed. This example preserves its
existing numerical configuration, so do not interpret it as a validated
year-2 atmosphere. An explicitly applied profile and its effect on the
trajectory need separate validation.

References are central-body-relative, geometric J2000 SPICE states in SI units
at the saved simulation times. Missing SPK coverage is reported and omitted.
The target is a one-based spacecraft index, not its public id. Separation is a
sampled model comparison, not mission reconstruction accuracy. Apollo uses an
illustrative parking orbit with SPICE orientation, not a flown LM reference.
Without integrated attitude these examples use the viewer's velocity-aligned
pose; a displayed NASA model is not a reconstructed attitude history.

Titan's static zero-tide 5x5 field uses fully normalized coefficients, a
2,575,000 m reference radius and GM 8,978,126,919,238.97 m^3/s^2 from NASA PGDA
product 91. No time-varying body tide is implied. Odyssey's optional
`SPACEAGORA_DEMO_ODYSSEY_ATMOSPHERE=accelerometer` requires a separately supplied
local density table; it is not silently substituted for GRAM.
`SPACEAGORA_DEMO_MAGELLAN_ANTENNA=forward|aft` controls the model pose.
`SPACEAGORA_DEMO_MESH_AERO=1` enables SpaceAGORA's built-in mesh aerodynamic model
and fits the selected geometry afresh; it does not reuse a stale mesh fit.

The canonical output is standalone `simulation_results_viewer.html`. Optional
`SPACEAGORA_DEMO_CDN=1` invokes the separate `build_cdn_page.py` tool when present.
That network-dependent sharing feature remains separate from these drivers.
`check_odyssey.jl [results_directory]` inspects saved articulation and altitude
from the basic energy-depletion Odyssey demonstration without propagating it.
Without an argument, it reads the `odyssey_two_orbits` results directory under
the configured viewer-demo output parent.

### Mission time and reference alignment

The event searches retain their estimated SPICE epoch. Each run converts its
start to the simulator's millisecond clock before sampling the initial state;
the saved scene and SPICE reference use that same epoch. Differences later in
the run measure the selected simulation against the navigation reconstruction,
not an initial sub-millisecond timestamp mismatch. This does not validate the
atmosphere's native time or coordinate conventions. Apollo's lunar orbit remains
a nominal example, not a reconstruction from an Apollo flight kernel.

Full-duration validation checks that these examples run and that their clocks,
initial states and reference tables agree. It does not establish flight accuracy.
For example, the default native-GRAM Odyssey run reaches about 552 km separation
from its navigation reference over 34.6 hours, with the second periapsis about
120 seconds late. Quantitative reconstruction requires further model validation.
