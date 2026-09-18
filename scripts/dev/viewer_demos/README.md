# Local viewer demonstrations

Run these scripts from the repository with its Julia environment instantiated.
They use analytic ephemerides and either no atmosphere or an exponential model.
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
