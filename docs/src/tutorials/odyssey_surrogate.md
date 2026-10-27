# Odyssey guidance and control with a frozen atmosphere

A frozen atmosphere gives a guidance, control or autonomy algorithm the same
atmospheric environment on every run. This example uses the named preset
`odyssey_p20_frozen_v1`, version `1.0.0`. It needs the Julia GRAMSuite package
but no native GRAM installation and no GRAM data directory. That package is the
public Julia wrapper, which the setup step below installs; it is not NASA's GRAM
Suite 2.0 software, which only the native GRAM backend needs.

The example propagates Odyssey's P20 aerobraking passage with the
maximum-energy-depletion guidance and its panel controller, once with a
90-degree panel-angle cap and once with a 30-degree cap. Both runs use the same
atmosphere, initial state, force models and solver settings. The example checks
that the controller applied each cap and that the trajectory and the
accumulated panel heating differ between the two runs.

## Run the example

All paths on this page are relative to the repository root. Blocks marked
"Terminal" are typed at the shell prompt there (PowerShell or Command Prompt on
Windows, a terminal on Linux or macOS); blocks marked "Julia REPL" go to the
`julia>` prompt. First install the example's own environment, which pins the
public Julia wrapper and uses this SpaceAGORA checkout, then run the
comparison. Terminal:

```sh
julia --project=examples/odyssey_surrogate_env examples/odyssey_surrogate_env/setup.jl
julia --project=examples/odyssey_surrogate_env examples/odyssey_surrogate.jl
```

The setup script changes two things in the checkout. It clones the EDG guidance
companion from GitHub into `packages/gnc/SpaceAGORAEDG.jl`, a Git submodule, so
the first setup needs network access. It then adds that companion to the
example environment, which rewrites the tracked file
`examples/odyssey_surrogate_env/Project.toml`; `git status` lists that file as
modified, which is expected. The resolved `Manifest.toml` beside it is ignored
by Git. You don't need to instantiate the root project first: this environment
resolves its own manifest, SpaceAGORA's dependencies included.

After pulling a newer SpaceAGORA, run the setup script again. It moves the
environment to the GRAMSuite revision the checkout now pins, which a plain
`Pkg.instantiate()` of the earlier setup would not do.

The first setup downloads and compiles several hundred packages, so it takes a
while; later setups reuse them. On a machine with an empty Julia depot, expect
downloads on the order of a gigabyte, and a depot (`~/.julia` by default) that
grows by a few gigabytes once the compiled caches and this example's assets are
in place.

Each command-line run starts a new Julia process, which compiles the example
before its first passage, so expect a pause. The first command in a new
environment, even `--help`, also loads and compiles its packages. The example
prints a reminder of this before the first passage. The run time it reports for
that passage, which is also `wall_s` in that case's `summary.toml`, includes the
compilation, so the baseline case looks much slower than the variant after it;
the difference is compilation, not the panel angle. An interactive session
compiles once and reuses the code for later runs. When the comparison finishes,
it prints the absolute path of its results directory and the files it wrote
there, which [Inspect the result](#Inspect-the-result) describes.

In a notebook or an interactive session, select the same environment before
including the example: start Julia from the repository root with
`julia --project=examples/odyssey_surrogate_env`, or run
`] activate examples/odyssey_surrogate_env` in a session started there. The
`include` path is relative to the repository root. Julia REPL:

```julia
include("examples/odyssey_surrogate.jl")
using .OdysseySurrogateExample

comparison = OdysseySurrogateExample.compare_panel_caps(output_dir="output/odyssey_notebook");
```

The trailing semicolon keeps the REPL from printing the returned tables. Give
each run a new directory; the examples on this page use directories under
`output/`, which Git ignores.

On the first run, the preset resolver installs the atmosphere grid (47 MB) and
the scenario helper installs the four public SPICE kernels and the Mars gravity
coefficients (176 MB). Each installation prints one line naming its source and
one confirming that the SHA256 checksums match. The files go into the
`artifacts/` folder of your Julia depot, by default `~/.julia/artifacts`, and
later runs reuse them. The atmosphere grid is published under CC BY 4.0: cite
it by preset name and version, and credit NASA's GRAM Suite (Mars-GRAM) as the
source model. The scenario assets keep the original terms of the NAIF kernels
and gravity coefficients (see the
[GRAMSuite.jl repository](https://github.com/Space-FALCON-Lab/GRAMSuite.jl)).
To list, prefetch or check presets without running the example, use the
commands in [Named surrogate data](../cli.md#Named-surrogate-data).

Once the assets are installed, a run can work offline. Terminal:

```sh
julia --project=examples/odyssey_surrogate_env examples/odyssey_surrogate.jl --offline --output=output/odyssey_offline
```

Each run needs a new output directory. An existing one is refused, so earlier
results are never overwritten; pass `--output=DIR` on the command line, or
`output_dir="DIR"` in Julia. A missing or incompatible asset is an error: the
example never starts native GRAM or substitutes another atmosphere.

The example replaces the process's SPICE kernel pool with its four pinned
kernels, so run it in a separate Julia session if your notebook also uses
another SPICE scenario.

## Change an active setting

For a single run, change `panel_cap_deg`. Julia REPL, in the example
environment:

```julia
run = OdysseySurrogateExample.run_case(
    panel_cap_deg=45.0,
    output_dir="output/odyssey_cap_45",
);
```

`run_case` returns a named tuple: `summary` is the dictionary written to
`summary.toml`, `rows` holds the samples written to `trajectory.csv`, `common`
is the sample at the 600-second comparison time and `endpoint` is the sample at
the outbound 250 km crossing.

To compare another cap with the 90-degree baseline, pass it on the command
line; `--help` lists all options. Terminal:

```sh
julia --project=examples/odyssey_surrogate_env examples/odyssey_surrogate.jl --help
julia --project=examples/odyssey_surrogate_env examples/odyssey_surrogate.jl --cap=60 --output=output/odyssey_cap_60
```

The second command runs the 90-degree baseline and a 60-degree variant and
writes `cap_90_deg/`, `cap_60_deg/` and `comparison.toml` in
`output/odyssey_cap_60/`. `--baseline-cap=DEG` changes the baseline. Both caps
must be in (0, 90] degrees and must differ. The same comparison from the Julia
REPL, in its own directory:

```julia
comparison = OdysseySurrogateExample.compare_panel_caps(
    variant_cap_deg=60.0,
    output_dir="output/odyssey_cap_60_julia",
);
```

The cap becomes `max_alpha_rad` in `AerobrakingEnergyDepletionConfig`. The
guidance runs every 3 seconds and selects maximum energy depletion; the
controller commands the cap on spacecraft links 2 and 3 every 0.1 seconds. EDG
computes the prediction, mode, switch times and constrained angle; only the
control callback moves the panels, so a guidance decision by itself changes
neither the panel geometry nor the heat loads. Aerodynamics and heating use the
articulated panel geometry through `fixed_attitude_incidence=:attitude`; the
spacecraft attitude itself is not integrated.

This is a prescribed panel-cap exercise. Its heat-rate limit is infinite, and it
does no targeting or heat-load prediction, so it does not demonstrate
closed-loop thermal protection, actuator rate limits or hardware dynamics. The
reported heating comes from the production thermal model and the integrated
heat loads, not from the controller's internal heat-rate estimate.

`configure_case` in the example shows how to replace the guidance and control
configuration while keeping the atmosphere and the scenario. The P20
regression case, with guidance and control disabled, stays a separate
reference; results from this active-control exercise don't inherit its
mission-accuracy claim.

## EDG state across runs and passes

Each `run_simulation` starts EDG from a fresh state, and the paired guidance
and controller share one `AerobrakingEnergyDepletionState` with one entry per
spacecraft and matching configurations. EDG refuses checkpoint writing and
resume, because a checkpoint does not hold its pass state. The example
registers its outbound 250 km exit event with `EDGGeodeticBoundaryCondition`,
so EDG closes the passage, records the exit heat and clears its plan before the
run stops; the `endpoint` sample's mode is `inactive`.

After an exit, the next control update uses the fresh-pass default: the minimum
angle when only targeting is configured, the maximum angle when maximum
depletion is configured. That command still matters between passes when
density, solar pressure or attitude coupling act above the entry interface.
Because EDG asks for atmosphere crossing events even when both solver phases
have the same settings, the added root-finding can change the integration path:
unchanged guidance formulas do not imply bit-identical trajectories.
[State between passages and repeated runs](../user/edg_guidance.md#State-between-passages-and-repeated-runs)
in the EDG guide gives the full rules.

## Inspect the result

Unless you pass `--output=DIR` (or `output_dir` in Julia), the comparison
writes to `odyssey_surrogate_results/` in the current directory, and `run_case`
writes to `odyssey_surrogate_<cap>/`, for example `odyssey_surrogate_45.0/`.
Both print the absolute path when they finish. The comparison directory
contains:

- `cap_90_deg/trajectory.csv` and `cap_30_deg/trajectory.csv`, named after the compared caps: sampled height, latitude, longitude, density, temperature, panel commands, drag, panel heat rates, accumulated heat loads and Cartesian states.
- A `summary.toml` beside each table: preset and scenario provenance, the resolved start time, the solver result, the measured run time and the state at the common comparison time.
- `comparison.toml`: the baseline and variant caps (`baseline_cap_deg`, `variant_cap_deg`) and the position, velocity and panel-heating differences at 600 seconds, while both trajectories are still inside the atmosphere's domain.

The three differences are lengths of difference vectors at 600 seconds, so they
are never negative and carry no direction. `position_difference_m` is the
distance in metres between the two inertial positions, `velocity_difference_m_s`
the magnitude in m/s of the difference of the two velocities, and
`panel_heat_load_difference_norm_J_cm2` the Euclidean norm in J/cm² of the
per-panel differences between the two runs' accumulated heat loads. To see which
run is ahead or hotter, compare the signed columns of the two `trajectory.csv`
files, for example by plotting them. The example's own environment does not
include Plots, so use the repository's examples environment, prepared as in
[Installation & Environment](../user/installation_environment.md#Examples-environment).
Julia REPL started with `julia --project=examples` from the repository root,
after the default comparison:

```julia
using CSV, DataFrames, Plots
dir = "odyssey_surrogate_results"
baseline = CSV.read(joinpath(dir, "cap_90_deg", "trajectory.csv"), DataFrame)
variant = CSV.read(joinpath(dir, "cap_30_deg", "trajectory.csv"), DataFrame)
drag = plot(baseline.elapsed_s, baseline.drag_N; label="cap 90 deg",
    xlabel="Elapsed time (s)", ylabel="Drag (N)")
plot!(drag, variant.elapsed_s, variant.drag_N; label="cap 30 deg")
height = plot(baseline.elapsed_s, baseline.height_km; label="cap 90 deg",
    xlabel="Elapsed time (s)", ylabel="Ellipsoidal height (km)")
plot!(height, variant.elapsed_s, variant.height_km; label="cap 30 deg")
savefig(plot(drag, height; layout=(2, 1), size=(800, 700)),
    joinpath(dir, "passages.png"))
```

For the default comparison the saved figure looks like this:

![Two plots against elapsed time for panel caps of 90 and 30 degrees. Drag peaks at periapsis, at about 14 N with the 90-degree cap and about 9 N with the 30-degree cap. The two height curves, from 250 km down to about 111 km and back, lie on top of each other.](../assets/figures/odyssey_panel_cap_comparison.png)

For another pair of caps, change the directory and the two `cap_*_deg` folder
names to the ones your run printed.

The Cartesian states are Mars-centred J2000 inertial position in metres and
velocity in m/s: the columns `x_m`, `y_m`, `z_m`, `vx_m_s`, `vy_m_s` and
`vz_m_s` in `trajectory.csv`, and `initial_state_m_m_s` and
`common_state_m_m_s` in `summary.toml`. This is the frame of the SPICE initial
state, and each `summary.toml` records it as `state_frame`. In the same table,
`height_km` is ellipsoidal height, `latitude_deg` geodetic latitude and
`longitude_deg` east longitude, the conventions of the preset's domain.

The atmosphere provenance in `summary.toml` describes how the published preset
was made. Its `generation_provenance.runtime_wrapper_revision` is the GRAMSuite
revision used to generate and publish the grid, not necessarily the one running
now. The revision your session runs is pinned in
`examples/odyssey_surrogate_env/Project.toml` and recorded in that folder's
`Manifest.toml`.

`summary.toml` carries three timestamps that are easy to confuse.
`requested_start_utc` and `start_utc`, at the top of the file, are the start of
the simulation: the instant at which the initial state is read from the Odyssey
kernel. Under `[atmosphere.atmosphere]`, `epoch_utc` is the frozen epoch of the
atmosphere snapshot, the instant whose Mars-GRAM state the grid stores;
advancing simulation time never changes it. Under
`[atmosphere.generation_provenance]`, `payload_created_at_utc` is when the grid
file was generated and packaged for publication. It describes the data product,
not the mission, and it is the only one of the three that is not a 2001 date.

Each run stops at its outbound 250 km ellipsoidal-height crossing, so
`solver_retcode = "Terminated"` in `summary.toml` is the normal result.
Comparing at a common time keeps a control effect apart from a difference in the
timing of the exit event. The example's checks for a small nonzero effect only
establish that the selected setting is active; they are not accuracy
tolerances or release requirements.

The cap is constant once active, so the force and heat-rate values sampled after
the fact use the same panel geometry as the run. A time-varying controller needs
its panel geometry recorded with each state before it can reuse this
postprocessing.

The run time in each summary covers the simulation and its sample recording; it
excludes package loading and the first asset resolution. Because the first run
can include compilation, the two run times are not a fair comparison of the
algorithm's speed.

To query the frozen atmosphere directly, for example along your own trajectory,
load both packages in the example environment, build the model and call
`getDensity`. Julia REPL:

```julia
using SpaceAGORA
import GRAMSuite  # the public wrapper that provides the grid atmosphere

model = surrogate_preset_model("odyssey_p20_frozen_v1"; version="1.0.0")
rho, temperature, wind = SpaceAGORA.getDensity(model, 150e3, deg2rad(60.0), deg2rad(30.0), 0.0, true)
```

The arguments are ellipsoidal height in metres, geodetic latitude and east
longitude in radians, elapsed time in seconds and a `Bool` wind flag. The frozen
snapshot ignores the elapsed time, and the flag does not remove its stored
winds. The call returns density in kg/m³, temperature in kelvin and the east,
north and up wind components in m/s.
[GRAM Grid Surrogates](../user/gram_grid_surrogates.md#Fixed-GRAM-grid-snapshot)
gives the domain policy, and [Extensibility](../extensibility.md) the full
interface.

## Atmospheric and scenario assumptions

The grid stores density, temperature and east, north and up winds. The
atmosphere is frozen at `2001-11-07T11:51:04.794789Z` (`epoch_utc` in
`summary.toml`), with Mars-GRAM configured for TES Mapping Year 2; advancing
simulation time does not change these fields. The example starts at
`2001-11-07T11:45:45.690115Z`, and the initial state is read from the Odyssey
kernel at the engine's resolved start time.

The stored domain is 100 to 260 km ellipsoidal altitude, 40 to 90 degrees
geodetic north latitude, and all east longitudes. The Mars equatorial and polar
radii are 3,396,190 and 3,376,200 metres. A query outside the height or latitude
bounds is an error. A named preset fixes this policy: `surrogate_preset_model`
rejects grid options such as `above_grid=:vacuum` with an `ArgumentError`.
[GRAM Grid Surrogates](../user/gram_grid_surrogates.md#Fixed-GRAM-grid-snapshot)
shows the generic grid model for other policies; a query below the floor fails
with either model. The example uses the stored nominal winds; `wind=false` is
not a supported way to remove winds from this preset.

At the exact pole, the stored winds keep the native longitude-label convention,
which does not define a unique physical horizontal vector there. The P20
passage avoids the exact pole, and the example records its maximum latitude. A
new algorithm must keep its trajectory inside the supported domain, and must
handle the pole winds if it uses them.

The force models are degree and order 20 Mars gravity harmonics, Sun, Earth,
Moon and Jupiter gravity, solar radiation pressure and free-molecular
aerodynamics. The spacecraft geometry and mass come from the P20 reference
configuration. Native-free does not mean asset-free: the SPICE kernels supply
the ephemerides and Mars orientation, and the gravity coefficients supply the
harmonic model.

The example selects the explicit `Tsit5` solver with a maximum step of
0.2 seconds, relative tolerance `1e-7` and absolute tolerance `1e-9` in both the
orbit and the atmosphere phases. It is a usability exercise in its own right,
not a replay of the reference case's automatic stiff-solver configuration.

The frozen snapshot suits repeatable algorithm development. It does not model
other dates, uncertain forcing, atmospheric variability or the accuracy of a
flight mission. For deeper atmospheric work, configure native GRAM explicitly
([Compare a surrogate with native GRAM](../user/gram_grid_surrogates.md#Compare-a-surrogate-with-native-GRAM)),
or generate and validate a new mission-specific snapshot before using it in
repeated runs
([Prepare another frozen snapshot](../user/gram_grid_surrogates.md#Prepare-another-frozen-snapshot)).
Accuracy and robustness requirements for a particular application are separate
from this example.
