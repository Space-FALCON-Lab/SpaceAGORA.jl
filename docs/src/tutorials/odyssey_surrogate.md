# Odyssey guidance and control with a frozen atmosphere

Use a surrogate when you want to work on your guidance, control or autonomy algorithm with a reproducible atmospheric environment. The example uses the named `odyssey_p20_frozen_v1` atmosphere, version `1.0.0`. It requires the Julia GRAMSuite package, but no native GRAM installation or GRAM data directory.

The example propagates Odyssey's P20 aerobraking passage and runs the existing maximum-energy-depletion guidance and panel controller. It compares a 90-degree panel-angle cap with a 30-degree cap. Both runs use the same atmosphere, initial state, force models and solver settings. The example checks that the controller applies each cap and that the trajectory and accumulated panel heating change.

## Run the example

From the SpaceAGORA checkout, install the dedicated example environment. It pins
the public Julia wrapper and uses the current SpaceAGORA checkout:

```sh
julia --project=examples/odyssey_surrogate_env examples/odyssey_surrogate_env/setup.jl
julia --project=examples/odyssey_surrogate_env examples/odyssey_surrogate.jl
```

For a notebook or an interactive session, select that same environment before
including the example:

```julia
include("examples/odyssey_surrogate.jl")
using .OdysseySurrogateExample

comparison = OdysseySurrogateExample.compare_panel_caps()
```


The preset resolver installs the identified atmosphere once, and the scenario helper supplies the four identified public SPICE kernels and Mars gravity coefficients. Subsequent runs reuse installed files. An offline run requires these assets to be installed already:

```sh
julia --project=examples/odyssey_surrogate_env examples/odyssey_surrogate.jl --offline --output=odyssey_offline_results
```

Each invocation needs a new output directory. A missing or incompatible asset produces an error. It does not start native GRAM or substitute another atmosphere.

The example replaces the process's SPICE kernel pool with its four pinned kernels.
Run it in a separate Julia session if your notebook also uses another SPICE scenario.

## Change an active setting

For a single run, change `panel_cap_deg`:

```julia
run = OdysseySurrogateExample.run_case(
    panel_cap_deg=45.0,
    output_dir="odyssey_cap_45_results",
)
```

The setting becomes `max_alpha_rad` in `AerobrakingEnergyDepletionConfig`. Guidance selects maximum energy depletion, and the controller commands the cap on spacecraft links 2 and 3 every 0.1 seconds. The guidance interval is 3 seconds. Aerodynamics and heating use the articulated panel geometry through `fixed_attitude_incidence=:attitude`; spacecraft attitude itself is not integrated.

This is a prescribed panel-cap exercise. Its thermal threshold is infinite, and it performs no targeting or heat-load prediction. It does not demonstrate closed-loop thermal protection, actuator rate limits or hardware dynamics. The reported heating comes from the production thermal model and integrated heat loads, rather than the controller's internal heat-rate estimate.

`configure_case` in the example shows how to replace the guidance and control configuration while retaining the atmosphere and scenario. The unchanged P20 regression remains a separate reference with guidance and control disabled. Results from this active-control exercise do not inherit a mission-accuracy claim from that reference.

## Inspect the result

The comparison writes:

- `cap_90_deg/trajectory.csv` and `cap_30_deg/trajectory.csv`: sampled height, latitude, longitude, density, temperature, panel commands, drag, panel heat rates, accumulated heat loads and Cartesian states.
- A `summary.toml` beside each table: preset and scenario provenance, resolved start time, solver result, measured runtime and state at a common elapsed time.
- `comparison.toml`: position, velocity and panel-heating differences at 600 seconds, while both trajectories remain in the atmosphere domain.

Each run stops at its outbound 250 km ellipsoidal-height crossing. The common-time comparison avoids confusing a control effect with a difference in the exit event's timing. The example's small nonzero effect checks establish that the selected setting is active. They are not accuracy tolerances or release requirements.

The cap is constant after activation. Historical force and heat-rate sampling in
this example therefore uses the same panel geometry. A time-varying controller
needs its panel geometry recorded alongside each state before applying this
postprocessing pattern.

Runtime in each summary covers the simulation and its sample recording. It excludes package loading and initial asset resolution. The first run may include compilation, so the two wall times are not a fair algorithm-performance comparison.

## Atmospheric and scenario assumptions

The grid stores density, temperature and east/north/up winds. Its atmosphere is frozen at `2001-11-07T11:51:04.794789Z`, with the recorded Mars-GRAM TES Mapping Year 2 configuration. Advancing simulation time does not change these fields. The example starts near `2001-11-07T11:45:45.690115Z`; the initial state is queried from the identified Odyssey kernel at the engine's resolved start time.

The stored domain is 100 to 260 km ellipsoidal altitude, 40 to 90 degrees geodetic north latitude, and periodic east longitude. The Mars equatorial and polar radii are 3,396,190 and 3,376,200 metres. A query outside the height or latitude bounds fails. The example uses nominal stored winds; `wind=false` is not a supported way to remove winds from this preset.

Exact-pole winds retain the native longitude-label convention, which does not define a unique physical horizontal vector there. The retained P20 passage avoids the exact pole, and the example records its maximum latitude. A new algorithm must keep its resulting trajectory inside the supported domain and address pole winds if it uses them.

Force models are degree/order 20 Mars gravity harmonics, Sun/Earth/Moon/Jupiter gravity, solar radiation pressure and free-molecular aerodynamics. The spacecraft geometry and mass come from the P20 reference configuration. Native-free does not mean asset-free: SPICE kernels supply ephemerides and Mars orientation, and the gravity coefficients supply the harmonic model.

The example selects the explicit `Tsit5` solver with a maximum step of 0.2 seconds,
orbit/atmosphere relative tolerance `1e-7` and absolute tolerance `1e-9`.
This is an independent usability exercise, not a replay of the reference's
automatic stiff-solver configuration.

This frozen snapshot is useful for repeatable algorithm development. It does not model other dates, uncertain forcing, atmospheric variability or the accuracy of a flight mission. For advanced atmospheric investigations, configure native GRAM explicitly, or generate and validate a new mission-specific snapshot before using it in repeated runs. Application-specific accuracy and robustness requirements remain separate from this example.
