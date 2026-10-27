# Atmosphere Models

This page helps you choose and configure a density model once you have run the
quickstart and want something other than its vacuum baseline. Every model
here except the native GRAM one works without NASA software. The two GRAM grid
models also need the public GRAMSuite Julia package and a payload file.

The quickstart itself uses `NoAtmosphereModel()`:

```text
julia --project=examples examples/AGORA_Basic_Quickstart.jl
```

`examples/AGORA_Mars_NoGRAM.jl` is the same kind of run with
`ExponentialAtmosphereModel(planet)` on Mars.

Related pages:

- [Assets & Modes](../assets.md)
- [Simulation Configuration](simulation_configuration.md)
- [GRAM Grid Surrogates](gram_grid_surrogates.md)
- [GRAMSuite Setup](gramsuite_setup.md)

The snippets below assume `using SpaceAGORA`.

## Choosing a model

| Model | Fidelity | What it needs | Notes |
|---|---|---|---|
| `NoAtmosphereModel()` | Vacuum | Nothing | Quickstart default |
| `ExponentialAtmosphereModel(planet)` | Low | Nothing | One scale height; fits near one altitude band |
| `PiecewiseExponentialAtmosphereModel(...)` | Low to medium | Nothing | Several layers; better shape over a wide altitude range |
| `NRLMSISE00AtmosphereModel(...)` | Medium | Nothing with fixed indices; internet access for live indices | Empirical model, documented for about 0 to 1000 km |
| `GRAMGridAtmosphereModel(...)` | Fixed snapshot | GRAMSuite with its grid API and a trusted grid payload | No native GRAM; evaluates within the grid's coverage |
| `GRAMNearSurfaceAtmosphereModel(...)` | Fixed snapshot | GRAMSuite with its near-surface API and a trusted near-surface payload | No native GRAM; Mars from 5 m above the surface to the payload's top. Format 1 payloads return zero wind, format 2 payloads their stored winds |
| `GRAMAtmosphereModel(...)` | High | The NASA GRAM Suite distribution, requested from NASA | Needs the GRAM setup described on [GRAMSuite Setup](gramsuite_setup.md) |

[GRAM Grid Surrogates](gram_grid_surrogates.md) describes the two grid models.
[GRAMSuite Setup](gramsuite_setup.md) covers the native model and how to set it
up.

## NoAtmosphereModel

The vacuum baseline. Density and wind are zero at every altitude, and the
temperature is the planet's reference temperature `T_ref`. Start with it for
orbit propagation and for checking an installation before you add
atmospheric effects.

```julia
density_model = NoAtmosphereModel()
```

## ExponentialAtmosphereModel

An analytic model with a single scale height, a constant temperature and zero
wind.

The planet form takes the reference density, reference altitude and scale
height from the planet object. `make_no_gram_planet` builds that object
without SPICE kernels and accepts `:earth`, `:mars` and `:venus`:

```julia
planet = make_no_gram_planet(:earth)
density_model = ExponentialAtmosphereModel(planet)
```

The general form takes the three parameters directly. Here it is with every
keyword at its default; put your own values in place of `ρ_ref`, `h_ref` and
`H`:

```julia
density_model = ExponentialAtmosphereModel(
    ρ_ref,   # reference density at h_ref, kg/m³
    h_ref,   # reference altitude, m
    H;       # scale height, m
    temperature_k = 200.0,          # constant temperature, K
    valid_min_altitude_m = h_ref,   # advisory lower bound, m
    valid_max_altitude_m = h_ref + 5*H  # advisory upper bound, m
)
```

The planet form accepts the same keywords. `H` must be positive, and the lower
advisory bound may not exceed the upper one.

The advisory band, `[h_ref, h_ref + 5H]` by default, is only a label. Outside
it the model keeps evaluating the same exponential; it neither clamps nor
raises an error. If your trajectory spans a wide altitude range, use
`PiecewiseExponentialAtmosphereModel` instead.

## PiecewiseExponentialAtmosphereModel

A multi-layer exponential atmosphere. Each layer has its own reference
density, reference altitude and scale height, which gives a better fit over a
wide altitude range and still needs no external data.

```julia
# Two layers; the values only illustrate the call and are not a fitted Earth profile
h_breaks_m = [80e3, 200e3, 600e3]    # N+1 breakpoints for N layers (strictly increasing), m
ρ_refs     = [6e-5,  2e-10]          # reference density per layer, kg/m³
Hs         = [7e3,   50e3]           # scale height per layer, m

density_model = PiecewiseExponentialAtmosphereModel(
    h_breaks_m,
    ρ_refs,
    Hs;
    h_refs = nothing,                  # optional; defaults to lower breakpoint of each layer
    temperature_k = 200.0,
    valid_min_altitude_m = 80e3,       # optional; defaults to h_breaks_m[1]
    valid_max_altitude_m = 600e3       # optional; defaults to h_breaks_m[end]
)
```

For `N` layers, give `N + 1` strictly increasing breakpoints, `N` reference
densities and `N` positive scale heights. If you pass `h_refs`, it also needs
`N` entries. The constructor rejects any other lengths.

Below the first breakpoint the bottom layer is extrapolated, and above the
last breakpoint the top layer is. Temperature is constant and wind is zero.

## NRLMSISE00AtmosphereModel

The NRLMSISE-00 empirical model, evaluated by
`SatelliteToolboxAtmosphericModels.jl`. Its documented range is about 0 to
1000 km; evaluation is not clamped to that range. Wind is always zero. The
model takes its solar and geomagnetic indices in one of three ways.

### Fixed geophysical indices

Fixed indices are the most reproducible choice, and the right one when a study
should hold solar activity constant. The values below are the defaults:

```julia
density_model = NRLMSISE00AtmosphereModel(
    f107a = 150.0,   # 81-day average solar flux (SFU)
    f107  = 150.0,   # previous-day solar flux (SFU)
    ap    = 4.0      # geomagnetic index (scalar or 7-element vector)
)
```

Solar-minimum conditions are roughly `f107a=70`, `f107=70`, `ap=2`, and
solar-maximum conditions roughly `f107a=220`, `f107=220`, `ap=15`.

### Live CelesTrak space indices

With `use_space_indices=true` the model looks up F10.7 and Ap from CelesTrak
data, through `SpaceIndices.jl`, for the date being simulated. The date comes
from the run's `initial_time` in `SimulationConfiguration` plus the elapsed
time. The first use needs internet access to fetch the dataset.

```julia
# Optional: prewarm the dataset before the solver starts
init_nrlmsise_space_indices!()

density_model = NRLMSISE00AtmosphereModel(use_space_indices=true)
```

Without the explicit call, the dataset is initialized on the first atmosphere
evaluation that needs it. Calling `init_nrlmsise_space_indices!()` before a
long run or a Monte Carlo campaign makes any download happen before the solver
starts. `init_nrlmsise_space_indices!(force_download=true)` refreshes the
dataset.

Below 80 km the built-in provider does not read the dataset at all; it returns
the standard fallback values `f107a=150`, `f107=150`, `ap=4`.

### Custom index provider

You can also pass any callable that returns `(f107a, f107, ap)` or a named
tuple with those keys. It must accept either `(instant)` or
`(instant, h, lat, lon)`; `instant` is a `DateTime`.

```julia
# Minimal provider: returns fixed indices regardless of time and position.
# Replace with a real lookup if solar activity varies in your study.
my_provider(::Any, ::Any, ::Any, ::Any) = (f107a=150.0, f107=150.0, ap=4.0)

density_model = NRLMSISE00AtmosphereModel(index_provider=my_provider)
```

`use_space_indices=true` and a custom `index_provider` cannot be combined; the
constructor raises an `ArgumentError`.

The model needs the scenario epoch to turn elapsed time into a date, so call
it through the seven-argument `getDensity(model, h, lat, lon, el_time, wind, p)`.
The six-argument form raises an `ArgumentError` for this model rather than
guess a date.

## GRAM-derived atmospheres

A frozen GRAM grid snapshot gives a GRAM-derived atmosphere without native
GRAM; see [GRAM Grid Surrogates](gram_grid_surrogates.md). The native
GRAM-backed model needs the NASA GRAM Suite; see
[GRAMSuite Setup](gramsuite_setup.md#GRAM-backed-atmosphere).
