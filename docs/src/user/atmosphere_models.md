# Atmosphere Models

Use this page when you need to choose and configure an atmosphere model for a
no-GRAM or open-data simulation.

This page is for users who have already completed the quickstart and want to
move beyond the vacuum baseline to a higher-fidelity open-data atmosphere.

Shortest successful command:

```text
julia --project=. examples/AGORA_Basic_Quickstart.jl
```

What to read next:

- [Assets & Modes](../assets.md)
- [Simulation Configuration](simulation_configuration.md)
- [GRAMSuite Setup](gramsuite_setup.md)

## Choosing a model

| Model | Fidelity | Assets required | Notes |
|---|---|---|---|
| `NoAtmosphereModel()` | Vacuum | None | Default quickstart baseline |
| `ExponentialAtmosphereModel(planet)` | Low | None | Single scale height; valid near one altitude band |
| `PiecewiseExponentialAtmosphereModel(...)` | Low–medium | None | Multi-layer; better altitude-shape fit |
| `NRLMSISE00AtmosphereModel(...)` | Medium | None (fixed indices) or internet (live indices) | Standard empirical model; ~0–1000 km |
| `GRAMGridAtmosphereModel(...)` | Fixed snapshot | GRAMSuite with its grid API and a trusted grid payload | Native-free evaluation within documented grid coverage |
| `GRAMAtmosphereModel(...)` | High | Licensed NASA GRAM | Requires GRAM asset setup |

For GRAM setup, see [GRAMSuite Setup](gramsuite_setup.md).

---

## NoAtmosphereModel

Vacuum baseline. Returns zero density and zero wind at every altitude. This is
the recommended starting model for orbit propagation and for validating
installation before introducing atmosphere effects.

```julia
density_model = NoAtmosphereModel()
```

---

## ExponentialAtmosphereModel

Single-scale-height analytic model with constant temperature and zero wind.

**Planet-default constructor** — uses reference values from the built-in planet
object:

```julia
planet = make_no_gram_planet(:earth)
density_model = ExponentialAtmosphereModel(planet)
```

Supported planets: `:earth`, `:mars`, `:venus`.

**Manual constructor** — specify the parameters directly:

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

**Important:** the model evaluates the same exponential function outside the
advisory valid band — it does not clamp or error. The default advisory range is
`[h_ref, h_ref + 5H]`. Use `PiecewiseExponentialAtmosphereModel` if your
trajectory spans a wide altitude range.

---

## PiecewiseExponentialAtmosphereModel

Multi-layer exponential atmosphere. Each layer has its own reference density,
reference altitude, and scale height. This gives a better fit over a wide
altitude range while still requiring no external data.

```julia
# Example: two-layer Earth model
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

Array lengths must satisfy: `length(h_breaks_m) == N+1`, `length(ρ_refs) == N`,
`length(Hs) == N` where `N` is the number of layers. Breakpoints must be
strictly increasing.

For altitudes below the first breakpoint or above the last, the nearest layer
is used for extrapolation.

---

## NRLMSISE00AtmosphereModel

Standard empirical atmosphere model backed by `SatelliteToolboxAtmosphericModels.jl`. Covers
approximately 0–1000 km altitude. Accepts three input modes:

### Fixed geophysical indices

Simplest and most reproducible — use for studies where solar activity should
be held constant:

```julia
density_model = NRLMSISE00AtmosphereModel(
    f107a = 150.0,   # 81-day average solar flux (SFU)
    f107  = 150.0,   # previous-day solar flux (SFU)
    ap    = 4.0      # geomagnetic index (scalar or 7-element vector)
)
```

Solar-minimum conditions are approximately `f107a=70`, `f107=70`, `ap=2`.
Solar-maximum conditions are approximately `f107a=220`, `f107=220`, `ap=15`.

### Live CelesTrak space indices

Fetches real F10.7 and Ap data from CelesTrak. Requires internet access on
first use. The `InitialTime` epoch in `SimulationConfiguration` is used to
look up the correct date:

```julia
# Optional: prewarm the dataset before the solver starts
init_nrlmsise_space_indices!()

density_model = NRLMSISE00AtmosphereModel(use_space_indices=true)
```

If you skip `init_nrlmsise_space_indices!()`, the dataset initializes lazily on
the first atmosphere evaluation. Call it explicitly before long runs or Monte
Carlo campaigns so any download happens before the solver starts.

Below 80 km, the built-in provider returns the standard fallback values
`f107a=150 / f107=150 / ap=4` without touching the dataset.

### Custom index provider

Provide a callable that returns `(f107a, f107, ap)` or a named tuple with
those keys. The callable must accept either `(instant)` or
`(instant, h, lat, lon)`:

```julia
# Minimal provider: returns fixed indices regardless of time and position.
# Replace with a real lookup if solar activity varies in your study.
my_provider(::Any, ::Any, ::Any, ::Any) = (f107a=150.0, f107=150.0, ap=4.0)

density_model = NRLMSISE00AtmosphereModel(index_provider=my_provider)
```

`use_space_indices=true` and a custom `index_provider` cannot be combined.

---

## Fixed GRAM grid snapshot

For a named, automatically retrieved atmosphere, start with the [Odyssey surrogate workflow](../tutorials/odyssey_surrogate.md). Use `surrogate_preset_model("odyssey_p20_frozen_v1"; version="1.0.0")` after loading `GRAMSuite`; retrieval and verification occur once before the solver.

`GRAMGridAtmosphereModel` connects GRAMSuite's existing offline interpolation
kernel to SpaceAGORA. It needs the Julia wrapper with its native-free grid API
and a trusted serialized grid payload. Construction and density evaluation use
no native GRAM installation or native fallback.

```julia
using SpaceAGORA
import GRAMSuite

density_model = GRAMGridAtmosphereModel(
    planet="earth",
    surrogate_file="/path/to/authorized/earth_surrogate.jls",
    above_grid=:error,
)
```

The constructor forwards grid-loader options, including `search_roots` and
`expected_sha256`. Missing files, undownloaded Git LFS pointers, and invalid
payloads produce errors. A checksum pins file contents; it does not establish
scientific provenance or permission to distribute the file. Grid distribution
and an independently reproducible public installation remain separate work.

Queries use altitude in metres and latitude/longitude in radians, on the grid's
documented reference surfaces. Results contain density in kg/m³, temperature in
kelvin, and local east/north/up winds in m/s. The stored atmosphere is frozen:
elapsed time and the wind selector do not alter its density or stored winds.
Select and validate the atmospheric epoch, coordinates, forcing, and coverage
for the mission before use. Legacy metadata may leave these facts unknown;
loading a file cannot recover them or establish physical accuracy.

The default policy rejects altitude and latitude outside the grid. Longitude is
periodic. Explicit `above_grid=:vacuum` returns zero density and wind above the
ceiling, with `vacuum_temperature` in kelvin; lower-bound extrapolation remains
an error. The adapter does not apply the native model's entry-interface
polynomial, fixed 2000 km cutoff, or lower-altitude clamp.

The grid model bypasses `SPACEAGORA_VACUUM_GRAM_CACHE` and
`SPACEAGORA_DENSITY_FREEZE_PER_STEP`, and it reevaluates coordinates even when a
buffered sample has the same timestamp. Native GRAM track caches, isolated pools,
and per-satellite native copies do not apply. Shared read-only queries support
threaded evaluation; keep the arrays and metadata unchanged during a run.
Ordinary `deepcopy` produces independent arrays, including when configuration
isolation copies an entire run. Managed process-worker startup has separate
native warm-up behavior and is outside this native-free adapter's scope.

---

## GRAM-backed atmosphere

When GRAMSuite assets are available, the GRAM-backed constructor is provided by
the `SpaceAGORAGRAMSuiteExt` extension and activated by importing `GRAMSuite`.
See [GRAMSuite Setup](gramsuite_setup.md) for the full setup walkthrough.

The basic usage once assets are in place:

```julia
setup_gram_example!()   # from examples/common.jl — loads GRAMSuite extension
density_model = GRAMAtmosphereModel(planet_name="earth")
```

`GRAMAtmosphereModel` is not exported from the root `SpaceAGORA` module. It is
available after `GRAMSuite` is loaded and is accessed through
`SpaceAGORA.SimulationModel.GRAMAtmosphereModel` (or via `setup_gram_example!`
in the examples).

### GRAM epoch alignment

A GRAM model carries its own epoch: the `initial_time` keyword of
`GRAMAtmosphereModel` (the GRAMSuite default when omitted). The simulation
carries another in `SimulationConfiguration.initial_time`. Before every run,
`run_simulation` passes the density model through
[`with_density_model_epoch`](@ref) so the two agree. This happens before the
configuration is copied for state isolation, and the caller's configuration is
left unchanged.

For a keyword-built `GRAMAtmosphereModel` the rule is:

- An unchanged epoch returns the same model object. Epochs are compared on
  their integer year, month, day, hour and minute and on the seconds as
  `Float32`, with no time tolerance.
- A changed epoch rebuilds a fresh native model from the recorded constructor
  keywords with only `initial_time` replaced; planet, paths, perturbation
  options and every other recorded keyword are preserved. The rebuild runs
  under the GRAM setup lock.
- A `GRAMAtmosphereModel` wrapped around a raw GRAMSuite core, whose
  constructor keywords are unknown, cannot be realigned and raises an
  `ArgumentError` when the epochs differ. Construct it with the keyword
  constructor at the run's `initial_time` instead.

A `GRAMAtmosphereModelSurrogate` built over a GRAM base keeps its object, file
and fallback setting for an unchanged epoch and rejects a changed one with an
`ArgumentError`: its table was generated for one epoch, and rebuilding the
native fallback would not re-epoch the table. Build and validate a surrogate
for the epoch you intend to run. Surrogates over other base models, the grid
model, `NRLMSISE00AtmosphereModel` and the analytic models are not touched by
the hook.

To avoid a rebuild, construct the model at the run's epoch:

```julia
density_model = GRAMAtmosphereModel(planet_name="earth", initial_time=args.initial_time)
```

The hook changes the construction epoch only. It does not convert between time
systems, validate atmospheric coordinate or datum conventions, certify cached
or surrogate data, or carry over a native random stream that was advanced or a
handle that was edited by hand before the run. The existing native cache and
environment policies apply to the rebuilt model as to any other.
