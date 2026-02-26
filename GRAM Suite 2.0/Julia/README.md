# Julia FFI (GRAM C Interface)

This folder provides a minimal Julia wrapper around the existing GRAM C API (`Gram_C.h`) using `ccall`.

## 1. Build GRAM as a shared library

1. Install GNU Make (macOS):
   - `brew install make`
2. Confirm `Build/makefile.defs` includes `../makefile.defs.macos` (this repo is now configured that way)
3. Ensure a macOS CSPICE static archive exists and set `SPICE_LIB` if needed:
   - Apple Silicon default: `common/cspice/lib/cspice_macos_arm64.a`
   - Intel default: `common/cspice/lib/cspice_macos_x86_64.a`
4. Build:
   - `cd Build`
   - `gmake clean`
   - `gmake shared -j`

Expected output:
- macOS: `Build/lib/libGRAM.dylib`
- Linux: `Build/lib/libGRAM.so`
- Windows (MinGW): `Build/lib/libGRAM.dll`

## 2. Use from Julia

From the repository root:

```julia
include("Julia/GRAM.jl")
using .GRAM

libext = Sys.iswindows() ? "dll" : (Sys.isapple() ? "dylib" : "so")
set_library!(joinpath("Build", "lib", "libGRAM.$libext"))
initialize!("SPICE")

atmos = create_atmosphere(BODY_MARS; data_path=joinpath("Mars", "data"))
set_start_time!(atmos; year=2020, month=3, day=15, hour=0, minute=0, seconds=0.0, scale=1, frame=1)
set_position!(atmos; height=50.0, latitude=22.0, longitude=48.0, elapsed_time=100.0)

err = update!(atmos)
if err != 0
    error(get_error_message())
end

state = get_dynamics_state(atmos)
println("Temperature: ", state.temperature)
println("Density: ", state.density)

# Optional: generate a full trajectory track in one C call
track = generate_trajectory(
    atmos;
    initial_height=50.0,
    initial_latitude=22.0,
    initial_longitude=48.0,
    initial_elapsed_time=0.0,
    delta_height=0.5,
    delta_latitude=0.02,
    delta_longitude=0.05,
    delta_elapsed_time=10.0,
    n_points=16
)
println("track points: ", length(track))
println("first density: ", track[1].dynamics.density)

# Monte Carlo trajectories (different seeds)
tracks_mc = generate_monte_carlo_trajectories(
    atmos;
    initial_height=50.0,
    initial_latitude=22.0,
    initial_longitude=48.0,
    delta_height=0.5,
    delta_latitude=0.02,
    delta_longitude=0.05,
    delta_elapsed_time=10.0,
    n_points=16,
    n_runs=5,
    initial_seed=1001
)
println("mc runs: ", length(tracks_mc))

# Mars-specific controls
set_map_year!(atmos, 2)
set_mgcm_dust_levels!(atmos; constant_level=2.0, min_level=0.0, max_level=0.0)
set_dust_storm!(atmos;
    longitude_sun=250.0,
    duration=30.0,
    intensity=1.0,
    max_radius=20.0,
    latitude=10.0,
    longitude=120.0
)

close!(atmos)
```

## Notes

- The wrapper is intentionally small and covers the common generic API surface.
- Exposed advanced APIs include:
  - Generic C controls/state: `set_delta!`, `set_ephemeris_state!`, `get_density_state`, `get_gases_state`, `get_constituent_gas`, `get_ephemeris_state`, `get_perturbation_state`
  - Trajectory batch + Monte Carlo: `generate_trajectory`, `generate_monte_carlo_trajectories`
  - Mars controls/state: `set_map_year!`, `set_mgcm_dust_levels!`, `set_dust_storm!`, `set_dust_density!`, `set_f107!`, `set_wind_scales!`, `set_mola_heights!`, `set_min_max!`, `get_daily_dynamics_state`, `get_mars_state`, `get_mars_gases_state`
