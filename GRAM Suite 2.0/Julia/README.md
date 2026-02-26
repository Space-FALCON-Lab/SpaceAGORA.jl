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

close!(atmos)
```

## Notes

- The wrapper is intentionally small and covers the common generic API surface.
- You can extend it by adding more structs/functions from `common/include/Gram_C.h` and the planet-specific `*_Atmosphere_C.h` headers.
