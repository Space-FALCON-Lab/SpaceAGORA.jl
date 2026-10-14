const STUDY_ROOT = @__DIR__
const STUDY_PROJECT = joinpath(STUDY_ROOT, "Project.toml")

if something(Base.active_project(), "") != STUDY_PROJECT
    import Pkg
    Pkg.activate(STUDY_ROOT; io=devnull)
end

using Dates
using Printf
using Statistics

using SpaceAGORA
import GRAMSuite

const SM = SpaceAGORA.SimulationModel
const RuntimeServices = SpaceAGORA.RuntimeServices

include(joinpath(STUDY_ROOT, "corridor.jl"))
include(joinpath(STUDY_ROOT, "gram_prior.jl"))
include(joinpath(STUDY_ROOT, "merra2.jl"))
include(joinpath(STUDY_ROOT, "merra2_native.jl"))

# Checks the native reader against real granules. The synthetic-granule tests in
# runtests.jl pin the reader's logic, but two things can only be settled on real
# data:
#
#   1. Whether `H` is geopotential or geometric height. The reader assumes
#      geopotential and converts; the difference is about 1 km at 80 km, which
#      matters against an 18 km vertical correlation scale. Check 3 below prints
#      both interpretations against GRAM's own geometric-altitude profile.
#   2. Whether the level ordering in a *subsetted* file survived. The MERRA-2
#      spec warns that some post-processors flip the vertical grid.
#
# Run after fetch_merra2_native.sh.

function main()
    dir = merra2_native_dir()
    isdir(dir) || error("No granule directory at $dir. Run fetch_merra2_native.sh first.")
    files = sort(filter(f -> endswith(f, ".nc4"), readdir(dir)))
    isempty(files) && error("No .nc4 granules in $dir. Run fetch_merra2_native.sh first.")
    @printf("granules in %s: %d\n", dir, length(files))

    case = first(default_study_cases())
    pair = build_trajectory_pair(case)
    t0 = first(pair.aerocapture).dt
    probe = pair.edl
    lats = [p.lat_deg for p in probe]; lons = [p.lon_deg for p in probe]

    window = load_merra2_native_span(
        first(probe).dt, last(probe).dt;
        dir, lat_min=minimum(lats) - 6, lat_max=maximum(lats) + 6,
        lon_min=minimum(lons) - 6, lon_max=maximum(lons) + 6,
    )
    @printf("\nwindow: %d levels x %d lat x %d lon x %d times, %s to %s\n",
            size(window.density, 1), length(window.lats), length(window.lons),
            length(window.times), first(window.times), last(window.times))

    lat, lon = 20.0, -65.0
    ceil_native = merra2_native_ceiling_m(window, lat, lon)
    ceil_clim = merra2_ceiling_m(load_merra2_grid(Dates.month(t0), 9), lat, lon)
    println("\n--- 1. Ceiling: native should reach ~80 km against the climatology's ~64 km ---")
    @printf("  native      %.2f km\n  climatology %.2f km\n  gain        %.2f km\n",
            ceil_native * 1e-3, ceil_clim * 1e-3, (ceil_native - ceil_clim) * 1e-3)
    ceil_native > ceil_clim + 8.0e3 ||
        @warn "Native ceiling is not meaningfully above the climatology's; check the level ordering."

    # Probe inside the loaded span: the window covers only the analyses bracketing
    # this case's corridors, so a fixed offset can fall outside it.
    probe_dt = first(window.times) + Millisecond(Dates.value(last(window.times) - first(window.times)) ÷ 2)
    @printf("  probing at %s\n", probe_dt)

    println("\n--- 2. Monotonicity: density must fall with altitude ---")
    zs = 5.0:5.0:75.0
    prof = [merra2_native_density(window, lat, lon, 1.0e3 * z, probe_dt) for z in zs]
    ok = filter(isfinite, prof)
    @printf("  %d/%d finite, monotone decreasing: %s\n", length(ok), length(prof),
            length(ok) >= 2 ? string(all(diff(ok) .< 0)) : "too few finite samples")
    length(ok) >= 2 || @warn "Almost nothing was finite — check the probe time is inside the loaded window."
    length(ok) >= 2 && !all(diff(ok) .< 0) &&
        @warn "Density is not monotone in altitude — the vertical grid is likely flipped."

    println("\n--- 3. Height convention: native vs GRAM nominal over 10-55 km ---")
    println("    GRAM is queried at geometric altitude. Whichever column is closer to")
    println("    zero is the correct reading of MERRA-2's H. Residual weather makes a")
    println("    few percent normal; a systematic ~1.5% ramp with altitude is the tell.")
    model = _build_gram_model(_to_initial_time(t0), "earth", 1)
    println("   alt_km   log(m2/gram) as-read    same, H taken as geometric")
    for z in 10.0:5.0:55.0
        p = TrajectoryPoint(probe_dt, 0.0, lat, lon, 1.0e3 * z)
        rg = _gram_density_at_point(model, p; elapsed_time_s=_elapsed_from_initial_s(p, t0))
        r_conv = merra2_native_density(window, lat, lon, 1.0e3 * z, p.dt)
        # Undo the conversion by asking for the altitude that the raw H would
        # have labelled this level.
        z_raw = EARTH_RADIUS_M * (1.0e3 * z) / (EARTH_RADIUS_M + 1.0e3 * z)
        r_raw = merra2_native_density(window, lat, lon, z_raw, p.dt)
        (isfinite(rg) && isfinite(r_conv) && isfinite(r_raw)) || continue
        @printf("   %5.1f   %+10.4f              %+10.4f\n", z, log(r_conv / rg), log(r_raw / rg))
    end

    println("\n--- 4. Time evolution: the gap must actually move the atmosphere ---")
    println("    This is what the climatology could not provide.")
    span_h = Dates.value(last(window.times) - first(window.times)) / 3.6e6
    for gap_h in (1, 3, 6)
        gap_h > span_h && continue
        a = merra2_native_density(window, lat, lon, 35.0e3, first(window.times))
        b = merra2_native_density(window, lat, lon, 35.0e3, first(window.times) + Hour(gap_h))
        (isfinite(a) && isfinite(b)) || continue
        @printf("   %d hr: log-density change %+.5f\n", gap_h, log(b / a))
    end
    @printf("   (loaded window spans %.1f hr; longer gaps need a case with a longer gap_s)\n", span_h)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
