# Direct density comparison for Figure 1 of the aerobraking reconstruction
# record: Odyssey accelerometer-derived (Tolson) density divided by the
# TES Mapping Year 2 MarsGRAM density at matched time and location.
#
# For every accelerometer altitude bin inside [alt-min, alt-max] (record: 110
# to 130 km), the matching point on the flight trajectory is located from the
# Mars Odyssey navigation kernel m01_ab_v2.bsp (the altitude crossing of the
# inbound or outbound leg around the tagged periapsis), and MarsGRAM is
# queried at that epoch, latitude, longitude, and altitude through the same
# point-density path the simulator uses. The record's comparison spans 2,938
# points over 315 passes.
#
# Interpretation caveat (record §6): this is NOT a direct-density measurement
# comparison — the Tolson density embeds a flight-side aerodynamic model.
# The panel exists to interpret the MarsGRAM comparator run.
#
# Requires the native GRAM Suite build. Usage:
#   julia --project=. benchmarks/studies/telemetry_validation/density_compare.jl \
#       [--alt-min=110] [--alt-max=130] [--map-year=2] [--window-s=1200] \
#       [--step-s=2.0] [--max-passes=0] [--out-root=DIR]

include(joinpath(@__DIR__, "common.jl"))

const DC_OPTS = parse_kv_args(copy(ARGS), ())
load_gramsuite!()

using SpaceAGORA
using SpaceAGORA.SimulationModel
using SpaceAGORA.SimulationModel: InitialTime
using SpaceAGORA.TelemetryVerification
using Arrow
using CSV
using DataFrames
using Printf
using StaticArrays
using Statistics
using SPICE: furnsh, str2et, spkgeo

const TV = SpaceAGORA.TelemetryVerification
const EM = SpaceAGORA.SimulationModel.EnvironmentModels

const DENSITY_FILE = joinpath(REPO_ROOT, "data", "telemetry", "Odyssey", "odyssey_accelerometer_density.feather")
const MISSION_KERNEL = joinpath(TV.SPICE_PATH, "spk", "missions", "m01_ab_v2.bsp")
const ODYSSEY_NAIF_ID = -53
const MARS_NAIF_ID = 499
# Scenario epoch of the benchmark of record; GRAM elapsed time is measured
# from here (pre-epoch walk-in passes get negative elapsed times).
const SCENARIO_EPOCH = InitialTime(year=2001, month=11, day=6, hour=19, minute=0, second=32.0)

function kernel_state_lla(et::Float64, planet)
    state, _ = spkgeo(ODYSSEY_NAIF_ID, et, "IAU_MARS", MARS_NAIF_ID)
    pos_m = SVector{3, Float64}(state[1], state[2], state[3]) .* 1e3
    return TV.rtolatlong(pos_m, planet)
end

# Locate the epoch at which the kernel trajectory crosses `target_alt_m` on
# the requested leg around the pass periapsis, by bisection on a bracketing
# sample pair. Returns nothing when the window does not cover the crossing.
function find_crossing(
    et_grid::Vector{Float64},
    alt_grid::Vector{Float64},
    idx_peri::Int,
    target_alt_m::Float64,
    leg::AbstractString,
    planet
)
    rng = leg == "in" ? (1:(idx_peri - 1)) : (idx_peri:(length(et_grid) - 1))
    for i in rng
        a, b = alt_grid[i] - target_alt_m, alt_grid[i + 1] - target_alt_m
        a * b > 0.0 && continue
        # Inbound legs must be descending, outbound ascending, so a pass with
        # multiple crossings (long window) keeps the physically correct one.
        descending = alt_grid[i + 1] < alt_grid[i]
        (leg == "in") == descending || continue
        lo, hi = et_grid[i], et_grid[i + 1]
        for _ in 1:40
            mid = 0.5 * (lo + hi)
            alt_mid = kernel_state_lla(mid, planet)[1]
            if (alt_mid - target_alt_m) * a <= 0.0
                hi = mid
            else
                lo = mid
                a = alt_mid - target_alt_m
            end
        end
        return 0.5 * (lo + hi)
    end
    return nothing
end

function main()
    alt_min_km = parse(Float64, get(DC_OPTS, "alt-min", "110.0"))
    alt_max_km = parse(Float64, get(DC_OPTS, "alt-max", "130.0"))
    map_year = parse(Int, get(DC_OPTS, "map-year", "2"))
    window_s = parse(Float64, get(DC_OPTS, "window-s", "1200.0"))
    step_s = parse(Float64, get(DC_OPTS, "step-s", "2.0"))
    max_passes = parse(Int, get(DC_OPTS, "max-passes", "0"))
    out_dir = joinpath(results_root(get(DC_OPTS, "out-root", "")), "density_compare")
    mkpath(out_dir)
    out_csv = joinpath(out_dir, "odyssey_accel_vs_marsgram_density.csv")

    isfile(DENSITY_FILE) || error("Missing accelerometer density archive: $DENSITY_FILE")
    isfile(MISSION_KERNEL) || error("Missing Odyssey navigation kernel: $MISSION_KERNEL")

    # Mars(...) furnishes the leapsecond/PCK/planetary kernels; the mission
    # kernel is furnished on top of them.
    planet = SimulationModel.Mars("", TV.SPICE_PATH)
    furnsh(MISSION_KERNEL)
    et0 = TV._initial_time_et(SCENARIO_EPOCH)

    truth = TV.AtmosphereTruthConfig(
        assumption_id="odyssey_marsgram_tes_my2_density_compare",
        atmosphere_model="GRAM",
        atmosphere_dataset="MarsGRAM TES Mapping Year 2 (MY25) climatology",
        gram_seed=1001,
        gram_perturbation_scales=(0.0, 0.0, 0.0, 0.0),
        gram_offline_surrogate="off",
        mars_map_year=map_year
    )
    model = TV._make_required_gram_density_model("mars", SCENARIO_EPOCH, truth)

    density = DataFrame(Arrow.Table(DENSITY_FILE))
    in_band = density[(density.alt_km .>= alt_min_km) .& (density.alt_km .<= alt_max_km), :]
    passes = sort(unique(in_band.P))
    max_passes > 0 && (passes = passes[1:min(end, max_passes)])
    n_selected_rows = count(p -> p in passes, in_band.P)
    println("Accelerometer archive: $(nrow(density)) rows, $(length(unique(density.P))) passes")
    println("Comparison band $(alt_min_km)-$(alt_max_km) km: $n_selected_rows candidate rows over $(length(passes)) selected passes")

    rows = NamedTuple[]
    skipped_no_crossing = 0
    skipped_pass = 0
    for pass_number in passes
        pass_rows = in_band[in_band.P .== pass_number, :]
        et_peri = try
            str2et(String(pass_rows.t_peri_utc[1]))
        catch err
            @warn "Skipping pass $pass_number: bad periapsis tag" error=sprint(showerror, err)
            skipped_pass += 1
            continue
        end
        et_grid = collect((et_peri - window_s):step_s:(et_peri + window_s))
        alt_grid = try
            [kernel_state_lla(et, planet)[1] for et in et_grid]
        catch err
            @warn "Skipping pass $pass_number: kernel coverage" error=sprint(showerror, err)
            skipped_pass += 1
            continue
        end
        idx_peri = argmin(alt_grid)
        for row in eachrow(pass_rows)
            et_x = find_crossing(et_grid, alt_grid, idx_peri, row.alt_km * 1e3, String(row.leg), planet)
            if et_x === nothing
                skipped_no_crossing += 1
                continue
            end
            lla = kernel_state_lla(et_x, planet)
            rho_gram, _, _ = Base.invokelatest(
                EM._gram_point_density, model, lla[1], lla[2], lla[3], et_x - et0, true
            )
            rho_gram > 0.0 || (skipped_no_crossing += 1; continue)
            push!(rows, (
                P=row.P,
                leg=String(row.leg),
                alt_km=row.alt_km,
                elapsed_s=et_x - et0,
                lat_deg=rad2deg(lla[2]),
                lon_deg=rad2deg(lla[3]),
                rho_accel_kgm3=row.rho_kgm3,
                rho_gram_kgm3=rho_gram,
                density_ratio=row.rho_kgm3 / rho_gram
            ))
        end
    end

    isempty(rows) && error("No matched density points produced; check GRAM build and kernel coverage.")
    out = DataFrame(rows)
    CSV.write(out_csv, out)

    ratios = out.density_ratio
    println(@sprintf(
        "Matched %d density points over %d passes (skipped: %d rows without a crossing, %d passes).",
        nrow(out), length(unique(out.P)), skipped_no_crossing, skipped_pass
    ))
    println(@sprintf(
        "density ratio accel/MarsGRAM: median %.3f  mean %.3f  p5 %.3f  p95 %.3f",
        median(ratios), mean(ratios), quantile(ratios, 0.05), quantile(ratios, 0.95)
    ))
    println("Saved: $out_csv")
    return nothing
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
