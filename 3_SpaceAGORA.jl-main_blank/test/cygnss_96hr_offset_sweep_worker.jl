# cygnss_96hr_offset_sweep_worker.jl
#
# Runs the CYGNSS 96hr scenario at a batch of initial-time offsets and records
# the RMSE for each.  Designed to be invoked in parallel across multiple Julia
# processes, each handling a different slice of the sweep.
#
# Required env vars:
#   SPACEAGORA_96HR_OFFSETS  – comma-separated integer second offsets, e.g. "-5,-4,-3"
#   SPACEAGORA_96HR_OUTPUT   – path for the output CSV (columns defined below)

using Arrow
using DataFrames
using TOML
using StaticArrays
using LinearAlgebra
using Statistics

const _REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

if !isdefined(Main, :SimulationModel)
    include(joinpath(_REPO_ROOT, "src", "core", "simulation_model.jl"))
end
if !isdefined(Main, :SimulationEngine)
    include(joinpath(_REPO_ROOT, "src", "simulation", "engine", "simulation_engine.jl"))
end
if !isdefined(Main, :TelemetryVerification)
    include(joinpath(_REPO_ROOT, "src", "analysis", "verification", "telemetry_verification.jl"))
end

const TV = TelemetryVerification

const _96HR_FEATHER   = joinpath(_REPO_ROOT, "data", "telemetry", "CYGNSS", "cyg01_nasa_pvt_96hr.feather")
const _HARMONICS_FILE = "data/Gravity_harmonics_data/EarthGGM05C.csv"

# ---------------------------------------------------------------------------
# Convert an integer second offset from the base epoch (2025-06-05 00:00:00 UTC)
# into a scenario initial_time dict, handling minute and hour rollover.
# ---------------------------------------------------------------------------
function _offset_to_epoch(offset_s::Int)::Dict{String, Any}
    year, month, day = 2025, 6, 5
    hour, minute = 0, 0
    second = Float64(offset_s)

    while second < 0.0;   second += 60.0; minute -= 1; end
    while second >= 60.0; second -= 60.0; minute += 1; end
    while minute < 0;     minute += 60;   hour   -= 1; end
    while minute >= 60;   minute -= 60;   hour   += 1; end
    while hour < 0;       hour   += 24;   day    -= 1; end
    while hour >= 24;     hour   -= 24;   day    += 1; end

    return Dict{String, Any}(
        "year" => year, "month" => month, "day" => day,
        "hour" => hour, "minute" => minute, "second" => second
    )
end

# ---------------------------------------------------------------------------
# Build the standardized Arrow telemetry file from the raw 96hr feather.
# This is done once per process and shared across all offset runs.
# ---------------------------------------------------------------------------
function _build_reference(outdir::String)::String
    df = DataFrame(Arrow.Table(_96HR_FEATHER))

    t_unix  = Float64.(df[!, "pvt_unix_seconds"])
    x_km    = Float64.(df[!, "sc_pos_x_pvt_m"])  .* 1e-3
    y_km    = Float64.(df[!, "sc_pos_y_pvt_m"])  .* 1e-3
    z_km    = Float64.(df[!, "sc_pos_z_pvt_m"])  .* 1e-3
    vx_kmps = Float64.(df[!, "sc_vel_x_pvt_mps"]) .* 1e-3
    vy_kmps = Float64.(df[!, "sc_vel_y_pvt_mps"]) .* 1e-3
    vz_kmps = Float64.(df[!, "sc_vel_z_pvt_mps"]) .* 1e-3

    perm = sortperm(t_unix)
    t_unix  = t_unix[perm]
    x_km    = x_km[perm];  y_km    = y_km[perm];  z_km    = z_km[perm]
    vx_kmps = vx_kmps[perm]; vy_kmps = vy_kmps[perm]; vz_kmps = vz_kmps[perm]

    t_rel  = t_unix .- t_unix[1]
    planet = TV._planet_from_name("earth")
    r_km   = sqrt.(x_km .^ 2 .+ y_km .^ 2 .+ z_km .^ 2)
    alt_km = r_km .- planet.Rp_e * 1e-3

    oe0 = TV.rvtoorbitalelement(
        SVector{3, Float64}(x_km[1], y_km[1], z_km[1]) .* 1e3,
        SVector{3, Float64}(vx_kmps[1], vy_kmps[1], vz_kmps[1]) .* 1e3,
        planet
    )

    n = length(t_rel)
    tdf = DataFrame(
        time_s      = t_rel,
        altitude_km = alt_km,
        x_km = x_km, y_km = y_km, z_km = z_km,
        sma_km   = fill(oe0[1] * 1e-3, n),
        ecc      = fill(oe0[2], n),
        inc_deg  = fill(rad2deg(oe0[3]), n),
        raan_deg = fill(rad2deg(oe0[4]), n),
        aop_deg  = fill(rad2deg(oe0[5]), n),
        ta_deg   = fill(rad2deg(oe0[6]), n),
        x_ic_km    = fill(x_km[1],    n),
        y_ic_km    = fill(y_km[1],    n),
        z_ic_km    = fill(z_km[1],    n),
        vx_ic_kmps = fill(vx_kmps[1], n),
        vy_ic_kmps = fill(vy_kmps[1], n),
        vz_ic_kmps = fill(vz_kmps[1], n)
    )

    path = joinpath(outdir, "cygnss_96hr_ref.arrow")
    Arrow.write(path, tdf)
    return path
end

# ---------------------------------------------------------------------------
# Run verification for one time offset; returns a NamedTuple of metrics.
# ---------------------------------------------------------------------------
function _run_offset(offset_s::Int, telemetry_path::String)
    initial_time = _offset_to_epoch(offset_s)

    scenario = Dict{String, Any}(
        "name"   => "cygnss_96hr_pvt",
        "kind"   => "time_aligned_state",
        "events" => Any["altitude_time", "state_x_time", "state_y_time", "state_z_time"],
        "telemetry" => telemetry_path,
        "telemetry_columns" => Dict{String, Any}(
            "time" => "time_s", "altitude" => "altitude_km",
            "x"    => "x_km",  "y"         => "y_km",  "z"  => "z_km",
            "x_ic" => "x_ic_km", "y_ic"    => "y_ic_km", "z_ic" => "z_ic_km",
            "vx_ic" => "vx_ic_kmps", "vy_ic" => "vy_ic_kmps", "vz_ic" => "vz_ic_kmps"
        ),
        "max_points_quick" => 500000,
        "max_points_full"  => 500000,
        "min_eval_points"  => 2,
        "units" => Dict{String, Any}(
            "x"            => "s",  "altitude_time" => "km",
            "state_x_time" => "km", "state_y_time"  => "km", "state_z_time" => "km"
        ),
        "tolerances_quick" => Dict{String, Any}(
            "altitude_time" => Dict("max_abs_km"=>1e6, "max_nmae"=>1e6, "max_rmse_km"=>1e6),
            "state_x_time"  => Dict("max_abs_km"=>1e6, "max_nmae"=>1e6, "max_rmse_km"=>1e6),
            "state_y_time"  => Dict("max_abs_km"=>1e6, "max_nmae"=>1e6, "max_rmse_km"=>1e6),
            "state_z_time"  => Dict("max_abs_km"=>1e6, "max_nmae"=>1e6, "max_rmse_km"=>1e6)
        ),
        "tolerances_full" => Dict{String, Any}(
            "altitude_time" => Dict("max_abs_km"=>1e6, "max_nmae"=>1e6, "max_rmse_km"=>1e6),
            "state_x_time"  => Dict("max_abs_km"=>1e6, "max_nmae"=>1e6, "max_rmse_km"=>1e6),
            "state_y_time"  => Dict("max_abs_km"=>1e6, "max_nmae"=>1e6, "max_rmse_km"=>1e6),
            "state_z_time"  => Dict("max_abs_km"=>1e6, "max_nmae"=>1e6, "max_rmse_km"=>1e6)
        ),
        "planet"                   => "earth",
        "gravity_model"            => "inverse_squared",
        "gravity_harmonics_degree" => 50,
        "gravity_harmonics_order"  => 50,
        "gravity_harmonics_file"   => _HARMONICS_FILE,
        "nbody_bodies"             => Any[],
        "orbit_altitude_mode"      => "oblate",
        "drag_enabled"             => false,
        "include_wind"             => false,
        "cartesian_ic_frame"       => "planet_fixed",
        "comparison_frame"         => "planet_fixed",
        "initial_time"             => initial_time,
        "EI_km"                    => 300.0,
        "comparison_mode"          => "time_aligned_state",
        "spacecraft" => Dict{String, Any}(
            "bus_dims_m"         => Any[2.05e-1, 3.7e-1, 0.8e-1],
            "panel_dims_m"       => Any[10e-3, 28.5e-3, 0.0001],
            "bus_mass_kg"        => 29.0,
            "panel_mass_each_kg" => 0.0,
            "panel_offset_y_m"   => 2.45,
            "prop_mass_kg"       => 0.0,
            "id"                 => 1002
        ),
        "atmosphere_truth" => Dict{String, Any}(
            "assumption_id"           => "earth_gmat_gram_deterministic_v1",
            "atmosphere_model"        => "GRAM",
            "atmosphere_dataset"      => "MERRA2 All Mean",
            "space_weather_model"     => "EarthGRAM MERRA2 climatology (deterministic)",
            "solar_flux_model"        => "EarthGRAM/MERRA2 epoch-fixed defaults",
            "gram_seed"               => 1001,
            "gram_perturbation_scales"=> Any[0.0, 0.0, 0.0, 0.0],
            "gram_offline_surrogate"  => "auto",
            "gram_static_grid"        => false,
            "gram_track_cache"        => false,
            "gram_global_lock"        => "on"
        ),
        "calibration" => Dict{String, Any}("enabled" => false)
    )

    mktempdir() do tmp
        manifest_path = joinpath(tmp, "manifest.toml")
        open(manifest_path, "w") do io
            TOML.print(io, Dict{String, Any}("scenarios" => Any[scenario]))
        end

        req = TV.VerificationRequest(
            profile        = :quick,
            out_summary    = joinpath(tmp, "summary.csv"),
            out_errors     = joinpath(tmp, "errors.csv"),
            manifest_path  = manifest_path,
            enforce        = false,
            generate_plots = false
        )

        result = withenv(
            "SPACEAGORA_TELEMETRY_SOLVER_MODE"   => "dp8",
            "SPACEAGORA_TELEMETRY_DT_MAX_ORBIT"  => "10.0",
            "SPACEAGORA_TELEMETRY_RELTOL_ORBIT"  => "1e-12",
            "SPACEAGORA_TELEMETRY_ABSTOL_ORBIT"  => "1e-12",
            "SPACEAGORA_TELEMETRY_RELTOL_ATM"    => "1e-12",
            "SPACEAGORA_TELEMETRY_ABSTOL_ATM"    => "1e-12"
        ) do
            TV.run_verification(req)
        end

        summary = result.summary
        retcode = let rows = summary[summary.scenario .== "cygnss_96hr_pvt", :]
            nrow(rows) > 0 ? String(rows.solver_retcode[1]) : "Unknown"
        end

        get_rmse(event) = begin
            rows = summary[summary.event .== event, :]
            nrow(rows) > 0 ? Float64(rows.rmse_km[1]) : NaN
        end

        pos_rows = summary[
            (summary.scenario .== "cygnss_96hr_pvt") .&
            in.(summary.event, Ref(["state_x_time", "state_y_time", "state_z_time"])), :]
        mean_pos_rmse = nrow(pos_rows) == 3 ? mean(Float64.(pos_rows.rmse_km)) : NaN

        return (
            offset_s          = offset_s,
            rmse_altitude_km  = get_rmse("altitude_time"),
            rmse_x_km         = get_rmse("state_x_time"),
            rmse_y_km         = get_rmse("state_y_time"),
            rmse_z_km         = get_rmse("state_z_time"),
            mean_pos_rmse_km  = mean_pos_rmse,
            solver_retcode    = retcode
        )
    end
end

# ---------------------------------------------------------------------------
# Main: parse env vars, build reference once, loop over offsets.
# ---------------------------------------------------------------------------
offsets_str = get(ENV, "SPACEAGORA_96HR_OFFSETS", "0")
output_path = get(ENV, "SPACEAGORA_96HR_OUTPUT", "/tmp/cygnss_96hr_sweep_result.csv")

offsets = [parse(Int, strip(s)) for s in split(offsets_str, ",") if !isempty(strip(s))]
isempty(offsets) && error("SPACEAGORA_96HR_OFFSETS is empty or missing")

println("Worker starting: offsets = $offsets → $output_path")
flush(stdout)

results = mktempdir() do tmp
    telemetry_path = _build_reference(tmp)
    println("Reference telemetry built: $telemetry_path")
    flush(stdout)

    map(offsets) do offset_s
        println("Running offset = $(offset_s)s ...")
        flush(stdout)
        r = _run_offset(offset_s, telemetry_path)
        println("  offset=$(r.offset_s)s  mean_pos_rmse=$(round(r.mean_pos_rmse_km, digits=2)) km  retcode=$(r.solver_retcode)")
        flush(stdout)
        r
    end
end

# Write CSV
mkpath(dirname(output_path))
open(output_path, "w") do io
    println(io, "offset_s,rmse_altitude_km,rmse_x_km,rmse_y_km,rmse_z_km,mean_pos_rmse_km,solver_retcode")
    for r in results
        println(io, "$(r.offset_s),$(r.rmse_altitude_km),$(r.rmse_x_km),$(r.rmse_y_km),$(r.rmse_z_km),$(r.mean_pos_rmse_km),$(r.solver_retcode)")
    end
end

println("Done. Results written to $output_path")
