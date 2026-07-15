##
# Position-accuracy replication of test/gmat_scenario_matrix.jl's
# "CYGNSS Legacy 48hr Entry Point" testset (the ~1.58 km axis-averaged RMSE
# result quoted throughout this study's README), using the SAME 48-hour
# dedicated PVT telemetry file the rest of this folder's 1-hour slew study
# does NOT use.
#
# This does not hand-reimplement the verification pipeline (an earlier
# version of this script did, and after matching every setting visible in
# the scenario dict -- data, N=5 SMA-averaged IC, 50x50 EGM05C harmonics,
# epoch, dp8 solver, no drag -- it still landed at ~4.5 km instead of
# ~1.58 km with no identified cause). Instead it copies the handful of
# helper functions test/gmat_scenario_matrix.jl uses to build and run this
# exact scenario (_build_cygnss_48hr_reference, _base_scenario_dict,
# _scenario_rmse, _telemetry_solver_env_overrides) verbatim and calls the
# real TV.VerificationRequest/TV.run_verification pipeline directly, so
# there is no hand-reimplementation left to diverge from the original.
#
# Usage: julia --project=. benchmarks/studies/cygnss_attitude_torque_reconstruction/run_position_accuracy_48hr.jl
##

using Test
using TOML
using Arrow
using CSV
using DataFrames
using Statistics
using LinearAlgebra
using StaticArrays

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :SimulationModel)
    include(joinpath(REPO_ROOT, "src", "core", "simulation_model.jl"))
end
if !isdefined(Main, :SimulationEngine)
    include(joinpath(REPO_ROOT, "src", "simulation", "engine", "simulation_engine.jl"))
end
if !isdefined(Main, :TelemetryVerification)
    include(joinpath(REPO_ROOT, "src", "analysis", "verification", "telemetry_verification.jl"))
end

const TV = TelemetryVerification
const SM = SimulationModel

const STUDY_DIR = @__DIR__
const PLOTS_DIR = joinpath(STUDY_DIR, "plots")
mkpath(PLOTS_DIR)

const _GMAT_HARMONICS_EARTH_FILE = joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv")
const _CYGNSS_48HR_TELEMETRY_FEATHER = joinpath(REPO_ROOT, "data", "telemetry", "CYGNSS", "cygnss_data_48hr.feather")
const TEST_MODE::Symbol = :quick

# ==============================================================================
# Copied verbatim from test/gmat_scenario_matrix.jl (see file header).
# ==============================================================================

@inline function _parse_bool_env(name::String, default::Bool)::Bool
    token = lowercase(strip(get(ENV, name, default ? "1" : "0")))
    if token in ("1", "true", "yes", "on")
        return true
    elseif token in ("0", "false", "no", "off")
        return false
    end
    return default
end

@inline function _gmat_parity_solver_mode_enabled()::Bool
    return _parse_bool_env("SPACEAGORA_GMAT_PARITY_SOLVER", false)
end

function _telemetry_solver_env_overrides()::Dict{String, String}
    if _gmat_parity_solver_mode_enabled()
        return Dict(
            "SPACEAGORA_TELEMETRY_SOLVER_MODE" => "auto_stiff",
            "SPACEAGORA_TELEMETRY_DT_MAX_ORBIT" => "1.0",
            "SPACEAGORA_TELEMETRY_RELTOL_ORBIT" => "1e-13",
            "SPACEAGORA_TELEMETRY_ABSTOL_ORBIT" => "1e-13",
            "SPACEAGORA_TELEMETRY_RELTOL_ATM" => "1e-13",
            "SPACEAGORA_TELEMETRY_ABSTOL_ATM" => "1e-13"
        )
    end
    return Dict(
        "SPACEAGORA_TELEMETRY_SOLVER_MODE" => "dp8",
        "SPACEAGORA_TELEMETRY_DT_MAX_ORBIT" => "10.0",
        "SPACEAGORA_TELEMETRY_RELTOL_ORBIT" => "1e-12",
        "SPACEAGORA_TELEMETRY_ABSTOL_ORBIT" => "1e-12",
        "SPACEAGORA_TELEMETRY_RELTOL_ATM" => "1e-12",
        "SPACEAGORA_TELEMETRY_ABSTOL_ATM" => "1e-12"
    )
end

function _required_column(df::DataFrame, candidates::Vector{String})::String
    for cname in candidates
        if cname in names(df)
            return cname
        end
    end
    throw(ArgumentError("None of the candidate columns were found: $(join(candidates, ", "))"))
end

function _base_scenario_dict(name::String, telemetry_path::String)
    return Dict{String, Any}(
        "name" => name,
        "kind" => "time_aligned_state",
        "events" => Any["altitude_time", "state_x_time", "state_y_time", "state_z_time"],
        "telemetry" => telemetry_path,
        "telemetry_columns" => Dict{String, Any}(
            "time" => "time_s",
            "altitude" => "altitude_km",
            "x" => "x_km",
            "y" => "y_km",
            "z" => "z_km",
            "x_ic" => "x_ic_km",
            "y_ic" => "y_ic_km",
            "z_ic" => "z_ic_km",
            "vx_ic" => "vx_ic_kmps",
            "vy_ic" => "vy_ic_kmps",
            "vz_ic" => "vz_ic_kmps"
        ),
        "max_points_quick" => 10000,
        "max_points_full" => 100000,
        "min_eval_points" => 2,
        "units" => Dict{String, Any}(
            "x" => "s",
            "altitude_time" => "km",
            "state_x_time" => "km",
            "state_y_time" => "km",
            "state_z_time" => "km"
        ),
        "tolerances_quick" => Dict{String, Any}(
            "altitude_time" => Dict("max_abs_km" => 1.0e6, "max_nmae" => 1.0e6, "max_rmse_km" => 1.0e6),
            "state_x_time" => Dict("max_abs_km" => 1.0e6, "max_nmae" => 1.0e6, "max_rmse_km" => 1.0e6),
            "state_y_time" => Dict("max_abs_km" => 1.0e6, "max_nmae" => 1.0e6, "max_rmse_km" => 1.0e6),
            "state_z_time" => Dict("max_abs_km" => 1.0e6, "max_nmae" => 1.0e6, "max_rmse_km" => 1.0e6)
        ),
        "tolerances_full" => Dict{String, Any}(
            "altitude_time" => Dict("max_abs_km" => 1.0e6, "max_nmae" => 1.0e6, "max_rmse_km" => 1.0e6),
            "state_x_time" => Dict("max_abs_km" => 1.0e6, "max_nmae" => 1.0e6, "max_rmse_km" => 1.0e6),
            "state_y_time" => Dict("max_abs_km" => 1.0e6, "max_nmae" => 1.0e6, "max_rmse_km" => 1.0e6),
            "state_z_time" => Dict("max_abs_km" => 1.0e6, "max_nmae" => 1.0e6, "max_rmse_km" => 1.0e6)
        ),
        "initial_time" => Dict{String, Any}(
            "year" => 2026,
            "month" => 1,
            "day" => 1,
            "hour" => 12,
            "minute" => 0,
            "second" => 0.0
        ),
        "spacecraft" => Dict{String, Any}(
            "bus_dims_m" => Any[2.05e-1, 3.7e-1, 0.8e-1],
            "panel_dims_m" => Any[10e-3, 28.5e-3, 0.0001],
            "bus_mass_kg" => 29.0,
            "panel_mass_each_kg" => 0.0,
            "panel_offset_y_m" => 2.45,
            "prop_mass_kg" => 0.0,
            "id" => 1002
        ),
        "atmosphere_truth" => Dict{String, Any}(
            "assumption_id" => "earth_gmat_gram_deterministic_v1",
            "atmosphere_model" => "GRAM",
            "atmosphere_dataset" => "MERRA2 All Mean",
            "space_weather_model" => "EarthGRAM MERRA2 climatology (deterministic)",
            "solar_flux_model" => "EarthGRAM/MERRA2 epoch-fixed defaults",
            "gram_seed" => 1001,
            "gram_perturbation_scales" => Any[0.0, 0.0, 0.0, 0.0],
            "gram_offline_surrogate" => "auto",
            "gram_static_grid" => false,
            "gram_track_cache" => false,
            "gram_global_lock" => "on"
        ),
        "calibration" => Dict{String, Any}("enabled" => false),
        "drag_enabled" => false,
        "EI_km" => 300.0,
        "comparison_mode" => "time_aligned_state"
    )
end

function _scenario_rmse(summary::DataFrame, name::String)::Float64
    rows = summary[(summary.scenario .== name) .& in.(summary.event, Ref(["state_x_time", "state_y_time", "state_z_time"])), :]
    @assert nrow(rows) == 3
    rmse = mean(Float64.(rows.rmse_km))
    @assert isfinite(rmse)
    return rmse
end

function _build_cygnss_48hr_reference(outdir::String, stem::String)
    df = DataFrame(Arrow.Table(_CYGNSS_48HR_TELEMETRY_FEATHER))
    @assert nrow(df) > 0

    time_col = _required_column(df, ["TIME OFFSET", "time"])
    x_col = _required_column(df, ["OBS4.ENG_PVT.DDMI_PVT_SCPOS_X (m)", "pos_ii_1"])
    y_col = _required_column(df, ["OBS4.ENG_PVT.DDMI_PVT_SCPOS_Y (m)", "pos_ii_2"])
    z_col = _required_column(df, ["OBS4.ENG_PVT.DDMI_PVT_SCPOS_Z (m)", "pos_ii_3"])
    vx_col = _required_column(df, ["OBS4.ENG_PVT.DDMI_PVT_SCVEL_X (m/s)", "vel_ii_1"])
    vy_col = _required_column(df, ["OBS4.ENG_PVT.DDMI_PVT_SCVEL_Y (m/s)", "vel_ii_2"])
    vz_col = _required_column(df, ["OBS4.ENG_PVT.DDMI_PVT_SCVEL_Z (m/s)", "vel_ii_3"])

    t = Float64.(df[!, time_col])
    x_km = Float64.(df[!, x_col]) .* 1.0e-3
    y_km = Float64.(df[!, y_col]) .* 1.0e-3
    z_km = Float64.(df[!, z_col]) .* 1.0e-3
    vx_kmps = Float64.(df[!, vx_col]) .* 1.0e-3
    vy_kmps = Float64.(df[!, vy_col]) .* 1.0e-3
    vz_kmps = Float64.(df[!, vz_col]) .* 1.0e-3

    perm = sortperm(t)
    t = t[perm]
    x_km = x_km[perm]; y_km = y_km[perm]; z_km = z_km[perm]
    vx_kmps = vx_kmps[perm]; vy_kmps = vy_kmps[perm]; vz_kmps = vz_kmps[perm]

    t_rel = t .- t[1]
    planet = TV._planet_from_name("earth")
    r_km = sqrt.(x_km .^ 2 .+ y_km .^ 2 .+ z_km .^ 2)
    alt_km = r_km .- planet.Rp_e * 1.0e-3

    _n_sma_avg = min(5, length(x_km))
    _a_samples = [
        1.0 / (2.0 / sqrt(x_km[k]^2 + y_km[k]^2 + z_km[k]^2) - (vx_kmps[k]^2 + vy_kmps[k]^2 + vz_kmps[k]^2) / (planet.μ * 1.0e-9))
        for k in 1:_n_sma_avg
    ]
    _a_target_km = mean(_a_samples)
    _r0_km = sqrt(x_km[1]^2 + y_km[1]^2 + z_km[1]^2)
    _v0_kmps = sqrt(vx_kmps[1]^2 + vy_kmps[1]^2 + vz_kmps[1]^2)
    _v_target_kmps = sqrt((planet.μ * 1.0e-9) * (2.0 / _r0_km - 1.0 / _a_target_km))
    _v_scale = _v_target_kmps / _v0_kmps

    x_ic_km, y_ic_km, z_ic_km = x_km[1], y_km[1], z_km[1]
    vx_ic_kmps, vy_ic_kmps, vz_ic_kmps = vx_kmps[1] * _v_scale, vy_kmps[1] * _v_scale, vz_kmps[1] * _v_scale

    oe0 = TV.rvtoorbitalelement(
        SVector{3, Float64}(x_ic_km, y_ic_km, z_ic_km) .* 1.0e3,
        SVector{3, Float64}(vx_ic_kmps, vy_ic_kmps, vz_ic_kmps) .* 1.0e3,
        planet
    )
    sma_km = oe0[1] * 1.0e-3
    ecc = oe0[2]
    inc_deg = rad2deg(oe0[3])
    raan_deg = rad2deg(oe0[4])
    aop_deg = rad2deg(oe0[5])
    ta_deg = rad2deg(oe0[6])

    telemetry_df = DataFrame(
        time_s=t_rel,
        altitude_km=alt_km,
        x_km=x_km, y_km=y_km, z_km=z_km,
        sma_km=fill(sma_km, length(t_rel)),
        ecc=fill(ecc, length(t_rel)),
        inc_deg=fill(inc_deg, length(t_rel)),
        aop_deg=fill(aop_deg, length(t_rel)),
        raan_deg=fill(raan_deg, length(t_rel)),
        ta_deg=fill(ta_deg, length(t_rel)),
        x_ic_km=fill(x_ic_km, length(t_rel)),
        y_ic_km=fill(y_ic_km, length(t_rel)),
        z_ic_km=fill(z_ic_km, length(t_rel)),
        vx_ic_kmps=fill(vx_ic_kmps, length(t_rel)),
        vy_ic_kmps=fill(vy_ic_kmps, length(t_rel)),
        vz_ic_kmps=fill(vz_ic_kmps, length(t_rel))
    )

    telemetry_path = joinpath(outdir, "$(stem)_time_aligned.arrow")
    Arrow.write(telemetry_path, telemetry_df)

    return (
        telemetry_path=telemetry_path, t0_s=t_rel[1], tf_s=t_rel[end],
        x_ic_km=x_ic_km, y_ic_km=y_ic_km, z_ic_km=z_ic_km,
        vx_ic_kmps=vx_ic_kmps, vy_ic_kmps=vy_ic_kmps, vz_ic_kmps=vz_ic_kmps,
        sma_km=sma_km, ecc=ecc, inc_deg=inc_deg, aop_deg=aop_deg, raan_deg=raan_deg, ta_deg=ta_deg
    )
end

# ==============================================================================
# Run the real verification pipeline (test/gmat_scenario_matrix.jl:_run_cygnss_48hr_result_once,
# copied verbatim minus test-framework caching/@test wrapping).
# ==============================================================================

println("=== Running the real TV.VerificationRequest/TV.run_verification pipeline for cygnss_48hr_pvt ===")

result = mktempdir() do tmp
    traj = _build_cygnss_48hr_reference(tmp, "cygnss_48hr_pvt")

    scenario = _base_scenario_dict("cygnss_48hr_pvt", traj.telemetry_path)
    merge!(scenario, Dict{String, Any}(
        "planet" => "earth",
        "gravity_model" => "inverse_squared",
        "gravity_harmonics_degree" => 50,
        "gravity_harmonics_order" => 50,
        "gravity_harmonics_file" => _GMAT_HARMONICS_EARTH_FILE,
        "nbody_bodies" => Any[],
        "orbit_altitude_mode" => "oblate",
        "drag_enabled" => false,
        "include_wind" => false,
        "max_points_quick" => 200000,
        "max_points_full" => 200000,
        "initial_time" => Dict{String, Any}(
            "year" => 2025, "month" => 6, "day" => 6,
            "hour" => 0, "minute" => 0, "second" => 0.0
        )
    ))

    manifest_path = joinpath(tmp, "cygnss_48hr_manifest.toml")
    open(manifest_path, "w") do io
        TOML.print(io, Dict{String, Any}("scenarios" => Any[scenario]))
    end

    req = TV.VerificationRequest(
        profile=TEST_MODE,
        out_summary=joinpath(tmp, "summary.csv"),
        out_errors=joinpath(tmp, "errors.csv"),
        manifest_path=manifest_path,
        enforce=false,
        generate_plots=false
    )

    solver_env = _telemetry_solver_env_overrides()
    return withenv(pairs(solver_env)...) do
        TV.run_verification(req)
    end
end

pos_rmse = _scenario_rmse(result.summary, "cygnss_48hr_pvt")
println()
println("cygnss_48hr_pvt mean position-axis RMSE [km]: $(pos_rmse)")
println("(reference documented value: ~1.58 km)")

CSV.write(joinpath(STUDY_DIR, "data", "cygnss_48hr_reference_summary.csv"), result.summary)
println("full summary written to: $(joinpath(STUDY_DIR, "data", "cygnss_48hr_reference_summary.csv"))")
# result.errors (one row per telemetry sample per event, ~170k rows/event) is
# used below for the plot but deliberately not persisted to disk -- it's
# large (~100MB+) and trivially regenerated by re-running this script.

# ==============================================================================
# Plot: per-axis position error vs time, IQR-filtered to drop isolated spikes
# before plotting (copied verbatim from test/gmat_scenario_matrix.jl:
# _extract_position_error_series / _iqr_inlier_mask / _filter_position_series_
# for_plot). The spikes are real telemetry-comparison samples (not a plotting
# bug), most likely isolated GPS-fix glitches in the 48hr telemetry itself --
# filtering is purely a *plotting* choice; _scenario_rmse above is computed
# from the unfiltered errors table, so the 1.578 km RMSE figure is unaffected
# by this filter.
# ==============================================================================

using Plots

Plots.default(left_margin=10Plots.mm, bottom_margin=8Plots.mm)

function _extract_position_error_series(errors::DataFrame, scenario_name::String)
    rows = errors[errors.scenario .== scenario_name, :]
    xrow = sort(rows[rows.event .== "state_x_time", :], :idx)
    yrow = sort(rows[rows.event .== "state_y_time", :], :idx)
    zrow = sort(rows[rows.event .== "state_z_time", :], :idx)

    n = min(nrow(xrow), nrow(yrow), nrow(zrow))
    n >= 2 || error("Need at least two aligned position-error samples for plotting scenario $(scenario_name).")

    t_s = Float64.(xrow.telemetry_axis[1:n])
    ex = Float64.(xrow.error_km[1:n])
    ey = Float64.(yrow.error_km[1:n])
    ez = Float64.(zrow.error_km[1:n])

    total_pos_err = sqrt.(ex .^ 2 .+ ey .^ 2 .+ ez .^ 2)
    running_total_rmse = sqrt.(cumsum(total_pos_err .^ 2) ./ collect(1:n))
    total_rmse = sqrt(mean(total_pos_err .^ 2))

    return (t_s=t_s, ex=ex, ey=ey, ez=ez, total_pos_err=total_pos_err, running_total_rmse=running_total_rmse, total_rmse=total_rmse)
end

function _iqr_inlier_mask(values::Vector{Float64}; whisker_scale::Float64=1.5)
    finite_idx = findall(isfinite, values)
    isempty(finite_idx) && return falses(length(values))

    finite_vals = values[finite_idx]
    q1 = quantile(finite_vals, 0.25)
    q3 = quantile(finite_vals, 0.75)
    iqr = q3 - q1
    if !isfinite(iqr) || iqr <= 0.0
        mask = falses(length(values))
        mask[finite_idx] .= true
        return mask
    end

    lo = q1 - whisker_scale * iqr
    hi = q3 + whisker_scale * iqr
    mask = falses(length(values))
    @inbounds for i in eachindex(values)
        v = values[i]
        mask[i] = isfinite(v) && lo <= v <= hi
    end
    return mask
end

function _filter_position_series_for_plot(series)
    mask = _iqr_inlier_mask(series.total_pos_err)
    sum(mask) >= 2 || (mask .= true)
    return (
        t_s=series.t_s[mask], ex=series.ex[mask], ey=series.ey[mask], ez=series.ez[mask],
        total_pos_err=series.total_pos_err[mask],
        running_total_rmse=sqrt.(cumsum(series.total_pos_err[mask] .^ 2) ./ collect(1:sum(mask))),
        total_rmse=series.total_rmse, inlier_count=sum(mask), total_count=length(mask),
    )
end

raw_series = _extract_position_error_series(result.errors, "cygnss_48hr_pvt")
series = _filter_position_series_for_plot(raw_series)
println("IQR filter kept $(series.inlier_count)/$(series.total_count) samples (dropped $(series.total_count - series.inlier_count) spikes)")

t_hr = series.t_s ./ 3600.0
err_fig = plot(
    t_hr, series.ex;
    label="x error", lw=1.2, alpha=0.9, color=:steelblue,
    title="CYGNSS 48hr Reference Replication: Per-Axis Position Error (IQR-Filtered)",
    xlabel="Time (hr)", ylabel="Error (km)",
    legend=:topleft, size=(1000, 550), titlefontsize=11,
)
plot!(err_fig, t_hr, series.ey; label="y error", lw=1.2, alpha=0.9, color=:darkorange)
plot!(err_fig, t_hr, series.ez; label="z error", lw=1.2, alpha=0.9, color=:seagreen)
plot!(err_fig, t_hr, series.running_total_rmse; label="total running RMSE", lw=2.0, color=:black)
hline!(err_fig, [series.total_rmse]; label="total RMSE (unfiltered)", ls=:dash, color=:black)
err_path = joinpath(PLOTS_DIR, "cygnss_48hr_reference_replication_error_timeseries.png")
savefig(err_fig, err_path)
println("error time series plot: $(err_path)")
