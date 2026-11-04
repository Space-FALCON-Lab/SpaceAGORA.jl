using SpaceAGORA
using Arrow
using DataFrames
using Dates
using LinearAlgebra
using Plots
using Printf
using StaticArrays
using Statistics

const TV = SpaceAGORA.TelemetryVerification

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const OUTDIR = joinpath(REPO_ROOT, "output", "cygnss_96hr_analysis")
const CYGNSS_48HR_PATH = joinpath(REPO_ROOT, "data", "telemetry", "CYGNSS", "cygnss_data_48hr.feather")
const CYG04_96HR_PATH = joinpath(REPO_ROOT, "data", "telemetry", "CYGNSS", "cyg04_nasa_pvt_96hr.feather")
const EARTH_HARMONICS_FILE = "data/Gravity_harmonics_data/EarthGGM05C.csv"

Base.@kwdef struct TrajectorySeries
    name::String
    initial_time::TV.InitialTime
    t_s::Vector{Float64}
    x_m::Vector{Float64}
    y_m::Vector{Float64}
    z_m::Vector{Float64}
    vx_mps::Vector{Float64}
    vy_mps::Vector{Float64}
    vz_mps::Vector{Float64}
end

@inline function _initial_time_48hr()::TV.InitialTime
    return TV.InitialTime(
        year=Int32(2025),
        month=Int16(6),
        day=Int16(6),
        hour=Int16(0),
        minute=Int16(0),
        second=Float32(0.0)
    )
end

function _initial_time_from_unix(unix_seconds::Float64)::TV.InitialTime
    dt = Dates.unix2datetime(unix_seconds)
    return TV.InitialTime(
        year=Int32(year(dt)),
        month=Int16(month(dt)),
        day=Int16(day(dt)),
        hour=Int16(hour(dt)),
        minute=Int16(minute(dt)),
        second=Float32(second(dt) + millisecond(dt) / 1.0e3)
    )
end

function _load_48hr_inertial_series()::TrajectorySeries
    df = DataFrame(Arrow.Table(CYGNSS_48HR_PATH))
    perm = sortperm(Float64.(df.time))
    df = df[perm, :]
    return TrajectorySeries(
        name="CYGNSS 48hr",
        initial_time=_initial_time_48hr(),
        t_s=Float64.(df.time),
        x_m=Float64.(df.pos_ii_1),
        y_m=Float64.(df.pos_ii_2),
        z_m=Float64.(df.pos_ii_3),
        vx_mps=Float64.(df.vel_ii_1),
        vy_mps=Float64.(df.vel_ii_2),
        vz_mps=Float64.(df.vel_ii_3)
    )
end

function _load_cyg04_inertial_series()::TrajectorySeries
    df = DataFrame(Arrow.Table(CYG04_96HR_PATH))
    perm = sortperm(Float64.(df.time))
    df = df[perm, :]
    return TrajectorySeries(
        name="CYG04 96hr",
        initial_time=_initial_time_from_unix(Float64(df[1, :pvt_unix_seconds])),
        t_s=Float64.(df.time) .- Float64(df[1, :time]),
        x_m=Float64.(df.pos_ii_1),
        y_m=Float64.(df.pos_ii_2),
        z_m=Float64.(df.pos_ii_3),
        vx_mps=Float64.(df.vel_ii_1),
        vy_mps=Float64.(df.vel_ii_2),
        vz_mps=Float64.(df.vel_ii_3)
    )
end

function _series_to_planet_fixed(series::TrajectorySeries; planet_name::String="earth")::TrajectorySeries
    et0 = TV._initial_time_et(series.initial_time)
    n = length(series.t_s)
    x_m = Vector{Float64}(undef, n)
    y_m = Vector{Float64}(undef, n)
    z_m = Vector{Float64}(undef, n)
    vx_mps = Vector{Float64}(undef, n)
    vy_mps = Vector{Float64}(undef, n)
    vz_mps = Vector{Float64}(undef, n)
    @inbounds for i in eachindex(series.t_s)
        r_pf, v_pf = TV._j2000_to_planet_fixed_state(
            planet_name,
            et0 + series.t_s[i],
            SVector{3, Float64}(series.x_m[i], series.y_m[i], series.z_m[i]),
            SVector{3, Float64}(series.vx_mps[i], series.vy_mps[i], series.vz_mps[i])
        )
        x_m[i] = r_pf[1]
        y_m[i] = r_pf[2]
        z_m[i] = r_pf[3]
        vx_mps[i] = v_pf[1]
        vy_mps[i] = v_pf[2]
        vz_mps[i] = v_pf[3]
    end
    return TrajectorySeries(
        name=series.name * " planet-fixed",
        initial_time=series.initial_time,
        t_s=copy(series.t_s),
        x_m=x_m,
        y_m=y_m,
        z_m=z_m,
        vx_mps=vx_mps,
        vy_mps=vy_mps,
        vz_mps=vz_mps
    )
end

function _interp_series_linear(x::Vector{Float64}, y::Vector{Float64}, xq::Vector{Float64})::Vector{Float64}
    n = length(x)
    n == length(y) || throw(ArgumentError("x/y length mismatch in interpolation"))
    n >= 2 || throw(ArgumentError("Need at least 2 samples for interpolation"))
    out = similar(xq, Float64)
    j = 1
    @inbounds for i in eachindex(xq)
        xi = xq[i]
        if xi <= x[1]
            out[i] = y[1]
            continue
        elseif xi >= x[end]
            out[i] = y[end]
            continue
        end
        while j < n - 1 && x[j + 1] < xi
            j += 1
        end
        x0 = x[j]
        x1 = x[j + 1]
        y0 = y[j]
        y1 = y[j + 1]
        α = (xi - x0) / (x1 - x0)
        out[i] = (1.0 - α) * y0 + α * y1
    end
    return out
end

function _position_error_against_reference(ref::TrajectorySeries, other::TrajectorySeries; tf_s::Float64)
    mask = ref.t_s .<= tf_s
    t = ref.t_s[mask]
    x_ref = ref.x_m[mask]
    y_ref = ref.y_m[mask]
    z_ref = ref.z_m[mask]
    x_other = _interp_series_linear(other.t_s, other.x_m, t)
    y_other = _interp_series_linear(other.t_s, other.y_m, t)
    z_other = _interp_series_linear(other.t_s, other.z_m, t)
    dx = x_other .- x_ref
    dy = y_other .- y_ref
    dz = z_other .- z_ref
    total_m = sqrt.(dx .^ 2 .+ dy .^ 2 .+ dz .^ 2)
    return (t_s=t, dx_m=dx, dy_m=dy, dz_m=dz, total_m=total_m)
end

function _downsample_indices(n::Int; max_points::Int=12000)::Vector{Int}
    if n <= max_points
        return collect(1:n)
    end
    step = ceil(Int, n / max_points)
    idx = collect(1:step:n)
    idx[end] == n || push!(idx, n)
    return idx
end

function _robust_limits(values...; lower_q::Float64=0.01, upper_q::Float64=0.99, padding_frac::Float64=0.08)
    data = Float64[]
    for v in values
        append!(data, vec(Float64.(v)))
    end
    finite_data = filter(isfinite, data)
    isempty(finite_data) && return nothing
    lo = quantile(finite_data, lower_q)
    hi = quantile(finite_data, upper_q)
    if !isfinite(lo) || !isfinite(hi)
        return nothing
    end
    if hi <= lo
        center = 0.5 * (lo + hi)
        halfspan = max(abs(center) * 0.1, 1.0)
        return (center - halfspan, center + halfspan)
    end
    span = hi - lo
    pad = max(span * padding_frac, 1.0e-9)
    return (lo - pad, hi + pad)
end

function _save_trajectory_overlay_plot(outpath::String, series_a::TrajectorySeries, series_b::TrajectorySeries; title::String)
    idx_a = _downsample_indices(length(series_a.t_s))
    idx_b = _downsample_indices(length(series_b.t_s))
    xlims_km = _robust_limits(series_a.x_m .* 1e-3, series_b.x_m .* 1e-3)
    ylims_km = _robust_limits(series_a.y_m .* 1e-3, series_b.y_m .* 1e-3)
    zlims_km = _robust_limits(series_a.z_m .* 1e-3, series_b.z_m .* 1e-3)
    p = plot(
        series_a.x_m[idx_a] .* 1e-3,
        series_a.y_m[idx_a] .* 1e-3,
        series_a.z_m[idx_a] .* 1e-3;
        label=series_a.name,
        lw=1.5,
        xlabel="x (km)",
        ylabel="y (km)",
        zlabel="z (km)",
        title=title,
        legend=:best,
        xlims=xlims_km,
        ylims=ylims_km,
        zlims=zlims_km
    )
    plot!(
        p,
        series_b.x_m[idx_b] .* 1e-3,
        series_b.y_m[idx_b] .* 1e-3,
        series_b.z_m[idx_b] .* 1e-3;
        label=series_b.name,
        lw=1.5
    )
    savefig(p, outpath)
    return outpath
end

function _save_error_plot(outpath::String, err_inertial, err_pf; title_prefix::String)
    t_hr = err_inertial.t_s ./ 3600.0
    inertial_axis_lims = _robust_limits(err_inertial.dx_m .* 1e-3, err_inertial.dy_m .* 1e-3, err_inertial.dz_m .* 1e-3)
    inertial_total_lims = _robust_limits(err_inertial.total_m .* 1e-3; lower_q=0.0, upper_q=0.99)
    pf_axis_lims = _robust_limits(err_pf.dx_m .* 1e-3, err_pf.dy_m .* 1e-3, err_pf.dz_m .* 1e-3)
    pf_total_lims = _robust_limits(err_pf.total_m .* 1e-3; lower_q=0.0, upper_q=0.99)
    p1 = plot(t_hr, err_inertial.dx_m .* 1e-3; label="dx", xlabel="Time (hr)", ylabel="km", title="$title_prefix Inertial Axis Error", lw=1.2, ylims=inertial_axis_lims)
    plot!(p1, t_hr, err_inertial.dy_m .* 1e-3; label="dy", lw=1.2)
    plot!(p1, t_hr, err_inertial.dz_m .* 1e-3; label="dz", lw=1.2)
    p2 = plot(t_hr, err_inertial.total_m .* 1e-3; label="total", xlabel="Time (hr)", ylabel="km", title="$title_prefix Inertial Total Error", lw=1.5, color=:black, ylims=inertial_total_lims)
    p3 = plot(t_hr, err_pf.dx_m .* 1e-3; label="dx", xlabel="Time (hr)", ylabel="km", title="$title_prefix Planet-Fixed Axis Error", lw=1.2, ylims=pf_axis_lims)
    plot!(p3, t_hr, err_pf.dy_m .* 1e-3; label="dy", lw=1.2)
    plot!(p3, t_hr, err_pf.dz_m .* 1e-3; label="dz", lw=1.2)
    p4 = plot(t_hr, err_pf.total_m .* 1e-3; label="total", xlabel="Time (hr)", ylabel="km", title="$title_prefix Planet-Fixed Total Error", lw=1.5, color=:black, ylims=pf_total_lims)
    fig = plot(p1, p2, p3, p4; layout=(2, 2), size=(1400, 900))
    savefig(fig, outpath)
    return outpath
end

function _telemetry_arrow_from_series(series::TrajectorySeries, initial_time::TV.InitialTime, outpath::String)
    planet = TV._planet_from_name("earth")
    r_km = sqrt.(series.x_m .^ 2 .+ series.y_m .^ 2 .+ series.z_m .^ 2) .* 1e-3
    alt_km = r_km .- planet.Rp_e * 1e-3
    oe0 = TV.rvtoorbitalelement(
        SVector{3, Float64}(series.x_m[1], series.y_m[1], series.z_m[1]),
        SVector{3, Float64}(series.vx_mps[1], series.vy_mps[1], series.vz_mps[1]),
        planet
    )
    df = DataFrame(
        time_s=series.t_s,
        altitude_km=alt_km,
        x_km=series.x_m .* 1e-3,
        y_km=series.y_m .* 1e-3,
        z_km=series.z_m .* 1e-3,
        sma_km=fill(oe0[1] * 1e-3, length(series.t_s)),
        ecc=fill(oe0[2], length(series.t_s)),
        inc_deg=fill(rad2deg(oe0[3]), length(series.t_s)),
        aop_deg=fill(rad2deg(oe0[5]), length(series.t_s)),
        raan_deg=fill(rad2deg(oe0[4]), length(series.t_s)),
        ta_deg=fill(rad2deg(oe0[6]), length(series.t_s)),
        x_ic_km=fill(series.x_m[1] * 1e-3, length(series.t_s)),
        y_ic_km=fill(series.y_m[1] * 1e-3, length(series.t_s)),
        z_ic_km=fill(series.z_m[1] * 1e-3, length(series.t_s)),
        vx_ic_kmps=fill(series.vx_mps[1] * 1e-3, length(series.t_s)),
        vy_ic_kmps=fill(series.vy_mps[1] * 1e-3, length(series.t_s)),
        vz_ic_kmps=fill(series.vz_mps[1] * 1e-3, length(series.t_s))
    )
    Arrow.write(outpath, df)
    return outpath
end

function _base_time_aligned_cfg(name::String, telemetry_path::String, initial_time::TV.InitialTime)::TV.TimeAlignedScenarioConfig
    return TV.TimeAlignedScenarioConfig(
        name=name,
        planet_name="earth",
        telemetry_path=telemetry_path,
        telemetry_time_col="time_s",
        telemetry_altitude_col="altitude_km",
        telemetry_x_col="x_km",
        telemetry_y_col="y_km",
        telemetry_z_col="z_km",
        telemetry_sma_col="sma_km",
        telemetry_ecc_col="ecc",
        telemetry_inc_col="inc_deg",
        telemetry_aop_col="aop_deg",
        telemetry_raan_col="raan_deg",
        telemetry_ta_col="ta_deg",
        telemetry_x_ic_col="x_ic_km",
        telemetry_y_ic_col="y_ic_km",
        telemetry_z_ic_col="z_ic_km",
        telemetry_vx_ic_col="vx_ic_kmps",
        telemetry_vy_ic_col="vy_ic_kmps",
        telemetry_vz_ic_col="vz_ic_kmps",
        max_points_quick=200000,
        max_points_full=200000,
        min_eval_points=2,
        units_x="s",
        units_y=Dict(
            "altitude_time" => "km",
            "state_x_time" => "km",
            "state_y_time" => "km",
            "state_z_time" => "km"
        ),
        tolerances_quick=Dict(
            "altitude_time" => (max_abs_km=1.0e6, max_nmae=1.0e6, max_rmse_km=1.0e6),
            "state_x_time" => (max_abs_km=1.0e6, max_nmae=1.0e6, max_rmse_km=1.0e6),
            "state_y_time" => (max_abs_km=1.0e6, max_nmae=1.0e6, max_rmse_km=1.0e6),
            "state_z_time" => (max_abs_km=1.0e6, max_nmae=1.0e6, max_rmse_km=1.0e6)
        ),
        tolerances_full=Dict(
            "altitude_time" => (max_abs_km=1.0e6, max_nmae=1.0e6, max_rmse_km=1.0e6),
            "state_x_time" => (max_abs_km=1.0e6, max_nmae=1.0e6, max_rmse_km=1.0e6),
            "state_y_time" => (max_abs_km=1.0e6, max_nmae=1.0e6, max_rmse_km=1.0e6),
            "state_z_time" => (max_abs_km=1.0e6, max_nmae=1.0e6, max_rmse_km=1.0e6)
        ),
        initial_time=initial_time,
        spacecraft=TV.SpacecraftConfig(
            bus_dims=(2.05e-1, 3.7e-1, 0.8e-1),
            panel_dims=(10e-3, 28.5e-3, 0.0001),
            bus_mass_kg=29.0,
            panel_mass_each_kg=0.0,
            panel_offset_y_m=2.45,
            prop_mass_kg=0.0,
            id=1002
        ),
        gravity_model=:inverse_squared,
        gravity_harmonics_degree=50,
        gravity_harmonics_order=50,
        gravity_harmonics_file=joinpath(REPO_ROOT, EARTH_HARMONICS_FILE),
        nbody_bodies=String[],
        srp_enabled=false,
        srp_cr=1.3,
        srp_area_m2=0.0,
        drag_enabled=false,
        include_wind=false,
        orbit_altitude_mode=:oblate,
        cartesian_ic_frame=:inertial,
        comparison_frame=:inertial,
        comparison_mode=:time_aligned_state,
        atmosphere_truth=TV.AtmosphereTruthConfig(
            assumption_id="earth_gmat_gram_deterministic_v1",
            atmosphere_model="GRAM",
            atmosphere_dataset="MERRA2 All Mean",
            space_weather_model="EarthGRAM MERRA2 climatology (deterministic)",
            solar_flux_model="EarthGRAM/MERRA2 epoch-fixed defaults",
            gram_seed=1001,
            gram_perturbation_scales=(0.0, 0.0, 0.0, 0.0),
            gram_offline_surrogate="auto",
            gram_static_grid=false,
            gram_track_cache=false,
            gram_global_lock="on"
        ),
        calibration=TV.CalibrationConfig(enabled=false),
        EI_km=300.0
    )
end

function _run_48hr_sim(cfg::TV.TimeAlignedScenarioConfig, ic, label::String)
    telemetry = TV._load_time_aligned_telemetry(cfg, cfg.max_points_quick)
    mission_time_s = max(telemetry.time_s[end] - telemetry.time_s[1], 1.0)
    args = TV._make_time_aligned_args(cfg, mission_time_s, ic)
    args = TV._with_study_settings(args; quick=true)
    run = withenv(
        "SPACEAGORA_TELEMETRY_SOLVER_MODE" => "dp8",
        "SPACEAGORA_TELEMETRY_DT_MAX_ORBIT" => "10.0",
        "SPACEAGORA_TELEMETRY_RELTOL_ORBIT" => "1e-12",
        "SPACEAGORA_TELEMETRY_ABSTOL_ORBIT" => "1e-12",
        "SPACEAGORA_TELEMETRY_RELTOL_ATM" => "1e-12",
        "SPACEAGORA_TELEMETRY_ABSTOL_ATM" => "1e-12"
    ) do
        TV._run_simulation_dataframe(args, label, cfg.atmosphere_truth, :quick)
    end
    rows, errors = TV._time_aligned_rows_errors(cfg, args, run.results_df, telemetry)
    return (args=args, run=run, rows=rows, errors=errors)
end

function _results_df_to_series(name::String, initial_time::TV.InitialTime, results_df::DataFrame)::TrajectorySeries
    return TrajectorySeries(
        name=name,
        initial_time=initial_time,
        t_s=Float64.(results_df.time),
        x_m=Float64.(results_df.sc1_pos_1),
        y_m=Float64.(results_df.sc1_pos_2),
        z_m=Float64.(results_df.sc1_pos_3),
        vx_mps=Float64.(results_df.sc1_vel_1),
        vy_mps=Float64.(results_df.sc1_vel_2),
        vz_mps=Float64.(results_df.sc1_vel_3)
    )
end

function _save_sim_overlay_plot(outpath::String, ref::TrajectorySeries, sim_a::TrajectorySeries, sim_b::TrajectorySeries; title::String)
    idx_ref = _downsample_indices(length(ref.t_s))
    idx_a = _downsample_indices(length(sim_a.t_s))
    idx_b = _downsample_indices(length(sim_b.t_s))
    xlims_km = _robust_limits(ref.x_m .* 1e-3, sim_a.x_m .* 1e-3, sim_b.x_m .* 1e-3)
    ylims_km = _robust_limits(ref.y_m .* 1e-3, sim_a.y_m .* 1e-3, sim_b.y_m .* 1e-3)
    zlims_km = _robust_limits(ref.z_m .* 1e-3, sim_a.z_m .* 1e-3, sim_b.z_m .* 1e-3)
    p = plot(
        ref.x_m[idx_ref] .* 1e-3,
        ref.y_m[idx_ref] .* 1e-3,
        ref.z_m[idx_ref] .* 1e-3;
        label=ref.name,
        lw=1.5,
        xlabel="x (km)",
        ylabel="y (km)",
        zlabel="z (km)",
        title=title,
        legend=:best,
        xlims=xlims_km,
        ylims=ylims_km,
        zlims=zlims_km
    )
    plot!(p, sim_a.x_m[idx_a] .* 1e-3, sim_a.y_m[idx_a] .* 1e-3, sim_a.z_m[idx_a] .* 1e-3; label=sim_a.name, lw=1.5)
    plot!(p, sim_b.x_m[idx_b] .* 1e-3, sim_b.y_m[idx_b] .* 1e-3, sim_b.z_m[idx_b] .* 1e-3; label=sim_b.name, lw=1.5)
    savefig(p, outpath)
    return outpath
end

function _save_sim_error_plot(outpath::String, err_inertial_a, err_inertial_b, err_pf_a, err_pf_b; title::String)
    t_hr_a = err_inertial_a.t_s ./ 3600.0
    t_hr_b = err_inertial_b.t_s ./ 3600.0
    inertial_total_lims = _robust_limits(err_inertial_a.total_m .* 1e-3, err_inertial_b.total_m .* 1e-3; lower_q=0.0, upper_q=0.99)
    inertial_axis_lims = _robust_limits(
        err_inertial_a.dx_m .* 1e-3, err_inertial_a.dy_m .* 1e-3, err_inertial_a.dz_m .* 1e-3,
        err_inertial_b.dx_m .* 1e-3, err_inertial_b.dy_m .* 1e-3, err_inertial_b.dz_m .* 1e-3
    )
    pf_total_lims = _robust_limits(err_pf_a.total_m .* 1e-3, err_pf_b.total_m .* 1e-3; lower_q=0.0, upper_q=0.99)
    pf_axis_lims = _robust_limits(
        err_pf_a.dx_m .* 1e-3, err_pf_a.dy_m .* 1e-3, err_pf_a.dz_m .* 1e-3,
        err_pf_b.dx_m .* 1e-3, err_pf_b.dy_m .* 1e-3, err_pf_b.dz_m .* 1e-3
    )
    p1 = plot(t_hr_a, err_inertial_a.total_m .* 1e-3; label="48hr IC", xlabel="Time (hr)", ylabel="km", title="$title Inertial Total Error", lw=1.5, ylims=inertial_total_lims)
    plot!(p1, t_hr_b, err_inertial_b.total_m .* 1e-3; label="CYG04 IC", lw=1.5)
    p2 = plot(t_hr_a, err_inertial_a.dx_m .* 1e-3; label="48hr IC dx", xlabel="Time (hr)", ylabel="km", title="$title Inertial Axis Errors", lw=1.2, ylims=inertial_axis_lims)
    plot!(p2, t_hr_a, err_inertial_a.dy_m .* 1e-3; label="48hr IC dy", lw=1.2)
    plot!(p2, t_hr_a, err_inertial_a.dz_m .* 1e-3; label="48hr IC dz", lw=1.2)
    plot!(p2, t_hr_b, err_inertial_b.dx_m .* 1e-3; label="CYG04 IC dx", lw=1.2, ls=:dash)
    plot!(p2, t_hr_b, err_inertial_b.dy_m .* 1e-3; label="CYG04 IC dy", lw=1.2, ls=:dash)
    plot!(p2, t_hr_b, err_inertial_b.dz_m .* 1e-3; label="CYG04 IC dz", lw=1.2, ls=:dash)
    p3 = plot(t_hr_a, err_pf_a.total_m .* 1e-3; label="48hr IC", xlabel="Time (hr)", ylabel="km", title="$title Planet-Fixed Total Error", lw=1.5, ylims=pf_total_lims)
    plot!(p3, t_hr_b, err_pf_b.total_m .* 1e-3; label="CYG04 IC", lw=1.5)
    p4 = plot(t_hr_a, err_pf_a.dx_m .* 1e-3; label="48hr IC dx", xlabel="Time (hr)", ylabel="km", title="$title Planet-Fixed Axis Errors", lw=1.2, ylims=pf_axis_lims)
    plot!(p4, t_hr_a, err_pf_a.dy_m .* 1e-3; label="48hr IC dy", lw=1.2)
    plot!(p4, t_hr_a, err_pf_a.dz_m .* 1e-3; label="48hr IC dz", lw=1.2)
    plot!(p4, t_hr_b, err_pf_b.dx_m .* 1e-3; label="CYG04 IC dx", lw=1.2, ls=:dash)
    plot!(p4, t_hr_b, err_pf_b.dy_m .* 1e-3; label="CYG04 IC dy", lw=1.2, ls=:dash)
    plot!(p4, t_hr_b, err_pf_b.dz_m .* 1e-3; label="CYG04 IC dz", lw=1.2, ls=:dash)
    fig = plot(p1, p2, p3, p4; layout=(2, 2), size=(1400, 900))
    savefig(fig, outpath)
    return outpath
end

function _mean_axis_rmse(rows)::Float64
    return mean(Float64.(getfield.(rows, :rmse_km)))
end

function main()
    mkpath(OUTDIR)

    # Standalone script entrypoints need to load the Earth SPICE kernels
    # before utc2et/sxform-based frame conversions.
    TV._planet_from_name("earth")

    ref48 = _load_48hr_inertial_series()
    cyg04 = _load_cyg04_inertial_series()
    ref48_pf = _series_to_planet_fixed(ref48)
    cyg04_pf = _series_to_planet_fixed(cyg04)

    overlap_tf_s = 48.0 * 3600.0
    file_err_inertial = _position_error_against_reference(ref48, cyg04; tf_s=overlap_tf_s)
    file_err_pf = _position_error_against_reference(ref48_pf, cyg04_pf; tf_s=overlap_tf_s)

    file_inertial_plot = _save_trajectory_overlay_plot(
        joinpath(OUTDIR, "cyg04_vs_48hr_inertial_trajectory.png"),
        ref48,
        cyg04;
        title="CYG04 96hr vs 48hr Trajectory (Inertial)"
    )
    file_pf_plot = _save_trajectory_overlay_plot(
        joinpath(OUTDIR, "cyg04_vs_48hr_planet_fixed_trajectory.png"),
        ref48_pf,
        cyg04_pf;
        title="CYG04 96hr vs 48hr Trajectory (Planet-Fixed)"
    )
    file_error_plot = _save_error_plot(
        joinpath(OUTDIR, "cyg04_vs_48hr_overlap_errors.png"),
        file_err_inertial,
        file_err_pf;
        title_prefix="CYG04 vs 48hr (First 48hr)"
    )

    telemetry_path = _telemetry_arrow_from_series(ref48, ref48.initial_time, joinpath(OUTDIR, "cygnss_48hr_reference.arrow"))
    cfg_48_ic = _base_time_aligned_cfg("cygnss_48hr_from_48_ic", telemetry_path, ref48.initial_time)
    cfg_cyg04_ic = _base_time_aligned_cfg("cygnss_48hr_from_cyg04_ic", telemetry_path, ref48.initial_time)
    telemetry = TV._load_time_aligned_telemetry(cfg_48_ic, cfg_48_ic.max_points_quick)
    ic48 = TV._initial_condition_from_time_aligned_telemetry(cfg_48_ic, telemetry)
    ic_cyg04 = TV.CartesianInitialCondition(
        [cyg04.x_m[1], cyg04.y_m[1], cyg04.z_m[1]],
        [cyg04.vx_mps[1], cyg04.vy_mps[1], cyg04.vz_mps[1]]
    )

    println("Running 48-hour simulation from 48-hour IC ...")
    sim48 = _run_48hr_sim(cfg_48_ic, ic48, cfg_48_ic.name)
    println("Running 48-hour simulation from CYG04 96-hour IC ...")
    sim_cyg04 = _run_48hr_sim(cfg_cyg04_ic, ic_cyg04, cfg_cyg04_ic.name)

    sim48_series = _results_df_to_series("Sim from 48hr IC", ref48.initial_time, sim48.run.results_df)
    sim_cyg04_series = _results_df_to_series("Sim from CYG04 IC", ref48.initial_time, sim_cyg04.run.results_df)
    sim48_pf = _series_to_planet_fixed(sim48_series)
    sim_cyg04_pf = _series_to_planet_fixed(sim_cyg04_series)
    sim48_err_inertial = _position_error_against_reference(ref48, sim48_series; tf_s=overlap_tf_s)
    sim_cyg04_err_inertial = _position_error_against_reference(ref48, sim_cyg04_series; tf_s=overlap_tf_s)
    sim48_err_pf = _position_error_against_reference(ref48_pf, sim48_pf; tf_s=overlap_tf_s)
    sim_cyg04_err_pf = _position_error_against_reference(ref48_pf, sim_cyg04_pf; tf_s=overlap_tf_s)

    sim_inertial_plot = _save_sim_overlay_plot(
        joinpath(OUTDIR, "sim_48hr_ic_vs_cyg04_ic_inertial.png"),
        ref48,
        sim48_series,
        sim_cyg04_series;
        title="48hr SpaceAGORA Runs vs 48hr Reference (Inertial)"
    )
    sim_pf_plot = _save_sim_overlay_plot(
        joinpath(OUTDIR, "sim_48hr_ic_vs_cyg04_ic_planet_fixed.png"),
        ref48_pf,
        sim48_pf,
        sim_cyg04_pf;
        title="48hr SpaceAGORA Runs vs 48hr Reference (Planet-Fixed)"
    )
    sim_error_plot = _save_sim_error_plot(
        joinpath(OUTDIR, "sim_48hr_ic_vs_cyg04_ic_errors.png"),
        sim48_err_inertial,
        sim_cyg04_err_inertial,
        sim48_err_pf,
        sim_cyg04_err_pf;
        title="48hr SpaceAGORA Runs"
    )

    sim48_rmse = _mean_axis_rmse(sim48.rows[2:4])
    sim_cyg04_rmse = _mean_axis_rmse(sim_cyg04.rows[2:4])
    file_overlap_rmse_km = sqrt(mean(file_err_inertial.total_m .^ 2)) * 1e-3

    println()
    println("Saved plots:")
    println("  $(file_inertial_plot)")
    println("  $(file_pf_plot)")
    println("  $(file_error_plot)")
    println("  $(sim_inertial_plot)")
    println("  $(sim_pf_plot)")
    println("  $(sim_error_plot)")
    println()
    println(@sprintf("CYG04 vs 48hr first-48hr mean position difference: %.6f km", file_overlap_rmse_km))
    println(@sprintf("48hr reference run mean position-axis RMSE: %.6f km", sim48_rmse))
    println(@sprintf("CYG04-IC run mean position-axis RMSE: %.6f km", sim_cyg04_rmse))
end

main()
