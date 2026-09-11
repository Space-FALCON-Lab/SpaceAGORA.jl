@inline _tolerances_for(cfg::TimeAlignedScenarioConfig, profile::Symbol) = profile == :quick ? cfg.tolerances_quick : cfg.tolerances_full

@inline function _default_plots_outdir(out_summary::String, profile::Symbol)::String
    return normpath(joinpath(dirname(out_summary), "telemetry_plots_$(String(profile))"))
end

function _generate_plots(summary_csv::String, errors_csv::String, profile::Symbol)::String
    plot_script = joinpath(REPO_ROOT, "scripts", "plotting", "telemetry_orbit_accuracy_plots.jl")
    isfile(plot_script) || throw(ErrorException("Missing plotting script: $plot_script"))
    outdir = _default_plots_outdir(summary_csv, profile)
    plot_project = let agora_project = joinpath(REPO_ROOT, ".AGORA")
        isdir(agora_project) ? agora_project : REPO_ROOT
    end
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$plot_project $plot_script --summary=$summary_csv --errors=$errors_csv --outdir=$outdir`
    run(cmd)
    return outdir
end

@inline _axis_units(cfg::OrbitEventsScenarioConfig) = cfg.units_x
@inline _axis_units(cfg::TimeAlignedScenarioConfig) = cfg.units_x

@inline function _value_units(cfg::OrbitEventsScenarioConfig, event::String)
    if occursin("speed", event) || occursin("vx", event) || occursin("vy", event) || occursin("vz", event)
        return get(cfg.units_y, event, "km/s")
    end
    return get(cfg.units_y, event, "km")
end
@inline function _value_units(cfg::TimeAlignedScenarioConfig, event::String)
    if occursin("speed", event) || occursin("vx", event) || occursin("vy", event) || occursin("vz", event)
        return get(cfg.units_y, event, "km/s")
    end
    return get(cfg.units_y, event, "km")
end

@inline function _display_value_units(value_units::AbstractString)::String
    token = lowercase(strip(String(value_units)))
    return token == "km/s" ? "m/s" : String(value_units)
end

@inline function _display_value_scale(value_units::AbstractString)::Float64
    token = lowercase(strip(String(value_units)))
    return token == "km/s" ? 1e3 : 1.0
end

function _append_display_metric_columns!(summary_df::DataFrame)::Nothing
    nrow(summary_df) == 0 && return nothing
    hasproperty(summary_df, :value_units) || return nothing

    units_raw = [String(v) for v in summary_df.value_units]
    scales = [_display_value_scale(u) for u in units_raw]
    units_display = [_display_value_units(u) for u in units_raw]
    summary_df.value_units_display = units_display

    function _add_scaled_column!(src::Symbol, dst::Symbol)
        hasproperty(summary_df, src) || return nothing
        src_values = summary_df[!, src]
        out = Vector{Float64}(undef, length(src_values))
        @inbounds for i in eachindex(src_values)
            out[i] = Float64(src_values[i]) * scales[i]
        end
        summary_df[!, dst] = out
        return nothing
    end

    _add_scaled_column!(:mae_km, :mae_display)
    _add_scaled_column!(:rmse_km, :rmse_display)
    _add_scaled_column!(:max_abs_km, :max_abs_display)
    _add_scaled_column!(:p95_abs_km, :p95_abs_display)
    _add_scaled_column!(:bias_km, :bias_display)
    _add_scaled_column!(:limit_max_abs_km, :limit_max_abs_display)
    _add_scaled_column!(:limit_max_rmse_km, :limit_max_rmse_display)
    return nothing
end

function _append_display_error_columns!(errors_df::DataFrame, summary_df::DataFrame)::Nothing
    nrow(errors_df) == 0 && return nothing
    hasproperty(errors_df, :scenario) || return nothing
    hasproperty(errors_df, :event) || return nothing

    unit_map = Dict{Tuple{String, String}, Tuple{Float64, String}}()
    if nrow(summary_df) > 0 && hasproperty(summary_df, :scenario) && hasproperty(summary_df, :event) && hasproperty(summary_df, :value_units)
        for row in eachrow(summary_df)
            key = (String(row.scenario), String(row.event))
            units = String(row.value_units)
            unit_map[key] = (_display_value_scale(units), _display_value_units(units))
        end
    end

    n = nrow(errors_df)
    scales = Vector{Float64}(undef, n)
    units_display = Vector{String}(undef, n)
    @inbounds for i in 1:n
        key = (String(errors_df.scenario[i]), String(errors_df.event[i]))
        scale, unit = get(unit_map, key, (1.0, "km"))
        scales[i] = scale
        units_display[i] = unit
    end
    errors_df.value_units_display = units_display

    function _add_scaled_column!(src::Symbol, dst::Symbol)
        hasproperty(errors_df, src) || return nothing
        src_values = errors_df[!, src]
        out = Vector{Float64}(undef, length(src_values))
        @inbounds for i in eachindex(src_values)
            out[i] = Float64(src_values[i]) * scales[i]
        end
        errors_df[!, dst] = out
        return nothing
    end

    _add_scaled_column!(:telemetry_value_km, :telemetry_value_display)
    _add_scaled_column!(:sim_interp_value_km, :sim_interp_value_display)
    _add_scaled_column!(:error_km, :error_display)
    return nothing
end

@inline _source_file(cfg::OrbitEventsScenarioConfig, event::String) = event == "peri" ? cfg.telemetry_peri_path : cfg.telemetry_apo_path
@inline _source_file(cfg::TimeAlignedScenarioConfig, event::String) = cfg.telemetry_path

@inline _orbit_altitude_mode(cfg::OrbitEventsScenarioConfig) = String(cfg.orbit_altitude_mode)
@inline function _orbit_altitude_mode(cfg::TimeAlignedScenarioConfig)
    return cfg.comparison_mode == :orbit_events ? String(cfg.orbit_altitude_mode) : "n/a"
end

@inline _maneuver_count(cfg::OrbitEventsScenarioConfig) = length(cfg.maneuver_orbit_numbers)
@inline _maneuver_count(::TimeAlignedScenarioConfig) = 0

@inline _maneuver_replay_scale_mode(cfg::OrbitEventsScenarioConfig) = cfg.maneuver_replay_scale_mode
@inline _maneuver_replay_scale_mode(::TimeAlignedScenarioConfig) = "delta_v"

@inline _min_eval_points(cfg::OrbitEventsScenarioConfig) = cfg.min_eval_points
@inline _min_eval_points(cfg::TimeAlignedScenarioConfig) = cfg.min_eval_points

@inline function _scenario_status_extra(cfg::OrbitEventsScenarioConfig)::String
    return ", altitude_mode=$(String(cfg.orbit_altitude_mode)), maneuvers=$(length(cfg.maneuver_orbit_numbers))"
end
@inline function _scenario_status_extra(cfg::TimeAlignedScenarioConfig)::String
    if cfg.comparison_mode == :orbit_events
        return ", comparison_mode=orbit_events, altitude_mode=$(String(cfg.orbit_altitude_mode))"
    end
    return ", comparison_mode=time_aligned_state"
end

function _evaluate_thresholds(row, cfg::AbstractScenarioConfig, profile::Symbol)
    tmap = _tolerances_for(cfg, profile)
    event_name = String(row.event)
    tol = get(tmap, event_name, nothing)
    if tol === nothing && endswith(event_name, "_speed")
        base_event = replace(event_name, "_speed" => "")
        tol = get(tmap, base_event, nothing)
    end
    tol === nothing && throw(ArgumentError("No tolerance configured for event '$(row.event)' in scenario '$(row.scenario)'"))

    pass_points = row.n_sim >= _min_eval_points(cfg)
    pass_abs = row.max_abs_km <= tol.max_abs_km
    pass_nmae = row.nmae <= tol.max_nmae
    pass_rmse = row.rmse_km <= tol.max_rmse_km
    pass = pass_points && pass_abs && pass_nmae && pass_rmse

    return (
        pass=pass,
        pass_points=pass_points,
        pass_abs=pass_abs,
        pass_nmae=pass_nmae,
        pass_rmse=pass_rmse,
        limit_max_abs_km=tol.max_abs_km,
        limit_nmae=tol.max_nmae,
        limit_max_rmse_km=tol.max_rmse_km,
        min_eval_points=_min_eval_points(cfg)
    )
end
