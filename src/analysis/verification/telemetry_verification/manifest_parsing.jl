@inline function _safe_parse_bool(raw::AbstractString, default::Bool)::Bool
    token = lowercase(strip(raw))
    if token in ("1", "true", "yes", "on")
        return true
    elseif token in ("0", "false", "no", "off")
        return false
    end
    return default
end

@inline function _parse_positive_int_env(name::String, default::Int)::Int
    raw = strip(get(ENV, name, ""))
    isempty(raw) && return default
    parsed = try
        parse(Int, raw)
    catch
        throw(ArgumentError("$name must be a positive integer, got '$raw'"))
    end
    parsed > 0 || throw(ArgumentError("$name must be a positive integer, got $parsed"))
    return parsed
end

@inline function _parse_positive_float_env(name::String)::Union{Nothing, Float64}
    raw = strip(get(ENV, name, ""))
    isempty(raw) && return nothing
    parsed = try
        parse(Float64, raw)
    catch
        throw(ArgumentError("$name must be a positive float, got '$raw'"))
    end
    parsed > 0.0 || throw(ArgumentError("$name must be > 0, got $parsed"))
    return parsed
end

@inline function _telemetry_solver_maxiters(profile::Symbol)::Int
    default_val = profile == :quick ? TELEMETRY_SOLVER_MAXITERS_QUICK_DEFAULT : TELEMETRY_SOLVER_MAXITERS_FULL_DEFAULT
    return _parse_positive_int_env("SPACEAGORA_TELEMETRY_SOLVER_MAXITERS", default_val)
end

@inline function _telemetry_solver_retry_maxiters(base_maxiters::Int)::Int
    default_retry = max(base_maxiters * 4, base_maxiters + 1_000_000)
    retry = _parse_positive_int_env("SPACEAGORA_TELEMETRY_SOLVER_MAXITERS_RETRY", default_retry)
    return max(retry, base_maxiters + 1)
end

@inline function _telemetry_solver_mode()::String
    mode = strip(get(ENV, "SPACEAGORA_TELEMETRY_SOLVER_MODE", ""))
    if isempty(mode)
        mode = strip(get(ENV, "SPACEAGORA_SOLVER_MODE", ""))
    end
    return isempty(mode) ? "auto_stiff" : mode
end

@inline _is_maxiters_error(err)::Bool = occursin("MaxIters", sprint(showerror, err))

@inline function _require_key(tbl, key::String, context::String)
    haskey(tbl, key) || throw(ArgumentError("Missing '$key' in $context"))
    return tbl[key]
end

@inline function _require_table(tbl, key::String, context::String)
    value = _require_key(tbl, key, context)
    value isa AbstractDict || throw(ArgumentError("Expected table '$key' in $context"))
    return value
end

@inline function _require_str(tbl, key::String, context::String)::String
    return String(_require_key(tbl, key, context))
end

@inline function _require_float(tbl, key::String, context::String)::Float64
    return Float64(_require_key(tbl, key, context))
end

@inline function _require_int(tbl, key::String, context::String)::Int
    return Int(_require_key(tbl, key, context))
end

@inline function _optional_int(tbl, key::String, default::Int)::Int
    return haskey(tbl, key) ? Int(tbl[key]) : default
end

@inline function _optional_float(tbl, key::String, default::Float64)::Float64
    return haskey(tbl, key) ? Float64(tbl[key]) : default
end

@inline function _optional_bool(tbl, key::String, default::Bool)::Bool
    return haskey(tbl, key) ? Bool(tbl[key]) : default
end

@inline function _optional_str(tbl, key::String, default::String)::String
    return haskey(tbl, key) ? String(tbl[key]) : default
end

@inline function _optional_symbol_vector(tbl, key::String, default::Vector{Symbol})::Vector{Symbol}
    if !haskey(tbl, key)
        return copy(default)
    end
    raw = tbl[key]
    raw isa AbstractVector || throw(ArgumentError("Expected array '$key' in manifest scenario table."))
    out = Symbol[]
    for token in raw
        push!(out, Symbol(lowercase(strip(String(token)))))
    end
    isempty(out) && throw(ArgumentError("Calibration profiles cannot be empty."))
    return out
end

@inline function _require_str_vector(tbl, key::String, context::String)::Vector{String}
    raw = _require_key(tbl, key, context)
    raw isa AbstractVector || throw(ArgumentError("Expected array '$key' in $context"))
    return [String(v) for v in raw]
end

@inline function _optional_str_vector(tbl, key::String)::Vector{String}
    if !haskey(tbl, key)
        return String[]
    end
    raw = tbl[key]
    raw isa AbstractVector || throw(ArgumentError("Expected array '$key' in manifest scenario table."))
    return [String(v) for v in raw]
end

@inline function _optional_int64_vector(tbl, key::String)::Vector{Int64}
    if !haskey(tbl, key)
        return Int64[]
    end
    raw = tbl[key]
    raw isa AbstractVector || throw(ArgumentError("Expected array '$key' in manifest scenario table."))
    return Int64[Int64(v) for v in raw]
end

@inline function _optional_float64_vector(tbl, key::String)::Vector{Float64}
    if !haskey(tbl, key)
        return Float64[]
    end
    raw = tbl[key]
    raw isa AbstractVector || throw(ArgumentError("Expected array '$key' in manifest scenario table."))
    return Float64[Float64(v) for v in raw]
end

@inline function _parse_orbit_altitude_mode(raw::String, context::String)::Symbol
    key = lowercase(strip(raw))
    if key in ("vacuum", "spherical", "rp_e")
        return :vacuum
    elseif key in ("true", "oblate", "geodetic")
        return :oblate
    end
    throw(ArgumentError(
        "Unsupported orbit_altitude_mode='$raw' in $context; use vacuum|true."
    ))
end

@inline function _parse_element_frame(raw::String, context::String)::Symbol
    key = lowercase(strip(raw))
    if key in ("j2000", "eme2000")
        return :j2000
    elseif key in ("body_equator_inertial", "body_equator", "planet_equator")
        return :body_equator_inertial
    end
    throw(ArgumentError(
        "Unsupported element_frame='$raw' in $context; use j2000|body_equator_inertial."
    ))
end

@inline function _parse_time_aligned_comparison_mode(raw::String, context::String)::Symbol
    key = lowercase(strip(raw))
    if key in ("time_aligned_state", "time_aligned", "state")
        return :time_aligned_state
    elseif key in ("orbit_events", "apo_peri", "extrema")
        return :orbit_events
    end
    throw(ArgumentError(
        "Unsupported comparison_mode='$raw' in $context; use time_aligned_state|orbit_events."
    ))
end

@inline function _parse_reference_frame(raw::String, context::String)::Symbol
    key = lowercase(strip(raw))
    if key in ("inertial", "j2000", "eci", "earthmj2000eq")
        return :inertial
    elseif key in ("planet_fixed", "fixed", "ecef", "itrf", "itrf93")
        return :planet_fixed
    end
    throw(ArgumentError(
        "Unsupported reference frame '$raw' in $context; use inertial|planet_fixed."
    ))
end

function _parse_maneuver_config(tbl, context::String)
    if !haskey(tbl, "maneuvers")
        return (
            orbit_numbers=Int64[],
            orbit_numbers_campaign=Int64[],
            delta_v_mps=Float64[],
            replay_scale_mode="delta_v",
            flight_apoapsis_alt_m=Float64[],
            thrust_n=0.0,
            isp_s=0.0,
            guidance_rate_s=30.0,
            control_rate_s=10.0
        )
    end
    mtbl = _require_table(tbl, "maneuvers", context)
    orbit_numbers = _optional_int64_vector(mtbl, "orbit_numbers")
    delta_v_mps = _optional_float64_vector(mtbl, "delta_v_mps")
    isempty(orbit_numbers) && throw(ArgumentError("maneuvers.orbit_numbers must be non-empty in $context"))
    length(delta_v_mps) == length(orbit_numbers) || throw(ArgumentError(
        "maneuvers.delta_v_mps length ($(length(delta_v_mps))) must match maneuvers.orbit_numbers length ($(length(orbit_numbers))) in $context"
    ))
    any(v -> v <= 0, orbit_numbers) && throw(ArgumentError("maneuvers.orbit_numbers must be positive integers in $context"))

    # Flight maneuver histories are often numbered from the campaign origin
    # (e.g. orbit insertion), while the replay fires at the Nth apoapsis after
    # the scenario epoch. orbit_number_offset converts campaign numbers to
    # epoch-relative ones; burns that executed before the epoch (shifted
    # number < 1) are already represented by the initial condition and are
    # dropped here rather than replayed a second time.
    offset = _optional_int(mtbl, "orbit_number_offset", 0)
    offset >= 0 || throw(ArgumentError("maneuvers.orbit_number_offset must be >= 0 in $context"))
    # The campaign-numbered list (including pre-epoch burns) is kept for
    # diagnostics that operate on campaign orbit axes — the truth curve carries
    # jumps at every campaign burn regardless of what the replay fires.
    # Diagnostic replay scaling (see types.jl): flight apoapsis altitudes are
    # required per burn in "flight_apoapsis_ratio" mode and follow the same
    # pre-epoch drop filter as the burns themselves.
    replay_scale_mode = _optional_str(mtbl, "replay_scale_mode", "delta_v")
    replay_scale_mode in ("delta_v", "flight_apoapsis_ratio") || throw(ArgumentError(
        "maneuvers.replay_scale_mode must be \"delta_v\" or \"flight_apoapsis_ratio\" in $context"
    ))
    flight_apo_alt_km = _optional_float64_vector(mtbl, "flight_apoapsis_alt_km")
    if replay_scale_mode == "flight_apoapsis_ratio"
        length(flight_apo_alt_km) == length(orbit_numbers) || throw(ArgumentError(
            "maneuvers.flight_apoapsis_alt_km length ($(length(flight_apo_alt_km))) must match maneuvers.orbit_numbers length ($(length(orbit_numbers))) in $context"
        ))
        any(v -> !(isfinite(v) && v > 0.0), flight_apo_alt_km) && throw(ArgumentError(
            "maneuvers.flight_apoapsis_alt_km entries must be positive and finite in $context"
        ))
    elseif !isempty(flight_apo_alt_km)
        throw(ArgumentError(
            "maneuvers.flight_apoapsis_alt_km requires replay_scale_mode = \"flight_apoapsis_ratio\" in $context"
        ))
    end

    orbit_numbers_campaign = copy(orbit_numbers)
    if offset > 0
        keep = findall(o -> o - offset >= 1, orbit_numbers)
        n_dropped = length(orbit_numbers) - length(keep)
        n_dropped > 0 && println(
            "maneuver_offset context=$context offset=$offset dropped_pre_epoch_burns=$n_dropped"
        )
        orbit_numbers = Int64[orbit_numbers[i] - offset for i in keep]
        delta_v_mps = delta_v_mps[keep]
        if !isempty(flight_apo_alt_km)
            flight_apo_alt_km = flight_apo_alt_km[keep]
        end
    end
    return (
        orbit_numbers=orbit_numbers,
        orbit_numbers_campaign=orbit_numbers_campaign,
        delta_v_mps=delta_v_mps,
        replay_scale_mode=replay_scale_mode,
        flight_apoapsis_alt_m=Float64[v * 1000.0 for v in flight_apo_alt_km],
        thrust_n=_optional_float(mtbl, "thrust_n", 4.0),
        isp_s=_optional_float(mtbl, "isp_s", 220.0),
        guidance_rate_s=_optional_float(mtbl, "guidance_rate_s", 30.0),
        control_rate_s=_optional_float(mtbl, "control_rate_s", 10.0)
    )
end

@inline function _optional_float_tuple(tbl, key::String, n::Int, context::String)
    if !haskey(tbl, key)
        return nothing
    end
    raw = tbl[key]
    raw isa AbstractVector || throw(ArgumentError("Expected array '$key' in $context"))
    length(raw) == n || throw(ArgumentError("Expected '$key' length=$n in $context, got $(length(raw))"))
    values = ntuple(i -> Float64(raw[i]), n)
    return values
end

function _parse_atmosphere_truth_config(tbl, context::String)::AtmosphereTruthConfig
    if !haskey(tbl, "atmosphere_truth")
        return AtmosphereTruthConfig()
    end
    t = _require_table(tbl, "atmosphere_truth", context)
    assumption_id = _optional_str(t, "assumption_id", "gram_default")
    atmosphere_model = _require_str(t, "atmosphere_model", "$context.atmosphere_truth")
    atmosphere_dataset = _require_str(t, "atmosphere_dataset", "$context.atmosphere_truth")
    space_weather_model = _require_str(t, "space_weather_model", "$context.atmosphere_truth")
    solar_flux_model = _require_str(t, "solar_flux_model", "$context.atmosphere_truth")
    seed = _optional_int(t, "gram_seed", 1001)
    scales_raw = _optional_float_tuple(t, "gram_perturbation_scales", 4, "$context.atmosphere_truth")
    scales = scales_raw === nothing ? (0.0, 0.0, 0.0, 0.0) : (
        scales_raw[1], scales_raw[2], scales_raw[3], scales_raw[4]
    )
    min_step = haskey(t, "gram_min_relative_step_size") ? _optional_float(t, "gram_min_relative_step_size", 0.0) : nothing
    offline_surrogate = lowercase(_optional_str(t, "gram_offline_surrogate", "off"))
    offline_surrogate in ("off", "on", "auto") || throw(ArgumentError(
        "Unsupported atmosphere_truth.gram_offline_surrogate='$offline_surrogate' in $context; use off|on|auto."
    ))
    global_lock = lowercase(_optional_str(t, "gram_global_lock", "on"))
    global_lock in ("on", "off") || throw(ArgumentError(
        "Unsupported atmosphere_truth.gram_global_lock='$global_lock' in $context; use on|off."
    ))
    mars_mgcm_raw = _optional_float_tuple(t, "mars_mgcm_dust_levels", 3, "$context.atmosphere_truth")
    mars_dust_raw = _optional_float_tuple(t, "mars_dust_storm", 6, "$context.atmosphere_truth")
    mars_wind_raw = _optional_float_tuple(t, "mars_wind_scales", 2, "$context.atmosphere_truth")
    return AtmosphereTruthConfig(
        assumption_id=assumption_id,
        atmosphere_model=atmosphere_model,
        atmosphere_dataset=atmosphere_dataset,
        space_weather_model=space_weather_model,
        solar_flux_model=solar_flux_model,
        gram_seed=seed,
        gram_perturbation_scales=scales,
        gram_min_relative_step_size=min_step,
        gram_offline_surrogate=offline_surrogate,
        gram_static_grid=_optional_bool(t, "gram_static_grid", false),
        gram_track_cache=_optional_bool(t, "gram_track_cache", false),
        gram_global_lock=global_lock,
        mars_map_year=haskey(t, "mars_map_year") ? _optional_int(t, "mars_map_year", 0) : nothing,
        mars_mgcm_dust_levels=mars_mgcm_raw === nothing ? nothing : (
            mars_mgcm_raw[1], mars_mgcm_raw[2], mars_mgcm_raw[3]
        ),
        mars_dust_storm=mars_dust_raw === nothing ? nothing : (
            mars_dust_raw[1], mars_dust_raw[2], mars_dust_raw[3], mars_dust_raw[4], mars_dust_raw[5], mars_dust_raw[6]
        ),
        mars_f107=haskey(t, "mars_f107") ? _optional_float(t, "mars_f107", 0.0) : nothing,
        mars_wind_scales=mars_wind_raw === nothing ? nothing : (mars_wind_raw[1], mars_wind_raw[2]),
        mars_mola_heights=haskey(t, "mars_mola_heights") ? _optional_bool(t, "mars_mola_heights", true) : nothing,
        mars_min_max=haskey(t, "mars_min_max") ? _optional_int(t, "mars_min_max", 0) : nothing
    )
end

function _parse_calibration_config(tbl, context::String)::CalibrationConfig
    if !haskey(tbl, "calibration")
        return CalibrationConfig()
    end
    c = _require_table(tbl, "calibration", context)
    profiles = _optional_symbol_vector(c, "profiles", Symbol[:full])
    objective = lowercase(_optional_str(c, "objective", "mean_nmae"))
    objective in ("mean_nmae", "mean_rmse_km", "max_nmae") || throw(ArgumentError(
        "Unsupported calibration.objective='$objective' in $context; use mean_nmae|mean_rmse_km|max_nmae."
    ))
    cd_steps = _optional_int(c, "cd_scale_steps", 3)
    cr_steps = _optional_int(c, "cr_steps", 3)
    cd_steps > 0 || throw(ArgumentError("calibration.cd_scale_steps must be > 0 in $context"))
    cr_steps > 0 || throw(ArgumentError("calibration.cr_steps must be > 0 in $context"))
    cd_min = _optional_float(c, "cd_scale_min", 0.85)
    cd_max = _optional_float(c, "cd_scale_max", 1.15)
    cr_min = _optional_float(c, "cr_min", 1.15)
    cr_max = _optional_float(c, "cr_max", 1.45)
    cd_max >= cd_min || throw(ArgumentError("calibration.cd_scale_max must be >= cd_scale_min in $context"))
    cr_max >= cr_min || throw(ArgumentError("calibration.cr_max must be >= cr_min in $context"))
    return CalibrationConfig(
        enabled=_optional_bool(c, "enabled", false),
        profiles=profiles,
        search_on_quick_subset=_optional_bool(c, "search_on_quick_subset", true),
        fit_cd_scale=_optional_bool(c, "fit_cd_scale", true),
        cd_scale_min=cd_min,
        cd_scale_max=cd_max,
        cd_scale_steps=cd_steps,
        fit_cr=_optional_bool(c, "fit_cr", true),
        cr_min=cr_min,
        cr_max=cr_max,
        cr_steps=cr_steps,
        fit_bias=_optional_bool(c, "fit_bias", true),
        bias_abs_max_km=_optional_float(c, "bias_abs_max_km", 500.0),
        objective=objective
    )
end

@inline function _resolve_repo_path(path::String)::String
    return isabspath(path) ? path : normpath(joinpath(REPO_ROOT, path))
end

@inline function _parse_gravity_model(raw::String, context::String)::Symbol
    key = lowercase(strip(raw))
    if key in ("inverse_squared", "is2", "point_mass")
        return :inverse_squared
    elseif key in ("inverse_squared_j2", "is2_j2", "j2")
        return :inverse_squared_j2
    end
    throw(ArgumentError("Unsupported gravity_model='$raw' in $context"))
end

function _parse_initial_time(tbl, context::String)::InitialTime
    return InitialTime(
        year=_require_int(tbl, "year", context),
        month=_require_int(tbl, "month", context),
        day=_require_int(tbl, "day", context),
        hour=_require_int(tbl, "hour", context),
        minute=_require_int(tbl, "minute", context),
        second=_require_float(tbl, "second", context)
    )
end

function _parse_vec3(tbl, key::String, context::String)::NTuple{3, Float64}
    raw = _require_key(tbl, key, context)
    raw isa AbstractVector || throw(ArgumentError("Expected 3-element array '$key' in $context"))
    length(raw) == 3 || throw(ArgumentError("Expected 3 values for '$key' in $context, got $(length(raw))"))
    return (Float64(raw[1]), Float64(raw[2]), Float64(raw[3]))
end

function _parse_spacecraft_config(tbl, context::String)::SpacecraftConfig
    stbl = _require_table(tbl, "spacecraft", context)
    return SpacecraftConfig(
        bus_dims=_parse_vec3(stbl, "bus_dims_m", "$context.spacecraft"),
        panel_dims=_parse_vec3(stbl, "panel_dims_m", "$context.spacecraft"),
        bus_mass_kg=_require_float(stbl, "bus_mass_kg", "$context.spacecraft"),
        panel_mass_each_kg=_require_float(stbl, "panel_mass_each_kg", "$context.spacecraft"),
        panel_offset_y_m=_require_float(stbl, "panel_offset_y_m", "$context.spacecraft"),
        prop_mass_kg=_require_float(stbl, "prop_mass_kg", "$context.spacecraft"),
        id=Int64(_require_int(stbl, "id", "$context.spacecraft"))
    )
end

function _parse_units(tbl, events::Vector{String}, context::String)::Tuple{String, Dict{String, String}}
    utbl = _require_table(tbl, "units", context)
    x_units = _require_str(utbl, "x", "$context.units")
    y_units = Dict{String, String}()
    for event in events
        y_units[event] = _require_str(utbl, event, "$context.units")
    end
    return x_units, y_units
end

function _parse_tolerances(tbl, key::String, events::Vector{String}, context::String)::Dict{String, EventTolerance}
    ttbl = _require_table(tbl, key, context)
    out = Dict{String, EventTolerance}()
    for event in events
        etbl = _require_table(ttbl, event, "$context.$key")
        out[event] = (
            max_abs_km=_require_float(etbl, "max_abs_km", "$context.$key.$event"),
            max_nmae=_require_float(etbl, "max_nmae", "$context.$key.$event"),
            max_rmse_km=_optional_float(etbl, "max_rmse_km", Inf)
        )
    end
    return out
end

function _load_scenarios_from_manifest(manifest_path::String)::Vector{AbstractScenarioConfig}
    doc = TOML.parsefile(manifest_path)
    raw = _require_key(doc, "scenarios", "manifest")
    raw isa AbstractVector || throw(ArgumentError("Manifest field 'scenarios' must be an array of tables."))

    scenarios = AbstractScenarioConfig[]
    for (idx, entry) in pairs(raw)
        entry isa AbstractDict || throw(ArgumentError("manifest.scenarios[$idx] must be a table."))
        tbl = entry
        context = "manifest.scenarios[$idx]"

        name = _require_str(tbl, "name", context)
        kind = lowercase(_require_str(tbl, "kind", context))
        planet_name = lowercase(_require_str(tbl, "planet", context))
        events = _require_str_vector(tbl, "events", context)
        units_x, units_y = _parse_units(tbl, events, context)
        tolerances_quick = _parse_tolerances(tbl, "tolerances_quick", events, context)
        tolerances_full = _parse_tolerances(tbl, "tolerances_full", events, context)
        min_eval_points = _require_int(tbl, "min_eval_points", context)
        initial_time = _parse_initial_time(_require_table(tbl, "initial_time", context), "$context.initial_time")
        spacecraft = _parse_spacecraft_config(tbl, context)
        gravity_model = _parse_gravity_model(_require_str(tbl, "gravity_model", context), context)
        EI_km = _require_float(tbl, "EI_km", context)
        gravity_harmonics_degree = _optional_int(tbl, "gravity_harmonics_degree", 0)
        gravity_harmonics_order = _optional_int(tbl, "gravity_harmonics_order", 0)
        raw_harmonics_file = _optional_str(tbl, "gravity_harmonics_file", "")
        gravity_harmonics_file = isempty(strip(raw_harmonics_file)) ? "" : _resolve_repo_path(raw_harmonics_file)
        nbody_bodies = _optional_str_vector(tbl, "nbody_bodies")
        srp_enabled = _optional_bool(tbl, "srp_enabled", false)
        srp_cr = _optional_float(tbl, "srp_cr", 1.3)
        srp_area_m2 = _optional_float(tbl, "srp_area_m2", 0.0)
        drag_enabled = _optional_bool(tbl, "drag_enabled", true)
        include_wind = _optional_bool(tbl, "include_wind", false)
        orbit_altitude_mode = _parse_orbit_altitude_mode(_optional_str(tbl, "orbit_altitude_mode", "vacuum"), context)
        maneuver = _parse_maneuver_config(tbl, context)
        atmosphere_truth = _parse_atmosphere_truth_config(tbl, context)
        calibration = _parse_calibration_config(tbl, context)

        if kind == "orbit_events"
            push!(scenarios, OrbitEventsScenarioConfig(
                name=name,
                planet_name=planet_name,
                telemetry_peri_path=_resolve_repo_path(_require_str(tbl, "telemetry_peri", context)),
                telemetry_apo_path=_resolve_repo_path(_require_str(tbl, "telemetry_apo", context)),
                target_orbits_quick=_require_int(tbl, "target_orbits_quick", context),
                target_orbits_full=_require_int(tbl, "target_orbits_full", context),
                compare_points_quick=_require_int(tbl, "compare_points_quick", context),
                compare_points_full=_require_int(tbl, "compare_points_full", context),
                min_eval_points=min_eval_points,
                units_x=units_x,
                units_y=units_y,
                tolerances_quick=tolerances_quick,
                tolerances_full=tolerances_full,
                initial_time=initial_time,
                ra_m=_require_float(tbl, "ra_m", context),
                rp_altitude_m=_require_float(tbl, "rp_altitude_m", context),
                i_deg=_require_float(tbl, "i_deg", context),
                aop_deg=_require_float(tbl, "aop_deg", context),
                raan_deg=_require_float(tbl, "raan_deg", context),
                ta_deg=_require_float(tbl, "ta_deg", context),
                element_frame=_parse_element_frame(_optional_str(tbl, "element_frame", "j2000"), context),
                initial_state_j2000_m=_optional_float_tuple(tbl, "initial_state_j2000_m", 6, context),
                epoch_orbit_offset=haskey(tbl, "epoch_orbit_offset") ?
                    _require_float(tbl, "epoch_orbit_offset", context) : nothing,
                spacecraft=spacecraft,
                gravity_model=gravity_model,
                gravity_harmonics_degree=gravity_harmonics_degree,
                gravity_harmonics_order=gravity_harmonics_order,
                gravity_harmonics_file=gravity_harmonics_file,
                nbody_bodies=nbody_bodies,
                srp_enabled=srp_enabled,
                srp_cr=srp_cr,
                srp_area_m2=srp_area_m2,
                drag_enabled=drag_enabled,
                include_wind=include_wind,
                orbit_altitude_mode=orbit_altitude_mode,
                maneuver_orbit_numbers=maneuver.orbit_numbers,
                maneuver_orbit_numbers_campaign=maneuver.orbit_numbers_campaign,
                maneuver_delta_v_mps=maneuver.delta_v_mps,
                maneuver_replay_scale_mode=maneuver.replay_scale_mode,
                maneuver_flight_apoapsis_alt_m=maneuver.flight_apoapsis_alt_m,
                maneuver_thrust_n=maneuver.thrust_n,
                maneuver_isp_s=maneuver.isp_s,
                maneuver_guidance_rate_s=maneuver.guidance_rate_s,
                maneuver_control_rate_s=maneuver.control_rate_s,
                atmosphere_truth=atmosphere_truth,
                calibration=calibration,
                EI_km=EI_km
            ))
        elseif kind == "time_aligned_state"
            ctbl = _require_table(tbl, "telemetry_columns", context)
            comparison_mode = _parse_time_aligned_comparison_mode(
                _optional_str(tbl, "comparison_mode", "time_aligned_state"),
                context
            )
            extrema_min_separation_s = _optional_float(tbl, "extrema_min_separation_s", 500.0)
            extrema_min_separation_s > 0.0 || throw(ArgumentError(
                "extrema_min_separation_s must be > 0.0 in $context"
            ))
            push!(scenarios, TimeAlignedScenarioConfig(
                name=name,
                planet_name=planet_name,
                telemetry_path=_resolve_repo_path(_require_str(tbl, "telemetry", context)),
                telemetry_time_col=_require_str(ctbl, "time", "$context.telemetry_columns"),
                telemetry_altitude_col=_require_str(ctbl, "altitude", "$context.telemetry_columns"),
                telemetry_x_col=_require_str(ctbl, "x", "$context.telemetry_columns"),
                telemetry_y_col=_require_str(ctbl, "y", "$context.telemetry_columns"),
                telemetry_z_col=_require_str(ctbl, "z", "$context.telemetry_columns"),
                telemetry_sma_col=haskey(ctbl, "sma") ? String(ctbl["sma"]) : nothing,
                telemetry_ecc_col=haskey(ctbl, "ecc") ? String(ctbl["ecc"]) : nothing,
                telemetry_inc_col=haskey(ctbl, "inc") ? String(ctbl["inc"]) : nothing,
                telemetry_aop_col=haskey(ctbl, "aop") ? String(ctbl["aop"]) : nothing,
                telemetry_raan_col=haskey(ctbl, "raan") ? String(ctbl["raan"]) : nothing,
                telemetry_ta_col=haskey(ctbl, "ta") ? String(ctbl["ta"]) : nothing,
                telemetry_x_ic_col=haskey(ctbl, "x_ic") ? String(ctbl["x_ic"]) : nothing,
                telemetry_y_ic_col=haskey(ctbl, "y_ic") ? String(ctbl["y_ic"]) : nothing,
                telemetry_z_ic_col=haskey(ctbl, "z_ic") ? String(ctbl["z_ic"]) : nothing,
                telemetry_vx_ic_col=haskey(ctbl, "vx_ic") ? String(ctbl["vx_ic"]) : nothing,
                telemetry_vy_ic_col=haskey(ctbl, "vy_ic") ? String(ctbl["vy_ic"]) : nothing,
                telemetry_vz_ic_col=haskey(ctbl, "vz_ic") ? String(ctbl["vz_ic"]) : nothing,
                max_points_quick=_require_int(tbl, "max_points_quick", context),
                max_points_full=_require_int(tbl, "max_points_full", context),
                min_eval_points=min_eval_points,
                units_x=units_x,
                units_y=units_y,
                tolerances_quick=tolerances_quick,
                tolerances_full=tolerances_full,
                initial_time=initial_time,
                spacecraft=spacecraft,
                gravity_model=gravity_model,
                gravity_harmonics_degree=gravity_harmonics_degree,
                gravity_harmonics_order=gravity_harmonics_order,
                gravity_harmonics_file=gravity_harmonics_file,
                nbody_bodies=nbody_bodies,
                srp_enabled=srp_enabled,
                srp_cr=srp_cr,
                srp_area_m2=srp_area_m2,
                drag_enabled=drag_enabled,
                include_wind=include_wind,
                orbit_altitude_mode=orbit_altitude_mode,
                cartesian_ic_frame=_parse_reference_frame(
                    _optional_str(tbl, "cartesian_ic_frame", "inertial"),
                    "$context.cartesian_ic_frame"
                ),
                comparison_frame=_parse_reference_frame(
                    _optional_str(tbl, "comparison_frame", "inertial"),
                    "$context.comparison_frame"
                ),
                comparison_mode=comparison_mode,
                extrema_min_separation_s=extrema_min_separation_s,
                atmosphere_truth=atmosphere_truth,
                calibration=calibration,
                EI_km=EI_km
            ))
        else
            throw(ArgumentError("Unsupported scenario kind '$kind' in $context"))
        end
    end

    isempty(scenarios) && throw(ArgumentError("No scenarios defined in manifest: $manifest_path"))
    return scenarios
end

function parse_cli(args::Vector{String})::StudyConfig
    profile = Symbol(lowercase(get(ENV, "SPACEAGORA_TELEMETRY_PROFILE", "quick")))
    out_summary = get(ENV, "SPACEAGORA_TELEMETRY_OUT_SUMMARY", joinpath(DEFAULT_OUTPUT_DIR, "telemetry_orbit_accuracy_summary.csv"))
    out_errors = get(ENV, "SPACEAGORA_TELEMETRY_OUT_ERRORS", joinpath(DEFAULT_OUTPUT_DIR, "telemetry_orbit_accuracy_errors.csv"))
    manifest_path = get(ENV, "SPACEAGORA_TELEMETRY_MANIFEST", DEFAULT_MANIFEST_PATH)
    enforce = _safe_parse_bool(get(ENV, "SPACEAGORA_TELEMETRY_ENFORCE", "0"), false)
    generate_plots = _safe_parse_bool(get(ENV, "SPACEAGORA_TELEMETRY_PLOTS", "1"), true)

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("quick", "full")
            profile = Symbol(arg)
        elseif startswith(arg, "--profile=")
            profile = Symbol(lowercase(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--out-summary=")
            out_summary = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--out-errors=")
            out_errors = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--manifest=")
            manifest_path = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--enforce=")
            enforce = _safe_parse_bool(split(arg, "=", limit=2)[2], enforce)
        elseif startswith(arg, "--plots=")
            generate_plots = _safe_parse_bool(split(arg, "=", limit=2)[2], generate_plots)
        else
            throw(ArgumentError("Unsupported argument '$arg'. Supported: [quick|full], --profile=..., --manifest=..., --out-summary=..., --out-errors=..., --enforce=true|false, --plots=true|false"))
        end
        i += 1
    end

    profile in (:quick, :full) || throw(ArgumentError("Unsupported profile '$profile'. Use quick or full."))
    return StudyConfig(
        profile=profile,
        out_summary=abspath(out_summary),
        out_errors=abspath(out_errors),
        manifest_path=abspath(manifest_path),
        enforce=enforce,
        generate_plots=generate_plots
    )
end

@inline function _study_config(request::VerificationRequest)::StudyConfig
    profile = Symbol(lowercase(String(request.profile)))
    profile in (:quick, :full) || throw(ArgumentError("Unsupported request.profile '$profile'. Use quick or full."))
    return StudyConfig(
        profile=profile,
        out_summary=abspath(request.out_summary),
        out_errors=abspath(request.out_errors),
        manifest_path=abspath(request.manifest_path),
        enforce=request.enforce,
        generate_plots=request.generate_plots
    )
end

@inline function _request_from_study_config(cfg::StudyConfig)::VerificationRequest
    return VerificationRequest(
        profile=cfg.profile,
        out_summary=cfg.out_summary,
        out_errors=cfg.out_errors,
        manifest_path=cfg.manifest_path,
        enforce=cfg.enforce,
        generate_plots=cfg.generate_plots
    )
end
