# Read settings only when called; importing reports does not alter the environment.
@inline function _parse_positive_int_env(name::String, default::Int)::Int
    raw = strip(get(ENV, name, string(default)))
    value = try
        parse(Int, raw)
    catch
        throw(ArgumentError("$name must be an integer, got '$raw'"))
    end
    return max(1, value)
end

@inline function _parse_nonnegative_int_env(name::String, default::Int)::Int
    raw = strip(get(ENV, name, string(default)))
    value = try
        parse(Int, raw)
    catch
        throw(ArgumentError("$name must be an integer, got '$raw'"))
    end
    return max(0, value)
end

@inline function _parse_nonnegative_float_env(name::String, default::Float64)::Float64
    raw = strip(get(ENV, name, string(default)))
    value = try
        parse(Float64, raw)
    catch
        throw(ArgumentError("$name must be a number, got '$raw'"))
    end
    return max(0.0, value)
end

@inline function _parse_bool_env(name::String, default::Bool)::Bool
    raw = lowercase(strip(get(ENV, name, default ? "1" : "0")))
    if raw in ("1", "true", "yes", "on")
        return true
    elseif raw in ("0", "false", "no", "off")
        return false
    end
    throw(ArgumentError("Invalid $name='$raw'. Use one of: 1/0, true/false, yes/no, on/off."))
end

@inline function _perf_smoke_mode()::Bool
    raw = lowercase(strip(get(ENV, "SPACEAGORA_PERF_SMOKE", "0")))
    return raw in ("1", "true", "yes", "on")
end

function _perf_case_name_filter()::Vector{String}
    raw = strip(get(ENV, "SPACEAGORA_PERF_CASES", ""))
    isempty(raw) && return String[]
    names = String[]
    seen = Set{String}()
    for part in split(raw, ",")
        name = strip(part)
        isempty(name) && continue
        if !(name in seen)
            push!(names, name)
            push!(seen, name)
        end
    end
    return names
end

@inline function _exclude_entry_scenarios()::Bool
    return _parse_bool_env("SPACEAGORA_PERF_EXCLUDE_ENTRY_SCENARIOS", false)
end

@inline function _split_rollout_enabled()::Bool
    return _parse_bool_env("SPACEAGORA_PERF_SPLIT_ROLLOUT_GATE", false)
end

@inline function _split_rollout_enforce()::Bool
    return _parse_bool_env("SPACEAGORA_PERF_SPLIT_ROLLOUT_ENFORCE", false)
end

@inline function _split_rollout_case_names()::Vector{String}
    raw = strip(get(
        ENV,
        "SPACEAGORA_PERF_SPLIT_ROLLOUT_CASES",
        "single_orientation_aero,proximity_2sat_orientation_fullstack_gnc_highrate,multi_4_gravity"
    ))
    tokens = String[]
    for token in split(raw, ",")
        t = strip(token)
        isempty(t) && continue
        push!(tokens, t)
    end
    return unique(tokens)
end

@inline function _split_rollout_solver_variants()::Vector{String}
    raw = lowercase(strip(get(ENV, "SPACEAGORA_PERF_SPLIT_ROLLOUT_SOLVERS", "kencarp4,kencarp47")))
    tokens = String[]
    for token in split(raw, ",")
        t = lowercase(strip(token))
        isempty(t) && continue
        if !(t in ("kencarp4", "kencarp47", "kencarp58"))
            throw(ArgumentError(
                "SPACEAGORA_PERF_SPLIT_ROLLOUT_SOLVERS token '$t' is unsupported. " *
                "Use comma-separated values from: kencarp4, kencarp47, kencarp58."
            ))
        end
        push!(tokens, t)
    end
    isempty(tokens) && return ["kencarp4"]
    return unique(tokens)
end

@inline function _split_rollout_max_slowdown_ratio()::Float64
    return max(eps(Float64), _parse_nonnegative_float_env("SPACEAGORA_PERF_SPLIT_MAX_SLOWDOWN_RATIO", 1.15))
end

@inline function _split_rollout_pos_rel_tol()::Float64
    return _parse_nonnegative_float_env("SPACEAGORA_PERF_SPLIT_POS_REL_TOL", 5e-4)
end

@inline function _split_rollout_vel_rel_tol()::Float64
    return _parse_nonnegative_float_env("SPACEAGORA_PERF_SPLIT_VEL_REL_TOL", 5e-4)
end

@inline function _split_rollout_q_angle_tol_rad()::Float64
    return _parse_nonnegative_float_env("SPACEAGORA_PERF_SPLIT_Q_ANGLE_TOL_RAD", 5e-4)
end

@inline function _split_rollout_omega_rel_tol()::Float64
    return _parse_nonnegative_float_env("SPACEAGORA_PERF_SPLIT_OMEGA_REL_TOL", 1e-3)
end

@inline function _split_rollout_sample_count()::Int
    return _parse_positive_int_env("SPACEAGORA_PERF_SPLIT_TRAJ_SAMPLES", 128)
end

@inline function _multirate_rollout_enabled()::Bool
    return _parse_bool_env("SPACEAGORA_PERF_MULTIRATE_ROLLOUT_GATE", false)
end

@inline function _multirate_rollout_enforce()::Bool
    return _parse_bool_env("SPACEAGORA_PERF_MULTIRATE_ROLLOUT_ENFORCE", false)
end

@inline function _multirate_rollout_case_names()::Vector{String}
    raw = strip(get(
        ENV,
        "SPACEAGORA_PERF_MULTIRATE_ROLLOUT_CASES",
        "single_orientation_aero,proximity_2sat_orientation_fullstack_gnc_highrate,multi_4_gravity"
    ))
    tokens = String[]
    for token in split(raw, ",")
        t = strip(token)
        isempty(t) && continue
        push!(tokens, t)
    end
    return unique(tokens)
end

@inline function _multirate_rollout_max_slowdown_ratio()::Float64
    return max(eps(Float64), _parse_nonnegative_float_env("SPACEAGORA_PERF_MULTIRATE_MAX_SLOWDOWN_RATIO", 1.25))
end

@inline function _multirate_rollout_pos_rel_tol()::Float64
    return _parse_nonnegative_float_env("SPACEAGORA_PERF_MULTIRATE_POS_REL_TOL", 7.5e-4)
end

@inline function _multirate_rollout_vel_rel_tol()::Float64
    return _parse_nonnegative_float_env("SPACEAGORA_PERF_MULTIRATE_VEL_REL_TOL", 7.5e-4)
end

@inline function _multirate_rollout_q_angle_tol_rad()::Float64
    return _parse_nonnegative_float_env("SPACEAGORA_PERF_MULTIRATE_Q_ANGLE_TOL_RAD", 7.5e-4)
end

@inline function _multirate_rollout_omega_rel_tol()::Float64
    return _parse_nonnegative_float_env("SPACEAGORA_PERF_MULTIRATE_OMEGA_REL_TOL", 1.5e-3)
end

@inline function _multirate_rollout_sample_count()::Int
    return _parse_positive_int_env("SPACEAGORA_PERF_MULTIRATE_TRAJ_SAMPLES", 128)
end

@inline function _multirate_env_float_or_missing(name::String)
    raw = strip(get(ENV, name, ""))
    isempty(raw) && return missing
    parsed = try
        parse(Float64, raw)
    catch
        return missing
    end
    return parsed
end

@inline function _multirate_env_int_or_missing(name::String)
    raw = strip(get(ENV, name, ""))
    isempty(raw) && return missing
    parsed = try
        parse(Int, raw)
    catch
        return missing
    end
    return parsed
end

@inline function _multirate_rollout_setting_snapshot()
    slow_solver = lowercase(strip(get(ENV, "SPACEAGORA_MULTIRATE_SLOW_SOLVER", "tsit5")))
    fast_solver = lowercase(strip(get(ENV, "SPACEAGORA_MULTIRATE_FAST_SOLVER", "auto_stiff")))
    slow_dt_s = _multirate_env_float_or_missing("SPACEAGORA_MULTIRATE_SLOW_DT_S")
    fast_substeps = _multirate_env_int_or_missing("SPACEAGORA_MULTIRATE_FAST_SUBSTEPS")
    return (
        slow_solver=slow_solver,
        fast_solver=fast_solver,
        slow_dt_s=slow_dt_s,
        fast_substeps=fast_substeps
    )
end

