@enum ParallelProfile begin
    R0
    R1_a
    R1_b
    R2
    R3
    R4
    R5
end

# Backward-compatible alias for historical profile naming.
const R4_full_auto = R5

Base.@kwdef struct ParallelProfileConfig
    profile::ParallelProfile
    label::String
    outer_backend::Symbol
    inner_adaptive::Bool
    outer_route_adaptive::Bool
    density_mode::String
    control_mode::String
    thermal_mode::String
    multibody_mode::String
    effector_mode::String
    inner_scheduler::String = "static"
    adaptive_window::Int = 8
    adaptive_control_tail_guard::Bool = false
    adaptive_measured_reward::Bool = false
    persistent_hints::Bool = false
    persistent_state_persist::Bool = false
end

@inline function parallel_profile_name(profile::ParallelProfile)::String
    if profile == R0
        return "R0"
    elseif profile == R1_a
        return "R1_a"
    elseif profile == R1_b
        return "R1_b"
    elseif profile == R2
        return "R2"
    elseif profile == R3
        return "R3"
    elseif profile == R4
        return "R4"
    end
    return "R5"
end

@inline function _normalize_profile_token(raw::AbstractString)::String
    token = lowercase(strip(String(raw)))
    token = replace(token, "-" => "_")
    return replace(token, " " => "")
end

function parse_parallel_profile(raw::AbstractString)::ParallelProfile
    token = _normalize_profile_token(raw)
    if token in ("r0", "serial", "r0_true_serial", "true_serial")
        return R0
    elseif token in ("r1_a", "r1a", "r1_a_outer_only", "outer_only", "threads", "outer_only_threads")
        return R1_a
    elseif token in ("r1_b", "r1b", "r1_b_outer_only_process", "outer_only_process", "process_outer")
        return R1_b
    elseif token in ("r2", "r2_inner_only", "inner_only")
        return R2
    elseif token in ("r3", "r3_outer_inner_static", "outer_inner_static", "process_static")
        return R3
    elseif token in ("r4", "r4_outer_inner_adaptive", "outer_inner_adaptive", "auto_adaptive")
        return R4
    elseif token in (
        "r5",
        "r5_full_auto",
        "r5fullauto",
        "r5_calibration_full_auto",
        "full_smart",
        # Legacy aliases:
        "r4_full_auto",
        "r4fullauto",
        "r4_calibration_full_auto"
    )
        return R5
    end
    throw(ArgumentError(
        "Unsupported parallel profile '$raw'. Use one of: R0, R1_a, R1_b, R2, R3, R4, R5."
    ))
end

parse_parallel_profile(raw::Symbol)::ParallelProfile = parse_parallel_profile(String(raw))
parse_parallel_profile(profile::ParallelProfile)::ParallelProfile = profile

function profile_config(profile_in)::ParallelProfileConfig
    profile = parse_parallel_profile(profile_in)
    if profile == R0
        return ParallelProfileConfig(
            profile=profile,
            label="r0_true_serial",
            outer_backend=:none,
            inner_adaptive=false,
            outer_route_adaptive=false,
            density_mode="off",
            control_mode="off",
            thermal_mode="off",
            multibody_mode="off",
            effector_mode="off"
        )
    elseif profile == R1_a
        return ParallelProfileConfig(
            profile=profile,
            label="r1_a_outer_only",
            outer_backend=:threads,
            inner_adaptive=false,
            outer_route_adaptive=false,
            density_mode="off",
            control_mode="off",
            thermal_mode="off",
            multibody_mode="off",
            effector_mode="off"
        )
    elseif profile == R1_b
        return ParallelProfileConfig(
            profile=profile,
            label="r1_b_outer_only_process",
            outer_backend=:process,
            inner_adaptive=false,
            outer_route_adaptive=false,
            density_mode="off",
            control_mode="off",
            thermal_mode="off",
            multibody_mode="off",
            effector_mode="off"
        )
    elseif profile == R2
        return ParallelProfileConfig(
            profile=profile,
            label="r2_inner_only",
            outer_backend=:none,
            inner_adaptive=false,
            outer_route_adaptive=false,
            density_mode="auto",
            control_mode="auto",
            thermal_mode="auto",
            multibody_mode="auto",
            effector_mode="auto"
        )
    elseif profile == R3
        return ParallelProfileConfig(
            profile=profile,
            label="r3_outer_inner_static",
            outer_backend=:auto,
            inner_adaptive=false,
            outer_route_adaptive=false,
            density_mode="auto",
            control_mode="auto",
            thermal_mode="auto",
            multibody_mode="auto",
            effector_mode="auto"
        )
    elseif profile == R4
        return ParallelProfileConfig(
            profile=profile,
            label="r4_outer_inner_adaptive",
            outer_backend=:auto,
            inner_adaptive=true,
            outer_route_adaptive=true,
            density_mode="auto",
            control_mode="auto",
            thermal_mode="auto",
            multibody_mode="auto",
            effector_mode="auto"
        )
    end
    return ParallelProfileConfig(
        profile=profile,
        label="r5",
        outer_backend=:auto,
        inner_adaptive=true,
        outer_route_adaptive=true,
        density_mode="auto",
        control_mode="auto",
        thermal_mode="on",
        multibody_mode="auto",
        effector_mode="auto",
        inner_scheduler="dynamic",
        adaptive_window=4,
        adaptive_control_tail_guard=true,
        adaptive_measured_reward=true,
        persistent_hints=true,
        persistent_state_persist=true
    )
end

