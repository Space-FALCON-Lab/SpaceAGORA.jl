struct Transition
    observation::Vector{Float32}
    action_index::Int
    reward::Float32
    next_observation::Vector{Float32}
    terminated::Bool
    truncated::Bool
    info_index::Int
end

Base.@kwdef mutable struct EpisodeSummary
    episode_index::Int = 0
    worker_id::Int = 1
    seed::Int = 0
    pass_count::Int = 0
    protected_passes::Int = 0
    episode_reward::Float64 = 0.0
    success::Bool = false
    impact::Bool = false
    out_of_drag_passage::Bool = false
    thermal_violations::Int = 0
    target_error_m::Float64 = NaN
    mission_duration_days::Float64 = 0.0
    total_delta_v_mps::Float64 = 0.0
    abm_delta_v_mps::Float64 = 0.0
    apoapsis_correction_delta_v_mps::Float64 = 0.0
    periapsis_raise_delta_v_mps::Float64 = 0.0
    total_mission_delta_v_mps::Float64 = 0.0
    maneuver_count::Int = 0
    solver_failures::Int = 0
    heat_rate_trace::Vector{Float64} = Float64[]
    action_trace::Vector{Float64} = Float64[]
    apoapsis_trace_m::Vector{Float64} = Float64[]
    periapsis_trace_m::Vector{Float64} = Float64[]
    omega_trace_rad::Vector{Float64} = Float64[]
    raan_trace_rad::Vector{Float64} = Float64[]
    inclination_trace_rad::Vector{Float64} = Float64[]
    reward_trace::Vector{Float64} = Float64[]
    protected_trace::Vector{Bool} = Bool[]
end

empty_episode_summary(; episode_index::Int=0, worker_id::Int=1, seed::Int=0) =
    EpisodeSummary(episode_index=episode_index, worker_id=worker_id, seed=seed)
