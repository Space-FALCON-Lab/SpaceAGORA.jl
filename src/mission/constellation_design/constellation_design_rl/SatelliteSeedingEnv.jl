module SatelliteSeedingEnv

using ReinforcementLearningCore
using ..CostFunction
using ..DeepSetPolicy
using LinearAlgebra
using Statistics

"""
    ConstellationRLSatelliteSeedingObservation

Observation structure for the satellite seeding constellation RL environment.
"""
struct ConstellationRLSatelliteSeedingObservation
    # Variable-sized sets (DeepSet inputs)
    satellite_orbitals::Matrix{Float64}     # (n_sats, 6) [a, e, inc, raan, arg_p, ta]
    client_trajectories::Matrix{Float64}     # (n_clients, 6) [a, e, inc, raan, arg_p, ta]
    
    # Fixed-size context
    orbital_bounds::Matrix{Float64}          # (6, 2) search bounds [min, max] per element
    remaining_budget::Int                     # satellites remaining to place
    current_step::Int                        # current placement step
    max_sats::Int                            # total constellation capacity
    
    # Certificate state (from cost function)
    total_deficit::Float64
    unsafe_cost::Float64
    safe_cost::Float64  
    pred_cost::Float64
    feasible::Bool
end

"""
    SatelliteSeedingEnv <: AbstractEnv

RL environment for satellite constellation seeding.
"""
mutable struct SatelliteSeedingEnv <: AbstractEnv
    config_dict::Dict{String,Any}
    current_sats::Matrix{Float64}  # (n_placed, 6)
    client_orbitals::Matrix{Float64} # (n_clients, 6)
    step::Int
    max_sats::Int
    previous_cost::Float64
    action_space::Box{Float64}
    observation_space::CompositeObservation
    orbital_bounds::Matrix{Float64}
end

"""
    SatelliteSeedingEnv(config_dict::Dict{String,Any})

Construct a SatelliteSeedingEnv from configuration dictionary.
"""
function SatelliteSeedingEnv(config_dict::Dict{String,Any})
    # Extract parameters
    opt_params = get(config_dict, "optimizer_params", Dict{String,Any}())
    constellation_rl_config = get(opt_params, "constellation_rl_config", Dict{String,Any}())
    
    max_sats = Int(get(constellation_rl_config, "constellation_rl_max_sats", 64))
    
    # Get client data
    client_orbitals = get(config_dict, "client_orbitals", Matrix{Float64}(undef, 0, 6))
    orbital_bounds = get(config_dict, "orbital_bounds", zeros(6, 2))
    
    # Action space: 6D continuous normalized to [-1, 1]
    action_space = Box(-1.0, 1.0, (6,))
    
    # Observation space (composite)
    observation_space = CompositeObservation(
        satellite_orbitals = ArraySpec(Float64, (0, 6)),  # Variable size
        client_trajectories = ArraySpec(Float64, (size(client_orbitals, 1), 6)),
        orbital_bounds = ArraySpec(Float64, (6, 2)),
        remaining_budget = NumberSpec(Int),
        current_step = NumberSpec(Int),
        max_sats = NumberSpec(Int),
        total_deficit = NumberSpec(Float64),
        unsafe_cost = NumberSpec(Float64),
        safe_cost = NumberSpec(Float64),
        pred_cost = NumberSpec(Float64),
        feasible = NumberSpec(Bool),
    )
    
    env = SatelliteSeedingEnv(
        config_dict,
        zeros(0, 6),  # current_sats
        client_orbitals,
        0,  # step
        max_sats,
        Inf,  # previous_cost
        action_space,
        observation_space,
        orbital_bounds,
    )
    
    return env
end

"""
    RLBase.reset!(env::SatelliteSeedingEnv) -> RLSatelliteSeedingObservation

Reset the environment to initial state (empty constellation).
"""
function RLBase.reset!(env::SatelliteSeedingEnv)
    env.current_sats = zeros(0, 6)
    env.step = 1
    env.previous_cost = Inf
    
    # Compute initial cost (empty constellation)
    cost = compute_stage0_cost(env.config_dict, env.current_sats)
    
    return _build_observation(env, cost)
end

"""
    (env::SatelliteSeedingEnv)(action::Vector{Float64}) -> Tuple

Execute action: place a satellite and return (observation, reward, done, info).
"""
function (env::SatelliteSeedingEnv)(action::Vector{Float64})
    # Denormalize action from [-1, 1] to orbital bounds
    new_satellite = _denormalize_action(action, env)
    
    # Add satellite to constellation
    env.current_sats = vcat(env.current_sats, new_satellite')
    
    # Compute cost using cost function
    cost = compute_stage0_cost(env.config_dict, env.current_sats)
    
    # Compute reward
    reward = _compute_reward(env, cost)
    
    # Update previous cost
    env.previous_cost = cost["total_cost"]
    
    # Check termination
    done = env.step >= env.max_sats || cost["feasible"]
    
    # Build observation
    obs = _build_observation(env, cost)
    
    # Info dictionary
    info = Dict{String,Any}(
        "n_sats" => size(env.current_sats, 1),
        "total_cost" => cost["total_cost"],
        "feasible" => cost["feasible"],
    )
    
    env.step += 1
    
    return obs, reward, done, info
end

"""
    _denormalize_action(action::Vector{Float64}, env::SatelliteSeedingEnv) -> Vector{Float64}

Denormalize action from [-1, 1] to orbital element bounds.
"""
function _denormalize_action(action::Vector{Float64}, env::SatelliteSeedingEnv)
    bounds = env.orbital_bounds
    denormalized = zeros(6)
    
    for i in 1:6
        min_val = bounds[i, 1]
        max_val = bounds[i, 2]
        # Map from [-1, 1] to [min, max]
        denormalized[i] = min_val + (action[i] + 1) / 2 * (max_val - min_val)
    end
    
    # Ensure eccentricity is non-negative
    denormalized[2] = max(0.0, denormalized[2])
    
    # Normalize angles to [0, 2π]
    for i in [3, 4, 5, 6]  # inc, raan, arg_p, ta
        denormalized[i] = mod(denormalized[i], 2π)
    end
    
    return denormalized
end

"""
    _compute_reward(env::SatelliteSeedingEnv, cost::Dict) -> Float64

Compute reward based on cost improvement and penalties.
"""
function _compute_reward(env::SatelliteSeedingEnv, cost::Dict)
    # Get penalty weights
    opt_params = get(env.config_dict, "optimizer_params", Dict{String,Any}())
    constellation_rl_config = get(opt_params, "constellation_rl_config", Dict{String,Any}())
    
    unsafe_weight = Float64(get(constellation_rl_config, "constellation_rl_unsafe_weight", 100.0))
    safe_weight = Float64(get(constellation_rl_config, "constellation_rl_safe_weight", 10.0))
    pred_weight = Float64(get(constellation_rl_config, "constellation_rl_pred_weight", 5.0))
    
    # Immediate reward: cost improvement
    if isfinite(env.previous_cost)
        immediate_reward = env.previous_cost - cost["total_cost"]
    else
        immediate_reward = -cost["total_cost"]  # First step
    end
    
    # Penalty for constraint violations
    penalty = 0.0
    if cost["unsafe_cost"] > 1e-8
        penalty -= unsafe_weight * cost["unsafe_cost"]
    end
    if cost["safe_cost"] > 1e-8
        penalty -= safe_weight * cost["safe_cost"]
    end
    if cost["pred_cost"] > 1e-8
        penalty -= pred_weight * cost["pred_cost"]
    end
    
    # Episode bonus for feasibility
    episode_bonus = cost["feasible"] ? 1000.0 : 0.0
    
    return immediate_reward + penalty + episode_bonus
end

"""
    _build_observation(env::SatelliteSeedingEnv, cost::Dict) -> ConstellationRLSatelliteSeedingObservation

Build observation from current environment state.
"""
function _build_observation(env::SatelliteSeedingEnv, cost::Dict)
    return ConstellationRLSatelliteSeedingObservation(
        satellite_orbitals = copy(env.current_sats),
        client_trajectories = copy(env.client_orbitals),
        orbital_bounds = copy(env.orbital_bounds),
        remaining_budget = env.max_sats - size(env.current_sats, 1),
        current_step = env.step,
        max_sats = env.max_sats,
        total_deficit = cost["total_deficit"],
        unsafe_cost = cost["unsafe_cost"],
        safe_cost = cost["safe_cost"],
        pred_cost = cost["pred_cost"],
        feasible = cost["feasible"],
    )
end

export SatelliteSeedingEnv, ConstellationRLSatelliteSeedingObservation

end # module SatelliteSeedingEnv
