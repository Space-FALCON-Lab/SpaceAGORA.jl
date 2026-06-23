module CapoIntegration

using ..SatelliteSeedingEnv
using ..DeepSetPolicy
using ..PPOTrainer
using ..Scenarios
using ..CostFunction
using JLD2
using Random
using Distributions

"""
    run_rl_stage0_seeding(config_dict::Dict{String,Any}) -> Dict{String,Any}

Run RL-based stage 0 seeding using a trained policy.
"""
function run_rl_stage0_seeding(config_dict::Dict{String,Any})
    # Get model path from config
    opt_params = get(config_dict, "optimizer_params", Dict{String,Any}())
    rl_config = get(opt_params, "rl_config", Dict{String,Any}())
    model_path = get(rl_config, "rl_model_path", "data/rl_models/latest_model.jld2")
    
    # Load trained policy
    policy, loaded_config, episode = load_checkpoint(model_path)
    
    # Create environment
    env = SatelliteSeedingEnv(config_dict)
    
    # Run inference
    constellation = run_policy_inference(policy, env)
    
    return Dict{String,Any}(
        "constellation_orbitals" => constellation,
        "n_sats" => size(constellation, 1),
        "model_path" => model_path,
        "training_episode" => episode,
    )
end

"""
    run_policy_inference(policy::DualDeepSetPolicy, env::SatelliteSeedingEnv) -> Matrix{Float64}

Run policy inference to place satellites.
Returns matrix of orbital elements (n_sats × 6).
"""
function run_policy_inference(policy::DualDeepSetPolicy, env::SatelliteSeedingEnv)
    # Reset environment
    obs = RLBase.reset!(env)
    
    constellation_orbitals = zeros(0, 6)
    
    # Place satellites until done
    while true
        # Get action from policy
        action = policy(obs.satellite_orbitals, obs.client_trajectories)
        
        # Step environment
        next_obs, reward, done, info = env(action)
        
        # Store satellite
        new_sat = env.current_sats[end:end, :]
        constellation_orbitals = vcat(constellation_orbitals, new_sat)
        
        obs = next_obs
        
        if done
            break
        end
    end
    
    return constellation_orbitals
end

"""
    load_trained_policy(model_path::String) -> DualDeepSetPolicy

Load a trained policy from disk.
"""
function load_trained_policy(model_path::String)
    policy, config, episode = load_checkpoint(model_path)
    return policy
end

"""
    run_stochastic_greedy_seeding(config_dict::Dict{String,Any}) -> Dict{String,Any}

Run stochastic greedy seeding (baseline for comparison).
Iteratively adds satellites that maximize cost improvement.
"""
function run_stochastic_greedy_seeding(config_dict::Dict{String,Any})
    # Get parameters
    opt_params = get(config_dict, "optimizer_params", Dict{String,Any}())
    rl_config = get(opt_params, "rl_config", Dict{String,Any}())
    max_sats = Int(get(rl_config, "max_sats", 64))
    greedy_variant = get(rl_config, "greedy_variant", "pure")  # "pure", "stochastic", "epsilon"
    top_k = Int(get(rl_config, "greedy_top_k", 5))
    epsilon = Float64(get(rl_config, "greedy_epsilon", 0.1))
    n_restarts = Int(get(rl_config, "greedy_restarts", 1))
    
    # Get orbital bounds
    orbital_bounds = get(config_dict, "orbital_bounds", zeros(6, 2))
    
    # Run multiple restarts and pick best
    best_constellation = nothing
    best_cost = Inf
    
    for restart in 1:n_restarts
        rng = MersenneTwister(42 + restart)
        constellation_orbitals = _run_greedy_algorithm(
            config_dict, orbital_bounds, max_sats, 
            greedy_variant, top_k, epsilon, rng
        )
        
        # Evaluate final constellation
        cost_result = compute_stage0_cost(config_dict, constellation_orbitals)
        
        if cost_result["total_cost"] < best_cost
            best_cost = cost_result["total_cost"]
            best_constellation = constellation_orbitals
        end
    end
    
    return Dict{String,Any}(
        "constellation_orbitals" => best_constellation,
        "n_sats" => size(best_constellation, 1),
        "method" => "stochastic_greedy",
        "final_cost" => best_cost,
    )
end

"""
    _run_greedy_algorithm(config_dict, orbital_bounds, max_sats, variant, top_k, epsilon, rng)

Internal function to run the greedy algorithm with specified variant.
"""
function _run_greedy_algorithm(
    config_dict::Dict{String,Any},
    orbital_bounds::Matrix{Float64},
    max_sats::Int,
    variant::String,
    top_k::Int,
    epsilon::Float64,
    rng::AbstractRNG
)
    constellation_orbitals = zeros(0, 6)
    
    for sat_idx in 1:max_sats
        # Sample candidate satellites
        n_candidates = min(100, 1000)  # Number of candidates to evaluate
        candidates = _sample_candidates(orbital_bounds, n_candidates, rng)
        
        # Evaluate each candidate
        best_candidate = nothing
        best_improvement = -Inf
        
        for candidate in eachrow(candidates)
            # Add candidate to constellation
            test_constellation = vcat(constellation_orbitals, candidate')
            
            # Compute cost
            cost_result = compute_stage0_cost(config_dict, test_constellation)
            
            # Compute improvement (negative cost)
            improvement = -cost_result["total_cost"]
            
            if improvement > best_improvement
                best_improvement = improvement
                best_candidate = copy(candidate)
            end
        end
        
        # Select based on variant
        if variant == "pure"
            # Always pick best
            selected = best_candidate
        elseif variant == "stochastic"
            # Pick from top-k randomly
            top_candidates = _get_top_candidates(
                orbital_bounds, config_dict, constellation_orbitals, 
                top_k, rng
            )
            selected = top_candidates[rand(rng, 1:length(top_candidates))]
        elseif variant == "epsilon"
            # Epsilon-greedy: explore with probability epsilon
            if rand(rng) < epsilon
                selected = _sample_candidates(orbital_bounds, 1, rng)[1, :]
            else
                selected = best_candidate
            end
        else
            selected = best_candidate
        end
        
        # Add selected satellite
        constellation_orbitals = vcat(constellation_orbitals, selected')
        
        # Check if feasible
        cost_result = compute_stage0_cost(config_dict, constellation_orbitals)
        if cost_result["feasible"]
            break
        end
    end
    
    return constellation_orbitals
end

"""
    _sample_candidates(orbital_bounds, n_candidates, rng) -> Matrix{Float64}

Sample candidate satellite orbital elements within bounds.
"""
function _sample_candidates(orbital_bounds::Matrix{Float64}, n_candidates::Int, rng::AbstractRNG)
    candidates = zeros(n_candidates, 6)
    
    for i in 1:n_candidates
        for j in 1:6
            min_val = orbital_bounds[j, 1]
            max_val = orbital_bounds[j, 2]
            candidates[i, j] = rand(rng, Uniform(min_val, max_val))
        end
        # Ensure eccentricity is non-negative
        candidates[i, 2] = max(0.0, candidates[i, 2])
        # Normalize angles
        for j in [3, 4, 5, 6]
            candidates[i, j] = mod(candidates[i, j], 2π)
        end
    end
    
    return candidates
end

"""
    _get_top_candidates(orbital_bounds, config_dict, constellation, top_k, rng)

Get top-k candidates by cost improvement (for stochastic variant).
"""
function _get_top_candidates(
    orbital_bounds::Matrix{Float64},
    config_dict::Dict{String,Any},
    constellation::Matrix{Float64},
    top_k::Int,
    rng::AbstractRNG
)
    n_candidates = min(100, 1000)
    candidates = _sample_candidates(orbital_bounds, n_candidates, rng)
    
    candidate_scores = []
    
    for candidate in eachrow(candidates)
        test_constellation = vcat(constellation, candidate')
        cost_result = compute_stage0_cost(config_dict, test_constellation)
        improvement = -cost_result["total_cost"]
        push!(candidate_scores, (copy(candidate), improvement))
    end
    
    # Sort by improvement and take top-k
    sort!(candidate_scores, by = x -> x[2], rev = true)
    top_candidates = [c[1] for c in candidate_scores[1:min(top_k, length(candidate_scores))]]
    
    return top_candidates
end

export run_rl_stage0_seeding, run_policy_inference, load_trained_policy, run_stochastic_greedy_seeding

end # module CapoIntegration
