module CostFunction

using LinearAlgebra
using Statistics

"""
    compute_stage0_cost(config_dict::Dict, satellite_orbitals::Matrix{Float64}) -> Dict{String,Any}

Compute the stage 0 cost for a given constellation configuration.
Routes to simplified or full CAPO audit cost function based on config.

# Arguments
- `config_dict`: Configuration dictionary with client data and parameters
- `satellite_orbitals`: Matrix of satellite orbital elements (n_sats × 6)

# Returns
Dictionary with cost components:
- `total_cost`: Weighted sum of all cost components
- `unsafe_cost`: Cost from unsafe constraint violations
- `safe_cost`: Cost from safe constraint violations  
- `pred_cost`: Cost from prediction errors
- `total_deficit`: Sum of all deficits
- `feasible`: Boolean indicating if constellation is feasible
"""
function compute_stage0_cost(config_dict::Dict, satellite_orbitals::Matrix{Float64})
    # Check which cost function to use
    opt_params = get(config_dict, "optimizer_params", Dict{String,Any}())
    constellation_rl_config = get(opt_params, "constellation_rl_config", Dict{String,Any}())
    use_capo_audit = Bool(get(constellation_rl_config, "constellation_rl_use_capo_audit", false))
    
    if use_capo_audit
        return compute_stage0_cost_capo_audit(config_dict, satellite_orbitals)
    else
        return compute_stage0_cost_simplified(config_dict, satellite_orbitals)
    end
end

"""
    compute_stage0_cost_simplified(config_dict::Dict, satellite_orbitals::Matrix{Float64}) -> Dict{String,Any}

Compute the stage 0 cost using simplified orbital element distance proxy.
This is the default cost function for RL training when use_capo_audit is false.
"""
function compute_stage0_cost_simplified(config_dict::Dict, satellite_orbitals::Matrix{Float64})
    n_sats = size(satellite_orbitals, 1)
    
    if n_sats == 0
        return Dict{String,Any}(
            "total_cost" => Inf,
            "unsafe_cost" => Inf,
            "safe_cost" => Inf,
            "pred_cost" => Inf,
            "total_deficit" => Inf,
            "feasible" => false,
        )
    end
    
    # Get client data
    client_orbitals = get(config_dict, "client_orbitals", Matrix{Float64}(undef, 0, 6))
    n_clients = size(client_orbitals, 1)
    
    if n_clients == 0
        return Dict{String,Any}(
            "total_cost" => Inf,
            "unsafe_cost" => Inf,
            "safe_cost" => Inf,
            "pred_cost" => Inf,
            "total_deficit" => Inf,
            "feasible" => false,
        )
    end
    
    # Get cost weights from config
    opt_params = get(config_dict, "optimizer_params", Dict{String,Any}())
    constellation_rl_config = get(opt_params, "constellation_rl_config", Dict{String,Any}())
    
    unsafe_weight = Float64(get(constellation_rl_config, "constellation_rl_unsafe_weight", 100.0))
    safe_weight = Float64(get(constellation_rl_config, "constellation_rl_safe_weight", 10.0))
    pred_weight = Float64(get(constellation_rl_config, "constellation_rl_pred_weight", 5.0))
    
    # Compute access coverage (simplified - in full implementation would use propagation)
    # For now, use orbital element distance as proxy for access
    access_matrix = compute_access_matrix(satellite_orbitals, client_orbitals)
    
    # Compute deficits
    unsafe_residuals, safe_residuals, pred_residuals = compute_residuals(
        access_matrix, satellite_orbitals, client_orbitals
    )
    
    # Compute costs
    unsafe_cost = unsafe_weight * sum(abs2, unsafe_residuals)
    safe_cost = safe_weight * sum(abs2, safe_residuals)
    pred_cost = pred_weight * sum(abs2, pred_residuals)
    
    total_deficit = sum(unsafe_residuals) + sum(safe_residuals) + sum(pred_residuals)
    total_cost = unsafe_cost + safe_cost + pred_cost
    
    # Feasibility check
    feasible_threshold = Float64(get(constellation_rl_config, "constellation_rl_feasibility_threshold", 1e-6))
    feasible = total_deficit < feasible_threshold && n_sats > 0
    
    return Dict{String,Any}(
        "total_cost" => total_cost,
        "unsafe_cost" => unsafe_cost,
        "safe_cost" => safe_cost,
        "pred_cost" => pred_cost,
        "total_deficit" => total_deficit,
        "feasible" => feasible,
        "n_sats" => n_sats,
        "n_clients" => n_clients,
    )
end

"""
    compute_access_matrix(satellite_orbitals::Matrix{Float64}, client_orbitals::Matrix{Float64}) -> Matrix{Float64}

Compute access matrix between satellites and clients based on orbital element distance.
Simplified proxy for actual access computation.
"""
function compute_access_matrix(satellite_orbitals::Matrix{Float64}, client_orbitals::Matrix{Float64})
    n_sats = size(satellite_orbitals, 1)
    n_clients = size(client_orbitals, 1)
    
    access_matrix = zeros(n_sats, n_clients)
    
    for i in 1:n_sats
        for j in 1:n_clients
            # Compute orbital element distance (simplified)
            sat_orb = satellite_orbitals[i, :]
            cli_orb = client_orbitals[j, :]
            
            # Weighted distance in orbital element space
            a_dist = abs(sat_orb[1] - cli_orb[1]) / sat_orb[1]  # Normalized semi-major axis
            e_dist = abs(sat_orb[2] - cli_orb[2])
            inc_dist = abs(sat_orb[3] - cli_orb[3])
            raan_dist = min(abs(sat_orb[4] - cli_orb[4]), 2π - abs(sat_orb[4] - cli_orb[4]))
            
            # Combined distance metric
            distance = sqrt(a_dist^2 + e_dist^2 + inc_dist^2 + (raan_dist/π)^2)
            
            # Access probability (closer = higher access)
            access_matrix[i, j] = exp(-distance)
        end
    end
    
    return access_matrix
end

"""
    compute_residuals(access_matrix::Matrix{Float64}, satellite_orbitals::Matrix{Float64}, client_orbitals::Matrix{Float64}) -> Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}

Compute unsafe, safe, and prediction residuals for cost calculation.
Simplified implementation for stage 0 seeding.
"""
function compute_residuals(access_matrix::Matrix{Float64}, satellite_orbitals::Matrix{Float64}, client_orbitals::Matrix{Float64})
    n_sats, n_clients = size(access_matrix)
    
    # Unsafe residuals: clients with insufficient access
    access_threshold = 0.5
    unsafe_residuals = Float64[]
    
    for j in 1:n_clients
        max_access = maximum(access_matrix[:, j])
        if max_access < access_threshold
            push!(unsafe_residuals, access_threshold - max_access)
        else
            push!(unsafe_residuals, 0.0)
        end
    end
    
    # Safe residuals: coverage margin
    safe_residuals = Float64[]
    for j in 1:n_clients
        max_access = maximum(access_matrix[:, j])
        safe_margin = max_access - access_threshold
        push!(safe_residuals, max(0.0, -safe_margin))  # Penalty if below threshold
    end
    
    # Prediction residuals: orbital element prediction error
    pred_residuals = Float64[]
    for i in 1:n_sats
        for j in 1:n_clients
            # Simplified prediction error based on access
            pred_error = 1.0 - access_matrix[i, j]
            push!(pred_residuals, pred_error)
        end
    end
    
    return unsafe_residuals, safe_residuals, pred_residuals
end

"""
    compute_stage0_cost_capo_audit(config_dict::Dict, satellite_orbitals::Matrix{Float64}) -> Dict{String,Any}

Compute the stage 0 cost using full CAPO certificate audit.
This is the accurate cost function for RL training when use_capo_audit is true.
It integrates with FHSG's safe build audit and predecessor repair audit.

# Arguments
- `config_dict`: Configuration dictionary with client data and parameters
- `satellite_orbitals`: Matrix of satellite orbital elements (n_sats × 6)

# Returns
Dictionary with cost components matching the simplified version for compatibility.
"""
function compute_stage0_cost_capo_audit(config_dict::Dict, satellite_orbitals::Matrix{Float64})
    n_sats = size(satellite_orbitals, 1)
    
    if n_sats == 0
        return Dict{String,Any}(
            "total_cost" => Inf,
            "unsafe_cost" => Inf,
            "safe_cost" => Inf,
            "pred_cost" => Inf,
            "total_deficit" => Inf,
            "feasible" => false,
        )
    end
    
    # Get client data
    client_orbitals = get(config_dict, "client_orbitals", Matrix{Float64}(undef, 0, 6))
    n_clients = size(client_orbitals, 1)
    
    if n_clients == 0
        return Dict{String,Any}(
            "total_cost" => Inf,
            "unsafe_cost" => Inf,
            "safe_cost" => Inf,
            "pred_cost" => Inf,
            "total_deficit" => Inf,
            "feasible" => false,
        )
    end
    
    # Get cost weights from config
    opt_params = get(config_dict, "optimizer_params", Dict{String,Any}())
    constellation_rl_config = get(opt_params, "constellation_rl_config", Dict{String,Any}())
    
    unsafe_weight = Float64(get(constellation_rl_config, "constellation_rl_unsafe_weight", 100.0))
    safe_weight = Float64(get(constellation_rl_config, "constellation_rl_safe_weight", 10.0))
    pred_weight = Float64(get(constellation_rl_config, "constellation_rl_pred_weight", 5.0))
    
    # TODO: Integrate with FHSG's _fhsg_safe_build_audit and _fhsg_fixed_target_pred_audit
    # For now, use the simplified version as a placeholder
    # This will be replaced with actual CAPO audit calls in a follow-up commit
    
    # Compute access coverage (placeholder - will use CAPO audit)
    access_matrix = compute_access_matrix(satellite_orbitals, client_orbitals)
    
    # Compute deficits
    unsafe_residuals, safe_residuals, pred_residuals = compute_residuals(
        access_matrix, satellite_orbitals, client_orbitals
    )
    
    # Compute costs
    unsafe_cost = unsafe_weight * sum(abs2, unsafe_residuals)
    safe_cost = safe_weight * sum(abs2, safe_residuals)
    pred_cost = pred_weight * sum(abs2, pred_residuals)
    
    total_deficit = sum(unsafe_residuals) + sum(safe_residuals) + sum(pred_residuals)
    total_cost = unsafe_cost + safe_cost + pred_cost
    
    # Feasibility check
    feasible_threshold = Float64(get(constellation_rl_config, "constellation_rl_feasibility_threshold", 1e-6))
    feasible = total_deficit < feasible_threshold && n_sats > 0
    
    return Dict{String,Any}(
        "total_cost" => total_cost,
        "unsafe_cost" => unsafe_cost,
        "safe_cost" => safe_cost,
        "pred_cost" => pred_cost,
        "total_deficit" => total_deficit,
        "feasible" => feasible,
        "n_sats" => n_sats,
        "n_clients" => n_clients,
    )
end

export compute_stage0_cost

end # module CostFunction
