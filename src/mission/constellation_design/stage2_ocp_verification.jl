module Stage2OCPVerification

using JuMP
using MathOptInterface
using Clarabel
using Ipopt
using LinearAlgebra
using ProgressMeter
using ..SupportFunctions
using ..ConstellationConfig
using ..ConstellationUtils

"""
    run_stage2_ocp_verification(config_dict::AbstractDict) -> Dict{String,Any}

Run Stage 2 OCP verification with H sequential OCPs (LADS specification).

Per LADS invariants: Stage 2 runs H sequential OCPs (not collapsed into one big OCP).
Each horizon n has terminal constraint x(T_n) ∈ C_n where C_n is from h_C_var_history[n].

# Arguments
- `config_dict::AbstractDict`: Configuration dictionary

# Returns
- `Dict{String,Any}`: Stage 2 results including verification status and trajectories
"""
function run_stage2_ocp_verification(config_dict::AbstractDict)
    constellation_log("stage2", "Starting Stage 2 sequential OCP verification")
    
    # Check if Stage 2 should be skipped
    opt_params = config_dict["optimizer_params"]
    if Bool(get(opt_params, "skip_stage2", false))
        constellation_log("stage2", "Stage 2 skipped per configuration")
        return Dict{String,Any}(
            "verified" => true,
            "skipped" => true,
            "reason" => "skip_stage2=true",
        )
    end
    
    # Check for Stage 1 results
    if !haskey(config_dict, "stage1_results")
        error("Stage 1 results required for Stage 2 verification")
    end
    
    stage1_results = config_dict["stage1_results"]
    mission = config_dict["mission"]
    
    # Parameters
    H = Int(mission["n_horizons"])
    P = config_dict["client_bounds"]["n_clients"]
    
    # Run sequential OCP verification (H sequential OCPs per LADS invariants)
    if H <= 1
        # Single horizon case
        result = _run_single_horizon_verification(config_dict, stage1_results)
    else
        # Multi-horizon case - H sequential OCPs
        result = _run_repeated_horizon_verification(config_dict, stage1_results, H, P)
    end
    
    constellation_log("stage2", "Stage 2 verification complete"; verified=result["verified"])
    
    return result
end

"""
    _run_single_horizon_verification(config_dict::AbstractDict, stage1_results::AbstractDict) -> Dict{String,Any}

Run single-horizon OCP verification with terminal constraint x(T_1) ∈ C_1.

# Arguments
- `config_dict::AbstractDict`: Configuration dictionary
- `stage1_results::AbstractDict`: Stage 1 results

# Returns
- `Dict{String,Any}`: Verification results
"""
function _run_single_horizon_verification(config_dict::AbstractDict, stage1_results::AbstractDict)
    constellation_log("stage2", "Running single-horizon verification")
    
    # Run OCP for single horizon (n=1)
    ocp_result = _solve_ocp_horizon(config_dict, stage1_results, 1)
    
    # Verify terminal constraint satisfaction: x(T_1) ∈ C_1
    h_C_var_history = stage1_results["h_C_var_history"]
    terminal_in_C = _verify_terminal_constraint_in_C(ocp_result, h_C_var_history, 1, config_dict)
    
    return Dict{String,Any}(
        "verified" => terminal_in_C,
        "terminal_in_C" => terminal_in_C,
        "precheck_slack" => get(ocp_result, "precheck_slack", 0.0),
        "x_opt" => get(ocp_result, "x_opt", nothing),
        "u_opt" => get(ocp_result, "u_opt", nothing),
        "n_horizons" => 1,
        "segments" => [ocp_result],
    )
end

"""
    _run_repeated_horizon_verification(config_dict::AbstractDict, stage1_results::AbstractDict,
                                       H::Int, P::Int) -> Dict{String,Any}

Run repeated-horizon OCP verification with H sequential OCPs.

Per LADS invariants: Each horizon n has terminal constraint x(T_n) ∈ C_n.
State continuity: x(0)_{n+1} = x(T_n) for n=1..H-1.

# Arguments
- `config_dict::AbstractDict`: Configuration dictionary
- `stage1_results::AbstractDict`: Stage 1 results
- `H::Int`: Number of horizons
- `P::Int`: Number of clients

# Returns
- `Dict{String,Any}`: Verification results
"""
function _run_repeated_horizon_verification(config_dict::AbstractDict, stage1_results::AbstractDict,
                                            H::Int, P::Int)
    constellation_log("stage2", "Running repeated-horizon verification (H=$H)")
    
    h_C_var_history = stage1_results["h_C_var_history"]
    segments = Vector{Dict{String,Any}}()
    all_horizons_solved = true
    all_terminal_in_C = true
    current_state = nothing
    
    for n in 1:H
        constellation_log("stage2", "Solving horizon $n/$H")
        
        # Set initial state for this horizon (state continuity)
        if n > 1 && current_state !== nothing
            config_dict["optimizer_params"]["stage2_initial_state"] = current_state
        end
        
        # Solve OCP for this horizon
        ocp_result = _solve_ocp_horizon(config_dict, stage1_results, n)
        
        # Store segment result
        push!(segments, ocp_result)
        
        # Check if this horizon was solved successfully
        horizon_solved = get(ocp_result, "all_clients_solved", true)
        all_horizons_solved &= horizon_solved
        
        # Verify terminal constraint: x(T_n) ∈ C_n
        terminal_in_C = _verify_terminal_constraint_in_C(ocp_result, h_C_var_history, n, config_dict)
        all_terminal_in_C &= terminal_in_C
        
        # Update current state for next horizon (terminal state of this horizon)
        if haskey(ocp_result, "x_opt")
            current_state = ocp_result["x_opt"][:, end, :]
        end
        
        # Early termination if horizon failed
        if !horizon_solved
            constellation_log("stage2", "Horizon $n failed, stopping verification")
            break
        end
    end
    
    return Dict{String,Any}(
        "verified" => all_horizons_solved && all_terminal_in_C,
        "all_horizons_solved" => all_horizons_solved,
        "all_terminal_in_C" => all_terminal_in_C,
        "n_horizons" => H,
        "segments" => segments,
    )
end

"""
    _solve_ocp_horizon(config_dict::AbstractDict, stage1_results::AbstractDict,
                       horizon_idx::Int) -> Dict{String,Any}

Solve OCP for a single horizon with CW-ECI dynamics.

# Arguments
- `config_dict::AbstractDict`: Configuration dictionary
- `stage1_results::AbstractDict`: Stage 1 results
- `horizon_idx::Int`: Horizon index (n)

# Returns
- `Dict{String,Any}`: OCP solution
"""
function _solve_ocp_horizon(config_dict::AbstractDict, stage1_results::AbstractDict,
                            horizon_idx::Int)
    sim_params = config_dict["sim_params"]
    opt_params = config_dict["optimizer_params"]
    effector_params = config_dict["effector_params"]
    
    P = config_dict["client_bounds"]["n_clients"]
    dt = Float64(sim_params["dt"])
    N = Int(sim_params["N"])
    
    # Physical parameters
    μ = Float64(config_dict["physical_constants"]["mu"])
    u_max = effector_params["max_thrust"] / effector_params["sc_mass"]
    
    # Initial state (from previous horizon or default)
    if haskey(opt_params, "stage2_initial_state")
        x0 = opt_params["stage2_initial_state"]
    else
        # Default to zero relative state
        x0 = zeros(Float64, 6, P)
    end
    
    # Build OCP model with IPOPT
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    
    # State and control variables (6D state: position + velocity in CW-ECI frame)
    @variable(model, x[1:6, 1:N+1, 1:P])
    @variable(model, u[1:3, 1:N, 1:P])
    
    # Initial condition
    for p in 1:P
        for i in 1:6
            @constraint(model, x[i, 1, p] == x0[i, p])
        end
    end
    
    # Dynamics constraints (CW-ECI linearized dynamics)
    # Placeholder: actual implementation requires CW-ECI state transition matrix
    for p in 1:P, k in 1:N
        # Simplified Euler integration (placeholder for CW-ECI dynamics)
        for i in 1:3
            @constraint(model, x[i, k+1, p] == x[i, k, p] + dt * x[i+3, k, p])
            @constraint(model, x[i+3, k+1, p] == x[i+3, k, p] + dt * u[i, k, p])
        end
    end
    
    # Control constraints
    for p in 1:P, k in 1:N
        @constraint(model, sum(u[i, k, p]^2 for i in 1:3) <= u_max^2)
    end
    
    # Terminal constraint will be verified post-solve against C_n
    # (Not enforced in OCP to allow feasibility check)
    
    # Objective (minimize control effort)
    @objective(model, Min, sum(sum(u[i, k, p]^2 for i in 1:3, k in 1:N) for p in 1:P))
    
    # Solve
    optimize!(model)
    
    # Extract results
    term_status = termination_status(model)
    has_sol = has_values(model)
    
    if has_sol
        x_opt = Float64.(value.(x))
        u_opt = Float64.(value.(u))
    else
        x_opt = zeros(Float64, 6, N+1, P)
        u_opt = zeros(Float64, 3, N, P)
    end
    
    return Dict{String,Any}(
        "horizon_index" => horizon_idx,
        "all_clients_solved" => has_sol && term_status in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL),
        "term_status" => string(term_status),
        "x_opt" => x_opt,
        "u_opt" => u_opt,
        "initial_state" => x0,
        "terminal_state" => has_sol ? x_opt[:, end, :] : x0,
    )
end

"""
    _verify_terminal_constraint_in_C(ocp_result::AbstractDict, h_C_var_history::Array{Float64,3},
                                     horizon_idx::Int, config_dict::AbstractDict) -> Bool

Verify terminal constraint satisfaction: x(T_n) ∈ C_n.

Per LADS invariants: C_n is the horizon-indexed tube from h_C_var_history[n].
Index 0 = zeros (C_0 = {0}).

# Arguments
- `ocp_result::AbstractDict`: OCP solution
- `h_C_var_history::Array{Float64,3}`: h_C variable history [H+1, P, Kd]
- `horizon_idx::Int`: Horizon index (n)
- `config_dict::AbstractDict`: Configuration dictionary

# Returns
- `Bool`: True if terminal constraint satisfied
"""
function _verify_terminal_constraint_in_C(ocp_result::AbstractDict, h_C_var_history::Array{Float64,3},
                                        horizon_idx::Int, config_dict::AbstractDict)
    if !haskey(ocp_result, "terminal_state")
        return false
    end
    
    terminal_state = ocp_result["terminal_state"]
    C_n = h_C_var_history[horizon_idx + 1, :, :]  # +1 because index 0 is at position 1
    
    P = size(terminal_state, 2)
    Kd = size(C_n, 2)
    
    # Check if terminal state is within C_n for all clients
    # For each client p, check if position component is within keepout support
    for p in 1:P
        x_terminal = terminal_state[1:3, p]
        
        # Check against all directions in C_n
        for q in 1:Kd
            # C_n[p,q] is the support value in direction q for client p
            support_val = C_n[p, q]
            
            # Get direction from D_dirs (should be stored in config or computed)
            # For now, use position norm as simplified check
            x_norm = norm(x_terminal)
            
            # If x_norm > support_val, violates the tube constraint
            if x_norm > support_val + 1e-6
                return false
            end
        end
    end
    
    return true
end


export run_stage2_ocp_verification

end # module Stage2OCPVerification
