module Stage1Controllable

using JuMP
using MathOptInterface
using Clarabel
using HiGHS
using GLPK
using Ipopt
using LinearAlgebra
using ProgressMeter
using ..SupportFunctions
using ..ConstellationConfig
using ..ConstellationUtils

"""
    run_stage1_controllable_optimization(config_dict::AbstractDict) -> Dict{String,Any}

Run Stage 1 constellation optimization with horizon-indexed tube certificate constraints.

# Arguments
- `config_dict::AbstractDict`: Configuration dictionary

# Returns
- `Dict{String,Any}`: Stage 1 results including β, z_tube, and solver status
"""
function run_stage1_controllable_optimization(config_dict::AbstractDict)
    constellation_log("stage1", "Starting Stage 1 horizon-indexed tube certificate optimization")
    
    # Check for cached Stage 0 results
    opt_params = config_dict["optimizer_params"]
    if !haskey(config_dict, "stage0_h_fwd_exact_coeffs")
        constellation_log("stage1", "Stage 0 coefficients not found, running Stage 0")
        using ..Stage0FHSG
        stage0_data = run_fhsg_stage0(config_dict)
        config_dict["stage0_h_fwd_exact_coeffs"] = stage0_data["h_fwd_exact_coeffs"]
        config_dict["stage0_h_Wcorr_coeffs"] = stage0_data["h_Wcorr_coeffs"]
        config_dict["stage0_support_lift_coeffs"] = stage0_data["support_lift_coeffs"]
        config_dict["stage0_backward_lift_coeffs"] = stage0_data["backward_lift_coeffs"]
    end
    
    # Extract coefficients from Stage 0
    h_fwd_exact_coeffs = config_dict["stage0_h_fwd_exact_coeffs"]
    h_Wcorr_coeffs = config_dict["stage0_h_Wcorr_coeffs"]
    support_lift_coeffs = config_dict["stage0_support_lift_coeffs"]
    
    # Build and solve optimization model
    result = _build_and_solve_stage1_model(config_dict, h_fwd_exact_coeffs, h_Wcorr_coeffs, 
                                           support_lift_coeffs)
    
    constellation_log("stage1", "Stage 1 optimization complete")
    
    return result
end

"""
    _build_and_solve_stage1_model(config_dict::AbstractDict, h_fwd_exact_coeffs::Array{<:Real,3}, 
                                   h_Wcorr_coeffs::Array{<:Real,3}, support_lift_coeffs::Array{<:Real,4}) -> Dict{String,Any}

Build and solve the Stage 1 optimization model with horizon-indexed tube certificate.

# Arguments
- `config_dict::AbstractDict`: Configuration dictionary
- `h_fwd_exact_coeffs::Array{<:Real,3}`: Forward support coefficients [Kd, M, P]
- `h_Wcorr_coeffs::Array{<:Real,3}`: Correction authority coefficients [Kd, M, P]
- `support_lift_coeffs::Array{<:Real,4}`: SL coefficients [Kd, Kd, P, H]

# Returns
- `Dict{String,Any}`: Optimization results
"""
function _build_and_solve_stage1_model(config_dict::AbstractDict, h_fwd_exact_coeffs::Array{<:Real,3}, 
                                       h_Wcorr_coeffs::Array{<:Real,3}, support_lift_coeffs::Array{<:Real,4})
    opt_params = config_dict["optimizer_params"]
    effector_params = config_dict["effector_params"]
    mission = config_dict["mission"]
    
    # Dimensions
    H = Int(mission["n_horizons"])
    Kd = Int(opt_params["controllable_set_dirs"])
    M = size(h_fwd_exact_coeffs, 2)
    P = size(h_fwd_exact_coeffs, 3)
    
    # Direction bank
    keepout_dirs = _build_keepout_position_dirs(Kd)
    dirs_mat = _build_polyhedral_dirs(Kd)
    
    # Parameters
    d_safe = effector_params["range"]
    use_predecessor = true
    enforced_horizons = H
    
    # LP warm start enabled?
    lp_warm_start = Bool(get(opt_params, "stage1_lp_warm_start_enabled", true))
    
    # Select solver
    solver_name = lowercase(opt_params["solver"])
    optimizer = _select_optimizer(solver_name)
    
    # Create model
    model = Model(optimizer)
    set_silent(model)
    
    # Variables
    @variable(model, 0 <= β[1:M] <= 1)
    @variable(model, z_tube[1:H, 1:P, 1:Kd] >= 0.0)
    
    # Add constraints
    _add_chained_predecessor_constraints!(model, β, z_tube, h_fwd_exact_coeffs, h_Wcorr_coeffs,
                                          support_lift_coeffs, dirs_mat, d_safe, H, P, Kd, M,
                                          use_predecessor, enforced_horizons)
    
    # Objective (convex - no β(1-β) term per LADS invariants)
    γ1 = Float64(opt_params["num_sats_weight"])
    @objective(model, Min, γ1 * sum(β[i] for i in 1:M))
    
    # Solve
    constellation_log("stage1", "Solving Stage 1 optimization"; solver=solver_name, lp_warm_start=lp_warm_start)
    optimize!(model)
    
    # Extract results
    term_status = termination_status(model)
    has_sol = has_values(model)
    β_val = has_sol ? Float64.(value.(β)) : zeros(Float64, M)
    z_val = has_sol ? Float64.(value.(z_tube)) : zeros(Float64, H, P, Kd)
    objective_val = has_sol ? try objective_value(model) catch; Inf end : Inf
    solve_time = has_sol ? solve_time(model) : NaN
    
    # Count active satellites
    beta_threshold = Float64(opt_params["beta_threshold"])
    num_active = count(>(beta_threshold), β_val)
    
    # Compute h_C_var_history for verification (index 0 = zeros for C_0={0})
    h_C_var_history = _compute_h_C_var_history(β_val, h_fwd_exact_coeffs, support_lift_coeffs, H, P, Kd, M)
    
    return Dict{String,Any}(
        "beta" => β_val,
        "z_tube" => z_val,
        "h_C_var_history" => h_C_var_history,
        "objective_value" => objective_val,
        "term_status" => string(term_status),
        "num_active" => num_active,
        "solve_time" => solve_time,
        "solver" => solver_name,
        "feasible" => has_sol && term_status in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL),
    )
end

"""
    _add_chained_predecessor_constraints!(model::Model, β::Vector{VariableRef}, z_tube::Array{VariableRef,3},
                                          h_fwd_exact_coeffs::Array{<:Real,3}, h_Wcorr_coeffs::Array{<:Real,3},
                                          support_lift_coeffs::Array{<:Real,4}, dirs_mat::Matrix{Float64},
                                          d_safe::Float64, H::Int, P::Int, Kd::Int, M::Int,
                                          use_predecessor::Bool, enforced_horizons::Int)

Add chained predecessor constraints to the Stage 1 model with horizon-indexed tube certificate.

# Arguments
- `model::Model`: JuMP model
- `β::Vector{VariableRef}`: Satellite selection variables
- `z_tube::Array{VariableRef,3}`: Tube size variables [H, P, Kd]
- `h_fwd_exact_coeffs::Array{<:Real,3}`: Forward support coefficients [Kd, M, P]
- `h_Wcorr_coeffs::Array{<:Real,3}`: Correction authority coefficients [Kd, M, P]
- `support_lift_coeffs::Array{<:Real,4}`: SL coefficients [Kd, Kd, P, H]
- `dirs_mat::Matrix{Float64}`: Direction matrix [6, Kd]
- `d_safe::Float64`: Safe distance
- `H::Int`: Number of horizons
- `P::Int`: Number of clients
- `Kd::Int`: Number of directions
- `M::Int`: Number of satellites
- `use_predecessor::Bool`: Whether to use predecessor constraints
- `enforced_horizons::Int`: Number of horizons to enforce keepout
"""
function _add_chained_predecessor_constraints!(model::Model, β::Vector{VariableRef}, z_tube::Array{VariableRef,3},
                                              h_fwd_exact_coeffs::Array{<:Real,3}, h_Wcorr_coeffs::Array{<:Real,3},
                                              support_lift_coeffs::Array{<:Real,4}, dirs_mat::Matrix{Float64},
                                              d_safe::Float64, H::Int, P::Int, Kd::Int, M::Int,
                                              use_predecessor::Bool, enforced_horizons::Int)
    for n in 1:H
        sl_n = support_lift_coeffs[:, :, :, n]   # [Kd, Kd, P]

        for p in 1:P, q in 1:Kd
            # Forward containment: z[n,p,q] ≤ h_fwd(β)[q,p]
            scale_contain = max(maximum(abs.(view(h_fwd_exact_coeffs, q, :, p))), 1.0)
            h_fwd_qp = AffExpr(0.0)
            for i in 1:M
                c = h_fwd_exact_coeffs[q, i, p]
                add_to_expression!(h_fwd_qp, c / scale_contain, β[i])
            end
            @constraint(model, z_tube[n, p, q] / scale_contain <= h_fwd_qp)

            # Predecessor chain (n≥2): SL·z[n-1] ≤ z[n] + h_neg(β)
            if use_predecessor && n >= 2
                scale_pred = max(
                    maximum(abs.(view(sl_n, :, q, p))),
                    maximum(abs.(view(h_Wcorr_coeffs, q, :, p))),
                    1.0,
                )
                h_Wcorr_qp = AffExpr(0.0)
                for i in 1:M
                    c = h_Wcorr_coeffs[q, i, p]
                    add_to_expression!(h_Wcorr_qp, c / scale_pred, β[i])
                end
                sl_z_prev = AffExpr(0.0)
                for j in 1:Kd
                    c = sl_n[j, q, p]
                    add_to_expression!(sl_z_prev, c / scale_pred, z_tube[n - 1, p, j])
                end
                @constraint(model, sl_z_prev <= z_tube[n, p, q] / scale_pred + h_Wcorr_qp)
            end
        end
    end

    # All-direction keepout support floor: h_{B(d_safe)}(d_q) ≤ z[n,p,q] for all q with nonzero position component
    for n in 1:enforced_horizons, p in 1:P, q in 1:Kd
        d_pos_norm = norm(view(dirs_mat, 1:3, q))
        d_pos_norm > 1e-10 || continue
        keepout_required = d_safe * d_pos_norm
        scale_keepout = max(abs(keepout_required), 1.0)
        @constraint(model, z_tube[n, p, q] / scale_keepout >= keepout_required / scale_keepout)
    end
end

"""
    _select_optimizer(solver_name::AbstractString) -> Optimizer

Select the appropriate optimizer based on configuration.

# Arguments
- `solver_name::AbstractString`: Solver name

# Returns
- Optimizer instance
"""
function _select_optimizer(solver_name::AbstractString)
    if solver_name == "clarabel"
        return Clarabel.Optimizer
    elseif solver_name == "highs"
        return HiGHS.Optimizer
    elseif solver_name == "glpk"
        return GLPK.Optimizer
    elseif solver_name == "ipopt"
        return Ipopt.Optimizer
    else
        @warn "Unknown solver $solver_name, defaulting to Clarabel"
        return Clarabel.Optimizer
    end
end

"""
    _compute_h_C_var_history(β::Vector{Float64}, h_fwd_exact_coeffs::Array{<:Real,3},
                             support_lift_coeffs::Array{<:Real,4}, H::Int, P::Int, Kd::Int, M::Int) -> Array{Float64,3}

Compute h_C variable history for Stage 2 verification with horizon-indexed tube.

# Arguments
- `β::Vector{Float64}`: Satellite selection variables
- `h_fwd_exact_coeffs::Array{<:Real,3}`: Forward support coefficients [Kd, M, P]
- `support_lift_coeffs::Array{<:Real,4}`: SL coefficients [Kd, Kd, P, H]
- `H::Int`: Number of horizons
- `P::Int`: Number of clients
- `Kd::Int`: Number of directions
- `M::Int`: Number of satellites

# Returns
- `Array{Float64,3}`: h_C variable history [H+1, P, Kd] (index 0 = zeros for C_0={0})
"""
function _compute_h_C_var_history(β::Vector{Float64}, h_fwd_exact_coeffs::Array{<:Real,3},
                                 support_lift_coeffs::Array{<:Real,4}, H::Int, P::Int, Kd::Int, M::Int)
    h_C_var = zeros(Float64, H + 1, P, Kd)
    
    # Index 0 = zeros (C_0 = {0})
    h_C_var[1, :, :] .= 0.0
    
    # Compute forward-lifted support for each horizon
    for p in 1:P
        # Initial condition from H
        for q in 1:Kd
            h_C_var[H + 1, p, q] = sum(h_fwd_exact_coeffs[q, i, p] * β[i] for i in 1:M)
        end
        
        # Backward lift through SL-chain
        for n in H:-1:1
            sl_n = support_lift_coeffs[:, :, :, n]  # [Kd, Kd, P]
            for q in 1:Kd
                h_C_var[n, p, q] = sum(sl_n[j, q, p] * h_C_var[n + 1, p, j] for j in 1:Kd)
            end
        end
    end
    
    return h_C_var
end

export run_stage1_controllable_optimization

end # module Stage1Controllable
