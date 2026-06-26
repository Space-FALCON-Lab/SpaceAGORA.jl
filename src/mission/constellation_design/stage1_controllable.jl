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
    """
    Stage 1 LADS tube optimization with complete CAPO formulation:
    1. LP warm-start (convex, HiGHS)
    2. Nonconvex optimization (Ipopt, optional β(1-β) penalty)
    3. Beta ladder multistart for robustness
    """
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
    
    # Nonconvex fractionality penalty
    use_nonconvex_fractionality_penalty = Bool(get(opt_params, "stage1_nonconvex_fractionality_penalty_enabled", false))
    stage1_solver = use_nonconvex_fractionality_penalty ? "ipopt" : lowercase(opt_params["solver"])
    
    # Beta ladder multistart configuration
    beta_ladder_enabled = Bool(get(opt_params, "stage1_beta_ladder_enabled", true))
    beta_ladder_start = Float64(get(opt_params, "stage1_beta_ladder_start", 0.01))
    beta_ladder_step = Float64(get(opt_params, "stage1_beta_ladder_step", 0.01))
    beta_ladder_stop = Float64(get(opt_params, "stage1_beta_ladder_stop", 1.0))
    
    # LP warm-start configuration
    lp_warm_start_enabled = Bool(get(opt_params, "stage1_lp_warm_start_enabled", true))
    lp_coeff_rel_drop = Float64(get(opt_params, "stage1_lp_warm_start_coeff_rel_drop", 1e-12))
    
    # Build keepout position directions for constraints
    using ..SupportFunctions
    keepout_dirs_pos = _build_keepout_position_dirs(Int(opt_params["controllable_keepout_dirs"]))
    
    # Compute Phi_n state transition matrices for horizon-indexed mapping
    using ..TwoBodyLinearization
    a_ref = config_dict["constellation_params"]["a_max"]
    μ = config_dict["physical_constants"]["mu"]
    N_ctrl = config_dict["sim_params"]["N"]
    dt = config_dict["sim_params"]["dt"]
    Phi_horizons = compute_horizon_phi_matrices(a_ref, μ, N_ctrl, dt, H, P)
    
    # Generate beta initialization ladder
    beta_init_levels = _generate_beta_ladder(beta_ladder_enabled, beta_ladder_start, 
                                            beta_ladder_step, beta_ladder_stop)
    
    # Step 1: LP warm-start (convex, HiGHS)
    lp_warm_start_result = nothing
    if lp_warm_start_enabled
        constellation_log("stage1", "Running LP warm-start with HiGHS")
        lp_warm_start_result = _solve_lp_warm_start(
            M, H, P, Kd, h_fwd_exact_coeffs, h_Wcorr_coeffs, 
            support_lift_coeffs, dirs_mat, d_safe, use_predecessor, 
            enforced_horizons, keepout_dirs_pos, Phi_horizons, 
            opt_params, lp_coeff_rel_drop
        )
    end
    
    # Step 2: Nonconvex optimization with multistart
    constellation_log("stage1", "Running nonconvex optimization with multistart"; 
                      solver=stage1_solver, multistart_count=length(beta_init_levels))
    
    best_result = nothing
    best_objective = Inf
    
    # Combine LP warm-start with beta ladder for multistart
    start_points = []
    start_labels = []
    
    if lp_warm_start_result !== nothing
        push!(start_points, lp_warm_start_result["beta"])
        push!(start_labels, "lp_warm_start")
    end
    
    for beta_init in beta_init_levels
        push!(start_points, fill(beta_init, M))
        push!(start_labels, "beta_ladder_$(beta_init)")
    end
    
    for (attempt_idx, beta_start) in enumerate(start_points)
        result = _solve_nonconvex_stage1(
            beta_start, M, H, P, Kd, h_fwd_exact_coeffs, h_Wcorr_coeffs,
            support_lift_coeffs, dirs_mat, d_safe, use_predecessor,
            enforced_horizons, keepout_dirs_pos, Phi_horizons,
            opt_params, use_nonconvex_fractionality_penalty, stage1_solver,
            label=start_labels[attempt_idx]
        )
        
        if result["feasible"] && result["objective_value"] < best_objective
            best_objective = result["objective_value"]
            best_result = result
            best_result["selected_start_label"] = start_labels[attempt_idx]
            best_result["selected_start_attempt"] = attempt_idx
        end
    end
    
    # Fallback to first attempt if none succeeded
    if best_result === nothing
        best_result = _solve_nonconvex_stage1(
            zeros(M), M, H, P, Kd, h_fwd_exact_coeffs, h_Wcorr_coeffs,
            support_lift_coeffs, dirs_mat, d_safe, use_predecessor,
            enforced_horizons, keepout_dirs_pos, Phi_horizons,
            opt_params, use_nonconvex_fractionality_penalty, stage1_solver,
            label="fallback"
        )
    end
    
    # Compute h_C_var_history for Stage 2
    h_C_var_history = _compute_h_C_var_history_from_z(best_result["z_tube"], H, P, Kd)
    
    return merge(best_result, Dict{String,Any}(
        "h_C_var_history" => h_C_var_history,
        "lp_warm_start_used" => lp_warm_start_result !== nothing,
        "lp_warm_start_objective" => lp_warm_start_result !== nothing ? lp_warm_start_result["objective_value"] : NaN,
        "multistart_attempts" => length(start_points),
        "selected_start_label" => get(best_result, "selected_start_label", "unknown"),
    ))
end

"""
    _add_chained_predecessor_constraints!(model::Model, β::Vector{VariableRef}, z_tube::Array{VariableRef,3},
                                          h_fwd_exact_coeffs::Array{<:Real,3}, h_Wcorr_coeffs::Array{<:Real,3},
                                          support_lift_coeffs::Array{<:Real,4}, dirs_mat::Matrix{Float64},
                                          d_safe::Float64, H::Int, P::Int, Kd::Int, M::Int,
                                          use_predecessor::Bool, enforced_horizons::Int, keepout_dirs_pos::Matrix{Float64},
                                          Phi_horizons::Vector{Matrix{Float64}})

Add chained predecessor constraints to the Stage 1 model with horizon-indexed tube certificate.
Uses Phi_n state space mapping (CAPO LADS tube formulation).

# Arguments
- `model::Model`: JuMP model
- `β::Vector{VariableRef}`: Satellite selection variables
- `z_tube::Array{VariableRef,3}`: Tube size variables [H, P, 2*Kd] (positive/negative splitting)
- `h_fwd_exact_coeffs::Array{<:Real,3}`: Forward support coefficients [Kd, M, P]
- `h_Wcorr_coeffs::Array{<:Real,3}`: Correction authority coefficients [Kd, M, P]
- `support_lift_coeffs::Array{<:Real,4}`: SL coefficients [Kd, Kd, P, H] (unused with Phi_n mapping)
- `dirs_mat::Matrix{Float64}`: Direction matrix [6, Kd]
- `d_safe::Float64`: Safe distance
- `H::Int`: Number of horizons
- `P::Int`: Number of clients
- `Kd::Int`: Number of directions
- `M::Int`: Number of satellites
- `use_predecessor::Bool`: Whether to use predecessor constraints
- `enforced_horizons::Int`: Number of horizons to enforce keepout
- `keepout_dirs_pos::Matrix{Float64}`: Keepout position directions [3, Kkeep]
- `Phi_horizons::Vector{Matrix{Float64}}: Horizon-indexed state transition matrices [6, 6, P] per horizon
"""
function _add_chained_predecessor_constraints!(model::Model, β::Vector{VariableRef}, z_tube::Array{VariableRef,3},
                                              h_fwd_exact_coeffs::Array{<:Real,3}, h_Wcorr_coeffs::Array{<:Real,3},
                                              support_lift_coeffs::Array{<:Real,4}, dirs_mat::Matrix{Float64},
                                              d_safe::Float64, H::Int, P::Int, Kd::Int, M::Int,
                                              use_predecessor::Bool, enforced_horizons::Int, keepout_dirs_pos::Matrix{Float64},
                                              Phi_horizons::Vector{Matrix{Float64}})
    for n in 1:H
        Phi_n = Phi_horizons[n]  # [6, 6, P]

        for p in 1:P, q in 1:Kd
            # Forward containment: z[n,p,2q-1] and z[n,p,2q] for positive/negative direction q
            scale_contain = max(maximum(abs.(view(h_fwd_exact_coeffs, q, :, p))), 1.0)
            h_fwd_qp = AffExpr(0.0)
            for i in 1:M
                c = h_fwd_exact_coeffs[q, i, p]
                add_to_expression!(h_fwd_qp, c / scale_contain, β[i])
            end
            @constraint(model, z_tube[n, p, 2*q-1] / scale_contain <= h_fwd_qp)
            @constraint(model, z_tube[n, p, 2*q] / scale_contain <= h_fwd_qp)

            # Predecessor chain with Phi_n mapping (CAPO LADS tube formulation)
            if use_predecessor && n >= 2
                scale_pred = max(maximum(abs.(view(h_Wcorr_coeffs, q, :, p))), 1.0)
                h_Wcorr_qp = AffExpr(0.0)
                for i in 1:M
                    c = h_Wcorr_coeffs[q, i, p]
                    add_to_expression!(h_Wcorr_qp, c / scale_pred, β[i])
                end
                
                # Phi_n mapping: mapped_dir = transpose(Phi_n) * dirs_mat[:, q]
                Phi_n_p = view(Phi_n, :, :, p)  # [6, 6]
                mapped_dir = transpose(Phi_n_p) * view(dirs_mat, :, q)  # [6]
                
                h_prev_forward = AffExpr(0.0)
                for s in 1:6
                    coeff_pos = max(mapped_dir[s], 0.0)
                    coeff_neg = max(-mapped_dir[s], 0.0)
                    if coeff_pos > 1e-12
                        add_to_expression!(h_prev_forward, coeff_pos / scale_pred, z_tube[n - 1, p, 2*s-1])
                    end
                    if coeff_neg > 1e-12
                        add_to_expression!(h_prev_forward, coeff_neg / scale_pred, z_tube[n - 1, p, 2*s])
                    end
                end
                
                @constraint(model, h_prev_forward <= z_tube[n, p, 2*q-1] / scale_pred + h_Wcorr_qp)
                @constraint(model, h_prev_forward <= z_tube[n, p, 2*q] / scale_pred + h_Wcorr_qp)
            end
        end
    end

    # Keepout constraints with positive/negative splitting (CAPO LADS tube formulation)
    # Uses 3D position directions embedded in 6D state space
    Kkeep = size(keepout_dirs_pos, 2)
    
    for n in 0:H, p in 1:P, q in 1:Kkeep
        @constraint(model, sum(
            max(keepout_dirs_pos[s, q], 0.0) * z_tube[n+1, p, 2*s-1] +
            max(-keepout_dirs_pos[s, q], 0.0) * z_tube[n+1, p, 2*s]
            for s in 1:3
        ) >= d_safe)
    end
end

"""
    _add_chained_predecessor_constraints!(model::Model, β::Vector{VariableRef}, z_tube::Array{VariableRef,3},
                                          h_fwd_exact_coeffs::Array{<:Real,3}, h_Wcorr_coeffs::Array{<:Real,3},
                                          support_lift_coeffs::Array{<:Real,4}, dirs_mat::Matrix{Float64},
                                          d_safe::Float64, H::Int, P::Int, Kd::Int, M::Int,
                                          use_predecessor::Bool, enforced_horizons::Int, keepout_dirs_pos::Matrix{Float64},
                                          Phi_horizons::Vector{Matrix{Float64}}; coeff_rel_drop::Float64=0.0)

Add chained predecessor constraints with optional coefficient dropping for numerical stability.
"""
function _add_chained_predecessor_constraints!(model::Model, β::Vector{VariableRef}, z_tube::Array{VariableRef,3},
                                              h_fwd_exact_coeffs::Array{<:Real,3}, h_Wcorr_coeffs::Array{<:Real,3},
                                              support_lift_coeffs::Array{<:Real,4}, dirs_mat::Matrix{Float64},
                                              d_safe::Float64, H::Int, P::Int, Kd::Int, M::Int,
                                              use_predecessor::Bool, enforced_horizons::Int, keepout_dirs_pos::Matrix{Float64},
                                              Phi_horizons::Vector{Matrix{Float64}}; coeff_rel_drop::Float64=0.0)
    for n in 1:H
        Phi_n = Phi_horizons[n]  # [6, 6, P]

        for p in 1:P, q in 1:Kd
            # Forward containment: z[n,p,2q-1] and z[n,p,2q] for positive/negative direction q
            scale_contain = max(maximum(abs.(view(h_fwd_exact_coeffs, q, :, p))), 1.0)
            contain_cutoff = coeff_rel_drop * scale_contain
            h_fwd_qp = AffExpr(0.0)
            for i in 1:M
                c = h_fwd_exact_coeffs[q, i, p]
                abs(c) > contain_cutoff || continue
                add_to_expression!(h_fwd_qp, c / scale_contain, β[i])
            end
            @constraint(model, z_tube[n, p, 2*q-1] / scale_contain <= h_fwd_qp)
            @constraint(model, z_tube[n, p, 2*q] / scale_contain <= h_fwd_qp)

            # Predecessor chain with Phi_n mapping (CAPO LADS tube formulation)
            if use_predecessor && n >= 2
                scale_pred = max(maximum(abs.(view(h_Wcorr_coeffs, q, :, p))), 1.0)
                pred_cutoff = coeff_rel_drop * scale_pred
                h_Wcorr_qp = AffExpr(0.0)
                for i in 1:M
                    c = h_Wcorr_coeffs[q, i, p]
                    abs(c) > pred_cutoff || continue
                    add_to_expression!(h_Wcorr_qp, c / scale_pred, β[i])
                end
                
                # Phi_n mapping: mapped_dir = transpose(Phi_n) * dirs_mat[:, q]
                Phi_n_p = view(Phi_n, :, :, p)  # [6, 6]
                mapped_dir = transpose(Phi_n_p) * view(dirs_mat, :, q)  # [6]
                
                h_prev_forward = AffExpr(0.0)
                for s in 1:6
                    coeff_pos = max(mapped_dir[s], 0.0)
                    coeff_neg = max(-mapped_dir[s], 0.0)
                    if coeff_pos > 1e-12
                        add_to_expression!(h_prev_forward, coeff_pos / scale_pred, z_tube[n - 1, p, 2*s-1])
                    end
                    if coeff_neg > 1e-12
                        add_to_expression!(h_prev_forward, coeff_neg / scale_pred, z_tube[n - 1, p, 2*s])
                    end
                end
                
                @constraint(model, h_prev_forward <= z_tube[n, p, 2*q-1] / scale_pred + h_Wcorr_qp)
                @constraint(model, h_prev_forward <= z_tube[n, p, 2*q] / scale_pred + h_Wcorr_qp)
            end
        end
    end

    # Keepout constraints with positive/negative splitting (CAPO LADS tube formulation)
    # Uses 3D position directions embedded in 6D state space
    Kkeep = size(keepout_dirs_pos, 2)
    
    for n in 0:H, p in 1:P, q in 1:Kkeep
        @constraint(model, sum(
            max(keepout_dirs_pos[s, q], 0.0) * z_tube[n+1, p, 2*s-1] +
            max(-keepout_dirs_pos[s, q], 0.0) * z_tube[n+1, p, 2*s]
            for s in 1:3
        ) >= d_safe)
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
    _solve_lp_warm_start(M::Int, H::Int, P::Int, Kd::Int,
                         h_fwd_exact_coeffs::Array{<:Real,3}, 
                         h_Wcorr_coeffs::Array{<:Real,3},
                         support_lift_coeffs::Array{<:Real,4},
                         dirs_mat::Matrix{Float64}, d_safe::Float64,
                         use_predecessor::Bool, enforced_horizons::Int,
                         keepout_dirs_pos::Matrix{Float64},
                         Phi_horizons::Vector{Matrix{Float64}},
                         opt_params::Dict{String,Any},
                         coeff_rel_drop::Float64) -> Union{Nothing, Dict{String,Any}}

Solve convex LP warm-start using HiGHS with coefficient dropping for numerical stability.
"""
function _solve_lp_warm_start(M::Int, H::Int, P::Int, Kd::Int,
                               h_fwd_exact_coeffs::Array{<:Real,3}, 
                               h_Wcorr_coeffs::Array{<:Real,3},
                               support_lift_coeffs::Array{<:Real,4},
                               dirs_mat::Matrix{Float64}, d_safe::Float64,
                               use_predecessor::Bool, enforced_horizons::Int,
                               keepout_dirs_pos::Matrix{Float64},
                               Phi_horizons::Vector{Matrix{Float64}},
                               opt_params::Dict{String,Any},
                               coeff_rel_drop::Float64)
    lp_model = Model(HiGHS.Optimizer)
    set_silent(lp_model)
    
    @variable(lp_model, 0 <= β_lp[1:M] <= 1)
    @variable(lp_model, z_tube_lp[1:H, 1:P, 1:2*Kd] >= 0.0)
    
    _add_chained_predecessor_constraints!(lp_model, β_lp, z_tube_lp, 
                                          h_fwd_exact_coeffs, h_Wcorr_coeffs,
                                          support_lift_coeffs, dirs_mat, d_safe, H, P, Kd, M,
                                          use_predecessor, enforced_horizons, 
                                          keepout_dirs_pos, Phi_horizons;
                                          coeff_rel_drop=coeff_rel_drop)
    
    γ1 = Float64(get(opt_params, "num_sats_weight", 1.0))
    @objective(lp_model, Min, γ1 * sum(β_lp[i] for i in 1:M))
    
    optimize!(lp_model)
    
    term_status = termination_status(lp_model)
    has_sol = has_values(lp_model)
    
    if has_sol && term_status in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL)
        β_val = Float64.(value.(β_lp))
        z_val = Float64.(value.(z_tube_lp))
        objective_val = objective_value(lp_model)
        
        constellation_log("stage1", "LP warm-start successful"; 
                          objective=objective_val, 
                          active_sats=count(>(0.1), β_val))
        
        return Dict{String,Any}(
            "beta" => β_val,
            "z_tube" => z_val,
            "objective_value" => objective_val,
            "term_status" => string(term_status),
        )
    else
        constellation_log_warn("stage1", "LP warm-start failed"; status=string(term_status))
        return nothing
    end
end

"""
    _solve_nonconvex_stage1(beta_start::Vector{Float64}, M::Int, H::Int, P::Int, Kd::Int,
                            h_fwd_exact_coeffs::Array{<:Real,3}, 
                            h_Wcorr_coeffs::Array{<:Real,3},
                            support_lift_coeffs::Array{<:Real,4},
                            dirs_mat::Matrix{Float64}, d_safe::Float64,
                            use_predecessor::Bool, enforced_horizons::Int,
                            keepout_dirs_pos::Matrix{Float64},
                            Phi_horizons::Vector{Matrix{Float64}},
                            opt_params::Dict{String,Any},
                            use_nonconvex_fractionality_penalty::Bool,
                            solver_name::String;
                            label::String="unknown") -> Dict{String,Any}

Solve nonconvex Stage 1 optimization with optional β(1-β) fractionality penalty.
"""
function _solve_nonconvex_stage1(beta_start::Vector{Float64}, M::Int, H::Int, P::Int, Kd::Int,
                                  h_fwd_exact_coeffs::Array{<:Real,3}, 
                                  h_Wcorr_coeffs::Array{<:Real,3},
                                  support_lift_coeffs::Array{<:Real,4},
                                  dirs_mat::Matrix{Float64}, d_safe::Float64,
                                  use_predecessor::Bool, enforced_horizons::Int,
                                  keepout_dirs_pos::Matrix{Float64},
                                  Phi_horizons::Vector{Matrix{Float64}},
                                  opt_params::Dict{String,Any},
                                  use_nonconvex_fractionality_penalty::Bool,
                                  solver_name::String;
                                  label::String="unknown")
    optimizer = _select_optimizer(solver_name)
    model = Model(optimizer)
    set_silent(model)
    
    # Configure solver-specific settings
    if solver_name == "ipopt"
        tol = Float64(get(opt_params, "tolerance", 1e-4))
        max_iter = Int(get(opt_params, "max_iter", 5000))
        set_optimizer_attribute(model, "tol", tol)
        set_optimizer_attribute(model, "max_iter", max_iter)
        set_optimizer_attribute(model, "print_level", 0)
        set_optimizer_attribute(model, "hessian_approximation", "limited-memory")
        set_optimizer_attribute(model, "warm_start_init_point", "yes")
    else
        tol = Float64(get(opt_params, "tolerance", 1e-4))
        max_iter = Int(get(opt_params, "max_iter", 200))
        set_optimizer_attribute(model, "verbose", false)
        set_optimizer_attribute(model, "max_iter", max_iter)
        set_optimizer_attribute(model, "tol_gap_abs", tol)
        set_optimizer_attribute(model, "tol_gap_rel", tol)
        set_optimizer_attribute(model, "tol_feas", tol)
    end
    
    # Variables
    @variable(model, 0 <= β[1:M] <= 1)
    @variable(model, z_tube[1:H, 1:P, 1:2*Kd] >= 0.0)
    
    # Set warm-start values
    for i in 1:M
        set_start_value(β[i], beta_start[i])
    end
    
    # Add constraints
    _add_chained_predecessor_constraints!(model, β, z_tube, h_fwd_exact_coeffs, h_Wcorr_coeffs,
                                          support_lift_coeffs, dirs_mat, d_safe, H, P, Kd, M,
                                          use_predecessor, enforced_horizons, keepout_dirs_pos, Phi_horizons)
    
    # Objective with optional nonconvex fractionality penalty
    γ1 = Float64(get(opt_params, "num_sats_weight", 1.0))
    γ2 = Float64(get(opt_params, "l1_reg_weight", 0.0))
    γ3 = Float64(get(opt_params, "tube_size_weight", 0.0))
    
    if use_nonconvex_fractionality_penalty
        @objective(model, Min,
            γ1 * sum(β[i] for i in 1:M) +
            γ2 * sum(β[i] * (1.0 - β[i]) for i in 1:M) -
            γ3 * sum(z_tube[n, p, q] for n in 1:H, p in 1:P, q in 1:2*Kd)
        )
    elseif γ3 > 0.0
        @objective(model, Min, γ1 * sum(β[i] for i in 1:M) - γ3 * sum(z_tube))
    else
        @objective(model, Min, γ1 * sum(β[i] for i in 1:M))
    end
    
    # Solve
    optimize!(model)
    
    # Extract results
    term_status = termination_status(model)
    has_sol = has_values(model)
    β_val = has_sol ? Float64.(value.(β)) : zeros(Float64, M)
    z_val = has_sol ? Float64.(value.(z_tube)) : zeros(Float64, H, P, 2*Kd)
    objective_val = has_sol ? try objective_value(model) catch; Inf end : Inf
    solve_time = has_sol ? solve_time(model) : NaN
    
    feasible = has_sol && term_status in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_OPTIMAL)
    
    beta_threshold = Float64(opt_params["beta_threshold"])
    num_active = count(>(beta_threshold), β_val)
    
    constellation_log("stage1", "Nonconvex solve complete"; 
                      label=label, status=string(term_status), 
                      feasible=feasible, objective=objective_val,
                      active_sats=num_active)
    
    return Dict{String,Any}(
        "beta" => β_val,
        "z_tube" => z_val,
        "objective_value" => objective_val,
        "term_status" => string(term_status),
        "num_active" => num_active,
        "solve_time" => solve_time,
        "solver" => solver_name,
        "feasible" => feasible,
        "beta_start_value" => beta_start[1],
    )
end

"""
    _generate_beta_ladder(enabled::Bool, start::Float64, step::Float64, stop::Float64) -> Vector{Float64}

Generate beta initialization ladder for multistart.
"""
function _generate_beta_ladder(enabled::Bool, start::Float64, step::Float64, stop::Float64)
    if !enabled
        return [0.0, 0.5, 1.0]
    end
    
    vals = Float64[]
    current = max(start, 0.01)
    stop = max(current, stop)
    step = max(step, 1e-9)
    
    while current <= stop + 1e-9
        push!(vals, clamp(round(current; digits=10), 0.0, 1.0))
        current += step
    end
    
    if isempty(vals) || vals[end] < stop - 1e-9
        push!(vals, stop)
    end
    
    unique!(vals)
    return vals
end

"""
    _compute_h_C_var_history_from_z(z_tube::Array{Float64,3}, H::Int, P::Int, Kd::Int) -> Array{Float64,3}

Compute h_C variable history from z_tube by taking max of positive/negative pairs.
"""
function _compute_h_C_var_history_from_z(z_tube::Array{Float64,3}, H::Int, P::Int, Kd::Int)
    h_C_var_history = zeros(Float64, H + 1, P, Kd)
    for n in 1:H, p in 1:P, q in 1:Kd
        h_C_var_history[n+1, p, q] = max(z_tube[n, p, 2*q-1], z_tube[n, p, 2*q])
    end
    return h_C_var_history
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
