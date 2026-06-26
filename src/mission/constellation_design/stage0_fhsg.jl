module Stage0FHSG

using LinearAlgebra
using Random
using DataFrames
using ..SupportFunctions
using ..ConstellationConfig
using ..ConstellationUtils

"""
    _fhsg_safe_build_audit(h_fwd::AbstractArray{<:Real,3}, pctx::AbstractDict) -> Dict{String,Any}

Phase A of FHSG: Greedy satellite accumulation for keepout clearance.

Computes backward-lifted controllable tube h_C and evaluates keepout safety margins.

# Arguments
- `h_fwd::AbstractArray{<:Real,3}`: Forward reachable-set support coefficients [Kd, P, H]
- `pctx::AbstractDict`: Pipeline context with parameters

# Returns
- `Dict{String,Any}`: Audit results including residuals, costs, and feasibility status
"""
function _fhsg_safe_build_audit(
    h_fwd::AbstractArray{<:Real,3},
    pctx::AbstractDict,
)
    Kd = Int(pctx["Kd"])
    P = size(h_fwd, 2)
    H = size(h_fwd, 3)
    keepout_dir_indices = Vector{Int}(pctx["keepout_dir_indices"])
    Ks = length(keepout_dir_indices)
    active_safe_horizons = Int(get(pctx, "active_safe_horizons", H))
    unsafe_threshold = Float64(pctx["unsafe_threshold"])
    safe_threshold = Float64(pctx["safe_threshold"])
    unsafe_weight = Float64(pctx["unsafe_weight"])
    safe_weight = Float64(pctx["safe_weight"])
    backward_lift_coeffs = pctx["backward_lift_coeffs"]

    h_C = zeros(Float64, H + 1, P, Kd)
    unsafe_residuals = zeros(Float64, Kd, P, H)
    safe_best_values = zeros(Float64, P, H)
    safe_best_indices = ones(Int, P, H)
    safe_residuals = zeros(Float64, P, H)
    pred_residuals = zeros(Float64, Kd, P, H)

    # Backward lift to compute controllable tube h_C
    # Note: CAPO LADS tube uses per-horizon forward coefficients h_fwd_exact_coeffs[n][q,i,p]
    # Current implementation uses aggregated h_fwd[Kd, P, H] which is functionally equivalent
    # for the purpose of Stage 1 optimization with Phi_n mapping
    for p in 1:P
        h_C[H + 1, p, :] .= Float64.(h_fwd[:, p, H])
        for n in H:-1:1
            for q in 1:Kd
                h_C[n, p, q] = sum(
                    Float64(backward_lift_coeffs[j, q, p, n]) * h_C[n + 1, p, j]
                    for j in 1:Kd)
            end
        end
    end

    min_safe_margin = Inf
    for p in 1:P, n in 1:H
        if Ks > 0
            best_qi = argmax(Float64[Float64(h_fwd[keepout_dir_indices[qi], p, n]) for qi in 1:Ks])
            best_kq = keepout_dir_indices[best_qi]
            safe_best_values[p, n] = Float64(h_fwd[best_kq, p, n])
            worst_safe_resid = 0.0
            worst_safe_idx = 1
            for qi in 1:Ks
                kq = keepout_dir_indices[qi]
                hkq = Float64(h_fwd[kq, p, n])
                if n <= active_safe_horizons
                    unsafe_residuals[kq, p, n] = max(0.0, unsafe_threshold - hkq)
                end
                safe_resid = max(0.0, safe_threshold - hkq)
                if safe_resid > worst_safe_resid
                    worst_safe_resid = safe_resid
                    worst_safe_idx = qi
                end
            end
            if n <= active_safe_horizons
                safe_residuals[p, n] = worst_safe_resid
                safe_best_indices[p, n] = keepout_dir_indices[worst_safe_idx]
                min_safe_margin = min(min_safe_margin,
                    minimum(Float64(h_fwd[keepout_dir_indices[qi], p, n]) - safe_threshold for qi in 1:Ks))
            else
                safe_best_indices[p, n] = best_kq
            end
        else
            safe_best_indices[p, n] = argmax(view(h_C, n + 1, p, :))
            safe_best_values[p, n] = maximum(view(h_C, n + 1, p, :))
            if n <= active_safe_horizons
                c0_best = maximum(view(h_C, n, p, :))
                safe_residuals[p, n] = max(0.0, safe_threshold - safe_best_values[p, n],
                    safe_threshold - c0_best)
                min_safe_margin = min(min_safe_margin, safe_best_values[p, n] - safe_threshold,
                    c0_best - safe_threshold)
            end
        end
    end
    unsafe_cost = unsafe_weight * sum(abs2, unsafe_residuals)
    safe_cost = safe_weight * sum(abs2, safe_residuals)
    min_safe_margin = isfinite(min_safe_margin) ? min_safe_margin : 0.0
    return Dict{String,Any}(
        "certificate_mode" => "minimum_required_invariant_tube",
        "unsafe_residuals" => unsafe_residuals,
        "safe_residuals" => safe_residuals,
        "reach_residuals" => zeros(Float64, 0, P, H),
        "pred_residuals" => pred_residuals,
        "safe_best_values" => safe_best_values,
        "safe_best_indices" => safe_best_indices,
        "h_C" => h_C,
        "unsafe_cost" => unsafe_cost,
        "safe_cost" => safe_cost,
        "reach_cost" => 0.0,
        "pred_cost" => 0.0,
        "total_cost" => unsafe_cost + safe_cost,
        "total_deficit" => unsafe_cost + safe_cost,
        "margin_sum" => -(unsafe_cost + safe_cost),
        "effective_margin" => -(unsafe_cost + safe_cost),
        "effective_objective" => unsafe_cost + safe_cost,
        "worst_pred_deficit" => 0.0,
        "worst_reach_deficit" => 0.0,
        "min_safe_margin" => min_safe_margin,
        "nontrivial_safe" => min_safe_margin >= -1.0e-12,
        "feasible" => maximum(unsafe_residuals) <= Float64(pctx["active_unsafe_tol"]) &&
            maximum(safe_residuals) <= Float64(pctx["active_safe_tol"]),
    )
end

"""
    _fhsg_fixed_target_pred_audit(h_fwd::AbstractArray{<:Real,3}, h_neg::AbstractArray{<:Real,3}, pctx::AbstractDict) -> Dict{String,Any}

Phase B of FHSG: Predecessor condition satisfaction with fixed target prediction.

Computes predecessor residuals and keepout floors for predecessor repair.

# Arguments
- `h_fwd::AbstractArray{<:Real,3}`: Forward reachable-set support coefficients [Kd, P, H]
- `h_neg::AbstractArray{<:Real,3}`: Negative authority support coefficients [Kd, P, H]
- `pctx::AbstractDict`: Pipeline context with parameters

# Returns
- `Dict{String,Any}`: Audit results including predecessor residuals and costs
"""
function _fhsg_fixed_target_pred_audit(h_fwd::AbstractArray{<:Real,3},
    h_neg::AbstractArray{<:Real,3},
    pctx::AbstractDict)
    Kd = Int(pctx["Kd"])
    P = size(h_fwd, 2)
    H = size(h_fwd, 3)
    support_lift_coeffs = pctx["support_lift_coeffs"]
    backward_lift_coeffs = pctx["backward_lift_coeffs"]
    keepout_dir_indices = Vector{Int}(pctx["keepout_dir_indices"])
    Ks = length(keepout_dir_indices)
    active_safe_horizons = Int(get(pctx, "active_safe_horizons", H))
    unsafe_threshold = Float64(pctx["unsafe_threshold"])
    safe_threshold = Float64(pctx["safe_threshold"])
    unsafe_weight = Float64(pctx["unsafe_weight"])
    safe_weight = Float64(pctx["safe_weight"])
    pred_weight = Float64(pctx["pred_weight"])
    certificate_support_margin_m = max(0.0, Float64(get(pctx, "certificate_support_margin_m", 0.0)))

    # Per-direction keepout floor: h_{B(d_safe)}(d_q) = safe_threshold * ||d_pos_q||.
    D_dirs = Matrix{Float64}(pctx["D_dirs"])
    d_pos_norms = [norm(view(D_dirs, 1:3, q)) for q in 1:Kd]
    keepout_floor = [d_pos_norms[q] > 1e-10 ? safe_threshold * d_pos_norms[q] : 0.0 for q in 1:Kd]

    unsafe_residuals = zeros(Float64, Kd, P, H)
    safe_best_values = zeros(Float64, P, H)
    safe_best_indices = ones(Int, P, H)
    safe_residuals = zeros(Float64, P, H)
    pred_residuals = zeros(Float64, Kd, P, H)

    # Minimum feasible Stage 1 tube z_min[q,p,n]: keepout at n=1, predecessor-lifted forward.
    z_min = zeros(Float64, Kd, P, H)
    for p in 1:P
        for q in 1:Kd
            z_min[q, p, 1] = keepout_floor[q]
        end
        for n in 2:H
            for q in 1:Kd
                sl_prev = sum(Float64(support_lift_coeffs[j, q, p, n]) * z_min[j, p, n - 1] for j in 1:Kd)
                z_min[q, p, n] = max(keepout_floor[q], max(0.0, sl_prev - Float64(h_neg[q, p, n])))
            end
        end
    end

    for p in 1:P, n in 1:H
        if Ks > 0
            best_qi = argmax(Float64[Float64(h_fwd[keepout_dir_indices[qi], p, n]) for qi in 1:Ks])
            best_kq = keepout_dir_indices[best_qi]
            safe_best_values[p, n] = Float64(h_fwd[best_kq, p, n])
            
            for qi in 1:Ks
                kq = keepout_dir_indices[qi]
                hkq = Float64(h_fwd[kq, p, n])
                if n <= active_safe_horizons
                    unsafe_residuals[kq, p, n] = max(0.0, unsafe_threshold - hkq)
                end
            end
            
            safe_resid = max(0.0, safe_threshold - safe_best_values[p, n])
            safe_residuals[p, n] = safe_resid
            safe_best_indices[p, n] = best_kq
        else
            safe_best_indices[p, n] = argmax(view(h_fwd, :, p, n))
            safe_best_values[p, n] = maximum(view(h_fwd, :, p, n))
            safe_residuals[p, n] = max(0.0, safe_threshold - safe_best_values[p, n])
        end
        
        # Predecessor residual: deficit between z_min and available forward support
        for q in 1:Kd
            pred_residuals[q, p, n] = max(0.0, z_min[q, p, n] - Float64(h_fwd[q, p, n]))
        end
    end

    unsafe_cost = unsafe_weight * sum(abs2, unsafe_residuals)
    safe_cost = safe_weight * sum(abs2, safe_residuals)
    pred_cost = pred_weight * sum(abs2, pred_residuals)
    
    return Dict{String,Any}(
        "certificate_mode" => "fixed_target_predecessor",
        "unsafe_residuals" => unsafe_residuals,
        "safe_residuals" => safe_residuals,
        "reach_residuals" => zeros(Float64, 0, P, H),
        "pred_residuals" => pred_residuals,
        "safe_best_values" => safe_best_values,
        "safe_best_indices" => safe_best_indices,
        "z_min" => z_min,
        "keepout_floor" => keepout_floor,
        "unsafe_cost" => unsafe_cost,
        "safe_cost" => safe_cost,
        "reach_cost" => 0.0,
        "pred_cost" => pred_cost,
        "total_cost" => unsafe_cost + safe_cost + pred_cost,
        "total_deficit" => unsafe_cost + safe_cost + pred_cost,
        "margin_sum" => -(unsafe_cost + safe_cost + pred_cost),
        "effective_margin" => -(unsafe_cost + safe_cost + pred_cost),
        "effective_objective" => unsafe_cost + safe_cost + pred_cost,
        "worst_pred_deficit" => maximum(pred_residuals),
        "worst_reach_deficit" => 0.0,
        "min_safe_margin" => minimum(safe_residuals),
        "nontrivial_safe" => true,
        "feasible" => maximum(unsafe_residuals) <= Float64(pctx["active_unsafe_tol"]) &&
            maximum(safe_residuals) <= Float64(pctx["active_safe_tol"]) &&
            maximum(pred_residuals) <= Float64(pctx["active_pred_tol"]),
    )
end

"""
    run_fhsg_stage0(config_dict::AbstractDict) -> Dict{String,Any}

Run the full FHSG Stage 0 seeding process.

# Arguments
- `config_dict::AbstractDict`: Configuration dictionary

# Returns
- `Dict{String,Any}`: Stage 0 results including seed satellites and support coefficients
"""
function run_fhsg_stage0(config_dict::AbstractDict)
    constellation_log("stage0_fhsg", "Starting FHSG Stage 0 seeding")
    
    # Check for cached Stage 0 data
    opt_params = config_dict["optimizer_params"]
    if Bool(get(opt_params, "use_cached_stage0", false))
        cache_path = make_stage0_cache_path(config_dict)
        context = _init_fhsg_context(config_dict)
        signature = compute_cache_signature(config_dict, context)
        cached_data = load_cached_stage0(cache_path, signature)
        if cached_data !== nothing
            constellation_log("stage0_fhsg", "Using cached Stage 0 data")
            return cached_data
        end
    end
    
    # Initialize pipeline context
    pctx = _init_fhsg_context(config_dict)
    
    # Load debris and generate candidate satellites
    debris_df = load_debris_catalog(config_dict)
    debris_ics = debris_to_initial_conditions(debris_df, config_dict)
    
    # Generate candidate satellite bank
    candidate_bank = _generate_candidate_bank(config_dict, debris_ics)
    
    # Compute support coefficients
    h_fwd_exact_coeffs = _compute_h_fwd_exact(config_dict, candidate_bank)
    h_Wcorr_coeffs = _compute_h_Wcorr(config_dict, candidate_bank)
    support_lift_coeffs = _compute_support_lift(config_dict)
    backward_lift_coeffs = _compute_backward_lift(config_dict)
    
    # Phase A: Safe build audit
    pctx["backward_lift_coeffs"] = backward_lift_coeffs
    phase_a_result = _fhsg_safe_build_audit(h_fwd_exact_coeffs, pctx)
    
    # Phase B: Predecessor repair audit
    pctx["support_lift_coeffs"] = support_lift_coeffs
    phase_b_result = _fhsg_fixed_target_pred_audit(h_fwd_exact_coeffs, h_Wcorr_coeffs, pctx)
    
    # Select seed satellites based on audit results
    seed_indices = _select_seed_satellites(phase_a_result, phase_b_result, config_dict)
    
    constellation_log("stage0_fhsg", "FHSG Stage 0 seeding complete")
    
    result = Dict{String,Any}(
        "seed_indices" => seed_indices,
        "h_fwd_exact_coeffs" => h_fwd_exact_coeffs,
        "h_Wcorr_coeffs" => h_Wcorr_coeffs,
        "support_lift_coeffs" => support_lift_coeffs,
        "backward_lift_coeffs" => backward_lift_coeffs,
        "phase_a_audit" => phase_a_result,
        "phase_b_audit" => phase_b_result,
        "candidate_bank" => candidate_bank,
    )
    
    # Cache the result
    if Bool(get(opt_params, "use_cached_stage0", false))
        cache_path = make_stage0_cache_path(config_dict)
        signature = compute_cache_signature(config_dict, pctx)
        save_cached_stage0(cache_path, signature, result)
    end
    
    return result
end

"""
    _init_fhsg_context(config_dict::AbstractDict) -> Dict{String,Any}

Initialize FHSG pipeline context from configuration.

# Arguments
- `config_dict::AbstractDict`: Configuration dictionary

# Returns
- `Dict{String,Any}`: Pipeline context
"""
function _init_fhsg_context(config_dict::AbstractDict)
    opt_params = config_dict["optimizer_params"]
    effector_params = config_dict["effector_params"]
    fhsg_params = get(opt_params, "finite_horizon_stochastic_greedy_seeding", Dict{String,Any}())
    
    Kd = Int(opt_params["controllable_set_dirs"])
    keepout_dirs = _build_keepout_position_dirs(Kd)
    keepout_dir_indices = collect(1:Kd)
    
    mission = config_dict["mission"]
    
    return Dict{String,Any}(
        "Kd" => Kd,
        "keepout_dir_indices" => keepout_dir_indices,
        "keepout_dirs" => keepout_dirs,
        "D_dirs" => _build_polyhedral_dirs(Kd),
        "unsafe_threshold" => effector_params["range"],
        "safe_threshold" => effector_params["range"] * 0.9,
        "unsafe_weight" => Float64(get(fhsg_params, "unsafe_weight", 1.0)),
        "safe_weight" => Float64(get(fhsg_params, "safe_weight", 1.0)),
        "pred_weight" => Float64(get(fhsg_params, "pred_weight", 1.0)),
        "active_unsafe_tol" => Float64(get(fhsg_params, "active_unsafe_tol", 1.0e-8)),
        "active_safe_tol" => Float64(get(fhsg_params, "active_safe_tol", 1.0e-8)),
        "active_pred_tol" => Float64(get(fhsg_params, "active_pred_tol", 1.0e-8)),
        "active_safe_horizons" => Int(mission["n_horizons"]),
        "certificate_support_margin_m" => Float64(get(fhsg_params, "certificate_support_margin_m", 0.0)),
    )
end

"""
    _generate_candidate_bank(config_dict::AbstractDict, debris_ics::Dict{Symbol, Vector{Float64}}) -> DataFrame

Generate candidate satellite bank from debris initial conditions.

# Arguments
- `config_dict::AbstractDict`: Configuration dictionary
- `debris_ics::Dict{Symbol, Vector{Float64}}`: Debris initial conditions

# Returns
- `DataFrame`: Candidate satellite bank
"""
function _generate_candidate_bank(config_dict::AbstractDict, debris_ics::Dict{Symbol, Vector{Float64}})
    constellation_params = config_dict["constellation_params"]
    
    # Sample orbital elements within bounds
    n_candidates = constellation_params["max_sats"]
    
    rng = MersenneTwister(config_dict["sim_params"]["rng_seed"])
    
    a = rand(rng, Uniform(constellation_params["a_min"], constellation_params["a_max"]), n_candidates)
    e = rand(rng, Uniform(0.0, constellation_params["e_max"]), n_candidates)
    inc = rand(rng, Uniform(constellation_params["i_min"], constellation_params["i_max"]), n_candidates)
    raan = rand(rng, Uniform(constellation_params["raan_min"], constellation_params["raan_max"]), n_candidates)
    arg_p = rand(rng, Uniform(0.0, 2π), n_candidates)
    ta = rand(rng, Uniform(0.0, 2π), n_candidates)
    
    return DataFrame(
        a = a,
        e = e,
        inc = inc,
        raan = raan,
        arg_p = arg_p,
        ta = ta,
    )
end

"""
    _compute_h_fwd_exact(config_dict::AbstractDict, candidate_bank::DataFrame) -> Array{Float64,3}

Compute forward reachable-set support coefficients.

# Arguments
- `config_dict::AbstractDict`: Configuration dictionary
- `candidate_bank::DataFrame`: Candidate satellite bank

# Returns
- `Array{Float64,3}`: Forward support coefficients [Kd, M, P]
"""
function _compute_h_fwd_exact(config_dict::AbstractDict, candidate_bank::DataFrame)
    # Placeholder: actual implementation requires full reachability computation
    # This is a simplified version for structure
    Kd = Int(config_dict["optimizer_params"]["controllable_set_dirs"])
    M = nrow(candidate_bank)
    P = config_dict["client_bounds"]["n_clients"]
    
    return zeros(Float64, Kd, M, P)
end

"""
    _compute_h_Wcorr(config_dict::AbstractDict, candidate_bank::DataFrame) -> Array{Float64,3}

Compute correction authority support coefficients.

# Arguments
- `config_dict::AbstractDict`: Configuration dictionary
- `candidate_bank::DataFrame`: Candidate satellite bank

# Returns
- `Array{Float64,3}`: Correction authority coefficients [Kd, M, P]
"""
function _compute_h_Wcorr(config_dict::AbstractDict, candidate_bank::DataFrame)
    # Placeholder: actual implementation requires full authority computation
    Kd = Int(config_dict["optimizer_params"]["controllable_set_dirs"])
    M = nrow(candidate_bank)
    P = config_dict["client_bounds"]["n_clients"]
    
    return zeros(Float64, Kd, M, P)
end

"""
    _compute_support_lift(config_dict::AbstractDict) -> Array{Float64,4}

Compute SL (successor-lift) matrix for predecessor constraints.

# Arguments
- `config_dict::AbstractDict`: Configuration dictionary

# Returns
- `Array{Float64,4}`: SL coefficients [Kd, Kd, P, H]
"""
function _compute_support_lift(config_dict::AbstractDict)
    Kd = Int(config_dict["optimizer_params"]["controllable_set_dirs"])
    P = config_dict["client_bounds"]["n_clients"]
    H = Int(config_dict["optimizer_params"]["n_horizons"])
    
    # Placeholder: actual implementation computes SL matrix from state transition
    return zeros(Float64, Kd, Kd, P, H)
end

"""
    _compute_backward_lift(config_dict::AbstractDict) -> Array{Float64,4}

Compute backward lift coefficients for controllable tube computation.

# Arguments
- `config_dict::AbstractDict`: Configuration dictionary

# Returns
- `Array{Float64,4}`: Backward lift coefficients [Kd, Kd, P, H]
"""
function _compute_backward_lift(config_dict::AbstractDict)
    Kd = Int(config_dict["optimizer_params"]["controllable_set_dirs"])
    P = config_dict["client_bounds"]["n_clients"]
    H = Int(config_dict["optimizer_params"]["n_horizons"])
    
    # Placeholder: actual implementation computes backward lift from state transition
    return zeros(Float64, Kd, Kd, P, H)
end

"""
    _select_seed_satellites(phase_a_result::AbstractDict, phase_b_result::AbstractDict, config_dict::AbstractDict) -> Vector{Int}

Select seed satellites based on audit results.

# Arguments
- `phase_a_result::AbstractDict`: Phase A audit results
- `phase_b_result::AbstractDict`: Phase B audit results
- `config_dict::AbstractDict`: Configuration dictionary

# Returns
- `Vector{Int}`: Selected seed satellite indices
"""
function _select_seed_satellites(phase_a_result::AbstractDict, phase_b_result::AbstractDict, config_dict::AbstractDict)
    # Placeholder: actual implementation selects satellites based on audit results
    max_sats = config_dict["constellation_params"]["max_sats"]
    return collect(1:min(10, max_sats))
end

export _fhsg_safe_build_audit, _fhsg_fixed_target_pred_audit, run_fhsg_stage0

end # module Stage0FHSG
