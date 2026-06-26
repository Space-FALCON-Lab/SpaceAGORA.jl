module ConstellationConfig

using ..ConstellationUtils
using YAML

"""
    parse_constellation_config(config_path::AbstractString) -> Dict{String,Any}

Parse a constellation configuration YAML file and return a validated configuration dictionary.

# Arguments
- `config_path::AbstractString`: Path to YAML configuration file

# Returns
- Dictionary containing validated constellation configuration
"""
function parse_constellation_config(config_path::AbstractString)
    config_dict = ingest_yaml(config_path)
    return validate_constellation_config(config_dict)
end

"""
    validate_constellation_config(config_dict::AbstractDict) -> Dict{String,Any}

Validate and normalize a constellation configuration dictionary.

# Arguments
- `config_dict::AbstractDict`: Configuration dictionary

# Returns
- Validated and normalized configuration dictionary
"""
function validate_constellation_config(config_dict::AbstractDict)
    # Ensure required top-level sections exist
    config = copy(config_dict)
    
    # Mission parameters
    if !haskey(config, "mission")
        config["mission"] = Dict{String,Any}()
    end
    mission = config["mission"]
    get!(mission, "type", "lads_debris_steering")
    get!(mission, "drm", "lads_publication")
    get!(mission, "n_horizons", 5)
    get!(mission, "planning_horizon_hours", 48.0)
    get!(mission, "mission_duration_hours", 240.0)
    
    # Physical constants
    if !haskey(config, "physical_constants")
        config["physical_constants"] = Dict{String,Any}()
    end
    phys = config["physical_constants"]
    get!(phys, "mu", 3.986004418e14)
    get!(phys, "J2", 1.08263e-3)
    get!(phys, "radius", 6371.0e3)
    
    # Debris parameters
    if !haskey(config, "debris_params")
        config["debris_params"] = Dict{String,Any}()
    end
    debris = config["debris_params"]
    get!(debris, "data_source", "sample")
    get!(debris, "file_path", nothing)
    get!(debris, "n_clients", 10)
    get!(debris, "synthetic_seed", 42)
    get!(debris, "altitude_range_km", [400, 800])
    get!(debris, "inclination_range_deg", [0, 90])
    
    # Constellation parameters
    if !haskey(config, "constellation_params")
        config["constellation_params"] = Dict{String,Any}()
    end
    constel = config["constellation_params"]
    get!(constel, "max_sats", 100)
    get!(constel, "a_min", 7.0e6)
    get!(constel, "a_max", 8.0e6)
    get!(constel, "e_max", 0.1)
    get!(constel, "i_min", 0.0)
    get!(constel, "i_max", deg2rad(90))
    get!(constel, "raan_min", 0.0)
    get!(constel, "raan_max", 2π)
    
    # Effector parameters
    if !haskey(config, "effector_params")
        config["effector_params"] = Dict{String,Any}()
    end
    effector = config["effector_params"]
    get!(effector, "range", 600e3)
    get!(effector, "max_thrust", 0.1)
    get!(effector, "sc_mass", 100.0)
    
    # Optimizer parameters
    if !haskey(config, "optimizer_params")
        config["optimizer_params"] = Dict{String,Any}()
    end
    opt = config["optimizer_params"]
    get!(opt, "mode", "full")
    get!(opt, "stage0_mode", "fhsg")
    get!(opt, "stage1_certificate_mode", "lads_tube")
    get!(opt, "stage2_mode", "optimal_control")
    
    # Stage 0: Finite Horizon Stochastic Greedy parameters
    if !haskey(opt, "finite_horizon_stochastic_greedy_seeding")
        opt["finite_horizon_stochastic_greedy_seeding"] = Dict{String,Any}()
    end
    fhsg = opt["finite_horizon_stochastic_greedy_seeding"]
    get!(fhsg, "run_count", 5)
    get!(fhsg, "max_greedy_steps", 64)
    get!(fhsg, "stage_a_max_steps", Int(get(fhsg, "max_greedy_steps", 64)))
    get!(fhsg, "stage_b_max_steps", Int(get(fhsg, "max_greedy_steps", 64)))
    get!(fhsg, "batch_size", 1000)
    get!(fhsg, "group_add_count", 5)
    get!(fhsg, "top_k", 8)
    get!(fhsg, "stage0_cost_mode", "minimum_required_invariant_tube")
    get!(fhsg, "two_stage_mode", true)
    get!(fhsg, "unsafe_weight", 1.0)
    get!(fhsg, "safe_weight", 1.0)
    get!(fhsg, "reach_weight", 1.0)
    get!(fhsg, "pred_weight", 1.0)
    get!(fhsg, "safe_margin_epsilon", 5000.0)
    get!(fhsg, "pred_margin_epsilon", 5000.0)
    get!(fhsg, "certificate_support_margin_m", 0.0)
    get!(fhsg, "min_standoff_distance_m", 5000.0)
    get!(fhsg, "active_unsafe_tol", 1.0e-8)
    get!(fhsg, "active_safe_tol", 1.0e-8)
    get!(fhsg, "active_reach_tol", 1.0e-8)
    get!(fhsg, "active_pred_tol", 1.0e-8)
    get!(fhsg, "directional_prefilter_eta", 0.0)
    get!(fhsg, "safe_prefilter_eta", 0.0)
    get!(fhsg, "sample_relaxation_gamma", 1.0)
    get!(fhsg, "safe_direction_count", Int(get(opt, "controllable_keepout_dirs", 72)))
    get!(fhsg, "correction_direction_count", Int(get(opt, "controllable_set_dirs", 72)))
    get!(fhsg, "reuse_existing_library", true)
    get!(fhsg, "write_selected_ic_arrow", true)
    get!(fhsg, "selected_ic_arrow_path", "")
    get!(fhsg, "nc_support_sharing", Int(get(fhsg, "max_greedy_steps", 64)))
    get!(fhsg, "strict_determinism", true)
    get!(fhsg, "parallel_candidate_sampling", false)
    get!(fhsg, "parallel_candidate_scoring", true)
    get!(fhsg, "parallel_runs", true)
    get!(fhsg, "stage1_probe_enabled", false)
    get!(fhsg, "stage1_probe_skip_stage2", true)
    get!(fhsg, "stage1_probe_require_solved", true)
    get!(fhsg, "merge_shells_across_runs", true)
    get!(fhsg, "thread_batch_chunk", 32)
    get!(fhsg, "show_progress_bars", true)
    get!(fhsg, "log_progress_csv", true)
    get!(fhsg, "progress_log_every_step", true)
    get!(fhsg, "fail_if_infeasible", true)
    get!(fhsg, "sampling_max_tries", 4096)
    get!(fhsg, "use_pred_repair", false)
    get!(fhsg, "repair_interval", 10)
    
    # Stage 1 parameters
    get!(opt, "controllable_set_dirs", 72)
    get!(opt, "controllable_keepout_dirs", 72)
    get!(opt, "stage1_keepout_margin", 0.0)
    get!(opt, "stage1_min_safe_support_margin", 0.0)
    get!(opt, "stage1_beta_ladder_enabled", true)
    get!(opt, "stage1_beta_ladder_start", 0.01)
    get!(opt, "stage1_beta_ladder_step", 0.01)
    get!(opt, "stage1_beta_ladder_stop", 1.0)
    get!(opt, "stage1_beta_ladder_max_attempts", 101)
    get!(opt, "constructive_multistart_count", 101)
    get!(opt, "constructive_predecessor_contraction_min", 0.5)
    get!(opt, "certificate_mode", "support_function")
    
    # Stage 2 parameters
    get!(opt, "skip_stage2", false)
    get!(opt, "witness_alpha_reward", 1.0)
    get!(opt, "witness_predecessor_reward", 0.25)
    get!(opt, "witness_productivity_prune_enabled", false)
    get!(opt, "witness_productivity_topk", 1)
    get!(opt, "witness_productivity_metric_basis", "constructive_dictionary")
    
    # Optimization weights
    get!(opt, "num_sats_weight", 100.0)
    get!(opt, "l1_reg_weight", 10.0)
    get!(opt, "tube_size_weight", 1.0e-5)
    get!(opt, "beta_threshold", 0.1)
    get!(opt, "beta_budget", 1.0)
    
    # Solver settings
    get!(opt, "solver", "clarabel")
    get!(opt, "tolerance", 1.0e-4)
    get!(opt, "max_iter", 200)
    get!(opt, "capo_log_solver_trace", false)
    
    # Relative dynamics (CW-ECI only for LADS benchmark)
    get!(opt, "relative_dynamics_mode", "cw_eci")
    
    # Stage 1 LP warm-start (CAPO LADS tube)
    get!(opt, "stage1_lp_warm_start_enabled", true)
    get!(opt, "stage1_lp_warm_start_coeff_rel_drop", 1e-12)
    get!(opt, "stage1_lp_warm_start_verbose", false)
    get!(opt, "stage1_lp_warm_start_plot_only", false)
    
    # Stage 1 nonconvex objective (β(1-β) fractionality penalty)
    get!(opt, "stage1_nonconvex_fractionality_penalty_enabled", false)  # Per LADS invariants, default to convex
    
    # Certificate mode
    get!(opt, "certificate_mode", "support_function")
    
    # Controllable set mode
    get!(opt, "controllable_set_mode", "polyhedral")
    
    # Repeated horizon (default enabled for LADS)
    get!(opt, "repeated_horizon_enabled", true)
    get!(opt, "repeated_horizon_client_cache_enabled", true)
    get!(opt, "compact_phase_access_single_horizon", false)
    
    # LP warm start
    get!(opt, "stage1_lp_warm_start_enabled", true)
    
    # Fidelity
    get!(opt, "fidelity_mode", "keplerian")
    
    # High-fidelity generator propagation (Stage 0.5)
    get!(opt, "use_high_fidelity_generators", false)
    
    # Fidelity level preset (0-6), overrides individual toggles if set
    get!(opt, "high_fidelity_preset", 0)  # 0=pure_keplerian, 1=keplerian_j2, 2=keplerian_j2_simple_drag, 3=keplerian_j2_nrlmsise00, 4=keplerian_j2_nrlmsise00_srp, 5=keplerian_j2_nrlmsise00_srp_third_body, 6=all_disturbances
    
    # Individual disturbance toggles (if preset not set)
    get!(opt, "disturbance_enable_j2", false)
    get!(opt, "disturbance_enable_harmonics", false)
    get!(opt, "disturbance_harmonics_degree", 4)
    get!(opt, "disturbance_harmonics_order", 4)
    get!(opt, "disturbance_enable_atmosphere", false)
    get!(opt, "disturbance_atmosphere_model", "nrlmsise00")  # Options: "none", "exponential", "nrlmsise00", "gram"
    get!(opt, "disturbance_enable_srp", false)
    get!(opt, "disturbance_enable_albedo", false)
    get!(opt, "disturbance_enable_ir", false)
    get!(opt, "disturbance_enable_third_body", false)
    get!(opt, "disturbance_third_body_bodies", ["sun", "moon"])  # Options: "sun", "moon", "mercury", "venus", "mars", "jupiter", "saturn"
    
    # Propagation settings
    get!(opt, "high_fidelity_num_threads", Threads.nthreads())
    get!(opt, "high_fidelity_cache_enabled", true)
    get!(opt, "high_fidelity_cache_dir", "cache/high_fidelity_generators")
    
    # Stage 0 caching
    get!(opt, "use_cached_stage0", false)
    get!(opt, "cached_stage0_data", nothing)
    
    # Simulation parameters
    if !haskey(config, "sim_params")
        config["sim_params"] = Dict{String,Any}()
    end
    sim = config["sim_params"]
    get!(sim, "dt", 10.0)
    get!(sim, "N", 1000)
    get!(sim, "rng_seed", 67)
    
    # Plotting parameters
    if !haskey(config, "plotting")
        config["plotting"] = Dict{String,Any}()
    end
    plot = config["plotting"]
    get!(plot, "save", true)
    get!(plot, "display", false)
    get!(plot, "base_outdir", "plots/constellation")
    get!(plot, "plot_stage1_lp_warm_start_beta_histogram", true)
    get!(plot, "plot_tube_3d", true)
    get!(plot, "plot_access_matrix", true)
    get!(plot, "plot_stage2_trajectories", true)
    
    # Client bounds (derived from debris_params)
    if !haskey(config, "client_bounds")
        config["client_bounds"] = Dict{String,Any}()
    end
    config["client_bounds"]["n_clients"] = debris["n_clients"]
    
    constellation_log("config", "Configuration validated successfully")
    
    return config
end

"""
    deg2rad(deg::Real) -> Float64

Convert degrees to radians.
"""
deg2rad(deg::Real) = deg * π / 180.0

export parse_constellation_config, validate_constellation_config, deg2rad

end # module ConstellationConfig
