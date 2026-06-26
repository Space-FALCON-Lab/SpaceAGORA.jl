module ConstellationDesign

include("constellation_utils.jl")
include("constellation_slots.jl")
include("constellation_plotting.jl")
include("constellation_config.jl")
include("stage0_seeding.jl")
include("stage0_fhsg.jl")
include("stage0p5_generator_propagation.jl")
include("stage1_design.jl")
include("stage1_controllable.jl")
include("stage2_verification.jl")
include("stage2_ocp_verification.jl")
include("math/support_functions.jl")
include("math/two_body_linearization.jl")
include("math/access_matrix.jl")
include("physics/dynamics_adapter.jl")
include("data/debris_catalog.jl")
include("visualization/plotly_plots.jl")
include("visualization/plot_config.jl")
include("test/stage0_cache_loader.jl")
include("test/parity_test_harness.jl")
include("constellation_design_rl/SatelliteSeedingEnv.jl")
include("constellation_design_rl/DeepSetPolicy.jl")
include("constellation_design_rl/PPOTrainer.jl")
include("constellation_design_rl/CostFunction.jl")
include("constellation_design_rl/Scenarios.jl")
include("constellation_design_rl/capo_integration.jl")
include("constellation_design_rl/src/ROEBoundsCalculator.jl")
include("constellation_design_rl/src/ClusterCombinations.jl")
include("constellation_design_rl/src/TrainingScenarioBuilder.jl")
include("constellation_design_rl/src/TrainingOrchestrator.jl")
include("constellation_design_rl/src/ComparisonHarness.jl")

# debris_controllable_sim functions are defined directly in this module

"""
    parse_debris_args(args::Vector{String}) -> Dict{String,Any}

Parse command-line arguments for debris_controllable_sim.
"""
function parse_debris_args(args::Vector{String})
    parsed = Dict{String,Any}()
    
    for arg in args
        if contains(arg, "=")
            key, value = split(arg, "="; limit=2)
            parsed[strip(key)] = strip(value)
        end
    end
    
    return parsed
end

"""
    run_debris_controllable_sim(config_path::AbstractString; kwargs...) -> Dict{String,Any}

Run the full debris controllable simulation pipeline with all three stages.
"""
function run_debris_controllable_sim(config_path::AbstractString; kwargs...)
    config_dict = ingest_yaml(config_path)
    for (key, value) in kwargs
        config_dict[key] = value
    end
    return run_debris_controllable_sim(config_dict)
end

"""
    run_debris_controllable_sim(config_dict::AbstractDict) -> Dict{String,Any}

Run the full debris controllable simulation pipeline with all three stages.
"""
function run_debris_controllable_sim(config_dict::AbstractDict)
    constellation_log_init!(config_dict; context="debris_controllable_sim")
    
    try
        constellation_log("pipeline", "Starting debris_controllable_sim pipeline")
        
        opt_params = get(config_dict, "optimizer_params", Dict{String,Any}())
        mode = lowercase(String(get(opt_params, "mode", "full")))
        
        results = Dict{String,Any}()
        
        if mode in ("full", "stochastic_greedy", "stage0", "rl")
            stage0_method = get(opt_params, "stage0_method", "stochastic")
            constellation_log("pipeline", "Running Stage 0: Seeding (method=$stage0_method)")
            results["stage0"] = run_stage0_seeding(config_dict; method=stage0_method)
        end
        
        # Stage 0.5: High-Fidelity Generator Propagation (optional)
        use_hf_generators = Bool(get(opt_params, "use_high_fidelity_generators", false))
        if use_hf_generators && mode in ("full", "stage0p5")
            constellation_log("pipeline", "Running Stage 0.5: High-Fidelity Generator Propagation")
            results["stage0p5"] = run_high_fidelity_generator_propagation(config_dict)
            
            # Update config_dict with time-varying generators
            config_dict["h_fwd_exact_coeffs"] = results["stage0p5"]["h_fwd_time_varying"]
            config_dict["h_Wcorr_coeffs"] = results["stage0p5"]["h_Wcorr_time_varying"]
            config_dict["generator_mode"] = "time_varying"
        else
            # Use static FHSG coefficients
            config_dict["generator_mode"] = "static"
        end
        
        if mode in ("full", "heuristic", "stage1")
            constellation_log("pipeline", "Running Stage 1: Constellation Design")
            results["stage1"] = run_constellation_design(config_dict)
        end
        
        if mode in ("full", "stage2")
            constellation_log("pipeline", "Running Stage 2: Control Verification")
            results["stage2"] = run_stage2_verification(config_dict)
        end
        
        constellation_log("pipeline", "Pipeline completed successfully")
        return results
    catch err
        constellation_log_exception("pipeline", err)
        rethrow(err)
    finally
        constellation_log_close!()
    end
end

"""
    run_capo_pipeline(config_dict::AbstractDict) -> Dict{String,Any}

Run the full CAPO pipeline (alias for run_debris_controllable_sim).
"""
function run_capo_pipeline(config_dict::AbstractDict)
    return run_debris_controllable_sim(config_dict)
end

function run_capo_pipeline(config_path::AbstractString)
    return run_debris_controllable_sim(config_path)
end

export run_stage0_seeding, run_constellation_design, run_stage2_verification, run_capo_pipeline
export run_debris_controllable_sim, parse_debris_args

# RL exports
export SatelliteSeedingEnv, ConstellationRLSatelliteSeedingObservation
export DualDeepSetPolicy, check_gpu_availability
export PPOConfig, train_ppo, setup_tensorboard_logging
export compute_stage0_cost
export load_training_scenarios, load_cluster_from_csv, load_superset_csv
export run_constellation_rl_stage0_seeding, run_policy_inference, load_trained_policy, run_stochastic_greedy_seeding
export compute_orbital_bounds_from_cluster, roe_to_orbital_elements, orbital_elements_to_roe, compute_reference_orbit, compute_shell_bounds_from_roe
export load_labeled_debris_csv, filter_cluster, generate_cluster_combinations, generate_specific_combinations
export convert_to_orbital_elements, write_client_csv, build_cluster_scenario, build_all_scenarios
export TrainingScenario, build_training_scenarios, generate_parameter_combinations
export compute_difficulty, save_training_scenarios, load_training_scenarios
export TrainingResult, train_scenario, train_campaign, generate_campaign_summary
export ComparisonResult, run_comparison, run_comparison_batch, generate_comparison_report, build_comparison_config

# Utility exports
using .ConstellationUtils
using .ConstellationSlots
using .ConstellationPlotting
using .ConstellationConfig
using .Stage0Seeding
using .Stage0FHSG
using .Stage0p5GeneratorPropagation
using .Stage1Design
using .Stage1Controllable
using .Stage2Verification
using .Stage2OCPVerification
using .SupportFunctions
using .TwoBodyLinearization
using .AccessMatrix
using .PhysicsAdapter
using .DebrisCatalog
using .PlotlyPlots
using .PlotConfig
using .Stage0CacheLoader
using .ParityTestHarness
using .SatelliteSeedingEnv
using .DeepSetPolicy
using .PPOTrainer
using .CostFunction
using .Scenarios
using .CapoIntegration
using .ROEBoundsCalculator
using .ClusterCombinations
using .TrainingScenarioBuilder
using .TrainingOrchestrator
using .ComparisonHarness

end # module ConstellationDesign
