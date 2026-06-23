#!/usr/bin/env julia

# RL Constellation Inference Example
# This script demonstrates how to use a trained RL policy for satellite seeding

using SpaceAGORA.ConstellationDesign
using CSV
using DataFrames

"""
    main()

Main inference function.
"""
function main()
    println("=== RL Constellation Inference ===")
    println("Starting at: $(Dates.now())")
    
    # Load a specific cluster for inference
    println("\n[1/4] Loading debris cluster...")
    cluster_id = "1"
    cluster_path = "data/debris_clusters/clients_cluster_$cluster_id.csv"
    
    if !isfile(cluster_path)
        error("Cluster file not found: $cluster_path")
    end
    
    config = load_cluster_from_csv(cluster_path)
    println("Loaded cluster $cluster_id with $(config["n_clients"]) clients")
    
    # Configure RL parameters
    println("\n[2/4] Configuring RL parameters...")
    config["optimizer_params"] = Dict{String,Any}(
        "stage0_method" => "rl",
        "rl_config" => Dict{String,Any}(
            "rl_model_path" => "data/rl_models/latest_model.jld2",
            "max_sats" => 64,
            "max_steps" => 100,
            "unsafe_weight" => 100.0,
            "safe_weight" => 10.0,
            "pred_weight" => 5.0,
            "feasibility_threshold" => 1e-6,
        )
    )
    
    # Check if model exists
    model_path = config["optimizer_params"]["rl_config"]["rl_model_path"]
    if !isfile(model_path)
        error("Trained model not found: $model_path. Please train a model first using rl_constellation_training.jl")
    end
    
    println("Model path: $model_path")
    
    # Run RL seeding
    println("\n[3/4] Running RL-based seeding...")
    @time result = run_rl_stage0_seeding(config)
    
    println("RL seeding completed:")
    println("  - Number of satellites: $(result["n_sats"])")
    println("  - Model used: $(result["model_path"])")
    println("  - Training episode: $(result["training_episode"])")
    
    # Compare with stochastic greedy baseline
    println("\n[4/4] Comparing with stochastic greedy baseline...")
    @time baseline_result = run_stochastic_greedy_seeding(config)
    
    println("Stochastic greedy seeding completed:")
    println("  - Number of satellites: $(baseline_result["n_sats"])")
    
    # Compute costs for comparison
    rl_cost = compute_stage0_cost(config, result["constellation_orbitals"])
    baseline_cost = compute_stage0_cost(config, baseline_result["constellation_orbitals"])
    
    println("\n=== Cost Comparison ===")
    println("RL Seeding:")
    println("  - Total cost: $(rl_cost["total_cost"])")
    println("  - Unsafe cost: $(rl_cost["unsafe_cost"])")
    println("  - Safe cost: $(rl_cost["safe_cost"])")
    println("  - Pred cost: $(rl_cost["pred_cost"])")
    println("  - Total deficit: $(rl_cost["total_deficit"])")
    println("  - Feasible: $(rl_cost["feasible"])")
    
    println("\nStochastic Greedy:")
    println("  - Total cost: $(baseline_cost["total_cost"])")
    println("  - Unsafe cost: $(baseline_cost["unsafe_cost"])")
    println("  - Safe cost: $(baseline_cost["safe_cost"])")
    println("  - Pred cost: $(baseline_cost["pred_cost"])")
    println("  - Total deficit: $(baseline_cost["total_deficit"])")
    println("  - Feasible: $(baseline_cost["feasible"])")
    
    println("\n=== Inference Completed ===")
    println("Finished at: $(Dates.now())")
    
    return result, baseline_result
end

# Run main function if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
