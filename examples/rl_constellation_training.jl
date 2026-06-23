#!/usr/bin/env julia

# RL Constellation Training Example
# This script demonstrates how to train a DeepSet-based RL agent for satellite seeding
# Usage: julia examples/rl_constellation_training.jl [--resume CHECKPOINT_PATH]

using SpaceAGORA.ConstellationDesign
using ProgressMeter
using Dates

"""
    parse_args(args::Vector{String}) -> Dict{String,Any}

Parse command-line arguments.
"""
function parse_args(args::Vector{String})
    parsed = Dict{String,Any}("resume_from" => "")
    
    i = 1
    while i <= length(args)
        arg = args[i]
        if startswith(arg, "--resume=")
            parsed["resume_from"] = replace(arg, "--resume=" => "")
        elseif arg == "--resume" && i < length(args)
            parsed["resume_from"] = args[i + 1]
            i += 1  # Skip next argument since we consumed it
        end
        i += 1
    end
    
    return parsed
end

"""
    main()

Main training function.
"""
function main(args::Vector{String}=ARGS)
    println("=== RL Constellation Training ===")
    println("Starting at: $(Dates.now())")
    
    # Parse command-line arguments
    cli_args = parse_args(args)
    resume_from = cli_args["resume_from"]
    
    if !isempty(resume_from)
        println("Resume mode: will resume from checkpoint: $resume_from")
    end
    
    # Load training scenarios
    println("\n[1/5] Loading training scenarios...")
    scenarios = load_training_scenarios()
    println("Loaded $(length(scenarios)) training scenarios")
    
    if isempty(scenarios)
        error("No training scenarios found. Please ensure debris cluster CSV files are in data/debris_clusters/")
    end
    
    # Configure PPO hyperparameters
    println("\n[2/5] Configuring PPO hyperparameters...")
    ppo_config = PPOConfig(
        # PPO parameters
        clip_range = 0.2,
        epochs = 10,
        batch_size = 64,
        gamma = 0.99,
        gae_lambda = 0.95,
        
        # Network parameters
        hidden_size = 32,
        embed_size = 16,
        learning_rate = 3e-4,
        
        # Training parameters
        max_episodes = 10000,
        max_steps_per_episode = 100,
        update_frequency = 2048,
        
        # Parallel training
        num_workers = 32,
        num_envs_per_worker = 4,
        
        # Logging
        tensorboard_log_dir = "data/rl_models/tensorboard",
        log_frequency = 10,
        
        # Reward weights
        unsafe_weight = 100.0,
        safe_weight = 10.0,
        pred_weight = 5.0,
        feasibility_threshold = 1e-6,
        
        # Checkpoint resume
        resume_from = resume_from,
    )
    
    println("PPO Configuration:")
    println("  - Max episodes: $(ppo_config.max_episodes)")
    println("  - Max steps per episode: $(ppo_config.max_steps_per_episode)")
    println("  - Hidden size: $(ppo_config.hidden_size)")
    println("  - Embed size: $(ppo_config.embed_size)")
    println("  - Learning rate: $(ppo_config.learning_rate)")
    println("  - Num workers: $(ppo_config.num_workers)")
    
    # Setup TensorBoard logging
    println("\n[3/5] Setting up TensorBoard logging...")
    lg = setup_tensorboard_logging(ppo_config.tensorboard_log_dir)
    println("TensorBoard logs will be written to: $(ppo_config.tensorboard_log_dir)")
    println("View logs with: tensorboard --logdir=$(ppo_config.tensorboard_log_dir)")
    
    # Check GPU availability
    println("\n[4/5] Checking GPU availability...")
    gpu_available, device = check_gpu_availability()
    println("Device: $(gpu_available ? "GPU (CUDA)" : "CPU")")
    
    # Train PPO agent
    println("\n[5/5] Starting PPO training...")
    println("This may take a while depending on hardware...")
    
    @time policy = train_ppo(ppo_config, scenarios)
    
    println("\n=== Training Completed ===")
    println("Finished at: $(Dates.now())")
    println("Final model saved to: data/rl_models/latest_model.jld2")
    println("\nTo view training progress:")
    println("  tensorboard --logdir=$(ppo_config.tensorboard_log_dir)")
    
    return policy
end

# Run main function if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
