module PPOTrainer

using ReinforcementLearningCore
using ReinforcementLearningZoo
using TensorBoardLogger
using Logging
using Flux
using ProgressMeter
using Dates
using ..SatelliteSeedingEnv
using ..DeepSetPolicy
using ..CostFunction
using ..Scenarios

"""
    PPOConfig

Configuration for PPO training.
"""
Base.@kwdef struct PPOConfig
    # PPO parameters
    clip_range::Float64 = 0.2
    epochs::Int = 10
    batch_size::Int = 64
    gamma::Float64 = 0.99
    gae_lambda::Float64 = 0.95
    
    # Network parameters
    hidden_size::Int = 32
    embed_size::Int = 16
    learning_rate::Float64 = 3e-4
    
    # Training parameters
    max_episodes::Int = 10000
    max_steps_per_episode::Int = 100
    update_frequency::Int = 2048
    
    # Parallel training
    num_workers::Int = 32
    num_envs_per_worker::Int = 4
    
    # Logging
    tensorboard_log_dir::String = "data/rl_models/tensorboard"
    log_frequency::Int = 10  # Log every N episodes
    
    # Reward weights
    unsafe_weight::Float64 = 100.0
    safe_weight::Float64 = 10.0
    pred_weight::Float64 = 5.0
    feasibility_threshold::Float64 = 1e-6
    
    # Checkpoint resume
    resume_from::String = ""  # Path to checkpoint file to resume from
end

"""
    setup_tensorboard_logging(log_dir::String) -> TBLogger

Setup TensorBoard logging for training.
"""
function setup_tensorboard_logging(log_dir::String)
    # Create timestamped directory
    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
    log_path = joinpath(log_dir, timestamp)
    mkpath(log_path)
    
    lg = TBLogger(log_path, min_level=Logging.Info)
    @info "TensorBoard logging initialized at: $log_path"
    @info "View logs with: tensorboard --logdir=$log_dir"
    
    return lg
end

"""
    train_ppo(config::PPOConfig, scenarios::Vector{Dict{String,Any}}) -> DualDeepSetPolicy

Train PPO agent on given scenarios. Can resume from checkpoint if config.resume_from is specified.
"""
function train_ppo(config::PPOConfig, scenarios::Vector{Dict{String,Any}})
    # Check GPU availability
    gpu_available, device = check_gpu_availability()
    
    # Setup TensorBoard logging
    lg = setup_tensorboard_logging(config.tensorboard_log_dir)
    
    # Check if resuming from checkpoint
    start_episode = 1
    if !isempty(config.resume_from)
        if isfile(config.resume_from)
            @info "Resuming from checkpoint: $(config.resume_from)"
            policy, loaded_config, loaded_episode = load_checkpoint(config.resume_from)
            start_episode = loaded_episode + 1
            @info "Resuming training from episode $start_episode"
            
            # Log resume information to TensorBoard
            with_logger(lg) do
                @info "training/resume" checkpoint_path=config.resume_from resumed_from_episode=loaded_episode
            end
        else
            @warn "Checkpoint file not found: $(config.resume_from). Starting fresh training."
            # Create fresh policy network
            policy = DualDeepSetPolicy(
                hidden_size=config.hidden_size,
                embed_size=config.embed_size,
                device=device
            )
        end
    else
        # Create fresh policy network
        policy = DualDeepSetPolicy(
            hidden_size=config.hidden_size,
            embed_size=config.embed_size,
            device=device
        )
    end
    
    # Setup optimizer
    opt = Adam(config.learning_rate)
    
    # Training loop
    @showprogress 1 "Training PPO..." for episode in start_episode:config.max_episodes
        # Select random scenario
        scenario = rand(scenarios)
        
        # Create environment
        env = SatelliteSeedingEnv(scenario)
        
        # Reset environment
        obs = RLBase.reset!(env)
        
        episode_reward = 0.0
        episode_length = 0
        
        # Episode loop
        for step in 1:config.max_steps_per_episode
            # Get action from policy
            action = policy(obs.satellite_orbitals, obs.client_trajectories)
            
            # Step environment
            next_obs, reward, done, info = env(action)
            
            episode_reward += reward
            episode_length += 1
            
            obs = next_obs
            
            if done
                break
            end
        end
        
        # Log metrics
        if episode % config.log_frequency == 0
            # Get embeddings for logging
            embeddings = get_embeddings(
                policy, 
                env.current_sats, 
                env.client_orbitals
            )
            
            with_logger(lg) do
                @info "training/episode" episode=episode total_reward=episode_reward episode_length=episode_length
                @info "training/hyperparams" lr=config.learning_rate clip_range=config.clip_range
                @info "constellation/state" n_sats=size(env.current_sats, 1) feasible=info["feasible"]
                @info "cost/total" total_cost=info["total_cost"]
                
                # Log embeddings for visualization
                @info "embeddings/satellite" sat_embeddings=embeddings.sat_φ_output
                @info "embeddings/client" client_embeddings=embeddings.client_φ_output
                @info "embeddings/cross_attention" cross_attn_output=embeddings.cross_attn_output
            end
        end
        
        # Save checkpoint periodically
        if episode % 100 == 0
            save_checkpoint(policy, config, episode)
        end
    end
    
    # Save final model
    save_checkpoint(policy, config, config.max_episodes)
    
    @info "Training completed"
    @info "Final model saved"
    
    return policy
end

"""
    save_checkpoint(policy::DualDeepSetPolicy, config::PPOConfig, episode::Int)

Save policy checkpoint to disk.
"""
function save_checkpoint(policy::DualDeepSetPolicy, config::PPOConfig, episode::Int)
    model_dir = "data/rl_models"
    mkpath(model_dir)
    
    filename = "policy_checkpoint_episode_$episode.jld2"
    filepath = joinpath(model_dir, filename)
    
    # Save using JLD2
    using JLD2
    jldsave(filepath; policy=policy, config=config, episode=episode)
    
    # Also save as latest
    latest_path = joinpath(model_dir, "latest_model.jld2")
    cp(filepath, latest_path; force=true)
    
    @info "Checkpoint saved: $filepath"
end

"""
    load_checkpoint(filepath::String) -> Tuple{DualDeepSetPolicy, PPOConfig, Int}

Load policy checkpoint from disk.
"""
function load_checkpoint(filepath::String)
    using JLD2
    data = JLD2.load(filepath)
    
    policy = data["policy"]
    config = data["config"]
    episode = data["episode"]
    
    @info "Checkpoint loaded: $filepath (episode $episode)"
    
    return policy, config, episode
end

export PPOConfig, train_ppo, setup_tensorboard_logging, save_checkpoint, load_checkpoint

end # module PPOTrainer
