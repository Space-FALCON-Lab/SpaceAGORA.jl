module TrainingOrchestrator

using ..PPOTrainer
using ..TrainingScenarioBuilder
using ..SatelliteSeedingEnv
using ..DeepSetPolicy
using ..Scenarios
using JLD2
using TensorBoardLogger
using Logging
using Dates

"""
    TrainingResult

Structure representing the result of training on a single scenario.
"""
struct TrainingResult
    scenario_id::String
    success::Bool
    final_episode::Int
    best_reward::Float64
    final_cost::Float64
    model_path::String
    tensorboard_dir::String
    training_time::Float64
    error_message::String
end

"""
    train_scenario(scenario::TrainingScenario, base_output_dir::String) -> TrainingResult

Train a single scenario and return the result.
"""
function train_scenario(scenario::TrainingScenario, base_output_dir::String)
    start_time = now()
    
    # Create output directory for this scenario
    scenario_dir = joinpath(base_output_dir, scenario.scenario_id)
    mkpath(scenario_dir)
    
    # Setup TensorBoard logging
    tensorboard_dir = joinpath(scenario_dir, "tensorboard")
    mkpath(tensorboard_dir)
    
    try
        # Build config dict from scenario
        config_dict = build_config_from_scenario(scenario, scenario_dir)
        
        # Setup PPO config
        ppo_config = build_ppo_config(scenario.training_config)
        
        # Setup TensorBoard logger
        logger = TBLogger(tensorboard_dir, tb_overwrite=true)
        
        # Load training scenarios
        training_scenarios = [scenario]
        
        # Train PPO
        policy, final_episode, best_reward = train_ppo(
            training_scenarios,
            ppo_config,
            logger = logger,
            checkpoint_dir = scenario_dir,
        )
        
        # Save final model
        model_path = joinpath(scenario_dir, "final_model.jld2")
        save_checkpoint(model_path, policy, scenario.training_config, final_episode)
        
        # Compute final cost
        env = SatelliteSeedingEnv(config_dict)
        obs = RLBase.reset!(env)
        constellation = zeros(0, 6)
        
        while true
            action = policy(obs.satellite_orbitals, obs.client_trajectories)
            next_obs, reward, done, info = env(action)
            constellation = vcat(constellation, env.current_sats[end:end, :])
            obs = next_obs
            if done
                break
            end
        end
        
        cost_result = compute_stage0_cost(config_dict, constellation)
        
        training_time = (now() - start_time) / Millisecond(1000)
        
        return TrainingResult(
            scenario.scenario_id,
            true,
            final_episode,
            best_reward,
            cost_result["total_cost"],
            model_path,
            tensorboard_dir,
            training_time,
            "",
        )
    catch err
        training_time = (now() - start_time) / Millisecond(1000)
        error_msg = sprint(showerror, err)
        
        return TrainingResult(
            scenario.scenario_id,
            false,
            0,
            -Inf,
            Inf,
            "",
            tensorboard_dir,
            training_time,
            error_msg,
        )
    end
end

"""
    build_config_from_scenario(scenario::TrainingScenario, output_dir::String) -> Dict{String,Any}

Build configuration dictionary from training scenario.
"""
function build_config_from_scenario(scenario::TrainingScenario, output_dir::String)
    config_dict = Dict{String,Any}(
        "client_orbitals" => scenario.client_orbitals,
        "orbital_bounds" => scenario.orbital_bounds,
        "optimizer_params" => Dict{String,Any}(
            "rl_config" => scenario.training_config,
        ),
    )
    
    # Add effector parameters
    for (key, value) in scenario.effector_params
        if !haskey(config_dict, "effector_params")
            config_dict["effector_params"] = Dict{String,Any}()
        end
        config_dict["effector_params"][key] = value
    end
    
    return config_dict
end

"""
    build_ppo_config(training_config::Dict{String,Any}) -> PPOConfig

Build PPOConfig from training configuration.
"""
function build_ppo_config(training_config::Dict{String,Any})
    return PPOConfig(
        max_episodes = Int(get(training_config, "max_episodes", 5000)),
        max_steps_per_episode = Int(get(training_config, "max_steps_per_episode", 100)),
        learning_rate = Float64(get(training_config, "learning_rate", 3e-4)),
        hidden_size = Int(get(training_config, "hidden_size", 32)),
        embed_size = Int(get(training_config, "embed_size", 16)),
        gamma = Float64(get(training_config, "gamma", 0.99)),
        gae_lambda = Float64(get(training_config, "gae_lambda", 0.95)),
        clip_epsilon = Float64(get(training_config, "clip_epsilon", 0.2)),
        value_loss_coef = Float64(get(training_config, "value_loss_coef", 0.5)),
        entropy_coef = Float64(get(training_config, "entropy_coef", 0.01)),
        max_grad_norm = Float64(get(training_config, "max_grad_norm", 0.5)),
        ppo_epochs = Int(get(training_config, "ppo_epochs", 10)),
        batch_size = Int(get(training_config, "batch_size", 64)),
    )
end

"""
    train_campaign(scenarios::Vector{TrainingScenario}, output_dir::String; parallel::Bool=true) -> Vector{TrainingResult}

Train multiple scenarios in a campaign.
"""
function train_campaign(scenarios::Vector{TrainingScenario}, output_dir::String; parallel::Bool=true)
    mkpath(output_dir)
    
    if parallel && length(scenarios) > 1
        return train_campaign_parallel(scenarios, output_dir)
    else
        return train_campaign_sequential(scenarios, output_dir)
    end
end

"""
    train_campaign_sequential(scenarios::Vector{TrainingScenario}, output_dir::String) -> Vector{TrainingResult}

Train scenarios sequentially.
"""
function train_campaign_sequential(scenarios::Vector{TrainingScenario}, output_dir::String)
    results = Vector{TrainingResult}()
    
    for scenario in scenarios
        @info "Training scenario: $(scenario.scenario_id)"
        result = train_scenario(scenario, output_dir)
        push!(results, result)
        
        if !result.success
            @warn "Scenario $(scenario.scenario_id) failed: $(result.error_message)"
        end
    end
    
    return results
end

"""
    train_campaign_parallel(scenarios::Vector{TrainingScenario}, output_dir::String) -> Vector{TrainingResult}

Train scenarios in parallel using Distributed.
"""
function train_campaign_parallel(scenarios::Vector{TrainingScenario}, output_dir::String)
    # Check if Distributed is available
    try
        using Distributed
    catch
        @warn "Distributed not available, falling back to sequential training"
        return train_campaign_sequential(scenarios, output_dir)
    end
    
    results = Vector{TrainingResult}()
    
    # Add workers if needed
    n_workers = min(length(scenarios), nprocs() - 1)
    if n_workers < 1
        n_workers = 1
    end
    
    @sync begin
        for scenario in scenarios
            @async begin
                result = train_scenario(scenario, output_dir)
                push!(results, result)
            end
        end
    end
    
    return results
end

"""
    generate_campaign_summary(results::Vector{TrainingResult}, output_path::String)

Generate a summary report of the training campaign.
"""
function generate_campaign_summary(results::Vector{TrainingResult}, output_path::String)
    mkpath(dirname(output_path))
    
    successful = filter(r -> r.success, results)
    failed = filter(r -> !r.success, results)
    
    summary = Dict{String,Any}(
        "total_scenarios" => length(results),
        "successful" => length(successful),
        "failed" => length(failed),
        "total_training_time" => sum(r.training_time for r in results),
        "average_training_time" => mean(r.training_time for r in results),
        "best_reward" => maximum(r.best_reward for r in successful; init=-Inf),
        "worst_cost" => minimum(r.final_cost for r in successful; init=Inf),
        "results" => [
            Dict{String,Any}(
                "scenario_id" => r.scenario_id,
                "success" => r.success,
                "final_episode" => r.final_episode,
                "best_reward" => r.best_reward,
                "final_cost" => r.final_cost,
                "training_time" => r.training_time,
                "error_message" => r.error_message,
            ) for r in results
        ],
    )
    
    # Save summary
    using YAML
    YAML.write_file(output_path, summary)
    
    return summary
end

export TrainingResult, train_scenario, train_campaign, generate_campaign_summary

end # module TrainingOrchestrator
