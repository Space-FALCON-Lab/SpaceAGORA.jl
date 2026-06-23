#!/usr/bin/env julia

"""
    run_rl_training_campaign.jl

Entry point script for running RL training campaigns.
Supports cluster combinations, parameter sweeps, and checkpoint/resume functionality.
Compatible with existing PPOTrainer checkpoint system and TensorBoard logging.
"""

using SpaceAGORA.ConstellationDesign
using YAML
using ArgParse
using Logging
using Dates

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--config", "-c"
            help = "Path to campaign configuration YAML file"
            arg_type = String
            required = true
        "--output-dir", "-o"
            help = "Output directory for training results"
            arg_type = String
            default = "data/rl_campaigns"
        "--resume"
            help = "Resume from checkpoint (scenario_id or 'all')"
            arg_type = String
            default = ""
        "--dry-run"
            help = "Generate scenarios without training"
            action = :store_true
        "--parallel"
            help = "Enable parallel training"
            action = :store_true
            default = true
        "--verbose", "-v"
            help = "Enable verbose logging"
            action = :store_true
    end

    return parse_args(s)
end

function main()
    args = parse_commandline()
    
    # Setup logging
    if args["verbose"]
        logger = ConsoleLogger(stdout, Logging.Debug)
        global_logger(logger)
    else
        logger = ConsoleLogger(stdout, Logging.Info)
        global_logger(logger)
    end
    
    @info "Starting RL training campaign"
    @info "Config: $(args["config"])"
    @info "Output directory: $(args["output-dir"])"
    
    # Load campaign configuration
    campaign_config = YAML.load_file(args["config"])
    
    # Build training scenarios
    @info "Building training scenarios..."
    scenarios = build_training_scenarios(campaign_config)
    @info "Generated $(length(scenarios)) training scenarios"
    
    # Save scenarios for reference
    scenarios_file = joinpath(args["output-dir"], "training_scenarios.yaml")
    mkpath(args["output-dir"])
    save_training_scenarios(scenarios, scenarios_file)
    @info "Saved scenarios to $scenarios_file"
    
    # Dry run: just generate scenarios
    if args["dry-run"]
        @info "Dry run complete. Scenarios generated but no training performed."
        return
    end
    
    # Resume from checkpoint if specified
    if !isempty(args["resume"])
        @info "Resume mode: $(args["resume"])"
        scenarios = filter_scenarios_for_resume(scenarios, args["resume"], args["output-dir"])
        @info "Filtered to $(length(scenarios)) scenarios for resume"
    end
    
    # Run training campaign
    @info "Starting training campaign..."
    results = train_campaign(scenarios, args["output-dir"]; parallel=args["parallel"])
    
    # Generate summary
    summary_file = joinpath(args["output-dir"], "campaign_summary.yaml")
    summary = generate_campaign_summary(results, summary_file)
    @info "Campaign summary saved to $summary_file"
    
    # Print summary
    @info "Campaign complete"
    @info "Total scenarios: $(summary["total_scenarios"])"
    @info "Successful: $(summary["successful"])"
    @info "Failed: $(summary["failed"])"
    @info "Total training time: $(summary["total_training_time"])s"
    @info "Best reward: $(summary["best_reward"])"
    @info "Worst cost: $(summary["worst_cost"])"
end

"""
    filter_scenarios_for_resume(scenarios::Vector{TrainingScenario}, resume_spec::String, output_dir::String) -> Vector{TrainingScenario}

Filter scenarios based on resume specification.
"""
function filter_scenarios_for_resume(scenarios::Vector{TrainingScenario}, resume_spec::String, output_dir::String)
    if resume_spec == "all"
        # Resume all incomplete scenarios
        return filter(s -> !is_scenario_complete(s, output_dir), scenarios)
    else
        # Resume specific scenario
        return filter(s -> s.scenario_id == resume_spec, scenarios)
    end
end

"""
    is_scenario_complete(scenario::TrainingScenario, output_dir::String) -> Bool

Check if a scenario has been completed (has final model).
"""
function is_scenario_complete(scenario::TrainingScenario, output_dir::String)
    model_path = joinpath(output_dir, scenario.scenario_id, "final_model.jld2")
    return isfile(model_path)
end

# Run main function
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
