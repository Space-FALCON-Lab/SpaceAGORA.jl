#!/usr/bin/env julia

using SpaceAGORA_RL

function main(args=ARGS)
    config_path = isempty(args) ? default_config_path() : args[1]
    config = resolve_config(config_path)
    episodes = length(args) >= 2 ?
               parse(Int, args[2]) :
               PAPER_GENERALIZATION_EVALUATION_EPISODES
    scenario = paper_evaluation_scenario(
        config.scenario;
        max_passes=max(1000, config.scenario.termination_config.max_passes),
    )
    outputs = Dict{String,Any}()
    for (name, suite_scenario) in generalization_suite_configs(scenario)
        outputs[name] = evaluate_baselines(
            suite_scenario;
            episodes=episodes,
            seed=config.training.seed,
        )
    end
    println("completed suites: ", join(sort(collect(keys(outputs))), ", "))
    return outputs
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
