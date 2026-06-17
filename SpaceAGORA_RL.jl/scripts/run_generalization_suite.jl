#!/usr/bin/env julia

using SpaceAGORA_RL

function main(args=ARGS)
    config_path = isempty(args) ? default_config_path() : args[1]
    config = resolve_config(config_path)
    outputs = Dict{String,Any}()
    for (name, scenario) in generalization_suite_configs(config.scenario)
        outputs[name] = evaluate_baselines(scenario; episodes=config.training.episodes, seed=config.training.seed)
    end
    println("completed suites: ", join(sort(collect(keys(outputs))), ", "))
    return outputs
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
