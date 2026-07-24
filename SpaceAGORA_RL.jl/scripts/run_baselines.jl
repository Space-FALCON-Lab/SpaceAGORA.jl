#!/usr/bin/env julia

using SpaceAGORA_RL

function main(args=ARGS)
    config_path = isempty(args) ? default_config_path() : args[1]
    config = resolve_config(config_path)
    episodes = length(args) >= 2 ?
               parse(Int, args[2]) :
               PAPER_IID_EVALUATION_EPISODES
    scenario = paper_evaluation_scenario(
        config.scenario;
        max_passes=max(250, config.scenario.termination_config.max_passes),
    )
    results = evaluate_baselines(scenario; episodes=episodes, seed=config.training.seed)
    paths = write_evaluation_artifacts(config.reports.output_dir, results)
    println("wrote baseline CSV artifacts: ", paths)
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
