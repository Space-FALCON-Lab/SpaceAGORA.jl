#!/usr/bin/env julia

using SpaceAGORA_RL

function main(args=ARGS)
    config_path = isempty(args) ? default_config_path() : args[1]
    config = resolve_config(config_path)
    results = evaluate_baselines(config.scenario; episodes=config.training.episodes, seed=config.training.seed)
    paths = write_evaluation_artifacts(config.reports.output_dir, results)
    println("wrote baseline CSV artifacts: ", paths)
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
