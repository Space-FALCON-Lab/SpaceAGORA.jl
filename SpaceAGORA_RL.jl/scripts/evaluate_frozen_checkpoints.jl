#!/usr/bin/env julia

using SpaceAGORA_RL

function main(args=ARGS)
    length(args) >= 2 || throw(ArgumentError(
        "usage: evaluate_frozen_checkpoints.jl CONFIG_PATH CHECKPOINT_DIR " *
        "[OUTPUT_DIR] [EPISODES]",
    ))
    config = resolve_config(args[1])
    checkpoint_dir = args[2]
    output_dir = length(args) >= 3 ?
                 args[3] :
                 joinpath(config.reports.output_dir, "frozen_checkpoints")
    episodes = length(args) >= 4 ?
               parse(Int, args[4]) :
               PAPER_IID_EVALUATION_EPISODES
    scenario = paper_evaluation_scenario(
        config.scenario;
        max_passes=max(1000, config.scenario.termination_config.max_passes),
    )
    results = evaluate_frozen_checkpoints(
        checkpoint_dir,
        scenario;
        episodes=episodes,
        seed=config.training.seed,
        output_dir=output_dir,
    )
    for checkpoint_path in sort(collect(keys(results)))
        for mode in PAPER_EVALUATION_MODES
            println(basename(checkpoint_path), " ", mode, ": ",
                    results[checkpoint_path][mode].aggregate)
        end
    end
    println("wrote frozen-checkpoint evaluation artifacts to ", output_dir)
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
