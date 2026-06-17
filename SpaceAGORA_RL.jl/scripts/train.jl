#!/usr/bin/env julia

using SpaceAGORA_RL

function main(args=ARGS)
    config_path = isempty(args) ? default_config_path() : args[1]
    config = resolve_config(config_path)
    session = build_training_session(config)
    result = train_parallel!(session)
    println("wrote run artifacts to ", result.output_dir)
    println("episodes: ", length(result.summaries),
            " transitions: ", length(result.transitions),
            " global_step: ", result.global_step)
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
