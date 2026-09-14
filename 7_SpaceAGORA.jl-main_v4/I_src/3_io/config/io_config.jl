module IOConfig

using Dates
using Random

@inline function _results_bundle_prefix(args)::String
    return joinpath(args.simulation_settings.results_directory, "simulation_results")
end

@inline function _results_csv_path(args)::String
    return joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
end

@inline function _collision_results_csv_path(args)::String
    stamp = Dates.format(now(UTC), dateformat"yyyymmddTHHMMSSsss")
    token = string(stamp, ".", getpid(), ".", rand(UInt))
    return joinpath(args.simulation_settings.results_directory, "simulation_results.$token.csv")
end

@inline function _checkpoint_directory(args)::String
    if isempty(args.simulation_settings.checkpoint_directory)
        return joinpath(args.simulation_settings.results_directory, "checkpoints")
    end
    return args.simulation_settings.checkpoint_directory
end

@inline function _checkpoint_paths(args)
    ckpt_dir = _checkpoint_directory(args)
    return (
        data=joinpath(ckpt_dir, "simulation_checkpoint.bin"),
        manifest=joinpath(ckpt_dir, "simulation_checkpoint.manifest.toml")
    )
end

export _results_bundle_prefix
export _results_csv_path
export _collision_results_csv_path
export _checkpoint_directory
export _checkpoint_paths

end # module IOConfig
