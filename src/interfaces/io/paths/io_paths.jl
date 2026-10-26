module IOConfig

using Dates
using Random

@inline function _results_bundle_prefix(args)::String
    return joinpath(args.simulation_settings.results_directory, "simulation_results")
end

@inline function _results_csv_path(args)::String
    return joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
end

# The link schedule of a run with interlinks, written beside its results.
@inline function _link_schedule_paths(args)
    prefix = _results_bundle_prefix(args)
    return (csv=prefix * "_link_schedule.csv", feather=prefix * "_link_schedule.feather")
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
        manifest=joinpath(ckpt_dir, "simulation_checkpoint.manifest.toml"),
        rows=joinpath(ckpt_dir, "simulation_checkpoint.results.bin"),
        interlinks=joinpath(ckpt_dir, "simulation_checkpoint.interlinks.bin"),
        # A long-layout run's streamed results tables (Arrow IPC streams).
        stream_spacecraft=joinpath(ckpt_dir, "simulation_checkpoint.results_spacecraft.arrows"),
        stream_run=joinpath(ckpt_dir, "simulation_checkpoint.results_run.arrows")
    )
end

export _results_bundle_prefix
export _results_csv_path
export _collision_results_csv_path
export _link_schedule_paths
export _checkpoint_directory
export _checkpoint_paths

end # module IOConfig
