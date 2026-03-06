@inline _checkpoint_directory(args) = SimulationModel.IOConfig._checkpoint_directory(args)
@inline _checkpoint_paths(args) = SimulationModel.IOConfig._checkpoint_paths(args)

@inline function _write_checkpoint!(args, t::Float64, u_state)
    return SimulationModel.IOSerialization._write_checkpoint!(args, t, u_state, CHECKPOINT_SCHEMA_VERSION)
end

@inline _load_checkpoint(args) = SimulationModel.IOSerialization._load_checkpoint(args)
@inline _clear_checkpoint!(args) = SimulationModel.IOSerialization._clear_checkpoint!(args)
