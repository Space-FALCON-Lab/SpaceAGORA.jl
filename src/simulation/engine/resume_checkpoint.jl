@inline _checkpoint_directory(args) = SimulationModel.IOConfig._checkpoint_directory(args)
@inline _checkpoint_paths(args) = SimulationModel.IOConfig._checkpoint_paths(args)

function _base_thruster_runtime_snapshot(control_model::SimulationModel.BaseThrusterModel)
    return (
        kind=:BaseThrusterModel,
        Δv=copy(control_model.Δv),
        direction=copy(control_model.direction),
        start_burn_time=copy(control_model.start_burn_time),
        stop_burn_time=copy(control_model.stop_burn_time)
    )
end

function _checkpoint_runtime_state(p)
    control_effectors = Any[]
    for control_model in p.args.control_model.control_effectors
        if control_model isa SimulationModel.BaseThrusterModel
            push!(control_effectors, _base_thruster_runtime_snapshot(control_model))
        else
            push!(control_effectors, nothing)
        end
    end
    return (
        orbit_counter=copy(p.orbit_counter),
        is_active=copy(p.is_active),
        maneuver_commands=deepcopy(p.shared_buffers.maneuver_commands),
        control_effectors=control_effectors
    )
end

function _restore_vector_field!(
    dest::AbstractVector{T},
    src,
    field_name::Symbol,
    effector_idx::Int
) where {T}
    length(dest) == length(src) || throw(ArgumentError(
        "Checkpoint runtime field $(field_name) for control effector $(effector_idx) has length $(length(src)); expected $(length(dest))."
    ))
    dest .= T.(src)
    return nothing
end

function _restore_base_thruster_runtime!(
    control_model::SimulationModel.BaseThrusterModel,
    snapshot,
    effector_idx::Int
)
    snapshot === nothing && return nothing
    hasproperty(snapshot, :kind) && snapshot.kind == :BaseThrusterModel || throw(ArgumentError(
        "Checkpoint runtime field for control effector $(effector_idx) is not a BaseThrusterModel snapshot."
    ))
    _restore_vector_field!(control_model.Δv, snapshot.Δv, :Δv, effector_idx)
    _restore_vector_field!(control_model.direction, snapshot.direction, :direction, effector_idx)
    _restore_vector_field!(control_model.start_burn_time, snapshot.start_burn_time, :start_burn_time, effector_idx)
    _restore_vector_field!(control_model.stop_burn_time, snapshot.stop_burn_time, :stop_burn_time, effector_idx)
    return nothing
end

function _restore_checkpoint_runtime_state!(p, runtime_state)
    runtime_state === nothing && return false

    if hasproperty(runtime_state, :orbit_counter)
        length(p.orbit_counter) == length(runtime_state.orbit_counter) || throw(ArgumentError(
            "Checkpoint orbit_counter length $(length(runtime_state.orbit_counter)) does not match simulation length $(length(p.orbit_counter))."
        ))
        p.orbit_counter .= Int64.(runtime_state.orbit_counter)
    end
    if hasproperty(runtime_state, :is_active)
        length(p.is_active) == length(runtime_state.is_active) || throw(ArgumentError(
            "Checkpoint is_active length $(length(runtime_state.is_active)) does not match simulation length $(length(p.is_active))."
        ))
        p.is_active .= Bool.(runtime_state.is_active)
    end
    if hasproperty(runtime_state, :maneuver_commands)
        length(p.shared_buffers.maneuver_commands) == length(runtime_state.maneuver_commands) || throw(ArgumentError(
            "Checkpoint maneuver_commands length $(length(runtime_state.maneuver_commands)) does not match simulation length $(length(p.shared_buffers.maneuver_commands))."
        ))
        p.shared_buffers.maneuver_commands .= deepcopy(runtime_state.maneuver_commands)
    end
    if hasproperty(runtime_state, :control_effectors)
        length(p.args.control_model.control_effectors) == length(runtime_state.control_effectors) || throw(ArgumentError(
            "Checkpoint control_effectors length $(length(runtime_state.control_effectors)) does not match simulation length $(length(p.args.control_model.control_effectors))."
        ))
        for (idx, control_model) in pairs(p.args.control_model.control_effectors)
            if control_model isa SimulationModel.BaseThrusterModel
                _restore_base_thruster_runtime!(control_model, runtime_state.control_effectors[idx], idx)
            end
        end
    end
    return true
end

@inline function _write_checkpoint!(args, t::Float64, u_state; runtime_state=nothing)
    return SimulationModel.IOSerialization._write_checkpoint!(
        args,
        t,
        u_state,
        CHECKPOINT_SCHEMA_VERSION;
        runtime_state=runtime_state
    )
end

@inline function _write_checkpoint!(args, t::Float64, u_state, p)
    return _write_checkpoint!(args, t, u_state; runtime_state=_checkpoint_runtime_state(p))
end

function _rewrite_checkpoint_state!(args, t::Float64, u_state)
    ckpt = _load_checkpoint(args)
    ckpt === nothing && throw(ArgumentError("Cannot rewrite checkpoint state because no checkpoint was found."))
    return _write_checkpoint!(args, t, u_state; runtime_state=ckpt.runtime_state)
end

@inline _load_checkpoint(args) = SimulationModel.IOSerialization._load_checkpoint(args)
@inline _clear_checkpoint!(args) = SimulationModel.IOSerialization._clear_checkpoint!(args)
