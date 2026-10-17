"""Cached RPO control command held between control updates."""
Base.@kwdef mutable struct RPOHeldActuation
    force_ii::SVector{3, Float64} = SVector{3, Float64}(0.0, 0.0, 0.0)
    torque_body::SVector{3, Float64} = SVector{3, Float64}(0.0, 0.0, 0.0)
    thruster_forces_n::SVector{6, Float64} = SVector{6, Float64}(zeros(6))
    rw_torque_body::SVector{3, Float64} = SVector{3, Float64}(0.0, 0.0, 0.0)
end

"""
Optional per-update record of an `RPOMPCControlModel`: update time, chaser
relative state in RTN, the LQ-MPC acceleration command (RTN), the QP status,
the desired body force, the allocated thruster forces, the chaser mass and
attitude. It is filled only when attached to the model; read it after
`run_simulation(args; isolate_state=false)`, since the default isolated run
advances a copy of the model.
"""
Base.@kwdef mutable struct RPOControlCommandLog
    t_s::Vector{Float64} = Float64[]
    x_rel_rtn::Vector{SVector{6, Float64}} = SVector{6, Float64}[]
    accel_cmd_rtn::Vector{SVector{3, Float64}} = SVector{3, Float64}[]
    qp_status::Vector{Symbol} = Symbol[]
    force_body_desired_n::Vector{SVector{3, Float64}} = SVector{3, Float64}[]
    thruster_forces_n::Vector{SVector{6, Float64}} = SVector{6, Float64}[]
    mass_kg::Vector{Float64} = Float64[]
    q_chaser::Vector{SVector{4, Float64}} = SVector{4, Float64}[]
end

"""
RPO control effector that tracks guidance references with LQ-MPC and allocates
actuators. `chaser_index` and `target_index` are the run positions of the chaser
and the target; `chaser_idx` and `target_idx`, their old names, still work with a
deprecation warning until SpaceAGORA 0.5.0.

`actuation = :ideal` (the default) applies the allocated force, rotated to the
inertial frame at the control tick and held there, with the jets' and the
wheels' torque. `actuation = :mounted` drives the chaser's mounted thrusters at
the six allocated forces instead, so the thrust turns with the vehicle between
ticks: the forces go to the first six mounted thrusters, so mount
`six_axis_thrusters(model.thrusters)` on the chaser first; any other mounted
actuator stays idle. The wheel torque is applied as in `:ideal`.

Each run starts the held command from its value before the model's first run,
the controller's warm start from zero with a fresh solver workspace, and an
attached command log empty (`initialize_control!`).
"""
mutable struct RPOMPCControlModel <: AbstractControlEffectorModel
    chaser_index::Int
    target_index::Int
    thrusters::SixAxisThrusterModel
    minimum_throttle_fraction::Float64
    controller::Any
    plan_buffer::RPOPlanBuffer
    held::RPOHeldActuation
    control_dt_s::Float64
    attitude_kp::Float64
    rate_kd::Float64
    max_rw_torque_nm::Float64
    command_log::Union{Nothing, RPOControlCommandLog}
    actuation::Symbol
    # `:mounted` actuation: one level per actuator mounted on the chaser, refreshed each tick.
    mounted_levels::Vector{Float64}
    run_start::Union{Nothing, RPOHeldActuation}   # the held command before the first run
end
# The deprecated names forward, as keywords and as properties.
function RPOMPCControlModel(; chaser_index=_NOT_PASSED, target_index=_NOT_PASSED, chaser_idx=_NOT_PASSED, target_idx=_NOT_PASSED,
        thrusters=SixAxisThrusterModel(), minimum_throttle_fraction=0.0, controller=nothing,
        plan_buffer=RPOPlanBuffer(), held=RPOHeldActuation(), control_dt_s=1.0, attitude_kp=0.0, rate_kd=0.0,
        max_rw_torque_nm=Inf, command_log=nothing, actuation=:ideal, mounted_levels=Float64[])
    T = RPOMPCControlModel
    return T(_renamed_keyword(T, :chaser_idx, chaser_idx, :chaser_index, chaser_index, 1),
        _renamed_keyword(T, :target_idx, target_idx, :target_index, target_index, 2),
        thrusters, minimum_throttle_fraction, controller, plan_buffer, held, control_dt_s, attitude_kp, rate_kd,
        max_rw_torque_nm, command_log, actuation, mounted_levels, nothing)
end
@inline _rpo_mpc_field(name::Symbol) =
    name === :chaser_idx ? _renamed_field_warning(RPOMPCControlModel, :chaser_idx, :chaser_index) :
    name === :target_idx ? _renamed_field_warning(RPOMPCControlModel, :target_idx, :target_index) : name
@inline Base.getproperty(m::RPOMPCControlModel, name::Symbol) = getfield(m, _rpo_mpc_field(name))
@inline function Base.setproperty!(m::RPOMPCControlModel, name::Symbol, x)
    field = _rpo_mpc_field(name)
    return setfield!(m, field, convert(fieldtype(RPOMPCControlModel, field), x))
end
