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

"""RPO control effector that tracks guidance references with LQ-MPC and allocates actuators."""
Base.@kwdef mutable struct RPOMPCControlModel <: AbstractControlEffectorModel
    chaser_idx::Int = 1
    target_idx::Int = 2
    thrusters::SixAxisThrusterModel = SixAxisThrusterModel()
    controller::Any = nothing
    plan_buffer::RPOPlanBuffer = RPOPlanBuffer()
    held::RPOHeldActuation = RPOHeldActuation()
    control_dt_s::Float64 = 1.0
    attitude_kp::Float64 = 0.0
    rate_kd::Float64 = 0.0
    max_rw_torque_nm::Float64 = Inf
    command_log::Union{Nothing, RPOControlCommandLog} = nothing
end
