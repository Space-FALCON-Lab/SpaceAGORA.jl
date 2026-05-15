Base.@kwdef mutable struct RPOHeldActuation
    force_ii::SVector{3, Float64} = SVector{3, Float64}(0.0, 0.0, 0.0)
    torque_body::SVector{3, Float64} = SVector{3, Float64}(0.0, 0.0, 0.0)
    thruster_forces_n::SVector{6, Float64} = SVector{6, Float64}(zeros(6))
end

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
end
