module ControlHooks
    using ..Structure

    using ..ConfigTypes: ODEParams
    using ..AbstractTypes: AbstractPlanet, AbstractControlEffectorModel, AbstractThrusterModel
    using ..ThrusterModels: BaseThrusterModel, SixAxisThrusterModel
    using ..CommandTypes: PropulsiveManeuverCommand, PropulsiveBurnPlan
    using ..GuidanceModels: RPOPlan, RPOPlanBuffer
    using ..GravityEffectors: aerobraking_gravity_force_ii
    using ..GuidanceHooks: AerobrakingGuidanceInput, dispatch_aerobraking_guidance
    using ..AerobrakingPolicy: AerobrakingPolicyConfig, DefaultAerobrakingPolicySelector
    using ..ReferenceSystems
    using ..LinearAlgebra
    using ..StaticArrays
    using ..Kinematics
    using SparseArrays
    using OSQP

    const config = Structure

    export calcControlForceTorque, calcControlEffect!, calcControlMassFlowRate
    export RpoLQMPCController, init_rpo_lqmpc, rpo_lqmpc_control
    export RPOHeldActuation, RPOMPCControlModel
    export rpo_allocate_six_axis_thrusters, rpo_thruster_wrench_body

    include(joinpath(@__DIR__, "..", "internal", "bridge_helpers.jl"))
    include(joinpath(@__DIR__, "..", "..", "core", "numerics", "quaternion_utils.jl"))
    include(joinpath(@__DIR__, "propulsive_maneuvers.jl"))
    include(joinpath(@__DIR__, "rpo_mpc", "lqmpc.jl"))
    include(joinpath(@__DIR__, "rpo_mpc", "rpo_control_types.jl"))
    include(joinpath(@__DIR__, "rpo_mpc", "thruster_allocator.jl"))
    include(joinpath(@__DIR__, "rpo_mpc", "reaction_wheel_allocator.jl"))
    include(joinpath(@__DIR__, "rpo_mpc", "rpo_mpc_control_model.jl"))
    include(joinpath(@__DIR__, "aerobraking", "control_commands.jl"))
    include(joinpath(@__DIR__, "aerobraking", "constraint_tracking.jl"))
    include(joinpath(@__DIR__, "aerobraking", "tracking_executor.jl"))
end
