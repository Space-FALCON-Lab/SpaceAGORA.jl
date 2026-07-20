module ControlHooks
    using ..Structure

    using ..ConfigTypes: ODEParams
    using ..AbstractTypes: AbstractPlanet, AbstractControlEffectorModel, AbstractThrusterModel
    using ..RobotArmPlanning: RobotArmPlan, robot_arm_plan_sample
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

    export calcControlForceTorque, calcControlEffect!, calcControlMassFlowRate, calcReactionWheelTorque
    export RpoLQMPCController, init_rpo_lqmpc, rpo_lqmpc_control
    export RPOHeldActuation, RPOMPCControlModel
    export RobotArmHeldActuation, RobotArmJointMPCController, RobotArmControlEffector
    export init_robot_arm_joint_mpc, robot_arm_joint_mpc_reference_preview
    export robot_arm_joint_mpc_control, robot_arm_measured_joint_state
    export rpo_allocate_six_axis_thrusters, rpo_thruster_wrench_body

    include(joinpath(@__DIR__, "..", "internal", "bridge_helpers.jl"))
    include(joinpath(@__DIR__, "..", "..", "core", "numerics", "quaternion_utils.jl"))
    include(joinpath(@__DIR__, "propulsive_maneuvers.jl"))
    include(joinpath(@__DIR__, "rpo_mpc", "lqmpc.jl"))
    include(joinpath(@__DIR__, "rpo_mpc", "rpo_control_types.jl"))
    include(joinpath(@__DIR__, "rpo_mpc", "thruster_allocator.jl"))
    include(joinpath(@__DIR__, "rpo_mpc", "reaction_wheel_allocator.jl"))
    include(joinpath(@__DIR__, "rpo_mpc", "rpo_mpc_control_model.jl"))
    include(joinpath(@__DIR__, "robot_arm_control.jl"))
    include(joinpath(@__DIR__, "aerobraking", "control_commands.jl"))
    include(joinpath(@__DIR__, "aerobraking", "constraint_tracking.jl"))
    include(joinpath(@__DIR__, "aerobraking", "tracking_executor.jl"))
end
