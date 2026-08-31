module ControlHooks
    using ..Structure

    using ..ConfigTypes: ODEParams
    using ..AbstractTypes: AbstractPlanet, AbstractControlEffectorModel, AbstractThrusterModel
    using ..RobotArmPlanning: RobotArmPlan, robot_arm_plan_sample
    using ..ThrusterModels: BaseThrusterModel, SixAxisThrusterModel
    using ..CommandTypes: PropulsiveManeuverCommand, PropulsiveBurnPlan
    using ..GuidanceModels: RPOPlan, RPOPlanBuffer
    using ..GravityEffectors: aerobraking_gravity_force_ii
    using ..AerodynamicEffectors: aerodynamic_coefficient_fM
    using ..GuidanceHooks: AerobrakingGuidanceInput, dispatch_aerobraking_guidance
    using ..GuidanceHooks: AerobrakingEnergyDepletionConfig, AerobrakingEnergyDepletionState
    using ..AerobrakingPolicy: AerobrakingPolicyConfig, DefaultAerobrakingPolicySelector
    using ..EnvironmentModels: getDensity
    using ..EphemeridesModels: ephemerides_requires_spice, planet_frame_lpi
    using ..ReferenceSystems
    using ..LinearAlgebra
    using ..StaticArrays
    using ..Kinematics
    using SparseArrays
    using OSQP

    const config = Structure

    export calcControlForceTorque, calcControlEffect!, calcControlMassFlowRate, calcReactionWheelTorque
    export AerobrakingEnergyDepletionControlModel, SolarPanelAngleOfAttackControlModel
    export AerobrakingMPCMode, TargetEnergyMode, MaxEnergyDepletionMode
    export AerobrakingMPCParams, AerobrakingMPCConfig, AerobrakingMPCProblem
    export AerobrakingMPCSolution, AerobrakingMPCState, AerobrakingMPCConstraintSet
    export AerobrakingMPCReferenceConfig, AerobrakingMPCControlModel
    export mpc_constraints, constraint_active, constraint_names, apply_constraints
    export mpc_params_from_spaceagora, mpc_prediction_gravity_model
    export spacecraft_mass_kg, spacecraft_reference_areas, mpc_config_from_spaceagora
    export density_and_gradient_from_spaceagora, density_function_from_spaceagora
    export build_reference_drag_pass, build_mpc_problem, solve_mpc_qp
    export objective_kind, objective_label, commanded_area_fraction
    export alpha_from_commanded_area, commanded_area_from_alpha, apply_commanded_area!
    export mpc_control_save_fields
    export RpoLQMPCController, init_rpo_lqmpc, rpo_lqmpc_control
    export RPOHeldActuation, RPOMPCControlModel
    export MagneticMomentumManagerModel
    export RobotArmHeldActuation, RobotArmJointMPCController, RobotArmControlEffector
    export init_robot_arm_joint_mpc, robot_arm_joint_mpc_reference_preview
    export robot_arm_joint_mpc_control, robot_arm_measured_joint_state
    export rpo_allocate_six_axis_thrusters, rpo_thruster_wrench_body

    include(joinpath(@__DIR__, "..", "internal", "bridge_helpers.jl"))
    include(joinpath(@__DIR__, "..", "..", "core", "numerics", "quaternion_utils.jl"))
    include(joinpath(@__DIR__, "propulsive_maneuvers.jl"))
    include(joinpath(@__DIR__, "heat_rate_control.jl"))
    include(joinpath(@__DIR__, "heat_load_control.jl"))
    include(joinpath(@__DIR__, "struct_load_control.jl"))
    include(joinpath(@__DIR__, "targeting_control.jl"))
    include(joinpath(@__DIR__, "rpo_mpc", "lqmpc.jl"))
    include(joinpath(@__DIR__, "rpo_mpc", "rpo_control_types.jl"))
    include(joinpath(@__DIR__, "rpo_mpc", "thruster_allocator.jl"))
    include(joinpath(@__DIR__, "rpo_mpc", "reaction_wheel_allocator.jl"))
    include(joinpath(@__DIR__, "rpo_mpc", "rpo_mpc_control_model.jl"))
    include(joinpath(@__DIR__, "robot_arm_control.jl"))
    include(joinpath(@__DIR__, "momentum_manager.jl"))
    include(joinpath(@__DIR__, "aerobraking_mpc", "aerobraking_mpc.jl"))
    include(joinpath(@__DIR__, "aerobraking", "control_commands.jl"))
    include(joinpath(@__DIR__, "aerobraking", "constraint_tracking.jl"))
    include(joinpath(@__DIR__, "aerobraking", "tracking_executor.jl"))
end
