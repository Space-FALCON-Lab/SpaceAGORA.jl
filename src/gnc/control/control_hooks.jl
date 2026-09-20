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
    using ..GuidanceHooks: ApolloDescentConfig, ApolloDescentState, _validate_descent_config, _descent_indices
    using ..AbstractTypes: AbstractTerrainModel
    using ..TerrainModels: NoTerrainModel, DEMTerrainModel
    using ..AerobrakingPolicy: AerobrakingPolicyConfig, DefaultAerobrakingPolicySelector
    using ..EnvironmentModels: getDensity
    using ..EphemeridesModels: ephemerides_requires_spice, planet_frame_lpi
    using ..ReferenceSystems
    using ..LinearAlgebra
    using ..StaticArrays
    using ..Kinematics
    using SparseArrays
    using AstroTime
    using SPICE
    using OSQP

    const config = Structure

    export calcControlForceTorque, calcControlEffect!, calcControlMassFlowRate, calcReactionWheelTorque
    export control_thruster_levels, touchdown_spec
    export AerobrakingEnergyDepletionControlModel, SolarPanelAngleOfAttackControlModel
    export ApolloDescentControlConfig, ApolloDescentControlModel, ApolloDescentControlState, attitude_error_vector
    export RpoLQMPCController, init_rpo_lqmpc, rpo_lqmpc_control
    export RPOHeldActuation, RPOMPCControlModel
    export MagneticMomentumManagerModel
    export RobotArmHeldActuation, RobotArmJointMPCController, RobotArmControlEffector
    export init_robot_arm_joint_mpc, robot_arm_joint_mpc_reference_preview
    export robot_arm_joint_mpc_control, robot_arm_measured_joint_state
    export rpo_allocate_six_axis_thrusters, rpo_thruster_wrench_body

    """
        control_thruster_levels(effector, i::Int) -> Union{Nothing, AbstractVector{<:Real}}

    Optional hook: the firing level (0 to 1) of every thruster of spacecraft `i`,
    in the order the visualization scene lists them (the spacecraft's links in
    order, each link's `thrusters` in order). Effectors that drive no named
    thruster return `nothing`, the default for any effector type that does not
    override this method; an effector that drives only some of them may return a
    shorter vector, and the remaining thrusters are idle.

    The engine reads this at save time, so an effector must compute the levels
    in `calcControlEffect!` and keep them in its own actuator state rather than
    recomputing them here.
    """
    control_thruster_levels(effector, i::Int) = nothing

    """
        touchdown_spec(effector, spacecraft_index::Int) -> Union{Nothing, NamedTuple}

    Optional terrain-contact event for a selected spacecraft's index in the run
    (not its user-defined ID). Return `nothing` for unselected spacecraft, or
    `(terrain, reference_radius_m, height_m, on_touchdown)`: an AbstractTerrainModel,
    its explicit positive reference-sphere radius in metres, a finite nonnegative
    clearance in metres, and a callable accepting `(t, r_p, v_p, spacecraft_index)`.
    The callback receives planet-fixed position and ground-relative velocity in
    metres and metres per second. It may record contact and turn off actuators.

    The first downward crossing deactivates only this spacecraft; the solve ends
    once all spacecraft are inactive. Initial contact/below-ground states and
    multiple specifications for one index are rejected. The usual 50 km impact
    stop remains active for every spacecraft with no touchdown specification.
    """
    touchdown_spec(effector, spacecraft_index::Int) = nothing

    include(joinpath(@__DIR__, "..", "internal", "bridge_helpers.jl"))
    using ..QuaternionMath
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
    include(joinpath(@__DIR__, "aerobraking", "control_commands.jl"))
    include(joinpath(@__DIR__, "aerobraking", "constraint_tracking.jl"))
    include(joinpath(@__DIR__, "aerobraking", "tracking_executor.jl"))
    include(joinpath(@__DIR__, "landing", "apollo_descent_control.jl"))
end
