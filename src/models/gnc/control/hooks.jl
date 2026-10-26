module ControlHooks
    import ...SimulationLifecycle
    import ...GNCRuntime
    using ...GNCRuntime: BoundSpacecraft, CoupledSpacecraft
    using ..Structure

    using ..ConfigTypes: ODEParams
    using ..EffectorSampling: StateSample
    using ..AbstractTypes: AbstractPlanet, AbstractControlEffectorModel, AbstractThrusterModel
    using ..RobotArmPlanning: RobotArmPlan, RobotArmPlanBuffer, robot_arm_plan_sample
    using ..Robotics: ClothArmModel
    using ..ThrusterModels: BaseThrusterModel, SixAxisThrusterModel
    import ..Actuators: actuation_mode, actuator_levels
    using ..Actuators: propellant_mass_rate, allocate_actuators
    using ..ActuatorMounting: mounted_actuator_count, mounted_thruster_count, mounted_wheel_count, mounted_magnetorquer_count, actuator_context, move_hinge!
    using ..Components: thrust_direction, PanelHinge, ThrusterGimbal
    using ..ActuatorMounting: mounted_actuators
    using ..CommandTypes: PropulsiveManeuverCommand, PropulsiveBurnPlan
    using ..GuidanceModels: RPOPlan, RPOPlanBuffer
    using ..GravityEffectors: aerobraking_gravity_force_ii
    using ..AerodynamicEffectors: aerodynamic_coefficient_fM
    import ..GuidanceHooks: _accepted_campaign_request
    using ..GuidanceHooks: AerobrakingEnergyDepletionConfig, AerobrakingEnergyDepletionState, AerobrakingEnergyDepletionGuidanceModel
    using ..GuidanceHooks: ApolloDescentConfig, ApolloDescentState, _validate_descent_config, _descent_indices
    using ..AbstractTypes: AbstractTerrainModel
    using ..TerrainModels: NoTerrainModel, DEMTerrainModel
    using ..EnvironmentModels: getDensity, density_altitude_derivative, _environment_wind
    using ..DynamicsKS: cartesian_to_ks_state, ks_state_to_cartesian, ks_position, ks_velocity
    using ..DynamicsKS: ks_rk4_step, ks_first_order_tangent_map
    using ..DynamicsKS: ks_kinematics_jacobians, ks_density_value_gradient, ks_rotation_cross_matrix
    using ..EphemeridesModels: ephemerides_requires_spice
    using ..ReferenceSystems
    using ..LinearAlgebra
    using ..StaticArrays
    using ..Kinematics
    using SparseArrays
    using AstroTime
    using SPICE
    using OSQP

    const config = Structure

    export update_control!, control_wrench, control_propellant_rate, control_wheel_torque
    export control_thruster_levels, touchdown_spec, refresh_burn!, mounted_burn_tick!, burn_edge
    export AerobrakingEnergyDepletionControlModel, SolarPanelAngleOfAttackControlModel
    export ApolloDescentControlConfig, ApolloDescentControlModel, ApolloDescentControlState, attitude_error_vector
    export RpoLQMPCController, init_rpo_lqmpc, rpo_lqmpc_control
    export RPOControlCommandLog, RPOHeldActuation, RPOMPCControlModel
    export MagneticMomentumManagerModel
    export AerobrakingMPCMode, TargetEnergyMode, MaxEnergyDepletionMode
    export AerobrakingMPCParams, AerobrakingMPCConfig, AerobrakingMPCProblem
    export AerobrakingMPCSolution, AerobrakingMPCState, AerobrakingMPCConstraintSet
    export AerobrakingMPCReferenceConfig, AerobrakingMPCControlModel
    export AerobrakingMPCCampaignState, AerobrakingMPCCampaignControlModel
    export mpc_constraints, constraint_active, constraint_names, apply_constraints
    export mpc_params_from_spaceagora, mpc_config_from_spaceagora
    export spacecraft_mass_kg, spacecraft_reference_areas
    export density_function_from_spaceagora
    export build_reference_drag_pass, build_mpc_problem, solve_mpc_qp
    export objective_kind, objective_label, commanded_area_fraction
    export alpha_from_commanded_area, commanded_area_from_alpha
    export interpolate_mpc_plan, interpolate_mpc_history
    export evaluate_cartesian_mpc_outputs, evaluate_ks_mpc_outputs
    export cumulative_mpc_heat_load, propagate_ks_mpc_plan
    export mpc_control_save_fields, mpc_campaign_save_fields
    export RobotArmHeldActuation, RobotArmJointMPCController, RobotArmControlEffector
    export init_robot_arm_joint_mpc, robot_arm_joint_mpc_reference_preview
    export robot_arm_joint_mpc_control, robot_arm_measured_joint_state
    export rpo_allocate_six_axis_thrusters, rpo_thruster_wrench_body

    """
        control_thruster_levels(effector, spacecraft_index::Int) -> Union{Nothing, AbstractVector{<:Real}}

    Optional hook: the firing level (0 to 1) of every thruster of the spacecraft
    at `spacecraft_index` (its position in the run, from 1), in the order the
    visualization scene lists them (the spacecraft's links in
    order, each link's `thrusters` in order). Effectors that drive no named
    thruster return `nothing`, the default for any effector type that does not
    override this method; an effector that drives only some of them may return a
    shorter vector, and the remaining thrusters are idle.

    The engine reads this at save time, so an effector must compute the levels
    in `update_control!` and keep them in its own actuator state rather than
    recomputing them here. During the run it asks only for the spacecraft the
    effector's `spacecraft_scope` names (for every spacecraft in SpaceAGORA 0.6
    and earlier). Before the run, when the result columns are sized, it asks
    every effector about every spacecraft with thrusters, so keep returning
    `nothing` for a spacecraft outside the scope.
    """
    control_thruster_levels(effector, spacecraft_index::Int) = nothing

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

    """
        control_wrench(model, x::StateSample, p, t::Float64) -> (force_ii, torque_body)

    The force and torque a control effector applies to its spacecraft, which the
    dynamics evaluate at every right-hand-side evaluation in `:ideal` actuation,
    the default `actuation_mode`. Every control effector whose actuation is not
    `:mounted` needs a method for its own type. There is no fallback:
    `run_simulation` checks for the method before integration starts, as it does
    for `update_control!`.

    Return a tuple of two 3-vectors:

    - `force_ii`: in the inertial frame, in N;
    - `torque_body`: in the bus frame, about the root body's centre of mass, in
      N·m. The dynamics apply it only when the run propagates attitude.

    An effector that applies neither returns two zero vectors. The arguments:

    - `model`: the control effector.
    - `x`: the [`StateSample`](@ref) of that spacecraft. `x.pos_ii` (m) and
      `x.vel_ii` (m/s) are relative to the planet's centre in the inertial frame
      and `x.mass_kg` is its mass. With attitude propagation, `x.q_ib` is the
      quaternion `q = (x, y, z, w)`, scalar last, with
      `v_body = rot(q) * v_inertial` (body means bus axes), and `x.ω_body` the bus
      angular velocity in rad/s; both are `nothing` otherwise.
      `x.spacecraft_index` is the spacecraft's position in the run (`Int`, from
      1), not its ID.
    - `p`: the run's parameters.
    - `t`: elapsed time since the run's start epoch, in seconds (`Float64`).

    Return the command held since the last `update_control!` tick rather than
    recomputing the control law here. A method typed `x::StateSample`, or with
    `x` untyped, is the one the dynamics call. The dynamics call it, and the
    other output hooks, only for the spacecraft the model's `spacecraft_scope`
    names (SpaceAGORA 0.6 and earlier called them for every spacecraft); return
    zeros for any other, which direct callers may still ask about.
    """
    function control_wrench end

    """
        control_propellant_rate(model, x::StateSample, p, t::Float64) -> Float64

    The propellant a control effector uses, in kg/s: **positive** while
    propellant is used, and `0.0` when none is. The engine subtracts it from the
    spacecraft's mass derivative. A nonfinite rate is ignored.

    The dynamics call it at every right-hand-side evaluation in `:ideal`
    actuation, so return the rate held since the last `update_control!` tick.
    An effector that uses no propellant needs no method: the fallback, for any
    model, returns `0.0`. The arguments are those of [`control_wrench`](@ref); the
    spacecraft's index in the run is `x.spacecraft_index`. The
    `BaseThrusterModel` method returns `norm(force) / (Isp * g0)` during an
    active valid burn and `0.0` otherwise.
    """
    function control_propellant_rate end
    control_propellant_rate(model, x, p, t) = 0.0

    """
        control_wheel_torque(model, x::StateSample, p, t::Float64) -> Union{Nothing, AbstractVector}

    Optional hook for a control effector that drives reaction wheels: the part
    of its commanded torque on the spacecraft that the wheels produce, a 3-vector
    in the bus frame in N·m. An effector that drives no wheel needs no method:
    the fallback, for any model, returns `nothing`.

    The momentum stored in the wheels (the `h_wheels` state) changes at the
    opposite rate, mapped onto the wheel axes. In `:ideal` actuation the engine
    applies this torque only to the wheels: the effector includes it in the
    torque its `control_wrench` returns. In `:mounted` actuation the engine adds
    it to the torque of the mounted actuators.

    The dynamics call it at every right-hand-side evaluation, so return the
    torque held since the last `update_control!` tick. The arguments are
    those of [`control_wrench`](@ref); the spacecraft's index in the run is
    `x.spacecraft_index`.
    """
    function control_wheel_torque end
    control_wheel_torque(model, x, p, t) = nothing

    """
        update_control!(model, u, p, t, spacecraft_index) -> nothing

    The control tick, which every effector in a `ControlModel` implements for its
    own type. The effector computes its command (a thrust, a torque, actuator
    levels, panel angles) and keeps it in its own state; the engine ignores the
    return value and applies the held command until the next tick:

    - in `:ideal` actuation, the default `actuation_mode`, through the effector's
      `control_wrench`, `control_propellant_rate` and `control_wheel_torque`,
      which the dynamics evaluate at every step;
    - in `:mounted` actuation, through the levels `actuator_levels` returns for
      the actuators mounted on the spacecraft; `control_wrench` and
      `control_propellant_rate` are not called.

    The arguments:

    - `model`: the control effector.
    - `u`: the full integrated state. `u.sc[spacecraft_index]` is the spacecraft, with
      `pos` (m) and `vel` (m/s) relative to the planet's centre in the inertial
      frame and `mass` (kg); with attitude propagation also `q = (x, y, z, w)`,
      scalar last, with `v_body = rot(q) * v_inertial` (body means bus axes), and
      `ω`, the bus angular velocity in rad/s.
    - `p`: the run's parameters; `p.args` is the `SimulationConfiguration`.
    - `t`: elapsed time since the run's start epoch, in seconds (`Float64`).
    - `spacecraft_index`: the spacecraft's position in the run (`Int`, from 1), not its ID.

    Effector `k` runs on the simulation's control clock every `control_rates[k]`
    seconds, once for each spacecraft in its declared scope. Checkpoint segments
    and resume retain that clock rather than starting a new period. The calls for
    different spacecraft run one after another unless the effector's type is
    declared thread-safe (`SimulationModel.SimulationCallbacks.control_model_threadsafe(::MyEffector) = true`)
    or the run sets `SPACEAGORA_CONTROL_CALLBACK_ASSUME_THREADSAFE=1`; then they
    may run at the same time on different threads, and state that several
    spacecraft share needs a lock. `BaseThrusterModel` considers new maneuver
    requests on its control ticks, after guidance due at the same time. Startup
    and atmosphere-exit bookkeeping do not force an extra request-processing tick;
    a control tick pending at resume is still taken once. Accepted burn windows
    and explicitly installed plans retain their
    ignition and stop times, including boundaries between control ticks.

    Before integration starts, `run_simulation` checks that every control
    effector has a method of this function and, unless its actuation is
    `:mounted`, of `control_wrench`, and stops with an `ArgumentError` naming the
    effector otherwise.
    """
    function update_control! end

    using ..QuaternionMath
    include(joinpath(@__DIR__, "maneuvers", "propulsive_maneuvers.jl"))
    include(joinpath(@__DIR__, "aerobraking", "targeting_control.jl"))
    include(joinpath(@__DIR__, "aerobraking", "edg_heat_accounting.jl"))
    include(joinpath(@__DIR__, "aerobraking", "energy_depletion_lifecycle.jl"))
    include(joinpath(@__DIR__, "rpo", "lqmpc.jl"))
    include(joinpath(@__DIR__, "rpo", "rpo_control_types.jl"))
    include(joinpath(@__DIR__, "..", "..", "vehicle", "actuators", "thruster_allocation.jl"))
    include(joinpath(@__DIR__, "rpo", "reaction_wheel_allocator.jl"))
    include(joinpath(@__DIR__, "rpo", "rpo_mpc_control_model.jl"))
    include(joinpath(@__DIR__, "robot_arm", "robot_arm_control.jl"))
    include(joinpath(@__DIR__, "attitude", "momentum_manager.jl"))
    include(joinpath(@__DIR__, "landing", "apollo_descent_control.jl"))
    # KS aerobraking MPC: data contracts, prediction model, QP, and runtime hooks.
    include(joinpath(@__DIR__, "aerobraking_mpc", "types.jl"))
    include(joinpath(@__DIR__, "aerobraking_mpc", "constraints.jl"))
    include(joinpath(@__DIR__, "aerobraking_mpc", "condensed_core.jl"))
    include(joinpath(@__DIR__, "aerobraking_mpc", "objectives.jl"))
    include(joinpath(@__DIR__, "aerobraking_mpc", "density_adapter.jl"))
    include(joinpath(@__DIR__, "aerobraking_mpc", "scenario_config.jl"))
    include(joinpath(@__DIR__, "aerobraking_mpc", "trajectory_reference.jl"))
    include(joinpath(@__DIR__, "aerobraking_mpc", "prediction_model.jl"))
    include(joinpath(@__DIR__, "aerobraking_mpc", "qp_solver.jl"))
    include(joinpath(@__DIR__, "aerobraking_mpc", "area_mapping.jl"))
    include(joinpath(@__DIR__, "aerobraking_mpc", "plan_validation.jl"))
    include(joinpath(@__DIR__, "aerobraking_mpc", "spaceagora_control_model.jl"))
    include(joinpath(@__DIR__, "aerobraking_mpc", "campaign_control_model.jl"))

    # Per-spacecraft GNC binding (see SimulationLifecycle.bind_spacecraft).
    function SimulationLifecycle.bind_spacecraft(m::MagneticMomentumManagerModel, sat_idx::Int)
        m.spacecraft_index == sat_idx && return m
        b = SimulationLifecycle._shallow_copy(m)
        b.spacecraft_index = sat_idx
        b.mounted_levels = Float64[]          # its own buffer, not the template's
        return b
    end
    function SimulationLifecycle.bind_spacecraft(m::RobotArmControlEffector, sat_idx::Int)
        m.spacecraft_index == sat_idx && return m
        b = SimulationLifecycle._shallow_copy(m)
        b.spacecraft_index = sat_idx
        b.controller = deepcopy(m.controller)  # its own MPC warm start, not the template's
        b.held = deepcopy(m.held)              # its own held command, starting from the template's
        b.run_start, b.run_written = deepcopy(m.run_start), deepcopy(m.run_written)
        return b
    end
    function SimulationLifecycle.bind_spacecraft(m::ApolloDescentControlModel, sat_idx::Int)
        ids = _descent_indices(m.state, (sat_idx,))
        m.spacecraft_indices == ids && return m
        return ApolloDescentControlModel(m.config, m.guidance, m.state, m.terrain, m.actuators, ids, m.actuation)
    end

    # The spacecraft each built-in control model acts on, and the command its
    # output hooks read between ticks (GNCRuntime). Each held command is returned
    # as a value or a named tuple of values, not the model's own storage. The
    # descent controller records its tick times; the others do not. Every built-in
    # control model stops ticking for an inactive spacecraft (`skip_inactive`, true by
    # default): its tick writes only that spacecraft's command and run state, which
    # the right-hand side no longer reads. The descent controller declares it,
    # since its `release_command!` turns the engine off.
    GNCRuntime.spacecraft_scope(m::RPOMPCControlModel) = CoupledSpacecraft(m.chaser_index, (m.target_index,))
    GNCRuntime.held_command(m::RPOMPCControlModel, spacecraft_index) =
        spacecraft_index == m.chaser_index ?
            (force_ii=m.held.force_ii, torque_body=m.held.torque_body, thruster_forces_n=m.held.thruster_forces_n,
             rw_torque_body=m.held.rw_torque_body) : nothing

    GNCRuntime.spacecraft_scope(m::RobotArmControlEffector) = BoundSpacecraft((m.spacecraft_index,))
    # A controller with a plan buffer reads the plans its planner installs there.
    GNCRuntime.gnc_links(m::RobotArmControlEffector) =
        m.plan_buffer === nothing ? () : (GNCRuntime.Consumes(m.plan_buffer),)
    GNCRuntime.held_command(m::RobotArmControlEffector, spacecraft_index) =
        spacecraft_index == m.spacecraft_index ?
            (joint_torque_nm=copy(m.held.joint_torque_nm), base_force_ii=m.held.base_force_ii,
             base_torque_body=m.held.base_torque_body) : nothing

    GNCRuntime.spacecraft_scope(m::MagneticMomentumManagerModel) = BoundSpacecraft((m.spacecraft_index,))
    GNCRuntime.held_command(m::MagneticMomentumManagerModel, spacecraft_index) =
        spacecraft_index == m.spacecraft_index ? (torque_body=m.held_torque_body, dipole_am2=m.held_dipole_am2) : nothing
    GNCRuntime.held_command_time(m::MagneticMomentumManagerModel, spacecraft_index) =
        spacecraft_index == m.spacecraft_index ? m.last_update_s : NaN

    GNCRuntime.spacecraft_scope(m::ApolloDescentControlModel) = BoundSpacecraft(m.spacecraft_indices)
    GNCRuntime.skip_inactive(::ApolloDescentControlModel) = true
    GNCRuntime.held_command(m::ApolloDescentControlModel, spacecraft_index) =
        spacecraft_index in m.spacecraft_indices ?
            (thrust_n=m.actuators.thrust_n[spacecraft_index], torque_nm=m.actuators.torque_nm[spacecraft_index]) :
            nothing
    GNCRuntime.held_command_time(m::ApolloDescentControlModel, spacecraft_index) =
        spacecraft_index in m.spacecraft_indices ? m.actuators.last_update_s[spacecraft_index] : NaN

    # The energy-depletion controller holds the panels' commanded angle of attack
    # (rad), NaN before its first command.
    GNCRuntime.held_command(m::AerobrakingEnergyDepletionControlModel, spacecraft_index) =
        1 <= spacecraft_index <= length(m.state.last_alpha_rad) ? m.state.last_alpha_rad[spacecraft_index] : nothing

    # The burn scheduler holds a schedule: the effective window, its thrust and
    # whether the burn is on (NaN window ends when none is planned).
    function GNCRuntime.held_command(m::BaseThrusterModel, spacecraft_index)
        1 <= spacecraft_index <= length(m.burning) || return nothing
        start_s, stop_s, thrust_n = m.burn_windows[spacecraft_index]
        return (window_start_s=start_s, window_stop_s=stop_s, thrust_n=thrust_n, burning=m.burning[spacecraft_index])
    end
    # Its release, at each tick while the spacecraft is inactive: in `:mounted`
    # actuation the engine and the wheels off, as between burns. The schedule stays.
    function GNCRuntime.release_command!(m::BaseThrusterModel, t, spacecraft_index)
        m.actuation === :mounted && 1 <= spacecraft_index <= length(m.mounted_levels) &&
            fill!(m.mounted_levels[spacecraft_index], 0.0)
        return nothing
    end
end
