module GuidanceHooks
    import ...SimulationLifecycle
    import ...GNCRuntime
    using ..Structure
    using ..Geodesy: geodetic_altitude, ellipsoid_surface_radius, radius_for_geodetic_altitude

    using ..ConfigTypes: ODEParams, Solution
    using ..AbstractTypes: AbstractTerrainModel
    using ..TerrainModels: NoTerrainModel, DEMTerrainModel, terrain_height, terrain_radius
    using ..QuaternionMath: dcm_to_quaternion
    using ..AbstractTypes: AbstractPlanet, AbstractControlEffectorModel, AbstractGuidanceModel
    using ..GuidanceModels: AerobrakingCampaignPropulsiveManeuverGuidanceModel
    using ..GuidanceModels: ApoapsisTargetPeriapsisRaiseGuidanceModel
    using ..GuidanceModels: RPOGuidanceModel, RPOPlan, RPOPlanBuffer, update_rpo_plan_buffer!
    using ..NavigationHooks: RPOReferenceGeometry, RPOStationGeometry
    using ..NavigationHooks: rpo_clearance_distance_to_station, rpo_clearance_to_station, rpo_path_clearance_stats
    using ..CommandTypes: PropulsiveManeuverCommand, AerobrakingControlCommand
    using ..EphemeridesModels: _ephemerides_frame_lpi
    using ..GravityEffectors: aerobraking_gravity_force_ii
    using ..ReferenceSystems
    using ..LinearAlgebra
    using ..StaticArrays
    using ..Kinematics
    using ..HYPRUtils
    using ..HookDeprecation: _hook_deprecation
    using ...RPOPlannerInterfaces: AbstractRPOPlanner, RPOPlanningHeadroom
    using Base.Threads: @threads, maxthreadid, threadid
    using Logging: Logging
    using Random
    using AstroTime
    using SPICE

    const config = Structure
    const _PARENT = parentmodule(@__MODULE__)

    @inline function _control_module()
        return getfield(_PARENT, :ControlHooks)
    end

    export update_guidance!, calcGuidanceEffect!
    export RPOPSOAdaptiveSettings, RPOPSOConfig, RPOPSOConfigurator, RPOPSOCullSettings
    export RPOPSOEarlyStoppingSettings, RPOPSOObjectiveSettings
    export RPOPSOProbeSettings, RPOPSOReexploreSettings, RPOPSORefinementSettings
    export RPOPSORetimingSettings, RPOPSOScheduleSettings, RPOPSOSwarmSettings
    export RPOPSORRTConnectWarmstartSettings, RPOPSOStagnationSettings, rpo_pso_config
    export rpo_pso_plan_path, rpo_reference_from_path
    export RPOReplanningConfig, RPOReplanningSphere
    export rpo_active_replanning_spheres, rpo_geometry_with_replanning_spheres
    export rpo_reference_tracking_error, rpo_remaining_reference_path
    export rpo_replanning_decision, rpo_replanning_sphere_center
    export RPOCHOMPSettings, RPORRTConnectSettings, RPORRTStarSettings, RPOSTOMPSettings, RPOTrajectoryOptimizerSettings
    export RPOLQMPCTrackingSettings, RPOPlannerComparisonCase, RPOPlannerComparisonConfig
    export normalize_rpo_comparison_planner_type, rpo_comparison_planner_label
    export HYPRComparisonPlanner, UnrefinedPSOPlanner, RRTConnectPlanner, RRTConnectBezierPlanner, RRTStarPlanner, CHOMPPlanner, STOMPPlanner
    export rpo_740_mpc_final_pso_config, rpo_chomp_plan_path, rpo_rrt_connect_plan_path, rpo_rrt_connect_bezier_plan_path, rpo_rrt_star_plan_path, rpo_stomp_plan_path
    export rpo_plan_comparison_path, rpo_run_planner_comparison_batch
    export rpo_track_retimed_path_lqmpc, rpo_flatten_planner_results
    export rpo_comparison_cost_iteration_plot, rpo_write_planner_comparison_outputs
    export AerobrakingControlCommand
    export AerobrakingEnergyDepletionConfig, AerobrakingEnergyDepletionState
    export AerobrakingEnergyDepletionGuidanceModel
    export DescentPhaseTargets, ApolloDescentConfig, ApolloDescentState, ApolloDescentGuidanceModel
    export apollo11_descent_targets, descent_site_frame, descent_accel_command, descent_time_to_go
    export descent_attitude_command, descent_throttle, descent_environment, quaternion_from_passive_dcm

    """
        update_guidance!(model, u, p, t, spacecraft_index) -> nothing

    The guidance tick, which every effector in a `GuidanceModel` implements for
    its own type. The effector updates its own state (a plan, a target, a mode)
    for control effectors to read; the engine ignores the return value.

    - `model`: the guidance effector.
    - `u`: the full integrated state. `u.sc[spacecraft_index]` is the spacecraft, with
      `pos` (m) and `vel` (m/s) relative to the planet's centre in the inertial
      frame and `mass` (kg); with attitude propagation also `q = (x, y, z, w)`,
      scalar last, with `v_body = rot(q) * v_inertial` (body means bus axes), and
      `ω`, the bus angular velocity in rad/s.
    - `p`: the run's parameters; `p.args` is the `SimulationConfiguration`.
    - `t`: elapsed time since the run's start epoch, in seconds (`Float64`).
    - `spacecraft_index`: the spacecraft's position in the run (`Int`, from 1), not its ID.

    Effector `k` runs every `guidance_rates[k]` seconds of simulated time, once
    for each spacecraft in index order, starting one period after the solve
    starts. When ticks coincide, navigation effectors run first, then control
    effectors, then guidance effectors; energy depletion guidance with targeting
    runs before the control effectors instead. The propulsive-maneuver thruster
    control effector (see `update_control!`) also runs every guidance
    effector at the start of the solve and when a spacecraft leaves the
    atmosphere, just before it schedules that spacecraft's burns.

    Before integration starts, `run_simulation` checks that every guidance
    effector has a method for these arguments (of this function or of the
    deprecated `calcGuidanceEffect!`) and stops with an `ArgumentError` naming
    the effector otherwise. This hook replaces the deprecated
    `calcGuidanceEffect!`.
    """
    function update_guidance! end

    """
        calcGuidanceEffect!(model, u, p, t, spacecraft_index) -> nothing

    Deprecated, use [`update_guidance!`](@ref); this name will be removed in
    SpaceAGORA 0.5.0. Only the name changes: the arguments, their order and when
    the engine calls the hook are those of `update_guidance!`.

    Until then the engine calls an effector's own method of this hook when its
    type has no `update_guidance!`, and `run_simulation` warns once per type
    before integration starts. A direct call on an effector that implements
    `update_guidance!` runs it, and warns once per type that the name is
    deprecated (under `julia --depwarn=error` it throws instead).
    """
    function calcGuidanceEffect! end

    using ..FrameTransforms
    # Shared calculations and compatibility types/bindings retain this module.
    include(joinpath(@__DIR__, "..", "planning", "rpo", "sampling_settings.jl"))
    include(joinpath(@__DIR__, "..", "planning", "adapters", "hypr", "pso_parameters.jl"))
    include(joinpath(@__DIR__, "..", "planning", "adapters", "hypr", "execution_contracts.jl"))
    include(joinpath(@__DIR__, "..", "planning", "rpo", "path_geometry.jl"))
    include(joinpath(@__DIR__, "..", "planning", "rpo", "profile_evaluation.jl"))
    include(joinpath(@__DIR__, "..", "planning", "rpo", "path_retiming.jl"))
    include(joinpath(@__DIR__, "..", "planning", "rpo", "path_sampling.jl"))
    include(joinpath(@__DIR__, "..", "planning", "rpo", "path_metrics.jl"))
    include(joinpath(@__DIR__, "..", "planning", "rpo", "rrt.jl"))
    include(joinpath(@__DIR__, "..", "planning", "rpo", "comparison", "trajectory_optimizers.jl"))
    include(joinpath(@__DIR__, "..", "planning", "rpo", "comparison", "planner_objects.jl"))
    include(joinpath(@__DIR__, "..", "planning", "rpo", "comparison", "planner_comparison.jl"))
    include(joinpath(@__DIR__, "..", "planning", "adapters", "hypr", "replanning.jl"))
    include(joinpath(@__DIR__, "rpo", "rpo_guidance_hooks.jl"))
    include(joinpath(@__DIR__, "aerobraking", "typed_edg", "guidance_models.jl"))
    include(joinpath(@__DIR__, "aerobraking", "energy_depletion_lifecycle.jl"))

    include(joinpath(@__DIR__, "thruster_guidance", "thruster_guidance_functions.jl"))
    include(joinpath(@__DIR__, "landing", "apollo_descent_guidance.jl"))

    # Per-spacecraft GNC binding (see SimulationLifecycle.bind_spacecraft).
    function SimulationLifecycle.bind_spacecraft(m::ApolloDescentGuidanceModel, sat_idx::Int)
        ids = _descent_indices(m.state, (sat_idx,))
        m.spacecraft_indices == ids && return m
        return ApolloDescentGuidanceModel(m.config, m.state, m.terrain, ids)
    end

    # The spacecraft each built-in guidance model acts on (GNCRuntime). The legacy
    # RPO guidance plans for the chaser from the target's state; the descent
    # guidance flies its selected spacecraft and stops at touchdown.
    GNCRuntime.spacecraft_scope(m::RPOGuidanceModel) = GNCRuntime.CoupledSpacecraft(m.chaser_idx, (m.target_idx,))
    GNCRuntime.spacecraft_scope(m::ApolloDescentGuidanceModel) = GNCRuntime.BoundSpacecraft(m.spacecraft_indices)
    GNCRuntime.skip_inactive(::ApolloDescentGuidanceModel) = true

    # Run state resets at each run's start (G-8): the descent guidance's entries of
    # the shared state, and the periapsis-raise guidance's command states.
    function SimulationLifecycle.initialize_guidance!(m::ApolloDescentGuidanceModel, u, p, t)
        _reset_descent_entries!(m.state, m.spacecraft_indices)
        return nothing
    end
    function SimulationLifecycle.initialize_guidance!(m::ApoapsisTargetPeriapsisRaiseGuidanceModel, u, p, t)
        empty!(m.command_state)
        return nothing
    end

    # The deprecated RPO guidance rewrites its goal, its plan buffer (shared with the
    # controller, so restored in place) and its replanning record during a run; each
    # run starts them from their values before the model's first run.
    const _RPO_GUIDANCE_RUN_FIELDS = (:goal_rtn, :force_replan, :replanning_events, :replan_count, :retime_count,
        :safe_hold_count, :replan_failure_count, :last_replanning_time_s, :last_replanning_signature,
        :replanning_persistence_count, :replanning_phase, :pending_replan)
    function SimulationLifecycle.initialize_guidance!(m::RPOGuidanceModel, u, p, t)
        buffer = m.plan_buffer
        if m.run_start === nothing
            m.run_start = deepcopy((fields=NamedTuple{_RPO_GUIDANCE_RUN_FIELDS}(map(f -> getfield(m, f), _RPO_GUIDANCE_RUN_FIELDS)),
                valid=buffer.valid, plan=buffer.plan, updated_at_s=buffer.updated_at_s))
        else
            start = deepcopy(m.run_start)
            foreach(f -> setfield!(m, f, getfield(start.fields, f)), _RPO_GUIDANCE_RUN_FIELDS)
            buffer.valid, buffer.plan, buffer.updated_at_s = start.valid, start.plan, start.updated_at_s
        end
        return nothing
    end
end
