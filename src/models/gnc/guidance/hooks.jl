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
    using ..GuidanceModels: RPOPlan, RPOPlanBuffer, update_rpo_plan_buffer!
    using ..NavigationHooks: RPOReferenceGeometry, RPOStationGeometry
    using ..NavigationHooks: rpo_clearance_distance_to_station, rpo_clearance_to_station, rpo_path_clearance_stats
    using ..CommandTypes: PropulsiveManeuverCommand, AerobrakingControlCommand
    using ..EphemeridesModels: ephemerides_frame_lpi
    using ..GravityEffectors: aerobraking_gravity_force_ii
    using ..ReferenceSystems
    using ..LinearAlgebra
    using ..StaticArrays
    using ..Kinematics
    using ..HYPRUtils
    using ...RPOPlannerInterfaces: AbstractRPOPlanner, RPOPlanningHeadroom
    using Base.Threads: @threads, maxthreadid, threadid
    using Logging: Logging
    using Random
    using AstroTime
    using SPICE

    const config = Structure
    # Guidance does not look up ControlHooks: what it needs from control comes by
    # dispatch on methods ControlHooks adds (`_is_edg_control`).

    export update_guidance!
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

    Effector `k` runs every `guidance_rates[k]` seconds of simulated time for
    the spacecraft in its declared scope, starting one period after the run
    starts. Checkpoint segments and resume retain that run clock. When ticks
    coincide, navigation runs first, followed by control and guidance; energy
    depletion guidance with targeting runs before control instead. Maneuver
    planning runs after all guidance due at that time. Startup and atmosphere
    exits do not force extra guidance ticks; a tick pending at resume is taken
    once on the retained clock.

    Before integration starts, `run_simulation` checks that every guidance
    effector has a method of this function for these arguments and stops with
    an `ArgumentError` naming the effector otherwise.
    """
    function update_guidance! end

    using ..FrameTransforms
    # Shared calculations and compatibility types/bindings retain this module.
    include(joinpath(@__DIR__, "..", "planning", "rpo", "sampling_settings.jl"))
    include(joinpath(@__DIR__, "..", "companions", "hypr", "pso_parameters.jl"))
    include(joinpath(@__DIR__, "..", "companions", "hypr", "execution_contracts.jl"))
    include(joinpath(@__DIR__, "..", "planning", "rpo", "path_geometry.jl"))
    include(joinpath(@__DIR__, "..", "planning", "rpo", "profile_evaluation.jl"))
    include(joinpath(@__DIR__, "..", "planning", "rpo", "path_retiming.jl"))
    include(joinpath(@__DIR__, "..", "planning", "rpo", "path_sampling.jl"))
    include(joinpath(@__DIR__, "..", "planning", "rpo", "path_metrics.jl"))
    include(joinpath(@__DIR__, "..", "planning", "rpo", "rrt.jl"))
    include(joinpath(@__DIR__, "..", "planning", "rpo", "comparison", "trajectory_optimizers.jl"))
    include(joinpath(@__DIR__, "..", "planning", "rpo", "comparison", "planner_objects.jl"))
    include(joinpath(@__DIR__, "..", "planning", "rpo", "comparison", "planner_comparison.jl"))
    include(joinpath(@__DIR__, "..", "companions", "hypr", "replanning.jl"))
    include(joinpath(@__DIR__, "aerobraking", "typed_edg", "guidance_models.jl"))
    include(joinpath(@__DIR__, "aerobraking", "energy_depletion_lifecycle.jl"))

    include(joinpath(@__DIR__, "maneuvers", "thruster_guidance_functions.jl"))
    include(joinpath(@__DIR__, "landing", "apollo_descent_guidance.jl"))

    # Per-spacecraft GNC binding (see SimulationLifecycle.bind_spacecraft).
    function SimulationLifecycle.bind_spacecraft(m::ApolloDescentGuidanceModel, sat_idx::Int)
        ids = _descent_indices(m.state, (sat_idx,))
        m.spacecraft_indices == ids && return m
        return ApolloDescentGuidanceModel(m.config, m.state, m.terrain, ids)
    end
    # The maneuver guidance models, bound to one spacecraft. A periapsis-raise copy
    # gets its own command states, so copies bound from one template share no run
    # state; a campaign copy shares the template's maneuver table, which a run only
    # reads.
    function SimulationLifecycle.bind_spacecraft(m::ApoapsisTargetPeriapsisRaiseGuidanceModel, sat_idx::Int)
        m.spacecraft_index == sat_idx && return m
        b = SimulationLifecycle._shallow_copy(m)
        b.spacecraft_index = sat_idx
        b.command_state = copy(m.command_state)
        return b
    end
    function SimulationLifecycle.bind_spacecraft(m::AerobrakingCampaignPropulsiveManeuverGuidanceModel, sat_idx::Int)
        m.spacecraft_index == sat_idx && return m
        return AerobrakingCampaignPropulsiveManeuverGuidanceModel(m.maneuver_orbit_number, m.maneuver_Δv,
            m.maneuver_flight_apoapsis_radius_m, sat_idx)
    end

    # The spacecraft each built-in guidance model acts on (GNCRuntime): the descent
    # guidance flies its selected spacecraft and stops at touchdown. Like every
    # built-in guidance model, it does not plan for an inactive spacecraft
    # (`skip_inactive`, true by default; the descent guidance declares it).
    GNCRuntime.spacecraft_scope(m::ApolloDescentGuidanceModel) = GNCRuntime.BoundSpacecraft(m.spacecraft_indices)
    GNCRuntime.skip_inactive(::ApolloDescentGuidanceModel) = true
    # A maneuver guidance model acts on every spacecraft, or on the one it is bound to.
    GNCRuntime.spacecraft_scope(m::Union{ApoapsisTargetPeriapsisRaiseGuidanceModel,
            AerobrakingCampaignPropulsiveManeuverGuidanceModel}) =
        m.spacecraft_index == 0 ? GNCRuntime.AllSpacecraft() : GNCRuntime.BoundSpacecraft((m.spacecraft_index,))

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
    # The same run state goes into checkpoint files, and a resumed run puts it back
    # after the reset: the descent guidance's entries of the shared state (which its
    # control reads) and the periapsis-raise guidance's command states.
    SimulationLifecycle.checkpoint_state(m::ApolloDescentGuidanceModel) =
        (entries=_descent_entries(m.state, m.spacecraft_indices),)
    function SimulationLifecycle.restore_checkpoint_state!(m::ApolloDescentGuidanceModel, data)
        _restore_descent_entries!(m.state, m.spacecraft_indices, data.entries)
        return nothing
    end
    # A periapsis-raise copy bound to a spacecraft keeps only that spacecraft's entry,
    # `nothing` before its first tick; the model for every spacecraft keeps them all.
    function SimulationLifecycle.checkpoint_state(m::ApoapsisTargetPeriapsisRaiseGuidanceModel)
        k = m.spacecraft_index
        k == 0 && return (command_state=copy(m.command_state),)
        return (entry=k <= length(m.command_state) ? m.command_state[k] : nothing,)
    end
    function SimulationLifecycle.restore_checkpoint_state!(m::ApoapsisTargetPeriapsisRaiseGuidanceModel, data)
        k = m.spacecraft_index
        if k == 0
            copy!(m.command_state, data.command_state)
        elseif data.entry !== nothing
            _ensure_apo_target_state!(m, k)
            m.command_state[k] = data.entry
        end
        return nothing
    end
end
