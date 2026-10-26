module RPOComparisonPlannerObjectTests
using Test, SpaceAGORA, Random, LinearAlgebra, StaticArrays, Logging
const P = SpaceAGORA.RPOPlannerInterfaces
const GA = SpaceAGORA.GNCAlgorithms
const S = SpaceAGORA.SimulationModel
const G = S.GuidanceHooks

const IDS = (:hypr, :pso_unrefined, :rrt_connect, :rrt_connect_bezier, :rrt_star, :chomp, :stomp)
const TYPES = (HYPRComparisonPlanner, UnrefinedPSOPlanner, RRTConnectPlanner, RRTConnectBezierPlanner,
    RRTStarPlanner, CHOMPPlanner, STOMPPlanner)
const ALIASES = (:pso => :hypr, :unrefined_hypr => :pso_unrefined, :rrtconnect => :rrt_connect,
    "RRT-Connect-Bezier" => :rrt_connect_bezier, :rrtstar => :rrt_star)

geometry(; keepout=0.6) = S.RPOReferenceGeometry(S.RPOStationGeometry([0.0 0.0; 0.0 2.0; 0.0 0.0]; keepout_radius_m=keepout);
    chaser=S.RPOCubeSatGeometry(dims_m=(0.1, 0.1, 0.1)))
const pso_cfg = G.rpo_740_mpc_final_pso_config(safe_distance_m=0.2; n_particles=8, n_iters=3, n_waypoints=2,
    adaptive_enable=false, reexplore_enable=false, cull_enable=false, sample_ds_m=0.25, retime_dt_s=0.5,
    retime_a_max_mps2=0.05, retime_max_steps=5_000, probe_candidates=4, refinement_rounds=2,
    refinement_waypoint_passes=2)
function comparison_config(; kwargs...)
    defaults = (pso_config=pso_cfg,
        rrt_connect=G.RPORRTConnectSettings(n_iters=120, step_size_m=0.5, shortcut_iters=6),
        rrt_star=G.RPORRTStarSettings(n_iters=90, step_size_m=0.5, shortcut_iters=6),
        chomp=G.RPOCHOMPSettings(n_iters=3), stomp=G.RPOSTOMPSettings(n_iters=3, n_rollouts=4),
        optimizer=G.RPOTrajectoryOptimizerSettings(match_hypr_runtime=false, runtime_limit_s=Inf, force_full_iters=true),
        tracking=G.RPOLQMPCTrackingSettings(dt_s=0.5, horizon=4, propellant_mass_kg=0.5, settle_time_s=2.0,
            final_position_tol_m=0.75, u_max_mps2=SVector(0.05, 0.05, 0.05)),
        safe_distance_m=0.2, write_plotly_outputs=false, write_failed_path_outputs=false, show_progress=false, rng_seed=741)
    return G.RPOPlannerComparisonConfig(; merge(defaults, (; kwargs...))...)
end
# The object each canonical identifier selects, written with the public constructors.
function object_for(id, cfg)
    id === :hypr && return HYPRComparisonPlanner(cfg.pso_config)
    id === :pso_unrefined && return UnrefinedPSOPlanner(cfg.pso_config)
    id === :rrt_connect && return RRTConnectPlanner(cfg.rrt_connect; config=cfg.pso_config, optimizer=cfg.optimizer)
    id === :rrt_connect_bezier && return RRTConnectBezierPlanner(cfg.rrt_connect; config=cfg.pso_config, optimizer=cfg.optimizer)
    id === :rrt_star && return RRTStarPlanner(cfg.rrt_star; config=cfg.pso_config, optimizer=cfg.optimizer)
    id === :chomp && return CHOMPPlanner(cfg.chomp; config=cfg.pso_config, optimizer=cfg.optimizer, seed=cfg.rrt_connect)
    id === :stomp && return STOMPPlanner(cfg.stomp; config=cfg.pso_config, optimizer=cfg.optimizer, seed=cfg.rrt_connect)
end
same_fields(a, b) = typeof(a) == typeof(b) && all(f -> isequal(getfield(a, f), getfield(b, f)), fieldnames(typeof(a)))

# Bit-for-bit equality of planner outputs, skipping measured wall-clock fields.
const TIMING = (:planner_compute_time, :planner_plan_time, :planner_refinement_time)
bits(x::Float64) = reinterpret(UInt64, x)
same_bits(a::Float64, b::Float64) = bits(a) == bits(b)
same_bits(a::AbstractArray{<:AbstractFloat}, b::AbstractArray{<:AbstractFloat}) =
    size(a) == size(b) && all(i -> same_bits(Float64(a[i]), Float64(b[i])), eachindex(a, b))
same_bits(a::Union{Tuple, AbstractVector}, b::Union{Tuple, AbstractVector}) =
    length(a) == length(b) && all(same_bits(x, y) for (x, y) in zip(a, b))
same_bits(a::NamedTuple, b::NamedTuple) = keys(a) == keys(b) &&
    all(k -> k in TIMING || same_bits(getfield(a, k), getfield(b, k)), keys(a))
same_bits(a::AbstractDict, b::AbstractDict) = keys(a) == keys(b) && all(k -> same_bits(a[k], b[k]), keys(a))
function same_bits(a, b)
    typeof(a) == typeof(b) || return false
    isstructtype(typeof(a)) && !isempty(fieldnames(typeof(a))) || return isequal(a, b)
    return all(f -> same_bits(getfield(a, f), getfield(b, f)), fieldnames(typeof(a)))
end
quietly(f) = with_logger(f, NullLogger())

const CASES = [G.RPOPlannerComparisonCase(start_rtn=SVector(-3.0, 0.0, 0.0), goal_rtn=SVector(3.0, 0.0, 0.0), label="blocked"),
    G.RPOPlannerComparisonCase(start_rtn=SVector(-3.0, 1.5, 0.5), goal_rtn=SVector(3.0, -1.0, 0.0), label="offset")]

@testset "Comparison planner objects: construction, identity and declarations" begin
    cfg = comparison_config()
    @test Tuple(typeof.(cfg.planners)) == TYPES
    @test Tuple(G.normalize_rpo_comparison_planner_type.(cfg.planners)) == IDS
    @test [G.rpo_comparison_planner_label(p) for p in cfg.planners] == [G.rpo_comparison_planner_label(id) for id in IDS]
    for (id, planner) in zip(IDS, cfg.planners)
        @test same_fields(planner, object_for(id, cfg))
        @test planner isa AbstractRPOPlanner
        @test GA.algorithm_family(planner) isa P.RPOReferenceFamily
        caps = planner_capabilities(planner)
        @test caps.state_sources == (:truth,) && caps.frames == (:target_rtn,)
        @test !caps.retiming && !caps.restart
        @test consumes_rng(planner) # every mode searches, swarms or seeds from the stream
    end
    # Pinning the comparison's safety distance keeps the selected objects, not names.
    pinned = G._rpo_comparison_config_with_fixed_safe_distance(cfg)
    @test all(pinned.planners[i] === cfg.planners[i] for i in eachindex(cfg.planners))
    @test pinned.pso_config.safe_distance_m == 0.5
    # Defaults follow the comparison convention; settings are owned copies.
    @test HYPRComparisonPlanner().config == G.rpo_740_mpc_final_pso_config(safe_distance_m=0.5)
    @test CHOMPPlanner(config=pso_cfg).settings == G.RPOCHOMPSettings(n_iters=pso_cfg.n_iters)
    @test STOMPPlanner(config=pso_cfg).settings == G.RPOSTOMPSettings(n_iters=pso_cfg.n_iters)
    @test RRTStarPlanner().max_reference_samples == 100_000
    @test_throws ArgumentError RRTConnectPlanner(; max_reference_samples=1)
    @test_throws ArgumentError HYPRComparisonPlanner(pso_cfg; max_reference_samples=typemax(Int))
end

@testset "Name selection was removed" begin
    cfg = comparison_config()
    # A Symbol or string is not a planner: the configuration and a direct plan refuse
    # it before HYPR is required or any planning starts.
    for selection in (:rrt_star, "rrt-star", :stopm)
        @test_throws MethodError comparison_config(planners=[selection])
        @test_throws MethodError G.rpo_plan_comparison_path(selection, CASES[1].start_rtn, CASES[1].goal_rtn,
            geometry(), cfg; rng=MersenneTwister(1))
    end
    @test_throws MethodError G.RPOPlannerComparisonConfig(planners=[:rrt_star, RRTStarPlanner()])
    # The identifiers and labels the reports use still accept names and objects, quietly.
    quiet, _ = Test.collect_test_logs() do
        for (name, id) in (map(id -> id => id, IDS)..., ALIASES...)
            @test G.normalize_rpo_comparison_planner_type(name) === id
        end
        foreach(G.rpo_comparison_planner_label, (IDS..., cfg.planners...))
        foreach(G.normalize_rpo_comparison_planner_type, cfg.planners)
    end
    @test isempty(quiet)
end

@testset "Object selections in batches" begin
    cfg = comparison_config()
    finite = comparison_config(optimizer=G.RPOTrajectoryOptimizerSettings(match_hypr_runtime=false,
        runtime_limit_s=1.0e6, runtime_max_iters=40, force_full_iters=true))
    # Batches: the default objects against the same objects listed, and repeats.
    batch(c) = quietly(() -> G.rpo_run_planner_comparison_batch(CASES, geometry(), c))
    # Every result and plan field but the measured planner times, the LQ-MPC closed
    # loop included.
    science(b) = (b.planner_types, b.results_by_planner, b.plans_by_planner)
    finite_kwargs = (optimizer=finite.optimizer,)
    defaults = batch(comparison_config(; finite_kwargs...))
    listed = batch(comparison_config(; planners=[object_for(id, finite) for id in IDS], finite_kwargs...))
    @test defaults.planner_types == collect(IDS)
    @test same_bits(science(defaults), science(listed))
    mixed = [:rrt_star, :hypr, :rrt_connect, :chomp, :rrt_star, :stomp]
    repeated = batch(comparison_config(; planners=[object_for(id, finite) for id in mixed], finite_kwargs...))
    @test repeated.planner_types == mixed
    # Runtime matching still moves the first HYPR selection to the front.
    matched = comparison_config(planners=[object_for(:rrt_connect, cfg), object_for(:hypr, cfg), object_for(:stomp, cfg),
            object_for(:hypr, cfg)],
        optimizer=G.RPOTrajectoryOptimizerSettings(match_hypr_runtime=true, runtime_max_iters=40))
    order = quietly(() -> G.rpo_run_planner_comparison_batch(CASES[1:1], geometry(), matched)).planner_types
    @test order == [:hypr, :rrt_connect, :stomp]
end

# Contract path fixtures, as in the HYPR adapter tests.
function request(; kwargs...)
    defaults = (request_id=1, chaser_id=101, target_id=201, epoch="fixture", time_s=0.0,
        x_rtn=(-3.0, 0.0, 0.0, 0.0, 0.0, 0.0), target_state_ii=(7e6, 0.0, 0.0, 0.0, 7500.0, 0.0),
        goal_rtn_m=(3.0, 0.0, 0.0), geometry=geometry(), geometry_revision="fixture-v1",
        constraints=P.RPOPlanningConstraints(clearance_m=0.15, max_speed_mps=0.5),
        reference_dt_s=0.5, preview_horizon_steps=2, valid_until_s=600.0)
    return P.RPOPlanningRequest(; merge(defaults, (; kwargs...))...)
end
validate(req, out) = P.validate_rpo_result(req, out; clearance_at=(p, g) -> S.rpo_clearance_distance_to_station(p, g))
# The unchanged kernels, called directly with the request mapping the contract documents.
function composed(planner, req, rng)
    cfg = planner.config; clearance = req.constraints.clearance_m
    budget = P.rpo_planning_budget(req, planner.headroom;
        position_scale_m=max(norm(collect(req.x_rtn[1:3])), norm(collect(req.goal_rtn_m))),
        velocity_scale_mps=max(something(req.constraints.max_speed_mps, cfg.retime_max_speed_mps), cfg.retime_initial_speed_mps))
    input = rpo_pso_config(cfg; safe_distance_m=clearance, retime_dt_s=req.reference_dt_s,
        retime_max_speed_mps=min(cfg.retime_max_speed_mps, budget.max_speed_mps), retime_a_max_mps2=cfg.retime_a_max_mps2,
        retime_max_steps=min(cfg.retime_max_steps, planner.max_reference_samples - 1))
    a, b = SVector{3, Float64}(req.x_rtn[1:3]), SVector{3, Float64}(req.goal_rtn_m)
    seed() = G.rpo_rrt_connect_plan_path(a, b, req.geometry, input; safe_distance_m=clearance,
        settings=planner.seed, rng=rng).path
    raw = if planner isa HYPRComparisonPlanner
        S.rpo_pso_plan_path(a, b, req.geometry, input; safe_distance_m=clearance, rng=rng)
    elseif planner isa UnrefinedPSOPlanner
        S.rpo_pso_plan_path(a, b, req.geometry, rpo_pso_config(input; refinement_enable=false, adaptive_enable=false,
            schedule_enable=false, stagnation_learning_enable=false, rrt_warmstart_enable=false, cull_enable=false,
            reexplore_enable=false); safe_distance_m=clearance, rng=rng)
    elseif planner isa RRTConnectPlanner
        G.rpo_rrt_connect_plan_path(a, b, req.geometry, input; safe_distance_m=clearance, settings=planner.settings, rng=rng)
    elseif planner isa RRTConnectBezierPlanner
        G.rpo_rrt_connect_bezier_plan_path(a, b, req.geometry, input; safe_distance_m=clearance, settings=planner.settings, rng=rng)
    elseif planner isa RRTStarPlanner
        G.rpo_rrt_star_plan_path(a, b, req.geometry, input; safe_distance_m=clearance, settings=planner.settings, rng=rng)
    elseif planner isa CHOMPPlanner
        G.rpo_chomp_plan_path(a, b, req.geometry, input; safe_distance_m=clearance, settings=planner.settings,
            optimizer=planner.optimizer, initial_path=seed())
    else
        G.rpo_stomp_plan_path(a, b, req.geometry, input; safe_distance_m=clearance, settings=planner.settings,
            optimizer=planner.optimizer, rng=rng, initial_path=seed())
    end
    return input, raw, G.rpo_reference_from_path(raw.path, req.geometry, raw.config; safe_distance_m=clearance)
end

@testset "plan_rpo! uses the request's limits and the unchanged kernels" begin
    # More RRT* iterations than the comparison fixture: with 90, this request's search finds no path.
    cfg = comparison_config(rrt_star=G.RPORRTStarSettings(n_iters=400, step_size_m=0.5, shortcut_iters=6))
    @test plan_rpo!(nothing, comparison_config().planners[5], request(), MersenneTwister(742)).termination ===
        :path_not_found
    for (id, planner) in zip(IDS, cfg.planners)
        req = request(); kept = deepcopy(req); settings = deepcopy(planner)
        rng = MersenneTwister(742); control = MersenneTwister(742)
        out = plan_rpo!(initialize_planner(planner, nothing), planner, req, rng)
        @test out.status === :candidate
        @test validate(req, out).accepted
        input, raw, (t, r, v) = composed(planner, req, control)
        @test out.reference.t_ref_s == t && out.reference.r_ref_rtn_m == r && out.reference.v_ref_rtn_mps == v
        @test rand(rng, UInt64, 2) == rand(control, UInt64, 2)
        d = out.diagnostics
        @test d.planner === id && d.input_config == input && d.effective_config == raw.config
        @test d.input_config.safe_distance_m == 0.15 # the request's clearance, not the comparison's 0.5 m
        @test d.requested_config == planner.config && d.path == raw.path && d.cost == raw.cost
        @test out.reference.r_ref_rtn_m[:, 1] == [-3.0, 0.0, 0.0] && out.reference.r_ref_rtn_m[:, end] == [3.0, 0.0, 0.0]
        @test out.reference.chaser_id == 101 && out.reference.target_id == 201
        @test out.reference.geometry_revision == "fixture-v1" && out.reference.frame === :target_rtn
        @test same_fields(planner, settings) && isequal(repr(req), repr(kept))
        # Repeating with a copied stream repeats the candidate.
        again = plan_rpo!(nothing, planner, deepcopy(req), MersenneTwister(742))
        @test again.reference.r_ref_rtn_m == out.reference.r_ref_rtn_m
    end
end

@testset "plan_rpo! refusals stay refusals" begin
    cfg = comparison_config()
    for planner in cfg.planners
        rng = MersenneTwister(5); before = copy(rng)
        plan(r; p=planner) = plan_rpo!(nothing, p, r, rng)
        @test plan(request(reason=:retime)).termination === :retiming_not_implemented
        @test plan(request(geometry=nothing)).termination === :unsupported_geometry
        accel = request(constraints=P.RPOPlanningConstraints(clearance_m=0.15, max_speed_mps=0.5, max_acceleration_mps2=0.1))
        @test plan(accel).termination === :acceleration_limited_retimer_required
        @test plan(request(x_rtn=(1e16, 0.0, 0.0, 0.0, 0.0, 0.0), goal_rtn_m=(nextfloat(1e16), 0.0, 0.0))).termination ===
            :insufficient_reference_precision
        @test rand(rng) == rand(before) # refused before any kernel draws
        # Unsupported retiming is refused, never planned again from scratch.
        @test P.validate_rpo_capabilities(planner, request(reason=:retime)).reason === :unsupported_retiming
        @test P.validate_rpo_capabilities(planner, request(); checkpoint=true).reason === :unsupported_restart
        active = P.RPOReference(t_ref_s=[0.0, 0.5], r_ref_rtn_m=zeros(3, 2), v_ref_rtn_mps=zeros(3, 2),
            origin_time_s=0.0, valid_until_s=600.0, chaser_id=101, target_id=201, geometry_revision="fixture-v1")
        @test retime_rpo!(nothing, planner, request(reason=:retime), active, rng).termination === :retiming_not_implemented
        @test plan(request(valid_until_s=1.0)).termination === :insufficient_reference_lifetime
        small = typeof(planner).name.wrapper === HYPRComparisonPlanner || planner isa UnrefinedPSOPlanner ?
            typeof(planner)(planner.config; max_reference_samples=3) :
            typeof(planner)(planner.settings; config=planner.config, optimizer=planner.optimizer,
                (planner isa Union{CHOMPPlanner, STOMPPlanner} ? (seed=planner.seed,) : ())..., max_reference_samples=3)
        @test plan(request(); p=small).termination === :reference_work_limit
    end
    # A search that finds no path fails; its direct-segment diagnostic is never a candidate.
    stuck = RRTConnectPlanner(G.RPORRTConnectSettings(n_iters=1, step_size_m=0.01); config=pso_cfg)
    lost = plan_rpo!(nothing, stuck, request(), MersenneTwister(1))
    @test lost.status === :failed && lost.termination === :path_not_found && lost.reference === nothing
    # With the acceleration-limited retimer, a curve meets an acceleration limit; a polyline's
    # corners do not, and the unchanged validator refuses it.
    accel_cfg = rpo_pso_config(pso_cfg; retime_accel_limit_enable=true, retime_a_max_mps2=0.02)
    accel = request(constraints=P.RPOPlanningConstraints(clearance_m=0.15, max_speed_mps=0.5, max_acceleration_mps2=0.025),
        valid_until_s=3000.0)
    search = G.RPORRTConnectSettings(n_iters=120, step_size_m=0.5, shortcut_iters=6)
    curve = plan_rpo!(nothing, RRTConnectBezierPlanner(search; config=accel_cfg), accel, MersenneTwister(742))
    @test curve.status === :candidate && validate(accel, curve).accepted
    @test curve.diagnostics.input_config.retime_a_max_mps2 == 0.02 && curve.diagnostics.path_kind === :bezier_control_polygon
    corners = plan_rpo!(nothing, RRTConnectPlanner(search; config=accel_cfg), accel, MersenneTwister(742))
    @test corners.status === :failed && corners.termination === :reference_rejected
    @test corners.diagnostics.validation.reason === :acceleration_limit
end

@testset "A comparison planner runs through the RPO planner binding" begin
    cfg = rpo_pso_config(pso_cfg; retime_accel_limit_enable=true, retime_max_speed_mps=0.5,
        retime_a_max_mps2=0.5 * 0.05 / 5.2, retime_dt_s=0.1)
    planner = RRTStarPlanner(G.RPORRTStarSettings(n_iters=60, step_size_m=0.5, shortcut_iters=4); config=cfg)
    args = make_rpo_configuration(planner=planner, mission_time_s=0.4)
    sol = run_simulation(args; return_solution=true)
    report = rpo_run_report(sol).planners[1]
    installed = filter(x -> x.event === :installed, report.records)
    @test length(installed) == 1 && only(installed).diagnostics.planner === :rrt_star
    @test report.planner == string(RRTStarPlanner)
    # A candidate the planner refuses is never installed: the run stops with the refusal.
    tiny = RRTStarPlanner(planner.settings; config=cfg, max_reference_samples=3)
    args = make_rpo_configuration(planner=tiny, mission_time_s=0.4)
    err = try run_simulation(args; return_solution=true); nothing catch e; e end
    @test err isa RPOPlanningError && err.reason === :failed
    @test isempty(filter(x -> x.event === :installed, err.report.records))
    @test !only(args.guidance_model.guidance_effectors).plan_buffer.valid
end
end
