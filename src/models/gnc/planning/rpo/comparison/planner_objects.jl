# Planner objects for the seven RPO comparison modes. Each object holds a copy of
# its settings. The comparison harness (planner_comparison.jl) runs it under the
# comparison's fixed safety distance; the contract adapter (planner_contract.jl)
# runs the same kernels under a request's own clearance and limits.

using ..SimConfig: _with_fields

"""Internal supertype of the seven RPO comparison planners."""
abstract type AbstractRPOComparisonPlanner <: AbstractRPOPlanner end

function _rpo_comparison_reference_samples(max_reference_samples::Integer)
    2 <= max_reference_samples < typemax(Int) || throw(ArgumentError("Invalid reference sample budget."))
    return Int(max_reference_samples)
end
_rpo_comparison_owned_config(config::RPOPSOConfig) = deepcopy(rpo_pso_config(config))
_rpo_comparison_default_config() =
    rpo_740_mpc_final_pso_config(safe_distance_m=RPO_PLANNER_COMPARISON_SAFE_DISTANCE_M)

"""
    HYPRComparisonPlanner(config=<comparison HYPR configuration>; headroom, max_reference_samples)

The `:hypr` mode of the RPO planner comparison: HYPR's swarm optimizer with every
feature `config` enables. It holds a copy of `config`, an `RPOPSOConfig`, which
also sets the cost and the retiming of the planned path. The default is the
comparison's tuned configuration with its 0.5 m safety distance.

As an `AbstractRPOPlanner`, `plan_rpo!` plans under the request's clearance,
reference step and limits; `headroom` and `max_reference_samples` have the meaning
they have for `HYPRRPOPlanner`. This is a comparison mode, not a replacement for
`HYPRRPOPlanner`: it does not retime and has no replan warm-start option. All seven
comparison planners need the optional HYPR companion and draw from the planner's
random stream.
"""
struct HYPRComparisonPlanner <: AbstractRPOComparisonPlanner
    config::RPOPSOConfig
    headroom::RPOPlanningHeadroom
    max_reference_samples::Int
    function HYPRComparisonPlanner(config::RPOPSOConfig=_rpo_comparison_default_config();
            headroom::RPOPlanningHeadroom=RPOPlanningHeadroom(), max_reference_samples::Integer=100_000)
        new(_rpo_comparison_owned_config(config), headroom, _rpo_comparison_reference_samples(max_reference_samples))
    end
end

"""
    UnrefinedPSOPlanner(config=<comparison HYPR configuration>; headroom, max_reference_samples)

The `:pso_unrefined` comparison mode: the HYPR swarm kernel of
`HYPRComparisonPlanner` with refinement, adaptation, scheduling, stagnation
learning, the RRT warm start, culling and re-exploration switched off when it
plans. The other settings come from its copy of `config`.
"""
struct UnrefinedPSOPlanner <: AbstractRPOComparisonPlanner
    config::RPOPSOConfig
    headroom::RPOPlanningHeadroom
    max_reference_samples::Int
    function UnrefinedPSOPlanner(config::RPOPSOConfig=_rpo_comparison_default_config();
            headroom::RPOPlanningHeadroom=RPOPlanningHeadroom(), max_reference_samples::Integer=100_000)
        new(_rpo_comparison_owned_config(config), headroom, _rpo_comparison_reference_samples(max_reference_samples))
    end
end

"""
    RRTConnectPlanner(settings=RPORRTConnectSettings(); config, optimizer, headroom, max_reference_samples)

The `:rrt_connect` comparison mode: bidirectional RRT-Connect with the search
`settings`, scored, refined and retimed with `config` (an `RPOPSOConfig`, default
the comparison's HYPR configuration). `optimizer`, an
`RPOTrajectoryOptimizerSettings`, supplies the comparison's iteration policy:
`match_hypr_iters` takes the iteration count from `config`, and `runtime_max_iters`
replaces it when the comparison sets a finite runtime limit.
"""
struct RRTConnectPlanner <: AbstractRPOComparisonPlanner
    config::RPOPSOConfig
    settings::RPORRTConnectSettings
    optimizer::RPOTrajectoryOptimizerSettings
    headroom::RPOPlanningHeadroom
    max_reference_samples::Int
    function RRTConnectPlanner(settings::RPORRTConnectSettings=RPORRTConnectSettings();
            config::RPOPSOConfig=_rpo_comparison_default_config(),
            optimizer::RPOTrajectoryOptimizerSettings=RPOTrajectoryOptimizerSettings(),
            headroom::RPOPlanningHeadroom=RPOPlanningHeadroom(), max_reference_samples::Integer=100_000)
        new(_rpo_comparison_owned_config(config), deepcopy(settings), deepcopy(optimizer), headroom,
            _rpo_comparison_reference_samples(max_reference_samples))
    end
end

"""
    RRTConnectBezierPlanner(settings=RPORRTConnectSettings(); config, optimizer, headroom, max_reference_samples)

The `:rrt_connect_bezier` comparison mode: the RRT-Connect search of
`RRTConnectPlanner`, with its path refitted as a Bezier curve under `config`.
The keywords are those of `RRTConnectPlanner`.
"""
struct RRTConnectBezierPlanner <: AbstractRPOComparisonPlanner
    config::RPOPSOConfig
    settings::RPORRTConnectSettings
    optimizer::RPOTrajectoryOptimizerSettings
    headroom::RPOPlanningHeadroom
    max_reference_samples::Int
    function RRTConnectBezierPlanner(settings::RPORRTConnectSettings=RPORRTConnectSettings();
            config::RPOPSOConfig=_rpo_comparison_default_config(),
            optimizer::RPOTrajectoryOptimizerSettings=RPOTrajectoryOptimizerSettings(),
            headroom::RPOPlanningHeadroom=RPOPlanningHeadroom(), max_reference_samples::Integer=100_000)
        new(_rpo_comparison_owned_config(config), deepcopy(settings), deepcopy(optimizer), headroom,
            _rpo_comparison_reference_samples(max_reference_samples))
    end
end

"""
    RRTStarPlanner(settings=RPORRTStarSettings(); config, optimizer, headroom, max_reference_samples)

The `:rrt_star` comparison mode: RRT* with the search `settings`, scored, refined
and retimed with `config`. The keywords are those of `RRTConnectPlanner`.
"""
struct RRTStarPlanner <: AbstractRPOComparisonPlanner
    config::RPOPSOConfig
    settings::RPORRTStarSettings
    optimizer::RPOTrajectoryOptimizerSettings
    headroom::RPOPlanningHeadroom
    max_reference_samples::Int
    function RRTStarPlanner(settings::RPORRTStarSettings=RPORRTStarSettings();
            config::RPOPSOConfig=_rpo_comparison_default_config(),
            optimizer::RPOTrajectoryOptimizerSettings=RPOTrajectoryOptimizerSettings(),
            headroom::RPOPlanningHeadroom=RPOPlanningHeadroom(), max_reference_samples::Integer=100_000)
        new(_rpo_comparison_owned_config(config), deepcopy(settings), deepcopy(optimizer), headroom,
            _rpo_comparison_reference_samples(max_reference_samples))
    end
end

"""
    CHOMPPlanner(settings=RPOCHOMPSettings(n_iters=config.n_iters); config, optimizer, seed, headroom, max_reference_samples)

The `:chomp` comparison mode: a CHOMP-like covariant gradient optimizer with
`settings`, under the cost of `config` and the obstacle, stopping and iteration
settings of `optimizer`. The optimizer is deterministic, but its initial path is
an RRT-Connect search with the `seed` settings (an `RPORRTConnectSettings`) drawn
from the planner's random stream, unless the comparison supplies a warm start.
"""
struct CHOMPPlanner <: AbstractRPOComparisonPlanner
    config::RPOPSOConfig
    settings::RPOCHOMPSettings
    optimizer::RPOTrajectoryOptimizerSettings
    seed::RPORRTConnectSettings
    headroom::RPOPlanningHeadroom
    max_reference_samples::Int
    function CHOMPPlanner(settings::Union{Nothing, RPOCHOMPSettings}=nothing;
            config::RPOPSOConfig=_rpo_comparison_default_config(),
            optimizer::RPOTrajectoryOptimizerSettings=RPOTrajectoryOptimizerSettings(),
            seed::RPORRTConnectSettings=RPORRTConnectSettings(),
            headroom::RPOPlanningHeadroom=RPOPlanningHeadroom(), max_reference_samples::Integer=100_000)
        chosen = settings === nothing ? RPOCHOMPSettings(n_iters=config.n_iters) : settings
        new(_rpo_comparison_owned_config(config), deepcopy(chosen), deepcopy(optimizer), deepcopy(seed),
            headroom, _rpo_comparison_reference_samples(max_reference_samples))
    end
end

"""
    STOMPPlanner(settings=RPOSTOMPSettings(n_iters=config.n_iters); config, optimizer, seed, headroom, max_reference_samples)

The `:stomp` comparison mode: a STOMP-like stochastic optimizer with `settings`,
started from an RRT-Connect path as for `CHOMPPlanner`. The seed search draws from
the planner's random stream before the rollouts do. The keywords are those of
`CHOMPPlanner`.
"""
struct STOMPPlanner <: AbstractRPOComparisonPlanner
    config::RPOPSOConfig
    settings::RPOSTOMPSettings
    optimizer::RPOTrajectoryOptimizerSettings
    seed::RPORRTConnectSettings
    headroom::RPOPlanningHeadroom
    max_reference_samples::Int
    function STOMPPlanner(settings::Union{Nothing, RPOSTOMPSettings}=nothing;
            config::RPOPSOConfig=_rpo_comparison_default_config(),
            optimizer::RPOTrajectoryOptimizerSettings=RPOTrajectoryOptimizerSettings(),
            seed::RPORRTConnectSettings=RPORRTConnectSettings(),
            headroom::RPOPlanningHeadroom=RPOPlanningHeadroom(), max_reference_samples::Integer=100_000)
        chosen = settings === nothing ? RPOSTOMPSettings(n_iters=config.n_iters) : settings
        new(_rpo_comparison_owned_config(config), deepcopy(chosen), deepcopy(optimizer), deepcopy(seed),
            headroom, _rpo_comparison_reference_samples(max_reference_samples))
    end
end

# Canonical identifiers stay the report keys and labels of each mode.
normalize_rpo_comparison_planner_type(::HYPRComparisonPlanner) = :hypr
normalize_rpo_comparison_planner_type(::UnrefinedPSOPlanner) = :pso_unrefined
normalize_rpo_comparison_planner_type(::RRTConnectPlanner) = :rrt_connect
normalize_rpo_comparison_planner_type(::RRTConnectBezierPlanner) = :rrt_connect_bezier
normalize_rpo_comparison_planner_type(::RRTStarPlanner) = :rrt_star
normalize_rpo_comparison_planner_type(::CHOMPPlanner) = :chomp
normalize_rpo_comparison_planner_type(::STOMPPlanner) = :stomp

# The comparison's iteration policy: a finite runtime limit selects
# `runtime_max_iters`, otherwise `match_hypr_iters` takes HYPR's count.
function _rpo_comparison_iterations(optimizer::RPOTrajectoryOptimizerSettings, base_cfg::RPOPSOConfig,
        runtime_limit_s::Real, default_iters::Integer)
    return isfinite(Float64(runtime_limit_s)) && optimizer.runtime_max_iters > 0 ?
        optimizer.runtime_max_iters :
        (optimizer.match_hypr_iters ? base_cfg.n_iters : default_iters)
end

"""Copy RRT-Connect settings with a caller-selected iteration cap."""
function _rpo_comparison_rrt_connect_settings(src::RPORRTConnectSettings, n_iters::Integer)
    return _with_fields(src; n_iters=Int(n_iters))
end

"""Copy RRT* settings with a caller-selected iteration cap."""
function _rpo_comparison_rrt_star_settings(src::RPORRTStarSettings, n_iters::Integer)
    return _with_fields(src; n_iters=Int(n_iters))
end

# The RRT-Connect path that warm-starts CHOMP and STOMP, drawn from `rng`.
function _rpo_comparison_seed_path(planner::Union{CHOMPPlanner, STOMPPlanner}, start_rtn, goal_rtn, geometry,
        base_cfg::RPOPSOConfig; safe_distance_m::Real, runtime_limit_s::Real, rng)
    n_iters = _rpo_comparison_iterations(planner.optimizer, base_cfg, runtime_limit_s, planner.seed.n_iters)
    plan = rpo_rrt_connect_plan_path(
        start_rtn,
        goal_rtn,
        geometry,
        base_cfg;
        safe_distance_m=safe_distance_m,
        settings=_rpo_comparison_rrt_connect_settings(planner.seed, n_iters),
        max_runtime_s=runtime_limit_s,
        rng=rng,
    )
    return plan.path
end

"""
    _rpo_comparison_planner_path(planner, start_rtn, goal_rtn, geometry, base_cfg;
        safe_distance_m, rng, runtime_limit_s, warmstart_path) -> (plan, iterations)

Run one comparison mode's kernels with `base_cfg`, the planner's configuration
already carrying the caller's safety distance (and, on the contract path, the
request's limits). Both the comparison harness and `plan_rpo!` use this, so the
kernels, settings and draw order are the same on both routes.
"""
function _rpo_comparison_planner_path(planner::HYPRComparisonPlanner, start_rtn, goal_rtn, geometry,
        base_cfg::RPOPSOConfig; safe_distance_m::Real, rng, runtime_limit_s::Real=Inf, warmstart_path=nothing)
    plan = rpo_pso_plan_path(start_rtn, goal_rtn, geometry, base_cfg; safe_distance_m=safe_distance_m, rng=rng)
    return plan, length(plan.cost_history)
end

function _rpo_comparison_planner_path(planner::UnrefinedPSOPlanner, start_rtn, goal_rtn, geometry,
        base_cfg::RPOPSOConfig; safe_distance_m::Real, rng, runtime_limit_s::Real=Inf, warmstart_path=nothing)
    unrefined_cfg = rpo_pso_config(base_cfg; refinement_enable=false,
        adaptive_enable=false, schedule_enable=false, stagnation_learning_enable=false,
        rrt_warmstart_enable=false,
        cull_enable=false, reexplore_enable=false)
    plan = rpo_pso_plan_path(start_rtn, goal_rtn, geometry, unrefined_cfg; safe_distance_m=safe_distance_m, rng=rng)
    return plan, length(plan.cost_history)
end

function _rpo_comparison_planner_path(planner::RRTConnectPlanner, start_rtn, goal_rtn, geometry,
        base_cfg::RPOPSOConfig; safe_distance_m::Real, rng, runtime_limit_s::Real=Inf, warmstart_path=nothing)
    n_iters = _rpo_comparison_iterations(planner.optimizer, base_cfg, runtime_limit_s, planner.settings.n_iters)
    plan = rpo_rrt_connect_plan_path(
        start_rtn,
        goal_rtn,
        geometry,
        base_cfg;
        safe_distance_m=safe_distance_m,
        settings=_rpo_comparison_rrt_connect_settings(planner.settings, n_iters),
        max_runtime_s=runtime_limit_s,
        rng=rng,
    )
    return plan, plan.iterations
end

function _rpo_comparison_planner_path(planner::RRTConnectBezierPlanner, start_rtn, goal_rtn, geometry,
        base_cfg::RPOPSOConfig; safe_distance_m::Real, rng, runtime_limit_s::Real=Inf, warmstart_path=nothing)
    n_iters = _rpo_comparison_iterations(planner.optimizer, base_cfg, runtime_limit_s, planner.settings.n_iters)
    plan = rpo_rrt_connect_bezier_plan_path(
        start_rtn,
        goal_rtn,
        geometry,
        base_cfg;
        safe_distance_m=safe_distance_m,
        settings=_rpo_comparison_rrt_connect_settings(planner.settings, n_iters),
        max_runtime_s=runtime_limit_s,
        rng=rng,
    )
    return plan, plan.iterations
end

function _rpo_comparison_planner_path(planner::RRTStarPlanner, start_rtn, goal_rtn, geometry,
        base_cfg::RPOPSOConfig; safe_distance_m::Real, rng, runtime_limit_s::Real=Inf, warmstart_path=nothing)
    n_iters = _rpo_comparison_iterations(planner.optimizer, base_cfg, runtime_limit_s, planner.settings.n_iters)
    plan = rpo_rrt_star_plan_path(
        start_rtn,
        goal_rtn,
        geometry,
        base_cfg;
        safe_distance_m=safe_distance_m,
        settings=_rpo_comparison_rrt_star_settings(planner.settings, n_iters),
        max_runtime_s=runtime_limit_s,
        rng=rng,
    )
    return plan, plan.iterations
end

function _rpo_comparison_planner_path(planner::CHOMPPlanner, start_rtn, goal_rtn, geometry,
        base_cfg::RPOPSOConfig; safe_distance_m::Real, rng, runtime_limit_s::Real=Inf, warmstart_path=nothing)
    src = planner.settings
    settings = _with_fields(src; n_iters=_rpo_comparison_iterations(planner.optimizer, base_cfg, runtime_limit_s, src.n_iters))
    initial_path = warmstart_path === nothing ? _rpo_comparison_seed_path(
        planner,
        start_rtn,
        goal_rtn,
        geometry,
        base_cfg;
        safe_distance_m=safe_distance_m,
        runtime_limit_s=runtime_limit_s,
        rng=rng,
    ) : Matrix{Float64}(warmstart_path)
    plan = rpo_chomp_plan_path(
        start_rtn,
        goal_rtn,
        geometry,
        base_cfg;
        safe_distance_m=safe_distance_m,
        settings=settings,
        optimizer=planner.optimizer,
        max_runtime_s=runtime_limit_s,
        initial_path=initial_path,
    )
    return plan, plan.iterations
end

function _rpo_comparison_planner_path(planner::STOMPPlanner, start_rtn, goal_rtn, geometry,
        base_cfg::RPOPSOConfig; safe_distance_m::Real, rng, runtime_limit_s::Real=Inf, warmstart_path=nothing)
    src = planner.settings
    settings = _with_fields(src; n_iters=_rpo_comparison_iterations(planner.optimizer, base_cfg, runtime_limit_s, src.n_iters))
    initial_path = warmstart_path === nothing ? _rpo_comparison_seed_path(
        planner,
        start_rtn,
        goal_rtn,
        geometry,
        base_cfg;
        safe_distance_m=safe_distance_m,
        runtime_limit_s=runtime_limit_s,
        rng=rng,
    ) : Matrix{Float64}(warmstart_path)
    plan = rpo_stomp_plan_path(
        start_rtn,
        goal_rtn,
        geometry,
        base_cfg;
        safe_distance_m=safe_distance_m,
        settings=settings,
        optimizer=planner.optimizer,
        max_runtime_s=runtime_limit_s,
        rng=rng,
        initial_path=initial_path,
    )
    return plan, plan.iterations
end
