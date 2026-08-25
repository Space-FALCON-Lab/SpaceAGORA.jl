# ParallelProfiles names arrive from the SimulationCampaigns module header so
# that monte_carlo.jl (included first) can annotate its adaptive keywords.
using ..SimulationModel: SimulationConfiguration, EnvironmentModels, ParallelPolicy

const _CAMPAIGN_OUTER_ROUTE_STATE = OuterRouteState()

# Cross-session persistence for the outer-route bandit.
#
# `save_outer_route_state` and `load_outer_route_state!` have existed and been
# exported since the state was written, and nothing in `src/` ever called them,
# so the bandit started cold in every process. That is not a small loss on the
# workloads it matters for: `select_outer_route!` refuses to exploit until every
# candidate has `adaptive_min_samples` observations, so a fresh process spends
# its first samples on `:none` and `:threads` before it will use `:process` --
# even when `default_outer_route` already answers `:process` correctly and the
# same machine learned that in the previous run. Measured on
# montecarlo_heavy_aerobraking at 64 samples, the adaptive profiles ran ~19%
# behind the pinned process route, and forced cold exploration is a direct
# contributor.
#
# Persistence is gated on the same knob the inner hints use, which R5 already
# sets (`persistent_state_persist=true` in profile_definitions.jl), so this
# turns on for exactly the profile that declares it and stays off for the
# static baselines that must not drift between runs.
const _CAMPAIGN_ROUTE_STATE_LOADED = Ref(false)
const _CAMPAIGN_ROUTE_STATE_ATEXIT = Ref(false)
const _CAMPAIGN_ROUTE_STATE_LOCK = ReentrantLock()

@inline function _campaign_route_state_persist_enabled()::Bool
    return ParallelPolicy.persistent_hints_persist_enabled()
end

"""
    campaign_route_state_path() -> String

File backing the persisted outer-route history.

Keyed by profile, machine label and thread count for the same reason the inner
policy state is: routing statistics gathered under a different profile or a
different thread budget describe a different machine as far as the bandit is
concerned, and replaying them would be worse than starting cold.
`SPACEAGORA_OUTER_ROUTE_STATE_PATH` overrides.
"""
function campaign_route_state_path()::String
    override = strip(get(ENV, "SPACEAGORA_OUTER_ROUTE_STATE_PATH", ""))
    if !isempty(override)
        return normpath(isabspath(override) ? override : joinpath(pwd(), override))
    end
    profile = ParallelPolicy._safe_token(get(ENV, "SPACEAGORA_PARALLEL_PROFILE", "default"))
    machine = ParallelPolicy._safe_token(get(ENV, "SPACEAGORA_PERF_MACHINE_LABEL", "default"))
    return normpath(joinpath(
        pwd(), "output", "parallel_policy_state",
        "outer_route_state_$(profile)_$(machine)_t$(Base.Threads.nthreads()).toml",
    ))
end

"""
    ensure_campaign_route_state_loaded!(state = campaign_outer_route_state())

Load persisted routing history into `state` once per process, and register an
`atexit` hook to write it back.

Loading is additive (`replace=false`): a session's own observations are merged
with the persisted ones rather than discarding either. Failures are swallowed
deliberately -- a missing, unreadable, or malformed state file must degrade the
router to cold-start behaviour, never break a campaign that would otherwise run.
"""
function ensure_campaign_route_state_loaded!(state::OuterRouteState = campaign_outer_route_state())::Nothing
    _campaign_route_state_persist_enabled() || return nothing
    lock(_CAMPAIGN_ROUTE_STATE_LOCK) do
        if !_CAMPAIGN_ROUTE_STATE_LOADED[]
            _CAMPAIGN_ROUTE_STATE_LOADED[] = true
            path = campaign_route_state_path()
            try
                load_outer_route_state!(state, path; replace=false)
            catch err
                @debug "Outer-route state could not be loaded; starting cold." path exception=err
            end
        end
        if !_CAMPAIGN_ROUTE_STATE_ATEXIT[]
            _CAMPAIGN_ROUTE_STATE_ATEXIT[] = true
            atexit(() -> save_campaign_route_state())
        end
        return nothing
    end
    return nothing
end

"""
    save_campaign_route_state(state = campaign_outer_route_state())

Write the accumulated routing history back to [`campaign_route_state_path`](@ref).
Returns without writing when persistence is disabled, and never throws: a state
file that cannot be written must cost the router its memory, not the run.
"""
function save_campaign_route_state(state::OuterRouteState = campaign_outer_route_state())::Nothing
    _campaign_route_state_persist_enabled() || return nothing
    try
        save_outer_route_state(
            state,
            campaign_route_state_path();
            metadata=Dict{String, Any}(
                "profile" => get(ENV, "SPACEAGORA_PARALLEL_PROFILE", "default"),
                "threads" => Base.Threads.nthreads(),
            ),
        )
    catch err
        @debug "Outer-route state could not be saved." exception=err
    end
    return nothing
end

"""
    reset_campaign_route_state_persistence!()

Forget that persisted state was loaded, so the next campaign reloads it. Exists
for tests; a normal process loads once and keeps the state in memory.
"""
function reset_campaign_route_state_persistence!()::Nothing
    lock(_CAMPAIGN_ROUTE_STATE_LOCK) do
        _CAMPAIGN_ROUTE_STATE_LOADED[] = false
    end
    return nothing
end

"""
    campaign_outer_route_state() -> OuterRouteState

Return the process-global [`OuterRouteState`](@ref) that the campaign runners
consult and update when called with `threads=:auto`.

Repeated campaigns in the same session share this state, so route selection
converges to the empirically fastest serial/threaded allocation per workload
signature. Inspect it with [`outer_route_stats_snapshot`](@ref), reset it with
[`reset_outer_route_state!`](@ref), or persist it across sessions with
`save_outer_route_state` and `load_outer_route_state!`. Pass `route_state` to
[`run_monte_carlo`](@ref) or [`run_constellation_ensemble`](@ref) to use an
isolated state instead.
"""
campaign_outer_route_state()::OuterRouteState = _CAMPAIGN_OUTER_ROUTE_STATE

function _campaign_density_family(density_model)::String
    if density_model isa EnvironmentModels.NoAtmosphereModel
        return "none"
    elseif density_model isa EnvironmentModels.GRAMAtmosphereModelSurrogate
        return "gram_surrogate"
    elseif density_model isa EnvironmentModels.GRAMAtmosphereModel
        return "gram_point"
    elseif density_model isa EnvironmentModels.ExponentialAtmosphereModel
        return "exponential"
    elseif density_model isa EnvironmentModels.PolynomialFitAtmosphereModel
        return "polyfit"
    elseif density_model isa EnvironmentModels.NRLMSISE00AtmosphereModel
        return "nrlmsise00"
    end
    return lowercase(string(nameof(typeof(density_model))))
end

"""
    campaign_route_features(; samples, n_sats=1, density_family="unknown",
                            mission_time_s=0.0) -> OuterRouteFeatures
    campaign_route_features(args::SimulationConfiguration; samples,
                            n_sats=length(args.dynamics_model.spacecraft)) -> OuterRouteFeatures

Build the [`OuterRouteFeatures`](@ref) describing a campaign workload for
adaptive outer-route selection (`threads=:auto` in [`run_monte_carlo`](@ref)
and [`run_constellation_ensemble`](@ref)).

`samples` is the campaign sample count (Monte Carlo seeds or ensemble members)
and `n_sats` the per-sample satellite count. The `SimulationConfiguration`
method additionally derives the density-model family, mission length,
orientation flag, and effector counts from `args`. Feedback recorded via
[`record_outer_route_feedback!`](@ref) is bucketed by
[`outer_route_signature`](@ref), so campaigns with the same shape share
empirical runtime statistics.

Campaign features always carry the `"montecarlo"` routing category: the other
categories' default-route rules answer for a single coupled simulation (inner
versus outer split), not for a sample fan-out, and can leave the adaptive
runner without a feasible default. The adaptive path also normalizes any
caller-provided `route_features` to this category before routing and recording.
"""
function campaign_route_features(;
    samples::Integer,
    n_sats::Integer=1,
    density_family::AbstractString="unknown",
    mission_time_s::Real=0.0
)::OuterRouteFeatures
    samples >= 0 || throw(ArgumentError("campaign_route_features samples must be >= 0; got $(samples)."))
    n_sats >= 1 || throw(ArgumentError("campaign_route_features n_sats must be >= 1; got $(n_sats)."))
    return OuterRouteFeatures(
        category="montecarlo",
        n_sats=Int(n_sats),
        mission_time_s=Float64(mission_time_s),
        density_family=String(density_family),
        montecarlo_samples=Int(samples)
    )
end

function campaign_route_features(
    args::SimulationConfiguration;
    samples::Integer,
    n_sats::Integer=length(args.dynamics_model.spacecraft)
)::OuterRouteFeatures
    samples >= 0 || throw(ArgumentError("campaign_route_features samples must be >= 0; got $(samples)."))
    n_sats >= 1 || throw(ArgumentError("campaign_route_features n_sats must be >= 1; got $(n_sats)."))
    spacecraft = args.dynamics_model.spacecraft
    per_sat_links = isempty(spacecraft) ? 1 : maximum(sc -> max(1, length(sc.links)), spacecraft)
    # A per-sample view (n_sats below the constellation size, e.g. ensemble
    # members) carries one satellite's links; the full-constellation view sums them.
    n_links = if Int(n_sats) >= length(spacecraft)
        max(1, sum(sc -> length(sc.links), spacecraft; init=0))
    else
        per_sat_links
    end
    control_effector_count = length(args.control_model.control_effectors)
    density_family = _campaign_density_family(args.environment_model.density_model)
    return OuterRouteFeatures(
        category="montecarlo",
        n_sats=Int(n_sats),
        n_links=n_links,
        max_links_per_sat=per_sat_links,
        mission_time_s=Float64(args.mission_configuration.mission_time),
        has_control=control_effector_count > 0,
        orientation_on=args.mission_configuration.orientation_sim,
        density_family=density_family,
        dt_max_orbit_s=Float64(args.integration_tolerances.dt_max_orbit),
        gram_surrogate_enabled=density_family == "gram_surrogate",
        control_effector_count=control_effector_count,
        dynamic_effector_count=length(args.dynamics_model.dynamic_effectors),
        montecarlo_samples=Int(samples)
    )
end

function _campaign_route_tuning()::OuterRouteTuning
    # Campaign runners can route to a real :process backend (ParallelProcess),
    # so the default thresholds apply unmodified.
    return OuterRouteTuning()
end

function _campaign_features_for_routing(f::OuterRouteFeatures, samples::Int)::OuterRouteFeatures
    if f.montecarlo_samples == samples && f.category == "montecarlo"
        return f
    end
    fields = NamedTuple{fieldnames(OuterRouteFeatures)}(
        ntuple(i -> getfield(f, i), fieldcount(OuterRouteFeatures))
    )
    # Force the montecarlo routing category: other categories' default-route
    # rules describe a single coupled simulation and can answer :process for a
    # shape whose candidate set is [:none, :threads], which select_outer_route!
    # then clamps to a serial default on medium/large machines.
    return OuterRouteFeatures(; fields..., category="montecarlo", montecarlo_samples=samples)
end

function _campaign_route_plan(
    features::OuterRouteFeatures,
    n_samples::Int;
    state::OuterRouteState,
    tuning::OuterRouteTuning
)
    if ParallelPolicy.outer_parallel_active()
        # A higher-level campaign (or benchmark harness) already owns the outer
        # split; running nested workers would oversubscribe the pool, and the
        # contended timings would poison the shared route statistics — run
        # serially and skip both selection and feedback.
        return (route=:none, threads=1, inner_thread_budget=1, record=false)
    end
    # Merge any persisted history before the first selection, so a repeat run on
    # the same machine exploits what the last one learned instead of re-paying
    # for cold exploration.
    ensure_campaign_route_state_loaded!(state)
    threads_available = Base.Threads.nthreads() > 1
    route = select_outer_route!(
        state,
        features;
        tuning=tuning,
        machine_class=ParallelProfiles._machine_parallel_class(),
        threads_available=threads_available,
        parallel_enabled=true
    )
    # Process workers run --threads=1 each and don't share the coordinator's
    # thread pool, so they're sized off tuning.process_max_workers
    # (Sys.CPU_THREADS by default), not Threads.nthreads() like the thread route.
    workers = if route === :process
        max(1, min(n_samples, tuning.process_max_workers))
    elseif route === :threads
        max(1, min(n_samples, Base.Threads.nthreads()))
    else
        1
    end
    inner_thread_budget = max(1, fld(Base.Threads.nthreads(), workers))
    return (route=route, threads=workers, inner_thread_budget=inner_thread_budget, record=true)
end

function _run_campaign_with_route_env(f, spec::MonteCarloSpec, plan)
    worker_count = min(spec.threads, length(spec.seeds))
    worker_count > 1 || return run_monte_carlo(f, spec)
    if plan.route === :process
        # Bypasses run_monte_carlo (which validates its thread count against
        # Threads.nthreads() -- not meaningful here, since process workers
        # aren't Julia threads) and dispatches straight to the process pool.
        pool = campaign_process_pool()
        # warmup_fn reuses f itself (the exact closure about to be dispatched
        # for real) so a newly-added worker's large one-time JIT/specialization
        # cost (see ensure_process_workers!'s docstring) is paid here, once,
        # rather than silently inside this call's own timed dispatch below.
        worker_ids = ensure_process_workers!(pool, worker_count; warmup_fn=() -> f(first(spec.seeds)))
        active_workers = worker_ids[1:min(worker_count, length(worker_ids))]
        start_ns = time_ns()
        samples = _run_monte_carlo_process(f, spec.seeds, spec, active_workers)
        elapsed_s = (time_ns() - start_ns) / 1.0e9
        spec.fail_fast && _throw_first_monte_carlo_failure(samples)
        return MonteCarloResult(samples, elapsed_s, length(active_workers))
    end
    env_pairs = Pair{String, String}["SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1"]
    if isempty(strip(get(ENV, "SPACEAGORA_INNER_THREAD_BUDGET", "")))
        # Split the thread pool between the outer workers and each sample's
        # inner parallelism instead of oversubscribing; an explicit user budget
        # always wins.
        push!(env_pairs, "SPACEAGORA_INNER_THREAD_BUDGET" => string(plan.inner_thread_budget))
    end
    return withenv(env_pairs...) do
        run_monte_carlo(f, spec)
    end
end

function _record_campaign_route_feedback!(
    state::OuterRouteState,
    features::OuterRouteFeatures,
    route::Symbol,
    result::MonteCarloResult;
    tuning::OuterRouteTuning
)::Nothing
    n_samples = length(result.samples)
    n_samples <= 0 && return nothing
    successes = length(result.successful)
    failures = length(result.failed)
    # Amortized campaign wall time per sample, not per-sample latency: threaded
    # samples overlap, so summing individual elapsed_s would charge the route
    # for concurrency instead of crediting its throughput.
    per_sample_s = result.elapsed_s / n_samples
    record_outer_route_feedback!(
        state,
        features;
        route=route,
        successes=successes,
        failures=failures,
        elapsed_success_s=per_sample_s * successes,
        elapsed_success_sq_sum_s=per_sample_s^2 * successes,
        tuning=tuning
    )
    return nothing
end

function _run_campaign_adaptive(
    f,
    seeds;
    fail_fast::Bool,
    features::OuterRouteFeatures,
    state::OuterRouteState,
    tuning::OuterRouteTuning
)::MonteCarloResult
    seed_values = collect(seeds)
    isempty(seed_values) && return MonteCarloResult(MonteCarloSampleResult[], 0.0, 0)
    routed_features = _campaign_features_for_routing(features, length(seed_values))
    plan = _campaign_route_plan(routed_features, length(seed_values); state=state, tuning=tuning)
    spec = MonteCarloSpec(seeds=seed_values, threads=plan.threads, fail_fast=fail_fast)
    result = _run_campaign_with_route_env(f, spec, plan)
    if plan.record
        _record_campaign_route_feedback!(state, routed_features, plan.route, result; tuning=tuning)
    end
    return result
end

function _run_monte_carlo_adaptive(
    f,
    seeds;
    fail_fast::Bool,
    route_features::Union{Nothing, OuterRouteFeatures},
    route_state::Union{Nothing, OuterRouteState},
    route_tuning::Union{Nothing, OuterRouteTuning}
)::MonteCarloResult
    features = route_features === nothing ? campaign_route_features(samples=0) : route_features
    state = route_state === nothing ? campaign_outer_route_state() : route_state
    tuning = route_tuning === nothing ? _campaign_route_tuning() : route_tuning
    return _run_campaign_adaptive(f, seeds; fail_fast=fail_fast, features=features, state=state, tuning=tuning)
end
