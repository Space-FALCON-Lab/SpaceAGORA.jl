# ParallelProfiles names arrive from the SimulationCampaigns module header so
# that monte_carlo.jl (included first) can annotate its adaptive keywords.
using ..SimulationModel: SimulationConfiguration, EnvironmentModels

const _CAMPAIGN_OUTER_ROUTE_STATE = OuterRouteState()

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
                            mission_time_s=0.0, category="montecarlo") -> OuterRouteFeatures
    campaign_route_features(args::SimulationConfiguration; samples,
                            n_sats=length(args.dynamics_model.spacecraft),
                            category="montecarlo") -> OuterRouteFeatures

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
"""
function campaign_route_features(;
    samples::Integer,
    n_sats::Integer=1,
    density_family::AbstractString="unknown",
    mission_time_s::Real=0.0,
    category::AbstractString="montecarlo"
)::OuterRouteFeatures
    samples >= 0 || throw(ArgumentError("campaign_route_features samples must be >= 0; got $(samples)."))
    n_sats >= 1 || throw(ArgumentError("campaign_route_features n_sats must be >= 1; got $(n_sats)."))
    return OuterRouteFeatures(
        category=String(category),
        n_sats=Int(n_sats),
        mission_time_s=Float64(mission_time_s),
        density_family=String(density_family),
        montecarlo_samples=Int(samples)
    )
end

function campaign_route_features(
    args::SimulationConfiguration;
    samples::Integer,
    n_sats::Integer=length(args.dynamics_model.spacecraft),
    category::AbstractString="montecarlo"
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
        category=String(category),
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
    # Campaign runners execute in-process (serial or threaded worker tasks); a
    # :process backend does not exist here, so suppress the process-preferring
    # Monte Carlo and constellation rules up front.
    return OuterRouteTuning(
        spice_constellation_process_enabled=false,
        mc_process_min_samples=typemax(Int),
        mc_process_min_mission_s=Inf
    )
end

function _with_montecarlo_samples(f::OuterRouteFeatures, samples::Int)::OuterRouteFeatures
    f.montecarlo_samples == samples && return f
    fields = NamedTuple{fieldnames(OuterRouteFeatures)}(
        ntuple(i -> getfield(f, i), fieldcount(OuterRouteFeatures))
    )
    return OuterRouteFeatures(; fields..., montecarlo_samples=samples)
end

function _campaign_route_plan(
    features::OuterRouteFeatures,
    n_samples::Int;
    state::OuterRouteState,
    tuning::OuterRouteTuning
)
    threads_available = Base.Threads.nthreads() > 1
    chosen = select_outer_route!(
        state,
        features;
        tuning=tuning,
        machine_class=ParallelProfiles._machine_parallel_class(),
        threads_available=threads_available,
        parallel_enabled=true
    )
    # The routing layer can still answer :process (native-GRAM point densities
    # bypass the tuning thresholds); degrade to the closest supported route and
    # record feedback under the route that actually ran.
    route = chosen === :process ? (threads_available ? :threads : :none) : chosen
    workers = route === :threads ? min(n_samples, Base.Threads.nthreads()) : 1
    workers = max(1, workers)
    inner_thread_budget = max(1, fld(Base.Threads.nthreads(), workers))
    return (route=route, threads=workers, inner_thread_budget=inner_thread_budget)
end

function _run_campaign_with_route_env(f, spec::MonteCarloSpec, plan)
    worker_count = min(spec.threads, length(spec.seeds))
    worker_count > 1 || return run_monte_carlo(f, spec)
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
    routed_features = _with_montecarlo_samples(features, length(seed_values))
    plan = _campaign_route_plan(routed_features, length(seed_values); state=state, tuning=tuning)
    spec = MonteCarloSpec(seeds=seed_values, threads=plan.threads, fail_fast=fail_fast)
    result = _run_campaign_with_route_env(f, spec, plan)
    _record_campaign_route_feedback!(state, routed_features, plan.route, result; tuning=tuning)
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
