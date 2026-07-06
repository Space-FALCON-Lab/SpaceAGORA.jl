using ..SimulationModel: SimulationConfiguration, SimulationSettings, DynamicsModel, SpacecraftModel
import ..SimulationEngine

function _ensemble_member_settings(settings::SimulationSettings, member_tag::String)::SimulationSettings
    # Checkpoints can be read or written whenever checkpointing OR resume is on,
    # and an empty checkpoint_directory falls back to results_directory/checkpoints —
    # so that fallback case must split results_directory even when results=false.
    checkpoint_active = settings.checkpoint_enabled || settings.resume_from_checkpoint
    explicit_checkpoint_dir = !isempty(strip(settings.checkpoint_directory))
    needs_results_split = settings.results || (checkpoint_active && !explicit_checkpoint_dir)
    needs_checkpoint_split = checkpoint_active && explicit_checkpoint_dir
    if !needs_results_split && !needs_checkpoint_split
        return settings
    end
    return SimulationSettings(
        results=settings.results,
        verbose=settings.verbose,
        results_directory=needs_results_split ?
            joinpath(settings.results_directory, member_tag) : settings.results_directory,
        generate_plots=settings.generate_plots,
        generate_filenames=settings.generate_filenames,
        normalize=settings.normalize,
        save_csv=settings.save_csv,
        checkpoint_enabled=settings.checkpoint_enabled,
        checkpoint_interval_s=settings.checkpoint_interval_s,
        checkpoint_directory=needs_checkpoint_split ?
            joinpath(settings.checkpoint_directory, member_tag) : settings.checkpoint_directory,
        resume_from_checkpoint=settings.resume_from_checkpoint
    )
end

function _ensemble_member_configuration(
    args::SimulationConfiguration,
    spacecraft::SpacecraftModel,
    member_tag::String
)::SimulationConfiguration
    return SimulationConfiguration(
        file_paths=args.file_paths,
        simulation_settings=_ensemble_member_settings(args.simulation_settings, member_tag),
        mission_configuration=args.mission_configuration,
        environment_model=args.environment_model,
        dynamics_model=DynamicsModel([spacecraft], args.dynamics_model.dynamic_effectors),
        guidance_model=args.guidance_model,
        navigation_model=args.navigation_model,
        control_model=args.control_model,
        initial_time=args.initial_time,
        integration_tolerances=args.integration_tolerances,
        solver_config=args.solver_config
    )
end

function _validate_ensemble_uncoupled(args::SimulationConfiguration, allow_gnc_effectors::Bool)
    allow_gnc_effectors && return nothing
    coupled_surfaces = String[]
    isempty(args.guidance_model.guidance_effectors) || push!(coupled_surfaces, "guidance_effectors")
    isempty(args.navigation_model.navigation_effectors) || push!(coupled_surfaces, "navigation_effectors")
    isempty(args.control_model.control_effectors) || push!(coupled_surfaces, "control_effectors")
    isempty(coupled_surfaces) && return nothing
    throw(ArgumentError(
        "run_constellation_ensemble requires an uncoupled constellation, but the configuration has " *
        "non-empty $(join(coupled_surfaces, ", ")). GNC effectors may couple satellites (for example " *
        "RPO planners or coordinated maneuvers), and each ensemble member propagates alone. If every " *
        "effector acts on a single satellite only, opt in with allow_gnc_effectors=true; otherwise use " *
        "the monolithic run_simulation path."
    ))
end

"""
    run_constellation_ensemble(args::SimulationConfiguration; threads=1, fail_fast=false,
                               allow_gnc_effectors=false, route_state=nothing,
                               route_tuning=nothing, run_kwargs...) -> MonteCarloResult

Propagate each spacecraft of a multi-satellite `SimulationConfiguration` as an
independent single-satellite simulation, parallelised across satellites with the
[`run_monte_carlo`](@ref) worker-task runner.

For constellations whose members do not interact dynamically, this outer-level split
replaces the monolithic coupled state vector, which has two costs: per-timestep thread
dispatch across satellites (whose overhead dominates when per-satellite work is light),
and a shared adaptive timestep that forces every satellite to the global minimum step.
Each ensemble member instead dispatches to a worker once for its entire propagation and
keeps its own adaptive step size.

Sample `i` of the returned [`MonteCarloResult`](@ref) holds the `run_simulation` return
value for spacecraft `i` (in `args.dynamics_model.spacecraft` order) in `value`. When
result saving, checkpointing, or checkpoint resume is enabled, each member uses
per-satellite `sat_<index>_id_<spacecraft id>` subdirectories so members do not read or
clobber each other's files.

`threads` caps worker tasks exactly as in [`run_monte_carlo`](@ref); Julia must be started
with at least that many threads. While more than one worker is active the runner sets
`SPACEAGORA_OUTER_PARALLEL_ACTIVE=1` so inner thread policies yield to the outer split.

Isolation: when more than one worker is active, each member solves a `deepcopy` of its
configuration taken inside its worker task, because model state such as
`environment_model.planet.L_PI` is mutated during a solve and must never be shared
between concurrent members. The copy replaces (not adds to) `run_simulation`'s own
`isolate_state` copy, so the total cost matches a plain isolated run; an explicit
`isolate_state=false` is ignored under parallel execution. With `threads=1` the
keyword is forwarded to [`run_simulation`](@ref) unchanged.

Pass `threads=:auto` to let the outer-route bandit choose serial or threaded execution
from empirical runtime history. Features are built with [`campaign_route_features`](@ref)
from the ensemble shape (member count as the sample count, one satellite per sample, the
configuration's density-model family and mission length); the decision comes from
[`select_outer_route!`](@ref) and per-member success and amortized wall-clock feedback is
recorded via [`record_outer_route_feedback!`](@ref), so repeated ensembles converge to the
fastest allocation. History lives in [`campaign_outer_route_state`](@ref) unless an
isolated [`OuterRouteState`](@ref) is passed as `route_state`; `route_tuning` overrides
the [`OuterRouteTuning`](@ref). When the adaptive route is threaded, the runner also sets
`SPACEAGORA_INNER_THREAD_BUDGET` (unless already set) so inner and outer parallelism
split the thread pool instead of oversubscribing it. Adaptive members are isolated the
same way as the fixed-thread parallel path, whichever route the bandit selects.

The configuration must be uncoupled: guidance, navigation, and control effector tuples
must be empty, because effectors that coordinate several satellites cannot act across
ensemble members. If every configured effector acts on one satellite only, pass
`allow_gnc_effectors=true` to opt in. Remaining keyword arguments are forwarded to
[`run_simulation`](@ref) (for example `return_solution=true`).

```julia
result = run_constellation_ensemble(args; threads=8, return_solution=true)
solutions = [sample.value for sample in result.successful]
```
"""
function run_constellation_ensemble(
    args::SimulationConfiguration;
    threads::Union{Integer, Symbol}=1,
    fail_fast::Bool=false,
    allow_gnc_effectors::Bool=false,
    route_state::Union{Nothing, OuterRouteState}=nothing,
    route_tuning::Union{Nothing, OuterRouteTuning}=nothing,
    run_kwargs...
)
    spacecraft = args.dynamics_model.spacecraft
    isempty(spacecraft) && throw(ArgumentError(
        "run_constellation_ensemble requires at least one spacecraft in args.dynamics_model.spacecraft."
    ))
    _validate_ensemble_uncoupled(args, allow_gnc_effectors)

    member_configs = [
        _ensemble_member_configuration(args, sc, "sat_$(idx)_id_$(sc.id)")
        for (idx, sc) in enumerate(spacecraft)
    ]

    if threads isa Symbol
        threads === :auto || throw(ArgumentError(
            "run_constellation_ensemble threads must be a positive integer or :auto; got :$(threads)."
        ))
        # The bandit may route this campaign to threaded execution, so members
        # must be isolated exactly as in the fixed-thread parallel path below:
        # a per-member deepcopy inside the worker replaces run_simulation's own
        # isolate_state copy (cost-equivalent when the route is serial, required
        # when it is threaded).
        adaptive_kwargs = Dict{Symbol, Any}(pairs(run_kwargs))
        adaptive_kwargs[:isolate_state] = false
        run_member_isolated = idx -> SimulationEngine.run_simulation(
            deepcopy(member_configs[idx]); adaptive_kwargs...
        )
        # Each ensemble member propagates one satellite, so the routed workload
        # is single-sat per sample regardless of the constellation size.
        features = campaign_route_features(args; samples=length(member_configs), n_sats=1)
        state = route_state === nothing ? campaign_outer_route_state() : route_state
        tuning = route_tuning === nothing ? _campaign_route_tuning() : route_tuning
        return _run_campaign_adaptive(
            run_member_isolated,
            1:length(member_configs);
            fail_fast=fail_fast,
            features=features,
            state=state,
            tuning=tuning
        )
    end
    if !(route_state === nothing && route_tuning === nothing)
        throw(ArgumentError(
            "run_constellation_ensemble route_state/route_tuning are only consulted with threads=:auto."
        ))
    end
    spec = MonteCarloSpec(seeds=1:length(member_configs), threads=threads, fail_fast=fail_fast)
    worker_count = min(spec.threads, length(member_configs))

    if worker_count > 1
        # Member splits alias the parent's environment/GNC models, and solves mutate
        # model state (for example planet.L_PI in the planet-frame callback). Copy per
        # member inside the worker so concurrent members share nothing, and disable
        # run_simulation's own isolate_state copy so isolation is paid exactly once.
        parallel_kwargs = Dict{Symbol, Any}(pairs(run_kwargs))
        parallel_kwargs[:isolate_state] = false
        run_member = idx -> SimulationEngine.run_simulation(
            deepcopy(member_configs[idx]); parallel_kwargs...
        )
        return withenv("SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1") do
            run_monte_carlo(run_member, spec)
        end
    end

    run_member = idx -> SimulationEngine.run_simulation(member_configs[idx]; run_kwargs...)
    return run_monte_carlo(run_member, spec)
end
