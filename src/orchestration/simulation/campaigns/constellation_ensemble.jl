import ..SimulationModel
using ..SimulationModel: SimulationConfiguration, SimulationSettings, DynamicsModel, SpacecraftModel
import ..SimulationEngine
using ..SimulationModel.SimConfig: _with_configuration, _with_fields
using ..SimulationModel.DynamicEffectors: PlumeSurfaceInteraction

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
    return _with_fields(settings;
        results_directory=needs_results_split ?
            joinpath(settings.results_directory, member_tag) : settings.results_directory,
        checkpoint_directory=needs_checkpoint_split ?
            joinpath(settings.checkpoint_directory, member_tag) : settings.checkpoint_directory,
    )
end

# The folder of ensemble member `index` of `count`: its index, zero-padded to the
# digits of `count` so the folders sort in member order, and its spacecraft id.
_ensemble_member_tag(index::Integer, count::Integer, id::Integer) =
    "sat_$(lpad(index, ndigits(count), '0'))_id_$(id)"

function _ensemble_member_configuration(
    args::SimulationConfiguration,
    spacecraft::SpacecraftModel,
    member_tag::String;
    effectors::Tuple=args.dynamics_model.dynamic_effectors,
)::SimulationConfiguration
    return _with_configuration(args;
        simulation_settings=_ensemble_member_settings(args.simulation_settings, member_tag),
        dynamics_model=DynamicsModel([spacecraft], effectors),
    )
end

# The dynamic effectors of ensemble member `index`, from a configuration whose
# effector scopes are resolved (SimulationEngine._resolve_effector_scopes): every
# effector without a table, otherwise the effectors whose scope holds the member's
# position in the constellation. The member runs them unscoped.
function _ensemble_member_effectors(args::SimulationConfiguration, index::Integer)::Tuple
    scopes = args.dynamics_model.effector_scopes
    effectors = args.dynamics_model.dynamic_effectors
    scopes === nothing && return effectors
    applies = SimulationModel.SpacecraftModels._effector_applies
    return Tuple(e for (j, e) in enumerate(effectors) if applies(scopes, j, Int(index)))
end

# The configuration of ensemble member `index` of `count`. Checkpoints of an
# ensemble of ten or more members written before the folders were zero-padded
# (SpaceAGORA 0.4 and earlier) are under the unpadded name, `sat_1_id_1`: on
# resume, a member whose padded folder has no checkpoint and whose unpadded one
# has resumes there, and goes on writing there.
function _ensemble_member_configuration(args::SimulationConfiguration, spacecraft::SpacecraftModel,
                                        index::Integer, count::Integer)::SimulationConfiguration
    tag = _ensemble_member_tag(index, count, spacecraft.id)
    effectors = _ensemble_member_effectors(args, index)
    config = _ensemble_member_configuration(args, spacecraft, tag; effectors=effectors)
    unpadded = "sat_$(index)_id_$(spacecraft.id)"
    (args.simulation_settings.resume_from_checkpoint && unpadded != tag) || return config
    isfile(SimulationEngine._checkpoint_paths(config).data) && return config
    earlier = _ensemble_member_configuration(args, spacecraft, unpadded; effectors=effectors)
    return isfile(SimulationEngine._checkpoint_paths(earlier).data) ? earlier : config
end

# The run positions a dynamic effector is bound to, by the fields the engine
# reads as spacecraft indices: `spacecraft_index` (the robot-arm reaction and
# reaction-wheel effectors) and `chaser_index`/`target_index`.
function _ensemble_effector_bindings(effector)
    bindings = Pair{Symbol, Int}[]
    index = SimulationEngine._effector_spacecraft_index(effector)
    index isa Integer && push!(bindings, :spacecraft_index => Int(index))
    for name in SimulationEngine._CROSS_SATELLITE_INDEX_FIELDS
        hasproperty(effector, name) || continue
        index = getproperty(effector, name)
        index isa Integer && push!(bindings, name => Int(index))
    end
    return bindings
end

# Every ensemble member runs its spacecraft alone, at run position 1, with the
# parent's whole dynamic effector tuple. An effector bound to a run position
# therefore keeps its meaning only in a one-member ensemble bound to position 1:
# bound to k > 1 it never acts in member k, and bound to 1 it acts on every
# member. A plume carries thrust and state for the parent's spacecraft count,
# while each member's plume check needs a count of one. Both are refused rather
# than remapped, so no effector silently moves to another spacecraft.
function _validate_ensemble_dynamic_effectors(args::SimulationConfiguration)
    spacecraft = args.dynamics_model.spacecraft
    n_members = length(spacecraft)
    # run_constellation_ensemble refuses an empty constellation before this check.
    n_members == 0 && return nothing
    member(k) = "member $(k) (spacecraft id $(spacecraft[k].id))"
    for (position, effector) in enumerate(args.dynamics_model.dynamic_effectors)
        problem = nothing
        for (field, index) in _ensemble_effector_bindings(effector)
            problem = if !(1 <= index <= n_members)
                "its $(field)=$(index) names no member of this $(n_members)-member ensemble"
            elseif index != 1
                "its $(field)=$(index) binds it to $(member(index)), whose spacecraft runs alone at " *
                "position 1, so it would never act there"
            elseif n_members > 1
                "its $(field)=1 binds it to $(member(1)), but every member's spacecraft runs alone at " *
                "position 1, so it would also act on $(member(2))"
            else
                nothing
            end
            problem === nothing || break
        end
        if problem === nothing && effector isa PlumeSurfaceInteraction.PlumeSurfaceInteractionModel
            thrusts = PlumeSurfaceInteraction._control_spacecraft_count(effector.control)
            problem = if thrusts != 1
                "it carries thrust and state for $(thrusts) spacecraft, but each member runs one, so " *
                "its plume check would fail in $(member(1))"
            elseif n_members > 1
                "it carries thrust and state for one spacecraft, but every member would run it, so it " *
                "would also act on $(member(2))"
            else
                nothing
            end
        end
        problem === nothing || throw(ArgumentError(
            "run_constellation_ensemble cannot split dynamic_effectors[$(position)] " *
            "($(nameof(typeof(effector)))) across its members: $(problem). Each member gets the whole " *
            "dynamic effector tuple, and spacecraft indices are not remapped. Use run_simulation for " *
            "the whole constellation, or one configuration per spacecraft with its own effectors."
        ))
    end
    return nothing
end

# With allow_gnc_effectors, every member receives the configuration-level guidance,
# navigation and control models unchanged and runs its spacecraft alone at run
# position 1. A model whose scope names run positions (a BoundSpacecraft, a
# CoupledSpacecraft or a SpacecraftSelection) would then act on every member, or
# fail in some, so it is refused, pointing to declaring it on its spacecraft. An
# ensemble of one runs its spacecraft at position 1 as the whole run would, so there
# a scope that names that position alone is kept, as for the dynamic effectors.
function _validate_ensemble_gnc_scopes(args::SimulationConfiguration)
    GNC = SimulationEngine.GNCRuntime
    spacecraft = args.dynamics_model.spacecraft
    groups = ((:guidance, args.guidance_model.guidance_effectors),
        (:navigation, args.navigation_model.navigation_effectors), (:control, args.control_model.control_effectors))
    for (group, models) in groups, (j, model) in enumerate(models)
        scope = GNC.spacecraft_scope(model)
        scope isa GNC.AllSpacecraft && continue
        length(spacecraft) == 1 && _ensemble_scope_positions(scope, spacecraft) == [1] && continue
        throw(ArgumentError(
            "run_constellation_ensemble: the configuration's $(SimulationEngine._gnc_scope_owner(group, j, model)) " *
            "declares spacecraft_scope $(SimulationModel.SimulationCallbacks._scope_text(scope)), but each " *
            "ensemble member runs its spacecraft alone, at run " *
            "position 1, with the configuration-level models unchanged, so it would act on every member or fail in " *
            "some. Declare it on the spacecraft it acts on (the SpacecraftModel's $(group) keyword), where each " *
            "member binds it to its own spacecraft, or give it AllSpacecraft() to act on every member."))
    end
    return nothing
end

# The run positions a GNC scope names among `spacecraft`, or `nothing` when it names
# one the run does not have.
_ensemble_scope_positions(scope::SimulationEngine.GNCRuntime.CoupledSpacecraft, spacecraft) =
    sort!(unique!([scope.primary; collect(scope.partners)]))
function _ensemble_scope_positions(scope, spacecraft)
    try
        return SimulationModel.resolve_spacecraft_scope(scope, spacecraft)
    catch err
        err isa ArgumentError || rethrow()
        return nothing
    end
end

function _validate_ensemble_uncoupled(args::SimulationConfiguration, allow_gnc_effectors::Bool)
    args.interlink_model === nothing || throw(ArgumentError(
        "run_constellation_ensemble cannot split an interlinked constellation into independent spacecraft. Use run_simulation instead."))
    _validate_ensemble_dynamic_effectors(args)
    allow_gnc_effectors && return _validate_ensemble_gnc_scopes(args)
    coupled_surfaces = String[]
    isempty(args.guidance_model.guidance_effectors) || push!(coupled_surfaces, "guidance_effectors")
    isempty(args.navigation_model.navigation_effectors) || push!(coupled_surfaces, "navigation_effectors")
    isempty(args.control_model.control_effectors) || push!(coupled_surfaces, "control_effectors")
    isempty(coupled_surfaces) && return nothing
    throw(ArgumentError(
        "run_constellation_ensemble requires an uncoupled constellation, but the configuration has " *
        "non-empty $(join(coupled_surfaces, ", ")). GNC effectors may couple satellites (for example " *
        "RPO planners or coordinated maneuvers), and each ensemble member propagates alone. If every " *
        "configuration-level effector acts on each satellite on its own (spacecraft_scope AllSpacecraft()), " *
        "opt in with allow_gnc_effectors=true. Declare GNC for one satellite on its SpacecraftModel instead, " *
        "which needs no flag; otherwise use the monolithic run_simulation path."
    ))
end

"""
    run_constellation_ensemble(args::SimulationConfiguration; threads=1, fail_fast=false,
                               allow_gnc_effectors=false, merge_results=false,
                               run_kwargs...) -> MonteCarloResult or ConstellationEnsembleResult

Propagate each spacecraft of a multi-satellite `SimulationConfiguration` as an
independent single-satellite simulation, parallelized across satellites with the
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
clobber each other's files. The index is zero-padded to the digits of the member count
(`sat_01` to `sat_12` for twelve members), so the folders sort in member order. On resume,
a member with a checkpoint only under its unpadded folder (written by SpaceAGORA 0.4 or
earlier, `sat_1_id_1`) resumes from that folder and keeps writing to it.

# Parallel execution

When `args.solver_config` is `SolverConfig(parallel=true)`, the ensemble lets the
predictive planner choose serial, threaded or process-worker execution for its
members, under the flag's settings for the duration of the call:

```julia
args = SimulationConfiguration(...; solver_config=SolverConfig(parallel=true))
result = run_constellation_ensemble(args; return_solution=true)
solutions = [sample.value for sample in result.successful]
```

The flag chooses the worker count itself, so combining it with an integer `threads`
throws an `ArgumentError`. Without the flag the ensemble is serial unless `threads`
gives a fixed number of worker tasks, exactly as in [`run_monte_carlo`](@ref); Julia
must be started with at least that many threads. While more than one worker is active
the runner sets `SPACEAGORA_OUTER_PARALLEL_ACTIVE=1` so inner thread policies yield to
the outer split, and an ensemble started inside another campaign's worker yields to
the enclosing split.

Isolation: when more than one worker is active, each member solves a `deepcopy` of its
configuration taken inside its worker task, because model state such as
`environment_model.planet.L_PI` is mutated during a solve and must never be shared
between concurrent members. The copy replaces (not adds to) `run_simulation`'s own
`isolate_state` copy, so the total cost matches a plain isolated run; an explicit
`isolate_state=false` is ignored under parallel execution. With `threads=1` the
keyword is forwarded to [`run_simulation`](@ref) unchanged.

`threads=:auto`, `route_state` and `route_tuning` drive the internal routing machinery
the benchmark harness measures (see `SpaceAGORA.ParallelProfiles`); they are not the
supported way to parallelize an ensemble.

The configuration must be uncoupled: `interlink_model` must be `nothing`, and configuration-level guidance, navigation, and control
effector tuples must be empty, because effectors that coordinate several satellites cannot
act across ensemble members. If every configuration-level effector acts on each satellite on
its own, with `spacecraft_scope` `AllSpacecraft()`, pass `allow_gnc_effectors=true` to opt in;
the effectors are passed to every member unchanged, so spacecraft indices are not remapped,
and a configuration-level model whose `spacecraft_scope` names run positions
(`BoundSpacecraft`, `CoupledSpacecraft` or a `SpacecraftSelection`) is refused with an
`ArgumentError` (since SpaceAGORA 0.7.0): it would act on every member, or fail in some,
unless the ensemble has one member and the scope names position 1 alone. A model that
picks its spacecraft by an index field of its own without declaring `spacecraft_scope`
is not detected. GNC declared on a `SpacecraftModel` (its `guidance`,
`navigation` and `control` keywords) always travels with that spacecraft, is bound to index 1
in its member run, and needs no flag.

Every member runs its spacecraft alone at run position 1 with the whole dynamic effector
tuple, so a dynamic effector bound to a spacecraft by run position cannot be split:
one with a `spacecraft_index`, `chaser_index` or `target_index` (such as
`RobotArmReactionEffector` and `ReactionWheelMomentumModel`), and the per-spacecraft
`PlumeSurfaceInteractionModel`. The ensemble refuses such an effector with an
`ArgumentError` naming it and the member, with or without `allow_gnc_effectors`, unless
the ensemble has one member and the effector is bound to position 1.

# Merged results

With `merge_results=true` the call returns a [`ConstellationEnsembleResult`](@ref):
the `MonteCarloResult` as `campaign`, plus a `members` table with every member's
saved rows and its identity (`spacecraft`, `spacecraft_id`, `plane`, `slot`, in
the long layout's columns, see [`long_table`](@ref)) and a `constellation` table
with one row per member: its identity, folder, status and last saved row. With
`results=true` it also writes `constellation_members.feather` (and `.csv` with
`save_csv`) and `constellation_summary.csv` in `results_directory`, beside the
member folders. Each member's own files, solve and sample `value` are what they
are without the keyword. The merged tables hold every saved row in memory;
`SPACEAGORA_RESULTS_THIN_STRIDE` and `SPACEAGORA_RESULTS_FINAL_ONLY` thin each
member's rows. `merge_results=true` with `return_solver_metadata=true` is
refused before any member runs, and so is a run-wide save field whose column
prefix is a per-spacecraft field's.

Save fields passed with `save_fields` run in each member, which has one
spacecraft: build them from a member's configuration, not from the whole
constellation's, whose default fields are sized for every spacecraft.

```julia
res = run_constellation_ensemble(args; merge_results=true)
plane_two = res.members[coalesce.(res.members.plane .== 2, false), :]
res.constellation[:, [:spacecraft, :spacecraft_id, :success, :time]]
```

Remaining keyword arguments are forwarded to [`run_simulation`](@ref) (for example
`return_solution=true`).
"""
function run_constellation_ensemble(
    args::SimulationConfiguration;
    threads::Union{Nothing, Integer, Symbol}=nothing,
    fail_fast::Bool=false,
    allow_gnc_effectors::Bool=false,
    merge_results::Bool=false,
    route_state::Union{Nothing, OuterRouteState}=nothing,
    route_tuning::Union{Nothing, OuterRouteTuning}=nothing,
    run_kwargs...
)
    isempty(args.external_propagators) || throw(ArgumentError(
        "run_constellation_ensemble does not support externally propagated spacecraft (args.external_propagators): " *
        "the owned spacecraft are coupled bodies of one scene and cannot be split into independent members."))
    spacecraft = args.dynamics_model.spacecraft
    isempty(spacecraft) && throw(ArgumentError(
        "run_constellation_ensemble requires at least one spacecraft in args.dynamics_model.spacecraft."
    ))
    # Scoped effectors are resolved against the whole constellation first, so the
    # checks below see the effectors themselves; each member then keeps the
    # effectors whose scope holds its position, unscoped.
    args = SimulationEngine._resolve_effector_scopes(args)
    _validate_ensemble_uncoupled(args, allow_gnc_effectors)
    # Refused here rather than by every member.
    args.simulation_settings.results_layout in (:wide, :long) || throw(ArgumentError(
        "SimulationSettings.results_layout must be :wide or :long; got :$(args.simulation_settings.results_layout)."))
    if merge_results
        # run_simulation refuses this pair with return_results, which every
        # member of a merged ensemble runs with.
        get(run_kwargs, :return_solver_metadata, false) === true && throw(ArgumentError(
            "run_constellation_ensemble with merge_results=true cannot return solver metadata: each member runs " *
            "with return_results=true, which run_simulation refuses with return_solver_metadata=true."))
        _validate_merged_save_fields(run_kwargs)
    end
    SimulationEngine._warn_repeated_spacecraft_ids(args)

    member_configs = [
        _ensemble_member_configuration(args, sc, idx, length(spacecraft))
        for (idx, sc) in enumerate(spacecraft)
    ]

    # SolverConfig(parallel=true) on the members: the whole ensemble runs under
    # the flag's scope with the planner choosing the route (threads=:auto).
    parallel = _campaign_parallel_flag(member_configs)
    threads_resolved = _campaign_threads(threads, parallel, "run_constellation_ensemble")
    run_members() = _run_constellation_members(args, member_configs, threads_resolved, fail_fast,
                                               route_state, route_tuning, run_kwargs, merge_results)
    campaign = parallel ? SimulationEngine._with_parallel_flag(run_members, true) : run_members()
    merge_results || return campaign
    return _constellation_ensemble_result(args, member_configs, campaign)
end

function _run_constellation_members(args::SimulationConfiguration, member_configs, threads, fail_fast::Bool,
                                    route_state, route_tuning, run_kwargs, merge_results::Bool=false)
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
        run_member_isolated = idx -> _constellation_member_run(
            deepcopy(member_configs[idx]), adaptive_kwargs, merge_results
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
        run_member = idx -> _constellation_member_run(
            deepcopy(member_configs[idx]), parallel_kwargs, merge_results
        )
        return withenv("SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1") do
            run_monte_carlo(run_member, spec)
        end
    end

    run_member = idx -> _constellation_member_run(member_configs[idx], run_kwargs, merge_results)
    return run_monte_carlo(run_member, spec)
end
