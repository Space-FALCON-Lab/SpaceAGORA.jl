using ..SimulationModel: SimulationConfiguration, SimulationSettings, DynamicsModel, SpacecraftModel
import ..SimulationEngine

function _ensemble_member_settings(settings::SimulationSettings, member_tag::String)::SimulationSettings
    needs_results_split = settings.results
    needs_checkpoint_split = settings.checkpoint_enabled && !isempty(strip(settings.checkpoint_directory))
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
                               allow_gnc_effectors=false, run_kwargs...) -> MonteCarloResult

Propagate each spacecraft of a multi-satellite [`SimulationConfiguration`](@ref) as an
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
result saving is enabled, each member writes under
`results_directory/sat_<index>_id_<spacecraft id>` so members do not clobber each other.

`threads` caps worker tasks exactly as in [`run_monte_carlo`](@ref); Julia must be started
with at least that many threads. While more than one worker is active the runner sets
`SPACEAGORA_OUTER_PARALLEL_ACTIVE=1` so inner thread policies yield to the outer split.

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
    threads::Integer=1,
    fail_fast::Bool=false,
    allow_gnc_effectors::Bool=false,
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

    run_member = idx -> SimulationEngine.run_simulation(member_configs[idx]; run_kwargs...)
    spec = MonteCarloSpec(seeds=1:length(member_configs), threads=threads, fail_fast=fail_fast)

    worker_count = min(spec.threads, length(member_configs))
    if worker_count > 1
        return withenv("SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1") do
            run_monte_carlo(run_member, spec)
        end
    end
    return run_monte_carlo(run_member, spec)
end
