# Monte Carlo campaign whose samples each write a viewer-ready results bundle,
# plus the ensemble manifest and page that merge them
# (docs/architecture/interactive_visualization_plan.md, phase 4).

using ..SimulationModel: SceneVisualization
using ..SimulationModel.SimConfig: _with_configuration
using ..SimulationModel.SceneVisualization: EnsembleSample, sample_results_directory, with_results_directory
using ..SimulationModel.SceneVisualization: default_sample_scalar, write_ensemble_manifest, export_ensemble_visualization
using ..SimulationEngine: run_simulation
using Arrow: Arrow
using DataFrames: DataFrame

"""
    run_monte_carlo_visualization(build_args, seeds, campaign_dir; scalar=default_sample_scalar, scalar_name="final spherical periapsis altitude (km)", threads=1, fail_fast=false, export_page=true, labels=nothing, nominal=nothing, page_kwargs...) -> (result, page)

Run one simulation per seed with `build_args(seed)::SimulationConfiguration`,
each writing its bundle and scene sidecar to `sample_results_directory(campaign_dir, index)`.
Checkpoint and resume files are isolated by the same sample tag, including
when the builder supplies a shared explicit checkpoint directory.
Then write the ensemble manifest and (unless `export_page=false`) build the
ensemble viewer page. `scalar(df)` reads one number per finished sample from
its results table for the page's color scale. Returns the
`MonteCarloResult` (each successful sample's `value` is its scalar) and the
page path or `nothing`. `page_kwargs` go to `export_ensemble_visualization`.
`nominal`, when given, is passed to `build_args` like a seed and run as one
more sample labeled "nominal"; the manifest records it and the page draws it
distinctly with the 3-sigma tube of the other samples around it.
"""
function run_monte_carlo_visualization(
    build_args,
    seeds,
    campaign_dir::AbstractString;
    scalar=default_sample_scalar,
    scalar_name::AbstractString="final spherical periapsis altitude (km)",
    threads::Union{Integer, Symbol}=1,
    fail_fast::Bool=false,
    export_page::Bool=true,
    labels::Union{Nothing, AbstractVector}=nothing,
    nominal=nothing,
    page_kwargs...
)
    campaign_dir = String(campaign_dir)
    mkpath(campaign_dir)
    seed_list = collect(seeds)
    labels === nothing || length(labels) == length(seed_list) || throw(ArgumentError("labels must have one entry per seed."))
    nominal_index = nominal === nothing ? nothing : length(seed_list) + 1
    nominal === nothing || push!(seed_list, nominal)
    indexed = collect(enumerate(seed_list))
    function sample((index, seed))
        base_args = build_args(seed)
        sample_dir = sample_results_directory(campaign_dir, index)
        # Share the existing ensemble checkpoint ownership policy. An explicit
        # checkpoint root must not let one seed overwrite or resume another.
        settings = _ensemble_member_settings(base_args.simulation_settings, basename(sample_dir))
        args = with_results_directory(
            _with_configuration(base_args; simulation_settings=settings),
            sample_dir; visualization=true,
        )
        run_simulation(args)
        df = DataFrame(Arrow.Table(joinpath(args.simulation_settings.results_directory, "simulation_results.feather")))
        return Float64(scalar(df))
    end
    result = run_monte_carlo(sample, indexed; threads=threads, fail_fast=fail_fast)
    samples = EnsembleSample[]
    for s in result.samples
        index, seed = s.seed
        value = s.success && s.value isa Real ? Float64(s.value) : NaN
        label = index == nominal_index ? "nominal" : (labels === nothing ? "sample $(index) (seed $(seed))" : String(labels[index]))
        push!(samples, EnsembleSample(Int(index), string(seed), s.success, value, label, sample_results_directory(campaign_dir, index)))
    end
    sort!(samples; by=s -> s.index)
    write_ensemble_manifest(campaign_dir, samples; scalar_name=scalar_name, nominal=nominal_index)
    page = export_page ? export_ensemble_visualization(campaign_dir; page_kwargs...) : nothing
    return result, page
end
