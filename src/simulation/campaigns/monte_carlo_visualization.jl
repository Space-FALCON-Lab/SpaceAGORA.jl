# Monte Carlo campaign whose samples each write a viewer-ready results bundle,
# plus the ensemble manifest and page that merge them
# (docs/architecture/interactive_visualization_plan.md, phase 4).

using ..SimulationModel: SceneVisualization
using ..SimulationModel.SceneVisualization: EnsembleSample, sample_results_directory, with_results_directory
using ..SimulationModel.SceneVisualization: default_sample_scalar, write_ensemble_manifest, export_ensemble_visualization
using ..SimulationEngine: run_simulation
using Arrow: Arrow
using DataFrames: DataFrame

"""
    run_monte_carlo_visualization(build_args, seeds, campaign_dir; scalar=default_sample_scalar, scalar_name="final periapsis altitude (km)", threads=1, fail_fast=false, export_page=true, labels=nothing, page_kwargs...) -> (result, page)

Run one simulation per seed with `build_args(seed)::SimulationConfiguration`,
each writing its bundle and scene sidecar to `sample_results_directory(campaign_dir, index)`,
then write the ensemble manifest and (unless `export_page=false`) build the
ensemble viewer page. `scalar(df)` reads one number per finished sample from
its results table for the page's colour scale. Returns the
`MonteCarloResult` (each successful sample's `value` is its scalar) and the
page path or `nothing`. `page_kwargs` go to `export_ensemble_visualization`.
"""
function run_monte_carlo_visualization(
    build_args,
    seeds,
    campaign_dir::AbstractString;
    scalar=default_sample_scalar,
    scalar_name::AbstractString="final periapsis altitude (km)",
    threads::Union{Integer, Symbol}=1,
    fail_fast::Bool=false,
    export_page::Bool=true,
    labels::Union{Nothing, AbstractVector}=nothing,
    page_kwargs...
)
    campaign_dir = String(campaign_dir)
    mkpath(campaign_dir)
    seed_list = collect(seeds)
    labels === nothing || length(labels) == length(seed_list) || throw(ArgumentError("labels must have one entry per seed."))
    indexed = collect(enumerate(seed_list))
    function sample((index, seed))
        args = with_results_directory(build_args(seed), sample_results_directory(campaign_dir, index); visualization=true)
        run_simulation(args)
        df = DataFrame(Arrow.Table(joinpath(args.simulation_settings.results_directory, "simulation_results.feather")))
        return Float64(scalar(df))
    end
    result = run_monte_carlo(sample, indexed; threads=threads, fail_fast=fail_fast)
    samples = EnsembleSample[]
    for s in result.samples
        index, seed = s.seed
        value = s.success && s.value isa Real ? Float64(s.value) : NaN
        label = labels === nothing ? "sample $(index) (seed $(seed))" : String(labels[index])
        push!(samples, EnsembleSample(Int(index), string(seed), s.success, value, label, sample_results_directory(campaign_dir, index)))
    end
    sort!(samples; by=s -> s.index)
    write_ensemble_manifest(campaign_dir, samples; scalar_name=scalar_name)
    page = export_page ? export_ensemble_visualization(campaign_dir; page_kwargs...) : nothing
    return result, page
end
