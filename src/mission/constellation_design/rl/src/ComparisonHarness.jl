module ComparisonHarness

using ..CapoIntegration
using ..CostFunction
using ..Scenarios
using ..ClusterCombinations
using Statistics
using DataFrames
using CSV
using YAML

"""
    ComparisonResult

Structure representing the result of comparing two seeding methods.
"""
struct ComparisonResult
    scenario_id::String
    cluster_set::Vector{String}
    n_clients::Int
    
    # Stochastic greedy results
    sg_cost::Float64
    sg_n_sats::Int
    sg_feasible::Bool
    sg_time::Float64
    
    # RL results
    rl_cost::Float64
    rl_n_sats::Int
    rl_feasible::Bool
    rl_time::Float64
    
    # Comparison metrics
    cost_improvement::Float64
    sat_reduction::Int
    time_speedup::Float64
    rl_better::Bool
end

"""
    run_comparison(config_dict::Dict{String,Any}, rl_model_path::String) -> ComparisonResult

Run a single comparison between stochastic greedy and RL seeding.
"""
function run_comparison(config_dict::Dict{String,Any}, rl_model_path::String)
    scenario_id = get(config_dict, "scenario_id", "unknown")
    cluster_set = get(config_dict, "cluster_set", ["unknown"])
    
    # Run stochastic greedy
    sg_start = time()
    sg_result = run_stochastic_greedy_seeding(config_dict)
    sg_time = time() - sg_start
    
    # Evaluate stochastic greedy
    sg_cost_result = compute_stage0_cost(config_dict, sg_result["constellation_orbitals"])
    
    # Run RL seeding
    rl_start = time()
    config_dict["optimizer_params"]["rl_config"]["rl_model_path"] = rl_model_path
    rl_result = run_rl_stage0_seeding(config_dict)
    rl_time = time() - rl_start
    
    # Evaluate RL
    rl_cost_result = compute_stage0_cost(config_dict, rl_result["constellation_orbitals"])
    
    # Compute comparison metrics
    cost_improvement = sg_cost_result["total_cost"] - rl_cost_result["total_cost"]
    sat_reduction = sg_result["n_sats"] - rl_result["n_sats"]
    time_speedup = sg_time / rl_time
    rl_better = rl_cost_result["total_cost"] < sg_cost_result["total_cost"]
    
    return ComparisonResult(
        scenario_id,
        cluster_set,
        sg_cost_result["n_clients"],
        
        sg_cost_result["total_cost"],
        sg_result["n_sats"],
        sg_cost_result["feasible"],
        sg_time,
        
        rl_cost_result["total_cost"],
        rl_result["n_sats"],
        rl_cost_result["feasible"],
        rl_time,
        
        cost_improvement,
        sat_reduction,
        time_speedup,
        rl_better,
    )
end

"""
    run_comparison_batch(scenarios::Vector{Dict{String,Any}}, rl_model_path::String) -> Vector{ComparisonResult}

Run comparisons for multiple scenarios.
"""
function run_comparison_batch(scenarios::Vector{Dict{String,Any}}, rl_model_path::String)
    results = Vector{ComparisonResult}()
    
    for scenario in scenarios
        @info "Running comparison for scenario: $(scenario["scenario_id"])"
        try
            result = run_comparison(scenario, rl_model_path)
            push!(results, result)
        catch err
            @warn "Scenario $(scenario["scenario_id"]) failed: $err"
        end
    end
    
    return results
end

"""
    generate_comparison_report(results::Vector{ComparisonResult}, output_path::String)

Generate a comprehensive comparison report.
"""
function generate_comparison_report(results::Vector{ComparisonResult}, output_path::String)
    mkpath(dirname(output_path))
    
    # Compute aggregate statistics
    n_results = length(results)
    rl_better_count = count(r -> r.rl_better, results)
    rl_better_pct = 100.0 * rl_better_count / n_results
    
    avg_cost_improvement = mean(r.cost_improvement for r in results)
    avg_sat_reduction = mean(r.sat_reduction for r in results)
    avg_time_speedup = mean(r.time_speedup for r in results)
    
    sg_feasible_count = count(r -> r.sg_feasible, results)
    rl_feasible_count = count(r -> r.rl_feasible, results)
    
    # Create DataFrame for detailed results
    df = DataFrame(
        scenario_id = [r.scenario_id for r in results],
        cluster_set = [join(r.cluster_set, "+") for r in results],
        n_clients = [r.n_clients for r in results],
        
        sg_cost = [r.sg_cost for r in results],
        sg_n_sats = [r.sg_n_sats for r in results],
        sg_feasible = [r.sg_feasible for r in results],
        sg_time = [r.sg_time for r in results],
        
        rl_cost = [r.rl_cost for r in results],
        rl_n_sats = [r.rl_n_sats for r in results],
        rl_feasible = [r.rl_feasible for r in results],
        rl_time = [r.rl_time for r in results],
        
        cost_improvement = [r.cost_improvement for r in results],
        sat_reduction = [r.sat_reduction for r in results],
        time_speedup = [r.time_speedup for r in results],
        rl_better = [r.rl_better for r in results],
    )
    
    # Save detailed results
    CSV.write(output_path, df)
    
    # Create summary
    summary = Dict{String,Any}(
        "total_scenarios" => n_results,
        "rl_better_count" => rl_better_count,
        "rl_better_percentage" => rl_better_pct,
        "avg_cost_improvement" => avg_cost_improvement,
        "avg_sat_reduction" => avg_sat_reduction,
        "avg_time_speedup" => avg_time_speedup,
        "sg_feasible_count" => sg_feasible_count,
        "rl_feasible_count" => rl_feasible_count,
        "detailed_results_path" => output_path,
    )
    
    # Save summary
    summary_path = replace(output_path, ".csv" => "_summary.yaml")
    YAML.write_file(summary_path, summary)
    
    @info "Comparison report generated"
    @info "  Total scenarios: $n_results"
    @info "  RL better: $rl_better_count ($rl_better_pct%)"
    @info "  Avg cost improvement: $avg_cost_improvement"
    @info "  Avg sat reduction: $avg_sat_reduction"
    @info "  Avg time speedup: $avg_time_speedup"
    
    return summary
end

"""
    build_comparison_config(labeled_csv_path::String, cluster_sets::Vector{Vector{String}}, output_dir::String) -> Vector{Dict{String,Any}}

Build configuration dictionaries for comparison scenarios.
"""
function build_comparison_config(labeled_csv_path::String, cluster_sets::Vector{Vector{String}}, output_dir::String)
    # Load labeled debris
    df = load_labeled_debris_csv(labeled_csv_path)
    
    # Build scenarios
    scenarios = Vector{Dict{String,Any}}()
    
    for (idx, cluster_set) in enumerate(cluster_sets)
        # Build cluster scenario
        cluster_scenario = build_cluster_scenario(df, cluster_set, output_dir)
        
        if cluster_scenario === nothing
            continue
        end
        
        # Build config dict
        config_dict = Dict{String,Any}(
            "scenario_id" => "comparison_$(lpad(idx, 3, '0'))",
            "cluster_set" => cluster_set,
            "client_orbitals" => cluster_scenario["client_orbitals"],
            "orbital_bounds" => cluster_scenario["orbital_bounds"],
            "optimizer_params" => Dict{String,Any}(
                "rl_config" => Dict{String,Any}(
                    "max_sats" => 64,
                    "greedy_variant" => "pure",
                    "greedy_top_k" => 5,
                    "greedy_epsilon" => 0.1,
                    "greedy_restarts" => 1,
                    "unsafe_weight" => 100.0,
                    "safe_weight" => 10.0,
                    "pred_weight" => 5.0,
                    "feasibility_threshold" => 1e-6,
                ),
            ),
        )
        
        push!(scenarios, config_dict)
    end
    
    return scenarios
end

export ComparisonResult, run_comparison, run_comparison_batch, generate_comparison_report, build_comparison_config

end # module ComparisonHarness
