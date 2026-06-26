#!/usr/bin/env julia

"""
    rl_vs_stochastic_comparison.jl

Example script for comparing RL-based seeding with stochastic greedy seeding.
This script demonstrates the test harness for performance comparison.
"""

using SpaceAGORA.ConstellationDesign
using ArgParse
using Logging

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--labeled-csv"
            help = "Path to labeled debris CSV"
            arg_type = String
            default = "analysis_tools/cluster_output/labeled_debris.csv"
        "--rl-model"
            help = "Path to trained RL model"
            arg_type = String
            default = "data/constellation_rl_models/latest_model.jld2"
        "--output-dir"
            help = "Output directory for comparison results"
            arg_type = String
            default = "data/constellation_rl_comparison"
        "--clusters"
            help = "Cluster IDs to test (comma-separated)"
            arg_type = String
            default = "1a,2a,3a"
        "--combinations"
            help = "Test cluster combinations (individual, pairs, all)"
            arg_type = String
            default = "individual"
        "--verbose", "-v"
            help = "Enable verbose logging"
            action = :store_true
    end

    return parse_args(s)
end

function main()
    args = parse_commandline()
    
    # Setup logging
    if args["verbose"]
        logger = ConsoleLogger(stdout, Logging.Debug)
        global_logger(logger)
    else
        logger = ConsoleLogger(stdout, Logging.Info)
        global_logger(logger)
    end
    
    @info "Starting constellation RL vs Stochastic Greedy comparison"
    @info "Labeled CSV: $(args["labeled-csv"])"
    @info "Constellation RL Model: $(args["rl-model"])"
    @info "Output directory: $(args["output-dir"])"
    
    # Parse cluster IDs
    cluster_ids = split(args["clusters"], ",")
    
    # Generate cluster combinations based on mode
    if args["combinations"] == "individual"
        cluster_sets = [[id] for id in cluster_ids]
    elseif args["combinations"] == "pairs"
        cluster_sets = []
        for i in 1:length(cluster_ids)
            for j in (i+1):length(cluster_ids)
                push!(cluster_sets, [cluster_ids[i], cluster_ids[j]])
            end
        end
    elseif args["combinations"] == "all"
        cluster_sets = generate_cluster_combinations(cluster_ids)
    else
        @error "Unknown combinations mode: $(args["combinations"])"
        return
    end
    
    @info "Generated $(length(cluster_sets)) cluster combinations"
    
    # Build comparison scenarios
    @info "Building comparison scenarios..."
    scenarios = build_comparison_config(
        args["labeled-csv"],
        cluster_sets,
        args["output-dir"]
    )
    
    @info "Built $(length(scenarios)) scenarios"
    
    # Run comparisons
    @info "Running comparisons..."
    results = run_comparison_batch(scenarios, args["rl-model"])
    
    @info "Completed $(length(results)) comparisons"
    
    # Generate report
    report_path = joinpath(args["output-dir"], "comparison_results.csv")
    summary = generate_comparison_report(results, report_path)
    
    @info "Comparison complete"
    @info "  Results saved to: $report_path"
    @info "  Summary saved to: $(replace(report_path, ".csv" => "_summary.yaml"))"
end

# Run main function
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
