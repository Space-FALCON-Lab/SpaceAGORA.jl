#!/usr/bin/env julia

# debris_controllable_sim.jl
#
# Mission entry point for the CAPO constellation design pipeline in SpaceAGORA.
# This script replicates the debris_controllable_sim pipeline from CAPOConstellations,
# using the default formulation (polyhedral set mode with constructive zonotope certificates).
#
# Usage:
#   julia --project src/mission/constellation_design/debris_controllable_sim.jl config=<path>
#   julia --project src/mission/constellation_design/debris_controllable_sim.jl config=<path> mode=stochastic_greedy

"""
    parse_debris_args(args::Vector{String}) -> Dict{String,Any}

Parse command-line arguments for debris_controllable_sim.
"""
function parse_debris_args(args::Vector{String})
    parsed = Dict{String,Any}()
    
    for arg in args
        if contains(arg, "=")
            key, value = split(arg, "="; limit=2)
            parsed[strip(key)] = strip(value)
        end
    end
    
    return parsed
end

"""
    run_debris_controllable_sim(config_path::AbstractString; kwargs...) -> Dict{String,Any}

Run the full debris controllable simulation pipeline with all three stages.

# Arguments
- `config_path::AbstractString`: Path to YAML configuration file
- `kwargs`: Additional configuration overrides

# Returns
- Dictionary containing results from all stages
"""
function run_debris_controllable_sim(config_path::AbstractString; kwargs...)
    # Load configuration
    config_dict = ingest_yaml(config_path)
    
    # Apply command-line overrides
    for (key, value) in kwargs
        config_dict[key] = value
    end
    
    return run_debris_controllable_sim(config_dict)
end

"""
    run_debris_controllable_sim(config_dict::AbstractDict) -> Dict{String,Any}

Run the full debris controllable simulation pipeline with all three stages.

This is the main pipeline integration point for the CAPOConstellation port.
Enforces Stage 0 → Stage 1 → Stage 2 ordering with LADS tube formulation.

# Arguments
- `config_dict::AbstractDict`: Configuration dictionary

# Returns
- Dictionary containing results from all stages
"""
function run_debris_controllable_sim(config_dict::AbstractDict)
    constellation_log_init!(config_dict; context="debris_controllable_sim")
    
    try
        constellation_log("pipeline", "Starting debris_controllable_sim pipeline")
        
        results = Dict{String,Any}()
        
        # Stage 0: FHSG Seeding (always runs unless cached)
        opt_params = get(config_dict, "optimizer_params", Dict{String,Any}())
        stage0_method = String(get(opt_params, "stage0_mode", "fhsg"))
        constellation_log("pipeline", "Running Stage 0: Seeding (method=$stage0_method)")
        results["stage0"] = run_stage0_seeding(config_dict; method=stage0_method)
        
        # Stage 1: LADS Tube Certificate (always runs after Stage 0)
        constellation_log("pipeline", "Running Stage 1: LADS Tube Certificate")
        results["stage1"] = run_constellation_design(config_dict)
        
        # Stage 2: Sequential OCP Verification (always runs after Stage 1)
        constellation_log("pipeline", "Running Stage 2: Sequential OCP Verification")
        results["stage2"] = run_stage2_verification(config_dict)
        
        constellation_log("pipeline", "Pipeline completed successfully")
        
        return results
    catch err
        constellation_log_exception("pipeline", err)
        rethrow(err)
    finally
        constellation_log_close!()
    end
end

"""
    main(args::Vector{String}=ARGS)

Main entry point for debris_controllable_sim.
"""
function main(args::Vector{String}=ARGS)
    # Parse arguments
    parsed_args = parse_debris_args(args)
    
    # Get config path
    config_path = get(parsed_args, "config", nothing)
    if config_path === nothing
        error("config=<path> argument is required")
    end
    
    if !isfile(config_path)
        error("Config file not found: $config_path")
    end
    
    # Run pipeline
    results = run_debris_controllable_sim(config_path; 
        mode=get(parsed_args, "mode", "full"),
        dry_run=get(parsed_args, "dry_run", "false") == "true",
    )
    
    println("Pipeline completed successfully")
    println("Results: $(keys(results))")
    
    return results
endbspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
