module Stage1Design

using ..ConstellationUtils
using ..ConstellationSlots
using ..ConstellationPlotting
using ..Stage0Seeding
using ..Stage1Controllable
using LinearAlgebra
using Statistics
using Random
using DataFrames, CSV
using Arrow
using Dates
using ProgressMeter

"""
    run_constellation_design(config_dict::AbstractString) -> Dict{String,Any}

Stage 1 constellation design entry point. Loads a YAML config file and runs
constellation design optimization.

# Arguments
- `config_path::AbstractString`: Path to YAML configuration file

# Returns
- Dictionary containing optimization results and metadata
"""
function run_constellation_design(config_path::AbstractString)
    config_dict = ingest_yaml(config_path)
    return run_constellation_design(config_dict)
end

"""
    run_constellation_design(config_dict::AbstractDict) -> Dict{String,Any}

Stage 1 constellation design entry point. Runs LADS tube constellation design
optimization using the provided configuration dictionary.

# Arguments
- `config_dict::AbstractDict`: Configuration dictionary

# Returns
- Dictionary containing optimization results and metadata
"""
function run_constellation_design(config_dict::AbstractDict)
    constellation_log_init!(config_dict; context="stage1_design")
    
    try
        constellation_log("stage1", "Starting Stage 1 LADS tube constellation design")
        
        # Run Stage 0 if seeds not already generated
        if !haskey(config_dict, "stage0_h_fwd_exact_coeffs")
            constellation_log("stage1", "Running Stage 0 seeding first")
            opt_params = get(config_dict, "optimizer_params", Dict{String,Any}())
            stage0_method = String(get(opt_params, "stage0_mode", "fhsg"))
            Stage0Seeding.run_stage0_seeding(config_dict; method=stage0_method)
        end
        
        # Run Stage 1 LADS tube optimization
        constellation_log("stage1", "Running Stage 1 LADS tube optimization")
        opt_result = run_stage1_controllable_optimization(config_dict)
        
        # Store results in config_dict for downstream stages
        config_dict["stage1_results"] = opt_result
        
        constellation_log("stage1", "Stage 1 completed successfully")
        
        return opt_result
    catch err
        constellation_log_exception("stage1", err)
        rethrow(err)
    finally
        constellation_log_close!()
    end
end

export run_constellation_design

end # module Stage1Design
