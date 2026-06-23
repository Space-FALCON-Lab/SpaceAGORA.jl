module ConstellationDesign

include("constellation_utils.jl")
include("constellation_slots.jl")
include("constellation_plotting.jl")
include("stage0_seeding.jl")
include("stage1_design.jl")
include("stage2_verification.jl")

# debris_controllable_sim functions are defined directly in this module

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
"""
function run_debris_controllable_sim(config_path::AbstractString; kwargs...)
    config_dict = ingest_yaml(config_path)
    for (key, value) in kwargs
        config_dict[key] = value
    end
    return run_debris_controllable_sim(config_dict)
end

"""
    run_debris_controllable_sim(config_dict::AbstractDict) -> Dict{String,Any}

Run the full debris controllable simulation pipeline with all three stages.
"""
function run_debris_controllable_sim(config_dict::AbstractDict)
    constellation_log_init!(config_dict; context="debris_controllable_sim")
    
    try
        constellation_log("pipeline", "Starting debris_controllable_sim pipeline")
        
        opt_params = get(config_dict, "optimizer_params", Dict{String,Any}())
        mode = lowercase(String(get(opt_params, "mode", "full")))
        
        results = Dict{String,Any}()
        
        if mode in ("full", "stochastic_greedy", "stage0")
            constellation_log("pipeline", "Running Stage 0: Stochastic Seeding")
            results["stage0"] = run_stage0_seeding(config_dict)
        end
        
        if mode in ("full", "heuristic", "stage1")
            constellation_log("pipeline", "Running Stage 1: Constellation Design")
            results["stage1"] = run_constellation_design(config_dict)
        end
        
        if mode in ("full", "stage2")
            constellation_log("pipeline", "Running Stage 2: Control Verification")
            results["stage2"] = run_stage2_verification(config_dict)
        end
        
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
    run_capo_pipeline(config_dict::AbstractDict) -> Dict{String,Any}

Run the full CAPO pipeline (alias for run_debris_controllable_sim).
"""
function run_capo_pipeline(config_dict::AbstractDict)
    return run_debris_controllable_sim(config_dict)
end

function run_capo_pipeline(config_path::AbstractString)
    return run_debris_controllable_sim(config_path)
end

export run_stage0_seeding, run_constellation_design, run_stage2_verification, run_capo_pipeline
export run_debris_controllable_sim, parse_debris_args

# Utility exports
using .ConstellationUtils
using .ConstellationSlots
using .ConstellationPlotting
using .Stage0Seeding
using .Stage1Design
using .Stage2Verification

end # module ConstellationDesign
