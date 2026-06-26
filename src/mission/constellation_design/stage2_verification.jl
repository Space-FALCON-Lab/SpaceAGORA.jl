module Stage2Verification

using ..ConstellationUtils
using ..ConstellationSlots
using ..ConstellationPlotting
using ..Stage0Seeding
using ..Stage1Design
using ..Stage2OCPVerification
using LinearAlgebra
using Statistics
using Random
using DataFrames, CSV
using Arrow
using Dates
using ProgressMeter

"""
    generate_verification_plots(config_dict::AbstractDict) -> Dict{String,String}

Generate visualization plots for control verification.
"""
function generate_verification_plots(config_dict::AbstractDict)
    constellation_log("stage2", "Generating verification plots")
    
    plotting_params = get(config_dict, "plotting", Dict{String,Any}())
    save_plots = Bool(get(plotting_params, "save", true))
    display_plots = Bool(get(plotting_params, "display", false))
    
    plot_paths = Dict{String,String}()
    
    if save_plots
        # Generate constellation plot
        stage1_results = get(config_dict, "stage1_results", Dict{String,Any}())
        active_sats = get(stage1_results, "active_satellites", Int[])
        
        if !isempty(active_sats)
            # Create synthetic positions for plotting
            positions = rand(3, length(active_sats)) * 1e6
            
            outdir = get(plotting_params, "controllable_outdir", "plots")
            mkpath(outdir)
            
            plot_path = joinpath(outdir, "constellation_verification.png")
            plot_constellation_3d(positions; 
                                  title_str="Constellation Verification",
                                  save_path=plot_path)
            
            plot_paths["constellation"] = plot_path
        end
    end
    
    constellation_log("stage2", "Generated $(length(plot_paths)) plots")
    
    return plot_paths
end

"""
    run_stage2_verification(config_dict::AbstractString) -> Dict{String,Any}

Stage 2 control verification entry point. Loads a YAML config file and runs
control verification and plotting.

# Arguments
- `config_path::AbstractString`: Path to YAML configuration file

# Returns
- Dictionary containing verification results and plot paths
"""
function run_stage2_verification(config_path::AbstractString)
    config_dict = ingest_yaml(config_path)
    return run_stage2_verification(config_dict)
end

"""
    run_stage2_verification(config_dict::AbstractDict) -> Dict{String,Any}

Stage 2 control verification entry point. Runs OCP verification
and plotting using the provided configuration dictionary.

# Arguments
- `config_dict::AbstractDict`: Configuration dictionary

# Returns
- Dictionary containing verification results and plot paths
"""
function run_stage2_verification(config_dict::AbstractDict)
    constellation_log_init!(config_dict; context="stage2_verification")
    
    try
        constellation_log("stage2", "Starting Stage 2 OCP verification")
        
        # Ensure Stage 1 has been run
        if !haskey(config_dict, "stage1_results")
            constellation_log_warn("stage2", "Stage 1 results not found, running Stage 1 first")
            Stage1Design.run_constellation_design(config_dict)
        end
        
        # Run Stage 2 OCP verification
        constellation_log("stage2", "Running Stage 2 OCP verification")
        verification_results = run_stage2_ocp_verification(config_dict)
        
        # Generate plots
        plot_paths = generate_verification_plots(config_dict)
        
        # Store results in config_dict
        config_dict["stage2_results"] = verification_results
        config_dict["stage2_plots"] = plot_paths
        
        constellation_log("stage2", "Stage 2 completed successfully")
        
        return merge(verification_results, Dict("plot_paths" => plot_paths))
    catch err
        constellation_log_exception("stage2", err)
        rethrow(err)
    finally
        constellation_log_close!()
    end
end

export run_stage2_verification, generate_verification_plots

end # module Stage2Verification
