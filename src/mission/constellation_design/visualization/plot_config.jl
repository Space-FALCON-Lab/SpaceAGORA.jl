module PlotConfig

using ..PlotlyPlots
using ..ConstellationConfig
using ..ConstellationUtils

"""
    generate_constellation_plots(config_dict::AbstractDict, results::AbstractDict) -> Dict{String,Any}

Generate all constellation design plots based on configuration and results.

# Arguments
- `config_dict::AbstractDict`: Configuration dictionary
- `results::AbstractDict`: Pipeline results dictionary

# Returns
- `Dict{String,Any}`: Dictionary of plot file paths
"""
function generate_constellation_plots(config_dict::AbstractDict, results::AbstractDict)
    plotting_params = config_dict["plotting"]
    save_plots = plotting_params["save"]
    display_plots = plotting_params["display"]
    base_outdir = plotting_params["base_outdir"]
    
    plot_paths = Dict{String,Any}()
    
    # Create output directory
    if save_plots
        mkpath(base_outdir)
    end
    
    # Stage 1 plots
    if haskey(results, "stage1") && plotting_params["plot_stage1_lp_warm_start_beta_histogram"]
        stage1_result = results["stage1"]
        if haskey(stage1_result, "beta")
            beta_plot = plot_beta_histogram(
                stage1_result["beta"];
                title="Stage 1 Satellite Selection",
                threshold=Float64(config_dict["optimizer_params"]["beta_threshold"]),
            )
            plot_path = joinpath(base_outdir, "stage1_beta_histogram.html")
            if save_plots
                save_plot(beta_plot, plot_path)
                plot_paths["stage1_beta_histogram"] = plot_path
            end
            if display_plots
                # Display plot (implementation depends on environment)
            end
        end
    end
    
    # Tube 3D plot
    if haskey(results, "stage1") && plotting_params["plot_tube_3d"]
        stage1_result = results["stage1"]
        if haskey(stage1_result, "z_tube")
            Kd = Int(config_dict["optimizer_params"]["controllable_set_dirs"])
            keepout_dirs = _build_keepout_position_dirs(Kd)
            tube_plot = plot_tube_3d(
                stage1_result["z_tube"],
                keepout_dirs;
                title="Controllable Tube",
            )
            plot_path = joinpath(base_outdir, "tube_3d.html")
            if save_plots
                save_plot(tube_plot, plot_path)
                plot_paths["tube_3d"] = plot_path
            end
        end
    end
    
    # Access matrix plot
    if haskey(results, "stage1") && plotting_params["plot_access_matrix"]
        stage1_result = results["stage1"]
        if haskey(stage1_result, "access_matrix")
            access_plot = plot_access_matrix(
                stage1_result["access_matrix"]["F"];
                title="Access Matrix",
            )
            plot_path = joinpath(base_outdir, "access_matrix.html")
            if save_plots
                save_plot(access_plot, plot_path)
                plot_paths["access_matrix"] = plot_path
            end
        end
    end
    
    # Stage 2 trajectory plots
    if haskey(results, "stage2") && plotting_params["plot_stage2_trajectories"]
        stage2_result = results["stage2"]
        if haskey(stage2_result, "segments") && !isempty(stage2_result["segments"])
            for (idx, segment) in enumerate(stage2_result["segments"])
                if haskey(segment, "x_opt") && haskey(segment, "u_opt")
                    traj_plot = plot_stage2_trajectories(
                        segment["x_opt"],
                        segment["u_opt"];
                        title="Stage 2 Trajectory - Horizon $idx",
                    )
                    plot_path = joinpath(base_outdir, "stage2_trajectory_horizon_$idx.html")
                    if save_plots
                        save_plot(traj_plot, plot_path)
                        plot_paths["stage2_trajectory_horizon_$idx"] = plot_path
                    end
                end
            end
        end
    end
    
    constellation_log("plotting", "Generated $(length(plot_paths)) plots")
    
    return plot_paths
end

"""
    _build_keepout_position_dirs(n::Int) -> Matrix{Float64}

Build keepout position direction bank (wrapper for SupportFunctions).

# Arguments
- `n::Int`: Number of directions

# Returns
- `Matrix{Float64}`: Direction matrix
"""
function _build_keepout_position_dirs(n::Int)
    using ..SupportFunctions
    return _build_keepout_position_dirs(n)
end

export generate_constellation_plots

end # module PlotConfig
