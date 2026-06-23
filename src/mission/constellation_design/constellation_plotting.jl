module ConstellationPlotting

using Plots: plot, plot!, scatter!, savefig, plot3d, plot3d!, xlabel!, ylabel!, title!
import Plots
using LinearAlgebra
using Statistics

const _PLOTTING_ENABLED = true

function _init_plotting_backend()
    if _PLOTTING_ENABLED
        ENV["GKS_WSTYPE"] = "100"
        Plots.gr()
    end
end

# Initialize plotting backend on module load
_init_plotting_backend()

function plot_constellation_3d(positions::Matrix{<:Real}; 
                               title_str::String="Constellation",
                               xlabel_str::String="X (m)",
                               ylabel_str::String="Y (m)",
                               zlabel_str::String="Z (m)",
                               markersize::Real=2.0,
                               save_path::Union{String,Nothing}=nothing)
    x = positions[1, :]
    y = positions[2, :]
    z = positions[3, :]
    
    plt = plot3d(x, y, z, 
                 markersize=markersize,
                 title=title_str,
                 xlabel=xlabel_str,
                 ylabel=ylabel_str,
                 zlabel=zlabel_str,
                 legend=false)
    
    if save_path !== nothing
        savefig(plt, save_path)
    end
    
    return plt
end

function plot_access_matrix(access_matrix::Matrix{<:Real};
                            title_str::String="Access Matrix",
                            save_path::Union{String,Nothing}=nothing)
    heatmap(access_matrix,
            title=title_str,
            xlabel="Satellite Index",
            ylabel="Client Index",
            colorbar=true)
    
    if save_path !== nothing
        savefig(save_path)
    end
end

export plot_constellation_3d, plot_access_matrix

end # module ConstellationPlotting
