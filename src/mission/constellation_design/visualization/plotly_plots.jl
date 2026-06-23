module PlotlyPlots

using PlotlyJS
using LinearAlgebra

"""
    plot_constellation_3d(satellite_positions::AbstractArray{<:Real,3}, client_positions::AbstractArray{<:Real,3};
                          title::AbstractString="Constellation 3D", 
                          satellite_names::AbstractVector{<:AbstractString}=String[],
                          client_names::AbstractVector{<:AbstractString}=String[]) -> Plot

Create a 3D plot of constellation satellites and client debris.

# Arguments
- `satellite_positions::AbstractArray{<:Real,3}`: Satellite positions [3, N, M]
- `client_positions::AbstractArray{<:Real,3}`: Client positions [3, N, P]
- `title::AbstractString="Constellation 3D"`: Plot title
- `satellite_names::AbstractVector{<:AbstractString}`: Satellite names
- `client_names::AbstractVector{<:AbstractString}`: Client names

# Returns
- `Plot`: PlotlyJS Plot object
"""
function plot_constellation_3d(satellite_positions::AbstractArray{<:Real,3}, 
                               client_positions::AbstractArray{<:Real,3};
                               title::AbstractString="Constellation 3D",
                               satellite_names::AbstractVector{<:AbstractString}=String[],
                               client_names::AbstractVector{<:AbstractString}=String[])
    traces = PlotlyJS.GenericTrace[]
    
    # Plot satellite trajectories
    M = size(satellite_positions, 3)
    for m in 1:M
        sat_path = satellite_positions[:, :, m]
        name = length(satellite_names) >= m ? satellite_names[m] : "Sat $m"
        
        push!(traces, scatter3d(
            x=sat_path[1, :],
            y=sat_path[2, :],
            z=sat_path[3, :],
            mode="lines+markers",
            name=name,
            marker=attr(size=4),
            line=attr(width=2),
        ))
    end
    
    # Plot client trajectories
    P = size(client_positions, 3)
    for p in 1:P
        client_path = client_positions[:, :, p]
        name = length(client_names) >= p ? client_names[p] : "Client $p"
        
        push!(traces, scatter3d(
            x=client_path[1, :],
            y=client_path[2, :],
            z=client_path[3, :],
            mode="lines",
            name=name,
            line=attr(width=2, dash="dash"),
        ))
    end
    
    layout = Layout(
        title=title,
        scene=attr(
            xaxis_title="X [m]",
            yaxis_title="Y [m]",
            zaxis_title="Z [m]",
            aspectmode="data",
        ),
        showlegend=true,
    )
    
    return Plot(traces, layout)
end

"""
    plot_tube_3d(z_tube::AbstractArray{<:Real,3}, keepout_dirs::AbstractMatrix{<:Real};
                 title::AbstractString="Controllable Tube",
                 scale::Real=1.0) -> Plot

Create a 3D plot of the controllable tube.

# Arguments
- `z_tube::AbstractArray{<:Real,3}`: Tube sizes [H, P, Kd]
- `keepout_dirs::AbstractMatrix{<:Real}`: Direction vectors [3, Kd]
- `title::AbstractString="Controllable Tube"`: Plot title
- `scale::Real=1.0`: Scale factor for visualization

# Returns
- `Plot`: PlotlyJS Plot object
"""
function plot_tube_3d(z_tube::AbstractArray{<:Real,3}, keepout_dirs::AbstractMatrix{<:Real};
                     title::AbstractString="Controllable Tube",
                     scale::Real=1.0)
    traces = PlotlyJS.GenericTrace[]
    
    H, P, Kd = size(z_tube)
    
    # Plot tube surface for each client and horizon
    for p in 1:P, n in 1:H
        # Create tube surface points
        tube_points = zeros(3, Kd)
        for q in 1:Kd
            d = keepout_dirs[:, q]
            tube_points[:, q] = d * z_tube[n, p, q] * scale
        end
        
        push!(traces, scatter3d(
            x=tube_points[1, :],
            y=tube_points[2, :],
            z=tube_points[3, :],
            mode="lines+markers",
            name="Client $p, Horizon $n",
            marker=attr(size=3),
        ))
    end
    
    layout = Layout(
        title=title,
        scene=attr(
            xaxis_title="X [m]",
            yaxis_title="Y [m]",
            zaxis_title="Z [m]",
            aspectmode="data",
        ),
        showlegend=true,
    )
    
    return Plot(traces, layout)
end

"""
    plot_access_matrix(F::AbstractArray{<:Real,3};
                       title::AbstractString="Access Matrix",
                       client_names::AbstractVector{<:AbstractString}=String[],
                       satellite_names::AbstractVector{<:AbstractString}=String[]) -> Plot

Create a heatmap of the access matrix.

# Arguments
- `F::AbstractArray{<:Real,3}`: Access matrix [M, P, N]
- `title::AbstractString="Access Matrix"`: Plot title
- `client_names::AbstractVector{<:AbstractString}`: Client names
- `satellite_names::AbstractVector{<:AbstractString}`: Satellite names

# Returns
- `Plot`: PlotlyJS Plot object
"""
function plot_access_matrix(F::AbstractArray{<:Real,3};
                            title::AbstractString="Access Matrix",
                            client_names::AbstractVector{<:AbstractString}=String[],
                            satellite_names::AbstractVector{<:AbstractString}=String[])
    M, P, N = size(F)
    
    # Aggregate access over time (sum over N)
    F_agg = sum(F, dims=3)[:, :, 1]
    
    # Create labels
    x_labels = length(client_names) >= P ? client_names : ["Client $p" for p in 1:P]
    y_labels = length(satellite_names) >= M ? satellite_names : ["Sat $m" for m in 1:M]
    
    trace = heatmap(
        z=F_agg,
        x=x_labels,
        y=y_labels,
        colorscale="Viridis",
    )
    
    layout = Layout(
        title=title,
        xaxis_title="Client",
        yaxis_title="Satellite",
    )
    
    return Plot(trace, layout)
end

"""
    plot_stage2_trajectories(x_opt::AbstractArray{<:Real,3}, u_opt::AbstractArray{<:Real,3};
                              title::AbstractString="Stage 2 Trajectories") -> Plot

Create plots of Stage 2 verification trajectories and controls.

# Arguments
- `x_opt::AbstractArray{<:Real,3}`: State trajectories [6, N+1, P]
- `u_opt::AbstractArray{<:Real,3}`: Control inputs [3, N, P]
- `title::AbstractString="Stage 2 Trajectories"`: Plot title

# Returns
- `Plot`: PlotlyJS Plot object with subplots
"""
function plot_stage2_trajectories(x_opt::AbstractArray{<:Real,3}, u_opt::AbstractArray{<:Real,3};
                                   title::AbstractString="Stage 2 Trajectories")
    traces = PlotlyJS.GenericTrace[]
    
    N = size(x_opt, 2) - 1
    P = size(x_opt, 3)
    t = collect(0:N)
    
    # Plot position trajectories
    for p in 1:P
        push!(traces, scatter(
            x=t,
            y=x_opt[1, :, p],
            mode="lines",
            name="Client $p X",
        ))
        push!(traces, scatter(
            x=t,
            y=x_opt[2, :, p],
            mode="lines",
            name="Client $p Y",
        ))
        push!(traces, scatter(
            x=t,
            y=x_opt[3, :, p],
            mode="lines",
            name="Client $p Z",
        ))
    end
    
    layout = Layout(
        title=title,
        xaxis_title="Time step",
        yaxis_title="Position [m]",
        showlegend=true,
    )
    
    return Plot(traces, layout)
end

"""
    plot_beta_histogram(beta::AbstractVector{<:Real};
                        title::AbstractString="Satellite Selection",
                        threshold::Real=0.1) -> Plot

Create a histogram of satellite selection variables β.

# Arguments
- `beta::AbstractVector{<:Real}`: Satellite selection variables
- `title::AbstractString="Satellite Selection"`: Plot title
- `threshold::Real=0.1`: Selection threshold

# Returns
- `Plot`: PlotlyJS Plot object
"""
function plot_beta_histogram(beta::AbstractVector{<:Real};
                             title::AbstractString="Satellite Selection",
                             threshold::Real=0.1)
    M = length(beta)
    
    trace = histogram(
        x=beta,
        nbinsx=50,
        name="β values",
    )
    
    # Add threshold line
    threshold_line = shape(
        type="line",
        x0=threshold,
        x1=threshold,
        y0=0,
        y1=M / 5,
        line=attr(color="red", dash="dash"),
    )
    
    layout = Layout(
        title=title,
        xaxis_title="β",
        yaxis_title="Count",
        shapes=[threshold_line],
    )
    
    return Plot(trace, layout)
end

"""
    save_plot(plot::Plot, path::AbstractString; format::AbstractString="html")

Save a PlotlyJS plot to file.

# Arguments
- `plot::Plot`: PlotlyJS Plot object
- `path::AbstractString`: Output file path
- `format::AbstractString="html"`: Output format ("html", "png", "svg", "pdf")
"""
function save_plot(plot::Plot, path::AbstractString; format::AbstractString="html")
    if format == "html"
        savefig(plot, path)
    elseif format == "png"
        savefig(plot, path, format="png")
    elseif format == "svg"
        savefig(plot, path, format="svg")
    elseif format == "pdf"
        savefig(plot, path, format="pdf")
    else
        error("Unknown format: $format")
    end
end

export plot_constellation_3d, plot_tube_3d, plot_access_matrix,
       plot_stage2_trajectories, plot_beta_histogram, save_plot

end # module PlotlyPlots
