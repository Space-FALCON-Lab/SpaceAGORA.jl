include("Misc.jl")

using PlotlyJS
using LaTeXStrings
using Arrow
using DataFrames

using PlotlyJS: circle, scatter, scatter3d, surface, plot, Layout, savefig, attr, traces, Trace

function plots_multi(state_dict, m_dict, name_dict, args, temp_name_dict)
    data_table = Dict{Any,DataFrame}()
    for (key, _) in temp_name_dict
        data_table[key] = DataFrame(Arrow.Table(temp_name_dict[key]))
    end

    if args[:keplerian] == true
        println("keplerian true")
        traces3D_all = PlotlyJS.GenericTrace[]  # Use PlotlyJS.GenericTrace instead of Trace
        traces2D_all = PlotlyJS.GenericTrace[]
        layout3D = nothing
        layout2D = nothing


        for (key, state) in pairs(state_dict)
            t3, L3 = traj_3D(state, m_dict[key], name_dict[key], args, data_table[key])
            t2, L2 = traj_2D(state, m_dict[key], name_dict[key], args, data_table[key])
            append!(traces3D_all, t3)
            append!(traces2D_all, t2)
            layout3D === nothing && (layout3D = L3)
            layout2D === nothing && (layout2D = L2)
        end

        display(PlotlyJS.plot(traces3D_all, layout3D))
        display(PlotlyJS.plot(traces2D_all, layout2D))

    else
        println("not keplerian")
        traces2D_all = PlotlyJS.GenericTrace[]
        for (key, state) in pairs(state_dict)
            t2, _ = traj_2D(state, m_dict[key], name_dict[key], args, data_table[key])
            append!(traces2D_all, t2)
        end
        layout = Layout(title="Combined 2D Trajectories", template="simple_white", showlegend=true)
        display(PlotlyJS.plot(traces2D_all, layout=layout))
    end
end
function traj_3D(state, m, name, args, data_table)
    x = data_table.pos_ii_1 * 1e-3
    y = data_table.pos_ii_2 * 1e-3
    z = data_table.pos_ii_3 * 1e-3

    r = m.planet.Rp_e*1e-3
    n = 100
    u = range(-pi, pi; length=n)
    v = range(0, pi; length=n)
    xs = [r * cos(ui) * sin(vj) for ui in u, vj in v]
    ys = [r * sin(ui) * sin(vj) for ui in u, vj in v]
    zs = [r * cos(vj) for ui in u, vj in v]
    sphere1 = surface(x=xs, y=ys, z=zs, opacity=0.9, showscale=false,
                      surfacecolor=xs.^2 .+ ys.^2 .+ zs.^2 .+ 100, name="Planet")

    r2 = r + args[:EI]
    xs2 = [r2 * cos(ui) * sin(vj) for ui in u, vj in v]
    ys2 = [r2 * sin(ui) * sin(vj) for ui in u, vj in v]
    zs2 = [r2 * cos(vj) for ui in u, vj in v]
    sphere2 = surface(x=xs2, y=ys2, z=zs2, opacity=0.2, showscale=false,
                      surfacecolor=xs2.^2 .+ ys2.^2 .+ zs2.^2 .- 100, name="EI shell")

    # split into continuous segments (optional)
    idxs = [1]
    for k in 1:length(data_table.number_of_passage)-1
        if data_table.number_of_passage[k+1] > data_table.number_of_passage[k]
            push!(idxs, k)
        end
    end
    push!(idxs, length(data_table.number_of_passage))

    traj_traces = PlotlyJS.GenericTrace[]
    if data_table.number_of_passage[end] == 1
        push!(traj_traces, scatter3d(x=x, y=y, z=z, mode="lines",
                                     line=attr(width=3), name=String(name)))
    else
        for s in 1:length(idxs)-1
            xs = x[idxs[s]:idxs[s+1]]
            ys = y[idxs[s]:idxs[s+1]]
            zs = z[idxs[s]:idxs[s+1]]
            push!(traj_traces, scatter3d(x=xs, y=ys, z=zs, mode="lines",
                                         line=attr(width=3), name=String(name)))
        end
    end

    axis_min = minimum((minimum(x), minimum(y), minimum(z), -(r2+250)))
    axis_max = maximum((maximum(x), maximum(y), maximum(z),  (r2+250)))

    layout = Layout(
        scene_aspectmode="cube",
        scene=attr(
            xaxis=attr(range=[axis_min, axis_max], title="x [km]"),
            yaxis=attr(range=[axis_min, axis_max], title="y [km]"),
            zaxis=attr(range=[axis_min, axis_max], title="z [km]")
        ),
        template="simple_white",
        showlegend=true
    )

    # include planet shells + trajectory traces
    return [sphere1, sphere2, traj_traces...], layout
end

function traj_2D(state, m, name, args, data_table)
    x = data_table.pos_ii_1 * 1e-3
    y = data_table.pos_ii_2 * 1e-3
    z = data_table.pos_ii_3 * 1e-3

    i = data_table.i[1]; Ω = data_table.OMEGA[1]; ω = data_table.omega[1]
    T_ijk = [cos(Ω)*cos(ω)-sin(Ω)*sin(ω)*cos(i)  sin(Ω)*cos(ω)+cos(Ω)*sin(ω)*cos(i)  sin(ω)*sin(i);
             -cos(Ω)*sin(ω)-sin(Ω)*cos(ω)*cos(i) -sin(Ω)*sin(ω)+cos(Ω)*cos(ω)*cos(i) cos(ω)*sin(i);
              sin(Ω)*sin(i)                       -cos(Ω)*sin(i)                       cos(i)]

    vecs = T_ijk * hcat(x, y, z)'
    xmin = minimum(vecs[1,:]); xmax = maximum(vecs[1,:])
    ymin = minimum(vecs[2,:]); ymax = maximum(vecs[2,:])
    lo = minimum((xmin, ymin)); hi = maximum((xmax, ymax))

    tr = scatter(x=vecs[1,:], y=vecs[2,:], mode="lines",
                 line=attr(width=2), name=String(name))  # give legend name

    layout = Layout(
        width=800, height=800,
        xaxis_title="x [km]", yaxis_title="y [km]",
        xaxis_range=[lo, hi], yaxis_range=[lo, hi],
        template="simple_white", showlegend=true,
        shapes=[
            circle(xref="x", yref="y", fillcolor="OrangeRed",
                   x0=-m.planet.Rp_e*1e-3, y0=-m.planet.Rp_p*1e-3,
                   x1= m.planet.Rp_e*1e-3, y1= m.planet.Rp_p*1e-3,
                   line_color="OrangeRed"),
            circle(xref="x", yref="y", fillcolor="Yellow", opacity=0.2,
                   x0=-(m.planet.Rp_e*1e-3 + args[:EI]),
                   y0=-(m.planet.Rp_e*1e-3 + args[:EI]),
                   x1=  m.planet.Rp_e*1e-3 + args[:EI],
                   y1=  m.planet.Rp_e*1e-3 + args[:EI],
                   line_color="Yellow")
        ]
    )
    return [tr], layout
end




# function plots_multi(state_dict, m_dict, name_dict, args, temp_name_dict)
#     data_table = Dict{Any,DataFrame}()
#     println("Keys in temp_name_dict: ", keys(temp_name_dict))
#     for (key, _) in temp_name_dict
#         println("Processing key for data table: $key")
#         data_table[key] = DataFrame(Arrow.Table(temp_name_dict[key]))
#     end

#     if args[:keplerian] == true

#         plots_3D = []
#         plots_2D = []
#         layout_3D = []
#         layout_2D = []

#         for (key, state) in state_dict
#             println("Processing state: $key")
#             plot_3D,layout3d = traj_3D(state, m_dict[key], name_dict[key], args, data_table[key])
#             plot_2D,layout2d = traj_2D(state, m_dict[key], name_dict[key], args, data_table[key])
#             push!(plots_3D, plot_3D)
#             push!(plots_2D, plot_2D)
#             push!(layout_3D, layout3d)
#             push!(layout_2D, layout2d)

#         end
#         traces_3D = [trace for plot_group in plots_3D for trace in plot_group]
#         traces_2D = [trace for plot_group in plots_2D for trace in plot_group]
#         combined_3D = PlotlyJS.plot(traces_3D, layout=layout_3D[1])
#         combined_2D = PlotlyJS.plot(traces_2D, layout=layout_2D[1])

#         display(combined_3D)
#         display(combined_2D)
#     else
#         # traj_3D(state, m, name, args)
        
#         plots_2D = []

#         for (key, state) in state_dict
#             println("Processing state: $key")
#             plot_2D = traj_2D(state, m_dict[key], name_dict[key], args, data_table[key])
#             push!(plots_2D, plot_2D)
#         end

#         combined_2D = PlotlyJS.plot(plots_2D..., layout=Layout(title="Combined 2D Trajectories"))

#         display(combined_2D)
#         # performance_plots(state, m, name, args, data_table)
#     end

# end

# function traj_2D(state, m, name, args, data_table)
#     x = data_table.pos_ii_1 * 1e-3
#     y = data_table.pos_ii_2 * 1e-3
#     z = data_table.pos_ii_3 * 1e-3
#     # y = [item*1e-3 for item in config.solution.orientation.pos_ii[2]]
#     # z = [item*1e-3 for item in config.solution.orientation.pos_ii[3]]
#     i = data_table.i[1]
#     Ω = data_table.OMEGA[1]
#     ω = data_table.omega[1]

#     T_ijk = [cos(Ω)*cos(ω)-sin(Ω)*sin(ω)*cos(i) sin(Ω)*cos(ω)+cos(Ω)*sin(ω)*cos(i) sin(ω)*sin(i);
#              -cos(Ω)*sin(ω)-sin(Ω)*cos(ω)*cos(i) -sin(Ω)*sin(ω)+cos(Ω)*cos(ω)*cos(i) cos(ω)*sin(i);
#              sin(Ω)*sin(i) -cos(Ω)*sin(i) cos(i)]
    
#     vector = T_ijk * hcat(x, y, z)'

#     x_min = minimum(vector[1,:])
#     x_max = maximum(vector[1,:])
#     y_min = minimum(vector[2,:])
#     y_max = maximum(vector[2,:])

#     min = minimum([x_min, y_min])
#     max = maximum([x_max, y_max])

#     x_labels = range(start=min, step=5000, stop=max)
#     y_labels = range(start=min, step=5000, stop=max)
    
#     plot_traces = scatter(x=vector[1,:], y=vector[2,:], mode="lines", line=attr(color="black"))
#     layout = Layout(width=800, height=800, xaxis_title="x [km]", yaxis_title="y [km]", xaxis_range=[min, max], yaxis_range=[min, max], 
#                     shapes=[circle(xref="x", yref="y", fillcolor="OrangeRed", x0=-m.planet.Rp_e*1e-3, y0=-m.planet.Rp_p*1e-3, x1=m.planet.Rp_e*1e-3, y1=m.planet.Rp_p*1e-3, line_color="OrangeRed"), 
#                     circle(xref="x", yref="y", fillcolor="Yellow", opacity=0.2, x0=-m.planet.Rp_e*1e-3-args[:EI], y0=-m.planet.Rp_e*1e-3-args[:EI], x1=m.planet.Rp_e*1e-3+args[:EI], y1=m.planet.Rp_e*1e-3+args[:EI], line_color="Yellow")], 
#                     template="simple_white", showlegend=false) # , xaxis_tickvals=x_labels, yaxis_tickvals=y_labels)
#     p = plot(plot_traces, layout)
#     # display(p)
#     # savefig(p, name * "_traj2D.pdf", width=800, height=800, format="pdf")
#     return plot_traces, layout
#     # display(p)
#     # savefig(p, name * "_traj2D.pdf", width=800, height=800, format="pdf")
# end

# function traj_3D(state, m, name, args, data_table)

#     x = data_table.pos_ii_1 * 1e-3
#     y = data_table.pos_ii_2 * 1e-3
#     z = data_table.pos_ii_3 * 1e-3
#     # y = [item*1e-3 for item in config.solution.orientation.pos_ii[2]]
#     # z = [item*1e-3 for item in config.solution.orientation.pos_ii[3]]

#     planet_radius = m.planet.Rp_e*1e-3
#     x_min = min(minimum(x), -(planet_radius + 250))
#     x_max = max(maximum(x), (planet_radius + 250))
#     y_min = min(minimum(y), -(planet_radius + 250))
#     y_max = max(maximum(y), (planet_radius + 250))
#     z_min = min(minimum(z), -(planet_radius + 250))
#     z_max = max(maximum(z), (planet_radius + 250))

#     axis_min = minimum([x_min, y_min, z_min])
#     axis_max = maximum([x_max, y_max, z_max])

#     r = m.planet.Rp_e*1e-3

#     n = 100
#     u = range(start=-pi, stop=pi, length=n)
#     v = range(start=0, stop=pi, length=n)
#     xs = r * cos.(u) * sin.(v)'
#     ys = r * sin.(u) * sin.(v)'
#     zs = r * ones(n) * cos.(v)'

#     sphere1 = surface(x=xs, y=ys, z=zs, opacity=0.9, showscale=false, surfacecolor=@. xs^2 + ys^2 + zs^2 + 100)

#     r += args[:EI]                       # AE alt 160 km    
#     xs = r * cos.(u) * sin.(v)'
#     ys = r * sin.(u) * sin.(v)'
#     zs = r * ones(n) * cos.(v)'

#     sphere2 = surface(x=xs, y=ys, z=zs, opacity=0.2, showscale=false, surfacecolor=@. xs^2 + ys^2 + zs^2 - 100) # "rgba(255, 255, 0, 0.25)")

#     index = [1]
#     for i in range(start=1, step=1, stop=length(data_table.number_of_passage)-1)
#         if data_table.number_of_passage[i+1] - data_table.number_of_passage[i] > 0
#             append!(index, i)
#         end
#     end

#     append!(index, length(data_table.number_of_passage))

#     x_lables = range(start=axis_min, step=5000, stop=axis_max)
#     y_labels = range(start=axis_min, step=5000, stop=axis_max)
#     z_labels = range(start=axis_min, step=5000, stop=axis_max)

#     traj_3D_traces = []
#     for i in range(start=1, step=1, stop=length(index)-1)
#         if data_table.number_of_passage[end] == 1 
#             x_s = x
#             y_s = y
#             z_s = z
#         else
#             x_s = x[index[i]:index[i+1]]
#             y_s = y[index[i]:index[i+1]]
#             z_s = z[index[i]:index[i+1]]
#         end

#         push!(traj_3D_traces, scatter3d(x=x_s, y=y_s, z=z_s, mode="lines", line=attr(color="black"), row=1, col=1))
#     end
#     if args[:type_of_mission] == "Entry"
#         push!(traj_3D_traces, scatter3d(x=[x[1]], y=[y[1]], z=[z[1]], mode="markers", marker=attr(color="black", marker="x"), row=1, col=1))
#     end
#     layout = Layout(scene_aspectmode="cube", scene_xaxis_range=[axis_min, axis_max], scene_yaxis_range=[axis_min, axis_max], scene_zaxis_range=[axis_min, axis_max], xaxis_title="x [km]", yaxis_title="y [km]", zaxis_title="z [km]", template="simple_white", showlegend=false)
#     p = plot([sphere1, sphere2, traj_3D_traces...], layout)
#     # display(p)
#     # savefig(p, name * "_traj3D.pdf", format="pdf")
#     return traj_3D_traces, layout
# end
