include("Misc.jl")

using PlotlyJS
using LaTeXStrings
using Arrow
using DataFrames
using Statistics

using PlotlyJS: circle, scatter, scatter3d, surface, plot, Layout, savefig, attr, traces, GenericTrace
using PlotlyBase: Frame

# using PlotlyJS
# using PlotlyBase: Frame

function plots_multi(state_dict, m_dict, name_dict, args, temp_name_dict)
    data_table = Dict{Any,DataFrame}()
    m = m_dict[1]
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

        access_traces, access_layout = access_graph(state_dict, m, name_dict[1], args, data_table[1])

        display(PlotlyJS.plot(traces3D_all, layout3D))
        display(PlotlyJS.plot(traces2D_all, layout2D))
        display(PlotlyJS.plot(access_traces, access_layout))




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

function traj_3D(state, m, name, args, data_table,anim=false)
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
    if anim
        return traj_traces, layout
    else
        # include planet shells + trajectory traces
        return [sphere1, sphere2, traj_traces...], layout
    end
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

function access_graph(state,m,name,args,data_table)
    #sort accesses by spacecraft
    access_dict = Dict()
    # println(args[:space_objects_dict][1].laser_effector)
    # println(args[:space_objects_dict][1].access_log)
    # println(args[:space_objects_dict][2].access_log)
    # if isa(args[:space_objects_dict][1].access_log, Vector{Tuple{Float64, Float64, Int64}})
    #     print("detected!")
    # end
    for (sc_id,sc) in args[:space_objects_dict]
        # if length(sc.access_log) == 0
        #     println("Empty access log detected for spacecraft $sc_id")
        #     continue
        # end
        println(sc_id)
        if args[:space_objects_dict][sc_id].laser_effector == true
            println("laser effector detected")
            continue
        end
        sc_accesses = sc.access_log #get accesses
        accesser_id_list = []
        for access in sc_accesses #get all spacecraft that accesed the target
            # println("iterating through accesses")
            accesser_id = access[3]
            if !(accesser_id in accesser_id_list)
                push!(accesser_id_list, accesser_id)
            end #add to the accesser_id_list if not already there
        end
        #create blank time series access vector
        access_vector = zeros(length(data_table.time))
        last_time = data_table.time[end]
        first_time = data_table.time[1]
        println(first_time)
        println(last_time)
        #check every access to sum them all up for each time_step
        for access in sc_accesses
            # println("iterating through")
            access_start = access[1]
            access_end = access[2]
            #fill access vector accordingly
            time_vector = data_table.time
            for (idx, time) in enumerate(time_vector)
                # print("iterating through access")
                # fill such that higher magnitudes indicate more spacecraft covering
                # compare actual times (not indices) to avoid mismatches and double counting
                if time >= access_start && time <= access_end
                    access_vector[idx] += 1 
                    # print("adding to access/")
                end

            end
                
        end
        # println(access_vector)
        access_dict[sc_id] = access_vector
    end

    #now plot each access vector over time for each spacecraft
    time_vector = data_table.time
    # print(args[:space_objects_dict][1].access_log)
    access_traces = PlotlyJS.GenericTrace[]
    for (sc_id, access_vector) in access_dict
        # Get spacecraft name from space_objects_dict
        # println(access_vector)
        if args[:space_objects_dict][sc_id].laser_effector
            sc_name = "laser_effector"
            push!(access_traces, scatter(x=time_vector, y=access_vector, mode="lines",
                                        line = attr(width=2, shape="hv"), name="$(sc_name) (ID: $(sc_id))"))
        else
            sc_name = "target_sc"
            push!(access_traces, scatter(x=time_vector, y=access_vector, mode="lines",
                                        line = attr(width=2, shape="hv"), name="$(sc_name) (ID: $(sc_id))"))
        end
    end
    layout = Layout(
        width=800, height=600,
        xaxis_title="Time [s]", yaxis_title="Number of Accesses", 
        template="simple_white", showlegend=true,
        title="Access Graph for $(name)",
        legend=attr(
            x=1.02,  # Position legend to the right of plot
            y=1,     # Align to top
            xanchor="left",
            yanchor="top"
        )
    )
    return access_traces, layout


end

