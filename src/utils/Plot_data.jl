include("Misc.jl")

using PlotlyJS
using LaTeXStrings

 # import .config

function plots(state, m, name, args)
    if args[:keplerian] == true
        traj_3D(state, m, name, args)
        traj_2D(state, m, name, args)
    else
        traj_3D(state, m, name, args)
        traj_2D(state, m, name, args)

        performance_plots(state, m, name, args)
    end

    if args[:body_shape] == "Spacecraft" && !config.cnf.impact
        if lowercase(m.planet.name) == "mars"
            closed_form_solution_plot(name, m)
        end
        angle_of_attack_plot(name, args)
    end

    if args[:Odyssey_sim] == 1 || args[:vex_sim] == 1 || args[:magellan_sim] == 1
        ABM_periapsis(name)
    end

end

function angle_of_attack_plot(name, args)
    alt_idx = findall(x -> x < args[:AE]*1e3, config.solution.orientation.alt)

    index_orbit = [1]
    time_0 = [config.solution.orientation.time[alt_idx[1]]]

    for i in range(start=2, step=1, stop=length(alt_idx))
        if alt_idx[i] - alt_idx[i - 1] > 2
            append!(index_orbit, i)
            append!(time_0, config.solution.orientation.time[alt_idx[i-1]])
        end
    end

    append!(index_orbit, length(alt_idx))

    # println(alt_idx)
    # println(alt_idx[index_orbit])
    # println(index_orbit)
    # println(length(alt_idx))

    if length(index_orbit) <= 2 # If onlly 1 orbit
        time = [config.solution.orientation.time[i] for i in alt_idx]
        aoa = [rad2deg(config.solution.physical_properties.α[i]) for i in alt_idx]
        trace1 = scatter(x=time, y=aoa, mode="lines", line=attr(color="black"))
        layout = Layout(xaxis_title="Time [s]", yaxis_title="α [deg]", template="simple_white", showlegend=false)
        p = plot(trace1, layout)
    else
        time_end = 0
        x_labels = []
        aoa_end = []
        plots_aoa_line = []
        plots_aoa_mark = []

        for i in range(start=1, step=1, stop=length(index_orbit)-1)
            time = [config.solution.orientation.time[j] for j in alt_idx[index_orbit[i]:index_orbit[i+1]-1]]
            aoa = [rad2deg(config.solution.physical_properties.α[j]) for j in alt_idx[index_orbit[i]:index_orbit[i+1]-1]]

            push!(plots_aoa_mark, scatter(x=[i], y=[aoa[end]], mode="markers", marker=attr(color="red")))
            push!(plots_aoa_line, scatter(x=(time .- time[1])./(time[end] - time[1]) .+ (i-1), y=aoa, mode="lines", line=attr(color="black")))
            time_end = time[end]
            # append!(x_labels, (time_end + time[1])/2)
            append!(x_labels, i)
            append!(aoa_end, aoa[end])
        end

        # push!(plots_aoa_line, scatter(x=x_labels, y=aoa_end, mode="lines", line=attr(color="black")))
        layout = Layout(xaxis_title="Orbits", yaxis_title="α [deg]", template="simple_white", showlegend=false)
        p = plot([plots_aoa_line..., plots_aoa_mark...], layout)
    end

    display(p)
    savefig(p, name * "_angle_of_attack_profile.pdf", format="pdf")
end

function closed_form_solution_plot(name, mission)
    println("alt: $(config.solution.orientation.alt)")

    alt = [item - mission.planet.Rp_e for item in config.solution.orientation.pos_ii_mag]
    alt_idx = findall(x -> x <= args[:EI]*1e3, alt)

    index_orbit = [1]
    time_0 = [config.solution.orientation.time[alt_idx[1]]]

    for i in range(start=2, step=1, stop=length(alt_idx))
        if alt_idx[i] - alt_idx[i - 1] > 50
            append!(index_orbit, i)
            append!(time_0, config.solution.orientation.time[alt_idx[i]])
        end
    end
    append!(index_orbit, length(alt_idx))

    alt_idx_cf = findall(x -> (x > 0) && (x <= args[:EI]*1e3), config.solution.closed_form.h_cf)
    index_orbit_cf = [1]
    time_0_cf = [config.solution.closed_form.t_cf[alt_idx_cf[1]]]
    
    for i in range(start=2, step=1, stop=length(alt_idx_cf))
        if alt_idx_cf[i] - alt_idx_cf[i - 1] > 50
            append!(index_orbit_cf, i)
            append!(time_0_cf, config.solution.closed_form.t_cf[alt_idx_cf[i]])
        end
    end
    append!(index_orbit_cf, length(alt_idx_cf))

    plot_traces_alt = []
    plot_traces_gamma = []
    plot_traces_v = []

    plot_traces_alt_cf = []
    plot_traces_gamma_cf = []
    plot_traces_v_cf = []

    for i in range(start=1, step=1, stop=length(index_orbit)-1)
        time_end = 0

        time = [config.solution.orientation.time[j] - time_0[i] + time_end for j in alt_idx[index_orbit[i]:index_orbit[i+1]-1]]
        alt = [(config.solution.orientation.pos_ii_mag[j] - mission.planet.Rp_e)/1e3 for j in alt_idx[index_orbit[i]:index_orbit[i+1]-1]]
        gamma = [rad2deg(config.solution.orientation.γ_ii[j]) for j in alt_idx[index_orbit[i]:index_orbit[i+1]-1]]
        v = [config.solution.orientation.vel_ii_mag[j] for j in alt_idx[index_orbit[i]:index_orbit[i+1]-1]]

        time_cf = [config.solution.closed_form.t_cf[j] - time_0_cf[i] + time_end for j in alt_idx_cf[index_orbit_cf[i]:index_orbit_cf[i+1]-1]]
        alt_cf = [config.solution.closed_form.h_cf[j]/1e3 for j in alt_idx_cf[index_orbit_cf[i]:index_orbit_cf[i+1]-1]]
        gamma_cf = [rad2deg(config.solution.closed_form.γ_cf[j]) for j in alt_idx_cf[index_orbit_cf[i]:index_orbit_cf[i+1]-1]]
        v_cf = [config.solution.closed_form.v_cf[j] for j in alt_idx_cf[index_orbit_cf[i]:index_orbit_cf[i+1]-1]]

        push!(plot_traces_alt, scatter(x=time, y=alt, mode="lines", line=attr(color="black")))
        push!(plot_traces_gamma, scatter(x=time, y=gamma, mode="lines", line=attr(color="black")))
        push!(plot_traces_v, scatter(x=time, y=v, mode="lines", line=attr(color="black")))

        push!(plot_traces_alt_cf, scatter(x=time_cf, y=alt_cf, mode="lines", line=attr(color="gray")))
        push!(plot_traces_gamma_cf, scatter(x=time_cf, y=gamma_cf, mode="lines", line=attr(color="gray")))
        push!(plot_traces_v_cf, scatter(x=time_cf, y=v_cf, mode="lines", line=attr(color="gray")))
    end

    layout_alt = Layout(xaxis_title="Time [s]", yaxis_title="Altitude [km]")
    layout_gamma = Layout(xaxis_title="Time [s]", yaxis_title="γ [deg]")
    layout_v = Layout(xaxis_title="Time [s]", yaxis_title="Velocity [km/s]")

    p_alt = plot([plot_traces_alt..., plot_traces_alt_cf...], layout_alt)
    p_gamma = plot([plot_traces_gamma..., plot_traces_gamma_cf...], layout_gamma)
    p_v = plot([plot_traces_v..., plot_traces_v_cf...], layout_v)

    p = [p_alt p_gamma p_v]
    relayout!(p, width=2200, height=1000, template="simple_white", showlegend=false)

    display(p)
    savefig(p, name * "_closed_form_solution.pdf", format="pdf")
end

function performance_plots(state, m, name, args)
    if args[:body_shape] == "Spacecraft"
        plot_traces_1 = scatter(x=[item/(60*60*24) for item in config.solution.orientation.time], y=[item * 1e-6 for item in config.solution.forces.energy], mode="lines", line=attr(color="black"))
        layout_1 = Layout(xaxis_title="Time [days]", yaxis_title="Energy [MJ/kg]")
    else
        index = findall(x -> x < 200*1e3, config.solution.orientation.alt)
        alt = config.solution.orientation.alt[index]

        index_orbit = [1]
        time_0 = [config.solution.orientation.time[index[1]]]

        for i in range(start=2, step=1, stop=length(index))
            if index[i] - index[i - 1] > 50
                append!(index_orbit, i)
                append!(time_0, config.solution.orientation.time[index[i]])
            end
        end
        append!(index_orbit, length(index))

        if length(index_orbit) <= 2
            time = [config.solution.orientation.time[i] for i in index]
            plot_traces_1 = scatter(x=time/(60*60*24), y=alt, mode="markers", marker=attr(color="black"))
            layout_1 = Layout(xaxis_title="Time [days]", yaxis_title="Altitude [km]")
        end
    end

    plot_1 = plot([plot_traces_1], layout_1)

    index = findall(x -> x > 0, config.solution.performance.heat_rate)
    heat_rate = config.solution.performance.heat_rate[index]
    index_orbit = [1]
    time_0 = [config.solution.orientation.time[index[1]]]

    for i in range(start=2, step=1, stop=length(index))
        if index[i] - index[i - 1] > 50
            append!(index_orbit, i)
            append!(time_0, config.solution.orientation.time[index[i]])
        end
    end

    # append!(index_orbit, length(index))
    plot_traces_heat_rate = []

    if length(index_orbit) <= 2
        time = [config.solution.orientation.time[i] for i in index]
        push!(plot_traces_heat_rate, scatter(x=time, y=heat_rate, mode="lines", line=attr(color="black")))
    else
        time_end = 0
        for i in range(start=1, step=1, stop=length(index_orbit)-1)
            time = [config.solution.orientation.time[j] - time_0[i] + time_end for j in index[index_orbit[i]:index_orbit[i+1]-1]]
            heat_rate = [config.solution.performance.heat_rate[j] for j in index[index_orbit[i]:index_orbit[i+1]-1]]
            max_value = maximum(heat_rate)
            max_index = findfirst(x -> x == max_value, heat_rate)

            # push!(plot_traces_heat_rate, scatter(x=[time[max_index]] , y=[heat_rate[max_index]], mode="markers", marker=attr(color="black"))) # uncomment for time instead
            push!(plot_traces_heat_rate, scatter(x=[i] , y=[heat_rate[max_index]], mode="markers", marker=attr(color="black")))

            time_end = time[end]
        end
    end

    layout_heat_rate = Layout(xaxis_title="Orbits", yaxis_title="Heat rate [W/cm^2]")
    plot_heat_rate = plot([plot_traces_heat_rate...], layout_heat_rate)

    plot_traces_heat_load = []

    if length(index_orbit) <= 2
        time = [config.solution.orientation.time[i] for i in index]
        heat_load = [config.solution.performance.heat_load[i] for i in index]
        push!(plot_traces_heat_load, scatter(x=time, y=heat_load, mode="markers", marker=attr(color="black")))
    else
        time_end = 0

        for i in range(start=1, step=1, stop=length(index_orbit)-1)
            time = [config.solution.orientation.time[j] - time_0[i] + time_end for j in index[index_orbit[i]:index_orbit[i+1]-1]]
            heat_load = [config.solution.performance.heat_load[j] for j in index[index_orbit[i]:index_orbit[i+1]-1]]

            # push!(plot_traces_heat_load, scatter(x=[time[end]] , y=[heat_load[end]], mode="markers", marker=attr(color="black")))
            push!(plot_traces_heat_load, scatter(x=[i] , y=[heat_load[end]], mode="markers", marker=attr(color="black")))

            time_end = time[end]
        end
    end

    layout_heat_load = Layout(xaxis_title="Orbits", yaxis_title="Heat load [J/cm^2]")
    plot_heat_load = plot([plot_traces_heat_load...], layout_heat_load)

    if args[:body_shape] == "Spacecraft"
        plot_traces_4 = scatter(x=[item/(60*60*24) for item in config.solution.orientation.time], y=[item - args[:dry_mass] for item in config.solution.performance.mass], mode="lines", line=attr(color="black"))
        layout_4 = Layout(xaxis_title="Time [days]", yaxis_title="Mass [kg]")
    else
        index = findall(x -> x > 0, config.solution.performance.q)
        q = config.solution.performance.q[index]
        index_orbit = [1]
        time_0 = [config.solution.orientation.time[index[1]]]

        for i in range(start=2, step=1, stop=length(index))
            if index[i] - index[i - 1] > 50
                append!(index_orbit, i-1)
                append!(time_0, config.solution.orientation.time[index[i-1]])
            end
        end

        append!(index_orbit, length(index))

        if length(index_orbit) <= 2
            time = [config.solution.orientation.time[i] for i in index]
            plot_traces_4 = scatter(x=time/(60*60*24), y=q, mode="markers", marker=attr(color="black"))
            layout_4 = Layout(xaxis_title="Time [days]", yaxis_title="Dynamic pressure [Pa]")
        end
    end

    plot_4 = plot([plot_traces_4], layout_4)

    p = [plot_1 plot_heat_rate; plot_heat_load plot_4]
    relayout!(p, width=2100, height=1200, template="simple_white", showlegend=false)

    display(p)
    savefig(p, name * "_performance.pdf", format="pdf")
end

function traj_2D(state, m, name, args)
    x = [item*1e-3 for item in config.solution.orientation.pos_ii[1]]
    y = [item*1e-3 for item in config.solution.orientation.pos_ii[2]]
    z = [item*1e-3 for item in config.solution.orientation.pos_ii[3]]
    i = config.solution.orientation.oe[3][1]
    Ω = config.solution.orientation.oe[4][1]
    ω = config.solution.orientation.oe[5][1]

    T_ijk = [cos(Ω)*cos(ω)-sin(Ω)*sin(ω)*cos(i) sin(Ω)*cos(ω)+cos(Ω)*sin(ω)*cos(i) sin(ω)*sin(i);
             -cos(Ω)*sin(ω)-sin(Ω)*cos(ω)*cos(i) -sin(Ω)*sin(ω)+cos(Ω)*cos(ω)*cos(i) cos(ω)*sin(i);
             sin(Ω)*sin(i) -cos(Ω)*sin(i) cos(i)]
    
    vector = T_ijk * hcat(x, y, z)'

    x_min = minimum(vector[1,:])
    x_max = maximum(vector[1,:])
    y_min = minimum(vector[2,:])
    y_max = maximum(vector[2,:])

    min = minimum([x_min, y_min])
    max = maximum([x_max, y_max])

    x_labels = range(start=min, step=5000, stop=max)
    y_labels = range(start=min, step=5000, stop=max)
    
    plot_traces = scatter(x=vector[1,:], y=vector[2,:], mode="lines", line=attr(color="black"))
    layout = Layout(width=800, height=800, xaxis_title="x [km]", yaxis_title="y [km]", xaxis_range=[min, max], yaxis_range=[min, max], 
                    shapes=[circle(xref="x", yref="y", fillcolor="OrangeRed", x0=-m.planet.Rp_e*1e-3, y0=-m.planet.Rp_p*1e-3, x1=m.planet.Rp_e*1e-3, y1=m.planet.Rp_p*1e-3, line_color="OrangeRed"), 
                    circle(xref="x", yref="y", fillcolor="Yellow", opacity=0.2, x0=-m.planet.Rp_e*1e-3-args[:EI], y0=-m.planet.Rp_e*1e-3-args[:EI], x1=m.planet.Rp_e*1e-3+args[:EI], y1=m.planet.Rp_e*1e-3+args[:EI], line_color="Yellow")], 
                    template="simple_white", showlegend=false) # , xaxis_tickvals=x_labels, yaxis_tickvals=y_labels)
    p = plot(plot_traces, layout)
    display(p)
    savefig(p, name * "_traj2D.pdf", width=800, height=800, format="pdf")
end

function traj_3D(state, m, name, args)
    # p = make_subplots(rows=3, cols=3, 
    #                   specs=[Spec(kind="scene", rowspan=3, colspan=2) missing Spec(kind="scene");
    #                          missing missing Spec(kind="scene"); 
    #                          missing missing Spec(kind="scene")])

    x = [item*1e-3 for item in config.solution.orientation.pos_ii[1]]
    y = [item*1e-3 for item in config.solution.orientation.pos_ii[2]]
    z = [item*1e-3 for item in config.solution.orientation.pos_ii[3]]

    x_min = minimum(x)
    x_max = maximum(x)
    y_min = minimum(y)
    y_max = maximum(y)
    z_min = minimum(z)
    z_max = maximum(z)

    min = minimum([x_min, y_min, z_min])
    max = maximum([x_max, y_max, z_max])

    r = m.planet.Rp_e*1e-3

    n = 100
    u = range(start=-pi, stop=pi, length=n)
    v = range(start=0, stop=pi, length=n)
    xs = r * cos.(u) * sin.(v)'
    ys = r * sin.(u) * sin.(v)'
    zs = r * ones(n) * cos.(v)'

    # println("")
    # println(size(xs))
    # println("")

    sphere1 = surface(x=xs, y=ys, z=zs, opacity=0.9, showscale=false, surfacecolor=@. xs^2 + ys^2 + zs^2 + 100)

    r += args[:EI]                       # AE alt 160 km    
    xs = r * cos.(u) * sin.(v)'
    ys = r * sin.(u) * sin.(v)'
    zs = r * ones(n) * cos.(v)'

    sphere2 = surface(x=xs, y=ys, z=zs, opacity=0.2, showscale=false, surfacecolor=@. xs^2 + ys^2 + zs^2 - 100) # "rgba(255, 255, 0, 0.25)")

    index = [1]
    for i in range(start=1, step=1, stop=length(config.solution.orientation.number_of_passage)-1)
        if config.solution.orientation.number_of_passage[i+1] - config.solution.orientation.number_of_passage[i] > 0
            append!(index, i)
        end
    end
    append!(index, length(config.solution.orientation.number_of_passage))

    x_lables = range(start=min, step=5000, stop=max)
    y_labels = range(start=min, step=5000, stop=max)
    z_labels = range(start=min, step=5000, stop=max)

    traj_3D_traces = []
    for i in range(start=1, step=1, stop=length(index)-1)
        x_s = x[index[i]:index[i+1]]
        y_s = y[index[i]:index[i+1]]
        z_s = z[index[i]:index[i+1]]
        
        push!(traj_3D_traces, scatter3d(x=x_s, y=y_s, z=z_s, mode="lines", line=attr(color="black"), row=1, col=1))
    end

    layout = Layout(scene_aspectmode="cube", scene_xaxis_range=[min, max], scene_yaxis_range=[min, max], scene_zaxis_range=[min, max], xaxis_title="x [km]", yaxis_title="y [km]", zaxis_title="z [km]", template="simple_white", showlegend=false)
    p = plot([sphere1, sphere2, traj_3D_traces...], layout)
    display(p)
    savefig(p, name * "_traj3D.pdf", format="pdf")

end

function ABM_periapsis(name)
    orbit_number = config.cnf.orbit_number_list .- 1
    periapsis_altitude = config.cnf.periapsis_list

    delta_v = [0.0]

    delta_v_raise = []
    delta_v_lower = []

    for j in range(start=1, step=1, stop=length(config.cnf.Δv_list)-1)
        append!(delta_v, config.cnf.Δv_list[j+1] - config.cnf.Δv_list[j])
    end

    for j in range(start=1, step=1, stop=length(delta_v))
        if j in config.cnf.raise_man_orbit
            append!(delta_v_raise, delta_v[j])
        elseif j in config.cnf.lower_man_orbit
            append!(delta_v_lower, delta_v[j])
        end
    end

    plot_traces_palt = scatter(x=orbit_number, y=periapsis_altitude, mode="markers", marker=attr(color="black"))

    plot_traces_abm1 = scatter(x=[item for item in config.cnf.raise_man_orbit], y=delta_v_raise, mode="markers", marker=attr(color="blue"), yaxis="y2") # , label="ABM to raise periapsis")
    plot_traces_abm2 = scatter(x=[item for item in config.cnf.lower_man_orbit], y=delta_v_lower, mode="markers", marker=attr(color="red"), yaxis="y2") # , label="ABM to lower periapsis")
    
    p = plot([plot_traces_palt, plot_traces_abm1, plot_traces_abm2], Layout(xaxis_title_text="Orbit", yaxis_title_text="Periapsis altitude [km]", yaxis2 = attr(title="ABM Magnitude [m/s]", overlaying="y", side="right"), template="simple_white", showlegend=false))

    display(p)
    savefig(p, name * "_Periapsis_alt_and_maneuvers.pdf", format="pdf")

end