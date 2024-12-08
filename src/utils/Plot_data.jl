include("Misc.jl")

using PlotlyJS
using LaTeXStrings

import .config

function plots(state, m, name, args)
    traj_3D(state, m, name)
    traj_2D(state, m , name)

    performance_plots(state, m, name, args)

    if args[:body_shape] == "Spacecraft" && !config.cnf.impact
        closed_form_solution_plot(name, m)
        angle_of_attack_plot(name, args)
    end

    if args[:Odyssey_sim] == 1
        ABM_periapsis(name)
    end

end

function angle_of_attack_plot(name, args)
    alt_idx = findall(x -> x < args[:AE]*1e3, config.solution.orientation.alt)

    index_orbit = [1]
    time_0 = [config.solution.orientation.time[alt_idx[1]]]

    for i in range(start=2, step=1, stop=length(alt_idx))
        if alt_idx[i] - alt_idx[i - 1] > 2
            append!(index_orbit, i-1)
            append!(time_0, config.solution.orientation.time[alt_idx[i-1]])
        end
    end

    append!(index_orbit, length(alt_idx))

    if length(index_orbit) == 2 # If onlly 1 orbit
        time = [config.solution.orientation.time[i] for i in alt_idx]
        aoa = [rad2deg(config.solution.physical_properties.α[i]) for i in alt_idx]
        trace1 = scatter(x=time, y=aoa, mode="lines", line=attr(color="black"))
        layout = Layout(xaxis_title="Time [s]", yaxis_title=L"$\alpha$ [$^\circ$]")
        plot(trace1, layout)
    else
        time_end = 0
        x_labels = []
        plots_aoa_line = []
        plots_aoa_mark = []
        

        for i in range(start=1, step=1, stop=length(index_orbit)-1)
            time = [config.solution.orientation.time[j] - time_0[i] + time_end for j in alt_idx[index_orbit[i]:index_orbit[i+1]-1]]
            aoa = [rad2deg(config.solution.physical_properties.α[j]) for j in alt_idx[index_orbit[i]:index_orbit[i+1]-1]]
            append!(plots_aoa_line, scatter(x=time, y=aoa, mode="lines", line=attr(color="black")))
            append!(plots_aoa_mark, scatter(x=time, y=aoa, mode="markers", marker=attr(color="red")))
            time_end = time[end]
            append!(x_labels, (time_end + time[1])/2)
        end

        layout = Layout(xaxis_title="Time [s]", yaxis_title=L"$\alpha$ [$^\circ$]", xaxis_tickvals=x_labels)
        plot([plots_aoa_line..., plots_aoa_mark...], layout)

        savefig(plot([plots_aoa_line..., plots_aoa_mark...], layout), name * "_angle_of_attack_profile.png")
    end
end

function closed_form_solution_plot(name, mission)
    alt = [item - mission.planet.Rp_e foe item in config.solution.orientation.pos_ii_mag]
    alt_idx = findall(x -> x <= 160*1e3, alt)

    index_orbit = [1]
    time_0 = [config.solution.orientation.time[alt_idx[1]]]

    for i in range(start=2, step=1, stop=length(alt_idx))
        if alt_idx[i] - alt_idx[i - 1] > 50
            append!(index_orbit, i-1)
            append!(time_0, config.solution.orientation.time[alt_idx[i-1]])
        end
    end
    append!(index_orbit, length(alt_idx))

    alt_idx_cf = findall(x -> (x > 0) && (x <= 160*1e3), config.solution.closed_form.h_cf)
    index_orbit_cf = [1]
    time_0_cf = [config.solution.closed_form.t_cf[alt_idx_cf[1]]]
    
    for i in range(start=2, step=1, stop=length(alt_idx_cf))
        if alt_idx_cf[i] - alt_idx_cf[i - 1] > 50
            append!(index_orbit_cf, i-1)
            append!(time_0_cf, config.solution.closed_form.t_cf[alt_idx_cf[i-1]])
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

        append!(plot_traces_alt, scatter(x=time, y=alt, mode="lines", line=attr(color="black")))
        append!(plot_traces_gamma, scatter(x=time, y=gamma, mode="lines", line=attr(color="black")))
        append!(plot_traces_v, scatter(x=time, y=v, mode="lines", line=attr(color="black")))

        append!(plot_traces_alt_cf, scatter(x=time_cf, y=alt_cf, mode="lines", line=attr(color="gray")))
        append!(plot_traces_gamma_cf, scatter(x=time_cf, y=gamma_cf, mode="lines", line=attr(color="gray")))
        append!(plot_traces_v_cf, scatter(x=time_cf, y=v_cf, mode="lines", line=attr(color="gray")))
    end

    layout_alt = Layout(xaxis_title="Time [s]", yaxis_title="Altitude [km]")
    layout_gamma = Layout(xaxis_title="Time [s]", yaxis_title=L"$\gamma$ [$^\circ$]")
    layout_v = Layout(xaxis_title="Time [s]", yaxis_title="Velocity [km/s]")

    p_alt = plot([plot_traces_alt..., plot_traces_alt_cf...], layout_alt)
    p_gamma = plot([plot_traces_gamma..., plot_traces_gamma_cf...], layout_gamma)
    p_v = plot([plot_traces_v..., plot_traces_v_cf...], layout_v)

    p = [p_alt p_gamma p_v]
    relayout!(p, width=2200, height=1000)

    savefig(p, name * "_closed_form_solution.png")
end

function performance_plots(state, m, name, args)
    if args[:body_shape] == "Spacecraft"
        plot_traces_1 = scatter(x=[item/(60*60*24) for item in config.solution.orientation.time], y=[item * 1e-6 for item in config.solution.forces.energy], mode="lines", line=attr(color="black"))
        layout_1 = Layout(xaxis_title="Time [s]", yaxis_title="Energy [MJ/kg]")
    else
        index = findall(x -> x < 200*1e3, config.solution.orientation.alt)
        alt = config.solution.orientation.alt[index]

        index_orbit = [1]
        time_0 = [config.solution.orientation.time[index[1]]]

        for i in range(start=2, step=1, stop=length(index))
            if index[i] - index[i - 1] > 50
                append!(index_orbit, i-1)
                append!(time_0, config.solution.orientation.time[index[i-1]])
            end
        end
        append!(index_orbit, length(index))

        if length(index_orbit) == 2
            time = [config.solution.orientation.time[i] for i in index]
            plot_traces_1 = scatter(x=time, y=alt, mode="markers", marker=attr(color="black"))
            layout_1 = Layout(xaxis_title="Time [s]", yaxis_title="Altitude [km]")
        end
    end

    plot_1 = plot([plot_traces_1], layout_1)

    index = findall(x -> x > 0, config.solution.performance.heat_rate)
    index_orbit = [1]
    time_0 = [config.solution.orientation.time[index[1]]]

    for i in range(start=2, step=1, stop=length(index))
        if index[i] - index[i - 1] > 50
            append!(index_orbit, i-1)
            append!(time_0, config.solution.orientation.time[index[i-1]])
        end
    end

    append!(index_orbit, length(index))
    plot_traces_heat_rate = []

    if length(index_orbit) == 2
        time = [config.solution.orientation.time[i] for i in index]
        append!(plot_traces_heat_rate, scatter(x=time, y=heat_rate, mode="lines", line=attr(color="black")))
    else
        time_end = 0
        for i in range(start=1, step=1, stop=length(index_orbit)-1)
            time = [config.solution.orientation.time[j] - time_0[i] + time_end for j in index[index_orbit[i]:index_orbit[i+1]-1]]
            heat_rate = [config.solution.performance.heat_rate[j] for j in index[index_orbit[i]:index_orbit[i+1]-1]]
            max_value = maximum(heat_rate)
            max_index = findfirst(x -> x == max_value, heat_rate)

            append!(plot_traces_heat_rate, scatter(x=time[max_index] , y=heat_rate[max_index], mode="markers", marker=attr(color="black")))
            time_end = time[end]
        end
    end

    layout_heat_rate = Layout(xaxis_title="Time [s]", yaxis_title=L"$\text{Heat rate [W/cm}^2 \text{}}]$")
    plot_heat_rate = plot(plot_traces_heat_rate, layout_heat_rate)

    plot_traces_heat_load = []

    if length(index_orbit) == 2
        time = [config.solution.orientation.time[i] for i in index]
        heat_load = [config.solution.performance.heat_load[i] for i in index]
        append!(plot_traces_heat_load, scatter(x=time, y=heat_load, mode="markers", marker=attr(color="black")))
    else
        time_end = 0

        for i in range(start=1, step=1, stop=length(index_orbit)-1)
            time = [config.solution.orientation.time[j] - time_0[i] + time_end for j in index[index_orbit[i]:index_orbit[i+1]-1]]
            heat_load = [config.solution.performance.heat_load[j] for j in index[index_orbit[i]:index_orbit[i+1]-1]]

            append!(plot_traces_heat_load, scatter(x=time[end] , y=heat_load[end], mode="markers", marker=attr(color="black")))
            time_end = time[end]
        end
    end

    layout_heat_load = Layout(xaxis_title="Time [s]", yaxis_title=L"$\text{Heat load [J/cm}^2 \text{}}]$")
    plot_heat_load = plot(plot_traces_heat_load, layout_heat_load)

    if args[:body_shape] == "Spacecraft"
        plot_traces_4 = scatter(x=[item/(60*60*24) for item in config.solution.orientation.time], y=[item - args[:dry_mass] for item in config.solution.performance.mass], mode="lines", line=attr(color="black"))
        layout_4 = Layout(xaxis_title="Time [s]", yaxis_title="Mass [kg]")
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

        if length(index_orbit) == 2
            time = [config.solution.orientation.time[i] for i in index]
            plot_traces_4 = scatter(x=time, y=q, mode="markers", marker=attr(color="black"))
            layout_4 = Layout(xaxis_title="Time [s]", yaxis_title="Dynamic pressure [Pa]")
        end
    end

    plot_4 = plot([plot_traces_4], layout_4)

    p = [plot_1 plot_heat_rate; plot_heat_load plot_4]
    relayout!(p, width=2100, height=1200)

    savefig(p, name * "_performance.png")
end

function traj_2D(state, m, name)
    x = [item*1e-3 for item in config.solution.orientation.pos_ii[1]]
    y = [item*1e-3 for item in config.solution.orientation.pos_ii[2]]
    z = [item*1e-3 for item in config.solution.orientation.pos_ii[3]]
    i = config.solution.orientation.oe[3][1]
    Ω = config.solution.orientation.oe[4][1]
    ω = config.solution.orientation.oe[5][1]

    T_ijk = [cos(Ω)*cos(ω)-sin(Ω)*sin(ω)*cos(i) sin(Ω)*cos(ω)+cos(Ω)*sin(ω)*cos(i) sin(ω)*sin(i);
             -cos(Ω)*sin(ω)-sin(Ω)*cos(ω)*cos(i) -sin(Ω)*sin(ω)+cos(Ω)*cos(ω)*cos(i) cos(ω)*sin(i);
             sin(Ω)*sin(i) -cos(Ω)*sin(i) cos(i)]
    
    vector = T_ijk .* [x; y; z]
    
    plot_traces = scatter(x=vector[1,:], y=vector[2,:], mode="lines", line=attr(color="black"))
    layout = Layout(xaxis_title="x [km]", yaxis_title="y [km]", shapes=[circle(xref="x", yref="y", fillcolor="OrangeRed", x0=0, y0=0, x1=0, y1=0, line_color="Red")])
    plot(plot_traces, layout)
    savefig(plot(plot_traces, layout), name * "_traj2D.png")
end

function traj_3D(state, m, name)

    p = make_subplots(rows=3, cols=3, 
                      specs=[Spec(kind="scene", rowspan=3, colspan=2) missing Spec(kind="scene")
                             missing missing Spec(kind="scene") 
                             missing missing Spec(kind="scene")])

    x = [item*1e-3 for item in config.solution.orientation.pos_ii[1]]
    y = [item*1e-3 for item in config.solution.orientation.pos_ii[2]]
    z = [item*1e-3 for item in config.solution.orientation.pos_ii[3]]

    r = m.planet.Rp_e*1e-3

    n = 100
    u = range(start=0, stop=2*pi, length=n)
    v = range(start=0, stop=pi, length=n)
    xs = r * cos.(u) * sin.(v)'
    ys = r * sin.(u) * sin.(v)'
    zs = r * ones(n) * cos.(v)'

    sphere1 = mesh3d(x=xs, y=ys, z=zs, color="rgba(168, 64, 50, 1.0)")

    r += 160                        # AE alt 160 km    
    xs = r * cos.(u) * sin.(v)'
    ys = r * sin.(u) * sin.(v)'
    zs = r * ones(n) * cos.(v)'

    sphere2 = mesh3d(x=xs, y=ys, z=zs, color="rgba(255, 255, 0, 0.25)")

    index = [1]
    for i in range(start=1, step=1, stop=length(config.solution.orientation.number_of_passage)-1)
        if config.solution.orientation.number_of_passage[i+1] - config.solution.orientation.number_of_passage[i] > 0
            append!(index, i)
        end
    end
    append!(length(config.solution.orientation.number_of_passage))

    for i in range(start=1, step=1, stop=length(index)-1)
        x_s = x[index[i]:index[i+1]]
        y_s = y[index[i]:index[i+1]]
        z_s = z[index[i]:index[i+1]]
        
        add_trace!(p, scatter3d(x=x_s, y=y_s, z=z_s, mode="lines", line=attr(color="black"), row=1, col=1))
        add_trace!(p, scatter3d(x=x_s, y=y_s, z=z_s, mode="lines", line=attr(color="black"), row=1, col=3))
        add_trace!(p, scatter3d(x=x_s, y=y_s, z=z_s, mode="lines", line=attr(color="black"), row=2, col=3))
        add_trace!(p, scatter3d(x=x_s, y=y_s, z=z_s, mode="lines", line=attr(color="black"), row=3, col=3))
    end

    relayout!(p, xaxis_title="x [km]", yaxis_title="y [km]", zaxis_title="z [km]")
    savefig(p, name * "_traj3D.png")

end

# function plot_visible(elev, azimuth, x, y, z)
#     a = azimuth*pi/180 - pi
#     e = elev*pi/180 - pi/2

#     X = [sin(e)*cos(a) sin(e)*sin(a) cos(e)]

#     Z = [x y z]

#     cond = (dot.(Z,X) >= 0)

#     x_ca = x[cond]
#     y_ca = y[cond]
#     z_ca = z[cond]

#     return x_ca, y_ca, z_ca
# end

# function plots_continuous(x, y, z, ax)
#     index_disc = [1]

#     for i in range(start=1, step=1, stop=length(x)-1)
#         if norm([x[i] - x[i-1], y[i] - y[i-1], z[i] - z[i-1]]) > 500
#             append!(index_disc, i-1)
#             append!(index_disc, i)
#         end
#     end
#     append!(index_disc, length(x))

#     for i in range(start=1, step=2, stop=length(index_disc)-1)
#         x_s = x[index_disc[i]:index_disc[i+1]]
#         y_s = y[index_disc[i]:index_disc[i+1]]
#         z_s = z[index_disc[i]:index_disc[i+1]]
#         plot_traces = 
#     end
# end

function ABM_periapsis(name)
    orbit_number = config.cnf.orbit_number_list
    periapsis_altitude = config.cnf.periapsis_list
    delta_v = [0]
    raise_man_orbit = [7,55,87,110,127,160,178,193,210,222]
    delta_v_raise = []
    lower_man_orbit = [13,25,30,36,48,68,72,80]
    delta_v_lower = []

    for j in range(start=1, step=1, stop=length(config.cnf.Δv_list)-1)
        append!(delta_v, config.cnf.Δv_list[j+1] - config.cnf.Δv_list[j])
    end

    for i in range(start=1, step=1, stop=length(delta_v))
        if j in raise_man_orbit
            append!(delta_v_raise, delta_v[i])
        elseif j in lower_man_orbit
            append!(delta_v_lower, delta_v[i])
        end
    end

    plot_traces_palt = scatter(x=orbit_number, y=periapsis_altitude, mode="lines", line=attr(color="black"))
    layout_palt = Layout(xaxis_title="Orbit number", yaxis_title="Periapsis altitude [km]")

    plot_traces_abm1 = scatter(x=[item + 19 for item in raise_man_orbit], delta_v_raise, mode="markers", marker=attr(color="black"), label="ABM to raise periapsis")
    plot_traces_abm2 = scatter(x=[item + 19 for item in lower_man_orbit], delta_v_lower, mode="markers", marker=attr(color="red"), label="ABM to lower periapsis")
    layout_abm = Layout(xaxis_title="Orbit number", yaxis_title="ABM Magnitude [m/s]")


end