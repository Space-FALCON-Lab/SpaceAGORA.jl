# ── Constants and helpers required by animate_all_satellites_3d_smooth_helper_target ──
using GeometryBasics: uv_normal_mesh, Tesselation, Sphere, Point3f
using FileIO: load

const R_EARTH  = 6.371e6  # m – mean Earth radius
const R_ATMDEF = 6.471e6  # m – top of nominal atmosphere (~100 km above surface)

@inline function rtn_basis(r, v)
    rhat = r / norm(r)
    nhat = cross(r, v)
    nhat = nhat / norm(nhat)
    that = cross(nhat, rhat)
    return rhat, that, nhat
end

"""
    Interactive 3D animation showing all N satellites in one plot using pure GLMakie.
    This version interpolates the solution to create a smooth animation at a fixed frame rate.

    Inputs:
        sol: ODE solution object
        p: parameter dictionary with keys:
            :N - number of satellites
            plus all keys required by laser_forces()
        helper_num: number of helper satellites (the last satellite is the target)
        tail: number of previous points to show in the trail (default 1000)
        Δt: time between frames in seconds (default 0.1s)
        show_earth: whether to show the Earth sphere (default true)
        earth_radius: radius of the Earth sphere in meters (default R_EARTH)
        markersize: size of the satellite markers (default 15)
        figure_size: size of the figure in pixels (default (1200, 800))
        trail_alpha: transparency of the trail lines (default 0.6)
        animation_fps: desired frames per second for smooth animation (default 30.0)

    Returns:
        nothing (displays the animation)

"""
function animate_all_satellites_3d_smooth_helper_target(sol, p, helper_num;
    tail = 1000,
    Δt = 0.1,
    show_earth = true,
    earth_radius = R_EARTH,
    markersize = 15,
    figure_size = (1200, 800),
    trail_alpha = 0.6,
    animation_fps = 30.0)

    println("Setting up 3D animation for all satellites...")

    # Determine N satellites
    N = haskey(p, :N) ? p[:N] : length(sol.u[1].sc)
    #println("Animating $N satellites")

    # Create regular time grid for smooth animation
    t_start = sol.t[1]
    t_end = sol.t[end]
    dt_animation = (t_end - t_start) / (animation_fps * 100)  # 10 seconds at 30fps = 300 points
    t_regular = t_start:dt_animation:t_end
    n_points = length(t_regular)

    # Interpolate solution at regular intervals
    positions = zeros(3, N, n_points)  # [xyz, sat_id, time]
    velocities = zeros(3, N, n_points)  # [vxyz, sat_id, time]
    for (k, t) in enumerate(t_regular)
        u_interp = sol(t)  # This interpolates automatically
        for i in 1:N
            sc = u_interp.sc[i]
            positions[1, i, k] = sc.pos[1]  # X in m
            positions[2, i, k] = sc.pos[2]  # Y in m
            positions[3, i, k] = sc.pos[3]  # Z in m
            velocities[1, i, k] = sc.vel[1]  # Vx in m/s
            velocities[2, i, k] = sc.vel[2]  # Vy in m/s
            velocities[3, i, k] = sc.vel[3]  # Vz in m/s
        end
    end

    # Set up GLMakie figure
    #GLMakie.activate!()
    fig = GLMakie.Figure(size = figure_size)
    ax = GLMakie.Axis3(fig[1, 1],
               xlabel = "X [m]", 
               ylabel = "Y [m]", 
               zlabel = "Z [m]",
               title = "Multi-Satellite Orbital Animation (Interpolated)",
               aspect = (1, 1, 1))  # Ensure equal aspect ratio

    # # Colors for each satellite
    # basic_colors = [:red, :blue, :green, :orange, :purple, :cyan, :magenta, :yellow, :brown, :pink]
    # colors = basic_colors[1:min(N, length(basic_colors))]
    # if N > length(basic_colors)
    #     # Repeat colors if more satellites than colors
    #     colors = [colors; basic_colors[1:(N - length(basic_colors))]]
    # end

    # Add Earth sphere with texture
    # method from https://beautiful.makie.org/dev/examples/3d/meshes/Earth_planes
    if show_earth
        # Create a tessellated sphere for Earth
        earth_sphere = uv_normal_mesh(Tesselation(Sphere(Point3f(0.0, 0.0, 0.0), earth_radius), 64))
        earth_img_path = joinpath(@__DIR__, "8k_earth_daymap.jpg")
        if isfile(earth_img_path)
            earth_img = load(earth_img_path)
            GLMakie.mesh!(ax, earth_sphere;
                color = circshift(earth_img, (0, 3000)),
                ssao = true,
                alpha = 0.5,
                transparency = true)
        else
            # Fallback: solid blue sphere when texture file is unavailable
            GLMakie.mesh!(ax, earth_sphere;
                color = :steelblue,
                alpha = 0.5,
                transparency = true)
        end
    end

    # Create observables for animation
    current_time = GLMakie.Observable(1)

    # Pre-compute axis limits
    all_pos = reshape(positions, 3, :)
    xlims = (minimum(all_pos[1, :]), maximum(all_pos[1, :]))
    ylims = (minimum(all_pos[2, :]), maximum(all_pos[2, :]))
    zlims = (minimum(all_pos[3, :]), maximum(all_pos[3, :]))

    # Add padding and ensure non-zero range
    pad = 0.1
    min_range = 1500e3  # 1500 km minimum range
    xlims = xlims[1] ≈ xlims[2] ? (xlims[1] - min_range, xlims[1] + min_range) : (xlims[1] * (1 + pad), xlims[2] * (1 + pad))
    ylims = ylims[1] ≈ ylims[2] ? (ylims[1] - min_range, ylims[1] + min_range) : (ylims[1] * (1 + pad), ylims[2] * (1 + pad))
    zlims = zlims[1] ≈ zlims[2] ? (zlims[1] - min_range, zlims[1] + min_range) : (zlims[1] * (1 + pad), zlims[2] * (1 + pad))

    # Calculate the maximum absolute range across all dimensions
    max_limit = maximum([xlims[2]; ylims[2]; zlims[2]])

    # Set symmetric limits for all axes
    xlims = (-max_limit, max_limit)
    ylims = (-max_limit, max_limit)
    zlims = (-max_limit, max_limit)

    # Apply the limits to the axis
    GLMakie.limits!(ax, xlims, ylims, zlims)

    # Colors and labels for satellites
    helper_color = :blue
    target_color = :red

    for i in 1:N
        # Determine if the satellite is a helper or a target
        if i <= helper_num
            # Helper satellites
            trail_label =  "Helper Sat"  # Add label only for the first helper satellite
            current_label =  "Helper Sat"
            color = helper_color
        else
            # Target satellites
            trail_label = "Target Sat" # Add label only for the first target satellite
            current_label = "Target Sat"
            color = target_color
        end

        # Trail observable
        trail_points = GLMakie.@lift begin
            t_idx = $current_time
            start_idx = max(1, t_idx - tail)
            trail_x = positions[1, i, start_idx:t_idx]
            trail_y = positions[2, i, start_idx:t_idx]
            trail_z = positions[3, i, start_idx:t_idx]
            GLMakie.Point3f.(trail_x, trail_y, trail_z)
        end

        # Current position observable
        current_pos = GLMakie.@lift begin
            t_idx = $current_time
            [GLMakie.Point3f(positions[1, i, t_idx], positions[2, i, t_idx], positions[3, i, t_idx])]
        end

        # Plot trail
        GLMakie.lines!(ax, trail_points,
            color = color,
            alpha = trail_alpha,
            linewidth = 2,
            label = trail_label)

        # Plot current position
        GLMakie.scatter!(ax, current_pos,
            color = color,
            markersize = markersize,
            strokewidth = 1,
            strokecolor = :black,
            label = current_label)
    end


    # Add RTN arrows for the target satellite
    target_sat = N  # Assuming the last satellite is the target

   # one line segment (er): exactly two points
    er_line = GLMakie.@lift begin
        t = $current_time
        r = @SVector [positions[1,N,t], positions[2,N,t], positions[3,N,t]]
        v = @SVector [velocities[1,N,t], velocities[2,N,t], velocities[3,N,t]]
        er, _, _ = rtn_basis(r, v)
        s = 1e6
        GLMakie.Point3f.(
            [r[1], r[1] + s*er[1]],
            [r[2], r[2] + s*er[2]],
            [r[3], r[3] + s*er[3]],
        )
    end
    GLMakie.lines!(ax, er_line, color = :red, linewidth = 2)

    # one line segment (et): exactly two points
    et_line = GLMakie.@lift begin
        t = $current_time
        r = @SVector [positions[1,N,t], positions[2,N,t], positions[3,N,t]]
        v = @SVector [velocities[1,N,t], velocities[2,N,t], velocities[3,N,t]]
        _, et, _ = rtn_basis(r, v)
        s = 1e6
        GLMakie.Point3f.(
            [r[1], r[1] + s*et[1]],
            [r[2], r[2] + s*et[2]],
            [r[3], r[3] + s*et[3]],
        )
    end
    GLMakie.lines!(ax, et_line, color = :green, linewidth = 2)

    # one line segment (en): exactly two points
    en_line = GLMakie.@lift begin
        t = $current_time
        r = @SVector [positions[1,N,t], positions[2,N,t], positions[3,N,t]]
        v = @SVector [velocities[1,N,t], velocities[2,N,t], velocities[3,N,t]]
        _, _, en = rtn_basis(r, v)
        s = 1e6
        GLMakie.Point3f.(
            [r[1], r[1] + s*en[1]],
            [r[2], r[2] + s*en[2]],
            [r[3], r[3] + s*en[3]],
        )
    end
    GLMakie.lines!(ax, en_line, color = :blue, linewidth = 2)

    # Record satellite pair points dynamically with LOS status
    sat_pair_points = Dict{Tuple{Int, Int}, Vector{Any}}()
    if haskey(p, :Pmatrix) && haskey(p, :cavity)
        #println("Recording satellite pair points and LOS status for single-pass and open-cavity pairs...")
        Pm = p[:Pmatrix]  # Matrix indicating single-pass pairs
        cavity = p[:cavity]  # Dictionary indicating open-cavity satellites

        # Read LOS parameters
        use_los = get(p, :use_los, false)  # Default to false if :use_los is not provided
        R_atm = get(p, :R_atm, R_ATMDEF)  # Atmosphere radius
        atm_cl = get(p, :atm_clearance, 0.0)  # Minimum clearance above atmosphere
        minR    = get(p, :min_range, 0.0)
        maxR    = get(p, :max_range, Inf)

    #     # current_helpers as an Observable (N×N zeros)
    #     current_helpers = GLMakie.Observable(fill(0, (N, N)))

    #     GLMakie.@lift begin
    #         _ = $current_time  # make reactive
    #         current_helpers[] = fill(0, (N, N))  # Reset after initialization
    #         nothing
    #     end

    #     # Iterate over all satellite pairs
    #     for i in 1:N
    #         for j in i+1:N
    #             # Observable for satellite pair points
    #             pair_points = GLMakie.@lift begin
    #                 t_idx = $current_time
    #                 pos_i = GLMakie.Point3f(positions[1, i, t_idx], positions[2, i, t_idx], positions[3, i, t_idx])
    #                 pos_j = GLMakie.Point3f(positions[1, j, t_idx], positions[2, j, t_idx], positions[3, j, t_idx]) 
    #                 [pos_i, pos_j]
    #             end
                
    #             if use_los
    #                 pair_los_transparency = GLMakie.@lift begin
    #                     t_idx = $current_time  
                        
    #                     ri = @SVector [positions[1, i, t_idx], positions[2, i, t_idx], positions[3, i, t_idx]]
    #                     rj = @SVector [positions[1, j, t_idx], positions[2, j, t_idx], positions[3, j, t_idx]]
    #                     met = los_metrics(ri, rj; R_atm=R_atm)
    #                     los_ok = (!use_los) || ((!met.blocked) && (met.clearance >= atm_cl))
    #                     range_ok = (met.slant_range >= minR) && (met.slant_range <= maxR)
    #                     if los_ok && range_ok
    #                         if reserve_link(i, j, $current_helpers)
    #                             1 # Line of sight is clear and within range
    #                         else
    #                             0
    #                         end
    #                     else
    #                         0 # Line of sight is blocked or out of range
    #                     end
    #                 end
    #             else
    #                 pair_los_transparency = GLMakie.@lift begin
    #                     t_idx = $current_time  
    #                     ri = @SVector [positions[1, i, t_idx], positions[2, i, t_idx], positions[3, i, t_idx]]
    #                     rj = @SVector [positions[1, j, t_idx], positions[2, j, t_idx], positions[3, j, t_idx]]
    #                     met = los_metrics(ri, rj; R_atm=R_atm)
    #                     range_ok = (met.slant_range >= minR) && (met.slant_range <= maxR)
    #                     if range_ok
    #                         if reserve_link(i, j, $current_helpers)
    #                             1 # Line of sight is clear and within range
    #                         else
    #                             0
    #                         end
    #                     else
    #                         0 # Out of range
    #                     end
    #                 end
    #             end

    #             # Check if satellites i and j are single-pass pairs
    #             if Pm[i, j] != 0
    #                 GLMakie.lines!(ax, pair_points, color = :green, alpha = pair_los_transparency, linewidth = 2, label = "Single-pass pair")
    #                 GLMakie.@lift begin
    #                     _ = $current_time  # make reactive
    #                     H = copy(current_helpers[])
    #                     H[i, :] .= 1 # 1 means single_pass # i is busy and not available to be helper
    #                     H[:, i] .= 1 # 1 means single_pass # i is busy and not available to be target
    #                     H[j, :] .= 1 # 1 means single_pass # j is busy and not available to be helper
    #                     H[:, j] .= 1 # 1 means single_pass # j is busy and not available to be target
    #                     current_helpers[] = H
    #                     nothing
    #                 end
    #             end

    #             # Check if satellites i and j are open-cavity pairs
    #             if haskey(cavity, (i, j))
    #                 GLMakie.lines!(ax, pair_points, color = :orange, alpha = pair_los_transparency, linewidth = 2, label = "Open-cavity pair")
    #                 GLMakie.@lift begin
    #                     _ = $current_time  # make reactive
    #                     if pair_los_transparency == 1
    #                         H = copy(current_helpers[])
    #                         H[i, :] .= 2 # 2 means open_cavity # i is busy and not available to be helper
    #                         H[:, i] .= 2 # 2 means open_cavity # i is busy and not available to be target
    #                         H[j, :] .= 2 # 2 means open_cavity # j is busy and not available to be helper
    #                         H[:, j] .= 2 # 2 means open_cavity # j is busy and not available to be target
    #                         current_helpers[] = H
    #                     else
    #                         H = copy(current_helpers[])
    #                         H[i, :] .= 0 # 0 means not busy # i is not busy and available to be helper
    #                         H[:, i] .= 0 # 0 means not busy # i is not busy and available to be target
    #                         H[j, :] .= 0 # 0 means not busy # j is not busy and available to be helper
    #                         H[:, j] .= 0 # 0 means not busy # j is not busy and available to be target
    #                         current_helpers[] = H
    #                     nothing
    #                     end
    #                 end
    #             end

    #         end
    #     end

    # end

    
        # current_helpers as an Observable (N×N zeros)
        current_helpers = GLMakie.Observable(fill(0, (N, N)))

        # Compute all link alphas and reservations once per frame (deterministic)
        frame_links = GLMakie.@lift begin
            t_idx = $current_time
            # initialize fresh busy matrix for this frame
            H = fill(0, (N, N)) # H means helper status: 0=free, 1=single_pass, 2=open_cavity
            # per-pair alpha (Float32)
            A = fill(0f0, (N, N)) # A means transparency alpha

            for i in 1:N
                for j in i+1:N
                    # LOS / range
                    ri = @SVector [positions[1, i, t_idx], positions[2, i, t_idx], positions[3, i, t_idx]]
                    rj = @SVector [positions[1, j, t_idx], positions[2, j, t_idx], positions[3, j, t_idx]]
                    met = los_metrics(ri, rj; R_atm = R_atm)
                    los_ok = use_los ? ((!met.blocked) && (met.clearance >= atm_cl)) : true
                    range_ok = (met.slant_range >= minR) && (met.slant_range <= maxR)

                    if los_ok && range_ok
                        # enforce exclusivity: i and j must be free
                        i_free = all(H[i, :] .== 0) && all(H[:, i] .== 0)
                        j_free = all(H[j, :] .== 0) && all(H[:, j] .== 0)
                        if i_free && j_free
                            # decide kind (single-pass vs cavity); if both true, pick one (e.g., prefer single-pass)
                            kind = 0
                            if Pm[i, j] != 0
                                kind = 1
                            elseif haskey(cavity, (i, j))
                                kind = 2
                            end

                            if kind > 0
                                # reserve i and j
                                H[i, :] .= kind
                                H[:, i] .= kind
                                H[j, :] .= kind
                                H[:, j] .= kind
                                A[i, j] = 1f0
                                A[j, i] = 1f0
                            end
                        end
                    end
                end
            end

            # publish the busy matrix for the frame
            current_helpers[] = H
            A # return the transparency A to frame_links
        end

        # Iterate over all satellite pairs (draw using precomputed alpha)
        for i in 1:N
            for j in i+1:N
                pair_points = GLMakie.@lift begin
                    t = $current_time
                    pos_i = GLMakie.Point3f(positions[1, i, t], positions[2, i, t], positions[3, i, t])
                    pos_j = GLMakie.Point3f(positions[1, j, t], positions[2, j, t], positions[3, j, t])
                    [pos_i, pos_j]
                end

                # alpha for this pair is just frame_links[i,j]
                pair_alpha = GLMakie.@lift $frame_links[i, j]

                if Pm[i, j] != 0
                    GLMakie.lines!(ax, pair_points, color = :green, alpha = pair_alpha, linewidth = 2, label = "Single-pass pair")
                end
                if haskey(cavity, (i, j))
                    GLMakie.lines!(ax, pair_points, color = :orange, alpha = pair_alpha, linewidth = 2, label = "Open-cavity pair")
                end
            end
        end
    end
    # Add legend
    GLMakie.axislegend(ax, position = :lt, unique = true)

    # Animation controls
    controls_grid = fig[2, 1] = GLMakie.GridLayout()

    # Play/Pause button
    is_playing = GLMakie.Observable(false)
    play_button = GLMakie.Button(controls_grid[1, 1], label = GLMakie.@lift($is_playing ? "|| Pause" : "=> Play"))

    # Speed slider
    GLMakie.Label(controls_grid[1, 2], "Speed:")
    speed_slider = GLMakie.Slider(controls_grid[1, 3], range = 0.1:0.1:5.0, value = 1.0)

    # Time slider
    GLMakie.Label(controls_grid[1, 5], "Time:")
    time_slider = GLMakie.Slider(controls_grid[1, 6], range = 1:length(t_regular), value = 1)

    # Reset button
    reset_button = GLMakie.Button(controls_grid[1, 7], label = "<<= Reset")

    # Time display
    # time_text = GLMakie.@lift("t = $(round(t_regular[$current_time], digits=1)) s")
    # GLMakie.Label(fig[0, 1], time_text, tellwidth = false, fontsize = 16)
        # Time display
    time_text = GLMakie.@lift begin
        t_idx = $current_time
        target_distance = norm([positions[1, N, t_idx], positions[2, N, t_idx], positions[3, N, t_idx]])
        "t = $(round(t_regular[t_idx], digits=1)) s, Target Distance = $(round(target_distance / 1e3, digits=1)) km"
    end
    GLMakie.Label(fig[0, 1], time_text, tellwidth = false, fontsize = 16)


    # Button callbacks
    GLMakie.on(play_button.clicks) do _
        is_playing[] = !is_playing[]
    end

    GLMakie.on(time_slider.value) do val
        current_time[] = Int(val)  # Ensure integer value
    end

    GLMakie.on(reset_button.clicks) do _
        current_time[] = 1
        GLMakie.set_close_to!(time_slider, 1)
        is_playing[] = false
    end

    # Animation loop with consistent timing
    animation_task = @async begin
        target_dt = 1.0 / animation_fps
        while true
            if is_playing[]
                if current_time[] < length(t_regular)
                    speed_val = GLMakie.to_value(speed_slider.value)
                    sleep(target_dt / max(speed_val, 0.1))

                    new_time = min(current_time[] + 1, length(t_regular))
                    current_time[] = new_time
                    GLMakie.set_close_to!(time_slider, new_time)
                else
                    is_playing[] = false
                end
            else
                sleep(0.1)
            end
        end
    end

    screen = GLMakie.display(fig)

    return fig, (current_time = current_time,
                is_playing = is_playing,
                animation_task = animation_task,
                screen = screen)
end