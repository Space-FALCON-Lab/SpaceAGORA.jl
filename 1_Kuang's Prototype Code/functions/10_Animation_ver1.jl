"""
    Interactive 3D animation showing all N satellites in one plot using pure GLMakie.
    This bypasses DynamicalSystems' 3D limitation.

    This function animates the status of satellites at each time step. Hence, due to adaptive step length in ODE solver,
    there may be some jitter in the animation speed.

    Inputs:
        sol: ODE solution object
        p: parameter dictionary with keys:
            :N - number of satellites
            plus all keys required by laser_forces()
        tail: number of previous points to show in the trail (default 1000)
        Δt: time between frames in seconds (default 0.1s)
        show_earth: whether to show the Earth sphere (default true)
        earth_radius: radius of the Earth sphere in meters (default R_EARTH)
        markersize: size of the satellite markers (default 15)
        figure_size: size of the figure in pixels (default (1200, 800))
        trail_alpha: transparency of the trail lines (default 0.6)

    Returns:
        nothing (displays the animation)
"""
function animate_all_satellites_3d(sol, p;
    tail = 1000,
    Δt = 0.1,
    show_earth = true,
    earth_radius = R_EARTH,
    markersize = 15,
    figure_size = (1200, 800),
    trail_alpha = 0.6)

    println("Setting up 3D animation for all satellites...")
    
    # Determine N satellites
    N = haskey(p, :N) ? p[:N] : (length(sol.u[1]) ÷ 6)
    println("Animating $N satellites")
    
    # Extract all satellite trajectories
    n_points = length(sol.t)
    positions = zeros(3, N, n_points)  # [xyz, sat_id, time]
    
    for (k, t) in enumerate(sol.t)
        u = sol.u[k]
        for i in 1:N
            positions[1, i, k] = u[idx(i,1)] / 1e3  # X in km
            positions[2, i, k] = u[idx(i,2)] / 1e3  # Y in km
            positions[3, i, k] = u[idx(i,3)] / 1e3  # Z in km
        end
    end
    
    # Set up GLMakie figure
    #GLMakie.activate!()
    fig = GLMakie.Figure(size = figure_size)  # Added GLMakie.
    ax = GLMakie.Axis3(fig[1, 1],  # Added GLMakie.
               xlabel = "X [km]", 
               ylabel = "Y [km]", 
               zlabel = "Z [km]",
               title = "Multi-Satellite Orbital Animation")
    
    # Colors for each satellite
    basic_colors = [:red, :blue, :green, :orange, :purple, :cyan, :magenta, :yellow, :brown, :pink]
    colors = basic_colors[1:min(N, length(basic_colors))]
    if N > length(basic_colors)
        # Repeat colors if more satellites than colors
        colors = [colors; basic_colors[1:(N-length(basic_colors))]]
    end
    
    # # Add Earth sphere
    # if show_earth
    #     θ = range(0, π, length=30)
    #     ϕ = range(0, 2π, length=60)
    #     x_earth = [earth_radius*sin(t)*cos(p) for t in θ, p in ϕ] ./ 1e3
    #     y_earth = [earth_radius*sin(t)*sin(p) for t in θ, p in ϕ] ./ 1e3
    #     z_earth = [earth_radius*cos(t) for t in θ, p in ϕ] ./ 1e3
        
    #     GLMakie.surface!(ax, x_earth, y_earth, z_earth,
    #             color = :lightblue, alpha = 0.3, shading = true)
    # end
    
    # Create observables for animation
    current_time = GLMakie.Observable(1)  # Added GLMakie.
    
    # Pre-compute axis limits
    all_pos = reshape(positions, 3, :)
    xlims = (minimum(all_pos[1,:]), maximum(all_pos[1,:]))
    ylims = (minimum(all_pos[2,:]), maximum(all_pos[2,:]))
    zlims = (minimum(all_pos[3,:]), maximum(all_pos[3,:]))
    
    # Add padding and ensure non-zero range
    pad = 0.1
    min_range = 1000.0  # 1000 km minimum range
    
    xlims = xlims[1] ≈ xlims[2] ? (xlims[1] - min_range, xlims[1] + min_range) : (xlims[1]*(1+pad), xlims[2]*(1+pad))
    ylims = ylims[1] ≈ ylims[2] ? (ylims[1] - min_range, ylims[1] + min_range) : (ylims[1]*(1+pad), ylims[2]*(1+pad))
    zlims = zlims[1] ≈ zlims[2] ? (zlims[1] - min_range, zlims[1] + min_range) : (zlims[1]*(1+pad), zlims[2]*(1+pad))
    
    GLMakie.limits!(ax, xlims, ylims, zlims)
    
    # Plot trails and current positions for each satellite
    for i in 1:N
        # Trail observable
        trail_points = GLMakie.@lift begin  # Added GLMakie.
            t_idx = $current_time
            start_idx = max(1, t_idx - tail)
            trail_x = positions[1, i, start_idx:t_idx]
            trail_y = positions[2, i, start_idx:t_idx]
            trail_z = positions[3, i, start_idx:t_idx]
            GLMakie.Point3f.(trail_x, trail_y, trail_z)  # Added GLMakie.
        end
        
        # Current position observable
        current_pos = GLMakie.@lift begin  # Added GLMakie.
            t_idx = $current_time
            [GLMakie.Point3f(positions[1, i, t_idx], positions[2, i, t_idx], positions[3, i, t_idx])]  # Added GLMakie.
        end
        
        # Plot trail
        GLMakie.lines!(ax, trail_points, 
               color = colors[i], 
               alpha = trail_alpha,
               linewidth = 2,
               label = "Sat $i trail")
        
        # Plot current position
        GLMakie.scatter!(ax, current_pos,
                color = colors[i],
                markersize = markersize,
                strokewidth = 1,
                strokecolor = :black,
                label = "Sat $i")
    end
    
    # Add legend
    GLMakie.axislegend(ax, position = :lt, unique = true)
    
    # Time display
    time_text = GLMakie.@lift("t = $(round(sol.t[$current_time], digits=1)) s")  # Added GLMakie.
    GLMakie.Label(fig[0, 1], time_text, tellwidth = false, fontsize = 16)
    
    # Animation controls
    controls_grid = fig[2, 1] = GLMakie.GridLayout()
    
    # Play/Pause button
    is_playing = GLMakie.Observable(false)  # Added GLMakie.
    play_button = GLMakie.Button(controls_grid[1, 1], label = GLMakie.@lift($is_playing ? "⏸️ Pause" : "▶️ Play"))  # Added GLMakie.
    
    # Speed slider - FIX the label Observable
    GLMakie.Label(controls_grid[1, 2], "Speed:")
    speed_slider = GLMakie.Slider(controls_grid[1, 3], range = 0.1:0.1:5.0, value = 1.0)
    # Fix this line - create a proper Observable for the speed label
    # speed_label = GLMakie.Label(controls_grid[1, 4], GLMakie.@lift("$(round(GLMakie.to_value($speed_slider.value), digits=1))x"))
    
    # Time slider
    GLMakie.Label(controls_grid[1, 5], "Time:")
    time_slider = GLMakie.Slider(controls_grid[1, 6], range = 1:length(sol.t), value = 1)
    
    # Reset button
    reset_button = GLMakie.Button(controls_grid[1, 7], label = "⏮️ Reset")
    
    # Button callbacks - fix slider access
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
    
    # Animation loop - corrected version
    animation_task = @async begin
        while true
            if is_playing[]
                if current_time[] < length(sol.t)
                    speed_val = GLMakie.to_value(speed_slider.value)
                    sleep(Δt / max(speed_val, 0.1))
                    current_time[] = min(current_time[] + 1, length(sol.t))
                    GLMakie.set_close_to!(time_slider, current_time[])
                else
                    is_playing[] = false
                end
            else
                sleep(0.1)
            end
        end
    end
    
    GLMakie.display(fig)  # Added GLMakie.
    
    return fig, (current_time = current_time, 
                is_playing = is_playing,
                animation_task = animation_task)
end

"""
    Interactive 3D animation showing all N satellites in one plot using pure GLMakie.
    This version interpolates the solution to create a smooth animation at a fixed frame rate.

    Inputs:
        sol: ODE solution object
        p: parameter dictionary with keys:
            :N - number of satellites
            plus all keys required by laser_forces()
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
function animate_all_satellites_3d_smooth(sol, p;
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
    N = haskey(p, :N) ? p[:N] : (length(sol.u[1]) ÷ 6)
    println("Animating $N satellites")

    # Create regular time grid for smooth animation
    t_start = sol.t[1]
    t_end = sol.t[end]
    dt_animation = (t_end - t_start) / (animation_fps * 10)  # 10 seconds at 30fps = 300 points
    t_regular = t_start:dt_animation:t_end
    n_points = length(t_regular)

    println("  Original time points: ", length(sol.t))
    println("  Interpolated points: ", n_points)
    println("  Animation dt: ", dt_animation, " seconds")

    # Interpolate solution at regular intervals
    positions = zeros(3, N, n_points)  # [xyz, sat_id, time]
    for (k, t) in enumerate(t_regular)
        u_interp = sol(t)  # This interpolates automatically
        for i in 1:N
            positions[1, i, k] = u_interp[idx(i, 1)]  # X in m
            positions[2, i, k] = u_interp[idx(i, 2)]  # Y in m
            positions[3, i, k] = u_interp[idx(i, 3)]  # Z in m
        end
    end

    # Set up GLMakie figure
    GLMakie.activate!()
    fig = GLMakie.Figure(size = figure_size)
    ax = GLMakie.Axis3(fig[1, 1],
               xlabel = "X [m]", 
               ylabel = "Y [m]", 
               zlabel = "Z [m]",
               title = "Multi-Satellite Orbital Animation (Interpolated)",
               aspect = (1, 1, 1))  # Ensure equal aspect ratio

    # Colors for each satellite
    basic_colors = [:red, :blue, :green, :orange, :purple, :cyan, :magenta, :yellow, :brown, :pink]
    colors = basic_colors[1:min(N, length(basic_colors))]
    if N > length(basic_colors)
        # Repeat colors if more satellites than colors
        colors = [colors; basic_colors[1:(N - length(basic_colors))]]
    end

    # Add Earth sphere with texture
    # method from https://beautiful.makie.org/dev/examples/3d/meshes/Earth_planes
    if show_earth
        println("Adding Earth sphere with texture...")
        # Download and load the Earth texture
        earth_img_path = "Kuang's Prototype Code/input/8k_earth_daymap.jpg"
        earth_img = load(earth_img_path)  # Load the image as an array

        # Create a tessellated sphere for Earth
        earth_sphere = uv_normal_mesh(Tesselation(Sphere(Point3f(0.0, 0.0, 0.0), earth_radius), 64))

        # Plot the Earth sphere with the texture
        GLMakie.mesh!(ax, earth_sphere; 
                color = circshift(earth_img, (0, 3000)), 
                ssao = true, 
                alpha = 0.5, 
                transparency = true)  # Set transparency level (0 = fully transparent, 1 = fully opaque)
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

    # Plot trails and current positions for each satellite
    for i in 1:N
        # Trail observable
        trail_points = GLMakie.@lift begin #@lift means this block is reactive to changes in current_time
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
               color = colors[i], 
               alpha = trail_alpha,
               linewidth = 2,
               label = "Sat $i trail")

        # Plot current position
        GLMakie.scatter!(ax, current_pos,
                color = colors[i],
                markersize = markersize,
                strokewidth = 1,
                strokecolor = :black,
                label = "Sat $i")
    end

    # Record satellite pair points dynamically with LOS status
    sat_pair_points = Dict{Tuple{Int, Int}, Vector{Any}}()
    if haskey(p, :Pmatrix) && haskey(p, :cavity)
        println("Recording satellite pair points and LOS status for single-pass and open-cavity pairs...")
        Pm = p[:Pmatrix]  # Matrix indicating single-pass pairs
        cavity = p[:cavity]  # Dictionary indicating open-cavity satellites

        # Read LOS parameters
        use_los = get(p, :use_los, false)  # Default to false if :use_los is not provided
        R_atm = get(p, :R_atm, R_ATMDEF)  # Atmosphere radius
        atm_cl = get(p, :atm_clearance, 0.0)  # Minimum clearance above atmosphere
        minR    = get(p, :min_range, 0.0)
        maxR    = get(p, :max_range, Inf)
        
        # Iterate over all satellite pairs
        for i in 1:N
            for j in i+1:N
                # Observable for satellite pair points
                pair_points = GLMakie.@lift begin
                    t_idx = $current_time
                    pos_i = GLMakie.Point3f(positions[1, i, t_idx], positions[2, i, t_idx], positions[3, i, t_idx])
                    pos_j = GLMakie.Point3f(positions[1, j, t_idx], positions[2, j, t_idx], positions[3, j, t_idx])
                    
                    [pos_i, pos_j]
                end
                
                if use_los
                    pair_los_transparency = GLMakie.@lift begin
                        t_idx = $current_time  
                        ri = @SVector [positions[1, i, t_idx], positions[2, i, t_idx], positions[3, i, t_idx]]
                        rj = @SVector [positions[1, j, t_idx], positions[2, j, t_idx], positions[3, j, t_idx]]
                        met = los_metrics(ri, rj; R_atm=R_atm)
                        los_ok = (!use_los) || ((!met.blocked) && (met.clearance >= atm_cl))
                        range_ok = (met.slant_range >= minR) && (met.slant_range <= maxR)
                        if los_ok && range_ok
                            1 # Line of sight is clear and within range
                        else
                            0 # Line of sight is blocked or out of range
                        end
                    end
                else
                    pair_los_transparency = GLMakie.@lift begin
                        t_idx = $current_time  
                        ri = @SVector [positions[1, i, t_idx], positions[2, i, t_idx], positions[3, i, t_idx]]
                        rj = @SVector [positions[1, j, t_idx], positions[2, j, t_idx], positions[3, j, t_idx]]
                        met = los_metrics(ri, rj; R_atm=R_atm)
                        range_ok = (met.slant_range >= minR) && (met.slant_range <= maxR)
                        if range_ok
                            1 # Within range
                        else
                            0 # Out of range
                        end
                    end
                end
                
                # print pair_los type
                #println("Pair ($i, $j) LOS type: ", typeof(pair_los), " value: ", pair_los[])

                # Check if satellites i and j are single-pass pairs
                
                if Pm[i, j] != 0
                    GLMakie.lines!(ax, pair_points, color = :blue, alpha = pair_los_transparency, linewidth = 2, label = "Single-pass pair")
                end

                # Check if satellites i and j are open-cavity pairs
                if haskey(cavity, (i, j))
                    GLMakie.lines!(ax, pair_points, color = :green, alpha = pair_los_transparency, linewidth = 2, label = "Open-cavity pair")
                end

            end
        end
    end

    # Add legend
    GLMakie.axislegend(ax, position = :lt)

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
    time_text = GLMakie.@lift("t = $(round(t_regular[$current_time], digits=1)) s")
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

    GLMakie.display(fig)

    return fig, (current_time = current_time, 
                is_playing = is_playing,
                animation_task = animation_task)
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
    N = haskey(p, :N) ? p[:N] : (length(sol.u[1]) ÷ 6)
    #println("Animating $N satellites")

    # Create regular time grid for smooth animation
    t_start = sol.t[1]
    t_end = sol.t[end]
    dt_animation = (t_end - t_start) / (animation_fps * 10)  # 10 seconds at 30fps = 300 points
    t_regular = t_start:dt_animation:t_end
    n_points = length(t_regular)

    # Interpolate solution at regular intervals
    positions = zeros(3, N, n_points)  # [xyz, sat_id, time]
    velocities = zeros(3, N, n_points)  # [vxyz, sat_id, time]
    for (k, t) in enumerate(t_regular)
        u_interp = sol(t)  # This interpolates automatically
        for i in 1:N
            positions[1, i, k] = u_interp[idx(i, 1)]  # X in m
            positions[2, i, k] = u_interp[idx(i, 2)]  # Y in m
            positions[3, i, k] = u_interp[idx(i, 3)]  # Z in m
            velocities[1, i, k] = u_interp[idx(i, 4)]  # Vx in m/s
            velocities[2, i, k] = u_interp[idx(i, 5)]  # Vy in m/s
            velocities[3, i, k] = u_interp[idx(i, 6)]  # Vz in m/s
        end
    end

    # Set up GLMakie figure
    GLMakie.activate!()
    fig = GLMakie.Figure(size = figure_size)
    ax = GLMakie.Axis3(fig[1, 1],
               xlabel = "X [m]", 
               ylabel = "Y [m]", 
               zlabel = "Z [m]",
               title = "Multi-Satellite Orbital Animation (Interpolated)",
               aspect = (1, 1, 1))  # Ensure equal aspect ratio

    # Colors for each satellite
    basic_colors = [:red, :blue, :green, :orange, :purple, :cyan, :magenta, :yellow, :brown, :pink]
    colors = basic_colors[1:min(N, length(basic_colors))]
    if N > length(basic_colors)
        # Repeat colors if more satellites than colors
        colors = [colors; basic_colors[1:(N - length(basic_colors))]]
    end

    # Add Earth sphere with texture
    # method from https://beautiful.makie.org/dev/examples/3d/meshes/Earth_planes
    if show_earth
        #println("Adding Earth sphere with texture...")
        # Download and load the Earth texture
        earth_img_path = "Kuang's Prototype Code/input/8k_earth_daymap.jpg"
        earth_img = load(earth_img_path)  # Load the image as an array

        # Create a tessellated sphere for Earth
        earth_sphere = uv_normal_mesh(Tesselation(Sphere(Point3f(0.0, 0.0, 0.0), earth_radius), 64))

        # Plot the Earth sphere with the texture
        GLMakie.mesh!(ax, earth_sphere; 
                color = circshift(earth_img, (0, 3000)), 
                ssao = true, 
                alpha = 0.5, 
                transparency = true)  # Set transparency level (0 = fully transparent, 1 = fully opaque)
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
        
        # Iterate over all satellite pairs
        for i in 1:N
            for j in i+1:N
                # Observable for satellite pair points
                pair_points = GLMakie.@lift begin
                    t_idx = $current_time
                    pos_i = GLMakie.Point3f(positions[1, i, t_idx], positions[2, i, t_idx], positions[3, i, t_idx])
                    pos_j = GLMakie.Point3f(positions[1, j, t_idx], positions[2, j, t_idx], positions[3, j, t_idx])
                    
                    [pos_i, pos_j]
                end
                
                if use_los
                    pair_los_transparency = GLMakie.@lift begin
                        t_idx = $current_time  
                        ri = @SVector [positions[1, i, t_idx], positions[2, i, t_idx], positions[3, i, t_idx]]
                        rj = @SVector [positions[1, j, t_idx], positions[2, j, t_idx], positions[3, j, t_idx]]
                        met = los_metrics(ri, rj; R_atm=R_atm)
                        los_ok = (!use_los) || ((!met.blocked) && (met.clearance >= atm_cl))
                        range_ok = (met.slant_range >= minR) && (met.slant_range <= maxR)
                        if los_ok && range_ok
                            1 # Line of sight is clear and within range
                        else
                            0 # Line of sight is blocked or out of range
                        end
                    end
                else
                    pair_los_transparency = GLMakie.@lift begin
                        t_idx = $current_time  
                        ri = @SVector [positions[1, i, t_idx], positions[2, i, t_idx], positions[3, i, t_idx]]
                        rj = @SVector [positions[1, j, t_idx], positions[2, j, t_idx], positions[3, j, t_idx]]
                        met = los_metrics(ri, rj; R_atm=R_atm)
                        range_ok = (met.slant_range >= minR) && (met.slant_range <= maxR)
                        if range_ok
                            1 # Within range
                        else
                            0 # Out of range
                        end
                    end
                end

                # Check if satellites i and j are single-pass pairs
                if Pm[i, j] != 0
                    GLMakie.lines!(ax, pair_points, color = :green, alpha = pair_los_transparency, linewidth = 2, label = "Single-pass pair")
                end

                # Check if satellites i and j are open-cavity pairs
                if haskey(cavity, (i, j))
                    GLMakie.lines!(ax, pair_points, color = :orange, alpha = pair_los_transparency, linewidth = 2, label = "Open-cavity pair")
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

    GLMakie.display(fig)

    return fig, (current_time = current_time, 
                is_playing = is_playing,
                animation_task = animation_task)
end