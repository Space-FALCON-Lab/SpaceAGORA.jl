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
        helper_trails: true for all helpers, false for none, or helper IDs such as [1, 3]
        target_trail: whether to show the target trail (default true)
        show_projections: project trails and position markers onto the XY/XZ/YZ boundary planes (default false)
        helper_projections: enable helper projection trails and markers (default true)
        target_projections: enable target projection trails and markers (default true)
        Δt: time between frames in seconds (default 0.1s)
        show_earth: whether to show the Earth sphere (default true)
        earth_radius: radius of the Earth sphere in meters (default R_EARTH)
        markersize: size of the satellite markers (default 15)
        figure_size: size of the figure in pixels (default (1200, 800))
        trail_alpha: transparency of the trail lines (default 0.6)
        animation_fps: desired frames per second for smooth animation (default 30.0)
        duration_seconds: playback duration (default 100.0)
        output_file: video path; nothing retains interactive playback
        radial_exaggeration: radial stretch above the reference altitude; 1 retains true scale
        reference_altitude_km: unscaled reference altitude; defaults to helper 1's initial altitude
        axis_limit_km: physical geocentric radius mapped to symmetric scene limits; nothing fits the scene
        radial_ticks_km: altitude references included in the scene extent; outer ticks are evenly spaced
        recorded_links: optional function of simulation time returning active, encounters, and kinds matrices

    Returns:
        figure and controls; saved videos have animation_task = nothing

"""
function animate_all_satellites_3d_smooth_helper_target(sol, p, helper_num;
    tail = 1000,
    helper_trails = true,
    target_trail::Bool = true,
    show_projections::Bool = false,
    helper_projections::Bool = true,
    target_projections::Bool = true,
    Δt = 0.1,
    show_earth = true,
    earth_radius = R_EARTH,
    markersize = 15,
    figure_size = (1200, 800),
    trail_alpha = 0.6,
    animation_fps = 30.0,
    duration_seconds = 100.0,
    output_file = nothing,
    radial_exaggeration = 1.0,
    reference_altitude_km = nothing,
    axis_limit_km = nothing,
    radial_ticks_km = nothing,
    scenario_caption = nothing,
    recorded_links = nothing)

    isfinite(animation_fps) && animation_fps > 0 || throw(ArgumentError("animation_fps must be finite and positive"))
    isfinite(duration_seconds) && duration_seconds > 0 || throw(ArgumentError("duration_seconds must be finite and positive"))
    isfinite(radial_exaggeration) && radial_exaggeration > 0 || throw(ArgumentError("radial_exaggeration must be finite and positive"))
    reference_altitude_km === nothing || (isfinite(reference_altitude_km) && reference_altitude_km >= 0) ||
        throw(ArgumentError("reference_altitude_km must be finite and nonnegative"))
    axis_limit_km === nothing || (isfinite(axis_limit_km) && axis_limit_km > earth_radius / 1e3) ||
        throw(ArgumentError("axis_limit_km must be finite and larger than Earth's radius in km"))
    radial_ticks_km === nothing || (radial_ticks_km isa AbstractVector{<:Real} && !isempty(radial_ticks_km) &&
        all(altitude -> isfinite(altitude) && altitude >= 0, radial_ticks_km)) ||
        throw(ArgumentError("radial_ticks_km must be a nonempty vector of finite, nonnegative altitudes"))
    if !(helper_trails isa Bool)
        helper_trails isa AbstractVector{<:Integer} || throw(ArgumentError("helper_trails must be true, false, or a vector of helper IDs"))
        all(helper -> !(helper isa Bool) && 1 <= helper <= helper_num, helper_trails) ||
            throw(ArgumentError("helper trail IDs must be between 1 and $helper_num"))
    end

    println("Setting up 3D animation for all satellites...")

    # Determine N satellites
    N = haskey(p, :N) ? p[:N] : (length(sol.u[1]) ÷ 6)
    #println("Animating $N satellites")

    # Create regular time grid for smooth animation
    t_start = sol.t[1]
    t_end = sol.t[end]
    t_regular = range(t_start, t_end; length=max(2, round(Int, animation_fps * duration_seconds)))
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

    physical_positions = positions
    reference_radius = reference_altitude_km === nothing ? norm(positions[:, 1, 1]) : earth_radius + reference_altitude_km * 1e3
    reference_radius >= earth_radius || throw(ArgumentError("reference altitude must be nonnegative"))
    display_radius(radius) = radial_exaggeration == 1 || radius <= reference_radius ? radius :
        reference_radius + radial_exaggeration * (radius - reference_radius)
    positions = copy(physical_positions)
    for frame in 1:n_points, satellite in 1:N
        radius = norm(physical_positions[:, satellite, frame])
        if radius > earth_radius
            positions[:, satellite, frame] .*= display_radius(radius) / radius
        end
    end
    tick_labels(values) = string.(round.(values ./ 1e3; digits=3))
    reference_km = round((reference_radius - earth_radius) / 1e3; digits=1)
    radial_altitudes = if radial_ticks_km === nothing
        highest_altitude = maximum(norm(physical_positions[:, satellite, frame]) for satellite in 1:N, frame in 1:n_points) / 1e3 - earth_radius / 1e3
        upper_step = max(10.0, 10.0 * ceil((highest_altitude - reference_km) / 50.0))
        sort!(unique(vcat([0.0, max(0.0, reference_km / 2)], reference_km .+ (0:5) .* upper_step)))
    else
        sort!(unique(Float64.(radial_ticks_km)))
    end
    radial_radii = display_radius.(earth_radius .+ radial_altitudes .* 1e3)
    scale_caption = radial_exaggeration == 1 ? "Linear axes; real km (Earth unscaled)" :
        "Radial scale: $(radial_exaggeration)x above $(reference_km) km; Earth and lower altitudes unscaled"

    # Set up GLMakie figure
    #GLMakie.activate!()
    fig = GLMakie.Figure(size = figure_size)
    GLMakie.Label(fig[0, 1], "Multi-Satellite Orbital Animation", fontsize=16, font=:bold,
        tellwidth=false)
    ax = GLMakie.Axis3(fig[1, 1],
               xlabel = radial_exaggeration == 1 ? "X [km]" : "X / geocentric reference [km]",
               ylabel = radial_exaggeration == 1 ? "Y [km]" : "Y / geocentric reference [km]",
               zlabel = radial_exaggeration == 1 ? "Z [km]" : "Z / geocentric reference [km]",
               xlabelsize = 14,
               ylabelsize = 14,
               zlabelsize = 14,
               xgridvisible = true,
               ygridvisible = true,
               zgridvisible = true,
               xtickformat = radial_exaggeration == 1 ? tick_labels : GLMakie.Makie.automatic,
               ytickformat = radial_exaggeration == 1 ? tick_labels : GLMakie.Makie.automatic,
               ztickformat = radial_exaggeration == 1 ? tick_labels : GLMakie.Makie.automatic,
               xticks = GLMakie.LinearTicks(7),
               yticks = GLMakie.LinearTicks(7),
               zticks = GLMakie.LinearTicks(7),
               xticklabelsize = 12,
               yticklabelsize = 12,
               zticklabelsize = 12,
               title = scenario_caption === nothing ? scale_caption : scenario_caption,
               titlefont = :regular,
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
        #println("Adding Earth sphere with texture...")
        # Download and load the Earth texture
        earth_img_path = joinpath(@__DIR__, "..", "input", "8k_earth_daymap.jpg")
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

    # Calculate the maximum absolute range across all dimensions
    radial_extent = radial_exaggeration == 1 ? earth_radius : maximum(radial_radii)
    max_limit = axis_limit_km === nothing ? 1.1 * max(radial_extent, maximum(abs, positions)) : display_radius(axis_limit_km * 1e3)
    maximum(abs, positions) <= max_limit || throw(ArgumentError("axis_limit_km clips the physical orbits; increase it"))
    radial_extent <= max_limit || throw(ArgumentError("axis_limit_km clips the radial ticks; increase it or reduce radial_ticks_km"))

    # Set symmetric limits for all axes
    xlims = (-max_limit, max_limit)
    ylims = (-max_limit, max_limit)
    zlims = (-max_limit, max_limit)

    # Apply the limits to the axis
    GLMakie.limits!(ax, xlims, ylims, zlims)

    if radial_exaggeration != 1
        tick_candidates = range(-0.9 * max_limit, 0.9 * max_limit; length=7)
        physical_radius(radius) = radius <= reference_radius ? radius :
            reference_radius + (radius - reference_radius) / radial_exaggeration
        largest_hundred = floor(Int, physical_radius(0.9 * max_limit) / 1e5)
        reference_kilometres = sort!(unique([100 * clamp(
            round(Int, sign(position) * physical_radius(abs(position)) / 1e5),
            -largest_hundred, largest_hundred) for position in tick_candidates]))
        tick_positions = [sign(value) * display_radius(abs(value) * 1e3) for value in reference_kilometres]
        reference_labels = string.(reference_kilometres)
        ax.xticks = (tick_positions, reference_labels)
        ax.yticks = (tick_positions, reference_labels)
        ax.zticks = (tick_positions, reference_labels)
    end

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
            visible = i <= helper_num ? (helper_trails isa Bool ? helper_trails : i in helper_trails) : target_trail,
            label = trail_label)

        # Plot current position
        GLMakie.scatter!(ax, current_pos,
            color = color,
            markersize = markersize,
            strokewidth = 1,
            strokecolor = :black,
            label = current_label)

        if show_projections && (i <= helper_num ? helper_projections : target_projections)
            for (plane, fixed_coordinate, plane_position) in (("XY", 3, -0.999 * max_limit),
                    ("XZ", 2, 0.999 * max_limit), ("YZ", 1, 0.999 * max_limit))
                project_point = point -> GLMakie.Point3f(ntuple(component ->
                    component == fixed_coordinate ? plane_position : point[component], 3))
                projected_trail = GLMakie.@lift project_point.($trail_points)
                projected_position = GLMakie.@lift project_point.($current_pos)
                GLMakie.lines!(ax, projected_trail; color=color, alpha=trail_alpha,
                    linewidth=1.5, linestyle=:dash,
                    visible=i <= helper_num ? (helper_trails isa Bool ? helper_trails : i in helper_trails) : target_trail,
                    label="$plane projection")
                GLMakie.scatter!(ax, projected_position; color=color, markersize=0.6 * markersize,
                    strokewidth=0.5, strokecolor=:black)
            end
        end
    end


    # Add RTN arrows for the target satellite
    target_sat = N  # Assuming the last satellite is the target

   # one line segment (er): exactly two points
    er_line = GLMakie.@lift begin
        t = $current_time
        r = @SVector [physical_positions[1,N,t], physical_positions[2,N,t], physical_positions[3,N,t]]
        v = @SVector [velocities[1,N,t], velocities[2,N,t], velocities[3,N,t]]
        er, _, _ = rtn_basis(r, v)
        r = @SVector [positions[1,N,t], positions[2,N,t], positions[3,N,t]]
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
        r = @SVector [physical_positions[1,N,t], physical_positions[2,N,t], physical_positions[3,N,t]]
        v = @SVector [velocities[1,N,t], velocities[2,N,t], velocities[3,N,t]]
        er, _, _ = rtn_basis(r, v)
        r = @SVector [positions[1,N,t], positions[2,N,t], positions[3,N,t]]
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
        r = @SVector [physical_positions[1,N,t], physical_positions[2,N,t], physical_positions[3,N,t]]
        v = @SVector [velocities[1,N,t], velocities[2,N,t], velocities[3,N,t]]
        _, et, _ = rtn_basis(r, v)
        r = @SVector [positions[1,N,t], positions[2,N,t], positions[3,N,t]]
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
        r = @SVector [physical_positions[1,N,t], physical_positions[2,N,t], physical_positions[3,N,t]]
        v = @SVector [velocities[1,N,t], velocities[2,N,t], velocities[3,N,t]]
        _, _, en = rtn_basis(r, v)
        r = @SVector [positions[1,N,t], positions[2,N,t], positions[3,N,t]]
        s = 1e6
        GLMakie.Point3f.(
            [r[1], r[1] + s*en[1]],
            [r[2], r[2] + s*en[2]],
            [r[3], r[3] + s*en[3]],
        )
    end
    GLMakie.lines!(ax, en_line, color = :blue, linewidth = 2)

    # Record satellite pair points dynamically with LOS status
    link_status = GLMakie.Observable("Laser links: OFF")
    laser_legend_names = ("Single-pass pair", "Open-cavity pair")
    laser_legend_labels = [GLMakie.Observable("$name: OFF") for name in laser_legend_names]
    laser_legend_labels[2][] = "Laser active"
    laser_legend_colors = [GLMakie.Observable(GLMakie.Makie.to_color(:gray)) for name in laser_legend_names]
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
            if recorded_links !== nothing
                links = recorded_links(t_regular[t_idx])
                A, encounters, H = links.active, links.encounters, links.kinds
            else
            # initialize fresh busy matrix for this frame
            forces, H = laser_forces(sol(t_regular[t_idx]), p)
            # per-pair alpha (Float32)
            A = fill(0f0, (N, N)) # A means transparency alpha
            encounters = fill(false, (N, N))

            for i in 1:N
                for j in i+1:N
                    # LOS / range
                    ri = @SVector [physical_positions[1, i, t_idx], physical_positions[2, i, t_idx], physical_positions[3, i, t_idx]]
                    rj = @SVector [physical_positions[1, j, t_idx], physical_positions[2, j, t_idx], physical_positions[3, j, t_idx]]
                    met = los_metrics(ri, rj; R_atm = R_atm)
                    los_ok = use_los ? ((!met.blocked) && (met.clearance >= atm_cl)) : true
                    range_ok = (met.slant_range >= minR) && (met.slant_range <= maxR)

                    configured = Pm[i, j] > 0 || Pm[j, i] > 0 ||
                        haskey(cavity, (i, j)) || haskey(cavity, (j, i))
                    active = (haskey(forces, (i, j)) && !iszero(forces[(i, j)])) ||
                        (haskey(forces, (j, i)) && !iszero(forces[(j, i)]))
                    A[i, j] = A[j, i] = active ? 1f0 : 0f0
                    encounters[i, j] = encounters[j, i] = (configured && los_ok && range_ok) || active
                end
            end
            end

            # publish the busy matrix for the frame
            current_helpers[] = H
            active_pairs = ["$helper-$target" for helper in 1:N for target in helper+1:N if A[helper, target] > 0]
            link_status[] = isempty(active_pairs) ? "Laser links: OFF" : "Laser links: ON (satellites $(join(active_pairs, ", ")))"
            for (kind, color) in enumerate((:green, :orange))
                pairs = ["$helper-$target" for helper in 1:N for target in helper+1:N
                    if A[helper, target] > 0 && H[helper, target] == kind]
                if kind == 1
                    laser_legend_labels[kind][] = isempty(pairs) ? "$(laser_legend_names[kind]): OFF" :
                        "$(laser_legend_names[kind]): satellites $(join(pairs, ", "))"
                end
                laser_legend_colors[kind][] = GLMakie.Makie.to_color(isempty(pairs) ? :gray : color)
            end
            (; active=A, encounters, kinds=H)
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
                pair_alpha = GLMakie.@lift $frame_links.active[i, j]
                pair_visible = GLMakie.@lift $pair_alpha > 0

                if Pm[i, j] > 0 || Pm[j, i] > 0 || haskey(cavity, (i, j)) || haskey(cavity, (j, i))
                    encounter_visible = GLMakie.@lift $frame_links.encounters[i, j]
                    GLMakie.lines!(ax, pair_points, color=:gray, linestyle=:dash,
                        visible=encounter_visible, linewidth=2, label="Encounter detected")
                end
                if Pm[i, j] > 0 || Pm[j, i] > 0
                    single_visible = GLMakie.@lift $pair_visible && $frame_links.kinds[i, j] == 1
                    GLMakie.lines!(ax, pair_points, color = :green, alpha = pair_alpha,
                        visible = single_visible, linewidth = 2, label = "Single-pass pair")
                end
                if haskey(cavity, (i, j)) || haskey(cavity, (j, i))
                    cavity_visible = GLMakie.@lift $pair_visible && $frame_links.kinds[i, j] == 2
                    GLMakie.lines!(ax, pair_points, color = :orange, alpha = pair_alpha,
                        visible = cavity_visible, linewidth = 2, label = "Open-cavity pair")
                end
            end
        end
    end
    # Add legend
    legend_plots, legend_labels = GLMakie.Makie.get_labeled_plots(ax; merge=false, unique=true)
    legend_entries = map(label -> !(label in ("XY projection", "XZ projection", "YZ projection")), legend_labels)
    legend_plots = legend_plots[legend_entries]
    legend_labels = legend_labels[legend_entries]
    legend_contents = Any[legend_plots...]
    reactive_labels = Any[legend_labels...]
    legend_active = Dict{String,Any}()
    for (kind, name) in enumerate(laser_legend_names)
        entry_index = findfirst(==(name), legend_labels)
        if entry_index !== nothing
            legend_contents[entry_index] = GLMakie.LineElement(color=laser_legend_colors[kind], linewidth=2)
            if kind == 2
                active = GLMakie.@lift $(laser_legend_colors[kind]) != GLMakie.Makie.to_color(:gray)
                legend_active["Laser active"] = active
                sample_color = GLMakie.@lift $active ? :orange : :white
                legend_contents[entry_index] = GLMakie.LineElement(color=sample_color, linewidth=2)
            end
            reactive_labels[entry_index] = laser_legend_labels[kind]
        end
    end
    encounter_entry = findfirst(==("Encounter detected"), legend_labels)
    if encounter_entry !== nothing
        reactive_labels[encounter_entry] = "Encounter detected"
        active = GLMakie.@lift any($frame_links.encounters)
        legend_active["Encounter detected"] = active
        sample_color = GLMakie.@lift $active ? :gray : :white
        legend_contents[encounter_entry] = GLMakie.LineElement(color=sample_color, linestyle=:dash, linewidth=2)
    end
    legend = GLMakie.axislegend(ax, legend_contents, reactive_labels; position=:lt)
    for entry in legend.entrygroups[][1][2]
        if haskey(legend_active, entry.label[])
            active = legend_active[entry.label[]]
            entry.attributes[:labelcolor] = GLMakie.@lift $active ? :black : :white
        end
    end
    GLMakie.notify(legend.entrygroups)
    entry_grid = only(item.content for item in legend.grid.content if item.content isa GLMakie.GridLayout)
    for (label, active) in legend_active
        text = only(item.content for item in entry_grid.content if item.content isa GLMakie.Label && item.content.text[] == label)
        bounds = GLMakie.lift(text.layoutobservables.computedbbox, legend.patchsize, legend.patchlabelgap) do text_bounds, patch_size, gap
            sample_width = patch_size[1] + gap
            GLMakie.Rect2f(text_bounds.origin[1] - sample_width - 2, text_bounds.origin[2] - 1,
                text_bounds.widths[1] + sample_width + 4, text_bounds.widths[2] + 2)
        end
        inactive = GLMakie.@lift !$active
        GLMakie.Box(legend.scene; bbox=bounds, color=:gray, visible=inactive, strokewidth=0, z=-6)
    end

    if output_file !== nothing
        time_text = GLMakie.@lift "t = $(round(t_regular[$current_time], digits=1)) s | $($link_status)"
        GLMakie.map!(ax.title, time_text) do status
            "$(scenario_caption === nothing ? scale_caption : scenario_caption)\n$status"
        end
        mkpath(dirname(abspath(output_file)))
        GLMakie.record(fig, output_file, eachindex(t_regular); framerate=animation_fps, visible=false) do frame
            current_time[] = frame
        end
        println("Saved animation: $output_file")
        return fig, (current_time=current_time, is_playing=GLMakie.Observable(false), animation_task=nothing)
    end

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
        target_distance = norm([physical_positions[1, N, t_idx], physical_positions[2, N, t_idx], physical_positions[3, N, t_idx]])
        "t = $(round(t_regular[t_idx], digits=1)) s, Target Distance = $(round(target_distance / 1e3, digits=1)) km"
    end
    GLMakie.map!(ax.title, time_text) do status
        "$(scenario_caption === nothing ? scale_caption : scenario_caption)\n$status"
    end


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