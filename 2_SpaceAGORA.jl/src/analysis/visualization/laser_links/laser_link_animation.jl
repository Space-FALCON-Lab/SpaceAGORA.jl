module LaserLinkAnimation

export animate_oracle_satellites_3d

using LinearAlgebra
using StaticArrays

# Self-contained RTN basis (not imported from SimulationModel — this module is
# loaded optionally from the example script, outside the package include tree).
@inline function _rtn_basis(r::SVector{3,Float64}, v::SVector{3,Float64})
    dcm = rtn_dcm_from_inertial(r, v)  # rtn_dcm_from_inertial is in scope via `using SpaceAGORA`
    return SVector{3,Float64}(dcm[:,1]), SVector{3,Float64}(dcm[:,2]), SVector{3,Float64}(dcm[:,3])
end

# Constants
const R_EARTH  = 6.371e6   # m
const R_ATMDEF = 6.471e6   # m

"""
    animate_oracle_satellites_3d(sol, helper_num; kwargs...)

Interactive 3D animation of the ORACLE open-cavity constellation.
Requires GLMakie to be loaded before calling (see example script).

# Arguments
- `sol`          : SpaceAGORA ODE solution (supports interpolation via `sol(t)`)
- `helper_num`   : number of helper satellites; the last satellite (index N) is the target
- `tail`         : number of previous time points to show in the trail (default 1000)
- `Δt`           : time between frames in seconds (default 0.1)
- `show_earth`   : whether to show the textured Earth sphere (default true)
- `earth_texture_path` : path to the earth texture jpg (default next to this file)
- `earth_radius` : radius of the Earth sphere in metres (default R_EARTH)
- `markersize`   : marker size for satellite dots (default 15)
- `figure_size`  : figure size in pixels (default (1200, 800))
- `trail_alpha`  : trail transparency (default 0.6)
- `animation_fps`: target frames per second (default 30.0)

Returns `(fig, controls)` where `controls` is a named tuple with
`current_time`, `is_playing`, `animation_task`, `screen`.
"""
function animate_oracle_satellites_3d(sol, helper_num;
    tail::Int         = 1000,
    Δt::Float64       = 0.1,
    show_earth::Bool  = true,
    earth_texture_path::String = joinpath(@__DIR__, "8k_earth_daymap.jpg"),
    earth_radius::Float64 = R_EARTH,
    markersize::Real  = 15,
    figure_size::Tuple{Int,Int} = (1200, 800),
    trail_alpha::Float64 = 0.6,
    animation_fps::Float64 = 30.0)

    # GLMakie must already be loaded by the caller; access via the module name.
    GLMakie = Base.require(Base.PkgId(Base.UUID("e9467ef8-e4e7-5192-8a1a-b1aee30e663a"), "GLMakie"))
    FileIO  = Base.require(Base.PkgId(Base.UUID("5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"), "FileIO"))
    GeometryBasics = Base.require(Base.PkgId(Base.UUID("5c1252a2-5f33-56bf-86c9-59e7332b4326"), "GeometryBasics"))

    N = haskey(sol, :sc) ? length(sol.u[1].sc) : (helper_num + 1)

    t_start = sol.t[1]; t_end = sol.t[end]
    dt_anim = (t_end - t_start) / (animation_fps * 100)
    t_regular = t_start:dt_anim:t_end
    n_points  = length(t_regular)

    positions  = zeros(3, N, n_points)
    velocities = zeros(3, N, n_points)
    for (k, t) in enumerate(t_regular)
        u_interp = sol(t)
        for i in 1:N
            sc = u_interp.sc[i]
            positions[1,i,k]  = sc.pos[1]; positions[2,i,k]  = sc.pos[2]; positions[3,i,k]  = sc.pos[3]
            velocities[1,i,k] = sc.vel[1]; velocities[2,i,k] = sc.vel[2]; velocities[3,i,k] = sc.vel[3]
        end
    end

    fig = GLMakie.Figure(size=figure_size)
    ax  = GLMakie.Axis3(fig[1,1],
        xlabel="X [m]", ylabel="Y [m]", zlabel="Z [m]",
        title="ORACLE Satellite Constellation (Interpolated)",
        aspect=(1,1,1))

    if show_earth
        sphere = GeometryBasics.uv_normal_mesh(
            GeometryBasics.Tesselation(
                GeometryBasics.Sphere(GeometryBasics.Point3f(0f0, 0f0, 0f0), Float32(earth_radius)),
                64))
        if isfile(earth_texture_path)
            img = FileIO.load(earth_texture_path)
            GLMakie.mesh!(ax, sphere; color=circshift(img, (0,3000)), ssao=true, alpha=0.5, transparency=true)
        else
            GLMakie.mesh!(ax, sphere; color=:steelblue, alpha=0.5, transparency=true)
        end
    end

    current_time = GLMakie.Observable(1)

    all_pos   = reshape(positions, 3, :)
    max_limit = maximum(abs.(all_pos)) * 1.1
    max_limit = max(max_limit, 1_500_000.0)
    GLMakie.limits!(ax, (-max_limit, max_limit), (-max_limit, max_limit), (-max_limit, max_limit))

    for i in 1:N
        color = i <= helper_num ? :blue : :red
        label_trail   = i <= helper_num ? "Helper Sat" : "Target Sat"
        label_current = label_trail

        trail_pts = GLMakie.@lift begin
            t_idx = $current_time
            s_idx = max(1, t_idx - tail)
            GLMakie.Point3f.(positions[1,i,s_idx:t_idx], positions[2,i,s_idx:t_idx], positions[3,i,s_idx:t_idx])
        end
        cur_pt = GLMakie.@lift [GLMakie.Point3f(positions[1,i,$current_time], positions[2,i,$current_time], positions[3,i,$current_time])]

        GLMakie.lines!(ax, trail_pts; color=color, alpha=trail_alpha, linewidth=2, label=label_trail)
        GLMakie.scatter!(ax, cur_pt; color=color, markersize=markersize, strokewidth=1, strokecolor=:black, label=label_current)
    end

    # RTN frame arrows for the target satellite (last in array = index N)
    for (axis_idx, axis_color) in ((1,:red), (2,:green), (3,:blue))
        line_pts = GLMakie.@lift begin
            t = $current_time
            r = @SVector [positions[1,N,t], positions[2,N,t], positions[3,N,t]]
            v = @SVector [velocities[1,N,t], velocities[2,N,t], velocities[3,N,t]]
        rhat, that, nhat = _rtn_basis(r, v)
            basis = (rhat, that, nhat)
            e = basis[axis_idx]; s = 1e6
            GLMakie.Point3f.([r[1], r[1]+s*e[1]], [r[2], r[2]+s*e[2]], [r[3], r[3]+s*e[3]])
        end
        GLMakie.lines!(ax, line_pts; color=axis_color, linewidth=2)
    end

    GLMakie.axislegend(ax; position=:lt, unique=true)

    ctrl_grid = fig[2,1] = GLMakie.GridLayout()
    is_playing  = GLMakie.Observable(false)
    play_btn    = GLMakie.Button(ctrl_grid[1,1], label=GLMakie.@lift($is_playing ? "|| Pause" : "=> Play"))
    GLMakie.Label(ctrl_grid[1,2], "Speed:")
    speed_slider = GLMakie.Slider(ctrl_grid[1,3], range=0.1:0.1:5.0, value=1.0)
    GLMakie.Label(ctrl_grid[1,5], "Time:")
    time_slider  = GLMakie.Slider(ctrl_grid[1,6], range=1:n_points, value=1)
    reset_btn    = GLMakie.Button(ctrl_grid[1,7], label="<<= Reset")

    time_text = GLMakie.@lift begin
        t_idx = $current_time
        d = norm([positions[1,N,t_idx], positions[2,N,t_idx], positions[3,N,t_idx]])
        "t = $(round(t_regular[t_idx], digits=1)) s, Target Distance = $(round(d/1e3, digits=1)) km"
    end
    GLMakie.Label(fig[0,1], time_text; tellwidth=false, fontsize=16)

    GLMakie.on(play_btn.clicks)  do _; is_playing[] = !is_playing[]; end
    GLMakie.on(time_slider.value) do val; current_time[] = Int(val); end
    GLMakie.on(reset_btn.clicks) do _
        current_time[] = 1
        GLMakie.set_close_to!(time_slider, 1)
        is_playing[] = false
    end

    target_dt = 1.0 / animation_fps
    animation_task = @async begin
        while true
            if is_playing[]
                if current_time[] < n_points
                    speed_val = GLMakie.to_value(speed_slider.value)
                    sleep(target_dt / max(speed_val, 0.1))
                    new_time = min(current_time[] + 1, n_points)
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
    return fig, (current_time=current_time, is_playing=is_playing,
                 animation_task=animation_task, screen=screen)
end

end # module LaserLinkAnimation
