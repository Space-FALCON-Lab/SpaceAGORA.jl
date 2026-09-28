include(joinpath(@__DIR__, "common.jl"))

using CSV
using ComponentArrays
using DataFrames
using LinearAlgebra
import PlotlyJS
using Plots
using Printf
using Random
using StaticArrays

"""
    Robot_Arm_Planner_Cloth_Demo.jl

Standalone smoke/demo case for the robotics arm port:

- builds the canonical `ClothArmModel`
- chooses a reachable inspection target and random spherical path obstacles
- solves IK and generates a HYPR joint-space arm plan around the obstacles
- realizes the plan through the compliant Cloth multibody dynamics backend
- saves joint, end-effector, and tracking-error CSV/PNG artifacts

Run with:

```text
julia --project=. examples/Robot_Arm_Planner_Cloth_Demo.jl
```
"""

const OUT_DIR = joinpath(REPO_ROOT, "output", "robot_arm_planner_cloth_demo")
mkpath(OUT_DIR)

const _CURRENT_OBSTACLES = Ref{Vector{RobotArmSphereObstacle}}(RobotArmSphereObstacle[])
const RPO_3U_CUBESAT_DIMS_M = SVector{3, Float64}(0.1, 0.1, 0.3)
const RPO_3U_CUBESAT_VISUAL_DIMS_M = SVector{3, Float64}(0.3, 0.1, 0.1)
const RPO_3U_CUBESAT_DRY_MASS_KG = 5.0
const RPO_3U_CUBESAT_PROP_MASS_KG = 0.2
const RPO_3U_CUBESAT_MAX_AXIS_THRUST_N = 0.05
const RPO_3U_CUBESAT_RW_MAX_TORQUE_NM = 0.25
const RPO_3U_CUBESAT_ISP_S = 60.0
const STANDARD_GRAVITY_MPS2 = 9.80665
const BASE_STABILIZATION_ENABLED = true

struct CoupledBaseSimulation
    t_s::Vector{Float64}
    states::Vector{ComponentVector{Float64}}
    end_effector_positions::Matrix{Float64}
    tracking_error_m::Vector{Float64}
    base_positions_m::Matrix{Float64}
    base_velocities_m_s::Matrix{Float64}
    base_omega_rad_s::Matrix{Float64}
    base_attitude_error_rad::Vector{Float64}
    reaction_force_n::Matrix{Float64}
    reaction_torque_nm::Matrix{Float64}
    stabilization_thrust_n::Matrix{Float64}
    stabilization_rw_torque_nm::Matrix{Float64}
    thrust_saturation_ratio::Vector{Float64}
    rw_saturation_ratio::Vector{Float64}
end

function _save_csv(path::String, df::DataFrame)
    CSV.write(path, df)
    println("Saved ", abspath(path))
    return path
end

function _save_plot(path::String, plt)
    savefig(plt, path)
    println("Saved ", abspath(path))
    return path
end

function _save_html_plot(path::String, plt)
    PlotlyJS.savefig(plt, path)
    println("Saved ", abspath(path))
    return path
end

function _env_int(name::AbstractString, default::Integer)
    token = strip(get(ENV, name, ""))
    isempty(token) && return Int(default)
    return parse(Int, token)
end

function _env_bool(name::AbstractString, default::Bool)
    token = lowercase(strip(get(ENV, name, "")))
    isempty(token) && return default
    token in ("1", "true", "yes", "on") && return true
    token in ("0", "false", "no", "off") && return false
    error("Environment variable $name must be boolean-like, got '$token'.")
end

function _joint_dataframe(plan::RobotArmPlan)
    df = DataFrame(time_s=plan.t_ref_s)
    for j in axes(plan.q_ref, 1)
        df[!, Symbol("q$j" * "_rad")] = plan.q_ref[j, :]
        df[!, Symbol("dq$j" * "_rad_s")] = plan.dq_ref[j, :]
        df[!, Symbol("ddq$j" * "_rad_s2")] = plan.ddq_ref[j, :]
    end
    return df
end

function _ee_dataframe(plan::RobotArmPlan, sim::ClothRobotArmSimulation)
    df = DataFrame(time_s=sim.trajectory.t_s)
    df.ee_ref_x_m = sim.reference_end_effector_positions[1, :]
    df.ee_ref_y_m = sim.reference_end_effector_positions[2, :]
    df.ee_ref_z_m = sim.reference_end_effector_positions[3, :]
    df.ee_cloth_x_m = sim.end_effector_positions[1, :]
    df.ee_cloth_y_m = sim.end_effector_positions[2, :]
    df.ee_cloth_z_m = sim.end_effector_positions[3, :]
    df.tracking_error_m = sim.tracking_error_m
    df.reference_obstacle_clearance_m = _obstacle_clearance_series(sim.reference_end_effector_positions, _CURRENT_OBSTACLES[])
    df.cloth_obstacle_clearance_m = _obstacle_clearance_series(sim.end_effector_positions, _CURRENT_OBSTACLES[])
    return df
end

function _base_dataframe(base_sim::CoupledBaseSimulation)
    df = DataFrame(time_s=base_sim.t_s)
    for (axis, label) in enumerate(("x", "y", "z"))
        df[!, Symbol("base_pos_" * label * "_m")] = base_sim.base_positions_m[axis, :]
        df[!, Symbol("base_vel_" * label * "_m_s")] = base_sim.base_velocities_m_s[axis, :]
        df[!, Symbol("base_omega_" * label * "_rad_s")] = base_sim.base_omega_rad_s[axis, :]
        df[!, Symbol("arm_reaction_force_" * label * "_n")] = base_sim.reaction_force_n[axis, :]
        df[!, Symbol("arm_reaction_torque_" * label * "_nm")] = base_sim.reaction_torque_nm[axis, :]
        df[!, Symbol("stabilization_thrust_" * label * "_n")] = base_sim.stabilization_thrust_n[axis, :]
        df[!, Symbol("stabilization_rw_torque_" * label * "_nm")] = base_sim.stabilization_rw_torque_nm[axis, :]
    end
    df.base_displacement_m = vec(sqrt.(sum(base_sim.base_positions_m .^ 2; dims=1)))
    df.base_speed_m_s = vec(sqrt.(sum(base_sim.base_velocities_m_s .^ 2; dims=1)))
    df.base_angular_rate_rad_s = vec(sqrt.(sum(base_sim.base_omega_rad_s .^ 2; dims=1)))
    df.base_attitude_error_rad = base_sim.base_attitude_error_rad
    df.coupled_ee_tracking_error_m = base_sim.tracking_error_m
    df.thrust_saturation_ratio = base_sim.thrust_saturation_ratio
    df.rw_saturation_ratio = base_sim.rw_saturation_ratio
    return df
end

function _obstacle_dataframe(obstacles::Vector{RobotArmSphereObstacle})
    return DataFrame(
        obstacle_id=collect(1:length(obstacles)),
        center_x_m=[obs.center[1] for obs in obstacles],
        center_y_m=[obs.center[2] for obs in obstacles],
        center_z_m=[obs.center[3] for obs in obstacles],
        radius_m=[obs.radius_m for obs in obstacles],
    )
end

function _obstacle_clearance_series(points::AbstractMatrix, obstacles::Vector{RobotArmSphereObstacle})
    isempty(obstacles) && return fill(Inf, size(points, 2))
    return [
        minimum(norm(SVector{3, Float64}(points[:, k]) - obs.center) - obs.radius_m for obs in obstacles)
        for k in axes(points, 2)
    ]
end

function _trapz(t, y)
    length(t) < 2 && return 0.0
    total = 0.0
    for i in 2:length(t)
        total += 0.5 * (Float64(y[i]) + Float64(y[i - 1])) * (Float64(t[i]) - Float64(t[i - 1]))
    end
    return total
end

function _stabilization_fuel_used_kg(
    base_sim::CoupledBaseSimulation;
    isp_s::Real=RPO_3U_CUBESAT_ISP_S,
    g0_mps2::Real=STANDARD_GRAVITY_MPS2,
)
    axis_thrust_sum = vec(sum(abs.(base_sim.stabilization_thrust_n); dims=1))
    return _trapz(base_sim.t_s, axis_thrust_sum) / max(Float64(isp_s) * Float64(g0_mps2), 1.0e-9)
end

function _fuel_used_pct(fuel_used_kg::Real, propellant_mass_kg::Real=RPO_3U_CUBESAT_PROP_MASS_KG)
    denom = Float64(propellant_mass_kg)
    return isfinite(denom) && denom > 0.0 ? 100.0 * Float64(fuel_used_kg) / denom : NaN
end

function _random_unit_vector(rng)
    v = SVector{3, Float64}(randn(rng), randn(rng), randn(rng))
    return v / max(norm(v), eps(Float64))
end

function _random_sphere_obstacles_on_nominal_path(
    plan::RobotArmPlan;
    rng=MersenneTwister(740),
    count::Int=5,
    radius_range_m=(0.018, 0.034),
    lateral_offset_m=0.008,
    x_offset_range_m=(-0.4, 0.1),
    endpoint_safe_distance_m=0.006,
    endpoint_exclusion_fraction=0.22,
)
    n = size(plan.ee_ref, 2)
    lo = max(2, floor(Int, endpoint_exclusion_fraction * n))
    hi = min(n - 1, ceil(Int, (1.0 - endpoint_exclusion_fraction) * n))
    lo < hi || throw(ArgumentError("Plan does not have enough samples to place interior obstacles."))

    start = SVector{3, Float64}(plan.ee_ref[:, 1])
    goal = SVector{3, Float64}(plan.ee_ref[:, end])
    obstacles = RobotArmSphereObstacle[]
    used_idxs = Set{Int}()
    attempts = 0
    while length(obstacles) < count && attempts < 250
        attempts += 1
        idx = rand(rng, lo:hi)
        idx in used_idxs && continue
        radius = rand(rng) * (radius_range_m[2] - radius_range_m[1]) + radius_range_m[1]
        x_offset = x_offset_range_m[1] + rand(rng) * (x_offset_range_m[2] - x_offset_range_m[1])
        center = SVector{3, Float64}(plan.ee_ref[:, idx]) +
            SVector{3, Float64}(x_offset, 0.0, 0.0) +
            lateral_offset_m * (2rand(rng) - 1.0) * _random_unit_vector(rng)
        endpoint_clearance = min(norm(center - start), norm(center - goal)) - radius
        endpoint_clearance > 0.04 || continue
        candidate = RobotArmSphereObstacle(center, radius)
        endpoint_q = hcat(plan.q_ref[:, 1], plan.q_ref[:, end])
        endpoint_stats = robot_arm_clearance_stats_from_samples(
            plan.model,
            plan.base_pose,
            endpoint_q,
            [candidate],
            endpoint_safe_distance_m,
        )
        endpoint_stats.violation_count == 0 || continue
        push!(used_idxs, idx)
        push!(obstacles, candidate)
    end
    length(obstacles) == count ||
        throw(ErrorException("Failed to place $count interior sphere obstacles without obstructing endpoints."))
    return obstacles
end

function _plot_joint_history(plan::RobotArmPlan)
    plt = plot(
        xlabel="Time (s)",
        ylabel="Joint angle (rad)",
        title="Robot Arm Joint Plan",
        lw=2,
        size=(950, 520),
    )
    for j in axes(plan.q_ref, 1)
        plot!(plt, plan.t_ref_s, plan.q_ref[j, :]; label="q$j")
    end
    return plt
end

function _plot_ee_path(plan::RobotArmPlan, sim::ClothRobotArmSimulation, obstacles::Vector{RobotArmSphereObstacle})
    lo, hi = _axis_limits_3d(plan, sim, obstacles)
    plt = plot(
        sim.reference_end_effector_positions[1, :],
        sim.reference_end_effector_positions[2, :],
        sim.reference_end_effector_positions[3, :];
        label="Reference EE",
        xlabel="x (m)",
        ylabel="y (m)",
        zlabel="z (m)",
        title="End-Effector Path",
        lw=3,
        camera=(35, 25),
        size=(850, 700),
        xlims=(lo[1], hi[1]),
        ylims=(lo[2], hi[2]),
        zlims=(lo[3], hi[3]),
        aspect_ratio=:equal,
    )
    plot!(
        plt,
        sim.end_effector_positions[1, :],
        sim.end_effector_positions[2, :],
        sim.end_effector_positions[3, :];
        label="Cloth realization",
        lw=2,
        ls=:dash,
    )
    scatter!(
        plt,
        [obs.center[1] for obs in obstacles],
        [obs.center[2] for obs in obstacles],
        [obs.center[3] for obs in obstacles];
        label="Sphere obstacle centers",
        ms=[max(5, 180 * obs.radius_m) for obs in obstacles],
        alpha=0.70,
    )
    scatter!(
        plt,
        [plan.target[1]],
        [plan.target[2]],
        [plan.target[3]];
        label="Target",
        ms=7,
    )
    return plt
end

function _plot_tracking_error(sim::ClothRobotArmSimulation)
    return plot(
        sim.trajectory.t_s,
        1.0e3 .* sim.tracking_error_m;
        label=false,
        xlabel="Time (s)",
        ylabel="EE tracking error (mm)",
        title="Cloth Realization Tracking Error",
        lw=2,
        size=(900, 420),
    )
end

function _plot_base_motion(base_sim::CoupledBaseSimulation)
    displacement_mm = 1.0e3 .* vec(sqrt.(sum(base_sim.base_positions_m .^ 2; dims=1)))
    attitude_mrad = 1.0e3 .* base_sim.base_attitude_error_rad
    plt = plot(
        base_sim.t_s,
        displacement_mm;
        label="hub displacement",
        xlabel="Time (s)",
        ylabel="Hub motion",
        title="3U CubeSat Base Motion from Arm Reaction",
        lw=2,
        size=(950, 460),
    )
    plot!(plt, base_sim.t_s, attitude_mrad; label="attitude error (mrad)", lw=2, ls=:dash)
    return plt
end

function _plot_stabilization_demand(base_sim::CoupledBaseSimulation)
    thrust_mn = 1.0e3 .* vec(sqrt.(sum(base_sim.stabilization_thrust_n .^ 2; dims=1)))
    torque_mnm = 1.0e3 .* vec(sqrt.(sum(base_sim.stabilization_rw_torque_nm .^ 2; dims=1)))
    plt = plot(
        base_sim.t_s,
        thrust_mn;
        label="ideal thrust demand",
        xlabel="Time (s)",
        ylabel="Command magnitude",
        title="Ideal Base Stabilization Demand (Not Applied)",
        lw=2,
        size=(950, 460),
    )
    plot!(plt, base_sim.t_s, torque_mnm; label="ideal RW torque demand", lw=2, ls=:dash)
    return plt
end

function _arm_polyline_from_plan(plan::RobotArmPlan, q)
    pose = cloth_fk(plan.model, plan.base_pose, q)
    pts = vcat([pose.joint_origins[1]], pose.link_tip_positions)
    return (
        x=[p[1] for p in pts],
        y=[p[2] for p in pts],
        z=[p[3] for p in pts],
    )
end

function _arm_polyline_from_cloth_state(plan::RobotArmPlan, x)
    n = length(plan.model.links)
    pts = Vector{SVector{3, Float64}}(undef, n + 1)
    pts[1] = plan.base_pose.position + _rot_local(plan.base_pose.quaternion) * plan.model.mount_offset_body
    for i in 1:n
        state = compliant_state_parts(x, i)
        link = plan.model.links[i]
        tip_body = link.vector_parent - link.com_offset_parent
        pts[i + 1] = state.r + _rot_local(state.q) * tip_body
    end
    return (
        x=[p[1] for p in pts],
        y=[p[2] for p in pts],
        z=[p[3] for p in pts],
    )
end

function _arm_polyline_from_tracked_motion(plan::RobotArmPlan, x, base_state)
    local_line = _arm_polyline_from_cloth_state(plan, x)
    R_base = _rot_local(base_state.q)
    r_base = SVector{3, Float64}(base_state.pos)
    pts = [
        r_base + R_base * SVector{3, Float64}(local_line.x[i], local_line.y[i], local_line.z[i])
        for i in eachindex(local_line.x)
    ]
    return (
        x=[p[1] for p in pts],
        y=[p[2] for p in pts],
        z=[p[3] for p in pts],
    )
end

function _base_prism_corners(position, quaternion)
    half = 0.5 .* RPO_3U_CUBESAT_VISUAL_DIMS_M
    local_corners = SVector{3, Float64}[
        SVector{3, Float64}(sx * half[1], sy * half[2], sz * half[3])
        for sx in (-1.0, 1.0) for sy in (-1.0, 1.0) for sz in (-1.0, 1.0)
    ]
    r = SVector{3, Float64}(position)
    R = _rot_local(quaternion)
    return [r + R * c for c in local_corners]
end

function _base_prism_extent_matrix(base_sim::CoupledBaseSimulation)
    isempty(base_sim.states) && return zeros(3, 0)
    pts = SVector{3, Float64}[]
    for idx in (1, length(base_sim.states))
        state = base_sim.states[idx]
        append!(pts, _base_prism_corners(state.pos, state.q))
    end
    return hcat(pts...)
end

function _base_prism_trace(state, name::AbstractString, color::AbstractString; opacity::Real=0.30)
    corners = _base_prism_corners(state.pos, state.q)
    # Corner order from _base_prism_corners:
    # 1=(-,-,-), 2=(-,-,+), 3=(-,+,-), 4=(-,+,+), 5=(+,-,-), 6=(+,-,+), 7=(+,+,-), 8=(+,+,+)
    faces = (
        (1, 2, 4), (1, 4, 3),
        (5, 7, 8), (5, 8, 6),
        (1, 5, 6), (1, 6, 2),
        (3, 4, 8), (3, 8, 7),
        (1, 3, 7), (1, 7, 5),
        (2, 6, 8), (2, 8, 4),
    )
    return PlotlyJS.mesh3d(
        x=[p[1] for p in corners],
        y=[p[2] for p in corners],
        z=[p[3] for p in corners],
        i=[f[1] - 1 for f in faces],
        j=[f[2] - 1 for f in faces],
        k=[f[3] - 1 for f in faces],
        name=name,
        color=color,
        opacity=opacity,
        flatshading=true,
        showscale=false,
    )
end

function _rot_local(q_raw)
    q = SVector{4, Float64}(q_raw) ./ max(norm(q_raw), eps(Float64))
    x, y, z, w = q
    return SMatrix{3, 3, Float64}(
        1 - 2y^2 - 2z^2, 2x * y + 2z * w, 2x * z - 2y * w,
        2x * y - 2z * w, 1 - 2x^2 - 2z^2, 2y * z + 2x * w,
        2x * z + 2y * w, 2y * z - 2x * w, 1 - 2x^2 - 2y^2,
    )
end

function _quat_mul_local(a, b)
    aq = SVector{4, Float64}(a)
    bq = SVector{4, Float64}(b)
    ax, ay, az, aw = aq
    bx, by, bz, bw = bq
    return SVector{4, Float64}(
        aw * bx + bw * ax + ay * bz - az * by,
        aw * by + bw * ay + az * bx - ax * bz,
        aw * bz + bw * az + ax * by - ay * bx,
        aw * bw - ax * bx - ay * by - az * bz,
    )
end

function _unit_quat_local(q_raw)
    q = SVector{4, Float64}(q_raw)
    nq = norm(q)
    return isfinite(nq) && nq > eps(Float64) ? q / nq : SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)
end

function _quat_derivative_local(q_raw, omega_body)
    q = _unit_quat_local(q_raw)
    ω = SVector{3, Float64}(omega_body)
    return 0.5 .* _quat_mul_local(q, SVector{4, Float64}(ω[1], ω[2], ω[3], 0.0))
end

function _quat_angle_error_local(q_raw)
    q = _unit_quat_local(q_raw)
    return 2.0 * atan(norm(SVector{3, Float64}(q[1], q[2], q[3])), abs(q[4]))
end

function _normalize_coupled_state_quaternions!(state)
    state.q .= _unit_quat_local(state.q)
    for i in axes(state.arm_q, 2)
        state.arm_q[:, i] .= _unit_quat_local(state.arm_q[:, i])
    end
    return state
end

function _coupled_end_effector_from_state(plan::RobotArmPlan, state)
    n = length(plan.model.links)
    link = plan.model.links[end]
    tip_body = link.vector_parent - link.com_offset_parent
    q = _unit_quat_local(state.arm_q[:, n])
    return SVector{3, Float64}(state.arm_r[:, n]) + _rot_local(q) * tip_body
end

function _rpo_3u_cubesat_inertia()
    root = Link(
        root=true,
        m=RPO_3U_CUBESAT_DRY_MASS_KG,
        dims=MVector{3, Float64}(RPO_3U_CUBESAT_DIMS_M),
        ref_area=0.03,
    )
    return root.inertia
end

function _coupled_base_initial_state(plan::RobotArmPlan)
    shape = merge((
        pos=zeros(3),
        vel=zeros(3),
        mass=RPO_3U_CUBESAT_DRY_MASS_KG + RPO_3U_CUBESAT_PROP_MASS_KG,
        heat_loads=zeros(1),
        q=Float64[0.0, 0.0, 0.0, 1.0],
        ω=zeros(3),
    ), coupled_cloth_robot_arm_state_shape(plan))
    state = ComponentVector(shape)
    initialize_coupled_cloth_robot_arm_state!(state, plan; t_s=0.0)
    return state
end

function _base_only_initial_state()
    return ComponentVector(
        pos=zeros(3),
        vel=zeros(3),
        q=Float64[0.0, 0.0, 0.0, 1.0],
        ω=zeros(3),
        mass=RPO_3U_CUBESAT_DRY_MASS_KG + RPO_3U_CUBESAT_PROP_MASS_KG,
    )
end

function _prescribed_cloth_state_for_reaction(plan::RobotArmPlan, x::AbstractVector)
    state = _coupled_base_initial_state(plan)
    n = length(plan.model.links)
    for i in 1:n
        part = compliant_state_parts(x, i)
        state.arm_r[:, i] .= part.r
        state.arm_q[:, i] .= part.q
        state.arm_v[:, i] .= part.v
        state.arm_ω[:, i] .= part.ω
    end
    return state
end

function _root_reaction_from_cloth_state(
    plan::RobotArmPlan,
    x::AbstractVector,
    t::Real;
    k_translation_n_m::Real,
    c_translation_n_s_m::Real,
    k_rotation_n_m_rad::Real,
    c_rotation_n_m_s_rad::Real,
)
    state = _prescribed_cloth_state_for_reaction(plan, x)
    du = zero(state)
    reaction_force = MVector{3, Float64}(0.0, 0.0, 0.0)
    reaction_torque = MVector{3, Float64}(0.0, 0.0, 0.0)
    assign_coupled_cloth_robot_arm_rhs!(
        du,
        state,
        plan,
        t,
        reaction_force,
        reaction_torque;
        k_translation_n_m=k_translation_n_m,
        c_translation_n_s_m=c_translation_n_s_m,
        k_rotation_n_m_rad=k_rotation_n_m_rad,
        c_rotation_n_m_s_rad=c_rotation_n_m_s_rad,
    )
    return SVector{3, Float64}(reaction_force), SVector{3, Float64}(reaction_torque)
end

function _advance_base_from_wrench(state, force, torque, dt::Real, inertia_tensor)
    next_state = deepcopy(state)
    h = Float64(dt)
    ω = SVector{3, Float64}(state.ω)
    τ = SVector{3, Float64}(torque)
    a = SVector{3, Float64}(force) / state.mass
    α = inertia_tensor \ (τ - cross(ω, inertia_tensor * ω))
    next_state.pos .= state.pos .+ h .* state.vel .+ 0.5h^2 .* a
    next_state.vel .= state.vel .+ h .* a
    next_state.q .= state.q .+ h .* _quat_derivative_local(state.q, ω)
    next_state.ω .= state.ω .+ h .* α
    next_state.q .= _unit_quat_local(next_state.q)
    return next_state
end

function _coupled_base_rhs!(
    du,
    state,
    plan::RobotArmPlan,
    t::Real,
    inertia_tensor;
    k_translation_n_m::Real,
    c_translation_n_s_m::Real,
    k_rotation_n_m_rad::Real,
    c_rotation_n_m_s_rad::Real,
)
    du .= 0.0
    reaction_force = MVector{3, Float64}(0.0, 0.0, 0.0)
    reaction_torque = MVector{3, Float64}(0.0, 0.0, 0.0)
    assign_coupled_cloth_robot_arm_rhs!(
        du,
        state,
        plan,
        t,
        reaction_force,
        reaction_torque;
        k_translation_n_m=k_translation_n_m,
        c_translation_n_s_m=c_translation_n_s_m,
        k_rotation_n_m_rad=k_rotation_n_m_rad,
        c_rotation_n_m_s_rad=c_rotation_n_m_s_rad,
    )
    ω = SVector{3, Float64}(state.ω)
    τ = SVector{3, Float64}(reaction_torque)
    du.pos .= state.vel
    du.vel .= SVector{3, Float64}(reaction_force) / state.mass
    du.mass = 0.0
    du.heat_loads .= 0.0
    du.q .= _quat_derivative_local(state.q, ω)
    du.ω .= inertia_tensor \ (τ - cross(ω, inertia_tensor * ω))
    return SVector{3, Float64}(reaction_force), τ
end

function _coupled_base_derivative(
    state,
    plan::RobotArmPlan,
    t::Real,
    inertia_tensor;
    k_translation_n_m::Real,
    c_translation_n_s_m::Real,
    k_rotation_n_m_rad::Real,
    c_rotation_n_m_s_rad::Real,
)
    du = zero(state)
    f, τ = _coupled_base_rhs!(
        du,
        state,
        plan,
        t,
        inertia_tensor;
        k_translation_n_m=k_translation_n_m,
        c_translation_n_s_m=c_translation_n_s_m,
        k_rotation_n_m_rad=k_rotation_n_m_rad,
        c_rotation_n_m_s_rad=c_rotation_n_m_s_rad,
    )
    return du, f, τ
end

function _advance_coupled_base_rk4(
    state,
    plan::RobotArmPlan,
    t::Real,
    dt::Real,
    inertia_tensor;
    k_translation_n_m::Real,
    c_translation_n_s_m::Real,
    k_rotation_n_m_rad::Real,
    c_rotation_n_m_s_rad::Real,
)
    h = Float64(dt)
    k1, _, _ = _coupled_base_derivative(
        state,
        plan,
        t,
        inertia_tensor;
        k_translation_n_m=k_translation_n_m,
        c_translation_n_s_m=c_translation_n_s_m,
        k_rotation_n_m_rad=k_rotation_n_m_rad,
        c_rotation_n_m_s_rad=c_rotation_n_m_s_rad,
    )
    s2 = _normalize_coupled_state_quaternions!(deepcopy(state) .+ 0.5h .* k1)
    k2, _, _ = _coupled_base_derivative(
        s2,
        plan,
        t + 0.5h,
        inertia_tensor;
        k_translation_n_m=k_translation_n_m,
        c_translation_n_s_m=c_translation_n_s_m,
        k_rotation_n_m_rad=k_rotation_n_m_rad,
        c_rotation_n_m_s_rad=c_rotation_n_m_s_rad,
    )
    s3 = _normalize_coupled_state_quaternions!(deepcopy(state) .+ 0.5h .* k2)
    k3, _, _ = _coupled_base_derivative(
        s3,
        plan,
        t + 0.5h,
        inertia_tensor;
        k_translation_n_m=k_translation_n_m,
        c_translation_n_s_m=c_translation_n_s_m,
        k_rotation_n_m_rad=k_rotation_n_m_rad,
        c_rotation_n_m_s_rad=c_rotation_n_m_s_rad,
    )
    s4 = _normalize_coupled_state_quaternions!(deepcopy(state) .+ h .* k3)
    k4, _, _ = _coupled_base_derivative(
        s4,
        plan,
        t + h,
        inertia_tensor;
        k_translation_n_m=k_translation_n_m,
        c_translation_n_s_m=c_translation_n_s_m,
        k_rotation_n_m_rad=k_rotation_n_m_rad,
        c_rotation_n_m_s_rad=c_rotation_n_m_s_rad,
    )
    next_state = deepcopy(state)
    next_state .+= (h / 6.0) .* (k1 .+ 2.0 .* k2 .+ 2.0 .* k3 .+ k4)
    return _normalize_coupled_state_quaternions!(next_state)
end

function simulate_coupled_base_motion(
    plan::RobotArmPlan;
    arm_sim::ClothRobotArmSimulation,
    dt_s::Real=0.001,
    duration_s::Real=plan.t_ref_s[end],
    k_translation_n_m::Real=8.0e3,
    c_translation_n_s_m::Real=40.0,
    k_rotation_n_m_rad::Real=25.0,
    c_rotation_n_m_s_rad::Real=1.2,
)
    times = arm_sim.trajectory.t_s
    inertia_tensor = _rpo_3u_cubesat_inertia()
    states = Vector{ComponentVector{Float64}}(undef, length(times))
    states[1] = _base_only_initial_state()
    ee = zeros(3, length(times))
    err = zeros(length(times))
    base_pos = zeros(3, length(times))
    base_vel = zeros(3, length(times))
    base_ω = zeros(3, length(times))
    base_att = zeros(length(times))
    reaction_force = zeros(3, length(times))
    reaction_torque = zeros(3, length(times))
    stab_force = zeros(3, length(times))
    stab_torque = zeros(3, length(times))
    thrust_sat = zeros(length(times))
    rw_sat = zeros(length(times))

    for k in eachindex(times)
        state = states[k]
        f, τ = _root_reaction_from_cloth_state(
            plan,
            arm_sim.trajectory.states[k],
            times[k],
            k_translation_n_m=k_translation_n_m,
            c_translation_n_s_m=c_translation_n_s_m,
            k_rotation_n_m_rad=k_rotation_n_m_rad,
            c_rotation_n_m_s_rad=c_rotation_n_m_s_rad,
        )
        local_ee = SVector{3, Float64}(arm_sim.end_effector_positions[:, k])
        ee_k = SVector{3, Float64}(state.pos) + _rot_local(state.q) * local_ee
        ee[:, k] .= ee_k
        err[k] = norm(ee_k - robot_arm_plan_sample(plan, times[k]).ee)
        base_pos[:, k] .= state.pos
        base_vel[:, k] .= state.vel
        base_ω[:, k] .= state.ω
        base_att[k] = _quat_angle_error_local(state.q)
        reaction_force[:, k] .= f
        reaction_torque[:, k] .= τ
        stab_force[:, k] .= -f
        stab_torque[:, k] .= -τ
        thrust_sat[k] = maximum(abs.(stab_force[:, k])) / RPO_3U_CUBESAT_MAX_AXIS_THRUST_N
        rw_sat[k] = maximum(abs.(stab_torque[:, k])) / RPO_3U_CUBESAT_RW_MAX_TORQUE_NM
        if k < length(times)
            step_dt = times[k + 1] - times[k]
            applied_force = BASE_STABILIZATION_ENABLED ? f + SVector{3, Float64}(stab_force[:, k]) : f
            applied_torque = BASE_STABILIZATION_ENABLED ? τ + SVector{3, Float64}(stab_torque[:, k]) : τ
            states[k + 1] = _advance_base_from_wrench(
                state,
                applied_force,
                applied_torque,
                step_dt,
                inertia_tensor,
            )
        end
    end

    return CoupledBaseSimulation(
        times,
        states,
        ee,
        err,
        base_pos,
        base_vel,
        base_ω,
        base_att,
        reaction_force,
        reaction_torque,
        stab_force,
        stab_torque,
        thrust_sat,
        rw_sat,
    )
end

function _obstacle_extent_matrix(obstacles::Vector{RobotArmSphereObstacle})
    isempty(obstacles) && return zeros(3, 0)
    pts = SVector{3, Float64}[]
    for obs in obstacles
        r = obs.radius_m
        push!(pts, obs.center + SVector{3, Float64}(r, 0.0, 0.0))
        push!(pts, obs.center - SVector{3, Float64}(r, 0.0, 0.0))
        push!(pts, obs.center + SVector{3, Float64}(0.0, r, 0.0))
        push!(pts, obs.center - SVector{3, Float64}(0.0, r, 0.0))
        push!(pts, obs.center + SVector{3, Float64}(0.0, 0.0, r))
        push!(pts, obs.center - SVector{3, Float64}(0.0, 0.0, r))
    end
    return hcat(pts...)
end

function _axis_limits_3d(
    plan::RobotArmPlan,
    sim::ClothRobotArmSimulation,
    obstacles::Vector{RobotArmSphereObstacle};
    base_sim::Union{Nothing, CoupledBaseSimulation}=nothing,
)
    all_xyz = hcat(
        plan.ee_ref,
        sim.end_effector_positions,
        base_sim === nothing ? zeros(3, 0) : base_sim.end_effector_positions,
        _obstacle_extent_matrix(obstacles),
        reshape(collect(plan.base_pose.position), 3, 1),
        base_sim === nothing ? zeros(3, 0) : base_sim.base_positions_m,
        base_sim === nothing ? zeros(3, 0) : _base_prism_extent_matrix(base_sim),
    )
    lo = vec(minimum(all_xyz; dims=2))
    hi = vec(maximum(all_xyz; dims=2))
    pad = max(0.08, 0.15 * maximum(hi .- lo))
    center = 0.5 .* (lo .+ hi)
    half_span = 0.5 * maximum(hi .- lo) + pad
    return (center .- half_span, center .+ half_span)
end

function _sphere_surface_trace(obs::RobotArmSphereObstacle, idx::Int)
    θ = collect(range(0.0, 2π, length=34))
    ϕ = collect(range(0.0, π, length=18))
    x = [obs.center[1] + obs.radius_m * cos(t) * sin(p) for p in ϕ, t in θ]
    y = [obs.center[2] + obs.radius_m * sin(t) * sin(p) for p in ϕ, t in θ]
    z = [obs.center[3] + obs.radius_m * cos(p) for p in ϕ, t in θ]
    return PlotlyJS.surface(
        x=x,
        y=y,
        z=z,
        name="sphere obstacle $idx",
        opacity=0.33,
        showscale=false,
        colorscale=[[0.0, "rgb(180,75,65)"], [1.0, "rgb(180,75,65)"]],
    )
end

function _plot_arm_motion_html(
    plan::RobotArmPlan,
    sim::ClothRobotArmSimulation,
    obstacles::Vector{RobotArmSphereObstacle},
    base_sim::CoupledBaseSimulation,
)
    traces = PlotlyJS.GenericTrace[]
    snapshot_count = min(14, length(plan.t_ref_s))
    idxs = unique(round.(Int, range(1, length(plan.t_ref_s), length=snapshot_count)))
    colors = range(0.15, 0.95, length=length(idxs))

    for (n, idx) in pairs(idxs)
        q = plan.q_ref[:, idx]
        arm_line = _arm_polyline_from_plan(plan, q)
        opacity = 0.25 + 0.55 * colors[n]
        push!(traces, PlotlyJS.scatter3d(
            x=arm_line.x,
            y=arm_line.y,
            z=arm_line.z,
            mode="lines+markers",
            name=@sprintf("reference arm t=%.2fs", plan.t_ref_s[idx]),
            line=PlotlyJS.attr(width=5, color=@sprintf("rgba(35,95,170,%.3f)", opacity)),
            marker=PlotlyJS.attr(size=4, color=@sprintf("rgba(35,95,170,%.3f)", opacity)),
            showlegend=(n == 1),
        ))
    end

    tracked_snapshot_idxs = unique(round.(Int, range(1, length(sim.trajectory.states), length=min(8, length(sim.trajectory.states)))))
    for (n, idx) in pairs(tracked_snapshot_idxs)
        cloth_line = _arm_polyline_from_tracked_motion(plan, sim.trajectory.states[idx], base_sim.states[idx])
        push!(traces, PlotlyJS.scatter3d(
            x=cloth_line.x,
            y=cloth_line.y,
            z=cloth_line.z,
            mode="lines+markers",
            name=idx == 1 ? "tracked arm start" : (idx == length(sim.trajectory.states) ? "tracked arm final" : @sprintf("tracked arm t=%.2fs", sim.trajectory.t_s[idx])),
            line=PlotlyJS.attr(
                width=idx in (1, length(sim.trajectory.states)) ? 6 : 4,
                color=idx == 1 ? "rgb(60,140,90)" : (idx == length(sim.trajectory.states) ? "rgb(210,85,55)" : @sprintf("rgba(90,90,90,%.3f)", 0.25 + 0.45 * n / length(tracked_snapshot_idxs))),
            ),
            marker=PlotlyJS.attr(size=idx in (1, length(sim.trajectory.states)) ? 5 : 3),
            showlegend=idx in (1, length(sim.trajectory.states)),
        ))
    end

    push!(traces, PlotlyJS.scatter3d(
        x=plan.ee_ref[1, :],
        y=plan.ee_ref[2, :],
        z=plan.ee_ref[3, :],
        mode="lines",
        name="reference EE path",
        line=PlotlyJS.attr(width=7, color="rgb(20,70,150)"),
    ))
    push!(traces, PlotlyJS.scatter3d(
        x=base_sim.end_effector_positions[1, :],
        y=base_sim.end_effector_positions[2, :],
        z=base_sim.end_effector_positions[3, :],
        mode="lines",
        name="tracked EE path with free base",
        line=PlotlyJS.attr(width=4, dash="dash", color="rgb(220,90,50)"),
    ))
    push!(traces, PlotlyJS.scatter3d(
        x=base_sim.base_positions_m[1, :],
        y=base_sim.base_positions_m[2, :],
        z=base_sim.base_positions_m[3, :],
        mode="lines",
        name="free hub path",
        line=PlotlyJS.attr(width=5, color="rgb(35,35,35)"),
    ))
    push!(traces, _base_prism_trace(base_sim.states[1], "3U base start", "rgb(80,120,170)"; opacity=0.26))
    push!(traces, _base_prism_trace(base_sim.states[end], "3U base final (free)", "rgb(230,130,40)"; opacity=0.38))
    for (idx, obs) in pairs(obstacles)
        push!(traces, _sphere_surface_trace(obs, idx))
        push!(traces, PlotlyJS.scatter3d(
            x=[obs.center[1]],
            y=[obs.center[2]],
            z=[obs.center[3]],
            mode="markers",
            name="obstacle center $idx",
            marker=PlotlyJS.attr(size=3, color="rgb(120,40,35)"),
            showlegend=false,
        ))
    end
    push!(traces, PlotlyJS.scatter3d(
        x=[plan.target[1]],
        y=[plan.target[2]],
        z=[plan.target[3]],
        mode="markers",
        name="inspection target",
        marker=PlotlyJS.attr(size=8, color="rgb(180,40,80)", symbol="diamond"),
    ))

    lo, hi = _axis_limits_3d(plan, sim, obstacles; base_sim=base_sim)
    layout = PlotlyJS.Layout(
        title="Robot Arm Motion with Free 3U Base Drift",
        scene=PlotlyJS.attr(
            aspectmode="cube",
            xaxis=PlotlyJS.attr(title="x (m)", range=[lo[1], hi[1]]),
            yaxis=PlotlyJS.attr(title="y (m)", range=[lo[2], hi[2]]),
            zaxis=PlotlyJS.attr(title="z (m)", range=[lo[3], hi[3]]),
        ),
        legend=PlotlyJS.attr(x=0.02, y=0.98),
        margin=PlotlyJS.attr(l=0, r=0, t=55, b=0),
    )
    return PlotlyJS.Plot(traces, layout)
end

function _planner_metrics_html(
    plan::RobotArmPlan,
    sim::ClothRobotArmSimulation,
    obstacles::Vector{RobotArmSphereObstacle},
    hypr::RobotArmHYPRResult,
    base_sim::CoupledBaseSimulation,
)
    traces = PlotlyJS.GenericTrace[]
    for j in axes(plan.q_ref, 1)
        push!(traces, PlotlyJS.scatter(x=plan.t_ref_s, y=plan.q_ref[j, :], mode="lines", name="q$j rad", xaxis="x1", yaxis="y1"))
        push!(traces, PlotlyJS.scatter(x=plan.t_ref_s, y=plan.dq_ref[j, :], mode="lines", name="dq$j rad/s", xaxis="x2", yaxis="y2"))
        push!(traces, PlotlyJS.scatter(x=plan.t_ref_s, y=plan.ddq_ref[j, :], mode="lines", name="ddq$j rad/s^2", xaxis="x3", yaxis="y3"))
    end
    ee_speed = [0.0; vec(sqrt.(sum(diff(plan.ee_ref; dims=2).^2; dims=1))) ./ diff(plan.t_ref_s)]
    ref_clearance_mm = 1.0e3 .* _obstacle_clearance_series(plan.ee_ref, obstacles)
    cloth_clearance_mm = 1.0e3 .* _obstacle_clearance_series(sim.end_effector_positions, obstacles)
    push!(traces, PlotlyJS.scatter(x=plan.t_ref_s, y=ee_speed, mode="lines", name="EE speed m/s", xaxis="x4", yaxis="y4"))
    push!(traces, PlotlyJS.scatter(x=sim.trajectory.t_s, y=1.0e3 .* sim.tracking_error_m, mode="lines", name="Cloth EE error mm", xaxis="x4", yaxis="y4"))
    push!(traces, PlotlyJS.scatter(x=plan.t_ref_s, y=ref_clearance_mm, mode="lines", name="Reference obstacle clearance mm", xaxis="x4", yaxis="y4"))
    push!(traces, PlotlyJS.scatter(x=sim.trajectory.t_s, y=cloth_clearance_mm, mode="lines", name="Cloth obstacle clearance mm", xaxis="x4", yaxis="y4"))
    push!(traces, PlotlyJS.scatter(
        x=collect(1:length(hypr.cost_history)),
        y=hypr.cost_history,
        mode="lines+markers",
        name="HYPR best cost",
        xaxis="x5",
        yaxis="y5",
    ))
    base_disp_mm = 1.0e3 .* vec(sqrt.(sum(base_sim.base_positions_m .^ 2; dims=1)))
    base_att_mrad = 1.0e3 .* base_sim.base_attitude_error_rad
    thrust_mn = 1.0e3 .* vec(sqrt.(sum(base_sim.stabilization_thrust_n .^ 2; dims=1)))
    rw_mnm = 1.0e3 .* vec(sqrt.(sum(base_sim.stabilization_rw_torque_nm .^ 2; dims=1)))
    push!(traces, PlotlyJS.scatter(x=base_sim.t_s, y=base_disp_mm, mode="lines", name="Base displacement mm", xaxis="x6", yaxis="y6"))
    push!(traces, PlotlyJS.scatter(x=base_sim.t_s, y=base_att_mrad, mode="lines", name="Base attitude error mrad", xaxis="x6", yaxis="y6"))
    push!(traces, PlotlyJS.scatter(x=base_sim.t_s, y=thrust_mn, mode="lines", name="Ideal cancellation thrust mN (not applied)", xaxis="x6", yaxis="y6"))
    push!(traces, PlotlyJS.scatter(x=base_sim.t_s, y=rw_mnm, mode="lines", name="Ideal cancellation RW torque mN-m (not applied)", xaxis="x6", yaxis="y6"))

    layout = PlotlyJS.Layout(
        title="Robot Arm HYPR Planner Metrics with 3U Base Reaction",
        grid=PlotlyJS.attr(rows=3, columns=2, pattern="independent"),
        xaxis=PlotlyJS.attr(title="time (s)"),
        yaxis=PlotlyJS.attr(title="joint angle (rad)"),
        xaxis2=PlotlyJS.attr(title="time (s)"),
        yaxis2=PlotlyJS.attr(title="joint velocity (rad/s)"),
        xaxis3=PlotlyJS.attr(title="time (s)"),
        yaxis3=PlotlyJS.attr(title="joint acceleration (rad/s^2)"),
        xaxis4=PlotlyJS.attr(title="time (s)"),
        yaxis4=PlotlyJS.attr(title="EE speed / error / clearance"),
        xaxis5=PlotlyJS.attr(title="HYPR iteration"),
        yaxis5=PlotlyJS.attr(title="best path cost", type="log"),
        xaxis6=PlotlyJS.attr(title="time (s)"),
        yaxis6=PlotlyJS.attr(title="free base / ideal cancellation"),
        legend=PlotlyJS.attr(orientation="h", y=-0.18),
        margin=PlotlyJS.attr(l=70, r=35, t=70, b=85),
    )
    return PlotlyJS.Plot(traces, layout)
end

function _planner_summary_html(
    plan::RobotArmPlan,
    sim::ClothRobotArmSimulation,
    obstacles::Vector{RobotArmSphereObstacle},
    hypr::RobotArmHYPRResult,
    base_sim::CoupledBaseSimulation,
)
    ref_clearance = minimum(_obstacle_clearance_series(plan.ee_ref, obstacles))
    cloth_clearance = minimum(_obstacle_clearance_series(sim.end_effector_positions, obstacles))
    base_disp = maximum(vec(sqrt.(sum(base_sim.base_positions_m .^ 2; dims=1))))
    base_att = maximum(base_sim.base_attitude_error_rad)
    thrust_peak = maximum(vec(sqrt.(sum(base_sim.stabilization_thrust_n .^ 2; dims=1))))
    rw_peak = maximum(vec(sqrt.(sum(base_sim.stabilization_rw_torque_nm .^ 2; dims=1))))
    metric_names = [
        "IK final error (mm)",
        "Cloth max EE error (mm)",
        "Cloth final EE error (mm)",
        "Min reference clearance (mm)",
        "Min Cloth clearance (mm)",
        "Plan duration (s)",
        "Reference samples",
        "HYPR best cost",
        "HYPR link violations",
        "Max base displacement (mm)",
        "Max base attitude (mrad)",
        "Peak ideal thrust (mN)",
        "Peak ideal RW torque (mN-m)",
        "Max thrust axis saturation",
        "Max RW torque saturation",
    ]
    values = [
        1.0e3 * plan.final_error_m,
        1.0e3 * maximum(sim.tracking_error_m),
        1.0e3 * sim.tracking_error_m[end],
        1.0e3 * ref_clearance,
        1.0e3 * cloth_clearance,
        plan.t_ref_s[end],
        length(plan.t_ref_s),
        hypr.cost,
        hypr.components.violation_count,
        1.0e3 * base_disp,
        1.0e3 * base_att,
        1.0e3 * thrust_peak,
        1.0e3 * rw_peak,
        maximum(base_sim.thrust_saturation_ratio),
        maximum(base_sim.rw_saturation_ratio),
    ]
    return PlotlyJS.Plot(
        [PlotlyJS.bar(x=metric_names, y=values, marker=PlotlyJS.attr(color="rgb(42,86,125)"))],
        PlotlyJS.Layout(
            title="Robot Arm HYPR Planner Summary",
            yaxis=PlotlyJS.attr(title="value"),
            margin=PlotlyJS.attr(l=70, r=30, t=70, b=130),
        ),
    )
end

function _single_arm_case_spec()
    return (
        case_id=1,
        seed=740,
        q_start=[0.08, 0.95, -0.85],
        q_target_truth=[-0.12, -0.08, 0.06],
        obstacle_count=5,
        obstacle_radius_range_m=(0.018, 0.034),
        lateral_offset_m=0.010,
        x_offset_range_m=(-0.4, 0.1),
    )
end

function _arm_batch_case_label(case_id::Integer)
    return lpad(case_id, 3, '0')
end

function generate_robot_arm_batch_cases(; n_cases::Integer=50, seed::Integer=740)
    templates = [
        ([0.08, 0.95, -0.85], [-0.12, -0.08, 0.06], 5, (0.018, 0.034), 0.010, (-0.4, 0.1)),
        ([-0.25, 0.78, -0.62], [0.35, -0.12, 0.10], 5, (0.016, 0.032), 0.011, (-0.4, 0.1)),
        ([0.35, -0.08, 0.05], [-0.35, 0.88, -0.70], 6, (0.014, 0.030), 0.012, (-0.4, 0.1)),
        ([-0.45, 0.68, -0.55], [0.25, -0.18, 0.18], 5, (0.018, 0.036), 0.010, (-0.4, 0.1)),
        ([0.20, -0.05, 0.02], [-0.50, 1.02, -0.82], 6, (0.015, 0.031), 0.012, (-0.4, 0.1)),
        ([-0.10, 1.08, -0.95], [0.55, -0.16, 0.12], 5, (0.017, 0.033), 0.011, (-0.4, 0.1)),
        ([0.50, 0.72, -0.66], [-0.20, -0.10, 0.08], 6, (0.014, 0.030), 0.012, (-0.4, 0.1)),
        ([-0.35, -0.08, 0.04], [0.40, 0.92, -0.78], 5, (0.018, 0.035), 0.010, (-0.4, 0.1)),
    ]
    cases = NamedTuple[]
    for case_id in 1:n_cases
        q_start, q_target, obstacle_count, radius_range, lateral_offset, x_offset_range =
            templates[mod1(case_id, length(templates))]
        push!(
            cases,
            (
                case_id=case_id,
                seed=seed + case_id - 1,
                q_start=copy(q_start),
                q_target_truth=copy(q_target),
                obstacle_count=obstacle_count,
                obstacle_radius_range_m=radius_range,
                lateral_offset_m=lateral_offset,
                x_offset_range_m=x_offset_range,
            ),
        )
    end
    return cases
end

function _robot_arm_case_metrics(case, plan, sim, obstacles, hypr, base_sim; case_dir::AbstractString, motion_html::AbstractString="")
    ref_clearance = minimum(_obstacle_clearance_series(plan.ee_ref, obstacles))
    cloth_clearance = minimum(_obstacle_clearance_series(sim.end_effector_positions, obstacles))
    base_disp = maximum(vec(sqrt.(sum(base_sim.base_positions_m .^ 2; dims=1))))
    base_att = maximum(base_sim.base_attitude_error_rad)
    thrust_peak = maximum(vec(sqrt.(sum(base_sim.stabilization_thrust_n .^ 2; dims=1))))
    rw_peak = maximum(vec(sqrt.(sum(base_sim.stabilization_rw_torque_nm .^ 2; dims=1))))
    fuel_used = _stabilization_fuel_used_kg(base_sim)
    return (
        case_id=case.case_id,
        seed=case.seed,
        start_q1_rad=case.q_start[1],
        start_q2_rad=case.q_start[2],
        start_q3_rad=case.q_start[3],
        target_truth_q1_rad=case.q_target_truth[1],
        target_truth_q2_rad=case.q_target_truth[2],
        target_truth_q3_rad=case.q_target_truth[3],
        target_x_m=plan.target[1],
        target_y_m=plan.target[2],
        target_z_m=plan.target[3],
        success=(plan.final_error_m <= 5.0e-3 && maximum(sim.tracking_error_m) <= 0.035 && cloth_clearance >= -0.005) ? 1 : 0,
        obstacle_count=length(obstacles),
        planned_duration_s=plan.t_ref_s[end],
        reference_samples=length(plan.t_ref_s),
        fuel_used_kg=fuel_used,
        fuel_used_pct=_fuel_used_pct(fuel_used),
        propellant_mass_kg=RPO_3U_CUBESAT_PROP_MASS_KG,
        thrust_impulse_ns=_trapz(base_sim.t_s, vec(sum(abs.(base_sim.stabilization_thrust_n); dims=1))),
        torque_impulse_nms=_trapz(base_sim.t_s, vec(sqrt.(sum(base_sim.stabilization_rw_torque_nm .^ 2; dims=1)))),
        peak_ideal_thrust_n=thrust_peak,
        peak_ideal_rw_torque_nm=rw_peak,
        max_thrust_saturation_ratio=maximum(base_sim.thrust_saturation_ratio),
        max_rw_saturation_ratio=maximum(base_sim.rw_saturation_ratio),
        thrust_saturation_fraction=count(x -> x >= 1.0, base_sim.thrust_saturation_ratio) / max(length(base_sim.thrust_saturation_ratio), 1),
        rw_saturation_fraction=count(x -> x >= 1.0, base_sim.rw_saturation_ratio) / max(length(base_sim.rw_saturation_ratio), 1),
        ik_final_error_m=plan.final_error_m,
        cloth_tracking_final_m=sim.tracking_error_m[end],
        cloth_tracking_mean_m=sum(sim.tracking_error_m) / max(length(sim.tracking_error_m), 1),
        cloth_tracking_max_m=maximum(sim.tracking_error_m),
        min_reference_clearance_m=ref_clearance,
        min_cloth_clearance_m=cloth_clearance,
        max_base_displacement_m=base_disp,
        max_base_attitude_error_rad=base_att,
        hypr_best_cost=hypr.cost,
        hypr_initial_cost=isempty(hypr.cost_history) ? hypr.cost : first(hypr.cost_history),
        hypr_iterations=length(hypr.cost_history),
        hypr_early_stopped=hypr.early_stopped ? 1 : 0,
        hypr_link_violations=hypr.components.violation_count,
        hypr_retime_base_force_ratio=hypr.components.retime_base_force_ratio,
        motion_3d_html=abspath(motion_html),
        case_dir=abspath(case_dir),
    )
end

function _run_robot_arm_case(case; out_dir::AbstractString, save_plots::Bool)
    mkpath(out_dir)
    arm = default_cloth_arm_model(
        mount_offset_body=(0.15, 0.0, 0.0),
    )
    base_pose = ClothArmBasePose([0.0, 0.0, 0.0])
    q_start = Float64.(case.q_start)

    # Use a known reachable FK pose as the synthetic inspection target. This keeps
    # the example deterministic and makes IK/planning errors meaningful.
    q_target_truth = Float64.(case.q_target_truth)
    target = cloth_fk(arm, base_pose, q_target_truth).end_effector_position

    planner_cfg = RobotArmPlannerConfig(
        dt_s=0.05,
        duration_s=2.0,
        ik_tol_m=5.0e-4,
        ik_max_iters=150,
        ik_damping=1.0e-3,
    )
    nominal_plan = plan_robot_arm_motion(
        arm,
        base_pose,
        q_start,
        target;
        config=planner_cfg,
        planner=:cloth_quintic,
    )
    obstacles = _random_sphere_obstacles_on_nominal_path(
        nominal_plan;
        rng=MersenneTwister(case.seed),
        count=case.obstacle_count,
        radius_range_m=case.obstacle_radius_range_m,
        lateral_offset_m=case.lateral_offset_m,
        x_offset_range_m=case.x_offset_range_m,
        endpoint_safe_distance_m=0.006,
        endpoint_exclusion_fraction=0.22,
    )
    _CURRENT_OBSTACLES[] = obstacles
    hypr_cfg = RobotArmHYPRConfig(
        n_waypoints=5,
        n_particles=120,
        n_iters=70,
        n_samples=80,
        safe_distance_m=0.006,
        w_len=1.0,
        w_smooth=0.30,
        w_obs=800.0,
        spread_scale=0.24,
        cull_start_iter=8,
        early_stopping_min_iters=20,
        early_stopping_patience=14,
        early_stopping_require_feasible=false,
        retime_enable=true,
        retime_max_joint_velocity_rad_s=0.55,
        retime_max_joint_acceleration_rad_s2=1.2,
        retime_reaction_time_s=0.0,
        retime_reaction_gain=1.0,
        retime_max_scale=100.0,
        retime_max_base_force_n=RPO_3U_CUBESAT_MAX_AXIS_THRUST_N,
        retime_max_base_torque_nm=RPO_3U_CUBESAT_RW_MAX_TORQUE_NM,
        retime_base_wrench_model_gain=1.0,
        retime_base_wrench_margin=1.15,
        retime_base_wrench_iters=6,
        retime_cloth_physics_enable=true,
        retime_cloth_dt_s=0.025,
        retime_cloth_k_translation_n_m=8.0e3,
        retime_cloth_c_translation_n_s_m=40.0,
        retime_cloth_k_rotation_n_m_rad=25.0,
        retime_cloth_c_rotation_n_m_s_rad=1.2,
    )
    hypr = plan_robot_arm_motion_hypr(
        arm,
        base_pose,
        q_start,
        target;
        planner_config=planner_cfg,
        hypr_config=hypr_cfg,
        obstacles=obstacles,
        rng=MersenneTwister(case.seed + 1),
    )
    plan = hypr.plan
    sim_duration_s = plan.t_ref_s[end]

    sim = simulate_cloth_robot_arm_plan(
        plan;
        dt_s=0.025,
        duration_s=sim_duration_s,
        integrator=:implicit_midpoint,
        k_translation_n_m=8.0e3,
        c_translation_n_s_m=40.0,
        k_rotation_n_m_rad=25.0,
        c_rotation_n_m_s_rad=1.2,
    )
    base_sim = simulate_coupled_base_motion(
        plan;
        arm_sim=sim,
        dt_s=0.001,
        duration_s=sim_duration_s,
        k_translation_n_m=8.0e3,
        c_translation_n_s_m=40.0,
        k_rotation_n_m_rad=25.0,
        c_rotation_n_m_s_rad=1.2,
    )

    joint_csv = _save_csv(joinpath(out_dir, "arm_joint_plan.csv"), _joint_dataframe(plan))
    ee_csv = _save_csv(joinpath(out_dir, "arm_end_effector_tracking.csv"), _ee_dataframe(plan, sim))
    obstacle_csv = _save_csv(joinpath(out_dir, "arm_sphere_obstacles.csv"), _obstacle_dataframe(obstacles))
    base_csv = _save_csv(joinpath(out_dir, "arm_base_motion_and_stabilization.csv"), _base_dataframe(base_sim))

    motion_html = _save_html_plot(joinpath(out_dir, "arm_motion_3d.html"), _plot_arm_motion_html(plan, sim, obstacles, base_sim))
    metrics_html = ""
    summary_html = ""
    if save_plots
        _save_plot(joinpath(out_dir, "arm_joint_plan.png"), _plot_joint_history(plan))
        _save_plot(joinpath(out_dir, "arm_end_effector_path.png"), _plot_ee_path(plan, sim, obstacles))
        _save_plot(joinpath(out_dir, "arm_tracking_error.png"), _plot_tracking_error(sim))
        _save_plot(joinpath(out_dir, "arm_base_motion.png"), _plot_base_motion(base_sim))
        _save_plot(joinpath(out_dir, "arm_base_stabilization_demand.png"), _plot_stabilization_demand(base_sim))
        metrics_html = _save_html_plot(joinpath(out_dir, "arm_hypr_planner_metrics.html"), _planner_metrics_html(plan, sim, obstacles, hypr, base_sim))
        summary_html = _save_html_plot(joinpath(out_dir, "arm_hypr_planner_summary.html"), _planner_summary_html(plan, sim, obstacles, hypr, base_sim))
    end

    max_error_mm = 1.0e3 * maximum(sim.tracking_error_m)
    final_error_mm = 1.0e3 * sim.tracking_error_m[end]
    min_ref_clearance_mm = 1.0e3 * minimum(_obstacle_clearance_series(plan.ee_ref, obstacles))
    min_cloth_clearance_mm = 1.0e3 * minimum(_obstacle_clearance_series(sim.end_effector_positions, obstacles))
    max_base_disp_mm = 1.0e3 * maximum(vec(sqrt.(sum(base_sim.base_positions_m .^ 2; dims=1))))
    max_base_att_mrad = 1.0e3 * maximum(base_sim.base_attitude_error_rad)
    peak_thrust_mn = 1.0e3 * maximum(vec(sqrt.(sum(base_sim.stabilization_thrust_n .^ 2; dims=1))))
    peak_rw_mnm = 1.0e3 * maximum(vec(sqrt.(sum(base_sim.stabilization_rw_torque_nm .^ 2; dims=1))))
    fuel_used = _stabilization_fuel_used_kg(base_sim)
    println("  case $(_arm_batch_case_label(case.case_id)) complete: ",
        "fuel=", @sprintf("%.6f kg", fuel_used),
        ", max EE error=", @sprintf("%.3f mm", max_error_mm),
        ", min clearance=", @sprintf("%.3f mm", min_cloth_clearance_mm),
        ", duration=", @sprintf("%.3f s", plan.t_ref_s[end]))

    return (
        case=case,
        plan=plan,
        nominal_plan=nominal_plan,
        hypr=hypr,
        simulation=sim,
        base_simulation=base_sim,
        obstacles=obstacles,
        joint_csv=joint_csv,
        ee_csv=ee_csv,
        obstacle_csv=obstacle_csv,
        base_csv=base_csv,
        motion_html=motion_html,
        metrics_html=metrics_html,
        summary_html=summary_html,
        metrics=_robot_arm_case_metrics(case, plan, sim, obstacles, hypr, base_sim; case_dir=out_dir, motion_html=motion_html),
    )
end

function run_robot_arm_planner_cloth_demo()
    result = _run_robot_arm_case(_single_arm_case_spec(); out_dir=OUT_DIR, save_plots=true)
    println()
    println("Robot arm planner + Cloth dynamics demo complete.")
    println("  planner              = ", result.plan.planner)
    println("  sphere obstacles      = ", length(result.obstacles))
    println("  fuel used             = ", @sprintf("%.6f kg (%.3f%%)", result.metrics.fuel_used_kg, result.metrics.fuel_used_pct))
    println("  HYPR best cost        = ", @sprintf("%.6g", result.hypr.cost))
    println("  HYPR link violations  = ", result.hypr.components.violation_count)
    println("  retimed duration      = ", @sprintf("%.3f s", result.plan.t_ref_s[end]))
    println("  max thrust saturation = ", @sprintf("%.3f x", maximum(result.base_simulation.thrust_saturation_ratio)))
    println("  max RW saturation     = ", @sprintf("%.3f x", maximum(result.base_simulation.rw_saturation_ratio)))
    println("  output directory      = ", abspath(OUT_DIR))
    return result
end

function _arm_batch_metrics_plot(df::DataFrame)
    specs = [
        (:fuel_used_pct, "Fuel used (%)"),
        (:cloth_tracking_max_m, "Max tracking error (m)"),
        (:min_cloth_clearance_m, "Min cloth clearance (m)"),
        (:max_thrust_saturation_ratio, "Max thrust sat. (x)"),
        (:max_rw_saturation_ratio, "Max RW sat. (x)"),
        (:planned_duration_s, "Plan duration (s)"),
    ]
    traces = PlotlyJS.GenericTrace[]
    for (idx, (metric, label)) in enumerate(specs)
        push!(traces, PlotlyJS.bar(
            x=df.case_id,
            y=Float64.(df[!, metric]),
            name=label,
            xaxis=idx == 1 ? "x" : "x$idx",
            yaxis=idx == 1 ? "y" : "y$idx",
        ))
    end
    return PlotlyJS.Plot(
        traces,
        PlotlyJS.Layout(
            title="Robot Arm Batch Metrics",
            grid=PlotlyJS.attr(rows=3, columns=2, pattern="independent"),
            xaxis=PlotlyJS.attr(title="case"),
            yaxis=PlotlyJS.attr(title="fuel used (%)"),
            xaxis2=PlotlyJS.attr(title="case"),
            yaxis2=PlotlyJS.attr(title="tracking error (m)"),
            xaxis3=PlotlyJS.attr(title="case"),
            yaxis3=PlotlyJS.attr(title="clearance (m)"),
            xaxis4=PlotlyJS.attr(title="case"),
            yaxis4=PlotlyJS.attr(title="thrust saturation"),
            xaxis5=PlotlyJS.attr(title="case"),
            yaxis5=PlotlyJS.attr(title="RW saturation"),
            xaxis6=PlotlyJS.attr(title="case"),
            yaxis6=PlotlyJS.attr(title="duration (s)"),
            showlegend=false,
            margin=PlotlyJS.attr(l=70, r=30, t=70, b=55),
        ),
    )
end

function _arm_batch_tracking_plot(results)
    traces = PlotlyJS.GenericTrace[]
    for result in results
        push!(traces, PlotlyJS.scatter(
            x=result.simulation.trajectory.t_s,
            y=result.simulation.tracking_error_m,
            mode="lines",
            name="case $(_arm_batch_case_label(result.case.case_id))",
        ))
    end
    return PlotlyJS.Plot(
        traces,
        PlotlyJS.Layout(
            title="Robot Arm Batch Cloth Tracking Error",
            xaxis_title="time (s)",
            yaxis_title="tracking error (m)",
        ),
    )
end

function run_robot_arm_planner_cloth_batch(;
    n_cases::Integer=_env_int("SPACEAGORA_ARM_BATCH_CASES", get(ENV, "SPACEAGORA_EXAMPLE_SMOKE", "0") == "1" ? 2 : 50),
    seed::Integer=_env_int("SPACEAGORA_ARM_BATCH_SEED", 740),
    save_case_plots::Bool=_env_bool("SPACEAGORA_ARM_BATCH_SAVE_CASE_PLOTS", false),
)
    batch_dir = joinpath(REPO_ROOT, "output", "robot_arm_planner_cloth_batch")
    mkpath(batch_dir)
    cases = generate_robot_arm_batch_cases(n_cases=n_cases, seed=seed)
    results = NamedTuple[]
    rows = NamedTuple[]
    println("Running robot arm batch with $(length(cases)) cases.")
    for case in cases
        label = _arm_batch_case_label(case.case_id)
        println("Starting arm case $label")
        result = _run_robot_arm_case(case; out_dir=joinpath(batch_dir, "case_$label"), save_plots=save_case_plots)
        push!(results, result)
        push!(rows, result.metrics)
    end
    summary_df = DataFrame(rows)
    cases_df = DataFrame(cases)
    summary_csv = _save_csv(joinpath(batch_dir, "arm_batch_metrics_summary.csv"), summary_df)
    cases_csv = _save_csv(joinpath(batch_dir, "arm_batch_case_manifest.csv"), cases_df)
    metrics_html = _save_html_plot(joinpath(batch_dir, "arm_batch_metrics.html"), _arm_batch_metrics_plot(summary_df))
    tracking_html = _save_html_plot(joinpath(batch_dir, "arm_batch_tracking_error.html"), _arm_batch_tracking_plot(results))

    success_pct = 100.0 * sum(summary_df.success) / max(nrow(summary_df), 1)
    println()
    println("Robot arm batch complete.")
    println("  cases                 = ", nrow(summary_df))
    println("  success               = ", @sprintf("%.1f%%", success_pct))
    println("  mean fuel used        = ", @sprintf("%.6f kg", sum(summary_df.fuel_used_kg) / max(nrow(summary_df), 1)))
    println("  max fuel used         = ", @sprintf("%.6f kg", maximum(summary_df.fuel_used_kg)))
    println("  max tracking error    = ", @sprintf("%.3f mm", 1.0e3 * maximum(summary_df.cloth_tracking_max_m)))
    println("  min cloth clearance   = ", @sprintf("%.3f mm", 1.0e3 * minimum(summary_df.min_cloth_clearance_m)))
    println("  summary CSV           = ", abspath(summary_csv))
    println("  case manifest CSV     = ", abspath(cases_csv))
    println("  metrics HTML          = ", abspath(metrics_html))
    println("  tracking HTML         = ", abspath(tracking_html))

    return (
        results=results,
        summary=summary_df,
        cases=cases_df,
        summary_csv=summary_csv,
        cases_csv=cases_csv,
        metrics_html=metrics_html,
        tracking_html=tracking_html,
    )
end

if _env_bool("SPACEAGORA_ARM_SINGLE", false)
    run_robot_arm_planner_cloth_demo()
else
    run_robot_arm_planner_cloth_batch()
end
