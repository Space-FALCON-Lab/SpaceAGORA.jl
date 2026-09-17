include(joinpath(@__DIR__, "common.jl"))
using CSV
using DataFrames
using PlotlyJS
using Random
using StaticArrays
using LinearAlgebra

"""
    build_rpo_cubesat_mpc_demo(; mission_time=180.0)

Build a complete two-spacecraft RPO scenario: a passive station target and a
six-axis CubeSat chaser using HYPR/PSO guidance and LQ-MPC tracking.
"""
function build_rpo_cubesat_mpc_demo(;
    mission_time=180.0,
    results_directory=joinpath(REPO_ROOT, "output", "rpo_single_case"),
    seed::Integer=741,
    start_rtn=SVector{3, Float64}(-8.0, -4.0, 2.0),
    goal_rtn=SVector{3, Float64}(5.0, 0.0, 0.0),
    pso_n_particles::Integer=120,
    pso_n_iters::Integer=35,
    pso_config=nothing,
    pso_configurator=nothing,
    n_station_points::Integer=10000,
    station_geometry_seed::Integer=seed,
    data_rate_s::Real=10.0,
    pso_iteration_runtime_limit_s=nothing,
    pso_iteration_callback=nothing,
    verbose::Bool=true,
)
    start_rtn = SVector{3, Float64}(start_rtn)
    goal_rtn = SVector{3, Float64}(goal_rtn)
    planet = Earth("", SPICE_PATH)
    orbit_radius = planet.Rp_e + 420e3
    μ = planet.μ
    n = sqrt(μ / orbit_radius^3)
    r_station_ii = SVector{3, Float64}(orbit_radius, 0.0, 0.0)
    v_station_ii = SVector{3, Float64}(0.0, sqrt(μ / orbit_radius), 0.0)
    r0_rel_rtn = start_rtn
    v0_rel_rtn = SVector{3, Float64}(0.0, 0.0, 0.0)
    r_chaser_ii, v_chaser_ii = SimulationModel.GuidanceHooks.rtn_to_inertial_relative_state(
        r0_rel_rtn,
        v0_rel_rtn,
        r_station_ii,
        v_station_ii,
    )

    station_points = SpaceAGORA.load_rpo_station_cad_pointcloud(:gateway; n_points=n_station_points, rng=MersenneTwister(station_geometry_seed))
    station_geometry = RPOStationGeometry(station_points; keepout_radius_m=0.25, name="gateway_core")
    chaser_geometry = RPOCubeSatGeometry(dims_m=(0.1, 0.1, 0.3))
    geometry = RPOReferenceGeometry(station_geometry; chaser=chaser_geometry)

    q_identity = SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)
    station_root = Link(
        root=true,
        m=500.0,
        dims=MVector{3, Float64}(4.0, 2.0, 2.0),
        ref_area=8.0,
    )
    station = SpacecraftModel(
        joints=Joint[],
        links=[station_root],
        root=station_root,
        prop_mass=0.0,
        inertia_tensor=station_root.inertia,
        initial_condition=CartesianInitialCondition(r_station_ii, v_station_ii; q=q_identity),
        id=201,
    )

    chaser_root = Link(
        root=true,
        m=5.0,
        dims=MVector{3, Float64}(0.1, 0.1, 0.3),
        ref_area=0.03,
    )
    chaser = SpacecraftModel(
        joints=Joint[],
        links=[chaser_root],
        root=chaser_root,
        prop_mass=0.2,
        inertia_tensor=chaser_root.inertia,
        initial_condition=CartesianInitialCondition(r_chaser_ii, v_chaser_ii; q=q_identity),
        id=101,
    )

    thrusters = SixAxisThrusterModel(
        max_thrust_n=SVector{6, Float64}(fill(0.05, 6)),
        isp_s=SVector{6, Float64}(fill(60.0, 6)),
    )
    chaser_initial_mass_kg = chaser_root.m + chaser.prop_mass
    max_axis_accel_mps2 = minimum(thrusters.max_thrust_n ./ chaser_initial_mass_kg)
    retime_accel_mps2 = 0.5 * max_axis_accel_mps2

    if pso_config !== nothing && pso_configurator !== nothing
        throw(ArgumentError("Pass either pso_config or pso_configurator, not both."))
    end
    default_pso_configurator = RPOPSOConfigurator(
        swarm=RPOPSOSwarmSettings(
            n_waypoints=5,
            n_particles=Int(pso_n_particles),
            n_iters=Int(pso_n_iters),
            curve_type=:bezier,
            sample_ds_m=0.05,
        ),
        adaptive=RPOPSOAdaptiveSettings(
            allow_downscale=true,
            n_particles_min=80,
            n_particles_max=160,
            n_iters_min=8,
            n_iters_max=35,
        ),
        early_stopping=RPOPSOEarlyStoppingSettings(
            enabled=true,
            patience=8,
            min_iters=12,
            min_rel_improvement=1.0e-4,
        ),
        objective=RPOPSOObjectiveSettings(
            cost_ref_distance_m=20.0,
            mass_kg=chaser_initial_mass_kg,
            tf_s=Float64(mission_time),
        ),
        retiming=RPOPSORetimingSettings(
            dt_s=0.1,
            a_max_mps2=retime_accel_mps2,
        ),
    )
    pso_cfg = if pso_config !== nothing
        rpo_pso_config(pso_config)
    elseif pso_configurator !== nothing
        rpo_pso_config(
            pso_configurator;
            mass_kg=chaser_initial_mass_kg,
            tf_s=Float64(mission_time),
            retime_a_max_mps2=retime_accel_mps2,
        )
    else
        rpo_pso_config(default_pso_configurator)
    end
    if pso_iteration_runtime_limit_s !== nothing
        pso_cfg = rpo_pso_config(pso_cfg; iteration_runtime_limit_s=Float64(pso_iteration_runtime_limit_s))
    end

    plan_buffer = RPOPlanBuffer()
    plan_result = rpo_pso_plan_path(
        r0_rel_rtn,
        goal_rtn,
        geometry,
        pso_cfg;
        safe_distance_m=0.1,
        rng=MersenneTwister(seed),
        iteration_callback=pso_iteration_callback,
    )
    t_ref, r_ref, v_ref = rpo_reference_from_path(
        plan_result.path,
        geometry,
        plan_result.config;
        safe_distance_m=0.1,
    )
    simulation_time_s = max(Float64(mission_time), t_ref[end] + 20.0)
    initial_plan = RPOPlan(
        valid=true,
        t_ref_s=t_ref,
        r_ref_rtn=r_ref,
        v_ref_rtn=v_ref,
        path_rtn=plan_result.path,
        cost=plan_result.cost,
        diagnostics=(
            components=plan_result.components,
            adaptive=plan_result.adaptive,
            refinement_improved=plan_result.refinement_improved,
            iteration_timed_out=plan_result.iteration_timed_out,
            iteration_timeout_iter=plan_result.iteration_timeout_iter,
            iteration_timeout_phase=plan_result.iteration_timeout_phase,
            iteration_timeout_events=plan_result.iteration_timeout_events,
            cost_history=plan_result.cost_history,
            planned_at_s=0.0,
        ),
    )
    update_rpo_plan_buffer!(plan_buffer, initial_plan, 0.0)
    guidance = RPOGuidanceModel(
        chaser_idx=1,
        target_idx=2,
        goal_rtn=goal_rtn,
        geometry=geometry,
        plan_buffer=plan_buffer,
        pso_config=pso_cfg,
        safe_distance_m=0.1,
    )

    Q = Diagonal([20.0, 20.0, 20.0, 2.0, 2.0, 2.0])
    R = 0.1I(3)
    Qf = 10.0 .* Matrix(Q)
    controller = init_rpo_lqmpc(
        n,
        pso_cfg.retime_dt_s,
        Matrix(Q),
        Matrix(R),
        Qf,
        12;
        u_min=fill(-max_axis_accel_mps2, 3),
        u_max=fill(max_axis_accel_mps2, 3),
    )
    control = RPOMPCControlModel(
        chaser_idx=1,
        target_idx=2,
        thrusters=thrusters,
        controller=controller,
        plan_buffer=plan_buffer,
        control_dt_s=pso_cfg.retime_dt_s,
        attitude_kp=0.02,
        rate_kd=0.08,
        max_rw_torque_nm=0.002,
    )

    args = SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=true,
            verbose=verbose,
            generate_plots=false,
            results_directory=results_directory,
            normalize=false,
        ),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=false,
            number_of_orbits=1,
            mission_time=simulation_time_s,
            orientation_sim=true,
            num_steps_to_save=400,
            data_rate=Float64(data_rate_s),
        ),
        environment_model=EnvironmentModel(
            planet=planet,
            EI=120.0,
            density_model=NoAtmosphereModel(),
            ephemerides_model=SimpleEphemeridesModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false,
            wind=false,
        ),
        dynamics_model=DynamicsModel([chaser, station], (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel(guidance_effectors=(guidance,), guidance_rates=[pso_cfg.retime_dt_s]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(control,), control_rates=[pso_cfg.retime_dt_s]),
        initial_time=InitialTime(
            year=2026,
            month=1,
            day=1,
            hour=0,
            minute=0,
            second=0.0,
        ),
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            reltol_quaternion=1e-8,
            abstol_quaternion=1e-8,
            dt_max_orbit=0.05,
        ),
    )

    return (
        args=args,
        geometry=geometry,
        guidance=guidance,
        control=control,
        pso_config=pso_cfg,
        plan_buffer=plan_buffer,
        initial_relative_state_rtn=vcat(r0_rel_rtn, v0_rel_rtn),
        goal_rtn=guidance.goal_rtn,
        initial_plan=initial_plan,
        plan_result=plan_result,
        seed=Int(seed),
    )
end

function _interp_ref_columns(values::Matrix{Float64}, t_ref::Vector{Float64}, t::Real)
    size(values, 2) == 0 && return zeros(3)
    tq = Float64(t)
    tq <= t_ref[1] && return copy(values[:, 1])
    tq >= t_ref[end] && return copy(values[:, end])
    hi = searchsortedfirst(t_ref, tq)
    lo = max(hi - 1, 1)
    denom = max(t_ref[hi] - t_ref[lo], 1.0e-9)
    α = (tq - t_ref[lo]) / denom
    return (1.0 - α) .* values[:, lo] .+ α .* values[:, hi]
end

function _earth_surface_trace(radius_m; n_lon=48, n_lat=24)
    lon = range(0.0, 2π; length=n_lon)
    lat = range(-π / 2, π / 2; length=n_lat)
    x = [radius_m * cos(φ) * cos(λ) for φ in lat, λ in lon]
    y = [radius_m * cos(φ) * sin(λ) for φ in lat, λ in lon]
    z = [radius_m * sin(φ) for φ in lat, λ in lon]
    return surface(
        x=x ./ 1000.0,
        y=y ./ 1000.0,
        z=z ./ 1000.0,
        opacity=0.35,
        colorscale=[[0.0, "rgb(40,90,170)"], [1.0, "rgb(70,130,190)"]],
        showscale=false,
        name="Earth",
    )
end

function _rpo_postprocess(csv_path::AbstractString, demo)
    df = CSV.read(csv_path, DataFrame)
    n_samples = nrow(df)
    actual_rtn = zeros(3, n_samples)
    ref_rtn = zeros(3, n_samples)
    tracking_error = zeros(n_samples)
    for k in 1:n_samples
        r_chaser = SVector{3, Float64}(df.sc1_pos_1[k], df.sc1_pos_2[k], df.sc1_pos_3[k])
        v_chaser = SVector{3, Float64}(df.sc1_vel_1[k], df.sc1_vel_2[k], df.sc1_vel_3[k])
        r_target = SVector{3, Float64}(df.sc2_pos_1[k], df.sc2_pos_2[k], df.sc2_pos_3[k])
        v_target = SVector{3, Float64}(df.sc2_vel_1[k], df.sc2_vel_2[k], df.sc2_vel_3[k])
        x_rel = SimulationModel.GuidanceHooks.inertial_to_rtn_relative_state(
            r_chaser,
            v_chaser,
            r_target,
            v_target,
        )
        actual_rtn[:, k] .= x_rel[1:3]
        ref_rtn[:, k] .= _interp_ref_columns(demo.initial_plan.r_ref_rtn, demo.initial_plan.t_ref_s, df.time[k])
        tracking_error[k] = norm(actual_rtn[:, k] - ref_rtn[:, k])
    end
    return df, actual_rtn, ref_rtn, tracking_error
end

function _save_plot(plot_obj, output_dir::AbstractString, filename::AbstractString)
    mkpath(output_dir)
    path = joinpath(output_dir, filename)
    PlotlyJS.savefig(plot_obj, path)
    return path
end

const _STATION_MESH_TRACE_CACHE = Ref{Any}(nothing)

function _build_station_mesh_trace()
    triangles = SpaceAGORA.load_rpo_station_cad_triangles(:gateway)
    ntri = size(triangles, 2) ÷ 3
    tri_idxs = collect(1:ntri)
    nfaces = length(tri_idxs)
    vertices = zeros(Float64, 3, 3 * nfaces)
    i = zeros(Int, nfaces)
    j = zeros(Int, nfaces)
    k = zeros(Int, nfaces)
    for (face_idx, tri_idx) in enumerate(tri_idxs)
        src = 3 * (tri_idx - 1)
        dst = 3 * (face_idx - 1)
        vertices[:, dst + 1] .= triangles[:, src + 1]
        vertices[:, dst + 2] .= triangles[:, src + 2]
        vertices[:, dst + 3] .= triangles[:, src + 3]
        i[face_idx] = dst
        j[face_idx] = dst + 1
        k[face_idx] = dst + 2
    end
    return mesh3d(
        x=vertices[1, :],
        y=vertices[2, :],
        z=vertices[3, :],
        i=i,
        j=j,
        k=k,
        color="rgb(150,160,170)",
        opacity=0.55,
        flatshading=true,
        name="Gateway mesh",
        hoverinfo="skip",
    )
end

function _station_mesh_trace(; refresh::Bool=false)
    if refresh || _STATION_MESH_TRACE_CACHE[] === nothing
        _STATION_MESH_TRACE_CACHE[] = _build_station_mesh_trace()
    end
    return deepcopy(_STATION_MESH_TRACE_CACHE[])
end

function _cuboid_mesh_trace(centers::Matrix{Float64}, half_extents; max_cubes::Integer=12)
    n = size(centers, 2)
    n > 0 || return mesh3d(x=Float64[], y=Float64[], z=Float64[], i=Int[], j=Int[], k=Int[], name="CubeSat geometry")
    idxs = unique(round.(Int, range(1, n; length=min(max_cubes, n))))
    hx, hy, hz = Tuple(Float64.(half_extents))
    offsets = [
        (-hx, -hy, -hz), ( hx, -hy, -hz), ( hx,  hy, -hz), (-hx,  hy, -hz),
        (-hx, -hy,  hz), ( hx, -hy,  hz), ( hx,  hy,  hz), (-hx,  hy,  hz),
    ]
    faces = [
        (1, 2, 3), (1, 3, 4),
        (5, 7, 6), (5, 8, 7),
        (1, 5, 6), (1, 6, 2),
        (2, 6, 7), (2, 7, 3),
        (3, 7, 8), (3, 8, 4),
        (4, 8, 5), (4, 5, 1),
    ]
    vertices = zeros(Float64, 3, 8 * length(idxs))
    face_i = Int[]
    face_j = Int[]
    face_k = Int[]
    for (cube_idx, sample_idx) in enumerate(idxs)
        center = centers[:, sample_idx]
        vbase = 8 * (cube_idx - 1)
        for (offset_idx, offset) in enumerate(offsets)
            vertices[:, vbase + offset_idx] .= center .+ collect(offset)
        end
        for face in faces
            push!(face_i, vbase + face[1] - 1)
            push!(face_j, vbase + face[2] - 1)
            push!(face_k, vbase + face[3] - 1)
        end
    end
    return mesh3d(
        x=vertices[1, :],
        y=vertices[2, :],
        z=vertices[3, :],
        i=face_i,
        j=face_j,
        k=face_k,
        color="rgb(38,120,180)",
        opacity=0.38,
        flatshading=true,
        name="CubeSat geometry samples",
        hoverinfo="skip",
    )
end

function _save_rpo_single_case_plots(csv_path::AbstractString, demo)
    df, actual_rtn, ref_rtn, tracking_error = _rpo_postprocess(csv_path, demo)
    plot_dir = joinpath(dirname(csv_path), "plots")
    mkpath(plot_dir)

    earth_plot = Plot(
        [
            _earth_surface_trace(demo.args.environment_model.planet.Rp_e),
            scatter3d(
                x=df.sc1_pos_1 ./ 1000.0,
                y=df.sc1_pos_2 ./ 1000.0,
                z=df.sc1_pos_3 ./ 1000.0,
                mode="lines+markers",
                name="CubeSat",
            ),
            scatter3d(
                x=df.sc2_pos_1 ./ 1000.0,
                y=df.sc2_pos_2 ./ 1000.0,
                z=df.sc2_pos_3 ./ 1000.0,
                mode="lines+markers",
                name="station",
            ),
        ],
        Layout(
            title="RPO Earth-Centered Trajectory",
            scene=attr(aspectmode="data", xaxis_title="x (km)", yaxis_title="y (km)", zaxis_title="z (km)"),
        ),
    )

    pso_path = SM.GuidanceHooks.rpo_sample_path(
        demo.initial_plan.path_rtn,
        demo.plan_result.config.sample_ds_m;
        curve_type=demo.plan_result.config.curve_type,
    )
    station_plot = Plot(
        [
            _station_mesh_trace(),
            scatter3d(
                x=pso_path[1, :],
                y=pso_path[2, :],
                z=pso_path[3, :],
                mode="lines",
                line=attr(width=6, color="rgb(210,85,45)"),
                name="PSO optimizer path",
            ),
            scatter3d(
                x=demo.initial_plan.r_ref_rtn[1, :],
                y=demo.initial_plan.r_ref_rtn[2, :],
                z=demo.initial_plan.r_ref_rtn[3, :],
                mode="lines",
                line=attr(dash="dot"),
                name="retimed reference",
            ),
            _cuboid_mesh_trace(actual_rtn, demo.geometry.chaser.half_extents_body),
            scatter3d(
                x=[demo.initial_relative_state_rtn[1]],
                y=[demo.initial_relative_state_rtn[2]],
                z=[demo.initial_relative_state_rtn[3]],
                mode="markers",
                marker=attr(size=7, color="rgb(45,150,90)", symbol="circle"),
                name="start",
            ),
            scatter3d(
                x=[demo.goal_rtn[1]],
                y=[demo.goal_rtn[2]],
                z=[demo.goal_rtn[3]],
                mode="markers",
                marker=attr(size=8, color="rgb(190,60,55)", symbol="diamond"),
                name="goal",
            ),
            scatter3d(
                x=actual_rtn[1, :],
                y=actual_rtn[2, :],
                z=actual_rtn[3, :],
                mode="lines",
                line=attr(width=4, color="rgb(20,90,150)"),
                name="actual CubeSat center",
            ),
        ],
        Layout(
            title="Station-Centered RPO Path and Tracking",
            scene=attr(aspectmode="data", xaxis_title="radial (m)", yaxis_title="along-track (m)", zaxis_title="cross-track (m)"),
        ),
    )

    tracking_plot = Plot(
        [
            scatter(x=df.time, y=tracking_error, mode="lines+markers", name="position tracking error"),
            scatter(x=df.time, y=actual_rtn[1, :] .- ref_rtn[1, :], mode="lines", name="radial error"),
            scatter(x=df.time, y=actual_rtn[2, :] .- ref_rtn[2, :], mode="lines", name="along-track error"),
            scatter(x=df.time, y=actual_rtn[3, :] .- ref_rtn[3, :], mode="lines", name="cross-track error"),
        ],
        Layout(title="RPO Tracking Error", xaxis_title="time (s)", yaxis_title="error (m)"),
    )

    cost_history = collect(get(demo.initial_plan.diagnostics, :cost_history, Float64[]))
    pso_cost_plot = Plot(
        [scatter(x=1:length(cost_history), y=cost_history, mode="lines+markers", name="best PSO cost")],
        Layout(title="PSO Cost History", xaxis_title="iteration", yaxis_title="normalized objective"),
    )

    components = demo.initial_plan.diagnostics.components
    component_plot = Plot(
        [bar(
            x=["length", "obstacle", "fuel"],
            y=[components.J_len_norm^2, components.J_obs, components.J_fuel_norm^2],
            name="normalized terms",
        )],
        Layout(title="PSO Normalized Cost Terms", yaxis_title="term value"),
    )

    mass_plot = Plot(
        [scatter(x=df.time, y=df.sc1_mass, mode="lines+markers", name="CubeSat mass")],
        Layout(title="CubeSat Mass History", xaxis_title="time (s)", yaxis_title="mass (kg)"),
    )

    paths = [
        _save_plot(earth_plot, plot_dir, "earth_centered_trajectory.html"),
        _save_plot(station_plot, plot_dir, "station_centered_tracking.html"),
        _save_plot(tracking_plot, plot_dir, "tracking_error.html"),
        _save_plot(pso_cost_plot, plot_dir, "pso_cost_history.html"),
        _save_plot(component_plot, plot_dir, "pso_cost_components.html"),
        _save_plot(mass_plot, plot_dir, "cubesat_mass_history.html"),
    ]
    return paths
end

function run_rpo_cubesat_mpc_demo(; mission_time=180.0, kwargs...)
    demo = build_rpo_cubesat_mpc_demo(; mission_time=mission_time, kwargs...)
    csv_path = run_and_report(demo.args)
    plot_paths = csv_path === nothing ? String[] : _save_rpo_single_case_plots(csv_path, demo)
    return merge(demo, (csv_path=csv_path, plot_paths=plot_paths))
end

if abspath(PROGRAM_FILE) == @__FILE__
    demo = run_rpo_cubesat_mpc_demo()
    println("RPO single-case run complete. Results: $(abspath(demo.csv_path))")
    println("Plots:")
    for path in demo.plot_paths
        println("  ", abspath(path))
    end
end
