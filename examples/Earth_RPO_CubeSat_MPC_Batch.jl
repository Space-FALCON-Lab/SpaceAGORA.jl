include(joinpath(@__DIR__, "Earth_RPO_CubeSat_MPC.jl"))

using CSV
using DataFrames
using PlotlyJS
using Printf
using Random
using Statistics
using StaticArrays
using LinearAlgebra

const RPO_NAV = SM.NavigationHooks

function _env_int(name::AbstractString, default::Integer)
    token = strip(get(ENV, name, ""))
    isempty(token) && return Int(default)
    return parse(Int, token)
end

function _env_float(name::AbstractString, default::Real)
    token = strip(get(ENV, name, ""))
    isempty(token) && return Float64(default)
    return parse(Float64, token)
end

function _env_bool(name::AbstractString, default::Bool)
    token = lowercase(strip(get(ENV, name, "")))
    isempty(token) && return default
    token in ("1", "true", "yes", "on") && return true
    token in ("0", "false", "no", "off") && return false
    error("Environment variable $name must be boolean-like, got '$token'.")
end

function _rpo_batch_smoke_mode()
    return get(ENV, "SPACEAGORA_EXAMPLE_SMOKE", "0") == "1"
end

function _rpo_batch_case_label(case_id::Integer)
    return lpad(case_id, 3, '0')
end

function _progress_bar_string(frac::Real; width::Integer=30)
    f = clamp(Float64(frac), 0.0, 1.0)
    filled = Int(floor(f * width))
    empty = width - filled
    return "[" * repeat("=", filled) * repeat(" ", empty) * "]"
end

function _print_progress_line!(label::AbstractString, frac::Real, t_s::Real, tf_s::Real)
    percent = 100.0 * clamp(Float64(frac), 0.0, 1.0)
    @printf("\r%s %s %5.1f%%  t=%7.1f/%7.1f s", label, _progress_bar_string(frac), percent, Float64(t_s), Float64(tf_s))
    flush(stdout)
    return nothing
end

function _rpo_triangle_areas_and_cdf(triangles::AbstractMatrix)
    ntri = size(triangles, 2) ÷ 3
    areas = zeros(Float64, ntri)
    for tri_idx in 1:ntri
        base = 3 * (tri_idx - 1)
        v1 = SVector{3, Float64}(triangles[:, base + 1])
        v2 = SVector{3, Float64}(triangles[:, base + 2])
        v3 = SVector{3, Float64}(triangles[:, base + 3])
        areas[tri_idx] = 0.5 * norm(cross(v2 - v1, v3 - v1))
    end
    total = sum(areas)
    total > 0.0 || error("Station mesh has zero triangle area; cannot sample batch endpoints.")
    return areas, cumsum(areas) ./ total
end

function _rpo_sample_triangle_point(rng, v1, v2, v3)
    u = rand(rng)
    v = rand(rng)
    if u + v > 1.0
        u = 1.0 - u
        v = 1.0 - v
    end
    return v1 + u * (v2 - v1) + v * (v3 - v1)
end

function _rpo_directional_nearest_distances(p, points)
    min_d2 = fill(Inf, 6)
    for i in 1:size(points, 2)
        dx = points[1, i] - p[1]
        dy = points[2, i] - p[2]
        dz = points[3, i] - p[3]
        d2 = dx * dx + dy * dy + dz * dz
        if dx > 0.0 && d2 < min_d2[1]
            min_d2[1] = d2
        end
        if dx < 0.0 && d2 < min_d2[2]
            min_d2[2] = d2
        end
        if dy > 0.0 && d2 < min_d2[3]
            min_d2[3] = d2
        end
        if dy < 0.0 && d2 < min_d2[4]
            min_d2[4] = d2
        end
        if dz > 0.0 && d2 < min_d2[5]
            min_d2[5] = d2
        end
        if dz < 0.0 && d2 < min_d2[6]
            min_d2[6] = d2
        end
    end
    return sqrt.(min_d2)
end

function _rpo_is_surrounded_endpoint(p, points; max_distance=nothing)
    dists = _rpo_directional_nearest_distances(p, points)
    max_distance === nothing && return all(isfinite, dists)
    return all(d -> d <= max_distance, dists)
end

function _sample_rpo_near_surface_endpoint(
    rng,
    geometry,
    triangles,
    tri_cdf;
    min_clearance_m::Real,
    max_clearance_m::Real,
    max_tries::Integer=4000,
)
    buffer = geometry.station.keepout_radius_m + maximum(geometry.chaser.half_extents_body)
    mesh_center = SVector{3, Float64}(vec(mean(triangles; dims=2)))
    ntri = size(triangles, 2) ÷ 3
    for _ in 1:max_tries
        tri_idx = clamp(searchsortedfirst(tri_cdf, rand(rng)), 1, ntri)
        base = 3 * (tri_idx - 1)
        v1 = SVector{3, Float64}(triangles[:, base + 1])
        v2 = SVector{3, Float64}(triangles[:, base + 2])
        v3 = SVector{3, Float64}(triangles[:, base + 3])
        surface_point = _rpo_sample_triangle_point(rng, v1, v2, v3)
        normal = cross(v2 - v1, v3 - v1)
        normal_norm = norm(normal)
        normal_norm < 1.0e-9 && continue
        normal = normal / normal_norm
        tri_center = (v1 + v2 + v3) / 3.0
        dot(normal, tri_center - mesh_center) < 0.0 && (normal = -normal)

        clearance = Float64(min_clearance_m) + rand(rng) * max(0.0, Float64(max_clearance_m - min_clearance_m))
        candidate = surface_point + normal * (buffer + clearance)
        # The old 740 runner samples a clearance in the requested shell and then
        # performs a final obstacle-clearance sanity check. Here the CAD normal
        # offset defines the shell distance; the point-cloud check rejects points
        # that are too close to nearby structure without making generation hinge
        # on an expensive all-triangle distance query.
        final_clearance = RPO_NAV.rpo_clearance_to_station(candidate, geometry).clearance
        if final_clearance >= min_clearance_m
            return candidate
        end
    end
    error("Failed to sample a near-surface RPO endpoint after $(max_tries) tries.")
end

function generate_rpo_seeded_batch_cases(;
    n_cases::Integer=200,
    seed::Integer=740,
    geometry_seed::Integer=740,
    n_station_points::Integer=10000,
    endpoint_min_clearance_m::Real=0.26,
    endpoint_max_clearance_m::Real=1.0,
    min_separation_m::Real=1.5,
    surrounded_max_distance_m=2.0,
    max_sampling_tries::Integer=4000,
)
    station_points = SpaceAGORA.load_rpo_station_cad_pointcloud(
        :gateway;
        n_points=n_station_points,
        rng=MersenneTwister(geometry_seed),
    )
    geometry = RPOReferenceGeometry(
        RPOStationGeometry(station_points; keepout_radius_m=0.25, name="gateway_core");
        chaser=RPOCubeSatGeometry(dims_m=(0.1, 0.1, 0.3)),
    )
    station_triangles = SpaceAGORA.load_rpo_station_cad_triangles(:gateway)
    _, tri_cdf = _rpo_triangle_areas_and_cdf(station_triangles)
    endpoint_max_clearance_m = max(Float64(endpoint_max_clearance_m), Float64(endpoint_min_clearance_m))
    rng = MersenneTwister(seed)
    cases = NamedTuple[]
    for case_id in 1:n_cases
        start = nothing
        for _ in 1:max_sampling_tries
            candidate = _sample_rpo_near_surface_endpoint(
                rng,
                geometry,
                station_triangles,
                tri_cdf;
                min_clearance_m=endpoint_min_clearance_m,
                max_clearance_m=endpoint_max_clearance_m,
                max_tries=max_sampling_tries,
            )
            if !_rpo_is_surrounded_endpoint(
                candidate,
                station_points;
                max_distance=surrounded_max_distance_m,
            )
                start = candidate
                break
            end
        end
        start === nothing && error("Failed to sample a start that is not surrounded after $(max_sampling_tries) tries.")

        goal = nothing
        for _ in 1:max_sampling_tries
            candidate = _sample_rpo_near_surface_endpoint(
                rng,
                geometry,
                station_triangles,
                tri_cdf;
                min_clearance_m=endpoint_min_clearance_m,
                max_clearance_m=endpoint_max_clearance_m,
                max_tries=max_sampling_tries,
            )
            if norm(candidate - start) >= min_separation_m &&
                !_rpo_is_surrounded_endpoint(
                    candidate,
                    station_points;
                    max_distance=surrounded_max_distance_m,
                )
                goal = candidate
                break
            end
        end
        goal === nothing && error("Failed to sample a goal with minimum separation and not surrounded after $(max_sampling_tries) tries.")
        push!(cases, (case_id=case_id, seed=seed + case_id, start_rtn=start, goal_rtn=goal))
    end
    return cases
end

function _rpo_control_effector_from_integrator(integrator)
    for effector in integrator.p.args.control_model.control_effectors
        effector isa RPOMPCControlModel && return effector
    end
    return nothing
end

function _rpo_demo_chaser_propellant_mass_kg(demo)::Float64
    spacecraft = demo.args.dynamics_model.spacecraft
    isempty(spacecraft) && return NaN
    prop_mass = Float64(spacecraft[1].prop_mass)
    return isfinite(prop_mass) && prop_mass > 0.0 ? prop_mass : NaN
end

function _rpo_fuel_used_pct(fuel_used_kg::Real, propellant_mass_kg::Real)::Float64
    denom = Float64(propellant_mass_kg)
    return isfinite(denom) && denom > 0.0 ? 100.0 * Float64(fuel_used_kg) / denom : NaN
end

function _rpo_batch_progress_save_field(args, case_label::AbstractString)
    tf_s = Float64(args.mission_configuration.mission_time)
    last_bucket = Ref(-1)
    label = "  sim case $case_label"
    return SM.SaveField(
        :rpo_progress_fraction,
        (u, t, integrator) -> begin
            frac = tf_s <= 0.0 ? 1.0 : clamp(Float64(t) / tf_s, 0.0, 1.0)
            bucket = Int(floor(100.0 * frac))
            if bucket > last_bucket[] || frac >= 1.0
                last_bucket[] = bucket
                _print_progress_line!(label, frac, t, tf_s)
                frac >= 1.0 && println("")
            end
            return frac
        end;
        column_prefix="rpo_progress_fraction",
    )
end

function _rpo_batch_save_fields(args; case_label::AbstractString="")
    fields = SM.default_save_fields(args)
    if !isempty(case_label)
        push!(fields, _rpo_batch_progress_save_field(args, case_label))
    end
    push!(
        fields,
        SM.SaveField(
            :rpo_thruster_forces,
            (u, t, integrator) -> begin
                effector = _rpo_control_effector_from_integrator(integrator)
                effector === nothing ? SVector{6, Float64}(zeros(6)) : effector.held.thruster_forces_n
            end;
            column_prefix="rpo_thruster_force",
        ),
    )
    push!(
        fields,
        SM.SaveField(
            :rpo_control_force_ii,
            (u, t, integrator) -> begin
                effector = _rpo_control_effector_from_integrator(integrator)
                effector === nothing ? SVector{3, Float64}(0.0, 0.0, 0.0) : effector.held.force_ii
            end;
            column_prefix="rpo_control_force_ii",
        ),
    )
    push!(
        fields,
        SM.SaveField(
            :rpo_control_torque_body,
            (u, t, integrator) -> begin
                effector = _rpo_control_effector_from_integrator(integrator)
                effector === nothing ? SVector{3, Float64}(0.0, 0.0, 0.0) : effector.held.torque_body
            end;
            column_prefix="rpo_control_torque_body",
        ),
    )
    return fields
end

function _trapz(t, y)
    length(t) < 2 && return 0.0
    total = 0.0
    for i in 2:length(t)
        total += 0.5 * (y[i] + y[i - 1]) * (t[i] - t[i - 1])
    end
    return total
end

function _column_or_zeros(df, name::Symbol)
    return name in propertynames(df) ? Float64.(df[!, name]) : zeros(Float64, nrow(df))
end

function _rpo_batch_case_metrics(case, demo, csv_path::AbstractString, planner_runtime_s::Real, simulation_runtime_s::Real)
    df, actual_rtn, ref_rtn, tracking_error = _rpo_postprocess(csv_path, demo)
    t = Float64.(df.time)
    mass = Float64.(df.sc1_mass)
    force_norm = sqrt.(
        _column_or_zeros(df, :rpo_control_force_ii_1).^2 .+
        _column_or_zeros(df, :rpo_control_force_ii_2).^2 .+
        _column_or_zeros(df, :rpo_control_force_ii_3).^2
    )
    torque_norm = sqrt.(
        _column_or_zeros(df, :rpo_control_torque_body_1).^2 .+
        _column_or_zeros(df, :rpo_control_torque_body_2).^2 .+
        _column_or_zeros(df, :rpo_control_torque_body_3).^2
    )
    thruster_cols = [_column_or_zeros(df, Symbol("rpo_thruster_force_$i")) for i in 1:6]
    thruster_max = maximum(Float64.(demo.control.thrusters.max_thrust_n))
    thruster_any_sat = [
        any(col[k] >= 0.995 * thruster_max for col in thruster_cols)
        for k in 1:nrow(df)
    ]
    rw_limit = demo.control.max_rw_torque_nm
    rw_sat = isfinite(rw_limit) && rw_limit > 0.0 ?
        count(x -> x >= 0.995 * rw_limit, torque_norm) / max(length(torque_norm), 1) :
        0.0
    clearances = [
        RPO_NAV.rpo_clearance_to_station(SVector{3, Float64}(actual_rtn[:, k]), demo.geometry).clearance
        for k in axes(actual_rtn, 2)
    ]

    final_goal_error = norm(SVector{3, Float64}(actual_rtn[:, end]) - SVector{3, Float64}(case.goal_rtn))
    fuel_used = mass[1] - mass[end]
    propellant_mass_kg = _rpo_demo_chaser_propellant_mass_kg(demo)
    success = final_goal_error <= 0.75 && minimum(clearances) >= -0.05 && fuel_used >= -1.0e-8

    return (
        case_id=case.case_id,
        seed=case.seed,
        start_x=case.start_rtn[1],
        start_y=case.start_rtn[2],
        start_z=case.start_rtn[3],
        goal_x=case.goal_rtn[1],
        goal_y=case.goal_rtn[2],
        goal_z=case.goal_rtn[3],
        success=success ? 1 : 0,
        planner_runtime_s=Float64(planner_runtime_s),
        simulation_runtime_s=Float64(simulation_runtime_s),
        total_runtime_s=Float64(planner_runtime_s + simulation_runtime_s),
        planned_duration_s=demo.initial_plan.t_ref_s[end],
        fuel_used_kg=fuel_used,
        fuel_used_pct=_rpo_fuel_used_pct(fuel_used, propellant_mass_kg),
        propellant_mass_kg=propellant_mass_kg,
        force_impulse_ns=_trapz(t, force_norm),
        torque_impulse_nms=_trapz(t, torque_norm),
        thrust_saturation_fraction=count(identity, thruster_any_sat) / max(nrow(df), 1),
        reaction_wheel_saturation_fraction=rw_sat,
        tracking_final_m=tracking_error[end],
        tracking_mean_m=mean(tracking_error),
        tracking_max_m=maximum(tracking_error),
        goal_error_final_m=final_goal_error,
        min_clearance_m=minimum(clearances),
        pso_final_cost=demo.plan_result.cost,
        pso_initial_cost=isempty(demo.plan_result.cost_history) ? demo.plan_result.cost : first(demo.plan_result.cost_history),
        pso_iterations=length(demo.plan_result.cost_history),
        pso_refinement_improved=demo.plan_result.refinement_improved ? 1 : 0,
        csv_path=abspath(csv_path),
    )
end

function _write_rpo_batch_summary(path::AbstractString, rows)
    df = DataFrame(rows)
    CSV.write(path, df)
    return df
end

function _rpo_seeded_cases_dataframe(cases)
    return DataFrame([
        (
            case_id=case.case_id,
            seed=case.seed,
            start_x=case.start_rtn[1],
            start_y=case.start_rtn[2],
            start_z=case.start_rtn[3],
            goal_x=case.goal_rtn[1],
            goal_y=case.goal_rtn[2],
            goal_z=case.goal_rtn[3],
        )
        for case in cases
    ])
end

function _save_plot(plot_obj, output_dir::AbstractString, filename::AbstractString)
    mkpath(output_dir)
    path = joinpath(output_dir, filename)
    PlotlyJS.savefig(plot_obj, path)
    return path
end

function _tracking_family_plot(case_traces)
    traces = PlotlyJS.GenericTrace[]
    for trace in case_traces
        push!(traces, scatter(
            x=trace.time_s,
            y=trace.tracking_error_m,
            mode="lines",
            name="case $(_rpo_batch_case_label(trace.case_id))",
        ))
    end
    return Plot(traces, Layout(title="RPO Batch Tracking Error Family", xaxis_title="time (s)", yaxis_title="tracking error (m)"))
end

function _family_metric_specs()
    return [
        (:success_pct, "Success", "success (%)"),
        (:total_runtime_s, "Runtime", "s"),
        (:fuel_used_pct, "Fuel", "%"),
        (:thrust_saturation_fraction, "Thruster Sat.", "fraction"),
        (:reaction_wheel_saturation_fraction, "RW Sat.", "fraction"),
        (:tracking_max_m, "Tracking Max", "m"),
        (:tracking_mean_m, "Tracking Mean", "m"),
        (:goal_error_final_m, "Goal Error", "m"),
        (:min_clearance_m, "Min Clearance", "m"),
    ]
end

function _metric_values(df::DataFrame, metric::Symbol)
    metric == :success_pct && return 100.0 .* Float64.(df.success)
    return Float64.(df[!, metric])
end

@inline _metric_rows(df::DataFrame, metric::Symbol) = metric == :success_pct ? df : df[df.success .== 1, :]

function _axis_name(prefix::AbstractString, idx::Integer)
    idx == 1 && return Symbol(prefix)
    return Symbol(prefix, idx)
end

function _trace_axis_name(prefix::AbstractString, idx::Integer)
    idx == 1 && return prefix
    return string(prefix, idx)
end

function _rpo_family_metric_summary_plot(df::DataFrame; rows::Integer=3, cols::Integer=3)
    specs = _family_metric_specs()
    traces = PlotlyJS.GenericTrace[]
    layout_kwargs = Dict{Symbol, Any}(
        :title => "RPO Batch Family Metrics",
        :height => 920,
        :width => 1200,
        :showlegend => false,
        :margin => attr(l=70, r=30, t=80, b=50),
        :plot_bgcolor => "rgb(250,250,250)",
    )
    xgap = 0.055
    ygap = 0.095
    xspan = (1.0 - xgap * (cols - 1)) / cols
    yspan = (1.0 - ygap * (rows - 1)) / rows

    for (idx, (metric, label, units)) in enumerate(specs)
        row = div(idx - 1, cols) + 1
        col = mod(idx - 1, cols) + 1
        x0 = (col - 1) * (xspan + xgap)
        x1 = x0 + xspan
        y1 = 1.0 - (row - 1) * (yspan + ygap)
        y0 = y1 - yspan
        xaxis = _trace_axis_name("x", idx)
        yaxis = _trace_axis_name("y", idx)
        metric_df = _metric_rows(df, metric)
        values = _metric_values(metric_df, metric)
        xvals = fill(label, length(values))
        case_labels = ["case $(_rpo_batch_case_label(id))" for id in metric_df.case_id]

        push!(traces, box(
            x=xvals,
            y=values,
            name=label,
            boxpoints=false,
            boxmean=true,
            fillcolor="rgba(64,120,165,0.28)",
            line=attr(color="rgb(42,86,125)", width=2),
            marker=attr(color="rgb(42,86,125)"),
            xaxis=xaxis,
            yaxis=yaxis,
            hovertemplate="$(label)<br>%{y:.5g} $(units)<extra></extra>",
        ))
        push!(traces, scatter(
            x=xvals,
            y=values,
            mode="markers",
            name="$(label) cases",
            marker=attr(size=8, color="rgb(205,92,55)", line=attr(color="white", width=1)),
            text=case_labels,
            hovertemplate="%{text}<br>$(label): %{y:.5g} $(units)<extra></extra>",
            xaxis=xaxis,
            yaxis=yaxis,
        ))

        layout_kwargs[_axis_name("xaxis", idx)] = attr(
            domain=[x0, x1],
            anchor=yaxis,
            showticklabels=false,
            ticks="",
            title=label,
        )
        layout_kwargs[_axis_name("yaxis", idx)] = attr(
            domain=[y0, y1],
            anchor=xaxis,
            title=units,
            zeroline=true,
            gridcolor="rgb(225,225,225)",
        )
    end
    return Plot(traces, Layout(; layout_kwargs...))
end

function _path_family_plot(case_traces)
    traces = PlotlyJS.GenericTrace[_station_mesh_trace()]
    for trace in case_traces
        label = _rpo_batch_case_label(trace.case_id)
        push!(traces, scatter3d(
            x=trace.actual_rtn[1, :],
            y=trace.actual_rtn[2, :],
            z=trace.actual_rtn[3, :],
            mode="lines",
            name="case $(label) actual",
        ))
        push!(traces, scatter3d(
            x=trace.ref_rtn[1, :],
            y=trace.ref_rtn[2, :],
            z=trace.ref_rtn[3, :],
            mode="lines",
            line=attr(dash="dot"),
            name="case $(label) ref",
        ))
        push!(traces, scatter3d(
            x=[trace.start_rtn[1]],
            y=[trace.start_rtn[2]],
            z=[trace.start_rtn[3]],
            mode="markers",
            marker=attr(size=6, color="rgb(45,150,90)", symbol="circle"),
            name="case $(label) start",
        ))
        push!(traces, scatter3d(
            x=[trace.goal_rtn[1]],
            y=[trace.goal_rtn[2]],
            z=[trace.goal_rtn[3]],
            mode="markers",
            marker=attr(size=7, color="rgb(190,60,55)", symbol="diamond"),
            name="case $(label) goal",
        ))
    end
    return Plot(
        traces,
        Layout(
            title="RPO Batch Path Family",
            scene=attr(aspectmode="data", xaxis_title="radial (m)", yaxis_title="along-track (m)", zaxis_title="cross-track (m)"),
        ),
    )
end

function _save_rpo_batch_family_plots(output_dir::AbstractString, df::DataFrame, case_traces)
    plot_paths = String[]
    push!(plot_paths, _save_plot(_rpo_family_metric_summary_plot(df), output_dir, "family_metric_summary.html"))
    push!(plot_paths, _save_plot(_path_family_plot(case_traces), output_dir, "family_paths.html"))
    return plot_paths
end

function run_rpo_cubesat_mpc_batch(;
    n_cases::Integer=_rpo_batch_smoke_mode() ? 1 : _env_int("SPACEAGORA_RPO_BATCH_N", 200),
    seed::Integer=_env_int("SPACEAGORA_RPO_BATCH_SEED", 740),
    mission_time::Real=_env_float("SPACEAGORA_EXAMPLE_SMOKE_MISSION_TIME", 180.0),
    results_directory::AbstractString=joinpath(REPO_ROOT, "output", "rpo_batch_cases"),
    pso_n_particles::Integer=_rpo_batch_smoke_mode() ? 10 : _env_int("SPACEAGORA_RPO_BATCH_PSO_PARTICLES", 120),
    pso_n_iters::Integer=_rpo_batch_smoke_mode() ? 4 : _env_int("SPACEAGORA_RPO_BATCH_PSO_ITERS", 35),
    pso_config=nothing,
    pso_configurator=nothing,
    n_station_points::Integer=_rpo_batch_smoke_mode() ? 800 : _env_int("SPACEAGORA_RPO_BATCH_STATION_POINTS", 10000),
    save_case_plots::Bool=!_rpo_batch_smoke_mode(),
)
    if pso_config !== nothing && pso_configurator !== nothing
        throw(ArgumentError("Pass either pso_config or pso_configurator, not both."))
    end
    mkpath(results_directory)
    cases = generate_rpo_seeded_batch_cases(
        n_cases=n_cases,
        seed=seed,
        geometry_seed=seed,
        n_station_points=n_station_points,
    )
    rows = NamedTuple[]
    case_traces = NamedTuple[]
    for case in cases
        label = _rpo_batch_case_label(case.case_id)
        case_dir = joinpath(results_directory, "case_$label")
        println("RPO batch case $label seed=$(case.seed)")
        print("  planning case $label ... ")
        flush(stdout)
        planner_runtime = @elapsed demo = build_rpo_cubesat_mpc_demo(
            mission_time=mission_time,
            results_directory=case_dir,
            seed=case.seed,
            station_geometry_seed=seed,
            start_rtn=case.start_rtn,
            goal_rtn=case.goal_rtn,
            pso_n_particles=pso_n_particles,
            pso_n_iters=pso_n_iters,
            pso_config=pso_config,
            pso_configurator=pso_configurator,
            n_station_points=n_station_points,
            data_rate_s=_rpo_batch_smoke_mode() ? 10.0 : 2.0,
            verbose=false,
        )
        @printf("done in %.2f s\n", planner_runtime)
        simulation_runtime = @elapsed run_simulation(demo.args; save_fields=_rpo_batch_save_fields(demo.args; case_label=label))
        csv_path = joinpath(case_dir, "simulation_results.csv")
        isfile(csv_path) || error("Missing simulation results for RPO batch case $label at $csv_path")
        row = _rpo_batch_case_metrics(case, demo, csv_path, planner_runtime, simulation_runtime)
        push!(rows, row)
        df_case, actual_rtn, ref_rtn, tracking_error = _rpo_postprocess(csv_path, demo)
        push!(case_traces, (
            case_id=case.case_id,
            time_s=Float64.(df_case.time),
            actual_rtn=actual_rtn,
            ref_rtn=ref_rtn,
            tracking_error_m=tracking_error,
            start_rtn=case.start_rtn,
            goal_rtn=case.goal_rtn,
        ))
        if save_case_plots
            _save_rpo_single_case_plots(csv_path, demo)
        end
        @printf("  success=%d fuel=%.6f kg max_tracking=%.3f m runtime=%.2f s\n", row.success, row.fuel_used_kg, row.tracking_max_m, row.total_runtime_s)
    end

    summary_path = joinpath(results_directory, "rpo_batch_summary.csv")
    summary_df = _write_rpo_batch_summary(summary_path, rows)
    plot_paths = _save_rpo_batch_family_plots(results_directory, summary_df, case_traces)
    println("RPO batch complete. Summary: $(abspath(summary_path))")
    println("Family plots:")
    for path in plot_paths
        println("  ", abspath(path))
    end
    return (
        cases=cases,
        summary=summary_df,
        summary_path=summary_path,
        plot_paths=plot_paths,
        output_dir=results_directory,
    )
end

function run_rpo_cubesat_mpc_planner_comparison_batch(;
    n_cases::Integer=_rpo_batch_smoke_mode() ? 1 : _env_int("SPACEAGORA_RPO_COMPARISON_N", 50),
    seed::Integer=_env_int("SPACEAGORA_RPO_COMPARISON_SEED", _env_int("SPACEAGORA_RPO_BATCH_SEED", 740)),
    results_directory::AbstractString=joinpath(REPO_ROOT, "output", "rpo_planner_comparison_cases"),
    pso_n_particles::Integer=_rpo_batch_smoke_mode() ? 8 : _env_int("SPACEAGORA_RPO_COMPARISON_PSO_PARTICLES", 100),
    pso_n_iters::Integer=_rpo_batch_smoke_mode() ? 2 : _env_int("SPACEAGORA_RPO_COMPARISON_PSO_ITERS", 60),
    pso_n_waypoints::Integer=_rpo_batch_smoke_mode() ? 1 : _env_int("SPACEAGORA_RPO_COMPARISON_PSO_WAYPOINTS", 5),
    rrt_connect_iters::Integer=_rpo_batch_smoke_mode() ? 25 : _env_int("SPACEAGORA_RPO_COMPARISON_RRT_CONNECT_ITERS", 1000),
    rrt_connect_step_size_m::Real=_env_float("SPACEAGORA_RPO_COMPARISON_RRT_CONNECT_STEP_SIZE", 0.75),
    rrt_star_iters::Integer=_rpo_batch_smoke_mode() ? 25 : _env_int("SPACEAGORA_RPO_COMPARISON_RRT_STAR_ITERS", rrt_connect_iters),
    rrt_star_step_size_m::Real=_env_float("SPACEAGORA_RPO_COMPARISON_RRT_STAR_STEP_SIZE", rrt_connect_step_size_m),
    rrt_star_neighbor_radius_m::Real=_env_float("SPACEAGORA_RPO_COMPARISON_RRT_STAR_NEIGHBOR_RADIUS", 2.0),
    chomp_iters::Integer=_rpo_batch_smoke_mode() ? 2 : _env_int("SPACEAGORA_RPO_COMPARISON_CHOMP_ITERS", pso_n_iters),
    stomp_iters::Integer=_rpo_batch_smoke_mode() ? 2 : _env_int("SPACEAGORA_RPO_COMPARISON_STOMP_ITERS", pso_n_iters),
    stomp_rollouts::Integer=_rpo_batch_smoke_mode() ? 3 : _env_int("SPACEAGORA_RPO_COMPARISON_STOMP_ROLLOUTS", 20),
    n_station_points::Integer=_rpo_batch_smoke_mode() ? 800 : _env_int("SPACEAGORA_RPO_COMPARISON_STATION_POINTS", 10000),
    safe_distance_m::Real=_env_float("SPACEAGORA_RPO_COMPARISON_SAFE_DISTANCE", 0.1),
    obstacle_weight_scale::Real=_env_float("SPACEAGORA_RPO_COMPARISON_OBS_WEIGHT_SCALE", 10.0),
    obstacle_margin_m::Real=_env_float("SPACEAGORA_RPO_COMPARISON_OBS_MARGIN", 0.5),
    match_hypr_runtime::Bool=_env_bool("SPACEAGORA_RPO_COMPARISON_MATCH_HYPR_RUNTIME", true),
    runtime_limit_s::Real=_env_float("SPACEAGORA_RPO_COMPARISON_RUNTIME_LIMIT_S", 30.0),
    force_full_iters::Bool=_env_bool("SPACEAGORA_RPO_COMPARISON_FORCE_FULL_ITERS", true),
    mpc_horizon::Integer=_rpo_batch_smoke_mode() ? 4 : _env_int("SPACEAGORA_RPO_COMPARISON_MPC_HORIZON", 80),
    mpc_u_max_mps2::Real=_rpo_batch_smoke_mode() ? 0.05 : _env_float("SPACEAGORA_RPO_COMPARISON_MPC_U_MAX", 0.01875),
    mpc_q_pos::Real=_rpo_batch_smoke_mode() ? 10.0 : _env_float("SPACEAGORA_RPO_COMPARISON_MPC_Q_POS", 30.0),
    mpc_q_vel::Real=_rpo_batch_smoke_mode() ? 1.0 : _env_float("SPACEAGORA_RPO_COMPARISON_MPC_Q_VEL", 4.0),
    mpc_r_accel::Real=_rpo_batch_smoke_mode() ? 0.1 : _env_float("SPACEAGORA_RPO_COMPARISON_MPC_R_ACCEL", 0.025),
    mpc_qf_pos::Real=_rpo_batch_smoke_mode() ? 50.0 : _env_float("SPACEAGORA_RPO_COMPARISON_MPC_QF_POS", 300.0),
    mpc_qf_vel::Real=_rpo_batch_smoke_mode() ? 5.0 : _env_float("SPACEAGORA_RPO_COMPARISON_MPC_QF_VEL", 25.0),
    mpc_settle_time_s::Real=_rpo_batch_smoke_mode() ? 2.0 : _env_float("SPACEAGORA_RPO_COMPARISON_MPC_SETTLE_TIME", 30.0),
    mpc_final_position_tol_m::Real=_rpo_batch_smoke_mode() ? 0.75 : _env_float("SPACEAGORA_RPO_COMPARISON_MPC_FINAL_TOL", 0.25),
    write_plotly_outputs::Bool=_env_bool("SPACEAGORA_RPO_COMPARISON_PLOTS", true),
    show_progress::Bool=_env_bool("SPACEAGORA_RPO_COMPARISON_PROGRESS", true),
)
    mkpath(results_directory)
    println("RPO planner comparison batch: HYPR vs PSO (unrefined) vs RRT-Connect vs RRT-Connect + Bezier vs RRT* vs CHOMP vs STOMP")
    println("  cases=$(n_cases), seed=$(seed), station_points=$(n_station_points)")
    println("  output=$(abspath(results_directory))")

    station_points = SpaceAGORA.load_rpo_station_cad_pointcloud(
        :gateway;
        n_points=n_station_points,
        rng=MersenneTwister(seed),
    )
    geometry = RPOReferenceGeometry(
        RPOStationGeometry(station_points; keepout_radius_m=0.25, name="gateway_core");
        chaser=RPOCubeSatGeometry(dims_m=(0.1, 0.1, 0.3)),
    )
    generated_cases = generate_rpo_seeded_batch_cases(
        n_cases=n_cases,
        seed=seed,
        geometry_seed=seed,
        n_station_points=n_station_points,
    )
    comparison_cases = [
        RPOPlannerComparisonCase(
            start_rtn=SVector{3, Float64}(case.start_rtn),
            goal_rtn=SVector{3, Float64}(case.goal_rtn),
            label=_rpo_batch_case_label(case.case_id),
        )
        for case in generated_cases
    ]

    pso_cfg = rpo_740_mpc_final_pso_config(
        safe_distance_m=safe_distance_m;
        n_particles=Int(pso_n_particles),
        n_iters=Int(pso_n_iters),
        n_waypoints=Int(pso_n_waypoints),
        adaptive_enable=!_rpo_batch_smoke_mode(),
        refinement_enable=true,
    )
    if _rpo_batch_smoke_mode()
        pso_cfg = rpo_pso_config(
            pso_cfg;
            cull_enable=false,
            schedule_enable=false,
            reexplore_enable=false,
            refinement_enable=false,
            sample_ds_m=0.5,
            retime_dt_s=0.5,
            retime_a_max_mps2=0.05,
            retime_max_steps=5_000,
        )
    end
    mpc_u_max = Float64(mpc_u_max_mps2)
    tracking = RPOLQMPCTrackingSettings(
        dt_s=pso_cfg.retime_dt_s,
        horizon=Int(mpc_horizon),
        settle_time_s=Float64(mpc_settle_time_s),
        final_position_tol_m=Float64(mpc_final_position_tol_m),
        u_max_mps2=SVector{3, Float64}(mpc_u_max, mpc_u_max, mpc_u_max),
        q_pos=Float64(mpc_q_pos),
        q_vel=Float64(mpc_q_vel),
        r_accel=Float64(mpc_r_accel),
        qf_pos=Float64(mpc_qf_pos),
        qf_vel=Float64(mpc_qf_vel),
    )
    cfg = RPOPlannerComparisonConfig(
        pso_config=pso_cfg,
        rrt_connect=RPORRTConnectSettings(
            n_iters=Int(rrt_connect_iters),
            step_size_m=Float64(rrt_connect_step_size_m),
            collision_sample_ds_m=max(pso_cfg.sample_ds_m, 0.10),
            shortcut_iters=_rpo_batch_smoke_mode() ? 4 : 80,
        ),
        rrt_star=RPORRTStarSettings(
            n_iters=Int(rrt_star_iters),
            step_size_m=Float64(rrt_star_step_size_m),
            collision_sample_ds_m=max(pso_cfg.sample_ds_m, 0.10),
            neighbor_radius_m=Float64(rrt_star_neighbor_radius_m),
            shortcut_iters=_rpo_batch_smoke_mode() ? 4 : 80,
        ),
        chomp=RPOCHOMPSettings(n_iters=Int(chomp_iters)),
        stomp=RPOSTOMPSettings(n_iters=Int(stomp_iters), n_rollouts=Int(stomp_rollouts)),
        optimizer=RPOTrajectoryOptimizerSettings(
            w_obs_scale=Float64(obstacle_weight_scale),
            obstacle_margin_m=Float64(obstacle_margin_m),
            match_hypr_runtime=match_hypr_runtime,
            runtime_limit_s=Float64(runtime_limit_s),
            force_full_iters=force_full_iters,
        ),
        tracking=tracking,
        safe_distance_m=Float64(safe_distance_m),
        output_dir=String(results_directory),
        write_plotly_outputs=write_plotly_outputs,
        rng_seed=Int(seed),
        show_progress=show_progress,
    )

    runtime = @elapsed batch = rpo_run_planner_comparison_batch(comparison_cases, geometry, cfg)
    batch = merge(batch, (station_triangles=SpaceAGORA.load_rpo_station_cad_triangles(:gateway),))
    outputs = rpo_write_planner_comparison_outputs(batch)
    rows = rpo_flatten_planner_results(batch)
    summary_df = DataFrame(rows)
    summary_path = outputs.csv
    cases_path = joinpath(results_directory, "rpo_planner_comparison_cases.csv")
    CSV.write(cases_path, _rpo_seeded_cases_dataframe(generated_cases))

    println("RPO planner comparison complete in $(round(runtime, digits=2)) s.")
    println("Summary CSV: ", abspath(summary_path))
    println("Cases CSV:   ", abspath(cases_path))
    if outputs.plotly_outputs
        println("Plots:")
        println("  ", abspath(outputs.metrics_plot))
        println("  ", abspath(outputs.path_plot))
        for planner in batch.planner_types
            path = outputs.method_path_plots[planner]
            println("  ", abspath(path))
        end
    else
        println("Plotly outputs: disabled")
    end
    planner_report_rows = NamedTuple[]
    for planner in cfg.planners
        planner_rows = filter(row -> row.planner == normalize_rpo_comparison_planner_type(planner), rows)
        isempty(planner_rows) && continue
        successful_rows = filter(row -> row.success, planner_rows)
        success_pct = 100.0 * mean([row.success ? 1.0 : 0.0 for row in planner_rows])
        mean_fuel_pct = isempty(successful_rows) ? NaN : mean([row.fuel_used_pct for row in successful_rows])
        fuel_values_pct = [Float64(row.fuel_used_pct) for row in successful_rows]
        std_fuel_pct = length(fuel_values_pct) > 1 ? std(fuel_values_pct) : (isempty(fuel_values_pct) ? NaN : 0.0)
        mean_runtime = mean([row.planner_compute_time for row in planner_rows])
        push!(planner_report_rows, (
            planner=planner,
            success_pct=success_pct,
            mean_fuel_pct=mean_fuel_pct,
            std_fuel_pct=std_fuel_pct,
            mean_runtime=mean_runtime,
        ))
    end
    hypr_success_fuel_by_case = Dict{Int, Float64}()
    for row in rows
        row.planner == :hypr && row.success || continue
        hypr_success_fuel_by_case[Int(row.case_id)] = Float64(row.fuel_used_pct)
    end
    for row in planner_report_rows
        planner_rows = filter(result -> result.planner == normalize_rpo_comparison_planner_type(row.planner), rows)
        paired_fuel_pct = Float64[]
        paired_hypr_fuel_pct = Float64[]
        for result in planner_rows
            result.success || continue
            hypr_fuel_pct = get(hypr_success_fuel_by_case, Int(result.case_id), NaN)
            isfinite(hypr_fuel_pct) || continue
            push!(paired_fuel_pct, Float64(result.fuel_used_pct))
            push!(paired_hypr_fuel_pct, hypr_fuel_pct)
        end
        paired_mean_fuel_pct = isempty(paired_fuel_pct) ? NaN : mean(paired_fuel_pct)
        paired_hypr_mean_fuel_pct = isempty(paired_hypr_fuel_pct) ? NaN : mean(paired_hypr_fuel_pct)
        fuel_diff_pct = isfinite(paired_mean_fuel_pct) && isfinite(paired_hypr_mean_fuel_pct) && paired_hypr_mean_fuel_pct > 0.0 ?
            100.0 * (paired_mean_fuel_pct - paired_hypr_mean_fuel_pct) / paired_hypr_mean_fuel_pct :
            NaN
        @printf("  %-21s success=%5.1f%% mean_fuel=%7.3f%% std_fuel=%6.3f%% fuel_vs_hypr=%+6.1f%% mean_planner=%.3f s\n",
            rpo_comparison_planner_label(row.planner), row.success_pct, row.mean_fuel_pct, row.std_fuel_pct, fuel_diff_pct, row.mean_runtime)
    end

    return (
        cases=generated_cases,
        batch=batch,
        summary=summary_df,
        summary_path=summary_path,
        cases_path=cases_path,
        outputs=outputs,
        output_dir=results_directory,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    mode = isempty(ARGS) ? get(ENV, "SPACEAGORA_RPO_BATCH_MODE", "batch") : lowercase(strip(ARGS[1]))
    if mode in ("batch", "hypr")
        run_rpo_cubesat_mpc_batch()
    elseif mode in ("comparison", "planner_comparison", "planners")
        run_rpo_cubesat_mpc_planner_comparison_batch()
    else
        error("Unknown RPO batch mode '$mode'. Use batch or comparison.")
    end
end
