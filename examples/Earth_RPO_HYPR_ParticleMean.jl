include(joinpath(@__DIR__, "Earth_RPO_CubeSat_MPC_PlannerComparison.jl"))
using PlotlyJS

"""Export a finite-cost companion figure while explicitly showing invalid particles."""
function write_rpo_hypr_finite_particle_mean_plot(plans, output_directory)
    batch = (plans_by_planner=Dict(:hypr => [
        (particle_mean_cost_history=plan.particle_finite_mean_cost_history,) for plan in plans
    ]),)
    plot = SM.GuidanceHooks.rpo_comparison_cost_iteration_plot(batch; planner=:hypr, history=:particle_mean)
    plot.data[3][:name] = "Mean finite particle cost"
    plot.data[3][:hovertemplate] = "iteration %{x}<br>mean finite particle cost %{y:.6g}<br>cases %{customdata}<extra></extra>"
    iters = collect(1:maximum(length(plan.particle_invalid_fraction_history) for plan in plans))
    invalid_pct = [100 * sum(plan.particle_invalid_fraction_history[
        min(iter, length(plan.particle_invalid_fraction_history))] for plan in plans) / length(plans)
        for iter in iters]
    PlotlyJS.addtraces!(plot, PlotlyJS.scatter(x=iters, y=invalid_pct,
        mode="lines", line=PlotlyJS.attr(color="black", width=2, dash="dash"),
        xaxis="x2", yaxis="y2", showlegend=false,
        hovertemplate="iteration %{x}<br>invalid particles %{y:.2f}%<extra></extra>"))
    PlotlyJS.relayout!(plot;
        height=440,
        yaxis_domain=[0.34, 1.0], yaxis_title_text="Finite particle cost",
        xaxis_title_text="", xaxis_showticklabels=false,
        xaxis2=PlotlyJS.attr(title="Iteration", anchor="y2", matches="x",
            ticks="outside", showline=true, mirror=true, linecolor="black", showgrid=false, zeroline=false),
        yaxis2=PlotlyJS.attr(title="Invalid (%)", domain=[0.0, 0.18], anchor="x2",
            ticks="outside", showline=true, mirror=true, linecolor="black",
            gridcolor="rgb(230,230,230)", rangemode="tozero", zeroline=false, nticks=3),
    )
    base = joinpath(output_directory, "rpo_planner_finite_particle_mean_cost_vs_iteration_hypr")
    PlotlyJS.savefig(plot, base * ".html")
    PlotlyJS.savefig(plot, base * ".pdf"; width=nothing, height=nothing)
    PlotlyJS.savefig(plot, base * ".png"; width=nothing, height=nothing, scale=6.25)
    open(base * "_caption.txt", "w") do io
        println(io, "HYPR finite particle cost versus iteration across $(length(plans)) cases. ",
            "Upper panel: each case averages only finite current particle costs; the solid line ",
            "averages those case means equally and shading spans their minimum and maximum. ",
            "Lower panel: the mean percentage of non-finite particle costs across cases. ",
            "Final values are carried forward after termination. Increases reflect swarm exploration. ",
            "This companion statistic differs from the literal all-particle mean when costs are infinite. ",
            "Insert at 7.16 inches wide.")
    end
    return plot
end

"""Rerun the saved comparison endpoints with complete swarm-cost logging."""
function run_rpo_hypr_particle_mean(;
    source_directory=joinpath(REPO_ROOT, "output", "rpo_planner_comparison_500_tol5mm"),
    output_directory=joinpath(source_directory, "hypr_particle_mean"),
)
    mkpath(output_directory)
    GH = SM.GuidanceHooks
    points = SpaceAGORA.load_rpo_station_cad_pointcloud(:gateway;
        n_points=10000, rng=MersenneTwister(740))
    geometry = RPOReferenceGeometry(
        RPOStationGeometry(points; keepout_radius_m=0.25, name="gateway_core");
        chaser=RPOCubeSatGeometry(dims_m=(0.1, 0.1, 0.3)))
    pso = rpo_740_mpc_final_pso_config(safe_distance_m=0.25,
        n_particles=100, n_iters=60, n_waypoints=5,
        sample_ds_m=0.5, refinement_sample_ds_m=0.5,
        clearance_feasibility_tol_m=0.005, adaptive_enable=true, refinement_enable=true)
    tracking = RPOLQMPCTrackingSettings(dt_s=pso.retime_dt_s, horizon=60,
        settle_time_s=20.0, final_position_tol_m=0.25,
        u_max_mps2=SVector(0.0125, 0.0125, 0.0125),
        q_pos=10.0, q_vel=1.0, r_accel=0.1, qf_pos=50.0, qf_vel=5.0)
    cfg = RPOPlannerComparisonConfig(planners=[:hypr], pso_config=pso,
        tracking=tracking, safe_distance_m=0.25, rng_seed=10000740)
    cases = collect(CSV.File(joinpath(source_directory, "rpo_planner_comparison_cases.csv")))
    plans = NamedTuple[]
    open(joinpath(output_directory, "run_settings.txt"), "w") do io
        println(io, "Source endpoints: ", abspath(source_directory))
        println(io, "Cases: ", length(cases), "; Julia threads: ", Threads.nthreads())
        println(io, "Current workspace implementation; complete costs with no objective cutoff.")
        println(io, "Mean of evaluated current particles before culling; Inf is retained.")
        println(io, "Planner runs only; no terminal tracking rerun.")
        show(io, MIME("text/plain"), cfg)
    end
    for (index, row) in enumerate(cases)
        case = RPOPlannerComparisonCase(
            start_rtn=SVector(row.start_x, row.start_y, row.start_z),
            goal_rtn=SVector(row.goal_x, row.goal_y, row.goal_z),
            label=lpad(row.case_id, 3, '0'))
        seed_path, _ = GH.rpo_pso_rrt_warmstart_path(case.start_rtn, case.goal_rtn,
            geometry, pso, 0.25, _rpo_50_of_50_planner_rng(:rrt_connect_seed, row.case_id, case))
        seed_path === nothing && (seed_path = hcat(case.start_rtn, case.goal_rtn))
        started = time()
        plan = GH.rpo_pso_plan_path(case.start_rtn, case.goal_rtn, geometry, pso;
            safe_distance_m=0.25,
            rng=_rpo_50_of_50_planner_rng(:hypr, row.case_id, case),
            initial_path=seed_path,
            objective_evaluator=_rpo_new_hypr_objective_evaluator(case, geometry, cfg),
            record_particle_mean=true)
        push!(plans, (cost_history=plan.cost_history,
            particle_mean_cost_history=plan.particle_mean_cost_history,
            particle_finite_mean_cost_history=plan.particle_finite_mean_cost_history,
            particle_invalid_fraction_history=plan.particle_invalid_fraction_history))
        CSV.write(joinpath(output_directory, "case_$(case.label).csv"), DataFrame(
            iteration=collect(eachindex(plan.cost_history)),
            best_cost=plan.cost_history,
            particle_mean_cost=plan.particle_mean_cost_history,
            finite_particle_mean_cost=plan.particle_finite_mean_cost_history,
            invalid_particle_fraction=plan.particle_invalid_fraction_history,
            particle_count=fill(plan.config.n_particles, length(plan.cost_history))))
        println("Case $(index)/$(length(cases)): $(length(plan.cost_history)) iterations, ",
            round(time() - started; digits=1), " s; nonfinite swarm means=",
            count(!isfinite, plan.particle_mean_cost_history))
        flush(stdout)
    end
    batch = (plans_by_planner=Dict(:hypr => plans),)
    plot = GH.rpo_comparison_cost_iteration_plot(batch; planner=:hypr, history=:particle_mean)
    base = joinpath(output_directory, "rpo_planner_particle_mean_cost_vs_iteration_hypr")
    PlotlyJS.savefig(plot, base * ".html")
    PlotlyJS.savefig(plot, base * ".pdf"; width=nothing, height=nothing)
    PlotlyJS.savefig(plot, base * ".png"; width=nothing, height=nothing, scale=6.25)
    open(base * "_caption.txt", "w") do io
        println(io, "HYPR swarm cost versus iteration across $(length(cases)) cases. ",
            "Each case records the arithmetic mean of all current particle costs before culling. ",
            "The line averages these case means equally; shading spans their minimum and maximum. ",
            "Final swarm means are carried forward after termination. ",
            "Non-finite means are retained and create gaps, rather than excluding particles or cases. ",
            "Increases can reflect particle exploration. Insert at 7.16 inches wide.")
    end
    println("Saved particle-mean figures: ", base)
    write_rpo_hypr_finite_particle_mean_plot(plans, output_directory)
    return batch
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_rpo_hypr_particle_mean()
end
