function _trajectory_delta_metrics(
    sol_ref,
    sol_cmp,
    n_sats::Int,
    orientation::Bool;
    n_samples::Int
)
    t_ref_start = Float64(first(sol_ref.t))
    t_ref_end = Float64(last(sol_ref.t))
    t_cmp_start = Float64(first(sol_cmp.t))
    t_cmp_end = Float64(last(sol_cmp.t))
    t_start = max(t_ref_start, t_cmp_start)
    t_end = min(t_ref_end, t_cmp_end)
    if !(isfinite(t_start) && isfinite(t_end)) || t_end < t_start
        throw(ArgumentError("Cannot compare trajectories: incompatible time spans [$(t_ref_start), $(t_ref_end)] vs [$(t_cmp_start), $(t_cmp_end)]."))
    end
    sample_count = max(2, n_samples)
    sample_times = collect(range(t_start, t_end; length=sample_count))

    pos_rel_max = 0.0
    vel_rel_max = 0.0
    q_angle_max_rad = 0.0
    omega_rel_max = 0.0
    pos_rel_sum2 = 0.0
    vel_rel_sum2 = 0.0
    q_angle_sum2 = 0.0
    omega_rel_sum2 = 0.0
    sample_count_total = 0
    orientation_count_total = 0

    for t in sample_times
        u_ref = _solution_state_at(sol_ref, t)
        u_cmp = _solution_state_at(sol_cmp, t)
        for sat_idx in 1:n_sats
            sc_ref = u_ref.sc[sat_idx]
            sc_cmp = u_cmp.sc[sat_idx]
            pos_rel = _relative_vector_delta(sc_ref.pos, sc_cmp.pos)
            vel_rel = _relative_vector_delta(sc_ref.vel, sc_cmp.vel)
            pos_rel_max = max(pos_rel_max, pos_rel)
            vel_rel_max = max(vel_rel_max, vel_rel)
            pos_rel_sum2 += pos_rel^2
            vel_rel_sum2 += vel_rel^2
            sample_count_total += 1
            if orientation
                q_angle = _quaternion_angle_delta_rad(sc_ref.q, sc_cmp.q)
                omega_rel = _relative_vector_delta(sc_ref.ω, sc_cmp.ω)
                q_angle_max_rad = max(q_angle_max_rad, q_angle)
                omega_rel_max = max(omega_rel_max, omega_rel)
                q_angle_sum2 += q_angle^2
                omega_rel_sum2 += omega_rel^2
                orientation_count_total += 1
            end
        end
    end

    return (
        t_start=t_start,
        t_end=t_end,
        sample_count=sample_count,
        compared_state_count=sample_count_total,
        pos_rel_max=pos_rel_max,
        pos_rel_rms=sample_count_total > 0 ? sqrt(pos_rel_sum2 / sample_count_total) : missing,
        vel_rel_max=vel_rel_max,
        vel_rel_rms=sample_count_total > 0 ? sqrt(vel_rel_sum2 / sample_count_total) : missing,
        q_angle_max_rad=orientation ? q_angle_max_rad : missing,
        q_angle_rms_rad=orientation && orientation_count_total > 0 ? sqrt(q_angle_sum2 / orientation_count_total) : missing,
        omega_rel_max=orientation ? omega_rel_max : missing,
        omega_rel_rms=orientation && orientation_count_total > 0 ? sqrt(omega_rel_sum2 / orientation_count_total) : missing
    )
end

function _run_split_gate_solution(
    case::BenchmarkCase,
    profile_name::String;
    solver_mode::String,
    split_solver::Union{Nothing, String}=nothing
)
    GC.gc()
    args_run = deepcopy(case.args_template)
    started_ns = time_ns()
    try
        run_once = () -> _run_perf_simulation(
            args_run;
            return_solution=true,
            return_solver_metadata=true,
            profile_name=profile_name,
            solver_mode_override=solver_mode,
            split_imex_solver_override=split_solver,
            entry_target_count_override=case.entry_target_count_override
        )
        result = isempty(case.env_overrides) ? run_once() : withenv(case.env_overrides...) do
            run_once()
        end
        elapsed_s = (time_ns() - started_ns) / 1e9
        sol = result.solution
        solver_trace = result.solver_trace
        solver_sequence = isempty(solver_trace) ? missing : join([meta.solver for meta in solver_trace], "->")
        return (
            ok=true,
            elapsed_s=elapsed_s,
            success=_solve_success_for_case(sol, case),
            retcode=string(sol.retcode),
            solver_mode=result.solver_mode,
            solver_sequence=solver_sequence,
            solution=sol,
            error_text=missing
        )
    catch err
        if err isa InterruptException
            rethrow()
        end
        elapsed_s = (time_ns() - started_ns) / 1e9
        retcode = _solve_retcode_from_error(err)
        if ismissing(retcode)
            retcode = "Exception"
        end
        return (
            ok=false,
            elapsed_s=elapsed_s,
            success=false,
            retcode=String(retcode),
            solver_mode=missing,
            solver_sequence=missing,
            solution=nothing,
            error_text=_perf_error_text(err)
        )
    end
end

function _write_split_rollout_gate_report(
    path::String,
    spec::ProfileSpec,
    gate_df::DataFrame,
    gate_csv_path::String
)
    open(path, "w") do io
        println(io, "# Split-IMEX Rollout Gate (`$(spec.name)` profile)")
        println(io)
        println(io, "- Generated (UTC): $(string(now(UTC)))")
        println(io, "- Cases requested: `$(join(_split_rollout_case_names(), ", "))`")
        println(io, "- Split solvers: `$(join(_split_rollout_solver_variants(), ", "))`")
        println(io, "- Runtime slowdown ceiling: `$(_split_rollout_max_slowdown_ratio())x`")
        println(io, "- Position relative tolerance: `$(_split_rollout_pos_rel_tol())`")
        println(io, "- Velocity relative tolerance: `$(_split_rollout_vel_rel_tol())`")
        println(io, "- Quaternion angle tolerance [rad]: `$(_split_rollout_q_angle_tol_rad())`")
        println(io, "- Angular-rate relative tolerance: `$(_split_rollout_omega_rel_tol())`")
        println(io, "- Trajectory samples per comparison: `$(_split_rollout_sample_count())`")
        println(io, "- Enforce mode: `$(_split_rollout_enforce())`")
        println(io)
        println(io, "- Gate CSV: `$(gate_csv_path)`")
        println(io)
        pass_count = (nrow(gate_df) == 0 || !(:pass_all in names(gate_df))) ? 0 : count(Bool.(gate_df.pass_all))
        println(io, "- Gate pass count: `$(pass_count)/$(nrow(gate_df))`")
        println(io)
        println(io, "| Scenario | Split Solver | Baseline Retcode | Split Retcode | Runtime Ratio | Pos Rel Max | Vel Rel Max | Q Angle Max [rad] | Omega Rel Max | Pass Runtime | Pass Trajectory | Pass All |")
        println(io, "|---|---|---|---|---:|---:|---:|---:|---:|---|---|---|")
        for row in eachrow(gate_df)
            pass_traj = Bool(row.pass_pos) && Bool(row.pass_vel) && Bool(row.pass_q) && Bool(row.pass_omega)
            println(
                io,
                "| $(row.scenario) | $(row.split_solver) | $(row.baseline_retcode) | $(row.split_retcode) | $(_fmt(row.runtime_ratio; digits=4)) | " *
                "$(_fmt(row.pos_rel_max; digits=4)) | $(_fmt(row.vel_rel_max; digits=4)) | $(_fmt(row.q_angle_max_rad; digits=4)) | " *
                "$(_fmt(row.omega_rel_max; digits=4)) | $(row.pass_runtime) | $(pass_traj) | $(row.pass_all) |"
            )
        end
    end
    return nothing
end

function evaluate_split_rollout_gate(
    spec::ProfileSpec,
    cases::Vector{BenchmarkCase},
    outdir::String
)
    requested_names = _split_rollout_case_names()
    split_solvers = _split_rollout_solver_variants()
    case_pool = selected_cases(spec, cases)
    case_by_name = Dict(c.name => c for c in case_pool)
    max_slowdown = _split_rollout_max_slowdown_ratio()
    pos_tol = _split_rollout_pos_rel_tol()
    vel_tol = _split_rollout_vel_rel_tol()
    q_tol = _split_rollout_q_angle_tol_rad()
    omega_tol = _split_rollout_omega_rel_tol()
    sample_count = _split_rollout_sample_count()
    rows = NamedTuple[]

    for scenario_name in requested_names
        if !haskey(case_by_name, scenario_name)
            @warn "[split-rollout] requested scenario '$scenario_name' was not found in profile=$(spec.name); skipping."
            continue
        end
        case = case_by_name[scenario_name]
        baseline_case = _case_with_solver(case; solver_mode_override="auto_stiff", split_imex_solver_override=nothing)

        # Warm up both solver paths so the gate compares runtime behavior, not first-call compilation.
        run_warmup(baseline_case, 1, spec.name)
        for split_solver in split_solvers
            split_case = _case_with_solver(case; solver_mode_override="split_imex", split_imex_solver_override=split_solver)
            run_warmup(split_case, 1, spec.name)

            baseline_run = _run_split_gate_solution(
                baseline_case,
                spec.name;
                solver_mode="auto_stiff",
                split_solver=nothing
            )
            split_run = _run_split_gate_solution(
                split_case,
                spec.name;
                solver_mode="split_imex",
                split_solver=split_solver
            )

            runtime_ratio = (baseline_run.elapsed_s > 0.0) ? split_run.elapsed_s / baseline_run.elapsed_s : Inf
            pass_runtime = baseline_run.success && split_run.success && isfinite(runtime_ratio) && runtime_ratio <= max_slowdown

            pos_rel_max = missing
            vel_rel_max = missing
            q_angle_max_rad = missing
            omega_rel_max = missing
            compared_t_start = missing
            compared_t_end = missing
            compared_samples = missing
            pass_pos = false
            pass_vel = false
            pass_q = !case.args_template.mission_configuration.orientation_sim
            pass_omega = !case.args_template.mission_configuration.orientation_sim

            if baseline_run.success && split_run.success
                metrics = _trajectory_delta_metrics(
                    baseline_run.solution,
                    split_run.solution,
                    length(case.args_template.dynamics_model.spacecraft),
                    case.args_template.mission_configuration.orientation_sim;
                    n_samples=sample_count
                )
                pos_rel_max = metrics.pos_rel_max
                vel_rel_max = metrics.vel_rel_max
                q_angle_max_rad = metrics.q_angle_max_rad
                omega_rel_max = metrics.omega_rel_max
                compared_t_start = metrics.t_start
                compared_t_end = metrics.t_end
                compared_samples = metrics.sample_count
                pass_pos = metrics.pos_rel_max <= pos_tol
                pass_vel = metrics.vel_rel_max <= vel_tol
                if case.args_template.mission_configuration.orientation_sim
                    pass_q = !(metrics.q_angle_max_rad isa Missing) && metrics.q_angle_max_rad <= q_tol
                    pass_omega = !(metrics.omega_rel_max isa Missing) && metrics.omega_rel_max <= omega_tol
                end
            end

            pass_all = pass_runtime && pass_pos && pass_vel && pass_q && pass_omega
            push!(rows, (
                profile=spec.name,
                scenario=scenario_name,
                split_solver=split_solver,
                satellites=length(case.args_template.dynamics_model.spacecraft),
                orientation=case.args_template.mission_configuration.orientation_sim,
                baseline_elapsed_s=baseline_run.elapsed_s,
                split_elapsed_s=split_run.elapsed_s,
                runtime_ratio=runtime_ratio,
                max_slowdown_ratio=max_slowdown,
                pass_runtime=pass_runtime,
                pos_rel_max=pos_rel_max,
                vel_rel_max=vel_rel_max,
                q_angle_max_rad=q_angle_max_rad,
                omega_rel_max=omega_rel_max,
                pos_rel_tol=pos_tol,
                vel_rel_tol=vel_tol,
                q_angle_tol_rad=q_tol,
                omega_rel_tol=omega_tol,
                compared_t_start_s=compared_t_start,
                compared_t_end_s=compared_t_end,
                compared_samples=compared_samples,
                pass_pos=pass_pos,
                pass_vel=pass_vel,
                pass_q=pass_q,
                pass_omega=pass_omega,
                pass_all=pass_all,
                baseline_solver_mode=baseline_run.solver_mode,
                baseline_solver_sequence=baseline_run.solver_sequence,
                baseline_retcode=baseline_run.retcode,
                baseline_error=baseline_run.error_text,
                split_solver_mode=split_run.solver_mode,
                split_solver_sequence=split_run.solver_sequence,
                split_retcode=split_run.retcode,
                split_error=split_run.error_text
            ))
        end
    end

    gate_df = DataFrame(rows)
    stamp = Dates.format(now(UTC), dateformat"yyyymmdd_HHMMSS")
    gate_csv_path = joinpath(outdir, "split_rollout_gate_$(spec.name)_$(stamp).csv")
    gate_report_path = joinpath(outdir, "split_rollout_gate_$(spec.name)_$(stamp).md")
    CSV.write(gate_csv_path, gate_df)
    _write_split_rollout_gate_report(gate_report_path, spec, gate_df, gate_csv_path)

    if _split_rollout_enforce() && nrow(gate_df) > 0 && (:pass_all in names(gate_df)) && any(.!Bool.(gate_df.pass_all))
        failing = gate_df[.!gate_df.pass_all, :]
        summary = join(["$(row.scenario):$(row.split_solver)" for row in eachrow(failing)], ", ")
        error("Split rollout gate failed for $(nrow(failing)) configuration(s): $summary")
    end

    return (df=gate_df, csv_path=gate_csv_path, report_path=gate_report_path)
end

function _write_multirate_rollout_gate_report(
    path::String,
    spec::ProfileSpec,
    gate_df::DataFrame,
    gate_csv_path::String
)
    settings = _multirate_rollout_setting_snapshot()
    open(path, "w") do io
        println(io, "# Multirate Rollout Gate (`$(spec.name)` profile)")
        println(io)
        println(io, "- Generated (UTC): $(string(now(UTC)))")
        println(io, "- Cases requested: `$(join(_multirate_rollout_case_names(), ", "))`")
        println(io, "- Runtime slowdown ceiling: `$(_multirate_rollout_max_slowdown_ratio())x`")
        println(io, "- Position relative tolerance: `$(_multirate_rollout_pos_rel_tol())`")
        println(io, "- Velocity relative tolerance: `$(_multirate_rollout_vel_rel_tol())`")
        println(io, "- Quaternion angle tolerance [rad]: `$(_multirate_rollout_q_angle_tol_rad())`")
        println(io, "- Angular-rate relative tolerance: `$(_multirate_rollout_omega_rel_tol())`")
        println(io, "- Trajectory samples per comparison: `$(_multirate_rollout_sample_count())`")
        println(io, "- Enforce mode: `$(_multirate_rollout_enforce())`")
        println(io, "- Multirate slow solver: `$(settings.slow_solver)`")
        println(io, "- Multirate fast solver: `$(settings.fast_solver)`")
        println(io, "- Multirate slow dt [s]: `$(_fmt(settings.slow_dt_s; digits=4))`")
        println(io, "- Multirate fast substeps: `$(_fmt(settings.fast_substeps; digits=0))`")
        println(io)
        println(io, "- Gate CSV: `$(gate_csv_path)`")
        println(io)
        pass_count = (nrow(gate_df) == 0 || !(:pass_all in names(gate_df))) ? 0 : count(Bool.(gate_df.pass_all))
        println(io, "- Gate pass count: `$(pass_count)/$(nrow(gate_df))`")
        println(io)
        println(io, "| Scenario | Baseline Retcode | Multirate Retcode | Runtime Ratio | Pos Rel Max | Vel Rel Max | Q Angle Max [rad] | Omega Rel Max | Pass Runtime | Pass Trajectory | Pass All |")
        println(io, "|---|---|---|---:|---:|---:|---:|---:|---|---|---|")
        for row in eachrow(gate_df)
            pass_traj = Bool(row.pass_pos) && Bool(row.pass_vel) && Bool(row.pass_q) && Bool(row.pass_omega)
            println(
                io,
                "| $(row.scenario) | $(row.baseline_retcode) | $(row.multirate_retcode) | $(_fmt(row.runtime_ratio; digits=4)) | " *
                "$(_fmt(row.pos_rel_max; digits=4)) | $(_fmt(row.vel_rel_max; digits=4)) | $(_fmt(row.q_angle_max_rad; digits=4)) | " *
                "$(_fmt(row.omega_rel_max; digits=4)) | $(row.pass_runtime) | $(pass_traj) | $(row.pass_all) |"
            )
        end
    end
    return nothing
end

function evaluate_multirate_rollout_gate(
    spec::ProfileSpec,
    cases::Vector{BenchmarkCase},
    outdir::String
)
    requested_names = _multirate_rollout_case_names()
    case_pool = selected_cases(spec, cases)
    case_by_name = Dict(c.name => c for c in case_pool)
    max_slowdown = _multirate_rollout_max_slowdown_ratio()
    pos_tol = _multirate_rollout_pos_rel_tol()
    vel_tol = _multirate_rollout_vel_rel_tol()
    q_tol = _multirate_rollout_q_angle_tol_rad()
    omega_tol = _multirate_rollout_omega_rel_tol()
    sample_count = _multirate_rollout_sample_count()
    settings = _multirate_rollout_setting_snapshot()
    rows = NamedTuple[]

    for scenario_name in requested_names
        if !haskey(case_by_name, scenario_name)
            @warn "[multirate-rollout] requested scenario '$scenario_name' was not found in profile=$(spec.name); skipping."
            continue
        end
        case = case_by_name[scenario_name]
        baseline_case = _case_with_solver(case; solver_mode_override="auto_stiff", split_imex_solver_override=nothing)
        multirate_case = _case_with_solver(case; solver_mode_override="multirate", split_imex_solver_override=nothing)

        # Warm up both solver paths so the gate compares runtime behavior, not first-call compilation.
        run_warmup(baseline_case, 1, spec.name)
        run_warmup(multirate_case, 1, spec.name)

        baseline_run = _run_split_gate_solution(
            baseline_case,
            spec.name;
            solver_mode="auto_stiff",
            split_solver=nothing
        )
        multirate_run = _run_split_gate_solution(
            multirate_case,
            spec.name;
            solver_mode="multirate",
            split_solver=nothing
        )

        runtime_ratio = (baseline_run.elapsed_s > 0.0) ? multirate_run.elapsed_s / baseline_run.elapsed_s : Inf
        pass_runtime = baseline_run.success && multirate_run.success && isfinite(runtime_ratio) && runtime_ratio <= max_slowdown

        pos_rel_max = missing
        vel_rel_max = missing
        q_angle_max_rad = missing
        omega_rel_max = missing
        compared_t_start = missing
        compared_t_end = missing
        compared_samples = missing
        pass_pos = false
        pass_vel = false
        pass_q = !case.args_template.mission_configuration.orientation_sim
        pass_omega = !case.args_template.mission_configuration.orientation_sim

        if baseline_run.success && multirate_run.success
            metrics = _trajectory_delta_metrics(
                baseline_run.solution,
                multirate_run.solution,
                length(case.args_template.dynamics_model.spacecraft),
                case.args_template.mission_configuration.orientation_sim;
                n_samples=sample_count
            )
            pos_rel_max = metrics.pos_rel_max
            vel_rel_max = metrics.vel_rel_max
            q_angle_max_rad = metrics.q_angle_max_rad
            omega_rel_max = metrics.omega_rel_max
            compared_t_start = metrics.t_start
            compared_t_end = metrics.t_end
            compared_samples = metrics.sample_count
            pass_pos = metrics.pos_rel_max <= pos_tol
            pass_vel = metrics.vel_rel_max <= vel_tol
            if case.args_template.mission_configuration.orientation_sim
                pass_q = !(metrics.q_angle_max_rad isa Missing) && metrics.q_angle_max_rad <= q_tol
                pass_omega = !(metrics.omega_rel_max isa Missing) && metrics.omega_rel_max <= omega_tol
            end
        end

        pass_all = pass_runtime && pass_pos && pass_vel && pass_q && pass_omega
        push!(rows, (
            profile=spec.name,
            scenario=scenario_name,
            satellites=length(case.args_template.dynamics_model.spacecraft),
            orientation=case.args_template.mission_configuration.orientation_sim,
            baseline_elapsed_s=baseline_run.elapsed_s,
            multirate_elapsed_s=multirate_run.elapsed_s,
            runtime_ratio=runtime_ratio,
            max_slowdown_ratio=max_slowdown,
            pass_runtime=pass_runtime,
            pos_rel_max=pos_rel_max,
            vel_rel_max=vel_rel_max,
            q_angle_max_rad=q_angle_max_rad,
            omega_rel_max=omega_rel_max,
            pos_rel_tol=pos_tol,
            vel_rel_tol=vel_tol,
            q_angle_tol_rad=q_tol,
            omega_rel_tol=omega_tol,
            compared_t_start_s=compared_t_start,
            compared_t_end_s=compared_t_end,
            compared_samples=compared_samples,
            pass_pos=pass_pos,
            pass_vel=pass_vel,
            pass_q=pass_q,
            pass_omega=pass_omega,
            pass_all=pass_all,
            multirate_slow_solver=settings.slow_solver,
            multirate_fast_solver=settings.fast_solver,
            multirate_slow_dt_s=settings.slow_dt_s,
            multirate_fast_substeps=settings.fast_substeps,
            baseline_solver_mode=baseline_run.solver_mode,
            baseline_solver_sequence=baseline_run.solver_sequence,
            baseline_retcode=baseline_run.retcode,
            baseline_error=baseline_run.error_text,
            multirate_solver_mode=multirate_run.solver_mode,
            multirate_solver_sequence=multirate_run.solver_sequence,
            multirate_retcode=multirate_run.retcode,
            multirate_error=multirate_run.error_text
        ))
    end

    gate_df = DataFrame(rows)
    stamp = Dates.format(now(UTC), dateformat"yyyymmdd_HHMMSS")
    gate_csv_path = joinpath(outdir, "multirate_rollout_gate_$(spec.name)_$(stamp).csv")
    gate_report_path = joinpath(outdir, "multirate_rollout_gate_$(spec.name)_$(stamp).md")
    CSV.write(gate_csv_path, gate_df)
    _write_multirate_rollout_gate_report(gate_report_path, spec, gate_df, gate_csv_path)

    if _multirate_rollout_enforce() && nrow(gate_df) > 0 && (:pass_all in names(gate_df)) && any(.!Bool.(gate_df.pass_all))
        failing = gate_df[.!gate_df.pass_all, :]
        summary = join([String(row.scenario) for row in eachrow(failing)], ", ")
        error("Multirate rollout gate failed for $(nrow(failing)) configuration(s): $summary")
    end

    return (df=gate_df, csv_path=gate_csv_path, report_path=gate_report_path)
end

@inline _plot_label(name::AbstractString) = replace(name, "_" => " ")
@inline _plot_axis_label(name::AbstractString) = replace(name, "_" => "\n")
@inline _plot_number(v) = v isa Missing ? NaN : Float64(v)

