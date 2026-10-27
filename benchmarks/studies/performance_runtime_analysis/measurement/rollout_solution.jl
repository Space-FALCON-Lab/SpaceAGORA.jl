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

