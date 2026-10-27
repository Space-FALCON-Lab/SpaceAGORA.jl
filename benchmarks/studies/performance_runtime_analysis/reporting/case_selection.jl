# Case records supply the named fields used here; no simulation package is loaded.
@inline function selected_cases(spec::ProfileSpec, cases::AbstractVector)
    case_filter = _perf_case_name_filter()
    selected = if isempty(case_filter)
        spec.name == "full" ? cases : [c for c in cases if c.run_in_quick]
    else
        by_name = Dict(c.name => c for c in cases)
        [by_name[name] for name in case_filter if haskey(by_name, name)]
    end
    if _exclude_entry_scenarios()
        selected = [c for c in selected if c.category != "entry"]
    end
    if _perf_smoke_mode()
        return selected[1:min(end, 3)]
    end
    return selected
end

@inline function _case_with_solver(
    case;
    solver_mode_override::Union{Nothing, String}=nothing,
    split_imex_solver_override::Union{Nothing, String}=nothing
)
    return typeof(case)(
        name=case.name,
        category=case.category,
        description=case.description,
        args_template=case.args_template,
        run_in_quick=case.run_in_quick,
        solver_mode_override=solver_mode_override,
        split_imex_solver_override=split_imex_solver_override,
        entry_target_count_override=case.entry_target_count_override,
        env_overrides=copy(case.env_overrides)
    )
end

