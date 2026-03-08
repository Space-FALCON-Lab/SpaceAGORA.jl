@inline function _resolve_component_tolerance(component_tol::Float64, fallback_tol::Float64, name::String)::Float64
    if component_tol < 0.0
        throw(ArgumentError("$name must be >= 0.0, got $component_tol"))
    end
    return component_tol == 0.0 ? fallback_tol : component_tol
end

@inline function _requires_componentwise_tolerances(args)::Bool
    tol = args.integration_tolerances
    return args.mission_configuration.orientation_sim ||
           tol.reltol_mass != 0.0 || tol.abstol_mass != 0.0 ||
           tol.reltol_heat_load != 0.0 || tol.abstol_heat_load != 0.0
end

function _build_solver_tolerances(u_state::ComponentVector, args)
    tol = args.integration_tolerances
    if !_requires_componentwise_tolerances(args)
        return tol.reltol_orbit, tol.abstol_orbit
    end

    reltol_mass = _resolve_component_tolerance(tol.reltol_mass, tol.reltol_orbit, "reltol_mass")
    abstol_mass = _resolve_component_tolerance(tol.abstol_mass, tol.abstol_orbit, "abstol_mass")
    reltol_heat = _resolve_component_tolerance(tol.reltol_heat_load, tol.reltol_orbit, "reltol_heat_load")
    abstol_heat = _resolve_component_tolerance(tol.abstol_heat_load, tol.abstol_orbit, "abstol_heat_load")
    reltol_ω = _resolve_component_tolerance(tol.reltol_angular_rate, tol.reltol_orbit, "reltol_angular_rate")
    abstol_ω = _resolve_component_tolerance(tol.abstol_angular_rate, tol.abstol_orbit, "abstol_angular_rate")

    reltol_state = copy(u_state)
    abstol_state = copy(u_state)
    reltol_state .= tol.reltol_orbit
    abstol_state .= tol.abstol_orbit
    @inbounds for i in eachindex(reltol_state.sc)
        reltol_state.sc[i].mass = reltol_mass
        abstol_state.sc[i].mass = abstol_mass
        reltol_state.sc[i].heat_loads .= reltol_heat
        abstol_state.sc[i].heat_loads .= abstol_heat
    end

    if args.mission_configuration.orientation_sim
        @inbounds for i in eachindex(reltol_state.sc)
            reltol_state.sc[i].ω .= reltol_ω
            abstol_state.sc[i].ω .= abstol_ω
        end
    end

    if args.mission_configuration.orientation_sim
        @inbounds for i in eachindex(reltol_state.sc)
        reltol_state.sc[i].q .= tol.reltol_quaternion
        abstol_state.sc[i].q .= tol.abstol_quaternion
        end
    end
    return reltol_state, abstol_state
end

@inline function _solver_policy_mode()::Symbol
    mode = lowercase(strip(_engine_env_get("SPACEAGORA_SOLVER_MODE", "tsit5")))
    if mode in ("tsit5", "default")
        return :tsit5
    elseif mode in ("symplectic", "kahanli8", "verlet")
        return :symplectic
    elseif mode in ("auto_stiff", "auto-stiff", "autostiff", "auto")
        return :auto_stiff
    elseif mode in ("rodas5p", "rodas", "stiff")
        return :rodas5p
    elseif mode in ("split_imex", "split-imex", "split", "imex")
        return :split_imex
    elseif mode in ("multirate", "multirate_split", "split_multirate", "mr")
        return :multirate
    end
    throw(ArgumentError(
        "Unsupported SPACEAGORA_SOLVER_MODE='$mode'. Use one of: tsit5, symplectic, auto_stiff, rodas5p, split_imex, multirate."
    ))
end

@inline function _retcode_is_stiff_symptom(retcode)::Bool
    rc = string(retcode)
    return rc in ("Unstable", "DtLessThanMin", "MaxIters", "InitialFailure")
end

@inline function _auto_stiff_switched(sol)::Bool
    hasproperty(sol, :alg_choice) || return false
    choices = getproperty(sol, :alg_choice)
    isempty(choices) && return false
    first_choice = first(choices)
    @inbounds for choice in choices
        if choice != first_choice
            return true
        end
    end
    return false
end

@inline function _solver_maxiters()::Union{Nothing, Int}
    raw = strip(_engine_env_get("SPACEAGORA_SOLVER_MAXITERS", ""))
    isempty(raw) && return nothing
    parsed = try
        parse(Int, raw)
    catch
        throw(ArgumentError("SPACEAGORA_SOLVER_MAXITERS must be an integer, got '$raw'."))
    end
    parsed > 0 || throw(ArgumentError("SPACEAGORA_SOLVER_MAXITERS must be > 0, got $parsed."))
    return parsed
end

@inline function _symplectic_fixed_dt_s(args)::Float64
    raw = strip(_engine_env_get("SPACEAGORA_SYMPLECTIC_DT_S", ""))
    dt = if isempty(raw)
        args.integration_tolerances.dt_max_orbit
    else
        parsed = try
            parse(Float64, raw)
        catch
            throw(ArgumentError("SPACEAGORA_SYMPLECTIC_DT_S must be a number, got '$raw'."))
        end
        parsed
    end
    dt > 0.0 || throw(ArgumentError("SPACEAGORA_SYMPLECTIC_DT_S must be > 0.0, got $dt."))
    return dt
end

@inline function _symplectic_conservative_eligible(args)::Bool
    args.mission_configuration.orientation_sim && return false
    length(args.dynamics_model.spacecraft) == 1 || return false
    isempty(args.control_model.control_effectors) || return false
    isempty(args.guidance_model.guidance_effectors) || return false
    isempty(args.navigation_model.navigation_effectors) || return false
    effectors = args.dynamics_model.dynamic_effectors
    length(effectors) == 1 || return false
    return first(effectors) isa SimulationModel.InverseSquaredGravityModel
end

@inline function _split_imex_solver_spec()
    mode = lowercase(strip(_engine_env_get("SPACEAGORA_SPLIT_IMEX_SOLVER", "kencarp4")))
    if mode in ("kencarp4", "ken4", "default")
        return (alg=KenCarp4(autodiff=AutoFiniteDiff()), label="KenCarp4")
    elseif mode in ("kencarp47", "ken47")
        return (alg=KenCarp47(autodiff=AutoFiniteDiff()), label="KenCarp47")
    elseif mode in ("kencarp58", "ken58")
        return (alg=KenCarp58(autodiff=AutoFiniteDiff()), label="KenCarp58")
    end
    throw(ArgumentError(
        "Unsupported SPACEAGORA_SPLIT_IMEX_SOLVER='$mode'. Use one of: kencarp4, kencarp47, kencarp58."
    ))
end

@inline function _multirate_fast_substeps()::Int
    raw = strip(_engine_env_get("SPACEAGORA_MULTIRATE_FAST_SUBSTEPS", "8"))
    parsed = try
        parse(Int, raw)
    catch
        throw(ArgumentError("SPACEAGORA_MULTIRATE_FAST_SUBSTEPS must be an integer, got '$raw'."))
    end
    parsed > 0 || throw(ArgumentError("SPACEAGORA_MULTIRATE_FAST_SUBSTEPS must be > 0, got $parsed."))
    return parsed
end

@inline function _multirate_slow_dt_s(args)::Float64
    default_dt = min(args.integration_tolerances.dt_max_orbit, 2.0)
    raw = strip(_engine_env_get("SPACEAGORA_MULTIRATE_SLOW_DT_S", ""))
    dt = if isempty(raw)
        default_dt
    else
        parsed = try
            parse(Float64, raw)
        catch
            throw(ArgumentError("SPACEAGORA_MULTIRATE_SLOW_DT_S must be a number, got '$raw'."))
        end
        parsed
    end
    dt > 0.0 || throw(ArgumentError("SPACEAGORA_MULTIRATE_SLOW_DT_S must be > 0.0, got $dt."))
    return min(dt, args.integration_tolerances.dt_max_orbit)
end

@inline function _multirate_solver_spec(env_name::String, default_mode::String)
    mode = lowercase(strip(_engine_env_get(env_name, default_mode)))
    if mode in ("tsit5", "tsit", "default")
        return (alg=Tsit5(), label="Tsit5", auto_switch_capable=false)
    elseif mode in ("auto_stiff", "auto-stiff", "autostiff", "auto")
        return (
            alg=AutoTsit5(Rodas5P(autodiff=AutoFiniteDiff())),
            label="AutoTsit5(Rodas5P)",
            auto_switch_capable=true
        )
    elseif mode in ("rodas5p", "rodas", "stiff")
        return (alg=Rodas5P(autodiff=AutoFiniteDiff()), label="Rodas5P", auto_switch_capable=false)
    elseif mode in ("kencarp4", "ken4")
        return (alg=KenCarp4(autodiff=AutoFiniteDiff()), label="KenCarp4", auto_switch_capable=false)
    end
    throw(ArgumentError(
        "Unsupported $(env_name)='$mode'. Use one of: tsit5, auto_stiff, rodas5p, kencarp4."
    ))
end

@inline _multirate_slow_solver_spec() = _multirate_solver_spec("SPACEAGORA_MULTIRATE_SLOW_SOLVER", "tsit5")
@inline _multirate_fast_solver_spec() = _multirate_solver_spec("SPACEAGORA_MULTIRATE_FAST_SOLVER", "auto_stiff")

@inline function _solve_with_explicit_solver(prob, args, alg, reltol_tol, abstol_tol; dtmax_override::Union{Nothing, Float64}=nothing)
    maxiters = _solver_maxiters()
    dtmax_use = isnothing(dtmax_override) ? args.integration_tolerances.dt_max_orbit : dtmax_override
    dtmax_use > 0.0 || throw(ArgumentError("Solver dtmax must be > 0.0, got $dtmax_use."))
    if maxiters === nothing
        return solve(
            prob,
            alg;
            reltol=reltol_tol,
            abstol=abstol_tol,
            dtmax=dtmax_use
        )
    end
    return solve(
        prob,
        alg;
        reltol=reltol_tol,
        abstol=abstol_tol,
        dtmax=dtmax_use,
        maxiters=maxiters
    )
end

@inline function _solve_with_fixed_step_solver(prob, alg, dt_s::Float64)
    maxiters = _solver_maxiters()
    if maxiters === nothing
        return solve(prob, alg; dt=dt_s)
    end
    return solve(prob, alg; dt=dt_s, maxiters=maxiters)
end

@inline function _is_partitioned_second_order_problem(prob)::Bool
    hasproperty(prob, :problem_type) || return false
    return getproperty(prob, :problem_type) isa SecondOrderODEProblem
end

@inline function _split_subproblem(prob, f, u, tspan)
    return ODEProblem(f, u, tspan, prob.p; prob.kwargs...)
end

function _solve_with_multirate_solver(prob, args, reltol_tol, abstol_tol)
    if !(hasproperty(prob.f, :f1) && hasproperty(prob.f, :f2))
        throw(ArgumentError("SPACEAGORA_SOLVER_MODE=multirate requires a split problem with f1/f2 components."))
    end

    t_start = Float64(first(prob.tspan))
    t_end = Float64(last(prob.tspan))
    if t_end <= t_start
        sol = _solve_with_explicit_solver(prob, args, Tsit5(), reltol_tol, abstol_tol)
        return sol, (
            slow_solver="Tsit5",
            fast_solver="Tsit5",
            macro_steps=0,
            fast_substeps=0,
            slow_dt_s=0.0,
            fast_dt_s=0.0,
            auto_switch_events=0
        )
    end

    slow_spec = _multirate_slow_solver_spec()
    fast_spec = _multirate_fast_solver_spec()
    fast_substeps = _multirate_fast_substeps()
    slow_dt_s = _multirate_slow_dt_s(args)
    fast_dt_s = slow_dt_s / fast_substeps

    t_cursor = t_start
    u_cursor = deepcopy(prob.u0)
    final_sol = nothing
    macro_steps = 0
    auto_switch_events = 0

    while t_cursor < t_end
        t_next = min(t_cursor + slow_dt_s, t_end)
        macro_steps += 1

        # Strang splitting: fast half-step -> slow full-step -> fast half-step.
        segment_dt = t_next - t_cursor
        half_dt = 0.5 * segment_dt
        t_half = t_cursor + half_dt

        if half_dt > 0.0
            fast_prob_pre = _split_subproblem(prob, prob.f.f2, u_cursor, (t_cursor, t_half))
            sol_fast_pre = _solve_with_explicit_solver(
                fast_prob_pre,
                args,
                fast_spec.alg,
                reltol_tol,
                abstol_tol;
                dtmax_override=min(fast_dt_s, half_dt)
            )
            if fast_spec.auto_switch_capable && _auto_stiff_switched(sol_fast_pre)
                auto_switch_events += 1
            end
            if !SciMLBase.successful_retcode(sol_fast_pre.retcode)
                return sol_fast_pre, (
                    slow_solver=slow_spec.label,
                    fast_solver=fast_spec.label,
                    macro_steps=macro_steps,
                    fast_substeps=fast_substeps,
                    slow_dt_s=slow_dt_s,
                    fast_dt_s=fast_dt_s,
                    auto_switch_events=auto_switch_events
                )
            end
            u_cursor = deepcopy(sol_fast_pre.u[end])
            final_sol = sol_fast_pre
        end

        slow_prob = _split_subproblem(prob, prob.f.f1, u_cursor, (t_cursor, t_next))
        sol_slow = _solve_with_explicit_solver(
            slow_prob,
            args,
            slow_spec.alg,
            reltol_tol,
            abstol_tol;
            dtmax_override=segment_dt
        )
        if slow_spec.auto_switch_capable && _auto_stiff_switched(sol_slow)
            auto_switch_events += 1
        end
        if !SciMLBase.successful_retcode(sol_slow.retcode)
            return sol_slow, (
                slow_solver=slow_spec.label,
                fast_solver=fast_spec.label,
                macro_steps=macro_steps,
                fast_substeps=fast_substeps,
                slow_dt_s=slow_dt_s,
                fast_dt_s=fast_dt_s,
                auto_switch_events=auto_switch_events
            )
        end
        u_cursor = deepcopy(sol_slow.u[end])
        final_sol = sol_slow

        if half_dt > 0.0
            fast_prob_post = _split_subproblem(prob, prob.f.f2, u_cursor, (t_half, t_next))
            sol_fast_post = _solve_with_explicit_solver(
                fast_prob_post,
                args,
                fast_spec.alg,
                reltol_tol,
                abstol_tol;
                dtmax_override=min(fast_dt_s, half_dt)
            )
            if fast_spec.auto_switch_capable && _auto_stiff_switched(sol_fast_post)
                auto_switch_events += 1
            end
            if !SciMLBase.successful_retcode(sol_fast_post.retcode)
                return sol_fast_post, (
                    slow_solver=slow_spec.label,
                    fast_solver=fast_spec.label,
                    macro_steps=macro_steps,
                    fast_substeps=fast_substeps,
                    slow_dt_s=slow_dt_s,
                    fast_dt_s=fast_dt_s,
                    auto_switch_events=auto_switch_events
                )
            end
            u_cursor = deepcopy(sol_fast_post.u[end])
            final_sol = sol_fast_post
        end

        t_cursor = t_next
    end

    return final_sol, (
        slow_solver=slow_spec.label,
        fast_solver=fast_spec.label,
        macro_steps=macro_steps,
        fast_substeps=fast_substeps,
        slow_dt_s=slow_dt_s,
        fast_dt_s=fast_dt_s,
        auto_switch_events=auto_switch_events
    )
end

function _solve_with_solver_policy(prob, args, reltol_tol, abstol_tol)
    mode = _solver_policy_mode()
    if mode == :symplectic
        _symplectic_conservative_eligible(args) || throw(ArgumentError(
            "SPACEAGORA_SOLVER_MODE=symplectic currently requires a single-spacecraft, translational-only inverse-squared gravity configuration with no control/guidance/navigation effectors."
        ))
        _is_partitioned_second_order_problem(prob) || throw(ArgumentError(
            "SPACEAGORA_SOLVER_MODE=symplectic requires a partitioned SecondOrderODEProblem. The typed run_simulation path still builds a first-order ODEProblem, so use tsit5/auto_stiff there until a partitioned runtime path is added."
        ))
        sol = _solve_with_fixed_step_solver(prob, KahanLi8(), _symplectic_fixed_dt_s(args))
        return sol, (
            solver="KahanLi8(Symplectic)",
            initial_solver="KahanLi8",
            fallback_used=false,
            trigger_retcode=missing
        )
    end

    if mode == :rodas5p
        sol = _solve_with_explicit_solver(prob, args, Rodas5P(autodiff=AutoFiniteDiff()), reltol_tol, abstol_tol)
        return sol, (
            solver="Rodas5P",
            initial_solver="Rodas5P",
            fallback_used=false,
            trigger_retcode=missing
        )
    end

    if mode == :auto_stiff
        # True stiffness-aware autoswitching handled internally by OrdinaryDiffEq.
        # This replaces the manual "retry with Rodas5P on Tsit5 failure" policy.
        autoswitch_alg = AutoTsit5(Rodas5P(autodiff=AutoFiniteDiff()))
        sol = _solve_with_explicit_solver(prob, args, autoswitch_alg, reltol_tol, abstol_tol)
        switched = _auto_stiff_switched(sol)
        return sol, (
            solver="AutoTsit5(Rodas5P)",
            initial_solver="AutoTsit5",
            fallback_used=switched,
            trigger_retcode=switched ? "internal_autoswitch" : missing
        )
    end

    if mode == :split_imex
        split_solver = _split_imex_solver_spec()
        sol = _solve_with_explicit_solver(prob, args, split_solver.alg, reltol_tol, abstol_tol)
        return sol, (
            solver="$(split_solver.label)(IMEX)",
            initial_solver=split_solver.label,
            fallback_used=false,
            trigger_retcode=missing
        )
    end

    if mode == :multirate
        sol, multirate_meta = _solve_with_multirate_solver(prob, args, reltol_tol, abstol_tol)
        switched = multirate_meta.auto_switch_events > 0
        return sol, (
            solver="Multirate(Strang; slow=$(multirate_meta.slow_solver), fast=$(multirate_meta.fast_solver))",
            initial_solver=multirate_meta.slow_solver,
            fallback_used=switched,
            trigger_retcode=switched ? "internal_autoswitch" : missing
        )
    end

    tsit_sol = _solve_with_explicit_solver(prob, args, Tsit5(), reltol_tol, abstol_tol)
    return tsit_sol, (
        solver="Tsit5",
        initial_solver="Tsit5",
        fallback_used=false,
        trigger_retcode=missing
    )
end
