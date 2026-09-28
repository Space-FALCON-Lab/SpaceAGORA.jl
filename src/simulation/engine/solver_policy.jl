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
    orientation_sim = args.mission_configuration.orientation_sim
    reltol_q = tol.reltol_quaternion
    abstol_q = tol.abstol_quaternion
    # Single pass over all spacecraft to set mass, heat load, and (when enabled)
    # orientation tolerances — avoids iterating the spacecraft array 2–3 times.
    @inbounds for i in eachindex(reltol_state.sc)
        reltol_state.sc[i].mass = reltol_mass
        abstol_state.sc[i].mass = abstol_mass
        reltol_state.sc[i].heat_loads .= reltol_heat
        abstol_state.sc[i].heat_loads .= abstol_heat
        if orientation_sim
            reltol_state.sc[i].ω .= reltol_ω
            abstol_state.sc[i].ω .= abstol_ω
            reltol_state.sc[i].q .= reltol_q
            abstol_state.sc[i].q .= abstol_q
        end
    end
    return reltol_state, abstol_state
end

@inline _solver_policy_mode(cfg::SolverConfig)::Symbol = cfg.solver_mode
@inline _solver_policy_mode()::Symbol = _solver_policy_mode(_active_solver_config())

@inline function _retcode_is_stiff_symptom(retcode)::Bool
    # Convert once to Symbol (zero-allocation for Symbol/ReturnCode inputs,
    # one allocation for String inputs — but String inputs are test-only).
    # This avoids the heap String allocation of `string(retcode) in (...)`.
    sym = Symbol(retcode)
    return sym === :Unstable     ||
           sym === :DtLessThanMin ||
           sym === :MaxIters      ||
           sym === :InitialFailure
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

@inline _solver_maxiters(cfg::SolverConfig)::Union{Nothing, Int} = cfg.maxiters
@inline _active_solver_config()::SolverConfig = begin
    active_config = _engine_active_config_ref[]
    active_config === nothing ? simulation_engine_config_from_env(; solver_strict=true).solver : active_config.solver
end
@inline _solver_maxiters()::Union{Nothing, Int} = _solver_maxiters(_active_solver_config())

@inline function _symplectic_fixed_dt_s(cfg::SolverConfig, args)::Float64
    dt = isnothing(cfg.symplectic_dt_s) ? args.integration_tolerances.dt_max_orbit : cfg.symplectic_dt_s
    dt > 0.0 || throw(ArgumentError("SolverConfig.symplectic_dt_s must be > 0.0, got $dt."))
    return dt
end
@inline _symplectic_fixed_dt_s(args)::Float64 = _symplectic_fixed_dt_s(_active_solver_config(), args)

@inline function _gravity_backbone_fixed_dt_s(cfg::SolverConfig, args)::Float64
    dt = isnothing(cfg.gravity_backbone_dt_s) ? args.integration_tolerances.dt_max_orbit : cfg.gravity_backbone_dt_s
    dt > 0.0 || throw(ArgumentError("SolverConfig.gravity_backbone_dt_s must be > 0.0, got $dt."))
    return dt
end
@inline _gravity_backbone_fixed_dt_s(args)::Float64 = _gravity_backbone_fixed_dt_s(_active_solver_config(), args)

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

@inline _auto_stiff_smooth_gravity_tsit5_enabled(cfg::SolverConfig)::Bool = cfg.auto_stiff_gravity_tsit5

@inline function _auto_stiff_smooth_gravity_effector(effector)::Bool
    return effector isa SimulationModel.InverseSquaredGravityModel ||
           effector isa SimulationModel.InverseSquaredJ2GravityModel ||
           effector isa SimulationModel.GravitationalHarmonicsModel ||
           effector isa SimulationModel.NBodyGravityModel
end

@inline function _auto_stiff_smooth_gravity_reject_reason(cfg::SolverConfig, args)::Union{Nothing, String}
    _auto_stiff_smooth_gravity_tsit5_enabled(cfg) || return "SolverConfig.auto_stiff_gravity_tsit5=false disables smooth-gravity Tsit5 routing."
    args.mission_configuration.orientation_sim && return "orientation_sim=true can couple attitude states into the RHS."
    isempty(args.control_model.control_effectors) || return "control effectors are active."
    isempty(args.guidance_model.guidance_effectors) || return "guidance effectors are active."
    isempty(args.navigation_model.navigation_effectors) || return "navigation effectors are active."

    effectors = args.dynamics_model.dynamic_effectors
    isempty(effectors) && return "no dynamic effectors are active."
    @inbounds for effector in effectors
        _auto_stiff_smooth_gravity_effector(effector) || return "$(nameof(typeof(effector))) is not a supported smooth-gravity effector."
        req = SimulationModel.environment_requirements(effector)
        req.atmosphere && return "$(nameof(typeof(effector))) requires atmosphere samples."
        req.solar && return "$(nameof(typeof(effector))) requires solar samples."
    end
    return nothing
end

@inline _auto_stiff_smooth_gravity_eligible(cfg::SolverConfig, args)::Bool = isnothing(_auto_stiff_smooth_gravity_reject_reason(cfg, args))
@inline _auto_stiff_smooth_gravity_eligible(args)::Bool = _auto_stiff_smooth_gravity_eligible(_active_solver_config(), args)

# Consecutive stiff-detection events required before AutoTsit5 commits to Rodas5P.
# The OrdinaryDiffEq default is 5, which can trigger during brief atmospheric passes
# where the entry/exit density gradient looks locally stiff.  A higher value requires
# sustained stiffness across many steps before switching, avoiding unnecessary implicit
# solves during aerobraking passages.
@inline _auto_stiff_switch_max(cfg::SolverConfig)::Int = cfg.auto_stiff_switch_max
@inline _auto_stiff_switch_max()::Int = _auto_stiff_switch_max(_active_solver_config())

@inline function _gravity_backbone_structure_validated(effector)::Symbol
    structure = SimulationModel.gravity_backbone_structure(effector)
    if structure === :unsupported || structure === :position_only_static_gravity
        return structure
    end
    throw(ArgumentError(
        "gravity_backbone_structure($(nameof(typeof(effector)))) must return :unsupported or :position_only_static_gravity, got $(repr(structure))."
    ))
end

@inline function _gravity_backbone_kick_structure_validated(effector)::Symbol
    structure = SimulationModel.gravity_backbone_kick_structure(effector)
    if structure === :unsupported || structure === :velocity_kick_explicit
        return structure
    end
    throw(ArgumentError(
        "gravity_backbone_kick_structure($(nameof(typeof(effector)))) must return :unsupported or :velocity_kick_explicit, got $(repr(structure))."
    ))
end

@inline function _gravity_backbone_has_kicks(args)::Bool
    @inbounds for effector in args.dynamics_model.dynamic_effectors
        if _gravity_backbone_kick_structure_validated(effector) === :velocity_kick_explicit
            return true
        end
    end
    return false
end

@inline function _gravity_backbone_time_reached(t_now::Float64, t_target::Float64)::Bool
    tol = max(1e-9, 64 * eps(Float64) * max(1.0, abs(t_target)))
    return abs(t_now - t_target) <= tol
end

@inline function _gravity_backbone_reject_reason(args)::Union{Nothing, String}
    args.mission_configuration.orientation_sim && return "SPACEAGORA_SOLVER_MODE=gravity_backbone_split requires orientation_sim=false."
    isempty(args.control_model.control_effectors) || return "SPACEAGORA_SOLVER_MODE=gravity_backbone_split requires no control effectors."
    isempty(args.guidance_model.guidance_effectors) || return "SPACEAGORA_SOLVER_MODE=gravity_backbone_split requires no guidance effectors."
    isempty(args.navigation_model.navigation_effectors) || return "SPACEAGORA_SOLVER_MODE=gravity_backbone_split requires no navigation effectors."
    effectors = args.dynamics_model.dynamic_effectors
    isempty(effectors) && return "SPACEAGORA_SOLVER_MODE=gravity_backbone_split requires at least one supported gravity effector."
    has_core = false
    @inbounds for effector in effectors
        core_structure = _gravity_backbone_structure_validated(effector)
        kick_structure = _gravity_backbone_kick_structure_validated(effector)
        if core_structure === :position_only_static_gravity
            has_core = true
        elseif kick_structure !== :velocity_kick_explicit
            return "SPACEAGORA_SOLVER_MODE=gravity_backbone_split does not support $(nameof(typeof(effector)))."
        end
        req = SimulationModel.environment_requirements(effector)
        if core_structure === :position_only_static_gravity
            req.atmosphere && return "SPACEAGORA_SOLVER_MODE=gravity_backbone_split rejects atmosphere-dependent effectors."
            req.solar && return "SPACEAGORA_SOLVER_MODE=gravity_backbone_split rejects solar/SRP-dependent gravity core effectors."
            isempty(req.third_body_names) || return "SPACEAGORA_SOLVER_MODE=gravity_backbone_split rejects third-body/ephemeris-dependent gravity core effectors."
        else
            req.atmosphere && return "SPACEAGORA_SOLVER_MODE=gravity_backbone_split rejects atmosphere-dependent effectors."
            req.planet_frame && return "SPACEAGORA_SOLVER_MODE=gravity_backbone_split does not support planet-frame-dependent perturbation kicks."
        end
    end
    has_core || return "SPACEAGORA_SOLVER_MODE=gravity_backbone_split requires at least one supported gravity core effector."
    return nothing
end

@inline _gravity_backbone_eligible(args)::Bool = isnothing(_gravity_backbone_reject_reason(args))

# KLUFactorization requires a sparse W, which only exists when the component
# function carries a jac_prototype; on a dense W it produces wrong Newton
# corrections rather than an error, so it is opt-in per problem exactly as in
# _rodas5p_alg below.
#
# A per-satellite dense-LU solver was built and measured against KLU here and
# was strictly worse, so it is deliberately not an option. Block LU does win the
# factorization (0.091 ms against 0.31 ms on a 256-block, n=2560 W), but it
# loses the solve badly -- 0.0495 ms against 0.0226 ms -- because 256 small
# LAPACK ldiv! calls cost more per call than one sparse triangular solve. This
# workload is solve-dominated by roughly 42:1 (nsolve=543 against nw=13 on a
# representative constellation run), so the solve is what counts and KLU wins
# it. Measuring factorize+solve as a unit is what makes block LU look 3.4x
# faster; that ratio does not survive contact with the integrator's actual mix.
@inline _sparse_linsolve_or_default(sparse_jac::Bool) = sparse_jac ? KLUFactorization() : nothing

"""Return whether a problem component function carries a sparse Jacobian prototype."""
@inline function _has_sparse_jac_prototype(f)::Bool
    hasproperty(f, :jac_prototype) || return false
    return getproperty(f, :jac_prototype) isa SparseMatrixCSC
end

"""Return whether split component `name` (`:f1`/`:f2`) of `prob` carries a sparse prototype."""
@inline _split_component_has_sparse_jac(prob, name::Symbol)::Bool =
    hasproperty(prob.f, name) && _has_sparse_jac_prototype(getproperty(prob.f, name))

@inline function _split_imex_solver_spec(cfg::SolverConfig, sparse_jac::Bool=false)
    mode = cfg.split_imex_solver
    ls = _sparse_linsolve_or_default(sparse_jac)
    mode === :kencarp4  && return (alg=KenCarp4(autodiff=AutoFiniteDiff(), linsolve=ls),  label="KenCarp4")
    mode === :kencarp47 && return (alg=KenCarp47(autodiff=AutoFiniteDiff(), linsolve=ls), label="KenCarp47")
    mode === :kencarp58 && return (alg=KenCarp58(autodiff=AutoFiniteDiff(), linsolve=ls), label="KenCarp58")
    throw(ArgumentError(
        "Unsupported SolverConfig.split_imex_solver=$(repr(mode)). Use one of: :kencarp4, :kencarp47, :kencarp58."
    ))
end
@inline _split_imex_solver_spec() = _split_imex_solver_spec(_active_solver_config())

@inline _multirate_fast_substeps(cfg::SolverConfig)::Int = cfg.multirate_fast_substeps
@inline _multirate_fast_substeps()::Int = _multirate_fast_substeps(_active_solver_config())

@inline function _multirate_slow_dt_s(cfg::SolverConfig, args)::Float64
    default_dt = min(args.integration_tolerances.dt_max_orbit, 2.0)
    dt = isnothing(cfg.multirate_slow_dt_s) ? default_dt : cfg.multirate_slow_dt_s
    dt > 0.0 || throw(ArgumentError("SolverConfig.multirate_slow_dt_s must be > 0.0, got $dt."))
    return min(dt, args.integration_tolerances.dt_max_orbit)
end
@inline _multirate_slow_dt_s(args)::Float64 = _multirate_slow_dt_s(_active_solver_config(), args)

@inline function _multirate_solver_spec_from_sym(mode::Symbol, field_name::String, sparse_jac::Bool=false)
    ls = _sparse_linsolve_or_default(sparse_jac)
    mode === :tsit5     && return (alg=Tsit5(), label="Tsit5", auto_switch_capable=false)
    mode === :auto_stiff && return (
        alg=AutoTsit5(Rodas5P(autodiff=AutoFiniteDiff(), linsolve=ls)),
        label="AutoTsit5(Rodas5P)",
        auto_switch_capable=true
    )
    mode === :rodas5p   && return (alg=Rodas5P(autodiff=AutoFiniteDiff(), linsolve=ls), label="Rodas5P", auto_switch_capable=false)
    mode === :kencarp4  && return (alg=KenCarp4(autodiff=AutoFiniteDiff(), linsolve=ls), label="KenCarp4", auto_switch_capable=false)
    mode === :dp8       && return (alg=DP8(), label="DP8", auto_switch_capable=false)
    throw(ArgumentError(
        "Unsupported SolverConfig.$field_name=$(repr(mode)). Use one of: :tsit5, :dp8, :auto_stiff, :rodas5p, :kencarp4."
    ))
end

@inline _multirate_slow_solver_spec(cfg::SolverConfig, sparse_jac::Bool=false) =
    _multirate_solver_spec_from_sym(cfg.multirate_slow_solver, "multirate_slow_solver", sparse_jac)
@inline _multirate_fast_solver_spec(cfg::SolverConfig, sparse_jac::Bool=false) =
    _multirate_solver_spec_from_sym(cfg.multirate_fast_solver, "multirate_fast_solver", sparse_jac)
@inline _multirate_slow_solver_spec() = _multirate_slow_solver_spec(_active_solver_config())
@inline _multirate_fast_solver_spec() = _multirate_fast_solver_spec(_active_solver_config())

mutable struct SolverIntegratorCache
    integrator::Any
    # Save options the cached integrator was init'ed with.  DiffEq bakes these
    # into the integrator, so reuse is only valid when the requesting call
    # resolves the same options; otherwise the cache is re-initialized (a cache
    # first used on a no-output run must not serve an endpoints-only integrator
    # to a later return_solution=true run).
    save_everystep::Bool
    save_on::Bool
    save_start::Bool
    save_end::Bool
    # dtmax is baked into the integrator at init and `reinit!` does not reset
    # it, so a cache init'ed at one dtmax must not serve a call asking for
    # another -- it would silently integrate with the wrong maximum step. The
    # multirate driver relies on this: its slow half asks for the macro-step
    # width, which shrinks on the final partial macro-step.
    dtmax::Float64
    SolverIntegratorCache() = new(nothing, true, true, true, true, NaN)
end

@inline function _solver_cache_options_match(
    solver_cache::SolverIntegratorCache,
    save_everystep::Bool,
    save_on::Bool,
    save_start::Bool,
    save_end::Bool,
    dtmax::Float64,
)::Bool
    return solver_cache.save_everystep == save_everystep &&
           solver_cache.save_on == save_on &&
           solver_cache.save_start == save_start &&
           solver_cache.save_end == save_end &&
           solver_cache.dtmax == dtmax
end

@inline function _cache_integrator!(
    solver_cache::SolverIntegratorCache,
    integ,
    save_everystep::Bool,
    save_on::Bool,
    save_start::Bool,
    save_end::Bool,
    dtmax::Float64,
)
    solver_cache.integrator = integ
    solver_cache.save_everystep = save_everystep
    solver_cache.save_on = save_on
    solver_cache.save_start = save_start
    solver_cache.save_end = save_end
    solver_cache.dtmax = dtmax
    return integ
end

# Solver save knobs go through the _engine_env_get adapter (overrides → ENV →
# default), keeping ENV access out of this file per the architecture contract.
@inline function _solver_save_everystep()::Bool
    raw = lowercase(strip(_engine_env_get("SPACEAGORA_SOLVER_SAVE_EVERYSTEP", "true")))
    return raw in ("1", "true", "yes", "on")
end

@inline function _solver_bool_env(name::String, default::Bool)::Bool
    raw = lowercase(strip(_engine_env_get(name, default ? "true" : "false")))
    return raw in ("1", "true", "yes", "on")
end

@inline function _solve_with_explicit_solver(prob, cfg::SolverConfig, args, alg, reltol_tol, abstol_tol;
    dtmax_override::Union{Nothing, Float64}=nothing,
    solver_cache::Union{Nothing, SolverIntegratorCache}=nothing,
    needs_full_solution::Bool=true)
    maxiters = _solver_maxiters(cfg)
    dtmax_use = isnothing(dtmax_override) ? args.integration_tolerances.dt_max_orbit : dtmax_override
    dtmax_use > 0.0 || throw(ArgumentError("Solver dtmax must be > 0.0, got $dtmax_use."))
    # When nothing reads the trajectory (return_solution=false, results=false, no
    # solver metadata), skip per-step solution/dense storage — it is the dominant
    # solver-side allocation in campaign runs.  save_on must be gated too: the
    # per-step DiscreteCallbacks save before/after states via save_positions
    # regardless of save_everystep.  save_start/save_end stay on so endpoints
    # and retcode remain available.  Explicitly set SPACEAGORA_SOLVER_SAVE_*
    # env vars still win (documented env semantics); DiffEq derives dense
    # output from save_everystep, so no dense kwarg is needed.
    save_everystep = _engine_env_haskey_with_env_fallback("SPACEAGORA_SOLVER_SAVE_EVERYSTEP") ?
        _solver_save_everystep() : needs_full_solution
    save_on = _engine_env_haskey_with_env_fallback("SPACEAGORA_SOLVER_SAVE_ON") ?
        _solver_bool_env("SPACEAGORA_SOLVER_SAVE_ON", true) : needs_full_solution
    save_start = _solver_bool_env("SPACEAGORA_SOLVER_SAVE_START", true)
    save_end = _solver_bool_env("SPACEAGORA_SOLVER_SAVE_END", true)

    # Reuse the cached integrator only when it was init'ed with the same save
    # options this call resolved; otherwise fall through and re-init the cache.
    if solver_cache !== nothing && solver_cache.integrator !== nothing &&
       _solver_cache_options_match(solver_cache, save_everystep, save_on, save_start, save_end, dtmax_use)
        integ = solver_cache.integrator
        integ.p = prob.p
        SciMLBase.reinit!(integ, prob.u0;
            t0=Float64(first(prob.tspan)),
            tf=Float64(last(prob.tspan)),
            erase_sol=true,
            reinit_callbacks=false)
        return DiffEqBase.solve!(integ)
    end

    if maxiters === nothing
        if solver_cache !== nothing
            integ = DiffEqBase.init(prob, alg; reltol=reltol_tol, abstol=abstol_tol, dtmax=dtmax_use, save_everystep=save_everystep, save_on=save_on, save_start=save_start, save_end=save_end)
            _cache_integrator!(solver_cache, integ, save_everystep, save_on, save_start, save_end, dtmax_use)
            return DiffEqBase.solve!(integ)
        end
        return solve(prob, alg; reltol=reltol_tol, abstol=abstol_tol, dtmax=dtmax_use, save_everystep=save_everystep, save_on=save_on, save_start=save_start, save_end=save_end)
    end
    if solver_cache !== nothing
        integ = DiffEqBase.init(prob, alg; reltol=reltol_tol, abstol=abstol_tol, dtmax=dtmax_use, maxiters=maxiters, save_everystep=save_everystep, save_on=save_on, save_start=save_start, save_end=save_end)
        _cache_integrator!(solver_cache, integ, save_everystep, save_on, save_start, save_end, dtmax_use)
        return DiffEqBase.solve!(integ)
    end
    return solve(prob, alg; reltol=reltol_tol, abstol=abstol_tol, dtmax=dtmax_use, maxiters=maxiters, save_everystep=save_everystep, save_on=save_on, save_start=save_start, save_end=save_end)
end

@inline function _solve_with_explicit_solver(prob, args, alg, reltol_tol, abstol_tol;
    dtmax_override::Union{Nothing, Float64}=nothing,
    solver_cache::Union{Nothing, SolverIntegratorCache}=nothing,
    needs_full_solution::Bool=true)
    return _solve_with_explicit_solver(
        prob,
        _active_solver_config(),
        args,
        alg,
        reltol_tol,
        abstol_tol;
        dtmax_override=dtmax_override,
        solver_cache=solver_cache,
        needs_full_solution=needs_full_solution,
    )
end

@inline function _solve_with_fixed_step_solver(prob, cfg::SolverConfig, alg, dt_s::Float64;
    needs_full_solution::Bool=true)
    maxiters = _solver_maxiters(cfg)
    save_everystep = _engine_env_haskey_with_env_fallback("SPACEAGORA_SOLVER_SAVE_EVERYSTEP") ?
        _solver_save_everystep() : needs_full_solution
    save_on = _engine_env_haskey_with_env_fallback("SPACEAGORA_SOLVER_SAVE_ON") ?
        _solver_bool_env("SPACEAGORA_SOLVER_SAVE_ON", true) : needs_full_solution
    save_start = _solver_bool_env("SPACEAGORA_SOLVER_SAVE_START", true)
    save_end = _solver_bool_env("SPACEAGORA_SOLVER_SAVE_END", true)
    if maxiters === nothing
        return solve(prob, alg; dt=dt_s, save_everystep=save_everystep, save_on=save_on, save_start=save_start, save_end=save_end)
    end
    return solve(prob, alg; dt=dt_s, maxiters=maxiters, save_everystep=save_everystep, save_on=save_on, save_start=save_start, save_end=save_end)
end
@inline _solve_with_fixed_step_solver(prob, alg, dt_s::Float64) = _solve_with_fixed_step_solver(prob, _active_solver_config(), alg, dt_s)

@inline function _is_partitioned_second_order_problem(prob)::Bool
    hasproperty(prob, :problem_type) || return false
    return getproperty(prob, :problem_type) isa SecondOrderODEProblem
end

@inline function _split_subproblem(prob, f, u, tspan)
    return ODEProblem(f, u, tspan, prob.p; prob.kwargs...)
end

function _solve_with_multirate_solver(prob, cfg::SolverConfig, args, reltol_tol, abstol_tol)
    if !(hasproperty(prob.f, :f1) && hasproperty(prob.f, :f2))
        throw(ArgumentError("SolverConfig.solver_mode=:multirate requires a split problem with f1/f2 components."))
    end

    t_start = Float64(first(prob.tspan))
    t_end = Float64(last(prob.tspan))
    if t_end <= t_start
        sol = _solve_with_explicit_solver(prob, cfg, args, Tsit5(), reltol_tol, abstol_tol)
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

    # Each half is handed to _split_subproblem as its own ODEProblem, so the
    # sparse linear solve is opt-in per component: whichever of f1/f2 carries a
    # prototype gets KLU, the other keeps the dense default.
    slow_spec = _multirate_slow_solver_spec(cfg, _split_component_has_sparse_jac(prob, :f1))
    fast_spec = _multirate_fast_solver_spec(cfg, _split_component_has_sparse_jac(prob, :f2))
    fast_substeps = _multirate_fast_substeps(cfg)
    slow_dt_s = _multirate_slow_dt_s(cfg, args)
    fast_dt_s = slow_dt_s / fast_substeps

    t_cursor = t_start
    # Carried across macro-steps and copied into, never rebound: the sub-solve
    # results can alias a cached integrator's own state buffer, which the next
    # reinit! overwrites, so the handoff state must live in a buffer the
    # integrators do not own.
    u_cursor = deepcopy(prob.u0)
    final_sol = nothing
    macro_steps = 0
    auto_switch_events = 0

    # One cached integrator per split component. They cannot share one: `reinit!`
    # rebinds the state and time span but never the RHS, and f1/f2 are different
    # functions. Both fast half-steps resolve the same dtmax, so a single fast
    # cache serves them; the slow cache re-inits itself on the final partial
    # macro-step, where the requested dtmax shrinks.
    #
    # Without this the driver built a complete integrator -- cache, jac config,
    # W, and KLU symbolic factorization -- three times per macro-step and threw
    # it away. Measured on a 256-satellite-shaped block-diagonal problem
    # (n=2560): init 0.834 ms / 4.3 MB against reinit! 0.005 ms / 0.1 MB, and a
    # fresh init costs more than four integration steps while each sub-solve
    # takes only a handful.
    slow_cache = SolverIntegratorCache()
    fast_cache = SolverIntegratorCache()

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
                cfg,
                args,
                fast_spec.alg,
                reltol_tol,
                abstol_tol;
                dtmax_override=min(fast_dt_s, half_dt),
                solver_cache=fast_cache
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
            copyto!(u_cursor, sol_fast_pre.u[end])
            final_sol = sol_fast_pre
        end

        slow_prob = _split_subproblem(prob, prob.f.f1, u_cursor, (t_cursor, t_next))
        sol_slow = _solve_with_explicit_solver(
            slow_prob,
            cfg,
            args,
            slow_spec.alg,
            reltol_tol,
            abstol_tol;
            dtmax_override=segment_dt,
            solver_cache=slow_cache
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
        copyto!(u_cursor, sol_slow.u[end])
        final_sol = sol_slow

        if half_dt > 0.0
            fast_prob_post = _split_subproblem(prob, prob.f.f2, u_cursor, (t_half, t_next))
            sol_fast_post = _solve_with_explicit_solver(
                fast_prob_post,
                cfg,
                args,
                fast_spec.alg,
                reltol_tol,
                abstol_tol;
                dtmax_override=min(fast_dt_s, half_dt),
                solver_cache=fast_cache
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
            copyto!(u_cursor, sol_fast_post.u[end])
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
_solve_with_multirate_solver(prob, args, reltol_tol, abstol_tol) =
    _solve_with_multirate_solver(prob, _active_solver_config(), args, reltol_tol, abstol_tol)

function _solve_with_gravity_backbone_solver(prob, cfg::SolverConfig, args)
    t_start = Float64(first(prob.tspan))
    t_end = Float64(last(prob.tspan))
    alg = KahanLi8()
    solver_label = _gravity_backbone_has_kicks(args) ? "KahanLi8(GravityBackbone+Kicks)" : "KahanLi8(GravityBackbone)"

    if t_end <= t_start
        sol = DiffEqBase.build_solution(
            prob,
            alg,
            Float64[t_start],
            [deepcopy(prob.u0)];
            retcode=SciMLBase.ReturnCode.Success,
        )
        return sol, solver_label
    end

    dt_s = _gravity_backbone_fixed_dt_s(cfg, args)
    t_cursor = t_start
    u_cursor = deepcopy(prob.u0)
    solution_ts = Float64[t_cursor]
    solution_us = [deepcopy(u_cursor)]
    final_retcode = SciMLBase.ReturnCode.Success

    while t_cursor < t_end
        t_next = min(t_cursor + dt_s, t_end)
        segment_dt = t_next - t_cursor
        if !(isfinite(segment_dt) && segment_dt > 0.0)
            break
        end
        half_dt = 0.5 * segment_dt

        if half_dt > 0.0
            _gravity_backbone_half_kick!(u_cursor, prob.p, t_cursor, half_dt)
        end

        core_prob = remake(prob; u0=u_cursor, tspan=(t_cursor, t_next))
        core_sol = _solve_with_fixed_step_solver(core_prob, cfg, alg, segment_dt)
        reached_t = Float64(core_sol.t[end])
        u_cursor = deepcopy(core_sol.u[end])
        final_retcode = core_sol.retcode

        if !SciMLBase.successful_retcode(core_sol.retcode)
            if reached_t > solution_ts[end]
                push!(solution_ts, reached_t)
                push!(solution_us, deepcopy(u_cursor))
            end
            break
        end

        if string(core_sol.retcode) != "Success"
            if reached_t > solution_ts[end]
                push!(solution_ts, reached_t)
                push!(solution_us, deepcopy(u_cursor))
            end
            break
        end

        if half_dt > 0.0
            _gravity_backbone_half_kick!(u_cursor, prob.p, reached_t, half_dt)
        end

        t_cursor = reached_t
        if t_cursor > solution_ts[end]
            push!(solution_ts, t_cursor)
            push!(solution_us, deepcopy(u_cursor))
        else
            solution_us[end] = deepcopy(u_cursor)
        end

        _gravity_backbone_time_reached(t_cursor, t_next) || break
    end

    sol = DiffEqBase.build_solution(
        prob,
        alg,
        solution_ts,
        solution_us;
        retcode=final_retcode,
    )
    return sol, solver_label
end

function _solve_with_solver_policy(prob, cfg::SolverConfig, args, reltol_tol, abstol_tol;
    solver_cache::Union{Nothing, SolverIntegratorCache}=nothing,
    needs_full_solution::Bool=true)
    mode = _solver_policy_mode(cfg)
    if mode == :symplectic
        _symplectic_conservative_eligible(args) || throw(ArgumentError(
            "SolverConfig.solver_mode=:symplectic currently requires a single-spacecraft, translational-only inverse-squared gravity configuration with no control/guidance/navigation effectors."
        ))
        _is_partitioned_second_order_problem(prob) || throw(ArgumentError(
            "SolverConfig.solver_mode=:symplectic requires a partitioned SecondOrderODEProblem. The typed run_simulation path still builds a first-order ODEProblem, so use :tsit5/:auto_stiff there until a partitioned runtime path is added."
        ))
        sol = _solve_with_fixed_step_solver(prob, cfg, KahanLi8(), _symplectic_fixed_dt_s(cfg, args); needs_full_solution=needs_full_solution)
        return sol, (
            solver="KahanLi8(Symplectic)",
            initial_solver="KahanLi8",
            fallback_used=false,
            trigger_retcode=missing
        )
    end

    if mode == :gravity_backbone_split
        reject_reason = _gravity_backbone_reject_reason(args)
        isnothing(reject_reason) || throw(ArgumentError(reject_reason))
        _is_partitioned_second_order_problem(prob) || throw(ArgumentError(
            "SolverConfig.solver_mode=:gravity_backbone_split requires a SecondOrderODEProblem built over translational position/velocity states."
        ))
        sol, solver_label = _solve_with_gravity_backbone_solver(prob, cfg, args)
        return sol, (
            solver=solver_label,
            initial_solver="KahanLi8",
            fallback_used=false,
            trigger_retcode=missing
        )
    end

    # KLUFactorization requires a sparse W matrix. Only use it when a sparse
    # jac_prototype was provided (multi-satellite path); fall back to dense LU
    # otherwise so that single-satellite runs don't get wrong Newton corrections
    # from KLU on a dense matrix. The check asks for a sparse prototype rather
    # than a non-nothing one so a dense prototype cannot route here either.
    _rodas5p_alg() = Rodas5P(
        autodiff=AutoFiniteDiff(),
        linsolve=_sparse_linsolve_or_default(_has_sparse_jac_prototype(prob.f)),
    )

    if mode == :rodas5p
        sol = _solve_with_explicit_solver(prob, cfg, args, _rodas5p_alg(), reltol_tol, abstol_tol; solver_cache=solver_cache, needs_full_solution=needs_full_solution)
        return sol, (
            solver="Rodas5P",
            initial_solver="Rodas5P",
            fallback_used=false,
            trigger_retcode=missing
        )
    end

    if mode == :auto_stiff
        if _auto_stiff_smooth_gravity_eligible(cfg, args)
            sol = _solve_with_explicit_solver(prob, cfg, args, Tsit5(), reltol_tol, abstol_tol; solver_cache=solver_cache, needs_full_solution=needs_full_solution)
            return sol, (
                solver="Tsit5",
                initial_solver="Tsit5",
                fallback_used=false,
                trigger_retcode=missing
            )
        end

        # True stiffness-aware autoswitching handled internally by OrdinaryDiffEq.
        # This replaces the manual "retry with Rodas5P on Tsit5 failure" policy.
        autoswitch_alg = AutoTsit5(_rodas5p_alg(); switch_max=_auto_stiff_switch_max(cfg))
        # Keep per-step storage here regardless of needs_full_solution:
        # _auto_stiff_switched inspects sol.alg_choice across saved steps.
        sol = _solve_with_explicit_solver(prob, cfg, args, autoswitch_alg, reltol_tol, abstol_tol; solver_cache=solver_cache)
        switched = _auto_stiff_switched(sol)
        return sol, (
            solver="AutoTsit5(Rodas5P)",
            initial_solver="AutoTsit5",
            fallback_used=switched,
            trigger_retcode=switched ? "internal_autoswitch" : missing
        )
    end

    if mode == :split_imex
        split_solver = _split_imex_solver_spec(cfg, _split_component_has_sparse_jac(prob, :f1))
        sol = _solve_with_explicit_solver(prob, cfg, args, split_solver.alg, reltol_tol, abstol_tol; solver_cache=solver_cache, needs_full_solution=needs_full_solution)
        return sol, (
            solver="$(split_solver.label)(IMEX)",
            initial_solver=split_solver.label,
            fallback_used=false,
            trigger_retcode=missing
        )
    end

    if mode == :multirate
        sol, multirate_meta = _solve_with_multirate_solver(prob, cfg, args, reltol_tol, abstol_tol)
        switched = multirate_meta.auto_switch_events > 0
        return sol, (
            solver="Multirate(Strang; slow=$(multirate_meta.slow_solver), fast=$(multirate_meta.fast_solver))",
            initial_solver=multirate_meta.slow_solver,
            fallback_used=switched,
            trigger_retcode=switched ? "internal_autoswitch" : missing
        )
    end

    if mode == :dp8
        sol = _solve_with_explicit_solver(prob, cfg, args, DP8(), reltol_tol, abstol_tol; solver_cache=solver_cache, needs_full_solution=needs_full_solution)
        return sol, (
            solver="DP8",
            initial_solver="DP8",
            fallback_used=false,
            trigger_retcode=missing
        )
    end

    tsit_sol = _solve_with_explicit_solver(prob, cfg, args, Tsit5(), reltol_tol, abstol_tol; solver_cache=solver_cache, needs_full_solution=needs_full_solution)
    return tsit_sol, (
        solver="Tsit5",
        initial_solver="Tsit5",
        fallback_used=false,
        trigger_retcode=missing
    )
end

function _solve_with_solver_policy(prob, args, reltol_tol, abstol_tol;
    solver_cache::Union{Nothing, SolverIntegratorCache}=nothing,
    needs_full_solution::Bool=true)
    return _solve_with_solver_policy(
        prob,
        _active_solver_config(),
        args,
        reltol_tol,
        abstol_tol;
        solver_cache=solver_cache,
        needs_full_solution=needs_full_solution,
    )
end
