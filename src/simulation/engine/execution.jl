"""
Block-diagonal Jacobian sparsity for a constellation state, one block per satellite.

`is_active` (when given) drops an inactive satellite's block to just its
diagonal. The RHS already zeroes those derivatives, so their off-diagonal
entries are structurally zero, and carrying them costs nnz in every W and one
colour in every finite-difference Jacobian for the rest of the run. This is only
sound because deactivation is permanent -- `p.is_active[idx] = false` in
event_callbacks.jl is the sole write and nothing sets it back -- so a pattern
built from a later snapshot can only be denser than the truth, never sparser.

Satellite blocks are contiguous index ranges in the ComponentArray layout, which
lets the CSC arrays be filled directly: `colptr` and `rowval` follow from each
column's owning block, with no COO triple and no `sparse()` sort. The general
path below still handles a layout where they are not contiguous.
"""
function _build_block_diagonal_jac_prototype(
    u0::ComponentVector,
    is_active::Union{Nothing, AbstractVector{Bool}}=nothing,
)::SparseMatrixCSC{Float64, Int}
    n_sats = length(u0.sc)
    n_total = length(u0)
    # Probe every satellite block in ONE pass: stamp each satellite's slots with
    # its own index, flatten once, then read the stamps back. This makes no
    # assumption about uniform block size (heterogeneous n_bodies / heat_loads
    # per satellite are recovered from the stamps), and costs one O(n_total)
    # copy plus one O(n_total) scan instead of the per-satellite `zero(u0)` +
    # full-length copy + full-length `findall` it used to, which made
    # construction O(n_sats * n_total) -- quadratic in constellation size.
    probe = zero(u0)
    @inbounds for i in 1:n_sats
        probe.sc[i] .= Float64(i)
    end
    flat = Vector{Float64}(probe)

    owner = zeros(Int, n_total)
    first_idx = fill(typemax(Int), n_sats)
    last_idx = zeros(Int, n_sats)
    counts = zeros(Int, n_sats)
    @inbounds for k in eachindex(flat)
        stamp = flat[k]
        stamp == 0.0 && continue
        s = Int(stamp)
        owner[k] = s
        k < first_idx[s] && (first_idx[s] = k)
        k > last_idx[s] && (last_idx[s] = k)
        counts[s] += 1
    end

    @inline _active(s::Int) = is_active === nothing || s > length(is_active) || is_active[s]

    contiguous = true
    @inbounds for s in 1:n_sats
        counts[s] == 0 && continue
        if last_idx[s] - first_idx[s] + 1 != counts[s]
            contiguous = false
            break
        end
    end

    if contiguous
        nnz_total = 0
        @inbounds for s in 1:n_sats
            counts[s] == 0 && continue
            nnz_total += _active(s) ? counts[s] * counts[s] : counts[s]
        end
        # An index owned by no satellite still needs its diagonal, or the column
        # is empty and W is structurally singular.
        @inbounds for k in 1:n_total
            owner[k] == 0 && (nnz_total += 1)
        end

        colptr = Vector{Int}(undef, n_total + 1)
        rowval = Vector{Int}(undef, nnz_total)
        pos = 1
        @inbounds for j in 1:n_total
            colptr[j] = pos
            s = owner[j]
            if s == 0 || !_active(s)
                rowval[pos] = j
                pos += 1
            else
                for i in first_idx[s]:last_idx[s]
                    rowval[pos] = i
                    pos += 1
                end
            end
        end
        colptr[n_total + 1] = pos
        return SparseMatrixCSC(n_total, n_total, colptr, rowval, ones(Float64, nnz_total))
    end

    # Non-contiguous layout: fall back to assembling from index lists.
    blocks = [Int[] for _ in 1:n_sats]
    @inbounds for k in eachindex(flat)
        owner[k] == 0 && continue
        push!(blocks[owner[k]], k)
    end
    nnz_total = 0
    @inbounds for (s, b) in enumerate(blocks)
        nnz_total += _active(s) ? length(b)^2 : length(b)
    end
    rows = Vector{Int}(undef, nnz_total)
    cols = Vector{Int}(undef, nnz_total)
    pos = 0
    @inbounds for (s, b) in enumerate(blocks)
        if _active(s)
            for j in b, k in b
                pos += 1
                rows[pos] = k
                cols[pos] = j
            end
        else
            for j in b
                pos += 1
                rows[pos] = j
                cols[pos] = j
            end
        end
    end
    return sparse(rows, cols, ones(Float64, nnz_total), n_total, n_total)
end

# Satellite-index fields an effector uses to name a satellite other than the one
# it is being evaluated for. An effector carrying two of these with different
# values reads across satellites, so its contribution to the Jacobian lands
# off-block.
const _CROSS_SATELLITE_INDEX_FIELDS = (:chaser_idx, :target_idx)

"""Return whether `effector` reads state belonging to more than one satellite."""
@inline function _effector_couples_satellites(effector)::Bool
    first_idx = nothing
    @inbounds for name in _CROSS_SATELLITE_INDEX_FIELDS
        hasproperty(effector, name) || continue
        idx = getproperty(effector, name)
        idx isa Integer || continue
        first_idx === nothing && (first_idx = idx; continue)
        first_idx == idx || return true
    end
    return false
end

"""Return whether any configured effector couples one satellite's dynamics to another's."""
function _config_has_cross_satellite_coupling(args)::Bool
    for (model_name, effectors_name) in (
        (:dynamics_model, :dynamic_effectors),
        (:guidance_model, :guidance_effectors),
        (:navigation_model, :navigation_effectors),
        (:control_model, :control_effectors),
    )
        hasproperty(args, model_name) || continue
        model = getproperty(args, model_name)
        hasproperty(model, effectors_name) || continue
        @inbounds for effector in getproperty(model, effectors_name)
            _effector_couples_satellites(effector) && return true
        end
    end
    return false
end

# Attach the block-diagonal sparsity to one component of a split problem.
#
# The prototype has to sit on the inner ODEFunction wrapping the component, NOT
# on the enclosing SplitFunction. Passing `jac_prototype` to SplitFunction alone
# builds a sparse W -- so KLUFactorization accepts it and nothing errors -- but
# the finite-difference Jacobian is still taken dense, column by column, because
# the sparsity that OrdinaryDiffEq colors is read off the component function it
# actually differentiates. Measured on a 100-block/600-state IMEX repro:
# SplitFunction-level prototype left the implicit-part evaluation count at 4351,
# identical to no prototype at all; moving it onto the inner ODEFunction dropped
# it to 228. The silent no-op is the whole reason this helper exists rather than
# a `jac_prototype=` keyword at the SplitODEProblem call sites.
@inline function _split_component_function(f, jac_prototype::Union{Nothing, SparseMatrixCSC{Float64, Int}})
    jac_prototype === nothing && return f
    return ODEFunction(f; jac_prototype=jac_prototype)
end

@inline function _build_typed_solver_problem(u0, tspan, p, callbacks, solver_mode::Symbol,
    jac_prototype::Union{Nothing, SparseMatrixCSC{Float64, Int}}=nothing)
    mode = solver_mode
    if mode == :gravity_backbone_split
        q0, dq0 = _gravity_backbone_initial_states(u0, p.args)
        return SecondOrderODEProblem(
            spacecraft_dynamics_gravity_backbone!,
            dq0,
            q0,
            tspan,
            p,
            callback=callbacks
        )
    elseif mode == :split_imex
        return SplitODEProblem(
            _split_component_function(spacecraft_dynamics_implicit_atmosphere!, jac_prototype),
            spacecraft_dynamics_explicit_remainder!,
            u0,
            tspan,
            p,
            callback=callbacks
        )
    elseif mode == :multirate
        return SplitODEProblem(
            _split_component_function(spacecraft_dynamics_slow!, jac_prototype),
            _split_component_function(spacecraft_dynamics_fast_control!, jac_prototype),
            u0,
            tspan,
            p,
            callback=callbacks
        )
    end
    if jac_prototype !== nothing
        f = ODEFunction(spacecraft_dynamics!; jac_prototype=jac_prototype)
        return ODEProblem(f, u0, tspan, p; callback=callbacks)
    end
    return ODEProblem(spacecraft_dynamics!, u0, tspan, p, callback=callbacks)
end

@inline function _build_typed_solver_problem(u0, tspan, p, callbacks,
    jac_prototype::Union{Nothing, SparseMatrixCSC{Float64, Int}}=nothing)
    return _build_typed_solver_problem(u0, tspan, p, callbacks, _solver_policy_mode(_active_solver_config()), jac_prototype)
end

@inline function _append_backbone_saved_segment!(
    times_acc::Vector{Float64},
    data_acc::Vector{SimulationModel.SaveData},
    sol,
    save_fields,
    p,
)
    length(sol.t) <= 1 && return nothing
    integrator_view = (p=p,)
    @inbounds for idx in 2:length(sol.t)
        t_sample = Float64(sol.t[idx])
        push!(times_acc, t_sample)
        push!(data_acc, SimulationModel.SimulationCallbacks._save_snapshot(save_fields, sol.u[idx], t_sample, integrator_view))
    end
    return nothing
end

@inline function _append_checkpoint_saved_segment!(
    times_acc::Vector{Float64},
    data_acc::Vector,
    saved_values,
    sol,
    save_fields,
    p,
)
    _append_saved_segment!(times_acc, data_acc, saved_values)
    isempty(sol.t) && return nothing

    t_final = Float64(sol.t[end])
    if isempty(times_acc) || !isapprox(times_acc[end], t_final; atol=0.0, rtol=0.0)
        integrator_view = (p=p,)
        push!(times_acc, t_final)
        push!(data_acc, SimulationModel.SimulationCallbacks._save_snapshot(save_fields, sol.u[end], t_final, integrator_view))
    end
    return nothing
end

function _save_simulation_results_if_enabled!(
    args,
    solver_mode::Symbol,
    checkpoint_active::Bool,
    save_fields_resolved,
    saved_values,
    checkpoint_saved_times,
    checkpoint_saved_data,
    backbone_saved_times,
    backbone_saved_data,
)
    args.simulation_settings.results || return nothing
    results_times = if solver_mode == :gravity_backbone_split
        checkpoint_active ? checkpoint_saved_times : backbone_saved_times
    else
        checkpoint_active ? checkpoint_saved_times : saved_values.t
    end
    results_data = if solver_mode == :gravity_backbone_split
        checkpoint_active ? checkpoint_saved_data : backbone_saved_data
    else
        checkpoint_active ? checkpoint_saved_data : saved_values.saveval
    end
    results_df = _build_results_dataframe(results_times, results_data, save_fields_resolved, args)
    csv_path = _write_results_csv!(results_df, args)
    if _typed_save_bundle_enabled()
        _write_results_bundle!(results_df, results_times, args; csv_path=csv_path)
    end
    return csv_path
end

function _try_save_simulation_results_if_enabled!(args...)
    try
        return _save_simulation_results_if_enabled!(args...)
    catch err
        @warn "Failed to save partial simulation results after solve failure." reason=sprint(showerror, err)
        return nothing
    end
end

function run_simulation(
    args::SimulationConfiguration;
    isolate_state::Bool=true,
    return_solution::Bool=false,
    return_solver_metadata::Bool=false,
    save_fields=nothing,
    extra_callbacks=(),
    solver_cache::Union{Nothing, SolverIntegratorCache}=nothing
)
    return SimulationModel.ParallelPolicy.with_policy_context() do
    # Isolate mutable campaign/model state by default so repeated/concurrent runs
    # do not alias shared in-memory objects.
    args = isolate_state ? deepcopy(args) : args
    # Resolve effective solver config: explicit field on args wins; otherwise read from env
    # (which respects any active SimulationEngineConfig overrides via _engine_env_get).
    solver_cfg = isnothing(args.solver_config) ? _solver_config_from_env() : args.solver_config
    solver_mode = solver_cfg.solver_mode

    # Typed pipeline is SI-native (meters, seconds, kilograms). The
    # `simulation_settings.normalize` field is legacy-only and rejected by default.
    _enforce_typed_normalize_policy!(args)
    _validate_orientation_inertia!(args)
    _validate_thermal_model_support!(args)
    _validate_ephemerides_support!(args)
    _warn_density_without_atmospheric_effector(args)
    try
        SimulationModel.ParallelPolicy.reset_policy_telemetry!()
        if SimulationModel.ParallelPolicy.persistent_hints_state_reset_requested()
            SimulationModel.ParallelPolicy.reset_persistent_hint_state!()
        end
    catch
    end

    # Set up the model and initial conditions
    initial_conditions = build_initial_conditions(args)
    if args.simulation_settings.verbose
        println("Initial conditions:")
        println(initial_conditions)
    end

    # Define the ODE parameters and callbacks
    p = SimulationModel.ODEParams(n_sats=length(args.dynamics_model.spacecraft), args=args) # Define the parameters for the ODE problem, including the shared buffers for the callbacks
    _initialize_in_atmosphere_flags!(p, initial_conditions)
    # Snapshot env-derived runtime config (policy/RHS-plan/callback knobs) once,
    # inside the active SimulationEngineConfig override scope, so hot paths read
    # typed struct fields instead of re-parsing process-global ENV per RHS call.
    _initialize_runtime_env_config!(p)
    _initialize_heat_rate_buffers!(p)
    _initialize_save_cache_buffers!(p)
    _initialize_density_model_instances!(p)
    _initialize_density_cache_buffers!(p)
    _initialize_gram_isolated_pool_buffers!(p)
    _initialize_harmonics_workspace_buffers!(p)
    _initialize_nbody_workspace_buffers!(p)
    _initialize_aero_workspace_buffers!(p)
    _initialize_nbody_ephemeris_cache_buffer!(p)
    _initialize_srp_sun_cache_buffer!(p)
    _initialize_planet_frame_cache_buffer!(p)
    _initialize_spice_rhs_memo_mode!(p)
    _reset_spice_runtime_counters!(p)
    _reset_spice_rhs_memo!(p)
    p.shared_buffers.debug_control[] = _engine_env_get("SPACEAGORA_DEBUG_CONTROL", "0") == "1"
    p.shared_buffers.debug_initial_derivative[] = _engine_env_get("SPACEAGORA_DEBUG_INITIAL_DERIVATIVE", "0") == "1"
    save_fields_resolved = isnothing(save_fields) ? SimulationModel.default_save_fields(args) : collect(save_fields)
    save_field_names = Symbol[field.name for field in save_fields_resolved]
    length(unique(save_field_names)) == length(save_field_names) || throw(ArgumentError("save_fields names must be unique. Got $(save_field_names)."))
    saved_values = SavedValues(Float64, SimulationModel.SaveData)
    callbacks = SimulationModel.get_callbacks(
        length(args.dynamics_model.spacecraft),
        args.dynamics_model.dynamic_effectors,
        args;
        saved_values=saved_values,
        save_fields=save_fields_resolved,
        extra_callbacks=extra_callbacks
    ) # Get the callbacks based on the number of satellites and the dynamic effectors being used in the simulation
    ephemerides_model = args.environment_model.ephemerides_model
    et_start = SimulationModel.ephemerides_time_seconds(args.initial_time, ephemerides_model)
    p.shared_buffers.et_start[] = et_start
    args.environment_model.planet.L_PI .= SimulationModel.planet_frame_lpi(args.environment_model.planet, et_start, ephemerides_model)
    if SimulationModel.ephemerides_requires_spice(ephemerides_model)
        Base.Threads.atomic_add!(p.shared_buffers.spice_runtime_counters.planet_pxform_runtime_calls, 1)
    end
    mission_end = args.mission_configuration.mission_time
    _initialize_nbody_ephemeris_cache!(p, et_start, mission_end)
    _initialize_srp_sun_ephemeris_cache!(p, et_start, mission_end)
    _initialize_planet_frame_ephemeris_cache!(p, et_start, mission_end)
    checkpoint_active = _typed_checkpoint_enabled(args)
    if checkpoint_active && args.simulation_settings.checkpoint_interval_s <= 0.0
        throw(ArgumentError("SimulationSettings.checkpoint_interval_s must be > 0 when checkpointing is enabled."))
    end

    t_start = 0.0
    u_start = initial_conditions
    if args.simulation_settings.resume_from_checkpoint
        ckpt = _load_checkpoint(args)
        if ckpt === nothing
            @warn "resume_from_checkpoint=true but no checkpoint file was found; starting from initial conditions."
        else
            if solver_mode == :gravity_backbone_split
                ckpt.solver_mode == "gravity_backbone_split" || throw(ArgumentError(
                    "gravity_backbone_split can only resume from checkpoints written by SPACEAGORA_SOLVER_MODE=gravity_backbone_split."
                ))
            elseif ckpt.solver_mode == "gravity_backbone_split"
                throw(ArgumentError(
                    "Checkpoint was written by gravity_backbone_split and cannot be resumed with solver mode $(solver_mode)."
                ))
            end
            t_start = ckpt.t
            u_start = ckpt.u
            if args.simulation_settings.verbose
                println("Resuming simulation from checkpoint at t=$(round(t_start, digits=6)) s")
            end
        end
    end

    # println("Initial conditions:")
    # println(initial_conditions)
    # println("ODE parameters:")
    # println(p)
    # println("args.mission_configuration.mission_time: $(args.mission_configuration.mission_time)")
    p.shared_buffers.solve_segment_end_time[] = mission_end
    prob_debug_state = solver_mode == :gravity_backbone_split ? initial_conditions : u_start
    prob_debug = ODEProblem(spacecraft_dynamics!, prob_debug_state, (t_start, mission_end), p, callback=callbacks)
    if p.shared_buffers.debug_initial_derivative[] && solver_mode != :gravity_backbone_split
        # 1. Manually evaluate the derivative at the start
        du_test = copy(prob_debug.u0)
        try
            prob_debug.f(du_test, prob_debug.u0, prob_debug.p, prob_debug.tspan[1])
        catch e
            @error "The derivative function itself crashed!" exception=e
            throw(ErrorException("Initial derivative evaluation failed; aborting solve in debug mode."))
        end

        # 2. Check for NaNs and print exactly where they are
        if any(isnan, du_test)
            println("--- INITIAL NaN DETECTED ---")

            # Check global parameters in p
            _debug_print_nan_parameter_paths!(prob_debug.p)

            # Check the state vector (u)
            # Assuming your u has a .sc field for satellites
            for (i, sat) in enumerate(du_test.sc)
                if any(isnan, sat.pos) || any(isnan, sat.vel)
                    println("NaN found in Satellite $i derivative!")
                    println("  Pos: $(sat.pos)")
                    println("  Vel: $(sat.vel)")
                end
            end

            throw(ErrorException("Initial derivative contains NaN; aborting solve in debug mode."))
        end
    elseif p.shared_buffers.debug_initial_derivative[] && solver_mode == :gravity_backbone_split
        @warn "Initial derivative debug is not supported for gravity_backbone_split; skipping first-order RHS debug probe."
    end

    reltol_tol, abstol_tol = if solver_mode == :gravity_backbone_split
        args.integration_tolerances.reltol_orbit, args.integration_tolerances.abstol_orbit
    else
        _build_solver_tolerances(u_start, args)
    end
    # Block-diagonal Jacobian sparsity for the implicit solver paths.
    #
    # :gravity_backbone_split is excluded because it integrates with an explicit
    # symplectic method and never forms a Jacobian. :split_imex and :multirate
    # are NOT excluded: both run implicit solvers (KenCarp*, Rodas5P) over the
    # full constellation state, and without a prototype they build a dense n x n
    # W and take the finite-difference Jacobian dense, costing n+1 RHS
    # evaluations per Jacobian and an O(n^3) factorization. Measured on a
    # 200-block/1200-state IMEX repro: 0.696 s and 8551 implicit-part
    # evaluations dense, against 0.002 s and 228 with this prototype.
    #
    # The pattern asserts that satellite i's derivative depends only on
    # satellite i's state. An effector that reads a second satellite (RPO
    # planners and the RPO MPC controller both read the target's pos/vel while
    # computing the chaser's wrench) breaks that assumption, and the entries it
    # needs sit off-block where this pattern has structural zeros. Supplying a
    # too-sparse Jacobian is not merely a Newton convergence penalty for
    # Rodas5P: Rosenbrock methods carry the supplied Jacobian in their order
    # conditions, so a truncated one reduces the order of the method. Fall back
    # to a dense Jacobian in that case -- correctness first; extending the
    # pattern with the coupling blocks would recover the sparsity and is the
    # better fix, but it needs the effectors to report which pairs they couple.
    jac_prototype = (
        solver_mode != :gravity_backbone_split &&
        u_start isa ComponentVector &&
        length(u_start.sc) > 1 &&
        !_config_has_cross_satellite_coupling(args)
    ) ? _build_block_diagonal_jac_prototype(u_start, p.is_active) : nothing
    # Satellites that finish early stop contributing off-diagonal structure.
    # Rebuild between checkpoint segments when the active set has shrunk, so a
    # long run does not keep paying peak Jacobian width after most of the
    # constellation has dropped out. Deactivation is permanent, so the count
    # alone identifies a change.
    jac_prototype_active_count = jac_prototype === nothing ? 0 : count(identity, p.is_active)

    # Auto-calibration: time candidate execution plans before the solve and pin the
    # fastest one for the duration.  No-ops when budget <= 1, SPACEAGORA_RHS_CALIBRATE=off,
    # or a cached result already exists for this machine + scenario signature.
    # The native-lock Amdahl bound for this solve.
    #
    # THE DENOMINATOR IS THE WHOLE DIFFICULTY. What bounds inner width is the
    # fraction of an RHS EVALUATION spent holding the shared native lock,
    # because that is the part a wider split cannot divide. Dividing instead by
    # the solve's wall clock folds in the solver spine, the callbacks and the
    # saving path -- work that has nothing to do with the RHS split -- and
    # understates it by two orders of magnitude. Measured on a live-GRAM solve:
    # hold over RHS time is 0.89, hold over solve wall is 0.008, and only the
    # first is a statement about what widening the RHS can achieve. A first
    # version of this took the wall-clock ratio and the cap never engaged on the
    # one workload it exists for.
    #
    # `probe_contention_inputs!` measures exactly that ratio over warm passes,
    # so the ceiling is derived from it rather than from the lock counters'
    # running totals. A failed or invalid probe applies no constraint: "not
    # measured" and "no constraint" must reach the same place, since a bound
    # invented from a failed measurement is worse than none.
    if _rhs_lock_width_cap_enabled()
        contention = try
            probe_contention_inputs!(
                p, u_start, args;
                k = 5,
                t_step = Float64(args.mission_configuration.mission_time) / 8)
        catch err
            @debug "Lock width cap: contention probe failed; no constraint." exception=err
            nothing
        end
        if contention !== nothing && contention.valid
            rho = SimulationModel.ParallelCost.lock_duty_cycle(contention, 1)
            if rho >= _rhs_lock_cap_floor_rho()
                p.shared_buffers.rhs_width_ceiling[] = max(1, ceil(Int, 1.0 / rho))
            end
        end
    end

    _calibrate_rhs_plan_if_needed!(p, u_start, args)

    # In-run width identification, AFTER the sweep and only where the sweep
    # produced no plan.
    #
    # Both mechanisms drive `rhs_plan_override`, and running them together does
    # not merely duplicate work -- it corrupts the sweep. The sweep measures a
    # candidate by pinning it and calling the RHS; with a trial installed, the
    # RHS wrapper replaces that pin with whichever arm the trial wants next, so
    # every reading the sweep takes belongs to a different plan than the one it
    # believes it is timing. It then caches that verdict, and the bad entry
    # outlives the run. Measured on the 256-satellite fullstack shape at a
    # 3600 s arc with a clean store: identification-before-sweep ran 57% slower
    # than the sweep alone, 15-0 and significant, and the penalty survived
    # warm-up because it was baked into the cached entry.
    #
    # Ordering it after the sweep also makes the division of labour the useful
    # one. The sweep pins a plan when it can separate the candidates by its
    # margin; identification is for the case where it cannot and abstains,
    # which is where a solve is currently left with the runtime heuristic and no
    # measurement behind it.
    if _rhs_identify_enabled() && p.shared_buffers.rhs_plan_override[] === nothing
        p.shared_buffers.rhs_width_trial[] =
            build_rhs_width_trial(p, args.dynamics_model.dynamic_effectors)
    end

    # Skip per-step solution/dense storage when nothing reads the trajectory.
    # gravity_backbone_split backfills save data from interior solution points,
    # and _auto_stiff_switched/solver metadata read per-step solver state, so
    # those cases keep full storage.  Explicit SPACEAGORA_SOLVER_SAVE_* env
    # settings still override inside the solve helpers.
    #
    # simulation_settings.results is deliberately NOT a reason to keep the
    # solution: the results pipeline builds its output from `saved_values`, the
    # SavingCallback, which records on its own cadence regardless of
    # save_everystep. The only solution access on that path is the final
    # endpoint, which `save_end` governs separately. Requiring full storage here
    # made every results-writing run retain a per-step trajectory nothing read.
    # Verified by writing the same scenario both ways at 256 and 1024
    # satellites: the CSV and feather outputs are byte-identical, while
    # allocation fell 15% (638.6 MB -> 541.0 MB at 1024) and GC went from 10.3%
    # to 0.9%.
    needs_full_solution = return_solution || return_solver_metadata ||
        solver_mode == :gravity_backbone_split

    last_sol = nothing
    solver_trace = NamedTuple[]
    checkpoint_saved_times = Float64[]
    checkpoint_saved_data = SimulationModel.SaveData[]
    backbone_saved_times = Float64[]
    backbone_saved_data = SimulationModel.SaveData[]

    if t_start < mission_end && checkpoint_active
        interval = args.simulation_settings.checkpoint_interval_s
        t_cursor = t_start
        u_cursor = deepcopy(u_start)

        while t_cursor < mission_end
            t_next = min(t_cursor + interval, mission_end)
            empty!(saved_values.t)
            empty!(saved_values.saveval)
            p.shared_buffers.solve_segment_end_time[] = t_next
            if jac_prototype !== nothing
                active_now = count(identity, p.is_active)
                if active_now < jac_prototype_active_count
                    jac_prototype = _build_block_diagonal_jac_prototype(u_cursor, p.is_active)
                    jac_prototype_active_count = active_now
                    # The cached integrator baked in the previous sparsity and
                    # reinit! does not rebuild W, so it would quietly keep using
                    # the old pattern and the narrower one would never take
                    # effect. Drop it and let the next segment re-init.
                    solver_cache === nothing || (solver_cache.integrator = nothing)
                end
            end
            prob = _build_typed_solver_problem(u_cursor, (t_cursor, t_next), p, callbacks, solver_mode, jac_prototype)
            seg_sol, solve_meta = try
                # Every segment resolves the same dtmax and save options, so the
                # cache hits from the second segment on and each one reuses the
                # integrator instead of rebuilding its cache, jac config, W and
                # symbolic factorization. Safe here because the segment state is
                # deepcopied into u_cursor below and _save_snapshot's getters
                # materialise values (fresh SVectors), so nothing retains a
                # reference into the integrator's own state buffer.
                _solve_with_solver_policy(prob, solver_cfg, args, reltol_tol, abstol_tol; solver_cache=solver_cache, needs_full_solution=needs_full_solution)
            catch err
                _try_save_simulation_results_if_enabled!(
                    args,
                    solver_mode,
                    checkpoint_active,
                    save_fields_resolved,
                    saved_values,
                    checkpoint_saved_times,
                    checkpoint_saved_data,
                    backbone_saved_times,
                    backbone_saved_data,
                )
                rethrow()
            end
            push!(solver_trace, solve_meta)
            if solver_mode == :gravity_backbone_split
                _append_backbone_saved_segment!(checkpoint_saved_times, checkpoint_saved_data, seg_sol, save_fields_resolved, p)
            else
                _append_checkpoint_saved_segment!(
                    checkpoint_saved_times,
                    checkpoint_saved_data,
                    saved_values,
                    seg_sol,
                    save_fields_resolved,
                    p,
                )
            end
            last_sol = seg_sol
            if !SciMLBase.successful_retcode(seg_sol.retcode)
                _try_save_simulation_results_if_enabled!(
                    args,
                    solver_mode,
                    checkpoint_active,
                    save_fields_resolved,
                    saved_values,
                    checkpoint_saved_times,
                    checkpoint_saved_data,
                    backbone_saved_times,
                    backbone_saved_data,
                )
                throw(ErrorException("Checkpointed solve failed with retcode=$(seg_sol.retcode)."))
            end
            t_cursor = Float64(seg_sol.t[end])
            u_cursor = deepcopy(seg_sol.u[end])
            _write_checkpoint!(args, t_cursor, u_cursor, string(solver_mode))
            if string(seg_sol.retcode) != "Success" || !_gravity_backbone_time_reached(t_cursor, t_next)
                break
            end
        end

    elseif t_start < mission_end
        p.shared_buffers.solve_segment_end_time[] = mission_end
        prob = _build_typed_solver_problem(u_start, (t_start, mission_end), p, callbacks, solver_mode, jac_prototype)
        sol, solve_meta = try
            _solve_with_solver_policy(prob, solver_cfg, args, reltol_tol, abstol_tol; solver_cache=solver_cache, needs_full_solution=needs_full_solution)
        catch err
            _try_save_simulation_results_if_enabled!(
                args,
                solver_mode,
                checkpoint_active,
                save_fields_resolved,
                saved_values,
                checkpoint_saved_times,
                checkpoint_saved_data,
                backbone_saved_times,
                backbone_saved_data,
            )
            rethrow()
        end
        push!(solver_trace, solve_meta)
        last_sol = sol
        if solver_mode == :gravity_backbone_split
            _append_backbone_saved_segment!(backbone_saved_times, backbone_saved_data, sol, save_fields_resolved, p)
        end
        if !SciMLBase.successful_retcode(sol.retcode)
            _try_save_simulation_results_if_enabled!(
                args,
                solver_mode,
                checkpoint_active,
                save_fields_resolved,
                saved_values,
                checkpoint_saved_times,
                checkpoint_saved_data,
                backbone_saved_times,
                backbone_saved_data,
            )
            throw(ErrorException("Solve failed with retcode=$(sol.retcode)."))
        end
    end

    # Time the solve, so the next process with this signature knows whether the
    # ~0.1 s sweep is worth paying. Placed here rather than at the end of the
    # function so results/telemetry work is not counted as solve cost.
    _rhs_calib_record_solve_time!()

    # Process and save results
    _save_simulation_results_if_enabled!(
        args,
        solver_mode,
        checkpoint_active,
        save_fields_resolved,
        saved_values,
        checkpoint_saved_times,
        checkpoint_saved_data,
        backbone_saved_times,
        backbone_saved_data,
    )

    if return_solution
        if checkpoint_active && args.simulation_settings.checkpoint_interval_s < mission_end
            @warn "return_solution=true with checkpointed integration returns the final segment ODESolution, not a stitched full-history ODESolution."
        end
        if return_solver_metadata
            parallel_policy = try
                SimulationModel.ParallelPolicy.policy_telemetry_snapshot()
            catch
                nothing
            end
            return (
                solution=last_sol,
                solver_mode=string(solver_mode),
                solver_trace=solver_trace,
                parallel_policy=parallel_policy,
                spice_counters=_spice_runtime_counters_snapshot(p)
            )
        end
        return last_sol
    end
    return nothing
    end
end
