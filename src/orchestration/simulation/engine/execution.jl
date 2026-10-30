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
const _CROSS_SATELLITE_INDEX_FIELDS = (:chaser_index, :target_index)

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
    hasproperty(args, :interlink_model) && args.interlink_model !== nothing && return true
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
    # Preserve the explicit sparse wrapper's v6 specialization. Unwrapped
    # problems continue to follow the solver library's normal defaults.
    jac_prototype === nothing && return f
    return ODEFunction{true, SciMLBase.FullSpecialize}(f; jac_prototype=jac_prototype)
end

@inline function _build_typed_solver_problem(u0, tspan, p, callbacks, solver_mode::Symbol,
    jac_prototype::Union{Nothing, SparseMatrixCSC{Float64, Int}}=nothing)
    mode = solver_mode
    if mode in (:gravity_backbone_split, :split_imex, :multirate) &&
            any(SimulationModel.has_battery, p.args.dynamics_model.spacecraft)
        throw(ArgumentError("A spacecraft with a battery needs a solver mode that integrates the full right-hand side; " *
            "solver_mode=$(repr(mode)) splits it. Use :tsit5, :dp8, :rodas5p or :auto_stiff."))
    end
    if mode in (:gravity_backbone_split, :split_imex, :multirate) &&
            any(SimulationModel.has_data_store, p.args.dynamics_model.spacecraft)
        throw(ArgumentError("A spacecraft with a data store needs a solver mode that integrates the full right-hand side; " *
            "solver_mode=$(repr(mode)) splits it. Use :tsit5, :dp8, :rodas5p or :auto_stiff."))
    end
    # The split modes carry the CMG rates in their control component only
    # (`actuator_rates` in `_assign_orientation_rhs!`); the gravity backbone has
    # no attitude state.
    if mode === :gravity_backbone_split &&
            any(sc -> !isempty(sc.root.rw_assembly.cmgs), p.args.dynamics_model.spacecraft)
        throw(ArgumentError("Control moment gyros need the attitude equation, which " *
            "solver_mode=:gravity_backbone_split does not integrate. Use another solver mode."))
    end
    # The gravity backbone runs no sensor sampling, and multirate subsolves
    # revisit overlapping intervals, which would sample twice or out of order.
    if mode in (:gravity_backbone_split, :multirate) &&
            any(sc -> !isempty(SimulationModel.mounted_sensors(sc)), p.args.dynamics_model.spacecraft)
        throw(ArgumentError("Mounted sensors need a solver mode that takes each step once with the run's callbacks; " *
            "solver_mode=$(repr(mode)) does not. Use :tsit5, :dp8, :rodas5p, :auto_stiff or :split_imex."))
    end
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
    elseif mode == :multirate_mri
        _any_robot_arm_effector(p.args) && throw(ArgumentError(
            "Coupled robot-arm effectors are not supported by :multirate_mri; use a full-state solver."))
        # SciML's MRI methods take the fast component first. The original
        # Strang path below retains its historical slow-first order.
        return SplitODEProblem(spacecraft_dynamics_mri_fast!, spacecraft_dynamics_mri_slow!,
            u0, tspan, p; callback=callbacks)
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
        f = ODEFunction{true, SciMLBase.FullSpecialize}(spacecraft_dynamics!; jac_prototype=jac_prototype)
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

"""
    _results_thin_stride()::Int

Opt-in results thinning stride, read from `SPACEAGORA_RESULTS_THIN_STRIDE`
(default `"1"`, meaning no thinning). Any value `<= 1` or unparseable is
treated as `1`. See `_thin_results_segment` for what a stride of `k` keeps.
"""
@inline function _results_thin_stride()::Int
    raw = strip(_engine_env_get("SPACEAGORA_RESULTS_THIN_STRIDE", "1"))
    parsed = tryparse(Int, raw)
    (parsed === nothing || parsed < 1) ? 1 : parsed
end

"""
    _results_final_only()::Bool

Opt-in "final state only" results thinning, read from
`SPACEAGORA_RESULTS_FINAL_ONLY` (default `"0"`/off). Takes priority over
`_results_thin_stride` when both are set.
"""
@inline function _results_final_only()::Bool
    lowercase(strip(_engine_env_get("SPACEAGORA_RESULTS_FINAL_ONLY", "0"))) in ("1", "true", "yes")
end

"""
    _thin_results_segment(times, data)

Opt-in state thinning for the *written results output* (CSV / results
bundle), NOT for the returned `ODESolution` -- `run_simulation`'s
`return_solution=true` path is unaffected by either setting.

With both settings at their default (off), this is the identity: the exact
`times`/`data` vectors already built by the caller are returned unchanged,
so the default output is byte-identical to before this function existed.

With `SPACEAGORA_RESULTS_FINAL_ONLY=1`, only the last saved state is kept.
Otherwise, with `SPACEAGORA_RESULTS_THIN_STRIDE=k` (k > 1), every k-th saved
state is kept, always including the final state (a checkpointed or
callback-driven save cadence need not land exactly on a multiple of k, and a
results file that silently dropped the mission-end value would be a
correctness hazard for anything reading it).

Every element of the returned vectors is an element of the input vectors --
no state is recomputed or copied element-wise -- so a retained state is
always byte-identical to what the default (unthinned) path would have
written for that same saved time.
"""
@inline function _thin_results_segment(
    times::Vector{Float64},
    data::Vector{SimulationModel.SaveData},
)::Tuple{Vector{Float64}, Vector{SimulationModel.SaveData}}
    isempty(times) && return (times, data)
    if _results_final_only()
        return ([times[end]], [data[end]])
    end
    stride = _results_thin_stride()
    stride <= 1 && return (times, data)
    n = length(times)
    idxs = collect(1:stride:n)
    idxs[end] == n || push!(idxs, n)
    return (times[idxs], data[idxs])
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
    resumed_from_s::Union{Nothing, Float64},
    rows_log=nothing;
    return_table::Bool=false,
    write_files::Bool=true,
)
    (args.simulation_settings.results || return_table) || return nothing
    # A fallback that fails before its first checkpoint saves only its new rows
    # in its own results directory; the old checkpoint streams remain untouched.
    if checkpoint_active && rows_log !== nothing && rows_log.stream !== nothing && !rows_log.keep_earlier
        return _save_streamed_results!(args, rows_log, save_fields_resolved, checkpoint_saved_times,
            checkpoint_saved_data, resumed_from_s; return_table=return_table, write_files=write_files)
    end
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
    results_times, results_data = _thin_results_segment(
        convert(Vector{Float64}, results_times),
        convert(Vector{SimulationModel.SaveData}, results_data),
    )
    sink = return_table ? _LONG_TABLES_SINK[] : nothing
    if args.simulation_settings.results_layout === :long
        tables = _build_long_results_tables(results_times, results_data, save_fields_resolved, args)
        sink === nothing || (sink[] = tables)
        csv_paths = (spacecraft=nothing, run=nothing)
        schedule_files = nothing
        if write_files && args.simulation_settings.results
            schedule_files = _save_link_schedule_if_enabled!(args)
            csv_paths = _write_long_results!(tables, results_times, args; write_bundle=_typed_save_bundle_enabled(),
                resumed_from_s=resumed_from_s, link_schedule=schedule_files)
        end
        return return_table ? (; csv_path=csv_paths.spacecraft, run_csv_path=csv_paths.run, table=tables.spacecraft,
            run_table=tables.run, schedule_files) : csv_paths.spacecraft
    end
    sink === nothing || (sink[] = _build_long_results_tables(results_times, results_data, save_fields_resolved, args))
    results_df = _build_results_dataframe(results_times, results_data, save_fields_resolved, args)
    csv_path = nothing
    schedule_files = nothing
    if write_files && args.simulation_settings.results
        if args.simulation_settings.save_csv
            csv_path = _write_results_csv!(results_df, args)
        end
        schedule_files = _save_link_schedule_if_enabled!(args)
        if _typed_save_bundle_enabled()
            _write_results_bundle!(results_df, results_times, args; csv_path=csv_path, resumed_from_s=resumed_from_s,
                link_schedule=schedule_files)
        end
    end
    return return_table ? (; csv_path, table=results_df, schedule_files) : csv_path
end

# Persist only the schedule written by this save, so stale files are never listed.
function _save_link_schedule_if_enabled!(args)
    args.interlink_model === nothing && return nothing
    schedule = link_schedule(args.interlink_model; policy=args.scheduling_policy_model)
    written = _write_link_schedule!(schedule, args; csv=args.simulation_settings.save_csv,
        feather=_typed_save_bundle_enabled())
    return (; written..., metadata=_link_schedule_manifest(schedule, args.scheduling_policy_model,
        args.interlink_model))
end

# The results of a run that streams them (`_stream_results_rows!`): the rows still
# in memory join the streams, then the feathers and the manifest are written from
# them, and the returned tables read the feathers without copying them. A resume
# that solves nothing writes nothing: it returns the rows its checkpoint records
# in the streams (the first `written` saved times) and the ones it keeps, read
# into memory. The one exception is a resume that took back the streams of a run
# killed while writing its files (`_start_results_rows!`; `args` are then that
# run's): it writes them.
function _save_streamed_results!(args, log, save_fields, times_acc::Vector{Float64}, data_acc::Vector,
                                 resumed_from_s::Union{Nothing, Float64}; return_table::Bool, write_files::Bool)
    stream = log.stream
    sink = return_table ? _LONG_TABLES_SINK[] : nothing
    in_memory(df) = DataFrame(df; copycols=true)
    if !write_files && (stream.fresh || !_streamed_results_unfinished(args))
        keep = Int[i for i in eachindex(times_acc) if _results_row_written(stream.released + i, i == length(times_acc))]
        kept = _build_long_results_tables(times_acc[keep], data_acc[keep], save_fields, args)
        per_time = (length(args.dynamics_model.spacecraft), 1)
        streamed(path, i) = stream.bytes[i] > 0 ?
            in_memory(SimulationModel.IOOutputs._read_long_results_table(path)[1:(stream.written * per_time[i]), :]) : nothing
        joined(a, b) = a === nothing ? b : b === nothing ? a : vcat(a, b)
        tables = (spacecraft=joined(streamed(stream.spacecraft, 1), kept.spacecraft),
                  run=joined(streamed(stream.run, 2), kept.run))
        sink === nothing || (sink[] = tables)
        return return_table ? (; csv_path=nothing, run_csv_path=nothing, table=tables.spacecraft, run_table=tables.run, schedule_files=nothing) :
            nothing
    end
    _stream_results_rows!(log, times_acc, data_acc, save_fields, args; final=true)
    # Final-write recovery may regenerate a schedule only from validated history.
    # Without it, preserve earlier schedule bytes without claiming ownership.
    has_history = args.interlink_model !== nothing && !isempty(args.interlink_model.history)
    schedule_files = (write_files || has_history) ? _save_link_schedule_if_enabled!(args) : nothing
    written = _finish_long_stream!((spacecraft=stream.spacecraft, run=stream.bytes[2] > 0 ? stream.run : nothing), args;
        steps=stream.written, first_time_s=isnan(stream.first_s) ? nothing : stream.first_s, resumed_from_s=resumed_from_s,
        link_schedule=schedule_files)
    mapped(path) = path === nothing ? nothing : SimulationModel.IOOutputs._read_long_results_table(path)
    tables = (spacecraft=mapped(written.feather.spacecraft), run=mapped(written.feather.run))
    sink === nothing || (sink[] = (spacecraft=in_memory(tables.spacecraft),
                                   run=tables.run === nothing ? nothing : in_memory(tables.run)))
    return return_table ? (; csv_path=written.csv.spacecraft, run_csv_path=written.csv.run, table=tables.spacecraft,
        run_table=tables.run, schedule_files) : written.csv.spacecraft
end

function _try_save_simulation_results_if_enabled!(args...)
    try
        return _save_simulation_results_if_enabled!(args...)
    catch err
        @warn "Failed to save partial simulation results after solve failure." reason=sprint(showerror, err)
        return nothing
    end
end

# One solve of the run, from `u0` over `tspan`, under the run's solver policy
# (`setup`: what every solve of the run shares), recorded in `solver_trace`. When
# the solve throws, the results so far are saved first (`results_inputs`) and the
# exception is rethrown as it was, with its backtrace.
function _solve_run_interval!(solver_trace, setup, u0, tspan, jac_prototype, results_inputs)
    prob = _build_typed_solver_problem(u0, tspan, setup.p, setup.callbacks, setup.solver_mode, jac_prototype)
    sol, solve_meta = try
        _solve_with_solver_policy(prob, setup.solver_cfg, setup.args, setup.reltol, setup.abstol;
            solver_cache=setup.solver_cache, needs_full_solution=setup.needs_full_solution)
    catch
        _try_save_simulation_results_if_enabled!(results_inputs...)
        rethrow()
    end
    push!(solver_trace, solve_meta)
    return sol
end

# A solve that ended without success: the results so far are saved, then an
# `ErrorException` names `what` and the solve's return code.
function _require_successful_solve(sol, results_inputs, what)
    SciMLBase.successful_retcode(sol.retcode) && return nothing
    _try_save_simulation_results_if_enabled!(results_inputs...)
    throw(ErrorException("$what failed with retcode=$(sol.retcode)."))
end

function _with_density_model_epoch(args::SimulationConfiguration)
    environment = args.environment_model
    density_model = SimulationModel.with_density_model_epoch(environment.density_model, args.initial_time)
    density_model === environment.density_model && return args
    aligned_environment = SimulationModel.SimConfig._with_fields(environment; density_model)
    return SimulationModel.SimConfig._with_configuration(args;
        environment_model=aligned_environment)
end

_declares_gnc(sc) = !(isempty(sc.guidance.guidance_effectors) &&
                      isempty(sc.navigation.navigation_effectors) &&
                      isempty(sc.control.control_effectors))

# Spacecraft copy without its per-spacecraft GNC (links, joints and the rest are shared).
_without_gnc(sc::SimulationModel.SpacecraftModel) = SimulationModel.SimConfig._with_fields(sc;
    guidance=SimulationModel.GuidanceModel((), Float64[]), navigation=SimulationModel.NavigationModel((), Float64[]),
    control=SimulationModel.ControlModel((), Float64[]))

"""
Fold per-spacecraft GNC declarations (`SpacecraftModel` `guidance`/`navigation`/`control`) into the
configuration-level GNC tuples, once, at the start of `run_simulation`: configuration-level effectors
first, then each spacecraft's effectors bound to its position by `bind_spacecraft`, in spacecraft order.
The returned configuration's spacecraft declare no GNC, so a second pass is a no-op, and the caller's
spacecraft are not modified. Returns `args` itself when no spacecraft declares GNC.
"""
function _flatten_spacecraft_gnc(args::SimulationConfiguration)
    spacecraft = args.dynamics_model.spacecraft
    any(_declares_gnc, spacecraft) || return args
    g = Any[args.guidance_model.guidance_effectors...]; gr = copy(args.guidance_model.guidance_rates)
    n = Any[args.navigation_model.navigation_effectors...]; nr = copy(args.navigation_model.navigation_rates)
    c = Any[args.control_model.control_effectors...]; cr = copy(args.control_model.control_rates)
    warned = Set{Any}()
    for (idx, sc) in enumerate(spacecraft)
        append!(g, _bind_declared_gnc(e, idx, sc.id, warned) for e in sc.guidance.guidance_effectors)
        append!(gr, sc.guidance.guidance_rates)
        append!(n, _bind_declared_gnc(e, idx, sc.id, warned) for e in sc.navigation.navigation_effectors)
        append!(nr, sc.navigation.navigation_rates)
        append!(c, _bind_declared_gnc(e, idx, sc.id, warned) for e in sc.control.control_effectors)
        append!(cr, sc.control.control_rates)
    end
    return SimulationModel.SimConfig._with_configuration(args;
        guidance_model=SimulationModel.GuidanceModel(Tuple(g), gr),
        navigation_model=SimulationModel.NavigationModel(Tuple(n), nr),
        control_model=SimulationModel.ControlModel(Tuple(c), cr),
        dynamics_model=SimulationModel.DynamicsModel(
            SimulationModel.SpacecraftModel[_declares_gnc(sc) ? _without_gnc(sc) : sc for sc in spacecraft],
            args.dynamics_model.dynamic_effectors, args.dynamics_model.effector_scopes))
end

# A model a spacecraft declares, bound to it. A bound copy whose scope still names
# every spacecraft acts on every spacecraft, as one declared at configuration level
# does, which is rarely what declaring it on one meant: warned once per type and run
# (each member of an ensemble is a run of its own).
function _bind_declared_gnc(effector, sat_idx::Int, id, warned::Set{Any})
    bound = SimulationLifecycle.bind_spacecraft(effector, sat_idx)
    T = typeof(bound)
    if GNCRuntime.spacecraft_scope(bound) isa GNCRuntime.AllSpacecraft && !(T in warned)
        push!(warned, T)
        @warn "$(nameof(T)) is declared on the spacecraft at run position $sat_idx (id $id) and bound to it by " *
            "bind_spacecraft, but its spacecraft_scope is AllSpacecraft(), so it acts on every spacecraft. Return " *
            "BoundSpacecraft(($sat_idx,)) from spacecraft_scope for the copy bound to run position $sat_idx."
    end
    return bound
end

# Opt-in factories own fresh callback state and capture the configuration that
# actually runs. The empty default is handled at the call site, preserving the
# existing extra_callbacks path and construction order without additional work.
function _callbacks_from_factories(extra_callbacks, factories, args::SimulationConfiguration)
    factories isa Union{Tuple, AbstractVector} || throw(ArgumentError(
        "extra_callback_factories must be a tuple or vector of callback factories."))
    for (i, factory) in enumerate(factories)
        applicable(factory, args) || throw(ArgumentError(
            "extra_callback_factories[$i] must accept the isolated SimulationConfiguration."))
    end
    built = map(enumerate(Tuple(factories))) do (i, factory)
        callback = factory(args)
        callback isa DiffEqBase.DECallback || throw(ArgumentError(
            "extra_callback_factories[$i] must return a solver callback or CallbackSet; got $(typeof(callback))."))
        callback
    end
    return (extra_callbacks..., built...)
end

# Print the initial state as labelled values (the state's own units: m, m/s, kg,
# J/cm^2 for heat loads), one line per spacecraft, rather than the state array's
# full type. Field names follow `u.sc[i]` (see the
# integrated-state documentation); other top-level components print by name.
function _print_initial_conditions(io::IO, u0::ComponentVector)
    println(io, "Initial conditions:")
    _plain(value) = value isa AbstractArray ? collect(value) : value
    for name in keys(u0)
        component = getproperty(u0, name)
        if name === :sc && component isa AbstractVector
            for (i, sc) in enumerate(component)
                fields = (string(field, " = ", _plain(getproperty(sc, field))) for field in keys(sc))
                println(io, "  spacecraft ", i, ": ", join(fields, ", "))
            end
        else
            println(io, "  ", name, " = ", _plain(component))
        end
    end
    return nothing
end

# The start of a run, in three phases that `run_simulation` calls in this order
# (see "Run state: owners, resets and restores" in
# docs/contributing/architecture.md). The engine's and the callbacks' run state
# comes back from a checkpoint before the models start, and the models' own run
# state after their initializers: two restores at two different times.

# Phase 1, a resumed run's engine and callback state: read the checkpoint, check
# that it fits the run, restore the spacecraft's activity, the interlinks' and the
# callbacks' run state, and the pending periodic ticks, and set the resume flags
# the callbacks (`callback_resume`) and the models' initializers (`model_resume`)
# read. The models' own state (`model_snapshot`) and their spacecraft scopes
# (`scope_entry`) are only checked here. Without a checkpoint to resume from, the
# run starts at 0 from `initial_conditions`.
function _resume_engine_state!(args, p, solver_mode, initial_conditions)
    t_start = 0.0
    u_start = initial_conditions
    resumed, model_snapshot, scope_entry = false, nothing, nothing
    pending_ticks = (at_s=NaN, from_s=-Inf, clocks=nothing)
    results_rows = nothing
    interlink_history = nothing
    link_attitudes = nothing
    callback_checkpoint = nothing
    forecast_restart = _CONTROL_FORECAST_RESTART[]
    if args.simulation_settings.resume_from_checkpoint || forecast_restart !== nothing
        ckpt = forecast_restart === nothing ? _load_checkpoint(args) : forecast_restart.checkpoint
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
            # A gravity_backbone_split run never has interlinks (_validate_interlinks!
            # refuses them), and its state has no `sc` components to check.
            solver_mode != :gravity_backbone_split && _check_interlink_checkpoint_state(args, ckpt)
            _check_effector_scope_checkpoint_state(args, ckpt)
            scope_entry = _gnc_scope_checkpoint_entry(args, ckpt)
            link_attitudes = _check_link_attitude_checkpoint_state(args, ckpt)
            callback_checkpoint = ckpt
            _restore_activity_checkpoint_state!(args, ckpt, p.is_active)
            _restore_interlink_checkpoint_state!(args, ckpt, p.is_active)
            _restore_callback_checkpoint_state!(args, ckpt, p)
            _restore_link_attitude_checkpoint_state!(args, link_attitudes)
            _restore_maneuver_buffer_checkpoint_state!(ckpt, p)
            _restore_orbit_counter_checkpoint_state!(ckpt, p)
            _restore_entry_counter_checkpoint_state!(ckpt, p)
            model_snapshot = _check_model_checkpoint_state(args, ckpt)
            pending_ticks = _pending_periodic_ticks(ckpt)
            results_rows = _results_rows_checkpoint_state(ckpt)
            interlink_history = _interlink_history_checkpoint_state(ckpt)
            resumed = true
            t_start = ckpt.t
            u_start = _resumed_checkpoint_state(args, ckpt.u)
            _initialize_in_atmosphere_flags!(p, u_start)
            p.shared_buffers.model_resume[] = model_snapshot !== nothing
            if args.simulation_settings.verbose
                println("Resuming simulation from checkpoint at t=$(round(t_start, digits=6)) s")
            end
        end
    end
    return (; resumed, model_snapshot, scope_entry, pending_ticks, t_start, u_start, results_rows, interlink_history, link_attitudes, callback_checkpoint)
end

# Phase 2, the run clock of the periodic navigation, guidance and control ticks
# (see callbacks/periodic_tick_callbacks.jl), and the ticks a resumed run takes at
# its start.
function _start_run_clock!(p, run_clock, mission_end, pending_ticks)
    forecast_restart = _CONTROL_FORECAST_RESTART[]
    origin = forecast_restart === nothing ? 0.0 : forecast_restart.clock_origin_s
    p.shared_buffers.periodic_clock_origin_s[] = run_clock ? origin : NaN
    p.shared_buffers.periodic_clock_horizon_s[] = run_clock ? Float64(mission_end) : NaN
    p.shared_buffers.periodic_clock_pending_s[] = pending_ticks.at_s
    p.shared_buffers.periodic_clock_pending_from_s[] = pending_ticks.from_s
    p.shared_buffers.periodic_clock_pending_clocks[] = pending_ticks.clocks
    return nothing
end

# Phase 3, the models: the links, the control clocks, the models' initializers,
# then, on a resumed run, callback state, model-owned state and exact link
# attitudes from the checkpoint replace the initialization values. Returns the models' spacecraft scopes as the initializers
# left them, which a checkpointed run writes with its checkpoints (`nothing` for a
# run that neither writes nor resumes one).
function _start_run_models!(args, p, u_start, t_start, run_clock, resumed, model_snapshot, scope_entry,
                            checkpoint_active; link_attitudes=nothing, callback_checkpoint=nothing)
    # The links first (each initialized once), then the control clocks, then
    # navigation, guidance and control: each stage may read what the previous one
    # set up.
    SimulationModel.SimulationCallbacks.initialize_gnc_links!(args, u_start, p, t_start)
    SimulationModel.SimulationCallbacks._report_gnc_clocks!(args,
        run_clock ? p.shared_buffers.periodic_clock_origin_s[] : t_start)
    for navigation in args.navigation_model.navigation_effectors
        SimulationLifecycle.initialize_navigation!(navigation, u_start, p, t_start)
    end
    for guidance in args.guidance_model.guidance_effectors
        SimulationLifecycle.initialize_guidance!(guidance, u_start, p, t_start)
    end
    for control in args.control_model.control_effectors
        SimulationLifecycle.initialize_control!(control, u_start, p, t_start)
    end
    # The scopes the initializers resolved (the run's GNC scopes, which the ticks,
    # the right-hand side and the readers follow from here on), compared with the
    # checkpoint's before the models' own run state replaces their reset on a
    # resumed run.
    scopes = _gnc_run_scopes(args)
    gnc_scopes = checkpoint_active || resumed ? _gnc_scope_record(args, scopes) : nothing
    resumed && _check_gnc_scope_baseline(args, scope_entry, gnc_scopes)
    p.shared_buffers.gnc_scopes[] = scopes
    # Initializers can command mounted hardware. Reinstate its saved state,
    # then restore model-owned state (including custom sensor connections),
    # preserving the generic-before-custom restoration order.
    resumed && callback_checkpoint !== nothing && _restore_callback_checkpoint_state!(args, callback_checkpoint, p)
    resumed && _restore_model_checkpoint_state!(args, model_snapshot)
    # Hinge restoration recomputes quaternions; the saved geometry is exact.
    resumed && _restore_link_attitude_checkpoint_state!(args, link_attitudes)
    checkpoint_active && _check_model_checkpoint_hooks(args)
    SimulationModel.ControlHooks._edg_initialize_heat_accounting!(args, u_start, Float64(t_start))
    return gnc_scopes
end

"""
    run_simulation(args...; isolate_state=true, kwargs...)

Stable package entrypoint for simulation execution used by calibration integrations.

Density models are first aligned to `args.initial_time` through
[`with_density_model_epoch`](@ref). Known GRAM construction settings are retained;
unknown raw-core realignment and fixed-surrogate epoch changes are rejected.

By default, `isolate_state=true` deep-copies the simulation configuration after epoch alignment, before execution.
This preserves correctness and reentrancy across repeated runs and concurrent callers by
preventing one run from mutating shared campaign or model state that another run still
references.

Set `isolate_state=false` only as an advanced performance lever when the caller owns the
configuration instance and will not reuse it concurrently or across runs that may mutate
shared state. This can reduce setup cost for large mission definitions or many short runs,
but it trades away the default isolation guarantee. When any spacecraft declares GNC, the
configuration that runs is a flattened copy. With `isolate_state=false`, configuration-level
effectors remain the caller's objects. Per-spacecraft binding may return the original effector
when its index already matches, or a shallow copy that shares mutable members. Scalar field
updates on a copied effector are not reflected in the original, while mutations of shared
members are. GNC-stripped spacecraft likewise share their links, joints and initial conditions.
Keep the default isolation when the caller needs its state preserved; use the returned
`SimulationResults.configuration` with `return_results=true` to inspect the objects that ran.

With `return_results=true` the call returns a [`SimulationResults`](@ref): the
results table built in memory (same columns and values as the CSV, even when
`simulation_settings.results=false`, in which case result files are not written), the
configuration that actually ran (the deep copy under `isolate_state=true`, so
controller state and logs are reachable), the files the run wrote, and the
solution when `return_solution=true` is also given (otherwise `nothing`).
`return_results=true` cannot be combined with `return_solver_metadata=true`. With
the default `return_results=false` the return value is unchanged.

`extra_callback_factories=()` optionally supplies functions of the isolated
configuration. Each is called once after epoch alignment, state isolation and
model preflight, before run output is created, and must return a fresh solver
callback or callback set. Factory-created callbacks are appended after
`extra_callbacks`; the usual built-in callback ordering is unchanged. Factories
must allocate their captured mutable state per call, avoid I/O, and refrain from
mutating configuration during construction. Their callback closures can then
refer to the same controller/model objects used by the integrator. The empty
default retains the existing callback path. `priority_callback_factories=()`
uses the same isolated factory contract but installs callbacks before normal
guidance, navigation, control and result saving. Engine impact handling and
external-state synchronization retain their ordering. Operational mission
segments use this hook; continuous events still follow solver event ordering.

A checkpointed run keeps its saved rows beside its checkpoint, so a resumed run's
results files and table hold the whole run; a resume that finds no usable rows
warns and writes its results to `resumed_at_<t>s` inside `results_directory`.

Explicitly enabled checkpoint writing remains independent of result-file output.
With `return_solver_metadata=true`, a checkpointed run returns one `solver_trace`
entry per segment solved, and the `retcode` (and `solution`, with
`return_solution=true`) of the last segment solved. A run resumed from a
checkpoint at or past the mission end solves nothing: its `retcode` and `solution`
are `nothing` and its `solver_trace` is empty, so test the mission rather than
these fields.

# Examples
```jldoctest
julia> args = spaceagora_no_gram_example_args();

julia> sol = run_simulation(args; return_solution=true);

julia> length(sol.t) > 1
true
```
"""
function run_simulation(
    args::SimulationConfiguration;
    isolate_state::Bool=true,
    return_solution::Bool=false,
    return_solver_metadata::Bool=false,
    return_results::Bool=false,
    save_fields=nothing,
    extra_callbacks=(),
    extra_callback_factories::Union{Tuple, AbstractVector}=(),
    priority_callback_factories::Union{Tuple, AbstractVector}=(),
    solver_cache::Union{Nothing, SolverIntegratorCache}=nothing,
    visualization::Bool=(_engine_env_get("SPACEAGORA_VISUALIZATION", "0") == "1")
)
    visualization && SimulationModel.SceneVisualization.VisualizerServices.require_available()
    return_results && return_solver_metadata &&
        throw(ArgumentError("return_results=true and return_solver_metadata=true cannot be combined."))
    args.simulation_settings.results_layout in (:wide, :long) || throw(ArgumentError(
        "SimulationSettings.results_layout must be :wide or :long; got :$(args.simulation_settings.results_layout)."))
    args = _flatten_spacecraft_gnc(args)
    # Scoped dynamic effectors become bare effectors and a table of run positions,
    # before any check reads the effectors by type (see effector_scopes.jl).
    args = _resolve_effector_scopes(args)
    for guidance in args.guidance_model.guidance_effectors
        SimulationLifecycle.preflight_guidance(guidance, args; isolate_state=isolate_state)
    end
    for control in args.control_model.control_effectors
        SimulationLifecycle.preflight_control(control, args; isolate_state=isolate_state)
    end
    SimulationModel.SimulationCallbacks.check_control_periods(args)
    # After the models' own checks, so that their more specific errors come first.
    SimulationModel.SimulationCallbacks.check_spacecraft_scopes(args)
    # SolverConfig(parallel=true): re-enter under the flag's scoped environment
    # (see `_with_parallel_flag`). Inside it the flag reads as resolved, so this
    # branch is taken once; with the flag off nothing here runs.
    if args.solver_config !== nothing && _parallel_flag_applies(args.solver_config.parallel)
        return _with_parallel_flag(true) do
            run_simulation(args; isolate_state=isolate_state, return_solution=return_solution,
                           return_solver_metadata=return_solver_metadata, return_results=return_results,
                           save_fields=save_fields, extra_callbacks=extra_callbacks,
                           extra_callback_factories=extra_callback_factories,
                           priority_callback_factories=priority_callback_factories, solver_cache=solver_cache,
                           visualization=visualization)
        end
    end
    # After the re-entry above, so that each link's `check_link` runs once per run.
    SimulationModel.SimulationCallbacks.check_gnc_links(args)
    try
    return SimulationModel.ParallelPolicy.with_policy_context() do
    # `visualization=true` (or SPACEAGORA_VISUALIZATION=1, so an unmodified
    # example script can opt in) turns the scene sidecar on for this run and
    # builds the viewer page once the results are written (see SceneVisualization).
    args = visualization ? SimulationModel.SceneVisualization.with_visualization_scene(args, true) : args
    args = _with_density_model_epoch(args)
    # Isolate mutable campaign/model state by default so repeated/concurrent runs
    # do not alias shared in-memory objects.
    args = isolate_state ? deepcopy(args) : args
    # Resolve effective solver config: explicit field on args wins; otherwise read from env
    # (which respects any active SimulationEngineConfig overrides via _engine_env_get).
    solver_cfg = _run_solver_config(args)
    # MRI phase callbacks need the effective macro-step cap, including when
    # it came from the environment rather than an explicit configuration.
    if solver_cfg.solver_mode === :multirate_mri && args.solver_config === nothing
        args = SimulationModel.SimConfig._with_configuration(args; solver_config=solver_cfg)
    end
    solver_mode = solver_cfg.solver_mode
    solver_mode === :ks && _validate_ks_config!(args)

    _validate_run_models!(args, solver_mode)
    if args.interlink_model !== nothing
        empty!(args.interlink_model.history)
        graph = args.interlink_model.graph
        fill!(graph.available, false)
        fill!(graph.active, false)
        args.scheduling_policy_model.matching = nothing
    end
    if !isempty(extra_callback_factories)
        extra_callbacks = _callbacks_from_factories(extra_callbacks, extra_callback_factories, args)
    end
    priority_callbacks = isempty(priority_callback_factories) ? () :
        _callbacks_from_factories((), priority_callback_factories, args)
    try
        SimulationModel.ParallelPolicy.reset_policy_telemetry!()
        if SimulationModel.ParallelPolicy.persistent_hints_state_reset_requested()
            SimulationModel.ParallelPolicy.reset_persistent_hint_state!()
        end
    catch
    end

    # Set up the model and initial conditions
    initial_conditions = build_initial_conditions(args)
    args.simulation_settings.verbose && _print_initial_conditions(stdout, initial_conditions)

    # Define the ODE parameters and callbacks
    p = _prepare_run_runtime!(args, initial_conditions, solver_mode)
    save_fields_resolved = isnothing(save_fields) ? SimulationModel.default_save_fields(args) : collect(save_fields)
    # Explicit save_fields built before the visualization flag was applied
    # (e.g. `vcat(default_save_fields(args), ...)` in an example run under
    # SPACEAGORA_VISUALIZATION=1) would silently miss the link poses.
    for extra in SimulationModel.SimulationCallbacks.visualization_save_fields(args)
        any(field -> field.name === extra.name, save_fields_resolved) && continue
        save_fields_resolved = vcat(save_fields_resolved, [extra])
    end
    save_field_names = Symbol[field.name for field in save_fields_resolved]
    length(unique(save_field_names)) == length(save_field_names) || throw(ArgumentError("save_fields names must be unique. Got $(save_field_names)."))
    saved_values = SavedValues(Float64, SimulationModel.SaveData)
    callbacks = SimulationModel.get_callbacks(
        length(args.dynamics_model.spacecraft),
        args.dynamics_model.dynamic_effectors,
        args;
        saved_values=saved_values,
        save_fields=save_fields_resolved,
        extra_callbacks=extra_callbacks,
        priority_callbacks=priority_callbacks,
        record_saved_values=return_results,
        solver_mode=solver_mode
    ) # Get the callbacks based on the number of satellites and the dynamic effectors being used in the simulation
    if args.interlink_model !== nothing
        initialization_active = Ref(copy(args.interlink_model.graph.active))
        callbacks = CallbackSet(SimulationModel.SimulationCallbacks.interlink_scheduler_initialization(; initialization_active), callbacks,
            SimulationModel.interlink_scheduler_callback(; initialization_active))
    end
    # Last, after every wrapper: building a CallbackSet from another flattens it
    # back into tuples.
    callbacks = _solver_callback_set(callbacks)
    mission_end = _prepare_run_ephemerides!(args, p)
    checkpoint_active = _typed_checkpoint_enabled(args)
    if checkpoint_active && args.simulation_settings.checkpoint_interval_s <= 0.0
        throw(ArgumentError("SimulationSettings.checkpoint_interval_s must be > 0 when checkpointing is enabled."))
    end

    # The run's start, in three phases (defined above): a resumed run's engine and
    # callback state, the run clock, then the models.
    (; resumed, model_snapshot, scope_entry, pending_ticks, t_start, u_start, results_rows, interlink_history, link_attitudes, callback_checkpoint) =
        _resume_engine_state!(args, p, solver_mode, initial_conditions)


    # println("Initial conditions:")
    # println(initial_conditions)
    # println("ODE parameters:")
    # println(p)
    # println("args.mission_configuration.mission_time: $(args.mission_configuration.mission_time)")
    # A checkpointed or resumed run keeps one tick grid for its periodic
    # navigation, guidance and control, from the run's start at 0 to its end, and
    # ticks at segment boundaries on it; a resumed run takes the ticks pending
    # at its start.
    run_clock = checkpoint_active || resumed
    _start_run_clock!(p, run_clock, mission_end, pending_ticks)
    gnc_scopes = _start_run_models!(args, p, u_start, t_start, run_clock, resumed, model_snapshot, scope_entry,
        checkpoint_active; link_attitudes, callback_checkpoint)
    p.shared_buffers.solve_segment_end_time[] = mission_end
    # prob_debug exists only to feed the NaN-probe below, which itself only
    # runs when SPACEAGORA_DEBUG_INITIAL_DERIVATIVE is set. Building it
    # unconditionally meant every solve -- debug flag on or off -- paid an
    # ODEProblem allocation whose only reader is a branch almost no run takes.
    # Nothing outside these two debug branches reads prob_debug, so deferring
    # its construction into the branch that needs it changes no observable
    # behavior.
    if p.shared_buffers.debug_initial_derivative[] && solver_mode != :gravity_backbone_split
        prob_debug_state = solver_mode == :gravity_backbone_split ? initial_conditions : u_start
        prob_debug = ODEProblem(spacecraft_dynamics!, prob_debug_state, (t_start, mission_end), p, callback=callbacks)
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

    _prepare_run_execution_policy!(p, u_start, args)

    # Skip per-step solution/dense storage when nothing reads the trajectory.
    # gravity_backbone_split backfills save data from interior solution points,
    # so it keeps full storage.  Explicit SPACEAGORA_SOLVER_SAVE_* env settings
    # still override inside the solve helpers.
    #
    # return_solver_metadata is deliberately NOT a reason to keep the solution
    # either. Every field of the metadata tuple is produced beside the solve
    # rather than read out of it: solver_mode and solver_trace come from
    # _solve_with_solver_policy's own return, parallel_policy from the policy
    # telemetry snapshot, spice_counters from the shared buffers, and retcode
    # from the ODESolution object, which exists at any storage setting. The one
    # metadata field that does read per-step solver state -- :auto_stiff's
    # fallback_used, from _auto_stiff_switched(sol).alg_choice -- is computed
    # inside the solver_policy branch that pins its own per-step storage on
    # regardless of this flag, so it is unaffected. Requiring full storage here
    # meant a caller that wanted a retcode and nothing else paid for a whole
    # trajectory, and (because the metadata was only reachable under
    # return_solution) got `nothing` back for it.
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
    needs_full_solution = return_solution || solver_mode == :gravity_backbone_split

    last_sol = nothing
    checkpoint_written = false
    solver_trace = NamedTuple[]
    checkpoint_saved_times = Float64[]
    checkpoint_saved_data = SimulationModel.SaveData[]
    backbone_saved_times = Float64[]
    backbone_saved_data = SimulationModel.SaveData[]
    # A checkpointed run keeps its saved rows beside its checkpoints; a resumed
    # one takes back the rows up to its checkpoint first, so its results hold the
    # whole run (see `_start_results_rows!`).
    rows_recorded = solver_mode == :gravity_backbone_split || args.simulation_settings.results || return_results
    # A resume that takes ticks at its start saves its first row again after them.
    # The backbone's rows come from the solution, which has no row of its own there.
    retake_last_row = resumed && t_start < mission_end && solver_mode != :gravity_backbone_split &&
        pending_ticks.clocks !== nothing && !isempty(pending_ticks.clocks)
    # Likewise the interlink history, so that a resumed run's link schedule holds
    # the whole run (see `_start_interlink_history!`).
    history_log = checkpoint_active ? _start_interlink_history!(args, resumed, t_start, interlink_history;
        solves=t_start < mission_end) : nothing
    rows_log = checkpoint_active ? _start_results_rows!(args, rows_recorded, save_fields_resolved, resumed, t_start,
        results_rows, checkpoint_saved_times, checkpoint_saved_data; solves=t_start < mission_end,
        returns_table=return_results, retake_last_row) : nothing
    results_args = _results_destination(args, rows_log)
    # The results so far, as both the partial saves after a failed solve and the
    # final save read them: containers changed only in place, never rebound.
    results_inputs = (results_args, solver_mode, checkpoint_active, save_fields_resolved, saved_values,
        checkpoint_saved_times, checkpoint_saved_data, backbone_saved_times, backbone_saved_data,
        resumed ? t_start : nothing, rows_log)
    # What every solve of the run shares (`_solve_run_interval!`).
    solve_setup = (; p, callbacks, solver_mode, solver_cfg, args, reltol=reltol_tol, abstol=abstol_tol, solver_cache,
        needs_full_solution)

    if t_start < mission_end && checkpoint_active
        interval = args.simulation_settings.checkpoint_interval_s
        t_cursor = t_start
        u_cursor = deepcopy(u_start)

        while t_cursor < mission_end
            t_next = _checkpoint_segment_end(t_cursor, interval, mission_end)
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
            _initialize_in_atmosphere_flags!(p, u_cursor)
            # Segments with the same initial phase and save options reuse
            # the cached integrator; a phase change rebuilds it. Reuse avoids
            # rebuilding its cache, jac config, W and symbolic factorization. Safe here because the segment state is
            # deepcopied into u_cursor below and _save_snapshot's getters
            # materialise values (fresh SVectors), so nothing retains a
            # reference into the integrator's own state buffer.
            seg_sol = _solve_run_interval!(solver_trace, solve_setup, u_cursor, (t_cursor, t_next), jac_prototype, results_inputs)
            if solver_mode == :gravity_backbone_split
                _append_backbone_saved_segment!(checkpoint_saved_times, checkpoint_saved_data, seg_sol, save_fields_resolved, p)
            else
                _retake_restored_row!(rows_log, checkpoint_saved_times, checkpoint_saved_data, saved_values.t)
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
            _require_successful_solve(seg_sol, results_inputs, "Checkpointed solve")
            t_cursor = Float64(seg_sol.t[end])
            u_cursor = deepcopy(seg_sol.u[end])
            # The rows first: a checkpoint never records rows that are not in the file.
            # A run that streams its results writes all but the last row now.
            _preserve_earlier_results!(rows_log)
            _stream_results_rows!(rows_log, checkpoint_saved_times, checkpoint_saved_data, save_fields_resolved,
                results_args)
            rows_entry = _append_results_rows!(rows_log, checkpoint_saved_times, checkpoint_saved_data)
            history_entry = _append_interlink_history!(history_log, args.interlink_model)
            _write_checkpoint!(args, t_cursor, u_cursor, string(solver_mode); is_active=p.is_active,
                callback_state=_callback_checkpoint_state(args), model_state=_model_checkpoint_state(args),
                periodic_ticks=_periodic_tick_checkpoint_state(p, t_cursor), gnc_scopes=gnc_scopes,
                results_rows=rows_entry, interlink_history=history_entry, maneuver_buffers=_maneuver_buffer_checkpoint_state(p),
                orbits=_orbit_counter_checkpoint_state(p), entries=_entry_counter_checkpoint_state(p))
            checkpoint_written = true
            if string(seg_sol.retcode) != "Success" || !_gravity_backbone_time_reached(t_cursor, t_next)
                break
            end
        end

    elseif t_start < mission_end
        p.shared_buffers.solve_segment_end_time[] = mission_end
        sol = _solve_run_interval!(solver_trace, solve_setup, u_start, (t_start, mission_end), jac_prototype, results_inputs)
        last_sol = sol
        if solver_mode == :gravity_backbone_split
            _append_backbone_saved_segment!(backbone_saved_times, backbone_saved_data, sol, save_fields_resolved, p)
        end
        _require_successful_solve(sol, results_inputs, "Solve")
    end

    # Time the solve, so the next process with this signature knows whether the
    # ~0.1 s sweep is worth paying. Placed here rather than at the end of the
    # function so results/telemetry work is not counted as solve cost.
    _rhs_calib_record_solve_time!()
    if !isempty(p.shared_buffers.rhs_region_tuners) &&
       lowercase(strip(_engine_env_get("SPACEAGORA_RHS_REGION_WIDTH_TRACE", "0"))) in ("1", "true", "yes", "on")
        println("[SpaceAGORA] RHS region widths: " * join(
            ("$(k)=$(t.width)/$(t.allotment)(best $(t.best_width) $(round(t.best_ns / 1e3; digits=1))us, " *
             "full $(round(t.full_ns / 1e3; digits=1))us$(t.searching ? ", searching" : ""))"
             for (k, t) in sort!(collect(p.shared_buffers.rhs_region_tuners); by=first)), " "))
    end
    sbt = p.shared_buffers.spin_barrier_tuner
    if sbt.enabled &&
       lowercase(strip(_engine_env_get("SPACEAGORA_RHS_REGION_WIDTH_TRACE", "0"))) in ("1", "true", "yes", "on")
        println("[SpaceAGORA] spin barrier auto: $(sbt.incumbent ? "on" : "off") " *
                "(probes $(sbt.probes), switches $(sbt.switches))")
    end

    # Process and save results. A resumed run that solved nothing (its checkpoint
    # is at or past the mission's end) returns its results but writes no files:
    # those already hold the run.
    solved = t_start < mission_end
    saved_results = _save_simulation_results_if_enabled!(results_inputs...; return_table=return_results,
        write_files=solved)
    scene_path = solved ? _write_visualization_scene_if_enabled!(results_args; density_params=p) : nothing
    export_path = nothing
    if solved && visualization && args.simulation_settings.results
        export_path = SimulationModel.SceneVisualization.export_visualization(results_args)
    end

    if return_solution && checkpoint_active && args.simulation_settings.checkpoint_interval_s < mission_end
        @warn "return_solution=true with checkpointed integration returns the final segment ODESolution, not a stitched full-history ODESolution."
    end

    if return_solver_metadata
        parallel_policy = try
            SimulationModel.ParallelPolicy.policy_telemetry_snapshot()
        catch
            nothing
        end
        # `solution` stays first so callers that read result.solution keep
        # working; it is `nothing` in the metadata-only case, which is the point
        # of that case. `retcode` is new, and is what a caller that wanted only
        # the outcome should read instead of reaching through the solution.
        #
        # Two things it does not do. It does not report failure: an unsuccessful
        # retcode throws above, before this tuple is built, so the values
        # reachable here are the successful ones and a callback termination. And
        # it is `nothing` when no solve ran at all, which happens when a resumed
        # checkpoint is at or past mission end — a caller comparing it against
        # "Success" reads that as failure, so test the mission, not the field,
        # if resuming is in play. Note also that appending `retcode` here shifts
        # the position of every later field, which matters only to a caller
        # destructuring the tuple positionally; none in this repository does.
        return (
            solution=return_solution ? last_sol : nothing,
            retcode=last_sol === nothing ? nothing : string(last_sol.retcode),
            solver_mode=string(solver_mode),
            solver_trace=solver_trace,
            parallel_policy=parallel_policy,
            spice_counters=_spice_runtime_counters_snapshot(p)
        )
    end
    if return_results
        # Files this run left behind: only paths that exist, so a disabled
        # output never shows up. The bundle and checkpoint paths are the
        # engine's own naming (IOConfig), not guesses.
        prefix = _results_bundle_prefix(results_args)
        long_layout = args.simulation_settings.results_layout === :long
        candidates = Any[saved_results.csv_path, scene_path, export_path]
        # The link schedule this run wrote, beside its results.
        if saved_results.schedule_files !== nothing
            push!(candidates, saved_results.schedule_files.csv, saved_results.schedule_files.feather)
        end
        long_layout && push!(candidates, saved_results.run_csv_path)
        # The bundle is this run's only when this run wrote it; a leftover
        # from an earlier run in the same directory must not be listed.
        if solved && args.simulation_settings.results && _typed_save_bundle_enabled()
            if long_layout
                push!(candidates, prefix * "_spacecraft.feather", prefix * ".manifest.toml")
                saved_results.run_table === nothing || push!(candidates, prefix * "_run.feather")
            else
                push!(candidates, prefix * ".feather", prefix * ".manifest.toml")
            end
        end
        if checkpoint_written
            ck = SimulationModel.IOConfig._checkpoint_paths(args)
            push!(candidates, ck.data, ck.manifest, ck.rows, ck.interlinks)
            rows_log !== nothing && rows_log.stream !== nothing && push!(candidates, ck.stream_spacecraft, ck.stream_run)
        end
        files = String[string(f) for f in candidates if f !== nothing && isfile(f)]
        # A resume without usable rows writes beside the earlier results; its
        # configuration names the folder it wrote to.
        return SimulationResults(saved_results.table, results_args, files, return_solution ? last_sol : nothing,
            long_layout ? saved_results.run_table : nothing)
    end
    return_solution && return last_sol
    return nothing
    end
    finally
        # A GRAM model this run realigned or deep-copied owns a native
        # atmosphere that only the collector's finalizer frees, and the
        # collector cannot see its size. Collect once enough of them have been
        # discarded (see collect_unreferenced_gram_atmospheres!). Changes when
        # native memory is released, not what the run computed.
        SimulationModel.EnvironmentModels.collect_unreferenced_gram_atmospheres!()
    end
end
