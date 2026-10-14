# Run-start measurement of per-effector cost, for the analytic cost model.
#
# Lives here rather than in ParallelCost because probing means calling
# `_evaluate_dynamic_effector`, which needs live `ODEParams`, populated prefill
# buffers, and the engine's wrench/calcForceTorque dispatch -- all of which sit
# above SimulationModel. ParallelCost owns the estimators and the result types;
# this file owns the driver that feeds them real effectors.
#
# Cost: bounded at roughly (k + 2) * target_window_ns per effector, because
# `timed_min` sizes its inner repeat count to fill a window rather than issuing
# a fixed number of calls. At the defaults that is well under a millisecond for
# a whole effector tuple -- against a pre-solve sweep that evaluated the entire
# RHS fifteen times for each of up to ten candidate plans.

"""
    probe_effector_costs!(p, u0, args; kwargs...) -> Dict{Int, EffectorProbe}

Measure each dynamic effector's per-satellite cost against live state, keyed by
its index in `args.dynamics_model.dynamic_effectors`.

This is what makes the cost model work for effectors it cannot count: a
user-defined effector needs no `effector_cost_terms` method, it is simply
measured. Effectors that *do* declare are probed as well unless
`validate_declared = false`, so the declaration can be checked against reality
rather than trusted; the resulting `EffectorProbe` carries `declared = true` and
`probe_ns_map` omits it from the cost path.

One full RHS evaluation runs first. That is not warm-up for its own sake --
wrench-based effectors read per-satellite environment samples out of
`shared_buffers.rhs_flat_state_*`, which only the RHS's prefill phase populates,
and n-body/solar effectors need `_prefill_shared_body_samples!` to have warmed
the SPICE memo. Probing before that would measure a cold path the solve never
takes, or throw.

Two distinct failure modes, both non-fatal:

  - The warm-up RHS throws. Nothing can be probed, so an empty table is
    returned. A solve whose RHS cannot be evaluated once was not going to run
    anyway, so this costs nothing; the model simply has no data.
  - An individual effector throws while being probed in isolation, having
    survived the warm-up -- a stateful effector that fails on re-evaluation, or
    one that needs a partition context the isolated call does not supply. It is
    recorded with `ns_per_sat = NaN`.

Either way the model degrades to "cannot predict this workload" and the caller
keeps its heuristic route. `constellation_work_counts` turns an unknown
effector with no usable measurement into `in_domain = false`, so a failed probe
cannot silently become a zero-cost effector.

`t_probe` defaults to the mission start. An effector whose cost varies strongly
with state -- aerodynamics in vacuum versus in the atmosphere -- is measured in
whatever regime that corresponds to, which is a real limitation of any
start-of-run measurement. The observed-cost EMA already maintained in
`_update_rhs_effector_cost_model!` is what catches divergence later.
"""
function probe_effector_costs!(
    p,
    u0,
    args;
    k::Int = 7,
    validate_declared::Bool = true,
    t_probe::Float64 = 0.0,
)::Dict{Int, SimulationModel.ParallelCost.EffectorProbe}
    PC = SimulationModel.ParallelCost
    probes = Dict{Int, PC.EffectorProbe}()

    effectors = args.dynamics_model.dynamic_effectors
    isempty(effectors) && return probes

    sat_idx = findfirst(identity, p.is_active)
    sat_idx === nothing && return probes

    # Populate prefill buffers and warm the SPICE memo; see docstring.
    try
        du = zero(u0)
        spacecraft_dynamics!(du, u0, p, t_probe)
    catch err
        @debug "Cost probe: RHS warm-up failed; skipping effector probes." exception=err
        return probes
    end

    sc_state = u0.sc
    spacecraft = args.dynamics_model.spacecraft
    orientation_sim = args.mission_configuration.orientation_sim

    # One reference reading for the whole probe pass, taken alongside the
    # measurements so contention is common-mode between them.
    reference_ns = PC.reference_kernel_ns()

    @inbounds for idx in eachindex(effectors)
        effector = effectors[idx]
        declared = PC.effector_cost_terms(effector) !== nothing
        if declared && !validate_declared
            continue
        end

        measured = try
            @views sc_view = sc_state[sat_idx]
            state_sample = _wrench_method_available(effector) ?
                _rhs_flat_state_sample_from_buffers(
                    p.shared_buffers, spacecraft, sat_idx, orientation_sim
                ) :
                nothing
            PC.timed_min(k = k) do
                force, torque = _evaluate_dynamic_effector(
                    effector, sc_view, state_sample, p, sat_idx, t_probe
                )
                force[1] + torque[1]
            end
        catch err
            @debug "Cost probe: effector evaluation failed; leaving it unmeasured." index=idx exception=err
            NaN
        end

        ratio = (isfinite(measured) && reference_ns > 0.0) ? measured / reference_ns : NaN
        probes[idx] = PC.EffectorProbe(measured, ratio, declared)
    end

    return probes
end

"""
    probe_contention_inputs!(p, u0, args; kwargs...) -> ContentionInputs

Measure how badly this workload will scale, from warm full-RHS passes.

Three readings, taken together so they describe one state of the machine:

  - allocation bytes and elapsed time for one pass, both minimised over `k`
    repeats. The minimum, not the mean: a pass that happened to catch a GC
    or a scheduler preemption measures the interruption, not the workload.
  - native-lock occupancy over the same window, differenced from the shared
    counters rather than reset, so a caller already measuring a wider region
    is not clobbered.

WARM PASSES ONLY. The caller must have run at least one RHS evaluation before
this, which `probe_effector_costs!` already does for its own reasons; a first
pass measures compilation, and in this codebase that has been worth three orders
of magnitude rather than a few percent.

`alloc_bytes_per_pass` comes from `Base.gc_bytes()`, which is process-wide. On a
probe running alone -- which is when this is called, at run start, before any
split exists -- that is the pass's own allocation. Called concurrently with
other work it would over-report, so it is not a hot-path measurement and must
not become one.

Returns `valid = false` on any failure, which the predictor must read as "no
scaling evidence" rather than "scales perfectly".
"""
function probe_contention_inputs!(
    p,
    u0,
    args;
    k::Int = 5,
    t_probe::Float64 = 0.0,
    t_step::Float64 = 0.0,
    workspace_bytes_per_sat::Real = 0.0,
)::SimulationModel.ParallelCost.ContentionInputs
    PC = SimulationModel.ParallelCost
    RS = SpaceAGORA_runtime_services()

    n_active = count(identity, p.is_active)

    du = try
        zero(u0)
    catch err
        @debug "Contention probe: could not allocate a derivative buffer." exception=err
        return PC.ContentionInputs()
    end

    # Warm-up pass, discarded. See the docstring: the first evaluation in a
    # process measures compilation.
    try
        spacecraft_dynamics!(du, u0, p, t_probe)
    catch err
        @debug "Contention probe: warm-up RHS failed; no scaling evidence." exception=err
        return PC.ContentionInputs()
    end

    before = RS.native_lock_stats_snapshot()
    best_ns = Inf
    best_bytes = Inf
    for rep in 1:max(1, k)
        # `t_step` exists because a fixed-epoch probe systematically
        # under-measures lock occupancy. Third-body ephemeris samples are
        # memoised on the epoch, so repeated passes at one `t` hit the memo and
        # take the native lock exactly never -- measured rho = 0 on a workload
        # whose solve holds the lock on every step. Stepping the epoch makes the
        # memo miss the way an advancing solve does. Zero keeps the
        # single-epoch behaviour for callers that want it.
        t = t_probe + Float64(rep - 1) * t_step
        bytes0 = Base.gc_bytes()
        t0 = time_ns()
        ok = try
            spacecraft_dynamics!(du, u0, p, t)
            true
        catch err
            @debug "Contention probe: RHS pass failed mid-measurement." exception=err
            false
        end
        ok || return PC.ContentionInputs()
        elapsed = Float64(time_ns() - t0)
        allocated = Float64(Base.gc_bytes() - bytes0)
        best_ns = min(best_ns, elapsed)
        # Allocation is deterministic per pass where it is real, so the minimum
        # is the pass's own figure with any GC bookkeeping that landed in a
        # particular repeat excluded.
        allocated >= 0.0 && (best_bytes = min(best_bytes, allocated))
    end
    after = RS.native_lock_stats_snapshot()

    isfinite(best_ns) || return PC.ContentionInputs()
    isfinite(best_bytes) || (best_bytes = 0.0)

    # Occupancy is charged against the whole measured window, so divide the
    # differenced totals by the repeat count to put them on a per-pass footing
    # with `pass_ns`.
    reps = Float64(max(1, k))
    return PC.ContentionInputs(
        alloc_bytes_per_pass = best_bytes,
        pass_ns = best_ns,
        workspace_bytes_per_sat = Float64(workspace_bytes_per_sat),
        n_active_sats = n_active,
        lock_acquisitions = after.acquisitions - before.acquisitions,
        lock_hold_ns = max(0.0, Float64(after.hold_ns - before.hold_ns)) / reps,
        lock_wait_ns = max(0.0, Float64(after.wait_ns - before.wait_ns)) / reps,
        valid = true,
    )
end

# RuntimeServices sits above SimulationModel in the include order, so it is
# reached through the module ancestry rather than imported -- the same pattern
# perturbations.jl and reference_system.jl use for the lock itself.
@inline function SpaceAGORA_runtime_services()
    mod = @__MODULE__
    while true
        if isdefined(mod, :RuntimeServices)
            return getproperty(mod, :RuntimeServices)
        end
        parent = parentmodule(mod)
        parent === mod && break
        mod = parent
    end
    error("RuntimeServices not found in module ancestry for rhs_cost_probe.jl")
end

"""
    warn_on_declaration_mismatch(probes, predicted_ns; factor = 2.0) -> Int

Compare declared effector costs against their probes and warn on gross
disagreement, returning the number of mismatches.

`predicted_ns` maps effector index to the cost the declaration implies under the
current machine constants. Only indices present in both, with a finite probe,
are checked. A declaration inside `factor` is accepted -- the model consumes
rankings, not absolute times, so a factor-of-two error in one effector among
several rarely flips a candidate ordering, while an order-of-magnitude error
usually does.

Warning rather than erroring is deliberate: a wrong declaration should be
visible and should cost the user the model's accuracy, not their run.
"""
function warn_on_declaration_mismatch(
    probes::Dict{Int, SimulationModel.ParallelCost.EffectorProbe},
    predicted_ns::Dict{Int, Float64};
    factor::Float64 = 2.0,
)::Int
    PC = SimulationModel.ParallelCost
    mismatches = 0
    for (idx, probe) in probes
        probe.declared || continue
        isfinite(probe.ns_per_sat) || continue
        predicted = get(predicted_ns, idx, nothing)
        predicted === nothing && continue
        if !PC.validate_declaration(probe.ns_per_sat, predicted; factor=factor)
            mismatches += 1
            @warn(
                "effector_cost_terms declares a cost that disagrees with measurement; " *
                "the cost model will route this workload less accurately. " *
                "Either correct the declaration or delete it and let the probe measure it.",
                effector_index = idx,
                measured_ns_per_sat = probe.ns_per_sat,
                declared_ns_per_sat = predicted,
                ratio = probe.ns_per_sat / predicted,
            )
        end
    end
    return mismatches
end
