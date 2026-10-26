# Run-start measurement of how a workload scales, for the native-lock width cap.
#
# Lives here rather than in ParallelCost because measuring means running
# `spacecraft_dynamics!` on live `ODEParams`, which sits above SimulationModel.
# ParallelCost owns the result type and its readers; this file owns the driver.

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

WARM PASSES ONLY. A first pass in a process measures compilation, and in this
codebase that has been worth three orders of magnitude rather than a few
percent, so one pass is run and discarded before any reading.

`alloc_bytes_per_pass` comes from `Base.gc_bytes()`, which is process-wide. On a
probe running alone -- which is when this is called, at run start, before any
split exists -- that is the pass's own allocation. Called concurrently with
other work it would over-report, so it is not a hot-path measurement and must
not become one.

Returns `valid = false` on any failure, which a caller must read as "no
scaling evidence" rather than "scales perfectly".
"""
function probe_contention_inputs!(
    p,
    u0,
    args;
    k::Int = 5,
    t_probe::Float64 = 0.0,
    t_step::Float64 = 0.0,
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
    error("RuntimeServices not found in module ancestry for contention_probe.jl")
end
