# One-shot per-machine measurement of the parallel cost constants.
#
# Runs once per machine, not per solve, which is what lets it be careful: every
# reading is a minimum over repeated windows. It measures what the predictive
# campaign planner reads: the serial fraction of the dispatch primitive
# (`usl_alpha_base`) and the extra contention that allocation adds
# (`usl_beta_alloc`), both from speedup ladders over the synthetic kernels in
# calibration_kernels.jl.

"""
    machine_fingerprint() -> String

Identity the constants were measured against.

Includes the effective CPU quota, not just `Sys.CPU_THREADS`. In a container or
under a cgroup the thread count reports the host's cores while the process is
allowed a fraction of them, so a fingerprint without the quota would happily
reuse a 64-core calibration inside a 4-core cgroup and predict throughput the
process cannot reach. Julia's version is included because codegen changes move
the arithmetic rates.
"""
function machine_fingerprint()::String
    override = strip(get(ENV, "SPACEAGORA_PERF_MACHINE_LABEL", ""))
    cpu_info = Sys.cpu_info()
    cpu_model = isempty(cpu_info) ? "unknown" : String(cpu_info[1].model)
    raw = string(
        isempty(override) ? cpu_model : override, "|",
        Sys.CPU_THREADS, "|",
        _cgroup_cpu_quota(), "|",
        VERSION,
    )
    return bytes2hex(SHA.sha256(codeunits(raw))[1:8])
end

# Effective CPU quota as a float number of cores, or -1.0 when unlimited or
# unreadable.
#
# Delegates to ParallelProfiles.cgroup_cpu_quota rather than reading the cgroup
# files again. This file had the only copy until the routing layer needed the
# same number to size a core budget, and a container's quota described two ways
# in one process is a divergence waiting to happen -- the fingerprint below
# depends on this value, so the two copies drifting would silently invalidate
# cached constants on one path and not the other. Value semantics are unchanged
# (cgroup v2 first, then v1, -1.0 for unlimited or unreadable), so fingerprints
# computed before this delegation still match.
function _cgroup_cpu_quota()::Float64
    mod = @__MODULE__
    while true
        if isdefined(mod, :ParallelProfiles)
            return getproperty(mod, :ParallelProfiles).cgroup_cpu_quota()
        end
        parent = parentmodule(mod)
        parent === mod && break
        mod = parent
    end
    error("ParallelProfiles.cgroup_cpu_quota not found in module ancestry for machine_calibration.jl")
end

"""
    calibrate_machine(; k = 15, verbose = false) -> MachineConstants

Measure this machine's parallel cost constants: the USL terms the predictive
campaign planner reads.
"""
function calibrate_machine(; k::Int = 15, verbose::Bool = false)::MachineConstants
    verbose && println("[calibrate] fingerprint = $(machine_fingerprint())")
    usl = _calibrate_usl(k, verbose)
    return MachineConstants(
        usl_alpha_base = usl.alpha_base,
        usl_beta_alloc = usl.beta_alloc,
        fingerprint = machine_fingerprint(),
        schema_version = CALIBRATION_SCHEMA_VERSION,
        measured_threads = Threads.nthreads(),
    )
end


# ── USL contention parameters ─────────────────────────────────────────────────

"""
    _fit_usl(workers, speedups) -> (alpha, beta)

Least-squares fit of `S(k) = k / (1 + α(k−1) + βk(k−1))` to a measured ladder.

Fitted in the linearised deficiency form rather than by nonlinear optimisation.
Rearranging gives `k/S(k) − 1 = α(k−1) + βk(k−1)`, and dividing by `(k−1)` makes
it a straight line in `k`:

    (k/S − 1)/(k − 1)  =  α  +  β·k

so `α` is the intercept and `β` the slope of an ordinary regression over the
rungs with `k > 1`. No solver, no starting guess, no convergence to check --
which matters because this runs inside calibration, where a fit that fails to
converge would have to be silently defaulted, and a silent default here reads as
"scales perfectly".

Both parameters are clamped at zero. A ladder can produce a small negative slope
from measurement noise, and a negative `β` would predict speedup rising without
bound.
"""
function _fit_usl(workers::Vector{Int}, speedups::Vector{Float64})::NTuple{2, Float64}
    xs = Float64[]
    ys = Float64[]
    for (w, s) in zip(workers, speedups)
        w > 1 || continue
        s > 0.0 || continue
        kf = Float64(w)
        deficiency = kf / s - 1.0
        push!(xs, kf)
        push!(ys, deficiency / (kf - 1.0))
    end
    isempty(xs) && return (0.0, 0.0)
    if length(xs) == 1
        # One rung cannot separate a serial fraction from a contention term.
        # Attribute it entirely to alpha, which is the conservative choice: it
        # saturates rather than inverting, so it never predicts a turning point
        # the ladder did not show.
        return (max(0.0, ys[1]), 0.0)
    end
    n = Float64(length(xs))
    sx = sum(xs); sy = sum(ys)
    sxx = sum(x -> x * x, xs); sxy = sum(i -> xs[i] * ys[i], eachindex(xs))
    denom = n * sxx - sx * sx
    if abs(denom) < 1e-12
        return (max(0.0, sy / n), 0.0)
    end
    beta = (n * sxy - sx * sy) / denom
    alpha = (sy - beta * sx) / n
    return (clamp(alpha, 0.0, 1.0), max(0.0, beta))
end

# Speedup ladder over an arbitrary per-worker body, shared by the two USL
# ladders so they differ only in what the workers do.
#
# FIXED TOTAL WORK, SPLIT W WAYS -- the one thing a speedup ladder must get
# right. The first version of this gave every worker the whole body, so `w`
# workers did `w` times the work in roughly the time one worker takes and the
# ladder reported a speedup of about one at every rung. Fitted, that came out as
# alpha = 0.83: a machine that cannot parallelise arithmetic at all. The body
# therefore receives the worker count and is responsible for taking `1/w` of
# the total.
function _speedup_ladder(body!, k::Int, budget::Int)
    workers = Int[]
    speedups = Float64[]
    serial = timed_min(k = k) do
        body!(1, 1)
    end
    serial > 0.0 || return (Int[], Float64[])
    for w in unique(Int[1, 2, 4, 8, budget])
        (w > budget || w < 1) && continue
        t = timed_min(k = k) do
            ParallelPolicy.threaded_foreach_worker_persistent(:cost_calib_usl, w, w) do wid, _i
                body!(wid, w)
            end
        end
        t > 0.0 || continue
        push!(workers, w)
        push!(speedups, max(1e-6, serial / t))
    end
    return (workers, speedups)
end

# The two USL terms the predictive planner reads, from two ladders over the same
# dispatch primitive.
#
# `alpha_base` is the arithmetic ladder's serial fraction: the dispatch
# primitive's own residual, present in every workload. The allocation term is
# the EXCESS beta an allocation ladder shows over the arithmetic one, because an
# allocation ladder alone would fold in whatever contention the dispatch
# primitive itself carries; taking the difference leaves the part due to the
# allocator. It is then expressed per GiB/s of allocation per worker, so a
# workload that allocates at half the calibration rate is predicted to contend
# half as much. The arithmetic ladder's own beta is only that baseline and is not
# saved.
function _calibrate_usl(k::Int, verbose::Bool)
    budget = max(1, Threads.nthreads())

    # 1. Arithmetic baseline: per-worker buffers, no allocation, cache-resident.
    tab = calib_table(52, 64)
    outs = [zeros(Float64, 2048, 16) for _ in 1:budget]
    total_touch = 4096
    arith_workers, arith_speedups = _speedup_ladder(k, budget) do wid, w
        calib_batch_kernel!(outs[min(wid, length(outs))], tab, 26, cld(total_touch, w), 4)
    end
    alpha_base, beta_arith = _fit_usl(arith_workers, arith_speedups)
    verbose && println("[calibrate] usl arith   ladder=", collect(zip(arith_workers,
                       round.(arith_speedups, digits = 2))))
    verbose && println("[calibrate] usl arith   alpha=$(round(alpha_base, digits=4)) " *
                       "beta=$(round(beta_arith, digits=6))")

    # 2. Allocation ladder.
    #
    # MANY SMALL SHORT-LIVED OBJECTS, not a few large ones, because that is the
    # shape the RHS actually produces. Iterating a heterogeneous effector tuple
    # whose element type is abstract boxes 144-240 bytes per satellite-effector
    # pair, immediately garbage; a first version of this ladder allocated 32 KiB
    # arrays instead and measured no contention at all (beta_total below the
    # arithmetic baseline at 23 GiB/s per worker), because large blocks and
    # short-lived boxes exercise different parts of the allocator.
    #
    # The rate one worker demands is measured here rather than assumed, so the
    # coefficient is per GiB/s and a workload allocating at a different rate
    # scales the prediction rather than being assumed to match.
    total_arrays, alloc_elems = 400_000, 4
    alloc_bytes = Float64(total_arrays * alloc_elems * 8)
    alloc_serial_ns = timed_min(k = k) do
        calib_alloc_kernel!(total_arrays, alloc_elems)
    end
    alloc_workers, alloc_speedups = _speedup_ladder(k, budget) do _wid, w
        calib_alloc_kernel!(cld(total_arrays, w), alloc_elems)
    end
    _, beta_alloc_total = _fit_usl(alloc_workers, alloc_speedups)
    alloc_gibs = alloc_serial_ns > 0.0 ?
        (alloc_bytes * 1.0e9 / alloc_serial_ns) / 1024^3 : 0.0
    beta_alloc = alloc_gibs > 0.0 ?
        max(0.0, beta_alloc_total - beta_arith) / alloc_gibs : 0.0
    verbose && println("[calibrate] usl alloc   ladder=", collect(zip(alloc_workers,
                       round.(alloc_speedups, digits = 2))))
    verbose && println("[calibrate] usl alloc   beta_total=$(round(beta_alloc_total, digits=6)) " *
                       "at $(round(alloc_gibs, digits=2)) GiB/s/worker -> " *
                       "coeff=$(round(beta_alloc, digits=6))")

    return (alpha_base = alpha_base, beta_arith = beta_arith, beta_alloc = beta_alloc)
end

"""
    machine_constants_path() -> String

Where this machine's constants live. `SPACEAGORA_COST_CONSTANTS_PATH` overrides;
otherwise a fingerprinted file under `output/parallel_policy_state/`, alongside
the RHS calibration and inner-policy state.
"""
function machine_constants_path()::String
    override = strip(get(ENV, "SPACEAGORA_COST_CONSTANTS_PATH", ""))
    if !isempty(override)
        return normpath(isabspath(override) ? override : joinpath(pwd(), override))
    end
    return normpath(joinpath(
        pwd(), "output", "parallel_policy_state",
        "cost_constants_$(machine_fingerprint()).toml",
    ))
end

function save_machine_constants(mc::MachineConstants, path::AbstractString = machine_constants_path())::String
    payload = Dict{String, Any}(
        "schema_version" => mc.schema_version,
        "fingerprint" => mc.fingerprint,
        "measured_threads" => mc.measured_threads,
        "usl_alpha_base" => mc.usl_alpha_base,
        "usl_beta_alloc" => mc.usl_beta_alloc,
    )
    path_s = String(path)
    # The `[campaign]` table holds costs measured from campaign traces, not by
    # this calibration (the predictive planner's pool terms; see
    # scripts/calibration/extract_campaign_cost_terms.py). Re-calibrating must not erase
    # them, so an existing table is carried across the rewrite.
    if isfile(path_s)
        previous = try
            TOML.parsefile(path_s)
        catch
            nothing
        end
        if previous !== nothing && get(previous, "campaign", nothing) isa AbstractDict
            payload["campaign"] = previous["campaign"]
        end
    end
    RuntimeServices.write_file_atomically(io -> TOML.print(io, payload), path_s)
    Threads.atomic_add!(_MACHINE_CONSTANTS_GENERATION, 1)
    return path_s
end

# Bumped by every `save_machine_constants` in this process, so a consumer that
# caches the loaded constants (the predictive planner does) can tell that the
# file it cached -- possibly as "absent" -- has been written since.
const _MACHINE_CONSTANTS_GENERATION = Threads.Atomic{Int}(0)

"""
    machine_constants_generation() -> Int

How many times this process has written machine constants. Caches of
[`load_machine_constants`](@ref) compare it to decide whether to re-read.
"""
machine_constants_generation()::Int = _MACHINE_CONSTANTS_GENERATION[]

# Paths `ensure_machine_constants!` has already settled in this process, so a
# campaign of many parallel runs checks the file once, not once per run.
const _ENSURED_CONSTANTS_PATHS = Set{String}()
const _ENSURE_CONSTANTS_LOCK = ReentrantLock()

"""
    ensure_machine_constants!(; path = machine_constants_path(),
                              calibrate = calibrate_machine,
                              threads = Threads.nthreads(),
                              lock_wait_s = 120.0, lock_stale_s = 900.0) -> Symbol

Make sure this machine's cost constants exist, calibrating them once if not.
This is what `SolverConfig(parallel=true)` calls before its first run, so the
predictive planner never runs without machine constants merely because nobody
ran `scripts/calibration/calibrate_machine.jl`.

Returns
  - `:present` when `path` already holds current constants (nothing is
    measured);
  - `:calibrated` when it did not and `calibrate()` was run and saved there
    (the same measurement and the same fingerprinted path as the script);
  - `:checked` when this process has already settled `path`;
  - `:failed` when the calibration threw (it is logged and not retried in this
    process; the planner then runs with whatever constants the file holds, or
    none, as it did before calibration was automatic);
  - `:deferred` when another process held the calibration lock for longer than
    `lock_wait_s` (this process does not measure while that one does; it uses
    whatever the file holds when the wait ends).
The check runs once per path per process.

Constants are current when the file has the current schema and was measured
with at least `threads` Julia threads. The USL ladders only reach the thread
count of the session that measured them, so a file written by a one-thread
session is re-measured the first time a wider session runs with the flag; a file
measured wider serves a narrower session as is, and is never overwritten by it.
A file that does not record the thread count, or records 0, is measured again.

Several processes started together in one working directory (a job array,
parallel test runners) share the fingerprinted path. Only one of them measures:
the others wait on the advisory lock file `path * ".lock"` and then read what
it saved, so the saved constants are not inflated by the processes contending
with each other. A lock older than `lock_stale_s` is taken to be left by a
process that died and is removed.
"""
function ensure_machine_constants!(;
    path::AbstractString = machine_constants_path(),
    calibrate = calibrate_machine,
    threads::Integer = Threads.nthreads(),
    lock_wait_s::Real = 120.0,
    lock_stale_s::Real = 900.0,
)::Symbol
    path_s = normpath(String(path))
    lock(_ENSURE_CONSTANTS_LOCK) do
        path_s in _ENSURED_CONSTANTS_PATHS && return :checked
        _constants_current(load_machine_constants(path_s), threads) &&
            return _settle_constants_path!(path_s, :present)
        lock_path = path_s * ".lock"
        if !_acquire_calibration_lock(lock_path, lock_stale_s)
            # Another process is measuring. Measuring alongside it would bake
            # the contention into both results; wait for its file instead.
            deadline = time() + Float64(lock_wait_s)
            while isfile(lock_path) && time() < deadline
                sleep(0.25)
            end
            _constants_current(load_machine_constants(path_s), threads) &&
                return _settle_constants_path!(path_s, :present)
            _acquire_calibration_lock(lock_path, lock_stale_s) ||
                return _settle_constants_path!(path_s, :deferred)
        end
        try
            # The process that held the lock may have saved while this one
            # was deciding to take it.
            _constants_current(load_machine_constants(path_s), threads) &&
                return _settle_constants_path!(path_s, :present)
            @info "Measuring this machine's parallel cost constants once (a few seconds); later runs reuse them." path = path_s threads = Int(threads)
            push!(_ENSURED_CONSTANTS_PATHS, path_s)
            try
                save_machine_constants(calibrate(), path_s)
            catch err
                err isa InterruptException && rethrow()
                # A run the user asked for must not fail because a measurement
                # the planner can do without did; it then predicts without
                # contention, exactly as it did before calibration was automatic.
                @warn "Machine calibration failed; the parallel planner will run without machine constants." exception = (err, catch_backtrace())
                return :failed
            end
            return :calibrated
        finally
            rm(lock_path; force = true)
        end
    end
end

function _settle_constants_path!(path_s::String, outcome::Symbol)::Symbol
    push!(_ENSURED_CONSTANTS_PATHS, path_s)
    return outcome
end

_constants_current(mc::Nothing, threads::Integer)::Bool = false
_constants_current(mc::MachineConstants, threads::Integer)::Bool =
    mc.measured_threads >= max(1, threads)

# Create `lock_path` exclusively (O_CREAT|O_EXCL), holding this process's id.
# False when another process holds it; a lock older than `stale_s` seconds is
# removed first (its owner died before its `finally` ran). Any other failure
# (a read-only or missing directory) returns true: the lock is advisory, and
# the save that follows reports the real problem through the `:failed` path,
# exactly as it did before there was a lock.
function _acquire_calibration_lock(lock_path::String, stale_s::Real)::Bool
    flags = Base.Filesystem.JL_O_CREAT | Base.Filesystem.JL_O_EXCL | Base.Filesystem.JL_O_WRONLY
    for _ in 1:2
        try
            mkpath(dirname(lock_path))
            io = Base.Filesystem.open(lock_path, flags, 0o644)
            try
                write(io, string(gethostname(), " ", getpid(), "\n"))
            finally
                close(io)
            end
            return true
        catch err
            err isa InterruptException && rethrow()
            (err isa Base.IOError && err.code == Base.UV_EEXIST) || return true
            age = time() - mtime(lock_path)   # mtime is 0 once the file is gone
            if isfile(lock_path) && age <= stale_s
                return false
            end
            isfile(lock_path) && rm(lock_path; force = true)
        end
    end
    return false
end

"""
    load_machine_constants(path = machine_constants_path()) -> Union{Nothing, MachineConstants}

Read cached constants, or `nothing` if absent, unparseable, written by a
different schema version, or missing one of the two USL terms or the measuring
thread count. A stale-schema
file is discarded rather than reinterpreted: term semantics changed by
definition when the version moved, so its numbers do not mean what the current
readers would take them to mean.
"""
function load_machine_constants(path::AbstractString = machine_constants_path())::Union{Nothing, MachineConstants}
    path_s = String(path)
    isfile(path_s) || return nothing
    parsed = try
        TOML.parsefile(path_s)
    catch
        return nothing
    end
    get(parsed, "schema_version", -1) == CALIBRATION_SCHEMA_VERSION || return nothing
    try
        return MachineConstants(
            usl_alpha_base = Float64(parsed["usl_alpha_base"]),
            usl_beta_alloc = Float64(parsed["usl_beta_alloc"]),
            fingerprint = String(get(parsed, "fingerprint", "")),
            schema_version = CALIBRATION_SCHEMA_VERSION,
            measured_threads = Int(parsed["measured_threads"]),
        )
    catch
        return nothing
    end
end
