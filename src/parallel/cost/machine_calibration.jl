# One-shot per-machine measurement of the cost model's rate constants.
#
# Runs once per machine, not per solve, which is what lets it be careful: every
# reading is a minimum over repeated windows, and the whole pass costs a few
# seconds against a per-solve sweep it removes entirely.
#
# The design matrix is synthetic on purpose. Measured across every configuration
# the real harmonics kernel can produce, `simd_terms` and `coeff_touches` stay
# 0.998 correlated, so a fit against real workloads cannot separate them -- and
# the separation is exactly what distinguishes the routing candidates. The
# synthetic kernels take touch count and arithmetic-per-touch as independent
# arguments, so the two columns can be made orthogonal by construction.
#
# The real kernel is then held out entirely and used to *validate* the result by
# cold prediction. Fitting on it would make the validation circular; holding it
# out makes a good prediction evidence that the synthetic constants transfer,
# and a bad one a localised, diagnosable residual.

# 4: adds the USL contention parameters (usl_alpha_base, usl_beta_arith,
# usl_beta_alloc, usl_beta_bw, llc_bytes). Bumped rather than defaulted because
# a file without them would silently predict zero contention -- which reads as
# "this machine scales perfectly" and is the most dangerous possible default.
const CALIBRATION_SCHEMA_VERSION = 4

# Footprint ladder for the coefficient-touch curve: square-ish tables, matching
# the (L+2) x (M+2) shape of a full-field harmonics model, spanning L1-resident
# to well past L3.
const _COEFF_KNOTS = ((22, 20), (52, 50), (102, 100), (202, 200),
                      (402, 400), (802, 800), (1602, 1600))

# (batch, depth) pairs for the SIMD-lane curve, which is indexed by WORKSPACE
# FOOTPRINT rather than by batch width alone.
#
# Both dimensions are swept because the real workspace is
# batch x (L+3) x (M+2): at batch 1024 with L=20 that is 4.1 MB and nowhere near
# cache, while a bare 1024-element vector is 8 KB and L1-resident. Calibrating
# against the vector measured arithmetic throughput and missed the memory
# traffic that dominates the real kernel at large batch. A wide-shallow and a
# narrow-deep workspace reaching the same footprint cost the same, and only
# footprint predicts that.
#
# The ladder still reaches batch 1: `satellite_batch` processes one satellite at
# a time and lives entirely at that end, where per-pass overhead cannot amortize.
const _LANE_KNOTS = ((1, 8), (4, 8), (16, 8), (64, 8), (256, 8),
                     (256, 64), (1024, 32), (1024, 128), (4096, 64),
                     (4096, 256), (16384, 128), (16384, 512))

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

Measure every rate constant on this machine.

Order matters in one place: the coefficient-touch curve is measured first so its
contribution can be subtracted from the SIMD-lane readings, which are taken with
the batch kernel and therefore include a small number of touches.
"""
function calibrate_machine(; k::Int = 15, verbose::Bool = false)::MachineConstants
    verbose && println("[calibrate] fingerprint = $(machine_fingerprint())")

    coeff_curve = _calibrate_coeff_touch(k, verbose)
    lane_curve = _calibrate_simd_lane(k, verbose, coeff_curve)
    scalar_ns = _calibrate_scalar_item(k, verbose)
    node_ns, pool_base, pool_per_worker, atomic_ns = _calibrate_dispatch(k, verbose)
    batch_base, batch_per_worker = _calibrate_polyester_dispatch(k, verbose)
    speedup_curve = _calibrate_parallel_speedup(k, verbose)
    usl = _calibrate_usl(k, verbose)

    return MachineConstants(
        simd_lane = lane_curve,
        coeff_touch = coeff_curve,
        parallel_speedup = speedup_curve,
        ns_per_scalar_item = scalar_ns,
        ns_per_queue_node = node_ns,
        dispatch_pool_ns_base = pool_base,
        dispatch_pool_ns_per_worker = pool_per_worker,
        dispatch_batch_ns_base = batch_base,
        dispatch_batch_ns_per_worker = batch_per_worker,
        ns_per_atomic = atomic_ns,
        reference_fma_ns = reference_kernel_ns(k = k),
        reference_mem_ns = reference_memory_kernel_ns(k = k),
        usl_alpha_base = usl.alpha_base,
        usl_beta_arith = usl.beta_arith,
        usl_beta_alloc = usl.beta_alloc,
        usl_beta_bw = usl.beta_bw,
        llc_bytes = usl.llc_bytes,
        fingerprint = machine_fingerprint(),
        schema_version = CALIBRATION_SCHEMA_VERSION,
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

# Last-level cache size in bytes, for normalising a working set.
#
# Read from sysfs rather than assumed: the footprint term is a ratio against
# this, so a wrong value rescales every bandwidth prediction on the machine. A
# conservative 8 MiB fallback is used where sysfs is unavailable, and callers
# treat `0.0` as "unknown" and drop the term entirely.
function _llc_bytes()::Float64
    best = 0.0
    base = "/sys/devices/system/cpu/cpu0/cache"
    isdir(base) || return 8.0 * 1024^2
    try
        for entry in readdir(base)
            startswith(entry, "index") || continue
            path = joinpath(base, entry, "size")
            isfile(path) || continue
            raw = strip(read(path, String))
            m = match(r"^(\d+)([KMG]?)$", raw)
            m === nothing && continue
            value = parse(Float64, m.captures[1])
            scale = m.captures[2] == "K" ? 1024.0 :
                    m.captures[2] == "M" ? 1024.0^2 :
                    m.captures[2] == "G" ? 1024.0^3 : 1.0
            best = max(best, value * scale)
        end
    catch
        return 8.0 * 1024^2
    end
    return best > 0.0 ? best : 8.0 * 1024^2
end

# Speedup ladder over an arbitrary per-worker body, shared by the three
# contention calibrations so they differ only in what the workers do.
#
# FIXED TOTAL WORK, SPLIT W WAYS -- the same shape `_calibrate_parallel_speedup`
# uses, and the one thing a speedup ladder must get right. The first version of
# this gave every worker the whole body, so `w` workers did `w` times the work
# in roughly the time one worker takes and the ladder reported a speedup of
# about one at every rung. Fitted, that came out as alpha = 0.83: a machine that
# cannot parallelise arithmetic at all. The body therefore receives the worker
# count and is responsible for taking `1/w` of the total.
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

# Contention coefficients, measured as the EXCESS beta a resource-bound ladder
# shows over an arithmetic one.
#
# Measuring beta from an allocation ladder alone would fold in whatever
# contention the dispatch primitive itself carries, and that is already present
# in every workload; taking the difference leaves the part attributable to the
# resource. The coefficients are then per unit of the workload quantity the
# predictor can count -- GiB/s of allocation per worker, and working set per
# worker as a multiple of the last-level cache -- so a workload that allocates
# at half the calibration rate is predicted to contend half as much.
function _calibrate_usl(k::Int, verbose::Bool)
    budget = max(1, Threads.nthreads())
    llc = _llc_bytes()

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

    # 3. Bandwidth ladder. Each worker streams its OWN buffer, so the per-worker
    #    footprint -- the axis the coefficient is expressed in -- is held fixed
    #    while the aggregate grows past the cache and the workers contend for
    #    DRAM. Splitting one shared buffer instead would shrink the footprint
    #    per worker as the ladder widens and confound the two effects.
    stream_elems = max(1024, Int(ceil(0.5 * llc / 8)))
    buffers = [fill(1.0, stream_elems) for _ in 1:budget]
    total_passes = 16
    stream_workers, stream_speedups = _speedup_ladder(k, budget) do wid, w
        calib_stream_kernel!(buffers[min(wid, length(buffers))], cld(total_passes, w))
    end
    _, beta_bw_total = _fit_usl(stream_workers, stream_speedups)
    footprint_ratio = llc > 0.0 ? (stream_elems * 8.0) / llc : 0.0
    beta_bw = footprint_ratio > 0.0 ?
        max(0.0, beta_bw_total - beta_arith) / footprint_ratio : 0.0
    verbose && println("[calibrate] usl bw      ladder=", collect(zip(stream_workers,
                       round.(stream_speedups, digits = 2))))
    verbose && println("[calibrate] usl bw      beta_total=$(round(beta_bw_total, digits=6)) " *
                       "at $(round(footprint_ratio, digits=2)) x LLC -> " *
                       "coeff=$(round(beta_bw, digits=6))")

    return (alpha_base = alpha_base, beta_arith = beta_arith,
            beta_alloc = beta_alloc, beta_bw = beta_bw, llc_bytes = llc)
end

function _calibrate_coeff_touch(k::Int, verbose::Bool)::RateCurve
    log2_size = Float64[]
    ns = Float64[]
    for (rows, cols) in _COEFF_KNOTS
        tab = calib_table(rows, cols)
        row = max(1, rows ÷ 2)
        bytes = rows * cols * 8
        per_touch = timed_min(() -> calib_touch_kernel!(tab, row, cols); k = k) / cols
        push!(log2_size, log2(Float64(bytes)))
        push!(ns, per_touch)
        verbose && println("[calibrate] coeff_touch  $(rows)x$(cols)  " *
                           "$(round(bytes/1024, digits=1)) KiB -> $(round(per_touch, digits=4)) ns")
    end
    return RateCurve(log2_size, ns)
end

function _calibrate_simd_lane(k::Int, verbose::Bool, coeff_curve::RateCurve)::RateCurve
    # Fixed, small touch count so the subtraction below stays a correction
    # rather than a substantial part of the reading.
    n_touch, spt = 16, 4
    rows, cols = 52, 64
    tab = calib_table(rows, cols)
    touch_ns = rate_at(coeff_curve, rows * cols * 8) * n_touch

    pairs = Tuple{Float64, Float64}[]
    for (B, depth) in _LANE_KNOTS
        ws = zeros(Float64, B, depth)
        total = timed_min(() -> calib_batch_kernel!(ws, tab, rows ÷ 2, n_touch, spt); k = k)
        # The synthetic pass is one muladd per lane = 2 flops, and simd_terms is
        # counted in flops, so the rate must be per flop for the two to compose.
        lane_terms = Float64(n_touch * spt * B) * 2.0
        per_lane = max(0.0, total - touch_ns) / lane_terms
        bytes = Float64(B * depth * 8)
        push!(pairs, (log2(bytes), per_lane))
        verbose && println("[calibrate] simd_lane    B=$(B) d=$(depth) " *
                           "$(round(bytes/1024, digits=1)) KiB -> $(round(per_lane, digits=6)) ns/flop-lane")
    end
    # Several (batch, depth) pairs can land on the same footprint; keep the
    # cheapest reading at each, since the minimum is the uncontended estimate.
    sort!(pairs; by = first)
    log2_size = Float64[]
    ns = Float64[]
    for (lg, v) in pairs
        if !isempty(log2_size) && isapprox(log2_size[end], lg; atol = 1e-9)
            ns[end] = min(ns[end], v)
        else
            push!(log2_size, lg)
            push!(ns, v)
        end
    end
    return RateCurve(log2_size, ns)
end

# Achieved parallel speedup against worker count, measured rather than assumed.
#
# The model divided per-worker work by the worker count, i.e. assumed perfect
# scaling. Measured on this box with a kernel shaped like the harmonics batch,
# achieved speedup saturates near 3x at twelve workers (efficiency 0.99 at one
# worker, 0.26 at twelve), and the real harmonics kernel shows the same ~2.3x
# ceiling. Assuming linear scaling over-credits every wide plan, and no
# correction to the other constants can compensate for a term that is absent.
#
# Batch width is held constant across worker counts so this measures contention
# alone; the batch-width effect is already carried by the lane curve and would
# otherwise be counted twice.
function _calibrate_parallel_speedup(k::Int, verbose::Bool)::RateCurve
    budget = Threads.nthreads()
    tab = calib_table(52, 64)
    outs = [zeros(Float64, 2048, 16) for _ in 1:max(1, budget)]
    total_touch = 512

    serial = timed_min(k = k) do
        calib_batch_kernel!(outs[1], tab, 26, total_touch, 4)
    end

    log2_w = Float64[]
    speedup = Float64[]
    for w in unique(Int[1, 2, 3, 4, 6, 8, budget])
        w > budget && continue
        per = cld(total_touch, w)
        t = timed_min(k = k) do
            ParallelPolicy.threaded_foreach_worker_persistent(:cost_calib_speedup, w, w) do wid, _i
                calib_batch_kernel!(outs[wid], tab, 26, per, 4)
            end
            @inbounds outs[1][1, 1]
        end
        s = t > 0.0 ? serial / t : 1.0
        push!(log2_w, log2(Float64(w)))
        push!(speedup, max(1.0, s))
        verbose && println("[calibrate] speedup      w=$(w) -> $(round(s, digits=2))x " *
                           "(efficiency $(round(s / w, digits=3)))")
    end
    return RateCurve(log2_w, speedup)
end

function _calibrate_scalar_item(k::Int, verbose::Bool)::Float64
    n = 4096
    state = [1.0 + i * 1e-6 for i in 1:n]
    out = zeros(Float64, n)
    # Per flop, matching how scalar_items is counted.
    per_item = timed_min(() -> calib_scalar_kernel!(out, state, n); k = k) / (n * CALIB_SCALAR_FLOPS)
    verbose && println("[calibrate] scalar_item  -> $(round(per_item, digits=5)) ns")
    return per_item
end

# Dispatch, per-node, and atomic costs all come from the real inner-dispatch
# primitive rather than a synthetic stand-in: they are properties of that
# implementation (persistent channel pool, spawn/join, the dynamic scheduler's
# atomic counter), not of the machine in isolation, so measuring anything else
# would calibrate the wrong thing.
function _calibrate_dispatch(k::Int, verbose::Bool)
    PP = ParallelPolicy
    budget = Threads.nthreads()

    # Every loop body below writes to a preallocated buffer rather than being
    # empty. An empty body is elided, and the elided measurement is not
    # obviously wrong -- the first version of this file reported
    # ns_per_queue_node = 0.021 ns, which is *below the harness noise floor* and
    # was timing a loop that no longer existed. A body that stores makes the
    # work observable to the compiler; the cost of that one store is included in
    # ns_per_queue_node, which is correct in kind since a real queue node also
    # writes its result.
    n_items = 8192
    sink = zeros(Float64, n_items)

    # Per-node cost from the SLOPE between two item counts, so the fixed
    # dispatch cost cancels instead of being amortised into the per-item number.
    per_node = let
        t_small = timed_min(k = k) do
            PP.threaded_foreach_worker(n_items ÷ 8, 1) do _w, i
                @inbounds sink[i] = i * 1.0000001
            end
            @inbounds sink[1]
        end
        t_large = timed_min(k = k) do
            PP.threaded_foreach_worker(n_items, 1) do _w, i
                @inbounds sink[i] = i * 1.0000001
            end
            @inbounds sink[1]
        end
        max(0.0, t_large - t_small) / (n_items - n_items ÷ 8)
    end
    verbose && println("[calibrate] queue_node   -> $(round(per_node, digits=5)) ns")

    floor_ns = timing_noise_floor(k = k)
    if per_node <= 2 * floor_ns / n_items
        @warn "Calibration: per-node cost is at the measurement floor; the loop body may have been optimised away." per_node floor_ns
    end

    if budget <= 1
        verbose && println("[calibrate] dispatch/atomic skipped: single-threaded session")
        return (per_node, NaN, NaN, NaN)
    end

    # Measured against threaded_foreach_worker_persistent, which is the
    # dispatcher the flat routes actually use. The first version measured
    # Threads.@spawn instead and applied the result to every route; the
    # persistent pool is 1.7-3.7x SLOWER at the same width, and Polyester
    # (satellite_batch) is cheaper than both, so one constant was wrong in
    # magnitude and in sign-of-error depending on which route it was applied to.
    #
    # Dispatch cost against worker count: one item per worker, so the reading is
    # dominated by fan-out and join rather than by the body. Least-squares line
    # through the knots gives base and per-worker.
    # Six knots rather than three, and a larger k. Thread dispatch is the
    # highest-variance thing measured here -- it pays OS wake latency on every
    # single call, so there is no uncontended window for `timed_min` to find,
    # and the minimum is much less effective than it is on pure compute. Fitting
    # a line through three noisy points amplified that: dispatch_ns_per_worker
    # reproduced to only 53% across three passes. More knots and more windows
    # is the cheap part of the fix; the honest part is that this constant is
    # measured less precisely than the arithmetic ones and the predictor should
    # not lean on it to break ties.
    dispatch_k = 4 * k
    ws = Float64[]
    ts = Float64[]
    knots = unique(Int[2, 3, 4, max(4, budget ÷ 3), max(4, budget ÷ 2), budget])
    for w in knots
        w > budget && continue
        t = timed_min(k = dispatch_k) do
            PP.threaded_foreach_worker_persistent(:cost_calib_dispatch, w, w) do _w, i
                @inbounds sink[i] = i * 1.0000001
            end
            @inbounds sink[1]
        end
        push!(ws, Float64(w))
        push!(ts, t)
    end
    slope, intercept = _least_squares_line(ws, ts)
    verbose && println("[calibrate] dispatch     base=$(round(intercept, digits=1)) ns  " *
                       "per_worker=$(round(slope, digits=1)) ns")

    # Atomic cost: identical work under the dynamic and static schedulers at
    # chunk 1, so the difference is one contended read-modify-write per item.
    t_static = timed_min(k = dispatch_k) do
        PP.threaded_foreach_worker_persistent(:cost_calib_sched, n_items, budget; scheduler = :static) do _w, i
            @inbounds sink[i] = i * 1.0000001
        end
        @inbounds sink[1]
    end
    t_dynamic = timed_min(k = dispatch_k) do
        PP.threaded_foreach_worker_persistent(:cost_calib_sched, n_items, budget; scheduler = :dynamic, chunk = 1) do _w, i
            @inbounds sink[i] = i * 1.0000001
        end
        @inbounds sink[1]
    end
    atomic_ns = max(0.0, t_dynamic - t_static) / n_items
    verbose && println("[calibrate] atomic       -> $(round(atomic_ns, digits=5)) ns " *
                       "(static $(round(t_static, digits=0)) ns, dynamic $(round(t_dynamic, digits=0)) ns)")

    return (per_node, intercept, slope, atomic_ns)
end

function _least_squares_line(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    n >= 2 || return (0.0, isempty(y) ? 0.0 : y[1])
    x̄ = sum(x) / n
    ȳ = sum(y) / n
    sxx = sum((xi - x̄)^2 for xi in x)
    sxx <= 0.0 && return (0.0, ȳ)
    slope = sum((x[i] - x̄) * (y[i] - ȳ) for i in 1:n) / sxx
    return (slope, ȳ - slope * x̄)
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
        "ns_per_scalar_item" => mc.ns_per_scalar_item,
        "ns_per_queue_node" => mc.ns_per_queue_node,
        "dispatch_pool_ns_base" => mc.dispatch_pool_ns_base,
        "dispatch_pool_ns_per_worker" => mc.dispatch_pool_ns_per_worker,
        "dispatch_batch_ns_base" => mc.dispatch_batch_ns_base,
        "dispatch_batch_ns_per_worker" => mc.dispatch_batch_ns_per_worker,
        "ns_per_atomic" => mc.ns_per_atomic,
        "reference_fma_ns" => mc.reference_fma_ns,
        "reference_mem_ns" => mc.reference_mem_ns,
        "simd_lane" => Dict{String, Any}(
            "log2_size" => mc.simd_lane.log2_size,
            "ns" => mc.simd_lane.ns,
        ),
        "coeff_touch" => Dict{String, Any}(
            "log2_size" => mc.coeff_touch.log2_size,
            "ns" => mc.coeff_touch.ns,
        ),
        "usl_alpha_base" => mc.usl_alpha_base,
        "usl_beta_arith" => mc.usl_beta_arith,
        "usl_beta_alloc" => mc.usl_beta_alloc,
        "usl_beta_bw" => mc.usl_beta_bw,
        "llc_bytes" => mc.llc_bytes,
        "parallel_speedup" => Dict{String, Any}(
            "log2_size" => mc.parallel_speedup.log2_size,
            "ns" => mc.parallel_speedup.ns,
        ),
    )
    path_s = String(path)
    mkpath(dirname(path_s))
    tmp = path_s * ".tmp"
    open(tmp, "w") do io
        TOML.print(io, payload)
    end
    mv(tmp, path_s; force = true)
    return path_s
end

"""
    load_machine_constants(path = machine_constants_path()) -> Union{Nothing, MachineConstants}

Read cached constants, or `nothing` if absent, unparseable, or written by a
different schema version. A stale-schema file is discarded rather than
reinterpreted: term semantics changed by definition when the version moved, so
its numbers do not mean what the current predictor would take them to mean.
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
        curve(key) = RateCurve(
            Vector{Float64}(parsed[key]["log2_size"]),
            Vector{Float64}(parsed[key]["ns"]),
        )
        return MachineConstants(
            simd_lane = curve("simd_lane"),
            coeff_touch = curve("coeff_touch"),
            parallel_speedup = curve("parallel_speedup"),
            usl_alpha_base = Float64(get(parsed, "usl_alpha_base", 0.0)),
            usl_beta_arith = Float64(get(parsed, "usl_beta_arith", 0.0)),
            usl_beta_alloc = Float64(get(parsed, "usl_beta_alloc", 0.0)),
            usl_beta_bw = Float64(get(parsed, "usl_beta_bw", 0.0)),
            llc_bytes = Float64(get(parsed, "llc_bytes", 0.0)),
            ns_per_scalar_item = Float64(parsed["ns_per_scalar_item"]),
            ns_per_queue_node = Float64(parsed["ns_per_queue_node"]),
            dispatch_pool_ns_base = Float64(parsed["dispatch_pool_ns_base"]),
            dispatch_pool_ns_per_worker = Float64(parsed["dispatch_pool_ns_per_worker"]),
            dispatch_batch_ns_base = Float64(parsed["dispatch_batch_ns_base"]),
            dispatch_batch_ns_per_worker = Float64(parsed["dispatch_batch_ns_per_worker"]),
            ns_per_atomic = Float64(parsed["ns_per_atomic"]),
            reference_fma_ns = Float64(parsed["reference_fma_ns"]),
            reference_mem_ns = Float64(parsed["reference_mem_ns"]),
            fingerprint = String(get(parsed, "fingerprint", "")),
            schema_version = CALIBRATION_SCHEMA_VERSION,
        )
    catch
        return nothing
    end
end

"""
    constants_are_current(mc; tolerance = 0.25, k = 9) -> (ok, measured_ns, drift)

The staleness canary. Re-measures the cheap arithmetic reference and compares it
against its value at calibration time.

This is the *only* freshness check the model needs, and that is a consequence of
the caching rule rather than an optimism: what is cached is expensive and stable
(these constants), and what is volatile is re-derived every run (per-effector
probes). So the single thing that can silently go stale is whether these numbers
still describe this machine -- and they stop doing so not when the hardware
changes but when its *effective* capacity does: a different cgroup quota, a
changed governor, a co-tenant, thermal state.

Costs microseconds. A drift beyond `tolerance` means abstain and keep the
heuristic route, not recalibrate mid-run: recalibrating on a machine that is
currently contended would bake the contention into the constants.

The default tolerance is deliberately loose. Measured on this box, ten busy
loops on twelve cores moved the reference about 5%, and a co-tenant that
degrades everything uniformly does not invalidate a model whose consumers use
ratios. What 25% catches is the categorical error -- constants from a different
machine, or a quota that makes the recorded throughput unreachable.
"""
function constants_are_current(
    mc::MachineConstants;
    tolerance::Float64 = 0.25,
    k::Int = 9,
)::NamedTuple{(:ok, :measured_ns, :drift), Tuple{Bool, Float64, Float64}}
    measured = reference_kernel_ns(k = k)
    if !(isfinite(mc.reference_fma_ns) && mc.reference_fma_ns > 0.0)
        return (ok = false, measured_ns = measured, drift = NaN)
    end
    drift = measured / mc.reference_fma_ns - 1.0
    return (ok = abs(drift) <= tolerance, measured_ns = measured, drift = drift)
end

# Polyester dispatch, measured separately because `satellite_batch` uses
# `@batch` rather than the channel pool. Polyester keeps its own parked worker
# set and hands out contiguous index ranges without going through a Channel, so
# its fan-out is substantially cheaper -- which is exactly why applying the
# pool's constant to `satellite_batch` mispriced it.
#
# `minbatch` is set the way dynamics_rhs.jl sets it, so the measured cost is the
# cost of the dispatch shape the RHS actually issues rather than of Polyester in
# the abstract.
function _calibrate_polyester_dispatch(k::Int, verbose::Bool)
    budget = Threads.nthreads()
    n = 4096
    sink = zeros(Float64, n)

    if budget <= 1 || Polyester.num_cores() <= 1
        verbose && println("[calibrate] polyester    skipped: single-threaded session")
        return (NaN, NaN)
    end

    ws = Float64[]
    ts = Float64[]
    for w in unique(Int[2, 3, 4, max(4, budget ÷ 3), max(4, budget ÷ 2), budget])
        w > budget && continue
        # minbatch chosen so the loop splits into exactly w chunks.
        mb = max(1, cld(n, w))
        t = timed_min(k = 4 * k) do
            @batch minbatch=mb for i in 1:n
                @inbounds sink[i] = i * 1.0000001
            end
            @inbounds sink[1]
        end
        push!(ws, Float64(w))
        push!(ts, t)
    end
    slope, intercept = _least_squares_line(ws, ts)
    verbose && println("[calibrate] polyester    base=$(round(intercept, digits=1)) ns  " *
                       "per_worker=$(round(slope, digits=1)) ns")
    return (intercept, slope)
end
