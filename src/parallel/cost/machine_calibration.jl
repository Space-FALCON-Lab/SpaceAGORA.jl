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

const CALIBRATION_SCHEMA_VERSION = 1

# Footprint ladder for the coefficient-touch curve: square-ish tables, matching
# the (L+2) x (M+2) shape of a full-field harmonics model, spanning L1-resident
# to well past L3.
const _COEFF_KNOTS = ((22, 20), (52, 50), (102, 100), (202, 200),
                      (402, 400), (802, 800), (1602, 1600))

# Batch widths for the SIMD-lane curve. The low end is where per-pass loop
# overhead dominates and cannot amortize -- a real regime, since batch width is
# satellites-per-worker and a small constellation on many threads lands there.
const _LANE_KNOTS = (16, 64, 256, 1024, 4096, 16384, 65536, 262144)

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
# unreadable. cgroup v2 first, then v1.
function _cgroup_cpu_quota()::Float64
    try
        v2 = "/sys/fs/cgroup/cpu.max"
        if isfile(v2)
            parts = split(strip(read(v2, String)))
            length(parts) == 2 || return -1.0
            parts[1] == "max" && return -1.0
            return parse(Float64, parts[1]) / parse(Float64, parts[2])
        end
        quota_p = "/sys/fs/cgroup/cpu/cpu.cfs_quota_us"
        period_p = "/sys/fs/cgroup/cpu/cpu.cfs_period_us"
        if isfile(quota_p) && isfile(period_p)
            q = parse(Float64, strip(read(quota_p, String)))
            p = parse(Float64, strip(read(period_p, String)))
            q <= 0 && return -1.0
            return q / p
        end
    catch
        return -1.0
    end
    return -1.0
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
    node_ns, dispatch_base, dispatch_per_worker, atomic_ns = _calibrate_dispatch(k, verbose)

    return MachineConstants(
        simd_lane = lane_curve,
        coeff_touch = coeff_curve,
        ns_per_scalar_item = scalar_ns,
        ns_per_queue_node = node_ns,
        dispatch_ns_base = dispatch_base,
        dispatch_ns_per_worker = dispatch_per_worker,
        ns_per_atomic = atomic_ns,
        reference_fma_ns = reference_kernel_ns(k = k),
        reference_mem_ns = reference_memory_kernel_ns(k = k),
        fingerprint = machine_fingerprint(),
        schema_version = CALIBRATION_SCHEMA_VERSION,
    )
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

    log2_size = Float64[]
    ns = Float64[]
    for B in _LANE_KNOTS
        out = zeros(Float64, B)
        total = timed_min(() -> calib_batch_kernel!(out, tab, rows ÷ 2, n_touch, spt); k = k)
        lane_terms = Float64(n_touch * spt * B)
        per_lane = max(0.0, total - touch_ns) / lane_terms
        push!(log2_size, log2(Float64(B)))
        push!(ns, per_lane)
        verbose && println("[calibrate] simd_lane    B=$(B) -> $(round(per_lane, digits=5)) ns/lane-term")
    end
    return RateCurve(log2_size, ns)
end

function _calibrate_scalar_item(k::Int, verbose::Bool)::Float64
    state = [1.0 + i * 1e-6 for i in 1:4096]
    n = 4096
    per_item = timed_min(() -> calib_scalar_kernel!(state, n); k = k) / n
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
            PP.threaded_foreach_worker(w, w) do _w, i
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
        PP.threaded_foreach_worker(n_items, budget; scheduler = :static) do _w, i
            @inbounds sink[i] = i * 1.0000001
        end
        @inbounds sink[1]
    end
    t_dynamic = timed_min(k = dispatch_k) do
        PP.threaded_foreach_worker(n_items, budget; scheduler = :dynamic, chunk = 1) do _w, i
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
        "dispatch_ns_base" => mc.dispatch_ns_base,
        "dispatch_ns_per_worker" => mc.dispatch_ns_per_worker,
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
            ns_per_scalar_item = Float64(parsed["ns_per_scalar_item"]),
            ns_per_queue_node = Float64(parsed["ns_per_queue_node"]),
            dispatch_ns_base = Float64(parsed["dispatch_ns_base"]),
            dispatch_ns_per_worker = Float64(parsed["dispatch_ns_per_worker"]),
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
