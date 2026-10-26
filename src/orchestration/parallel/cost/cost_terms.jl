# Machine-constant and contention types for the parallel planners.
#
# `MachineConstants` holds the two USL terms `calibrate_machine` measures once
# per machine, which the predictive campaign planner reads.
# `ContentionInputs` holds what one warm RHS pass shows about how a workload
# scales; the native-lock width cap reads its lock occupancy.

"""
    ContentionInputs

What one warm RHS pass reveals about how this workload will behave under a
split, as opposed to how much work it contains.

It answers "how badly does this work scale", which needs different evidence
from how much work there is. It measures two of the mechanisms found to cap
thread scaling in this codebase:

- `alloc_bytes_per_pass` / `pass_ns`  allocation pressure. Shared-heap
                    allocation and the GC it drives are process-wide, which is
                    why separate processes escape a ceiling that more threads
                    only worsen. The measured contrast is stark and mechanical:
                    a vacuum constellation iterates a homogeneous effector
                    1-tuple and allocates 0 bytes per pass, while the same
                    physics plus atmosphere reaches tens of GiB over a solve and
                    inverts past four threads.
- `lock_*`          native-lock occupancy, read straight off the shared
                    GRAM/SPICE critical section. Unlike allocation this is
                    not a proxy: `lock_hold_ns` over the window is an Amdahl
                    serial fraction with a hard `1/rho` speedup ceiling, and
                    `lock_wait_ns / lock_hold_ns` reports whether the width in
                    use is already past what the lock admits.

`valid` is false when the warm pass could not be completed, in which case every
field is meaningless rather than zero -- a workload that allocates nothing and
one that could not be measured must not read the same.
"""
Base.@kwdef struct ContentionInputs
    alloc_bytes_per_pass::Float64 = 0.0
    pass_ns::Float64 = 0.0
    n_active_sats::Int = 0
    lock_acquisitions::Int64 = 0
    lock_hold_ns::Float64 = 0.0
    lock_wait_ns::Float64 = 0.0
    valid::Bool = false
end

"""
    alloc_bytes_per_second(ci) -> Float64

Allocation rate of one worker running this workload. Zero for an invalid or
zero-length measurement, never `Inf`: a caller sizing a width from this must get
"no pressure observed" rather than a poisoned number.
"""
@inline function alloc_bytes_per_second(ci::ContentionInputs)::Float64
    (ci.valid && ci.pass_ns > 0.0) || return 0.0
    return ci.alloc_bytes_per_pass * 1.0e9 / ci.pass_ns
end

"""
    lock_duty_cycle(ci, workers = 1) -> Float64

Fraction rho of available worker-time spent inside the native lock. `1/rho`
bounds achievable speedup at any width.
"""
@inline function lock_duty_cycle(ci::ContentionInputs, workers::Integer = 1)::Float64
    window = ci.pass_ns * max(1.0, Float64(workers))
    (ci.valid && window > 0.0) || return 0.0
    return clamp(ci.lock_hold_ns / window, 0.0, 1.0)
end

"""
    lock_wait_hold_ratio(ci) -> Float64

Time blocked outside the native lock against time spent inside it. Much greater
than one means the split is already wider than the lock admits, and the finding
needs no history to reach: it is visible in the window that measured it.
"""
@inline function lock_wait_hold_ratio(ci::ContentionInputs)::Float64
    (ci.valid && ci.lock_hold_ns > 0.0) || return 0.0
    return ci.lock_wait_ns / ci.lock_hold_ns
end

# 4: adds the USL contention parameters (usl_alpha_base, usl_beta_arith,
# usl_beta_alloc, usl_beta_bw, llc_bytes). Bumped rather than defaulted because
# a file without them would silently predict zero contention -- which reads as
# "this machine scales perfectly" and is the most dangerous possible default.
#
# 5: the USL arithmetic ladder and the parallel speedup curve no longer read
# past their table (`calib_batch_kernel!` wraps its column). Constants saved at
# 4 timed the memory after the table, so they are measured again, not kept.
#
# 6: only the constants something reads are kept: `usl_alpha_base`,
# `usl_beta_alloc`, the fingerprint and the measuring thread count. The other
# rates were measured for the RHS-plan predictor, which 0.6.0 removed. A file at
# 5 is measured again rather than read for its two kept values, so the loader
# reads one version only.
const CALIBRATION_SCHEMA_VERSION = 6

"""
    MachineConstants

This machine's parallel cost constants, measured by [`calibrate_machine`](@ref)
and cached in a fingerprinted TOML file (`machine_constants_path`). The
predictive campaign planner reads the two USL terms.

- `usl_alpha_base`    serial fraction of the dispatch primitive itself, from the
                      arithmetic speedup ladder: the USL `α` every workload pays.
- `usl_beta_alloc`    contention per GiB/s of allocation per worker: the excess
                      USL `β` an allocation ladder shows over the arithmetic one.
- `fingerprint`       machine identity these constants were measured on.
- `schema_version`    bumped when the file's meaning changes, so stale files are
                      measured again rather than reinterpreted. Defaults to the
                      current schema, so a value built here saves a file the
                      loader reads.
- `measured_threads`  `Threads.nthreads()` of the session that measured them, the
                      widest rung of its ladders; 0 when unknown, which is never
                      current.
"""
Base.@kwdef struct MachineConstants
    usl_alpha_base::Float64 = 0.0
    usl_beta_alloc::Float64 = 0.0
    fingerprint::String = ""
    schema_version::Int = CALIBRATION_SCHEMA_VERSION
    measured_threads::Int = 0
end

"""
    usl_speedup(α, β, k) -> Float64

Universal Scalability Law: `S(k) = k / (1 + α(k−1) + βk(k−1))`.

Two parameters rather than one because the measured curves do two different
things. `α` alone (Amdahl) can only saturate, and the workloads that matter here
do not saturate -- they **invert**: the aero RHS at N=1024 runs 4957 µs serial,
1637 at four threads and 4755 at twelve. Only a term that grows faster than the
worker count can bend a speedup curve back down, which is what `βk(k−1)` does,
and it puts the turning point at `k* = √((1−α)/β)` in closed form rather than
requiring that the inverting width have been measured.

Returns at least `1.0`: a plan can always decline to split, so a predicted
speedup below one means "do not widen", not "widening makes the work itself
slower".
"""
@inline function usl_speedup(alpha::Real, beta::Real, k::Integer)::Float64
    kf = Float64(max(1, k))
    kf <= 1.0 && return 1.0
    a = max(0.0, Float64(alpha))
    b = max(0.0, Float64(beta))
    denom = 1.0 + a * (kf - 1.0) + b * kf * (kf - 1.0)
    denom > 0.0 || return 1.0
    return max(1.0, kf / denom)
end

"""
    usl_peak_workers(α, β) -> Float64

The width at which [`usl_speedup`](@ref) turns over, `k* = √((1−α)/β)`. `Inf`
when `β` is zero, which is the saturating case: more workers never help, but
they never hurt either.
"""
@inline function usl_peak_workers(alpha::Real, beta::Real)::Float64
    b = max(0.0, Float64(beta))
    b > 0.0 || return Inf
    a = clamp(Float64(alpha), 0.0, 1.0)
    return sqrt(max(0.0, 1.0 - a) / b)
end
