module RuntimeServices

# Despite living in separate files (`libGRAM.dylib`'s statically-linked CSPICE
# copy vs. SPICE.jl's own libcspice), the two are not isolated at the OS
# symbol-resolution level: `nm -gU libGRAM.dylib` shows it globally exports
# CSPICE internals (`chkin_`, `chkout_`, `trcpkg_`, `subslr_c`, ...) under the
# same symbol names SPICE.jl's libcspice also exports. Concurrently invoking
# a native-GRAM call (previously serialized only by GRAM_LOCK) on one thread
# and a SpaceAGORA ephemerides/frame-transform call (serialized only by
# SPICE_LOCK) on another corrupts CSPICE's internal call-trace stack
# (SPICE(NAMESDONOTMATCH)/CHKOUT errors), because both locks were guarding
# what is actually one shared C-level critical section. GRAM_LOCK and
# SPICE_LOCK must therefore be the same lock, not merely same-shaped.
const SPICE_LOCK = ReentrantLock()
const GRAM_LOCK = SPICE_LOCK

# ── Native-lock occupancy ─────────────────────────────────────────────────────
#
# How much of a solve is spent inside the shared native critical section, and
# how much is spent waiting to get into it.
#
# This is the one contention term in the parallel cost model that does not have
# to be inferred from proxies. Allocation pressure and memory footprint are
# predicted from countable workload quantities; lock occupancy can simply be
# read, and the two quantities it yields answer different questions:
#
#   hold / (elapsed * workers)   the serial fraction. A lock caps achievable
#                                speedup at 1/rho no matter how wide the split,
#                                which is an Amdahl term rather than the
#                                superlinear one the cost model currently
#                                declines to represent.
#   wait / hold                  over-width, directly. Much more waiting than
#                                holding means more workers than the lock
#                                admits, so the inner width should shrink. No
#                                sweep, no arms and no history are needed to
#                                read it.
#
# ATTRIBUTED BY SITE, because GRAM_LOCK === SPICE_LOCK means one lock serves two
# subsystems. A live-GRAM constellation with SPICE third bodies contends against
# itself across both, and a workload's density-model family does not predict its
# occupancy: a vacuum + SPICE third-body run has real occupancy and no GRAM at
# all, while a surrogate run names GRAM and takes the lock almost never.
#
# All bookkeeping state below is written only while the lock is held, so it
# needs no atomics and no lock of its own. `_stats_snapshot` takes the lock to
# read.

const NATIVE_LOCK_SITES = (
    :gram_density,      # native GRAM point/profile density sampling
    :gram_cache,        # GRAM track-cache refresh
    :gram_setup,        # GRAM model construction, warm-up, deepcopy
    :spice_ephemeris,   # body states, third-body positions
    :spice_frame,       # frame transforms, reference-system conversion
    :spice_body,        # body constants: radii, GM, frame names
    :other,
)

Base.@kwdef mutable struct NativeLockSiteStats
    acquisitions::Int64 = 0
    contended::Int64 = 0
    hold_ns::Int64 = 0
    wait_ns::Int64 = 0
end

const _NATIVE_LOCK_STATS =
    NativeLockSiteStats[NativeLockSiteStats() for _ in 1:length(NATIVE_LOCK_SITES)]

# Reentrancy bookkeeping. A ReentrantLock lets the owning task re-acquire, and
# charging hold time once per nesting level would over-count occupancy by
# whatever the call depth happens to be. Only the outermost acquisition is
# timed; the site recorded is the outermost one, since that is the span the
# other tasks are actually excluded for.
const _native_lock_owner = Ref{Union{Nothing, Task}}(nothing)
const _native_lock_depth = Ref{Int}(0)
const _native_lock_entry_ns = Ref{Int64}(0)
const _native_lock_entry_wait_ns = Ref{Int64}(0)
const _native_lock_entry_site = Ref{Int}(length(NATIVE_LOCK_SITES))

# Two `time_ns()` calls and a `trylock` per outermost acquisition, against a
# critical section that exists to guard microsecond-scale FFI. Off restores the
# bare `lock(f, SPICE_LOCK)` path exactly, and exists so the overhead can be
# measured rather than asserted.
const _NATIVE_LOCK_STATS_ENABLED = Ref{Bool}(true)

# When the counters were last zeroed, so occupancy can be expressed against a
# window rather than against an unknown span. A duty cycle is a fraction of
# something; without the denominator the counters are only comparable to
# themselves.
const _NATIVE_LOCK_RESET_NS = Ref{UInt64}(time_ns())

function __init__()
    raw = lowercase(strip(get(ENV, "SPACEAGORA_NATIVE_LOCK_STATS", "1")))
    _NATIVE_LOCK_STATS_ENABLED[] = !(raw in ("0", "false", "no", "off"))
    return nothing
end

@inline function native_lock_site_index(site::Symbol)::Int
    site === :gram_density    && return 1
    site === :gram_cache      && return 2
    site === :gram_setup      && return 3
    site === :spice_ephemeris && return 4
    site === :spice_frame     && return 5
    site === :spice_body      && return 6
    return 7
end

"""
    TrackedLock(site)

A view of [`SPICE_LOCK`](@ref) that records occupancy against `site`.

It is a view and not a lock: every instance locks the one shared
`ReentrantLock`, so two `TrackedLock`s with different sites still exclude each
other, exactly as the underlying native critical section requires. The site only
decides which counter the time lands in.

Constructing one is free (it holds a single `Int`), which matters because the
GRAM hot path builds one per density call. It subtypes `Base.AbstractLock` so
`lock(f, l)` and `@lock` work unchanged, and so it can be handed to GRAMSuite as
a `lock_obj` — that keyword is untyped, whereas `GRAMSuite._GRAM_DEFAULT_LOCK_HOOK`
is a `Ref{Union{Nothing, ReentrantLock}}` and must keep receiving the plain lock.
"""
struct TrackedLock <: Base.AbstractLock
    site::Int
end

TrackedLock(site::Symbol) = TrackedLock(native_lock_site_index(site))

"""
    tracked_lock(site) -> TrackedLock

The shared native lock, attributed to `site`. See [`NATIVE_LOCK_SITES`](@ref)
for the accepted names; anything else is charged to `:other` rather than
rejected, so a new call site cannot break a run by naming itself.
"""
@inline tracked_lock(site::Symbol)::TrackedLock = TrackedLock(site)

@inline function _enter_tracked!(site::Int, wait_ns::Int64)::Nothing
    if _native_lock_depth[] > 0 && _native_lock_owner[] === current_task()
        _native_lock_depth[] += 1
    else
        _native_lock_owner[] = current_task()
        _native_lock_depth[] = 1
        _native_lock_entry_site[] = site
        _native_lock_entry_wait_ns[] = wait_ns
        _native_lock_entry_ns[] = time_ns()
    end
    return nothing
end

function Base.lock(l::TrackedLock)
    if !_NATIVE_LOCK_STATS_ENABLED[]
        lock(SPICE_LOCK)
        return nothing
    end
    wait_ns = Int64(0)
    if !trylock(SPICE_LOCK)
        t0 = time_ns()
        lock(SPICE_LOCK)
        wait_ns = Int64(time_ns() - t0)
    end
    _enter_tracked!(l.site, wait_ns)
    return nothing
end

function Base.trylock(l::TrackedLock)::Bool
    if !_NATIVE_LOCK_STATS_ENABLED[]
        return trylock(SPICE_LOCK)
    end
    trylock(SPICE_LOCK) || return false
    _enter_tracked!(l.site, Int64(0))
    return true
end

function Base.unlock(::TrackedLock)
    # Branch on the bookkeeping state rather than on the enable flag: a flag
    # flipped between lock and unlock (tests do this) must not leave a
    # half-recorded acquisition or skip the release.
    if _native_lock_depth[] > 0 && _native_lock_owner[] === current_task()
        if _native_lock_depth[] > 1
            _native_lock_depth[] -= 1
        else
            st = @inbounds _NATIVE_LOCK_STATS[_native_lock_entry_site[]]
            st.acquisitions += 1
            st.hold_ns += Int64(time_ns() - _native_lock_entry_ns[])
            wait = _native_lock_entry_wait_ns[]
            if wait > 0
                st.wait_ns += wait
                st.contended += 1
            end
            _native_lock_depth[] = 0
            _native_lock_owner[] = nothing
        end
    end
    unlock(SPICE_LOCK)
    return nothing
end

Base.islocked(::TrackedLock) = islocked(SPICE_LOCK)

"""
    with_native_lock(f, site)

Run `f` holding the shared native lock, charging the time to `site`. The
call-site form of [`tracked_lock`](@ref); prefer it over `lock(SPICE_LOCK)` in
`src/` so that occupancy is attributable.
"""
@inline function with_native_lock(f::F, site::Symbol) where {F}
    return lock(f, tracked_lock(site))
end

"""
    native_lock_stats_enabled() -> Bool
    set_native_lock_stats!(on) -> Bool

Whether occupancy is being recorded. `SPACEAGORA_NATIVE_LOCK_STATS=0` disables
it at startup; the setter exists for tests and for the overhead A/B.
"""
@inline native_lock_stats_enabled()::Bool = _NATIVE_LOCK_STATS_ENABLED[]

function set_native_lock_stats!(on::Bool)::Bool
    previous = _NATIVE_LOCK_STATS_ENABLED[]
    _NATIVE_LOCK_STATS_ENABLED[] = on
    return previous
end

"""
    reset_native_lock_stats!()

Zero every site's counters. Called at the start of a measured region; occupancy
is only meaningful relative to a known elapsed window.
"""
function reset_native_lock_stats!()::Nothing
    lock(SPICE_LOCK) do
        _NATIVE_LOCK_RESET_NS[] = time_ns()
        for st in _NATIVE_LOCK_STATS
            st.acquisitions = 0
            st.contended = 0
            st.hold_ns = 0
            st.wait_ns = 0
        end
    end
    return nothing
end

"""
    native_lock_stats_snapshot() -> NamedTuple

Per-site counters plus totals. `hold_ns` is time inside the critical section,
`wait_ns` is time blocked outside it, and `contended` counts the acquisitions
that had to wait at all.

Divide `hold_ns` by the wall time of the region and the worker count to get the
duty cycle rho; `wait_ns / hold_ns` reports whether the current width is already
past what the lock admits.
"""
function native_lock_stats_snapshot()
    return lock(SPICE_LOCK) do
        sites = NamedTuple(
            site => (
                acquisitions = st.acquisitions,
                contended = st.contended,
                hold_ns = st.hold_ns,
                wait_ns = st.wait_ns,
            )
            for (site, st) in zip(NATIVE_LOCK_SITES, _NATIVE_LOCK_STATS)
        )
        acquisitions = sum(st.acquisitions for st in _NATIVE_LOCK_STATS; init = Int64(0))
        hold_ns = sum(st.hold_ns for st in _NATIVE_LOCK_STATS; init = Int64(0))
        wait_ns = sum(st.wait_ns for st in _NATIVE_LOCK_STATS; init = Int64(0))
        contended = sum(st.contended for st in _NATIVE_LOCK_STATS; init = Int64(0))
        return (
            enabled = _NATIVE_LOCK_STATS_ENABLED[],
            acquisitions = acquisitions,
            contended = contended,
            hold_ns = hold_ns,
            wait_ns = wait_ns,
            wait_hold_ratio = hold_ns > 0 ? wait_ns / hold_ns : 0.0,
            sites = sites,
        )
    end
end

"""
    native_lock_occupancy(; workers = 1) -> NamedTuple

Occupancy over the window since the counters were last reset: `rho`, the
`wait_hold` ratio, the acquisition count and the window itself.

This is the per-solve form. `run_simulation` resets at the start of a solve, so
a caller asking mid-solve or just after gets that solve's numbers rather than
whatever the process has accumulated since it started -- which is the difference
between a workload property and a process-lifetime statistic.

`rho` is an Amdahl serial fraction: speedup is bounded by `1/rho` at any width.
`wait_hold` says whether the width in use is already past what the lock admits.
"""
function native_lock_occupancy(; workers::Integer = 1)
    snap = native_lock_stats_snapshot()
    window_ns = Float64(time_ns() - _NATIVE_LOCK_RESET_NS[])
    denom = window_ns * max(1.0, Float64(workers))
    return (
        rho = denom > 0.0 ? clamp(snap.hold_ns / denom, 0.0, 1.0) : 0.0,
        wait_hold = snap.wait_hold_ratio,
        acquisitions = snap.acquisitions,
        window_s = window_ns / 1.0e9,
    )
end

"""
    lock_width_ceiling(; workers = 1, floor_rho = 0.02) -> Int

The widest inner split the native lock admits: `ceil(1/rho)`, or `typemax(Int)`
when occupancy is below `floor_rho`.

A lock caps achievable speedup at `1/rho` however wide the split, so any width
above that ceiling adds contention and no throughput. Measured on live GRAM
density, `rho = 0.89` gives a ceiling of 2 -- threads cannot help it, which is
why process isolation is the only route there.

`floor_rho` exists so the ceiling is not computed from noise: below it the
implied ceiling exceeds any real core budget anyway, so returning "no
constraint" is both cheaper and more honest than returning a large number
derived from a handful of acquisitions.
"""
function lock_width_ceiling(; workers::Integer = 1, floor_rho::Float64 = 0.02)::Int
    occ = native_lock_occupancy(workers = workers)
    occ.rho >= floor_rho || return typemax(Int)
    return max(1, ceil(Int, 1.0 / occ.rho))
end

"""
    native_lock_duty_cycle(elapsed_ns, workers) -> Float64

Occupancy rho: the fraction of the available worker-time that was spent inside
the native critical section. `1/rho` is the ceiling on achievable speedup for
any width, which is what makes this an Amdahl term rather than a contention one.

Returns `0.0` for a non-positive window, never `Inf` or `NaN`: a caller sizing a
width from this must get "no constraint observed", not a poisoned number.
"""
function native_lock_duty_cycle(elapsed_ns::Real, workers::Integer = 1)::Float64
    elapsed = Float64(elapsed_ns) * max(1.0, Float64(workers))
    elapsed > 0.0 || return 0.0
    snap = native_lock_stats_snapshot()
    return clamp(snap.hold_ns / elapsed, 0.0, 1.0)
end

end # module RuntimeServices
