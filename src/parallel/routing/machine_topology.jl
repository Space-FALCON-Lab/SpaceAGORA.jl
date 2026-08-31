# How many cores this process may actually use.
#
# Everything in `src/` currently sizes parallelism from `Sys.CPU_THREADS`, which
# answers a different question than the one routing asks. It reports the logical
# CPUs the kernel exposes on the host, not the cores this process is allowed to
# run on and not the cores that add throughput:
#
#   - SMT siblings are counted. On the 12-physical/24-logical reference box every
#     workload measured regressed from 16 to 24 threads, so a ladder topped out
#     at CPU_THREADS puts its worst point last. `process_max_workers` defaults to
#     CPU_THREADS and is therefore twice the useful width on that machine.
#   - An affinity mask is ignored. Under `taskset -c 0-3` the process sees 24 and
#     may use 4.
#   - A cgroup quota is ignored. In a container the host's core count is visible
#     while the process is throttled to a fraction of it.
#
# The physical-core ladder exists today only in the benchmark harness
# (`SPACEAGORA_PPC_PHYSICAL_CORES`), which is the wrong side of the boundary: the
# harness measures the runtime, so a fact the runtime needs cannot live there.
#
# NOTHING HERE CHANGES AN EXISTING DECISION. `usable_core_budget` is computed,
# cached and reported; no current profile consults it. R4/R5 keep sizing from
# `Sys.CPU_THREADS` so that every measurement taken against them stays
# reproducible, and the allocator that will consult this is the only thing that
# should move because of it.

# Readings are cached because they cannot change within a process for the two
# that matter (topology, quota) and because the third (affinity) is only ever
# changed deliberately. `refresh_machine_topology!` exists for the deliberate
# case and for tests.
const _TOPOLOGY_CACHE = Ref{Union{Nothing, NamedTuple}}(nothing)
const _TOPOLOGY_LOCK = ReentrantLock()

@inline function _topology_env(name::String)::String
    return strip(get(ENV, name, ""))
end

function _positive_int_env(name::String)::Int
    raw = _topology_env(name)
    isempty(raw) && return -1
    value = tryparse(Int, raw)
    if value === nothing || value <= 0
        throw(ArgumentError("$name must be a positive integer, got '$raw'"))
    end
    return value
end

# ── Physical cores ────────────────────────────────────────────────────────────

# (package, core) pairs keyed by logical CPU id, from sysfs. Empty when sysfs is
# unavailable (non-Linux, or a kernel built without topology export).
function _sysfs_core_map()::Dict{Int, Tuple{Int, Int}}
    out = Dict{Int, Tuple{Int, Int}}()
    base = "/sys/devices/system/cpu"
    isdir(base) || return out
    for entry in readdir(base)
        m = match(r"^cpu(\d+)$", entry)
        m === nothing && continue
        cpu = parse(Int, m.captures[1])
        core_path = joinpath(base, entry, "topology", "core_id")
        pkg_path = joinpath(base, entry, "topology", "physical_package_id")
        (isfile(core_path) && isfile(pkg_path)) || continue
        try
            core = parse(Int, strip(read(core_path, String)))
            pkg = parse(Int, strip(read(pkg_path, String)))
            out[cpu] = (pkg, core)
        catch
            # A single unreadable CPU must not void the whole map; it only makes
            # the count conservative.
            continue
        end
    end
    return out
end

# Fallback for kernels that expose /proc/cpuinfo but not the sysfs topology
# tree. Same (physical id, core id) pairing, parsed from the flat record format.
function _proc_cpuinfo_core_map()::Dict{Int, Tuple{Int, Int}}
    out = Dict{Int, Tuple{Int, Int}}()
    isfile("/proc/cpuinfo") || return out
    processor = -1
    pkg = 0
    core = -1
    try
        for line in eachline("/proc/cpuinfo")
            if isempty(strip(line))
                if processor >= 0 && core >= 0
                    out[processor] = (pkg, core)
                end
                processor = -1; pkg = 0; core = -1
                continue
            end
            parts = split(line, ':'; limit=2)
            length(parts) == 2 || continue
            key = strip(parts[1]); val = strip(parts[2])
            if key == "processor"
                processor = something(tryparse(Int, val), -1)
            elseif key == "physical id"
                pkg = something(tryparse(Int, val), 0)
            elseif key == "core id"
                core = something(tryparse(Int, val), -1)
            end
        end
        if processor >= 0 && core >= 0
            out[processor] = (pkg, core)
        end
    catch
        return Dict{Int, Tuple{Int, Int}}()
    end
    return out
end

function _darwin_physical_cores()::Int
    try
        raw = read(`sysctl -n hw.physicalcpu`, String)
        value = tryparse(Int, strip(raw))
        return value === nothing ? -1 : max(1, value)
    catch
        return -1
    end
end

"""
    physical_core_count() -> Int

Number of distinct physical cores visible to this process, ignoring SMT
siblings. Falls back to `Sys.CPU_THREADS` where topology cannot be read, which
is conservative in the wrong direction but never zero.

`SPACEAGORA_PHYSICAL_CORES` overrides. The benchmark harness has needed that
override under `taskset` since the B-series runs; this is the same knob, read by
the runtime rather than duplicated in the harness.
"""
function physical_core_count()::Int
    override = _positive_int_env("SPACEAGORA_PHYSICAL_CORES")
    override > 0 && return override

    if Sys.islinux()
        map = _sysfs_core_map()
        isempty(map) && (map = _proc_cpuinfo_core_map())
        if !isempty(map)
            return max(1, length(Set(values(map))))
        end
    elseif Sys.isapple()
        n = _darwin_physical_cores()
        n > 0 && return n
    end
    return max(1, Sys.CPU_THREADS)
end

# ── Affinity ──────────────────────────────────────────────────────────────────

# Parse a kernel CPU list ("0-3,8,12-15") into the set of logical CPU ids.
function _parse_cpu_list(raw::AbstractString)::Set{Int}
    out = Set{Int}()
    for token in split(strip(raw), ',')
        t = strip(token)
        isempty(t) && continue
        if occursin('-', t)
            bounds = split(t, '-'; limit=2)
            lo = tryparse(Int, strip(bounds[1]))
            hi = tryparse(Int, strip(bounds[2]))
            (lo === nothing || hi === nothing) && continue
            for c in lo:hi
                push!(out, c)
            end
        else
            v = tryparse(Int, t)
            v === nothing || push!(out, v)
        end
    end
    return out
end

"""
    allowed_cpus() -> Set{Int}

Logical CPU ids this process may be scheduled on, from
`/proc/self/status`'s `Cpus_allowed_list`. Empty when unreadable, which callers
must treat as "unknown", never as "none".

Read from `/proc` rather than through `sched_getaffinity` so this needs no
`ccall` and no libc struct layout assumption.
"""
function allowed_cpus()::Set{Int}
    isfile("/proc/self/status") || return Set{Int}()
    try
        for line in eachline("/proc/self/status")
            if startswith(line, "Cpus_allowed_list:")
                return _parse_cpu_list(split(line, ':'; limit=2)[2])
            end
        end
    catch
        return Set{Int}()
    end
    return Set{Int}()
end

# ── cgroup quota ──────────────────────────────────────────────────────────────

"""
    cgroup_cpu_quota() -> Float64

Effective CPU quota as a fractional core count, or `-1.0` when unlimited or
unreadable. cgroup v2 (`cpu.max`) first, then v1
(`cpu.cfs_quota_us` / `cpu.cfs_period_us`).

This is the canonical implementation; `ParallelCost.machine_calibration` reads
it through here so a container's quota cannot be described two ways in one
process. Its value semantics are unchanged from that earlier copy, so machine
fingerprints computed against it stay stable and cached constants remain valid.
"""
function cgroup_cpu_quota()::Float64
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

# ── Combined budget ───────────────────────────────────────────────────────────

# Physical cores reachable under the affinity mask.
#
# Counting the mask's entries directly would answer in logical CPUs and mix
# units with the physical count: pinned to 4 SMT threads of 2 cores, the mask
# says 4 and the machine delivers 2. Mapping the mask through the topology and
# counting distinct (package, core) pairs answers in the same unit as
# `physical_core_count`, which is the unit the budget is expressed in.
function _affinity_physical_cores(total_physical::Int)::Int
    allowed = allowed_cpus()
    isempty(allowed) && return -1
    map = _sysfs_core_map()
    isempty(map) && (map = _proc_cpuinfo_core_map())
    if isempty(map)
        # No topology to map through: the mask's own size is the only bound
        # available, and it is an upper bound in physical terms.
        return max(1, length(allowed))
    end
    cores = Set{Tuple{Int, Int}}()
    for cpu in allowed
        entry = get(map, cpu, nothing)
        entry === nothing && continue
        push!(cores, entry)
    end
    isempty(cores) && return -1
    return min(max(1, length(cores)), max(1, total_physical))
end

"""
    machine_topology(; refresh = false) -> NamedTuple

Everything known about this machine's usable parallelism, measured once and
cached:

  - `cpu_threads`       `Sys.CPU_THREADS`, the number every current caller uses.
  - `physical_cores`    distinct physical cores, SMT siblings collapsed.
  - `affinity_cores`    physical cores reachable under the affinity mask, or
                        `-1` when no mask could be read.
  - `cgroup_quota`      fractional core quota, or `-1.0` when unlimited.
  - `usable_cores`      the binding minimum of the above — the budget.
  - `smt_ratio`         `cpu_threads / physical_cores`; `1.0` means no SMT.
  - `source`            `:override` when `SPACEAGORA_CORE_BUDGET` set the
                        budget, otherwise which term bound it (`:physical`,
                        `:affinity`, `:quota`).
"""
function machine_topology(; refresh::Bool = false)::NamedTuple
    if !refresh
        cached = _TOPOLOGY_CACHE[]
        cached === nothing || return cached
    end
    return lock(_TOPOLOGY_LOCK) do
        if !refresh
            cached = _TOPOLOGY_CACHE[]
            cached === nothing || return cached
        end
        threads = max(1, Sys.CPU_THREADS)
        physical = physical_core_count()
        affinity = _affinity_physical_cores(physical)
        quota = cgroup_cpu_quota()

        usable = physical
        source = :physical
        if affinity > 0 && affinity < usable
            usable = affinity
            source = :affinity
        end
        if quota > 0.0
            quota_cores = max(1, Int(floor(quota)))
            if quota_cores < usable
                usable = quota_cores
                source = :quota
            end
        end

        override = _positive_int_env("SPACEAGORA_CORE_BUDGET")
        if override > 0
            usable = override
            source = :override
        end

        snapshot = (
            cpu_threads = threads,
            physical_cores = physical,
            affinity_cores = affinity,
            cgroup_quota = quota,
            usable_cores = max(1, usable),
            smt_ratio = physical > 0 ? threads / physical : 1.0,
            source = source,
        )
        _TOPOLOGY_CACHE[] = snapshot
        return snapshot
    end
end

"""
    usable_core_budget() -> Int

The number of cores an allocation may spend, as opposed to the number the host
advertises. This is the quantity `w · k ≤ P` is written against.
"""
@inline usable_core_budget()::Int = machine_topology().usable_cores

"""
    refresh_machine_topology!() -> NamedTuple

Re-read topology, affinity and quota. Affinity is the only one that can move
within a process (a deliberate re-pin), and tests need to observe env overrides
taking effect.
"""
refresh_machine_topology!()::NamedTuple = machine_topology(refresh = true)
