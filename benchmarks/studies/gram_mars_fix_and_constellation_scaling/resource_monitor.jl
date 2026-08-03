# Process-scoped CPU/RAM sampling shared by leo_constellation_size_scaling.jl
# (self-monitoring: samples the current process, since all (n_sats, mode)
# points run in one process) and leo_thread_scaling_by_mode.jl (samples a
# spawned worker subprocess by pid).
#
# Sampling runs as a genuinely separate OS process (a small `bash` polling
# loop around `ps -p <pid>`) rather than a Julia Task, so it keeps sampling
# on its own schedule even while every Julia thread -- including thread 1,
# where Julia's own task scheduler would otherwise have to find an opening --
# is pinned busy inside a solve (e.g. Polyester-batched harmonics). `ps -p
# <pid>` reports only that one process's RSS/%CPU, so other load sharing the
# machine (background browser tabs, another terminal, etc.) never leaks into
# the numbers -- this deliberately does NOT read system-wide CPU/memory.

struct ResourceUsage
    peak_rss_mb::Float64
    mean_cpu_pct::Float64
    peak_cpu_pct::Float64
    n_samples::Int
end

const _EMPTY_RESOURCE_USAGE = ResourceUsage(NaN, NaN, NaN, 0)

struct ResourceMonitor
    proc::Base.Process
    log_path::String
end

# Starts the background sampler for `pid`, polling every `interval_s` seconds
# until the caller calls `stop_and_collect!`. `pid` defaults to the current
# process (`getpid()`), the right choice for a script that runs many
# (n_sats, mode) points back-to-back in one process; pass a worker
# subprocess's pid instead when monitoring it externally.
function start_resource_monitor(pid::Integer=getpid(); interval_s::Real=0.2)::ResourceMonitor
    log_path = tempname()
    touch(log_path)
    script = """
    while kill -0 $(pid) 2>/dev/null; do
        ps -o rss=,pcpu= -p $(pid) >> "$(log_path)" 2>/dev/null
        sleep $(interval_s)
    done
    """
    proc = run(pipeline(`bash -c $(script)`; stdout=devnull, stderr=devnull); wait=false)
    return ResourceMonitor(proc, log_path)
end

# Stops sampling and parses the accumulated log into a ResourceUsage. Safe to
# call even if no samples were ever collected (e.g. a point that finished
# faster than one sampling interval) -- returns NaN fields with n_samples=0
# rather than throwing, since "too fast to sample" is a real, expected case
# for small n_sats.
function stop_and_collect!(mon::ResourceMonitor)::ResourceUsage
    try
        Base.process_running(mon.proc) && kill(mon.proc)
        wait(mon.proc)
    catch
    end

    rss_samples = Float64[]
    cpu_samples = Float64[]
    if isfile(mon.log_path)
        for line in eachline(mon.log_path)
            parts = split(strip(line))
            length(parts) == 2 || continue
            rss_kb = tryparse(Float64, parts[1])
            cpu_pct = tryparse(Float64, parts[2])
            (rss_kb === nothing || cpu_pct === nothing) && continue
            push!(rss_samples, rss_kb / 1024)
            push!(cpu_samples, cpu_pct)
        end
        rm(mon.log_path; force=true)
    end

    isempty(rss_samples) && return _EMPTY_RESOURCE_USAGE
    return ResourceUsage(
        maximum(rss_samples),
        sum(cpu_samples) / length(cpu_samples),
        maximum(cpu_samples),
        length(rss_samples),
    )
end

# Convenience wrapper: runs `f()` while sampling `pid`'s resource usage,
# returning `(f(), ResourceUsage)`. Used where the monitored process is the
# caller itself (pid defaults to getpid()).
function measure_resource_usage(f, pid::Integer=getpid(); interval_s::Real=0.2)
    mon = start_resource_monitor(pid; interval_s=interval_s)
    local result
    local usage
    try
        result = f()
    finally
        usage = stop_and_collect!(mon)
    end
    return result, usage
end
