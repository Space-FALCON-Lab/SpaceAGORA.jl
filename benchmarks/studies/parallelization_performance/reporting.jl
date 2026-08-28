ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
using Plots
using Plots.PlotMeasures: mm

function ppc_hardware_snapshot()
    return (
        timestamp_utc=string(now(UTC)),
        machine=gethostname(),
        julia_version=string(VERSION),
        cpu_threads=Sys.CPU_THREADS,
        julia_threads=Threads.nthreads(),
        os=string(Sys.KERNEL),
        arch=string(Sys.ARCH),
        git_commit=ppc_git_commit(),
        # Machine load at the moment this row was produced, and how much of it is
        # not ours. A benchmark that cannot get the cores it was given measures
        # contention, not routing -- and nothing in the recorded row previously
        # said so, which is how a 26-hour job from another user on a shared
        # 64-core box silently turned an entire run into noise (the adaptive
        # profiles, which dispatch the most concurrent work, degraded ~4x while
        # the narrow static ones degraded ~2x, which reads exactly like a routing
        # failure). Recording it per row makes that visible after the fact;
        # ppc_assert_machine_quiet! refuses to start when it is already true.
        load_avg_1min=ppc_load_average(),
        load_headroom=ppc_load_headroom()
    )
end

"""
    ppc_load_average() -> Float64

1-minute load average, or `NaN` where the platform does not report one.
"""
function ppc_load_average()::Float64
    try
        return Float64(Sys.loadavg()[1])
    catch
        return NaN
    end
end

# Fraction of the machine's physical cores not already claimed by the 1-minute
# load. 1.0 is idle; <= 0.0 means the machine is oversubscribed before this
# benchmark starts.
function ppc_load_headroom()::Float64
    la = ppc_load_average()
    isnan(la) && return NaN
    cores = _ppc_physical_core_count()
    return (cores - la) / max(1, cores)
end

"""
    ppc_cpu_busy_fraction(sample_s=2.0) -> Float64

Fraction of total CPU capacity in use, sampled directly from `/proc/stat` over
`sample_s` seconds. `NaN` where /proc/stat is unavailable.

Distinct from the load average, and needed because the load average could not see
the failure this exists to catch. On 2026-08-24 a leftover probe process and a
desktop file indexer together held ~1.8 of 12 cores for ten hours; the 1-minute
load average read 4.09, i.e. 66% headroom, comfortably inside the guard's 50%
threshold, so every run started and every timing came back 8-34% slow -- serial
baselines included, which is the signature of stolen cores rather than a code
change. Load average answers "how many runnable tasks", which a benchmark's own
workers dominate; this answers "how much of the machine is already gone", which
is the question that matters before any worker exists.
"""
function ppc_cpu_busy_fraction(sample_s::Float64=2.0; cpus::Vector{Int}=Int[])::Float64
    # With `cpus` empty this reads the aggregate "cpu" line (whole machine). With
    # a CPU list it sums only those per-core lines, which is the capacity a run
    # confined by taskset/cgroup can actually obtain -- see ppc_assert_machine_quiet!.
    wanted = isempty(cpus) ? nothing : Set("cpu" * string(c) for c in cpus)
    read_stat() = try
        if wanted === nothing
            line = first(eachline("/proc/stat"))
            vals = parse.(Int64, split(line)[2:end])
            idle = vals[4] + (length(vals) >= 5 ? vals[5] : 0)   # idle + iowait
            (sum(vals), idle)
        else
            total = 0; idle = 0; seen = 0
            for line in eachline("/proc/stat")
                parts = split(line)
                (isempty(parts) || !(parts[1] in wanted)) && continue
                vals = parse.(Int64, parts[2:end])
                total += sum(vals)
                idle += vals[4] + (length(vals) >= 5 ? vals[5] : 0)
                seen += 1
            end
            seen == 0 ? nothing : (total, idle)
        end
    catch
        nothing
    end
    a = read_stat()
    a === nothing && return NaN
    sleep(max(0.1, sample_s))
    b = read_stat()
    b === nothing && return NaN
    dtotal = b[1] - a[1]
    didle = b[2] - a[2]
    dtotal <= 0 && return NaN
    return clamp(1.0 - didle / dtotal, 0.0, 1.0)
end

"""
    ppc_pinned_cpus() -> Vector{Int}

CPUs this process is confined to, or empty when it can use the whole machine.

Read from `/proc/self/status` `Cpus_allowed_list`, so it reflects whatever set
the run actually got -- `taskset`, a cgroup cpuset, or a scheduler placement --
rather than anything the harness was told. Empty is the "not confined" answer,
which keeps every caller on the pre-existing whole-machine path.
"""
function ppc_pinned_cpus()::Vector{Int}
    Sys.islinux() || return Int[]
    try
        for line in eachline("/proc/self/status")
            startswith(line, "Cpus_allowed_list:") || continue
            cpus = _ppc_parse_cpu_list(strip(split(line, ":"; limit=2)[2]))
            # A mask covering everything is not a confinement; treat it as absent
            # so the whole-machine checks below apply unchanged.
            length(cpus) >= Sys.CPU_THREADS && return Int[]
            return cpus
        end
    catch
    end
    return Int[]
end

"""
    ppc_assert_machine_quiet!(; required_headroom, allow_override)

Refuse to start a timing run on a machine that is already busy.

Timing benchmarks are only meaningful when the process can actually obtain the
cores it asks for. Set `SPACEAGORA_PPC_ALLOW_BUSY=1` to proceed anyway (recorded
in the row either way through `load_avg_1min`).
"""
function ppc_assert_machine_quiet!(; required_headroom::Float64=0.5, wait_s::Int=600)
    allow_busy = lowercase(strip(get(ENV, "SPACEAGORA_PPC_ALLOW_BUSY", "0"))) in ("1", "true", "yes", "on")
    # Announce before any blocking wait. This check runs ahead of the controller's
    # first println, so without this a run that is legitimately waiting for the
    # machine to settle is indistinguishable from a hung process: no output, no
    # worker, no CPU. That ambiguity cost two killed runs on 2026-08-28.
    println("[parallelization-performance] machine-quiet check (waits up to $(wait_s)s)")
    flush(stdout)

    pinned = ppc_pinned_cpus()
    if !isempty(pinned)
        # Confined to a subset of the machine, so judge only that subset.
        #
        # The 1-minute load average is machine-wide and has no per-cpuset
        # equivalent, so dividing it by the pinned core count compares a whole
        # machine's runnable tasks against a slice of its capacity. On a 64-core
        # box shared with another user, a run pinned to 12 idle cores saw load
        # 13.4 -> "headroom -11%" and refused to start, while the cores it had
        # were completely free. Sampling the pinned cores answers the question
        # the guard is actually asking: is the capacity this run can obtain
        # already spoken for?
        deadline = time() + max(0, wait_s)
        busy = ppc_cpu_busy_fraction(; cpus=pinned)
        while !isnan(busy) && busy > 0.15 && time() < deadline
            sleep(15)
            busy = ppc_cpu_busy_fraction(; cpus=pinned)
        end
        pinned_desc = "$(length(pinned)) pinned core(s) [$(first(pinned))-$(last(pinned))]"
        if isnan(busy) || busy <= 0.15
            println("[parallelization-performance] $(pinned_desc) " *
                    "$(isnan(busy) ? "unreadable" : string(round(100 * busy; digits=0)) * "% busy") — proceeding")
            flush(stdout)
            return nothing
        end
        offenders = try
            readchomp(pipeline(`ps -eo pcpu,etime,comm --sort=-pcpu`, `head -4`))
        catch
            "(process list unavailable)"
        end
        msg = "Pinned cores are $(round(100 * busy; digits=0))% busy before this run starts any " *
              "worker ($(pinned_desc)). Timing results taken now would measure contention. " *
              "Top consumers:\n$(offenders)"
        allow_busy || error(msg)
        @warn msg
        return nothing
    end

    la = ppc_load_average()
    isnan(la) && return nothing
    # Wait for the machine to settle before refusing.
    #
    # The 1-minute load average is a decaying average, so it stays elevated for
    # minutes after work stops. A controller starting immediately after a long
    # phase therefore sees the *previous* phase's load and aborts, which is what
    # happened to B14 on 2026-08-24: it failed with 0 s elapsed at 7.6 load
    # seconds after B11 finished a 5 h 20 m phase, on an otherwise idle machine.
    # Refusing to measure on a busy machine is right; refusing to wait 90 s for
    # a decaying average is not.
    deadline = time() + max(0, wait_s)
    while ppc_load_headroom() < required_headroom && time() < deadline
        sleep(15)
    end
    # Foreign-CPU check, run before any worker of this run exists: at this point
    # anything busy is somebody else's, so a simple total-utilisation sample is
    # an honest reading. Threshold at 15% of the machine, above desktop idle
    # noise and below the ~15% a single stolen core represents on 12 cores.
    busy = ppc_cpu_busy_fraction()
    if !isnan(busy) && busy > 0.15
        offenders = try
            readchomp(pipeline(`ps -eo pcpu,etime,comm --sort=-pcpu`, `head -4`))
        catch
            "(process list unavailable)"
        end
        msg = "Machine has $(round(100 * busy; digits=0))% of its CPU already in use before this run " *
              "starts any worker. Timing results taken now would measure contention. Top consumers:\n$(offenders)"
        if allow_busy
            @warn msg
        else
            error(msg)
        end
    end
    headroom = ppc_load_headroom()
    headroom >= required_headroom && return nothing
    cores = _ppc_physical_core_count()
    msg = "Machine is busy: 1-min load average $(round(la; digits=1)) on $(cores) " *
          "physical cores (headroom $(round(100 * headroom; digits=0))%, need " *
          "$(round(100 * required_headroom; digits=0))%). Timing results taken now " *
          "would measure contention rather than parallel routing. Wait for the " *
          "machine to clear (waited $(wait_s) s), or set SPACEAGORA_PPC_ALLOW_BUSY=1 to proceed anyway."
    if allow_busy
        @warn msg
        return nothing
    end
    error(msg)
end

function ppc_git_commit()::String
    try
        return strip(read(`git -C $(PPC_REPO_ROOT) rev-parse HEAD`, String))
    catch
        return "unknown"
    end
end

function ppc_write_rows(path::String, rows::Vector{NamedTuple})
    mkpath(dirname(path))
    df = DataFrame(rows)
    CSV.write(path, df)
    return path
end

function ppc_read_optional(path::String)::DataFrame
    isfile(path) || return DataFrame()
    return CSV.read(path, DataFrame)
end

function ppc_summarize(raw::DataFrame, parity::DataFrame)::DataFrame
    nrow(raw) == 0 && return DataFrame()
    grouped = groupby(raw, [:case, :family, :mode, :thread_count, :mc_samples])
    summary = combine(
        grouped,
        :success => (x -> count(identity, skipmissing(x))) => :success_count,
        :success => length => :sample_count,
        :wall_time_s => mean => :wall_time_mean_s,
        :wall_time_s => std => :wall_time_std_s,
        :wall_time_s => (x -> length(x) <= 1 ? 0.0 : 100.0 * std(x) / max(abs(mean(x)), eps(Float64))) => :wall_time_cv_pct,
        :throughput_samples_per_s => mean => :throughput_samples_per_s_mean,
        :execution_scope => (x -> first(skipmissing(x))) => :execution_scope,
        :outer_backend_actual => (x -> first(skipmissing(x))) => :outer_backend_actual
    )
    summary[!, :wall_time_ci95_low_s] = summary.wall_time_mean_s .- 1.96 .* coalesce.(summary.wall_time_std_s, 0.0) ./ sqrt.(summary.sample_count)
    summary[!, :wall_time_ci95_high_s] = summary.wall_time_mean_s .+ 1.96 .* coalesce.(summary.wall_time_std_s, 0.0) ./ sqrt.(summary.sample_count)

    serial = summary[summary.mode .== "serial", [:case, :mc_samples, :wall_time_mean_s]]
    rename!(serial, :wall_time_mean_s => :serial_wall_time_mean_s)
    summary = leftjoin(summary, serial; on=[:case, :mc_samples])
    summary[!, :speedup_vs_serial] = [
        (ismissing(s) || t <= 0.0) ? missing : s / t
        for (s, t) in zip(summary.serial_wall_time_mean_s, summary.wall_time_mean_s)
    ]
    summary[!, :parallel_efficiency] = [
        (v isa Missing || th <= 0) ? missing : v / max(1, th)
        for (v, th) in zip(summary.speedup_vs_serial, summary.thread_count)
    ]

    if nrow(parity) > 0 && "pass" in names(parity)
        parity_ok = combine(groupby(parity, [:case, :mode]), :pass => (x -> all(Bool.(x))) => :parity_pass)
        summary = leftjoin(summary, parity_ok; on=[:case, :mode])
    else
        summary[!, :parity_pass] = fill(missing, nrow(summary))
    end
    summary[!, :headline_eligible] = [
        ok && (p === missing || p === true)
        for (ok, p) in zip(summary.success_count .== summary.sample_count, summary.parity_pass)
    ]
    return summary
end

function ppc_write_report(path::String, cfg::PPCConfig, raw::DataFrame, summary::DataFrame, parity::DataFrame)
    open(path, "w") do io
        println(io, "# SpaceAGORA Parallelization Performance Study")
        println(io)
        println(io, "- Generated UTC: $(now(UTC))")
        println(io, "- Profile: `$(cfg.profile)`")
        println(io, "- Solver mode: `$(cfg.solver_mode)`")
        println(io, "- Repeats: `$(cfg.repeats)`, warmup: `$(cfg.warmup)`")
        println(io, "- Modes: `$(join(cfg.modes, ", "))`")
        println(io, "- Thread counts: `$(join(cfg.threads, ", "))`")
        println(io)
        println(io, "## Outputs")
        println(io)
        println(io, "- Raw rows: `$(nrow(raw))`")
        println(io, "- Summary rows: `$(nrow(summary))`")
        println(io, "- Trajectory parity rows: `$(nrow(parity))`")
        println(io)
        println(io, "## Headline Summary")
        println(io)
        if nrow(summary) == 0
            println(io, "No summary rows were generated.")
        else
            cols = intersect([:case, :mode, :thread_count, :mc_samples, :execution_scope, :outer_backend_actual, :wall_time_mean_s, :throughput_samples_per_s_mean, :speedup_vs_serial, :headline_eligible], Symbol.(names(summary)))
            show(io, MIME"text/plain"(), first(select(summary, cols), min(30, nrow(summary))))
            println(io)
        end
        println(io)
        println(io, "## Trajectory Correctness")
        println(io)
        if nrow(parity) == 0
            println(io, "No trajectory parity rows were generated.")
        else
            cols = intersect([:case, :mode, :pass, :pos_rel_rms, :pos_rel_max, :vel_rel_rms, :vel_rel_max, :event_count_equal, :event_time_abs_max_s], Symbol.(names(parity)))
            show(io, MIME"text/plain"(), select(parity, cols))
            println(io)
        end
    end
    return path
end

function ppc_write_plots(outdir::String, summary::DataFrame, parity::DataFrame)::Vector{String}
    paths = String[]
    try
        ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
        @eval using Plots
        if nrow(summary) > 0 && "speedup_vs_serial" in names(summary)
            df = summary[summary.mode .!= "serial", :]
            if nrow(df) > 0
                p = Plots.scatter(
                    string.(df.case),
                    coalesce.(df.speedup_vs_serial, NaN);
                    group=string.(df.mode),
                    xrotation=35,
                    xlabel="case",
                    ylabel="speedup vs R0",
                    legend=:outertopright,
                    tickfont=Plots.font(8),
                    guidefont=Plots.font(12),
                    legend_font=Plots.font(9),
                    size=(1500, 850),
                    left_margin=18mm,
                    right_margin=34mm,
                    top_margin=10mm,
                    bottom_margin=36mm
                )
                path = joinpath(outdir, "parallelization_speedup_by_case.png")
                Plots.savefig(p, path)
                push!(paths, path)
            end
        end
        if nrow(summary) > 0 && all(name in names(summary) for name in ("speedup_vs_serial", "execution_scope"))
            inner = summary[(summary.mode .!= "serial") .& (summary.execution_scope .== "single_simulation"), :]
            if nrow(inner) > 0
                p = Plots.scatter(
                    string.(inner.case),
                    coalesce.(inner.speedup_vs_serial, NaN);
                    group=string.(inner.mode),
                    xrotation=35,
                    xlabel="case",
                    ylabel="single-simulation speedup vs R0",
                    legend=:outertopright,
                    tickfont=Plots.font(8),
                    guidefont=Plots.font(12),
                    legend_font=Plots.font(9),
                    size=(1500, 850),
                    left_margin=18mm,
                    right_margin=34mm,
                    top_margin=10mm,
                    bottom_margin=36mm
                )
                path = joinpath(outdir, "parallelization_inner_speedup_by_case.png")
                Plots.savefig(p, path)
                push!(paths, path)
            end
        end
        if nrow(summary) > 0 && all(name in names(summary) for name in ("throughput_samples_per_s_mean", "outer_backend_actual"))
            outer = summary[(summary.mode .!= "serial") .& (summary.mc_samples .> 1) .& (summary.outer_backend_actual .!= "serial"), :]
            if nrow(outer) > 0
                p = Plots.scatter(
                    outer.mc_samples,
                    outer.throughput_samples_per_s_mean;
                    group=string.(outer.mode) .* " / " .* string.(outer.case),
                    xscale=:log10,
                    xlabel="Monte Carlo samples",
                    ylabel="throughput (samples/s)",
                    legend=:outertopright,
                    tickfont=Plots.font(8),
                    guidefont=Plots.font(12),
                    legend_font=Plots.font(8),
                    size=(1500, 850),
                    left_margin=18mm,
                    right_margin=48mm,
                    top_margin=10mm,
                    bottom_margin=18mm
                )
                path = joinpath(outdir, "parallelization_outer_throughput_by_samples.png")
                Plots.savefig(p, path)
                push!(paths, path)
            end
        end
        if nrow(parity) > 0 && "pos_rel_max" in names(parity)
            p = Plots.scatter(
                string.(parity.case),
                parity.pos_rel_max;
                group=string.(parity.mode),
                xrotation=45,
                yscale=:log10,
                ylabel="max relative position error",
                legend=:outertopright,
                size=(1200, 700)
            )
            path = joinpath(outdir, "parallelization_parity_position_error.png")
            Plots.savefig(p, path)
            push!(paths, path)
        end
    catch err
        @warn "Plot generation failed" exception=(err, catch_backtrace())
    end
    return paths
end
