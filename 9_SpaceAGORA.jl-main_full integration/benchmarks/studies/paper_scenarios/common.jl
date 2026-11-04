# Controller-side helpers shared by the paper_scenarios suite (s1_*.jl .. s5_*.jl).
#
# Deliberately does NOT load SpaceAGORA: controllers only build environment-variable
# grids, spawn `julia --threads=N scenario_worker.jl` subprocesses (thread count is
# fixed at Julia startup, so every (point, thread-count) pair must be its own
# process), parse the worker's single PS_RESULT line, and write CSV rows. All
# simulation work happens in scenario_worker.jl.
#
# Environment knobs honored by every controller:
#   PS_THREADS    parallel thread budget per point (default Sys.CPU_THREADS; set this
#                 to the machine's PHYSICAL core count for paper-quality numbers)
#   PS_REPEATS    timed repeats per point (default 3)
#   PS_WARMUP     warmup runs per point (default 1)
#   PS_TIMEOUT_S  per-point wall-clock kill deadline (default 3600)
#   PS_MISSION_S  mission length override (per-scenario defaults otherwise)

using Printf
using Statistics

const PS_DIR = @__DIR__
const PS_REPO_ROOT = normpath(joinpath(PS_DIR, "..", "..", ".."))
const PS_WORKER_SCRIPT = joinpath(PS_DIR, "scenario_worker.jl")

ps_env_int(name::String, default::Int)::Int = parse(Int, get(ENV, name, string(default)))
ps_env_float(name::String, default::Float64)::Float64 = parse(Float64, get(ENV, name, string(default)))
ps_env_str(name::String, default::String)::String = get(ENV, name, default)

ps_threads_default()::Int = ps_env_int("PS_THREADS", Sys.CPU_THREADS)
ps_repeats()::Int = ps_env_int("PS_REPEATS", 3)
ps_warmup()::Int = ps_env_int("PS_WARMUP", 1)
ps_timeout_s()::Float64 = ps_env_float("PS_TIMEOUT_S", 3600.0)

"""Results directory for this host, optionally suffixed via PS_RESULTS_SUFFIX
(e.g. "_thread_ladder") so an exploratory run's CSVs land in their own
pseudo-host directory -- auto-discovered by plot_results.jl as a separate
host -- instead of overwriting the canonical per-host baseline data that
lives under the bare hostname."""
function ps_results_dir()::String
    dir = joinpath(PS_DIR, "results", gethostname() * ps_env_str("PS_RESULTS_SUFFIX", ""))
    mkpath(dir)
    return dir
end

"""Env pairs for one constellation-workload point (mirrors the proven
leo_2048_constellation_gram_scaling.jl / mc_multisat_thread_allocation_sample.jl
combinations, including the GRAM safety knobs those studies established)."""
function ps_constellation_env(;
    rhs_mode::String,             # "serial" | "flat"
    density::String,              # "none" | "gram_standard" | "gram_lookahead" | "gram_surrogate"
    mission_s::Float64,
    inner_budget::Int,
)::Vector{Pair{String, String}}
    pairs = Pair{String, String}[
        "SPACEAGORA_RHS_EXECUTION_MODE" => rhs_mode,
        "SPACEAGORA_HARMONICS_BATCH_ENABLED" => (rhs_mode == "serial" ? "0" : "1"),
        "SPACEAGORA_EFFECTOR_PARALLEL" => (rhs_mode == "serial" ? "off" : "auto"),
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0",
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "0",
        "SPACEAGORA_RHS_CALIBRATE" => "off",
        "SPACEAGORA_INNER_THREAD_BUDGET" => string(inner_budget),
    ]
    if density == "gram_standard"
        # Native GRAM without the look-ahead cache: freeze density per accepted step
        # (per-RK-stage perturbation noise otherwise collapses the adaptive dt) and
        # keep density sampling serial — concurrent direct native calls across
        # satellites serialize on the process-wide lock at best and are the exact
        # hazard the surrogate/look-ahead paths exist to avoid.
        append!(pairs, [
            "SPACEAGORA_DENSITY_FREEZE_PER_STEP" => "1",
            "SPACEAGORA_VACUUM_GRAM_CACHE" => "0",
            "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "off",
        ])
    elseif density == "gram_lookahead"
        # Vacuum-predicted trajectory cache. Horizon + deviation loosened past what
        # the mission can reach so only the (proven-safe) initial build ever runs —
        # workaround for the known native rebuild hang with 2+ satellites (see
        # [[project_gram_multisat_rebuild_bug]] / leo_2048 header).
        append!(pairs, [
            "SPACEAGORA_DENSITY_FREEZE_PER_STEP" => "0",
            "SPACEAGORA_VACUUM_GRAM_CACHE" => "1",
            "SPACEAGORA_VACUUM_GRAM_CACHE_NPOINTS" => "20",
            "SPACEAGORA_VACUUM_GRAM_CACHE_HORIZON_S" => string(mission_s + 500.0),
            "SPACEAGORA_VACUUM_GRAM_CACHE_DEVIATION_M" => "1e8",
            "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => (rhs_mode == "serial" ? "off" : "auto"),
        ])
    elseif density == "gram_surrogate"
        # Lock-free offline surrogate. Force density-callback threading "on": the
        # generic auto-gate's thread floor would otherwise leave sampling serial at
        # small thread counts (see mc_multisat_thread_allocation_sample.jl).
        push!(pairs, "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => (rhs_mode == "serial" ? "off" : "on"))
    end
    return pairs
end

"""Env pairs for one Monte Carlo point (outer campaign × inner RHS split)."""
function ps_mc_env(;
    density::String,
    mission_s::Float64,
    outer_backend::String,        # "serial" | "threads" | "process"
    outer_workers::Int,
    inner_budget::Int,
)::Vector{Pair{String, String}}
    pairs = ps_constellation_env(
        rhs_mode=(inner_budget > 1 ? "flat" : "serial"),
        density=density,
        mission_s=mission_s,
        inner_budget=inner_budget,
    )
    # Under an outer *thread* split, tell inner policies they share the pool. The
    # process backend must NOT see this flag pre-set — its own dispatch checks it
    # to detect being nested inside an enclosing split and would collapse to serial.
    outer_active = (outer_backend == "threads" && outer_workers > 1) ? "1" : "0"
    return [p.first == "SPACEAGORA_OUTER_PARALLEL_ACTIVE" ? (p.first => outer_active) : p for p in pairs]
end

"""Spawn one scenario_worker.jl subprocess, stream output to a log file, parse the
PS_RESULT line. Returns a Dict{String,String} (empty on failure/timeout)."""
function ps_run_worker(;
    label::String,
    julia_threads::Int,
    env::Vector{Pair{String, String}},
    log_dir::String,
    timeout_s::Float64=ps_timeout_s(),
)::Dict{String, String}
    julia_bin = Base.julia_cmd().exec[1]
    cmd = Cmd([
        julia_bin,
        "--project=$(PS_REPO_ROOT)",
        "--threads=$(julia_threads)",
        PS_WORKER_SCRIPT,
    ])
    cmd = addenv(cmd, env...)
    log_path = joinpath(log_dir, "$(label).log")
    result = Dict{String, String}()
    open(log_path, "w") do log_io
        proc = run(pipeline(cmd; stdout=log_io, stderr=log_io); wait=false)
        deadline = time() + timeout_s
        while process_running(proc) && time() < deadline
            sleep(1.0)
        end
        if process_running(proc)
            # Reap the subprocess's own children first (orphaned Distributed pool
            # workers hold the stdout pipe open past the parent's death, which was
            # observed to hang wait(proc) indefinitely in earlier studies).
            try
                run(`pkill -9 -P $(getpid(proc))`; wait=false)
            catch
            end
            kill(proc, Base.SIGKILL)
            println("  -> TIMEOUT after $(timeout_s)s (log: $(log_path))")
            return result
        end
    end
    for line in eachline(log_path)
        if startswith(line, "PS_RESULT ")
            for token in split(chopprefix(line, "PS_RESULT "))
                kv = split(token, "="; limit=2)
                length(kv) == 2 && (result[kv[1]] = kv[2])
            end
        end
    end
    isempty(result) && println("  -> FAILED (no PS_RESULT line; log: $(log_path))")
    return result
end

"""Append rows (Vector of Pairs col=>value, uniform keys) to a CSV file."""
function ps_write_csv(path::String, rows::Vector{<:AbstractDict{String, String}}, columns::Vector{String})
    isempty(rows) && return nothing
    open(path, "w") do io
        println(io, join(columns, ","))
        for row in rows
            println(io, join([get(row, c, "") for c in columns], ","))
        end
    end
    println("Wrote $(length(rows)) row(s) to $(path)")
    return nothing
end

"""Merge PS_RESULT dict with the point's own metadata columns."""
function ps_row(meta::Vector{Pair{String, String}}, result::Dict{String, String})::Dict{String, String}
    row = Dict{String, String}(meta...)
    for (k, v) in result
        row[k] = v
    end
    haskey(row, "ok") || (row["ok"] = "false")
    return row
end

ps_hostname_meta()::Vector{Pair{String, String}} = [
    "hostname" => gethostname(),
    "cpu_threads" => string(Sys.CPU_THREADS),
]
