# Controller for the GRAM memory-footprint study.
#
# Why this study exists: the outer-route memory model in
# `src/parallel/routing/machine_topology.jl` charges a per-spacecraft term to
# every process worker on a native-GRAM constellation, and that term decides
# whether the process route is offered at all. It has to be a measured number
# per density path, not one number for "GRAM".
#
# One subprocess per (density path, constellation size) point -- a shared
# process would carry the previous point's heap into the next point's baseline,
# which is precisely the quantity being measured. Every point runs
# `--threads=1`, the configuration a process-pool worker is started in.
#
# Usage:
#   julia --project=. benchmarks/studies/gram_memory_footprint/run_footprint.jl
#
# Environment knobs:
#   GMF_SIZES       comma-separated constellation sizes (default 16,64,256,1024,4096)
#   GMF_PATHS       comma-separated paths (default none,gram_point,gram_freeze_step,
#                                                  gram_lookahead,gram_surrogate)
#   GMF_MISSION_S   mission length per point [s] (default 300)
#   GMF_GRAVITY     invsq | l20 | l50 (default l20)
#   GMF_TIMEOUT_S   per-point kill deadline (default 2400)
#   GMF_REPEATS     repeats of the whole grid, appended as separate rows (default 1)
#   GMF_MIN_FREE_GB skip a point when the kernel reports less available memory
#                   than this (default 12)
#   GMF_SOLVE_REPEATS  solves of the same size per point before the peak is
#                   read (default 1; the paper_scenarios workers run four)
#   GMF_ISOLATE     1 = run_simulation's production `isolate_state` deep copy
#                   (default), 0 = the paper_scenarios convention
#   GMF_SUFFIX      results file-name suffix, for an exploratory grid

using Printf

const GMF_DIR = @__DIR__
const GMF_REPO_ROOT = normpath(joinpath(GMF_DIR, "..", "..", ".."))
const GMF_WORKER = joinpath(GMF_DIR, "footprint_worker.jl")

gmf_env_str(name::String, default::String)::String = get(ENV, name, default)
gmf_env_int(name::String, default::Int)::Int = parse(Int, get(ENV, name, string(default)))
gmf_env_float(name::String, default::Float64)::Float64 = parse(Float64, get(ENV, name, string(default)))

const SIZES = parse.(Int, split(gmf_env_str("GMF_SIZES", "16,64,256,1024,4096"), ","))
const PATHS = String.(split(gmf_env_str(
    "GMF_PATHS", "none,gram_point,gram_freeze_step,gram_lookahead,gram_surrogate"), ","))
const MISSION_S = gmf_env_float("GMF_MISSION_S", 300.0)
const GRAVITY = gmf_env_str("GMF_GRAVITY", "l20")
const TIMEOUT_S = gmf_env_float("GMF_TIMEOUT_S", 2400.0)
const REPEATS = gmf_env_int("GMF_REPEATS", 1)
const MIN_FREE_GB = gmf_env_float("GMF_MIN_FREE_GB", 12.0)
const ISOLATE = gmf_env_str("GMF_ISOLATE", "1")
const SOLVE_REPEATS = gmf_env_str("GMF_SOLVE_REPEATS", "1")
const SUFFIX = gmf_env_str("GMF_SUFFIX", "")

"""Memory the kernel could hand out right now, in GB (`MemAvailable`)."""
function available_gb()::Float64
    Sys.islinux() && isfile("/proc/meminfo") || return Sys.free_memory() / (1 << 30)
    for line in eachline("/proc/meminfo")
        startswith(line, "MemAvailable:") || continue
        kb = tryparse(Int, first(split(strip(split(line, ":"; limit=2)[2]))))
        kb === nothing || return kb * 1024 / (1 << 30)
    end
    return Sys.free_memory() / (1 << 30)
end

"""SPACEAGORA_* switches that select one native-GRAM calling method. Mirrors
`ps_constellation_env` in benchmarks/studies/paper_scenarios/common.jl, so a
path here is the same path the S2 scenario measured."""
function path_env(path::String)::Vector{Pair{String, String}}
    base = Pair{String, String}[
        "SPACEAGORA_RHS_EXECUTION_MODE" => "serial",
        "SPACEAGORA_HARMONICS_BATCH_ENABLED" => "0",
        "SPACEAGORA_EFFECTOR_PARALLEL" => "off",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0",
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "0",
        "SPACEAGORA_RHS_CALIBRATE" => "off",
        "SPACEAGORA_INNER_THREAD_BUDGET" => "1",
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "off",
    ]
    if path == "gram_point"
        # Direct native point density: no per-step freeze, no cache. This is the
        # path `_is_native_gram_point_density` is named for.
        append!(base, [
            "SPACEAGORA_DENSITY_FREEZE_PER_STEP" => "0",
            "SPACEAGORA_VACUUM_GRAM_CACHE" => "0",
        ])
    elseif path == "gram_freeze_step"
        append!(base, [
            "SPACEAGORA_DENSITY_FREEZE_PER_STEP" => "1",
            "SPACEAGORA_VACUUM_GRAM_CACHE" => "0",
        ])
    elseif path == "gram_lookahead"
        # Horizon and deviation loosened past what the mission can reach, so only
        # the initial (proven-safe) cache build runs -- the workaround the S2
        # scenario uses for the known multi-satellite native rebuild hang.
        append!(base, [
            "SPACEAGORA_DENSITY_FREEZE_PER_STEP" => "0",
            "SPACEAGORA_VACUUM_GRAM_CACHE" => "1",
            "SPACEAGORA_VACUUM_GRAM_CACHE_NPOINTS" => "20",
            "SPACEAGORA_VACUUM_GRAM_CACHE_HORIZON_S" => string(MISSION_S + 500.0),
            "SPACEAGORA_VACUUM_GRAM_CACHE_DEVIATION_M" => "1e8",
        ])
    end
    return base
end

function run_point(path::String, n::Int, rep::Int, log_dir::String)::Dict{String, String}
    julia_bin = Base.julia_cmd().exec[1]
    cmd = Cmd([julia_bin, "--project=$(GMF_REPO_ROOT)", "--threads=1", GMF_WORKER])
    cmd = addenv(cmd, vcat(path_env(path), [
        "GMF_PATH" => path,
        "GMF_N_SATS" => string(n),
        "GMF_MISSION_S" => string(MISSION_S),
        "GMF_GRAVITY" => GRAVITY,
        "GMF_ISOLATE" => ISOLATE,
        "GMF_SOLVE_REPEATS" => SOLVE_REPEATS,
    ])...)
    log_path = joinpath(log_dir, "$(path)_n$(n)_iso$(ISOLATE)_x$(SOLVE_REPEATS)_r$(rep).log")
    result = Dict{String, String}()
    open(log_path, "w") do log_io
        proc = run(pipeline(cmd; stdout=log_io, stderr=log_io); wait=false)
        deadline = time() + TIMEOUT_S
        while process_running(proc) && time() < deadline
            sleep(1.0)
        end
        if process_running(proc)
            try
                run(`pkill -9 -P $(getpid(proc))`; wait=false)
            catch
            end
            kill(proc, Base.SIGKILL)
            println("  -> TIMEOUT after $(TIMEOUT_S) s (log: $(log_path))")
            result["note"] = "timeout"
            return result
        end
    end
    for line in eachline(log_path)
        startswith(line, "GMF_RESULT ") || continue
        for token in split(chopprefix(line, "GMF_RESULT "))
            kv = split(token, "="; limit=2)
            length(kv) == 2 && (result[kv[1]] = kv[2])
        end
    end
    if isempty(result)
        println("  -> FAILED (no GMF_RESULT line; log: $(log_path))")
        result["note"] = "failed"
    end
    return result
end

const COLUMNS = [
    "hostname", "cpu_threads", "path", "n_sats", "repeat", "mission_s", "gravity", "isolate", "solve_repeats",
    "julia_threads", "ok", "rss_base_mb", "rss_built_mb", "rss_solved_mb", "rss_retained_mb",
    "maxrss_mb", "vmrss_base_mb", "vmrss_retained_mb", "solve_s", "retcode", "note",
]

function main()
    out_dir = joinpath(GMF_DIR, "results")
    log_dir = joinpath(out_dir, "logs")
    mkpath(log_dir)
    host = gethostname()
    csv_path = joinpath(out_dir, "gram_memory_footprint_$(host)$(SUFFIX).csv")

    println("GRAM memory footprint: paths=$(PATHS) sizes=$(SIZES) mission_s=$(MISSION_S) " *
            "gravity=$(GRAVITY) repeats=$(REPEATS)")
    println("Available memory now: $(round(available_gb(); digits=1)) GB")

    rows = Dict{String, String}[]
    for rep in 1:REPEATS, path in PATHS, n in SIZES
        free = available_gb()
        row = Dict{String, String}(
            "hostname" => host,
            "cpu_threads" => string(Sys.CPU_THREADS),
            "path" => path,
            "n_sats" => string(n),
            "repeat" => string(rep),
            "mission_s" => string(MISSION_S),
            "gravity" => GRAVITY,
            "isolate" => ISOLATE,
            "solve_repeats" => SOLVE_REPEATS,
            "julia_threads" => "1",
            "ok" => "false",
            "note" => "",
        )
        if free < MIN_FREE_GB
            println("[$(path) n=$(n) r=$(rep)] SKIPPED: only $(round(free; digits=1)) GB available")
            row["note"] = "skipped_low_memory"
            push!(rows, row)
            continue
        end
        @printf("[%s n=%d r=%d] available=%.1f GB ... ", path, n, rep, free)
        flush(stdout)
        t0 = time()
        result = run_point(path, n, rep, log_dir)
        merge!(row, result)
        @printf("ok=%s rss_retained=%s MB (%.0f s wall)\n",
                get(row, "ok", "false"), get(row, "rss_retained_mb", "-"), time() - t0)
        flush(stdout)
        push!(rows, row)
        # Write after every point: a long grid that is interrupted still leaves
        # everything it did measure on disk.
        write_csv(csv_path, rows)
    end
    println("Wrote $(length(rows)) row(s) to $(csv_path)")
    return nothing
end

function write_csv(path::String, rows::Vector{Dict{String, String}})
    open(path, "w") do io
        println(io, join(COLUMNS, ","))
        for row in rows
            println(io, join([get(row, c, "") for c in COLUMNS], ","))
        end
    end
    return nothing
end

main()
