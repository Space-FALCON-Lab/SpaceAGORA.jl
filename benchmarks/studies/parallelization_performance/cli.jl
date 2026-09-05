const PPC_STUDY_DIR = @__DIR__
const PPC_REPO_ROOT = normpath(joinpath(PPC_STUDY_DIR, "..", "..", ".."))
const PPC_DEFAULT_OUTDIR = joinpath(PPC_REPO_ROOT, "output", "performance", "parallelization_performance")
const PPC_LAUNCHER = joinpath(PPC_REPO_ROOT, "benchmarks", "studies", "parallelization_performance.jl")

using CSV
using DataFrames
using Dates
using Distributed
using LinearAlgebra
using Random
using SHA
using Sockets
using StaticArrays
using Statistics

Base.@kwdef struct PPCConfig
    profile::String = "smoke"
    outdir::String = PPC_DEFAULT_OUTDIR
    modes::Vector{String} = String[]
    cases::Vector{String} = String[]
    parity_cases::Vector{String} = String[]
    threads::Vector{Int} = Int[]
    repeats::Int = 1
    warmup::Int = 0
    seed::Int = 20260615
    solver_mode::String = "auto_stiff"
    process_workers::Int = 2
    mc_samples::Vector{Int} = Int[]
    # True when mc_samples came from --mc-samples or SPACEAGORA_PPC_MC_SAMPLES
    # rather than from the profile's default ladder. Only the joint_routing
    # family consults it: those cases carry their own sample count (a rung is
    # N spacecraft x S samples and the pair is what holds its total work fixed),
    # so they use it unless the caller asked for something specific.
    mc_samples_explicit::Bool = false
    parity_samples::Int = 128
    cpu_pinning::Vector{Int} = Int[]
    worker::Bool = false
    worker_case::String = ""
    worker_mode::String = ""
    worker_threads::Int = 1
    worker_repeat::Int = 1
    # How many timed repeats the worker runs inside its own process, starting at
    # worker_repeat. One Julia subprocess spends ~80 s on startup plus JIT of the
    # RHS/solver stack before it can time anything (measured: 79 s per row for a
    # case whose solve is 0.016 s), so launching a fresh subprocess per repeat
    # spends that cost N times to collect N samples of the same point. Running the
    # repeats in one process pays it once. The repeats stay independently timed
    # and all of them are post-warm-up, which is what the measurement requires;
    # what changes is only that repeats 2..N inherit a process that has already
    # warmed its caches -- if anything a better model of steady state than repeat
    # 1, which is the sample every previous run of this harness relied on.
    worker_repeats::Int = 1
    worker_seed::Int = 20260615
    worker_mc_samples::Int = 1
    worker_outfile::String = ""
    worker_parity::Bool = false
end

# Copy of `cfg` with the named fields replaced. Base.@kwdef gives a keyword
# constructor but not an update-an-existing-instance one, and the worker needs a
# per-repeat variant of its config (repeat index and seed) without restating the
# other two dozen fields.
function _ppc_with(cfg::PPCConfig; kwargs...)
    overrides = Dict{Symbol, Any}(kwargs)
    return PPCConfig(; (f => get(overrides, f, getfield(cfg, f)) for f in fieldnames(PPCConfig))...)
end

@inline function _ppc_bool(raw::AbstractString)::Bool
    token = lowercase(strip(String(raw)))
    token in ("1", "true", "yes", "on") && return true
    token in ("0", "false", "no", "off") && return false
    throw(ArgumentError("Invalid boolean token '$raw'."))
end

function _ppc_csv(raw::AbstractString)::Vector{String}
    token = strip(String(raw))
    isempty(token) && return String[]
    values = String[]
    seen = Set{String}()
    for part in split(token, ",")
        value = strip(part)
        isempty(value) && continue
        if !(value in seen)
            push!(values, value)
            push!(seen, value)
        end
    end
    return values
end

function _ppc_int_csv(raw::AbstractString)::Vector{Int}
    values = Int[]
    for token in _ppc_csv(raw)
        value = parse(Int, token)
        value > 0 || throw(ArgumentError("Expected positive integer, got $value."))
        push!(values, value)
    end
    return values
end

@inline function _ppc_arg_value(arg::String)::String
    occursin("=", arg) || throw(ArgumentError("Expected --key=value argument, got '$arg'."))
    return split(arg, "=", limit=2)[2]
end

# Parses a taskset-style CPU list, e.g. "0-7,16,20-23", into a flat, ordered
# list of unique CPU ids. Used to reserve a fixed pool of physical cores that
# worker subprocesses get pinned to via `taskset` (see ppc_worker_cmd).
function _ppc_parse_cpu_list(raw::AbstractString)::Vector{Int}
    token = strip(String(raw))
    isempty(token) && return Int[]
    cpus = Int[]
    for part in split(token, ",")
        piece = strip(part)
        isempty(piece) && continue
        if occursin("-", piece)
            bounds = split(piece, "-")
            length(bounds) == 2 || throw(ArgumentError("Invalid CPU range '$piece'."))
            lo, hi = parse(Int, bounds[1]), parse(Int, bounds[2])
            lo <= hi || throw(ArgumentError("Invalid CPU range '$piece': start must be <= end."))
            append!(cpus, lo:hi)
        else
            push!(cpus, parse(Int, piece))
        end
    end
    return unique(cpus)
end

# Physical (not logical) core count. Julia's Sys.CPU_THREADS counts SMT
# siblings, and the RHS hot path is FP/SIMD-bound spherical-harmonics work that
# gains nothing from SMT while contending for the same execution ports: on this
# repo's 12-core/24-thread reference box every constellation workload measured
# regressed 1.5-2x going from 16 to 24 Julia threads, and on a 64-core/128-thread
# box the ladder's own top rung would therefore be its *worst* data point. The
# ladder below is built from this instead of Sys.CPU_THREADS so "max threads"
# means "one thread per physical core".
#
# Override with SPACEAGORA_PPC_PHYSICAL_CORES when the detection is wrong or the
# run is confined to a subset of the machine (e.g. under taskset/cgroups).
function _ppc_physical_core_count()::Int
    override = strip(get(ENV, "SPACEAGORA_PPC_PHYSICAL_CORES", ""))
    if !isempty(override)
        parsed = tryparse(Int, override)
        parsed !== nothing && parsed > 0 && return parsed
        @warn "Ignoring non-positive/unparseable SPACEAGORA_PPC_PHYSICAL_CORES." value=override
    end
    if Sys.islinux()
        try
            # /proc/cpuinfo lists one stanza per logical CPU; a physical core is
            # a unique (physical id, core id) pair, so SMT siblings collapse.
            cores = Set{Tuple{String, String}}()
            package = ""
            core = ""
            for line in eachline("/proc/cpuinfo")
                if startswith(line, "physical id")
                    package = strip(split(line, ":", limit=2)[2])
                elseif startswith(line, "core id")
                    core = strip(split(line, ":", limit=2)[2])
                elseif isempty(strip(line)) && !isempty(core)
                    push!(cores, (package, core))
                    package = ""
                    core = ""
                end
            end
            isempty(core) || push!(cores, (package, core))
            isempty(cores) || return length(cores)
        catch err
            @warn "Physical-core detection failed; falling back to Sys.CPU_THREADS." reason=sprint(showerror, err)
        end
    end
    return max(1, Sys.CPU_THREADS)
end

function _ppc_full_thread_ladder(cpu_threads::Int=_ppc_physical_core_count())::Vector{Int}
    max_threads = max(1, cpu_threads)
    if max_threads < 7
        return collect(1:max_threads)
    elseif max_threads >= 64
        return [1, 2, 4, 8, 16, 32, 64]
    elseif max_threads >= 24
        return unique(sort(Int[
            1,
            2,
            4,
            8,
            max(1, round(Int, max_threads / 2)),
            max(1, round(Int, 2 * max_threads / 3)),
            max_threads
        ]))
    end

    return unique(sort(Int[
        1,
        2,
        4,
        max(1, round(Int, max_threads / 2)),
        max(1, round(Int, 3 * max_threads / 4)),
        max_threads
    ]))
end

function _ppc_defaults(profile::String)::NamedTuple
    if profile == "test"
        return (
            modes=["serial", "outer_threads", "outer_process", "inner_only", "full_smart"],
            cases=["montecarlo_multi_sat"],
            parity_cases=["montecarlo_multi_sat"],
            threads=[max(1, min(2, Sys.CPU_THREADS))],
            repeats=1,
            warmup=0,
            mc_samples=[2],
            parity_samples=16
        )
    elseif profile == "smoke"
        return (
            modes=["serial", "full_smart"],
            cases=["single_inverse_square_vacuum", "gravity_4sat_inverse_square_vacuum", "montecarlo_high_accuracy"],
            parity_cases=["single_harmonics_l20_vacuum", "gravity_16sat_l20_vacuum", "montecarlo_high_accuracy"],
            threads=[1, max(1, min(2, Sys.CPU_THREADS))],
            repeats=1,
            warmup=0,
            mc_samples=[1],
            parity_samples=128
        )
    elseif profile == "full"
        return (
            modes=["serial", "outer_threads", "outer_process", "inner_only", "outer_inner_static", "outer_inner_adaptive", "full_smart"],
            cases=[
                "single_inverse_square_vacuum",
                "gravity_4sat_inverse_square_vacuum",
                "gravity_16sat_inverse_square_vacuum",
                "gravity_64sat_inverse_square_vacuum",
                "gravity_256sat_inverse_square_vacuum",
                "gravity_1024sat_inverse_square_vacuum",
                "gravity_2048sat_inverse_square_vacuum",
                "single_harmonics_l20_vacuum",
                "gravity_4sat_l20_vacuum",
                "gravity_16sat_l20_vacuum",
                "gravity_64sat_l20_vacuum",
                "gravity_256sat_l20_vacuum",
                "gravity_1024sat_l20_vacuum",
                "gravity_2048sat_l20_vacuum",
                "single_harmonics_l50_vacuum",
                "srp_heavy_high_area",
                "articulated_1sat_fullstack",
                "multi_16_aero_surrogate_cached",
                "multi_64_high_fidelity",
                "multi_128_high_fidelity",
                "callback_128_aero_thermal",
                "montecarlo_high_accuracy",
                "montecarlo_multi_sat",
                "montecarlo_mars_aerobraking",
            ],
            parity_cases=["single_harmonics_l20_vacuum", "gravity_16sat_l20_vacuum", "articulated_1sat_fullstack", "multi_16_aero_surrogate_cached", "montecarlo_high_accuracy", "montecarlo_multi_sat", "montecarlo_mars_aerobraking"],
            threads=_ppc_full_thread_ladder(),
            repeats=3,
            warmup=1,
            mc_samples=[1, 4, 16, 64, 256, 1024],
            parity_samples=512
        )
    end
    throw(ArgumentError("Unsupported profile '$profile'. Use test, smoke, or full."))
end

function parse_parallelization_performance_cli(args::Vector{String}=ARGS)::PPCConfig
    profile = lowercase(strip(get(ENV, "SPACEAGORA_PPC_PROFILE", "smoke")))
    outdir = get(ENV, "SPACEAGORA_PPC_OUTDIR", PPC_DEFAULT_OUTDIR)
    modes = _ppc_csv(get(ENV, "SPACEAGORA_PPC_MODES", ""))
    cases = _ppc_csv(get(ENV, "SPACEAGORA_PPC_CASES", ""))
    parity_cases = _ppc_csv(get(ENV, "SPACEAGORA_PPC_PARITY_CASES", ""))
    threads = _ppc_int_csv(get(ENV, "SPACEAGORA_PPC_THREADS", ""))
    repeats = parse(Int, get(ENV, "SPACEAGORA_PPC_REPEATS", "0"))
    # -1 (not 0) is the "unset" sentinel: the fallback below is `warmup < 0`, so
    # defaulting this to 0 silently pinned every worker at zero warm-up runs and
    # made the profile defaults (`full` => 1) dead code. With no warm-up the
    # timed run IS the first run, so every reported wall time was dominated by
    # Julia's JIT compilation of the RHS/solver stack rather than by the
    # simulation: measured 21.4 s first solve vs 0.017 s steady-state solve for
    # gravity_16sat_l20_vacuum, i.e. >99.9% compilation. 0 remains selectable
    # explicitly via --warmup=0 / SPACEAGORA_PPC_WARMUP=0.
    warmup = parse(Int, get(ENV, "SPACEAGORA_PPC_WARMUP", "-1"))
    seed = parse(Int, get(ENV, "SPACEAGORA_PPC_SEED", "20260615"))
    solver_mode = lowercase(strip(get(ENV, "SPACEAGORA_PPC_SOLVER_MODE", "auto_stiff")))
    process_workers = parse(Int, get(ENV, "SPACEAGORA_PPC_PROCESS_WORKERS", "2"))
    mc_samples = _ppc_int_csv(get(ENV, "SPACEAGORA_PPC_MC_SAMPLES", ""))
    parity_samples = parse(Int, get(ENV, "SPACEAGORA_PPC_PARITY_SAMPLES", "0"))
    cpu_pinning = _ppc_parse_cpu_list(get(ENV, "SPACEAGORA_PPC_CPU_LIST", ""))
    worker = false
    worker_case = ""
    worker_mode = ""
    worker_threads = 1
    worker_repeat = 1
    worker_repeats = 1
    worker_seed = seed
    worker_mc_samples = 1
    worker_outfile = ""
    worker_parity = false

    for arg in args
        if arg in ("test", "smoke", "full")
            profile = arg
        elseif arg == "--worker"
            worker = true
        elseif startswith(arg, "--profile=")
            profile = lowercase(strip(_ppc_arg_value(arg)))
        elseif startswith(arg, "--outdir=")
            outdir = _ppc_arg_value(arg)
        elseif startswith(arg, "--modes=")
            modes = _ppc_csv(_ppc_arg_value(arg))
        elseif startswith(arg, "--cases=")
            cases = _ppc_csv(_ppc_arg_value(arg))
        elseif startswith(arg, "--parity-cases=")
            parity_cases = _ppc_csv(_ppc_arg_value(arg))
        elseif startswith(arg, "--threads=")
            threads = _ppc_int_csv(_ppc_arg_value(arg))
        elseif startswith(arg, "--repeats=")
            repeats = parse(Int, _ppc_arg_value(arg))
        elseif startswith(arg, "--warmup=")
            warmup = parse(Int, _ppc_arg_value(arg))
        elseif startswith(arg, "--seed=")
            seed = parse(Int, _ppc_arg_value(arg))
        elseif startswith(arg, "--solver-mode=")
            solver_mode = lowercase(strip(_ppc_arg_value(arg)))
        elseif startswith(arg, "--process-workers=")
            process_workers = parse(Int, _ppc_arg_value(arg))
        elseif startswith(arg, "--mc-samples=")
            mc_samples = _ppc_int_csv(_ppc_arg_value(arg))
        elseif startswith(arg, "--parity-samples=")
            parity_samples = parse(Int, _ppc_arg_value(arg))
        elseif startswith(arg, "--cpu-list=")
            cpu_pinning = _ppc_parse_cpu_list(_ppc_arg_value(arg))
        elseif startswith(arg, "--case=")
            worker_case = _ppc_arg_value(arg)
        elseif startswith(arg, "--mode=")
            worker_mode = _ppc_arg_value(arg)
        elseif startswith(arg, "--thread-count=")
            worker_threads = parse(Int, _ppc_arg_value(arg))
        elseif startswith(arg, "--repeat=")
            worker_repeat = parse(Int, _ppc_arg_value(arg))
        elseif startswith(arg, "--worker-repeats=")
            worker_repeats = parse(Int, _ppc_arg_value(arg))
        elseif startswith(arg, "--worker-seed=")
            worker_seed = parse(Int, _ppc_arg_value(arg))
        elseif startswith(arg, "--worker-mc-samples=")
            worker_mc_samples = parse(Int, _ppc_arg_value(arg))
        elseif startswith(arg, "--outfile=")
            worker_outfile = _ppc_arg_value(arg)
        elseif startswith(arg, "--parity=")
            worker_parity = _ppc_bool(_ppc_arg_value(arg))
        else
            throw(ArgumentError("Unknown argument '$arg'."))
        end
    end

    defaults = _ppc_defaults(profile)
    isempty(modes) && (modes = defaults.modes)
    isempty(cases) && (cases = defaults.cases)
    isempty(parity_cases) && (parity_cases = defaults.parity_cases)
    isempty(threads) && (threads = defaults.threads)
    repeats <= 0 && (repeats = defaults.repeats)
    warmup < 0 && (warmup = defaults.warmup)
    # Record this before the default ladder fills the gap: afterwards the two
    # cases are indistinguishable.
    mc_samples_explicit = !isempty(mc_samples)
    isempty(mc_samples) && (mc_samples = defaults.mc_samples)
    parity_samples <= 0 && (parity_samples = defaults.parity_samples)
    process_workers = max(1, process_workers)

    return PPCConfig(
        profile=profile,
        outdir=abspath(outdir),
        modes=modes,
        cases=cases,
        parity_cases=parity_cases,
        threads=threads,
        repeats=repeats,
        warmup=warmup,
        seed=seed,
        solver_mode=solver_mode,
        process_workers=process_workers,
        mc_samples=mc_samples,
        mc_samples_explicit=mc_samples_explicit,
        parity_samples=parity_samples,
        cpu_pinning=cpu_pinning,
        worker=worker,
        worker_case=worker_case,
        worker_mode=worker_mode,
        worker_threads=worker_threads,
        worker_repeat=worker_repeat,
        worker_repeats=worker_repeats,
        worker_seed=worker_seed,
        worker_mc_samples=worker_mc_samples,
        worker_outfile=worker_outfile,
        worker_parity=worker_parity
    )
end
