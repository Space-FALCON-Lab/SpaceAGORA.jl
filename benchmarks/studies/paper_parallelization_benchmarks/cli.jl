const PPB_STUDY_DIR    = @__DIR__
const PPB_DEFAULT_OUTDIR = joinpath(PPC_REPO_ROOT, "output", "performance", "paper_benchmarks")

# ── Phase definition ─────────────────────────────────────────────────────────

# thread_mode controls how the thread ladder is applied for each phase:
#   :full_ladder  — use the full machine-scaled ladder (B1, B2, B3)
#   :max_only     — fix at the machine maximum for a fair profile comparison (B5)
#   :single       — force 1 thread per Julia instance; parallelism comes from
#                   process count only (B4 outer_process runs)
#
# worker_ladder: when non-empty, the phase is run once per entry, each time
# varying process_workers to that value.  This is used for B4 to
# produce throughput-vs-workers curves.  Serial and full_smart modes run
# identically regardless of worker count, so their rows will be duplicated
# in the output CSV; the analysis layer should filter on mode when plotting.

Base.@kwdef struct PPBPhase
    id::String
    label::String
    cases::Vector{String}
    parity_cases::Vector{String}
    modes::Vector{String}
    mc_samples::Vector{Int}
    repeats::Int
    warmup::Int
    thread_mode::Symbol        = :full_ladder
    worker_ladder::Vector{Int} = Int[]
end

Base.@kwdef struct PPBConfig
    phases::Vector{String}  = String[]
    outdir::String          = PPB_DEFAULT_OUTDIR
    threads::Vector{Int}    = Int[]
    process_workers::Int    = 32
    seed::Int               = 20260615
    solver_mode::String     = "auto_stiff"
    cpu_pinning::Vector{Int} = Int[]
    dry_run::Bool           = false
    preview::Bool           = false
end

# Preview mode: caps N_sat at 64, MC samples at 16, workers at 4, repeats at 2.
const PPB_PREVIEW_MAX_N_SAT   = 64
const PPB_PREVIEW_MAX_SAMPLES = 16
const PPB_PREVIEW_MAX_WORKERS = 4
const PPB_PREVIEW_REPEATS     = 2
const PPB_PREVIEW_WARMUP      = 1
const PPB_PREVIEW_SKIP_PHASES = Set{String}()

function _ppb_preview_phase(phase::PPBPhase)::PPBPhase
    cases = filter(c -> _ppb_n_sat(c) <= PPB_PREVIEW_MAX_N_SAT, phase.cases)
    isempty(cases) && (cases = phase.cases[1:min(2, end)])
    parity = filter(c -> _ppb_n_sat(c) <= PPB_PREVIEW_MAX_N_SAT, phase.parity_cases)
    isempty(parity) && (parity = phase.parity_cases[1:min(1, end)])
    samples = filter(s -> s <= PPB_PREVIEW_MAX_SAMPLES, phase.mc_samples)
    isempty(samples) && (samples = [1])
    workers = filter(w -> w <= PPB_PREVIEW_MAX_WORKERS, phase.worker_ladder)
    isempty(workers) && !isempty(phase.worker_ladder) && (workers = [1, 2])
    return PPBPhase(
        id            = phase.id,
        label         = phase.label * " [preview]",
        cases         = cases,
        parity_cases  = parity,
        modes         = phase.modes,
        mc_samples    = samples,
        repeats       = PPB_PREVIEW_REPEATS,
        warmup        = PPB_PREVIEW_WARMUP,
        thread_mode   = phase.thread_mode,
        worker_ladder = workers,
    )
end

# ── Phase catalog ─────────────────────────────────────────────────────────────

const PAPER_BENCHMARK_PHASES = PPBPhase[

    PPBPhase(
        id    = "B1",
        label = "Constellation Scaling — Gravity-Only Baseline",
        cases = [
            "single_inverse_square_vacuum",
            "gravity_4sat_inverse_square_vacuum",
            "gravity_16sat_inverse_square_vacuum",
            "gravity_64sat_inverse_square_vacuum",
            "gravity_256sat_inverse_square_vacuum",
            "gravity_1024sat_inverse_square_vacuum",
        ],
        parity_cases = [
            "single_inverse_square_vacuum",
            "gravity_16sat_inverse_square_vacuum",
        ],
        modes      = ["serial", "outer_threads", "full_smart"],
        mc_samples = [1],
        repeats    = 5,
        warmup     = 2,
        thread_mode = :full_ladder,
    ),

    PPBPhase(
        id    = "B2",
        label = "Constellation Scaling — L20 Harmonics",
        cases = [
            "single_harmonics_l20_vacuum",
            "gravity_4sat_l20_vacuum",
            "gravity_16sat_l20_vacuum",
            "gravity_64sat_l20_vacuum",
            "gravity_256sat_l20_vacuum",
            "gravity_1024sat_l20_vacuum",
        ],
        parity_cases = [
            "single_harmonics_l20_vacuum",
            "gravity_16sat_l20_vacuum",
        ],
        modes      = ["serial", "outer_threads", "outer_inner_adaptive", "full_smart"],
        mc_samples = [1],
        repeats    = 5,
        warmup     = 2,
        thread_mode = :full_ladder,
    ),

    PPBPhase(
        id    = "B3",
        label = "Constellation Scaling — Atmosphere and GRAMTrackCache Surrogate",
        cases = [
            "multi_16_aero_surrogate_cached",
            "multi_64_high_fidelity",
        ],
        parity_cases = [
            "multi_16_aero_surrogate_cached",
        ],
        modes      = ["serial", "outer_threads", "full_smart"],
        mc_samples = [1],
        repeats    = 5,
        warmup     = 2,
        thread_mode = :full_ladder,
    ),

    PPBPhase(
        id    = "B4",
        label = "Monte Carlo Throughput Scaling",
        cases = [
            "montecarlo_mars_aerobraking",
            "montecarlo_high_accuracy",
        ],
        parity_cases = [
            "montecarlo_mars_aerobraking",
            "montecarlo_high_accuracy",
        ],
        modes      = ["serial", "outer_process", "full_smart"],
        mc_samples = [1, 4, 16, 64, 256, 1024],
        repeats    = 3,
        warmup     = 1,
        # 1 thread per process: parallelism comes from process count.
        # 256 GB RAM allows up to 32 GRAM+SPICE processes safely.
        thread_mode   = :single,
        worker_ladder = [1, 2, 4, 8, 16, 32],
    ),

    PPBPhase(
        id    = "B5",
        label = "Routing Profile Comparison R0–R5",
        cases = [
            "gravity_16sat_l20_vacuum",
            "multi_16_aero_surrogate_cached",
            "articulated_1sat_fullstack",
        ],
        parity_cases = [
            "gravity_16sat_l20_vacuum",
            "multi_16_aero_surrogate_cached",
            "articulated_1sat_fullstack",
        ],
        modes = [
            "serial",
            "outer_threads",
            "inner_only",
            "outer_inner_static",
            "outer_inner_adaptive",
            "full_smart",
        ],
        mc_samples  = [1],
        repeats     = 5,
        warmup      = 2,
        # Fix thread count at machine maximum so all profiles are compared at
        # the same total thread budget.
        thread_mode = :max_only,
    ),

]

# ── CLI parsing ───────────────────────────────────────────────────────────────

function ppb_parse_cli(args::Vector{String}=ARGS)::PPBConfig
    phases          = _ppc_csv(get(ENV, "SPACEAGORA_PPB_PHASES", ""))
    outdir          = get(ENV, "SPACEAGORA_PPB_OUTDIR", PPB_DEFAULT_OUTDIR)
    threads         = _ppc_int_csv(get(ENV, "SPACEAGORA_PPB_THREADS", ""))
    process_workers = parse(Int, get(ENV, "SPACEAGORA_PPB_PROCESS_WORKERS", "32"))
    seed            = parse(Int, get(ENV, "SPACEAGORA_PPB_SEED", "20260615"))
    solver_mode     = lowercase(strip(get(ENV, "SPACEAGORA_PPB_SOLVER_MODE", "auto_stiff")))
    cpu_pinning     = _ppc_parse_cpu_list(get(ENV, "SPACEAGORA_PPB_CPU_LIST", ""))
    dry_run         = _ppc_bool(get(ENV, "SPACEAGORA_PPB_DRY_RUN", "0"))
    preview         = _ppc_bool(get(ENV, "SPACEAGORA_PPB_PREVIEW", "0"))

    valid_phases = Set(p.id for p in PAPER_BENCHMARK_PHASES)
    for arg in args
        if arg in valid_phases
            push!(phases, arg)
        elseif startswith(arg, "--phases=")
            append!(phases, _ppc_csv(_ppc_arg_value(arg)))
        elseif startswith(arg, "--outdir=")
            outdir = _ppc_arg_value(arg)
        elseif startswith(arg, "--threads=")
            threads = _ppc_int_csv(_ppc_arg_value(arg))
        elseif startswith(arg, "--process-workers=")
            process_workers = parse(Int, _ppc_arg_value(arg))
        elseif startswith(arg, "--seed=")
            seed = parse(Int, _ppc_arg_value(arg))
        elseif startswith(arg, "--solver-mode=")
            solver_mode = lowercase(strip(_ppc_arg_value(arg)))
        elseif startswith(arg, "--cpu-list=")
            cpu_pinning = _ppc_parse_cpu_list(_ppc_arg_value(arg))
        elseif arg == "--dry-run"
            dry_run = true
        elseif arg == "--preview"
            preview = true
        else
            throw(ArgumentError("Unknown argument '$arg'. Valid phases: $(join(sort(collect(valid_phases)), ", "))."))
        end
    end

    unique!(phases)
    unknown = [p for p in phases if p ∉ valid_phases]
    isempty(unknown) || throw(ArgumentError("Unknown phase(s): $(join(unknown, ", "))."))

    return PPBConfig(
        phases          = phases,
        outdir          = abspath(outdir),
        threads         = threads,
        process_workers = max(1, process_workers),
        seed            = seed,
        solver_mode     = solver_mode,
        cpu_pinning     = cpu_pinning,
        dry_run         = dry_run,
        preview         = preview,
    )
end
