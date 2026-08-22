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
    mc_samples_max::Union{Nothing, Int} = nothing
    seed::Int               = 20260615
    solver_mode::String     = "auto_stiff"
    cpu_pinning::Vector{Int} = Int[]
    dry_run::Bool           = false
    preview::Bool           = false
    resume::String          = ""
end

# Preview mode: caps N_sat at 64, MC samples at 16, workers at 4, repeats at 2.
const PPB_PREVIEW_MAX_N_SAT   = 64
const PPB_PREVIEW_MAX_SAMPLES = 16
const PPB_PREVIEW_MAX_WORKERS = 4
const PPB_PREVIEW_REPEATS     = 2
const PPB_PREVIEW_WARMUP      = 1
# B7/B8 are skipped under --preview. Preview exists to smoke-test the phase
# structure on a laptop and caps N_sat at 64, but B7's whole premise is
# workloads large enough for scaling to be observable -- and its case filter
# would drop every case (all are >64 satellites) and then fall back to running
# the two largest anyway. B8's samples are seconds each by design. Run them on
# the benchmark machine or not at all.
const PPB_PREVIEW_SKIP_PHASES = Set{String}(["B7", "B8"])

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

    PPBPhase(
        id    = "B6",
        label = "Expanded Router Evaluation — Spacecraft Count, Atmosphere/GRAM, Actuators, Propagation Mode, Duration/Cadence",
        # B5 evaluates the full routing-profile ladder on 3 fixed workloads at
        # a single thread budget. This phase is point 8's response to that
        # narrowness: it holds the mode ladder to the routes that matter most
        # for router *regret* (serial/outer_threads/outer_process as static
        # baselines, outer_inner_adaptive/full_smart as the routes being
        # judged against them -- see PPB_STATIC_MODES/PPB_ADAPTIVE_MODES in
        # reporting.jl) and instead sweeps the workload axes B5 held fixed:
        #   - spacecraft count:      4 / 16 / 64 / 256 (many_sat_high_fidelity ladder)
        #   - atmosphere/GRAM usage: montecarlo_mars_gram_live (native
        #                            GRAMAtmosphereModel, not the analytic
        #                            ExponentialAtmosphereModel every other
        #                            atmosphere case here uses). multi_16_gram_live
        #                            (the interacting/16-satellite counterpart) is
        #                            defined in cases.jl but deliberately NOT run
        #                            here: it triggers an unbounded per-call memory
        #                            leak in the vendored GRAMSuite.jl native
        #                            binding (reproduced with solver_mode=tsit5 too,
        #                            so it isn't the auto_stiff/Rodas5P path) --
        #                            pre-existing, out of scope for this phase,
        #                            being investigated separately. Do not add it
        #                            back to this case list until that's resolved.
        #   - force/actuator count:  articulated_1sat_fullstack (harmonics+aero+
        #                            thermal+attitude) and multi_8sat_magnetorquer_attitude
        #                            (the only case in this whole catalog with a
        #                            real control_effector -- every other case has
        #                            control_effectors=())
        #   - interacting vs.
        #     independent:          constellation cases vs. montecarlo_mars_gram_live/
        #                           montecarlo_multi_sat
        #   - thread/process budget: thread_mode=:low_high (1 and machine-max,
        #                            rather than B5's single fixed max) plus
        #                           outer_process in the mode ladder
        #   - duration/cadence:     gravity_16sat_l20_vacuum_longmission (~4x
        #                            mission time) and _sparse_output (15x fewer
        #                            saved steps), both vs. the unmodified
        #                            gravity_16sat_l20_vacuum baseline
        cases = [
            "gravity_16sat_l20_vacuum",
            "multi_4_aero_surrogate_cached",
            "multi_16_aero_surrogate_cached",
            "multi_64_high_fidelity",
            "multi_256_high_fidelity",
            "montecarlo_mars_gram_live",
            "montecarlo_multi_sat",
            "articulated_1sat_fullstack",
            "multi_8sat_magnetorquer_attitude",
            "gravity_16sat_l20_vacuum_longmission",
            "gravity_16sat_l20_vacuum_sparse_output",
        ],
        parity_cases = [
            "multi_4_aero_surrogate_cached",
            "multi_256_high_fidelity",
            "montecarlo_mars_gram_live",
            "multi_8sat_magnetorquer_attitude",
            "gravity_16sat_l20_vacuum_longmission",
            "gravity_16sat_l20_vacuum_sparse_output",
        ],
        modes = [
            "serial",
            "outer_threads",
            "outer_process",
            "outer_inner_adaptive",
            "full_smart",
        ],
        mc_samples  = [1, 8],
        repeats     = 3,
        warmup      = 1,
        thread_mode = :low_high,
    ),

    PPBPhase(
        id    = "B7",
        label = "Heavy Constellation Scaling — Thread Ladder on Workloads Large Enough to Scale",
        # B1-B6 all draw from cases whose post-warm-up solve time is a fraction of
        # a second on a modern machine (0.017 s for gravity_16sat_l20_vacuum,
        # 0.18 s for multi_64_high_fidelity). At that size the measurement is
        # dominated by fixed per-solve setup and the serial spine of the
        # integrator, so the thread ladder is flat by construction and no routing
        # profile can distinguish itself. This phase runs the same ladder against
        # workloads whose serial baseline is ~8 s and whose parallelisable RHS is
        # the majority of that, which is the only regime where the scaling claim
        # is testable.
        #
        # The case list is deliberately a contrast, not a sweep:
        #   heavy_1024sat_l50_6hr       single heavy effector, no callbacks --
        #                               the cleanest read on outer-loop scaling.
        #   heavy_4096sat_l50_1hr       same physics, 4x the satellites per RHS
        #                               call. Reference measurement on 12 physical
        #                               cores: 1024 sats saturate at ~4 workers
        #                               (3.9x) while 4096 keeps scaling to 12+
        #                               (5.1x), so the pair localises where
        #                               per-RHS fork/join overhead stops being
        #                               amortised instead of reporting one number.
        #   heavy_1024sat_fullstack_1hr heterogeneous effector mix plus density
        #                               and thermal callbacks at scale -- the case
        #                               where inner/callback parallelism has work.
        #   heavy_256sat_coupled6dof_2hr attitude propagation on for every
        #                               satellite; the only many-satellite case in
        #                               the catalog that is actually 6-DOF.
        cases = [
            "heavy_1024sat_l50_6hr",
            "heavy_4096sat_l50_1hr",
            "heavy_1024sat_fullstack_1hr",
            "heavy_256sat_coupled6dof_2hr",
        ],
        parity_cases = [
            "heavy_1024sat_l50_6hr",
            "heavy_1024sat_fullstack_1hr",
            "heavy_256sat_coupled6dof_2hr",
        ],
        modes = [
            "serial",
            "outer_threads",
            "outer_inner_adaptive",
            "full_smart",
        ],
        mc_samples  = [1],
        repeats     = 3,
        warmup      = 1,
        thread_mode = :full_ladder,
    ),

    PPBPhase(
        id    = "B8",
        label = "Heavy Monte Carlo Process Throughput",
        # B4's process-throughput curve uses samples that take milliseconds each,
        # so its wall time is mostly per-worker process startup and JIT rather
        # than the trials, and adding workers can make the curve go the wrong way.
        # This phase uses a 12x longer aerobraking arc: ~1.05 s of integration per
        # sample (measured), against which pmap dispatch is ~1% rather than the
        # dominant term. (The harness change that provisions and warms the
        # Distributed pool before the clock starts removes the startup component;
        # this removes what remained of the signal-to-noise problem.)
        #
        # The sample ladder stops at 64 on purpose. Every worker-ladder entry is a
        # separate controller run that re-executes the full mode x sample x repeat
        # grid, and serial/full_smart produce identical numbers at every worker
        # count (both run the batch in-process at thread_mode=:single), so their
        # cost is paid seven times over. At [16, 64] that is ~13 min per entry;
        # adding 256 quadruples it for no additional shape in the curve, since 64
        # samples across 64 workers is already exactly one dispatch round.
        cases        = ["montecarlo_heavy_aerobraking"],
        parity_cases = ["montecarlo_heavy_aerobraking"],
        modes        = ["serial", "outer_process", "full_smart"],
        mc_samples   = [16, 64],
        repeats      = 3,
        warmup       = 1,
        thread_mode   = :single,
        worker_ladder = [1, 2, 4, 8, 16, 32, 64],
    ),

]

# ── CLI parsing ───────────────────────────────────────────────────────────────

function ppb_parse_cli(args::Vector{String}=ARGS)::PPBConfig
    phases          = _ppc_csv(get(ENV, "SPACEAGORA_PPB_PHASES", ""))
    outdir          = get(ENV, "SPACEAGORA_PPB_OUTDIR", PPB_DEFAULT_OUTDIR)
    threads         = _ppc_int_csv(get(ENV, "SPACEAGORA_PPB_THREADS", ""))
    process_workers = parse(Int, get(ENV, "SPACEAGORA_PPB_PROCESS_WORKERS", "32"))
    mc_samples_max  = let raw = strip(get(ENV, "SPACEAGORA_PPB_MC_SAMPLES_MAX", ""))
        isempty(raw) ? nothing : parse(Int, raw)
    end
    seed            = parse(Int, get(ENV, "SPACEAGORA_PPB_SEED", "20260615"))
    solver_mode     = lowercase(strip(get(ENV, "SPACEAGORA_PPB_SOLVER_MODE", "auto_stiff")))
    cpu_pinning     = _ppc_parse_cpu_list(get(ENV, "SPACEAGORA_PPB_CPU_LIST", ""))
    dry_run         = _ppc_bool(get(ENV, "SPACEAGORA_PPB_DRY_RUN", "0"))
    preview         = _ppc_bool(get(ENV, "SPACEAGORA_PPB_PREVIEW", "0"))
    resume          = get(ENV, "SPACEAGORA_PPB_RESUME", "")

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
        elseif startswith(arg, "--mc-samples-max=")
            mc_samples_max = parse(Int, _ppc_arg_value(arg))
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
        elseif startswith(arg, "--resume=")
            resume = _ppc_arg_value(arg)
        else
            throw(ArgumentError("Unknown argument '$arg'. Valid phases: $(join(sort(collect(valid_phases)), ", "))."))
        end
    end

    unique!(phases)
    unknown = [p for p in phases if p ∉ valid_phases]
    isempty(unknown) || throw(ArgumentError("Unknown phase(s): $(join(unknown, ", "))."))
    (mc_samples_max === nothing || mc_samples_max >= 1) ||
        throw(ArgumentError("--mc-samples-max must be >= 1, got $(mc_samples_max)."))
    resume = strip(resume)
    (isempty(resume) || isdir(resume)) ||
        throw(ArgumentError("--resume path does not exist: $(resume)"))

    return PPBConfig(
        phases          = phases,
        outdir          = abspath(outdir),
        threads         = threads,
        process_workers = max(1, process_workers),
        mc_samples_max  = mc_samples_max,
        seed            = seed,
        solver_mode     = solver_mode,
        cpu_pinning     = cpu_pinning,
        dry_run         = dry_run,
        preview         = preview,
        resume          = isempty(resume) ? "" : abspath(resume),
    )
end
