const PPB_STUDY_DIR    = @__DIR__
const PPB_DEFAULT_OUTDIR = joinpath(PPC_REPO_ROOT, "output", "performance", "paper_benchmarks")

# ── Phase definition ─────────────────────────────────────────────────────────

# thread_mode controls how the thread ladder is applied for each phase:
#   :full_ladder   — use the full machine-scaled ladder (B1, B2, B3, B7)
#   :max_only      — fix at the machine maximum for a fair profile comparison (B5)
#   :single        — force 1 thread per Julia instance; parallelism comes from
#                    process count only (B4, B8 outer_process runs)
#   :low_high      — 1 and the machine maximum only (B6)
#   :router_ladder — geometric [1, 8, 32] clipped to physical cores and to
#                    PPB_ROUTER_LADDER_MAX_THREADS, the
#                    standard budget axis for the expanded router evaluation
#                    (B9-B14). See _ppb_thread_ladder for why.
#
# worker_ladder: when non-empty, the phase is run once per entry, each time
# varying process_workers to that value.  This is used for B4 to
# produce throughput-vs-workers curves.  Serial and full_smart modes run
# identically regardless of worker count, so their rows will be duplicated
# in the output CSV; the analysis layer should filter on mode when plotting.
#
# budget_grid: when non-empty, the phase is run once per (process_workers,
# threads) pair with the thread ladder pinned to that pair's thread count.  Used
# by B13 to split one fixed total core budget between process and thread
# parallelism without running the full cross product (which would also vary the
# total budget, the one thing that phase has to hold constant).  Mutually
# exclusive with worker_ladder.

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
    budget_grid::Vector{Tuple{Int, Int}} = Tuple{Int, Int}[]
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
# Ceiling on the router ladder's top rung, independent of how many physical
# cores the machine has.
#
# These phases are about which route the router picks, not about how far the
# machine scales -- that is B1/B2/B7's job, and they run the full ladder. Taking
# the router ladder all the way to a 64-core box's top rung buys a data point
# that is mostly a measurement of oversubscription: §9 records every workload
# regressing past the physical-core count, the benchmark machine is shared with
# other users' jobs, and a rung that spends every core is the one most likely to
# be contending with something else rather than measuring routing. 32 keeps the
# ladder spanning a 32x range while leaving headroom on a 64-core host.
const PPB_ROUTER_LADDER_MAX_THREADS = 32

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
# B9-B14 join B7/B8 for the same reason: their whole premise is workloads large
# enough for routing to be observable, and preview's 64-satellite / 16-sample caps
# would either drop every case or shrink them back under the measurability floor
# the phases exist to stay above. Smoke-test them with an explicit
# --phases=B10 --threads=1 run against the `test` profile instead.
const PPB_PREVIEW_SKIP_PHASES =
    Set{String}(["B7", "B8", "B9", "B10", "B11", "B12", "B13", "B14"])

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
        budget_grid   = filter(p -> p[1] * p[2] <= PPB_PREVIEW_MAX_WORKERS, phase.budget_grid),
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
        label = "Small-Workload Control — Router Evaluation Below the Measurability Floor",
        # RETAINED AS A CONTROL, NOT AS THE ROUTER EVALUATION. This phase was
        # point 8's first response to B5's narrowness, and its numbers turned out
        # not to be reportable: every case it draws on resolves in 0.017-2.9 s
        # post-warm-up, so fixed per-solve setup dominates, no route can
        # distinguish itself, and the run of 2026-08-18 measured 1.0x speedup on 9
        # of 11 cases with regret of 0.0-3.3% -- all of it inside the ~16%
        # run-to-run variance §7.4 of the findings doc records for measurements
        # this size. B9-B14 below are the actual expanded evaluation, sized so
        # routing is observable.
        #
        # It is kept, rather than deleted, because "below roughly a 3 s serial
        # baseline no routing profile is distinguishable" is a real and useful
        # boundary on the paper's claim, and this phase is the evidence for it.
        # Its regret rows are flagged `below_noise_floor` in reporting so they are
        # never quoted as router performance.
        #
        # Original rationale, still accurate as a description of what it sweeps:
        # it holds the mode ladder to the routes that matter most
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
        #   - duration:             gravity_16sat_l20_vacuum_longmission (~4x
        #                            mission time) vs. the unmodified
        #                            gravity_16sat_l20_vacuum baseline. The
        #                            cadence half of this axis used to be
        #                            _sparse_output; it measured nothing at all
        #                            (num_steps_to_save never reaches the solver;
        #                            §10 of the findings doc) and is removed. B14
        #                            carries the real cadence ladder, through
        #                            simulation_settings.results / data_rate.
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
        ],
        parity_cases = [
            "multi_4_aero_surrogate_cached",
            "multi_256_high_fidelity",
            "montecarlo_mars_gram_live",
            "multi_8sat_magnetorquer_attitude",
            "gravity_16sat_l20_vacuum_longmission",
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

    # ── B9-B14: the expanded router evaluation (review point 8) ───────────────
    #
    # Point 8 asks for the router to be evaluated across six workload axes and for
    # router regret, (T_selected - T_best_static)/T_best_static, to be reported
    # against the best tested static route. B6 attempted that and produced
    # unreportable numbers (see its own comment); these six phases are the
    # replacement, one axis each, and they differ from B6 in three ways that all
    # matter for whether a regret figure means anything:
    #
    #   1. Sized to be measurable. Every case here targets a serial baseline of
    #      several seconds, following B7's finding that sub-second workloads are
    #      fixed-setup-dominated and flatten every routing profile onto the same
    #      number. Reporting flags any point whose serial baseline lands under 3 s
    #      as below_noise_floor rather than silently quoting it.
    #   2. The full static ladder. B6 omitted inner_only (R2) and
    #      outer_inner_static (R3), so its `best_static` was a minimum over a
    #      subset -- which biases regret *downward*, flattering the router. Every
    #      phase here runs all four static routes (plus outer_process where the
    #      workload is independent), so the bar the adaptive routes are held to is
    #      the real one.
    #   3. A budget axis with more than two points (thread_mode=:router_ladder),
    #      because the best static route is itself budget-dependent: regret
    #      measured at one thread count says nothing about regret at another.
    #
    # Each phase varies exactly one workload property across its case list, with
    # constellation size, mission duration and physics otherwise held fixed, so a
    # regret difference between two rungs is attributable to the axis rather than
    # to problem size. The per-case rationale lives on the ppc_single_config
    # branches in parallelization_performance/cases.jl.

    PPBPhase(
        id    = "B9",
        label = "Router Evaluation — Spacecraft Count",
        # Axis: satellites per RHS call, at fixed physics (L50 harmonics, vacuum)
        # and fixed 1 h duration. One heavy effector and no callbacks, so outer-loop
        # routing is the only thing that can move the number.
        #
        # The ladder deliberately straddles the measurability floor: 64 satellites
        # solve in ~0.15 s serial and 4096 in ~8 s, so the low rungs are the
        # evidence for *where* routing stops being decidable rather than wasted
        # points. Reuses the existing gravity_<N>sat_l50_vacuum_1hr family.
        # 256 and 1024 spacecraft resolve in 0.22 s and 1.26 s serially, under the
        # 3 s floor, so their regret was never quotable and they are not re-run.
        # The evidence they carried -- that routing stops being decidable below
        # the floor -- is established and does not need re-measuring.
        cases = [
            "gravity_4096sat_l50_vacuum_1hr",
        ],
        parity_cases = [
            "gravity_4096sat_l50_vacuum_1hr",
        ],
        modes = [
            "serial", "outer_threads", "inner_only", "outer_inner_static",
            "outer_inner_adaptive", "full_smart",
        ],
        mc_samples  = [1],
        repeats     = 3,
        warmup      = 2,
        thread_mode = :router_ladder,
    ),

    PPBPhase(
        id    = "B10",
        label = "Router Evaluation — Atmosphere and GRAM Usage",
        # Axis: how density is obtained, at fixed N=256, fixed 600 s duration and
        # fixed L50 harmonics. Five rungs of increasing native-library
        # serialisation: none -> analytic exponential -> precomputed GRAM
        # surrogate -> live native GRAM -> live GRAM plus SPICE third-body.
        #
        # This is the load-bearing phase. Review point 7 reports that R5 runs ~25%
        # slower than the best tested route on the GRAM workload and slower than
        # serial, which is the single strongest argument against keeping adaptive
        # routing as a contribution; nothing in the existing harness sweeps GRAM
        # fidelity as an isolated variable, so that finding currently cannot be
        # localised to a cause. 600 s rather than 1 h because that is the duration
        # the heavy_<N>sat_gram_nbody_l50 family is already validated at up to
        # N=256, and live GRAM at 256 satellites is the slowest thing in the
        # catalog.
        cases = [
            # atmo256_vacuum_10min dropped: 0.14 s serial, far below the floor.
            "atmo256_exponential_10min",
            "atmo256_gram_surrogate_10min",
            "atmo256_gram_live_10min",
            "atmo256_gram_live_nbody_10min",
        ],
        parity_cases = [
            "atmo256_exponential_10min",
            "atmo256_gram_live_10min",
        ],
        modes = [
            "serial", "outer_threads", "inner_only", "outer_inner_static",
            "outer_inner_adaptive", "full_smart",
        ],
        mc_samples  = [1],
        repeats     = 3,
        warmup      = 2,
        thread_mode = :router_ladder,
    ),

    PPBPhase(
        id    = "B11",
        label = "Router Evaluation — Number of Force and Actuator Models",
        # Axis: model count, at fixed N=256 and fixed 1 h duration, each rung
        # adding exactly one model to the one before it (harmonics -> +SRP ->
        # +aero -> +third-body -> +attitude propagation -> +actuators).
        #
        # §7.1 of the findings doc identifies the heterogeneous effector tuple as a
        # first-order cost the flat work queue has to schedule around, but no
        # existing case sweeps model count in isolation -- the multi-effector cases
        # each pick one arbitrary stack. e5 and e6 are split so that turning on
        # attitude propagation and adding the actuator/control-callback path are
        # measured separately rather than confounded.
        # Two sub-ladders, because the actuator rung is quadratic in satellite
        # count and no single N puts every rung in a measurable band (measured
        # per 10 s of simulated mission, serial: e1 0.0009 s / e5 0.0157 s /
        # e6 0.583 s at N=32, and e5 0.13 s / e6 33.4 s at N=256 -- a full-length
        # e6 at N=256 would run into the hours).
        #
        #   model count   e1..e5 at N=256 over 1 h, spanning ~3-47 s serial
        #   actuators     e5 vs e6 at N=32 over 600 s, where e6 lands at ~35 s
        #
        # Each sub-ladder is internally single-variable, which is what makes the
        # comparison mean anything; stack32_e5_6dof appears in both roles as the
        # actuator pair's control. It sits under the noise floor at that size and
        # is flagged accordingly -- the measurable half of that pair is e6, which
        # is the rung carrying the axis.
        # stack256_e1_harm (0.17 s) and stack32_e5_6dof (1.24 s) are below the
        # floor and dropped. Note the cost: stack32_e5_6dof was the control for
        # the actuator pair, so the single-variable e5-vs-e6 comparison now rests
        # on the earlier measurement rather than on this run.
        cases = [
            "stack256_e3_aero",
            "stack256_e4_nbody",
            "stack256_e5_6dof",
            "stack32_e6_actuated",
        ],
        parity_cases = [
            "stack256_e4_nbody",
            "stack32_e6_actuated",
        ],
        modes = [
            "serial", "outer_threads", "inner_only", "outer_inner_static",
            "outer_inner_adaptive", "full_smart",
        ],
        mc_samples  = [1],
        repeats     = 3,
        warmup      = 2,
        # Four rungs rather than five: this phase has the largest case list and the
        # most expensive cases in the expanded set.
        thread_mode = :router_ladder,
    ),

    PPBPhase(
        id    = "B12",
        label = "Router Evaluation — Interacting vs. Independent Propagation",
        # Axis: how a fixed amount of work is arranged. interact_<N>sat_1hr
        # propagates N satellites as one coupled simulation; independent_1sat_1hr
        # run at mc_samples=N propagates the same N satellite-hours as N fully
        # independent solves. Same planet, orbit, effectors, density model,
        # tolerances and duration on both sides, so coupling is the only variable.
        #
        # B6 gestured at this axis by putting constellation and Monte Carlo cases
        # in one phase, but those differ in physics, duration and size as well as
        # in coupling, so nothing there isolated coupling. This pair is the direct
        # test of the manuscript's third contribution -- that the router
        # distinguishes shared-state multi-spacecraft propagation from
        # process-isolable independent campaigns -- at matched total work.
        #
        # outer_process joins the mode ladder (it is only meaningful on the
        # independent side, and its presence there is the point); inner_only drops
        # out, since with no outer work to contend with it duplicates
        # outer_inner_static on the interacting side.
        cases = [
            # interact_16sat_1hr dropped: 1.38 s serial, below the floor.
            "interact_64sat_1hr",
            "interact_256sat_1hr",
            "independent_1sat_1hr",
        ],
        parity_cases = [
            "interact_64sat_1hr",
            "independent_1sat_1hr",
        ],
        modes = [
            "serial", "outer_threads", "outer_process", "outer_inner_static",
            "outer_inner_adaptive", "full_smart",
        ],
        # Only the independent case consumes these (ppc_run_controller runs
        # non-Monte-Carlo cases once regardless), so each entry is one more
        # independent-side point against the same three interacting cases. 64 is
        # the matched partner of interact_64sat_1hr and carries the comparison;
        # 16 and 256 bracket it so the independent side has its own size trend.
        # mc=64 gives a 2.63 s baseline, below the floor; only 256 is retained.
        mc_samples  = [256],
        repeats     = 3,
        warmup      = 2,
        thread_mode = :router_ladder,
    ),

    PPBPhase(
        id    = "B13",
        label = "Router Evaluation — Thread vs. Process Budget Split",
        # Axis: where a fixed total core budget is spent. Every entry in the grid
        # below multiplies out to 32 cores, split six ways between process workers
        # and threads per worker, so the sweep varies the *shape* of the
        # parallelism and not its amount. 32 rather than the machine's 64 for the
        # same reason the router ladder stops there -- see
        # PPB_ROUTER_LADDER_MAX_THREADS.
        #
        # Nothing in the harness tests this today: B4 and B8 both pin
        # thread_mode=:single and vary workers alone, and every constellation phase
        # runs one process. The nested-split case is precisely what
        # SPACEAGORA_OUTER_PARALLEL_ACTIVE / SPACEAGORA_INNER_THREAD_BUDGET exist
        # to arbitrate (see §3 of the findings doc, where asserting outer_active
        # unconditionally pinned the RHS to one worker), which makes it the
        # configuration the router is most likely to get wrong and the least
        # covered by evidence.
        #
        # montecarlo_heavy_aerobraking because at ~1.05 s per sample the trials
        # dominate dispatch overhead, and because it is the only case that can be
        # arranged either way -- 64 independent samples can be spread across
        # processes, threads, or any mix of the two.
        cases        = ["montecarlo_heavy_aerobraking"],
        parity_cases = ["montecarlo_heavy_aerobraking"],
        modes        = ["outer_process", "outer_threads", "outer_inner_adaptive", "full_smart"],
        mc_samples   = [64],
        repeats      = 3,
        warmup       = 2,
        # Stops at 32 workers rather than 64: each Distributed worker is a full
        # Julia process with SpaceAGORA and its Mars/SPICE state loaded (~2 GB),
        # and the target box had ~140 GB free at sizing time, so a 64-worker pool
        # would run within a few GB of the limit for the sake of one more point.
        # (32, 2) already brackets the process-heavy end of the axis.
        budget_grid  = [(1, 32), (2, 16), (4, 8), (8, 4), (16, 2), (32, 1)],
    ),

    PPBPhase(
        id    = "B14",
        label = "Router Evaluation — Mission Duration and Output Cadence",
        # Two sub-ladders at N=1024, L50, vacuum, sharing the 1 h rung:
        #
        #   duration  15 min / 1 h / 6 h at no trajectory output
        #   cadence   at 1 h: no output / saveat 60 s / 10 s / 1 s
        #
        # The cadence half is new work rather than a re-run. B6's cadence case
        # varied num_steps_to_save, which nothing outside the telemetry-
        # verification code reads -- the solve runs save_everystep with no saveat,
        # so that case was byte-identical work to its own baseline (§10 of the
        # findings doc). The real mechanism is simulation_settings.results gating
        # get_data_saving_callback's SavingCallback(...; saveat=data_rate), and
        # this harness had output switched off entirely. See ppc_build_config.
        #
        # At 1024 satellites the saving callback's per-satellite snapshot is a
        # genuinely serial spine, so this ladder measures an Amdahl term the router
        # has no visibility into at all -- which is the result either way it comes
        # out.
        # cadence_1024sat_none is the same configuration as
        # gravity_1024sat_l50_vacuum_1hr (1 h, L50, 1024 satellites, no output),
        # so the duration ladder's middle rung doubles as the cadence ladder's
        # zero point rather than being run twice.
        # The 15 min and 1 h rungs are 0.37 s and 1.39 s serially, both below the
        # floor, and are dropped. That leaves the duration ladder with a single
        # measurable point, so this phase now carries the cadence axis only: at
        # 1024 spacecraft in vacuum, mission duration is not separable from noise
        # until 6 h, which is itself the result for that axis.
        cases = [
            "heavy_1024sat_l50_6hr",
            "cadence_1024sat_10s",
            "cadence_1024sat_1s",
        ],
        parity_cases = [
            "cadence_1024sat_10s",
        ],
        modes = [
            "serial", "outer_threads", "inner_only", "outer_inner_static",
            "outer_inner_adaptive", "full_smart",
        ],
        mc_samples  = [1],
        repeats     = 3,
        warmup      = 2,
        thread_mode = :router_ladder,
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
