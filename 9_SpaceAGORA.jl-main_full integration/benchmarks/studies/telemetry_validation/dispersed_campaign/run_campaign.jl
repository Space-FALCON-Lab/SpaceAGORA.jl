# Dispersed Odyssey MarsGRAM campaign: N perturbed members (design B, per-pass
# look-ahead, per-pass reseed, GRAM density scale 1.0, winds nominal) run as ONE
# campaign under the adaptive execution policy (profile R7, route chosen by
# run_monte_carlo(threads=:auto)), or the nominal member alone.
#
#   julia --project=. --threads=32 .../run_campaign.jl --out=DIR --orbits=42 --seeds=101:132
#   julia --project=. .../run_campaign.jl --out=DIR --orbits=42 --member=nominal
#
# --orbits=N+2 runs to the apoapsis after pass N+1, which gives the N+1 apoapses
# and N passes an N-pass analysis needs (analyze_campaign.py; N = 40 by default).
#
# The pool size is the policy's (SPACEAGORA_PERF_PROCS caps it; the remote
# launcher sets it from --threads). Writes DIR/samples.csv (one row per member,
# including each member's own solve time and the dispatcher's per-sample
# elapsed time), DIR/campaign.toml (route, split, wall time), and per-member
# per_orbit.csv / per_pass_r.csv under DIR/<member>/. Each member's outcome is
# classified (member_status.jl): dispatch failure, solver failure, early
# termination before the horizon, or complete; campaign.toml counts each.

include(joinpath(@__DIR__, "..", "common.jl"))
const OPTS = parse_kv_args(copy(ARGS))
const OUT = abspath(OPTS["out"])
const ORBITS = parse(Int, OPTS["orbits"])
const MEMBER = get(OPTS, "member", "campaign")
const SEEDS = let r = split(get(OPTS, "seeds", "101:132"), ":")
    collect(parse(Int, r[1]):parse(Int, r[2]))
end
# The analysis horizon these orbit events reach: N passes for --orbits=N+2.
const HORIZON_PASSES = ORBITS - 2
HORIZON_PASSES >= 1 || error("--orbits=$ORBITS reaches no complete pass; use N+2 for an N-pass horizon")
include(joinpath(@__DIR__, "member_status.jl"))
using .DispersedMemberStatus
mkpath(OUT)

# ── Run-wide environment, set before anything spawns so pool workers inherit it.
# The runner's solver/GRAM settings (as TelemetryVerification's runner applies
# them per solve) plus the perturbation mode. No per-member ENV change happens
# anywhere, which is what makes the members safe on any route.
const PERTURBED = MEMBER != "nominal"
for (k, v) in [
    "SPACEAGORA_WARN_NORMALIZE" => "0",
    "SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0",
    "SPACEAGORA_GRAM_OFFLINE_SURROGATE" => "off",
    "SPACEAGORA_GRAM_STATIC_GRID" => "off",
    "SPACEAGORA_GRAM_TRACK_CACHE" => "off",
    "SPACEAGORA_GRAM_GLOBAL_LOCK" => "on",
    "SPACEAGORA_GRAM_DENSITY_PERTURBATION" => PERTURBED ? "pass" : "off",
    "SPACEAGORA_GRAM_DENSITY_PERTURBATION_PASS_RESEED" => "1",
    "SPACEAGORA_GRAM_DENSITY_PERTURBATION_PASS_DT_S" => "1.0",
    "SPACEAGORA_GRAM_DENSITY_PERTURBATION_LOG" => "results_directory",
    "SPACEAGORA_CAMPAIGN_DISPATCH_TRACE" => "1",
]
    ENV[k] = v
end

load_gramsuite!()
using SpaceAGORA
using SpaceAGORA.TelemetryVerification
using DataFrames, CSV, TOML, Dates, Printf

const TV = SpaceAGORA.TelemetryVerification
const SC = SpaceAGORA.SimulationCampaigns
const SAMPLE_FILE = joinpath(@__DIR__, "odyssey_sample.jl")

const SOLVER_MODE = TV._telemetry_solver_mode()
ENV["SPACEAGORA_SOLVER_MODE"] = SOLVER_MODE
ENV["SPACEAGORA_SOLVER_MAXITERS"] = string(TV._telemetry_solver_maxiters(:full))
ENV["SPACEAGORA_SOLVER_SAVE_EVERYSTEP"] = TV._telemetry_solver_save_env("SPACEAGORA_SOLVER_SAVE_EVERYSTEP", SOLVER_MODE)
ENV["SPACEAGORA_SOLVER_SAVE_ON"] = TV._telemetry_solver_save_env("SPACEAGORA_SOLVER_SAVE_ON", SOLVER_MODE)
# The adaptive execution policy, as the paper's full-budget runs use it (R7).
for (k, v) in SpaceAGORA.ParallelProfiles.profile_env_pairs("R7"; preserve_existing=false)
    ENV[k] = v
end

TV._planet_from_name("mars")                      # furnish kernels on the coordinator
Base.include(Main, SAMPLE_FILE)

# Per-pass seeds must not collide across members (_gram_pass_seed).
let s = [SpaceAGORA.SimulationModel.SimulationCallbacks._gram_pass_seed(m, k) for m in SEEDS for k in 1:ORBITS]
    length(unique(s)) == length(s) || error("per-pass GRAM seeds collide across members")
end

member_dir(seed) = PERTURBED ? joinpath(OUT, "sample_seed$(seed)") : joinpath(OUT, "nominal")

const ENV_SNAPSHOT = Dict(k => v for (k, v) in ENV if startswith(k, "SPACEAGORA_") || k in ("JULIA_NUM_THREADS",))

if MEMBER == "nominal"
    r = Base.invokelatest(Main.OdysseyDispersedSample.run_member, 1001, false, ORBITS, member_dir(1001))
    r = merge(r, (status=member_status(merge(r, (dispatch_success=true,)); horizon_passes=HORIZON_PASSES),))
    CSV.write(joinpath(OUT, "nominal_summary.csv"), DataFrame([r]))
    println(r)
else
    # The sample closure: on a pool worker the module is included on first use
    # (workers bootstrap only SpaceAGORA and GRAMSuite), after furnishing the
    # Mars kernels there.
    sample = let sample_file = SAMPLE_FILE, orbits = ORBITS, out = OUT
        seed -> begin
            if !isdefined(Main, :OdysseyDispersedSample)
                SpaceAGORA.TelemetryVerification._planet_from_name("mars")
                Base.include(Main, sample_file)
            end
            Base.invokelatest(Main.OdysseyDispersedSample.run_member, seed, true, orbits,
                              joinpath(out, "sample_seed$(seed)"))
        end
    end
    probe_args = TV._with_study_settings(
        TV._make_orbit_args(Main.OdysseyDispersedSample.member_cfg(SEEDS[1], true), ORBITS); quick=false)
    features = SC.campaign_route_features(probe_args; samples=length(SEEDS))
    tuning = SpaceAGORA.ParallelProfiles.OuterRouteTuning(trace=true)

    # Provision and warm the pool before the campaign, concurrently, as the
    # paper harness does (parallelization_performance/execution.jl, "Provision
    # and warm the Distributed pool *before* the clock starts"). Without it the
    # process route's own ensure_process_workers! warms each new worker by
    # running the REAL first sample, one worker after another: with ~11-minute
    # members that is ~6 h before the dispatch starts (job
    # 20261002-115416-4125610, killed). The warm member is 2 orbits (one drag
    # pass, so the pass-mode path is compiled) with a seed outside SEEDS, into a
    # scratch directory. The route and split are still the policy's: this only
    # sizes the pool the policy may use; if the policy does not choose the
    # process route the pool simply idles.
    pool_n = parse(Int, get(ENV, "SPACEAGORA_PERF_PROCS", string(length(SEEDS))))
    warm_seed = minimum(SEEDS) - 1
    warm = let sample_file = SAMPLE_FILE, out = joinpath(OUT, "_warmup"), s = warm_seed
        () -> begin
            if !isdefined(Main, :OdysseyDispersedSample)
                SpaceAGORA.TelemetryVerification._planet_from_name("mars")
                Base.include(Main, sample_file)
            end
            Base.invokelatest(Main.OdysseyDispersedSample.run_member, s, true, 2,
                              joinpath(out, "pid$(getpid())"))
            nothing
        end
    end
    pool = SC.campaign_process_pool()
    provision_s = @elapsed (pool_ids = SC.ensure_process_workers!(pool, pool_n))
    warm_s = @elapsed begin
        @sync begin
            Threads.@spawn warm()                             # the coordinator (its own thread), for any local slots
            for w in pool_ids
                @async SpaceAGORA.SimulationCampaigns.Distributed.remotecall_wait(warm, w)
            end
        end
    end
    @printf("prewarm pool=%d provision=%.1f s warm=%.1f s\n", length(pool_ids), provision_s, warm_s)
    flush(stdout)

    started = now(UTC)
    result = SC.run_monte_carlo(sample, SEEDS; threads=:auto, route_features=features, route_tuning=tuning)
    finished = now(UTC)

    rows = map(result.samples) do s
        v = s.success ? s.value : nothing
        base = (index=s.index, seed=s.seed, dispatch_success=s.success, dispatch_elapsed_s=s.elapsed_s,
                error=s.success ? "" : sprint(showerror, s.error))
        v === nothing ? base : merge(base, v)
    end
    df = vcat([DataFrame([r]) for r in rows]...; cols=:union)
    df.status = [member_status(r; horizon_passes=HORIZON_PASSES) for r in rows]
    CSV.write(joinpath(OUT, "samples.csv"), df)
    counts = status_counts(rows; horizon_passes=HORIZON_PASSES)
    execution = member_execution_totals(rows)
    summary = Dict{String, Any}(
        "route" => string(result.route), "consumers" => result.threads, "local_slots" => result.local_slots,
        "wall_s" => result.elapsed_s, "started_utc" => string(started), "finished_utc" => string(finished),
        "n_samples" => length(result.samples), "horizon_passes" => HORIZON_PASSES,
        "n_complete" => counts["complete"], "n_early_terminated" => counts["early_termination"],
        "n_solver_failed" => counts["solver_failure"], "n_dispatch_failed" => counts["dispatch_failure"],
        "sum_member_solve_s" => execution.sum_member_solve_s,
        "sum_dispatch_elapsed_s" => sum(df.dispatch_elapsed_s),
        "distinct_member_pids" => execution.distinct_member_pids,
        "coordinator_pid" => getpid(), "julia_threads" => Threads.nthreads(),
        "seeds" => SEEDS, "orbits_requested" => ORBITS,
        "prewarm_pool_workers" => length(pool_ids), "prewarm_provision_s" => provision_s,
        "prewarm_warm_s" => warm_s, "prewarm_seed" => warm_seed, "prewarm_orbits" => 2,
        "commit" => get(ENV, "SPACEAGORA_RUN_COMMIT", "unknown"), "host" => gethostname(),
        "env" => ENV_SNAPSHOT)
    open(io -> TOML.print(io, summary), joinpath(OUT, "campaign.toml"), "w")
    @printf("campaign route=%s consumers=%d local_slots=%d wall=%.1f s complete=%d early_terminated=%d solver_failed=%d dispatch_failed=%d sum_solve=%.1f s pids=%d\n",
            summary["route"], result.threads, result.local_slots, result.elapsed_s,
            counts["complete"], counts["early_termination"], counts["solver_failure"], counts["dispatch_failure"],
            summary["sum_member_solve_s"], summary["distinct_member_pids"])
end
