# ── RHS execution plan auto-calibration ─────────────────────────────────────
#
# When multiple Julia threads are available, a short pre-solve sweep times a
# small set of candidate execution plans and pins the fastest one for the
# duration of the solve.  Results are persisted to a machine-specific TOML file
# so subsequent runs reuse the result without re-sweeping.
#
# The winning plan is written to SharedBuffers.rhs_plan_override[].
# _rhs_execution_plan in setup.jl reads it first and returns immediately when
# it is non-nothing, bypassing all heuristic routing logic.
#
# Environment controls:
#   SPACEAGORA_RHS_CALIBRATE=auto     (default) calibrate if no cached result
#   SPACEAGORA_RHS_CALIBRATE=force    always re-calibrate
#   SPACEAGORA_RHS_CALIBRATE=off      never calibrate
#   SPACEAGORA_RHS_CALIBRATION_PATH   override the TOML file path
#   SPACEAGORA_RHS_CALIBRATE_N_WARMUP number of discarded warm-up calls (default 5)
#   SPACEAGORA_RHS_CALIBRATE_N_TIMED  number of timed calls per candidate (default 10)
#   SPACEAGORA_RHS_CALIBRATE_SCHEDULERS  static,dynamic (default) -- which inner
#                                     schedulers the flat ladder is swept over
#
# The sweep covers (scheduler x allotment) for the flat route, not allotment
# alone. The scheduler used to come from SPACEAGORA_PARALLEL_POLICY_INNER_SCHEDULER
# and was therefore a profile constant the sweep held fixed while it optimised
# around it: R5 declares `dynamic`, which costs an atomic RMW per chunk on work
# items that are uniform by construction, and calibration could not discover
# that `static` was faster for the same allotment.

const _CALIB_MACHINE_LABEL    = Ref{String}("")
const _rhs_calib_cache        = Dict{String, Dict{String, Any}}()
const _rhs_calib_lock         = ReentrantLock()
const _rhs_calib_loaded       = Ref{Bool}(false)

@inline function _rhs_calibration_mode()::Symbol
    raw = lowercase(strip(_engine_env_get("SPACEAGORA_RHS_CALIBRATE", "auto")))
    raw in ("auto", "")    && return :auto
    raw in ("force", "always") && return :force
    raw in ("off", "false", "0", "none") && return :off
    throw(ArgumentError(
        "SPACEAGORA_RHS_CALIBRATE must be auto, force, or off; got '$raw'"
    ))
end

@inline function _rhs_calibrate_n_warmup()::Int
    n = try parse(Int, strip(_engine_env_get("SPACEAGORA_RHS_CALIBRATE_N_WARMUP", "5"))) catch; 5 end
    return max(1, n)
end

@inline function _rhs_calibrate_n_timed()::Int
    n = try parse(Int, strip(_engine_env_get("SPACEAGORA_RHS_CALIBRATE_N_TIMED", "10"))) catch; 10 end
    return max(1, n)
end

# Which schedulers the flat ladder is swept over. Restricting this to a single
# value reproduces the pre-change sweep (one scheduler, allotment only) and is
# the escape hatch if the doubled sweep cost ever matters on a short solve.
# How much faster a swept plan must be before it displaces the heuristic.
#
# Not zero: the sweep's readings carry real spread, and overriding on a
# difference inside that spread trades a plan that adapts per call for one that
# cannot, on no evidence. Ten percent is above the repeat-to-repeat spread
# measured on this harness and well below the margins that motivate calibration
# in the first place -- the wins it exists to capture are 3-5x, not 3%.
@inline function _rhs_calibrate_override_margin()::Float64
    raw = strip(_engine_env_get("SPACEAGORA_RHS_CALIBRATE_OVERRIDE_MARGIN", "0.10"))
    v = try
        parse(Float64, raw)
    catch
        throw(ArgumentError("SPACEAGORA_RHS_CALIBRATE_OVERRIDE_MARGIN must be a float, got '$raw'"))
    end
    return clamp(v, 0.0, 0.9)
end

# Persist a "the heuristic won" verdict, so it is reached once per machine
# rather than re-derived once per solve.
#
# The no-regret floor's retain-the-heuristic outcome used to return `nothing`
# and store nothing, and `_rhs_calib_lookup` had no representation for it -- so
# the next solve on the same signature saw a cache MISS and ran the entire
# successive-halving sweep again, to reach the same conclusion, forever. On the
# shapes where the floor actually fires (small constellations, three-effector
# queues -- the light workloads, where a fixed pre-solve cost is largest
# relative to the solve) that is roughly 110 discarded RHS evaluations per run
# that the floor guarantees cannot buy anything.
#
# Measured with the paired probe in its cold-calibration mode (which drops the
# in-process memo before each sample, modelling the fresh-process regime the
# benchmark harness and every real user actually run in): light_16_harm was
# +30.7% and light_64_aero +14.9% cold against warm, and neither stores an entry
# even under SPACEAGORA_RHS_CALIBRATE=force.
#
# Off restores the previous behaviour, which is what the A/B probe reverts on
# the B side to isolate this mechanism.
@inline function _rhs_calibrate_cache_heuristic()::Bool
    return SimulationModel.ParallelPolicy.parse_bool_env(
        "SPACEAGORA_RHS_CALIBRATE_CACHE_HEURISTIC", true
    )
end

# Sentinel `mode` for a cached retain-the-heuristic verdict. Deliberately not a
# plan: nothing is pinned, rhs_plan_override stays unset and the runtime
# heuristic runs, exactly as it does after a live sweep reaches the same answer.
# Older builds parse this row, fail both plan-mode comparisons in
# `_rhs_calib_lookup` and fall through to a miss, so the file stays
# backward-compatible.
const _CALIB_HEURISTIC_MODE = "heuristic"

function _rhs_calibrate_schedulers()::Vector{Symbol}
    raw = lowercase(strip(_engine_env_get("SPACEAGORA_RHS_CALIBRATE_SCHEDULERS", "static,dynamic")))
    out = Symbol[]
    for token in split(raw, ',')
        t = strip(token)
        if t in ("static", "strided")
            :static in out || push!(out, :static)
        elseif t == "dynamic"
            :dynamic in out || push!(out, :dynamic)
        elseif !isempty(t)
            throw(ArgumentError(
                "SPACEAGORA_RHS_CALIBRATE_SCHEDULERS entries must be static or dynamic; got '$t'"
            ))
        end
    end
    return isempty(out) ? Symbol[:static, :dynamic] : out
end

# ── Machine fingerprint ───────────────────────────────────────────────────────

function _calib_machine_label()::String
    if isempty(_CALIB_MACHINE_LABEL[])
        label_override = strip(_engine_env_get("SPACEAGORA_PERF_MACHINE_LABEL", ""))
        if !isempty(label_override)
            _CALIB_MACHINE_LABEL[] = SimulationModel.ParallelPolicy._safe_token(label_override)
        else
            cpu_info = Sys.cpu_info()
            cpu_model = isempty(cpu_info) ? "unknown" : String(cpu_info[1].model)
            raw = "$(cpu_model)_$(Sys.CPU_THREADS)"
            _CALIB_MACHINE_LABEL[] = bytes2hex(SHA.sha256(codeunits(raw))[1:8])
        end
    end
    return _CALIB_MACHINE_LABEL[]
end

@inline function _calib_sat_bucket(n::Int)::String
    n <= 1   && return "1"
    n <= 4   && return "2_4"
    n <= 8   && return "5_8"
    n <= 16  && return "9_16"
    n <= 32  && return "17_32"
    n <= 64  && return "33_64"
    n <= 128 && return "65_128"
    n <= 256 && return "129_256"
    return "257p"
end

# ── Signature ─────────────────────────────────────────────────────────────────

# Stable fingerprint of WHICH effectors are in the queue, not just how many.
#
# The signature keyed on `effs=<count>` plus a has-harmonics flag, so a
# 2-effector harmonics+drag constellation and a 2-effector inverse-square+SRP
# one shared a cache entry and replayed each other's pinned plan. effs=2 is the
# most populated bucket in this machine's own state file and it is the light
# workload bucket, so the collision lands where it is least affordable. Type
# names are stable across runs and cheap to hash; the concrete values inside an
# effector are deliberately excluded, since a plan is a function of the work
# SHAPE, not of a drag coefficient.
function _rhs_calib_effector_token(dynamic_effectors)::String
    names = sort!([string(nameof(typeof(e))) for e in dynamic_effectors])
    return bytes2hex(SHA.sha256(codeunits(join(names, ",")))[1:4])
end

# The density model's identity, because it decides both the cost of a pass and
# whether the cost model is willing to rank plans at all: model_in_cost_domain
# rejects GRAMAtmosphereModel, whose process-wide lock has no representation as
# a sum of per-unit rates. Plain type name rather than a hash -- there is only
# ever one, and a readable cache file is worth more than four bytes.
_rhs_calib_density_token(density_model)::String = string(nameof(typeof(density_model)))

function _rhs_calib_signature(p, dynamic_effectors, density_model)::String
    budget      = SimulationModel.ParallelPolicy.effective_inner_thread_budget()
    active_sats = count(identity, p.is_active)
    n_eff       = length(dynamic_effectors)
    has_harmonics = n_eff == 1 &&
        dynamic_effectors[1] isa SimulationModel.GravitationalHarmonicsModel
    # An enclosing outer split changes which plan the heuristic will pick --
    # _rhs_execution_plan_uncached clamps the flat routes' allotment to 1 under
    # `outer_serialized` specifically to avoid nesting a thread split inside an
    # already-blocked outer worker. A signature blind to it lets an entry
    # calibrated in a single-simulation run be replayed under a campaign's outer
    # split at full width (the nested-oversubscription hazard) and an entry
    # calibrated under an outer split be replayed single-simulation at
    # allotment 1 (throughput left on the table). Budget does not separate the
    # two: the benchmark harness asserts SPACEAGORA_OUTER_PARALLEL_ACTIVE=1 even
    # for single-simulation constellation cases, where there is no outer split
    # at all.
    outer_active = SimulationModel.ParallelPolicy.outer_parallel_active()
    return join([
        # v3: entries predating the no-regret floor must not be replayed.
        #
        # The floor only runs inside a sweep. A cached entry short-circuits the
        # sweep entirely (`_rhs_calib_lookup` returns and calibration returns),
        # so a machine carrying pre-floor entries would keep pinning plans the
        # floor exists to reject -- indefinitely, since nothing re-sweeps a
        # signature that already has an answer. Invalidating them is the only
        # way the fix reaches a machine that has already been calibrated.
        #
        # v2 pinned the inner scheduler, which v1 left to ENV; a v1 entry
        # replayed under v2 would have silently reinstated whatever scheduler the
        # current profile declared.
        # v4: v3 entries key on effector COUNT and are blind to the outer split,
        # so replaying one is exactly the collision the eff= and outer= terms
        # below exist to prevent -- a v3 entry has no way to say which of the
        # colliding shapes produced it. Nothing re-sweeps a signature that
        # already has an answer, so invalidating them is the only way the fix
        # reaches an already-calibrated machine.
        # v5: v4 entries are blind to the DENSITY MODEL, and that is not a
        # near-miss collision -- it silently defeats the cost model's own
        # abstention.
        #
        # select_plan declines when counts.in_domain is false, and
        # model_in_cost_domain clears it for GRAMAtmosphereModel precisely
        # because the native lock's cost is superlinear in concurrency. But a
        # cached entry short-circuits calibration before any of that runs, so a
        # plan calibrated against a cheap analytic atmosphere is replayed
        # verbatim against live GRAM -- the one workload the model refuses to
        # rank. Measured on B10, where three cases share an effector set and
        # differ only in density model: atmo256_exponential_10min calibrates to
        # satellite_batch and caches it, then atmo256_gram_surrogate_10min
        # replays it at +79% regret and atmo256_gram_live_10min at +119%, both
        # slower than serial. atmo256_gram_live_nbody_10min escapes only by
        # accident -- its extra effectors change the eff= token, so it misses
        # the cache, recomputes, and correctly abstains.
        #
        # Nothing re-sweeps a signature that already has an answer, so
        # invalidating v4 is the only way this reaches an already-calibrated
        # machine.
        # v6: v5 entries carry no evidence -- just a verdict -- so the
        # confirm-on-hit in `_calibrate_rhs_plan_if_needed!` has nothing to
        # re-measure against, and their `elapsed_mean_ns` was produced under the
        # pre-2026-08-30 asymmetric comparison (candidates interleaved, the
        # heuristic contiguous) so it is not comparable to a contiguous confirm
        # anyway. Nothing re-sweeps a signature that already has an answer, so
        # invalidating is the only way this reaches an already-calibrated
        # machine -- the same argument v3, v4 and v5 each made.
        "v6",
        "machine=$(_calib_machine_label())",
        "budget=$(budget)",
        "sats=$(_calib_sat_bucket(active_sats))",
        "effs=$(n_eff)",
        "harm=$(has_harmonics ? "1" : "0")",
        "eff=$(_rhs_calib_effector_token(dynamic_effectors))",
        "dens=$(_rhs_calib_density_token(density_model))",
        "outer=$(outer_active ? "1" : "0")",
    ], "|")
end

# ── Persistence ───────────────────────────────────────────────────────────────

function _rhs_calib_path()::String
    override = strip(_engine_env_get("SPACEAGORA_RHS_CALIBRATION_PATH", ""))
    if !isempty(override)
        return normpath(isabspath(override) ? override : joinpath(pwd(), override))
    end
    return normpath(joinpath(
        pwd(), "output", "parallel_policy_state",
        "rhs_calibration_$(_calib_machine_label()).toml"
    ))
end

function _rhs_calib_load!()::Nothing
    lock(_rhs_calib_lock) do
        _rhs_calib_loaded[] && return nothing
        _rhs_calib_loaded[] = true
        path = _rhs_calib_path()
        isfile(path) || return nothing
        parsed = try TOML.parsefile(path) catch; return nothing end
        rows = get(parsed, "calibrations", Any[])
        rows isa AbstractVector || return nothing
        for row in rows
            row isa AbstractDict || continue
            sig = get(row, "signature", "")
            isempty(sig) && continue
            entry = Dict{String, Any}(
                "mode"           => String(get(row, "mode", "")),
                "allotment"      => Int(get(row, "allotment", 1)),
                "scheduler"      => String(get(row, "scheduler", "auto")),
                "elapsed_mean_ns"=> Float64(get(row, "elapsed_mean_ns", 0.0)),
                "solve_ns"       => Float64(get(row, "solve_ns", 0.0)),
            )
            _rhs_calib_cache[sig] = entry
        end
        return nothing
    end
end

function _rhs_calib_save!()::Nothing
    lock(_rhs_calib_lock) do
        isempty(_rhs_calib_cache) && return nothing
        path = _rhs_calib_path()
        rows = Dict{String, Any}[]
        for sig in sort!(collect(keys(_rhs_calib_cache)))
            e = _rhs_calib_cache[sig]
            row = Dict{String, Any}(
                "signature"       => sig,
                "mode"            => get(e, "mode", ""),
                "allotment"       => Int(get(e, "allotment", 1)),
                "scheduler"       => get(e, "scheduler", "auto"),
                "elapsed_mean_ns" => Float64(get(e, "elapsed_mean_ns", 0.0)),
                # How long the SOLVE this verdict was formed for actually took.
                # The gate below reads it; 0.0 means "never measured", which is
                # treated as "sweep", not as "short".
                "solve_ns"        => Float64(get(e, "solve_ns", 0.0)),
            )
            push!(rows, row)
        end
        payload = Dict{String, Any}(
            "schema_version" => 1,
            "calibrations"   => rows,
        )
        try
            mkpath(dirname(path))
            tmp = path * ".tmp"
            open(tmp, "w") do io
                TOML.print(io, payload)
            end
            mv(tmp, path; force=true)
        catch e
            @warn "RHS calibration: failed to save to $(path)" exception=e
        end
        return nothing
    end
end

# Three outcomes, not two: `nothing` is a MISS (sweep), a NamedTuple is a plan to
# pin, and `:heuristic` is a cached decision to pin nothing and not sweep.
function _rhs_calib_lookup(sig::String)::Union{Nothing, Symbol, NamedTuple}
    _rhs_calib_load!()
    entry = lock(_rhs_calib_lock) do
        get(_rhs_calib_cache, sig, nothing)
    end
    entry === nothing && return nothing
    mode_str  = get(entry, "mode", "")
    allotment = max(1, Int(get(entry, "allotment", 1)))
    scheduler = Symbol(get(entry, "scheduler", "auto"))
    if mode_str == _CALIB_HEURISTIC_MODE
        return _rhs_calibrate_cache_heuristic() ? :heuristic : nothing
    end
    mode_str == "satellite_batch"                  && return _make_calib_satellite_batch_plan()
    mode_str == "flat_constellation_effector_queue" && return _make_calib_flat_plan(allotment, scheduler)
    return nothing
end

function _rhs_calib_store_heuristic!(sig::String, heuristic_ns::Float64 = 0.0)::Nothing
    _rhs_calib_load!()
    lock(_rhs_calib_lock) do
        entry = Dict{String, Any}(
            "mode"            => _CALIB_HEURISTIC_MODE,
            "allotment"       => 1,
            "scheduler"       => "auto",
            "elapsed_mean_ns" => heuristic_ns,
        )
        _rhs_calib_cache[sig] = entry
    end
    return nothing
end

function _rhs_calib_store!(sig::String, plan, elapsed_mean_ns::Float64)::Nothing
    # Settle the lazy one-time disk load first: otherwise a later first lookup
    # would load persisted entries over fresher in-process stores.
    _rhs_calib_load!()
    lock(_rhs_calib_lock) do
        entry = Dict{String, Any}(
            "mode"            => String(plan.mode),
            "allotment"       => Int(plan.allotment),
            "scheduler"       => String(plan.scheduler),
            "elapsed_mean_ns" => elapsed_mean_ns,
        )
        _rhs_calib_cache[sig] = entry
    end
    return nothing
end

# ── Plan constructors ─────────────────────────────────────────────────────────

const _CALIB_SERIAL_EFFECTOR_DECISION = (
    use_threads  = false,
    allotment    = 1,
    mode         = :off,
    policy_applied = false,
)

@inline function _make_calib_satellite_batch_plan()
    return (
        mode           = :satellite_batch,
        allotment      = 1,
        scheduler      = :static,
        dominant_axis  = :satellite,
        policy_applied = true,
        effector_decision = _CALIB_SERIAL_EFFECTOR_DECISION,
    )
end

@inline function _make_calib_flat_plan(allotment::Int, scheduler::Symbol=:dynamic)
    return (
        mode           = :flat_constellation_effector_queue,
        allotment      = max(1, allotment),
        # Honoured by the RHS dispatch sites now, rather than being overridden
        # by SPACEAGORA_PARALLEL_POLICY_INNER_SCHEDULER; :auto defers to that
        # env var, which is what a legacy cache entry without this field means.
        scheduler      = (scheduler === :static || scheduler === :dynamic) ? scheduler : :auto,
        dominant_axis  = :flat_effector,
        policy_applied = true,
        effector_decision = _CALIB_SERIAL_EFFECTOR_DECISION,
    )
end

# ── Candidate enumeration ─────────────────────────────────────────────────────

function _rhs_plan_candidates(p, dynamic_effectors)
    budget      = SimulationModel.ParallelPolicy.effective_inner_thread_budget()
    active_sats = count(identity, p.is_active)
    min_sats_floor = SimulationModel.ParallelPolicy.harmonics_batch_spin_barrier_enabled() ?
        1 : _rhs_harmonics_batch_min_sats_per_worker()
    viable_workers = fld(active_sats, max(1, min_sats_floor))

    candidates = Any[_make_calib_satellite_batch_plan()]

    if viable_workers >= 2 && _rhs_flat_supported(dynamic_effectors)
        # Ladder is geometric in the THREAD BUDGET, not in viable_workers.
        #
        # viable_workers is a SIMD batch-sizing quantity (active_sats /
        # min_sats_per_worker) and is routinely far larger than the budget: at
        # 1024 satellites with the default floor of 4 it is 256, against a
        # 12-thread budget. The old ladder was built from it -- 2,
        # viable_workers/2, viable_workers -- and every entry was then clamped by
        # min(a, budget), so [2, 128, 256] collapsed to [2, 12, 12]. The sweep
        # therefore only ever compared width 2 against the full budget and could
        # not discover anything in between.
        #
        # That matters because the optimum is in between. On this repo's
        # 12-physical-core reference box, every multi-effector constellation case
        # measured peaks at 4 workers and is *slower than serial* at 12
        # (heavy_1024sat_fullstack_1hr: 13.7 s at 4 threads vs 25.1 s at 12,
        # against 24.0 s serial; same shape at 256 satellites and in 6-DOF). A
        # ladder that skips 4 cannot find it.
        #
        # 1 is included so the sweep can conclude "do not thread this at all",
        # which is the right answer for workloads whose curve inverts.
        max_workers = max(1, min(budget, viable_workers))
        allotments = Int[1]
        a = 2
        while a < max_workers
            push!(allotments, a)
            a *= 2
        end
        push!(allotments, max_workers)
        sort!(unique!(allotments))
        # (scheduler x allotment), not allotment alone. The two axes interact:
        # dynamic's atomic-per-chunk cost scales with the item count while its
        # load-balancing benefit scales with per-item cost variance, so the best
        # allotment under one scheduler is not the best under the other. At
        # allotment 1 the dispatch degenerates to a serial loop and the
        # scheduler is unobservable, so it is only crossed for allotment >= 2.
        for scheduler in _rhs_calibrate_schedulers()
            for a in allotments
                a >= 2 || continue
                push!(candidates, _make_calib_flat_plan(a, scheduler))
            end
        end
        if 1 in allotments
            push!(candidates, _make_calib_flat_plan(1, :static))
        end
    end

    return candidates
end

# ── Sweep ─────────────────────────────────────────────────────────────────────

# Which statistic reduces a candidate's timed samples to one score.
#
#   trimmed  (default) drop the slowest sample, then take the MEAN of the rest
#   min                the previous behaviour, minimum over all samples
#
# Minimum was wrong, and specifically wrong in a direction that biases the
# sweep. Its justification -- "interference is strictly additive and one-sided,
# a preempted call is longer and never shorter" -- is correct for EXTERNAL
# interference and false for a parallel plan's own straggler variance. That
# variance is not noise to be filtered out: it is part of what the plan costs,
# it grows with allotment, and the solve pays it on every call. Scoring on the
# minimum therefore reads a wide plan at its luckiest and a narrow plan at very
# nearly its average, and systematically prefers the wide one.
#
# Trimming the single slowest sample keeps the one-sided-interference filter the
# minimum was reaching for -- a 2 ms preemption anywhere in the block is still
# discarded -- while the mean over what remains is the statistic the solve
# actually pays. At fewer than three samples there is nothing to trim and this
# degenerates to the minimum, which is the old behaviour.
#
# MEASURED EFFECT: none that this harness can resolve. Paired probe, 21 pairs,
# cold, isolated stores, trimmed against min: gravity_4096 +1.6%, heavy_1024
# -0.9%, interact_256 -1.2%, light_64_aero -0.8%, all p >= 0.38. Verdict
# stability over 16 forced sweeps on gravity_4096 moved from 12/16 to 13/16
# retain-the-heuristic, which is inside its own spread. The bias argued above is
# real in principle and too small to detect at this sample size; this is kept
# because it is the more defensible statistic at no cost, NOT because it was
# shown to be faster. Revert the default to `min` freely if it ever complicates
# anything.
# Off restores contiguous per-candidate blocks, which is what the B side of the
# isolating A/B reverts to.
# Off keeps the sweep's oversized flat partials buffer for the rest of the solve,
# which is the previous behaviour and what the isolating A/B reverts on B.
# How close two swept candidates must be before the narrower one wins.
#
# DEFAULT OFF (0.0), and that is a measured decision rather than caution. The
# argument for it is sound -- see the tie-break block in _run_rhs_sweep! -- but
# on the one path with enough near-ties to exercise it (gravity_4096 under an
# enclosing outer split, where the pinned allotment is unstable at 12/8/4 even
# with interleaving) a 5% margin over 16 forced sweeps did not stabilise the
# verdict and WIDENED the outcome set, pulling in satellite_batch results the
# pure-minimum chooser never produced:
#
#   margin 0.05   flat(12)x5  flat(8)x4  flat(4)x4  satellite_batch(1)x3
#   margin 0.00   flat(12)x7  flat(8)x5  flat(4)x4
#
# So the mechanism ships and the default does not. Set it if a workload is found
# where near-ties are both frequent and genuinely equal; do not set it on the
# reasoning alone, which is what this comment previously amounted to.
@inline function _rhs_calibrate_tie_margin()::Float64
    raw = strip(_engine_env_get("SPACEAGORA_RHS_CALIBRATE_TIE_MARGIN", "0.0"))
    v = try
        parse(Float64, raw)
    catch
        throw(ArgumentError("SPACEAGORA_RHS_CALIBRATE_TIE_MARGIN must be a float, got '$raw'"))
    end
    return clamp(v, 0.0, 0.5)
end

@inline function _rhs_calibrate_release_scratch()::Bool
    return SimulationModel.ParallelPolicy.parse_bool_env(
        "SPACEAGORA_RHS_CALIBRATE_RELEASE_SCRATCH", true
    )
end

# DEFAULT OFF since 2026-08-30, and the reasoning that put it on is still
# sound -- it was applied to one side of a two-sided comparison.
#
# Interleaving does cancel drift, which is why the mechanism stays. But
# `_measure_round!` interleaves the swept candidates while the no-regret floor
# measures the heuristic as a CONTIGUOUS block, so the accept rule compares two
# different measurement regimes. Measured on interact_256sat_1hr, same process,
# contiguous block against round-robin passes:
#
#     flat(static,  8)    0.225 -> 0.874 ms/call   (3.9x)
#     flat(dynamic, 8)    0.230 -> 0.644           (2.8x)
#     satellite_batch     0.314 -> 0.329           (1.05x)
#
# The penalty lands on the flat plans and not on satellite_batch, so the
# comparison was biased against pinning any flat plan -- in the direction of
# retaining the heuristic, on evidence that is an artefact of how the candidates
# were sampled rather than of how they run.
#
# Contiguous is also simply the right model: the solve runs ONE plan for ~16,000
# consecutive calls. Round-robin measures a regime that never occurs.
#
# Set this to 1 on a machine with real drift over the sweep's duration (thermal
# ramp, a noisy neighbour) where the ordering bias is the larger error; the
# comparison is then at least symmetric because the heuristic is measured the
# same way -- see `_measure_heuristic!`.
@inline function _rhs_calibrate_interleave()::Bool
    return SimulationModel.ParallelPolicy.parse_bool_env(
        "SPACEAGORA_RHS_CALIBRATE_INTERLEAVE", false
    )
end

@inline function _rhs_calibrate_score_statistic()::Symbol
    raw = lowercase(strip(_engine_env_get("SPACEAGORA_RHS_CALIBRATE_SCORE", "trimmed")))
    raw in ("trimmed", "", "trimmed_mean") && return :trimmed
    raw in ("min", "minimum", "best")       && return :min
    throw(ArgumentError(
        "SPACEAGORA_RHS_CALIBRATE_SCORE must be trimmed or min; got '$raw'"
    ))
end

@inline function _rhs_reduce_samples(samples::Vector{Float64}, statistic::Symbol)::Float64
    isempty(samples) && return Inf
    (statistic === :min || length(samples) < 3) && return minimum(samples)
    sort!(samples)
    kept = @view samples[1:(length(samples) - 1)]   # drop the slowest
    return sum(kept) / length(kept)
end

function _rhs_sweep_measure!(p, du, u0, candidate, reps::Int)::Float64
    p.shared_buffers.rhs_plan_override[] = candidate
    statistic = _rhs_calibrate_score_statistic()
    samples = Float64[]
    sizehint!(samples, max(1, reps))
    try
        for _ in 1:max(1, reps)
            t0 = time_ns()
            spacecraft_dynamics!(du, u0, p, 0.0)
            push!(samples, Float64(time_ns() - t0))
        end
        return _rhs_reduce_samples(samples, statistic)
    catch e
        @warn "RHS calibration: candidate skipped due to error" mode=candidate.mode allotment=candidate.allotment exception=e
        return Inf
    finally
        p.shared_buffers.rhs_plan_override[] = nothing
    end
end

function _run_rhs_sweep!(p, u0, dynamic_effectors, verbose::Bool, args = nothing)
    n_warmup   = _rhs_calibrate_n_warmup()
    n_timed    = _rhs_calibrate_n_timed()
    candidates = _rhs_plan_candidates(p, dynamic_effectors)

    # NO-REGRET FLOOR: the runtime heuristic competes as a candidate.
    #
    # Without this the sweep can only pick the best of the plans it happens to
    # enumerate, and it enumerates satellite_batch plus a flat ladder -- never
    # the plan `_rhs_execution_plan_uncached` would have chosen on its own. So
    # calibration could confidently pin something worse than doing nothing, and
    # nothing in the measurement could detect it.
    #
    # That was not hypothetical. On gravity_4096sat_l50_vacuum_1hr and
    # heavy_1024sat_l50_6hr at eight threads, both the static profiles and R5
    # record policy_decisions_total = 0 -- no policy consultation happens at all
    # -- and the entire 10-12% gap is that R5 pins a calibrated
    # flat(allotment=8) while the static profiles run the heuristic. The
    # heuristic is better on those workloads because it re-derives the plan
    # against the live satellite count and outer-split state on every call,
    # where a pinned plan is fixed at whatever the pre-solve probe saw.
    #
    # The outer-route bandit has always had this property: default_outer_route's
    # answer is in its candidate set and ranked first. Plan selection did not.
    heuristic = nothing
    if args !== nothing
        heuristic = try
            _rhs_execution_plan_uncached(args, p, dynamic_effectors, length(args.dynamics_model.spacecraft))
        catch
            nothing
        end
    end

    # Nothing to compare: only the baseline candidate exists.
    length(candidates) <= 1 && return nothing, 0.0, :aborted, 0.0, nothing

    du = zero(u0)

    # SUCCESSIVE HALVING rather than an exhaustive equal-budget sweep.
    #
    # The old sweep gave every candidate the full warm-up plus timed block --
    # 10 candidates x (5 + 15) = 200 RHS evaluations, all discarded -- and spent
    # exactly as much measuring an obviously terrible plan as the eventual
    # winner. Most candidates are separable after one call; only the last two or
    # three need precision.
    #
    # So each round times the survivors, keeps the better half, and doubles the
    # repetitions for the next round. Bad candidates are eliminated cheaply and
    # the budget concentrates where the decision is actually close:
    #
    #     10 x 1 + 5 x 2 + 3 x 4 + 2 x 8 + 1 x 15  =  63 timed calls
    #
    # against 150, for the same final precision on the winner. This is the
    # standard fixed-budget bandit allocation, and it needs no cost model -- it
    # is strictly cheaper than what it replaces with the same answer whenever the
    # eliminated candidates were genuinely worse, which one call is enough to
    # establish for all but near-ties.
    #
    # Warm-up is per-sweep plus one call per candidate, not the full block per
    # candidate. The expensive part of warm-up is JIT and thread-pool spin-up,
    # which is shared: candidates differ by allotment and scheduler, both runtime
    # values, so only the two plan *modes* need separate compilation. That takes
    # warm-up from 50 calls to 15.
    p.shared_buffers.rhs_plan_override[] = first(candidates)
    try
        for _ in 1:n_warmup
            spacecraft_dynamics!(du, u0, p, 0.0)
        end
    catch e
        @warn "RHS calibration: warm-up failed; skipping calibration." exception=e
        p.shared_buffers.rhs_plan_override[] = nothing
        return nothing, 0.0, :aborted, 0.0, nothing
    finally
        p.shared_buffers.rhs_plan_override[] = nothing
    end
    for candidate in candidates
        _rhs_sweep_measure!(p, du, u0, candidate, 1)
    end

    if verbose
        println(
            "[SpaceAGORA] RHS calibration: successive halving over " *
            "$(length(candidates)) candidates (budget $(n_timed) timed calls)"
        )
    end

    # Interleave the round's samples across candidates instead of running each
    # candidate's whole block contiguously.
    #
    # This file's own measurement methodology says never to compare in blocks:
    # the same thermal-callback comparison read +11.1% block-ordered and +1.7%
    # paired, and the 11.1% was entirely ordering. That lesson was applied to
    # the external A/B instrument and not to this sweep, which is the internal
    # one and makes an irreversible decision. Each candidate's reps ran
    # back-to-back, satellite_batch always ran first in every round, and any
    # drift across the round -- frequency ramp, thermal throttle, a GC pause, a
    # background process -- landed on whichever candidate happened to be
    # executing.
    #
    # Round-robin over reps fixes it the same way the paired probe does: every
    # candidate is sampled once per pass, so slow drift scales all of them
    # together and cancels in the comparison rather than accumulating on one.
    # The sweep order also reverses on alternate passes, so a systematic
    # first-position or last-position bias cancels across passes instead of
    # always favouring the candidate that happens to be enumerated first.
    #
    # Same total number of RHS calls; only their order changes.
    @inline function _measure_round!(cands::Vector, reps::Int, out::Dict{Any, Float64})
        if !_rhs_calibrate_interleave() || length(cands) <= 1
            for c in cands
                out[c] = _rhs_sweep_measure!(p, du, u0, c, reps)
            end
            return reps * length(cands)
        end
        statistic = _rhs_calibrate_score_statistic()
        samples = Dict{Any, Vector{Float64}}(c => Float64[] for c in cands)
        failed = Set{Any}()
        for pass in 1:max(1, reps)
            order = isodd(pass) ? cands : reverse(cands)
            for c in order
                c in failed && continue
                # One rep at a time, so the plan override is re-armed per sample.
                # That is the same Ref write the block version does once per
                # block; at ~10 ns against an RHS call it is not a measurable
                # addition to what is being timed.
                v = _rhs_sweep_measure!(p, du, u0, c, 1)
                isfinite(v) ? push!(samples[c], v) : push!(failed, c)
            end
        end
        for c in cands
            out[c] = c in failed ? Inf : _rhs_reduce_samples(samples[c], statistic)
        end
        return reps * length(cands)
    end

    survivors = collect(candidates)
    scores = Dict{Any, Float64}()
    # Two repetitions in the first round, not one.
    #
    # Elimination is irreversible, so the round that discards the most
    # candidates is the one that can least afford a bad reading. A single timing
    # is a min-of-1, which is just the raw sample and fully exposed to the
    # one-sided interference the minimum is supposed to filter -- so the true
    # best could be eliminated in round one on a single unlucky call. Two
    # samples costs ten extra calls out of roughly sixty and removes the
    # single-sample failure mode; candidates surviving to later rounds get
    # geometrically more evidence anyway.
    reps = min(2, n_timed)
    total_calls = 0
    round_index = 0

    # EARLY-OUT ON CANDIDATE SPREAD.
    #
    # The heuristic joins round one instead of being measured only at the very
    # end, so the sweep can find out cheaply whether the decision it is about to
    # spend ~100 ms resolving carries any value at all.
    #
    # Solve duration bounds how much a wrong verdict can COST; it does not
    # predict whether calibration will FIND anything, and those come apart
    # across machines. Measured on interact_256sat_1hr full profile, 20 samples:
    # on space-falcon-1 at 8 threads every candidate sits within a few percent of
    # the heuristic, so sweeping is pure overhead and costs 17.5 % (374 -> 453
    # ms/solve); on trx50 at 32 threads satellite_batch runs 0.108 ms/call
    # against flat plans at 0.5-1.4, so sweeping pays for itself and saves 2.4 %.
    # A gate on solve length gets one of those right and the other slightly
    # wrong. Candidate spread is the quantity that actually differs.
    #
    # When the heuristic is already inside the accept margin after round one,
    # the remaining halving rounds cannot change the verdict -- they only refine
    # the ranking among plans that will lose to it -- so the sweep jumps straight
    # to the final paired round, which is the measurement that actually decides.
    # ~40 calls instead of ~112. The decision is NOT taken on round-one data;
    # round one only decides whether to keep refining.
    while true
        round_index += 1
        empty!(scores)
        total_calls += _measure_round!(survivors, reps, scores)
        viable = [c for c in survivors if isfinite(get(scores, c, Inf))]
        isempty(viable) && return nothing, 0.0, :aborted, 0.0, nothing

        if verbose
            for candidate in sort(viable; by = c -> scores[c])
                label = candidate.mode == :satellite_batch ?
                    "satellite_batch                 " :
                    "flat($(rpad(candidate.scheduler, 7)) allotment=$(lpad(candidate.allotment, 3)))"
                println("  [x$(lpad(reps, 2))] $(label) → $(round(scores[candidate] / 1e6, digits=3)) ms/call")
            end
        end

        if length(viable) == 1 || reps >= n_timed
            survivors = viable
            break
        end

        # Round one only: if the heuristic is already within the accept margin of
        # the best candidate, stop refining and let the final paired round
        # decide. Deliberately compared at the SAME resolution -- both scores
        # come out of this round -- so this is a like-for-like read, just a
        # coarse one, and it never decides the verdict itself.
        # EARLY-OUT: is this decision worth more rounds?
        #
        # The heuristic's own plan is almost always already in the ladder -- it
        # differs from a swept candidate only in bookkeeping fields, so it is
        # matched on (mode, allotment, scheduler) and its ROUND-ONE CANDIDATE
        # SCORE is reused. It is deliberately not measured as an extra entry:
        # appending it puts its block last in the round and it then reads ~12 %
        # slow purely from position (0.259 against the same plan's 0.231 as a
        # candidate, measured), which is the same ordering contamination that
        # made the confirm-on-hit unusable.
        #
        # Rank, not ratio: `min` over N noisy two-sample scores is biased low, so
        # a ratio test flatters the candidates and never fires. Top-3 is the
        # record's own resolution limit -- the true best sits there 85 % of the
        # time at 0.0 % median regret -- so a heuristic inside it will not be
        # beaten by the accept margin however the rounds below reshuffle.
        #
        # This skips the REFINEMENT rounds, never the decision: the round-one
        # leader still goes to the final paired round and is still compared
        # against the heuristic at full resolution. On a workload where the
        # ladder is a near-tie that is ~40 calls instead of ~112; where a plan
        # genuinely wins (trx50: satellite_batch 0.108 against flat 0.5-1.4) the
        # heuristic ranks far down, nothing fires, and the full budget is spent.
        if round_index == 1 && heuristic !== nothing
            hmatch = findfirst(
                c -> c.mode === heuristic.mode &&
                     c.allotment == heuristic.allotment &&
                     c.scheduler === heuristic.scheduler,
                viable
            )
            if hmatch !== nothing
                h = scores[viable[hmatch]]
                rank = 1 + count(c -> scores[c] < h, viable)
                if verbose
                    println("  [x$(lpad(reps, 2))] heuristic plan ranks $(rank) of " *
                            "$(length(viable)) ($(round(h / 1e6, digits=3)) ms/call)")
                end
                if rank <= 3
                    verbose && println(
                        "  → heuristic already in the top 3; further rounds cannot " *
                        "change the verdict, skipping to the decisive comparison"
                    )
                    sort!(viable; by = c -> scores[c])
                    survivors = viable[1:1]
                    break
                end
            end
        end

        keep = max(1, cld(length(viable), 2))
        sort!(viable; by = c -> scores[c])
        survivors = viable[1:keep]
        reps = min(n_timed, reps * 2)
    end

    best_plan = survivors[1]
    best_elapsed = scores[best_plan]
    for candidate in survivors
        if scores[candidate] < best_elapsed
            best_elapsed = scores[candidate]
            best_plan = candidate
        end
    end

    # Among survivors that are indistinguishable from the winner, take the
    # NARROWEST rather than the nominal minimum.
    #
    # Cold validation of the cost model established that the true best sits in
    # the top-2 of a ranked candidate set 65% of the time and the top-3 85%,
    # both at 0.0% median regret -- that is, near-ties are unresolvable by
    # measurement AND do not matter, because the candidates involved cost the
    # same. That evidence was collected and never fed back into the chooser,
    # which still resolves a near-tie by raw minimum and so lets noise pick
    # between plans that are equally good on the sweep but not equally good
    # afterwards.
    #
    # Narrower is the right way to break it. A wide plan and a narrow plan that
    # time the same on an idle pre-solve probe are not the same bet during a
    # solve: the wide one has more straggler exposure, holds more of the pool
    # against anything else that wants it, and is the one that inverts on the
    # workloads whose scaling curve turns over. Preferring the narrow one costs
    # nothing when they really are equal and protects the tail when they are
    # not.
    tie_margin = _rhs_calibrate_tie_margin()
    if tie_margin > 0.0 && isfinite(best_elapsed)
        threshold_ns = best_elapsed * (1.0 + tie_margin)
        for candidate in survivors
            isfinite(scores[candidate]) || continue
            scores[candidate] <= threshold_ns || continue
            if candidate.allotment < best_plan.allotment
                best_plan = candidate
                best_elapsed = scores[candidate]
            end
        end
    end

    heuristic_ns = 0.0
    # Measure the heuristic on the same footing and keep it unless the swept
    # winner clears it by more than the sweep's own resolution. Returning
    # `nothing` leaves rhs_plan_override unset, so the runtime heuristic runs --
    # which is strictly better than pinning a copy of it, since it re-derives
    # per call rather than freezing the pre-solve answer.
    if heuristic !== nothing
        # The swept winner and the heuristic are re-measured TOGETHER, in one
        # round, through the same `_measure_round!` the halving used.
        #
        # Previously the winner's score came out of the halving rounds (which
        # interleave) while the heuristic was measured as a contiguous block, so
        # the accept rule compared two numbers produced by different measurement
        # regimes -- worth up to 3.9x on a flat candidate (see
        # `_rhs_calibrate_interleave`). Symmetry is the property that matters,
        # not which regime wins: routing both through `_measure_round!` keeps
        # them comparable whichever way that knob is set.
        #
        # Re-measuring the winner is not redundant. It is the decisive
        # comparison of the whole sweep, and it is the one place worth spending
        # a paired, equal-budget measurement: everything before it only has to
        # rank candidates well enough to pick a finalist.
        if best_plan == heuristic
            # The sweep's winner IS what the runtime heuristic would pick.
            # Pinning a copy of it is strictly worse than pinning nothing, since
            # the heuristic re-derives per call against the live satellite count
            # and outer-split state while a pinned plan freezes the pre-solve
            # answer. No measurement can separate them, so do not spend one.
            verbose && println("  → heuristic retained (swept winner is the heuristic's own plan)")
            # No rival: nothing was compared, so there is no evidence to store
            # and a later confirm has nothing to re-run. Signalled by the zeros.
            return nothing, 0.0, :heuristic, 0.0, nothing
        end
        final_reps = max(2, n_timed)
        final_scores = Dict{Any, Float64}()
        total_calls += _measure_round!(Any[best_plan, heuristic], final_reps, final_scores)
        best_final  = get(final_scores, best_plan, best_elapsed)
        heuristic_ns = get(final_scores, heuristic, Inf)
        margin = _rhs_calibrate_override_margin()
        # An unmeasurable heuristic must not win by default: `Inf` here means the
        # floor could not be evaluated, and retaining it would hand the solve a
        # plan the probe could not even run once.
        if !isfinite(best_final) ||
           (isfinite(heuristic_ns) && best_final > heuristic_ns * (1.0 - margin))
            if verbose
                println("  → heuristic retained ($(round(heuristic_ns / 1e6, digits=3)) ms/call " *
                        "vs best swept $(round(best_final / 1e6, digits=3)); " *
                        "margin $(round(100 * margin, digits=1))% not cleared)")
            end
            # :heuristic, not :aborted. This is a MEASURED verdict -- the sweep
            # ran, the floor compared, and the heuristic won -- so it is worth
            # caching. The :aborted returns above are failures to measure and
            # must stay uncached, or one transient error would permanently
            # suppress calibration for that signature.
            return nothing, heuristic_ns, :heuristic, best_final, best_plan
        end
        # Carry the paired measurement forward: it is the one taken on the same
        # footing as the number it beat, so it is what should be cached.
        best_elapsed = best_final
    end

    if verbose && best_plan !== nothing
        label = best_plan.mode == :satellite_batch ?
            "satellite_batch" :
            "flat($(best_plan.scheduler), allotment=$(best_plan.allotment))"
        println("  → best: $(label)  ($(round(best_elapsed / 1e6, digits=3)) ms/call, " *
                "$(total_calls) timed calls)")
    end

    return best_plan, best_elapsed, :pinned, heuristic_ns, heuristic
end

# ── Public entry point ────────────────────────────────────────────────────────


# ── Ahead-of-time plan-mode compilation ──────────────────────────────────────
#
# Precompilation runs single-threaded (`Threads.nthreads() == 1` in the
# precompile subprocess), so calibration can never compile the plan modes on its
# own however the workload is configured: `_calibrate_rhs_plan_if_needed!`
# returns at the `effective_inner_thread_budget() <= 1` gate, and even past it
# `_rhs_plan_candidates` collapses its flat ladder to `allotment = [1]` because
# the ladder is geometric in the thread budget.
#
# So the modes are driven directly instead. `allotment` and `scheduler` are
# runtime values carried in the plan NamedTuple -- that is the whole reason the
# sweep can vary them without recompiling -- so the values used here do not
# restrict what the compiled code can run later. Only the two plan *modes* and
# the two schedulers need separate compilation.
#
# Coverage is per type-combination, not per plan: `spacecraft_dynamics!`
# specialises on `ODEParams{SimulationConfiguration{P, D, E, T, DM}}`, so this
# hook only helps a configuration whose planet, density model, ephemerides
# model, thermal model and dynamic-effector tuple all match one the precompile
# workload actually ran. Measured: 0.048 s first solve for a covered
# combination against 31.5 s for an otherwise identical one differing in a
# single effector type. Satellite count is not a type parameter.
@inline function _rhs_precompile_plans_requested()::Bool
    return SimulationModel.ParallelPolicy.parse_bool_env(
        "SPACEAGORA_RHS_PRECOMPILE_PLANS", false
    )
end

function _precompile_rhs_plans!(p, u0)::Nothing
    du = zero(u0)
    plans = Any[
        _make_calib_satellite_batch_plan(),
        _make_calib_flat_plan(2, :static),
        _make_calib_flat_plan(2, :dynamic),
    ]
    for plan in plans
        p.shared_buffers.rhs_plan_override[] = plan
        try
            spacecraft_dynamics!(du, u0, p, 0.0)
        catch
            # A plan this configuration cannot run is not an error here. The
            # point is to compile whatever compiles; the sweep applies its own
            # viability rules at run time.
        finally
            p.shared_buffers.rhs_plan_override[] = nothing
        end
    end
    return nothing
end

# ── Solve-cost gate ───────────────────────────────────────────────────────────
#
# Whether calibrating is worth it at all is a function of how long the solve is,
# and the crossover is measurable rather than a matter of taste.
#
# Measured on this machine, five repeat solves in one warm process so only the
# calibration path varies (interact_256sat_1hr, 8 threads):
#
#     calibration off          15.3 ms min
#     cache hit, no sweep      15.0 ms min
#     always sweep             91.6 ms min      -> the sweep costs 76-118 ms
#
# The sweep therefore costs ~0.10 s and exists to avoid ~13 % of regret (the
# recorded swing on this workload is 17 points; the individual mistakes are
# +24 % and +34 %). It pays for itself when 0.10 < 0.13 x T_solve, i.e. once the
# solve is longer than roughly a second.
#
# Above that line the right policy is to sweep EVERY process and never replay a
# cached verdict: a verdict is a reading of the process that formed it (§2.7 of
# the record -- satellite_batch alone moves 1.76x between processes on one
# machine at one signature), and 4 % of a long solve is a cheap price for not
# betting on a stale one. Below it the sweep costs more than the regret it
# prevents -- on the 15 ms test-profile solve it is 500-800 % overhead, which is
# the regime the cached-heuristic path was added for (light_16_harm +30.7 %,
# light_64_aero +14.9 %) -- so the cache is honoured and the sweep runs once.
#
# This replaces a confirm-on-hit probe that re-measured the cached verdict
# before honouring it. That probe was implemented and measured and did not work:
# it caught a deliberately stale entry 0/5 times, because measuring a flat plan
# next to a satellite_batch block reproducibly inflates the flat plan ~3x
# (0.72 against its true 0.23 in the same process). Reproducing the sweep's
# measurement context in miniature turned out to be the hard part; this gate
# needs no such assumption, because above the threshold it simply runs the
# sweep.
@inline function _rhs_calibrate_min_solve_seconds()::Float64
    raw = strip(_engine_env_get("SPACEAGORA_RHS_CALIBRATE_MIN_SOLVE_S", "1.0"))
    v = try
        parse(Float64, raw)
    catch
        throw(ArgumentError("SPACEAGORA_RHS_CALIBRATE_MIN_SOLVE_S must be a float, got '$raw'"))
    end
    return max(0.0, v)
end

# Signature -> when calibration finished, so the solve that follows can be timed.
# Last-writer-wins under concurrent solves of the same signature in one process;
# this feeds a threshold comparison, not a measurement, so a lost sample costs a
# redundant sweep rather than a wrong answer.
const _rhs_calib_solve_start = Dict{String, UInt64}()

# True when this signature has been seen to solve for longer than the threshold,
# i.e. we are in the always-sweep regime. An unmeasured signature answers TRUE:
# the first encounter sweeps, which is also the only conservative choice, since
# nothing is yet known about what the cache would be replaying into.
function _rhs_calib_solve_exceeds_threshold(sig::String)::Bool
    _rhs_calib_load!()
    entry = lock(_rhs_calib_lock) do
        get(_rhs_calib_cache, sig, nothing)
    end
    entry === nothing && return true
    solve_ns = Float64(get(entry, "solve_ns", 0.0))
    (isfinite(solve_ns) && solve_ns > 0.0) || return true
    return solve_ns >= _rhs_calibrate_min_solve_seconds() * 1e9
end

# Called once the solve is done. See the call in simulation/engine/execution.jl.
function _rhs_calib_record_solve_time!()::Nothing
    now = time_ns()
    lock(_rhs_calib_lock) do
        isempty(_rhs_calib_solve_start) && return nothing
        for (sig, started) in collect(_rhs_calib_solve_start)
            delete!(_rhs_calib_solve_start, sig)
            now > started || continue
            entry = get(_rhs_calib_cache, sig, nothing)
            entry === nothing && continue
            sample = Float64(now - started)
            # MINIMUM across samples, not the latest.
            #
            # The first solve in a process is dominated by compilation -- 24.6 s
            # measured against a 15 ms solve -- so a single cold sample would
            # classify every workload as long and the cache would never be used.
            # Later solves in the same process are warm and are the ones that
            # describe the solve itself. The minimum picks those out without
            # needing to know which sample was cold.
            #
            # A process that only ever solves once therefore never contributes a
            # warm sample and stays in the always-sweep regime. That is the right
            # answer for it: it spent tens of seconds compiling, so the sweep's
            # ~0.1 s is a rounding error there anyway.
            prev = Float64(get(entry, "solve_ns", 0.0))
            entry["solve_ns"] = (isfinite(prev) && prev > 0.0) ? min(prev, sample) : sample
        end
        return nothing
    end
    _rhs_calib_save!()
    return nothing
end

function _calibrate_rhs_plan_if_needed!(p, u0, args)
    # Ahead of every gate below: precompilation is single-threaded, so the
    # budget gate would otherwise return before either plan mode is exercised.
    if _rhs_precompile_plans_requested()
        _precompile_rhs_plans!(p, u0)
        return
    end
    _rhs_calibration_mode() == :off && return
    SimulationModel.ParallelPolicy.effective_inner_thread_budget() <= 1 && return

    dynamic_effectors = args.dynamics_model.dynamic_effectors
    isempty(dynamic_effectors) && return
    count(identity, p.is_active) < 2 && return

    sig     = _rhs_calib_signature(
        p, dynamic_effectors, args.environment_model.density_model
    )
    verbose = args.simulation_settings.verbose

    # Note the start unconditionally: whichever way the gate goes, the solve that
    # follows is the sample that decides the regime next time.
    lock(_rhs_calib_lock) do
        _rhs_calib_solve_start[sig] = time_ns()
    end

    if _rhs_calibration_mode() != :force && !_rhs_calib_solve_exceeds_threshold(sig)
        cached = _rhs_calib_lookup(sig)
        if cached === :heuristic
            # Cached retain-the-heuristic verdict: pin nothing, and -- the point --
            # do not sweep. rhs_plan_override stays unset, so _rhs_execution_plan
            # runs the per-call heuristic exactly as it does when a live sweep
            # reaches this answer. Recorded as a plan selection so the outcome is
            # still visible in policy_telemetry_snapshot rather than looking like
            # calibration never ran.
            SimulationModel.ParallelPolicy.record_rhs_plan_selection!(
                :cache, :heuristic, 0, :none
            )
            verbose && println(
                "[SpaceAGORA] RHS calibration: cached verdict → heuristic retained (no sweep)"
            )
            return
        elseif cached !== nothing
            p.shared_buffers.rhs_plan_override[] = cached
            SimulationModel.ParallelPolicy.record_rhs_plan_selection!(
                :cache, cached.mode, cached.allotment, cached.scheduler
            )
            if verbose
                label = cached.mode == :satellite_batch ?
                    "satellite_batch" :
                    "flat($(cached.scheduler), allotment=$(cached.allotment))"
                println("[SpaceAGORA] RHS calibration: loaded cached plan → $(label)")
            end
            return
        end
    end

    # The two trailing values are the sweep's own record of what the winner was
    # measured against. Nothing consumes them any more -- the confirm-on-hit that
    # did was removed, see the solve-cost gate above -- but the sweep still
    # computes them and they are the natural place to hang outcome feedback off
    # when phase 6 lands.
    best_plan, best_elapsed, verdict, _rival_ns, _rival_plan =
        _run_rhs_sweep!(p, u0, dynamic_effectors, verbose, args)

    if best_plan === nothing
        # The sweep ran and grew the buffer even though nothing was pinned, so
        # the heuristic now inherits a partials buffer sized to the widest
        # candidate tried. Release it for the same reason as the pinned path.
        _rhs_calibrate_release_scratch() &&
            _release_oversized_flat_scratch!(p.shared_buffers, 1)
        if verdict === :heuristic && _rhs_calibrate_cache_heuristic()
            SimulationModel.ParallelPolicy.record_rhs_plan_selection!(
                :sweep, :heuristic, 0, :none
            )
            _rhs_calib_store_heuristic!(sig, best_elapsed)
            _rhs_calib_save!()
        end
        return
    end

    p.shared_buffers.rhs_plan_override[] = best_plan
    # The sweep grew the flat partials buffer to its widest candidate; the solve
    # runs at best_plan.allotment. Leaving it oversized costs a strided zeroing
    # instead of a memset on every RHS call for the rest of the run.
    _rhs_release_oversized_scratch(p, best_plan)
    SimulationModel.ParallelPolicy.record_rhs_plan_selection!(
        :sweep, best_plan.mode, best_plan.allotment, best_plan.scheduler
    )
    _rhs_calib_store!(sig, best_plan, best_elapsed)
    _rhs_calib_save!()

    return nothing
end

@inline function _rhs_release_oversized_scratch(p, plan)::Nothing
    _rhs_calibrate_release_scratch() || return nothing
    _release_oversized_flat_scratch!(p.shared_buffers, max(1, plan.allotment))
    return nothing
end

# ── In-run width identification ───────────────────────────────────────────────
#
# The pre-solve sweep answers "which plan is fastest" by stopping the world and
# timing candidates before the solve starts. That has three costs the solve
# itself does not have to pay:
#
#   - It is charged entirely to setup. On light shapes -- small constellations,
#     short effector queues -- the fixed cost of roughly 110 discarded RHS
#     evaluations is large against the solve, which is why the no-regret floor
#     exists and why a cached "the heuristic won" verdict had to be added.
#   - It measures ONE state. The sweep runs at a single epoch, so a workload
#     whose cost varies over the mission -- an aerobraking pass against a coast
#     arc, the atmosphere entering and leaving the picture -- is ranked on
#     whichever regime the probe happened to sit in.
#   - It cannot revisit. A plan pinned before the solve is pinned for all of it.
#
# Identifying in-run removes all three: the observations are RHS evaluations the
# solve was going to perform anyway, they are spread across the mission because
# the solve is, and the trial can be restarted.
#
# What it costs instead is that the observations are noisy -- they carry
# whatever the solver, the callbacks and the machine were doing at the time --
# which is precisely what the interleaved, order-rotating, sign-tested design of
# StreamingPairedTrial is for. It can return "too close to call", and on this
# workload that is usually the right answer.

"""
    RhsWidthTrial

An in-run identification of the RHS execution plan: the candidate plans, the
streaming trial that ranks them, and the widths they correspond to so the
observations can also be fitted to a scalability law.
"""
mutable struct RhsWidthTrial
    trial::SimulationModel.ParallelCost.StreamingPairedTrial
    plans::Vector{Any}
    widths::Vector{Int}
    signature::String
    committed::Bool
end

# Identification verdicts share the calibration store but not its namespace.
#
# Same file, same machine keying, same one-time load and atexit save -- there is
# no reason for a second cache with the same lifecycle. A separate PREFIX,
# though, because the two mechanisms answer the same question by different means
# and a sweep entry must never be read as an identification verdict or the
# reverse: the sweep times candidates in isolation before the solve, the trial
# times them inside it, and where they disagree that disagreement is information
# rather than something to average away. This mirrors the `split|` prefix the
# outer-route selector uses, whose comment records what sharing a bucket cost.
const _IDENTIFY_SIGNATURE_PREFIX = "identify|"

@inline function _rhs_identify_signature(p, dynamic_effectors)::String
    return _IDENTIFY_SIGNATURE_PREFIX * _rhs_calib_signature(
        p, dynamic_effectors, p.args.environment_model.density_model)
end

@inline function _rhs_identify_persist_enabled()::Bool
    return SimulationModel.ParallelPolicy.parse_bool_env(
        "SPACEAGORA_RHS_IDENTIFY_PERSIST", true)
end

# Whether the native-lock Amdahl bound clamps inner width.
#
# OFF BY DEFAULT, and that is a measured result rather than caution.
#
# The premise was sound and the arithmetic is right: a lock caps achievable
# speedup at 1/rho however wide the split, so widths above that ceiling add
# contention and no throughput. Live GRAM density measures rho = 0.89 under a
# pure-RHS probe -- a ceiling of 2, meaning threads cannot help it at all.
#
# It does not survive contact with a solve. On a live-GRAM solve at N=16 the
# probe-derived rho came out below the 0.02 floor, so no ceiling was applied,
# while the probe needed to reach that conclusion cost 13%: paired,
# order-alternating, 11 pairs, cap on lost 1-10 at a median ratio of 1.13,
# p = 0.012. A mechanism that measures a workload in order to decline to
# constrain it has paid for nothing.
#
# Two things would have to change before this is worth enabling. The rho a
# solve exhibits has to be established across configurations -- 0.89 under the
# probe and below 0.02 in a solve of the same density model is not a
# discrepancy that a floor constant resolves. And the measurement has to come
# free, from evaluations the solve performs anyway, rather than from five extra
# warm passes bolted onto setup.
@inline function _rhs_lock_width_cap_enabled()::Bool
    return SimulationModel.ParallelPolicy.parse_bool_env(
        "SPACEAGORA_RHS_LOCK_WIDTH_CAP", false)
end

# Occupancy below which no ceiling is computed. Below it the implied ceiling
# exceeds any real core budget anyway, so "no constraint" is both cheaper and
# more honest than a large number derived from a handful of acquisitions.
@inline function _rhs_lock_cap_floor_rho()::Float64
    raw = strip(_engine_env_get("SPACEAGORA_RHS_LOCK_CAP_FLOOR_RHO", "0.02"))
    v = tryparse(Float64, raw)
    return (v === nothing || v <= 0.0) ? 0.02 : v
end

@inline function _rhs_identify_enabled()::Bool
    return SimulationModel.ParallelPolicy.parse_bool_env("SPACEAGORA_RHS_IDENTIFY", false)
end

@inline function _rhs_identify_rounds()::Int
    n = try
        parse(Int, strip(_engine_env_get("SPACEAGORA_RHS_IDENTIFY_ROUNDS", "9")))
    catch
        9
    end
    return max(2, n)
end

"""
    build_rhs_width_trial(p, dynamic_effectors) -> Union{Nothing, RhsWidthTrial}

Assemble the candidate plans for in-run identification, or `nothing` when there
is nothing to identify.

The candidate set is the sweep's own -- `_rhs_plan_candidates` -- so the two
mechanisms rank the same alternatives and a comparison between them is about
*how* the answer is reached rather than about which answers were available. Arm
one is the width-1 flat plan, the incumbent: a trial that never reaches
significance therefore leaves the solve running serially inside the RHS, which
is the conservative outcome rather than an arbitrary one.
"""
function build_rhs_width_trial(p, dynamic_effectors)::Union{Nothing, RhsWidthTrial}
    budget = SimulationModel.ParallelPolicy.effective_inner_thread_budget()
    budget > 1 || return nothing

    # A verdict already reached on this shape, on this machine, is the whole
    # point of persisting one: the exploration is paid once rather than by every
    # solve. Three outcomes, the same three the sweep's cache has -- a plan to
    # pin, a decision to pin nothing, or a miss that means go and find out.
    if _rhs_identify_persist_enabled()
        cached = _rhs_calib_lookup(_rhs_identify_signature(p, dynamic_effectors))
        if cached === :heuristic
            # Identified once and could not separate the arms. Running the trial
            # again would re-pay the exploration to reach the same answer.
            p.shared_buffers.rhs_plan_override[] = nothing
            return nothing
        elseif cached !== nothing
            p.shared_buffers.rhs_plan_override[] = cached
            return nothing
        end
    end

    candidates = _rhs_plan_candidates(p, dynamic_effectors)
    length(candidates) >= 2 || return nothing

    # `_rhs_plan_candidates` returns pinnable plans, not descriptors, so they
    # are used as they are. Arm one is the narrowest, which makes the trial's
    # ratio read directly as "how much did widening buy".
    # ONE ARM PER (mode, width). The sweep crosses the width ladder with both
    # schedulers, which is right for a sweep -- it times every candidate to
    # convergence and can afford the width. A sign test cannot: this repo has
    # already measured that opening a five-rung width ladder to the outer-route
    # selector left it in a WORSE steady state than three rungs, because ranking
    # many noisy arms from short samples is unreliable in a way that ranking few
    # is not. The scheduler axis is the one to drop, since `static` is the
    # measured default and the sweep remains available to cross both.
    seen = Set{Tuple{Symbol, Int}}()
    plans = Any[]
    widths = Int[]
    for cand in sort(collect(candidates); by = c -> (_rhs_plan_width(c), c.scheduler !== :static))
        key = (cand.mode, _rhs_plan_width(cand))
        key in seen && continue
        push!(seen, key)
        push!(plans, cand)
        push!(widths, _rhs_plan_width(cand))
    end
    length(plans) >= 2 || return nothing

    # A LENGTH GATE, and it is not optional -- without it identification is a
    # measured regression rather than a cost.
    #
    # The trial spends `n_arms * (rounds + 1)` evaluations deliberately running
    # arms it expects to be bad, and only then commits. That pays when the solve
    # has many evaluations left afterwards and does not when it does not.
    # Measured on the 256-satellite fullstack shape at a 600 s mission: a warm
    # solve is ~73 ms at ~0.9 ms per RHS pass, so roughly 80 evaluations, of
    # which a six-arm nine-round trial wants 60. Paired against identification
    # off, 21 pairs in one process, it lost every pair at a median of 1.67x --
    # it spent the solve exploring and committed with nothing left to exploit.
    #
    # So the trial is built only when the solve is long enough to amortise it.
    # The estimate is deliberately crude: the decision is order-of-magnitude
    # ("is this a hundred evaluations or a hundred thousand"), and a precise
    # step count is not available before the adaptive solver has run.
    rounds = _rhs_identify_rounds()
    trial_calls = length(plans) * (rounds + 1)
    estimated_calls = _rhs_estimated_evaluations(p)
    if estimated_calls < _rhs_identify_min_ratio() * trial_calls
        return nothing
    end

    return RhsWidthTrial(
        SimulationModel.ParallelCost.StreamingPairedTrial(
            length(plans); rounds = rounds, warmup_rounds = 1),
        plans, widths, _rhs_identify_signature(p, dynamic_effectors), false)
end

# Rough count of RHS evaluations a solve will perform: accepted steps times
# stages, with rejections ignored. Only the magnitude matters -- see the length
# gate above -- and anything more precise would need the adaptive controller's
# own step history, which does not exist yet when this is called.
function _rhs_estimated_evaluations(p)::Float64
    args = p.args
    mission_s = Float64(args.mission_configuration.mission_time)
    dt = Float64(args.integration_tolerances.dt_max_orbit)
    (isfinite(mission_s) && mission_s > 0.0) || return 0.0
    (isfinite(dt) && dt > 0.0) || (dt = 60.0)
    steps = mission_s / dt
    # Stages per step for the explicit methods this engine defaults to. A stiff
    # solver performs more, so this under-estimates rather than over-estimates,
    # which errs toward not identifying.
    return max(0.0, steps * 7.0)
end

# How many times the trial's own cost the remaining solve must be worth before
# identifying is allowed. Four, not one: breaking even is not a reason to run a
# mechanism, and the arms the trial explores are measurably worse than the one
# it will pick -- on the ladder measured above the mean arm costs roughly twice
# the best, so a ratio of one would still lose.
@inline function _rhs_identify_min_ratio()::Float64
    raw = strip(_engine_env_get("SPACEAGORA_RHS_IDENTIFY_MIN_RATIO", "4.0"))
    v = tryparse(Float64, raw)
    return (v === nothing || v <= 0.0) ? 4.0 : v
end

# The width a plan actually runs at, so the fitted scalability parameters are
# indexed by something physical rather than by arm number.
#
# `satellite_batch` takes its width from Polyester's own pool and honours
# neither `allotment` nor the inner thread budget, so it is reported at the full
# budget rather than at whatever its allotment field happens to say.
@inline function _rhs_plan_width(plan)::Int
    plan.mode === :satellite_batch &&
        return max(1, SimulationModel.ParallelPolicy.effective_inner_thread_budget())
    return max(1, Int(plan.allotment))
end

"""
    rhs_width_trial_step!(du, u, p, t, wt, dispatch!) -> Any

Run one RHS evaluation under the trial: pin the arm's plan, time the call, and
record it.

Timing wraps the dispatch and nothing else -- not the plan lookup, not the
bookkeeping below -- because what the trial is comparing is the cost of the
evaluation under each plan.

When the trial finishes it writes its answer into `rhs_plan_override` and clears
itself from `shared_buffers`, so every subsequent call takes the ordinary path
with no trial branch and no per-call switching. An undecided trial writes
nothing, which leaves the runtime heuristic in charge.
"""
function rhs_width_trial_step!(du, u, p, t::Float64, wt::RhsWidthTrial, dispatch!)
    PC = SimulationModel.ParallelCost
    if !PC.trial_active(wt.trial)
        _rhs_width_trial_commit!(p, wt)
        return dispatch!(du, u, p, t)
    end

    arm = PC.next_arm(wt.trial)
    @inbounds p.shared_buffers.rhs_plan_override[] = wt.plans[arm]
    t0 = time_ns()
    result = dispatch!(du, u, p, t)
    PC.observe!(wt.trial, time_ns() - t0)

    if !PC.trial_active(wt.trial)
        _rhs_width_trial_commit!(p, wt)
    end
    return result
end

# Fit the scalability law to the WIDTH LADDER ONLY.
#
# The arm set mixes two dispatch mechanisms: a flat queue swept across widths,
# and `satellite_batch`, which takes its width from Polyester's own pool. A
# single USL curve through both is not a curve through anything -- the
# functional form describes how one mechanism responds to being widened, and
# fitting it across a mechanism change produced alpha = 0.861 on a ladder whose
# widest flat arm was slower than its narrowest. Arms sharing arm one's mode are
# the ones that vary only in width, so they are the only ones fitted; fewer than
# three of them means there is no curve to fit.
function _rhs_identify_fit(wt::RhsWidthTrial, speedups::Vector{Float64})
    PC = SimulationModel.ParallelCost
    base_mode = wt.plans[1].mode
    idx = [i for i in eachindex(wt.plans) if wt.plans[i].mode === base_mode]
    length(idx) >= 3 || return (0.0, 0.0)
    return PC._fit_usl(wt.widths[idx], speedups[idx])
end

function _rhs_width_trial_commit!(p, wt::RhsWidthTrial)::Nothing
    wt.committed && return nothing
    wt.committed = true
    PC = SimulationModel.ParallelCost
    verdict = PC.trial_verdict(wt.trial)

    if verdict.significant && verdict.arm != 1
        @inbounds plan = wt.plans[verdict.arm]
        p.shared_buffers.rhs_plan_override[] = plan
        if _rhs_identify_persist_enabled()
            _rhs_calib_store!(wt.signature, plan, 0.0)
            _rhs_calib_save!()
        end
    else
        # Undecided, or the incumbent won. Clear the override rather than
        # pinning arm one: the runtime heuristic adapts per call and a trial
        # that could not separate the arms has no evidence for replacing it.
        p.shared_buffers.rhs_plan_override[] = nothing
        # Persist that verdict too. "The arms could not be separated" is a
        # result, and one this shape will reach again at the same cost every
        # run if it is not written down -- which is exactly the defect the
        # sweep's own cached-heuristic entry was added to fix.
        if _rhs_identify_persist_enabled()
            _rhs_calib_store_heuristic!(wt.signature, 0.0)
            _rhs_calib_save!()
        end
    end

    if lowercase(strip(_engine_env_get("SPACEAGORA_RHS_IDENTIFY_TRACE", "0"))) in ("1", "true", "yes", "on")
        speedups = PC.trial_speedups(wt.trial)
        alpha, beta = _rhs_identify_fit(wt, speedups)
        # Arms named by MODE and width, not by index. A trace that printed
        # widths alone was misread once already: `satellite_batch` reports the
        # full budget as its width, so it sits in the same slot as the widest
        # flat plan and the two are indistinguishable in a bare width list. The
        # two routes differ by 5x on this shape, in opposite directions, so
        # which one an index refers to is not a detail.
        labels = [string(pl.mode === :satellite_batch ? "batch" : "flat",
                         "@", w) for (pl, w) in zip(wt.plans, wt.widths)]
        println("[SpaceAGORA] RHS identification: rounds=$(verdict.rounds) " *
                "winner=$(labels[verdict.arm]) significant=$(verdict.significant) " *
                "ratio=$(round(verdict.median_ratio; digits=3)) " *
                "arms=$(labels) speedups=$(round.(speedups; digits=2)) " *
                "alpha=$(round(alpha; digits=4)) beta=$(round(beta; digits=6))")
    end
    p.shared_buffers.rhs_width_trial[] = nothing
    return nothing
end
