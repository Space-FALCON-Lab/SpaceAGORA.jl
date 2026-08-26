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

function _rhs_calib_signature(p, dynamic_effectors)::String
    budget      = SimulationModel.ParallelPolicy.effective_inner_thread_budget()
    active_sats = count(identity, p.is_active)
    n_eff       = length(dynamic_effectors)
    has_harmonics = n_eff == 1 &&
        dynamic_effectors[1] isa SimulationModel.GravitationalHarmonicsModel
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
        "v3",
        "machine=$(_calib_machine_label())",
        "budget=$(budget)",
        "sats=$(_calib_sat_bucket(active_sats))",
        "effs=$(n_eff)",
        "harm=$(has_harmonics ? "1" : "0")",
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
            _rhs_calib_cache[sig] = Dict{String, Any}(
                "mode"           => String(get(row, "mode", "")),
                "allotment"      => Int(get(row, "allotment", 1)),
                "scheduler"      => String(get(row, "scheduler", "auto")),
                "elapsed_mean_ns"=> Float64(get(row, "elapsed_mean_ns", 0.0)),
            )
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
            push!(rows, Dict{String, Any}(
                "signature"       => sig,
                "mode"            => get(e, "mode", ""),
                "allotment"       => Int(get(e, "allotment", 1)),
                "scheduler"       => get(e, "scheduler", "auto"),
                "elapsed_mean_ns" => Float64(get(e, "elapsed_mean_ns", 0.0)),
            ))
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

function _rhs_calib_lookup(sig::String)::Union{Nothing, NamedTuple}
    _rhs_calib_load!()
    entry = lock(_rhs_calib_lock) do
        get(_rhs_calib_cache, sig, nothing)
    end
    entry === nothing && return nothing
    mode_str  = get(entry, "mode", "")
    allotment = max(1, Int(get(entry, "allotment", 1)))
    scheduler = Symbol(get(entry, "scheduler", "auto"))
    mode_str == "satellite_batch"                  && return _make_calib_satellite_batch_plan()
    mode_str == "flat_constellation_effector_queue" && return _make_calib_flat_plan(allotment, scheduler)
    return nothing
end

function _rhs_calib_store!(sig::String, plan, elapsed_mean_ns::Float64)::Nothing
    # Settle the lazy one-time disk load first: otherwise a later first lookup
    # would load persisted entries over fresher in-process stores.
    _rhs_calib_load!()
    lock(_rhs_calib_lock) do
        _rhs_calib_cache[sig] = Dict{String, Any}(
            "mode"            => String(plan.mode),
            "allotment"       => Int(plan.allotment),
            "scheduler"       => String(plan.scheduler),
            "elapsed_mean_ns" => elapsed_mean_ns,
        )
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

# Measure one candidate as the MINIMUM over `reps` single-call timings.
#
# Minimum, not mean. Interference is strictly additive and one-sided -- a
# preempted call is longer, never shorter -- so the minimum estimates the
# uncontended cost while the mean estimates "this plus whatever else the machine
# was doing". The sweep previously took a mean over its whole timed block, where
# a single 2 ms preemption anywhere inside moved a candidate by 200 microseconds
# per call and could silently change which plan won.
function _rhs_sweep_measure!(p, du, u0, candidate, reps::Int)::Float64
    p.shared_buffers.rhs_plan_override[] = candidate
    try
        best = Inf
        for _ in 1:max(1, reps)
            t0 = time_ns()
            spacecraft_dynamics!(du, u0, p, 0.0)
            sample = Float64(time_ns() - t0)
            sample < best && (best = sample)
        end
        return best
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
    length(candidates) <= 1 && return nothing, 0.0

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
        return nothing, 0.0
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

    while true
        empty!(scores)
        for candidate in survivors
            scores[candidate] = _rhs_sweep_measure!(p, du, u0, candidate, reps)
            total_calls += reps
        end
        viable = [c for c in survivors if isfinite(scores[c])]
        isempty(viable) && return nothing, 0.0

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

    # Measure the heuristic on the same footing and keep it unless the swept
    # winner clears it by more than the sweep's own resolution. Returning
    # `nothing` leaves rhs_plan_override unset, so the runtime heuristic runs --
    # which is strictly better than pinning a copy of it, since it re-derives
    # per call rather than freezing the pre-solve answer.
    if heuristic !== nothing
        heuristic_ns = _rhs_sweep_measure!(p, du, u0, heuristic, max(2, n_timed))
        total_calls += max(2, n_timed)
        margin = _rhs_calibrate_override_margin()
        if !isfinite(best_elapsed) || best_elapsed > heuristic_ns * (1.0 - margin)
            if verbose
                println("  → heuristic retained ($(round(heuristic_ns / 1e6, digits=3)) ms/call " *
                        "vs best swept $(round(best_elapsed / 1e6, digits=3)); " *
                        "margin $(round(100 * margin, digits=1))% not cleared)")
            end
            return nothing, 0.0
        end
    end

    if verbose && best_plan !== nothing
        label = best_plan.mode == :satellite_batch ?
            "satellite_batch" :
            "flat($(best_plan.scheduler), allotment=$(best_plan.allotment))"
        println("  → best: $(label)  ($(round(best_elapsed / 1e6, digits=3)) ms/call, " *
                "$(total_calls) timed calls)")
    end

    return best_plan, best_elapsed
end

# ── Public entry point ────────────────────────────────────────────────────────

function _calibrate_rhs_plan_if_needed!(p, u0, args)
    _rhs_calibration_mode() == :off && return
    SimulationModel.ParallelPolicy.effective_inner_thread_budget() <= 1 && return

    dynamic_effectors = args.dynamics_model.dynamic_effectors
    isempty(dynamic_effectors) && return
    count(identity, p.is_active) < 2 && return

    sig     = _rhs_calib_signature(p, dynamic_effectors)
    verbose = args.simulation_settings.verbose

    if _rhs_calibration_mode() != :force
        cached = _rhs_calib_lookup(sig)
        if cached !== nothing
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

    best_plan, best_elapsed = _run_rhs_sweep!(p, u0, dynamic_effectors, verbose, args)
    best_plan === nothing && return

    p.shared_buffers.rhs_plan_override[] = best_plan
    SimulationModel.ParallelPolicy.record_rhs_plan_selection!(
        :sweep, best_plan.mode, best_plan.allotment, best_plan.scheduler
    )
    _rhs_calib_store!(sig, best_plan, best_elapsed)
    _rhs_calib_save!()

    return nothing
end
