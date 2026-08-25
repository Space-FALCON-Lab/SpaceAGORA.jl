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
        # v2: entries now also pin the inner scheduler, which v1 left to ENV.
        # A v1 entry replayed under v2 would silently reinstate whatever
        # scheduler the *current* profile declares, which is exactly the
        # coupling this change removes -- so the version bump invalidates them.
        "v2",
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

function _run_rhs_sweep!(p, u0, dynamic_effectors, verbose::Bool)
    n_warmup   = _rhs_calibrate_n_warmup()
    n_timed    = _rhs_calibrate_n_timed()
    candidates = _rhs_plan_candidates(p, dynamic_effectors)

    # Nothing to compare: only the baseline candidate exists.
    length(candidates) <= 1 && return nothing, 0.0

    if verbose
        println(
            "[SpaceAGORA] RHS calibration: sweeping $(length(candidates)) candidates " *
            "($(n_warmup) warmup + $(n_timed) timed calls each)"
        )
    end

    du         = zero(u0)
    best_plan    = nothing
    best_elapsed = Inf

    for candidate in candidates
        p.shared_buffers.rhs_plan_override[] = candidate
        elapsed_mean = Inf
        try
            for _ in 1:n_warmup
                spacecraft_dynamics!(du, u0, p, 0.0)
            end
            t0 = time_ns()
            for _ in 1:n_timed
                spacecraft_dynamics!(du, u0, p, 0.0)
            end
            elapsed_mean = Float64(time_ns() - t0) / n_timed
        catch e
            @warn "RHS calibration: candidate skipped due to error" mode=candidate.mode allotment=candidate.allotment exception=e
        finally
            p.shared_buffers.rhs_plan_override[] = nothing
        end

        elapsed_mean < Inf || continue

        if elapsed_mean < best_elapsed
            best_elapsed = elapsed_mean
            best_plan    = candidate
        end

        if verbose
            label = candidate.mode == :satellite_batch ?
                "satellite_batch                 " :
                "flat($(rpad(candidate.scheduler, 7)) allotment=$(lpad(candidate.allotment, 3)))"
            println("  $(label) → $(round(elapsed_mean / 1e6, digits=3)) ms/call")
        end
    end

    if verbose && best_plan !== nothing
        label = best_plan.mode == :satellite_batch ?
            "satellite_batch" :
            "flat($(best_plan.scheduler), allotment=$(best_plan.allotment))"
        println("  → best: $(label)  ($(round(best_elapsed / 1e6, digits=3)) ms/call)")
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

    best_plan, best_elapsed = _run_rhs_sweep!(p, u0, dynamic_effectors, verbose)
    best_plan === nothing && return

    p.shared_buffers.rhs_plan_override[] = best_plan
    SimulationModel.ParallelPolicy.record_rhs_plan_selection!(
        :sweep, best_plan.mode, best_plan.allotment, best_plan.scheduler
    )
    _rhs_calib_store!(sig, best_plan, best_elapsed)
    _rhs_calib_save!()

    return nothing
end
