const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

function _read(rel)
    path = joinpath(REPO_ROOT, rel)
    isfile(path) || error("Missing required file: $rel")
    return read(path, String)
end

calib_src   = _read(joinpath("src", "simulation", "engine", "rhs_calibration.jl"))
setup_src   = _read(joinpath("src", "simulation", "engine", "setup.jl"))
engine_src  = _read(joinpath("src", "simulation", "engine", "simulation_engine.jl"))
exec_src    = _read(joinpath("src", "simulation", "engine", "execution.jl"))
types_src   = _read(joinpath("src", "core", "types", "runtime_types.jl"))

# ── §1: All required functions present in rhs_calibration.jl ─────────────────

for fn in (
    "function _rhs_calibration_mode()",
    "function _rhs_calibrate_n_warmup()",
    "function _rhs_calibrate_n_timed()",
    "function _calib_machine_label()",
    "function _calib_sat_bucket(n::Int)",
    "function _rhs_calib_signature(",
    "function _rhs_calib_path()",
    "function _rhs_calib_load!()",
    "function _rhs_calib_save!()",
    "function _rhs_calib_lookup(sig::String)",
    "function _rhs_calib_store!(sig::String",
    "function _make_calib_satellite_batch_plan()",
    # Scheduler is a swept axis, not a constant: the flat plan constructor takes
    # it as a second argument and the sweep crosses it with the allotment ladder.
    "function _make_calib_flat_plan(allotment::Int, scheduler::Symbol",
    "function _rhs_calibrate_schedulers()",
    "function _rhs_plan_candidates(",
    "function _run_rhs_sweep!(",
    "function _calibrate_rhs_plan_if_needed!(",
)
    occursin(fn, calib_src) ||
        error("Required calibration function '$(fn)' is missing from rhs_calibration.jl")
end

# ── §2: Environment variable keys wired up correctly ─────────────────────────

for key in (
    "SPACEAGORA_RHS_CALIBRATE",
    "SPACEAGORA_RHS_CALIBRATION_PATH",
    "SPACEAGORA_RHS_CALIBRATE_N_WARMUP",
    "SPACEAGORA_RHS_CALIBRATE_N_TIMED",
    "SPACEAGORA_RHS_CALIBRATE_SCHEDULERS",
    "SPACEAGORA_PERF_MACHINE_LABEL",
)
    occursin(key, calib_src) ||
        error("Env var '$(key)' is not referenced in rhs_calibration.jl")
end

# ── §3: TOML schema keys present in persistence helpers ───────────────────────

for key in ("\"calibrations\"", "\"signature\"", "\"mode\"", "\"allotment\"", "\"scheduler\"", "\"elapsed_mean_ns\"", "\"schema_version\"")
    occursin(key, calib_src) ||
        error("TOML schema key $(key) is missing from rhs_calibration.jl persistence helpers")
end

# ── §4: Plan constructors return the two expected modes ───────────────────────

occursin(":satellite_batch", calib_src) ||
    error("satellite_batch plan constructor is missing from rhs_calibration.jl")
occursin(":flat_constellation_effector_queue", calib_src) ||
    error("flat_constellation_effector_queue plan constructor is missing from rhs_calibration.jl")

# Both plan constructors must wire a serial effector_decision (no inner threading during calibrated sweep).
n_satellite_batch_plans = length(collect(eachmatch(r"mode\s*=\s*:satellite_batch", calib_src)))
n_flat_plans            = length(collect(eachmatch(r"mode\s*=\s*:flat_constellation_effector_queue", calib_src)))
n_serial_decisions      = length(collect(eachmatch(r"_CALIB_SERIAL_EFFECTOR_DECISION", calib_src)))
n_serial_decisions >= n_satellite_batch_plans + n_flat_plans ||
    error("Fewer _CALIB_SERIAL_EFFECTOR_DECISION references ($n_serial_decisions) than plan modes " *
          "($(n_satellite_batch_plans) satellite_batch + $(n_flat_plans) flat) in rhs_calibration.jl")

# ── §5: Guard conditions present — calibration skips when budget <= 1 or empty effectors ──

occursin("effective_inner_thread_budget() <= 1 && return", calib_src) ||
    error("Single-thread short-circuit guard is missing from _calibrate_rhs_plan_if_needed! in rhs_calibration.jl")
occursin("isempty(dynamic_effectors) && return", calib_src) ||
    error("Empty effectors guard is missing from _calibrate_rhs_plan_if_needed! in rhs_calibration.jl")
occursin("count(identity, p.is_active) < 2 && return", calib_src) ||
    error("Single-satellite guard is missing from _calibrate_rhs_plan_if_needed! in rhs_calibration.jl")

# ── §6: rhs_plan_override field present in SharedBuffers (runtime_types.jl) ──

occursin("rhs_plan_override", types_src) ||
    error("rhs_plan_override field is missing from SharedBuffers in runtime_types.jl")
# Concretely typed Ref: a Ref{Any} here made the plan infer as Any in the RHS,
# boxing every plan access and re-boxing ODEParams per dynamics call (see
# RhsExecutionPlan in runtime_types.jl).
occursin("rhs_plan_override::Base.RefValue{Union{Nothing, RhsExecutionPlan}}", types_src) ||
    error("rhs_plan_override is not typed as Base.RefValue{Union{Nothing, RhsExecutionPlan}} in runtime_types.jl")

# ── §7: Override check present at top of _rhs_execution_plan (setup.jl) ──────

occursin("rhs_plan_override", setup_src) ||
    error("rhs_plan_override override check is missing from setup.jl")
occursin("p.shared_buffers.rhs_plan_override[]", setup_src) ||
    error("p.shared_buffers.rhs_plan_override[] dereference is missing from setup.jl")

# ── §8: rhs_calibration.jl included in simulation_engine.jl ─────────────────

occursin("rhs_calibration.jl", engine_src) ||
    error("rhs_calibration.jl is not included in simulation_engine.jl")

# Must be included after dynamics_rhs.jl (calibration calls spacecraft_dynamics! through the RHS).
dynamics_pos = findfirst("dynamics_rhs.jl", engine_src)
calib_pos    = findfirst("rhs_calibration.jl", engine_src)
dynamics_pos !== nothing && calib_pos !== nothing ||
    error("Could not locate dynamics_rhs.jl or rhs_calibration.jl in simulation_engine.jl include list")
first(dynamics_pos) < first(calib_pos) ||
    error("rhs_calibration.jl is included before dynamics_rhs.jl in simulation_engine.jl — it must come after")

# ── §9: _calibrate_rhs_plan_if_needed! called from execution.jl ──────────────

occursin("_calibrate_rhs_plan_if_needed!(p, u_start, args)", exec_src) ||
    error("_calibrate_rhs_plan_if_needed! is not called in execution.jl")

# The calibration call must appear after jac_prototype construction and before the solve loop.
jac_pos   = findfirst("jac_prototype", exec_src)
calib_call_pos = findfirst("_calibrate_rhs_plan_if_needed!", exec_src)
jac_pos !== nothing && calib_call_pos !== nothing ||
    error("Could not locate jac_prototype or _calibrate_rhs_plan_if_needed! in execution.jl")
first(jac_pos) < first(calib_call_pos) ||
    error("_calibrate_rhs_plan_if_needed! is called before jac_prototype construction in execution.jl")

println("ci_rhs_calibration_gate_ok")
