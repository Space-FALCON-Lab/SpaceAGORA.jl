const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

function _read(rel)
    path = joinpath(REPO_ROOT, rel)
    isfile(path) || error("Missing required file: $rel")
    return read(path, String)
end

setup_src    = _read(joinpath("src", "simulation", "engine", "setup.jl"))
dynamics_src = _read(joinpath("src", "simulation", "engine", "dynamics_rhs.jl"))
assembly_src = _read(joinpath("src", "simulation", "callbacks", "density_callbacks", "assembly.jl"))
types_src    = _read(joinpath("src", "parallel", "policy", "types.jl"))
obs_src      = _read(joinpath("src", "parallel", "policy", "observation_tracking.jl"))
telemetry_src = _read(joinpath("src", "parallel", "policy", "policy_telemetry.jl"))
pp_src       = _read(joinpath("src", "parallel", "policy", "parallel_policy.jl"))

# ── §1/§2: Satellite-batch route contract ───────────────────────────────────
# _satellite_batch_saturates_pool must exist and be used in auto routing.
occursin("function _satellite_batch_saturates_pool(active_sats::Int, budget::Int)::Bool", setup_src) ||
    error("_satellite_batch_saturates_pool helper is missing from setup.jl")
occursin("_satellite_batch_saturates_pool(active_sats, budget)", setup_src) ||
    error("_satellite_batch_saturates_pool is not called in the auto routing section of setup.jl")

# _with_serial_effector_decision must exist and be used.
occursin("function _with_serial_effector_decision(effector_decision)", setup_src) ||
    error("_with_serial_effector_decision helper is missing from setup.jl")
occursin("_with_serial_effector_decision(effector_decision)", setup_src) ||
    error("_with_serial_effector_decision is not called in setup.jl")

# No old-style inline serial effector decision tuples (they were replaced by the helper).
!occursin("effector_decision=(use_threads=false, allotment=1, mode=", setup_src) ||
    error("Old-style inline serial effector_decision tuple found in setup.jl — use _with_serial_effector_decision instead")

# Every satellite_batch return must use the serial coercion helper (count parity).
n_satellite_batch = length(collect(eachmatch(r"mode=:satellite_batch", setup_src)))
n_serial_coercion = length(collect(eachmatch(r"_with_serial_effector_decision\(effector_decision\)", setup_src)))
n_serial_coercion >= n_satellite_batch ||
    error("Fewer _with_serial_effector_decision calls ($n_serial_coercion) than mode=:satellite_batch returns ($n_satellite_batch) in setup.jl — a satellite_batch route may be missing the serial coercion")

# ── §2: dominant_axis field ──────────────────────────────────────────────────
occursin("dominant_axis=:serial", setup_src) ||
    error("dominant_axis=:serial is missing from _rhs_execution_plan in setup.jl")
occursin("dominant_axis=:satellite", setup_src) ||
    error("dominant_axis=:satellite is missing from _rhs_execution_plan in setup.jl")
occursin("dominant_axis=:flat_effector", setup_src) ||
    error("dominant_axis=:flat_effector is missing from _rhs_execution_plan in setup.jl")
occursin("dominant_axis=:per_satellite_inner_effector", setup_src) ||
    error("dominant_axis=:per_satellite_inner_effector is missing from _rhs_execution_plan in setup.jl")

# ── §5: Per-worker work threshold gate ──────────────────────────────────────
occursin("function _rhs_flat_work_per_worker_ns_threshold()", setup_src) ||
    error("_rhs_flat_work_per_worker_ns_threshold is missing from setup.jl")
occursin("SPACEAGORA_RHS_FLAT_WORK_PER_WORKER_NS_THRESHOLD", setup_src) ||
    error("SPACEAGORA_RHS_FLAT_WORK_PER_WORKER_NS_THRESHOLD env var is not wired up in setup.jl")
# The thresholds are read from the run-scoped env snapshot (RhsPlanEnvConfig)
# rather than re-parsed from ENV per RHS call.
occursin("estimated_work_per_worker >= env.flat_work_per_worker_ns_threshold", setup_src) ||
    error("Per-worker work threshold gate is not used in flat route selection in setup.jl")
# The per-worker gate must be part of the compound condition that selects the flat route,
# not a standalone guard — verify both total and per-worker gates appear together.
occursin("estimated_total_work_ns   >= env.flat_work_ns_threshold &&\n        estimated_work_per_worker >= env.flat_work_per_worker_ns_threshold", setup_src) ||
    error("Total and per-worker work threshold gates are not paired correctly in setup.jl")

# ── §3: Staged density callback split ───────────────────────────────────────
occursin("function _requires_density_for_rhs(", assembly_src) ||
    error("_requires_density_for_rhs is missing from assembly.jl")
occursin("function _requires_staged_density_callback(", assembly_src) ||
    error("_requires_staged_density_callback is missing from assembly.jl")
# _requires_density_callback must still exist (backward-compat alias for tests).
occursin("function _requires_density_callback(", assembly_src) ||
    error("_requires_density_callback backward-compat alias is missing from assembly.jl")
# get_callbacks must use the staged variant, not the raw density check.
occursin("_requires_staged_density_callback(effectors, args)", assembly_src) ||
    error("get_callbacks in assembly.jl is not gating density callback installation with _requires_staged_density_callback")
# Staged variant must check at least one non-RHS consumer.
occursin("ParallelPolicy.parse_bool_env(\"SPACEAGORA_FORCE_DENSITY_CALLBACK\"", assembly_src) ||
    error("_requires_staged_density_callback does not check SPACEAGORA_FORCE_DENSITY_CALLBACK env var in assembly.jl")
occursin("_requires_thermal_callback(effectors, args)", assembly_src) ||
    error("_requires_staged_density_callback does not check _requires_thermal_callback in assembly.jl")
occursin("_requires_drag_state_callback(effectors, args)", assembly_src) ||
    error("_requires_staged_density_callback does not check _requires_drag_state_callback in assembly.jl")

# ── §6: Batchable effector trait ────────────────────────────────────────────
occursin("@inline _batchable_effector(::Any)::Bool = false", dynamics_src) ||
    error("Fallback _batchable_effector(::Any) = false trait is missing from dynamics_rhs.jl")
occursin("@inline _batchable_effector(::SimulationModel.NBodyGravityModel)::Bool = true", dynamics_src) ||
    error("NBodyGravityModel is not marked as batchable in dynamics_rhs.jl")
occursin("@inline _batchable_effector(::SimulationModel.SolarRadiationPressureModel)::Bool = true", dynamics_src) ||
    error("SolarRadiationPressureModel is not marked as batchable in dynamics_rhs.jl")

# GRAM/aero must NOT be marked batchable.
!occursin("_batchable_effector(::SimulationModel.AerodynamicCoefficientfM)::Bool = true", dynamics_src) ||
    error("AerodynamicCoefficientfM is incorrectly marked as batchable — GRAM density is per-satellite")

# Batchable pre-pass: only runs on the unpartitioned path.
occursin("partition === nothing && _has_any_batchable_effector(dynamic_effectors)", dynamics_src) ||
    error("Batchable pre-pass is not gated by 'partition === nothing' in dynamics_rhs.jl")

# Early exit when all effectors are batchable.
occursin("_count_non_batchable_effectors(dynamic_effectors) == 0 && return nothing", dynamics_src) ||
    error("Early exit for all-batchable effector set is missing from dynamics_rhs.jl")

# Flat queue workers must skip batchable AND harmonics effectors to prevent double-counting.
occursin("partition === nothing && (_batchable_effector(effector) || _harmonics_prepass_effector(effector)) && return nothing", dynamics_src) ||
    error("Non-packet flat queue worker is missing combined batchable+harmonics skip guard in dynamics_rhs.jl")
occursin("partition === nothing && (_batchable_effector(effector) || _harmonics_prepass_effector(effector)) && continue", dynamics_src) ||
    error("Packet flat queue worker is missing combined batchable+harmonics skip guard in dynamics_rhs.jl")

# needs_state_sample must be true whenever batchable effectors are present so that
# pos_buffers/mass_buffers are prefilled before the batchable pre-pass runs.
occursin("_has_any_batchable_effector(dynamic_effectors)", dynamics_src) ||
    error("needs_state_sample guard does not include _has_any_batchable_effector check in dynamics_rhs.jl")

# Batch kernels for NBody and SRP must exist.
occursin("function _accumulate_nbody_flat_batch!", dynamics_src) ||
    error("_accumulate_nbody_flat_batch! kernel is missing from dynamics_rhs.jl")
occursin("function _accumulate_srp_flat_batch!", dynamics_src) ||
    error("_accumulate_srp_flat_batch! kernel is missing from dynamics_rhs.jl")
occursin("function _accumulate_batchable_effector_flat!", dynamics_src) ||
    error("_accumulate_batchable_effector_flat! dispatch function is missing from dynamics_rhs.jl")

# ── §6b: Harmonics SIMD pre-pass ────────────────────────────────────────────
# The harmonics pre-pass runs the SIMD batch kernel in the multi-effector flat
# queue path, writing directly into totals without re-zeroing the matrix.
occursin("@inline _harmonics_prepass_effector(::Any)::Bool = false", dynamics_src) ||
    error("Fallback _harmonics_prepass_effector(::Any) = false trait is missing from dynamics_rhs.jl")
occursin("@inline _harmonics_prepass_effector(::SimulationModel.GravitationalHarmonicsModel)::Bool = true", dynamics_src) ||
    error("GravitationalHarmonicsModel is not marked as harmonics-prepass in dynamics_rhs.jl")
occursin("function _has_any_harmonics_effector(", dynamics_src) ||
    error("_has_any_harmonics_effector helper is missing from dynamics_rhs.jl")
occursin("function _count_flat_queue_only_effectors(", dynamics_src) ||
    error("_count_flat_queue_only_effectors helper is missing from dynamics_rhs.jl")

# Harmonics pre-pass must be gated by partition === nothing.
occursin("partition === nothing && _has_any_harmonics_effector(dynamic_effectors)", dynamics_src) ||
    error("Harmonics pre-pass is not gated by 'partition === nothing' in dynamics_rhs.jl")

# Early exit when no flat queue work remains after both pre-passes.
occursin("_count_flat_queue_only_effectors(dynamic_effectors) == 0 && return nothing", dynamics_src) ||
    error("Early exit for all-pre-passed effector set is missing from dynamics_rhs.jl")

# _accumulate_harmonics_flat_batch! must accept init_scratch keyword so the pre-pass
# can skip re-zeroing totals (which would overwrite batchable pre-pass contributions).
occursin("init_scratch::Bool=true", dynamics_src) ||
    error("_accumulate_harmonics_flat_batch! is missing init_scratch::Bool=true keyword in dynamics_rhs.jl")
occursin("_accumulate_harmonics_flat_batch!(sc_state, p, t, effector, plan; init_scratch=false)", dynamics_src) ||
    error("Harmonics pre-pass does not call _accumulate_harmonics_flat_batch! with init_scratch=false in dynamics_rhs.jl")

# _prepare_rhs_flat_work_items! must exclude pre-pass effectors from the work list.
occursin("_batchable_effector(effector) || _harmonics_prepass_effector(effector)) && continue", dynamics_src) ||
    error("_prepare_rhs_flat_work_items! does not exclude pre-pass effectors in dynamics_rhs.jl")

# ── §8: PolicyTelemetry proposed / dispatched / discarded counters ──────────
occursin("policy_threading_proposed_total::Int64", types_src) ||
    error("policy_threading_proposed_total is missing from PolicyTelemetry in types.jl")
occursin("policy_threading_dispatched_total::Int64", types_src) ||
    error("policy_threading_dispatched_total is missing from PolicyTelemetry in types.jl")
occursin("policy_discarded_by_route_total::Int64", types_src) ||
    error("policy_discarded_by_route_total is missing from PolicyTelemetry in types.jl")

occursin("policy_threading_proposed_total += use_threads ? 1 : 0", telemetry_src) ||
    error("policy_threading_proposed_total is not incremented in policy_telemetry.jl")
occursin("policy_threading_dispatched_total += 1", obs_src) ||
    error("policy_threading_dispatched_total is not incremented in record_policy_observation! in observation_tracking.jl")
occursin("function record_route_discard!()", obs_src) ||
    error("record_route_discard! is missing from observation_tracking.jl")
occursin("policy_discarded_by_route_total += 1", obs_src) ||
    error("policy_discarded_by_route_total is not incremented in record_route_discard! in observation_tracking.jl")

# record_route_discard! must be exported from the ParallelPolicy module.
occursin("record_route_discard!", pp_src) ||
    error("record_route_discard! is not referenced in parallel_policy.jl (check export list)")

println("ci_rhs_parallel_route_gate_ok")
