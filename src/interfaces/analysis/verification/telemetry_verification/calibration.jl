@inline function _calibration_active(cal::CalibrationConfig, profile::Symbol)::Bool
    return cal.enabled && (profile in cal.profiles)
end

"""
    _single_point_calibration(use_calibration, cd_candidates, cr_candidates, eval_profile, profile)

True when calibration is active but the candidate grid has exactly one point
and the eval and final solves share a profile. The eval solve then has the
same configuration as the final solve, so the runner reuses it instead of
solving twice (the bias fit is applied to the final rows either way).
"""
@inline function _single_point_calibration(
    use_calibration::Bool,
    cd_candidates::AbstractVector,
    cr_candidates::AbstractVector,
    eval_profile::Symbol,
    profile::Symbol
)::Bool
    return use_calibration && length(cd_candidates) == 1 && length(cr_candidates) == 1 && eval_profile == profile
end

@inline function _grid_values(min_v::Float64, max_v::Float64, steps::Int)::Vector{Float64}
    if steps <= 1 || isapprox(min_v, max_v; rtol=0.0, atol=1e-12)
        return [min_v]
    end
    return collect(range(min_v, max_v, length=steps))
end

# The bias exists to absorb a constant datum/measurement offset, so it is
# estimated robustly from the earliest compared orbits, where accumulated
# model drift is smallest. A bias at or beyond the configured cap is a
# model-mismatch signal (see the drag_decay_ratio diagnostic), not a datum
# offset — applying it would silently shift accurate apsides — so it is
# reported and NOT applied.
const _BIAS_ESTIMATE_POINTS = 10

function _estimate_event_biases(error_tables::Vector{DataFrame}, bias_abs_max_km::Float64)::Dict{String, Float64}
    out = Dict{String, Float64}()
    for df in error_tables
        nrow(df) == 0 && continue
        event = String(df.event[1])
        errs = _to_float_vector(df.error_km, "error_km:$event")
        k = min(_BIAS_ESTIMATE_POINTS, length(errs))
        bias = -median(errs[1:k])
        if abs(bias) >= bias_abs_max_km
            println(
                "calibration_bias_saturated event=$event bias_km=$(round(bias, digits=3)) " *
                "cap_km=$bias_abs_max_km -- bias NOT applied (model-mismatch signal)"
            )
            bias = 0.0
        end
        out[event] = bias
    end
    return out
end

function _calibration_score(rows::AbstractVector{<:NamedTuple}, objective::String)::Float64
    isempty(rows) && return Inf
    if objective == "mean_nmae"
        return mean([Float64(r.nmae) for r in rows])
    elseif objective == "mean_rmse_km"
        return mean([Float64(r.rmse_km) for r in rows])
    elseif objective == "max_nmae"
        return maximum([Float64(r.nmae) for r in rows])
    end
    throw(ArgumentError("Unsupported calibration objective '$objective'"))
end

function _annotate_calibration_rows(
    rows::AbstractVector{<:NamedTuple},
    cd_scale::Float64,
    cr_value::Float64,
    bias_by_event::Dict{String, Float64},
    score::Float64,
    selected_runtime_s::Float64,
    dt_max_orbit_s::Float64,
    calibration_runtime_s::Float64,
    calibration_used::Bool,
    solver_info
)::Vector{NamedTuple}
    out = NamedTuple[]
    for row in rows
        push!(
            out,
            merge(
                row,
                (
                    calibration_used=calibration_used,
                    calibrated_cd_scale=cd_scale,
                    calibrated_cr=cr_value,
                    calibrated_bias_km=get(bias_by_event, String(row.event), 0.0),
                    calibration_score=score,
                    selected_simulation_runtime_s=selected_runtime_s,
                    dt_max_orbit_s=dt_max_orbit_s,
                    calibration_runtime_s=calibration_runtime_s,
                    solver_mode=solver_info.solver_mode,
                    solver_sequence=solver_info.solver_sequence,
                    solver_fallback_used=solver_info.solver_fallback_used,
                    solver_fallback_count=solver_info.solver_fallback_count,
                    solver_fallback_trigger=solver_info.solver_fallback_trigger,
                    solver_retcode=solver_info.solver_retcode,
                    solver_maxiters=solver_info.solver_maxiters,
                    solver_maxiters_retry_used=solver_info.solver_maxiters_retry_used
                )
            )
        )
    end
    return out
end

