function _interp_linear(x::Vector{Float64}, y::Vector{Float64}, xq::Vector{Float64})
    length(x) == length(y) || throw(ArgumentError("x/y length mismatch: $(length(x)) vs $(length(y))"))
    n = length(x)
    n >= 1 || throw(ArgumentError("Interpolation domain is empty."))
    n >= 2 || return fill(y[1], length(xq))

    out = Vector{Float64}(undef, length(xq))
    j = 1
    @inbounds for i in eachindex(xq)
        q = xq[i]
        q <= x[1] && (out[i] = y[1]; continue)
        q >= x[end] && (out[i] = y[end]; continue)
        while j < n - 1 && x[j + 1] < q
            j += 1
        end
        x0, x1 = x[j], x[j + 1]
        y0, y1 = y[j], y[j + 1]
        w = x1 == x0 ? 0.0 : (q - x0) / (x1 - x0)
        out[i] = y0 + w * (y1 - y0)
    end
    return out
end

@inline function _normalized_axis(n::Int)
    n <= 1 && return [0.0]
    return collect(range(0.0, 1.0, length=n))
end

# Placeholder decay-rate fields so every summary row shares one schema; the
# orbit-events apoapsis rows overwrite them via _apo_decay_diagnostic.
const _DECAY_DIAGNOSTIC_EMPTY = (
    drag_decay_ratio_median=NaN,
    drag_decay_ratio_total=NaN,
    drag_decay_n=0
)

# Non-compounding drag-fidelity diagnostic: per-orbit apoapsis decay rates for
# telemetry and simulation, compared at telemetry interval midpoints. Absolute
# apsis errors accumulate any per-pass density bias over the whole campaign
# (and inherit orbit-index misalignment once the period drifts); the ratio of
# decay RATES measures per-pass drag directly. Intervals bracketing a maneuver
# orbit are excluded because the impulse, not drag, dominates them. Ratios are
# only formed where the telemetry rate is resolvable; the median is robust to
# the telemetry altitude quantization, and the total is the aggregate
# sim-to-truth decay over the overlapping span.
function _apo_decay_diagnostic(
    tele_orbit::Vector{Float64},
    tele_alt::Vector{Float64},
    sim_axis::Vector{Float64},
    sim_alt::Vector{Float64},
    maneuver_orbits::Vector{Float64}
)
    (length(tele_orbit) >= 3 && length(sim_axis) >= 3) || return _DECAY_DIAGNOSTIC_EMPTY

    _rates(axis, val) = (
        (val[2:end] .- val[1:(end - 1)]) ./ (axis[2:end] .- axis[1:(end - 1)]),
        (axis[2:end] .+ axis[1:(end - 1)]) ./ 2.0
    )
    d_tele, mid_tele = _rates(tele_orbit, tele_alt)
    d_sim, mid_sim = _rates(sim_axis, sim_alt)

    keep = trues(length(mid_tele))
    @inbounds for k in eachindex(mid_tele)
        if mid_tele[k] < mid_sim[1] || mid_tele[k] > mid_sim[end]
            keep[k] = false
            continue
        end
        for m in maneuver_orbits
            if tele_orbit[k] <= m <= tele_orbit[k + 1]
                keep[k] = false
                break
            end
        end
    end
    any(keep) || return _DECAY_DIAGNOSTIC_EMPTY

    d_sim_interp = _interp_linear(mid_sim, d_sim, mid_tele[keep])
    d_tele_kept = d_tele[keep]
    spans_kept = (tele_orbit[2:end] .- tele_orbit[1:(end - 1)])[keep]

    ratios = Float64[]
    @inbounds for i in eachindex(d_tele_kept)
        abs(d_tele_kept[i]) > 1e-9 && push!(ratios, d_sim_interp[i] / d_tele_kept[i])
    end
    # Rates are km/orbit; the aggregate must weight each interval by its orbit
    # span (equivalently, sum the altitude deltas), or gaps in the telemetry
    # sampling would count the same as single-orbit intervals.
    tele_total = sum(d_tele_kept .* spans_kept)
    ratio_total = abs(tele_total) > 1e-9 ? sum(d_sim_interp .* spans_kept) / tele_total : NaN
    return (
        drag_decay_ratio_median=isempty(ratios) ? NaN : median(ratios),
        drag_decay_ratio_total=ratio_total,
        drag_decay_n=length(ratios)
    )
end

function _compare_orbit_curve(
    scenario::String,
    event::String,
    telemetry_axis::Vector{Float64},
    telemetry_values::Vector{Float64},
    sim_values::Vector{Float64};
    sim_axis::Union{Nothing, Vector{Float64}}=nothing
)
    n_tel = length(telemetry_values)
    n_sim = length(sim_values)
    n_tel == 0 && error("Telemetry series for $scenario/$event is empty.")

    if n_sim == 0
        return (
            scenario=scenario,
            event=event,
            n_telemetry=n_tel,
            n_sim=0,
            coverage=0.0,
            mae_km=Inf,
            rmse_km=Inf,
            max_abs_km=Inf,
            p95_abs_km=Inf,
            bias_km=Inf,
            nmae=Inf,
            nrmse=Inf,
            _DECAY_DIAGNOSTIC_EMPTY...,
            telemetry_axis_start=telemetry_axis[1],
            telemetry_axis_end=telemetry_axis[end]
        ), DataFrame(
            scenario=String[],
            event=String[],
            idx=Int[],
            telemetry_axis=Float64[],
            telemetry_value_km=Float64[],
            sim_interp_value_km=Float64[],
            error_km=Float64[]
        )
    end

    sim_interp = if sim_axis === nothing
        u_tel = _normalized_axis(n_tel)
        u_sim = _normalized_axis(n_sim)
        _interp_linear(u_sim, sim_values, u_tel)
    else
        length(sim_axis) == n_sim || throw(ArgumentError(
            "sim_axis length ($(length(sim_axis))) must match sim_values length ($n_sim) for $scenario/$event"
        ))
        _interp_linear(sim_axis, sim_values, telemetry_axis)
    end

    err = sim_interp .- telemetry_values
    abs_err = abs.(err)
    tel_range = max(maximum(telemetry_values) - minimum(telemetry_values), 1e-9)

    return (
        scenario=scenario,
        event=event,
        n_telemetry=n_tel,
        n_sim=n_sim,
        coverage=min(n_tel, n_sim) / n_tel,
        mae_km=mean(abs_err),
        rmse_km=sqrt(mean(err .^ 2)),
        max_abs_km=maximum(abs_err),
        p95_abs_km=quantile(abs_err, 0.95),
        bias_km=mean(err),
        nmae=mean(abs_err) / tel_range,
        nrmse=sqrt(mean(err .^ 2)) / tel_range,
        _DECAY_DIAGNOSTIC_EMPTY...,
        telemetry_axis_start=telemetry_axis[1],
        telemetry_axis_end=telemetry_axis[end]
    ), DataFrame(
        scenario=fill(scenario, n_tel),
        event=fill(event, n_tel),
        idx=collect(1:n_tel),
        telemetry_axis=telemetry_axis,
        telemetry_value_km=telemetry_values,
        sim_interp_value_km=sim_interp,
        error_km=err
    )
end

function _compare_time_series(
    scenario::String,
    event::String,
    telemetry_time::Vector{Float64},
    telemetry_values::Vector{Float64},
    sim_time::Vector{Float64},
    sim_values::Vector{Float64}
)
    n_tel = length(telemetry_values)
    n_sim = length(sim_values)
    n_tel == 0 && error("Telemetry series for $scenario/$event is empty.")

    if n_sim == 0
        return (
            scenario=scenario,
            event=event,
            n_telemetry=n_tel,
            n_sim=0,
            coverage=0.0,
            mae_km=Inf,
            rmse_km=Inf,
            max_abs_km=Inf,
            p95_abs_km=Inf,
            bias_km=Inf,
            nmae=Inf,
            nrmse=Inf,
            _DECAY_DIAGNOSTIC_EMPTY...,
            telemetry_axis_start=telemetry_time[1],
            telemetry_axis_end=telemetry_time[end]
        ), DataFrame(
            scenario=String[],
            event=String[],
            idx=Int[],
            telemetry_axis=Float64[],
            telemetry_value_km=Float64[],
            sim_interp_value_km=Float64[],
            error_km=Float64[]
        )
    end

    sim_interp = _interp_linear(sim_time, sim_values, telemetry_time)
    err = sim_interp .- telemetry_values
    abs_err = abs.(err)
    tel_range = max(maximum(telemetry_values) - minimum(telemetry_values), 1e-9)

    return (
        scenario=scenario,
        event=event,
        n_telemetry=n_tel,
        n_sim=n_sim,
        coverage=min(n_tel, n_sim) / n_tel,
        mae_km=mean(abs_err),
        rmse_km=sqrt(mean(err .^ 2)),
        max_abs_km=maximum(abs_err),
        p95_abs_km=quantile(abs_err, 0.95),
        bias_km=mean(err),
        nmae=mean(abs_err) / tel_range,
        nrmse=sqrt(mean(err .^ 2)) / tel_range,
        _DECAY_DIAGNOSTIC_EMPTY...,
        telemetry_axis_start=telemetry_time[1],
        telemetry_axis_end=telemetry_time[end]
    ), DataFrame(
        scenario=fill(scenario, n_tel),
        event=fill(event, n_tel),
        idx=collect(1:n_tel),
        telemetry_axis=telemetry_time,
        telemetry_value_km=telemetry_values,
        sim_interp_value_km=sim_interp,
        error_km=err
    )
end

