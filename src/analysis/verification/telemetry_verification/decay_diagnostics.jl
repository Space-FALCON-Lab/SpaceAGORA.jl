# Decay/energy diagnostics for telemetry-verification series.
#
# The one rule of this module (see issue #52): a secular orbit-energy slope
# estimated from any finite arc is NOT physical on its own. The J2 secular
# apsidal precession slowly modulates the amplitudes and phases of the
# osculating-SMA short-period harmonics, and ANY fixed-amplitude regression
# leaks that modulation into its linear term (+11–12 m/day on a drag-free
# 48-hour LEO arc — measured, config/grid/tolerance-independent, present for
# pure J2). The leakage is common-mode between a simulation and flight data
# of the same arc, so it cancels exactly under zero-reference subtraction:
# always difference against a drag-free reference propagation scored on the
# same comparison grid ([`zero_referenced_decay`](@ref)). The corrected decay
# is estimator-invariant (validated at the 1% level between fixed- and
# drifting-amplitude estimators on the CYGNSS 48-hour campaign arc).

"""
    visviva_sma(r_m, v_mps, mu) -> a_m

Osculating semi-major axis from position magnitude `r_m` [m] and velocity
magnitude `v_mps` [m/s] via vis-viva. Broadcasts over vectors.
"""
@inline visviva_sma(r_m::Real, v_mps::Real, mu::Real) = 1.0 / (2.0 / r_m - v_mps^2 / mu)
visviva_sma(r_m::AbstractVector{<:Real}, v_mps::AbstractVector{<:Real}, mu::Real) =
    visviva_sma.(r_m, v_mps, mu)

"""
    secular_sma_slope(t_s, sma_m; period_s, n_harmonics=3, drifting_amplitudes=false)
        -> slope_m_per_day

Secular SMA slope by harmonic regression: a linear trend fit jointly with
explicit 1..`n_harmonics`/rev sin/cos terms at the orbital period `period_s`.
This is robust to sampling gaps (unlike rolling-mean pre-smoothers, whose
sample-counted windows stretch past one period at gaps and alias ±50 m/day of
phase-dependent leakage into the slope). With `drifting_amplitudes=true` the
fit also carries t·sin/t·cos columns (first-order secular amplitude drift).

The returned slope still contains the J2-precession window leakage described
in the module header — difference it against a drag-free reference via
[`zero_referenced_decay`](@ref) before physical interpretation.
"""
function secular_sma_slope(
    t_s::AbstractVector{<:Real},
    sma_m::AbstractVector{<:Real};
    period_s::Real,
    n_harmonics::Int=3,
    drifting_amplitudes::Bool=false
)::Float64
    n = length(t_s)
    n == length(sma_m) || throw(ArgumentError(
        "secular_sma_slope: time and SMA vectors must have equal length, got $n vs $(length(sma_m))."
    ))
    min_cols = 2 + 2 * n_harmonics * (drifting_amplitudes ? 2 : 1)
    n > min_cols || throw(ArgumentError(
        "secular_sma_slope: need more samples ($n) than fit columns ($min_cols)."
    ))
    period_s > 0.0 || throw(ArgumentError("secular_sma_slope: period_s must be > 0."))
    n_harmonics >= 1 || throw(ArgumentError("secular_sma_slope: n_harmonics must be >= 1."))
    t = Float64.(t_s)
    tc = t .- (sum(t) / n)
    cols = Vector{Vector{Float64}}()
    push!(cols, tc)
    push!(cols, ones(n))
    for k in 1:n_harmonics
        w = 2.0 * pi * k / Float64(period_s)
        s = sin.(w .* tc)
        c = cos.(w .* tc)
        push!(cols, s)
        push!(cols, c)
        if drifting_amplitudes
            push!(cols, tc .* s)
            push!(cols, tc .* c)
        end
    end
    A = reduce(hcat, cols)
    coef = A \ Float64.(sma_m)
    return coef[1] * 86400.0
end

"""
    zero_referenced_decay(t_s, sma_m, t_ref_s, sma_ref_m; period_s, kwargs...)
        -> (decay_m_per_day, raw_m_per_day, reference_m_per_day)

The physical secular decay of a measured SMA series: the identical estimator
is applied to the measurement and to a DRAG-FREE reference propagation of the
same arc (same gravity field, same window, same gap structure), and the
reference slope — pure estimator/window leakage, since a conservative field
cannot change SMA secularly — is subtracted. Keyword arguments are forwarded
to [`secular_sma_slope`](@ref) for both series.
"""
function zero_referenced_decay(
    t_s::AbstractVector{<:Real},
    sma_m::AbstractVector{<:Real},
    t_ref_s::AbstractVector{<:Real},
    sma_ref_m::AbstractVector{<:Real};
    period_s::Real,
    kwargs...
)
    raw = secular_sma_slope(t_s, sma_m; period_s=period_s, kwargs...)
    ref = secular_sma_slope(t_ref_s, sma_ref_m; period_s=period_s, kwargs...)
    return (decay_m_per_day=raw - ref, raw_m_per_day=raw, reference_m_per_day=ref)
end

"""
    flight_density_table(t_s, sma_m, t_ref_s, sma_ref_m;
                         mu, cd_area_m2, mass_kg, period_s,
                         window_s=21600.0, step_s=3600.0, kwargs...) -> DataFrame

Flight-inferred along-track density history from windowed, zero-referenced
SMA decay rates. Each window's corrected decay is converted through the
near-circular drag relation `da/dt = -rho * (Cd*A/m) * sqrt(mu * a)`, so the
returned density absorbs the supplied effective `cd_area_m2` — scale errors
in Cd*A rescale rho uniformly without changing its shape. The result is a
`DataFrame(time_s, rho_kgm3)` in exactly the `tabulated_time` scenario-source
format (`atmosphere_truth.atmosphere_model = "tabulated_time"`), i.e. ready
to fly a digital-twin replay. `time_s` are window centers in the same
timebase as `t_s`; windows with fewer samples than fit columns are skipped.
Keyword arguments are forwarded to [`secular_sma_slope`](@ref).
"""
function flight_density_table(
    t_s::AbstractVector{<:Real},
    sma_m::AbstractVector{<:Real},
    t_ref_s::AbstractVector{<:Real},
    sma_ref_m::AbstractVector{<:Real};
    mu::Real,
    cd_area_m2::Real,
    mass_kg::Real,
    period_s::Real,
    window_s::Real=21600.0,
    step_s::Real=3600.0,
    n_harmonics::Int=3,
    kwargs...
)::DataFrame
    cd_area_m2 > 0.0 || throw(ArgumentError("flight_density_table: cd_area_m2 must be > 0."))
    mass_kg > 0.0 || throw(ArgumentError("flight_density_table: mass_kg must be > 0."))
    window_s > 0.0 && step_s > 0.0 || throw(ArgumentError(
        "flight_density_table: window_s and step_s must be > 0."
    ))
    t = Float64.(t_s)
    tr = Float64.(t_ref_s)
    min_pts = 2 + 2 * n_harmonics + 1
    centers = Float64[]
    rhos = Float64[]
    c = minimum(t) + window_s / 2.0
    t_max = maximum(t)
    while c <= t_max - window_s / 2.0 + step_s / 2.0
        lo = c - window_s / 2.0
        hi = c + window_s / 2.0
        m = (t .>= lo) .& (t .< hi)
        mr = (tr .>= lo) .& (tr .< hi)
        if count(m) > min_pts && count(mr) > min_pts
            adot = zero_referenced_decay(
                t[m], Float64.(sma_m)[m], tr[mr], Float64.(sma_ref_m)[mr];
                period_s=period_s, n_harmonics=n_harmonics, kwargs...
            ).decay_m_per_day / 86400.0
            a_mean = sum(Float64.(sma_m)[m]) / count(m)
            rho = -adot * Float64(mass_kg) / (Float64(cd_area_m2) * sqrt(Float64(mu) * a_mean))
            if isfinite(rho) && rho > 0.0
                push!(centers, c)
                push!(rhos, rho)
            end
        end
        c += step_s
    end
    isempty(centers) && throw(ArgumentError(
        "flight_density_table: no window produced a valid (finite, positive) density."
    ))
    return DataFrame(time_s=centers, rho_kgm3=rhos)
end
