struct NormalizationBounds
    names::Vector{Symbol}
    lower::Vector{Float64}
    upper::Vector{Float64}
end

function normalize_value(value::Real, lower::Real, upper::Real)
    lower_f = Float64(lower)
    upper_f = Float64(upper)
    upper_f == lower_f && return Float32(0)
    return Float32((Float64(value) - lower_f) / (upper_f - lower_f))
end

function paper_normalization_bounds(phase::AbstractString)
    constants = mars_odyssey_phase_constants(phase)
    names = paper_observation_names()
    lower = [400.0, constants.mars_radius_m, 85e3, 0.0, 1.57, 1.39, 730486.0, 0.0, 0.0]
    upper = [1200.0, constants.r_norm_m, 145e3, 2pi, 3.14, 1.75, 731216.0, 1.5e-7, 0.5]
    return NormalizationBounds(names, lower, upper)
end

function normalize_observation(values::AbstractVector{<:Real}, bounds::NormalizationBounds)
    length(values) == length(bounds.lower) || throw(DimensionMismatch("observation length does not match bounds"))
    return [normalize_value(values[i], bounds.lower[i], bounds.upper[i]) for i in eachindex(values)]
end
