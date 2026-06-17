Base.@kwdef struct AerobrakingRandomizationConfig
    nominal::Bool = true
    apoapsis_jitter_m::Float64 = 2.5e3
    periapsis_jitter_m::Float64 = 1.0e3
    angle_jitter_deg::Float64 = 0.25
    nonnominal_inclination_low_deg::Float64 = 88.6
    nonnominal_inclination_high_deg::Float64 = 98.6
    nonnominal_angle_low_deg::Float64 = 110.0
    nonnominal_angle_high_deg::Float64 = 120.0
    process_noise::Bool = false
    process_noise_scale::Float64 = 0.0
end

uniform_jitter(rng::AbstractRNG, span::Real) = (2rand(rng) - 1) * Float64(span)
