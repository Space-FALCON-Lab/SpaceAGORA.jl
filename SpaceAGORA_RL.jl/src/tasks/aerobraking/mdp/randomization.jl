Base.@kwdef struct AerobrakingRandomizationConfig
    nominal::Bool = true
    apoapsis_jitter_m::Float64 = 2.5e3
    periapsis_jitter_m::Float64 = 2.5e3
    angle_jitter_deg::Float64 = 0.25
    nonnominal_inclination_low_deg::Float64 = 88.6
    nonnominal_inclination_high_deg::Float64 = 98.6
    nonnominal_aop_low_deg::Float64 = 60.0
    nonnominal_aop_high_deg::Float64 = 90.0
    nonnominal_raan_low_deg::Float64 = 110.0
    nonnominal_raan_high_deg::Float64 = 120.0
    initial_date_start::Union{Nothing,Date} = nothing
    initial_date_days::Int = 0
    randomize_initial_time_of_day::Bool = true
    initial_true_anomaly_jitter_deg::Float64 = 0.0
    process_noise::Bool = false
    process_noise_scale::Float64 = 0.0
    aerodynamic_coefficient_dispersion::Bool = false
    aerodynamic_coefficient_span::Float64 = 0.10
    aerodynamic_cd_span::Float64 = NaN
    aerodynamic_cl_span::Float64 = NaN
    marsgram_perturbation_scale::Float64 = 1.0
    marsgram_seed_base::Int = 1001
end

uniform_jitter(rng::AbstractRNG, span::Real) = (2rand(rng) - 1) * Float64(span)
randomized_cd_span(config::AerobrakingRandomizationConfig) =
    isfinite(config.aerodynamic_cd_span) ? config.aerodynamic_cd_span : config.aerodynamic_coefficient_span
randomized_cl_span(config::AerobrakingRandomizationConfig) =
    isfinite(config.aerodynamic_cl_span) ? config.aerodynamic_cl_span : config.aerodynamic_coefficient_span
