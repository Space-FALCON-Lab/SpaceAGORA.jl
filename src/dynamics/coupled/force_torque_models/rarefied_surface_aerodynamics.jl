"""
    free_molecular_surface_coefficients(
        speed_ratio, incidence, wall_to_freestream_temperature_ratio,
        normal_accommodation, tangential_accommodation,
    ) -> (pressure_coefficient, tangential_multiplier)

Evaluate the Schaaf-Chambre free-molecular flat-surface pressure and shear
law reproduced on page 133 of NASA CR-182076. `incidence` is the
dot product of the freestream direction and the inward unit normal. The
returned tangential multiplier multiplies the tangential projection of the
unit freestream vector; this retains the complete shear vector without a
grazing-incidence division.

Both accommodation coefficients use the conventional range `[0, 1]`: zero
is specular/no tangential accommodation and one is fully accommodated.
"""
function free_molecular_surface_coefficients(
    speed_ratio::Real,
    incidence::Real,
    wall_to_freestream_temperature_ratio::Real,
    normal_accommodation::Real,
    tangential_accommodation::Real,
)
    S = Float64(speed_ratio)
    μ = Float64(incidence)
    temperature_ratio = Float64(wall_to_freestream_temperature_ratio)
    fn = Float64(normal_accommodation)
    ft = Float64(tangential_accommodation)
    isfinite(S) && S > 0.0 || throw(DomainError(
        speed_ratio,
        "molecular speed ratio must be finite and positive.",
    ))
    isfinite(μ) && -1.0 - 8 * eps(Float64) <= μ <= 1.0 + 8 * eps(Float64) || throw(
        DomainError(incidence, "surface incidence must be finite and lie in [-1, 1]."),
    )
    isfinite(temperature_ratio) && temperature_ratio > 0.0 || throw(DomainError(
        wall_to_freestream_temperature_ratio,
        "wall/freestream temperature ratio must be finite and positive.",
    ))
    isfinite(fn) && 0.0 <= fn <= 1.0 || throw(DomainError(
        normal_accommodation,
        "normal accommodation must be finite and lie in [0, 1].",
    ))
    isfinite(ft) && 0.0 <= ft <= 1.0 || throw(DomainError(
        tangential_accommodation,
        "tangential accommodation must be finite and lie in [0, 1].",
    ))

    μ = clamp(μ, -1.0, 1.0)
    normal_speed_ratio = S * μ
    exponential = exp(-normal_speed_ratio^2)
    cumulative = erfc(-normal_speed_ratio)
    sqrt_temperature_ratio = sqrt(temperature_ratio)

    first_pressure_factor = (
        (2.0 - fn) * normal_speed_ratio / sqrt(π) +
        0.5 * fn * sqrt_temperature_ratio
    )
    second_pressure_factor = (
        (2.0 - fn) * (normal_speed_ratio^2 + 0.5) +
        0.5 * fn * sqrt(π) * sqrt_temperature_ratio * normal_speed_ratio
    )
    pressure_numerator = if normal_speed_ratio < 0.0
        exponential * (
            first_pressure_factor +
            second_pressure_factor * erfcx(-normal_speed_ratio)
        )
    else
        first_pressure_factor * exponential + second_pressure_factor * cumulative
    end
    pressure_coefficient = pressure_numerator / S^2

    molecular_flux_factor = if normal_speed_ratio < 0.0
        # This form avoids subtracting two nearly equal, separately rounded
        # exponentially small quantities on a leeward surface.
        exponential * (
            1.0 + sqrt(π) * normal_speed_ratio * erfcx(-normal_speed_ratio)
        )
    else
        exponential + sqrt(π) * normal_speed_ratio * cumulative
    end
    tangential_multiplier = ft * molecular_flux_factor / (sqrt(π) * S)
    return pressure_coefficient, tangential_multiplier
end

"""
    free_molecular_aerodynamic_coefficients(
        geometry, freestream_unit_body;
        speed_ratio, temperature_inf_k, wall_temperature_k,
        normal_accommodation=1, tangential_accommodation=1,
    )

Integrate the complete free-molecular pressure and tangential shear law over
every element of an [`AerodynamicGeometry`](@ref), including force and moment
about its configured reference point. Unlike Newtonian theory, leeward and
base surfaces are evaluated because the molecular velocity distribution has a
finite thermal width.

The geometry represents one gas-facing side per surface sample. A thin plate
must therefore be represented by two oppositely oriented surfaces when both
sides are exposed.
"""
function free_molecular_aerodynamic_coefficients(
    geometry::AerodynamicGeometry,
    freestream_unit_body::AbstractVector{<:Real};
    speed_ratio::Real,
    temperature_inf_k::Real,
    wall_temperature_k::Real,
    normal_accommodation::Real=1.0,
    tangential_accommodation::Real=1.0,
)
    freestream = SVector{3, Float64}(freestream_unit_body)
    all(isfinite, freestream) || throw(DomainError(
        freestream_unit_body,
        "freestream direction must be finite.",
    ))
    freestream_magnitude = norm(freestream)
    freestream_magnitude > eps(Float64) || throw(DomainError(
        freestream_unit_body,
        "freestream direction must be nonzero.",
    ))
    freestream /= freestream_magnitude
    temperature_inf = Float64(temperature_inf_k)
    wall_temperature = Float64(wall_temperature_k)
    isfinite(temperature_inf) && temperature_inf > 0.0 || throw(DomainError(
        temperature_inf_k,
        "freestream temperature must be finite and positive.",
    ))
    isfinite(wall_temperature) && wall_temperature > 0.0 || throw(DomainError(
        wall_temperature_k,
        "wall temperature must be finite and positive.",
    ))
    temperature_ratio = wall_temperature / temperature_inf

    force_integral = MVector{3, Float64}(0.0, 0.0, 0.0)
    moment_integral = MVector{3, Float64}(0.0, 0.0, 0.0)
    surface = geometry.surface
    @inbounds for index in eachindex(surface.area_weights_m2)
        normal = surface.inward_normals_body[index]
        incidence = dot(freestream, normal)
        pressure_coefficient, tangential_multiplier =
            free_molecular_surface_coefficients(
                speed_ratio,
                incidence,
                temperature_ratio,
                normal_accommodation,
                tangential_accommodation,
            )
        tangential_freestream = freestream - incidence * normal
        local_force_coefficient =
            pressure_coefficient * normal +
            tangential_multiplier * tangential_freestream
        weighted_force = surface.area_weights_m2[index] * local_force_coefficient
        force_integral .+= weighted_force
        moment_integral .+= cross(surface.positions_body_m[index], weighted_force)
    end

    force_body = SVector{3, Float64}(force_integral) / geometry.reference_area_m2
    moment_body = SVector{3, Float64}(moment_integral) /
                  (geometry.reference_area_m2 * geometry.reference_length_m)
    return _surface_aerodynamic_coefficients(
        force_body,
        moment_body,
        freestream,
        1.0,
    )
end

function free_molecular_aerodynamic_coefficients(
    geometry::AerodynamicGeometry,
    alpha_rad::Real,
    beta_rad::Real;
    kwargs...,
)
    alpha = Float64(alpha_rad)
    beta = Float64(beta_rad)
    isfinite(alpha) || throw(DomainError(alpha_rad, "angle of attack must be finite."))
    isfinite(beta) || throw(DomainError(beta_rad, "sideslip must be finite."))
    freestream = SVector{3, Float64}(
        -cos(alpha) * cos(beta),
        -sin(beta),
        -sin(alpha) * cos(beta),
    )
    return free_molecular_aerodynamic_coefficients(geometry, freestream; kwargs...)
end

"""
    gas_mean_free_path(density_kg_m3, temperature_k, gas_constant_j_kg_k,
                       dynamic_viscosity_pa_s)

Return the kinetic-theory mean free path
`lambda = mu/rho * sqrt(pi/(2*R*T))`. Dynamic viscosity is an explicit input;
the function does not silently impose a temperature-viscosity correlation.
"""
function gas_mean_free_path(
    density_kg_m3::Real,
    temperature_k::Real,
    gas_constant_j_kg_k::Real,
    dynamic_viscosity_pa_s::Real,
)
    density = Float64(density_kg_m3)
    temperature = Float64(temperature_k)
    gas_constant = Float64(gas_constant_j_kg_k)
    viscosity = Float64(dynamic_viscosity_pa_s)
    isfinite(density) && density > 0.0 || throw(DomainError(
        density_kg_m3,
        "gas density must be finite and positive.",
    ))
    isfinite(temperature) && temperature > 0.0 || throw(DomainError(
        temperature_k,
        "gas temperature must be finite and positive.",
    ))
    isfinite(gas_constant) && gas_constant > 0.0 || throw(DomainError(
        gas_constant_j_kg_k,
        "specific gas constant must be finite and positive.",
    ))
    isfinite(viscosity) && viscosity > 0.0 || throw(DomainError(
        dynamic_viscosity_pa_s,
        "dynamic viscosity must be finite and positive.",
    ))
    return viscosity / density * sqrt(π / (2.0 * gas_constant * temperature))
end

"""
    gas_knudsen_number(density_kg_m3, temperature_k, gas_constant_j_kg_k,
                       dynamic_viscosity_pa_s, characteristic_length_m)

Compute `Kn = lambda/L` using [`gas_mean_free_path`](@ref).
"""
function gas_knudsen_number(
    density_kg_m3::Real,
    temperature_k::Real,
    gas_constant_j_kg_k::Real,
    dynamic_viscosity_pa_s::Real,
    characteristic_length_m::Real,
)
    length_scale = Float64(characteristic_length_m)
    isfinite(length_scale) && length_scale > 0.0 || throw(DomainError(
        characteristic_length_m,
        "Knudsen characteristic length must be finite and positive.",
    ))
    return gas_mean_free_path(
        density_kg_m3,
        temperature_k,
        gas_constant_j_kg_k,
        dynamic_viscosity_pa_s,
    ) / length_scale
end

"""
    transitional_free_molecular_weight(
        knudsen_number; continuum_limit=1e-3, free_molecular_limit=10,
    )

Return the Wilmoth-Blanchard-Moss sine-squared bridge weight on a logarithmic
Knudsen-number coordinate. The result is exactly zero below the configured
continuum anchor and exactly one above the free-molecular anchor, with zero
slope at both anchors. The defaults reproduce their Shuttle bridging form;
the limits remain configurable because the paper explicitly finds no universal
transition correlation.
"""
function transitional_free_molecular_weight(
    knudsen_number::Real;
    continuum_limit::Real=1.0e-3,
    free_molecular_limit::Real=10.0,
)
    Kn = Float64(knudsen_number)
    Kn_continuum = Float64(continuum_limit)
    Kn_free_molecular = Float64(free_molecular_limit)
    !isnan(Kn) && Kn >= 0.0 || throw(DomainError(
        knudsen_number,
        "Knudsen number must be nonnegative and not NaN.",
    ))
    isfinite(Kn_continuum) && Kn_continuum > 0.0 || throw(DomainError(
        continuum_limit,
        "continuum bridge limit must be finite and positive.",
    ))
    isfinite(Kn_free_molecular) && Kn_free_molecular > Kn_continuum || throw(DomainError(
        free_molecular_limit,
        "free-molecular bridge limit must be finite and greater than the continuum limit.",
    ))
    Kn <= Kn_continuum && return 0.0
    Kn >= Kn_free_molecular && return 1.0
    logarithmic_fraction = (
        log10(Kn) - log10(Kn_continuum)
    ) / (
        log10(Kn_free_molecular) - log10(Kn_continuum)
    )
    return sinpi(0.5 * logarithmic_fraction)^2
end

"""
    AerodynamicRegimeResult

Surface-integrated aerodynamic coefficients together with the actual regime,
Knudsen number, and free-molecular bridge weight used to obtain them.
"""
struct AerodynamicRegimeResult
    coefficients::SurfaceAerodynamicCoefficients
    regime::Symbol
    knudsen_number::Union{Nothing, Float64}
    free_molecular_weight::Float64
end

@inline function _blend_surface_coefficients(
    continuum::SurfaceAerodynamicCoefficients,
    free_molecular::SurfaceAerodynamicCoefficients,
    freestream::SVector{3, Float64},
    free_molecular_weight::Float64,
)
    continuum_weight = 1.0 - free_molecular_weight
    return _surface_aerodynamic_coefficients(
        continuum_weight * continuum.force_body +
            free_molecular_weight * free_molecular.force_body,
        continuum_weight * continuum.moment_body +
            free_molecular_weight * free_molecular.moment_body,
        freestream,
        continuum_weight * continuum.exposed_area_fraction +
            free_molecular_weight * free_molecular.exposed_area_fraction,
    )
end

"""
    aerodynamic_regime_coefficients(
        geometry, freestream_unit_body;
        flow_regime=:automatic, knudsen_number=nothing,
        continuum_limit=1e-3, free_molecular_limit=10,
        pressure_model=:regular_newtonian, mach=nothing, gamma=nothing,
        speed_ratio=nothing, temperature_inf_k=nothing,
        wall_temperature_k=nothing, normal_accommodation=1,
        tangential_accommodation=1,
    )

Evaluate continuum Newtonian, free-molecular Schaaf-Chambre, or a transitional
blend on one shared geometry. Transitional force and moment vectors are
blended before body- and wind-axis coefficients are derived, so all reported
coefficients remain mutually consistent.
"""
function aerodynamic_regime_coefficients(
    geometry::AerodynamicGeometry,
    freestream_unit_body::AbstractVector{<:Real};
    flow_regime::Symbol=:automatic,
    knudsen_number::Union{Nothing, Real}=nothing,
    continuum_limit::Real=1.0e-3,
    free_molecular_limit::Real=10.0,
    pressure_model::Symbol=:regular_newtonian,
    mach::Union{Nothing, Real}=nothing,
    gamma::Union{Nothing, Real}=nothing,
    speed_ratio::Union{Nothing, Real}=nothing,
    temperature_inf_k::Union{Nothing, Real}=nothing,
    wall_temperature_k::Union{Nothing, Real}=nothing,
    normal_accommodation::Real=1.0,
    tangential_accommodation::Real=1.0,
)
    flow_regime in (:automatic, :continuum, :transitional, :free_molecular) || throw(
        ArgumentError(
            "unsupported aerodynamic flow regime $(repr(flow_regime)); expected :automatic, :continuum, :transitional, or :free_molecular.",
        ),
    )
    freestream = SVector{3, Float64}(freestream_unit_body)
    all(isfinite, freestream) || throw(DomainError(
        freestream_unit_body,
        "freestream direction must be finite.",
    ))
    freestream_magnitude = norm(freestream)
    freestream_magnitude > eps(Float64) || throw(DomainError(
        freestream_unit_body,
        "freestream direction must be nonzero.",
    ))
    freestream /= freestream_magnitude

    Kn = if knudsen_number === nothing
        flow_regime === :continuum || flow_regime === :free_molecular || throw(
            ArgumentError("automatic and transitional flow selection require `knudsen_number`."),
        )
        nothing
    else
        value = Float64(knudsen_number)
        !isnan(value) && value >= 0.0 || throw(DomainError(
            knudsen_number,
            "Knudsen number must be nonnegative and not NaN.",
        ))
        value
    end

    free_molecular_weight = if flow_regime === :continuum
        0.0
    elseif flow_regime === :free_molecular
        1.0
    else
        transitional_free_molecular_weight(
            something(Kn);
            continuum_limit=continuum_limit,
            free_molecular_limit=free_molecular_limit,
        )
    end

    if free_molecular_weight == 0.0
        coefficients = newtonian_aerodynamic_coefficients(
            geometry,
            freestream;
            pressure_model=pressure_model,
            mach=mach,
            gamma=gamma,
        )
        return AerodynamicRegimeResult(coefficients, :continuum, Kn, 0.0)
    end

    speed_ratio === nothing && throw(ArgumentError(
        "free-molecular and transitional aerodynamics require `speed_ratio`.",
    ))
    temperature_inf_k === nothing && throw(ArgumentError(
        "free-molecular and transitional aerodynamics require `temperature_inf_k`.",
    ))
    wall_temperature_k === nothing && throw(ArgumentError(
        "free-molecular and transitional aerodynamics require `wall_temperature_k`.",
    ))
    free_molecular_coefficients = free_molecular_aerodynamic_coefficients(
        geometry,
        freestream;
        speed_ratio=speed_ratio,
        temperature_inf_k=temperature_inf_k,
        wall_temperature_k=wall_temperature_k,
        normal_accommodation=normal_accommodation,
        tangential_accommodation=tangential_accommodation,
    )
    free_molecular_weight == 1.0 && return AerodynamicRegimeResult(
        free_molecular_coefficients,
        :free_molecular,
        Kn,
        1.0,
    )

    continuum_coefficients = newtonian_aerodynamic_coefficients(
        geometry,
        freestream;
        pressure_model=pressure_model,
        mach=mach,
        gamma=gamma,
    )
    coefficients = _blend_surface_coefficients(
        continuum_coefficients,
        free_molecular_coefficients,
        freestream,
        free_molecular_weight,
    )
    return AerodynamicRegimeResult(
        coefficients,
        :transitional,
        Kn,
        free_molecular_weight,
    )
end

function aerodynamic_regime_coefficients(
    geometry::AerodynamicGeometry,
    alpha_rad::Real,
    beta_rad::Real;
    kwargs...,
)
    alpha = Float64(alpha_rad)
    beta = Float64(beta_rad)
    isfinite(alpha) || throw(DomainError(alpha_rad, "angle of attack must be finite."))
    isfinite(beta) || throw(DomainError(beta_rad, "sideslip must be finite."))
    freestream = SVector{3, Float64}(
        -cos(alpha) * cos(beta),
        -sin(beta),
        -sin(alpha) * cos(beta),
    )
    return aerodynamic_regime_coefficients(geometry, freestream; kwargs...)
end
