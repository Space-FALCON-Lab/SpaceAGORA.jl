@inline function acceleration_from_force(
    net_force::AbstractVector{<:Real},
    mass::Real,
)::SVector{3, Float64}
    mass_f64 = Float64(mass)
    if !isfinite(mass_f64) || abs(mass_f64) <= eps(Float64)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    inv_mass = inv(mass_f64)
    return SVector{3, Float64}(
        Float64(net_force[1]) * inv_mass,
        Float64(net_force[2]) * inv_mass,
        Float64(net_force[3]) * inv_mass,
    )
end

@inline function mass_derivative(mass_rate::Real)::Float64
    return isfinite(mass_rate) ? Float64(mass_rate) : 0.0
end

@inline function assign_full_translational_rhs!(
    du_view,
    sc_view,
    net_force::AbstractVector{<:Real},
    mass_rate::Real,
)
    du_view.pos .= position_derivative(sc_view.vel)
    du_view.vel .= acceleration_from_force(net_force, sc_view.mass)
    du_view.mass = mass_derivative(mass_rate)
    return nothing
end

@inline function assign_slow_translational_rhs!(
    du_view,
    sc_view,
    net_force::AbstractVector{<:Real},
)
    du_view.pos .= position_derivative(sc_view.vel)
    du_view.vel .= acceleration_from_force(net_force, sc_view.mass)
    du_view.mass = 0.0
    return nothing
end

@inline function assign_control_only_translational_rhs!(
    du_view,
    sc_view,
    net_force::AbstractVector{<:Real},
    mass_rate::Real,
)
    du_view.pos .= zero_position_derivative()
    du_view.vel .= acceleration_from_force(net_force, sc_view.mass)
    du_view.mass = mass_derivative(mass_rate)
    return nothing
end

@inline function assign_force_only_translational_rhs!(
    du_view,
    sc_view,
    net_force::AbstractVector{<:Real},
)
    du_view.pos .= zero_position_derivative()
    du_view.vel .= acceleration_from_force(net_force, sc_view.mass)
    du_view.mass = 0.0
    return nothing
end
