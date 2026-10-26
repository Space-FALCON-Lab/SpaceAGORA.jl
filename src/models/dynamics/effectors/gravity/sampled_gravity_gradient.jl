@inline function _gravity_gradient_torque_body(
    model,
    x::StateSample,
    planet
)::SVector{3, Float64}
    if !model.gravity_gradient || x.q_ib === nothing || x.spacecraft === nothing
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    r = norm(x.pos_ii)
    if !isfinite(r) || r <= 0.0
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    inertia_tensor = x.spacecraft.inertia_tensor
    r_body = rot(x.q_ib) * x.pos_ii
    return gravity_gradient(inertia_tensor, r_body, planet.μ)
end
