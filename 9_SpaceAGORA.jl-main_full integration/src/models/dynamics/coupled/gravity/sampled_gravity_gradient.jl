@inline function _gravity_gradient_torque_body(
    model,
    pos_ii::SVector{3, Float64},
    x,
    param::ODEParams,
    i::Int64
)::SVector{3, Float64}
    if !model.gravity_gradient || !param.args.mission_configuration.orientation_sim
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    if i < 1 || i > length(param.args.dynamics_model.spacecraft)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    if !hasproperty(x, :q)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    r = norm(pos_ii)
    if !isfinite(r) || r <= 0.0
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    inertia_tensor = param.args.dynamics_model.spacecraft[i].inertia_tensor
    q_ib = getproperty(x, :q)
    q_body = SVector{4, Float64}(Float64(q_ib[1]), Float64(q_ib[2]), Float64(q_ib[3]), Float64(q_ib[4]))
    r_body = rot(q_body) * pos_ii
    return gravity_gradient(inertia_tensor, r_body, param.args.environment_model.planet.μ)
end

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
