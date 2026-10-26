function aerobraking_gravity_force_ii(
    gm_code::Integer,
    mass::Float64,
    pos_ii::SVector{3, Float64},
    vel_ii::SVector{3, Float64},
    pos_pp::SVector{3, Float64},
    lat::Float64,
    lon::Float64,
    alt::Float64,
    planet,
    et::Float64,
    args,
    gram_atmosphere,
    gram,
    n_bodies_list,
)::SVector{3, Float64}
    gravity_ii = if gm_code == 2
        _inverse_squared_j2_gravity_accel(pos_ii, planet)
    else
        # The legacy gm=0/1/3 aerobraking paths all reduced to inverse-squared gravity.
        _inverse_squared_gravity_accel(pos_ii, planet)
    end

    force_ii = mass * gravity_ii

    n_bodies = _gravity_runtime_field(args, :n_bodies, ())
    if !isempty(n_bodies)
        for k in eachindex(n_bodies)
            force_ii += mass * PerturbationEffectors.gravity_n_bodies(et, pos_ii, planet, n_bodies_list[k])
        end
    end

    if _gravity_runtime_field(args, :gravity_harmonics, 0) == 1
        L = Int(_gravity_runtime_field(args, :L, 0))
        M = Int(_gravity_runtime_field(args, :M, 0))
        force_ii += mass * planet.L_PI' * PerturbationEffectors.acc_gravity_pines!(pos_pp, planet.Clm, planet.Slm, L, M, planet.μ, planet.Rp_e, planet)
    end

    return force_ii
end

@inline function wrench(
    model::ConstantGravityModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    force_ii = _inverse_squared_force_ii(x.pos_ii, x.mass_kg, env.planet)
    torque_body = _gravity_gradient_torque_body(model, x, env.planet)
    return force_ii, torque_body
end

@inline gravity_backbone_structure(::ConstantGravityModel) = :position_only_static_gravity

@inline function gravity_backbone_acceleration_ii(
    model::ConstantGravityModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::SVector{3, Float64}
    return _inverse_squared_gravity_accel(x.pos_ii, env.planet)
end

@inline function gravity_backbone_relative_acceleration_ii(
    model::ConstantGravityModel, x_base::StateSample, ρ, x_far::StateSample, env_base::EnvironmentSample, env_far::EnvironmentSample, t::Float64,
)::SVector{3, Float64}
    return encke_point_mass_difference(Float64(env_base.planet.μ), x_base.pos_ii, ρ)
end

@inline function wrench(
    model::InverseSquaredGravityModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    force_ii = _inverse_squared_force_ii(x.pos_ii, x.mass_kg, env.planet)
    torque_body = _gravity_gradient_torque_body(model, x, env.planet)
    return force_ii, torque_body
end

@inline gravity_backbone_structure(::InverseSquaredGravityModel) = :position_only_static_gravity

@inline function gravity_backbone_acceleration_ii(
    model::InverseSquaredGravityModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::SVector{3, Float64}
    return _inverse_squared_gravity_accel(x.pos_ii, env.planet)
end

@inline function gravity_backbone_relative_acceleration_ii(
    model::InverseSquaredGravityModel, x_base::StateSample, ρ, x_far::StateSample, env_base::EnvironmentSample, env_far::EnvironmentSample, t::Float64,
)::SVector{3, Float64}
    return encke_point_mass_difference(Float64(env_base.planet.μ), x_base.pos_ii, ρ)
end

@inline environment_requirements(::InverseSquaredJ2GravityModel) = EffectorEnvironmentRequirements(planet_frame=true)

@inline function wrench(
    model::InverseSquaredJ2GravityModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    planet_frame = env.planet_frame
    planet_frame === nothing && throw(ArgumentError("InverseSquaredJ2GravityModel wrench requires env.planet_frame."))
    force_ii = _inverse_squared_j2_force_ii(x.pos_ii, x.mass_kg, planet_frame.l_pi, env.planet)
    torque_body = _gravity_gradient_torque_body(model, x, env.planet)
    return force_ii, torque_body
end

@inline gravity_backbone_structure(::InverseSquaredJ2GravityModel) = :position_only_static_gravity

@inline function gravity_backbone_acceleration_ii(
    model::InverseSquaredJ2GravityModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::SVector{3, Float64}
    planet_frame = env.planet_frame
    planet_frame === nothing && throw(ArgumentError("InverseSquaredJ2GravityModel gravity-backbone acceleration requires env.planet_frame."))
    gravity_pp = _inverse_squared_j2_gravity_accel(planet_frame.pos_pp, env.planet)
    return planet_frame.l_pi' * gravity_pp
end

@inline function gravity_backbone_relative_acceleration_ii(
    model::InverseSquaredJ2GravityModel, x_base::StateSample, ρ, x_far::StateSample, env_base::EnvironmentSample, env_far::EnvironmentSample, t::Float64,
)::SVector{3, Float64}
    pb = env_base.planet_frame
    pf = env_far.planet_frame
    (pb === nothing || pf === nothing) && throw(ArgumentError("InverseSquaredJ2GravityModel relative gravity requires env.planet_frame."))
    planet = env_base.planet
    point_mass = encke_point_mass_difference(Float64(planet.μ), x_base.pos_ii, ρ)
    # J2 is a small term (1e-3 of the point mass), so its plain difference is accurate enough.
    return point_mass + pb.l_pi' * (_j2_only_gravity_accel(pf.pos_pp, planet) - _j2_only_gravity_accel(pb.pos_pp, planet))
end
