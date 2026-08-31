#=
"""
    KS propagation model.

    Inverse-square gravity is represented through the KS energy parameter. J2
    and drag enter as perturbing accelerations.
"""
=#
function _skew_rotation(params::AerobrakingMPCParams)
    Ω = params.Ω
    return @SMatrix [0.0 -Ω 0.0; Ω 0.0 0.0; 0.0 0.0 0.0]
end

function j2_acceleration_si(rvec, params::AerobrakingMPCParams)
    r = norm(rvec)
    r <= eps(Float64) && return @SVector [0.0, 0.0, 0.0]
    scale = 3.0 * params.J2 * params.μ * params.Re^2 / (2.0 * r^5)
    z2_r2 = rvec[3]^2 / r^2
    return scale * @SVector [
        rvec[1] * (5.0 * z2_r2 - 1.0),
        rvec[2] * (5.0 * z2_r2 - 1.0),
        rvec[3] * (5.0 * z2_r2 - 3.0),
    ]
end

function drag_acceleration_si(
    rvec,
    vvec,
    params::AerobrakingMPCParams,
    area_m2::Real;
    density_kg_m3::Real,
    drag_coefficient::Real,
    mass_kg::Real,
)
    ρ = max(0.0, Float64(density_kg_m3))
    m = Float64(mass_kg)
    m > 0.0 || throw(ArgumentError("mass_kg must be positive."))
    Vrel = vvec - _skew_rotation(params) * rvec
    Vn = norm(Vrel)
    Vn <= eps(Float64) && return @SVector [0.0, 0.0, 0.0]
    return -0.5 * ρ * Float64(drag_coefficient) * Float64(area_m2) / m * Vn * Vrel
end

function cartesian_to_ks_state(position_ii_m, velocity_ii_m, params::AerobrakingMPCParams; elapsed_time_s::Real=0.0)
    p = SVector{4, Float64}(u_from_x(position_ii_m)...)
    q = SVector{4, Float64}(uprime_from_xdot(velocity_ii_m, p)...)
    r = Λ(p) * p
    v = xdot_from_u(p, q)
    energy = 0.5 * dot(v, v) - params.μ / norm(r)
    h = -2.0 * energy
    return collect(vcat(p, q, h, Float64(elapsed_time_s)))
end

function ks_state_to_cartesian(state)
    p = SVector{4, Float64}(state[1:4])
    q = SVector{4, Float64}(state[5:8])
    return (
        position_ii_m=Λ(p) * p,
        velocity_ii_m=xdot_from_u(p, q),
        energy_parameter=Float64(state[9]),
        elapsed_time_s=Float64(state[10]),
    )
end

function ks_rhs(
    state::AbstractVector,
    params::AerobrakingMPCParams,
    area_m2::Real;
    config::Union{Nothing, AerobrakingMPCConfig}=nothing,
    density_kg_m3::Real=0.0,
    use_drag::Bool,
)
    p = SVector{4, Float64}(state[1:4])
    q = SVector{4, Float64}(state[5:8])
    h = Float64(state[9])
    rvec = Λ(p) * p
    r = dot(p, p)
    vvec = xdot_from_u(p, q)
    a = j2_acceleration_si(rvec, params)
    if use_drag
        config === nothing && throw(ArgumentError("ks_rhs requires config with user-specified mass_kg when use_drag=true."))
        a += drag_acceleration_si(
            rvec,
            vvec,
            params,
            area_m2;
            density_kg_m3=density_kg_m3,
            drag_coefficient=config.drag_coefficient,
            mass_kg=config.mass_kg,
        )
    end
    ap4 = SVector(a[1], a[2], a[3], 0.0)
    dp = q
    dq = -0.25 * h .* p + 0.5 * r .* (transpose(ks_L(p)) * ap4)
    dh = -2.0 * r * dot(vvec, a)
    dt = r
    return collect(vcat(dp, dq, dh, dt))
end

function rk4_step(
    state::AbstractVector,
    params::AerobrakingMPCParams,
    area_m2::Real,
    Δs::Real;
    config::Union{Nothing, AerobrakingMPCConfig}=nothing,
    density_kg_m3::Real=0.0,
    use_drag::Bool,
)
    x = Float64.(collect(state))
    step = Float64(Δs)
    f(z) = ks_rhs(
        z,
        params,
        area_m2;
        config=config,
        density_kg_m3=density_kg_m3,
        use_drag=use_drag,
    )
    k1 = f(x)
    k2 = f(x .+ 0.5 * step .* k1)
    k3 = f(x .+ 0.5 * step .* k2)
    k4 = f(x .+ step .* k3)
    return x .+ (step / 6.0) .* (k1 .+ 2.0 .* k2 .+ 2.0 .* k3 .+ k4)
end
