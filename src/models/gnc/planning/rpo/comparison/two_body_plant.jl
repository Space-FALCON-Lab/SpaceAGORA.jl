using ..Planets: Earth
# Recovered from SpaceAGORA_RL_NEW 5335ae7c, adapted without the RL dependency.
"""Circular Earth reference orbit matching the HCW prediction mean motion."""
function rpo_two_body_reference_orbit(mean_motion_radps::Real)
    n = Float64(mean_motion_radps)
    isfinite(n) && n > 0.0 || throw(ArgumentError("RPO mean motion must be finite and positive."))
    μ = Earth().μ
    radius = cbrt(μ / n^2)
    return (μ=μ, radius_m=radius,
        r_target=SVector(radius, 0.0, 0.0),
        v_target=SVector(0.0, sqrt(μ / radius), 0.0))
end

"""Inertial target/chaser states for the nonlinear comparison plant."""
mutable struct RPOTwoBodyPlant
    state::SVector{12, Float64}
    μ::Float64
end

function rpo_init_two_body_plant(relative_state, mean_motion_radps::Real)
    orbit = rpo_two_body_reference_orbit(mean_motion_radps)
    r_chaser, v_chaser = rtn_to_inertial_relative_state(
        relative_state[1:3], relative_state[4:6], orbit.r_target, orbit.v_target,
    )
    return RPOTwoBodyPlant(
        SVector{12, Float64}(vcat(orbit.r_target, orbit.v_target, r_chaser, v_chaser)), orbit.μ,
    )
end

@inline function _rpo_two_body_derivative(state::SVector{12, Float64}, acceleration_rtn, μ, inertial::Bool=false)
    r_target = SVector(state[1], state[2], state[3])
    v_target = SVector(state[4], state[5], state[6])
    r_chaser = SVector(state[7], state[8], state[9])
    v_chaser = SVector(state[10], state[11], state[12])
    a_target = -μ * r_target / norm(r_target)^3
    a_chaser = -μ * r_chaser / norm(r_chaser)^3 +
        (inertial ? acceleration_rtn : rtn_dcm_from_inertial(r_target, v_target) * acceleration_rtn)
    return SVector{12, Float64}(vcat(v_target, a_target, v_chaser, a_chaser))
end

"""
Advance both inertial orbits with RK4, then return rotating-frame RTN feedback.
Acceleration defaults to a held RTN vector rotated at every RK stage.
`acceleration_frame=:inertial` instead holds its inertial direction, matching
the full simulator actuator zero-order hold. Actuator evaluators supply their interval-averaged
realized acceleration. The target geometry remains fixed in RTN.
"""
function rpo_step_two_body!(plant::RPOTwoBodyPlant, acceleration_rtn, dt_s::Real; max_step_s::Real=0.05, acceleration_frame::Symbol=:rtn)
    dt_s > 0.0 && max_step_s > 0.0 || throw(ArgumentError("RPO integration steps must be positive."))
    acceleration_frame in (:rtn,:inertial) || throw(ArgumentError("Unknown RPO acceleration frame."))
    inertial = acceleration_frame === :inertial
    steps = max(1, ceil(Int, dt_s / max_step_s))
    h = Float64(dt_s) / steps
    u = SVector{3, Float64}(acceleration_rtn)
    y = plant.state
    for _ in 1:steps
        k1 = _rpo_two_body_derivative(y, u, plant.μ, inertial)
        k2 = _rpo_two_body_derivative(y + (h / 2) * k1, u, plant.μ, inertial)
        k3 = _rpo_two_body_derivative(y + (h / 2) * k2, u, plant.μ, inertial)
        k4 = _rpo_two_body_derivative(y + h * k3, u, plant.μ, inertial)
        y += (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
    end
    plant.state = y
    return inertial_to_rtn_relative_state(y[7:9], y[10:12], y[1:3], y[4:6])
end
