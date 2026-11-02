## Splitting integration
#
# The KS right-hand side is split into the unperturbed oscillator with t′ = |u|²
# (flow A) and a position-only perturbation (flow B). Both flows are exact, so
# the step error comes only from their noncommutativity. With p frozen, B moves
# u′ by the constant g = (|u|²/2) Λ(u)ᵀ a and h_KS by -2 aᵀΛ(u) ∫ u′ ds.

const KS_SPLITTING_WEIGHTS = (
    order2=(1.0,),
    order4=let w1 = 1.0 / (2.0 - cbrt(2.0))
        (w1, 1.0 - 2.0 * w1, w1)
    end,
)

function ks_splitting_weights(order::Integer)
    order == 2 && return KS_SPLITTING_WEIGHTS.order2
    order == 4 && return KS_SPLITTING_WEIGHTS.order4
    throw(ArgumentError("KS splitting order must be 2 or 4."))
end

"""
    ks_oscillator_flow(u, u_prime, h_ks, delta_s) -> (u, u_prime, delta_t)

Exact flow of `u′′ = -(h_KS/2) u`, `t′ = |u|²` over `delta_s` at fixed `h_KS > 0`.
"""
@inline function ks_oscillator_flow(u, u_prime, h_ks, delta_s)
    h_ks > 0.0 || throw(ArgumentError("The KS oscillator flow requires h_KS > 0."))
    omega = sqrt(0.5 * h_ks)
    sine, cosine = sincos(omega * delta_s)
    u_new = cosine * u + (sine / omega) * u_prime
    u_prime_new = (-omega * sine) * u + cosine * u_prime
    half_sin_double = sine * cosine / omega
    delta_t = 0.5 * dot(u, u) * (delta_s + half_sin_double) +
        dot(u_prime, u_prime) / (2.0 * omega^2) * (delta_s - half_sin_double) +
        dot(u, u_prime) / omega^2 * sine^2
    return u_new, u_prime_new, delta_t
end

"""
    ks_j2_kick_terms(u, params) -> (Λᵀa, g)

Evaluate the J2 acceleration `a` once at the KS coordinate `u` and return
`Λ(u)ᵀ a` and `g = (|u|²/2) Λ(u)ᵀ a`.
"""
@inline function ks_j2_kick_terms(u, params)
    u4 = _ks_coordinate(u)
    radius = dot(u4, u4)
    lambda_t_a = transpose(ks_lambda_matrix(u4)) *
        ks_j2_acceleration_with_radius(ks_position(u4), radius, params)
    return lambda_t_a, (0.5 * radius) * lambda_t_a
end

"""
    ks_perturbation_kick(u_prime, h_ks, kick_terms, tau) -> (u_prime, h_ks)

Exact flow over `tau` of `u′′ = g`, `h_KS′ = -2 aᵀΛ(u) u′` at frozen `u`.
"""
@inline function ks_perturbation_kick(u_prime, h_ks, kick_terms, tau)
    lambda_t_a, g = kick_terms
    return u_prime + tau * g,
        h_ks - 2.0 * dot(lambda_t_a, tau * u_prime + (0.5 * tau^2) * g)
end

"""
    ks_splitting_step(state, kick_terms, params, delta_s; order=2) -> (state, kick_terms)

Advance the KS state by one J2-only splitting step: the Strang composition
`B(Δs/2) ∘ A(Δs) ∘ B(Δs/2)` (`order=2`) or the Yoshida triple jump of three Strang
steps (`order=4`). `kick_terms` must belong to `state` (see
[`ks_j2_kick_terms`](@ref)); the returned pair belongs to the new state, so a
step costs one (order 2) or three (order 4) perturbation evaluations.
"""
function ks_splitting_step(state::AbstractVector, kick_terms, params,
    delta_s::Real; order::Integer=2)
    step = Float64(delta_s)
    step > 0.0 || throw(ArgumentError("delta_s must be positive."))
    u = SVector{4,Float64}(state[1], state[2], state[3], state[4])
    u_prime = SVector{4,Float64}(state[5], state[6], state[7], state[8])
    h_ks = Float64(state[9])
    t = Float64(state[10])
    for weight in ks_splitting_weights(order)
        stage = weight * step
        u_prime, h_ks = ks_perturbation_kick(u_prime, h_ks, kick_terms, 0.5 * stage)
        u, u_prime, delta_t = ks_oscillator_flow(u, u_prime, h_ks, stage)
        t += delta_t
        kick_terms = ks_j2_kick_terms(u, params)
        u_prime, h_ks = ks_perturbation_kick(u_prime, h_ks, kick_terms, 0.5 * stage)
    end
    return SVector{10,Float64}(u..., u_prime..., h_ks, t), kick_terms
end

"""
    ks_splitting_step_size(state, steps_per_orbit) -> Float64

Fictitious-time step `π/(ω N)`, with `ω² = h_KS/2`; one orbit spans `π/ω` in `s`.
"""
ks_splitting_step_size(state::AbstractVector, steps_per_orbit::Real) =
    pi / (sqrt(0.5 * Float64(state[9])) * Float64(steps_per_orbit))

"""
    ks_splitting_propagate(state, params, final_time_s, steps_per_orbit;
        order=2, save_history=true) -> (states, work)

Propagate with fixed J2-only splitting steps until the physical time reaches
`final_time_s`; the last step is taken in full. `states` holds every KS state
(only the initial one without `save_history`). `work` holds the step count,
`delta_s`, the perturbation evaluations including the initial one, and the
final state.
"""
function ks_splitting_propagate(state::AbstractVector, params, final_time_s::Real,
    steps_per_orbit::Real; order::Integer=2, save_history::Bool=true)
    x = SVector{10,Float64}(state)
    delta_s = ks_splitting_step_size(x, steps_per_orbit)
    stages = length(ks_splitting_weights(order))
    kick_terms = ks_j2_kick_terms(x, params)
    states = [x]
    steps = 0
    while x[10] < final_time_s
        x, kick_terms = ks_splitting_step(x, kick_terms, params, delta_s; order)
        steps += 1
        save_history && push!(states, x)
    end
    return states, (; steps, delta_s, evaluations=1 + stages * steps, final_state=x)
end
