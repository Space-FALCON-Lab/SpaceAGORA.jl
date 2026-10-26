"""
    GravityGradientTorqueModel(; gravity_gradient=true)

Central-field gravity-gradient torque as a separate effector, with zero
translational force. The body-frame torque is `3μ/r³ (r̂_body × J r̂_body)`,
using the spacecraft inertia `J`, planet parameter `μ`, and body-frame radial
unit vector. This reuses the same calculation as `gravity_gradient=true` on
an analytic gravity model; it does not include nonspherical-field corrections.

Use with `orientation_sim=true`. The state-sample hook returns zero torque
when the attitude or spacecraft model is absent. The optional
`gravity_gradient=false` setting disables the torque.

Pair this with a gravity force model, including `GravitationalHarmonicsModel`,
without adding another gravitational force. Do not also enable
`gravity_gradient=true` on an analytic gravity effector in the same run:
that would count the central-field torque twice.
"""
@kwdef struct GravityGradientTorqueModel <: AbstractForceTorqueModel
    # Match the existing analytic-gravity option so both hooks reuse its owner.
    gravity_gradient::Bool = true
end

@inline function wrench(
    model::GravityGradientTorqueModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    torque_body = _gravity_gradient_torque_body(model, x, env.planet)
    return SVector{3, Float64}(0.0, 0.0, 0.0), torque_body
end

# ---------------------------------------------------------------------------
# Per-body gravity gradient for articulated spacecraft (live-pose loads)
# ---------------------------------------------------------------------------

"""
    gravity_gradient_requested(effector) -> Bool

`true` when `effector` asks for a gravity-gradient torque: `gravity_gradient=true` on an analytic gravity
model, or an enabled [`GravityGradientTorqueModel`](@ref). On articulated spacecraft with live-pose loads
any such request becomes the per-body torque of [`body_gravity_gradient_torque_ii`](@ref).
"""
@inline gravity_gradient_requested(::Any)::Bool = false
@inline gravity_gradient_requested(
    model::Union{ConstantGravityModel, InverseSquaredGravityModel, InverseSquaredJ2GravityModel, GravityGradientTorqueModel},
)::Bool = model.gravity_gradient

"""True for the stand-alone [`GravityGradientTorqueModel`](@ref), which the live-pose path replaces by per-body torques."""
@inline is_gravity_gradient_torque_effector(::Any)::Bool = false
@inline is_gravity_gradient_torque_effector(::GravityGradientTorqueModel)::Bool = true

"""
    body_gravity_gradient_torque_ii(R_b, I_b, r_ii, μ) -> torque_ii

Central-field gravity-gradient torque about the COM of one dynamic body: `3μ/r³ r̂ × (R_b I_b R_bᵀ) r̂`, with
`I_b` the body inertia about its COM in the body frame, `R_b` the body-to-inertial rotation and `r_ii` the
body COM position. The result is inertial.
"""
@inline body_gravity_gradient_torque_ii(R_b, I_b, r_ii::SVector{3, Float64}, μ::Float64)::SVector{3, Float64} =
    gravity_gradient(SMatrix{3, 3, Float64, 9}(R_b * I_b * R_b'), r_ii, μ)
