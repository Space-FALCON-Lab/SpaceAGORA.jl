"""
    gravity_gradient(J::SMatrix{3,3,Float64}, rVec::SVector{3,Float64}, μ::Float64)

Calculate the body-frame gravity-gradient torque for a rigid spacecraft in a
central field.

# Args
- `J`: The spacecraft inertia tensor expressed in the body frame [kg·m²].
- `rVec`: The spacecraft position vector expressed in the same body frame [m].
- `μ`: The gravitational parameter of the central body [m³/s²].

# Returns
- A 3-element `SVector` representing the gravity-gradient torque in the body
  frame `[τ_x, τ_y, τ_z]` [N·m].
"""
function gravity_gradient(J::SMatrix{3,3,Float64}, rVec::SVector{3,Float64}, μ::Float64)
    r = norm(rVec)
    if !isfinite(r) || r <= 0.0
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    r_hat = rVec / r
    return 3*μ/r^3 * cross(r_hat, J * r_hat)
end
