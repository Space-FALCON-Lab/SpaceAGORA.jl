@inline function angular_acceleration(
    omega_body::SVector{3, Float64},
    inertia_tensor::AbstractMatrix{<:Real},
    tau_body::SVector{3, Float64};
    include_gyroscopic::Bool=true,
)::SVector{3, Float64}
    rhs = include_gyroscopic ? (tau_body - cross(omega_body, inertia_tensor * omega_body)) : tau_body
    return SVector{3, Float64}(inertia_tensor \ rhs)
end
