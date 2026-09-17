@inline function angular_acceleration(
    omega_body::SVector{3, Float64},
    inertia_tensor::AbstractMatrix{<:Real},
    tau_body::SVector{3, Float64};
    include_gyroscopic::Bool=true,
    # Accepts any 3-vector rather than strictly an SVector: the reaction-wheel
    # Jacobian is a plain Matrix (see ReactionWheelAssembly), so `J_rw * h_wheels`
    # yields a Vector at every call site. Matching the permissiveness
    # inertia_tensor already has keeps that a callee detail instead of making
    # every caller convert.
    h_wheel_body::AbstractVector{<:Real}=SVector{3, Float64}(0.0, 0.0, 0.0),
)::SVector{3, Float64}
    h_wheel = SVector{3, Float64}(h_wheel_body)
    rhs = include_gyroscopic ?
        (tau_body - cross(omega_body, inertia_tensor * omega_body + h_wheel)) :
        tau_body
    return SVector{3, Float64}(inertia_tensor \ rhs)
end
