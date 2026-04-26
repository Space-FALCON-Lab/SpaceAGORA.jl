@inline function quaternion_derivative(
    omega_body::SVector{3, Float64},
    q::AbstractVector{<:Real}
)::SVector{4, Float64}
    qv = SVector{3, Float64}(Float64(q[1]), Float64(q[2]), Float64(q[3]))
    qs = Float64(q[4])

    vector_part = qs * omega_body + cross(omega_body, qv)
    scalar_part = -dot(omega_body, qv)

    return 0.5 * SVector{4, Float64}(
        vector_part[1],
        vector_part[2],
        vector_part[3],
        scalar_part,
    )
end
