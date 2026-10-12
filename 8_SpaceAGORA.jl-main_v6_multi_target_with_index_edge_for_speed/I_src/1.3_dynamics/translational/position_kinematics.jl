@inline function position_derivative(
    velocity::AbstractVector{<:Real}
)::SVector{3, Float64}
    return SVector{3, Float64}(Float64(velocity[1]), Float64(velocity[2]), Float64(velocity[3]))
end

@inline function zero_position_derivative()::SVector{3, Float64}
    return SVector{3, Float64}(0.0, 0.0, 0.0)
end
