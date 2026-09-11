@inline function body_torque(torque::AbstractVector{<:Real})::SVector{3, Float64}
    return SVector{3, Float64}(Float64(torque[1]), Float64(torque[2]), Float64(torque[3]))
end

@inline function body_angular_velocity(omega::AbstractVector{<:Real})::SVector{3, Float64}
    return SVector{3, Float64}(Float64(omega[1]), Float64(omega[2]), Float64(omega[3]))
end

@inline function combine_torques(
    dynamic_torque::SVector{3, Float64},
    control_torque::SVector{3, Float64}
)::SVector{3, Float64}
    return dynamic_torque + control_torque
end
