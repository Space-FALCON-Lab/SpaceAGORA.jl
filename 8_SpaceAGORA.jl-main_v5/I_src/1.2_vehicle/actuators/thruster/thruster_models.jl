@kwdef struct BaseThrusterModel <: AbstractThrusterModel
    thrust::Vector{Float64}
    direction::Vector{Float64}
    Δv::Vector{Float64}
    start_burn_time::Vector{Float64}
    stop_burn_time::Vector{Float64}
    Isp::Vector{Float64}

    function BaseThrusterModel(
        thrust::Vector{Float64},
        direction::Vector{Float64},
        Δv::Vector{Float64},
        start_burn_time::Vector{Float64},
        stop_burn_time::Vector{Float64},
        Isp::Vector{Float64}
    )
        n = length(thrust)
        if length(direction) != n || length(Δv) != n || length(start_burn_time) != n || length(stop_burn_time) != n || length(Isp) != n
            throw(ArgumentError("BaseThrusterModel vector lengths must match; got (thrust=$(length(thrust)), direction=$(length(direction)), Δv=$(length(Δv)), start_burn_time=$(length(start_burn_time)), stop_burn_time=$(length(stop_burn_time)), Isp=$(length(Isp)))"))
        end

        return new(thrust, direction, Δv, start_burn_time, stop_burn_time, Isp)
    end
end

"""
    SixAxisThrusterModel

Fixed-direction six-thruster CubeSat actuator description. Directions and
locations are expressed in the spacecraft body frame. Control logic should use
this as hardware data and perform allocation in a control effector.
"""
struct SixAxisThrusterModel <: AbstractThrusterModel
    directions_body::SMatrix{3, 6, Float64}
    locations_body::SMatrix{3, 6, Float64}
    max_thrust_n::SVector{6, Float64}
    isp_s::SVector{6, Float64}
    min_firing_time_s::Float64
end

function SixAxisThrusterModel(;
    directions_body = SMatrix{3, 6, Float64}(
         1.0,  0.0,  0.0,
        -1.0,  0.0,  0.0,
         0.0,  1.0,  0.0,
         0.0, -1.0,  0.0,
         0.0,  0.0,  1.0,
         0.0,  0.0, -1.0,
    ),
    locations_body = SMatrix{3, 6, Float64}(zeros(3, 6)),
    max_thrust_n = SVector{6, Float64}(ones(6)),
    isp_s = SVector{6, Float64}(fill(60.0, 6)),
    min_firing_time_s::Float64 = 0.0,
)
    dirs = SMatrix{3, 6, Float64}(directions_body)
    locs = SMatrix{3, 6, Float64}(locations_body)
    max_thrust = SVector{6, Float64}(max_thrust_n)
    isp = SVector{6, Float64}(isp_s)

    cols = ntuple(j -> begin
        d = SVector{3, Float64}(dirs[:, j])
        d_norm = norm(d)
        d_norm <= eps(Float64) && throw(ArgumentError("SixAxisThrusterModel direction $j is zero."))
        d / d_norm
    end, 6)
    normalized_dirs = SMatrix{3, 6, Float64}(hcat(cols...))

    any(x -> !isfinite(x) || x < 0.0, max_thrust) &&
        throw(ArgumentError("SixAxisThrusterModel max_thrust_n values must be finite and nonnegative."))
    any(x -> !isfinite(x) || x <= 0.0, isp) &&
        throw(ArgumentError("SixAxisThrusterModel isp_s values must be finite and positive."))
    min_firing_time_s < 0.0 &&
        throw(ArgumentError("SixAxisThrusterModel min_firing_time_s must be nonnegative."))

    return SixAxisThrusterModel(normalized_dirs, locs, max_thrust, isp, min_firing_time_s)
end
