@kwdef mutable struct BaseThrusterModel <: AbstractThrusterModel
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
